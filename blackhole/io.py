"""
blackhole.io — FITS loading and inspection utilities.

WHY THIS MODULE EXISTS
======================
FITS (Flexible Image Transport System) is the standard astronomical data format,
formalized by NASA/IAU in 1981 and still in active use 45 years later. Every
NASA mission from Einstein (1978) through JWST (2021) writes FITS. This module
exists because the format hides three traps that bite every new user:

  TRAP 1: A FITS file is NOT a single image. It's a list of HDUs
          (Header/Data Units), each with its own header and payload. You
          always start by listing the HDUs, then index into the one you want.

  TRAP 2: X-ray "data" from Chandra, NuSTAR, XMM-Newton is NOT a 2D image.
          It is a table of individual detected photons (an "event list"),
          each row a single photon's arrival time, sky position, energy, and
          quality flags. To see a picture, you bin the photons into pixels
          yourself. Confusing event files for images is the #1 new-user error.

  TRAP 3: A FITS file's data block uses big-endian byte order regardless of
          machine architecture. Astropy handles this for you on access, but
          if you cache numpy arrays and operate on them with mixed code paths
          you can get silent corruption. Always work through the astropy API
          unless you have a reason not to.

REFERENCES
==========
- FITS Standard v4.0 (IAU FITS Working Group, 2018):
  https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf
- Astropy FITS I/O guide:
  https://docs.astropy.org/en/stable/io/fits/
- Chandra event-file structure (CIAO docs):
  https://cxc.harvard.edu/ciao/dictionary/eventlist.html
- HEASARC FITS file format primer:
  https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/primer/node5.html
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

# ---------------------------------------------------------------------------
# Inspection
# ---------------------------------------------------------------------------

@dataclass
class HDUSummary:
    """One-line summary of a single HDU. Returned by `inspect`."""

    index: int
    name: str
    hdu_type: str                    # 'PrimaryHDU' | 'ImageHDU' | 'BinTableHDU' | etc.
    shape: tuple[int, ...] | None    # pixel dims for images; (nrows,) for tables
    columns: list[str] | None        # column names for tables, else None
    telescope: str | None
    instrument: str | None
    exposure_s: float | None


def inspect(path: str | Path) -> list[HDUSummary]:
    """List every HDU in a FITS file with the metadata you actually care about.

    This is your first call on any unfamiliar FITS file. It answers:
      - How many HDUs are there?
      - Which are images, which are tables?
      - What instrument and exposure time?

    Example:
        >>> for hdu in inspect("fits_data/ngc1068_chandra.evt2.fits"):
        ...     print(hdu)

    The header keys we look for (TELESCOP, INSTRUME, EXPOSURE) are FITS
    *reserved* keywords with standard meanings — see FITS Standard v4.0 §4.4.2.
    Missing keys return None rather than raising, because real-world FITS
    files from the 1990s sometimes lack them.
    """
    path = Path(path)
    out: list[HDUSummary] = []
    with fits.open(path) as hdul:
        for i, hdu in enumerate(hdul):
            hdr = hdu.header
            is_table = isinstance(hdu, (fits.BinTableHDU, fits.TableHDU))
            shape = None if is_table else (None if hdu.data is None else tuple(hdu.data.shape))
            cols = [c.name for c in hdu.columns] if is_table else None
            if is_table and hdu.data is not None:
                # For tables, "shape" is more useful as row count.
                shape = (len(hdu.data),)
            out.append(HDUSummary(
                index=i,
                name=hdu.name or f"HDU{i}",
                hdu_type=type(hdu).__name__,
                shape=shape,
                columns=cols,
                telescope=hdr.get("TELESCOP"),
                instrument=hdr.get("INSTRUME"),
                exposure_s=hdr.get("EXPOSURE") or hdr.get("ONTIME") or hdr.get("LIVETIME"),
            ))
    return out


def header_dict(path: str | Path, hdu_index: int = 0) -> dict[str, object]:
    """Return one HDU's header as an ordinary dict.

    The native astropy Header object is dict-like but carries comment cards
    and HISTORY entries — useful for forensic work, but noise when you just
    want to dump metadata into a UI table or JSON file.
    """
    with fits.open(path) as hdul:
        h = hdul[hdu_index].header
        # Strip card-level metadata. astropy gives us comments as a parallel
        # mapping; we discard them here on purpose.
        return {k: h[k] for k in h if k and k not in ("COMMENT", "HISTORY", "")}


# ---------------------------------------------------------------------------
# Image loading (Image HDUs)
# ---------------------------------------------------------------------------

@dataclass
class ImageData:
    """A 2D astronomical image plus its World Coordinate System."""

    array: np.ndarray
    wcs: WCS | None
    header: fits.Header
    source_path: Path
    hdu_index: int


def load_image(path: str | Path, hdu_index: int | None = None) -> ImageData:
    """Load a 2D image and its WCS.

    If `hdu_index` is None we pick the first HDU that has a 2D data array.
    The Primary HDU is often empty (NAXIS=0) when secondary extensions hold
    the actual image — a common convention since the 1990s.

    PITFALL: We `.copy()` the array because astropy memory-maps FITS data
    by default. When the `with` block closes the file, the memmap is
    invalidated and any further access on the unwrapped numpy array
    produces undefined behavior (sometimes garbage, sometimes a segfault).
    Copying is cheap for typical-size images (< 100 MB) and removes the
    foot-gun entirely. See:
    https://docs.astropy.org/en/stable/io/fits/appendix/faq.html#i-m-getting-a-segfault
    """
    path = Path(path)
    with fits.open(path) as hdul:
        if hdu_index is None:
            for i, h in enumerate(hdul):
                if h.data is not None and h.data.ndim == 2:
                    hdu_index = i
                    break
            if hdu_index is None:
                raise ValueError(f"No 2D image HDU found in {path.name}")
        hdu = hdul[hdu_index]
        if hdu.data is None or hdu.data.ndim != 2:
            raise ValueError(
                f"HDU {hdu_index} of {path.name} is not a 2D image "
                f"(shape={None if hdu.data is None else hdu.data.shape})"
            )
        array = hdu.data.copy()  # see PITFALL note above
        header = hdu.header.copy()

    # Build WCS. Some FITS files lack WCS keywords entirely (small cutouts
    # from old missions); we return None rather than raising so the caller
    # can decide whether to plot in pixel coordinates.
    #
    # PITFALL (#3 in docs/PITFALLS.md): astropy.wcs.WCS(header) returns a
    # valid 2-axis object even when the header has no CTYPE/CRVAL keywords,
    # silently using the identity mapping. We must verify the WCS actually
    # encodes a celestial projection via has_celestial.
    try:
        wcs = WCS(header)
        if wcs.naxis < 2 or not wcs.has_celestial:
            wcs = None
    except Exception:
        wcs = None

    return ImageData(array=array, wcs=wcs, header=header,
                     source_path=path, hdu_index=hdu_index)


# ---------------------------------------------------------------------------
# Event lists (BinTable HDUs from X-ray missions)
# ---------------------------------------------------------------------------

@dataclass
class EventList:
    """An X-ray photon event list. Each row = one detected photon.

    Standard columns across Chandra/XMM/NuSTAR (column names vary slightly
    by mission, hence the per-mission accessor methods):

      Chandra ACIS evt2:   TIME, X, Y, ENERGY (eV), CCD_ID, GRADE
      XMM EPIC PN evt:     TIME, X, Y, PI (eV), PATTERN, FLAG
      NuSTAR evt:          TIME, X, Y, PI (channels), DET_ID

    The mission and column-name lookup is handled by `bin_to_image`.
    """

    times: np.ndarray         # seconds, mission reference frame
    x: np.ndarray             # detector or sky pixel coordinate
    y: np.ndarray
    energies: np.ndarray | None  # eV (Chandra) or channel-units (others)
    energy_unit: str          # "eV", "PI", etc.
    header: fits.Header
    source_path: Path
    mission: str


def load_events(path: str | Path) -> EventList:
    """Load an X-ray event list.

    We look for the EVENTS extension first (Chandra/XMM standard), then
    fall back to scanning for any BinTable with a TIME column.

    PITFALL: Don't try to plt.imshow() the result. There are no pixels yet.
    Call `bin_to_image` first.
    """
    path = Path(path)
    with fits.open(path) as hdul:
        # Prefer the explicitly-named EVENTS extension if present.
        ev_hdu = None
        for h in hdul:
            if h.name.upper() == "EVENTS" and isinstance(h, fits.BinTableHDU):
                ev_hdu = h
                break
        if ev_hdu is None:
            # Fallback: first BinTable with a TIME column.
            for h in hdul:
                if isinstance(h, fits.BinTableHDU) and "TIME" in [c.name for c in h.columns]:
                    ev_hdu = h
                    break
        if ev_hdu is None:
            raise ValueError(f"No event-list extension found in {path.name}")

        data = ev_hdu.data
        header = ev_hdu.header.copy()

        cols_upper = {c.name.upper(): c.name for c in ev_hdu.columns}
        mission = (header.get("TELESCOP") or "").upper()

        times = np.asarray(data[cols_upper["TIME"]])
        # X/Y might be called X,Y or DETX,DETY or RAWX,RAWY.
        # Sky-frame X,Y is what we want for imaging.
        x_col = cols_upper.get("X") or cols_upper.get("DETX") or cols_upper.get("RAWX")
        y_col = cols_upper.get("Y") or cols_upper.get("DETY") or cols_upper.get("RAWY")
        if not (x_col and y_col):
            raise ValueError(f"No X/Y position columns in {path.name}")
        x = np.asarray(data[x_col]).astype(float)
        y = np.asarray(data[y_col]).astype(float)

        # Energy column: ENERGY (Chandra, eV) or PI (XMM/NuSTAR, channel/eV-ish)
        energies = None
        energy_unit = "unknown"
        if "ENERGY" in cols_upper:
            energies = np.asarray(data[cols_upper["ENERGY"]]).astype(float)
            energy_unit = "eV"
        elif "PI" in cols_upper:
            energies = np.asarray(data[cols_upper["PI"]]).astype(float)
            energy_unit = "PI"

    return EventList(
        times=times, x=x, y=y, energies=energies, energy_unit=energy_unit,
        header=header, source_path=path, mission=mission or "UNKNOWN",
    )


def bin_to_image(
    events: EventList,
    bins: int = 512,
    energy_range_ev: tuple[float, float] | None = None,
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    """Bin an event list into a 2D image (photon count per pixel).

    Returns (image_array, extent) where extent = (xmin, xmax, ymin, ymax)
    is suitable for matplotlib's `imshow(extent=...)`.

    `energy_range_ev` filters by energy before binning. For Chandra in eV:
      Soft:  500-2000
      Hard:  2000-7000
      Broad: 500-7000

    Passing energy_range_ev with PI-channel data (XMM/NuSTAR) will silently
    misbehave — channels are not eV. For those missions, convert PI to eV
    first via the RMF response file, or pass channel ranges that you've
    computed yourself.

    PITFALL: A long Chandra observation might have ~10^6 events spread over
    a wide field. Default 512x512 is a reasonable balance. Push to 1024 for
    detail at the cost of noisier per-pixel statistics.
    """
    x = events.x
    y = events.y
    if energy_range_ev is not None and events.energies is not None and events.energy_unit == "eV":
        mask = (events.energies >= energy_range_ev[0]) & (events.energies <= energy_range_ev[1])
        x = x[mask]
        y = y[mask]

    if x.size == 0:
        # Empty selection — return a zero image so downstream code doesn't crash.
        return np.zeros((bins, bins)), (0.0, 1.0, 0.0, 1.0)

    h, xedges, yedges = np.histogram2d(x, y, bins=bins)
    # histogram2d returns image transposed relative to imshow conventions.
    # We transpose so that pixel (i,j) lines up with matplotlib display.
    image = h.T
    extent = (float(xedges[0]), float(xedges[-1]),
              float(yedges[0]), float(yedges[-1]))
    return image, extent


# ---------------------------------------------------------------------------
# Iterators
# ---------------------------------------------------------------------------

def iter_fits(directory: str | Path) -> Iterator[Path]:
    """Yield every .fits / .fits.gz / .fit file in a directory, sorted."""
    d = Path(directory)
    patterns = ("*.fits", "*.fits.gz", "*.fit", "*.fit.gz", "*.fz")
    found: list[Path] = []
    for pat in patterns:
        found.extend(d.glob(pat))
    yield from sorted(found)
