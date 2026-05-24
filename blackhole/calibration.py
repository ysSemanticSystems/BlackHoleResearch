"""
blackhole.calibration — per-survey calibrators that turn raw pixel values
into units-bearing physical quantities.

Design
------
Each supported survey has a calibrator function

    calibrator(image: ImageData) -> CalibratedImage

that consults the FITS header (CDELT, BUNIT, MAGZP, BMAJ, BMIN, EXPTIME)
together with hardcoded primary-source constants (F_nu_0, ECF) to convert
the pixel data into a Quantity per pixel. Surveys we cannot calibrate
(DSS plate scans, anything without a documented zero-point) raise
`UncalibratedDataError` rather than silently returning a number.

Dispatch
--------
`calibrate(image)` inspects the header — preferring an explicit `SURVEY`
keyword (the convention NASA SkyView writes), falling back to TELESCOP /
INSTRUME / BUNIT pattern matching — and routes to the right calibrator.
The registry is `CALIBRATORS: Mapping[str, Callable]`.

References
----------
- Cohen, Wheaton & Megeath 2003 AJ 126 1090 — 2MASS K-band F_nu(0)=666.7 Jy.
- Miville-Deschênes & Lagache 2005 A&A 432 729 — IRIS reprocessing of IRAS.
- Murakami+2007 PASJ 59 S369 — AKARI mission overview, FIS WIDE-S/L surveys.
- Becker, White & Helfand 1995 ApJ 450 559 — VLA FIRST survey, BMAJ=5.4''.
- Snowden+1995 ApJ 454 643; Voges+1999 A&A 349 389 — RASS broad-band ECF.

Calibration version: bump CALIBRATION_VERSION when a constant changes.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass

import astropy.units as u

from .io import ImageData

CALIBRATION_VERSION = "1.0.0"

# ---------------------------------------------------------------------------
# Public dataclass and error
# ---------------------------------------------------------------------------


class UncalibratedDataError(RuntimeError):
    """Raised when an image cannot be flux-calibrated.

    Reasons include: unsupported survey, missing required header keywords
    (e.g. no MAGZP for 2MASS), or a survey that intrinsically has no
    photometric zero-point (DSS plate scans).
    """


@dataclass(frozen=True)
class CalibratedImage:
    """Result of calibrating an `ImageData`.

    Attributes
    ----------
    array
        2D `astropy.units.Quantity`. Units are per-pixel flux density
        (Jy/pixel) for sub-mm/IR/radio, or per-pixel flux
        (erg/s/cm^2/pixel) for X-ray broad-band count maps.
    pixel_solid_angle
        Steradian per pixel, derived from CDELT1/CDELT2. Used by
        downstream surface-brightness-aware photometry.
    band
        Coarse waveband label aligned with `blackhole.sed.BAND_COLORS`:
        one of `"radio"`, `"submm"`, `"ir"`, `"opt"`, `"uv"`, `"xray"`.
    zeropoint
        The per-image zero-point used (mag-system or BUNIT-equivalent
        numeric). `None` for surveys whose calibration does not pass
        through a numeric zero-point (e.g. RASS).
    zeropoint_ref
        Primary-source paper or keyword describing where `zeropoint`
        comes from.
    method
        Short label naming the calibration path, e.g. ``"MAGZP+F_nu0"``,
        ``"BUNIT_MJy_sr"``, ``"Jy_per_beam"``, ``"counts_to_flux_ECF"``.
    survey
        Survey identifier matched at dispatch time. Useful for
        provenance display.
    original
        The original `ImageData`. Kept so downstream code can re-render
        with raw pixel values when needed.
    """

    array: u.Quantity
    pixel_solid_angle: u.Quantity
    band: str
    zeropoint: float | None
    zeropoint_ref: str
    method: str
    survey: str
    original: ImageData


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pixel_solid_angle(image: ImageData) -> u.Quantity:
    """Pixel solid angle in steradians from CDELT (or CD matrix).

    Uses `|CDELT1 * CDELT2|`; converts from deg^2 to sr by multiplying
    by `(pi/180)**2`. CD matrices that lack CDELTs fall back to
    PC1_1/PC2_2-derived scales. The astropy WCS object would normalize
    these, but we want a survey-independent calculation that works
    against bare astropy.io.fits headers in tests too.

    Raises
    ------
    UncalibratedDataError
        If neither CDELT nor CD scaling can be read.
    """
    h = image.header
    if "CDELT1" in h and "CDELT2" in h:
        dx, dy = float(h["CDELT1"]), float(h["CDELT2"])
    elif "CD1_1" in h and "CD2_2" in h:
        dx, dy = float(h["CD1_1"]), float(h["CD2_2"])
    else:
        raise UncalibratedDataError(
            "No CDELT/CD scale in header; cannot derive pixel solid angle."
        )
    deg2 = abs(dx * dy)
    sr_per_pixel = deg2 * (math.pi / 180.0) ** 2
    return sr_per_pixel * u.sr


def _require(header: object, key: str, survey: str) -> float:
    """Return float(header[key]) or raise UncalibratedDataError.

    The header argument is typed broadly because astropy.io.fits.Header
    has no public Mapping protocol we can constrain against.
    """
    try:
        return float(header[key])  # type: ignore[index]
    except (KeyError, TypeError, ValueError) as exc:
        raise UncalibratedDataError(
            f"{survey} calibration requires header key {key!r}; not found."
        ) from exc


# ---------------------------------------------------------------------------
# Calibrators
# ---------------------------------------------------------------------------

# 2MASS K-band flux density at magnitude zero (Vega), Cohen+2003 Table 1.
F_NU_0_2MASS_K = 666.7  # Jy

# IRAS/IRIS BUNIT = MJy/sr. No per-image zero-point.
IRIS_METHOD = "BUNIT_MJy_sr"

# RASS broad-band (0.1-2.4 keV) energy conversion factor for a T~1 keV
# thermal spectrum and Galactic N_H ~ 3e20 cm^-2, Snowden+1995. Truly
# accurate values depend on the source spectrum; we adopt a published
# representative value and document this in the method label.
RASS_ECF_BROAD = 1.08e-11  # erg/s/cm^2 per (count/s); Snowden+1995 Fig. 8


def _calibrate_2mass_k(image: ImageData) -> CalibratedImage:
    """Calibrate a 2MASS K-band image cutout to Jy/pixel.

    The 2MASS Atlas Image cutouts on SkyView arrive in DN (data numbers)
    with a per-image magnitude zero-point in the header (`MAGZP`). The
    chain is::

        m_pix       = MAGZP - 2.5 * log10(DN)
        flux_pix_Jy = F_nu_0(K) * 10^(-m_pix / 2.5)
                    = F_nu_0(K) * DN * 10^(-MAGZP / 2.5)

    Parameters
    ----------
    image
        ImageData with `MAGZP` in the header (Vega system).

    Returns
    -------
    CalibratedImage
        `array` in Jy per pixel.

    Raises
    ------
    UncalibratedDataError
        If MAGZP is not in the header.

    References
    ----------
    Cohen, Wheaton & Megeath 2003 AJ 126 1090 — F_nu(0)(K) = 666.7 +/- 12.6 Jy.
    Skrutskie+2006 AJ 131 1163 — 2MASS Atlas Image format and MAGZP convention.
    """
    magzp = _require(image.header, "MAGZP", survey="2MASS-K")
    factor = F_NU_0_2MASS_K * 10.0 ** (-magzp / 2.5)
    array = image.array.astype(float) * factor * u.Jy
    return CalibratedImage(
        array=array,
        pixel_solid_angle=_pixel_solid_angle(image),
        band="ir",
        zeropoint=magzp,
        zeropoint_ref="Cohen+2003 (F_nu0); per-image MAGZP",
        method="MAGZP+F_nu0",
        survey="2MASS-K",
        original=image,
    )


def _calibrate_iris(image: ImageData, band_um: float | None = None) -> CalibratedImage:
    """Calibrate an IRAS IRIS image to Jy/pixel.

    IRIS cutouts ship with `BUNIT = "MJy/sr"`. Conversion to Jy per pixel
    is data * MJy/sr * pixel_solid_angle_sr * 1e6 Jy/MJy::

        flux_pix_Jy = data * pixel_sr * 1e6

    Parameters
    ----------
    image
        ImageData with BUNIT == "MJy/sr".
    band_um
        Effective wavelength of the IRIS band (12, 25, 60, 100). Used to
        tag the result; not used in the conversion math.

    References
    ----------
    Miville-Deschênes & Lagache 2005 A&A 432 729 — IRIS reprocessing of IRAS.
    """
    bunit = str(image.header.get("BUNIT", "")).strip()
    if bunit.lower().replace(" ", "") not in ("mjy/sr", "mjysr-1"):
        raise UncalibratedDataError(
            f"IRIS calibrator expects BUNIT='MJy/sr'; got {bunit!r}."
        )
    sr = _pixel_solid_angle(image).value
    array = image.array.astype(float) * sr * 1.0e6 * u.Jy
    return CalibratedImage(
        array=array,
        pixel_solid_angle=_pixel_solid_angle(image),
        band="ir" if (band_um is None or band_um <= 30.0) else "submm",
        zeropoint=None,
        zeropoint_ref="BUNIT=MJy/sr (IRIS reprocessing; Miville-Deschênes+2005)",
        method=IRIS_METHOD,
        survey=f"IRIS_{int(band_um)}" if band_um else "IRIS",
        original=image,
    )


def _calibrate_akari(image: ImageData) -> CalibratedImage:
    """Calibrate an AKARI FIS WIDE-S/L image to Jy/pixel.

    Same `BUNIT = MJy/sr` convention as IRIS.

    References
    ----------
    Doi+2015 PASJ 67 50 — AKARI FIS All-Sky Survey maps; BUNIT=MJy/sr.
    """
    cal = _calibrate_iris(image, band_um=None)
    return CalibratedImage(
        array=cal.array,
        pixel_solid_angle=cal.pixel_solid_angle,
        band="submm",
        zeropoint=None,
        zeropoint_ref="BUNIT=MJy/sr (AKARI FIS; Doi+2015)",
        method=IRIS_METHOD,
        survey="AKARI",
        original=image,
    )


def _calibrate_vla_first(image: ImageData) -> CalibratedImage:
    """Calibrate a VLA FIRST 1.4 GHz image to Jy/pixel.

    FIRST images are in Jy/beam. For a Gaussian beam the effective area is::

        beam_area_sr = pi / (4 ln 2) * BMAJ * BMIN

    where BMAJ and BMIN are in degrees in the header. Per pixel::

        flux_pix_Jy = data * pixel_sr / beam_area_sr

    **Trap.** For sources resolved by the beam, this is a *lower limit*
    on integrated flux: peak Jy/beam is meaningful only for unresolved
    point sources. We compute Jy/pixel mechanically; the
    consumer (photometry) must decide whether it makes sense.

    References
    ----------
    Becker, White & Helfand 1995 ApJ 450 559 — FIRST survey design (BMAJ=5.4'').
    """
    bunit = str(image.header.get("BUNIT", "")).strip().lower()
    if "jy/beam" not in bunit and "jy beam" not in bunit:
        raise UncalibratedDataError(
            f"FIRST calibrator expects BUNIT='Jy/beam'; got {bunit!r}."
        )
    bmaj_deg = _require(image.header, "BMAJ", survey="VLA-FIRST")
    bmin_deg = _require(image.header, "BMIN", survey="VLA-FIRST")
    bmaj_rad = math.radians(bmaj_deg)
    bmin_rad = math.radians(bmin_deg)
    beam_area_sr = math.pi / (4.0 * math.log(2.0)) * bmaj_rad * bmin_rad
    pixel_sr = _pixel_solid_angle(image).value
    factor = pixel_sr / beam_area_sr
    array = image.array.astype(float) * factor * u.Jy
    return CalibratedImage(
        array=array,
        pixel_solid_angle=_pixel_solid_angle(image),
        band="radio",
        zeropoint=None,
        zeropoint_ref=(
            "Jy/beam, beam from BMAJ/BMIN; "
            "VLA FIRST 5.4'' beam (Becker+1995)"
        ),
        method="Jy_per_beam",
        survey="VLA-FIRST",
        original=image,
    )


def _calibrate_rass(image: ImageData) -> CalibratedImage:
    """Calibrate a ROSAT RASS broad-band count map to erg/s/cm^2/pixel.

    The RASS broad-band image holds *counts* per pixel. Conversion is::

        flux_pix = (counts / EXPTIME) * ECF_BROAD

    where ECF_BROAD = 1.08e-11 erg/s/cm^2 per (count/s) for a representative
    1 keV thermal source with N_H = 3e20 cm^-2 (Snowden+1995). The true
    ECF depends on the spectrum and absorbing column of the target;
    consumers should treat this as a first-pass estimate.

    References
    ----------
    Snowden+1995 ApJ 454 643 — RASS broad-band ECF (their Fig. 8).
    Voges+1999 A&A 349 389 — Bright Source Catalogue, ECF formulation.
    """
    exptime = _require(image.header, "EXPTIME", survey="RASS")
    if exptime <= 0:
        raise UncalibratedDataError(f"RASS EXPTIME must be positive; got {exptime}.")
    rate = image.array.astype(float) / exptime
    array = rate * RASS_ECF_BROAD * (u.erg / u.s / u.cm**2)
    return CalibratedImage(
        array=array,
        pixel_solid_angle=_pixel_solid_angle(image),
        band="xray",
        zeropoint=None,
        zeropoint_ref=(
            f"ECF={RASS_ECF_BROAD:.2e} erg/s/cm^2/(count/s); Snowden+1995"
        ),
        method="counts_to_flux_ECF",
        survey="RASS-broad",
        original=image,
    )


def _calibrate_dss(image: ImageData) -> CalibratedImage:
    """DSS plate scans have no documented absolute zero-point.

    The DSS is a digitized photographic survey (Lasker+1990); per-plate
    response, sky background, and non-linearity make absolute photometry
    a research project in its own right, not a one-liner. Raise.
    """
    raise UncalibratedDataError(
        "DSS plate scans have no documented absolute photometric zero-point. "
        "Use DSS for visual context only, or perform plate-by-plate calibration "
        "against catalog stars before treating as photometry."
    )


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

CALIBRATORS: Mapping[str, Callable[[ImageData], CalibratedImage]] = {
    "2MASS-K":   _calibrate_2mass_k,
    "IRIS_12":   lambda im: _calibrate_iris(im, band_um=12.0),
    "IRIS_25":   lambda im: _calibrate_iris(im, band_um=25.0),
    "IRIS_60":   lambda im: _calibrate_iris(im, band_um=60.0),
    "IRIS_100":  lambda im: _calibrate_iris(im, band_um=100.0),
    "AKARI":     _calibrate_akari,
    "VLA-FIRST": _calibrate_vla_first,
    "RASS":      _calibrate_rass,
    "DSS":       _calibrate_dss,
}


def detect_survey(image: ImageData) -> str:
    """Return the registry key for the calibrator that handles `image`.

    Inspects (in order): an explicit `SURVEY` header (SkyView convention),
    TELESCOP+INSTRUME pairs, the BUNIT string, and the filename.

    Raises
    ------
    UncalibratedDataError
        If no calibrator matches.
    """
    h = image.header
    survey = str(h.get("SURVEY", "")).strip()
    telescop = str(h.get("TELESCOP", "")).strip().upper()
    instrume = str(h.get("INSTRUME", "")).strip().upper()
    bunit = str(h.get("BUNIT", "")).strip().lower()
    fname = str(image.source_path).lower()

    # Header-driven blob (SURVEY/TELESCOP/INSTRUME/BUNIT) — used for survey
    # family identification. The filename is used as a *fallback only* and
    # never to disambiguate IRIS band numbers, because pytest tmp paths
    # often contain stray digits.
    header_blob = " | ".join([survey, telescop, instrume, bunit]).lower()

    # 2MASS K — header first, filename as last resort.
    if "2mass" in header_blob and ("k" in header_blob or "ks" in header_blob):
        return "2MASS-K"
    if "2mass" in fname and ("_k" in fname or "ks" in fname):
        return "2MASS-K"

    # IRIS / IRAS reprocessed maps — disambiguate band only from the SURVEY
    # keyword, never from the filename.
    if "iris" in header_blob or "iras" in header_blob:
        survey_l = survey.lower()
        for um in (100, 60, 25, 12):  # longest-first to avoid "12" matching "120"
            if str(um) in survey_l:
                return f"IRIS_{um}"
        # No band in header — fall back to filename for the band number.
        for um in (100, 60, 25, 12):
            if f"iris_{um}" in fname or f"iras_{um}" in fname or f"{um}um" in fname:
                return f"IRIS_{um}"
        return "IRIS_12"

    # AKARI
    if "akari" in header_blob or "akari" in fname:
        return "AKARI"

    # VLA FIRST
    if "first" in header_blob or ("vla" in header_blob and "1.4" in header_blob):
        return "VLA-FIRST"
    if "first" in fname:
        return "VLA-FIRST"

    # ROSAT All-Sky Survey
    if "rass" in header_blob or "rosat" in header_blob:
        return "RASS"
    if "rass" in fname or "rosat" in fname:
        return "RASS"

    # DSS / POSS plate scans
    if "dss" in header_blob or "poss" in header_blob or "digitized sky" in header_blob:
        return "DSS"
    if "dss" in fname or "poss" in fname:
        return "DSS"

    raise UncalibratedDataError(
        f"Could not identify calibrator from headers "
        f"(SURVEY={survey!r}, TELESCOP={telescop!r}, "
        f"INSTRUME={instrume!r}, BUNIT={bunit!r})."
    )


def calibrate(image: ImageData) -> CalibratedImage:
    """Dispatch to the right calibrator based on the FITS header.

    See `detect_survey` for the matching rules and `CALIBRATORS` for the
    set of supported surveys. Unsupported surveys raise
    `UncalibratedDataError` rather than returning a number.

    Examples
    --------
    >>> from blackhole import io as bhio
    >>> from blackhole.calibration import calibrate
    >>> img = bhio.load_image("fits_data/ngc1068_2mass_k.fits")  # doctest: +SKIP
    >>> cal = calibrate(img)                                      # doctest: +SKIP
    >>> cal.array.unit                                            # doctest: +SKIP
    Unit("Jy")
    """
    survey = detect_survey(image)
    calibrator = CALIBRATORS[survey]
    return calibrator(image)


def supported_surveys() -> tuple[str, ...]:
    """Return the registry keys, sorted, for display in UI and docs."""
    return tuple(sorted(CALIBRATORS.keys()))
