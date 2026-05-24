"""
tests/conftest.py — shared synthetic-FITS fixtures.

Provides pytest fixtures that build well-formed FITS files in a temporary
directory so unit tests run offline and deterministically. Each fixture
returns a filesystem path; reading is the test's responsibility (so we
exercise the same load path users hit).

Fixtures
========
tiny_image_fits
    64x64 image HDU with a valid celestial WCS (TAN projection) and a
    bright synthetic source near the centre. TELESCOP/INSTRUME populated
    so io.inspect() has something to report.

tiny_events_fits
    1000-row EVENTS BinTable HDU (Chandra-like ENERGY column in eV) plus
    a GTI extension with a single 1 ks interval. Headers carry TELESCOP,
    INSTRUME, EXPOSURE.

tiny_pha_fits
    OGIP SPECTRUM extension carrying a synthetic power-law-like channel
    distribution (CHANNEL, COUNTS, STAT_ERR), 256 channels.

References
==========
- FITS Standard v4.0, IAU FITS Working Group 2018.
- OGIP/92-007 (PHA file format), George & Yusaf 1992.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS


@pytest.fixture(scope="session")
def rng() -> np.random.Generator:
    """Deterministic NumPy RNG for any randomized test."""
    return np.random.default_rng(42)


@pytest.fixture
def tiny_image_fits(tmp_path: Path, rng: np.random.Generator) -> Path:
    """Build a 64x64 image FITS with WCS and a centred Gaussian source."""
    nx = ny = 64
    yy, xx = np.mgrid[0:ny, 0:nx]
    cx, cy = nx / 2.0, ny / 2.0
    sigma = 3.0
    source = 1000.0 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sigma**2))
    bg = rng.normal(loc=5.0, scale=1.0, size=(ny, nx))
    data = (source + bg).astype(np.float32)

    w = WCS(naxis=2)
    w.wcs.crpix = [cx + 1, cy + 1]   # 1-indexed FITS convention
    w.wcs.cdelt = [-1.0 / 3600.0, 1.0 / 3600.0]  # ~1 arcsec/pixel
    w.wcs.crval = [40.66962, -0.01328]           # NGC 1068 approx ICRS
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    hdr = w.to_header()
    hdr["TELESCOP"] = "SYNTHETIC"
    hdr["INSTRUME"] = "TEST"
    hdr["EXPOSURE"] = 1000.0
    hdr["BUNIT"] = "count"

    primary = fits.PrimaryHDU(data=data, header=hdr)
    path = tmp_path / "tiny_image.fits"
    primary.writeto(path, overwrite=True)
    return path


@pytest.fixture
def tiny_events_fits(tmp_path: Path, rng: np.random.Generator) -> Path:
    """Build a Chandra-like event list with EVENTS + GTI extensions."""
    n_events = 1000
    t_start = 0.0
    t_end = 1000.0   # 1 ks
    times = np.sort(rng.uniform(t_start, t_end, size=n_events))

    # Sky-frame X, Y in detector pixels, centred at (4096, 4096).
    x = rng.normal(loc=4096.0, scale=8.0, size=n_events).astype(np.float32)
    y = rng.normal(loc=4096.0, scale=8.0, size=n_events).astype(np.float32)
    # Energies in eV (Chandra convention), spread across 0.5-8 keV.
    energies = rng.uniform(500.0, 8000.0, size=n_events).astype(np.float32)

    col_time = fits.Column(name="TIME",   array=times,    format="D", unit="s")
    col_x    = fits.Column(name="X",      array=x,        format="E")
    col_y    = fits.Column(name="Y",      array=y,        format="E")
    col_e    = fits.Column(name="ENERGY", array=energies, format="E", unit="eV")
    events_hdu = fits.BinTableHDU.from_columns([col_time, col_x, col_y, col_e],
                                               name="EVENTS")
    events_hdu.header["TELESCOP"] = "CHANDRA"
    events_hdu.header["INSTRUME"] = "ACIS"
    events_hdu.header["EXPOSURE"] = t_end - t_start
    events_hdu.header["TSTART"]   = t_start
    events_hdu.header["TSTOP"]    = t_end

    gti_start = np.array([t_start], dtype=np.float64)
    gti_stop  = np.array([t_end],   dtype=np.float64)
    gti_hdu = fits.BinTableHDU.from_columns([
        fits.Column(name="START", array=gti_start, format="D", unit="s"),
        fits.Column(name="STOP",  array=gti_stop,  format="D", unit="s"),
    ], name="GTI")

    hdul = fits.HDUList([fits.PrimaryHDU(), events_hdu, gti_hdu])
    path = tmp_path / "tiny_events.fits"
    hdul.writeto(path, overwrite=True)
    return path


@pytest.fixture
def gapped_events_fits(tmp_path: Path, rng: np.random.Generator) -> Path:
    """Events split across two GTIs with a 1 ks gap between them.

    Layout
    ------
        GTI 1:     0      ->   1000 s
        gap:    1000      ->   2000 s   (NO events here)
        GTI 2:  2000      ->   3000 s

    Used by M4 tests to confirm bin_events_to_lightcurve drops the
    in-gap bin rather than rendering it as a zero-count bin.
    """
    n_events = 500
    n_a = n_events // 2
    n_b = n_events - n_a
    times = np.concatenate([
        np.sort(rng.uniform(0.0,    1000.0, size=n_a)),
        np.sort(rng.uniform(2000.0, 3000.0, size=n_b)),
    ])
    x = rng.normal(loc=4096.0, scale=8.0, size=times.size).astype(np.float32)
    y = rng.normal(loc=4096.0, scale=8.0, size=times.size).astype(np.float32)
    energies = rng.uniform(500.0, 8000.0, size=times.size).astype(np.float32)

    events_hdu = fits.BinTableHDU.from_columns([
        fits.Column(name="TIME",   array=times,    format="D", unit="s"),
        fits.Column(name="X",      array=x,        format="E"),
        fits.Column(name="Y",      array=y,        format="E"),
        fits.Column(name="ENERGY", array=energies, format="E", unit="eV"),
    ], name="EVENTS")
    events_hdu.header["TELESCOP"] = "CHANDRA"
    events_hdu.header["INSTRUME"] = "ACIS"
    events_hdu.header["TSTART"]   = 0.0
    events_hdu.header["TSTOP"]    = 3000.0
    events_hdu.header["EXPOSURE"] = 2000.0

    gti_hdu = fits.BinTableHDU.from_columns([
        fits.Column(name="START", array=np.array([0.0,    2000.0]), format="D", unit="s"),
        fits.Column(name="STOP",  array=np.array([1000.0, 3000.0]), format="D", unit="s"),
    ], name="GTI")

    hdul = fits.HDUList([fits.PrimaryHDU(), events_hdu, gti_hdu])
    path = tmp_path / "gapped_events.fits"
    hdul.writeto(path, overwrite=True)
    return path


@pytest.fixture
def tiny_pha_fits(tmp_path: Path, rng: np.random.Generator) -> Path:
    """Build a 256-channel OGIP PHA-like SPECTRUM HDU with a power-law shape."""
    n_chan = 256
    channels = np.arange(n_chan, dtype=np.int32)
    # N(ch) ~ (ch + 1)^(-1.8), Poisson-realized
    expected = 1.0e5 * np.power(channels + 1.0, -1.8)
    counts = rng.poisson(expected).astype(np.int32)
    errors = np.sqrt(np.maximum(counts, 1.0)).astype(np.float32)

    cols = [
        fits.Column(name="CHANNEL",  array=channels, format="J"),
        fits.Column(name="COUNTS",   array=counts,   format="J"),
        fits.Column(name="STAT_ERR", array=errors,   format="E"),
    ]
    spec_hdu = fits.BinTableHDU.from_columns(cols, name="SPECTRUM")
    spec_hdu.header["TELESCOP"] = "SYNTHETIC"
    spec_hdu.header["INSTRUME"] = "TEST"
    spec_hdu.header["EXPOSURE"] = 1000.0
    spec_hdu.header["HDUCLAS1"] = "SPECTRUM"
    spec_hdu.header["POISSERR"] = True

    hdul = fits.HDUList([fits.PrimaryHDU(), spec_hdu])
    path = tmp_path / "tiny_spectrum.pha"
    hdul.writeto(path, overwrite=True)
    return path
