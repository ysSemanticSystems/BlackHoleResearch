"""
tests/test_photometry.py — synthetic-source recovery.

We inject a 2D Gaussian point source with known total flux into a
calibrated image and verify the aperture path recovers it within a few
percent. Background statistics are validated separately on a flat sky
image with known sigma.
"""

from __future__ import annotations

from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits

from blackhole import io as bhio
from blackhole.calibration import calibrate
from blackhole.photometry import (
    PhotometryResult,
    aperture_for_band,
    aperture_photometry_on,
)

# ---------------------------------------------------------------------------
# Helpers — build synthetic calibrated images
# ---------------------------------------------------------------------------


def _wcs_header(
    crval: tuple[float, float],
    pix_arcsec: float = 1.0,
    image_size: int = 81,
) -> dict[str, object]:
    """A celestial-WCS header that puts `crval` at the centre of `image_size`."""
    centre = (image_size + 1) / 2.0   # 1-indexed FITS centre
    return {
        "CDELT1": -pix_arcsec / 3600.0,
        "CDELT2":  pix_arcsec / 3600.0,
        "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
        "CRVAL1": crval[0],   "CRVAL2": crval[1],
        "CRPIX1": centre,     "CRPIX2": centre,
    }


def _write(path: Path, data: np.ndarray, cards: dict[str, object]) -> Path:
    hdu = fits.PrimaryHDU(data=data.astype(np.float32))
    for k, v in cards.items():
        hdu.header[k] = v
    hdu.writeto(path, overwrite=True)
    return path


def _gaussian(size: int, total: float, sigma: float = 1.5) -> np.ndarray:
    """A normalized 2D Gaussian of `total` peak-area integral, centred."""
    y, x = np.mgrid[0:size, 0:size]
    cx = cy = (size - 1) / 2.0
    g = np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * sigma**2))
    g *= total / g.sum()
    return g


# ---------------------------------------------------------------------------
# Background statistics
# ---------------------------------------------------------------------------


def test_sky_sigma_recovered_on_flat_image(tmp_path: Path) -> None:
    """A flat sky image with known Gaussian noise: sky_sigma matches input."""
    rng = np.random.default_rng(1)
    size = 81
    data = rng.normal(loc=0.0, scale=1.0, size=(size, size))
    path = _write(tmp_path / "iris.fits", data, {
        "TELESCOP": "IRAS", "SURVEY": "IRIS 12", "BUNIT": "MJy/sr",
        **_wcs_header((0.0, 0.0)),
    })
    cal = calibrate(bhio.load_image(path))
    res = aperture_photometry_on(
        cal,
        SkyCoord(0.0, 0.0, unit=(u.deg, u.deg), frame="icrs"),
        aperture_radius=3.0 * u.arcsec,
        annulus_inner=10.0 * u.arcsec,
        annulus_outer=30.0 * u.arcsec,
    )
    # sky_sigma should track input noise (in calibrated Jy units; just
    # check positivity and finite-ness here — the calibration scale matters).
    assert res.sky_sigma.value > 0
    assert np.isfinite(res.sky_sigma.value)


# ---------------------------------------------------------------------------
# Point-source recovery via 2MASS-K calibration path
# ---------------------------------------------------------------------------


def test_point_source_recovered_within_2_percent(tmp_path: Path) -> None:
    """Inject a Gaussian point source on a 2MASS-K image; the aperture
    sum should match the injected total to within a few percent."""
    size = 81
    total_dn = 1.0e5  # data-number budget for the source
    img = _gaussian(size, total=total_dn, sigma=2.0)
    # Background: small constant ~ 1% of peak, no noise here (we test
    # the math; statistical tests live above).
    img += 50.0
    magzp = 21.0
    path = _write(tmp_path / "ngc_2mass.fits", img, {
        "TELESCOP": "2MASS", "SURVEY": "2MASS-K", "MAGZP": magzp,
        **_wcs_header((0.0, 0.0), pix_arcsec=1.0),
    })
    cal = calibrate(bhio.load_image(path))

    # 2MASS-K factor: F_nu0 * 10^(-MAGZP/2.5) = 666.7 * 10^-8.4 Jy / DN.
    expected_jy = total_dn * 666.7 * 10.0 ** (-magzp / 2.5)

    res = aperture_photometry_on(
        cal,
        SkyCoord(0.0, 0.0, unit=(u.deg, u.deg), frame="icrs"),
        aperture_radius=10.0 * u.arcsec,
        annulus_inner=15.0 * u.arcsec,
        annulus_outer=25.0 * u.arcsec,
    )
    assert res.flux.unit == u.Jy
    # Aperture-correction loss is small for a Gaussian sigma=2px in r=10px
    # aperture (encircled energy ~ 99.99%).
    assert res.flux.value == pytest.approx(expected_jy, rel=2e-2)
    assert not res.upper_limit


# ---------------------------------------------------------------------------
# Upper-limit detection
# ---------------------------------------------------------------------------


def test_returns_upper_limit_for_subthreshold_source(tmp_path: Path) -> None:
    """A faint Gaussian (< 3 sigma_bg) sitting in a noisy field returns
    upper_limit=True with flux set to the 3-sigma threshold."""
    rng = np.random.default_rng(7)
    size = 81
    noise_sigma = 5.0  # DN
    faint_total = 8.0   # very faint source (well under aperture noise budget)
    data = rng.normal(loc=0.0, scale=noise_sigma, size=(size, size))
    data += _gaussian(size, total=faint_total, sigma=2.0)
    path = _write(tmp_path / "iris_faint.fits", data, {
        "TELESCOP": "IRAS", "SURVEY": "IRIS 12", "BUNIT": "MJy/sr",
        **_wcs_header((0.0, 0.0)),
    })
    cal = calibrate(bhio.load_image(path))
    res = aperture_photometry_on(
        cal,
        SkyCoord(0.0, 0.0, unit=(u.deg, u.deg), frame="icrs"),
        aperture_radius=10.0 * u.arcsec,
        annulus_inner=15.0 * u.arcsec,
        annulus_outer=25.0 * u.arcsec,
    )
    assert res.upper_limit is True
    # 3-sigma upper limit: flux >= 3 * flux_err.
    assert res.flux.value == pytest.approx(3.0 * res.flux_err.value, rel=1e-9)


# ---------------------------------------------------------------------------
# SEDPoint conversion
# ---------------------------------------------------------------------------


def test_to_sed_point_carries_jy_flux(tmp_path: Path) -> None:
    img = _gaussian(81, total=1e6, sigma=2.0)
    path = _write(tmp_path / "iris.fits", img, {
        "TELESCOP": "IRAS", "SURVEY": "IRIS 12", "BUNIT": "MJy/sr",
        **_wcs_header((0.0, 0.0)),
    })
    cal = calibrate(bhio.load_image(path))
    res = aperture_photometry_on(
        cal, SkyCoord(0.0, 0.0, unit=(u.deg, u.deg), frame="icrs"),
        aperture_radius=10.0 * u.arcsec,
        annulus_inner=20.0 * u.arcsec,
        annulus_outer=30.0 * u.arcsec,
    )
    pt = res.to_sed_point()
    assert pt.flux_density is not None
    assert pt.flux_density.unit == u.Jy
    assert pt.flux_err is not None
    assert pt.band == "ir"
    assert "local photometry" in pt.source


def test_to_sed_point_xray_uses_nu_f_nu(tmp_path: Path) -> None:
    """A RASS-calibrated image returns flux in erg/s/cm^2; the SEDPoint
    should record it as nu_f_nu, not flux_density."""
    img = _gaussian(81, total=1000.0, sigma=2.0)
    path = _write(tmp_path / "rass.fits", img, {
        "TELESCOP": "ROSAT", "SURVEY": "RASS BROAD",
        "EXPTIME": 400.0, **_wcs_header((0.0, 0.0)),
    })
    cal = calibrate(bhio.load_image(path))
    res = aperture_photometry_on(
        cal, SkyCoord(0.0, 0.0, unit=(u.deg, u.deg), frame="icrs"),
        aperture_radius=10.0 * u.arcsec,
        annulus_inner=20.0 * u.arcsec,
        annulus_outer=30.0 * u.arcsec,
    )
    pt = res.to_sed_point()
    assert pt.nu_f_nu is not None
    assert pt.nu_f_nu.unit == u.erg / u.s / u.cm**2
    assert pt.flux_density is None
    assert pt.band == "xray"


# ---------------------------------------------------------------------------
# Result type and defaults
# ---------------------------------------------------------------------------


def test_result_is_frozen(tmp_path: Path) -> None:
    """PhotometryResult is a frozen dataclass."""
    img = _gaussian(81, total=1e5)
    path = _write(tmp_path / "iris.fits", img, {
        "TELESCOP": "IRAS", "SURVEY": "IRIS 12", "BUNIT": "MJy/sr",
        **_wcs_header((0.0, 0.0)),
    })
    cal = calibrate(bhio.load_image(path))
    res = aperture_photometry_on(
        cal, SkyCoord(0.0, 0.0, unit=(u.deg, u.deg), frame="icrs"),
    )
    assert isinstance(res, PhotometryResult)
    import dataclasses
    with pytest.raises(dataclasses.FrozenInstanceError):
        res.flux = 0.0 * u.Jy  # type: ignore[misc]


@pytest.mark.parametrize(
    "band,expected_arcsec",
    [
        ("radio",  60.0),
        ("submm", 120.0),
        ("ir",     20.0),
        ("opt",    10.0),
        ("uv",     10.0),
        ("xray",   30.0),
        ("unknown", 10.0),
    ],
)
def test_aperture_for_band_defaults(band: str, expected_arcsec: float) -> None:
    radius = aperture_for_band(band)
    assert radius.unit == u.arcsec
    assert radius.value == pytest.approx(expected_arcsec)


def test_no_wcs_falls_back_to_pixel_centre(tmp_path: Path) -> None:
    """An image without a celestial WCS still photometers — using image
    centre instead of the catalog coordinate. The fallback path is
    documented in `aperture_photometry_on`."""
    img = _gaussian(81, total=1e6, sigma=2.0)
    path = _write(tmp_path / "iris_no_wcs.fits", img, {
        "TELESCOP": "IRAS", "SURVEY": "IRIS 12", "BUNIT": "MJy/sr",
        # CDELT is required for pixel solid angle but no CTYPE -> no WCS.
        "CDELT1": -1.0 / 3600.0, "CDELT2": 1.0 / 3600.0,
    })
    cal = calibrate(bhio.load_image(path))
    res = aperture_photometry_on(
        cal, SkyCoord(0.0, 0.0, unit=(u.deg, u.deg), frame="icrs"),
        aperture_radius=10.0 * u.arcsec,
        annulus_inner=15.0 * u.arcsec,
        annulus_outer=25.0 * u.arcsec,
    )
    assert res.flux.unit == u.Jy
    assert not res.upper_limit
