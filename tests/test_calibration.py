"""
tests/test_calibration.py — one synthetic FITS per supported survey.

Each test builds a 2-pixel image with header keywords representative of
the real survey, runs `calibrate(image)`, and asserts the resulting
`Quantity` matches a hand-computed expected value to within 1e-4.

The synthetic images are intentionally tiny (2x2) so the math is
inspectable; pixel-solid-angle is computed from CDELT1/CDELT2.
"""

from __future__ import annotations

import math
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits

from blackhole import io as bhio
from blackhole.calibration import (
    CALIBRATION_VERSION,
    CALIBRATORS,
    F_NU_0_2MASS_K,
    RASS_ECF_BROAD,
    CalibratedImage,
    UncalibratedDataError,
    calibrate,
    detect_survey,
    supported_surveys,
)

# ---------------------------------------------------------------------------
# Synthetic-FITS factories
# ---------------------------------------------------------------------------


def _write_image(
    tmp_path: Path,
    name: str,
    data: np.ndarray,
    cards: dict[str, object],
) -> Path:
    """Write a 1-HDU FITS image with the given header cards."""
    hdu = fits.PrimaryHDU(data=data.astype(np.float32))
    for k, v in cards.items():
        hdu.header[k] = v
    path = tmp_path / name
    hdu.writeto(path, overwrite=True)
    return path


def _common_cdelt() -> dict[str, object]:
    """A standard 1-arcsec/pixel WCS scaffold used by most tests."""
    return {
        "CDELT1": -1.0 / 3600.0,
        "CDELT2":  1.0 / 3600.0,
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "CRVAL1": 0.0,
        "CRVAL2": 0.0,
        "CRPIX1": 1.0,
        "CRPIX2": 1.0,
    }


# ---------------------------------------------------------------------------
# Registry sanity
# ---------------------------------------------------------------------------


def test_calibration_version_is_pinned() -> None:
    assert CALIBRATION_VERSION == "1.0.0"


def test_supported_surveys_includes_phase1_set() -> None:
    surveys = set(supported_surveys())
    assert {"2MASS-K", "IRIS_12", "IRIS_25", "IRIS_60", "IRIS_100",
            "AKARI", "VLA-FIRST", "RASS", "DSS"} <= surveys


# ---------------------------------------------------------------------------
# 2MASS-K
# ---------------------------------------------------------------------------


def test_2mass_k_calibration_matches_hand_computed(tmp_path: Path) -> None:
    """A pixel with DN=1, MAGZP=21 gives F_nu_0 * 10^(-21/2.5) Jy.

    666.7 * 10^(-8.4) = 2.6541e-6 Jy.
    """
    data = np.array([[1.0, 2.0]])
    path = _write_image(tmp_path, "2mass.fits", data, {
        "TELESCOP": "2MASS",
        "INSTRUME": "K",
        "SURVEY":   "2MASS-K",
        "MAGZP":    21.0,
        **_common_cdelt(),
    })
    img = bhio.load_image(path)
    cal = calibrate(img)
    assert isinstance(cal, CalibratedImage)
    assert cal.survey == "2MASS-K"
    assert cal.band == "ir"
    assert cal.method == "MAGZP+F_nu0"
    expected = F_NU_0_2MASS_K * 10.0 ** (-21.0 / 2.5)
    assert cal.array.unit == u.Jy
    assert cal.array[0, 0].value == pytest.approx(expected, rel=1e-9)
    assert cal.array[0, 1].value == pytest.approx(2.0 * expected, rel=1e-9)


def test_2mass_k_requires_magzp(tmp_path: Path) -> None:
    path = _write_image(tmp_path, "2mass_no_magzp.fits", np.zeros((2, 2)), {
        "TELESCOP": "2MASS",
        "SURVEY":   "2MASS-K",
        **_common_cdelt(),
    })
    img = bhio.load_image(path)
    with pytest.raises(UncalibratedDataError, match="MAGZP"):
        calibrate(img)


# ---------------------------------------------------------------------------
# IRIS (IRAS reprocessing)
# ---------------------------------------------------------------------------


def test_iris_12_calibration_matches_hand_computed(tmp_path: Path) -> None:
    """MJy/sr * pixel_sr * 1e6 = Jy/pixel.

    With 1 arcsec/pixel: pixel_sr = (1/3600)^2 * (pi/180)^2 = 2.35044e-11 sr.
    A pixel value of 100 MJy/sr → 100 * 2.35044e-11 * 1e6 = 2.35044e-3 Jy.
    """
    data = np.array([[100.0, 200.0]])
    path = _write_image(tmp_path, "iris12.fits", data, {
        "TELESCOP": "IRAS",
        "SURVEY":   "IRIS 12",
        "BUNIT":    "MJy/sr",
        **_common_cdelt(),
    })
    img = bhio.load_image(path)
    cal = calibrate(img)
    sr = (1.0 / 3600.0) ** 2 * (math.pi / 180.0) ** 2
    expected = 100.0 * sr * 1e6
    assert cal.survey == "IRIS_12"
    assert cal.band == "ir"
    assert cal.array.unit == u.Jy
    assert cal.array[0, 0].value == pytest.approx(expected, rel=1e-4)
    assert cal.array[0, 1].value == pytest.approx(2.0 * expected, rel=1e-4)


def test_iris_60_classified_as_submm(tmp_path: Path) -> None:
    path = _write_image(tmp_path, "iris60.fits", np.array([[1.0, 1.0]]), {
        "TELESCOP": "IRAS", "SURVEY": "IRIS 60", "BUNIT": "MJy/sr",
        **_common_cdelt(),
    })
    img = bhio.load_image(path)
    cal = calibrate(img)
    assert cal.band == "submm"


def test_iris_rejects_wrong_bunit(tmp_path: Path) -> None:
    path = _write_image(tmp_path, "iris_wrong.fits", np.zeros((2, 2)), {
        "TELESCOP": "IRAS", "SURVEY": "IRIS 12", "BUNIT": "Jy/beam",
        **_common_cdelt(),
    })
    img = bhio.load_image(path)
    with pytest.raises(UncalibratedDataError, match="MJy/sr"):
        calibrate(img)


# ---------------------------------------------------------------------------
# AKARI
# ---------------------------------------------------------------------------


def test_akari_calibration_uses_iris_path(tmp_path: Path) -> None:
    data = np.array([[10.0, 20.0]])
    path = _write_image(tmp_path, "akari.fits", data, {
        "TELESCOP": "AKARI", "SURVEY": "AKARI WIDE-S",
        "BUNIT": "MJy/sr", **_common_cdelt(),
    })
    img = bhio.load_image(path)
    cal = calibrate(img)
    assert cal.survey == "AKARI"
    assert cal.band == "submm"
    assert cal.array.unit == u.Jy
    sr = (1.0 / 3600.0) ** 2 * (math.pi / 180.0) ** 2
    assert cal.array[0, 0].value == pytest.approx(10.0 * sr * 1e6, rel=1e-4)


# ---------------------------------------------------------------------------
# VLA FIRST
# ---------------------------------------------------------------------------


def test_vla_first_calibration_matches_hand_computed(tmp_path: Path) -> None:
    """For a 5.4 arcsec beam (BMAJ=BMIN=1.5e-3 deg) and 1.8 arcsec/pixel:

        beam_area_sr = pi/(4 ln 2) * (1.5e-3 deg)^2 * (pi/180)^2 sr
        pixel_sr     = (1.8/3600)^2 * (pi/180)^2 sr
        factor       = pixel_sr / beam_area_sr

    A pixel value of 0.01 Jy/beam converts to 0.01 * factor Jy/pixel.
    """
    bmaj_deg = 1.5e-3
    pix_deg = 1.8 / 3600.0
    data = np.array([[0.01, 0.02]])
    path = _write_image(tmp_path, "first.fits", data, {
        "TELESCOP": "VLA",
        "SURVEY":   "VLA FIRST 1.4 GHz",
        "BUNIT":    "Jy/beam",
        "BMAJ":     bmaj_deg,
        "BMIN":     bmaj_deg,
        "CDELT1":  -pix_deg,
        "CDELT2":   pix_deg,
        "CTYPE1":  "RA---TAN",
        "CTYPE2":  "DEC--TAN",
        "CRVAL1":   0.0, "CRVAL2": 0.0,
        "CRPIX1":   1.0, "CRPIX2": 1.0,
    })
    img = bhio.load_image(path)
    cal = calibrate(img)
    assert cal.survey == "VLA-FIRST"
    assert cal.band == "radio"
    assert cal.array.unit == u.Jy
    pixel_sr = (pix_deg ** 2) * (math.pi / 180.0) ** 2
    beam_sr = math.pi / (4 * math.log(2)) * (math.radians(bmaj_deg)) ** 2
    factor = pixel_sr / beam_sr
    assert cal.array[0, 0].value == pytest.approx(0.01 * factor, rel=1e-6)


def test_vla_first_requires_bmaj_bmin(tmp_path: Path) -> None:
    path = _write_image(tmp_path, "first_nobeam.fits", np.zeros((2, 2)), {
        "TELESCOP": "VLA", "SURVEY": "VLA FIRST 1.4 GHz",
        "BUNIT": "Jy/beam", **_common_cdelt(),
    })
    img = bhio.load_image(path)
    with pytest.raises(UncalibratedDataError, match="BMAJ"):
        calibrate(img)


# ---------------------------------------------------------------------------
# RASS
# ---------------------------------------------------------------------------


def test_rass_calibration_matches_hand_computed(tmp_path: Path) -> None:
    """data (counts) / EXPTIME * ECF = erg/s/cm^2/pixel."""
    data = np.array([[10.0, 20.0]])
    exptime = 400.0  # s
    path = _write_image(tmp_path, "rass.fits", data, {
        "TELESCOP": "ROSAT", "SURVEY": "RASS BROAD",
        "EXPTIME": exptime, **_common_cdelt(),
    })
    img = bhio.load_image(path)
    cal = calibrate(img)
    expected_0 = (10.0 / exptime) * RASS_ECF_BROAD
    assert cal.survey == "RASS-broad"
    assert cal.band == "xray"
    assert cal.array.unit == u.erg / u.s / u.cm**2
    assert cal.array[0, 0].value == pytest.approx(expected_0, rel=1e-9)


def test_rass_rejects_zero_exptime(tmp_path: Path) -> None:
    path = _write_image(tmp_path, "rass_zero.fits", np.zeros((2, 2)), {
        "TELESCOP": "ROSAT", "SURVEY": "RASS BROAD",
        "EXPTIME": 0.0, **_common_cdelt(),
    })
    img = bhio.load_image(path)
    with pytest.raises(UncalibratedDataError, match="EXPTIME must be positive"):
        calibrate(img)


# ---------------------------------------------------------------------------
# DSS — must raise
# ---------------------------------------------------------------------------


def test_dss_raises_uncalibrated(tmp_path: Path) -> None:
    path = _write_image(tmp_path, "dss.fits", np.zeros((2, 2)), {
        "TELESCOP": "Palomar", "SURVEY": "Digitized Sky Survey (DSS2)",
        **_common_cdelt(),
    })
    img = bhio.load_image(path)
    with pytest.raises(UncalibratedDataError, match="DSS"):
        calibrate(img)


# ---------------------------------------------------------------------------
# Dispatch: unknown survey raises
# ---------------------------------------------------------------------------


def test_unknown_survey_raises(tmp_path: Path) -> None:
    path = _write_image(tmp_path, "alien.fits", np.zeros((2, 2)), {
        "TELESCOP": "ALIEN-OBS-2099",
        "SURVEY":   "Unknown Survey",
        **_common_cdelt(),
    })
    img = bhio.load_image(path)
    with pytest.raises(UncalibratedDataError, match="Could not identify"):
        calibrate(img)


def test_detect_survey_uses_filename_fallback(tmp_path: Path) -> None:
    """If the SURVEY/TELESCOP keys are absent but the filename mentions
    a survey, detect_survey should still pick the right one."""
    path = _write_image(tmp_path, "obj_2mass_K.fits", np.zeros((2, 2)), {
        "MAGZP": 21.0, **_common_cdelt(),
    })
    img = bhio.load_image(path)
    assert detect_survey(img) == "2MASS-K"


# ---------------------------------------------------------------------------
# All calibrators are in the registry and callable
# ---------------------------------------------------------------------------


def test_registry_completeness() -> None:
    for key, fn in CALIBRATORS.items():
        assert callable(fn), key


def test_pixel_solid_angle_uses_cd_matrix(tmp_path: Path) -> None:
    """When CDELT keys are absent but CD1_1/CD2_2 are present, scale
    still resolves. Uses the IRIS calibrator path."""
    path = _write_image(tmp_path, "iris_cdmatrix.fits",
                        np.array([[1.0, 1.0]]),
                        {
                            "TELESCOP": "IRAS", "SURVEY": "IRIS 12",
                            "BUNIT": "MJy/sr",
                            "CD1_1": -1.0 / 3600.0, "CD2_2": 1.0 / 3600.0,
                            "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
                            "CRVAL1": 0.0, "CRVAL2": 0.0,
                            "CRPIX1": 1.0, "CRPIX2": 1.0,
                        })
    img = bhio.load_image(path)
    cal = calibrate(img)
    assert cal.array.unit == u.Jy
    # 1 arcsec/pixel → ~2.35e-5 Jy/pixel for 1 MJy/sr.
    sr = (1.0 / 3600.0) ** 2 * (math.pi / 180.0) ** 2
    assert cal.array[0, 0].value == pytest.approx(sr * 1e6, rel=1e-6)


def test_pixel_solid_angle_raises_with_no_scale(tmp_path: Path) -> None:
    path = _write_image(tmp_path, "iris_no_scale.fits",
                        np.array([[1.0, 1.0]]),
                        {"TELESCOP": "IRAS", "SURVEY": "IRIS 12",
                         "BUNIT": "MJy/sr"})
    img = bhio.load_image(path)
    with pytest.raises(UncalibratedDataError, match="pixel solid angle"):
        calibrate(img)
