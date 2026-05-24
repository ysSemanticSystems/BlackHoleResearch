"""
tests/test_spectra.py — seed tests for blackhole.spectra.

Exercises OGIP PHA loading and the descriptive channel-space power-law
fit on a synthetic spectrum with known input slope.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from blackhole import spectra as sp


def test_load_pha_spectrum_columns(tiny_pha_fits: Path) -> None:
    spec = sp.load_pha_spectrum(tiny_pha_fits)
    assert spec.channels.size == 256
    assert spec.counts.size == 256
    assert spec.errors is not None
    assert spec.exposure_s == pytest.approx(1000.0)
    assert spec.mission == "SYNTHETIC"
    assert spec.instrument == "TEST"


def test_load_pha_raises_without_spectrum_hdu(tmp_path: Path) -> None:
    from astropy.io import fits
    path = tmp_path / "empty.fits"
    fits.PrimaryHDU().writeto(path)
    with pytest.raises(ValueError, match="No SPECTRUM extension"):
        sp.load_pha_spectrum(path)


def test_fit_power_law_recovers_input_slope_within_uncertainty(
    tiny_pha_fits: Path,
) -> None:
    """Synthetic spectrum was generated with α=1.8 in channel space.

    The descriptive fit should land within ~3 sigma of that value.
    """
    spec = sp.load_pha_spectrum(tiny_pha_fits)
    gamma, gamma_err, fit_curve = sp.fit_power_law(spec)
    assert np.isfinite(gamma)
    assert np.isfinite(gamma_err)
    assert fit_curve.shape == spec.channels.shape
    # The slope is descriptive (channel-space), but the synthetic process is
    # known so we expect roughly 1.8 ± a few sigma.
    assert abs(gamma - 1.8) < 0.5, f"gamma={gamma}, err={gamma_err}"
    assert gamma_err > 0


def test_fit_power_law_with_channel_range(tiny_pha_fits: Path) -> None:
    spec = sp.load_pha_spectrum(tiny_pha_fits)
    g_full, _, _ = sp.fit_power_law(spec)
    g_sub, _, _ = sp.fit_power_law(spec, channel_range=(0, 128))
    assert np.isfinite(g_full)
    assert np.isfinite(g_sub)


def test_fit_power_law_raises_with_insufficient_data(tmp_path: Path) -> None:
    from astropy.io import fits
    cols = [
        fits.Column(name="CHANNEL",  array=np.array([0, 1, 2]), format="J"),
        fits.Column(name="COUNTS",   array=np.array([0, 0, 0]), format="J"),
        fits.Column(name="STAT_ERR", array=np.array([1.0, 1.0, 1.0]), format="E"),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name="SPECTRUM")
    path = tmp_path / "starved.pha"
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)
    spec = sp.load_pha_spectrum(path)
    with pytest.raises(ValueError, match="Not enough"):
        sp.fit_power_law(spec)


def test_render_spectrum_smoke(tiny_pha_fits: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    spec = sp.load_pha_spectrum(tiny_pha_fits)
    fit = sp.fit_power_law(spec)
    fig = sp.render_spectrum(spec, fit=fit)
    assert fig is not None
    assert len(fig.axes) >= 1


def test_render_spectrum_without_fit(tiny_pha_fits: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    spec = sp.load_pha_spectrum(tiny_pha_fits)
    fig = sp.render_spectrum(spec, fit=None)
    assert fig is not None


def test_render_spectrum_without_errors(tmp_path: Path) -> None:
    """When STAT_ERR is missing entirely, render should still produce a figure."""
    import matplotlib
    matplotlib.use("Agg")
    from astropy.io import fits

    channels = np.arange(64, dtype=np.int32)
    counts = np.round(1000.0 * np.power(channels + 1, -1.8)).astype(np.int32)
    cols = [
        fits.Column(name="CHANNEL", array=channels, format="J"),
        fits.Column(name="RATE",    array=counts.astype(np.float32),
                    format="E"),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name="SPECTRUM")
    path = tmp_path / "noerr.pha"
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)

    spec = sp.load_pha_spectrum(path)
    fig = sp.render_spectrum(spec, fit=None)
    assert fig is not None
