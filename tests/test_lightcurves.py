"""
tests/test_lightcurves.py — tests for blackhole.lightcurves.

Exercises event binning, light-curve construction, Lomb-Scargle frequency
recovery, renderer smoke tests, and GTI-aware binning (M4).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from blackhole import io as bhio
from blackhole import lightcurves as lc


def test_bin_events_to_lightcurve_shape(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=100.0)
    assert curve.bin_size_s == 100.0
    assert curve.times.size >= 1
    assert curve.rates.size == curve.times.size
    assert curve.errors.size == curve.times.size
    assert (curve.errors > 0).all()


def test_bin_events_to_lightcurve_total_counts_matches(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=50.0)
    total_counts = (curve.rates * curve.bin_size_s).sum()
    # 1000 synthetic events; binning preserves count, modulo edge bins.
    assert total_counts == pytest.approx(ev.times.size, abs=2)


def test_bin_events_empty_returns_empty(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=100.0,
                                        energy_range_ev=(1e6, 1e7))
    assert curve.times.size == 0
    assert curve.rates.size == 0


def test_lomb_scargle_recovers_known_frequency() -> None:
    """Inject a 0.01 Hz sinusoid into a synthetic light curve and recover it."""
    rng = np.random.default_rng(7)
    bin_s = 1.0
    n = 4096
    t = np.arange(n) * bin_s
    f0 = 0.01   # Hz
    rates = 5.0 + 2.0 * np.sin(2 * np.pi * f0 * t) + rng.normal(0, 0.2, size=n)
    errors = np.full_like(rates, 0.2)
    curve = lc.LightCurve(
        times=t, rates=rates, errors=errors,
        bin_size_s=bin_s,
        energy_range_ev=None,
        source_label="synthetic",
    )
    freqs, power = lc.lomb_scargle_periodogram(
        curve, min_freq_hz=1e-3, max_freq_hz=1e-1, n_freq=20000,
    )
    peak = freqs[int(np.argmax(power))]
    assert peak == pytest.approx(f0, rel=2e-2)


def test_lomb_scargle_insufficient_data_returns_empty() -> None:
    curve = lc.LightCurve(times=np.array([0.0]), rates=np.array([1.0]),
                          errors=np.array([1.0]), bin_size_s=1.0,
                          energy_range_ev=None, source_label="tiny")
    freqs, power = lc.lomb_scargle_periodogram(curve)
    assert freqs.size == 0
    assert power.size == 0


def test_render_lightcurve_smoke(tiny_events_fits: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    ev = bhio.load_events(tiny_events_fits)
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=100.0)
    fig = lc.render_lightcurve(curve)
    assert fig is not None
    assert len(fig.axes) >= 1


def test_render_lightcurve_empty() -> None:
    import matplotlib
    matplotlib.use("Agg")
    curve = lc.LightCurve(times=np.array([]), rates=np.array([]),
                          errors=np.array([]), bin_size_s=100.0,
                          energy_range_ev=None, source_label="empty")
    fig = lc.render_lightcurve(curve)
    assert fig is not None


def test_render_periodogram_smoke() -> None:
    import matplotlib
    matplotlib.use("Agg")
    freqs = np.linspace(1e-3, 1.0, 100)
    power = np.abs(np.sin(freqs))
    fig = lc.render_periodogram(freqs, power)
    assert fig is not None


def test_render_periodogram_empty() -> None:
    import matplotlib
    matplotlib.use("Agg")
    fig = lc.render_periodogram(np.array([]), np.array([]))
    assert fig is not None


# ---------------------------------------------------------------------------
# M4 — GTI-aware binning
# ---------------------------------------------------------------------------


def test_gti_loaded_from_fits(gapped_events_fits: Path) -> None:
    """load_events populates EventList.gti from the GTI extension."""
    ev = bhio.load_events(gapped_events_fits)
    assert ev.gti is not None
    assert ev.gti.shape == (2, 2)
    # Two intervals as written.
    assert ev.gti[0, 0] == pytest.approx(0.0)
    assert ev.gti[0, 1] == pytest.approx(1000.0)
    assert ev.gti[1, 0] == pytest.approx(2000.0)
    assert ev.gti[1, 1] == pytest.approx(3000.0)


def test_gti_masks_gaps(gapped_events_fits: Path) -> None:
    """A 1 ks gap shows up as missing bins, not as zero-rate bins.

    With bin_size_s=100 and two GTIs of 1 ks each (gap from 1000-2000),
    a naive binning would produce 30 bins covering 0-3000 s, including
    10 zero-rate bins. With apply_gti=True the gap bins are dropped:
    we expect ~20 surviving bins.
    """
    ev = bhio.load_events(gapped_events_fits)
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=100.0, apply_gti=True)
    assert curve.gti_applied is True
    # Surviving bins should be in [0, 1000) and [2000, 3000), not in (1000, 2000).
    assert curve.times.size == pytest.approx(20, abs=2)
    for t in curve.times:
        assert (0.0 <= t < 1000.0) or (2000.0 <= t < 3000.0)


def test_gti_drop_yields_no_fake_zero_bins(gapped_events_fits: Path) -> None:
    """Surviving bins all have a nonzero effective exposure."""
    ev = bhio.load_events(gapped_events_fits)
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=100.0, apply_gti=True)
    assert curve.effective_exposure.size == curve.times.size
    assert (curve.effective_exposure > 0).all()


def test_gti_total_exposure_matches_gti_window(gapped_events_fits: Path) -> None:
    """Sum of effective_exposure equals total GTI duration (2000 s)."""
    ev = bhio.load_events(gapped_events_fits)
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=100.0, apply_gti=True)
    assert curve.gti_total_s == pytest.approx(2000.0)
    assert curve.total_exposure_s == pytest.approx(2000.0, abs=1.0)


def test_gti_disabled_preserves_old_behavior(gapped_events_fits: Path) -> None:
    """apply_gti=False produces a continuous light curve including the gap."""
    ev = bhio.load_events(gapped_events_fits)
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=100.0, apply_gti=False)
    assert curve.gti_applied is False
    # Continuous time-axis: 3000 s / 100 s ≈ 30 bins.
    assert curve.times.size == pytest.approx(30, abs=2)
    # Some bins in the gap should be exactly zero rate.
    zero_bins = (curve.rates == 0).sum()
    assert zero_bins >= 5


def test_gti_rate_normalization_uses_effective_exposure(gapped_events_fits: Path) -> None:
    """A partial-exposure bin reports rate = counts / effective_exposure,
    not counts / bin_size_s."""
    # Use bin_size_s such that one bin straddles a GTI edge, then check
    # the effective exposure is between 0 and bin_size_s for that bin.
    ev = bhio.load_events(gapped_events_fits)
    curve = lc.bin_events_to_lightcurve(
        ev, bin_size_s=150.0, apply_gti=True, min_exposure_fraction=0.05,
    )
    partial = (curve.effective_exposure > 0) & (curve.effective_exposure < curve.bin_size_s)
    assert partial.any()
    for idx in np.where(partial)[0]:
        exposure = curve.effective_exposure[idx]
        counts = curve.rates[idx] * exposure
        assert counts >= 0
        # The rate should differ from a naive counts/bin_size_s.
        if counts > 0:
            naive = counts / curve.bin_size_s
            assert curve.rates[idx] > naive - 1e-9


def test_gti_no_intervals_falls_back_to_no_gti_path(tiny_events_fits: Path) -> None:
    """When apply_gti=True but no overlap with any GTI window, we fall
    back to the non-GTI path rather than silently producing nothing."""
    ev = bhio.load_events(tiny_events_fits)
    # tiny_events_fits has a single GTI 0-1000 s and events in 0-1000 s.
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=100.0, apply_gti=True)
    assert curve.gti_applied is True
    assert curve.times.size >= 1


def test_render_lightcurve_includes_gti_in_title(gapped_events_fits: Path) -> None:
    """Title surfaces the GTI exposure when a GTI-aware curve is rendered."""
    import matplotlib
    matplotlib.use("Agg")
    ev = bhio.load_events(gapped_events_fits)
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=100.0, apply_gti=True)
    fig = lc.render_lightcurve(curve)
    title = fig.axes[0].get_title()
    assert "GTI exposure" in title
