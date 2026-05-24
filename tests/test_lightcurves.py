"""
tests/test_lightcurves.py — seed tests for blackhole.lightcurves.

Exercises event binning, light-curve construction, Lomb-Scargle frequency
recovery, and renderer smoke tests. GTI-aware behaviour is deferred to
M4 (PHASE2_PLAN.md).
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
