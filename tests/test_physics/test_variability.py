"""
Physics tests for blackhole.physics.variability — F_var, excess variance.

References
----------
- Vaughan et al. 2003, MNRAS 345, 1271 — definitive treatment.
"""

from __future__ import annotations

import numpy as np
import pytest

from blackhole.physics import variability as v


def test_excess_variance_pure_noise_is_near_zero() -> None:
    """A flat signal with the right per-bin error should have ~zero excess variance."""
    rng = np.random.default_rng(0)
    rates = 10.0 + rng.normal(0, 1.0, size=1024)
    errors = np.full_like(rates, 1.0)
    xs = v.excess_variance(rates, errors)
    # Sample variance ~ 1.0, mean error² ~ 1.0; the difference is small relative to 1.0.
    assert abs(xs) < 0.2


def test_fractional_rms_on_pure_noise_returns_nan() -> None:
    rng = np.random.default_rng(1)
    rates = 10.0 + rng.normal(0, 1.0, size=1024)
    errors = np.full_like(rates, 1.0)
    fvar = v.fractional_rms(rates, errors)
    # Excess variance ≤ 0 → NaN by convention.
    assert np.isnan(fvar) or fvar < 0.05


def test_fractional_rms_on_known_sinusoid_recovers_amplitude() -> None:
    """A pure sinusoid of amplitude A around mean μ has σ_rms = A/√2.

    F_var should ≈ (A/√2) / μ when noise is much smaller than A.
    """
    n = 4096
    mu = 10.0
    A = 2.0
    t = np.arange(n)
    rates = mu + A * np.sin(2 * np.pi * t / 100.0)
    errors = np.full_like(rates, 0.01)   # tiny errors so excess variance ≈ var(rates)
    fvar = v.fractional_rms(rates, errors)
    expected = (A / np.sqrt(2)) / mu
    assert fvar == pytest.approx(expected, rel=2e-2)


def test_fractional_rms_error_is_positive_when_fvar_finite() -> None:
    rng = np.random.default_rng(2)
    rates = 10.0 + 2.0 * np.sin(np.linspace(0, 20, 500))
    errors = rng.uniform(0.05, 0.1, size=500)
    fvar = v.fractional_rms(rates, errors)
    err = v.fractional_rms_error(rates, errors)
    if np.isfinite(fvar):
        assert err > 0


def test_fractional_rms_zero_mean_returns_nan() -> None:
    rates = np.zeros(100)
    errors = np.ones(100)
    assert np.isnan(v.fractional_rms(rates, errors))


def test_fractional_rms_short_input_returns_nan() -> None:
    assert np.isnan(v.fractional_rms(np.array([1.0]), np.array([0.1])))
    assert np.isnan(v.excess_variance(np.array([1.0]), np.array([0.1])))


def test_fractional_rms_error_nan_when_fvar_nan() -> None:
    rates = np.zeros(10)
    errors = np.ones(10)
    assert np.isnan(v.fractional_rms_error(rates, errors))
