"""
Variability metrics — fractional rms, excess variance.

FRACTIONAL RMS VARIABILITY AMPLITUDE F_var
==========================================
A model-independent measure of how variable a light curve is, accounting
for Poisson noise:

  σ_xs²   = σ²_obs - <σ_err²>             (excess variance)
  F_var   = sqrt(σ_xs² / <x>²)             (rms fractional variability)

If F_var = 0.1, the source varies by 10% rms above what counting noise
would produce. AGN typically show 5-30% F_var on hour-day timescales;
X-ray binaries can hit ~50% during bright hard states.

REFERENCES
==========
- Vaughan et al. 2003, MNRAS 345, 1271 — the definitive treatment of F_var
  uncertainties.
- Edelson et al. 2002, ApJ 568, 610 — earlier excess-variance work.
"""

from __future__ import annotations

import numpy as np


def excess_variance(rates: np.ndarray, errors: np.ndarray) -> float:
    """σ²_xs = σ²(rates) - <σ²(errors)>. Can be negative if data is
    pure noise."""
    if rates.size < 2:
        return float("nan")
    sigma_obs_sq = float(np.var(rates, ddof=1))
    mean_err_sq = float(np.mean(errors**2))
    return sigma_obs_sq - mean_err_sq


def fractional_rms(rates: np.ndarray, errors: np.ndarray) -> float:
    """F_var = sqrt(σ²_xs) / <rates>. Returns NaN if the excess variance is
    negative (i.e. consistent with no intrinsic variability)."""
    if rates.size < 2:
        return float("nan")
    mean_rate = float(np.mean(rates))
    if mean_rate <= 0:
        return float("nan")
    xs = excess_variance(rates, errors)
    if xs <= 0:
        return float("nan")
    return float(np.sqrt(xs) / mean_rate)


def fractional_rms_error(rates: np.ndarray, errors: np.ndarray) -> float:
    """Approximate 1-σ error on F_var (Vaughan et al. 2003 eq. B2).

    err(F_var) ≈ sqrt[ (1/2N) · <σ_err²>² / (F_var · <x>²)²
                     + (1/N) · <σ_err²> / <x>² ]^(1/2)
    """
    F = fractional_rms(rates, errors)
    if not np.isfinite(F):
        return float("nan")
    N = rates.size
    mean_rate = float(np.mean(rates))
    mean_err_sq = float(np.mean(errors**2))
    term1 = (1 / (2 * N)) * (mean_err_sq**2) / ((F * mean_rate**2)**2)
    term2 = (1 / N) * mean_err_sq / mean_rate**2
    return float(np.sqrt(term1 + term2))
