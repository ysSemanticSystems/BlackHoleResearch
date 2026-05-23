"""
spectral_xray — X-ray spectral diagnostics.

HARDNESS RATIO
==============
The simplest model-independent spectral diagnostic. Counts in a hard band
divided by counts in a soft band, normalized:

  HR = (H - S) / (H + S)

  HR ranges in [-1, +1]. Positive = hard (Compton-thick AGN, accretion-disk
  corona dominated). Negative = soft (unobscured AGN, thermal disk).

PHOTON INDEX Γ
==============
For a power-law spectrum N(E) ∝ E^(-Γ):

  Γ ≈ 1.7-2.1   Typical type-1 (unobscured) AGN
  Γ ≈ 1.4-1.7   Obscured / type-2 AGN with hard reflection-dominated spectrum
  Γ > 2.5       Very steep — sometimes ultra-soft tidal disruption events
  Γ < 1.4       Reflection-dominated; absorbed continuum

REFERENCES
==========
- Ricci et al. 2017, ApJS 233, 17 — Swift/BAT AGN survey, full population
  distributions of Γ and HR.
- Hickox & Alexander 2018, ARA&A 56, 625 — review of obscured AGN selection.
"""

from __future__ import annotations

import numpy as np


def hardness_ratio(
    soft_counts: float | np.ndarray,
    hard_counts: float | np.ndarray,
) -> float | np.ndarray:
    """HR = (H - S) / (H + S). Returns NaN where total is zero."""
    s = np.asarray(soft_counts, dtype=float)
    h = np.asarray(hard_counts, dtype=float)
    total = h + s
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(total > 0, (h - s) / total, np.nan)


def hardness_ratio_error(
    soft_counts: float | np.ndarray,
    hard_counts: float | np.ndarray,
) -> float | np.ndarray:
    """1-sigma error on HR via Poisson error propagation.

    σ_HR = 2 sqrt(S²·H + H²·S) / (H + S)²
    Assumes independent Poisson statistics in each band.
    """
    s = np.asarray(soft_counts, dtype=float)
    h = np.asarray(hard_counts, dtype=float)
    total = h + s
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(
            total > 0,
            2 * np.sqrt(s**2 * h + h**2 * s) / total**2,
            np.nan,
        )


def classify_photon_index(gamma: float) -> str:
    """Heuristic verbal label for a photon index value. Population ranges
    from Ricci+2017 and related Swift/BAT survey work."""
    if gamma < 1.0:
        return "unphysically hard — check fit"
    if gamma < 1.4:
        return "very hard — reflection-dominated or Compton-thick"
    if gamma < 1.7:
        return "hard — likely obscured AGN"
    if gamma < 2.1:
        return "typical type-1 AGN range"
    if gamma < 2.5:
        return "soft — possible thermal contribution"
    return "very soft — TDE-like or super-Eddington"
