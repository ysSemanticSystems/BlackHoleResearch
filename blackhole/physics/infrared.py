"""
Infrared diagnostics — dust torus blackbody, WISE color cuts.

WHY MID-IR IDENTIFIES AGN
=========================
In the AGN unification picture (Antonucci 1993, Urry & Padovani 1995),
the central engine is surrounded by a dusty torus on parsec scales. UV/
optical/X-ray light from the disk and corona heats the torus to T ~ few
hundred to ~1500 K (above which dust sublimates). The torus re-radiates
as a thermal IR bump peaking at 10-30 µm — observable by WISE, Spitzer,
JWST/MIRI.

This is the cleanest single-band AGN selection because the torus signal:
  - Is roughly isotropic (unlike the obscured optical/X-ray light)
  - Has no real galaxy/stellar contaminant at those temperatures
  - Persists for obscured (type-2) AGN that look quiet in optical

WISE COLOR CUTS — STERN AND DONLEY WEDGES
==========================================
WISE bands: W1 (3.4 µm), W2 (4.6 µm), W3 (12 µm), W4 (22 µm).
On a W1-W2 vs W2-W3 color-color diagram, AGN occupy a distinct region.
- Stern et al. 2012, ApJ 753, 30: simple W1-W2 > 0.8 cut.
- Donley et al. 2012, ApJ 748, 142: refined wedge using all 4 WISE bands,
  much higher purity but lower completeness.

PITFALL: WISE colors are in *Vega magnitudes* by convention (not AB). The
zero-point conversions (Wright et al. 2010, AJ 140, 1868):
  W1: AB - Vega = 2.699
  W2: AB - Vega = 3.339
  W3: AB - Vega = 5.174
  W4: AB - Vega = 6.620

LITTLE RED DOTS BREAK THIS TOO
==============================
LRDs in JWST/COSMOS-Web are NOT detected in mid-IR even when stacked
(Akins et al. 2025, arXiv 2406.10341). They violate Stern/Donley because
they lack the hot dust torus. Why is still debated — possibilities
include intrinsically dust-poor accretion, super-Eddington outflows
disrupting the torus, or these objects being something other than
classical AGN entirely.
"""

from __future__ import annotations

import astropy.constants as const
import numpy as np


def blackbody_nu_fnu(
    frequency_hz: np.ndarray | float, temperature_k: float,
) -> np.ndarray | float:
    """νB_ν(T) in arbitrary units (for SED overplotting).

    B_ν = (2 h ν³ / c²) · 1/(exp(hν/kT) - 1)
    """
    nu = np.asarray(frequency_hz, dtype=float)
    h = const.h.cgs.value
    k = const.k_B.cgs.value
    c = const.c.cgs.value
    with np.errstate(over="ignore"):
        Bnu = (2 * h * nu**3 / c**2) / np.expm1(h * nu / (k * temperature_k))
    result: np.ndarray | float = nu * Bnu
    return result


def stern_2012_agn(w1_w2_vega: float | np.ndarray) -> bool | np.ndarray:
    """Stern et al. 2012 single-color WISE AGN cut.

    Returns True for sources with W1-W2 > 0.8 (Vega), the canonical AGN
    color selection.
    """
    return np.asarray(w1_w2_vega) > 0.8


def donley_2012_agn(
    w1_w2_vega: float | np.ndarray, w2_w3_vega: float | np.ndarray,
) -> bool | np.ndarray:
    """Donley et al. 2012 four-band WISE AGN wedge (Vega magnitudes).

    Applies all four cuts that define the AGN wedge:
      W1 - W2 > 0.8 - 0.1·(W2 - W3)     (lower-left boundary)
      W1 - W2 > -0.7                      (faint AGN cutoff)
      W1 - W2 < 1.7·(W2 - W3) - 1.4       (upper-left)
      W2 - W3 < 2.2                       (upper-right boundary)

    Higher purity than Stern at the cost of completeness — recommended
    for catalog AGN selection.
    """
    c12 = np.asarray(w1_w2_vega)
    c23 = np.asarray(w2_w3_vega)
    return (
        (c12 > (0.8 - 0.1 * c23)) &
        (c12 > -0.7) &
        (c12 < (1.7 * c23 - 1.4)) &
        (c23 < 2.2)
    )
