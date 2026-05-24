"""
Physics tests for blackhole.physics.infrared — blackbody, Stern+2012 cut,
Donley+2012 wedge.

References
----------
- Stern et al. 2012, ApJ 753, 30 — W1-W2 > 0.8 (Vega) AGN cut.
- Donley et al. 2012, ApJ 748, 142 — refined four-band WISE wedge.
"""

from __future__ import annotations

import numpy as np

from blackhole.physics import infrared as ir


def test_blackbody_nu_fnu_positive_for_realistic_temperatures() -> None:
    # 1 micron in Hz -> ~3e14 Hz, T=5000K is a reasonable star.
    nu = np.array([3e14, 1e15])
    out = ir.blackbody_nu_fnu(nu, temperature_k=5000.0)
    assert np.all(out > 0)
    assert out.shape == nu.shape


def test_blackbody_nu_fnu_higher_temperature_shifts_peak_to_higher_nu() -> None:
    """Wien's law: λ_max T = const. Hotter -> higher nu."""
    nu = np.logspace(12, 17, 1000)
    bb_cool = ir.blackbody_nu_fnu(nu, temperature_k=300.0)
    bb_hot = ir.blackbody_nu_fnu(nu, temperature_k=10000.0)
    nu_peak_cool = nu[int(np.argmax(bb_cool))]
    nu_peak_hot = nu[int(np.argmax(bb_hot))]
    assert nu_peak_hot > nu_peak_cool


def test_blackbody_nu_fnu_scalar_input() -> None:
    val = ir.blackbody_nu_fnu(1e14, temperature_k=1000.0)
    assert val > 0


def test_stern_2012_cut_above_threshold() -> None:
    assert bool(ir.stern_2012_agn(1.0)) is True
    assert bool(ir.stern_2012_agn(0.5)) is False


def test_stern_2012_at_threshold() -> None:
    # Spec says > 0.8 strictly.
    assert bool(ir.stern_2012_agn(0.8)) is False
    assert bool(ir.stern_2012_agn(0.81)) is True


def test_stern_2012_vectorized() -> None:
    arr = np.array([0.3, 0.8, 1.0, 1.5])
    result = ir.stern_2012_agn(arr)
    np.testing.assert_array_equal(result, np.array([False, False, True, True]))


def test_donley_2012_classic_agn_corner() -> None:
    """A clear-AGN colour pair from the Donley+2012 Fig. 1 wedge interior."""
    # W2-W3 = 1.5, W1-W2 = 1.0 sits inside the wedge.
    assert bool(ir.donley_2012_agn(w1_w2_vega=1.0, w2_w3_vega=1.5)) is True


def test_donley_2012_rejects_blue_stellar_colours() -> None:
    """Main-sequence stars cluster near W1-W2 ~ 0, W2-W3 ~ 0; should not be AGN."""
    assert bool(ir.donley_2012_agn(w1_w2_vega=0.0, w2_w3_vega=0.0)) is False


def test_donley_2012_rejects_red_galaxies_outside_wedge() -> None:
    """Star-forming galaxies often sit at large W2-W3 but small W1-W2 — outside."""
    # W2-W3 = 3.0 > 2.2 upper limit -> reject.
    assert bool(ir.donley_2012_agn(w1_w2_vega=1.0, w2_w3_vega=3.0)) is False


def test_donley_2012_vectorized_input() -> None:
    c12 = np.array([0.0, 1.0, 1.0])
    c23 = np.array([0.0, 1.5, 3.0])
    out = ir.donley_2012_agn(c12, c23)
    assert out.shape == (3,)
    assert not bool(out[0])
    assert bool(out[1])
    assert not bool(out[2])


def test_donley_lower_left_boundary() -> None:
    """W1-W2 = 0.8 - 0.1*(W2-W3): at W2-W3 = 0, threshold is 0.8."""
    # Just above the lower-left boundary at W2-W3 = 1.0 → threshold 0.7.
    assert bool(ir.donley_2012_agn(w1_w2_vega=0.71, w2_w3_vega=1.0)) is False  # upper cut
    assert bool(ir.donley_2012_agn(w1_w2_vega=0.71, w2_w3_vega=1.3)) is True
