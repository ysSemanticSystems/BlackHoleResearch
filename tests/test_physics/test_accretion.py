"""
Physics regression tests for blackhole.physics.accretion.

These reproduce textbook values to specified precision and are intentionally
strict — a drift in any of these numbers is a sign of a physics-affecting
change and should be visible in CI.

References
----------
- Frank, King & Raine 2002, "Accretion Power in Astrophysics", 3rd ed.
  Cambridge University Press. Eq. (1.5): L_Edd ~ 1.26e38 erg/s for 1 M_sun.
- Shakura & Sunyaev 1973, A&A 24, 337.
"""

from __future__ import annotations

import astropy.units as u
import pytest

from blackhole.physics import accretion as a


def test_eddington_luminosity_solar_textbook() -> None:
    """L_Edd(1 M_sun) ≈ 1.26e38 erg/s (Frank+2002 Eq. 1.5)."""
    L = a.eddington_luminosity(1.0)
    assert L.unit == u.erg / u.s
    assert L.value == pytest.approx(1.26e38, rel=5e-3)


def test_eddington_luminosity_scales_linearly_with_mass() -> None:
    L1 = a.eddington_luminosity(1.0).value
    L1e8 = a.eddington_luminosity(1e8).value
    assert L1e8 / L1 == pytest.approx(1e8, rel=1e-10)


def test_eddington_ratio_at_limit_is_one() -> None:
    L_Edd = a.eddington_luminosity(1.0).value
    lambda_Edd = a.eddington_ratio(L_Edd, 1.0)
    assert lambda_Edd == pytest.approx(1.0, rel=1e-12)


def test_eddington_ratio_at_half_is_one_half() -> None:
    L_Edd = a.eddington_luminosity(1e7).value
    half = a.eddington_ratio(0.5 * L_Edd, 1e7)
    assert half == pytest.approx(0.5, rel=1e-12)


def test_accretion_rate_at_eddington_for_1e8_msun() -> None:
    """At L = L_Edd with η = 0.1, Ṁ ~ 2.21 M_sun/yr for M = 1e8 M_sun.

    L_Edd(1e8 M_sun) = 1.26e46 erg/s; Ṁ = L / (0.1 * c²) ~ 1.40e26 g/s
    ~ 2.22 M_sun/yr.
    """
    M = 1e8
    L_Edd = a.eddington_luminosity(M).value
    Mdot = a.accretion_rate_msun_yr(L_Edd, efficiency=0.1)
    assert Mdot == pytest.approx(2.22, rel=5e-2)


def test_accretion_rate_higher_efficiency_lower_mdot() -> None:
    L = 1e44
    m1 = a.accretion_rate_msun_yr(L, efficiency=0.1)
    m2 = a.accretion_rate_msun_yr(L, efficiency=0.4)
    assert m2 < m1
    # 4x efficiency -> 4x less mass to feed.
    assert m1 / m2 == pytest.approx(4.0, rel=1e-12)


def test_shakura_sunyaev_temperature_decreases_with_radius() -> None:
    """T(r) ∝ r^(-3/4) far from inner boundary."""
    M = 1e7
    lam = 0.1
    T_inner = a.shakura_sunyaev_temperature(20.0, M, lam).value
    T_outer = a.shakura_sunyaev_temperature(200.0, M, lam).value
    assert T_inner > T_outer
    # Ratio of T at r1 vs r2 should follow (r1/r2)^(-3/4).
    ratio_expected = (20.0 / 200.0) ** (-0.75)
    ratio_actual = T_inner / T_outer
    assert ratio_actual == pytest.approx(ratio_expected, rel=1e-6)


def test_shakura_sunyaev_temperature_increases_with_eddington_ratio() -> None:
    """At fixed r and M, T ∝ Ṁ^(1/4) ∝ λ_Edd^(1/4)."""
    M = 1e7
    r = 50.0
    T_low = a.shakura_sunyaev_temperature(r, M, 0.01).value
    T_high = a.shakura_sunyaev_temperature(r, M, 1.0).value
    assert T_high > T_low
    # Ratio should match (1.0/0.01)^(1/4) = 100^0.25 ≈ 3.162.
    assert (T_high / T_low) == pytest.approx(100 ** 0.25, rel=1e-6)


def test_shakura_sunyaev_temperature_units_are_kelvin() -> None:
    T = a.shakura_sunyaev_temperature(20.0, 1e7, 0.1)
    assert T.unit == u.K
    assert T.value > 0.0
