"""
Accretion physics — the formulas that turn a BH mass into observable
luminosities and temperatures.

DERIVATIONS AND SOURCES
=======================
Eddington luminosity: L_Edd = 4π G M m_p c / σ_T
  - Sets the maximum steady-state luminosity from spherical accretion
    where radiation pressure on free electrons balances gravity on protons.
  - Numerically: L_Edd ≈ 1.26 × 10^38 (M/M_sun) erg/s
  - Eddington 1926, "The Internal Constitution of the Stars"
  - Caveats: applies to fully ionized hydrogen; assumes spherical symmetry;
    real accretion disks can briefly exceed L_Edd ("super-Eddington" — see
    Little Red Dot models, e.g. Pacucci & Narayan 2025).

Accretion rate from luminosity: Ṁ = L / (η c²)
  - η = radiative efficiency. For a Schwarzschild BH η ≈ 0.057.
    For a maximally spinning Kerr BH η ≈ 0.42. Standard "thin disk"
    assumption: η ≈ 0.1.
  - Novikov & Thorne 1973, "Black Holes (Les Astres Occlus)" eds.
    DeWitt & DeWitt — original GR thin disk derivation.

Shakura-Sunyaev disk temperature: T(r) ∝ (M Ṁ)^(1/4) r^(-3/4)
  - Shakura & Sunyaev 1973, A&A 24, 337. Defines the "α-disk" framework
    that still underpins almost all accretion disk modeling.
"""

from __future__ import annotations

import astropy.constants as const
import astropy.units as u
import numpy as np


def eddington_luminosity(mass_msun: float) -> u.Quantity:
    """L_Edd in erg/s for a black hole of mass M (in solar masses).

    L_Edd = 4π G M m_p c / σ_T  (for ionized hydrogen)
    """
    M = mass_msun * const.M_sun
    L = (4 * np.pi * const.G * M * const.m_p * const.c / const.sigma_T)
    return L.to(u.erg / u.s)


def eddington_ratio(luminosity_erg_s: float, mass_msun: float) -> float:
    """λ_Edd = L_bol / L_Edd. Dimensionless. Order 1 means radiation-pressure
    limit; > 1 indicates super-Eddington accretion (the LRD regime)."""
    return float(luminosity_erg_s / eddington_luminosity(mass_msun).value)


def accretion_rate_msun_yr(luminosity_erg_s: float, efficiency: float = 0.1) -> float:
    """Ṁ in M_sun/yr from luminosity, assuming efficiency η.

    L = η Ṁ c²  =>  Ṁ = L / (η c²)
    Default η = 0.1 is the standard thin-disk assumption.
    """
    L = luminosity_erg_s * u.erg / u.s
    Mdot = L / (efficiency * const.c**2)
    return float(Mdot.to(u.M_sun / u.yr).value)


def shakura_sunyaev_temperature(
    radius_rg: float, mass_msun: float, eddington_ratio_value: float,
    efficiency: float = 0.1,
) -> u.Quantity:
    """Disk temperature at radius r (in units of gravitational radii R_g = GM/c²).

    Standard Shakura-Sunyaev 1973 thin disk profile, ignoring the inner
    boundary correction factor [1 - (R_in/r)^(1/2)] for simplicity. For
    realistic spectra include that factor and integrate.

    Returns temperature in Kelvin.
    """
    # R_g in cm
    Rg_cm = (const.G * mass_msun * const.M_sun / const.c**2).to(u.cm).value
    r_cm = radius_rg * Rg_cm
    # Ṁ in g/s from the Eddington ratio
    L_bol = eddington_ratio_value * eddington_luminosity(mass_msun).value  # erg/s
    Mdot_gs = L_bol / (efficiency * const.c.cgs.value**2)
    # SS73 temperature: T(r) = [3 G M Ṁ / (8 π σ r³)]^(1/4)
    G_cgs = const.G.cgs.value
    M_cgs = mass_msun * const.M_sun.cgs.value
    sigma_cgs = const.sigma_sb.cgs.value
    T = (3 * G_cgs * M_cgs * Mdot_gs / (8 * np.pi * sigma_cgs * r_cm**3))**0.25
    return T * u.K
