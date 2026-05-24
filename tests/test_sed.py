"""
tests/test_sed.py — seed tests for blackhole.sed.

Validates the unit-safe conversions that the SED layer depends on: frequency
resolution from wavelength/energy, and nu_f_nu conversion from flux_density
in Jy. Hand-computed expected values from astropy constants.
"""

from __future__ import annotations

import astropy.constants as const
import astropy.units as u
import numpy as np
import pytest

from blackhole import sed as sedmod


def test_to_frequency_from_wavelength() -> None:
    p = sedmod.SEDPoint(label="2MASS K", wavelength=2.2 * u.micron, band="ir")
    nu = sedmod.to_frequency(p)
    expected_hz = (const.c / (2.2 * u.micron)).to(u.Hz).value
    assert nu.unit == u.Hz
    assert nu.value == pytest.approx(expected_hz, rel=1e-12)


def test_to_frequency_from_energy() -> None:
    p = sedmod.SEDPoint(label="2 keV", energy=2.0 * u.keV, band="xray")
    nu = sedmod.to_frequency(p)
    expected_hz = ((2.0 * u.keV) / const.h).to(u.Hz).value
    assert nu.value == pytest.approx(expected_hz, rel=1e-12)


def test_to_frequency_from_frequency_passthrough() -> None:
    p = sedmod.SEDPoint(label="1.4 GHz", frequency=1.4 * u.GHz, band="radio")
    nu = sedmod.to_frequency(p)
    assert nu.to(u.GHz).value == pytest.approx(1.4, rel=1e-12)


def test_to_frequency_raises_without_axis_value() -> None:
    p = sedmod.SEDPoint(label="bad", band="unknown")
    with pytest.raises(ValueError, match="wavelength, frequency, or energy"):
        sedmod.to_frequency(p)


def test_to_nu_fnu_one_jansky_at_1ghz_matches_hand_value() -> None:
    """νF_ν at 1 GHz, 1 Jy.

    F_ν = 1 Jy = 1e-23 erg/s/cm²/Hz.
    νF_ν = 1e9 Hz * 1e-23 erg/s/cm²/Hz = 1e-14 erg/s/cm².
    """
    p = sedmod.SEDPoint(label="1 Jy at 1 GHz",
                        frequency=1.0 * u.GHz,
                        flux_density=1.0 * u.Jy,
                        band="radio")
    nufnu = sedmod.to_nu_fnu(p)
    assert nufnu.unit == u.erg / u.s / u.cm**2
    assert nufnu.value == pytest.approx(1e-14, rel=1e-12)


def test_to_nu_fnu_passthrough_when_already_set() -> None:
    p = sedmod.SEDPoint(label="x-ray", energy=2.0 * u.keV,
                        nu_f_nu=3.5e-12 * (u.erg / u.s / u.cm**2),
                        band="xray")
    out = sedmod.to_nu_fnu(p)
    assert out.value == pytest.approx(3.5e-12, rel=1e-12)


def test_to_nu_fnu_raises_when_no_flux_supplied() -> None:
    p = sedmod.SEDPoint(label="nope", wavelength=1.0 * u.micron, band="ir")
    with pytest.raises(ValueError, match="flux_density or nu_f_nu"):
        sedmod.to_nu_fnu(p)


def test_sed_add_and_iter() -> None:
    s = sedmod.SED(target_name="NGC 1068")
    s.add(sedmod.SEDPoint(label="a", wavelength=2.2 * u.micron,
                          flux_density=1.0 * u.Jy, band="ir"))
    s.add(sedmod.SEDPoint(label="b", energy=2.0 * u.keV,
                          nu_f_nu=1e-12 * (u.erg / u.s / u.cm**2),
                          band="xray"))
    assert len(s.points) == 2
    assert s.target_name == "NGC 1068"
    assert s.points[0].band == "ir"
    assert s.points[1].band == "xray"


def test_render_sed_runs_and_returns_figure() -> None:
    """Smoke test: rendering an SED with one point returns a Matplotlib Figure."""
    import matplotlib
    matplotlib.use("Agg")
    s = sedmod.SED(target_name="Test")
    s.add(sedmod.SEDPoint(label="K", wavelength=2.2 * u.micron,
                          flux_density=1.0 * u.Jy, band="ir"))
    fig = sedmod.render_sed(s)
    assert fig is not None
    # Two scatter calls or one — at minimum, one axes with data.
    assert len(fig.axes) >= 1


def test_render_sed_empty_does_not_crash() -> None:
    import matplotlib
    matplotlib.use("Agg")
    s = sedmod.SED(target_name="Empty")
    fig = sedmod.render_sed(s)
    assert fig is not None


def test_band_colors_cover_all_used_bands() -> None:
    for band in ("radio", "submm", "ir", "opt", "uv", "xray", "gamma", "unknown"):
        assert band in sedmod.BAND_COLORS


def test_to_nu_fnu_at_2_keV_matches_energy_density() -> None:
    """νF_ν at 2 keV given F_ν directly (in CGS units)."""
    nu = ((2.0 * u.keV) / const.h).to(u.Hz)
    fnu = 1e-30 * u.erg / u.s / u.cm**2 / u.Hz   # arbitrary
    p = sedmod.SEDPoint(label="ck", energy=2.0 * u.keV,
                        flux_density=fnu, band="xray")
    nufnu = sedmod.to_nu_fnu(p)
    expected = (nu.value * fnu.value)
    assert nufnu.value == pytest.approx(expected, rel=1e-12)


def test_to_frequency_keV_consistent_with_x_ray_band() -> None:
    p = sedmod.SEDPoint(label="hard", energy=10.0 * u.keV, band="xray")
    nu = sedmod.to_frequency(p)
    # 10 keV ~ 2.4e18 Hz
    assert nu.value == pytest.approx(2.4e18, rel=2e-2)


def test_unit_round_trip_jy_to_cgs_to_jy() -> None:
    fnu = 1.0 * u.Jy
    cgs = fnu.to(u.erg / u.s / u.cm**2 / u.Hz)
    back = cgs.to(u.Jy)
    assert back.value == pytest.approx(1.0, rel=1e-12)
    assert cgs.value == pytest.approx(1e-23, rel=1e-12)


@pytest.mark.parametrize("wl_um,fnu_jy", [(2.2, 1.0), (12.0, 18.0), (24.0, 20.0)])
def test_to_nu_fnu_parametrized_ir_points(wl_um: float, fnu_jy: float) -> None:
    p = sedmod.SEDPoint(label=f"{wl_um}um", wavelength=wl_um * u.micron,
                        flux_density=fnu_jy * u.Jy, band="ir")
    nufnu = sedmod.to_nu_fnu(p).value
    nu = sedmod.to_frequency(p).value
    fnu_cgs = fnu_jy * 1e-23
    assert nufnu == pytest.approx(nu * fnu_cgs, rel=1e-10)


def test_array_wavelength_and_energy_compatible() -> None:
    """One point with wavelength and another with energy should both convert
    cleanly to a common frequency-space representation."""
    a = sedmod.SEDPoint(label="IR", wavelength=22.0 * u.micron,
                        flux_density=10.0 * u.Jy, band="ir")
    b = sedmod.SEDPoint(label="X-ray", energy=5.0 * u.keV,
                        nu_f_nu=1e-12 * (u.erg / u.s / u.cm**2), band="xray")
    nu_a = sedmod.to_frequency(a).to(u.Hz).value
    nu_b = sedmod.to_frequency(b).to(u.Hz).value
    assert nu_b > nu_a   # X-ray higher frequency than mid-IR


def test_no_negative_nufnu_for_positive_inputs() -> None:
    p = sedmod.SEDPoint(label="pos", wavelength=1.0 * u.micron,
                        flux_density=0.5 * u.Jy, band="ir")
    assert sedmod.to_nu_fnu(p).value > 0


def test_overplot_template_helper_does_not_crash() -> None:
    """The Elvis+1994 schematic overplot helper is private but exercised here."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    sedmod._overplot_elvis_template(ax)
    plt.close(fig)


def test_render_sed_with_overplot_runs() -> None:
    import matplotlib
    matplotlib.use("Agg")
    s = sedmod.SED(target_name="t")
    s.add(sedmod.SEDPoint(label="K", wavelength=2.2 * u.micron,
                          flux_density=np.float64(1.0) * u.Jy, band="ir"))
    fig = sedmod.render_sed(s, overplot_quasar_template=True)
    assert fig is not None
