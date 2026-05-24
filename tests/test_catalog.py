"""
tests/test_catalog.py — pin every published value in blackhole.catalog
to its primary source.

A test failure here means a catalog entry drifted from the paper it
cites. If the drift is intentional (e.g. you updated to a newer paper),
update both the catalog entry's `*_ref` field and the test in the same
commit.

References
----------
- NGC 1068:  Lodato & Bertin 2003 A&A 408 1015; Tully+2013 AJ 146 86.
- M87:       EHT Collaboration 2019 ApJL 875 L6; Bird+2010 A&A 524 A71.
- Cyg X-1:   Miller-Jones+2021 Science 371 1046.
- L_Edd:     Frank, King & Raine 2002, Eq. 1.5.
"""

from __future__ import annotations

import dataclasses

import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord

from blackhole import catalog as cat

# ---------------------------------------------------------------------------
# Pinned values — change only if you also update the catalog reference.
# ---------------------------------------------------------------------------


def test_catalog_has_phase1_targets() -> None:
    """Every Phase-1 target must be present and resolvable by short_id."""
    assert cat.by_short_id("ngc1068") is cat.NGC1068
    assert cat.by_short_id("m87") is cat.M87
    assert cat.by_short_id("cygx1") is cat.CYG_X1
    assert len(cat.CATALOG) >= 3


# --- NGC 1068 --------------------------------------------------------------


def test_ngc1068_mass_matches_lodato_bertin_2003() -> None:
    s = cat.NGC1068
    assert s.m_bh_msun == pytest.approx(8.0e6, rel=1e-9)
    assert s.m_bh_err_msun == pytest.approx(0.3e6, rel=1e-9)
    assert "Lodato" in s.m_bh_ref and "2003" in s.m_bh_ref


def test_ngc1068_distance_and_redshift() -> None:
    s = cat.NGC1068
    assert s.distance_mpc == pytest.approx(14.4, rel=1e-9)
    assert s.redshift == pytest.approx(0.003793, rel=1e-6)
    assert "Tully" in s.distance_ref


def test_ngc1068_classification_and_coord() -> None:
    s = cat.NGC1068
    assert s.type == "seyfert2"
    assert isinstance(s.coord, SkyCoord)
    assert s.coord.ra.deg == pytest.approx(40.66962, abs=1e-3)
    assert s.coord.dec.deg == pytest.approx(-0.01328, abs=1e-3)


# --- M87 -------------------------------------------------------------------


def test_m87_mass_matches_eht_2019() -> None:
    s = cat.M87
    assert s.m_bh_msun == pytest.approx(6.5e9, rel=1e-9)
    assert s.m_bh_err_msun == pytest.approx(0.7e9, rel=1e-9)
    assert "EHT" in s.m_bh_ref and "2019" in s.m_bh_ref


def test_m87_distance() -> None:
    s = cat.M87
    assert s.distance_mpc == pytest.approx(16.8, rel=1e-9)
    assert "Bird" in s.distance_ref


def test_m87_classification_and_coord() -> None:
    s = cat.M87
    assert s.type == "llagn"
    assert s.coord.ra.deg == pytest.approx(187.70593, abs=1e-3)
    assert s.coord.dec.deg == pytest.approx(12.39112, abs=1e-3)


# --- Cyg X-1 ---------------------------------------------------------------


def test_cygx1_mass_matches_miller_jones_2021() -> None:
    s = cat.CYG_X1
    assert s.m_bh_msun == pytest.approx(21.2, rel=1e-9)
    assert s.m_bh_err_msun == pytest.approx(2.2, rel=1e-9)
    assert "Miller-Jones" in s.m_bh_ref and "2021" in s.m_bh_ref


def test_cygx1_distance_in_mpc_corresponds_to_2_22_kpc() -> None:
    s = cat.CYG_X1
    assert s.distance_mpc == pytest.approx(2.22e-3, rel=1e-9)
    # Converted to kpc, it should be 2.22.
    d = cat.distance_to(s).to(u.kpc).value
    assert d == pytest.approx(2.22, rel=1e-9)


def test_cygx1_is_galactic_no_redshift() -> None:
    assert cat.CYG_X1.redshift is None
    assert cat.CYG_X1.type == "xrb_hmxb"


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------


def test_eddington_luminosity_ngc1068_textbook_value() -> None:
    """L_Edd(8e6 Msun) = 1.26e38 * 8e6 = 1.008e45 erg/s.

    Hand-check from Frank, King & Raine 2002 Eq. 1.5.
    """
    L = cat.eddington_luminosity_of("ngc1068")
    assert L.unit == u.erg / u.s
    assert L.value == pytest.approx(1.008e45, rel=5e-3)


def test_eddington_luminosity_m87_textbook_value() -> None:
    """L_Edd(6.5e9 Msun) ~ 8.19e47 erg/s."""
    L = cat.eddington_luminosity_of(cat.M87)
    assert L.value == pytest.approx(8.19e47, rel=5e-3)


def test_eddington_luminosity_cygx1_textbook_value() -> None:
    """L_Edd(21.2 Msun) ~ 2.67e39 erg/s.

    Reference: Frank, King & Raine 2002 Eq. 1.5, 1.26e38 * 21.2.
    """
    L = cat.eddington_luminosity_of("cygx1")
    assert L.value == pytest.approx(2.67e39, rel=5e-3)


def test_eddington_luminosity_raises_when_no_mass() -> None:
    """If we ever add a source with m_bh_msun=None, L_Edd must raise."""
    placeholder = cat.Source(
        name="placeholder",
        short_id="pl",
        aliases=(),
        coord=SkyCoord(ra=0.0, dec=0.0, unit=(u.deg, u.deg), frame="icrs"),
        redshift=None,
        distance_mpc=None,
        distance_ref="",
        m_bh_msun=None,
        m_bh_err_msun=None,
        m_bh_ref="",
        type="lrd",
        notes="",
    )
    with pytest.raises(ValueError, match="no catalogued M_BH"):
        cat.eddington_luminosity_of(placeholder)


# ---------------------------------------------------------------------------
# Lookups
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected",
    [
        ("NGC 1068",  "NGC 1068"),
        ("ngc 1068",  "NGC 1068"),
        ("NGC1068",   "NGC 1068"),
        ("M77",       "NGC 1068"),
        ("Messier 77","NGC 1068"),
        ("M87",       "M87"),
        ("NGC 4486",  "M87"),
        ("Virgo A",   "M87"),
        ("Cyg X-1",   "Cyg X-1"),
        ("Cygnus X-1","Cyg X-1"),
        ("CygX-1",    "Cyg X-1"),
    ],
)
def test_by_name_resolves_aliases(name: str, expected: str) -> None:
    s = cat.by_name(name)
    assert s is not None
    assert s.name == expected


def test_by_name_returns_none_for_unknown() -> None:
    assert cat.by_name("a galaxy we never observed") is None


@pytest.mark.parametrize(
    "filename,expected",
    [
        ("ngc1068_chandra.fits",     "NGC 1068"),
        ("NGC1068_HST.fits",         "NGC 1068"),
        ("ngc_1068_xmm.fits",        "NGC 1068"),
        ("m87_radio_5GHz.fits",      "M87"),
        ("ngc4486_optical.fits",     "M87"),
        ("M-87_chandra.fits",        "M87"),
        ("cygx1_2mass_k.fits",       "Cyg X-1"),
        ("Cyg_X-1_RXTE.fits",        "Cyg X-1"),
        ("cyg-x1_evt.fits",          "Cyg X-1"),
    ],
)
def test_by_filename_resolves_real_filenames(filename: str, expected: str) -> None:
    s = cat.by_filename(filename)
    assert s is not None
    assert s.name == expected


def test_by_filename_returns_none_for_unknown() -> None:
    assert cat.by_filename("totally_unrelated_file.fits") is None


def test_distance_to_returns_quantity() -> None:
    d = cat.distance_to("m87")
    assert d.unit == u.Mpc
    assert d.value == pytest.approx(16.8, rel=1e-9)


def test_distance_to_raises_when_missing() -> None:
    placeholder = cat.Source(
        name="pl", short_id="pl", aliases=(),
        coord=SkyCoord(0, 0, unit=(u.deg, u.deg), frame="icrs"),
        redshift=None, distance_mpc=None, distance_ref="",
        m_bh_msun=None, m_bh_err_msun=None, m_bh_ref="",
        type="lrd", notes="",
    )
    with pytest.raises(ValueError, match="no catalogued distance"):
        cat.distance_to(placeholder)


def test_redshift_of() -> None:
    assert cat.redshift_of("ngc1068") == pytest.approx(0.003793, rel=1e-6)
    assert cat.redshift_of("cygx1") is None


def test_resolve_raises_keyerror_for_unknown_string() -> None:
    with pytest.raises(KeyError):
        cat._resolve("nothing_here")


def test_source_is_frozen() -> None:
    """The dataclass is immutable so callers can't accidentally edit values."""
    s = cat.NGC1068
    with pytest.raises(dataclasses.FrozenInstanceError):
        s.m_bh_msun = 0.0  # type: ignore[misc]


def test_catalog_is_a_tuple_not_a_list() -> None:
    """CATALOG must be a tuple so it cannot be appended to at runtime."""
    assert isinstance(cat.CATALOG, tuple)
