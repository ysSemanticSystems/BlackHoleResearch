"""
blackhole.catalog — typed catalog of studied black-hole sources.

This module is the single source of truth for source-level metadata
(mass, distance, redshift, classification). The rest of the codebase
reads from `CATALOG` rather than rederiving or hardcoding these values.

Design notes
------------
- `Source` is frozen; the `CATALOG` tuple is immutable at import time.
  Mutation at runtime is not supported — fork the module if you need a
  scenario-specific override.
- Every numeric field has a primary-source reference string. The test
  suite (`tests/test_catalog.py`) pins each value to its citation.
- The `short_id` is the filename prefix used by `scripts/download_data.py`
  and the in-page UI selector. Resolving a FITS filename to a catalog
  entry should go through `by_filename()`, not ad hoc string matching.

References
----------
- SIMBAD source resolver: https://simbad.u-strasbg.fr/simbad/
- NED extragalactic database: https://ned.ipac.caltech.edu/
- Eddington luminosity: Frank, King & Raine 2002, "Accretion Power in
  Astrophysics", 3rd ed., §1.2 (Eq. 1.5).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import astropy.units as u
from astropy.coordinates import SkyCoord

from .physics.accretion import eddington_luminosity

SourceType = Literal[
    "seyfert1",
    "seyfert2",
    "llagn",
    "xrb_hmxb",
    "xrb_lmxb",
    "lrd",
    "quasar",
    "blazar",
    "tde",
]


@dataclass(frozen=True)
class Source:
    """A single studied black-hole source.

    Every numeric field carries an explicit `*_ref` companion field
    pointing to its primary source. The convention is "First-author+YEAR
    Journal volume page" so a reader can copy-paste into ADS.

    Attributes
    ----------
    name
        Canonical display name, e.g. "NGC 1068".
    short_id
        Lower-case identifier with no whitespace; used as the filename
        prefix by `scripts/download_data.py` and matched against FITS
        filenames in `by_filename()`.
    aliases
        Tuple of alternative spellings/IDs (SIMBAD-resolvable, NED IDs,
        common abbreviations, filename variants).
    coord
        ICRS SkyCoord with explicit `unit=(u.deg, u.deg)`. Used for any
        future SIMBAD/NED query layer; for the UI it's display-only.
    redshift
        Heliocentric redshift z. None for Galactic sources.
    distance_mpc
        Luminosity distance in Mpc. For Galactic XRBs we still report a
        distance in Mpc (converted from kpc) so a single number suffices
        downstream; tag the source `type='xrb_*'` to indicate scale.
    distance_ref
        Primary reference string for the distance value.
    m_bh_msun
        Black-hole mass in M_sun. None if no robust dynamical / reverberation
        / parallax measurement exists.
    m_bh_err_msun
        Symmetric 1-sigma uncertainty on M_BH (Msun). Asymmetric errors
        are documented in `notes` and we store the larger side here.
    m_bh_ref
        Primary reference for the M_BH value.
    type
        Classification literal; one of `SourceType`.
    notes
        Free-form caveats: alternate values in the literature, asymmetric
        errors, anything a downstream consumer needs to know.

    References
    ----------
    Each instance documents its own measurements; see the `CATALOG`
    definitions in this module for the field-by-field provenance.
    """

    name: str
    short_id: str
    aliases: tuple[str, ...]
    coord: SkyCoord
    redshift: float | None
    distance_mpc: float | None
    distance_ref: str
    m_bh_msun: float | None
    m_bh_err_msun: float | None
    m_bh_ref: str
    type: SourceType
    notes: str = field(default="")


# ---------------------------------------------------------------------------
# Catalog
# ---------------------------------------------------------------------------
#
# Values below are pinned to specific primary-source papers. When you
# update an entry, also update the corresponding test in
# tests/test_catalog.py so silent drift is impossible.
# ---------------------------------------------------------------------------

NGC1068 = Source(
    name="NGC 1068",
    short_id="ngc1068",
    aliases=("NGC1068", "M77", "Messier 77", "3C 71", "NGC_1068", "n1068"),
    coord=SkyCoord(ra=40.66962, dec=-0.01328, unit=(u.deg, u.deg), frame="icrs"),
    redshift=0.003793,
    distance_mpc=14.4,
    distance_ref="Tully+2013 AJ 146 86 (CosmicFlows-2; redshift-independent)",
    # H2O-maser disk dynamics, corrected for disk self-gravity.
    # See notes for the pre-2003 maser-only value, which is still widely cited.
    m_bh_msun=8.0e6,
    m_bh_err_msun=0.3e6,
    m_bh_ref="Lodato & Bertin 2003 A&A 408 1015",
    type="seyfert2",
    notes=(
        "Compton-thick (N_H > 1.5e24 cm^-2; Marinucci+2016 NuSTAR). "
        "Earlier maser-only analysis (Greenhill+1996 ApJ 472 L21) gave "
        "1.7e7 Msun; the Lodato+2003 self-gravity correction is preferred."
    ),
)

M87 = Source(
    name="M87",
    short_id="m87",
    aliases=("Messier 87", "NGC 4486", "NGC4486", "Virgo A", "3C 274", "m_87"),
    coord=SkyCoord(ra=187.70593, dec=12.39112, unit=(u.deg, u.deg), frame="icrs"),
    redshift=0.00428,
    distance_mpc=16.8,
    distance_ref=(
        "Bird+2010 A&A 524 A71 (surface-brightness fluctuations); "
        "consistent with EHT Coll. 2019 Paper VI"
    ),
    m_bh_msun=6.5e9,
    m_bh_err_msun=0.7e9,
    m_bh_ref="EHT Collaboration 2019 ApJL 875 L6",
    type="llagn",
    notes=(
        "First direct event-horizon image (EHT 2019). "
        "Gas-dynamical (Walsh+2013) and stellar-dynamical (Gebhardt+2011) "
        "M_BH estimates differ by ~2x; EHT shadow size is the cleanest constraint."
    ),
)

CYG_X1 = Source(
    name="Cyg X-1",
    short_id="cygx1",
    aliases=("Cygnus X-1", "CygX-1", "Cyg X1", "HDE 226868", "V1357 Cyg", "cyg_x1"),
    coord=SkyCoord(ra=299.59032, dec=35.20162, unit=(u.deg, u.deg), frame="icrs"),
    redshift=None,
    # 2.22 +/- 0.18 kpc from VLBI parallax. Stored as Mpc for type uniformity.
    distance_mpc=2.22e-3,
    distance_ref="Miller-Jones+2021 Science 371 1046 (VLBI parallax)",
    m_bh_msun=21.2,
    m_bh_err_msun=2.2,
    m_bh_ref="Miller-Jones+2021 Science 371 1046",
    type="xrb_hmxb",
    notes=(
        "High-mass X-ray binary; companion HDE 226868 is an O9.7 Iab supergiant. "
        "Earlier mass estimates were ~15 Msun; the upward revision in "
        "Miller-Jones+2021 follows from the new parallax-driven distance."
    ),
)


CATALOG: tuple[Source, ...] = (NGC1068, M87, CYG_X1)


# ---------------------------------------------------------------------------
# Lookups
# ---------------------------------------------------------------------------

def by_short_id(short_id: str) -> Source | None:
    """Return the catalog entry whose `short_id` matches (case-insensitive).

    Examples
    --------
    >>> by_short_id("ngc1068").name
    'NGC 1068'
    >>> by_short_id("missing") is None
    True
    """
    sid = short_id.lower()
    for s in CATALOG:
        if s.short_id == sid:
            return s
    return None


def by_name(name: str) -> Source | None:
    """Return the catalog entry matching `name` exactly or via aliases.

    Case- and whitespace-tolerant; tries exact name match, then aliases.
    Does *not* perform fuzzy matching — callers needing fuzzy lookup
    should use SIMBAD via the (future) `resolve_simbad` helper.
    """
    target = _normalize(name)
    for s in CATALOG:
        if _normalize(s.name) == target:
            return s
        for alias in s.aliases:
            if _normalize(alias) == target:
                return s
    return None


def by_filename(filename: str) -> Source | None:
    """Resolve a FITS filename (or basename) to a catalog entry.

    Matches against `short_id` first (the convention used by
    `scripts/download_data.py`), then against each entry's aliases.
    Punctuation is stripped before comparison so `cygx1`, `cyg-x1`,
    `cyg_x1` all resolve.

    Returns
    -------
    Source | None
        The matching entry, or None if no catalog entry matches.
    """
    needle = _normalize(filename)
    for s in CATALOG:
        if s.short_id in needle:
            return s
        for alias in s.aliases:
            if _normalize(alias) and _normalize(alias) in needle:
                return s
    return None


def _normalize(s: str) -> str:
    """Lowercase, strip whitespace/punctuation used inconsistently across IDs."""
    return (
        s.lower()
        .replace(" ", "")
        .replace("_", "")
        .replace("-", "")
        .replace(".", "")
    )


# ---------------------------------------------------------------------------
# Derived physical quantities
# ---------------------------------------------------------------------------

def eddington_luminosity_of(source: Source | str) -> u.Quantity:
    """Compute the Eddington luminosity from the catalog M_BH.

    Parameters
    ----------
    source
        Either a `Source` instance or its `short_id` / canonical name.

    Returns
    -------
    Quantity
        L_Edd in erg/s. Raises if the source has no catalogued M_BH.

    Notes
    -----
    Delegates to `blackhole.physics.accretion.eddington_luminosity`, which
    uses the standard pure-hydrogen Thomson-opacity expression
    L_Edd = 4 pi G M m_p c / sigma_T  ~  1.26e38 (M / Msun) erg/s
    (Frank, King & Raine 2002, Eq. 1.5).
    """
    s = _resolve(source)
    if s.m_bh_msun is None:
        raise ValueError(
            f"Source {s.name!r} has no catalogued M_BH; "
            "Eddington luminosity is undefined."
        )
    return eddington_luminosity(s.m_bh_msun)


def distance_to(source: Source | str) -> u.Quantity:
    """Luminosity distance as an astropy Quantity (Mpc).

    Raises if the source has no distance entry. For Galactic XRBs we
    return Mpc (1 kpc = 1e-3 Mpc); callers should `.to(u.kpc)` for
    readability when the value is < 0.1 Mpc.
    """
    s = _resolve(source)
    if s.distance_mpc is None:
        raise ValueError(f"Source {s.name!r} has no catalogued distance.")
    return s.distance_mpc * u.Mpc


def redshift_of(source: Source | str) -> float | None:
    """Heliocentric redshift z, or None for Galactic sources."""
    return _resolve(source).redshift


def _resolve(source: Source | str) -> Source:
    """Internal: accept either a Source or a name/short_id."""
    if isinstance(source, Source):
        return source
    candidate = by_short_id(source) or by_name(source) or by_filename(source)
    if candidate is None:
        raise KeyError(f"No catalog entry resolves {source!r}.")
    return candidate
