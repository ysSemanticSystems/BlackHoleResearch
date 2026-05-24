"""
docs/literature_seds.py — published SED reference points.

This module exists so that the literature comparison values are
*data*, not *logic*, and so that `app.py` does not embed numerical
photometry in code. Aperture photometry from the local FITS cutouts
(via `blackhole.photometry.aperture_photometry_on`) is the primary
data path; the values here are an opt-in overlay for visual comparison.

Numbers are illustrative, drawn from NED compilations, AllWISE,
2MASS XSC, and the named papers. They are NOT measurement-quality
values for new science — they exist as a calibration check for the
Phase-2 photometry pipeline.

Schema
------
SEDS:  {target_name -> tuple[(label, wavelength_um, flux_density_jy, band, source), ...]}
XRAY:  {target_name -> tuple[(energy_keV, nu_f_nu_cgs, label, source), ...]}

Bands match `blackhole.sed.BAND_COLORS`.
"""

from __future__ import annotations

LiteraturePoint = tuple[str, float, float, str, str]
XRayLiteraturePoint = tuple[float, float, str, str]


SEDS: dict[str, tuple[LiteraturePoint, ...]] = {
    "NGC 1068": (
        # (label, wavelength_um, flux_density_jy, band, source)
        ("VLA 1.4 GHz",        214000.0, 2.0,    "radio", "NED"),
        ("VLA 5 GHz",           60000.0, 1.4,    "radio", "NED"),
        ("Spitzer 24 µm",          24.0, 20.0,   "ir",    "Bendo+2012 MNRAS 419 1833"),
        ("WISE W4 22 µm",          22.0, 18.0,   "ir",    "AllWISE"),
        ("WISE W3 12 µm",          12.0, 18.0,   "ir",    "AllWISE"),
        ("Spitzer 8 µm",            8.0,  7.0,   "ir",    "Bendo+2012 MNRAS 419 1833"),
        ("WISE W2 4.6 µm",          4.6,  2.6,   "ir",    "AllWISE"),
        ("WISE W1 3.4 µm",          3.4,  1.6,   "ir",    "AllWISE"),
        ("2MASS K 2.2 µm",          2.2,  1.0,   "ir",    "2MASS XSC"),
        ("2MASS J 1.25 µm",        1.25,  0.55,  "ir",    "2MASS XSC"),
        ("DSS R 0.66 µm",          0.66,  0.07,  "opt",   "DSS"),
        ("Swift UVOT UVW1",        0.26,  0.005, "uv",    "Swift UVOT"),
    ),
    "M87": (
        ("VLA 1.4 GHz core",  214000.0, 4.0,    "radio", "NED"),
        ("VLA 5 GHz core",     60000.0, 2.9,    "radio", "NED"),
        ("Spitzer 24 µm",         24.0, 0.4,    "ir",    "Shi+2007 ApJ 655 781"),
        ("WISE W3 12 µm",         12.0, 0.5,    "ir",    "AllWISE"),
        ("WISE W1 3.4 µm",         3.4, 0.5,    "ir",    "AllWISE"),
        ("2MASS K 2.2 µm",         2.2, 1.5,    "ir",    "2MASS XSC"),
        ("DSS R 0.66 µm",         0.66, 0.6,    "opt",   "DSS"),
    ),
    "Cyg X-1": (
        ("2MASS K 2.2 µm",   2.2,  0.8,  "ir",  "2MASS"),
        ("2MASS J 1.25 µm",  1.25, 0.55, "ir",  "2MASS"),
        ("DSS R 0.66 µm",    0.66, 0.4,  "opt", "DSS"),
    ),
}


XRAY: dict[str, tuple[XRayLiteraturePoint, ...]] = {
    "NGC 1068": (
        (2.0,  5e-13, "Chandra 2 keV",   "Bauer+2015 ApJ 812 116"),
        (10.0, 1e-12, "NuSTAR 10 keV",   "Marinucci+2016 MNRAS 456 L94"),
        (30.0, 5e-12, "NuSTAR 30 keV",   "Marinucci+2016 MNRAS 456 L94"),
    ),
    "M87": (
        (1.0,  5e-12, "Chandra 1 keV",   "Wilson & Yang 2002 ApJ 568 133"),
        (5.0,  2e-12, "Chandra 5 keV",   "Wilson & Yang 2002 ApJ 568 133"),
    ),
    "Cyg X-1": (
        (2.0,  3e-9,  "RXTE 2 keV",       "Wilms+2006 A&A 447 245"),
        (10.0, 2e-9,  "RXTE 10 keV",      "Wilms+2006 A&A 447 245"),
        (50.0, 5e-10, "INTEGRAL 50 keV",  "Wilms+2006 A&A 447 245"),
    ),
}


def literature_for(target_name: str) -> tuple[
    tuple[LiteraturePoint, ...], tuple[XRayLiteraturePoint, ...]
]:
    """Return (SEDS, XRAY) for the named target. Returns empty tuples if
    the target is absent."""
    return SEDS.get(target_name, ()), XRAY.get(target_name, ())
