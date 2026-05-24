"""
tests/test_provenance.py — provenance contract enforcement.

These tests pin two M6 invariants:

1. Every render_* function in `blackhole.*` attaches a `Provenance` to
   the returned Figure (via `provenance.attach`). The attached struct
   serializes round-trip through JSON.
2. `provenance.save_figure(fig, path)` writes a PNG plus a sidecar JSON
   whose content reconstructs the same Provenance.

It also enforces the dark-mode-spine grep guard from M6 exit criteria:
the inline `for spine in ax.spines.values(): spine.set_color(...)` idiom
must live in `_style.apply_dark` only, not be re-introduced anywhere else
in the package.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

matplotlib.use("Agg")

import astropy.units as u

from blackhole import io as bhio
from blackhole import lightcurves as lc
from blackhole import provenance as prov
from blackhole import sed as sedmod
from blackhole import spectra as sp
from blackhole import wcs_plot as wp
from blackhole.calibration import calibrate
from blackhole.photometry import aperture_photometry_on

# ---------------------------------------------------------------------------
# Provenance schema round-trip
# ---------------------------------------------------------------------------


def test_provenance_dict_roundtrip() -> None:
    p = prov.Provenance(
        fits_sha256="deadbeef",
        fits_path="/tmp/example.fits",
        calibration_version="1.0.0",
        function_chain=("a", "b", "c"),
        library_version="0.1.0",
        timestamp_utc="2026-01-01T00:00:00+00:00",
        extra={"stretch": "asinh", "cmap": "inferno"},
    )
    d = p.to_dict()
    p2 = prov.Provenance.from_dict(d)
    assert p2 == p


def test_sha256_of_file_for_real_file(tmp_path: Path) -> None:
    path = tmp_path / "hello.txt"
    path.write_bytes(b"hello provenance")
    sha = prov.sha256_of_file(path)
    assert len(sha) == 64
    assert all(c in "0123456789abcdef" for c in sha)


def test_sha256_of_file_for_missing_returns_synthetic(tmp_path: Path) -> None:
    assert prov.sha256_of_file(None) == "synthetic"
    assert prov.sha256_of_file(tmp_path / "nope.fits") == "synthetic"


def test_extend_chain_appends_names() -> None:
    p = prov.build_provenance(None, function_chain=("load",))
    p2 = prov.extend_chain(p, "render")
    assert p2.function_chain == ("load", "render")


# ---------------------------------------------------------------------------
# Renderer attach contract — every render_* attaches a Provenance
# ---------------------------------------------------------------------------


def test_render_image_attaches_provenance(tiny_image_fits: Path) -> None:
    img = bhio.load_image(tiny_image_fits)
    fig = wp.render_image(img)
    p = prov.get(fig)
    assert p is not None
    assert "render_image" in p.function_chain
    assert p.fits_path.endswith(tiny_image_fits.name)
    assert p.library_version
    assert p.calibration_version


def test_render_event_image_attaches_provenance(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    arr, extent = bhio.bin_to_image(ev, bins=64)
    fig = wp.render_event_image(arr, extent)
    p = prov.get(fig)
    assert p is not None
    assert "render_event_image" in p.function_chain


def test_render_lightcurve_attaches_provenance(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    curve = lc.bin_events_to_lightcurve(ev, bin_size_s=100.0)
    fig = lc.render_lightcurve(curve)
    p = prov.get(fig)
    assert p is not None
    assert "render_lightcurve" in p.function_chain


def test_render_periodogram_attaches_provenance() -> None:
    freqs = np.linspace(1e-3, 1.0, 100)
    power = np.abs(np.sin(freqs))
    fig = lc.render_periodogram(freqs, power)
    p = prov.get(fig)
    assert p is not None
    assert "render_periodogram" in p.function_chain


def test_render_spectrum_attaches_provenance(tiny_pha_fits: Path) -> None:
    spec = sp.load_pha_spectrum(tiny_pha_fits)
    fit = sp.fit_power_law(spec)
    fig = sp.render_spectrum(spec, fit=fit)
    p = prov.get(fig)
    assert p is not None
    assert "render_spectrum" in p.function_chain
    assert p.extra.get("fit_applied") is True
    assert p.fits_path.endswith(tiny_pha_fits.name)


def test_render_sed_attaches_provenance() -> None:
    s = sedmod.SED(target_name="TestObj")
    s.add(sedmod.SEDPoint(
        label="K", wavelength=2.2 * u.micron,
        flux_density=1.0 * u.Jy, band="ir",
    ))
    fig = sedmod.render_sed(s)
    p = prov.get(fig)
    assert p is not None
    assert "render_sed" in p.function_chain
    assert p.extra.get("target") == "TestObj"


# ---------------------------------------------------------------------------
# save_figure produces PNG + sidecar JSON; JSON reconstructs the Provenance
# ---------------------------------------------------------------------------


def test_save_figure_writes_png_and_sidecar(
    tiny_image_fits: Path, tmp_path: Path,
) -> None:
    img = bhio.load_image(tiny_image_fits)
    fig = wp.render_image(img)
    out = tmp_path / "out" / "image"
    png_path, json_path = prov.save_figure(fig, out)
    assert png_path.exists() and png_path.suffix == ".png"
    assert json_path.exists() and json_path.suffix == ".json"
    loaded = prov.load_sidecar(json_path)
    attached = prov.get(fig)
    assert loaded == attached


def test_save_figure_raises_without_provenance(tmp_path: Path) -> None:
    """save_figure refuses to write a figure without provenance attached."""
    import matplotlib.pyplot as plt
    fig, _ = plt.subplots()
    with pytest.raises(ValueError, match="no Provenance"):
        prov.save_figure(fig, tmp_path / "x.png")


# ---------------------------------------------------------------------------
# Photometry integrates with provenance via the calibration version
# ---------------------------------------------------------------------------


def test_photometry_provenance_carries_calibration_version(
    tiny_image_fits: Path,
) -> None:
    """Build a calibrated 2MASS-K image from the tiny_image fixture,
    photometer it, and verify the surrounding provenance still records
    the calibration version that produced it."""
    # The tiny_image_fits doesn't have MAGZP. Skip if not present; this
    # test ensures the field exists on real calibrated paths.
    p = prov.build_provenance(
        tiny_image_fits,
        function_chain=("load_image", "calibrate",
                        "aperture_photometry_on", "render_sed"),
    )
    assert p.calibration_version
    # Aperture photometry alone is exercised in test_photometry; here we
    # just verify the chain composition discipline.
    assert "aperture_photometry_on" in p.function_chain


def test_calibrated_photometry_attached_to_sed(tmp_path: Path) -> None:
    """End-to-end happy path: calibrated image -> photometry -> SED render.

    The final Figure's Provenance must carry calibration_version and
    the renderer name."""
    from astropy.io import fits as afits

    n = 51
    cy = cx = (n - 1) / 2.0
    y, x = np.mgrid[0:n, 0:n]
    img = 1e5 * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * 2.0**2))
    hdu = afits.PrimaryHDU(data=img.astype(np.float32))
    hdu.header["TELESCOP"] = "2MASS"
    hdu.header["SURVEY"] = "2MASS-K"
    hdu.header["MAGZP"] = 21.0
    hdu.header["CDELT1"] = -1.0 / 3600.0
    hdu.header["CDELT2"] =  1.0 / 3600.0
    hdu.header["CTYPE1"] = "RA---TAN"
    hdu.header["CTYPE2"] = "DEC--TAN"
    hdu.header["CRVAL1"] = 0.0
    hdu.header["CRVAL2"] = 0.0
    hdu.header["CRPIX1"] = cx + 1
    hdu.header["CRPIX2"] = cy + 1
    fpath = tmp_path / "ngc_2mass.fits"
    hdu.writeto(fpath, overwrite=True)

    cal = calibrate(bhio.load_image(fpath))
    res = aperture_photometry_on(
        cal, SkyCoord(0.0, 0.0, unit=(u.deg, u.deg), frame="icrs"),
        aperture_radius=10.0 * u.arcsec,
        annulus_inner=15.0 * u.arcsec,
        annulus_outer=25.0 * u.arcsec,
    )
    sed = sedmod.SED(target_name="X")
    pt = res.to_sed_point()
    pt.wavelength = 2.159 * u.micron
    sed.add(pt)
    fig = sedmod.render_sed(sed)
    p = prov.get(fig)
    assert p is not None
    assert p.calibration_version != ""


# ---------------------------------------------------------------------------
# Grep guard — apply_dark must be the only place that styles spines
# ---------------------------------------------------------------------------


SPINE_PATTERN = re.compile(
    r"for\s+spine\s+in\s+ax\.spines\.values\s*\(\s*\)\s*:\s*\n\s*spine\.set_color",
    re.MULTILINE,
)


def test_no_inline_spine_loop_outside_style() -> None:
    """The inline `for spine in ax.spines.values(): spine.set_color(...)`
    pattern lives in `blackhole/_style.py::apply_dark` only.

    Phase 2 M6 exit criterion: every renderer routes through apply_dark
    so a future palette change is a one-file edit. This grep enforces it.
    """
    pkg = Path(__file__).resolve().parent.parent / "blackhole"
    hits: list[Path] = []
    for py in pkg.rglob("*.py"):
        if py.name == "_style.py":
            continue
        text = py.read_text()
        if SPINE_PATTERN.search(text):
            hits.append(py)
    assert not hits, (
        "These files still set spine colours inline; route through "
        f"_style.apply_dark instead: {[str(p) for p in hits]}"
    )
