"""
tests/test_wcs_plot.py — renderer smoke tests for blackhole.wcs_plot.

These verify the render functions return a usable Matplotlib Figure with
the expected number of axes, across all supported stretches. No baseline
image regression (planned as optional pytest-mpl in later milestone).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from blackhole import io as bhio
from blackhole import wcs_plot as wp


@pytest.mark.parametrize("stretch", ["linear", "sqrt", "log", "asinh", "zscale"])
def test_render_image_each_stretch(tiny_image_fits: Path, stretch: str) -> None:
    img = bhio.load_image(tiny_image_fits)
    fig = wp.render_image(img, stretch=stretch, cmap="inferno")
    assert fig is not None
    assert len(fig.axes) >= 1


def test_render_image_without_wcs(tmp_path: Path) -> None:
    """When WCS is absent, fall back to pixel axes."""
    from astropy.io import fits
    data = np.zeros((16, 16), dtype=np.float32)
    path = tmp_path / "nowcs.fits"
    fits.PrimaryHDU(data=data).writeto(path)
    img = bhio.load_image(path)
    assert img.wcs is None
    fig = wp.render_image(img, stretch="asinh")
    assert fig is not None


def test_render_event_image_smoke(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    arr, extent = bhio.bin_to_image(ev, bins=64)
    fig = wp.render_event_image(arr, extent, stretch="asinh",
                                cmap="inferno", title="test",
                                energy_band_label="0.5-7 keV")
    assert fig is not None


def test_render_event_image_no_band_label(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    arr, extent = bhio.bin_to_image(ev, bins=32)
    fig = wp.render_event_image(arr, extent)
    assert fig is not None
