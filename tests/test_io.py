"""
tests/test_io.py — seed tests for blackhole.io.

Covers inspect, header_dict, load_image (with and without WCS), load_events,
and bin_to_image. The fixtures live in tests/conftest.py.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from blackhole import io as bhio


def test_inspect_image(tiny_image_fits: Path) -> None:
    hdus = bhio.inspect(tiny_image_fits)
    assert len(hdus) == 1
    primary = hdus[0]
    assert primary.hdu_type == "PrimaryHDU"
    assert primary.telescope == "SYNTHETIC"
    assert primary.instrument == "TEST"
    assert primary.exposure_s == pytest.approx(1000.0)
    assert primary.shape == (64, 64)


def test_inspect_events(tiny_events_fits: Path) -> None:
    hdus = bhio.inspect(tiny_events_fits)
    assert len(hdus) == 3
    primary, events, gti = hdus
    assert primary.hdu_type == "PrimaryHDU"
    assert events.name == "EVENTS"
    assert events.hdu_type == "BinTableHDU"
    assert events.shape == (1000,)
    assert events.telescope == "CHANDRA"
    assert events.instrument == "ACIS"
    assert gti.name == "GTI"
    assert events.columns is not None and "TIME" in events.columns


def test_header_dict_strips_history(tiny_image_fits: Path) -> None:
    d = bhio.header_dict(tiny_image_fits, hdu_index=0)
    assert "COMMENT" not in d
    assert "HISTORY" not in d
    assert d.get("TELESCOP") == "SYNTHETIC"


def test_load_image_returns_array_and_wcs(tiny_image_fits: Path) -> None:
    img = bhio.load_image(tiny_image_fits)
    assert img.array.shape == (64, 64)
    assert img.array.dtype.kind == "f"
    assert img.wcs is not None
    assert isinstance(img.wcs, WCS)
    assert img.wcs.has_celestial


def test_load_image_copies_array(tmp_path: Path) -> None:
    """Verify the loaded array survives file closure (PITFALLS #9)."""
    data = np.arange(16, dtype=np.float32).reshape(4, 4)
    path = tmp_path / "mini.fits"
    fits.PrimaryHDU(data=data).writeto(path)

    img = bhio.load_image(path)
    # Force-mutating the file under us must not corrupt the loaded copy.
    (tmp_path / "mini.fits").unlink()
    assert img.array.shape == (4, 4)
    assert float(img.array[3, 3]) == pytest.approx(15.0)


def test_load_image_no_wcs_returns_none(tmp_path: Path) -> None:
    """Headers without CTYPE keywords must yield img.wcs == None."""
    data = np.zeros((8, 8), dtype=np.float32)
    path = tmp_path / "nowcs.fits"
    fits.PrimaryHDU(data=data).writeto(path)
    img = bhio.load_image(path)
    assert img.wcs is None


def test_load_events_columns_and_mission(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    assert ev.times.size == 1000
    assert ev.x.size == 1000
    assert ev.y.size == 1000
    assert ev.energies is not None
    assert ev.energy_unit == "eV"
    assert ev.mission == "CHANDRA"
    assert ev.times.min() >= 0.0
    assert ev.times.max() <= 1000.0


def test_bin_to_image_returns_correct_shape(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    arr, extent = bhio.bin_to_image(ev, bins=128)
    assert arr.shape == (128, 128)
    assert int(arr.sum()) == ev.times.size
    xmin, xmax, ymin, ymax = extent
    assert xmin < xmax and ymin < ymax


def test_bin_to_image_energy_filter_drops_counts(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    arr_full, _ = bhio.bin_to_image(ev, bins=64)
    arr_hard, _ = bhio.bin_to_image(ev, bins=64, energy_range_ev=(2000, 8000))
    assert arr_hard.sum() < arr_full.sum()
    assert arr_hard.sum() > 0


def test_bin_to_image_empty_returns_zeros(tiny_events_fits: Path) -> None:
    ev = bhio.load_events(tiny_events_fits)
    arr, _ = bhio.bin_to_image(ev, bins=32, energy_range_ev=(1e6, 1e7))
    assert arr.shape == (32, 32)
    assert int(arr.sum()) == 0


def test_iter_fits_finds_extensions(tmp_path: Path) -> None:
    (tmp_path / "a.fits").write_bytes(b"")
    (tmp_path / "b.fit").write_bytes(b"")
    (tmp_path / "c.txt").write_bytes(b"")
    found = sorted(p.name for p in bhio.iter_fits(tmp_path))
    assert found == ["a.fits", "b.fit"]
