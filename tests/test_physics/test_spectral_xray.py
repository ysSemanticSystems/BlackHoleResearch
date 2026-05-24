"""
Physics tests for blackhole.physics.spectral_xray — hardness ratios and
photon-index classification.
"""

from __future__ import annotations

import numpy as np
import pytest

from blackhole.physics import spectral_xray as sx


@pytest.mark.parametrize(
    "soft,hard,expected",
    [
        (100, 100, 0.0),
        (0, 100, 1.0),
        (100, 0, -1.0),
        (50, 150, 0.5),
        (150, 50, -0.5),
    ],
)
def test_hardness_ratio_known_values(soft: float, hard: float, expected: float) -> None:
    assert sx.hardness_ratio(soft, hard) == pytest.approx(expected, abs=1e-12)


def test_hardness_ratio_zero_total_returns_nan() -> None:
    result = sx.hardness_ratio(0, 0)
    assert np.isnan(float(result))


def test_hardness_ratio_vectorized() -> None:
    s = np.array([100, 50, 0])
    h = np.array([100, 150, 100])
    hr = sx.hardness_ratio(s, h)
    assert hr.shape == (3,)
    assert hr[0] == pytest.approx(0.0)
    assert hr[1] == pytest.approx(0.5)
    assert hr[2] == pytest.approx(1.0)


def test_hardness_ratio_error_positive_for_positive_counts() -> None:
    err = float(sx.hardness_ratio_error(100, 100))
    assert err > 0
    assert err < 1


def test_hardness_ratio_error_zero_when_zero_total() -> None:
    err = sx.hardness_ratio_error(0, 0)
    assert np.isnan(float(err))


@pytest.mark.parametrize(
    "gamma,expected_keyword",
    [
        (0.5, "unphysically hard"),
        (1.2, "reflection"),
        (1.5, "obscured"),
        (1.9, "type-1"),
        (2.3, "soft"),
        (3.0, "TDE"),
    ],
)
def test_classify_photon_index_returns_meaningful_label(
    gamma: float, expected_keyword: str,
) -> None:
    label = sx.classify_photon_index(gamma)
    assert isinstance(label, str)
    assert expected_keyword.lower() in label.lower()


def test_classify_photon_index_boundary_values() -> None:
    # Boundary conditions in the implementation.
    assert "unphysically" in sx.classify_photon_index(0.99).lower()
    assert "reflection" in sx.classify_photon_index(1.39).lower()
    assert "obscured" in sx.classify_photon_index(1.69).lower()
    assert "type-1" in sx.classify_photon_index(2.09).lower()
    assert "soft" in sx.classify_photon_index(2.49).lower()
