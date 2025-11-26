
import math
import numpy as np

from cgeom import clib


def _elliptic_E_numeric(phi, k, n=2000):
    """Simple Simpson integrator for E(phi, k) as an independent oracle."""
    sign = 1 if phi >= 0 else -1
    a, b = 0.0, abs(phi)
    h = (b - a) / n
    total = 0.0
    for i in range(n + 1):
        x = a + i * h
        fx = math.sqrt(1 - k * k * math.sin(x) * math.sin(x))
        coeff = 1 if i == 0 or i == n else 4 if i % 2 == 1 else 2
        total += coeff * fx
    return sign * total * h / 3.0


def test_geom2d_norm():
    assert clib.geom2d_norm(3, 4) == 5.0


def test_geom2d_points_distance():
    assert clib.geom2d_points_distance(1, 2, 4, 6) == 5.0


def test_elliptic_E_matches_numeric_integrator():
    phi = 1.234
    k = 0.6
    expected = _elliptic_E_numeric(phi, k, n=2000)
    assert math.isclose(clib.geom2d_elliptic_E(phi, k), expected, rel_tol=1e-6, abs_tol=1e-6)


def test_elliptic_E_complete_matches_phi_over_two():
    k = 0.5
    expected = _elliptic_E_numeric(0.5 * math.pi, k, n=2000)
    assert math.isclose(clib.geom2d_elliptic_E_complete(k), expected, rel_tol=1e-6, abs_tol=1e-6)


def test_elliptic_E_symmetry_and_edge_cases():
    phi = 0.9
    k = 0.4
    assert math.isclose(clib.geom2d_elliptic_E(-phi, k), -clib.geom2d_elliptic_E(phi, k), rel_tol=1e-9)
    assert clib.geom2d_elliptic_E(0.0, k) == 0.0


def test_elliptic_E_invalid_k_returns_nan():
    assert math.isnan(clib.geom2d_elliptic_E(1.0, -0.1))
    assert math.isnan(clib.geom2d_elliptic_E(1.0, 1.0))
    assert math.isnan(clib.geom2d_elliptic_E_complete(-0.1))
    assert math.isnan(clib.geom2d_elliptic_E_complete(1.0))


def test_merge_sorted():
    a = np.array([1.0, 2.0, 2.0, 4.0], dtype=np.float64)
    b = np.array([2.0, 3.0], dtype=np.float64)
    out_len = np.zeros(1, dtype=np.int32)
    merged, out_len_val = clib.merge_sorted(a, b, out_len=out_len)
    np.testing.assert_allclose(merged, np.array([1.0, 2.0, 3.0, 4.0]))
    assert out_len_val == 4