import numpy as np
from cgeom import Path2D, clib


def test_poly_boundary_distance_bruteforce():
    a = [
        (1, 3),
        (3, 1),
        (3, -1),
        (1, -3),
        (-1, -3),
        (-3, -1),
        (-3, 1),
        (-1, 3),
        (1, 3),
    ]
    b = [(2, 0), (2, -1), (1, -1), (1, 0), (2, 0)]

    points_a = np.array(a, dtype=clib.G2DPoint.dtype)
    points_b = np.array(b, dtype=clib.G2DPoint.dtype)

    result = clib.geom2d_poly_boundary_distance_bruteforce(points_a, points_b)
    expected = np.sqrt(2) / 2

    assert np.isclose(result, expected, rtol=1e-10, atol=1e-10)


def test_poly_boundary_distance_convex():
    a = [
        (1, 3),
        (3, 1),
        (3, -1),
        (1, -3),
        (-1, -3),
        (-3, -1),
        (-3, 1),
        (-1, 3),
        (1, 3),
    ]
    b = [(2, 0), (2, -1), (1, -1), (1, 0), (2, 0)]

    points_a = np.array(a, dtype=clib.G2DPoint.dtype)
    points_b = np.array(b, dtype=clib.G2DPoint.dtype)

    result = clib.geom2d_polygon_boundary_distance_convex_nested(points_a, points_b)
    expected = np.sqrt(2) / 2

    assert np.isclose(result, expected, rtol=1e-10, atol=1e-10)
