import numpy as np
from cgeom import Path2D, clib


def test_is_point_inside_polygon_ellipse():
    rx = 2
    ry = 3
    ellipse = Path2D.from_ellipse(rx=2, ry=3)
    ellipse_poly = ellipse.get_points(ds_min=0.01)

    @np.vectorize
    def in_ellipse(x, y):
        point = clib.G2DPoint(x=x, y=y)
        return clib.geom2d_is_point_inside_polygon(point, ellipse_poly)

    extent = np.linspace(-10, 10, 100)
    xs, ys = np.meshgrid(extent, extent)

    result = in_ellipse(xs, ys).view(bool)
    expected = (xs ** 2 / rx ** 2 + ys ** 2 / ry ** 2 - 1) < 0

    assert not np.all(result) and np.any(result)  # sanity check
    assert np.all(result == expected)


def test_is_point_inside_polygon_path():
    # Define a shape that is a rectangle spanning (-1, -1) through (3, 2) minus
    # a rectangle (1, 0.5) through (2, 2)

    poly = [
        (1, 2),
        (1, .5),
        (2, .5),
        (2, 2),
        (3, 2),
        (3, -1),
        (-1, -1),
        (-1, 2),
        (1, 2),
    ]
    points = np.array(poly, dtype=clib.G2DPoint.dtype)

    @np.vectorize
    def in_poly(x, y):
        point = clib.G2DPoint(x=x, y=y)
        return clib.geom2d_is_point_inside_polygon(point, points)

    extent = np.linspace(-5, 5, 100)
    xs, ys = np.meshgrid(extent, extent)

    result = in_poly(xs, ys).view(bool)

    in_rec1 = (-1 < xs) & (xs < 3) & (-1 < ys) & (ys < 2)
    in_rec2 = (1 < xs) & (xs < 2) & (0.5 < ys) & (ys < 2)
    expected = in_rec1 & ~in_rec2

    assert not np.all(result) and np.any(result)  # sanity check
    assert np.all(result == expected)
