#ifndef CGEOM_POLYGON_DIST_H
#define CGEOM_POLYGON_DIST_H

#include "base.h"

double geom2d_poly_boundary_distance_bruteforce(
    const G2DPoint* a,
    size_t len_a,
    const G2DPoint* b,
    size_t len_b
);

double geom2d_polygon_boundary_distance_convex_nested(
    G2DPoint *a, size_t len_a,
    G2DPoint *b, size_t len_b
);

#endif
