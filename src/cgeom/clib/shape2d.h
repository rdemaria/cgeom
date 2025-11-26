#ifndef CGEOM_SHAPE2D_H
#define CGEOM_SHAPE2D_H

#include "base.h"
#include "path.h"


void geom2d_circle_get_n_points(double cx, double cy, double r, int len_points, G2DPoint *out_points);
void geom2d_circle_get_points(double cx, double cy, double r, G2DPoint *out_points);
void geom2d_rect_get_corners(double cx, double cy, double halfw, double halfh, G2DPoint *out_corners);
void geom2d_rect_get_points(double cx, double cy, double halfw, double halfh, G2DPoint *out_points);



#endif /* CGEOM_SHAPE2D_H */