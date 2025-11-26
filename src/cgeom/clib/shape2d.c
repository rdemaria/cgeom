#include "base.h"
#include "path.h"
#include <math.h>


void geom2d_circle_get_n_points(double cx, double cy, double r, int len_points, G2DPoint *out_points)
/* Get n points on the circumference of a circle

Contract: len(out_points)=len_points
*/
{
    if (len_points <= 0)
    {
        return;
    }

    double angle_step = 2.0 * M_PI / (len_points - 1);

    for (int i = 0; i < len_points; i++)
    {
        double angle = i * angle_step;
        out_points[i].x = cx + r * cos(angle);
        out_points[i].y = cy + r * sin(angle);
    }
}


void geom2d_circle_get_points(double cx, double cy, double r, G2DPoint *out_points)
/* Get n points on the circumference of a circle between two angles

Contract: len(out_points)=101
*/
{
    geom2d_circle_get_n_points(cx, cy, r, 101, out_points);
}

void geom2d_rect_get_corners(double cx, double cy, double halfw, double halfh, G2DPoint *out_corners)
/* Get the 4 corners of a rectangle

Contract: len(out_corners)=4
*/
{
    out_corners[0].x = cx - halfw;
    out_corners[0].y = cy - halfh;
    out_corners[1].x = cx + halfw;
    out_corners[1].y = cy - halfh;
    out_corners[2].x = cx + halfw;
    out_corners[2].y = cy + halfh;
    out_corners[3].x = cx - halfw;
    out_corners[3].y = cy + halfh;
}

void geom2d_rect_get_points(double cx, double cy, double halfw, double halfh, G2DPoint *out_points)
/* Get the 5 points of a rectangle (last point = first point)

Contract: len(out_points)=5
*/
{
    out_points[0].x = cx - halfw;
    out_points[0].y = cy - halfh;
    out_points[1].x = cx + halfw;
    out_points[1].y = cy - halfh;
    out_points[2].x = cx + halfw;
    out_points[2].y = cy + halfh;
    out_points[3].x = cx - halfw;
    out_points[3].y = cy + halfh;
    out_points[4].x = cx - halfw;
    out_points[4].y = cy - halfh;
}