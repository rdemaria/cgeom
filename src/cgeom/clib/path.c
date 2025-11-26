#include <math.h>
#include <stdlib.h>

#include "base.h"
#include "path.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* 
geom2d_<objec>_<operation>....

get_n_points : get n points of the profiles
get_points : get sufficient points for a good representation of the profile
get_steps: define steps along the profile that garantee a minimum resolution and the feautures of the profile
*/


void geom2d_line_segment_from_start_length(double x0, double y0, double dx, double dy, double length, G2DSegment *out)
/* Get line data from starting point, direction (assuming dx,dy have norm=1) and length

*/
{
    double ux = dx;
    double uy = dy;
    out->data[0] = x0;
    out->data[1] = y0;
    out->data[2] = x0 + ux * length;
    out->data[3] = y0 + uy * length;
    out->type = 0; /* line */
}

double geom2d_line_segment_get_length(const G2DSegment *seg)
/* Get length of a line segment */
{
    double x1 = seg->data[0];
    double y1 = seg->data[1];
    double x2 = seg->data[2];
    double y2 = seg->data[3];
    return geom2d_norm(x2 - x1, y2 - y1);
}

double geom2d_line_segment_get_points_at_steps(const G2DSegment *seg, const double *steps, int len_points, G2DPoint *out_points)
/* Get points along a line segment at specified steps

Contract: len(steps)=len_points; len(out_points)=len_points
*/
{
    double x1 = seg->data[0];
    double y1 = seg->data[1];
    double x2 = seg->data[2];
    double y2 = seg->data[3];
    double line_length = geom2d_norm(x2 - x1, y2 - y1);
    double ux = (x2 - x1) / line_length;
    double uy = (y2 - y1) / line_length;

    for (int i = 0; i < len_points; i++)
    {
        double at = steps[i];
        out_points[i].x = x1 + ux * at;
        out_points[i].y = y1 + uy * at;
    }
    return line_length;
}



void geom2d_arc_segment_from_start_length_angle(double x0, double y0, double dx, double dy, double length, double angle, G2DSegment *out)
/* Get arc data from starting point, direction (assuming dx,dy have norm=1), length and angle

*/
{
    double norm = geom2d_norm(dx, dy);
    if (angle == 0.0)
    {
        geom2d_line_segment_from_start_length(x0, y0, dx, dy, length, out);
        return;
    }

    double r = length / fabs(angle);
    double cx = x0 - dy * r / norm;
    double cy = y0 + dx * r / norm;
    double start_angle = atan2(y0 - cy, x0 - cx);
    double end_angle = start_angle + angle;
    out->data[0] = cx;
    out->data[1] = cy;
    out->data[2] = r;
    out->data[3] = start_angle;
    out->data[4] = end_angle;
    out->type = 1; /* arc */
    return;
}

void geom2d_arc_segment_get_pd_at_length(const G2DSegment *seg, double at, double *out_x, double *out_y, double *out_dx, double *out_dy)
/* Get point and direction at length 'at' along an arc segment */
{
    /* arc */
    double cx = seg->data[0];
    double cy = seg->data[1];
    double r = seg->data[2];
    double start_angle = seg->data[3];
    double end_angle = seg->data[4];
    double total_angle = end_angle - start_angle;
    double arc_length = fabs(total_angle) * r;
    double angle_at = start_angle + (at / arc_length) * total_angle;
    *out_x = cx + r * cos(angle_at);
    *out_y = cy + r * sin(angle_at);
    *out_dx = -sin(angle_at);
    *out_dy = cos(angle_at);
}

double geom2d_arc_segment_get_length(const G2DSegment *seg)
/* Get length of an arc segment */
{
    double r = seg->data[2];
    double start_angle = seg->data[3];
    double end_angle = seg->data[4];
    double total_angle = end_angle - start_angle;
    return fabs(total_angle) * r;
}

double geom2d_arc_segment_get_points_at_steps(const G2DSegment *seg, const double *steps, int len_points, G2DPoint *out_points)
/* Get points along an arc segment at specified steps

Contract: len(steps)=len_points; len(out_points)=len_points
*/
{
    double cx = seg->data[0];
    double cy = seg->data[1];
    double r = seg->data[2];
    double start_angle = seg->data[3];
    double end_angle = seg->data[4];
    double total_angle = end_angle - start_angle;
    double arc_length = fabs(total_angle) * r;

    for (int i = 0; i < len_points; i++)
    {
        double at = steps[i];
        double angle_at = start_angle + (at / arc_length) * total_angle;
        out_points[i].x = cx + r * cos(angle_at);
        out_points[i].y = cy + r * sin(angle_at);
    }
    return arc_length;
}

static double geom2d_ellipse_cumulative_length(double angle, double rx, double ry)
/* Cumulative length of ellipse from 0 to angle */
{
    if (rx <= 0.0 || ry <= 0.0)
        return 0.0;
    if (rx == ry)
        return rx * fabs(angle);

    if (rx >= ry)
    {
        double k = sqrt(1.0 - (ry * ry) / (rx * rx));
        double complete_E = geom2d_elliptic_E_complete(k);
        return rx * (complete_E - geom2d_elliptic_E(0.5 * M_PI - angle, k));
    }

    double k = sqrt(1.0 - (rx * rx) / (ry * ry));
    return ry * geom2d_elliptic_E(angle, k);
}

static double geom2d_ellipse_arc_length_between(double start, double end, double rx, double ry)
/* Arc length of ellipse between start and end angles */
{
    return fabs(geom2d_ellipse_cumulative_length(end, rx, ry) - geom2d_ellipse_cumulative_length(start, rx, ry));
}

static double geom2d_ellipse_angle_at_length(double start, double end, double rx, double ry, double target)
/* Find ellipse parameter angle at a given arc length from start towards end */
{
    if (rx <= 0.0 || ry <= 0.0)
        return start;
    double total_len = geom2d_ellipse_arc_length_between(start, end, rx, ry);
    if (total_len == 0.0)
        return start;
    if (target <= 0.0)
        return start;
    if (target >= total_len)
        return end;

    double dir = (end >= start) ? 1.0 : -1.0;
    double theta_low = 0.0;
    double theta_high = fabs(end - start);

    if (rx == ry)
        return start + dir * (target / rx);

    for (int i = 0; i < 60; i++)
    {
        double theta_mid = 0.5 * (theta_low + theta_high);
        double angle_mid = start + dir * theta_mid;
        double len_mid = geom2d_ellipse_arc_length_between(start, angle_mid, rx, ry);
        if (len_mid < target)
            theta_low = theta_mid;
        else
            theta_high = theta_mid;
    }
    double theta_mid = 0.5 * (theta_low + theta_high);
    return start + dir * theta_mid;
}

double geom2d_ellipse_segment_get_length(const G2DSegment *seg)
/* Get length of an ellipse arc segment */
{
    double rx = seg->data[2];
    double ry = seg->data[3];
    double start_angle = seg->data[5];
    double end_angle = seg->data[6];
    return geom2d_ellipse_arc_length_between(start_angle, end_angle, rx, ry);
}

void geom2d_ellipse_segment_get_points_at_steps(const G2DSegment *seg, const double *steps, int len_points, G2DPoint *out_points)
/* Get points along an ellipse arc segment at specified steps

Contract: len(steps)=len_points; len(out_points)=len_points
*/
{
    double cx = seg->data[0];
    double cy = seg->data[1];
    double rx = seg->data[2];
    double ry = seg->data[3];
    double rotation = seg->data[4];
    double start_angle = seg->data[5];
    double end_angle = seg->data[6];

    double cos_rot = cos(rotation);
    double sin_rot = sin(rotation);

    for (int i = 0; i < len_points; i++)
    {
        double at = steps[i];
        double angle_at = geom2d_ellipse_angle_at_length(start_angle, end_angle, rx, ry, at);
        double x_ellipse = rx * cos(angle_at);
        double y_ellipse = ry * sin(angle_at);
        // Apply rotation
        out_points[i].x = cx + (x_ellipse * cos_rot - y_ellipse * sin_rot);
        out_points[i].y = cy + (x_ellipse * sin_rot + y_ellipse * cos_rot);
    }
}

double geom2d_segment_get_length(const G2DSegment *seg)
/* Get length of a segment */
{
    switch (seg->type)
    {
    case 0: /* line */
        return geom2d_line_segment_get_length(seg);
    case 1: /* arc */
        return geom2d_arc_segment_get_length(seg);
    case 2: /* ellipse arc */
        return geom2d_ellipse_segment_get_length(seg);
    default:
        return 0.0;
    }
}


void geom2d_rectangle_to_path(double halfwidth, double halfheight, G2DSegment *out_segments)
/* Create a path for a rectangle centered at (0,0)

Contract: len(out_segments)=4
*/
{
    out_segments[0].type = 0; /* line */
    out_segments[0].data[0] = -halfwidth;
    out_segments[0].data[1] = -halfheight;
    out_segments[0].data[2] = halfwidth;
    out_segments[0].data[3] = -halfheight;

    out_segments[1].type = 0; /* line */
    out_segments[1].data[0] = halfwidth;
    out_segments[1].data[1] = -halfheight;
    out_segments[1].data[2] = halfwidth;
    out_segments[1].data[3] = halfheight;

    out_segments[2].type = 0; /* line */
    out_segments[2].data[0] = halfwidth;
    out_segments[2].data[1] = halfheight;
    out_segments[2].data[2] = -halfwidth;
    out_segments[2].data[3] = halfheight;

    out_segments[3].type = 0; /* line */
    out_segments[3].data[0] = -halfwidth;
    out_segments[3].data[1] = halfheight;
    out_segments[3].data[2] = -halfwidth;
    out_segments[3].data[3] = -halfheight;
}

void geom2d_ellipse_to_path(double rx, double ry, G2DSegment *out_segment)
/* Create a path for an ellipse centered at (0,0)

Contract: len(out_segments) must be at least 1
*/
{
    out_segment[0].type = 3;      /* arc */
    out_segment[0].data[0] = 0.0; /* cx */
    out_segment[0].data[1] = 0.0; /* cy */
    out_segment[0].data[2] = rx;  /* rx */
    out_segment[0].data[3] = ry;  /* ry */
    out_segment[0].data[4] = 0.0; /* rotation */
    out_segment[0].data[5] = 0.0; /* start_angle */
    out_segment[0].data[6] = 2.0 * M_PI; /* end_angle */
}

void geom2d_rectellipse_to_path(double halfwidth, double halfheight, double rx, double ry, G2DSegment *out_segments, int *out_len)
/* Create a path for the intersection between a reactangle and an ellipse

Contract: len(out_segments) must be at least 8
Post-contract: len(out_segments)=out_len
*/
{
    double arg1 = halfwidth / rx;
    double arg2 = halfheight / ry;

    if (arg1 > 1.0) arg1 = 1.0;
    if (arg2 > 1.0) arg2 = 1.0;

    if (arg1 == 1.0 && arg2 == 1.0) //ellipse inside rectangle
    {  geom2d_ellipse_to_path(rx, ry, out_segments);
       *out_len = 1;
       return;
    }

    double angle1 = acos(halfwidth / rx);
    double angle2 = asin(halfheight / ry);

    if (angle2 < angle1)
    { // rectangle is inside ellipse
        geom2d_rectangle_to_path(halfwidth, halfheight, out_segments);
        *out_len = 4;
        return;
    }

    double iy=ry*sin(angle1);
    double ix=rx*cos(angle2);

    if (arg1 < 1 && arg2 < 1){ /* 8 segments */
        out_segments[0].type = 0; /* line */
        out_segments[0].data[0] = halfwidth;
        out_segments[0].data[1] = -iy;
        out_segments[0].data[2] = halfwidth;
        out_segments[0].data[3] = iy;
        out_segments[1].type = 3;      /* arc */
        out_segments[1].data[0] = 0.0; /* cx */
        out_segments[1].data[1] = 0.0; /* cy */
        out_segments[1].data[2] = rx;  /* rx */
        out_segments[1].data[3] = ry;  /* ry */
        out_segments[1].data[4] = 0.0; /* rotation */
        out_segments[1].data[5] = angle1;
        out_segments[1].data[6] = angle2;
        out_segments[2].type = 0; /* line */
        out_segments[2].data[0] = ix;
        out_segments[2].data[1] = halfheight;
        out_segments[2].data[2] = -ix;
        out_segments[2].data[3] = halfheight;
        out_segments[3].type = 3;      /* arc */
        out_segments[3].data[0] = 0.0; /* cx */
        out_segments[3].data[1] = 0.0; /* cy */
        out_segments[3].data[2] = rx;  /* rx */
        out_segments[3].data[3] = ry;  /* ry */
        out_segments[3].data[4] = 0.0; /* rotation */
        out_segments[3].data[5] = M_PI/2 + angle1;
        out_segments[3].data[6] = M_PI/2 + angle2;
        out_segments[4].type = 0; /* line */
        out_segments[4].data[0] = -halfwidth;
        out_segments[4].data[1] = iy;
        out_segments[4].data[2] = -halfwidth;
        out_segments[4].data[3] = -iy;
        out_segments[5].type = 3;      /* arc */
        out_segments[5].data[0] = 0.0; /* cx */
        out_segments[5].data[1] = 0.0; /* cy */
        out_segments[5].data[2] = rx;  /* rx */
        out_segments[5].data[3] = ry;  /* ry */
        out_segments[5].data[4] = 0.0; /* rotation */
        out_segments[5].data[5] = M_PI + angle1;
        out_segments[5].data[6] = M_PI + angle2;
        out_segments[6].type = 0; /* line */
        out_segments[6].data[0] = -ix;
        out_segments[6].data[1] = -halfheight;
        out_segments[6].data[2] = ix;
        out_segments[6].data[3] = -halfheight;
        out_segments[7].type = 3;      /* arc */
        out_segments[7].data[0] = 0.0; /* cx */
        out_segments[7].data[1] = 0.0; /* cy */
        out_segments[7].data[2] = rx;  /* rx */
        out_segments[7].data[3] = ry;  /* ry */
        out_segments[7].data[4] = 0.0; /* rotation */
        out_segments[7].data[5] = 3*M_PI/2 + angle1;
        out_segments[7].data[6] = 3*M_PI/2 + angle2;
        *out_len = 8;
        return;
    }
    if (angle1 <1){ /*  arc top and bottom */
        out_segments[0].type = 0; /* line */
        out_segments[0].data[0] = halfwidth;
        out_segments[0].data[1] = -iy;
        out_segments[0].data[2] = halfwidth;
        out_segments[0].data[3] = iy;
        out_segments[1].type = 3; /* arc */
        out_segments[1].data[0] = 0.0; /* cx */
        out_segments[1].data[1] = 0.0; /* cy */
        out_segments[1].data[2] = rx;  /* rx */
        out_segments[1].data[3] = ry;  /* ry */
        out_segments[1].data[4] = 0.0; /* rotation */
        out_segments[1].data[5] = angle1;
        out_segments[1].data[6] = M_PI - angle1;
        out_segments[2].type = 0; /* line */
        out_segments[2].data[0] = -halfwidth;
        out_segments[2].data[1] = iy;
        out_segments[2].data[2] = -halfwidth;;
        out_segments[2].data[3] = -iy;
        out_segments[3].type = 3; /* arc */
        out_segments[3].data[0] = 0.0; /* cx */
        out_segments[3].data[1] = 0.0; /* cy */
        out_segments[3].data[2] = rx;  /* rx */
        out_segments[3].data[3] = ry;  /* ry */
        out_segments[3].data[4] = 0.0; /* rotation */
        out_segments[3].data[5] = M_PI + angle1;;
        out_segments[3].data[6] = 2*M_PI - angle1;
        *out_len = 4;
        return;
    }
    if (angle2 <1){ /* arc right and left */
        out_segments[0].type = 3; /* arc */
        out_segments[0].data[0] = 0.0; /* cx */
        out_segments[0].data[1] = 0.0; /* cy */
        out_segments[0].data[2] = rx;  /* rx */
        out_segments[0].data[3] = ry;  /* ry */
        out_segments[0].data[4] = 0.0; /* rotation */
        out_segments[0].data[5] = -angle2;
        out_segments[0].data[6] = angle2;
        out_segments[1].type = 0; /* line */
        out_segments[1].data[0] = ix;
        out_segments[1].data[1] = halfheight;
        out_segments[1].data[2] = -ix;
        out_segments[1].data[3] = halfheight;
        out_segments[2].type = 3; /* arc */
        out_segments[2].data[0] = 0.0; /* cx */
        out_segments[2].data[1] = 0.0; /* cy */
        out_segments[2].data[2] = rx;  /* rx */
        out_segments[2].data[3] = ry;  /* ry */
        out_segments[2].data[4] = 0.0; /* rotation */
        out_segments[2].data[5] = M_PI - angle2;
        out_segments[2].data[6] = M_PI + angle2;
        out_segments[3].type = 0; /* line */
        out_segments[3].data[0] = -ix;
        out_segments[3].data[1] = -halfheight;
        out_segments[3].data[2] = ix;
        out_segments[3].data[3] = -halfheight;
        *out_len = 4;
        return;
    }
}

int geom2d_path_get_len_steps(const G2DSegment *segments, int len_segments, double ds_min){
/* Get the number of steps needed to represent a path defined by segments

Contract: len_segments=len(segments)
*/
    int total_steps = 1;
    int nsteps;
    for (int i = 0; i < len_segments; i++)
    {
        switch (segments[i].type)
        {
        case 0: /* line */
            total_steps += ceil(geom2d_line_segment_get_length(&segments[i])/ds_min);
            break;
        case 1: /* arc */
            nsteps = ceil(geom2d_arc_segment_get_length(&segments[i])/ds_min);
            total_steps += nsteps>11 ? nsteps : 11;
            break;
        case 2: /* ellipse arc */
            nsteps = ceil(geom2d_ellipse_segment_get_length(&segments[i])/ds_min);
            total_steps += nsteps>11 ? nsteps : 11;
            break;
        default:
            break;
        }
    }
    return total_steps;
}

double geom2d_path_get_length(const G2DSegment *segments, int len_segments){
/* Get length of a path defined by segments

Contract: len_segments=len(segments)
*/
    double total_length = 0.0;
    for (int i = 0; i < len_segments; i++)
    {
        total_length += geom2d_segment_get_length(&segments[i]);
    }
    return total_length;
}

void geom2d_path_get_steps(const G2DSegment *segments, int len_segments, double ds_min, double *out_steps){
/* Get steps along a path defined by segments

Contract: len_segments=len(segments); len(out_steps)=geom2d_path_get_len_steps(segments,len_segments,ds_min)
*/
    int idx = 0;
    out_steps[idx++] = 0.0;
    double seg_length;
    double nsteps_d;
    int nsteps;
    double ds;
    double length_acc = 0.0;

    for (int i = 0; i < len_segments; i++)
    {
        switch (segments[i].type)
        {
        case 0: /* line */
            seg_length = geom2d_line_segment_get_length(&segments[i]);
            nsteps_d = seg_length / ds_min;
            nsteps = (int)ceil(nsteps_d);
            ds = seg_length / nsteps;
            for (int j = 1; j <= nsteps; j++)
            {
                length_acc += ds;
                out_steps[idx++] = length_acc;
            }
            break;
        case 1: /* arc */
            seg_length = geom2d_arc_segment_get_length(&segments[i]);
            nsteps_d = seg_length / ds_min;
            nsteps = (int)ceil(nsteps_d);
            if (nsteps < 11)
                nsteps = 11;
            ds = seg_length / nsteps;
            for (int j = 1; j <= nsteps; j++)
            {
                length_acc += ds;
                out_steps[idx++] = length_acc;
            }
            break;
        case 2: /* ellipse arc */
            seg_length = geom2d_ellipse_segment_get_length(&segments[i]);
            nsteps_d = seg_length / ds_min;
            nsteps = (int)ceil(nsteps_d);
            if (nsteps < 11)
                nsteps = 11;
            ds = seg_length / nsteps;
            for (int j = 1; j <= nsteps; j++)
            {
                length_acc += ds;
                out_steps[idx++] = length_acc;
            }
            break;
        default:
            break;
        }
    }
}

void geom2d_path_get_points_at_steps(const G2DSegment *segments, int len_segments, const double *steps, int len_points, G2DPoint *out_points)
/* Get points along a path defined by segments at specified steps

Contract: len_segments=len(segments); len(steps)=len_points; len(out_points)=len_points
*/
{
    int seg_idx = 0;
    double seg_start_length = 0.0;
    double seg_end_length = 0.0;
    double seg_length;
    for (int i = 0; i < len_segments; i++)
    {
        seg_length = geom2d_segment_get_length(&segments[i]);
        seg_end_length = seg_start_length + seg_length;

        while (seg_idx < len_points && steps[seg_idx] <= seg_end_length)
        {
            double at = steps[seg_idx] - seg_start_length;
            switch (segments[i].type)
            {
            case 0: /* line */
                geom2d_line_segment_get_points_at_steps(&segments[i], &at, 1, &out_points[seg_idx]);
                break;
            case 1: /* arc */
                geom2d_arc_segment_get_points_at_steps(&segments[i], &at, 1, &out_points[seg_idx]);
                break;
            case 2: /* ellipse arc */
                // Implement ellipse arc point retrieval if needed
                break;
            default:
                break;
            }
            seg_idx++;
        }
        seg_start_length = seg_end_length;
    }
}
