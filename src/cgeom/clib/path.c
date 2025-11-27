#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "base.h"
#include "path.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
geom2d_<objec>_<operation>....

*/

/* ===== Line segment functions ===== */

void geom2d_line_segment_from_start_end(double x0, double y0, double x1, double y1, G2DSegment *out)
/* Get line data from starting and ending points

*/
{
    out->data[0] = x0;
    out->data[1] = y0;
    out->data[2] = x1;
    out->data[3] = y1;
    out->type = 0; /* line */
}

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

void geom2d_line_segment_get_points_at_steps(const G2DSegment *seg, const double *steps, int len_points, G2DPoint *out_points)
/* Get points along a line segment at specified steps

Contract: len_points=len(steps); len(out_points)=len_points
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
}

/* ===== End line segment functions ===== */

/* ===== Arc segment functions ===== */

void geom2d_arc_segment_from_center_radius_angles(double cx, double cy, double r, double start_angle, double end_angle, G2DSegment *out)
/* Get arc data from center, radius and start/end angles

*/
{
    out->data[0] = cx;
    out->data[1] = cy;
    out->data[2] = r;
    out->data[3] = start_angle;
    out->data[4] = end_angle;
    out->type = 1; /* arc */
}

void geom2d_arc_segment_from_ref_length_angle(double x0, double y0, double dx, double dy, double length, double angle, G2DSegment *out)
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

void geom2d_arc_segment_get_ref_at_length(const G2DSegment *seg, double at, double *out_x, double *out_y, double *out_dx, double *out_dy)
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

void geom2d_arc_segment_get_points_at_steps(const G2DSegment *seg, const double *steps, int len_points, G2DPoint *out_points)
/* Get points along an arc segment at specified steps

Contract: len_points=len(steps); len(out_points)=len_points
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
}

/* ===== End arc segment functions ===== */

/* ===== Ellipse arc segment functions ===== */

void geom2d_ellipse_arc_segment_from_center_radii_rotation_angles(double cx, double cy, double rx, double ry, double rotation, double start_angle, double end_angle, G2DSegment *out)
{
    out->data[0] = cx;
    out->data[1] = cy;
    out->data[2] = rx;
    out->data[3] = ry;
    out->data[4] = rotation;
    out->data[5] = start_angle;
    out->data[6] = end_angle;
    out->type = CGEOM_ELLIPSE_ARC_SEGMENT_TYPE;
}

void geom2d_maybe_ellipse_arc_segment_from_center_radii_rotation_angles(double cx, double cy, double rx, double ry, double rotation, double start_angle, double end_angle, G2DSegment *out)
{
    if (rx == ry)
    {
        double r = rx;
        geom2d_arc_segment_from_center_radius_angles(cx, cy, r, start_angle, end_angle, out);
        return;
    }
    geom2d_ellipse_arc_segment_from_center_radii_rotation_angles(cx, cy, rx, ry, rotation, start_angle, end_angle, out);
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

/* ===== End ellipse arc segment functions ===== */

/* ===== Racetrack segment functions ===== */
void geom2d_segments_from_racetrack(double halfwidth, double halfheight, double rx, double ry, G2DSegment *out_segments, int *out_len)
/* Create a path for a racetrack shape centered at (0,0)

Contract: len(out_segments)=8
Post-Contract: len(out_segments)=out_len
*/
{
    *out_len = 0;
    // Straight right
    if (halfheight > 0){
        geom2d_line_segment_from_start_end(halfwidth, -halfheight+ry, halfwidth, halfheight-ry, &out_segments[0]);
        (*out_len)++;
    }
    if (rx>0 && ry>0){
         // Top-right arc
         geom2d_maybe_ellipse_arc_segment_from_center_radii_rotation_angles(halfwidth-rx, halfheight-ry, rx, ry, 0.0, 0, M_PI/2, &out_segments[*out_len]);
         (*out_len)++;
    }
    if (halfwidth > 0){
        // Straight top
        geom2d_line_segment_from_start_end(halfwidth - rx, halfheight , -halfwidth + rx, halfheight, &out_segments[*out_len]);
        (*out_len)++;
    }
    if (rx>0 && ry>0){
        // Top-left arc
        geom2d_maybe_ellipse_arc_segment_from_center_radii_rotation_angles(-halfwidth+rx, halfheight-ry, rx, ry, 0.0, M_PI/2, M_PI, &out_segments[*out_len]);
        (*out_len)++;
    }
    if (halfheight > 0){
        // Straight left
        geom2d_line_segment_from_start_end(-halfwidth, halfheight - ry, -halfwidth, -halfheight + ry, &out_segments[*out_len]);
        (*out_len)++;
    }
    if (rx>0 && ry>0){
        // Bottom-left arc
        geom2d_maybe_ellipse_arc_segment_from_center_radii_rotation_angles(-halfwidth+rx, -halfheight+ry, rx, ry, 0.0, M_PI, 3.0*M_PI/2.0, &out_segments[*out_len]);
        (*out_len)++;
    }
    if (halfheight > 0){
        // Straight bottom
        geom2d_line_segment_from_start_end(-halfwidth + rx, -halfheight, halfwidth - rx, -halfheight, &out_segments[*out_len]);
        (*out_len)++;
    }
    if (rx>0 && ry>0){
        // Bottom-right arc
        geom2d_maybe_ellipse_arc_segment_from_center_radii_rotation_angles(halfwidth-rx, -halfheight+ry, rx, ry, 0.0, 3.0*M_PI/2.0, 2.0*M_PI, &out_segments[*out_len]);
        (*out_len)++;
    }
}



/* ===== End racetrack segment functions ===== */


/* ===== Segment functions ===== */

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

/* ===== End segment functions ===== */

/* ===== Segments from shapes ===== */

void geom2d_segments_from_rectangle(double halfwidth, double halfheight, G2DSegment *out_segments)
/* Create a path for a rectangle centered at (0,0)

Contract: len(out_segments)=4
*/
{
    geom2d_line_segment_from_start_end(-halfwidth, -halfheight, halfwidth, -halfheight, &out_segments[0]);
    geom2d_line_segment_from_start_end(halfwidth, -halfheight, halfwidth, halfheight, &out_segments[1]);
    geom2d_line_segment_from_start_end(halfwidth, halfheight, -halfwidth, halfheight, &out_segments[2]);
    geom2d_line_segment_from_start_end(-halfwidth, halfheight, -halfwidth, -halfheight, &out_segments[3]);
}

void geom2d_segments_from_circle(double r, G2DSegment *out_segments)
/* Create a path for a circle centered at (0,0)

Contract: len(out_segments)=1
*/
{
    geom2d_arc_segment_from_center_radius_angles(0.0, 0.0, r, 0.0, 2.0 * M_PI, &out_segments[0]);
}

void geom2d_segments_from_ellipse(double rx, double ry, G2DSegment *out_segments)
/* Create a path for an ellipse centered at (0,0)

Contract: len(out_segments)=1
*/
{
    if (rx == ry)
    {
        geom2d_segments_from_circle(rx, out_segments);
        return;
    }
    geom2d_ellipse_arc_segment_from_center_radii_rotation_angles(0.0, 0.0, rx, ry, 0.0, 0.0, 2.0 * M_PI, &out_segments[0]);
}

void geom2d_segments_from_rectellipse(double halfwidth, double halfheight, double rx, double ry, G2DSegment *out_segments, int *out_len)
/* Create a path for the intersection between a rectangle and an ellipse

Contract: len(out_segments)=8
Post-contract: len(out_segments)=out_len
*/
{
    double arg1 = halfwidth / rx;
    double arg2 = halfheight / ry;
 
    if (arg1 > 1.0)
        arg1 = 1.0;
    if (arg2 > 1.0)
        arg2 = 1.0;

    if (arg1 == 1.0 && arg2 == 1.0)
    {
        // printf("ellipse inside rectangle\n");
        geom2d_segments_from_ellipse(rx, ry, out_segments);
        *out_len = 1;
        return;
    }

    double angle1 = acos(halfwidth / rx);
    double angle2 = asin(halfheight / ry);
    // printf("angle1=%f angle2=%f\n",angle1,angle2);

    if (angle2 < angle1)
    {
        // printf("rectangle inside ellipse\n");
        geom2d_segments_from_rectangle(halfwidth, halfheight, out_segments);
        *out_len = 4;
        return;
    }

    double iy = ry * sin(angle1);
    double ix = rx * cos(angle2);

    if (arg1 < 1 && arg2 < 1)
    { // printf("8 segments\n");
        geom2d_line_segment_from_start_end(halfwidth, -iy, halfwidth, iy, &out_segments[0]);
        geom2d_maybe_ellipse_segment_from_center_radius_angles(0.0, 0.0, rx, ry, 0.0, angle1, angle2, &out_segments[1]);
        geom2d_line_segment_from_start_end(ix, halfheight, -ix, halfheight, &out_segments[2]);
        geom2d_maybe_ellipse_segment_from_center_radius_angles(0.0, 0.0, rx, ry, M_PI / 2 + angle1, M_PI / 2 + angle2, &out_segments[3]);
        geom2d_line_segment_from_start_end(-halfwidth, iy, -halfwidth, -iy, &out_segments[4]);
        geom2d_maybe_ellipse_segment_from_center_radius_angles(0.0, 0.0, rx, ry, M_PI + angle1, M_PI + angle2, &out_segments[5]);
        geom2d_line_segment_from_start_end(-ix, -halfheight, ix, -halfheight, &out_segments[6]);
        geom2d_maybe_ellipse_segment_from_center_radius_angles(0.0, 0.0, rx, ry, 3 * M_PI / 2 + angle1, 3 * M_PI / 2 + angle2, &out_segments[7]);
        *out_len = 8;
        return;
    }
    if (arg1 < 1)
    { // printf("flat sides\n");
        geom2d_line_segment_from_start_end(halfwidth, -iy, halfwidth, iy, &out_segments[0]);
        geom2d_maybe_ellipse_segment_from_center_radius_angles(0.0, 0.0, rx, ry, angle1, M_PI - angle1, &out_segments[1]);
        geom2d_line_segment_from_start_end(-halfwidth, iy, -halfwidth, -iy, &out_segments[2]);
        geom2d_maybe_ellipse_segment_from_center_radius_angles(0.0, 0.0, rx, ry, M_PI + angle1, 2 * M_PI - angle1, &out_segments[3]);
        *out_len = 4;
        return;
    }
    if (arg2 < 1)
    { // printf("flat top/bottom\n");
        geom2d_maybe_ellipse_segment_from_center_radius_angles(0.0, 0.0, rx, ry, -angle2, angle2, &out_segments[0]);
        geom2d_line_segment_from_start_end(ix, halfheight, -ix, halfheight, &out_segments[1]);
        geom2d_maybe_ellipse_segment_from_center_radius_angles(0.0, 0.0, rx, ry, M_PI - angle2, M_PI + angle2, &out_segments[2]);
        geom2d_line_segment_from_start_end(-ix, -halfheight, ix, -halfheight, &out_segments[3]);
        *out_len = 4;
        return;
    }
}

/* ===== End segments from shapes ===== */

/* ===== Path functions ===== */

int geom2d_path_get_len_steps(const G2DPath *path, double ds_min)
{
    /* Get the number of steps needed to represent a path defined by segments

    Contract: len_segments=len(segments)
    */
    int total_steps = 1;
    int nsteps;
    for (int i = 0; i < path->len_segments; i++)
    {
        switch (path->segments[i].type)
        {
        case 0: /* line */
            total_steps += ceil(geom2d_line_segment_get_length(&path->segments[i]) / ds_min);
            break;
        case 1: /* arc */
            nsteps = ceil(geom2d_arc_segment_get_length(&path->segments[i]) / ds_min);
            total_steps += nsteps > 10 ? nsteps : 10;
            break;
        case 2: /* ellipse arc */
            nsteps = ceil(geom2d_ellipse_segment_get_length(&path->segments[i]) / ds_min);
            total_steps += nsteps > 10 ? nsteps : 10;
            break;
        default:
            break;
        }
    }
    return total_steps;
}

double geom2d_path_get_length(const G2DPath *path)
{
    /* Get length of a path defined by segments

    Contract: len_segments=len(segments)
    */
    double total_length = 0.0;
    for (int i = 0; i < path->len_segments; i++)
    {
        total_length += geom2d_segment_get_length(&path->segments[i]);
    }
    return total_length;
}

void geom2d_path_get_steps(const G2DPath *path, double ds_min, double *out_steps)
{
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

    for (int i = 0; i < path->len_segments; i++)
    {
        double seg_start = length_acc;
        switch (path->segments[i].type)
        {
        case 0: /* line */
            seg_length = geom2d_line_segment_get_length(&path->segments[i]);
            nsteps_d = seg_length / ds_min;
            nsteps = (int)ceil(nsteps_d);
            ds = seg_length / nsteps;
            for (int j = 1; j <= nsteps; j++)
            {
                if (j == nsteps)
                    length_acc = seg_start + seg_length;
                else
                    length_acc += ds;
                out_steps[idx++] = length_acc;
            }
            break;
        case 1: /* arc */
            seg_length = geom2d_arc_segment_get_length(&path->segments[i]);
            nsteps_d = seg_length / ds_min;
            nsteps = (int)ceil(nsteps_d);
            if (nsteps < 10)
                nsteps = 10;
            ds = seg_length / nsteps;
            for (int j = 1; j <= nsteps; j++)
            {
                if (j == nsteps)
                    length_acc = seg_start + seg_length;
                else
                    length_acc += ds;
                out_steps[idx++] = length_acc;
            }
            break;
        case 2: /* ellipse arc */
            seg_length = geom2d_ellipse_segment_get_length(&path->segments[i]);
            nsteps_d = seg_length / ds_min;
            nsteps = (int)ceil(nsteps_d);
            if (nsteps < 10)
                nsteps = 10;
            ds = seg_length / nsteps;
            for (int j = 1; j <= nsteps; j++)
            {
                if (j == nsteps)
                    length_acc = seg_start + seg_length;
                else
                    length_acc += ds;
                out_steps[idx++] = length_acc;
            }
            break;
        default:
            break;
        }
    }; // Ensure last step is exact length
}

void geom2d_path_get_points_at_steps(const G2DPath *path, const double *steps, int len_points, G2DPoint *out_points)
/* Get points along a path defined by segments at specified steps

Contract: len_points=len(steps); len(out_points)=len_points
*/
{
    int seg_idx = 0;
    double seg_start_length = 0.0;
    double seg_end_length = 0.0;
    double seg_length;
    for (int i = 0; i < path->len_segments; i++)
    {
        seg_length = geom2d_segment_get_length(&path->segments[i]);
        seg_end_length = seg_start_length + seg_length;

        while (seg_idx < len_points && steps[seg_idx] <= seg_end_length)
        {
            double at = steps[seg_idx] - seg_start_length;
            switch (path->segments[i].type)
            {
            case 0: /* line */
                geom2d_line_segment_get_points_at_steps(&path->segments[i], &at, 1, &out_points[seg_idx]);
                break;
            case 1: /* arc */
                geom2d_arc_segment_get_points_at_steps(&path->segments[i], &at, 1, &out_points[seg_idx]);
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

int geom2d_path_get_len_points(const G2DPath *path)
/*Return the number of points needed for a good representation of the path

*/
{
    int total_points = 1;
    for (int i = 0; i < path->len_segments; i++)
    {
        if (path->segments[i].type == 0)
        {
            total_points += 1;
        }
        else if (path->segments[i].type == 1)
        {
            total_points += 10;
        }
    }
    return total_points;
}

void geom2d_path_get_points(const G2DPath *path, G2DPoint *out_points)
/* Get points along a path defined by segments for a good representation of the path

Contract: len(out_points)=geom2d_path_get_len_points(path)
*/
{
    int idx = 0;
    out_points[idx++].x = path->segments[0].data[0];
    out_points[idx++].y = path->segments[0].data[1];
    for (int i = 0; i < path->len_segments; i++)
    {
        if (path->segments[i].type == 0)
        {
            geom2d_line_segment_get_points_at_steps(&path->segments[i], NULL, 1, &out_points[idx]);
            idx += 1;
        }
        else if (path->segments[i].type == 1)
        {
            double seg_length = geom2d_arc_segment_get_length(&path->segments[i]);
            double ds = seg_length / 10.0;
            double steps[10];
            for (int j = 1; j <= 10; j++)
            {
                steps[j - 1] = j * ds;
            }
            geom2d_arc_segment_get_points_at_steps(&path->segments[i], steps, 10, &out_points[idx]);
            idx += 10;
        }
    }
}

int geom2d_path_get_len_corners(const G2DPath *path)
{
    /*Return number of corners in path

    */
    return path->len_segments + 1;
}

void geom2d_path_get_corner_steps(const G2DPath *path, double *out_steps)
{
    /*Get steps at corners of path

    Contract: len(out_steps)=geom2d_path_get_len_corners(path)
    */
    double length_acc = 0.0;
    out_steps[0] = 0.0;
    for (int i = 0; i < path->len_segments; i++)
    {
        length_acc += geom2d_segment_get_length(&path->segments[i]);
        out_steps[i + 1] = length_acc;
    }
}

/* ===== End path functions ===== */
