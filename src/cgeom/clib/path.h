#ifndef CGEOM_PATH_H
#define CGEOM_PATH_H

#include "base.h"

/* 2D PATH made of segments  
type: 0=line, 1=arc, 2=ellipse arc, 3=quadratic bezier, 4=cubic bezier
data:
 line: x1,y1,x2,y2
 arc: cx,cy,start_angle,end_angle
 ellipse arc: cx,cy,rx,ry,rotation,start_angle,end_angle
 quadratic bezier: x1,y1,x2,y2,cx,cy
 cubic bezier: x1,y1,x2,y2,cx1,cy1,cx2,cy2

ISSUES:
 - no bezier segments yet
 - no splines yet
 - path assumed to be continuous but not enforced by the structure


*/

#define CGEOM_LINE_SEGMENT_TYPE 0
#define CGEOM_ARC_SEGMENT_TYPE 1
#define CGEOM_ELLIPSE_ARC_SEGMENT_TYPE 2
#define CGEOM_QUADRATIC_BEZIER_SEGMENT_TYPE 3
#define CGEOM_CUBIC_BEZIER_SEGMENT_TYPE 4


typedef struct {
    int type; /* 0=line, 1=arc, 2=ellipse arc */
    double data[8]; 
} G2DSegment;

typedef struct {
    G2DSegment *segments;
    int len_segments;
} G2DPath;

/* Line segment functions */
void geom2d_line_segment_from_start_end(double x0, double y0, double x1, double y1, G2DSegment *out);
void geom2d_line_segment_from_start_length(double x0, double y0, double dx, double dy, double length, G2DSegment *out);
double geom2d_line_segment_get_length(const G2DSegment *seg);
void geom2d_line_segment_get_points_at_steps(const G2DSegment *seg, const double *steps, int len_points, G2DPoint *out_points);

/* Arc segment functions */
void geom2d_arc_segment_from_center_radius_angles(double cx, double cy, double r, double start_angle, double end_angle, G2DSegment *out);
void geom2d_arc_segment_from_ref_length_angle(double x0, double y0, double dx, double dy, double length, double angle, G2DSegment *out);
void geom2d_arc_segment_get_ref_at_length(const G2DSegment *seg, double at, double *out_x, double *out_y, double *out_dx, double *out_dy);
double geom2d_arc_segment_get_length(const G2DSegment *seg);
void geom2d_arc_segment_get_points_at_steps(const G2DSegment *seg, const double *steps, int len_points, G2DPoint *out_points);

/* Ellipse arc segment functions */
void geom2d_ellipse_arc_segment_from_center_radii_rotation_angles(double cx, double cy, double rx, double ry, double rotation, double start_angle, double end_angle, G2DSegment *out);
void geom2d_maybe_ellipse_arc_segment_from_center_radii_rotation_angles(double cx, double cy, double rx, double ry, double rotation, double start_angle, double end_angle, G2DSegment *out);
double geom2d_ellipse_segment_get_length(const G2DSegment *seg);
void geom2d_ellipse_segment_get_points_at_steps(const G2DSegment *seg, const double *steps, int len_points, G2DPoint *out_points);

/* Segment functions */
double geom2d_segment_get_length(const G2DSegment *seg);

/* Segments from shapes */
void geom2d_segments_from_rectangle(double halfwidth, double halfheight, G2DSegment *out_segments);
void geom2d_segments_from_circle(double radius, G2DSegment *out_segments);
void geom2d_segments_from_ellipse(double rx, double ry, G2DSegment *out_segments);
void geom2d_segments_from_rectellipse(double halfwidth, double halfheight, double rx, double ry, G2DSegment *out_segments, int *len_segments);
void geom2d_segments_from_racetrack(double halfhside, double halfvside, double rx, double ry, G2DSegment *out_segments, int *out_len);

/* Path functions */
int geom2d_path_get_len_steps(const G2DPath *path, double ds_min);
double geom2d_path_get_length(const G2DPath *path);
void geom2d_path_get_steps(const G2DPath *path, double ds_min, double *out_steps);
void geom2d_path_get_points_at_steps(const G2DPath *path, const double *steps, int len_points, G2DPoint *out_points);
int geom2d_path_get_len_corners(const G2DPath *path);
void geom2d_path_get_corner_steps(const G2DPath *path, double *out_steps);


#endif /* CGEOM_PATH_H */
