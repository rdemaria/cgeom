#ifndef CGEOM_BASE_H
#define CGEOM_BASE_H
/* Basic geometric types and functions */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    double x;
    double y;
} G2DPoint;

typedef struct {
    double x;
    double y;
    double z;
} G3DPoint;

typedef struct {
    double mat[4*4];
} G3DTransform;



double geom2d_norm(double x, double y);
double geom2d_dot(G2DPoint a, G2DPoint b);
double geom2d_cross(G2DPoint a, G2DPoint b);
G2DPoint geom2d_sub(G2DPoint a, G2DPoint b);
double geom2d_clamp(double t, double lo, double hi);
double geom2d_points_distance(double x1, double y1, double x2, double y2);
void geom2d_points_translate(double dx, double dy, G2DPoint* points, const int len_points);
double geom2d_elliptic_E(double phi, double k);
double geom2d_elliptic_E_complete(double k);

#endif /* CGEOM_BASE_H */
