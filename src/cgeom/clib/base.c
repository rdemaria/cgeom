#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>

#include "base.h"


double geom2d_norm(double x, double y)
/* Get the norm of a 2D vector */
{
    return sqrt(x * x + y * y);
}

double geom2d_points_distance(double x1, double y1, double x2, double y2)
/* Get the distance between two 2D points */
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    return geom2d_norm(dx, dy);
}

static double geom2d_elliptic_E_adaptive(double a, double b, double eps, int depth, double m)
/* Incomplete elliptic integral of the second kind E(phi, k)
   Uses libm ellint_2/comp_ellint_2 when available (glibc),
   otherwise falls back to a small adaptive Simpson integrator.
*/
{
    double c = 0.5 * (a + b);
    double fa = sqrt(1.0 - m * sin(a) * sin(a));
    double fb = sqrt(1.0 - m * sin(b) * sin(b));
    double fc = sqrt(1.0 - m * sin(c) * sin(c));
    double h = b - a;
    double S = (fa + 4.0 * fc + fb) * h / 6.0;
    double left_c = 0.5 * (a + c);
    double right_c = 0.5 * (c + b);
    double f_left_c = sqrt(1.0 - m * sin(left_c) * sin(left_c));
    double f_right_c = sqrt(1.0 - m * sin(right_c) * sin(right_c));
    double Sleft = (fa + 4.0 * f_left_c + fc) * (h / 2.0) / 6.0;
    double Sright = (fc + 4.0 * f_right_c + fb) * (h / 2.0) / 6.0;
    if (depth <= 0 || fabs(Sleft + Sright - S) < 15.0 * eps)
        return Sleft + Sright + (Sleft + Sright - S) / 15.0;
    return geom2d_elliptic_E_adaptive(a, c, eps / 2.0, depth - 1, m) +
           geom2d_elliptic_E_adaptive(c, b, eps / 2.0, depth - 1, m);
}

static double geom2d_elliptic_E_numeric(double phi, double k)
{
    double m = k * k;
    double sign = (phi >= 0.0) ? 1.0 : -1.0;
    double abs_phi = fabs(phi);
    double period = M_PI; /* integrand is periodic in pi */
    double half_pi = 0.5 * M_PI;
    double base_complete = geom2d_elliptic_E_adaptive(0.0, half_pi, 1e-10, 12, m);
    double per_value = 2.0 * base_complete; /* integral over one period */
    long periods = (long)(abs_phi / period);
    double remainder = abs_phi - periods * period;
    double total = periods * per_value;
    if (remainder > 0.0)
        total += geom2d_elliptic_E_adaptive(0.0, remainder, 1e-10, 12, m);
    return sign * total;
}

double geom2d_elliptic_E(double phi, double k)
{
    if (k < 0.0 || k >= 1.0)
        return NAN;
    return geom2d_elliptic_E_numeric(phi, k);
}

double geom2d_elliptic_E_complete(double k)
{
    if (k < 0.0 || k >= 1.0)
        return NAN;
    return geom2d_elliptic_E_numeric(0.5 * M_PI, k);
}
