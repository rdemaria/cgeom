#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "base.h"
#include "polygon_dist.h"


static double point_segment_dist_sq(G2DPoint p, G2DPoint a, G2DPoint b)
{
    G2DPoint ab = geom2d_sub(b, a);
    G2DPoint ap = geom2d_sub(p, a);
    double denominator = geom2d_dot(ab, ab);
    if (denominator == 0.0) {
        G2DPoint d = geom2d_sub(p, a);
        return geom2d_dot(d, d);
    }
    double t = geom2d_dot(ap, ab) / denominator;
    t = geom2d_clamp(t, 0.0, 1.0);
    G2DPoint proj = (G2DPoint){ a.x + ab.x * t, a.y + ab.y * t };
    G2DPoint d = geom2d_sub(p, proj);
    return geom2d_dot(d, d);
}


static double orient(G2DPoint a, G2DPoint b, G2DPoint c)
{
    return geom2d_cross(geom2d_sub(b, a), geom2d_sub(c, a));
}


static int on_segment(G2DPoint a, G2DPoint b, G2DPoint p)
{
    return (fmin(a.x, b.x) <= p.x && p.x <= fmax(a.x, b.x) &&
        fmin(a.y, b.y) <= p.y && p.y <= fmax(a.y, b.y));
}


static int segments_intersect(G2DPoint a, G2DPoint b, G2DPoint c, G2DPoint d)
{
    double o1 = orient(a, b, c);
    double o2 = orient(a, b, d);
    double o3 = orient(c, d, a);
    double o4 = orient(c, d, b);

    if ((o1 > 0) != (o2 > 0) && (o3 > 0) != (o4 > 0)) {
        return 1;
    }
    const double eps = 1e-12;
    if (fabs(o1) < eps && on_segment(a, b, c)) return 1;
    if (fabs(o2) < eps && on_segment(a, b, d)) return 1;
    if (fabs(o3) < eps && on_segment(c, d, a)) return 1;
    if (fabs(o4) < eps && on_segment(c, d, b)) return 1;
    return 0;
}


static double segment_segment_dist_sq(G2DPoint a, G2DPoint b, G2DPoint c, G2DPoint d)
{
    if (segments_intersect(a,b,c,d)) return 0.0;

    double d2 = point_segment_dist_sq(a, c, d);
    double t2 = point_segment_dist_sq(b, c, d);
    if (t2 < d2) d2 = t2;

    t2 = point_segment_dist_sq(c, a, b);
    if (t2 < d2) d2 = t2;

    t2 = point_segment_dist_sq(d, a, b);
    if (t2 < d2) d2 = t2;

    return d2;
}


double geom2d_poly_boundary_distance_bruteforce(
    const G2DPoint* a,
    size_t len_a,
    const G2DPoint* b,
    size_t len_b
)
/* Compute minimal distance between two polygons.

Contract: len_a=len(a); len_b=len(b)
*/
{
    // polygons are closed => edges = count-1
    if (len_a < 2 || len_b < 2) return 0.0;

    size_t a_edges = len_a - 1;
    size_t b_edges = len_b - 1;

    double best2 = INFINITY;
    for (size_t i = 0; i < a_edges; ++i) {
        G2DPoint a1 = a[i];
        G2DPoint a2 = a[i+1];
        for (size_t j = 0; j < b_edges; ++j) {
            G2DPoint b1 = b[j];
            G2DPoint b2 = b[j+1];
            double d2 = segment_segment_dist_sq(a1, a2, b1, b2);
            if (d2 < best2) {
                best2 = d2;
                if (best2 == 0.0) return 0.0;
            }
        }
    }
    return sqrt(best2);
}


static double polygon_signed_area_closed(const G2DPoint *p, size_t count)
{
    /* count includes repeated last point */
    double area = 0.0;
    for (size_t i = 0; i + 1 < count; ++i) {
        area += p[i].x * p[i + 1].y - p[i + 1].x * p[i].y;
    }
    return 0.5 * area;
}


static void reverse_closed_polygon(G2DPoint *p, size_t count)
{
    size_t n = count - 1; /* unique vertices */
    for (size_t i = 0; i < n / 2; ++i) {
        G2DPoint tmp = p[i];
        p[i] = p[n - 1 - i];
        p[n - 1 - i] = tmp;
    }
    p[count - 1] = p[0]; /* re-close */
}


static size_t advance_support_index(
    const G2DPoint *p,
    size_t count,
    size_t idx,
    G2DPoint direction
) {
    size_t n = count - 1; /* unique vertices */
    double best = p[idx].x * direction.x + p[idx].y * direction.y;

    for (;;) {
        size_t next = idx + 1;
        if (next == n) next = 0;

        double value = p[next].x * direction.x + p[next].y * direction.y;
        if (value > best + 1e-12) {
            idx = next;
            best = value;
        } else {
            break;
        }
    }
    return idx;
}


static size_t support_index_linear(
    const G2DPoint *p,
    size_t count,
    G2DPoint direction
) {
    size_t n = count - 1; /* unique vertices */
    size_t best_i = 0;
    double best = p[0].x * direction.x + p[0].y * direction.y;

    for (size_t i = 1; i < n; ++i) {
        double v = p[i].x * direction.x + p[i].y * direction.y;
        if (v > best) {
            best = v;
            best_i = i;
        }
    }
    return best_i;
}


double geom2d_polygon_boundary_distance_convex_nested(
    G2DPoint *a, size_t len_a,
    G2DPoint *b, size_t len_b
)
/* Compute minimal distance between two polygons (`a` must be convex).

Contract: len_a=len(a); len_b=len(b)
*/
{
    if (len_a < 4 || len_b < 4) return 0.0;

    /* enforce CCW */
    if (polygon_signed_area_closed(a, len_a) < 0.0)
        reverse_closed_polygon(a, len_a);
    if (polygon_signed_area_closed(b, len_b) < 0.0)
        reverse_closed_polygon(b, len_b);

    size_t a_edges = len_a - 1;
    size_t b_edges = len_b - 1;

    /* start indices: minimum-angle outward normal on each polygon */
    size_t a_start = 0, b_start = 0;
    double best_angle = INFINITY;

    for (size_t i = 0; i < a_edges; ++i) {
        G2DPoint e = { a[i + 1].x - a[i].x, a[i + 1].y - a[i].y };
        G2DPoint n = { e.y, -e.x };
        double ang = atan2(n.y, n.x);
        if (ang < 0) ang += 2.0 * M_PI;
        if (ang < best_angle) { best_angle = ang; a_start = i; }
    }

    best_angle = INFINITY;
    for (size_t i = 0; i < b_edges; ++i) {
        G2DPoint e = { b[i + 1].x - b[i].x, b[i + 1].y - b[i].y };
        G2DPoint n = { e.y, -e.x };
        double ang = atan2(n.y, n.x);
        if (ang < 0) ang += 2.0 * M_PI;
        if (ang < best_angle) { best_angle = ang; b_start = i; }
    }

    size_t i = 0, j = 0;

    G2DPoint ea0 = { a[a_start + 1].x - a[a_start].x, a[a_start + 1].y - a[a_start].y };
    G2DPoint eb0 = { b[b_start + 1].x - b[b_start].x, b[b_start + 1].y - b[b_start].y };
    G2DPoint na0 = { ea0.y, -ea0.x };
    G2DPoint nb0 = { eb0.y, -eb0.x };

    /* pick the earlier direction in CCW order (same rule as merge) */
    double cr0 = na0.x * nb0.y - na0.y * nb0.x;
    G2DPoint d0 = (cr0 > 1e-15) ? na0 : nb0;

    /* true argmax init */
    size_t ia = support_index_linear(a, len_a, d0);
    size_t ib = support_index_linear(b, len_b, d0);

    double best = INFINITY;

    for (size_t steps = 0; steps < a_edges + b_edges; ++steps) {
        size_t ai = (a_start + i) % a_edges;
        size_t bj = (b_start + j) % b_edges;

        G2DPoint ea = { a[ai + 1].x - a[ai].x, a[ai + 1].y - a[ai].y };
        G2DPoint eb = { b[bj + 1].x - b[bj].x, b[bj + 1].y - b[bj].y };

        G2DPoint na = { ea.y, -ea.x };
        G2DPoint nb = { eb.y, -eb.x };

        double cr = na.x * nb.y - na.y * nb.x;
        G2DPoint direction;

        if (cr > 1e-15) {
            direction = na;
            i = (i + 1) % a_edges;
        } else if (cr < -1e-15) {
            direction = nb;
            j = (j + 1) % b_edges;
        } else {
            direction = na;
            i = (i + 1) % a_edges;
            j = (j + 1) % b_edges;
        }

        ia = advance_support_index(a, len_a, ia, direction);
        ib = advance_support_index(b, len_b, ib, direction);

        double gap =
            (a[ia].x * direction.x + a[ia].y * direction.y) -
            (b[ib].x * direction.x + b[ib].y * direction.y);

        double len = hypot(direction.x, direction.y);
        if (len > 0.0) {
            double dist = gap / len;
            if (dist < best) best = dist;
        }
    }

    if (best < 0.0 && best > -1e-9) best = 0.0;
    return best;
}
