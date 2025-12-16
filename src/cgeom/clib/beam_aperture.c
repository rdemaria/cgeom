#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "base.h"
#include "path.h"
#include "beam_aperture.h"

void geom2d_get_beam_envelope(const G2DBeamData *beam_data, const G2DTwissData *twiss_data, const G2DBeamApertureData *aperture_data, int len_points, G2DPoint *out_points)
/* Create beam envelope based on beam data, twiss data, aperture data is needed to get tolerances on the shape

See pyoptics/aperture.py: get_halo

Contract: len(out_points)=len_points
*/
{
    double x0 = twiss_data->x;  // assuming closed orbit relative to aperture center
    double y0 = twiss_data->y;  // assuming closed orbit relative to aperture center
    double betx = twiss_data->betx;
    double bety = twiss_data->bety;
    double dx = twiss_data->dx;
    double dy = twiss_data->dy;
    double gamma = twiss_data->gamma;

    double emitx_norm = beam_data->emitx_norm;
    double emity_norm = beam_data->emity_norm;
    double delta_rms = beam_data->delta_rms;

    double tol_r = aperture_data->tol_r;
    double tol_x = aperture_data->tol_x;
    double tol_y = aperture_data->tol_y;

    double hr = beam_data->halo_r;
    double hx = beam_data->halo_x;
    double hy = beam_data->halo_y;


    double ex = emitx_norm / gamma;
    double ey = emity_norm / gamma;

    double sigma_x = sqrt(ex * betx + dx * dx * delta_rms * delta_rms);
    double sigma_y = sqrt(ey * bety + dy * dy * delta_rms * delta_rms);

    /*
        We describe the beam of the shape described by hx, hy, and hr as a
        racetrack, leading to the following equations, where sh and sv are the
        width and height of a rectangle, which, when convolved with a circle of
        radius sr (Minkowski sum) yields our racetrack:

        { hx = sh + sr  (width),
        { hy = sv + sr  (height),
        { hr = sqrt(sh ** 2 + sv ** 2) + sr  (radial maximum: through (sh, sv)).

        Solving these for sh, sv, and sr yields the following equations:
    */
    double tmp = sqrt(2) * sqrt((hr - hx) * (hr - hy));
    double sh = hr - hy + tmp;
    double sv = hr - hx + tmp;
    double sr = hx + hy - hr - tmp;

    /*
        Remembering that hx, hy, and hr are specified in sigmas, we convolve the
        beam racetrack with the aperture tolerance racetrack, to get our beam
        envelope racetrack:
    */
    double h = tol_x + sh * sigma_x;
    double v = tol_y + sv * sigma_y;
    double a = tol_r + sr * sigma_x;
    double b = tol_r + sr * sigma_y;

    G2DSegment segments[8];
    G2DPath path;
    path.segments = segments;
    path.len_segments=8;
    geom2d_segments_from_racetrack(h + a, v + b, a, b, path.segments, &path.len_segments);
    geom2d_path_get_n_uniform_points(&path, len_points, out_points);
}
