#include "base.h"
#include "path.h"

#ifndef CGEOM_BEAM_APERTURE_H
#define CGEOM_BEAM_APERTURE_H

typedef struct twiss_data
{
    double x;     // closed orbit x
    double y;     // closed orbit y
    double betx;  // beta x
    double bety;  // beta y
    double dx;    // dispersion x
    double dy;    // dispersion y
    double delta; // relative energy deviation
    double gamma; // relativistic gamma
} G2DTwissData;

typedef struct beam_data
{
    double emitx_norm;        // normalized emittance x
    double emity_norm;        // normalized emittance y
    double delta_rms;         // rms energy spread
    double tol_co;            // tolerance for closed orbit
    double tol_disp;          // tolerance for normalized dispersion
    double tol_disp_ref_dx;   // tolerance for reference dispersion derivative
    double tol_disp_ref_beta; // tolerance for reference dispersion beta
    double tol_energy;        // tolerance for energy error
    double tol_betabeating;   // tolerance for betabeating in sigma
    double halo_x ;           // n sigma of horizontal halo
    double halo_y ;           // n sigma of vertical halo
    double halo_r ;           // n sigma of 45 degree halo
    double halo_primary;      // n sigma of primary halo
} G2DBeamData;

typedef struct aperture_data
{
    G2DPoint *points; // points defining the aperture shape
    int n_points;     // number of points defining the aperture shape
    double tol_r;     // radial tolerance for point-in-aperture check
    double tol_x;     // horizontal tolerance for point-in-aperture check
    double tol_y;     // vertical tolerance for point-in-aperture check
} G2DBeamApertureData;

void generate_beam_envelope_from_sigma_xy(G2DBeamData *beam_data, G2DTwissData *twiss_data, G2DBeamApertureData *aperture_data, G2DSegment *segments, int n_segments);

#endif // CGEOM_BEAM_APERTURE_H
