#ifndef CGEOM_BEAM_PIPES_H
#define CGEOM_BEAM_PIPES_H

#include "base.h"
#include "path.h"

typedef struct profile
{
    G2DPath path; // Cross-sectional shape of the beam pipe
    G3DTransform transform; // Position and orientation in 3D space of the profile
} G2DBeamPipeProfile;

typedef struct
{
    G2DBeamPipeProfile *profiles; // Array of beam pipe profiles along its length
    int n_profiles;               // Number of profiles
    G2DPath centerline;      // Centerline path of the beam pipe in 2D
} G2DBeamPipe;


/* Potential extensions to 3D
typedef struct profile
{
    G3DPath path; // Cross-sectional shape of the beam pipe
    G3DTransform transform; // Position and orientation in 3D space of the profile
} G3DBeamPipeProfile;

typedef struct
{
    G3DBeamPipeProfile *profiles; // Array of beam pipe profiles along its length
    int n_profiles;               // Number of profiles
    G3DPath centerline;      // Centerline path of the beam pipe in 2D
} G3DBeamPipe;
 */



#endif // CGEOM_BEAM_PIPES_H


