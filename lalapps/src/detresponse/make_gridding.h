/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
#ifndef _DETRESPONSE_MAKE_GRIDDING_H
#define _DETRESPONSE_MAKE_GRIDDING_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>

#include <lal/AVFactories.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/Velocity.h>
#include <lal/Units.h>
#include <lal/VectorOps.h>
#include <lal/SkyCoordinates.h>

#include <lal/PrintFTSeries.h>
#include <lal/StreamOutput.h>

#include "util.h"

/*
 * Per axis/coord. grid geometry
 *
 * Regular grid:   regular division of the coordinate; total
 *                 grid cells == num_ra or num_dec
 * Irregular grid: irregular division of the coordinate;
 *                 total grid points == num_ra or num_dec
 * Variable grid:  variable grid; number of grid points varies
 *                 as a function of the other coordinate;
 *                 we only support RA grid vary as cos(dec)
 *
 * NB: only the RA grid supports variable gridding
 */
typedef enum 
{
  DETRESP_REGGRID, 
  DETRESP_IRRGRID, 
  DETRESP_VARGRID,
} 
gridding_geom_t;

/*
 * Struct that will contain the coordinates of every 
 * grid cell center
 */
typedef 
struct gridding_tag 
{
  gridding_geom_t ra_geom;
  gridding_geom_t dec_geom;
  REAL8Vector *ra;
  REAL8Vector **ra_irr;  /* irregular grid, Dec dependent */
  REAL8Vector *dec;   
}
gridding_t;


void init_gridding(gridding_t *gridding);

void make_gridding(LALStatus *status, gridding_t *gridding, 
                   UINT4 num_ra, gridding_geom_t ra_geom, 
                   UINT4 num_dec, gridding_geom_t dec_geom,
                   LIGOTimeGPS *gps);

void cleanup_gridding(LALStatus *status, gridding_t *gridding);

void zero_gridding(LALStatus *status, gridding_t *gridding);

void print_gridding(gridding_t *gridding, char *filename);

#endif
