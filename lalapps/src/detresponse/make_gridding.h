/*
*  Copyright (C) 2007 David Chin
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
 * Mollweide: Mollweide projection 
 *            (see http://mathworld.wolfram.com/MollweideProjection.html
 *
 * NB: only the RA grid supports variable gridding
 */
typedef enum 
{
  DETRESP_REGGRID, 
  DETRESP_IRRGRID, 
  DETRESP_VARGRID,
  DETRESP_MOLLWEIDEGRID,
} 
gridding_geom_t;

typedef enum
{
  DETRESP_HUMANREAD,
  DETRESP_XYPAIRS_ASCII,
  DETRESP_XYPAIRS_BIN,
}
gridding_printmode_t;

/*
 * Struct that will contain the coordinates of every 
 * grid cell center
 */
typedef 
struct tag_gridding_t
{
  LIGOTimeGPS        gps;
  gridding_geom_t ra_geom;
  gridding_geom_t dec_geom;
  REAL8Vector *ra;
  REAL8Vector **ra_irr;  /* irregular grid, Dec dependent */
  REAL8Vector *dec;   
}
gridding_t;


void init_gridding(gridding_t *p_gridding);

void make_gridding(LALStatus *status, gridding_t *p_gridding, 
                   UINT4 num_ra, gridding_geom_t ra_geom, 
                   UINT4 num_dec, gridding_geom_t dec_geom,
                   EphemerisData *p_ephem, LIGOTimeGPS *p_gps);

void cleanup_gridding(LALStatus *status, gridding_t *p_gridding);

void zero_gridding(LALStatus *status, gridding_t *p_gridding);

void print_gridding(gridding_t *p_gridding, char *filename,
                    gridding_printmode_t mode);

void print_ra_grid(gridding_t *p_gridding, char *filename);

void print_dec_grid(gridding_t *p_gridding, char *filename);

#endif
