/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */

#ifndef _DETRESPONSE_SKYGRID_H
#define _DETRESPONSE_SKYGRID_H
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
#include <lal/SeqFactories.h>
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
#include "cmdline.h"

/* number of grid points in declination and right ascension */
#define NUM_DEC 11
#define NUM_RA  25

extern struct gengetopt_args_info args_info;

typedef REAL4 skygrid_t[NUM_RA * NUM_DEC];

REAL4 skygrid_avg(const skygrid_t response);
void  skygrid_square(skygrid_t square, const skygrid_t input);
REAL4 skygrid_rms(const skygrid_t input);
void  skygrid_sqrt(skygrid_t result, const skygrid_t input);
INT4  skygrid_copy(skygrid_t dest, const skygrid_t src);
void  skygrid_print(const LIGOTimeGPS * gps, const skygrid_t input,
                    const char * filename);
void  skygrid_fabs(skygrid_t absgrid, const skygrid_t input);
void  skygrid_add(skygrid_t sum, const skygrid_t a, const skygrid_t b);
void  skygrid_subtract(skygrid_t sum, const skygrid_t a, const skygrid_t b);
void  skygrid_scalar_mult(skygrid_t result, const skygrid_t a, REAL4 b);
void  skygrid_zero(skygrid_t a);

void multiply_vectors(REAL4Vector * out,
                      const REAL4Vector * a, const REAL4Vector * b);

#endif
