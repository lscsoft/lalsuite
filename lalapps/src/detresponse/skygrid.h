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

#define DR_ANGULAR_OFFSET 0.005
#define DR_TWOPI ((1. - DR_ANGULAR_OFFSET) * LAL_TWOPI)

extern struct gengetopt_args_info args_info;

typedef REAL8Array * skygrid_t;

void init_ephemeris(LALStatus *status, EphemerisData *ephemeris_data);
void cleanup_ephemeris(LALStatus *status, EphemerisData *ephemeris_data);
void init_skygrid(LALStatus *status);
skygrid_t * alloc_skygrid(LALStatus *status, skygrid_t *g);
void free_skygrid(LALStatus *status, skygrid_t *skygrid);
void cleanup_skygrid(LALStatus *status);
REAL8 skygrid_avg(LALStatus *status, const skygrid_t response);
void  skygrid_square(LALStatus *status, skygrid_t square, const skygrid_t input);
REAL8 skygrid_rms(LALStatus *status, const skygrid_t input);
void  skygrid_sqrt(LALStatus *status, skygrid_t result, const skygrid_t input);
INT4  skygrid_copy(LALStatus *status, skygrid_t dest, const skygrid_t src);
void  skygrid_print(LALStatus *status, const LIGOTimeGPS * gps, const skygrid_t input,
                    const char * filename);
void  skygrid_fabs(LALStatus *status, skygrid_t absgrid, const skygrid_t input);
void  skygrid_add(LALStatus *status, skygrid_t sum, const skygrid_t a, const skygrid_t b);
void  skygrid_subtract(LALStatus *status, skygrid_t sum, const skygrid_t a, const skygrid_t b);
void  skygrid_scalar_mult(LALStatus *status, skygrid_t result, const skygrid_t a, REAL8 b);
void  skygrid_zero(LALStatus *status, skygrid_t a);

void multiply_vectors(LALStatus *status, REAL4Vector * out,
                      const REAL4Vector * a, const REAL4Vector * b);

#endif
