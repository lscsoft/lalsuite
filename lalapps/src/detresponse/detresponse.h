/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */

#ifndef _DETRESPONSE_DETRESPONSE_H
#define _DETRESPONSE_DETRESPONSE_H
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


/* macro for minimum of two arguments */
#define DETRESPONSE_MIN(a, b) (((a) < (b)) ? (a) : (b))

/* number of grid points in declination and right ascension */
#define NUM_DEC 11
#define NUM_RA  25

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

/* wrap fopen(3) and fclose(3) */
FILE *xfopen(const char *path, const char *mode);
int   xfclose(FILE *stream);

/* wrap strncpy(3) */
char *mystrlcpy(char *dst, const char *src, size_t len);

/* print detector parameters */
void PrintLALDetector(LALDetector * const detector);

/* print source params */
void print_source_maybe(const LALSource * source);

REAL8 deg_to_rad(REAL8 degrees);

void  set_source_params(LALSource * source, const char * name, REAL8 ra_rad,
                        REAL8 dec_rad, REAL8 orien_rad);

void generate_timeseries_response(LALStatus * status);
void compute_skygrid(LALStatus * status);

void multiply_vectors(REAL4Vector * out,
                      const REAL4Vector * a, const REAL4Vector * b);

/* globals */
extern int lalDebugLevel;
extern int verbosity_level;
extern struct gengetopt_args_info args_info;

#endif
