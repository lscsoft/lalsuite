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
#include <lal/Units.h>

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/Units.h>

#include <lal/PrintFTSeries.h>
#include <lal/StreamOutput.h>

#include "util.h"

/* macro for minimum of two arguments */
#define DETRESPONSE_MIN(a, b) (((a) < (b)) ? (a) : (b))

/* wrap fopen(3) and fclose(3) */
FILE *xfopen(const char *path, const char *mode);
int   xfclose(FILE *stream);

/* wrap strncpy(3) */
char *strlcpy(char *dst, const char *src, size_t len);

/* print detector parameters */
void PrintLALDetector(LALDetector * const detector);

/* print source params */
void print_source_maybe(const LALSource * source);

REAL8 deg_to_rad(REAL8 degrees);

void  set_source_params(LALSource * source, const char * name, REAL8 ra_rad,
                        REAL8 dec_rad, REAL8 orien_rad);


/* globals */
extern int lalDebugLevel;
extern int verbose_level;

#endif
