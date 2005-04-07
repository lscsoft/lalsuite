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
#include "skygrid.h"
#include "cmdline.h"

#include "set_source_params.h"


REAL8 deg_to_rad(REAL8 degrees);

void generate_timeseries_response(LALStatus * status);
void compute_skygrid(LALStatus * status, EphemerisData *p_ephemeris_data,
                     skygrid_t grid_cros_sq, skygrid_t grid_plus_sq,
                     skygrid_t grid_sum_sq, skygrid_t grid_relfreq);

/* globals */
extern int lalDebugLevel;
extern int verbosity_level;
extern struct gengetopt_args_info args_info;

#endif
