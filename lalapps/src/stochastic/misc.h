/*
 * misc.h - SGWB Standalone Analysis Pipeline
 *        - Miscellaneous Function Prototypes
 * 
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 *
 * $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/LALStdio.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <lalapps.h>

/* can't use round() as its not C89 */
double myround(double x);

/* cut a time series between given start and end times */
REAL4TimeSeries *cut_time_series(LALStatus *status,
    REAL4TimeSeries *input,
    LIGOTimeGPS start,
    UINT4 duration);

/*
 * vim: et
 */
