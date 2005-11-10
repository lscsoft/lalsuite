/*
 * data_input.h - SGWB Standalone Analysis Pipeline
 *              - Data Input Function Prototypes
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

#include <FrameL.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/FrameStream.h>
#include <lal/LALStdio.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <lalapps.h>

/* read a LIGO time series */
REAL4TimeSeries *get_ligo_data(LALStatus *status,
    FrStream *stream,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end);

/* read and high pass filter a GEO time series */
REAL4TimeSeries *get_geo_data(LALStatus *status,
    FrStream *stream,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end,
    INT4 order,
    REAL8 frequency,
    REAL8 attenuation);

/* read a time series */
REAL4TimeSeries *get_time_series(LALStatus *status,
    CHAR *ifo,
    CHAR *cache_file,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end,
    INT4 resample_rate,
    INT4 hpf_order,
    REAL8 hpf_frequency,
    REAL8 hpf_attenuation,
    INT4 geo_hpf_order,
    REAL8 geo_hpf_frequency,
    REAL8 geo_hpf_attenuation,
    INT4 buffer);

/*
 * vim: et
 */
