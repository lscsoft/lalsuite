/*
 * data_input.h - SGWB Standalone Analysis Pipeline
 *              - Data Input Function Prototypes
 *
 * Copyright (C) 2002-2006,2010 Adam Mercer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 *
 */

#ifndef DATA_INPUT_H
#define DATA_INPUT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/LALFrStream.h>
#include <lal/LALStdio.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <lalapps.h>

/* read a LIGO time series */
REAL4TimeSeries *get_ligo_data(LALFrStream *stream,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end);

/* read and high pass filter a GEO time series */
REAL4TimeSeries *get_geo_data(LALFrStream *stream,
    CHAR *channel,
    LIGOTimeGPS start,
    LIGOTimeGPS end,
    INT4 hpf_order,
    REAL8 hpf_frequency,
    REAL8 hpf_attenuation);

/* read a time series */
REAL4TimeSeries *get_time_series(CHAR *ifo,
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

#endif /* DATA_INPUT_H */

/*
 * vim: et
 */
