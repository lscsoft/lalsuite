/*
 * misc.h - SGWB Standalone Analysis Pipeline
 *        - Miscellaneous Function Prototypes
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

#ifndef MISC_H
#define MISC_H

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

/* cut a time series between given start and end times */
REAL4TimeSeries *cut_time_series(REAL4TimeSeries *input,
    LIGOTimeGPS start,
    UINT4 duration);

#endif /* MISC_H */

/*
 * vim: et
 */
