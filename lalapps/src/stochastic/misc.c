/*
 * misc.c - SGWB Standalone Analysis Pipeline
 *        - Miscellaneous Functions
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

#include "misc.h"


/* cut a time series between given start and end times */
REAL4TimeSeries *cut_time_series(REAL4TimeSeries *input,
    LIGOTimeGPS start,
    UINT4 duration)
{
  /* variables */
  REAL4TimeSeries *series;
  INT4 length;
  INT4 first;

  /* calculate length of segment to cut */
  length = floor((duration / input->deltaT) + 0.5);

  /* get first bin */
  first = (INT4)((start.gpsSeconds - input->epoch.gpsSeconds) / input->deltaT);

  /* allocate memory */
  series = XLALCreateREAL4TimeSeries(input->name, &start, input->f0, \
      input->deltaT, &input->sampleUnits, length);

  /* cut time series */
  series = XLALCutREAL4TimeSeries(input, first, length);

  return(series);
}

/*
 * vim: et
 */
