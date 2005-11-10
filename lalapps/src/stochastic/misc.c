/*
 * misc.c - SGWB Standalone Analysis Pipeline
 *        - Miscellaneous Functions
 * 
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 *
 * $Id$
 */

#include "misc.h"

NRCSID(DATAINPUTC, "$Id$");
RCSID("$Id$");

/* can't use round() as its not C89 */
double myround(double x)
{
  return(x < 0 ? ceil(x - 0.5): floor(x + 0.5));
}

/* cut a time series between given start and end times */
REAL4TimeSeries *cut_time_series(LALStatus *status,
    REAL4TimeSeries *input,
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
  LAL_CALL(LALCreateREAL4TimeSeries(status, &series, input->name, start, \
        input->f0, input->deltaT, input->sampleUnits, length), status);

  /* cut time series */
  LAL_CALL(LALCutREAL4TimeSeries(status, &series, input, first, length), \
      status);

  return(series);
}

/*
 * vim: et
 */
