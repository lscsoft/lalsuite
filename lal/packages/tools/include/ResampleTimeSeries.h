/**** <lalVerbatim file="ResampleTimeSeriesHV">
 * Author: Brown, D. A.
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{ResampleTimeSeries.h}}
 *
 * High-level routines for resampleing time series data.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/ResampleTimeSeries.h>
 * \end{verbatim}
 *
 * Provides a high level interface for resampling time series data.
 *
 **** </lalLaTeX> */

#include <lal/LALDatatypes.h>

#ifndef _RESAMPLETIMESERIES_H
#define _RESAMPLETIMESERIES_H

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( RESAMPLETIMESERIESH, "$Id$" );

/**** <lalLaTeX>
 *
 * \subsection*{Error conditions}
 *
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define RESAMPLETIMESERIESH_ENULL 1
#define RESAMPLETIMESERIESH_ENNUL 2
#define RESAMPLETIMESERIESH_EZERO 3

#define RESAMPLETIMESERIESH_MSGENULL "No calibration generated: Null pointer"
#define RESAMPLETIMESERIESH_MSGENNUL "No calibration generated: Non-null pointer"
#define RESAMPLETIMESERIESH_MSGEZERO "Length of input time series is zero"
/**** </lalErrTable> */

/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 *
 **** </lalLaTeX> */

typedef struct
tagResampleTimeSeriesParams
{
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  UINT4         length;
}
ResampleTimeSeriesParams;

/**** <lalLaTeX>
 * 
 * \vfill{\footnotesize\input{ResampleTimeSeriesHV}}
 * \newpage\input{ResampleTimeSeriesC}
 *
 **** </lalLaTeX> */

void
LALResampleREAL4TimeSeries(
    LALStatus                  *status,
    REAL4TimeSeries            *data, 
    ResampleTimeSeriesParams   *params
    );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _RESAMPLETIMESERIES_H */

