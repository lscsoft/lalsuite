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
#include <lal/BandPassTimeSeries.h>

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
#define RESAMPLETIMESERIESH_ERATE 4
#define RESAMPLETIMESERIESH_EUPSM 5
#define RESAMPLETIMESERIESH_EHIGH 6
#define RESAMPLETIMESERIESH_ELOG2 7
#define RESAMPLETIMESERIESH_EFILT 8
#define RESAMPLETIMESERIESH_EINVD 9
#define RESAMPLETIMESERIESH_ELDAS 10

#define RESAMPLETIMESERIESH_MSGENULL "Null pointer"
#define RESAMPLETIMESERIESH_MSGENNUL "Non-null pointer"
#define RESAMPLETIMESERIESH_MSGEZERO "Length of input time series is zero"
#define RESAMPLETIMESERIESH_MSGERATE "Sample rate is zero"
#define RESAMPLETIMESERIESH_MSGEUPSM "Cannot upsample"
#define RESAMPLETIMESERIESH_MSGEHIGH "Input sample rate is greater than 32kHz"
#define RESAMPLETIMESERIESH_MSGELOG2 "Only power-of-two resampling is avaliable"
#define RESAMPLETIMESERIESH_MSGEFILT "Unknown filter type"
#define RESAMPLETIMESERIESH_MSGEINVD "Invalid or non-integer resample factor"
#define RESAMPLETIMESERIESH_MSGELDAS "Input resample factor with LDAS FIR"
/**** </lalErrTable> */

/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 *
 **** </lalLaTeX> */

typedef enum
{
  defaultButterworth,
  GDSfirLSOne,
  GDSfirLSTwo,
  GDSfirLSThree,
  GDSfirPMOne,
  LDASorderTen,
  firLSOne,
  firLSTwo,
  firLSThree,
  firEqOne
}
ResampleTSFilter;

typedef union
tagResampleTimeSeriesFilterPars
{
  PassBandParamStruc    butterworth;
  REAL8IIRFilter        iirfilter;
}
ResampleTSFilterParams;

typedef struct
tagResampleTimeSeriesParams
{
  REAL8                   deltaT;
  ResampleTSFilter        filterType;
  ResampleTSFilterParams  filterParams;
}
ResampleTSParams;

/**** <lalLaTeX>
 * 
 * \vfill{\footnotesize\input{ResampleTimeSeriesHV}}
 * \newpage\input{ResampleTimeSeriesC}
 * \newpage\input{ResampleTimeSeriesTestC}
 *
 **** </lalLaTeX> */

void
LALResampleREAL4TimeSeries(
    LALStatus          *status,
    REAL4TimeSeries    *ts, 
    ResampleTSParams   *params
    );

void
LALDecimateREAL4TimeSeries(
    LALStatus          *status,
    REAL4TimeSeries    *ts, 
    ResampleTSParams   *params
    );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _RESAMPLETIMESERIES_H */

