/**** <lalVerbatim file="ResampleTimeSeriesHV">
 * Author: Brown, D. A.
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{ResampleTimeSeries.h}}
 *
 * Provides routines to resample a time series. At present only integer
 * downsampling or \verb|REAL4TimeSeries| by a power of two is supported.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/ResampleTimeSeries.h>
 * \end{verbatim}
 *
 * \noindent This header covers routines that resample time series by applying
 * a low pass filter and decimating the resulting time series. Further
 * documentation is given in the individual routines' modules.
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

/************************************************************* <lalLaTeX>
\subsection*{Types}
\idx[Type]{ResampleTSFilter}
\idx[Type]{ResampleTSFilterParams}
\idx[Type]{ResampleTSParams}

\subsubsection*{Enum \texttt{ResampleTSFilter}}
This enum type contains the different low pass filters available to
prevent power above the new Nyquist frequency entering the resampled 
time series due to aliasing. The options are:

\begin{description}
\item[\texttt{defaultButterworth}] An IIR butterwoth filter of order 20 
with attenuation 0.1 at the new Nyquist frequency. See the package tdfilters
for documentation of butterworth filters in LAL.

\item[\texttt{LDASfirLP}] For downsampling by a factor of 2, 4 or 8 an
implementation of the FIR filter used by the LDAS datacondAPI resample().
This is provided for testing the result of standalone codes and codes 
running under LDAS. The LDAS filter provided here has filter order parameter
10, so the order of the filter is $2 \times 10 \times q$ where $q$ is the
resampling ratio.
\end{description}

\subsubsection*{Union \texttt{ResampleTSFilterParams}}
This union is provided so that the code can store the parameters of the 
filter in a place accessible by the user for user designed low pass filters.
This is not presently implemented and this structure may be ignored.
\begin{description}
\item[\texttt{butterworth}] A structure of type \verb|PassBandParamStruc| used
to store the parameters of the butterworth filter used to perform low pass
filtering.

\item[\texttt{iirfilter}] A structure of type \verb|REAL8IIRFilter| used
to store the parameters of the IIR or FIR filter used to perform low pass
filtering.
\end{description}

\subsubsection*{Structure \texttt{tagResampleTimeSeriesParams}}
This structure controls the behaviour of the resampling function.
\begin{description}
\item[\texttt{deltaT}] The sample interval desired in the down sampled time
series.

\item[\texttt{filterType}] The type of filter with which to perform the low
pass filtering.

\item[\texttt{filterParams}] Filter parameters for the low pass filter. 
Presently ignored.
\end{description}

************************************************************* </lalLaTeX> */

typedef enum
{
  defaultButterworth,
  LDASfirLP
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

int XLALResampleREAL4TimeSeries( REAL4TimeSeries *series, REAL8 dt );
int XLALResampleREAL8TimeSeries( REAL8TimeSeries *series, REAL8 dt );

void
LALResampleREAL4TimeSeries(
    LALStatus          *status,
    REAL4TimeSeries    *ts, 
    ResampleTSParams   *params
    );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _RESAMPLETIMESERIES_H */

