/*************************** <lalVerbatim file="BandPassTimeSeriesHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{BandPassTimeSeries.h}}
\label{s:BandPassTimeSeries.h}

Provides routines to low- or high-pass filter a time series.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/BandPassTimeSeries.h>
\end{verbatim}

\noindent This header covers routines that apply a time-domain low- or
high-pass filter to a data series of type \verb@<datatype>TimeSeries@.
Further documentation is given in the individual routines' modules.

******************************************************* </lalLaTeX> */

#ifndef _BANDPASSTIMESERIES_H
#define _BANDPASSTIMESERIES_H

#include <lal/LALStdlib.h>
#include <lal/IIRFilter.h>
#include <lal/ZPGFilter.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID(BANDPASSTIMESERIESH,"$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define BANDPASSTIMESERIESH_ENUL 1
#define BANDPASSTIMESERIESH_EBAD 2

#define BANDPASSTIMESERIESH_MSGENUL "Unexpected null pointer in arguments"
#define BANDPASSTIMESERIESH_MSGEBAD "Bad filter parameters"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{PassBandParamStruc}}
\idx[Type]{PassBandParamStruc}

This structure stores data used for constructing a low- or high-pass
filter: either the order and characteristic frequency of the filter,
or the frequencies and desired attenuations at the ends of some
transition band.  In the latter case, a nonzero filter order parameter
\verb@n@ indicates a maximum allowed order.  The fields are:

\begin{description}
\item[\texttt{CHAR *name}] A user-assigned name.

\item[\texttt{INT4 n}] The maximum desired filter order (actual order
  may be less if specified attenuations do not require a high order).

\item[\texttt{REAL8 f1}, \texttt{f2}] The reference frequencies of the
  transition band.

\item[\texttt{REAL8 a1}, \texttt{a2}] The minimal desired attenuation
  factors at the reference frequencies.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagPassBandParamStruc{
  CHAR *name;
  INT4 nMax;
  REAL8 f1;
  REAL8 f2;
  REAL8 a1;
  REAL8 a2;
} PassBandParamStruc;

/* <lalLaTeX>
\vfill{\footnotesize\input{BandPassTimeSeriesHV}}
</lalLaTeX> */

/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{ButterworthTimeSeriesC}
</lalLaTeX> */

int XLALButterworthREAL4TimeSeries( REAL4TimeSeries *series, PassBandParamStruc *params );
int XLALButterworthREAL8TimeSeries( REAL8TimeSeries *series, PassBandParamStruc *params );
int XLALLowPassREAL4TimeSeries( REAL4TimeSeries *series,
    REAL8 frequency, REAL8 amplitude, INT4 filtorder );
int XLALLowPassREAL8TimeSeries( REAL8TimeSeries *series,
    REAL8 frequency, REAL8 amplitude, INT4 filtorder );
int XLALHighPassREAL4TimeSeries( REAL4TimeSeries *series,
    REAL8 frequency, REAL8 amplitude, INT4 filtorder );
int XLALHighPassREAL8TimeSeries( REAL8TimeSeries *series,
    REAL8 frequency, REAL8 amplitude, INT4 filtorder );



void
LALButterworthREAL4TimeSeries( LALStatus          *status,
			       REAL4TimeSeries    *series,
			       PassBandParamStruc *params );

void
LALButterworthREAL8TimeSeries( LALStatus          *status,
			       REAL8TimeSeries    *series,
			       PassBandParamStruc *params );

void
LALDButterworthREAL4TimeSeries( LALStatus          *status,
				REAL4TimeSeries    *series,
				PassBandParamStruc *params );

/* Chebyshev filters should also be added, but I'm too busy to write
   the routines now. */

/* Test program. */

/* <lalLaTeX>
\newpage\input{BandPassTestC}
</lalLaTeX> */

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _BANDPASSTIMESERIES_H */
