/*----------------------------------------------------------------------- 
 * 
 * File Name: BandPassTimeSeries.h
 * 
 * Author: Creighton, T. D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------*/

/* <lalLaTeX>

\section{Header \texttt{BandPassTimeSeries.h}}

Provides routines to low- or high-pass filter a time series.

\subsection{Synopsis}
\begin{verbatim}
#include "BandPassTimeSeries.h"
\end{verbatim}

\noindent This header covers routines that apply a time-domain low- or
high-pass filter to a data series of type \verb@<datatype>TimeSeries@.
Further documentation is given in the individual routines' modules.

</lalLaTeX> */

#ifndef _BANDPASSTIMESERIES_H
#define _BANDPASSTIMESERIES_H

#include "LALStdlib.h"
#include "IIRFilter.h"
#include "ZPGFilter.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID(BANDPASSTIMESERIESH,"$Id$");

/* <lalLaTeX>

\subsection{Error conditions}
\begin{tabular}{|c|l|l|}
\hline
status & status                    & Explanation                           \\
 code  & description               &                                       \\
\hline
\tt 1  & \tt Null pointer          & Missing a required pointer.           \\
\tt 2  & \tt Bad filter parameters & Filter creation parameters outside of \\
       &                           & acceptable ranges.                    \\
\hline
\end{tabular}

</lalLaTeX> */

#define BANDPASSTIMESERIES_ENUL 1
#define BANDPASSTIMESERIES_EBAD 2

#define BANDPASSTIMESERIES_MSGENUL "Null pointer"
#define BANDPASSTIMESERIES_MSGEBAD "Bad filter parameters"

/* <lalLaTeX>

\subsection{Structures}

\begin{verbatim}
struct PassBandParamStruc
\end{verbatim}

\noindent This structure stores data used for constructing a low- or
high-pass filter: either the order and characteristic frequency of the
filter, or the frequencies and desired attenuations at the ends of
some transition band.  In the latter case, a nonzero filter order
parameter n indicates a maximum allowed order.  The fields are:

\begin{description}
\item[\texttt{CHAR *name}] A user-assigned name.

\item[\texttt{INT4 n}] The maximum desired filter order (actual order
  may be less if specified attenuations do not require a high order).

\item[\texttt{REAL8 f1}, \texttt{f2}] The reference frequencies of the
  transition band.

\item[\texttt{REAL8 a1}, \texttt{a2}] The minimal desired attenuation
  factors at the reference frequencies.
\end{description}

</lalLaTeX> */

typedef struct tagPassBandParamStruc{
  CHAR *name;
  INT4 nMax;
  REAL8 f1;
  REAL8 f2;
  REAL8 a1;
  REAL8 a2;
} PassBandParamStruc;

/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{ButterworthTimeSeriesC}
</lalLaTeX> */
void ButterworthREAL4TimeSeries(Status             *stat,
				REAL4TimeSeries    *series,
				PassBandParamStruc *params);

void ButterworthREAL8TimeSeries(Status             *stat,
				REAL8TimeSeries    *series,
				PassBandParamStruc *params);

/* Chebyshev filters should also be added, but I'm too busy to write
   the routines now. */

/* <lalLaTeX>
\newpage\input{BandPassTestC}
</lalLaTeX> */


#ifdef  __cplusplus
}
#endif

#endif /* _BANDPASSTIMESERIES_H */
