/************************************ <lalVerbatim file="PrintFTSeriesHV">
Author: J. T. Whelan <whelan@oates.utb.edu>
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{PrintFTSeries.h}}
\label{s:PrintFTSeries.h}

This is a simple utility to print time and frequency series into a
file.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/PrintFTSeries.h>
\end{verbatim}

\noindent Provides prototype information for the routines in 
\verb+PrintTimeSeries.c+ and \verb+PrintFrequencySeries.c+.

\vfill{\footnotesize\input{PrintFTSeriesHV}}
\newpage\input{PrintTimeSeriesC}
\newpage\input{PrintFTSeriesTestC}
\newpage\input{PrintFrequencySeriesC}

</lalLaTeX> */

#ifndef _PRINTFTSERIES_H
#define _PRINTFTSERIES_H

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( PRINTFTSERIESH, "$Id$" );

void LALI2PrintTimeSeries( INT2TimeSeries *series , const CHAR *filename );
void LALI4PrintTimeSeries( INT4TimeSeries *series , const CHAR *filename );
void LALI8PrintTimeSeries( INT8TimeSeries *series , const CHAR *filename );
void LALU2PrintTimeSeries( UINT2TimeSeries *series , const CHAR *filename );
void LALU4PrintTimeSeries( UINT4TimeSeries *series , const CHAR *filename );
void LALU8PrintTimeSeries( UINT8TimeSeries *series , const CHAR *filename );
void LALPrintTimeSeries( REAL4TimeSeries *series , const CHAR *filename );
void LALSPrintTimeSeries( REAL4TimeSeries *series , const CHAR *filename );
void LALDPrintTimeSeries( REAL8TimeSeries *series , const CHAR *filename );
void LALCPrintTimeSeries( COMPLEX8TimeSeries *series , const CHAR *filename );
void LALZPrintTimeSeries( COMPLEX16TimeSeries *series , const CHAR *filename );

void LALI2PrintFrequencySeries( INT2FrequencySeries *series , const CHAR *filename );
void LALI4PrintFrequencySeries( INT4FrequencySeries *series , const CHAR *filename );
void LALI8PrintFrequencySeries( INT8FrequencySeries *series , const CHAR *filename );
void LALU2PrintFrequencySeries( UINT2FrequencySeries *series , const CHAR *filename );
void LALU4PrintFrequencySeries( UINT4FrequencySeries *series , const CHAR *filename );
void LALU8PrintFrequencySeries( UINT8FrequencySeries *series , const CHAR *filename );
void LALPrintFrequencySeries( REAL4FrequencySeries *series , const CHAR *filename );
void LALSPrintFrequencySeries( REAL4FrequencySeries *series , const CHAR *filename );
void LALDPrintFrequencySeries( REAL8FrequencySeries *series , const CHAR *filename );
void LALCPrintFrequencySeries( COMPLEX8FrequencySeries *series , const CHAR *filename );
void LALZPrintFrequencySeries( COMPLEX16FrequencySeries *series , const CHAR *filename );

#ifdef  __cplusplus
}
#endif

#endif /* _PRINTFTSERIES_H */
