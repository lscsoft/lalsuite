/*************** <lalVerbatim file="PrintFTSeriesTestCV"> *******
$Id$
**************** </lalVerbatim> ***********************************/

/* <lalLaTeX>

\subsection{Program \texttt{PrintFTSeriesTest.c}}
\label{s:PrintFTSeriesTest.c}

Tests the routines in \verb@PrintTimeSeries.c@ and
\verb@PrintFrequenceSeries.c@.

\subsubsection*{Usage}
\begin{verbatim}
PrintFTSeriesTest
\end{verbatim}

\subsubsection*{Description}

This program generates and prints a sequence of frequency and time
series; the program itself always returns success, so the testing
function is actually served by examinaton of the output files.

The program as written generates and prints single and double
precision real and complex time and frequency series.  The routines
for integers are not tested (although the formats are the same as in
the \verb+PrintVector+ module).

\subsubsection*{Exit codes}
\begin{tabular}{|c|l|}
\hline
 Code & Explanation                   \\
\hline
\tt 0 & Always returned.              \\
\hline
\end{tabular}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{PrintFTSeriesTestCV}}

</lalLaTeX> */

#ifndef _PRINTVECTOR_H
#include <lal/PrintFTSeries.h>
#ifndef _PRINTVECTOR_H
#define _PRINTVECTOR_H
#endif
#endif

#ifndef _MATH_H
#include <math.h>
#ifndef _MATH_H
#define _MATH_H
#endif
#endif

#ifndef _STDIO_H
#include <stdio.h>
#ifndef _STDIO_H
#define _STDIO_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include <lal/LALStdlib.h>
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _AVFACTORIES_H
#include <lal/AVFactories.h>
#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H
#endif
#endif

NRCSID( PRINTVECTORTESTC, "$Id$" );

INT4 lalDebugLevel = 3;

int main( void )
{
  LALStatus                status = {0};

  COMPLEX16Sequence   *zSequence;
  COMPLEX8Sequence    *cSequence;
  REAL8Sequence       *dSequence;
  REAL4Sequence       *sSequence;

  COMPLEX16TimeSeries   zTimeSeries;
  COMPLEX8TimeSeries    cTimeSeries;
  REAL8TimeSeries       dTimeSeries;
  REAL4TimeSeries       sTimeSeries;

  COMPLEX16FrequencySeries   zFrequencySeries;
  COMPLEX8FrequencySeries    cFrequencySeries;
  REAL8FrequencySeries       dFrequencySeries;
  REAL4FrequencySeries       sFrequencySeries;

  COMPLEX16        *z;
  COMPLEX8         *c;
  REAL8            *d;
  REAL4            *s;

  INT2             n;
  LIGOTimeGPS      t80;
  LIGOTimeGPS      t00;
  LIGOTimeGPS      t10;

  t80.gpsSeconds = 0;
  t80.gpsNanoSeconds = 0;

  t00.gpsSeconds = 3600 * 24 * (15 * 365 + 5 * 366);
  t00.gpsNanoSeconds = 0;

  t10.gpsSeconds = 3600 * 24 * (22 * 365 + 8 * 366);
  t10.gpsNanoSeconds = 0;

  zSequence = NULL;

  LALZCreateVector( &status, &zSequence, 8 );
  for ( n=zSequence->length, z=zSequence->data; n > 0 ; --n, ++z ) {
    z->re = sinh(90.0*(4-n));
    z->im = - 1 / (1e-300 + z->re);
  }
  zTimeSeries.deltaT = 1e-3;
  zTimeSeries.epoch = t80;
  zTimeSeries.data = zSequence;

  LALZPrintTimeSeries(&zTimeSeries, "zTS.dat");

  zFrequencySeries.deltaF = 1;
  zFrequencySeries.epoch = t80;
  zFrequencySeries.f0 = 0;
  zFrequencySeries.data = zSequence;

  LALZPrintFrequencySeries(&zFrequencySeries, "zFS.dat");

  cSequence = NULL;

  LALCCreateVector( &status, &cSequence, 8 );
  for ( n=cSequence->length, c=cSequence->data; n > 0 ; --n, ++c ) {
    c->re = sinh(9.0*(4-n));
    c->im = - 1 / (1e-30 + c->re);
  } 
  cTimeSeries.deltaT = 1.0/1024.0;
  cTimeSeries.epoch = t00;
  cTimeSeries.data = cSequence;
 
  LALCPrintTimeSeries(&cTimeSeries, "cTS.dat");

  cFrequencySeries.deltaF = 1;
  cFrequencySeries.epoch = t80;
  cFrequencySeries.f0 = 0;
  cFrequencySeries.data = cSequence;

  LALCPrintFrequencySeries(&cFrequencySeries, "cFS.dat");

  dSequence = NULL;

  LALDCreateVector( &status, &dSequence, 8 );
  for ( n=dSequence->length, d=dSequence->data; n > 0 ; --n, ++d ) {
    *d = sinh(90.0*(4-n));
  }
  dTimeSeries.deltaT = 1.0/1024.0;
  dTimeSeries.epoch = t00;
  dTimeSeries.data = dSequence;

  LALDPrintTimeSeries(&dTimeSeries, "dTS.dat");
  for ( n=dSequence->length, d=dSequence->data; n > 0 ; --n, ++d ) {
    *d = 1 / (1e-300 + *d);
  }
  dFrequencySeries.f0 = -3*128;
  dFrequencySeries.deltaF = 128;
  dFrequencySeries.epoch = t00;
  dFrequencySeries.data = dSequence;
  LALDPrintFrequencySeries(&dFrequencySeries, "dFS.dat");

  sSequence = NULL;

  LALSCreateVector( &status, &sSequence, 8 );
  for ( n=sSequence->length, s=sSequence->data; n > 0 ; --n, ++s ) {
    *s = sinh(9.0*(4-n));
  }
  sTimeSeries.deltaT = 1.0/1024.0;
  sTimeSeries.epoch = t10;
  sTimeSeries.data = sSequence;
  LALSPrintTimeSeries(&sTimeSeries, "sTS.dat");

  for ( n=sSequence->length, s=sSequence->data; n > 0 ; --n, ++s ) {
    *s = 1 / (1e-30 + *s);
  }
  sFrequencySeries.f0 = 0;
  sFrequencySeries.deltaF = 128;
  sFrequencySeries.epoch = t10;
  sFrequencySeries.data = sSequence;
  LALSPrintFrequencySeries(&sFrequencySeries, "sFS.dat");

  return 0;
}
