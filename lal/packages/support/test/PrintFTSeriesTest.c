/*************** <lalVerbatim file="PrintFTSeriesTestCV"> *******
Author: Whelan, J. T.
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

#include <lal/Units.h>
#include <lal/PrintFTSeries.h>
#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>

NRCSID( PRINTVECTORTESTC, "$Id$" );

INT4 lalDebugLevel = 3;

int main( void )
{
  static LALStatus status;

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
  int i;

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
  strncpy(zTimeSeries.name,"Complex time series",LALNameLength);
  zTimeSeries.sampleUnits = lalDimensionlessUnit;
  zTimeSeries.deltaT = 1e-3;
  zTimeSeries.f0 = 0;
  zTimeSeries.epoch = t80;
  zTimeSeries.data = zSequence;

  LALZPrintTimeSeries(&zTimeSeries, "zTS.dat");

  strncpy(zFrequencySeries.name,"Complex frequency series",LALNameLength);
  zFrequencySeries.sampleUnits = lalDimensionlessUnit; 
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
  strncpy(cTimeSeries.name,"Complex time series",LALNameLength);
  cTimeSeries.sampleUnits = lalDimensionlessUnit;
  cTimeSeries.deltaT = 1.0/1024.0;
  cTimeSeries.f0 = 0;
  cTimeSeries.epoch = t00;
  cTimeSeries.data = cSequence;
 
  LALCPrintTimeSeries(&cTimeSeries, "cTS.dat");

  strncpy(cFrequencySeries.name,"Complex frequency series",LALNameLength);
  cFrequencySeries.sampleUnits = lalDimensionlessUnit;
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
  strncpy(dTimeSeries.name,"Real time series",LALNameLength);
  dTimeSeries.sampleUnits = lalDimensionlessUnit;
  dTimeSeries.f0 = 0;
  dTimeSeries.deltaT = 1.0/1024.0;
  dTimeSeries.epoch = t00;
  dTimeSeries.data = dSequence;

  LALDPrintTimeSeries(&dTimeSeries, "dTS.dat");
/*    for ( n=dSequence->length, d=dSequence->data; n > 0 ; --n, ++d ) { */
/*      *d = 1 / (1e-300 + *d); */
/*    } */
  for ( n=dSequence->length, d=dSequence->data; n > 0 ; --n, ++d ) {
    *d = 1 / (1e-300 + *d);
  }
  strncpy(dFrequencySeries.name,"Real frequency series",LALNameLength);
  dFrequencySeries.sampleUnits = lalDimensionlessUnit;
  dFrequencySeries.f0 = 0     ;
  dFrequencySeries.deltaF = 128;
  dFrequencySeries.epoch = t00;
  dFrequencySeries.data = dSequence;
  LALDPrintFrequencySeries(&dFrequencySeries, "dFS.dat");

  sSequence = NULL;

  LALSCreateVector( &status, &sSequence, 8 );
  for ( n=sSequence->length, s=sSequence->data; n > 0 ; --n, ++s ) {
    *s = sinh(9.0*(4-n));
  }
  strncpy(sFrequencySeries.name,"Real frequency series",LALNameLength);
  sTimeSeries.sampleUnits = lalDimensionlessUnit;
  sTimeSeries.deltaT = 1.0/1024.0;
  sTimeSeries.f0 = 0;
  sTimeSeries.epoch = t10;
  sTimeSeries.data = sSequence;
  LALSPrintTimeSeries(&sTimeSeries, "sTS.dat");

  for ( n=sSequence->length, s=sSequence->data; n > 0 ; --n, ++s ) {
    *s = 1 / (1e-30 + *s);
  }
  sFrequencySeries.sampleUnits = lalDimensionlessUnit;
  sFrequencySeries.f0 = 0;
  sFrequencySeries.deltaF = 128;
  sFrequencySeries.epoch = t10;
  sFrequencySeries.data = sSequence;
  LALSPrintFrequencySeries(&sFrequencySeries, "sFS.dat");

  strncpy(cFrequencySeries.name,"Heterodyned frequency series",LALNameLength);
  cFrequencySeries.sampleUnits = lalDimensionlessUnit;
  cFrequencySeries.f0 = 500.0;
  cFrequencySeries.deltaF = 50.0;
  cFrequencySeries.epoch = t00;
  cFrequencySeries.data = cSequence;

  for (i=0; i<(int)cSequence->length; ++i) {
    cSequence->data[i].re = 1.0*i;
    cSequence->data[i].im = cFrequencySeries.f0 + cFrequencySeries.deltaF * i;
  }

  LALCPrintFrequencySeries(&cFrequencySeries, "hFSe.dat");

  LALCDestroyVector( &status, &cSequence );
  LALCCreateVector( &status, &cSequence, 9 );

  cFrequencySeries.sampleUnits = lalDimensionlessUnit;
  cFrequencySeries.f0 = 500.0;
  cFrequencySeries.deltaF = 50.0;
  cFrequencySeries.epoch = t00;
  cFrequencySeries.data = cSequence;

  for (i=0; i<(int)cSequence->length; ++i) {
    cSequence->data[i].re = 1.0*i;
    cSequence->data[i].im = cFrequencySeries.f0 + cFrequencySeries.deltaF*i;
  }

  LALCPrintFrequencySeries(&cFrequencySeries, "hFSo.dat");

  return 0;
}
