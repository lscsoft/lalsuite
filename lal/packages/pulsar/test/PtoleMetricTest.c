/************************************ <lalVerbatim file="PtoleMetricTestCV">
Author: Owen, B. J.
$Id$
********************************************************** </lalVerbatim> */

/**************************************************************** <lalLaTeX>

\subsection{Program \texttt{PtoleMetricTest}}
\label{ss:PtoleMetricTest}

Tests the \texttt{PtoleMetric()} function.

\subsubsection*{Usage}
\begin{verbatim}
PtoleMetricTest
\end{verbatim}

\subsubsection*{Description}

This program computes examples of metric components for two points in
parameter space. With \texttt{lalDebugLevel} nonzero, it showcases the
error checking in \texttt{PtoleMetric()}. The ordering of the components
shown is $(f_0, \alpha, \delta, f_1, \ldots)$ for the unprojected metric,
and $(\alpha, \delta, f_1, \ldots)$ for the projected metric.

If \texttt{lalDebugLevel} is nonzero, the program will showcase a list of
error traps in \texttt{PtoleMetric()}. If not, only some valid results are
displayed.

\subsubsection*{Exit Codes}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\begin{verbatim}
lalDebugLevel
LALProjectMetric()
LALPtoleMetric()
\end{verbatim}

\subsubsection*{Notes}

If you want to see a plot and have \texttt{xmgrace} installed on your
system, search for XMGRACE in the source code and uncomment the big chunk of
code following it. The graph shows that the patches' overall area is
independent of right ascension but that those near the equator rotate,
which adds a new complication to the tiling.

\vfill{\footnotesize\input{PtoleMetricTestCV}}

************************************************************* </lalLaTeX> */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <lal/AVFactories.h>
#include <lal/LALMalloc.h>
#include <lal/LALStdlib.h>
#include <lal/PtoleMetric.h>
#include <lal/StackMetric.h>

NRCSID( PTOLEMETRICTESTC, "$Id" );

#define NUM_SPINDOWN 1 /* Number of spindown parameters */
#define SPOKES 30

int lalDebugLevel = 0;

int main( int argc, char *argv[] ) {
  static LALStatus status;          /* Status structure */
  PtoleMetricIn    in;              /* PtoleMetric() input structure */
  REAL4Vector     *spindown = NULL; /* Spindown parameters */
  REAL8Vector     *metric = NULL;   /* Parameter-space metric */
  INT4             j, k;            /* Loop counters */
  int ra, dec, i;
  FILE *pvc;

  if (lalDebugLevel) {
    printf("\nTesting bad I/O structures...\n");
    LALPtoleMetric( &status, metric, NULL );
    printf("\nTesting bad sky position...\n");
    in.position.longitude = 1e10;
    LALPtoleMetric( &status, metric, &in );
  }
  /* Good sky position: ra=19h13m, dec=+16deg, sound familiar? */
  in.position.system = COORDINATESYSTEM_EQUATORIAL;
  in.position.longitude = 288.25*LAL_PI/180;
  in.position.latitude = 16*LAL_PI/180;

  if (lalDebugLevel){
    printf("\nTesting bad spindown parameters...\n");
    LALPtoleMetric( &status, metric, &in );
  }
  /* Good spindown parameters: all zero */
  LALCreateVector( &status, &spindown, NUM_SPINDOWN );
  for (j=0; j<spindown->length; j++)
    spindown->data[j] = 0;
  in.spindown = spindown;

  /* Good epoch: simple as possible, need to find those constants */
  in.epoch.gpsSeconds = 0;
  in.epoch.gpsNanoSeconds = 0;

  if (lalDebugLevel) {
    printf("\nTesting bad duration...\n");
    in.duration = -1;
    LALPtoleMetric( &status, metric, &in );
  }
  /* Good duration: typical first-stage coherent integration */
  in.duration = 1e5;

  if (lalDebugLevel) {
    printf("\nTesting bad maximum frequency...\n");
    in.maxFreq = 0;
    LALPtoleMetric( &status, metric, &in );
  }
  /* Good maximum frequency */
  in.maxFreq = 1e3;

  if (lalDebugLevel) {
    printf("\nTesting bad detector site...\n");
    in.site.vertexLatitudeDegrees = 100;
    LALPtoleMetric( &status, metric, &in );
  }
  /* Use GEO600 site. */
  in.site = lalCachedDetectors[LALDetectorIndexGEO600DIFF].frDetector;

  if (lalDebugLevel) {
    printf("\nTesting bad output contents...\n");
    LALPtoleMetric( &status, metric, &in );
  }
  /* Allocate storage for output metric. */
  LALDCreateVector( &status, &metric, (4+NUM_SPINDOWN)*(5+NUM_SPINDOWN)/2 );

  printf( "\nValid results for duration 1e5 seconds:\n" );
  in.duration = 1e5;
  LALPtoleMetric( &status, metric, &in );
  for (j=0; j<=2+NUM_SPINDOWN; j++) {
    for (k=0; k<=j; k++)
      printf( "  %+.3e", metric->data[k+j*(j+1)/2] );
    printf("\n");
  }
  LALProjectMetric( &status, metric, 0 );
  printf( "With f0 projected out:\n" );
  for (j=1; j<=2+NUM_SPINDOWN; j++) {
    for (k=1; k<=j; k++)
      printf( "  %+.3e", metric->data[k+j*(j+1)/2] );
    printf( "\n" );
  }

  printf( "\nValid results for duration 1e7 seconds:\n" );
  in.duration = 1e7;
  LALPtoleMetric( &status, metric, &in );
  for (j=0; j<=2+NUM_SPINDOWN; j++) {
    for (k=0; k<=j; k++)
      printf( "  %+.3e", metric->data[k+j*(j+1)/2] );
    printf("\n");
  }
  printf( "With f0 projected out:\n" );
  for (j=1; j<=2+NUM_SPINDOWN; j++) {
    for (k=1; k<=j; k++)
      printf( "  %+.3e", metric->data[k+j*(j+1)/2] );
    printf( "\n" );
  }

/* If you want xmgrace plots, then uncomment the next line and recompile: */
/* #define XMGRACE 1 */
#if XMGRACE
  /* Take care of header stuff. */
  in.duration = 1e5;
  pvc = fopen( "pvc", "w" );
  fprintf( pvc, "@xaxis label \"Right ascension (degrees)\"\n" );
  fprintf( pvc, "@yaxis label \"Declination (degrees)\"\n" );
  /* Step around the sky. */
  j = 0;
  for (dec=80; dec>0; dec-=10) {
    for (ra=0; ra<=90; ra+=15) {
      float gaa, gad, gdd, angle, smaj, smin;
 
      /* Get the metric at this ra, dec. */
      in.position.longitude = ra*LAL_PI_180;
      in.position.latitude = dec*LAL_PI_180;
      LALPtoleMetric( &status, metric, &in );
      /* Rename \gamma_{\alpha\alpha}. */
      gaa = metric->data[2];
      /* Rename \gamma_{\alpha\delta}. */
      gad = metric->data[4];
      /* Rename \gamma_{\delta\delta}. */
      gdd = metric->data[5];
      /* Larger eigenvalue of metric. */
      smin = gaa+gdd + sqrt( pow(gaa-gdd,2) + pow(2*gad,2) );
      smin = sqrt(.02/smin);
      /* Smaller eigenvalue of metric. */
      smaj = gaa+gdd - sqrt( pow(gaa-gdd,2) + pow(2*gad,2) );
      smaj = sqrt(.02/smaj);
      /* Angle of semimajor axis with "horizontal" (equator). */
      angle = atan2( gad, .02/smaj/smaj-gdd );
      if (angle <= -M_PI/2) angle += M_PI;
      if (angle > M_PI/2) angle -= M_PI;
 
      /* Print set header. */
      fprintf( pvc, "@s%d color (0,0,0)\n", j );
      fprintf( pvc, "@target G0.S%d\n@type xy\n", j++ );
      /* Print center of patch. */
      fprintf( pvc, "%16.8g %16.8g\n", (float)ra, (float)dec );
      /* Loop around patch ellipse. */
      for (i=0; i<=SPOKES; i++) {
        float c, r;
        c = 2*M_PI*i/SPOKES;
        r = 10*LAL_180_PI*smaj*smin/sqrt( pow(smaj*sin(c),2) + pow(smin*cos(c),2) );
        fprintf( pvc, "%e %e\n", ra+r*cos(angle-c), dec+r*sin(angle-c) );
      } /* for (a...) */
 
    } /* for (ra...) */
  } /* for (dec...) */
  fclose( pvc );
  system( "xmgrace pvc &" );
  system( "sleep 2; rm -f pvc" );
#endif

  printf("\nCleaning up and leaving...\n");
  LALDestroyVector( &status, &spindown );
  LALDDestroyVector( &status, &metric );
  LALCheckMemoryLeaks();

  /* Get while the getting's good. */
  return 0;
} /* main() */
