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

This program computes metric components using \texttt{PtoleMetric()}. The
ordering of the components is $(f_0, \alpha, \delta, f_1, \ldots)$ for the
unprojected metric, and $(\alpha, \delta, f_1, \ldots)$ for the metric with
$f_0$ projected out.

With no options, this program displays components for a single point in
parameter space for two different integration times.

The \texttt{-t} option causes the program to showcase error messages when
given bad parameter values, etc.

The \texttt{-x} option produces a graph of the 2% power mismatch contours on
the sky. Dimensions of the ellipses have been exaggerated by 10 for
legibility. The purpose of the graph is to get a feel for how ellipses are
varying across parameter space. Note that this option makes a system call to
the xmgrace program, and will not work if that program is not installed on
your system.

\subsubsection*{Exit Codes}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\begin{verbatim}
lalDebugLevel
LALCheckMemoryLeaks()
LALCreateVector()
LALDestroyVector()
LALDCreateVector()
LALDDestroyVector()
LALProjectMetric()
LALPtoleMetric()
\end{verbatim}

\subsubsection*{Notes}

The graph shows that the patches' overall area is independent of right
ascension but that those near the equator rotate, which adds a new
complication to the tiling.

\vfill{\footnotesize\input{PtoleMetricTestCV}}

************************************************************* </lalLaTeX> */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
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
  int              j, k;            /* Parameter-space indices */
  int              opt;             /* Command-line option. */
  BOOLEAN          test = 0;        /* Whether we showcase error messages */
  BOOLEAN          grace = 0;       /* Whether or not we use xmgrace */
  int              ra, dec, i;      /* Loop variables for xmgrace option */
  FILE            *pvc;             /* Temporary file for xmgrace option */

  /* Parse options. */
  while ((opt = getopt( argc, argv, "tx" )) != -1) {
    switch (opt) {
    case 't':
      test = 1;
      lalDebugLevel = 1;
      break;
    case 'x':
      grace = 1;
      break;
    }
  }

  if (test) {
    printf("\nTesting bad I/O structures...\n");
    LALPtoleMetric( &status, metric, NULL );
    printf("\nTesting bad sky position...\n");
    in.position.longitude = 1e10;
    LALPtoleMetric( &status, metric, &in );
  }
  /* Good sky position: ra=19h13m, dec=+16deg, sound familiar? */
  in.position.system = COORDINATESYSTEM_EQUATORIAL;
  in.position.longitude = 288.25*LAL_PI_180;
  in.position.latitude = 16*LAL_PI_180;

  if (test){
    printf("\nTesting bad spindown parameters...\n");
    LALPtoleMetric( &status, metric, &in );
  }
  /* Good spindown parameters: all zero */
  LALCreateVector( &status, &spindown, NUM_SPINDOWN );
  for (j=0; (UINT4)j<spindown->length; j++)
    spindown->data[j] = 0;
  in.spindown = spindown;

  /* Good epoch: simple as possible, need to find those constants */
  in.epoch.gpsSeconds = 0;
  in.epoch.gpsNanoSeconds = 0;

  if (test) {
    printf("\nTesting bad duration...\n");
    in.duration = -1;
    LALPtoleMetric( &status, metric, &in );
  }
  /* Good duration: typical first-stage coherent integration */
  in.duration = 1e5;

  if (test) {
    printf("\nTesting bad maximum frequency...\n");
    in.maxFreq = 0;
    LALPtoleMetric( &status, metric, &in );
  }
  /* Good maximum frequency */
  in.maxFreq = 1e3;

  if (test) {
    printf("\nTesting bad detector site...\n");
    in.site.vertexLatitudeDegrees = 100;
    LALPtoleMetric( &status, metric, &in );
  }
  /* Use GEO600 site. */
  in.site = lalCachedDetectors[LALDetectorIndexGEO600DIFF].frDetector;

  if (test) {
    printf("\nTesting bad output contents...\n");
    LALPtoleMetric( &status, metric, &in );
  }
  /* Allocate storage for output metric. */
  LALDCreateVector( &status, &metric, (4+NUM_SPINDOWN)*(5+NUM_SPINDOWN)/2 );

  /* Print results if no options. */
  if (argc == 1) {

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
  } /* if (argc...) */

  /* Here is the code that uses xmgrace with the -x option. */
  if (grace) {

    /* Take care of preliminaries. */
    in.duration = 1e5;
    pvc = fopen( "pvc", "w" );
    fprintf( pvc, "@xaxis label \"Right ascension (degrees)\"\n" );
    fprintf( pvc, "@yaxis label \"Declination (degrees)\"\n" );

    /* Step around the sky: a grid in ra and dec. */
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
        /* Semiminor axis from larger eigenvalue of metric. */
        smin = gaa+gdd + sqrt( pow(gaa-gdd,2) + pow(2*gad,2) );
        smin = sqrt(.02/smin);
        /* Semiminor axis from smaller eigenvalue of metric. */
        smaj = gaa+gdd - sqrt( pow(gaa-gdd,2) + pow(2*gad,2) );
        smaj = sqrt(.02/smaj);
        /* Angle of semimajor axis with "horizontal" (equator). */
        angle = atan2( gad, .02/smaj/smaj-gdd );
        if (angle <= -LAL_PI_2) angle += LAL_PI;
        if (angle > LAL_PI_2) angle -= LAL_PI;
 
        /* Print set header. */
        fprintf( pvc, "@s%d color (0,0,0)\n", j );
        fprintf( pvc, "@target G0.S%d\n@type xy\n", j++ );
        /* Print center of patch. */
        fprintf( pvc, "%16.8g %16.8g\n", (float)ra, (float)dec );
        /* Loop around patch ellipse. */
        for (i=0; i<=SPOKES; i++) {
          float c, r;
          c = LAL_TWOPI*i/SPOKES;
          r = 10*LAL_180_PI*smaj*smin/sqrt( pow(smaj*sin(c),2)
              + pow(smin*cos(c),2) );
          fprintf( pvc, "%e %e\n", ra+r*cos(angle-c), dec+r*sin(angle-c) );
        } /* for (a...) */
 
      } /* for (ra...) */
    } /* for (dec...) */
    fclose( pvc );
    /* Yes, system() calls are dirty, but we're only doing two. */
    system( "xmgrace pvc &" );
    system( "sleep 2; rm -f pvc" );
  } /* if (grace) */

  printf("\nCleaning up and leaving...\n");
  LALDestroyVector( &status, &spindown );
  LALDDestroyVector( &status, &metric );
  LALCheckMemoryLeaks();

  /* Get while the getting's good. */
  return 0;
} /* main() */
