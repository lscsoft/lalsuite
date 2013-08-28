/*
*  Copyright (C) 2007 Jolien Creighton, Ian Jones, Benjamin Owen
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**
 * \author Owen, B. J.,   Jones, D. I.
 * \file
 * \ingroup PtoleMetric_h
 * \brief Tests the LALPtoleMetric() function.
 *
 * ### Program \c PtoleMetricTest ###
 *
 *
 * ### Usage ###
 *
 * \code
 * PtoleMetricTest
 * \endcode
 *
 * ### Description ###
 *
 * This program computes metric components using LALPtoleMetric(). The
 * ordering of the components is \f$(f_0, \alpha, \delta, f_1, \ldots)\f$ for the
 * unprojected metric, and \f$(\alpha, \delta, f_1, \ldots)\f$ for the metric with
 * \f$f_0\f$ projected out.
 *
 * With no options, this program displays metric components for a single point
 * in parameter space for the default integration time (see the <tt>-t</tt>
 * option).
 *
 * The <tt>-b</tt> option sets the beginning time of integration to the option
 * argument. (Default is \f$0\f$ seconds)
 *
 * The <tt>-e</tt> option causes the program to showcase error messages when
 * given bad parameter values, etc.
 *
 * The <tt>-m</tt> option sets the mismatch (default is \f$0.02\f$).
 *
 * The <tt>-p</tt> option is provided for users who wish to view the
 * power mismatch contours provided by the <tt>-x</tt> option (see below)
 * but don't have xmgrace installed.  All necessary data is simply
 * written to a file ``nongrace.data''; it's probably best to look at the
 * code to see the exact format.  The user should then write a small
 * script to convert the data into the format appropriate to their
 * favorite graphics package.
 *
 * The <tt>-t</tt> option sets the duration of integration in seconds. The default
 * is \f$10^5\f$ seconds, which is chosen because it is about one day but not an
 * integer multiple of any astronomically important timescale.
 *
 * The <tt>-x</tt> option produces a graph of the 2\% power mismatch
 * contours on the sky. Dimensions of the ellipses have been exaggerated
 * by a factor of \c MAGNIFY (specified within the code) for
 * legibility. The purpose of the graph is to get a feel for how ellipses
 * are varying across parameter space. Note that this option makes a
 * system call to the \c xmgrace program, and will not work if that
 * program is not installed on your system.
 *
 * ### Algorithm ###
 *
 *
 * ### Uses ###
 *
 * \code
 * lalDebugLevel
 * LALCheckMemoryLeaks()
 * LALCreateVector()
 * LALDestroyVector()
 * LALDCreateVector()
 * LALDDestroyVector()
 * LALProjectMetric()
 * LALPtoleMetric()
 * xmgrace
 * \endcode
 *
 * ### Notes ###
 *
 * The graph shows that the patches' overall area is independent of right
 * ascension but that those near the equator rotate, which adds a new
 * complication to the tiling.
 *
 */

/** \name Error Codes */ /*@{*/
#define PTOLEMETRICTESTC_EMEM 1
#define PTOLEMETRICTESTC_ESUB 2
#define PTOLEMETRICTESTC_ESYS 3

#define PTOLEMETRICTESTC_MSGEMEM "memory (de)allocation error"
#define PTOLEMETRICTESTC_MSGESUB "subroutine failed"
#define PTOLEMETRICTESTC_MSGESYS "system call failed"
/*@}*/

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

/** \cond DONT_DOXYGEN */
extern char *optarg;

#define DEFAULT_DURATION 1e5 /* seconds */
#define NUM_SPINDOWN 0       /* Number of spindown parameters */
#define SPOKES 30
#define MAGNIFY 1.0           /* Magnification factor of ellipses */


int main( int argc, char *argv[] ) {
  static LALStatus status;          /* Status structure */
  PtoleMetricIn    in;              /* PtoleMetric() input structure */
  REAL4            mismatch;        /* mismatch threshold of mesh */
  REAL4Vector     *spindown = NULL; /* Spindown parameters */
  REAL8Vector     *metric = NULL;   /* Parameter-space metric */
  int              j, k;            /* Parameter-space indices */
  int              opt;             /* Command-line option. */
  BOOLEAN          test = 0;        /* Whether we showcase error messages */
  BOOLEAN          grace = 0;       /* Whether or not we use xmgrace */
  BOOLEAN          nongrace = 0;    /* Whether or not to output data to file*/
  int              ra, dec, i;      /* Loop variables for xmgrace option */
  FILE            *pvc=NULL;        /* Temporary file for xmgrace option */
  FILE            *fnongrace=NULL;  /* File contaning ellipse coordinates */


  /* Default values. */
  in.duration = DEFAULT_DURATION;
  in.epoch.gpsSeconds = 0.0;
  in.epoch.gpsNanoSeconds = 0.0;
  mismatch = 0.02;


  /* Parse options. */
  while ((opt = getopt( argc, argv, "b:em:pt:x" )) != -1) {
    switch (opt) {
    case 'b':
      in.epoch.gpsSeconds = atof( optarg );
      break;
    case 'e':
      test = 1;
      break;
    case 'm':
      mismatch = atof( optarg );
      break;
    case 'p':
      nongrace = 1;
      break;
    case 't':
      in.duration = atof( optarg );
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
  /* Good spindown parameters: all zero (if we have any) */
  if( NUM_SPINDOWN > 0 )
  {
    LALCreateVector( &status, &spindown, NUM_SPINDOWN );
    if( status.statusCode )
    {
      printf( "%s line %d: %s\n", __FILE__, __LINE__,
              PTOLEMETRICTESTC_MSGEMEM );
      return PTOLEMETRICTESTC_EMEM;
    }
    for (j=0; (UINT4)j<spindown->length; j++)
      spindown->data[j] = 0;
    in.spindown = spindown;
  }
  else
    in.spindown = NULL;



  if (test) {
    REAL4 old_duration = in.duration;
    printf("\nTesting bad duration...\n");
    in.duration = -1;
    LALPtoleMetric( &status, metric, &in );
    in.duration = old_duration;
  }

  if (test) {
    printf("\nTesting bad maximum frequency...\n");
    in.maxFreq = 0;
    LALPtoleMetric( &status, metric, &in );
  }
  /* Good maximum frequency */
  in.maxFreq = 1e3;

  /* JC: DISABLE THIS
  if (test) {
    printf("\nTesting bad detector site...\n");
    in.site->frDetector.vertexLatitudeRadians = 100 * LAL_PI / 180;
    LALPtoleMetric( &status, metric, &in );
  }
  */
  /* Use GEO600 site. */
  in.site = &lalCachedDetectors[LALDetectorIndexGEO600DIFF];

  if (test) {
    printf("\nTesting bad output contents...\n");
    LALPtoleMetric( &status, metric, &in );
  }
  /* Allocate storage for output metric. */
  LALDCreateVector( &status, &metric, (3+NUM_SPINDOWN)*(4+NUM_SPINDOWN)/2 );
  if( status.statusCode )
  {
    printf( "%s line %d: %s\n", __FILE__, __LINE__, PTOLEMETRICTESTC_MSGEMEM );
    return PTOLEMETRICTESTC_EMEM;
  }

  /* Print results if no options. */
  if (argc == 1) {

    printf( "\nValid results for duration %e seconds:\n", in.duration );
    LALPtoleMetric( &status, metric, &in );
    if( status.statusCode )
    {
      printf( "%s line %d: %s\n", __FILE__, __LINE__,
              PTOLEMETRICTESTC_MSGESUB );
      return PTOLEMETRICTESTC_ESUB;
    }
    for (j=0; j<=2+NUM_SPINDOWN; j++) {
      for (k=0; k<=j; k++)
        printf( "  %+.3e", metric->data[k+j*(j+1)/2] );
      printf("\n");
    }
    LALProjectMetric( &status, metric, 0 );
    if( status.statusCode )
    {
      printf( "%s line %d: %s\n", __FILE__, __LINE__,
              PTOLEMETRICTESTC_MSGESUB );
      return PTOLEMETRICTESTC_ESUB;
    }
    printf( "With f0 projected out:\n" );
    for (j=1; j<=2+NUM_SPINDOWN; j++) {
      for (k=1; k<=j; k++)
        printf( "  %+.3e", metric->data[k+j*(j+1)/2] );
      printf( "\n" );
    }

#if 0
    printf( "\nValid results for duration 1e7 seconds:\n" );
    in.duration = 1e7;
    LALPtoleMetric( &status, metric, &in );
    if( status.statusCode )
    {
      printf( "%s line %d: %s\n", __FILE__, __LINE__,
              PTOLEMETRICTESTC_MSGESUB );
      return PTOLEMETRICTESTC_ESUB;
    }
    for (j=0; j<=2+NUM_SPINDOWN; j++) {
      for (k=0; k<=j; k++)
        printf( "  %+.3e", metric->data[k+j*(j+1)/2] );
      printf("\n");
    }
    LALProjectMetric( &status, metric, 0 );
    if( status.statusCode )
    {
      printf( "%s line %d: %s\n", __FILE__, __LINE__,
              PTOLEMETRICTESTC_MSGESUB );
      return PTOLEMETRICTESTC_ESUB;
    }
    printf( "With f0 projected out:\n" );
    for (j=1; j<=2+NUM_SPINDOWN; j++) {
      for (k=1; k<=j; k++)
        printf( "  %+.3e", metric->data[k+j*(j+1)/2] );
      printf( "\n" );
    }
#endif
  } /* if (argc...) */

  /* Here is the code that uses xmgrace with the -x option, */
  /* and outputs data to a file with the -t option. */
  if (grace || nongrace) {

    /* Take care of preliminaries. */
    if(grace)
      {
	pvc = popen( "xmgrace -pipe", "w" );
	if( !pvc )
	  {
	    printf( "%s line %d: %s\n", __FILE__, __LINE__,
		    PTOLEMETRICTESTC_MSGESYS );
	    return PTOLEMETRICTESTC_ESYS;
	  }
	fprintf( pvc, "@xaxis label \"Right ascension (degrees)\"\n" );
	fprintf( pvc, "@yaxis label \"Declination (degrees)\"\n" );
      }
    if(nongrace)
      {
	fnongrace = fopen( "nongrace.data", "w" );
	if( !fnongrace )
	  {
	    printf( "%s line %d: %s\n", __FILE__, __LINE__,
		    PTOLEMETRICTESTC_MSGESYS );
	    return PTOLEMETRICTESTC_ESYS;
	  }
      }

    /* Step around the sky: a grid in ra and dec. */
    j = 0;
    for (dec=80; dec>0; dec-=10) {
      for (ra=0; ra<=90; ra+=15) {
        float gaa, gad, gdd, angle, smaj, smin;

        /* Get the metric at this ra, dec. */
        in.position.longitude = ra*LAL_PI_180;
        in.position.latitude = dec*LAL_PI_180;
        LALPtoleMetric( &status, metric, &in );
        if( status.statusCode )
        {
          printf( "%s line %d: %s\n", __FILE__, __LINE__,
                  PTOLEMETRICTESTC_MSGESUB );
          return PTOLEMETRICTESTC_ESUB;
        }

	/*  Project metric: */
	LALProjectMetric( &status, metric, 0 );
	if( status.statusCode )
	  {
          printf( "%s line %d: %s\n", __FILE__, __LINE__,
                  PTOLEMETRICTESTC_MSGESUB );
          return PTOLEMETRICTESTC_ESUB;
	  }

        /* Rename \gamma_{\alpha\alpha}. */
        gaa = metric->data[2];
        /* Rename \gamma_{\alpha\delta}. */
        gad = metric->data[4];
        /* Rename \gamma_{\delta\delta}. */
        gdd = metric->data[5];
        /* Semiminor axis from larger eigenvalue of metric. */
        smin = gaa+gdd + sqrt( pow(gaa-gdd,2) + pow(2*gad,2) );
        smin = sqrt(2*mismatch/smin);
        /* Semiminor axis from smaller eigenvalue of metric. */
        smaj = gaa+gdd - sqrt( pow(gaa-gdd,2) + pow(2*gad,2) );
        smaj = sqrt(2*mismatch/smaj);
        /* Angle of semimajor axis with "horizontal" (equator). */
        angle = atan2( gad, mismatch/smaj/smaj-gdd );
        if (angle <= -LAL_PI_2) angle += LAL_PI;
        if (angle > LAL_PI_2) angle -= LAL_PI;

        if(grace)
	  {
	    /* Print set header. */
	    fprintf( pvc, "@s%d color (0,0,0)\n", j );
	    fprintf( pvc, "@target G0.S%d\n@type xy\n", j++ );
	    /* Print center of patch. */
	    fprintf( pvc, "%16.8g %16.8g\n", (float)ra, (float)dec );
	  }
	if(nongrace)
	  /* Print center of patch. */
	  fprintf( fnongrace, "%16.8g %16.8g\n", (float)ra, (float)dec );
	/* Loop around patch ellipse. */
        for (i=0; i<=SPOKES; i++) {
          float c, r;
          c = LAL_TWOPI*i/SPOKES;
          r = MAGNIFY*LAL_180_PI*smaj*smin/sqrt( pow(smaj*sin(c),2)
              + pow(smin*cos(c),2) );
	  if(grace)
	    fprintf( pvc, "%e %e\n", ra+r*cos(angle-c), dec+r*sin(angle-c) );
	  if(nongrace)
	    fprintf( fnongrace, "%e %e\n", ra+r*cos(angle-c),
		     dec+r*sin(angle-c) );

        } /* for (a...) */

      } /* for (ra...) */
    } /* for (dec...) */
    if(grace)
      fclose( pvc );
    if(nongrace)
      fclose( fnongrace );
  } /* if (grace || nongrace) */

  printf("\nCleaning up and leaving...\n");
  if( spindown )
  {
    LALDestroyVector( &status, &spindown );
    if( status.statusCode )
    {
      printf( "%s line %d: %s\n", __FILE__, __LINE__,
              PTOLEMETRICTESTC_MSGEMEM );
      return PTOLEMETRICTESTC_EMEM;
    }
  }
  LALDDestroyVector( &status, &metric );
  if( status.statusCode )
  {
    printf( "%s line %d: %s\n", __FILE__, __LINE__,
            PTOLEMETRICTESTC_MSGEMEM );
    return PTOLEMETRICTESTC_EMEM;
  }
  LALCheckMemoryLeaks();
  return 0;
} /* main() */

/** \endcond */
