/*
*  Copyright (C) 2007 David M. Whitbeck, Jolien Creighton, Ian Jones, Benjamin Owen
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \author Jones, D. I.,     Owen, B. J.
 * \file
 * \ingroup PtoleMetric_h
 *
 * \brief Tests the various LAL metric functions, by outputting the metric at a
 * point in parameter space, and also producing an array of ellipses of
 * constant mismatch.
 *
 * ### Program <tt>GeneralMetricTest.c</tt> ###
 *
 *
 * ### Usage ###
 *
 * \code
 * GeneralMetricTest
 * \endcode
 *
 * ### Description ###
 *
 * This program computes metric components using a metric function of the
 * user's specification.  The ordering of the components is \f$ (f_0,
 * \alpha, \delta, f_1\ldots) \f$ for the unprojected metric, and \f$ (\alpha,
 * \delta, f_1\ldots) \f$ for the metric with \f$ f_0 \f$ projected out.
 *
 * With no options, this program displays metric components for a single point
 * in parameter space for the default parameter values.
 *
 * The <b>-c</b> option determines the point on the sky where the metric
 * is evaluated.  This option is hard-coded to use equatorial coordinates
 * and the argument should be given in hh:mm:ss:dd:mm:ss format.
 * (Default is the center of the globular cluster 47 Tuc).
 *
 * The <b>-f</b> option sets the maximum frequency (in Hz) to search. (The
 * default is 1000.)
 *
 * The <b>-l</b> option determines the limits in right ascension and
 * declination of the rectangular region over which the mismatch contours
 * are computed.  The argument should be given in degrees as
 * %RA(min):%RA(max):dec(min):dec(max).  (The default is the octant of the
 * sky defined by \f$ 0 < \textrm{RA} < 90 \f$ and \f$ 0< \textrm{dec} < 85 \f$ ; this avoids the
 * coordinate singularity at the poles.)
 *
 * The <b>-m</b> option sets the mismatch (default is \f$ 0.02 \f$ ).
 *
 * The <b>-n</b> option sets the number of spindown parameters (default 0).
 *
 * The <b>-p</b> option is provided for users who wish to view the
 * power mismatch contours provided by the <b>-x</b> option (see below)
 * but don't have xmgrace installed.  All necessary data is simply
 * written to a file ``nongrace.data''; it's probably best to look at the
 * code to see the exact format.  The user should then write a small
 * script to convert the data into the format appropriate to their
 * favorite graphics package.
 *
 * The <b>-t</b> option sets the duration of integration in seconds. (The
 * default is \f$ 39600 \f$ seconds \f$ = 11 \f$ hours, which is chosen because it is of
 * the right size for S2 analyses).
 *
 * The <b>-x</b> option produces a graph of the 2\% power mismatch
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
 * The code does not yet really work with more than one spindown parameter.
 *
 */


/** \name Error Codes */ /** @{ */
#define GENERALMETRICTESTC_EMEM 1
#define GENERALMETRICTESTC_ESUB 2
#define GENERALMETRICTESTC_ESYS 3
#define GENERALMETRICTESTC_EMET 4

#define GENERALMETRICTESTC_MSGEMEM "memory (de)allocation error"
#define GENERALMETRICTESTC_MSGESUB "subroutine failed"
#define GENERALMETRICTESTC_MSGESYS "system call failed"
#define GENERALMETRICTESTC_MSGEMET "determinant of projected metric negative"

/** @} */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <lal/AVFactories.h>
#include <lal/LALgetopt.h>
#include <lal/LALMalloc.h>
#include <lal/LALStdlib.h>
#include <lal/PtoleMetric.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>

/** \cond DONT_DOXYGEN */

#define DEFAULT_DURATION 39600 /* seconds */
#define SPOKES 30
#define MAGNIFY 1.0            /* Magnification factor of ellipses */


int main( int argc, char *argv[] )
{
  static LALStatus status;          /* Status structure */
  PtoleMetricIn    in;              /* PtoleMetric() input structure */
  REAL8            mismatch;        /* mismatch threshold of mesh */
  REAL8Vector     *metric;          /* Parameter-space metric */
  int              j, k;            /* Parameter-space indices */
  int              opt;             /* Command-line option. */
  BOOLEAN          grace;           /* Whether or not we use xmgrace */
  BOOLEAN          nongrace;        /* Whether or not to output data to file*/
  int              ra, dec, i;      /* Loop variables for xmgrace option */
  FILE            *pvc = NULL;      /* Temporary file for xmgrace option */
  FILE            *fnongrace = NULL;/* File contaning ellipse coordinates */
  REAL8           ra_point;         /* RA at which metric is evaluated */
  REAL8           dec_point;        /* dec at which metric is evaluated */
  float           a, b, c, d, e, f; /* To input point in standard format */
  int             ra_min, ra_max;   /* Min and max RA for ellipse plot */
  int             dec_min, dec_max; /* Min and max dec for ellipse plot */
  float           c_ellipse;        /* Centers of ellipses */
  float           r_ellipse;        /* Radii of ellipses */
  REAL8           determinant;      /* Determinant of projected metric */
  REAL4           f0;               /* carrier frequency */
  UINT2           numSpindown;      /* Number of spindowns */

  /* Defaults that can be overwritten: */
  in.epoch.gpsSeconds = 731265908;
  in.epoch.gpsNanoSeconds = 0.0;
  mismatch = 0.02;
  nongrace = 0;
  in.duration = DEFAULT_DURATION;
  grace = 0;
  ra_point  = ( 24.1 / 60 ) * LAL_PI_180; /* 47 Tuc */
  dec_point = -( 72 + 5. / 60 ) * LAL_PI_180;
  ra_min = 0;
  ra_max = 90;
  dec_min = 0;
  dec_max = 85;
  f0 = 1000;
  numSpindown = 0;

  /* Parse options. */
  while ( ( opt = LALgetopt( argc, argv, "a:b:c:d:ef:l:m:n:pt:s:x" ) ) != -1 ) {
    switch ( opt ) {
    case 'b':
      in.epoch.gpsSeconds = atoi( LALoptarg );
      break;
    case 'c':
      if ( sscanf( LALoptarg, "%f:%f:%f:%f:%f:%f", &a, &b, &c, &d, &e, &f ) != 6 ) {
        fprintf( stderr, "coordinates should be hh:mm:ss:dd:mm:ss\n" );
      }
      ra_point = ( 15 * a + b / 4 + c / 240 ) * LAL_PI_180;
      dec_point = ( d + e / 60 + f / 3600 ) * LAL_PI_180;
      break;
    case 'f':
      f0 = atof( LALoptarg );
      break;
    case 'l':
      if ( sscanf( LALoptarg, "%d:%d:%d:%d",
                   &ra_min, &ra_max, &dec_min, &dec_max ) != 4 ) {
        fprintf( stderr, "coordinates should be ra_min, ra_max, dec_min, dec_max all in degrees" );
      }
      break;
    case 'm':
      mismatch = atof( LALoptarg );
      break;
    case 'n':
      numSpindown = atoi( LALoptarg );
      break;
    case 'p':
      nongrace = 1;
      break;
    case 's':
      break;
    case 'x':
      grace = 1;
      break;
    }
  }

  /* Allocate storage. */
  metric = NULL;
  LALDCreateVector( &status, &metric, ( 3 + numSpindown ) * ( 4 + numSpindown ) / 2 );
  if ( status.statusCode ) {
    printf( "%s line %d: %s\n", __FILE__, __LINE__,
            GENERALMETRICTESTC_MSGEMEM );
    return GENERALMETRICTESTC_EMEM;
  }

  /* Position in parameter space (sky, frequency, spindowns) */
  in.position.system = COORDINATESYSTEM_EQUATORIAL;
  in.position.longitude = ra_point;
  in.position.latitude = dec_point;
  in.maxFreq = f0;
  in.spindown = NULL;
  if ( numSpindown > 0 ) {
    LALCreateVector( &status, &( in.spindown ), numSpindown );
    if ( status.statusCode ) {
      printf( "%s line %d: %s\n", __FILE__, __LINE__,
              GENERALMETRICTESTC_MSGEMEM );
      return GENERALMETRICTESTC_EMEM;
    }
    for ( i = 0; i < numSpindown; i++ ) {
      in.spindown->data[i] = 0;
    }
  }

  /* Detector site */
  in.site = &lalCachedDetectors[LALDetectorIndexGEO600DIFF];

  /* Evaluate metric components. */
  LALPtoleMetric( &status, metric, &in );
  if ( status.statusCode ) {
    printf( "%s line %d: %s\n", __FILE__, __LINE__,
            GENERALMETRICTESTC_MSGESUB );
    return GENERALMETRICTESTC_ESUB;
  }

  /* Print metric. */
  printf( "\nmetric (f0, alpha, delta, ...) at the requested point\n" );
  for ( j = 0; j <= 2 + numSpindown; j++ ) {
    for ( k = 0; k <= j; k++ ) {
      printf( "  %+.4e", metric->data[k + j * ( j + 1 ) / 2] );
    }
    printf( "\n" );
  }

  /* Print determinants. */
  determinant = metric->data[5] * metric->data[2] - pow( metric->data[4], 2 );
  printf( "\nSky-determinant %e\n", determinant );
  if ( numSpindown == 1 ) {
    determinant = metric->data[2] * metric->data[5] * metric->data[9]
                  - metric->data[2] * metric->data[8] * metric->data[8]
                  + metric->data[4] * metric->data[8] * metric->data[7]
                  - metric->data[4] * metric->data[4] * metric->data[9]
                  + metric->data[7] * metric->data[4] * metric->data[8]
                  - metric->data[7] * metric->data[7] * metric->data[5];
    printf( "S&S determinant %e\n", determinant );
  }

  /* Project carrier frequency out of metric. */
  LALProjectMetric( &status, metric, 0 );
  if ( status.statusCode ) {
    printf( "%s line %d: %s\n", __FILE__, __LINE__,
            GENERALMETRICTESTC_MSGESUB );
    return GENERALMETRICTESTC_ESUB;
  }

  /* Print projected metric. */
  printf( "\nf-projected metric (alpha, delta, ...) at the requested point\n" );
  for ( j = 1; j <= 2 + numSpindown; j++ ) {
    for ( k = 1; k <= j; k++ ) {
      printf( "  %+.4e", metric->data[k + j * ( j + 1 ) / 2] );
    }
    printf( "\n" );
  }

  /* Print determinants. */
  determinant = metric->data[5] * metric->data[2] - pow( metric->data[4], 2 );
  printf( "\nSky-determinant %e\n", determinant );
  if ( numSpindown == 1 ) {
    determinant = metric->data[2] * metric->data[5] * metric->data[9]
                  - metric->data[2] * metric->data[8] * metric->data[8]
                  + metric->data[4] * metric->data[8] * metric->data[7]
                  - metric->data[4] * metric->data[4] * metric->data[9]
                  + metric->data[7] * metric->data[4] * metric->data[8]
                  - metric->data[7] * metric->data[7] * metric->data[5];
    printf( "S&S determinant %e\n", determinant );
  }

  /* Here is the code that uses xmgrace with the -x option, */
  /* and outputs data to a file with the -t option. */
  if ( grace || nongrace ) {

    /* Take care of preliminaries. */
    if ( grace ) {
      pvc = popen( "xmgrace -pipe", "w" );
      if ( !pvc ) {
        printf( "%s line %d: %s\n", __FILE__, __LINE__,
                GENERALMETRICTESTC_MSGESYS );
        return GENERALMETRICTESTC_ESYS;
      }
      fprintf( pvc, "@xaxis label \"Right ascension (degrees)\"\n" );
      fprintf( pvc, "@yaxis label \"Declination (degrees)\"\n" );
    }
    if ( nongrace ) {
      fnongrace = fopen( "nongrace.data", "w" );
      if ( !fnongrace ) {
        printf( "%s line %d: %s\n", __FILE__, __LINE__,
                GENERALMETRICTESTC_MSGESYS );
        return GENERALMETRICTESTC_ESYS;
      }
    }

    /* Step around the sky: a grid in ra and dec. */
    j = 0;
    for ( dec = dec_max; dec >= dec_min; dec -= 10 ) {
      for ( ra = ra_min; ra <= ra_max; ra += 15 ) {
        REAL8 gaa, gad, gdd, angle, smaj, smin;

        /* Get the metric at this ra, dec. */
        in.position.longitude = ra * LAL_PI_180;
        in.position.latitude  = dec * LAL_PI_180;

        /* Evaluate metric: */
        LALPtoleMetric( &status, metric, &in );
        if ( status.statusCode ) {
          printf( "%s line %d: %s\n", __FILE__, __LINE__,
                  GENERALMETRICTESTC_MSGESUB );
          return GENERALMETRICTESTC_ESUB;
        }

        /*  Project metric: */
        LALProjectMetric( &status, metric, 0 );
        if ( status.statusCode ) {
          printf( "%s line %d: %s\n", __FILE__, __LINE__,
                  GENERALMETRICTESTC_MSGESUB );
          return GENERALMETRICTESTC_ESUB;
        }
        determinant = metric->data[5] * metric->data[2] - pow( metric->data[4], 2.0 );
        if ( determinant < 0.0 ) {
          printf( "%s line %d: %s\n", __FILE__, __LINE__,
                  GENERALMETRICTESTC_MSGEMET );
          return GENERALMETRICTESTC_EMET;
        }



        /* Rename \gamma_{\alpha\alpha}. */
        gaa = metric->data[2];
        /* Rename \gamma_{\alpha\delta}. */
        gad = metric->data[4];
        /* Rename \gamma_{\delta\delta}. */
        gdd = metric->data[5];
        /* Semiminor axis from larger eigenvalue of metric. */
        smin = gaa + gdd + sqrt( pow( gaa - gdd, 2 ) + pow( 2 * gad, 2 ) );
        smin = sqrt( 2 * mismatch / smin );
        /* Semiminor axis from smaller eigenvalue of metric. */
        smaj = gaa + gdd - sqrt( pow( gaa - gdd, 2 ) + pow( 2 * gad, 2 ) );
        /*printf("ra = %d, dec = %d, temp = %g\n", ra, dec, smaj);*/
        smaj = sqrt( 2 * mismatch / smaj );
        /* Angle of semimajor axis with "horizontal" (equator). */
        angle = atan2( gad, mismatch / smaj / smaj - gdd );
        if ( angle <= -LAL_PI_2 ) {
          angle += LAL_PI;
        }
        if ( angle > LAL_PI_2 ) {
          angle -= LAL_PI;
        }

        if ( grace ) {
          /* Print set header. */
          fprintf( pvc, "@s%d color (0,0,0)\n", j );
          fprintf( pvc, "@target G0.S%d\n@type xy\n", j++ );
          /* Print center of patch. */
          fprintf( pvc, "%16.8g %16.8g\n", ( float )ra, ( float )dec );
        }
        if ( nongrace )
          /* Print center of patch. */
        {
          fprintf( fnongrace, "%16.8g %16.8g\n", ( float )ra, ( float )dec );
        }
        /* Loop around patch ellipse. */
        for ( i = 0; i <= SPOKES; i++ ) {
          c_ellipse = LAL_TWOPI * i / SPOKES;
          r_ellipse = MAGNIFY * LAL_180_PI * smaj * smin /
                      sqrt( pow( smaj * sin( c_ellipse ), 2 ) + pow( smin * cos( c_ellipse ), 2 ) );
          if ( grace )
            fprintf( pvc, "%e %e\n", ra + r_ellipse * cos( angle - c_ellipse ),
                     dec + r_ellipse * sin( angle - c_ellipse ) );
          if ( nongrace )
            fprintf( fnongrace, "%e %e\n", ra + r_ellipse * cos( angle - c_ellipse ),
                     dec + r_ellipse * sin( angle - c_ellipse ) );

        } /* for (a...) */

      } /* for (ra...) */
    } /* for (dec...) */
    if ( grace ) {
      pclose( pvc );
    }
    if ( nongrace ) {
      fclose( fnongrace );
    }
  } /* if (grace || nongrace) */

  printf( "\nCleaning up and leaving...\n" );

  LALDDestroyVector( &status, &metric );
  if ( status.statusCode ) {
    printf( "%s line %d: %s\n", __FILE__, __LINE__,
            GENERALMETRICTESTC_MSGEMEM );
    return GENERALMETRICTESTC_EMEM;
  }
  if ( in.spindown ) {
    LALDestroyVector( &status, &( in.spindown ) );
  }

  if ( status.statusCode ) {
    printf( "%s line %d: %s\n", __FILE__, __LINE__,
            GENERALMETRICTESTC_MSGEMEM );
    return GENERALMETRICTESTC_EMEM;
  }
  LALCheckMemoryLeaks();
  return 0;
} /* main() */

/** \endcond */
