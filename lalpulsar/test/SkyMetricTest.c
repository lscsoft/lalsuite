/*
*  Copyright (C) 2007 Teviet Creighton
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
\author Creighton, T. D.
\file
\ingroup StackMetric_h

\heading{Program <tt>SkyMetricTest.c</tt>}
\latexonly\label{ss_SkyMetricTest_c}\endlatexonly

\brief Computes the sky-position metric for a coherent or semicoherent pulsar search.

\heading{Usage}
\code
SkyMetricTest [-o metricfile rangefile] [-p n dt t0 f0] [-l lat lon]
              [-r ra1 ra2 dec1 dec2] [-d debuglevel]
\endcode

\heading{Description}

This test program computes the parameter metric for the
two-dimensional sky-position pulsar stack search, over a grid of
points on the sky.  The following option flags are accepted:
<ul>
<li>[<tt>-o</tt>] Prints the metric grid and a grid specifying the
parameter ranges to the files \c metricfile and \c rangefile,
respectively, using the standard formatting routine
<tt>LALSWriteGrid()</tt> in the \c support package.  See below for a
description of these two grids.  If absent, the routines are
exercised, but no output is written.</li>
<li>[<tt>-p</tt>] Sets the search parameters: the number of stacks
\c n, the length of each stack \c dt (in seconds), and the
start time of the first stack \c t0 (in seconds of GPS time), and
the maximum source frequency \c f0 (in Hz).  If absent,
<tt>-t 1 86400 0 1000</tt> is assumed.</li>
<li>[<tt>-l</tt>] Sets the detector latitude to \c lat (in
degrees north from the equator) and longitude to \c lon (in
degrees east of the prime meridian).  If absent,
<tt>-l 52.247 9.822</tt> (GEO600) is assumed.</li>
<li>[<tt>-r</tt>] Sets the ``box'' on the sky to be covered: with
right ascension in the range [\c ra1,\c ra2] and declination
in the range [\c ra1,\c ra2], in degrees.  If absent,
<tt>-r 188.8594813 196.8594813 23.1282511 31.1282511</tt> is assumed (an
\f$8^\circ\f$ square centred on the Galactic core).</li>
<li>[<tt>-d</tt>] Sets the debug level to \c debuglevel; if
absent, <tt>-d 0</tt> is assumed.</li>
</ul>

\heading{Grid formats:} The metric grid is stored as a
\c REAL4Grid with physical dimension 2 and data dimesnion 3: for
each point in the two-dimensional sky grid \f$(\delta,\alpha)\f$,
representing declination and right ascension, it stores a 3-element
vector consisting of the three metric components \f$g_{\delta\delta}\f$,
\f$g_{\alpha\alpha}\f$, and \f$g_{\alpha\delta}\f$, in that order.  The
coordinates are given in degrees, and the metric components in
\f$\mathrm{degrees}^{-2}\f$.

The range grid stores the parameter space boundaries as a
\c REAL4Grid with physical dimension 1 and data dimension 2: for
each grid point in the right ascension, it stores a 2-element vector
consisting of the lower and upper values of declination for that right
ascension.  Both coordinates are given in degrees.


\heading{Algorithm}

After parsing the command-line inputs, this program defines a grid
covering the desired space, and determines the metric using the
routines LALStackMetric() and LALProjectMetric() with
the canonical time routine LALDTBaryPtolemaic().  These
routines return uncertainties along with their best-estimate
components, so we adopt a (perhaps overly) conservative approach: each
component is adjusted up or down by up to \f$1\sigma\f$ of estimated
uncertainty, in such a way as to maximize the metric determinant.
This will tend to give \e overcoverage of the parameter space.  If
the uncertainty in a component is larger than the component itself,
then a warning is generated.  If the final metric has non-positive
diagonal components or determinant, then the parameter space range is
shortened to exclude these points.  These degeneracies usually occur
at points near the poles for very short observation times, where there
is almost no modulation of the signal; usually in these cases an
unmodulated template will suffice to cover the excluded regions of the
parameter space.

The metric and range grids are written to the output file by the
routine <tt>LALSWriteGrid()</tt>, from which they can be read using
<tt>LALSReadGrid()</tt> in template-placement programs such as
\c TwoDMeshTest.c

\heading{Note:} This program, like the LALDTBaryPtolemaic()
phase model it relies on, uses right ascension and declination as its
parameter basis.  Since these are polar coordinates, they will go
singular at the poles regardless of the actual Doppler phase
modulation: specifically, the \f$g_{\alpha\alpha}\f$ and
\f$g_{\alpha\delta}\f$ metric components will go to zero, leading to
equimatch ellipses that are highly elongated in the \f$\alpha\f$
direction.

This gives rise to problems in template placement using the routines
in TwoDMesh.h, \e if columns are arranged along the
declination axis, since ellipse size and shape varies greatly along
the length of a column.  The problems largely go away if the columns
are arranged along the right ascension axis, where ellipse variations
are less severe.  This routine therefore orients the sky parameter
space so that declination runs along the first (\f$x\f$) axis and right
ascension along the second (\f$y\f$) axis.

\heading{Uses}
\code
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALMalloc()                     LALFree()
LALDCreateVector()              LALDDestroyVector()
LALU4CreateVector()             LALU4DestroyVector()
LALSCreateGrid()                LALSDestroyGrid()
LALSWriteGrid()                 snprintf()
LALStackMetric()                LALProjectMetric()
LALDTBaryPtolemaic()            LALGetEarthTimes()
\endcode

\heading{Notes}

*/

/** \name Error Codes */ /*@{*/
#define SKYMETRICTESTC_ENORM 0
#define SKYMETRICTESTC_ESUB  1
#define SKYMETRICTESTC_EARG  2
#define SKYMETRICTESTC_EVAL  3
#define SKYMETRICTESTC_EMEM  4
#define SKYMETRICTESTC_ERNG  5
#define SKYMETRICTESTC_EFILE 6

#define SKYMETRICTESTC_MSGENORM "Normal exit"
#define SKYMETRICTESTC_MSGESUB  "Subroutine failed"
#define SKYMETRICTESTC_MSGEARG  "Error parsing arguments"
#define SKYMETRICTESTC_MSGEVAL  "Input argument out of valid range"
#define SKYMETRICTESTC_MSGEMEM  "Memory allocation error"
#define SKYMETRICTESTC_MSGERNG  "Could not find area with valid metric"
#define SKYMETRICTESTC_MSGEFILE "Could not open file"
/*@}*/


/** \cond DONT_DOXYGEN */
#include <math.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Grid.h>
#include <lal/StreamOutput.h>
#include <lal/StackMetric.h>
#include <lal/PulsarTimes.h>

/* Default parameter settings. */
#define NSTACKS 1
#define STACKLENGTH 86400.0 /* arbitrary */
#define STARTTIME 0.0       /* arbitrary */
#define LATITUDE  52.247    /* GEO600 location */
#define LONGITUDE 9.822     /* GEO600 location */
#define FREQUENCY 1000.0    /* arbitrary */
#define RA_MIN 188.8594813  /* Galactic core search */
#define RA_MAX 196.8594813  /* Galactic core search */
#define DEC_MIN 23.1282511  /* Galactic core search */
#define DEC_MAX 31.1282511  /* Galactic core search */
#define MISMATCH 0.25       /* arbitrary but reasonably optimal */

/* Usage format string. */
#define USAGE "Usage: %s [-o metricfile rangefile] [-p n dt t0 f0]\n" \
"\t[-l lat lon] [-r ra1 ra2 dec1 dec2] [-d debuglevel]\n"

/* Input error checking: Some accepted parameter ranges. */
#define NMAX  10000 /* 1 <= number of stacks <= NMAX */
#define DTMAX  3e10 /* 1/f_0 < stack length <= DTMAX */
#define F0MAX  1e4  /* 0 < f_0 <= FOMAX */
/* Also: |latitude|  will be restricted to <= LAL_PI,
         |longitude| will be restricted to <= LAL_TWOPI
         mismatch    will be restricted to (0,1] */

/* Other internal constants. */
#define SPACING 5.0 /* maximum spacing between grid points (degrees) */
#define MSGLEN 1024 /* maximum warning message length */


/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, "$Id$", statement ? statement :   \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  XLALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
                 "        %s\n", *argv, __FILE__, __LINE__,          \
                 "$Id$", (statement) );                      \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 "$Id$", (statement) );                      \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( SKYMETRICTESTC_ESUB, SKYMETRICTESTC_MSGESUB,                \
         "Function call \"" #func "\" failed:" );                    \
  return SKYMETRICTESTC_ESUB;                                        \
}                                                                    \
while (0)

#define CHECKVAL( val, lower, upper )                                \
do                                                                   \
if ( ( (val) <= (lower) ) || ( (val) > (upper) ) )                   \
{                                                                    \
  ERROR( SKYMETRICTESTC_EVAL, SKYMETRICTESTC_MSGESUB,                \
         "Value of " #val " out of range:" );                        \
  XLALPrintError( #val " = %f, range = (%f,%f]\n", (REAL8)(val),      \
                 (REAL8)(lower), (REAL8)(upper) );                   \
  return SKYMETRICTESTC_EVAL;                                        \
}                                                                    \
while (0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main(int argc, char **argv)
{
  /* Variables for parsing command-line arguments. */
  INT4  arg;               /* argument list counter */
  CHAR *metricfile = NULL; /* metric output filename */
  CHAR *rangefile = NULL;  /* parameter range output filename */
  UINT4 n = NSTACKS;       /* number of stacks */
  REAL8 dt = STACKLENGTH;  /* length of each stack (s) */
  REAL8 t0 = STARTTIME;    /* start time of first stack (GPS s) */
  REAL8 f0 = FREQUENCY;    /* maximum frequency of search */
  REAL8 lat = LATITUDE;    /* latitude of detector (degrees N) */
  REAL8 lon = LONGITUDE;   /* longitude of detector (degrees E) */
  REAL8 ra[2], dec[2];     /* RA and dec ranges (degrees) */

  /* Other variables. */
  UINT4 i, j, k;                 /* indecies */
  UINT4 nRA, nDec;               /* sizes of metric grid in RA and dec */
  REAL8 dRA, dDec;               /* grid intervals in RA and dec */
  UINT4 warnings = 0;            /* points with large uncertainties */
  UINT4 errors = 0;              /* points with bad metrics */
  static LALStatus stat;         /* top-level status structure */
  static LIGOTimeGPS start;      /* GPS start time of first stack */
  REAL4Grid *metricGrid = NULL;  /* grid of metric values */
  REAL4 *gData;                  /* pointer to metricGrid->data->data */
  REAL4Grid *rangeGrid = NULL;   /* grid of range delimiters */
  REAL4 *rData;                  /* pointer to rangeGrid->data->data */
  REAL4 propVol = 0.0;           /* average proper volume element */
  UINT4 nVol = 0;                /* number of points used in average */
  UINT4Vector *dimLength = NULL; /* used to allocate metricGrid */
  REAL8Vector *lambda = NULL;    /* metric parameter space vector */
  REAL8Vector *metric = NULL;    /* metric components */
  REAL8 *mData;                  /* pointer to metric->data */
  INT4 *topRange, *botRange;     /* preliminary grid range limits */
  INT4 *topRange2, *botRange2;   /* adjusted grid range limits */
  static MetricParamStruc params; /* metric computation parameters */
  static PulsarTimesParamStruc baryParams; /* barycentring parameters */


  /* Some more initializations. */
  ra[0] = RA_MIN; dec[0] = DEC_MIN;
  ra[1] = RA_MAX; dec[1] = DEC_MAX;

  /******************************************************************
   * ARGUMENT PARSING                                               *
   ******************************************************************/

  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	metricfile = argv[arg++];
	rangefile = argv[arg++];
      } else {
	ERROR( SKYMETRICTESTC_EARG, SKYMETRICTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return SKYMETRICTESTC_EARG;
      }
    }
    /* Parse search parameter option. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 4 ) {
	arg++;
	n = atoi( argv[arg++] );
	dt = atof( argv[arg++] );
	t0 = atof( argv[arg++] );
	f0 = atof( argv[arg++] );
      } else {
	ERROR( SKYMETRICTESTC_EARG, SKYMETRICTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return SKYMETRICTESTC_EARG;
      }
    }
    /* Parse detector position option. */
    else if ( !strcmp( argv[arg], "-l" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	lat = atof( argv[arg++] );
	lon = atof( argv[arg++] );
      } else {
	ERROR( SKYMETRICTESTC_EARG, SKYMETRICTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return SKYMETRICTESTC_EARG;
      }
    }
    /* Parse parameter space range option. */
    else if ( !strcmp( argv[arg], "-r" ) ) {
      if ( argc > arg + 4 ) {
	arg++;
	ra[0] = atof( argv[arg++] );
	ra[1] = atof( argv[arg++] );
	dec[0] = atof( argv[arg++] );
	dec[1] = atof( argv[arg++] );
      } else {
	ERROR( SKYMETRICTESTC_EARG, SKYMETRICTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return SKYMETRICTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
      } else {
	ERROR( SKYMETRICTESTC_EARG, SKYMETRICTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return SKYMETRICTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( SKYMETRICTESTC_EARG, SKYMETRICTESTC_MSGEARG, 0 );
      XLALPrintError( USAGE, *argv );
      return SKYMETRICTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Do error trapping on input parameters. */
  if ( lalDebugLevel & LALERROR ) {
    CHECKVAL( n, 0, NMAX );
    CHECKVAL( f0, 0.0, F0MAX );
    CHECKVAL( dt, 1.0/f0, DTMAX );
    CHECKVAL( lat, -90.0, 90.0 );
    CHECKVAL( lon, -360.0, 360.0 );
  }

  /******************************************************************
   * METRIC COMPUTATION                                             *
   ******************************************************************/

  /* Set up start time. */
  start.gpsSeconds = (INT4)t0;
  start.gpsNanoSeconds = (INT4)( 1.0e9*( t0 - start.gpsSeconds ) );
  t0 = 0.0;

  /* Set up constant parameters for barycentre transformation. */
  baryParams.epoch = start;
  baryParams.latitude = LAL_PI_180*lat;
  baryParams.longitude = LAL_PI_180*lon;
  SUB( LALGetEarthTimes( &stat, &baryParams ), &stat );

  /* Set up constant parameters for metric calculation. */
  params.dtCanon = LALDTBaryPtolemaic;
  params.constants = &baryParams;
  params.start = t0;
  params.deltaT = dt;
  params.n = n;
  params.errors = 1;

  /* Set up the metric vector, metric grid, and parameter vector. */
  if ( params.errors )
    SUB( LALDCreateVector( &stat, &metric, 12 ), &stat );
  else
    SUB( LALDCreateVector( &stat, &metric, 6 ), &stat );
  mData = metric->data;
  SUB( LALU4CreateVector( &stat, &dimLength, 3 ), &stat );
  nRA = (UINT4)( ( ra[1] - ra[0] )/SPACING ) + 2;
  nDec = (UINT4)( ( dec[1] - dec[0] )/SPACING ) + 2;
  dimLength->data[0] = nDec;
  dimLength->data[1] = nRA;
  dimLength->data[2] = 3;
  SUB( LALSCreateGrid( &stat, &metricGrid, dimLength, 2 ), &stat );
  SUB( LALU4DestroyVector( &stat, &dimLength ), &stat );
  metricGrid->offset->data[0] = dec[0];
  metricGrid->offset->data[1] = ra[0];
  metricGrid->interval->data[0] = dDec = ( dec[1] - dec[0] )/( nDec - 1.0 );
  metricGrid->interval->data[1] = dRA = ( ra[1] - ra[0] )/( nRA - 1.0 );
  gData = metricGrid->data->data;
  SUB( LALDCreateVector( &stat, &lambda, 3 ), &stat );
  lambda->data[0] = f0;

  /* Set up preliminary range delimiters. */
  topRange = (INT4 *)LALMalloc( nRA*sizeof( INT4 ) );
  botRange = (INT4 *)LALMalloc( nRA*sizeof( INT4 ) );
  if ( !topRange || !botRange ) {
    ERROR( SKYMETRICTESTC_EMEM, SKYMETRICTESTC_MSGEMEM, 0 );
    return SKYMETRICTESTC_EMEM;
  }
  for ( i = 0; i < nRA; i++ ) {
    topRange[i] = 0;
    botRange[i] = nRA;
  }

  /* Fill the metric grid with the projected metric. */
  k = 0;
  fprintf( stdout, "Computing metric on %ix%i grid:\n", nDec, nRA );
  for ( i = 0; i < nDec; i++ ) {
    BOOLEAN bottom = 1, top = 0; /* whether we're at the bottom or top
                                    of the grid column */
    lambda->data[2] = LAL_PI_180*( dec[0] + i*dDec );
    for ( j = 0; j < nRA; j++ ) {
      if ( top )
	fprintf( stdout, "o" );
      else {
	REAL4 gxx, gyy, gxy; /* metric components */
	lambda->data[1] = LAL_PI_180*( ra[0] + j*dRA );
	SUB( LALStackMetric( &stat, metric, lambda, &params ), &stat );
	SUB( LALProjectMetric( &stat, metric, params.errors ), &stat );

	/* If uncertainties are given, generate conservative metric. */
	if ( params.errors ) {
	  INT2 sign = 1; /* sign of xy metric component */
	  gxx = mData[10] + mData[11];
	  gyy = mData[4] + mData[5];
	  gxy = mData[8];
	  if ( gxy < 0.0 )
	    sign = -1;
	  if ( fabs( gxy ) <= fabs( mData[9] ) )
	    gxy = 0.0;
	  else
	    gxy = sign*( fabs( gxy ) - fabs( mData[9] ) );
	  if ( gxx < 0.0 || gyy < 0.0 || gxx*gyy - gxy*gxy < 0.0 ) {
	    errors++;
	    fprintf( stdout, "o" );
	  } else if ( fabs( mData[4]*mData[11] ) +
		      fabs( mData[5]*mData[10] ) +
		      fabs( mData[8]*mData[9] )*2.0 >=
		      mData[4]*mData[10] - mData[8]*mData[8] ) {
	    warnings++;
	    fprintf( stdout, "*" );
	  } else
	    fprintf( stdout, "+" );
	}

	/* If uncertainties are not given, just use best estimate. */
	else {
	  gxx = mData[5];
	  gxy = mData[4];
	  gyy = mData[2];
	  if ( gxx < 0.0 || gyy < 0.0 || gxx*gyy - gxy*gxy < 0.0 ) {
	    errors++;
	    fprintf( stdout, "o" );
	  } else
	    fprintf( stdout, "+" );
	}

	/* See whether we've crossed a range boundary. */
	if ( gxx > 0.0 && gyy > 0.0 && gxx*gyy - gxy*gxy > 0.0 ) {
	  if ( bottom ) {
	    bottom = 0;
	    botRange[i] = j - 1;
	  }
	  propVol += LAL_PI_180*LAL_PI_180*sqrt( gxx*gyy - gxy*gxy );
	  nVol++;
	} else {
	  if ( !bottom ) {
	    top = 1;
	    topRange[i] = j;
	  }
	}
	*(gData++) = LAL_PI_180*LAL_PI_180*gxx;
	*(gData++) = LAL_PI_180*LAL_PI_180*gyy;
	*(gData++) = LAL_PI_180*LAL_PI_180*gxy;
      }
      fflush( stdout );
    }
    if ( !top )
      topRange[i] = j;
    fprintf( stdout, "\n" );
  }

  /* Free memory that's no longer needed. */
  SUB( LALDDestroyVector( &stat, &metric ), &stat );
  SUB( LALDDestroyVector( &stat, &lambda ), &stat );

  /* Warn if any points had bad values or large uncertainties. */
  if ( errors ) {
    CHAR msg[MSGLEN];
    snprintf( msg, MSGLEN, "%i of %i points had"
		 " non-positive-definite metric", warnings,
		 nRA*nDec );
    WARNING( msg );
  }
  if ( warnings ) {
    CHAR msg[MSGLEN];
    snprintf( msg, MSGLEN, "%i of %i points had metric component"
		 " errors larger than values", warnings, nRA*nDec );
    WARNING( msg );
  }

  /* Write the proper volume element. */
  fprintf( stdout, "Average proper volume element: %10.3e per square"
	   " degree\n", propVol /= nVol );

  /******************************************************************
   * RANGE COMPUTATION                                              *
   ******************************************************************/

  /* Each bad metric point affects the four squares adjoining it, so
     top and bottom ranges need to be adjusted down or up one point
     based on the adjoining columns. */
  topRange2 = (INT4 *)LALMalloc( nDec*sizeof(INT4) );
  botRange2 = (INT4 *)LALMalloc( nDec*sizeof(INT4) );
  if ( !topRange2 || !botRange2 ) {
    ERROR( SKYMETRICTESTC_EMEM, SKYMETRICTESTC_MSGEMEM, 0 );
    return SKYMETRICTESTC_EMEM;
  }
  topRange2[0] = topRange[1] < topRange[0] ?
    topRange[1] - 1 : topRange[0] - 1;
  botRange2[0] = botRange[1] > botRange[0] ?
    botRange[1] + 1 : botRange[0] + 1;
  topRange2[nDec-1] = topRange[nDec] < topRange[nDec-1] ?
    topRange[nDec] - 1 : topRange[nDec-1] - 1;
  botRange2[nDec-1] = botRange[nDec] > botRange[nDec-1] ?
    botRange[nDec] + 1 : botRange[nDec-1] + 1;
  for ( i = 1; i < nDec - 1; i++ ) {
    INT4 topMin = topRange[i-1], botMax = botRange[i-1];
    if ( topRange[i] < topMin )
      topMin = topRange[i];
    if ( topRange[i+1] < topMin )
      topMin = topRange[i+1];
    if ( botRange[i] > botMax )
      botMax = botRange[i];
    if ( botRange[i+1] > botMax )
      botMax = botRange[i+1];
    topRange2[i] = topMin - 1;
    botRange2[i] = botMax + 1;
  }
  LALFree( topRange );
  LALFree( botRange );

  /* Determine the contiguous domain with open range intervals. */
  if ( topRange2[0] > botRange2[0] && topRange2[1] > botRange2[1] )
    i = 0;
  else {
    i = 1;
    while ( i < nDec - 1 && topRange2[i] <= botRange2[i] &&
	    topRange2[i+1] <= botRange2[i+1] )
      i++;
  }
  j = i;
  while ( j < nDec - 1 && topRange2[j+1] > botRange2[j+1] )
    j++;
  if ( j == i ) {
    ERROR( SKYMETRICTESTC_ERNG, SKYMETRICTESTC_MSGERNG, 0 );
    return SKYMETRICTESTC_ERNG;
  }

  /* Set up range grid. */
  j = j - i + 1;
  SUB( LALU4CreateVector( &stat, &dimLength, 2 ), &stat );
  dimLength->data[0] = j;
  dimLength->data[1] = 2;
  SUB( LALSCreateGrid( &stat, &rangeGrid, dimLength, 1 ), &stat );
  SUB( LALU4DestroyVector( &stat, &dimLength ), &stat );
  rangeGrid->offset->data[0] = dec[0] + i*dDec;
  rangeGrid->interval->data[0] = dDec;
  rData = rangeGrid->data->data;
  nVol = 0;
  for ( k = 0; k < j; k++ ) {
    *(rData++) = ra[0] + botRange2[k]*dRA;
    *(rData++) = ra[0] + topRange2[k]*dRA;
    nVol += topRange2[k] - botRange2[k] - 1;
  }
  nVol -= topRange2[j-1] - botRange2[j-1] - 1;
  LALFree( topRange2 );
  LALFree( botRange2 );

  /* Write the total proper volume of the covered range. */
  fprintf( stdout, "Total proper volume:           %10.3e\n",
	   propVol*nVol*dRA*dDec );

  /******************************************************************
   * FINISH                                                         *
   ******************************************************************/

  /* Write the output files, if requested. */
  if ( metricfile ) {
    FILE *fp = fopen( metricfile, "w" ); /* output file pointer */
    if ( !fp ) {
      ERROR( SKYMETRICTESTC_EFILE, "- " SKYMETRICTESTC_MSGEFILE,
	     metricfile );
      return SKYMETRICTESTC_EFILE;
    }
    SUB( LALSWriteGrid( &stat, fp, metricGrid ), &stat );
    fclose( fp );
  }
  if ( rangefile ) {
    FILE *fp = fopen( rangefile, "w" ); /* output file pointer */
    if ( !fp ) {
      ERROR( SKYMETRICTESTC_EFILE, "- " SKYMETRICTESTC_MSGEFILE,
	     rangefile );
      return SKYMETRICTESTC_EFILE;
    }
    SUB( LALSWriteGrid( &stat, fp, rangeGrid ), &stat );
    fclose( fp );
  }

  /* Clean up and exit. */
  SUB( LALSDestroyGrid( &stat, &metricGrid ), &stat );
  SUB( LALSDestroyGrid( &stat, &rangeGrid ), &stat );
  LALCheckMemoryLeaks();
  INFO( SKYMETRICTESTC_MSGENORM );
  return SKYMETRICTESTC_ENORM;
}

/** \endcond */
