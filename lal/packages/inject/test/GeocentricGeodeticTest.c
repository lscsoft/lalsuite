/*
*  Copyright (C) 2008 Teviet Creighton
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

/*********************** <lalVerbatim file="GeocentricGeodeticTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{GeocentricGeodeticTest.c}}
\label{ss:GeocentricGeodeticTest.c}

Tests geocentric to geodetic conversion.

\subsubsection*{Usage}
\begin{verbatim}
GeocentricGeodeticTest [-x xmin xmax nx] [-y ymin ymax ny]
                       [-z zmin zmax nz] [-v] [-d debuglevel]
\end{verbatim}

\subsubsection*{Description}

This program converts a point or set of points from geocentric to
geodetic coordinates and back, using the routines
\verb@LALGeocentricToGeodetic()@ and \verb@LALGeodeticToGeocentric()@.
The reconverted position is compared with the original, and the
maximum difference (in metres) reported to \verb@stderr@.  The
following option flags are accepted:
\begin{itemize}
\item[\texttt{-x}] Specifies a range and number of points to span in
the geocentric Cartesian $x$-coordinate.  Range limits are given in
Earth radii, and the number of points must be at least 1.  If not
specified, a single value of 0 is assumed.
\item[\texttt{-y}] As \verb@-x@, above, but for the $y$-coordinate.
\item[\texttt{-z}] As \verb@-x@, above, but for the $z$-coordinate.
\item[\texttt{-v}] Specifies verbosity level.  The default is level 0,
printing only the maximum difference to \verb@stderr@.  Level 1 in
addition prints to \verb@stdout@ the difference measured at each
point.  Level 2 prints the elevation of the point, followed by the
difference, both in metres (this is to facilitate shaded diagrams such
as Fig.~\ref{fig:geodetictest}).  Level 3 prints the $x$, $y$, and $z$
coordinates of each point followed by the difference measured at that
point (all in metres).  Level 4 prints the geocentric Cartesian
coordinates, above, folowed by the geodetic elevation, latitude and
longitude, followed by the difference (in metres or degrees as
appropriate).
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@.  If not
specified, level 0 is assumed.
\end{itemize}
If neither \verb@-x@, \verb@-y@, or \verb@-z@ were specified, then the
program will test a single randomly generated position between 0.5 and
2 Earth radii, and return an error if the conversion produces a
difference greater than a micron.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define GEOCENTRICGEODETICTESTC_ENORM 0
#define GEOCENTRICGEODETICTESTC_ESUB  1
#define GEOCENTRICGEODETICTESTC_EARG  2
#define GEOCENTRICGEODETICTESTC_EMEM  3
#define GEOCENTRICGEODETICTESTC_ETEST 4

#define GEOCENTRICGEODETICTESTC_MSGENORM "Normal exit"
#define GEOCENTRICGEODETICTESTC_MSGESUB  "Subroutine failed"
#define GEOCENTRICGEODETICTESTC_MSGEARG  "Error parsing arguments"
#define GEOCENTRICGEODETICTESTC_MSGEMEM  "Out of memory"
#define GEOCENTRICGEODETICTESTC_MSGETEST "Test case failed"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

\begin{wrapfigure}{r}{0.6\textwidth}
\vspace{-2ex}
\begin{center}
\resizebox{0.55\textwidth}{!}{\includegraphics{inject_geodetictest}}
\\ \parbox{0.55\textwidth}{\caption{\label{fig:geodetictest} Precision
of geocentric-geodetic conversion algorithm.  Shaded ellipse is the
reference ellipsoid.  The wedges of (comparatively) lower precision
occur near the equator, where the series expansion in $B$ is
required.}}
\end{center}
\vspace{-4ex}
\end{wrapfigure}
See \verb@TerrestrialCoordinates.c()@ for documentation about the
geocentric/geodetic conversion algorithm.  Since
\verb@LALGeodeticToGeocentric()@ is fairly straightforward and
numerically robust, this program is basically a test of the
\verb@LALGeocentricToGeodetic()@ algorithm.

Running with verbosity level 2 gives error data that can be used to
generate figures such as Fig.~\ref{fig:geodetictest}.  First we note
that for points near the surface of the Earth, the position error is
never greater than a micron, or about one part in $10^{12}$ --- a
reasonable expectation for double-precision arithmetic.  The largest
errors occur near the equator, where the expression for $v$
experiences loss of precision.  In this limit we replace the ``exact''
expression for $v$ with a power series, as described in
\verb@TerrestrialCoordinates.c()@.  This restores precision near the
equatorial plane.  The point of transition has been tuned by hand, so
that the errors due to series truncation just within the wedge are
about the same as the numerical loss of precision just outside of it.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALGeocentricToGeodetic()       LALGeodeticToGeocentric()
LALCreateRandomParams()         LALDestroyRandomParams()
LALUniformDeviate()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GeocentricGeodeticTestCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>

NRCSID( GEOCENTRICGEODETICTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;

/* Usage format string. */
#define USAGE "Usage: %s [-x xmin xmax nx] [-y ymin ymax ny]\n" \
"\t[-z zmin zmax nz] [-o outfle] [-d debuglevel]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, GEOCENTRICGEODETICTESTC, statement ?      \
                 statement : "", (msg) );                            \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 GEOCENTRICGEODETICTESTC, (statement) );             \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 GEOCENTRICGEODETICTESTC, (statement) );             \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( GEOCENTRICGEODETICTESTC_ESUB,                               \
	 GEOCENTRICGEODETICTESTC_MSGESUB,                            \
         "Function call \"" #func "\" failed:" );                    \
  return GEOCENTRICGEODETICTESTC_ESUB;                               \
}                                                                    \
while (0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main( int argc, char **argv )
{
  int arg;                      /* command-line argument counter */
  BOOLEAN xyz = 0;              /* whether -x, -y, or -z options were given */
  INT4 verbosity = 0;           /* verbosity level */
  REAL8 x_1 = 0.0, x_2 = 0.0, dx; /* range and increment in x */
  REAL8 y_1 = 0.0, y_2 = 0.0, dy; /* range and increment in y */
  REAL8 z_1 = 0.0, z_2 = 0.0, dz; /* range and increment in z */
  INT4 nx = 1, ny = 1, nz = 1;  /* number of steps in each direction */
  INT4 i, j, k;                 /* index in each direction */
  REAL8 x, y, z, ddx, ddy, ddz; /* position and error in each direction */
  REAL8 ddr, ddmax = 0.0;       /* overall and maximum position error */
  static LALStatus stat;        /* status structure */
  EarthPosition earth;          /* terrestrial coordinates */

  /*******************************************************************
   * PARSE ARGUMENTS (arg stores the current position)               *
   *******************************************************************/

  arg = 1;
  while ( arg < argc ) {

    /* Parse range options. */
    if ( !strcmp( argv[arg], "-x" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	xyz = 1;
	x_1 = atof( argv[arg++] );
	x_2 = atof( argv[arg++] );
	nx = atoi( argv[arg++] );
      } else {
	ERROR( GEOCENTRICGEODETICTESTC_EARG,
	       GEOCENTRICGEODETICTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GEOCENTRICGEODETICTESTC_EARG;
      }
    } else if ( !strcmp( argv[arg], "-y" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	xyz = 1;
	y_1 = atof( argv[arg++] );
	y_2 = atof( argv[arg++] );
	ny = atoi( argv[arg++] );
      } else {
	ERROR( GEOCENTRICGEODETICTESTC_EARG,
	       GEOCENTRICGEODETICTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GEOCENTRICGEODETICTESTC_EARG;
      }
    } else if ( !strcmp( argv[arg], "-z" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	xyz = 1;
	z_1 = atof( argv[arg++] );
	z_2 = atof( argv[arg++] );
	nz = atoi( argv[arg++] );
      } else {
	ERROR( GEOCENTRICGEODETICTESTC_EARG,
	       GEOCENTRICGEODETICTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GEOCENTRICGEODETICTESTC_EARG;
      }
    }

    /* Parse verbosity option. */
    else if ( !strcmp( argv[arg], "-v" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	verbosity = atoi( argv[arg++] );
      } else {
	ERROR( GEOCENTRICGEODETICTESTC_EARG,
	       GEOCENTRICGEODETICTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GEOCENTRICGEODETICTESTC_EARG;
      }
    }

    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      } else {
	ERROR( GEOCENTRICGEODETICTESTC_EARG,
	       GEOCENTRICGEODETICTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return GEOCENTRICGEODETICTESTC_EARG;
      }
    }

    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( GEOCENTRICGEODETICTESTC_EARG,
	     GEOCENTRICGEODETICTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return GEOCENTRICGEODETICTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /*******************************************************************
   * SET PARAMETERS                                                  *
   *******************************************************************/

  /* If none of -x, -y, -z were given, generate a random position. */
  if ( !xyz ) {
    REAL4 lon, coslat, sinlat, rad; /* polar coordinates */
    RandomParams *rparams = NULL;   /* pseudorandom sequence parameters */
    SUB( LALCreateRandomParams( &stat, &rparams, 0 ), &stat );
    SUB( LALUniformDeviate( &stat, &lon, rparams ), &stat );
    lon *= LAL_TWOPI;
    SUB( LALUniformDeviate( &stat, &sinlat, rparams ), &stat );
    coslat = sqrt( 1.0 - sinlat*sinlat );
    SUB( LALUniformDeviate( &stat, &rad, rparams ), &stat );
    rad = 1.5*rad + 0.5;
    x_1 = x_2 = rad*coslat*cos( lon );
    y_1 = y_2 = rad*coslat*sin( lon );
    z_1 = z_2 = rad*sinlat;
    SUB( LALDestroyRandomParams( &stat, &rparams ), &stat );
  }

  /* Compute stepsizes. */
  dx = dy = dz = 0.0;
  if ( nx > 1 )
    dx = ( x_2 - x_1 )/( nx - 1.0 );
  if ( ny > 1 )
    dy = ( y_2 - y_1 )/( ny - 1.0 );
  if ( nz > 1 )
    dz = ( z_2 - z_1 )/( nz - 1.0 );

  /*******************************************************************
   * PERFORM TEST                                                    *
   *******************************************************************/

  /* Loop over each direction. */
  for ( i = 0; i < nx; i++ ) {
    x = LAL_REARTH_SI*( x_1 + i*dx );
    for ( j = 0; j < ny; j++ ) {
      y = LAL_REARTH_SI*( y_1 + j*dy );
      for ( k = 0; k < nz; k++ ) {
	z = LAL_REARTH_SI*( z_1 + k*dz );

	/* Do transformation. */
	earth.x = x;
	earth.y = y;
	earth.z = z;
	SUB( LALGeocentricToGeodetic( &stat, &earth ), &stat );
 	SUB( LALGeodeticToGeocentric( &stat, &earth ), &stat );

	/* Compute difference. */
	ddx = x - earth.x;
	ddy = y - earth.y;
	ddz = z - earth.z;
	ddr = sqrt( ddx*ddx + ddy*ddy + ddz*ddz );
	if ( ddr > ddmax )
	  ddmax = ddr;
	if ( verbosity == 1 )
	  fprintf( stdout, "%+23.16e\n", ddr );
	else if ( verbosity == 2 )
	  fprintf( stdout, "%+23.16e %+23.16e\n", earth.elevation, ddr );
	else if ( verbosity == 3 )
	  fprintf( stdout, "%+23.16e %+23.16e %+23.16e %+23.16e\n",
		   x, y, z, ddr );
	else if ( verbosity == 4 )
	  fprintf( stdout, "%+23.16e %+23.16e %+23.16e"
		   " %+23.16e %+23.16e %+23.16e %+23.16e\n", x, y, z,
		   earth.elevation, LAL_180_PI*earth.geodetic.latitude,
		   LAL_180_PI*earth.geodetic.longitude, ddr );
      }
    }
  }

  /* Print maximum. */
  fprintf( stdout, "Maximum error: %.16em\n", ddmax );
  if ( !xyz && ddmax > 1.0e-6 ) {
    ERROR( GEOCENTRICGEODETICTESTC_ETEST,
	   GEOCENTRICGEODETICTESTC_MSGETEST, 0 );
    return GEOCENTRICGEODETICTESTC_ETEST;
  }
  LALCheckMemoryLeaks();
  INFO( GEOCENTRICGEODETICTESTC_MSGENORM );
  return GEOCENTRICGEODETICTESTC_ENORM;
}
