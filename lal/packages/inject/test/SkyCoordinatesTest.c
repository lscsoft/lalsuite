/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix, Teviet Creighton
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

/*************************** <lalVerbatim file="SkyCoordinatesTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{SkyCoordinatesTest.c}}
\label{ss:SkyCoordinatesTest.c}

Transforms coordinates among various systems.

\subsubsection*{Usage}
\begin{verbatim}
SkyCoordinatesTest [-i system lat lon] [-o system] [-z lat lon]
                   [-a altitude] [-c lat lon rad] [-t sec nsec] [-d debuglevel]
\end{verbatim}

\subsubsection*{Description}

This program converts between any two celestial coordinate systems, or
between geocentric and geodetic terrestrial coordinates, using the
routines in \verb@SkyCoordinates.h@.  The following option flags are
accepted:
\begin{itemize}
\item[\texttt{-i}] Sets the input coordinate system and coordinate
values for a celestial coordinate trasformation: \verb@system@ may be
one of \verb@horizon@, \verb@geographic@, \verb@equatorial@,
\verb@ecliptic@, or \verb@galactic@; \verb@lat@ and \verb@lon@ are the
latitude and longitude coordinates in that system (in degrees).  If
the \verb@-i@ option is not given, then no celestial coordinate
transformation will be performed (although a terrestrial coordinate
transformation may still occur; see below).
\item[\texttt{-o}] Sets the output coordinate system for a cellestial
coodinate transformation: \verb@system@ may be any of the above.  If the
\verb@-o@ option is not given, then no celestial coordinate
transformation will be performed (although a terrestrial coordinate
transformation may still occur; see below).
\item[\texttt{-z}] Sets the \emph{geodetic} latitude and longitude of
the observer to \verb@lat@ and \verb@lon@, respectively (in degrees).
Either this or the \verb@-c@ option (below) is required for a
celestial coordinate transformation involving the horizon system.
\item[\texttt{-a}] Sets the elevation of the observer above the
Earth's reference ellipsoid to \verb@altitude@ (in metres).  If given
along with the \verb@-z@ option, above, the program will compute and
print out the geocentric coordinates of the observer as well.
\item[\texttt{-c}] Sets the \emph{geocentric} latitude and longitude
of the observer to \verb@lat@ and \verb@lon@, respectively (in
degrees), and the distance from the geocentre to \verb@rad@ (in
metres).  The program will convert and print out the geodetic
coordinates of the observer.  Either this or the \verb@-z@ option
(below) is required for a celestial coordinate transformation
involving the horizon system; if both are given, this option is
ignored.
\item[\texttt{-t}] Sets the GPS time of the conversion to \verb@sec@
seconds plus \verb@nsec@ nanoseconds.  The time will be printed in
various other formats.  This option is required for any transformation
between Earth-fixed and sky-fixed coordinate systems.
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@.  If not
specified, level 0 is assumed.
\end{itemize}
If no option flags are specified at all, then the routine will
randomly generate a sky position in Galactic coordinates, convert it
to ecliptic coordinates and back again, and return an error if the
result disagrees by more than a milliradian.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define SKYCOORDINATESTESTC_ENORM 0
#define SKYCOORDINATESTESTC_ESUB  1
#define SKYCOORDINATESTESTC_EARG  2
#define SKYCOORDINATESTESTC_EMEM  3
#define SKYCOORDINATESTESTC_ETEST 4

#define SKYCOORDINATESTESTC_MSGENORM "Normal exit"
#define SKYCOORDINATESTESTC_MSGESUB  "Subroutine failed"
#define SKYCOORDINATESTESTC_MSGEARG  "Error parsing arguments"
#define SKYCOORDINATESTESTC_MSGEMEM  "Out of memory"
#define SKYCOORDINATESTESTC_MSGETEST "Test case failed"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALMalloc()                     LALFree()
LALGeocentricToGeodetic()       LALGeodeticToGeocentric()
LALGPStoGMST1()
LALCHARCreateVector()           LALCHARDestroyVector()
LALCreateRandomParams()         LALDestroyRandomParams()
LALUniformDeviate()             LALConvertSkyCoordinates()
LALNormalizeSkyPosition()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{SkyCoordinatesTestCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>

NRCSID( SKYCOORDINATESTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;

/* Usage format string. */
#define USAGE "Usage: %s [-i system lat lon] [-o system] [-z lat lon]\n" \
"\t[-a altitude] [-c lat lon rad] [-t sec nsec] [-d debuglevel]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, SKYCOORDINATESTESTC, statement ?          \
                 statement : "", (msg) );                            \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 SKYCOORDINATESTESTC, (statement) );                 \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 SKYCOORDINATESTESTC, (statement) );                 \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( SKYCOORDINATESTESTC_ESUB, SKYCOORDINATESTESTC_MSGESUB,      \
         "Function call \"" #func "\" failed:" );                    \
  return SKYCOORDINATESTESTC_ESUB;                                   \
}                                                                    \
while (0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif


int
main( int argc, char **argv )
{
  int arg;                 /* command-line argument counter */
  BOOLEAN in = 0, out = 0; /* whether -i or -o options were given */
  BOOLEAN gd = 0, alt = 0; /* whether -z or -a options were given */
  BOOLEAN gc = 0, t = 0;   /* whether -c or -t options were given */
  static LALStatus stat;   /* status structure */
  LIGOTimeGPS gpsTime;     /* time of transformation */
  SkyPosition sky;         /* celestial coordinates */
  EarthPosition earth;     /* terrestrial coordinates */
  ConvertSkyParams params; /* additional parameters for conversion */
  memset( &params, 0, sizeof(ConvertSkyParams ) );

  /*******************************************************************
   * PARSE ARGUMENTS (arg stores the current position)               *
   *******************************************************************/

  arg = 1;
  while ( arg < argc ) {

    /* Parse input option. */
    if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	in = 1;
        if ( !strncmp( argv[arg], "h", 1 ) )
	  sky.system = COORDINATESYSTEM_HORIZON;
        else if ( !strncmp( argv[arg], "ge", 2 ) )
	  sky.system = COORDINATESYSTEM_GEOGRAPHIC;
        else if ( !strncmp( argv[arg], "eq", 2 ) )
	  sky.system = COORDINATESYSTEM_EQUATORIAL;
        else if ( !strncmp( argv[arg], "ec", 2 ) )
	  sky.system = COORDINATESYSTEM_ECLIPTIC;
        else if ( !strncmp( argv[arg], "ga", 2 ) )
	  sky.system = COORDINATESYSTEM_GALACTIC;
	else {
	  ERROR( SKYCOORDINATESTESTC_EARG, SKYCOORDINATESTESTC_MSGEARG, 0 );
	  LALPrintError( USAGE, *argv );
	  return SKYCOORDINATESTESTC_EARG;
	}
	arg++;
	sky.latitude = LAL_PI_180*atof( argv[arg++] );
	sky.longitude = LAL_PI_180*atof( argv[arg++] );
      } else {
	ERROR( SKYCOORDINATESTESTC_EARG, SKYCOORDINATESTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SKYCOORDINATESTESTC_EARG;
      }
    }

    /* Parse output option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	out = 1;
        if ( !strncmp( argv[arg], "h", 1 ) )
	  params.system = COORDINATESYSTEM_HORIZON;
        else if ( !strncmp( argv[arg], "ge", 2 ) )
	  params.system = COORDINATESYSTEM_GEOGRAPHIC;
        else if ( !strncmp( argv[arg], "eq", 2 ) )
	  params.system = COORDINATESYSTEM_EQUATORIAL;
        else if ( !strncmp( argv[arg], "ec", 2 ) )
	  params.system = COORDINATESYSTEM_ECLIPTIC;
        else if ( !strncmp( argv[arg], "ga", 2 ) )
	  params.system = COORDINATESYSTEM_GALACTIC;
	else {
	  ERROR( SKYCOORDINATESTESTC_EARG, SKYCOORDINATESTESTC_MSGEARG, 0 );
	  LALPrintError( USAGE, *argv );
	  return SKYCOORDINATESTESTC_EARG;
	}
	arg++;
      } else {
	ERROR( SKYCOORDINATESTESTC_EARG, SKYCOORDINATESTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SKYCOORDINATESTESTC_EARG;
      }
    }

    /* Parse zenith position option. */
    else if ( !strcmp( argv[arg], "-z" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	gd = 1;
	earth.geodetic.system = COORDINATESYSTEM_GEOGRAPHIC;
	earth.geodetic.latitude = LAL_PI_180*atof( argv[arg++] );
	earth.geodetic.longitude = LAL_PI_180*atof( argv[arg++] );
      } else {
	ERROR( SKYCOORDINATESTESTC_EARG, SKYCOORDINATESTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SKYCOORDINATESTESTC_EARG;
      }
    }

    /* Parse observer altitude option. */
    else if ( !strcmp( argv[arg], "-a" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	alt = 1;
	earth.elevation = atof( argv[arg++] );
      } else {
	ERROR( SKYCOORDINATESTESTC_EARG, SKYCOORDINATESTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SKYCOORDINATESTESTC_EARG;
      }
    }

    /* Parse geocentric coordinates option. */
    else if ( !strcmp( argv[arg], "-c" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	gc = 1;
	earth.geocentric.system = COORDINATESYSTEM_GEOGRAPHIC;
	earth.geocentric.latitude = LAL_PI_180*atof( argv[arg++] );
	earth.geocentric.longitude = LAL_PI_180*atof( argv[arg++] );
	earth.radius = atof( argv[arg++] );
      } else {
	ERROR( SKYCOORDINATESTESTC_EARG, SKYCOORDINATESTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SKYCOORDINATESTESTC_EARG;
      }
    }

    /* Parse GPS time option. */
    else if ( !strcmp( argv[arg], "-t" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	t = 1;
	gpsTime.gpsSeconds = atof( argv[arg++] );
	gpsTime.gpsNanoSeconds = atof( argv[arg++] );
      } else {
	ERROR( SKYCOORDINATESTESTC_EARG, SKYCOORDINATESTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SKYCOORDINATESTESTC_EARG;
      }
    }

    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      }else{
	ERROR( SKYCOORDINATESTESTC_EARG, SKYCOORDINATESTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SKYCOORDINATESTESTC_EARG;
      }
    }

    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( SKYCOORDINATESTESTC_EARG, SKYCOORDINATESTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return SKYCOORDINATESTESTC_EARG;
    }
  } /* End of argument parsing loop. */


  /*******************************************************************
   * TEST CASE                                                       *
   *******************************************************************/

  /* If no arguments were given, convert a random position and check
     the result. */
  if ( !in && !out && !gd && !gc && !alt && !t ) {
    REAL4 deviate;    /* random number between 0 and 1 */
    REAL8 diff;       /* difference in converted longitude */
    SkyPosition sky2; /* converted celestial coordinates */
    RandomParams *rparams = NULL; /* pseudorandom sequence parameters */


    printf ("Testing LALConvertSkyCoordinates()\n");
    /* Set up random Galactic position. */
    SUB( LALCreateRandomParams( &stat, &rparams, 0 ), &stat );
    SUB( LALUniformDeviate( &stat, &deviate, rparams ), &stat );
    sky.latitude = LAL_PI*deviate - LAL_PI_2;
    SUB( LALUniformDeviate( &stat, &deviate, rparams ), &stat );
    sky.longitude = LAL_TWOPI*deviate;
    SUB( LALDestroyRandomParams( &stat, &rparams ), &stat );
    sky.system = COORDINATESYSTEM_GALACTIC;

    /* Convert to Ecliptic and back. */
    fprintf( stdout, "Galactic   (%7.3f,%7.3f)\n",
	     LAL_180_PI*sky.latitude, LAL_180_PI*sky.longitude );
    params.system = COORDINATESYSTEM_EQUATORIAL;
    SUB( LALConvertSkyCoordinates( &stat, &sky2, &sky, &params ),
	 &stat );
    fprintf( stdout, "Equatorial (%7.3f,%7.3f)\n",
	     LAL_180_PI*sky2.latitude, LAL_180_PI*sky2.longitude );
    params.system = COORDINATESYSTEM_ECLIPTIC;
    SUB( LALConvertSkyCoordinates( &stat, &sky2, &sky2, &params ),
	 &stat );
    fprintf( stdout, "Ecliptic   (%7.3f,%7.3f)\n",
	     LAL_180_PI*sky2.latitude, LAL_180_PI*sky2.longitude );
    params.system = COORDINATESYSTEM_EQUATORIAL;
    SUB( LALConvertSkyCoordinates( &stat, &sky2, &sky2, &params ),
	 &stat );
    fprintf( stdout, "Equatorial (%7.3f,%7.3f)\n",
	     LAL_180_PI*sky2.latitude, LAL_180_PI*sky2.longitude );
    params.system = COORDINATESYSTEM_GALACTIC;
    SUB( LALConvertSkyCoordinates( &stat, &sky2, &sky2, &params ),
	 &stat );
    fprintf( stdout, "Galactic   (%7.3f,%7.3f)\n",
	     LAL_180_PI*sky2.latitude, LAL_180_PI*sky2.longitude );

    /* Make sure conversion is consistent. */
    diff = sky2.longitude - sky.longitude;
    while ( diff < -LAL_PI )
      diff += LAL_TWOPI;
    while ( diff > LAL_PI )
      diff -= LAL_TWOPI;
    if ( fabs( sky2.latitude - sky.latitude ) > 0.001 ||
	 fabs( diff ) > 0.001 ) {
      ERROR( SKYCOORDINATESTESTC_ETEST, SKYCOORDINATESTESTC_MSGETEST, 0 );
      return SKYCOORDINATESTESTC_ETEST;
    }


    /***********************************************************************
     * Test LALNormalizeSkyPosition()
     * we completely ignore the coordinate-system here, they are all treated
     * the same (provided we can assume they are spherical and their coordinates
     * lie within the range [0,2pi)x[-pi/2,pi/2]
     ************************************************************************/
    {
      SkyPosition testIn, testOut;
      SkyPosition correct = { 3.4, 0.9, 0 };
      printf ("Testing LALNormalizeSkyPosition()\n");
      /* now add some cycles of 2pi*/
      testIn.longitude = correct.longitude + 4 * LAL_TWOPI;
      testIn.latitude  = correct.latitude - 7 * LAL_TWOPI;
      SUB (LALNormalizeSkyPosition (&stat, &testOut, &testIn), &stat);
      if ( (fabs(testOut.longitude - correct.longitude) > 1e-14)
	   || (fabs(testOut.latitude-correct.latitude)>1e-14))
	{
	  printf ( "1.) LALNormalizeSkyPosition failed: got (%f,%f) instead of (%f,%f)\n",
		   testOut.longitude, testOut.latitude, correct.longitude, correct.latitude);
	  ERROR( SKYCOORDINATESTESTC_ETEST, SKYCOORDINATESTESTC_MSGETEST, 0 );
	  return SKYCOORDINATESTESTC_ETEST;
	}
      /* try going over the pole */
      testIn.latitude = LAL_PI - testIn.latitude;
      testIn.longitude += LAL_PI + 2 * LAL_TWOPI;
      SUB (LALNormalizeSkyPosition (&stat, &testOut, &testIn), &stat);
      if ( (fabs(testOut.longitude - correct.longitude) > 1e-14)
	   || (fabs(testOut.latitude-correct.latitude)>1e-14))
	{
	  printf ( "2.) LALNormalizeSkyPosition failed: got (%f,%f) instead of (%f,%f)\n",
		   testOut.longitude, testOut.latitude, correct.longitude, correct.latitude);
	  ERROR( SKYCOORDINATESTESTC_ETEST, SKYCOORDINATESTESTC_MSGETEST, 0 );
	  return SKYCOORDINATESTESTC_ETEST;
	}


    } /* testing LALNormalizeSkyPosition() */
    /***********************************************************************/


    /* Everything's fine, and nothing should have been allocated. */
    LALCheckMemoryLeaks();
    INFO( SKYCOORDINATESTESTC_MSGENORM );
    return SKYCOORDINATESTESTC_ENORM;
  }


  /*******************************************************************
   * TERRESTRIAL COORDINATES                                         *
   *******************************************************************/

  /* Convert geocentric to geodetic, if required. */
  fprintf( stdout, "\n" );
  if ( gc && !gd ) {
    fprintf( stdout, "TERRESTRIAL COORDINATES\n"
	     "Geocentric: latitude = %6.2f deg, longitude = %6.2f deg,"
	     " radius = %.0fm\n", LAL_180_PI*earth.geocentric.latitude,
	     LAL_180_PI*earth.geocentric.longitude, 1.0*earth.radius );
    earth.x = earth.radius*cos( earth.geocentric.latitude )
      *cos( earth.geocentric.longitude );
    earth.y = earth.radius*cos( earth.geocentric.latitude )
      *sin( earth.geocentric.longitude );
    earth.z = earth.radius*sin( earth.geocentric.latitude );
    fprintf( stdout,
	     "            x = %.0fm, y = %.0fm, z = %.0fm\n",
	     1.0*earth.x, 1.0*earth.y, 1.0*earth.z );
    SUB( LALGeocentricToGeodetic( &stat, &earth ), &stat );
    fprintf( stdout,
	     "Geodetic:   latitude = %6.2f deg, longitude = %6.2f"
	     " deg, elevation = %.0fm\n",
	     LAL_180_PI*earth.geodetic.latitude,
	     LAL_180_PI*earth.geodetic.longitude,
	     1.0*earth.elevation );
    params.zenith = &(earth.geodetic);
    fprintf( stdout, "\n" );
  }

  /* Convert geodetic to geocentric, if required. */
  else if ( gd ) {
    fprintf( stdout, "TERRESTRIAL COORDINATES\n"
	     "Geodetic:   latitude = %6.2f deg, longitude = %6.2f"
	     " deg", LAL_180_PI*earth.geodetic.latitude,
	     LAL_180_PI*earth.geodetic.longitude );
    if ( alt ) {
      fprintf( stdout, ", elevation = %.0fm\n", 1.0*earth.elevation );
      SUB( LALGeodeticToGeocentric( &stat, &earth ), &stat );
      fprintf( stdout,
	       "Geocentric: latitude = %6.2f deg, longitude = %6.2f"
	       " deg, radius = %.0fm\n"
	       "            x = %.0fm, y = %.0fm, z = %.0fm\n",
	       LAL_180_PI*earth.geocentric.latitude,
	       LAL_180_PI*earth.geocentric.longitude, 1.0*earth.radius,
	       1.0*earth.x, 1.0*earth.y, 1.0*earth.z );
    } else
      fprintf( stdout, "\n" );

    /* Assign strcuture for other location-dependent conversions. */
    params.zenith = &(earth.geodetic);
    fprintf( stdout, "\n" );
  }


  /*******************************************************************
   * TIME COORDINATE                                                 *
   *******************************************************************/

  /* Print the time in various formats. */
  if ( t ) {
    INT8 nsec;    /* time as INT8 nanoseconds */
    struct tm date;	/* UTC */
    REAL8 gmst;   /* Greenwich mean sidereal time */

    /* Convert to INT8 seconds and back, just to test things (and get
       the LIGOTimeGPS structure into standard form). */
    nsec = XLALGPSToINT8NS(&gpsTime);
    XLALINT8NSToGPS(&gpsTime, nsec);
    fprintf( stdout, "TIME COORDINATE\n" );

    /* Convert to UTC timestamp. */
    if ( gpsTime.gpsSeconds >= 0 ) {
      CHARVector *timeStamp = NULL; /* date string */
      /* don't bomb if leap second table isn't up to date */
      LALLeapSecAccuracy acc = LALLEAPSEC_LOOSE;
      LALMSTUnitsAndAcc uAcc;
      uAcc.units = MST_DEG;
      uAcc.accuracy = acc;

      fprintf( stdout, "GPS time: %i.%09is\n", gpsTime.gpsSeconds,
	       gpsTime.gpsNanoSeconds );
      XLALGPSToUTC(&date, gpsTime.gpsSeconds);
      SUB( LALCHARCreateVector( &stat, &timeStamp, 32 ), &stat );
      strftime(timeStamp->data, timeStamp->length, "%F %T UTC %a", &date);
      fprintf( stdout, "UTC time: %s\n", timeStamp->data );
      SUB( LALCHARDestroyVector( &stat, &timeStamp ), &stat );

      /* Convert to Greenwich mean sidereal time. */
      SUB( LALGPStoGMST1( &stat, &gmst, &gpsTime, &uAcc ), &stat );
      fprintf( stdout, "Greenwich mean sidereal time: %6.2f deg\n",
	       1.0*gmst );
    } else
      fprintf( stdout, "GPS time: %i.%09is\n", gpsTime.gpsSeconds,
	       -gpsTime.gpsNanoSeconds );

    /* Assign strcuture for other time-dependent conversions. */
    params.gpsTime = &gpsTime;
    fprintf( stdout, "\n" );
  }


  /*******************************************************************
   * CELESTIAL COORDINATES                                           *
   *******************************************************************/

  /* Print the input coordinates. */
  if ( in ) {
    fprintf( stdout, "CELESTIAL COORDINATES\n" );
    switch ( sky.system ) {
    case COORDINATESYSTEM_HORIZON: fprintf( stdout, "Horizon   " );
      break;
    case COORDINATESYSTEM_GEOGRAPHIC: fprintf( stdout, "Geographic" );
      break;
    case COORDINATESYSTEM_EQUATORIAL: fprintf( stdout, "Equatorial" );
      break;
    case COORDINATESYSTEM_ECLIPTIC: fprintf( stdout, "Ecliptic  " );
      break;
    case COORDINATESYSTEM_GALACTIC: fprintf( stdout, "Galactic  " );
      break;
    default: fprintf( stdout, "Unknown   " );
    }
    fprintf( stdout, ": latitude = %6.2f deg, longitude = %6.2f deg\n",
	     LAL_180_PI*sky.latitude, LAL_180_PI*sky.longitude );

    /* Print the output coordinates. */
    if ( out ) {
      SUB( LALConvertSkyCoordinates( &stat, &sky, &sky, &params ),
	   &stat );
      switch ( sky.system ) {
      case COORDINATESYSTEM_HORIZON: fprintf( stdout, "Horizon   " );
	break;
      case COORDINATESYSTEM_GEOGRAPHIC: fprintf( stdout, "Geographic" );
	break;
      case COORDINATESYSTEM_EQUATORIAL: fprintf( stdout, "Equatorial" );
	break;
      case COORDINATESYSTEM_ECLIPTIC: fprintf( stdout, "Ecliptic  " );
	break;
      case COORDINATESYSTEM_GALACTIC: fprintf( stdout, "Galactic  " );
	break;
      default: fprintf( stdout, "Unknown   " );
      }
      fprintf( stdout, ": latitude = %6.2f deg, longitude = %6.2f deg\n",
	       LAL_180_PI*sky.latitude, LAL_180_PI*sky.longitude );
    }
    fprintf( stdout, "\n" );
  }

  /* Clean up and exit. */
  LALCheckMemoryLeaks();
  INFO( SKYCOORDINATESTESTC_MSGENORM );
  return SKYCOORDINATESTESTC_ENORM;
}
