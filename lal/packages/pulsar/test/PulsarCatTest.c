/******************************** <lalVerbatim file="PulsarCatTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{PulsarCatTest.c}}
\label{ss:PulsarCatTest.c}

Tests routines to manipulate pulsar data.

\subsubsection*{Usage}
\begin{verbatim}
PulsarCatTest [-p posepoch ra dec pmra pmdec] [-l site earthfile sunfile] [-h]
              [-t newepoch] [-o outfile] [-d debuglevel] [fepoch f0 [f1 ...]]
\end{verbatim}

\subsubsection*{Description}

This program reads in or randomly generates pulsar parameters, stores
them in a pulsar catalogue structure, and manipulates them based on
command-line arguments.  The following option flags are accepted (in
any order):
\begin{itemize}
\item[\texttt{-p}] The pulsar position at time \verb@posepoch@ is set
to \verb@ra@ radians right ascension and \verb@dec@ radians
declination, with proper motions of \verb@pmra@ and \verb@pmdec@
milliarcseconds per year, respectively.  See below for parsing formats
for \verb@posepoch@.  If the \verb@-p@ option is not specified, a
random source position is generated.
\item[\texttt{-l}] Sets the detector location.  \verb@site@ must be
one of the following character strings: \verb@LHO@, \verb@LLO@,
\verb@VIRGO@, \verb@GEO600@, \verb@TAMA300@, or \verb@CIT40@.
\verb@earthfile@ and \verb@sunfile@ are ephemeris files of the format
expected by \verb@LALInitBarycenter()@.  If the \verb@-l@ option is
not specified, the detector is placed at the solar system barycentre.
\item[\texttt{-h}] Prints usage information and then exits.
\item[\texttt{-t}] Sets the new epoch to which the pulsar data will be
updated.  See below for parsing formats for \verb@newepoch@.  If the
\verb@-t@ option is not given, \verb@-t J2000.0@ is assumed.
\item[\texttt{-o}] Writes the generated time series to the file
\verb@outfile@.  If absent, the routines are exercised, but no output
is written.
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@.  If
absent, level 0 is assumed.
\end{itemize}

Once all valid options are read, the remaining command-line arguments
are taken to be the epoch \verb@fepoch@ at which the pulsar spin
frequency \verb@f@ (in Hz) was measured, plus zero or more frequency
derivatives \verb@f1@$\ldots$ (in Hz${}^{k+1}$ for the $k^\mathrm{th}$
derivative).  If no additional arguments are given, spin timing
information will be omitted.

Meaurement epoch \verb@posepoch@, \verb@newepoch@, and \verb@fepoch@
may be specified either as a \verb@REAL8@ Julian epoch preceded by a
\verb@J@ character (e.g.\ \verb@JD2000.0@), a \verb@REAL8@ number of
Julian days preceded by \verb@JD@ (e.g.\ \verb@JD2451545.0@), or as an
\verb@INT8@ number of GPS nanoseconds with no prefix (e.g.\
\verb@630763213000000000@).  Note that the preceding examples all
refer to noon UTC, January 1, 2000.  Also, note that each Julian epoch
is assumed to be exactly 365.25 Julian days, so J2001.0 corresponds to
18:00 UTC, January 1, 2001.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define PULSARCATTESTC_ENORM 0
#define PULSARCATTESTC_ESUB  1
#define PULSARCATTESTC_EARG  2
#define PULSARCATTESTC_EVAL  3
#define PULSARCATTESTC_EMEM  4
#define PULSARCATTESTC_EFILE 5

#define PULSARCATTESTC_MSGENORM "Normal exit"
#define PULSARCATTESTC_MSGESUB  "Subroutine failed"
#define PULSARCATTESTC_MSGEARG  "Error parsing arguments"
#define PULSARCATTESTC_MSGEVAL  "Input argument out of valid range"
#define PULSARCATTESTC_MSGEMEM  "Out of memory"
#define PULSARCATTESTC_MSGEFILE "Could not open file"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

This routine simply parses the input arguments, stuffs the data into a
\verb@PulsarCatNode@ structure, and then calls
\verb@LAUpdatePulsarCatNode()@ to update it to the new epoch.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALMalloc()                     LALFree()
LALDCreateVector()              LALDDestroyVector()
LALCreateRandomParams()         LALDestroyRandomParams()
LALUniformDeviate()             LALInitBarycenter()
LALGPStoINT8()                  LALINT8toGPS()
LALStringToD()                  LALStringToI8()
LALLeapSec()                    LALUpdatePulsarCat()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{PulsarCatTestCV}}

******************************************************* </lalLaTeX> */

#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/DetectorSite.h>
#include <lal/PulsarCat.h>

NRCSID( PULSARCATTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;
#define J2000GPS 630763213 /* J2000.0 epoch in GPS seconds */
#define J2000JD    2451545 /* J2000.0 epoch in Julian days */

/* Usage format string. */
#define USAGE "Usage: %s [-p posepoch ra dec pmra pmdec]\n" \
"\t[-l site earthfile sunfile] [-h]\n" \
"\t[-t newepoch] [-o outfile] [-d debuglevel] [fepoch f0 [f1 ...]]\n"

/* Maximum output message length. */
#define MSGLEN (1024)

/* Upper cutoff frequency for default detector response function. */
#define FSTOP (16384.0)

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, PULSARCATTESTC,                           \
		 statement ? statement : "", (msg) );                \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 PULSARCATTESTC, (statement) );                      \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 PULSARCATTESTC, (statement) );                      \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( PULSARCATTESTC_ESUB, PULSARCATTESTC_MSGESUB,                \
	 "Function call \"" #func "\" failed:" );                    \
  return PULSARCATTESTC_ESUB;                                        \
}                                                                    \
while (0)

#define CHECKVAL( val, lower, upper )                                \
do                                                                   \
if ( ( (val) < (lower) ) || ( (val) > (upper) ) )                    \
{                                                                    \
  ERROR( PULSARCATTESTC_EVAL, PULSARCATTESTC_MSGEVAL,                \
         "Value of " #val " out of range:" );                        \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( #val " = %f, range = [%f,%f]\n", (REAL8)(val),    \
                   (REAL8)(lower), (REAL8)(upper) );                 \
  return PULSARCATTESTC_EVAL;                                        \
}                                                                    \
while (0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

/* Prototype for a routine to parse epoch strings. */
static void
ParseEpoch( LALStatus *stat, LIGOTimeGPS *epoch, const CHAR *string );


int
main(int argc, char **argv)
{
  static LALStatus stat;       /* status structure */
  int arg;                     /* command-line argument counter */
  BOOLEAN posGiven = 0;        /* whether position was given */
  CHAR *outfile = NULL;        /* name of outfile */
  CHAR *earthfile = NULL;      /* name of earthfile */
  CHAR *sunfile = NULL;        /* name of sunfile */
  CHAR *site = NULL;           /* name of detector site */
  PulsarCatNode node;          /* pulsar data structure */
  LIGOTimeGPS epoch;           /* epoch of update */
  LALPlaceAndGPS detectorTime; /* epoch and detector site */
  EphemerisData *edat = NULL;  /* detector ephemerides */

  /* First set up some defaults. */
  memset( &node, 0, sizeof(PulsarCatNode) );
  epoch.gpsSeconds = J2000GPS;
  epoch.gpsNanoSeconds = 0;
  detectorTime.p_gps = &epoch;
  detectorTime.p_detector = NULL;

  /*******************************************************************
   * ARGUMENT PARSING (arg stores the current position)              *
   *******************************************************************/

  /* Start reading options. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse position option. */
    if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 5 ) {
	arg++;
	SUB( ParseEpoch( &stat, &(node.posepoch), argv[arg++] ),
	     &stat );
	node.pos.system = COORDINATESYSTEM_EQUATORIAL;
	node.pos.longitude = atof( argv[arg++] );
	node.pos.latitude = atof( argv[arg++] );
	node.pm.system = COORDINATESYSTEM_EQUATORIAL;
	node.pm.longitude = atof( argv[arg++] );
	node.pm.latitude = atof( argv[arg++] );
	posGiven = 1;
      } else {
	ERROR( PULSARCATTESTC_EARG, PULSARCATTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PULSARCATTESTC_EARG;
      }
    }
    /* Parse detector location option. */
    else if ( !strcmp( argv[arg], "-l" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	site = argv[arg++];
	earthfile = argv[arg++];
	sunfile = argv[arg++];
      } else {
	ERROR( PULSARCATTESTC_EARG, PULSARCATTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PULSARCATTESTC_EARG;
      }
    }
    /* Parse output time option. */
    else if ( !strcmp( argv[arg], "-t" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	SUB( ParseEpoch( &stat, &epoch, argv[arg++] ), &stat );
      } else {
	ERROR( PULSARCATTESTC_EARG, PULSARCATTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PULSARCATTESTC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      } else {
	ERROR( PULSARCATTESTC_EARG, PULSARCATTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PULSARCATTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      } else {
	ERROR( PULSARCATTESTC_EARG, PULSARCATTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PULSARCATTESTC_EARG;
      }
    }
    /* Parse help option. */
    else if ( !strcmp( argv[arg], "-h" ) ) {
      INFO( PULSARCATTESTC_MSGENORM );
      LALPrintError( USAGE, *argv );
      return PULSARCATTESTC_ENORM;
    }
    /* Any additional command-line options are timing parameters. */
    else if ( argc >= arg + 2 ) {
      REAL8 *f; /* pointer to timing data */
      SUB( ParseEpoch( &stat, &(node.fepoch), argv[arg++] ), &stat );
      SUB( LALDCreateVector( &stat, &(node.f), argc - arg ), &stat );
      for ( f = node.f->data; arg < argc; f++ )
	*f = atof( argv[arg++] );
    }
    /* Unknown option. */
    else {
      ERROR( PULSARCATTESTC_EARG, PULSARCATTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return PULSARCATTESTC_EARG;
    }
  }

  /*******************************************************************
   * SETUP                                                           *
   *******************************************************************/

  /* If position was not given, generate one randomly. */
  if ( !posGiven ) {
    RandomParams *params = NULL; /* random number parameters */
    REAL4 x;                     /* random deviate */
    SUB( LALCreateRandomParams( &stat, &params, 0 ), &stat );
    node.pos.system = COORDINATESYSTEM_EQUATORIAL;
    SUB( LALUniformDeviate( &stat, &x, params ), &stat );
    node.pos.longitude = LAL_TWOPI*x;
    SUB( LALUniformDeviate( &stat, &x, params ), &stat );
    node.pos.latitude = LAL_PI*( x - 0.5 );
    SUB( LALDestroyRandomParams( &stat, &params ), &stat );
    node.posepoch.gpsSeconds = J2000GPS;
    node.posepoch.gpsNanoSeconds = 0;
  }

  /* If the detector was specified, set up the detector position. */
  if ( site ) {
    UINT4 i; /* site index */
    if ( !strcmp( site, "LHO" ) )
      i = LALDetectorIndexLHODIFF;
    else if ( !strcmp( site, "LLO" ) )
      i = LALDetectorIndexLLODIFF;
    else if ( !strcmp( site, "VIRGO" ) )
      i = LALDetectorIndexVIRGODIFF;
    else if ( !strcmp( site, "GEO600" ) )
      i = LALDetectorIndexGEO600DIFF;
    else if ( !strcmp( site, "TAMA300" ) )
      i = LALDetectorIndexTAMA300DIFF;
    else if ( !strcmp( site, "CIT40" ) )
      i = LALDetectorIndexCIT40DIFF;
    else {
      ERROR( PULSARCATTESTC_EVAL, PULSARCATTESTC_MSGEVAL,
	     "Unrecognized site:" );
      if ( lalDebugLevel & LALERROR )
	LALPrintError( "%s", site );
      return PULSARCATTESTC_EVAL;
    }
    detectorTime.p_detector =
      (LALDetector *)LALMalloc( sizeof(LALDetector) );
    if ( !(detectorTime.p_detector) ) {
      ERROR( PULSARCATTESTC_EMEM, PULSARCATTESTC_MSGEMEM, 0 );
      return PULSARCATTESTC_EMEM;
    }
    *(detectorTime.p_detector) = lalCachedDetectors[i];

    /* Read ephemerides. */
    edat = (EphemerisData *)LALMalloc( sizeof(EphemerisData) );
    if ( !(edat) ) {
      ERROR( PULSARCATTESTC_EMEM, PULSARCATTESTC_MSGEMEM, 0 );
      return PULSARCATTESTC_EMEM;
    }
    edat->ephiles.earthEphemeris = earthfile;
    edat->ephiles.sunEphemeris = sunfile;
    SUB( LALInitBarycenter( &stat, edat ), &stat );
  }


  /*******************************************************************
   * CONVERSION                                                      *
   *******************************************************************/

  /* Perform update. */
  SUB( LALUpdatePulsarCatNode( &stat, &node, &detectorTime, edat ),
       &stat );

  /* If output was requested, print output. */
  if ( outfile ) {
    PulsarCatNode *here = &node;
    FILE *fp = fopen( outfile, "w" );
    if ( fp == NULL ) {
      ERROR( PULSARCATTESTC_EFILE, PULSARCATTESTC_MSGEFILE, outfile );
      return PULSARCATTESTC_EFILE;
    }
    while ( here ) {
      CompanionNode *companion = here->companion; /* companion data */
      UINT4 compNo = 0; /* companion number */
      INT8 ep; /* an epoch */
      if ( here->jname[0] != '\0' ) {
	fprintf( fp, "PULSAR %s", here->jname );
	if ( here->bname[0] != '\0' )
	  fprintf( fp, " (%s)", here->bname );
	fprintf( fp, "\n" );
      } else if ( here->bname[0] != '\0' )
	fprintf( fp, "PULSAR %s\n", here->bname );
      else
	fprintf( fp, "PULSAR (Unknown)\n" );
      fprintf( fp, "ra    = %23.16e rad\n", here->pos.longitude );
      fprintf( fp, "dec   = %23.16e rad\n", here->pos.latitude );
      fprintf( fp, "pmra  = %23.16e rad/s\n", here->pm.longitude );
      fprintf( fp, "pmdec = %23.16e rad/s\n", here->pm.latitude );
      SUB( LALGPStoINT8( &stat, &ep, &(here->posepoch) ), &stat );
      fprintf( fp, "posepoch = %lli\n", ep );
      if ( here->f ) {
	UINT4 i; /* an index */
	fprintf( fp, "f = %23.16e Hz\n", here->f->data[0] );
	for ( i = 1; i < here->f->length; i++ )
	  fprintf( fp, "    %23.16e Hz^%i\n", here->f->data[i], i+1 );
	SUB( LALGPStoINT8( &stat, &ep, &(here->fepoch) ), &stat );
	fprintf( fp, "fepoch = %lli\n", ep );
      }
      if ( here->dist > 0.0 )
	fprintf( fp, "dist = %15.8e m\n", here->dist );
      if ( here->dmin > 0.0 )
	fprintf( fp, "dmin = %15.8e m\n", here->dmin );
      if ( here->dmax > 0.0 )
	fprintf( fp, "dmax = %15.8e m\n", here->dmax );
      if ( here->lcode >= 'a' && here->lcode <= 'd' )
	fprintf( fp, "lcode = %c\n", here->lcode );
      if ( here->ucode >= 'a' && here->ucode <= 'd' )
	fprintf( fp, "ucode = %c\n", here->ucode );
      fprintf( fp, "typecode = %i\n", here->typecode );
      while ( companion ) {
	fprintf( fp, "Companion %i\n", ++compNo );
	SUB( LALGPStoINT8( &stat, &ep, &(companion->epoch) ), &stat );
	fprintf( fp, "  epoch = %lli\n", ep );
	fprintf( fp, "  x     = %23.16e s\n", companion->x );
	fprintf( fp, "  p     = %23.16e s\n", companion->p );
	fprintf( fp, "  pDot  = %23.16e\n", companion->pDot );
	fprintf( fp, "  w     = %23.16e rad\n", companion->w );
	fprintf( fp, "  wDot  = %23.16e rad/s\n", companion->wDot );
	fprintf( fp, "  ecc   = %23.16e\n", companion->ecc );
	fprintf( fp, "  gamma = %23.16e s\n", companion->gamma );
	fprintf( fp, "  sin   = %23.16e\n", companion->sin );
	fprintf( fp, "  r     = %23.16e\n", companion->r );
	companion = companion->next;
      }
      here = here->next;
    }
    fclose( fp );
  }


  /*******************************************************************
   * CLEANUP                                                         *
   *******************************************************************/

  /* Destroy any allocated memory. */
  if ( node.f ) {
    SUB( LALDDestroyVector( &stat, &(node.f) ), &stat );
  }
  if ( edat ) {
    if ( edat->ephemE )
      LALFree( edat->ephemE );
    if ( edat->ephemS )
      LALFree( edat->ephemS );
    LALFree( edat );
  }
  if ( detectorTime.p_detector )
    LALFree( detectorTime.p_detector );

  /* Done! */
  LALCheckMemoryLeaks();
  INFO( PULSARCATTESTC_MSGENORM );
  return PULSARCATTESTC_ENORM;
}


/* A routine to parse epoch strings. */
static void
ParseEpoch( LALStatus *stat, LIGOTimeGPS *epoch, const CHAR *string )
{
  INT8 gpsNan;  /* number of GPS nanoseconds */
  CHAR *endptr; /* pointer to end of parsed data */

  INITSTATUS( stat, "ParseEpoch", PULSARCATTESTC );
  ATTATCHSTATUSPTR( stat );

  /* Parse Julian days, or Julian epochs (converted to days). */
  if ( string[0] == 'J' ) {
    REAL8 julianDay;            /* Julian date */
    INT4 leap1, leap2;          /* number of leap seconds to date */
    LALLeapSecFormatAndAcc acc; /* accuracy of leap second computation */
    if ( string[1] == 'D' ) {
      TRY( LALStringToD( stat->statusPtr, &julianDay, string+2,
			 &endptr ), stat );
      if ( endptr == string+2 ) {
	ABORT( stat, PULSARCATTESTC_EARG, PULSARCATTESTC_MSGEARG );
      }
    } else {
      TRY( LALStringToD( stat->statusPtr, &julianDay, string+1,
			 &endptr ), stat );
      if ( endptr == string+1 ) {
	ABORT( stat, PULSARCATTESTC_EARG, PULSARCATTESTC_MSGEARG );
      }
      julianDay -= 2000.0;
      julianDay *= 365.25;
      julianDay += J2000JD;
    }

    /* Convert Julian days to GPS nanoseconds. */
    acc.accuracy = LALLEAPSEC_STRICT;
    acc.format = LALLEAPSEC_GPSUTC;
    gpsNan = (INT8)( ( julianDay - 2444244.5 )*(8.64e13L) );
    TRY( LALINT8toGPS( stat->statusPtr, epoch, &gpsNan ), stat );
    leap2 = 0;
    do {
      leap1 = leap2;
      TRY( LALLeapSecs( stat->statusPtr, &leap2, epoch, &acc ),
	   stat );
      epoch->gpsSeconds += leap2 - leap1;
    } while ( leap2 != leap1 );
  }

  /* Parse GPS times directly. */
  else {
    TRY( LALStringToI8( stat->statusPtr, &gpsNan, string, &endptr ),
	 stat );
    if ( endptr == string ) {
      ABORT( stat, PULSARCATTESTC_EARG, PULSARCATTESTC_MSGEARG );
    }
    TRY( LALINT8toGPS( stat->statusPtr, epoch, &gpsNan ), stat );
  }

  /* That's all. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
