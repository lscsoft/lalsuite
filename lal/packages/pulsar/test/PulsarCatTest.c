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
              [-t newepoch] [-i infile] [-o outfile] [-d debuglevel]
              [fepoch f0 [f1 ...]]
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
radians per second, respectively.  See below for parsing formats for
\verb@posepoch@.  If the \verb@-p@ option is not specified, a random
source position is generated.
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
\item[\texttt{-i}] Reads the pulsar data from the file \verb@infile@,
whose format is described below.  This overrides the position and spin
information read from the command line.  If the name \verb@stdin@ is
given, it will read from standard input, \emph{not} from a file named
\verb@stdin@.
\item[\texttt{-o}] Writes the pulsar data to the file \verb@outfile@.
If the name \verb@stdout@ or \verb@stderr@ is given, it will write to
standard output or standard error (respectively), \emph{not} to a file
named \verb@stdout@ or \verb@stderr@.  If the \verb@-o@ option is not
given, the routines are exercised, but no output is written.
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

If an input file is specified, it should consist of a header line
that, when tokenized, can be parsed by \verb@LALReadPulsarCatHead()@,
followed by one or more lines of pulsar data parseable (when
tokenized) by \verb@LALReadPulsarCatLine()@.  Blank lines (with no
tokens) or divider lines (with only one token) will be skipped.

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

If the \verb@-i@ option is given, the corresponding file is opened and
read by \verb@LALCHARReadVectorSequence()@, then each line is
tokenized by \verb@LALCreateTokenList()@.

Output via the \verb@-o@ option is in a custom human-readable format,
which should be easy to figure out.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALMalloc()                     LALFree()
LALDCreateVector()              LALDDestroyVector()
LALCreateRandomParams()         LALDestroyRandomParams()
LALUniformDeviate()             LALInitBarycenter()
LALGPStoINT8()                  LALINT8toGPS()
LALCHARReadVectorSequence()     LALCHARDestroyVectorSequence()
LALCreateTokenList()            LALDestroyTokenList()
LALReadPulsarCatHead()          LALReadPulsarCatLine()
LALStringToD()                  LALStringToI8()
LALLeapSec()                    LALUpdatePulsarCat()
LALSnprintf()
\end{verbatim}

\subsubsection*{Notes}

At present the routine is kludged up to ignore pulsar position and
frequency information from the command line, using hardwired
parameters for \mbox{PSR J0034-0534} instead.  It can still override
this with the \verb@-i@ option.

\vfill{\footnotesize\input{PulsarCatTestCV}}

******************************************************* </lalLaTeX> */

#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/StringInput.h>
#include <lal/StreamInput.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
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
"\t[-t newepoch] [-i infile] [-o outfile] [-d debuglevel]\n" \
"\t[fepoch f0 [f1 ...]]\n"

/* List of whitespace characters used for tokenizing input. */
#define WHITESPACES " \b\f\n\r\t\v"

/* Maximum length of fprintderr() format or error/info message
   strings. */
#define MAXLEN 1024

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

/* Prototype for a routine to print floating-point numbers with
   uncertainties. */
int
fprintderr( FILE *fp, REAL8 x, REAL8 dx );

/* Prototype for a routine to print GPS epochs. */
int
fprintepoch( FILE *fp, LIGOTimeGPS epoch );


int
main(int argc, char **argv)
{
  static LALStatus stat;       /* status structure */
  int arg;                     /* command-line argument counter */
  BOOLEAN posGiven = 0;        /* whether position was given */
  CHAR *infile = NULL;         /* name of infile */
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
    /* Parse input file option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	infile = argv[arg++];
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

    /* Instead, kludge it up for PSR J1939+2134, B1937+21.
    memcpy( node.jname, "J1939+2134", 11*sizeof(CHAR) );
    memcpy( node.bname, "B1937+21", 9*sizeof(CHAR) );
    node.pos.system  = COORDINATESYSTEM_EQUATORIAL;
    node.dpos.system = COORDINATESYSTEM_EQUATORIAL;
    node.pm.system   = COORDINATESYSTEM_EQUATORIAL;
    node.dpm.system  = COORDINATESYSTEM_EQUATORIAL;
    node.pos.latitude   = 0.3766960555;
    node.dpos.latitude  = 0.0000000034;
    node.pos.longitude  = 5.1471621406;
    node.dpos.longitude = 0.0000000015;
    node.pm.latitude   = -7.13E-17;
    node.dpm.latitude  =  0.14E-17;
    node.pm.longitude  = -2.00E-17;
    node.dpm.longitude =  0.12E-17;
    SUB( ParseEpoch( &stat, &(node.posepoch), "JD2447500.0" ), &stat );
    SUB( LALDCreateVector( &stat, &(node.f), 2 ), &stat );
    SUB( LALDCreateVector( &stat, &(node.df), 2 ), &stat );
    node.f->data[0]  = 641.92825287201;
    node.df->data[0] =   0.00000000017;
    node.f->data[1]  =  -4.3313E-14;
    node.df->data[1] =   0.0013E-14;
    SUB( ParseEpoch( &stat, &(node.fepoch), "JD2450100.0" ), &stat ); */

    /* Instead, kludge it up for PSR J0034-0534.
    memcpy( node.jname, "J0034-0534", 11*sizeof(CHAR) );
    node.pos.system  = COORDINATESYSTEM_EQUATORIAL;
    node.dpos.system = COORDINATESYSTEM_EQUATORIAL;
    node.pm.system   = COORDINATESYSTEM_EQUATORIAL;
    node.dpm.system  = COORDINATESYSTEM_EQUATORIAL;
    node.pos.latitude   = -0.097334152;
    node.dpos.latitude  =  0.000000097;
    node.pos.longitude  =  0.149940232;
    node.dpos.longitude =  0.000000036;
    node.pm.latitude   = -0.5E-15;
    node.dpm.latitude  =  3.5E-15;
    node.pm.longitude  =  2.3E-15;
    node.dpm.longitude =  1.7E-15;
    SUB( ParseEpoch( &stat, &(node.posepoch), "JD2449550.0" ), &stat );
    SUB( LALDCreateVector( &stat, &(node.f), 2 ), &stat );
    SUB( LALDCreateVector( &stat, &(node.df), 2 ), &stat );
    node.f->data[0]  = 532.71343832081;
    node.df->data[0] =   0.00000000015;
    node.f->data[1]  =  -1.436E-15;
    node.df->data[1] =   0.009E-15;
    SUB( ParseEpoch( &stat, &(node.fepoch), "JD2449550.0" ), &stat ); */
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

  /* If an input file was specified, it overrides the preceding
     position and frequency information. */
  if ( infile ) {
    UINT4 i = 0;                     /* line number */
    UINT4 n = 0;                     /* number of pulsars read */
    CHAR msg[MAXLEN];                /* error/info message string */
    CHARVectorSequence *file = NULL; /* input file */
    TokenList *list = NULL;          /* input line parsed into tokens */
    INT4 indx[PULSARCATINDEX_NUM];   /* ordering of tokens in line */
    PulsarCatNode *here = &node;     /* structure for current line */
    FILE *fp = NULL;                 /* input file pointer */

    /* Read in file. */
    if ( !strcmp( infile, "stdin" ) )
      fp = stdin;
    else
      fp = fopen( infile, "r" );
    if ( fp == NULL ) {
      ERROR( PULSARCATTESTC_EFILE, PULSARCATTESTC_MSGEFILE, infile );
      return PULSARCATTESTC_EFILE;
    }
    SUB( LALCHARReadVectorSequence( &stat, &file, fp ), &stat );
    fclose( fp );

    /* Read header, skipping blank or divider lines. */
    SUB( LALCreateTokenList( &stat, &list, file->data, WHITESPACES ),
	 &stat );
    while ( list->nTokens < 3 && ++i < file->length ) {
      SUB( LALDestroyTokenList( &stat, &list ), &stat );
      SUB( LALCreateTokenList( &stat, &list, file->data +
			       i*file->vectorLength, WHITESPACES ),
	   &stat );
    }
    if ( i >= file->length ) {
      WARNING( "No header found in input file" );
    } else {
      SUB( LALReadPulsarCatHead( &stat, indx, list ), &stat );
    }
    SUB( LALDestroyTokenList( &stat, &list ), &stat );

    /* Read body, skipping blank or divider lines. */
    while ( ++i < file->length ) {
      SUB( LALCreateTokenList( &stat, &list, file->data +
			       i*file->vectorLength, WHITESPACES ),
	   &stat );
      while ( list->nTokens < 3 && ++i < file->length ) {
	SUB( LALDestroyTokenList( &stat, &list ), &stat );
	SUB( LALCreateTokenList( &stat, &list, file->data +
				 i*file->vectorLength, WHITESPACES ),
	     &stat );
      }
      if ( i < file->length ) {
	if ( !( here->next = (PulsarCatNode *)
		LALMalloc( sizeof(PulsarCatNode) ) ) ) {
	  ERROR( PULSARCATTESTC_EMEM, PULSARCATTESTC_MSGEMEM, 0 );
	  return PULSARCATTESTC_EMEM;
	}
	memset( here->next, 0, sizeof(PulsarCatNode) );

	/* Ignore lines that don't parse, moving on to the next. */
	LALReadPulsarCatLine( &stat, here->next, list, indx );
	if ( stat.statusCode ) {
	  LALSnprintf( msg, MAXLEN, "Error reading line %i of %s:"
		       " skipping", i, infile );
	  ERROR( 0, msg, 0 );
	  if ( stat.statusPtr ) {
	    FREESTATUSPTR( &stat );
	  }
	  memset( &stat, 0, sizeof(LALStatus) );
	  SUB( LALDestroyPulsarCat( &stat, &(here->next) ), &stat );
	} else {
	  here = here->next;
	  n++;
	}
      }
      SUB( LALDestroyTokenList( &stat, &list ), &stat );
    }
    LALSnprintf( msg, MAXLEN, "File %s had %i lines and %i parseable"
		 " pulsars", infile, i, n );
    INFO( msg );
    SUB( LALCHARDestroyVectorSequence( &stat, &file ), &stat );
  }


  /*******************************************************************
   * CONVERSION                                                      *
   *******************************************************************/

  /* Perform update. */
  SUB( LALUpdatePulsarCat( &stat, &node, &detectorTime, edat ),
       &stat );

  /* If output was requested, print output. */
  if ( outfile ) {
    PulsarCatNode *here = &node; /* current position in catalogue */
    FILE *fp = NULL;             /* output file pointer */
    if ( !strcmp( outfile, "stdout" ) )
      fp = stdout;
    else if ( !strcmp( outfile, "stderr" ) )
      fp = stderr;
    else
      fp = fopen( outfile, "w" );
    if ( fp == NULL ) {
      ERROR( PULSARCATTESTC_EFILE, PULSARCATTESTC_MSGEFILE, outfile );
      return PULSARCATTESTC_EFILE;
    }
    if ( infile )
      here = here->next;

    /* First print epoch information. */
    if ( site ) {
      fprintf( fp, "%s time = ", site );
      fprintepoch( fp, epoch );
      fprintf( fp, "\n" );
    }

    /* Now print pulsar information. */
    while ( here ) {
      CompanionNode *companion = here->companion; /* companion data */
      UINT4 compNo = 0; /* companion number */
      fprintf( fp, "\n" );
      if ( here->jname[0] != '\0' ) {
	fprintf( fp, "PULSAR %s", here->jname );
	if ( here->bname[0] != '\0' )
	  fprintf( fp, " (%s)", here->bname );
	fprintf( fp, "\n" );
      } else if ( here->bname[0] != '\0' )
	fprintf( fp, "PULSAR %s\n", here->bname );
      else
	fprintf( fp, "PULSAR (Unknown)\n" );
      fprintf( fp, "epoch = " );
      fprintepoch( fp, here->posepoch );
      fprintf( fp, "\n" );
      fprintf( fp, "ra    = " );
      fprintderr( fp, here->pos.longitude, here->dpos.longitude );
      fprintf( fp, " rad\n" );
      fprintf( fp, "dec   = " );
      fprintderr( fp, here->pos.latitude, here->dpos.latitude );
      fprintf( fp, " rad\n" );
      fprintf( fp, "pmra  = " );
      fprintderr( fp, here->pm.longitude, here->dpm.longitude );
      fprintf( fp, " rad/s\n" );
      fprintf( fp, "pmdec = " );
      fprintderr( fp, here->pm.latitude, here->dpm.latitude );
      fprintf( fp, " rad/s\n" );
      /* SUB( LALGPStoINT8( &stat, &ep, &(here->posepoch) ), &stat );
	 fprintf( fp, "posepoch = %lli ns\n", ep ); */
      if ( here->f ) {
	UINT4 i; /* an index */
	fprintf( fp, "f0    = " );
	fprintderr( fp, 2.0*here->f->data[0], 2.0*here->df->data[0] );
	fprintf( fp, " Hz\n" );
	for ( i = 1; i < here->f->length; i++ ) {
	  fprintf( fp, "f%i    = ", i );
	  fprintderr( fp, 2.0*here->f->data[i], 2.0*here->df->data[i] );
	  fprintf( fp, " Hz^%i\n", i+1 );
	}
	/* SUB( LALGPStoINT8( &stat, &ep, &(here->fepoch) ), &stat );
	   fprintf( fp, "fepoch = %lli ns\n", ep ); */
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
      /* fprintf( fp, "typecode = %i\n", here->typecode ); */
      while ( companion ) {
	fprintf( fp, "Companion %i\n", ++compNo );
	/* SUB( LALGPStoINT8( &stat, &ep, &(companion->epoch) ), &stat );
	   fprintf( fp, "  epoch = %lli ns\n", ep ); */
	fprintf( fp, "  epoch = " );
	fprintepoch( fp, companion->epoch );
	fprintf( fp, "\n" );
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
  if ( node.next ) {
    SUB( LALDestroyPulsarCat( &stat, &(node.next) ), &stat );
  }
  if ( node.f ) {
    SUB( LALDDestroyVector( &stat, &(node.f) ), &stat );
  }
  if ( node.df ) {
    SUB( LALDDestroyVector( &stat, &(node.df) ), &stat );
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


int
fprintderr( FILE *fp, REAL8 x, REAL8 dx ) {
  CHAR format[MAXLEN]; /* format string for fprintf() */
  INT4 gsd = 0;        /* place index of greatest significant digit */
  INT4 lsd = 0;        /* place index of least significant digit */
  REAL8 norm;          /* normalization factor */

  /* Compute gsd, lsd, y, and dy. */
  if ( dx < LAL_REAL8_EPS*fabs( x ) )
    dx = 0.0;
  if ( dx > 0.0 ) {
    REAL8 lsdd = log( 0.5*dx )/log( 10.0 );
    if ( lsdd >= 0.0 )
      lsd = (INT4)( lsdd );
    else
      lsd = (INT4)( lsdd ) - 1;
  }
  if ( x != 0.0 ) {
    REAL8 gsdd = log( fabs( x ) )/log( 10.0 );
    if ( gsdd >= 0.0 )
      gsd = (INT4)( gsdd );
    else
      gsd = (INT4)( gsdd ) - 1;
  }

  /* If x is zero, format is determined entirely by dx. */
  if ( x == 0.0 ) {
    if ( dx <= 0.0 )
      return fprintf( fp, "0" );
    if ( abs( lsd ) > 3 ) {
      norm = pow( 10.0, -lsd );
      return fprintf( fp, "( 0 +/- %.0f )e%+i", dx*norm, lsd );
    }
    if ( lsd <= 0 ) {
      LALSnprintf( format, MAXLEN, "%%.%if +/- %%.%if", -lsd, -lsd );
      return fprintf( fp, format, 0.0, dx );
    }
    norm = pow( 10.0, -lsd );
    LALSnprintf( format, MAXLEN, "0 +/- %%.0f%%0%ii", lsd );
    return fprintf( fp, format, dx*norm, 0 );
  }

  /* If number is exact to 8-byte precision, print it as such. */
  if ( dx <= 0.0 ) {
    if ( abs( gsd ) > 3 )
      return fprintf( fp, "%.16e", x );
    LALSnprintf( format, MAXLEN, "%%.%if", 16 - gsd );
    return fprintf( fp, format, x );
  }

  /* Otherwise, format depends on x and dx. */
  if ( gsd < lsd )
    gsd = lsd;
  if ( lsd > 3 || gsd < -3 ) {
    norm = pow( 10.0, -gsd );
    LALSnprintf( format, MAXLEN, "( %%.%if +/- %%.%if )e%+i",
		 gsd - lsd, gsd - lsd, gsd );
    return fprintf( fp, format, x*norm, dx*norm );
  }
  if ( lsd <= 0 ) {
    LALSnprintf( format, MAXLEN, "%%.%if +/- %%.%if", -lsd, -lsd );
    return fprintf( fp, format, x, dx );
  }
  norm = pow( 10.0, -lsd );
  LALSnprintf( format, MAXLEN, "%%.0f%%0%ii +/- %%.0f%%0%ii", lsd,
	       lsd );
  return fprintf( fp, format, x*norm, 0, dx*norm, 0 );
}


int
fprintepoch( FILE *fp, LIGOTimeGPS epoch ) {
  while ( epoch.gpsNanoSeconds >= 1000000000 ) {
    epoch.gpsSeconds += 1;
    epoch.gpsNanoSeconds -= 1000000000;
  }
  while ( epoch.gpsNanoSeconds < 0 ) {
    epoch.gpsSeconds -= 1;
    epoch.gpsNanoSeconds += 1000000000;
  }
  if ( epoch.gpsSeconds < 0 ) {
    epoch.gpsSeconds += 1;
    epoch.gpsNanoSeconds -= 1000000000;
  }
  return fprintf( fp, "%i s %09i ns", epoch.gpsSeconds,
		  abs( epoch.gpsNanoSeconds ) );
}
