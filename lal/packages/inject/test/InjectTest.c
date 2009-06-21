/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

/*********************************** <lalVerbatim file="InjectTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Program \texttt{InjectTest.c}}
\label{ss:InjectTest.c}

Injects an inspiral signal into detector noise.

\subsubsection*{Usage}
\begin{verbatim}
InjectTest [-s sourcefile] [-e earthfile sunfile] [-r randomseed]
           [-d debuglevel] [-h] indxfile
\end{verbatim}

\subsubsection*{Description}

This program generates Galactic inspiral waveform signals and injects
them into ADC data.  The following option flags are accepted:
\begin{itemize}
\item[\texttt{-s}] Reads source data from the file \verb@sourcefile@,
whose format is given below.  If not specified, no injections are
performed.
\item[\texttt{-e}] Specifies ephemeris files giving the location of
the Earth and Sun with respect to a barycentric reference point for
arrival times, whose format is as required by
\verb@LALInitBarycenter()@ in the \verb@support@ package.  If not
specified, the Earth's centre is treated as the barycentric reference
point (this is not consistent for long signal durations).
\item[\texttt{-r}] Sets the random number seed to \verb@randomseed@.
If not specified, the seed is gerenated from the current time.
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@.  If not
specified, level 0 is assumed.
\item[\texttt{-h}] Prints out the usage message and exits.
\end{itemize}
After parsing these recognized options, there must be one remaining
argument: the name of an index file \verb@indxfile@ specifying the
data input and output files, along with detector information
associated with each file; the format is given below.

\paragraph{Format for \texttt{sourcefile}:} The source file consists
of one or more lines, each representing a particular source to be
injected.  Each line consists of the name of a routine for generating
a waveform, followed by a set of numerical arguments required for
setting the parameters for that generator.  The generators currently
supported are as follows:
\begin{verbatim}
LALGeneratePPNInspiral tc m1 m2 d inc ra dec psi phic fi ff
LALGenerateTaylorCW    t0 t1 t2 a1 a2 ra dec psi phi0 f0 [f1 [...]]
LALGenerateSpinOrbitCW t0 t1 t2 a1 a2 ra dec psi phi0 f0 [f1 [...]] arg udot rp e
\end{verbatim}
where \verb@tc@ is the time of coalescence, \verb@t0@ is a reference
time where system properties are specified, \verb@t1@ and \verb@t2@
are start and stop times for the signal (all given as \verb@INT8@ GPS
nanoseconds), \verb@m1@ and \verb@m2@ are the component masses of a
binary system (\verb@REAL4@ solar masses), \verb@d@ is the distance to
the system (\verb@REAL4@ Mpc), \verb@inc@ is the inclination of the
system (\verb@REAL4@ degrees), \verb@a1@ and \verb@a2@ are intrinsic
GW amplitudes of the + and $\times$ polarizations (\verb@REAL4@
strain), \verb@ra@ and \verb@dec@ are the right ascension and
declination of the system (\verb@REAL8@ degrees), \verb@psi@ is the
polarization angle of the system (\verb@REAL4@ degrees), \verb@phic@
is the wave phase at coalescence (\verb@REAL4@ degrees), \verb@phi0@
is the wave phase at the reference time (\verb@REAL8@ degrees),
\verb@fi@ and \verb@ff@ are the start and stop frequencies of the
waveform (\verb@REAL4@ Hz), \verb@f0@ is the wave frequency at the
reference time (\verb@REAL8@ Hz),
\verb@f1@$,\ldots,$\verb@f@${}_k,\ldots$ are the (optional)
frequency-normalized spindown coefficience (\verb@REAL8@ s${}^{-k}$),
\verb@arg@ is the argument of the source's orbital periapsis
(\verb@REAL4@ degrees), \verb@udot@ is the orbital angular speed of
the source at periapsis (\verb@REAL4@ Hz), \verb@rp@ is the normalized
projected orbital periapsis distance $(r_p/c)\sin i$ of the source
(\verb@REAL4@ s), and \verb@e@ is the orbital eccentricity of the
source (\verb@REAL4@),

\paragraph{Format for \texttt{indxfile}:} The index file consists
of one or more lines, each representing a stretch of data to be
injected.  Each line consists of four strings: The site name, the
response function file name, the input file name, and the output file
name (which may be the same as the input).  These are treated as
follows:

The detector name may be one of \verb@LHO@, \verb@LLO@, \verb@VIRGO@,
\verb@GEO600@, \verb@TAMA300@, or \verb@CIT40@.  Additionally, a name
of \verb@-@ represents a fictitious detector with a fixed location in
the barycentric frame, and responding purely to the signal's +
polarization.

The response function file name should give the location of a readable
file relative to the current execution directory, which will be read
using the routine \verb@LALCReadFSeries()@; see \verb@StreamInput.h@
for discussion of the file format.  It should have a
\verb@sampleUnits@ field set to \verb@"strain count^-1"@, but
incorrect units will generate only a warning.  A response file name of
\verb@-@ specifies a unit response (i.e.\ raw strain is injected into
the data); note that this means that actual response files named
\verb@-@ are not permitted.

The input file name should give the location of a readable file
relative to the current execution directory, which will be read using
the routine \verb@LALSReadTimeSeries()@; see \verb@StreamInput.h@ for
discussion of the file format.  It should have a \verb@sampleUnits@
field set to \verb@"count"@, but incorrect units will generate only a
warning.  The input file name may be replaced with a \verb@-@ followed
by 4 whitespace-delimited numerical tokens
\verb@- epoch npt dt sigma@, indicating that a time series will be
generated starting at (\verb@INT8@) \verb@epoch@ GPS nanoseconds,
consisting of (\verb@UINT4@) \verb@npt@ points sampled at
(\verb@REAL8@) \verb@dt@ second intervals having Gaussian random noise
with rms value (\verb@REAL4@) \verb@sigma@ ADC counts.  Again, an
actual input file named \verb@-@ is therefore not permitted.

The output file name should give the location relative to the current
execution directory, where the injected data will be written using the
routine \verb@LALSWriteTSeries()@; see \verb@StreamInput.h@ for
discussion of the file format.  A name of \verb@-@ specifies writing
to standard output (\emph{not} to a file named \verb@-@).  To suppress
output, specify \verb@/dev/null@ as the output file.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define INJECTTESTC_ENORM 0
#define INJECTTESTC_ESUB  1
#define INJECTTESTC_EARG  2
#define INJECTTESTC_EVAL  3
#define INJECTTESTC_EFILE 4
#define INJECTTESTC_EMEM  5

#define INJECTTESTC_MSGENORM "Normal exit"
#define INJECTTESTC_MSGESUB  "Subroutine failed"
#define INJECTTESTC_MSGEARG  "Error parsing arguments"
#define INJECTTESTC_MSGEVAL  "Input argument out of valid range"
#define INJECTTESTC_MSGEFILE "Could not open file"
#define INJECTTESTC_MSGEMEM  "Out of memory"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

The program first reads \verb@sourcefile@ using
\verb@LALCHARReadVector()@, parses its arguments with
\verb@LALCreateTokenList()@, and generates a linked list of
\verb@CoherentGW@ structures containing the waveforms to be injected:
see below for a discussion of how the various generators' input lines
are parsed.  Blank lines are ignored, as are comments (a \verb@#@ or
\verb@%@ character causes the remainder of the line to be ignored).
The waveforms will be held in memory together for the rest of the
injection procedure; if this taxes the available system memory, break
up \verb@sourcefile@ and run the program on each separate file, using
the output from each run as the input for the next.

Once all the waveforms have been generated, the program starts reading
\verb@indxfile@: for each stretch of data specified, it reads in the
input file and response function file using \verb@LALSReadTSeries()@
and \verb@LALSReadFSeries()@, inverts the response function to get a
transfer function, simulates and injects the waveforms using
\verb@LALSimulateCoherentGW()@ and \verb@LALSSInjectTimeSeries()@, and
writes the resulting time series to the output file using
\verb@LALSWriteTSeries()@.  If a given stretch of data has the same
response function file name as the previous one, then the previous
transfer function is reused rather than re-reading and re-inverting
the data; thus, one should naturally arrange \verb@indxfile@ so that
data stretches using the same response function are listed
consecutively.

\paragraph{\texttt{LALGeneratePPNInspiral}:} Most arguments are parsed
from the corresponding tokens using the functions in
\verb@StringConvert.c@.  However, two input parameters require some
nontrivial computation.  First, the generator requires a start time,
while \verb@sourcefile@ specifies a coalescence time.  The solution is
to generate a waveform with an arbitrary start time (in this case zero
GPS seconds), take the returned coalescence time, and adjust the start
time accordingly.  Second, \verb@sourcefile@ does not specify a
sampling interval to be used in the generated waveform.  Here the
solution is to estimate an appropriate interval from the requested
termination frequency or from the point of post-Newtonian breakdown
for the specified masses.  If a nonzero termination frequency $f$ is
specified (either positive or negative), then the maximum rate of
change of frequency in the (post${}^0$-)Newtonian approximation is
$\dot{f}_\mathrm{max}=1.8\pi^{8/3}\tau^{5/3}f^{11/3}$, where
$\tau=Gm_\mathrm{tot}/c^3$ is the relativistic minimum timescale of
the system.  As explained in \verb@GeneratePPNInspiral.c@, for linear
interpolation of the waveform to be accurate to within $\pi/2$
radians, we would like $\Delta f\Delta t\lessim2$ over any sampling
interval.  This implies that $\Delta
t\approx\sqrt{2/\dot{f}_\mathrm{max}}$, or:
$$
\Delta t \approx 0.1403\tau^{-5/6}f^{-11/6} \;.
$$
If the maximum frequency is given as zero (i.e.\ unspecified) or
positive, then an additional constraint is imposed by the guaranteed
post-Newtonian breakdown at or before $r=2Gm_\mathrm{tot}/c^2$, giving
$\dot{f}_\mathrm{max}=(2\sqrt{2}/40\pi)\tau^{-2}$.  This implies that:
$$
\Delta t_\mathrm{min} \approx 7.697\tau \;,
$$
and we can write the previous expression as $\Delta t = \max\{\Delta
t_\mathrm{min}, 0.7685(\Delta t_\mathrm{min})^{-5/6}f^{-11/6}\}$.
When the waveform is actually generated, the maximum value of $\Delta
f\Delta t$ is returned, and, if greater than 2, $\Delta t$ is reduced
accordingly and the waveform regenerated.

\paragraph{\texttt{LALGenerateTaylorCW}:} Most arguments are parsed
from the corresponding tokens using the routines in
\verb@StringConvert.c@.  However, two input parameters require
computation.  First, the base frequency $f_0$ and spindown terms
$f_1,\ldots f_N$ need to be transformed from the reference time to the
start time, by the following formula:
\begin{eqnarray}
f_{0\mathrm{(start)}} & = & f_0\left( 1+\sum_{k=1}^N f_k t^k \right)
	\;,\nonumber\\
f_{j\mathrm{(start)}} & = & \frac{\sum_{k=j}^N{j \choose k}f_k t^{k-j}}
	{1+\sum_{k=1}^N f_k t^k} \;,\nonumber
\end{eqnarray}
where $t=t_\mathrm{start}-t_\mathrm{ref}$ is the time shift.  These
calculations are done to double precision; even so, one is advised to
specify the frequency and spindown terms at a reference time that is
not too far from the times being considered.  Roughly speaking, for
phases to be accurate within a fraction of a cycle, this means that
$1+\sum_k|f_k\tau^k|\ll10^{15}/|f_0\tau|$, where $\tau$ is the larger
of $t_\mathrm{start}-t_\mathrm{ref}$ and
$t_\mathrm{stop}-t_\mathrm{ref}$.

Second, the program must choose an appropriate sampling interval for
the waveforms, such that the maximum $\Delta f\Delta t$ over any
interval is less than 2.  This is fairly straightforward; one simply
takes:
$$
\Delta t \approx \left( 0.5 f_0 \sum_{k=1}^N k|f_k T^{k-1}|
	\right)^{-1/2} \;,
$$
where $T=t_\mathrm{stop}-t_\mathrm{start}$ is the total signal length.
If this gives $\Delta t>T$ (as for instance when there are no nonzero
spindown terms), then $\Delta t$ is set equal to $T$.

\paragraph{\texttt{LALGenerateSpinOrbitCW}:} Most arguments are parsed
from the corresponding tokens using the routines in
\verb@StringConvert.c@.  The reference time is assumed to correspond
to the orbital epoch of the system, and the conversion from reference
time to start time discussed above is performed automatically within
the \verb@LALGenerateSpinOrbitCW()@ function.  However, the program
still needs to choose an appropriate sampling interval for the
waveforms, such that the maximum $\Delta f\Delta t$ over any interval
is less than 2.  We simply take the formula for $\Delta t$ above and
add to the spindown terms a Doppler modulation at least as large as
the maximum possible orbital acceleration, to get:
$$
\Delta t \approx \left( 0.5 f_0 \left[ \dot{\upsilon}_p^2 r_p\sin i/c +
	\sum_{k=1}^N k|f_k T^{k-1}| \right] \right)^{-1/2} \;,
$$
where $\dot{\upsilon}_p$ is the angular speed at periapsis, $r_p\sin
i/c$ is the projected, normalized periapsis distance, and
$T=t_\mathrm{stop}-t_\mathrm{start}$ is the total signal length.  If
this gives $\Delta t>T$ then $\Delta t$ is set equal to $T$.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALMalloc()                     LALFree()
LALCreateRandomParams()         LALDestroyRandomParams()
LALCreateTokenList()            LALDestroyTokenList()
LALCHARCreateVector()           LALCHARDestroyVector()
LALCCreateVector()              LALCDestroyVector()
LALDCreateVector()              LALDDestroyVector()
LALSDestroyVector()             LALSDestroyVectorSequence()
LALCReadFSeries()               LALSReadTSeries()
LALCHARReadVector()             LALSWriteTSeries()
LALCCVectorDivide()             LALNormalDeviates()
LALUnitRaise()                  LALUnitMultiply()
LALUnitCompare()                LALStringToI8()
LALStringToS()                  LALStringToD()
LALGeneratePPNInspiral()        LALGenerateTaylorCW()
LALGenerateSpinOrbitCW()        LALSimulateCoherentGW()
LALSSInjectTimeSeries()         LALInitBarycenter()
snprintf()                   LALPrintError()
LALCheckMemoryLeaks()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{InjectTestCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>
#include <lal/Sort.h>
#include <lal/Random.h>
#include <lal/VectorOps.h>
#include <lal/StringInput.h>
#include <lal/StreamInput.h>
#include <lal/StreamOutput.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateTaylorCW.h>
#include <lal/GenerateSpinOrbitCW.h>

NRCSID( INJECTTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;

#define MSGLEN 1024  /* Max. length of warning/info messages */
#define FMAX (1.0e9) /* Max. frequency we will ever use (Hz) */
#define WHITESPACE " \f\n\r\t\v" /* list of whitespace characters */

/* Usage format string. */
#define USAGE "Usage: %s [-s sourcefile] [-e earthfile sunfile]\n"   \
"\t[-r randomseed] [-d debuglevel] [-h] indxfile\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), program, __FILE__,       \
		 __LINE__, INJECTTESTC, statement ? statement :      \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
		 "        %s\n", program, __FILE__, __LINE__,        \
		 INJECTTESTC, (statement) );                         \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", program, __FILE__, __LINE__,        \
		 INJECTTESTC, (statement) );                         \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( INJECTTESTC_ESUB, INJECTTESTC_MSGESUB,                      \
         "Function call \"" #func "\" failed:" );                    \
  exit( INJECTTESTC_ESUB );                                          \
}                                                                    \
while (0)

#define CHECKVAL( val, lower, upper )                                \
do                                                                   \
if ( ( (val) <= (lower) ) || ( (val) > (upper) ) )                   \
{                                                                    \
  ERROR( INJECTTESTC_EVAL, INJECTTESTC_MSGEVAL,                      \
         "Value of " #val " out of range:" );                        \
  LALPrintError( #val " = %f, range = (%f,%f]\n", (REAL8)(val),      \
                 (REAL8)(lower), (REAL8)(upper) );                   \
  exit( INJECTTESTC_EVAL );                                          \
}                                                                    \
while (0)

/* A macro for parsing a numerical token in a string.  The variables
   token, endptr, and ok must already exist within the scope of the
   macro call. */
#define PARSETOKEN( parser, variable )                               \
do {                                                                 \
  if ( ok ) {                                                        \
    if ( *(++token) ) {                                              \
      CHAR *endptr;                                                  \
      SUB( (parser)( stat, (variable), *token, &endptr ), stat );    \
      ok = ( *endptr == '\0' );                                      \
    } else                                                           \
      ok = 0;                                                        \
  }                                                                  \
} while (0)


/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

/* Some more global constants (set once and used at all levels). */
char *program;
LALUnit countsPerStrainUnit, strainPerCountUnit;

/* A linked list of CoherentGW structures. */
typedef struct tagCoherentGWNode {
  CoherentGW waveform;
  struct tagCoherentGWNode *next;
} CoherentGWNode;


/* A function to convert INT8 nanoseconds to LIGOTimeGPS. */
void
I8ToLIGOTimeGPS( LIGOTimeGPS *output, INT8 input );

/* Function to parse source lines to generate various waveform. */
INT4
GenerateWaveform( LALStatus      *stat,
		  CoherentGWNode *head,
		  TokenList      *tokens,
		  const CHAR     *tag,
		  RandomParams   *randParams );

/* Function to inject sources into each input file. */
INT4
InjectWaveforms( LALStatus        *stat,
		 CoherentGWNode   *head,
		 TokenList        *tokens,
		 const CHAR       *tag,
		 RandomParams     *randParams,
		 DetectorResponse *detector,
		 CHARVector       **prevResponse );

/* Function to compute the combination of two numbers. */
UINT4
choose( UINT4 a, UINT4 b );


int
main(int argc, char **argv)
{
  /* Command-line parsing variables. */
  int arg;                         /* command-line argument counter */
  static LALStatus stat;           /* status structure */
  CHAR *sourcefile = NULL;         /* name of source file */
  CHAR *earthfile = NULL;          /* name of earth ephemeris file */
  CHAR *sunfile = NULL;            /* name of sun ephemeris file */
  CHAR *indxfile = NULL;           /* name of index file */
  INT4 seed = 0;                   /* random number seed */

  CoherentGWNode head;             /* head of list of waveforms */
  CHAR msg[MSGLEN];                /* warning or info messages */
  DetectorResponse detector;       /* detector site geometry */
  CHARVector *prevResponse = NULL; /* previous response filename */
  RandomParams *randParams = NULL; /* random parameters */

  UINT4 i;                         /* current position in line */
  UINT4 lineno;                    /* current line number */
  UINT4 seriesno;                  /* number of time series written */
  FILE *fp;                        /* generic file pointer */

  /*******************************************************************
   *                                                                 *
   * Parse command-line arguments.  arg stores the current position. *
   *                                                                 *
   *******************************************************************/

  program = *argv;
  arg = 1;
  while ( arg < argc ) {
    /* Parse source file option. */
    if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	sourcefile = argv[arg++];
      } else {
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Parse ephemeris file option. */
    else if ( !strcmp( argv[arg], "-e" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	earthfile = argv[arg++];
	sunfile = argv[arg++];
      } else {
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Parse random seed option. */
    else if ( !strcmp( argv[arg], "-r" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	seed = atoi( argv[arg++] );
      } else {
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      } else {
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Parse help option. */
    else if ( !strcmp( argv[arg], "-h" ) ) {
      LALPrintError( USAGE, *argv );
      INFO( INJECTTESTC_MSGENORM );
      return INJECTTESTC_ENORM;
    }
    /* Get index file name. */
    else {
      if ( argc == arg + 1 ) {
	indxfile = argv[arg++];
      } else {
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
  } /* End of argument parsing loop. */

  /* Make sure that the index file was specified.  If it isn't, print
     a warning but exit normally (since all test programs should
     return normally when run without arguments). */
  if ( !indxfile ) {
    WARNING( "No index file given (exiting)" );
    INFO( INJECTTESTC_MSGENORM );
    return INJECTTESTC_ENORM;
  }

  /* Initialize random parameters here, in case any future waveform
     generators will require it (none do at present). */
  SUB( LALCreateRandomParams( &stat, &randParams, seed ), &stat );


  /*******************************************************************
   *                                                                 *
   * Read sourcefile, and create waveforms.                          *
   *                                                                 *
   *******************************************************************/

  /* Set up list and start reading file. */
  memset( &head, 0, sizeof(CoherentGWNode) );
  if ( sourcefile ) {
    UINT4 sourceno = 0; /* source number */
    if ( ( fp = fopen( sourcefile, "r" ) ) == NULL ) {
      ERROR( INJECTTESTC_EFILE, "- " INJECTTESTC_MSGEFILE, sourcefile );
      return INJECTTESTC_EFILE;
    }

    /* Read a line, skipping over comments. */
    lineno = 0;
    while ( !feof( fp ) ) {
      CHARVector *line = NULL;  /* input line */
      TokenList *tokens = NULL; /* input line parsed into tokens */
      SUB( LALCHARReadVector( &stat, &line, fp ), &stat );
      lineno++;
      for ( i = 0; i < line->length; i++ )
	if ( line->data[i] == '%' || line->data[i] == '#' ) {
	  line->data[i] = '\0';
	  i = line->length;
	}
      SUB( LALCreateTokenList( &stat, &tokens, line->data,
			       WHITESPACE ), &stat );
      SUB( LALCHARDestroyVector( &stat, &line ), &stat );

      /* If line isn't blank, generate a waveform with it. */
      if ( tokens->nTokens > 0 ) {
	snprintf( msg, MSGLEN, "%s line %i:", sourcefile, lineno );
	if ( !GenerateWaveform( &stat, &head, tokens, msg,
				randParams ) )
	  sourceno++;
      }

      /* Done parsing line. */
      SUB( LALDestroyTokenList( &stat, &tokens ), &stat );
    }

    /* Done parsing sourcefile. */
    fclose( fp );
    snprintf( msg, MSGLEN, "%s: Parsed %i lines, generated %i"
		 " sources", sourcefile, lineno, sourceno );
    INFO( msg );
  }


  /*******************************************************************
   *                                                                 *
   * Read index file, and perform injections.                        *
   *                                                                 *
   *******************************************************************/

  /* First, set up ephemeris data. */
  memset( &detector, 0, sizeof(DetectorResponse) );
  if ( earthfile && sunfile ) {
    if ( !( detector.ephemerides = (EphemerisData *)
	    LALMalloc( sizeof(EphemerisData) ) ) ) {
      ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
      return INJECTTESTC_EMEM;
    }
    memset( detector.ephemerides, 0, sizeof(EphemerisData) );
    detector.ephemerides->ephiles.earthEphemeris = earthfile;
    detector.ephemerides->ephiles.sunEphemeris = sunfile;
    SUB( LALInitBarycenter( &stat, detector.ephemerides ), &stat );
  }

  /* Set up (dummy) initial response function and filename, and
     transfer function units. */
  if ( !( detector.transfer = (COMPLEX8FrequencySeries *)
	  LALMalloc( sizeof(COMPLEX8FrequencySeries) ) ) ) {
    ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
    return INJECTTESTC_EMEM;
  }
  memset( detector.transfer, 0, sizeof(COMPLEX8FrequencySeries) );
  SUB( LALCCreateVector( &stat, &(detector.transfer->data), 1 ),
       &stat );
  {
    RAT4 negOne;
    LALUnit unit;
    LALUnitPair pair;
    negOne.numerator = -1;
    negOne.denominatorMinusOne = 0;
    SUB( LALUnitRaise( &stat, &unit, &lalStrainUnit, &negOne ),
	 &stat );
    pair.unitOne = &lalADCCountUnit;
    pair.unitTwo = &unit;
    SUB( LALUnitMultiply( &stat, &countsPerStrainUnit, &pair ),
	 &stat );
    SUB( LALUnitRaise( &stat, &strainPerCountUnit,
		       &countsPerStrainUnit, &negOne ), &stat );
  }
  SUB( LALCHARCreateVector( &stat, &prevResponse, 1 ), &stat );
  prevResponse->data[0] = '\0';


  /* Now start reading index file. */
  if ( ( fp = fopen( indxfile, "r" ) ) == NULL ) {
    ERROR( INJECTTESTC_EFILE, "- " INJECTTESTC_MSGEFILE, indxfile );
    return INJECTTESTC_EFILE;
  }

  /* Read a line, skipping over comments. */
  lineno = 0;
  seriesno = 0;
  while ( !feof( fp ) ) {
    CHARVector *line = NULL;  /* input line */
    TokenList *tokens = NULL; /* input line parsed into tokens */

    SUB( LALCHARReadVector( &stat, &line, fp ), &stat );
    lineno++;
    for ( i = 0; i < line->length; i++ )
      if ( line->data[i] == '%' || line->data[i] == '#' ) {
	line->data[i] = '\0';
	i = line->length;
      }
    SUB( LALCreateTokenList( &stat, &tokens, line->data, WHITESPACE ),
	 &stat );
    SUB( LALCHARDestroyVector( &stat, &line ), &stat );

    /* If line isn't blank, use it to perform injections. */
    if ( tokens->nTokens ) {
      snprintf( msg, MSGLEN, "%s line %i:", indxfile, lineno );
      if ( !InjectWaveforms( &stat, &head, tokens, msg, randParams,
			     &detector, &prevResponse ) )
	seriesno++;
    }
    SUB( LALDestroyTokenList( &stat, &tokens ), &stat );
  }

  fclose( fp );
  snprintf( msg, MSGLEN, "%s: Parsed %i lines, wrote %i time"
	       " series out", indxfile, lineno, seriesno );
  INFO( msg );


  /*******************************************************************
   *                                                                 *
   * Clean up.                                                       *
   *                                                                 *
   *******************************************************************/

  /* Clean up memory. */
  SUB( LALCHARDestroyVector( &stat, &prevResponse ), &stat );
  SUB( LALCDestroyVector( &stat, &(detector.transfer->data) ),
       &stat );
  LALFree( detector.transfer );
  if ( detector.ephemerides )
    LALFree( detector.ephemerides );
  SUB( LALDestroyRandomParams( &stat, &randParams ), &stat );
  {
    CoherentGWNode *here = head.next;
    while ( here ) {
      CoherentGWNode *last = here;
      here = here->next;
      SUB( LALSDestroyVectorSequence( &stat, &(last->waveform.a->data) ),
	   &stat );
      SUB( LALSDestroyVector( &stat, &(last->waveform.f->data) ),
	   &stat );
      SUB( LALDDestroyVector( &stat, &(last->waveform.phi->data) ),
	   &stat );
      LALFree( last->waveform.a );
      LALFree( last->waveform.f );
      LALFree( last->waveform.phi );
      LALFree( last );
    }
  }

  /* Exit. */
  LALCheckMemoryLeaks();
  INFO( INJECTTESTC_MSGENORM );
  return INJECTTESTC_ENORM;
}


/* A function to convert INT8 nanoseconds to LIGOTimeGPS. */
void
I8ToLIGOTimeGPS( LIGOTimeGPS *output, INT8 input )
{
  INT8 s = input / 1000000000LL;
  output->gpsSeconds = (INT4)( s );
  output->gpsNanoSeconds = (INT4)( input - 1000000000LL*s );
  return;
}


/* Function to parse source lines to generate various waveform. */
INT4
GenerateWaveform( LALStatus      *stat,
		  CoherentGWNode *head,
		  TokenList      *tokens,
		  const CHAR     *tag,
		  RandomParams   *randParams )
{
  CoherentGWNode *here = head;   /* current waveform in list */
  CHAR msg[MSGLEN];              /* warning messages */
  BOOLEAN ok = 1;                /* whether current line is okay */
  CHAR **token = tokens->tokens; /* handle to current token */

  /* At present none of the generators use random data.  This line is
     to keep gcc from complaining about unused variables. */
  randParams = randParams;

  /* Jump to last source in linked list. */
  while ( here->next )
    here = here->next;

  /* Generate PPN inspiral waveform. */
  if ( !strcmp( *token, "LALGeneratePPNInspiral" ) ) {
    INT8 t1, t2;  /* times to coalescence (nanoseconds) */
    REAL8 dt;     /* sample interval (s) */
    REAL4 m1, m2; /* component binary masses (solar) */
    PPNParamStruc params; /* PPN inspiral parameters */
    memset( &params, 0, sizeof(PPNParamStruc) );

    if ( tokens->nTokens < 12 ) {
      snprintf( msg, MSGLEN, "%s Need at least 11 arguments for"
		   " generator %s, got %i (skipping line)", tag,
		   *token, tokens->nTokens - 1 );
      WARNING( msg );
      return 1;
    }

    /* Parse input line. */
    PARSETOKEN( LALStringToI8, &t1 );
    PARSETOKEN( LALStringToS, &m1 );
    PARSETOKEN( LALStringToS, &m2 );
    PARSETOKEN( LALStringToS, &(params.d) );
    PARSETOKEN( LALStringToS, &(params.inc) );
    PARSETOKEN( LALStringToD, &(params.position.longitude) );
    PARSETOKEN( LALStringToD, &(params.position.latitude) );
    PARSETOKEN( LALStringToS, &(params.psi) );
    PARSETOKEN( LALStringToS, &(params.phi) );
    PARSETOKEN( LALStringToS, &(params.fStartIn) );
    PARSETOKEN( LALStringToS, &(params.fStopIn) );
    if ( !ok ) {
      snprintf( msg, MSGLEN, "%s Error parsing argument %i"
		   " (skipping line)", tag,
		   (UINT4)( token - tokens->tokens ) );
      WARNING( msg );
      return 1;
    }

    /* Perform necessary input conversions. */
    params.mTot = m1 + m2;
    params.eta = m1*m2/params.mTot/params.mTot;
    params.d *= 1000000.0*LAL_PC_SI;
    params.psi *= LAL_PI/180.0;
    params.inc *= LAL_PI/180.0;
    params.phi *= LAL_PI/180.0;
    params.position.longitude *= LAL_PI/180.0;
    params.position.latitude *= LAL_PI/180.0;
    params.position.system = COORDINATESYSTEM_EQUATORIAL;

    /* Estimate sampling interval. */
    dt = 7.697*params.mTot / LAL_MTSUN_SI;
    if ( params.fStop != 0.0 )
      params.deltaT = 0.7685*pow( dt, -5.0/6.0 )
	*pow( fabs( params.fStop ), -11.0/6.0 );
    if ( params.fStop < 0.0 || params.deltaT < dt )
      params.deltaT = dt;

    /* Generate waveform.  Repeat if dfdt is too large. */
    if ( !( here->next = (CoherentGWNode *)
	    LALMalloc( sizeof(CoherentGWNode) ) ) ) {
      ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
      exit( INJECTTESTC_EMEM );
    }
    here = here->next;
    do {
      memset( here, 0, sizeof(CoherentGWNode) );
      SUB( LALGeneratePPNInspiral( stat, &(here->waveform), &params ),
	   stat );
      if ( params.dfdt > 2.0 ) {
	dt = sqrt( 2.0/params.dfdt );
	if ( dt > 0.5 )
	  dt = 0.5;
	params.deltaT *= dt;
	SUB( LALSDestroyVectorSequence( stat, &(here->waveform.a->data) ),
	     stat );
	SUB( LALSDestroyVector( stat, &(here->waveform.f->data) ),
	     stat );
	SUB( LALDDestroyVector( stat, &(here->waveform.phi->data) ),
	     stat );
	LALFree( here->waveform.a );
	LALFree( here->waveform.f );
	LALFree( here->waveform.phi );
      }
    } while ( params.dfdt > 2.0 );

    /* Match the coalescence time. */
    t2 = (INT8)( params.tc );
    params.tc -= t2;
    t2 *= 1000000000LL;
    t2 += (INT8)( ( 1.0e9 )*params.tc );
    I8ToLIGOTimeGPS( &(here->waveform.a->epoch), t1 - t2 );
    here->waveform.f->epoch = here->waveform.phi->epoch
      = here->waveform.a->epoch;
  }


  /* Generate Taylor-parametrized continuous waveform. */
  else if ( !strcmp( *token, "LALGenerateTaylorCW" ) ) {
    UINT4 i, nSpin;   /* index over and number of spindown terms */
    INT8 t0, t1, t2;  /* reference, start, and stop times (ns) */
    REAL8 t, dt;      /* series length and sample interval (s) */
    REAL8 dfFac;      /* estimated maximum fdot/f0 (Hz) */
    TaylorCWParamStruc params; /* Taylor CW parameters */
    memset( &params, 0, sizeof(PPNParamStruc) );

    if ( tokens->nTokens < 11 ) {
      snprintf( msg, MSGLEN, "%s Need at least 10 arguments for"
		   " generator %s, got %i (skipping line)", tag,
		   *token, tokens->nTokens - 1 );
      WARNING( msg );
      return 1;
    }

    /* Allocate spindown parameters. */
    nSpin = tokens->nTokens - 11;
    if ( nSpin ) {
      SUB( LALDCreateVector( stat, &(params.f), nSpin ), stat );
    }

    /* Parse input line. */
    PARSETOKEN( LALStringToI8, &t0 );
    PARSETOKEN( LALStringToI8, &t1 );
    PARSETOKEN( LALStringToI8, &t2 );
    PARSETOKEN( LALStringToS, &(params.aPlus) );
    PARSETOKEN( LALStringToS, &(params.aCross) );
    PARSETOKEN( LALStringToD, &(params.position.longitude) );
    PARSETOKEN( LALStringToD, &(params.position.latitude) );
    PARSETOKEN( LALStringToS, &(params.psi) );
    PARSETOKEN( LALStringToD, &(params.phi0) );
    PARSETOKEN( LALStringToD, &(params.f0) );
    for ( i = 0; i < nSpin; i++ ) {
      PARSETOKEN( LALStringToD, params.f->data + i );
    }
    if ( !ok ) {
      snprintf( msg, MSGLEN, "%s Error parsing argument %i"
		   " (skipping line)", tag,
		   (UINT4)( token - tokens->tokens ) );
      WARNING( msg );
      return 1;
    }
    if ( t2 <= t1 ) {
      snprintf( msg, MSGLEN, "%s Stop time precedes start time"
		   " (skipping line)", tag );
      WARNING( msg );
      return 1;
    }

    /* Perform necessary input conversions. */
    params.psi *= LAL_PI/180.0;
    params.phi0 *= LAL_PI/180.0;
    params.position.longitude *= LAL_PI/180.0;
    params.position.latitude *= LAL_PI/180.0;
    params.position.system = COORDINATESYSTEM_EQUATORIAL;

    /* Shift reference time to start of waveform. */
    I8ToLIGOTimeGPS( &(params.epoch), t1 );
    t = (1.0e-9)*(REAL8)( t1 - t0 );
    if ( nSpin ) {
      REAL8 tN = 1.0;   /* t raised to various powers */
      REAL8 fFac = 1.0; /* fractional change in frequency */
      REAL8 tFac = 1.0; /* time integral of fFac */
      REAL8 *fData = params.f->data; /* spindown coeficients */
      for ( i = 0; i < nSpin; i++ ) {
	UINT4 j;        /* index over remaining spindown terms */
	REAL8 tM = 1.0; /* t raised to various powers */
	fFac += fData[i]*( tN *= t );
	tFac += fData[i]*tN/( i + 2.0 );
	for ( j = i + 1; j < nSpin; j++ )
	  fData[i] += choose( j + 1, i + 1 )*fData[j]*( tM *= t );
      }
      params.phi0 += LAL_TWOPI*params.f0*t*tFac;
      params.f0 *= fFac;
      for ( i = 0; i < nSpin; i++ )
	fData[i] /= fFac;
    }

    /* Estimate sampling interval. */
    t = params.deltaT = (1.0e-9)*(REAL8)( t2 - t1 );
    params.length = 2;
    if ( nSpin ) {
      REAL8 tN = 1.0;    /* t raised to various powers */
      REAL8 *fData = params.f->data; /* spindown coeficients */
      dfFac = 0.0;
      for ( i = 0; i < nSpin; i++ )
	dfFac += i*fabs( fData[i]*( tN *= t ) );
      if ( dfFac > 0.0 ) {
	dt = sqrt( 2.0/( params.f0*dfFac ) );
	if ( dt < params.deltaT ) {
	  params.deltaT = dt;
	  params.length = (UINT4)( t/params.deltaT ) + 1;
	}
      }
    }

    /* Generate waveform.  Repeat if dfdt is too large. */
    if ( !( here->next = (CoherentGWNode *)
	    LALMalloc( sizeof(CoherentGWNode) ) ) ) {
      ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
      exit( INJECTTESTC_EMEM );
    }
    here = here->next;
    do {
      memset( here, 0, sizeof(CoherentGWNode) );
      SUB( LALGenerateTaylorCW( stat, &(here->waveform), &params ),
	   stat );
      if ( params.dfdt > 2.0 ) {
	dt = sqrt( 2.0/params.dfdt );
	if ( dt > 0.5 )
	  dt = 0.5;
	params.deltaT *= dt;
	params.length = (UINT4)( t/params.deltaT ) + 1;
	SUB( LALSDestroyVectorSequence( stat, &(here->waveform.a->data) ),
	     stat );
	SUB( LALSDestroyVector( stat, &(here->waveform.f->data) ),
	     stat );
	SUB( LALDDestroyVector( stat, &(here->waveform.phi->data) ),
	     stat );
	LALFree( here->waveform.a );
	LALFree( here->waveform.f );
	LALFree( here->waveform.phi );
      }
    } while ( params.dfdt > 2.0 );

    /* Deallocate spindown parameters. */
    if ( nSpin ) {
      SUB( LALDDestroyVector( stat, &(params.f) ), stat );
    }
  }


  /* Generate Taylor-parametrized and Doppler-modulated
     continuous waveform. */
  else if ( !strcmp( *token, "LALGenerateSpinOrbitCW" ) ) {
    UINT4 i, nSpin;   /* index over and number of spindown terms */
    INT8 t0, t1, t2;  /* reference, start, and stop times (ns) */
    REAL8 t, dt;      /* series length and sample interval (s) */
    REAL8 dfFac;      /* estimated maximum fdot/f0 (Hz). */
    SpinOrbitCWParamStruc params; /* Taylor CW parameters */
    memset( &params, 0, sizeof(PPNParamStruc) );

    if ( tokens->nTokens < 15 ) {
      snprintf( msg, MSGLEN, "%s Need at least 14 arguments for"
		   " generator %s, got %i (skipping line)", tag,
		   *token, tokens->nTokens - 1 );
      WARNING( msg );
      return 1;
    }

    /* Allocate spindown parameters. */
    nSpin = tokens->nTokens - 15;
    if ( nSpin ) {
      SUB( LALDCreateVector( stat, &(params.f), nSpin ), stat );
    }

    /* Parse input line. */
    PARSETOKEN( LALStringToI8, &t0 );
    PARSETOKEN( LALStringToI8, &t1 );
    PARSETOKEN( LALStringToI8, &t2 );
    PARSETOKEN( LALStringToS, &(params.aPlus) );
    PARSETOKEN( LALStringToS, &(params.aCross) );
    PARSETOKEN( LALStringToD, &(params.position.longitude) );
    PARSETOKEN( LALStringToD, &(params.position.latitude) );
    PARSETOKEN( LALStringToS, &(params.psi) );
    PARSETOKEN( LALStringToD, &(params.phi0) );
    PARSETOKEN( LALStringToD, &(params.f0) );
    for ( i = 0; i < nSpin; i++ ) {
      PARSETOKEN( LALStringToD, params.f->data + i );
    }
    PARSETOKEN( LALStringToD, &(params.omega) );
    PARSETOKEN( LALStringToD, &(params.angularSpeed) );
    PARSETOKEN( LALStringToD, &(params.rPeriNorm) );
    PARSETOKEN( LALStringToD, &(params.oneMinusEcc) );
    if ( !ok ) {
      snprintf( msg, MSGLEN, "%s Error parsing argument %i"
		   " (skipping line)", tag,
		   (UINT4)( token - tokens->tokens ) );
      WARNING( msg );
      return 1;
    }
    if ( t2 <= t1 ) {
      snprintf( msg, MSGLEN, "%s Stop time precedes start time"
		   " (skipping line)", tag );
      WARNING( msg );
      return 1;
    }

    /* Perform necessary input conversions. */
    params.psi *= LAL_PI/180.0;
    params.phi0 *= LAL_PI/180.0;
    params.position.longitude *= LAL_PI/180.0;
    params.position.latitude *= LAL_PI/180.0;
    params.position.system = COORDINATESYSTEM_EQUATORIAL;
    params.oneMinusEcc = 1.0 - params.oneMinusEcc;
    I8ToLIGOTimeGPS( &(params.epoch), t1 );
    I8ToLIGOTimeGPS( &(params.spinEpoch), t0 );
    params.orbitEpoch = params.spinEpoch;

    /* Estimate sampling interval. */
    t = params.deltaT = (1.0e-9)*(REAL8)( t2 - t1 );
    params.length = 2;
    dfFac = fabs( params.angularSpeed*params.angularSpeed*
		  params.rPeriNorm );
    if ( nSpin ) {
      REAL8 tN = 1.0;    /* t raised to various powers */
      REAL8 *fData = params.f->data; /* spindown coeficients */
      for ( i = 0; i < nSpin; i++ )
	dfFac += i*fabs( fData[i]*( tN *= t ) );
    }
    if ( dfFac > 0.0 ) {
      dt = sqrt( 2.0/( params.f0*dfFac ) );
      if ( dt < params.deltaT ) {
	params.deltaT = dt;
	params.length = (UINT4)( t/params.deltaT ) + 1;
      }
    }

    /* Generate waveform.  Repeat if dfdt is too large. */
    if ( !( here->next = (CoherentGWNode *)
	    LALMalloc( sizeof(CoherentGWNode) ) ) ) {
      ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
      exit( INJECTTESTC_EMEM );
    }
    here = here->next;
    do {
      memset( here, 0, sizeof(CoherentGWNode) );
      SUB( LALGenerateSpinOrbitCW( stat, &(here->waveform), &params ),
	   stat );
      if ( params.dfdt > 2.0 ) {
	dt = sqrt( 2.0/params.dfdt );
	if ( dt > 0.5 )
	  dt = 0.5;
	params.deltaT *= dt;
	params.length = (UINT4)( t/params.deltaT ) + 1;
	SUB( LALSDestroyVectorSequence( stat, &(here->waveform.a->data) ),
	     stat );
	SUB( LALSDestroyVector( stat, &(here->waveform.f->data) ),
	     stat );
	SUB( LALDDestroyVector( stat, &(here->waveform.phi->data) ),
	     stat );
	LALFree( here->waveform.a );
	LALFree( here->waveform.f );
	LALFree( here->waveform.phi );
      }
    } while ( params.dfdt > 2.0 );

    /* Deallocate spindown parameters. */
    if ( nSpin ) {
      SUB( LALDDestroyVector( stat, &(params.f) ), stat );
    }
  }


  /* Generator method does not match one of the known types. */
  else {
    snprintf( msg, MSGLEN, "%s Unknown generator %s"
		 " (skipping line)", tag, *token );
    WARNING( msg );
    return 1;
  }

  return 0;
}


/* Function to inject sources into each input file. */
INT4
InjectWaveforms( LALStatus        *stat,
		 CoherentGWNode   *head,
		 TokenList        *tokens,
		 const CHAR       *tag,
		 RandomParams     *randParams,
		 DetectorResponse *detector,
		 CHARVector       **prevResponse )
{
  UINT4 sourceno = 0;       /* number of sources injected */
  CHAR name[LALNameLength]; /* name of time series */
  CHAR msg[MSGLEN];         /* warning messages */
  CHAR **token;             /* handle to current token */
  BOOLEAN ok;               /* whether current line is okay */
  CoherentGWNode *here;     /* current source in list */
  REAL4TimeSeries input;    /* input data */
  static LALDetector site;
  memset( &input, 0, sizeof(REAL4TimeSeries) );

  /* Line must have 4 arguments iff input is read from a file, and
     8 arguments iff input is generated randomly. */
  token = tokens->tokens + 2;
  ok = ( tokens->nTokens == 4 && strcmp( *token, "-" ) )
    || ( tokens->nTokens == 8 && !strcmp( *token, "-" ) );
  if ( !ok ) {
    snprintf( msg, MSGLEN, "%s Unreognized number of tokens"
		 " (skipping line)", tag );
    WARNING( msg );
    return 1;
  }

  /* Set detector location. */
  token = tokens->tokens;
  detector->site = &site;
  if ( !strcmp( *token, "-" ) )
    detector->site = NULL;
  else if ( !strcmp( *token, "LHO" ) )
    site = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else if ( !strcmp( *token, "LLO" ) )
    site = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if ( !strcmp( *token, "VIRGO" ) )
    site = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  else if ( !strcmp( *token, "GEO600" ) )
    site = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  else if ( !strcmp( *token, "TAMA300" ) )
    site = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  else if ( !strcmp( *token, "CIT40" ) )
    site = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
  else {
    snprintf( msg, MSGLEN, "%s Unreognized detector %s"
		 " (skipping line)", tag, *token );
    WARNING( msg );
    return 1;
  }

  /* Get transfer function, if necessary. */
  token = tokens->tokens + 1;
  if ( strcmp( *token, (*prevResponse)->data ) ) {
    COMPLEX8FrequencySeries transfer; /* transfer function */
    memset( &transfer, 0, sizeof(COMPLEX8FrequencySeries) );

    /* Check if a unit response was requested. */
    if ( !strcmp( *token, "-" ) ) {
      transfer.f0 = 0.0;
      transfer.deltaF = FMAX;
      SUB( LALCCreateVector( stat, &(transfer.data), 2 ), stat );
      transfer.data->data[0].re = transfer.data->data[1].re = 1.0;
      transfer.data->data[0].im = transfer.data->data[1].im = 0.0;
      transfer.sampleUnits = countsPerStrainUnit;
    }

    /* Otherwise, read response from a file. */
    else {
      BOOLEAN unitsOK;                 /* whether read units are right */
      UINT4 i;                         /* generic index over data */
      LALUnitPair pair;                /* input for unit comparison */
      COMPLEX8Vector *unity = NULL;    /* complex 1's */
      FILE *fp = fopen( *token, "r" ); /* response file */
      if ( !fp ) {
	snprintf( msg, MSGLEN, "%s Response file %s not found"
		     " (skipping line)", tag, *token );
	WARNING( msg );
	return 1;
      }
      SUB( LALCReadFSeries( stat, &transfer, fp ), stat );
      fclose( fp );
      if ( transfer.data == NULL ) {
	snprintf( msg, MSGLEN, "%s Response file %s has no data"
		     " (skipping line)", tag, *token );
	WARNING( msg );
	return 1;
      }
      if ( transfer.deltaF <= 0.0 ) {
	snprintf( msg, MSGLEN, "%s Response file %s has bad"
		     " metadata (skipping line)", tag, *token );
	WARNING( msg );
	SUB( LALCDestroyVector( stat, &(transfer.data) ), stat );
	return 1;
      }

      /* Invert response function to get transfer function. */
      SUB( LALCCreateVector( stat, &unity, transfer.data->length ),
	   stat );
      for ( i = 0; i < unity->length; i++ ) {
	unity->data[i].re = 1.0;
	unity->data[i].im = 0.0;
      }
      SUB( LALCCVectorDivide( stat, transfer.data, unity,
			      transfer.data ), stat );
      SUB( LALCDestroyVector( stat, &unity ), stat );

      /* Check input file units. */
      pair.unitOne = &strainPerCountUnit;
      pair.unitOne = &(transfer.sampleUnits);
      SUB( LALUnitCompare( stat, &unitsOK, &pair ), stat );
      if ( !unitsOK ) {
	snprintf( msg, MSGLEN, "%s Response file %s has incorrect"
		     " units (ignoring units)", tag, *token );
	WARNING( msg );
      }
      transfer.sampleUnits = countsPerStrainUnit;
    }

    /* Copy new transfer function (and filename). */
    SUB( LALCDestroyVector( stat, &(detector->transfer->data) ),
	 stat );
    memcpy( detector->transfer, &transfer,
	    sizeof(COMPLEX8FrequencySeries) );
    SUB( LALCHARDestroyVector( stat, prevResponse ), stat );
    SUB( LALCHARCreateVector( stat, prevResponse,
			      strlen( *token ) + 1 ), stat );
    memcpy( (*prevResponse)->data, *token, strlen( *token ) + 1 );
  }

  /* Check if simulated noise was requested. */
  token = tokens->tokens + 2;
  if ( !strcmp( *token, "-" ) ) {
    INT8 epoch;  /* time series epoch */
    UINT4 npt;   /* time series length */
    UINT4 i;     /* index over data */
    REAL4 sigma; /* RMS time series value */

    PARSETOKEN( LALStringToI8, &epoch );
    PARSETOKEN( LALStringToU4, &npt );
    PARSETOKEN( LALStringToD, &(input.deltaT) );
    PARSETOKEN( LALStringToS, &sigma );
    if ( !ok ) {
      snprintf( msg, MSGLEN, "%s Error parsing noise argument %i"
		   " (skipping line)", tag,
		   (UINT4)( token - tokens->tokens - 2 ) );
      WARNING( msg );
      return 1;
    }
    I8ToLIGOTimeGPS( &(input.epoch), epoch );
    SUB( LALSCreateVector( stat, &(input.data), npt ), stat );
    SUB( LALNormalDeviates( stat, input.data, randParams ), stat );
    for ( i = 0; i < npt; i++ )
      input.data->data[i] *= sigma;
    input.sampleUnits = lalADCCountUnit;
    if ( sigma != 0.0 )
      snprintf( input.name, LALNameLength, "Simulated noise" );
    else
      input.name[0] = '\0';
  }

  /* Otherwise, read data from input file. */
  else {
    BOOLEAN unitsOK;  /* whether read units are right */
    LALUnitPair pair; /* input for unit comparison */
    FILE *fp = fopen( *token, "r" );
    if ( !fp ) {
      snprintf( msg, MSGLEN, "%s Input file %s not found"
		   " (skipping line)", tag, *token );
      WARNING( msg );
      return 1;
    }
    SUB( LALSReadTSeries( stat, &input, fp ), stat );
    fclose( fp );
    if ( input.data == NULL ) {
      snprintf( msg, MSGLEN, "%s Input file %s has no data"
		   " (skipping line)", tag, *token );
      WARNING( msg );
      return 1;
    }
    if ( input.deltaT <= 0.0 ) {
      snprintf( msg, MSGLEN, "%s Input file %s has bad metadata"
		   " (skipping line)", tag, *token );
      WARNING( msg );
      SUB( LALSDestroyVector( stat, &(input.data) ), stat );
      return 1;
    }

    /* Check units. */
    pair.unitOne = &lalADCCountUnit;
    pair.unitTwo = &(input.sampleUnits);
    SUB( LALUnitCompare( stat, &unitsOK, &pair ), stat );
    if ( !unitsOK ) {
      snprintf( msg, MSGLEN, "%s Input file %s has incorrect units"
		   " (ignoring units)", tag, *token );
      WARNING( msg );
      input.sampleUnits = lalADCCountUnit;
      input.sampleUnits = strainPerCountUnit;
    }
  }

  /* Simulate and inject signals. */
  here = head->next;
  snprintf( name, LALNameLength, "%s", input.name );
  while ( here ) {
    REAL4TimeSeries output = input;
    output.data = NULL;
    SUB( LALSCreateVector( stat, &(output.data), input.data->length ),
	 stat );
    SUB( LALSimulateCoherentGW( stat, &output, &(here->waveform),
				detector ), stat );
    SUB( LALSSInjectTimeSeries( stat, &input, &output ), stat );
    SUB( LALSDestroyVector( stat, &(output.data) ), stat );
    here = here->next;
    sourceno++;
  }

  /* Give the output a name. */
  if ( sourceno == 0 ) {
    snprintf( input.name, LALNameLength, "%s", name );
  } else if ( sourceno == 1 ) {
    if ( strlen( name ) )
      snprintf( input.name, LALNameLength,
		   "%s + 1 injected signal", name );
    else
      snprintf( input.name, LALNameLength, "1 injected signal" );
  } else {
    if ( strlen( name ) )
      snprintf( input.name, LALNameLength,
		   "%s + %i injected signals", name, sourceno );
    else
      snprintf( input.name, LALNameLength, "%i injected signals",
		   sourceno );
  }

  /* Write output to stdout or to requested output file. */
  token = tokens->tokens + tokens->nTokens - 1;
  if ( !strcmp( *token, "-" ) ) {
    SUB( LALSWriteTSeries( stat, stdout, &input ), stat );
  } else {
    FILE *fp = fopen( *token, "w" );
    if ( !fp ) {
      snprintf( msg, MSGLEN, "%s Output file %s could not be"
		   " opened (skipping line)", tag, *token );
      WARNING( msg );
      SUB( LALSDestroyVector( stat, &(input.data) ), stat );
      return 1;
    }
    SUB( LALSWriteTSeries( stat, fp, &input ), stat );
    fclose( fp );
  }

  /* We're done for this line; free input data. */
  SUB( LALSDestroyVector( stat, &(input.data) ), stat );
  return 0;
}


/* Function to compute the combination of two numbers. */
UINT4
choose( UINT4 a, UINT4 b )
{
  UINT4 numer = 1;
  UINT4 denom = 1;
  UINT4 index = b + 1;
  while ( --index ) {
    numer *= a - b + index;
    denom *= index;
  }
  return numer/denom;
}
