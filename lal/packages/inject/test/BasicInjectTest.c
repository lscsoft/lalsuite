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

/****************************** <lalVerbatim file="BasicInjectTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{BasicInjectTest.c}}
\label{ss:BasicInjectTest.c}

Injects inspiral signals into detector noise.

\subsubsection*{Usage}
\begin{verbatim}
BasicInjectTest [-s sourcefile] [-r respfile] [-o outfile] [-e seed]
                [-i infile | -n sec nsec npt dt sigma] [-d debuglevel]
\end{verbatim}

\subsubsection*{Description}

This program generates inspiral waveforms with specified parameters,
and injects them into ADC data.  The following option flags are
accepted:
\begin{itemize}
\item[\texttt{-s}] Reads source information from the file
\verb@sourcefile@.  If absent, it injects a single
1.4$M_\odot$--1.4$M_\odot$ inspiral, optimally oriented, at a distance
of $10^{-5}$ solar Schwarzschild radii ($0.00002GM_\odot/c^2$).
\item[\texttt{-r}] Reads a detector response function from the file
\verb@respfile@.  If absent, it generates raw dimensionless strain.
\item[\texttt{-i}] Reads ADC input from the file \verb@infile@.  This
takes precedence over the \verb@-n@ option, below.
\item[\texttt{-n}] Generates random ADC input data starting from a GPS
epoch of \verb@sec@ seconds plus \verb@nsec@ nanoseconds, with
\verb@npt@ data sampled at \verb@dt@ second intervals, with white
Gaussian noise having standard deviation \verb@sigma@.  If neither
\verb@-i@ (above) nor \verb@-n@ are given, the program assumes
\verb@-n 0 0 1048576 9.765625e-4 0.0@.
\item[\texttt{-o}] Writes injected ADC data to the file
\verb@outfile@.  If absent, the routines are exercised, but no output
is written.
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@.  If not
specified, level 0 is assumed.
\item[\texttt{-r}] Sets the random number seed to \verb@randomseed@.
If not specified, the seed is gerenated from the current time.
\end{itemize}

\paragraph{Format for \texttt{sourcefile}:} The source file consists
of any number of lines of data, each specifying a chirp waveform.
Each line must begin with a character code (\verb@CHAR@ equal to one
of \verb@'i'@, \verb@'f'@, or \verb@'c'@), followed by 6
whitespace-delimited numerical fields: the GPS epoch of the chirp
(\verb@INT8@ nanoseconds), the two binary masses (\verb@REAL4@
$M_\odot$), the distance to the source (\verb@REAL4@ kpc), and the
source's inclination and phase at coalescence (\verb@REAL4@ degrees).
The character codes have the following meanings:
\begin{itemize}
\item[\texttt{'i'}] The epoch represents the GPS time of the start of
the chirp waveform.
\item[\texttt{'f'}] The epoch represents the GPS time of the end of
the chirp waveform.
\item[\texttt{'c'}] The epoch represents the GPS time when the
binaries would coalesce in the point-mass approximation.
\end{itemize}
Thus a typical input line for two $1.4M_\odot$ objects at 11\,000\,kpc
inclined $30^\circ$ with an initial phase of $45^\circ$, coalescing at
315\,187\,245 GPS seconds, will have the following line in the input
file:
\begin{verbatim}
c 315187245000000000 1.4 1.4 11000.0 30.0 45.0
\end{verbatim}

\paragraph{Format for \texttt{respfile}:} The response function $R(f)$
gives the real and imaginary components of the transformation
\emph{from} ADC output $o$ \emph{to} tidal strain $h$ via
$\tilde{h}(f)=R(f)\tilde{o}(f)$.  It is inverted internally to give
the detector \emph{transfer function} $T(f)=1/R(f)$.  The format
\verb@respfile@ is a header specifying the GPS epoch $t_0$ at which
the response was taken (\verb@INT8@ nanoseconds), the lowest frequency
$f_0$ at which the response is given (\verb@REAL8@ Hz), and the
frequency sampling interval $\Delta f$ (\verb@REAL8@ Hz):

\medskip
\begin{tabular}{l}
\verb@# epoch = @$t_0$ \\
\verb@# f0 = @$f_0$ \\
\verb@# deltaF = @$\Delta f$
\end{tabular}
\medskip

\noindent followed by two columns of \verb@REAL4@ data giving the real
and imaginary components of $R(f_0+k\Delta f)$.

\paragraph{Format for \texttt{infile}:} The input file consists of a
header giving the GPS epoch $t_0$ of the first time sample
(\verb@INT8@ nanoseconds) and the sampling interval $\Delta t$
(\verb@REAL8@ seconds):

\medskip
\begin{tabular}{l}
\verb@# epoch = @$t_0$ \\
\verb@# deltaT = @$\Delta t$
\end{tabular}
\medskip

\noindent followed by a single column of ADC data.  The ADC data
should be integers in the range of an \verb@INT2@ (from $-32768$ to
$32767$), but is assumed to be written in floating-point notation in
accordance with frame format.

The output file \verb@outfile@ containing injected data is written in
the same format.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define BASICINJECTTESTC_ENORM  0
#define BASICINJECTTESTC_ESUB   1
#define BASICINJECTTESTC_EARG   2
#define BASICINJECTTESTC_EVAL   3
#define BASICINJECTTESTC_EFILE  4
#define BASICINJECTTESTC_EINPUT 5
#define BASICINJECTTESTC_EMEM   6

#define BASICINJECTTESTC_MSGENORM  "Normal exit"
#define BASICINJECTTESTC_MSGESUB   "Subroutine failed"
#define BASICINJECTTESTC_MSGEARG   "Error parsing arguments"
#define BASICINJECTTESTC_MSGEVAL   "Input argument out of valid range"
#define BASICINJECTTESTC_MSGEFILE  "Could not open file"
#define BASICINJECTTESTC_MSGEINPUT "Error reading file"
#define BASICINJECTTESTC_MSGEMEM   "Out of memory"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALMalloc()                     LALFree()
LALCreateRandomParams()         LALDestroyRandomParams()
LALI2CreateVector()             LALI2DestroyVector()
LALSCreateVector()              LALSDestroyVector()
LALCCreateVector()              LALCDestroyVector()
LALSReadVectorSequence()        LALSDestroyVectorSequence()
LALCCVectorDivide()             LALGeneratePPNInspiral()
LALSimulateCoherentGW()         LALSI2InjectTimeSeries()
LALNormalDeviates()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{BasicInjectTestCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/VectorOps.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/StreamInput.h>

NRCSID( BASICINJECTTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;
#define EPOCH (0)
#define DIST  (0.00002*LAL_MRSUN_SI )
#define M1    (1.4)
#define M2    (1.4)
#define INC   (0.0)
#define PHIC  (0.0)
#define SEC   (0)
#define NSEC  (0)
#define NPT   (1048576)
#define DT    (1.0/1024.0)
#define SIGMA (0.0)

/* Global constants. */
#define MSGLEN (256)    /* maximum length of warning/info messages */
#define FSTART (25.0)   /* initial frequency of waveform */
#define FSTOP  (500.0) /* termination frequency of waveform */
#define DELTAT (0.01)   /* sampling interval of amplitude and phase */

/* Usage format string. */
#define USAGE "Usage: %s [-s sourcefile] [-r respfile] [-o outfile] [-e seed]\n                [-i infile | -n sec nsec npt dt sigma] [-d debuglevel]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, BASICINJECTTESTC, statement ? statement : \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 BASICINJECTTESTC, (statement) );                    \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 BASICINJECTTESTC, (statement) );                    \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( BASICINJECTTESTC_ESUB, BASICINJECTTESTC_MSGESUB,            \
         "Function call \"" #func "\" failed:" );                    \
  return BASICINJECTTESTC_ESUB;                                      \
}                                                                    \
while (0)

#define CHECKVAL( val, lower, upper )                                \
do                                                                   \
if ( ( (val) <= (lower) ) || ( (val) > (upper) ) )                   \
{                                                                    \
  ERROR( BASICINJECTTESTC_EVAL, BASICINJECTTESTC_MSGEVAL,            \
         "Value of " #val " out of range:" );                        \
  LALPrintError( #val " = %f, range = (%f,%f]\n", (REAL8)(val),      \
                 (REAL8)(lower), (REAL8)(upper) );                   \
  return BASICINJECTTESTC_EVAL;                                      \
}                                                                    \
while (0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

/* A function to convert INT8 nanoseconds to LIGOTimeGPS. */
void
I8ToLIGOTimeGPS( LIGOTimeGPS *output, INT8 input );


int
main(int argc, char **argv)
{
  /* Command-line parsing variables. */
  int arg;                   /* command-line argument counter */
  static LALStatus stat;     /* status structure */
  CHAR *sourcefile = NULL;   /* name of sourcefile */
  CHAR *respfile = NULL;     /* name of respfile */
  CHAR *infile = NULL;       /* name of infile */
  CHAR *outfile = NULL;      /* name of outfile */
  INT4 seed = 0;             /* random number seed */
  INT4 sec = SEC;            /* ouput epoch.gpsSeconds */
  INT4 nsec = NSEC;          /* ouput epoch.gpsNanoSeconds */
  INT4 npt = NPT;            /* number of output points */
  REAL8 dt = DT;             /* output sampling interval */
  REAL4 sigma = SIGMA;       /* noise amplitude */

  /* File reading variables. */
  FILE *fp = NULL; /* generic file pointer */
  BOOLEAN ok = 1;  /* whether input format is correct */
  UINT4 i;         /* generic index over file lines */
  INT8 epoch;      /* epoch stored as an INT8 */

  /* Other global variables. */
  RandomParams *params = NULL; /* parameters of pseudorandom sequence */
  DetectorResponse detector;   /* the detector in question */
  INT2TimeSeries output;       /* detector ACD output */


  /*******************************************************************
   * PARSE ARGUMENTS (arg stores the current position)               *
   *******************************************************************/

  /* Exit gracefully if no arguments were given.
  if ( argc <= 1 ) {
    INFO( "No testing done." );
    return 0;
  } */

  arg = 1;
  while ( arg < argc ) {
    /* Parse source file option. */
    if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	sourcefile = argv[arg++];
      }else{
	ERROR( BASICINJECTTESTC_EARG, BASICINJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BASICINJECTTESTC_EARG;
      }
    }
    /* Parse response file option. */
    else if ( !strcmp( argv[arg], "-r" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	respfile = argv[arg++];
      }else{
	ERROR( BASICINJECTTESTC_EARG, BASICINJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BASICINJECTTESTC_EARG;
      }
    }
    /* Parse input file option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	infile = argv[arg++];
      }else{
	ERROR( BASICINJECTTESTC_EARG, BASICINJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BASICINJECTTESTC_EARG;
      }
    }
    /* Parse noise output option. */
    else if ( !strcmp( argv[arg], "-n" ) ) {
      if ( argc > arg + 5 ) {
	arg++;
	sec = atoi( argv[arg++] );
	nsec = atoi( argv[arg++] );
	npt = atoi( argv[arg++] );
	dt = atof( argv[arg++] );
	sigma = atof( argv[arg++] );
      }else{
	ERROR( BASICINJECTTESTC_EARG, BASICINJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BASICINJECTTESTC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      }else{
	ERROR( BASICINJECTTESTC_EARG, BASICINJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BASICINJECTTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      }else{
	ERROR( BASICINJECTTESTC_EARG, BASICINJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BASICINJECTTESTC_EARG;
      }
    }
    /* Parse random seed option. */
    else if ( !strcmp( argv[arg], "-e" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	seed = atoi( argv[arg++] );
      }else{
	ERROR( BASICINJECTTESTC_EARG, BASICINJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BASICINJECTTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( BASICINJECTTESTC_EARG, BASICINJECTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return BASICINJECTTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Check for redundant or bad argument values. */
  CHECKVAL( npt, 0, 2147483647 );
  CHECKVAL( dt, 0, LAL_REAL4_MAX );


  /*******************************************************************
   * SETUP                                                           *
   *******************************************************************/

  /* Set up output, detector, and random parameter structures. */
  output.data = NULL;
  detector.transfer = (COMPLEX8FrequencySeries *)
    LALMalloc( sizeof(COMPLEX8FrequencySeries) );
  if ( !(detector.transfer) ) {
    ERROR( BASICINJECTTESTC_EMEM, BASICINJECTTESTC_MSGEMEM, 0 );
    return BASICINJECTTESTC_EMEM;
  }
  detector.transfer->data = NULL;
  detector.site = NULL;
  SUB( LALCreateRandomParams( &stat, &params, seed ), &stat );

  /* Set up units. */
  {
    RAT4 negOne = { -1, 0 };
    LALUnit unit;
    LALUnitPair pair;
    output.sampleUnits = lalADCCountUnit;
    pair.unitOne = &lalADCCountUnit;
    pair.unitTwo = &lalStrainUnit;
    SUB( LALUnitRaise( &stat, &unit, pair.unitTwo,
		       &negOne ), &stat );
    pair.unitTwo = &unit;
    SUB( LALUnitMultiply( &stat, &(detector.transfer->sampleUnits),
			  &pair ), &stat );
  }

  /* Read response function. */
  if ( respfile ) {
    REAL4VectorSequence *resp = NULL; /* response as vector sequence */
    COMPLEX8Vector *response = NULL;  /* response as complex vector */
    COMPLEX8Vector *unity = NULL;     /* vector of complex 1's */

    if ( ( fp = fopen( respfile, "r" ) ) == NULL ) {
      ERROR( BASICINJECTTESTC_EFILE, BASICINJECTTESTC_MSGEFILE,
	     respfile );
      return BASICINJECTTESTC_EFILE;
    }

    /* Read header. */
    ok &= ( fscanf( fp, "# epoch = %lld\n", &epoch ) == 1 );
    I8ToLIGOTimeGPS( &( detector.transfer->epoch ), epoch );
    ok &= ( fscanf( fp, "# f0 = %lf\n", &( detector.transfer->f0 ) )
	    == 1 );
    ok &= ( fscanf( fp, "# deltaF = %lf\n",
		    &( detector.transfer->deltaF ) ) == 1 );
    if ( !ok ) {
      ERROR( BASICINJECTTESTC_EINPUT, BASICINJECTTESTC_MSGEINPUT,
	     respfile );
      return BASICINJECTTESTC_EINPUT;
    }

    /* Read and convert body to a COMPLEX8Vector. */
    SUB( LALSReadVectorSequence( &stat, &resp, fp ), &stat );
    fclose( fp );
    if ( resp->vectorLength != 2 ) {
      ERROR( BASICINJECTTESTC_EINPUT, BASICINJECTTESTC_MSGEINPUT,
	     respfile );
      return BASICINJECTTESTC_EINPUT;
    }
    SUB( LALCCreateVector( &stat, &response, resp->length ), &stat );
    memcpy( response->data, resp->data, 2*resp->length*sizeof(REAL4) );
    SUB( LALSDestroyVectorSequence( &stat, &resp ), &stat );

    /* Convert response function to a transfer function. */
    SUB( LALCCreateVector( &stat, &unity, response->length ), &stat );
    for ( i = 0; i < response->length; i++ ) {
      unity->data[i].re = 1.0;
      unity->data[i].im = 0.0;
    }
    SUB( LALCCreateVector( &stat, &( detector.transfer->data ),
			   response->length ), &stat );
    SUB( LALCCVectorDivide( &stat, detector.transfer->data, unity,
			    response ), &stat );
    SUB( LALCDestroyVector( &stat, &response ), &stat );
    SUB( LALCDestroyVector( &stat, &unity ), &stat );
  }

  /* No response file, so generate a unit response. */
  else {
    I8ToLIGOTimeGPS( &( detector.transfer->epoch ), EPOCH );
    detector.transfer->f0 = 0.0;
    detector.transfer->deltaF = 1.5*FSTOP;
    SUB( LALCCreateVector( &stat, &( detector.transfer->data ), 2 ),
	 &stat );
    detector.transfer->data->data[0].re = 1.0;
    detector.transfer->data->data[1].re = 1.0;
    detector.transfer->data->data[0].im = 0.0;
    detector.transfer->data->data[1].im = 0.0;
  }


  /* Read input data. */
  if ( infile ) {
    REAL4VectorSequence *input = NULL; /* input as vector sequence */
    if ( ( fp = fopen( infile, "r" ) ) == NULL ) {
      ERROR( BASICINJECTTESTC_EFILE, BASICINJECTTESTC_MSGEFILE,
	     infile );
      return BASICINJECTTESTC_EFILE;
    }

    /* Read header. */
    ok &= ( fscanf( fp, "# epoch = %lld\n", &epoch ) == 1 );
    I8ToLIGOTimeGPS( &( output.epoch ), epoch );
    ok &= ( fscanf( fp, "# deltaT = %lf\n", &( output.deltaT ) )
	    == 1 );
    if ( !ok ) {
      ERROR( BASICINJECTTESTC_EINPUT, BASICINJECTTESTC_MSGEINPUT,
	     infile );
      return BASICINJECTTESTC_EINPUT;
    }

    /* Read and convert body. */
    SUB( LALSReadVectorSequence( &stat, &input, fp ), &stat );
    fclose( fp );
    if ( input->vectorLength != 1 ) {
      ERROR( BASICINJECTTESTC_EINPUT, BASICINJECTTESTC_MSGEINPUT,
	     infile );
      return BASICINJECTTESTC_EINPUT;
    }
    SUB( LALI2CreateVector( &stat, &( output.data ), input->length ),
	 &stat );
    for ( i = 0; i < input->length; i++ )
      output.data->data[i] = (INT2)( input->data[i] );
    SUB( LALSDestroyVectorSequence( &stat, &input ), &stat );
  }

  /* No input file, so generate one randomly. */
  else {
    output.epoch.gpsSeconds = sec;
    output.epoch.gpsNanoSeconds = nsec;
    output.deltaT = dt;
    SUB( LALI2CreateVector( &stat, &( output.data ), npt ), &stat );
    if ( sigma == 0 ) {
      memset( output.data->data, 0, npt*sizeof(INT2) );
    } else {
      REAL4Vector *deviates = NULL; /* unit Gaussian deviates */
      SUB( LALSCreateVector( &stat, &deviates, npt ), &stat );
      SUB( LALNormalDeviates( &stat, deviates, params ), &stat );
      for ( i = 0; i < (UINT4)( npt ); i++ )
	output.data->data[i] = (INT2)
	  floor( sigma*deviates->data[i] + 0.5 );
      SUB( LALSDestroyVector( &stat, &deviates ), &stat );
    }
  }


  /*******************************************************************
   * INJECTION                                                       *
   *******************************************************************/

  /* Open sourcefile. */
  if ( sourcefile ) {
    if ( ( fp = fopen( sourcefile, "r" ) ) == NULL ) {
      ERROR( BASICINJECTTESTC_EFILE, BASICINJECTTESTC_MSGEFILE,
	     sourcefile );
      return BASICINJECTTESTC_EFILE;
    }
  }

  /* For each line in the sourcefile... */
  while ( ok ) {
    PPNParamStruc ppnParams;       /* wave generation parameters */
    REAL4 m1, m2, dist, inc, phic; /* unconverted parameters */
    CoherentGW waveform;           /* amplitude and phase structure */
    REAL4TimeSeries signalvec;     /* GW signal */
    REAL8 time;                    /* length of GW signal */
    CHAR timeCode;                 /* code for signal time alignment */
    CHAR message[MSGLEN];          /* warning/info messages */

    /* Read and convert input line. */
    if ( sourcefile ) {
      ok &= ( fscanf( fp, "%c %lli %f %f %f %f %f\n", &timeCode,
		      &epoch, &m1, &m2, &dist, &inc, &phic ) == 7 );
      ppnParams.mTot = m1 + m2;
      ppnParams.eta = m1*m2/( ppnParams.mTot*ppnParams.mTot );
      ppnParams.d = dist*LAL_PC_SI*1000.0;
      ppnParams.inc = inc*LAL_PI/180.0;
      ppnParams.phi = phic*LAL_PI/180.0;
    } else {
      timeCode = 'i';
      ppnParams.mTot = M1 + M2;
      ppnParams.eta = M1*M2/( ppnParams.mTot*ppnParams.mTot );
      ppnParams.d = DIST;
      ppnParams.inc = INC;
      ppnParams.phi = PHIC;
      epoch = EPOCH;
    }

    if ( ok ) {
      /* Set up other parameter structures. */
      ppnParams.epoch.gpsSeconds = ppnParams.epoch.gpsNanoSeconds = 0;
      ppnParams.position.latitude = ppnParams.position.longitude = 0.0;
      ppnParams.position.system = COORDINATESYSTEM_EQUATORIAL;
      ppnParams.psi = 0.0;
      ppnParams.fStartIn = FSTART;
      ppnParams.fStopIn = FSTOP;
      ppnParams.lengthIn = 0;
      ppnParams.ppn = NULL;
      ppnParams.deltaT = DELTAT;
      memset( &waveform, 0, sizeof(CoherentGW) );

      /* Generate waveform at zero epoch. */
      SUB( LALGeneratePPNInspiral( &stat, &waveform, &ppnParams ),
	   &stat );
      snprintf( message, MSGLEN, "%d: %s", ppnParams.termCode,
		   ppnParams.termDescription );
      INFO( message );
      if ( ppnParams.dfdt > 2.0 ) {
	snprintf( message, MSGLEN,
		     "Waveform sampling interval is too large:\n"
		     "\tmaximum df*dt = %f", ppnParams.dfdt );
	WARNING( message );
      }

      /* Compute epoch for waveform. */
      time = waveform.a->data->length*DELTAT;
      if ( timeCode == 'f' )
	epoch -= (INT8)( 1000000000.0*time );
      else if ( timeCode == 'c' )
	epoch -= (INT8)( 1000000000.0*ppnParams.tc );
      I8ToLIGOTimeGPS( &( waveform.a->epoch ), epoch );
      waveform.f->epoch = waveform.phi->epoch = waveform.a->epoch;

      /* Generate and inject signal. */
      signalvec.epoch = waveform.a->epoch;
      signalvec.epoch.gpsSeconds -= 1;
      signalvec.deltaT = output.deltaT/4.0;
      signalvec.f0 = 0;
      signalvec.data = NULL;
      time = ( time + 2.0 )/signalvec.deltaT;
      SUB( LALSCreateVector( &stat, &( signalvec.data ), (UINT4)time ),
	   &stat );
      SUB( LALSimulateCoherentGW( &stat, &signalvec, &waveform,
				  &detector ), &stat );
      SUB( LALSI2InjectTimeSeries( &stat, &output, &signalvec, params ),
	   &stat );
      SUB( LALSDestroyVectorSequence( &stat, &( waveform.a->data ) ),
	   &stat );
      SUB( LALSDestroyVector( &stat, &( waveform.f->data ) ), &stat );
      SUB( LALDDestroyVector( &stat, &( waveform.phi->data ) ), &stat );
      LALFree( waveform.a );
      LALFree( waveform.f );
      LALFree( waveform.phi );
      SUB( LALSDestroyVector( &stat, &( signalvec.data ) ), &stat );
    }

    /* If there is no source file, inject only one source. */
    if ( !sourcefile )
      ok = 0;
  }

  /* Input file is exhausted (or has a badly-formatted line ). */
  if ( sourcefile )
    fclose( fp );


  /*******************************************************************
   * CLEANUP                                                         *
   *******************************************************************/

  /* Print output file. */
  if ( outfile ) {
    if ( ( fp = fopen( outfile, "w" ) ) == NULL ) {
      ERROR( BASICINJECTTESTC_EFILE, BASICINJECTTESTC_MSGEFILE,
	     outfile );
      return BASICINJECTTESTC_EFILE;
    }
    epoch = 1000000000LL*(INT8)( output.epoch.gpsSeconds );
    epoch += (INT8)( output.epoch.gpsNanoSeconds );
    fprintf( fp, "# epoch = %lld\n", epoch );
    fprintf( fp, "# deltaT = %23.16e\n", output.deltaT );
    for ( i = 0; i < output.data->length; i++ )
      fprintf( fp, "%8.1f\n", (REAL4)( output.data->data[i] ) );
    fclose( fp );
  }

  /* Destroy remaining memory. */
  SUB( LALDestroyRandomParams( &stat, &params ), &stat );
  SUB( LALI2DestroyVector( &stat, &( output.data ) ), &stat );
  SUB( LALCDestroyVector( &stat, &( detector.transfer->data ) ),
       &stat );
  LALFree( detector.transfer );

  /* Done! */
  LALCheckMemoryLeaks();
  INFO( BASICINJECTTESTC_MSGENORM );
  return BASICINJECTTESTC_ENORM;
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
