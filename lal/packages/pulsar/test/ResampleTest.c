/********************************* <lalVerbatim file="ResampleTestCV">
Author: Creighton, T. D.
Revision: $Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Program \texttt{ResampleTest.c}}
\label{ss:ResampleTest.c}

Tests the routines in \verb@Resample.h@.

\subsubsection*{Usage}
\begin{verbatim}
ResampleTest [-d debuglevel] [-p psfile] [-t tfile] [-c n a f] [-m dec df fm]
\end{verbatim}

This program generates a quasiperiodic time series having a sinusoidal
phase modulation, and generates a piecewise-polynomial fit to the
phase function.  It then generates and applies stroboscopic resampling
rules to the time series to produce a monochromatic signal.  The
following option flags are accepted:
\begin{itemize}
\item[\texttt{-d}] Sets the global \verb@lalDebugLevel@ to the
specified \verb@debuglevel@.
\item[\texttt{-p}] Power spectra of the time series before and after
demodulation will be written to the file \verb@psfile@.
\item[\texttt{-t}] The timing difference function $(\tau-t)/\Delta t$
computed in three ways (analytically, from a polynomial fit, and from
the resampling rules), will be written to the file \verb@tfile@.  See
below for notation.
\item[\texttt{-c}] Sets parameters for the ``carrier'' signal: the
number of points \verb@n@, the amplitude \verb@a@, and the frequency
\verb@f@.
\item[\texttt{-m}] Sets parameters for the signal modulation and
resampling: the decimation factor \verb@dec@, the maximum change in
signal frequency \verb@df@, and the frequency of the modulation
\verb@fm@.
\end{itemize}
All frequencies are in units of the sampling rate, which is an
arbitrary scale.  With no options, the program runs with
\verb@lalDebugLevel@=0, produces no output, and uses internally
\verb@#define@d signal parameters.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define RESAMPLETESTC_ENORM 0
#define RESAMPLETESTC_ESUB  1
#define RESAMPLETESTC_EARG  2
#define RESAMPLETESTC_EBAD  3
#define RESAMPLETESTC_EFILE 4

#define RESAMPLETESTC_MSGENORM "Success, normal exit"
#define RESAMPLETESTC_MSGESUB  "Recursive error"
#define RESAMPLETESTC_MSGEARG  "Error parsing arguments"
#define RESAMPLETESTC_MSGEBAD  "Bad argument value"
#define RESAMPLETESTC_MSGEFILE "Error opening or writing to output file"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

The modulated signal is of the form $s(t)=A\sin[\phi(t)]$, where $A$
is a rather arbitrary amplitude, and $\phi(t)$ is a phase function of
the form:
$$
\phi(t) = 2\pi f_c t + \frac{\Delta f}{f_m}\sin(2\pi f_m t) \; .
$$
Here $f_c$ is the average ``carrier'' frequency, $f_m$ is the
frequency of the modulation, and $\Delta f$ is the maximum change in
the ``instantaneous'' frequency of the signal (or, equivalently, the
separation between the carrier and sidebands in the power spectral
density).  The canonical (demodulated) time coordinate for this phase
function is $\tau=\phi/2\pi f_c$.  The demodulation routines require
quadratic fits to the function $\tau-t$ at various times $t_0$:
\begin{eqnarray}
\tau - t & = & \frac{(\Delta f/f_c)}{2\pi f_m}\sin(2\pi f_m t_0)
		\nonumber\\
         & + & (\Delta f/f_c)\cos(2\pi f_m t_0)(t-t_0) \nonumber\\
         & - & \pi f_m (\Delta f/f_c)\sin(2\pi f_m t_0)(t-t_0)^2 \; ,
\label{eq:polyco-formulae}
\end{eqnarray}
with residuals less than $(2/3)\pi^2 f_m^2(\Delta f/f_c)(t-t_0)^3$.
We require this residual to be always less than one sample interval
$\Delta t$.  This means that a piecewise-quadratic fit to the phase
function must be evaluated at times $t_0$ separated by no more than:
\begin{equation}
\Delta t_0 \lessim \sqrt[3]{\frac{12f_c\Delta t}{\pi^2f_m^2\Delta f}}
	\; ,
\label{eq:polyco-interval}
\end{equation}
noting that each piecewise fit is good for a time interval
$t_0\pm\Delta t_0/2$ about each central time $t_0$.

Thus to create a piecewise-polynomial fit defined by
\verb@PolycoStruc@, this program simply define a set of fitting times
\verb@t0[@$k$\verb@]@$=(2k+1)\Delta t_0/2$, and computes the
appropriate components of \verb@polyco@ from
Eq.~(\ref{eq:polyco-formulae}), above.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()
LALSCreateVector()
LALSDestroyVector()
LALCreateResampleRules()
LALApplyResampleRules()
LALDestroyResampleRules()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ResampleTestCV}}

******************************************************* </lalLaTeX> */

#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/Resample.h>

NRCSID( RESAMPLETESTC, "$Id$" );

/* Default parameter settings. */
INT4 lalDebugLevel = 0;
#define NPTS 4096
#define AMP 1
#define FREQ 0.1
#define DEC 4
#define DF 0.02
#define FM 0.001

/* Output files for timing information. */
#define OUTFILE0 "out0.dat"
#define OUTFILE1 "out1.dat"
#define OUTFILE2 "out2.dat"

/* Other constants. */
#define NPOLY 3    /* Number of polynomial coefficients (quadratic) */
#define BUFFER 5.0 /* Polynomial fit covers this many seconds before
                      and after requested timespan */
REAL8 df_f;    /* Maximum fractional change in carrier frequency */
REAL8 twoPiFm; /* 2*Pi times the frequency of modulation */

/* Usage format string. */
#define USAGE "Usage: %s [-d debuglevel] [-p psfile] [-t tfile] [-c n a f] [-m dec df fm]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, RESAMPLETESTC, statement ? statement : "",\
		 (msg) );                                            \
}                                                                    \
else (void)(0)

#define INFO( statement )                                            \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 RESAMPLETESTC, (statement) );                       \
}                                                                    \
else (void)(0)

#define SUB( func, statusptr )                                       \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( RESAMPLETESTC_ESUB, RESAMPLETESTC_MSGESUB,                  \
         "Function call \"" #func "\" failed:" );                    \
  return RESAMPLETESTC_ESUB;                                         \
}                                                                    \
else (void)(0)


/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif


/* Local subroutines. */
REAL8
TDiff( REAL8 t );
/* Returns the value tau-t as a function of t. */

REAL8
DTDiff( REAL8 t );
/* Returns the derivative of tau-t as a function of t. */

REAL8
DDTDiff( REAL8 t );
/* Returns the double derivative of tau-t as a function of t. */

REAL8
DDDTDiffMax( void );
/* Returns the triple derivative of tau-t as a function of t,
   maximized over all t. */

static /* const */ REAL4TimeSeries emptyREAL4TimeSeries;
static /* const */ CreateVectorSequenceIn emptyCreateVectorSequenceIn;
static /* const */ ResampleParamStruc emptyResampleParamStruct;

int
main( int argc, char **argv )
{
  static LALStatus stat;/* Head of status list */
  CHAR *tfile = NULL;   /* Timing difference output filename */
  CHAR *psfile = NULL;  /* Power spectrum output filename */
  INT4 arg;
  UINT4 i;              /* Index counter */
  UINT4 n = NPTS;       /* Number of points in time series */
  UINT4 m = NPTS/DEC;   /* Number of points in decimated series */
  UINT4 dec = DEC;      /* Decimation factor */
  REAL4 a = AMP;        /* Amplitude of carrier wave */
  REAL4 f = FREQ;       /* Frequency of carrier wave */
  REAL4 df = DF;        /* Amplitude of frequency modulation */
  REAL4 fm = FM;        /* Frequency of modulation */
  REAL4TimeSeries input = emptyREAL4TimeSeries;  /* Modulated time series */
  REAL4TimeSeries output = emptyREAL4TimeSeries; /* Decimated time series */
  static PolycoStruc polyco;     /* Polynomial coefficient structure */
  ResampleRules *rules = NULL; /* Resampling rules structure */

  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse carrier parameters option. */
    if ( !strcmp( argv[arg], "-c" ) ) {
      if ( argc > arg + 3 ) {
        arg++;
        n = atoi( argv[arg++] );
        a = atof( argv[arg++] );
        f = atof( argv[arg++] );
      } else {
	ERROR( RESAMPLETESTC_EARG, RESAMPLETESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return RESAMPLETESTC_EARG;
      }
    }
    /* Parse modulation parameters option. */
    else if ( !strcmp( argv[arg], "-m" ) ) {
      if ( argc > arg + 3 ) {
        arg++;
        dec = atoi( argv[arg++] );
        df = atof( argv[arg++] );
        fm = atof( argv[arg++] );
      } else {
	ERROR( RESAMPLETESTC_EARG, RESAMPLETESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return RESAMPLETESTC_EARG;
      }
    }
    /* Parse debuglevel option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        lalDebugLevel = atoi( argv[arg++] );
      } else {
	ERROR( RESAMPLETESTC_EARG, RESAMPLETESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return RESAMPLETESTC_EARG;
      }
    }
    /* Parse power spectrum output option. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        psfile = argv[arg++];
      } else {
	ERROR( RESAMPLETESTC_EARG, RESAMPLETESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return RESAMPLETESTC_EARG;
      }
    }
    /* Parse timing difference output option. */
    else if ( !strcmp( argv[arg], "-t" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        tfile = argv[arg++];
      } else {
	ERROR( RESAMPLETESTC_EARG, RESAMPLETESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return RESAMPLETESTC_EARG;
      }
    }
    /* Unrecognized option. */
    else {
      ERROR( RESAMPLETESTC_EARG, RESAMPLETESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return RESAMPLETESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Make sure that parameter values won't crash the system. */
  if ( n*dec*f*df*fm == 0 ) {
    ERROR( RESAMPLETESTC_EBAD, RESAMPLETESTC_MSGEBAD, 0 );
    return RESAMPLETESTC_EBAD;
  }
  m = n/dec;
  df_f = df/f;
  twoPiFm = LAL_TWOPI*fm;

  /* Generate modulated time series. */
  /* input.name = "Modulated sinusoid"; */
  input.deltaT = 1.0;
  SUB( LALSCreateVector( &stat, &(input.data), n ), &stat );
  {
    REAL4 phi_c = LAL_TWOPI*f;      /* Factor preceding tau */
    REAL4 *data = input.data->data; /* Generic pointer to data */
    for ( i=0; i < n; i++ )
      *(data++) = cos( phi_c*( i + TDiff( i ) ) );
  }

  /* Create polynomial coefficients. */
  {
    REAL8 t0; /* Fitting points */
    REAL8 dtMax = 2.0*pow( 6.0/DDDTDiffMax(), 1.0/3.0 );
    /* Maximum separation between fitting points */
    UINT4 n0 = (UINT4)( ( n + 2*BUFFER)/dtMax ) + 1;
    /* Number of fits required to span the desired input time */
    REAL4 *data1, *data2, *data3; /* Generic pointers to data */
    CreateVectorSequenceIn in = emptyCreateVectorSequenceIn;

    /* Set fields of polyco. */
    polyco.start.gpsSeconds = (INT4)( -BUFFER );
    polyco.start.gpsNanoSeconds = (INT4)
      ( 1.0e9*( -BUFFER - polyco.start.gpsSeconds ) );
    SUB( LALSCreateVector( &stat, &(polyco.tBound), n0 ), &stat );
    SUB( LALSCreateVector( &stat, &(polyco.t0), n0 ), &stat );
    in.length = n0;
    in.vectorLength = NPOLY;
    SUB( LALSCreateVectorSequence( &stat, &(polyco.polyco), &in ),
	 &stat );
    memset( polyco.polyco->data, 0, in.length*in.vectorLength );

    /* Fill in data. */
    t0=0.5*dtMax-BUFFER;
    data1=polyco.polyco->data;
    data2=polyco.tBound->data;
    data3=polyco.t0->data;
    while( n0-- ) {
      data1[0] = TDiff( t0 );
      if ( NPOLY > 0 )
	data1[1] = DTDiff( t0 );
      if ( NPOLY > 1 )
	data1[2] = 0.5*DDTDiff( t0 );
      data1 += NPOLY;
      *(data2++) = t0 + 0.5*dtMax + BUFFER;
      *(data3++) = t0 + BUFFER;
      t0 += dtMax;
    }
  }

  /* Generate demodulated time series. */
  /* output.name = "Decimated sinusoid"; */
  output.deltaT = (REAL8)( dec );
  SUB( LALSCreateVector( &stat, &(output.data), m ), &stat );
  {
    ResampleParamStruc params = emptyResampleParamStruct;
    params.start.gpsSeconds = -1;
    params.start.gpsNanoSeconds = 0;
    params.stop.gpsSeconds = n + 1;
    params.stop.gpsNanoSeconds = 0;
    params.deltaT = 1.0;
    params.decimate = dec;
    SUB( LALCreateResampleRules( &stat, &rules, &polyco, &params ),
	 &stat );
    SUB( LALApplyResampleRules( &stat, &output, &input, rules ),
	 &stat );
  }

  /* Print modulation function, if requested. */
  if ( tfile ) {
    REAL4 *data0, *data1, *data2; /* Pointers to data */
    REAL4TimeSeries diffs0 = emptyREAL4TimeSeries;  /* Analytic timing difference */
    REAL4TimeSeries diffs1 = emptyREAL4TimeSeries;  /* Polynomial timing difference */
    REAL4TimeSeries diffs2 = emptyREAL4TimeSeries;  /* Resampled timing difference */
    FILE *fp = fopen( tfile, "w" ); /* Output file pointer */

    /* Make sure output file was created. */
    if ( !fp ) {
      ERROR( RESAMPLETESTC_EFILE, RESAMPLETESTC_MSGEFILE, 0 );
      return RESAMPLETESTC_EFILE;
    }

    /* Create the three time series. */
    diffs0.deltaT = diffs1.deltaT = diffs2.deltaT = 1.0;
    SUB( LALSCreateVector( &stat, &(diffs0.data), n ), &stat );
    SUB( LALSCreateVector( &stat, &(diffs1.data), n ), &stat );
    SUB( LALSCreateVector( &stat, &(diffs2.data), n ), &stat );

    /* Fill the three time series. */
    data0 = diffs0.data->data;
    for ( i=0; i < n; i++ )
      *(data0++) = TDiff( i );
    SUB( LALPolycoToTimingDifference( &stat, &diffs1, &polyco ),
	 &stat );
    SUB( LALRulesToTimingDifference( &stat, &diffs2, rules ), &stat );

    /* Print the result. */
    i = n;
    data0 = diffs0.data->data;
    data1 = diffs1.data->data;
    data2 = diffs2.data->data;
    for ( i=0; i < n; i++ )
      if ( fprintf( fp, "%10i %10.3e %10.3e %10.3e\n", i, *(data0++),
		    *(data1++), *(data2++) ) < 44 ) {
	ERROR( RESAMPLETESTC_EFILE, RESAMPLETESTC_MSGEFILE, 0 );
	return RESAMPLETESTC_EFILE;
      }
    fclose( fp );

    /* Destroy the three time series. */
    SUB( LALSDestroyVector( &stat, &(diffs0.data) ), &stat );
    SUB( LALSDestroyVector( &stat, &(diffs1.data) ), &stat );
    SUB( LALSDestroyVector( &stat, &(diffs2.data) ), &stat );
  }

  /* Destroy polynomial coefficients and resampling rules. */
  SUB( LALSDestroyVectorSequence( &stat, &(polyco.polyco) ), &stat );
  SUB( LALSDestroyVector( &stat, &(polyco.t0) ), &stat );
  SUB( LALSDestroyVector( &stat, &(polyco.tBound) ), &stat );
  SUB( LALDestroyResampleRules( &stat, &rules ), &stat );

  /* Write power spectra to output file. */
  if ( psfile ) {
    REAL4 *data1, *data2;      /* Generic pointers to data */
    REAL4Vector *PSIn = NULL;  /* Power spectrum of input */
    REAL4Vector *PSOut = NULL; /* Power spectrum of output */
    RealFFTPlan *plan = NULL;  /* Plan for computing power spectra */
    FILE *fp = fopen( psfile, "w" ); /* Output file pointer */

    /* Make sure output file was created. */
    if ( !fp ) {
      ERROR( RESAMPLETESTC_EFILE, RESAMPLETESTC_MSGEFILE, 0 );
      return RESAMPLETESTC_EFILE;
    }

    /* Generate power spectrum of output. */
    SUB( LALCreateForwardRealFFTPlan( &stat, &plan, m, 0 ), &stat );
    SUB( LALSCreateVector( &stat, &PSOut, m/2 + 1 ), &stat );
    SUB( LALRealPowerSpectrum( &stat, PSOut, output.data, plan ),
	 &stat );
    /* CHANGE BACK TO ORIGINAL NORMALIZATION -- JC */
    {
      REAL4Vector *myvector = PSOut;
      UINT4 mybin;
      for ( mybin = 1; mybin < myvector->length - 1; ++mybin )
        myvector->data[mybin] *= 0.5;
    }


    /* Create a new output that is a decimated but not demodulated
       data set, then discard the initial input. */
    for ( i=0, data1=input.data->data, data2=output.data->data; i < m;
	  i++, data1+=dec, data2++ )
      *data2 = *data1;
    SUB( LALSDestroyVector( &stat, &(input.data) ), &stat );

    /* Generate its power spectrum, then discard time series data (and
       the FFT plan). */
    SUB( LALSCreateVector( &stat, &PSIn, m/2 + 1 ), &stat );
    SUB( LALRealPowerSpectrum( &stat, PSIn, output.data, plan ),
	 &stat );
    SUB( LALSDestroyVector( &stat, &(output.data) ), &stat );
    SUB( LALDestroyRealFFTPlan( &stat, &plan ), &stat );

    /* Print both power spectra to the output file and close it. */
    data1 = PSIn->data;
    data2 = PSOut->data;
    for ( i = 0; i <= m/2; i++ )
      if ( fprintf( fp, "%10i %10.3e %10.3e\n", i, *(data1++),
		    *(data2++) ) < 33 ) {
	ERROR( RESAMPLETESTC_EFILE, RESAMPLETESTC_MSGEFILE, 0 );
	return RESAMPLETESTC_EFILE;
      }
    fclose( fp );

    /* Free local memory. */
    SUB( LALSDestroyVector( &stat, &PSIn ), &stat );
    SUB( LALSDestroyVector( &stat, &PSOut ), &stat );
  }

  /* Since the ``if'' branch freed the input and output time series,
     this needs to be repeated even if output was not requested.  (The
     reason for not collecting these operation after the two branches
     rejoined is that I wanted to keep as much memory free as possible
     when creating the power spectra, in case this program is run with
     huge datasets.) */
  else {
    SUB( LALSDestroyVector( &stat, &(input.data) ), &stat );
    SUB( LALSDestroyVector( &stat, &(output.data) ), &stat );
  }

  /* Test program executed successfully!  (User will have to examine
     the output file, if any, to determine that the algorithm achieved
     the correct result.) */
  LALCheckMemoryLeaks();
  INFO( RESAMPLETESTC_MSGENORM );
  return RESAMPLETESTC_ENORM;
}


REAL8
TDiff( REAL8 t )
{
  return df_f*sin( twoPiFm*t )/twoPiFm;
}


REAL8
DTDiff( REAL8 t )
{
  return df_f*cos( twoPiFm*t );
}


REAL8
DDTDiff( REAL8 t )
{
  return -df_f*twoPiFm*sin( twoPiFm*t );
}


REAL8
DDDTDiffMax( void )
{
  return df_f*twoPiFm*twoPiFm;
}
