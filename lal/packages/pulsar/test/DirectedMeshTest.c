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

/***************************** <lalVerbatim file="DirectedMeshTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{DirectedMeshTest.c}}
\label{ss:DirectedMeshTest.c}

Computes the sky-position metric for a coherent or semicoherent pulsar
search.

\subsubsection*{Usage}
\begin{verbatim}
DirectedMeshTest [-o outfile] [-d debuglevel] [-p n dt t0 f0] [-l lat lon]
                 [-s ra dec] [-r dra ddec] [-t tau] [-m mismatch]
\end{verbatim}

\subsubsection*{Description}

This test program computes template meshes for directed pulsar
searches with spindown, where it is assumed that the parameter metric
is constant over the search space.  The following option flags are
accepted:
\begin{itemize}
\item[\texttt{-o}] Prints the template mesh to the file
\verb@outfile@: each line consists of a sequence of
whitespace-separated numbers representing the coordinates of the
template.  If absent, the routines are exercised, but no output is
written.
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@; if
absent, \verb@-d 0@ is assumed.
\item[\texttt{-p}] Sets the search parameters: the number of stacks
\verb@n@, the length of each stack \verb@dt@ (in seconds), and the
start time of the first stack \verb@t0@ (in seconds of GPS time), and
the maximum source frequency \verb@f0@ (in Hz).  If absent,
\verb@-t 1 86400 0 1000@ is assumed.
\item[\texttt{-l}] Sets the detector latitude to \verb@lat@ (in
degrees north from the equator) and longitude to \verb@lon@ (in
degrees east of the prime meridian).  If absent,
\verb@-l 52.247 9.822@ (GEO600) is assumed.
\item[\texttt{-s}] Sets the right ascension and declination of the
target to \verb@ra@ and \verb@dec@ degrees, respectively.  If absent,
\verb@-r 192.8594813 27.1282511@ is assumed (the Galactic core).
\item[\texttt{-r}] Sets the range of the sky search to $\pm$\verb@dra@
degrees in right ascension and $\pm$\verb@ddec@ degrees in declination
about the target point.  If absent, \verb@-s 0 0@ is assumed (no
search over sky position).
\item[\texttt{-t}] Sets the range of the spindown search according to
the spindown timescale \verb@tau@, in seconds: the spindown parameter
$f_k$ is constrained by $|f_k|\leq$\verb@tau@${}^{-k}$.  If absent,
\verb@-t 3.16e9@ (century-long spindown) is assumed.
\item[\texttt{-m}] Sets the maximum mismatch threshold of the mesh to
\verb@mismatch@.  If absent, \verb@-m 0.25@ is assumed.
\end{itemize}

The program automatically determines how many spindown terms are
required to cover the parameter space, by starting with none (or one
if the search is \emph{only} over spindown), computing the local
template density from the parameter metric, and estimating the number
of templates required to cover the search volume.  The number of
spindown parameters is then increased and the number of templates
re-estimated.  Eventually the estimated number of templates will start
to \emph{decrease}, as the proper width of the parameter space in the
new dimensions is less than one template width.  The last dimension
before that happens is the correct dimension to use.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define DIRECTEDMESHTESTC_ENORM 0
#define DIRECTEDMESHTESTC_ESUB  1
#define DIRECTEDMESHTESTC_EARG  2
#define DIRECTEDMESHTESTC_EVAL  3
#define DIRECTEDMESHTESTC_EMEM  4
#define DIRECTEDMESHTESTC_EDET  5
#define DIRECTEDMESHTESTC_EFILE 6

#define DIRECTEDMESHTESTC_MSGENORM "Normal exit"
#define DIRECTEDMESHTESTC_MSGESUB  "Subroutine failed"
#define DIRECTEDMESHTESTC_MSGEARG  "Error parsing arguments"
#define DIRECTEDMESHTESTC_MSGEVAL  "Input argument out of valid range"
#define DIRECTEDMESHTESTC_MSGEMEM  "Memory allocation error"
#define DIRECTEDMESHTESTC_MSGEDET  "Non-positive metric determinant"
#define DIRECTEDMESHTESTC_MSGEFILE "Could not open file"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

The program is fairly straightforward.  It uses
\verb@LALStackMetric()@ to compute the parameter metric, passing it a
canonical time function \verb@LALDTComp()@ that composites the
\verb@LALDTSpin()@ and \verb@LALDTBaryPtolemaic()@ time
transformations.  It starts with a single spindown parameter.  (If a
\emph{search} over sky position is also indicated, it will start with
\emph{no} spindown parameters, using only the barycentric time
transformation rather than the composite transformation.)  The
determinant of the metric is computed, and the number of patches
estimated as:
\begin{equation}
N_\mathrm{patches} \approx \frac{\sqrt{|\mathsf{g}_{ab}|}}
	{(\mathrm{mismatch}/n)^{n/2}}
	\left\{\Delta\alpha\Delta\delta\right\}
	\tau^{-n_s(n_s+1)/2}
\end{equation}
where $\mathsf{g}_{ab}$ is the parameter metric (spindown sector only
if the sky search space is a single point), $n$ is the number of
dimensions in the search, $n_s$ is the number of spindown terms, and
$\Delta\alpha$ and $\Delta\delta$ are the half-ranges of the search in
right ascension and declination, respectively (assuming these ranges
are nonzero).  For the first round we are considering a search
\emph{only} over sky position ($n=2$, $n_s=0$) or a single spindown
search ($n=n_s=1$, ignore the term in braces).  The determinant is
computed using \verb@LALDMatrixDeterminantErr()@, repacking into
\verb@REAL8Array@s the metric components and uncertainties returned by
\verb@LALStackMetric()@.  An error is generated if the determinant is
non-positive, or a warning if it is smaller than its estimated
uncertainty.

In subsequent trials, we increase $n_s$ successivlely by 1, and
recompute $\mathsf{g}_{ab}$ and $N_\mathrm{patches}$.  Eventually,
when the width of the added dimension is less than one patch witdh,
$N_\mathrm{patches}$ will decrease.  When this happend, we back up to
the value of $n_s$ that gave the largest number of patches, and use
that parameter metric.

The program then uses \verb@LALDSymmetricEigenVectors()@ to compute
the eigenvalues and eigenvectors of the metric; these are combined and
repacked into a \verb@REAL4VectorSequence@ used by
\verb@LALFlatMesh()@, as described in \verb@FlatMesh.h@.  The inverse
transformation is computed using \verb@LALDMatrixInverse()@, and again
repacked into a \verb@REAL4VectorSequence@.  The search area is taken
to be a rectangular space controled by \verb@LALRectIntersect()@,
covering the sky area $|\alpha-\mathtt{ra}|\leq\mathtt{dra}$ and
$|\alpha-\mathtt{ra}|\leq\mathtt{dra}$ (provided these ranges are
nonzero), and the spindown volume $|f_k|\leq\mathtt{tau}^{-k}$ for
$k=1,\ldots,n_s$.  The volume boundaries are increased by half the
maximum patch size in each direction, to ensure total coverage of the
edges, as described in \verb@FlatMeshTest.c@.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALCalloc()                     LALFree()
LALU4CreateVector()             LALU4DestroyVector()
LALSCreateVector()              LALSDestroyVector()
LALDCreateVector()              LALDDestroyVector()
LALSCreateVectorSequence()      LALSDestroyVectorSequence()
LALDCreateArray()               LALDDestroyArray()
LALDTBaryPtolemaic()            LALTBaryPtolemaic()
LALDTSpin()                     LALTSpin()
LALDTComp()                     LALGetEarthTimes()
LALCreateFlatMesh()             LALRectIntersect()
LALStackMetric()                LALProjectMetric()
LALDMatrixDeterminantErr()      LALDMatrixInverse()
LALDSymmetricEigenVectors()     snprintf()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DirectedMeshTestCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/MatrixUtils.h>
#include <lal/StackMetric.h>
#include <lal/PulsarTimes.h>
#include <lal/FlatMesh.h>

NRCSID( DIRECTEDMESHTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;
#define NSTACKS 1
#define STACKLENGTH 86400.0 /* arbitrary */
#define STARTTIME 0.0       /* arbitrary */
#define LATITUDE  52.247    /* GEO600 location */
#define LONGITUDE 9.822     /* GEO600 location */
#define FREQUENCY 1000.0    /* arbitrary */
#define RA 192.8594813      /* Galactic core */
#define DEC 27.1282511      /* Galactic core */
#define TAU 3.16e9          /* century-scale spindown */
#define MISMATCH 0.25       /* arbitrary but reasonably optimal */

/* Usage format string. */
#define USAGE "Usage: %s [-o outfile] [-d debuglevel] [-p n dt t0 f0]\n" \
"\t[-l lat lon] [-s ra dec] [-r dra ddec] [-t tau] [-m mismatch]\n"

/* Input error checking: Some accepted parameter ranges. */
#define NMAX  10000 /* 1 <= number of stacks <= NMAX */
#define DTMAX  3e10 /* 1/f_0 < stack length <= DTMAX */
#define F0MAX  1e4  /* 0 < f_0 <= FOMAX */
/* Also: |latitude| and |dec| will be restricted to <= 90 degrees,
         |longitude| and |ra| will be restricted to <= 360 degrees */

/* Other internal constants. */
#define MAXLEN 1024 /* maximum format or warning message length */


/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, DIRECTEDMESHTESTC, statement ? statement :\
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
                 "        %s\n", *argv, __FILE__, __LINE__,          \
                 DIRECTEDMESHTESTC, (statement) );                   \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 DIRECTEDMESHTESTC, (statement) );                   \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( DIRECTEDMESHTESTC_ESUB, DIRECTEDMESHTESTC_MSGESUB,          \
         "Function call \"" #func "\" failed:" );                    \
  return DIRECTEDMESHTESTC_ESUB;                                     \
}                                                                    \
while (0)

#define CHECKVAL( val, lower, upper )                                \
do                                                                   \
if ( ( (val) <= (lower) ) || ( (val) > (upper) ) )                   \
{                                                                    \
  ERROR( DIRECTEDMESHTESTC_EVAL, DIRECTEDMESHTESTC_MSGESUB,          \
         "Value of " #val " out of range:" );                        \
  LALPrintError( #val " = %f, range = (%f,%f]\n", (REAL8)(val),      \
                 (REAL8)(lower), (REAL8)(upper) );                   \
  return DIRECTEDMESHTESTC_EVAL;                                     \
}                                                                    \
while (0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

/* Prototype for a routine to print floating-point numbers with
   uncertainties. */
int
fprintderr( FILE *fp, REAL8 x, REAL8 dx );

int
main(int argc, char **argv)
{
  /* Variables for parsing command-line arguments. */
  INT4  arg;                   /* argument list counter */
  CHAR *outfile = NULL;        /* output filename */
  UINT4 n = NSTACKS;           /* number of stacks or dimensions */
  REAL8 dt = STACKLENGTH;      /* length of each stack (s) */
  REAL8 t0 = STARTTIME;        /* start time of first stack (GPS s) */
  REAL8 f0 = FREQUENCY;        /* maximum frequency of search */
  REAL8 lat = LATITUDE;        /* latitude of detector (degrees N) */
  REAL8 lon = LONGITUDE;       /* longitude of detector (degrees E) */
  REAL8 ra = RA, dec = DEC;    /* sky position (degrees) */
  REAL8 dra = 0.0, ddec = 0.0; /* sky position range (degrees) */
  REAL8 tau = TAU;             /* spindown timescale (s) */
  REAL8 mismatch = MISMATCH;   /* maximum mismatch threshold */

  /* Other variables. */
  UINT4 i, j, ij, ji;           /* indecies */
  UINT4 nSpin, nSky;            /* dimensions of search */
  UINT4 err = 1;                /* 1 if no errors, 2 if errors. */
  REAL8Vector *lambda = NULL;   /* parameter/eigenvalue vector */
  REAL8Vector *metric = NULL;   /* current metric as a REAL8Vector */
  REAL8Array *current = NULL;   /* current metric as a REAL8Array */
  REAL8Array *dMatrix = NULL;   /* current metric uncertainties */
  REAL8Array *previous = NULL;  /* previous metric as a REAL8Array */
  REAL8 *mData, *cData, *dData; /* (metric, current, dMatrix)->data */
  REAL4 *sData, *width;         /* data pointers for edge calculations */
  REAL8 det[2];                 /* metric determinant and uncertainty */
  REAL8 vol, np;                /* variables for volumes, patch numbers */
  UINT4Vector dimLength;        /* structure to create arrays */
  UINT4 dimData[2];             /* dimLength.data storage */
  CreateVectorSequenceIn in;    /* structure to create vector sequence */
  REAL4VectorSequence *mesh = NULL; /* template list */
  static LALStatus stat;            /* top-level status structure */
  static LIGOTimeGPS start;         /* GPS start time of first stack */
  static MetricParamStruc params;   /* metric computation parameters */
  static FlatMeshParamStruc meshParams; /* mesh creation parameters */
  static PulsarTimesParamStruc baryParams; /* barycentring parameters */
  static PulsarTimesParamStruc spinParams; /* spindown parameters */
  static PulsarTimesParamStruc compParams; /* composite parameters */

  /* Set up array creation structure. */
  dimLength.length = 2;
  dimLength.data = dimData;

  /******************************************************************
   * ARGUMENT PARSING                                               *
   ******************************************************************/

  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      } else {
	ERROR( DIRECTEDMESHTESTC_EARG, DIRECTEDMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return DIRECTEDMESHTESTC_EARG;
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
	ERROR( DIRECTEDMESHTESTC_EARG, DIRECTEDMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return DIRECTEDMESHTESTC_EARG;
      }
    }
    /* Parse detector position option. */
    else if ( !strcmp( argv[arg], "-l" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	lat = atof( argv[arg++] );
	lon = atof( argv[arg++] );
      } else {
	ERROR( DIRECTEDMESHTESTC_EARG, DIRECTEDMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return DIRECTEDMESHTESTC_EARG;
      }
    }
    /* Parse sky position option. */
    else if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	ra = atof( argv[arg++] );
	dec = atof( argv[arg++] );
      } else {
	ERROR( DIRECTEDMESHTESTC_EARG, DIRECTEDMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return DIRECTEDMESHTESTC_EARG;
      }
    }
    /* Parse sky range option. */
    else if ( !strcmp( argv[arg], "-r" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	dra = atof( argv[arg++] );
	ddec = atof( argv[arg++] );
      } else {
	ERROR( DIRECTEDMESHTESTC_EARG, DIRECTEDMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return DIRECTEDMESHTESTC_EARG;
      }
    }
    /* Parse spindown timescale option. */
    else if ( !strcmp( argv[arg], "-t" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	tau = atof( argv[arg++] );
      } else {
	ERROR( DIRECTEDMESHTESTC_EARG, DIRECTEDMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return DIRECTEDMESHTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      } else {
	ERROR( DIRECTEDMESHTESTC_EARG, DIRECTEDMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return DIRECTEDMESHTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else {
      ERROR( DIRECTEDMESHTESTC_EARG, DIRECTEDMESHTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return DIRECTEDMESHTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Do error trapping on input parameters. */
  if ( lalDebugLevel & LALERROR ) {
    CHECKVAL( n, 0, NMAX );
    CHECKVAL( f0, 0.0, F0MAX );
    CHECKVAL( dt, 1.0/f0, DTMAX );
    CHECKVAL( lat, -90.0, 90.0 );
    CHECKVAL( lon, -360.0, 360.0 );
    CHECKVAL( dec, -90.0, 90.0 );
    CHECKVAL( ra, -360.0, 360.0 );
    CHECKVAL( tau, 0.0, LAL_REAL8_MAX );
  }
  if ( tau < n*dt )
    WARNING( "Spindown timescale less than observation time!" );

  /* Convert degrees to radians. */
  lat *= LAL_PI_180;
  lon *= LAL_PI_180;
  ra *= LAL_PI_180;
  dec *= LAL_PI_180;
  dra *= LAL_PI_180;
  ddec *= LAL_PI_180;

  /******************************************************************
   * METRIC COMPUTATION                                             *
   ******************************************************************/

  /* Set up start time. */
  start.gpsSeconds = (INT4)t0;
  start.gpsNanoSeconds = (INT4)( 1.0e9*( t0 - start.gpsSeconds ) );
  t0 = 0.0;

  /* Set up constant parameters for barycentre transformation. */
  baryParams.epoch = start;
  baryParams.latitude = lat;
  baryParams.longitude = lon;
  SUB( LALGetEarthTimes( &stat, &baryParams ), &stat );

  /* Set up constant parameters for spindown transformation. */
  spinParams.epoch = start;
  spinParams.t0 = t0;

  /* Set up constant parameters for composed transformation. */
  compParams.epoch = start;
  compParams.t1 = LALTBaryPtolemaic;
  compParams.t2 = LALTSpin;
  compParams.dt1 = LALDTBaryPtolemaic;
  compParams.dt2 = LALDTSpin;
  compParams.constants1 = &baryParams;
  compParams.constants2 = &spinParams;
  compParams.nArgs = 2;

  /* Set up constant parameters for initial metric calculation. */
  if ( dra == 0.0 || ddec == 0.0 ) {
    params.dtCanon = LALDTComp;
    params.constants = &compParams;
    nSpin = 1;
    nSky = 0;
  } else {
    params.dtCanon = LALDTBaryPtolemaic;
    params.constants = &baryParams;
    nSpin = 0;
    nSky = 2;
  }
  params.start = t0;
  params.deltaT = dt;
  params.n = n;
  params.errors = 1;
  err = 2;

  /* Set up the initial metric and parameter vectors. */
  n = nSpin + 3;
  SUB( LALDCreateVector( &stat, &metric, err*n*(n+1)/2 ), &stat );
  mData = metric->data;
  SUB( LALDCreateVector( &stat, &lambda, n ), &stat );
  lambda->data[0] = f0;
  lambda->data[1] = ra;
  lambda->data[2] = dec;
  if ( nSpin )
    memset( lambda->data + 3, 0, nSpin*sizeof(REAL8) );

  /* Compute the initial metric and its projection. */
  SUB( LALStackMetric( &stat, metric, lambda, &params ), &stat );
  SUB( LALDDestroyVector( &stat, &lambda ), &stat );
  SUB( LALProjectMetric( &stat, metric, params.errors ), &stat );

  /* Copy appropriate metric components into a matrx. */
  dimData[0] = dimData[1] = n = nSpin + nSky;
  SUB( LALDCreateArray( &stat, &current, &dimLength ), &stat );
  cData = current->data;
  for ( i = 0; i < n; i++ ) {
    cData[i*(n+1)] = mData[err*(i+3-nSky)*(i+6-nSky)/2];
    for ( j = 0; j < i; j++ )
      cData[i*n+j] = cData[j*n+i]
	= mData[err*(j+3-nSky) + err*(i+3-nSky)*(i+4-nSky)/2];
  }
  if ( params.errors ) {
    SUB( LALDCreateArray( &stat, &dMatrix, &dimLength ), &stat );
    dData = dMatrix->data;
    for ( i = 0; i < n; i++ ) {
      dData[i*(n+1)] = mData[err*(i+3-nSky)*(i+6-nSky)/2 + 1];
      for ( j = 0; j < i; j++ )
	dData[i*n+j] = dData[j*n+i]
	  = mData[err*(j+3-nSky) + err*(i+3-nSky)*(i+4-nSky)/2 + 1];
    }
  }
  SUB( LALDDestroyVector( &stat, &metric ), &stat );


  for ( i = 0; i < n; i++ ) {
    fprintf( stderr, "%25.16e", cData[i*n] );
    for ( j = 1; j < n; j++ )
      fprintf( stderr, "% 25.16e", cData[i*n+j] );
    fprintf( stderr, "\n" );
  }



  /* Estimate the number of templates required. */
  SUB( LALDMatrixDeterminantErr( &stat, det, current, dMatrix ),
       &stat );
  if ( params.errors ) {
    SUB( LALDDestroyArray( &stat, &dMatrix ), &stat );
  }
  if ( det[0] <= 0.0 ) {
    ERROR( DIRECTEDMESHTESTC_EDET, DIRECTEDMESHTESTC_MSGEDET, 0 );
    return DIRECTEDMESHTESTC_EDET;
  }
  if ( det[0] <= det[1] )
    WARNING( "Uncertainty in determinant greater than value" );
  det[0] = sqrt( det[0] );
  det[1] /= 2.0*det[0];
  vol = pow( mismatch/( (REAL8)( n ) ), 0.5*( (REAL8)( n ) ) );
  vol *= pow( tau, -0.5*( (REAL8)( nSpin*( nSpin + 1 ) ) ) );
  if ( nSky )
    vol *= dra*ddec;
  det[0] *= vol;
  det[1] *= vol;
  fprintf( stdout, "Estimated number of templates:\n" );
  fprintf( stdout, "%i spindown: number of templates = ", nSpin );
  fprintderr( stdout, det[0], det[1] );
  fprintf( stdout, "\n" );

  /* Increase number of spindown terms, and re-estimate number of
     templates. */
  params.dtCanon = LALDTComp;
  params.constants = &compParams;
  do {
    nSpin++;

    /* Store old results. */
    np = det[0];
    if ( previous ) {
      SUB( LALDDestroyArray( &stat, &previous ), &stat );
    }
    previous = current;
    current = NULL;

    /* Set up the current metric and parameter vectors. */
    n = nSpin + 3;
    SUB( LALDCreateVector( &stat, &metric, err*n*(n+1)/2 ), &stat );
    mData = metric->data;
    SUB( LALDCreateVector( &stat, &lambda, n ), &stat );
    lambda->data[0] = f0;
    lambda->data[1] = ra;
    lambda->data[2] = dec;
    memset( lambda->data + 3, 0, nSpin*sizeof(REAL8) );

    /* Compute the current metric and its projection. */
    SUB( LALStackMetric( &stat, metric, lambda, &params ), &stat );
    SUB( LALDDestroyVector( &stat, &lambda ), &stat );
    SUB( LALProjectMetric( &stat, metric, params.errors ), &stat );

    for ( i = 0; i < metric->length/2; i++ )
      fprintf( stderr, "%25.16e\n", metric->data[2*i] );


    /* Copy appropriate metric components into a matrx. */
    dimData[0] = dimData[1] = n = nSpin + nSky;
    SUB( LALDCreateArray( &stat, &current, &dimLength ), &stat );
    cData = current->data;
    for ( i = 0; i < n; i++ ) {
      cData[i*(n+1)] = mData[err*(i+3-nSky)*(i+6-nSky)/2];
      for ( j = 0; j < i; j++ )
	cData[i*n+j] = cData[j*n+i]
	  = mData[err*(j+3-nSky) + err*(i+3-nSky)*(i+4-nSky)/2];
    }
    if ( params.errors ) {
      SUB( LALDCreateArray( &stat, &dMatrix, &dimLength ), &stat );
      dData = dMatrix->data;
      for ( i = 0; i < n; i++ ) {
	dData[i*(n+1)] = mData[err*(i+3-nSky)*(i+6-nSky)/2 + 1];
	for ( j = 0; j < i; j++ )
	  dData[i*n+j] = dData[j*n+i]
	    = mData[err*(j+3-nSky) + err*(i+3-nSky)*(i+4-nSky)/2 + 1];
      }
    }
    SUB( LALDDestroyVector( &stat, &metric ), &stat );


    fprintf( stderr, "\n" );
    for ( i = 0; i < n; i++ ) {
      fprintf( stderr, "%25.16e", cData[i*n] );
      for ( j = 1; j < n; j++ )
	fprintf( stderr, "% 25.16e", cData[i*n+j] );
      fprintf( stderr, "\n" );
    }



    /* Estimate the number of templates required. */
    SUB( LALDMatrixDeterminantErr( &stat, det, current, dMatrix ),
	 &stat );
    if ( params.errors ) {
      SUB( LALDDestroyArray( &stat, &dMatrix ), &stat );
    }
    if ( det[0] <= 0.0 ) {
      ERROR( DIRECTEDMESHTESTC_EDET, DIRECTEDMESHTESTC_MSGEDET, 0 );
      return DIRECTEDMESHTESTC_EDET;
    }
    if ( det[0] <= det[1] )
      WARNING( "Uncertainty in determinant greater than value" );
    det[0] = sqrt( det[0] );
    det[1] /= 2.0*det[0];
    vol = pow( mismatch/( (REAL8)( n ) ), 0.5*( (REAL8)( n ) ) );
    vol *= pow( tau, -0.5*( (REAL8)( nSpin*( nSpin + 1 ) ) ) );
    if ( nSky )
      vol *= dra*ddec;
    det[0] *= vol;
    det[1] *= vol;
    fprintf( stdout, "%i spindown: number of templates = ", nSpin );
    fprintderr( stdout, det[0], det[1] );
    fprintf( stdout, "\n" );
  } while ( det[0] > np );

  /******************************************************************
   * MESH PLACEMENT                                                 *
   ******************************************************************/

  /* Use next-to-last number of spindown terms. */
  nSpin--;
  dimData[0] = dimData[1] = n = nSpin + nSky;
  SUB( LALDDestroyArray( &stat, &current ), &stat );
  SUB( LALDCreateArray( &stat, &current, &dimLength ), &stat );
  SUB( LALDCreateArray( &stat, &dMatrix, &dimLength ), &stat );
  SUB( LALDCreateVector( &stat, &lambda, n ), &stat );

  /* Compute transformation matrices: current stores transformation
     matrix, dMatrix stores its inverse, previous gets mangled. */
  SUB( LALDSymmetricEigenVectors( &stat, lambda, previous ), &stat );
  memcpy( current->data, previous->data, n*n*sizeof(REAL8) );
  for ( i = 0; i < n; i++ ) {
    REAL8 factor = 2.0*mismatch*sqrt( n*lambda->data[i] );
    for ( j = 0, ji = j; j < n; j++, ji += n )
      current->data[ji] *= factor;
  }
  SUB( LALDDestroyVector( &stat, &lambda ), &stat );
  SUB( LALDMatrixInverse( &stat, det, previous, dMatrix ), &stat );
  SUB( LALDDestroyArray( &stat, &previous ), &stat );

  /* Create flat mesh parameter structure fields. */
  in.length = in.vectorLength = n;
  SUB( LALSCreateVectorSequence( &stat, &(meshParams.matrix), &in ),
       &stat );
  SUB( LALSCreateVectorSequence( &stat, &(meshParams.matrixInv),
				 &in ), &stat );
  SUB( LALSCreateVector( &stat, &(meshParams.xMax), n ), &stat );
  SUB( LALSCreateVector( &stat, &(meshParams.xMin), n ), &stat );
  in.length = 2;
  SUB( LALSCreateVectorSequence( &stat, &(meshParams.controlPoints),
				 &in ), &stat );

  /* Transpose matrices to agree with convention in FlatMesh.h. */
  for ( i = 0; i < n; i++ )
    for ( j = 0, ij = i*n, ji = j; j < n; j++, ij++, ji += n ) {
      meshParams.matrix->data[ij] = current->data[ji];
      meshParams.matrixInv->data[ij] = dMatrix->data[ji];
    }
  SUB( LALDDestroyArray( &stat, &current ), &stat );
  SUB( LALDDestroyArray( &stat, &dMatrix ), &stat );

  /* Set up boundaries. */
  sData = meshParams.controlPoints->data;
  if ( nSky ) {
    sData[0] = ra - dra;
    sData[n] = ra + dra;
    sData[1] = dec - ddec;
    sData[n+1] = dec + ddec;
  }
  vol = 1.0;
  for ( i = nSky; i < n; i++ ) {
    sData[i] = -( vol /= tau );
    sData[n+i] = vol;
  }

  /* Expand boundaries to ensure edge coverage. */
  width = (REAL4 *)LALCalloc( n, sizeof(REAL4) );
  if ( !width ) {
    ERROR( DIRECTEDMESHTESTC_EMEM, DIRECTEDMESHTESTC_MSGEMEM, 0 );
    return DIRECTEDMESHTESTC_EMEM;
  }
  for ( sData = meshParams.matrix->data, i = 0; i < n; i++ )
    for ( j = 0; j < n; j++, sData++ )
      width[j] += fabs( *sData );
  for ( sData = meshParams.controlPoints->data, i = 0; i < n;
	i++, sData++ ) {
    INT2 direct = ( sData[0] < sData[n] ) ? -1 : 1;
    sData[0] += 0.5*direct*width[i];
    sData[n] -= 0.5*direct*width[i];
  }
  LALFree( width );

  /* Copy final boundaries into mesh parameters. */
  memcpy( meshParams.xMin->data, meshParams.controlPoints->data,
	  n*sizeof(REAL4) );
  memcpy( meshParams.xMax->data, meshParams.controlPoints->data + n,
	  n*sizeof(REAL4) );

  /* Compute the mesh, and then clean up memory. */
  SUB( LALCreateFlatMesh( &stat, &mesh, &meshParams ), &stat );
  SUB( LALSDestroyVector( &stat, &(meshParams.xMax) ), &stat );
  SUB( LALSDestroyVector( &stat, &(meshParams.xMin) ), &stat );
  SUB( LALSDestroyVectorSequence( &stat, &(meshParams.matrix) ),
       &stat );
  SUB( LALSDestroyVectorSequence( &stat, &(meshParams.matrixInv) ),
       &stat );
  SUB( LALSDestroyVectorSequence( &stat, &(meshParams.controlPoints) ),
       &stat );

  /* Write output if requested. */
  if ( outfile ) {
    FILE *fp = fopen( outfile, "w" );
    if ( !fp ) {
      ERROR( DIRECTEDMESHTESTC_EFILE, "- " DIRECTEDMESHTESTC_MSGEFILE,
	     outfile );
      return DIRECTEDMESHTESTC_EFILE;
    }
    sData = mesh->data;
    i = mesh->length;
    while ( i-- ) {
      fprintf( fp, "%16.9e", *(sData++) );
      for ( j = 1; j < n; j++ )
	fprintf( fp, " %16.9e", *(sData++) );
      fprintf( fp, "\n" );
    }
    fclose( fp );
  }

  /* Done. */
  SUB( LALSDestroyVectorSequence( &stat, &mesh ), &stat );
  LALCheckMemoryLeaks();
  INFO( DIRECTEDMESHTESTC_MSGENORM );
  return DIRECTEDMESHTESTC_ENORM;
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
      snprintf( format, MAXLEN, "%%.%if +/- %%.%if", -lsd, -lsd );
      return fprintf( fp, format, 0.0, dx );
    }
    norm = pow( 10.0, -lsd );
    snprintf( format, MAXLEN, "0 +/- %%.0f%%0%ii", lsd );
    return fprintf( fp, format, dx*norm, 0 );
  }

  /* If number is exact to 8-byte precision, print it as such. */
  if ( dx <= 0.0 ) {
    if ( abs( gsd ) > 3 )
      return fprintf( fp, "%.16e", x );
    snprintf( format, MAXLEN, "%%.%if", 16 - gsd );
    return fprintf( fp, format, x );
  }

  /* Otherwise, format depends on x and dx. */
  if ( gsd < lsd )
    gsd = lsd;
  if ( lsd > 3 || gsd < -3 ) {
    norm = pow( 10.0, -gsd );
    snprintf( format, MAXLEN, "( %%.%if +/- %%.%if )e%+i",
		 gsd - lsd, gsd - lsd, gsd );
    return fprintf( fp, format, x*norm, dx*norm );
  }
  if ( lsd <= 0 ) {
    snprintf( format, MAXLEN, "%%.%if +/- %%.%if", -lsd, -lsd );
    return fprintf( fp, format, x, dx );
  }
  norm = pow( 10.0, -lsd );
  snprintf( format, MAXLEN, "%%.0f%%0%ii +/- %%.0f%%0%ii", lsd,
	       lsd );
  return fprintf( fp, format, x*norm, 0, dx*norm, 0 );
}
