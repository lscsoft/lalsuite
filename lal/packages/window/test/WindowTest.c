/*
*  Copyright (C) 2007 Bruce Allen, Duncan Brown, Jolien Creighton, Kipp Cannon, Teviet Creighton
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

/******************************** <lalVerbatim file="WindowTestCV">
Authora:Allen, B., Brown, D. A., and Creighton, T.
Revision: $Id$
**************************************************** </lalVerbatim> */
/********************************************************** <lalLaTeX>
\subsection{Program \texttt{WindowTest.c}}
\label{s:WindowTest.c}


Tests the routines in \verb@Window.h@.


\subsubsection*{Usage}
\begin{verbatim}
WindowTest [-d debuglevel] [-w flag] [-p] [-b beta] [-i infile] [-n width npts]
\end{verbatim}

\subsubsection*{Description}

This program writes output files \verb@PrintVector.@$nnn$, where $nnn$
is an integer from 0 to \verb@NumberWindowTypes@$-1$, each containing
the corresponding window function (along with a header of metadata
lines beginning with \verb@#@).  The following option flags are
accepted:
\begin{itemize}
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@ (default
0).
\item[\texttt{-w}] Restricts which window is used according to the
integer \verb@flag@, which is interpreted as a bit mask for turning on
each individual window type: the value of \verb@flag@ is a sum of
powers $2^{nnn}$, where $nnn$ is the number of the desired window in
the enumeration of window types (starting with rectangular at
$nnn=0$).  The default behaviour is equivalent to \verb@-w -1@, i.e.\
all available windows are used.  \textbf{NOTE:} This will need to
change if the number of windows grows to more than 31.
\item[\texttt{-p}] Switches to double-precision windowing.
\item[\texttt{-b}] Sets the window-dependent shape parameter $\beta$
equal to \verb@beta@.  By default, $\beta=6$ for a Kaiser window, 2
for a Creighton window, and is ignored by other windows.
\item[\texttt{-i}] Instead of simply printing out the window function,
the window is \emph{applied} to data read from \verb@infile@, which is
read using \verb@LALSReadSequence()@ or \verb@LALDReadSequence()@.
\item[\texttt{-n}] Sets the width of the window to \verb@width@, and
zero-pads the output out to a total length of \verb@npts@.  This is
useful to generate oversampled frequency spectra such as in
Fig.~\ref{fig:window-pectra}.  If not specified, the width and length
are taken to be 1024, or the length of data read with the \verb@-i@
option.
\end{itemize}
The output files simply columns of window(ed) data.  The window files
can be viewed for example by using the public domain graphing program
\verb@xmgr@ by typing:
\begin{verbatim}
xmgr PrintVector.*
\end{verbatim}

\verb@WindowTest@ also tests all error conditions, and, if the default
window size (1024 points) and $\beta$ parameter (6 for Kaiser, 2 for
Creighton) are used, checks that the sum of these windows squared add
to the correct values.  If there is an error in execution, the
corresponding error message is printed.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define WINDOWTESTC_ENORM 0
#define WINDOWTESTC_ESUB  1
#define WINDOWTESTC_EARG  2
#define WINDOWTESTC_ETEST 3
#define WINDOWTESTC_EFILE 4

#define WINDOWTESTC_MSGENORM "Normal exit"
#define WINDOWTESTC_MSGESUB  "Subroutine failed"
#define WINDOWTESTC_MSGEARG  "Error parsing arguments"
#define WINDOWTESTC_MSGETEST "Window creation function failed a test"
#define WINDOWTESTC_MSGEFILE "Could not open file"
/******************************************** </lalErrTable><lalLaTeX>

\vfill{\footnotesize\input{WindowTestCV}}

********************************************************** </lalLaTeX> */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/Window.h>
#include <lal/PrintVector.h>
#include <lal/StreamInput.h>
#include <lal/LALStdio.h>

NRCSID ( WINDOWTESTC, "$Id$");

/* modify this value to get a stack trace of test */
int lalDebugLevel = 0;

/* modify this to turn on printing windows into files for checking */
#define PRINT 1

/* default number of points */
#define NPOINTS 1024

/* Usage format string. */
#define USAGE "Usage: %s [-d debuglevel] [-w flag] [-p] [-i infile] [-n width npts]"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
                   __LINE__, WINDOWTESTC, statement ? statement :    \
                   "", (msg) );                                      \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( WINDOWTESTC_ESUB, WINDOWTESTC_MSGESUB,                    \
           "Function call \"" #func "\" failed:" );                  \
    return WINDOWTESTC_ESUB;                                         \
  }                                                                  \
} while (0)



int
main( int argc, char **argv )
{
  int arg = 1;                 /* counter over input arguments */
  FILE *fp;                    /* input/output file pointer */
  BOOLEAN pOption = 0;         /* whether to use REAL8 windows */
  BOOLEAN bOption = 0;         /* whether beta was user-specified */
  CHAR *infile = NULL;         /* input filename */
  CHAR outfile[LALNameLength]; /* output filename */
  UINT4 i;                     /* index over data */
  UINT4 npts = NPOINTS;        /* number of points in output */
  UINT4 width = NPOINTS;       /* width of window function */
  UINT4 flag = (UINT4)( -1 );  /* mask of windows to apply */
  REAL8 beta;                  /* window shape parameter */
  static LALStatus status;     /* top-level status structure */
  REAL4Vector *sVector = NULL; /* input data vector (single-precision) */
  REAL8Vector *dVector = NULL; /* input data vector (double-precision) */
  REAL4Window *sWindow = NULL; /* window to apply (single-precision) */
  REAL8Window *dWindow = NULL; /* window to apply (double-precision) */
  LALWindowParams params;      /* window creation parameters */
  WindowType wintype;          /* window type */
  REAL8 testsquares[] =        /* sum of squares for NPOINTS=1024: */
  { 1024.0,                /* rectangular */
    384.0,                 /* Hann */
    546.0+2.0/15.0,        /* Welch */
    341.333984375,         /* Bartlett */
    276.1142857152779,     /* Parzen */
    300.357781729967622,   /* Papoulis */
    406.9376,              /* Hamming */
    375.544713725875234,   /* Kaiser */
    393.028878331734330};  /* Creighton */


  /* Read command line arguments. */
  while ( arg < argc ) {

    /* Parse debug level option. */
    if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      } else {
	ERROR( WINDOWTESTC_EARG, WINDOWTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return WINDOWTESTC_EARG;
      }
    }

    /* Parse window flags option. */
    else if ( !strcmp( argv[arg], "-w" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	flag = (UINT4)( atoi( argv[arg++] ) );
      } else {
	ERROR( WINDOWTESTC_EARG, WINDOWTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return WINDOWTESTC_EARG;
      }
    }

    /* Parse precision option. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      arg++;
      pOption = 1;
    }

    /* Parse window size option. */
    else if ( !strcmp( argv[arg], "-b" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	beta = atof( argv[arg++] );
	bOption = 1;
      } else {
	ERROR( WINDOWTESTC_EARG, WINDOWTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return WINDOWTESTC_EARG;
      }
    }

    /* Parse window size option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	infile = argv[arg++];
      } else {
	ERROR( WINDOWTESTC_EARG, WINDOWTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return WINDOWTESTC_EARG;
      }
    }

    /* Parse window size option. */
    else if ( !strcmp( argv[arg], "-n" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	width = atoi( argv[arg++] );
	npts = atoi( argv[arg++] );
      } else {
	ERROR( WINDOWTESTC_EARG, WINDOWTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return WINDOWTESTC_EARG;
      }
    }

    /* Check for unrecognized arguments. */
    else {
      ERROR( WINDOWTESTC_EARG, WINDOWTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return WINDOWTESTC_EARG;
    }
  }

  /* Do argument-list-level tests. */
#ifndef LAL_NDEBUG
  if ( ! lalNoDebug ) {
    params.type = Rectangular;
    params.length = NPOINTS;
    params.beta = 1.0;

    /* Test behavior for non-positive length  */
    params.length=0;
    sWindow = XLALCreateREAL4Window(params.length, params.type, params.beta);
    if (sWindow) {
      ERROR(0, "XLALCreateREAL4Window() failed to detect 0 length", NULL);
      return 1;
    }

    /* Test failures for undefined window type on lower and upper bounds */
    params.type=-1;
    sWindow = XLALCreateREAL4Window(params.length, params.type, params.beta);
    if (sWindow) {
      ERROR(0, "XLALCreateREAL4Window() failed to detect out-of-bound type", NULL);
      return 1;
    }
    params.type=NumberWindowTypes;
    sWindow = XLALCreateREAL4Window(params.length, params.type, params.beta);
    if (sWindow) {
      ERROR(0, "XLALCreateREAL4Window() failed to detect out-of-bound type", NULL);
      return 1;
    }
  }
#endif


  /* Read input file, if any. */
  if ( infile ) {
    FILE *fp = fopen( infile, "r" );
    if ( !fp ) {
      ERROR( WINDOWTESTC_EFILE, WINDOWTESTC_MSGEFILE, infile );
      return WINDOWTESTC_EFILE;
    }
    if ( pOption )
      SUB( LALDReadSequence( &status, &dVector, fp ), &status );
    else
      SUB( LALSReadSequence( &status, &sVector, fp ), &status );
    fclose( fp );
  }


  /* Compute (and apply) windows. */
  params.length = width;
  for ( wintype = 0; wintype < NumberWindowTypes; wintype++ ) {
    if ( flag & ( 1 << wintype ) ) {

      /* Create window. */
      params.type = wintype;
      if ( !bOption ) {
	if ( wintype == Kaiser )
	  params.beta = 6.0;
	else if ( wintype == Creighton )
	  params.beta = 2.0;
      }
      if ( pOption ) {
      	dWindow = XLALCreateREAL8Window(params.length, params.type, params.beta);
	params.sumofsquares = dWindow->sumofsquares;
	params.windowname = dWindow->windowname;
      } else {
      	sWindow = XLALCreateREAL4Window(params.length, params.type, params.beta);
	params.sumofsquares = sWindow->sumofsquares;
	params.windowname = sWindow->windowname;
      }

      /* Check sum of squares. */
      if ( width == NPOINTS )
	if ( ( wintype != Kaiser || params.beta == 6.0 ) &&
	     ( wintype != Creighton || params.beta == 2.0 ) &&
	     ( wintype != Tukey ) ) {
	  if ( fabs( params.sumofsquares - testsquares[(int)wintype] )
	       > 1.e-5 ) {
	    printf("FAIL: Window %s appears incorrect.\n",
		   params.windowname );
	    printf("Expected %16.12f, got %16.12f\n",
		   testsquares[(int)wintype], params.sumofsquares );
	    return 1;
	  }
	}

      /* Apply/print window. */
      LALSnprintf( outfile, LALNameLength, "PrintVector.%03u",
		   (UINT4)( wintype ) );
      if ( !( fp = fopen( outfile, "w" ) ) ) {
	ERROR( WINDOWTESTC_EFILE, WINDOWTESTC_MSGEFILE, outfile );
	return WINDOWTESTC_EFILE;
      }

      /* Double-precision window: */
      if ( pOption ) {
	if ( infile ) {
	  for ( i = 0; i < dVector->length && i < width; i++ )
	    fprintf( fp, "%23.16e\n",
		     dVector->data[i]*dWindow->data->data[i] );
	  if ( i++ < npts && i < dVector->length )
	    fprintf( fp, "%23.16e\n",
		     dVector->data[i]*dWindow->data->data[0] );
	  SUB( LALDDestroyVector( &status, &dVector ), &status );
	} else {
	  for ( i = 0; i < width; i++ )
	    fprintf( fp, "%23.16e\n", dWindow->data->data[i] );
	  if ( i++ < npts )
	    fprintf( fp, "%23.16e\n", dWindow->data->data[0] );
	}
	XLALDestroyREAL8Window(dWindow);
	for ( ; i < npts; i++ )
	  fprintf( fp, "%23.16e\n", 0.0 );
      }

      /* Single-precision window: */
      else {
	if ( infile ) {
	  for ( i = 0; i < sVector->length && i < width; i++ )
	    fprintf( fp, "%16.9e\n",
		     sVector->data[i]*sWindow->data->data[i] );
	  if ( i++ < npts && i < sVector->length )
	    fprintf( fp, "%16.9e\n",
		     sVector->data[i]*sWindow->data->data[0] );
	  SUB( LALSDestroyVector( &status, &sVector ), &status );
	} else {
	  for ( i = 0; i < width; i++ )
	    fprintf( fp, "%16.9e\n", sWindow->data->data[i] );
	  if ( i++ < npts )
	    fprintf( fp, "%16.9e\n", sWindow->data->data[0] );
	}
	XLALDestroyREAL4Window(sWindow);
	for ( ; i < npts; i++ )
	  fprintf( fp, "%16.9e\n", 0.0 );
      }

      /* Done. */
      fclose( fp );
    }
  }
  return 0;
}
