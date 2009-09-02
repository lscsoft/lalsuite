/*
*  Copyright (C) 2007 John Whelan
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

/******************************** <lalVerbatim file="StochasticOmegaGWTestCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{StochasticOmegaGWTest.c}}
\label{stochastic:ss:StochasticOmegaGWTest.c}

A program to test \texttt{LALStochasticOmegaGW()}.

\subsubsection*{Usage}

\begin{verbatim}
./StochasticOmegaGWTest [options]
Options:
  -h             print usage message
  -q             quiet: run silently
  -v             verbose: print extra information
  -d level       set lalDebugLevel to level
  -a alpha       set power law exponent to alpha
  -O omegaRef    set amplitude to omegaRef
  -F fRef        set normalization reference frequency to fRef
  -f f0          set start frequency to f0
  -e deltaF      set frequency spacing to deltaF
  -n length      set number of points in frequency series to length
  -o filename    print gravitational-wave spectrum to file filename
\end{verbatim}

\subsubsection*{Description}

This program tests the function {\tt LALStochasticOmegaGW()\/}, which outputs a
power law spectrum
\begin{equation}
h_{100}^2\Omega_{\scriptstyle{\rm GW}}(f)
=\Omega_{\scriptstyle{\rm R}}
\left(
  \frac{f}{f_{\scriptstyle{\rm R}}}
\right)^\alpha
\end{equation}

First, it tests that the correct error codes
(\textit{cf.}\ Sec.~\ref{stochastic:s:StochasticCrossCorrelation.h})
are generated for the following error conditions (tests in
\textit{italics} are not performed if \verb+LAL_NDEBUG+ is set, as
the corresponding checks in the code are made using the ASSERT macro):
\begin{itemize}
\item \textit{null pointer to output series}
\item \textit{null pointer to parameter structure}
\item \textit{null pointer to data member of output series}
\item \textit{null pointer to data member of data member of output series}
\item \textit{zero length parameter}
\item \textit{negative frequency spacing}
\item \textit{zero frequency spacing}
\item mismatch between length of output series and length parameter
\item zero reference frequency $f_{\scriptstyle{\rm R}}$
% \item reference frequency $f_{\scriptstyle{\rm R}}$
% smaller than lowest positive output frequency
\item negative amplitude parameter $\Omega_{\scriptstyle{\rm R}}$
\item zero amplitude parameter $\Omega_{\scriptstyle{\rm R}}$
\end{itemize}

It then verifies that the correct frequency series are generated for
two simple test cases: $\alpha=2.5$ and $\alpha=0$.  For each
successful test (both of these valid data and the invalid ones
described above), it prints ``\texttt{PASS}'' to standard output; if a
test fails, it prints ``\texttt{FAIL}''.

If the \texttt{filename} argument is present, it also calculates a
spectrum based on user-specified data.
Figure~\ref{stochastic:f:quadOmega} illustrates the output of the
command with the following arguments:
\begin{verbatim}
StochasticOmegaGWTest -e 1 -n 1000 -F 100 -O 1e-6 -a 2 -o OmegaGW.dat
\end{verbatim}

\begin{figure}[htb!]
\begin{center}
\noindent
\includegraphics[width=4in,angle=-90]{stochasticOmegaGWQuadratic}
\caption{\label{stochastic:f:quadOmega}
A quadratic stochastic gravitational-wave background spectrum.}
\end{center}
\end{figure}

\subsubsection*{Exit codes}
\input{StochasticOmegaGWTestCE}

\subsubsection*{Uses}

\begin{verbatim}
lalDebugLevel
getopt()
LALSCreateVector()
LALStochasticOmegaGW()
LALSPrintFrequencySeries
LALSDestroyVector()
LALCheckMemoryLeaks()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}
  \item No specific error checking is done on user-specified data.  If
\texttt{deltaF} or \texttt{length} are missing, the resulting defaults
will cause a bad data error.  If other arguments are unspecified, the
following defaults are used:
\begin{description}
\item[\texttt{alpha}] 0
\item[\texttt{f0}] 0
\item[\texttt{fRef}] 1\,Hz
\item[\texttt{omegaRef}] 1
\end{description}
\item The routine \texttt{LALStochasticOmegaGW()} will eventually be generalized to
include ``broken'' power law spectra
\begin{equation}
h_{100}^2\Omega_{\scriptstyle{\rm GW}}
= \left\{
\begin{array}{cc}
\Omega_1 f^{\alpha_1} & f\le f_c\\
\Omega_2 f^{\alpha_2} & f\ge f_c
\end{array}
\right.
\end{equation}
\end{itemize}

\vfill{\footnotesize\input{StochasticOmegaGWTestCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <config.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/StochasticCrossCorrelation.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/Units.h>

#include "CheckStatus.h"

NRCSID(STOCHASTICOMEGAGWTESTC, "$Id$");

#define STOCHASTICOMEGAGWTESTC_OMEGAREF  1e-6
#define STOCHASTICOMEGAGWTESTC_FREF      100.0
#define STOCHASTICOMEGAGWTESTC_ALPHA     2.5
#define STOCHASTICOMEGAGWTESTC_LENGTH    8
#define STOCHASTICOMEGAGWTESTC_F0        0.0
#define STOCHASTICOMEGAGWTESTC_DELTAF    80.0
#define STOCHASTICOMEGAGWTESTC_TOL       1e-6

#define STOCHASTICOMEGAGWTESTC_TRUE     1
#define STOCHASTICOMEGAGWTESTC_FALSE    0

extern char *optarg;
extern int   optind;

/* int lalDebugLevel = LALMSGLVL3; */
int lalDebugLevel = LALNDEBUG;
BOOLEAN optVerbose    = STOCHASTICOMEGAGWTESTC_FALSE;
REAL8 optDeltaF     = -1.0;
UINT4 optLength     = 0;
REAL4 optAlpha    = 0.0;
REAL8 optF0       = 0.0;
REAL4 optFR       = 1.0;
REAL4 optOmegaR   = 1.0;
CHAR optFile[LALNameLength] = "";

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

/***************************** <lalErrTable file="StochasticOmegaGWTestCE"> */
#define STOCHASTICOMEGAGWTESTC_ENOM 0
#define STOCHASTICOMEGAGWTESTC_EARG 1
#define STOCHASTICOMEGAGWTESTC_ECHK 2
#define STOCHASTICOMEGAGWTESTC_EFLS 3
#define STOCHASTICOMEGAGWTESTC_EUSE 4
#define STOCHASTICOMEGAGWTESTC_MSGENOM "Nominal exit"
#define STOCHASTICOMEGAGWTESTC_MSGEARG "Error parsing command-line arguments"
#define STOCHASTICOMEGAGWTESTC_MSGECHK "Error checking failed to catch bad data"
#define STOCHASTICOMEGAGWTESTC_MSGEFLS "Incorrect answer for valid data"
#define STOCHASTICOMEGAGWTESTC_MSGEUSE "Bad user-entered data"
/***************************** </lalErrTable> */

int main( int argc, char *argv[] )
{
  static LALStatus                status;

  StochasticOmegaGWParameters   parameters;
  REAL4FrequencySeries     omegaGW;

  REAL4FrequencySeries     dummyOutput;

  REAL4                *tempPtr;

  UINT4 i;
  REAL4 omega, f;
  INT4 code;

  /* define valid parameters */

  parameters.alpha    = STOCHASTICOMEGAGWTESTC_ALPHA;
  parameters.fRef     = STOCHASTICOMEGAGWTESTC_FREF;
  parameters.omegaRef = STOCHASTICOMEGAGWTESTC_OMEGAREF;
  parameters.length   = STOCHASTICOMEGAGWTESTC_LENGTH;
  parameters.f0       = STOCHASTICOMEGAGWTESTC_F0;
  parameters.deltaF   = STOCHASTICOMEGAGWTESTC_DELTAF;

  omegaGW.data = NULL;

  dummyOutput.data = NULL;

  ParseOptions( argc, argv );

  LALSCreateVector(&status, &(omegaGW.data), STOCHASTICOMEGAGWTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICOMEGAGWTESTC_EFLS,
			    STOCHASTICOMEGAGWTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /* TEST INVALID DATA HERE ------------------------------------------ */
#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to real frequency series for output */
    LALStochasticOmegaGW(&status, NULL, &parameters);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOMEGAGWTESTC_ECHK,
			      STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to output series results in error:       \n\"%s\"\n",
	   STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);


    /* test behavior for null pointer to parameter structure */
    LALStochasticOmegaGW(&status, &omegaGW, NULL);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOMEGAGWTESTC_ECHK,
			      STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to parameter structure results in error:       \n\"%s\"\n",
	   STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);


    /* test behavior for null pointer to data member of real frequency
       series for output */
    LALStochasticOmegaGW(&status, &dummyOutput, &parameters);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOMEGAGWTESTC_ECHK,
			      STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of output series results in error:       \n\"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);


    /* Create a vector for testing null data-data pointer */
    LALSCreateVector(&status, &(dummyOutput.data), STOCHASTICOMEGAGWTESTC_LENGTH);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOMEGAGWTESTC_EFLS,
			      STOCHASTICOMEGAGWTESTC_MSGEFLS) ) )
    {
      return code;
    }
    tempPtr = dummyOutput.data->data;
    dummyOutput.data->data = NULL;

    /* test behavior for null pointer to data member of data member of
       real frequency series for output */
    LALStochasticOmegaGW(&status, &dummyOutput, &parameters);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICOMEGAGWTESTC_ECHK,
			      STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of output series results in error:       \n\"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* clean up */

    dummyOutput.data->data = tempPtr;
    LALSDestroyVector(&status, &(dummyOutput.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOMEGAGWTESTC_EFLS,
			      STOCHASTICOMEGAGWTESTC_MSGEFLS) ) )
    {
      return code;
    }

    /* test behavior for length parameter equal to zero */
    parameters.length = 0;
    LALStochasticOmegaGW(&status, &omegaGW, &parameters);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EZEROLEN,
			      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN,
			      STOCHASTICOMEGAGWTESTC_ECHK,
			      STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero length parameter results in error:       \n\"%s\"\n",
	   STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);
    /* assign valid length parameter */
    parameters.length = STOCHASTICOMEGAGWTESTC_LENGTH;

    /* test behavior for frequency spacing less than or equal to zero */
    parameters.deltaF = -1;
    LALStochasticOmegaGW(&status, &omegaGW, &parameters);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
			      STOCHASTICOMEGAGWTESTC_ECHK,
			      STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: negative frequency spacing results in error:       \n\"%s\"\n",
	   STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

    parameters.deltaF = 0;
    LALStochasticOmegaGW(&status, &omegaGW, &parameters);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
			      STOCHASTICOMEGAGWTESTC_ECHK,
			      STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
    {
        return code;
    }
    printf("  PASS: zero frequency spacing results in error:       \n\"%s\"\n",
	   STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* assign valid frequency spacing */
    parameters.deltaF = STOCHASTICOMEGAGWTESTC_DELTAF;
  }

#endif /* LAL_NDEBUG */

  /* test behavior for negative start frequency */
  parameters.f0 = -20.0;
  LALStochasticOmegaGW(&status, &omegaGW, &parameters);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN,
			    STOCHASTICOMEGAGWTESTC_ECHK,
			    STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: negative start frequency results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);

  /* reassign valid start frequency */
  parameters.f0 = STOCHASTICOMEGAGWTESTC_F0;

  /* test behavior for length of data member of real frequency series
     for output not equal to length specified in input parameters */
  parameters.length += 1;
  LALStochasticOmegaGW(&status, &omegaGW, &parameters);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICOMEGAGWTESTC_ECHK,
			    STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: mismatch between length of output series and length parameter results in error:       \n\"%s\"\n",
	 STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  /* reassign valid length to data member of dummy output */
  parameters.length -= 1;

  /* test behavior for fRef < deltaf */
  parameters.fRef = 0;
  LALStochasticOmegaGW(&status, &omegaGW, &parameters);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EOORFREF,
                          STOCHASTICCROSSCORRELATIONH_MSGEOORFREF,
			  STOCHASTICOMEGAGWTESTC_ECHK,
			    STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: zero reference frequency results in error:       \n\"%s\"\n",
	 STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  parameters.fRef = STOCHASTICOMEGAGWTESTC_FREF;

  /* test behavior for omegaRef <=0 */
  parameters.omegaRef = -1.0;
  LALStochasticOmegaGW(&status, &omegaGW, &parameters);
  if ( ( code = CheckStatus(&status,  STOCHASTICCROSSCORRELATIONH_ENONPOSOMEGA,
                          STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA,
			  STOCHASTICOMEGAGWTESTC_ECHK,
			    STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: negative amplitude parameter results in error:       \n\"%s\"\n",
	 STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA);

  parameters.omegaRef = 0.0;
  LALStochasticOmegaGW(&status, &omegaGW, &parameters);
  if ( ( code = CheckStatus(&status,  STOCHASTICCROSSCORRELATIONH_ENONPOSOMEGA,
			    STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA,
			    STOCHASTICOMEGAGWTESTC_ECHK,
			    STOCHASTICOMEGAGWTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: zero amplitude parameter results in error:       \n\"%s\"\n",
	 STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA);
  parameters.omegaRef = STOCHASTICOMEGAGWTESTC_OMEGAREF;

  /* TEST VALID DATA HERE -------------------------------------------- */

  /* generate omegaGW */
  LALStochasticOmegaGW(&status, &omegaGW, &parameters);
  if ( ( code = CheckStatus(&status,0, "",
			    STOCHASTICOMEGAGWTESTC_EFLS,
			    STOCHASTICOMEGAGWTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /* test values */

  for (i=0; i<STOCHASTICOMEGAGWTESTC_LENGTH; ++i)
  {
    f = i * STOCHASTICOMEGAGWTESTC_DELTAF;
    omega = STOCHASTICOMEGAGWTESTC_OMEGAREF
      * pow(f/STOCHASTICOMEGAGWTESTC_FREF,STOCHASTICOMEGAGWTESTC_ALPHA);
    if (optVerbose)
    {
      printf("Omega(%f Hz)=%g, should be %g\n",
             f, omegaGW.data->data[i], omega);
    }
    if ( (omegaGW.data->data[i] - omega) &&
         abs((omegaGW.data->data[i] - omega)/omega) > STOCHASTICOMEGAGWTESTC_TOL )
    {
      printf("  FAIL: Valid data test #1 (alpha=%f)\n",STOCHASTICOMEGAGWTESTC_ALPHA);
      return STOCHASTICOMEGAGWTESTC_EFLS;
    }
  }
  printf("  PASS: Valid data test #1 (alpha=%f)\n",STOCHASTICOMEGAGWTESTC_ALPHA);

  /* change parameters */
  parameters.alpha = 0.0;

  /* generate omegaGW */
  LALStochasticOmegaGW(&status, &omegaGW, &parameters);
  if ( ( code = CheckStatus(&status,0, "",
			    STOCHASTICOMEGAGWTESTC_EFLS,
			    STOCHASTICOMEGAGWTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /* test values */

  for (i=0; i<STOCHASTICOMEGAGWTESTC_LENGTH; ++i)
  {
    f = i * STOCHASTICOMEGAGWTESTC_DELTAF;
    if (optVerbose) {
      printf("Omega(%f Hz)=%g, should be %g\n",
             f, omegaGW.data->data[i], STOCHASTICOMEGAGWTESTC_OMEGAREF);
    }
    if ( abs(omegaGW.data->data[i] - STOCHASTICOMEGAGWTESTC_OMEGAREF)
         / STOCHASTICOMEGAGWTESTC_OMEGAREF > STOCHASTICOMEGAGWTESTC_TOL )
    {
      printf("  FAIL: Valid data test #2 (alpha=0)\n");
      return STOCHASTICOMEGAGWTESTC_EFLS;
    }
  }
  printf("  PASS: Valid data test #2 (alpha=0)\n");

  /* clean up valid data */
  LALSDestroyVector(&status, &(omegaGW.data));
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICOMEGAGWTESTC_EFLS,
			    STOCHASTICOMEGAGWTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCheckMemoryLeaks();

  printf("PASS: all tests\n");

  if (optFile[0]) {
    parameters.alpha = optAlpha;
    parameters.length = optLength;
    parameters.deltaF = optDeltaF;
    parameters.f0 = optF0;
    parameters.omegaRef = optOmegaR;
    parameters.fRef = optFR;
    LALSCreateVector(&status, &(omegaGW.data), optLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOMEGAGWTESTC_EUSE,
			      STOCHASTICOMEGAGWTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALStochasticOmegaGW(&status, &omegaGW, &parameters);
    if ( ( code = CheckStatus(&status,0, "",
			      STOCHASTICOMEGAGWTESTC_EUSE,
			      STOCHASTICOMEGAGWTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALSPrintFrequencySeries( &omegaGW, optFile );

    printf("=== Stochastic Gravitational-wave Spectrum Written to File %s ===\n", optFile);


    LALSDestroyVector(&status, &(omegaGW.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICOMEGAGWTESTC_EUSE,
			      STOCHASTICOMEGAGWTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCheckMemoryLeaks();
  }

  return STOCHASTICOMEGAGWTESTC_ENOM;
}

/*
 * Usage ()
 *
 * Prints a usage message for program program and exits with code exitcode.
 *
 */
static void
Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h             print this message\n");
  fprintf (stderr, "  -q             quiet: run silently\n");
  fprintf (stderr, "  -v             verbose: print extra information\n");
  fprintf (stderr, "  -d level       set lalDebugLevel to level\n");
  fprintf (stderr, "  -a alpha       set power law exponent to alpha\n");
  fprintf (stderr, "  -O omegaRef    set amplitude to omegaRef\n");
  fprintf (stderr, "  -F fRef        set normalization reference frequency to fRef\n");
  fprintf (stderr, "  -f f0          set start frequency to f0\n");
  fprintf (stderr, "  -e deltaF      set frequency spacing to deltaF\n");
  fprintf (stderr, "  -n length      set number of points in frequency series to length\n");
  fprintf (stderr, "  -o filename    print gravitational-wave spectrum to file filename\n");
  exit (exitcode);
}

/*
 * ParseOptions ()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */
static void
ParseOptions (int argc, char *argv[])
{
  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvd:a:O:F:e:f:n:o:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'o': /* specify output file */
        strncpy (optFile, optarg, LALNameLength);
        break;

      case 'n': /* specify number of points in frequency series */
        optLength = atoi (optarg);
        break;

      case 'e': /* specify frequency resolution */
        optDeltaF = atof (optarg);
        break;

      case 'f': /* specify start frequency */
        optF0 = atof (optarg);
        break;

      case 'O': /* specify amplitude at reference frequency */
        optOmegaR = atof (optarg);
        break;

      case 'F': /* specify reference frequency */
        optFR = atof (optarg);
        break;

      case 'a': /* specify power law exponent */
        optAlpha = atof (optarg);
        break;

      case 'd': /* set debug level */
        lalDebugLevel = atoi (optarg);
        break;

      case 'v': /* optVerbose */
        optVerbose = STOCHASTICOMEGAGWTESTC_TRUE;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        freopen ("/dev/null", "w", stderr);
        freopen ("/dev/null", "w", stdout);
        break;

      case 'h':
        Usage (argv[0], 0);
        break;

      default:
        Usage (argv[0], 1);
    }

  }

  if (optind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}
