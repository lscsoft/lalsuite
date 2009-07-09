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

/********* <lalVerbatim file="StochasticCrossCorrelationSpectrumTestCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu (original by S. Drasco)
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{StochasticCrossCorrelationSpectrumTest.c}}
\label{stochastic:ss:StochasticCrossCorrelationSpectrumTest.c}

A program to test \texttt{LALStochasticCrossCorrelationSpectrum()}.

\subsubsection*{Usage}

\begin{verbatim}
./StochasticCrossCorrelationSpectrumTest [options]
Options:
  -h             print usage message
  -q             quiet: run silently
  -v             verbose: print extra information
  -d level       set lalDebugLevel to level
  -o filename    write spectrum to file filename
  -i filename    read first data stream from file filename
  -j filename    read second data stream from file filename
  -k filename    read optimal filter from file filename
  -m length      optimal filter contains length points
  -n length      data streams contain length points
  -t             epochs need not match
\end{verbatim}

This program tests the function
\texttt{LALStochasticCrossCorrelationSpectrum()}, which calculates
the cross-correlation spectrum given two zero-padded and
Fourier-transformed data streams and a (frequency domain) optimal
filter.

First, it tests that the correct error codes
(\textit{cf.}\ Sec.~\ref{stochastic:s:StochasticCrossCorrelation.h})
are generated for the following error conditions (tests in
\textit{italics} are not performed if \verb+LAL_NDEBUG+ is set, as
the corresponding checks in the code are made using the ASSERT macro):
\begin{itemize}
\item \textit{null pointer to output series}
\item \textit{null pointer to input structure}
\item \textit{null pointer to first data stream}
\item \textit{null pointer to second data stream}
\item \textit{null pointer to optimal filter}
\item \textit{null pointer to data member of first data stream}
\item \textit{null pointer to data member of second data stream}
\item \textit{null pointer to data member of optimal filter}
\item \textit{null pointer to data member of data member of first data stream}
\item \textit{null pointer to data member of data member of second data stream}
\item \textit{null pointer to data member of data member of optimal filter}
\item \textit{zero length}
\item \textit{negative frequency spacing}
\item \textit{zero frequency spacing}
\item negative start frequency
\item length mismatch between data streams
\item frequency spacing mismatch between data streams
\item start frequency mismatch between data streams
\item mismatch between epochs of data streams
\end{itemize}

It then verifies that the correct cross-correlation statistic (value
and units) is generated for the following simple test case:
$\widetilde{Q}(f) = x(1-x)$; $\widetilde{\bar{h}}_1(f)=x^2+ix$,
  $\widetilde{\bar{h}}_2(f)=x^{-2}-ix^{-1}$, with $x=f/400\,\textrm{Hz}$.
  The expected result in this case is $Y(f)=-i(2-x)(x^2+2)$.
For each successful test
(this valid data test and the invalid ones described above), it
prints ``\texttt{PASS}'' to standard output; if a test fails, it
prints ``\texttt{FAIL}''.

If the \texttt{filename} arguments are present, it also reads in the
optimal filter and the two data streams from the specified files and
use the specified parameters to calculate the cross-correlation
statistic.  The result is printed to the specified output file.

\subsubsection*{Exit codes}
\input{StochasticCrossCorrelationSpectrumTestCE}

\subsubsection*{Uses}

\begin{verbatim}
LALStochasticCrossCorrelationSpectrum()
LALCheckMemoryLeaks()
LALCReadFrequencySeries()
LALCCreateVector()
LALCDestroyVector()
LALCHARCreateVector()
LALCHARDestroyVector()
LALUnitAsString()
LALUnitCompare()
getopt()
printf()
fprintf()
freopen()
fabs()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}
  \item No specific error checking is done on user-specified data.  If
    \texttt{length} is missing, the resulting default will cause a bad
    data error.
  \item The length of the user-provided series must be specified, even
    though it could in principle be deduced from the input file,
    because the data sequences must be allocated before the
    \texttt{LALCReadFrequencySeries()} function is called.
  \item If some, but not all, of the \texttt{filename} arguments are
    present, the user-specified data will be silently ignored.
\end{itemize}

\vfill{\footnotesize\input{StochasticCrossCorrelationSpectrumTestCV}}

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
#include <lal/ReadFTSeries.h>
#include <lal/PrintFTSeries.h>
#include <lal/Units.h>

#include "CheckStatus.h"

NRCSID (STOCHASTICCROSSCORRELATIONSPECTRUMTESTC, "$Id$");

#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH    9
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_F0        80.0
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_DELTAF    80.0
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TOL       1e-6

#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_WINMIN   300.0
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_WINMAX   500.0
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_FLIM     800.0
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EXP2     116.8

#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE     1
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_FALSE    0

extern char *optarg;
extern int   optind;

/* int lalDebugLevel = LALMSGLVL3; */
int lalDebugLevel = LALNDEBUG;
BOOLEAN optVerbose    = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_FALSE;
BOOLEAN optMatch    = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE;
UINT4 optStreamLength     = 0;
UINT4 optFilterLength     = 0;
CHAR optData1File[LALNameLength] = "";
CHAR optData2File[LALNameLength] = "";
CHAR optFilterFile[LALNameLength] = "";
CHAR optOutputFile[LALNameLength] = "";

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

/*************** <lalErrTable file="StochasticCrossCorrelationSpectrumTestCE"> */
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ENOM 0
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EARG 1
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK 2
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS 3
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EUSE 4

#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGENOM "Nominal exit"
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEARG "Error parsing command-line arguments"
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK "Error checking failed to catch bad data"
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS "Incorrect answer for valid data"
#define STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEUSE "Bad user-entered data"
/***************************** </lalErrTable> */

int main( int argc, char *argv[] )
{
  static LALStatus                status;

  StochasticCrossCorrelationInput           input;

  COMPLEX8FrequencySeries  goodData1;
  COMPLEX8FrequencySeries  goodData2;
  COMPLEX8FrequencySeries  goodFilter;
  COMPLEX8FrequencySeries  goodOutput;

  COMPLEX8FrequencySeries  badData1;
  COMPLEX8FrequencySeries  badData2;
  COMPLEX8FrequencySeries  badFilter;
  COMPLEX8FrequencySeries  badOutput;

  COMPLEX8                *tempPtr;
  LIGOTimeGPS              epoch0 = {0,0};
  LIGOTimeGPS              epoch1 = {630720000,123456789};
  LIGOTimeGPS              epoch2 = {630720000,987654321};
  LIGOTimeGPS              epoch3 = {630722222,123456789};

  UINT4 i;
  REAL4 f, x;
  INT4  code;
  REAL4 expIm;

  ParseOptions( argc, argv );

  /* define valid parameters */

  goodFilter.f0     = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_F0;
  goodFilter.deltaF = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_DELTAF;
  goodFilter.epoch  = epoch0;
  goodFilter.data   = NULL;

  goodData1 = goodFilter;
  goodOutput = goodFilter;

  goodData1.epoch = epoch1;
  badFilter = goodFilter;
  badData2 = badData1 = goodData2 = goodData1;
  badOutput = goodOutput;

  LALCCreateVector(&status, &(goodData1.data),
                          STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCCreateVector(&status, &(goodData2.data),
                          STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCCreateVector(&status, &(goodFilter.data),
                          STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCCreateVector(&status, &(goodOutput.data),
                          STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
  {
    return code;
  }

  input.hBarTildeOne  = &goodData1;
  input.hBarTildeTwo  = &goodData2;
  input.optimalFilter = &goodFilter;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to output series */
    LALStochasticCrossCorrelationSpectrum(&status, NULL, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
      {
        return code;
      }
    printf("  PASS: null pointer to output series results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to input structure */
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, NULL, STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to input structure results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to first data stream */
    input.hBarTildeOne = NULL;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to first data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to second data stream */
    input.hBarTildeOne = &goodData1;

    /* test behavior for null pointer to second data stream */
    input.hBarTildeTwo = NULL;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to second data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to second data stream */
    input.hBarTildeTwo = &goodData2;

    /* test behavior for null pointer to optimal filter */
    input.optimalFilter = NULL;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to optimal filter results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to optimal filter */
    input.optimalFilter = &goodFilter;

    /* test behavior for null pointer to data member of output series */
    LALStochasticCrossCorrelationSpectrum(&status, &badOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of output series results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to data member of first data stream */
    input.hBarTildeOne = &badData1;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of first data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of first data stream */
    input.hBarTildeOne = &goodData1;

    /* test behavior for null pointer to data member of second data stream */
    input.hBarTildeTwo = &badData2;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of second data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of second data stream */
    input.hBarTildeTwo = &goodData2;

    /* test behavior for null pointer to data member of optimal filter */
    input.optimalFilter = &badFilter;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of optimal filter results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of optimal filter */
    input.optimalFilter = &goodFilter;

    /* Create a vector for testing null data-data pointers */
    LALCCreateVector(&status, &(badFilter.data),
                          STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
    {
      return code;
    }
    tempPtr = badFilter.data->data;
    badFilter.data->data = NULL;
    badOutput.data = badData1.data = badData2.data = badFilter.data;

    /* test behavior for null pointer to data member of data member of output series */
    LALStochasticCrossCorrelationSpectrum(&status, &badOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of output series results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to data member of data member of first data stream */
    input.hBarTildeOne = &badData1;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of first data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of data member of first data stream */
    input.hBarTildeOne = &goodData1;

    /* test behavior for null pointer to data member of data member of second data stream */
    input.hBarTildeTwo = &badData2;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of second data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of data member of second data stream */
    input.hBarTildeTwo = &goodData2;

    /* test behavior for null pointer to data member of data member of optimal filter */
    input.optimalFilter = &badFilter;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of optimal filter results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of data member of optimal filter */
    input.optimalFilter = &goodFilter;

    /* clean up */

    badFilter.data->data = tempPtr;
    LALCDestroyVector(&status, &(badFilter.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
    {
      return code;
    }
    badData1.data = badData2.data = badFilter.data;

    /* test behavior for zero length */
    goodData1.data->length = goodData2.data->length
      = goodFilter.data->length = 0;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EZEROLEN,
			      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero length results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);
    /* reassign valid length */
    goodData1.data->length = goodData2.data->length
      = goodFilter.data->length = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH;

    /* test behavior for negative frequency spacing */
    goodData1.deltaF = goodData2.deltaF
      = goodFilter.deltaF = -STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_DELTAF;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

    /* test behavior for zero frequency spacing */
    goodData1.deltaF = goodData2.deltaF
      = goodFilter.deltaF = 0;
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* reassign valid frequency spacing */
    goodData1.deltaF = goodData2.deltaF
      = goodFilter.deltaF = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_DELTAF;
   } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

  /* test behavior for negative start frequency */
  goodData1.f0 = goodData2.f0
    = goodFilter.f0 = -20.0;
  LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: negative start frequency results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);

  /* reassign valid start frequency */
  goodData1.f0 = goodData2.f0
    = goodFilter.f0 = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_F0;

  /* test behavior for length mismatch between data streams */
  goodData2.data->length = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH - 1;
  LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: length mismatch between data streams results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);

  /* reassign correct length */
  goodData2.data->length = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH;

  /* test behavior for frequency spacing mismatch between data streams */
  goodData2.deltaF = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_DELTAF * 2.0;
  LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: frequency spacing mismatch between data streams results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);

  /* reassign correct frequency spacing */
  goodData2.deltaF = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_DELTAF;

  /* test behavior for start frequency mismatch between data streams */
  goodData2.f0 = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_F0 + 2.0;
  LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between data streams results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);

  /* reassign correct start frequency */
  goodData2.f0 = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_F0;

  /* test behavior for mismatch between epochs of data streams */
  goodData2.epoch = epoch2;
  LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMTIME,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMTIME,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
  {
    return code;
  }
  goodData2.epoch = epoch3;
  LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMTIME,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMTIME,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ECHK,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: mismatch between epochs of data streams results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMTIME);

  /* reassign correct epoch */
  goodData2.epoch = epoch1;

  /******************** Test Valid Data Case #1 ***********************/
  goodData1.sampleUnits = lalDimensionlessUnit;
  goodData1.sampleUnits.unitNumerator[LALUnitIndexStrain] = 1;
  goodData1.sampleUnits.unitNumerator[LALUnitIndexSecond] = 1;
  goodData2.sampleUnits = goodData1.sampleUnits;
  goodFilter.sampleUnits = lalDimensionlessUnit;
  goodFilter.sampleUnits.unitNumerator[LALUnitIndexStrain] = -1;

  goodData1.f0 = goodData2.f0 = goodFilter.f0 = 0.0;

  goodData1.data->data[0].re = goodData1.data->data[0].im
    = goodData2.data->data[0].re = goodData2.data->data[0].im
    = goodFilter.data->data[0].re = goodFilter.data->data[0].im
    = 0.0;

  for (i=1; i<STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH; ++i)
  {
    f = i * STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_DELTAF;
    x = f / (STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_FLIM / 2.0);
    goodData1.data->data[i].re = x*x;
    goodData1.data->data[i].im = x;
    goodData2.data->data[i].re = 1.0/goodData1.data->data[i].re;
    goodData2.data->data[i].im = -1.0/goodData1.data->data[i].im;
    goodFilter.data->data[i].re = x * (2-x);
    goodFilter.data->data[i].im = 0.0;
  }

  LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input,  STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
  {
    return code;
  }

   /* check output values */

  if (optVerbose)
  {
    printf("Y(0)=%g + %g i, should be 0\n",
           goodOutput.data->data[0].re, goodOutput.data->data[0].im);
  }
  if ( ( fabs(goodOutput.data->data[0].re)
         > STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TOL )
       || ( fabs(goodOutput.data->data[0].im)
            > STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TOL ) )
  {
    printf("  FAIL: Valid data test\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS);
    }
    return STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS;
  }

  for (i=1; i<STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_LENGTH; ++i)
  {
    f = i * STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_DELTAF;
    x = f / (STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_FLIM / 2.0);
    expIm = -1 * (2-x) * (1 + x*x);

    if (optVerbose)
    {
      printf("Y(%f Hz)=%g + %g i, should be %g i\n",
             f, goodOutput.data->data[i].re, goodOutput.data->data[i].im,
             expIm);
     }
     if ( fabs(goodOutput.data->data[i].re)
          > STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TOL
          || fabs(goodOutput.data->data[i].im - expIm)
          > STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TOL )
     {
       printf("  FAIL: Valid data test\n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n", STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS);
       }
       return STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS;
     }
   }


  /* clean up */
  LALCDestroyVector(&status, &(goodOutput.data));
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCDestroyVector(&status, &(goodFilter.data));
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCDestroyVector(&status, &(goodData1.data));
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCDestroyVector(&status, &(goodData2.data));
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			    STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
  {
    return code;
  }

  printf("PASS: all tests\n");
  LALCheckMemoryLeaks();

  if (optData1File[0] && optData2File[0] && optFilterFile[0]
      && optOutputFile[0])
  {

    /* Allocate Memory */
    LALCCreateVector(&status, &(goodOutput.data), optFilterLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EUSE,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(goodFilter.data), optFilterLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EUSE,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(goodData1.data), optStreamLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EUSE,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(goodData2.data), optStreamLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EUSE,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEUSE) ) )
    {
      return code;
    }
    /* Read Data From Files */
    LALCReadFrequencySeries(&status, &(goodFilter), optFilterFile);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EUSE,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCReadFrequencySeries(&status, &(goodData1), optData1File);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EUSE,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCReadFrequencySeries(&status, &(goodData2), optData2File);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EUSE,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEUSE) ) )
    {
      return code;
    }
    /* Calculate CC Spectrum */
    LALStochasticCrossCorrelationSpectrum(&status, &goodOutput, &input, optMatch);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EUSE,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEUSE) ) )
    {
      return code;
    }

    /* Print result to file */

    LALCPrintFrequencySeries(&goodOutput, optOutputFile);
    printf("=== Cross-Correlation Spectrum for User-Specified Data Written to File %s ===\n",
         optOutputFile);

    /* Deallocate Memory */
    LALCDestroyVector(&status, &(goodOutput.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(goodFilter.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(goodData1.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(goodData2.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_EFLS,
			      STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_MSGEFLS) ) )
    {
      return code;
    }
  }


  /* normal exit */
  LALCheckMemoryLeaks();
  return STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_ENOM;
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
  fprintf (stderr, "  -o filename    write spectrum to file filename\n");
  fprintf (stderr, "  -i filename    read first data stream from file filename\n");
  fprintf (stderr, "  -j filename    read second data stream from file filename\n");
  fprintf (stderr, "  -k filename    read optimal filter from file filename\n");
  fprintf (stderr, "  -m length      optimal filter contains length points\n");
  fprintf (stderr, "  -n length      data streams contain length points\n");
  fprintf (stderr, "  -t             epochs need not match\n");
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

    c = getopt (argc, argv, "hqvd:i:j:k:m:n:o:t");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 't': /* epochs need not match */
        optMatch = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_FALSE;
        break;

      case 'o': /* specify output file */
        strncpy (optOutputFile, optarg, LALNameLength);
        break;

      case 'i': /* specify file containing first data stream */
        strncpy (optData1File, optarg, LALNameLength);
        break;

      case 'j': /* specify file containing second data stream */
        strncpy (optData2File, optarg, LALNameLength);
        break;

      case 'k': /* specify file containing optimal filter */
        strncpy (optFilterFile, optarg, LALNameLength);
        break;

      case 'm': /* specify number of points in optimal filter */
        optFilterLength = atoi (optarg);
        break;

      case 'n': /* specify number of points in data streams */
        optStreamLength = atoi (optarg);
        break;

      case 'd': /* set debug level */
        lalDebugLevel = atoi (optarg);
        break;

      case 'v': /* optVerbose */
        optVerbose = STOCHASTICCROSSCORRELATIONSPECTRUMTESTC_TRUE;
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
