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

/****************** <lalVerbatim file="CCoarseGrainFrequencySeriesTestCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{CCoarseGrainFrequencySeriesTest.c}}
\label{utilities:ss:CCoarseGrainFrequencySeriesTest.c}

Test suite for \texttt{LALCCoarseGrainFrequencySeries()}.

\subsubsection*{Usage}
\begin{verbatim}
./CCoarseGrainFrequencySeriesTest
Options:
  -h             print usage message
  -q             quiet: run silently
  -v             verbose: print extra information
  -d level       set lalDebugLevel to level
  -i filename    read fine grained series from file filename
  -o filename    print coarse grained  series to file filename
  -n length      input series contains length points
  -m length      output series contains length points
  -e deltaF      set coarse grained frequency spacing to deltaF
  -f f0          set start frequency of output to f0
\end{verbatim}

\subsubsection*{Description}

This program tests the routine
\texttt{LALCCoarseGrainFrequencySeries()}, which coarse-grains a
frequency series.

First, it tests that the correct error codes
(\textit{cf.}\ Sec.~\ref{utilities:s:CoarseGrainFrequencySeries.h})
are generated for the following error conditions (tests in
\textit{italics} are not performed if \verb+LAL_NDEBUG+ is set, as
the corresponding checks in the code are made using the ASSERT macro):
\begin{itemize}
\item \textit{null pointer to output series}
\item \textit{null pointer to input series}
\item \textit{null pointer to data member of output series}
\item \textit{null pointer to data member of input series}
\item \textit{null pointer to data member of data member of input series}
\item \textit{null pointer to data member of data member of output series}
% \item \textit{duplicate pointers to input and output series}
% \item \textit{duplicate pointers to data members of input and output series}
% \item \textit{duplicate pointers to data members of data members of input and output series}
\item \textit{zero length}
\item \textit{negative frequency spacing}
\item \textit{zero frequency spacing}
\end{itemize}

It then verifies that the correct values are obtained for some simple
test cases
\begin{itemize}
\item $\{h_\ell'\}=\{0,1,2,3,4,5,6,7\}$, $f'_0=f_0$, $\delta f'=1$, $\delta
f=2$, $N=3$; the expected output is $\{h_k\}=\{1/2,2,4,6\}$.
\item $\{h_\ell'\}=\{0,1,2,3,4,5,6,7\}$, $f'_0=f_0$, $\delta f'=1$, $\delta
f=3$, $N=3$; the expected output is $\{h_k\}=\{2/3,3,6\}$.
\item $f_0'=40$, $\delta f'= 1$,
  $\{h_k\}=\{f_k+i\,f_k^{-1}|k=0,\ldots,4\}$, $f_0=41$,
$f_0=f'_0$, $\delta f=2$.
$\delta f'=3$, $N=$ ; the expected output is
$$
\{h'_\ell\}=\left\{41+i\left(\frac{1}{40}+\frac{2}{41}+\frac{1}{42}\right),
  43+i\left(\frac{1}{42}+\frac{2}{43}+\frac{1}{44}\right)
\right\}
$$
% '
\end{itemize}
For each successful test (both of these valid data and the invalid
ones described above), it prints ``\texttt{PASS}'' to standard output;
if a test fails, it prints ``\texttt{FAIL}''.

If the \texttt{filename} arguments are present, it also reads a
frequency series from a file, calls
\texttt{LALCCoarseGrainFrequencySeries()}, and writes the results to
the specified output file.

\subsubsection*{Exit codes}
\input{CCoarseGrainFrequencySeriesTestCE}

\subsubsection*{Uses}
\begin{verbatim}
LALCCoarseGrainFrequencySeries()
LALCheckMemoryLeaks()
LALCReadFrequencySeries()
LALCPrintFrequencySeries()
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
\item In addition to the error checks tested in this routine, the
  function checks for errors related to inconsistency of coarse
  graining parameters, as well as duplicate input and output pointers.
  Tests of these error checks are still to be
  added to this test program.
\item No specific error checking is done on user-specified data.  If
  \texttt{length} is missing, the resulting default will cause a bad
  data error.
\item The length of the user-provided series must be specified, even
  though it could in principle be deduced from the input file, because
  the data sequences must be allocated before the
  \texttt{LALCReadFrequencySeries()} function is called.
\item If one \texttt{filename} argument, but not both, is present,
  the user-specified data will be silently ignored.
\end{itemize}

\vfill{\footnotesize\input{CCoarseGrainFrequencySeriesTestCV}}

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

#include <lal/CoarseGrainFrequencySeries.h>
#include <lal/AVFactories.h>
#include <lal/ReadFTSeries.h>
#include <lal/PrintFTSeries.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/TimeFreqFFT.h>

#include "CheckStatus.h"

NRCSID(CCOARSEGRAINFREQUENCYSERIESTESTC, "$Id$");

#define CCOARSEGRAINFREQUENCYSERIESTESTC_TOL           1e-6

#define CCOARSEGRAINFREQUENCYSERIESTESTC_RANTOL           1e-5

#define CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC      1234
#define CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS       56789

#define CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0    1.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_F00        0.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0    8

#define CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1    2.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_F01        0.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1    4

#define CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2    3.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_F02        0.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2    3

#define CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF3    1.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_F03        40.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3    5

#define CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4    2.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_F04        41.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4    2

#define CCOARSEGRAINFREQUENCYSERIESTESTC_TSDT     (1.0 / 1024.0)
#define CCOARSEGRAINFREQUENCYSERIESTESTC_TSLEN     16384

#define CCOARSEGRAINFREQUENCYSERIESTESTC_FMIN     50
#define CCOARSEGRAINFREQUENCYSERIESTESTC_FMAX     500

#define CCOARSEGRAINFREQUENCYSERIESTESTC_RESRATIO 9.14

#define CCOARSEGRAINFREQUENCYSERIESTESTC_TRUE     1
#define CCOARSEGRAINFREQUENCYSERIESTESTC_FALSE    0

extern char *optarg;
extern int   optind;

/* int lalDebugLevel = LALMSGLVL3; */
int lalDebugLevel  = LALNDEBUG;
BOOLEAN optVerbose = CCOARSEGRAINFREQUENCYSERIESTESTC_FALSE;
UINT4 optInLength    = 0;
UINT4 optOutLength   = 0;
REAL8 optDeltaF     = -1.0;
REAL8 optF0       = 0.0;

INT4 optSeed = 1682;

CHAR optInputFile[LALNameLength] = "";
CHAR optOutputFile[LALNameLength] = "";
INT4 code;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

/***************************** <lalErrTable file="CCoarseGrainFrequencySeriesTestCE"> */
#define CCOARSEGRAINFREQUENCYSERIESTESTC_ENOM 0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_EARG 1
#define CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK 2
#define CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS 3
#define CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE 4

#define CCOARSEGRAINFREQUENCYSERIESTESTC_MSGENOM "Nominal exit"
#define CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEARG "Error parsing command-line arguments"
#define CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK "Error checking failed to catch bad data"
#define CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS "Incorrect answer for valid data"
#define CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE "Bad user-entered data"
/***************************** </lalErrTable> */

int
main( int argc, char *argv[] )
{

   static LALStatus         status;

   UINT4      i;
   REAL8      f;

   COMPLEX8                   *cPtr;

   const COMPLEX8  testInputDataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0]
     = {{0.0,0.0}, {1.0,0.0}, {2.0,0.0}, {3.0,0.0},
        {4.0,0.0}, {5.0,0.0}, {6.0,0.0}, {7.0,0.0}};

   const COMPLEX8
     expectedOutput1DataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1]
     = {{0.5,0.0}, {2.0,0.0}, {4.0,0.0}, {6.0,0.0}};

   const COMPLEX8
     expectedOutput2DataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2]
     = {{(2.0/3.0),0.0}, {3.0,0.0}, {6.0,0.0}};

   const COMPLEX8
     testInput3DataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3]
     = {{40.0,1.0/40.0}, {41.0,1.0/41.0}, {42.0,1.0/42.0},
        {43.0,1.0/43.0}, {44.0,1.0/44.0}};

   const COMPLEX8
     expectedOutput4DataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4]
     = {{41.0, (1.0/40.0+2.0/41.0+1.0/42.0) / 4.0},
        {43.0, (1.0/42.0+2.0/43.0+1.0/44.0) / 4.0}};

   COMPLEX8FrequencySeries             goodInput, badInput;
   COMPLEX8FrequencySeries     goodOutput, badOutput;

   BOOLEAN                result;
   LALUnitPair            unitPair;

   CHARVector             *unitString;

   FrequencySamplingParams     params;

   COMPLEX8               coarseTotal, fineTotal;
   COMPLEX8               cError;
   const COMPLEX8         cZero = {0,0};
   REAL4TimeSeries        timeSeries;

   RealFFTPlan            *fftPlan;

   RandomParams           *randomParams;

   ParseOptions( argc, argv );

   /* TEST INVALID DATA HERE ------------------------------------------- */

   /* define valid parameters */
   goodInput.f0                   = CCOARSEGRAINFREQUENCYSERIESTESTC_F00;
   goodInput.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   goodInput.epoch.gpsSeconds     = CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC;
   goodInput.epoch.gpsNanoSeconds = CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS;
   goodInput.data                 = NULL;
   goodOutput.data                = NULL;

   params.f0                      = CCOARSEGRAINFREQUENCYSERIESTESTC_F00;
   params.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   params.length               = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;

   badInput = goodInput;
   badOutput = goodOutput;

   /* allocate input and output vectors */
   LALCCreateVector(&status, &(goodInput.data),
                    CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }
   LALCCreateVector(&status, &(goodOutput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

#ifndef LAL_NDEBUG
   if ( ! lalNoDebug )
   {
     /* test behavior for null pointer to output series */
     LALCCoarseGrainFrequencySeries(&status, NULL, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to input series */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, NULL, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to parameter structure */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, NULL);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to parameter structure results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of output series */
     LALCCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of input series */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of data member of output series */
     LALCCreateVector(&status, &(badOutput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     cPtr = badOutput.data->data;
     badOutput.data->data = NULL;
     LALCCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = cPtr;
     LALCDestroyVector(&status, &(badOutput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for null pointer to data member of data member of output series */
     LALCCreateVector(&status, &(badInput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     cPtr = badInput.data->data;
     badInput.data->data = NULL;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badInput.data->data = cPtr;
     LALCDestroyVector(&status, &(badInput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for duplicate pointers */

     /* input and output series */
     LALCCoarseGrainFrequencySeries(&status, &goodInput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     badOutput = goodInput;
     badOutput.data = goodInput.data;

     /* data members of input and output series */
     LALCCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to data members of input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     badOutput.data = NULL;

     /* data members of data members of input and output series */
     LALCCreateVector(&status, &(badOutput.data),
		      CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     cPtr = badOutput.data->data;
     badOutput.data->data = goodInput.data->data;
     LALCCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to data members of data members of input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = cPtr;
     LALCDestroyVector(&status, &(badOutput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for zero length */
     /* input */

     goodInput.data->length = 0;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
			       COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero length in input results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

     goodInput.data->length = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;
     goodOutput.data->length = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

     /* output */

     goodOutput.data->length = params.length = 0;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
			       COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero length in output results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

     goodOutput.data->length = params.length
       = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;

     /* test behavior for negative frequency spacing */
     goodInput.deltaF = params.deltaF
       = -CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status,
			       COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
			       COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

     /* test behavior for zero frequency spacing */
     goodInput.deltaF = params.deltaF = 0;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status,
			       COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
			       COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

     /* reassign valid frequency spacing */
     goodInput.deltaF = params.deltaF
       = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;

   } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

   LALCDestroyVector(&status, &(goodOutput.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* TEST VALID DATA HERE --------------------------------------------- */

   params.f0                      = CCOARSEGRAINFREQUENCYSERIESTESTC_F01;
   params.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
   params.length               = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

   /* allocate input and output vectors */
   LALCCreateVector(&status, &(goodOutput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* fill input time-series parameters */
   strncpy(goodInput.name,"Dummy test data",LALNameLength);
   goodInput.sampleUnits  = lalDimensionlessUnit;

     /* fill input data */
   for (i=0; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0; ++i)
   {
     goodInput.data->data[i] = testInputDataData[i];
   }

   /* coarse grain */
   LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            CCOARSEGRAINFREQUENCYSERIESTESTC_F01);
   }
   if (goodOutput.f0)
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1);
   }
   if ( fabs(goodOutput.deltaF-CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1)
        / CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1
        > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   unitPair.unitOne = &(goodInput.sampleUnits);
   unitPair.unitTwo = &(goodOutput.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
   }

   if (!result)
   {
     printf("  FAIL: Valid data test #1\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%g + %g i, should be %g + %g i\n",
            goodOutput.data->data[0].re, goodOutput.data->data[0].im,
            expectedOutput1DataData[0].re, expectedOutput1DataData[0].im);
   }
   if ((fabs(goodOutput.data->data[0].re - expectedOutput1DataData[0].re)
        /expectedOutput1DataData[0].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
        ||
       (fabs(goodOutput.data->data[0].im - expectedOutput1DataData[0].im)
        /expectedOutput1DataData[0].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
   {
     printf("  FAIL: Valid data test #1\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1; ++i)
   {
     f = i * CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%g + %g i, should be %g + %g i\n", f,
              goodOutput.data->data[i].re, goodOutput.data->data[i].im,
              expectedOutput1DataData[i].re, expectedOutput1DataData[i].im);
     }
     if ((fabs(goodOutput.data->data[i].re - expectedOutput1DataData[i].re)
          /expectedOutput1DataData[i].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         ||
         (fabs(goodOutput.data->data[i].im - expectedOutput1DataData[i].im)
          /expectedOutput1DataData[i].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         )
     {
       printf("  FAIL: Valid data test #1\n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   LALCDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }


   /*-------Test #2-------*/

   params.f0                      = CCOARSEGRAINFREQUENCYSERIESTESTC_F02;
   params.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
   params.length               = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2;

   LALCCreateVector(&status, &(goodOutput.data),
                    CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* coarse grain */
   LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            CCOARSEGRAINFREQUENCYSERIESTESTC_F02);
   }
   if (goodOutput.f0)
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2);
   }
   if ( fabs(goodOutput.deltaF-CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2)
        / CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2 > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   unitPair.unitOne = &(goodInput.sampleUnits);
   unitPair.unitTwo = &(goodOutput.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
   }

   if (!result)
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%g + %g i, should be %g + %g i\n",
            goodOutput.data->data[0].re, goodOutput.data->data[0].im,
            expectedOutput2DataData[0].re, expectedOutput2DataData[0].im);
   }
   if ((fabs(goodOutput.data->data[0].re - expectedOutput2DataData[0].re)
        /expectedOutput2DataData[0].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
        ||
       (fabs(goodOutput.data->data[0].im - expectedOutput2DataData[0].im)
        /expectedOutput2DataData[0].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2; ++i)
   {
     f = i * CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%g + %g i, should be %g + %g i\n", f,
              goodOutput.data->data[i].re, goodOutput.data->data[i].im,
              expectedOutput2DataData[i].re, expectedOutput2DataData[i].im);
     }
     if ((fabs(goodOutput.data->data[i].re - expectedOutput2DataData[i].re)
          /expectedOutput2DataData[i].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         ||
         (fabs(goodOutput.data->data[i].im - expectedOutput2DataData[i].im)
          /expectedOutput2DataData[i].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
     {
       printf("  FAIL: Valid data test #2 \n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   /* clean up valid data */
   LALCDestroyVector(&status, &goodInput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCheckMemoryLeaks();

   /*-------Test #3-------*/

   goodInput.f0               = CCOARSEGRAINFREQUENCYSERIESTESTC_F03;
   goodInput.deltaF           = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF3;

   LALCCreateVector(&status, &(goodInput.data),
                    CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* fill input data */
   for (i=0; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3; ++i)
   {
     goodInput.data->data[i] = testInput3DataData[i];
   }

   params.f0                  = CCOARSEGRAINFREQUENCYSERIESTESTC_F04;
   params.deltaF              = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4;
   params.length              = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4;

   LALCCreateVector(&status, &(goodOutput.data),
                    CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* coarse grain */
   LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            CCOARSEGRAINFREQUENCYSERIESTESTC_F04);
   }
   if ( fabs(goodOutput.f0-CCOARSEGRAINFREQUENCYSERIESTESTC_F04)
        / CCOARSEGRAINFREQUENCYSERIESTESTC_F04 > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4);
   }
   if ( fabs(goodOutput.deltaF-CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4)
        / CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4 > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   unitPair.unitOne = &(goodInput.sampleUnits);
   unitPair.unitTwo = &(goodOutput.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
   }

   if (!result)
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%g + %g i, should be %g + %g i\n",
            goodOutput.data->data[0].re, goodOutput.data->data[0].im,
            expectedOutput4DataData[0].re, expectedOutput4DataData[0].im);
   }
   if ((fabs(goodOutput.data->data[0].re - expectedOutput4DataData[0].re)
        /expectedOutput4DataData[0].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
        ||
       (fabs(goodOutput.data->data[0].im - expectedOutput4DataData[0].im)
        /expectedOutput4DataData[0].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4; ++i)
   {
     f = i * CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%g + %g i, should be %g + %g i\n", f,
              goodOutput.data->data[i].re, goodOutput.data->data[i].im,
              expectedOutput4DataData[i].re, expectedOutput4DataData[i].im);
     }
     if ((fabs(goodOutput.data->data[i].re - expectedOutput4DataData[i].re)
          /expectedOutput4DataData[i].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         ||
         (fabs(goodOutput.data->data[i].im - expectedOutput4DataData[i].im)
          /expectedOutput4DataData[i].re > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
     {
       printf("  FAIL: Valid data test #3 \n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   /* clean up valid data */
   LALCDestroyVector(&status, &goodInput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCheckMemoryLeaks();

   /*-------Test #4-------*/

   strncpy(timeSeries.name,"Random test data",LALNameLength);
   timeSeries.deltaT = CCOARSEGRAINFREQUENCYSERIESTESTC_TSDT;
   timeSeries.f0 = 0.0;
   timeSeries.sampleUnits = lalDimensionlessUnit;
   timeSeries.epoch = goodOutput.epoch;

   timeSeries.data = NULL;

   LALSCreateVector(&status, &(timeSeries.data),
		    CCOARSEGRAINFREQUENCYSERIESTESTC_TSLEN );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   goodInput.data = NULL;

   LALCCreateVector(&status, &(goodInput.data),
		    ( CCOARSEGRAINFREQUENCYSERIESTESTC_TSLEN / 2 + 1 ) );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* construct plan */
   fftPlan = NULL;

   LALCreateForwardRealFFTPlan(&status, &fftPlan,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_TSLEN,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_FALSE);
   if ( ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS ) ) )
   {
     return code;
   }

   randomParams = NULL;

   LALCreateRandomParams( &status, &randomParams, optSeed );

   /* fill time series with normal deviates */

   LALNormalDeviates( &status, timeSeries.data, randomParams );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALTimeFreqRealFFT( &status, &goodInput, &timeSeries, fftPlan );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   for ( i = 0 ;
	 i * goodInput.deltaF <= CCOARSEGRAINFREQUENCYSERIESTESTC_FMIN ;
	 ++i )
   {
     goodInput.data->data[i] = cZero;
   }

   for ( i = CCOARSEGRAINFREQUENCYSERIESTESTC_FMAX / goodInput.deltaF ;
	 i < goodInput.data->length ;
	 ++i )
   {
     goodInput.data->data[i] = cZero;
   }

   fineTotal = cZero;
   for ( i = 0 ; i < goodInput.data->length ; ++i )
   {
     fineTotal.re += goodInput.data->data[i].re;
     fineTotal.im += goodInput.data->data[i].im;
   }
   fineTotal.re *= goodInput.deltaF;
   fineTotal.im *= goodInput.deltaF;

   params.deltaF
     = goodInput.deltaF * CCOARSEGRAINFREQUENCYSERIESTESTC_RESRATIO;
   params.f0  = CCOARSEGRAINFREQUENCYSERIESTESTC_FMIN - params.deltaF;
   params.length = (CCOARSEGRAINFREQUENCYSERIESTESTC_FMAX - params.f0)
                   / params.deltaF + 3;

   goodOutput.data = NULL;
   LALCCreateVector(&status, &(goodOutput.data),
		    params.length );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* coarse grain */
   LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   coarseTotal = cZero;
   for ( i = 0 ; i < goodOutput.data->length ; ++i )
   {
     coarseTotal.re += goodOutput.data->data[i].re;
     coarseTotal.im += goodOutput.data->data[i].im;
   }
   coarseTotal.re *= goodOutput.deltaF;
   coarseTotal.im *= goodOutput.deltaF;

   cError.re = ( (coarseTotal.re - fineTotal.re) * fineTotal.re
		 + (coarseTotal.im - fineTotal.im) * fineTotal.im )
               / ( fineTotal.re * fineTotal.re + fineTotal.im * fineTotal.im );
   cError.im = ( (coarseTotal.im - fineTotal.im) * fineTotal.re
		 - (coarseTotal.re - fineTotal.re) * fineTotal.im )
               / ( fineTotal.re * fineTotal.re + fineTotal.im * fineTotal.im );

   if (optVerbose)
   {
     printf("Integral of fine-grained frequency series = %g + %g i\n",
	    fineTotal.re, fineTotal.im);
     printf("Integral of coarse-grained frequency series = %g + %g i\n",
	    coarseTotal.re, coarseTotal.im);
     printf( "Fractional error is %e + % e i\n",
	     cError.re, cError.im );

     LALSPrintTimeSeries( &timeSeries, "ccg_TimeSeries.dat" );
     LALCPrintFrequencySeries( &goodInput, "ccg_fgFreqSeries.dat" );
     LALCPrintFrequencySeries( &goodOutput, "ccg_cgFreqSeries.dat" );
   }

   if ( cError.re * cError.re + cError.im * cError.im
	> ( CCOARSEGRAINFREQUENCYSERIESTESTC_RANTOL
	    * CCOARSEGRAINFREQUENCYSERIESTESTC_RANTOL )
      )
   {
     printf("  FAIL: Valid data test #4 \n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
	      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   LALDestroyRandomParams( &status, &randomParams );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALDestroyRealFFTPlan(&status, &fftPlan);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALSDestroyVector(&status, &(timeSeries.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCDestroyVector(&status, &(goodInput.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCDestroyVector(&status, &(goodOutput.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCheckMemoryLeaks();

   printf("PASS: all tests\n");

   /**************** Process User-Entered Data, If Any **************/

   if (optInputFile[0] && optOutputFile[0]) {

     params.f0 = optF0;
     params.length = optOutLength;
     params.deltaF = optDeltaF;

     goodInput.data  = NULL;
     goodOutput.data = NULL;

     LALCCreateVector(&status, &goodInput.data, optInLength);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALCCreateVector(&status, &goodOutput.data, optOutLength);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     /* Read input file */
     LALCReadFrequencySeries(&status, &goodInput, optInputFile);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     /* coarse grain */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     LALCPrintFrequencySeries(&goodOutput, optOutputFile);

     printf("===== Coarse-Graining of User-Specified Series Written to File %s =====\n", optOutputFile);

     /* clean up valid data */
     LALCDestroyVector(&status, &goodInput.data);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALCDestroyVector(&status, &goodOutput.data);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALCheckMemoryLeaks();
   }
   return CCOARSEGRAINFREQUENCYSERIESTESTC_ENOM;
}

/*------------------------------------------------------------------------*/

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
  fprintf (stderr, "  -s ranseed        use random number seed ranseed\n");
  fprintf (stderr, "  -i filename    read fine grained series from file filename\n");
  fprintf (stderr, "  -o filename    print coarse grained series to file filename\n");
  fprintf (stderr, "  -n length      input series contains length points\n");
  fprintf (stderr, "  -m length      output series contains length points\n");
  fprintf (stderr, "  -e deltaF      set coarse grained frequency spacing to deltaF\n");
  fprintf (stderr, "  -f f0          set start frequency of output to f0\n");
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

    c = getopt (argc, argv, "hqvd:s:i:o:n:m:e:f:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'i': /* specify input file */
        strncpy (optInputFile, optarg, LALNameLength);
        break;

      case 'o': /* specify output file */
        strncpy (optOutputFile, optarg, LALNameLength);
        break;

      case 'n': /* specify number of points in input series */
        optInLength = atoi (optarg);
        break;

      case 'm': /* specify number of points in output series */
        optOutLength = atoi (optarg);
        break;

      case 'e': /* specify frequency resolution */
        optDeltaF = atof (optarg);
        break;

      case 'f': /* specify start frequency */
        optF0 = atof (optarg);
        break;

      case 's': /* set random number seed */
        optSeed = atoi (optarg);
        break;

      case 'd': /* set debug level */
        lalDebugLevel = atoi (optarg);
        break;

      case 'v': /* optVerbose */
        optVerbose = CCOARSEGRAINFREQUENCYSERIESTESTC_TRUE;
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
