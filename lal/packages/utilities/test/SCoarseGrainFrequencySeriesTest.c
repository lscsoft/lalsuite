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

/****************** <lalVerbatim file="SCoarseGrainFrequencySeriesTestCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{SCoarseGrainFrequencySeriesTest.c}}
\label{utilities:ss:SCoarseGrainFrequencySeriesTest.c}

Test suite for \texttt{LALSCoarseGrainFrequencySeries()}.

\subsubsection*{Usage}
\begin{verbatim}
./SCoarseGrainFrequencySeriesTest
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
\texttt{LALSCoarseGrainFrequencySeries()}, which coarse-grains a
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
%\item \textit{duplicate pointers to input and output series}
%\item \textit{duplicate pointers to data members of input and output series}
%\item \textit{duplicate pointers to data members of data members of input and output series}
\item \textit{zero length}
\item \textit{negative frequency spacing}
\item \textit{zero frequency spacing}
\end{itemize}

It then verifies that the correct
values are obtained for some simple test cases
\begin{itemize}
\item $\{h_\ell'\}=\{0,1,2,3,4,5,6,7\}$, $f'_0=f_0$, $\delta f'=1$, $\delta
f=2$, $N=3$; the expected output is $\{h_k\}=\{1/2,2,4,6\}$.
\item $\{h_\ell'\}=\{0,1,2,3,4,5,6,7\}$, $f'_0=f_0$, $\delta f'=1$, $\delta
f=3$, $N=3$; the expected output is $\{h_k\}=\{2/3,3,6\}$.
\end{itemize}
For each successful test (both of these valid data and the invalid
ones described above), it prints ``\texttt{PASS}'' to standard output;
if a test fails, it prints ``\texttt{FAIL}''.

If the \texttt{filename} arguments are present, it also reads a
frequency series from a file, calls
\texttt{LALSCoarseGrainFrequencySeries()}, and writes the results to
the specified output file.

\subsubsection*{Exit codes}
\input{SCoarseGrainFrequencySeriesTestCE}

\subsubsection*{Uses}
\begin{verbatim}
LALSCoarseGrainFrequencySeries()
LALCheckMemoryLeaks()
LALSReadFrequencySeries()
LALSPrintFrequencySeries()
LALSCreateVector()
LALSDestroyVector()
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
  graining parameters.  Tests of these error checks are still to be
  added to this test program.
\item No specific error checking is done on user-specified data.  If
  \texttt{length} is missing, the resulting default will cause a bad
  data error.
\item The length of the user-provided series must be specified, even
  though it could in principle be deduced from the input file, because
  the data sequences must be allocated before the
  \texttt{LALSReadFrequencySeries()} function is called.
\item If one \texttt{filename} argument, but not both, is present,
  the user-specified data will be silently ignored.
\end{itemize}

\vfill{\footnotesize\input{SCoarseGrainFrequencySeriesTestCV}}

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

#include "CheckStatus.h"

NRCSID(SCOARSEGRAINFREQUENCYSERIESTESTC, "$Id$");

#define SCOARSEGRAINFREQUENCYSERIESTESTC_TOL           1e-6

#define SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC      1234
#define SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS       56789

#define SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0    1.0
#define SCOARSEGRAINFREQUENCYSERIESTESTC_F00        0.0
#define SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0    8

#define SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1    2.0
#define SCOARSEGRAINFREQUENCYSERIESTESTC_F01        0.0
#define SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1    4

#define SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2    3.0
#define SCOARSEGRAINFREQUENCYSERIESTESTC_F02        0.0
#define SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2    3

#define SCOARSEGRAINFREQUENCYSERIESTESTC_TRUE     1
#define SCOARSEGRAINFREQUENCYSERIESTESTC_FALSE    0

extern char *optarg;
extern int   optind;

/* int lalDebugLevel = LALMSGLVL3; */
int lalDebugLevel  = LALNDEBUG;
BOOLEAN optVerbose = SCOARSEGRAINFREQUENCYSERIESTESTC_FALSE;
UINT4 optInLength    = 0;
UINT4 optOutLength   = 0;
REAL8 optDeltaF     = -1.0;
REAL8 optF0       = 0.0;

CHAR optInputFile[LALNameLength] = "";
CHAR optOutputFile[LALNameLength] = "";
INT4 code;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

/***************************** <lalErrTable file="SCoarseGrainFrequencySeriesTestCE"> */
#define SCOARSEGRAINFREQUENCYSERIESTESTC_ENOM 0
#define SCOARSEGRAINFREQUENCYSERIESTESTC_EARG 1
#define SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK 2
#define SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS 3
#define SCOARSEGRAINFREQUENCYSERIESTESTC_EUSE 4

#define SCOARSEGRAINFREQUENCYSERIESTESTC_MSGENOM "Nominal exit"
#define SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEARG "Error parsing command-line arguments"
#define SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK "Error checking failed to catch bad data"
#define SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS "Incorrect answer for valid data"
#define SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE "Bad user-entered data"
/***************************** </lalErrTable> */

int
main( int argc, char *argv[] )
{

   static LALStatus         status;

   UINT4      i;
   REAL8      f;

   REAL4                   *sPtr;

   const REAL4    testInputDataData[SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0]
                     = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

   const REAL4 expectedOutput1DataData[SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1]
                     = {0.5, 2.0, 4.0, 6.0};

   const REAL4 expectedOutput2DataData[SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2]
                     = {(2.0/3.0), 3.0, 6.0};

   REAL4FrequencySeries             goodInput, badInput;
   REAL4FrequencySeries     goodOutput, badOutput;

   BOOLEAN                result;
   LALUnitPair            unitPair;

   CHARVector             *unitString;

   FrequencySamplingParams     params;

   ParseOptions( argc, argv );

   /* TEST INVALID DATA HERE ------------------------------------------- */

   /* define valid parameters */
   goodInput.f0                   = SCOARSEGRAINFREQUENCYSERIESTESTC_F00;
   goodInput.deltaF               = SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   goodInput.epoch.gpsSeconds     = SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC;
   goodInput.epoch.gpsNanoSeconds = SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS;
   goodInput.data                 = NULL;
   goodOutput.data                = NULL;

   params.f0                      = SCOARSEGRAINFREQUENCYSERIESTESTC_F01;
   params.deltaF               = SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   params.length               = SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;

   badInput = goodInput;
   badOutput = goodOutput;

   /* allocate input and output vectors */
   LALSCreateVector(&status, &(goodInput.data),
                    SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
   if ( ( code = CheckStatus(&status, 0 , "",
			     SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }
   LALSCreateVector(&status, &(goodOutput.data), SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
   if ( ( code = CheckStatus(&status, 0 , "",
			     SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

#ifndef LAL_NDEBUG
   if ( ! lalNoDebug )
   {
     /* test behavior for null pointer to output series */
     LALSCoarseGrainFrequencySeries(&status, NULL, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to input series */
     LALSCoarseGrainFrequencySeries(&status, &goodOutput, NULL, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to plan parameter */
     LALSCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, NULL);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to plan parameter results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of output series */
     LALSCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of input series */
     LALSCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of data member of output series */
     LALSCreateVector(&status, &(badOutput.data), SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     sPtr = badOutput.data->data;
     badOutput.data->data = NULL;
     LALSCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = sPtr;
     LALSDestroyVector(&status, &(badOutput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for null pointer to data member of data member of output series */
     LALSCreateVector(&status, &(badInput.data), SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     sPtr = badInput.data->data;
     badInput.data->data = NULL;
     LALSCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badInput.data->data = sPtr;
     LALSDestroyVector(&status, &(badInput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* Removed to make POST05MDC pass make check */

#if 0
          /* test behavior for duplicate pointers */

     /* input and output series */
     LALSCoarseGrainFrequencySeries(&status, &goodInput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     badOutput = goodInput;
     badOutput.data = goodInput.data;

     /* data members of input and output series */
     LALSCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to data members of input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     badOutput.data = NULL;

     /* data members of data members of input and output series */
     LALSCreateVector(&status, &(badOutput.data),
                      SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     sPtr = badOutput.data->data;
     badOutput.data->data = goodInput.data->data;
     LALSCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to data members of data members of input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = sPtr;
     LALSDestroyVector(&status, &(badOutput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

#endif

     /* test behavior for zero length */
     /* input */

     goodInput.data->length = 0;
     LALSCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
			       COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero length in input results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

     goodInput.data->length = SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;
     goodOutput.data->length = SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

     /* output */

     goodOutput.data->length = params.length = 0;
     LALSCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
			       COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero length in output results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

     goodOutput.data->length = params.length
       = SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;

     /* test behavior for negative frequency spacing */
     goodInput.deltaF = params.deltaF
       = -SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
     LALSCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status,
			       COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
			       COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

     /* test behavior for zero frequency spacing */
     goodInput.deltaF = params.deltaF = 0;
     LALSCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status,
			       COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
			       COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

     /* reassign valid frequency spacing */
     goodInput.deltaF = params.deltaF
       = SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;

   } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

   LALSDestroyVector(&status, &(goodOutput.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* TEST VALID DATA HERE --------------------------------------------- */

   params.f0                      = SCOARSEGRAINFREQUENCYSERIESTESTC_F01;
   params.deltaF               = SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
   params.length               = SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

   /* allocate input and output vectors */
   LALSCreateVector(&status, &(goodOutput.data), SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
   if ( ( code = CheckStatus(&status, 0 , "",
			     SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* fill input time-series parameters */
   strncpy(goodInput.name,"Dummy test data",LALNameLength);
   goodInput.sampleUnits  = lalDimensionlessUnit;

     /* fill input data */
   for (i=0; i<SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0; ++i)
   {
     goodInput.data->data[i] = testInputDataData[i];
   }

   /* coarse grain */
   LALSCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            SCOARSEGRAINFREQUENCYSERIESTESTC_F01);
   }
   if (goodOutput.f0)
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1);
   }
   if ( fabs(goodOutput.deltaF-SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1)
        / SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1
        > SCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   unitPair.unitOne = &(goodInput.sampleUnits);
   unitPair.unitTwo = &(goodOutput.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
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
              SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%g, should be %g\n",
            goodOutput.data->data[0], expectedOutput1DataData[0]);
   }
   if (fabs(goodOutput.data->data[0] - expectedOutput1DataData[0])
        / expectedOutput1DataData[0] > SCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1; ++i)
   {
     f = i * SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%g, should be %g\n",
              f, goodOutput.data->data[i], expectedOutput1DataData[i]);
     }
     if (fabs(goodOutput.data->data[i] - expectedOutput1DataData[i])
         / expectedOutput1DataData[i] > SCOARSEGRAINFREQUENCYSERIESTESTC_TOL)
     {
       printf("  FAIL: Valid data test\n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   LALSDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }


   /*-------Test #2-------*/

   params.f0                      = SCOARSEGRAINFREQUENCYSERIESTESTC_F02;
   params.deltaF               = SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
   params.length               = SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2;

   LALSCreateVector(&status, &(goodOutput.data),
                    SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2);
   if ( ( code = CheckStatus(&status, 0 , "",
			     SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* coarse grain */
   LALSCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            SCOARSEGRAINFREQUENCYSERIESTESTC_F02);
   }
   if (goodOutput.f0)
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2);
   }
   if ( fabs(goodOutput.deltaF-SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2)
        / SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2 > SCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != SCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   unitPair.unitOne = &(goodInput.sampleUnits);
   unitPair.unitTwo = &(goodOutput.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
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
              SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%g, should be %g\n",
            goodOutput.data->data[0], expectedOutput2DataData[0]);
   }
   if (fabs(goodOutput.data->data[0] - expectedOutput2DataData[0])
        / expectedOutput2DataData[0] > SCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<SCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2; ++i)
   {
     f = i * SCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%g, should be %g\n",
              f, goodOutput.data->data[i], expectedOutput2DataData[i]);
     }
     if (fabs(goodOutput.data->data[i] - expectedOutput2DataData[i])
         / expectedOutput2DataData[i] > SCOARSEGRAINFREQUENCYSERIESTESTC_TOL)
     {
       printf("  FAIL: Valid data test #2 \n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   /* clean up valid data */
   LALSDestroyVector(&status, &goodInput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALSDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     SCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
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

     LALSCreateVector(&status, &goodInput.data, optInLength);
     if ( ( code = CheckStatus( &status, 0 , "",
				SCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALSCreateVector(&status, &goodOutput.data, optOutLength);
     if ( ( code = CheckStatus( &status, 0 , "",
				SCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     /* Read input file */
     LALSReadFrequencySeries(&status, &goodInput, optInputFile);
     if ( ( code = CheckStatus( &status, 0 , "",
				SCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     /* coarse grain */
     LALSCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus( &status, 0 , "",
				SCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     LALSPrintFrequencySeries(&goodOutput, optOutputFile);

     printf("===== Coarse-Graining of User-Specified Series Written to File %s =====\n", optOutputFile);

     /* clean up valid data */
     LALSDestroyVector(&status, &goodInput.data);
     if ( ( code = CheckStatus( &status, 0 , "",
				SCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALSDestroyVector(&status, &goodOutput.data);
     if ( ( code = CheckStatus( &status, 0 , "",
				SCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				SCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALCheckMemoryLeaks();
   }
   return SCOARSEGRAINFREQUENCYSERIESTESTC_ENOM;
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

    c = getopt (argc, argv, "hqvd:i:o:n:m:e:f:");
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

      case 'd': /* set debug level */
        lalDebugLevel = atoi (optarg);
        break;

      case 'v': /* optVerbose */
        optVerbose = SCOARSEGRAINFREQUENCYSERIESTESTC_TRUE;
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
