/****************** <lalVerbatim file="CCoarseGrainFrequencySeriesTestCV">
Author: UTB Relativity Group; contact J. T. Whelan
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{CCoarseGrainFrequencySeriesTest.c}}
\label{stochastic:ss:CCoarseGrainFrequencySeriesTest.c}

{\bf {\Large WARNING} The functionality of this module has been expanded
and modified, so the tests and  documentation are not yet complete
and/or correct.}

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
(\textit{cf.}\ Sec.~\ref{stochastic:s:CoarseGrainFrequencySeries.h})
are generated for the following error conditions (tests in
\textit{italics} are not performed if \verb+LAL_NEDEBUG+ is set, as
the corresponding checks in the code are made using the ASSERT macro):
\begin{itemize}
\item \textit{null pointer to output series}
\item \textit{null pointer to input series}
\item \textit{null pointer to data member of output series}
\item \textit{null pointer to data member of input series}
\item \textit{null pointer to data member of data member of input series}
\item \textit{null pointer to data member of data member of output series}
\item \textit{zero length}
\item \textit{negative frequency spacing}
\item \textit{zero frequency spacing}
\end{itemize}

It then verifies that the correct
values are obtained for some simple test cases \textbf{:TODO:}.  For each
successful test (both of these valid data and the invalid ones
described above), it prints ``\texttt{PASS}'' to standard output; if a
test fails, it prints ``\texttt{FAIL}''.

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

#include "CheckStatus.h"

NRCSID(CCOARSEGRAINFREQUENCYSERIESTESTC, "$Id$");

#define CCOARSEGRAINFREQUENCYSERIESTESTC_TOL           1e-8

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
     = {{1.0,0.0}, {2.0,0.0}, {3.0,0.0}, {4.0,0.0}, {5.0,0.0}, {6.0,0.0},
	{7.0,0.0}};

   const COMPLEX8 expectedOutput1DataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1] 
                     = {{0.5,0.0}, {2.0,0.0}, {4.0,0.0}, {6.0,0.0}};

   const COMPLEX8 expectedOutput2DataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2] 
                     = {{(2.0/3.0),0.0}, {3.0,0.0}, {6.0,0.0}};

   COMPLEX8FrequencySeries             goodInput, badInput;
   COMPLEX8FrequencySeries     goodOutput, badOutput;

   BOOLEAN                result;
   LALUnitPair            unitPair;
   LALUnit                dimensionless = { 0 };
   CHARVector             *unitString;

   FrequencySamplingParams     params;

   ParseOptions( argc, argv );

   /* TEST INVALID DATA HERE ------------------------------------------- */

   /* define valid parameters */
   goodInput.f0                   = CCOARSEGRAINFREQUENCYSERIESTESTC_F00;
   goodInput.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   goodInput.epoch.gpsSeconds     = CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC;
   goodInput.epoch.gpsNanoSeconds = CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS;
   goodInput.data                 = NULL;
   goodOutput.data                = NULL;

   params.f0                      = CCOARSEGRAINFREQUENCYSERIESTESTC_F01;
   params.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
   params.length               = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

   badInput = goodInput;
   badOutput = goodOutput;

   /* allocate input and output vectors */
   LALCCreateVector(&status, &(goodInput.data),
		    CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
   if ( code = CheckStatus(&status, 0 , "",
                           CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                           CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
   {
     return code;
   }
   LALCCreateVector(&status, &(goodOutput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
   if ( code = CheckStatus(&status, 0 , "",
                           CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                           CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
   {
     return code;
   }

#ifndef LAL_NDEBUG
   if ( ! lalNoDebug )
   {
     /* test behavior for null pointer to output series */
     LALCCoarseGrainFrequencySeries(&status, NULL, &goodInput, &params);
     if ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR, 
                             COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK )) 
     {
       return code;
     }
     printf("  PASS: null pointer to output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to input series */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, NULL, &params);
     if ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR, 
                             COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK )) 
     {
       return code;
     }
     printf("  PASS: null pointer to input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
   
     /* test behavior for null pointer to plan parameter */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, NULL);
     if ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR, 
                             COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK )) 
     {
       return code;
     }
     printf("  PASS: null pointer to plan parameter results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
   
     /* test behavior for null pointer to data member of output series */
     LALCCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR, 
                             COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK )) 
     {
       return code;
     }
     printf("  PASS: null pointer to data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of input series */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR, 
                             COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK )) 
     {
       return code;
     }
     printf("  PASS: null pointer to data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of data member of output series */
     LALCCreateVector(&status, &(badOutput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
     if ( code = CheckStatus(&status, 0 , "",
                             CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
     {
       return code;
     }
     cPtr = badOutput.data->data;
     badOutput.data->data = NULL;
     LALCCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR, 
                             COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK )) 
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = cPtr;
     LALCDestroyVector(&status, &(badOutput.data));
     if ( code = CheckStatus(&status, 0 , "",
                             CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
     {
       return code;
     }
     
     /* test behavior for null pointer to data member of data member of output series */
     LALCCreateVector(&status, &(badInput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( code = CheckStatus(&status, 0 , "",
                             CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
     {
       return code;
     }
     cPtr = badInput.data->data;
     badInput.data->data = NULL;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR, 
                             COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK )) 
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badInput.data->data = cPtr;
     LALCDestroyVector(&status, &(badInput.data));
     if ( code = CheckStatus(&status, 0 , "",
                             CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
     {
       return code;
     }

     /* test behavior for zero length
     goodInput.data->length = goodOutput.data->length = 0;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
                             COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK)) 
     {
       return code;
     }
     printf("  PASS: zero length results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);
     */
     /* reassign valid length 
     goodInput.data->length = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;
     goodOutput.data->length = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;
*/
     /* test behavior for negative frequency spacing 
     goodInput.deltaT = -CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAT;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( code = CheckStatus(&status,
                             COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAT,
                             COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAT,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK)) 
     {
       return code;
     }
     printf("  PASS: negative time spacing results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAT);
    */
     /* test behavior for zero time spacing 
     goodInput.deltaT = 0;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( code = CheckStatus(&status,
                             COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAT,
                             COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAT,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK)) 
     {
       return code;
     }
     printf("  PASS: zero time spacing results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAT);
*/
     /* reassign valid time spacing 
     goodInput.deltaT = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAT;
     */

   } /* if ( ! lalNoDebug ) */
#endif /* NDEBUG */

   /* TEST VALID DATA HERE --------------------------------------------- */

   /* fill input time-series parameters */
   strncpy(goodInput.name,"Dummy test data",LALNameLength);
   goodInput.sampleUnits  = dimensionless;

     /* fill input data */
   for (i=0; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0; ++i)
   {
     goodInput.data->data[i] = testInputDataData[i];
   }

   /* coarse grain */
   LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( code = CheckStatus( &status, 0 , "",
			    CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                            CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) )
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
   unitPair.unitOne = goodInput.sampleUnits;
   unitPair.unitTwo = goodOutput.sampleUnits;
   LALUnitCompare(&status, &result, &unitPair);
   if ( code = CheckStatus(&status, 0 , "",
                           CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                           CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
   {
     return code;
   }

   if (optVerbose) 
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( code = CheckStatus(&status, 0 , "",
                             CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
     {
       return code;
     }
    
     LALUnitAsString( &status, unitString, &(unitPair.unitTwo) );
     if ( code = CheckStatus(&status, 0 , "",
                            CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                            CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );
     
     LALUnitAsString( &status, unitString, &(unitPair.unitOne) );
     if ( code = CheckStatus(&status, 0 , "",
                             CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );
     
     LALCHARDestroyVector(&status, &unitString);
     if ( code = CheckStatus(&status, 0 , "",
                             CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
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
	&&
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
	 &&
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
   if ( code = CheckStatus(&status, 0 , "",
			   CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                            CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
   {
     return code;
   }


   /*-------Test #2-------*/

   params.f0                      = CCOARSEGRAINFREQUENCYSERIESTESTC_F02;
   params.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
   params.length               = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2;

   LALCCreateVector(&status, &(goodOutput.data),
		    CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2);
   if ( code = CheckStatus(&status, 0 , "",
                           CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                           CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
   {
     return code;
   }

   /* coarse grain */
   LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( code = CheckStatus( &status, 0 , "",
			    CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                            CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) )
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
   unitPair.unitOne = goodInput.sampleUnits;
   unitPair.unitTwo = goodOutput.sampleUnits;
   LALUnitCompare(&status, &result, &unitPair);
   if ( code = CheckStatus(&status, 0 , "",
                           CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                           CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
   {
     return code;
   }

   if (optVerbose) 
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( code = CheckStatus(&status, 0 , "",
                             CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
     {
       return code;
     }
    
     LALUnitAsString( &status, unitString, &(unitPair.unitTwo) );
     if ( code = CheckStatus(&status, 0 , "",
                            CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                            CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );
     
     LALUnitAsString( &status, unitString, &(unitPair.unitOne) );
     if ( code = CheckStatus(&status, 0 , "",
                             CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );
     
     LALCHARDestroyVector(&status, &unitString);
     if ( code = CheckStatus(&status, 0 , "",
                             CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                             CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
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
	&&
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
	 &&
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
   if ( code = CheckStatus(&status, 0 , "",
			   CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
                            CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) 
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
     if ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
                              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) 
     {
       return code;
     }
     LALCCreateVector(&status, &goodOutput.data, optOutLength);
     if ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
                              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) 
     {
       return code;
     }

     /* Read input file */
     LALCReadFrequencySeries(&status, &goodInput, optInputFile);
     if ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
                              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) 
     {
       return code;
     }
     
     /* coarse grain */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
                              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) 
     {
       return code;
     }

     LALCPrintFrequencySeries(&goodOutput, optOutputFile);
     
     printf("===== Coarse-Graining of User-Specified Series Written to File %s =====\n", optOutputFile);
     
     /* clean up valid data */
     LALCDestroyVector(&status, &goodInput.data);
     if ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
                              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) 
     {
       return code;
     }
     LALCDestroyVector(&status, &goodOutput.data);
     if ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
                              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) 
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
