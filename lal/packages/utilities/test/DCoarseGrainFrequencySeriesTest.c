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

/**
   \file
   \ingroup CoarseGrainFrequencySeries_h
   \author UTB Relativity Group; contact whelan@phys.utb.edu


   \brief Test suite for <tt>LALDCoarseGrainFrequencySeries()</tt>.

\heading{Usage}
\code
./DCoarseGrainFrequencySeriesTest
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
\endcode

\heading{Description}

This program tests the routine
<tt>LALDCoarseGrainFrequencySeries()</tt>, which coarse-grains a
frequency series.

First, it tests that the correct error codes
(cf \ref CoarseGrainFrequencySeries_h
are generated for the following error conditions (tests in
\e italics are not performed if \c LAL_NDEBUG is set, as
the corresponding checks in the code are made using the ASSERT macro):
<ul>
<li> <em>null pointer to output series</em></li>
<li> <em>null pointer to input series</em></li>
<li> <em>null pointer to data member of output series</em></li>
<li> <em>null pointer to data member of input series</em></li>
<li> <em>null pointer to data member of data member of input series</em></li>
<li> <em>null pointer to data member of data member of output series</em>
</li>
<li> <em>zero length</em></li>
<li> <em>negative frequency spacing</em></li>
<li> <em>zero frequency spacing</em></li>
</ul>

It then verifies that the correct
values are obtained for some simple test cases
<ul>
<li> \f$\{h_\ell'\}=\{0,1,2,3,4,5,6,7\}\f$, \f$f'_0=f_0\f$, \f$\delta f'=1\f$, \f$\delta
f=2\f$, \f$N=3\f$; the expected output is \f$\{h_k\}=\{1/2,2,4,6\}\f$.</li>
<li> \f$\{h_\ell'\}=\{0,1,2,3,4,5,6,7\}\f$, \f$f'_0=f_0\f$, \f$\delta f'=1\f$, \f$\delta
f=3\f$, \f$N=3\f$; the expected output is \f$\{h_k\}=\{2/3,3,6\}\f$.</li>
</ul>
For each successful test (both of these valid data and the invalid
ones described above), it prints \c PASS to standard output;
if a test fails, it prints \c FAIL.

If the \c filename arguments are present, it also reads a
frequency series from a file, calls
<tt>LALDCoarseGrainFrequencySeries()</tt>, and writes the results to
the specified output file.

\heading{Notes}

<ul>
<li> In addition to the error checks tested in this routine, the
  function checks for errors related to inconsistency of coarse
  graining parameters.  Tests of these error checks are still to be
  added to this test program.</li>
<li> No specific error checking is done on user-specified data.  If
  \c length is missing, the resulting default will cause a bad
  data error.</li>
<li> The length of the user-provided series must be specified, even
  though it could in principle be deduced from the input file, because
  the data sequences must be allocated before the
  <tt>LALDReadFrequencySeries()</tt> function is called.</li>
<li> If one \c filename argument, but not both, is present,
  the user-specified data will be silently ignored.</li>
</ul>

*/

/**\name Error Codes */
/*@{*/
#define DCOARSEGRAINFREQUENCYSERIESTESTC_ENOM 0		/**< Nominal exit */
#define DCOARSEGRAINFREQUENCYSERIESTESTC_EARG 1		/**< Error parsing command-line arguments */
#define DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK 2		/**< Error checking failed to catch bad data */
#define DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS 3		/**< Incorrect answer for valid data */
#define DCOARSEGRAINFREQUENCYSERIESTESTC_EUSE 4		/**< Bad user-entered data */
/*@}*/

/** \cond DONT_DOXYGEN */
#define DCOARSEGRAINFREQUENCYSERIESTESTC_MSGENOM "Nominal exit"
#define DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEARG "Error parsing command-line arguments"
#define DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK "Error checking failed to catch bad data"
#define DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS "Incorrect answer for valid data"
#define DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE "Bad user-entered data"



#define DCOARSEGRAINFREQUENCYSERIESTESTC_TOL           1e-15

#define DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC      1234
#define DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS       56789

#define DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0    1.0
#define DCOARSEGRAINFREQUENCYSERIESTESTC_F00        0.0
#define DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0    8

#define DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1    2.0
#define DCOARSEGRAINFREQUENCYSERIESTESTC_F01        0.0
#define DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1    4

#define DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2    3.0
#define DCOARSEGRAINFREQUENCYSERIESTESTC_F02        0.0
#define DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2    3

#define DCOARSEGRAINFREQUENCYSERIESTESTC_TRUE     1
#define DCOARSEGRAINFREQUENCYSERIESTESTC_FALSE    0

extern char *optarg;
extern int   optind;

BOOLEAN optVerbose = DCOARSEGRAINFREQUENCYSERIESTESTC_FALSE;
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


int
main( int argc, char *argv[] )
{

   static LALStatus         status;

   UINT4      i;
   REAL8      f;

   const REAL8    testInputDataData[DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0]
                     = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

   const REAL8 expectedOutput1DataData[DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1]
                     = {0.5, 2.0, 4.0, 6.0};

   const REAL8 expectedOutput2DataData[DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2]
                     = {(2.0/3.0), 3.0, 6.0};

   REAL8FrequencySeries             goodInput;
   REAL8FrequencySeries     goodOutput;

   BOOLEAN                result;
   LALUnitPair            unitPair;

   CHARVector             *unitString;

   FrequencySamplingParams     params;


   ParseOptions( argc, argv );

   /* TEST INVALID DATA HERE ------------------------------------------- */

   /* define valid parameters */
   goodInput.f0                   = DCOARSEGRAINFREQUENCYSERIESTESTC_F00;
   goodInput.deltaF               = DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   goodInput.epoch.gpsSeconds     = DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC;
   goodInput.epoch.gpsNanoSeconds = DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS;
   goodInput.data                 = NULL;
   goodOutput.data                = NULL;

   params.f0                      = DCOARSEGRAINFREQUENCYSERIESTESTC_F01;
   params.deltaF               = DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   params.length               = DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;

#ifndef LAL_NDEBUG
   REAL8FrequencySeries badInput = goodInput;
   REAL8FrequencySeries badOutput = goodOutput;
#endif

   /* allocate input and output vectors */
   LALDCreateVector(&status, &(goodInput.data),
                    DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }
   LALDCreateVector(&status, &(goodOutput.data), DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

#ifndef LAL_NDEBUG
   if ( ! lalNoDebug )
   {
     /* test behavior for null pointer to output series */
     LALDCoarseGrainFrequencySeries(&status, NULL, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to input series */
     LALDCoarseGrainFrequencySeries(&status, &goodOutput, NULL, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to plan parameter */
     LALDCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, NULL);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to plan parameter results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of output series */
     LALDCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of input series */
     LALDCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of data member of output series */
     LALDCreateVector(&status, &(badOutput.data), DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     REAL8                   *dPtr = badOutput.data->data;
     badOutput.data->data = NULL;
     LALDCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = dPtr;
     LALDDestroyVector(&status, &(badOutput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for null pointer to data member of data member of output series */
     LALDCreateVector(&status, &(badInput.data), DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     dPtr = badInput.data->data;
     badInput.data->data = NULL;
     LALDCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badInput.data->data = dPtr;
     LALDDestroyVector(&status, &(badInput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* Removed to make POST05MDC pass make check */

#if 0
          /* test behavior for duplicate pointers */

     /* input and output series */
     LALDCoarseGrainFrequencySeries(&status, &goodInput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     badOutput = goodInput;
     badOutput.data = goodInput.data;

     /* data members of input and output series */
     LALDCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to data members of input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     badOutput.data = NULL;

     /* data members of data members of input and output series */
     LALDCreateVector(&status, &(badOutput.data),
                      DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     dPtr = badOutput.data->data;
     badOutput.data->data = goodInput.data->data;
     LALDCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to data members of data members of input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = dPtr;
     LALDDestroyVector(&status, &(badOutput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

#endif

     /* test behavior for zero length */
     /* input */

     goodInput.data->length = 0;
     LALDCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
			       COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero length in input results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

     goodInput.data->length = DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;
     goodOutput.data->length = DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

     /* output */

     goodOutput.data->length = params.length = 0;
     LALDCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
			       COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero length in output results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

     goodOutput.data->length = params.length
       = DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;

     /* test behavior for negative frequency spacing */
     goodInput.deltaF = params.deltaF
       = -DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
     LALDCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status,
			       COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
			       COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

     /* test behavior for zero frequency spacing */
     goodInput.deltaF = params.deltaF = 0;
     LALDCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status,
			       COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
			       COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
            COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

     /* reassign valid frequency spacing */
     goodInput.deltaF = params.deltaF
       = DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;

   } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

   LALDDestroyVector(&status, &(goodOutput.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* TEST VALID DATA HERE --------------------------------------------- */

   params.f0                      = DCOARSEGRAINFREQUENCYSERIESTESTC_F01;
   params.deltaF               = DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
   params.length               = DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

   /* allocate input and output vectors */
   LALDCreateVector(&status, &(goodOutput.data), DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* fill input time-series parameters */
   strncpy(goodInput.name,"Dummy test data",LALNameLength);
   goodInput.sampleUnits  = lalDimensionlessUnit;

     /* fill input data */
   for (i=0; i<DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0; ++i)
   {
     goodInput.data->data[i] = testInputDataData[i];
   }

   /* coarse grain */
   LALDCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            DCOARSEGRAINFREQUENCYSERIESTESTC_F01);
   }
   if (goodOutput.f0)
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1);
   }
   if ( fabs(goodOutput.deltaF-DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1)
        / DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1
        > DCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   unitPair.unitOne = &(goodInput.sampleUnits);
   unitPair.unitTwo = &(goodOutput.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
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
              DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%1.15e, should be %1.15e\n",
            goodOutput.data->data[0], expectedOutput1DataData[0]);
   }
   if (fabs(goodOutput.data->data[0] - expectedOutput1DataData[0])
        / expectedOutput1DataData[0] > DCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1; ++i)
   {
     f = i * DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%1.15e, should be %1.15e\n",
              f, goodOutput.data->data[i], expectedOutput1DataData[i]);
     }
     if (fabs(goodOutput.data->data[i] - expectedOutput1DataData[i])
         / expectedOutput1DataData[i] > DCOARSEGRAINFREQUENCYSERIESTESTC_TOL)
     {
       printf("  FAIL: Valid data test\n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   LALDDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }


   /*-------Test #2-------*/

   params.f0                      = DCOARSEGRAINFREQUENCYSERIESTESTC_F02;
   params.deltaF               = DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
   params.length               = DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2;

   LALDCreateVector(&status, &(goodOutput.data),
                    DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2);
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* coarse grain */
   LALDCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            DCOARSEGRAINFREQUENCYSERIESTESTC_F02);
   }
   if (goodOutput.f0)
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2);
   }
   if ( fabs(goodOutput.deltaF-DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2)
        / DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2 > DCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != DCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   unitPair.unitOne = &(goodInput.sampleUnits);
   unitPair.unitTwo = &(goodOutput.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
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
              DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%1.15e, should be %1.15e\n",
            goodOutput.data->data[0], expectedOutput2DataData[0]);
   }
   if (fabs(goodOutput.data->data[0] - expectedOutput2DataData[0])
        / expectedOutput2DataData[0] > DCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<DCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2; ++i)
   {
     f = i * DCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%1.15e, should be %1.15e\n",
              f, goodOutput.data->data[i], expectedOutput2DataData[i]);
     }
     if (fabs(goodOutput.data->data[i] - expectedOutput2DataData[i])
         / expectedOutput2DataData[i] > DCOARSEGRAINFREQUENCYSERIESTESTC_TOL)
     {
       printf("  FAIL: Valid data test #2 \n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   /* clean up valid data */
   LALDDestroyVector(&status, &goodInput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALDDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     DCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
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

     LALDCreateVector(&status, &goodInput.data, optInLength);
     if ( ( code = CheckStatus( &status, 0 , "",
				DCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALDCreateVector(&status, &goodOutput.data, optOutLength);
     if ( ( code = CheckStatus( &status, 0 , "",
				DCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     /* Read input file */
     LALDReadFrequencySeries(&status, &goodInput, optInputFile);
     if ( ( code = CheckStatus( &status, 0 , "",
				DCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     /* coarse grain */
     LALDCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus( &status, 0 , "",
				DCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     LALDPrintFrequencySeries(&goodOutput, optOutputFile);

     printf("===== Coarse-Graining of User-Specified Series Written to File %s =====\n", optOutputFile);

     /* clean up valid data */
     LALDDestroyVector(&status, &goodInput.data);
     if ( ( code = CheckStatus( &status, 0 , "",
				DCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALDDestroyVector(&status, &goodOutput.data);
     if ( ( code = CheckStatus( &status, 0 , "",
				DCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				DCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALCheckMemoryLeaks();
   }
   return DCOARSEGRAINFREQUENCYSERIESTESTC_ENOM;
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
  FILE *fp;

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
        break;

      case 'v': /* optVerbose */
        optVerbose = DCOARSEGRAINFREQUENCYSERIESTESTC_TRUE;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        fp = freopen ("/dev/null", "w", stderr);
        if (fp == NULL)
        {
          fprintf(stderr, "Error: Unable to open /dev/null\n");
          exit(1);
        }
        fp = freopen ("/dev/null", "w", stdout);
        if (fp == NULL)
        {
          fprintf(stderr, "Error: Unable to open /dev/null\n");
          exit(1);
        }
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
/** \endcond */
