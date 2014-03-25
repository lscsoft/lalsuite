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

#include <complex.h>
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
 * \file
 * \ingroup CoarseGrainFrequencySeries_h
 * \author UTB Relativity Group; contact whelan@phys.utb.edu
 *
 * \brief Test suite for <tt>LALZCoarseGrainFrequencySeries()</tt>.
 *
 * ### Usage ###
 *
 * \code
 * ./ZCoarseGrainFrequencySeriesTest
 * Options:
 * -h             print usage message
 * -q             quiet: run silently
 * -v             verbose: print extra information
 * -d level       set lalDebugLevel to level
 * -i filename    read fine grained series from file filename
 * -o filename    print coarse grained  series to file filename
 * -n length      input series contains length points
 * -m length      output series contains length points
 * -e deltaF      set coarse grained frequency spacing to deltaF
 * -f f0          set start frequency of output to f0
 * \endcode
 *
 * ### Description ###
 *
 * This program tests the routine
 * <tt>LALZCoarseGrainFrequencySeries()</tt>, which coarse-grains a
 * frequency series.
 *
 * First, it tests that the correct error codes
 * (cf \ref CoarseGrainFrequencySeries_h
 * are generated for the following error conditions (tests in
 * \e italics are not performed if \c LAL_NDEBUG is set, as
 * the corresponding checks in the code are made using the ASSERT macro):
 * <ul>
 * <li> <em>null pointer to output series</em></li>
 * <li> <em>null pointer to input series</em></li>
 * <li> <em>null pointer to data member of output series</em></li>
 * <li> <em>null pointer to data member of input series</em></li>
 * <li> <em>null pointer to data member of data member of input series</em></li>
 * <li> <em>null pointer to data member of data member of output series</em>
 * <li> <em>zero length</em></li>
 * <li> <em>negative frequency spacing</em></li>
 * <li> <em>zero frequency spacing</em></li>
 * </ul>
 *
 * It then verifies that the correct values are obtained for some simple
 * test cases
 * <ul>
 * <li> \f$\{h_\ell'\}=\{0,1,2,3,4,5,6,7\}\f$, \f$f'_0=f_0\f$, \f$\delta f'=1\f$, \f$\delta
 * f=2\f$, \f$N=3\f$; the expected output is \f$\{h_k\}=\{1/2,2,4,6\}\f$.</li>
 * <li> \f$\{h_\ell'\}=\{0,1,2,3,4,5,6,7\}\f$, \f$f'_0=f_0\f$, \f$\delta f'=1\f$, \f$\delta
 * f=3\f$, \f$N=3\f$; the expected output is \f$\{h_k\}=\{2/3,3,6\}\f$.</li>
 * <li> \f$f_0'=40\f$, \f$\delta f'= 1\f$,
 * \f$\{h_k\}=\{f_k+i\,f_k^{-1}|k=0,\ldots,4\}\f$, \f$f_0=41\f$,
 * \f$f_0=f'_0\f$, \f$\delta f=2\f$.
 * \f$\delta f'=3\f$, \f$N=\f$ ; the expected output is
 * \f[
 * \{h'_\ell\}=\left\{41+i\left(\frac{1}{40}+\frac{2}{41}+\frac{1}{42}\right),
 * 43+i\left(\frac{1}{42}+\frac{2}{43}+\frac{1}{44}\right)
 * \right\}
 * \f]</li>
 * </ul>
 * For each successful test (both of these valid data and the invalid
 * ones described above), it prints \c PASS to standard output;
 * if a test fails, it prints \c FAIL.
 *
 * If the \c filename arguments are present, it also reads a
 * frequency series from a file, calls
 * <tt>LALZCoarseGrainFrequencySeries()</tt>, and writes the results to
 * the specified output file.
 *
 * ### Notes ###
 *
 * <ul>
 * <li> In addition to the error checks tested in this routine, the
 * function checks for errors related to inconsistency of coarse
 * graining parameters, as well as duplicate input and output pointers.
 * Tests of these error checks are still to be
 * added to this test program.</li>
 * <li> No specific error checking is done on user-specified data.  If
 * \c length is missing, the resulting default will cause a bad
 * data error.</li>
 * <li> The length of the user-provided series must be specified, even
 * though it could in principle be deduced from the input file, because
 * the data sequences must be allocated before the
 * <tt>LALZReadFrequencySeries()</tt> function is called.</li>
 * <li> If one \c filename argument, but not both, is present,
 * the user-specified data will be silently ignored.</li>
 * </ul>
 *
 */
/**\name Error Codes */
/*@{*/
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_ENOM 0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_EARG 1
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK 2
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS 3
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_EUSE 4
/*@}*/

/** \cond DONT_DOXYGEN */
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGENOM "Nominal exit"
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEARG "Error parsing command-line arguments"
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK "Error checking failed to catch bad data"
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS "Incorrect answer for valid data"
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE "Bad user-entered data"



#define ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL           1e-15

#define ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC      1234
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS       56789

#define ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0    1.0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_F00        0.0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0    8

#define ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1    2.0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_F01        0.0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1    4

#define ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2    3.0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_F02        0.0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2    3

#define ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF3    1.0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_F03        40.0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3    5

#define ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4    2.0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_F04        41.0
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4    2

#define ZCOARSEGRAINFREQUENCYSERIESTESTC_TRUE     1
#define ZCOARSEGRAINFREQUENCYSERIESTESTC_FALSE    0

extern char *optarg;
extern int   optind;

BOOLEAN optVerbose = ZCOARSEGRAINFREQUENCYSERIESTESTC_FALSE;
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

   const COMPLEX16  testInputDataData[ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0]
     = {0.0, 1.0, 2.0, 3.0,
        4.0, 5.0, 6.0, 7.0};

   const COMPLEX16
     expectedOutput1DataData[ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1]
     = {0.5, 2.0, 4.0, 6.0};

   const COMPLEX16
     expectedOutput2DataData[ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2]
     = {(2.0/3.0), 3.0, 6.0};

   const COMPLEX16
     testInput3DataData[ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3]
     = {40.0 + I * 1.0/40.0, 41.0 + I * 1.0/41.0, 42.0 + I * 1.0/42.0,
        43.0 + I * 1.0/43.0, 44.0 + I * 1.0/44.0};

   const COMPLEX16
     expectedOutput4DataData[ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4]
     = {41.0 + I *  (1.0/40.0+2.0/41.0+1.0/42.0) / 4.0,
        43.0 + I *  (1.0/42.0+2.0/43.0+1.0/44.0) / 4.0};

   COMPLEX16FrequencySeries             goodInput;
   COMPLEX16FrequencySeries     goodOutput;

   int                    result;

   FrequencySamplingParams     params;


   ParseOptions( argc, argv );

   /* TEST INVALID DATA HERE ------------------------------------------- */

   /* define valid parameters */
   goodInput.f0                   = ZCOARSEGRAINFREQUENCYSERIESTESTC_F00;
   goodInput.deltaF               = ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   goodInput.epoch.gpsSeconds     = ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC;
   goodInput.epoch.gpsNanoSeconds = ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS;
   goodInput.data                 = NULL;
   goodOutput.data                = NULL;

   params.f0                      = ZCOARSEGRAINFREQUENCYSERIESTESTC_F00;
   params.deltaF               = ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   params.length               = ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;

#ifndef LAL_NDEBUG
   COMPLEX16FrequencySeries badInput = goodInput;
   COMPLEX16FrequencySeries badOutput = goodOutput;
#endif

   /* allocate input and output vectors */
   LALZCreateVector(&status, &(goodInput.data),
                    ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }
   LALZCreateVector(&status, &(goodOutput.data), ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

#ifndef LAL_NDEBUG
   if ( ! lalNoDebug )
   {
     /* test behavior for null pointer to output series */
     LALZCoarseGrainFrequencySeries(&status, NULL, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to input series */
     LALZCoarseGrainFrequencySeries(&status, &goodOutput, NULL, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to parameter structure */
     LALZCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, NULL);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to parameter structure results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of output series */
     LALZCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of input series */
     LALZCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of data member of output series */
     LALZCreateVector(&status, &(badOutput.data), ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     COMPLEX16 *zPtr = badOutput.data->data;
     badOutput.data->data = NULL;
     LALZCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = zPtr;
     LALZDestroyVector(&status, &(badOutput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for null pointer to data member of data member of output series */
     LALZCreateVector(&status, &(badInput.data), ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     zPtr = badInput.data->data;
     badInput.data->data = NULL;
     LALZCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badInput.data->data = zPtr;
     LALZDestroyVector(&status, &(badInput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for duplicate pointers */

     /* input and output series */
     LALZCoarseGrainFrequencySeries(&status, &goodInput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     badOutput = goodInput;
     badOutput.data = goodInput.data;

     /* data members of input and output series */
     LALZCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to data members of input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     badOutput.data = NULL;

     /* data members of data members of input and output series */
     LALZCreateVector(&status, &(badOutput.data),
		      ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     zPtr = badOutput.data->data;
     badOutput.data->data = goodInput.data->data;
     LALZCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to data members of data members of input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = zPtr;
     LALZDestroyVector(&status, &(badOutput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for zero length */
     /* input */

     goodInput.data->length = 0;
     LALZCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
			       COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero length in input results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

     goodInput.data->length = ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;
     goodOutput.data->length = ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

     /* output */

     goodOutput.data->length = params.length = 0;
     LALZCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
			       COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero length in output results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

     goodOutput.data->length = params.length
       = ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;

     /* test behavior for negative frequency spacing */
     goodInput.deltaF = params.deltaF
       = -ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
     LALZCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status,
			       COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
			       COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

     /* test behavior for zero frequency spacing */
     goodInput.deltaF = params.deltaF = 0;
     LALZCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status,
			       COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
			       COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

     /* reassign valid frequency spacing */
     goodInput.deltaF = params.deltaF
       = ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;

   } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

   LALZDestroyVector(&status, &(goodOutput.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* TEST VALID DATA HERE --------------------------------------------- */

   params.f0                      = ZCOARSEGRAINFREQUENCYSERIESTESTC_F01;
   params.deltaF               = ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
   params.length               = ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

   /* allocate input and output vectors */
   LALZCreateVector(&status, &(goodOutput.data), ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* fill input time-series parameters */
   strncpy(goodInput.name,"Dummy test data",LALNameLength);
   goodInput.sampleUnits  = lalDimensionlessUnit;

     /* fill input data */
   for (i=0; i<ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0; ++i)
   {
     goodInput.data->data[i] = testInputDataData[i];
   }

   /* coarse grain */
   LALZCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_F01);
   }
   if (goodOutput.f0)
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1);
   }
   if ( fabs(goodOutput.deltaF-ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1)
        / ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1
        > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   result = XLALUnitCompare( &(goodInput.sampleUnits), &(goodOutput.sampleUnits) );

   if (optVerbose)
   {
     CHAR unitString[LALUnitTextSize];
     if ( XLALUnitAsString( unitString, LALUnitTextSize, &(goodOutput.sampleUnits)) == NULL ) {
       return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
     printf( "Units are \"%s\", ", unitString );
     if ( XLALUnitAsString( unitString, LALUnitTextSize, &(goodInput.sampleUnits)) == NULL ) {
       return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
     printf( "should be \"%s\"\n", unitString );
   }

   if (result != 0)
   {
     printf("  FAIL: Valid data test #1\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%1.15e + %1.15e i\n  should be %1.15e + %1.15e i\n",
            creal(goodOutput.data->data[0]), cimag(goodOutput.data->data[0]),
            creal(expectedOutput1DataData[0]), cimag(expectedOutput1DataData[0]));
   }
   if ((fabs(creal(goodOutput.data->data[0] - expectedOutput1DataData[0]))
        /creal(expectedOutput1DataData[0]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
        ||
       (fabs(cimag(goodOutput.data->data[0] - expectedOutput1DataData[0]))
        /creal(expectedOutput1DataData[0]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
   {
     printf("  FAIL: Valid data test #1\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1; ++i)
   {
     f = i * ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%1.15e + %1.15e i\n  should be %1.15e + %1.15e i\n", f,
              creal(goodOutput.data->data[i]), cimag(goodOutput.data->data[i]),
              creal(expectedOutput1DataData[i]), cimag(expectedOutput1DataData[i]));
     }
     if ((fabs(creal(goodOutput.data->data[i] - expectedOutput1DataData[i]))
          /creal(expectedOutput1DataData[i]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         ||
         (fabs(cimag(goodOutput.data->data[i] - expectedOutput1DataData[i]))
          /creal(expectedOutput1DataData[i]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         )
     {
       printf("  FAIL: Valid data test #1\n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   LALZDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }


   /*-------Test #2-------*/

   params.f0                      = ZCOARSEGRAINFREQUENCYSERIESTESTC_F02;
   params.deltaF               = ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
   params.length               = ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2;

   LALZCreateVector(&status, &(goodOutput.data),
                    ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* coarse grain */
   LALZCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_F02);
   }
   if (goodOutput.f0)
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2);
   }
   if ( fabs(goodOutput.deltaF-ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2)
        / ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2 > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   result = XLALUnitCompare( &(goodInput.sampleUnits), &(goodOutput.sampleUnits) );

   if (optVerbose)
   {
     CHAR unitString[LALUnitTextSize];
     if ( XLALUnitAsString( unitString, LALUnitTextSize, &(goodOutput.sampleUnits)) == NULL ) {
       return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
     printf( "Units are \"%s\", ", unitString );
     if ( XLALUnitAsString( unitString, LALUnitTextSize, &(goodInput.sampleUnits)) == NULL ) {
       return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
     printf( "should be \"%s\"\n", unitString );
   }

   if (result != 0)
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%1.15e + %1.15e i\n  should be %1.15e + %1.15e i\n",
            creal(goodOutput.data->data[0]), cimag(goodOutput.data->data[0]),
            creal(expectedOutput2DataData[0]), cimag(expectedOutput2DataData[0]));
   }
   if ((fabs(creal(goodOutput.data->data[0] - expectedOutput2DataData[0]))
        /creal(expectedOutput2DataData[0]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
        ||
       (fabs(cimag(goodOutput.data->data[0] - expectedOutput2DataData[0]))
        /creal(expectedOutput2DataData[0]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2; ++i)
   {
     f = i * ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%1.15e + %1.15e i\n  should be %1.15e + %1.15e i\n", f,
              creal(goodOutput.data->data[i]), cimag(goodOutput.data->data[i]),
              creal(expectedOutput2DataData[i]), cimag(expectedOutput2DataData[i]));
     }
     if ((fabs(creal(goodOutput.data->data[i] - expectedOutput2DataData[i]))
          /creal(expectedOutput2DataData[i]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         ||
         (fabs(cimag(goodOutput.data->data[i] - expectedOutput2DataData[i]))
          /cimag(expectedOutput2DataData[i]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
     {
       printf("  FAIL: Valid data test #2 \n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   /* clean up valid data */
   LALZDestroyVector(&status, &goodInput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALZDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCheckMemoryLeaks();

   /*-------Test #3-------*/

   goodInput.f0               = ZCOARSEGRAINFREQUENCYSERIESTESTC_F03;
   goodInput.deltaF           = ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF3;

   LALZCreateVector(&status, &(goodInput.data),
                    ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* fill input data */
   for (i=0; i<ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3; ++i)
   {
     goodInput.data->data[i] = testInput3DataData[i];
   }

   params.f0                  = ZCOARSEGRAINFREQUENCYSERIESTESTC_F04;
   params.deltaF              = ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4;
   params.length              = ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4;

   LALZCreateVector(&status, &(goodOutput.data),
                    ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* coarse grain */
   LALZCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_F04);
   }
   if ( fabs(goodOutput.f0-ZCOARSEGRAINFREQUENCYSERIESTESTC_F04)
        / ZCOARSEGRAINFREQUENCYSERIESTESTC_F04 > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4);
   }
   if ( fabs(goodOutput.deltaF-ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4)
        / ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4 > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != ZCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   result = XLALUnitCompare( &(goodInput.sampleUnits), &(goodOutput.sampleUnits) );

   if (optVerbose)
   {
     CHAR unitString[LALUnitTextSize];
     if ( XLALUnitAsString( unitString, LALUnitTextSize, &(goodOutput.sampleUnits)) == NULL ) {
       return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
     printf( "Units are \"%s\", ", unitString );
     if ( XLALUnitAsString( unitString, LALUnitTextSize, &(goodInput.sampleUnits)) == NULL ) {
       return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
     printf( "should be \"%s\"\n", unitString );
   }

   if (result != 0)
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%1.15e + %1.15e i\n  should be %1.15e + %1.15e i\n",
            creal(goodOutput.data->data[0]), cimag(goodOutput.data->data[0]),
            creal(expectedOutput4DataData[0]), cimag(expectedOutput4DataData[0]));
   }
   if ((fabs(creal(goodOutput.data->data[0] - expectedOutput4DataData[0]))
        /creal(expectedOutput4DataData[0]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
        ||
       (fabs(cimag(goodOutput.data->data[0] - expectedOutput4DataData[0]))
        /creal(expectedOutput4DataData[0]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<ZCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4; ++i)
   {
     f = i * ZCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%1.15e + %1.15e i\n  should be %1.15e + %1.15e i\n", f,
              creal(goodOutput.data->data[i]), cimag(goodOutput.data->data[i]),
              creal(expectedOutput4DataData[i]), cimag(expectedOutput4DataData[i]));
     }
     if ((fabs(creal(goodOutput.data->data[i] - expectedOutput4DataData[i]))
          /creal(expectedOutput4DataData[i]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         ||
         (fabs(cimag(goodOutput.data->data[i] - expectedOutput4DataData[i]))
          /creal(expectedOutput4DataData[i]) > ZCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
     {
       printf("  FAIL: Valid data test #3 \n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   /* clean up valid data */
   LALZDestroyVector(&status, &goodInput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALZDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
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

     LALZCreateVector(&status, &goodInput.data, optInLength);
     if ( ( code = CheckStatus( &status, 0 , "",
				ZCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALZCreateVector(&status, &goodOutput.data, optOutLength);
     if ( ( code = CheckStatus( &status, 0 , "",
				ZCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     /* Read input file */
     LALZReadFrequencySeries(&status, &goodInput, optInputFile);
     if ( ( code = CheckStatus( &status, 0 , "",
				ZCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     /* coarse grain */
     LALZCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus( &status, 0 , "",
				ZCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     LALZPrintFrequencySeries(&goodOutput, optOutputFile);

     printf("===== Coarse-Graining of User-Specified Series Written to File %s =====\n", optOutputFile);

     /* clean up valid data */
     LALZDestroyVector(&status, &goodInput.data);
     if ( ( code = CheckStatus( &status, 0 , "",
				ZCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALZDestroyVector(&status, &goodOutput.data);
     if ( ( code = CheckStatus( &status, 0 , "",
				ZCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				ZCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALCheckMemoryLeaks();
   }
   return ZCOARSEGRAINFREQUENCYSERIESTESTC_ENOM;
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
        optVerbose = ZCOARSEGRAINFREQUENCYSERIESTESTC_TRUE;
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
