/*
*  Copyright (C) 2007 Jolien Creighton
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

/**
 * \author Mendell, Greg A.
 * \file
 * \ingroup FoldAmplitudes_h
 *
 * ### Program FoldAmplitudesTest.c ###
 *
 * %[One-line description of test program]
 *
 * The test program test each of the error conditions, and then test the output
 * of known input with the expected output.
 *
 * ### Usage ###
 *
 * \code
 * FoldAmplitudesTest
 * \endcode
 *
 * ### Uses ###
 *
 * \code
 * LALFoldAmplitudes()
 * \endcode
 *
 */

/* ****** INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */
#include <math.h>

/* ****** INCLUDE ANY LDAS LIBRARY HEADERS ************/

/* ****** INCLUDE ANY LAL HEADERS ************/
#include <lal/LALStdlib.h>
#include <lal/FoldAmplitudes.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>

/* ****** DEFINE LOCAL CONSTANTS AND MACROS ************/

/**\name Error Codes */
/*@{*/
#define FOLDAMPLITUDESTESTC_ENOM 0
#define FOLDAMPLITUDESTESTC_ECHK 1
#define FOLDAMPLITUDESTESTC_EFLS 2

#define FOLDAMPLITUDESTESTC_MSGENOM "Nominal exit"
#define FOLDAMPLITUDESTESTC_MSGECHK "Error checking failed to catch bad data"
#define FOLDAMPLITUDESTESTC_MSGEFLS "Incorrect answer for valid data"
/*@}*/

/* Define parameters and expected results for test cases here. */
#define FOLDAMPLITUDESTESTC_TOL                 1.0e-3
#define FOLDAMPLITUDESTESTC_LENGTH              10
#define FOLDAMPLITUDESTESTC_BADLENGTH           1
#define FOLDAMPLITUDESTESTC_TIMESTEP            0.01
#define FOLDAMPLITUDESTESTC_NUMBINS             10
#define FOLDAMPLITUDESTESTC_FREQ                33.2
#define FOLDAMPLITUDESTESTC_FREQDOT             2.5
/* ****** DECLARE AND SET GLOBAL lalDebugLevel ************/

/** \cond DONT_DOXYGEN */

/* See the section (currently 7.4.1) of the LSD
 * on "Status-reporting objects" for list of predefined debug levels */

int main( void )
{

  /* Declare inputs of LAL function */
  static LALStatus     status;
  REAL4Vector		*output;
  FoldAmplitudesInput   input;
  FoldAmplitudesInput   badinput;
  FoldAmplitudesParams  param;

  /* Declare other variables */
  REAL4                 f = FOLDAMPLITUDESTESTC_FREQ;  /* A frequency for producing test amplitudes */
  REAL4                 fDot = FOLDAMPLITUDESTESTC_FREQDOT;    /* A frequency for producing test phases */
  REAL4			delT = FOLDAMPLITUDESTESTC_TIMESTEP;       /* A time step for producing test data */
  REAL4			twoPi = (REAL4) LAL_TWOPI;                 /* For phase bins between 0 an 2*pi    */
  REAL4			binRange;                                  /* binMax - binMin */
  INT4                  lengthAmpVec = FOLDAMPLITUDESTESTC_LENGTH; /* length of the vector of amplitudes */
  INT4			i;                                         /* generic integer index */
  INT4			j;                                         /* generic integer index */
  INT4			k;                                         /* generic integer index */
  INT2			gotError = 0;                                /* Set nonzero if error condition occurs */


  /* Allocate memory */
  input.phaseVec = NULL;
  input.amplitudeVec = NULL;
  badinput.phaseVec = NULL;
  badinput.amplitudeVec = NULL;
  output = NULL;
  LALSCreateVector( &status, &input.phaseVec, FOLDAMPLITUDESTESTC_LENGTH  );
  LALSCreateVector( &status, &input.amplitudeVec, FOLDAMPLITUDESTESTC_LENGTH  );
  LALSCreateVector( &status, &badinput.phaseVec, FOLDAMPLITUDESTESTC_BADLENGTH  );
  LALSCreateVector( &status, &badinput.amplitudeVec, FOLDAMPLITUDESTESTC_LENGTH  );
  LALSCreateVector( &status, &output, FOLDAMPLITUDESTESTC_NUMBINS );

  /* ******  TEST RESPONSE TO INVALID DATA  ************/

  /* Test that all the error conditions are correctly detected by the function */

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {

    /* Test NULL output */
    param.numBins = FOLDAMPLITUDESTESTC_NUMBINS;
    param.binMin = 0.0;
    param.binMax = twoPi;

    LALFoldAmplitudes( &status, NULL, &input, &param );

    if ( status.statusCode != FOLDAMPLITUDESH_ENULLP
        || strcmp(status.statusDescription, FOLDAMPLITUDESH_MSGENULLP) )
    {
      printf( "Got error code %d and message %s\n",
          status.statusCode, status.statusDescription );
      printf( "Expected error code %d and message %s\n",
          FOLDAMPLITUDESH_ENULLP, FOLDAMPLITUDESH_MSGENULLP );
      return FOLDAMPLITUDESTESTC_ECHK;
    }

    /* Test input vectors of different length */
    param.numBins = FOLDAMPLITUDESTESTC_NUMBINS;
    param.binMin = 0.0;
    param.binMax = twoPi;

    LALFoldAmplitudes( &status, output, &badinput, &param );

    if ( status.statusCode != FOLDAMPLITUDESH_EVECSIZE
        || strcmp(status.statusDescription, FOLDAMPLITUDESH_MSGEVECSIZE) )
    {
      printf( "Got error code %d and message %s\n",
          status.statusCode, status.statusDescription );
      printf( "Expected error code %d and message %s\n",
          FOLDAMPLITUDESH_EVECSIZE, FOLDAMPLITUDESH_MSGEVECSIZE );
      return FOLDAMPLITUDESTESTC_ECHK;
    }

    /* Test number of bins < 1 */
    param.numBins = 0;
    param.binMin = 0.0;
    param.binMax = twoPi;

    LALFoldAmplitudes( &status, output, &input, &param );

    if ( status.statusCode != FOLDAMPLITUDESH_ENUMBINS
        || strcmp(status.statusDescription, FOLDAMPLITUDESH_MSGENUMBINS) )
    {
      printf( "Got error code %d and message %s\n",
          status.statusCode, status.statusDescription );
      printf( "Expected error code %d and message %s\n",
          FOLDAMPLITUDESH_ENUMBINS, FOLDAMPLITUDESH_MSGENUMBINS );
      return FOLDAMPLITUDESTESTC_ECHK;
    }

    /* Test bin max <= bin min */
    param.numBins = FOLDAMPLITUDESTESTC_NUMBINS;
    param.binMin = 0.0;
    param.binMax = 0.0;

    LALFoldAmplitudes( &status, output, &input, &param );

    if ( status.statusCode != FOLDAMPLITUDESH_EBINSIZE
        || strcmp(status.statusDescription, FOLDAMPLITUDESH_MSGEBINSIZE) )
    {
      printf( "Got error code %d and message %s\n",
          status.statusCode, status.statusDescription );
      printf( "Expected error code %d and message %s\n",
          FOLDAMPLITUDESH_EBINSIZE, FOLDAMPLITUDESH_MSGEBINSIZE );
      return FOLDAMPLITUDESTESTC_ECHK;
    }

    /* Test bin min != 0 */
    param.numBins = FOLDAMPLITUDESTESTC_NUMBINS;
    param.binMin = -0.5;
    param.binMax = 0.5;

    LALFoldAmplitudes( &status, output, &input, &param );

    if ( status.statusCode != FOLDAMPLITUDESH_EBINMIN
        || strcmp(status.statusDescription, FOLDAMPLITUDESH_MSGEBINMIN) )
    {
      printf( "Got error code %d and message %s\n",
          status.statusCode, status.statusDescription );
      printf( "Expected error code %d and message %s\n",
          FOLDAMPLITUDESH_EBINMIN, FOLDAMPLITUDESH_MSGEBINMIN );
      return FOLDAMPLITUDESTESTC_ECHK;
    }

  }
#endif /* LAL_NDEBUG */

  /* ******  TEST RESPONSE TO VALID DATA  ************/

  /* Test that valid data generate the correct answers */

  /* Test 1: Simple test with constant amplitudes and constant phases  */
  /* All the amplitudes should end up in one bin; specifically the (j + 8)/2 % param.numBins bin.  */

  for (k = 0; k < 2; ++k) {
    param.numBins = 4;
    param.binMin = 0.0;
    if (k == 0) {
        param.binMax = 1.0;      /* test phase measured in cycles */
    } else {
        param.binMax = twoPi;   /* test phase measured in radians */
    }
    binRange = param.binMax - param.binMin;

    for (j = -8; j <= 8; ++j) {

        for ( i = 0 ; i < (INT4)input.amplitudeVec->length ; ++i )
        {
                input.amplitudeVec->data[i] = 1.0;
                /* Set the phase:  Add in a multiple of twoPi to check that values bin correctly */
                input.phaseVec->data[i] = j*(binRange/8.0) + j*j*j*binRange + .001;
        }

        /* Initialize the output vector */
        for ( i = 0 ; i < (INT4)output->length ; ++i )
        {
                output->data[i] = 0.0;
        }

        /*
        printf("\n");
        for ( i = 0 ; i < input.amplitudeVec->length ; ++i )
        {
                printf("Index %i  Input Amplitude %f \n",i, input.amplitudeVec->data[i]);
        }

        printf("\n");
        for ( i = 0 ; i < input.phaseVec->length ; ++i )
        {
                printf("Index %i  Input Phase %f \n",i, input.phaseVec->data[i]);
        }
        */

        LALFoldAmplitudes( &status, output, &input, &param );

        if ( status.statusCode )
        {
                printf( "Unexpectedly got error code %d and message %s\n",
                        status.statusCode, status.statusDescription );
                return FOLDAMPLITUDESTESTC_EFLS;
        }

        printf("\n");
        for ( i = 0 ; i < param.numBins ; ++i )
        {
                printf("Constant phase test binRange %g, index %i, Bin %i, Output Amplitude %f \n",binRange,j,i,output->data[i]);
                if (i == (j + 8)/2 % param.numBins) {
                    if (output->data[i] != lengthAmpVec*1.0) {
                        printf ("For binRange %g expected all the amplitudes to fold into bin %i but instead got %g in this bin. \n",binRange, i,output->data[i]);
                        gotError = 1;
                    }
                } else  {
                    if (output->data[i] != 0.0) {
                        printf ("For binRange %g expected no amplitudes to fold into bin %i but instead got %g in this bin. \n",binRange, i,output->data[i]);
                        gotError = 1;
                    }
                }
        }

    }  /* end for (j = -8; j < 8; ++j) */

  } /* end for (k = 0; k < 1; ++k) */

  printf("\n");
  if  (gotError > 0) {
        printf("Test 1 in FoldAmplitudesTest.c Failed. \n");
        return FOLDAMPLITUDESTESTC_EFLS;
  } else {
        printf("Test 1 in FoldAmplitudesTest.c Passed. \n");
  }

  /* Test 2: Constant amplitudes, simple phases */
  /* One should end up in each bin. */

  param.numBins = 10;
  param.binMin = 0.0;
  param.binMax = twoPi;   /* test phase measured in radians */
  binRange = param.binMax - param.binMin;

  for ( i = 0 ; i < (INT4)input.amplitudeVec->length ; ++i )
  {
        input.amplitudeVec->data[i] = 1;
        input.phaseVec->data[i] = twoPi*i/10.0 + .0001;
  }

  /* Initialize the output vector */
  for ( i = 0 ; i < (INT4)output->length ; ++i )
  {
        output->data[i] = 0.0;
  }

  /*
  printf("\n");
  for ( i = 0 ; i < input.amplitudeVec->length ; ++i )
  {
        printf("Index %i  Input Amplitude %f \n",i, input.amplitudeVec->data[i]);
  }

  printf("\n");
  for ( i = 0 ; i < input.phaseVec->length ; ++i )
  {
        printf("Index %i  Input Phase %f \n",i, input.phaseVec->data[i]);
  }
  */

  LALFoldAmplitudes( &status, output, &input, &param );

  if ( status.statusCode )
  {
    printf( "Unexpectedly got error code %d and message %s\n",
            status.statusCode, status.statusDescription );
    return FOLDAMPLITUDESTESTC_EFLS;
  }

  printf("\n");
  for ( i = 0 ; i < param.numBins ; ++i )
  {
        printf("Bin %i  Output Amplitude %f \n",i,output->data[i]);
        if (output->data[i] != 1.0) {
                printf ("In bin %i expected 1 but instead got %g in this bin. \n",i,output->data[i]);
                gotError = 1;
        }

  }

  printf("\n");
  if  (gotError > 0) {
        printf("Test 2 in FoldAmplitudesTest.c Failed. \n");
        return FOLDAMPLITUDESTESTC_EFLS;
  } else {
        printf("Test 2 in FoldAmplitudesTest.c Passed. \n");
  }


  /* Test 3: One amplitude goes into each bin, so that sin(\PHI) should return sin(\PHI) */

  param.numBins = 10;
  param.binMin = 0.0;
  param.binMax = twoPi;   /* test phase measured in radians */
  binRange = param.binMax - param.binMin;

  for ( i = 0 ; i < (INT4)input.amplitudeVec->length ; ++i )
  {
        input.amplitudeVec->data[i] = sin(twoPi*i/10.0 + .0001);
        input.phaseVec->data[i] = twoPi*i/10.0 + .0001;
  }

  /* Initialize the output vector */
  for ( i = 0 ; i < (INT4)output->length ; ++i )
  {
        output->data[i] = 0.0;
  }

  /*
  printf("\n");
  for ( i = 0 ; i < input.amplitudeVec->length ; ++i )
  {
        printf("Index %i  Input Amplitude %f \n",i, input.amplitudeVec->data[i]);
  }

  printf("\n");
  for ( i = 0 ; i < input.phaseVec->length ; ++i )
  {
        printf("Index %i  Input Phase %f \n",i, input.phaseVec->data[i]);
  }
  */

  LALFoldAmplitudes( &status, output, &input, &param );

  if ( status.statusCode )
  {
    printf( "Unexpectedly got error code %d and message %s\n",
            status.statusCode, status.statusDescription );
    return FOLDAMPLITUDESTESTC_EFLS;
  }

  printf("\n");
  for ( i = 0 ; i < param.numBins ; ++i )
  {
        printf("Bin %i  Output Amplitude %f \n",i,output->data[i]);
        if (fabs(output->data[i] - sin(twoPi*i/10.0 + .0001)) > fabs(output->data[i]*FOLDAMPLITUDESTESTC_TOL)) {
                printf ("In bin %i expected %g but instead got %g in this bin. \n",i,sin(twoPi*i/10.0 + .0001),output->data[i]);
                gotError = 1;
        }

  }

  printf("\n");
  if  (gotError > 0) {
        printf("Test 3 in FoldAmplitudesTest.c Failed. \n");
        return FOLDAMPLITUDESTESTC_EFLS;
  } else {
        printf("Test 3 in FoldAmplitudesTest.c Passed. \n");
  }


  /* Test 4: Two amplitudes go into each bin, such that sin(2*\PHI) should return 2*sin(\PHI) */
  /* Also, phases are input as cycles. */

  param.numBins = 5;
  param.binMin = 0.0;
  param.binMax = 1.0;   /* test phase measured in radians */
  binRange = param.binMax - param.binMin;

  for ( i = 0 ; i < (INT4)input.amplitudeVec->length ; ++i )
  {
        input.amplitudeVec->data[i] = sin(2.0*twoPi*i/5.0 + .001);
        input.phaseVec->data[i] = i/5.0 + .001;
  }

  /* Initialize the output vector */
  for ( i = 0 ; i < (INT4)output->length ; ++i )
  {
        output->data[i] = 0.0;
  }

  printf("\n");
  for ( i = 0 ; i < (INT4)input.amplitudeVec->length ; ++i )
  {
        printf("Index %i  Input Amplitude %f \n",i, input.amplitudeVec->data[i]);
  }

  printf("\n");
  for ( i = 0 ; i < (INT4)input.phaseVec->length ; ++i )
  {
        printf("Index %i  Input Phase in Cycles %f \n",i, input.phaseVec->data[i]);
  }

  LALFoldAmplitudes( &status, output, &input, &param );

  if ( status.statusCode )
  {
    printf( "Unexpectedly got error code %d and message %s\n",
            status.statusCode, status.statusDescription );
    return FOLDAMPLITUDESTESTC_EFLS;
  }

  printf("\n");
  for ( i = 0 ; i < param.numBins ; ++i )
  {
        printf("Bin %i  Output Amplitude %f \n",i,output->data[i]);
        if (fabs(output->data[i] - 2.0*sin(2.0*twoPi*i/5.0 + .001)) > fabs(output->data[i]*FOLDAMPLITUDESTESTC_TOL)) {
                printf ("In bin %i expected %g but instead got %g in this bin. \n",i,2.0*sin(2.0*twoPi*i/5.0 + .001),output->data[i]);
                gotError = 1;
        }

  }

  printf("\n");
  if  (gotError > 0) {
        printf("Test 4 in FoldAmplitudesTest.c Failed. \n");
        return FOLDAMPLITUDESTESTC_EFLS;
  } else {
        printf("Test 4 in FoldAmplitudesTest.c Passed. \n");
  }


  /* Test 5: More realistic data; check output by hand */

  param.numBins = FOLDAMPLITUDESTESTC_NUMBINS;
  param.binMin = 0.0;
  param.binMax = twoPi;

  printf("\n");
  printf("Test 5 check by hand: numBins, binMax = %i, %g \n",param.numBins,param.binMax);

  for ( i = 0 ; i < (INT4)input.amplitudeVec->length ; ++i )
  {
        input.amplitudeVec->data[i] = sin(twoPi*f*i*delT + 0.5*fDot*i*delT*i*delT - 10.001);
        input.phaseVec->data[i] = twoPi*f*i*delT + 0.5*fDot*i*delT*i*delT - 10.001;
  }

  for ( i = 0 ; i < (INT4)output->length ; ++i )
  {
        output->data[i] = 0.0;
  }

  printf("\n");
  for ( i = 0 ; i < (INT4)input.amplitudeVec->length ; ++i )
  {
        printf("Index %i  Input Amplitude %f \n",i, input.amplitudeVec->data[i]);
  }

  printf("\n");
  for ( i = 0 ; i < (INT4)input.phaseVec->length ; ++i )
  {
        printf("Index %i  Input Phase %f \n",i, input.phaseVec->data[i]);
  }

  LALFoldAmplitudes( &status, output, &input, &param );

  if ( status.statusCode )
  {
    printf( "Unexpectedly got error code %d and message %s\n",
            status.statusCode, status.statusDescription );
    return FOLDAMPLITUDESTESTC_EFLS;
  }

  printf("\n");
  for ( i = 0 ; i < param.numBins ; ++i )
  {
        printf("Bin %i  Output Amplitude %f \n",i,output->data[i]);
  }

  printf("\n");
  if  (gotError > 0) {
        printf("Test 5in FoldAmplitudesTest.c Failed. \n");
        return FOLDAMPLITUDESTESTC_EFLS;
  } else {
        printf("Test 5 in FoldAmplitudesTest.c Passed. \n");
  }


  /* some check on the contents of output */

  /* ******  CLEAN UP  ************/
  LALSDestroyVector( &status, &input.phaseVec );
  LALSDestroyVector( &status, &input.amplitudeVec );
  LALSDestroyVector( &status, &badinput.phaseVec );
  LALSDestroyVector( &status, &badinput.amplitudeVec );
  LALSDestroyVector( &status, &output );

  LALCheckMemoryLeaks();

  printf("\n");
  printf("PASS: All tests \n");
  printf("\n");

  return FOLDAMPLITUDESTESTC_ENOM;
}
/** \endcond */
