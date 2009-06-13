/*
*  Copyright (C) 2007 Edward Daw, Jolien Creighton
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

/************************* <lalVerbatim file="SlopeDetectorFilterTestCV">
Author: Daw, E. J.
$Id$
************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{SlopeDetectorFilterTest.c}}

This test program ensures that all errors are flagged correctly.

\subsubsection*{Usage}
\begin{verbatim}
SlopeDetectorFilterTest
\end{verbatim}

\subsubsection*{Description}

\subsubsection*{Exit codes}
\input{SlopeDetectorFilterTestCE}

\subsubsection*{Uses}
\begin{verbatim}
LALSlopeDetectorFilter()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{SlopeDetectorFilterTestCV}}
******************************************************* </lalLaTeX> */

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */

#include <math.h>

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/LALStdlib.h>
#include <lal/SlopeDetectorFilter.h>
#include <lal/AVFactories.h>

/******* DEFINE RCS ID STRING ************/

NRCSID( SLOPEDETECTORFILTERTESTC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

/********************* <lalErrTable file="SlopeDetectorFilterTestCE"> */
#define SLOPEDETECTORFILTERTESTC_EOK 0
#define SLOPEDETECTORFILTERTESTC_ENOM 1
#define SLOPEDETECTORFILTERTESTC_ECHK 2
#define SLOPEDETECTORFILTERTESTC_EFLS 3
#define SLOPEDETECTORFILTERTESTC_ECOR 4
#define SLOPEDETECTORFILTERTESTC_EBADVECTOR 5
#define SLOPEDETECTORFILTERTESTC_EBADFREE 6
#define SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT 7

#define SLOPEDETECTORFILTERTESTC_MSGEOK "Test finished OK"
#define SLOPEDETECTORFILTERTESTC_MSGENOM "Error checking failed for bad input data"
#define SLOPEDETECTORFILTERTESTC_MSGECHK "Error checking failed for bad output data pointer"
#define SLOPEDETECTORFILTERTESTC_MSGEFLS "Error checking failed for input data shorter than number of filter taps"
#define SLOPEDETECTORFILTERTESTC_MSGECOR "Returned unexpected error upon calling function with good data"
#define SLOPEDETECTORFILTERTESTC_MSGEBADVECTOR "Error making input vector"
#define SLOPEDETECTORFILTERTESTC_MSGEBADFREE "Error freeing allocated vector"
#define SLOPEDETECTORFILTERTESTC_MSGEINCOMPATIBLEOUTPUT "Filter output differs from expected"

/***************************** </lalErrTable> */

/* might also wish to define parameters and expected results for test
   cases here, for example */

#define SLOPEDETECTORFILTERTESTC_INPUTVECTORLENGTH 8
#define SLOPEDETECTORFILTERTESTC_OUTPUTVECTORLENGTH 5
#define SLOPEDETECTORFILTERTESTC_NTAPS 4
#define SLOPEDETECTORFILTERTESTC_REDUCELENGTH 5
#define SLOPEDETECTORFILTERTESTC_EQUALITYTOLERANCE 0.00001

/******* DECLARE AND SET GLOBAL lalDebugLevel ************/

int lalDebugLevel = LALMSGLVL3;

/* See the section (currently 7.4.1) of the LSD
 * on "Status-reporting objects" for list of predefined debug levels */

int main( void )
{
  static LALStatus     status;
  REAL4Vector          *output;
  REAL4Vector          *input;

  /** VARIABLES AND DUMMY INPUT DATA FOR TESTING THE ALGORITHM **/
  REAL4  testdata[SLOPEDETECTORFILTERTESTC_INPUTVECTORLENGTH] =
  {-7.3, 2.6, 27.2, 13.0, 12.9, -100.4, -13.2, 4.4};
  REAL4  slopedetectexpout[SLOPEDETECTORFILTERTESTC_OUTPUTVECTORLENGTH] =
  {8.55, 1.67, -38.29, -19.19, 6.17};
  REAL4  offsetdetectexpout[SLOPEDETECTORFILTERTESTC_OUTPUTVECTORLENGTH] =
  {-3.9500, 11.4200, 45.6100, 6.8600, -33.3300};
  REAL4  alfdetectexpout[SLOPEDETECTORFILTERTESTC_OUTPUTVECTORLENGTH] =
  {1580.035278, 986.760803, 17376.673828, 8603.057617, 1472.700439};
  REAL4  boxcartapsexpout[SLOPEDETECTORFILTERTESTC_NTAPS] = {0.25, 0.25, 0.25, 0.25};
  REAL4  gaussiantapsexpout[SLOPEDETECTORFILTERTESTC_NTAPS] =
  {0.011109, 0.606531, 0.606531, 0.011109};
  REAL4  sinetapsexpout[SLOPEDETECTORFILTERTESTC_NTAPS] = {0, -0.866025, 0.866025, 0};
  REAL4  convfilterexpout[SLOPEDETECTORFILTERTESTC_OUTPUTVECTORLENGTH] =
  {18.137936, 24.554724,  14.895966, -53.073654, -68.709694};
  REAL4  sampper = 0.00390625;  /* 1/256 */
  UINT4  i;
  UINT4  ntaps = (UINT4)SLOPEDETECTORFILTERTESTC_NTAPS;
  /* UINT4  statustest; */
  SLOPEFilterParams fparams;

  /*******  CREATE A TEST VECTOR FOR DATA INPUT  ************/

  input = NULL;
  /* create an input vector */
  LALSCreateVector(&status, &input,
		   (UINT4)SLOPEDETECTORFILTERTESTC_INPUTVECTORLENGTH);
  /* check that the input vector was created OK */
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode, status.statusDescription);
    return SLOPEDETECTORFILTERTESTC_EBADVECTOR;
  }

  /********  FILL THE INPUT DATA WITH A TEST SET *******/

  input->length = (UINT4)SLOPEDETECTORFILTERTESTC_INPUTVECTORLENGTH ;

  for(i=0; i < (UINT4)SLOPEDETECTORFILTERTESTC_INPUTVECTORLENGTH ; ++i) {
    input->data[i] = testdata[i];
  }

  /*******  CREATE A TEST VECTOR FOR DATA OUTPUT  **********/

  output = NULL;
  /* create an output vector */
  LALSCreateVector(&status, &output,
		   (UINT4)SLOPEDETECTORFILTERTESTC_OUTPUTVECTORLENGTH);
  /* check that the output vector was created OK */
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode, status.statusDescription);
    return SLOPEDETECTORFILTERTESTC_EBADVECTOR;
  }

  /*******************************************************************************/
  /****** TESTS OF ERROR HANDLING BY THE LALSlopeDetectorFilter() ALGORITHM ******/
  /*******************************************************************************/

  /*******  TEST THAT APPLYING A FILTER WITH NULL INPUT DATA RETURNS AN ERROR  *********/

  /* check error test for bad input data */
  LALSlopeDetectorFilter(&status, output, NULL, 1);
  if( (status.statusCode != SLOPEDETECTORFILTERH_EINPUTNULLP) ||
      strcmp(status.statusDescription, SLOPEDETECTORFILTERH_MSGEINPUTNULLP) ) {
    printf("Error trap on null pointer for input data failed\n");
    return SLOPEDETECTORFILTERTESTC_ENOM;
  }

  /*******  TEST THAT APPLYING A FILTER WITH NULL OUTPUT DATA RETURNS AN ERROR  *******/

  /* check error test for bad output data */
  LALSlopeDetectorFilter(&status, NULL, input, 1);
  if( (status.statusCode != SLOPEDETECTORFILTERH_EOUTPUTNULLP) ||
      strcmp(status.statusDescription, SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP) ){
    printf("Error trap on null pointer for output data failed\n");
    return SLOPEDETECTORFILTERTESTC_ECHK;
  }

  /*******  TEST THAT IF THE INPUT DATA IS SHORTER THAN THE NUMBER OF INPUT TAPS THEN AN ERROR IS RETURNED ****/

  /* check error test for input data shorter or equal to filter taps */
  input->length -= (UINT4)SLOPEDETECTORFILTERTESTC_REDUCELENGTH ;
  LALSlopeDetectorFilter(&status, output, input, ntaps );
  if( (status.statusCode != SLOPEDETECTORFILTERH_EDATATOOSHORT) ||
      strcmp(status.statusDescription, SLOPEDETECTORFILTERH_MSGEDATATOOSHORT)) {
    printf("Error trap on input data vector shorter than filter failed\n");
    return SLOPEDETECTORFILTERTESTC_EFLS;
  }
  input->length += (UINT4)SLOPEDETECTORFILTERTESTC_REDUCELENGTH ;

  /*******  TEST THAT THE NUMBER OF TAPS IS NOT ZERO  *******/

  LALSlopeDetectorFilter(&status, output, input, (UINT4)0);
  if( (status.statusCode != SLOPEDETECTORFILTERH_EDIVBYZERO) ||
      strcmp(status.statusDescription, SLOPEDETECTORFILTERH_MSGEDIVBYZERO)) {
    printf("Error trap for division by zero failed\n");
    return SLOPEDETECTORFILTERTESTC_EFLS;
  }

  /*******  TEST RESPONSE TO VALID DATA  ************/

  LALSlopeDetectorFilter(&status, output, input, ntaps);
  if ( status.statusCode )
  {
    printf( "Unexpectedly got error code %d and message %s\n",
	    status.statusCode, status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_ECOR;
  }

  /*******  CHECK THAT OUTPUT DATA AGREES WITH THAT CALCULATED IN MATLAB  *******/

  for(i=0;i<output->length;++i) {

    /*this line here for debugging*/
    /*printf("%lf\t%lf\n",fabs( ((double)output->data[i]) ),fabs((double)slopedetectexpout[i]) );*/

    if ( (fabs( (double)output->data[i] / (double)slopedetectexpout[i] - 1.0)) >
	 (double)SLOPEDETECTORFILTERTESTC_EQUALITYTOLERANCE ) {
      printf("Error: In LALSlopeDetectorFilter(), output differs from expected for element %u.\n",i);
      return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
    }
  }
  printf("Output from LALSlopeDetectorFilter() passed tests.\n");

  /*******************************************************************************/
  /************** Tests of LALSlopeConvolutionFilter() ***************************/
  /*******************************************************************************/

  /* Check that boxcar filter taps are correctly assigned */

  fparams.forder = ntaps;
  fparams.tap = LALMalloc(ntaps*sizeof(REAL4));
  fparams.taps_set = LALMalloc(sizeof(UINT4));
  *(fparams.taps_set) = 0;
  fparams.function_select = FILTER_OUTPUT_TAPS_SET_BOXCAR;
  fparams.waveform_offset = 0;

  /* Set taps to boxcar */
  LALSlopeConvolutionFilter(&status,NULL,NULL,fparams);
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode,
	   status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_EBADFREE;
  }

  /* Test that taps come out as expected */
  if(*(fparams.taps_set) != 1) {
    printf("After setting boxcar taps, taps bit not set to one.\n");
    return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
  }
  for(i=0; i<ntaps; ++i) {
    if(fparams.tap[i] != boxcartapsexpout[i]) {
      printf("Boxcar tap %u not as expected.\n",i);
      return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
    }
  }
  printf("Boxcar taps test passed OK.\n");

  /*Check that sinewave filter taps are correctly assigned */

  *(fparams.taps_set) = 0;
  fparams.function_select = FILTER_OUTPUT_TAPS_SET_SINE;

  /* Set taps to sine */
  LALSlopeConvolutionFilter(&status,NULL,NULL,fparams);
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode,
	   status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_EBADFREE;
  }

  /* check that taps come out as expected */
  if(*(fparams.taps_set) != 1) {
    printf("After setting sinewave taps, taps bit not set to one.\n");
    return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
  }
  for(i=0;i<fparams.forder;++i) {
    if(fabs((REAL8)fparams.tap[i] - (REAL8)sinetapsexpout[i]) >
       SLOPEDETECTORFILTERTESTC_EQUALITYTOLERANCE) {
      printf("Sine tap %u error exceeds tolerance.\n",i);
      return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
    }
  }
  printf("Sine taps test passed OK.\n");

  /* Check that Gaussian filter taps are correctly assigned */
  *(fparams.taps_set) = 0;
  fparams.function_select = FILTER_OUTPUT_TAPS_SET_GAUSSIAN;

  /* Set taps to gaussian */
  LALSlopeConvolutionFilter(&status,NULL,NULL,fparams);
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode,
	   status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_EBADFREE;
  }

  /* check that taps come out as expected */
  if(*(fparams.taps_set) != 1) {
    printf("After setting Gaussian taps, taps bit not set to one.\n");
    return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
  }
  for(i=0; i<ntaps; ++i) {
    if(fabs((REAL8)fparams.tap[i]/(REAL8)gaussiantapsexpout[i]-1) >
       SLOPEDETECTORFILTERTESTC_EQUALITYTOLERANCE) {
      printf("Gaussian tap %u error exceeds tolerance.\n",i);
      return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
    }
  }
  printf("Gaussian taps test passed OK.\n");

  /*** Test convolution with preallocated (Gaussian) taps ***/
  fparams.history_allocated = LALMalloc(sizeof(UINT4));
  *(fparams.history_allocated) = 0;
  fparams.history = LALMalloc((fparams.forder - 1)*sizeof(REAL4));
  fparams.function_select = FILTER_OUTPUT_CONVOLVE;
  fparams.sampling_period_s = sampper;
  /* apply convolution filter to well formed data */
  LALSlopeConvolutionFilter(&status,output,input,fparams);
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode,
	   status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_EBADFREE;
  }
  /* print the output data */
  for(i=0;i<output->length;++i) {
    if(fabs((REAL8)output->data[i]/(REAL8)convfilterexpout[i]-1) >
       SLOPEDETECTORFILTERTESTC_EQUALITYTOLERANCE) {
      printf("Convolution result error exceeds tolerance at LALSlopeConvolutionFilter. %d\n",i);
      return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
    }
  }
  printf("Output from LALSlopeConvolutionFilter() passed tests.\n");

  /******************************************************************************/
  /************** Tests of LALSlopeLineFitFilter()     **************************/
  /******************************************************************************/

  /*** Test application of slope filter ***/
  fparams.function_select = FILTER_OUTPUT_SLOPE;
  LALSlopeLineFitFilter(&status,output,input,fparams);
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode,
	   status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_EBADFREE;
  }

  /* print filter outputs */
  for(i=0;i<output->length;++i) {
    /* printf("Slope output %u is %f.\n",i,output->data[i]); */
    if(fabs((REAL8)sampper*(REAL8)output->data[i]/(REAL8)slopedetectexpout[i]-1) >
       SLOPEDETECTORFILTERTESTC_EQUALITYTOLERANCE) {
      printf("Convolution result error exceeds tolerance at LALSlopeLineFitFilter (slope). %d\n",i);
      return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
    }
  }
  printf("Slope fit test passed OK.\n");

  /*** Test application of offset filter ***/
  fparams.function_select = FILTER_OUTPUT_OFFSET;
  LALSlopeLineFitFilter(&status,output,input,fparams);
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode,
	   status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_EBADFREE;
  }

  for(i=0;i<output->length;++i) {
    /* printf("Offset output %u is %f.\n",i,output->data[i]); */
    if(fabs((REAL8)output->data[i]/(REAL8)offsetdetectexpout[i]-1) >
       SLOPEDETECTORFILTERTESTC_EQUALITYTOLERANCE) {
      printf("Convolution error exceeds tolerance at LALSlopeLineFitFilter (offset).\n %d",i);
      return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
      }
  }
  printf("Offset fit test passed OK.\n");

  /*** Test application of alf filter ***/
  fparams.function_select = FILTER_OUTPUT_ALF;
  LALSlopeLineFitFilter(&status,output,input,fparams);
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode,
	   status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_EBADFREE;
  }

  for(i=0;i<output->length;++i) {
    /* printf("ALF output %u is %f.\n",i,output->data[i]); */
    if(fabs((REAL8)output->data[i]/(REAL8)alfdetectexpout[i]-1) >
       SLOPEDETECTORFILTERTESTC_EQUALITYTOLERANCE) {
      printf("Convolution error exceeds tolerance at LALSlopeLineFitFilter (alf).\n %d",i);
      return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
      }
  }
  printf("ALF filter test passed OK.\n");
  printf("Output from LALSlopeLineFitFilter() passed tests.\n");

  /*******  CLEAN UP  ************/

  LALFree(fparams.history);
  LALFree(fparams.history_allocated);
  LALFree(fparams.tap);
  LALFree(fparams.taps_set);

  LALSDestroyVector( &status, &input );
  if( status.statusCode ) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode,
	   status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_EBADFREE;
  }

  LALSDestroyVector( &status, &output );
  if( status.statusCode ) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode,
	   status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_EBADFREE;
  }

  printf("Finished freeing. run memory check.\n");

  LALCheckMemoryLeaks();

  printf("PASS: All tests\n");

  return SLOPEDETECTORFILTERTESTC_EOK;
}










