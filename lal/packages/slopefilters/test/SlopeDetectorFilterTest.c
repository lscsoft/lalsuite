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

int main( int argc, char *argv[] )
{
  static LALStatus     status;
  REAL4Vector          *output;
  REAL4Vector          *input;

  /** THIS ARRAY CONTAINS DUMMY INPUT DATA FOR TESTING THE ALGORITHM **/
  REAL4  testdata[SLOPEDETECTORFILTERTESTC_INPUTVECTORLENGTH] = {-7.3, 2.6, 27.2, 13.0, 12.9, -100.4, -13.2, 4.4};
  REAL4  slopedetectexpout[SLOPEDETECTORFILTERTESTC_OUTPUTVECTORLENGTH] = {8.55, 1.67, -38.29, -19.19, 6.17};
  UINT4  i;
  UINT4  ntaps = (UINT4)SLOPEDETECTORFILTERTESTC_NTAPS;
  UINT4  statustest;

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

  LALSlopeDetectorFilter(&status, output, input, 0);
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

    if ( (fabs( (double)output->data[i] / (double)slopedetectexpout[i] ) - 1.0) >
	 (double)SLOPEDETECTORFILTERTESTC_EQUALITYTOLERANCE ) {
      printf("Error: In LALSlopeDetectorFilter(), output differs from expected for element %u.\n",i);
      return SLOPEDETECTORFILTERTESTC_EINCOMPATIBLEOUTPUT;
    }
  }
  printf("Output from LALSlopeDetectorFilter() passed tests.\n");

  /*******  CLEAN UP  ************/

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



