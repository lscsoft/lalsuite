/************************* <lalVerbatim file="SlopeDetectorFilterTestCV">
Author: Daw, E. J.
$Id$
************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{SlopeDetectorFilterTest.c}}

[One-line description of test program]

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

#define SLOPEDETECTORFILTERTESTC_MSGEOK "Test finished OK"
#define SLOPEDETECTORFILTERTESTC_MSGENOM "Error checking failed for bad input data"
#define SLOPEDETECTORFILTERTESTC_MSGECHK "Error checking failed for bad output data pointer"
#define SLOPEDETECTORFILTERTESTC_MSGEFLS "Error checking failed for input data shorter than number of filter taps"
#define SLOPEDETECTORFILTERTESTC_MSGECOR "Returned unexpected error upon calling function with good data"
#define SLOPEDETECTORFILTERTESTC_MSGEBADVECTOR "Error making input vector"
#define SLOPEDETECTORFILTERTESTC_MSGEBADFREE "Error freeing allocated vector"
/***************************** </lalErrTable> */

/* might also wish to define parameters and expected results for test
   cases here, for example */

#define SLOPEDETECTORFILTERTESTC_INPUTVECTORLENGTH 100
#define SLOPEDETECTORFILTERTESTC_OUTPUTVECTORLENGTH 100
#define SLOPEDETECTORFILTERTESTC_REDUCELENGTH 20

/******* DECLARE AND SET GLOBAL lalDebugLevel ************/

int lalDebugLevel = LALMSGLVL3;

/* See the section (currently 7.4.1) of the LSD 
 * on "Status-reporting objects" for list of predefined debug levels */

int main( int argc, char *argv[] )
{
  static LALStatus     status;
  REAL4Vector          *output;
  REAL4Vector          *input;

  /*******  CREATE A TEST VECTOR INPUT AND OUTPUT  ************/

  input = NULL;
  /* create an input vector */
  LALSCreateVector(&status, &input, 
		   SLOPEDETECTORFILTERTESTC_INPUTVECTORLENGTH);
  /* check that the input vector was created OK */
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode, status.statusDescription);
    return SLOPEDETECTORFILTERTESTC_EBADVECTOR;
  }

  output = NULL;
  /* create an output vector */
  LALSCreateVector(&status, &output, 
		   SLOPEDETECTORFILTERTESTC_OUTPUTVECTORLENGTH);
  /* check that the output vector was created OK */
  if(status.statusCode) {
    printf("Unexpectedly got error code %d and message %s.\n",
	   status.statusCode, status.statusDescription);
    return SLOPEDETECTORFILTERTESTC_EBADVECTOR;
  }

  /*******  TEST RESPONSE TO INVALID INPUT DATA  *********/

  /* check error test for bad input data */
  LALSlopeDetectorFilter(&status, output, NULL, 1);
  if( (status.statusCode != SLOPEDETECTORFILTERH_EINPUTNULLP) ||
      strcmp(status.statusDescription, SLOPEDETECTORFILTERH_MSGEINPUTNULLP) ) {
    printf("Error trap on null pointer for input data failed");
    return SLOPEDETECTORFILTERTESTC_ENOM;
  }

  /* check error test for bad output data */
  LALSlopeDetectorFilter(&status, NULL, input, 1);
  if( (status.statusCode != SLOPEDETECTORFILTERH_EOUTPUTNULLP) ||
      strcmp(status.statusDescription, SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP) ){
    printf("Error trap on null pointer for output data failed");
    return SLOPEDETECTORFILTERTESTC_ECHK;
  }

  /* check error test for input data shorter or equal to filter taps */
  input->length -= (UINT4)SLOPEDETECTORFILTERTESTC_REDUCELENGTH ;
  LALSlopeDetectorFilter(&status, output, input, 
	     (UINT4)(SLOPEDETECTORFILTERTESTC_INPUTVECTORLENGTH) - 1);
  if( (status.statusCode != SLOPEDETECTORFILTERH_EDATATOOSHORT) ||
      strcmp(status.statusDescription, SLOPEDETECTORFILTERH_MSGEDATATOOSHORT)){
    printf("Error trap on input data vector shorter than filter failed");
    return SLOPEDETECTORFILTERTESTC_EFLS;
  }
  input->length += SLOPEDETECTORFILTERTESTC_REDUCELENGTH ;

  /*******  TEST RESPONSE TO VALID DATA  ************/

  LALSlopeDetectorFilter(&status, output, input, 1);
  if ( status.statusCode ) 
  {
    printf( "Unexpectedly got error code %d and message %s\n",
	    status.statusCode, status.statusDescription );
    return SLOPEDETECTORFILTERTESTC_ECOR;
  }

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

  LALCheckMemoryLeaks();

  printf("PASS: All tests\n");

  return SLOPEDETECTORFILTERTESTC_EOK;
}



