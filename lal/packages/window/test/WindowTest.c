/*----------------------------------------------------------------------- 
 * 
 * File Name: WindowTest.c
 * 
 * Author: Allen, Bruce ballen@dirac.phys.uwm.edu 
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 *   main()
 *
 * SYNOPSIS 
 * 
 * DESCRIPTION 
 *   Test suite for Window generator
 * 
 * DIAGNOSTICS
 *   Writes PASS, FAIL to stdout as tests are passed or failed. 
 *   Returns 0 if all tests passed and -1 if any test failed. 
 * 
 * CALLS
 *   Window
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "LALDatatypes.h"
#include "AVFactories.h"
#include "Window.h"
#include "PrintVector.h"

NRCSID (MAIN, "$Id$");

/* modify this value to get a stack trace of test */
int debuglevel=2;

/* modify this to turn on printing windows into files for checking */
#define PRINT 1

int
check(Status *status,INT4 code,CHAR * message)
{
  if (status->statusCode != code) {
    printf("FAIL: did not recognize %s\n",message);
    return 1;
  }
  else if (strcmp(message,status->statusDescription)) {
    printf("FAIL: incorrect warning message %s not %s\n",status->statusDescription,message);
    return 1;
  }
  return 0;
}


int main()
{

  static Status status;     
  REAL4Vector *vector = NULL;
  REAL4Vector dummy;
  LALWindowParams params;
  WindowType wintype;
  REAL8 testsquares[]=
    {1024.0,     /* rectangular */
    384.0,       /* Hann */
    546.0+2.0/15.0,   /* Welch */
    341.333984375,     /* Bartlett */
    276.1142857152779,   /* Parzen */
    300.357781729967622,   /* Papoulis */
    406.9376};     /* Hamming */

  CreateVector (&status, &vector, 1024);

  /* Test behavior for null parameter block */
  LALWindow(&status,vector,NULL );
  if (check(&status,WINDOW_NULLPARAM,WINDOW_MSGNULLPARAM)) return 1;

  /* Test behavior for null vector block */
  LALWindow(&status,NULL,&params );
  if (check(&status,WINDOW_NULLHANDLE,WINDOW_MSGNULLHANDLE)) return 1;

  /* Test behavior for non-positive length  */
  params.length=0;
  LALWindow(&status,vector,&params);
  if (check(&status,WINDOW_ELENGTH,WINDOW_MSGELENGTH)) return 1;

  /* Test failures for undefined window type on lower and upper bounds */
  params.length=1024;
  params.type=-1;
  LALWindow( &status, vector, &params );
  if (check(&status,WINDOW_TYPEUNKNOWN,WINDOW_MSGTYPEUNKNOWN)) return 1;
  params.type=NumberWindowTypes;
  LALWindow( &status, vector, &params );
  if (check(&status,WINDOW_TYPEUNKNOWN,WINDOW_MSGTYPEUNKNOWN)) return 1;

  params.type=Rectangular;

  /* test that we get an error if the wrong vector length is present */
  dummy.length=1234;
  dummy.data=NULL;
  LALWindow( &status, &dummy, &params );
  if (check(&status,WINDOW_WRONGLENGTH,WINDOW_MSGWRONGLENGTH)) return 1;

  /* test that we get an error if the vector data area null */
  dummy.length=params.length;
  LALWindow( &status, &dummy, &params );
  if (check(&status,WINDOW_NULLDATA,WINDOW_MSGNULLDATA)) return 1;


  /* Test normalizations */
  for (wintype=Rectangular;wintype<=Hamming;wintype++)
  {
    params.type=wintype;
    LALWindow(&status,vector,&params);
    if (fabs(params.sumofsquares-testsquares[(int)wintype])>1.e-5)
    {
      printf("FAIL: Window %s appears incorrect.\n",params.windowname);
      return 1;
    }
    if (PRINT) PrintVector(vector);
  }


  printf("PASS Window()\n");

  DestroyVector (&status, &vector);

  return 0;
}
