/******************************** <lalVerbatim file="WindowTestCV">
Authora:Allen, B. and Brown, D. A.
Revision: $Id$
**************************************************** </lalVerbatim> */
/********************************************************** <lalLaTeX>
\subsection{Program \texttt{WindowTest.c}}
\label{s:WindowTest.c}


Tests the routines in \verb@Window.h@.


\subsubsection*{Usage}
\begin{verbatim}
WindowTest
\end{verbatim}

\subsubsection*{Description}
This program outputs 7 text files, with names
\texttt{PrintVector.000} to \texttt{PrintVector.007}.
Each of these contain one of the window functions show in
Fig.~\ref{f:window} for $N=1024$.   These files
may be viewed for example by using the public domain graphing
program xmgr, and typing:
\begin{verbatim}
xmgr PrintVector.*
\end{verbatim}

The program also tests all error conditions, and checks that the sum of these
windows squared add to the correct values.  If there is an error in execution,
the corresponding error message is printed.

Upon successful completion, the program prints:
\texttt{PASS Window()}. 


\vfill{\footnotesize\input{WindowTestCV}}

********************************************************** </lalLaTeX> */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/Window.h>
#include <lal/PrintVector.h>

NRCSID (MAIN, "$Id$");

/* modify this value to get a stack trace of test */
int lalDebugLevel=2;

/* modify this to turn on printing windows into files for checking */
#define PRINT 1

static int
check(LALStatus *status,INT4 code,const CHAR * message)
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


int main( void )
{

  static LALStatus status;     
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
    406.9376,     /* Hamming */
    375.178192052468 }; /* Kaiser */

    LALCreateVector (&status, &vector, 1024);

#ifndef LAL_NDEBUG
    if ( ! lalNoDebug )
    {
      /* Test behavior for null parameter block */
      LALWindow(&status,vector,NULL );
      if (check(&status,WINDOWH_ENULLPARAM,WINDOWH_MSGENULLPARAM)) return 1;

      /* Test behavior for null vector block */
      LALWindow(&status,NULL,&params );
      if (check(&status,WINDOWH_ENULLHANDLE,WINDOWH_MSGENULLHANDLE)) return 1;

      /* Test behavior for non-positive length  */
      params.length=0;
      LALWindow(&status,vector,&params);
      if (check(&status,WINDOWH_EELENGTH,WINDOWH_MSGEELENGTH)) return 1;

      /* Test failures for undefined window type on lower and upper bounds */
      params.length=1024;
      params.type=-1;
      LALWindow( &status, vector, &params );
      if (check(&status,WINDOWH_ETYPEUNKNOWN,WINDOWH_MSGETYPEUNKNOWN)) return 1;
      params.type=NumberWindowTypes;
      LALWindow( &status, vector, &params );
      if (check(&status,WINDOWH_ETYPEUNKNOWN,WINDOWH_MSGETYPEUNKNOWN)) return 1;

      params.type=Rectangular;

      /* test that we get an error if the wrong vector length is present */
      dummy.length=1234;
      dummy.data=NULL;
      LALWindow( &status, &dummy, &params );
      if (check(&status,WINDOWH_EWRONGLENGTH,WINDOWH_MSGEWRONGLENGTH)) return 1;

      /* test that we get an error if the vector data area null */
      dummy.length=params.length;
      LALWindow( &status, &dummy, &params );
      if (check(&status,WINDOWH_ENULLDATA,WINDOWH_MSGENULLDATA)) return 1;
    }
#endif



    /* Test normalizations */
    for (wintype=Rectangular;wintype<=Kaiser;wintype++)
    {
      params.length=vector->length;
      params.type=wintype;
      params.beta = 6.0;
      LALWindow(&status,vector,&params);
      if (fabs(params.sumofsquares-testsquares[(int)wintype])>1.e-5)
      {
        printf("FAIL: Window %s appears incorrect.\n",params.windowname);
        printf("Expected %16.12f, got %16.12f\n", testsquares[(int)wintype],
            params.sumofsquares );
        return 1;
      }
      if (PRINT) LALPrintVector(vector);
    }


    printf("PASS Window()\n");

    LALDestroyVector (&status, &vector);

    return 0;
}

