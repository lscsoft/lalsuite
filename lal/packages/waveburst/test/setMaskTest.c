/******************************** <lalVerbatim file="setMaskTestCV">
Author: Klimenko, Sergei and Yakushin, Igor
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{setMaskTest.c}}

[One-line description of test program]

\subsubsection*{Usage}
\begin{verbatim}
setMaskTest
\end{verbatim}

\subsubsection*{Description}

\subsubsection*{Exit codes}
\input{setMaskTestCE}

\subsubsection*{Uses}
\begin{verbatim}
setMaskTest()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{setMaskTestCV}}
******************************************************* </lalLaTeX> */

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */
#include <math.h>

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/LALStdlib.h>
#include <lal/LALWavelet.h>
#include <lal/AVFactories.h>
#include <lal/FileIO.h>

/******* DEFINE RCS ID STRING ************/

NRCSID( SETMASKTESTC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

/***************************** <lalErrTable file="setMaskTestCE"> */
#define SETMASKTESTC_ENOM 0
#define SETMASKTESTC_ECHK 1
#define SETMASKTESTC_EFLS 2
#define SETMASKTESTC_EIO  3
#define SETMASKTESTC_EUP  4

#define SETMASKTESTC_MSGENOM "Nominal exit"
#define SETMASKTESTC_MSGECHK "Error checking failed to catch bad data"
#define SETMASKTESTC_MSGEFLS "Incorrect answer for valid data"
#define SETMASKTESTC_MSGEIO  "Input/output error"
#define SETMASKTESTC_MSGEUP "Unexpected parameter values"

/***************************** </lalErrTable> */

char testsFileName []= "test_setMask.dat";

int lalDebugLevel=3;

#include "wavelet_static.h"
#include "wavelet_test_static.h"

typedef struct
tagMaskParameters
{
  int nc;
  BOOLEAN aura;
  int nF;
} maskParameters;

static void readParams(maskParameters *p, char *fileName);
static void releaseMemory(ClusterWavelet **is, 
			  ClusterWavelet **should,
			  maskParameters **p);
static void readClusterWaveletIs(ClusterWavelet **wavelet, char *fileName);
static void readClusterWaveletShould(ClusterWavelet **wavelet, char *fileName);

int main( int argc, char *argv[] )
{
  static LALStatus status;
  int errors=0;

  ClusterWavelet *is;
  ClusterWavelet *should;
  maskParameters *p;

  Test testRecord;
  FILE *tests;

  tests=LALOpenDataFile(testsFileName);
  if(tests==NULL)
    {
      printf("Cannot open %s\n",testsFileName);
      exit(1);
    }
  readTestRecord(tests,&testRecord);
  while(!feof(tests))
    {
      printf("Test = setMask, waveletFile=%s, paramsFile=%s, shouldFile=%s\n",
	     testRecord.waveletFileName, testRecord.paramsFileName, testRecord.shouldFileName);
      fflush(stdout);

      p=(maskParameters*)LALMalloc(sizeof(maskParameters));
      if(p==NULL)
	{
	  fprintf(stderr,"Cannot allocate memory\n");
	  exit(1);
	}

      readClusterWaveletIs(&is,testRecord.waveletFileName);
      readClusterWaveletShould(&should,testRecord.shouldFileName);
      readParams(p,testRecord.paramsFileName);

      if(!p->aura)
	{
	  is->nonZeroFractionAfterCoincidence=1.2/p->nF;
	}
      else
	{
	  is->nonZeroFractionAfterCoincidence=1.0;
	}

      _setMask(is,p->nc,p->aura,is->wavelet);
 
     
      testRecord.result=compareWavelets_MD_TS_PM(is,should,FALSE) > 0 ? TRUE : FALSE;
      errors+=testRecord.result;

/*        printf("errors=%d\n",errors); */

      releaseMemory(&is,&should,&p);
      readTestRecord(tests,&testRecord);
    }
  LALFclose(tests);

  LALCheckMemoryLeaks();

  return errors;
}


static void readParams(maskParameters *p, char *fileName)
{
  FILE *in;
  in=LALOpenDataFile(fileName);
  if(in==NULL)
    {
      printf("Cannot open %s\n",fileName);
      exit(1);
    }

  fscanf(in,"%d %d %d",&p->nc,&p->aura,&p->nF);
  LALFclose(in);  
}


static void releaseMemory(ClusterWavelet **input, 
			  ClusterWavelet **output, maskParameters **p)
{
  _freeClusterWavelet(input);
  _freeClusterWavelet(output);
  LALFree(*p);
}

static void readClusterWaveletShould(ClusterWavelet **w, char *fileName)
{
  FILE *in;
  in=LALOpenDataFile(fileName);
  _createClusterWavelet(w);
  readWavelet_TS_PM(w,in);

  LALFclose(in);
}

static void readClusterWaveletIs(ClusterWavelet **w, char *fileName)
{
  FILE *in;
  in=LALOpenDataFile(fileName);
  _createClusterWavelet(w);
  readWavelet_TS(&(*w)->wavelet,in);
  LALFclose(in);
}
