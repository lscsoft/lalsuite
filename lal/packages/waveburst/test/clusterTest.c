/******************************** <lalVerbatim file="clusterTestCV">
Author: Klimenko, Sergei and Yakushin, Igor
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{clusterTest.c}}

[One-line description of test program]

\subsubsection*{Usage}
\begin{verbatim}
clusterTest
\end{verbatim}

\subsubsection*{Description}

\subsubsection*{Exit codes}
\input{clusterTestCE}

\subsubsection*{Uses}
\begin{verbatim}
clusterTest()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{clusterTestCV}}
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

NRCSID( CLUSTERTESTC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

/***************************** <lalErrTable file="clusterTestCE"> */
#define CLUSTERTESTC_ENOM 0
#define CLUSTERTESTC_ECHK 1
#define CLUSTERTESTC_EFLS 2
#define CLUSTERTESTC_EIO  3
#define CLUSTERTESTC_EUP  4

#define CLUSTERTESTC_MSGENOM "Nominal exit"
#define CLUSTERTESTC_MSGECHK "Error checking failed to catch bad data"
#define CLUSTERTESTC_MSGEFLS "Incorrect answer for valid data"
#define CLUSTERTESTC_MSGEIO  "Input/output error"
#define CLUSTERTESTC_MSGEUP "Unexpected parameter values"

/***************************** </lalErrTable> */

char testsFileName []= "test_cluster.dat";

int lalDebugLevel=1;

#include "wavelet_static.h"
#include "wavelet_test_static.h"


typedef struct
tagClusterParameters
{
  int dummyForNow;
} clusterParameters;

static void readParams(clusterParameters *p, char *fileName);
static void releaseMemory(ClusterWavelet **is, 
			  ClusterWavelet **should,
			  clusterParameters **p);
static void readClusterWaveletIs(ClusterWavelet **w, char *fileName);
static void readClusterWaveletShould(ClusterWavelet **w, char *fileName);

int main( int argc, char *argv[] )
{
  static LALStatus status;
  int errors=0;
  int i;

  ClusterWavelet *is;
  ClusterWavelet *should;
  clusterParameters *p;

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
      printf("Test = cluster, waveletFile=%s, paramsFile=%s, shouldFile=%s\n",
	     testRecord.waveletFileName, testRecord.paramsFileName, testRecord.shouldFileName);
      fflush(stdout);

      p=(clusterParameters*)LALMalloc(sizeof(clusterParameters));
      if(p==NULL)
	{
	  fprintf(stderr,"Cannot allocate memory\n");
	  exit(1);
	}

      readClusterWaveletIs(&is,testRecord.waveletFileName);
      readClusterWaveletShould(&should,testRecord.shouldFileName);
      readParams(p,testRecord.paramsFileName);

      _clusterMain(is);

      testRecord.result=compareWavelets_ALL(is,should,FALSE) > 0 ? TRUE : FALSE;

      errors+=testRecord.result;

      printf("errors=%d\n",errors);fflush(stdout);

      releaseMemory(&is,&should,&p);

      LALCheckMemoryLeaks();

      readTestRecord(tests,&testRecord);
    }
  LALFclose(tests);


/*    for(i=0;i<32;i++) */
/*      { */
/*        printf("index=%d l=%d\n",i,_f2l(5,i)); */
/*      } */

  return errors;
}


static void readParams(clusterParameters *p, char *fileName)
{
/*    FILE *in; */
/*    in=LALOpenDataFile(fileName); */
/*    if(in==NULL) */
/*      { */
/*        printf("Cannot open %s\n",fileName); */
/*        exit(1); */
/*      } */

/*    fscanf(in,"%d %d %d",&p->nc,&p->aura,&p->nF); */
/*    LALFclose(in);   */
}


static void releaseMemory(ClusterWavelet **input, 
			  ClusterWavelet **output, clusterParameters **p)
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
  readWavelet_ALL(w,in);
  LALFclose(in);
}

static void readClusterWaveletIs(ClusterWavelet **w, char *fileName)
{
  FILE *in;
  in=LALOpenDataFile(fileName);
  _createClusterWavelet(w);
  readWavelet_TS_PM(w,in);
  LALFclose(in);
}
