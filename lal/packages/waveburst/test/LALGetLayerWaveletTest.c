/******************************** <lalVerbatim file="LALGetLayerWaveletTestCV">
Author: Klimenko, Sergei and Yakushin, Igor
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{LALGetLayerWaveletTest.c}}

[One-line description of test program]

\subsubsection*{Usage}
\begin{verbatim}
LALGetLayerWaveletTest
\end{verbatim}

\subsubsection*{Description}

\subsubsection*{Exit codes}
\input{LALGetLayerWaveletTestCE}

\subsubsection*{Uses}
\begin{verbatim}
LALGetLayerWaveletTest()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALGetLayerWaveletTestCV}}
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

NRCSID( LALGETLAYERWAVELETTESTC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

/***************************** <lalErrTable file="LALGetLayerWaveletTestCE"> */
#define LALGETLAYERWAVELETTESTC_ENOM 0
#define LALGETLAYERWAVELETTESTC_ECHK 1
#define LALGETLAYERWAVELETTESTC_EFLS 2
#define LALGETLAYERWAVELETTESTC_EIO  3
#define LALGETLAYERWAVELETTESTC_EUP  4

#define LALGETLAYERWAVELETTESTC_MSGENOM "Nominal exit"
#define LALGETLAYERWAVELETTESTC_MSGECHK "Error checking failed to catch bad data"
#define LALGETLAYERWAVELETTESTC_MSGEFLS "Incorrect answer for valid data"
#define LALGETLAYERWAVELETTESTC_MSGEIO  "Input/output error"
#define LALGETLAYERWAVELETTESTC_MSGEUP "Unexpected parameter values"

/***************************** </lalErrTable> */

char testsFileName []= "test_getLayer.dat";

int lalDebugLevel=3;

#include "wavelet_static.h"
#include "wavelet_test_static.h"


static void readParams(InputLayerWavelet *w, char *fileName);
static void releaseMemory(InputLayerWavelet **input, 
			  OutputLayerWavelet **output,
			  REAL4TimeSeries **layer);
static void readWavelet(Wavelet **wavelet, char *fileName);

int main( int argc, char *argv[] )
{
  static LALStatus status;
  int errors=0;

  REAL4TimeSeries *LayerShould;
  InputLayerWavelet *input;
  OutputLayerWavelet *output;
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
      printf("Test = LALGetLayerWavelet, waveletFile=%s, paramsFile=%s, shouldFile=%s\n",
	     testRecord.waveletFileName, testRecord.paramsFileName, testRecord.shouldFileName);

      input=(InputLayerWavelet*)LALMalloc(sizeof(InputLayerWavelet));

      readWavelet(&input->wavelet,testRecord.waveletFileName);

      readParams(input,testRecord.paramsFileName);
      readREAL4TimeSeries(&LayerShould,testRecord.shouldFileName);

      LALGetLayerWavelet(&status, &output, input);
      
      if(status.statusCode){
	printf("Error code=%d, message=%s\n", status.statusCode,status.statusDescription);
	return  LALGETLAYERWAVELETTESTC_EFLS;
      }

      testRecord.result=compareREAL4TimeSeries(output->layer,LayerShould,FALSE) > 0 ? TRUE : FALSE;
      errors+=testRecord.result;
      releaseMemory(&input,&output,&LayerShould);
      readTestRecord(tests,&testRecord);
    }
  LALFclose(tests);

  LALCheckMemoryLeaks();

  return errors;
}


static void readParams(InputLayerWavelet *w, char *fileName)
{
  FILE *in;
  in=LALOpenDataFile(fileName);
  if(in==NULL)
    {
      printf("Cannot open %s\n",fileName);
      exit(1);
    }

  fscanf(in,"%d",&w->index);
  LALFclose(in);  
}


static void releaseMemory(InputLayerWavelet **input, 
			  OutputLayerWavelet **output, REAL4TimeSeries **layer)
{
  _freeWavelet(&(*input)->wavelet);
  LALFree((*input));

  _freeREAL4TimeSeries(&(*output)->layer);
  LALFree((*output));

  _freeREAL4TimeSeries(layer);
}

static void readWavelet(Wavelet **wavelet, char *fileName)
{
  FILE *in;
  in=LALOpenDataFile(fileName);
  readWavelet_TS(wavelet,in);
  LALFclose(in);
}
