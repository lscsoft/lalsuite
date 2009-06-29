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

/******************************** <lalVerbatim file="LALCorrelationTestCV">
Author: Yakushin, Igor
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{LALCorrelationTest.c}}

[One-line description of test program]

\subsubsection*{Usage}
\begin{verbatim}
LALCorrelationTest
\end{verbatim}

\subsubsection*{Description}

\subsubsection*{Exit codes}
\input{LALCorrelationTestCE}

\subsubsection*{Uses}
\begin{verbatim}
LALCorrelationTest()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALCorrelationTestCV}}
******************************************************* </lalLaTeX> */

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */
#include <math.h>

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/LALStdlib.h>
#include <lal/LALCorrelation.h>
#include <lal/AVFactories.h>

/******* DEFINE RCS ID STRING ************/

NRCSID( LALCORRELATIONTESTC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

/***************************** <lalErrTable file="LALCorrelationTestCE"> */
#define LALCORRELATIONTESTC_ENOM 0
#define LALCORRELATIONTESTC_ECHK 1
#define LALCORRELATIONTESTC_EFLS 2
#define LALCORRELATIONTESTC_EIO  3
#define LALCORRELATIONTESTC_EUP  4

#define LALCORRELATIONTESTC_MSGENOM "Nominal exit"
#define LALCORRELATIONTESTC_MSGECHK "Error checking failed to catch bad data"
#define LALCORRELATIONTESTC_MSGEFLS "Incorrect answer for valid data"
#define LALCORRELATIONTESTC_MSGEIO  "Input/output error"
#define LALCORRELATIONTESTC_MSGEUP "Unexpected parameter values"

/***************************** </lalErrTable> */


#define LALCORRELATIONTESTC_TOL    1.0e-6

static INT4 readData(InputCorrelation **in, CorrelationParams **p, OutputCorrelation **should);
static INT4 compareOutputs(OutputCorrelation *is, OutputCorrelation *should);

int lalDebugLevel = LALMSGLVL3;
const CHAR fileName[] = "data.txt";

int main( void )
{
  static LALStatus     status;

  OutputCorrelation       *outputIs=NULL;
  OutputCorrelation       *outputShould=NULL;
  InputCorrelation        *input=NULL;
  CorrelationParams       *params=NULL;
  INT4 io;
  INT4 comparison;

  printf("\n\n---- Test of LALCorrelation -----\n\n");

  if((io=readData(&input,&params,&outputShould))!=0){
    printf("%s\n",LALCORRELATIONTESTC_MSGEIO);
    printf("io=%d\n",io);
    return io;
  }


  LALCorrelation(&status, &outputIs, input, params);

  if ( status.statusCode )
  {
    printf( "Unexpectedly got error code %d and message %s\n",
	    status.statusCode, status.statusDescription );
    return LALCORRELATIONTESTC_EFLS;
  }

  if(status.statusCode){
    printf("Error code=%d, message=%s\n", status.statusCode,status.statusDescription);
    return  LALCORRELATIONTESTC_EFLS;
  }

  comparison=compareOutputs(outputIs,outputShould);
  if(comparison){
    printf("Bad computation\n");
    return comparison;             ;
  }

  /*******  CLEAN UP  ************/

  LALFree(params);
  LALFree(input->one->data->data);
  LALFree(input->one->data);
  LALFree(input->one);
  LALFree(input->two->data->data);
  LALFree(input->two->data);
  LALFree(input->two);
  LALFree(input);
  LALFree(outputIs->timeShiftedCorrelation);
  LALFree(outputIs);
  LALFree(outputShould->timeShiftedCorrelation);
  LALFree(outputShould);


  LALCheckMemoryLeaks();

  printf("PASS: All tests\n");

  return 0;
}


static INT4 readData(InputCorrelation **input, CorrelationParams **p, OutputCorrelation **s)
{
  INT4 length=15;
  INT4 shift=3;
  REAL8 deltaT=1.0e7;
  REAL4 dataOne [] = {-0.252, 0.113, 0.307, -0.096, -0.137, -0.142, 0.058, 0.158, 0.256, 0.208, -0.212, 0.331, 0.383, -0.078, 0.034};
  REAL4 dataTwo [] = {-0.357, -0.341, 0.277, -0.218, 0.494, -0.092, 0.436, 0.220, 0.245, -0.187, -0.355, -0.193, -0.494, -0.209, 0.026};
  REAL4 corr [] = {-0.242256, -0.176615, -0.297668, 0.274943, -0.2567, -0.189163, 0.230703};
  INT4 i;


  InputCorrelation *in;
  CorrelationParams *params;
  OutputCorrelation *should;

  in=(InputCorrelation*)LALMalloc(sizeof(InputCorrelation));
  in->one=(REAL4TimeSeries*)LALMalloc(sizeof(REAL4TimeSeries));
  in->two=(REAL4TimeSeries*)LALMalloc(sizeof(REAL4TimeSeries));
  in->one->data=(REAL4Sequence*)LALMalloc(sizeof(REAL4Sequence));
  in->two->data=(REAL4Sequence*)LALMalloc(sizeof(REAL4Sequence));
  params=(CorrelationParams*)LALMalloc(sizeof(CorrelationParams));
  should=(OutputCorrelation*)LALMalloc(sizeof(OutputCorrelation));

  should->length=length;
  in->one->data->length=length;
  in->two->data->length=length;
  in->one->data->data=(REAL4*)LALMalloc(length*sizeof(REAL4));
  in->two->data->data=(REAL4*)LALMalloc(length*sizeof(REAL4));
  strcpy(in->one->name,"one");
  strcpy(in->two->name,"two");

  should->shift=shift;
  should->timeShiftedCorrelation=(REAL4*)LALMalloc((2*shift+1)*sizeof(REAL4));


  should->deltaT=deltaT;
  in->one->deltaT=deltaT;
  in->two->deltaT=deltaT;
  params->maxTimeShiftNan=ceil(deltaT*shift*pow(10,9));

  should->start.gpsSeconds=700000000;
  in->one->epoch.gpsSeconds=should->start.gpsSeconds;
  in->two->epoch.gpsSeconds=should->start.gpsSeconds;

  should->start.gpsNanoSeconds=394730;
  in->one->epoch.gpsNanoSeconds=should->start.gpsNanoSeconds;
  in->two->epoch.gpsNanoSeconds=should->start.gpsNanoSeconds;

  for(i=0;i<length;i++){
    in->one->data->data[i]=dataOne[i];
    in->two->data->data[i]=dataTwo[i];
  }

  for(i=0;i<2*shift+1;i++){
    should->timeShiftedCorrelation[i]=corr[i];
  }

  should->minCorrelationTimeShift=-1;

  should->minCorrelationValue=-0.297668;

  should->maxCorrelationTimeShift=0;

  should->maxCorrelationValue=0.274943;

  *input=in;
  *p=params;
  *s=should;

  return 0;
}


static INT4 compareOutputs(OutputCorrelation *is, OutputCorrelation *should)
{
  INT4 diff=0;
  INT4 i;

  printf("Is:     maxCorrelationTimeShift=%d\n", is->maxCorrelationTimeShift);
  printf("Should: maxCorrelationTimeShift=%d\n\n", should->maxCorrelationTimeShift);

  printf("Is:     maxCorrelationValue=%f\n",is->maxCorrelationValue);
  printf("Should: maxCorrelationValue=%f\n\n",should->maxCorrelationValue);

  printf("Is:     minCorrelationTimeShift=%d\n", is->minCorrelationTimeShift);
  printf("Should: minCorrelationTimeShift=%d\n\n", should->minCorrelationTimeShift);

  printf("Is:     minCorrelationValue=%f\n",is->minCorrelationValue);
  printf("Should: minCorrelationValue=%f\n\n",should->minCorrelationValue);

  printf("Is:     start.gpsSeconds=%d\n",is->start.gpsSeconds);
  printf("Should: start.gpsSeconds=%d\n\n\n",should->start.gpsSeconds);

  printf("Is:     start.gpsNanoSeconds=%d\n",is->start.gpsNanoSeconds);
  printf("Should: start.gpsNanoSeconds=%d\n\n",should->start.gpsNanoSeconds);

  printf("Is:     length=%d\n",is->length);
  printf("Should: length=%d\n\n",should->length);

  printf("Is:     deltaT=%g\n",is->deltaT);
  printf("Should: deltaT=%g\n\n",should->deltaT);

  printf("Is:     shift=%d\n",is->shift);
  printf("Should: shift=%d\n\n",should->shift);

  printf("Correlation comparison\n");
  for(i = 0; i <= 2*is->shift; i++){
    printf("i=%d\n",i);fflush(stdout);
    printf("%f = %f\n" , is->timeShiftedCorrelation[i], should->timeShiftedCorrelation[i]);
    fflush(stdout);
  }

  if(is->maxCorrelationTimeShift!=should->maxCorrelationTimeShift ||
     fabs(is->maxCorrelationValue-should->maxCorrelationValue)>LALCORRELATIONTESTC_TOL ||
     is->minCorrelationTimeShift!=should->minCorrelationTimeShift ||
     fabs(is->minCorrelationValue-should->minCorrelationValue)>LALCORRELATIONTESTC_TOL ||
     is->start.gpsSeconds!=should->start.gpsSeconds ||
     is->start.gpsNanoSeconds!=should->start.gpsNanoSeconds ||
     is->length!=should->length ||
     is->deltaT!=should->deltaT ||
     is->shift!=should->shift){
    return LALCORRELATIONTESTC_EFLS;
  }

  printf("Correlation:\n");
  for(i=0;i<2*is->shift+1;i++)
    {
      if(fabs(is->timeShiftedCorrelation[i]-should->timeShiftedCorrelation[i])>LALCORRELATIONTESTC_TOL)
	{
	  diff++;
	}
    }

  printf("diff=%d\n",diff);

  if(diff>0) return LALCORRELATIONTESTC_EFLS;
  return 0;
}



