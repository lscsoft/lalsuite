/*************** <lalVerbatim file="IIRFilterTestCV"> ***************
Author: Creighton, T. D.
$Id$
**************** </lalVerbatim> ***********************************/

/* <lalLaTeX>

\subsection{Program \texttt{IIRFilterTest.c}}
\label{s:IIRFilterTest.c}

Tests the routines in \verb@IIRFilter.h@.

\subsubsection*{Usage}
\begin{verbatim}
IIRFilterTest [-o] [-d [debug-level]]
\end{verbatim}

\subsubsection*{Description}

This program generates a time series vector with an impulse in it, and
passes it through a third-order Butterworth low-pass filter.  By
default, running this program with no arguments simply tests the
subroutines, producing no output.  All filter parameters are set from
\verb@#define@d constants.

The \verb@-o@ flag tells the program to print the input and output
vectors to data files: \verb@out.0@ stores the initial impulse,
\verb@out.1@ the response computed using \verb@IIRFilterREAL4()@,
\verb@out.2@ the response computed using \verb@SIIRFilter()@,
\verb@out.3@ the response from \verb@IIRFilterREAL4Vector()@, and
\verb@out.4@ the response from \verb@IIRFilterREAL4VectorR()@.  The
\verb@-d@ option increases the default debug level from 0 to 1, or
sets it to the specified value \verb@debug-level@.

\subsubsection*{Exit codes}
\begin{tabular}{|c|l|}
\hline
 Code & Explanation                   \\
\hline
\tt 0 & Success, normal exit.         \\
\tt 1 & Subroutine failed.            \\
\tt 2 & Could not create output file. \\
\hline
\end{tabular}

\subsubsection*{Uses}
\begin{verbatim}
debuglevel
LALPrintError()
SCreateVector()
SDestroyVector()
CreateCOMPLEX8ZPGFilter()
DestroyCOMPLEX8ZPGFilter()
WToZCOMPLEX8ZPGFilter()
CreateREAL4IIRFilter()
DestroyREAL4IIRFilter()
SIIRFilter()
IIRFilterREAL4()
IIRFilterREAL4Vector()
IIRFilterREAL4VectorR()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{IIRFilterTestCV}}

</lalLaTeX> */

#include "LALStdlib.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "AVFactories.h"
#include "IIRFilter.h"
#include "ZPGFilter.h"

NRCSID(IIRFILTERTESTC,"$Id$");

/* Default parameters. */
#define NPOINTS 4096 /* Length of time series. */
#define OFFSET  1024 /* Offset of the impulse from the start. */
#define WC      0.01 /* Characteristic frequency in w-plane. */
#define OUTFILE0 "out.0"
#define OUTFILE1 "out.1"
#define OUTFILE2 "out.2"
#define OUTFILE3 "out.3"
#define OUTFILE4 "out.4"

#define IIRFILTERTEST_ESUB  1
#define IIRFILTERTEST_EFILE 2

#define IIRFILTERTEST_MSGESUB  "Subroutine returned error"
#define IIRFILTERTEST_MSGEFILE "File creation error"

#define CHECKSTAT(statusptr)                                         \
do{                                                                  \
  if((statusptr).statusCode){                                        \
    if(debuglevel>0){                                                \
      LALPrintError("%s: %s\n",argv[0],IIRFILTERTEST_MSGESUB);       \
      REPORTSTATUS(&stat);                                           \
    }                                                                \
    return IIRFILTERTEST_ESUB;                                       \
  }                                                                  \
} while(0)

static void PrintVector(FILE *fp, REAL4Vector *vector);

INT4 debuglevel=0;

INT4 main(INT4 argc, CHAR **argv)
{
  Status stat={0};          /* Status pointer for subroutines. */
  BOOLEAN printout=0;       /* Whether output will be printed. */
  INT4 i;                   /* Index counter. */
  REAL4Vector *input1=NULL; /* A time series input vector. */
  REAL4Vector *input2=NULL; /* Another time series input vector. */
  REAL4Vector *output=NULL; /* A time series output vector. */
  REAL4 *data1;             /* Time series data. */
  REAL4 *data2;             /* More time series data. */
  REAL4IIRFilter *iirFilter=NULL;    /* The IIR filter. */
  COMPLEX8ZPGFilter *zpgFilter=NULL; /* The ZPG filter used to
                                        construct the IIR filter. */
  FILE *fp=NULL; /* The output file. */

  /* Parse the input parameters. */
  for(i=1;argc>1;i++,argc--){
    if(!strcmp(argv[i],"-o")){
      printout=1;
    }else if(!strcmp(argv[i],"-d")){
      if((argc>2)&&(argv[i+1][0]!='-')){
	debuglevel=atoi(argv[++i]);
	argc--;
      }else
	debuglevel=1;
    }else
      LALPrintError("%s: Ignoring argument: %s\n",argv[0],argv[i]);
  }

  /* Create the time-domain filter. */
  CreateCOMPLEX8ZPGFilter(&stat,&zpgFilter,0,3);
  CHECKSTAT(stat);
  zpgFilter->poles->data[0].re=WC*sqrt(0.5);
  zpgFilter->poles->data[0].im=WC*sqrt(0.5);
  zpgFilter->poles->data[1].re=0.0;
  zpgFilter->poles->data[1].im=WC;
  zpgFilter->poles->data[2].re=-WC*sqrt(0.5);
  zpgFilter->poles->data[2].im=WC*sqrt(0.5);
  zpgFilter->gain.re=0.0;
  zpgFilter->gain.im=WC*WC*WC;
  WToZCOMPLEX8ZPGFilter(&stat,zpgFilter);
  CHECKSTAT(stat);
  CreateREAL4IIRFilter(&stat,&iirFilter,zpgFilter);
  CHECKSTAT(stat);
  DestroyCOMPLEX8ZPGFilter(&stat,&zpgFilter);
  CHECKSTAT(stat);

  /* Allocate memory for the time series. */
  SCreateVector(&stat,&input1,NPOINTS);
  CHECKSTAT(stat);
  SCreateVector(&stat,&input2,NPOINTS);
  CHECKSTAT(stat);
  SCreateVector(&stat,&output,NPOINTS);
  CHECKSTAT(stat);

  /* Create the input time series. */
  data1=input1->data;
  data2=input2->data;
  for(i=0;i<NPOINTS;i++)
    *(data1++)=*(data2++)=0.0;
  data1=input1->data;
  data2=input2->data;
  data1[OFFSET]=data2[OFFSET]=1.0;
  if(printout){
    if(!(fp=fopen(OUTFILE0,"w"))){
      LALPrintError("%s: %s\n",argv[0],IIRFILTERTEST_MSGEFILE);
      return IIRFILTERTEST_EFILE;
    }
    PrintVector(fp,input1);
  }

  /* Filter the time series using IIRFilterREAL4(). */
  data1=input1->data;
  data2=output->data;
  for(i=0;i<NPOINTS;i++){
    IIRFilterREAL4(&stat,data2++,*(data1++),iirFilter);
    CHECKSTAT(stat);
  }
  if(printout){
    if(!(fp=fopen(OUTFILE1,"w"))){
      LALPrintError("%s: %s\n",argv[0],IIRFILTERTEST_MSGEFILE);
      return IIRFILTERTEST_EFILE;
    }
    PrintVector(fp,output);
  }

  /* Filter the time series using SIIRFilter(). */
  data1=input1->data;
  data2=output->data;
  for(i=0;i<NPOINTS;i++)
    *(data2++)=SIIRFilter(*(data1++),iirFilter);
  if(printout){
    if(!(fp=fopen(OUTFILE2,"w"))){
      LALPrintError("%s: %s\n",argv[0],IIRFILTERTEST_MSGEFILE);
      return IIRFILTERTEST_EFILE;
    }
    PrintVector(fp,output);
  }

  /* Filter the time series using IIRFilterREAL4Vector(). */
  IIRFilterREAL4Vector(&stat,input1,iirFilter);
  CHECKSTAT(stat);
  if(printout){
    if(!(fp=fopen(OUTFILE3,"w"))){
      LALPrintError("%s: %s\n",argv[0],IIRFILTERTEST_MSGEFILE);
      return IIRFILTERTEST_EFILE;
    }
    PrintVector(fp,input1);
  }

  /* Filter the time series using IIRFilterREAL4VectorR(). */
  IIRFilterREAL4VectorR(&stat,input2,iirFilter);
  CHECKSTAT(stat);
  if(printout){
    if(!(fp=fopen(OUTFILE4,"w"))){
      LALPrintError("%s: %s\n",argv[0],IIRFILTERTEST_MSGEFILE);
      return IIRFILTERTEST_EFILE;
    }
    PrintVector(fp,input2);
  }

  /* Free memory and exit. */
  SDestroyVector(&stat,&input1);
  CHECKSTAT(stat);
  SDestroyVector(&stat,&input2);
  CHECKSTAT(stat);
  SDestroyVector(&stat,&output);
  CHECKSTAT(stat);
  DestroyREAL4IIRFilter(&stat,&iirFilter);
  CHECKSTAT(stat);
  return 0;
}


static void PrintVector(FILE *fp, REAL4Vector *vector)
{
  INT4 i=vector->length;
  REAL4 *data=vector->data;

  while(i--)
    fprintf(fp,"%10.3e\n",*(data++));

  return;
}
