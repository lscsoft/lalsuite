/*----------------------------------------------------------------------- 
 * 
 * File Name: BandPassTest.c
 * 
 * Author: Creighton, T. D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------*/

/* <lalLaTeX>

\subsection{Program \texttt{BandPassTest.c}}

Tests time-domain high- and low-pass filters.

\subsubsection{Usage}
\begin{verbatim}
BandPassTest [-o [outfile]] [-d [debug-level]]
\end{verbatim}

\subsubsection{Description}

This program generates a time series with an impulse in it, and passes
it through a time-domain low-pass or high-pass filter.  By default,
running this program with no arguments simply tests the subroutines,
producing no output.  All filter parameters are set from
\verb@#define@d constants.  The program returns a value of 0 upon
successful completion, 1 if any of its function calls failed, or 2 if
output file creation failed.

The \verb@-o@ flag tells the program to print the impulse response to
a data file; if \verb@outfile@ is not specified, it will write to the
file \verb@out.dat@.  The \verb@-d@ option increases the default debug
level from 0 to 1, or sets it to the specified value
\verb@debug-level@.

\subsubsection{Algorithm}

\subsubsection{Uses}
\begin{verbatim}
debuglevel
SCreateVector()
SDestroyVector()
ButterworthREAL4TimeSeries()
LALPrintError()
\end{verbatim}

\subsubsection{Notes}

</lalLaTeX> */

#include "LALStdlib.h"
#include <string.h>
#include <stdlib.h>
#include "AVFactories.h"
#include "BandPassTimeSeries.h"

NRCSID(BANDPASSTESTC,"$Id$");

/* Default parameters. */
#define NPOINTS 4096 /* Length of time series. */
#define OFFSET  1024 /* Offset of the impulse from the start. */
#define F1 0.01      /* Lower frequency of transition band. */
#define F2 0.015     /* Upper frequency of transition band. */
#define A1 0.9       /* Desired attenuation at F1. */
#define A2 0.1       /* Desired attenuation at F2. */
#define NMAX 20      /* Maximum filter order. */
#define OUTFILE "out.dat" /* Default output filename. */

#define BANDPASSTEST_ESUB  1
#define BANDPASSTEST_EFILE 2

#define BANDPASSTEST_MSGESUB  "Subroutine returned error"
#define BANDPASSTEST_MSGEFILE "File creation error"

INT4 debuglevel=0;

INT4 main(INT4 argc, CHAR **argv)
{
  Status stat = {0};         /* Status pointer for subroutines. */
  CHAR *fname=NULL;          /* The output filename. */
  INT4 i;                    /* Index counter. */
  REAL4TimeSeries series;    /* Time series. */
  REAL4 *data;               /* Time series data. */
  PassBandParamStruc params; /* Filter parameters. */
  FILE *fp=NULL;             /* Output file. */

  /* Parse the input parameters. */
  for(i=1;argc>1;i++,argc--){
    if(!strcmp(argv[i],"-o")){
      if((argc>2)&&(argv[i+1][0]!='-')){
	fname=argv[++i];
	argc--;
      }else
	fname=OUTFILE;
    }else if(!strcmp(argv[i],"-d")){
      if((argc>2)&&(argv[i+1][0]!='-')){
	debuglevel=atoi(argv[++i]);
	argc--;
      }else
	debuglevel=1;
    }else
      LALPrintError("%s: Ignoring argument: %s\n",argv[0],argv[i]);
  }

  /* Set up the filter parameters. */
  params.f1=F1;
  params.f2=F2;
  params.a1=A1;
  params.a2=A2;
  params.nMax=NMAX;

  /* Create the time series. */
  series.name="Impulse";
  series.epoch.gpsSeconds=0.0;
  series.epoch.gpsNanoSeconds=0.0;
  series.deltaT=1.0;
  series.f0=0.0;
  series.sampleUnits=NULL;
  series.data=NULL;
  SCreateVector(&stat,&(series.data),NPOINTS);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],BANDPASSTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return BANDPASSTEST_ESUB;
  }
  for(data=series.data->data,i=0;i<NPOINTS;data++,i++)
    *data=0.0;
  data=series.data->data;
  data[OFFSET]=1.0;

  /* Filter the time series. */
  ButterworthREAL4TimeSeries(&stat,&series,&params);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],BANDPASSTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return BANDPASSTEST_ESUB;
  }

  /* Print the output, if the -o option was given. */
  if(!fname){
    fp=fopen(OUTFILE,"w");
    if(!fp){
      LALPrintError("%s: %s\n",argv[0],BANDPASSTEST_MSGEFILE);
      return BANDPASSTEST_EFILE;
    }
    for(data=series.data->data,i=0;i<NPOINTS;data++,i++)
      fprintf(fp,"%8.3e\n",*data);
    fclose(fp);
  }

  /* Free memory. */
  SDestroyVector(&stat,&(series.data));
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],BANDPASSTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return BANDPASSTEST_ESUB;
  }

  /* Normal exit. */
  return 0;
}
