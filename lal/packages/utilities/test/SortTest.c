/*----------------------------------------------------------------------- 
 * 
 * File Name: SortTest.c
 * 
 * Author: Creighton, T. D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------*/

/* <lalLaTeX>

\subsection{Program \texttt{SortTest.c}}

A program to test sorting routines.

\subsubsection{Usage}
\begin{verbatim}
SortTest [-s seed] [-d [debug-level]]
\end{verbatim}

\subsubsection{Description}

This test program creates rank and index arrays for an unordered list
of numbers, and then sorts the list.  The data for the list are
generated randomly, and the output is to \verb@stdout@ unless
redirected.  \verb@SortTest@ returns 0 if it executes successfully,
and 1 if any of the subroutines fail.

The \verb@-s@ option sets the seed for the random number generator; if
\verb@seed@ is set to zero (or if no \verb@-s@ option is given) then
the seed is taken from the processor clock.  The \verb@-d@ option
increases the default debug level from 0 to 1, or sets it to the
specified value \verb@debug-level@.


\subsubsection{Algorithm}

\subsubsection{Uses}
\begin{verbatim}
lalDebugLevel
CreateI4Vector()
CreateSVector()
DestroyI4Vector()
DestroySVector()
LALCreateRandomParams()
LALDestroyRandomParams()
LALUniformDeviate()
LALPrintError()
\end{verbatim}

\subsubsection{Notes}

</lalLaTeX> */

#include <lal/LALStdlib.h>
#include <stdlib.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/Sort.h>

NRCSID(SORTTESTC,"$Id$");

#define NPTS 5

#define SORTTEST_ESUB 1
#define SORTTEST_MSGESUB "Subroutine returned error"

INT4 lalDebugLevel=0;

int main(int argc, char **argv)
{
  static LALStatus stat;
  INT4         i;
  INT4         seed=0;
  INT4Vector   *index=NULL;
  REAL4Vector  *data=NULL;
  RandomParams *params=NULL;

  /* Set debug level to the input argument, if any. */
  for(i=1;argc>1;i++,argc--){
    if(!strcmp(argv[i],"-s")){
      if(argc>2){
	seed=atoi(argv[++i]);
	argc--;
      }else
	LALPrintError("%s: Ignoring argument: -s\n",argv[0]);
    }else if(!strcmp(argv[i],"-d")){
      if((argc>2)&&(argv[i+1][0]!='-')){
	lalDebugLevel=atoi(argv[++i]);
	argc--;
      }else
	lalDebugLevel=1;
    }else
      LALPrintError("%s: Ignoring argument: %s\n",argv[0],argv[i]);
  }

  /* Create vectors and random parameters. */
  LALI4CreateVector(&stat,&index,NPTS);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  LALSCreateVector(&stat,&data,NPTS);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  LALCreateRandomParams(&stat,&params,seed);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }

  /* Initialize data. */
  for(i=0;i<NPTS;i++){
    LALUniformDeviate(&stat,data->data+i,params);
    if(stat.statusCode){
      LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
      REPORTSTATUS(&stat);
      return SORTTEST_ESUB;
    }
    data->data[i]*=9.99;
  }

  /* Rank data; print data and ranks. */
  LALSHeapRank(&stat,index,data);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  fprintf(stdout,  " Data    Rank\n");
  for(i=0;i<NPTS;i++)
    fprintf(stdout," %4.2f     %2i\n",data->data[i],index->data[i]);
  fprintf(stdout,"\n");

  /* Index and sort data; print sorted data and indecies. */
  LALSHeapIndex(&stat,index,data);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  LALSHeapSort(&stat,data);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  fprintf(stdout,  "Sorted  Index\n");
  for(i=0;i<NPTS;i++)
    fprintf(stdout," %4.2f     %2i\n",data->data[i],index->data[i]);

  /* Free and clear. */
  LALI4DestroyVector(&stat,&index);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  LALSDestroyVector(&stat,&data);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  LALDestroyRandomParams(&stat,&params);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  return 0;
}
