/******************************** <lalVerbatim file="GEO600InspiralSearch">
Author: Sathyaprakash, B. S.
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{GEO600InspiralSearch.c}}
\label{ss:GEO600InspiralSearch.c}

LAL code for inspiral searches.

\subsubsection*{Usage}
\begin{verbatim}
GEO600InspiralSearch
\end{verbatim}

\subsubsection*{Description}

This code filters the GEO600 data through a template bank.
\subsubsection*{Exit codes}
\input{GEO600InspiralSearchCE}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GEO600InspiralSearchCV}}
******************************************************* </lalLaTeX> */

/***************************** <lalErrTable file="GEO600InspiralSearchCE"> */
/***************************** </lalErrTable> */
#include <config.h>
#include <stdio.h>

#if !defined HAVE_GSL_GSL_FFT_REAL_H || !defined HAVE_LIBGSL \
 || !defined HAVE_MYSQL_H || !defined HAVE_LIBMYSQLCLIENT \
 || !defined HAVE_FRAMEL_H || !defined HAVE_LIBFRAME \
 || !defined HAVE_MPI_H

 int main( void ) { return fputs( "disabled\n", stderr ), 77; }

#else

#include <mpi.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>
#include <lal/RealFFT.h>
#include <GEO600InspiralSearch.h>

INT4 lalDebugLevel=1;
     
int 
main (int argc, char * argv[]) 
{
   int myrank;
   static LALStatus status;
   
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if (myrank == 0)
   {
	   GEO600InspiralSearchDoubleMasterDB(&status);
   } else {
	   GEO600InspiralSearchSlave(&status);
   }
   MPI_Finalize();
   LALCheckMemoryLeaks();
   return(0);
}
#endif
