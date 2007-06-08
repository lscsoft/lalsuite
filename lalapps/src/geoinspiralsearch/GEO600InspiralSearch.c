/*
*  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash
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
