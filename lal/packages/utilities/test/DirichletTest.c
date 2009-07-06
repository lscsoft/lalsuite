/*
*  Copyright (C) 2007 John Whelan
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

/******************************** <lalVerbatim file="DirichletTestCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{DirichletTest.c}}
\label{utilities:ss:DirichletTest.c}

Test suite for LALDirichlet().

\subsubsection*{Usage}
\begin{verbatim}
./DirichletTest
\end{verbatim}

\subsubsection*{Description}

This program tests the function {\tt LALDirichlet()}.
It tests all error conditions listed in the Error codes table.
It also writes to files the values of the Dirichlet kernel for three
different valid test cases.
See Figs.~\ref{utilities:f:dirichlet_fig1}-\ref{utilities:f:dirichlet_fig3}.

% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\begin{figure}[htbp!]
\begin{center}
\noindent\includegraphics[angle=-90,width=4in]{utilitiesDirichletFig1}
\caption{\label{utilities:f:dirichlet_fig1}
Dirichlet kernel for $N=10$, $\Delta x =.01$, and $0\le x\le 1$.}
\end{center}
\end{figure}
%
\begin{figure}[htbp!]
\begin{center}
\noindent\includegraphics[angle=-90,width=4in]{utilitiesDirichletFig2}
\caption{\label{utilities:f:dirichlet_fig2}
Dirichlet kernel for $N=11$, $\Delta x =.01$, and $0\le x\le 1$.}
\end{center}
\end{figure}
%
\begin{figure}[htbp!]
\begin{center}
\noindent\includegraphics[angle=-90,width=4in]{utilitiesDirichletFig3}
\caption{\label{utilities:f:dirichlet_fig3}
Dirichlet kernel for $N=10$, $\Delta x =.01$, and $0\le x\le 2$.}
\end{center}
\end{figure}
%

\subsubsection*{Exit codes}
\begin{tabular}{|c|l|}
\hline
 Code & Explanation                   \\
\hline
\tt 0 & Success, normal exit.         \\
\tt 1 & Subroutine failed.            \\
\hline
\end{tabular}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALDirichlet()
LALPrintVector()
LALSCreateVector()
LALSDestroyVector()
LALCheckMemoryLeaks()
\end{verbatim}

\subsubsection*{Notes}
None.

\vfill{\footnotesize\input{DirichletTestCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <math.h>
#include <string.h>
#include <lal/AVFactories.h>
#include <lal/PrintVector.h>
#include <lal/Dirichlet.h>

/* INT4 lalDebugLevel = LALMSGLVL3; */
int lalDebugLevel = LALNDEBUG;

int check ( LALStatus*, INT4, const CHAR* );

NRCSID(MAIN, "$Id$");

int
main( void )
{

  static LALStatus            status;
  REAL4Vector*                poutput = NULL;
  REAL4Vector                 dummy;
  DirichletParameters         parameters;

  UINT4                        length;
  length = 10;


#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    parameters.length = 11 ;
    dummy.length = 11;
    /* test behavior for null pointer to input parameters */
    LALDirichlet( &status, &dummy, NULL );

    if ( check( &status, DIRICHLETH_ENULLPIN, DIRICHLETH_MSGENULLPIN ) )
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGENULLPIN);

    /* test behavior for LALDirichlet parameter N <= 0  */

    parameters.n = 0;
    LALDirichlet( &status, &dummy, &parameters);
    if ( check( &status, DIRICHLETH_ENVALUE, DIRICHLETH_MSGENVALUE ) )
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGENVALUE);
}
#endif /* LAL_NDEBUG */

  /* define valid value of N */
  parameters.n = 10 ;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    parameters.length = 0;
    LALDirichlet( &status, &dummy, &parameters);
    if ( check( &status, DIRICHLETH_ESIZE, DIRICHLETH_MSGESIZE ) )
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGESIZE);
  }
#endif /* LAL_NDEBUG */

  /* define valid value for specified length of output vector */
  parameters.length = 11 ;
  dummy.length = 11;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for x spacing <= 0 */
    parameters.deltaX = -4.0;
    LALDirichlet( &status, &dummy, &parameters);
    if ( check( &status, DIRICHLETH_EDELTAX, DIRICHLETH_MSGEDELTAX ) )
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGEDELTAX );
    parameters.deltaX = 0.0;
    LALDirichlet( &status, &dummy, &parameters);
    if ( check( &status, DIRICHLETH_EDELTAX, DIRICHLETH_MSGEDELTAX ) )
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGEDELTAX );
  }
#endif  /* LAL_NDEBUG */

  /* define valid delta x */
  parameters.deltaX = 0.1;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to output vector */
    LALDirichlet( &status, NULL, &parameters );
    if ( check( &status, DIRICHLETH_ENULLPOUT, DIRICHLETH_MSGENULLPOUT))
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGENULLPOUT);

    /* test behavior for length of output vector not equal to length  */
    /* specified in input parameters */
    dummy.length = 10;
    LALDirichlet( &status, &dummy, &parameters );
    if ( check( &status, DIRICHLETH_ESIZEMM, DIRICHLETH_MSGESIZEMM ) )
    {
      return 1;
    }
    printf( "PASS: %s\n", DIRICHLETH_MSGESIZEMM );
  }
#endif  /* LAL_NDEBUG */

  /* assign valid output vector length */
  dummy.length = parameters.length;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to data member of output vector */
    dummy.data = NULL;
    LALDirichlet( &status, &dummy, &parameters );
    if ( check( &status, DIRICHLETH_ENULLPDOUT, DIRICHLETH_MSGENULLPDOUT))
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGENULLPDOUT);
  }
# endif  /* LAL_NDEBUG */

  /* VALID TEST DATA #1 */
  /* call Dirichet() with valid data (N=even) */
  parameters.n      = 10;
  parameters.length = 101;
  parameters.deltaX = 0.01;


  LALSCreateVector (&status, &poutput, parameters.length);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALDirichlet( &status, poutput, &parameters );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALPrintVector(poutput);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALSDestroyVector( &status, &poutput );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }

  /* VALID TEST DATA #2 */
  /* call Dirichet() with valid data (N=odd) */
  parameters.n      = 11;
  parameters.length = 101;
  parameters.deltaX = 0.01;


  LALSCreateVector(&status, &poutput, parameters.length);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALDirichlet( &status, poutput, &parameters );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALPrintVector(poutput);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALSDestroyVector( &status, &poutput );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  /* VALID TEST DATA #3 */
  /* call Dirichet() with valid data (x=0 to 2) */
  parameters.n      = 10;
  parameters.length = 201;
  parameters.deltaX = 0.01;


  LALSCreateVector(&status, &poutput, parameters.length);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALDirichlet( &status, poutput, &parameters );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALPrintVector(poutput);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALSDestroyVector( &status, &poutput );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  return 0;
}
/*------------------------------------------------------------------------*/

int
check( LALStatus* status, INT4 code, const CHAR* message )
{
  if ( status->statusCode!= code )
  {
    printf ( "FAIL: did not recognize \"%s\"\n", message );
    return 1;
  }
  else if (code && strcmp( message, status->statusDescription)) {
    printf("FAIL: incorrect warning message \"%s\" not \"%s\"\n",
	   status->statusDescription, message);

    return 1;
  }

  return 0;
}
