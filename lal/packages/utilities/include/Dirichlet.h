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

/************************************ <lalVerbatim file="DirichletHV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
*********************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{Dirichlet.h}}
\label{utilities:s:Dirichlet.h}

Provides prototype and error code information for \texttt{LALDirichlet()},
a routine which calculates the values of the Dirichlet kernel
${\cal D}_N(x)$.

\subsection*{Synopsis}
\begin{verbatim}
#include "Dirichlet.h"
\end{verbatim}

\noindent
% Add more information here, if necessary.

\subsection*{Error conditions}
\input{DirichletHErrTable}

\subsection*{Structures}
\begin{verbatim}
struct DirichletParameters
\end{verbatim}
\idx[Type]{DirichletParameters}

\noindent
Contains parameters that specify the Dirichlet kernel $\mathcal{D}_N(x)$.
The fields are:

\begin{description}
\item[\texttt{UINT4 n}] Dirichlet parameter $N$.

\item[\texttt{UINT4 length}] Specified length of output vector.

\item[\texttt{REAL8 deltaX}] Spacing of $x$ values.
\end{description}

\vfill{\footnotesize\input{DirichletHV}}

\newpage\input{DirichletC}
\newpage\input{DirichletTestC}

*********************************************************** </lalLaTeX> */

#ifndef  _DIRICHLET_H
#define  _DIRICHLET_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (DIRICHLETH, "$Id$");

/******************************** <lalErrTable file="DirichletHErrTable"> */

#define DIRICHLETH_ENULLPIN    1
#define DIRICHLETH_ENVALUE     2
#define DIRICHLETH_ESIZE       3
#define DIRICHLETH_EDELTAX     4
#define DIRICHLETH_ENULLPOUT   5
#define DIRICHLETH_ESIZEMM     6
#define DIRICHLETH_ENULLPDOUT  7

#define DIRICHLETH_MSGENULLPIN   "Null pointer to input parameters"
#define DIRICHLETH_MSGENVALUE    "Dirichlet parameter N less than or equal to zero"
#define DIRICHLETH_MSGESIZE      "Length parameter less than or equal to zero"
#define DIRICHLETH_MSGEDELTAX    "Spacing of x values less than or equal to zero"
#define DIRICHLETH_MSGENULLPOUT  "Null pointer to ouput vector"
#define DIRICHLETH_MSGESIZEMM    "Length of data member of output vector does not equal length specified in input parameters"
#define DIRICHLETH_MSGENULLPDOUT "Null pointer to data member of output vector"

/************************************ </lalErrTable> */

typedef struct tagDirichletParameters{
  UINT4	 n;       /* LALDirichlet parameter N */
  UINT4	 length;  /* specified length of output vector */
  REAL8	 deltaX;  /* spacing of x values */
} DirichletParameters;

void
LALDirichlet(LALStatus*,
	     REAL4Vector*,
	     const DirichletParameters*);

#ifdef  __cplusplus
}
#endif /* C++ protection */
#endif /* _DIRICHLET_H */
