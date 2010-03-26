/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton
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

/*  <lalVerbatim file="LALInverse3CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInverse3.c}}


Uses $g^{ij} = 1/(2g) e^{ikl} e^{jab} g_{ka} g_{lb}$ to compute the inverse.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInverse3CP}
\idx{LALInverse3()}

\begin{itemize}
   \item \texttt{inverse,} Output, inverse of the $(3\times 3)$ input matrix
   \item \texttt{*matrix,} Input, matrix whose inverse is required
\end{itemize}
\subsubsection*{Description}


Uses $g^{ij} = 1/(2g) e^{ikl} e^{jab} g_{ka} g_{lb}$ to compute
the inverse; though not efficient, it is good enough for the
3-d matrix that we have. Prevents the need for having a new library.

\subsubsection*{Algorithm}


\subsubsection*{Uses}

\begin{verbatim}
LALDeterminant3
\end{verbatim}

\subsubsection*{Notes}
Do not generalise to more than 3 dimensions.

\vfill{\footnotesize\input{LALInverse3CV}}

</lalLaTeX>  */

#include <lal/LALInspiralBank.h>

NRCSID(LALINVERSE3C, "$Id$");
#define Dim 3

/*  <lalVerbatim file="LALInverse3CP"> */

void LALInverse3(LALStatus *status,
                 REAL8     **inverse,
                 REAL8     **matrix)
{ /* </lalVerbatim> */

   REAL8 epsilon[Dim][Dim][Dim] = {{
                              { 0, 0, 0},
                              { 0, 0, 1},
                              { 0,-1, 0}},
                             {{ 0, 0,-1},
                              { 0, 0, 0},
                              { 1, 0, 0}},
                             {{ 0, 1, 0},
                              {-1, 0, 0},
                              { 0, 0, 0}}};
   REAL8 det,x;
   int i,j,p,q,a,b;

   INITSTATUS(status, "LALInverse3", LALINVERSE3C);
   ATTATCHSTATUSPTR(status);
   ASSERT (inverse,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (matrix,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   LALDeterminant3(status->statusPtr, &det, matrix); CHECKSTATUSPTR(status);

   ASSERT (det!=0,  status, LALINSPIRALBANKH_EDIV0, LALINSPIRALBANKH_MSGEDIV0);

   for (i=0; i<Dim; i++) {
   for (j=0; j<Dim; j++) {
      x = 0;
      for (a=0; a<Dim; a++) { for (b=0; b<Dim; b++) {
      for (p=0; p<Dim; p++) { for (q=0; q<Dim; q++) {
         x+=epsilon[j][a][p] * epsilon[i][b][q] * matrix[a][b] * matrix[p][q];
      }}}}
      inverse[i][j] = x/(2.*det);
   }}
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

#undef Dim
