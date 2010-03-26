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

/*  <lalVerbatim file="LALDeterminant3CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALDeterminant.c}}

Module to calculate the determinant of a 3-dimensional matrix $g_{ij}$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALDeterminant3CP}
\idx{LALDeterminant() 3D determinant}
\begin{itemize}
   \item \texttt{determinant,} Output, determinant of the matrix.
   \item \texttt{matrix,} Input, the input $(3\times 3)$ matrix whose determinant is required.
\end{itemize}

\subsubsection*{Description}

This code computes the determinant of a 3-dimensional matrix.

\subsubsection*{Algorithm}
Given a matrix $g_{ij}$ its determinant
is computed using  the formula $g = \epsilon^{ijk} g_{i1} g_{j2} g_{k3},$
where $\epsilon$ is the totally anti-symmetric tensor in 3-dimensions.

\subsubsection*{Uses}
None.

\subsubsection*{Notes}
Don't ever generalise this to higher dimensions since this
would take many more operations than some of the standard routines.

\vfill{\footnotesize\input{LALDeterminant3CV}}

</lalLaTeX>  */

#include <lal/LALInspiralBank.h>

NRCSID(LALDETERMINANT3C, "$Id$");

/*  <lalVerbatim file="LALDeterminant3CP"> */

void LALDeterminant3(LALStatus *status,
                     REAL8     *determinant,
                     REAL8     **matrix)
{ /* </lalVerbatim> */

   REAL8 epsilon[3][3][3] = {{
                              { 0, 0, 0},
                              { 0, 0, 1},
                              { 0,-1, 0}},
                             {{ 0, 0,-1},
                              { 0, 0, 0},
                              { 1, 0, 0}},
                             {{ 0, 1, 0},
                              {-1, 0, 0},
                              { 0, 0, 0}}};
   INT4 Dim=3,i,j,k;
   INITSTATUS(status, "LALDeterminant3", LALDETERMINANT3C);
   ATTATCHSTATUSPTR(status);

   ASSERT (matrix,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   *determinant = 0;
   for (i=0; i<Dim; i++) {
   for (j=0; j<Dim; j++) {
   for (k=0; k<Dim; k++) {
      *determinant+=epsilon[i][j][k]*matrix[0][i]*matrix[1][j]*matrix[2][k];
   }}}

   ASSERT (*determinant != 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
