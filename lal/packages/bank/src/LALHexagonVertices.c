/*
*  Copyright (C) 2007 Thomas Cokelaer
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

/*  <lalVerbatim file="LALHexagonVerticesCV">
Author: Cokelaer Thomas.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALHexagonVertices.c}}

Module to find the vertices of an hexagon inscribed in an ellipse
given its centre, half side-lengths and orientation angle.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALRectangleVerticesCP}
\index{\verb&LALHexagonVertices()&}
\begin{itemize}
   \item \texttt{out,} Output.
   \item \texttt{in,} Input.
\end{itemize}

\subsubsection*{Description}

This code computes the vertices of an hexagon for plotting
a grid of templates with xmgr, useful when looking at the
minimal-match-Hexagons around mesh points in a template bank.
Used by SpaceCovering in the test directory.

\subsubsection*{Algorithm}
Given the centre $(x_0,y_0)$ and half-sides $(dx,dy),$
the vertices of a Hexagon in a {\it diagonal} coordinate
system are given by
\begin{eqnarray}
x_1 & = & x_0 - dx, \ \ y_1 = y_0 - dy, \nonumber \\
x_2 & = & x_0 + dx, \ \ y_2 = y_0 - dy, \nonumber \\
x_3 & = & x_0 + dx, \ \ y_3 = y_0 + dy, \nonumber \\
x_4 & = & x_0 - dx, \ \ y_4 = y_0 + dy. \nonumber
\end{eqnarray}
The coordinates of a Hexagon oriented at an angle $\theta$ is
found by using the formulas
\begin{eqnarray}
x' = x \cos(\theta) - y \sin(\theta),\nonumber \\
y' = y \cos(\theta) + x \sin(\theta).\nonumber
\end{eqnarray}
The function returns 7 coordinate points (1,2,3,4,5,6,1),
and not just the 6 verticies, to help a plotting programme
to complete the Hexagon.

\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALHexagonVerticesCV}}

</lalLaTeX>  */

#include <lal/LALInspiralBank.h>

NRCSID(LALHEXAGONVERTICESC, "$Id$");

/*  <lalVerbatim file="LALHexagonVerticesCP"> */

void
LALHexagonVertices(
   LALStatus *status,
   HexagonOut *out,
   RectangleIn *in
)
{ /* </lalVerbatim> */

   REAL4 x1, x2, x3, x4, y1, y2, y3, y4, x5, y5, x6, y6, x7, y7;
   REAL4 ctheta,stheta, sca;
   INITSTATUS(status, "LALHexagonVertices", LALHEXAGONVERTICESC);
   ATTATCHSTATUSPTR(status);

   ASSERT (out,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (in,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   sca = sqrt(3);

   x1 = -in->dx/2;
   y1 = -in->dy/sca/2;
   x2 = 0;
   y2 = -in->dy/sqrt(3);
   x3 = in->dx/2;
   y3 = -in->dy/sca/2;
   x4 = in->dx/2;
   y4 = in->dy/sca/2;
   x5 = 0;
   y5 = in->dy/sqrt(3);
   x6 = -in->dx/2;
   y6 = in->dy/sca/2;

   ctheta=cos(in->theta);
   stheta=sin(in->theta);

   out->x1 = in->x0 + x1 * ctheta - y1 * stheta;
   out->y1 = in->y0 + y1 * ctheta + x1 * stheta;
   out->x2 = in->x0 + x2 * ctheta - y2 * stheta;
   out->y2 = in->y0 + y2 * ctheta + x2 * stheta;
   out->x3 = in->x0 + x3 * ctheta - y3 * stheta;
   out->y3 = in->y0 + y3 * ctheta + x3 * stheta;
   out->x4 = in->x0 + x4 * ctheta - y4 * stheta;
   out->y4 = in->y0 + y4 * ctheta + x4 * stheta;
   out->x5 = in->x0 + x5 * ctheta - y5 * stheta;
   out->y5 = in->y0 + y5 * ctheta + x5 * stheta;
   out->x6 = in->x0 + x6 * ctheta - y6 * stheta;
   out->y6 = in->y0 + y6 * ctheta + x6 * stheta;

   out->x7 = out->x1;
   out->y7 = out->y1;
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
