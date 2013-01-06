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

/**
\author Cokelaer Thomas.
\file

\brief Module to find the vertices of an hexagon inscribed in an ellipse
given its centre, half side-lengths and orientation angle.

\heading{Prototypes}

<tt>LALHexagonVertices()</tt>
<ul>
   <li> <tt>out,</tt> Output.
   </li><li> <tt>in,</tt> Input.</li>
</ul>

\heading{Description}

This code computes the vertices of an hexagon for plotting
a grid of templates with xmgr, useful when looking at the
minimal-match-Hexagons around mesh points in a template bank.
Used by SpaceCovering in the test directory.

\heading{Algorithm}
Given the centre \f$(x_0,y_0)\f$ and half-sides \f$(dx,dy),\f$
the vertices of a Hexagon in a \e diagonal coordinate
system are given by
\f{eqnarray}{
x_1 & = & x_0 - dx, \quad y_1 = y_0 - dy, \nonumber \\
x_2 & = & x_0 + dx, \quad y_2 = y_0 - dy, \nonumber \\
x_3 & = & x_0 + dx, \quad y_3 = y_0 + dy, \nonumber \\
x_4 & = & x_0 - dx, \quad y_4 = y_0 + dy. \nonumber
\f}
The coordinates of a Hexagon oriented at an angle \f$\theta\f$ is
found by using the formulas
\f{eqnarray}{
x' = x \cos(\theta) - y \sin(\theta),\nonumber \\
y' = y \cos(\theta) + x \sin(\theta).\nonumber
\f}
The function returns 7 coordinate points (1,2,3,4,5,6,1),
and not just the 6 verticies, to help a plotting programme
to complete the Hexagon.

\heading{Uses}
None.

\heading{Notes}

*/

#include <lal/LALInspiralBank.h>

void
LALHexagonVertices(
   LALStatus *status,
   HexagonOut *out,
   RectangleIn *in
)
{

   REAL4 x_1, x_2, x_3, x_4, x_5, x_6;
   REAL4 y_1, y_2, y_3, y_4, y_5, y_6;
   REAL4 ctheta,stheta, sca;
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   ASSERT (out,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (in,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   sca = sqrt(3);

   x_1 = -in->dx/2;
   y_1 = -in->dy/sca/2;
   x_2 = 0;
   y_2 = -in->dy/sqrt(3);
   x_3 = in->dx/2;
   y_3 = -in->dy/sca/2;
   x_4 = in->dx/2;
   y_4 = in->dy/sca/2;
   x_5 = 0;
   y_5 = in->dy/sqrt(3);
   x_6 = -in->dx/2;
   y_6 = in->dy/sca/2;

   ctheta=cos(in->theta);
   stheta=sin(in->theta);

   out->x1 = in->x0 + x_1 * ctheta - y_1 * stheta;
   out->y1 = in->y0 + y_1 * ctheta + x_1 * stheta;
   out->x2 = in->x0 + x_2 * ctheta - y_2 * stheta;
   out->y2 = in->y0 + y_2 * ctheta + x_2 * stheta;
   out->x3 = in->x0 + x_3 * ctheta - y_3 * stheta;
   out->y3 = in->y0 + y_3 * ctheta + x_3 * stheta;
   out->x4 = in->x0 + x_4 * ctheta - y_4 * stheta;
   out->y4 = in->y0 + y_4 * ctheta + x_4 * stheta;
   out->x5 = in->x0 + x_5 * ctheta - y_5 * stheta;
   out->y5 = in->y0 + y_5 * ctheta + x_5 * stheta;
   out->x6 = in->x0 + x_6 * ctheta - y_6 * stheta;
   out->y6 = in->y0 + y_6 * ctheta + x_6 * stheta;

   out->x7 = out->x1;
   out->y7 = out->y1;
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
