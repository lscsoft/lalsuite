/*
*  Copyright (C) 2007 David Churches, Jolien Creighton, Thomas Cokelaer
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

#include <lal/LALInspiralBank.h>

/**
\ingroup LALInspiralBank_h
\brief Function to find the vertices of a rectangle given its centre, half side-lengths and orientation angle.
\author Sathyaprakash, B. S.

This code computes the vertices of a rectangle for plotting
a grid of templates with xmgr, useful when looking at the
minimal-match-rectangles around mesh points in a template bank.

\heading{Algorithm}
Given the centre \f$(x_0,y_0)\f$ and half-sides \f$(dx,dy),\f$
the vertices of a rectangle in a \e diagonal coordinate
system are given by
\f{eqnarray}{
x_1 & = & x_0 - dx, \quad y_1 = y_0 - dy, \nonumber \\
x_2 & = & x_0 + dx, \quad y_2 = y_0 - dy, \nonumber \\
x_3 & = & x_0 + dx, \quad y_3 = y_0 + dy, \nonumber \\
x_4 & = & x_0 - dx, \quad y_4 = y_0 + dy. \nonumber
\f}
The coordinates of a rectangle oriented at an angle \f$\theta\f$ is
found by using the formulas
\f{eqnarray}{
x' = x \cos(\theta) - y \sin(\theta),\nonumber \\
y' = y \cos(\theta) + x \sin(\theta).\nonumber
\f}
The function returns five coordinate points (1,2,3,4,1),
and not just the four verticies, to help
a plotting programme to complete the rectangle.
 */
void
LALRectangleVertices(
   LALStatus *status,	/**< LAL status pointer */
   RectangleOut *out,	/**< Output */
   RectangleIn *in	/**< Input */
)
{

   REAL4 x1, x2, x3, x4, myy1, y2, y3, y4;
   REAL4 ctheta,stheta;
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   ASSERT (out,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (in,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

/*
   x1 = out->x1 = in->x0-in->dx/2.;
   myy1 = out->y1 = in->y0-in->dy/2.;
   x2 = out->x2 = in->x0+in->dx/2.;
   y2 = out->y2 = in->y0-in->dy/2.;
   x3 = out->x3 = in->x0+in->dx/2.;
   y3 = out->y3 = in->y0+in->dy/2.;
   x4 = out->x4 = in->x0-in->dx/2.;
   y4 = out->y4 = in->y0+in->dy/2.;

   out->x1 = x1 * cos(in->theta) - myy1 * sin(in->theta);
   out->y1 = myy1 * cos(in->theta) + x1 * sin(in->theta);
   out->x2 = x2 * cos(in->theta) - y2 * sin(in->theta);
   out->y2 = y2 * cos(in->theta) + x2 * sin(in->theta);
   out->x3 = x3 * cos(in->theta) - y3 * sin(in->theta);
   out->y3 = y3 * cos(in->theta) + x3 * sin(in->theta);
   out->x4 = x4 * cos(in->theta) - y4 * sin(in->theta);
   out->y4 = y4 * cos(in->theta) + x4 * sin(in->theta);
*/
   x1 = -in->dx/2.;
   myy1 = -in->dy/2.;
   x2 =  in->dx/2.;
   y2 = -in->dy/2.;
   x3 =  in->dx/2.;
   y3 =  in->dy/2.;
   x4 = -in->dx/2.;
   y4 =  in->dy/2.;
   ctheta=cos(in->theta);
   stheta=sin(in->theta);

   out->x1 = in->x0 + x1 * ctheta - myy1 * stheta;
   out->y1 = in->y0 + myy1 * ctheta + x1 * stheta;
   out->x2 = in->x0 + x2 * ctheta - y2 * stheta;
   out->y2 = in->y0 + y2 * ctheta + x2 * stheta;
   out->x3 = in->x0 + x3 * ctheta - y3 * stheta;
   out->y3 = in->y0 + y3 * ctheta + x3 * stheta;
   out->x4 = in->x0 + x4 * ctheta - y4 * stheta;
   out->y4 = in->y0 + y4 * ctheta + x4 * stheta;

   out->x5 = out->x1;
   out->y5 = out->y1;

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
