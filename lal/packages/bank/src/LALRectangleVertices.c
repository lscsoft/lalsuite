/*  <lalVerbatim file="LALRectangleVerticesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALRectangleVertices.c}}

Module to find the vertices of a rectangle given its centre,
half side-lengths and orientation angle.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALRectangleVerticesCP}
\index{\verb&LALRectangleVertices()&}

\subsubsection*{Description}

This code computes the vertices of a rectangle for plotting
a grid of templates with xmgr.

\subsubsection*{Algorithm}
Given the centre $(x_0,y_0)$ and half-sides $(dx,dy),$ 
the vertices of a rectangle in a {\it diangle} coordinate 
system are given by
\begin{eqnarray}
x_1 & = & x_0 - dx, \ \ y_1 = y_0 - dy, \nonumber \\
x_2 & = & x_0 + dx, \ \ y_2 = y_0 - dy, \nonumber \\
x_3 & = & x_0 + dx, \ \ y_3 = y_0 + dy, \nonumber \\
x_4 & = & x_0 - dx, \ \ y_4 = y_0 + dy. \nonumber
\end{eqnarray}
The coordinates of a rectangle oriented at an angle $\theta$ is
found by using the formulas
\begin{eqnarray}
x' = x \cos(\theta) - y \sin(\theta),\nonumber \\
y' = y \cos(\theta) + x \sin(\theta).\nonumber 
\end{eqnarray}

\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALRectangleVerticesCV}}

</lalLaTeX>  */

#include <lal/LALInspiralBank.h>

/*
typedef struct
tagRectangleIn 
   {REAL8 x0, y0, dx, dy, theta;}
RectangleIn;

typedef struct
tagRectangleOut 
   {REAL8 x1, y1, x2, y2, x3, y3, x4, y4, x5, y5;}
RectangleOut;

void 
LALRectangleVertices(
   LALStatus *status, 
   RectangleOut *out,
   RectangleIn *in
);
*/

NRCSID(LALRECTANGLEVERTICESC, "$Id$");

/*  <lalVerbatim file="LALRectangleVerticesCP"> */

void 
LALRectangleVertices(
   LALStatus *status, 
   RectangleOut *out,
   RectangleIn *in
) 
{ /* </lalVerbatim> */

   REAL8 x1, x2, x3, x4, y1, y2, y3, y4;
   INITSTATUS(status, "LALRectangleVertices", LALRECTANGLEVERTICESC);
   ATTATCHSTATUSPTR(status);

   ASSERT (out,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (in,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

/*
   x1 = out->x1 = in->x0-in->dx/2.;
   y1 = out->y1 = in->y0-in->dy/2.;
   x2 = out->x2 = in->x0+in->dx/2.;
   y2 = out->y2 = in->y0-in->dy/2.;
   x3 = out->x3 = in->x0+in->dx/2.;
   y3 = out->y3 = in->y0+in->dy/2.;
   x4 = out->x4 = in->x0-in->dx/2.;
   y4 = out->y4 = in->y0+in->dy/2.;

   out->x1 = x1 * cos(in->theta) - y1 * sin(in->theta);
   out->y1 = y1 * cos(in->theta) + x1 * sin(in->theta);
   out->x2 = x2 * cos(in->theta) - y2 * sin(in->theta);
   out->y2 = y2 * cos(in->theta) + x2 * sin(in->theta);
   out->x3 = x3 * cos(in->theta) - y3 * sin(in->theta);
   out->y3 = y3 * cos(in->theta) + x3 * sin(in->theta);
   out->x4 = x4 * cos(in->theta) - y4 * sin(in->theta);
   out->y4 = y4 * cos(in->theta) + x4 * sin(in->theta);
*/
   x1 = -in->dx/2.;
   y1 = -in->dy/2.;
   x2 =  in->dx/2.;
   y2 = -in->dy/2.;
   x3 =  in->dx/2.;
   y3 =  in->dy/2.;
   x4 = -in->dx/2.;
   y4 =  in->dy/2.;

   out->x1 = in->x0 + x1 * cos(in->theta) - y1 * sin(in->theta);
   out->y1 = in->y0 + y1 * cos(in->theta) + x1 * sin(in->theta);
   out->x2 = in->x0 + x2 * cos(in->theta) - y2 * sin(in->theta);
   out->y2 = in->y0 + y2 * cos(in->theta) + x2 * sin(in->theta);
   out->x3 = in->x0 + x3 * cos(in->theta) - y3 * sin(in->theta);
   out->y3 = in->y0 + y3 * cos(in->theta) + x3 * sin(in->theta);
   out->x4 = in->x0 + x4 * cos(in->theta) - y4 * sin(in->theta);
   out->y4 = in->y0 + y4 * cos(in->theta) + x4 * sin(in->theta);

   out->x5 = out->x1;
   out->y5 = out->y1;

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
