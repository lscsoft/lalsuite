/*  <lalVerbatim file="LALInsidePolygonCV">
Author: Cokelaer. T
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInsidePolygon.c}}

Module to check whether a point with coordinate (x,y) in inside 
a polygon defined by the vector (vx, vy). It return valid = 1 if 
the point is inside or valid = 0 if outside. 

The polygon must be defined by its n points (not n+1 ). The module
 uses n+1 points (append the first point to the end of the arrays)
 in order to close the polygon. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInsidePolygonCP}
\idx{LALInsidePolygon()}
\begin{itemize}
   \item \texttt{vx, vy} Input, two arrays of floats defining the polygon. 
   \item \texttt{x, y} Input, coordinate of the point to look at.
   \item \texttt{minimalmatch} Input, the minimal match
\end{itemize}

\subsubsection*{Description/Algorithm}
The algorithm computes the angle between the point(x,y) and each segment
 of the polygon. If the pint is inside the polygon, the sum of the angles 
should be 2 pi. If ouside it is less than 2 pi. 




\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInsidePolygonCV}}

</lalLaTeX>  */



#include <lal/LALInspiralBank.h>

NRCSID (LALINSIDEPOLYGONC, "Id: $");

/*  <lalVerbatim file="LALInsidePolygonCP"> */
void LALInsidePolygon(LALStatus          *status,
                             REAL4              *inputx,
                             REAL4              *inputy,
                             INT4               n,
                             REAL4              x0,
                             REAL4              y0,
                             INT4               *valid)

                             
{  /*  </lalVerbatim>  */


   INITSTATUS (status, "LALInsidePolygon", LALINSIDEPOLYGONC);
   ATTATCHSTATUSPTR(status);
   ASSERT (n>2,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   
   {
     int i, j, c = 0;
     for (i = 0, j = n-1; i < n; j = i++) {
       if ((((inputy[i] <= y0) && (y0 < inputy[j])) ||
	    ((inputy[j] <= y0) && (y0 < inputy[i]))) &&
	   (x0 < (inputx[j] - inputx[i]) * (y0 - inputy[i]) / (inputy[j] - inputy[i]) + inputx[i]))
	 c = !c;
     }
     *valid = c;
   }

   

   
   DETATCHSTATUSPTR(status);
   RETURN(status);

}
