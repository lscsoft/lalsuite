/*  <lalVerbatim file="LALInsidePolygonCV">
Author: Cokelaer. T
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInsidePolygon.c}}

Module to check whether a point with coordinates (x,y) in inside 
a polygon defined by the vectors (vx, vy). It return valid = 1 if 
the point is inside or valid = 0 if outside. 

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

\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInsidePolygonCV}}

</lalLaTeX>  */



#include <lal/LALInspiralBank.h>

NRCSID (LALINSIDEPOLYGONC, "Id: $");

/*  <lalVerbatim file="LALInsidePolygonCP"> */
void LALInsidePolygon(  LALStatus          *status,
                        REAL4              *inputx,
                        REAL4              *inputy,
                        INT4               n,
                        REAL4              x0,
                        REAL4              y0,
                        INT4               *valid)

                             
{  /*  </lalVerbatim>  */


   INITSTATUS (status, "LALInsidePolygon", LALINSIDEPOLYGONC);
   ATTATCHSTATUSPTR(status);
   ASSERT (n>=3,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   
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
