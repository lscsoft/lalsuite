/*  <lalVerbatim file="LALMatrixTransformCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALMatrixTransform.c}}


\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALMatrixTransformCP}
\index{\texttt{LALMatrixTransform()}}

\subsubsection*{Description}


\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALMatrixTransformCV}}

</lalLaTeX>  */



#include <lal/LALInspiralBank.h>

/*  <lalVerbatim file="LALMatrixTransformCP"> */
void LALMatrixTransform (LALStatus *status, 
                         INT4 n, 
                         REAL8 **data1, 
                         REAL8 **data2, 
                         REAL8 **data3)
{ /* </lalVerbatim> */

   INT4 i, j, l, m;
   status = NULL;

   for (i=0; i<n; i++) {
   for (j=0; j<n; j++) {
      data3[i][j] = 0.0;
      for (l=0; l<n; l++) {
      for (m=0; m<n; m++) {
	 data3[i][j] += data1[i][m]*data2[m][l]*data1[j][l];
      }}
   }}
}
