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
\index{\texttt{LALInverse3()}}

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

\vfill{\footnotesize\input{LALInverse3CV}}

</lalLaTeX>  */





#include <lal/LALInspiralBank.h>

NRCSID(LALINVERSE3C, "$Id$");
#define Dim 3

/*  <lalVerbatim file="LALInverse3CP"> */
void LALInverse3(LALStatus *status, 
                 REAL8 **inverse, 
                 REAL8 **matrix) 
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

   LALDeterminant3(status->statusPtr, &det, matrix);
   CHECKSTATUSPTR(status);
   if(det==0) {
      fprintf(stderr, "LALInverse3.c: Singular matrix\n");
      exit(0);
   }
   for (i=0; i<Dim; i++) {
   for (j=0; j<Dim; j++) {
      x = 0;
      for (a=0; a<Dim; a++) { for (b=0; b<Dim; b++) {
      for (p=0; p<Dim; p++) { for (q=0; q<Dim; q++) { 
         x+=epsilon[j][a][p] * epsilon[i][b][q] * matrix[a][b] * matrix[p][q];
      }}}}
      inverse[i][j] = x/(2.*det);
/*
      fprintf(stderr, "In LALInverse3.c %d %d %e %e\n", i, j, matrix[i][j], inverse[i][j]);
*/
   }}
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

#undef Dim
