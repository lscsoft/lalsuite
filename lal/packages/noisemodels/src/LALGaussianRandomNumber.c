/*  <lalVerbatim file="LALGaussianRandomNumberCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Module \texttt{LALGaussianRandomNumber.c}}
Module to create a Gaussian random number.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALGaussianRandomNumberCP}
\index{\verb&LALGaussianRandomNumber()&}

\subsubsection*{Description}
Using the c internal random number generates this routine
transforms a uniform random number into a Gaussian random number.
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALGaussianRandomNumberCV}}
</lalLaTeX>  */
#include <stdlib.h>
#include <lal/LALNoiseModels.h>

#define random() rand()
#define srandom( seed ) srand( seed )

NRCSID (LALGAUSSIANRANDOMNUMBERC, "$Id$");
/*  <lalVerbatim file="LALGaussianRandomNumberCP"> */
void
LALGaussianRandomNumber(
   LALStatus *status, 
   REAL4 *randnum)
{  /*  </lalVerbatim>  */
   static int iset=0;
   static REAL8 gset;
   REAL8 fac,rsq,v1,v2;

   INITSTATUS (status, "LALGaussianRandomNumber", LALGAUSSIANRANDOMNUMBERC);
   ATTATCHSTATUSPTR(status);

   if  (iset == 0) {
      do {
         v1 = 2.*(random()/(REAL8)RAND_MAX) - 1.0;
         v2 = 2.*(random()/(REAL8)RAND_MAX) - 1.0;
         rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      *randnum = (REAL4) (v2*fac);
   } else {
      iset=0;
      *randnum = (REAL4) (gset);
   }

   DETATCHSTATUSPTR(status);
   RETURN(status); 
}
