/*  <lalVerbatim file="LALGaussianNoiseCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/* <lalLaTeX>
\subsection{Module \texttt{LALGaussianNoise.c}}
Module to create GaussianNoise.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALGaussianNoiseCP}
\index{\verb&LALGaussianNoise()&}

\subsubsection*{Description}
Using a random number seed (unsigned integer) this module returns a
Gaussian random number of lenght defined by the \texttt{noisy} structure.
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
LALGaussianRandomNumber()
random()
srandom()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALGaussianNoiseCV}}
</lalLaTeX>  */

#include <stdlib.h>
#include <lal/LALNoiseModels.h>

#define random() rand()
#define srandom( seed ) srand( seed )

NRCSID (LALGAUSSIANNOISEC, "$Id$");
/*  <lalVerbatim file="LALGaussianNoiseCP"> */
void 
LALGaussianNoise (
   LALStatus   *status,
   REAL4Vector *noisy, 
   UINT8       *seed) 
{  /*  </lalVerbatim>  */
   int i;

   INITSTATUS (status, "LALGaussianNoise", LALGAUSSIANNOISEC);
   ATTATCHSTATUSPTR(status);
   ASSERT (noisy,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (noisy->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (*seed > 0, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

   srandom(*seed);
   *seed = random();
   for (i=0; i<(int)noisy->length; i++){
      LALGaussianRandomNumber(status->statusPtr,noisy->data+i);
      CHECKSTATUSPTR(status);
   }
   DETATCHSTATUSPTR(status);
   RETURN(status); 
}

