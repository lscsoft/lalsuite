/*  <lalVerbatim file="LALCorrelateCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/* <lalLaTeX>
\subsection{Module \texttt{LALCorrelate.c}}
Module to compute the correlation of two data sets.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALCorrelateCP}
\index{\verb&LALCorrelate()&}

\subsubsection*{Description}
The module expects two inputs \texttt{signal1, signal2}
in the Fourier-domain, computes the correlation
and returns the correlated output in the time-domain weighted by the
noise \texttt{psd}.
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALCorrelateCV}}
</lalLaTeX>  */

#include <lal/LALNoiseModels.h>
#include <lal/RealFFT.h>

NRCSID (LALCORRELATEC, "$Id$");

/*  <lalVerbatim file="LALCorrelateCP"> */
void
LALCorrelate (
   LALStatus *status,
   REAL4Vector *output,
   CorrelateIn corrin)
{  /*  </lalVerbatim>  */
  INT4 n, nby2, i, k;
  REAL8 psd;
  REAL4Vector buff;


  INITSTATUS (status, "LALCorrelate", LALCORRELATEC);
  ATTATCHSTATUSPTR(status);

  if (corrin.psd.length != corrin.signal1.length/2+1) {
      fprintf(stderr, "LALCorrelate: Incompatible lengths\n");
      exit(1);
  }
  if (corrin.signal1.length != corrin.signal2.length) {
      fprintf(stderr, "LALCorrelate: Incompatible array lengths\n");
      exit(1);
  }

  n = corrin.signal1.length;

  buff.length = corrin.signal1.length;
  buff.data = (REAL4*) LALMalloc(sizeof(REAL4)*buff.length);

  nby2 = n/2;
  for (i=1; i<nby2; i++) {
    k=n-i;
    psd = corrin.psd.data[i+1];
    if (psd) {
       output->data[i] = (corrin.signal1.data[i]*corrin.signal2.data[i] 
                    +  corrin.signal1.data[k]*corrin.signal2.data[k]) / psd;
       output->data[k] = (corrin.signal1.data[k]*corrin.signal2.data[i] 
                    -  corrin.signal1.data[i]*corrin.signal2.data[k]) / psd;
    } else {
       output->data[i] = output->data[k] = 0;
    }
  }
  psd = corrin.psd.data[0];
  if (psd) 
     output->data[0] = corrin.signal1.data[0]*corrin.signal2.data[0] / psd;
  else
     output->data[0] = 0;
  psd = corrin.psd.data[nby2];
  if (psd) 
     output->data[nby2] = corrin.signal1.data[nby2]*corrin.signal2.data[nby2] / psd;
  else
     output->data[nby2] = 0;
  LALREAL4VectorFFT(status->statusPtr,&buff,output,corrin.revp);
  CHECKSTATUSPTR(status);
  for (i=0; i< (int) buff.length; i++) output->data[i] = buff.data[i]/2.;
  LALFree(buff.data);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
