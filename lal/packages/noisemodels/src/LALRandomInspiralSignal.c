/*  <lalVerbatim file="LALRandomInspiralSignalCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Module \texttt{LALRandomInspiralSignal.c}}
Module to generate signals with random masses in the parameter space
and to add a noisy component expected in a given detector. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALRandomInspiralSignalCP}
\index{\verb&LALRandomInspiralSignal()&}

\subsubsection*{Description}
Depending on the value of the parameter \texttt{RandomIn.type=0,1 or 2} this
code returns a pure signal, a pure noise or signal+noise.
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALRandomInspiralSignalCV}}
</lalLaTeX>  */
#include <stdlib.h>
#include <lal/LALNoiseModels.h>

#define random() rand()
#define srandom( seed ) srand( seed )

NRCSID (LALRANDOMINSPIRALSIGNALC, "$Id$");
/*  <lalVerbatim file="LALRandomInspiralSignalCP"> */
void
LALRandomInspiralSignal(
   LALStatus *status, 
   REAL4Vector *signal,
   RandomInspiralSignalIn *randIn)
{  /*  </lalVerbatim>  */
   REAL8 e1, e2, norm;
   REAL4Vector noisy, buff;
   AddVectorsIn addIn;
   
   INITSTATUS (status, "LALRandomInspiralSignal", LALRANDOMINSPIRALSIGNALC);
   ATTATCHSTATUSPTR(status);
   buff.length = signal->length;
   buff.data = (REAL4*) LALMalloc(sizeof(REAL4)*buff.length);
   srandom(randIn->useed);
   randIn->useed = random();
   e1 = random()/(float)RAND_MAX;
   e2 = random()/(float)RAND_MAX;
   randIn->param.mass1 = randIn->mMin + (randIn->MMax - 2.*randIn->mMin) * e1;
   randIn->param.mass2 = randIn->mMin + (randIn->MMax - randIn->param.mass1 - randIn->mMin) * e2;
   switch (randIn->type) {
      case 0:
         LALInspiralWave(status->statusPtr, &buff, &randIn->param);
         CHECKSTATUSPTR(status);
         LALREAL4VectorFFT(status->statusPtr, signal, &buff, randIn->fwdp);
         CHECKSTATUSPTR(status);
         LALNormalise(status->statusPtr, signal, &norm, randIn->psd);
         CHECKSTATUSPTR(status);
         break;
      case 1:
         LALGaussianNoise(status->statusPtr, &buff, &randIn->useed);
         CHECKSTATUSPTR(status);
         LALREAL4VectorFFT(status->statusPtr, signal, &buff, randIn->fwdp);
         CHECKSTATUSPTR(status);
         LALColoredNoise(status->statusPtr, signal, randIn->psd);
         CHECKSTATUSPTR(status);
         break;
      default:
         noisy.length = signal->length;
         noisy.data = (REAL4*) LALMalloc(sizeof(REAL4)*noisy.length);
         LALGaussianNoise(status->statusPtr, &buff, &randIn->useed);
         CHECKSTATUSPTR(status);
         LALREAL4VectorFFT(status->statusPtr, &noisy, &buff, randIn->fwdp);
         CHECKSTATUSPTR(status);
         LALColoredNoise(status->statusPtr, &noisy, randIn->psd);
         CHECKSTATUSPTR(status);

         LALInspiralWave(status->statusPtr, signal, &randIn->param);
         CHECKSTATUSPTR(status);
         LALREAL4VectorFFT(status->statusPtr, &buff, signal, randIn->fwdp);
         CHECKSTATUSPTR(status);
         LALNormalise(status->statusPtr, &buff, &norm, randIn->psd);
         CHECKSTATUSPTR(status);

         addIn.v1 = &buff;
         addIn.a1 = randIn->SignalAmp;
         addIn.v2 = &noisy;
         addIn.a2 = randIn->NoiseAmp;
         LALAddVectors(status->statusPtr, signal, addIn);
         CHECKSTATUSPTR(status);
         LALFree(noisy.data);
         break;
   }
   LALFree(buff.data);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

