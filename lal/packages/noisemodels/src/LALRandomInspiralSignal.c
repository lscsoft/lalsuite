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
\idx{LALRandomInspiralSignal()}

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
#include <lal/Random.h>

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
   INT4 valid;
   static RandomParams *randomparams;
   
   INITSTATUS (status, "LALRandomInspiralSignal", LALRANDOMINSPIRALSIGNALC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (randIn->psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (randIn->mMin > 0, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (randIn->MMax > 2*randIn->mMin, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (randIn->type >= 0, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (randIn->type <= 2, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

   buff.length = signal->length;
   if (!(buff.data = (REAL4*) LALMalloc(sizeof(REAL4)*buff.length))) {
      ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
   }
   srandom(randIn->useed);
   randIn->useed = random();
   valid = 0;
   while (!valid) 
   {
      e1 = random()/(float)RAND_MAX;
      e2 = random()/(float)RAND_MAX;
      switch (randIn->param.massChoice) 
      {
         case m1Andm2: 
            randIn->param.mass1 = randIn->mMin 
               + (randIn->MMax - 2.*randIn->mMin) * e1;
            randIn->param.mass2 = randIn->mMin 
               + (randIn->MMax - randIn->param.mass1 - randIn->mMin) * e2;
            break;
         case t02: 
            randIn->param.t0 = randIn->t0Min+(randIn->t0Max-randIn->t0Min)*e1;
            randIn->param.t2 = randIn->tnMin+(randIn->tnMax-randIn->tnMin)*e2;
            break;
         case t03: 
            randIn->param.t0 = randIn->t0Min+(randIn->t0Max-randIn->t0Min)*e1;
            randIn->param.t3 = randIn->tnMin+(randIn->tnMax-randIn->tnMin)*e2;
            break;
         default:
            ABORT (status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
            break;
      }
      LALInspiralParameterCalc(status->statusPtr, &(randIn->param));

      if (randIn->param.mass1 > randIn->mMin &&
          randIn->param.mass2 > randIn->mMin &&
          randIn->param.totalMass < randIn->MMax &&
          randIn->param.eta <= 0.25 &&
          randIn->param.eta > randIn->etaMin)
      {
           valid = 1;
      }

   }
   switch (randIn->type) 
   {
      case 0:
         LALInspiralWave(status->statusPtr, &buff, &randIn->param);
         CHECKSTATUSPTR(status);
         LALREAL4VectorFFT(status->statusPtr, signal, &buff, randIn->fwdp);
         CHECKSTATUSPTR(status);
         LALInspiralWaveNormalise(status->statusPtr, signal, &norm, randIn->psd);
         CHECKSTATUSPTR(status);
         break;
      case 1:
/*
         LALGaussianNoise(status->statusPtr, &buff, &randIn->useed);
*/
         LALCreateRandomParams(status->statusPtr, &randomparams, randIn->useed);
         CHECKSTATUSPTR(status);
         LALNormalDeviates(status->statusPtr, &buff, randomparams);
         CHECKSTATUSPTR(status);
         LALDestroyRandomParams(status->statusPtr, &randomparams);
         CHECKSTATUSPTR(status);
         LALREAL4VectorFFT(status->statusPtr, signal, &buff, randIn->fwdp);
         CHECKSTATUSPTR(status);
         LALColoredNoise(status->statusPtr, signal, randIn->psd);
         CHECKSTATUSPTR(status);
         break;
      default:
         noisy.length = signal->length;
         if (!(noisy.data = (REAL4*) LALMalloc(sizeof(REAL4)*noisy.length))) 
         {
            if (buff.data != NULL) LALFree(buff.data);
            buff.data = NULL;
            ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
         }
/*
         LALGaussianNoise(status->statusPtr, &buff, &randIn->useed);
*/
         LALCreateRandomParams(status->statusPtr, &randomparams, randIn->useed);
         CHECKSTATUSPTR(status);
         LALNormalDeviates(status->statusPtr, &buff, randomparams);
         CHECKSTATUSPTR(status);
         LALDestroyRandomParams(status->statusPtr, &randomparams);
         CHECKSTATUSPTR(status);
         LALREAL4VectorFFT(status->statusPtr, &noisy, &buff, randIn->fwdp);
         CHECKSTATUSPTR(status);
         LALColoredNoise(status->statusPtr, &noisy, randIn->psd);
         CHECKSTATUSPTR(status);

         LALInspiralWave(status->statusPtr, signal, &randIn->param);
         CHECKSTATUSPTR(status);
         LALREAL4VectorFFT(status->statusPtr, &buff, signal, randIn->fwdp);
         CHECKSTATUSPTR(status);
         LALInspiralWaveNormalise(status->statusPtr, &buff, &norm, randIn->psd);
         CHECKSTATUSPTR(status);

         addIn.v1 = &buff;
         addIn.a1 = randIn->SignalAmp;
         addIn.v2 = &noisy;
         addIn.a2 = randIn->NoiseAmp;
         LALAddVectors(status->statusPtr, signal, addIn);
         CHECKSTATUSPTR(status);
         if (noisy.data != NULL) LALFree(noisy.data);
         break;
   }
   if (buff.data != NULL) LALFree(buff.data);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

