/*  <lalVerbatim file="LALInspiralWaveOverlapCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Module \texttt{LALInspiralWaveOverlap.c}}
Module to compute the overlap of a given data set with
two orthogonal inspiral signals of specified parameters 
with a weight specified in a psd array. The code also returns
in a parameter structure the maximum of the overlap, the bin
where the maximum occured and the phase at the maximum. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWaveOverlapCP}
\index{\verb&LALInspiralWaveOverlap()&}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWaveOverlapCV}}
</lalLaTeX>  */
#include <lal/LALNoiseModels.h>

NRCSID (LALOVERLAPC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveOverlapCP"> */
void
LALInspiralWaveOverlap (
   LALStatus   *status,
   REAL4Vector *output,
   OverlapOut  *overlapout,
   OverlapIn   *overlapin)
{  /*  </lalVerbatim>  */
   REAL4Vector filter1, filter2, output1, output2;
   CorrelateIn corrin;
   REAL8 norm, x, y, z, phase;
   INT4 i;

   INITSTATUS (status, "LALOverlap", LALOVERLAPC);
   ATTATCHSTATUSPTR(status);
   output1.length = output2.length = overlapin->signal.length;
   filter1.length = filter2.length = overlapin->signal.length;
   output1.data = (REAL4*) LALMalloc(sizeof(REAL4)*output1.length);
   output2.data = (REAL4*) LALMalloc(sizeof(REAL4)*output2.length);
   filter1.data = (REAL4*) LALMalloc(sizeof(REAL4)*filter1.length);
   filter2.data = (REAL4*) LALMalloc(sizeof(REAL4)*filter2.length);
   overlapin->param.nStartPad = 0;
   overlapin->param.startPhase = LAL_PI_2;
   LALInspiralWave(status->statusPtr, &output2, &overlapin->param);
   CHECKSTATUSPTR(status);
   overlapin->param.startPhase = 0.;
   LALInspiralWave(status->statusPtr, &output1, &overlapin->param);
   CHECKSTATUSPTR(status);
   LALREAL4VectorFFT(status->statusPtr, &filter1, &output1, overlapin->fwdp);
   CHECKSTATUSPTR(status);
   LALREAL4VectorFFT(status->statusPtr, &filter2, &output2, overlapin->fwdp);
   CHECKSTATUSPTR(status);
   LALNormalise(status->statusPtr, &filter1, &norm, overlapin->psd);
   CHECKSTATUSPTR(status);
   LALNormalise(status->statusPtr, &filter2, &norm, overlapin->psd);
   CHECKSTATUSPTR(status);
   corrin.psd = overlapin->psd;
   corrin.signal1 = overlapin->signal;
   corrin.signal2 = filter1;
   corrin.revp = overlapin->revp;
   LALCorrelate(status->statusPtr, &output1, corrin);
   CHECKSTATUSPTR(status);
   corrin.signal2 = filter2;
   LALCorrelate(status->statusPtr, &output2, corrin);
   CHECKSTATUSPTR(status);

   x = output1.data[0];
   y = output2.data[0];
   overlapout->max = x*x + y*y;
   overlapout->bin = 0;
   overlapout->phase = atan2(y,x);
   for (i=1; i<(int)output1.length; i++) {
       x = output1.data[i];
       y = output2.data[i];
       z = x*x + y*y;
       if (z>overlapout->max) {
          overlapout->max = z;
          overlapout->bin = i;
          overlapout->phase = atan2(y,x);
       }
   }

   phase = overlapout->phase;
   for (i=0; i<(int)output->length; i++) 
      output->data[i] = cos(phase) * output1.data[i] 
                      + sin(phase) * output2.data[i];
   overlapout->max = sqrt(overlapout->max);
   LALFree(filter1.data);
   LALFree(filter2.data);
   LALFree(output1.data);
   LALFree(output2.data);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
