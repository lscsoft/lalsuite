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
\idx{LALInspiralWaveOverlap()}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
LALInspiralWave
LALREAL4VectorFFT
LALInspiralWaveNormaliseLSO
LALInspiralWaveCorrelate
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWaveOverlapCV}}
</lalLaTeX>  */
#include <lal/LALNoiseModels.h>

NRCSID (LALINSPIRALWAVEOVERLAPC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveOverlapCP"> */
void
LALInspiralWaveOverlap 
   (
   LALStatus               *status,
   REAL4Vector             *output,
   InspiralWaveOverlapOut  *overlapout,
   InspiralWaveOverlapIn   *overlapin
   )
{  /*  </lalVerbatim>  */
   REAL4Vector filter1, filter2, output1, output2;
   InspiralWaveCorrelateIn corrin;
   REAL8 norm, x, y, z, phase;
   INT4 i, nBegin, nEnd;
   InspiralWaveNormaliseIn normin;

   INITSTATUS (status, "LALInspiralWaveOverlap", LALINSPIRALWAVEOVERLAPC);
   ATTATCHSTATUSPTR(status);

   ASSERT (output->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (overlapin->psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (overlapin->signal.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   output1.length = output2.length = overlapin->signal.length;
   filter1.length = filter2.length = overlapin->signal.length;

   if (!(output1.data = (REAL4*) LALMalloc(sizeof(REAL4)*output1.length))) {
      ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
   }
   if (!(output2.data = (REAL4*) LALMalloc(sizeof(REAL4)*output2.length))) {
      LALFree(output1.data);
      output1.data = NULL;
      ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
   }
   if (!(filter1.data = (REAL4*) LALMalloc(sizeof(REAL4)*filter1.length))) {
      LALFree(output1.data);
      LALFree(output2.data);
      output1.data = NULL;
      output2.data = NULL;
      ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
   }
   if (!(filter2.data = (REAL4*) LALMalloc(sizeof(REAL4)*filter2.length))) {
      LALFree(output1.data);
      LALFree(output2.data);
      LALFree(filter1.data);
      output1.data = NULL;
      output2.data = NULL;
      filter1.data = NULL;
      ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
   }

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
   normin.psd = &(overlapin->psd);
   normin.df = overlapin->param.tSampling / (REAL8) overlapin->signal.length;
   normin.fCutoff = overlapin->param.fCutoff;
   normin.samplingRate = overlapin->param.tSampling;
   LALInspiralWaveNormaliseLSO(status->statusPtr, &filter1, &norm, &normin);
   CHECKSTATUSPTR(status);
   LALInspiralWaveNormaliseLSO(status->statusPtr, &filter2, &norm, &normin);
   CHECKSTATUSPTR(status);
   corrin.fCutoff = overlapin->param.fCutoff;
   corrin.samplingRate = overlapin->param.tSampling;
   corrin.df = overlapin->param.tSampling / (REAL8) overlapin->signal.length;
   corrin.psd = overlapin->psd;
   corrin.signal1 = overlapin->signal;
   corrin.signal2 = filter1;
   corrin.revp = overlapin->revp;
   LALInspiralWaveCorrelate(status->statusPtr, &output1, corrin);
   CHECKSTATUSPTR(status);
   corrin.signal2 = filter2;
   LALInspiralWaveCorrelate(status->statusPtr, &output2, corrin);
   CHECKSTATUSPTR(status);

   nBegin = overlapin->nBegin;
   nEnd = output1.length - overlapin->nEnd;
   x = output1.data[nBegin];
   y = output2.data[nBegin];
   overlapout->max = x*x + y*y;
   overlapout->bin = nBegin;
   overlapout->phase = atan2(y,x);
   for (i=nBegin; i<nEnd; i++) {
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
   if (filter1.data!=NULL) LALFree(filter1.data);
   if (filter2.data!=NULL) LALFree(filter2.data);
   if (output1.data!=NULL) LALFree(output1.data);
   if (output2.data!=NULL) LALFree(output2.data);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
