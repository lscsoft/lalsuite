/*  <lalVerbatim file="LALInspiralFindLoudestEventCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Module \texttt{LALInspiralFindLoudestEvent.c}}
Module to find events in a given data set with an SNR 
larger than a pre-specified threshold. The module uses
two orthogonal inspiral signals of specified parameters 
with a weight specified in a psd array. The code returns
the number of events found, and for each event the snr, 
the bin number and the phase of the template at that bin.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralFindLoudestEventCP}
\index{\verb&LALInspiralFindLoudestEvent()&}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
LALInspiralWave
LALREAL4VectorFFT
LALInspiralWaveNormalise
LALInspiralWaveCorrelate
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralFindLoudestEventCV}}
</lalLaTeX>  */
#include <lal/LALNoiseModels.h>

NRCSID (LALINSPIRALFINDEVENTSC, "$Id$");

/*  <lalVerbatim file="LALInspiralFindLoudestEventCP"> */
void
LALInspiralFindLoudestEvent 
   (
   LALStatus            *status,
   INT4                 *nEvents,
   InspiralEventsList   *eventlist,
   InspiralFindEventsIn *findeventsin
   )
{  /*  </lalVerbatim>  */

   /* 
    * We shall assume that the number of events to be found in 
    * each template is no more than 100000.
   */
   INT4 i;
   INT4 nBegin;
   INT4 nEnd;

   REAL8 dt;
   REAL8 norm;
   REAL8 x;
   REAL8 y;
   REAL8 z;

   REAL4Vector filter1;
   REAL4Vector filter2;
   REAL4Vector output1;
   REAL4Vector output2;
   REAL4Vector buffer;

   InspiralWaveCorrelateIn corrin;

   INITSTATUS (status, "LALInspiralFindLoudestEvent", LALINSPIRALFINDEVENTSC);
   ATTATCHSTATUSPTR(status);

   ASSERT (findeventsin->psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (findeventsin->signal.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   output1.length = output2.length = findeventsin->signal.length;
   filter1.length = filter2.length = findeventsin->signal.length;

   ASSERT (findeventsin->nEnd >= 0,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
   ASSERT (findeventsin->nEnd <= (INT4)output1.length,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
   ASSERT (findeventsin->nBegin >= 0,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
   ASSERT (findeventsin->nBegin <= (INT4)output1.length,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);


   dt = 1./findeventsin->param.tSampling;
   nBegin = findeventsin->nBegin;
   nEnd = findeventsin->signal.length - findeventsin->nEnd;

   buffer.length = nEnd - nBegin;
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
   if (!(buffer.data = (REAL4*) LALMalloc(sizeof(REAL4)*buffer.length))) {
      LALFree(output1.data);
      LALFree(output2.data);
      LALFree(filter1.data);
      LALFree(filter2.data);
      output1.data = NULL;
      output2.data = NULL;
      filter1.data = NULL;
      filter2.data = NULL;
      ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
   }

   findeventsin->param.nStartPad = 0;
   findeventsin->param.startPhase = LAL_PI_2;
   LALInspiralWave(status->statusPtr, &output2, &findeventsin->param);
   CHECKSTATUSPTR(status);
   findeventsin->param.startPhase = 0.;
   LALInspiralWave(status->statusPtr, &output1, &findeventsin->param);
   CHECKSTATUSPTR(status);
   if (findeventsin->displayTemplates)
   {
      for (i=0;i<output1.length;i++) 
         printf("%e %e\n", i*dt, output1.data[i]);printf("&\n");
      for (i=0;i<output1.length;i++) 
         printf("%e %e\n", i*dt, output2.data[i]);printf("&\n");
   }
   LALREAL4VectorFFT(status->statusPtr, &filter1, &output1, findeventsin->fwdp);
   CHECKSTATUSPTR(status);
   LALREAL4VectorFFT(status->statusPtr, &filter2, &output2, findeventsin->fwdp);
   CHECKSTATUSPTR(status);
   LALInspiralWaveNormalise(status->statusPtr, &filter1, &norm, findeventsin->psd);
   CHECKSTATUSPTR(status);
   LALInspiralWaveNormalise(status->statusPtr, &filter2, &norm, findeventsin->psd);
   CHECKSTATUSPTR(status);
   corrin.psd = findeventsin->psd;
   corrin.signal1 = findeventsin->signal;
   corrin.signal2 = filter1;
   corrin.revp = findeventsin->revp;
   LALInspiralWaveCorrelate(status->statusPtr, &output1, corrin);
   CHECKSTATUSPTR(status);
   corrin.signal2 = filter2;
   LALInspiralWaveCorrelate(status->statusPtr, &output2, corrin);
   CHECKSTATUSPTR(status);

   if (findeventsin->displayCorrelationStats)
   {
	   StatsREAL4VectorOut statsout1;
	   StatsREAL4VectorOut statsout2;
   
	   for (i=nBegin;i<nEnd;i++) 
		   buffer.data[i-nBegin] = output1.data[i];
	   LALStatsREAL4Vector(status->statusPtr, &statsout1, &buffer);
	   CHECKSTATUSPTR(status);
	   fprintf(stderr, "mean=%e std=%e max=%e\n", statsout1.mean, statsout1.stddev, statsout1.max);   
   
	   for (i=nBegin;i<nEnd;i++) 
		   buffer.data[i-nBegin] = output2.data[i];
	   LALStatsREAL4Vector(status->statusPtr, &statsout2, &buffer);
	   CHECKSTATUSPTR(status);
	   fprintf(stderr, "mean=%e std=%e max=%e\n", statsout2.mean, statsout2.stddev, statsout2.max);   
   }
   
   if (findeventsin->displayCorrelation)
   {
      for (i=nBegin;i<nEnd;i++) 
         {
            x = pow ( pow( output1.data[i], 2.) + pow( output2.data[i], 2.), 0.5); 
            printf("%e %e\n",i*dt, x);
         }
         printf("&\n");   
   }

   *nEvents = 0;
   
   x = output1.data[nBegin];
   y = output2.data[nBegin];
   z = sqrt(x*x + y*y);
   eventlist->max = z;
   eventlist->bin = nBegin;
   eventlist->phase = atan2(y,x);

   if (z>findeventsin->Threshold) (*nEvents)++;

   for (i=nBegin; i<nEnd; i++) 
   {
       x = output1.data[i];
       y = output2.data[i];
       z = sqrt(x*x + y*y);
       if (z>eventlist->max) 
       { 
          eventlist->max = z;
          eventlist->bin = i;
          eventlist->eventTime = (double) i * dt;
	  /* The following line needs to be amended appropriately after
	     confirming what eventTime_ns stands for
	  */
          eventlist->eventTime_ns = (double) i * dt * 2.;
          eventlist->phase = atan2(y,x);
       }

       if (z>findeventsin->Threshold) (*nEvents)++;
   }

   /* eventlist->max /= filter1.length; */
   LALFree(filter1.data);
   LALFree(filter2.data);
   LALFree(output1.data);
   LALFree(output2.data);
   LALFree(buffer.data);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

