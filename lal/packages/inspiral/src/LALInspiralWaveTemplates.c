/*  <lalVerbatim file="LALInspiralWaveTemplatesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWaveTemplates.c}}

Module to calculate two inspiral waveforms which differ only in their phases. The phase difference between
the waveforms is $\pi/2$. This is so that when we maximise over phase by correlating a signal seperately
with a zero and $\pi/2$ phase template, we do not need to call the waveform generation routine twice.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWaveTemplatesCP}
\idx{LALInspiralWaveTemplates()}

\subsubsection*{Description}

This module is identical to \texttt{LALInspiralWave} except that it has an additional output structure, which
corresponds to the $\pi/2$ waveform.

\subsubsection*{Algorithm}


\subsubsection*{Uses}

Depending on the user inputs, any one of the following:

\texttt{TimeDomain2Templates}
\texttt{TappRpnTdomTimeTemplates}
\texttt{TAppRpnTdomFreqTemplates}

\subsubsection*{Notes}

See the documentation for the module \texttt{LALInspiralWave} for further details.

\vfill{\footnotesize\input{LALInspiralWaveTemplatesCV}}

</lalLaTeX>  */


#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALWAVEC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveTemplatesCP"> */
void LALInspiralWaveTemplates(LALStatus *status,
		              REAL4Vector *signal0,
		              REAL4Vector *signal1,
		              InspiralTemplate *params)
{ /* </lalVerbatim>  */

   INITSTATUS(status, "LALInspiralWave", LALINSPIRALWAVEC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal0,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal0->length >= 2, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (signal1->length >= 2, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (signal0->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal1,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal1->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   ASSERT(params->domain >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->domain <= 1, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   switch (params->domain) {
      case TimeDomain:
         ASSERT(params->method >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->method <= 3, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         switch (params->method) {

            case one:
            case best:
               LALTimeDomain2Templates(status->statusPtr, signal0, signal1, params);
               CHECKSTATUSPTR(status);
            break;
            case two:
               LALTappRpnTdomFreqTemplates(status->statusPtr, signal0, signal1, params);
               CHECKSTATUSPTR(status);
            break;
            case three:
               LALTappRpnTdomTimeTemplates(status->statusPtr, signal0, signal1, params);
               CHECKSTATUSPTR(status);
            break;
         }
      break;
      case FrequencyDomain:
         fprintf(stderr,"LALInspiralWaveTemplates: No frequency domain waveforms yet \n");
         exit(0);
   }
   DETATCHSTATUSPTR(status);
   RETURN (status);

}
