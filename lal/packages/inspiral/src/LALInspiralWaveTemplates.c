/*  <lalVerbatim file="LALInspiralWaveTemplatesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWaveTemplates.c}}
This is exactly the same as \texttt{LALInspiralWave1,} except that
it generates two templates one for which the starting phase is 
\texttt{params.startPhase} and the other for which the phase is
\texttt{params.startPhase+$\pi/2$}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWaveTemplatesCP}
\idx{LALInspiralWaveTemplates()}

\subsubsection*{Description}

\subsubsection*{Algorithm}


\subsubsection*{Uses}
Depending on the user inputs one of the following functions is called:

\texttt{LALInspiralWave1Templates}
\texttt{LALInspiralWave2Templates}
\texttt{LALInspiralWave3Templates}
\texttt{LALEOBWaveformTemplates}
\texttt{LALDJSWaveformTemplates}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWaveTemplatesCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALWAVETEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveTemplatesCP"> */
void LALInspiralWaveTemplates(LALStatus *status,
		     REAL4Vector *signal1,
		     REAL4Vector *signal2,
		     InspiralTemplate *params)
{ /* </lalVerbatim>  */

   INITSTATUS(status, "LALInspiralWaveTemplates", LALINSPIRALWAVETEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal1->length >= 2, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (signal2->length >= 2, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (signal1,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal2,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal1->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal2->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   ASSERT(params->approximant >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->approximant <= 7, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->order >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->order <= 7, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   switch (params->approximant) 
   {
      case TaylorT1:
      case PadeT1:
           LALInspiralWave1Templates(status->statusPtr, signal1, signal2, params);
           CHECKSTATUSPTR(status);
      break;
      case TaylorT2:
           LALInspiralWave2Templates(status->statusPtr, signal1, signal2, params);
           CHECKSTATUSPTR(status);
      break;
      case TaylorT3:
           LALInspiralWave3Templates(status->statusPtr, signal1, signal2, params);
           CHECKSTATUSPTR(status);
      break;
      case EOB:
           LALEOBWaveformTemplates(status->statusPtr, signal1, signal2, params);
           CHECKSTATUSPTR(status);
      break;
      case TaylorF1:
      case TaylorF2:
      case PadeF1:
      case INSPA:
      case IRSPA:
           ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
      break;
   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}
