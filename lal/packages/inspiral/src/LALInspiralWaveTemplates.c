/*  <lalVerbatim file="LALInspiralWaveTemplatesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWaveTemplates.c}}
This is exactly the same as \texttt{LALInspiralWave,} except that
it generates two templates one for which the starting phase is 
\texttt{params.startPhase} and a second for which the starting phase is
\texttt{params.startPhase + $\pi/2$}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWaveTemplatesCP}
\idx{LALInspiralWaveTemplates()}
\begin{itemize}
\item {\tt signal1:} Output containing the 0-phase inspiral waveform.
\item {\tt signal2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}
{\tt *signla1} contains `0-phase' inspiral template and {\tt *signal2} contains 
a signal that is $\pi/2$ out of phase with respect to {\tt *signal1.}
Currently, a template pair is generated only for the following {\tt approximants:}
{\tt TaylorT1, TaylorT2, TaylorT3, PadeT1, EOB.}

\subsubsection*{Algorithm}


\subsubsection*{Uses}
Depending on the user inputs one of the following functions is called:\\
\texttt{LALInspiralWave1Templates\\}
\texttt{LALInspiralWave2Templates\\}
\texttt{LALInspiralWave3Templates\\}
\texttt{LALEOBWaveformTemplates}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWaveTemplatesCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALWAVETEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveTemplatesCP"> */
void 
LALInspiralWaveTemplates(
   LALStatus        *status,
   REAL4Vector      *signal1,
   REAL4Vector      *signal2,
   InspiralTemplate *params
   )
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


   ASSERT((INT4)params->approximant >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->approximant <= 7, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->order >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->order <= 7, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

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
      case BCV:
      case SpinTaylorT3:
           ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
      break;
   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}
