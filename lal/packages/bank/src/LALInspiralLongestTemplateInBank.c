/*  <lalVerbatim file="LALInspiralLongestTemplateInBankCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralLongestTemplateInBank.c}}
To find the longest template in a template bank.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralLongestTemplateInBankCP}
\idx{LALInspiralLongestTemplateInBank()}

\subsubsection*{Description}
Given the parameters of a template bank find the longest template
in the bank. This is done by looking at the duration for which
a signal corresponding to smallest masses lasts. One simply calls
the {\tt LALInspiralWaveLength} code for a system consisting
of two stars each of mass {\tt mMin.}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALInspiralParameterCalc 
LALInspiralWaveLength 
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralLongestTemplateInBankCV}}

</lalLaTeX>  */

#include <lal/LALInspiralBank.h>
 
NRCSID (INSPIRALSEARCHC, "$Id$");

/*  <lalVerbatim file="LALInspiralLongestTemplateInBankCP"> */
void
LALInspiralLongestTemplateInBank
   (
   LALStatus            *status, 
   UINT4                *templateLength,
   InspiralCoarseBankIn *coarseIn
   )
{  /*  </lalVerbatim>  */

   InspiralTemplate param;
   INITSTATUS (status, "LALInspiralLongestTemplateInBank", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);
   ASSERT (coarseIn,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   param.startTime = 0.0;
   param.startPhase = 0.0;
   param.nStartPad = 0;
   param.nEndPad = 0;
   param.ieta = 1;
   param.Theta = 0.;
   param.OmegaS = 0.;
   param.mass1 = coarseIn->mMin;
   param.mass2 = coarseIn->mMin;
   param.fLower = coarseIn->fLower;
   param.fCutoff = coarseIn->fUpper;
   param.tSampling = coarseIn->tSampling;
   param.signalAmplitude = 1.;
   param.order = coarseIn->order;
   param.approximant = coarseIn->approximant;
   param.massChoice = m1Andm2;
   LALInspiralParameterCalc (status->statusPtr, &param);
   CHECKSTATUSPTR(status);
   *templateLength = 0;
   LALInspiralWaveLength (status->statusPtr, templateLength, param);
   CHECKSTATUSPTR(status);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

