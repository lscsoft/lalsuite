/*  <lalVerbatim file="LALInspiralSetParamsCV">
Author: Churches, D. K and Sathyaprakash, B.S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralSetParams.c}}

A routine that fills an \texttt{InspiralTemplate} structure
based on the values in the \texttt{InspiralCoarseBankIn} structure.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralSetParamsCP}
\idx{LALInspiralSetParams()}

\subsubsection*{Description}

This function takes as an input a structure of type 
\texttt{InspiralCoarseBankIn} and it fills up the
elements of a structure of type \texttt{InspiralTemplate}.
The function sets the fields 
\texttt{massChoice}, \texttt{ieta}, \texttt{signalAmplitude},
\texttt{tSampling}, \texttt{fLower}, \texttt{fCutoff}, 
\texttt{method}, \texttt{order},
\texttt{approximant}, \texttt{domain}, \texttt{nStartPad}, \texttt{nEndPad}.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
\begin{verbatim}
None
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralSetParamsCV}}

</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALInspiralBank.h>

NRCSID (LALINSPIRALSETPARAMSC, "Id: $");

/*  <lalVerbatim file="LALInspiralSetParamsCP"> */

void LALInspiralSetParams(LALStatus            *status, 
                          InspiralTemplate     *tempPars, 
                          InspiralCoarseBankIn coarseIn) 
{  /*  </lalVerbatim>  */

   INITSTATUS (status, "LALInspiralSetParams", LALINSPIRALSETPARAMSC);
   ATTATCHSTATUSPTR(status);
   ASSERT (tempPars,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (tempPars->massChoice >= 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (tempPars->massChoice <= 6, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

   switch (coarseIn.space) {
      case Tau0Tau2:
         tempPars->massChoice = t02;
         break;
      case Tau0Tau3:
         tempPars->massChoice = t03;
         break;
   }

   tempPars->ieta = 1;
   tempPars->signalAmplitude = 1.;
   tempPars->tSampling = coarseIn.tSampling;
   tempPars->fLower = coarseIn.fLower;
   tempPars->fCutoff = coarseIn.fUpper;
   tempPars->method = coarseIn.method;
   tempPars->order = coarseIn.order;
   tempPars->approximant = coarseIn.approximant;
   tempPars->domain = coarseIn.domain;
   tempPars->nStartPad = 0;
   tempPars->nEndPad = 0;

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
