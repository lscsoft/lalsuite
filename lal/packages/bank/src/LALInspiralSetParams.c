/*  <lalVerbatim file="LALInspiralSetParamsCV">
Author: Churches, D. K and Sathyaprakash, B.S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralSetParams.c}}

Module to create a coarse grid of templates.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralSetParamsCP}
\index{\verb&LALInspiralSetParams()&}

\subsubsection*{Description}

This function takes as an input a structure of type \texttt{InspiralCoarseBankIn} and it fills up the
elements of a structure of type \texttt{InspiralTemplate}.
The function sets the fields \texttt{massChoice}, \texttt{ieta}, \texttt{signalAmplitude},
\texttt{tSampling}, \texttt{fLower}, \texttt{fCutoff}, \texttt{method}, \texttt{order},
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
void LALInspiralSetParams(LALStatus               *status, 
                          InspiralTemplate     *tempPars, 
                          InspiralCoarseBankIn coarseIn) 
{  /*  </lalVerbatim>  */

   INITSTATUS (status, "LALInspiralSetParams", LALINSPIRALSETPARAMSC);
   switch (coarseIn.space) {
      case Tau0Tau2:
         tempPars->massChoice = t02;
         break;
      case Tau0Tau3:
         tempPars->massChoice = t03;
         break;
      default:
         fprintf(stderr, "InspiralSetParams: Invalid coordinate space choice\n");
         exit(0);
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
   RETURN(status);
}
