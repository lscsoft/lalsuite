/*  <lalVerbatim file="LALInspiralWaveLengthCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWaveLength.c}}

Module to calculate the number of data points needed to make a given waveform.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWaveLengthCP}
\idx{LALInspiralWaveLength()}

\subsubsection*{Description}

This module calculates the length of a waveform by evaluating the following integral
\begin{equation}
t(v) = t_{0} - m \int^{v}_{v_{0}} \frac{E^{\prime}(v)}{\mathcal{F}(v)} \, dv \,\,,
\end{equation}
where $t_{0}$ is the starting time for the waveform, $v_{0}$ is its frequency at time $t_{0}$, and $v$
is the frequency of the waveform at time $t(v)$. The upper limit is choosen to be 
$v=v_{lso}$ or a user specified cutoff \texttt{fCutoff} in the parameter structure \texttt{InspiralTemplate}
if that is lower than the LSO.

This gives us the length of the waveform in seconds. If we then multiply this by the sampling
rate then we get the number of samples needed. We then add on the number of zeros which the user
wants to add at the start and end of the waveform, to get the final number of data points needed.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
\begin{verbatim}
InspiralSetup
ChooseModel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWaveLengthCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALWAVELENGTHC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveLengthCP"> */
void LALInspiralWaveLength(
   LALStatus        *status, 
   UINT4            *length,
   InspiralTemplate params) 
{ /* </lalVerbatim>  */

   INT4 ndx;
   REAL8 x;

   expnCoeffs ak;
   expnFunc func;

   INITSTATUS (status, "LALInspiralWaveLength", LALINSPIRALWAVELENGTHC);
   ATTATCHSTATUSPTR(status);

   ASSERT (params.nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params.nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params.tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (&params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

   LALInspiralSetup (status->statusPtr, &ak, &params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, &params);
   CHECKSTATUSPTR(status);

   x = (ak.tn) * params.tSampling + params.nStartPad + params.nEndPad;

   ndx = ceil(log10(x)/log10(2.));
   *length = pow(2, ndx);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
