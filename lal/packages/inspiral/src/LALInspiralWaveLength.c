/*  <lalVerbatim file="LALInspiralWaveLengthCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWaveLength.c}}

Module to calculate the number of data points (to the nearest power of 2)
needed to store a waveform.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWaveLengthCP}
\idx{LALInspiralWaveLength()}
\begin{itemize}
\item {\tt length:} output, number of bins to the nearest power of 2
greater than the minimum length required to store a wave of parameters as in {\tt params}.
\item {\tt params:} input, parameters of the binary system. 
\end{itemize}

\subsubsection*{Description}


This module first calls {\tt LALInspiralChooseModel,} which gives the length of the 
waveform in seconds. Multiplying this by the sampling rate {\tt params.tSampling}
gives the minimum number of samples needed to hold the waveform. To this are added the
number of bins of leading and trailing zeroes requested by the user in 
{\tt params.nStartPad} and {\tt params.nEndPad.} The resulting number is rounded to
an upward power of 2 and returned in {\tt length}. Note that this {\tt length} is good enough
for all waveforms except EOB.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
This function calls:\\
\texttt{
LALInspiralSetup\\
LALInspiralChooseModel
}

\subsubsection*{Notes}
\begin{itemize}
\item Sometimes waveforms of type {\tt EOB} might require an array larger than
the length returned by this function. It is desirable to use a vector
of size 2 times greater than returned by this code for EOB waveforms.

\item It would be safer if a call to the {\tt LALInspiralTemplate} is made
before calling this function.
\end{itemize}
\vfill{\footnotesize\input{LALInspiralWaveLengthCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

#define INSPIRAL_MINIMALWAVELENGTH 64

NRCSID (LALINSPIRALWAVELENGTHC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveLengthCP"> */
void 
LALInspiralWaveLength(
   LALStatus        *status, 
   UINT4            *length,
   InspiralTemplate  params
   ) 
{ /* </lalVerbatim>  */

   INT4 ndx;
   REAL8 x;
   CHAR message[256];
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
   /* not checkstatus before but we check if everything is fine, 
      then we return the length otherwise length is zero*/
   if (status->statusPtr->statusCode == 0){
	   /*we add a minimal value and 10 % of overestimation */
      	   x	= 1.1 * (ak.tn) * params.tSampling + params.nStartPad 
   		   + params.nEndPad + INSPIRAL_MINIMALWAVELENGTH;
	   ndx 	= ceil(log10(x)/log10(2.));
      	   *length = pow(2, ndx);
     
	DETATCHSTATUSPTR(status);
	RETURN(status);
   }
   else { 
	sprintf(message,
	       	"size is zero for the following waveform: totalMass = %lf, fLower = %lf, approximant = %d @ %lfPN"
	       , params.mass1 + params.mass2, params.fLower, params.approximant, params.order/2.);
       	LALWarning(status, message);
       	*length = 0;
       
       DETATCHSTATUSPTR(status);
       RETURN(status);
     }
}
