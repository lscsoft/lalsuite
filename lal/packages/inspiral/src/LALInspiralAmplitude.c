/*  <lalVerbatim file="LALInspiralAmplitudeCV">
Author: Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralAmplitude.c}}
Given an inspiral template structure containing the binary distance 
and a set of mass parameters, that module provides functions to compute 
the related amplitude.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralRestrictedAmplitudeCP}
\idx{LALInspiralRestrictedAmplitude()}

\subsubsection*{Description}
The inspiral template structure can stored (1) the distance of the binary (2)
a set of binary masses such as the two masses or the total mass and eta and (3) 
an amplitude which is arbitrary fixed to unity when templates are computed. 

\noindent However we might need to have a template with the physical amplitude (for instance 
to deal with injections). The function \texttt{LALInspiralRestrictedAmplitude}
takes anInspiralTemplate structure as input/output to return the restricted 
Newtonian amplitude by using the following formula. 

\begin{equation}
A = \frac{4c}{d \eta} M^{-5./3.}
\end{equation}

\noindent where $d$ is in Mpc and $M$ in solar mass. The result is stored in the signalAmplitude
variable of the inspiral template structure.



\subsubsection*{Uses}
When appropriate this function calls:\\
\texttt{LALInspiralParameterCalc}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralAmplitudeCV}}
</lalLaTeX>  */



#include <lal/LALInspiral.h>


NRCSID (LALINSPIRALAMPLITUDEC, "$Id$");

/*  <lalVerbatim file="LALInspiralRestrictedAmplitudeCP"> */
void 
LALInspiralRestrictedAmplitude (LALStatus        *status, 
				InspiralTemplate *params )
{ /* </lalVerbatim> */
  
  
  INITSTATUS (status, "LALInspiralAmplitude", LALINSPIRALAMPLITUDEC );
  ATTATCHSTATUSPTR(status);
  
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT((INT4)params->massChoice >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT((INT4)params->massChoice <= 8, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  
  if (params->massChoice != totalMassAndEta) 
    {
      LALInspiralParameterCalc(status->statusPtr, params );
      CHECKSTATUSPTR(status);
    }
  
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT((INT4)params->massChoice >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT((INT4)params->massChoice <= 8, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  

  params->signalAmplitude = 4. * params->totalMass  * params->eta   /  (LAL_PC_SI * 1e6 *params->distance / LAL_MRSUN_SI); 
    
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

