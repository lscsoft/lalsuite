/*  <lalVerbatim file="LALInspiralWaveCV">
Author: Churches, D. K and Sathyaprakash, B.S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave.c}}

Interface routine needed to generate a T- or a P-approximant.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWaveCP}
\idx{LALInspiralWave()}

\subsubsection*{Description}

The code \texttt{LALInspiralWave} is the user interface to the inspiral codes. It takes from the user all
the physical parameters which specify the binary, and calls the relevent code generation function.
We have plans to produce eleven different approximants of which six are fully implemented. 
Each approximant can be generated at seven different post-Newtonian orders, from Newtonian to 3.5 PN
order.  The approximant and the order are set up the enums \texttt{Approximant} and \texttt{Order,}
respectively.

The choices available represent the following possibilities.
 
Theoretical calculations have given us post-Newtonian expansions (i.e.\ a series in terms of ascending
powers of $x\equiv v^{2}$) of an energy function $E(x=v^2)$ and a GW luminosity function $\mathcal{F}(v)$.
Then there is the gravitational wave phasing formula at the restricted post-Newtonian approximation:
\begin{eqnarray}
t(v) & = & t_{\rm ref} + m \int_v^{v_{\rm ref}} \,
\frac{E'(v)}{{\cal F}(v)} \, dv, \nonumber \\
\phi (v) & = & \phi_{\rm ref} + 2 \int_v^{v_{\rm ref}}  v^3 \,
\frac{E'(v)}{{\cal F}(v)} \, dv,
\label{InspiralWavePhasingFormula}
\end{eqnarray}
where $v=(\pi m F)^{1/3}$ is an invariantly defined velocity, $F$ is the instantaneous GW frequency, and
$m$ is the total mass of the binary.
 
The variable $t$ is in our hands, and the phasing formula is solved by finding the value of $v$ (for a
given $t$) which satisfies the first of the equations in Eq. (\ref{InspiralWavePhasingFormula}). This value of $v$ is
then substituted into the second equation to yield $\phi(t)$.
 
There are basically three ways of solving the problem:

\begin{enumerate}
\item Leave  $E^{\prime}(v)/\mathcal{F}(v)$ as it is and integrate the equations numerically. This
corresponds to \texttt{TaylorT1} and \texttt{PadeT1} approximants in the enum \texttt{Approximant}.
 
\item Re-expand $E^{\prime}(v)/\mathcal{F}(v)$ in which case the integrals can be done analytically to
obtain a {\it parametric} representation of the phasing formula in terms of
polynomial expressions in the auxiliary variable $v$
\begin {eqnarray}
\phi^{(2)}(v)&=& \phi_{\rm ref} +
\phi^v_N (v)\sum_{k=0}^{n} \hat{\phi}^v_k v^k, \nonumber\\
t^{(2)}(v)&=& t_{\rm ref} +t^v_N(v) \sum_{k=0}^{n} \hat{t}^v_k v^k,
\label{InspiralWavePhase2}
\end {eqnarray}
This corresponds to approximant \texttt{TaylorT2} in the enum \texttt{Approximant}.
 
\item The second of the polynomials in Eq.~(\ref{InspiralWavePhase2}) can
be inverted and the resulting polynomial for $v$ in terms of
$t$ can be substituted in $\phi^{(2)}(v)$ to arrive at an explicit  time-domain
phasing formula
\begin{equation}
\phi^{(3)}(t)=\phi_{\rm ref}+\phi_N^t \sum_{k=0}^{n}
\hat{\phi}^t_k\theta^k
\end{equation}
\begin{equation}
F^{(3)}(t)= F_N^t \sum_{k=0}^{n} \hat{F}^t_k \theta^k,
\end{equation}         
where \\$\theta=[\eta (t_{\rm ref}-t)/(5m)]^{-1/8}$, \\
$F \equiv d \phi/ 2 \pi dt =v^3/(\pi m)$ is the instantaneous GW frequency.
This corresponds to method \texttt{TaylorT3} in the enum \texttt{approximant} above.
\end{enumerate}

Approximants \texttt{EOB} and \texttt{DJS} correspond to effective one-body approach
where a set of four ordinary differential equations are solved by making an anzatz for the
radiation reaction force. The four ODEs correspond to the evolution of the radial and angular
coordinates and the corresponding momenta.

See the test code in this module for an example of how to generate a waveform.

The parameter \texttt{Order}, which is of type \texttt{enum Order}, lets the user choose to which order of
post-Newtonian expansion they would like to go.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
Depending on the user inputs one of the following functions is called:

\texttt{LALInspiralWave1}
\texttt{LALInspiralWave2}
\texttt{LALInspiralWave3}
\texttt{LALEOBWaveform}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWaveCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALWAVEC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveCP"> */
void LALInspiralWave(LALStatus *status,
		     REAL4Vector *signal,
		     InspiralTemplate *params)
{ /* </lalVerbatim>  */

   INITSTATUS(status, "LALInspiralWave", LALINSPIRALWAVEC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal->length >= 2, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   ASSERT((INT4)params->approximant >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->approximant <= 11, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->order >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->order <= 7, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   switch (params->approximant) 
   {
      case TaylorT1:
      case PadeT1:
           LALInspiralWave1(status->statusPtr, signal, params);
           CHECKSTATUSPTR(status);
	   break;
      case TaylorT2:
           LALInspiralWave2(status->statusPtr, signal, params);
           CHECKSTATUSPTR(status);
	   break;
      case TaylorT3:
           LALInspiralWave3(status->statusPtr, signal, params);
           CHECKSTATUSPTR(status);
	   break;
      case EOB:
           LALEOBWaveform(status->statusPtr, signal, params);
           CHECKSTATUSPTR(status);
	   break;
      case TaylorF1:
	   LALInspiralStationaryPhaseApprox1(status->statusPtr, signal, params); 
           CHECKSTATUSPTR(status);
	   break;
      case TaylorF2:
	   LALInspiralStationaryPhaseApprox2(status->statusPtr, signal, params); 
           CHECKSTATUSPTR(status);
	   break;
      case PadeF1:
      case INSPA:
      case IRSPA:
      case DJS:
           ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
	   break;
      case SpinTaylorT3:
	   LALInspiralSpinModulatedWave(status->statusPtr, signal, params); 
           CHECKSTATUSPTR(status);
	   break;
   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}
