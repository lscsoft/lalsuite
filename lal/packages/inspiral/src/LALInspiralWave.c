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
\index{\texttt{LALInspiralWave()}}

\subsubsection*{Description}

The code \texttt{LALInspiralWave} is the user interface to the inspiral codes. It takes from the user all
the physical parameters which specify the binary, and calls the relevent code generation function.
At the moment we have three methods of generating waveforms. These are set by the enum \texttt{Method}, and
the choices available represent the following possibilities.

The choices available represent the following possibilities.
 
Theoretical calculations have given us post--Newtonian expansions (i.e.\ a series in terms of ascending
powers of $x\equiv v^{2}$) of an energy function $E(x)$ and a GW luminosity function $\mathcal{F}(v)$,
together which define the gravitational wave phasing formula:
 
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
corresponds to method \texttt{one} in the enum \texttt{Method} above.
 
\item Re--expand $E^{\prime}(v)/\mathcal{F}(v)$ in which case the integrals can be done analytically to
obtain
a {\it parametric} representation of the phasing formula in terms of
polynomial expressions in the auxiliary variable $v$
\begin {eqnarray}
\phi^{(2)}(v)&=& \phi_{\rm ref} +
\phi^v_N (v)\sum_{k=0}^{n} \hat{\phi}^v_k v^k, \nonumber\\
t^{(2)}(v)&=& t_{\rm ref} +t^v_N(v) \sum_{k=0}^{n} \hat{t}^v_k v^k,
\label{InspiralWavePhase2}
\end {eqnarray}
This corresponds to method \texttt{two} in the enum \texttt{Method} above.
 
\item the second of the polynomials in Eq.~(\ref{InspiralWavePhase2}) can
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
This corresponds to method \texttt{three} in the enum \texttt{Method} above.
 
\end{enumerate}

If the user chooses method \texttt{two} or \texttt{three}, then they also need to choose one of the
following ways of inputting the mass of the binary, which is handled by the parameter \texttt{MassChoice},
which is of type \texttt{enum InputMasses}.

The choice \texttt{m1Andm2} corresponds to choosing to declare the masses $m_{1}$ and $m_{2}$ of the two
objects in the binary. The choice \texttt{totalMassAndEta} corresponds to declaring $m=m_{1}+m_{2}$ and
$\eta=m_{1}m_{2}/(m_{1}+m_{2})^{2}$, and the choice \texttt{totalMassAndMu} corresponds to declaring
$m=m_{1}+m_{2}$ and $\mu=m_{1}m_{2}/(m_{1}+m_{2})$ as inputs. Each choice is completely equivalent, so user
may choose the most convenient pair of inputs.
 
The parameter \texttt{order}, which is of type \texttt{enum Order}, lets the user choose to which order of
post--Newtonian expansion they would like to go.

The parameter \texttt{domain}, which is of type \texttt{enum Domain}, lets the user to choose to generate
the waveform in the time or frequency domain.

If the user has chosen method \texttt{one}, (i.e.\ they have chosen to integrate
$E^{\prime}(v)/\mathcal{F}(v)$ numerically) then they have the additional choice of specifying the form of
the expansion which defines $E^{\prime}(v)/\mathcal{F}(v)$. This ratio can be expressed in the form of a
\emph{Taylor series} (T--Approximants) or by using \emph{P--Approximants}. This choice is handled by the
parameter \texttt{approximant}, which is of type \texttt{enum Approximant}.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
Depending on the user inputs, and one of the following functions:

\texttt{TimeDomain2}
\texttt{TappRpnTdomFreq}
\texttt{TappRpnTdomTime}

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

   ASSERT (signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   switch (params->domain) {
	case TimeDomain:
		switch (params->method) {

			case one:
			case best:
		        	LALTimeDomain2(status->statusPtr, signal, params);
   	                	CHECKSTATUSPTR(status);
				break;
			case two:
				LALTappRpnTdomFreq(status->statusPtr, signal, params);
   	                	CHECKSTATUSPTR(status);
				break;
			case three:
				LALTappRpnTdomTime(status->statusPtr, signal, params);
   	                	CHECKSTATUSPTR(status);
				break;
			default:
		                fprintf(stderr,"LALInspiralWave: You haven't chosen a method ... exiting\n");
                                exit(0);
				}
		break;

	case FrequencyDomain:
		fprintf(stderr,"LALInspiralWave: We don't have frequency domain waveforms yet \n");
                exit(0);
   }						

   DETATCHSTATUSPTR(status);
   RETURN (status);

}
