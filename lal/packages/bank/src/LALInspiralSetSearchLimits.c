/* <lalVerbatim file="LALInspiralSetSearchLimitsCV">
Author: Churches, D. K.
$Id$
</lalVerbatim>  */


/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralSetSearchLimits.c}}

Function which calculates the minimum and maximum values of $\tau_{0}$ and $\tau_{3}$,
as determined by the total mass of the binary $m$ and the symmetric 
mass ratio $\eta$.  The function also calulates the coordinates of the
first template in the bank. These coordinates are $\tau_{0}=\tau_{0min}$,
$\tau_{3}=\tau_{3min}$. 


\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralSetSearchLimitsCP}
\index{\texttt{LALInspiralSetSearchLimits()}}

\subsubsection*{Description}

We start with the definition of the chirp times $\tau_{0}$ and $\tau_{3}$,
\begin{equation}
\tau_{0} = \frac{5}{256 (\pi f_{a} )^{8/3} m^{5/3} \eta}
\end{equation}

and

\begin{equation}
\tau_{3} = \frac{1}{8 (\pi^{2} f_{a}^{3} )^{1/3} m^{2/3} \eta}
\end{equation}

$\tau_{0}$ is minimised when $\eta=1/4$ and $\mathtt{m=MMax}$.
$\tau_{0}$ is maximised when $\eta=1/4$ and $\mathtt{m=2mMin}$.
$\tau_{3}$ is minimised when $\eta=1/4$ and $\mathtt{m=MMax}$.
$\tau_{3}$ is maximised when
$\eta=\mathtt{ mMin(MMax-mMin)/MMax^{2} }$.


\subsubsection*{Algorithm}


\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralSetSearchLimitsCV}}

</lalLaTeX>  */



#include <lal/LALInspiralBank.h>
#include <lal/LALStdlib.h>


NRCSID (LALINSPIRALSETSEARCHLIMITSC, "$Id$");

/*  <lalVerbatim file="LALInspiralSetSearchLimitsCP">  */
void LALInspiralSetSearchLimits (LALStatus               *status,
                                 InspiralBankParams   *bankParams,
                                 InspiralCoarseBankIn coarseIn)
{  /*  </lalVerbatim>  */


   InspiralTemplate Pars1, Pars2, Pars3;

   INITSTATUS (status, "LALInspiralSetSearchLimits", LALINSPIRALSETSEARCHLIMITSC);
   ATTATCHSTATUSPTR(status);

/* Initiate three parameter vectors consistent with the coarseIn structure */

   LALInspiralSetParams(status->statusPtr, &Pars1, coarseIn); CHECKSTATUSPTR(status);
   LALInspiralSetParams(status->statusPtr, &Pars2, coarseIn); CHECKSTATUSPTR(status);
   LALInspiralSetParams(status->statusPtr, &Pars3, coarseIn); CHECKSTATUSPTR(status);
   Pars1.massChoice = Pars2.massChoice = Pars3.massChoice = m1Andm2;

/* Calculate the value of the parameters at the three corners of the search space */

   Pars1.mass1 = coarseIn.MMax/2.;
   Pars1.mass2 = coarseIn.MMax/2.;
   LALInspiralParameterCalc(status->statusPtr, &Pars1); CHECKSTATUSPTR(status);
   Pars2.mass1 = coarseIn.mMin;
   Pars2.mass2 = coarseIn.mMin;
   LALInspiralParameterCalc(status->statusPtr, &Pars2); CHECKSTATUSPTR(status);
   Pars3.mass1 = coarseIn.mMin;
   Pars3.mass2 = (coarseIn.MMax - coarseIn.mMin);
   LALInspiralParameterCalc(status->statusPtr, &Pars3); CHECKSTATUSPTR(status);


/* Find the minimum and maximum values of the parameters and set 
   the search space.  (The minimum values of chirp times are those 
   corresponding to m1 = m2 = MMax/2, i.e., Pars1 structure.
*/

   bankParams->x0 = bankParams->x0Min = Pars1.t0;
   bankParams->x0Max = Pars2.t0;

  switch (coarseIn.space) {
     case Tau0Tau2:
        bankParams->x1 = bankParams->x1Min = Pars1.t2;
        bankParams->x1Max = (Pars2.t2 > Pars3.t2) ? Pars2.t2 : Pars3.t2;
        break;
     case Tau0Tau3:
        bankParams->x1 = bankParams->x1Min = Pars1.t3;
        bankParams->x1Max = (Pars2.t3 > Pars3.t3) ? Pars2.t3 : Pars3.t3;
        break;
     default:
        fprintf(stderr, "LALInspiralSetSearchLimits: Invalid coordinate choice\n");
        exit(0);
        break;
  }
  DETATCHSTATUSPTR(status);
  RETURN (status);
}
