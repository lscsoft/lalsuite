/*  <lalVerbatim file="LALInspiralVelocityCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralVelocity.c}}

The function \texttt{LALInspiralVelocity} calculates the velocity $v$ which corresponds to a time $t$ in
the inspiralling binary system.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralVelocityCP}
\idx{LALInspiralVelocity()}

\subsubsection*{Description}

The function \texttt{LALInspiralVelocity} calculates the velocity $v$ corresponding to a time $t$ 
in the evolution of an inspiralling binary system.  It does this by iteratively solving
\begin{equation}
t(v) =  t_{0} - m \int_{v_{0}}^{v} \frac{E'(v)}{{\cal F}(v)} \, dv \,\,.
\label{tofv}
\end{equation}
\texttt{LALInspiralVelocity} calculates $v$, given $t(v)$, 
$t_{0}$, $m$, $v_{0}$, $E^{\prime}(v)$ and $\mathcal{F}(v)$.

\subsubsection*{Algorithm}


\subsubsection*{Uses}

\texttt{LALDBisectionFindRoot}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralVelocityCV}}

</lalLaTeX>  */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>

NRCSID (LALINSPIRALVELOCITYC, "$Id$"); 

/*  <lalVerbatim file="LALInspiralVelocityCP"> */
void 
LALInspiralVelocity(
   LALStatus *status,
   REAL8     *v,
   TofVIn    *ak
   )
{ /* </lalVerbatim>  */

  DFindRootIn rootIn;
  void *funcParams;


  INITSTATUS (status, "LALInspiralVelocity", LALINSPIRALVELOCITYC);
  ATTATCHSTATUSPTR(status);

  ASSERT (v, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  rootIn.function = LALInspiralTofV;
  rootIn.xmax = ak->vlso;
  rootIn.xmin = ak->v0/2.;
  rootIn.xacc = 1.0e-8;

  funcParams = (void *) ak;


  if (ak->t==ak->t0) 
  {
     *v = ak->v0;
     DETATCHSTATUSPTR(status);
     RETURN(status);
  }

  LALDBisectionFindRoot(status->statusPtr, v, &rootIn, funcParams);
  CHECKSTATUSPTR(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
