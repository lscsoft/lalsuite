/*  <lalVerbatim file="LALInspiralBankUtilsCV">
Author: Cokelaer Thomas
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralBankUtils.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralBankUtilsCP}
\idx{LALInspiralBankUtils()}

\subsubsection*{Description}
In a parameter space defined by m1 and m2, the equal mass line is used quite often to lay down template. It is particulary
important since most of the algorithm placement won't allow a placement below this line. This  functions returns the value of
$\tau_3$, given the value of $\tau_0$ and the lower cut-off frequency.

\subsubsection*{Algorithm}
We know that 
\begin{equation}
\tau_0 = \frac{A_0}{\eta} M^{-5/2}, 
\end{equation}
and
 \begin{equation}
\tau_3 = \frac{A_3}{\eta} M^{-2/3}, 
\end{equation}
where
\begin{equation}
A_0 = \frac{5}{256 (\pi *f_L)^{8/3}}, 
\end{equation}
and
\begin{equation}
A_3 = \frac{\pi}{8 (\pi *f_L)^{5/3}}, 
\end{equation}

Therefore, if $\eta=0.25$ on the equal-mass line, then it is straightforward to express $\tau_3$ as a function of $\tau_0$ :
\begin{equation}
\tau_3 = 4 A3 \left( \frac{\tau_0}{4 A_0} \right)^{2/5}
\end{equation}

\subsubsection*{Uses}
\begin{verbatim}
LALInspiralParameterCalc()
LALInspiralComputeMetric()
\end{verbatim}

\subsubsection*{Notes}

\vspace{0.1in}
\vfill{\footnotesize\input{LALInspiralBankUtilsCV}}
</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FindRoot.h>


NRCSID(LALINSPIRALBANKUTILSC, "$Id$");


/*  <lalVerbatim file="LALInspiralBankUtilsCP"> */
REAL4
XLALInspiralTau3FromTau0AndEqualMassLine(
    REAL4               tau0,
    REAL4               fL
    )
{   /*  </lalVerbatim>  */
  REAL4 A0, A3, pifL, tau3; 
  
  pifL = LAL_PI * fL;
  A0 = (5.0 / 256.0) * pow(pifL, (-8.0/3.0));
  A3  = LAL_PI / (8.0 * pow(pifL, (5.0/3.0)));

  tau3 = 4 * A3 * pow(tau0/4/A0, 2./5.);
  
  return tau3;
}
  

