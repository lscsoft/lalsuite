/*  <lalVerbatim file="LALEtaTau02CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALEtaTau02.c}}
Given $\tau_0$ and $\tau_2$ one can determine $\eta$ by solving 
\begin{equation}
-\eta^{2/5} \tau_2 + A_2 \left ( \frac {\tau_0}{A_0} \right )^{3/5}  
\left (1 + B_2\eta \right )  = 0,
\end{equation}
where $A_0 = 5/[256 (\pi f_{s} )^{8/3}],$ $A_2 = 3715 / [64512 (\pi f_s)^2],$
$B_2 = 4620/3715.$  
This function returns the LHS of the above
equation in \texttt{x} for a given \texttt{eta}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALEtaTau02CP}
\index{\verb&LALEtaTau02()&}

\subsubsection*{Description}
None.

\subsubsection*{Algorithm}
None.

\subsubsection*{Uses}
None.
</lalLaTeX> */


#include <lal/LALInspiral.h>

/*  <lalVerbatim file="LALEtaTau02CP"> */
void LALEtaTau02(LALStatus *status, 
                 REAL8 *x, 
                 REAL8 eta, 
                 void *p) 
{ /* </lalVerbatim> */
   EtaTau02In *q;
   status = NULL;
   q = (EtaTau02In *) p;
   *x = -q->t2 + q->A2/pow(eta,0.4) * (1. + q->B2*eta);
}
