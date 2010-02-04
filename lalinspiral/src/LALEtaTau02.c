/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, B.S. Sathyaprakash
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*  <lalVerbatim file="LALEtaTau02CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALEtaTau02.c}}
Given $\tau_0$ and $\tau_2$ compute the mass ratio $\eta.$
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALEtaTau02CP}
\idx{LALEtaTau02()}

\subsubsection*{Description}
Given $\tau_0$ and $\tau_2$ one can determine $\eta$ by solving
\begin{equation}
-\eta^{2/5} \tau_2 + A_2 \left ( \frac {\tau_0}{A_0} \right )^{3/5}
\left (1 + B_2\eta \right )  = 0,
\end{equation}
where $A_0 = 5/[256 (\pi f_{s} )^{8/3}],$ $A_2 = 3715 / [64512 (\pi f_s)^2],$
$B_2 = 4620/3715.$
This function returns the LHS of the above
equation in \texttt{x} for a given \texttt{eta}.

\subsubsection*{Algorithm}
None.

\subsubsection*{Uses}
None.

\subsubsection*{Notes}
The {\tt void pointer} {\tt *p} should point to a {\tt struct}
of type {\tt EtaTau02In:}\\[10pt]
{\tt
void *p;\\
EtaTau02In q;\\[5pt]

$\ldots$\\
p = (void *) \&q;\\
}

</lalLaTeX> */


#include <lal/LALInspiral.h>

NRCSID (LALETATAU02C, "$Id$");
/*  <lalVerbatim file="LALEtaTau02CP"> */
void
LALEtaTau02(
   LALStatus *status,
   REAL8     *x,
   REAL8     eta,
   void      *p
   )
{ /* </lalVerbatim> */
   EtaTau02In *q;

   INITSTATUS(status, "LALEtaTau02", LALETATAU02C);
   ATTATCHSTATUSPTR(status);
   ASSERT (p,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(eta > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   q = (EtaTau02In *) p;
   *x = -q->t2 + q->A2/pow(eta,0.4) * (1. + q->B2*eta);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
