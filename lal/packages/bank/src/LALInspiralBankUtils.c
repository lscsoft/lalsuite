/*
*  Copyright (C) 2007 Thomas Cokelaer
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

/*  <lalVerbatim file="LALInspiralBankUtilsCV">
Author: Cokelaer Thomas
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralBankUtils.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{XLALInspiralTau3FromTau0AndEqualMassLineCP}
\input{XLALInspiralTau3FromNonEqualMassCP}
\input{XLALInspiralTau0FromMEtaCP}
\input{XLALInspiralMFromTau0AndNonEqualMassCP}

\idx{LALInspiralBankUtils()}

\subsubsection*{Description}
In a parameter space defined by $m_1$ and $m_2$, or equivalently, $M=m_1+m_2$ and $\eta=\frac{m_1 m_2}{M^2}$, the conversion
to chirp-time parameter such as $\tau_0$ and $\tau_3$ si quite common. In particular, it is interesting to get the value of
$\tau_3$ when only $\tau_0$ is known, and a constraint on the masses exists (e.g., $m_1=m_2$ or one of the mass equals mMin or mMax.
This modules contains a few functions to perform these conversion.
\subsubsection*{Algorithm}
We know that
\begin{equation}
\tau_0 = \frac{A_0}{\eta} M^{-5/2},
\label{eq:tau0a}
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

Therefore, it is straightforward to express $\tau_3$ as a function of $\tau_0$ amd $\eta$:
\begin{equation}
\tau_3 = \frac{A3}{\eta} \left( \frac{\tau_0 \eta}{ A_0} \right)^{2/5}
\label{eq:tau3b}
\end{equation}
if $\eta=0.25$ on the equal-mass line, then
\begin{equation}
\tau_3 = 4 A3 \left( \frac{\tau_0}{ 4 A_0} \right)^{2/5}
\label{eq:tau3a}
\end{equation}


\noindent Equation \ref{eq:tau3b} returns $\tau_3$ given in $M, \eta$ and $f_L$ and is defined
in\texttt{XLALInspiralTau3FromNonEqualMassLine}.
\\\\
Equation \ref{eq:tau3a} returns tau3 in the particular case $m_1=m_2$, given
$\tau_0$ only, and is defined in \texttt{XLALInspiralTau3FromTau0AndEqualMassLine}.
 \\\\
Equation \ref{eq:tau0a} returns $tau_0$ given $M, \eta$ and $f_L$, and is defined
\texttt{XLALInspiralTau0FromMEta}.
\\\\
Finally, \texttt{XLALInspiralMFromTau0AndNonEqualMass} returns $M$ when $\tau_0$ is known
 and a constraint exists on one of the individual mass (e.g., $m_1={\rm mMax}$ or
 $m_1={\rm mMin}$). This functions requires a little more algebra and is used in the
HybridHexagonal placement. The function LALInspiralHybridHexagonal describes this algebra.

\vspace{0.1in}
\vfill{\footnotesize\input{LALInspiralBankUtilsCV}}
</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FindRoot.h>
#include <math.h>


NRCSID(LALINSPIRALBANKUTILSC, "$Id$");



/*  <lalVerbatim file="XLALInspiralTau3FromTau0AndEqualMassLineCP"> */
REAL4
XLALInspiralTau3FromTau0AndEqualMassLine(
    REAL4               tau0,
    REAL4               fL
    )
{   /*  </lalVerbatim>  */
  REAL4 A0, A3, tau3=0;


  A0 = (5.0 / 256.0) * pow(LAL_PI * fL, (-8.0/3.0));
  A3  = LAL_PI / (8.0 * pow(LAL_PI * fL, (5.0/3.0)));

  tau3 = 4 * A3 * pow(tau0/4/A0, 2./5.);

  return tau3;
}


/*  <lalVerbatim file="XLALInspiralTau3FromNonEqualMassCP"> */
REAL4
XLALInspiralTau3FromNonEqualMass(
    REAL4              	M,
    REAL4 		eta,
    REAL4		fL
 )
{ /*  </lalVerbatim>  */
  REAL4 A3;
  REAL4 tau3 = 0;

  A3  = LAL_PI / (8.0 * pow(LAL_PI*fL, (5.0/3.0)));
  tau3 = A3 * pow(M * LAL_MTSUN_SI, -2.0/3.0) / eta;

  return tau3;
}

/*  <lalVerbatim file="XLALInspiralTau0FromMEtaCP"> */
REAL4
XLALInspiralTau0FromMEta(
    REAL4              	M,
    REAL4 		eta,
    REAL4		fL
 )
{/*  </lalVerbatim>  */

/* This function returns tau3, computed from M and eta*/

  REAL4 A0;
  REAL4 tau0 = 0;

  A0 = (5.0 / 256.0) * pow( LAL_PI * fL, (-8.0/3.0));
  tau0 = A0 * pow(M*LAL_MTSUN_SI, -5.0/3.0) / eta;

  return tau0;
}

/*  <lalVerbatim file="XLALInspiralMFromTau0AndNonEqualMassCP"> */
REAL8
XLALInspiralMFromTau0AndNonEqualMass(
  REAL8 tau0,
  REAL8 extremMass,
  REAL8 fL)
{/*  </lalVerbatim>  */
  REAL8 result, A0, p, q, x;

  A0 = (5.0 / 256.0) * pow( LAL_PI * fL, (-8.0/3.0));

  /* from tau0, and M, we can get a poylomial expression where M is the
  unknowm of the form x^3+px+q =0 where x = M^(1/3) and p and q as follows :*/
  p = -A0/tau0/extremMass/LAL_MTSUN_SI;
  q = -extremMass * LAL_MTSUN_SI;

  x = pow((-q/2-0.5*sqrt((27*q*q + 4*p*p*p)/27)), 1./3.);
  x += pow((-q/2+0.5*sqrt((27*q*q + 4*p*p*p)/27)), 1./3.);

  /* This is a real solution and M is simply */
  result = x*x*x/LAL_MTSUN_SI;

  return result;
}
