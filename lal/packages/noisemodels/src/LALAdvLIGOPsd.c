/*  <lalVerbatim file="LALAdvLIGOPsdCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALAdvLIGOPsd.c}}

Module to calculate the noise power spectral density for the initial LIGO detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALAdvLIGOPsdCP}
\idx{LALAdvLIGOPsd()}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in Hz, and it 
calculates the noise spectral density (per Hz) $S_{h}(f)$ 
for that frequency. The noise PSD is based on data provided by
K. Blackburn (see T. Damour, B.R. Iyer and B.S. Sathyaprakash,
Phys. Rev. D 63, 044023 (2001)) and is approximated by
the following:
\begin{equation}
   S_h(f) = \left ( \frac {4.49 f}{f_0} \right )^{-56} + 
            0.16 \left ( \frac{f}{f_0} \right )^{-4.52} + 0.52 + 
            0.32 \left ( \frac {f}{f_0} \right )^2
\end{equation}
The returned value is scaled up by $s_0 = 10^{46}/9.$ In otherwords, 
the expected noise PSD is $9 \times 10^{-46}$ times the returned value.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALAdvLIGOPsdCV}}

</lalLaTeX>  */

#include <lal/LALNoiseModels.h>

/*  <lalVerbatim file="LALAdvLIGOPsdCP"> */
void
LALAdvLIGOPsd (LALStatus *status, REAL8 *psd, REAL8 f) 
{ /* </lalVerbatim> */

   REAL8 x2,x;
   x = f/215.;
   status = NULL;
   x2 = x*x;
   *psd = pow(x,-4.14) - 5./x2 + 111. * (1. - x2 + 0.5 * x2*x2)/(1. + 0.5*x2);
}
