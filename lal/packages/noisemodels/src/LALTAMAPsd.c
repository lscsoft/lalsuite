/*  <lalVerbatim file="LALTAMAPsdCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTAMAPsd.c}}

Module to calculate the noise power spectral density for the TAMA detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALTAMAPsdCP}
\index{\texttt{LALTAMAPsd()}}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in units of 400Hz, and it calculates the noise spectral density (per Hz) $S_{h}(f)$ for that frequency $f$.

\begin{equation}
   S_h(f) = \left ( \frac{f}{f_0} \right )^{-5} + 13 \frac{f_0}{f} +
   9 \left [1 + \left( \frac{f}{f_0} \right)^2 \right ].
\end{equation}
\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALTAMAPsdCV}}

</lalLaTeX>  */



#include <lal/LALNoiseModels.h>

/*  <lalVerbatim file="LALTAMAPsdCP"> */
REAL8  LALTAMAPsd(REAL8 x) 
{ /* </lalVerbatim> */

   REAL8 psd, seismic, thermal, shot;

   seismic = pow(x,-5);
   thermal = 13. / x;
   shot = 9. * (1. + x*x);
   psd = seismic + thermal + shot;
   return(psd);
}
