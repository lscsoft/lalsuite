/*  <lalVerbatim file="LALLIGOIPsdCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALLIGOIPsd.c}}

Module to calculate the noise power spectral density for the initial LIGO detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALLIGOIPsdCP}
\index{\verb&LALLIGOIPsd()&}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in units of 175 Hz, 
and it calculates the noise spectral density (per Hz)
$S_{h}(f)$ for that frequency $f$.
\begin{equation}
   S_h(f) = \left ( \frac {4.49 f}{f_0} \right )^{-56} + 
            0.16 \left ( \frac{f}{f_0} \right )^{-4.52} + 0.52 + 
            0.32 \left ( \frac {f}{f_0} \right )^2
\end{equation}

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALLIGOIPsdCV}}

</lalLaTeX>  */

#include <lal/LALNoiseModels.h>

/*  <lalVerbatim file="LALLIGOIPsdCP"> */
void
LALLIGOIPsd (LALStatus *status, REAL8 *psd, REAL8 f) 
{ /* </lalVerbatim> */

   REAL8 x2,x;
   x = f/175.;
   status = NULL;
   x2 = x*x;
   *psd = pow(4.49*x,-56.) + 0.16 * pow(x,-4.52) + 0.52 + 0.32 * x2; 
}
