/*  <lalVerbatim file="LALGEOPsdCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALGEOPsd.c}}

Module to calculate the noise power spectral density for the GEO600 detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALGEOPsdCP}
\index{\texttt{LALGEOPsd()}}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in units of 150Hz, and it calculates the noise spectral density (per Hz) $S_{h}(f)$ for that frequency $f$.

\begin{equation}
   S_h(f) = 10^{-16} \left ( \frac{f}{f_0} \right)^{-30} + 
            34 \frac{f_0 }{ f } +
   \frac{20 \left [1 - (f/f_0)^2 + 0.5 (f/f_0)^4 \right ] }{ 1 + 0.5 (f/f_0)^2}
\end{equation}
\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALGEOPsdCV}}

</lalLaTeX>  */



#include <lal/LALNoiseModels.h>



/* <lalVerbatim file="LALGEOPsdCP"> */
REAL8  LALGEOPsd(REAL8 x) 
{ /* </lalVerbatim> */

   REAL8 psd, seismic, thermal, shot;
   seismic = pow(10.,-16.) *  pow(x,-30.);
   thermal = 34. / x;
   shot = 20. * (1 - pow(x,2.) + 0.5 * pow(x,4.)) / (1. + 0.5 * pow(x,2.));
   psd = seismic + thermal + shot;
   return(psd);
}
