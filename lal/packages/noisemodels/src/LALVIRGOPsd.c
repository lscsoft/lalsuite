/*  <lalVerbatim file="LALVIRGOPsdCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALVIRGOPsd.c}}

Module to calculate the noise power spectral density for the VIRGO detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALVIRGOPsdCP}
\index{\texttt{LALVIRGOPsd()}}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in units of 500Hz, and it calculates the noise spectral density (per Hz) $S_{h}(f)$ for that frequency $f$.

\begin{equation}
   S_h(f) = 
   s_1 \left ( \frac {10f}{f_0} \right )^{-5} + s_2 \frac{f_0}{f}
   + s_3 \left [1 + \left (\frac {f}{f_0} \right)^2 \right ],
\end{equation}
where $s_1=34.6,$ $s_2=6.6,$ and $s_3=3.24.$
\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALVIRGOPsdCV}}

</lalLaTeX>  */



#include <lal/LALNoiseModels.h>

/*  <lalVerbatim file="LALVIRGOPsdCP"> */
REAL8 LALVIRGOPsd (REAL8 x) 
{ /* </lalVerbatim> */

   REAL8 s1, s2, s3, psd;

   s1 = 34.6;
   s2 = 6.60;
   s3 = 3.24;
   psd = s1*pow(10.*x,-5.) + s2/x + s3 * (1. + x*x);
   return psd;
}
