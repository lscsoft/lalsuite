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
\index{\verb&LALVIRGOPsd()&}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in units of 500Hz, and it calculates the noise spectral density (per Hz) $S_{h}(f)$ for that frequency $f$.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALVIRGOPsdCV}}

</lalLaTeX>  */



#include <lal/LALInspiralBank.h>

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
