/*  <lalVerbatim file="LALEGOPsdCV">
Author: Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALEGOPsd.c}}

Module to calculate the noise power spectral density for the EGO detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALEGOPsdCP}
\idx{LALEGOPsd()}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in Hz, and it 
calculates the noise spectral density (per Hz) $S_{h}(f)$ 
for that frequency. The noise PSD is based on data provided by
grqc/0607092
\begin{equation}
   S_h(f) = 
   s_0 \left (  to be completed \right ),
\end{equation}
where $s_0=1.61e-51$
\subsubsection*{Algorithm}
\subsubsection*{Uses}
None.
\subsubsection*{Notes}
\vfill{\footnotesize\input{LALEGOPsdCV}}
</lalLaTeX>  */

#include <lal/LALNoiseModels.h>

/*  <lalVerbatim file="LALEGOPsdCP"> */
void LALEGOPsd (LALStatus *status, REAL8 *psd, REAL8 f) 
{ /* </lalVerbatim> */
   REAL8 s0, x, x2, x3, x4, x5, x6;
   REAL8 a1, a2, p1, p2, c1, c2, c3, c4, b1, b2, b3, b4, b5, b6;
   REAL8 num,den;

   status = NULL;
   x = f/200.;
  
   a1 = 185.62;
   a2 = 232.56;
   p1 = -4.05;
   p2 = -0.69;
   b1 = 31.18;
   b2 = -64.72;
   b3 = 52.24;
   b4 = -42.16;
   b5 = 10.17;
   b6 = 11.53;
   c1 = 13.58;
   c2 = -36.46;
   c3 = 18.56;
   c4 = 27.43;

   x2 = x * x;
   x3 = x2 * x;
   x4 = x2 * x2;
   x5 = x3 * x2;
   x6 = x3 * x3;
   
   num = 1 + b1*x + b2*x2 + b3*x3 + b4*x4 + x5*b5 + x6*b6;
   den = 1 + c1*x + c2*x2 + c3*x3 + c4*x4;
   
   /*new psds from fitted on the Design sensitivity curve from virgo web site*/
   s0 = 1.62e-51;
   *psd = s0*( pow(x,p1) +a1*pow(x,p2) +a2*num/den);

}
