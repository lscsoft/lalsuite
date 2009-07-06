/*
*  Copyright (C) 2007 Bernd Machenschalk, David Churches, Duncan Brown, Jolien Creighton
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

/*  <lalVerbatim file="LALGEOPsdCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALGEOPsd.c}}

Module to calculate the expected
noise power spectral density for the GEO600 detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALGEOPsdCP}
\idx{LALGEOPsd()}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in Hz, and it
calculates the noise spectral density (per Hz) $S_{h}(f)$
for that frequency. The noise PSD is based on data provided by
J. Hough and G. Cagnoli (see T. Damour, B.R. Iyer and B.S. Sathyaprakash,
Phys. Rev. D 63, 044023 (2001)) and is approximated by
the following:
\begin{equation}
   S_h(f) = 10^{-16} \left ( \frac{f}{f_0} \right)^{-30} +
            34 \frac{f_0 }{ f } +
   \frac{20 \left [1 - (f/f_0)^2 + 0.5 (f/f_0)^4 \right ] }{ 1 + 0.5 (f/f_0)^2}
\end{equation}
The returned value is scaled up by $s_0 = 10^{46}.$ In otherwords,
the expected noise PSD is a factor $10^{46}$ lower.
\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALGEOPsdCV}}

</lalLaTeX>  */


#include <lal/LALNoiseModels.h>

NRCSID (LALGEOPSDC,"$Id$");

/* <lalVerbatim file="LALGEOPsdCP"> */
void
LALGEOPsd(LALStatus *status, REAL8 *psd, REAL8 f)
{ /* </lalVerbatim> */

   REAL8 x, seismic, thermal, shot;

   status = NULL;
   x = f/150.;
   seismic = pow(10.,-16.) *  pow(x,-30.);
   thermal = 34. / x;
   shot = 20. * (1 - pow(x,2.) + 0.5 * pow(x,4.)) / (1. + 0.5 * pow(x,2.));
   *psd = seismic + thermal + shot;
}
