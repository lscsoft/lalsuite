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
\idx{LALTAMAPsd()}

\subsubsection*{Description}

The module takes as an input a frequency $f$ in Hz, and it
calculates the noise spectral density (per Hz) $S_{h}(f)$
for that frequency. The noise PSD is based on data provided by
M.-K Fujimoto (see T. Damour, B.R. Iyer and B.S. Sathyaprakash,
Phys. Rev. D 63, 044023 (2001)) and is approximated by
the following:
\begin{equation}
   S_h(f) = \left ( \frac{f}{f_0} \right )^{-5} + 13 \frac{f_0}{f} +
   9 \left [1 + \left( \frac{f}{f_0} \right)^2 \right ].
\end{equation}
The returned value is scaled up by $s_0 = 10^{46}/75.$ In otherwords,
the expected noise PSD is $75 \times 10^{-46}$ times the returned value.
\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALTAMAPsdCV}}

</lalLaTeX>  */



#include <lal/LALNoiseModels.h>

NRCSID (LALTAMAPSDC,"$Id$");

/*  <lalVerbatim file="LALTAMAPsdCP"> */
void  LALTAMAPsd(LALStatus *status, REAL8 *psd, REAL8 f)
{ /* </lalVerbatim> */

   REAL8 seismic, thermal, shot, x;

   status = NULL;
   x = f/400.;
   seismic = pow(x,-5);
   thermal = 13. / x;
   shot = 9. * (1. + x*x);
   *psd = seismic + thermal + shot;
}
