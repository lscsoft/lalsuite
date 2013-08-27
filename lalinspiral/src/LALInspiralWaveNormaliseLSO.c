/*
*  Copyright (C) 2007 B.S. Sathyaprakash
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

/**
\author Sathyaprakash, B. S.
\file

\brief Module to find the norm of a signal and to return a normaliseLSOd
array. The original signal is left untouched.

\heading{Prototypes}

<tt>LALInspiralWaveNormaliseLSO()</tt>

\heading{Description}
Given the positive frequency Fourier components
\f$H_k,\f$ \f$k=0,\ldots,n-1,\f$ of a vector
and the noise PSD \f$S_m,\f$ \f$m=0,\ldots,n/2,\f$
this module first computes the norm \f$H\f$ of the vector treating
\f$S_m\f$ as the measure:
(note that in {\em fftw} notation, the zeroth frequency
component is \f$H_0,\f$ Nyquist
is \f$H_{n/2},\f$ \f$H_k,\f$ \f$k \ne 0,n/2,\f$ (\f$H_{n-k})\f$ is the real (imaginary)
part of the \f$k\f$th harmonic)
\anchor eq_inspiralnorm \f{equation}{
H = \sum_{k=1}^{n/2-1} \frac{H_k^2 + H^2_{n-k}}{S_k}.
\tag{eq_inspiralnorm}
\f}
<tt>The above sum is limited to frequency</tt> <tt>in->fCutoff.</tt>
Also, note that the zeroth and Nyquist frequency components
are ignored in the computation of the norm.
Moreover, <tt>array elements of</tt> \c filter corresponding
to frequencies greater than <tt>in->fCutoff</tt> are <tt>set to zero</tt>.
That is, the code replaces the original vector \f$H_k\f$ with <em>normalized
vector</em> using:
\f{eqnarray}{
\widehat H_k & = & \frac {H_k}{\sqrt H},
\ \ \ \mathrm{k \times in\rightarrow\mathrm df} \le \mathrm{in\rightarrow fCutoff},\nonumber \\
& = & 0, \ \ \ \mathrm{k \times in\rightarrow df} > \mathrm{in\rightarrow fCutoff}.
\f}
In addition, the 0th and Nyquist frequency components are also set to zero.
\heading{Algorithm}
\heading{Uses}
\code
none.
\endcode

\heading{Notes}


*/
#include <lal/LALNoiseModelsInspiral.h>

void
LALInspiralWaveNormaliseLSO
   (
   LALStatus   *status,
   REAL4Vector *filter,
   REAL8       *norm,
   InspiralWaveNormaliseIn *in
   )
{

  INT4 i, n, nby2;
  REAL8 psd, f;

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT (filter->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (in->psd->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (in->psd->length == filter->length/2+1, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

  n = filter->length;
  *norm = 0;
  nby2 = n/2;

  for (i=1; i<nby2; i++)
  {
	  f = i*in->df;
	  if (f>in->fCutoff)
	  {
		  /* Since the normalisation is done using power only up to fCutoff
		   * it is better to terminate the frequency-domain waveform beyond
		   * fCutoff so that the signal doesn't contribute to the correlation
		   * outside this frequency band
		   */
		  filter->data[i] = filter->data[n-i] = 0.;
	  }
	  else
	  {
		  psd = in->psd->data[i];
		  if (psd)
		  {
			  *norm += (filter->data[i]*filter->data[i] + filter->data[n-i]*filter->data[n-i])/(psd*0.5);
		  }
	  }
  }

  /* If we want to add the negative frequency components as well now is
   * the time to do that before 0th and Nyquist contributions are added
   */

  (*norm) *= 2.;

  /* Set the 0th and Nyquist frequency bins to be zero. */
  filter->data[0] = filter->data[nby2] = 0.;

  *norm /= ((double) n * in->samplingRate);
  *norm = sqrt(*norm);

  for (i=0; i<n; i++)
  {
	  *(filter->data+i) /= *norm;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

