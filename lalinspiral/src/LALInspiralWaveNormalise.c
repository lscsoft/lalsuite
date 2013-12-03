/*
*  Copyright (C) 2007 David Churches, B.S. Sathyaprakash
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
 * \author Sathyaprakash, B. S.
 * \file
 *
 * \brief Module to find the norm of a signal and to return a normalised
 * array. The original signal is left untouched.
 *
 * ### Prototypes ###
 *
 * <tt>LALInspiralWaveNormalise()</tt>
 *
 * ### Description ###
 *
 * Given the positive frequency Fourier components
 * \f$H_k,\f$ \f$k=0,\ldots,n-1,\f$ of a vector
 * and the noise PSD \f$S_m,\f$ \f$m=0,\ldots,n/2,\f$
 * this module first computes the norm \f$H\f$ of the vector treating
 * \f$S_m\f$ as the measure:
 * (note that in {\em fftw} notation, the zeroth frequency
 * component is \f$H_0,\f$ Nyquist
 * is \f$H_{n/2},\f$ \f$H_k,\f$ \f$k \ne 0,n/2,\f$ (\f$H_{n-k})\f$ is the real (imaginary)
 * part of the \f$k\f$th harmonic)
 * \f{equation}{
 * \label{eq_inspiralnorm}
 * H = \sum_{k=1}^{n/2-1} \frac{H_k^2 + H^2_{n-k}}{S_k}.
 * \f}
 * (Note that the zeroth and Nyquist components are ignored in the
 * computation of the norm.)
 * It then replaces the original vector \f$H_k\f$ with normalized
 * vector using:
 * \f{equation}{
 * \widehat H_k = \frac {H_k}{\sqrt H},\ \ k=0,\ldots n-1.
 * \f}
 *
 * ### Algorithm ###
 *
 *
 * ### Uses ###
 *
 * \code
 * none.
 * \endcode
 *
 * ### Notes ###
 *
 */
#include <lal/LALNoiseModelsInspiral.h>

void
LALInspiralWaveNormalise
   (
   LALStatus   *status,
   REAL4Vector *in,
   REAL8       *norm,
   REAL8Vector psd
   )
{

  INT4 i, n, nby2, k;
  REAL8 psdvalue;

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT (in->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (psd.length == in->length/2+1, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

  n = in->length;
  *norm = 0;
  nby2 = n/2;

  for (i=1; i<nby2; i++)
  {
     k = n-i;
     psdvalue = psd.data[i];
     if (psdvalue)
     {
	*norm += (in->data[i]*in->data[i] + in->data[k]*in->data[k])/psdvalue;
     }
  }

  /* If we want to add the negative frequency components as well now is
   * the time to do that before 0th and Nyquist contributions are added
   */

  *norm *= 2.;

  /*
  if (psd.data[nby2]) *norm += pow(in->data[nby2], 2.)/psd.data[nby2];
  if (psd.data[0]) *norm += pow(in->data[0], 2.)/psd.data[0];
  */
  /* Set the 0th and Nyquist frequency bins to be zero. */
  in->data[0] = in->data[nby2] = 0.;

  *norm = sqrt(*norm);

  for (i=0; i<n; i++)
  {
	  *(in->data+i) /= *norm;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
