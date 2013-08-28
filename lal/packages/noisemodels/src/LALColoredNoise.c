/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash
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

#include <lal/LALNoiseModels.h>

/**
 * \author Sathyaprakash, B. S.
 * \ingroup LALNoiseModels_h
 * \brief This function colors a given white noise input into a colored noise
 * of power spectral density \c psd.
 *
 * \heading{Description}
 * Given the Fourier transform \f$N(f)\f$ of  white noise, the
 * Fourier transform of noise of power spectral density \f$S(f)\f$ is
 * given by \f${\cal N}(f) = N(f) \times \sqrt{S(f)}.\f$
 * In the discrete version there is an additional normalisation:
 * \f[{\cal N}_k = N_k \times \sqrt{\frac{2 S_k}{n}},\ \
 * {\cal N}_{n-k} = N_{n-k} \times \sqrt{\frac{2 S_k}{n}},\ \
 * k=1, \ldots, \frac{n}{2}.\f]
 *
 */
void
LALColoredNoise
   (
   LALStatus   *status,
   REAL4Vector *noisy,
   REAL8Vector  psd
   )
{

   INT4 i, j, n, nby2;
   REAL8 x, length;


   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   ASSERT (noisy,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (noisy->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (psd.length == noisy->length/2+1, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

   n = length = noisy->length;
   nby2 = n/2;

   for (i=1; i<nby2; i++)
   {
      j = n-i;
      /* Since fftw requires n and NOT n/2, I presume we
         don't need the factor 2 in the normalisation
         x = sqrt(2. * psd.data[i] / length);
      */
      x = sqrt(2.*psd.data[i]);
      noisy->data[i] *= x;
      noisy->data[j] *= x;
   }
   x = sqrt(2.*psd.data[0]);
   noisy->data[0] *= x;
   x = sqrt(2.*psd.data[nby2]);
   noisy->data[nby2] *= x;

   DETATCHSTATUSPTR(status);
   RETURN (status);
}
