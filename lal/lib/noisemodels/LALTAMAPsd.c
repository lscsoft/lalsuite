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


#include <lal/LALNoiseModels.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * \author Sathyaprakash, B. S.
 * \ingroup LALNoiseModels_h
 * \brief Function to calculate the noise power spectral density for the TAMA detector.
 *
 * ### Description ###
 *
 * The module takes as an input a frequency \f$f\f$ in Hz, and it
 * calculates the noise spectral density (per Hz) \f$S_{h}(f)\f$
 * for that frequency. The noise PSD is based on data provided by
 * M.-K Fujimoto (see T. Damour, B.R. Iyer and B.S. Sathyaprakash,
 * Phys. Rev. D 63, 044023 (2001)) and is approximated by
 * the following:
 * \f{equation}{
 * S_h(f) = \left ( \frac{f}{f_0} \right )^{-5} + 13 \frac{f_0}{f} +
 * 9 \left [1 + \left( \frac{f}{f_0} \right)^2 \right ].
 * \f}
 * The returned value is scaled up by \f$s_0 = 10^{46}/75.\f$ In otherwords,
 * the expected noise PSD is \f$75 \times 10^{-46}\f$ times the returned value.
 */
void  LALTAMAPsd(LALStatus UNUSED *status, REAL8 *psd, REAL8 f)
{

   REAL8 seismic, thermal, shot, x;

   x = f/400.;
   seismic = pow(x,-5);
   thermal = 13. / x;
   shot = 9. * (1. + x*x);
   *psd = seismic + thermal + shot;
}
