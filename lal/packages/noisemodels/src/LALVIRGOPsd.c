/*
*  Copyright (C) 2007 Bernd Machenschalk, David Churches, Duncan Brown, Jolien Creighton, Thomas Cokelaer
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
 * \author Sathyaprakash, B. S., Cokelaer T.
 * \ingroup LALNoiseModels_h
 * \brief Module to calculate the noise power spectral density for the VIRGO detector.
 *
 * \heading{Description}
 *
 * The module takes as an input a frequency \f$f\f$ in Hz, and it
 * calculates the noise spectral density (per Hz) \f$S_{h}(f)\f$
 * for that frequency. The noise PSD is based on data provided by
 * J-Y. Vinet and is approximated by
 * the following:
 * \f{equation}{
 * S_h(f) =
 * s_0 \left ( \frac {7.87f}{f_0} \right )^{-4.8} + \frac{6}{17} \frac{f_0}{f}
 * + \left [1 + \left (\frac {f}{f_0} \right)^2 \right ],
 * \f}
 * where \f$s_0=10.2e-46\f$
 *
 */
void LALVIRGOPsd (LALStatus UNUSED *status, REAL8 *psd, REAL8 f)
{
   REAL8 s0, x;

   x = f/500.;
/*
 * s1 = 34.6;
   s2 = 6.60;
   s3 = 3.24;
   *psd = pow(6.23*x,-5.) + 2.04/x + 1. + x*x;
   */

   /*new psds from fitted on the Design sensitivity curve from virgo web site*/
   s0 = 10.2e-46;
   *psd = s0*( pow(7.87*x,-4.8) + 6./17./x + 1. + x*x);

}
