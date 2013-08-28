/*
*  Copyright (C) 2007 Bernd Machenschalk, B.S. Sathyaprakash, Thomas Cokelaer
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
 * \brief Module to calculate the noise power spectral density for the initial LIGO detector.
 *
 * \heading{Description}
 *
 * This function module takes as an input a frequency \f$f\f$ in Hz, and it
 * calculates the noise spectral density (per Hz) \f$S_{h}(f)\f$
 * for that frequency. The noise PSD is based on data provided by
 * Kip Thorne, and the fit by B.S.Sathyaprakash
 *
 * \f{equation}{
 * S_h(f) = S_0\left\{  \left(\frac{f}{f_0}\right)^{-4.14} - 5\left(\frac{f_0}{f}\right)^2 + 111  \left(\frac{1. -
 * \frac{f}{f_0}^2 + 0.5  \frac{f}{f_0}^4}{1. + 0.5\frac{f}{f_0}^2} \right)\right\};
 * \f}
 * where, \f$f_0=215\f$Hz
 * The returned value is scaled up by \f$S_0 = 10^{49}.\f$
 *
 */
void
LALAdvLIGOPsd (LALStatus UNUSED *status, REAL8 *psd, REAL8 f)
{

   REAL8 x2,x;
   x = f/215.;
   x2 = x*x;
   *psd = pow(x,-4.14) - 5./x2 + 111. * (1. - x2 + 0.5 * x2*x2)/(1. + 0.5*x2);
}
