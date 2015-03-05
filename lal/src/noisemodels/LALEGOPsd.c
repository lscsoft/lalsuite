/*
*  Copyright (C) 2007 Bernd Machenschalk, Thomas Cokelaer
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
 * \author Cokelaer T.
 * \ingroup LALNoiseModels_h
 * \brief Function to calculate the noise power spectral density for the EGO detector.
 *
 * ### Description ###
 *
 * The module takes as an input a frequency \f$f\f$ in Hz, and it
 * calculates the noise spectral density (per Hz) \f$S_{h}(f)\f$
 * for that frequency. The noise PSD is based on data provided by
 * grqc/0607092
 * \f{equation}{
 * S_h(f) =
 * s_0 \left\{  x^{p_1} + a_1x^{p_2} +a_2 \frac{1+b_1x +b_2x^2+b_3x^3+b_4x^4+b_5x^5+b_6x^6}{1+c_1x+c_2x^2+c_3x^3+c_4x^4} \right\}
 * \f}
 * where \f$S_0=1.61e-51\f$\\
 * \f$p_1=-4.05, p_2=-0.69\f$\\
 * \f$a_1=185.62, a_2=232.56\f$\\
 * \f$b_1 = 31.18, b_2=-64.72, b_3=52.24, b_4=-42.16, b_5=10.17, b_6=11.53\f$\\
 * and \f$c_1=13.58, c_2 = -36.46, c_3=18.56, c_4=27.43\f$
 *
 */
void LALEGOPsd (LALStatus UNUSED *status, REAL8 *psd, REAL8 f)
{
   REAL8 s0, x, x2, x3, x4, x5, x6;
   REAL8 a1, a2, p1, p2, c1, c2, c3, c4, b1, b2, b3, b4, b5, b6;
   REAL8 num,den;

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
