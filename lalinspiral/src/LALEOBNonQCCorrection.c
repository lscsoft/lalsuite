/*
*  Copyright (C) 2010 Craig Robinson 
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
 * \author Craig Robinson
 *
 * \brief Function to compute the factorized flux as uses in the new EOBNR_PP
 * model. Flux function given by Phys.Rev.D79:064004,2009.
 */

#include <complex.h>

#include <lal/LALComplex.h>
#include <lal/LALInspiral.h>

int  XLALEOBNonQCCorrection(
                      COMPLEX16             *nqc,
                      REAL8Vector           *values,
                      REAL8Vector           *dvalues,
                      EOBNonQCCoeffs        *coeffs
                     )

{

  REAL8 rOmega, rOmegaSq;
  REAL8 r, p;

  REAL8 mag, phase;


  r = values->data[0];
  p = values->data[1];

  rOmega = r * dvalues->data[1];
  rOmegaSq = rOmega*rOmega;

  mag = 1. + (p*p / rOmegaSq) * ( coeffs->a1
     + coeffs->a2 / r + coeffs->a3 / (r*sqrt(r))
     + coeffs->a4 / (r*r) );

  phase = coeffs->b1 * p / rOmega + coeffs->b2 * p*p*p/rOmega;

  nqc->re = mag * cos(phase);
  nqc->im = mag * sin(phase);

  return XLAL_SUCCESS;

}

