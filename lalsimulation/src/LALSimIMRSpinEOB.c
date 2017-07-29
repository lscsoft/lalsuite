/*
*  Copyright (C) 2011 Craig Robinson, Enrico Barausse
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
 * \ file
 *
 * \brief Functions for producing EOB waveforms for
 * spinning binaries, as described in Barausse and Buonanno ( arXiv 0912.3517 ).
 */


#include <math.h>
#include <complex.h>
#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
#include <lal/SphericalHarmonics.h>
#include <gsl/gsl_sf_gamma.h>

#include "LALSimIMRSpinEOB.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

SEOBHCoeffConstants XLALEOBSpinPrecCalcSEOBHCoeffConstants(REAL8 eta){

    SEOBHCoeffConstants CoeffConsts;

    REAL8 c20  = 1.712;
    REAL8 c21  = -1.803949138004582;
    REAL8 c22  = -39.77229225266885;
    REAL8 c23  = 103.16588921239249;

    /*
    coeffsb3  = 0.;
    coeffsbb3 = 0.;
    */

    REAL8 coeffsKK =  c20 + c21*eta + c22*eta*eta + c23*eta*eta*eta;
    REAL8 m1PlusEtaKK = -1. + eta*coeffsKK;

    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_coeffs-kC-parsedfinal.h"
    CoeffConsts.a0k2 = kC0;
    CoeffConsts.a1k2 = kC1;
    CoeffConsts.a0k3 = kC2;
    CoeffConsts.a1k3 = kC3;
    CoeffConsts.a0k4 = kC4;
    CoeffConsts.a1k4 = kC5;
    CoeffConsts.a2k4 = kC6;
    CoeffConsts.a0k5 = kC7;
    CoeffConsts.a1k5 = kC8;
    CoeffConsts.a2k5 = kC9;

    return CoeffConsts;
}
