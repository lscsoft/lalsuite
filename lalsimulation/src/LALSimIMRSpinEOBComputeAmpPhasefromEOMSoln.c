/*
*  Copyright (C) 2011 Craig Robinson, Enrico Barausse, Yi Pan
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
 * \author Craig Robinson, Yi Pan
 *
 * Functions for calculating the effective one-body Hamiltonian for
 * spinning binaries, as described in
 * Taracchini et al. ( PRD 86, 024011 (2012), arXiv 1202.0790 ).
 * All equation numbers in this file refer to equations of this paper,
 * unless otherwise specified.
 * This code borrows hugely from a C implementation originally written
 * by Enrico Barausse, following Barausse and Buonanno
 * PRD 81, 084024 (2010) and PRD 84, 104027 (2011), henceforth BB1 and BB2
 */

#ifndef _LALSIMIMRSPINEOBCOMPUTEAMPPHASEFROMEOMSOLN_C
#define _LALSIMIMRSPINEOBCOMPUTEAMPPHASEFROMEOMSOLN_C

#include <stdio.h>
#include <math.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */


static double GenericSpinAlignedHamiltonianWrapperOptimized( REAL8* values, SpinEOBParams *params );
static void GenerateAmpPhaseFromEOMSoln(UINT4 retLen,REAL8 *inputdata,SpinEOBParams *seobParams);

/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

/* This function simply reads in the four inputs from values[] and outputs the Optimized SEOBNRv2 Hamiltonian */
static double GenericSpinAlignedHamiltonianWrapperOptimized( REAL8* values, SpinEOBParams *params )
{
  EOBParams *eobParams = params->eobParams;
  REAL8 tmpVec[6];

  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p;
  REAL8Vector *s1Vec = params->s1Vec;
  REAL8Vector *s2Vec = params->s2Vec;
  REAL8Vector *sigmaKerr = params->sigmaKerr;
  REAL8Vector *sigmaStar = params->sigmaStar;

  /* Use a temporary vector to avoid corrupting the main function */
  memcpy(tmpVec, values, sizeof(tmpVec));

  /* Set the LAL-style vectors to point to the appropriate things */
  r.length   = p.length = 3;
  r.data     = tmpVec;
  p.data     = tmpVec+3;

  /* Notice that the Hamiltonian is not divided by eta. */
  REAL8 xlalham = XLALSimIMRSpinEOBHamiltonianOptimized( eobParams->eta, &r, &p, s1Vec, s2Vec, sigmaKerr, sigmaStar, params->tortoise, params->seobCoeffs );
  return xlalham;
}

static void GenerateAmpPhaseFromEOMSoln(UINT4 retLen,REAL8 *inputdata,SpinEOBParams *seobParams) {

  for(UINT4 jj=0;jj<retLen;jj++) {
    REAL8 yForOmegaCirc[4],cartValues[6],yy[4];
    for(INT4 kk=1;kk<=4;kk++) yy[kk-1] = inputdata[kk*retLen+jj];
    yForOmegaCirc[0] = yy[0];
    yForOmegaCirc[1] = yy[1];
    yForOmegaCirc[2] = 0.0;
    yForOmegaCirc[3] = yy[3];

    cartValues[0] = yy[0];
    cartValues[1] = 0.0;
    cartValues[2] = 0.0;
    cartValues[3] = yy[2];
    cartValues[4] = yy[3]/yy[0];
    cartValues[5] = 0.0;

    REAL8 ham       = GenericSpinAlignedHamiltonianWrapperOptimized(cartValues,    seobParams);
    REAL8 omega     = XLALSimIMRSpinAlignedEOBCalcOmegaOptimized   (yy,            seobParams);
    REAL8 omegaCirc = XLALSimIMRSpinAlignedEOBCalcOmegaOptimized   (yForOmegaCirc, seobParams);

    //CALEBS: Consolidated variable declarations:
    /* REAL8 rOmegaSq,rOmega,sqrtR,v,v2,vh,vh3,eulerlogxabs,deltalm,rholm,rholmPwrl,hlm0,hlmPhase,hlmMag,nqcMag,nqcPhase; */

    /* REAL8 legendre; */
    /* REAL8 vPhi2; */


    EOBNonQCCoeffs *nqcCoeffs = seobParams->nqcCoeffs;
    FacWaveformCoeffs *hCoeffs = seobParams->eobParams->hCoeffs;
    REAL8 eta = seobParams->eobParams->eta;
    REAL8 r = yy[0];
    REAL8 phi = yy[1];
    REAL8 p = yy[2];
    //REAL8 pp = yy[3];

    /* OPTIMIZATION COMMENT: BEGIN: STEP 7 */

    /* Calculate the Tail term, 3rd term in Eq. 17, given by Eq. A6 */
    // OPTIMIZATION COMMENT: Since l=m=2
    REAL8 k = 2 * omega; /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform */
    REAL8 hathatk = ham * k; /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform */

    //	XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
    //	XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( 3.0, -2.0*hathatk, &lnr1, &arg1 ) );
    //if(zr <= 0.5) //OPTIMIZATION COMMENT: zr = 3.0 so, proceeding to else...
    //else

    gsl_sf_result lnr1, arg1;
    INT4 UNUSED status; /* UNUSED is needed to avoid compiler warning about a variable being declared but unused. */
    XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( 3.0, -2.0*hathatk, &lnr1, &arg1 ) ); /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform */
    //	XLAL_CALLGSL( status = gsl_sf_fact_e( l, &z2 ) ); //OPTIMIZATION COMMENT: Since l=2, l! = 2 = z2

    REAL8 TlmPhase = arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)); /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform */
    REAL8 TlmMag = exp(lnr1.val + LAL_PI * hathatk)*0.5; /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform */

    /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5 */
    //if ( ( (l+m)%2 ) == 0) //OPTIMIZATION COMMENT: Since l=m=2...
    REAL8 Slm = (ham*ham - 1.)/(2.*eta) + 1.; /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform*/

    REAL8 v = cbrt(omega);          /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform*/
    //pr	= values->data[2];
    REAL8 v2 = v * v;               /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform*/
    //omega   = v2 * v;
    REAL8 vh3     = ham * omega;    /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform*/
    REAL8 vh = cbrt(vh3);           /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform*/

    //eulerlogxabs = LAL_GAMMA + log( 2.0 * (REAL8)m * v ); //OPTIMIZATION COMMENT: m=2
    REAL8 eulerlogxabs = LAL_GAMMA + log( 4.0 * v ); /*Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform*/

    REAL8 vPhi2 = (omega*omega/omegaCirc)*cbrt(1.0/omegaCirc); /*Location of this line of code:1,2:XLALSimIMRSpinEOBGetSpinFactorizedWaveform,XLALSimIMRSpinAlignedEOBNonKeplerCoeff*/

    REAL8 legendre = 0.75*sqrt( 5.0 / (6*LAL_PI)); /*Location of this line of code: XLALScalarSphHarmThetaPiBy2,XLALAssociatedLegendreXIsZero */

    REAL8 hlm0 = vPhi2*seobParams->eobParams->prefixes->values[2][2];/*Location of this line of code: XLALSimIMRSpinEOBCalculateNewtonianMultipole */


    /* Calculate the non-Keplerian velocity */
    //	  vPhi = XLALSimIMRSpinAlignedEOBNonKeplerCoeff( values->data, &seobParams ); //OPTIMIZATION COMMENT: Replaced with its output

    /* Calculate the residue phase and amplitude terms */
    /* deltalm is the 4th term in Eq. 17, delta 22 given by Eq. A15, others  */
    /* rholm is the 5th term in Eq. 17, given by Eqs. A8 - A14 */
    /* auxflm is a special part of the 5th term in Eq. 17, given by Eq. A15 */
    /* Actual values of the coefficients are defined in the next function of this file */

    //OPTIMIZATION COMMENT: Since l=m=2...
    REAL8 deltalm = vh3*(hCoeffs->delta22vh3 + vh3*(hCoeffs->delta22vh6 /* Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform */
                                              + vh*vh*(hCoeffs->delta22vh9*vh)))
      + hCoeffs->delta22v5 *v*v2*v2 + hCoeffs->delta22v6 *v2*v2*v2 + hCoeffs->delta22v8 *v2*v2*v2*v2;
    REAL8 rholm	= 1. + v2*(hCoeffs->rho22v2 + v*(hCoeffs->rho22v3
                                                 + v*(hCoeffs->rho22v4
                                                      + v*(hCoeffs->rho22v5 + v*(hCoeffs->rho22v6
                                                                                 + hCoeffs->rho22v6l*eulerlogxabs + v*(hCoeffs->rho22v7
                                                                                                                       + v*(hCoeffs->rho22v8 + hCoeffs->rho22v8l*eulerlogxabs
                                                                                                                            + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs)*v2)))))));
    /* Raise rholm to the lth power */ //OPTIMIZATION COMMENT: since l = 2...
    REAL8 rholmPwrl = rholm*rholm; /* Location of this line of code: XLALSimIMRSpinEOBGetSpinFactorizedWaveform */

    /* In the equal-mass odd m case, there is no contribution from nonspin terms,
     * and the only contribution comes from the auxflm term that is proportional to chiA (asymmetric spins).
     * In this case, we must ignore the nonspin terms directly, since the leading term defined by
     * CalculateThisMultipolePrefix in LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
     */

    //OPTIMIZATION COMMENT: Since m = 2, we can comment out:
    /*if (eta == 0.25 && m % 2)
      {
      rholmPwrl = auxflm;
      }
      else
      {*/
    //rholmPwrl += auxflm; /OPTIMIZATION COMMENT: ignoring because auxflm = 0.0;
    //} //OPTIMIZATION COMMENT: m=2.
    /*if (r > 8.5)
      {
      printf("YP::dynamics variables in waveform: %i, %i, %e, %e, %e\n",l,m,r,pp,values->data[1]);
      printf( "rholm^l = %.16e, Tlm = %.16e + i %.16e, |Tlm| = %.16e \nSlm = %.16e, hNewton = %.16e + i %.16e, delta = %.16e\n", rholmPwrl, creal(Tlm), cimag(Tlm), cabs(Tlm), Slm, creal(hNewton), cimag(hNewton), deltalm );
      printf("delta22 coeffs: vh3C = %.16e, vh6C = %.16e, vh9C = %.16e, v5C = %.16e, v6C = %.16e, v8C = %.16e\n",hCoeffs->delta22vh3,hCoeffs->delta22vh6,hCoeffs->delta22vh9,hCoeffs->delta22v5,hCoeffs->delta22v6,hCoeffs->delta22v8);
      printf("hNewt amp = %.16e, arg = %.16e\n",cabs(hNewton),carg(hNewton));
      }*/

    /*
      if (r > 8.5)
      {
      printf("YP::FullWave: Reh = %.16e, Imh = %.16e, hAmp = %.16e, hPhi = %.16e\n",creal(*hLM),cimag(*hLM),cabs(*hLM),carg(*hLM));
      }
    */

    REAL8 sqrtR = sqrt(r); /*Location of this line of code: XLALSimIMREOBNonQCCorrection*/
    REAL8 rOmega = r * omega; /*Location of this line of code: XLALSimIMREOBNonQCCorrection*/
    REAL8 rOmegaSq = rOmega*rOmega; /*Location of this line of code: XLALSimIMREOBNonQCCorrection*/

    /*printf("a1 = %.16e, a2 = %.16e, a3 = %.16e, a3S = %.16e, a4 = %.16e, a5 = %.16e\n",coeffs->a1,coeffs->a2,coeffs->a3,coeffs->a3S, coeffs->a4,coeffs->a5);
      printf("b1 = %.16e, b2 = %.16e, b3 = %.16e, b4 = %.16e\n",coeffs->b1,coeffs->b2,coeffs->b3,coeffs->b4);*/
    /* In EOBNRv2, coeffs->a3S, coeffs->a4 and coeffs->a5 are set to zero */
    /* through XLALSimIMREOBGetCalibratedNQCCoeffs() */
    /* and XLALSimIMREOBCalculateNQCCoefficients() */

    REAL8 nqcMag= 1. + (p*p / rOmegaSq) * ( nqcCoeffs->a1 /*Location of this line of code: XLALSimIMREOBNonQCCorrection*/
                                      + nqcCoeffs->a2 / r + ( nqcCoeffs->a3 + nqcCoeffs->a3S) / (r*sqrtR)
                                      + nqcCoeffs->a4 / (r*r) + nqcCoeffs->a5 / (r*r*sqrtR));

    REAL8 nqcPhase = nqcCoeffs->b1 * p / rOmega + p*p*p/rOmega * ( nqcCoeffs->b2 + nqcCoeffs->b3 / sqrtR + nqcCoeffs->b4 / r ); /*Location of this line of code: XLALSimIMREOBNonQCCorrection*/

    REAL8 hlmMag = hlm0 * legendre ; /*Location of this line of code: XLALScalarSphHarmThetaPiBy2 */
    REAL8 hlmPhase = -2.0*phi; /*Location of this line of code: XLALScalarSphHarmThetaPiBy2 */

    // *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl; //OPTIMIZATION COMMENT: From XLALSimIMRSpinEOBGetSpinFactorizedWaveform
    // OPTIMIZATION COMMENT: (the final) hlm = hlm*hNQC, (from XLALSimIMRSpinEOBGetSpinFactorizedWaveform) so, combine to get:

    /* Amplitude: */
    inputdata[(4+1)*retLen + jj] = TlmMag * hlmMag * (Slm*rholmPwrl) * nqcMag;

    /* Phase: */
    inputdata[(4+2)*retLen + jj] = (hlmPhase+TlmPhase+nqcPhase+deltalm);

  }
}


#endif /* _LALSIMIMRSPINEOBCOMPUTEAMPPHASEFROMEOMSOLN_C */
