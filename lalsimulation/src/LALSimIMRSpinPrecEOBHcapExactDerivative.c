/** Static function file **/
#ifndef _LALSIMIMRSPINPRECEOBHCAPEXACTDERIVATIVE_C
#define _LALSIMIMRSPINPRECEOBHCAPEXACTDERIVATIVE_C

#include "LALSimIMRCalculateSpinPrecEOBHCoeffs.c" /* SEOBNRv3_opt */

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

static INT4 XLALSEOBNRv3_opt_HDerivs_for_Omega(
            const REAL8 * valuestort,
            const REAL8 mass1, const REAL8 mass2, const REAL8 eta,
            SpinEOBHCoeffs * coeffs,
            REAL8 * derivs,
            SpinEOBParams * params);

static INT4 XLALSEOBNRv3_opt_ComputeHamiltonianDerivatives(
            const REAL8 * valuestort1, const REAL8 * valuestort2,
            SpinEOBHCoeffs * coeffs,
            REAL8 * derivs,
            SpinEOBParams * params,
            REAL8 * HReal);

/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

static INT4 XLALSEOBNRv3_opt_HDerivs_for_Omega(
						      const REAL8 * valuestort,
						      const REAL8 mass1, const REAL8 mass2, const REAL8 eta,
						      SpinEOBHCoeffs * coeffs,
						      REAL8 * derivs,
						      SpinEOBParams * params){

  REAL8 sigmaKerrData[3];
  REAL8 sigmaStarData[3];
  REAL8 s1VecData[3];
  REAL8 s2VecData[3];

  /* Spins are independent of tortoise value */
  memcpy(s1VecData, valuestort + 6, 3 * sizeof(REAL8));
  memcpy(s2VecData, valuestort + 9, 3 * sizeof(REAL8));

  for(int i=0; i<3; i++){
    sigmaKerrData[i] = s1VecData[i]+s2VecData[i];
    sigmaStarData[i] = (mass2/mass1)*s1VecData[i]+(mass1/mass2)*s2VecData[i];
  }

  REAL8Vector dsigmaKerr;
  REAL8Vector dsigmaStar;
  REAL8Vector ds1Vec;
  REAL8Vector ds2Vec;

  dsigmaKerr.data = sigmaKerrData;
  dsigmaStar.data = sigmaStarData;
  ds1Vec.data = s1VecData;
  ds2Vec.data = s2VecData;

  REAL8Vector * sigmaKerr;
  REAL8Vector * sigmaStar;
  REAL8Vector * s1Vec;
  REAL8Vector * s2Vec;

  s1Vec = &ds1Vec;
  s2Vec = &ds2Vec;
  sigmaKerr = &dsigmaKerr;
  sigmaStar = &dsigmaStar;

  SEOBHCoeffConstants * seobCoeffConsts = params->seobCoeffConsts;

  if(coeffs->updateHCoeffs){
    SpinEOBHCoeffs tmpCoeffs;
    REAL8 tmpa;

    tmpa = sqrt(sigmaKerr->data[0]*sigmaKerr->data[0]
                + sigmaKerr->data[1]*sigmaKerr->data[1]
                + sigmaKerr->data[2]*sigmaKerr->data[2]);

    if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs( &tmpCoeffs, eta, tmpa, coeffs->SpinAlignedEOBversion ) == XLAL_FAILURE )
      {
	XLAL_ERROR( XLAL_EFUNC );
      }

    tmpCoeffs.SpinAlignedEOBversion = coeffs->SpinAlignedEOBversion;
    tmpCoeffs.updateHCoeffs = coeffs->updateHCoeffs;

    coeffs = &tmpCoeffs;
  }

  REAL8 a2_prederiv = sigmaKerr->data[0] * sigmaKerr->data[0] + sigmaKerr->data[1] * sigmaKerr->data[1] +  sigmaKerr->data[2] * sigmaKerr->data[2];
  REAL8 a_prederiv = sqrt(a2_prederiv);

  INT4 divby0 = 0;

  REAL8 e3_x, e3_y, e3_z;

  if(a_prederiv !=0.)
    {
      const REAL8 inva = 1./a_prederiv;
      e3_x = sigmaKerr->data[0] * inva;
      e3_y = sigmaKerr->data[1] * inva;
      e3_z = sigmaKerr->data[2] * inva;
    }
  else
    {
      /* SEOBNRv3_opt: Since spin=0, we are free to choose the "spin direction". */
      e3_x = 1./sqrt(3.);
      e3_y = 1./sqrt(3.);
      e3_z = 1./sqrt(3.);

      divby0=1;
    }

  REAL8Vector xVec, pVec;
  xVec.length = pVec.length= 3;

  REAL8 xData[3] = {0.}, pData[3] = {0.};
  xVec.data = xData;
  pVec.data = pData;

  REAL8Vector * x;
  REAL8Vector * p;
  x=&xVec;
  p=&pVec;

  memcpy(xVec.data,valuestort,3*sizeof(REAL8));
  memcpy(pVec.data,valuestort+3,3*sizeof(REAL8));

  const REAL8 invr = 1./sqrt(xData[0]*xData[0]+xData[1]*xData[1]+xData[2]*xData[2]);

  if (1. - fabs(e3_x*(xData[0]*invr) + e3_y*(xData[1]*invr) + e3_z*(xData[2]*invr)) <= 1.e-8) {
    e3_x = e3_x+0.1;
    e3_y = e3_y+0.1;
    const REAL8 invnorm = 1./sqrt(e3_x*e3_x + e3_y*e3_y + e3_z*e3_z);
    e3_x = e3_x*invnorm;
    e3_y = e3_y*invnorm;
    e3_z = e3_z*invnorm;
    divby0 = 1;
  }

  if(divby0) {
    /* s1 & s2Vec's cannot all be zero when taking spin derivatives, because naturally
       some s1Vec's & s2Vec's appear in denominators of exact deriv expressions.
    */
    const double epsilon_spin=1e-14;
    if(fabs(s1Vec->data[0] + s2Vec->data[0])<epsilon_spin &&
       fabs(s1Vec->data[1] + s2Vec->data[1])<epsilon_spin &&
       fabs(s1Vec->data[2] + s2Vec->data[2])<epsilon_spin) {
      s1Vec->data[0] = epsilon_spin;
      s1Vec->data[1] = epsilon_spin;
      s1Vec->data[2] = epsilon_spin;
      s2Vec->data[0] = epsilon_spin;
      s2Vec->data[1] = epsilon_spin;
      s2Vec->data[2] = epsilon_spin;
    }
  }
  const REAL8 etainv = 1.0/eta;

  /* SEOBNRv3_opt: Must define the seob coeffs constants for spin derivatives */
  UNUSED REAL8 c0k2 = seobCoeffConsts->a0k2;
  UNUSED REAL8 c1k2 = seobCoeffConsts->a1k2;
  UNUSED REAL8 c0k3 = seobCoeffConsts->a0k3;
  UNUSED REAL8 c1k3 = seobCoeffConsts->a1k3;
  UNUSED REAL8 c0k4 = seobCoeffConsts->a0k4;
  UNUSED REAL8 c1k4 = seobCoeffConsts->a1k4;
  UNUSED REAL8 c2k4 = seobCoeffConsts->a2k4;
  UNUSED REAL8 c0k5 = seobCoeffConsts->a0k5;
  UNUSED REAL8 c1k5 = seobCoeffConsts->a1k5;
  UNUSED REAL8 c2k5 = seobCoeffConsts->a2k5;

  #include "seobnrv3_opt_exactderivs/exact_derivatives-HrealHeader.c"

  /* Momentum derivatives */
  INT4 tortoise = 1;
  {
    #include "seobnrv3_opt_exactderivs/exact_derivatives-Hreal.c"
    {
      #include "seobnrv3_opt_exactderivs/exact_derivatives-px.c"
      derivs[0]=Hrealprm*etainv;}
    {
      #include "seobnrv3_opt_exactderivs/exact_derivatives-py.c"
      derivs[1]=Hrealprm*etainv;}
    {
      #include "seobnrv3_opt_exactderivs/exact_derivatives-pz.c"
      derivs[2]=Hrealprm*etainv;}
  }

  if(isnan(derivs[0]*derivs[1]*derivs[2])){
    XLALPrintError("XLALSEOBNRv3_opt_HDerivs_for_Omega failed!\n");
    XLALPrintError("NAN in derivative! derivs0,1,2 = %e %e %e | divby0 = %d\n",derivs[0],derivs[1],derivs[2],divby0);
    XLAL_ERROR( XLAL_EFAILED );
  }

  return 0;
}

UNUSED static INT4 XLALSEOBNRv3_opt_ComputeHamiltonianDerivatives(const REAL8 * valuestort1 , const REAL8 * valuestort2, SpinEOBHCoeffs * coeffs, REAL8 * derivs, SpinEOBParams * params,REAL8 * HReal){

  REAL8 mass1 = params->eobParams->m1;
  REAL8 mass2 = params->eobParams->m2;
  REAL8 eta = params->eobParams->eta;

  REAL8 sigmaKerrData[3];
  REAL8 sigmaStarData[3];
  REAL8 s1VecData[3];
  REAL8 s2VecData[3];

  REAL8Vector * sigmaKerr;
  REAL8Vector * sigmaStar;
  REAL8Vector * s1Vec;
  REAL8Vector * s2Vec;

  REAL8Vector dsigmaKerr;
  REAL8Vector dsigmaStar;
  REAL8Vector ds1Vec;
  REAL8Vector ds2Vec;

  s1Vec = &ds1Vec;
  s2Vec = &ds2Vec;
  sigmaKerr = &dsigmaKerr;
  sigmaStar = &dsigmaStar;

  /* Spins are independent of tortoise value */
  memcpy(s1VecData, valuestort1 + 6, 3 * sizeof(REAL8));
  memcpy(s2VecData, valuestort1 + 9, 3 * sizeof(REAL8));

  const REAL8 m1divm2 = mass1/mass2;
  const REAL8 m2divm1 = mass2/mass1;

  for(int i=0; i<3; i++){
    sigmaKerrData[i] = s1VecData[i]+s2VecData[i];
    sigmaStarData[i] = m2divm1*s1VecData[i]+m1divm2*s2VecData[i];
  }

  dsigmaKerr.data = sigmaKerrData;
  dsigmaStar.data = sigmaStarData;
  ds1Vec.data = s1VecData;
  ds2Vec.data = s2VecData;

  SEOBHCoeffConstants * seobCoeffConsts = params->seobCoeffConsts;

  REAL8 a2_prederiv = sigmaKerr->data[0] * sigmaKerr->data[0] + sigmaKerr->data[1] * sigmaKerr->data[1] +  sigmaKerr->data[2] * sigmaKerr->data[2];
  REAL8 a_prederiv = sqrt(a2_prederiv);

  if(coeffs->updateHCoeffs){
    SpinEOBHCoeffs tmpCoeffs;

    if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs( &tmpCoeffs, eta, a_prederiv, coeffs->SpinAlignedEOBversion ) == XLAL_FAILURE ){
      XLAL_ERROR( XLAL_EFUNC );
    }

    tmpCoeffs.SpinAlignedEOBversion = coeffs->SpinAlignedEOBversion;
    tmpCoeffs.updateHCoeffs = coeffs->updateHCoeffs;

    coeffs = &tmpCoeffs;
  }

  INT4 divby0 = 0;

  REAL8 e3_x, e3_y, e3_z, e3x, e3y, e3z;
  const REAL8 inva = 1./a_prederiv;
  const REAL8 inva3 = inva * inva * inva;

  if(a_prederiv !=0.){
    e3x = sigmaKerr->data[0] * inva;
    e3y = sigmaKerr->data[1] * inva;
    e3z = sigmaKerr->data[2] * inva;

    e3_x = e3x;
    e3_y = e3y;
    e3_z = e3z;
  }
  else{
    /* SEOBNRv3_opt: Since spin=0, we are free to choose the "spin direction". */
    e3x = 1./sqrt(3.);
    e3y = 1./sqrt(3.);
    e3z = 1./sqrt(3.);

    e3_x = e3x;
    e3_y = e3y;
    e3_z = e3z;

    divby0=1;
  }
  REAL8Vector xVec, pVec;
  xVec.length = pVec.length = 3;

  REAL8 xData[3] = {0.}, pData[3] = {0.};
  xVec.data = xData;
  pVec.data = pData;

  REAL8Vector * x;
  REAL8Vector * p;

  x=&xVec;
  p=&pVec;

  memcpy(xVec.data,valuestort1,3*sizeof(REAL8));
  memcpy(pVec.data,valuestort1+3,3*sizeof(REAL8));

  const REAL8 invr = 1./sqrt(xData[0]*xData[0]+xData[1]*xData[1]+xData[2]*xData[2]);

  REAL8 invnorm = 0.;
  if (1. - fabs(e3_x*(xData[0]*invr) + e3_y*(xData[1]*invr) + e3_z*(xData[2]*invr)) <= 1.e-8) {
    e3_x = e3_x+0.1;
    e3_y = e3_y+0.1;
    invnorm = 1./sqrt(e3_x*e3_x + e3_y*e3_y + e3_z*e3_z);
    e3_x = e3_x*invnorm;
    e3_y = e3_y*invnorm;
    e3_z = e3_z*invnorm;
    divby0 = 1;
  }

  const REAL8 invnorm3 = invnorm * invnorm * invnorm;

  if(divby0) {
    /* s1 & s2Vec's cannot all be zero when taking spin derivatives, because naturally
       some s1Vec's & s2Vec's appear in denominators of exact deriv expressions.
    */
    const double epsilon_spin=1e-14;
    if(fabs(s1Vec->data[0] + s2Vec->data[0])<epsilon_spin &&
       fabs(s1Vec->data[1] + s2Vec->data[1])<epsilon_spin &&
       fabs(s1Vec->data[2] + s2Vec->data[2])<epsilon_spin) {
      s1Vec->data[0] = epsilon_spin;
      s1Vec->data[1] = epsilon_spin;
      s1Vec->data[2] = epsilon_spin;
      s2Vec->data[0] = epsilon_spin;
      s2Vec->data[1] = epsilon_spin;
      s2Vec->data[2] = epsilon_spin;
    }
  }

  const REAL8 etainv = 1.0/eta;

  /* SEOBNRv3_opt: Must define the seob coeffs constants for spin derivatives */
  UNUSED REAL8 c0k2 = seobCoeffConsts->a0k2;
  UNUSED REAL8 c1k2 = seobCoeffConsts->a1k2;
  UNUSED REAL8 c0k3 = seobCoeffConsts->a0k3;
  UNUSED REAL8 c1k3 = seobCoeffConsts->a1k3;
  UNUSED REAL8 c0k4 = seobCoeffConsts->a0k4;
  UNUSED REAL8 c1k4 = seobCoeffConsts->a1k4;
  UNUSED REAL8 c2k4 = seobCoeffConsts->a2k4;
  UNUSED REAL8 c0k5 = seobCoeffConsts->a0k5;
  UNUSED REAL8 c1k5 = seobCoeffConsts->a1k5;
  UNUSED REAL8 c2k5 = seobCoeffConsts->a2k5;

#include "seobnrv3_opt_exactderivs/exact_derivatives-HrealHeader.c"

  /* Position derivatives */
  INT4 tortoise = 2;
  {
  #include "seobnrv3_opt_exactderivs/exact_derivatives-Hreal.c"
  {
  #include "seobnrv3_opt_exactderivs/exact_derivatives-x.c"
  derivs[0]=Hrealprm*etainv;}
  {
  #include "seobnrv3_opt_exactderivs/exact_derivatives-y.c"
  derivs[1]=Hrealprm*etainv;}
  {
  #include "seobnrv3_opt_exactderivs/exact_derivatives-z.c"
  derivs[2]=Hrealprm*etainv;}
  }

  memcpy(xVec.data,valuestort2,3*sizeof(REAL8));
  memcpy(pVec.data,valuestort2+3,3*sizeof(REAL8));

  /* Momentum and spin derivatives  */
 tortoise = 1;
  {
  #include "seobnrv3_opt_exactderivs/exact_derivatives-Hreal.c"
  *HReal = Hreal;
  {
  #include "seobnrv3_opt_exactderivs/exact_derivatives-px.c"
  derivs[3]=Hrealprm*etainv;}
  {
  #include "seobnrv3_opt_exactderivs/exact_derivatives-py.c"
  derivs[4]=Hrealprm*etainv;}
  {
  #include "seobnrv3_opt_exactderivs/exact_derivatives-pz.c"
  derivs[5]=Hrealprm*etainv;}
  {
    double e3USCORExprm = inva - sigmaKerr->data[0] * sigmaKerr->data[0] * inva3;
    double e3USCOREyprm =      - sigmaKerr->data[1] * sigmaKerr->data[0] * inva3;
    double e3USCOREzprm =      - sigmaKerr->data[2] * sigmaKerr->data[0] * inva3;

    if(divby0){
      double e3xprm = e3USCORExprm;
      double e3yprm = e3USCOREyprm;
      double e3zprm = e3USCOREzprm;

      e3USCORExprm = e3xprm * invnorm - (e3x + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREyprm = e3yprm * invnorm - (e3y + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREzprm = e3zprm * invnorm -  e3z        * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
    }

    #include "seobnrv3_opt_exactderivs/exact_derivatives-s1x.c"
    derivs[6]=Hrealprm*etainv;}
  {
    double e3USCORExprm =      - sigmaKerr->data[0] * sigmaKerr->data[1] * inva3;
    double e3USCOREyprm = inva - sigmaKerr->data[1] * sigmaKerr->data[1] * inva3;
    double e3USCOREzprm =      - sigmaKerr->data[2] * sigmaKerr->data[1] * inva3;

    if(divby0){
      double e3xprm = e3USCORExprm;
      double e3yprm = e3USCOREyprm;
      double e3zprm = e3USCOREzprm;

      e3USCORExprm = e3xprm * invnorm - (e3x + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREyprm = e3yprm * invnorm - (e3y + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREzprm = e3zprm * invnorm -  e3z        * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
    }

    #include "seobnrv3_opt_exactderivs/exact_derivatives-s1y.c"
    derivs[7]=Hrealprm*etainv;}
  {
    double e3USCORExprm =      - sigmaKerr->data[0] * sigmaKerr->data[2] * inva3;
    double e3USCOREyprm =      - sigmaKerr->data[1] * sigmaKerr->data[2] * inva3;
    double e3USCOREzprm = inva - sigmaKerr->data[2] * sigmaKerr->data[2] * inva3;

    if(divby0){
      double e3xprm = e3USCORExprm;
      double e3yprm = e3USCOREyprm;
      double e3zprm = e3USCOREzprm;

      e3USCORExprm = e3xprm * invnorm - (e3x + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREyprm = e3yprm * invnorm - (e3y + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREzprm = e3zprm * invnorm -  e3z        * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
    }

    #include "seobnrv3_opt_exactderivs/exact_derivatives-s1z.c"
    derivs[8]=Hrealprm*etainv;}
  {
    double e3USCORExprm = inva - sigmaKerr->data[0] * sigmaKerr->data[0] * inva3;
    double e3USCOREyprm =      - sigmaKerr->data[1] * sigmaKerr->data[0] * inva3;
    double e3USCOREzprm =      - sigmaKerr->data[2] * sigmaKerr->data[0] * inva3;

    if(divby0){
      double e3xprm = e3USCORExprm;
      double e3yprm = e3USCOREyprm;
      double e3zprm = e3USCOREzprm;

      e3USCORExprm = e3xprm * invnorm - (e3x + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREyprm = e3yprm * invnorm - (e3y + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREzprm = e3zprm * invnorm -  e3z        * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
    }
    #include "seobnrv3_opt_exactderivs/exact_derivatives-s2x.c"
    derivs[9]=Hrealprm*etainv;}
  {
    double e3USCORExprm =      - sigmaKerr->data[0] * sigmaKerr->data[1] * inva3;
    double e3USCOREyprm = inva - sigmaKerr->data[1] * sigmaKerr->data[1] * inva3;
    double e3USCOREzprm =      - sigmaKerr->data[2] * sigmaKerr->data[1] * inva3;

    if(divby0){
      double e3xprm = e3USCORExprm;
      double e3yprm = e3USCOREyprm;
      double e3zprm = e3USCOREzprm;

      e3USCORExprm = e3xprm * invnorm - (e3x + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREyprm = e3yprm * invnorm - (e3y + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREzprm = e3zprm * invnorm -  e3z        * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
    }
    #include "seobnrv3_opt_exactderivs/exact_derivatives-s2y.c"
    derivs[10]=Hrealprm*etainv;}
  {
    double e3USCORExprm =      - sigmaKerr->data[0] * sigmaKerr->data[2] * inva3;
    double e3USCOREyprm =      - sigmaKerr->data[1] * sigmaKerr->data[2] * inva3;
    double e3USCOREzprm = inva - sigmaKerr->data[2] * sigmaKerr->data[2] * inva3;

    if(divby0){
      double e3xprm = e3USCORExprm;
      double e3yprm = e3USCOREyprm;
      double e3zprm = e3USCOREzprm;

      e3USCORExprm = e3xprm * invnorm - (e3x + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREyprm = e3yprm * invnorm - (e3y + 0.1) * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
      e3USCOREzprm = e3zprm * invnorm -  e3z        * (((e3x + 0.1) * e3xprm + (e3y + 0.1) * e3yprm + e3z * e3zprm) * invnorm3);
    }
    #include "seobnrv3_opt_exactderivs/exact_derivatives-s2z.c"
    derivs[11]=Hrealprm*etainv;}
  }

  if(isnan(derivs[0]*derivs[1]*derivs[2]*derivs[3]*derivs[4]*derivs[5]*derivs[6]*derivs[7]*derivs[8]*derivs[9]*derivs[10]*derivs[11])){
    XLALPrintError("XLALSEOBNRv3_opt_ComputeHamiltonianDerivatives failed!\n");
    XLALPrintError("NAN in derivative! derivs0,1,2,3,4,5,6,7,8,9,10,11 = %e %e %e %e %e %e %e %e %e %e %e %e | divby0 = %d\n",derivs[0],derivs[1],derivs[2],derivs[3],derivs[4],derivs[5],derivs[6],derivs[7],derivs[8],derivs[9],derivs[10],derivs[11],divby0);
    XLAL_ERROR( XLAL_EFAILED );
  }

 return 0;
}

#endif // _LALSIMIMRSPINPRECEOBHCAPEXACTDERIVATIVE_C
