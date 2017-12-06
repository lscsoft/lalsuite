#ifndef _LALSPINPRECHCAPRVECDERIVATIVE_V3OPT_C
#define _LALSPINPRECHCAPRVECDERIVATIVE_V3OPT_C

#include <stdio.h>
#include <math.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficientsPrec.c"

/**
 * Function to calculate numerical derivatives of the spin EOB Hamiltonian,
 * to get \f$dr/dt\f$, as decribed in Eqs. A4 of PRD 81, 084041 (2010)
 * This function is not used by the spin-aligned SEOBNRv1 model.
 */
static int XLALSpinPrecHcapRvecDerivative_exact(
                 double UNUSED     t,         /**<< UNUSED */
                 const  REAL8      values[],  /**<< Dynamical variables */
                 REAL8             dvalues[], /**<< Time derivatives of variables (returned) */
                 void             *funcParams /**<< EOB parameters */
                               )
{
  //UNUSED int debugPK = 1;
  //if (debugPK){
    for(int i =0; i < 12; i++){
      if( isnan(values[i]) ) {
        XLAL_PRINT_INFO("XLALSpinPrecHcapRvecDerivative_exact::values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]);
          XLALPrintError( "XLAL Error - %s: nan in input values \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
        }
    }
  //}
    SpinEOBParams * seobParams = (SpinEOBParams *)funcParams;

    REAL8 tmpDValues[12] = {0.}, Tmatrix[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    seobParams->seobCoeffs->updateHCoeffs = 1;

    /* OPTv3: First, find csi and rMag2 for Tmatrix */
    REAL8 csi=1.0, rMag2;
    rMag2 = values[0]*values[0] + values[1]*values[1] + values[2]*values[2];
    if ( seobParams->tortoise ){
        REAL8 a2, rMag, u, u2, u3, u4, u5, w2, D, m1PlusetaKK, bulk, logTerms, deltaU, deltaT, deltaR, eta;
        REAL8 sKerr[3] = {0.};

        eta=seobParams->eobParams->eta;
        SpinEOBHCoeffs * coeffs = seobParams->seobCoeffs;

        a2=0.;
        for(int i=0; i<3; i++){
            sKerr[i] = values[i+6] + values[i+9];
            a2 += sKerr[i]*sKerr[i];
        }

        rMag = sqrt(rMag2);

        u  = 1./rMag;
        u2 = 1./rMag2;
        u3 = u2*u;
        u4 = u2*u2;
        u5 = u4*u;
        w2 = rMag2 + a2;

        /* Eq. 5.83 of BB1, inverse */
        D = 1. + log(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);

        m1PlusetaKK = -1. + eta * coeffs->KK;

        /* Eq. 5.75 of BB1 */
        bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a2*u2;

        /* Eq. 5.73 of BB1 */
        logTerms = 1. + eta*coeffs->k0 + eta*log(fabs(1. + coeffs->k1*u
          + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
          + coeffs->k5*u5 + coeffs->k5l*u5*log(u)));

        /* Eq. 5.73 of BB1 */
        deltaU = fabs(bulk*logTerms);

        /* Eq. 5.71 of BB1 */
        deltaT = rMag2*deltaU;

        /* Eq. 5.38 of BB1 */
        deltaR = deltaT*D;

        csi = sqrt( fabs(deltaT * deltaR) )/ w2; /* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */

        XLALSimIMRCalculateSpinPrecEOBHCoeffs( seobParams->seobCoeffs, seobParams->eobParams->eta, sqrt(a2), seobParams->seobCoeffs->SpinAlignedEOBversion );
    }

    /* Calculate the T-matrix, required to convert P from tortoise to
    * non-tortoise coordinates, and/or vice-versa. This is given explicitly
    * in Eq. A3 of 0912.3466 */
    for( int i = 0; i < 3; i++ ){
        for( int j = 0; j <= i; j++ )
            Tmatrix[i][j] = Tmatrix[j][i] = (values[i]*values[j]/rMag2)* (csi - 1.);
    Tmatrix[i][i]++;
    }

    for ( int i = 3; i < 6; i++ )
        tmpDValues[i]=XLALSpinPrecHcapExactDerivativeWrapper(values,seobParams,i);

    {
    //OPTV3: The following updates hcoeffs
        REAL8 mass1 = seobParams->eobParams->m1;
        REAL8 mass2 = seobParams->eobParams->m2;
        REAL8 s1Data[3],s2Data[3],/*magL,polData[3],*/rcrossrDot[3];
        REAL8 /*Lx,Ly,Lz,*/rcrossrDotMag, s1dotLN, s2dotLN, chiS,chiA, tplspin;
        memcpy(s1Data,values+6,3*sizeof(REAL8));
        memcpy(s2Data,values+9,3*sizeof(REAL8));
        for ( int i = 0; i < 3; i++ )
        {
            s1Data[i] *= (mass1+mass2)*(mass1+mass2);
            s2Data[i] *= (mass1+mass2)*(mass1+mass2);
        }

        /*Compute \vec{L_N} = \vec{r} \times \.{\vec{r}}, \vec{S_i} \dot \vec{L_N} and chiS and chiA */
        rcrossrDot[0] = values[1]*tmpDValues[5] - values[2]*tmpDValues[4];
        rcrossrDot[1] = values[2]*tmpDValues[3] - values[0]*tmpDValues[5];
        rcrossrDot[2] = values[0]*tmpDValues[4] - values[1]*tmpDValues[3];
        rcrossrDotMag = sqrt( rcrossrDot[0]*rcrossrDot[0] + rcrossrDot[1]*rcrossrDot[1] + rcrossrDot[2]*rcrossrDot[2] );
        rcrossrDot[0] /= rcrossrDotMag;
        rcrossrDot[1] /= rcrossrDotMag;
        rcrossrDot[2] /= rcrossrDotMag;
        s1dotLN = (s1Data[0]*rcrossrDot[0] + s1Data[1]*rcrossrDot[1] + s1Data[2]*rcrossrDot[2])/ (mass1*mass1);
        s2dotLN = (s2Data[0]*rcrossrDot[0] + s2Data[1]*rcrossrDot[1] + s2Data[2]*rcrossrDot[2])/ (mass2*mass2) ;
        chiS = 0.5 * (s1dotLN + s2dotLN);
        chiA = 0.5 * (s1dotLN - s2dotLN);
        tplspin = (1.-2.*seobParams->eobParams->eta) * chiS + (mass1 - mass2)/(mass1 + mass2) * chiA;

        XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(seobParams->eobParams->hCoeffs, mass1, mass2, seobParams->eobParams->eta, tplspin,chiS, chiA, 3);
    }

    /* Now make the conversion */
    /* rVectorDot */
    for( int i = 0; i < 3; i++ ){
      dvalues[i]=0.;
      for( int j = 0.; j < 3; j++ )
          dvalues[i] += tmpDValues[j+3]*Tmatrix[i][j];
    }

    /* Now check for NaNs in the dvalues[] output array; only elements 0, 1, and 2 are set */
    for(int i = 0; i < 3; i++ ){
      if( isnan(dvalues[i]) ) {
        XLAL_PRINT_INFO("XLALSpinPrecHcapRvecDerivative_exact::dvalues %3.10f %3.10f %3.10f\n", dvalues[0], dvalues[1], dvalues[2]);
	XLALPrintError( "XLAL Error - %s: nan in the output dvalues \n", __func__);
	XLAL_ERROR( XLAL_EINVAL );
      }
    }


    return XLAL_SUCCESS;
}

#endif //_LALSPINPRECHCAPRVECDERIVATIVE_C
