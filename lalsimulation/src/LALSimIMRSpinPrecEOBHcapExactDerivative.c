/** Static function file **/
#ifndef _LALSIMIMRSPINPRECEOBHCAPEXACTDERIVATIVE_C
#define _LALSIMIMRSPINPRECEOBHCAPEXACTDERIVATIVE_C

//#include "LALSimIMRSpinEOBHamiltonian.c"
//#include "LALSimIMRSpinEOBHamiltonianPrec.c"

#include "LALSimIMRCalculateSpinPrecEOBHCoeffs.c" /* OPTV3 */

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

static REAL8 XLALSpinPrecHcapExactDerivativeWrapper(
            const REAL8 * values, /**<< Dynamical variables */
            SpinEOBParams * params, /**<< EOB Parameters */
            INT4 variedParam /**<< Index of variable to be differentiated w.r.t.*/
);

static REAL8 XLALSpinPrecHcapExactDerivativeNoWrap(
              const REAL8 * values,     /**<< Dynamical variables */
              INT4 variedParam,         /**<< Index of variable to be differentiated w.r.t.*/
              INT4 tortoise,            /**<< Tortoise coordinate indicator */
              REAL8 eta,                /**<< Symmetric mass ratio */
              REAL8Vector * sigmaKerr,  /**<< Spin vector sigma_kerr */
              REAL8Vector * sigmaStar,  /**<< Spin vector sigma_star */
              REAL8Vector * s1Vec,      /**<< Spin 1 vector */
              REAL8Vector * s2Vec,      /**<< Spin 2 vector */
              SpinEOBHCoeffs *coeffs,   /**<< Pre-computed coefficients which appear in the function */
              SEOBHCoeffConstants * seobCoeffConsts, /**<< Pre-computed Constants of SEOB coefficient calculation */
              REAL8 mass1,              /**<< mass 1 */
              REAL8 mass2               /**<< mass 2 */
              );

/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */


static REAL8 XLALSpinPrecHcapExactDerivativeWrapper(const REAL8 * values, SpinEOBParams * params, INT4 variedParam){
  REAL8 mass1 = params->eobParams->m1;
  REAL8 mass2 = params->eobParams->m2;
  REAL8 eta = params->eobParams->eta;
  INT4 tortoise = params->tortoise;
  REAL8 deriv;

  SpinEOBHCoeffs *coeffs = params->seobCoeffs;

  REAL8 sigmaKerrData[3];
  REAL8 sigmaStarData[3];
  REAL8 s1VecData[3];
  REAL8 s2VecData[3];

  REAL8Vector sigmaKerr;
  REAL8Vector sigmaStar;
  REAL8Vector s1Vec;
  REAL8Vector s2Vec;

  s1Vec.data = s1VecData;
  s2Vec.data = s2VecData;
  sigmaKerr.data = sigmaKerrData;
  sigmaStar.data = sigmaStarData;

  memcpy(s1VecData, values + 6, 3 * sizeof(REAL8));
  memcpy(s2VecData, values + 9, 3 * sizeof(REAL8));

  /**
   * Loop to calculate normalized spin of the deformed-Kerr background (sigmaKerr)
   * and normalized spin of the test particle (sigmaStarr) in SEOBNR.
   * Eqs. 5.2 and 5.3 of Barausse and Buonanno PRD 81, 084024 (2010).
   */

  /* Constants defined to reduce the number of floating point operations
     needed to compute sigmaStarData */
  REAL8 m1m2_inv = 1.0/(mass1*mass2);
  REAL8 m2overm1 = mass2*mass2*m1m2_inv;
  REAL8 m1overm2 = mass1*mass1*m1m2_inv;

  for(int i=0; i<3; i++){
    sigmaKerrData[i] = s1VecData[i] + s2VecData[i];
    sigmaStarData[i] = m2overm1*s1VecData[i] + m1overm2*s2VecData[i];
  }

  SEOBHCoeffConstants * seobCoeffConsts = params->seobCoeffConsts;
  deriv = XLALSpinPrecHcapExactDerivativeNoWrap(values, variedParam, tortoise, eta, &sigmaKerr, &sigmaStar, &s1Vec, &s2Vec, coeffs, seobCoeffConsts, mass1, mass2)/eta;
  return deriv;
}

/**
 * Call a particular derivative
 */
static REAL8 XLALSpinPrecHcapExactDerivativeNoWrap(
              const REAL8 * values,                 /**<< Dynamical variables */
              INT4 variedParam,                     /**<< Index of variable to be differentiated w.r.t.*/
              INT4 tortoise,                        /**<< Tortoise coordinate indicator */
              REAL8 eta,                            /**<< Symmetric mass ratio */
              REAL8Vector * sigmaKerr,              /**<< Spin vector sigma_kerr */
              REAL8Vector * sigmaStar,              /**<< Spin vector sigma_star */
              REAL8Vector * s1Vec,                  /**<< Spin 1 vector */
              REAL8Vector * s2Vec,                  /**<< Spin 2 vector */
              SpinEOBHCoeffs * coeffs,              /**<< Pre-computed coefficients which appear in the function */
              SEOBHCoeffConstants * seobCoeffConsts, /**<< SpinEOBH Coefficient Constants */
              REAL8 mass1,                          /**<< mass 1 */
              REAL8 mass2                           /**<< mass 2 */
              )
{
    REAL8 a2 = sigmaKerr->data[0] * sigmaKerr->data[0] + sigmaKerr->data[1] * sigmaKerr->data[1] +  sigmaKerr->data[2] * sigmaKerr->data[2];
    REAL8 a = sqrt(a2);

    if(coeffs->updateHCoeffs){
      SpinEOBHCoeffs tmpCoeffs;

      if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs( &tmpCoeffs, eta, a, coeffs->SpinAlignedEOBversion ) == XLAL_FAILURE )
	{
	  XLAL_ERROR( XLAL_EFUNC );
	}

      tmpCoeffs.SpinAlignedEOBversion = coeffs->SpinAlignedEOBversion;
      tmpCoeffs.updateHCoeffs = coeffs->updateHCoeffs;

      coeffs = &tmpCoeffs;
    }

    REAL8 deriv, e3_x, e3_y, e3_z;

    REAL8Vector xVec, pVec;
    xVec.length = pVec.length= 3;

    REAL8 xData[3] = {0.}, pData[3] = {0.};
    xVec.data = xData;
    pVec.data = pData;

    REAL8Vector * x;
    REAL8Vector * p;

    x=&xVec;
    p=&pVec;

    memcpy(xVec.data,values,3*sizeof(REAL8));
    memcpy(pVec.data,values+3,3*sizeof(REAL8));

    INT4 divby0 = 0;

    if(a !=0.)
    {
      const REAL8 inva = 1./a;
      e3_x = sigmaKerr->data[0] * inva;
      e3_y = sigmaKerr->data[1] * inva;
      e3_z = sigmaKerr->data[2] * inva;
    }
    else
    {
      /*OPTV3: Since spin=0, we are free to choose the "spin direction".*/
      e3_x = 1./sqrt(3.);
      e3_y = 1./sqrt(3.);
      e3_z = 1./sqrt(3.);

      divby0=1;
    }

    const REAL8 invr = 1./sqrt(xData[0]*xData[0]+xData[1]*xData[1]+xData[2]*xData[2]);

    if (1. - fabs(e3_x*(xData[0]*invr) + e3_y*(xData[1]*invr) + e3_z*(xData[2]*invr)) <= 1.e-8) {
      /** Prevent the source frame basis vectors from becoming too short, which may
       * cause denominators to become zero in exact derivative expressions.
       */
      e3_x = e3_x+0.1;
      e3_y = e3_y+0.1;
      const REAL8 invnorm = 1./sqrt(e3_x*e3_x + e3_y*e3_y + e3_z*e3_z);
      e3_x = e3_x*invnorm;
      e3_y = e3_y*invnorm;
      e3_z = e3_z*invnorm;
      divby0 = 1;
    }

    if(divby0) {
      /** Vector sum of s1 & s2 cannot be zero when taking spin derivatives, because
       * vector sum components appear in denominators of exact derivative expressions.
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

    /** KK defined in Barausse and Buonanno PRD 81, 084024 (2010)
     * Equations 6.11 and 6.13.
     */
    REAL8 m1PlusEtaKK = coeffs->KK*eta-1.0;

    //OPTV3: Must define the seob coeffs constants for spin derivatives
    REAL8 c0k2 = seobCoeffConsts->a0k2;
    REAL8 c1k2 = seobCoeffConsts->a1k2;
    REAL8 c0k3 = seobCoeffConsts->a0k3;
    REAL8 c1k3 = seobCoeffConsts->a1k3;
    REAL8 c0k4 = seobCoeffConsts->a0k4;
    REAL8 c1k4 = seobCoeffConsts->a1k4;
    REAL8 c2k4 = seobCoeffConsts->a2k4;
    REAL8 c0k5 = seobCoeffConsts->a0k5;
    REAL8 c1k5 = seobCoeffConsts->a1k5;
    REAL8 c2k5 = seobCoeffConsts->a2k5;

    /** variedParam specifies which partial derivative to compute;
     * the following switch structure sorts derivatives both by tortoise
     * value and divby0 status.  Both tortoise and divby0 act as
     * parameters that affect partial derivative computations.
     */
    if(variedParam <= 5 || !divby0 ){
        switch(tortoise){
        case 1:
            switch(variedParam){
                case 0:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dx0_Tortoise-1.h"
                    deriv=d100000;}
                    break;
                case 1:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dx1_Tortoise-1.h"
                    deriv=d010000;}
                    break;
                case 2:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dx2_Tortoise-1.h"
                    deriv=d001000;}
                    break;
                case 3:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dp0_Tortoise-1.h"
                    deriv=d000100;}
                    break;
                case 4:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dp1_Tortoise-1.h"
                    deriv=d000010;}
                    break;
                case 5:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dp2_Tortoise-1.h"
                    deriv=d000001;}
                    break;
                case 6:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-0_Tortoise-1.h"
                    deriv=ds100000;}
                    break;
                case 7:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-1_Tortoise-1.h"
                    deriv=ds010000;}
                    break;
                case 8:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-2_Tortoise-1.h"
                    deriv=ds001000;}
                    break;
                case 9:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-0_Tortoise-1.h"
                    deriv=ds000100;}
                    break;
                case 10:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-1_Tortoise-1.h"
                    deriv=ds000010;}
                    break;
                case 11:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-2_Tortoise-1.h"
                    deriv=ds000001;}
                    break;
                default:
                    printf("Err: Unexpected Index %d\n",variedParam);
                    XLAL_ERROR( XLAL_EFUNC );
                    exit(1);
                    break;

            }
            break;
        case 0:
            switch(variedParam){
                case 0:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dx0_Tortoise-0.h"
                    deriv=d100000;}
                    break;
                case 1:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dx1_Tortoise-0.h"
                    deriv=d010000;}
                    break;
                case 2:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dx2_Tortoise-0.h"
                    deriv=d001000;}
                    break;
                case 3:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dp0_Tortoise-0.h"
                    deriv=d000100;}
                    break;
                case 4:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dp1_Tortoise-0.h"
                    deriv=d000010;}
                    break;
                case 5:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dp2_Tortoise-0.h"
                    deriv=d000001;}
                    break;
                case 6:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-0_Tortoise-0.h"
                    deriv=ds100000;}
                    break;
                case 7:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-1_Tortoise-0.h"
                    deriv=ds010000;}
                    break;
                case 8:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-2_Tortoise-0.h"
                    deriv=ds001000;}
                    break;
                case 9:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-0_Tortoise-0.h"
                    deriv=ds000100;}
                    break;
                case 10:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-1_Tortoise-0.h"
                    deriv=ds000010;}
                    break;
                case 11:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-2_Tortoise-0.h"
                    deriv=ds000001;}
                    break;
                default:
                    printf("Err: Unexpected Index %d\n",variedParam);
                    XLAL_ERROR( XLAL_EFUNC );
                    exit(1);
                    break;
            }
        break;
        case 2:
            switch(variedParam){
                case 0:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dx0_Tortoise-2.h"
                    deriv=d100000;}
                    break;
                case 1:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dx1_Tortoise-2.h"
                    deriv=d010000;}
                    break;
                case 2:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dx2_Tortoise-2.h"
                    deriv=d001000;}
                    break;
                case 3:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dp0_Tortoise-2.h"
                    deriv=d000100;}
                    break;
                case 4:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dp1_Tortoise-2.h"
                    deriv=d000010;}
                    break;
                case 5:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_dp2_Tortoise-2.h"
                    deriv=d000001;}
                    break;
                case 6:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-0_Tortoise-2.h"
                    deriv=ds100000;}
                    break;
                case 7:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-1_Tortoise-2.h"
                    deriv=ds010000;}
                    break;
                case 8:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-2_Tortoise-2.h"
                    deriv=ds001000;}
                    break;
                case 9:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-0_Tortoise-2.h"
                    deriv=ds000100;}
                    break;
                case 10:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-1_Tortoise-2.h"
                    deriv=ds000010;}
                    break;
                case 11:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-2_Tortoise-2.h"
                    deriv=ds000001;}
                    break;
                default:
                    printf("Err: Unexpected Index %d\n",variedParam);
                    XLAL_ERROR( XLAL_EFUNC );
                    exit(1);
                    break;
            }
            break;
        default:
            printf("STRANGE TORTOISE: %d\n",tortoise);
            XLAL_ERROR( XLAL_EFUNC );
            exit(1);
            break;
        }
    }else{

        switch(tortoise){
        case 1:
            switch(variedParam){
                case 6:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-0_Tortoise-1_D-e3.h"
                    deriv=ds100000;}
                    break;
                case 7:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-1_Tortoise-1_D-e3.h"
                    deriv=ds010000;}
                    break;
                case 8:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-2_Tortoise-1_D-e3.h"
                    deriv=ds001000;}
                    break;
                case 9:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-0_Tortoise-1_D-e3.h"
                    deriv=ds000100;}
                    break;
                case 10:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-1_Tortoise-1_D-e3.h"
                    deriv=ds000010;}
                    break;
                case 11:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-2_Tortoise-1_D-e3.h"
                    deriv=ds000001;}
                    break;
                default:
                    printf("Err: Unexpected Index %d\n",variedParam);
                    XLAL_ERROR( XLAL_EFUNC );
                    exit(1);
                    break;

            }
            break;
        case 0:
            switch(variedParam){
                case 6:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-0_Tortoise-0_D-e3.h"
                    deriv=ds100000;}
                    break;
                case 7:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-1_Tortoise-0_D-e3.h"
                    deriv=ds010000;}
                    break;
                case 8:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-2_Tortoise-0_D-e3.h"
                    deriv=ds001000;}
                    break;
                case 9:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-0_Tortoise-0_D-e3.h"
                    deriv=ds000100;}
                    break;
                case 10:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-1_Tortoise-0_D-e3.h"
                    deriv=ds000010;}
                    break;
                case 11:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-2_Tortoise-0_D-e3.h"
                    deriv=ds000001;}
                    break;
                default:
                    printf("Err: Unexpected Index %d\n",variedParam);
                    XLAL_ERROR( XLAL_EFUNC );
                    exit(1);
                    break;
            }
        break;
        case 2:
            switch(variedParam){
                case 6:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-0_Tortoise-2_D-e3.h"
                    deriv=ds100000;}
                    break;
                case 7:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-1_Tortoise-2_D-e3.h"
                    deriv=ds010000;}
                    break;
                case 8:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds1-2_Tortoise-2_D-e3.h"
                    deriv=ds001000;}
                    break;
                case 9:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-0_Tortoise-2_D-e3.h"
                    deriv=ds000100;}
                    break;
                case 10:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-1_Tortoise-2_D-e3.h"
                    deriv=ds000010;}
                    break;
                case 11:{
                    #include "mathematica_codes/SEOBNRv3_opt/SEOBNRv3_opt_ds2-2_Tortoise-2_D-e3.h"
                    deriv=ds000001;}
                    break;
                default:
                    printf("Err: Unexpected Index %d\n",variedParam);
                    XLAL_ERROR( XLAL_EFUNC );
                    exit(1);
                    break;
            }
            break;
        default:
            printf("STRANGE TORTOISE: %d\n",tortoise);
            XLAL_ERROR( XLAL_EFUNC );
            exit(1);
            break;
        }
    }
    if(isnan(deriv)){
        printf("NAN in derivative! tortoise=%d; variedParam=%d; divby0=%d\n",tortoise,variedParam,divby0);
        XLAL_ERROR( XLAL_EFUNC );
        exit(1);
    }
	return deriv;
}

#endif // _LALSIMIMRSPINPRECEOBHCAPEXACTDERIVATIVE_C
