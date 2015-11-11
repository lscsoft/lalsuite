/*
*  Copyright (C) 2010 Craig Robinson, Yi Pan
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
 * \brief In newer versions of the EOBNR approximant, we
 * do not have an analytic expression for the derivative of the waveform.
 * As such, it is necessary to calculate the derivatives numerically. This
 * function provides the means to do just that.
 *
 */

#ifndef _LALSIMIMRSPINEOBHCAPEXACTDERIVATIVE_C
#define _LALSIMIMRSPINEOBHCAPEXACTDERIVATIVE_C

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <math.h>
#include <gsl/gsl_deriv.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinEOBHamiltonianOptimized.c"

//int UsePrec = 1;

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */


//static REAL8 GSLSpinHamiltonianWrapperOptimized( REAL8 x, void *params );

static REAL8 XLALSpinHcapExactDerivWRTParam(
                       const INT4 paramIdx,
                       const REAL8 values[],
                       SpinEOBParams *params
                       );

static REAL8 XLALSpinHcapHybDerivWRTParam(
                                          const INT4 paramIdx,      /**<< Index of the parameters */
                                          const REAL8 values[],     /**<< Dynamical variables */
                                          SpinEOBParams *funcParams /**<< EOB Parameters */
                                          );

static REAL8 GSLSpinAlignedHamiltonianWrapper_derivs_allatonce( REAL8 output[6], const REAL8 input[6], void *params );

static double GSLSpinAlignedHamiltonianWrapper_ExactDeriv( double x, void *params );

static REAL8 XLALSimIMRSpinEOBHamiltonian_ExactDeriv(
                                                     INT4 which_to_vary,                  /**<< Take a derivative with respect to "which_to_vary" variable */
                                                     const REAL8    eta,                  /**<< Symmetric mass ratio */
                                                     REAL8Vector    * restrict x,         /**<< Position vector */
                                                     REAL8Vector    * restrict p,         /**<< Momentum vector (tortoise radial component pr*) */
                                                     REAL8Vector    * restrict s1Vec,     /**<< Spin vector 1 */
                                                     REAL8Vector    * restrict s2Vec,     /**<< Spin vector 2 */
                                                     REAL8Vector    * restrict sigmaKerr, /**<< Spin vector sigma_kerr */
                                                     REAL8Vector    * restrict sigmaStar, /**<< Spin vector sigma_star */
                                                     INT4                      tortoise,  /**<< flag to state whether the momentum is the tortoise co-ord */
                                                     SpinEOBHCoeffs *coeffs               /**<< Structure containing various coefficients */
                                                     );

static REAL8 XLALSimIMRSpinEOBHamiltonian_derivs_allatonce(
							   REAL8 output[6],                     /**<< Output vector (contains all derivatives, though WARNING: known unused derivatives may be set to zero) */
							   const REAL8    eta,                  /**<< Symmetric mass ratio */
							   REAL8Vector    * restrict x,         /**<< Position vector */
							   REAL8Vector    * restrict p,	        /**<< Momentum vector (tortoise radial component pr*) */
							   REAL8Vector    * restrict s1Vec,     /**<< Spin vector 1 */
							   REAL8Vector    * restrict s2Vec,     /**<< Spin vector 2 */
							   REAL8Vector    * restrict sigmaKerr, /**<< Spin vector sigma_kerr */
							   REAL8Vector    * restrict sigmaStar, /**<< Spin vector sigma_star */
							   INT4                      tortoise,  /**<< flag to state whether the momentum is the tortoise co-ord */
							   SpinEOBHCoeffs *coeffs               /**<< Structure containing various coefficients */
							   );
/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

static REAL8 XLALSpinHcapHybDerivWRTParam(
                                          const INT4 paramIdx,      /**<< Index of the parameters */
                                          const REAL8 values[],     /**<< Dynamical variables */
                                          SpinEOBParams *funcParams /**<< EOB Parameters */
                                          )
{
  HcapDerivParams params;
  /* Set up pointers for GSL */
  params.values  = values;
  params.params  = funcParams;
  params.varyParam = paramIdx;

  REAL8 result;
  if(paramIdx>=0 && paramIdx<6) {
    REAL8 output[6];
    GSLSpinAlignedHamiltonianWrapper_derivs_allatonce( output, values, &params );
    result=output[paramIdx];
  } else {
    XLAL_ERROR( XLAL_EFUNC );
  }

  return result;
}

/* Wrapper for computing exact derivatives of the Hamiltonian */
static double GSLSpinAlignedHamiltonianWrapper_ExactDeriv( double x, void *params )
{
  HcapDerivParams *dParams = (HcapDerivParams *)params;

  EOBParams *eobParams = dParams->params->eobParams;

  REAL8 tmpVec[6];

  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p;
  REAL8Vector *s1Vec = dParams->params->s1Vec;
  REAL8Vector *s2Vec = dParams->params->s2Vec;
  REAL8Vector *sigmaKerr = dParams->params->sigmaKerr;
  REAL8Vector *sigmaStar = dParams->params->sigmaStar;

  /* Use a temporary vector to avoid corrupting the main function */
  memcpy( tmpVec, dParams->values,
	  sizeof(tmpVec) );

  /* Set the relevant entry in the vector to the correct value */
  tmpVec[dParams->varyParam] = x;

  /* Set the LAL-style vectors to point to the appropriate things */
  r.length = p.length = 3;
  r.data     = tmpVec;
  p.data     = tmpVec+3;

  return XLALSimIMRSpinEOBHamiltonian_ExactDeriv(dParams->varyParam, eobParams->eta, &r, &p, s1Vec, s2Vec, sigmaKerr, sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs ) / eobParams->eta;
}

/* Wrapper for GSL to call the Hamiltonian function */
static REAL8 GSLSpinAlignedHamiltonianWrapper_derivs_allatonce( REAL8 output[6], const REAL8 input[6], void *params )
{
  HcapDerivParams *dParams = (HcapDerivParams *)params;

  EOBParams *eobParams = dParams->params->eobParams;

  REAL8 tmpVec[6];
  REAL8 returnval=0;
  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p;
  REAL8Vector *s1Vec = dParams->params->s1Vec;
  REAL8Vector *s2Vec = dParams->params->s2Vec;
  REAL8Vector *sigmaKerr = dParams->params->sigmaKerr;
  REAL8Vector *sigmaStar = dParams->params->sigmaStar;

  /* Use a temporary vector to avoid corrupting the main function */
  memcpy( tmpVec, dParams->values,
	  sizeof(tmpVec) );

  /* Set the relevant entry in the vector to the correct value */
  for(int i=0;i<6;i++) tmpVec[i] = input[i];
  //tmpVec[dParams->varyParam] = x;

  /* Set the LAL-style vectors to point to the appropriate things */
  r.length = p.length = 3;
  r.data     = tmpVec;
  p.data     = tmpVec+3;

  XLALSimIMRSpinEOBHamiltonian_derivs_allatonce(output, eobParams->eta, &r, &p, s1Vec, s2Vec, sigmaKerr, sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs );

  for(int i=0;i<6;i++) output[i] /= eobParams->eta;

  return returnval;
}

/**
 * Calculate the derivative of the Hamiltonian w.r.t. a specific parameter
 * Used by generic spin EOB model, including initial conditions solver.
 */
static REAL8 XLALSpinHcapExactDerivWRTParam(
                                            const INT4 paramIdx,      /**<< Index of the parameters */
                                            const REAL8 values[],     /**<< Dynamical variables */
                                            SpinEOBParams *funcParams /**<< SEOB Parameters */
                                            )
{
  REAL8 result;

  REAL8Vector *sigmaKerr = funcParams->sigmaKerr;
  REAL8Vector *sigmaStar = funcParams->sigmaStar;
  SpinEOBHCoeffs *coeffs = funcParams->seobCoeffs;
  REAL8Vector *s1Vec = funcParams->s1Vec;
  REAL8Vector *s2Vec = funcParams->s2Vec;
  REAL8Vector *x=NULL;
  REAL8Vector *p=NULL;
  /* Note that this function is called a limited number of times in the initial condition setup,
   * so no worries about the below memory allocations slowing the code... at this point. */
  x=XLALCreateREAL8Vector(3);
  p=XLALCreateREAL8Vector(3);
  memcpy(x->data,&values[0],3*sizeof(REAL8));
  memcpy(p->data,&values[3],3*sizeof(REAL8));

  REAL8 eta = funcParams->eobParams->eta;
  REAL8 e3z = (0.0 <= sigmaKerr->data[2]) - (sigmaKerr->data[2] < 0.0); // This is a modified sign function: e3z = 1 if sigmaKerr->data[2]>=0, -1 otherwise

  /* Now calculate derivatives w.r.t. the required parameter */
  if (funcParams->tortoise==1){
    if (paramIdx==4){
#include "mathematica_codes/SEOBNRv2_opt_dp1tortoise.h"
      result = dHdp1/eta;
    } else if (paramIdx==3){
#include "mathematica_codes/SEOBNRv2_opt_dp0tortoise.h"
      result = dHdp0/eta;
    } else if (paramIdx==5){
#include "mathematica_codes/SEOBNRv2_opt_dp2tortoise.h"
      result = dHdp2/eta;
    } else if (paramIdx==0){
#include "mathematica_codes/SEOBNRv2_opt_dx0tortoise.h"
      result = dHdx0/eta;
    } else {
      XLAL_ERROR( XLAL_EFUNC );
    }
  }else if(funcParams->tortoise==0) {
    if (paramIdx==4){
#include "mathematica_codes/SEOBNRv2_opt_dp1.h"
      result = dHdp1/eta;
    } else if (paramIdx==3){
#include "mathematica_codes/SEOBNRv2_opt_dp0.h"
      result = dHdp0/eta;
    } else if (paramIdx==5){
#include "mathematica_codes/SEOBNRv2_opt_dp2.h"
      result = dHdp2/eta;
    } else if (paramIdx==0){
#include "mathematica_codes/SEOBNRv2_opt_dx0.h"
      result = dHdx0/eta;
    } else {
      XLAL_ERROR( XLAL_EFUNC );
    }
  } else {
      XLAL_ERROR( XLAL_EFUNC );
  }
  XLALDestroyREAL8Vector(x);
  XLALDestroyREAL8Vector(p);

  return result;
}

/* Evaluates derivative of Hamiltonian with respect to py.
 * We've made it extensible in case we'd like to take other derivs. */
static REAL8 XLALSimIMRSpinEOBHamiltonian_ExactDeriv(
                                                     INT4 which_to_vary,
                                                     const REAL8    eta,                  /**<< Symmetric mass ratio */
                                                     REAL8Vector    * restrict x,         /**<< Position vector */
                                                     REAL8Vector    * restrict p,         /**<< Momentum vector (tortoise radial component pr*) */
                                                     REAL8Vector    * restrict s1Vec,     /**<< Spin vector 1 */
                                                     REAL8Vector    * restrict s2Vec,     /**<< Spin vector 2 */
                                                     REAL8Vector    * restrict sigmaKerr, /**<< Spin vector sigma_kerr */
                                                     REAL8Vector    * restrict sigmaStar, /**<< Spin vector sigma_star */
                                                     INT4                      tortoise,  /**<< flag to state whether the momentum is the tortoise co-ord */
                                                     SpinEOBHCoeffs * coeffs              /**<< Structure containing various coefficients */
                                                     )
{
  REAL8 e3z = (0.0 <= sigmaKerr->data[2]) - (sigmaKerr->data[2] < 0.0); // This is a modified sign function: e3z = 1 if sigmaKerr->data[2]>=0, -1 otherwise
  if(tortoise==1) {
    switch(which_to_vary) {
    case 4:
      {
#include "mathematica_codes/SEOBNRv2_opt_pdata1tortoise.h"
	return Hreal;
      }
      break;
    default:
      {
	printf("Option not supported: %d!\n",which_to_vary); exit(1);
	break;
      }
    }
  } else {
    switch(which_to_vary) {
    case 4:
      {
#include "mathematica_codes/SEOBNRv2_opt_pdata1.h"
	return Hreal;
      }
      break;
    default:
      {
	printf("Option not supported: %d!\n",which_to_vary); exit(1);
	break;
      }
    }
  }
  return -1e-300;
}


/**
 * Wrapper for GSL to call the Hamiltonian function
 */
/* static REAL8 GSLSpinHamiltonianWrapperOptimized( REAL8 x, void *params ) */
/* { */
/*   HcapDerivParams *dParams = (HcapDerivParams *)params; */

/*   EOBParams *eobParams = dParams->params->eobParams; */

/*   REAL8 tmpVec[12]; */
/*   REAL8 s1normData[3], s2normData[3], sKerrData[3], sStarData[3]; */

/*   /\* These are the vectors which will be used in the call to the Hamiltonian *\/ */
/*   REAL8Vector r, p, spin1, spin2, spin1norm, spin2norm; */
/*   REAL8Vector sigmaKerr, sigmaStar; */

/*   int i; */
/*   REAL8 a; */
/*   REAL8 m1 = eobParams->m1; */
/*   REAL8 m2 = eobParams->m2; */
/*   REAL8 mT2 = (m1+m2)*(m1+m2); */

/*   /\* Use a temporary vector to avoid corrupting the main function *\/ */
/*   memcpy( tmpVec, dParams->values, sizeof(tmpVec) ); */

/*   /\* Set the relevant entry in the vector to the correct value *\/ */
/*   tmpVec[dParams->varyParam] = x; */

/*   /\* Set the LAL-style vectors to point to the appropriate things *\/ */
/*   r.length = p.length = spin1.length = spin2.length = spin1norm.length = spin2norm.length = 3; */
/*   sigmaKerr.length = sigmaStar.length = 3; */
/*   r.data     = tmpVec; */
/*   p.data     = tmpVec+3; */
/*   spin1.data = tmpVec+6; */
/*   spin2.data = tmpVec+9; */
/*   spin1norm.data = s1normData; */
/*   spin2norm.data = s2normData; */
/*   sigmaKerr.data = sKerrData; */
/*   sigmaStar.data = sStarData; */

/*   memcpy( s1normData, tmpVec+6, 3*sizeof(REAL8) ); */
/*   memcpy( s2normData, tmpVec+9, 3*sizeof(REAL8) ); */

/*   for ( i = 0; i < 3; i++ ) */
/*   { */
/*      s1normData[i] /= mT2; */
/*      s2normData[i] /= mT2; */
/*   } */

/*   /\* Calculate various spin parameters *\/ */
/*   XLALSimIMRSpinEOBCalculateSigmaKerr( &sigmaKerr, eobParams->m1, */
/* 				eobParams->m2, &spin1, &spin2 ); */
/*   XLALSimIMRSpinEOBCalculateSigmaStar( &sigmaStar, eobParams->m1, */
/* 				eobParams->m2, &spin1, &spin2 ); */
/*   a = sqrt( sigmaKerr.data[0]*sigmaKerr.data[0] + sigmaKerr.data[1]*sigmaKerr.data[1] */
/*               + sigmaKerr.data[2]*sigmaKerr.data[2] ); */
/*   //printf( "a = %e\n", a ); */
/*   //printf( "aStar = %e\n", sqrt( sigmaStar.data[0]*sigmaStar.data[0] + sigmaStar.data[1]*sigmaStar.data[1] + sigmaStar.data[2]*sigmaStar.data[2]) ); */
/*   if ( isnan( a ) ) */
/*   { */
/*     printf( "a is nan!!\n"); */
/*   } */
/*   //XLALSimIMRCalculateSpinEOBHCoeffs( dParams->params->seobCoeffs, eobParams->eta, a ); */
/*   if (UsePrec) */
/*   { */
/*     /\* Set up structures and calculate necessary PN parameters *\/ */
/*     /\* Due to precession, these need to get calculated in every step *\/ */
/*     /\* TODO: Only calculate non-spinning parts once *\/ */
/*     memset( dParams->params->seobCoeffs, 0, sizeof(SpinEOBHCoeffs) ); */

/*     /\* Update the Z component of the Kerr background spin */
/*      * Pre-compute the Hamiltonian coefficients            *\/ */
/*     REAL8Vector *delsigmaKerr 	    = dParams->params->sigmaKerr; */
/*     dParams->params->sigmaKerr 		= &sigmaKerr; */
/*     dParams->params->a 				= a; */
/*     //tmpsigmaKerr->data[2]; */
/*     if ( XLALSimIMRCalculateSpinEOBHCoeffs( dParams->params->seobCoeffs, */
/* 			eobParams->eta, a, */
/* 			dParams->params->seobCoeffs->SpinAlignedEOBversion ) == XLAL_FAILURE ) */
/*     { */
/*       XLALDestroyREAL8Vector( dParams->params->sigmaKerr ); */
/*       XLAL_ERROR( XLAL_EFUNC ); */
/*     } */

/*     /\* Release the old memory *\/ */
/*     XLALDestroyREAL8Vector( delsigmaKerr ); */
/*   } */

/*   //printf( "Hamiltonian = %e\n", XLALSimIMRSpinEOBHamiltonian( eobParams->eta, &r, &p, &sigmaKerr, &sigmaStar, dParams->params->seobCoeffs ) ); */
/*   return XLALSimIMRSpinEOBHamiltonianOptimized( eobParams->eta, &r, &p, &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs ) / eobParams->eta; */
/* } */

static REAL8 XLALSimIMRSpinEOBHamiltonian_derivs_allatonce(
							   REAL8 output[6],
							   const REAL8    eta,                  /**<< Symmetric mass ratio */
							   REAL8Vector    * restrict x,         /**<< Position vector */
							   REAL8Vector    * restrict p,	    /**<< Momentum vector (tortoise radial component pr*) */
							   REAL8Vector    * restrict s1Vec,     /**<< Spin vector 1 */
							   REAL8Vector    * restrict s2Vec,     /**<< Spin vector 2 */
							   REAL8Vector    * restrict sigmaKerr, /**<< Spin vector sigma_kerr */
							   REAL8Vector    * restrict sigmaStar, /**<< Spin vector sigma_star */
							   INT4                      tortoise,  /**<< flag to state whether the momentum is the tortoise co-ord */
							   SpinEOBHCoeffs *coeffs               /**<< Structure containing various coefficients */
							   )
{
  REAL8 returnval=0;
  REAL8 e3z = (0.0 <= sigmaKerr->data[2]) - (sigmaKerr->data[2] < 0.0); // This is a modified sign function: e3z = 1 if sigmaKerr->data[2]>=0, -1 otherwise
  if(tortoise==1) {
    // FASTEST OPTION. Note that it sets g2=[xdata2 derivative]=0.
#include "mathematica_codes/SEOBNRv2_opt_3derivstortoise.h"
    output[0]=g0;
    output[1]=g1;
    output[2]=g2;
    output[3]=g3;
    output[4]=g4;
    output[5]=g5;
  } else {
#include "mathematica_codes/SEOBNRv2_opt_3derivs.h"
    output[0]=g0;
    output[1]=g1;
    output[2]=g2;
    output[3]=g3;
    output[4]=g4;
    output[5]=g5;//    printf("ERROR: NOT EXPECTING THAT %d.\n",tortoise); exit(1);
  }

  return returnval;
}
#endif /* _LALSIMIMRSPINEOBHCAPEXACTDERIVATIVE_C */
