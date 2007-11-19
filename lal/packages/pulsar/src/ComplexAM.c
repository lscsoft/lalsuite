/*
 * Copyright (C) 2007 John T. Whelan
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

/** \author J. T. Whelan
 * \file 
 * \brief
 * Functions related to F-statistic calculation when the AM coefficients are complex.
 *                                                                          
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include <lal/ExtrapolatePulsarSpins.h>

/* GSL includes */
#include <lal/LALGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>


#include <lal/AVFactories.h>
#include "ComplexAM.h"

NRCSID( COMPLEXAMC, "$Id$");

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/* empty initializers  */
const Fcomponents empty_Fcomponents;
const MultiDetectorStateSeries empty_MultiDetectorStateSeries;

/*==================== FUNCTION DEFINITIONS ====================*/



/** Function to compute (multi-IFO) F-statistic for given parameter-space point ::psPoint,
 *  normalized SFT-data (normalized by <em>double-sided</em> PSD Sn), noise-weights
 *  and detector state-series 
 *
 * NOTE: for better efficiency some quantities that need to be recomputed only for different 
 * sky-positions are buffered in \a cfBuffer if given. 
 * - In order to 'empty' this buffer (at the end) use XLALEmptyComputeFBuffer()
 * - You CAN pass NULL for the \a cfBuffer if you don't want to use buffering (slower).
 *
 * NOTE2: there's a spaceholder for binary-pulsar parameters in \a psPoint, but this 
 * it not implemented yet.
 *
 * \note THIS VERSION IS FOR THE CASE WHERE THE DETECTOR TENSOR IS COMPLEX
 * AND THE DESCRIPTION NEEDS TO BE MODIFIED!
 *
 */
void
CmplxComputeFStat ( LALStatus *status, 
		    Fcomponents *Fstat,                 /**< [out] Fstatistic + Fa, Fb */
		    const PulsarDopplerParams *doppler, /**< parameter-space point to compute F for */
		    const MultiSFTVector *multiSFTs,    /**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
		    const MultiNoiseWeights *multiWeights,/**< noise-weights of all SFTs */
		    const MultiCmplxDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
		    const ComputeFParams *params,       /**< addition computational params */
		    CmplxComputeFBuffer *cfBuffer            /**< CF-internal buffering structure */
	       )
{
  Fcomponents retF = empty_Fcomponents;
  UINT4 X, numDetectors;	
  MultiSSBtimes *multiSSB = NULL;
  MultiCmplxAMCoeffs *multiAMcoef = NULL;
  REAL8 Ad, Bd, Cd, Dd_inv;

  INITSTATUS( status, "CmplxComputeFStat", COMPLEXAMC );
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT ( Fstat, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL );
  ASSERT ( multiSFTs, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL );
  ASSERT ( doppler, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL );
  ASSERT ( params, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL );

  numDetectors = multiSFTs->length;
  ASSERT ( multiDetStates->length == numDetectors, status, COMPLEXAMC_EINPUT, COMPLEXAMC_MSGEINPUT );
  if ( multiWeights ) {
    ASSERT ( multiWeights->length == numDetectors , status, COMPLEXAMC_EINPUT, COMPLEXAMC_MSGEINPUT );
  }

  if ( doppler->orbit ) {
    LALPrintError ("\nSorry, binary-pulsar search not yet implemented in CmplxComputeFStat()\n\n");
    ABORT ( status, COMPLEXAMC_EINPUT, COMPLEXAMC_MSGEINPUT );
  }

  /* check if that skyposition SSB+AMcoef were already buffered */
  if ( cfBuffer 
       && ( cfBuffer->multiDetStates == multiDetStates )
       && ( cfBuffer->Alpha == doppler->Alpha )
       && ( cfBuffer->Delta == doppler->Delta ) 
       && cfBuffer->multiSSB
       && cfBuffer->multiAMcoef )
    { /* yes ==> reuse */
      multiSSB = cfBuffer->multiSSB;
      multiAMcoef = cfBuffer -> multiAMcoef;
    }
  else 
    {
      SkyPosition skypos;
      /* use a dummy real detector state series for SSB times */
      MultiDetectorStateSeries multiRealDetStates
	= empty_MultiDetectorStateSeries;
      multiRealDetStates.length = multiDetStates.length;
      multiRealDetStates.startTime = multiDetStates.startTime;
      multiRealDetStates.Tspan = multiDetStates.Tspan;
      skypos.system =   COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;
      /* compute new AM-coefficients and SSB-times */
      TRY ( LALGetMultiSSBtimes ( status->statusPtr, &multiSSB, multiRealDetStates, skypos, doppler->refTime, params->SSBprec ), status );

      LALGetMultiCmplxAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, skypos );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      /* noise-weigh Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	LALPrintError("\nXLALWeighMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPLEXAMC_EXLAL, COMPLEXAMC_MSGEXLAL );
      }

      /* store these in buffer if available */
      if ( cfBuffer )
	{
	  XLALEmptyComputeFBuffer ( *cfBuffer );
	  cfBuffer->multiSSB = multiSSB;
	  cfBuffer->multiAMcoef = multiAMcoef;
	  cfBuffer->Alpha = doppler->Alpha;
	  cfBuffer->Delta = doppler->Delta;
	  cfBuffer->multiDetStates = multiDetStates ;
	} /* if cfBuffer */

    } /* if no buffer, different skypos or different detStates */

  Ad = multiAMcoef->Mmunu.Ad;
  Bd = multiAMcoef->Mmunu.Bd;
  Cd = multiAMcoef->Mmunu.Cd;
  Dd_inv = 1.0 / (Ad * Bd - Cd * Cd );

  /* ----- loop over detectors and compute all detector-specific quantities ----- */
  for ( X=0; X < numDetectors; X ++)
    {
      Fcomponents FcX = empty_Fcomponents;	/* for detector-specific FaX, FbX */

      if ( params->upsampling > 1) 
	{
	  if ( XLALComputeFaFbXavie (&FcX, multiSFTs->data[X], doppler->fkdot, multiSSB->data[X], multiAMcoef->data[X], params) != 0)
	    {
	      LALPrintError ("\nXALComputeFaFbXavie() failed\n");
	      ABORT ( status, COMPLEXAMC_EXLAL, COMPLEXAMC_MSGEXLAL );
	    }
	}
      else
	{
	  if ( XLALComputeFaFb (&FcX, multiSFTs->data[X], doppler->fkdot, multiSSB->data[X], multiAMcoef->data[X], params) != 0)
	    {
	      LALPrintError ("\nXALComputeFaFb() failed\n");
	      ABORT ( status, COMPLEXAMC_EXLAL, COMPLEXAMC_MSGEXLAL );
	    }
	}

#ifndef LAL_NDEBUG
      if ( !finite(FcX.Fa.re) || !finite(FcX.Fa.im) || !finite(FcX.Fb.re) || !finite(FcX.Fb.im) ) {
	LALPrintError("XLALComputeFaFb() returned non-finite: Fa=(%f,%f), Fb=(%f,%f)\n", 
		      FcX.Fa.re, FcX.Fa.im, FcX.Fb.re, FcX.Fb.im );
	ABORT (status,  COMPLEXAMC_EIEEE,  COMPLEXAMC_MSGEIEEE);
      }
#endif
 		 
      /* Fa = sum_X Fa_X */
      retF.Fa.re += FcX.Fa.re;
      retF.Fa.im += FcX.Fa.im;

      /* Fb = sum_X Fb_X */ 		  
      retF.Fb.re += FcX.Fb.re;
      retF.Fb.im += FcX.Fb.im;
  		  
    } /* for  X < numDetectors */
 
  /* ----- compute final Fstatistic-value ----- */

  /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
   * therefore there is a factor of 2 difference with respect to the equations in JKS, which 
   * where based on the single-sided PSD.
   */ 
 		       
  retF.F = Dd_inv * (  Bd * (SQ(retF.Fa.re) + SQ(retF.Fa.im) ) 
                     + Ad * ( SQ(retF.Fb.re) + SQ(retF.Fb.im) )
                     - 2.0 * Cd *( retF.Fa.re * retF.Fb.re + retF.Fa.im * retF.Fb.im )  
		   );

  (*Fstat) = retF;

  /* free memory if no buffer was available */
  if ( !cfBuffer )
    {
      XLALDestroyMultiSSBtimes ( multiSSB );
      XLALDestroyMultiCmplxAMCoeffs ( multiAMcoef );
    } /* if !cfBuffer */

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* CmplxComputeFStat() */

/** Compute the 'amplitude coefficients' \f$a(t)\sin\zeta\f$,
 * \f$b(t)\sin\zeta\f$ as defined in \ref JKS98 for a series of
 * timestamps.
 * 
 * The input consists of the DetectorState-timeseries, which contains
 * the detector-info and the LMST's corresponding to the different times.
 * 
 * In order to allow re-using the output-structure AMCoeffs for subsequent
 * calls, we require the REAL4Vectors a and b to be allocated already and 
 * to have the same length as the DetectoStates-timeseries.
 *
 * \note This is an alternative implementation to both LALComputeAM()
 * and LALGetAMCoeffs(), which uses the geometrical definition of
 * \f$a\sin\zeta\f$ and \f$b\sin\zeta\f$ as detector response
 * coefficients in a preferred polarization basis.  (It is thereby
 * more general than the JKS expressions and could be used e.g., with
 * the response tensor of a bar detector with no further modification
 * needed.)
 *
 * \note THIS VERSION IS FOR THE CASE WHERE THE DETECTOR TENSOR IS COMPLEX
 * AND THE DESCRIPTION NEEDS TO BE MODIFIED!
 *
 */ 
void
LALGetCmplxAMCoeffs(LALStatus *status,
		    CmplxAMCoeffs *coeffs,			/**< [out] amplitude-coeffs {a(f_0,t_i), b(f_0,t_i)} */
		    const CmplxDetectorStateSeries *DetectorStates,	/**< timeseries of detector states (note these may depend on f_0 and the sky position) */
		    SkyPosition skypos			/**< {alpha,delta} of the source */
	       )
{
  REAL4 delta, alpha;
  REAL4 sin1delta, cos1delta;
  REAL4 sin1alpha, cos1alpha;

  REAL4 xi1, xi2;
  REAL4 eta1, eta2, eta3;
  REAL4 norm;
  UINT4 i, numSteps;

  INITSTATUS (status, "LALGetCmplxAMCoeffs", COMPLEXAMC);

  /*---------- check input ---------- */
  ASSERT ( DetectorStates, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL);

  numSteps = DetectorStates->length;

  /* require the coeffients-vectors to be allocated and consistent with timestamps */
  ASSERT ( coeffs, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL);
  ASSERT ( coeffs->a && coeffs->b, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL);
  ASSERT ( (coeffs->a->length == numSteps) && (coeffs->b->length == numSteps), status,
	   COMPLEXAMC_EINPUT,  COMPLEXAMC_MSGEINPUT);

  /* require sky-pos to be in equatorial coordinates */
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, 
	   SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );

  /*---------- We write components of xi and eta vectors in SSB-fixed coords */
  alpha = skypos.longitude;
  delta = skypos.latitude;

  sin_cos_LUT (&sin1delta, &cos1delta, delta );
  sin_cos_LUT (&sin1alpha, &cos1alpha, alpha );
  xi1 = - sin1alpha;
  xi2 =  cos1alpha;
  eta1 = sin1delta * cos1alpha;
  eta2 = sin1delta * sin1alpha;
  eta3 = - cos1delta;

  /*---------- Compute the a(f_0,t_i) and b(f_0,t_i) ---------- */
  coeffs->A = 0;
  coeffs->B = 0;
  coeffs->C = 0;
  coeffs->E = 0;
  coeffs->D = 0;
  for ( i=0; i < numSteps; i++ )
    {
      COMPLEX8 ai, bi;

      CmplxDetectorTensor *d = &(DetectorStates->data[i].detT);
      
      ai.re = d->d11.re * ( xi1 * xi1 - eta1 * eta1 )
	+ 2 * d->d12.re * ( xi1*xi2 - eta1*eta2 )
	- 2 * d->d13.re *             eta1 * eta3
	+     d->d22.re * ( xi2*xi2 - eta2*eta2 )
	- 2 * d->d23.re *             eta2 * eta3
	-     d->d33.re *             eta3*eta3;

      ai.im = d->d11.im * ( xi1 * xi1 - eta1 * eta1 )
	+ 2 * d->d12.im * ( xi1*xi2 - eta1*eta2 )
	- 2 * d->d13.im *             eta1 * eta3
	+     d->d22.im * ( xi2*xi2 - eta2*eta2 )
	- 2 * d->d23.im *             eta2 * eta3
	-     d->d33.im *             eta3*eta3;

      bi.re = d->d11.re * 2 * xi1 * eta1
	+ 2 * d->d12.re *   ( xi1 * eta2 + xi2 * eta1 )
	+ 2 * d->d13.re *     xi1 * eta3
	+     d->d22.re * 2 * xi2 * eta2
	+ 2 * d->d23.re *     xi2 * eta3;

      bi.im = d->d11.im * 2 * xi1 * eta1
	+ 2 * d->d12.im *   ( xi1 * eta2 + xi2 * eta1 )
	+ 2 * d->d13.im *     xi1 * eta3
	+     d->d22.im * 2 * xi2 * eta2
	+ 2 * d->d23.im *     xi2 * eta3;

      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;

      /* sum A, B, C on the fly */
      coeffs->A += ai.re * ai.re + ai.im * ai.im;
      coeffs->B += bi.re * bi.re + bi.im * bi.im;
      coeffs->C += ai.re * bi.re + ai.im * bi.im;
      coeffs->E += ai.re * bi.im - ai.im * bi.re;

    } /* for i < numSteps */

  /* finish calculation of A,B,C,E,D */
  norm = 2.0f / numSteps;
  coeffs->A *= norm;
  coeffs->B *= norm;
  coeffs->C *= norm;
  coeffs->E *= norm;

  coeffs->D = coeffs->A * coeffs->B - coeffs->C * coeffs->C
    - coeffs->E * coeffs->E;

  RETURN(status);

} /* LALGetCmplxAMCoeffs() */

/** Multi-IFO version of LALGetCmplxAMCoeffs(). 
 * Get all antenna-pattern coefficients for all input detector-series.
 *
 * NOTE: contrary to LALGetCmplxAMCoeffs(), this functions *allocates* the output-vector,
 * use XLALDestroyMultiCmplxAMCoeffs() to free this.
 */
void
LALGetMultiCmplxAMCoeffs (LALStatus *status, 
		     MultiCmplxAMCoeffs **multiAMcoef,	/**< [out] AM-coefficients for all input detector-state series */
		     const MultiCmplxDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
		     SkyPosition skypos			/**< source sky-position [in equatorial coords!] */
		     )
{
  UINT4 X, numDetectors;
  MultiCmplxAMCoeffs *ret = NULL;

  INITSTATUS( status, "LALGetMultiCmplxAMCoeffs", COMPLEXAMC);
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT (multiDetStates, status,COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL);
  ASSERT (multiDetStates->length, status,COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL);
  ASSERT (multiAMcoef, status,COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL);
  ASSERT ( *multiAMcoef == NULL, status,COMPLEXAMC_ENONULL, COMPLEXAMC_MSGENONULL);
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, COMPLEXAMC_EINPUT, COMPLEXAMC_MSGEINPUT );

  numDetectors = multiDetStates->length;

  if ( ( ret = LALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    ABORT (status, COMPLEXAMC_EMEM, COMPLEXAMC_MSGEMEM);    
  }
  ret->length = numDetectors;
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, COMPLEXAMC_EMEM, COMPLEXAMC_MSGEMEM);
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      CmplxAMCoeffs *amcoeX = NULL;
      UINT4 numStepsX = multiDetStates->data[X]->length;

      ret->data[X] = LALCalloc ( 1, sizeof ( *(ret->data[X]) ) );
      amcoeX = ret->data[X];
      amcoeX->a = XLALCreateCOMPLEX8Vector ( numStepsX );
      if ( (amcoeX->b = XLALCreateCOMPLEX8Vector ( numStepsX )) == NULL ) {
	LALPrintError ("\nOut of memory!\n\n");
	goto failed;
      }

      LALGetCmplxAMCoeffs (status->statusPtr, amcoeX, multiDetStates->data[X], skypos );
      if ( status->statusPtr->statusCode ) 
	{
	  LALPrintError ( "\nCall to LALGetCmplxAMCoeffs() has failed ... \n\n");
	  goto failed;
	}
 
    } /* for X < numDetectors */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiCmplxAMCoeffs ( ret );
  ABORT ( status, -1, "LALGetMultiCmplxAMCoeffs() failed" );

 success:
  (*multiAMcoef) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetMultiCmplxAMCoeffs() */



/* ===== Object creation/destruction functions ===== */

/** Destroy a MultiCmplxAMCoeffs structure. 
 * Note, this is "NULL-robust" in the sense that it will not crash 
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs 
 */
void
XLALDestroyMultiCmplxAMCoeffs ( MultiCmplxAMCoeffs *multiAMcoef )
{
  UINT4 X;
  CmplxAMCoeffs *tmp;

  if ( ! multiAMcoef )
    return;

  if ( multiAMcoef->data )
    {
      for ( X=0; X < multiAMcoef->length; X ++ ) 
	{
	  if ( (tmp = multiAMcoef->data[X]) != NULL )
	    {
	      if ( tmp->a )
		XLALDestroyCOMPLEX8Vector ( tmp->a );
	      if ( tmp->b )
		XLALDestroyCOMPLEX8Vector ( tmp->b );
	      LALFree ( tmp );
	    } /* if multiAMcoef->data[X] */
	} /* for X < numDetectors */
      LALFree ( multiAMcoef->data );
    }
  LALFree ( multiAMcoef );

  return;

} /* XLALDestroyMultiCmplxAMCoeffs() */


/** Multiply AM-coeffs \f$a_{X\alpha}, b_{X\alpha}\f$ by weights \f$\sqrt(w_{X\alpha})\f$ and 
 * compute the resulting \f$A_d, B_d, C_d, E_d\f$ by simply *SUMMING* them, i.e.
 * \f$A_d \equiv \sum_{X,\alpha} \sqrt{w_{X\alpha} a_{X\alpha}^2\f$ etc.
 * 
 * NOTE: this function modifies the CmplxAMCoeffs *in place* !
 * NOTE2: if the weights = NULL, we assume unit-weights.
 */
int
XLALWeighMultiCmplxAMCoeffs (  MultiCmplxAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights )
{
  UINT4 numDetectors, X;
  REAL8 Ad, Bd, Cd, Ed;
  UINT4 alpha;

  if ( !multiAMcoef )
    XLAL_ERROR( "XLALWeighMultiCmplxAMCoeffs", XLAL_EINVAL );

  numDetectors = multiAMcoef->length;

  if ( multiWeights && ( multiWeights->length != numDetectors ) )
    {
      LALPrintError("\nmultiWeights must have same length as mulitAMcoef!\n\n");
      XLAL_ERROR( "XLALWeighMultiCmplxAMCoeffs", XLAL_EINVAL );
    }
  
  /* noise-weight Antenna-patterns and compute A,B,C,E */
  Ad = Bd = Cd = Ed = 0;

  if ( multiWeights  )
    {
      for ( X=0; X < numDetectors; X ++)
	{
	  CmplxAMCoeffs *amcoeX = multiAMcoef->data[X];
	  UINT4 numSteps = amcoeX->a->length;
	  
	  REAL8Vector *weightsX = multiWeights->data[X];;
	  if ( weightsX->length != numSteps ) 
	    {
	      LALPrintError("\nmultiWeights must have same length as mulitAMcoef!\n\n");
	      XLAL_ERROR( "XLALWeighMultiCmplxAMCoeffs", XLAL_EINVAL );
	    }
	  
	  for(alpha = 0; alpha < numSteps; alpha++)
	    {
	      REAL8 Sqwi = sqrt ( weightsX->data[alpha] );
	      COMPLEX16 ahat;
	      COMPLEX16 bhat;
	      ahat.re = Sqwi * amcoeX->a->data[alpha].re;
	      ahat.im = Sqwi * amcoeX->a->data[alpha].im;
	      bhat.re= Sqwi * amcoeX->b->data[alpha].re;
	      bhat.im= Sqwi * amcoeX->b->data[alpha].im;
	      
	      /* *replace* original a(t), b(t) by noise-weighed version! */
	      amcoeX->a->data[alpha].re = ahat.re;
	      amcoeX->a->data[alpha].im = ahat.im;
	      amcoeX->b->data[alpha].re = bhat.re;
	      amcoeX->b->data[alpha].im = bhat.im;
	      
	      /* sum A, B, C, E on the fly */
	      Ad += ahat.re * ahat.re + ahat.im * ahat.im;
	      Bd += bhat.re * bhat.re + bhat.im * bhat.im;
	      Cd += ahat.re * bhat.re + ahat.im * bhat.im;
	      Ed += ahat.re * bhat.re - ahat.im * bhat.im;
	    } /* for alpha < numSFTsX */
	} /* for X < numDetectors */
      multiAMcoef->Mmunu.Sinv_Tsft = multiWeights->Sinv_Tsft;
    }
  else /* if no noise-weights: simply add to get A,B,C */
    {
      for ( X=0; X < numDetectors; X ++)
	{
	  CmplxAMCoeffs *amcoeX = multiAMcoef->data[X];
	  UINT4 numSteps = amcoeX->a->length;
	  
	  for(alpha = 0; alpha < numSteps; alpha++)
	    {
	      COMPLEX16 ahat;
	      COMPLEX16 bhat;
	      ahat.re = amcoeX->a->data[alpha].re;
	      ahat.im = amcoeX->a->data[alpha].im;
	      bhat.re = amcoeX->b->data[alpha].re;
	      bhat.im = amcoeX->b->data[alpha].im;
	    
	      /* sum A, B, C, E on the fly */
	      Ad += ahat.re * ahat.re + ahat.im * ahat.im;
	      Bd += bhat.re * bhat.re + bhat.im * bhat.im;
	      Cd += ahat.re * bhat.re + ahat.im * bhat.im;
	      Ed += ahat.re * bhat.re - ahat.im * bhat.im;
	    } /* for alpha < numSFTsX */
	} /* for X < numDetectors */

    } /* if multiWeights == NULL */

  multiAMcoef->Mmunu.Ad = Ad;
  multiAMcoef->Mmunu.Bd = Bd;
  multiAMcoef->Mmunu.Cd = Cd;


  return XLAL_SUCCESS;

} /* XLALWeighMultiCmplxAMCoefs() */
