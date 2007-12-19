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
#include <lal/LISAspecifics.h>
#include <lal/ComplexAM.h>

NRCSID( COMPLEXAMC, "$Id$");

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*==================== FUNCTION DEFINITIONS ====================*/

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
		    const DetectorStateSeries *DetectorStates,	/**< timeseries of detector states */
		    PulsarDopplerParams doppler			/**< {alpha,delta} of the source */
	       )
{
  REAL4 sin1delta, cos1delta;
  REAL4 sin1alpha, cos1alpha;

  REAL4 xi[3];
  REAL4 eta[3];
  REAL4 norm;
  UINT4 i, numSteps;
  SymmTensor3 etaT, xiT, ePlus, eCross;

  CHAR channelNum;

  INITSTATUS (status, "LALGetCmplxAMCoeffs", COMPLEXAMC);

  /*---------- check input ---------- */
  ASSERT ( DetectorStates, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL);

  numSteps = DetectorStates->length;

  /* require the coeffients-vectors to be allocated and consistent with timestamps */
  ASSERT ( coeffs, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL);
  ASSERT ( coeffs->a && coeffs->b, status, COMPLEXAMC_ENULL, COMPLEXAMC_MSGENULL);
  ASSERT ( (coeffs->a->length == numSteps) && (coeffs->b->length == numSteps), status,
	   COMPLEXAMC_EINPUT,  COMPLEXAMC_MSGEINPUT);

  ASSERT ( DetectorStates->detector.frDetector.prefix[0] == 'Z', status,
	   COMPLEXAMC_ERAALISA, COMPLEXAMC_MSGERAALISA);

  /* need to know TDI channel number to calculate complex detector tensor */
  channelNum = DetectorStates->detector.frDetector.prefix[1];

  /*---------- We write components of xi and eta vectors in SSB-fixed coords */
  sin_cos_LUT (&sin1delta, &cos1delta, doppler.Delta );
  sin_cos_LUT (&sin1alpha, &cos1alpha, doppler.Alpha );
  xi[0] = - sin1alpha;
  xi[1] =   cos1alpha;
  xi[2] =   0.0f;

  eta[0] = sin1delta * cos1alpha;
  eta[1] = sin1delta * sin1alpha;
  eta[2] = - cos1delta;

  /* compute e+ */
  XLALTensorSquareVector3 ( &xiT,   xi );	  /* the tensor xi x xi */
  XLALTensorSquareVector3 ( &etaT, eta );	  /* the tensor eta x eta */
  XLALSubtractSymmTensor3s ( &ePlus, &xiT, &etaT );	/* e+ = xi x xi - eta x eta */

  /* compute ex */
  XLALSymmetricTensorProduct3 ( &eCross, xi, eta );	/* ex = xi x eta + eta x xi */


  /*---------- Compute the a(f_0,t_i) and b(f_0,t_i) ---------- */
  coeffs->A = 0;
  coeffs->B = 0;
  coeffs->C = 0;
  coeffs->E = 0;
  coeffs->D = 0;

  for ( i=0; i < numSteps; i++ )
    {
      COMPLEX8 ai, bi;
      CmplxDetectorTensor d;

      if ( XLALgetLISADetectorTensorRAA (&d, DetectorStates->data[i].detArms, doppler, channelNum) != 0 ) {
	LALPrintError ( "\nXLALgetCmplxLISADetectorTensor() failed ... errno = %d\n\n", xlalErrno );
	ABORT ( status, COMPLEXAMC_EXLAL, COMPLEXAMC_MSGEXLAL );
      }

      ai.re = XLALContractSymmTensor3s ( &d.re, &ePlus );
      ai.im = XLALContractSymmTensor3s ( &d.im, &ePlus );

      bi.re = XLALContractSymmTensor3s ( &d.re, &eCross );
      bi.im = XLALContractSymmTensor3s ( &d.im, &eCross );

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

  coeffs->D = coeffs->A * coeffs->B - coeffs->C * coeffs->C - coeffs->E * coeffs->E;

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
		     const MultiDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
		     PulsarDopplerParams doppler		     /**< source sky-position [in equatorial coords!], freq etc. */
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

      LALGetCmplxAMCoeffs (status->statusPtr, amcoeX, multiDetStates->data[X], doppler );
      if ( status->statusPtr->statusCode )
	{
	  LALPrintError ( "\nCall to LALGetCmplxAMCoeffs() has failed ... \n\n");
	  REPORTSTATUS ( status->statusPtr );
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
	      Ed += ahat.re * bhat.im - ahat.im * bhat.re;
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
	      Ed += ahat.re * bhat.im - ahat.im * bhat.re;
	    } /* for alpha < numSFTsX */
	} /* for X < numDetectors */

    } /* if multiWeights == NULL */

  multiAMcoef->Mmunu.Ad = Ad;
  multiAMcoef->Mmunu.Bd = Bd;
  multiAMcoef->Mmunu.Cd = Cd;
  multiAMcoef->Mmunu.Ed = Ed;

  return XLAL_SUCCESS;

} /* XLALWeighMultiCmplxAMCoefs() */
