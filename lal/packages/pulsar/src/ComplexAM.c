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
