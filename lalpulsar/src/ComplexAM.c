/*
 * Copyright (C) 2007 John T. Whelan, Reinhard Prix
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
 * \author J. T. Whelan, Reinhard Prix
 * \ingroup pulsarTODO
 * \file
 * \brief
 * Functions related to F-statistic calculation when the AM coefficients are complex.
 *
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

/* GSL includes */
#include <lal/LALGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

/* LAL includes */
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/AVFactories.h>
#include <lal/ComplexAM.h>
#include <lal/LISAspecifics.h>
#include <lal/CWFastMath.h>

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * Compute the 'amplitude coefficients' \f$a(t)\sin\zeta\f$,
 * \f$b(t)\sin\zeta\f$ as defined in \cite JKS98 for a series of
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
LALGetCmplxAMCoeffs(LALStatus *status,				/**< pointer to LALStatus structure */
		    CmplxAMCoeffs *coeffs,			/**< [out] amplitude-coeffs {a(f_0,t_i), b(f_0,t_i)} */
		    const DetectorStateSeries *DetectorStates,	/**< timeseries of detector states */
		    const FreqSkypos_t *freq_skypos		/**< Frequency and skyposition information */
		    )
{
  UINT4 i, numSteps;
  CHAR channelNum;

  INITSTATUS(status);

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

  /*---------- Compute the a(f_0,t_i) and b(f_0,t_i) ---------- */
  for ( i=0; i < numSteps; i++ )
    {
      COMPLEX8 ai, bi;
      CmplxDetectorTensor d;

      if ( XLALgetLISADetectorTensorRAA (&d, DetectorStates->data[i].detArms, channelNum, freq_skypos ) != 0 ) {
	XLALPrintError ( "\nXLALgetCmplxLISADetectorTensor() failed ... errno = %d\n\n", xlalErrno );
	ABORT ( status, COMPLEXAMC_EXLAL, COMPLEXAMC_MSGEXLAL );
      }

      ai = crectf( XLALContractSymmTensor3s ( &d.re, &(freq_skypos->ePlus) ), XLALContractSymmTensor3s ( &d.im, &(freq_skypos->ePlus) ) );

      bi = crectf( XLALContractSymmTensor3s ( &d.re, &(freq_skypos->eCross) ), XLALContractSymmTensor3s ( &d.im, &(freq_skypos->eCross) ) );

      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;

    } /* for i < numSteps */

  RETURN(status);

} /* LALGetCmplxAMCoeffs() */

/**
 * Multi-IFO version of LALGetCmplxAMCoeffs().
 * Get all antenna-pattern coefficients for all input detector-series.
 *
 * NOTE: contrary to LALGetCmplxAMCoeffs(), this functions *allocates* the output-vector,
 * use XLALDestroyMultiCmplxAMCoeffs() to free this.
 */
void
LALGetMultiCmplxAMCoeffs (LALStatus *status,				/**< pointer to LALStatus structure */
		     MultiCmplxAMCoeffs **multiAMcoef,			/**< [out] AM-coefficients for all input detector-state series */
		     const MultiDetectorStateSeries *multiDetStates, 	/**< [in] detector-states at timestamps t_i */
		     PulsarDopplerParams doppler		     	/**< source sky-position [in equatorial coords!], freq etc. */
		     )
{
  UINT4 X, numDetectors;
  MultiCmplxAMCoeffs *ret = NULL;
  REAL4 sin1Delta, cos1Delta;
  REAL4 sin1Alpha, cos1Alpha;
  FreqSkypos_t freq_skypos;

  INITSTATUS(status);
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

  if( XLALSinCosLUT (&sin1Delta, &cos1Delta, doppler.Delta ) != XLAL_SUCCESS )
    ABORT( status->statusPtr, LAL_EXLAL, "XLALSinCosLUT (&sin1Delta, &cos1Delta, doppler.Delta ) failed" );

  if( XLALSinCosLUT (&sin1Alpha, &cos1Alpha, doppler.Alpha ) != XLAL_SUCCESS )
    ABORT( status->statusPtr, LAL_EXLAL, "XLALSinCosLUT (&sin1Alpha, &cos1Alpha, doppler.Alpha ) failed" );


  freq_skypos.skyposV[0] = cos1Delta * cos1Alpha;
  freq_skypos.skyposV[1] = cos1Delta * sin1Alpha;
  freq_skypos.skyposV[2] = sin1Delta;

  /*---------- compute components of xi and eta vectors in SSB-fixed coords */
  {
    REAL4 xi[3];
    REAL4 eta[3];
    SymmTensor3 etaT, xiT;

    xi[0] = - sin1Alpha;
    xi[1] =   cos1Alpha;
    xi[2] =   0.0f;

    eta[0] = sin1Delta * cos1Alpha;
    eta[1] = sin1Delta * sin1Alpha;
    eta[2] = - cos1Delta;

    /* compute e+ */
    XLALTensorSquareVector3 ( &xiT,   xi );	  /* the tensor xi x xi */
    XLALTensorSquareVector3 ( &etaT, eta );	  /* the tensor eta x eta */
    XLALSubtractSymmTensor3s ( &(freq_skypos.ePlus), &xiT, &etaT );	/* e+ = xi x xi - eta x eta */

    /* compute ex */
    XLALSymmetricTensorProduct3 ( &(freq_skypos.eCross), xi, eta );	/* ex = xi x eta + eta x xi */
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      CmplxAMCoeffs *amcoeX = NULL;
      UINT4 numStepsX = multiDetStates->data[X]->length;

      ret->data[X] = LALCalloc ( 1, sizeof ( *(ret->data[X]) ) );
      amcoeX = ret->data[X];
      amcoeX->a = XLALCreateCOMPLEX8Vector ( numStepsX );
      if ( (amcoeX->b = XLALCreateCOMPLEX8Vector ( numStepsX )) == NULL ) {
	XLALPrintError ("\nOut of memory!\n\n");
	goto failed;
      }

      freq_skypos.Freq = doppler.fkdot[0];

      LALGetCmplxAMCoeffs (status->statusPtr, amcoeX, multiDetStates->data[X], &freq_skypos );
      if ( status->statusPtr->statusCode )
	{
	  XLALPrintError ( "\nCall to LALGetCmplxAMCoeffs() has failed ... \n\n");
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

/**
 * Destroy a MultiCmplxAMCoeffs structure.
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


/**
 * Multiply AM-coeffs \f$a_{X\alpha}, b_{X\alpha}\f$ by weights \f$\sqrt(w_{X\alpha})\f$ and
 * compute the resulting \f$\widehat{A}, \widehat{B}, \widehat{C}, \widehat{E}\f$ by simply *SUMMING* them, i.e.
 * \f$\widehat{A} \equiv \sum_{X,\alpha} w_{X\alpha} a_{X\alpha}^2\f$ etc.
 *
 * NOTE: this function modifies the CmplxAMCoeffs *in place* !
 * NOTE2: if the weights = NULL, we assume unit-weights.
 */
int
XLALWeightMultiCmplxAMCoeffs (  MultiCmplxAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights )
{
  UINT4 numDetectors, X;
  REAL8 Ad, Bd, Cd, Ed;
  UINT4 alpha;

  if ( !multiAMcoef )
    XLAL_ERROR( XLAL_EINVAL );

  numDetectors = multiAMcoef->length;

  if ( multiWeights && ( multiWeights->length != numDetectors ) )
    {
      XLALPrintError("\nmultiWeights must have same length as mulitAMcoef!\n\n");
      XLAL_ERROR( XLAL_EINVAL );
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
	      XLALPrintError("\nmultiWeights must have same length as mulitAMcoef!\n\n");
	      XLAL_ERROR( XLAL_EINVAL );
	    }

	  for(alpha = 0; alpha < numSteps; alpha++)
	    {
	      REAL8 Sqwi = sqrt ( weightsX->data[alpha] );
	      COMPLEX8 ahat;
	      COMPLEX8 bhat;
	      ahat = crect( Sqwi * crealf(amcoeX->a->data[alpha]), Sqwi * cimagf(amcoeX->a->data[alpha]) );
	      bhat = crect( Sqwi * crealf(amcoeX->b->data[alpha]), Sqwi * cimagf(amcoeX->b->data[alpha]) );

	      /* *replace* original a(t), b(t) by noise-weighed version! */
	      amcoeX->a->data[alpha] = crectf( creal(ahat), cimag(ahat) );
	      amcoeX->b->data[alpha] = crectf( creal(bhat), cimag(bhat) );

	      /* sum A, B, C, E on the fly */
	      Ad += creal(ahat) * creal(ahat) + cimag(ahat) * cimag(ahat);
	      Bd += creal(bhat) * creal(bhat) + cimag(bhat) * cimag(bhat);
	      Cd += creal(ahat) * creal(bhat) + cimag(ahat) * cimag(bhat);
	      Ed += creal(ahat) * cimag(bhat) - cimag(ahat) * creal(bhat);
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
	      COMPLEX8 ahat;
	      COMPLEX8 bhat;
	      ahat = crect( crealf(amcoeX->a->data[alpha]), cimagf(amcoeX->a->data[alpha]) );
	      bhat = crect( crealf(amcoeX->b->data[alpha]), cimagf(amcoeX->b->data[alpha]) );

	      /* sum A, B, C, E on the fly */
	      Ad += creal(ahat) * creal(ahat) + cimag(ahat) * cimag(ahat);
	      Bd += creal(bhat) * creal(bhat) + cimag(bhat) * cimag(bhat);
	      Cd += creal(ahat) * creal(bhat) + cimag(ahat) * cimag(bhat);
	      Ed += creal(ahat) * cimag(bhat) - cimag(ahat) * creal(bhat);
	    } /* for alpha < numSFTsX */
	} /* for X < numDetectors */

    } /* if multiWeights == NULL */

  multiAMcoef->Mmunu.Ad = Ad;
  multiAMcoef->Mmunu.Bd = Bd;
  multiAMcoef->Mmunu.Cd = Cd;
  multiAMcoef->Mmunu.Ed = Ed;
  multiAMcoef->Mmunu.Dd = Ad * Bd - Cd * Cd - Ed * Ed;

  return XLAL_SUCCESS;

} /* XLALWeightMultiCmplxAMCoefs() */


#if 0
/* Revamped version of XLALComputeFaFb() for the case where a and b
 * are complex.
 * Compute JKS's Fa and Fb, which are ingredients for
 * calculating the F-statistic.
 */
static int
XLALComputeFaFbCmplx ( Fcomponents *FaFb,               /* [out] Fa,Fb (and possibly atoms) returned */
                  const SFTVector *sfts,                /* [in] input SFTs */
                  const PulsarSpins fkdot,              /* [in] frequency and derivatives fkdot = d^kf/dt^k */
                  const SSBtimes *tSSB,                 /* [in] SSB timing series for particular sky-direction */
                  const CmplxAMCoeffs *amcoe,           /* [in] antenna-pattern coefficients for this sky-direction */
                  const ComputeFParams *params)         /* addition computational params */
{
  UINT4 alpha;                  /* loop index over SFTs */
  UINT4 spdnOrder;              /* maximal spindown-orders */
  UINT4 numSFTs;                /* number of SFTs (M in the Notes) */
  COMPLEX8 Fa, Fb;
  REAL8 Tsft;                   /* length of SFTs in seconds */
  INT4 freqIndex0;              /* index of first frequency-bin in SFTs */
  INT4 freqIndex1;              /* index of last frequency-bin in SFTs */

  COMPLEX8 *a_al, *b_al;        /* pointer to alpha-arrays over a and b */
  REAL8 *DeltaT_al, *Tdot_al;   /* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;              /* SFT alpha  */
  UINT4 Dterms = params->Dterms;

  REAL8 norm = OOTWOPI;

  /* ----- check validity of input */
#ifndef LAL_NDEBUG
  if ( !FaFb ) {
    XLALPrintError ("\nOutput-pointer is NULL !\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    XLALPrintError ("\nInput SFTs are NULL!\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !tSSB || !tSSB->DeltaT || !tSSB->Tdot || !amcoe || !amcoe->a || !amcoe->b || !params)
    {
      XLALPrintError ("\nIllegal NULL in input !\n\n");
      XLAL_ERROR ( XLAL_EINVAL);
    }

  if ( PULSAR_MAX_SPINS > LAL_FACT_MAX )
    {
      XLALPrintError ("\nInverse factorials table only up to order s=%d, can't handle %d spin-order\n\n",
                     LAL_FACT_MAX, PULSAR_MAX_SPINS - 1 );
      XLAL_ERROR ( XLAL_EINVAL);
    }

  if ( params->returnAtoms )
    {
      XLALPrintError ("%s: using the option 'returnAtoms' is not supported in this function!\n", __func__ );
      XLAL_ERROR ( XLAL_EINVAL);
    }
#endif

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;
  {
    REAL8 dFreq = sfts->data[0].deltaF;
    freqIndex0 = (UINT4) ( sfts->data[0].f0 / dFreq + 0.5); /* lowest freqency-index */
    freqIndex1 = freqIndex0 + sfts->data[0].data->length;
  }

  /* find highest non-zero spindown-entry */
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot[spdnOrder] )
      break;

  Fa = 0.0f;
  Fb = 0.0f;

  a_al = amcoe->a->data;        /* point to beginning of alpha-arrays */
  b_al = amcoe->b->data;
  DeltaT_al = tSSB->DeltaT->data;
  Tdot_al = tSSB->Tdot->data;
  SFT_al = sfts->data;

  /* Loop over all SFTs  */
  for ( alpha = 0; alpha < numSFTs; alpha++ )
    {
      COMPLEX8 a_alpha, b_alpha;

      INT4 kstar;               /* central frequency-bin k* = round(xhat_alpha) */
      INT4 k0, k1;

      COMPLEX8 *Xalpha = SFT_al->data->data; /* pointer to current SFT-data */
      COMPLEX8 *Xalpha_l;       /* pointer to frequency-bin k in current SFT */
      REAL4 s_alpha=0, c_alpha=0;/* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
      REAL4 realQ, imagQ;       /* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realXP, imagXP;     /* Re/Im of sum_k X_ak * P_ak */
      REAL4 realQXP, imagQXP;   /* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha, kappa_max, kappa_star;

      /* ----- calculate kappa_max and lambda_alpha */
      {
        UINT4 s;                /* loop-index over spindown-order */
        REAL8 phi_alpha, Dphi_alpha, DT_al;
        REAL8 Tas;      /* temporary variable to calculate (DeltaT_alpha)^s */

        /* init for s=0 */
        phi_alpha = 0.0;
        Dphi_alpha = 0.0;
        DT_al = (*DeltaT_al);
        Tas = 1.0;              /* DeltaT_alpha ^ 0 */

        for (s=0; s <= spdnOrder; s++)
          {
            REAL8 fsdot = fkdot[s];
            Dphi_alpha += fsdot * Tas * LAL_FACT_INV[s];        /* here: DT^s/s! */
            Tas *= DT_al;                               /* now: DT^(s+1) */
            phi_alpha += fsdot * Tas * LAL_FACT_INV[s+1];
          } /* for s <= spdnOrder */

        /* Step 3: apply global factors to complete Dphi_alpha */
        Dphi_alpha *= Tsft * (*Tdot_al);                /* guaranteed > 0 ! */

        lambda_alpha = phi_alpha - 0.5 * Dphi_alpha;

        /* real- and imaginary part of e^{-i 2 pi lambda_alpha } */
        if ( XLALSinCos2PiLUT ( &imagQ, &realQ, - lambda_alpha ) ) {
          XLAL_ERROR ( XLAL_EFUNC);
        }

        kstar = (INT4) (Dphi_alpha);    /* k* = floor(Dphi_alpha) for positive Dphi */
        kappa_star = Dphi_alpha - 1.0 * kstar;  /* remainder of Dphi_alpha: >= 0 ! */
        kappa_max = kappa_star + 1.0 * Dterms - 1.0;

        /* ----- check that required frequency-bins are found in the SFTs ----- */
        k0 = kstar - Dterms + 1;
        k1 = k0 + 2 * Dterms - 1;
        if ( (k0 < freqIndex0) || (k1 > freqIndex1) )
          {
            XLALPrintError ("Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n\n",
                           k0, k1, freqIndex0, freqIndex1 );
            XLAL_ERROR(XLAL_EDOM);
          }

      } /* compute kappa_star, lambda_alpha */

      /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
       * the trig-functions need to be calculated only once!
       * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the
       * closest to zero and will pose no numerical difficulties !
       */
      XLAL_CHECK( XLALSinCos2PiLUT ( &s_alpha, &c_alpha, kappa_star ) == XLAL_SUCCESS, XLAL_EFUNC );
      c_alpha -= 1.0f;

      /* ---------- calculate the (truncated to Dterms) sum over k ---------- */

      /* ---------- ATTENTION: this the "hot-loop", which will be
       * executed many millions of times, so anything in here
       * has a HUGE impact on the whole performance of the code.
       *
       * DON'T touch *anything* in here unless you really know
       * what you're doing !!
       *------------------------------------------------------------
       */

      Xalpha_l = Xalpha + k0 - freqIndex0;  /* first frequency-bin in sum */

      realXP = 0;
      imagXP = 0;

      /* if no danger of denominator -> 0 */
      if ( ( kappa_star > LD_SMALL4 ) && (kappa_star < 1.0 - LD_SMALL4) )
        {
          /* improved hotloop algorithm by Fekete Akos:
           * take out repeated divisions into a single common denominator,
           * plus use extra cleverness to compute the nominator efficiently...
           */
          REAL4 Sn = crealf(*Xalpha_l);
          REAL4 Tn = cimagf(*Xalpha_l);
          REAL4 pn = kappa_max;
          REAL4 qn = pn;
          REAL4 U_alpha, V_alpha;

          /* recursion with 2*Dterms steps */
          UINT4 l;
          for ( l = 1; l < 2*Dterms; l ++ )
            {
              Xalpha_l ++;

              pn = pn - 1.0f;                   /* p_(n+1) */
              Sn = pn * Sn + qn * crealf(*Xalpha_l);    /* S_(n+1) */
              Tn = pn * Tn + qn * cimagf(*Xalpha_l);    /* T_(n+1) */
              qn *= pn;                         /* q_(n+1) */
            } /* for l <= 2*Dterms */

          U_alpha = Sn / qn;
          V_alpha = Tn / qn;

#ifndef LAL_NDEBUG
          if ( !isfinite(U_alpha) || !isfinite(V_alpha) || !isfinite(pn) || !isfinite(qn) || !isfinite(Sn) || !isfinite(Tn) ) {
            XLALPrintError("XLALComputeFaFbCmplx() returned non-finite: U_alpha=%f, V_alpha=%f, pn=%f, qn=%f, Sn=%f, Tn=%f\n",
                           U_alpha, V_alpha, pn, qn, Sn, Tn);
            XLAL_ERROR (XLAL_EFPINVAL);
          }
#endif

          realXP = s_alpha * U_alpha - c_alpha * V_alpha;
          imagXP = c_alpha * U_alpha + s_alpha * V_alpha;

        } /* if |remainder| > LD_SMALL4 */
      else
        { /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
          UINT4 ind0;
          if ( kappa_star <= LD_SMALL4 ) ind0 = Dterms - 1;
          else ind0 = Dterms;
          realXP = TWOPI_FLOAT * crealf(Xalpha_l[ind0]);
          imagXP = TWOPI_FLOAT * cimagf(Xalpha_l[ind0]);
        } /* if |remainder| <= LD_SMALL4 */

      realQXP = realQ * realXP - imagQ * imagXP;
      imagQXP = realQ * imagXP + imagQ * realXP;

      /* we're done: ==> combine these into Fa and Fb */
      a_alpha = (*a_al);
      b_alpha = (*b_al);

      /* Fa contains complex conjugate of a */
      Fa += crect( crealf(a_alpha) * realQXP + cimagf(a_alpha) * imagQXP, crealf(a_alpha) * imagQXP - cimagf(a_alpha) * realQXP );

      /* Fb contains complex conjugate of b */
      Fb += crect( crealf(b_alpha) * realQXP + cimagf(b_alpha) * imagQXP, crealf(b_alpha) * imagQXP - cimagf(b_alpha) * realQXP );

      /* advance pointers over alpha */
      a_al ++;
      b_al ++;
      DeltaT_al ++;
      Tdot_al ++;
      SFT_al ++;

    } /* for alpha < numSFTs */

  /* return result */
  FaFb->Fa = norm * Fa;
  FaFb->Fb = norm * Fb;

  return XLAL_SUCCESS;

} /* XLALComputeFaFbCmplx() */
#endif
