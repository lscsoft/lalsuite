/*
 * Copyright (C) 2005 Reinhard Prix
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

/** \author R. Prix
 * \file 
 * \brief
 * Functions to calculate the so-called F-statistic for a given point in parameter-space, 
 * following the equations in \ref JKS98.
 *                                                                          
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include <lal/AVFactories.h>
#include "ComputeFstat.h"

NRCSID( COMPUTEFSTATC, "$Id$");

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)


#define LD_SMALL4       (1.0e-9)		/**< "small" number for REAL4*/
#define OOTWOPI         (1.0 / LAL_TWOPI)	/**< 1/2pi */

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */ 


/*----- Macros ----- */
/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

#define SQ(x) ( (x) * (x) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */
static const BarycenterInput empty_BarycenterInput;
static const LALStatus empty_status;
static const Fcomponents empty_Fcomponents;

/*---------- Global variables ----------*/
#define NUM_FACT 6
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0) };

/*---------- internal prototypes ----------*/

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
 */
void
ComputeFStat ( LALStatus *status, 
	       Fcomponents *Fstat, 		/**< [out] Fstatistic + Fa, Fb */
	       const CWParamSpacePoint *psPoint,/**< parameter-space point to compute F for */
	       const MultiSFTVector *multiSFTs, /**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
	       const MultiNoiseWeights *multiWeights,	/**< noise-weights of all SFTs */
	       const MultiDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
	       const ComputeFParams *params,	/**< addition computational params */
	       ComputeFBuffer *cfBuffer		/**< CF-internal buffering structure */
	       )
{
  Fcomponents retF = empty_Fcomponents;
  UINT4 X, numDetectors;	
  MultiSSBtimes *multiSSB = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  REAL8 A, B, C, Dinv;

  INITSTATUS( status, "ComputeFStat", COMPUTEFSTATC );
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT ( Fstat, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiSFTs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiWeights, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( psPoint, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( psPoint->fkdot, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  numDetectors = multiSFTs->length;
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  ASSERT ( multiWeights->length == numDetectors , status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );

  if ( psPoint->binary ) {
    LALPrintError ("\nSorry, binary-pulsar search not yet implemented in LALComputeFStat()\n\n");
    ABORT ( status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  }

  /* check if that skyposition SSB+AMcoef were already buffered */
  if ( cfBuffer 
       && ( cfBuffer->skypos.longitude == psPoint->skypos.longitude)
       && ( cfBuffer->skypos.latitude == psPoint->skypos.latitude) 
       && cfBuffer->multiSSB
       && cfBuffer->multiAMcoef )
    { /* yes ==> reuse */
      multiSSB = cfBuffer->multiSSB;
      multiAMcoef = cfBuffer -> multiAMcoef;
      A = cfBuffer->A;
      B = cfBuffer->B;
      C = cfBuffer->C;
      Dinv = cfBuffer->Dinv;
    }
  else 
    {
      /* compute new AM-coefficients and SSB-times */
      TRY ( LALGetMultiSSBtimes ( status->statusPtr, &multiSSB, multiDetStates, psPoint->skypos, psPoint->refTime, params->SSBprec ), status );

      LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, psPoint->skypos );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      /* noise-weight Antenna-patterns and compute A,B,C */
      A = B = C = Dinv = 0.0;
      for ( X=0; X < numDetectors; X ++)
	{
	  UINT4 alpha;
	  UINT4 numSFTsX = multiSFTs->data[X]->length;
	  AMCoeffs *amcoeX = multiAMcoef->data[X];
	  REAL8Vector *weightsX = multiWeights->data[X];

	  for(alpha = 0; alpha < numSFTsX; alpha++)
	    {
	      REAL8 Sqwi = sqrt ( weightsX->data[alpha] );
	      REAL8 ahat = Sqwi * amcoeX->a->data[alpha] ;
	      REAL8 bhat = Sqwi * amcoeX->b->data[alpha] ;
	  
	      /* *replace* original a(t), b(t) by noise-weighed version! */
	      amcoeX->a->data[alpha] = ahat;
	      amcoeX->b->data[alpha] = bhat;
 
	      /* sum A, B, C on the fly */
	      A += ahat * ahat;
	      B += bhat * bhat;
	      C += ahat * bhat;
	    } /* for alpha < numSFTsX */
	} /* for X < numDetectors */
      Dinv = 1.0 / (A * B - C * C );

      /* store these in buffer if available */
      if ( cfBuffer )
	{
	  XLALEmptyComputeFBuffer ( *cfBuffer );
	  cfBuffer->multiSSB = multiSSB;
	  cfBuffer->multiAMcoef = multiAMcoef;
	  cfBuffer->skypos = psPoint->skypos;
	  cfBuffer->A = A;
	  cfBuffer->B = B;
	  cfBuffer->C = C;
	  cfBuffer->Dinv = Dinv;
	} /* if cfBuffer */

    } /* if no buffer or different skypos */

  /* ----- loop over detectors and compute all detector-specific quantities ----- */
  for ( X=0; X < numDetectors; X ++)
    {
      Fcomponents FcX = empty_Fcomponents;	/* for detector-specific FaX, FbX */
 		  
      if ( XLALComputeFaFb (&FcX, multiSFTs->data[X], psPoint->fkdot, multiSSB->data[X], multiAMcoef->data[X], params->Dterms) != 0)
	{
	  LALPrintError ("\nXALNewLALDemod() failed\n");
	  ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
	}
 		  
      retF.Fa.re += FcX.Fa.re;
      retF.Fa.im += FcX.Fa.im;
 		  
      retF.Fb.re += FcX.Fb.re;
      retF.Fb.im += FcX.Fb.im;
  		  
    } /* for  X < numDetectors */
 
  /* ----- compute final Fstatistic-value ----- */

  /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
   * therefore there is a factor of 2 difference with respect to the equations in JKS, which 
   * where based on the single-sided PSD.
   */ 
 		       
  retF.F = Dinv * (B * (SQ(retF.Fa.re) + SQ(retF.Fa.im) ) 
		   + A * ( SQ(retF.Fb.re) + SQ(retF.Fb.im) )
		   - 2.0 * C *( retF.Fa.re * retF.Fb.re + retF.Fa.im * retF.Fb.im )  
		   );

  (*Fstat) = retF;

  /* free memory if no buffer was available */
  if ( !cfBuffer )
    {
      XLALDestroyMultiSSBtimes ( multiSSB );
      XLALDestroyMultiAMCoeffs ( multiAMcoef );
    } /* if !cfBuffer */

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStat() */


/** Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 */
int
XLALComputeFaFb ( Fcomponents *FaFb,
		  const SFTVector *sfts,
		  const REAL8Vector *fkdot,
		  const SSBtimes *tSSB,
		  const AMCoeffs *amcoe,
		  UINT4 Dterms) 
{ 
  UINT4 alpha;                 	/* loop index over SFTs */
  UINT4 spdnOrder;		/* maximal spindown-orders */
  UINT4 numSFTs;		/* number of SFTs (M in the Notes) */
  COMPLEX16 Fa, Fb;
  REAL8 f;			/* !! MUST be REAL8, or precision breaks down !! */
  REAL8 Tsft; 			/* length of SFTs in seconds */
  INT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  INT4 freqIndex1;		/* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;		/* pointer to alpha-arrays over a and b */
  REAL8 *DeltaT_al, *Tdot_al;	/* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;		/* SFT alpha  */

  /* ----- check validity of input */
#ifndef LAL_NDEBUG
  if ( !FaFb ) {
    LALPrintError ("\nOutput-pointer is NULL !\n\n");
    XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    LALPrintError ("\nInput SFTs are NULL!\n\n");
    XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
  }
  
  if ( !fkdot || !tSSB || !tSSB->DeltaT || !tSSB->Tdot || !amcoe || !amcoe->a || !amcoe->b )
    {
      LALPrintError ("\nIllegal NULL in input !\n\n");
      XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
    }

  if ( fkdot->length >= NUM_FACT )
    {
      LALPrintError ("\nInverse factorials table only up to order s=%d, can't handle %d spin-order\n\n",
		     NUM_FACT, fkdot->length );
      XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
    }
#endif

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;

  freqIndex0 = (UINT4) ( sfts->data[0].f0 / sfts->data[0].deltaF + 0.5); /* lowest freqency-index */
  freqIndex1 = freqIndex0 + sfts->data[0].data->length;

  /* find highest non-zero spindown-entry */
  for ( spdnOrder = fkdot->length - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot->data[spdnOrder] )
      break;

  f = fkdot->data[0];

  Fa.re = 0.0f;
  Fa.im = 0.0f;
  Fb.re = 0.0f;
  Fb.im = 0.0f;

  a_al = amcoe->a->data;	/* point to beginning of alpha-arrays */
  b_al = amcoe->b->data;
  DeltaT_al = tSSB->DeltaT->data;
  Tdot_al = tSSB->Tdot->data;
  SFT_al = sfts->data;

  /* Loop over all SFTs  */
  for ( alpha = 0; alpha < numSFTs; alpha++ )
    {
      REAL4 a_alpha, b_alpha;

      INT4 kstar;		/* central frequency-bin k* = round(xhat_alpha) */
      INT4 k0, k1;

      COMPLEX8 *Xalpha = SFT_al->data->data; /* pointer to current SFT-data */
      COMPLEX8 *Xalpha_l; 	/* pointer to frequency-bin k in current SFT */
      REAL4 s_alpha, c_alpha;	/* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
      REAL4 realQ, imagQ;	/* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realXP, imagXP;	/* Re/Im of sum_k X_ak * P_ak */
      REAL4 realQXP, imagQXP;	/* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha, kappa_alpha;
      REAL8 remainder;

      /* ----- calculate kappa_alpha and lambda_alpha */
      {
	UINT4 s; 		/* loop-index over spindown-order */
	REAL8 phi_alpha, Dphi_alpha, DT_al;
	REAL8 Tas;	/* temporary variable to calculate (DeltaT_alpha)^s */
	
	/* init for s=0 */
	phi_alpha = 0.0;
	Dphi_alpha = 0.0;
	DT_al = (*DeltaT_al);
	Tas = 1.0;		/* DeltaT_alpha ^ 0 */

	for (s=0; s <= spdnOrder; s++)
	  {
	    REAL8 fsdot = fkdot->data[s];
	    Dphi_alpha += fsdot * Tas * inv_fact[s]; 	/* here: DT^s/s! */
	    Tas *= DT_al;				/* now: DT^(s+1) */
	    phi_alpha += fsdot * Tas * inv_fact[s+1];
	  } /* for s <= spdnOrder */

	/* Step 3: apply global factors to complete Dphi_alpha */
	Dphi_alpha *= Tsft * (*Tdot_al);		/* guaranteed > 0 ! */

	lambda_alpha = phi_alpha - 0.5 * Dphi_alpha;
	
	/* real- and imaginary part of e^{-i 2 pi lambda_alpha } */
	if ( sin_cos_2PI_LUT ( &imagQ, &realQ, - lambda_alpha ) ) {
	  XLAL_ERROR ( "XLALComputeFaFb", XLAL_EFUNC);
	}

	kstar = (INT4) (Dphi_alpha + 0.5);	/* k* = round(Dphi_alpha) for positive Dphi */
	remainder = (Dphi_alpha - kstar);	/* rem(Dphi_alpha) */
	kappa_alpha = remainder + Dterms;

	k0 = kstar - Dterms;	
	k1 = k0 + 2 * Dterms;
	if ( (k0 < freqIndex0) || (k1 > freqIndex1) ) 
	  {
	    LALPrintError ("Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n\n",
			   k0, k1, freqIndex0, freqIndex1 );
	    XLAL_ERROR("XLALComputeFaFb", XLAL_EDOM);
	  }

      } /* compute kappa_alpha, lambda_alpha */

      /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
       * the trig-functions need to be calculated only once!
       * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the 
       * closest to zero and will pose no numerical difficulties !
       */
      sin_cos_2PI_LUT ( &s_alpha, &c_alpha, remainder );
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
      if ( ( remainder > LD_SMALL4 ) || (remainder < -LD_SMALL4) )	
	{ 
	  /* improved hotloop algorithm by Fekete Akos: 
	   * take out repeated divisions into a single common denominator,
	   * plus use extra cleverness to compute the nominator efficiently...
	   */
	  COMPLEX8 Xal = *Xalpha_l;
	  REAL4 Sn = Xal.re;
	  REAL4 Tn = Xal.im;
	  REAL4 pn = kappa_alpha;
	  REAL4 qn = pn;
	  REAL4 U_alpha, V_alpha;
	  
	  /* recursion with 2*Dterms steps */
	  UINT4 l;
	  for ( l = 1; l <= 2*Dterms; l ++ )
	    {
	      Xalpha_l ++;
	      Xal = *Xalpha_l;
	      
	      pn = pn - 1.0f;			/* p_(n+1) */
	      Sn = pn * Sn + qn * Xal.re;	/* S_(n+1) */
	      Tn = pn * Tn + qn * Xal.im;	/* T_(n+1) */
	      qn *= pn;				/* q_(n+1) */
	    } /* for l <= 2*Dterms */
	  
	  U_alpha = Sn / qn;
	  V_alpha = Tn / qn;
	  
	  realXP = s_alpha * U_alpha - c_alpha * V_alpha;
	  imagXP = c_alpha * U_alpha + s_alpha * V_alpha;
	  
	} /* if |remainder| > LD_SMALL4 */
      else
	{ /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
	  realXP = TWOPI_FLOAT * Xalpha_l[Dterms].re;
	  imagXP = TWOPI_FLOAT * Xalpha_l[Dterms].im;
	} /* if |remainder| <= LD_SMALL4 */
      
      realQXP = realQ * realXP - imagQ * imagXP;
      imagQXP = realQ * imagXP + imagQ * realXP;
      
      /* we're done: ==> combine these into Fa and Fb */
      a_alpha = (*a_al);
      b_alpha = (*b_al);

      Fa.re += a_alpha * realQXP;
      Fa.im += a_alpha * imagQXP;
      
      Fb.re += b_alpha * realQXP;
      Fb.im += b_alpha * imagQXP;

      /* advance pointers over alpha */
      a_al ++;
      b_al ++;
      DeltaT_al ++;
      Tdot_al ++;
      SFT_al ++;

    } /* for alpha < numSFTs */
      
  /* return result */
  FaFb->Fa.re = Fa.re * OOTWOPI;
  FaFb->Fa.im = Fa.im * OOTWOPI;
  FaFb->Fb.re = Fb.re * OOTWOPI;
  FaFb->Fb.im = Fb.im * OOTWOPI;

  return XLAL_SUCCESS;

} /* XLALComputeFaFb() */

/** Compute the 'amplitude coefficients' \f$a(t), b(t)\f$ as defined in 
 * \ref JKS98 for a series of timestamps.
 * 
 * The input consists of the DetectorState-timeseries, which contains
 * the detector-info and the LMST's corresponding to the different times.
 * 
 * In order to allow re-using the output-structure AMCoeffs for subsequent
 * calls, we require the REAL4Vectors a and b to be allocated already and 
 * to have the same length as the DetectoStates-timeseries.
 *
 * \note This is an alternative implementation to LALComputeAM() with 
 * the aim to be both simpler and faster.
 * The difference being that we don't implicitly re-derive the final expression
 * here but simply try to implement the final expressions (12), (13) in \ref JKS98
 * in the most economical way possible.
 */ 
void
LALGetAMCoeffs(LALStatus *status,
	       AMCoeffs *coeffs,			/**< [out] amplitude-coeffs {a(t_i), b(t_i)} */
	       const DetectorStateSeries *DetectorStates,	/**< timeseries of detector states */
	       SkyPosition skypos			/**< {alpha,delta} of the source */
	       )
{
  REAL4 ah1, ah2, ah3, ah4, ah5;
  REAL4 a1, a2, a3, a4, a5;
  
  REAL4 bh1, bh2, bh3, bh4;
  REAL4 b1, b2, b3, b4;

  REAL4 delta, alpha;
  REAL4 sin1delta, cos1delta, sin2delta, cos2delta;

  REAL4 gamma, lambda;
  REAL4 norm;
  UINT4 i, numSteps;

  INITSTATUS (status, "LALGetAMCoeffs", COMPUTEFSTATC);

  /*---------- check input ---------- */
  ASSERT ( DetectorStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  numSteps = DetectorStates->length;

  /* require the coeffients-vectors to be allocated and consistent with timestamps */
  ASSERT ( coeffs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( coeffs->a && coeffs->b, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( (coeffs->a->length == numSteps) && (coeffs->b->length == numSteps), status,
	   COMPUTEFSTATC_EINPUT,  COMPUTEFSTATC_MSGEINPUT);

  /* require sky-pos to be in equatorial coordinates */
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, 
	   SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );

  /*---------- detector paramters: lambda, L, gamma */
  {
    /* FIXME: put into DetectorStateSeries */
    /* orientation of detector arms */
    REAL8 xAzi = DetectorStates->detector.frDetector.xArmAzimuthRadians;
    REAL8 yAzi = DetectorStates->detector.frDetector.yArmAzimuthRadians;

    /* get detector orientation gamma */
    gamma = LAL_PI_2 - 0.5 * (xAzi + yAzi);
    /* get detector position latitude (lambda) */
    lambda = DetectorStates->detector.frDetector.vertexLatitudeRadians;
  }

  /*---------- coefficient ahN, bhN dependent ONLY on detector-position  ---------- */
  /* FIXME: put these coefficients into DetectorStateSeries */
  {
    REAL4 sin2gamma, cos2gamma;
    REAL4 sin1lambda, cos1lambda;
    REAL4 sin2lambda, cos2lambda;

    sin_cos_LUT (&sin2gamma, &cos2gamma, 2.0f * gamma );
    sin_cos_LUT (&sin1lambda, &cos1lambda, lambda );

    sin2lambda = 2.0f * sin1lambda * cos1lambda;
    cos2lambda = cos1lambda * cos1lambda - sin1lambda * sin1lambda;

    /* coefficients for a(t) */
    ah1 = 0.0625f * sin2gamma * (3.0f - cos2lambda);	/* 1/16 = 0.0625 */
    ah2 = - 0.25f * cos2gamma * sin1lambda;
    ah3 =   0.25f * sin2gamma * sin2lambda;
    ah4 =  -0.5f  * cos2gamma * cos1lambda;
    ah5 =  0.75f  * sin2gamma * cos1lambda * cos1lambda;

    /* coefficients for b(t) */
    bh1 =           cos2gamma * sin1lambda;
    bh2 =   0.25f * sin2gamma * (3.0f - cos2lambda);
    bh3 =           cos2gamma * cos1lambda;
    bh4 =   0.5f  * sin2gamma * sin2lambda;
  }

  /*---------- coefficients aN, bN dependent ONLY on {ahN, bhN} and source-latitude delta */
  alpha = skypos.longitude;
  delta = skypos.latitude;

  sin_cos_LUT (&sin1delta, &cos1delta, delta );
  sin2delta = 2.0f * sin1delta * cos1delta;
  cos2delta = cos1delta * cos1delta - sin1delta * sin1delta;

  /* coefficients for a(t) */
  a1 = ah1 * ( 3.0f - cos2delta );
  a2 = ah2 * ( 3.0f - cos2delta );
  a3 = ah3 * sin2delta;
  a4 = ah4 * sin2delta;
  a5 = ah5 * cos1delta * cos1delta;

  /* coefficients for b(t) */
  b1 = bh1 * sin1delta;
  b2 = bh2 * sin1delta;
  b3 = bh3 * cos1delta;
  b4 = bh4 * cos1delta;


  /*---------- Compute the a(t_i) and b(t_i) ---------- */
  coeffs->A = 0;
  coeffs->B = 0;
  coeffs->C = 0;
  coeffs->D = 0;
  for ( i=0; i < numSteps; i++ )
    {
      REAL4 ah;
      REAL4 cos1ah, sin1ah, cos2ah, sin2ah;
      REAL4 ai, bi;

      ah = alpha - DetectorStates->data[i].LMST;

      sin_cos_LUT ( &sin1ah, &cos1ah, ah );
      sin2ah = 2.0f * sin1ah * cos1ah;
      cos2ah = cos1ah * cos1ah - sin1ah * sin1ah;

      ai = a1 * cos2ah + a2 * sin2ah + a3 * cos1ah + a4 * sin1ah + a5;
      bi = b1 * cos2ah + b2 * sin2ah + b3 * cos1ah + b4 * sin1ah;
      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;

      /* sum A, B, C on the fly */
      coeffs->A += ai * ai;
      coeffs->B += bi * bi;
      coeffs->C += ai * bi;

    } /* for i < numSteps */

  /* finish calculation of A,B,C, D */
  norm = 2.0f / numSteps;
  coeffs->A *= norm;
  coeffs->B *= norm;
  coeffs->C *= norm;

  coeffs->D = coeffs->A * coeffs->B - coeffs->C * coeffs->C;

  RETURN(status);

} /* LALGetAMCoeffs() */


/** For a given DetectorStateSeries, calculate the time-differences
 *  \f$\Delta T_\alpha\equiv T(t_\alpha) - T_0\f$, and their
 *  derivatives \f$Tdot_\alpha \equiv d T / d t (t_\alpha)\f$.
 * 
 *  \note The return-vectors \a DeltaT and \a Tdot must be allocated already
 *  and have the same length as the input time-series \a DetStates.
 *
 */
void
LALGetSSBtimes (LALStatus *status, 
		SSBtimes *tSSB,			/**< [out] DeltaT_alpha = T(t_alpha) - T_0; and Tdot(t_alpha) */
		const DetectorStateSeries *DetectorStates,/**< [in] detector-states at timestamps t_i */
		SkyPosition pos,		/**< source sky-location */
		LIGOTimeGPS refTime,		/**< SSB reference-time T_0 of pulsar-parameters */
		SSBprecision precision		/**< relativistic or Newtonian SSB transformation? */
		)
{
  UINT4 numSteps, i;
  REAL8 vn[3];		/* unit-vector pointing to source in Cart. coord. */
  REAL8 alpha, delta;	/* source position */
  REAL8 refTimeREAL8;

  INITSTATUS( status, "LALGetSSBtimes", COMPUTEFSTATC);
  ATTATCHSTATUSPTR (status);

  ASSERT (DetectorStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  numSteps = DetectorStates->length;		/* number of timestamps */

  ASSERT (tSSB, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (tSSB->DeltaT, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (tSSB->Tdot, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  ASSERT (tSSB->DeltaT->length == numSteps, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT (tSSB->Tdot->length == numSteps, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);

  ASSERT (precision < SSBPREC_LAST, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT ( pos.system == COORDINATESYSTEM_EQUATORIAL, status,
	   COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  

  /* convenience variables */
  alpha = pos.longitude;
  delta = pos.latitude;
  refTimeREAL8 = GPS2REAL8(refTime);

  /*----- now calculate the SSB transformation in the precision required */
  switch (precision)
    {
    case SSBPREC_NEWTONIAN:	/* use simple vr.vn to calculate time-delay */

      /*----- get the cartesian source unit-vector */
      vn[0] = cos(alpha) * cos(delta);
      vn[1] = sin(alpha) * cos(delta);
      vn[2] = sin(delta);

      for (i=0; i < numSteps; i++ )
	{
	  LIGOTimeGPS *ti = &(DetectorStates->data[i].tGPS);
	  /* DeltaT_alpha */
	  tSSB->DeltaT->data[i]  = GPS2REAL8 ( (*ti) );
	  tSSB->DeltaT->data[i] += SCALAR(vn, DetectorStates->data[i].rDetector);
	  tSSB->DeltaT->data[i] -= refTimeREAL8;

	  /* Tdot_alpha */
	  tSSB->Tdot->data[i] = 1.0 + SCALAR(vn, DetectorStates->data[i].vDetector);
	  
	} /* for i < numSteps */

      break;

    case SSBPREC_RELATIVISTIC:	/* use LALBarycenter() to get SSB-times and derivative */
      for (i=0; i < numSteps; i++ )
	{
	  BarycenterInput baryinput = empty_BarycenterInput;
	  EmissionTime emit;
	  DetectorState *state = &(DetectorStates->data[i]);

	  baryinput.tgps = state->tGPS;
	  baryinput.site = DetectorStates->detector;
	  /* ARGHHH!!! */
	  baryinput.site.location[0] /= LAL_C_SI;
	  baryinput.site.location[1] /= LAL_C_SI;
	  baryinput.site.location[2] /= LAL_C_SI;

	  baryinput.alpha = alpha;
	  baryinput.delta = delta;
	  baryinput.dInv = 0;

	  TRY ( LALBarycenter(status->statusPtr, &emit, &baryinput, &(state->earthState)), status);

	  tSSB->DeltaT->data[i] = GPS2REAL8 ( emit.te ) - refTimeREAL8;
	  tSSB->Tdot->data[i] = emit.tDot;

	} /* for i < numSteps */

      break;
    default:
      LALPrintError ("\n?? Something went wrong.. this should never be called!\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
      break;
    } /* switch precision */

  /* finally: store the reference-time used into the output-structure */
  tSSB->refTime = refTime;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetSSBtimes() */

/** Get the 'detector state' (ie position, velocity, etc) for the given
 * vector of timestamps, shifted by a common time-shift \a tOffset.
 *
 * This function just calls LALBarycenterEarth() and LALBarycenter() for the
 * given vector of timestamps (shifted by tOffset) and returns the positions, 
 * velocities and LMSTs of the detector, stored in a DetectorStateSeries. 
 * There is also an entry containing the EarthState at each timestamp, which 
 * can be used as input for subsequent calls to LALBarycenter().
 *
 * \a tOffset allows one to easily use the midpoints of SFT-timestamps, for example.
 *
 * \note the DetectorStateSeries is allocated here and should be free'ed with
 * LALDestroyDetectorStateSeries().
 *
 */
void
LALGetDetectorStates (LALStatus *status,
		      DetectorStateSeries **DetectorStates,	/**< [out] series of DetectorStates */
		      const LIGOTimeGPSVector *timestamps,	/**< array of GPS timestamps t_i */
		      const LALDetector *detector,		/**< detector info */
		      const EphemerisData *edat,		/**< ephemeris file data */	
		      REAL8 tOffset
		      )
{
  UINT4 i, j, numSteps;
  DetectorStateSeries *ret = NULL;

  INITSTATUS( status, "LALGetDetectorStates", COMPUTEFSTATC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( DetectorStates != NULL, status,  COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *DetectorStates == NULL, status,  COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);

  ASSERT ( timestamps, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( detector, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( edat, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  numSteps = timestamps->length;

  TRY ( LALCreateDetectorStateSeries (status->statusPtr, &ret, numSteps), status);

  /* enter detector-info into the head of the state-vector */
  ret->detector = (*detector);
  
  /* now fill all the vector-entries corresponding to different timestamps */
  for ( i=0; i < numSteps; i++ )
    {
      BarycenterInput baryinput;
      EmissionTime emit;
      DetectorState *state = &(ret->data[i]);
      EarthState *earth = &(state->earthState);
      LIGOTimeGPS tgps;

      /* shift timestamp by tOffset */
      TRY ( LALAddFloatToGPS(status->statusPtr, &tgps, &timestamps->data[i], tOffset ), status);

      /*----- first get earth-state */
      LALBarycenterEarth (status->statusPtr, earth, &tgps, edat );
      BEGINFAIL(status){
	TRY ( LALDestroyDetectorStateSeries(status->statusPtr, &ret), status);
      }ENDFAIL(status);
      /*----- then get detector-specific info */
      baryinput.tgps = tgps;			/* irrelevant here! */
      baryinput.site = (*detector);
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;
      baryinput.alpha = baryinput.delta = 0;	/* irrelevant */
      baryinput.dInv = 0;

      LALBarycenter (status->statusPtr, &emit, &baryinput, earth);
      BEGINFAIL(status) {
	TRY ( LALDestroyDetectorStateSeries(status->statusPtr, &ret), status);
      }ENDFAIL(status);

      /*----- extract the output-data from this */
      for (j=0; j < 3; j++)	/* copy detector's position and velocity */
	{
	  state->rDetector[j] = emit.rDetector[j];
	  state->vDetector[j] = emit.vDetector[j];
	} /* for j < 3 */

      /* local mean sidereal time = GMST + longitude */
      state->LMST = earth->gmstRad + detector->frDetector.vertexLongitudeRadians;
      state->LMST = fmod (state->LMST, LAL_TWOPI );	/* normalize */

      /* insert timestamp */
      state->tGPS = tgps;

    } /* for i < numSteps */

  /* return result */
  (*DetectorStates) = ret;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALGetDetectorStates() */


/* ===== Multi-IFO versions of some of the above functions ===== */

/** Get the detector-time series for the given MultiSFTVector. 
 * (see LALGetDetectorStates for more comments).
 */
void
LALGetMultiDetectorStates( LALStatus *status, 
			   MultiDetectorStateSeries **mdetStates, 	/**< [out] multi-IFO detector-states */
			   const MultiSFTVector *multiSFTs, 		/**< [in] multi-IFO SFTs */
			   const EphemerisData *edat )			/**< ephemeris files data */				   
{
  UINT4 X, numDetectors;
  MultiDetectorStateSeries *ret = NULL;

  INITSTATUS (status, "LALGetMultiDetectorStates", COMPUTEFSTATC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( mdetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( multiSFTs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *mdetStates == NULL, status, COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);

  numDetectors = multiSFTs->length;

  /* prepare return-structure */
  if ( ( ret = LALCalloc ( 1, sizeof( *ret ) )) == NULL ) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof( *(ret->data) ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  ret->length = numDetectors;

  /* loop over detectors */
  for ( X=0; X < numDetectors; X ++ )
    {
      LIGOTimeGPSVector *ts = NULL;
      LALDetector *det = NULL;

      SFTVector *this_sftvect = multiSFTs->data[X];
      REAL8 tOffs = 0.5 / this_sftvect->data[0].deltaF;	/* Tsft / 2 */

      /* timestamps from SFTVector  of detector X */
      LALGetSFTtimestamps ( status->statusPtr, &ts, this_sftvect );
      if ( status->statusPtr->statusCode ) 
	{
	  LALPrintError ( "\nCall to LALGetSFTtimestamps() has failed ... \n\n");
	  goto failed;
	}
      /* LALDetector struct for this detector */
      if ( (det = XLALGetSiteInfo ( this_sftvect->data[0].name )) == NULL ) 
	{
	  LALPrintError ("\nCall to XLALGetSiteInfo() has failed ... \n\n");
	  XLALDestroyTimestampVector ( ts );
	  goto failed;
	}
      /* fill in the detector-state series for this detector */
      LALGetDetectorStates (status->statusPtr, &(ret->data[X]), ts, det, edat, tOffs );
      if ( status->statusPtr->statusCode ) 
	{
	  LALPrintError ( "\nCall to LALGetDetectorStates() has failed ... \n\n");
	  XLALDestroyTimestampVector ( ts );
	  LALFree ( det );
	  goto failed;
	}

      /* free temporary mem */
      XLALDestroyTimestampVector ( ts );
      ts = NULL;
      LALFree ( det );

    } /* for X < numDetectors */

  goto success;

 failed:
  /* free complete MultiDetectorStateSeries built up so far */
  XLALDestroyMultiDetectorStateSeries ( ret );	/* NOTE: this function is "NULL-robust" */
  ABORT ( status, -1, "LALGetMultiDetectorStates failed" );
  
 success:

  (*mdetStates) = ret;
  
  DETATCHSTATUSPTR (status);
  RETURN ( status );

} /* LALGetMultiDetectorStates() */


/** Multi-IFO version of LALGetSSBtimes(). 
 * Get all SSB-timings for all input detector-series.
 *
 * NOTE: contrary to LALGetSSBtimes(), this functions *allocates* the output-vector,
 * use XLALDestroyMultiSSBtimes() to free this.
 */
void
LALGetMultiSSBtimes (LALStatus *status, 
		     MultiSSBtimes **multiSSB,		/**< [out] SSB-timings for all input detector-state series */
		     const MultiDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
		     SkyPosition skypos,		/**< source sky-position [in equatorial coords!] */
		     LIGOTimeGPS refTime,		/**< SSB reference-time T_0 for SSB-timing */
		     SSBprecision precision		/**< use relativistic or Newtonian SSB timing?  */
		     )
{
  UINT4 X, numDetectors;
  MultiSSBtimes *ret = NULL;

  INITSTATUS( status, "LALGetMultiSSBtimes", COMPUTEFSTATC);
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT (multiDetStates, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (multiDetStates->length, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (multiSSB, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *multiSSB == NULL, status,COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );

  numDetectors = multiDetStates->length;

  if ( ( ret = LALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);    
  }
  ret->length = numDetectors;
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      SSBtimes *SSBtimesX = NULL;
      UINT4 numStepsX = multiDetStates->data[X]->length;

      ret->data[X] = LALCalloc ( 1, sizeof ( *(ret->data[X]) ) );
      SSBtimesX = ret->data[X];
      SSBtimesX->DeltaT = XLALCreateREAL8Vector ( numStepsX );
      if ( (SSBtimesX->Tdot = XLALCreateREAL8Vector ( numStepsX )) == NULL ) {
	LALPrintError ("\nOut of memory!\n\n");
	goto failed;
      }

      LALGetSSBtimes (status->statusPtr, SSBtimesX, multiDetStates->data[X], skypos, refTime, precision );
      if ( status->statusPtr->statusCode ) 
	{
	  LALPrintError ( "\nCall to LALGetSSBtimes() has failed ... \n\n");
	  goto failed;
	}
 
    } /* for X < numDet */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiSSBtimes ( ret );
  ABORT ( status, -1, "LALGetMultiSSBtimes failed" );

 success:
  (*multiSSB) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetMultiSSBtimes() */

/** Multi-IFO version of LALGetAMCoeffs(). 
 * Get all antenna-pattern coefficients for all input detector-series.
 *
 * NOTE: contrary to LALGetAMCoeffs(), this functions *allocates* the output-vector,
 * use XLALDestroyMultiAMCoeffs() to free this.
 */
void
LALGetMultiAMCoeffs (LALStatus *status, 
		     MultiAMCoeffs **multiAMcoef,	/**< [out] AM-coefficients for all input detector-state series */
		     const MultiDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
		     SkyPosition skypos			/**< source sky-position [in equatorial coords!] */
		     )
{
  UINT4 X, numDetectors;
  MultiAMCoeffs *ret = NULL;

  INITSTATUS( status, "LALGetMultiAMCoeffs", COMPUTEFSTATC);
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT (multiDetStates, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (multiDetStates->length, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (multiAMcoef, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *multiAMcoef == NULL, status,COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );

  numDetectors = multiDetStates->length;

  if ( ( ret = LALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);    
  }
  ret->length = numDetectors;
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      AMCoeffs *amcoeX = NULL;
      UINT4 numStepsX = multiDetStates->data[X]->length;

      ret->data[X] = LALCalloc ( 1, sizeof ( *(ret->data[X]) ) );
      amcoeX = ret->data[X];
      amcoeX->a = XLALCreateREAL4Vector ( numStepsX );
      if ( (amcoeX->b = XLALCreateREAL4Vector ( numStepsX )) == NULL ) {
	LALPrintError ("\nOut of memory!\n\n");
	goto failed;
      }

      LALGetAMCoeffs (status->statusPtr, amcoeX, multiDetStates->data[X], skypos );
      if ( status->statusPtr->statusCode ) 
	{
	  LALPrintError ( "\nCall to LALGetAMCoeffs() has failed ... \n\n");
	  goto failed;
	}
 
    } /* for X < numDetectors */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiAMCoeffs ( ret );
  ABORT ( status, -1, "LALGetMultiSSBtimes failed" );

 success:
  (*multiAMcoef) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetMultiAMCoeffs() */



/* ===== Object creation/destruction functions ===== */

/** Create a DetectorStateSeries */
void
LALCreateDetectorStateSeries (LALStatus *status, 
			      DetectorStateSeries **vect,	/**< output vector */
			      UINT4 length )			/**< number of entries */
{
  DetectorStateSeries *ret = NULL;

  INITSTATUS (status, "LALCreateDetectorStateSeries", COMPUTEFSTATC );

  ASSERT ( vect, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *vect == NULL, status, COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);

  if ( (ret = LALCalloc(1, sizeof(DetectorStateSeries) )) == NULL ) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  if ( (ret->data = LALCalloc (length, sizeof(DetectorState) )) == NULL ) {
    LALFree (ret);
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  ret->length = length;

  /* return result */
  (*vect) = ret;

  RETURN (status);

} /* LALCreateDetectorStateSeries() */

/* Get rid of a DetectorStateSeries */
void
XLALDestroyDetectorStateSeries ( DetectorStateSeries *detStates )
{
  if ( !detStates )
    return;

  if ( detStates->data ) LALFree ( detStates->data );
  LALFree ( detStates );

  return;

} /* XLALDestroyDetectorStateSeries() */

/** Destroy a DetectorStateSeries (and set it to NULL) */
void
LALDestroyDetectorStateSeries (LALStatus *status, 
			       DetectorStateSeries **detStates ) /**< pointer to vector to be destroyed */
{
  INITSTATUS (status, "LALDestroyDetectorStateSeries", COMPUTEFSTATC );

  ASSERT ( detStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  XLALDestroyDetectorStateSeries ( (*detStates) );

  (*detStates) = NULL;

  RETURN (status);
} /* LALDestroyDetectorStateSeries() */

/** Helper function to get rid of a multi-IFO DetectorStateSeries 
 * Note, this is "NULL-robust" in the sense that it will not crash 
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs.
 */
void
XLALDestroyMultiDetectorStateSeries ( MultiDetectorStateSeries *mdetStates )
{
  UINT4 X, numDet;

  if ( !mdetStates )
    return;

  numDet = mdetStates->length;
  if ( mdetStates->data )
    {
      for ( X=0; X < numDet ; X ++ )
	XLALDestroyDetectorStateSeries ( mdetStates->data[X] );

      LALFree ( mdetStates->data );
    }

  LALFree ( mdetStates );

  return;

} /* XLALDestroyMultiDetectorStateSeries() */

/** Destroy a MultiSSBtimes structure. 
 * Note, this is "NULL-robust" in the sense that it will not crash 
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs 
 */
void
XLALDestroyMultiSSBtimes ( MultiSSBtimes *multiSSB )
{
  UINT4 X;
  SSBtimes *tmp;

  if ( ! multiSSB )
    return;

  if ( multiSSB->data )
    {
      for ( X=0; X < multiSSB->length; X ++ ) 
	{
	  if ( (tmp = multiSSB->data[X]) != NULL )
	    {
	      if ( tmp->DeltaT )
		XLALDestroyREAL8Vector ( tmp->DeltaT );
	      if ( tmp->Tdot )
		XLALDestroyREAL8Vector ( tmp->Tdot );
	      LALFree ( tmp );
	    } /* if multiSSB->data[X] */
	} /* for X < numDetectors */
      LALFree ( multiSSB->data );
    }
  LALFree ( multiSSB );

  return;

} /* XLALDestroyMultiSSBtimes() */

/** Destroy a MultiAMCoeffs structure. 
 * Note, this is "NULL-robust" in the sense that it will not crash 
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs 
 */
void
XLALDestroyMultiAMCoeffs ( MultiAMCoeffs *multiAMcoef )
{
  UINT4 X;
  AMCoeffs *tmp;

  if ( ! multiAMcoef )
    return;

  if ( multiAMcoef->data )
    {
      for ( X=0; X < multiAMcoef->length; X ++ ) 
	{
	  if ( (tmp = multiAMcoef->data[X]) != NULL )
	    {
	      if ( tmp->a )
		XLALDestroyREAL4Vector ( tmp->a );
	      if ( tmp->b )
		XLALDestroyREAL4Vector ( tmp->b );
	      LALFree ( tmp );
	    } /* if multiAMcoef->data[X] */
	} /* for X < numDetectors */
      LALFree ( multiAMcoef->data );
    }
  LALFree ( multiAMcoef );

  return;

} /* XLALDestroyMultiAMCoeffs() */


/** Destruction of a ComputeFBuffer *contents*, 
 * i.e. the multiSSB and multiAMcoeff, while the 
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...) */
void
XLALEmptyComputeFBuffer ( ComputeFBuffer cfb )
{
  XLALDestroyMultiSSBtimes ( cfb.multiSSB );
  cfb.multiSSB = NULL;
  XLALDestroyMultiAMCoeffs ( cfb.multiAMcoef );
  cfb.multiAMcoef = NULL;

  return;
} /* XLALDestroyComputeFBuffer() */

/* ===== General internal helper functions ===== */

/** Calculate sin(x) and cos(x) to roughly 1e-7 precision using 
 * a lookup-table and Tayler-expansion.
 *
 * NOTE: this function will fail for arguments larger than
 * |x| > INT4_MAX = 2147483647 ~ 2e9 !!!
 *
 * return = 0: OK, nonzero=ERROR
 */
int
sin_cos_LUT (REAL4 *sinx, REAL4 *cosx, REAL8 x)
{
  return sin_cos_2PI_LUT ( sinx, cosx, x * OOTWOPI );
} /* sin_cos_LUT() */

#define LUT_RES         64      /* resolution of lookup-table */
#define LUT_RES_F	(1.0 * LUT_RES)
#define OO_LUT_RES	(1.0 / LUT_RES)

#define X_TO_IND	(1.0 * LUT_RES * OOTWOPI )
#define IND_TO_X	(LAL_TWOPI * OO_LUT_RES)
int
sin_cos_2PI_LUT (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x)
{
  REAL8 xt;
  INT4 i0;
  REAL8 d, d2;
  REAL8 ts, tc;

  static BOOLEAN firstCall = TRUE;
  static REAL4 sinVal[LUT_RES+1], cosVal[LUT_RES+1];

  /* the first time we get called, we set up the lookup-table */
  if ( firstCall )
    {
      UINT4 k;
      for (k=0; k <= LUT_RES; k++)
        {
          sinVal[k] = sin( LAL_TWOPI * k * OO_LUT_RES );
          cosVal[k] = cos( LAL_TWOPI * k * OO_LUT_RES );
        }
      firstCall = FALSE;
    }

  xt = x - (INT4)x;		/* xt in (-1, 1) */
  if ( xt < 0.0 )
    xt += 1.0;			/* xt in [0, 1 ) */
#ifndef LAL_NDEBUG
  if ( xt < 0.0 || xt > 1.0 )
    {
      LALPrintError("\nFailed numerica in sin_cos_2PI_LUT(): xt = %f not in [0,1)\n\n", xt );
      return XLAL_FAILURE;
    }
#endif

  i0 = (INT4)( xt * LUT_RES_F + 0.5 );	/* i0 in [0, LUT_RES ] */
  d = d2 = LAL_TWOPI * (xt - OO_LUT_RES * i0);
  d2 *= 0.5 * d;

  ts = sinVal[i0];
  tc = cosVal[i0];
   
  /* use Taylor-expansions for sin/cos around LUT-points */
  (*sin2pix) = ts + d * tc - d2 * ts;
  (*cos2pix) = tc - d * ts - d2 * tc;

  return XLAL_SUCCESS;
} /* sin_cos_2PI_LUT() */

