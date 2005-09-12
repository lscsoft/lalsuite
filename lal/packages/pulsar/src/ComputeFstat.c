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

#include "ComputeFstat.h"

NRCSID( COMPUTEFSTATC, "$Id$");

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Macros ----- */
/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])


/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */
static const BarycenterInput empty_BarycenterInput;

/*---------- Global variables ----------*/


/*---------- internal prototypes ----------*/
int sin_cos_LUT (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x); /* LUT-calculation of sin/cos */



/*==================== FUNCTION DEFINITIONS ====================*/


#define LD_SMALL4       (1.0e-6)		/**< "small" number for REAL4*/
#define OOTWOPI         (1.0 / LAL_TWOPI)	/**< 1/2pi */

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */ 

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
  UINT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  UINT4 freqIndex1;		/* index of last frequency-bin in SFTs */
  COMPLEX16 Fa, Fb;
  REAL8 f;		/* !! MUST be REAL8, or precision breaks down !! */
  REAL8 Tsft; 			/* length of SFTs in seconds */

  /* ----- check validity of input */
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

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;

  freqIndex0 = (UINT4) ( sfts->data[0].f0 / sfts->data[0].deltaF + 0.5); /* lowest freqency-index */
  freqIndex1 = freqIndex0 + sfts->data[0].data->length;

  spdnOrder = fkdot->length - 1;

  f = fkdot->data[0];

  Fa.re = 0.0;
  Fa.im = 0.0;
  Fb.re = 0.0;
  Fb.im = 0.0;


  /* Loop over all SFTs  */
  for ( alpha = 0; alpha < numSFTs; alpha++ )
    {
      REAL4 a = amcoe->a->data[alpha];
      REAL4 b = amcoe->b->data[alpha];

      REAL8 xhat_alpha, y_alpha;	/* xhat(alpha), y(alpha): need to be REAL8 !! */
      REAL4 x0;
      UINT4 k;			/* loop index over frequency-bins */
      UINT4 kstar;		/* central frequency-bin k* = round(xhat_alpha) */

      COMPLEX8 *Xalpha = sfts->data[alpha].data->data; /* pointer to current SFT-data */
      COMPLEX8 *Xalpha_k; 	/* pointer to frequency-bin k in current SFT */
      REAL4 sinx, cosxm1;	/* sin(x_alpha) and (cos(x_alpha)-1) */
      REAL4 realXP, imagXP;	/* the sum_k X_alpha_k P_alpha_k */
      REAL4 realQ, imagQ;	/* Re and Im of Q = e^{-i y} */
      REAL4 realQXP, imagQXP;	/* Re/Im of Q_alpha XP_alpha */
      UINT4 k0, k1;
      REAL4 x0_remainder;

      /* ----- calculate x(alpha,0) and y(alpha) */
      {
	UINT4 s; 		/* loop-index over spindown-order */
	REAL8 Tas; 		/* temporary variable to calculate (DeltaT_alpha)^2 */
	UINT4 sfact = 1;	/* store for s! */
	REAL8 DeltaTalpha = tSSB->DeltaT->data[alpha];
	Tas = 1.0; 	/* DeltaT_alpha = T^1 */

	/* Step 1: s = 0 */
	xhat_alpha = f * Tas;	/* f^{0) T^0 / 0! */
	Tas *= DeltaTalpha;
	y_alpha = f * Tas;	/* f^{0} T^1 / 1! */

	/* Step 2: sum s >= 1 */
	for (s=1; s <= spdnOrder; s++)
	  {
	    REAL8 fsdot = fkdot->data[s];
	    xhat_alpha += fsdot * Tas / sfact; 	/* Tas = T^s here, sfact=s! */
	    Tas *= DeltaTalpha; 		/* T^(s+1) */
	    sfact *= (s+1);			/* (s+1)! */	  
	    y_alpha += fsdot * Tas / sfact; 
	  } /* for s <= spdnOrder */

	/* Step 3: apply global factors and complete y_alpha */
	xhat_alpha *= Tsft * tSSB->Tdot->data[alpha];	/* guaranteed > 0 ! */
	y_alpha -= 0.5 * xhat_alpha;
	
	/* real- and imaginary part of e^{-i 2 pi y } */
	if ( sin_cos_LUT ( &imagQ, &realQ, y_alpha ) ) {
	  XLAL_ERROR ( "XLALComputeFaFb", XLAL_EFUNC);
	}
	imagQ = -imagQ;
      }
      /* ---------------------------------------- */

      /* xhat_alpha determines the 'central' frequency-bin k* in the sum */
      kstar = (UINT4) (xhat_alpha + 0.5);	/* k* = round(xhat_alpha) */

      /* Trick: sin[ 2pi (xhat - k) ] = sin [ 2pi xhat ], therefore
       * the trig-functions need to be calculated only once!
       * We choose the value sin[ 2pi(xhat - kstar) ] because it is the 
       * smallest and will pose no numerical difficulties !
       */

      /*-------------------- calculate sin(x), cos(x) */
      sin_cos_LUT ( &sinx, &cosxm1, xhat_alpha );
      cosxm1 -= 1.0f; 
      /*-------------------- */

      realXP = 0;
      imagXP = 0;

      k0 = kstar - Dterms;
      k1 = k0 + 2 * Dterms;
      if ( (k0 < freqIndex0) || (k1 > freqIndex1) ) {
	LALPrintError ("Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n\n",
		       k0, k1, freqIndex0, freqIndex1 );
	XLAL_ERROR("XLALComputeFaFb", XLAL_EDOM);
      }

      /* ---------- calculate the (truncated to Dterms) sum over k ---------- */

      /* ---------- ATTENTION: this the "hot-loop", which will be 
       * executed many millions of times, so anything in here 
       * has a HUGE impact on the whole performance of the code.
       * 
       * DON'T touch *anything* in here unless you really know 
       * what you're doing !!
       *------------------------------------------------------------
       */

      Xalpha_k = Xalpha + k0 - freqIndex0;  /* first frequency-bin in sum */
      x0 = (REAL4)(xhat_alpha - (REAL8)k0);	/* first xhat-value in the loop */

      /* we branch now (instead of inside the central loop)
       * depending on wether x0 can ever become SMALL in the loop or not, 
       * because it requires special treatment in the Dirichlet kernel
       * [We use that fact here that xhat_alpha > 0 ! (see above)] 
       */
      x0_remainder = (REAL4)( xhat_alpha - kstar );
      if ( x0_remainder < LD_SMALL4 ) /* too close to an integer? */
	{
	  /* count down 2*Dterms values */
	  for ( k = 2 * Dterms; k != 0;  k -- )
	    {
	      REAL4 realP, imagP;	/* real and imaginary parts of Dirichlet-kernel P_alpha_k */
	      COMPLEX8 Xa = *Xalpha_k;
	      REAL4 xinv;
	      
	      /* calculate Dirichlet-kernel: P_alpha_k */
	      if( fabs(x0) <  LD_SMALL4 ) /* If x0 is small: correct x->0 limit : P_apha_k = 1 */
		{
		  realXP += Xa.re;
		  imagXP += Xa.im;
		} /* x0 too near zero */      
	      else 	
		{ /* safe to invert x0 */
		  xinv = OOTWOPI_FLOAT / x0;
		  realP = sinx * xinv;
		  imagP = cosxm1 * xinv;
		  
		  /* calculate P_alpha_k * X_alpha_k */
		  realXP += realP * Xa.re - imagP * Xa.im;
		  imagXP += imagP * Xa.re + realP * Xa.im;
		} /* x0 not near zero */
	      
	      Xalpha_k ++;	/* point to next frequency-bin */
	      x0 -- ;	/* x0-value for next iteration */
	      
	    } /* for k=kstar-Dterms to kstar+Dterms */
	  
	} /* if x could become close to 0 */
      else
	{ /* normal loop: no danger of x0 becoming zero.. */
	  
	  /* count down 2*Dterms values */
	  for ( k = 2 * Dterms; k != 0;  k -- )
	    {
	      REAL4 realP, imagP;	/* real and imaginary parts of Dirichlet-kernel P_alpha_k */
	      COMPLEX8 Xa = *Xalpha_k;
	      REAL4 xinv = OOTWOPI_FLOAT / x0;
	      
	      /* calculate P_alpha_k */
	      realP = sinx * xinv;
	      imagP = cosxm1 * xinv;
	      
	      /* calculate P_alpha_k * X_alpha_k */
	      realXP += realP * Xa.re - imagP * Xa.im;
	      imagXP += imagP * Xa.re + realP * Xa.im;
	      
	      Xalpha_k ++;	/* point to next frequency-bin */
	      x0 -- ;	/* x0-value for next iteration */
	      
	    } /* for k=kstar-Dterms to kstar+Dterms */
	  
	} /* normal loop: no danger of x0 becoming zero */

      realQXP = realQ * realXP - imagQ * imagXP;
      imagQXP = realQ * imagXP + imagQ * realXP;
      
      /* we're done: ==> combine these into Fa and Fb */
      Fa.re += a * realQXP;
      Fa.im += a * imagQXP;
      
      Fb.re += b * realQXP;
      Fb.im += b * imagQXP;
      
#ifdef HAVE_ISFINITE
      if ( !isfinite(Fa.re) || !isfinite(Fa.im) || !isfinite(Fb.re) || !isfinite(Fb.im) )
	{
	  LALPrintError("\nnon-normal number encountered in SFT-loop alpha=%d!\n", alpha );
	  LALPrintError("Fa = %f + i %f, Fb = %f + i %f\n\n", Fa.re, Fa.im, Fb.re, Fb.im );
	  XLAL_ERROR("XLALComputeFaFb", XLAL_ERANGE);
	}
#endif

    } /* for alpha < numSFTs */
      
  /* return result */
  FaFb->Fa = Fa;
  FaFb->Fb = Fb;

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
    REAL4 sin2gamma = sinf ( 2.0f * gamma );
    REAL4 cos2gamma = cosf ( 2.0f * gamma );
    REAL4 sin1lambda = sinf ( lambda );
    REAL4 cos1lambda = cosf ( lambda );
    REAL4 sin2lambda = 2.0f * sin1lambda * cos1lambda;
    REAL4 cos2lambda = cos1lambda * cos1lambda - sin1lambda * sin1lambda;

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
  sin1delta = sinf (delta);
  cos1delta = cosf (delta);
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

      sin1ah = sinf ( ah );
      cos1ah = cosf ( ah );
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


/** For a given vector of GPS-times, calculate the time-differences
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
  
  /*----- get the cartesian source unit-vector */
  alpha = pos.longitude;
  delta = pos.latitude;
  vn[0] = cos(alpha) * cos(delta);
  vn[1] = sin(alpha) * cos(delta);
  vn[2] = sin(delta);

  /*----- now calculate the SSB transformation in the precision required */
  switch (precision)
    {
    case SSBPREC_NEWTONIAN:	/* use simple vr.vn to calculate time-delay */
      /*----- first figure out reference-time in SSB */

      for (i=0; i < numSteps; i++ )
	{
	  LIGOTimeGPS *ti = &(DetectorStates->data[i].tGPS);
	  /* DeltaT_alpha */
	  tSSB->DeltaT->data[i]  = GPS2REAL8 ( (*ti) );
	  tSSB->DeltaT->data[i] += SCALAR(vn, DetectorStates->data[i].rDetector);
	  tSSB->DeltaT->data[i] -= GPS2REAL8 ( refTime );

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

	  tSSB->DeltaT->data[i] = GPS2REAL8 ( emit.te ) - GPS2REAL8 ( refTime );
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
		      const EphemerisData *edat,		/**< ephemeris-files */	
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


/** Destroy a DetectorStateSeries (and set it to NULL) */
void
LALDestroyDetectorStateSeries (LALStatus *status, 
			       DetectorStateSeries **vect ) /**< pointer to vector to be destroyed */
{
  INITSTATUS (status, "LALDestroyDetectorStateSeries", COMPUTEFSTATC );

  ASSERT ( vect, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  if ( *vect != NULL ) 
    {
      LALFree ( (*vect)->data );
      LALFree ( *vect );
      *vect = NULL;
    }

  RETURN (status);
} /* LALDestroyDetectorStateSeries() */


/** Calculate sin(2 pi x) and cos(2 pi x) to roughly 1e-7 error using 
 * a lookup-table and Tayler-expansion.
 * This is meant to be fast, so we don't even check the input-pointers...
 *
 * However, for numerical sanity&safty, we *DO* check if the resulting
 * index is within bounds, which can fail in case the argument x is too large..
 *
 * return = 0: OK, nonzero=ERROR
 */
int
sin_cos_LUT (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x)
{
#define LUT_RES         64      /* resolution of lookup-table */
  UINT4 ind; 
  REAL8 rem;
  static BOOLEAN firstCall = TRUE;
  static REAL4 sinVal[LUT_RES+1], cosVal[LUT_RES+1];


  if ( firstCall )
    {
      UINT4 k;
      for (k=0; k <= LUT_RES; k++)
        {
          sinVal[k] = sin( (LAL_TWOPI*k)/LUT_RES );
          cosVal[k] = cos( (LAL_TWOPI*k)/LUT_RES );
        }
      firstCall = FALSE;
    }

  rem = x - (INT8)x;	/* rem in (-1, 1) */
  if ( rem < 0 )
    rem += 1.0;		/* rem in [0, 1) */

  /* security check if we didn't overstretch the numerics here (can happen for x too large) */
  if ( (rem < 0) || (rem > 1) )
    {
      LALPrintError ("\nLUT-index out of bounds. Input argument was probably too large!\n\n");
      XLAL_ERROR ( "sin_cos_LUT", XLAL_EDOM);
    }
  
			   
  ind = (UINT4)( rem * LUT_RES + 0.5 );   /* closest LUT-entry */
  {
    REAL8 d = LAL_TWOPI *(rem - (REAL8)ind/(REAL8)LUT_RES);
    REAL8 d2 = 0.5 * d * d;
    REAL8 ts = sinVal[ind];
    REAL8 tc = cosVal[ind];
                
    (*sin2pix) = ts + d * tc - d2 * ts;
    (*cos2pix) = tc - d * ts - d2 * tc;
  }
  
  return XLAL_SUCCESS;
} /* sin_cos_LUT() */

