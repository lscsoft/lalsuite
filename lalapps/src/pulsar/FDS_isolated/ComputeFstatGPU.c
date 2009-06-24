/*
 * Copyright (C) 2009 Reinhard Prix
 * Copyright (C) 2007 Chris Messenger
 * Copyright (C) 2006 John T. Whelan, Badri Krishnan
 * Copyright (C) 2005, 2006, 2007 Reinhard Prix
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

/** \author R. Prix, J. T. Whelan
 * \ingroup pulsarCoherent
 * \file
 * \brief
 * Functions to calculate the so-called F-statistic for a given point in parameter-space,
 * following the equations in \ref JKS98.
 *
 * This code is partly a descendant of an earlier implementation found in
 * LALDemod.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens, Bruce Allen
 * ComputSky.[ch] by Jolien Creighton, Reinhard Prix, Steve Berukoff
 * LALComputeAM.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens
 *
 * NOTE: this file contains specialized versions of the central CFSv2 routines, that are aimed at GPU optimization.
 * At the moment this only means they are internally using only single precision, but still aggree to within 1%
 * for Tobs ~ 1day and fmax ~ 1kHz.
 *
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/ComputeFstat.h>

#include "ComputeFstatGPU.h"

NRCSID( COMPUTEFSTATC, "$Id$");

/*---------- local DEFINES ----------*/
#define LD_SMALL4       ((REAL4)2.0e-4)		/**< "small" number for REAL4*/

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */

/*----- Macros ----- */
#define SQ(x) ( (x) * (x) )
#define REM(x) ( (x) - (INT4)(x) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
static const REAL4 inv_fact[PULSAR_MAX_SPINS] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* empty initializers  */
static const LALStatus empty_LALStatus;
static const AMCoeffs empty_AMCoeffs;

const SSBtimes empty_SSBtimes;
const MultiSSBtimes empty_MultiSSBtimes;
const AntennaPatternMatrix empty_AntennaPatternMatrix;
const MultiAMCoeffs empty_MultiAMCoeffs;
const Fcomponents empty_Fcomponents;
const ComputeFBuffer empty_ComputeFBuffer;
const PulsarSpins_REAL4 empty_PulsarSpins_REAL4;
const ComputeFBuffer_REAL4 empty_ComputeFBuffer_REAL4;
const Fcomponents_REAL4 empty_Fcomponents_REAL4;

/*---------- internal prototypes ----------*/
int finite(double x);

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
 * NOTE3: This is a special "single-precision" version aimed at GPU optimization
 *
 *
 */
void
ComputeFStat_REAL4 ( LALStatus *status,
                     Fcomponents *Fstat,                 		/**< [out] Fstatistic + Fa, Fb */
                     const PulsarDopplerParams *doppler, 		/**< parameter-space point to compute F for */
                     const MultiSFTVector *multiSFTs,    		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
                     const MultiNoiseWeights *multiWeights,		/**< noise-weights of all SFTs */
                     const MultiDetectorStateSeries *multiDetStates,	/**< 'trajectories' of the different IFOs */
                     UINT4 Dterms,       				/**< number of terms to use in Dirichlet kernel interpolation */
                     ComputeFBuffer_REAL4 *cfBuffer            		/**< CF-internal buffering structure */
                     )
{
  static const CHAR *fn = "ComputeFStat_REAL4()";

  Fcomponents_REAL4 retF = empty_Fcomponents_REAL4;
  UINT4 X, numDetectors;
  MultiSSBtimes_REAL4 *multiSSB4 = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  REAL4 Ad, Bd, Cd, Dd_inv;
  SkyPosition skypos;

  INITSTATUS( status, fn, COMPUTEFSTATC );
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT ( Fstat, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiSFTs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  numDetectors = multiSFTs->length;
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  if ( multiWeights ) {
    ASSERT ( multiWeights->length == numDetectors , status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  }

  skypos.system =   COORDINATESYSTEM_EQUATORIAL;
  skypos.longitude = doppler->Alpha;
  skypos.latitude  = doppler->Delta;

  /* check if that skyposition SSB+AMcoef were already buffered */
  if ( cfBuffer
       && ( cfBuffer->multiDetStates == multiDetStates )
       && ( cfBuffer->Alpha == doppler->Alpha )
       && ( cfBuffer->Delta == doppler->Delta )
       && cfBuffer->multiSSB
       && cfBuffer->multiAMcoef
       )
    { /* yes ==> reuse */
      multiSSB4 = cfBuffer->multiSSB;
      multiAMcoef = cfBuffer->multiAMcoef;

    } /* if have buffered stuff to reuse */
  else
    {
      /* compute new SSB timings */
      if ( (multiSSB4 = XLALGetMultiSSBtimes_REAL4 ( multiDetStates, doppler->Alpha, doppler->Delta, doppler->refTime)) == NULL ) {
        XLALPrintError ( "%s: XLALGetMultiSSBtimes_REAL4() failed. xlalErrno = %d.\n", fn, xlalErrno );
	ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
      }

      /* compute new AM-coefficients */
      LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, skypos );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes_REAL4 ( multiSSB4 );
      } ENDFAIL (status);

      /* noise-weight Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	LALPrintError("\nXLALWeighMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
      }

      /* store these in buffer if available */
      if ( cfBuffer )
	{
	  cfBuffer->Alpha = doppler->Alpha;
	  cfBuffer->Delta = doppler->Delta;
	  cfBuffer->multiDetStates = multiDetStates ;

	  XLALDestroyMultiAMCoeffs ( cfBuffer->multiAMcoef );
	  XLALDestroyMultiSSBtimes_REAL4 ( cfBuffer->multiSSB );

	  cfBuffer->multiSSB = multiSSB4;
	  cfBuffer->multiAMcoef = multiAMcoef;

	} /* buffer new SSB times */

    } /* could NOT reuse previously buffered quantites */

  /* ---------- prepare data in pure REAL4 quantities to be passed into
   * the central single-precision routine XLALComputeFaFb_REAL4()
   */
  UINT4 s, numSpins = PULSAR_MAX_SPINS;
  PulsarSpins_REAL4 spins = empty_PulsarSpins_REAL4;
  spins.FreqMain = (INT4)doppler->fkdot[0];
  /* fkdot[0] now only carries the *remainder* of Freq wrt FreqMain */
  spins.fkdot[0] = (REAL4)( doppler->fkdot[0] - (REAL8)spins.FreqMain );
  /* the remaining spins are simply down-cast to REAL4 */
  for ( s=1; s < numSpins; s ++ )
    spins.fkdot[s] = (REAL4) doppler->fkdot[s];

  /* ----- loop over detectors and compute all detector-specific quantities ----- */
  for ( X=0; X < numDetectors; X ++)
    {
      Fcomponents_REAL4 FcX = empty_Fcomponents_REAL4;	/* for detector-specific FaX, FbX */

      if ( XLALComputeFaFb_REAL4 (&FcX, multiSFTs->data[X], spins, multiSSB4->data[X], multiAMcoef->data[X], Dterms) != 0)
        {
          XLALPrintError ("%s: XALComputeFaFb_REAL4() failed\n", fn );
          ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
        }

#ifndef LAL_NDEBUG
      if ( !finite(FcX.Fa.re) || !finite(FcX.Fa.im) || !finite(FcX.Fb.re) || !finite(FcX.Fb.im) ) {
	XLALPrintError("%s: XLALComputeFaFb_REAL4() returned non-finite: Fa=(%f,%f), Fb=(%f,%f)\n",
                       fn, FcX.Fa.re, FcX.Fa.im, FcX.Fb.re, FcX.Fb.im );
	ABORT (status,  COMPUTEFSTATC_EIEEE,  COMPUTEFSTATC_MSGEIEEE);
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
  Ad = multiAMcoef->Mmunu.Ad;
  Bd = multiAMcoef->Mmunu.Bd;
  Cd = multiAMcoef->Mmunu.Cd;
  Dd_inv = 1.0f / multiAMcoef->Mmunu.Dd;

  /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
   * therefore there is a factor of 2 difference with respect to the equations in JKS, which
   * where based on the single-sided PSD.
   */
  retF.F = Dd_inv * (  Bd * (SQ(retF.Fa.re) + SQ(retF.Fa.im) )
		       + Ad * ( SQ(retF.Fb.re) + SQ(retF.Fb.im) )
		       - 2.0 * Cd *( retF.Fa.re * retF.Fb.re + retF.Fa.im * retF.Fb.im )
		       );

  Fstat->F = (REAL8) retF.F;
  Fstat->Fa.re = (REAL8) retF.Fa.re;
  Fstat->Fa.im = (REAL8) retF.Fa.im;
  Fstat->Fb.re = (REAL8) retF.Fb.re;
  Fstat->Fb.im = (REAL8) retF.Fb.im;

  /* free memory if no buffer was available */
  if ( !cfBuffer )
    {
      XLALDestroyMultiSSBtimes_REAL4 ( multiSSB4 );
      XLALDestroyMultiAMCoeffs ( multiAMcoef );
    } /* if !cfBuffer */

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStat_REAL4() */


/** Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 *
 * Note: this is a single-precision version aimed for GPU parallelization.
 *
 */
int
XLALComputeFaFb_REAL4 ( Fcomponents_REAL4 *FaFb,
                        const SFTVector *sfts,
                        const PulsarSpins_REAL4 spins,
                        const SSBtimes_REAL4 *tSSB,
                        const AMCoeffs *amcoe,
                        UINT4 Dterms)
{
  static const char *fn = "XLALComputeFaFb_REAL4()";

  UINT4 alpha;                 	/* loop index over SFTs */
  UINT4 spdnOrder;		/* maximal spindown-orders */
  UINT4 numSFTs;		/* number of SFTs (M in the Notes) */
  COMPLEX8 Fa, Fb;
  REAL4 Tsft; 			/* length of SFTs in seconds */
  INT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  INT4 freqIndex1;		/* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;		/* pointer to alpha-arrays over a and b */

  REAL4 *DeltaT_int_al, *DeltaT_rem_al, *TdotM1_al;	/* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;		/* SFT alpha  */

  REAL4 norm = OOTWOPI_FLOAT;

  /* ----- check validity of input */
#ifndef LAL_NDEBUG
  if ( !FaFb ) {
    XLALPrintError ("%s: Output-pointer is NULL !\n", fn);
    XLAL_ERROR ( fn, XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    XLALPrintError ("%s: Input SFTs are NULL!\n", fn);
    XLAL_ERROR ( fn, XLAL_EINVAL);
  }

  if ( !tSSB || !tSSB->DeltaT_int || !tSSB->DeltaT_rem || !tSSB->TdotM1 || !amcoe || !amcoe->a || !amcoe->b )
    {
      XLALPrintError ("%s: Illegal NULL in input !\n", fn);
      XLAL_ERROR ( fn, XLAL_EINVAL);
    }
#endif

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = (REAL4)( 1.0 / sfts->data[0].deltaF );

  REAL4 dFreq = sfts->data[0].deltaF;
  freqIndex0 = (UINT4) ( sfts->data[0].f0 / dFreq + 0.5); /* lowest freqency-index */
  freqIndex1 = freqIndex0 + sfts->data[0].data->length;

  REAL4 f0 = spins.FreqMain;
  REAL4 df = spins.fkdot[0];

  /* find highest non-zero spindown-entry */
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( spins.fkdot[spdnOrder] )
      break;

  Fa.re = 0.0f;
  Fa.im = 0.0f;
  Fb.re = 0.0f;
  Fb.im = 0.0f;

  /* convenient shortcuts, pointers to beginning of alpha-arrays */
  a_al = amcoe->a->data;
  b_al = amcoe->b->data;
  DeltaT_int_al = tSSB->DeltaT_int->data;
  DeltaT_rem_al = tSSB->DeltaT_rem->data;
  TdotM1_al = tSSB->TdotM1->data;

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

      REAL4 lambda_alpha;
      REAL4 kappa_max, kappa_star;

      /* ----- calculate kappa_max and lambda_alpha */
      {
	UINT4 s; 		/* loop-index over spindown-order */

	REAL4 Tas;		/* temporary variable to calculate (DeltaT_alpha)^s */
        REAL4 T0, dT, deltaT;
        REAL4 Dphi_alpha_int, Dphi_alpha_rem;
        REAL4 phi_alpha_rem;
        REAL4 TdotM1;		/* defined as Tdot_al - 1 */

        /* 1st oder: s = 0 */
        TdotM1 = *TdotM1_al;
        T0 = *DeltaT_int_al;
        dT = *DeltaT_rem_al;
        deltaT = T0 + dT;

        /* phi_alpha = f * Tas; */
        phi_alpha_rem = f0 * dT + df * T0 + df * dT;
        phi_alpha_rem = REM ( phi_alpha_rem );
        Dphi_alpha_int = f0;
        Dphi_alpha_rem = df * (1.0f + TdotM1) + f0 * TdotM1;

        /* higher-order spindowns */
        Tas = deltaT;
	for (s=1; s <= spdnOrder; s++)
	  {
	    REAL4 fsdot = spins.fkdot[s];
	    Dphi_alpha_rem += fsdot * Tas * inv_fact[s]; 	/* here: Tas = DT^s */
	    Tas *= deltaT;					/* now:  Tas = DT^(s+1) */
	    phi_alpha_rem += fsdot * Tas * inv_fact[s+1];
	  } /* for s <= spdnOrder */

	/* Step 3: apply global factor of Tsft to complete Dphi_alpha */
        Dphi_alpha_int *= Tsft;
	Dphi_alpha_rem *= Tsft;

        REAL4 tmp = REM( 0.5f * Dphi_alpha_int ) + REM ( 0.5f * Dphi_alpha_rem );
	lambda_alpha = phi_alpha_rem - tmp;

	/* real- and imaginary part of e^{-i 2 pi lambda_alpha } */
	if ( sin_cos_2PI_LUT ( &imagQ, &realQ, - lambda_alpha ) ) {
          XLALPrintError ("%s: sin_cos_2PI_LUT() failed.!\n", fn );
	  XLAL_ERROR ( fn, XLAL_EFUNC);
	}

        kstar = (INT4)Dphi_alpha_int + (INT4)Dphi_alpha_rem;
	kappa_star = REM(Dphi_alpha_int) + REM(Dphi_alpha_rem);
	kappa_max = kappa_star + 1.0f * Dterms - 1.0f;

	/* ----- check that required frequency-bins are found in the SFTs ----- */
	k0 = kstar - Dterms + 1;
	k1 = k0 + 2 * Dterms - 1;
	if ( (k0 < freqIndex0) || (k1 > freqIndex1) )
	  {
	    XLALPrintError ("%s: Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n\n",
                            fn, k0, k1, freqIndex0, freqIndex1 );
	    XLAL_ERROR( fn, XLAL_EDOM);
	  }

      } /* compute kappa_star, lambda_alpha */

      /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
       * the trig-functions need to be calculated only once!
       * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the
       * closest to zero and will pose no numerical difficulties !
       */
      sin_cos_2PI_LUT ( &s_alpha, &c_alpha, kappa_star );
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
	  REAL4 Sn = (*Xalpha_l).re;
	  REAL4 Tn = (*Xalpha_l).im;
	  REAL4 pn = kappa_max;
	  REAL4 qn = pn;
	  REAL4 U_alpha, V_alpha;

	  /* recursion with 2*Dterms steps */
	  UINT4 l;
	  for ( l = 1; l < 2*Dterms; l ++ )
	    {
	      Xalpha_l ++;

	      pn = pn - 1.0f; 			/* p_(n+1) */
	      Sn = pn * Sn + qn * (*Xalpha_l).re;	/* S_(n+1) */
	      Tn = pn * Tn + qn * (*Xalpha_l).im;	/* T_(n+1) */
	      qn *= pn;				/* q_(n+1) */
	    } /* for l <= 2*Dterms */

	  U_alpha = Sn / qn;
	  V_alpha = Tn / qn;

#ifndef LAL_NDEBUG
	  if ( !finite(U_alpha) || !finite(V_alpha) || !finite(pn) || !finite(qn) || !finite(Sn) || !finite(Tn) ) {
	    XLAL_ERROR (fn, XLAL_EFPINVAL );
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
	  realXP = TWOPI_FLOAT * Xalpha_l[ind0].re;
	  imagXP = TWOPI_FLOAT * Xalpha_l[ind0].im;
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

      DeltaT_int_al ++;
      DeltaT_rem_al ++;
      TdotM1_al ++;

      SFT_al ++;

    } /* for alpha < numSFTs */

  /* return result */
  FaFb->Fa.re = norm * Fa.re;
  FaFb->Fa.im = norm * Fa.im;
  FaFb->Fb.re = norm * Fb.re;
  FaFb->Fb.im = norm * Fb.im;

  return XLAL_SUCCESS;

} /* XLALComputeFaFb_REAL4() */



/** Destroy a MultiSSBtimes_REAL4 structure.
 *  Note, this is "NULL-robust" in the sense that it will not crash
 *  on NULL-entries anywhere in this struct, so it can be used
 *  for failure-cleanup even on incomplete structs
 */
void
XLALDestroyMultiSSBtimes_REAL4 ( MultiSSBtimes_REAL4 *multiSSB )
{
  UINT4 X;

  if ( ! multiSSB )
    return;

  if ( multiSSB->data )
    {
      for ( X=0; X < multiSSB->length; X ++ )
	{
          XLALDestroySSBtimes_REAL4 ( multiSSB->data[X] );
	} /* for X < numDetectors */
      LALFree ( multiSSB->data );
    }

  LALFree ( multiSSB );

  return;

} /* XLALDestroyMultiSSBtimes_REAL4() */


/** Destroy a SSBtimes_REAL4 structure.
 *  Note, this is "NULL-robust" in the sense that it will not crash
 *  on NULL-entries anywhere in this struct, so it can be used
 *  for failure-cleanup even on incomplete structs
 */
void
XLALDestroySSBtimes_REAL4 ( SSBtimes_REAL4 *tSSB )
{
  if ( ! tSSB )
    return;

  if ( tSSB->DeltaT_int )
    XLALDestroyREAL4Vector ( tSSB->DeltaT_int );

  if ( tSSB->DeltaT_rem )
    XLALDestroyREAL4Vector ( tSSB->DeltaT_rem );

  if ( tSSB->TdotM1 )
    XLALDestroyREAL4Vector ( tSSB->TdotM1 );

  LALFree ( tSSB );

  return;

} /* XLALDestroySSBtimes_REAL4() */


/** Multi-IFO version of LALGetSSBtimes_REAL4().
 * Get all SSB-timings for all input detector-series in REAL4 representation.
 *
 */
MultiSSBtimes_REAL4 *
XLALGetMultiSSBtimes_REAL4 ( const MultiDetectorStateSeries *multiDetStates, 	/**< [in] detector-states at timestamps t_i */
                             REAL8 Alpha, REAL8 Delta,				/**< source sky-location, in equatorial coordinates  */
                             LIGOTimeGPS refTime
                             )
{
  static const char *fn = "XLALGetMultiSSBtimes_REAL4()";

  UINT4 X, numDetectors;
  MultiSSBtimes_REAL4 *ret;

  /* check input */
  if ( !multiDetStates || multiDetStates->length==0 ) {
    XLALPrintError ("%s: illegal NULL or empty input 'multiDetStates'.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  numDetectors = multiDetStates->length;

  if ( ( ret = XLALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    XLALPrintError ("%s: XLALCalloc( 1, %d ) failed.\n", fn, sizeof( *ret ) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  ret->length = numDetectors;
  if ( ( ret->data = XLALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    XLALFree ( ret );
    XLALPrintError ("%s: XLALCalloc( %d, %d ) failed.\n", fn, numDetectors, sizeof( *ret->data ) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      if ( (ret->data[X] = XLALGetSSBtimes_REAL4 (multiDetStates->data[X], Alpha, Delta, refTime)) == NULL ) {
        XLALPrintError ("%s: XLALGetSSBtimes_REAL4() failed. xlalErrno = %d\n", fn, xlalErrno );
        goto failed;
      }

    } /* for X < numDet */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiSSBtimes_REAL4 ( ret );
  XLAL_ERROR_NULL (fn, XLAL_EFAILED );

 success:
  return ret;

} /* LALGetMultiSSBtimes() */



/** XLAL REAL4-version of LALGetSSBtimes()
 *
 */
SSBtimes_REAL4 *
XLALGetSSBtimes_REAL4 ( const DetectorStateSeries *DetectorStates,	/**< [in] detector-states at timestamps t_i */
                        REAL8 Alpha, REAL8 Delta,			/**< source sky-location, in equatorial coordinates  */
                        LIGOTimeGPS refTime
                        )
{
  static const char *fn = "XLALGetSSBtimes_REAL4()";

  UINT4 numSteps, i;
  REAL8 refTimeREAL8;
  SSBtimes_REAL4 *ret;

  /* check input consistency */
  if ( !DetectorStates || DetectorStates->length==0 ) {
    XLALPrintError ("%s: illegal NULL or empty input 'DetectorStates'.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  numSteps = DetectorStates->length;		/* number of timestamps */

  /* convenience variables */
  refTimeREAL8 = XLALGPSGetREAL8( &refTime );

  /* allocate return container */
  if ( ( ret = XLALMalloc( sizeof(*ret))) == NULL ) {
    XLALPrintError ("%s: XLALMalloc(%d) failed.\n", fn, sizeof(*ret) );
    goto failed;
  }
  if ( (ret->DeltaT_int  = XLALCreateREAL4Vector ( numSteps )) == NULL ||
       (ret->DeltaT_rem  = XLALCreateREAL4Vector ( numSteps )) == NULL ||
       (ret->TdotM1      = XLALCreateREAL4Vector ( numSteps )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateREAL4Vector ( %d ) failed.\n", fn, numSteps );
      goto failed;
    }

  /* loop over all SFTs and compute timing info */
  for (i=0; i < numSteps; i++ )
    {
      BarycenterInput baryinput = empty_BarycenterInput;
      EmissionTime emit;
      DetectorState *state = &(DetectorStates->data[i]);
      LALStatus status = empty_LALStatus;
      REAL8 deltaT;
      REAL4 deltaT_int;

      baryinput.tgps = state->tGPS;
      baryinput.site = DetectorStates->detector;
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;

      baryinput.alpha = Alpha;
      baryinput.delta = Delta;
      baryinput.dInv = 0;

      LALBarycenter(&status, &emit, &baryinput, &(state->earthState) );
      if ( status.statusCode ) {
        XLALPrintError ("%s: LALBarycenter() failed with status = %d, '%s'\n", fn, status.statusCode, status.statusDescription );
        goto failed;
      }

      deltaT = XLALGPSGetREAL8 ( &emit.te ) - refTimeREAL8;
      deltaT_int = (INT4)deltaT;

      ret->DeltaT_int->data[i] = deltaT_int;
      ret->DeltaT_rem->data[i]  = (REAL4)( deltaT - (REAL8)deltaT_int );
      ret->TdotM1->data[i] = (REAL4) ( emit.tDot - 1.0 );

    } /* for i < numSteps */

  /* finally: store the reference-time used into the output-structure */
  ret->refTime = refTime;

  goto success;

 failed:
  XLALDestroySSBtimes_REAL4 ( ret );
  XLAL_ERROR_NULL ( fn, XLAL_EFAILED );

 success:
  return ret;

} /* XLALGetSSBtimes_REAL4() */



/** Destruction of a ComputeFBuffer_REAL4 *contents*,
 * i.e. the multiSSB and multiAMcoeff, while the
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...) */
void
XLALEmptyComputeFBuffer_REAL4 ( ComputeFBuffer_REAL4 *cfb)
{
  if ( !cfb )
    return;

  XLALDestroyMultiSSBtimes_REAL4 ( cfb->multiSSB );
  cfb->multiSSB = NULL;

  XLALDestroyMultiAMCoeffs ( cfb->multiAMcoef );
  cfb->multiAMcoef = NULL;

  return;

} /* XLALDestroyComputeFBuffer_REAL4() */

