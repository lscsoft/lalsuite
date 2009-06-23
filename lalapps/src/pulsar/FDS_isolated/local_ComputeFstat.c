/*
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
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/FindRoot.h>

/* GSL includes */
#include <lal/LALGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>


#include <lal/AVFactories.h>
#include <lal/ComputeFstat.h>
#include <lal/ComplexAM.h>

NRCSID( COMPUTEFSTATC, "$Id$");

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)


#define LD_SMALL4       (2.0e-4)		/**< "small" number for REAL4*/
#define OOTWOPI         (1.0 / LAL_TWOPI)	/**< 1/2pi */

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */

#define EA_ACC          1E-9                    /* the timing accuracy of LALGetBinaryTimes in seconds */

/*----- Macros ----- */
#define SQ(x) ( (x) * (x) )
#define REM(x) ( (x) - (INT4)(x) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
#define NUM_FACT 7
static const REAL4 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* empty initializers  */
static const LALStatus empty_status;
static const AMCoeffs empty_AMCoeffs;

const SSBtimes empty_SSBtimes;
const MultiSSBtimes empty_MultiSSBtimes;
const AntennaPatternMatrix empty_AntennaPatternMatrix;
const MultiAMCoeffs empty_MultiAMCoeffs;
const Fcomponents empty_Fcomponents;
const ComputeFParams empty_ComputeFParams;
const ComputeFBuffer empty_ComputeFBuffer;

/*---------- internal prototypes ----------*/
int finite(double x);

int local_XLALComputeFaFb ( Fcomponents *FaFb,
                            const SFTVector *sfts,
                            const PulsarSpins fkdot,
                            const SSBtimes *tSSB,
                            const AMCoeffs *amcoe,
                            const ComputeFParams *params);


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
local_ComputeFStat ( LALStatus *status,
	       Fcomponents *Fstat,                 		/**< [out] Fstatistic + Fa, Fb */
	       const PulsarDopplerParams *doppler, 		/**< parameter-space point to compute F for */
	       const MultiSFTVector *multiSFTs,    		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
	       const MultiNoiseWeights *multiWeights,		/**< noise-weights of all SFTs */
	       const MultiDetectorStateSeries *multiDetStates,	/**< 'trajectories' of the different IFOs */
	       const ComputeFParams *params,       		/**< addition computational params */
	       ComputeFBuffer *cfBuffer            		/**< CF-internal buffering structure */
	       )
{
  Fcomponents retF = empty_Fcomponents;
  UINT4 X, numDetectors;
  MultiSSBtimes *multiSSB = NULL;
  MultiSSBtimes *multiBinary = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  MultiCmplxAMCoeffs *multiCmplxAMcoef = NULL;
  REAL8 Ad, Bd, Cd, Dd_inv, Ed;
  SkyPosition skypos;

  INITSTATUS( status, "ComputeFStat", COMPUTEFSTATC );
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT ( Fstat, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiSFTs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  numDetectors = multiSFTs->length;
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  if ( multiWeights ) {
    ASSERT ( multiWeights->length == numDetectors , status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  }

  /* check if that skyposition SSB+AMcoef were already buffered */
  if ( cfBuffer
       && ( cfBuffer->multiDetStates == multiDetStates )
       && ( cfBuffer->Alpha == doppler->Alpha )
       && ( cfBuffer->Delta == doppler->Delta )
       && cfBuffer->multiSSB )
    { /* yes ==> reuse */
      multiSSB = cfBuffer->multiSSB;

      /* re-use (LWL) AM coefficients whenever available */
      if ( cfBuffer->multiAMcoef )
	multiAMcoef = cfBuffer->multiAMcoef;

      /* re-use RAA AM coefficients *only* if bufferedRAA is TRUE !*/
      if ( params->bufferedRAA && cfBuffer->multiCmplxAMcoef  )
	multiCmplxAMcoef = cfBuffer->multiCmplxAMcoef;

    } /* if have buffered stuff to reuse */
  else
    {
      skypos.system =   COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;
      TRY ( LALGetMultiSSBtimes ( status->statusPtr, &multiSSB, multiDetStates, skypos, doppler->refTime, params->SSBprec ), status );
      if ( cfBuffer )
	{
	  XLALDestroyMultiSSBtimes ( cfBuffer->multiSSB );
	  cfBuffer->multiSSB = multiSSB;
	  cfBuffer->Alpha = doppler->Alpha;
	  cfBuffer->Delta = doppler->Delta;
	  cfBuffer->multiDetStates = multiDetStates ;
	} /* buffer new SSB times */

    } /* could not reuse previously buffered quantites */

    /* new orbital parameter corrections if not already buffered */
  if ( doppler->orbit )
    {
      /* if already buffered */
      if ( cfBuffer && cfBuffer->multiBinary )
	{ /* yes ==> reuse */
	  multiBinary = cfBuffer->multiBinary;
	}
      else
	{
	  /* compute binary time corrections to the SSB time delays and SSB time derivitive */
	  TRY ( LALGetMultiBinarytimes ( status->statusPtr, &multiBinary, multiSSB, multiDetStates, doppler->orbit, doppler->refTime ), status );

	  /* store these in buffer if available */
	  if ( cfBuffer )
	    {
	      XLALDestroyMultiSSBtimes ( cfBuffer->multiBinary );
	      cfBuffer->multiBinary = multiBinary;
	    } /* if cfBuffer */
 	}
    }
  else multiBinary = multiSSB;
  /*
  printf("multiSSB = %6.12f %6.12f multiBinary = %6.12f %6.12f\n",
	 multiSSB->data[0]->DeltaT->data[0],
	 multiSSB->data[0]->Tdot->data[0],
	 multiBinary->data[0]->DeltaT->data[0],
	 multiBinary->data[0]->Tdot->data[0]);
  */

  /* special treatment of AM coefficients */
  if ( params->useRAA && !multiCmplxAMcoef )
    {
      /* compute new RAA AM-coefficients */
      LALGetMultiCmplxAMCoeffs ( status->statusPtr, &multiCmplxAMcoef, multiDetStates, *doppler );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      /* noise-weight Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiCmplxAMCoeffs ( multiCmplxAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	LALPrintError("\nXLALWeighMultiCmplxAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
      }

      /* store in buffer if available */
      if ( cfBuffer )
	{
	  XLALDestroyMultiCmplxAMCoeffs ( cfBuffer->multiCmplxAMcoef );
	  cfBuffer->multiCmplxAMcoef = multiCmplxAMcoef;
	}

    } /* if RAA AM coefficients need to be computed */

  if ( !params->useRAA && !multiAMcoef )
    {
      /* compute new AM-coefficients */
      LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, skypos );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      /* noise-weight Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	LALPrintError("\nXLALWeighMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
      }

      /* store these in buffer if available */
      if ( cfBuffer )
	{
	  XLALDestroyMultiAMCoeffs ( cfBuffer->multiAMcoef );
	  cfBuffer->multiAMcoef = multiAMcoef;
	} /* if cfBuffer */

    } /* if LWL AM coefficient need to be computed */

  if ( multiAMcoef )
    {
      Ad = multiAMcoef->Mmunu.Ad;
      Bd = multiAMcoef->Mmunu.Bd;
      Cd = multiAMcoef->Mmunu.Cd;
      Dd_inv = 1.0 / multiAMcoef->Mmunu.Dd;
      Ed = 0;
    }
  else if ( multiCmplxAMcoef )
    {
      Ad = multiCmplxAMcoef->Mmunu.Ad;
      Bd = multiCmplxAMcoef->Mmunu.Bd;
      Cd = multiCmplxAMcoef->Mmunu.Cd;
      Ed = multiCmplxAMcoef->Mmunu.Ed;
      Dd_inv = 1.0 / multiCmplxAMcoef->Mmunu.Dd;
    }
  else
    {
      LALPrintError ( "Programming error: neither 'multiAMcoef' nor 'multiCmplxAMcoef' are available!\n");
      ABORT ( status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
    }

  /* ----- loop over detectors and compute all detector-specific quantities ----- */
  for ( X=0; X < numDetectors; X ++)
    {
      Fcomponents FcX = empty_Fcomponents;	/* for detector-specific FaX, FbX */

      if ( params->useRAA )
	{
	  if ( XLALComputeFaFbCmplx (&FcX, multiSFTs->data[X], doppler->fkdot, multiSSB->data[X], multiCmplxAMcoef->data[X], params) != 0)
	    {
	      LALPrintError ("\nXALComputeFaFbCmplx() failed\n");
	      ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
	    }
	}
      else if ( params->upsampling > 1)
	{
	  if ( XLALComputeFaFbXavie (&FcX, multiSFTs->data[X], doppler->fkdot, multiBinary->data[X], multiAMcoef->data[X], params) != 0)
	    {
	      LALPrintError ("\nXALComputeFaFbXavie() failed\n");
	      ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
	    }
	}
      else
	{
	  if ( local_XLALComputeFaFb (&FcX, multiSFTs->data[X], doppler->fkdot, multiBinary->data[X], multiAMcoef->data[X], params) != 0)
	    {
	      LALPrintError ("\nXALComputeFaFb() failed\n");
	      ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
	    }
	}

#ifndef LAL_NDEBUG
      if ( !finite(FcX.Fa.re) || !finite(FcX.Fa.im) || !finite(FcX.Fb.re) || !finite(FcX.Fb.im) ) {
	LALPrintError("XLALComputeFaFb() returned non-finite: Fa=(%f,%f), Fb=(%f,%f)\n",
		      FcX.Fa.re, FcX.Fa.im, FcX.Fb.re, FcX.Fb.im );
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

  /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
   * therefore there is a factor of 2 difference with respect to the equations in JKS, which
   * where based on the single-sided PSD.
   */
  retF.F = Dd_inv * (  Bd * (SQ(retF.Fa.re) + SQ(retF.Fa.im) )
		       + Ad * ( SQ(retF.Fb.re) + SQ(retF.Fb.im) )
		       - 2.0 * Cd *( retF.Fa.re * retF.Fb.re + retF.Fa.im * retF.Fb.im )
		       );
  if ( Ed != 0 ) /* extra term in RAA case */
    retF.F += - 2.0 * Dd_inv * Ed *( - retF.Fa.re * retF.Fb.im + retF.Fa.im * retF.Fb.re ); /* -2 E Im(Fa Fb^* ) / D */

  (*Fstat) = retF;

  /* free memory if no buffer was available */
  if ( !cfBuffer )
    {
      XLALDestroyMultiSSBtimes ( multiSSB );
      XLALDestroyMultiSSBtimes ( multiBinary );
      XLALDestroyMultiAMCoeffs ( multiAMcoef );
      XLALDestroyMultiCmplxAMCoeffs ( multiCmplxAMcoef );
    } /* if !cfBuffer */

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStat() */


/** Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 *
 * Note: this is a locally optimized version of XLALComputeFaFb(), targeted at GPU parallelization.
 *
 */
int
local_XLALComputeFaFb ( Fcomponents *FaFb,
		  const SFTVector *sfts,
		  const PulsarSpins fkdot,
		  const SSBtimes *tSSB,
		  const AMCoeffs *amcoe,
		  const ComputeFParams *params)       /**< addition computational params */
{
  UINT4 alpha;                 	/* loop index over SFTs */
  UINT4 spdnOrder;		/* maximal spindown-orders */
  UINT4 numSFTs;		/* number of SFTs (M in the Notes) */
  COMPLEX16 Fa, Fb;
  REAL4 Tsft; 			/* length of SFTs in seconds */
  INT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  INT4 freqIndex1;		/* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;		/* pointer to alpha-arrays over a and b */
  REAL8 *DeltaT_al, *Tdot_al;	/* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;		/* SFT alpha  */
  UINT4 Dterms = params->Dterms;

  REAL4 norm = OOTWOPI;

  REAL4 f0, df;

  f0 = (INT4)fkdot[0];
  df = fkdot[0] - (REAL8)f0;

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

  if ( !tSSB || !tSSB->DeltaT || !tSSB->Tdot || !amcoe || !amcoe->a || !amcoe->b || !params)
    {
      LALPrintError ("\nIllegal NULL in input !\n\n");
      XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
    }

  if ( PULSAR_MAX_SPINS > NUM_FACT )
    {
      LALPrintError ("\nInverse factorials table only up to order s=%d, can't handle %d spin-order\n\n",
		     NUM_FACT, PULSAR_MAX_SPINS - 1 );
      XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
    }
#endif

  if ( params->upsampling > 1 ) {
    fprintf (stderr, "\n===== WARNING: XLALComputeFaFb() should not be used with upsampled-SFTs!\n");
    XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
  }

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = (REAL4)( 1.0 / sfts->data[0].deltaF );

  REAL4 dFreq = sfts->data[0].deltaF;
  freqIndex0 = (UINT4) ( sfts->data[0].f0 / dFreq + 0.5); /* lowest freqency-index */
  freqIndex1 = freqIndex0 + sfts->data[0].data->length;

  /* find highest non-zero spindown-entry */
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot[spdnOrder] )
      break;

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

      REAL4 lambda_alpha;
      REAL4 kappa_max, kappa_star;

      /* ----- calculate kappa_max and lambda_alpha */
      {
	REAL8 DT_al = (*DeltaT_al);

	UINT4 s; 		/* loop-index over spindown-order */

	REAL4 Tas;	/* temporary variable to calculate (DeltaT_alpha)^s */
        REAL4 T0, dT;
        REAL4 Dphi_alpha_int, Dphi_alpha_rem;
        REAL4 phi_alpha_int, phi_alpha_rem;
        REAL4 Tdot_al_frac;	/* defined as Tdot_al - 1, *not* the remainder wrt floor ! */

        Tdot_al_frac = (*Tdot_al) - 1.0f;

        /* 1st oder: s = 0 */
        T0 = (INT4)DT_al;
        dT = DT_al - (REAL8)T0;

        /* phi_alpha = f * Tas; */
        phi_alpha_rem = f0 * dT + df * T0 + df * dT;
        phi_alpha_rem = REM ( phi_alpha_rem );
        Dphi_alpha_int = f0;
        Dphi_alpha_rem = df * (1.0f + Tdot_al_frac) + f0 * Tdot_al_frac;

        /* higher-order spindowns */
        Tas = DT_al;
	for (s=1; s <= spdnOrder; s++)
	  {
	    REAL4 fsdot = fkdot[s];
	    Dphi_alpha_rem += fsdot * Tas * inv_fact[s]; 	/* here: DT^s/s! */
	    Tas *= DT_al;					/* now: DT^(s+1) */
	    phi_alpha_rem += fsdot * Tas * inv_fact[s+1];
	  } /* for s <= spdnOrder */

	/* Step 3: apply global factor of Tsft to complete Dphi_alpha */
        Dphi_alpha_int *= Tsft;
	Dphi_alpha_rem *= Tsft;

        REAL4 tmp = REM( 0.5f * Dphi_alpha_int ) + REM ( 0.5f * Dphi_alpha_rem );
	lambda_alpha = phi_alpha_rem - tmp;

        kstar = (INT4)Dphi_alpha_int + (INT4)Dphi_alpha_rem;
	kappa_star = REM(Dphi_alpha_int) + REM(Dphi_alpha_rem);
	kappa_max = kappa_star + 1.0f * Dterms - 1.0f;

	/* real- and imaginary part of e^{-i 2 pi lambda_alpha } */
	if ( sin_cos_2PI_LUT ( &imagQ, &realQ, - lambda_alpha ) ) {
	  XLAL_ERROR ( "XLALComputeFaFb", XLAL_EFUNC);
	}

	/* ----- check that required frequency-bins are found in the SFTs ----- */
	k0 = kstar - Dterms + 1;
	k1 = k0 + 2 * Dterms - 1;
	if ( (k0 < freqIndex0) || (k1 > freqIndex1) )
	  {
	    LALPrintError ("Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n\n",
			   k0, k1, freqIndex0, freqIndex1 );
	    XLAL_ERROR("XLALComputeFaFb", XLAL_EDOM);
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
	    XLAL_ERROR ("XLALComputeFaFb()", COMPUTEFSTATC_EIEEE);
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
      DeltaT_al ++;
      Tdot_al ++;
      SFT_al ++;

    } /* for alpha < numSFTs */

  /* return result */
  FaFb->Fa.re = norm * Fa.re;
  FaFb->Fa.im = norm * Fa.im;
  FaFb->Fb.re = norm * Fb.re;
  FaFb->Fb.im = norm * Fb.im;

  return XLAL_SUCCESS;

} /* XLALComputeFaFb() */
