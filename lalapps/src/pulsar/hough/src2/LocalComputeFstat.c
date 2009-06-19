/*
*  Copyright (C) 2007 Bernd Machenschalk
*  Copyright (C) 2006 John T. Whelan, Badri Krishnan
*  Copyright (C) 2005, 2006 Reinhard Prix 
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

/*
 * Stripped together and modified from ComputeFStat.c in LAL
 * by Bernd Machenschalk for Einstein@Home
 * $Id$
 */


/*---------- INCLUDES ----------*/

#define __USE_ISOC99 1
#include <math.h>

#include <lal/ExtrapolatePulsarSpins.h>

#include <lal/AVFactories.h>
#include <lal/ComputeFstat.h>
#include <lal/LogPrintf.h>

#include "LocalOptimizationFlags.h"

#if defined(_MSC_VER) && (EAH_HOTLOOP_VARIANT == EAH_HOTLOOP_VARIANT_SSE)
#include "xmmintrin.h"
#endif

/* necessary for gcc on linux-powerpc, shouldn't hurt on Macs */
#if (EAH_HOTLOOP_VARIANT == EAH_HOTLOOP_VARIANT_ALTIVEC)
#include <altivec.h>
#endif

NRCSID( LOCALCOMPUTEFSTATC, "$Id$");


/*---------- local DEFINES ----------*/
#define TRUE  (1==1)
#define FALSE (1==0)

#define LD_SMALL4       (2.0e-4)		/**< "small" number for REAL4*/ 
#define OOTWOPI         (1.0 / LAL_TWOPI)	/**< 1/2pi */

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */ 

/* don't use finite() from win_lib.cpp here for performance reasons */
#ifdef _MSC_VER
#define finite _finite
#endif

/*---------- optimization dependant switches ----------*/


/*----- Macros ----- */
/** fixed DTERMS to allow for loop unrolling */
#define DTERMS 8

/** square */
#define SQ(x) ( (x) * (x) )

#ifndef __GNUC__
/** somehow the branch prediction of gcc-4.1.2 terribly failes
    with the current case distinction in the hot-loop,
    having a severe impact on runtime of the E@H Linux App.
    So let's allow to give gcc a hint which path has a higher probablility */
#define __builtin_expect(a,b) a

/** currently interleaving the kernel loop doesn't work with MSC */
#undef EAH_HOTLOOP_INTERLEAVED
#endif


#if EAH_HOUGH_PREFETCH > EAH_HOUGH_PREFETCH_NONE
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
#include "xmmintrin.h"
#define PREFETCH(a) _mm_prefetch((char *)(void *)(a),_MM_HINT_T0)
#elif defined(__GNUC__)
#define PREFETCH(a) __builtin_prefetch(a)
#else
#define PREFETCH(a) a
#endif
#else
#define PREFETCH(a) a
#endif

#include "EinsteinAtHome/sincos.ci"

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
#define NUM_FACT 6
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0) };

/* empty initializers  */
static const LALStatus empty_status;


/*---------- internal prototypes ----------*/
static void
LocalComputeFStat ( LALStatus*, Fcomponents*, const PulsarDopplerParams*,
		    const MultiSFTVector*, const MultiNoiseWeights*,
		    const MultiDetectorStateSeries*, const ComputeFParams*,
		    ComputeFBuffer*);

static int
LocalXLALComputeFaFb (Fcomponents*, const SFTVector*, const PulsarSpins,
		      const SSBtimes*, const AMCoeffs*, const ComputeFParams*);


/*==================== FUNCTION DEFINITIONS ====================*/


/** Function to compute a vector of Fstatistic values for a number of frequency bins.
    This function is simply a wrapper for LocalComputeFstat() which is called repeatedly for
    every frequency value.  The output, i.e. fstatVector must be properly allocated
    before this function is called.  The values of the start frequency, the step size
    in the frequency and the number of frequency values for which the Fstatistic is 
    to be calculated are read from fstatVector.  The other parameters are not checked and 
    they must be correctly set outside this function. 
*/
void LocalComputeFStatFreqBand ( LALStatus *status, 
				 REAL8FrequencySeries *fstatVector, /**< [out] Vector of Fstat values */
				 const PulsarDopplerParams *doppler,/**< parameter-space point to compute F for */
				 const MultiSFTVector *multiSFTs, /**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
				 const MultiNoiseWeights *multiWeights,	/**< noise-weights of all SFTs */
				 const MultiDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
				 const ComputeFParams *params	/**< addition computational params */
				 )
{

  UINT4 numDetectors, numBins, k;	
  REAL8 deltaF;
  Fcomponents Fstat;
  PulsarDopplerParams thisPoint;
  ComputeFBuffer cfBuffer = empty_ComputeFBuffer;

  INITSTATUS( status, "LocalComputeFStatFreqBand", LOCALCOMPUTEFSTATC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( multiSFTs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  numDetectors = multiSFTs->length;
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  ASSERT ( fstatVector, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( fstatVector->data, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( fstatVector->data->data, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( fstatVector->data->length > 0, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  ASSERT ( DTERMS == params->Dterms, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  {
    static int firstcall = TRUE;
    if (firstcall) {
      /* init sin/cos lookup tables */
      local_sin_cos_2PI_LUT_init();

      /* make sure Dterms is what we expect */
      if (DTERMS != params->Dterms) {
	LogPrintf(LOG_CRITICAL, "LocalComputeFstat has been compiled with fixed DTERMS (%d) != params->Dtems (%d)\n",DTERMS, params->Dterms);
	ABORT ( status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
      }

      /* write out optimization settings */
      fprintf(stderr,
	      "\n$Revision$ REV:%s, OPT:%d, "
	      "SCVAR:%d, SCTRIM:%d, "
	      "HOTVAR:%d, HOTDIV:%d, "
	      "HGHPRE:%d, HGHBAT:%d\n",
	      EAH_OPTIMIZATION_REVISION, EAH_OPTIMIZATION,
	      EAH_SINCOS_VARIANT,  EAH_SINCOS_ROUND,
	      EAH_HOTLOOP_VARIANT, EAH_HOTLOOP_DIVS,
	      EAH_HOUGH_PREFETCH,  EAH_HOUGH_BATCHSIZE_LOG2);
      firstcall = FALSE;
    }
  }

  /** something to improve/cleanup -- the start frequency is available both 
      from the fstatvector and from the input doppler point -- they could be inconsistent
      or the user of this function could misunderstand */

  /* copy values from 'doppler' to local variable 'thisPoint' */
  thisPoint = *doppler;

  numBins = fstatVector->data->length;
  deltaF = fstatVector->deltaF;

  /* loop over frequency values and fill up values in fstatVector */
  for ( k = 0; k < numBins; k++) {
 
    TRY (LocalComputeFStat ( status->statusPtr, &Fstat, &thisPoint, multiSFTs, multiWeights, 
			     multiDetStates, params, &cfBuffer ), status);

    fstatVector->data->data[k] = Fstat.F;
      
    thisPoint.fkdot[0] += deltaF;
  }

  XLALEmptyComputeFBuffer ( &cfBuffer );

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStatFreqBand() */




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
static void
LocalComputeFStat ( LALStatus *status, 
		    Fcomponents *Fstat,                 /**< [out] Fstatistic + Fa, Fb */
		    const PulsarDopplerParams *doppler, /**< parameter-space point to compute F for */
		    const MultiSFTVector *multiSFTs,    /**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
		    const MultiNoiseWeights *multiWeights,/**< noise-weights of all SFTs */
		    const MultiDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
		    const ComputeFParams *params,       /**< addition computational params */
		    ComputeFBuffer *cfBuffer            /**< CF-internal buffering structure */
		    )
{
  Fcomponents retF = empty_Fcomponents;
  UINT4 X, numDetectors;	
  MultiSSBtimes *multiSSB = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  REAL8 Ad, Bd, Cd, Dd_inv;

  INITSTATUS( status, "LocalComputeFStat", LOCALCOMPUTEFSTATC );
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

  if ( doppler->orbit ) {
    XLALPrintError ("\nSorry, binary-pulsar search not yet implemented in LALComputeFStat()\n\n");
    ABORT ( status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
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
      skypos.system =   COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;
      /* compute new AM-coefficients and SSB-times */
      TRY ( LALGetMultiSSBtimes ( status->statusPtr, &multiSSB, multiDetStates, skypos, doppler->refTime, params->SSBprec ), status );

      LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, skypos );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      /* noise-weigh Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	XLALPrintError("\nXLALWeighMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
      }

      /* store these in buffer if available */
      if ( cfBuffer )
	{
	  XLALEmptyComputeFBuffer ( cfBuffer );
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
	      XLALPrintError ("\nXALComputeFaFbXavie() failed\n");
	      ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
	    }
	}
      else
	{
	  if ( LocalXLALComputeFaFb (&FcX, multiSFTs->data[X], doppler->fkdot, multiSSB->data[X], multiAMcoef->data[X], params) != 0)
	    {
	      XLALPrintError ("\nLocalXALComputeFaFb() failed\n");
	      ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
	    }
	}

#ifndef LAL_NDEBUG
      if ( !finite(FcX.Fa.re) || !finite(FcX.Fa.im) || !finite(FcX.Fb.re) || !finite(FcX.Fb.im) ) {
	XLALPrintError("LocalXLALComputeFaFb() returned non-finite: Fa=(%f,%f), Fb=(%f,%f)\n", 
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

  (*Fstat) = retF;

  /* free memory if no buffer was available */
  if ( !cfBuffer )
    {
      XLALDestroyMultiSSBtimes ( multiSSB );
      XLALDestroyMultiAMCoeffs ( multiAMcoef );
    } /* if !cfBuffer */

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LocalComputeFStat() */


/** Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 */
static int
LocalXLALComputeFaFb ( Fcomponents *FaFb,
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
  REAL8 f;			/* !! MUST be REAL8, or precision breaks down !! */
  REAL8 Tsft; 			/* length of SFTs in seconds */
  INT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  INT4 freqIndex1;		/* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;		/* pointer to alpha-arrays over a and b */
  REAL8 *DeltaT_al, *Tdot_al;	/* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;		/* SFT alpha  */
  
  REAL8 norm = OOTWOPI; 

#ifndef LAL_NDEBUG
  /* ----- check validity of input */
  if ( !FaFb ) {
    XLALPrintError ("\nOutput-pointer is NULL !\n\n");
    XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    XLALPrintError ("\nInput SFTs are NULL!\n\n");
    XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EINVAL);
  }
  
  if ( !tSSB || !tSSB->DeltaT || !tSSB->Tdot || !amcoe || !amcoe->a || !amcoe->b || !params)
    {
      XLALPrintError ("\nIllegal NULL in input !\n\n");
      XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EINVAL);
    }

  if ( PULSAR_MAX_SPINS > NUM_FACT )
    {
      XLALPrintError ("\nInverse factorials table only up to order s=%d, can't handle %d spin-order\n\n",
		     NUM_FACT, PULSAR_MAX_SPINS - 1 );
      XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EINVAL);
    }
#endif

  if ( params->upsampling > 1 ) {
    fprintf (stderr, "\n===== WARNING: LocalXLALComputeFaFb() should not be used with upsampled-SFTs!\n");
    XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EINVAL);
  }

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
    if ( fkdot[spdnOrder] != 0.0 )
      break;

  f = fkdot[0];

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
      REAL4 s_alpha, c_alpha;	/* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
      REAL4 realQ, imagQ;	/* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realXP, imagXP;     /* re/im of sum_k X_ak * P_ak */
      REAL4 realQXP, imagQXP;	/* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha, kappa_max, kappa_star;

      /* ----- calculate kappa_max and lambda_alpha */
      {
	UINT4 s; 		/* loop-index over spindown-order */
	REAL8 phi_alpha, Dphi_alpha, DT_al;
	REAL8 Tas;	/* temporary variable to calculate (DeltaT_alpha)^s */
        REAL8 TAS_invfact_s;
	static const REAL8 Dterms_1f = DTERMS - 1; /* Dterms - 1 as a double constant */

	/* init for s=0 */
	phi_alpha = 0.0;
	Dphi_alpha = 0.0;
	DT_al = (*DeltaT_al);
	Tas = 1.0;		/* DeltaT_alpha ^ 0 */
	TAS_invfact_s=1.0;    /* TAS / s! */

	for (s=0; s <= spdnOrder; s++) {
	  REAL8 fsdot = fkdot[s];
	  Dphi_alpha += fsdot * TAS_invfact_s; 	/* here: DT^s/s! */
#ifdef EAH_CHECK_FINITE_DPHI
	  if (!finite(Dphi_alpha)) {
	    LogPrintf(LOG_CRITICAL, "non-finite Dphi_alpha:%e, alpha:%d, spind#:%d, fkdot:%e, Tas:%e, inv_fact[s]:%e, inv_fact[s+1]:%e, phi_alpha:%e. DT_al:%e\n",
		      Dphi_alpha, alpha, s, fkdot[s], Tas, inv_fact[s], inv_fact[s+1], phi_alpha, DT_al);
	    XLAL_ERROR("LocalXLALComputeFaFb", XLAL_EDOM);
	  }
#endif
	  Tas *= DT_al;				/* now: DT^(s+1) */
          TAS_invfact_s= Tas * inv_fact[s+1];
	  phi_alpha += fsdot * TAS_invfact_s;
	} /* for s <= spdnOrder */

	/* Step 3: apply global factors to complete Dphi_alpha */
	Dphi_alpha *= Tsft * (*Tdot_al);		/* guaranteed > 0 ! */
#ifndef EAH_NO_CHECK_FINITE_DPHI
	if (!finite(Dphi_alpha)) {
	  LogPrintf(LOG_CRITICAL, "non-finite Dphi_alpha:%e, alpha:%d, Tsft:%e, Tkdot_al:%e Tas:%e, DT_al:%e\n",
		    Dphi_alpha, alpha, Tsft, (*Tdot_al), Tas, DT_al);
	  XLAL_ERROR("LocalXLALComputeFaFb", XLAL_EDOM);
	}
#endif
	lambda_alpha = phi_alpha - 0.5 * Dphi_alpha;
	
	/* FIXME: that might be possible to do faster */
	kstar = (INT4) (Dphi_alpha);	/* k* = floor(Dphi_alpha*chi) for positive Dphi */
	kappa_star = Dphi_alpha - 1.0 * kstar;	/* remainder of Dphi_alpha: >= 0 ! */
	kappa_max  = kappa_star + Dterms_1f;

	/* ----- check that required frequency-bins are found in the SFTs ----- */
	k0 = kstar - DTERMS + 1;
	/* 
	   original:
	   k1 = k0 + 2 * DTERMS - 1;
	   inserted k0:
	   k1 = kstar - DTERMS + 1 + 2 * DTERMS - 1;
	   shortened:
	*/
	k1 = kstar + DTERMS;

	if ( (k0 < freqIndex0) || (k1 > freqIndex1) ) 
	  {
	    LogPrintf(LOG_CRITICAL,
		      "Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n"
		      "\t\t[Parameters: alpha:%d, Dphi_alpha:%e, Tsft:%e, *Tdot_al:%e]\n",
		      k0, k1, freqIndex0, freqIndex1,
		      alpha, Dphi_alpha, Tsft, *Tdot_al);
	    XLAL_ERROR("LocalXLALComputeFaFb", XLAL_EDOM);
	  }

      } /* compute kappa_star, lambda_alpha */

      /* ---------- calculate the (truncated to DTERMS) sum over k ---------- */

      /* ---------- ATTENTION: this the "hot-loop", which will be 
       * executed many millions of times, so anything in here 
       * has a HUGE impact on the whole performance of the code.
       * 
       * DON'T touch *anything* in here unless you really know 
       * what you're doing !!
       *------------------------------------------------------------
       */

      {
	COMPLEX8 *Xalpha_l = Xalpha + k0 - freqIndex0;  /* first frequency-bin in sum */

	/* if no danger of denominator -> 0 */
#ifdef __GNUC__
	/** somehow the branch prediction of gcc-4.1.2 terribly failes
	    with the current case distinction in the hot-loop,
	    having a severe impact on runtime of the E@H Linux App.
	    So let's allow to give gcc a hint which path has a higher probablility */
	if (__builtin_expect((kappa_star > LD_SMALL4) && (kappa_star < 1.0 - LD_SMALL4), (0==0)))
#else
	if ((kappa_star > LD_SMALL4) && (kappa_star < 1.0 - LD_SMALL4)
#endif
	  {
	  /* WARNING: all current optimized loops rely on current implementation of COMPLEX8 type */
#if (EAH_HOTLOOP_VARIANT == EAH_HOTLOOP_VARIANT_SSE_AKOS08)
#include "EinsteinAtHome/hotloop_precalc.ci"
#elif (EAH_HOTLOOP_VARIANT == EAH_HOTLOOP_VARIANT_SSE)
#include "EinsteinAtHome/hotloop_sse.ci"
#elif (EAH_HOTLOOP_VARIANT == EAH_HOTLOOP_VARIANT_ALTIVEC)
#include "EinsteinAtHome/hotloop_altivec.ci"
#elif (EAH_HOTLOOP_VARIANT == EAH_HOTLOOP_VARIANT_AUTOVECT)
#include "EinsteinAtHome/hotloop_autovect.ci"
#else /* EAH_HOTLOOP_VARIANT */
#include "EinsteinAtHome/hotloop_generic.ci"
#endif /* EAH_HOTLOOP_VARIANT */
	 
	  } /* if |remainder| > LD_SMALL4 */
	else
	  { /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
	    UINT4 ind0;
	    if ( kappa_star <= LD_SMALL4 )
	      ind0 = DTERMS - 1;
	    else
	      ind0 = DTERMS;
	    realXP = TWOPI_FLOAT * Xalpha_l[ind0].re;
	    imagXP = TWOPI_FLOAT * Xalpha_l[ind0].im;
	    REAL8 _lambda_alpha = -lambda_alpha;
	    SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
	    SINCOS_2PI_TRIMMED( &imagQ, &realQ, _lambda_alpha );
	  } /* if |remainder| <= LD_SMALL4 */
      }

      /* real- and imaginary part of e^{-i 2 pi lambda_alpha } */
      
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

} /* LocalXLALComputeFaFb() */
