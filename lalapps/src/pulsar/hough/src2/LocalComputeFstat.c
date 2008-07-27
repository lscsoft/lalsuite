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
#if defined(__INTEL_COMPILER) ||  defined(_MSC_VER)
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



/** the way of trimming x to the interval [0..1) for the sin_cos_LUT functions
    give significant differences in speed, so we provide various ways here.
    We also record the way we are using for logging */

#if EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_PLUS2
/* this only makes sense for the linear sin/cos approximation */
#ifdef _MSC_VER /* no C99 */
#define SINCOS_TRIM_X(y,x) \
  { \
    __asm FLD     QWORD PTR x 	\
    __asm FRNDINT             	\
    __asm FSUBR   QWORD PTR x 	\
    __asm FLD1                	\
    __asm FADDP   ST(1),ST	\
    __asm FSTP    QWORD PTR y 	\
    }
#else
#define SINCOS_TRIM_X(y,x) \
  y = x - rint(x) + 1.0;
#endif
#elif EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_FLOOR 
#define SINCOS_TRIM_X(y,x) \
  y = x - floor(x);
#elif EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_INT4
#define SINCOS_TRIM_X(y,x) \
  y = x-(INT4)x; \
  if ( y < 0.0 ) { y += 1.0; }
#elif EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_INT8
#define SINCOS_TRIM_X(y,x) \
  y = x-(INT8)x; \
  if ( y < 0.0 ) { y += 1.0; }
#elif EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_MODF
#define SINCOS_TRIM_X(y,x) \
{ \
  REAL8 dummy; \
  y = modf(x, &dummy); \
  if ( y < 0.0 ) { y += 1.0; } \
}
#endif

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
#define NUM_FACT 6
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0) };

/* empty initializers  */
static const LALStatus empty_status;

/* sin/cos Lookup tables */
#if (EAH_SINCOS_VARIANT == EAH_SINCOS_VARIANT_LAL)
#define SINCOS_LUT_RES 64 /* resolution of lookup-table */
static REAL4 sinLUT[SINCOS_LUT_RES+1];
static REAL4 cosLUT[SINCOS_LUT_RES+1];
#elif (EAH_SINCOS_VARIANT == EAH_SINCOS_VARIANT_LINEAR)
#define SINCOS_LUT_RES 1024 /* should be multiple of 4 */
static REAL4 sincosLUTbase[SINCOS_LUT_RES+SINCOS_LUT_RES/4];
static REAL4 sincosLUTdiff[SINCOS_LUT_RES+SINCOS_LUT_RES/4];
#endif

/*---------- internal prototypes ----------*/
extern int finite(double x);

static void
LocalComputeFStat ( LALStatus*, Fcomponents*, const PulsarDopplerParams*,
		    const MultiSFTVector*, const MultiNoiseWeights*,
		    const MultiDetectorStateSeries*, const ComputeFParams*,
		    ComputeFBuffer*);

static int
LocalXLALComputeFaFb (Fcomponents*, const SFTVector*, const PulsarSpins,
		      const SSBtimes*, const AMCoeffs*, const ComputeFParams*);

static int local_sin_cos_2PI_LUT_trimmed (REAL4 *sinx, REAL4 *cosx, REAL8 x); 
static void local_sin_cos_2PI_LUT_init (void);

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

      /* make sure Dterms is waht we expect */
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

  /* ----- check validity of input */
#ifndef LAL_NDEBUG
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
      COMPLEX8 *Xalpha_l; 	/* pointer to frequency-bin k in current SFT */
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
	kappa_max = kappa_star + Dterms_1f;

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

      Xalpha_l = Xalpha + k0 - freqIndex0;  /* first frequency-bin in sum */
/* not needed
      realXP = 0;
      imagXP = 0;
*/

      /* if no danger of denominator -> 0 */
      if (__builtin_expect((kappa_star > LD_SMALL4) && (kappa_star < 1.0 - LD_SMALL4), (0==0)))

	/* -------- NOTE: to understand the variants, read them in the order:
	 * - generic  (bottom-most after #else)
	 * - AUTOVECT (still plain C)
	 * - ALTIVEC  (uses vector intrinsics)
	 * - SSE      (Assembler)
	 * (in the code roughly from bottom to top)
	*/

	/* WARNING: all current optimized loops rely on current implementation of COMPLEX8 type */
#if (EAH_HOTLOOP_VARIANT == EAH_HOTLOOP_VARIANT_SSE_AKOS08)

	/** SSE version from Akos */

#ifdef __GNUC__
        {

          COMPLEX8 XSums __attribute__ ((aligned (16))); /* sums of Xa.re and Xa.im for SSE */
	  REAL4 kappa_s = kappa_star; /* single precision version of kappa_star */

	  static REAL4 *scd 		  =  &(sincosLUTdiff[0]);
	  static REAL4 *scb 		  =  &(sincosLUTbase[0]);
	  static REAL4 M1 = -1.0f;
	  static REAL8 sincos_adds = 402653184.0;
	  REAL8 tmp;
          REAL8 _lambda_alpha = -lambda_alpha;
          /* vector constants */
          /* having these not aligned will crash the assembler code */
          static  REAL4 D2222[4] __attribute__ ((aligned (16))) = { 2,2,2,2 };
	  
	  static  REAL4 D1100[4] __attribute__ ((aligned (16))) = {1.0f, 1.0f, 0.0f, 0.0f};
	  static  REAL4 D3322[4] __attribute__ ((aligned (16))) = {3.0f, 3.0f, 2.0f, 2.0f};
	  static  REAL4 D5544[4] __attribute__ ((aligned (16))) = {5.0f, 5.0f, 4.0f, 4.0f};
	  static  REAL4 D7766[4] __attribute__ ((aligned (16))) = {7.0f, 7.0f, 6.0f, 6.0f};

	  static  REAL4 Daabb[4] __attribute__ ((aligned (16))) = {-1.0f, -1.0f, -2.0f, -2.0f};
	  static  REAL4 Dccdd[4] __attribute__ ((aligned (16))) = {-3.0f, -3.0f, -4.0f, -4.0f};
	  static  REAL4 Deeff[4] __attribute__ ((aligned (16))) = {-5.0f, -5.0f, -6.0f, -6.0f};
	  static  REAL4 Dgghh[4] __attribute__ ((aligned (16))) = {-7.0f, -7.0f, -8.0f, -8.0f};

	  /* hand-coded SSE version from Akos */

	  /* one loop iteration as a macro */

#ifdef EAH_HOTLOOP_INTERLEAVED
/* Macros to interleave linear sin/cos calculation (in x87 opcodes)
   with SSE hotloop.*/

/* Version 1 : with trimming of input argument 
   to [0,2) */ 
#define LIN_SIN_COS_TRIM_P0A(alpha) \
		"fldl %[" #alpha "] \n\t"   /* st: alpha */ \
		"fistpll %[tmp] \n\t"	    /* tmp=(INT8)(round((alpha)) */ \
		"fld1 \n\t" 	            /* st: 1.0 */ \
		"fildll %[tmp] \n\t"        /* st: 1.0;(round((alpha))*/ 

#define LIN_SIN_COS_TRIM_P0B(alpha)\
		"fsubrp %%st,%%st(1) \n\t"  /* st: 1.0 -round(alpha) */  \
		"faddl %[" #alpha "] \n\t"  /* st: alpha -round(alpha)+1.0*/ \
		"faddl  %[sincos_adds]  \n\t" /* ..continue lin. sin/cos as lebow */ \
		"fstpl  %[tmp]    \n\t" 
/* Version 2 : assumes input argument is already trimmed */ 
		
#define LIN_SIN_COS_P0(alpha) \
		"fldl %[" #alpha "] \n\t"     /*st:alpha */\
		"faddl  %[sincos_adds]  \n\t" /*st:alpha+A */\
		"fstpl  %[tmp]    \n\t"
#define LIN_SIN_COS_P1 \
		"mov  %[tmp],%%eax \n\t"      /* alpha +A ->eax (ix)*/ \
                "mov  %%eax,%%edx  \n\t"      /* n  = ix & SINCOS_MASK2 */\
		"and  $0x3fff,%%eax \n\t"     	
#define LIN_SIN_COS_P2 \
		"mov  %%eax,%[tmp] \n\t"     \
		"mov  %[scd], %%eax \n\t"    \
		"and  $0xffffff,%%edx \n\t"   /* i  = ix & SINCOS_MASK1;*/
#define LIN_SIN_COS_P3 \
		"fildl %[tmp]\n\t" \
		"sar $0xe,%%edx \n\t"        /*  i  = i >> SINCOS_SHIFT;*/\
		"fld %%st  \n\t"   	     /* st: n; n; */
#define LIN_SIN_COS_P4 \
		"fmuls (%%eax,%%edx,4)   \n\t" \
		"mov  %[scb], %%edi \n\t" \
		"fadds (%%edi,%%edx,4)   \n\t" /*st:sincosLUTbase[i]+n*sincosLUTdiff[i]; n*/
#define LIN_SIN_COS_P5(sin)\
		"add $0x100,%%edx \n\t"   /*edx+=SINCOS_LUT_RES/4*/\
		"fstps %[" #sin "] \n\t"  /*(*sin)=sincosLUTbase[i]+n*sincosLUTdiff[i]*/\
		"fmuls (%%eax,%%edx,4)   \n\t"
#define LIN_SIN_COS_P6(cos) \
		"fadds (%%edi,%%edx,4)   \n\t" \
		"fstps %[" #cos "] \n\t" /*(*cos)=cosbase[i]+n*cosdiff[i];*/

	
#else
#define LIN_SIN_COS_TRIM_P0A(alpha) ""
#define LIN_SIN_COS_TRIM_P0B(alpha) ""
#define LIN_SIN_COS_P0(alpha) ""
#define LIN_SIN_COS_P1 ""
#define LIN_SIN_COS_P2 "" 
#define LIN_SIN_COS_P3 ""
#define LIN_SIN_COS_P4(sin) ""
#define LIN_SIN_COS_P5(cos) ""		

#ifndef LAL_NDEBUG
	  if ( local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star ) ) {
	    XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EFUNC);
	  }
#else
	  local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star );
#endif	  

          SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
#endif /* EAH_HOTLOOP_INTERLEAVED */
          __asm __volatile
	    (
		"movaps %[D7766],%%xmm0\n\t"
		"movaps %[D5544],%%xmm1	\n\t"
		"movups (%[Xa]),%%xmm2	\n\t"
		"movups 0x10(%[Xa]),%%xmm3 \n\t"
		"movss  %[kappa_s],%%xmm7\n\t"
		"shufps $0x0,%%xmm7,%%xmm7\n\t"
	LIN_SIN_COS_P0(kappa_star)
		"addps  %%xmm7,%%xmm0\n\t"
		"addps  %%xmm7,%%xmm1\n\t"
		"rcpps  %%xmm0,%%xmm0\n\t"
		"rcpps  %%xmm1,%%xmm1\n\t"
		"mulps  %%xmm2,%%xmm0\n\t"
		"mulps  %%xmm3,%%xmm1\n\t"
	LIN_SIN_COS_P1
		"addps  %%xmm1,%%xmm0\n\t"
		"movaps %[D3322],%%xmm2\n\t"
		"movaps %[Dccdd],%%xmm3\n\t"
		"movups 0x20(%[Xa]),%%xmm4\n\t"
		"movups 0x50(%[Xa]),%%xmm5\n\t"
	LIN_SIN_COS_P2
		"addps  %%xmm7,%%xmm2\n\t"
		"addps  %%xmm7,%%xmm3\n\t"
		"rcpps  %%xmm2,%%xmm2\n\t"
		"rcpps  %%xmm3,%%xmm3\n\t"
		"mulps  %%xmm4,%%xmm2\n\t"
		"mulps  %%xmm5,%%xmm3\n\t"
	LIN_SIN_COS_P3
		"addps  %%xmm3,%%xmm2\n\t"
		"movaps %[Deeff],%%xmm4\n\t"
		"movaps %[Dgghh],%%xmm5\n\t"
		"movups 0x60(%[Xa]),%%xmm1\n\t"
		"movups 0x70(%[Xa]),%%xmm6\n\t"
	LIN_SIN_COS_P4
		"addps  %%xmm7,%%xmm4\n\t"
		"addps  %%xmm7,%%xmm5\n\t"
		"rcpps  %%xmm4,%%xmm4\n\t"
		"rcpps  %%xmm5,%%xmm5\n\t"
		"mulps  %%xmm1,%%xmm4\n\t"
		"mulps  %%xmm6,%%xmm5\n\t"
	LIN_SIN_COS_P5(sin)
		"addps  %%xmm2,%%xmm0\n\t"
		"addps  %%xmm5,%%xmm4\n\t"
		"movaps %[D1100],%%xmm1\n\t"
		"movaps %[Daabb],%%xmm2\n\t"
	LIN_SIN_COS_P6(cos)
		"addps  %%xmm7,%%xmm1\n\t"
		"addps  %%xmm7,%%xmm2\n\t"
		"rcpps  %%xmm1,%%xmm5\n\t"
		"rcpps  %%xmm2,%%xmm6\n\t"
	LIN_SIN_COS_TRIM_P0A(_lambda_alpha)
		"addps  %%xmm4,%%xmm0\n\t"
		"movaps %[D2222],%%xmm3\n\t"
		"movaps %[D2222],%%xmm4\n\t"
		"mulps  %%xmm5,%%xmm1\n\t"
		"mulps  %%xmm6,%%xmm2\n\t"
	LIN_SIN_COS_TRIM_P0B(_lambda_alpha)
		"subps  %%xmm1,%%xmm3\n\t"
		"subps  %%xmm2,%%xmm4\n\t"
		"mulps  %%xmm3,%%xmm5\n\t"
		"mulps  %%xmm4,%%xmm6\n\t"
		"movups 0x30(%[Xa]),%%xmm1\n\t"
		"movups 0x40(%[Xa]),%%xmm2\n\t"
	LIN_SIN_COS_P1
		"mulps  %%xmm5,%%xmm1\n\t"
		"mulps  %%xmm6,%%xmm2\n\t"
		"addps  %%xmm1,%%xmm0\n\t"
		"addps  %%xmm2,%%xmm0\n\t"
	LIN_SIN_COS_P2
		"movhlps %%xmm0,%%xmm1\n\t"
		"addps  %%xmm1,%%xmm0\n\t"

/*	  
        c_alpha-=1.0f;
	  realXP = s_alpha * XSums.re - c_alpha * XSums.im;
	  imagXP = c_alpha * XSums.re + s_alpha * XSums.im;
*/

		"movss %[M1],%%xmm5 \n\t"
		"movaps %%xmm0,%%xmm3 \n\t"
		"shufps $1,%%xmm3,%%xmm3 \n\t"	
	LIN_SIN_COS_P3
		"movss %[cos],%%xmm2 \n\t"
		"movss %[sin],%%xmm1 \n\t"
		"addss %%xmm5,%%xmm2 \n\t"	
		"movss %%xmm2,%%xmm6 \n\t"	
	LIN_SIN_COS_P4
		"movss %%xmm1,%%xmm5  \n\t"
		"mulss %%xmm0,%%xmm1 \n\t"		
		"mulss %%xmm0,%%xmm2 \n\t"

	LIN_SIN_COS_P5(Qimag)
		"mulss %%xmm3,%%xmm5 \n\t"
		"mulss %%xmm3,%%xmm6 \n\t"
		"addss %%xmm5,%%xmm2 \n\t"
		"subss %%xmm6,%%xmm1 \n\t"
	LIN_SIN_COS_P6(Qreal)
		"MOVSS	%%xmm2,%[XPimag]   	\n\t"	/*  */
		"MOVSS	%%xmm1,%[XPreal]   	\n\t"	/*  */

	     /* interface */
	     :
	     /* output  (here: to memory)*/
	     [XPreal]      "=m" (realXP),
	     [XPimag]      "=m" (imagXP),
	     [Qreal]      "=m" (realQ),
	     [Qimag]      "=m" (imagQ),
	     [sin]	  "=m" (s_alpha),
	     [cos]	  "=m" (c_alpha),
	     [tmp]        "=m" (tmp)

	     :
	     /* input */
	     [Xa]          "r"  (Xalpha_l),
	     [kappa_s]     "m"  (kappa_s),
	     [kappa_star]  "m"  (kappa_star),
	     [_lambda_alpha] "m" (_lambda_alpha),
	     [scd]	   "m"  (scd),
	     [scb]	   "m"  (scb),
	     [sincos_adds] "m"  (sincos_adds),
	     [M1]	  "m" (M1),


	     /* vector constants */
	     [D2222]       "m"  (D2222[0]),
	     [D1100]       "m"  (D1100[0]),
	     [D3322]       "m"  (D3322[0]),
	     [D5544]       "m"  (D5544[0]),
	     [D7766]       "m"  (D7766[0]),
	     [Daabb]       "m"  (Daabb[0]),
	     [Dccdd]       "m"  (Dccdd[0]),
	     [Deeff]       "m"  (Deeff[0]),
	     [Dgghh]       "m"  (Dgghh[0])

	     :
	     /* clobbered registers */
	     "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","st","st(1)","st(2)","eax","edx","edi","cc"
	     );

	  /* moved the sin/cos call down here to avoid the store/forward stall of Core2s */

	  /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
	   * the trig-functions need to be calculated only once!
	   * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the 
	   * closest to zero and will pose no numerical difficulties !
	   */

	}
#else /* __GNUC__ */
	{
	  __declspec(align(16)) static struct { REAL4 a,b,c,d; } v0011 = {0.0, 0.0, 1.0, 1.0};
	  __declspec(align(16)) static struct { REAL4 a,b,c,d; } v2222 = {2.0, 2.0, 2.0, 2.0};
  	  __declspec(align(16)) COMPLEX8 STn; 
	
	  REAL4 kappa_m = kappa_max; /* single precision version of kappa_max */
	      
	  /* prelude */
	  __asm {
 	      mov      esi , Xalpha_l                 /* Xal = Xalpha_l         */
	      movss    xmm2, kappa_m                  /* pn[0] = kappa_max      */
	      movlps   xmm1, MMWORD PTR [esi]         /* STnV = Xal ...         */
	      movhps   xmm1, MMWORD PTR [esi+8]       /* ... continued          */
	      shufps   xmm2, xmm2, 0                  /* pn[3]=pn[2]=pn[1]=pn[0]*/
	      movaps   xmm4, XMMWORD PTR v2222        /* xmm4 = V2222           */
	      subps    xmm2, XMMWORD PTR v0011        /* pn[2]-=1.0; pn[3]-=1.0 */
	      movaps   xmm0, xmm2                     /* qn = pn                */
	      };

	  /* one loop iteration as a macro */
#define VEC_LOOP_AV(a,b)\
	  { \
	      __asm movlps   xmm3, MMWORD PTR [esi+a] /* Xai = Xal[a]  ...*/\
	      __asm movhps   xmm3, MMWORD PTR [esi+b] /* ... continued    */\
	      __asm subps    xmm2, xmm4		      /* pn   -= V2222    */\
	      __asm mulps    xmm3, xmm0		      /* Xai  *= qn       */\
	      __asm mulps    xmm1, xmm2		      /* STnV *= pn       */\
	      __asm mulps    xmm0, xmm2		      /* qn   *= pn       */\
	      __asm addps    xmm1, xmm3		      /* STnV += Xai      */\
	      }

	  /* seven macro calls i.e. loop iterations */
	  VEC_LOOP_AV(16,24);
	  VEC_LOOP_AV(32,40);
	  VEC_LOOP_AV(48,56);
	  VEC_LOOP_AV(64,72);
	  VEC_LOOP_AV(80,88);
	  VEC_LOOP_AV(96,104);
	  VEC_LOOP_AV(112,120);

	  /* four divisions and summing in SSE, then write out the result */
	  __asm {
	      divps    xmm1, xmm0                     /* STnV      /= qn       */
	      movhlps  xmm4, xmm1                     /* / STnV[0] += STnV[2] \ */
	      addps    xmm4, xmm1                     /* \ STnV[1] += STnV[3] / */
	      movlps   STn, xmm4                      /* STn = STnV */
	      };

	  /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
	   * the trig-functions need to be calculated only once!
	   * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the 
	   * closest to zero and will pose no numerical difficulties !
	   */
#ifndef LAL_NDEBUG
	  if ( local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star ) ) {
	    XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EFUNC);
	  }
#else
	  local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star );
#endif
	  c_alpha -= 1.0f;

	  realXP = s_alpha * STn.re - c_alpha * STn.im;
	  imagXP = c_alpha * STn.re + s_alpha * STn.im;
	     
	}
#endif /* __GNUC__ */




#elif (EAH_HOTLOOP_VARIANT == EAH_HOTLOOP_VARIANT_SSE)

	/** SSE version from Akos */

#ifdef __GNUC__
        {

          COMPLEX8 XSums __attribute__ ((aligned (16))); /* sums of Xa.re and Xa.im for SSE */
	  REAL4 kappa_m = kappa_max; /* single precision version of kappa_max */

	  static REAL4 *scd 		  =  &(sincosLUTdiff[0]);
	  static REAL4 *scb 		  =  &(sincosLUTbase[0]);
	  static REAL4 M1 = -1.0f;
	  static REAL8 sincos_adds = 402653184.0;
	  REAL8 tmp;
          REAL8 _lambda_alpha = -lambda_alpha;
          /* vector constants */
          /* having these not aligned will crash the assembler code */
          static REAL4 V0011[4] __attribute__ ((aligned (16))) = { 0,0,1,1 };
          static REAL4 V2222[4] __attribute__ ((aligned (16))) = { 2,2,2,2 };
	  
	  /* hand-coded SSE version from Akos */

	  /* one loop iteration as a macro */

#define VEC_LOOP_AV(a)\
	     "MOVUPS " #a "(%[Xa]),%%xmm3   	\n\t" \
	     "SUBPS	%%xmm4,%%xmm2   	\n\t" \
	     "MULPS	%%xmm0,%%xmm3   	\n\t" \
	     "MULPS	%%xmm2,%%xmm1   	\n\t" \
	     "MULPS	%%xmm2,%%xmm0   	\n\t" \
	     "ADDPS	%%xmm3,%%xmm1   	\n\t" 


#ifdef EAH_HOTLOOP_INTERLEAVED

#define LIN_SIN_COS_TRIM_P0A(alpha) \
		"fldl %[" #alpha "] \n\t" \
		"fistpll %[tmp] \n\t" \
		"fld1 \n\t" 	\
		"fildll %[tmp] \n\t" 

#define LIN_SIN_COS_TRIM_P0B(alpha)\
		"fsubrp %%st,%%st(1) \n\t" \
		"faddl %[" #alpha "] \n\t" \
		"faddl  %[sincos_adds]  \n\t" \
		"fstpl  %[tmp]    \n\t" 
		
#define LIN_SIN_COS_P0(alpha) \
		"fldl %[" #alpha "] \n\t" \
		"faddl  %[sincos_adds]  \n\t" \
		"fstpl  %[tmp]    \n\t" 

#define LIN_SIN_COS_P1 \
		"mov  %[tmp],%%eax \n\t"  \
                "mov  %%eax,%%edx  \n\t" \
		"and  $0x3fff,%%eax \n\t" 
#define LIN_SIN_COS_P2 \
		"mov  %%eax,%[tmp] \n\t"   \
		"mov  %[scd], %%eax \n\t"\
		"and  $0xffffff,%%edx \n\t" \
		"fildl %[tmp]\n\t" 
#define LIN_SIN_COS_P3 \
		"sar $0xe,%%edx \n\t" \
		"fld %%st  \n\t"   \
		"fmuls (%%eax,%%edx,4)   \n\t" \
		"mov  %[scb], %%edi \n\t"
#define LIN_SIN_COS_P4(sin)\
		"fadds (%%edi,%%edx,4)   \n\t" \
		"add $0x100,%%edx \n\t"   \
		"fstps %[" #sin "] \n\t" \
		"fmuls (%%eax,%%edx,4)   \n\t"
#define LIN_SIN_COS_P5(cos) \
		"fadds (%%edi,%%edx,4)   \n\t" \
		"fstps %[" #cos "] \n\t"		

	
#else
#define LIN_SIN_COS_TRIM_P0A(alpha) ""
#define LIN_SIN_COS_TRIM_P0B(alpha) ""
#define LIN_SIN_COS_P0(alpha) ""
#define LIN_SIN_COS_P1 ""
#define LIN_SIN_COS_P2 "" 
#define LIN_SIN_COS_P3 ""
#define LIN_SIN_COS_P4(sin) ""
#define LIN_SIN_COS_P5(cos) ""		

#ifndef LAL_NDEBUG
	  if ( local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star ) ) {
	    XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EFUNC);
	  }
#else
	  local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star );
#endif	  

          SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
#endif
          __asm __volatile
	    (
	     /* -------------------------------------------------------------------; */
	     /* Prepare common divisor method for 4 values ( two Re-Im pair ) */
	     /*  Im1, Re1, Im0, Re0 */
	     "MOVAPS	%[V0011],%%xmm6   	\n\t"	
	     "MOVSS	%[kappa_m],%%xmm2   	\n\t"	/* X2:  -   -   -   C */
	     "MOVUPS	(%[Xa]),%%xmm1   	\n\t"	/* X1: Y01 X01 Y00 X00 */
	     "SHUFPS	$0,%%xmm2,%%xmm2   	\n\t"	/* X2:  C   C   C   C */
	     "MOVAPS	%[V2222],%%xmm4   	\n\t"	/* X7:  2   2   2   2 */
	     "SUBPS	%%xmm6,%%xmm2   	\n\t"	/* X2: C-1 C-1  C   C */
	     /* -------------------------------------------------------------------; */
	     "MOVAPS	%%xmm2,%%xmm0   	\n\t"	/* X0: C-1 C-1  C   C */
	     /* -------------------------------------------------------------------; */
	     /* xmm0: collected denumerators -> a new element will multiply by this */
	     /* xmm1: collected numerators -> we will divide it by the denumerator last */
	     /* xmm2: current denumerator ( counter type ) */
	     /* xmm3: current numerator ( current Re,Im elements ) */
	     /* -------------------------------------------------------------------; */

	     /* seven "loop iterations" (unrolled) */
LIN_SIN_COS_P0(kappa_star)
	     VEC_LOOP_AV(16)
LIN_SIN_COS_P1
	     VEC_LOOP_AV(32)
LIN_SIN_COS_P2
	     VEC_LOOP_AV(48)
LIN_SIN_COS_P3
	     VEC_LOOP_AV(64)
LIN_SIN_COS_P4(sin)
	     VEC_LOOP_AV(80)
LIN_SIN_COS_P5(cos)
	     VEC_LOOP_AV(96)

LIN_SIN_COS_TRIM_P0A(_lambda_alpha)


#if EAH_HOTLOOP_DIVS == EAH_HOTLOOP_DIVS_RECIPROCAL
	     "movhlps %%xmm6,%%xmm6     \n\t"
#endif
	     "movss %[M1] , %%xmm5 \n\t"
	     VEC_LOOP_AV(112)


#if EAH_HOTLOOP_DIVS == EAH_HOTLOOP_DIVS_RECIPROCAL
	     "movhlps %%xmm0,%%xmm2   	\n\t"
	     "mulss %%xmm0,%%xmm2	\n\t"

#ifdef EAH_HOTLOOP_RENR 
	     "RCPSS %%xmm2,%%xmm6	\n\t"
	     "MULSS %%xmm6,%%xmm2	\n\t"
	     "MULSS %%xmm6,%%xmm2	\n\t"
	LIN_SIN_COS_TRIM_P0B(_lambda_alpha)
	     "ADDSS %%xmm6,%%xmm6	\n\t"
	     "SUBSS %%xmm2,%%xmm6	\n\t"
#else
	     "divss %%xmm2,%%xmm6	\n\t"
	LIN_SIN_COS_TRIM_P0B(_lambda_alpha)

#endif

	     "shufps $78,%%xmm0,%%xmm0  \n\t"
	     "mulps  %%xmm1,%%xmm0      \n\t"
	     "movhlps %%xmm0,%%xmm4	\n\t"
	     "addps %%xmm0,%%xmm4	\n\t"

	     "shufps $160,%%xmm6,%%xmm6   \n\t"
	     "mulps %%xmm6,%%xmm4	\n\t"		
#else
	     /* -------------------------------------------------------------------; */
	     /* Four divisions at once ( two for real parts and two for imaginary parts ) */
#ifdef EAH_HOTLOOP_RENR 
	     "RCPPS %%xmm0,%%xmm6	\n\t"
	     "MULPS %%xmm6,%%xmm0	\n\t"
	     "MULPS %%xmm6,%%xmm0	\n\t"
	LIN_SIN_COS_TRIM_P0B(_lambda_alpha)
	     "ADDPS %%xmm6,%%xmm6	\n\t"
	     "SUBPS %%xmm0,%%xmm6	\n\t"
	     "MULPS %%xmm6,%%xmm1	\n\t"
#else
	     "DIVPS	%%xmm0,%%xmm1   	\n\t"	/* X1: Y0G X0G Y1F X1F */
	LIN_SIN_COS_TRIM_P0B(_lambda_alpha)
#endif
	     /* -------------------------------------------------------------------; */
	     /* So we have to add the two real and two imaginary parts */
	     "MOVHLPS   %%xmm1,%%xmm4	        \n\t"	/* X4:  -   -  Y0G X0G */
	     "ADDPS	%%xmm1,%%xmm4   	\n\t"	/* X4:  -   -  YOK XOK */
	     /* -------------------------------------------------------------------; */
#endif

/*	  
        c_alpha-=1.0f;
	  realXP = s_alpha * XSums.re - c_alpha * XSums.im;
	  imagXP = c_alpha * XSums.re + s_alpha * XSums.im;
*/

		"movaps %%xmm4,%%xmm3 \n\t"
		"shufps $1,%%xmm3,%%xmm3 \n\t"
		"movss %[cos],%%xmm2 \n\t"
	LIN_SIN_COS_P1	
		"movss %[sin],%%xmm1 \n\t"
		"addss %%xmm5,%%xmm2\n\t"	
		"movss %%xmm2,%%xmm6 \n\t"
	LIN_SIN_COS_P2	
		"movss %%xmm1,%%xmm5  \n\t"
		"mulss %%xmm4,%%xmm1 \n\t"		
		"mulss %%xmm4,%%xmm2 \n\t"
	LIN_SIN_COS_P3
		"mulss %%xmm3,%%xmm5 \n\t"
		"mulss %%xmm3,%%xmm6 \n\t"
	LIN_SIN_COS_P4(Qimag)		
		"addss %%xmm5,%%xmm2 \n\t"
		"subss %%xmm6,%%xmm1 \n\t"
	LIN_SIN_COS_P5(Qreal)
		"MOVss	%%xmm2,%[XPimag]   	\n\t"	/*  */
		"MOVss	%%xmm1,%[XPreal]   	\n\t"	/*  */

	     /* interface */
	     :
	     /* output  (here: to memory)*/
	     [XPreal]      "=m" (realXP),
	     [XPimag]      "=m" (imagXP),
	     [Qreal]      "=m" (realQ),
	     [Qimag]      "=m" (imagQ),
	     [tmp]        "=m" (tmp),
	     [sin]	  "=m" (s_alpha),
	     [cos]	  "=m" (c_alpha)

	     :
	     /* input */
	     [Xa]          "r"  (Xalpha_l),
	     [kappa_m]     "m"  (kappa_m),
	     [kappa_star]  "m"  (kappa_star),
	     [_lambda_alpha] "m" (_lambda_alpha),
	     [scd]	   "m"  (scd),
	     [scb]	   "m"  (scb),
	     [sincos_adds]       "m"  (sincos_adds),
	     [M1]	  "m" (M1),


	     /* vector constants */
	     [V0011]       "m"  (V0011[0]),
	     [V2222]       "m"  (V2222[0])

#ifndef IGNORE_XMM_REGISTERS
	     :
	     /* clobbered registers */
	     "xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","st","st(1)","st(2)","eax","edx","edi","cc"
#endif
	     );

	  /* moved the sin/cos call down here to avoid the store/forward stall of Core2s */

	  /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
	   * the trig-functions need to be calculated only once!
	   * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the 
	   * closest to zero and will pose no numerical difficulties !
	   */

	}
#else /* __GNUC__ */
	{
	  __declspec(align(16)) static struct { REAL4 a,b,c,d; } v0011 = {0.0, 0.0, 1.0, 1.0};
	  __declspec(align(16)) static struct { REAL4 a,b,c,d; } v2222 = {2.0, 2.0, 2.0, 2.0};
  	  __declspec(align(16)) COMPLEX8 STn; 
	
	  REAL4 kappa_m = kappa_max; /* single precision version of kappa_max */
	      
	  /* prelude */
	  __asm {
 	      mov      esi , Xalpha_l                 /* Xal = Xalpha_l         */
	      movss    xmm2, kappa_m                  /* pn[0] = kappa_max      */
	      movlps   xmm1, MMWORD PTR [esi]         /* STnV = Xal ...         */
	      movhps   xmm1, MMWORD PTR [esi+8]       /* ... continued          */
	      shufps   xmm2, xmm2, 0                  /* pn[3]=pn[2]=pn[1]=pn[0]*/
	      movaps   xmm4, XMMWORD PTR v2222        /* xmm4 = V2222           */
	      subps    xmm2, XMMWORD PTR v0011        /* pn[2]-=1.0; pn[3]-=1.0 */
	      movaps   xmm0, xmm2                     /* qn = pn                */
	      };

	  /* one loop iteration as a macro */
#define VEC_LOOP_AV(a,b)\
	  { \
	      __asm movlps   xmm3, MMWORD PTR [esi+a] /* Xai = Xal[a]  ...*/\
	      __asm movhps   xmm3, MMWORD PTR [esi+b] /* ... continued    */\
	      __asm subps    xmm2, xmm4		      /* pn   -= V2222    */\
	      __asm mulps    xmm3, xmm0		      /* Xai  *= qn       */\
	      __asm mulps    xmm1, xmm2		      /* STnV *= pn       */\
	      __asm mulps    xmm0, xmm2		      /* qn   *= pn       */\
	      __asm addps    xmm1, xmm3		      /* STnV += Xai      */\
	      }

	  /* seven macro calls i.e. loop iterations */
	  VEC_LOOP_AV(16,24);
	  VEC_LOOP_AV(32,40);
	  VEC_LOOP_AV(48,56);
	  VEC_LOOP_AV(64,72);
	  VEC_LOOP_AV(80,88);
	  VEC_LOOP_AV(96,104);
	  VEC_LOOP_AV(112,120);

	  /* four divisions and summing in SSE, then write out the result */
	  __asm {
	      divps    xmm1, xmm0                     /* STnV      /= qn       */
	      movhlps  xmm4, xmm1                     /* / STnV[0] += STnV[2] \ */
	      addps    xmm4, xmm1                     /* \ STnV[1] += STnV[3] / */
	      movlps   STn, xmm4                      /* STn = STnV */
	      };

	  /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
	   * the trig-functions need to be calculated only once!
	   * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the 
	   * closest to zero and will pose no numerical difficulties !
	   */
#ifndef LAL_NDEBUG
	  if ( local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star ) ) {
	    XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EFUNC);
	  }
#else
	  local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star );
#endif
	  c_alpha -= 1.0f;

	  realXP = s_alpha * STn.re - c_alpha * STn.im;
	  imagXP = c_alpha * STn.re + s_alpha * STn.im;
	     
	}
#endif /* __GNUC__ */

#elif (EAH_HOTLOOP_VARIANT == EAH_HOTLOOP_VARIANT_ALTIVEC)

	{
#ifndef LAL_NDEBUG
	  if ( local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star ) ) {
	  XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EFUNC);
	  }
#else
	  local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star );
#endif
	  c_alpha -= 1.0f;
	  {
	    REAL4 *Xalpha_kR4 = (REAL4*)(Xalpha_l);
	    
	    float STn[4] __attribute__ ((aligned (16))); /* aligned for vector output */
	    /* the vectors actually become registers in the AVUnit */
	    vector unsigned char perm;      /* permutation pattern for unaligned memory access */
	    vector float load0, load1, load2;    /* temp registers for unaligned memory access */
	    vector float XaiV   /* xmm3 */;                     /* SFT data loaded from memory */
	    vector float STnV   /* xmm1 */;                            /* sums up the dividend */
	    vector float V0000             = {0,0,0,0};                /* zero vector constant */
	    vector float V2222  /* xmm4 */ = {2,2,2,2};                     /* vector constant */
	    vector float pnV    /* xmm2 */ = {((float)(kappa_max)),
					      ((float)(kappa_max)),
					      ((float)(kappa_max - 1)),
					      ((float)(kappa_max - 1)) };
	    vector float qnV    /* xmm0 */ = pnV;  /* common divisor, initally = 1.0 * pnV */
	    /*    this column above (^) lists the corresponding register in the SSE version */
	    
	    vector float tV;  /* temporary vector used for Newton-Rhapson iterarion */

	    /* init the memory access (load0,load1) */
	    load0   = vec_ld  (0,(Xalpha_kR4));
	    perm    = vec_lvsl(0,(Xalpha_kR4));
	    load1   = vec_ld  (0,(Xalpha_kR4+4));

	    /* first "iteration" & initialization */
	    XaiV    = vec_perm(load0,load1,perm);
	    qnV     = vec_re(pnV);
	    STnV    = vec_madd(XaiV, qnV, V0000);

	    /* use a reciprocal estimate as a replacement for a division.
	       in our case this is only valid for the "outer" elements of the kernel loop */
#define VEC_LOOP_RE(n,a,b)\
	    pnV     = vec_sub(pnV,V2222);\
	    perm    = vec_lvsl(0,(Xalpha_kR4+(n)));\
            load##b = vec_ld(0,(Xalpha_kR4+(n)+4));\
	    XaiV    = vec_perm(load##a,load##b,perm);\
	    qnV     = vec_re(pnV);\
	    STnV    = vec_madd(XaiV, qnV, STnV);  /* STnV = XaiV * qnV + STnV */

	    /* refine the reciprocal estimate to by a Newton-Rhapson iteration.
	       re1(x) = re0(x) * (2 - x * re0(x))
	       (see http://en.wikipedia.org/wiki/Division_(digital)#Newton-Raphson_division)
	       this should give as much precision as a normal float division */
#define VEC_LOOP_RE_NR(n,a,b)\
	    pnV     = vec_sub(pnV,V2222);\
	    perm    = vec_lvsl(0,(Xalpha_kR4+(n)));\
            load##b = vec_ld(0,(Xalpha_kR4+(n)+4));\
	    XaiV    = vec_perm(load##a,load##b,perm);\
	    qnV     = vec_re(pnV);\
            tV      = vec_madd(qnV,pnV,V0000);\
            tV      = vec_sub(V2222,tV);\
            qnV     = vec_madd(qnV,tV,V0000);\
	    STnV    = vec_madd(XaiV, qnV, STnV);

	    VEC_LOOP_RE(4,1,2);
	    VEC_LOOP_RE_NR(8,2,0);
	    VEC_LOOP_RE_NR(12,0,1);
	    VEC_LOOP_RE_NR(16,1,2);
	    VEC_LOOP_RE_NR(20,2,0);
	    VEC_LOOP_RE(24,0,1);
	    VEC_LOOP_RE(28,1,0);

	    /* output the vector */
	    vec_st(STnV,0,STn);

	    /* combine the sums */
	    {
	      REAL4 U_alpha = STn[0] + STn[2];
	      REAL4 V_alpha = STn[1] + STn[3];

	      realXP = s_alpha * U_alpha - c_alpha * V_alpha;
	      imagXP = c_alpha * U_alpha + s_alpha * V_alpha;
	    }
	  }
	}

#elif (EAH_HOTLOOP_VARIANT == EAH_HOTLOOP_VARIANT_AUTOVECT)

	/* designed for four vector elemens (ve) as there are in SSE and AltiVec */
	/* vectorizes with gcc-4.2.3 and gcc-4.1.3 */

	{
	  /* the initialization already handles the first elements,
	     thus there are only 7 loop iterations left */
	  UINT4 l;
	  UINT4 ve;
	  REAL4 *Xal   = (REAL4*)Xalpha_l;
	  REAL4 STn[4] = {Xal[0],Xal[1],Xal[2],Xal[3]};
	  REAL4 pn[4]  = {kappa_max, kappa_max, kappa_max-1.0f, kappa_max-1.0f};
	  REAL4 qn[4];
	  
	  for ( ve = 0; ve < 4; ve++)
	    qn[ve] = pn[ve];
	  
	  for ( l = 1; l < DTERMS; l ++ ) {
	    Xal += 4;
	    for ( ve = 0; ve < 4; ve++) {
	      pn[ve] -= 2.0f;
	      STn[ve] = pn[ve] * STn[ve] + qn[ve] * Xal[ve];
	      qn[ve] *= pn[ve];
	    }
	  }
	  
	  {
#ifndef LAL_NDEBUG
	    if ( local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star ) ) {
	      XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EFUNC);
	    }
#else
	    local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star );
#endif
	    c_alpha -= 1.0f;
	  }

	  /* combine the partial sums: */
	  {
#if EAH_HOTLOOP_DIVS == EAH_HOTLOOP_DIVS_RECIPROCAL
	    /* if the division is to be done outside the SIMD unit */
	    
	    REAL4 r_qn  = 1.0 / (qn[0] * qn[2]);
	    REAL4 U_alpha = (STn[0] * qn[2] + STn[2] * qn[0]) * r_qn;
	    REAL4 V_alpha = (STn[1] * qn[3] + STn[3] * qn[1]) * r_qn;
	    
	    realXP = s_alpha * U_alpha - c_alpha * V_alpha;
	    imagXP = c_alpha * U_alpha + s_alpha * V_alpha;
	    
#else /* EAH_HOTLOOP_DIVS */
	    /* if the division can and should be done inside the SIMD unit */
	    
	    REAL4 U_alpha, V_alpha;
	    
	    for ( ve = 0; ve < 4; ve++)
	      STn[ve] /= qn[ve];
	    
	    U_alpha = (STn[0] + STn[2]);
	    V_alpha = (STn[1] + STn[3]);
	    
	    realXP = s_alpha * U_alpha - c_alpha * V_alpha;
	    imagXP = c_alpha * U_alpha + s_alpha * V_alpha;

#endif /* EAH_HOTLOOP_DIVS */	    
	  }
	}

#else /* EAH_HOTLOOP_VARIANT */

	/* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
	 * the trig-functions need to be calculated only once!
	 * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the 
	 * closest to zero and will pose no numerical difficulties !
	 */
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

	  /* 2*DTERMS iterations */
	  UINT4 l;
	  for ( l = 1; l < 2*DTERMS; l ++ )
	    {
	      Xalpha_l ++;
	      
	      pn = pn - 1.0f;   		  /* p_(n+1) */
	      Sn = pn * Sn + qn * (*Xalpha_l).re; /* S_(n+1) */
	      Tn = pn * Tn + qn * (*Xalpha_l).im; /* T_(n+1) */
	      qn *= pn; 			  /* q_(n+1) */
	    } /* for l < 2*DTERMS */

#if EAH_HOTLOOP_DIVS == EAH_HOTLOOP_DIVS_RECIPROCAL
	  { /* could hardly be slower than two divisions */
	    REAL4 r_qn = 1.0 / qn;
	    U_alpha = Sn * r_qn;
	    V_alpha = Tn * r_qn;
	  }
#else /* EAH_HOTLOOP_DIVS */
	  U_alpha = Sn / qn;
	  V_alpha = Tn / qn;
#endif /* EAH_HOTLOOP_DIVS */

#ifndef LAL_NDEBUG
	  if ( local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star ) ) {
	    XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EFUNC);
	  }
#else
	  local_sin_cos_2PI_LUT_trimmed ( &s_alpha, &c_alpha, kappa_star );
#endif
	  c_alpha -= 1.0f;
	
 	  realXP = s_alpha * U_alpha - c_alpha * V_alpha;
 	  imagXP = c_alpha * U_alpha + s_alpha * V_alpha;
	}

#endif /* EAH_HOTLOOP_VARIANT */

      /* if |remainder| > LD_SMALL4 */
      else
	{ /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
 	  UINT4 ind0;
   	  if ( kappa_star <= LD_SMALL4 )
	    ind0 = DTERMS - 1;
   	  else
	    ind0 = DTERMS;
 	  realXP = TWOPI_FLOAT * Xalpha_l[ind0].re;
 	  imagXP = TWOPI_FLOAT * Xalpha_l[ind0].im;
#ifdef EAH_HOTLOOP_INTERLEAVED
	REAL8 _lambda_alpha = -lambda_alpha;
	SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
#ifndef LAL_NDEBUG
	if ( local_sin_cos_2PI_LUT_trimmed ( &imagQ, &realQ, _lambda_alpha ) ) {
	  XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EFUNC);
	}
#else
	local_sin_cos_2PI_LUT_trimmed ( &imagQ, &realQ, _lambda_alpha );
#endif
#endif
	} /* if |remainder| <= LD_SMALL4 */

      {
#ifndef EAH_HOTLOOP_INTERLEAVED
	REAL8 _lambda_alpha = -lambda_alpha;
	SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
#ifndef LAL_NDEBUG
	if ( local_sin_cos_2PI_LUT_trimmed ( &imagQ, &realQ, _lambda_alpha ) ) {
	  XLAL_ERROR ( "LocalXLALComputeFaFb", XLAL_EFUNC);
	}
#else
	local_sin_cos_2PI_LUT_trimmed ( &imagQ, &realQ, _lambda_alpha );
#endif
#endif
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



#define SINCOS_LUT_RES_F	(1.0 * SINCOS_LUT_RES)
#define OO_SINCOS_LUT_RES	(1.0 / SINCOS_LUT_RES)

#define X_TO_IND	(1.0 * SINCOS_LUT_RES * OOTWOPI )
#define IND_TO_X	(LAL_TWOPI * OO_SINCOS_LUT_RES)

#if (EAH_SINCOS_VARIANT == EAH_SINCOS_VARIANT_LAL)
/* "traditional" version with 2 LUT */

static void local_sin_cos_2PI_LUT_init (void)
{
  UINT4 k;
  static REAL8 const oo_lut_res = OO_SINCOS_LUT_RES;
  for (k=0; k <= SINCOS_LUT_RES; k++) {
    sinLUT[k] = sin( LAL_TWOPI * k * oo_lut_res );
    cosLUT[k] = cos( LAL_TWOPI * k * oo_lut_res );
  }
}

static int local_sin_cos_2PI_LUT_trimmed (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x)
{
  INT4 i0;
  REAL8 d, d2;
  REAL8 ts, tc;
  static REAL8 const oo_lut_res = OO_SINCOS_LUT_RES;

#ifndef LAL_NDEBUG
  if ( x < 0.0 || x >= 1.0 )
    {
      XLALPrintError("\nFailed numerica in local_sin_cos_2PI_LUT(): x = %f not in [0,1)\n\n", x );
      return XLAL_FAILURE;
    }
#endif

  i0 = (INT4)( x * SINCOS_LUT_RES_F + 0.5 );	/* i0 in [0, SINCOS_LUT_RES ] */
  d = d2 = LAL_TWOPI * (x - oo_lut_res * i0);
  d2 *= 0.5 * d;

  ts = sinLUT[i0];
  tc = cosLUT[i0];
   
  /* use Taylor-expansions for sin/cos around LUT-points */
  (*sin2pix) = ts + d * tc - d2 * ts;
  (*cos2pix) = tc - d * ts - d2 * tc;

  return XLAL_SUCCESS;
} /* local_sin_cos_2PI_LUT() */



#elif (EAH_SINCOS_VARIANT == EAH_SINCOS_VARIANT_LINEAR)

/* linear approximation developed with Akos */

#define SINCOS_ADDS  402653184.0
#define SINCOS_MASK1 0xFFFFFF
#define SINCOS_MASK2 0x003FFF
#define SINCOS_SHIFT 14

static void local_sin_cos_2PI_LUT_init (void)
{
  static const REAL8 step = LAL_TWOPI / (REAL8)SINCOS_LUT_RES;
  static const REAL8 div  = 1.0 / ( 1 << SINCOS_SHIFT );
  REAL8 start, end, true_mid, linear_mid;
  int i;

  start = 0.0; /* sin(0 * step) */
  for( i = 0; i < SINCOS_LUT_RES + SINCOS_LUT_RES/4; i++ ) {
    true_mid = sin( ( i + 0.5 ) * step );
    end = sin( ( i + 1 ) * step );
    linear_mid = ( start + end ) * 0.5;
    sincosLUTbase[i] = start + ( ( true_mid - linear_mid ) * 0.5 );
    sincosLUTdiff[i] = ( end - start ) * div;
    start = end;
  }
}


static int local_sin_cos_2PI_LUT_trimmed (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x) {
  INT4  i, n, ix;
  union {
    REAL8 asreal;
    struct {
      INT4 high;
      INT4 low;
    } as2int;
  } ux;

  static const REAL4* cosbase = sincosLUTbase + (SINCOS_LUT_RES/4);
  static const REAL4* cosdiff = sincosLUTdiff + (SINCOS_LUT_RES/4);

#ifndef LAL_NDEBUG
  if(x > SINCOS_ADDS) {
    LogPrintf(LOG_DEBUG,"sin_cos_LUT: x too large: %22f > %f\n",x,SINCOS_ADDS);
    return XLAL_FAILURE;
  } else if(x < -SINCOS_ADDS) {
    LogPrintf(LOG_DEBUG,"sin_cos_LUT: x too small: %22f < %f\n",x,-SINCOS_ADDS);
    return XLAL_FAILURE;
  }
#endif

#if EAH_SINCOS_F2IBITS == EAH_SINCOS_F2IBITS_MEMCPY
  ux.asreal = x + SINCOS_ADDS;
  memcpy(&(ix), &ux.asreal, sizeof(ix));
#elif EAH_SINCOS_F2IBITS == EAH_SINCOS_F2IBITS_UNION
  ux.asreal = x + SINCOS_ADDS;
#ifdef __BIG_ENDIAN__
  ix = ux.as2int.low;
#else /* BIG_ENDIAN */
  ix = ux.as2int.high;
#endif /* BIG_ENDIAN */
#endif /* EAH_SINCOS_F2IBITS */

  i  = ix & SINCOS_MASK1;
  n  = ix & SINCOS_MASK2;
  i  = i >> SINCOS_SHIFT;
  
  (*sin2pix) = sincosLUTbase[i]  + n * sincosLUTdiff[i];
  (*cos2pix) = cosbase[i]        + n * cosdiff[i];

  return XLAL_SUCCESS;
}

#else  /* EAH_SINCOS_VARIANT */
#error no valid EAH_SINCOS_VARIANT specified
#endif /* EAH_SINCOS_VARIANT */
