//
// Copyright (C) 2012, 2013, 2014 Karl Wette
// Copyright (C) 2007, 2008, 2009, 2010, 2012 Bernd Machenschalk
// Copyright (C) 2007 Chris Messenger
// Copyright (C) 2006 John T. Whelan, Badri Krishnan
// Copyright (C) 2005, 2006, 2007, 2009, 2010, 2012 Reinhard Prix
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

// This file implements the F-statistic demodulation algorithm. It is not compiled directly, but
// included from ComputeFstat.c

///////////////////////////////////////////////////////////////////
//////////////////// Old demodulation API code ////////////////////
///////////////////////////////////////////////////////////////////

#include <lal/LogPrintf.h>

/*----- SIN/COS Lookup table code inserted here -----*/
/* In particular this defines
  - a macro SINCOS_TRIM_X(y,x) which trims the value x to interval [0..2)
  - a function void local_sin_cos_2PI_LUT_init(void) that inits the lookup tables
  - a function local_sin_cos_2PI_LUT_trimmed(*sin,*cos,x) that uses the
    lookup tables to evaluate sin and cos values of 2*Pi*x
  - macros to be used in "hotloop" variants to interleave hotloop & sincos calculations
*/
#include "SinCosLUT.i"

/* Macro to check return value of local_sin_cos_2PI_LUT_trimmed only in DEBUG mode */
#ifndef LAL_NDEBUG
#define SINCOS_2PI_TRIMMED(s,c,x)               \
  if ( local_sin_cos_2PI_LUT_trimmed(s,c,x) ) { \
    XLAL_ERROR ( XLAL_EFUNC);                   \
  }
#else
#define SINCOS_2PI_TRIMMED(s,c,x)               \
  local_sin_cos_2PI_LUT_trimmed(s,c,x)
#endif

/* LocalComputeFaFb: fixed DTERMS to allow for loop unrolling */
#define DTERMS 8
const UINT4 OptimisedHotloopDterms = DTERMS;

/* select optimized hotloop variant to be included and used in LocalXLALComputeFaFb() */
#ifdef AUTOVECT_HOTLOOP
#define OPT_DEMOD_SOURCE "ComputeFstat_Demod_autovect.i"
#elif __ALTIVEC__
#define OPT_DEMOD_SOURCE "ComputeFstat_Demod_altivec.i"
#elif __SSE__ && defined(_MSC_VER)
#define OPT_DEMOD_SOURCE "ComputeFstat_Demod_sse_msc.i"
#elif __SSE__ && defined(__OPTIMIZE__)
#define OPT_DEMOD_SOURCE "ComputeFstat_Demod_precalc.i"
#else
#define OPT_DEMOD_SOURCE "ComputeFstat_Demod_generic.i"
#endif

#if defined(__SSE__) && defined(_MSC_VER)
const int DEMODHL_BEST = DEMODHL_SIMD_SSE_MSC;
#elif defined(__SSE__)
const int DEMODHL_BEST = DEMODHL_SIMD_SSE_PRECALC;
#elif defined(__ALTIVEC__)
const int DEMODHL_BEST = DEMODHL_SIMD_ALTIVEC;
#else
const int DEMODHL_BEST = DEMODHL_AKOSGENERIC20;
#endif

/* record which optimized hotloop variant was selected for use in LocalXLALComputeFaFb() */
const char *const OptimisedHotloopSource = OPT_DEMOD_SOURCE;

#define LD_SMALL4       (2.0e-4)                /* "small" number for REAL4*/
#define OOTWOPI         (1.0 / LAL_TWOPI)       /* 1/2pi */
#define TWOPI_FLOAT     6.28318530717958f       /* single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)    /* single-precision 1 / (2pi) */

/* Type containing F-statistic proper plus the two complex amplitudes Fa and Fb (for ML-estimators) */
typedef struct tagFcomponents {
  REAL8 F;                              /* F-statistic value */
  REAL8 FX[PULSAR_MAX_DETECTORS];               /* vector of single-detector F-statistic values (array of fixed size) */
  UINT4 numDetectors;                   /* number of detectors = effective vector length. numDetectors=0 should make all code ignore the FX field. */
  LIGOTimeGPS refTime;                  /* 'internal' refTime used to compute the F-statistic: only relevant for phase of complex amplitudes {Fa,Fb} */
  COMPLEX16 Fa;                         /* complex amplitude Fa */
  COMPLEX16 Fb;                         /* complex amplitude Fb */
  MultiFstatAtomVector *multiFstatAtoms;/* per-IFO, per-SFT arrays of F-stat 'atoms', ie quantities required to compute F-stat */
} Fcomponents;

/* [opaque] type holding a ComputeFBuffer for use in the resampling F-stat codes */
typedef struct tagComputeFBuffer_RS ComputeFBuffer_RS;

/* Extra parameters controlling the actual computation of F */
typedef struct tagComputeFParams {
  UINT4 Dterms;         /* how many terms to keep in the Dirichlet kernel (~16 is usually fine) */
  SSBprecision SSBprec; /* whether to use full relativist SSB-timing, or just simple Newtonian */
  ComputeFBuffer_RS *buffer; /* buffer for storing pre-resampled timeseries (used for resampling implementation) */
  const EphemerisData *edat;   /* ephemeris data for re-computing multidetector states */
  BOOLEAN returnAtoms;  /* whether or not to return the 'FstatAtoms' used to compute the F-statistic */
  BOOLEAN returnSingleF; /* in multi-detector case, whether or not to also return the single-detector Fstats computed from the atoms */
  DemodHLType demodHL; 	//!< which hotloop variant to use in Demod: 'vanilla', 'generic', 'autovect', 'SIMD', ...
} ComputeFParams;

/* Struct holding buffered ComputeFStat()-internal quantities to avoid unnecessarily
 * recomputing things that depend ONLY on the skyposition and detector-state series (but not on the spins).
 * For the first call of ComputeFStat() the pointer-entries should all be NULL.
 */
typedef struct tagComputeFBuffer {
  const MultiDetectorStateSeries *multiDetStates;/* buffer for each detStates (store pointer) and skypos */
  REAL8 Alpha, Delta;                           /* skyposition of candidate */
  MultiSSBtimes *multiSSB;
  MultiSSBtimes *multiBinary;
  MultiAMCoeffs *multiAMcoef;
} ComputeFBuffer;

// ----- local prototypes ----------
static int LocalXLALComputeFaFb ( Fcomponents *FaFb, const SFTVector *sfts, const PulsarSpins fkdot, const SSBtimes *tSSB, const AMCoeffs *amcoe, const ComputeFParams *params);
static int XLALComputeFaFb ( Fcomponents *FaFb, const SFTVector *sfts, const PulsarSpins fkdot, const SSBtimes *tSSB, const AMCoeffs *amcoe, const ComputeFParams *params );
static int ComputeFStat ( Fcomponents *Fstat, const PulsarDopplerParams *doppler, const MultiSFTVector *multiSFTs, const MultiNoiseWeights *multiWeights, const MultiDetectorStateSeries *multiDetStates, const ComputeFParams *params, ComputeFBuffer *cfBuffer );
static inline void DestroyFstatInput_Demod ( FstatInput_Demod* demod );

// ----- NEW F-stat Demod API ----------
static int ComputeFstat_Demod ( FstatResults* Fstats, const FstatInput_Common *common, FstatInput_Demod* demod );

static void XLALEmptyComputeFBuffer ( ComputeFBuffer *cfb);

// ----- function definitions ----------
/* Destruction of a ComputeFBuffer *contents*,
 * i.e. the multiSSB and multiAMcoeff, while the
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...) */
static void
XLALEmptyComputeFBuffer ( ComputeFBuffer *cfb)
{
  XLALDestroyMultiSSBtimes ( cfb->multiSSB );
  cfb->multiSSB = NULL;
  XLALDestroyMultiSSBtimes ( cfb->multiBinary );
  cfb->multiBinary = NULL;
  XLALDestroyMultiAMCoeffs ( cfb->multiAMcoef );
  cfb->multiAMcoef = NULL;
  return;
} // XLALDestroyComputeFBuffer()

/* Function to compute (multi-IFO) F-statistic for given parameter-space point \a doppler,
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
static int
ComputeFStat ( Fcomponents *Fstat,                              /* [out] Fstatistic + Fa, Fb */
               const PulsarDopplerParams *doppler,              /* parameter-space point to compute F for */
               const MultiSFTVector *multiSFTs,                 /* normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
               const MultiNoiseWeights *multiWeights,           /* noise-weights of all SFTs */
               const MultiDetectorStateSeries *multiDetStates,  /* 'trajectories' of the different IFOs */
               const ComputeFParams *params,                    /* addition computational params */
               ComputeFBuffer *cfBuffer                         /* CF-internal buffering structure */
               )
{
  /* check input */
  XLAL_CHECK ( Fstat != NULL, XLAL_EINVAL );
  XLAL_CHECK ( multiSFTs != NULL, XLAL_EINVAL );
  XLAL_CHECK ( doppler != NULL, XLAL_EINVAL );
  XLAL_CHECK ( multiDetStates != NULL, XLAL_EINVAL );
  XLAL_CHECK ( params != NULL, XLAL_EINVAL );

  UINT4 numDetectors = multiSFTs->length;
  XLAL_CHECK ( multiDetStates->length == numDetectors, XLAL_EINVAL );
  XLAL_CHECK ( multiWeights==NULL || (multiWeights->length == numDetectors), XLAL_EINVAL );

  Fcomponents XLAL_INIT_DECL(retF);
  MultiSSBtimes *multiSSB = NULL;
  MultiSSBtimes *multiBinary = NULL;
  const MultiSSBtimes *multiSSBTotal = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;

  /* ----- prepare return of 'FstatAtoms' if requested */
  if ( params->returnAtoms )
    {
      XLAL_CHECK ( (retF.multiFstatAtoms = XLALMalloc ( sizeof(*retF.multiFstatAtoms) )) != NULL, XLAL_ENOMEM );
      retF.multiFstatAtoms->length = numDetectors;
      XLAL_CHECK ( (retF.multiFstatAtoms->data = XLALMalloc ( numDetectors * sizeof(*retF.multiFstatAtoms->data) )) != NULL, XLAL_ENOMEM );
    } /* if returnAtoms */

  /* ----- check if that skyposition SSB+AMcoef were already buffered */
  if ( cfBuffer
       && ( cfBuffer->multiDetStates == multiDetStates )
       && ( cfBuffer->Alpha == doppler->Alpha )
       && ( cfBuffer->Delta == doppler->Delta )
       && cfBuffer->multiSSB
       && cfBuffer->multiAMcoef )
    { /* yes ==> reuse */
      multiSSB = cfBuffer->multiSSB;
      multiAMcoef = cfBuffer->multiAMcoef;
    } /* if have buffered stuff to reuse */
  else
    {
      SkyPosition skypos;
      skypos.system =   COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;
      /* compute new AM-coefficients and SSB-times */
      XLAL_CHECK ( (multiSSB = XLALGetMultiSSBtimes ( multiDetStates, skypos, doppler->refTime, params->SSBprec )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (multiAMcoef = XLALComputeMultiAMCoeffs ( multiDetStates, multiWeights, skypos )) != NULL, XLAL_EFUNC );

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

    } /* could not reuse previously buffered quantites */

  /* new orbital parameter corrections if not already buffered */
  if ( doppler->asini > 0 )
    {
      /* compute binary time corrections to the SSB time delays and SSB time derivitive */
      XLAL_CHECK ( XLALAddMultiBinaryTimes ( &multiBinary, multiSSB, doppler ) == XLAL_SUCCESS, XLAL_EFUNC );
      multiSSBTotal = multiBinary;
    }
  else {
    multiSSBTotal = multiSSB;
  }

  REAL8 Ad = multiAMcoef->Mmunu.Ad;
  REAL8 Bd = multiAMcoef->Mmunu.Bd;
  REAL8 Cd = multiAMcoef->Mmunu.Cd;
  REAL8 Dd_inv = 1.0 / multiAMcoef->Mmunu.Dd;
  REAL8 Ed = 0;

  /* if requested, prepare for returning single-IFO F-stat vector */
  if ( params->returnSingleF )
    {
      retF.numDetectors = numDetectors;
      XLAL_CHECK ( numDetectors <= PULSAR_MAX_DETECTORS, XLAL_EINVAL, "numDetectors = %d exceeds currently allowed upper value (%d) for returnSingleF=TRUE\n", numDetectors, PULSAR_MAX_DETECTORS );
    }

  /* ----- loop over detectors and compute all detector-specific quantities ----- */
  for ( UINT4 X=0; X < numDetectors; X ++)
    {
      Fcomponents XLAL_INIT_DECL(FcX);      /* for detector-specific FaX, FbX */

      if ( (params->Dterms != DTERMS) || params->returnAtoms ) {
        XLAL_CHECK ( XLALComputeFaFb      (&FcX, multiSFTs->data[X], doppler->fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], params) == XLAL_SUCCESS, XLAL_EFUNC );
      } else {
        XLAL_CHECK ( LocalXLALComputeFaFb (&FcX, multiSFTs->data[X], doppler->fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], params) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      if ( params->returnAtoms )
        {
          retF.multiFstatAtoms->data[X] = FcX.multiFstatAtoms->data[0];     /* copy pointer to IFO-specific Fstat-atoms 'contents' */
          /* free 'container', but not *contents*, which have been linked above */
          XLALFree ( FcX.multiFstatAtoms->data );
          XLALFree ( FcX.multiFstatAtoms );
        }

      XLAL_CHECK ( isfinite(creal(FcX.Fa)) && isfinite(cimag(FcX.Fa)) && isfinite(creal(FcX.Fb)) && isfinite(cimag(FcX.Fb)), XLAL_EFPOVRFLW );

      /* compute single-IFO F-stats, if requested */
      if ( params->returnSingleF )
        {
         REAL8 AdX = multiAMcoef->data[X]->A;
         REAL8 BdX = multiAMcoef->data[X]->B;
         REAL8 CdX = multiAMcoef->data[X]->C;
         REAL8 DdX_inv = 1.0 / multiAMcoef->data[X]->D;
         REAL8 EdX = 0;

         REAL8 FXa_re = creal(FcX.Fa);
         REAL8 FXa_im = cimag(FcX.Fa);
         REAL8 FXb_re = creal(FcX.Fb);
         REAL8 FXb_im = cimag(FcX.Fb);

         /* compute final single-IFO F-stat */
         retF.FX[X] = DdX_inv * (  BdX * ( SQ(FXa_re) + SQ(FXa_im) )
                                   + AdX * ( SQ(FXb_re) + SQ(FXb_im) )
                                   - 2.0 * CdX * (   FXa_re * FXb_re + FXa_im * FXb_im )
                                   - 2.0 * EdX * ( - FXa_re * FXb_im + FXa_im * FXb_re )		// nonzero only in RAA case where Ed!=0
                                   );
        } /* if returnSingleF */

      /* Fa = sum_X Fa_X */
      retF.Fa += FcX.Fa;

      /* Fb = sum_X Fb_X */
      retF.Fb += FcX.Fb;

    } /* for  X < numDetectors */

  /* ----- compute final Fstatistic-value ----- */

  /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
   * therefore there is a factor of 2 difference with respect to the equations in JKS, which
   * where based on the single-sided PSD.
   */
  REAL8 Fa_re = creal(retF.Fa);
  REAL8 Fa_im = cimag(retF.Fa);
  REAL8 Fb_re = creal(retF.Fb);
  REAL8 Fb_im = cimag(retF.Fb);

  retF.F = Dd_inv * (  Bd * ( SQ(Fa_re) + SQ(Fa_im) )
                       + Ad * ( SQ(Fb_re) + SQ(Fb_im) )
                       - 2.0 * Cd * ( Fa_re * Fb_re + Fa_im * Fb_im )
                       - 2.0 * Ed * ( - Fa_re * Fb_im + Fa_im * Fb_re )		// nonzero only in RAA case where Ed!=0
                       );

  /* set correct F-stat reference time (taken from template 'doppler') [relevant only for phase of {Fa,Fb}] */
  retF.refTime = doppler->refTime;

  /* free memory if no buffer was available */
  if ( !cfBuffer )
    {
      XLALDestroyMultiSSBtimes ( multiSSB );
      XLALDestroyMultiAMCoeffs ( multiAMcoef );
    } /* if !cfBuffer */

  /* this always needs to be free'ed, as it's no longer buffered */
  XLALDestroyMultiSSBtimes ( multiBinary );

  /* return final Fstat result */
  (*Fstat) = retF;

  return XLAL_SUCCESS;

} // ComputeFStat()


/* Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 */
static int
XLALComputeFaFb ( Fcomponents *FaFb,                    /* [out] Fa,Fb (and possibly atoms) returned */
                  const SFTVector *sfts,                /* [in] input SFTs */
                  const PulsarSpins fkdot,              /* [in] frequency and derivatives fkdot = d^kf/dt^k */
                  const SSBtimes *tSSB,                 /* [in] SSB timing series for particular sky-direction */
                  const AMCoeffs *amcoe,                /* [in] antenna-pattern coefficients for this sky-direction */
                  const ComputeFParams *params )        /* addition computational params */
{
  UINT4 alpha;                  /* loop index over SFTs */
  UINT4 spdnOrder;              /* maximal spindown-orders */
  UINT4 numSFTs;                /* number of SFTs (M in the Notes) */
  COMPLEX16 Fa, Fb;
  REAL8 Tsft;                   /* length of SFTs in seconds */
  INT4 freqIndex0;              /* index of first frequency-bin in SFTs */
  INT4 freqIndex1;              /* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;           /* pointer to alpha-arrays over a and b */
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
#endif

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;
  {
    REAL8 dFreq = sfts->data[0].deltaF;
    freqIndex0 = lround ( sfts->data[0].f0 / dFreq ); /* lowest freqency-index */
    freqIndex1 = freqIndex0 + sfts->data[0].data->length;
  }

  /* ----- prepare return of 'FstatAtoms' if requested */
  if ( params->returnAtoms )
    {
      if ( (FaFb->multiFstatAtoms = LALMalloc ( sizeof(*FaFb->multiFstatAtoms) )) == NULL ){
        XLAL_ERROR ( XLAL_ENOMEM );
      }
      FaFb->multiFstatAtoms->length = 1;        /* in this function: single-detector only */
      if ( (FaFb->multiFstatAtoms->data = LALMalloc ( 1 * sizeof( *FaFb->multiFstatAtoms->data) )) == NULL ){
        LALFree (FaFb->multiFstatAtoms);
        XLAL_ERROR ( XLAL_ENOMEM );
      }
      if ( (FaFb->multiFstatAtoms->data[0] = XLALCreateFstatAtomVector ( numSFTs )) == NULL ) {
        LALFree ( FaFb->multiFstatAtoms->data );
        LALFree ( FaFb->multiFstatAtoms );
        XLAL_ERROR( XLAL_ENOMEM );
      }

      FaFb->multiFstatAtoms->data[0]->TAtom = Tsft;     /* time-baseline of returned atoms is Tsft */

    } /* if returnAtoms */

  /* ----- find highest non-zero spindown-entry */
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
      REAL4 a_alpha, b_alpha;

      INT4 kstar;               /* central frequency-bin k* = round(xhat_alpha) */
      INT4 k0, k1;

      COMPLEX8 *Xalpha = SFT_al->data->data; /* pointer to current SFT-data */
      COMPLEX8 *Xalpha_l;       /* pointer to frequency-bin k in current SFT */
      REAL4 s_alpha=0, c_alpha=0;/* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
      REAL4 realQ, imagQ;       /* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realXP, imagXP;     /* Re/Im of sum_k X_ak * P_ak */
      REAL4 realQXP, imagQXP;   /* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha, kappa_max, kappa_star;
      COMPLEX8 Fa_alpha, Fb_alpha;

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
            XLALPrintError("XLALComputeFaFb() returned non-finite: U_alpha=%f, V_alpha=%f, pn=%f, qn=%f, Sn=%f, Tn=%f\n",
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

      Fa_alpha = crectf( a_alpha * realQXP, a_alpha * imagQXP );
      Fa += Fa_alpha;

      Fb_alpha = crectf( b_alpha * realQXP, b_alpha * imagQXP );
      Fb += Fb_alpha;

      /* store per-SFT F-stat 'atoms' for transient-CW search */
      if ( params->returnAtoms )
        {
          FaFb->multiFstatAtoms->data[0]->data[alpha].timestamp = (UINT4)XLALGPSGetREAL8( &SFT_al->epoch );
          FaFb->multiFstatAtoms->data[0]->data[alpha].a2_alpha   = a_alpha * a_alpha;
          FaFb->multiFstatAtoms->data[0]->data[alpha].b2_alpha   = b_alpha * b_alpha;
          FaFb->multiFstatAtoms->data[0]->data[alpha].ab_alpha   = a_alpha * b_alpha;
          FaFb->multiFstatAtoms->data[0]->data[alpha].Fa_alpha   = norm * Fa_alpha;
          FaFb->multiFstatAtoms->data[0]->data[alpha].Fb_alpha   = norm * Fb_alpha;
        }

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

} // XLALComputeFaFb()

/* Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 */
static int
LocalXLALComputeFaFb ( Fcomponents *FaFb,
                       const SFTVector *sfts,
                       const PulsarSpins fkdot,
                       const SSBtimes *tSSB,
                       const AMCoeffs *amcoe,
                       const ComputeFParams *params)       /* addition computational params */
{
  UINT4 alpha;                  /* loop index over SFTs */
  UINT4 spdnOrder;              /* maximal spindown-orders */
  UINT4 numSFTs;                /* number of SFTs (M in the Notes) */
  COMPLEX16 Fa, Fb;
  REAL8 Tsft;                   /* length of SFTs in seconds */
  INT4 freqIndex0;              /* index of first frequency-bin in SFTs */
  INT4 freqIndex1;              /* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;           /* pointer to alpha-arrays over a and b */
  REAL8 *DeltaT_al, *Tdot_al;   /* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;              /* SFT alpha  */

  REAL8 norm = OOTWOPI;

  XLAL_CHECK ( params->Dterms == DTERMS, XLAL_EDOM );
#ifndef LAL_NDEBUG
  /* ----- check validity of input */
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

#endif

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;
  {
    REAL8 dFreq = sfts->data[0].deltaF;
    freqIndex0 = lround ( sfts->data[0].f0 / dFreq ); /* lowest freqency-index */
    freqIndex1 = freqIndex0 + sfts->data[0].data->length;
  }

  static int firstcall = 1; /* for sin/cos lookup table initialization */
  if (firstcall)
    {
      /* init sin/cos lookup tables */
      local_sin_cos_2PI_LUT_init();

      /* make sure Dterms is what we expect */
      XLAL_CHECK (DTERMS == params->Dterms, XLAL_EINVAL, "LocalXLALComputeFaFb() has been compiled with fixed DTERMS (%d) != params->Dtems (%d)\n", DTERMS, params->Dterms );
      firstcall = 0;
    }

  /* find highest non-zero spindown-entry */
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot[spdnOrder] != 0.0 )
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
      REAL4 a_alpha, b_alpha;

      INT4 kstar;               /* central frequency-bin k* = round(xhat_alpha) */
      INT4 k0, k1;

      COMPLEX8 *Xalpha = SFT_al->data->data; /* pointer to current SFT-data */
      REAL4 s_alpha, c_alpha;   /* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
      REAL4 realQ, imagQ;       /* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realXP, imagXP;     /* re/im of sum_k X_ak * P_ak */
      REAL4 realQXP, imagQXP;   /* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha, kappa_star;

      /* ----- calculate lambda_alpha */
      {
        UINT4 s;                /* loop-index over spindown-order */
        REAL8 phi_alpha, Dphi_alpha, DT_al;
        REAL8 Tas;      /* temporary variable to calculate (DeltaT_alpha)^s */
        REAL8 TAS_invfact_s;

        /* init for s=0 */
        phi_alpha = 0.0;
        Dphi_alpha = 0.0;
        DT_al = (*DeltaT_al);
        Tas = 1.0;              /* DeltaT_alpha ^ 0 */
        TAS_invfact_s=1.0;    /* TAS / s! */

        for (s=0; s <= spdnOrder; s++) {
          REAL8 fsdot = fkdot[s];
          Dphi_alpha += fsdot * TAS_invfact_s;  /* here: DT^s/s! */
#ifdef EAH_CHECK_FINITE_DPHI
          if (!finite(Dphi_alpha)) {
            LogPrintf(LOG_CRITICAL, "non-finite Dphi_alpha:%e, alpha:%d, spind#:%d, fkdot:%e, Tas:%e, LAL_FACT_INV[s]:%e, LAL_FACT_INV[s+1]:%e, phi_alpha:%e. DT_al:%e\n",
                      Dphi_alpha, alpha, s, fkdot[s], Tas, LAL_FACT_INV[s], LAL_FACT_INV[s+1], phi_alpha, DT_al);
            XLAL_ERROR("LocalXLALComputeFaFb", XLAL_EDOM);
          }
#endif
          Tas *= DT_al;                         /* now: DT^(s+1) */
          TAS_invfact_s= Tas * LAL_FACT_INV[s+1];
          phi_alpha += fsdot * TAS_invfact_s;
        } /* for s <= spdnOrder */

        /* Step 3: apply global factors to complete Dphi_alpha */
        Dphi_alpha *= Tsft * (*Tdot_al);                /* guaranteed > 0 ! */
#ifndef EAH_NO_CHECK_FINITE_DPHI
        if (!finite(Dphi_alpha)) {
          LogPrintf(LOG_CRITICAL, "non-finite Dphi_alpha:%e, alpha:%d, Tsft:%e, Tkdot_al:%e Tas:%e, DT_al:%e\n",
                    Dphi_alpha, alpha, Tsft, (*Tdot_al), Tas, DT_al);
          XLAL_ERROR( XLAL_EDOM );
        }
#endif
        lambda_alpha = 0.5 * Dphi_alpha - phi_alpha;

        /* FIXME: that might be possible to do faster */
        kstar = (INT4) (Dphi_alpha);    /* k* = floor(Dphi_alpha*chi) for positive Dphi */
        kappa_star = Dphi_alpha - 1.0 * kstar;  /* remainder of Dphi_alpha: >= 0 ! */

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
            XLAL_ERROR(XLAL_EDOM);
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
        /* somehow the branch prediction of gcc-4.1.2 terribly failes
            with the current case distinction in the hot-loop,
            having a severe impact on runtime of the E@H Linux App.
            So let's allow to give gcc a hint which path has a higher probablility */
        if (__builtin_expect((kappa_star > LD_SMALL4) && (kappa_star < 1.0 - LD_SMALL4), (0==0)))
#else
        if ((kappa_star > LD_SMALL4) && (kappa_star < 1.0 - LD_SMALL4))
#endif
          {
          /* WARNING: all current optimized loops rely on current implementation of COMPLEX8 and DTERMS == 8 */

#include OPT_DEMOD_SOURCE

          } /* if |remainder| > LD_SMALL4 */
        else
          { /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
            UINT4 ind0;

            /* realQ/imagQ are calculated in the hotloop in the other case; in this case we have to do it too */
            SINCOS_TRIM_X (lambda_alpha,lambda_alpha);
            SINCOS_2PI_TRIMMED( &imagQ, &realQ, lambda_alpha );

            if ( kappa_star <= LD_SMALL4 ) {
              ind0 = DTERMS - 1;
            } else {
              ind0 = DTERMS;
            }

            realXP = TWOPI_FLOAT * crealf(Xalpha_l[ind0]);
            imagXP = TWOPI_FLOAT * cimagf(Xalpha_l[ind0]);

          } /* if |remainder| <= LD_SMALL4 */
      }

      /* real- and imaginary part of e^{-i 2 pi lambda_alpha } */

      realQXP = realQ * realXP - imagQ * imagXP;
      imagQXP = realQ * imagXP + imagQ * realXP;


      /* we're done: ==> combine these into Fa and Fb */

      a_alpha = (*a_al);
      b_alpha = (*b_al);

      Fa += crect( a_alpha * realQXP, a_alpha * imagQXP );

      Fb += crect( b_alpha * realQXP, b_alpha * imagQXP );

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

} // LocalXLALComputeFaFb()

///////////////////////////////////////////////////////////////////
//////////////////// New demodulation API code ////////////////////
///////////////////////////////////////////////////////////////////

struct tagFstatInput_Demod {
  UINT4 Dterms;                                 // Number of terms to keep in Dirichlet kernel
  DemodHLType demodHL; 	                        // Which hotloop variant to use
  MultiSFTVector *multiSFTs;                    // Input multi-detector SFTs
  ComputeFParams params;                        // Additional parameters for ComputeFStat()
  ComputeFBuffer buffer;                        // Internal buffer for ComputeFStat()
};

static inline void
DestroyFstatInput_Demod(
  FstatInput_Demod* demod
  )
{
  XLALDestroyMultiSFTVector(demod->multiSFTs);
  XLALEmptyComputeFBuffer(&demod->buffer);
  XLALFree(demod);
} // DestroyFstatInput_Demod()

FstatInput*
XLALCreateFstatInput_Demod(
  const UINT4 Dterms,
  const DemodHLType demodHL
  )
{

  // Check input
  XLAL_CHECK_NULL(Dterms > 0, XLAL_EINVAL);

  // Allocate input data struct
  FstatInput* input = XLALCalloc(1, sizeof(FstatInput));
  XLAL_CHECK_NULL(input != NULL, XLAL_ENOMEM);
  input->demod = XLALCalloc(1, sizeof(FstatInput_Demod));
  XLAL_CHECK_NULL(input->demod != NULL, XLAL_ENOMEM);

  // Save parameters
  input->demod->Dterms = Dterms;
  input->demod->demodHL = demodHL;

  return input;

}

static int
SetupFstatInput_Demod(
  FstatInput_Demod *demod,
  const FstatInput_Common *common,
  MultiSFTVector *multiSFTs
  )
{

  // Check input
  XLAL_CHECK(common != NULL, XLAL_EFAULT);
  XLAL_CHECK(demod != NULL, XLAL_EFAULT);
  XLAL_CHECK(multiSFTs != NULL, XLAL_EFAULT);

  // Save pointer to SFTs
  demod->multiSFTs = multiSFTs;

  // Set parameters to pass to ComputeFStat()
  demod->params.Dterms = demod->Dterms;
  demod->params.SSBprec = common->SSBprec;
  demod->params.buffer = NULL;
  demod->params.edat = common->ephemerides;
  demod->params.demodHL = demod->demodHL;

  return XLAL_SUCCESS;

}

static int
GetFstatExtraBins_Demod(
  FstatInput_Demod* demod
  )
{


  // Check input
  XLAL_CHECK(demod != NULL, XLAL_EFAULT);

  // Demodulation requires 'Dterms' extra frequency bins
  return demod->Dterms;

} // XLALSetupFstat_Demod()

static int
ComputeFstat_Demod(
  FstatResults* Fstats,
  const FstatInput_Common *common,
  FstatInput_Demod* demod
  )
{

  // Check input
  XLAL_CHECK(Fstats != NULL, XLAL_EFAULT);
  XLAL_CHECK(common != NULL, XLAL_EFAULT);
  XLAL_CHECK(demod != NULL, XLAL_EFAULT);

  // Get which F-statistic quantities to compute
  const FstatQuantities whatToCompute = Fstats->whatWasComputed;

  // Check which quantities can be computed
  XLAL_CHECK(!(whatToCompute & FSTATQ_FAFB_PER_DET), XLAL_EINVAL, "Demodulation does not currently support Fa & Fb per detector");

  // Set parameters to pass to ComputeFStat()
  demod->params.returnSingleF = (whatToCompute & FSTATQ_2F_PER_DET);
  demod->params.returnAtoms = (whatToCompute & FSTATQ_ATOMS_PER_DET);

  // Save local copy of doppler point and starting frequency
  PulsarDopplerParams thisPoint = Fstats->doppler;
  const REAL8 fStart = thisPoint.fkdot[0];

  // Call ComputeFStat() for each frequency bin
  for (UINT4 k = 0; k < Fstats->numFreqBins; ++k) {

    // Set frequency to search at
    thisPoint.fkdot[0] = fStart + k * Fstats->dFreq;

    // Call ComputeFStat()
    Fcomponents Fcomp;
    XLAL_CHECK ( ComputeFStat( &Fcomp, &thisPoint, demod->multiSFTs, common->noiseWeights, common->detectorStates, &demod->params, &demod->buffer) == XLAL_SUCCESS, XLAL_EFUNC );

    // Return multi-detector 2F
    if (whatToCompute & FSTATQ_2F) {
      Fstats->twoF[k] = 2.0 * Fcomp.F;   // *** Return value of 2F ***
    }

    // Return multi-detector Fa & Fb
    if (whatToCompute & FSTATQ_FAFB) {
      Fstats->FaFb[k].Fa = Fcomp.Fa;
      Fstats->FaFb[k].Fb = Fcomp.Fb;
    }

    // Return 2F per detector
    if (whatToCompute & FSTATQ_2F_PER_DET) {
      for (UINT4 X = 0; X < Fstats->numDetectors; ++X) {
        Fstats->twoFPerDet[X][k] = 2.0 * Fcomp.FX[X];   // *** Return value of 2F ***
      }
    }

    // Return Fa & Fb per detector
    if (whatToCompute & FSTATQ_FAFB_PER_DET) {
      XLAL_ERROR(XLAL_EFAILED, "Unimplemented!");
    }

    // Return F-atoms per detector
    if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
      XLALDestroyMultiFstatAtomVector(Fstats->multiFatoms[k]);
      Fstats->multiFatoms[k] = Fcomp.multiFstatAtoms;
    }

  } // for k < Fstats->numFreqBins

  // Return amplitude modulation coefficients
  if (demod->buffer.multiAMcoef != NULL) {
    Fstats->Mmunu.Ad = demod->buffer.multiAMcoef->Mmunu.Ad;
    Fstats->Mmunu.Bd = demod->buffer.multiAMcoef->Mmunu.Bd;
    Fstats->Mmunu.Cd = demod->buffer.multiAMcoef->Mmunu.Cd;
    Fstats->Mmunu.Ed = 0;
    Fstats->Mmunu.Dd = demod->buffer.multiAMcoef->Mmunu.Dd;
    Fstats->Mmunu.Sinv_Tsft = demod->buffer.multiAMcoef->Mmunu.Sinv_Tsft;
  }

  return XLAL_SUCCESS;

} // ComputeFstat_Demod()
