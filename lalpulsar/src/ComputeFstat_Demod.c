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
#include <lal/CWFastMath.h>
#include "config.h"

/* somehow the branch prediction of gcc-4.1.2 terribly failes
   with the current case distinction in the hot-loop,
   having a severe impact on runtime of the E@H Linux App.
   So let's allow to give gcc a hint which path has a higher probablility
*/
#ifdef __GNUC__
#define likely(x)       __builtin_expect((x),1)
#else
#define likely(x)       (x)
#endif

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
  FstatMethodType FstatMethod; //!< which Fstat-algorithm to use: see documentation for FstatMethodType
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
static int ComputeFStat ( Fcomponents *Fstat, const PulsarDopplerParams *doppler, const MultiSFTVector *multiSFTs, const MultiNoiseWeights *multiWeights, const MultiDetectorStateSeries *multiDetStates, const ComputeFParams *params, ComputeFBuffer *cfBuffer );
static inline void DestroyFstatInput_Demod ( FstatInput_Demod* demod );

// ----- NEW F-stat Demod API ----------
static int ComputeFstat_Demod ( FstatResults* Fstats, const FstatInput_Common *common, FstatInput_Demod* demod );
static void XLALEmptyComputeFBuffer ( ComputeFBuffer *cfb);

// ----- define various variants of ComputeFaFb() using different hotloops
// ComputeFaFb: DTERMS define used for loop unrolling in some hotloop variants
#define DTERMS 8
#include "SinCosLUT.i"

// ----- old (pre-Akos) LALDemod hotloop variant (unrestricted Dterms) ----------
#define FUNC XLALComputeFaFb_Generic
#define HOTLOOP_SOURCE "ComputeFstat_DemodHL_Generic.i"
#define RUNTIME_CHECK XLAL_CHECK ( Dterms > 0, XLAL_EINVAL );
#include "ComputeFstat_Demod_ComputeFaFb.c"
// ----- Akos generic hotloop code (Dterms <= 20) ----------
#define FUNC XLALComputeFaFb_OptC
#define HOTLOOP_SOURCE "ComputeFstat_DemodHL_OptC.i"
#define RUNTIME_CHECK XLAL_CHECK ( Dterms <= 20, XLAL_EINVAL, "Selected Hotloop variant 'OptC' only works for Dterms <= 20, got %d\n", Dterms );
#include "ComputeFstat_Demod_ComputeFaFb.c"
// ----- Akos hotloop precalc SSE code (Dterms=8) ----------
#if defined(HAVE_SSE)
#define FUNC XLALComputeFaFb_SSE
#define HOTLOOP_SOURCE "ComputeFstat_DemodHL_SSE.i"
#define RUNTIME_CHECK XLAL_CHECK ( Dterms == 8, XLAL_EINVAL, "Selected Hotloop variant 'SSE' only works for Dterms == 8, got %d\n", Dterms );
#include "ComputeFstat_Demod_ComputeFaFb.c"
#else
#define XLALComputeFaFb_SSE(...) (XLALPrintError("Selected Hotloop variant 'SSE' unavailable\n") || XLAL_EFAILED)
#endif
// ----- Akos hotloop Altivec code (Dterms=8) ----------
#if (defined(HAVE_ALTIVEC) || defined(__ALTIVEC__))
#include <altivec.h>
#define FUNC XLALComputeFaFb_Altivec
#define HOTLOOP_SOURCE "ComputeFstat_DemodHL_Altivec.i"
#define RUNTIME_CHECK XLAL_CHECK ( Dterms == 8, XLAL_EINVAL, "Selected Hotloop variant 'Altivec' only works for Dterms == 8, got %d\n", Dterms );
#include "ComputeFstat_Demod_ComputeFaFb.c"
#else
#define XLALComputeFaFb_Altivec(...) (XLALPrintError("Selected Hotloop variant 'Altivec' unavailable\n") || XLAL_EFAILED)
#endif
// ------------------------------------------------------------

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

      // chose ComputeFaFb main function depending on selected hotloop variant
      switch ( params->FstatMethod )
        {
        case  FMETHOD_DEMOD_GENERIC:
          XLAL_CHECK ( XLALComputeFaFb_Generic ( &FcX, multiSFTs->data[X], doppler->fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], params) == XLAL_SUCCESS, XLAL_EFUNC );
          break;
        case FMETHOD_DEMOD_OPTC:
          XLAL_CHECK ( XLALComputeFaFb_OptC ( &FcX, multiSFTs->data[X], doppler->fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], params) == XLAL_SUCCESS, XLAL_EFUNC );
          break;
        case FMETHOD_DEMOD_SSE:
          XLAL_CHECK ( XLALComputeFaFb_SSE ( &FcX, multiSFTs->data[X], doppler->fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], params) == XLAL_SUCCESS, XLAL_EFUNC );
          break;
        case FMETHOD_DEMOD_ALTIVEC:
          XLAL_CHECK ( XLALComputeFaFb_Altivec ( &FcX, multiSFTs->data[X], doppler->fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], params) == XLAL_SUCCESS, XLAL_EFUNC );
          break;
        default:
          XLAL_ERROR ( XLAL_EINVAL, "Invalid Fstat-method %d!\n", params->FstatMethod );
          break;
        } // switch ( demodHL )

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


///////////////////////////////////////////////////////////////////
//////////////////// New demodulation API code ////////////////////
///////////////////////////////////////////////////////////////////

struct tagFstatInput_Demod {
  UINT4 Dterms;                                 // Number of terms to keep in Dirichlet kernel
  FstatMethodType FstatMethod;                      // Which Fstat algorithm to use
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
  const FstatMethodType FstatMethod
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
  input->demod->FstatMethod = FstatMethod;

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
  demod->params.FstatMethod = demod->FstatMethod;

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
