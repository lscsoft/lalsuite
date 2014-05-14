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

#include <lal/LogPrintf.h>
#include <lal/CWFastMath.h>
#include "config.h"

// ========== Demod internals ==========

// ----- local types ----------
struct tagFstatInput_Demod {
  UINT4 Dterms;                                 // Number of terms to keep in Dirichlet kernel
  FstatMethodType FstatMethod;                  // Which Fstat algorithm to use
  MultiSFTVector *multiSFTs;                    // Input multi-detector SFTs
  REAL8 prevAlpha, prevDelta;                   // previous skyposition computed
  MultiSSBtimes *prevMultiSSBtimes;		// previous multiSSB times, unique to skypos + SFTs
  MultiAMCoeffs *prevMultiAMcoef;		// previous AM-coeffs, unique to skypos + SFTs
};

// ----- local prototypes ----------
static void DestroyFstatInput_Demod ( FstatInput_Demod* demod );
static int ComputeFstat_Demod ( FstatResults* Fstats, const FstatInput_Common *common, FstatInput_Demod* demod );


// ========== define various hotloop variants of ComputeFaFb() ==========
// ComputeFaFb: DTERMS define used for loop unrolling in some hotloop variants
#define DTERMS 8
#define LD_SMALL4       (2.0e-4)                /* "small" number for REAL4*/
#define OOTWOPI         (1.0 / LAL_TWOPI)       /* 1/2pi */
#define TWOPI_FLOAT     6.28318530717958f       /* single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)    /* single-precision 1 / (2pi) */
/* somehow the branch prediction of gcc-4.1.2 terribly fails
 * So let's give gcc a hint which path has a higher probablility
 */
#ifdef __GNUC__
#define likely(x)       __builtin_expect((x),1)
#else
#define likely(x)       (x)
#endif

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
// ======================================================================


// ----- local function definitions ----------
static int
ComputeFstat_Demod ( FstatResults* Fstats,
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

  // handy shortcuts
  BOOLEAN returnSingleF = (whatToCompute & FSTATQ_2F_PER_DET);
  BOOLEAN returnAtoms = (whatToCompute & FSTATQ_ATOMS_PER_DET);
  UINT4 Dterms = demod->Dterms;
  PulsarDopplerParams thisPoint = Fstats->doppler;
  const REAL8 fStart = thisPoint.fkdot[0];
  const MultiSFTVector *multiSFTs = demod->multiSFTs;
  const MultiNoiseWeights *multiWeights = common->noiseWeights;
  const MultiDetectorStateSeries *multiDetStates = common->detectorStates;

  UINT4 numDetectors = multiSFTs->length;
  XLAL_CHECK ( multiDetStates->length == numDetectors, XLAL_EINVAL );
  XLAL_CHECK ( multiWeights==NULL || (multiWeights->length == numDetectors), XLAL_EINVAL );

  MultiSSBtimes *multiSSB = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  // ----- check if we have buffered SSB+AMcoef for current sky-position
  if ( (demod->prevAlpha == thisPoint.Alpha) && (demod->prevDelta == thisPoint.Delta ) && demod->prevMultiSSBtimes && demod->prevMultiAMcoef )
    { // if yes ==> reuse
      multiSSB    = demod->prevMultiSSBtimes;
      multiAMcoef = demod->prevMultiAMcoef;
    }
  else
    { // if not, compute SSB + AMcoef for this skyposition
      SkyPosition skypos;
      skypos.system = COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = thisPoint.Alpha;
      skypos.latitude  = thisPoint.Delta;
      XLAL_CHECK ( (multiSSB = XLALGetMultiSSBtimes ( multiDetStates, skypos, thisPoint.refTime, common->SSBprec )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (multiAMcoef = XLALComputeMultiAMCoeffs ( multiDetStates, multiWeights, skypos )) != NULL, XLAL_EFUNC );

      // store these for possible later re-use in buffer
      XLALDestroyMultiSSBtimes ( demod->prevMultiSSBtimes );
      demod->prevMultiSSBtimes = multiSSB;
      XLALDestroyMultiAMCoeffs ( demod->prevMultiAMcoef );
      demod->prevMultiAMcoef = multiAMcoef;
      demod->prevAlpha = thisPoint.Alpha;
      demod->prevDelta = thisPoint.Delta;
    } // if could not reuse previously buffered quantites

  MultiSSBtimes *multiBinary = NULL;
  MultiSSBtimes *multiSSBTotal = NULL;
  // handly binary-orbital timing corrections, if applicable
  if ( thisPoint.asini > 0 )
    {
      // compute binary time corrections to the SSB time delays and SSB time derivitive
      XLAL_CHECK ( XLALAddMultiBinaryTimes ( &multiBinary, multiSSB, &thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );
      multiSSBTotal = multiBinary;
    }
  else
    {
      multiSSBTotal = multiSSB;
    }

  // ---------- Compute F-stat for each frequency bin ----------
  for ( UINT4 k = 0; k < Fstats->numFreqBins; k++ )
    {
      // Set frequency to search at
      thisPoint.fkdot[0] = fStart + k * Fstats->dFreq;

      REAL8 F;                       	// F-statistic value
      REAL8 FX[PULSAR_MAX_DETECTORS];	// vector of single-detector F-statistic values (array of fixed size)
      COMPLEX16 Fa = 0;       		// complex amplitude Fa
      COMPLEX16 Fb = 0;                 // complex amplitude Fb
      MultiFstatAtomVector *multiFstatAtoms = NULL;	// per-IFO, per-SFT arrays of F-stat 'atoms', ie quantities required to compute F-stat

      // prepare return of 'FstatAtoms' if requested
      if ( returnAtoms )
        {
          XLAL_CHECK ( (multiFstatAtoms = XLALMalloc ( sizeof(*multiFstatAtoms) )) != NULL, XLAL_ENOMEM );
          multiFstatAtoms->length = numDetectors;
          XLAL_CHECK ( (multiFstatAtoms->data = XLALMalloc ( numDetectors * sizeof(*multiFstatAtoms->data) )) != NULL, XLAL_ENOMEM );
        } // if returnAtoms

      // loop over detectors and compute all detector-specific quantities
      for ( UINT4 X=0; X < numDetectors; X ++)
        {
          COMPLEX16 FaX, FbX;
          FstatAtomVector *FstatAtoms = NULL;
          FstatAtomVector **FstatAtoms_p = returnAtoms ? (&FstatAtoms) : NULL;

          // chose ComputeFaFb main function depending on selected hotloop variant
          switch ( demod->FstatMethod )
            {
            case  FMETHOD_DEMOD_GENERIC:
              XLAL_CHECK ( XLALComputeFaFb_Generic ( &FaX, &FbX, FstatAtoms_p, multiSFTs->data[X], thisPoint.fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], Dterms ) == XLAL_SUCCESS, XLAL_EFUNC );
              break;
            case FMETHOD_DEMOD_OPTC:
              XLAL_CHECK ( XLALComputeFaFb_OptC ( &FaX, &FbX, FstatAtoms_p, multiSFTs->data[X], thisPoint.fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], Dterms ) == XLAL_SUCCESS, XLAL_EFUNC );
              break;
            case FMETHOD_DEMOD_SSE:
              XLAL_CHECK ( XLALComputeFaFb_SSE ( &FaX, &FbX, FstatAtoms_p, multiSFTs->data[X], thisPoint.fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], Dterms) == XLAL_SUCCESS, XLAL_EFUNC );
              break;
            case FMETHOD_DEMOD_ALTIVEC:
              XLAL_CHECK ( XLALComputeFaFb_Altivec ( &FaX, &FbX, FstatAtoms_p, multiSFTs->data[X], thisPoint.fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], Dterms) == XLAL_SUCCESS, XLAL_EFUNC );
              break;
            default:
              XLAL_ERROR ( XLAL_EINVAL, "Invalid Fstat-method %d!\n", demod->FstatMethod );
              break;
            } // switch ( FstatMethod )

          if ( returnAtoms ) {
            multiFstatAtoms->data[X] = FstatAtoms;     // copy pointer to IFO-specific Fstat-atoms 'contents'
          }

          XLAL_CHECK ( isfinite(creal(FaX)) && isfinite(cimag(FaX)) && isfinite(creal(FbX)) && isfinite(cimag(FbX)), XLAL_EFPOVRFLW );

          // compute single-IFO F-stats, if requested
          if ( returnSingleF )
            {
              REAL8 AdX = multiAMcoef->data[X]->A;
              REAL8 BdX = multiAMcoef->data[X]->B;
              REAL8 CdX = multiAMcoef->data[X]->C;
              REAL8 DdX_inv = 1.0 / multiAMcoef->data[X]->D;
              REAL8 EdX = 0;

              /* compute final single-IFO F-stat */
              FX[X] = XLALComputeFstatFromFaFb ( FaX, FbX, AdX, BdX, CdX, EdX, DdX_inv );

            } /* if returnSingleF */

          /* Fa = sum_X Fa_X */
          Fa += FaX;

          /* Fb = sum_X Fb_X */
          Fb += FbX;

        } // for  X < numDetectors

      // ----- compute final Fstatistic-value -----
      /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
       * therefore there is a factor of 2 difference with respect to the equations in JKS, which
       * where based on the single-sided PSD.
       */
      REAL8 Ad = multiAMcoef->Mmunu.Ad;
      REAL8 Bd = multiAMcoef->Mmunu.Bd;
      REAL8 Cd = multiAMcoef->Mmunu.Cd;
      REAL8 Dd_inv = 1.0 / multiAMcoef->Mmunu.Dd;
      REAL8 Ed = 0;

      F = XLALComputeFstatFromFaFb ( Fa, Fb, Ad, Bd, Cd, Ed, Dd_inv );

      // this needs to be free'ed, as it's currently not buffered
      XLALDestroyMultiSSBtimes ( multiBinary );
      // --------------------------------------------------

      // Return multi-detector 2F
      if (whatToCompute & FSTATQ_2F) {
        Fstats->twoF[k] = 2.0 * F;   // *** Return value of 2F ***
      }

      // Return multi-detector Fa & Fb
      if (whatToCompute & FSTATQ_FAFB) {
        Fstats->FaFb[k].Fa = Fa;
        Fstats->FaFb[k].Fb = Fb;
      }

      // Return 2F per detector
      if (whatToCompute & FSTATQ_2F_PER_DET) {
        for (UINT4 X = 0; X < Fstats->numDetectors; ++X) {
          Fstats->twoFPerDet[X][k] = 2.0 * FX[X];   // *** Return value of 2F ***
        }
      }

      // Return Fa & Fb per detector
      if (whatToCompute & FSTATQ_FAFB_PER_DET) {
        XLAL_ERROR(XLAL_EFAILED, "Unimplemented!");
      }

      // Return F-atoms per detector
      if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
        XLALDestroyMultiFstatAtomVector(Fstats->multiFatoms[k]);
        Fstats->multiFatoms[k] = multiFstatAtoms;
      }

    } // for k < Fstats->numFreqBins

  // Return amplitude modulation coefficients
  Fstats->Mmunu = demod->prevMultiAMcoef->Mmunu;

  return XLAL_SUCCESS;

} // ComputeFstat_Demod()


static void
DestroyFstatInput_Demod ( FstatInput_Demod* demod )
{

  XLALDestroyMultiSFTVector ( demod->multiSFTs);
  XLALDestroyMultiSSBtimes  ( demod->prevMultiSSBtimes );
  XLALDestroyMultiAMCoeffs  ( demod->prevMultiAMcoef );
  XLALFree ( demod );

} // DestroyFstatInput_Demod()

FstatInput*
XLALCreateFstatInput_Demod( const UINT4 Dterms,
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

} // XLALCreateFstatInput_Demod()

static int
SetupFstatInput_Demod ( FstatInput_Demod *demod,
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

  return XLAL_SUCCESS;

} // SetupFstatInput_Demod()


static int
GetFstatExtraBins_Demod ( FstatInput_Demod* demod )
{
  // Check input
  XLAL_CHECK(demod != NULL, XLAL_EFAULT);

  // Demodulation requires 'Dterms' extra frequency bins
  return demod->Dterms;

} // GetFstatExtraBins_Demod()
