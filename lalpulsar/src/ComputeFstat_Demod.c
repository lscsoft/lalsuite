//
// Copyright (C) 2012--2015 Karl Wette
// Copyright (C) 2005--2007, 2009, 2010, 2012, 2014 Reinhard Prix
// Copyright (C) 2007--2010, 2012 Bernd Machenschalk
// Copyright (C) 2007 Chris Messenger
// Copyright (C) 2006 John T. Whelan, Badri Krishnan
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ComputeFstat_internal.h"

#include <lal/Factorial.h>
#include <lal/LogPrintf.h>
#include <lal/SinCosLUT.h>

// ========== Demod internals ==========

// ----- local types ----------
typedef struct tagDemodTimingInfo
{ // NOTE: all times refer to a single-detector timing case
  BOOLEAN collectTiming;	// turn on/off the collection of F-stat-method-specific timing-data (stored in method-data 'demod')

  UINT4 numFreqBins;		// number of frequency bins to compute F-stat for
  UINT4 numSFTs;		// total number of SFTs used
  UINT4 numDetectors;		// number of detectors

  REAL8 tauTotal;		// total time spent in XLALComputeFstatDemod()
  REAL8 tauBary;		// time spent in barycentring
  REAL8 tauF1Buf;		// Demod timing 'constant': Fstat time per template per detector for a 'buffered' case (same skypos, same numFreqBins)
  REAL8 tauF1NoBuf;		// Demod timing 'constant': Fstat time per template per detector for an 'unbuffered' usage (different skypos and numFreqBins)
} DemodTimingInfo;

typedef struct {
  int (*computefafb_func) (			// XLALComputeFaFb_...() function for the selected Demod hotloop
    COMPLEX8 *, COMPLEX8 *, FstatAtomVector **, const SFTVector *, const PulsarSpins, const SSBtimes *, const AMCoeffs *, const UINT4 Dterms
    );
  UINT4 Dterms;					// Number of terms to keep in Dirichlet kernel
  MultiSFTVector *multiSFTs;			// Input multi-detector SFTs
  REAL8 prevAlpha, prevDelta;			// buffering: previous skyposition computed
  LIGOTimeGPS prevRefTime;			// buffering: keep track of previous refTime for SSBtimes buffering
  MultiSSBtimes *prevMultiSSBtimes;		// buffering: previous multiSSB times, unique to skypos + SFTs
  MultiAMCoeffs *prevMultiAMcoef;		// buffering: previous AM-coeffs, unique to skypos + SFTs

  DemodTimingInfo timingInfo;			// temporary storage for collecting timing data
} DemodMethodData;

// ----- local prototypes ----------

int XLALSetupFstatDemod  ( void **method_data, FstatCommon *common, FstatMethodFuncs* funcs, MultiSFTVector *multiSFTs, const FstatOptionalArgs *optArgs );

int XLALComputeFaFb_Generic ( COMPLEX8 *Fa, COMPLEX8 *Fb, FstatAtomVector **FstatAtoms, const SFTVector *sfts,
                              const PulsarSpins fkdot, const SSBtimes *tSSB, const AMCoeffs *amcoe, const UINT4 Dterms );

int XLALComputeFaFb_OptC    ( COMPLEX8 *Fa, COMPLEX8 *Fb, FstatAtomVector **FstatAtoms, const SFTVector *sfts,
                              const PulsarSpins fkdot, const SSBtimes *tSSB, const AMCoeffs *amcoe, const UINT4 Dterms );

#ifdef HAVE_ALTIVEC
int XLALComputeFaFb_Altivec ( COMPLEX8 *Fa, COMPLEX8 *Fb, FstatAtomVector **FstatAtoms, const SFTVector *sfts,
                              const PulsarSpins fkdot, const SSBtimes *tSSB, const AMCoeffs *amcoe, const UINT4 Dterms );
#endif

#ifdef HAVE_SSE_COMPILER
int XLALComputeFaFb_SSE     ( COMPLEX8 *Fa, COMPLEX8 *Fb, FstatAtomVector **FstatAtoms, const SFTVector *sfts,
                              const PulsarSpins fkdot, const SSBtimes *tSSB, const AMCoeffs *amcoe, const UINT4 Dterms );
#endif

// ----- local function definitions ----------
static int
XLALComputeFstatDemod ( FstatResults* Fstats,
                        const FstatCommon *common,
                        void *method_data
                      )
{
  // Check input
  XLAL_CHECK(Fstats != NULL, XLAL_EFAULT);
  XLAL_CHECK(common != NULL, XLAL_EFAULT);
  XLAL_CHECK(method_data != NULL, XLAL_EFAULT);

  DemodMethodData *demod = (DemodMethodData*) method_data;
  // get internal timing info
  DemodTimingInfo *ti = &(demod->timingInfo);
  REAL8 tic = 0, toc = 0;

  // Get which F-statistic quantities to compute
  const FstatQuantities whatToCompute = Fstats->whatWasComputed;

  // handy shortcuts
  BOOLEAN returnAtoms = (whatToCompute & FSTATQ_ATOMS_PER_DET);
  PulsarDopplerParams thisPoint = Fstats->doppler;
  const REAL8 fStart = thisPoint.fkdot[0];
  const MultiSFTVector *multiSFTs = demod->multiSFTs;
  const MultiNoiseWeights *multiWeights = common->multiNoiseWeights;
  const MultiDetectorStateSeries *multiDetStates = common->multiDetectorStates;

  UINT4 numDetectors = multiSFTs->length;
  XLAL_CHECK ( multiDetStates->length == numDetectors, XLAL_EINVAL );
  XLAL_CHECK ( multiWeights==NULL || (multiWeights->length == numDetectors), XLAL_EINVAL );
  UINT4 numSFTs = 0;
  for ( UINT4 X = 0; X < numDetectors; X ++ ) {
    numSFTs += multiDetStates->data[X]->length;
  }

  // initialize timing info struct
  if ( ti->collectTiming )
    {
      XLAL_INIT_MEM ( (*ti) );
      ti->collectTiming = 1;

      ti->numDetectors = numDetectors;
      ti->numFreqBins = Fstats->numFreqBins;
      ti->numSFTs = numSFTs;

      tic = XLALGetCPUTime();
    }

  MultiSSBtimes *multiSSB = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  // ----- check if we have buffered SSB+AMcoef for current sky-position
  if ( (demod->prevAlpha == thisPoint.Alpha) && (demod->prevDelta == thisPoint.Delta ) &&
       (demod->prevMultiSSBtimes != NULL) && ( XLALGPSDiff(&demod->prevRefTime, &thisPoint.refTime) == 0 ) &&	// have SSB times for same reftime?
       (demod->prevMultiAMcoef != NULL)
       )
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
      demod->prevRefTime = thisPoint.refTime;
      XLALDestroyMultiAMCoeffs ( demod->prevMultiAMcoef );
      demod->prevMultiAMcoef = multiAMcoef;
      demod->prevAlpha = thisPoint.Alpha;
      demod->prevDelta = thisPoint.Delta;
    } // if could not reuse previously buffered quantites

  MultiSSBtimes *multiBinary = NULL;
  MultiSSBtimes *multiSSBTotal = NULL;
  // handle binary-orbital timing corrections, if applicable
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

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->tauBary = (toc - tic);
  }

  // ----- compute final Fstatistic-value -----
  REAL4 Ad = multiAMcoef->Mmunu.Ad;
  REAL4 Bd = multiAMcoef->Mmunu.Bd;
  REAL4 Cd = multiAMcoef->Mmunu.Cd;
  REAL4 Ed = multiAMcoef->Mmunu.Ed;;
  REAL4 Dd_inv = 1.0 / multiAMcoef->Mmunu.Dd;

  // ---------- Compute F-stat for each frequency bin ----------
  for ( UINT4 k = 0; k < Fstats->numFreqBins; k++ )
    {
      // Set frequency to search at
      thisPoint.fkdot[0] = fStart + k * Fstats->dFreq;

      COMPLEX8 Fa = 0;       		// complex amplitude Fa
      COMPLEX8 Fb = 0;                 // complex amplitude Fb
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
          COMPLEX8 FaX, FbX;
          FstatAtomVector *FstatAtoms = NULL;
          FstatAtomVector **FstatAtoms_p = returnAtoms ? (&FstatAtoms) : NULL;

          // call XLALComputeFaFb_...() function for the user-requested hotloop variant
          XLAL_CHECK ( (demod->computefafb_func) ( &FaX, &FbX, FstatAtoms_p, multiSFTs->data[X], thisPoint.fkdot,
                                                   multiSSBTotal->data[X], multiAMcoef->data[X], demod->Dterms ) == XLAL_SUCCESS, XLAL_EFUNC );

          if ( returnAtoms ) {
            multiFstatAtoms->data[X] = FstatAtoms;     // copy pointer to IFO-specific Fstat-atoms 'contents'
          }

          XLAL_CHECK ( isfinite(creal(FaX)) && isfinite(cimag(FaX)) && isfinite(creal(FbX)) && isfinite(cimag(FbX)), XLAL_EFPOVRFLW );

          if ( whatToCompute & FSTATQ_FAFB_PER_DET )
            {
              Fstats->FaPerDet[X][k] = FaX;
              Fstats->FbPerDet[X][k] = FbX;
            }

          // compute single-IFO F-stats, if requested
          if ( whatToCompute & FSTATQ_2F_PER_DET )
            {
              REAL4 AdX = multiAMcoef->data[X]->A;
              REAL4 BdX = multiAMcoef->data[X]->B;
              REAL4 CdX = multiAMcoef->data[X]->C;
              REAL4 EdX = 0;
              REAL4 DdX_inv = 1.0 / multiAMcoef->data[X]->D;

              // compute final single-IFO F-stat
              Fstats->twoFPerDet[X][k] = XLALComputeFstatFromFaFb ( FaX, FbX, AdX, BdX, CdX, EdX, DdX_inv );

            } // if FSTATQ_2F_PER_DET

          /* Fa = sum_X Fa_X */
          Fa += FaX;

          /* Fb = sum_X Fb_X */
          Fb += FbX;

        } // for  X < numDetectors

      if ( whatToCompute & FSTATQ_2F )
        {
          Fstats->twoF[k] = XLALComputeFstatFromFaFb ( Fa, Fb, Ad, Bd, Cd, Ed, Dd_inv );
        }

      // Return multi-detector Fa & Fb
      if ( whatToCompute & FSTATQ_FAFB )
        {
          Fstats->Fa[k] = Fa;
          Fstats->Fb[k] = Fb;
        }

      // Return F-atoms per detector
      if ( whatToCompute & FSTATQ_ATOMS_PER_DET )
        {
          XLALDestroyMultiFstatAtomVector ( Fstats->multiFatoms[k] );
          Fstats->multiFatoms[k] = multiFstatAtoms;
        }

    } // for k < Fstats->numFreqBins

  // this needs to be free'ed, as it's currently not buffered
  XLALDestroyMultiSSBtimes ( multiBinary );

  // Return amplitude modulation coefficients
  Fstats->Mmunu = demod->prevMultiAMcoef->Mmunu;

  // return per-detector antenna-pattern matrices
  for ( UINT4 X=0; X < numDetectors; X ++ )
    {
      Fstats->MmunuX[X].Ad = multiAMcoef->data[X]->A;
      Fstats->MmunuX[X].Bd = multiAMcoef->data[X]->B;
      Fstats->MmunuX[X].Cd = multiAMcoef->data[X]->C;
      Fstats->MmunuX[X].Dd = multiAMcoef->data[X]->D;
      Fstats->MmunuX[X].Ed = 0;
    }

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->tauTotal = (toc - tic);
    ti->tauF1NoBuf = ti->tauTotal / ( Fstats->numFreqBins * numDetectors );
    ti->tauF1Buf   = (ti->tauTotal - ti->tauBary) / ( Fstats->numFreqBins * numDetectors );
  }

  return XLAL_SUCCESS;

} // XLALComputeFstatDemod()


static void
XLALDestroyDemodMethodData ( void* method_data )
{

  DemodMethodData *demod = (DemodMethodData*) method_data;

  XLALDestroyMultiSFTVector ( demod->multiSFTs);
  XLALDestroyMultiSSBtimes  ( demod->prevMultiSSBtimes );
  XLALDestroyMultiAMCoeffs  ( demod->prevMultiAMcoef );
  XLALFree ( demod );

} // XLALDestroyDemodMethodData()

int
XLALSetupFstatDemod ( void **method_data,
                      FstatCommon *common,
                      FstatMethodFuncs* funcs,
                      MultiSFTVector *multiSFTs,
                      const FstatOptionalArgs *optArgs
                    )
{
  // Check input
  XLAL_CHECK ( method_data != NULL, XLAL_EFAULT );
  XLAL_CHECK ( common != NULL, XLAL_EFAULT );
  XLAL_CHECK ( funcs != NULL, XLAL_EFAULT );
  XLAL_CHECK ( multiSFTs != NULL, XLAL_EFAULT );
  XLAL_CHECK ( optArgs != NULL, XLAL_EFAULT );

  // Allocate method data
  DemodMethodData *demod = *method_data = XLALCalloc( 1, sizeof(*demod) );
  XLAL_CHECK( demod != NULL, XLAL_ENOMEM );

  // Set method function pointers
  funcs->compute_func = XLALComputeFstatDemod;
  funcs->method_data_destroy_func = XLALDestroyDemodMethodData;
  funcs->workspace_destroy_func = NULL;

  // Save pointer to SFTs
  demod->multiSFTs = multiSFTs;

  // Save Dterms
  demod->Dterms = optArgs->Dterms;

  // Save flag about collecting timing data
  demod->timingInfo.collectTiming = optArgs->collectTiming;

  // Select XLALComputeFaFb_...() function for the user-requested hotloop variant
  switch ( optArgs->FstatMethod ) {
  case  FMETHOD_DEMOD_GENERIC:
    demod->computefafb_func = XLALComputeFaFb_Generic;
    break;
  case FMETHOD_DEMOD_OPTC:
    demod->computefafb_func = XLALComputeFaFb_OptC;
    break;
#ifdef HAVE_ALTIVEC
  case FMETHOD_DEMOD_ALTIVEC:
    demod->computefafb_func = XLALComputeFaFb_Altivec;
    break;
#endif
#ifdef HAVE_SSE_COMPILER
  case FMETHOD_DEMOD_SSE:
    demod->computefafb_func = XLALComputeFaFb_SSE;
    break;
#endif
  default:
    XLAL_ERROR ( XLAL_EINVAL, "Invalid Demod hotloop optArgs->FstatMethod='%d'", optArgs->FstatMethod );
    break;
  }

  return XLAL_SUCCESS;

} // XLALSetupFstatDemod()


// export timing constants from internally-stored 'method data' struct
int
XLALGetFstatTiming_Demod ( const void* method_data, REAL8 *tauF1Buf, REAL8 *tauF1NoBuf )
{
  XLAL_CHECK ( method_data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( (tauF1Buf != NULL) && (tauF1NoBuf != NULL), XLAL_EINVAL );

  const DemodMethodData *demod = (const DemodMethodData *)method_data;

  (*tauF1Buf) = demod->timingInfo.tauF1Buf;
  (*tauF1NoBuf) = demod->timingInfo.tauF1NoBuf;

  return XLAL_SUCCESS;

} // XLALGetFstatTiming_Demod()

// append detailed (method-specific) timing info to given file
int
AppendFstatTimingInfo2File_Demod ( const void* method_data, FILE *fp, BOOLEAN printHeader )
{
  XLAL_CHECK ( method_data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( fp != NULL, XLAL_EINVAL );

  const DemodMethodData *demod = (const DemodMethodData *)method_data;
  if ( printHeader ) {
    fprintf (fp, "%%%%%8s %4s %10s ", "Nfreq", "Ndet", "Nsft" );
    fprintf (fp, "%10s %10s %10s %10s %10s\n",
             "tauTotal", "tauBary", "tauF1NoBuf", "tauF1Buf", "tauF0Demod" );
  }

  const DemodTimingInfo *ti = &(demod->timingInfo);
  fprintf (fp, "%10d %4d %10d", ti->numFreqBins, ti->numDetectors, ti->numSFTs );
  REAL8 numSFTsPerDet = ti->numSFTs / ti->numDetectors;
  fprintf (fp, "%10.1e %10.1e %10.1e %10.1e %10.1e\n", ti->tauTotal, ti->tauBary, ti->tauF1NoBuf, ti->tauF1Buf, ti->tauF1Buf / numSFTsPerDet );

  return XLAL_SUCCESS;
} // AppendFstatTimingInfo2File_Demod()
