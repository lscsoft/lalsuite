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

// ---------- BEGIN: Demod-specific timing model data ----------
typedef struct tagFstatTimingDemod
{
  REAL4 Nsft;			   //  (average) number of SFTs per detector
  REAL4 tau0_coreLD;		   // 'fundamental' timing coefficient: fully-buffered time to compute F-stat for one SFT
  REAL4 tau0_bufferLD;		   // 'fundamental' timing coefficient: time to re-compute buffered quantities for one SFT
} FstatTimingDemod;

static char FstatTimingDemodHelp[] =
  "%%%% ----- Demod-specific F-statistic timing model -----\n"
  "%%%% Nsft:          (average) number of SFTs per detector\n"
  "%%%% tau0_coreLD:    timing coefficient for core Demod F-stat time\n"
  "%%%% tau0_bufferLD:  timing coefficient for computation of buffered quantities\n"
  "%%%%\n"
  "%%%% Demod F-statistic timing model:\n"
  "%%%% tauF_core      =  Nsft * tau0_coreLD\n"
  "%%%% tauF_buffer    =  Nsft * tau0_bufferLD / NFbin\n"
  "%%%%"
  "";
// ---------- END: Demod-specific timing model data ----------


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

  // ----- timing -----
  BOOLEAN collectTiming;			// flag whether or not to collect timing information
  FstatTimingGeneric timingGeneric;		// measured generic F-statistic timing values
  FstatTimingDemod   timingDemod;		// storage for measured Demod-specific timing model

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
  BOOLEAN collectTiming = demod->collectTiming;
  BOOLEAN BufferRecomputed = 0;	// keep track if buffer was recomputed in this call or not
  REAL8 Tau_buffer = 0;
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

  // initialize timing info struct
  if ( collectTiming ) {
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
      BufferRecomputed = 0;
    }
  else
    { // if not, compute SSB + AMcoef for this skyposition
      BufferRecomputed = 1;
      demod->timingGeneric.NBufferMisses ++;

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

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    Tau_buffer = (toc - tic);
  }

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

  if ( collectTiming )
    {
      toc = XLALGetCPUTime();
      REAL8 Tau_total = (toc - tic);

      FstatTimingGeneric *tiGen = &(demod->timingGeneric);	// generic part
      FstatTimingDemod   *tiLD  = &(demod->timingDemod);	// Demod-method specific part

      XLAL_CHECK ( numDetectors == tiGen->Ndet, XLAL_EINVAL, "Inconsistent number of detectors between XLALCreateSetup() [%d] and XLALComputeFstat() [%d]\n", tiGen->Ndet, numDetectors );

      // rescale all relevant timings to per-detector
      Tau_total  /= numDetectors;
      Tau_buffer /= numDetectors;

      // compute generic F-stat timing model contributions
      UINT4 NFbin     = Fstats->numFreqBins;
      REAL8 tauF_eff  = Tau_total / NFbin;
      REAL8 tauF_core = (Tau_total - Tau_buffer) / NFbin;

      // compute Demod timing model coefficients
      REAL8 NsftPerDet    = tiLD->Nsft;
      REAL8 tau0_coreLD   = tauF_core / NsftPerDet;

      // update the averaged timing-model quantities
      tiGen->NCalls ++;	// keep track of number of Fstat-calls for timing
#define updateAvgF(q) tiGen->q = ((tiGen->q *(tiGen->NCalls-1) + q)/(tiGen->NCalls))
      updateAvgF(tauF_eff);
      updateAvgF(tauF_core);
      // we also average NFbin, which can be different betnween different calls to XLALComputeFstat() (contrary to Ndet)
      updateAvgF(NFbin);

#define updateAvgLD(q) tiLD->q = ((tiLD->q *(tiGen->NCalls-1) + q)/(tiGen->NCalls))
      updateAvgLD(tau0_coreLD);

      // buffer-quantities only updated if buffer was actually recomputed
      if ( BufferRecomputed )
        {
          REAL8 tau0_bufferLD = Tau_buffer / NsftPerDet;
          REAL8 tauF_buffer   = Tau_buffer / NFbin;

          updateAvgF(tauF_buffer);
          updateAvgLD(tau0_bufferLD);
        } // if BufferRecomputed

    } // if collect timing

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

  // turn on timing collection if requested
  demod->collectTiming = optArgs->collectTiming;

  // initialize struct for collecting timing data, store invariant 'meta' quantities about this setup
  if ( demod->collectTiming )
    {
      XLAL_INIT_MEM(demod->timingGeneric);
      XLAL_INIT_MEM(demod->timingDemod);
      UINT4 numDetectors = demod->multiSFTs->length;
      UINT4 numSFTs = 0;
      for ( UINT4 X = 0; X < numDetectors; X ++ ) {
        numSFTs += demod->multiSFTs->data[X]->length;
      }
      demod->timingGeneric.Ndet = numDetectors;
      demod->timingDemod.Nsft	= 1.0 * numSFTs / numDetectors;	// average number of sfts *per detector*
    } // if collectTiming

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


int
XLALGetFstatTiming_Demod ( const void *method_data, FstatTimingGeneric *timingGeneric, FstatTimingModel *timingModel )
{
  XLAL_CHECK ( method_data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( timingGeneric != NULL, XLAL_EINVAL );
  XLAL_CHECK ( timingModel != NULL, XLAL_EINVAL );

  const DemodMethodData *demod = (const DemodMethodData*) method_data;
  XLAL_CHECK ( demod != NULL, XLAL_EINVAL );

  (*timingGeneric) = demod->timingGeneric; // struct-copy generic timing measurements

  const FstatTimingDemod *tiLD = &(demod->timingDemod);

  // return method-specific timing model values
  XLAL_INIT_MEM( (*timingModel) );

  UINT4 i = 0;
  timingModel->names[i]     = "Nsft";
  timingModel->values[i]    = tiLD->Nsft;

  i++;
  timingModel->names[i]     = "tau0_coreLD";
  timingModel->values[i]    = tiLD->tau0_coreLD;

  i++;
  timingModel->names[i]     = "tau0_bufferLD";
  timingModel->values[i]    = tiLD->tau0_bufferLD;

  timingModel->numVariables = i+1;
  timingModel->help         = FstatTimingDemodHelp;

  return XLAL_SUCCESS;
} // XLALGetFstatTiming_Demod()

// Generates an FstatInput Timeslice based on minStartGPS and maxStartGPS according to XLALCWGPSinRange().
void *
XLALFstatInputTimeslice_Demod ( const void *method_data,
                                const UINT4 iStart[PULSAR_MAX_DETECTORS],
                                const UINT4 iEnd[PULSAR_MAX_DETECTORS]
                                )
{
  XLAL_CHECK_NULL ( method_data != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( iStart != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( iEnd != NULL, XLAL_EINVAL );

  const DemodMethodData *demod_input = (const DemodMethodData *)method_data;

  // allocate memory and copy the input method_data struct
  DemodMethodData *demod_slice;
  XLAL_CHECK_NULL ( ( demod_slice = XLALCalloc ( 1, sizeof(*demod_input) ) ) != NULL, XLAL_ENOMEM );
  memcpy ( demod_slice, demod_input, sizeof(*demod_input) );

  UINT4 numIFOs = demod_input->multiSFTs->length;

  // prepare MultiSFTs struct
  MultiSFTVector *multiSFTs;
  XLAL_CHECK_NULL ( ( multiSFTs = XLALCalloc( 1, sizeof(*multiSFTs) ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( ( multiSFTs->data = XLALCalloc( numIFOs, sizeof(*multiSFTs->data) ) ) != NULL, XLAL_ENOMEM );
  multiSFTs->length = numIFOs;

  // loop over all detectors
  for ( UINT4 X = 0; X < numIFOs; X ++ )
    {
      XLAL_CHECK_NULL ( ( multiSFTs->data[X] = XLALCalloc ( 1, sizeof(*multiSFTs->data[X]) ) ) !=NULL, XLAL_EFUNC );

      if ( iStart[X] > iEnd[X] ) {
        continue; 	// empty slice
      }
      UINT4 sliceLength =  iEnd[X] - iStart[X] + 1; // guaranteed >= 1
      multiSFTs->data[X]->length = sliceLength;
      multiSFTs->data[X]->data   = &(demod_input->multiSFTs->data[X]->data[iStart[X]]);

    } // for loop over all detectors

  demod_slice->multiSFTs = multiSFTs;

  // empty all buffering quantities
  demod_slice->prevAlpha = 0;
  demod_slice->prevDelta = 0;
  XLAL_INIT_MEM(demod_slice->prevRefTime);
  demod_slice->prevMultiSSBtimes = NULL;
  demod_slice->prevMultiAMcoef = NULL;

  // reset timing counters
  XLAL_INIT_MEM(demod_slice->timingGeneric);
  XLAL_INIT_MEM(demod_slice->timingDemod);

  return demod_slice;

} // XLALFstatInputTimeslice_Demod()

///
/// Free all memory not needed by the orginal FstatInput structure
///
void
XLALDestroyFstatInputTimeslice_Demod ( void *method_data )
{
  if ( !method_data ) {
    return;
  }

  DemodMethodData *demod = (DemodMethodData*) method_data;

  XLALDestroyMultiSSBtimes  ( demod->prevMultiSSBtimes );
  XLALDestroyMultiAMCoeffs  ( demod->prevMultiAMcoef );

  for ( UINT4 X=0; X < demod->multiSFTs->length; X ++ ) {
    XLALFree ( demod->multiSFTs->data[X] );
  }

  XLALFree ( demod->multiSFTs->data );
  XLALFree ( demod->multiSFTs );
  XLALFree ( demod );

  return;

} // XLALDestroyFstatInputTimeslice_Demod()
