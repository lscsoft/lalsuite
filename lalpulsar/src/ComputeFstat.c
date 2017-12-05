//
// Copyright (C) 2012--2015 Karl Wette
// Copyright (C) 2005--2007, 2010, 2014 Reinhard Prix
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

#define _COMPUTE_FSTAT_C
#include "ComputeFstat_internal.h"

#include <lal/LALString.h>
#include <lal/LALSIMD.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/VectorMath.h>

// ---------- Internal struct definitions ---------- //

// Internal definition of input data structure
struct tagFstatInput {
  int singleFreqBin;					// True if XLALComputeFstat() can only compute a single frequency bin, due to zero dFreq being passed to XLALCreateFstatInput()
  FstatMethodType method;				// Method to use for computing the F-statistic
  FstatCommon common;					// Common input data
  int *workspace_refcount;				// Reference counter for the shared workspace 'common.workspace'
  FstatMethodFuncs method_funcs;			// Function pointers for F-statistic method
  void *method_data;					// F-statistic method data
};

// ---------- Internal prototypes ---------- //

static int XLALSelectBestFstatMethod ( FstatMethodType *method );

int XLALSetupFstatDemod  ( void **method_data, FstatCommon *common, FstatMethodFuncs* funcs, MultiSFTVector *multiSFTs, const FstatOptionalArgs *optArgs );
int XLALSetupFstatResamp ( void **method_data, FstatCommon *common, FstatMethodFuncs* funcs, MultiSFTVector *multiSFTs, const FstatOptionalArgs *optArgs );

// ---------- Constant variable definitions ---------- //

static const char *const FstatMethodNames[FMETHOD_END] = {
  [FMETHOD_DEMOD_GENERIC]	= "DemodGeneric",
  [FMETHOD_DEMOD_OPTC]		= "DemodOptC",
  [FMETHOD_DEMOD_ALTIVEC]	= "DemodAltivec",
  [FMETHOD_DEMOD_SSE]		= "DemodSSE",
  [FMETHOD_DEMOD_BEST]		= "DemodBest",

  [FMETHOD_RESAMP_GENERIC]	= "ResampGeneric",
  [FMETHOD_RESAMP_BEST]		= "ResampBest",
};


const FstatOptionalArgs FstatOptionalArgsDefaults = {
  .randSeed = 0,
  .SSBprec = SSBPREC_RELATIVISTICOPT,
  .Dterms = 8,
  .runningMedianWindow = 50,
  .FstatMethod = FMETHOD_DEMOD_BEST,
  .injectSources = NULL,
  .injectSqrtSX = NULL,
  .assumeSqrtSX = NULL,
  .prevInput = NULL,
  .collectTiming = 0,
  .resampFFTPowerOf2 = 1
};

// ==================== Function definitions =================== //

///
/// Create a #FstatInputVector of the given length, for example for setting up
/// F-stat searches over several segments.
///
FstatInputVector*
XLALCreateFstatInputVector ( const UINT4 length            ///< [in] Length of the #FstatInputVector.
                             )
{
  // Allocate and initialise vector container
  FstatInputVector* inputs;
  XLAL_CHECK_NULL ( (inputs = XLALCalloc ( 1, sizeof(*inputs))) != NULL, XLAL_ENOMEM );
  inputs->length = length;

  // Allocate and initialise vector data
  if (inputs->length > 0) {
    XLAL_CHECK_NULL ( (inputs->data = XLALCalloc ( inputs->length, sizeof(inputs->data[0]) )) != NULL, XLAL_ENOMEM );
  }

  return inputs;

} // XLALCreateFstatInputVector()

///
/// Free all memory associated with a #FstatInputVector structure.
///
void
XLALDestroyFstatInputVector ( FstatInputVector* inputs        ///< [in] #FstatInputVector structure to be freed.
                              )
{
  if ( inputs == NULL ) {
    return;
  }

  if ( inputs->data )
    {
      for ( UINT4 i = 0; i < inputs->length; ++i ) {
        XLALDestroyFstatInput ( inputs->data[i] );
      }
      XLALFree ( inputs->data );
    }

  XLALFree ( inputs );

  return;

} // XLALDestroyFstatInputVector()

///
/// Create a #FstatAtomVector of the given length.
///
FstatAtomVector*
XLALCreateFstatAtomVector ( const UINT4 length ///< [in] Length of the #FstatAtomVector.
                            )
{
  // Allocate and initialise vector container
  FstatAtomVector* atoms;
  XLAL_CHECK_NULL ( (atoms = XLALCalloc ( 1, sizeof(*atoms) )) != NULL, XLAL_ENOMEM );
  atoms->length = length;

  // Allocate and initialise vector data
  if (atoms->length > 0) {
    XLAL_CHECK_NULL ( (atoms->data = XLALCalloc (atoms->length, sizeof(atoms->data[0]) )) != NULL, XLAL_ENOMEM );
  }

  return atoms;

} // XLALCreateFstatAtomVector()

///
/// Free all memory associated with a #FstatAtomVector structure.
///
void
XLALDestroyFstatAtomVector ( FstatAtomVector *atoms      ///< [in] #FstatAtomVector structure to be freed.
                             )
{
  if ( atoms == NULL ) {
    return;
  }

  if ( atoms->data ) {
    XLALFree ( atoms->data );
  }
  XLALFree ( atoms );

  return;

} // XLALDestroyFstatAtomVector()

///
/// Create a #MultiFstatAtomVector of the given length.
///
MultiFstatAtomVector*
XLALCreateMultiFstatAtomVector ( const UINT4 length   ///< [in] Length of the #MultiFstatAtomVector.
                                 )
{
  // Allocate and initialise vector container
  MultiFstatAtomVector* multiAtoms;
  XLAL_CHECK_NULL ( (multiAtoms = XLALCalloc(1, sizeof(*multiAtoms))) != NULL, XLAL_ENOMEM );
  multiAtoms->length = length;

  // Allocate and initialise vector data
  if ( multiAtoms->length > 0 ) {
    XLAL_CHECK_NULL ( (multiAtoms->data = XLALCalloc ( multiAtoms->length, sizeof(multiAtoms->data[0]) )) != NULL, XLAL_ENOMEM );
  }

  return multiAtoms;

} // XLALCreateMultiFstatAtomVector()

///
/// Free all memory associated with a #MultiFstatAtomVector structure.
///
void
XLALDestroyMultiFstatAtomVector ( MultiFstatAtomVector *multiAtoms  ///< [in] #MultiFstatAtomVector structure to be freed.
                                  )
{
  if ( multiAtoms == NULL ) {
    return;
  }

  for ( UINT4 X = 0; X < multiAtoms->length; ++X ) {
    XLALDestroyFstatAtomVector ( multiAtoms->data[X] );
  }
  XLALFree ( multiAtoms->data );
  XLALFree ( multiAtoms );

  return;

} // XLALDestroyMultiFstatAtomVector()

///
/// Create a fully-setup \c FstatInput structure for computing the \f$\mathcal{F}\f$-statistic using XLALComputeFstat().
///
FstatInput *
XLALCreateFstatInput ( const SFTCatalog *SFTcatalog,              ///< [in] Catalog of SFTs to either load from files, or generate in memory.
                                                                  ///< The \c locator field of each ::SFTDescriptor must be \c !=NULL for SFT loading, and \c ==NULL for SFT generation.
                       const REAL8 minCoverFreq,                  ///< [in] Minimum instantaneous frequency which will be covered over the SFT time span.
                       const REAL8 maxCoverFreq,                  ///< [in] Maximum instantaneous frequency which will be covered over the SFT time span.
                       const REAL8 dFreq,                         ///< [in] Requested spacing of \f$\mathcal{F}\f$-statistic frequency bins. May be zero \e only for single-frequency searches.
                       const EphemerisData *ephemerides,          ///< [in] Ephemerides for the time-span of the SFTs.
                       const FstatOptionalArgs *optionalArgs      ///< [in] Optional 'advanced-level' and method-specific extra arguments; NULL: use defaults from FstatOptionalArgsDefaults.
                       )
{
  // Check catalog
  XLAL_CHECK_NULL ( SFTcatalog != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( SFTcatalog->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL ( SFTcatalog->data != NULL, XLAL_EINVAL );
  for ( UINT4 i = 1; i < SFTcatalog->length; ++i ) {
    XLAL_CHECK_NULL ( (SFTcatalog->data[0].locator == NULL) == (SFTcatalog->data[i].locator == NULL), XLAL_EINVAL,
                      "All 'locator' fields of SFTDescriptors in 'SFTcatalog' must be either NULL or !NULL." );
  }

  // Check remaining required parameters
  XLAL_CHECK_NULL ( isfinite(minCoverFreq) && ( minCoverFreq > 0 ) && isfinite(maxCoverFreq) && ( maxCoverFreq > 0 ), XLAL_EINVAL, "Check failed: minCoverFreq=%f and maxCoverFreq=%f must be finite and positive!", minCoverFreq, maxCoverFreq );
  XLAL_CHECK_NULL ( maxCoverFreq > minCoverFreq, XLAL_EINVAL, "Check failed: maxCoverFreq>minCoverFreq (%f<=%f)!", maxCoverFreq, minCoverFreq );
  XLAL_CHECK_NULL ( ephemerides != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( dFreq >= 0, XLAL_EINVAL);

  // Create local copy of optional arguments, or use defaults if not given
  FstatOptionalArgs optArgs;
  if ( optionalArgs != NULL ) {
    optArgs = *optionalArgs;
  } else {
    optArgs = FstatOptionalArgsDefaults;
  }

  // Check optional arguments sanity
  XLAL_CHECK_NULL ( (optArgs.injectSources == NULL) || ((optArgs.injectSources->length > 0) && (optArgs.injectSources->data != NULL)), XLAL_EINVAL );
  XLAL_CHECK_NULL ( (optArgs.injectSqrtSX == NULL) || (optArgs.injectSqrtSX->length > 0), XLAL_EINVAL );
  XLAL_CHECK_NULL ( (optArgs.assumeSqrtSX == NULL) || (optArgs.assumeSqrtSX->length > 0), XLAL_EINVAL );
  XLAL_CHECK_NULL ( optArgs.SSBprec < SSBPREC_LAST, XLAL_EINVAL );

  // Check optional Fstat method type argument
  XLAL_CHECK_NULL ( ( FMETHOD_START < optArgs.FstatMethod ) && ( optArgs.FstatMethod < FMETHOD_END ), XLAL_EINVAL );
  XLAL_CHECK_NULL ( FstatMethodNames[optArgs.FstatMethod] != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL ( XLALSelectBestFstatMethod( &optArgs.FstatMethod ) == XLAL_SUCCESS, XLAL_EFAULT );

  //
  // Parse which F-statistic method to use, and set these variables:
  // - extraBinsMethod:   any extra SFT frequency bins required by the method
  // - setupFuncMethod:   method setup function, called at end of XLALCreateFstatInput()
  //
  int extraBinsMethod = 0;
  int (*setupFuncMethod) ( void **, FstatCommon *, FstatMethodFuncs*, MultiSFTVector *, const FstatOptionalArgs * );
  switch (optArgs.FstatMethod) {

  case FMETHOD_DEMOD_GENERIC:		// Demod: generic C hotloop
    XLAL_CHECK_NULL ( optArgs.Dterms > 0, XLAL_EINVAL );
    extraBinsMethod = optArgs.Dterms;
    setupFuncMethod = XLALSetupFstatDemod;
    break;

  case FMETHOD_DEMOD_OPTC:		// Demod: gptimized C hotloop using Akos' algorithm
    XLAL_CHECK_NULL ( optArgs.Dterms <= 20, XLAL_EINVAL, "Selected Hotloop variant 'OptC' only works for Dterms <= 20, got %d\n", optArgs.Dterms );
    extraBinsMethod = optArgs.Dterms;
    setupFuncMethod = XLALSetupFstatDemod;
    break;

  case FMETHOD_DEMOD_ALTIVEC:		// Demod: Altivec hotloop variant
    XLAL_CHECK_NULL ( optArgs.Dterms == 8, XLAL_EINVAL, "Selected Hotloop variant 'Altivec' only works for Dterms == 8, got %d\n", optArgs.Dterms );
    extraBinsMethod = optArgs.Dterms;
    setupFuncMethod = XLALSetupFstatDemod;
    break;

  case FMETHOD_DEMOD_SSE:		// Demod: SSE hotloop with precalc divisors
    XLAL_CHECK_NULL ( optArgs.Dterms == 8, XLAL_EINVAL, "Selected Hotloop variant 'SSE' only works for Dterms == 8, got %d\n", optArgs.Dterms );
    extraBinsMethod = optArgs.Dterms;
    setupFuncMethod = XLALSetupFstatDemod;
    break;

  case FMETHOD_RESAMP_GENERIC:		// Resamp: generic implementation
    extraBinsMethod = 8;   // use 8 extra bins to give better agreement with Demod(w Dterms=8) near the boundaries
    setupFuncMethod = XLALSetupFstatResamp;
    break;

  default:
    XLAL_ERROR_NULL ( XLAL_EFAILED, "Missing switch case for optArgs.FstatMethod='%d'\n", optArgs.FstatMethod );
  }
  XLAL_CHECK_NULL ( extraBinsMethod >= 0, XLAL_EFAILED );
  XLAL_CHECK_NULL ( setupFuncMethod != NULL, XLAL_EFAILED );

  // Determine whether to load and/or generate SFTs
  const BOOLEAN loadSFTs = (SFTcatalog->data[0].locator != NULL);
  const BOOLEAN generateSFTs = (optArgs.injectSources != NULL) || (optArgs.injectSqrtSX != NULL);
  XLAL_CHECK_NULL ( loadSFTs || generateSFTs, XLAL_EINVAL, "Can neither load nor generate SFTs with given parameters" );

  // Create top-level input data struct
  FstatInput* input;
  XLAL_CHECK_NULL ( (input = XLALCalloc ( 1, sizeof(*input) )) != NULL, XLAL_ENOMEM );
  input->method = optArgs.FstatMethod;
  FstatCommon *common = &input->common;      // handy shortcut

  // Determine whether we can re-used workspace from a previous call to XLALCreateFstatInput()
  if ( optArgs.prevInput != NULL ) {

    // Check that F-stat method being used agrees with 'prevInput'
    XLAL_CHECK_NULL( optArgs.prevInput->method == input->method, XLAL_EFAILED, "Cannot use workspace from 'prevInput' with different FstatMethod '%d'!='%d'", optArgs.prevInput->method, input->method );

    // Get pointers to workspace and workspace reference counter in 'prevInput'
    common->workspace = optArgs.prevInput->common.workspace;
    input->workspace_refcount = optArgs.prevInput->workspace_refcount;

    // Increment reference count
    ++(*input->workspace_refcount);
    XLALPrintInfo( "%s: re-using workspace from 'optionalArgs.prevInput', reference count = %i\n", __func__, *input->workspace_refcount );

  } else {

    // Workspace must be allocated by method setup function
    common->workspace = NULL;

    // Allocate memory for reference counter; when reference count reaches 0, memory must be destroyed
    XLAL_CHECK_NULL ( ( input->workspace_refcount = XLALCalloc ( 1, sizeof(*input->workspace_refcount) ) ) != NULL, XLAL_ENOMEM );

    // Initialise reference counter to 1
    (*input->workspace_refcount) = 1;
    XLALPrintInfo( "%s: allocating new workspace, reference count = %i\n", __func__, *input->workspace_refcount );

  }

  // Determine the time baseline of an SFT
  const REAL8 Tsft = 1.0 / SFTcatalog->data[0].header.deltaF;

  // Compute the mid-time and time-span of the SFTs
  double Tspan = 0;
  {
    const LIGOTimeGPS startTime = SFTcatalog->data[0].header.epoch;
    const LIGOTimeGPS endTime   = SFTcatalog->data[SFTcatalog->length - 1].header.epoch;
    common->midTime = startTime;
    Tspan = Tsft + XLALGPSDiff( &endTime, &startTime );
    XLALGPSAdd ( &common->midTime, 0.5 * Tspan );
  }

  // Determine the frequency spacing: if dFreq==0, default to 1.0/Tspan and set singleFreqBin==true
  input->singleFreqBin = (dFreq == 0);
  common->dFreq = input->singleFreqBin ? 1.0/Tspan : dFreq;

  // Determine the frequency band required by each method 'minFreqMethod',
  // as well as the frequency band required to load or generate initial SFTs for 'minFreqFull'
  // the difference being that for noise-floor estimation, we need extra frequency-bands for the
  // running median
  REAL8 minFreqMethod, maxFreqMethod;
  REAL8 minFreqFull, maxFreqFull;
  {
    // Number of extra frequency bins required by: F-stat method, and running median
    int extraBinsFull = extraBinsMethod + optArgs.runningMedianWindow/2 + 1; // NOTE: running-median window needed irrespective of assumeSqrtSX!

    // Extend frequency range by number of extra bins times SFT bin width
    const REAL8 extraFreqMethod = extraBinsMethod / Tsft;
    minFreqMethod = minCoverFreq - extraFreqMethod;
    maxFreqMethod = maxCoverFreq + extraFreqMethod;

    const REAL8 extraFreqFull = extraBinsFull / Tsft;
    minFreqFull = minCoverFreq - extraFreqFull;
    maxFreqFull = maxCoverFreq + extraFreqFull;

  } // end: block to determine frequency-bins range

  // Load SFTs, if required, and extract detectors and timestamps
  MultiSFTVector *multiSFTs = NULL;
  if (loadSFTs)
    {
      // Load all SFTs at once
      XLAL_CHECK_NULL ( ( multiSFTs = XLALLoadMultiSFTs(SFTcatalog, minFreqFull, maxFreqFull) ) != NULL, XLAL_EFUNC );

      // Extract detectors and timestamps from SFTs
      XLAL_CHECK_NULL ( XLALMultiLALDetectorFromMultiSFTs ( &common->detectors, multiSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL ( ( common->multiTimestamps = XLALExtractMultiTimestampsFromSFTs ( multiSFTs ) ) != NULL,  XLAL_EFUNC );

    }
  else
    {
      // Create a multi-view of SFT catalog
      MultiSFTCatalogView *multiSFTcatalog;
      XLAL_CHECK_NULL ( (multiSFTcatalog = XLALGetMultiSFTCatalogView(SFTcatalog)) != NULL, XLAL_EFUNC );

      // Extract detectors and timestamps from multi-view of SFT catalog
      XLAL_CHECK_NULL ( XLALMultiLALDetectorFromMultiSFTCatalogView ( &common->detectors, multiSFTcatalog ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL ( ( common->multiTimestamps = XLALTimestampsFromMultiSFTCatalogView ( multiSFTcatalog ) ) != NULL,  XLAL_EFUNC );

      // Cleanup
      XLALDestroyMultiSFTCatalogView ( multiSFTcatalog );
    } // end: if !loadSFTs

  // Check length of multi-noise floor arrays
  XLAL_CHECK_NULL ( (optArgs.injectSqrtSX == NULL) || (optArgs.injectSqrtSX->length == common->detectors.length), XLAL_EINVAL );
  XLAL_CHECK_NULL ( (optArgs.assumeSqrtSX == NULL) || (optArgs.assumeSqrtSX->length == common->detectors.length), XLAL_EINVAL );

  // Generate SFTs with injections and noise, if required
  if (generateSFTs)
    {
      // Initialise parameters struct for XLALCWMakeFakeMultiData()
      CWMFDataParams XLAL_INIT_DECL(MFDparams);
      MFDparams.fMin = minFreqFull;
      MFDparams.Band = maxFreqFull - minFreqFull;
      MFDparams.multiIFO = common->detectors;
      MFDparams.multiTimestamps = *(common->multiTimestamps);
      MFDparams.randSeed = optArgs.randSeed;

      // Set noise floors if sqrtSX is given; otherwise noise floors are zero
      if ( optArgs.injectSqrtSX != NULL ) {
        MFDparams.multiNoiseFloor = *(optArgs.injectSqrtSX);
      } else {
        MFDparams.multiNoiseFloor.length = common->detectors.length;
      }

      // Generate SFTs with injections
      MultiSFTVector *fakeMultiSFTs = NULL;
      XLAL_CHECK_NULL ( XLALCWMakeFakeMultiData ( &fakeMultiSFTs, NULL, optArgs.injectSources, &MFDparams, ephemerides ) == XLAL_SUCCESS, XLAL_EFUNC );

      // If SFTs were loaded, add generated SFTs to then, otherwise just used generated SFTs
      if (multiSFTs != NULL) {
        XLAL_CHECK_NULL ( XLALMultiSFTVectorAdd ( multiSFTs, fakeMultiSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLALDestroyMultiSFTVector ( fakeMultiSFTs );
      } else {
        multiSFTs = fakeMultiSFTs;
      }

    } // if generateSFTs

  // Check that no single-SFT input vectors are given to avoid returning singular results
  for ( UINT4 X = 0; X < common->detectors.length; ++X ) {
    XLAL_CHECK_NULL ( multiSFTs->data[X]->length > 1, XLAL_EINVAL, "Need more than 1 SFTs per Detector!\n" );
  }

  // Normalise SFTs using either running median or assumed PSDs
  MultiPSDVector *runningMedian;
  XLAL_CHECK_NULL ( (runningMedian = XLALNormalizeMultiSFTVect ( multiSFTs, optArgs.runningMedianWindow, optArgs.assumeSqrtSX )) != NULL, XLAL_EFUNC );

  // Calculate SFT noise weights from PSD
  XLAL_CHECK_NULL ( (common->multiNoiseWeights = XLALComputeMultiNoiseWeights ( runningMedian, optArgs.runningMedianWindow, 0 )) != NULL, XLAL_EFUNC );

  // for use within this modul: remove the normalization from the noise-weights: the normalisation factor is saved separately
  // in the struct, and this allows us to use and extract subsets of the noise-weights without having to re-normalize them
  for ( UINT4 X = 0; X < common->detectors.length; ++X )
    {
      XLAL_CHECK_NULL ( XLALVectorScaleREAL8 ( common->multiNoiseWeights->data[X]->data, common->multiNoiseWeights->Sinv_Tsft, common->multiNoiseWeights->data[X]->data, common->multiNoiseWeights->data[X]->length) == XLAL_SUCCESS, XLAL_EFUNC);
    }
  common->multiNoiseWeights->isNotNormalized = (1 == 1);

  // at this point we're done with running-median noise estimation and can 'trim' the SFTs back to
  // the width actually required by the Fstat-methods *methods*.
  // NOTE: this is especially relevant for resampling, where the frequency-band determines the sampling
  // rate, and the number of samples that need to be FFT'ed
  XLAL_CHECK_NULL ( XLALMultiSFTVectorResizeBand ( multiSFTs, minFreqMethod, maxFreqMethod - minFreqMethod ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Get detector states, with a timestamp shift of Tsft/2
  const REAL8 tOffset = 0.5 * Tsft;
  XLAL_CHECK_NULL ( (common->multiDetectorStates = XLALGetMultiDetectorStates ( common->multiTimestamps, &common->detectors, ephemerides, tOffset )) != NULL, XLAL_EFUNC );

  // Save ephemerides and SSB precision
  common->ephemerides = ephemerides;
  common->SSBprec = optArgs.SSBprec;

  // Call the appropriate method function to setup their input data structures
  // - The method input data structures are expected to take ownership of the
  //   SFTs, which is why 'input->common' does not retain a pointer to them
  FstatMethodFuncs *funcs = &input->method_funcs;
  XLAL_CHECK_NULL( (setupFuncMethod) ( &input->method_data, common, funcs, multiSFTs, &optArgs ) == XLAL_SUCCESS, XLAL_EFUNC );

  // If setup function allocated a workspace, check that it also supplied a destructor function
  XLAL_CHECK_NULL( common->workspace == NULL || funcs->workspace_destroy_func != NULL, XLAL_EFAILED );

  // Cleanup
  XLALDestroyMultiPSDVector ( runningMedian );

  return input;

} // XLALCreateFstatInput()

///
/// Returns the human-readable name of the \f$\mathcal{F}\f$-statistic method being used by a \c FstatInput structure.
///
const CHAR *
XLALGetFstatInputMethodName ( const FstatInput* input    ///< [in] \c FstatInput structure.
                             )
{
  // Check input
  XLAL_CHECK_NULL ( input != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( FstatMethodNames[input->method] != NULL, XLAL_EFAULT );

  return FstatMethodNames[input->method];

} // XLALGetFstatMethodName()

///
/// Returns the detector information stored in a \c FstatInput structure.
///
const MultiLALDetector*
XLALGetFstatInputDetectors ( const FstatInput* input    ///< [in] \c FstatInput structure.
                             )
{
  // Check input
  XLAL_CHECK_NULL ( input != NULL, XLAL_EINVAL );

  return &input->common.detectors;

} // XLALGetFstatInputDetectors()

///
/// Returns the SFT timestamps stored in a \c FstatInput structure.
///
const MultiLIGOTimeGPSVector*
XLALGetFstatInputTimestamps ( const FstatInput* input   ///< [in] \c FstatInput structure.
                              )
{
  // Check input
  XLAL_CHECK_NULL ( input != NULL, XLAL_EINVAL );

  return input->common.multiTimestamps;

} // XLALGetFstatInputTimestamps()

///
/// Returns the multi-detector noise weights stored in a \c FstatInput structure.
/// Note: the returned object is allocated here and needs to be free'ed by caller.
///
MultiNoiseWeights*
XLALGetFstatInputNoiseWeights ( const FstatInput* input     ///< [in] \c FstatInput structure.
                                )
{
  // Check input
  XLAL_CHECK_NULL ( input != NULL, XLAL_EINVAL );

  // copy and rescale data to a new allocated MultiNoiseWeights structure
  UINT4 numIFOs = input->common.multiNoiseWeights->length;
  MultiNoiseWeights *ret = NULL;
  XLAL_CHECK_NULL ( ( ret = XLALCalloc(1, sizeof(*ret) ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( ( ret->data = XLALCalloc ( numIFOs, sizeof(ret->data[0]) ) ) != NULL, XLAL_ENOMEM );
  ret->length = numIFOs;
  ret->Sinv_Tsft = input->common.multiNoiseWeights->Sinv_Tsft;
  REAL8 ooSinv_Tsft = 1.0 / input->common.multiNoiseWeights->Sinv_Tsft;
  for ( UINT4 X = 0; X < numIFOs; X++ )
    {
      UINT4 length = input->common.multiNoiseWeights->data[X]->length;
      XLAL_CHECK_NULL ( ( ret->data[X] = XLALCreateREAL8Vector ( length ) ) != NULL, XLAL_EFUNC );
      // ComputeFstat stores *un-normalized* weights internally, so re-normalize them for returning them
      XLAL_CHECK_NULL ( input->common.multiNoiseWeights->isNotNormalized, XLAL_EERR, "BUG: noise weights stored in FstatInput must be unnormalized!\n" );
      XLAL_CHECK_NULL ( XLALVectorScaleREAL8 ( ret->data[X]->data, ooSinv_Tsft, input->common.multiNoiseWeights->data[X]->data, length) == XLAL_SUCCESS, XLAL_EFUNC);
    }

  return ret;

} // XLALGetFstatInputNoiseWeights()

///
/// Returns the multi-detector state series stored in a \c FstatInput structure.
///
const MultiDetectorStateSeries*
XLALGetFstatInputDetectorStates ( const FstatInput* input       ///< [in] \c FstatInput structure.
                                  )
{
  // Check input
  XLAL_CHECK_NULL ( input != NULL, XLAL_EINVAL );

  return input->common.multiDetectorStates;

} // XLALGetFstatInputDetectorStates()

///
/// Compute the \f$\mathcal{F}\f$-statistic over a band of frequencies.
///
int
XLALComputeFstat ( FstatResults **Fstats,               ///< [in/out] Address of a pointer to a #FstatResults results structure; if \c NULL, allocate here.
                   FstatInput *input,                   ///< [in] Input data structure created by one of the setup functions.
                   const PulsarDopplerParams *doppler,  ///< [in] Doppler parameters, including starting frequency, at which to compute \f$2\mathcal{F}\f$
                   const UINT4 numFreqBins,             ///< [in] Number of frequencies at which the \f$2\mathcal{F}\f$ are to be computed. Must be 1 if XLALCreateFstatInput() was passed zero \c dFreq.
                   const FstatQuantities whatToCompute  ///< [in] Bit-field of which \f$\mathcal{F}\f$-statistic quantities to compute.
                   )
{
  // Check input
  XLAL_CHECK ( Fstats != NULL, XLAL_EINVAL);
  XLAL_CHECK ( input != NULL, XLAL_EINVAL);
  XLAL_CHECK ( doppler != NULL, XLAL_EINVAL);
  XLAL_CHECK ( doppler->asini >= 0, XLAL_EINVAL);
  XLAL_CHECK ( numFreqBins > 0, XLAL_EINVAL);
  XLAL_CHECK ( !input->singleFreqBin || numFreqBins == 1, XLAL_EINVAL, "numFreqBins must be 1 if XLALCreateFstatInput() was passed zero dFreq" );
  XLAL_CHECK ( whatToCompute < FSTATQ_LAST, XLAL_EINVAL);

  // Allocate results struct, if needed
  if ( (*Fstats) == NULL ) {
    XLAL_CHECK ( ((*Fstats) = XLALCalloc ( 1, sizeof(**Fstats) )) != NULL, XLAL_ENOMEM );
  }

  // Get constant pointer to common input data
  const FstatCommon *common = &input->common;
  const UINT4 numDetectors = common->detectors.length;

  // Enlarge result arrays if they are too small
  const BOOLEAN moreFreqBins = (numFreqBins > (*Fstats)->internalalloclen);
  const BOOLEAN moreDetectors = (numDetectors > (*Fstats)->numDetectors);
  if (moreFreqBins || moreDetectors)
    {
      // Enlarge multi-detector 2F array
      if ( (whatToCompute & FSTATQ_2F) && moreFreqBins )
        {
          (*Fstats)->twoF = XLALRealloc ( (*Fstats)->twoF, numFreqBins*sizeof((*Fstats)->twoF[0]) );
          XLAL_CHECK ( (*Fstats)->twoF != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->twoF to length %u", numFreqBins );
        }

      // Enlarge multi-detector Fa & Fb array
      if ( (whatToCompute & FSTATQ_FAFB) && moreFreqBins )
        {
          (*Fstats)->Fa = XLALRealloc( (*Fstats)->Fa, numFreqBins * sizeof((*Fstats)->Fa[0]) );
          XLAL_CHECK ( (*Fstats)->Fa != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->Fa to length %u", numFreqBins );
          (*Fstats)->Fb = XLALRealloc( (*Fstats)->Fb, numFreqBins * sizeof((*Fstats)->Fb[0]) );
          XLAL_CHECK ( (*Fstats)->Fb != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->Fb to length %u", numFreqBins );
        }

      // Enlarge 2F per detector arrays
      if ( (whatToCompute & FSTATQ_2F_PER_DET) && (moreFreqBins || moreDetectors) )
        {
          for ( UINT4 X = 0; X < numDetectors; ++X )
            {
              (*Fstats)->twoFPerDet[X] = XLALRealloc ( (*Fstats)->twoFPerDet[X], numFreqBins * sizeof((*Fstats)->twoFPerDet[X][0]) );
              XLAL_CHECK ( (*Fstats)->twoFPerDet[X] != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->twoFPerDet[%u] to length %u", X, numFreqBins );
            }
        }

      // Enlarge Fa & Fb per detector arrays
      if ( ( whatToCompute & FSTATQ_FAFB_PER_DET) && (moreFreqBins || moreDetectors) )
        {
          for ( UINT4 X = 0; X < numDetectors; ++X )
            {
              (*Fstats)->FaPerDet[X] = XLALRealloc ( (*Fstats)->FaPerDet[X], numFreqBins*sizeof((*Fstats)->FaPerDet[X][0]) );
              XLAL_CHECK( (*Fstats)->FaPerDet[X] != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->FaPerDet[%u] to length %u", X, numFreqBins );
              (*Fstats)->FbPerDet[X] = XLALRealloc ( (*Fstats)->FbPerDet[X], numFreqBins*sizeof((*Fstats)->FbPerDet[X][0]) );
              XLAL_CHECK( (*Fstats)->FbPerDet[X] != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->FbPerDet[%u] to length %u", X, numFreqBins );
            }
        }

      // Enlarge F-atoms per detector arrays, and initialise to NULL
      if ( (whatToCompute & FSTATQ_ATOMS_PER_DET) && moreFreqBins )
        {
          UINT4 kPrev = 0;
          if ( (*Fstats)->multiFatoms != NULL ) {
            kPrev = (*Fstats)->internalalloclen; // leave previously-used frequency-bins untouched
          }

          (*Fstats)->multiFatoms = XLALRealloc ( (*Fstats)->multiFatoms, numFreqBins*sizeof((*Fstats)->multiFatoms[0]) );
          XLAL_CHECK ( (*Fstats)->multiFatoms != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->multiFatoms to length %u", numFreqBins );

          for ( UINT4 k = kPrev; k < numFreqBins; ++k ) {
            (*Fstats)->multiFatoms[k] = NULL;
          }

        } // if Atoms_per_det to enlarge

      // Update allocated length of arrays
      (*Fstats)->internalalloclen = numFreqBins;

    } // if (moreFreqBins || moreDetectors)

  // Extrapolate parameters in 'doppler' to SFT mid-time
  PulsarDopplerParams midDoppler = (*doppler);
  {
    const REAL8 dtau = XLALGPSDiff ( &common->midTime, &doppler->refTime );
    XLAL_CHECK ( XLALExtrapolatePulsarSpins ( midDoppler.fkdot, midDoppler.fkdot, dtau ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  midDoppler.refTime = common->midTime;

  // Initialise result struct parameters
  (*Fstats)->doppler      = midDoppler;
  (*Fstats)->dFreq        = input->singleFreqBin ? 0 : common->dFreq;
  (*Fstats)->numFreqBins  = numFreqBins;
  (*Fstats)->numDetectors = numDetectors;
  XLAL_INIT_MEM ( (*Fstats)->detectorNames);
  for (UINT4 X = 0; X < numDetectors; ++X) {
    strncpy ( (*Fstats)->detectorNames[X], common->detectors.sites[X].frDetector.prefix, 2 );
  }
  (*Fstats)->whatWasComputed = whatToCompute;

  // Call the appropriate method function to compute the F-statistic
  XLAL_CHECK ( (input->method_funcs.compute_func) ( *Fstats, common, input->method_data ) == XLAL_SUCCESS, XLAL_EFUNC );

  (*Fstats)->doppler = (*doppler);
  // Record the internal reference time used, which is required to compute a correct global signal phase
  (*Fstats)->refTimePhase = midDoppler.refTime;

  return XLAL_SUCCESS;

} // XLALComputeFstat()

///
/// Free all memory associated with a \c FstatInput structure.
///
void
XLALDestroyFstatInput ( FstatInput* input       ///< [in] \c FstatInput structure to be freed.
                        )
{
  if ( input == NULL ) {
    return;
  }
  if ( input->common.isTimeslice )
    {
      XLAL_CHECK_VOID ( input->method < FMETHOD_RESAMP_GENERIC, XLAL_EINVAL,
                        "Something is wrong: 'isTimeslice==TRUE' for non-LALDemod F-stat method '%s' is not supported!\n", XLALGetFstatInputMethodName(input));
      XLALDestroyFstatInputTimeslice_common ( &input->common );
      XLALDestroyFstatInputTimeslice_Demod ( input->method_data);
      XLALFree ( input );
      return;
    }

  XLALDestroyMultiTimestamps ( input->common.multiTimestamps );
  XLALDestroyMultiNoiseWeights ( input->common.multiNoiseWeights );
  XLALDestroyMultiDetectorStateSeries ( input->common.multiDetectorStates );

  // Release a reference to 'common.workspace'; if there are no more outstanding references ...
  if ( --(*input->workspace_refcount) == 0 ) {
    XLALPrintInfo( "%s: workspace reference count = %i, freeing workspace\n", __func__, *input->workspace_refcount );

    // If allocated, free method-specific workspace using destructor function
    if ( input->common.workspace != NULL ) {
      (input->method_funcs.workspace_destroy_func) ( input->common.workspace );
    }

    // Free memory for workspace reference count
    XLALFree ( input->workspace_refcount );

  } else {
    XLALPrintInfo( "%s: workspace reference count = %i\n", __func__, *input->workspace_refcount );
  }

  if ( input->method_data != NULL ) {
    // Free method-specific data using destructor function
    (input->method_funcs.method_data_destroy_func) ( input->method_data );
  }

  XLALFree ( input );

  return;
} // XLALDestroyFstatInput()

///
/// Free all memory associated with a #FstatResults structure.
///
void
XLALDestroyFstatResults ( FstatResults* Fstats  ///< [in] #FstatResults structure to be freed.
                          )
{
  if ( Fstats == NULL ) {
    return;
  }

  XLALFree ( Fstats->twoF );
  XLALFree ( Fstats->Fa );
  XLALFree ( Fstats->Fb );
  for ( UINT4 X = 0; X < PULSAR_MAX_DETECTORS; ++X )
    {
      XLALFree ( Fstats->twoFPerDet[X] );
      XLALFree ( Fstats->FaPerDet[X] );
      XLALFree ( Fstats->FbPerDet[X] );
    }
  if ( Fstats->multiFatoms != NULL )
    {
      for ( UINT4 n = 0; n < Fstats->internalalloclen; ++n )
        {
          XLALDestroyMultiFstatAtomVector ( Fstats->multiFatoms[n] );
        }
      XLALFree ( Fstats->multiFatoms );
    }

  XLALFree ( Fstats );

  return;
} // XLALDestroyFstatResults()

///
/// Add +4 to any multi-detector or per-detector 2F values computed by XLALComputeFstat().
/// This is for compatibility with programs which expect this normalisation if SFTs do not
/// contain noise, e.g. \c lalapps_ComputeFstatistic with the \c --SignalOnly option.
///
int
XLALAdd4ToFstatResults ( FstatResults* Fstats    ///< [in/out] #FstatResults structure.
                         )
{
  // Check input
  XLAL_CHECK( Fstats != NULL, XLAL_EINVAL );

  // Add +4 to multi-detector 2F array
  if ( Fstats->whatWasComputed & FSTATQ_2F )
    {
      for ( UINT4 k = 0; k < Fstats->numFreqBins; ++k ) {
        Fstats->twoF[k] += 4;
      }
    }

  // Add +4 to 2F per detector arrays
  if ( Fstats->whatWasComputed & FSTATQ_2F_PER_DET )
    {
      for ( UINT4 X = 0; X < Fstats->numDetectors; ++X ) {
        for ( UINT4 k = 0; k < Fstats->numFreqBins; ++k ) {
          Fstats->twoFPerDet[X][k] += 4;
        }
      }
    }

  return XLAL_SUCCESS;

} // XLALAdd4ToFstatResults()

///
/// Compute single-or multi-IFO Fstat '2F' from multi-IFO 'atoms'
///
REAL4
XLALComputeFstatFromAtoms ( const MultiFstatAtomVector *multiFstatAtoms,   ///< [in] Multi-detector atoms
                            const INT4                 X                   ///< [in] Detector number, give -1 for multi-Fstat
                            )
{
  // ----- check input parameters and report errors
  XLAL_CHECK_REAL4 ( multiFstatAtoms && multiFstatAtoms->data && multiFstatAtoms->data[0]->data, XLAL_EINVAL, "Empty pointer as input parameter." );
  XLAL_CHECK_REAL4 ( multiFstatAtoms->length > 0, XLAL_EBADLEN, "Input MultiFstatAtomVector has zero length. (i.e., no detectors)" );
  XLAL_CHECK_REAL4 ( X >= -1, XLAL_EDOM, "Invalid detector number X=%d. Only nonnegative numbers, or -1 for multi-F, are allowed.", X );
  XLAL_CHECK_REAL4 ( ( X < 0 ) || ( (UINT4)(X) <= multiFstatAtoms->length-1 ), XLAL_EDOM, "Requested X=%d, but FstatAtoms only have length %d.", X, multiFstatAtoms->length );

  // internal detector index Y to do both single- and multi-F case
  UINT4 Y, Ystart, Yend;
  if ( X == -1 ) { /* loop through all detectors to get multi-Fstat */
    Ystart = 0;
    Yend   = multiFstatAtoms->length-1;
  }
  else { /* just compute single-Fstat for 1 IFO */
    Ystart = X;
    Yend   = X;
  }

  // set up temporary Fatoms and matrix elements for summations
  REAL4 mmatrixA = 0.0, mmatrixB = 0.0, mmatrixC = 0.0;
  REAL4 twoF = 0.0;
  COMPLEX8 Fa, Fb;
  Fa = 0.0;
  Fb = 0.0;

  for (Y = Ystart; Y <= Yend; Y++ ) /* loop through detectors */
    {
      UINT4 alpha, numSFTs;
      numSFTs = multiFstatAtoms->data[Y]->length;
      XLAL_CHECK_REAL4 ( numSFTs > 0, XLAL_EDOM, "Input FstatAtomVector has zero length. (i.e., no timestamps for detector X=%d)", Y );

      for ( alpha = 0; alpha < numSFTs; alpha++ ) /* loop through SFTs */
        {
          FstatAtom *thisAtom = &multiFstatAtoms->data[Y]->data[alpha];
          /* sum up matrix elements and Fa, Fb */
          mmatrixA += thisAtom->a2_alpha;
          mmatrixB += thisAtom->b2_alpha;
          mmatrixC += thisAtom->ab_alpha;
          Fa += thisAtom->Fa_alpha;
          Fb += thisAtom->Fb_alpha;
        } /* loop through SFTs */

    } // loop through detectors

  // compute determinant and final Fstat (not twoF!)
  REAL4 Dinv = 1.0 / ( mmatrixA * mmatrixB - SQ(mmatrixC) );

  twoF = XLALComputeFstatFromFaFb ( Fa, Fb, mmatrixA, mmatrixB, mmatrixC, 0, Dinv );

  return twoF;

} // XLALComputeFstatFromAtoms()

///
/// If user asks for a 'best' #FstatMethodType, find and select it
///
static int
XLALSelectBestFstatMethod ( FstatMethodType *method )
{
  switch ( *method ) {

  case FMETHOD_DEMOD_BEST:
  case FMETHOD_RESAMP_BEST:
    // If user asks for a 'best' method:
    //   Decrement the current method, then check for the first available Fstat method. This assumes the FstatMethodType enum is ordered as follows:
    //     FMETHOD_..._GENERIC,      (must **always** avaiable)
    //     FMETHOD_..._OPTIMISED,    (always avaiable)
    //     FMETHOD_..._SUPERFAST     (not always available; requires special hardware)
    //     FMETHOD_..._BEST          (must **always** avaiable)
    XLALPrintInfo( "%s: trying to find best available Fstat method for '%s'\n", __func__, FstatMethodNames[*method] );
    while ( !XLALFstatMethodIsAvailable( --( *method ) ) ) {
      XLAL_CHECK ( FMETHOD_START < *method, XLAL_EFAILED );
      XLALPrintInfo( "%s: Fstat method '%s' is unavailable\n",  __func__, FstatMethodNames[*method] );
    }
    XLALPrintInfo( "%s: Fstat method '%s' is available; selected as best method\n", __func__, FstatMethodNames[*method] );
    break;

  default:
    // If user asks for a specific method:
    //   Check that method is available
    XLAL_CHECK ( XLALFstatMethodIsAvailable(*method), XLAL_EINVAL, "Fstat method '%s' is unavailable", FstatMethodNames[*method] );
    break;
  }
  return XLAL_SUCCESS;
}

///
/// Return true if given #FstatMethodType corresponds to a valid and *available* Fstat method, false otherwise
///
int
XLALFstatMethodIsAvailable ( FstatMethodType method )
{
  switch (method) {

  case FMETHOD_DEMOD_GENERIC:
  case FMETHOD_DEMOD_BEST:
  case FMETHOD_RESAMP_GENERIC:
  case FMETHOD_RESAMP_BEST:
    // 'Generic' and 'Best' methods must **always** be available
    return 1;

  case FMETHOD_DEMOD_OPTC:
    // This method is always avalable
    return 1;

  case FMETHOD_DEMOD_ALTIVEC:
    // This method is available only if compiled with Altivec support
#ifdef HAVE_ALTIVEC
    return 1;
#else
    return 0;
#endif

  case FMETHOD_DEMOD_SSE:
    // This method is available only if compiled with SSE support,
    // and SSE is available on the current execution machine
#ifdef HAVE_SSE_COMPILER
    return LAL_HAVE_SSE_RUNTIME();
#else
    return 0;
#endif

  default:
    return 0;

  }
} // XLALFstatMethodIsAvailable()

///
/// Return pointer to a static string giving the name of the #FstatMethodType \p method
///
const CHAR *
XLALFstatMethodName ( FstatMethodType method )
{
  XLAL_CHECK_NULL( method < XLAL_NUM_ELEM(FstatMethodNames) && FstatMethodNames[method] != NULL,
                   XLAL_EINVAL, "FstatMethodType = %i is invalid", method );
  return FstatMethodNames[method];
}

///
/// Return pointer to a static array of all (available) #FstatMethodType choices.
/// This data is used by the UserInput module to parse a user enumeration.
///
const UserChoices *
XLALFstatMethodChoices ( void )
{
  static int firstCall = 1;
  static UserChoices choices;
  if ( firstCall )
    {
      XLAL_INIT_MEM ( choices );
      size_t i = 0;
      for (FstatMethodType m = FMETHOD_START + 1; m < FMETHOD_END; m++ )
        {
          if ( ! XLALFstatMethodIsAvailable(m) ) {
            continue;
          }
          FstatMethodType bm = m;
          XLAL_CHECK_NULL ( XLALSelectBestFstatMethod( &bm ) == XLAL_SUCCESS, XLAL_EFAULT );
          if ( bm == m ) {
            // add an Fstat method
            choices[i].name = FstatMethodNames[m];
            choices[i].val = m;
            i++;
          } else {
            // add a "best" Fstat method twice, so that:
            // - name of "best" method will appear as equal to the actual method it selects,
            //   in the help string, e.g. "...|DemodSuperFast=DemodBest|..."
            // - FMETHOD_..._BEST will print as "...Best" in e.g. help default argument, log output
            choices[i].name = FstatMethodNames[m];
            choices[i].val = bm;
            i++;
            choices[i].name = FstatMethodNames[m];
            choices[i].val = m;
            i++;
          }
        } // for m < FMETHOD_LAST

      firstCall = 0;

    } // if firstCall

  return (const UserChoices *) &choices;
} // XLALFstatMethodChoices()

/// Return measured values and details about generic F-statistic timing and method-specific timing model,
// including static pointers to help-string documenting the generic F-stat timing, and the method-specific timing model.
/// See also https://dcc.ligo.org/LIGO-T1600531-v4 for details about the timing model.
int
XLALGetFstatTiming ( const FstatInput* input, FstatTimingGeneric *timingGeneric, FstatTimingModel *timingModel )
{
  XLAL_CHECK ( input != NULL, XLAL_EINVAL );
  XLAL_CHECK ( timingGeneric != NULL, XLAL_EINVAL );
  XLAL_CHECK ( timingModel != NULL, XLAL_EINVAL );

  if ( input->method >= FMETHOD_DEMOD_GENERIC && input->method <= FMETHOD_DEMOD_BEST)
    {
      XLAL_CHECK ( XLALGetFstatTiming_Demod ( input->method_data, timingGeneric, timingModel ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  else if ( input->method >= FMETHOD_RESAMP_GENERIC && input->method <= FMETHOD_RESAMP_BEST)
    {
      XLAL_CHECK ( XLALGetFstatTiming_Resamp ( input->method_data, timingGeneric, timingModel ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  else
    {
      XLAL_ERROR ( XLAL_EINVAL, "Unsupported F-stat method '%s'\n", FstatMethodNames [ input->method ] );
    }

  timingGeneric->help = FstatTimingGenericHelp;	// set static help-string pointer (not used or set otherwise)

  return XLAL_SUCCESS;

} // XLALGetFstatTiming()

int
XLALAppendFstatTiming2File ( const FstatInput* input, FILE *fp, BOOLEAN printHeader )
{
  XLAL_CHECK ( input != NULL, XLAL_EINVAL );
  XLAL_CHECK ( fp != NULL, XLAL_EINVAL );

  FstatTimingGeneric XLAL_INIT_DECL ( tiGen);
  FstatTimingModel XLAL_INIT_DECL ( tiModel );
  XLAL_CHECK ( XLALGetFstatTiming ( input, &tiGen, &tiModel ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( printHeader )
    {
      fprintf ( fp, "%s\n", tiGen.help );
      fprintf ( fp, "%s\n", tiModel.help );
      // generic F-stat timing header line
      fprintf ( fp, "%%%%%8s %10s %4s %10s %10s %11s %10s ", "NCalls", "NFbin", "Ndet", "tauF_eff", "tauF_core", "tauF_buffer", "b");
      // method-specific F-stat timing model header line
      fprintf (fp, "|");
      for ( UINT4 i = 0; i < tiModel.numVariables; i ++ ) {
        fprintf (fp, " %13s", tiModel.names[i] );
      }
      fprintf (fp, "\n");
    } // if (printHeader)

  // generic F-stat timing values
  fprintf (fp, "%10.0f %10d %4d %10.2e %10.2e %11.2e %10.2e ",
           tiGen.NCalls, tiGen.NFbin, tiGen.Ndet, tiGen.tauF_eff, tiGen.tauF_core, tiGen.tauF_buffer, tiGen.NBufferMisses/tiGen.NCalls );

  // method-specific F-stat timing model values
  fprintf (fp, " "); // for '|'
  for ( UINT4 i = 0; i < tiModel.numVariables; i ++ ) {
    fprintf (fp, " %13g", tiModel.values[i] );
  }
  fprintf (fp, "\n");

  return XLAL_SUCCESS;

} // AppendFstatTiming2File()

///
/// Extracts the resampled timeseries from a given FstatInput structure (must have been initialized previously by XLALCreateFstatInput() and a call to XLALComputeFstat()).
///
/// \note This function only returns a pointer to the time-series that are still part of the FstatInput structure, and will by managed by F-statistic API calls.
/// Therefore do *not* free the resampled timeseries with XLALDestroyMultiCOMPLEX8TimeSeries(), only
/// free the complete Fstat input structure with XLALDestroyFstatInputVector() at the end once its no longer needed.
///
int
XLALExtractResampledTimeseries ( MultiCOMPLEX8TimeSeries **multiTimeSeries_SRC_a, ///< [out] \c multi-detector SRC-frame timeseries, multiplied by AM function a(t).
                                 MultiCOMPLEX8TimeSeries **multiTimeSeries_SRC_b, ///< [out] \c multi-detector SRC-frame timeseries, multiplied by AM function b(t).
                                 const FstatInput *input ///< [in] \c FstatInput structure.
                                 )
{
  XLAL_CHECK ( input != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ( multiTimeSeries_SRC_a != NULL ) && ( multiTimeSeries_SRC_b != NULL ) , XLAL_EINVAL );
  XLAL_CHECK ( (input->method >= FMETHOD_RESAMP_GENERIC) && (input->method <= FMETHOD_RESAMP_BEST), XLAL_EINVAL,
               "%s() only works for resampling-Fstat methods, not with '%s'\n", __func__, XLALGetFstatInputMethodName ( input ) );

  XLAL_CHECK ( XLALExtractResampledTimeseries_intern ( multiTimeSeries_SRC_a, multiTimeSeries_SRC_b, input->method_data ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;
} // XLALExtractResampledTimeseries()

///
/// Create and return an FstatInput 'timeslice' for given input FstatInput object [must be using LALDemod Fstat method!]
/// which covers only the time interval  ['minStartGPS','maxStartGPS'), using conventions of XLALCWGPSinRange().
///
/// The returned FstatInput structure is equivalent to one that would have been produced for the given time-interval
/// in the first place, except that it uses "containers" and references the data in the original FstatInput object.
/// A special flag is set in the FstatInput object to notify the destructor to only free the data-containers.
///
/// Note: if a timeslice is empty for any detector, the corresponding timeseries 'container' will be allocated and
/// will have length==0 and data==NULL.
///
int
XLALFstatInputTimeslice ( FstatInput ** slice,                ///< [out] Address of a pointer to a \c FstatInput structure
                          const FstatInput* input,            ///< [in] Input data structure
                          const LIGOTimeGPS *minStartGPS,     ///< [in] Minimum starting GPS time
                          const LIGOTimeGPS *maxStartGPS      ///< [in] Maximum starting GPS time
                        )
{
  XLAL_CHECK ( input != NULL, XLAL_EINVAL );
  XLAL_CHECK ( slice != NULL && (*slice) == NULL, XLAL_EINVAL );
  XLAL_CHECK ( minStartGPS != NULL, XLAL_EINVAL );
  XLAL_CHECK ( maxStartGPS != NULL, XLAL_EINVAL );
  XLAL_CHECK( XLALGPSCmp( minStartGPS, maxStartGPS ) < 1 , XLAL_EINVAL , "minStartGPS (%"LAL_GPS_FORMAT") is greater than maxStartGPS (%"LAL_GPS_FORMAT")\n",
              LAL_GPS_PRINT(*minStartGPS), LAL_GPS_PRINT(*maxStartGPS) );

  // only supported for 'LALDemod' Fstat methods
  XLAL_CHECK ( input->method < FMETHOD_RESAMP_GENERIC, XLAL_EINVAL, "This function is not avavible for the chosen FstatMethod '%s'!", XLALGetFstatInputMethodName ( input ) );

  const FstatCommon *common = &(input->common);
  UINT4 numIFOs = common->detectors.length;

  // ---------- Find start- and end- indices of the requested timeslice for all detectors
  const MultiLIGOTimeGPSVector *mTS = common->multiTimestamps;
  UINT4 XLAL_INIT_DECL ( iStart, [PULSAR_MAX_DETECTORS] );
  UINT4 XLAL_INIT_DECL ( iEnd, [PULSAR_MAX_DETECTORS] );
  for ( UINT4 X = 0; X < numIFOs; X ++ ) {
    XLAL_CHECK ( XLALFindTimesliceBounds ( &(iStart[X]), &(iEnd[X]), input->common.multiTimestamps->data[X], minStartGPS, maxStartGPS ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // ----- prepare new top-level 'containers' for the timeslices in the 'common' part of FstatInput:
  MultiLIGOTimeGPSVector *multiTimestamps = NULL;
  XLAL_CHECK ( ( multiTimestamps = XLALCalloc ( 1, sizeof(*multiTimestamps) ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( multiTimestamps->data = XLALCalloc ( numIFOs, sizeof(*multiTimestamps->data) ) ) != NULL, XLAL_ENOMEM );
  multiTimestamps->length = numIFOs;

  MultiNoiseWeights *multiNoiseWeights = NULL;
  XLAL_CHECK ( ( multiNoiseWeights = XLALCalloc ( 1, sizeof(*multiNoiseWeights) ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( multiNoiseWeights->data = XLALCalloc ( numIFOs, sizeof(*multiNoiseWeights->data) ) ) != NULL, XLAL_ENOMEM );
  multiNoiseWeights->length = numIFOs;

  MultiDetectorStateSeries *multiDetectorStates = NULL;
  XLAL_CHECK ( ( multiDetectorStates = XLALCalloc ( 1, sizeof(*multiDetectorStates) ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( multiDetectorStates->data = XLALCalloc ( numIFOs, sizeof(*multiDetectorStates->data) ) ) != NULL, XLAL_ENOMEM );
  multiDetectorStates->length = numIFOs;

  // determine earliest and latest SFT times of timeslice over all IFOs, in order to correctly set 'FstatCommon->midTime' value
  LIGOTimeGPS earliest_time = {LAL_INT4_MAX,0};
  LIGOTimeGPS latest_time = {0,0};

  // for updating noise-normalization factor Sinv_Tsft over timesclie
  UINT4 numSFTsSlice = 0;
  REAL8 sumWeightsSlice = 0;

  // ----- loop over detectors, create per-detector containers and set them to the requested timeslice
  for ( UINT4 X = 0; X < numIFOs; X ++ )
    {
      // prepare multiTimestamps
      XLAL_CHECK ( ( multiTimestamps->data[X] = XLALCalloc ( 1, sizeof(*multiTimestamps->data[X]) ) ) != NULL, XLAL_EFUNC );
      multiTimestamps->data[X]->deltaT = common->multiTimestamps->data[X]->deltaT;

      // prepare multiNoiseWeights
      multiNoiseWeights->data[X] = NULL;
      XLAL_CHECK ( ( multiNoiseWeights->data[X] = XLALCalloc ( 1, sizeof(*multiNoiseWeights->data[X]) ) ) !=NULL, XLAL_EFUNC );

      // prepare multiDetectorStates, copy struct values
      XLAL_CHECK ( ( multiDetectorStates->data[X] = XLALCalloc ( 1, sizeof(*multiDetectorStates->data[X]) ) ) != NULL, XLAL_EFUNC );
      multiDetectorStates->data[X]->detector = common->multiDetectorStates->data[X]->detector;
      multiDetectorStates->data[X]->system   = common->multiDetectorStates->data[X]->system;
      multiDetectorStates->data[X]->deltaT   = common->multiDetectorStates->data[X]->deltaT;

      //---- set the timeslices
      if ( iStart[X] > iEnd[X] ) {
        continue;	// empty slice
      }

      UINT4 sliceLength =  iEnd[X] - iStart[X] + 1;	// guaranteed >= 1

      // slice multiTimestamps
      multiTimestamps->data[X]->length = sliceLength;
      multiTimestamps->data[X]->data   = &( mTS->data[X]->data[iStart[X]] );

      // keep track of earliest and latest SFT times included in timeslice (for computing 'midTime')
      if ( ( XLALGPSCmp ( &(multiTimestamps->data[X]->data[0]), &earliest_time ) == -1 ) ) {
        earliest_time = multiTimestamps->data[X]->data[0];
      }
      if ( XLALGPSCmp ( &(multiTimestamps->data[X]->data[sliceLength - 1]), &latest_time ) == 1 ) {
        latest_time = multiTimestamps->data[X]->data[sliceLength - 1];
      }

      // slice multiNoiseWeights
      multiNoiseWeights->data[X]->length = sliceLength;
      multiNoiseWeights->data[X]->data   = &( common->multiNoiseWeights->data[X]->data[iStart[X]] );

      // sum the weights for updating the overall noise-normalisation factor Sinv_Tsft
      numSFTsSlice += sliceLength;
      for ( UINT4 alpha = 0; alpha < sliceLength; alpha++ ) {
        sumWeightsSlice += multiNoiseWeights->data[X]->data[alpha];
      }

      // slice multiDetectorStates
      multiDetectorStates->data[X]->length = sliceLength;
      multiDetectorStates->data[X]->data   = &( common->multiDetectorStates->data[X]->data[iStart[X]] );

    } // for X < numIFOs

  // calculate new data midTime for timeslice
  XLALGPSAdd ( &latest_time, multiDetectorStates->data[0]->deltaT );
  REAL8 TspanSlice = XLALGPSDiff( &latest_time, &earliest_time );
  LIGOTimeGPS midTimeSlice = earliest_time;
  XLALGPSAdd ( &midTimeSlice, 0.5 * TspanSlice );

  // update noise-weight normalisation factor
  multiNoiseWeights->Sinv_Tsft = sumWeightsSlice / numSFTsSlice ;
  multiNoiseWeights->isNotNormalized = (1 == 1);

  // allocate memory and copy the orginal FstatInput struct
  XLAL_CHECK ( ( (*slice) = XLALCalloc ( 1 , sizeof(*input) ) ) != NULL, XLAL_ENOMEM );
  memcpy ( (*slice), input, sizeof ( *input ) );

  (*slice)->common.isTimeslice         = (1==1); // This is a timeslice
  (*slice)->common.midTime             = midTimeSlice;
  (*slice)->common.multiTimestamps     = multiTimestamps;
  (*slice)->common.multiDetectorStates = multiDetectorStates;
  (*slice)->common.multiNoiseWeights   = multiNoiseWeights;

  (*slice)->method_data = XLALFstatInputTimeslice_Demod ( input->method_data, iStart, iEnd );
  XLAL_CHECK ( (*slice)->method_data != NULL, XLAL_EFUNC );

  return XLAL_SUCCESS;

} // XLALFstatInputTimeslice()


void
XLALDestroyFstatInputTimeslice_common ( FstatCommon *common )
{
  if ( !common ) {
    return;
  }

  for ( UINT4 X=0; X < common->multiTimestamps->length; X ++ )
    {
      XLALFree ( common->multiTimestamps->data[X] );
      XLALFree ( common->multiDetectorStates->data[X] );
      XLALFree ( common->multiNoiseWeights->data[X] );
    }

  XLALFree ( common->multiTimestamps->data );
  XLALFree ( common->multiDetectorStates->data );
  XLALFree ( common->multiNoiseWeights->data );
  XLALFree ( common->multiTimestamps );
  XLALFree ( common->multiNoiseWeights );
  XLALFree ( common->multiDetectorStates );

  return;

} // XLALDestroyFstatInputTimeslice_common()
