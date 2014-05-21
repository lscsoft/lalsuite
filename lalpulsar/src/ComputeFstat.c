//
// Copyright (C) 2012, 2013, 2014 Karl Wette
// Copyright (C) 2007 Chris Messenger
// Copyright (C) 2006 John T. Whelan, Badri Krishnan
// Copyright (C) 2005, 2006, 2007, 2010 Reinhard Prix
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

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include <lal/ComputeFstat.h>
#include <lal/Factorial.h>
#include <lal/ComplexFFT.h>
#include <lal/TimeSeries.h>
#include <lal/LALComputeAM.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/LFTandTSutils.h>
#include <lal/CWFastMath.h>
#include <lal/NormalizeSFTRngMed.h>

// ----- macro definitions
#define MYSIGN(x) ( ((x) < 0) ? (-1.0):(+1.0) )
#define SQ(x) ( (x) * (x) )

// ---------- Internal struct definitions ---------- //

// Common input data for F-statistic methods
typedef struct {
  MultiLALDetector detectors;                           // List of detectors
  MultiLIGOTimeGPSVector *timestamps;                   // Multi-detector list of SFT timestamps
  MultiNoiseWeights *noiseWeights;                      // Multi-detector noise weights
  MultiDetectorStateSeries *detectorStates;             // Multi-detector state series
  const EphemerisData *ephemerides;                     // Ephemerides for the time-span of the SFTs
  SSBprecision SSBprec;                                 // Barycentric transformation precision
  FstatMethodType FstatMethod;                          // Method to use for computing the F-statistic
} FstatInput_Common;

// Input data specific to F-statistic methods
typedef struct tagFstatInput_Demod FstatInput_Demod;
typedef struct tagFstatInput_Resamp FstatInput_Resamp;

// Internal definition of input data structure
struct tagFstatInput {
  FstatInput_Common* common;                        // Common input data
  FstatInput_Demod* demod;                          // Demodulation input data
  FstatInput_Resamp* resamp;                        // Resampling input data
};

// ----- internal prototypes
static REAL8 ComputeFstatFromFaFb ( COMPLEX16 Fa, COMPLEX16 Fb, REAL8 A, REAL8 B, REAL8 C, REAL8 E, REAL8 Dinv );

// ---------- Check for various computer architectures ---------- //

#if defined(HAVE_SSE) || defined(__SSE__)
#define CFS_HAVE_SSE 1
#else
#define CFS_HAVE_SSE 0
#endif

#if defined(HAVE_ALTIVEC) || defined(__ALTIVEC__)
#define CFS_HAVE_ALTIVEC 1
#else
#define CFS_HAVE_ALTIVEC 0
#endif

// ---------- Include F-statistic method implementations ---------- //

#if CFS_HAVE_SSE
const int FMETHOD_DEMOD_BEST = FMETHOD_DEMOD_SSE;
#elif CFS_HAVE_ALTIVEC
const int FMETHOD_DEMOD_BEST = FMETHOD_DEMOD_ALTIVEC;
#else
const int FMETHOD_DEMOD_BEST = FMETHOD_DEMOD_OPTC;
#endif

const int FMETHOD_RESAMP_BEST = FMETHOD_RESAMP_GENERIC;

static const struct {
  const char *const name;
  BOOLEAN available;
} FstatMethodNames[FMETHOD_END] = {
  [FMETHOD_DEMOD_GENERIC]	= {"DemodGeneric",	1 },
  [FMETHOD_DEMOD_OPTC] 		= {"DemodOptC", 	1 },
  [FMETHOD_DEMOD_SSE] 		= {"DemodSSE", 		CFS_HAVE_SSE},
  [FMETHOD_DEMOD_ALTIVEC] 	= {"DemodAltivec", 	CFS_HAVE_ALTIVEC},

  [FMETHOD_RESAMP_GENERIC]	= {"ResampGeneric", 	1 }
} ;

#include "ComputeFstat_Demod.c"
#include "ComputeFstat_Resamp.c"

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
  FstatInputVector* inputs = XLALCalloc ( 1, sizeof(*inputs) );
  XLAL_CHECK_NULL ( inputs != NULL, XLAL_ENOMEM);
  inputs->length = length;

  // Allocate and initialise vector data
  if (inputs->length > 0)
    {
      inputs->data = XLALCalloc ( inputs->length, sizeof(inputs->data[0]) );
      XLAL_CHECK_NULL(inputs->data != NULL, XLAL_ENOMEM);
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
  FstatAtomVector* atoms = XLALCalloc ( 1, sizeof(*atoms) );
  XLAL_CHECK_NULL ( atoms != NULL, XLAL_ENOMEM );
  atoms->length = length;

  // Allocate and initialise vector data
  if (atoms->length > 0)
    {
      atoms->data = XLALCalloc (atoms->length, sizeof(atoms->data[0]) );
      XLAL_CHECK_NULL ( atoms->data != NULL, XLAL_ENOMEM );
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
  MultiFstatAtomVector* multiAtoms = XLALCalloc(1, sizeof(*multiAtoms));
  XLAL_CHECK_NULL ( multiAtoms != NULL, XLAL_ENOMEM );
  multiAtoms->length = length;

  // Allocate and initialise vector data
  if ( multiAtoms->length > 0 )
    {
      multiAtoms->data = XLALCalloc ( multiAtoms->length, sizeof(multiAtoms->data[0]) );
      XLAL_CHECK_NULL ( multiAtoms->data != NULL, XLAL_ENOMEM );
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
XLALCreateFstatInput ( const SFTCatalog *SFTcatalog,		  ///< [in] Catalog of SFTs to either load from files, or generate in memory.
                                                                  ///< The \c locator field of each ::SFTDescriptor must be \c !NULL for SFT loading,
                                                                  ///< and \c NULL for SFT generation.

                       const REAL8 minCoverFreq,		  ///< [in] Minimum instantaneous frequency which will be covered over the SFT time span.
                       const REAL8 maxCoverFreq,		  ///< [in] Maximum instantaneous frequency which will be covered over the SFT time span.
                       const PulsarParamsVector *injectSources,	  ///< [in] Optional vector of parameters of CW signals to simulate and inject.

                       const MultiNoiseFloor *injectSqrtSX,	  ///< [in] Optional array of single-sided PSD values governing fake Gaussian noise generation.
                                                                  ///< If supplied, then fake Gaussian noise with the given PSD values will be added to the SFTs.

                       const MultiNoiseFloor *assumeSqrtSX,	  ///< [in] Optional array of single-sided PSD values governing the calculation of SFT noise weights.
                                                                  ///< If supplied, then SFT noise weights are calculated from constant spectra with the given PSD
                                                                  ///< values; otherwise, SFT noise weights are calculated from PSDs computed from a running median
                                                                  ///< of the SFTs themselves.
                       const UINT4 runningMedianWindow,		  ///< [in] If SFT noise weights are calculated from the SFTs, the running median window length to use.

                       const EphemerisData *ephemerides,	  ///< [in] Ephemerides for the time-span of the SFTs.

                       const FstatMethodType FstatMethod,	  ///< [in] Method to use for computing the \f$\mathcal{F}\f$-statistic.
                       const FstatExtraParams *extraParams        ///< [in] Minor tuning or method-specific parameters.
                       )
{
  // Check catalog
  XLAL_CHECK_NULL ( SFTcatalog != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL ( SFTcatalog->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL ( SFTcatalog->data != NULL, XLAL_EFAULT );
  for ( UINT4 i = 1; i < SFTcatalog->length; ++i )
    {
      XLAL_CHECK_NULL ( (SFTcatalog->data[0].locator == NULL) == (SFTcatalog->data[i].locator == NULL), XLAL_EINVAL,
                  "All 'locator' fields of SFTDescriptors in 'SFTcatalog' must be either NULL or !NULL." );
    }

  // Create a multi-view of SFT catalog
  MultiSFTCatalogView *multiSFTcatalog = XLALGetMultiSFTCatalogView(SFTcatalog);
  XLAL_CHECK_NULL ( multiSFTcatalog != NULL, XLAL_EFUNC );
  const UINT4 numDetectors = multiSFTcatalog->length;

  // check that no single-SFT input vectors are given to avoid returning singular results:
  for ( UINT4 X=0; X < numDetectors; X ++ ) {
    XLAL_CHECK_NULL ( multiSFTcatalog->data[X].length > 1, XLAL_EINVAL, "Need more than 1 SFTs per Detector!\n" );
  }

  // Check remaining parameters
  XLAL_CHECK_NULL ( isfinite(minCoverFreq) && minCoverFreq > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL ( isfinite(maxCoverFreq) && maxCoverFreq > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL ( maxCoverFreq > minCoverFreq, XLAL_EINVAL );
  XLAL_CHECK_NULL ( injectSources == NULL || injectSources->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL ( injectSources == NULL || injectSources->data != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL ( injectSqrtSX == NULL || injectSqrtSX->length == numDetectors, XLAL_EINVAL );
  XLAL_CHECK_NULL ( assumeSqrtSX == NULL || assumeSqrtSX->length == numDetectors, XLAL_EINVAL );
  XLAL_CHECK_NULL ( ephemerides != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL ( extraParams != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL ( extraParams->SSBprec < SSBPREC_LAST, XLAL_EINVAL );

  // Determine whether to load and/or generate SFTs
  const BOOLEAN loadSFTs = (SFTcatalog->data[0].locator != NULL);
  const BOOLEAN generateSFTs = (injectSources != NULL) || (injectSqrtSX != NULL);
  XLAL_CHECK_NULL ( loadSFTs || generateSFTs, XLAL_EINVAL, "Can neither load nor generate SFTs with given parameters" );

  // Create top-level input data struct
  FstatInput* input;
  XLAL_CHECK_NULL ( (input = XLALCalloc ( 1, sizeof(*input) )) != NULL, XLAL_ENOMEM );

  // create common input data
  XLAL_CHECK_NULL ( (input->common = XLALCalloc ( 1, sizeof(*input->common))) != NULL, XLAL_ENOMEM );
  FstatInput_Common *const common = input->common;	// handy shortcut

  // create method-specific input data
  if ( XLALFstatMethodClassIsDemod ( FstatMethod ) )
    {
      XLAL_CHECK_NULL ( (input->demod = XLALCalloc ( 1, sizeof(FstatInput_Demod) )) != NULL, XLAL_ENOMEM );
      input->demod->Dterms = extraParams->Dterms;
    }
  else if ( XLALFstatMethodClassIsResamp ( FstatMethod ) )
    {
      XLAL_CHECK_NULL ( (input->resamp = XLALCalloc(1, sizeof(FstatInput_Resamp))) != NULL, XLAL_ENOMEM );
    }
  else
    {
      XLAL_ERROR_NULL ( XLAL_EINVAL, "Received invalid Fstat method enum '%d'\n", FstatMethod );
    }
  common->FstatMethod = FstatMethod;

  // Determine the time baseline of an SFT
  const REAL8 Tsft = 1.0 / SFTcatalog->data[0].header.deltaF;

  // Determine the frequency band to load or generate SFTs over
  REAL8 minFreq = minCoverFreq, maxFreq = maxCoverFreq;
  {
    // Determine whether the method being used requires extra frequency bins
    int extraBins = 0;
    if ( input->demod != NULL ) {
      extraBins = GetFstatExtraBins_Demod ( input->demod );
    } else if ( input->resamp != NULL ) {
      extraBins = GetFstatExtraBins_Resamp ( input->resamp );
    } else {
      XLAL_ERROR_NULL ( XLAL_EFAILED, "Invalid FstatInput struct passed to %s()", __func__);
    }
    XLAL_CHECK_NULL ( extraBins >= 0, XLAL_EFAILED );

    // Add number of extra frequency bins required by running median
    extraBins += runningMedianWindow/2 + 1;	// NOTE: this is currently needed irrespective of assumeSqrtSX usage!

    // Extend frequency range by number of extra bins times SFT bin width
    const REAL8 extraFreq = extraBins / Tsft;
    minFreq -= extraFreq;
    maxFreq += extraFreq;
  } // end: block to determine frequency-bins range

  // Extract detectors and timestamps from multi-view of SFT catalog
  XLAL_CHECK_NULL ( XLALMultiLALDetectorFromMultiSFTCatalogView( &common->detectors, multiSFTcatalog ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (common->timestamps = XLALTimestampsFromMultiSFTCatalogView(multiSFTcatalog)) != NULL,  XLAL_EFUNC );

  // Generate SFTs with injections and noise, if required
  MultiSFTVector *multiSFTs = NULL;
  if (generateSFTs)
    {
      // Initialise parameters struct for XLALCWMakeFakeMultiData()
      CWMFDataParams XLAL_INIT_DECL(MFDparams);
      MFDparams.fMin = minFreq;
      MFDparams.Band = maxFreq - minFreq;
      MFDparams.multiIFO = common->detectors;
      MFDparams.multiTimestamps = *(common->timestamps);
      MFDparams.randSeed = extraParams->randSeed;

      // Set noise floors if sqrtSX is given; otherwise noise floors are zero
      if ( injectSqrtSX != NULL ) {
        MFDparams.multiNoiseFloor = (*injectSqrtSX);
      } else {
        MFDparams.multiNoiseFloor.length = numDetectors;
      }

      // Generate SFTs with injections
      XLAL_CHECK_NULL ( XLALCWMakeFakeMultiData ( &multiSFTs, NULL, injectSources, &MFDparams, ephemerides ) == XLAL_SUCCESS, XLAL_EFUNC );

    } // if generateSFTs

  // Load SFTs, if required
  if (loadSFTs)
    {
      // If no SFTs have been generated, then load all SFTs at once
      if (multiSFTs == NULL)
        {
          XLAL_CHECK_NULL ( ( multiSFTs = XLALLoadMultiSFTs(SFTcatalog, minFreq, maxFreq) ) != NULL, XLAL_EFUNC );
        }
      else
        {
          // Otherwise, to minimize memory usage, loop over detectors, load SFTs and add them to generated SFTs
          for ( UINT4 X = 0; X < numDetectors; ++X )
            {
              SFTVector *loadedSFTs = NULL;
              XLAL_CHECK_NULL ( ( loadedSFTs = XLALLoadSFTs ( &multiSFTcatalog->data[X], minFreq, maxFreq) ) != NULL, XLAL_EFUNC );
              XLAL_CHECK_NULL ( XLALSFTVectorAdd ( multiSFTs->data[X], loadedSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );
              XLALDestroySFTVector ( loadedSFTs );
            } // for X < numDetectors
        } // if multiSFTs != NULL

    } // if loadSFTs

  // Normalise SFTs using either running median or assumed PSDs
  MultiPSDVector *runningMedian = XLALNormalizeMultiSFTVect ( multiSFTs, runningMedianWindow, assumeSqrtSX );
  XLAL_CHECK_NULL ( runningMedian != NULL, XLAL_EFUNC );

  // Calculate SFT noise weights from PSD
  common->noiseWeights = XLALComputeMultiNoiseWeights ( runningMedian, runningMedianWindow, 0 );
  XLAL_CHECK_NULL ( common->noiseWeights != NULL, XLAL_EFUNC );

  // Get detector states, with a timestamp shift of Tsft/2
  const REAL8 tOffset = 0.5 * Tsft;
  common->detectorStates = XLALGetMultiDetectorStates ( common->timestamps, &common->detectors, ephemerides, tOffset );
  XLAL_CHECK_NULL ( common->detectorStates != NULL, XLAL_EFUNC );

  // Save ephemerides and SSB precision
  common->ephemerides = ephemerides;
  common->SSBprec = extraParams->SSBprec;

  // Call the appropriate method function to setup their input data structures
  // - The method input data structures are expected to take ownership of the
  //   SFTs, which is why 'input->common' does not retain a pointer to them
  if ( input->demod != NULL )
    {
      XLAL_CHECK_NULL ( SetupFstatInput_Demod ( input->demod, common, multiSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  else if ( input->resamp != NULL )
    {
      XLAL_CHECK_NULL ( SetupFstatInput_Resamp ( input->resamp, common, multiSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  else
    {
      XLAL_ERROR_NULL ( XLAL_EFAILED, "Invalid FstatInput struct passed to %s()", __func__ );
    }
  multiSFTs = NULL;

  // Cleanup
  XLALDestroyMultiSFTCatalogView ( multiSFTcatalog );
  XLALDestroyMultiPSDVector ( runningMedian );

  return input;

} // XLALCreateFstatInput()

///
/// Returns the detector information stored in a \c FstatInput structure.
///
const MultiLALDetector*
XLALGetFstatInputDetectors ( const FstatInput* input    ///< [in] \c FstatInput structure.
                             )
{
  // Check input
  XLAL_CHECK_NULL ( input != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL ( input->common != NULL, XLAL_EINVAL, "'input' has not yet been set up" );

  return &input->common->detectors;

} // XLALGetFstatInputDetectors()

///
/// Returns the SFT timestamps stored in a \c FstatInput structure.
///
const MultiLIGOTimeGPSVector*
XLALGetFstatInputTimestamps ( const FstatInput* input	///< [in] \c FstatInput structure.
                              )
{
  // Check input
  XLAL_CHECK_NULL ( input != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL ( input->common != NULL, XLAL_EINVAL, "'input' has not yet been set up" );

  return input->common->timestamps;

} // XLALGetFstatInputTimestamps()

///
/// Returns the multi-detector noise weights stored in a \c FstatInput structure.
///
const MultiNoiseWeights*
XLALGetFstatInputNoiseWeights ( const FstatInput* input     ///< [in] \c FstatInput structure.
                                )
{
  // Check input
  XLAL_CHECK_NULL ( input != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL ( input->common != NULL, XLAL_EINVAL, "'input' has not yet been set up" );

  return input->common->noiseWeights;

} // XLALGetFstatInputNoiseWeights()

///
/// Returns the multi-detector state series stored in a \c FstatInput structure.
///
const MultiDetectorStateSeries*
XLALGetFstatInputDetectorStates ( const FstatInput* input	///< [in] \c FstatInput structure.
                                  )
{
  // Check input
  XLAL_CHECK_NULL ( input != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL ( input->common != NULL, XLAL_EINVAL, "'input' has not yet been set up" );

  return input->common->detectorStates;

} // XLALGetFstatInputDetectorStates()

///
/// Compute the \f$\mathcal{F}\f$-statistic over a band of frequencies.
///
int
XLALComputeFstat ( FstatResults **Fstats,	  	///< [in/out] Address of a pointer to a #FstatResults results structure; if \c NULL, allocate here.
                   FstatInput *input,		  	///< [in] Input data structure created by one of the setup functions.
                   const PulsarDopplerParams *doppler,  ///< [in] Doppler parameters, including starting frequency, at which to compute \f$2\mathcal{F}\f$
                   const REAL8 dFreq,	  		///< [in] Required spacing in frequency between each \f$\mathcal{F}\f$-statistic.
                   const UINT4 numFreqBins,		///< [in] Number of frequencies at which the \f$2\mathcal{F}\f$ are to be computed.
                   const FstatQuantities whatToCompute	///< [in] Bit-field of which \f$\mathcal{F}\f$-statistic quantities to compute.
                   )
{
  // Check input
  XLAL_CHECK(Fstats != NULL, XLAL_EFAULT);
  XLAL_CHECK(input != NULL, XLAL_EFAULT);
  XLAL_CHECK(input->common != NULL, XLAL_EINVAL, "'input' has not yet been set up");
  XLAL_CHECK(doppler != NULL, XLAL_EFAULT);
  XLAL_CHECK(doppler->asini >= 0, XLAL_EINVAL);
  XLAL_CHECK(dFreq > 0 || (numFreqBins == 1 && dFreq >= 0), XLAL_EINVAL);
  XLAL_CHECK(numFreqBins > 0, XLAL_EINVAL);
  XLAL_CHECK(0 < whatToCompute && whatToCompute < FSTATQ_LAST, XLAL_EINVAL);

  // Allocate results struct, if needed
  if ( (*Fstats) == NULL )
    {
      (*Fstats) = XLALCalloc ( 1, sizeof(**Fstats) );
      XLAL_CHECK ( (*Fstats) != NULL, XLAL_ENOMEM );
    }

  // Get constant pointer to common input data
  const FstatInput_Common *common = input->common;
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
          (*Fstats)->FaFb = XLALRealloc( (*Fstats)->FaFb, numFreqBins * sizeof((*Fstats)->FaFb[0]) );
          XLAL_CHECK ( (*Fstats)->FaFb != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->FaFb to length %u", numFreqBins );
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
              (*Fstats)->FaFbPerDet[X] = XLALRealloc ( (*Fstats)->FaFbPerDet[X], numFreqBins*sizeof((*Fstats)->FaFbPerDet[X][0]) );
              XLAL_CHECK( (*Fstats)->FaFbPerDet[X] != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->FaFbPerDet[%u] to length %u", X, numFreqBins );
            }
        }

      // Enlarge F-atoms per detector arrays, and initialise to NULL
      if ( (whatToCompute & FSTATQ_ATOMS_PER_DET) && (moreFreqBins || moreDetectors) )
        {
          for ( UINT4 X = 0; X < numDetectors; ++X )
            {
              (*Fstats)->multiFatoms = XLALRealloc ( (*Fstats)->multiFatoms, numFreqBins*sizeof((*Fstats)->multiFatoms[0]) );
              XLAL_CHECK ( (*Fstats)->multiFatoms != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->multiFatoms to length %u", numFreqBins );

              // If more detectors are needed, destroy multi-F-atom vectors so they can be re-allocated later
              if ( moreDetectors )
                {
                  for ( UINT4 k = 0; k < numFreqBins; ++k )
                    {
                      XLALDestroyMultiFstatAtomVector ( (*Fstats)->multiFatoms[k] );
                      (*Fstats)->multiFatoms[k] = NULL;
                    }
                }
              else
                {
                  for ( UINT4 k = (*Fstats)->internalalloclen; k < numFreqBins; ++k ) {
                    (*Fstats)->multiFatoms[k] = NULL;
                  }
                }
            } // for X < numDetectors

        } // if Atoms_per_det to enlarge

      // Update allocated length of arrays
      (*Fstats)->internalalloclen = numFreqBins;

    } // if (moreFreqBins || moreDetectors)

  // Initialise result struct parameters
  (*Fstats)->doppler = *doppler;
  (*Fstats)->dFreq = dFreq;
  (*Fstats)->numFreqBins = numFreqBins;
  (*Fstats)->numDetectors = numDetectors;
  memset ( (*Fstats)->detectorNames, 0, sizeof((*Fstats)->detectorNames) );
  for (UINT4 X = 0; X < numDetectors; ++X) {
    strncpy ( (*Fstats)->detectorNames[X], common->detectors.sites[X].frDetector.prefix, 2 );
  }
  (*Fstats)->whatWasComputed = whatToCompute;

  // Call the appropriate method function to compute the F-statistic
  if ( input->demod != NULL )
    {
      XLAL_CHECK ( ComputeFstat_Demod(*Fstats, common, input->demod) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  else if ( input->resamp != NULL )
    {
      XLAL_CHECK ( ComputeFstat_Resamp(*Fstats, common, input->resamp) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  else
    {
      XLAL_ERROR(XLAL_EFAILED, "Invalid FstatInput struct passed to %s()", __func__);
    }

  return XLAL_SUCCESS;

} // XLALComputeFstat()

///
/// Free all memory associated with a \c FstatInput structure.
///
void
XLALDestroyFstatInput ( FstatInput* input	///< [in] \c FstatInput structure to be freed.
                        )
{
  if ( input == NULL ) {
    return;
  }

  if (input->common != NULL)
    {
      XLALDestroyMultiTimestamps ( input->common->timestamps );
      XLALDestroyMultiNoiseWeights ( input->common->noiseWeights );
      XLALDestroyMultiDetectorStateSeries ( input->common->detectorStates );
      XLALFree ( input->common );
    }
  if (input->demod != NULL)
    {
      DestroyFstatInput_Demod ( input->demod );
    }
  else if ( input->resamp != NULL )
    {
      DestroyFstatInput_Resamp ( input->resamp );
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
  XLALFree ( Fstats->FaFb );
  for ( UINT4 X = 0; X < PULSAR_MAX_DETECTORS; ++X )
    {
      XLALFree ( Fstats->twoFPerDet[X] );
      XLALFree ( Fstats->FaFbPerDet[X] );
      if ( Fstats->multiFatoms != NULL )
        {
          for ( UINT4 n = 0; n < Fstats->internalalloclen; ++n )
            {
              XLALDestroyMultiFstatAtomVector ( Fstats->multiFatoms[n] );
            }
          XLALFree ( Fstats->multiFatoms );
        }
    }

  XLALFree ( Fstats );

  return;
} // XLALDestroyFstatResults()

///
/// Add +4 to any multi-detector or per-detector 2F values computed by XLALComputeFstat().
/// This is for compatibility with programs which expect this normalisation if SFTs do not
/// contain noise, e.g. \c lalapps_ComputeFStatistic with the \c --SignalOnly option.
///
int
XLALAdd4ToFstatResults ( FstatResults* Fstats    ///< [in/out] #FstatResults structure.
                         )
{
  // Check input
  XLAL_CHECK( Fstats != NULL, XLAL_EFAULT );

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
/// Estimate the amplitude parameters of a pulsar CW signal, given its phase parameters,
/// constituent parts of the \f$\mathcal{F}\f$-statistic, and antenna pattern matrix.
///
/// \note Parameter-estimation based on large parts on Yousuke's notes and implemention (in CFSv1),
/// extended for error-estimation.
///
int
XLALEstimatePulsarAmplitudeParams ( PulsarCandidate *pulsarParams,	///< [in,out] Pulsar candidate parameters.
                                    const LIGOTimeGPS* FaFb_refTime,	///< [in] Reference time of \f$F_a\f$ and \f$F_b\f$, may differ from pulsar candidate reference time.
                                    const COMPLEX16 Fa,			///< [in] Complex \f$\mathcal{F}\f$-statistic amplitude \f$F_a\f$.
                                    const COMPLEX16 Fb,			///< [in] Complex \f$\mathcal{F}\f$-statistic amplitude \f$F_b\f$.
                                    const AntennaPatternMatrix *Mmunu	///< [in] Antenna pattern matrix \f$M_{\mu\nu}\f$.
                                    )
{

  REAL8 A1h, A2h, A3h, A4h;
  REAL8 Ad, Bd, Cd, Dd, Ed;
  REAL8 normAmu;
  REAL8 A1check, A2check, A3check, A4check;

  REAL8 Asq, Da, disc;
  REAL8 aPlus, aCross;
  REAL8 Ap2, Ac2;
  REAL8 beta;
  REAL8 phi0, psi;
  REAL8 b1, b2, b3;
  REAL8 h0, cosi;

  REAL8 cosphi0, sinphi0, cos2psi, sin2psi;

  REAL8 tolerance = LAL_REAL4_EPS;

  gsl_vector *x_mu, *A_Mu;
  gsl_matrix *M_Mu_Nu;
  gsl_matrix *Jh_Mu_nu;
  gsl_permutation *permh;
  gsl_matrix *tmp, *tmp2;
  int signum;

  XLAL_CHECK ( pulsarParams != NULL, XLAL_EFAULT );
  XLAL_CHECK ( FaFb_refTime != NULL, XLAL_EFAULT );
  XLAL_CHECK ( isfinite(creal(Fa)) && isfinite(cimag(Fb)), XLAL_EDOM,
               "Fa = (%g, %g) is not finite", creal(Fa), cimag(Fa) );
  XLAL_CHECK ( isfinite(creal(Fb)) && isfinite(cimag(Fb)), XLAL_EDOM,
               "Fb = (%g, %g) is not finite", creal(Fb), cimag(Fb) );
  XLAL_CHECK ( Mmunu != NULL, XLAL_EFAULT );

  Ad = Mmunu->Ad;
  Bd = Mmunu->Bd;
  Cd = Mmunu->Cd;
  Ed = Mmunu->Ed;
  Dd = Ad * Bd - Cd * Cd - Ed * Ed;

  normAmu = 2.0 / sqrt(2.0 * Mmunu->Sinv_Tsft); // generally *very* small!!

  // ----- GSL memory allocation -----
  XLAL_CHECK ( ( x_mu = gsl_vector_calloc (4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( A_Mu = gsl_vector_calloc (4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( M_Mu_Nu = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( Jh_Mu_nu = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );

  XLAL_CHECK ( ( permh = gsl_permutation_calloc ( 4 ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( tmp = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( tmp2 = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );

  // ----- fill vector x_mu
  gsl_vector_set (x_mu, 0,   creal(Fa) );       // x_1
  gsl_vector_set (x_mu, 1,   creal(Fb) );       // x_2
  gsl_vector_set (x_mu, 2, - cimag(Fa) );       // x_3
  gsl_vector_set (x_mu, 3, - cimag(Fb) );       // x_4

  // ----- fill matrix M^{mu,nu} [symmetric: use UPPER HALF ONLY!!]
  gsl_matrix_set (M_Mu_Nu, 0, 0,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 1, 1,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 0, 1, - Cd / Dd );

  gsl_matrix_set (M_Mu_Nu, 0, 3, - Ed / Dd );
  gsl_matrix_set (M_Mu_Nu, 1, 2,   Ed / Dd );

  gsl_matrix_set (M_Mu_Nu, 2, 2,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 3, 3,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 2, 3, - Cd / Dd );

  // get (un-normalized) MLE's for amplitudes A^mu  = M^{mu,nu} x_nu

  /* GSL-doc: int gsl_blas_dsymv (CBLAS_UPLO_t Uplo, double alpha, const gsl_matrix * A,
   *                              const gsl_vector * x, double beta, gsl_vector * y )
   *
   * compute the matrix-vector product and sum: y = alpha A x + beta y
   * for the symmetric matrix A. Since the matrix A is symmetric only its
   * upper half or lower half need to be stored. When Uplo is CblasUpper
   * then the upper triangle and diagonal of A are used, and when Uplo
   * is CblasLower then the lower triangle and diagonal of A are used.
   */
  XLAL_CHECK( gsl_blas_dsymv (CblasUpper, 1.0, M_Mu_Nu, x_mu, 0.0, A_Mu) == 0, XLAL_EERR );

  A1h = gsl_vector_get ( A_Mu, 0 );
  A2h = gsl_vector_get ( A_Mu, 1 );
  A3h = gsl_vector_get ( A_Mu, 2 );
  A4h = gsl_vector_get ( A_Mu, 3 );

  Asq = SQ(A1h) + SQ(A2h) + SQ(A3h) + SQ(A4h);
  Da = A1h * A4h - A2h * A3h;
  disc = sqrt ( SQ(Asq) - 4.0 * SQ(Da) );

  Ap2  = 0.5 * ( Asq + disc );
  aPlus = sqrt(Ap2);            // not yet normalized

  Ac2 = 0.5 * ( Asq - disc );
  aCross = sqrt( Ac2 );
  aCross *= MYSIGN ( Da );      // not yet normalized

  beta = aCross / aPlus;

  b1 =   A4h - beta * A1h;
  b2 =   A3h + beta * A2h;
  b3 = - A1h + beta * A4h ;

  psi  = 0.5 * atan ( b1 /  b2 );       // in [-pi/4,pi/4] (gauge used also by TDS)
  phi0 =       atan ( b2 / b3 );        // in [-pi/2,pi/2]

  // Fix remaining sign-ambiguity by checking sign of reconstructed A1
  A1check = aPlus * cos(phi0) * cos(2.0*psi) - aCross * sin(phi0) * sin(2*psi);
  if ( A1check * A1h <  0 )
    phi0 += LAL_PI;

  cosphi0 = cos(phi0);
  sinphi0 = sin(phi0);
  cos2psi = cos(2*psi);
  sin2psi = sin(2*psi);

  // check numerical consistency of estimated Amu and reconstructed
  A1check =   aPlus * cosphi0 * cos2psi - aCross * sinphi0 * sin2psi;
  A2check =   aPlus * cosphi0 * sin2psi + aCross * sinphi0 * cos2psi;
  A3check = - aPlus * sinphi0 * cos2psi - aCross * cosphi0 * sin2psi;
  A4check = - aPlus * sinphi0 * sin2psi + aCross * cosphi0 * cos2psi;

  if ( ( fabs( (A1check - A1h)/A1h ) > tolerance ) ||
       ( fabs( (A2check - A2h)/A2h ) > tolerance ) ||
       ( fabs( (A3check - A3h)/A3h ) > tolerance ) ||
       ( fabs( (A4check - A4h)/A4h ) > tolerance ) )
  {
    if ( lalDebugLevel )
      XLALPrintError ( "WARNING %s(): Difference between estimated and reconstructed Amu exceeds tolerance of %g\n",
                       __func__, tolerance );
  }

  // translate A_{+,x} into {h_0, cosi}
  h0 = aPlus + sqrt ( disc );  // not yet normalized !
  cosi = aCross / h0;

  // ========== Estimate the errors ==========

  // ----- compute derivatives \partial A^\mu / \partial B^\nu, where
  // we consider the output-variables B^\nu = (h0, cosi, phi0, psi)
  // where aPlus = 0.5 * h0 * (1 + cosi^2)  and aCross = h0 * cosi
  { // Ahat^mu is defined as A^mu with the replacements: A_+ --> A_x, and A_x --> h0
    REAL8 A1hat =   aCross * cosphi0 * cos2psi - h0 * sinphi0 * sin2psi;
    REAL8 A2hat =   aCross * cosphi0 * sin2psi + h0 * sinphi0 * cos2psi;
    REAL8 A3hat = - aCross * sinphi0 * cos2psi - h0 * cosphi0 * sin2psi;
    REAL8 A4hat = - aCross * sinphi0 * sin2psi + h0 * cosphi0 * cos2psi;

    // ----- A1 =   aPlus * cosphi0 * cos2psi - aCross * sinphi0 * sin2psi; -----
    gsl_matrix_set (Jh_Mu_nu, 0, 0,   A1h / h0 );       /* dA1/h0 */
    gsl_matrix_set (Jh_Mu_nu, 0, 1,   A1hat );          /* dA1/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 0, 2,   A3h );            /* dA1/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 0, 3, - 2.0 * A2h );      /* dA1/dpsi */

    // ----- A2 =   aPlus * cosphi0 * sin2psi + aCross * sinphi0 * cos2psi; -----
    gsl_matrix_set (Jh_Mu_nu, 1, 0,   A2h / h0 );       /* dA2/h0 */
    gsl_matrix_set (Jh_Mu_nu, 1, 1,   A2hat );          /* dA2/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 1, 2,   A4h );            /* dA2/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 1, 3,   2.0 * A1h );      /* dA2/dpsi */

    // ----- A3 = - aPlus * sinphi0 * cos2psi - aCross * cosphi0 * sin2psi; -----
    gsl_matrix_set (Jh_Mu_nu, 2, 0,   A3h / h0 );       /* dA3/h0 */
    gsl_matrix_set (Jh_Mu_nu, 2, 1,   A3hat );          /* dA3/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 2, 2, - A1h );            /* dA3/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 2, 3, - 2.0 * A4h );      /* dA3/dpsi */

    // ----- A4 = - aPlus * sinphi0 * sin2psi + aCross * cosphi0 * cos2psi; -----
    gsl_matrix_set (Jh_Mu_nu, 3, 0,   A4h / h0 );       /* dA4/h0 */
    gsl_matrix_set (Jh_Mu_nu, 3, 1,   A4hat );          /* dA4/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 3, 2, - A2h );            /* dA4/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 3, 3,   2.0 * A3h );      /* dA4/dpsi */
  }

  // ----- compute inverse matrices Jh^{-1} by LU-decomposition -----
  XLAL_CHECK( gsl_linalg_LU_decomp (Jh_Mu_nu, permh, &signum ) == 0, XLAL_EERR );

  // inverse matrix
  XLAL_CHECK( gsl_linalg_LU_invert (Jh_Mu_nu, permh, tmp ) == 0, XLAL_EERR );
  gsl_matrix_memcpy ( Jh_Mu_nu, tmp );

  // ----- compute Jh^-1 . Minv . (Jh^-1)^T -----

  /* GSL-doc: gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha,
   *                          const gsl_matrix *A, const gsl_matrix *B, double beta, gsl_matrix *C)
   * These functions compute the matrix-matrix product and sum
   * C = \alpha op(A) op(B) + \beta C
   * where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans
   * and similarly for the parameter TransB.
   */

  // first tmp = Minv . (Jh^-1)^T
  XLAL_CHECK( gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, M_Mu_Nu, Jh_Mu_nu, 0.0, tmp ) == 0, XLAL_EERR );
  // then J^-1 . tmp , store result in tmp2
  XLAL_CHECK( gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Jh_Mu_nu, tmp, 0.0, tmp2 ) == 0, XLAL_EERR );
  gsl_matrix_memcpy ( Jh_Mu_nu, tmp2 );

  // ===== debug-output resulting matrices =====
  // propagate initial-phase from Fstat-reference-time to refTime of Doppler-params
  XLAL_CHECK( XLALExtrapolatePulsarPhase ( &phi0, pulsarParams->Doppler.fkdot, pulsarParams->Doppler.refTime, phi0, *FaFb_refTime )
              == XLAL_SUCCESS, XLAL_EFUNC );

  if ( phi0 < 0 )             /* make sure phi0 in [0, 2*pi] */
    phi0 += LAL_TWOPI;
  phi0 = fmod ( phi0, LAL_TWOPI );

  // fill candidate-struct with the obtained signal-parameters and error-estimations
  pulsarParams->Amp.h0     = normAmu * h0;
  pulsarParams->Amp.cosi   = cosi;
  pulsarParams->Amp.phi0   = phi0;
  pulsarParams->Amp.psi    = psi;

  // read out principal estimation-errors from diagonal elements of inverse Fisher-matrix
  pulsarParams->dAmp.h0     = normAmu * sqrt( gsl_matrix_get (Jh_Mu_nu, 0, 0 ) );
  pulsarParams->dAmp.cosi   = sqrt( gsl_matrix_get (Jh_Mu_nu, 1, 1 ) );
  pulsarParams->dAmp.phi0   = sqrt( gsl_matrix_get (Jh_Mu_nu, 2, 2 ) );
  pulsarParams->dAmp.psi    = sqrt( gsl_matrix_get (Jh_Mu_nu, 3, 3 ) );
  // return also the full Amplitude-Fisher matrix: invert Jh_Mu_nu
  XLAL_CHECK( gsl_linalg_LU_decomp (Jh_Mu_nu, permh, &signum ) == 0, XLAL_EERR );
  XLAL_CHECK( gsl_linalg_LU_invert (Jh_Mu_nu, permh, tmp ) == 0, XLAL_EERR );
  pulsarParams->AmpFisherMatrix = tmp;

  // ----- free GSL memory -----
  gsl_vector_free ( x_mu );
  gsl_vector_free ( A_Mu );
  gsl_matrix_free ( M_Mu_Nu );
  gsl_matrix_free ( Jh_Mu_nu );
  gsl_permutation_free ( permh );
  gsl_matrix_free ( tmp2 );

  return XLAL_SUCCESS;

} // XLALEstimatePulsarAmplitudeParams()

///
/// Convert amplitude params from 'physical' coordinates \f$(h_0, \cos\iota, \psi, \phi_0)\f$ into
/// 'canonical' coordinates \f$A^\mu = (A_1, A_2, A_3, A_4)\f$. The equations can be found in
/// \cite JKS98 or \cite Prix07 Eq.(2).
///
int
XLALAmplitudeParams2Vect ( PulsarAmplitudeVect A_Mu,		///< [out] Canonical amplitude coordinates \f$A^\mu = (A_1, A_2, A_3, A_4)\f$.
                           const PulsarAmplitudeParams Amp	///< [in] Physical amplitude params \f$(h_0, \cos\iota, \psi, \phi_0)\f$.
                           )
{

  REAL8 aPlus = 0.5 * Amp.h0 * ( 1.0 + SQ(Amp.cosi) );
  REAL8 aCross = Amp.h0 * Amp.cosi;
  REAL8 cos2psi = cos ( 2.0 * Amp.psi );
  REAL8 sin2psi = sin ( 2.0 * Amp.psi );
  REAL8 cosphi0 = cos ( Amp.phi0 );
  REAL8 sinphi0 = sin ( Amp.phi0 );

  XLAL_CHECK( A_Mu != NULL, XLAL_EFAULT );

  A_Mu[0] =  aPlus * cos2psi * cosphi0 - aCross * sin2psi * sinphi0;
  A_Mu[1] =  aPlus * sin2psi * cosphi0 + aCross * cos2psi * sinphi0;
  A_Mu[2] = -aPlus * cos2psi * sinphi0 - aCross * sin2psi * cosphi0;
  A_Mu[3] = -aPlus * sin2psi * sinphi0 + aCross * cos2psi * cosphi0;

  return XLAL_SUCCESS;

} // XLALAmplitudeParams2Vect()

///
/// Compute amplitude params \f$(h_0, \cos\iota, \psi, \phi_0)\f$ from amplitude-vector \f$A^\mu = (A_1, A_2, A_3, A_4)\f$.
/// Adapted from algorithm in XLALEstimatePulsarAmplitudeParams().
///
int
XLALAmplitudeVect2Params ( PulsarAmplitudeParams *Amp,		///< [out] Physical amplitude params \f$(h_0, \cos\iota, \psi, \phi_0)\f$.
                           const PulsarAmplitudeVect A_Mu	///< [in] Canonical amplitude coordinates \f$A^\mu = (A_1, A_2, A_3, A_4)\f$.
                           )
{

  REAL8 h0Ret, cosiRet, psiRet, phi0Ret;

  REAL8 A1, A2, A3, A4, Asq, Da, disc;
  REAL8 Ap2, Ac2, aPlus, aCross;
  REAL8 beta, b1, b2, b3;

  XLAL_CHECK( A_Mu != NULL, XLAL_EFAULT );
  XLAL_CHECK( Amp != NULL, XLAL_EFAULT );

  A1 = A_Mu[0];
  A2 = A_Mu[1];
  A3 = A_Mu[2];
  A4 = A_Mu[3];

  Asq = SQ(A1) + SQ(A2) + SQ(A3) + SQ(A4);
  Da = A1 * A4 - A2 * A3;

  disc = sqrt ( SQ(Asq) - 4.0 * SQ(Da) );

  Ap2  = 0.5 * ( Asq + disc );
  aPlus = sqrt(Ap2);

  Ac2 = 0.5 * ( Asq - disc );
  aCross = MYSIGN(Da) * sqrt( Ac2 );

  beta = aCross / aPlus;

  b1 =   A4 - beta * A1;
  b2 =   A3 + beta * A2;
  b3 = - A1 + beta * A4 ;

  // amplitude params in LIGO conventions
  psiRet  = 0.5 * atan2 ( b1,  b2 );  /* [-pi/2,pi/2] */
  phi0Ret =       atan2 ( b2,  b3 );  /* [-pi, pi] */

  // Fix remaining sign-ambiguity by checking sign of reconstructed A1
  REAL8 A1check = aPlus * cos(phi0Ret) * cos(2.0*psiRet) - aCross * sin(phi0Ret) * sin(2*psiRet);
  if ( A1check * A1 < 0 ) {
    phi0Ret += LAL_PI;
  }

  h0Ret = aPlus + sqrt ( disc );
  cosiRet = aCross / h0Ret;

  // make unique by fixing the gauge to be psi in [-pi/4, pi/4], phi0 in [0, 2*pi]
  while ( psiRet > LAL_PI_4 ) {
    psiRet  -= LAL_PI_2;
    phi0Ret -= LAL_PI;
  }
  while ( psiRet < - LAL_PI_4 ) {
    psiRet  += LAL_PI_2;
    phi0Ret += LAL_PI;
  }
  while ( phi0Ret < 0 ) {
    phi0Ret += LAL_TWOPI;
  }

  while ( phi0Ret > LAL_TWOPI ) {
    phi0Ret -= LAL_TWOPI;
  }

  // Return final answer
  Amp->h0   = h0Ret;
  Amp->cosi = cosiRet;
  Amp->psi  = psiRet;
  Amp->phi0 = phi0Ret;

  return XLAL_SUCCESS;

}


///
/// Compute single-or multi-IFO Fstat from multi-IFO 'atoms'
///
REAL8
XLALComputeFstatFromAtoms ( const MultiFstatAtomVector *multiFstatAtoms,   ///< [in] Multi-detector atoms
                            const INT4                 X                   ///< [in] Detector number, give -1 for multi-Fstat
                            )
{
  // ----- check input parameters and report errors
  XLAL_CHECK_REAL8 ( multiFstatAtoms && multiFstatAtoms->data && multiFstatAtoms->data[0]->data, XLAL_EFAULT, "Empty pointer as input parameter." );
  XLAL_CHECK_REAL8 ( multiFstatAtoms->length > 0, XLAL_EBADLEN, "Input MultiFstatAtomVector has zero length. (i.e., no detectors)" );
  XLAL_CHECK_REAL8 ( X >= -1, XLAL_EDOM, "Invalid detector number X=%d. Only nonnegative numbers, or -1 for multi-F, are allowed.", X );
  XLAL_CHECK_REAL8 ( ( X < 0 ) || ( (UINT4)(X) <= multiFstatAtoms->length-1 ), XLAL_EDOM, "Requested X=%d, but FstatAtoms only have length %d.", X, multiFstatAtoms->length );

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
  REAL8 mmatrixA = 0.0, mmatrixB = 0.0, mmatrixC = 0.0;
  REAL8 F = 0.0;
  COMPLEX8 Fa, Fb;
  Fa = 0.0;
  Fb = 0.0;

  for (Y = Ystart; Y <= Yend; Y++ ) /* loop through detectors */
    {
      UINT4 alpha, numSFTs;
      numSFTs = multiFstatAtoms->data[Y]->length;
      XLAL_CHECK_REAL8 ( numSFTs > 0, XLAL_EDOM, "Input FstatAtomVector has zero length. (i.e., no timestamps for detector X=%d)", Y );

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
  REAL8 Dinv = 1.0 / ( mmatrixA * mmatrixB - SQ(mmatrixC) );

  F = ComputeFstatFromFaFb ( Fa, Fb, mmatrixA, mmatrixB, mmatrixC, 0, Dinv );

  return F;

} // XLALComputeFstatFromAtoms()


///
/// Simple helper function which computes \f$\mathcal{F}\f$ from given \f$F_a\f$ and \f$F_b\f$, and antenna-pattern
/// coefficients \f$(A,B,C,E)\f$ with inverse determinant \f$\text{Dinv} = 1/D\f$ where \f$D = A * B - C^2 - E^2\f$.
///
static REAL8
ComputeFstatFromFaFb ( COMPLEX16 Fa, COMPLEX16 Fb, REAL8 A, REAL8 B, REAL8 C, REAL8 E, REAL8 Dinv )
{
  REAL8 Fa_re = creal(Fa);
  REAL8 Fa_im = cimag(Fa);
  REAL8 Fb_re = creal(Fb);
  REAL8 Fb_im = cimag(Fb);

  REAL8 F = Dinv * (  B * ( SQ(Fa_re) + SQ(Fa_im) )
                      + A * ( SQ(Fb_re) + SQ(Fb_im) )
                      - 2.0 * C * (   Fa_re * Fb_re + Fa_im * Fb_im )
                      - 2.0 * E * ( - Fa_re * Fb_im + Fa_im * Fb_re )		// nonzero only in RAA case where Ed!=0
                      );
  return F;

} // ComputeFstatFromFaFb()

///
/// Provide human-readable names for the different \f$\mathcal{F}\f$-statistic method variants in #FstatMethodType.
///
const CHAR *
XLALGetFstatMethodName ( FstatMethodType i )
{
  XLAL_CHECK_NULL ( (i > FMETHOD_START) && (i < FMETHOD_END) && FstatMethodNames[i].name!=NULL, XLAL_EDOM, "Invalid FstatMethodType = %d\n", i );
  return FstatMethodNames[i].name;
} // XLALGetFstatMethodName()

///
/// Return pointer to a static help string enumerating all (available) #FstatMethodType options.
/// Also indicates which is the (guessed) available 'best' method available.
///
const CHAR *
XLALFstatMethodHelpString ( void )
{
  static int firstCall = 1;
  static CHAR helpstr[1024];
  if ( firstCall )
    {
      CHAR buf[1024];
      strncpy (helpstr, "Available methods: (", sizeof(helpstr));
      UINT4 len = strlen(helpstr);
      const CHAR *separator = "";
      for (int i = FMETHOD_START + 1; i < FMETHOD_END; i++ )
        {
          if ( ! FstatMethodNames[i].available ) {
            continue;
          }
          snprintf ( buf, sizeof(buf), "%s%s", separator, FstatMethodNames[i].name );
          separator="|";
          if ( i == FMETHOD_DEMOD_BEST ) {
            strncat ( buf, "=DemodBest", sizeof(buf) - strlen(buf) - 1 );
          }
          if ( i == FMETHOD_RESAMP_BEST ) {
            strncat ( buf, "=ResampBest", sizeof(buf) - strlen(buf) - 1 );
          }
          len += strlen(buf);
          XLAL_CHECK_NULL ( len < sizeof(helpstr), XLAL_EBADLEN, "FstatMethod help-string exceeds buffer length (%d)\n", sizeof(helpstr) );
          strcat ( helpstr, buf );
        } // for i < FMETHOD_LAST

      strcat(helpstr, ") ");
      firstCall = 0;

    } // if firstCall

  return helpstr;
} // XLALFstatMethodHelpString()

///
/// Parse a given string into an #FstatMethodType number if valid and available,
/// return error otherwise.
///
int
XLALParseFstatMethodString ( FstatMethodType *Fmethod, 	//!< [out] Parsed #FstatMethodType enum
                             const char *s		//!< [in] String to parse
                             )
{
  XLAL_CHECK ( s != NULL, XLAL_EINVAL );
  XLAL_CHECK ( Fmethod != NULL, XLAL_EINVAL );

  // handle special user-input strings to select respective (guessed) best method
  if ( strcmp ( s, "DemodBest" ) == 0 )
    {
      (*Fmethod) = FMETHOD_DEMOD_BEST;
      return XLAL_SUCCESS;
    }
  if ( strcmp ( s, "ResampBest" ) == 0 )
    {
      (*Fmethod) = FMETHOD_RESAMP_BEST;
      return XLAL_SUCCESS;
    }

  // find matching FstatMethod string
  for (int i = FMETHOD_START + 1; i < FMETHOD_END; i++ )
    {
      if ( (FstatMethodNames[i].name != NULL) && (strcmp ( s, FstatMethodNames[i].name ) == 0) )
        {
          if ( FstatMethodNames[i].available )
            {
              (*Fmethod) = i;
              return XLAL_SUCCESS;
            }
          else
            {
              XLAL_ERROR ( XLAL_EINVAL, "Chosen FstatMethod '%s' valid but unavailable in this binary\n", s );
            }
        } // if found matching FstatMethod
    } // for i < FMETHOD_LAST

  XLAL_ERROR ( XLAL_EINVAL, "Unknown FstatMethod '%s'\n", s );

} // XLALParseFstatMethodString()
