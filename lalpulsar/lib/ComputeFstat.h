//
// Copyright (C) 2012--2014 David Keitel, Bernd Machenschalk, Reinhard Prix, Karl Wette
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
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301  USA
//

#ifndef _COMPUTEFSTAT_H
#define _COMPUTEFSTAT_H

#include <lal/LALStdlib.h>
#include <lal/UserInputParse.h>
#include <lal/PulsarDataTypes.h>
#include <lal/LALComputeAM.h>
#include <lal/LALComputeAM.h>
#include <lal/SSBtimes.h>
#include <lal/CWMakeFakeData.h>
#include <lal/LFTandTSutils.h>

#ifdef  __cplusplus
extern "C" {
#endif

///
/// \defgroup ComputeFstat_h Header ComputeFstat.h
/// \ingroup lalpulsar_coh
/// \authors Badri Krishnan, Bernd Machenschalk, Chris Messenger, David Keitel, Holger Pletsch,
///          John T. Whelan, Karl Wette, Pinkesh Patel, Reinhard Prix, Xavier Siemens
///
/// \brief The \f$\mathcal{F}\f$-statistic.
///
/// This module provides a API for computing the \f$\mathcal{F}\f$-statistic \cite JKS98 , using
/// various different methods.  All data required to compute the \f$\mathcal{F}\f$-statistic are
/// contained in the opaque structure \c FstatInput, which is shared by all methods. A function
/// XLALCreateFstatInput() is provided for creating an \c FstatInput structure configured
/// for the particular method.  The \c FstatInput structure is passed to the function
/// XLALComputeFstat(), which computes the \f$\mathcal{F}\f$-statistic using the chosen method, and
/// fills a \c FstatResults structure with the results.
///
/// \note The \f$\mathcal{F}\f$-statistic method codes are partly descended from earlier
/// implementations found in:
/// - <tt>LALDemod.[ch]</tt> by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve
///   Berukoff, Xavier Siemens, Bruce Allen
/// - <tt>ComputeSky.[ch]</tt> by Jolien Creighton, Reinhard Prix, Steve Berukoff
/// - <tt>LALComputeAM.[ch]</tt> by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve
///   Berukoff, Xavier Siemens
/// - <tt>ComputeFStatistic_resamp.c</tt> by Pinkesh Patel, Xavier Siemens, Reinhard Prix, Iraj
///   Gholami, Yousuke Itoh, Maria Alessandra Papa
///

/// @{

/// default maximum allowed F-stat mismatch from SFTs being too long,
/// to be used in XLALFstatCheckSFTLengthMismatch()
/// - this value allows a search using 1800-second SFTs for isolated CW signals with frequencies <~ 2054 Hz,
///   while preventing a search using 1800-second SFTs for a Sco X-1-like binary CW signal
///   (which would require <~ 200-second SFTs)
#define DEFAULT_MAX_MISMATCH_FROM_SFT_LENGTH 0.05

///
/// XLALComputeFstat() input data structure. Encapsulates all data, buffers, etc. used by the
/// \f$\mathcal{F}\f$-statistic methods.
///
typedef struct tagFstatInput FstatInput;

///
/// A vector of XLALComputeFstat() input data structures, for e.g. computing the
/// \f$\mathcal{F}\f$-statistic for multiple segments.
///
typedef struct tagFstatInputVector {
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(FstatInputVector, FstatInput*, data, UINT4, length));
#endif // SWIG
  UINT4 length;                     ///< Number of elements in array.
  FstatInput **data;                ///< Pointer to the data array.
} FstatInputVector;

///
/// Bit-field of \f$\mathcal{F}\f$-statistic quantities which can be computed by XLALComputeFstat().
/// Not all options are supported by all \f$\mathcal{F}\f$-statistic methods.
///
typedef enum tagFstatQuantities {
  FSTATQ_NONE           = 0x00,         ///< Do not compute \f$\mathcal{F}\f$-statistic, still compute buffered quantities
  FSTATQ_2F             = 0x01,         ///< Compute multi-detector \f$2\mathcal{F}\f$.
  FSTATQ_FAFB           = 0x02,         ///< Compute multi-detector \f$F_a\f$ and \f$F_b\f$.
  FSTATQ_2F_PER_DET     = 0x04,         ///< Compute \f$2\mathcal{F}\f$ for each detector.
  FSTATQ_FAFB_PER_DET   = 0x08,         ///< Compute \f$F_a\f$ and \f$F_b\f$ for each detector.
  FSTATQ_ATOMS_PER_DET  = 0x10,         ///< Compute per-SFT \f$\mathcal{F}\f$-statistic atoms for each detector (\a Demod only).
  FSTATQ_2F_CUDA        = 0x20,         ///< Compute multi-detector \f$2\mathcal{F}\f$, store results on CUDA GPU (CUDA implementation of \a Resamp only).
  FSTATQ_LAST           = 0x80
} FstatQuantities;

///
/// Different methods available to compute the \f$\mathcal{F}\f$-statistic, falling into two broad classes:
/// * \a Demod: Dirichlet kernel-based demodulation \cite Williams1999 .
/// * \a Resamp: FFT-based resampling \cite JKS98 \cite Prix2022 \cite DunnEtAl2022 .
///
typedef enum tagFstatMethodType {

  /// \cond DONT_DOXYGEN
  FMETHOD_START = 1,
  /// \endcond

  FMETHOD_DEMOD_GENERIC,	///< \a Demod: generic C hotloop, works for any number of Dirichlet kernel terms \f$\text{Dterms}\f$
  FMETHOD_DEMOD_OPTC,		///< \a Demod: gptimized C hotloop using Akos' algorithm, only works for \f$\text{Dterms} \lesssim 20\f$
  FMETHOD_DEMOD_ALTIVEC,	///< \a Demod: Altivec hotloop variant, uses fixed \f$\text{Dterms} = 8\f$
  FMETHOD_DEMOD_SSE,		///< \a Demod: SSE hotloop with precalc divisors, uses fixed \f$\text{Dterms} = 8\f$
  FMETHOD_DEMOD_BEST,		///< \a Demod: best guess of the fastest available hotloop

  FMETHOD_RESAMP_GENERIC,	///< \a Resamp: generic implementation \cite Prix2022
  FMETHOD_RESAMP_CUDA,		///< \a Resamp: CUDA resampling \cite DunnEtAl2022
  FMETHOD_RESAMP_BEST,		///< \a Resamp: best guess of the fastest available implementation

  /// \cond DONT_DOXYGEN
  FMETHOD_END
  /// \endcond

} FstatMethodType;

///
/// Struct of optional 'advanced level' and (potentially method-specific) arguments to be passed to the
/// \f$\mathcal{F}\f$-statistic setup function XLALCreateFstatInput().
///
typedef struct tagFstatOptionalArgs {
  UINT4 randSeed;			///< Random-number seed value used in case of fake Gaussian noise generation (\c injectSqrtSX)
  SSBprecision SSBprec;			///< Barycentric transformation precision.
  UINT4 Dterms;                  	///< Number of Dirichlet kernel terms, used by some \a Demod methods; see \c FstatMethodType.
  UINT4 runningMedianWindow;	  	///< If SFT noise weights are calculated from the SFTs, the running median window length to use.
  FstatMethodType FstatMethod;	  	///< Method to use for computing the \f$\mathcal{F}\f$-statistic.
  PulsarParamsVector *injectSources;	///< Vector of parameters of CW signals to simulate and inject.
  MultiNoiseFloor *injectSqrtSX;  	///< Single-sided PSD values for fake Gaussian noise to be added to SFT data.
  MultiNoiseFloor *assumeSqrtSX;  	///< Single-sided PSD values to be used for computing SFT noise weights instead of from a running median of the SFTs themselves.
  FstatInput *prevInput;		///< An \c FstatInput structure from a previous call to XLALCreateFstatInput(); may contain common workspace data than can be re-used to save memory.
  BOOLEAN collectTiming;		///< a flag to turn on/off the collection of F-stat-method-specific timing-data
  BOOLEAN resampFFTPowerOf2;		///< \a Resamp: round up FFT lengths to next power of 2; see \c FstatMethodType.
  REAL8 allowedMismatchFromSFTLength;      ///<  Optional override for XLALFstatCheckSFTLengthMismatch().
  REAL8 sourceDeltaT;			///< Optional source-frame sampling period for XLALCWMakeFakeData(); if zero, use the previous internal defaults.
} FstatOptionalArgs;

///
/// Global initializer for setting \c FstatOptionalArgs to default values
///
extern const FstatOptionalArgs FstatOptionalArgsDefaults;

///
/// An \f$\mathcal{F}\f$-statistic 'atom', i.e. the elementary per-SFT quantities required to compute the
/// \f$\mathcal{F}\f$-statistic, for one detector X.
///
typedef struct tagFstatAtom {
  UINT4 timestamp;                      ///< SFT GPS timestamp \f$t_i\f$ in seconds.
  REAL4 a2_alpha;                       ///< Antenna-pattern factor \f$a^2(X,t_i)\f$.
  REAL4 b2_alpha;                       ///< Antenna-pattern factor \f$b^2(X,t_i)\f$.
  REAL4 ab_alpha;                       ///< Antenna-pattern factor \f$a*b(X,t_i)\f$.
  COMPLEX8 Fa_alpha;                    ///< \f$Fa^X(t_i)\f$.
  COMPLEX8 Fb_alpha;                    ///< \f$Fb^X(t_i)\f$.
} FstatAtom;

///
/// A vector of \f$\mathcal{F}\f$-statistic 'atoms', i.e. all per-SFT quantities required to compute
/// the \f$\mathcal{F}\f$-statistic, for one detector X.
///
typedef struct tagFstatAtomVector {
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(FstatAtomVector, FstatAtom, data, UINT4, length));
#endif // SWIG
  UINT4 length;                         ///< Number of per-SFT 'atoms'.
  FstatAtom *data;                      ///< Array of \c FstatAtom pointers of given length.
  UINT4 TAtom;                          ///< Time-baseline of 'atoms', typically \f$T_{\mathrm{sft}}\f$.
} FstatAtomVector;

///
/// A multi-detector vector of \c FstatAtomVector.
///
typedef struct tagMultiFstatAtomVector {
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(MultiFstatAtomVector, FstatAtomVector*, data, UINT4, length));
#endif // SWIG
  UINT4 length;                         ///< Number of detectors.
  FstatAtomVector **data;               ///< Array of \c FstatAtomVector pointers, one for each detector X.
} MultiFstatAtomVector;

///
/// XLALComputeFstat() computed results structure.
///
#ifdef SWIG // SWIG interface directives
SWIGLAL(IGNORE_MEMBERS(tagFstatResults, internalalloclen));
SWIGLAL(ARRAY_MULTIPLE_LENGTHS(tagFstatResults, numFreqBins, numDetectors));
#endif // SWIG
typedef struct tagFstatResults {

  /// Doppler parameters, including the starting frequency, at which the \f$2\mathcal{F}\f$ were
  /// computed.
  PulsarDopplerParams doppler;

  /// For performance reasons the global phase of all returned 'Fa' and 'Fb' quantities (#Fa,#Fb,#FaPerDet,#FbPerDet, #multiFatoms),
  /// refers to this 'phase reference time' instead of the (#doppler).refTime.
  /// Use this to compute initial phase of signals at doppler.refTime, or if you need the correct global Fa,Fb-phase!
  LIGOTimeGPS refTimePhase;

  /// Spacing in frequency between each computed \f$\mathcal{F}\f$-statistic.
  REAL8 dFreq;

  /// Number of frequencies at which the \f$2\mathcal{F}\f$ were computed.
  UINT4 numFreqBins;

  /// Number of detectors over which the \f$2\mathcal{F}\f$ were computed.  Valid range is 1 to
  /// #PULSAR_MAX_DETECTORS.
  UINT4 numDetectors;

  /// Names of detectors over which the \f$2\mathcal{F}\f$ were computed.  Valid range is 1 to
  /// #PULSAR_MAX_DETECTORS.
  CHAR detectorNames[PULSAR_MAX_DETECTORS][3];

  /// Antenna pattern matrix \f$M_{\mu\nu}\f$, used in computing \f$2\mathcal{F}\f$.
  AntennaPatternMatrix Mmunu;

  /// Per detector antenna pattern matrix \f$M_{\mu\nu}^X\f$, used in computing \f$2\mathcal{F}^X\f$.
  AntennaPatternMatrix MmunuX[PULSAR_MAX_DETECTORS];

  /// Bit-field of which \f$\mathcal{F}\f$-statistic quantities were computed.
  FstatQuantities whatWasComputed;

  /// If #whatWasComputed & FSTATQ_2F is true, the multi-detector \f$2\mathcal{F}\f$ values computed
  /// at #numFreqBins frequencies spaced #dFreq apart.  This array should not be accessed if
  /// #whatWasComputed & FSTATQ_2F is false.
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(FstatResults, REAL4, twoF, UINT4, numFreqBins));
#endif // SWIG
  REAL4 *twoF;

#ifndef SWIG // exclude from SWIG interface
  /// If #whatWasComputed & FSTATQ_2F_CUDA is true, the multi-detector \f$2\mathcal{F}\f$ values
  /// as for #twoF, but stored in CUDA device memory.  This array should not be accessed if
  /// #whatWasComputed & FSTATQ_2F_CUDA is false.
  REAL4 *twoF_CUDA;
#endif

  /// If #whatWasComputed & FSTATQ_PARTS is true, the multi-detector \f$F_a\f$ and \f$F_b\f$
  /// computed at #numFreqBins frequencies spaced #dFreq apart.  This array should not be accessed
  /// if #whatWasComputed & FSTATQ_PARTS is false.
  /// \note global phase refers to #refTimePhase, not (#doppler).refTime
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(FstatResults, COMPLEX8, Fa, UINT4, numFreqBins));
  SWIGLAL(ARRAY_1D(FstatResults, COMPLEX8, Fb, UINT4, numFreqBins));
#endif // SWIG
  COMPLEX8 *Fa;
  COMPLEX8 *Fb;

  /// If #whatWasComputed & FSTATQ_2F_PER_DET is true, the \f$2\mathcal{F}\f$ values computed at
  /// #numFreqBins frequencies spaced #dFreq apart, and for #numDetectors detectors.  Only the first
  /// #numDetectors entries will be valid.  This array should not be accessed if #whatWasComputed &
  /// FSTATQ_2F_PER_DET is false.
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D_PTR_1D(FstatResults, REAL4, twoFPerDet, UINT4, numDetectors, numFreqBins));
#endif // SWIG
  REAL4 *twoFPerDet[PULSAR_MAX_DETECTORS];

  /// If #whatWasComputed & FSTATQ_PARTS_PER_DET is true, the \f$F_a\f$ and \f$F_b\f$ values
  /// computed at #numFreqBins frequencies spaced #dFreq apart, and for #numDetectors detectors.
  /// This array should not be accessed if #whatWasComputed & FSTATQ_PARTS_PER_DET is false.
  /// \note global phase refers to #refTimePhase, not (#doppler).refTime
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D_PTR_1D(FstatResults, COMPLEX8, FaPerDet, UINT4, numDetectors, numFreqBins));
  SWIGLAL(ARRAY_1D_PTR_1D(FstatResults, COMPLEX8, FbPerDet, UINT4, numDetectors, numFreqBins));
#endif // SWIG
  COMPLEX8 *FaPerDet[PULSAR_MAX_DETECTORS];
  COMPLEX8 *FbPerDet[PULSAR_MAX_DETECTORS];

  /// If #whatWasComputed & FSTATQ_ATOMS_PER_DET is true, the per-SFT \f$\mathcal{F}\f$-statistic
  /// multi-atoms computed at #numFreqBins frequencies spaced #dFreq apart.  This array should not
  /// be accessed if #whatWasComputed & FSTATQ_ATOMS_PER_DET is false.
  /// \note global phase of atoms refers to #refTimePhase, not (#doppler).refTime
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(FstatResults, MultiFstatAtomVector*, multiFatoms, UINT4, numFreqBins));
#endif // SWIG
  MultiFstatAtomVector** multiFatoms;

  /// \cond DONT_DOXYGEN
  UINT4 internalalloclen;
  /// \endcond

} FstatResults;

/// Generic F-stat timing coefficients (times in seconds)
/// [see https://dcc.ligo.org/LIGO-T1600531-v4 for details]
/// tauF_eff = tauF_core + b * tauF_buffer
/// where 'b' denotes the fraction of times per call to XLALComputeFstat() the buffer needed to be recomputed
/// ==> b = NBufferMisses / NCalls\n"
///
/// Note: All timing numbers are averaged over repeated calls to XLALComputeFstat(for fixed-setup)
///
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagFstatTimingGeneric, help));
#endif /* SWIG */
typedef struct tagFstatTimingGeneric
{
  REAL4 tauF_eff;	//< effective total time per F-stat call to XLALComputeFstat() per detector per output frequency
  REAL4 tauF_core;	//< core time per output frequency bin per detector, excluding time to compute the buffer
  REAL4 tauF_buffer;	//< time per detector per frequency bin to re-compute all buffered quantities once
  UINT4 NFbin;		//< (average over F-stat calls) number of F-stat output frequency bins
  UINT4 Ndet;		//< number of detectors
  REAL4 NCalls;		//< number of F-stat calls we average over
  REAL4 NBufferMisses;	//< number of times the buffer needed to be recomputed
  const char *help;	//< (static) string documenting the generic F-stat timing values
} FstatTimingGeneric;

#define TIMING_MODEL_MAX_VARS 10
/// Struct to carry the \f$\mathcal{F}\f$-statistic method-specific timing *model* in terms of
/// a list of variable names and corresponding REAL4 values, including a help-string documenting
/// the timing model variables.
/// See https://dcc.ligo.org/LIGO-T1600531-v4 for a more detailed discussion of the F-stat timing model.
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagFstatTimingModel, help));
#endif /* SWIG */
typedef struct tagFstatTimingModel
{
  UINT4 numVariables;                           //< number of (method-specific) timing model variables
  const char* names[TIMING_MODEL_MAX_VARS];     //< array of (static) timing-model variable names
  REAL4 values[TIMING_MODEL_MAX_VARS];          //< array of timing-model variable values
  const char *help;                            	//< (static) help string documenting the (method-specific) Fstat timing model
} FstatTimingModel;

// ---------- API function prototypes ----------
REAL8 XLALFstatMaximumSFTLength ( const REAL8 maxFreq, const REAL8 binaryMaxAsini, const REAL8 binaryMinPeriod, const REAL8 mu_SFT
);
int XLALFstatCheckSFTLengthMismatch ( const REAL8 Tsft, const REAL8 maxFreq, const REAL8 binaryMaxAsini, const REAL8 binaryMinPeriod, const REAL8 allowedMismatch );

int XLALFstatMethodIsAvailable ( FstatMethodType method );
const CHAR *XLALFstatMethodName ( FstatMethodType method );
const UserChoices *XLALFstatMethodChoices ( void );

FstatInputVector* XLALCreateFstatInputVector ( const UINT4 length );
void XLALDestroyFstatInputVector ( FstatInputVector* input );
FstatAtomVector* XLALCreateFstatAtomVector ( const UINT4 length );
void XLALDestroyFstatAtomVector ( FstatAtomVector *atoms );
MultiFstatAtomVector* XLALCreateMultiFstatAtomVector ( const UINT4 length );
void XLALDestroyMultiFstatAtomVector ( MultiFstatAtomVector *atoms );

FstatInput *
XLALCreateFstatInput ( const SFTCatalog *SFTcatalog, const REAL8 minCoverFreq, const REAL8 maxCoverFreq, const REAL8 dFreq,
                       const EphemerisData *ephemerides, const FstatOptionalArgs *optionalArgs );

int XLALGetFstatInputSFTBand ( const FstatInput *input, REAL8 *minFreqFull, REAL8 *maxFreqFull );
const CHAR *XLALGetFstatInputMethodName ( const FstatInput* input );
const MultiLALDetector* XLALGetFstatInputDetectors ( const FstatInput* input );
const MultiLIGOTimeGPSVector* XLALGetFstatInputTimestamps ( const FstatInput* input );
MultiNoiseWeights* XLALGetFstatInputNoiseWeights ( const FstatInput* input );
const MultiDetectorStateSeries* XLALGetFstatInputDetectorStates ( const FstatInput* input );
int  XLALExtractResampledTimeseries ( MultiCOMPLEX8TimeSeries **multiTimeSeries_SRC_a, MultiCOMPLEX8TimeSeries **multiTimeSeries_SRC_b, FstatInput *input );
int XLALGetFstatTiming ( const FstatInput* input, FstatTimingGeneric *timingGeneric, FstatTimingModel *timingModel );
int XLALAppendFstatTiming2File ( const FstatInput* input, FILE *fp, BOOLEAN printHeader );
int XLALFstatInputTimeslice ( FstatInput** slice, const FstatInput* input, const LIGOTimeGPS *minStartGPS, const LIGOTimeGPS *maxStartGPS);

#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(FstatResults**, Fstats));
#endif
int XLALComputeFstat ( FstatResults **Fstats, FstatInput *input, const PulsarDopplerParams *doppler,
                       const UINT4 numFreqBins, const FstatQuantities whatToCompute );

void XLALDestroyFstatInput ( FstatInput* input );
void XLALDestroyFstatResults ( FstatResults* Fstats );

REAL4 XLALComputeFstatFromFaFb ( COMPLEX8 Fa, COMPLEX8 Fb, REAL4 A, REAL4 B, REAL4 C, REAL4 E, REAL4 Dinv );

REAL4 XLALComputeFstatFromAtoms ( const MultiFstatAtomVector *multiFstatAtoms, const INT4 X );

/// @}

#ifdef  __cplusplus
}
#endif

#endif // _COMPUTEFSTAT_H
