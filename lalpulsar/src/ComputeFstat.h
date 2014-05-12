//
// Copyright (C) 2012, 2013, 2014 David Keitel, Bernd Machenschalk, Reinhard Prix, Karl Wette
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

#ifndef _COMPUTEFSTAT_H
#define _COMPUTEFSTAT_H

#include <lal/LALStdlib.h>
#include <lal/PulsarDataTypes.h>
#include <lal/LALComputeAM.h>
#include <lal/ComplexAM.h>
#include <lal/SSBtimes.h>
#include <lal/CWMakeFakeData.h>

#ifdef  __cplusplus
extern "C" {
#endif

///
/// \defgroup ComputeFstat_h Header ComputeFstat.h
/// \ingroup pkg_pulsarCoh
/// \authors Badri Krishnan, Bernd Machenschalk, Chris Messenger, David Keitel, Holger Pletsch,
///          John T. Whelan, Karl Wette, Pinkesh Patel, Reinhard Prix, Xavier Siemens
///
/// \brief The \f$\mathcal{F}\f$-statistic.
///
/// This module provides a API for computing the \f$\mathcal{F}\f$-statistic \cite JKS98, using
/// various different algorithms.  All data required to compute the \f$\mathcal{F}\f$-statistic are
/// contained in the opaque structure \c FstatInput, which is shared by all algorithms. A
/// function <tt>XLALCreateFstatInput_...()</tt> is provided by each algorithm for creating an
/// \c FstatInput structure configured for the particular algorithm.  The \c FstatInput
/// structure is then passed to the function XLALSetupFstatInput(), which performs general
/// initialisation tasks.  Finally, the \c FstatInput structure is passed to the function
/// XLALComputeFstat(), which computes the \f$\mathcal{F}\f$-statistic using the chosen algorithm,
/// and fills a \c FstatResults structure with the results.
///
/// \note The \f$\mathcal{F}\f$-statistic algorithm codes are partly descended from earlier
/// implementations found in:
/// - LALDemod.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff,
///   Xavier Siemens, Bruce Allen
/// - ComputeSky.[ch] by Jolien Creighton, Reinhard Prix, Steve Berukoff
/// - LALComputeAM.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff,
///   Xavier Siemens
/// - ComputeFStatistic_resamp.c by Pinkesh Patel, Xavier Siemens, Reinhard Prix, Iraj Gholami,
///   Yousuke Itoh, Maria Alessandra Papa
///

// @{

///
/// XLALComputeFstat() input data structure. Encapsulates all data, buffers, etc. used by the
/// \f$\mathcal{F}\f$-statistic algorithms.
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
  UINT4 length;                         ///< Number of elements in array.
  FstatInput **data;                ///< Pointer to the data array.
} FstatInputVector;

///
/// Bit-field of \f$\mathcal{F}\f$-statistic quantities which can be computed by XLALComputeFstat().
/// Not all options are supported by all \f$\mathcal{F}\f$-statistic algorithms.
///
typedef enum tagFstatQuantities {
  FSTATQ_NONE           = 0x00,
  FSTATQ_2F             = 0x01,         ///< Compute multi-detector \f$2\mathcal{F}\f$.
  FSTATQ_FAFB           = 0x02,         ///< Compute multi-detector \f$F_a\f$ and \f$F_b\f$.
  FSTATQ_2F_PER_DET     = 0x04,         ///< Compute \f$2\mathcal{F}\f$ for each detector.
  FSTATQ_FAFB_PER_DET   = 0x08,         ///< Compute \f$F_a\f$ and \f$F_b\f$ for each detector.
  FSTATQ_ATOMS_PER_DET  = 0x10,         ///< Compute per-SFT \f$\mathcal{F}\f$-statistic atoms for each detector (demodulation only).
  FSTATQ_LAST           = 0x20
} FstatQuantities;

///
/// Different algorithms available to compute the F-statitistc, falling into two broad classes:
/// 'Demod' = Dirichel-Kernel based demodulation, following Williams&Schutz's method
/// 'Resamp' = FFT-based resampling
///
typedef enum tagFstatMethodType {
  FMETHOD_START = 1,		///< overall start marker (set to 1 to allow range-check without warnings)

  FMETHOD_DEMOD_START,		///< demod start marker
  FMETHOD_DEMOD_GENERIC,	///< Old 'generic' C hotloop, works for all values of Dterms
  FMETHOD_DEMOD_OPTC,		///< Optimized C hotloop using Akos' algorithm, only works for Dterms <~20
  FMETHOD_DEMOD_SSE,		///< SSE hotloop (with precalc divisors) (Dterms=8)
  FMETHOD_DEMOD_ALTIVEC,	///< Altivec hotloop variant (Dterms=8)
  FMETHOD_DEMOD_END,		///< demod end marker

  FMETHOD_RESAMP_START,		///< resamp start marker
  FMETHOD_RESAMP_GENERIC,	///< 'generic' resampling implementation
  FMETHOD_RESAMP_END,		///< resamp end marker

  FMETHOD_END			///< overall end marker
} FstatMethodType;

/// Some helpful range macros
#define XLALFstatMethodClassIsDemod(x)  ( ((x) > FMETHOD_DEMOD_START )  && ((x) < FMETHOD_DEMOD_END ) )
#define XLALFstatMethodClassIsResamp(x) ( ((x) > FMETHOD_RESAMP_START ) && ((x) < FMETHOD_RESAMP_END ) )

///
/// Provide a 'best guess' heuristic as to which available demodulation hotloop variant will be fastest.
/// Can be used as a user default value.
///
extern const int FMETHOD_DEMOD_BEST;
extern const int FMETHOD_RESAMP_BEST;

const CHAR *XLALGetFstatMethodName ( FstatMethodType i);
const CHAR *XLALFstatMethodHelpString ( void );
int XLALParseFstatMethodString ( FstatMethodType *Fmethod, const char *s );

///
/// Complex \f$\mathcal{F}\f$-statistic amplitudes \f$F_a\f$ and \f$F_b\f$.
///
typedef struct tagFstatFaFb {
  COMPLEX16 Fa;                         ///< Complex amplitude \f$F_a\f$.
  COMPLEX16 Fb;                         ///< Complex amplitude \f$F_b\f$.
} FstatFaFb;

///
/// An \f$\mathcal{F}\f$-statistic 'atom', i.e. the elementary per-SFT quantities required to compute the
/// \f$\mathcal{F}\f$-statistic, for one detector X.
///
typedef struct tagFstatAtom {
  UINT4 timestamp;                      ///< SFT GPS timestamp \f$t_i\f$ in seconds.
  REAL8 a2_alpha;                       ///< Antenna-pattern factor \f$a^2(X,t_i)\f$.
  REAL8 b2_alpha;                       ///< Antenna-pattern factor \f$b^2(X,t_i)\f$.
  REAL8 ab_alpha;                       ///< Antenna-pattern factor \f$a*b(X,t_i)\f$.
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
  FstatAtom *data;                      ///< Array of #FstatAtom pointers of given length.
  UINT4 TAtom;                          ///< Time-baseline of 'atoms', typically \f$T_{\mathrm{sft}}\f$.
} FstatAtomVector;

///
/// A multi-detector vector of #FstatAtomVector.
///
typedef struct tagMultiFstatAtomVector {
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(MultiFstatAtomVector, FstatAtomVector*, data, UINT4, length));
#endif // SWIG
  UINT4 length;                         ///< Number of detectors.
  FstatAtomVector **data;               ///< Array of #FstatAtomVector pointers, one for each detector X.
} MultiFstatAtomVector;

///
/// XLALComputeFstat() computed results structure.
///
#ifdef SWIG // SWIG interface directives
SWIGLAL(IMMUTABLE_MEMBERS(tagFstatResults, internalalloclen));
SWIGLAL(ARRAY_MULTIPLE_LENGTHS(tagFstatResults, numFreqBins, numDetectors));
#endif // SWIG
typedef struct tagFstatResults {

  /// Doppler parameters, including the starting frequency, at which the \f$2\mathcal{F}\f$ were
  /// computed.
  PulsarDopplerParams doppler;

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

  /// Bit-field of which \f$\mathcal{F}\f$-statistic quantities were computed.
  FstatQuantities whatWasComputed;

  /// If #whatWasComputed & FSTATQ_2F is true, the multi-detector \f$2\mathcal{F}\f$ values computed
  /// at #numFreqBins frequencies spaced #dFreq apart.  This array should not be accessed if
  /// #whatWasComputed & FSTATQ_2F is false.
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(FstatResults, REAL4, twoF, UINT4, numFreqBins));
#endif // SWIG
  REAL4 *twoF;

  /// If #whatWasComputed & FSTATQ_PARTS is true, the multi-detector \f$F_a\f$ and \f$F_b\f$
  /// computed at #numFreqBins frequencies spaced #dFreq apart.  This array should not be accessed
  /// if #whatWasComputed & FSTATQ_PARTS is false.
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(FstatResults, FstatFaFb, FaFb, UINT4, numFreqBins));
#endif // SWIG
  FstatFaFb *FaFb;

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
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D_PTR_1D(FstatResults, FstatFaFb, FaFb, UINT4, numDetectors, numFreqBins));
#endif // SWIG
  FstatFaFb *FaFbPerDet[PULSAR_MAX_DETECTORS];

  /// If #whatWasComputed & FSTATQ_ATOMS_PER_DET is true, the per-SFT \f$\mathcal{F}\f$-statistic
  /// multi-atoms computed at #numFreqBins frequencies spaced #dFreq apart.  This array should not
  /// be accessed if #whatWasComputed & FSTATQ_ATOMS_PER_DET is false.
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(FstatResults, MultiFstatAtomVector*, multiFatoms, UINT4, numFreqBins));
#endif // SWIG
  MultiFstatAtomVector** multiFatoms;

  /// \cond DONT_DOXYGEN
  UINT4 internalalloclen;
  /// \endcond

} FstatResults;

///
/// Create a #FstatInputVector of the given length.
///
FstatInputVector*
XLALCreateFstatInputVector(
  const UINT4 length                            ///< [in] Length of the #FstatInputVector.
  );

///
/// Free all memory associated with a #FstatInputVector structure.
///
void
XLALDestroyFstatInputVector(
  FstatInputVector* input                   ///< [in] #FstatInputVector structure to be freed.
  );

///
/// Create a #FstatAtomVector of the given length.
///
FstatAtomVector*
XLALCreateFstatAtomVector(
  const UINT4 length                            ///< [in] Length of the #FstatAtomVector.
  );

///
/// Free all memory associated with a #FstatAtomVector structure.
///
void
XLALDestroyFstatAtomVector(
  FstatAtomVector *atoms                        ///< [in] #FstatAtomVector structure to be freed.
  );

///
/// Create a #MultiFstatAtomVector of the given length.
///
MultiFstatAtomVector*
XLALCreateMultiFstatAtomVector(
  const UINT4 length                            ///< [in] Length of the #MultiFstatAtomVector.
  );

///
/// Free all memory associated with a #MultiFstatAtomVector structure.
///
void
XLALDestroyMultiFstatAtomVector(
  MultiFstatAtomVector *atoms                   ///< [in] #MultiFstatAtomVector structure to be freed.
  );

///
/// Create a \c FstatInput structure which will compute the \f$\mathcal{F}\f$-statistic using demodulation \cite Williams1999.
///
FstatInput*
XLALCreateFstatInput_Demod(

  /// [in] Number of terms to keep in the Dirichlet kernel.
  const UINT4 Dterms,

  /// [in] Which Fstat algorithm/method to use: see the documentation for #FstatMethodType
  const FstatMethodType FstatMethod

  );

///
/// Create a \c FstatInput structure which will compute the \f$\mathcal{F}\f$-statistic using resampling \cite JKS98.
///
FstatInput*
XLALCreateFstatInput_Resamp(
  void
  );

///
/// Setup a \c FstatInput structure for computing the \f$\mathcal{F}\f$-statistic.
///
int
XLALSetupFstatInput(

  /// [in] Input data structure created by one of the setup functions.
  FstatInput *input,

  /// [in] Catalog of SFTs to either load from files, or generate in memory.  The \c locator field
  /// of each ::SFTDescriptor must be \c !NULL for SFT loading, and \c NULL for SFT generation.
  const SFTCatalog *SFTcatalog,

  /// [in] Minimum instantaneous frequency which will be covered over the SFT time span.
  const REAL8 minCoverFreq,

  /// [in] Maximum instantaneous frequency which will be covered over the SFT time span.
  const REAL8 maxCoverFreq,

  /// [in] Optional vector of parameters of CW signals to simulate and inject.
  const PulsarParamsVector *injectSources,

  /// [in] Optional array of single-sided PSD values governing fake Gaussian noise generation.  If
  /// supplied, then fake Gaussian noise with the given PSD values will be added to the SFTs.
  const MultiNoiseFloor *injectSqrtSX,

  /// [in] Optional array of single-sided PSD values governing the calculation of SFT noise weights.
  /// If supplied, then SFT noise weights are calculated from constant spectra with the given PSD
  /// values; otherwise, SFT noise weights are calculated from PSDs computed from a running median
  /// of the SFTs themselves.
  const MultiNoiseFloor *assumeSqrtSX,

  /// [in] If SFT noise weights are calculated from the SFTs, the running median window length to use.
  const UINT4 runningMedianWindow,

  /// [in] Ephemerides for the time-span of the SFTs.
  const EphemerisData *ephemerides,

  /// [in] Barycentric transformation precision.
  const SSBprecision SSBprec,

  /// [in] Seed value used for random number generation, if required.
  const UINT4 randSeed

  );

///
/// Returns the detector information stored in a \c FstatInput structure.
///
const MultiLALDetector*
XLALGetFstatInputDetectors(
  const FstatInput* input                   ///< [in] \c FstatInput structure.
  );

///
/// Returns the SFT timestamps stored in a \c FstatInput structure.
///
const MultiLIGOTimeGPSVector*
XLALGetFstatInputTimestamps(
  const FstatInput* input                   ///< [in] \c FstatInput structure.
  );

///
/// Returns the multi-detector noise weights stored in a \c FstatInput structure.
///
const MultiNoiseWeights*
XLALGetFstatInputNoiseWeights(
  const FstatInput* input                   ///< [in] \c FstatInput structure.
  );

///
/// Returns the multi-detector state series stored in a \c FstatInput structure.
///
const MultiDetectorStateSeries*
XLALGetFstatInputDetectorStates(
  const FstatInput* input                   ///< [in] \c FstatInput structure.
  );

#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(FstatResults**, Fstats));
#endif

///
/// Compute the \f$\mathcal{F}\f$-statistic over a band of frequencies.
///
int
XLALComputeFstat(

  /// [in/out] Address of a pointer to a #FstatResults results structure.  If the pointer is
  /// \c NULL, this function will allocate the structure.
  FstatResults **Fstats,

  /// [in] Input data structure created by one of the setup functions.
  FstatInput *input,

  /// [in] Doppler parameters, including the starting frequency, at which the \f$2\mathcal{F}\f$
  /// are to be computed.
  const PulsarDopplerParams *doppler,

  /// [in] Required spacing in frequency between each \f$\mathcal{F}\f$-statistic.
  const REAL8 dFreq,

  /// [in] Number of frequencies at which the \f$2\mathcal{F}\f$ are to be computed.
  const UINT4 numFreqBins,

  /// [in] Bit-field of which \f$\mathcal{F}\f$-statistic quantities to compute.
  const FstatQuantities whatToCompute

  );

///
/// Free all memory associated with a \c FstatInput structure.
///
void
XLALDestroyFstatInput(
  FstatInput* input                         ///< [in] \c FstatInput structure to be freed.
  );

///
/// Free all memory associated with a #FstatResults structure.
///
void
XLALDestroyFstatResults(
  FstatResults* Fstats                          ///< [in] #FstatResults structure to be freed.
  );

///
/// Add +4 to any multi-detector or per-detector 2F values computed by XLALComputeFstat().
/// This is for compatibility with programs which expect this normalisation if SFTs do not
/// contain noise, e.g. \c lalapps_ComputeFStatistic with the \c --SignalOnly option.
///
int
XLALAdd4ToFstatResults(
  FstatResults* Fstats                          ///< [in] #FstatResults structure.
  );

///
/// Estimate the amplitude parameters of a pulsar CW signal, given its phase parameters,
/// constituent parts of the \f$\mathcal{F}\f$-statistic, and antenna pattern matrix.
///
/// \note Parameter-estimation based on large parts on Yousuke's notes and implemention (in CFSv1),
/// extended for error-estimation.
///
int
XLALEstimatePulsarAmplitudeParams(
  PulsarCandidate *pulsarParams,                ///< [in,out] Pulsar candidate parameters.
  const LIGOTimeGPS* FaFb_refTime,              ///< [in] Reference time of \f$F_a\f$ and \f$F_b\f$, may differ from pulsar candidate reference time.
  const COMPLEX16 Fa,                           ///< [in] Complex \f$\mathcal{F}\f$-statistic amplitude \f$F_a\f$.
  const COMPLEX16 Fb,                           ///< [in] Complex \f$\mathcal{F}\f$-statistic amplitude \f$F_b\f$.
  const AntennaPatternMatrix *Mmunu             ///< [in] Antenna pattern matrix \f$M_{\mu\nu}\f$.
  );

///
/// Convert amplitude params from 'physical' coordinates \f$(h_0, \cos\iota, \psi, \phi_0)\f$ into
/// 'canonical' coordinates \f$A^\mu = (A_1, A_2, A_3, A_4)\f$. The equations can be found in
/// \cite JKS98 or \cite Prix07 Eq.(2).
///
int
XLALAmplitudeParams2Vect(
  PulsarAmplitudeVect A_Mu,                     ///< [out] Canonical amplitude coordinates \f$A^\mu = (A_1, A_2, A_3, A_4)\f$.
  const PulsarAmplitudeParams Amp               ///< [in] Physical amplitude params \f$(h_0, \cos\iota, \psi, \phi_0)\f$.
  );

///
/// Compute amplitude params \f$(h_0, \cos\iota, \psi, \phi_0)\f$ from amplitude-vector \f$A^\mu = (A_1, A_2, A_3, A_4)\f$.
/// Adapted from algorithm in XLALEstimatePulsarAmplitudeParams().
///
int
XLALAmplitudeVect2Params(
  PulsarAmplitudeParams *Amp,                   ///< [out] Physical amplitude params \f$(h_0, \cos\iota, \psi, \phi_0)\f$.
  const PulsarAmplitudeVect A_Mu                ///< [in] Canonical amplitude coordinates \f$A^\mu = (A_1, A_2, A_3, A_4)\f$.
  );

REAL8
XLALComputeFstatFromAtoms ( const MultiFstatAtomVector *multiFstatAtoms,
			    const INT4 X );

// @}

#ifdef  __cplusplus
}
#endif

#endif // _COMPUTEFSTAT_H
