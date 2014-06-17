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

// @{

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
  FSTATQ_NONE           = 0x00,
  FSTATQ_2F             = 0x01,         ///< Compute multi-detector \f$2\mathcal{F}\f$.
  FSTATQ_FAFB           = 0x02,         ///< Compute multi-detector \f$F_a\f$ and \f$F_b\f$.
  FSTATQ_2F_PER_DET     = 0x04,         ///< Compute \f$2\mathcal{F}\f$ for each detector.
  FSTATQ_FAFB_PER_DET   = 0x08,         ///< Compute \f$F_a\f$ and \f$F_b\f$ for each detector.
  FSTATQ_ATOMS_PER_DET  = 0x10,         ///< Compute per-SFT \f$\mathcal{F}\f$-statistic atoms for each detector (\a Demod only).
  FSTATQ_LAST           = 0x20
} FstatQuantities;

///
/// Different methods available to compute the F-statistic, falling into two broad classes:
/// * \a Demod: Dirichlet kernel-based demodulation \cite Williams1999
/// * \a Resamp: FFT-based resampling \cite JKS98
///
typedef enum tagFstatMethodType {

  /// \cond DONT_DOXYGEN
  // Overall/demod start marker; set to 1 to allow range-check without warnings
  FMETHOD_START = 1, FMETHOD_DEMOD_START = FMETHOD_START,
  /// \endcond

  FMETHOD_DEMOD_GENERIC,	///< \a Demod: generic C hotloop, works for any number of Dirichlet kernel terms \f$\text{Dterms}\f$
  FMETHOD_DEMOD_OPTC,		///< \a Demod: gptimized C hotloop using Akos' algorithm, only works for \f$\text{Dterms} \lesssim 20\f$
  FMETHOD_DEMOD_SSE,		///< \a Demod: SSE hotloop with precalc divisors, uses fixed \f$\text{Dterms} = 8\f$
  FMETHOD_DEMOD_ALTIVEC,	///< \a Demod: Altivec hotloop variant, uses fixed \f$\text{Dterms} = 8\f$

  /// \cond DONT_DOXYGEN
  // Demod end marker, resamp start marker
  FMETHOD_DEMOD_END, FMETHOD_RESAMP_START = FMETHOD_DEMOD_END,
  /// \endcond

  FMETHOD_RESAMP_GENERIC,	///< \a Resamp: generic implementation

  /// \cond DONT_DOXYGEN
  // Resamp/overall end marker
  FMETHOD_RESAMP_END, FMETHOD_END = FMETHOD_RESAMP_END
  /// \endcond

} FstatMethodType;

///
/// Determine if #FstatMethodType value \c x belongs to the \a Demod class
///
#define XLALFstatMethodClassIsDemod(x)  ( ((x) > FMETHOD_DEMOD_START )  && ((x) < FMETHOD_DEMOD_END ) )

///
/// Determine if #FstatMethodType value \c x belongs to the \a Resamp class
///
#define XLALFstatMethodClassIsResamp(x) ( ((x) > FMETHOD_RESAMP_START ) && ((x) < FMETHOD_RESAMP_END ) )

///
/// Provide a 'best guess' of the fastest available \a Demod method variant. Useful as a user default value.
///
extern const int FMETHOD_DEMOD_BEST;

///
/// Provide a 'best guess' of the fastest available \a Resamp method variant. Useful as a user default value.
///
extern const int FMETHOD_RESAMP_BEST;

///
/// Struct for collecting 'lower-level' tuning and method-specific parameters to be passed to the
/// \f$\mathcal{F}\f$-statistic setup function.
///
typedef struct tagFstatExtraParams {
  UINT4 randSeed;                 ///< Random-number seed value used for fake Gaussian noise generation.
  SSBprecision SSBprec;           ///< Barycentric transformation precision.
  UINT4 Dterms;                   ///< Number of Dirichlet kernel terms, used by some \a Demod methods; see #FstatMethodType
} FstatExtraParams;

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
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D_PTR_1D(FstatResults, COMPLEX8, FaPerDet, UINT4, numDetectors, numFreqBins));
  SWIGLAL(ARRAY_1D_PTR_1D(FstatResults, COMPLEX8, FaPerDet, UINT4, numDetectors, numFreqBins));
#endif // SWIG
  COMPLEX8 *FaPerDet[PULSAR_MAX_DETECTORS];
  COMPLEX8 *FbPerDet[PULSAR_MAX_DETECTORS];

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

// ---------- API function prototypes ----------
const CHAR *XLALGetFstatMethodName ( FstatMethodType i);
int XLALFstatMethodIsAvailable ( FstatMethodType i );
const CHAR *XLALFstatMethodHelpString ( void );
int XLALParseFstatMethodString ( FstatMethodType *Fmethod, const char *s );

FstatInputVector* XLALCreateFstatInputVector ( const UINT4 length );
void XLALDestroyFstatInputVector ( FstatInputVector* input );
FstatAtomVector* XLALCreateFstatAtomVector ( const UINT4 length );
void XLALDestroyFstatAtomVector ( FstatAtomVector *atoms );
MultiFstatAtomVector* XLALCreateMultiFstatAtomVector ( const UINT4 length );
void XLALDestroyMultiFstatAtomVector ( MultiFstatAtomVector *atoms );

FstatInput *
XLALCreateFstatInput ( const SFTCatalog *SFTcatalog, const REAL8 minCoverFreq, const REAL8 maxCoverFreq,
                       const PulsarParamsVector *injectSources, const MultiNoiseFloor *injectSqrtSX,
                       const MultiNoiseFloor *assumeSqrtSX, const UINT4 runningMedianWindow,
                       const EphemerisData *ephemerides,
                       const FstatMethodType FstatMethod, const FstatExtraParams *extraParams );

const MultiLALDetector* XLALGetFstatInputDetectors ( const FstatInput* input );
const MultiLIGOTimeGPSVector* XLALGetFstatInputTimestamps ( const FstatInput* input );
const MultiNoiseWeights* XLALGetFstatInputNoiseWeights ( const FstatInput* input );
const MultiDetectorStateSeries* XLALGetFstatInputDetectorStates ( const FstatInput* input );

#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(FstatResults**, Fstats));
#endif
int XLALComputeFstat ( FstatResults **Fstats, FstatInput *input, const PulsarDopplerParams *doppler,
                       const REAL8 dFreq, const UINT4 numFreqBins, const FstatQuantities whatToCompute );

void XLALDestroyFstatInput ( FstatInput* input );
void XLALDestroyFstatResults ( FstatResults* Fstats );
int XLALAdd4ToFstatResults ( FstatResults* Fstats );

int XLALEstimatePulsarAmplitudeParams ( PulsarCandidate *pulsarParams, const LIGOTimeGPS* FaFb_refTime,
                                        const COMPLEX16 Fa, const COMPLEX16 Fb, const AntennaPatternMatrix *Mmunu );

int XLALAmplitudeParams2Vect ( PulsarAmplitudeVect A_Mu, const PulsarAmplitudeParams Amp );
int XLALAmplitudeVect2Params( PulsarAmplitudeParams *Amp, const PulsarAmplitudeVect A_Mu );

REAL8 XLALComputeFstatFromAtoms ( const MultiFstatAtomVector *multiFstatAtoms, const INT4 X );
REAL8 XLALComputeFstatFromFaFb ( COMPLEX16 Fa, COMPLEX16 Fb, REAL8 A, REAL8 B, REAL8 C, REAL8 E, REAL8 Dinv );

// @}

#ifdef  __cplusplus
}
#endif

#endif // _COMPUTEFSTAT_H
