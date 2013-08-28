//
// Copyright (C) 2012, 2013 David Keitel, Bernd Machenschalk, Reinhard Prix, Karl Wette
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
#include <lal/ComputeFstat_OldDemodAPI.h>
#include <lal/ComputeFstat_OldResampAPI.h>

#ifdef  __cplusplus
extern "C" {
#endif

///
/// \defgroup ComputeFstat_h Header ComputeFstat.h
/// \ingroup pkg_pulsarCoh
/// \authors David Keitel, Bernd Machenschalk, Reinhard Prix, Karl Wette
///
/// \brief Unified API for computing the \f$\mathcal{F}\f$-statistic.
///
/// This module provides a unified API for computing the \f$\mathcal{F}\f$-statistic \cite JKS98
/// using different algorithms, e.g. demodulation \cite Williams1999 or resampling \cite JKS98.
/// Each algorithm provides a setup function, <tt>XLALSetupFstat_...()</tt>, which performs all
/// initialisation tasks required for the chosen algorithm.  A setup function must accept the
/// following input arguments:
///
/// - <tt>[in/out]</tt> #MultiSFTVector **<b>multiSFTs</b>: the address of a multi-detector SFT
///   array of type #MultiSFTVector*.  The setup function takes ownership of the SFTs (since,
///   depending on the algorithm, they may be modified and/or destroyed), and sets the supplied
///   pointer to \c NULL to indicate this. The SFTs should therfore not be accessed after the setup
///   function is called.
///
/// - <tt>[in/out]</tt> #MultiNoiseWeights **<b>multiWeights</b>: the address of a multi-detector
///   array containing noise weights.  If the address pointed to is NULL, the
///   \f$\mathcal{F}\f$-statistic is calculated assuming unity noise weights, and is then normalised
///   appropriately.  Otherwise, the setup function takes ownership of the noise weights, and sets
///   the supplied pointer to \c NULL to indicate this. The noise weights should therfore not be
///   accessed after the setup function is called.
///
/// - <tt>[in]</tt> const #EphemerisData *<b>edat</b>: ephemerides for the time-span of the SFTs.
///
/// - <tt>[in]</tt> const #SSBprecision <b>SSBprec</b>: precision of barycentric transformation.
///
/// The setup function returns a pointer to an opaque structure, #FstatInputData, which is shared by
/// all algorithms.  After the initial setup is performed, the calling code passes the
/// #FstatInputData pointer to the function #XLALComputeFstat(), which computes the
/// \f$\mathcal{F}\f$-statistic using the chosen algorithm, and fills a #FstatResults structure with
/// the results.
///
/// \note The \f$\mathcal{F}\f$-statistic algorithm codes are partly descended from earlier
/// implementations found in:
///
/// - LALDemod.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff,
///   Xavier Siemens, Bruce Allen
///
/// - ComputeSky.[ch] by Jolien Creighton, Reinhard Prix, Steve Berukoff
///
/// - LALComputeAM.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff,
///   Xavier Siemens
///
/// - ComputeFStatistic_resamp.c by Pinkesh Patel, Xavier Siemens, Reinhard Prix, Iraj Gholami,
///   Yousuke Itoh, Maria Alessandra Papa
///

// @{

///
/// XLALComputeFstat() input data structure. Encapsulates all data, buffers, etc. used by the
/// \f$\mathcal{F}\f$-statistic algorithms.
///
typedef struct tagFstatInputData FstatInputData;

///
/// A vector of XLALComputeFstat() input data structures, for e.g. computing the
/// \f$\mathcal{F}\f$-statistic for multiple segments.
///
typedef struct tagFstatInputDataVector {
  UINT4 length;				///< Number of elements in array.
  FstatInputData **data;		///< Pointer to the data array.
} FstatInputDataVector;

///
/// Bit-field of \f$\mathcal{F}\f$-statistic quantities which can be computed by XLALComputeFstat().
/// Not all options are supported by all \f$\mathcal{F}\f$-statistic algorithms.
///
typedef enum tagFstatQuantities {
  FSTATQ_NONE		= 0x00,
  FSTATQ_2F		= 0x01,		///< Compute multi-detector \f$2\mathcal{F}\f$.
  FSTATQ_FAFB		= 0x02,		///< Compute multi-detector \f$F_a\f$ and \f$F_b\f$.
  FSTATQ_2F_PER_DET	= 0x04,		///< Compute \f$2\mathcal{F}\f$ for each detector.
  FSTATQ_FAFB_PER_DET	= 0x08,		///< Compute \f$F_a\f$ and \f$F_b\f$ for each detector.
  FSTATQ_ATOMS_PER_DET	= 0x10,		///< Compute per-SFT \f$\mathcal{F}\f$-statistic atoms for each detector (demodulation only).
  FSTATQ_LAST		= 0x20
} FstatQuantities;

///
/// Amplitude modulation coefficient type to use when computing the \f$\mathcal{F}\f$-statistic.
/// Not all options are supported by all \f$\mathcal{F}\f$-statistic algorithms.
///
typedef enum tagDemodAMType {
  DEMODAM_LONG_WAVELENGTH,		///< Long-wavelength limit approximation.
  DEMODAM_RIGID_ADIABATIC,		///< Frequency- and sky-position-dependent rigid adiabatic response tensor (demodulation only).
  DEMODAM_BUFFERED_RIGID_ADIABATIC,	///< Approximated rigid adiabatic by assuming constant response over (small) frequency band.
  DEMODAM_LAST
} DemodAMType;

///
/// Complex \f$\mathcal{F}\f$-statistic amplitudes \f$F_a\f$ and \f$F_b\f$.
///
typedef struct tagFstatFaFb {
  COMPLEX16 Fa;				///< Complex amplitude \f$F_a\f$.
  COMPLEX16 Fb;				///< Complex amplitude \f$F_b\f$.
} FstatFaFb;

///
/// XLALComputeFstat() computed results structure.
///
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
  CmplxAntennaPatternMatrix Mmunu;

  /// Bit-field of which \f$\mathcal{F}\f$-statistic quantities were computed.
  FstatQuantities whatWasComputed;

  /// If #whatWasComputed & FSTATQ_2F is true, the multi-detector \f$2\mathcal{F}\f$ values computed
  /// at #numFreqBins frequencies spaced #dFreq apart.  This array should not be accessed if
  /// #whatWasComputed & FSTATQ_2F is false.
  REAL4 *twoF;

  /// If #whatWasComputed & FSTATQ_PARTS is true, the multi-detector \f$F_a\f$ and \f$F_b\f$
  /// computed at #numFreqBins frequencies spaced #dFreq apart.  This array should not be accessed
  /// if #whatWasComputed & FSTATQ_PARTS is false.
  FstatFaFb *FaFb;

  /// If #whatWasComputed & FSTATQ_2F_PER_DET is true, the \f$2\mathcal{F}\f$ values computed at
  /// #numFreqBins frequencies spaced #dFreq apart, and for #numDetectors detectors.  Only the first
  /// #numDetectors entries will be valid.  This array should not be accessed if #whatWasComputed &
  /// FSTATQ_2F_PER_DET is false.
  REAL4 *twoFPerDet[PULSAR_MAX_DETECTORS];

  /// If #whatWasComputed & FSTATQ_PARTS_PER_DET is true, the \f$F_a\f$ and \f$F_b\f$ values
  /// computed at #numFreqBins frequencies spaced #dFreq apart, and for #numDetectors detectors.
  /// This array should not be accessed if #whatWasComputed & FSTATQ_PARTS_PER_DET is false.
  FstatFaFb *FaFbPerDet[PULSAR_MAX_DETECTORS];

  /// If #whatWasComputed & FSTATQ_ATOMS_PER_DET is true, the per-SFT \f$\mathcal{F}\f$-statistic
  /// multi-atoms computed at #numFreqBins frequencies spaced #dFreq apart.  This array should not
  /// be accessed if #whatWasComputed & FSTATQ_ATOMS_PER_DET is false.
  MultiFstatAtomVector** multiFatoms;

  /// \cond DONT_DOXYGEN
  UINT4 internalalloclen;
  /// \endcond

} FstatResults;

///
/// Create a #FstatInputDataVector of the given length.
///
FstatInputDataVector*
XLALCreateFstatInputDataVector(
  const UINT4 length				///< [in] Length of the #FstatInputDataVector.
  );

///
/// Free all memory associated with a #FstatInputDataVector structure.
///
void
XLALDestroyFstatInputDataVector(
  FstatInputDataVector* input			///< [in] #FstatInputDataVector structure to be freed.
  );

///
/// Setup function for computing the \f$\mathcal{F}\f$-statistic using demodulation.  See description
/// in \ref ComputeFstat_h for further information on XLALComputeFstat() setup functions.
///
FstatInputData*
XLALSetupFstat_Demod(
  MultiSFTVector **multiSFTs,			///< [in/out] Address of multi-detector SFT array.
  MultiNoiseWeights **multiWeights,		///< [in/out] Address of multi-detector noise weights array.
  const EphemerisData *edat,			///< [in] Ephemerides over SFT time-span.
  const SSBprecision SSBprec,			///< [in] Barycentric transformation precision.
  const DemodAMType demodAM,			///< [in] Amplitude modulation coefficient type to use.
  const UINT4 Dterms				///< [in] Number of terms to keep in Dirichlet kernel.
  );

///
/// Setup function for computing the \f$\mathcal{F}\f$-statistic using resampling.  See description
/// in \ref ComputeFstat_h for further information on XLALComputeFstat() setup functions.
///
FstatInputData*
XLALSetupFstat_Resamp(
  MultiSFTVector **multiSFTs,			///< [in/out] Address of multi-detector SFT array.
  MultiNoiseWeights **multiWeights,		///< [in/out] Address of multi-detector noise weights array.
  const EphemerisData *edat,			///< [in] Ephemerides over SFT time-span.
  const SSBprecision SSBprec			///< [in] Barycentric transformation precision.
  );

///
/// Compute the \f$\mathcal{F}\f$-statistic over a band of frequencies.
///
int
XLALComputeFstat(
  FstatResults **Fstats,			///< [in/out] Address of a pointer to a #FstatResults results structure.  If the pointer is NULL, this function will allocate the structure.
  FstatInputData *input,			///< [in] Input data structure created by one of the setup functions.
  const PulsarDopplerParams *doppler,		///< [in] Doppler parameters, including the starting frequency, at which the \f$2\mathcal{F}\f$ are to be computed.
  const REAL8 dFreq,				///< [in] Required spacing in frequency between each \f$\mathcal{F}\f$-statistic.
  const UINT4 numFreqBins,			///< [in] Number of frequencies at which the \f$2\mathcal{F}\f$ are to be computed.
  const FstatQuantities whatToCompute		///< [in] Bit-field of which \f$\mathcal{F}\f$-statistic quantities were computed.
  );

///
/// Free all memory associated with a #FstatInputData structure.
///
void
XLALDestroyFstatInputData(
  FstatInputData* input				///< [in] #FstatInputData structure to be freed.
  );

///
/// Free all memory associated with a #FstatResults structure.
///
void
XLALDestroyFstatResults(
  FstatResults* Fstats				///< [in] #FstatResults structure to be freed.
  );

// @}

#ifdef  __cplusplus
}
#endif

#endif // _COMPUTEFSTAT_H
