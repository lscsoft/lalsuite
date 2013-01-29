/*
 * Copyright (C) 2004, 2005 Reinhard Prix
 * Copyright (C) 2004, 2005 Greg Mendell
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/* NOTES: */
/* 07/14/04 gam; add functions LALFastGeneratePulsarSFTs and LALComputeSkyAndZeroPsiAMResponse */
/* 10/08/04 gam; fix indexing into trig lookup tables (LUTs) by having table go from -2*pi to 2*pi */
/* 09/07/05 gam; Add Dterms parameter to LALFastGeneratePulsarSFTs; use this to fill in SFT bins with fake data as per LALDemod else fill in bin with zero */

#ifndef _GENERATEPULSARSIGNAL_H  /* Double-include protection. */
#define _GENERATEPULSARSIGNAL_H

#include <lal/LALDatatypes.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/GenerateSpinOrbitCW.h>
#include <lal/Date.h>
#include <lal/LALBarycenter.h>
#include <lal/PulsarDataTypes.h>
#include <lal/ComputeSky.h>
#include <lal/ComputeSkyBinary.h>
#include <lal/Window.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup GeneratePulsarSignal_h
 * \author Reinhard Prix, Greg Mendell
 * \date 2005
 *
 * \brief Pulsar signal-generation routines for hardware- and software-injections.
 *
 * \heading{Description}
 *
 * - The main function LALGeneratePulsarSignal() generates a fake
 * pulsar-signal, either for an isolated or a binary pulsar. It returns a
 * time-series with the generated signal as received by the detector.
 *
 * - The time-series can then be turned into a vector of short time-base FFTs
 * (the so-called "SFTs") by using the function LALSignalToSFTs().
 * These SFTs are the data-format used by most frequency-domain pulsar codes,
 * therefore these functions can be directly used in a Monte-Carlo
 * injection driver.
 *
 * This module also contains a few more general-purpose helper-functions:
 *
 *
 * - Namely, LALConvertSSB2GPS() and LALConvertGPS2SSB()
 * which convert arrival times for a given source (not necessarily a
 * pulsar!) the detector ("GPS") and the solar-system barycenter ("SSB").
 * NOTE: only the source-location (<tt>params-\>pulsar.position</tt>), the
 * detector-site (<tt>params-\>site</tt>) and the ephemeris-data
 * (<tt>params-\>ephemerides</tt>)are used from the
 * PulsarSignalParams structure.
 *
 *
 * \heading{Algorithm}
 *
 * LALGeneratePulsarSignal() is basically a wrapper for the two
 * LAL-functions GenerateSpinOrbitCW() to produce the source-signal,
 * and LALPulsarSimulateCoherentGW() to turn it into a time-series at the detector.
 *
 * LALSignalToSFTs() uses LALForwardRealFFT() appropriately on the input-timeseries to
 * produce the required output-SFTs ( v2-normalization! ).
 *
 *
 * \note
 *
 * LALSignalToSFTs() currently enforces the constraint of
 * <tt>Tsft * Band</tt> being an integer, so that the number of
 * time-samples per SFT is even. This follows <tt>makefakedata_v2</tt>.
 *
 * Furthermore, if the timestamps for SFT-creation passed to
 * LALSignalToSFTs() do not fit exactly on a time-step of the
 * input time-series, it will be "nudged" to the closest one.
 * If <tt>lalDebugLevel>0</tt> a warning will be printed about this.
 * The user could also see this effect in the actual timestamps of the
 * SFTs returned.
 *
 *
 * The FFTW-"plan" is currently created using the \c ESTIMATE flag,
 * which is fast but only yields an approximate plan. Better results
 * might be achieved by using \c MEASURE and an appropriate buffering
 * of the resulting plan (which doesnt change provided the SFT-length is
 * the same). Measuring the plan takes longer but might lead to
 * substantial speedups for many FFTs, which seems the most likely situation.
 *
 *
 *
 * \heading{Use of LALFastGeneratePulsarSFTs()}
 *
The functions LALComputeSkyAndZeroPsiAMResponse() and LALFastGeneratePulsarSFTs()
use approximate analytic formulas to generate SFTs.  This should be significantly
faster than LALGeneratePulsarSignal() and LALSignalToSFTs(), which generate the
time series data and then FFT it.  Simple tests performed by the code in
GeneratePulsarSignalTest.c indicate that the maximum modulus of the SFTs output
by the approximate and exact codes differs by less that 10\%.  Since the tests
are not exhaustive, the user should use caution and conduct their own test
to compare LALComputeSkyAndZeroPsiAMResponse() and LALFastGeneratePulsarSFTs() with
LALGeneratePulsarSignal() and LALSignalToSFTs().

The strain of a periodic signal at the detector is given by
\f[
h(t) = F_+(t) A_+ {\rm cos}\Phi(t) + F_\times(t) A_\times {\rm sin}\Phi(t),
\f]
where \f$F_+\f$ and \f$F_\times\f$ are the usual beam pattern response functions,
\f$A_+\f$ and \f$A_\times\f$ are the amplitudes of the gravitational wave for the
plus and cross polarizations, and \f$\Phi\f$ is the phase.  The phase contains modulations
from doppler shifts due to the relative motion between the source and the
detector and the spin evolution of the source.  (The functions discussed here
support both isolated sources and those in binary systems. The binary case
has not been tested.)

If we Taylor expand the phase out to first order about the time at the midpoint of
an SFT and approximate \f$F_+\f$ and \f$F_\times\f$ as constants, for one SFT we can write
\f[
\Phi(t) \approx \Phi_{1/2} + 2\pi f_{1/2}(t - t_{1/2}).
\f]
The strain at discrete time \f$t_j\f$, measured from the start of the SFT, can
thus be approximated as
\f[
h_j \approx F_{+ 1/2} A_+ {\rm cos} [\Phi_{1/2} + 2\pi f_{1/2}(t_0 + t_j - t_{1/2})]
+ F_{\times 1/2} A_\times {\rm sin} [\Phi_{1/2} + 2\pi f_{1/2}(t_0 + t_j - t_{1/2})],
\f]
where \f$t_0\f$ is the time as the start of the SFT, and \f$t_{1/2} - t_0 = T_{\rm sft}/2\f$,
where \f$T_{\rm sft}\f$ is the duration of one SFT.  This simplifies to
\f[
h_j \approx F_{+ 1/2} A_+ {\rm cos} (\Phi_0 + 2\pi f_{1/2}t_j)
+ F_{\times 1/2} A_\times {\rm sin} (\Phi_0 + 2\pi f_{1/2}t_j),
\f]
where \f$\Phi_0\f$ is the phase at the start of the SFT
(not the initial phase at the start of the observation), i.e.,
\f[
\Phi_0 = \Phi_{1/2} - 2 \pi f_{1/2} (T_{\rm sft} / 2).
\f]
Note that these are the same approximations used by LALDemod().

One can show that the Discrete Fourier Transform (DFT) of \f$h_j\f$ above is:
\f[
\tilde{h}_k = e^{i\Phi_0}  \frac{1}{2} ( F_{+ 1/2} A_+ - i F_{\times 1/2} A_\times)
\frac{ 1 - e^{2\pi i (\kappa - k)}}{1 - e^{2\pi i (\kappa - k)/N} }
\\
+ e^{-i\Phi_0}  \frac{1}{2} ( F_{+ 1/2} A_+ + i F_{\times 1/2} A_\times)
\frac{ 1 - e^{-2\pi i (\kappa + k)}}{ 1 - e^{-2\pi i (\kappa + k)/N} }
\f]
where \f$N\f$ is the number of time samples used to find the
DFT (i.e., the sample rate times \f$T_{\rm sft}\f$), and
\f[
\kappa \equiv f_{1/2} T_{\rm sft},
\f]
is usually not an integer.

Note that the factor \f$e^{\pm 2\pi i k}\f$ in the numerators of the equation for \f$\tilde{h}_k\f$
equals 1.  Furthermore, for \f$0 < \kappa < N/2\f$ and \f$|\kappa - k| << N\f$ the first term
dominates and can be Taylor expanded to give:
\f[
\tilde{h}_k = N e^{i\Phi_0} \frac{1}{2} ( F_{+ 1/2} A_+ - i F_{\times 1/2} A_\times)
\left [ \, \frac{ {\rm sin} (2\pi\kappa)}{2 \pi (\kappa - k) } \,
+ \, i \frac{ 1 - {\rm cos} (2\pi\kappa)}{2 \pi (\kappa - k) } \, \right ]
\f]
Note that the last factor in square brackets is \f$P_{\alpha k}^*\f$ and
\f$e^{i\Phi_0} = Q_{\alpha}^*\f$ used by LALDemod.

\heading{Example pseudocode}

The structs used by LALComputeSkyAndZeroPsiAMResponse() and LALFastGeneratePulsarSFTs()
are given in previous sections, and make use of those used by
LALGeneratePulsarSignal() and LALSignalToSFTs() plus a small number of
additional parameters.  Thus it is fairly easy to change between the above
approximate routines the exact routines. See GeneratePulsarSignalTest.c for
an example implementation of the code.

Note that one needs to call LALComputeSkyAndZeroPsiAMResponse() once per sky position,
and then call LALFastGeneratePulsarSFTs() for each set of signal parameters at that
sky position.  Thus, one could perform a Monte Carlo simulation, as shown
by the pseudo code:

\code
loop over sky positions {
   ...
   LALComputeSkyAndZeroPsiAMResponse();
   ...
   loop over spindown {
      ...
      loop over frequencies {
         ...
         LALFastGeneratePulsarSFTs();
         ...
      }
      ...
   }
   ...
}

\endcode

\heading{Notes on LALFastGeneratePulsarSFTs()}

-#  If \c *outputSFTs sent to LALFastGeneratePulsarSFTs() is \c NULL then
LALFastGeneratePulsarSFTs() allocates memory for the output SFTs; otherwise it assumes
memory has already been allocated.  Thus, the user does not have to deallocate
memory for the SFTs until all calls to LALFastGeneratePulsarSFTs() are completed.

-# \c fHeterodyne and <tt>0.5 * samplingRate</tt> set in the PulsarSignalParams struct
give the start frequency and frequency band of the SFTs output from LALFastGeneratePulsarSFTs().

-# If \c resTrig is set to zero in the SFTandSignalParams struct, then
the C math libary \c cos() \c sin() functions are called, else lookup tables (LUTs) are used
for calls to trig functions.  There may be a slight speedup in using LUTs.

-# To maximize the speed of SFT generations, LALFastGeneratePulsarSFTs() only generates
values for the bins in the band <tt>2*Dterms</tt> centered on the signal frequency in each SFT. Dterms must be
greater than zero and less than or equal to the number of frequency bins in the output SFTs. Note that
Dterms is used the same way here as it is in LALDemod(). Nothing is done to the other bins, unless
\c *outputSFTs is \c NULL; then, since memory is allocates for the output SFTs, the bins
not in the <tt>2*Dterms</tt> band are initialized to zero.

*/
/*@{*/

/** \name Error codes */
/*@{*/
#define GENERATEPULSARSIGNALH_ENULL 		1	/**< Arguments contained an unexpected null pointer */
#define GENERATEPULSARSIGNALH_ENONULL		2	/**< Output pointer is not NULL */
#define GENERATEPULSARSIGNALH_EMEM		3	/**< Out of memory */
#define GENERATEPULSARSIGNALH_ESAMPLING		4	/**< Waveform sampling interval too large. */
#define GENERATEPULSARSIGNALH_ESSBCONVERT	5	/**< SSB-\>GPS iterative conversion failed */
#define GENERATEPULSARSIGNALH_ESYS		6	/**< System error, probably while File I/O */
#define GENERATEPULSARSIGNALH_ETIMEBOUND	7	/**< Timestamp outside of allowed time-interval */
#define GENERATEPULSARSIGNALH_ENUMSFTS		8	/**< Inconsistent number of SFTs in timestamps and noise-SFTs */
#define GENERATEPULSARSIGNALH_EINCONSBAND	9	/**< Inconsistent values of sampling-rate and Tsft */
#define GENERATEPULSARSIGNALH_ENOISEDELTAF	10	/**< Frequency resolution of noise-SFTs inconsistent with signal */
#define GENERATEPULSARSIGNALH_ENOISEBAND	11	/**< Frequency band of noise-SFTs inconsistent with signal */
#define GENERATEPULSARSIGNALH_ENOISEBINS	12	/**< Frequency bins of noise-SFTs inconsistent with signal */
#define GENERATEPULSARSIGNALH_EBADCOORDS	13	/**< Current code requires sky position in equatorial coordinates */
#define GENERATEPULSARSIGNALH_ELUTS		14	/**< Lookup tables (LUTs) for trig functions must be defined on domain -2pi to 2pi inclusive */
#define GENERATEPULSARSIGNALH_EDTERMS		15	/**< Dterms must be greater than zero and less than or equal to half the number of SFT bins */
#define GENERATEPULSARSIGNALH_EINPUT		16	/**< Invalid input-arguments to function */
#define GENERATEPULSARSIGNALH_EDETECTOR		17	/**< Unknown detector-name */
/*@}*/

/** \cond DONT_DOXYGEN */
#define GENERATEPULSARSIGNALH_MSGENULL 		"Arguments contained an unexpected null pointer"
#define GENERATEPULSARSIGNALH_MSGENONULL	"Output pointer is not NULL"
#define GENERATEPULSARSIGNALH_MSGEMEM		"Out of memory"
#define GENERATEPULSARSIGNALH_MSGESAMPLING	"Waveform sampling interval too large."
#define GENERATEPULSARSIGNALH_MSGESSBCONVERT	"SSB->GPS iterative conversion failed"
#define GENERATEPULSARSIGNALH_MSGESYS		"System error, probably while File I/O"
#define GENERATEPULSARSIGNALH_MSGETIMEBOUND	"Timestamp outside of allowed time-interval"
#define GENERATEPULSARSIGNALH_MSGENUMSFTS	"Inconsistent number of SFTs in timestamps and noise-SFTs"
#define GENERATEPULSARSIGNALH_MSGEINCONSBAND	"Inconsistent values of sampling-rate and Tsft"
#define GENERATEPULSARSIGNALH_MSGENOISEDELTAF	"Frequency resolution of noise-SFTs inconsistent with signal"
#define GENERATEPULSARSIGNALH_MSGENOISEBAND	"Frequency band of noise-SFTs inconsistent with signal"
#define GENERATEPULSARSIGNALH_MSGENOISEBINS	"Frequency bins of noise-SFTs inconsistent with signal"
#define GENERATEPULSARSIGNALH_MSGEBADCOORDS	"Current code requires sky position in equatorial coordinates"
#define GENERATEPULSARSIGNALH_MSGELUTS		"Lookup tables (LUTs) for trig functions must be defined on domain -2pi to 2pi inclusive"
#define GENERATEPULSARSIGNALH_MSGEDTERMS	"Dterms must be greater than zero and less than or equal to half the number of SFT bins"
#define GENERATEPULSARSIGNALH_MSGEINPUT		"Invalid input-arguments to function"
#define GENERATEPULSARSIGNALH_MSGEDETECTOR	"Unknown detector-name"
/** \endcond */

/** Input parameters to GeneratePulsarSignal(), defining the source and the time-series
 */
typedef struct tagPulsarSignalParams {
  /* source-parameters */
  PulsarSourceParams pulsar;	/**< the actual pulsar-source */
  BinaryOrbitParams *orbit;	/**< and its binary orbit (NULL if isolated pulsar) */

  /* characterize the detector */
  COMPLEX8FrequencySeries *transfer; /**< detector transfer function (NULL if not used) */
  LALDetector *site;		/**< detector location and orientation */
  EphemerisData *ephemerides;	/**< Earth and Sun ephemerides */

  /* characterize the output time-series */
  LIGOTimeGPS startTimeGPS;     /**< start time of output time series */
  UINT4 duration;           	/**< length of time series in seconds */
  REAL8 samplingRate;		/**< sampling rate of time-series (= 2 * frequency-Band) */
  REAL8 fHeterodyne;		/**< heterodyning frequency for output time-series */
  UINT4 dtDelayBy2; 		/**< half-interval for the Doppler delay look-up table for LALPulsarSimulateCoherentGW() */
  UINT4 dtPolBy2; 		/**< half-interval for the polarisation response look-up table for LALPulsarSimulateCoherentGW() */
} PulsarSignalParams;

/** Parameters defining the SFTs to be returned from LALSignalToSFTs().
 */
typedef struct tagSFTParams {
  REAL8 Tsft;			 /**< length of each SFT in seconds */
  LIGOTimeGPSVector *timestamps; /**< timestamps to produce SFTs for (can be NULL) */
  SFTVector *noiseSFTs;		 /**< noise SFTs to be added (can be NULL) */
  REAL4Window *window;		 /**< window function for the time series (can be NULL, which is equivalent to a rectangular window) */
} SFTParams;

/** Parameters defining the pulsar signal and SFTs used by LALFastGeneratePulsarSFTs().  Lookup tables (LUTs) are
 * used for trig functions if \code resTrig > 0 \endcode the user must then initialize \c trigArg, \c sinVal, and
 * \c cosVal on the domain \f$[-2\pi, 2\pi]\f$ inclusive.  See GeneratePulsarSignalTest.c for an example.
 */
typedef struct tagSFTandSignalParams {
   PulsarSignalParams *pSigParams;
   SFTParams *pSFTParams;
   INT4  nSamples;  /**< nsample from noise SFT header; 2x this equals effective number of time samples  */
   INT4  Dterms;    /**< use this to fill in SFT bins with fake data as per LALDemod else fill in bin with zero */
   INT4  resTrig;   /**< length sinVal, cosVal; domain: -2pi to 2pi; resolution = 4pi/resTrig */
   REAL8 *trigArg;  /**< array of arguments to hold lookup table (LUT) values for doing trig calls */
   REAL8 *sinVal;   /**< sinVal holds lookup table (LUT) values for doing trig sin calls */
   REAL8 *cosVal;   /**< cosVal holds lookup table (LUT) values for doing trig cos calls */
} SFTandSignalParams;


/** Sky Constants and beam pattern response functions used by LALFastGeneratePulsarSFTs().
 * These are output from LALComputeSkyAndZeroPsiAMResponse().
 */
typedef struct tagSkyConstAndZeroPsiAMResponse {
      REAL8  *skyConst;      /**< vector of A and B sky constants */
      REAL4  *fPlusZeroPsi;  /**< vector of Fplus values for psi = 0 at midpoint of each SFT */
      REAL4  *fCrossZeroPsi; /**< vector of Fcross values for psi = 0 at midpoint of each SFT */
} SkyConstAndZeroPsiAMResponse;

/*---------- Global variables ----------*/
/** \name Empty init-structs for the types defined in here */
/*@{*/
extern const PulsarSignalParams empty_PulsarSignalParams;
extern const SFTParams empty_SFTParams;
extern const SFTandSignalParams empty_SFTandSignalParams;
/*@}*/

/* ---------- Function prototypes ---------- */
REAL4TimeSeries *XLALGeneratePulsarSignal ( const PulsarSignalParams *params );
REAL4TimeSeries *XLALGenerateLineFeature ( const PulsarSignalParams *params );
SFTVector *XLALSignalToSFTs ( const REAL4TimeSeries *signalvec, const SFTParams *params );
int XLALConvertGPS2SSB ( LIGOTimeGPS *SSBout, LIGOTimeGPS GPSin, const PulsarSignalParams *params );
int XLALConvertSSB2GPS ( LIGOTimeGPS *GPSout, LIGOTimeGPS GPSin, const PulsarSignalParams *params );

// ----- obsolete and deprecated LAL interface
void LALGeneratePulsarSignal (LALStatus *, REAL4TimeSeries **signalvec, const PulsarSignalParams *params);
void LALSignalToSFTs (LALStatus *, SFTVector **outputSFTs, const REAL4TimeSeries *signalvec, const SFTParams *params);

void LALComputeSkyAndZeroPsiAMResponse (LALStatus *, SkyConstAndZeroPsiAMResponse *output, const SFTandSignalParams *params);
void LALFastGeneratePulsarSFTs (LALStatus *, SFTVector **outputSFTs, const SkyConstAndZeroPsiAMResponse *input, const SFTandSignalParams *params);

/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
