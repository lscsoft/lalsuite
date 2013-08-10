/*
 * Copyright (C) 2005 Reinhard Prix
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
#ifndef _COMPUTEFSTAT_OLDDEMODAPI_H  /* Double-include protection. */
#define _COMPUTEFSTAT_OLDDEMODAPI_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup ComputeFstat_OldDemodAPI_h Header ComputeFstat_OldDemodAPI.h
 * \ingroup pkg_pulsarCoh
 * \author Reinhard Prix
 * \date 2005
 * \brief Header-file defining the API for the F-statistic functions.
 *
 * This code is (partly) a descendant of an earlier implementation found in
 * LALDemod.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens, Bruce Allen
 * ComputSky.[ch] by Jolien Creighton, Reinhard Prix, Steve Berukoff
 * LALComputeAM.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens
 *
 */
/*@{*/

/*---------- exported INCLUDES ----------*/
#include <lal/LALComputeAM.h>
#include <lal/ComplexAM.h>
#include <lal/SSBtimes.h>
#include <lal/PulsarDataTypes.h>
#include <lal/DetectorStates.h>
#include <gsl/gsl_vector.h>

/*---------- exported DEFINES ----------*/

/** \name Error codes */
/*@{*/
#define COMPUTEFSTATC_ENULL 		1
#define COMPUTEFSTATC_ENONULL 		2
#define COMPUTEFSTATC_EINPUT   		3
#define COMPUTEFSTATC_EMEM   		4
#define COMPUTEFSTATC_EXLAL		5
#define COMPUTEFSTATC_EIEEE		6

#define COMPUTEFSTATC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATC_MSGENONULL 	"Output pointer is non-NULL"
#define COMPUTEFSTATC_MSGEINPUT   	"Invalid input"
#define COMPUTEFSTATC_MSGEMEM   	"Out of memory. Bad."
#define COMPUTEFSTATC_MSGEXLAL		"XLAL function call failed"
#define COMPUTEFSTATC_MSGEIEEE		"Floating point failure"
/*@}*/

/*---------- exported types ----------*/

#ifndef SWIG
struct tagMultiFstatAtomVector;
#endif

/** Type containing F-statistic proper plus the two complex amplitudes Fa and Fb (for ML-estimators) */
typedef struct tagFcomponents {
  REAL8 F;				/**< F-statistic value */
  REAL8 FX[PULSAR_MAX_DETECTORS];		/**< vector of single-detector F-statistic values (array of fixed size) */
  UINT4 numDetectors;			/**< number of detectors = effective vector length. numDetectors=0 should make all code ignore the FX field. */
  LIGOTimeGPS refTime;			/**< 'internal' refTime used to compute the F-statistic: only relevant for phase of complex amplitudes {Fa,Fb} */
  COMPLEX16 Fa;				/**< complex amplitude Fa */
  COMPLEX16 Fb;				/**< complex amplitude Fb */
  struct tagMultiFstatAtomVector *multiFstatAtoms;/**< per-IFO, per-SFT arrays of F-stat 'atoms', ie quantities required to compute F-stat */
} Fcomponents;

/** [opaque] type holding a ComputeFBuffer for use in the resampling F-stat codes */
typedef struct tagComputeFBuffer_RS ComputeFBuffer_RS;

/** Extra parameters controlling the actual computation of F */
typedef struct tagComputeFParams {
  UINT4 Dterms;		/**< how many terms to keep in the Dirichlet kernel (~16 is usually fine) */
  REAL8 upsampling;	/**< frequency-upsampling applied to SFTs ==> dFreq != 1/Tsft ... */
  SSBprecision SSBprec; /**< whether to use full relativist SSB-timing, or just simple Newtonian */
  BOOLEAN useRAA;        /**< whether to use the frequency- and sky-position-dependent rigid adiabatic response tensor and not just the long-wavelength approximation */
  BOOLEAN bufferedRAA;	/**< approximate RAA by assuming constant response over (small) frequency band */
  ComputeFBuffer_RS *buffer; /**< buffer for storing pre-resampled timeseries (used for resampling implementation) */
  const EphemerisData *edat;   /**< ephemeris data for re-computing multidetector states */
  BOOLEAN returnAtoms;	/**< whether or not to return the 'FstatAtoms' used to compute the F-statistic */
  BOOLEAN returnSingleF; /**< in multi-detector case, whether or not to also return the single-detector Fstats computed from the atoms */
} ComputeFParams;


/**
 * Struct holding buffered ComputeFStat()-internal quantities to avoid unnecessarily
 * recomputing things that depend ONLY on the skyposition and detector-state series (but not on the spins).
 * For the first call of ComputeFStat() the pointer-entries should all be NULL.
 */
typedef struct tagComputeFBuffer {
  const MultiDetectorStateSeries *multiDetStates;/**< buffer for each detStates (store pointer) and skypos */
  REAL8 Alpha, Delta;				/**< skyposition of candidate */
  MultiSSBtimes *multiSSB;
  MultiSSBtimes *multiBinary;
  MultiAMCoeffs *multiAMcoef;
  MultiCmplxAMCoeffs *multiCmplxAMcoef;
} ComputeFBuffer;

  /** Struct containing vectors of multi- and single-IFO F-stats over a frequency range and full search parameter info in dopplerParams */
typedef struct tagMultiFstatFrequencySeries {
  PulsarDopplerParams doppler;	/**< full info about {sky position, fkdot, refTime, .. and *frequency band*} for which these F values are computed */
  REAL4Vector *F;		/**< 1D array of multi-IFO  F-stat values over {frequencies} */
  REAL4VectorSequence *FX;	/**< 2D array of single-IFO F-stat values over {detectors, frequencies}, ordered as (det1bin1,det1bin2,..,det1binN,det2bin1,...detMbinN) */
} MultiFstatFrequencySeries;

/* macro to index arrays in the MultiFstatFrequencySeries->FX structure */
#define FX_INDEX(FX, iDet, iFreq)          \
  ( ( (iDet) * (FX)->vectorLength ) + (iFreq) )

/*---------- exported Global variables ----------*/
/* empty init-structs for the types defined in here */
extern const Fcomponents empty_Fcomponents;
extern const ComputeFParams empty_ComputeFParams;
extern const ComputeFBuffer empty_ComputeFBuffer;

/*---------- exported prototypes [API] ----------*/
int
XLALComputeFaFb ( Fcomponents *FaFb,
		  const SFTVector *sfts,
		  const PulsarSpins fkdot,
		  const SSBtimes *tSSB,
		  const AMCoeffs *amcoe,
		  const ComputeFParams *params);

int
XLALComputeFaFbXavie ( Fcomponents *FaFb,
		       const SFTVector *sfts,
		       const PulsarSpins fkdot,
		       const SSBtimes *tSSB,
		       const AMCoeffs *amcoe,
		       const ComputeFParams *params);
int
XLALComputeFaFbCmplx ( Fcomponents *FaFb,
		       const SFTVector *sfts,
		       const PulsarSpins fkdot,
		       const SSBtimes *tSSB,
		       const CmplxAMCoeffs *amcoe,
		       const ComputeFParams *params);


void ComputeFStat ( LALStatus *, Fcomponents *Fstat,
		    const PulsarDopplerParams *doppler,
		    const MultiSFTVector *multiSFTs,
		    const MultiNoiseWeights *multiWeights,
		    const MultiDetectorStateSeries *multiDetStates,
		    const ComputeFParams *params,
		    ComputeFBuffer *cfBuffer );

void ComputeFStatFreqBand ( LALStatus *status,
			    REAL4FrequencySeries *FstatVector,
			    const PulsarDopplerParams *doppler,
			    const MultiSFTVector *multiSFTs,
			    const MultiNoiseWeights *multiWeights,
			    const MultiDetectorStateSeries *multiDetStates,
			    const ComputeFParams *params);

int XLALAmplitudeParams2Vect ( PulsarAmplitudeVect A_Mu, const PulsarAmplitudeParams Amp );
int XLALAmplitudeVect2Params ( PulsarAmplitudeParams *Amp, const PulsarAmplitudeVect A_Mu );

/* destructors */
void XLALEmptyComputeFBuffer ( ComputeFBuffer *cfb );

/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
