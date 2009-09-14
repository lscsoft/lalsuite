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

/**
 * \author Reinhard Prix
 * \date 2005
 * \ingroup pulsarCoherent
 * \file
 * \brief Header-file defining the API for the F-statistic functions.
 *
 * This code is (partly) a descendant of an earlier implementation found in
 * LALDemod.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens, Bruce Allen
 * ComputSky.[ch] by Jolien Creighton, Reinhard Prix, Steve Berukoff
 * LALComputeAM.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens
 *
 * $Id$
 *
 */

#ifndef _COMPUTEFSTAT_H  /* Double-include protection. */
#define _COMPUTEFSTAT_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

#include <lal/LALRCSID.h>
NRCSID( COMPUTEFSTATH, "$Id$" );

/*---------- exported INCLUDES ----------*/
#include <lal/LALComputeAM.h>
#include <lal/PulsarDataTypes.h>
#include <lal/DetectorStates.h>

/*---------- exported DEFINES ----------*/

/*----- Error-codes -----*/
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

/*---------- exported types ----------*/

/** Simple container for two REAL8-vectors, namely the SSB-timings DeltaT_alpha  and Tdot_alpha,
 * with one entry per SFT-timestamp. These are required input for XLALNewDemod().
 * We also store the SSB reference-time tau0.
 */
typedef struct {
  LIGOTimeGPS refTime;
  REAL8Vector *DeltaT;		/**< Time-difference of SFT-alpha - tau0 in SSB-frame */
  REAL8Vector *Tdot;		/**< dT/dt : time-derivative of SSB-time wrt local time for SFT-alpha */
} SSBtimes;

/** Multi-IFO container for SSB timings */
typedef struct {
  UINT4 length;		/**< number of IFOs */
  SSBtimes **data;	/**< array of SSBtimes (pointers) */
} MultiSSBtimes;


/** Struct holding the "antenna-pattern" matrix \f$\mathcal{M}_{\mu\nu} \equiv \left( \mathbf{h}_\mu|\mathbf{h}_\nu\right)\f$,
 * in terms of the multi-detector scalar product. This matrix can be shown to be expressible as
 * \f{equation}
 * \mathcal{M}_{\mu\nu} = \mathcal{S}^{-1}\,T_\mathrm{SFT}\,\left( \begin{array}{c c c c} A_d & C_d & 0 & 0 \\ C_d & B_d & 0 & 0 \\ 0 & 0 & A_d & C_d \\ 0 & 0 & C_d & B_d \\ \end{array}\right)\,,
 * \f}
 * where (here) \f$\mathcal{S} \equiv \frac{1}{N_\mathrm{SFT}}\sum_{X,\alpha} S_{X\alpha}\f$ characterizes the (single-sided!)
 * multi-detector noise-floor, and
 * \f{equation}
 * A_d \equiv \sum_{X,\alpha} \widehat{a}^X_\alpha \widehat{a}^X_\alpha\,,\quad
 * B_d \equiv \sum_{X,\alpha} \widehat{b}^X_\alpha \widehat{b}^X_\alpha \,,\quad
 * C_d \equiv \sum_{X,\alpha} \widehat{a}^X_\alpha \widehat{b}^X_\alpha \,,
 * \f}
 * and the noise-weighted atenna-functions \f$\widehat{a}^X_\alpha = \sqrt{w^X_\alpha}\,a^X_\alpha\f$,
 * \f$\widehat{b}^X_\alpha = \sqrt{w^X_\alpha}\,b^X_\alpha\f$, and noise-weights
 * \f$w^X_\alpha \equiv {S^{-1}_{X\alpha}/{\mathcal{S}^{-1}}\f$.
 *
 * \note One reason for storing the un-normalized \a Ad, \a Bd, \a Cd and the normalization-factor \a Sinv_Tsft separately
 * is that the former are of order unity, while \a Sinv_Tsft is very large, and it has numerical advantages for parameter-estimation
 * to use that fact.
 */
typedef struct {
  REAL8 Ad; 		/**<  \f$A_d \equiv \sum_{X,\alpha} \widehat{a}^X_\alpha \widehat{a}^X_\alpha\f$ */
  REAL8 Bd; 		/**<  \f$B_d \equiv \sum_{X,\alpha} \widehat{b}^X_\alpha \widehat{b}^X_\alpha\f$ */
  REAL8 Cd; 		/**<  \f$C_d \equiv \sum_{X,\alpha} \widehat{a}^X_\alpha \widehat{b}^X_\alpha\f$ */
  REAL8 Dd; 		/**<  determinant \f$D_d \equiv A_d B_d - C_d^2 \f$ */
  REAL8 Sinv_Tsft;	/**< normalization-factor \f$\mathcal{S}^{-1}\,T_\mathrm{SFT}\f$ (wrt single-sided PSD!) */
} AntennaPatternMatrix;


/** Multi-IFO container for antenna-pattern coefficients a^X(t), b^X(t) and atenna-pattern matrix M_mu_nu */
typedef struct {
  UINT4 length;		/**< number of IFOs */
  AMCoeffs **data;	/**< noise-weighted am-coeffs \f$\widehat{a}_{X\alpha}\f$, and \f$\widehat{b}_{X\alpha}\f$ */
  AntennaPatternMatrix Mmunu;	/**< antenna-pattern matrix \f$\mathcal{M}_{\mu\nu}\f$ */
} MultiAMCoeffs;

/** Struct holding the "antenna-pattern" matrix \f$\mathcal{M}_{\mu\nu} \equiv \left( \mathbf{h}_\mu|\mathbf{h}_\nu\right)\f$,
 * in terms of the multi-detector scalar product. This matrix can be shown to be expressible, in the case of complex AM co\"{e}fficients, as
 * \f{equation}
 * \mathcal{M}_{\mu\nu} = \mathcal{S}^{-1}\,T_\mathrm{SFT}\,\left( \begin{array}{c c c c} A_d & C_d & 0 & -E_d \\ C_d & B_d & E_d & 0 \\ 0 & E_d & A_d & C_d \\ -E_d & 0 & C_d & B_d \\ \end{array}\right)\,,
 * \f}
 * where (here) \f$\mathcal{S} \equiv \frac{1}{N_\mathrm{SFT}}\sum_{X,\alpha} S_{X\alpha}\f$ characterizes the multi-detector noise-floor, and
 * \f{equation}
 * A_d \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{a}^X_\alpha\,,\quad
 * B_d \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{b}^X_\alpha{}^* \widehat{b}^X_\alpha \,,\quad
 * C_d \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha \,,
 * E_d \equiv \mathrm{Im} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha \,,
 * \f}
 * and the noise-weighted atenna-functions \f$\widehat{a}^X_\alpha = \sqrt{w^X_\alpha}\,a^X_\alpha\f$,
 * \f$\widehat{b}^X_\alpha = \sqrt{w^X_\alpha}\,b^X_\alpha\f$, and noise-weights
 * \f$w^X_\alpha \equiv {S^{-1}_{X\alpha}/{\mathcal{S}^{-1}}\f$.
 *
 * \note One reason for storing the un-normalized \a Ad, \a Bd, \a Cd, \a Ed and the normalization-factor \a Sinv_Tsft separately
 * is that the former are of order unity, while \a Sinv_Tsft is very large, and it has numerical advantages for parameter-estimation
 * to use that fact.
 */
typedef struct {
  REAL8 Ad; 		/**<  \f$A_d \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{a}^X_\alpha\f$ */
  REAL8 Bd; 		/**<  \f$B_d \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{b}^X_\alpha{}^* \widehat{b}^X_\alpha\f$ */
  REAL8 Cd; 		/**<  \f$C_d \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha\f$ */
  REAL8 Ed; 		/**<  \f$E_d \equiv \mathrm{Im} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha\f$ */
  REAL8 Dd; 		/**<  determinant \f$D_d \equiv A_d B_d - C_d^2 -E_d^2 \f$ */
  REAL8 Sinv_Tsft;	/**< normalization-factor \f$\mathcal{S}^{-1}\,T_\mathrm{SFT}\f$ (wrt single-sided PSD!) */
} CmplxAntennaPatternMatrix;

/** Multi-IFO container for antenna-pattern coefficients a^X(t), b^X(t) and atenna-pattern matrix M_mu_nu */
typedef struct {
  UINT4 length;		/**< number of IFOs */
  CmplxAMCoeffs **data;	/**< noise-weighted am-coeffs \f$\widehat{a}_{X\alpha}\f$, and \f$\widehat{b}_{X\alpha}\f$ */
  CmplxAntennaPatternMatrix Mmunu;	/**< antenna-pattern matrix \f$\mathcal{M}_{\mu\nu}\f$ */
} MultiCmplxAMCoeffs;

/** Type containing F-statistic proper plus the two complex amplitudes Fa and Fb (for ML-estimators) */
typedef struct {
  REAL8 F;		/**< F-statistic value */
  COMPLEX16 Fa;		/**< complex amplitude Fa */
  COMPLEX16 Fb;		/**< complex amplitude Fb */
} Fcomponents;

/** The precision in calculating the barycentric transformation */
typedef enum {
  SSBPREC_NEWTONIAN,		/**< simple Newtonian: \f$\tau = t + \vec{r}\cdot\vec{n}/c\f$ */
  SSBPREC_RELATIVISTIC,		/**< detailed relativistic: \f$\tau=\tau(t; \vec{n}, \vec{r})\f$ */
  SSBPREC_LAST			/**< end marker */
} SSBprecision;


/** Extra parameters controlling the actual computation of F */
typedef struct {
  UINT4 Dterms;		/**< how many terms to keep in the Dirichlet kernel (~16 is usually fine) */
  REAL8 upsampling;	/**< frequency-upsampling applied to SFTs ==> dFreq != 1/Tsft ... */
  SSBprecision SSBprec; /**< whether to use full relativist SSB-timing, or just simple Newtonian */
  BOOLEAN useRAA;        /**< whether to use the frequency- and sky-position-dependent rigid adiabatic response tensor and not just the long-wavelength approximation */
  BOOLEAN bufferedRAA;	/**< approximate RAA by assuming constant response over (small) frequency band */
} ComputeFParams;

/** Struct holding buffered ComputeFStat()-internal quantities to avoid unnecessarily
 * recomputing things that depend ONLY on the skyposition and detector-state series (but not on the spins).
 * For the first call of ComputeFStat() the pointer-entries should all be NULL.
 */
typedef struct {
  const MultiDetectorStateSeries *multiDetStates;/**< buffer for each detStates (store pointer) and skypos */
  REAL8 Alpha, Delta;				/**< skyposition of candidate */
  MultiSSBtimes *multiSSB;
  MultiSSBtimes *multiBinary;
  MultiAMCoeffs *multiAMcoef;
  MultiCmplxAMCoeffs *multiCmplxAMcoef;
} ComputeFBuffer;


/*---------- exported Global variables ----------*/
/* empty init-structs for the types defined in here */
extern const SSBtimes empty_SSBtimes;
extern const MultiSSBtimes empty_MultiSSBtimes;
extern const AntennaPatternMatrix empty_AntennaPatternMatrix;
extern const MultiAMCoeffs empty_MultiAMCoeffs;
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

void
LALGetBinarytimes (LALStatus *,
		   SSBtimes *tBinary,
		   const SSBtimes *tSSB,
		   const DetectorStateSeries *DetectorStates,
		   const BinaryOrbitParams *binaryparams,
		   LIGOTimeGPS refTime);

void
LALGetMultiBinarytimes (LALStatus *status,
			MultiSSBtimes **multiBinary,
			const MultiSSBtimes *multiSSB,
			const MultiDetectorStateSeries *multiDetStates,
			const BinaryOrbitParams *binaryparams,
			LIGOTimeGPS refTime);

void
LALGetSSBtimes (LALStatus *,
		SSBtimes *tSSB,
		const DetectorStateSeries *DetectorStates,
		SkyPosition pos,
		LIGOTimeGPS refTime,
		SSBprecision precision);

void
LALGetAMCoeffs(LALStatus *,
	       AMCoeffs *coeffs,
	       const DetectorStateSeries *DetectorStates,
	       SkyPosition skypos);

void
LALNewGetAMCoeffs(LALStatus *,
		  AMCoeffs *coeffs,
		  const DetectorStateSeries *DetectorStates,
		  SkyPosition skypos);


int
XLALComputeAntennaPatternCoeffs ( REAL8 *ai,
				  REAL8 *bi,
				  const SkyPosition *skypos,
				  const LIGOTimeGPS *tGPS,
				  const LALDetector *site,
				  const EphemerisData *edat
				  );

void
LALGetMultiSSBtimes (LALStatus *,
		     MultiSSBtimes **multiSSB,
		     const MultiDetectorStateSeries *multiDetStates,
		     SkyPosition pos,
		     LIGOTimeGPS refTime,
		     SSBprecision precision );

void
LALGetMultiAMCoeffs (LALStatus *,
		     MultiAMCoeffs **multiAMcoef,
		     const MultiDetectorStateSeries *multiDetStates,
		     SkyPosition pos );


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

int
XLALWeighMultiAMCoeffs (  MultiAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights );

void
LALEstimatePulsarAmplitudeParams (LALStatus * status,
				  PulsarCandidate *pulsarParams,
				  const Fcomponents *Fstat,
				  const LIGOTimeGPS *FstatRefTime,
				  const CmplxAntennaPatternMatrix *Mmunu
				  );

/* destructors */
void XLALDestroyMultiSSBtimes ( MultiSSBtimes *multiSSB );
void XLALDestroyMultiAMCoeffs ( MultiAMCoeffs *multiAMcoef );
void XLALDestroyAMCoeffs ( AMCoeffs *amcoef );

void XLALEmptyComputeFBuffer ( ComputeFBuffer *cfb );


/* helpers */
int sin_cos_LUT (REAL4 *sinx, REAL4 *cosx, REAL8 x);
int sin_cos_2PI_LUT (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x);


#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
