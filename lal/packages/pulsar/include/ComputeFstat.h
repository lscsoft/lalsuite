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
 * \file 
 * \ingroup pulsar
 * \brief Header-file defining the API for the F-statistic functions.
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

/** Multi-IFO container for antenna-pattern coefficients a^X(t), b^X(t) */
typedef struct {
  UINT4 length;		/**< number of IFOs */
  AMCoeffs **data;	/**< array of (inverse-noise weighed) amcoeffs 
			   \f$\sqrt{w_{Xi}}\,a_{Xi}\f$, and \f$\sqrt{w_{Xi}}\,b_{Xi}\f$ */
  REAL8 A;		/**< multi-IFO inverse-noise weighed antenna-pattern coefficient 
			   \f$A = \sum_{X,i} w_{Xi} \, a_{Xi}^2\f$ */
  REAL8 B;		/**<  \f$B = \sum_{X,i} w_{Xi} \, b_{Xi}^2\f$ */
  REAL8 C;		/**<  \f$C = \sum_{X,i} w_{Xi} \, a_{Xi}\, b_{Xi}\f$ */
  REAL8 D;		/**< determinant \f$D = A\,B - C^2\f$ */
} MultiAMCoeffs;

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
  SSBprecision SSBprec; /**< wether to use full relativist SSB-timing, or just simple Newtonian */
} ComputeFParams;

/** Struct holding buffered ComputeFStat()-internal quantities to avoid unnecessarily 
 * recomputing things that depend ONLY on the skyposition and detector-state series (but not on the spins). 
 * For the first call of ComputeFStat() the pointer-entries should all be NULL.
 */
typedef struct {
  const MultiDetectorStateSeries *multiDetStates;/**< buffer for each detStates (store pointer) and skypos */
  REAL8 Alpha, Delta;				/**< skyposition of candidate */
  MultiSSBtimes *multiSSB;	
  MultiAMCoeffs *multiAMcoef;
} ComputeFBuffer;

/*---------- exported Global variables ----------*/

/*---------- exported prototypes [API] ----------*/
int
XLALComputeFaFb ( Fcomponents *FaFb,
		  const SFTVector *sfts, 
		  const PulsarSpins fkdot,
		  const SSBtimes *tSSB,
		  const AMCoeffs *amcoe,
		  UINT4 Dterms);

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
			    REAL8FrequencySeries *FstatVector,
			    const PulsarDopplerParams *doppler,
			    const MultiSFTVector *multiSFTs, 
			    const MultiNoiseWeights *multiWeights,
			    const MultiDetectorStateSeries *multiDetStates,
			    const ComputeFParams *params);

int
XLALWeighMultiAMCoeffs (  MultiAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights );

void
LALEstimatePulsarAmplitudeParams (LALStatus *,
				  PulsarAmplitudeParams *Amp,
				  PulsarAmplitudeParams *dAmp,
				  const Fcomponents *Fstat,
				  const MultiAMCoeffs *multiAMcoef,
				  REAL8 TsftShat);

/* destructors */
void XLALDestroyMultiSSBtimes ( MultiSSBtimes *multiSSB );
void XLALDestroyMultiAMCoeffs ( MultiAMCoeffs *multiAMcoef );

void XLALEmptyComputeFBuffer ( ComputeFBuffer cfb );


/* helpers */
int sin_cos_LUT (REAL4 *sinx, REAL4 *cosx, REAL8 x); 
int sin_cos_2PI_LUT (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x);


#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
