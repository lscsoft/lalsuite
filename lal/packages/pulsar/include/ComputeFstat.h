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

/*---------- exported DEFINES ----------*/

/*----- Error-codes -----*/
#define COMPUTEFSTATC_ENULL 		1
#define COMPUTEFSTATC_ENONULL 		2
#define COMPUTEFSTATC_EINPUT   		3
#define COMPUTEFSTATC_EMEM   		4
#define COMPUTEFSTATC_EXLAL		5

#define COMPUTEFSTATC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATC_MSGENONULL 	"Output pointer is non-NULL"
#define COMPUTEFSTATC_MSGEINPUT   	"Invalid input"
#define COMPUTEFSTATC_MSGEMEM   	"Out of memory. Bad."
#define COMPUTEFSTATC_MSGEXLAL		"XLAL function call failed"

/*---------- exported types ----------*/

/* ----- Output types for LALGetDetectorStates() */
/** State-info about position, velocity and LMST of a detector together 
 * with corresponding EarthState.
 */
typedef struct
{
  LIGOTimeGPS tGPS;		/**< GPS timestamps corresponding to this entry */
  REAL8 rDetector[3];		/**< Cartesian coords of detector position in ICRS J2000. Units=sec */
  REAL8 vDetector[3];		/**< Cart. coords. of detector velocity, in dimensionless units (v/c)*/
  REAL8 LMST;			/**< local mean sidereal time at the detector-location in radians */
  EarthState earthState;	/**< pointer to EarthState information */
} DetectorState;


/** Timeseries of DetectorState's, representing the detector-info at different timestamps.
 * In addition to the standard 'vector'-fields we also store the detector-info in here.
 */
typedef struct
{
  UINT4 length;			/**< total number of entries */
  DetectorState *data;		/**< array of DetectorState entries */
  LALDetector detector;		/**< detector-info corresponding to this timeseries */
} DetectorStateSeries;

/** Multi-IFO time-series of DetectorStates */
typedef struct
{
  UINT4 length;			/**< number of detectors */
  DetectorStateSeries **data;	/**< vector of pointers to DetectorStateSeries */
} MultiDetectorStateSeries;


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
  AMCoeffs **data;	/**< array of amcoeffs (pointers) */
} MultiAMCoeffs;

/** Type containing F-statistic proper plus the two complex amplitudes Fa and Fb (for ML-estimators) */
typedef struct {
  REAL8 F;		/**< F-statistic value */
  COMPLEX16 Fa;		/**< complex amplitude Fa */
  COMPLEX16 Fb;		/**< complex amplitude Fb */
  REAL8 Bstat;		/**< experimental: 'B-statistic' */
} Fcomponents; 

/** The precision in calculating the barycentric transformation */
typedef enum {
  SSBPREC_NEWTONIAN,		/**< simple Newtonian: \f$\tau = t + \vec{r}\cdot\vec{n}/c\f$ */
  SSBPREC_RELATIVISTIC,		/**< detailed relativistic: \f$\tau=\tau(t; \vec{n}, \vec{r})\f$ */
  SSBPREC_LAST			/**< end marker */
} SSBprecision;


/** One point in the general continuous-waves parameters-space */
typedef struct {
  LIGOTimeGPS refTime;		/**< reference-time for these parameters */
  SkyPosition skypos;		/**< skyposition longitude/latitude [in equatorial-coordinates!!] */
  REAL8Vector *fkdot;		/**< *intrinsic* pulsar frequency + spindowns */
  BinaryOrbitParams *binary; 	/**< PLACEHOLDER for future extension: binary-pulsar params */
} CWParamSpacePoint;

/** Extra parameters controlling the actual computation of F */
typedef struct {
  UINT4 Dterms;		/**< how many terms to keep in the Dirichlet kernel (~16 is usually fine) */
  SSBprecision SSBprec; /**< wether to use full relativist SSB-timing, or just simple Newtonian */
} ComputeFParams;

/** Struct holding buffered ComputeFStat()-internal quantities to avoid unnecessarily 
 * recomputing things. For the first call of ComputeFStat() the pointer-entries should all be NULL.
 */
typedef struct {
  SkyPosition skypos;
  MultiSSBtimes *multiSSB;
  MultiAMCoeffs *multiAMcoef;
  REAL8 A;
  REAL8 B;
  REAL8 C;
  REAL8 Dinv;
} ComputeFBuffer;

/*---------- exported Global variables ----------*/

/*---------- exported prototypes [API] ----------*/
int
XLALComputeFaFb ( Fcomponents *FaFb,
		  const SFTVector *sfts, 
		  const REAL8Vector *fkdot,
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
LALGetDetectorStates (LALStatus *, 
		      DetectorStateSeries **DetectorStates,
		      const LIGOTimeGPSVector *timestamps,
		      const LALDetector *detector,
		      const EphemerisData *edat,
		      REAL8 tOffset);

void 
LALGetMultiDetectorStates( LALStatus *, 
			   MultiDetectorStateSeries **mdetStates, 
			   const MultiSFTVector *multiSFTs, 
			   const EphemerisData *edat );

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
		    const CWParamSpacePoint *psPoint,
		    const MultiSFTVector *multiSFTs,
		    const MultiNoiseWeights *multiWeights,
		    const MultiDetectorStateSeries *multiDetStates,
		    const ComputeFParams *params,
		    ComputeFBuffer *cfBuffer );


void LALCreateDetectorStateSeries (LALStatus *, DetectorStateSeries **vect, UINT4 length );

/* destructors */
void XLALDestroyDetectorStateSeries ( DetectorStateSeries *detStates );
void LALDestroyDetectorStateSeries(LALStatus *, DetectorStateSeries **vect );
void XLALDestroyMultiDetectorStateSeries ( MultiDetectorStateSeries *mdetStates );
void XLALDestroyMultiSSBtimes ( MultiSSBtimes *multiSSB );
void XLALDestroyMultiAMCoeffs ( MultiAMCoeffs *multiAMcoef );

void XLALEmptyComputeFBuffer ( ComputeFBuffer cfb );

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
