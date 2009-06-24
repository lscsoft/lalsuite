/*
*  Copyright (C) 2009 Reinhard Prix
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

#ifndef _COMPUTEFSTATGPU_H
#define _COMPUTEFSTATGPU_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/* ---------- includes ---------- */
#include <lal/LALDatatypes.h>
#include <lal/DetectorStates.h>

/* ---------- exported defines and macros ---------- */

/* ---------- exported types ---------- */

/** REAL4 version of pulsar spins fkdot[] */
typedef struct {
  REAL4 FreqMain;			/**< "main" part of frequency fkdot[0], normally just the integral part */
  REAL4 fkdot[PULSAR_MAX_SPINS];	/**< remaining spin-parameters, including *fractional* part of Freq = fkdot[0] */
} PulsarSpins_REAL4;


/** Simple container for REAL4-vectors, holding the SSB-timings DeltaT_alpha  and Tdot_alpha,
 *  with one entry per SFT-timestamp. We also store the SSB reference-time tau0.
 * NOTE: this is a REAL4 version of SSBtimes, preserving the required precision by appropriate
 * 'splitting' of REAL8's into pairs of REAL4s.
 */
typedef struct {
  LIGOTimeGPS refTime;		/**< Reference time wrt to which the time-differences DeltaT are computed */
  REAL4Vector *DeltaT_int;	/**< Integral part of time-difference of SFT-alpha - tau0 in SSB-frame */
  REAL4Vector *DeltaT_rem;	/**< Remainder of time-difference of SFT-alpha - tau0 in SSB-frame */
  REAL4Vector *TdotM1;		/**< dT/dt - 1 : time-derivative of SSB-time wrt local time for SFT-alpha, of order O(1e-4) */
} SSBtimes_REAL4;


/** Multi-IFO container for SSB timings in REAL4-representation */
typedef struct {
  UINT4 length;			/**< number of IFOs */
  SSBtimes_REAL4 **data;	/**< array of SSBtimes (pointers) */
} MultiSSBtimes_REAL4;

/** Struct holding buffered ComputeFStat_REAL4()-internal quantities to avoid unnecessarily
 * recomputing things that depend ONLY on the skyposition and detector-state series (but not on the spins).
 * For the first call of ComputeFStat() the pointer-entries should all be NULL.
 */
typedef struct {
  REAL8 Alpha, Delta;				/**< skyposition of candidate */
  MultiDetectorStateSeries *multiDetStates;	/**< buffer for each detStates (store pointer) and skypos */
  MultiSSBtimes_REAL4 *multiSSB;
  MultiAMCoeffs *multiAMcoef;
} ComputeFBuffer_REAL4;


/** Type containing F-statistic proper plus the two complex amplitudes Fa and Fb (for ML-estimators).
 * NOTE: this is simply a REAL4 version of Fcomponents.
 */
typedef struct {
  REAL4 F;		/**< F-statistic value */
  COMPLEX8 Fa;		/**< complex amplitude Fa */
  COMPLEX8 Fb;		/**< complex amplitude Fb */
} Fcomponents_REAL4;


/* ---------- exported global variables ---------- */
extern const ComputeFBuffer_REAL4 empty_ComputeFBuffer_REAL4;

/* ---------- exported API prototypes ---------- */
void ComputeFStat_REAL4 ( LALStatus *status,
                          Fcomponents *Fstat,
                          const PulsarDopplerParams *doppler,
                          const MultiSFTVector *multiSFTs,
                          const MultiNoiseWeights *multiWeights,
                          const MultiDetectorStateSeries *multiDetStates,
                          UINT4 Dterms,
                          ComputeFBuffer_REAL4 *cfBuffer );

int XLALComputeFaFb_REAL4 ( Fcomponents_REAL4 *FaFb,
                            const SFTVector *sfts,
                            const PulsarSpins_REAL4 spins,
                            const SSBtimes_REAL4 *tSSB,
                            const AMCoeffs *amcoe,
                            UINT4 Dterms);

MultiSSBtimes_REAL4 *
XLALGetMultiSSBtimes_REAL4 ( const MultiDetectorStateSeries *multiDetStates,
                             REAL8 Alpha, REAL8 Delta,
                             LIGOTimeGPS refTime
                             );

SSBtimes_REAL4 *
XLALGetSSBtimes_REAL4 ( const DetectorStateSeries *DetectorStates,
                        REAL8 Alpha, REAL8 Delta,
                        LIGOTimeGPS refTime
                        );

void XLALDestroySSBtimes_REAL4 ( SSBtimes_REAL4 *tSSB );
void XLALDestroyMultiSSBtimes_REAL4 ( MultiSSBtimes_REAL4 *multiSSB );

void XLALEmptyComputeFBuffer_REAL4 ( ComputeFBuffer_REAL4 *cfb);


#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
