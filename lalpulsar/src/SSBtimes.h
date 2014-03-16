/*
 *  Copyright (C) 2009 Chris Messenger
 *  Copyright (C) 2005 Reinhard Prix
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
#ifndef _SSBTIMES_H  /* Double-include protection. */
#define _SSBTIMES_H

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup SSBtimes_h Header SSBtimes.h
 * \ingroup pkg_pulsarCoh
 * \author Reinhard Prix
 * \date 2005
 * \brief Functions for working with SSB times.
 */
/*@{*/

/*---------- exported INCLUDES ----------*/
#include <lal/LALStdlib.h>
#include <lal/PulsarDataTypes.h>
#include <lal/DetectorStates.h>

/*---------- exported types ----------*/

/** The precision in calculating the barycentric transformation */
typedef enum {
  SSBPREC_NEWTONIAN,		/**< simple Newtonian: \f$\tau = t + \vec{r}\cdot\vec{n}/c\f$ */
  SSBPREC_RELATIVISTIC,		/**< detailed relativistic: \f$\tau=\tau(t; \vec{n}, \vec{r})\f$ */
  SSBPREC_RELATIVISTICOPT,  	/**< optimized relativistic, numerically equivalent to #SSBPREC_RELATIVISTIC, but faster */
  SSBPREC_LAST			/**< end marker */
} SSBprecision;

/** Simple container for two REAL8-vectors, namely the SSB-timings DeltaT_alpha  and Tdot_alpha,
 * with one entry per SFT-timestamp. These are required input for XLALNewDemod().
 * We also store the SSB reference-time tau0.
 */
typedef struct tagSSBtimes {
  LIGOTimeGPS refTime;		/**< reference-time 'tau0' */
  REAL8Vector *DeltaT;		/**< Time-difference of SFT-alpha - tau0 in SSB-frame */
  REAL8Vector *Tdot;		/**< dT/dt : time-derivative of SSB-time wrt local time for SFT-alpha */
} SSBtimes;

/** Multi-IFO container for SSB timings */
typedef struct tagMultiSSBtimes {
  UINT4 length;		/**< number of IFOs */
  SSBtimes **data;	/**< array of SSBtimes (pointers) */
} MultiSSBtimes;

/*---------- exported Global variables ----------*/
/* empty init-structs for the types defined in here */
extern const SSBtimes empty_SSBtimes;
extern const MultiSSBtimes empty_MultiSSBtimes;

/*---------- exported prototypes [API] ----------*/
int XLALAddBinaryTimes ( SSBtimes **tSSBOut, const SSBtimes *tSSBIn, const PulsarDopplerParams *Doppler );
int XLALAddMultiBinaryTimes ( MultiSSBtimes **multiSSBOut, const MultiSSBtimes *multiSSBIn, const PulsarDopplerParams *Doppler );
SSBtimes *XLALDuplicateSSBtimes ( const SSBtimes *tSSB );
MultiSSBtimes *XLALDuplicateMultiSSBtimes ( const MultiSSBtimes *multiSSB );

SSBtimes *XLALGetSSBtimes ( const DetectorStateSeries *DetectorStates, SkyPosition pos, LIGOTimeGPS refTime, SSBprecision precision );
MultiSSBtimes *XLALGetMultiSSBtimes ( const MultiDetectorStateSeries *multiDetStates, SkyPosition skypos, LIGOTimeGPS refTime, SSBprecision precision);

int XLALEarliestMultiSSBtime ( LIGOTimeGPS *out, const MultiSSBtimes *multiSSB, const REAL8 Tsft );
int XLALLatestMultiSSBtime ( LIGOTimeGPS *out, const MultiSSBtimes *multiSSB,  const REAL8 Tsft );

/* destructors */
void XLALDestroyMultiSSBtimes ( MultiSSBtimes *multiSSB );

/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
