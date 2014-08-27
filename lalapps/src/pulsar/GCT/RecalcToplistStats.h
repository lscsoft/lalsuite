/*
 * Copyright (C) 2011 David Keitel
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
 * \defgroup RecalcToplistStats_h Header RecalcToplistStats.h
 * \ingroup pkg_pulsarCommon
 * \author David Keitel
 *
 * \brief Functions to recompute statistics for GCT/hough toplist entries
 */
/*@{*/

#ifndef _RECALCTOPLISTSTATS_H  /* Double-include protection. */
#define _RECALCTOPLISTSTATS_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- exported INCLUDES ----------*/

/* lal includes */
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/ComputeFstat.h>
#include <lal/StringVector.h>
#include <lal/LineRobustStats.h>

/* additional includes */
#include "GCTtoplist.h"
#include "../hough/src2/HoughFStatToplist.h"

/*---------- exported DEFINES ----------*/

/*---------- exported types ----------*/

/** Type containing multi- and single-detector \f$ \mathcal{F} \f$-statistics and line-robust statistic */
typedef struct tagRecalcStatsComponents {
  REAL4 avTwoF;				/**< multi-detector \f$ \mathcal{F} \f$-statistic, averaged over segments */
  REAL4 avTwoFX[PULSAR_MAX_DETECTORS];	/**< fixed-size array of single-detector \f$ \mathcal{F} \f$-statistics, averaged over segments */
  UINT4 numDetectors;			/**< number of detectors, numDetectors=0 should make all code ignore the TwoFX field. */
  REAL4 log10BSGL;			/**< line-robust statistic \f$ \log_{10}B_{\mathrm{SGL}} \f$ */
  REAL4 avTwoFWithoutLoudestSeg;	/**< average \f$ \mathcal{F} \f$-stat with loudest segment removed */
  INT4 loudestSeg;			/**< index of the segment with the largest contribution to avTwoF */
} RecalcStatsComponents;

/** Type containing input arguments for XLALComputeExtraStatsForToplist() */
typedef struct tagRecalcStatsParams {
  const char *listEntryTypeName;		/**< type of toplist entries, give as name string */
  FstatInputVector *Fstat_in_vec;	/**< vector of input data for XLALComputeFstat() */
  LALStringVector *detectorIDs;		/**< detector name vector with all detectors present in any data sements */
  LIGOTimeGPSVector *startTstack;	/**< starting GPS time of each stack */
  LIGOTimeGPS refTimeGPS;		/**< reference time for fkdot values in toplist */
  BSGLSetup *BSGLsetup;			/**< pre-computed setup for line-robust statistic BSGL */
  BOOLEAN loudestSegOutput;		/**< return extra info about loudest segment */
} RecalcStatsParams;

/*---------- exported Global variables ----------*/
/* empty init-structs for the types defined in here */

/*---------- exported prototypes [API] ----------*/
int
XLALComputeExtraStatsForToplist ( toplist_t *list,
				  const RecalcStatsParams *recalcParams
				);

int
XLALComputeExtraStatsSemiCoherent ( RecalcStatsComponents *stats,
				    const PulsarDopplerParams *dopplerParams,
				    const RecalcStatsParams *recalcParams
				  );

// @}

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
