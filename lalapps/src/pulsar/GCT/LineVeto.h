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
 * \author David Keitel
 * \date 2011
 * \ingroup pulsarCoherent
 * \file
 * \brief Header-file defining functions related to GCT Line Veto followups
 *
 * This code is partly based on work done by
 * Reinhard Prix, Maria Alessandra Papa, M. Siddiqi
 *
 */

#ifndef _LINEVETO_H  /* Double-include protection. */
#define _LINEVETO_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- exported INCLUDES ----------*/

/* lal includes */
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/ComputeFstat.h>
#include <lal/StringVector.h>

/* lalapps includes */
#include <lalapps.h>

/* additional includes */
#include "GCTtoplist.h"
#include "../hough/src2/HoughFStatToplist.h"

/*---------- exported DEFINES ----------*/

/*---------- exported types ----------*/

/** Type containing multi- and single-detector F statistics and Line Veto statistic */
typedef struct tagLVcomponents {
   REAL4 TwoF;                           /**< multi-detector F-statistic value */
   REAL4 LV;                             /**< multi-detector Line Veto statistic value */
   REAL4Vector *TwoFX;                   /**< vector of single-detector F-statistic values */
} LVcomponents;

/*---------- exported Global variables ----------*/
/* empty init-structs for the types defined in here */
extern const LVcomponents empty_LVcomponents;

/*---------- exported prototypes [API] ----------*/
int
XLALComputeExtraStatsForToplist ( toplist_t *list,
				  const char *listEntryTypeName,
				  const FstatInputDataVector *Fstat_in_vec,
				  const LALStringVector *detectorIDs,
				  const LIGOTimeGPSVector *startTstack,
				  const LIGOTimeGPS refTimeGPS,
				  const char* outputSingleSegStats );

int
XLALComputeExtraStatsSemiCoherent ( LVcomponents *lineVeto,
				    const PulsarDopplerParams *dopplerParams,
				    const FstatInputDataVector *Fstat_in_vec,
				    const LALStringVector *detectorIDs,
				    const LIGOTimeGPSVector *startTstack,
				    FILE *singleSegStatsFile );

REAL4
XLALComputeLineVeto ( const REAL4 TwoF,
		      const REAL4Vector *TwoFX,
		      const REAL8 rhomaxline,
		      const REAL8Vector *lX,
		      const BOOLEAN useAllTerms );

REAL4
XLALComputeLineVetoArray ( const REAL4 TwoF,
                           const UINT4 numDetectors,
                           const REAL4 *TwoFX,
                           const REAL8 logRhoTerm,
                           const REAL8 *loglX,
                           const BOOLEAN useAllTerms );

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
