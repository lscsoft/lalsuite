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
 * $Id$
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

/* lalapps includes */
#include <lalapps.h>

/* additional includes */
#include "HierarchSearchGCT.h"

/*---------- exported DEFINES ----------*/

/** \name Error codes */
/*@{*/
#define LINEVETOC_ENULL 	1
#define LINEVETOC_ENONULL 	2
#define LINEVETOC_EINPUT   	3
#define LINEVETOC_EMEM   	4
#define LINEVETOC_EXLAL		5
#define LINEVETOC_EIEEE		6

#define LINEVETOC_MSGENULL 	"Arguments contained an unexpected null pointer"
#define LINEVETOC_MSGENONULL 	"Output pointer is non-NULL"
#define LINEVETOC_MSGEINPUT   	"Invalid input"
#define LINEVETOC_MSGEMEM   	"Out of memory. Bad."
#define LINEVETOC_MSGEXLAL	"XLAL function call failed"
#define LINEVETOC_MSGEIEEE	"Floating point failure"
/*@}*/

/*---------- exported types ----------*/

  /** Type containing multi- and single-detector F statistics and Line Veto statistic */
  typedef struct tagLVcomponents {
    REAL8 TwoF;                           /**< multi-detector F-statistic value */
    REAL8 LV;                             /**< multi-detector Line Veto statistic value */
    REAL8Vector *TwoFX;                   /**< vector of single-detector F-statistic values */
  } LVcomponents;


/*---------- exported Global variables ----------*/
/* empty init-structs for the types defined in here */
extern const LVcomponents empty_LVcomponents;

/*---------- exported prototypes [API] ----------*/
int
XLALComputeExtraStatsForToplist ( toplist_t *list,
				  const MultiSFTVectorSequence *multiSFTs,
				  const MultiNoiseWeightsSequence *multiNoiseWeights,
				  const MultiDetectorStateSeriesSequence *multiDetStates,
				  const ComputeFParams *CFparams,
				  const LIGOTimeGPS refTimeGPS,
				  const LIGOTimeGPS tMidGPS,
				  const BOOLEAN SignalOnly );

int
XLALComputeExtraStatsSemiCoherent ( LVcomponents *lineVeto,
				    const PulsarDopplerParams *dopplerParams,
				    const MultiSFTVectorSequence *multiSFTs,
				    const MultiNoiseWeightsSequence *multiNoiseWeights,
				    const MultiDetectorStateSeriesSequence *multiDetStates,
				    const ComputeFParams *CFparams,
				    const BOOLEAN SignalOnly);

REAL8
XLALComputeFstatFromAtoms ( const Fcomponents *Fstat,
			    const UINT4 X );

REAL8
XLALComputeLineVeto ( const REAL8 TwoF,
		      const REAL8Vector *TwoFX,
		      const REAL8 rhomax,
		      const REAL8Vector *priorX );

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
