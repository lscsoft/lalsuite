/*
 * Copyright (C) 2011-2014 David Keitel
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
 * \defgroup LineRobustStats_h Header LineRobustStats.h
 * \ingroup pkg_pulsarCommon
 * \author David Keitel
 *
 * \brief Functions to compute line-robust CW statistics
 */
/*@{*/

#ifndef _LINEROBUSTSTATS_H  /* Double-include protection. */
#define _LINEROBUSTSTATS_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- exported INCLUDES ----------*/

/* lal includes */
#include <lal/LALStdlib.h>
#include <lal/PulsarDataTypes.h>
#include <lal/StringVector.h>
#include <lal/LALConstants.h>
#include <math.h>

/* additional includes */

/*---------- exported DEFINES ----------*/

/*---------- exported types ----------*/

/** Type containing multi- and single-detector F statistics and Line Veto statistic */
typedef struct tagLRScomponents {
  REAL4 TwoF;				/**< multi-detector F-statistic value */
  REAL4 TwoFX[PULSAR_MAX_DETECTORS];	/**< fixed-size array of single-detector F-statistic values */
  UINT4 numDetectors;			/**< number of detectors, numDetectors=0 should make all code ignore the TwoFX field. */
  REAL4 LRS;				/**< multi-detector line-robust statistic value */
} LRScomponents;

/*---------- exported Global variables ----------*/

/*---------- exported prototypes [API] ----------*/
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

// @}

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
