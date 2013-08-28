/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes
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
#ifndef _STATISTICS_H
#define _STATISTICS_H

#ifdef  __cplusplus
extern "C" {
#endif


/**
 * \author Krishnan, B., Sintes, A.M.
 * \defgroup Statistics_h Statistics
 * \ingroup pkg_pulsarHough
 * \brief Computes statistics of the Hough maps.
 *
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/Statistics.h>
 * \endcode
 *
 * Given a total Hough map, this calculates the maximum number count, minimum
 * number count, average and standard deviation and produces a histogram of the
 * number counts.
 *
 */
/*@{*/


/* *************
 *    Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */
#include<lal/Date.h>
#include<lal/LALDatatypes.h>
#include<lal/HoughMap.h>
#include<lal/LALStdlib.h>
#include<lal/LALConstants.h>


/** \name Error Codes */
/*@{*/
#define STATISTICSH_ENULL 1
#define STATISTICSH_EVAL 2
#define STATISTICSH_MSGENULL "Null Pointer"
#define STATISTICSH_MSGEVAL "Invalid Value"
/*@}*/


/* *****************************************************
 *   Structure, enum, union, etc., typdefs.
 */

/** Structure for storing statistical information about a Hough map */
typedef struct tagHoughStats {
  HoughTT    maxCount;    /**< maximum number count */
  UINT2      maxIndex[2]; /**< loctaion of maximum number count */
  HoughTT    minCount;    /**< minimum number count */
  UINT2      minIndex[2]; /**< location of minimum number count */
  REAL8      avgCount;    /**< average number count */
  REAL8      stdDev;      /**< standard deviation of number counts */
} HoughStats;

/*
 *  Extern Global variables. (discouraged)
 */

/* ***************************************************
 *  Functions Declarations (i.e., prototypes).
 */
/** Calculates max, min, average and standard deviation of Hough number counts */
void LALHoughStatistics(LALStatus      *status,			/**< pointer to LALStatus structure */
		        HoughStats     *out, /**< output containing statistics */
		        HOUGHMapTotal  *in /**< hough map */);

/** Calculates number count histogram */
void LALHoughHistogram(LALStatus       *status,			/**< pointer to LALStatus structure */
		       UINT8Vector     *out,  /**< histogram */
		       HOUGHMapTotal   *in  /**< hough map*/ );

void LALHoughmapMeanVariance( LALStatus     *status,
			      REAL8         *mean,
			      REAL8         *variance,
			      HOUGHMapTotal *in);

/*@}*/

/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif  /* end of double inclusion protection */













