/*
*  Copyright (C) 2007 Badri Krishnan, Alicia Sintes Olives
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
 * \file PeakSelect.h
 * \author Sintes, A.M., and Krishnan, B.
 * \brief Header file for  PeakSelect.c
 *
 * History:   Created by Sintes May 21, 2003
 * Modified by Krishnan Oct 2005
 *
 * ### Header \ref PeakSelect.h ###
 *
 * From periodogram to peakgram
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/PeakSelect.h>
 * \endcode
 */

/*
 * 4.  Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _PEAKSELECT_H
#define _PEAKSELECT_H

/*
 * 5. Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALRunningMedian.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/PHMD.h>
#include "./SFTbin.h"

/*
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

/*
 * 7. Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */

/**\name Error Codes */ /*@{*/

#define PEAKSELECTH_ENULL 1
#define PEAKSELECTH_EVAL 5

#define PEAKSELECTH_MSGENULL "Null pointer"
#define PEAKSELECTH_MSGEVAL  "Invalid value"

/*@}*/


/* ******************************************************
 * 8. Macros. But, note that macros are deprecated.
 *    They could be moved to the modules where are needed
 */


/* *******************************************************
 * 9. Constant Declarations. (discouraged)
 */



/* **************************************************************
 * 10. Structure, enum, union, etc., typdefs.
 */

/** Explicit peakgram structure -- 1 if power in bin is above threshold and 0 if below */  
typedef struct tagUCHARPeakGram{ 
  LIGOTimeGPS  epoch; /**< epoch of first series sample */
  REAL8        timeBase; /**< coherent time baseline used to construct peakgram */
  INT4         fminBinIndex; /**< first frequency bin of peakgram */
  INT4         length; /**< number of elements in data */
  INT4         nPeaks; /**< number of peaks selected in data */
  UCHAR       *data;  /**< pointer to the data {0,1}*/
}UCHARPeakGram; 
  
/** structure containing psd and periodogram of a sft -- obsolete -- use LAL functions */
typedef struct tagREAL8PeriodoPSD{ 
  REAL8Periodogram1    psd;
  REAL8Periodogram1    periodogram;
}REAL8PeriodoPSD; 
  
/*
 * 11. Extern Global variables. (discouraged)
 */


/*
 * 12. Functions Declarations (i.e., prototypes).
 */

/** to compute mean power from a periodogram -- obsolete -- use LAL functions in NormalizeSFTRngMed.c */
void LALComputeMeanPower (LALStatus  *status,		/**< pointer to LALStatus structure */
			  REAL8                *mean, /**< mean power */
			  REAL8Periodogram1    *peri /**< periodogram */);

/** select peakgram in white noise -- obsolete -- use LAL functions in NormalizeSFTRngMed.c */
void LALSelectPeakWhiteNoise(LALStatus  *status,
			     UCHARPeakGram        *pg,
			     REAL8                *thr, /*absolute threshold */
			     REAL8Periodogram1    *peri);

/** Compress explicit peak gram */	     
void LALUCHAR2HOUGHPeak(LALStatus  *status,		/**< pointer to LALStatus structure */
			HOUGHPeakGram        *pgOut, /**< compressed peakgram */
			UCHARPeakGram        *pgIn /**< explicit peakgram -- collection of 1s and 0s*/ 
			);

/** Wrapper for LALRunningMedian code -- obsolete -- use LAL functions in NormalizeSFTRngMed.c */
void LALPeriodo2PSDrng (LALStatus  *status,		/**< pointer to LALStatus structure */
			REAL8Periodogram1    *psd, /**< output psd */
			REAL8Periodogram1    *peri, /**< input periodogram */
			INT4                *blocksRNG /**< running median block size */
			);

/** Function for selecting peaks in colored noise -- obsolete -- use LAL functions in NormalizeSFTRngMed.c */
void LALSelectPeakColorNoise(LALStatus  *status,		/**< pointer to LALStatus structure */
			     UCHARPeakGram        *pg,  /**< output peakgram */
			     REAL8                *thr, /**< threshold reltive to psd */
			     REAL8PeriodoPSD      *in /**< input psd and periodogram */
			     );

/** Constructs peakgram from a normalized SFT -- uses standard pulsar data types */
void SFTtoUCHARPeakGram(LALStatus        *status,		/**< pointer to LALStatus structure */
			UCHARPeakGram    *pg, /**< output peakgram */
			const SFTtype    *sft, /**< standard pulsar sft type */  
			REAL8            thr /**< sft power threshold for peak selection */
			);

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _PEAKSELECT_H */

