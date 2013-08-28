/*  
 *  Copyright (C) 2005-2008 Badri Krishnan, Alicia Sintes, Bernd Machenschalk
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
 * 
 */



/**
 * \file
 * \brief Header file for DriveHoughFStat.c
 * \author Badri Krishnan, Alicia Sintes
 *
 */


#ifndef _DRIVEHOUGHCOLOR_H
#define _DRIVEHOUGHCOLOR_H

/* standard includes */
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h> 

/* lal includes */
#include <lal/UserInput.h>
#include <lal/LALStdlib.h>
#include <lal/PulsarDataTypes.h>
#include <lal/SFTfileIO.h>
#include <lal/AVFactories.h>
#include <lal/RngMedBias.h>
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Velocity.h>
#include <lal/LALDemod.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/Date.h>
#include <lal/LALHough.h> 
#include <lal/NormalizeSFTRngMed.h>
#include <lal/ComputeFstat.h>
#include <lal/Statistics.h>
#include <lal/GeneratePulsarSignal.h> 
#include <lal/LogPrintf.h>
#include <lal/DopplerScan.h>

/* lalapps includes */
#include <lalapps.h>

/* more efficient toplist using heaps */
#include "HoughFStatToplist.h"

/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif


/******************************************************
 *  Error codes and messages.
 */
 
#define HIERARCHICALSEARCH_ENORM 0
#define HIERARCHICALSEARCH_ESUB  1
#define HIERARCHICALSEARCH_EARG  2
#define HIERARCHICALSEARCH_EBAD  3
#define HIERARCHICALSEARCH_EFILE 4
#define HIERARCHICALSEARCH_ENULL 5
#define HIERARCHICALSEARCH_EVAL  6
#define HIERARCHICALSEARCH_ENONULL 7
#define HIERARCHICALSEARCH_EDLOPEN 8
#define HIERARCHICALSEARCH_EWORKER 9
#define HIERARCHICALSEARCH_ECHECKPT 10
#define HIERARCHICALSEARCH_EMEM 11
#define HIERARCHICALSEARCH_ESFT 12


#define HIERARCHICALSEARCH_MSGENORM    "Normal exit"
#define HIERARCHICALSEARCH_MSGESUB     "Subroutine failed"
#define HIERARCHICALSEARCH_MSGEARG     "Error parsing arguments"
#define HIERARCHICALSEARCH_MSGEBAD     "Bad argument values"
#define HIERARCHICALSEARCH_MSGEFILE    "Could not create output file"
#define HIERARCHICALSEARCH_MSGENULL    "Null pointer"
#define HIERARCHICALSEARCH_MSGEVAL     "Invalid value"
#define HIERARCHICALSEARCH_MSGENONULL  "Pointer not null"
#define HIERARCHICALSEARCH_MSGECHECKPT "Could not resume from checkpoint"
#define HIERARCHICALSEARCH_MSGEMEM     "Out of memory"
#define HIERARCHICALSEARCH_MSGESFT     "SFT validity check failed"


/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */


  /** sequence of SFT catalogs -- for each segment */
  typedef struct tagSFTCatalogSequence {
    UINT4 length;     		/**< the number of segments */
    SFTCatalog *data; 		/**< the catalogs */
  } SFTCatalogSequence;



  /** parameters for the semicoherent stage -- hough or stackslide */
  typedef struct tagSemiCoherentParams {
    LIGOTimeGPSVector *tsMid;  /**< timestamps of mid points of stacks */
    LIGOTimeGPS refTime;       /**< reference time for f, fdot definition */
    REAL8VectorSequence *vel;  /**< detector velocity for each stack */
    REAL8VectorSequence *pos;  /**< detector position for each stack */
    REAL8 alpha;               /**< right ascension of demodulation point */
    REAL8 delta;               /**< declination of demodulation point*/
    REAL8 pixelFactor;         /**< Resolution of semicoherent sky-grid */
    REAL8 patchSizeX;          /**< Size of semicoherent sky-patch */
    REAL8 patchSizeY;          /**< Size of semicoherent sky-patch */
    REAL8 fdot;                /**< spindown value of demodulation point */
    UINT4 nfdot;               /**< number of fdot values to search over */ 
    REAL8 dfdot;               /**< resolution in residual spindowns */
    CHAR *outBaseName;         /**< file for writing output -- if chosen */
    BOOLEAN useToplist;        /**< Use a toplist for producing candidates? */
    REAL8  threshold;          /**< Threshold for candidate selection */
    REAL8Vector *weightsV;     /**< Vector of weights for each stack */
    UINT4 extraBinsFstat;      /**< Extra bins required for Fstat calculation */
  } SemiCoherentParams;

  /** one hough or stackslide candidate */
  typedef struct tagSemiCohCandidate {
    REAL8 freq;        /**< frequency */
    REAL8 alpha;       /**< right ascension */
    REAL8 delta;       /**< declination */
    REAL8 fdot;        /**< spindown */
    REAL8 dFreq;       /**< frequency error */
    REAL8 dAlpha;      /**< alpha error */
    REAL8 dDelta ;     /**< delta error */
    REAL8 dFdot;       /**< fdot error */
    REAL8 significance;/**< significance */
    REAL8 alphaBest;   /**< alpha for best candidate in hough map */
    REAL8 deltaBest;   /**< delta for best candidate in hough map */
    REAL8 meanSig;     /**< mean of significance values in hough map */
    REAL8 varianceSig; /**< variance of significance values in Hough map */
  } SemiCohCandidate;  

  /** structure for storing candidates produced by Hough search */
  typedef struct tagSemiCohCandidateList {
    LIGOTimeGPS refTime;       /**< reference time for candidates */
    INT4 length;               /**< maximum allowed length of vectors */
    INT4 nCandidates;          /**< number of candidates -- must be less than length */
    SemiCohCandidate *list;    /**> list of candidates */
  } SemiCohCandidateList;



  /* function prototypes */

  void ComputeFstatHoughMap (LALStatus *status,
			     SemiCohCandidateList *out,
			     HOUGHPeakGramVector *pgV,
			     SemiCoherentParams *params);

  void FstatVectToPeakGram (LALStatus *status,
			    HOUGHPeakGramVector *pgV,
			    REAL4FrequencySeriesVector *FstatVect,
			    REAL4  thr);

  void SetUpStacks(LALStatus *status, 
		 SFTCatalogSequence  *out,  
		 REAL8 tStack,
		 SFTCatalog  *in,
		 UINT4 nStacks);

  void PrintHmap2file(LALStatus *status,
		      HOUGHMapTotal *ht, 
		      CHAR *fnameOut, 
		      INT4 iHmap);

  extern
  void GetHoughCandidates_threshold(LALStatus            *status,
				    SemiCohCandidateList *out,
				    HOUGHMapTotal        *ht,
				    HOUGHPatchGrid       *patch,
				    HOUGHDemodPar        *parDem,
				    REAL8                threshold);

  extern
  void GetHoughCandidates_toplist(LALStatus      *status,
				  toplist_t      *list,
				  HOUGHMapTotal  *ht,
				  HOUGHPatchGrid *patch,
				  HOUGHDemodPar  *parDem);

  void GetFstatCandidates_toplist(LALStatus *status,
				  toplist_t *list,
				  REAL8FrequencySeries   *FstatVec,
				  REAL8 alpha,
				  REAL8 delta,
				  REAL8 fdot);

  void GetChkPointIndex( LALStatus *status,
			 INT4 *loopindex, 
			 const CHAR *fnameChkPoint);


  
#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _HIERARCHICALSEARCH_H */

