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
 * 
 */



/**
 * \file DriveHoughFStat.h
 * \brief Header file for DriveHoughFStat.c 
 * \author Badri Krishnan, Alicia Sintes
 * \date $Date$
 * 
 ****/


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
#include <glob.h>
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

#include <lalapps.h>

/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif


/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID( DRIVEHOUGHFSTATH, "$Id$" );

/******************************************************
 *  Error codes and messages.
 */
 
#define DRIVEHOUGHFSTAT_ENORM 0
#define DRIVEHOUGHFSTAT_ESUB  1
#define DRIVEHOUGHFSTAT_EARG  2
#define DRIVEHOUGHFSTAT_EBAD  3
#define DRIVEHOUGHFSTAT_EFILE 4
#define DRIVEHOUGHFSTAT_ENULL 5

#define DRIVEHOUGHFSTAT_MSGENORM "Normal exit"
#define DRIVEHOUGHFSTAT_MSGESUB  "Subroutine failed"
#define DRIVEHOUGHFSTAT_MSGEARG  "Error parsing arguments"
#define DRIVEHOUGHFSTAT_MSGEBAD  "Bad argument values"
#define DRIVEHOUGHFSTAT_MSGEFILE "Could not create output file"
#define DRIVEHOUGHFSTAT_MSGENULL "Null pointer"




/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */


  /** sequence of SFT vectors -- for each stack */
  typedef struct tagSFTVectorSequence {
    UINT4 length;     /**< number of stacks */
    SFTVector *data; /**< the SFT vectors */
  } SFTVectorSequence;

  /** parameters for calculating Fstatistic for multiple stacks 
      --- this assumes that the output will be used in a Hough or perhaps 
      other semi-coherent method; therefore it also requires some
      Hough parameters */ 
  typedef struct tagFstatStackParams {
    INT4 nStacks;           /**< number of stacks */
    REAL8 tStackAvg;        /**< average timebaseline of each stack */
    REAL8 fBand;            /**< search band -- fstat is actually calculated in a larger band */
    REAL8 fStart;           /**< start search frequency -- fstat calculation actually starts at a somewhat smaller frequency */
    INT4 nfSizeCylinder;    /**< Hough nfSizeCylinder -- used for calculating extra frequency bins */
    INT4 *mCohSft;          /**< number of sfts in each stack */
    LIGOTimeGPSVector *tsStack; /**< mid times of stacks */
    REAL8 refTime;          /**< reference time for pulsar frequency and spndn. */
    INT4 SSBprecision;      /**< precision for transformation from detector to SSB times*/
    INT4 Dterms;            /**< value of Dterms for LALDemod */
    LALDetector detector;   /**< detector */
    EphemerisData *edat;    /**< ephemeris info */ 
    LIGOTimeGPSVector *ts;  /**< timestamp vector for each sft */
    REAL8 alpha;            /**< sky-location -- right acsension */
    REAL8 delta;            /**< sky-location -- declination */
    REAL8Vector *fdot;      /**< vector containing spindown values */
  } FstatStackParams;

  /** parameters for calculating Hough Maps */
  typedef struct tagHoughParams {
    REAL8 houghThr;            /**< number count threshold */
    REAL8 tStart;              /**< start time of first stack */
    INT8 fBinIni;              /**< start frequency */
    INT8 fBinFin;              /**< frequency band */
    INT4 nfSizeCylinder;       /**< cylinder of LUTs */
    LALDetector detector;      /**< detector */
    LIGOTimeGPSVector *ts;     /**< timestamps of mid points of stacks */
    REAL8VectorSequence *vel;  /**< detector velocity for each stack */
    REAL8VectorSequence *pos;  /**< detector position for each stack */
    REAL8 alpha;               /**< right ascension */
    REAL8 delta;               /**< declination */
    REAL8Vector *fdot;         /**< spindown parameters */
    CHAR *outBaseName;         /**< file for writing output -- if chosen */
  } HoughParams;
  

  /** structure for storing candidates produced by Hough search */
  typedef struct tagHoughCandidates {
    INT4 length;        /**< maximum allowed length of vectors */
    INT4 nCandidates;   /**< number of candidates -- must be less than length */
    REAL8 *freq;        /**< frequency */
    REAL8 *alpha;       /**< right ascension */
    REAL8 *delta;       /**< declination */
    REAL8 *fdot;        /**< spindown vector */
    REAL8 *dFreq;       /**< frequency error */
    REAL8 *dAlpha;      /**< alpha error */
    REAL8 *dDelta ;     /**< delta error */
    REAL8 *dFdot;       /**< fdot error */
  } HoughCandidates;



  void ComputeFstatStack (LALStatus *status, 
			  REAL8FrequencySeriesVector *out, 
			  const SFTVectorSequence *sfts, 
			  FstatStackParams *params);

  void ComputeFstatHoughMap (LALStatus *status,
			     HoughCandidates *out,
			     const HOUGHPeakGramVector *pgV,
			     HoughParams *params);

  void FstatVectToPeakGram (LALStatus *status,
			    HOUGHPeakGramVector *pgV,
			    const REAL8FrequencySeriesVector *FstatVect,
			    REAL8  thr);

  void SetUpStacks1( LALStatus *status,
		    SFTVectorSequence *out,
		    const SFTVector *sftVect,
		    INT4 nStacks);

  void SetUpStacks2( LALStatus *status,
		     SFTVectorSequence *out,
		     const SFTVector *sftVect,
		     const LIGOTimeGPSVector *ts,
		     INT4 nStacks);

  void PrintFstat( LALStatus *status,
		   REAL8FrequencySeries *Fstat, 
		   CHAR *fname,
		   INT4 stackIndex);

  void PrintHmap2file(LALStatus *status,
		      HOUGHMapTotal *ht, 
		      CHAR *fnameOut, 
		      INT4 iHmap);

  void GetHoughCandidates(LALStatus *status,
			  HoughCandidates *houghCand,
			  const HOUGHMapTotal *ht,
			  const HOUGHPatchGrid  *patch,
			  const HOUGHDemodPar   *parDem,
			  REAL8 houghThreshold);

  void PrintHoughCandidates(LALStatus *status,
			    HoughCandidates *in,
			    CHAR *fname);

  void GetLoudestFstat(LALStatus *status,
		       REAL8 *max,
		       REAL8FrequencySeries *Fstat);

  void GetFstatCandidates( LALStatus *status,
			   HoughCandidates *cand,
			   const REAL8FrequencySeries *in,
			   REAL8 FstatThr);

  void PrintHoughHistogram(LALStatus *status,
		      UINT4Vector *hist, 
		      CHAR *fnameOut);
  
#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _DRIVEHOUGHFSTAT_H */



