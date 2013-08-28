/*
*  Copyright (C) 2007 Gregory Mendell
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
 * \author G. Mendell
 */

/* REVISIONS: */
/* 04/12/05 gam; Move from using StackSlideOld to using StackSlide function. */
/* 04/12/05 gam; add RunStackSlideIsolatedMonteCarloSimulation to StackSlideIsolated.c */
/* 05/24/05 gam; make maxPower and totalEventCount part of params; change finishSUMs to finishPeriodicTable; end xml in FinalizeSearch */
/* 07/15/2005 gam; Change RunStackSlideIsolatedMonteCarloSimulation to inject at a */
/*                 random point in the parameters space a random number of times.  */
/* 07/29/05 gam; if (params->testFlag & 64) > 0 set searchSurroundingPts == 1 and  */
/*               search surrounding parameters space pts; else search nearest only */
/* 09/01/05 gam; If params->numMCRescalings > 0 use this and params->rescaleMCFractionSFTs to */
/*               rescale SFTs to run numMCRescalings Monte Carlo simulations in parallel.     */
/* 01/12/06 gam; Add function WriteStackSlideLEsToPriorResultsFile; if ((outputEventFlag & 8) > 0) write loudest events to params->priorResultsFile */
/* 01/12/06 gam; if ( (outputEventFlag & 8) > 0 ) and running MC, produce estimated UL based on loudest event from priorResultsFile */
/* 02/23/06 gam; add function EstimateStackSlidePower */

#ifndef _STACKSLIDEISOLATED_H
#define _STACKSLIDEISOLATED_H

/*********************************************/
/*                                           */
/* START SECTION: include header files       */
/*                                           */
/*********************************************/
#include <stdio.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALDemod.h>
/* next two are for xml I/O */
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLHeaders.h>
/* next is needed for tables defined in LAL */
#include <lal/LIGOMetadataTables.h>
/* 04/12/05 gam; next two are needed to inject signals for Monte Carlo simulations. */
#include <lal/GeneratePulsarSignal.h>
#include <lal/Random.h>
#include <gsl/gsl_sort.h>
#include "DriveStackSlide.h"
#include "StackSlide.h"
/* #include <lal/LALStackSlide.h> Will need to switch to this version when StackSlide is in LAL. */
/*********************************************/
/*                                           */
/* END SECTION: include header files         */
/*                                           */
/*********************************************/

#ifdef __cplusplus
extern "C" {
#endif
  

/*********************************************/
/*                                           */
/* START SECTION: define constants           */
/*                                           */
/*********************************************/
#define STACKSLIDEISOLATEDH_ENULL 1
#define STACKSLIDEISOLATEDH_ENNUL 2
#define STACKSLIDEISOLATEDH_ENEGA 3
#define STACKSLIDEISOLATEDH_EBADRESULTSFILE 4
#define STACKSLIDEISOLATEDH_EMCEVENTTHRESHOLDFLAGS 5
#define STACKSLIDEISOLATEDH_EIFO 6
#define STACKSLIDEISOLATEDH_ENOFREQBINS 7
#define STACKSLIDEISOLATEDH_ENUMMCRESCALINGS 8
#define STACKSLIDEISOLATEDH_E2NUMMCRESCALINGS 9
#define STACKSLIDEISOLATEDH_EMAXMC 10
#define STACKSLIDEISOLATEDH_EMAXMCERR 11
#define STACKSLIDEISOLATEDH_EESTUL 12
#define STACKSLIDEISOLATEDH_MSGENULL "Null Pointer"
#define STACKSLIDEISOLATEDH_MSGENNUL "Non-Null Pointer"
#define STACKSLIDEISOLATEDH_MSGENEGA "Bad Negative Value"
#define STACKSLIDEISOLATEDH_MSGEBADRESULTSFILE "Could not open priorResultsFile"
#define STACKSLIDEISOLATEDH_MSGEMCEVENTTHRESHOLDFLAGS "Monte Carlo needs (outputEventFlag and 2) > 0 and thresholdFlag <= 0"
#define STACKSLIDEISOLATEDH_MSGEIFO "Invalid or null IFO"
#define STACKSLIDEISOLATEDH_MSGENOFREQBINS "Avoiding or cleaning all freq bins! Monte Carlo cannot inject any signals!"
#define STACKSLIDEISOLATEDH_MSGENUMMCRESCALINGS "The command parameter numMCRescalings cannot be negative."
#define STACKSLIDEISOLATEDH_MSGE2NUMMCRESCALINGS "The command parameter numMCRescalings cannot be > 0 if not reporting MC results."
#define STACKSLIDEISOLATEDH_MSGEMAXMC "The number of MC iterations must be at least 1."
#define STACKSLIDEISOLATEDH_MSGEMAXMCERR "The absolute maximum MC error must be a positive number."
#define STACKSLIDEISOLATEDH_MSGEESTUL "Finding estimated UL failed because d2 for false dismisal bin or LE was less than 1"
/*********************************************/
/*                                           */
/* END SECTION: define constants             */
/*                                           */
/*********************************************/

void StackSlideIsolated (
    LALStatus                        *status,
    SnglStackSlidePeriodicTable      *loudestPeaksArray,
    LALFindStackSlidePeakOutputs     *pLALFindStackSlidePeakOutputs,
    LALFindStackSlidePeakParams      *pLALFindStackSlidePeakParams,
    LALUpdateLoudestStackSlideParams *pLALUpdateLoudestStackSlideParams,
    LALDetector                      *cachedDetector,
    StackSlideParams                 *stksldParams,
    StackSlideSearchParams           *params
);

/* 07/15/05 gam; rewrite to choose random points in the parameter space */
void RunStackSlideIsolatedMonteCarloSimulation(LALStatus *status, StackSlideSearchParams *params, INT4 nSamples);

/* 01/12/06 gam; use the powers from the first MC to generate an estimated UL. */
void FindStackSlideEstimatedUL(
                LALStatus *status,
                REAL8 *h0UpperLimit, 
                REAL8 *uncertainty, 
                REAL8 *power, 
                INT4 powerLength,
                REAL8 confidence,
                REAL8 loudestEvent,
                REAL8 injectedH0, 
                INT4 numBLKs);

/* 01/12/06 gam; Add function WriteStackSlideLEsToPriorResultsFile; write loudest events to params->priorResultsFile */
void WriteStackSlideLEsToPriorResultsFile(
     LALStatus *status, 
     CHAR *priorResultsFile, 
     SnglStackSlidePeriodicTable *loudestPeaksArray, 
     INT4 keepThisNumber,
     CHAR *IFO, 
     REAL8 startFreq, 
     REAL8 searchBand,
     REAL8 confidence,
     REAL8 h0UpperLimit,
     REAL8 uncertainty
);

/* 05/24/05 gam; Function that reads results from previous jobs in the pipeline */
void getStackSlidePriorResults(LALStatus *status,
                               REAL4 *priorLoudestEvent,
                               REAL8 *priorStartFreq,
                               REAL8 *priorBand,
                               REAL8 *priorConfidence,
                               REAL8 *priorUL,
                               REAL8 *priorUncertainty,
                               CHAR  *priorResultsFile);

/* 07/15/2005 gam; Change RunStackSlideIsolatedMonteCarloSimulation */
void StackSlideGetUniformDeviate(LALStatus *status, REAL8 *returnVal, REAL8 startVal, REAL8 range, RandomParams *randPar);

/* 09/01/05 gam */
void StackSlideULLeastSquaresLinFit(REAL8 *interpolatedUL, const REAL8 *arrayULs, const REAL8 *arrayConfs, REAL8 desiredConf, INT4 startIndex, INT4 numVals);

void EstimateStackSlidePower(LALStatus *status,
			REAL4FrequencySeries **SUMData,
			REAL4FrequencySeries **STKData,
			TdotsAndDeltaTs *pTdotsAndDeltaTs,
			StackSlideSearchParams *searchParams,
			LALDetector          *cachedDetector,
			INT4 iSky,
			StackSlideParams *params);

#ifdef __cplusplus
}
#endif

#endif /* _STACKSLIDEISOLATED_H */
