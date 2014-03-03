/*
*  Copyright (C) 2007 Gregory Mendell, Virginia Re
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
 * \author  Virginia Re
 */

/* Revisions: */
/* 04/12/05 gam; Add StackSlideSearchParams *params to StackSlideBinary. Need to include "DriveStackSlide.h" */
/* 04/12/05 gam; Remove from StackSlideParams *stksldParams, those already in StackSlideSearchParams *params */


#ifndef _STACKSLIDEBINARY_H
#define _STACKSLIDEBINARY_H

#define STACKSLIDEBINARY_DEBUG 1
#define STACKSLIDEBINARYH_EIFO 6
#define STACKSLIDEBINARYH_ENOFREQBINS 7
#define STACKSLIDEBINARYH_EMAXMC 8
#define STACKSLIDEBINARYH_EDELTASMAnSMA 9
#define STACKSLIDEBINARYH_MSGEIFO "Invalid or null IFO"
#define STACKSLIDEBINARYH_MSGENOFREQBINS "Avoiding or cleaning all freq bins! Monte Carlo cannot inject any signals!"
#define STACKSLIDEBINARYH_MSGEMAXMC "The number of MC iterations must be at least 1."
#define STACKSLIDEBINARYH_MSGDELTASMAnSMA "If params->deltaSMA > 0 you need nMaxSMA > 1"
/*********************************************/
/*                                           */
/* START SECTION: include header files       */
/*                                           */
/*********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <math.h>
#include <string.h> 
#include <lal/LALConstants.h>
#include <lal/StreamInput.h>
#include <lal/SeqFactories.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Random.h>
#include <getopt.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALDatatypes.h>
#include <lal/FindRoot.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/VectorOps.h>

#include <lal/GeneratePulsarSignal.h>

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
  

/* 04/12/05 gam; add StackSlideSearchParams *params */
void StackSlideBinary(  LALStatus *status,
                        StackSlideParams *stksldParams,
                        StackSlideSearchParams *params
                      );

/* 04/12/05 gam */
void FindBinaryLoudest(REAL8 *LoudestEvent, REAL8 *peakFreq, REAL4FrequencySeries **SUMData, StackSlideParams *stksldParams, REAL8 SemimajorAxis, FILE *fpSavedEvents);

/*05/08/07 vir: add function to run MC simulation in the binary case */
void RunStackSlideBinaryMonteCarloSimulation(LALStatus *status, StackSlideSearchParams *params, INT4 nSamples);

void getStackSlideBinaryPriorResults(LALStatus *status,
                               REAL4 *priorLoudestEvent,
                               REAL8 *priorStartFreq,
                               REAL8 *priorBand,
                               CHAR  *priorResultsFile);

void ComputeConfidence(LALStatus *status, REAL4 priorLoudestEvent, REAL4 maxPower, REAL8 *Confidence, REAL8 *conf_err);

void getStackSlideBinarySearchResults(LALStatus *status, StackSlideSearchParams *params, REAL8 *SearchLoudestEvent);

/*void ValidateMCResults(LALStatus *status, const REAL4FrequencySeries *oneSUM, StackSlideSearchParams *params, REAL4 *SNR);*/
void ValidateMCResults(LALStatus *status,
			REAL4FrequencySeries **SUMData,
			REAL4FrequencySeries **STKData,
			TdotsAndDeltaTs *pTdotsAndDeltaTs,
			StackSlideSearchParams *searchParams,
			INT4 iSky,
			StackSlideParams *params);

/*void ComputeUpperLimit(LALStatus *status, REAL8 *interpolatedUL, const REAL8 *arrayULs, const REAL8 *arrayConfs, REAL8 desiredConf, INT4 startIndex, INT4 numVals);*/
void ComputeUpperLimit(LALStatus *status, const REAL8 *arrayULs, const REAL8 *arrayConfs, REAL8 desiredConf);

void Readh0File(LALStatus *status, FILE *fph0, INT4 *N, REAL8 **h0);

#ifdef __cplusplus
}
#endif

#endif /* _STACKSLIDEBINARY_H */
