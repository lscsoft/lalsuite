/************************************ <lalVerbatim file="StackSlideIsolatedHV">
Author:  Mendell, G.
$Id$
************************************* </lalVerbatim> */

/* REVISIONS: */
/* 04/12/05 gam; Move from using StackSlideOld to using StackSlide function. */
/* 04/12/05 gam; add RunStackSlideIsolatedMonteCarloSimulation to StackSlideIsolated.c */
/* 05/24/05 gam; make maxPower and totalEventCount part of params; change finishSUMs to finishPeriodicTable; end xml in FinalizeSearch */
/* 07/15/2005 gam; Change RunStackSlideIsolatedMonteCarloSimulation to: */
/*                 inject at a random point in the parameters space a random number of times */
/* 07/29/05 gam; if (params->testFlag & 64) > 0 set searchSurroundingPts == 1 and   */
/*               search surrounding parameters space pts; else search nearest only  */
/* 07/29/05 gam; rescale SFTs for inject given number of times for repeat of search */

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
  
NRCSID (STACKSLIDEISOLATEDH, "$Id$");

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
#define STACKSLIDEISOLATEDH_MSGENULL "Null Pointer"
#define STACKSLIDEISOLATEDH_MSGENNUL "Non-Null Pointer"
#define STACKSLIDEISOLATEDH_MSGENEGA "Bad Negative Value"
#define STACKSLIDEISOLATEDH_MSGEBADRESULTSFILE "Could not open priorResultsFile"
#define STACKSLIDEISOLATEDH_MSGEMCEVENTTHRESHOLDFLAGS "Monte Carlo needs (outputEventFlag and 2) > 0 and thresholdFlag <= 0"
#define STACKSLIDEISOLATEDH_MSGEIFO "Invalid or null IFO"
#define STACKSLIDEISOLATEDH_MSGENOFREQBINS "Avoiding or cleaning all freq bins! Monte Carlo cannot inject any signals!"
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

#ifdef INCLUDE_RUNSTACKSLIDEISOLATEDMONTECARLO_CODEOLD
/* 04/12/05 gam */
void RunStackSlideIsolatedMonteCarloSimulationOld(LALStatus *status, StackSlideSearchParams *params, INT4 nSamples);
#endif

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

#ifdef __cplusplus
}
#endif

#endif /* _STACKSLIDEISOLATED_H */
