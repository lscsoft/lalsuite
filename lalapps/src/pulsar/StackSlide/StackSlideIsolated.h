/************************************ <lalVerbatim file="StackSlideIsolatedHV">
Author:  Mendell, G.
$Id$
************************************* </lalVerbatim> */

/* REVISIONS: */

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
#define STACKSLIDEISOLATEDH_ENEGA 4
#define STACKSLIDEISOLATEDH_MSGENULL "Null Pointer"
#define STACKSLIDEISOLATEDH_MSGENNUL "Non-Null Pointer"
#define STACKSLIDEISOLATEDH_MSGENEGA "Bad Negative Value"
/*********************************************/
/*                                           */
/* END SECTION: define constants             */
/*                                           */
/*********************************************/

void StackSlideIsolated (
    LALStatus                        *status,
    REAL4                            *maxPower,
    INT4                             *totalEventCount,
    SnglStackSlidePeriodicTable      *loudestPeaksArray,
    LALFindStackSlidePeakOutputs     *pLALFindStackSlidePeakOutputs,
    LALFindStackSlidePeakParams      *pLALFindStackSlidePeakParams,
    LALUpdateLoudestStackSlideParams *pLALUpdateLoudestStackSlideParams,
    LALDetector                      *cachedDetector,
    StackSlideParams                 *stksldParams,
    StackSlideSearchParams           *params
);

void StackSlideOld(	LALStatus *status, 
			REAL4FrequencySeries **SUMData, 
			REAL4FrequencySeries **STKData, 
			StackSlideParams *params);

#ifdef __cplusplus
}
#endif

#endif /* _STACKSLIDEISOLATED_H */
