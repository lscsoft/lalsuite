/************************************ <lalVerbatim file="ComputeSkyHV">
Author:  Virginia Re
$Id$
************************************* </lalVerbatim> */

/* Revisions: */
/* 04/12/05 gam; Add StackSlideSearchParams *params to StackSlideBinary. Need to include "DriveStackSlide.h" */
/* 04/12/05 gam; Remove from StackSlideParams *stksldParams, those already in StackSlideSearchParams *params */

#ifndef _STACKSLIDEBINARY_H
#define _STACKSLIDEBINARY_H

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
  
NRCSID (STACKSLIDEBINARYH, "$Id$");

/* 04/12/05 gam; add StackSlideSearchParams *params */
void StackSlideBinary(  LALStatus *status,
                        StackSlideParams *stksldParams,
                        StackSlideSearchParams *params
                      );

/* 04/12/05 gam */
void FindBinaryLoudest(REAL8 *LoudestEvent, REAL8 *peakFreq, REAL4FrequencySeries **SUMData, StackSlideParams *stksldParams);
#ifdef __cplusplus
}
#endif

#endif /* _STACKSLIDEBINARY_H */
