/********************************** <lalVerbatim file="ExcessPowerHV">
Author: Flanagan, E
$Id$
**************************************************** </lalVerbatim> */

#ifndef _EPSEARCH_H
#define _EPSEARCH_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/TFTransform.h>
#include <lal/ExcessPower.h>
#include <lal/BurstSearch.h>
#include <lal/LALRCSID.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID (EPSEARCHH, "$Id$");

/******** <lalErrTable file="ExcessPowerHErrTab"> ********/
#define EPSEARCHH_ENULLP       1
#define EPSEARCHH_EPOSARG      2
#define EPSEARCHH_EPOW2        4
#define EPSEARCHH_EMALLOC      8
#define EPSEARCHH_EINCOMP      16
#define EPSEARCHH_EORDER       32
#define EPSEARCHH_ENONNULL     64
#define EPSEARCHH_ETILES       65
#define EPSEARCHH_EDELF        128


#define EPSEARCHH_MSGENULLP    "Null pointer"
#define EPSEARCHH_MSGEPOSARG   "Arguments must be non-negative"
#define EPSEARCHH_MSGEPOW2     "Length of supplied data must be a power of 2"
#define EPSEARCHH_MSGEMALLOC   "Malloc failure"
#define EPSEARCHH_MSGEINCOMP   "Incompatible arguments"
#define EPSEARCHH_MSGEORDER    "Routines called in illegal order"
#define EPSEARCHH_MSGENONNULL  "Null pointer expected"
#define EPSEARCHH_MSGETILES    "Malloc failed while assigning memory for a tile"
#define EPSEARCHH_MSGEDELF     "Inconsistent deltaF in spectrum and data"
/******** </lalErrTable> ********/

typedef struct
tagEPInitParams
{
  UINT4                         numPoints;
  UINT4                         numSegments;
  UINT4                         segDutyCycle;
}
EPInitParams;

typedef struct
tagEPSearchParams
{
  BOOLEAN                       searchMaster;
  BOOLEAN                       haveData;
  UINT4                        *numSlaves;          
  UINT4                         simType;
  UINT4                         currentSegment;
  INT4                          ntotT;
  INT4                          numEvents;
  UINT4                         ovrlap;
  REAL8                         lambda;
  REAL8                         alphaThreshold;
  INT4                          events2Master;     /* Added Erik Kats */
  CHAR                         *channelName;       /* Added Erik Kats */
  EPInitParams                 *initParams;
  struct tagEPDataSegmentVector          *epSegVec;
  CreateTFTilingIn             *tfTilingInput;
  TFTiling                     *tfTiling;
  ComputeExcessPowerIn         *compEPInput;
  LALWindowParams               winParams;
}
EPSearchParams;

void
EPSearch (
		   LALStatus               *status,
                   EPSearchParams          *params,
                   BurstEvent             **burstEvent,
                   UINT4                    tmpDutyCyle
		   );



#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif




