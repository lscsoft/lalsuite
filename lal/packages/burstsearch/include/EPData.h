/**** <lalVerbatim file="EPDataHV"> ***********
Author: Brady, P
$Id$
***** </lalVerbatim> ***********************************/

#ifndef _EPDATA_H
#define _EPDATA_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/ExcessPower.h>
#include <lal/LALRCSID.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID (EPDATAH, "$Id$");

/******** <lalErrTable file="EPDataHErrTab"> ********/
#define EPDATA_ENULL    1
#define EPDATA_ENNUL    2
#define EPDATA_EALOC    3
#define EPDATA_ESEGZ    4
#define EPDATA_ENUMZ    5

#define EPDATA_MSGENULL "Null pointer"
#define EPDATA_MSGENNUL "Non-null pointer"
#define EPDATA_MSGEALOC "Memory allocation error"   
#define EPDATA_MSGESEGZ "Invalid number of segments"
#define EPDATA_MSGENUMZ "Invalid number of points in segment"
/******** </lalErrTable> ********/


typedef struct
tagEPDataSegment
{
  REAL4TimeSeries              *data;
  REAL4FrequencySeries         *spec;
  COMPLEX8FrequencySeries      *resp;
  INT4                          endOfData;
  INT4                          newLock;
  INT4                          newCal;
  INT4                          number;
}
EPDataSegment;


typedef struct
tagEPDataSegmentVector
{
  UINT4                         length;
  EPDataSegment                  *data;
}
EPDataSegmentVector;


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
  EPDataSegmentVector          *epSegVec;
  CreateTFTilingIn             *tfTilingInput;
  TFTiling                     *tfTiling;
  ComputeExcessPowerIn         *compEPInput;
}
EPSearchParams;



void
LALCreateEPDataSegmentVector (
    LALStatus                  *status,
    EPDataSegmentVector       **vector,
    EPInitParams               *params
    );

void
LALDestroyEPDataSegmentVector (
    LALStatus                  *status,
    EPDataSegmentVector       **vector
    );

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif


