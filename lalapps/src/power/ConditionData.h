/**** <lalVerbatim file="PowerConditionDataHV"> ***********
Author: Brady, P and Brown, D
$Id$
***** </lalVerbatim> ***********************************/

#ifndef _CONDITIONDATA_H
#define _CONDITIONDATA_H

#include "LALWrapperInterface.h"
#include "Power.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (CONDITIONDATAH, "power $Id$");

/******** <lalErrTable file="PowerConditionDataHErrTab"> ********/
#define CONDITIONDATAH_ENULL     1
#define CONDITIONDATAH_ENNUL     2
#define CONDITIONDATAH_EALOC     3
#define CONDITIONDATAH_EARGS     4
#define CONDITIONDATAH_EINPUT    5
#define CONDITIONDATAH_EDT       6
#define CONDITIONDATAH_EDTZ      7
#define CONDITIONDATAH_ENUMZ     8
#define CONDITIONDATAH_ESEGZ     9
#define CONDITIONDATAH_EDATZ     10
#define CONDITIONDATAH_ESPEC     11

#define CONDITIONDATAH_MSGENULL   "Null pointer"
#define CONDITIONDATAH_MSGENNUL   "Non-null pointer"
#define CONDITIONDATAH_MSGEALOC   "Memory allocation error"
#define CONDITIONDATAH_MSGEARGS   "Wrong number of arguments"
#define CONDITIONDATAH_MSGEINPUT  "Wrong data names in InPut structure"
#define CONDITIONDATAH_MSGEDT     "Mismatch between deltaT in wrapperParams and inPut"
#define CONDITIONDATAH_MSGEDTZ    "wrapperParams->deltaT is zero or negative"
#define CONDITIONDATAH_MSGENUMZ   "Data segment length is zero"
#define CONDITIONDATAH_MSGESEGZ   "Number of data segments is zero"
#define CONDITIONDATAH_MSGEDATZ   "Got less data than expected"
#define CONDITIONDATAH_MSGESPEC   "Unrecognized power spectrum type"
/******** </lalErrTable> ********/

/* the names of the multidim data channels */
#define INPUTNAME_CHANNEL   "channel"
#define INPUTNAME_SPECTRUM  "spectrum"
#define INPUTNAME_RESPONSE  "response"

#ifdef  __cplusplus
}
#endif

#endif /* _CONDITIONDATA_H */
