/**** <lalVerbatim file="PowerApplySearchHV"> ***********
Author: Brady, P and Brown, D
 $Id$
***** </lalVerbatim> ***********************************/

#ifndef _APPLYSEARCH_H
#define _APPLYSEARCH_H

#include <lal/EPData.h>
#include "Power.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (APPLYSEARCHH, "power $Id$");

/******** <lalErrTable file="PowerApplySearchHErrTab"> ********/
#define APPLYSEARCHH_ENULL 1
#define APPLYSEARCHH_ENNUL 2
#define APPLYSEARCHH_EALOC 3
#define APPLYSEARCHH_ENUMZ 4
#define APPLYSEARCHH_EDELT 8
#define APPLYSEARCHH_ERANK 9
#define APPLYSEARCHH_EUEXT 10
#define APPLYSEARCHH_EDELF 11

#define APPLYSEARCHH_MSGENULL "Null pointer"
#define APPLYSEARCHH_MSGENNUL "Non-null pointer"
#define APPLYSEARCHH_MSGEALOC "Memory allocation error"
#define APPLYSEARCHH_MSGENUMZ "Data segment length is zero"
#define APPLYSEARCHH_MSGEDELT "deltaT is zero or negative"
#define APPLYSEARCHH_MSGERANK "Search node has incorrect rank"
#define APPLYSEARCHH_MSGEUEXT "Unrecognised exchange type"
#define APPLYSEARCHH_MSGEDELF "Inconsistent deltaF in spectrum and data"
/******** </lalErrTable> ********/

#define APPLYSEARCHH_SEARCHNAME "power"
#define APPLYSEARCHH_STRLENGTH 2

typedef struct
tagPowerOutPutParams
{
  INT4                          numEvents;
  INT8                          tstart;
  REAL8                         flow;
}
PowerOutPutParams;

void
BuildPowerOutPut (
    LALStatus          *status,
    LALSearchOutput    *output,
    BurstEventList     *burstEventList,
    PowerOutPutParams  *params
                    );

#ifdef  __cplusplus
}
#endif

#endif /* _APPLYSEARCH_H */
