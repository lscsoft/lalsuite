/**** <lalVerbatim file="EPExchHV"> ***********
Author: Allen, B., Brady, P., Brown, D and Creighton, J. D. E.
$Id$
***** </lalVerbatim> ***********************************/

#ifndef _EPEXCH_H
#define _EPEXCH_H

#include <lal/LALDatatypes.h>
#include <lal/Comm.h>
#include <lal/DataBuffer.h>
#include <lal/FindChirp.h>
#include <lal/ExcessPower.h>
#include <lal/BurstSearch.h>
#include "EPData.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (EPEXCHH, "$Id$");

/******** <lalErrTable file="EPExchHErrTab"> ********/
#define EPEXCH_ENULL 1
#define EPEXCH_ENNUL 2
#define EPEXCH_ENOBJ 4
#define EPEXCH_EHAND 8
#define EPEXCH_EMPIE 16

#define EPEXCH_MSGENULL "Null pointer"
#define EPEXCH_MSGENNUL "Non-null pointer"
#define EPEXCH_MSGENOBJ "Invalid number of objects"
#define EPEXCH_MSGEHAND "Wrong handshake"
#define EPEXCH_MSGEMPIE "Problem exchanging event list"
/******** </lalErrTable> ********/

enum
{
  ExchDataSegment,
  ExchDataList,
  ExchEPEvent,
  ExchFinished
}
ExchObjectType;

void
LALExchangeEPEvent (
    LALStatus        *status,
    TFTile           *event,
    ExchParams       *exchParams
    );

void
LALExchangeEPEventList (
    LALStatus     *status,
    BurstEvent   **eventHead,
    INT4           numEvents,
    ExchParams    *exchParams
    );

#ifdef  __cplusplus
}
#endif

#endif /* _EPEXCH_H */
