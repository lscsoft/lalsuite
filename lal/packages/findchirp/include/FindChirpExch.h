/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpExch.h
 *
 * Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _FINDCHIRPEXCH_H
#define _FINDCHIRPEXCH_H

#include <lal/LALDatatypes.h>
#include <lal/Comm.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (FINDCHIRPEXCHH, "$Id$");

#define FINDCHIRPEXCH_ENULL 1
#define FINDCHIRPEXCH_ENNUL 2
#define FINDCHIRPEXCH_ENOBJ 4
#define FINDCHIRPEXCH_EHAND 8

#define FINDCHIRPEXCH_MSGENULL "Null pointer"
#define FINDCHIRPEXCH_MSGENNUL "Non-null pointer"
#define FINDCHIRPEXCH_MSGENOBJ "Invalid number of objects"
#define FINDCHIRPEXCH_MSGEHAND "Wrong handshake"


void
LALExchangeDataSegment (
    LALStatus      *status,
    DataSegment *segment,
    ExchParams  *exchParams
    );

void
LALExchangeInspiralBankIn (
    LALStatus         *status,
    InspiralBankIn *bankIn,
    ExchParams     *exchParams
    );

void
LALExchangeInspiralTemplate (
    LALStatus           *status,
    InspiralTemplate *tmplt,
    ExchParams       *exchParams
    );

void
LALExchangeInspiralEvent (
    LALStatus        *status,
    InspiralEvent *event,
    ExchParams    *exchParams
    );

void
LALExchangeInspiralEventList (
    LALStatus     *status,
    InspiralEvent **eventHead,
    ExchParams    *exchParams
    );

void
LALExchangeTemplateBank (
    LALStatus         *status,
    InspiralTemplate **tmpltHead,
    ExchParams        *exchParms
                 );


#ifdef  __cplusplus
}
#endif

#endif /* _FINDCHIRPEXCH_H */
