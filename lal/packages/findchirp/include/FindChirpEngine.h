/*----------------------------------------------------------------------- 
 * 
 * File Name: ApplySearch.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */


#ifndef _FINDCHIRPENGINEH_H
#define _FINDCHIRPENGINEH_H

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <lal/LALStdlib.h>
#include <lal/DataBuffer.h>
#include <lal/Comm.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpExch.h>
#include <lal/FindChirpSP.h>


#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (FINDCHIRPENGINEHH, "$Id$");

#define FINDCHIRPENGINEH_ENULL 1
#define FINDCHIRPENGINEH_ENNUL 2
#define FINDCHIRPENGINEH_ENUMZ 3
#define FINDCHIRPENGINEH_EDELT 4
#define FINDCHIRPENGINEH_ERHOZ 5
#define FINDCHIRPENGINEH_ECHIZ 6
#define FINDCHIRPENGINEH_ETMPL 7
#define FINDCHIRPENGINEH_EALOC 8
#define FINDCHIRPENGINEH_ERANK 9
#define FINDCHIRPENGINEH_EUEXT 10


#define FINDCHIRPENGINEH_MSGENULL "Null pointer"
#define FINDCHIRPENGINEH_MSGENNUL "Non-null pointer"
#define FINDCHIRPENGINEH_MSGENUMZ "Data segment length is zero"
#define FINDCHIRPENGINEH_MSGEDELT "deltaT is zero or negative"
#define FINDCHIRPENGINEH_MSGERHOZ "snr squared threshold is zero or negative"
#define FINDCHIRPENGINEH_MSGECHIZ "chi squared threshold is zero or negative"
#define FINDCHIRPENGINEH_MSGETMPL "linked list of templates to filter in null"
#define FINDCHIRPENGINEH_MSGEALOC "Memory allocation error"
#define FINDCHIRPENGINEH_MSGERANK "Search node has incorrect rank"
#define FINDCHIRPENGINEH_MSGEUEXT "Unrecognised exchange type"


enum ExchObjectType
{
  ExchDataSegment,
  ExchFindChirpSegment,
  ExchInspiralTemplate,
  ExchInspiralEvent,
  ExchFinished
};

typedef struct
tagFindChirpMasterParams
{
  UINT4                         numCoarseExch;
  MPI_Comm                     *mpiComm;
  UINT4                        *numSlaves;
  InspiralTemplate             *currentTmplt;
  UINT4                         numTmplts;
  BOOLEAN                      *notFinished;
  REAL4                        *fracRemaining;
}
FindChirpMasterParams;

typedef struct
tagFindChirpSlaveParams
{
  MPI_Comm                     *mpiComm;
  FindChirpSPDataParams        *dataParams;
  FindChirpSPTmpltParams       *tmpltParams;
  FindChirpFilterParams        *filterParams;
  FindChirpFilterInput         *filterInput;
}
FindChirpSlaveParams;

typedef struct
tagInspiralTemplateNode
{
  UINT4                                 level;
  struct tagInspiralTemplateNode       *next;
  struct tagInspiralTemplateNode       *prev;
  InspiralTemplate                     *tmpltPtr;
}
InspiralTemplateNode;


void
LALFindChirpMaster (
    LALStatus                  *status, 
    InspiralEvent             **eventList,
    FindChirpMasterParams       *params 
    );

void
LALFindChirpSlave (
    LALStatus                  *status, 
    BOOLEAN                    *notFinished,
    FindChirpSegmentVector     *fcSegVec,
    FindChirpSlaveParams        *params 
    );

void
LALFindChirpCreateInspiralBank (
    LALStatus                  *status,
    InspiralCoarseBankIn       *bankIn,
    InspiralTemplate          **head
    );

void
LALFindChirpDestroyInspiralBank (
    LALStatus           *status,
    InspiralTemplate   **head
    );

#ifdef  __cplusplus
}
#endif

#endif /* _FINDCHIRPENGINEH_H */
