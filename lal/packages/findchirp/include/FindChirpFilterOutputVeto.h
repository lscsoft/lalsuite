/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpFilterOutputVeto.h
 *
 * Author: 
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _FINDCHIRPFILTEROUTPUTVETOH_H
#define _FINDCHIRPFILTEROUTPUTVETOH_H

#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/FindChirpDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPFILTEROUTPUTVETOH, "$Id$");


#define FINDCHIRPFILTEROUTPUTVETOH_ENULL 1
#define FINDCHIRPFILTEROUTPUTVETOH_ENNUL 2
#define FINDCHIRPFILTEROUTPUTVETOH_MSGENULL "Null pointer"
#define FINDCHIRPFILTEROUTPUTVETOH_MSGENNUL "Non-null pointer"

typedef struct
tagFindChirpFilterOutputVetoParams
{
  UINT4         window;
  UINT4         length;
  REAL4         cutoff;
}
FindChirpFilterOutputVetoParams;

void LALFindChirpFilterOutputVeto( 
    LALStatus                          *status,
    SnglInspiralTable                 **eventList, 
    COMPLEX8Vector                     *qVec,
    REAL4                               qNorm,
    FindChirpFilterOutputVetoParams    *params
    );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPFILTEROUTPUTVETO_H */

