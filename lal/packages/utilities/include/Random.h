/*----------------------------------------------------------------------- 
 * 
 * File Name: Random.h
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _RANDOM_H
#define _RANDOM_H

#include "LALDatatypes.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (RANDOMH, "$Id$");

#define RANDOM_ENULL 1
#define RANDOM_ENNUL 2
#define RANDOM_ESIZE 4

#define RANDOM_MSGENULL "Null pointer"
#define RANDOM_MSGENNUL "Non-null pointer"
#define RANDOM_MSGESIZE "Invalid size"

typedef struct
tagRandomParams
{
  INT4 i;
  INT4 y;
  INT4 v[32];
}
RandomParams;

void
CreateRandomParams (
    Status        *status,
    RandomParams **params,
    INT4           seed
    );

void
DestroyRandomParams (
    Status        *status,
    RandomParams **params
    );

void
UniformDeviate (
    Status       *status,
    REAL4        *deviate,
    RandomParams *params
    );

void
NormalDeviates (
    Status       *status,
    REAL4Vector  *deviates,
    RandomParams *params
    );


#ifdef  __cplusplus
}
#endif

#endif /* _RANDOM_H */
