/**** <lalVerbatim file="PowerFinalizeSearchHV"> ***********
Author:  Brown, D
$Id$ 
***** </lalVerbatim> ***********************************/

#ifndef _FINALIZESEARCH_H
#define _FINALIZESEARCH_H

#include "LALWrapperInterface.h"
#include "Power.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (FINALIZESEARCHH, "power $Id$");

/******** <lalErrTable file="PowerFinalizeSearchHErrTab"> ********/
#define FINALIZESEARCHH_ENULL 1
#define FINALIZESEARCHH_ENNUL 2
#define FINALIZESEARCHH_EALOC 4
#define FINALIZESEARCHH_EARGS 8
#define FINALIZESEARCHH_ENUMZ 16

#define FINALIZESEARCHH_MSGENULL "Null pointer"
#define FINALIZESEARCHH_MSGENNUL "Non-null pointer"
#define FINALIZESEARCHH_MSGEALOC "Memory allocation error"
#define FINALIZESEARCHH_MSGEARGS "Wrong number of arguments"
#define FINALIZESEARCHH_MSGENUMZ "Number of data segments is zero"
/******** </lalErrTable> ********/


#ifdef  __cplusplus
}
#endif

#endif /* _FINALIZESEARCH_H */
