/**** <lalVerbatim file="EPDataHV"> ***********
Author: Brady, P
$Id$
***** </lalVerbatim> ***********************************/

#ifndef _EPDATA_H
#define _EPDATA_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/EPSearch.h>
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


