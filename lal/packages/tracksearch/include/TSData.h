/**** <lalVerbatim file="TSDataHV"> ***********
Author: Torres, C
$Id$
***** </lalVerbatim> ***********************************/

#ifndef _TSDATA_H
#define _TSDATA_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/TSSearch.h>
#include <lal/LALRCSID.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID (TSDATAH, "$Id$");

/******** <lalErrTable file="TSDataHErrTab"> ********/
#define TSDATA_ENULL    1
#define TSDATA_ENNUL    2
#define TSDATA_EALOC    3
#define TSDATA_ESEGZ    4
#define TSDATA_ENUMZ    5

#define TSDATA_MSGENULL "Null pointer"
#define TSDATA_MSGENNUL "Non-null pointer"
#define TSDATA_MSGEALOC "Memory allocation error"   
#define TSDATA_MSGESEGZ "Invalid number of segments"
#define TSDATA_MSGENUMZ "Invalid number of points in segment"
/******** </lalErrTable> ********/

typedef struct
tagTSCreateParams
{
  UINT4      dataSegmentPoints;     /* Fixed Length Varies Values  */
  UINT4      responseSegmentPoints; /* Fixed Copy for each segment */
  UINT4      spectraSegmentPoints;  /* Fixed Copy for each segment */
  UINT4      numberDataSegments;    /* Number of data segments to analyze */
}TSCreateParams;

void
LALCreateTSDataSegmentVector (
    LALStatus                  *status,
    TSSegmentVector           **vector,
    TSCreateParams             *params
    );

void
LALDestroyTSDataSegmentVector (
    LALStatus                  *status,
    TSSegmentVector           **vector
    );

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif


