/**** <lalVerbatim file="TSDatgenHV"> *********
Author: Torres. C
$ID: tracksearch.h,v 1.0 2004/04/14 02:00:00 charlie Exp $
***** </lalVerbatim> **********************************/

#ifndef TRACKSEARCHAVERAGER_H
#define TRACKSEARCHAVERAGER_H

#include "tracksearch.h"
#include "tracksearchToolbox.h"

/*
 * Error Codes:
 */

/* Non-Error Code Define Statements */
#define TRACKSEARCHAVERAGERC_NARGS 2
/* End N-E Defines */
/* Define Error Codes */
#define TRACKSEARCHAVERAGERC_ENORM              0
#define TRACKSEARCHAVERAGERC_ESUB               1
#define TRACKSEARCHAVERAGERC_EARGS              2
#define TRACKSEARCHAVERAGERC_EVAL               4
#define TRACKSEARCHAVERAGERC_EFILE              8
#define TRACKSEARCHAVERAGERC_EREAD              16
#define TRACKSEARCHAVERAGERC_EMEM               32
#define TRACKSEARCHAVERAGERC_EMISC              64
#define TRACKSEARCHAVERAGERC_EWRITE             128
#define TRACKSEARCHAVERAGERC_EDIMS              256 
#define TRACKSEARCHAVERAGERC_ENULL              512 
#define TRACKSEARCHAVERAGERC_ENPTR              1024
#define TRACKSEARCHAVERAGERC_EMAPO              2048

#define TRACKSEARCHAVERAGERC_MSGENORM          "Normal Exit"
#define TRACKSEARCHAVERAGERC_MSGESUB           "Subroutine Fail"
#define TRACKSEARCHAVERAGERC_MSGEARGS          "Arguement Parse Error"
#define TRACKSEARCHAVERAGERC_MSGEVAL           "Invalid Argument(s)"
#define TRACKSEARCHAVERAGERC_MSGEFILE          "Error Opening File"
#define TRACKSEARCHAVERAGERC_MSGEREAD          "Error Reading File"
#define TRACKSEARCHAVERAGERC_MSGEMEM           "Memory Error"
#define TRACKSEARCHAVERAGERC_MSGEMISC          "Unknown Error"
#define TRACKSEARCHAVERAGERC_MSGEWRITE         "Error Writing File"
#define TRACKSEARCHAVERAGERC_MSGEDIMS          "Map Scale Error"
#define TRACKSEARCHAVERAGERC_MSGENULL          "NULL pointer expected"
#define TRACKSEARCHAVERAGERC_MSGENPTR          "Unexpected NULL pointer"
#define TRACKSEARCHAVERAGERC_MSGEMAPO          "Input maps do not intersect properly"

#define TRUE     1
#define FALSE    0

/*
 * Averager data structures and functions 
 */
#define NumberAveragerOperations 2
typedef enum averagerOperationType {
  Unknown,Merge,Collapse
} averagerOperationType;

/*
 * Averager variable prototypes
 */
typedef struct tagTSAMap {
  CreateTimeFreqIn              imageCreateParams;
  UINT4                         clippedWith; /*1=another map*/
  LIGOTimeGPS                   clipperMapStart;
  TrackSearchMapMarkingParams   imageBorders;
  TimeFreqRep                  *imageRep;
} TSAMap;

typedef struct tagTSAcache {
  UINT4           numMapFilenames;
  CHARVector    **filename;
  REAL8          *mapStartTime;
}
TSAcache;

typedef struct tagTSACollapseParams {
  /* 
   * Should be UINT4 but timefreq package uses INT4
   * So we do this to avoid compiler warnings
   */
  INT4                  averageTBins;
  INT4                  averageFBins;
} TSACollapseParams;

typedef struct tagTSAparams {
  averagerOperationType  operation;
  CHARVector            *cacheFilename;
  CHARVector            *multiCacheFilename;
  TSACollapseParams      colParams;
  TSDiagnosticType       verbosity;
} TSAparams;

/*
 * Subroutine Prototypes
 */
void
LALappsTSACollapseMap(LALStatus*,
		      TSAMap**,
		      TSACollapseParams);

void
LALappsTSAMergeMap(LALStatus  *status,
		   TSAMap  **output,
		   TSAMap    inputA,
		   TSAMap    inputB);
void
LALappsTSAInitialize(LALStatus     *status,
		     int            argc,
		     char*          argv[],
		     TSAparams     *params);

void
LALappsTSAProcessSingleCache(LALStatus    *status,
			     TSAparams    *params);

void
tsaTest(LALStatus  *status);

#endif
