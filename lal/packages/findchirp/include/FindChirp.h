/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirp.h
 *
 * Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _FINDCHIRP_H
#define _FINDCHIRP_H

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/DataBuffer.h>
#include <lal/Comm.h>
#include <lal/Inspiral.h>
#include <lal/FindChirpChisq.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (FINDCHIRPH, "$Id$");

#define FINDCHIRP_ENULL 1
#define FINDCHIRP_ENNUL 2
#define FINDCHIRP_EALOC 3
#define FINDCHIRP_ENUMZ 5
#define FINDCHIRP_ESEGZ 6
#define FINDCHIRP_ECHIZ 7
#define FINDCHIRP_EDTZO 8
#define FINDCHIRP_ETRNC 10
#define FINDCHIRP_EFLOW 11
#define FINDCHIRP_EFREE 12
#define FINDCHIRP_ERHOT 15
#define FINDCHIRP_ECHIT 16


#define FINDCHIRP_MSGENULL "Null pointer"
#define FINDCHIRP_MSGENNUL "Non-null pointer"
#define FINDCHIRP_MSGEALOC "Memory allocation error"
#define FINDCHIRP_MSGENUMZ "Invalid number of points in segment"
#define FINDCHIRP_MSGESEGZ "Invalid number of segments"
#define FINDCHIRP_MSGECHIZ "Invalid number of chi squared bins"
#define FINDCHIRP_MSGEDTZO "deltaT is zero or negative"
#define FINDCHIRP_MSGETRNC "Duration of inverse spectrum in time domain is negative"
#define FINDCHIRP_MSGEFLOW "Inverse spectrum low frequency cutoff is negative"
#define FINDCHIRP_MSGEFREE "Memory free error"
#define FINDCHIRP_MSGERHOT "Rhosq threshold is zero or negative"
#define FINDCHIRP_MSGECHIT "Chisq threshold is zero or negative"



/*
 *
 * typedefs of structures used by the findchip functions
 *
 */


/* the various interferometer codes */
typedef enum
{
  Caltech40m, 
  Hanford4km, 
  Hanford2km, 
  Livingston4km, 
  GEO600m, 
  TAMA300m
}
Detector;

/* input for specifying a template bank */
typedef struct
tagInspiralBankIn
{
  REAL4                         mMin;
  REAL4                         mMax;
  REAL4                         ffCoarse;
  REAL4                         ffFine;
  Detector                      detector;
  Method                        method;
}
InspiralBankIn;

/* structure for describing a binary insipral event */
typedef struct
tagInspiralEvent
{
  UINT4                         id;
  LIGOTimeGPS                   time;
  UINT4                         timeIndex;
  InspiralTemplate              tmplt;
  REAL4                         snrsq;
  REAL4                         chisq;
  REAL4                         sigma;
  REAL4                         effDist;
  struct tagInspiralEvent      *next;
}
InspiralEvent;

/* vector of DataSegment, as defined the framedata package */
typedef struct
tagDataSegmentVector
{
  UINT4                         length;
  DataSegment                  *data;
}
DataSegmentVector;

/* processed data segment used by FindChirp filter routine */
typedef struct
tagFindChirpSegment
{
  COMPLEX8FrequencySeries      *data;
  UINT4Vector                  *chisqBinVec;
  REAL8                         deltaT;
  REAL4                         segNorm;
  UINT4                         number;
}
FindChirpSegment;

/* vector of FindChirpSegment defined above */
typedef struct
tagFindChirpSegmentVector
{
  UINT4                         length;
  FindChirpSegment             *data;
}
FindChirpSegmentVector;

/* structure to contain an inspiral template */
typedef struct
tagFindChirpTemplate
{
  COMPLEX8Vector               *data;
  REAL4                         tmpltNorm;
}
FindChirpTemplate;


/*
 *
 * typedefs of parameter structures used by functions in findchirp
 *
 */


/* parameter structure for all init funtions */
typedef struct
tagFindChirpInitParams
{
  UINT4                         numSegments;
  UINT4                         numPoints;
  UINT4                         numChisqBins;
  BOOLEAN                       createRhosqVec;
}
FindChirpInitParams;

/* parameter structure for the filtering function */
typedef struct
tagFindChirpFilterParams
{
  REAL4                         deltaT;
  REAL4                         rhosqThresh;
  REAL4                         chisqThresh;
  BOOLEAN                       computeNegFreq;
  COMPLEX8Vector               *qVec;
  COMPLEX8Vector               *qtildeVec;
  ComplexFFTPlan               *invPlan;
  REAL4Vector                  *rhosqVec;
  REAL4Vector                  *chisqVec;
  FindChirpChisqParams         *chisqParams;
  FindChirpChisqInput          *chisqInput;
}
FindChirpFilterParams;


/*
 *
 * typedefs of input structures used by functions in findchirp
 *
 */


/* input to the filtering functions */
typedef struct
tagFindChirpFilterInput
{
  InspiralTemplate             *tmplt;
  FindChirpTemplate            *fcTmplt;
  FindChirpSegment             *segment;
}
FindChirpFilterInput;


/*
 *
 * function prototypes for memory management functions
 *
 */


void
LALCreateDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector,
    FindChirpInitParams        *params
    );

void
LALDestroyDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector
    );

void
LALCreateFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector,
    FindChirpInitParams        *params
    );

void
LALDestroyFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector
    );


/*
 *
 * function prototypes for initialization and finalization functions
 *
 */


void
LALFindChirpFilterInit (
    LALStatus                  *status,
    FindChirpFilterParams     **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpFilterFinalize (
    LALStatus                  *status,
    FindChirpFilterParams     **output
    );

void
LALCreateFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output,
    FindChirpInitParams        *params
    );

void
LALDestroyFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output
    );


/*
 *
 * function prototype for the filtering function
 *
 */


void
LALFindChirpFilterSegment (
    LALStatus                  *status,
    InspiralEvent             **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );


#ifdef  __cplusplus
}
#endif

#endif /* _FINDCHIRP_H */
