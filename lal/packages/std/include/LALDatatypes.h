/*-----------------------------------------------------------------------
 *
 * File Name: LALDatatypes.h
 *
 * Author: Finn, L. S.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * LALDatatypes.h
 *
 * SYNOPSIS
 * #include "LALDatatypes.h"
 *
 * DESCRIPTION
 * Defines structured data types.
 *
 * DIAGNOSTICS
 *
 *----------------------------------------------------------------------
 */

#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H

#include "LALAtomicDatatypes.h"
#include "LALRCSID.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALDATATYPESH, "$Id$");


/* Vectors */

typedef struct
tagCHARVector
{
  UINT4  length;
  CHAR  *data;
}
CHARVector;

typedef struct
tagINT2Vector
{
  UINT4  length;
  INT2  *data;
}
INT2Vector;

typedef struct
tagUINT2Vector
{
  UINT4  length;
  UINT2 *data;
}
UINT2Vector;

typedef struct
tagINT4Vector
{
  UINT4  length;
  INT4  *data;
}
INT4Vector;

typedef struct
tagUINT4Vector
{
  UINT4  length;
  UINT4  *data;
}
UINT4Vector;

typedef struct
tagINT8Vector
{
  UINT4  length;
  INT8  *data;
}
INT8Vector;

typedef struct
tagUINT8Vector
{
  UINT4  length;
  UINT8 *data;
}
UINT8Vector;

typedef struct
tagREAL4Vector
{
  UINT4  length;
  REAL4 *data;
}
REAL4Vector;

typedef struct tagREAL8Vector
{
  UINT4  length;
  REAL8 *data;
}
REAL8Vector;

typedef struct tagCOMPLEX8Vector
{
  UINT4     length;
  COMPLEX8 *data;
}
COMPLEX8Vector;

typedef struct tagCOMPLEX16Vector
{
  UINT4      length;
  COMPLEX16 *data;
}
COMPLEX16Vector;

/* Arrays */

typedef struct
tagINT2Array
{
  UINT4Vector *dimLength;
  INT2        *data;
}
INT2Array;

typedef struct
tagUINT2Array
{
  UINT4Vector *dimLength;
  UINT2       *data;
}
UINT2Array;

typedef struct
tagINT4Array
{
  UINT4Vector *dimLength;
  INT4        *data;
}
INT4Array;

typedef struct
tagUINT4Array
{
  UINT4Vector *dimLength;
  UINT4       *data;
}
UINT4Array;

typedef struct
tagINT8Array
{
  UINT4Vector *dimLength;
  INT8        *data;
}
INT8Array;

typedef struct
tagUINT8Array
{
  UINT4Vector *dimLength;
  UINT8       *data;
}
UINT8Array;

typedef struct
tagREAL4Array
{
  UINT4Vector *dimLength;
  REAL4       *data;
}
REAL4Array;

typedef struct
tagREAL8Array
{
  UINT4Vector *dimLength;
  REAL8       *data;
}
REAL8Array;

typedef struct
tagCOMPLEX8Array
{
  UINT4Vector *dimLength;
  COMPLEX8    *data;
}
COMPLEX8Array;

typedef struct
tagCOMPLEX16Array
{
  UINT4Vector *dimLength;
  COMPLEX16   *data;
}
COMPLEX16Array;


/* Sequences */

typedef CHARVector      CHARSequence;
typedef INT2Vector      INT2Sequence;
typedef UINT2Vector     UINT2Sequence;
typedef INT4Vector      INT4Sequence;
typedef UINT4Vector     UINT4Sequence;
typedef INT8Vector      INT8Sequence;
typedef UINT8Vector     UINT8Sequence;
typedef REAL4Vector     REAL4Sequence;
typedef REAL8Vector     REAL8Sequence;
typedef COMPLEX8Vector  COMPLEX8Sequence;
typedef COMPLEX16Vector COMPLEX16Sequence;

/* Vector Sequences */

typedef struct
tagCHARVectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  CHAR  *data;
}
CHARVectorSequence;

typedef struct
tagINT2VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  INT2  *data;
}
INT2VectorSequence;

typedef struct
tagUINT2VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  UINT2 *data;
}
UINT2VectorSequence;

typedef struct
tagINT4VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  INT4  *data;
}
INT4VectorSequence;

typedef struct
tagUINT4VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  UINT4 *data;
}
UINT4VectorSequence;

typedef struct
tagINT8VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  INT8  *data;
}
INT8VectorSequence;

typedef struct
tagUINT8VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  UINT8 *data;
}
UINT8VectorSequence;

typedef struct
tagREAL4VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  REAL4 *data;
}
REAL4VectorSequence;

typedef struct
tagREAL8VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  REAL8 *data;
}
REAL8VectorSequence;

typedef struct
tagCOMPLEX8VectorSequence
{
  UINT4     length;
  UINT4     vectorLength;
  COMPLEX8 *data;
}
COMPLEX8VectorSequence;

typedef struct
tagCOMPLEX16VectorSequence
{
  UINT4      length;
  UINT4      vectorLength;
  COMPLEX16 *data;
}
COMPLEX16VectorSequence;

/* 
 * Array Sequences 
 */

typedef struct
tagINT2ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  INT2        *data;
}
INT2ArraySequence;

typedef struct
tagUINT2ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  UINT2       *data;
}
UINT2ArraySequence;

typedef struct
tagINT4ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  INT4        *data;
}
INT4ArraySequence;

typedef struct
tagUINT4ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  UINT4       *data;
}
UINT4ArraySequence;

typedef struct
tagINT8ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  INT8        *data;
}
INT8ArraySequence;

typedef struct
tagUINT8ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  UINT8       *data;
}
UINT8ArraySequence;

typedef struct
tagREAL4ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  REAL4       *data;
}
REAL4ArraySequence;

typedef struct
tagREAL8ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  REAL8       *data;
}
REAL8ArraySequence;

typedef struct
tagCOMPLEX8ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  COMPLEX8    *data;
}
COMPLEX8ArraySequence;

typedef struct
tagCOMPLEX16ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  COMPLEX16   *data;
}
COMPLEX16ArraySequence;


/* 
 * STRUCTURED DATA TYPES
 */

/* Time object */ 

typedef struct
tagLIGOTimeGPS
{
  INT4 gpsSeconds;
  INT4 gpsNanoSeconds;
}
LIGOTimeGPS;

/* Time Series */

typedef struct
tagINT2TimeSeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  REAL8         f0;
  CHARVector   *sampleUnits;
  INT2Sequence *data;
}
INT2TimeSeries;

typedef struct
tagUINT2TimeSeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  CHARVector    *sampleUnits;
  UINT2Sequence *data;
}
UINT2TimeSeries;

typedef struct
tagINT4TimeSeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  REAL8         f0;
  CHARVector   *sampleUnits;
  INT4Sequence *data;
}
INT4TimeSeries;

typedef struct
tagUINT4TimeSeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  CHARVector    *sampleUnits;
  UINT4Sequence *data;
}
UINT4TimeSeries;

typedef struct
tagINT8TimeSeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  REAL8         f0;
  CHARVector   *sampleUnits;
  INT8Sequence *data;
}
INT8TimeSeries;

typedef struct
tagUINT8TimeSeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  CHARVector    *sampleUnits;
  UINT8Sequence *data;
}
UINT8TimeSeries;

typedef struct
tagREAL4TimeSeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  CHARVector    *sampleUnits;
  REAL4Sequence *data;
}
REAL4TimeSeries;

typedef struct
tagREAL8TimeSeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  CHARVector    *sampleUnits;
  REAL8Sequence *data;
}
REAL8TimeSeries;

typedef struct
tagCOMPLEX8TimeSeries
{
  CHAR             *name;
  LIGOTimeGPS       epoch;
  REAL8             deltaT;
  REAL8             f0;
  CHARVector       *sampleUnits;
  COMPLEX8Sequence *data;
}
COMPLEX8TimeSeries;

typedef struct
tagCOMPLEX16TimeSeries
{
  CHAR              *name;
  LIGOTimeGPS        epoch;
  REAL8              deltaT;
  REAL8              f0;
  CHARVector        *sampleUnits;
  COMPLEX16Sequence *data;
}
COMPLEX16TimeSeries;

/* Vector Time Series */

typedef struct
tagINT2TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  INT2VectorSequence  *data;
}
INT2TimeVectorSeries;

typedef struct
tagUINT2TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  UINT2VectorSequence *data;
}
UINT2TimeVectorSeries;

typedef struct
tagINT4TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  INT4VectorSequence  *data;
}
INT4TimeVectorSeries;

typedef struct
tagUINT4TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  UINT4VectorSequence *data;
}
UINT4TimeVectorSeries;

typedef struct
tagINT8TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  INT8VectorSequence  *data;
}
INT8TimeVectorSeries;

typedef struct
tagUINT8TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  UINT8VectorSequence *data;
}
UINT8TimeVectorSeries;

typedef struct
tagREAL4TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  REAL4VectorSequence *data;
}
REAL4TimeVectorSeries;

typedef struct
tagREAL8TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  REAL8VectorSequence *data;
}
REAL8TimeVectorSeries;

typedef struct
tagCOMPLEX8TimeVectorSeries
{
  CHAR                    *name;
  LIGOTimeGPS              epoch;
  REAL8                    deltaT;
  REAL8                    f0;
  CHARVector              *sampleUnits;
  COMPLEX8VectorSequence  *data;
}
COMPLEX8TimeVectorSeries;

typedef struct
tagCOMPLEX16TimeVectorSeries
{
  CHAR                     *name;
  LIGOTimeGPS               epoch;
  REAL8                     deltaT;
  REAL8                     f0;
  CHARVector               *sampleUnits;
  COMPLEX16VectorSequence  *data;
}
COMPLEX16TimeVectorSeries;

/* Array Time Series */

typedef struct
tagINT2TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  INT2ArraySequence  *data;
}
INT2TimeArraySeries;

typedef struct
tagUINT2TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  UINT2ArraySequence *data;
}
UINT2TimeArraySeries;

typedef struct
tagINT4TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  INT4ArraySequence  *data;
}
INT4TimeArraySeries;

typedef struct
tagUINT4TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  UINT4ArraySequence *data;
}
UINT4TimeArraySeries;

typedef struct
tagINT8TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  INT8ArraySequence  *data;
}
INT8TimeArraySeries;

typedef struct
tagUINT8TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  UINT8ArraySequence *data;
}
UINT8TimeArraySeries;

typedef struct
tagREAL4TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  REAL4ArraySequence *data;
}
REAL4TimeArraySeries;

typedef struct
tagREAL8TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  REAL8ArraySequence *data;
}
REAL8TimeArraySeries;

typedef struct
tagCOMPLEX8TimeArraySeries
{
  CHAR                  *name;
  LIGOTimeGPS            epoch;
  REAL8                  deltaT;
  REAL8                  f0;
  CHARVector            *sampleUnits;
  COMPLEX8ArraySequence *data;
}
COMPLEX8TimeArraySeries;

typedef struct
tagCOMPLEX16TimeArraySeries
{
  CHAR                   *name;
  LIGOTimeGPS             epoch;
  REAL8                   deltaT;
  REAL8                   f0;
  CHARVector             *sampleUnits;
  COMPLEX16ArraySequence *data;
}
COMPLEX16TimeArraySeries;

/* Frequency Series */

typedef struct
tagINT2FrequencySeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         f0;
  REAL8         deltaF;
  CHARVector   *sampleUnits;
  INT2Sequence *data;
}
INT2FrequencySeries;

typedef struct
tagUINT2FrequencySeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          f0;
  REAL8          deltaF;
  CHARVector    *sampleUnits;
  UINT2Sequence *data;
}
UINT2FrequencySeries;

typedef struct
tagINT4FrequencySeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         f0;	
  REAL8         deltaF;
  CHARVector   *sampleUnits;
  INT4Sequence *data;
}
INT4FrequencySeries;

typedef struct
tagUINT4FrequencySeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          f0;	
  REAL8          deltaF;
  CHARVector    *sampleUnits;
  UINT4Sequence *data;
}
UINT4FrequencySeries;

typedef struct
tagINT8FrequencySeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         f0;
  REAL8         deltaF;
  CHARVector   *sampleUnits;
  INT8Sequence *data;
}
INT8FrequencySeries;

typedef struct
tagUINT8FrequencySeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          f0;
  REAL8          deltaF;
  CHARVector    *sampleUnits;
  UINT8Sequence *data;
}
UINT8FrequencySeries;

typedef struct
tagREAL4FrequencySeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          f0;	
  REAL8          deltaF;
  CHARVector    *sampleUnits;
  REAL4Sequence *data;
}
REAL4FrequencySeries;

typedef struct
tagREAL8FrequencySeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          f0;
  REAL8          deltaF;
  CHARVector    *sampleUnits;
  REAL8Sequence *data;
}
REAL8FrequencySeries;

typedef struct
tagCOMPLEX8FrequencySeries
{
  CHAR             *name;
  LIGOTimeGPS       epoch;
  REAL8             f0;	
  REAL8             deltaF;
  CHARVector       *sampleUnits;
  COMPLEX8Sequence *data;
}
COMPLEX8FrequencySeries;

typedef struct
tagCOMPLEX16FrequencySeries
{
  CHAR              *name;
  LIGOTimeGPS        epoch;
  REAL8              f0;	
  REAL8              deltaF;
  CHARVector        *sampleUnits;
  COMPLEX16Sequence *data;
}
COMPLEX16FrequencySeries;

/* Transfer functions */ 

typedef struct
tagCOMPLEX8ZPGFilter
{
  CHAR           *name;
  COMPLEX8Vector *zeros;	
  COMPLEX8Vector *poles;
  COMPLEX8        gain;
}
COMPLEX8ZPGFilter;

typedef struct
tagCOMPLEX16ZPGFilter
{
  CHAR            *name;
  COMPLEX16Vector *zeros;	
  COMPLEX16Vector *poles;
  COMPLEX16        gain;
}
COMPLEX16ZPGFilter;

/* LIGO Function Support: Universal Status Structure */

typedef struct
tagStatus
{
  INT4                 statusCode;
  const CHAR          *statusDescription;
  volatile const CHAR *Id;
  const CHAR          *function;
  const CHAR          *file;
  INT4                 line;
  struct tagStatus    *statusPtr;
  INT4                 level;
}
Status;


#ifdef  __cplusplus
}
#endif

#endif /* _LALDATATYPES_H */
