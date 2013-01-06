#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/StreamInput.h>

/* Define linked-list of pointers to vectors of arbitrary type. */
typedef union tagVector {
  CHARVector *CHV;
  INT2Vector *I2V;
  INT4Vector *I4V;
  INT8Vector *I8V;
  UINT2Vector *U2V;
  UINT4Vector *U4V;
  UINT8Vector *U8V;
  REAL4Vector *SV;
  REAL8Vector *DV;
} Vector;
typedef struct tagVectorList {
  Vector vector;
  struct tagVectorList *next;
} VectorList;

static const VectorList empty;


#define FREECHARVECTORLIST                                           \
do {                                                                 \
  if ( head.vector.CHV ) {                                           \
    TRY( LALCHARDestroyVector( stat->statusPtr,                      \
			     &(head.vector.CHV) ), stat );           \
  }                                                                  \
  here = head.next;                                                  \
  while ( here ) {                                                   \
    VectorList *nextPtr = here->next;                                \
    if ( here->vector.CHV ) {                                        \
      TRY( LALCHARDestroyVector( stat->statusPtr,                    \
			       &(here->vector.CHV) ), stat );        \
    }                                                                \
    LALFree( here );                                                 \
    here = nextPtr;                                                  \
  }                                                                  \
} while (0)


void
LALCHARReadVectorSequence( LALStatus          *stat,
			   CHARVectorSequence **sequence,
			   FILE               *stream )
{ 
  VectorList head = empty;   /* head of linked list of vectors */
  VectorList *here;          /* pointer to current position in list */
  CHAR *data;                /* pointer to vector data */
  UINT4 nRows, nCols;        /* number and length of lines */
  CreateVectorSequenceIn in; /* parameters for creating sequence */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Read the first line. */
  if ( !feof( stream ) ) {
    TRY( LALCHARReadVector( stat->statusPtr, &(head.vector.CHV), stream ),
	 stat );
  }
  here = &head;

  /* As long as lines remain... */
  while ( !feof( stream ) ) {

    /* Allocate space for next line. */
    here->next = (VectorList *)LALCalloc( 1, sizeof(VectorList) );
    if ( !here->next ) {
      FREECHARVECTORLIST;
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }

    /* Read in next line. */
    here = here->next;
    LALCHARReadVector( stat->statusPtr, &(here->vector.CHV), stream );
    BEGINFAIL( stat ) {
      FREECHARVECTORLIST;
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    } ENDFAIL( stat );
  }

  /* Lines have been read.  Now determine the maximum line length, and
     allocate the vector sequence.  Ignore lines containing only a
     single '\0' character. */
  nRows = nCols = 0;
  here = &head;
  while ( here ) {
    if ( here->vector.CHV->length > 1 ) {
      nRows++;
      if ( here->vector.CHV->length > nCols )
	nCols = here->vector.CHV->length;
    }
    here = here->next;
  }
  in.length = nRows;
  in.vectorLength = nCols;
  LALCHARCreateVectorSequence( stat->statusPtr, sequence, &in );
  BEGINFAIL( stat ) {
    FREECHARVECTORLIST;
    ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
  } ENDFAIL( stat );

  /* Now assign data to the sequence, padding with zeros as
     necessary. */
  here = &head;
  data = (*sequence)->data;
  while ( here ) {
    UINT4 length = here->vector.CHV->length;
    if ( length > 1 ) {
      memcpy( data, here->vector.CHV->data, length );
      if ( nCols - length > 0 )
	memset( data + length, 0, nCols - length );
      data += nCols;
    }
    here = here->next;
  }

  /* Free memory and exit. */
  FREECHARVECTORLIST;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}

#define TYPECODE I2
#define TYPE INT2
#define SIZE 2
#include "StreamVectorSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE I4
#define TYPE INT4
#define SIZE 4
#include "StreamVectorSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE I8
#define TYPE INT8
#define SIZE 8
#include "StreamVectorSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U2
#define TYPE UINT2
#define SIZE 2
#include "StreamVectorSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U4
#define TYPE UINT4
#define SIZE 4
#include "StreamVectorSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U8
#define TYPE UINT8
#define SIZE 8
#include "StreamVectorSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE S
#define TYPE REAL4
#define SIZE 4
#include "StreamVectorSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE D
#define TYPE REAL8
#define SIZE 8
#include "StreamVectorSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
