dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')dnl
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')dnl
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')dnl
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')dnl
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')dnl
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')dnl
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')dnl
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')dnl
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')dnl
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')dnl
define(`VTYPE',`format(`%sVectorSequence',TYPE)')dnl
define(`VTYPECODE',`format(`%sV',TYPECODE)')dnl
define(`FUNC',`format(`LAL%sReadVectorSequence',TYPECODE)')dnl
define(`VFUNC',`format(`LAL%sReadVector',TYPECODE)')dnl
define(`CREATEFUNC',`format(`LAL%sCreateVectorSequence',TYPECODE)')dnl
define(`DESTROYFUNC',`format(`LAL%sDestroyVector',TYPECODE)')dnl
define(`FREEMACRO',`format(`FREE%sVECTORLIST',TYPECODE)')dnl
ifelse(TYPECODE,`Z',`define(`SIZE',`16')')dnl
ifelse(TYPECODE,`C',`define(`SIZE',`8')')dnl
ifelse(TYPECODE,`D',`define(`SIZE',`8')')dnl
ifelse(TYPECODE,`S',`define(`SIZE',`4')')dnl
ifelse(TYPECODE,`I2',`define(`SIZE',`2')')dnl
ifelse(TYPECODE,`I4',`define(`SIZE',`4')')dnl
ifelse(TYPECODE,`I8',`define(`SIZE',`8')')dnl
ifelse(TYPECODE,`U2',`define(`SIZE',`2')')dnl
ifelse(TYPECODE,`U4',`define(`SIZE',`4')')dnl
ifelse(TYPECODE,`U8',`define(`SIZE',`8')')dnl

#define \
FREEMACRO \
do { \
  if ( head.vector.VTYPECODE ) { \
    TRY( DESTROYFUNC ( stat->statusPtr, &(head.vector.VTYPECODE) ), stat ); \
  } \
  here = head.next; \
  while ( here ) { \
    VectorList *nextPtr = here->next; \
    if ( here->vector.VTYPECODE ) { \
      TRY( DESTROYFUNC ( stat->statusPtr, &(here->vector.VTYPECODE) ), stat ); \
    } \
    LALFree( here ); \
    here = nextPtr; \
  } \
} while (0)

/* <lalVerbatim file="StreamVectorSequenceInputCP"> */
void
FUNC ( LALStatus  *stat, VTYPE **sequence, FILE *stream )
{ /* </lalVerbatim> */
  VectorList head = empty; /* head of linked list of vectors */
  VectorList *here;     /* pointer to current position in list */
  TYPE *data;           /* pointer to vector data */
  UINT4 nRows, nCols;   /* number and length of lines */
  CreateVectorSequenceIn in; /* parameters for creating sequence */

  INITSTATUS( stat, "FUNC", STREAMVECTORSEQUENCEINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Read the first line. */
  if ( !feof( stream ) ) {
    TRY( VFUNC ( stat->statusPtr, &(head.vector.VTYPECODE), stream, 0 ),
	 stat );
  }
  here = &head;

  /* As long as lines remain... */
  while ( !feof( stream ) ) {

    /* Allocate space for next line. */
    here->next = (VectorList *)LALCalloc( 1, sizeof(VectorList) );
    here = here->next;
    if ( !here ) {
      FREEMACRO;
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }

    /* Read in next line. */
    VFUNC ( stat->statusPtr, &(here->vector.VTYPECODE), stream, 0 );
    BEGINFAIL( stat ) {
      FREEMACRO;
    } ENDFAIL( stat );
  }

  /* Lines have been read.  Now determine minimum common line length,
     and allocate the vector sequence. */
  nRows = nCols = 0;
  here = &head;
  while ( here ) {
    if ( here->vector.VTYPECODE ) {
      if ( here->vector.VTYPECODE->length > nCols )
	nCols = here->vector.VTYPECODE->length;
      nRows++;
    }
    here = here->next;
  }
  if ( ( nRows == 0 ) || ( nCols == 0 ) ) {
    FREEMACRO;
    ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
  }
  in.length = nRows;
  in.vectorLength = nCols;
  CREATEFUNC ( stat->statusPtr, sequence, &in );
  BEGINFAIL( stat ) {
    FREEMACRO;
    ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
  } ENDFAIL( stat );

  /* Now assign data to the sequence, padding with zeros as
     necessary. */
  here = &head;
  data = (*sequence)->data;
  while ( here ) {
    if ( here->vector.VTYPECODE ) {
      UINT4 i;
      UINT4 length = here->vector.VTYPECODE->length;
      TYPE *hereData = here->vector.VTYPECODE->data;
      for ( i = 0; i < length; i++ )
	*(data++) = *(hereData++);
      if ( nCols - length > 0 )
	memset( data, 0, SIZE*( nCols - length ) );
      data += nCols - length;
    }
    here = here->next;
  }

  /* Free memory and exit. */
  FREEMACRO;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
