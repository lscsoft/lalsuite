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
define(`VTYPE',`format(`%sVector',TYPE)')dnl
define(`FUNC',`format(`LAL%sReadVector',TYPECODE)')dnl
define(`PARSEFUNC',`format(`LALStringTo%s',TYPECODE)')dnl
define(`CREATEFUNC',`format(`LAL%sCreateVector',TYPECODE)')dnl
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

/* <lalVerbatim file="StreamVectorInputCP"> */
void
FUNC ( LALStatus *stat, VTYPE **vector, FILE *stream, BOOLEAN strict )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *start, *end;       /* pointers to start and end of a token */
  TYPE *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "FUNC", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Read into first buffer, and see if more needs to be read. */
  start = end = line->data;
  for ( i = 0; more && ( i < BUFFSIZE/SIZE ); i++ ) {
    PARSEFUNC ( stat->statusPtr, head.buf.TYPECODE + i, start, &end );
    BEGINFAIL( stat )
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    ENDFAIL( stat );
    if ( start == end )
      more = 0;
    else {
      nTot++;
      head.size++;
      start = end;
    }
  }

  /* Read into remaining buffers. */
  here = &head;
  while ( more ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
    here = here->next;
    if ( !here ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here->size = 0;
    here->next = NULL;

    /* Read into the next buffer, and see if more needs to be read. */
    for ( i = 0; more && ( i < BUFFSIZE/SIZE ); i++ ) {
      PARSEFUNC ( stat->statusPtr, here->buf.TYPECODE + i, start, &end );
      BEGINFAIL( stat ) {
	FREEBUFFERLIST( head.next );
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      } ENDFAIL( stat );
      if ( start == end )
	more = 0;
      else {
	nTot++;
	here->size++;
	start = end;
      }
    }
  }

  /* Check for formatting problems, if required, and free the line. */
  if ( strict ) {
    if ( nTot == 0 ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
    }
    while ( isspace( *end ) )
      end++;
    if ( *end != '\0' ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EFMT, STREAMINPUTH_MSGEFMT );
    }
  }
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );

  /* Allocate **vector. */
  if ( nTot > 0 ) {
    CREATEFUNC ( stat->statusPtr, vector, nTot );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head.next );
    ENDFAIL( stat );

    /* Copy buffer list into (*vector)->data, and set
       (*vector)->length. */
    here = &head;
    data = (*vector)->data;
    while ( here ) {
      UINT4 j = here->size;
      TYPE *hereData = here->buf.TYPECODE;
      while ( j-- )
	*(data++) = *(hereData++);
      here = here->next;
    }
  }

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
