dnl $Id$
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')dnl
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')dnl
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')dnl
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')dnl
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')dnl
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')dnl
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')dnl
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')dnl
define(`VTYPE',`format(`%sSequence',TYPE)')dnl
define(`FUNC',`format(`LAL%sReadSequence',TYPECODE)')dnl
define(`CREATEFUNC',`format(`LAL%sCreateVector',TYPECODE)')dnl
define(`FMT',`format(`LAL_%s_FORMAT',TYPE)')dnl
define(`PARSEDATA',`done = ( fscanf( stream, "%" FMT, data++ ) != 1 )')dnl
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

/* <lalVerbatim file="StreamSequenceInputCP"> */
void
FUNC ( LALStatus *stat, VTYPE **sequence, FILE *stream )
{ /* </lalVerbatim> */
  BufferList head;    /* head of linked list of buffers */
  BufferList *here;   /* pointer to current position in list */
  BOOLEAN done;       /* whether to stop reading */
  TYPE *data;         /* pointer to vector data */
  size_t nTot = 0;    /* total number of values read */

  INITSTATUS( stat, "FUNC", STREAMSEQUENCEINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( sequence, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*sequence, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read file into linked list of buffers. */
  here = &head;
  here->next = NULL;
  done = feof( stream );
  while ( !done ) {
    size_t n = BUFFSIZE/SIZE + 1;
    data = here->buf.TYPECODE;
    while ( !done && --n )
      PARSEDATA;
    here->size = BUFFSIZE/SIZE - n;
    nTot += here->size;
    if ( !done ) {
      here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
      if ( !(here->next) ) {
	FREEBUFFERLIST( head.next );
	ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
      }
      here = here->next;
      here->next = NULL;
    }
  }

  /* If anything was read, allocate **sequence. */
  if ( !nTot ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
  }
  CREATEFUNC ( stat->statusPtr, sequence, nTot );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*sequence)->data. */
  here = &head;
  data = (*sequence)->data;
  while ( here ) {
    memcpy( data, here->buf.TYPECODE, SIZE*here->size );
    data += here->size;
    here = here->next;
  }

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
