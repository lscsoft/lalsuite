dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')dnl
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')dnl
ifelse(TYPECODE,`Z',`define(`FMT',`LAL_REAL8_FORMAT')')dnl
ifelse(TYPECODE,`C',`define(`FMT',`LAL_REAL4_FORMAT')')dnl
define(`VTYPE',`format(`%sSequence',TYPE)')dnl
define(`FUNC',`format(`LAL%sReadSequence',TYPECODE)')dnl
define(`CREATEFUNC',`format(`LAL%sCreateVector',TYPECODE)')dnl
ifelse(TYPECODE,`Z',`define(`SIZE',`16')')dnl
ifelse(TYPECODE,`C',`define(`SIZE',`8')')dnl
define(`PARSEDATA',`done = ( ( fscanf( stream, "%" FMT, &(data->re) ) != 1 ) ||
                ( fscanf( stream, "%" FMT, &(data->im) ) != 1 ) ); data++')dnl

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
    while ( !done && --n ) {
      PARSEDATA;
    }
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
