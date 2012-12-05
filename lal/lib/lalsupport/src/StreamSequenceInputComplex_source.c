#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define VTYPE CONCAT2(TYPE,Sequence)
#define FUNC CONCAT3(LAL,TYPECODE,ReadSequence)
#define CREATEFUNC CONCAT3(LAL,TYPECODE,CreateVector)


void
FUNC ( LALStatus *stat, VTYPE **sequence, FILE *stream )
{ 
  BufferList head;    /* head of linked list of buffers */
  BufferList *here;   /* pointer to current position in list */
  BOOLEAN done;       /* whether to stop reading */
  TYPE *data;         /* pointer to vector data */
  size_t nTot = 0;    /* total number of values read */

  INITSTATUS(stat);
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
      DATA re = 0.0, im = 0.0;
      done = ( ( fscanf( stream, "%" FORMAT, &re ) != 1 ) ||
          ( fscanf( stream, "%" FORMAT, &im ) != 1 ) );
      *data = re + I * im;
      data++;
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
