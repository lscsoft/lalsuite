#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>
#include <lal/StreamInput.h>

/* Define linked-list of buffers for storing an arbitrary number of
   arbitrary datatypes. */
#define BUFFSIZE 16
typedef union tagBuffer {
  CHAR CH[BUFFSIZE];
  INT2 I2[BUFFSIZE/2];
  INT4 I4[BUFFSIZE/4];
  INT8 I8[BUFFSIZE/8];
  UINT2 U2[BUFFSIZE/2];
  UINT4 U4[BUFFSIZE/4];
  UINT8 U8[BUFFSIZE/8];
  REAL4 S[BUFFSIZE/4];
  REAL8 D[BUFFSIZE/8];
} Buffer;
typedef struct tagBufferList {
  Buffer buf;
  size_t size;
  struct tagBufferList *next;
} BufferList;

/* Define a macro for freeing the linked list. */
#define FREEBUFFERLIST( headPtr )                                    \
if ( headPtr ) {                                                     \
  BufferList *herePtr = headPtr;                                     \
  while ( herePtr ) {                                                \
    BufferList *nextPtr = herePtr->next;                             \
    LALFree( herePtr );                                              \
    herePtr = nextPtr;                                               \
  }                                                                  \
} else (void)(0)

static const BufferList empty;


void
LALCHARReadVector( LALStatus *stat, CHARVector **vector, FILE *stream )
{ 
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;     /* pointer to current position in list */
  CHAR *data;           /* pointer to vector data */
  BOOLEAN done = 0;     /* whether or not to read more buffers */
  size_t nTot;          /* total number of characters read */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read into the first buffer at the head of the list, and see if
     more needs to be read. */
  if ( !( fgets( head.buf.CH, BUFFSIZE, stream ) ) )
    done = 1;
  nTot = head.size = strlen( head.buf.CH );
  done |= ( head.size < BUFFSIZE - 1 );
  done |= ( head.buf.CH[BUFFSIZE-2] == '\n' );
  here = &head;

  /* If we haven't yet reached the end of the line or file... */
  while ( !done ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
    here = here->next;
    if ( !here ) {
      FREEBUFFERLIST( head.next );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    memset( here, 0, sizeof(BufferList) );

    /* Read into the next buffer, and see if more needs to be read. */
    if ( !( fgets( here->buf.CH, BUFFSIZE, stream ) ) )
      done = 1;
    fflush( stream );
    nTot += here->size = strlen( here->buf.CH );
    done |= ( here->size < BUFFSIZE - 1 );
    done |= ( here->buf.CH[BUFFSIZE-2] == '\n' );
  }

  /* Finished reading the line.  Now allocate **vector.  Include space
     for a terminating '\0'. */
  LALCHARCreateVector( stat->statusPtr, vector, nTot+1 );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*vector)->data, adding final '\0'. */
  here = &head;
  data = (*vector)->data;
  while ( here ) {
    memcpy( data, here->buf.CH, here->size );
    data += here->size;
    here = here->next;
  }
  (*vector)->data[nTot] = '\0';

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}

#define TYPECODE I2
#define TYPE INT2
#define SIZE 2
#include "StreamVectorInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE I4
#define TYPE INT4
#define SIZE 4
#include "StreamVectorInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE I8
#define TYPE INT8
#define SIZE 8
#include "StreamVectorInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U2
#define TYPE UINT2
#define SIZE 2
#include "StreamVectorInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U4
#define TYPE UINT4
#define SIZE 4
#include "StreamVectorInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U8
#define TYPE UINT8
#define SIZE 8
#include "StreamVectorInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE S
#define TYPE REAL4
#define SIZE 4
#include "StreamVectorInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE D
#define TYPE REAL8
#define SIZE 8
#include "StreamVectorInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
