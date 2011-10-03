/**
\author Creighton, T. D.
\file
*/

/**

\heading{Module \ref StreamSequenceInput.c}
\latexonly\label{ss_StreamSequenceInput_c}\endlatexonly

Converts an input stream into a data sequence.

\heading{Prototypes}














\heading{Description}

These routines read data from the I/O stream <tt>*stream</tt> until the
end-of-input is reached.  (The input can be of arbitrary length; the
data is temporarily stored in a linked list of buffers.)  Once read, a
LAL sequence structure <tt>**sequence</tt> is created and the data
stored in it.  The routine passes back a pointer to the new structure.

The routine <tt>LALCHARReadSequence()</tt> simply stores the entire
remaining contents of the I/O stream in a \c CHARSequence,
including whitespace, newline <tt>'\n'</tt>, null <tt>'\0'</tt>, or other
special characters.  (It can in principle be used to read and store
binary data as a sequence of bytes.  Note that the end-of-transmission
byte <tt>'\004'</tt> does \e not necessarily mark the end-of-input,
which is instead determined using the <tt>feof()</tt> function.)

The other routines in this module interpret the input as a sequence of
whitespace-separated numbers, which are parsed directly from the I/O
stream using <tt>fscanf()</tt>.  The sequence is terminated at the
end-of-input or at any point where <tt>fscanf()</tt> is unable to parse
the input.

For the complex input routines <tt>LALCReadSequence()</tt> and
<tt>LALZReadSequence()</tt>, each pair of numbers read are interpreted
as the real and imaginary parts of a complex number.  The usual input
format is for each line to contain a pair of numbers, but
<tt>fscanf()</tt> does not distinguish between newline and other
whitespace characters, so neither do these routines.

Unlike the numerical routines in other \ref StreamInput.h modules,
these routines have no mechanism to deal with comments; every
whitespace-delimited substring will be treated as a number.

\heading{Algorithm}

These routines read data into a linked list of buffers, to allow
memory allocation to occur in batches for improved efficiency.  The
numerical routines also use <tt>fscanf()</tt> directly on the I/O stream
to avoid the inefficiency of storing and parsing intermediate
character strings, as is done by the corresponding vector sequence
input routines.  This reduces robustness and versatility (as
indicated, for instance, by the inability of dealing with comments),
and increases the number of potential points-of-failure (by requiring
a consistent implementation across platforms of <tt>getc()</tt> and
<tt>fscanf()</tt>, rather than the single function <tt>fgets()</tt> used
by other stream input routines).  However, these sacrifices are
necessary to allow LAL applications to ingest large quantities of
numerical data efficiently.

\heading{Uses}
\code
LALMalloc()                     LALFree()
LALWarning()                    LALCHARCreateVector()
LALI2CreateVector()             LALU2CreateVector()
LALI4CreateVector()             LALU4CreateVector()
LALI8CreateVector()             LALU8CreateVector()
LALSCreateVector()              LALDCreateVector()
LALCCreateVector()              LALZCreateVector()
\endcode

\heading{Notes}



*/

#include <stdio.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>
#include <lal/StreamInput.h>

NRCSID( STREAMSEQUENCEINPUTC, "$Id$" );

/* Define linked-list of buffers for storing an arbitrary number of
   arbitrary datatypes. */
#define BUFFSIZE 1024
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
  COMPLEX8 C[BUFFSIZE/8];
  COMPLEX16 Z[BUFFSIZE/16];
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

/* static const BufferList empty; */


void
LALCHARReadSequence( LALStatus *stat, CHARSequence **sequence, FILE *stream )
{ 
  BufferList head;  /* head of linked list of buffers */
  BufferList *here; /* pointer to current position in list */
  CHAR *data;       /* pointer to vector data */
  size_t nTot = 0;  /* total number of characters read */

  INITSTATUS( stat, "LALCHARReadSequence", STREAMSEQUENCEINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( sequence, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*sequence, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read file into linked list of buffers. */
  here = &head;
  here->next = NULL;
  while ( !feof( stream ) ) {
    size_t n = BUFFSIZE;
    data = here->buf.CH;
    while ( !feof( stream ) && n )
	{
	   *(data++) = (CHAR)getc( stream );
	    n--;
	}
    /* The very last value returned by getc() is EOF, which is not a
       character and should not be stored. */
    if ( feof( stream ) ) {
      data--;
      n++;
    }
    here->size = BUFFSIZE - n;
    nTot += here->size;
    if ( !feof( stream ) ) {
      here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
      if ( !(here->next) ) {
	FREEBUFFERLIST( head.next );
	ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
      }
      here = here->next;
      here->next = NULL;
    }
  }

  /* Allocate **sequence, include space for a terminating '\0'. */
  LALCHARCreateVector( stat->statusPtr, sequence, nTot+1 );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*sequence)->data, adding final '\0'. */
  here = &head;
  data = (*sequence)->data;
  while ( here ) {
    memcpy( data, here->buf.CH, here->size );
    data += here->size;
    here = here->next;
  }
  (*sequence)->data[nTot] = '\0';

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}



int
XLALCHARReadSequence( CHARSequence **sequence, FILE *stream );
int
XLALCHARReadSequence( CHARSequence **sequence, FILE *stream )
{ 
  BufferList head;  /* head of linked list of buffers */
  BufferList *here; /* pointer to current position in list */
  CHAR *data;       /* pointer to vector data */
  size_t nTot = 0;  /* total number of characters read */

  /* Check for valid input arguments. */
  if ( !stream || !sequence ) {
    fprintf(stderr, STREAMINPUTH_MSGENUL );
    return STREAMINPUTH_ENUL;
  }
  if ( *sequence != NULL ) {
    fprintf(stderr, STREAMINPUTH_MSGEOUT );
    return STREAMINPUTH_EOUT;
  }

  /* Read file into linked list of buffers. */
  here = &head;
  here->next = NULL;
  while ( !feof( stream ) ) {
    size_t n = BUFFSIZE;
    data = here->buf.CH;
    while ( !feof( stream ) && n )
	{
	   *(data++) = (CHAR)getc( stream );
	    n--;
	}
    /* The very last value returned by getc() is EOF, which is not a
       character and should not be stored. */
    if ( feof( stream ) ) {
      data--;
      n++;
    }
    here->size = BUFFSIZE - n;
    nTot += here->size;
    if ( !feof( stream ) ) {
      here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
      if ( !(here->next) ) {
	FREEBUFFERLIST( head.next );
	fprintf(stderr, STREAMINPUTH_MSGEMEM );
	return STREAMINPUTH_EMEM;
      }
      here = here->next;
      here->next = NULL;
    }
  }

  /* Allocate **sequence, include space for a terminating '\0'. */
  *sequence = XLALCreateCHARVector( nTot+1 );
  if ( ! *sequence ) {
    FREEBUFFERLIST( head.next );
  }

  /* Copy buffer list into (*sequence)->data, adding final '\0'. */
  here = &head;
  data = (*sequence)->data;
  while ( here ) {
    memcpy( data, here->buf.CH, here->size );
    data += here->size;
    here = here->next;
  }
  (*sequence)->data[nTot] = '\0';

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  
  return 0;
}


/* tell the GNU compiler to ignore issues with the `ll' length modifier */
#ifdef __GNUC__
#define fscanf __extension__ fscanf
#endif

#define TYPECODE I2
#define TYPE INT2
#define SIZE 2
#include "StreamSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE I4
#define TYPE INT4
#define SIZE 4
#include "StreamSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE I8
#define TYPE INT8
#define SIZE 8
#include "StreamSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U2
#define TYPE UINT2
#define SIZE 2
#include "StreamSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U4
#define TYPE UINT4
#define SIZE 4
#include "StreamSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U8
#define TYPE UINT8
#define SIZE 8
#include "StreamSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE S
#define TYPE REAL4
#define SIZE 4
#include "StreamSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE D
#define TYPE REAL8
#define SIZE 8
#include "StreamSequenceInput_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE Z
#define TYPE COMPLEX16
#define SIZE 16
#define FORMAT LAL_REAL8_FORMAT
#include "StreamSequenceInputComplex_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef FORMAT

#define TYPECODE C
#define TYPE COMPLEX8
#define SIZE 8
#define FORMAT LAL_REAL4_FORMAT
#include "StreamSequenceInputComplex_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef FORMAT
