/**************************** <lalVerbatim file="StreamSequenceInputCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{StreamSequenceInput.c}}
\label{ss:StreamSequenceInput.c}

Converts an input stream into a data sequence.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{StreamSequenceInputCP}
\idx{LALCHARReadSequence()}
\idx{LALI2ReadSequence()}
\idx{LALI4ReadSequence()}
\idx{LALI8ReadSequence()}
\idx{LALU2ReadSequence()}
\idx{LALU4ReadSequence()}
\idx{LALU8ReadSequence()}
\idx{LALSReadSequence()}
\idx{LALDReadSequence()}
\idx{LALCReadSequence()}
\idx{LALZReadSequence()}

\subsubsection*{Description}

These routines read data from the I/O stream \verb@*stream@ until the
end-of-input is reached.  (The input can be of arbitrary length; the
data is temporarily stored in a linked list of buffers.)  Once read, a
LAL sequence structure \verb@**sequence@ is created and the data
stored in it.  The routine passes back a pointer to the new structure.

The routine \verb@LALCHARReadSequence()@ simply stores the entire
remaining contents of the I/O stream in a \verb@CHARSequence@,
including whitespace, newline \verb@'\n'@, null \verb@'\0'@, or other
special characters.  (It can in principle be used to read and store
binary data as a sequence of bytes.  Note that the end-of-transmission
byte \verb@'\004'@ does \emph{not} necessarily mark the end-of-input,
which is instead determined using the \verb@feof()@ function.)

The other routines in this module interpret the input as a sequence of
whitespace-separated numbers, which are parsed directly from the I/O
stream using \verb@fscanf()@.  The sequence is terminated at the
end-of-input or at any point where \verb@fscanf()@ is unable to parse
the input.

For the complex input routines \verb@LALCReadSequence()@ and
\verb@LALZReadSequence()@, each pair of numbers read are interpreted
as the real and imaginary parts of a complex number.  The usual input
format is for each line to contain a pair of numbers, but
\verb@fscanf()@ does not distinguish between newline and other
whitespace characters, so neither do these routines.

Unlike the numerical routines in other \verb@StreamInput.h@ modules,
these routines have no mechanism to deal with comments; every
whitespace-delimited substring will be treated as a number.

\subsubsection*{Algorithm}

These routines read data into a linked list of buffers, to allow
memory allocation to occur in batches for improved efficiency.  The
numerical routines also use \verb@fscanf()@ directly on the I/O stream
to avoid the inefficiency of storing and parsing intermediate
character strings, as is done by the corresponding vector sequence
input routines.  This reduces robustness and versatility (as
indicated, for instance, by the inability of dealing with comments),
and increases the number of potential points-of-failure (by requiring
a consistent implementation across platforms of \verb@getc()@ and
\verb@fscanf()@, rather than the single function \verb@fgets()@ used
by other stream input routines).  However, these sacrifices are
necessary to allow LAL applications to ingest large quantities of
numerical data efficiently.

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                     LALFree()
LALWarning()                    LALCHARCreateVector()
LALI2CreateVector()             LALU2CreateVector()
LALI4CreateVector()             LALU4CreateVector()
LALI8CreateVector()             LALU8CreateVector()
LALSCreateVector()              LALDCreateVector()
LALCCreateVector()              LALZCreateVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StreamSequenceInputCV}}

******************************************************* </lalLaTeX> */

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

/* <lalVerbatim file="StreamSequenceInputCP"> */
void
LALCHARReadSequence( LALStatus *stat, CHARSequence **sequence, FILE *stream )
{ /* </lalVerbatim> */
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


/* <lalVerbatim file="StreamSequenceInputCP"> */
int
XLALCHARReadSequence( CHARSequence **sequence, FILE *stream );
int
XLALCHARReadSequence( CHARSequence **sequence, FILE *stream )
{ /* </lalVerbatim> */
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

define(`TYPECODE',`I2')dnl
include(`LALReadSequence.m4')dnl

define(`TYPECODE',`I4')dnl
include(`LALReadSequence.m4')dnl

define(`TYPECODE',`I8')dnl
include(`LALReadSequence.m4')dnl

define(`TYPECODE',`U2')dnl
include(`LALReadSequence.m4')dnl

define(`TYPECODE',`U4')dnl
include(`LALReadSequence.m4')dnl

define(`TYPECODE',`U8')dnl
include(`LALReadSequence.m4')dnl

define(`TYPECODE',`S')dnl
include(`LALReadSequence.m4')dnl

define(`TYPECODE',`D')dnl
include(`LALReadSequence.m4')dnl

define(`TYPECODE',`C')dnl
include(`LALReadComplexSequence.m4')dnl

define(`TYPECODE',`Z')dnl
include(`LALReadComplexSequence.m4')dnl
