/**************************** <lalVerbatim file="StreamVectorInputCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{StreamVectorInput.c}}
\label{ss:StreamVectorInput.c}

Reads data from a single line in an input stream.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{StreamVectorInputCP}
\idx{LALCHARReadVector()}
\idx{LALI2ReadVector()}
\idx{LALI4ReadVector()}
\idx{LALI8ReadVector()}
\idx{LALU2ReadVector()}
\idx{LALU4ReadVector()}
\idx{LALU8ReadVector()}
\idx{LALSReadVector()}
\idx{LALDReadVector()}

\subsubsection*{Description}

These routines read ASCII data from the I/O stream \verb@*stream@
until a newline or the end-of-input is reached.  (The line can be of
arbitrary length; the data is temporarily stored in a linked list of
buffers.)  Once read, a LAL vector structure \verb@**vector@ is
created and the data stored in it.  The routine passes back a pointer
to the new structure.  For the numerical routines, the \verb@strict@
parameter determines whether the routine will do strict error checking
based on the contents of the input stream (see below).

The basic routine in this module is \verb@LALCHARReadVector()@, which
simply stores bytes read from \verb@*stream@ until the next newline
character \verb@'\n'@, null character \verb@'\0'@, or the end of the
input as determined by the \verb@feof()@ function.  The vector
includes the newline (if present), and also an explicit \verb@'\0'@ at
the end, if one was not already present.  This routine should
\emph{not} be used to read a binary data stream, which are not
logically divided into `lines'.  Unless it aborts due to invalid
arguments or failed memory allocation, \verb@LALCHARReadVector()@ will
always return successfully regardless of the contents of the input
stream; \verb@*vector@ will created containing at least a single
\verb@'\0'@ terminator, if nothing else.

The other routines in this module use \verb@LALCHARReadVector()@ to
read a line, and then parse it into numerical datatypes using the
corresponding routine in the \verb@StringConvert.c@ module.
Conversion stops when the routine encounters a character that cannot
be parsed as part of a number.  If \verb@strict@ is 0, the routine
will fail only due to invalid arguments or memory allocation failure,
not from a poorly-formatted input stream; if no numbers are read,
\verb@*vector@ will remain \verb@NULL@, but no error will be reported.
(In this mode, the calling routine should always test the output
before trying to dereference it, in order to avoid segmentation
violations.)  If \verb@strict@ is nonzero, the routine will report an
error if the input stream was poorly formatted, either an \verb@ELEN@
error if no numbers were read, or \verb@EFMT@ if a character was
encountered that was neither part of a parseable number nor
whitespace.

Note that \verb@strict@=0 allows an input stream to contain blank
lines or comments.  A comment begins with any character that cannot
occur in a valid number, which will cause the numerical parser to skip
the rest of the line.  The usual comment delimiters are \verb@'#'@ and
\verb@'%'@, but any character except \verb@'+'@ \verb@'-'@,
\verb@'e'@, \verb@'E'@, \verb@'.'@, digits, and whitespace will work.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                     LALFree()
LALCHARCreateVector()           LALCHARDestroyVector()
LALI2CreateVector()             LALU2CreateVector()
LALI4CreateVector()             LALU4CreateVector()
LALI8CreateVector()             LALU8CreateVector()
LALSCreateVector()              LALDCreateVector()
LALWarning()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StreamVectorInputCV}}

******************************************************* </lalLaTeX> */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>
#include <lal/StreamInput.h>

NRCSID( STREAMVECTORINPUTC, "$Id$" );

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

/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALCHARReadVector( LALStatus *stat, CHARVector **vector, FILE *stream )
{ /* </lalVerbatim> */
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;     /* pointer to current position in list */
  CHAR *data;           /* pointer to vector data */
  BOOLEAN done = 0;     /* whether or not to read more buffers */
  size_t nTot;          /* total number of characters read */

  INITSTATUS( stat, "LALCHARReadVector", STREAMVECTORINPUTC );
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


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALI2ReadVector ( LALStatus *stat, INT2Vector **vector, FILE *stream, BOOLEAN strict )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *start, *end;       /* pointers to start and end of a token */
  INT2 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALI2ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Read into first buffer, and see if more needs to be read. */
  start = end = line->data;
  for ( i = 0; more && ( i < BUFFSIZE/2 ); i++ ) {
    LALStringToI2 ( stat->statusPtr, head.buf.I2 + i, start, &end );
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
    for ( i = 0; more && ( i < BUFFSIZE/2 ); i++ ) {
      LALStringToI2 ( stat->statusPtr, here->buf.I2 + i, start, &end );
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
    LALI2CreateVector ( stat->statusPtr, vector, nTot );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head.next );
    ENDFAIL( stat );

    /* Copy buffer list into (*vector)->data, and set
       (*vector)->length. */
    here = &head;
    data = (*vector)->data;
    while ( here ) {
      UINT4 j = here->size;
      INT2 *hereData = here->buf.I2;
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


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALI4ReadVector ( LALStatus *stat, INT4Vector **vector, FILE *stream, BOOLEAN strict )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *start, *end;       /* pointers to start and end of a token */
  INT4 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALI4ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Read into first buffer, and see if more needs to be read. */
  start = end = line->data;
  for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ ) {
    LALStringToI4 ( stat->statusPtr, head.buf.I4 + i, start, &end );
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
    for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ ) {
      LALStringToI4 ( stat->statusPtr, here->buf.I4 + i, start, &end );
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
    LALI4CreateVector ( stat->statusPtr, vector, nTot );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head.next );
    ENDFAIL( stat );

    /* Copy buffer list into (*vector)->data, and set
       (*vector)->length. */
    here = &head;
    data = (*vector)->data;
    while ( here ) {
      UINT4 j = here->size;
      INT4 *hereData = here->buf.I4;
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


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALI8ReadVector ( LALStatus *stat, INT8Vector **vector, FILE *stream, BOOLEAN strict )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *start, *end;       /* pointers to start and end of a token */
  INT8 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALI8ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Read into first buffer, and see if more needs to be read. */
  start = end = line->data;
  for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ ) {
    LALStringToI8 ( stat->statusPtr, head.buf.I8 + i, start, &end );
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
    for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ ) {
      LALStringToI8 ( stat->statusPtr, here->buf.I8 + i, start, &end );
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
    LALI8CreateVector ( stat->statusPtr, vector, nTot );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head.next );
    ENDFAIL( stat );

    /* Copy buffer list into (*vector)->data, and set
       (*vector)->length. */
    here = &head;
    data = (*vector)->data;
    while ( here ) {
      UINT4 j = here->size;
      INT8 *hereData = here->buf.I8;
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


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALU2ReadVector ( LALStatus *stat, UINT2Vector **vector, FILE *stream, BOOLEAN strict )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *start, *end;       /* pointers to start and end of a token */
  UINT2 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALU2ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Read into first buffer, and see if more needs to be read. */
  start = end = line->data;
  for ( i = 0; more && ( i < BUFFSIZE/2 ); i++ ) {
    LALStringToU2 ( stat->statusPtr, head.buf.U2 + i, start, &end );
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
    for ( i = 0; more && ( i < BUFFSIZE/2 ); i++ ) {
      LALStringToU2 ( stat->statusPtr, here->buf.U2 + i, start, &end );
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
    LALU2CreateVector ( stat->statusPtr, vector, nTot );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head.next );
    ENDFAIL( stat );

    /* Copy buffer list into (*vector)->data, and set
       (*vector)->length. */
    here = &head;
    data = (*vector)->data;
    while ( here ) {
      UINT4 j = here->size;
      UINT2 *hereData = here->buf.U2;
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


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALU4ReadVector ( LALStatus *stat, UINT4Vector **vector, FILE *stream, BOOLEAN strict )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *start, *end;       /* pointers to start and end of a token */
  UINT4 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALU4ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Read into first buffer, and see if more needs to be read. */
  start = end = line->data;
  for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ ) {
    LALStringToU4 ( stat->statusPtr, head.buf.U4 + i, start, &end );
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
    for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ ) {
      LALStringToU4 ( stat->statusPtr, here->buf.U4 + i, start, &end );
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
    LALU4CreateVector ( stat->statusPtr, vector, nTot );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head.next );
    ENDFAIL( stat );

    /* Copy buffer list into (*vector)->data, and set
       (*vector)->length. */
    here = &head;
    data = (*vector)->data;
    while ( here ) {
      UINT4 j = here->size;
      UINT4 *hereData = here->buf.U4;
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


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALU8ReadVector ( LALStatus *stat, UINT8Vector **vector, FILE *stream, BOOLEAN strict )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *start, *end;       /* pointers to start and end of a token */
  UINT8 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALU8ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Read into first buffer, and see if more needs to be read. */
  start = end = line->data;
  for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ ) {
    LALStringToU8 ( stat->statusPtr, head.buf.U8 + i, start, &end );
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
    for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ ) {
      LALStringToU8 ( stat->statusPtr, here->buf.U8 + i, start, &end );
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
    LALU8CreateVector ( stat->statusPtr, vector, nTot );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head.next );
    ENDFAIL( stat );

    /* Copy buffer list into (*vector)->data, and set
       (*vector)->length. */
    here = &head;
    data = (*vector)->data;
    while ( here ) {
      UINT4 j = here->size;
      UINT8 *hereData = here->buf.U8;
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


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALSReadVector ( LALStatus *stat, REAL4Vector **vector, FILE *stream, BOOLEAN strict )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *start, *end;       /* pointers to start and end of a token */
  REAL4 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALSReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Read into first buffer, and see if more needs to be read. */
  start = end = line->data;
  for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ ) {
    LALStringToS ( stat->statusPtr, head.buf.S + i, start, &end );
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
    for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ ) {
      LALStringToS ( stat->statusPtr, here->buf.S + i, start, &end );
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
    LALSCreateVector ( stat->statusPtr, vector, nTot );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head.next );
    ENDFAIL( stat );

    /* Copy buffer list into (*vector)->data, and set
       (*vector)->length. */
    here = &head;
    data = (*vector)->data;
    while ( here ) {
      UINT4 j = here->size;
      REAL4 *hereData = here->buf.S;
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


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALDReadVector ( LALStatus *stat, REAL8Vector **vector, FILE *stream, BOOLEAN strict )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty; /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *start, *end;       /* pointers to start and end of a token */
  REAL8 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALDReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Read into first buffer, and see if more needs to be read. */
  start = end = line->data;
  for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ ) {
    LALStringToD ( stat->statusPtr, head.buf.D + i, start, &end );
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
    for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ ) {
      LALStringToD ( stat->statusPtr, here->buf.D + i, start, &end );
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
    LALDCreateVector ( stat->statusPtr, vector, nTot );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head.next );
    ENDFAIL( stat );

    /* Copy buffer list into (*vector)->data, and set
       (*vector)->length. */
    here = &head;
    data = (*vector)->data;
    while ( here ) {
      UINT4 j = here->size;
      REAL8 *hereData = here->buf.D;
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
