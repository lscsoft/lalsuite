/**************************** <lalVerbatim file="StreamVectorInputCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{StreamVectorInput.c}}
\label{ss:StreamVectorInput.c}

Reads data from a single line in a file.

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

These routines read data from the I/O stream \verb@*stream@ until an
end-of-line or end-of-file is reached.  (The line can be of arbitrary
length; the data is temporarily stored in a linked list of buffers.)
Once read, a LAL vector structure \verb@**vector@ is created and the
data stored in it.  The routine passes back a pointer to the new
structure.  The pointer \verb@*end@, if non-null, stores the last
character read; this is useful for discerning the end of I/O streams
(such as standard input) that do not explicitly set an end-of-file
flag.

The basic routine in this module is \verb@LALCHARReadVector()@, which
simply stores bytes read from \verb@*stream@ until the next
end-of-line byte, or the end of the file.  The vector includes the
end-of-line byte (if present), and also an explicit \verb@'\0'@ byte
at the end.  This routine \emph{can} be used to read a binary data
stream, but this is not particularly meaningful since binary files are
not logically divided into ``lines''.

The other routines in this module use \verb@LALCHARReadVector()@ to
read a line, and then parse it into numerical datatypes using the C
routine \verb@sscanf()@.  The stream must be an ASCII text stream
containing whitespace-separated numbers in a format recognized by
\verb@sscanf()@.  A \verb@#@ sign at the beginning of a line, or a
\verb@%@ sign anywhere in the line, indicates that the remainder of
the line is a comment and will be ignored.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()                     LALFree()
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
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include "StreamInput.h"

NRCSID(STREAMVECTORINPUTC,"$Id$");

/* Define macros for determining the sscanf() format specifier. */
#define INTFORMAT( format, size )                                    \
do {                                                                 \
  if ( (size) == sizeof( short ) )                                   \
    (format) = "%hi";                                                \
  else if ( (size) == sizeof( int ) )                                \
    (format) = "%i";                                                 \
  else if ( (size) == sizeof( long ) )                               \
    (format) = "%li";                                                \
  else if ( (size) == sizeof( long long ) )                          \
    (format) = "%Li";                                                \
  else {                                                             \
    CHAR msg[64];                                                    \
    sprintf( msg, "Using default format %%i for type INT%u",         \
	     (size) );                                               \
    LALWarning( stat, msg );                                         \
    (format) = "%i";                                                 \
  }                                                                  \
} while(0)

#define UINTFORMAT( format, size )                                   \
do {                                                                 \
  if ( (size) == sizeof( unsigned short ) )                          \
    (format) = "%hu";                                                \
  else if ( (size) == sizeof( unsigned int ) )                       \
    (format) = "%u";                                                 \
  else if ( (size) == sizeof( unsigned long ) )                      \
    (format) = "%lu";                                                \
  else if ( (size) == sizeof( unsigned long long ) )                 \
    (format) = "%Lu";                                                \
  else {                                                             \
    CHAR msg[64];                                                    \
    sprintf( msg, "Using default format %%u for type INT%u",         \
	     (size) );                                               \
    LALWarning( stat, msg );                                         \
    (format) = "%u";                                                 \
  }                                                                  \
} while(0)

#define REALFORMAT( format, size )                                   \
do {                                                                 \
  if ( (size) == sizeof( float ) )                                   \
    (format) = "%f";                                                 \
  else if ( (size) == sizeof( double ) )                             \
    (format) = "%lf";                                                \
  else if ( (size) == sizeof( long double ) )                        \
    (format) = "%Lf";                                                \
  else {                                                             \
    CHAR msg[64];                                                    \
    sprintf( msg, "Using default format %%f for type REAL%u",        \
	     (size) );                                               \
    LALWarning( stat, msg );                                         \
    (format) = "%f";                                                 \
  }                                                                  \
} while(0)

/* Define linked-list of buffers for storing an arbitrary number of
   arbitrary datatypes. */
#define BUFFSIZE 16
typedef union tagBuffer {
  CHAR ch[BUFFSIZE];
  INT2 i2[BUFFSIZE/2];
  INT4 i4[BUFFSIZE/4];
  INT8 i8[BUFFSIZE/8];
  UINT2 u2[BUFFSIZE/2];
  UINT4 u4[BUFFSIZE/4];
  UINT8 u8[BUFFSIZE/8];
  REAL4 s[BUFFSIZE/4];
  REAL8 d[BUFFSIZE/8];
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

/* Define the set of whitespace characters. */
static char
whitespace[] = { ' ', '\f', '\n', '\r', '\t', '\v', EOF, '\0' };

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
  if ( !( fgets( head.buf.ch, BUFFSIZE, stream ) ) )
    done = 1;
  nTot = head.size = strlen( head.buf.ch );
  done |= ( head.size < BUFFSIZE - 1 );
  done |= ( head.buf.ch[BUFFSIZE-2] == '\n' );
  here = &head;

  /* If we haven't yet read the EOL or EOF character... */
  while ( !done ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALCalloc( 1, sizeof(BufferList) );
    if ( !here->next ) {
      FREEBUFFERLIST( head.next );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here = here->next;

    /* Read into the next buffer, and see if more needs to be read. */
    if ( !( fgets( here->buf.ch, BUFFSIZE, stream ) ) )
      done = 1;
    fflush( stream );
    nTot += here->size = strlen( here->buf.ch );
    done |= ( here->size < BUFFSIZE - 1 );
    done |= ( here->buf.ch[BUFFSIZE-2] == '\n' );
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
    memcpy( data, here->buf.ch, here->size );
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
LALI2ReadVector( LALStatus *stat, INT2Vector **vector, FILE *stream )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  const CHAR *format;      /* the input format specifier */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty;    /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *token;             /* pointer to a number to be converted */
  INT2 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALI2ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Determine the format for sscanf(). */
  INTFORMAT( format, 2 );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Scan for comment delimiters. */
  if ( line->data[0] == '#' )
    line->data[0] = '\0';
  else
    for ( i = 0; i < line->length; i++ )
      if ( line->data[i] == '%' )
	line->data[i] = '\0';

  /* Read into the first buffer. */
  token = strtok( line->data, whitespace );
  for ( i = 0; more && ( i < BUFFSIZE/2 ); i++ )
    if ( token )
    {
      if ( sscanf( token, format, head.buf.i2 + i ) != 1 )
	more = 0;
      else {
	nTot++;
	head.size++;
	token = strtok( NULL, whitespace );
      }
    }
    else
      more = 0;

  /* Read into remaining buffers. */
  here = &head;
  while ( more ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALCalloc( 1, sizeof(BufferList) );
    if ( !here->next ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here = here->next;

    /* Read into the next buffer, and see if more needs to be read. */
    for ( i = 0; more && ( i < BUFFSIZE/2 ); i++ )
      if ( token )
      {
	if ( sscanf( token, format, here->buf.i2 + i ) != 1 )
	  more = 0;
	else {
	  nTot++;
	  here->size++;
	  token = strtok( NULL, whitespace );
	}
      }
      else
	more = 0;
  }

  /* Line has been parsed into numbers, so free it and return an error
     if no numbers were parsed. */
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
  if ( nTot == 0 ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
  }

  /* Allocate **vector. */
  LALI2CreateVector( stat->statusPtr, vector, nTot );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*vector)->data, and set
     (*vector)->length. */
  here = &head;
  data = (*vector)->data;
  while ( here ) {
    UINT4 j = here->size;
    INT2 *hereData = here->buf.i2;
    while ( j-- )
      *(data++) = *(hereData++);
    here = here->next;
  }

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALI4ReadVector( LALStatus *stat, INT4Vector **vector, FILE *stream )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  const CHAR *format;      /* the input format specifier */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty;    /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *token;             /* pointer to a number to be converted */
  INT4 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALI4ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Determine the format for sscanf(). */
  INTFORMAT( format, 4 );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Scan for comment delimiters. */
  if ( line->data[0] == '#' )
    line->data[0] = '\0';
  else
    for ( i = 0; i < line->length; i++ )
      if ( line->data[i] == '%' )
	line->data[i] = '\0';

  /* Read into the first buffer. */
  token = strtok( line->data, whitespace );
  for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ )
    if ( token )
    {
      if ( sscanf( token, format, head.buf.i4 + i ) != 1 )
	more = 0;
      else {
	nTot++;
	head.size++;
	token = strtok( NULL, whitespace );
      }
    }
    else
      more = 0;

  /* Read into remaining buffers. */
  here = &head;
  while ( more ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALCalloc( 1, sizeof(BufferList) );
    if ( !here->next ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here = here->next;

    /* Read into the next buffer, and see if more needs to be read. */
    for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ )
      if ( token )
      {
	if ( sscanf( token, format, here->buf.i4 + i ) != 1 )
	  more = 0;
	else {
	  nTot++;
	  here->size++;
	  token = strtok( NULL, whitespace );
	}
      }
      else
	more = 0;
  }

  /* Line has been parsed into numbers, so free it and return an error
     if no numbers were parsed. */
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
  if ( nTot == 0 ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
  }

  /* Allocate **vector. */
  LALI4CreateVector( stat->statusPtr, vector, nTot );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*vector)->data, and set
     (*vector)->length. */
  here = &head;
  data = (*vector)->data;
  while ( here ) {
    UINT4 j = here->size;
    INT4 *hereData = here->buf.i4;
    while ( j-- )
      *(data++) = *(hereData++);
    here = here->next;
  }

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALI8ReadVector( LALStatus *stat, INT8Vector **vector, FILE *stream )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  const CHAR *format;      /* the input format specifier */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty;    /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *token;             /* pointer to a number to be converted */
  INT8 *data;              /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALI8ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Determine the format for sscanf(). */
  INTFORMAT( format, 8 );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Scan for comment delimiters. */
  if ( line->data[0] == '#' )
    line->data[0] = '\0';
  else
    for ( i = 0; i < line->length; i++ )
      if ( line->data[i] == '%' )
	line->data[i] = '\0';

  /* Read into the first buffer. */
  token = strtok( line->data, whitespace );
  for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ )
    if ( token )
    {
      if ( sscanf( token, format, head.buf.i8 + i ) != 1 )
	more = 0;
      else {
	nTot++;
	head.size++;
	token = strtok( NULL, whitespace );
      }
    }
    else
      more = 0;

  /* Read into remaining buffers. */
  here = &head;
  while ( more ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALCalloc( 1, sizeof(BufferList) );
    if ( !here->next ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here = here->next;

    /* Read into the next buffer, and see if more needs to be read. */
    for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ )
      if ( token )
      {
	if ( sscanf( token, format, here->buf.i8 + i ) != 1 )
	  more = 0;
	else {
	  nTot++;
	  here->size++;
	  token = strtok( NULL, whitespace );
	}
      }
      else
	more = 0;
  }

  /* Line has been parsed into numbers, so free it and return an error
     if no numbers were parsed. */
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
  if ( nTot == 0 ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
  }

  /* Allocate **vector. */
  LALI8CreateVector( stat->statusPtr, vector, nTot );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*vector)->data, and set
     (*vector)->length. */
  here = &head;
  data = (*vector)->data;
  while ( here ) {
    UINT4 j = here->size;
    INT8 *hereData = here->buf.i8;
    while ( j-- )
      *(data++) = *(hereData++);
    here = here->next;
  }

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALU2ReadVector( LALStatus *stat, UINT2Vector **vector, FILE *stream )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  const CHAR *format;      /* the input format specifier */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty;    /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *token;             /* pointer to a number to be converted */
  UINT2 *data;             /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALU2ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Determine the format for sscanf(). */
  INTFORMAT( format, 2 );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Scan for comment delimiters. */
  if ( line->data[0] == '#' )
    line->data[0] = '\0';
  else
    for ( i = 0; i < line->length; i++ )
      if ( line->data[i] == '%' )
	line->data[i] = '\0';

  /* Read into the first buffer. */
  token = strtok( line->data, whitespace );
  for ( i = 0; more && ( i < BUFFSIZE/2 ); i++ )
    if ( token )
    {
      if ( sscanf( token, format, head.buf.u2 + i ) != 1 )
	more = 0;
      else {
	nTot++;
	head.size++;
	token = strtok( NULL, whitespace );
      }
    }
    else
      more = 0;

  /* Read into remaining buffers. */
  here = &head;
  while ( more ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALCalloc( 1, sizeof(BufferList) );
    if ( !here->next ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here = here->next;

    /* Read into the next buffer, and see if more needs to be read. */
    for ( i = 0; more && ( i < BUFFSIZE/2 ); i++ )
      if ( token )
      {
	if ( sscanf( token, format, here->buf.u2 + i ) != 1 )
	  more = 0;
	else {
	  nTot++;
	  here->size++;
	  token = strtok( NULL, whitespace );
	}
      }
      else
	more = 0;
  }

  /* Line has been parsed into numbers, so free it and return an error
     if no numbers were parsed. */
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
  if ( nTot == 0 ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
  }

  /* Allocate **vector. */
  LALU2CreateVector( stat->statusPtr, vector, nTot );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*vector)->data, and set
     (*vector)->length. */
  here = &head;
  data = (*vector)->data;
  while ( here ) {
    UINT4 j = here->size;
    UINT2 *hereData = here->buf.u2;
    while ( j-- )
      *(data++) = *(hereData++);
    here = here->next;
  }

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALU4ReadVector( LALStatus *stat, UINT4Vector **vector, FILE *stream )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  const CHAR *format;      /* the input format specifier */
  CHARVector *line = NULL; /* the line of text stored as a CHARVector */
  BufferList head = empty;    /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *token;             /* pointer to a number to be converted */
  UINT4 *data;             /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALU4ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Determine the format for sscanf(). */
  INTFORMAT( format, 4 );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Scan for comment delimiters. */
  if ( line->data[0] == '#' )
    line->data[0] = '\0';
  else
    for ( i = 0; i < line->length; i++ )
      if ( line->data[i] == '%' )
	line->data[i] = '\0';

  /* Read into the first buffer. */
  token = strtok( line->data, whitespace );
  for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ )
    if ( token )
    {
      if ( sscanf( token, format, head.buf.u4 + i ) != 1 )
	more = 0;
      else {
	nTot++;
	head.size++;
	token = strtok( NULL, whitespace );
      }
    }
    else
      more = 0;

  /* Read into remaining buffers. */
  here = &head;
  while ( more ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALCalloc( 1, sizeof(BufferList) );
    if ( !here->next ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here = here->next;

    /* Read into the next buffer, and see if more needs to be read. */
    for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ )
      if ( token )
      {
	if ( sscanf( token, format, here->buf.u4 + i ) != 1 )
	  more = 0;
	else {
	  nTot++;
	  here->size++;
	  token = strtok( NULL, whitespace );
	}
      }
      else
	more = 0;
  }

  /* Line has been parsed into numbers, so free it and return an error
     if no numbers were parsed. */
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
  if ( nTot == 0 ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
  }

  /* Allocate **vector. */
  LALU4CreateVector( stat->statusPtr, vector, nTot );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*vector)->data, and set
     (*vector)->length. */
  here = &head;
  data = (*vector)->data;
  while ( here ) {
    UINT4 j = here->size;
    UINT4 *hereData = here->buf.u4;
    while ( j-- )
      *(data++) = *(hereData++);
    here = here->next;
  }

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALU8ReadVector( LALStatus *stat, UINT8Vector **vector, FILE *stream )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  const CHAR *format;      /* the input format specifier */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty;    /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *token;             /* pointer to a number to be converted */
  UINT8 *data;             /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALU8ReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Determine the format for sscanf(). */
  INTFORMAT( format, 8 );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Scan for comment delimiters. */
  if ( line->data[0] == '#' )
    line->data[0] = '\0';
  else
    for ( i = 0; i < line->length; i++ )
      if ( line->data[i] == '%' )
	line->data[i] = '\0';

  /* Read into the first buffer. */
  token = strtok( line->data, whitespace );
  for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ )
    if ( token )
    {
      if ( sscanf( token, format, head.buf.u8 + i ) != 1 )
	more = 0;
      else {
	nTot++;
	head.size++;
	token = strtok( NULL, whitespace );
      }
    }
    else
      more = 0;

  /* Read into remaining buffers. */
  here = &head;
  while ( more ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALCalloc( 1, sizeof(BufferList) );
    if ( !here->next ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here = here->next;

    /* Read into the next buffer, and see if more needs to be read. */
    for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ )
      if ( token )
      {
	if ( sscanf( token, format, here->buf.u8 + i ) != 1 )
	  more = 0;
	else {
	  nTot++;
	  here->size++;
	  token = strtok( NULL, whitespace );
	}
      }
      else
	more = 0;
  }

  /* Line has been parsed into numbers, so free it and return an error
     if no numbers were parsed. */
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
  if ( nTot == 0 ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
  }

  /* Allocate **vector. */
  LALU8CreateVector( stat->statusPtr, vector, nTot );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*vector)->data, and set
     (*vector)->length. */
  here = &head;
  data = (*vector)->data;
  while ( here ) {
    UINT4 j = here->size;
    UINT8 *hereData = here->buf.u8;
    while ( j-- )
      *(data++) = *(hereData++);
    here = here->next;
  }

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALSReadVector( LALStatus *stat, REAL4Vector **vector, FILE *stream )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  const CHAR *format;      /* the input format specifier */
  CHARVector *line = NULL; /* the line of text stored as a CHARVector */
  BufferList head = empty;    /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *token;             /* pointer to a number to be converted */
  REAL4 *data;             /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALSReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Determine the format for sscanf(). */
  REALFORMAT( format, 4 );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Scan for comment delimiters. */
  if ( line->data[0] == '#' )
    line->data[0] = '\0';
  else
    for ( i = 0; i < line->length; i++ )
      if ( line->data[i] == '%' )
	line->data[i] = '\0';

  /* Read into the first buffer. */
  token = strtok( line->data, whitespace );
  for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ )
    if ( token )
    {
      if ( sscanf( token, format, head.buf.s + i ) != 1 )
	more = 0;
      else {
	nTot++;
	head.size++;
	token = strtok( NULL, whitespace );
      }
    }
    else
      more = 0;

  /* Read into remaining buffers. */
  here = &head;
  while ( more ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALCalloc( 1, sizeof(BufferList) );
    if ( !here->next ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here = here->next;

    /* Read into the next buffer, and see if more needs to be read. */
    for ( i = 0; more && ( i < BUFFSIZE/4 ); i++ )
      if ( token )
      {
	if ( sscanf( token, format, here->buf.s + i ) != 1 )
	  more = 0;
	else {
	  nTot++;
	  here->size++;
	  token = strtok( NULL, whitespace );
	}
      }
      else
	more = 0;
  }

  /* Line has been parsed into numbers, so free it and return an error
     if no numbers were parsed. */
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
  if ( nTot == 0 ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
  }

  /* Allocate **vector. */
  LALSCreateVector( stat->statusPtr, vector, nTot );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*vector)->data, and set
     (*vector)->length. */
  here = &head;
  data = (*vector)->data;
  while ( here ) {
    UINT4 j = here->size;
    REAL4 *hereData = here->buf.s;
    while ( j-- )
      *(data++) = *(hereData++);
    here = here->next;
  }

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="StreamVectorInputCP"> */
void
LALDReadVector( LALStatus *stat, REAL8Vector **vector, FILE *stream )
{ /* </lalVerbatim> */
  UINT4 i;                 /* an index */
  const CHAR *format;      /* the input format specifier */
  CHARVector *line = NULL; /* a line of text stored as a CHARVector */
  BufferList head = empty;    /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHAR *token;             /* pointer to a number to be converted */
  REAL8 *data;             /* pointer to converted data */
  BOOLEAN more = 1;        /* whether or not to read more numbers */
  UINT4 nTot = 0;          /* total number of numbers read */

  INITSTATUS( stat, "LALSReadVector", STREAMVECTORINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( vector, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !*vector, stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  /* Determine the format for sscanf(). */
  REALFORMAT( format, 8 );

  /* Read the line of text as a CHARVector. */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );

  /* Scan for comment delimiters. */
  if ( line->data[0] == '#' )
    line->data[0] = '\0';
  else
    for ( i = 0; i < line->length; i++ )
      if ( line->data[i] == '%' )
	line->data[i] = '\0';

  /* Read into the first buffer. */
  token = strtok( line->data, whitespace );
  for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ )
    if ( token )
    {
      if ( sscanf( token, format, head.buf.d + i ) != 1 )
	more = 0;
      else {
	nTot++;
	head.size++;
	token = strtok( NULL, whitespace );
      }
    }
    else
      more = 0;

  /* Read into remaining buffers. */
  here = &head;
  while ( more ) {

    /* Allocate next buffer. */
    here->next = (BufferList *)LALCalloc( 1, sizeof(BufferList) );
    if ( !here->next ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here = here->next;

    /* Read into the next buffer, and see if more needs to be read. */
    for ( i = 0; more && ( i < BUFFSIZE/8 ); i++ )
      if ( token )
      {
	if ( sscanf( token, format, here->buf.d + i ) != 1 )
	  more = 0;
	else {
	  nTot++;
	  here->size++;
	  token = strtok( NULL, whitespace );
	}
      }
      else
	more = 0;
  }

  /* Line has been parsed into numbers, so free it and return an error
     if no numbers were parsed. */
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
  if ( nTot == 0 ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ELEN, STREAMINPUTH_MSGELEN );
  }

  /* Allocate **vector. */
  LALDCreateVector( stat->statusPtr, vector, nTot );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Copy buffer list into (*vector)->data, and set
     (*vector)->length. */
  here = &head;
  data = (*vector)->data;
  while ( here ) {
    UINT4 j = here->size;
    REAL8 *hereData = here->buf.d;
    while ( j-- )
      *(data++) = *(hereData++);
    here = here->next;
  }

  /* Free buffer list and exit. */
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
