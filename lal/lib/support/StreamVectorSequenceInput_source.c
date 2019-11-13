#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define VTYPE CONCAT2(TYPE,VectorSequence)
#define VTYPECODE CONCAT2(TYPECODE,V)
#define FUNC CONCAT3(LAL,TYPECODE,ReadVectorSequence)
#define VFUNC CONCAT3(LAL,TYPECODE,ReadVector)
#define CREATEFUNC CONCAT3(LAL,TYPECODE,CreateVectorSequence)
#define DESTROYFUNC CONCAT3(LAL,TYPECODE,DestroyVector)

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


void
FUNC ( LALStatus  *stat, VTYPE **sequence, FILE *stream )
{ 
  VectorList head = empty; /* head of linked list of vectors */
  VectorList *here;     /* pointer to current position in list */
  TYPE *data;           /* pointer to vector data */
  UINT4 nRows, nCols;   /* number and length of lines */
  CreateVectorSequenceIn in; /* parameters for creating sequence */

  INITSTATUS(stat);
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
