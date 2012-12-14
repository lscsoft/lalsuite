#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define GTYPE CONCAT2(TYPE,Grid)
#define FUNC CONCAT3(LAL,TYPECODE,WriteGrid)

void
FUNC ( LALStatus *stat, FILE *stream, GTYPE *grid )
{
  UINT4 i, j, length; /* indecies, and line/string length */
  UINT4 pLength, np;  /* length and number of data ``paragraphs'' */
  TYPE *data;         /* pointer to grid->data->data */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( grid, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( grid->dimUnits, stat, STREAMOUTPUTH_ENUL,
	  STREAMOUTPUTH_MSGENUL );
  ASSERT( grid->offset, stat, STREAMOUTPUTH_ENUL,
	  STREAMOUTPUTH_MSGENUL );
  ASSERT( grid->offset->data, stat, STREAMOUTPUTH_ENUL,
	  STREAMOUTPUTH_MSGENUL );
  ASSERT( grid->interval, stat, STREAMOUTPUTH_ENUL,
	  STREAMOUTPUTH_MSGENUL );
  ASSERT( grid->interval->data, stat, STREAMOUTPUTH_ENUL,
	  STREAMOUTPUTH_MSGENUL );
  ASSERT( grid->data, stat, STREAMOUTPUTH_ENUL,
	  STREAMOUTPUTH_MSGENUL );
  ASSERT( grid->data->data, stat, STREAMOUTPUTH_ENUL,
	  STREAMOUTPUTH_MSGENUL );
  ASSERT( grid->data->dimLength, stat, STREAMOUTPUTH_ENUL,
	  STREAMOUTPUTH_MSGENUL );
  ASSERT( grid->data->dimLength->data, stat, STREAMOUTPUTH_ENUL,
	  STREAMOUTPUTH_MSGENUL );

  /*******************************************************************
   * PRINT METADATA HEADER                                           *
   *******************************************************************/

  /* Print the datatype. */
  if ( fprintf( stream, "# datatype = " STRING(GTYPE) "\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the name. */
  if ( fprintf( stream, "# name = " ) < 0 ||
       LALWriteLiteral( stream, grid->name ) ||
       fprintf( stream, "\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the sample units and dimension units. */
  {
    CHARVector *unit = NULL; /* vector to store unit string */
    CHAR *cData;             /* pointer to unit data */
    CHAR c;                  /* value of *cData */
    int code = 0;            /* return code from putc() */
    BOOLEAN repeat;          /* whether a longer vector is needed */

    /* First, use LALUnitAsString() to generate unit string.  Repeat
       if a longer vector is required. */
    length = 64;
    TRY( LALCHARCreateVector( stat->statusPtr, &unit, length ), stat );
    memset( unit->data, 0, length*sizeof(CHAR) );
    do {
      LALUnitAsString( stat->statusPtr, unit, &(grid->sampleUnits) );
      if ( stat->statusPtr->statusCode == UNITSH_ESTRINGSIZE ) {
	repeat = 1;
	length *= 2;
	TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
	TRY( LALCHARCreateVector( stat->statusPtr, &unit, length ), stat );
	memset( unit->data, 0, length*sizeof(CHAR) );
#ifndef NDEBUG
	if ( lalDebugLevel & LALERROR )
	  LALPrintError( "\tCONTINUE: Dealt with preceding error\n" );
#endif
      } else {
	repeat = 0;
	BEGINFAIL( stat )
	  TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
	ENDFAIL( stat );
      }
    } while ( repeat );

    /* Write the resulting unit string enclosed in quotes. */
    cData = unit->data;
    if ( fprintf( stream, "# sampleUnits = \"" ) < 0 ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
    while ( code != EOF && ( c = *(cData++) ) != '\0' && length-- )
      code = putc( (int)c, stream );
    if ( code == EOF || fprintf( stream, "\"\n" ) < 0 ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }

    /* Repeat above for all the units in grid->dimUnits. */
    if ( fprintf( stream, "# dimUnits =" ) < 0 ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
    for ( i = 0; i < grid->offset->length; i++ ) {
      memset( unit->data, 0, length*sizeof(CHAR) );
      do {
	LALUnitAsString( stat->statusPtr, unit, &(grid->sampleUnits) );
	if ( stat->statusPtr->statusCode == UNITSH_ESTRINGSIZE ) {
	  repeat = 1;
	  length *= 2;
	  TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
	  TRY( LALCHARCreateVector( stat->statusPtr, &unit, length ), stat );
	  memset( unit->data, 0, length*sizeof(CHAR) );
#ifndef NDEBUG
	  if ( lalDebugLevel & LALERROR )
	    LALPrintError( "\tCONTINUE: Dealt with preceding error\n" );
#endif
	} else {
	  repeat = 0;
	  BEGINFAIL( stat )
	    TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
	  ENDFAIL( stat );
	}
      } while ( repeat );
      cData = unit->data;
      if ( fprintf( stream, " \"" ) < 0 ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
	ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
      }
      while ( code != EOF && ( c = *(cData++) ) != '\0' && length-- )
	code = putc( (int)c, stream );
      if ( code == EOF || fprintf( stream, "\"" ) < 0 ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
	ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
      }
    }
    TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
    if ( fprintf( stream, "\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }

  /* Print the offset, interval, and dimLength. */
  if ( fprintf( stream, "# offset =" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }
  for ( i = 0; i < grid->offset->length; i++ ) {
    if ( fprintf( stream, " %.16e", grid->offset->data[i] ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }
  if ( fprintf( stream, "\n# interval =" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }
  for ( i = 0; i < grid->interval->length; i++ ) {
    if ( fprintf( stream, " %.16e", grid->interval->data[i] ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }
  if ( fprintf( stream, "\n# dimLength =" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }
  for ( i = 0; i < grid->data->dimLength->length; i++ ) {
    if ( fprintf( stream, " %" LAL_UINT4_FORMAT ,
		  grid->data->dimLength->data[i] ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }
  if ( fprintf( stream, "\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /*******************************************************************
   * PRINT DATA                                                      *
   *******************************************************************/

  /* Compute number of data per line, per ``paragraph'' (grid point),
     and number of grid points.  length equals one *less* than the
     number of data per line. */
  pLength = np = 1;
  for ( i = 0; i < grid->offset->length; i++ )
    np *= grid->data->dimLength->data[i];
  for ( ; i < grid->data->dimLength->length - 1; i++ )
    pLength *= grid->data->dimLength->data[i];
  if ( i >= grid->offset->length )
    length = grid->data->dimLength->data[i] - 1;
  else
    length = 0;

  /* If each grid point is a single line, don't bother with
     ``paragraph'' breaks. */
  if ( pLength == 1 ) {
    pLength = np;
    np = 1;
  }
  data = grid->data->data;

  /* Print data. */
  while ( np-- ) {
    fprintf( stream, "\n" );
    for ( i = 0; i < pLength; i++ ) {
#if COMPLEX
      for ( j = 0; j < length; j++ ) {
	if ( fprintf( stream, FMT " " FMT " ", creal(*data), cimag(*data) ) < 0 ) {
	  ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
	}
	data++;
      }
      if ( fprintf( stream, FMT " " FMT "\n", creal(*data), cimag(*data) ) < 0 ) {
	ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
      }
      data++;
#else
      for ( j = 0; j < length; j++ ) {
	if ( fprintf( stream, FMT " ", *(data++) ) < 0 ) {
	  ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
	}
      }
      if ( fprintf( stream, FMT "\n", *(data++) ) < 0 ) {
	ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
      }
#endif
    }
  }
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
