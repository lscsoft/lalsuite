dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')define(`FMT',`"% .16e"')')dnl
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')define(`FMT',`"% .8e"')')dnl
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')define(`FMT',`"% .16e"')')dnl
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')define(`FMT',`"% .8e"')')dnl
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')define(`FMT',`"% " LAL_INT2_FORMAT')')dnl
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')define(`FMT',`"% " LAL_INT4_FORMAT')')dnl
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')define(`FMT',`"% " LAL_INT8_FORMAT')')dnl
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')define(`FMT',`"%" LAL_UINT2_FORMAT')')dnl
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')define(`FMT',`"%" LAL_UINT4_FORMAT')')dnl
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')define(`FMT',`"%" LAL_UINT8_FORMAT')')dnl
define(`COMPLEX',`0')dnl
ifelse(TYPECODE,`Z',`define(`COMPLEX',`1')')dnl
ifelse(TYPECODE,`C',`define(`COMPLEX',`1')')dnl
define(`GTYPE',`format(`%sGrid',TYPE)')dnl
define(`FUNC',`format(`LAL%sWriteGrid',TYPECODE)')dnl
dnl
void
FUNC ( LALStatus *stat, FILE *stream, GTYPE *grid )
{
  UINT4 i, j, length; /* indecies, and line/string length */
  UINT4 pLength, np;  /* length and number of data ``paragraphs'' */
  TYPE *data;         /* pointer to grid->data->data */

  INITSTATUS( stat, "FUNC", STREAMGRIDOUTPUTC );
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
  if ( fprintf( stream, "# datatype = GTYPE\n" ) < 0 ) {
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
	if ( fprintf( stream, FMT " " FMT " ", data->re, data->im ) < 0 ) {
	  ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
	}
	data++;
      }
      if ( fprintf( stream, FMT " " FMT "\n", data->re, data->im ) < 0 ) {
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
