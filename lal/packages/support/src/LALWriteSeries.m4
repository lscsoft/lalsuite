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
define(`STYPE',`format(`%sTimeSeries',TYPE)')dnl
define(`VTYPE',`format(`%sTimeVectorSeries',TYPE)')dnl
define(`ATYPE',`format(`%sTimeArraySeries',TYPE)')dnl
define(`FTYPE',`format(`%sFrequencySeries',TYPE)')dnl
define(`SFUNC',`format(`LAL%sWriteTSeries',TYPECODE)')dnl
define(`VFUNC',`format(`LAL%sWriteTVectorSeries',TYPECODE)')dnl
define(`AFUNC',`format(`LAL%sWriteTArraySeries',TYPECODE)')dnl
define(`FFUNC',`format(`LAL%sWriteFSeries',TYPECODE)')dnl
dnl
/* <lalVerbatim file="StreamSeriesInputCP"> */
void
SFUNC ( LALStatus *stat, FILE *stream, STYPE *series )
{ /* </lalVerbatim> */
  UINT4 length; /* length of data sequence */
  TYPE *data;   /* pointer to data in sequence */

  INITSTATUS( stat, "SFUNC", STREAMSERIESOUTPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series->data, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series->data->data, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );


  /*******************************************************************
   * PRINT METADATA HEADER                                           *
   *******************************************************************/

  /* Print the datatype. */
  if ( fprintf( stream, "# datatype = STYPE\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the name. */
  if ( fprintf( stream, "# name = " ) < 0 ||
       LALWriteLiteral( stream, series->name ) ||
       fprintf( stream, "\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the epoch, deltaT, and f0. */
  if ( fprintf( stream, "# epoch = %" LAL_INT8_FORMAT "\n",
                XLALGPSToINT8NS( &series->epoch ) ) < 0 ||
       fprintf( stream, "# deltaT = %.16e\n", series->deltaT ) < 0 ||
       fprintf( stream, "# f0 = %.16e\n", series->f0 ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the sample units. */
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
      LALUnitAsString( stat->statusPtr, unit, &(series->sampleUnits) );
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
    TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
    if ( code == EOF || fprintf( stream, "\"\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }

  /* Print the length. */
  if ( fprintf( stream, "# length = %" LAL_UINT4_FORMAT "\n",
		series->data->length ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }


  /*******************************************************************
   * PRINT DATA                                                      *
   *******************************************************************/

  length = series->data->length;
  data = series->data->data;
  while ( length-- ) {
#if COMPLEX
    if ( fprintf( stream, FMT " " FMT "\n", data->re, data->im ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
    data++;
#else
    if ( fprintf( stream, FMT "\n", *(data++) ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
#endif
  }
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}

/* <lalVerbatim file="StreamSeriesInputCP"> */
void
VFUNC ( LALStatus *stat, FILE *stream, VTYPE *series )
{ /* </lalVerbatim> */
  UINT4 length;       /* length of data sequence */
  UINT4 vectorLength; /* length of each element in sequence */
  TYPE *data; /* pointer to data in sequence */

  INITSTATUS( stat, "VFUNC", STREAMSERIESOUTPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series->data, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series->data->data, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );


  /*******************************************************************
   * PRINT METADATA HEADER                                           *
   *******************************************************************/

  /* Print the datatype. */
  if ( fprintf( stream, "# datatype = VTYPE\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the name. */
  if ( fprintf( stream, "# name = " ) < 0 ||
       LALWriteLiteral( stream, series->name ) ||
       fprintf( stream, "\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the epoch, deltaT, and f0. */
  if ( fprintf( stream, "# epoch = %" LAL_INT8_FORMAT "\n",
                XLALGPSToINT8NS( &series->epoch ) ) < 0 ||
       fprintf( stream, "# deltaT = %.16e\n", series->deltaT ) < 0 ||
       fprintf( stream, "# f0 = %.16e\n", series->f0 ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the sample units. */
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
      LALUnitAsString( stat->statusPtr, unit, &(series->sampleUnits) );
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
    TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
    if ( code == EOF || fprintf( stream, "\"\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }

  /* Print the length and vectorLength. */
  if ( fprintf( stream, "# length = %" LAL_UINT4_FORMAT "\n",
		series->data->length ) < 0 ||
       fprintf( stream, "# vectorLength = %" LAL_UINT4_FORMAT "\n",
		series->data->vectorLength ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }


  /*******************************************************************
   * PRINT DATA                                                      *
   *******************************************************************/

  length = series->data->length;
  vectorLength = series->data->vectorLength;
  data = series->data->data;
  while ( length-- ) {
    UINT4 n = vectorLength - 1;
#if COMPLEX
    if ( fprintf( stream, FMT " " FMT, data->re, data->im ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
    data++;
#else
    if ( fprintf( stream, FMT, *(data++) ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
#endif
    while ( n-- ) {
#if COMPLEX
      if ( fprintf( stream, " " FMT " " FMT, data->re, data->im ) < 0 ) {
	ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
      }
      data++;
#else
      if ( fprintf( stream, " " FMT, *(data++) ) < 0 ) {
	ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
      }
#endif
    }
    if ( fprintf( stream, "\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}

/* <lalVerbatim file="StreamSeriesInputCP"> */
void
AFUNC ( LALStatus *stat, FILE *stream, ATYPE *series )
{ /* </lalVerbatim> */
  UINT4 length;   /* length of data sequence */
  UINT4 arrayDim; /* length of each element in sequence */
  UINT4 *dimData; /* pointer to dimLength data */
  TYPE *data; /* pointer to data in sequence */

  INITSTATUS( stat, "AFUNC", STREAMSERIESOUTPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series->data, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series->data->data, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series->data->dimLength, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series->data->dimLength->data, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );


  /*******************************************************************
   * PRINT METADATA HEADER                                           *
   *******************************************************************/

  /* Print the datatype. */
  if ( fprintf( stream, "# datatype = ATYPE\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the name. */
  if ( fprintf( stream, "# name = " ) < 0 ||
       LALWriteLiteral( stream, series->name ) ||
       fprintf( stream, "\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the epoch, deltaT, and f0. */
  if ( fprintf( stream, "# epoch = %" LAL_INT8_FORMAT "\n",
                XLALGPSToINT8NS( &series->epoch ) ) < 0 ||
       fprintf( stream, "# deltaT = %.16e\n", series->deltaT ) < 0 ||
       fprintf( stream, "# f0 = %.16e\n", series->f0 ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the sample units. */
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
      LALUnitAsString( stat->statusPtr, unit, &(series->sampleUnits) );
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
    TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
    if ( code == EOF || fprintf( stream, "\"\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }

  /* Print the length, dimLength, and arrayDim. */
  if ( fprintf( stream, "# length = %" LAL_UINT4_FORMAT "\n",
		series->data->length ) < 0 ||
       fprintf( stream, "# arrayDim = %" LAL_UINT4_FORMAT "\n",
		series->data->arrayDim ) < 0 ||
       fprintf( stream, "# dimLength =" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }
  arrayDim = series->data->dimLength->length;
  dimData = series->data->dimLength->data;
  while ( arrayDim-- )
    if ( fprintf( stream, " %" LAL_UINT4_FORMAT, *(dimData++) ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  if ( fprintf( stream, "\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }


  /*******************************************************************
   * PRINT DATA                                                      *
   *******************************************************************/

  length = series->data->length;
  arrayDim = series->data->arrayDim;
  data = series->data->data;
  while ( length-- ) {
    UINT4 n = arrayDim - 1;
#if COMPLEX
    if ( fprintf( stream, FMT " " FMT, data->re, data->im ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
    data++;
#else
    if ( fprintf( stream, FMT, *(data++) ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
#endif
    while ( n-- ) {
#if COMPLEX
      if ( fprintf( stream, " " FMT " " FMT, data->re, data->im ) < 0 ) {
	ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
      }
      data++;
#else
      if ( fprintf( stream, " " FMT, *(data++) ) < 0 ) {
	ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
      }
#endif
    }
    if ( fprintf( stream, "\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}

/* <lalVerbatim file="StreamSeriesInputCP"> */
void
FFUNC ( LALStatus *stat, FILE *stream, FTYPE *series )
{ /* </lalVerbatim> */
  UINT4 length; /* length of data sequence */
  TYPE *data;   /* pointer to data in sequence */

  INITSTATUS( stat, "FFUNC", STREAMSERIESOUTPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series->data, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );
  ASSERT( series->data->data, stat, STREAMOUTPUTH_ENUL, STREAMOUTPUTH_MSGENUL );


  /*******************************************************************
   * PRINT METADATA HEADER                                           *
   *******************************************************************/

  /* Print the datatype. */
  if ( fprintf( stream, "# datatype = FTYPE\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the name. */
  if ( fprintf( stream, "# name = " ) < 0 ||
       LALWriteLiteral( stream, series->name ) ||
       fprintf( stream, "\n" ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the epoch, deltaF, and f0. */
  if ( fprintf( stream, "# epoch = %" LAL_INT8_FORMAT "\n",
                XLALGPSToINT8NS( &series->epoch ) ) < 0 ||
       fprintf( stream, "# deltaF = %.16e\n", series->deltaF ) < 0 ||
       fprintf( stream, "# f0 = %.16e\n", series->f0 ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }

  /* Print the sample units. */
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
      LALUnitAsString( stat->statusPtr, unit, &(series->sampleUnits) );
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
    TRY( LALCHARDestroyVector( stat->statusPtr, &unit ), stat );
    if ( code == EOF || fprintf( stream, "\"\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }

  /* Print the length. */
  if ( fprintf( stream, "# length = %" LAL_UINT4_FORMAT "\n",
		series->data->length ) < 0 ) {
    ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
  }


  /*******************************************************************
   * PRINT DATA                                                      *
   *******************************************************************/

  length = series->data->length;
  data = series->data->data;
  while ( length-- ) {
#if COMPLEX
    if ( fprintf( stream, FMT " " FMT "\n", data->re, data->im ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
    data++;
#else
    if ( fprintf( stream, FMT "\n", *(data++) ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
#endif
  }
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
