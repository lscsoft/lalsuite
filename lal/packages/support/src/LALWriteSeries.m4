/* m4_dnl $Id$
m4_ifelse(TYPECODE,`Z',`m4_define(`TYPE',`COMPLEX16')m4_define(`FMT',`"% .17e"')')m4_dnl
m4_ifelse(TYPECODE,`C',`m4_define(`TYPE',`COMPLEX8')m4_define(`FMT',`"% .9e"')')m4_dnl
m4_ifelse(TYPECODE,`D',`m4_define(`TYPE',`REAL8')m4_define(`FMT',`"% .17e"')')m4_dnl
m4_ifelse(TYPECODE,`S',`m4_define(`TYPE',`REAL4')m4_define(`FMT',`"% .9e"')')m4_dnl
m4_ifelse(TYPECODE,`I2',`m4_define(`TYPE',`INT2')m4_define(`FMT',`"% " LAL_INT2_FORMAT')')m4_dnl
m4_ifelse(TYPECODE,`I4',`m4_define(`TYPE',`INT4')m4_define(`FMT',`"% " LAL_INT4_FORMAT')')m4_dnl
m4_ifelse(TYPECODE,`I8',`m4_define(`TYPE',`INT8')m4_define(`FMT',`"% " LAL_INT8_FORMAT')')m4_dnl
m4_ifelse(TYPECODE,`U2',`m4_define(`TYPE',`UINT2')m4_define(`FMT',`"%" LAL_UINT2_FORMAT')')m4_dnl
m4_ifelse(TYPECODE,`U4',`m4_define(`TYPE',`UINT4')m4_define(`FMT',`"%" LAL_UINT4_FORMAT')')m4_dnl
m4_ifelse(TYPECODE,`U8',`m4_define(`TYPE',`UINT8')m4_define(`FMT',`"%" LAL_UINT8_FORMAT')')m4_dnl
m4_ifelse(TYPECODE,`Z',`m4_define(`COMPLEX',`1')')m4_dnl
m4_ifelse(TYPECODE,`C',`m4_define(`COMPLEX',`1')')m4_dnl
m4_define(`STYPE',`m4_format(`%sTimeSeries',TYPE)')m4_dnl
m4_define(`VTYPE',`m4_format(`%sTimeVectorSeries',TYPE)')m4_dnl
m4_define(`ATYPE',`m4_format(`%sTimeArraySeries',TYPE)')m4_dnl
m4_define(`FTYPE',`m4_format(`%sFrequencySeries',TYPE)')m4_dnl
m4_define(`SFUNC',`m4_format(`LAL%sWriteTSeries',TYPECODE)')m4_dnl
m4_define(`VFUNC',`m4_format(`LAL%sWriteTVectorSeries',TYPECODE)')m4_dnl
m4_define(`AFUNC',`m4_format(`LAL%sWriteTArraySeries',TYPECODE)')m4_dnl
m4_define(`FFUNC',`m4_format(`LAL%sWriteFSeries',TYPECODE)')m4_dnl
m4_dnl
<lalVerbatim file="StreamSeriesInputCP"> */
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

  /* Print the name. */
  {
    CHAR *cData = series->name; /* pointer to name data */
    CHAR c;                     /* value of *cData */
    int code = 0;               /* return code from putc() */
    length = LALNameLength;
    if ( fprintf( stream, "# name = \"" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
    while ( code != EOF && ( c = *(cData++) ) != '"' && c != '\n'
	    && c != '\0' && length-- )
      code = putc( (int)c, stream );
    if ( code == EOF || fprintf( stream, "\"\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }

  /* Print the epoch, deltaT, and f0. */
  if ( fprintf( stream, "# epoch = %" LAL_INT8_FORMAT "\n",
		(INT8)( series->epoch.gpsSeconds )*1000000000LL +
		(INT8)( series->epoch.gpsNanoSeconds ) ) < 0 ||
       fprintf( stream, "# deltaT = %.15e\n", series->deltaT ) < 0 ||
       fprintf( stream, "# f0 = %.15e\n", series->f0 ) < 0 ) {
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
`#'if COMPLEX
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

  /* Print the name. */
  {
    CHAR *cData = series->name; /* pointer to name data */
    CHAR c;                     /* value of *cData */
    int code = 0;               /* return code from putc() */
    length = LALNameLength;
    if ( fprintf( stream, "# name = \"" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
    while ( code != EOF && ( c = *(cData++) ) != '"' && c != '\n'
	    && c != '\0' && length-- )
      code = putc( (int)c, stream );
    if ( code == EOF || fprintf( stream, "\"\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }

  /* Print the epoch, deltaT, and f0. */
  if ( fprintf( stream, "# epoch = %" LAL_INT8_FORMAT "\n",
		(INT8)( series->epoch.gpsSeconds )*1000000000LL +
		(INT8)( series->epoch.gpsNanoSeconds ) ) < 0 ||
       fprintf( stream, "# deltaT = %.15e\n", series->deltaT ) < 0 ||
       fprintf( stream, "# f0 = %.15e\n", series->f0 ) < 0 ) {
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
`#'if COMPLEX
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
`#'if COMPLEX
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

  /* Print the name. */
  {
    CHAR *cData = series->name; /* pointer to name data */
    CHAR c;                     /* value of *cData */
    int code = 0;               /* return code from putc() */
    length = LALNameLength;
    if ( fprintf( stream, "# name = \"" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
    while ( code != EOF && ( c = *(cData++) ) != '"' && c != '\n'
	    && c != '\0' && length-- )
      code = putc( (int)c, stream );
    if ( code == EOF || fprintf( stream, "\"\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }

  /* Print the epoch, deltaT, and f0. */
  if ( fprintf( stream, "# epoch = %" LAL_INT8_FORMAT "\n",
		(INT8)( series->epoch.gpsSeconds )*1000000000LL +
		(INT8)( series->epoch.gpsNanoSeconds ) ) < 0 ||
       fprintf( stream, "# deltaT = %.15e\n", series->deltaT ) < 0 ||
       fprintf( stream, "# f0 = %.15e\n", series->f0 ) < 0 ) {
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
`#'if COMPLEX
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
`#'if COMPLEX
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

  /* Print the name. */
  {
    CHAR *cData = series->name; /* pointer to name data */
    CHAR c;                     /* value of *cData */
    int code = 0;               /* return code from putc() */
    length = LALNameLength;
    if ( fprintf( stream, "# name = \"" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
    while ( code != EOF && ( c = *(cData++) ) != '"' && c != '\n'
	    && c != '\0' && length-- )
      code = putc( (int)c, stream );
    if ( code == EOF || fprintf( stream, "\"\n" ) < 0 ) {
      ABORT( stat, STREAMOUTPUTH_EPRN, STREAMOUTPUTH_MSGEPRN );
    }
  }

  /* Print the epoch, deltaF, and f0. */
  if ( fprintf( stream, "# epoch = %" LAL_INT8_FORMAT "\n",
		(INT8)( series->epoch.gpsSeconds )*1000000000LL +
		(INT8)( series->epoch.gpsNanoSeconds ) ) < 0 ||
       fprintf( stream, "# f0 = %.15e\n", series->f0 ) < 0 ||
       fprintf( stream, "# deltaF = %.15e\n", series->deltaF ) < 0 ) {
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
`#'if COMPLEX
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
