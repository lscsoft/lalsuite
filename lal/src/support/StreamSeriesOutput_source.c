#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,TimeSeries)
#define VTYPE CONCAT2(TYPE,TimeVectorSeries)
#define ATYPE CONCAT2(TYPE,TimeArraySeries)
#define FTYPE CONCAT2(TYPE,FrequencySeries)
#define SFUNC CONCAT3(LAL,TYPECODE,WriteTSeries)
#define VFUNC CONCAT3(LAL,TYPECODE,WriteTVectorSeries)
#define AFUNC CONCAT3(LAL,TYPECODE,WriteTArraySeries)
#define FFUNC CONCAT3(LAL,TYPECODE,WriteFSeries)


void
SFUNC ( LALStatus *stat, FILE *stream, STYPE *series )
{ 
  UINT4 length; /* length of data sequence */
  TYPE *data;   /* pointer to data in sequence */

  INITSTATUS(stat);
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
    CHAR unitString[LALUnitTextSize];

    /* First, use XLALUnitAsString() to generate unit string. */
    if ( XLALUnitAsString( unitString, LALUnitTextSize, &(series->sampleUnits) ) == NULL ) {
      ABORTXLAL(stat);
    }

    /* Write the resulting unit string enclosed in quotes. */
    if ( fprintf( stream, "# sampleUnits = \"%s\"\n", unitString ) < 0 ) {
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
    if ( fprintf( stream, FMT " " FMT "\n", creal(*data), cimag(*data) ) < 0 ) {
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


void
VFUNC ( LALStatus *stat, FILE *stream, VTYPE *series )
{ 
  UINT4 length;       /* length of data sequence */
  UINT4 vectorLength; /* length of each element in sequence */
  TYPE *data; /* pointer to data in sequence */

  INITSTATUS(stat);
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
    CHAR unitString[LALUnitTextSize];

    /* First, use XLALUnitAsString() to generate unit string. */
    if ( XLALUnitAsString( unitString, LALUnitTextSize, &(series->sampleUnits) ) == NULL ) {
      ABORTXLAL(stat);
    }

    /* Write the resulting unit string enclosed in quotes. */
    if ( fprintf( stream, "# sampleUnits = \"%s\"\n", unitString ) < 0 ) {
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
    if ( fprintf( stream, FMT " " FMT, creal(*data), cimag(*data) ) < 0 ) {
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
      if ( fprintf( stream, " " FMT " " FMT, creal(*data), cimag(*data) ) < 0 ) {
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


void
AFUNC ( LALStatus *stat, FILE *stream, ATYPE *series )
{ 
  UINT4 length;   /* length of data sequence */
  UINT4 arrayDim; /* length of each element in sequence */
  UINT4 *dimData; /* pointer to dimLength data */
  TYPE *data; /* pointer to data in sequence */

  INITSTATUS(stat);
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
    CHAR unitString[LALUnitTextSize];

    /* First, use XLALUnitAsString() to generate unit string. */
    if ( XLALUnitAsString( unitString, LALUnitTextSize, &(series->sampleUnits) ) == NULL ) {
      ABORTXLAL(stat);
    }

    /* Write the resulting unit string enclosed in quotes. */
    if ( fprintf( stream, "# sampleUnits = \"%s\"\n", unitString ) < 0 ) {
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
    if ( fprintf( stream, FMT " " FMT, creal(*data), cimag(*data) ) < 0 ) {
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
      if ( fprintf( stream, " " FMT " " FMT, creal(*data), cimag(*data) ) < 0 ) {
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


void
FFUNC ( LALStatus *stat, FILE *stream, FTYPE *series )
{ 
  UINT4 length; /* length of data sequence */
  TYPE *data;   /* pointer to data in sequence */

  INITSTATUS(stat);
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
    CHAR unitString[LALUnitTextSize];

    /* First, use XLALUnitAsString() to generate unit string. */
    if ( XLALUnitAsString( unitString, LALUnitTextSize, &(series->sampleUnits) ) == NULL ) {
      ABORTXLAL(stat);
    }

    /* Write the resulting unit string enclosed in quotes. */
    if ( fprintf( stream, "# sampleUnits = \"%s\"\n", unitString ) < 0 ) {
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
    if ( fprintf( stream, FMT " " FMT "\n", creal(*data), cimag(*data) ) < 0 ) {
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
