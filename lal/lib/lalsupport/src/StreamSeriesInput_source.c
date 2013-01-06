#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,TimeSeries)
#define VTYPE CONCAT2(TYPE,TimeVectorSeries)
#define ATYPE CONCAT2(TYPE,TimeArraySeries)
#define FTYPE CONCAT2(TYPE,FrequencySeries)
#define SFUNC CONCAT3(LAL,TYPECODE,ReadTSeries)
#define VFUNC CONCAT3(LAL,TYPECODE,ReadTVectorSeries)
#define AFUNC CONCAT3(LAL,TYPECODE,ReadTArraySeries)
#define FFUNC CONCAT3(LAL,TYPECODE,ReadFSeries)
#define SCREATE CONCAT3(LAL,TYPECODE,CreateVector)
#define VCREATE CONCAT3(LAL,TYPECODE,CreateVectorSequence)
#define ACREATE CONCAT3(LAL,TYPECODE,CreateArraySequence)
#define SDESTROY CONCAT3(LAL,TYPECODE,DestroyVector)
#define VDESTROY CONCAT3(LAL,TYPECODE,DestroyVectorSequence)
#define ADESTROY CONCAT3(LAL,TYPECODE,DestroyArraySequence)
#define FMT CONCAT3(LAL_,DATA,_FORMAT)
#define STRINGTODATA CONCAT2(LALStringTo,DATACODE)


void
SFUNC ( LALStatus *stat, STYPE *series, FILE *stream )
{ 
  BufferList *head = NULL;  /* pointer to head of linked list of buffers */
  BufferList *here = NULL;  /* pointer to current position in list */
  CHARVector *line = NULL; /* current line being read */
  CHAR *start, *end; /* start and end of a token on a line */
  UINT4 length = 0;  /* number of sequence elements to be read */
  UINT4 n;           /* countdown index over elements read */
  DATA *data = NULL;        /* pointer to data in buffers */
  TYPE *sData = NULL;       /* pointer to data in output sequence */
  STYPE sCopy; /* internal copy of series */
  int numRead = 0; /* number of values read by parsing subroutine */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( series, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !(series->data), stat, STREAMINPUTH_EOUT,
	  STREAMINPUTH_MSGEOUT );
#if COMPLEX
  ASSERT( 2*(BUFFSIZE/SIZE/2) == BUFFSIZE/SIZE, stat,
	  STREAMINPUTH_EBUF, STREAMINPUTH_MSGEBUF );
#endif

  /*******************************************************************
   * PARSE METADATA HEADER                                           *
   *******************************************************************/

  sCopy = *series;

  /* Skip over blank lines; start points to the first non-whitespace
     character (or '\0' if there are none). */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
  start = line->data;
  while ( isspace( *start ) ) {
    start++;
    if ( *start == '\0' ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
      start = line->data;
    }
  }

  /* As long as '#' is the first non-whitespace character of any
     nonblank lines... */
  while ( *start == '#' ) {
    CHAR *startValue; /* start of substring giving metadata value */
    CHAR *endValue;   /* end of substring giving metadata value */

    /* Skip to the start of the metadata field tag. */
    do
      start++;
    while ( isspace( *start ) );

    /* Mark the end of the tag and the start of the metadata value. */
    end = start;
    while ( !isspace( *end ) && *end != '=' && *end != '\0' )
      end++;
    startValue = end;
    while ( isspace( *startValue ) || *startValue == '=' )
      startValue++;
    if ( startValue != end ) {
      *end = '\0';

      /* Parse name field. */
      if ( !strcmp( start, "name" ) ) {
	if ( *startValue != '"' )
	  LALWarning( stat, LALREADSERIESC_HEADER "name" );
	else
	  LALLiteralToString( stat->statusPtr, sCopy.name, startValue,
			      LALNameLength );
      }

      /* Parse epoch field. */
      else if ( !strcmp( start, "epoch" ) ) {
	INT8 sec, nsec;
	LALStringToI8( stat->statusPtr, &sec, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "epoch" );
	else {
	  CHAR *temp = endValue;
	  LALStringToI8( stat->statusPtr, &nsec, endValue, &endValue );
	  BEGINFAIL( stat ) {
	    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  } ENDFAIL( stat );
	  if ( endValue == temp ) {
	    sCopy.epoch.gpsSeconds = (INT4)( sec/1000000000LL );
	    sCopy.epoch.gpsNanoSeconds = (INT4)
	      ( sec - 1000000000LL*sCopy.epoch.gpsSeconds );
	  } else {
	    sCopy.epoch.gpsSeconds = (INT4)( sec );
	    sCopy.epoch.gpsNanoSeconds = (INT4)( nsec );
	  }
	}
      }

      /* Parse deltaT field. */
      else if ( !strcmp( start, "deltaT" ) ) {
	REAL8 deltaT;
	LALStringToD( stat->statusPtr, &deltaT, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "deltaT" );
	else
	  sCopy.deltaT = deltaT;
      }

      /* Parse f0 field. */
      else if ( !strcmp( start, "f0" ) ) {
	REAL8 f0;
	LALStringToD( stat->statusPtr, &f0, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "f0" );
	else
	  sCopy.f0 = f0;
      }

      /* Parse sampleUnits field. */
      else if ( !strcmp( start, "sampleUnits" ) ) {
	LALUnit unit;
	CHARVector unitString;
	endValue = ++startValue;
	while ( *endValue != '"' && *endValue != '\n' &&
		*endValue != '\0' )
	  endValue++;
	*endValue = '\0';
	unitString.length = strlen( startValue ) + 1;
	unitString.data = startValue;
	LALParseUnitString( stat->statusPtr, &unit, &unitString );
	if ( stat->statusPtr->statusCode == UNITSH_EPARSE ) {
#ifndef NDEBUG
	  if ( lalDebugLevel & LALERROR ) {
	    LALPrintError( "\tCONTINUE: Ignoring preceding error\n" );
	    DETATCHSTATUSPTR( stat );
	    ATTATCHSTATUSPTR( stat );
	  }
#endif
	  LALWarning( stat, LALREADSERIESC_HEADER "sampleUnits" );
	} else {
	  BEGINFAIL( stat ) {
	    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  } ENDFAIL( stat );
	  sCopy.sampleUnits = unit;
	}
      }

      /* Parse length field. */
      else if ( !strcmp( start, "length" ) ) {
	UINT4 tempLength;
	LALStringToU4( stat->statusPtr, &tempLength, startValue,
		       &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue || tempLength == 0 )
	  LALWarning( stat, LALREADSERIESC_HEADER "length" );
	else
	  length = tempLength;
      }

      /* Parse datatype field. */
#ifndef NDEBUG
      else if ( !strcmp( start, "datatype" ) ) {
	if ( lalDebugLevel & LALWARNING ) {
	  endValue = startValue;
	  while ( !isspace( *endValue ) && *endValue != '\0' )
	    endValue++;
	  *endValue = '\0';
	  if ( strcmp( startValue, "STYPE" ) ) {
	    LALWarning( stat, "STYPE data expected" );
	    LALPrintError( "\t%s data being read\n", startValue );
	  }
	}
      }
#endif

      /* No other recognized tags; ignore anything else. */
    }

    /* Read in next line, skipping over blank lines. */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
    start = line->data;
    while ( isspace( *start ) ) {
      start++;
      if ( *start == '\0' ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
	start = line->data;
      }
    }
  }


  /*******************************************************************
   * PARSE DATA                                                      *
   *******************************************************************/

  end = start;

  /* If length is given, we can simply allocate the required block of
     memory and parse data into it. */
  if ( length > 0 ) {
    SCREATE ( stat->statusPtr, &(sCopy.data), length );
    BEGINFAIL( stat ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    } ENDFAIL( stat );
    sData = sCopy.data->data;

    /* Begin with line in memory. */
    n = length;
    numRead = 0;
#if COMPLEX
    do {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      STRINGTODATA ( stat->statusPtr, &(zData.x->re), start = end,
		     &end );
      BEGINFAIL( stat ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
      } ENDFAIL( stat );
      if ( start != end ) {
	numRead = 1;
	STRINGTODATA ( stat->statusPtr, &(zData.x->im), start = end,
		       &end );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	} ENDFAIL( stat );
	if ( start != end ) {
	  numRead = 0;
	  sData++;
	}
      }
    } while ( end != start && --n );
#else
    do {
      STRINGTODATA ( stat->statusPtr, sData++, start = end, &end );
      BEGINFAIL( stat ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
      } ENDFAIL( stat );
    } while ( end != start && --n );
    sData--;
#endif

    /* Read remaining data with fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
#if COMPLEX
    if ( numRead == 1 ) {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      if ( fscanf( stream, "%" FMT, &(zData.x->im) ) != 1 ) {
	TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      sData++;
      n--;
    }
    while ( n-- ) {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      if ( ( fscanf( stream, "%" FMT, &(zData.x->re) ) != 1 ||
	   fscanf( stream, "%" FMT, &(zData.x->im) ) != 1 ) ) {
	TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      sData++;
    }
#else
    while ( n-- )
      if ( fscanf( stream, "%" FMT, sData++ ) != 1 ) {
	TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
#endif
  }


  /* If length was not specified, we must read data into buffers until
     we know how many there are. */
  else {
    here = head = (BufferList *)LALMalloc( sizeof(BufferList) );
    if ( !here ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here->next = NULL;
    data = here->buf.DATACODE;
    n = BUFFSIZE/SIZE;

    /* Read from the line currently in memory. */
    do {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	here = here->next;
	if ( !here ) {
	  FREEBUFFERLIST( head );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here->next = NULL;
	data = here->buf.DATACODE;
	length += n = BUFFSIZE/SIZE;
      } else
	data--;
    } while ( end != start );

    /* Read the remaining data using fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    do {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      length += BUFFSIZE/SIZE - n;
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	here = here->next;
	if ( !here ) {
	  FREEBUFFERLIST( head );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here->next = NULL;
	data = here->buf.DATACODE;
	n = BUFFSIZE/SIZE;
      }
    } while ( numRead == 1 );

    /* Create the data sequence, and copy the data into it. */
#if COMPLEX
    length /= 2;
#endif
    if ( !length ) {
      FREEBUFFERLIST( head );
      ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
    }
    SCREATE ( stat->statusPtr, &(sCopy.data), length );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head );
    ENDFAIL( stat );

    here = head;
    data = here->buf.DATACODE;
    sData = sCopy.data->data;
    n = BUFFSIZE/SIZE;
    while ( length ) {
#if COMPLEX
      n /= 2;
      if ( n > length )
	n = length;
      length -= n;
      while ( n-- ) {
	*sData = *(data++);
	*(sData++) += *(data++) * I;
      }
      here = here->next;
      if ( here )
	data = here->buf.DATACODE;
      n = BUFFSIZE/SIZE;
#else
      if ( n > length )
	n = length;
      memcpy( sData, here->buf.DATACODE, n*sizeof( TYPE ) );
      sData += n;
      length -= n;
      here = here->next;
      if ( here )
	data = here->buf.DATACODE;
      n = BUFFSIZE/SIZE;
#endif
    }
    FREEBUFFERLIST( head );
  }

  /* Data have been stored successfully.  So, clean up and exit. */
  *series = sCopy;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


void
VFUNC ( LALStatus *stat, VTYPE *series, FILE *stream )
{ 
  BufferList *head = NULL;  /* pointer to head of linked list of buffers */
  BufferList *here = NULL;  /* pointer to current position in list */
  CHARVector *line = NULL; /* current line being read */
  CHAR *start, *end; /* start and end of a token on a line */
  UINT4 length = 0;  /* number of sequence elements to be read */
  UINT4 vectorLength = 0; /* number of components per element */
  UINT4 nTot = 0, n; /* number of data read, and countdown index */
  DATA *data = NULL;        /* pointer to data in buffers */
  TYPE *sData = NULL;       /* pointer to data in output sequence */
  VTYPE sCopy; /* internal copy of series */
  CreateVectorSequenceIn in; /* structure to create sequence */
  int numRead = 0; /* number of values read per call of fscanf() */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( series, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !(series->data), stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );
#if COMPLEX
  ASSERT( 2*(BUFFSIZE/SIZE/2) == BUFFSIZE/SIZE, stat,
	  STREAMINPUTH_EBUF, STREAMINPUTH_MSGEBUF );
#endif

  /*******************************************************************
   * PARSE METADATA HEADER                                           *
   *******************************************************************/

  sCopy = *series;

  /* Skip over blank lines; start points to the first non-whitespace
     character (or '\0' if there are none). */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
  start = line->data;
  while ( isspace( *start ) ) {
    start++;
    if ( *start == '\0' ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
      start = line->data;
    }
  }

  /* As long as '#' is the first non-whitespace character of any
     nonblank lines... */
  while ( *start == '#' ) {
    CHAR *startValue; /* start of substring giving metadata value */
    CHAR *endValue;   /* end of substring giving metadata value */

    /* Skip to the start of the metadata field tag. */
    do
      start++;
    while ( isspace( *start ) );

    /* Mark the end of the tag and the start of the metadata value. */
    end = start;
    while ( !isspace( *end ) && *end != '=' && *end != '\0' )
      end++;
    startValue = end;
    while ( isspace( *startValue ) || *startValue == '=' )
      startValue++;
    if ( startValue != end ) {
      *end = '\0';

      /* Parse name field. */
      if ( !strcmp( start, "name" ) ) {
	if ( *startValue != '"' )
	  LALWarning( stat, LALREADSERIESC_HEADER "name" );
	else
	  LALLiteralToString( stat->statusPtr, sCopy.name, startValue,
			      LALNameLength );
      }

      /* Parse epoch field. */
      else if ( !strcmp( start, "epoch" ) ) {
	INT8 sec, nsec;
	LALStringToI8( stat->statusPtr, &sec, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "epoch" );
	else {
	  CHAR *temp = endValue;
	  LALStringToI8( stat->statusPtr, &nsec, endValue, &endValue );
	  BEGINFAIL( stat ) {
	    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  } ENDFAIL( stat );
	  if ( endValue == temp ) {
	    sCopy.epoch.gpsSeconds = (INT4)( sec/1000000000LL );
	    sCopy.epoch.gpsNanoSeconds = (INT4)
	      ( sec - 1000000000LL*sCopy.epoch.gpsSeconds );
	  } else {
	    sCopy.epoch.gpsSeconds = (INT4)( sec );
	    sCopy.epoch.gpsNanoSeconds = (INT4)( nsec );
	  }
	}
      }

      /* Parse deltaT field. */
      else if ( !strcmp( start, "deltaT" ) ) {
	REAL8 deltaT;
	LALStringToD( stat->statusPtr, &deltaT, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "deltaT" );
	else
	  sCopy.deltaT = deltaT;
      }

      /* Parse f0 field. */
      else if ( !strcmp( start, "f0" ) ) {
	REAL8 f0;
	LALStringToD( stat->statusPtr, &f0, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "f0" );
	else
	  sCopy.f0 = f0;
      }

      /* Parse sampleUnits field. */
      else if ( !strcmp( start, "sampleUnits" ) ) {
	LALUnit unit;
	CHARVector unitString;
	endValue = ++startValue;
	while ( *endValue != '"' && *endValue != '\n' &&
		*endValue != '\0' )
	  endValue++;
	*endValue = '\0';
	unitString.length = strlen( startValue ) + 1;
	unitString.data = startValue;
	LALParseUnitString( stat->statusPtr, &unit, &unitString );
	if ( stat->statusPtr->statusCode == UNITSH_EPARSE ) {
#ifndef NDEBUG
	  if ( lalDebugLevel & LALERROR ) {
	    LALPrintError( "\tCONTINUE: Ignoring preceding error\n" );
	    DETATCHSTATUSPTR( stat );
	    ATTATCHSTATUSPTR( stat );
	  }
#endif
	  LALWarning( stat, LALREADSERIESC_HEADER "sampleUnits" );
	} else {
	  BEGINFAIL( stat ) {
	    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  } ENDFAIL( stat );
	  sCopy.sampleUnits = unit;
	}
      }

      /* Parse length field. */
      else if ( !strcmp( start, "length" ) ) {
	UINT4 tempLength;
	LALStringToU4( stat->statusPtr, &tempLength, startValue,
		       &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue || tempLength == 0 )
	  LALWarning( stat, LALREADSERIESC_HEADER "length" );
	else
	  length = tempLength;
      }

      /* Parse vectorLength field. */
      else if ( !strcmp( start, "vectorLength" ) ) {
	UINT4 tempLength;
	LALStringToU4( stat->statusPtr, &tempLength, startValue,
		       &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue || tempLength == 0 )
	  LALWarning( stat, LALREADSERIESC_HEADER "vectorLength" );
	else
	  vectorLength = tempLength;
      }

      /* Parse datatype field. */
#ifndef NDEBUG
      else if ( !strcmp( start, "datatype" ) ) {
	if ( lalDebugLevel & LALWARNING ) {
	  endValue = startValue;
	  while ( !isspace( *endValue ) && *endValue != '\0' )
	    endValue++;
	  *endValue = '\0';
	  if ( strcmp( startValue, "VTYPE" ) ) {
	    LALWarning( stat, "VTYPE data expected" );
	    LALPrintError( "\t%s data being read\n", startValue );
	  }
	}
      }
#endif

      /* No other recognized tags; ignore anything else. */
    }

    /* Read in next line, skipping over blank lines. */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
    start = line->data;
    while ( isspace( *start ) ) {
      start++;
      if ( *start == '\0' ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
	start = line->data;
      }
    }
  }


  /*******************************************************************
   * PARSE DATA                                                      *
   *******************************************************************/

  end = start;
  in.length = length;
  in.vectorLength = vectorLength;

  /* If length and vectorLength were both specified, we know we will
     need length*vectorLength data, so we can just allocate this space
     from the start. */
  if ( length > 0 && vectorLength > 0 ) {
    VCREATE ( stat->statusPtr, &(sCopy.data), &in );
    BEGINFAIL( stat ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    } ENDFAIL( stat );
    sData = sCopy.data->data;

    /* Read the line in memory. */
    n = length*vectorLength;
    numRead = 0;
#if COMPLEX
    do {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      STRINGTODATA ( stat->statusPtr, &(zData.x->re), start = end,
		     &end );
      BEGINFAIL( stat ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( VDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
      } ENDFAIL( stat );
      if ( start != end ) {
	numRead = 1;
	STRINGTODATA ( stat->statusPtr, &(zData.x->im), start = end,
		       &end );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  TRY( VDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	} ENDFAIL( stat );
	if ( start != end ) {
	  numRead = 0;
	  sData++;
	}
      }
    } while ( end != start && --n );
#else
    do {
      STRINGTODATA ( stat->statusPtr, sData++, start = end, &end );
      BEGINFAIL( stat ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( VDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
      } ENDFAIL( stat );
    } while ( end != start && --n );
    sData--;
#endif
  }

  /* Otherwise, we will have to read at least the entire first line
     currently in memory. */
  else {
    here = head = (BufferList *)LALMalloc( sizeof(BufferList) );
    if ( !here ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here->next = NULL;
    data = here->buf.DATACODE;
    n = BUFFSIZE/SIZE;
    do {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	here = here->next;
	if ( !here ) {
	  FREEBUFFERLIST( head );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here->next = NULL;
	data = here->buf.DATACODE;
	nTot += n = BUFFSIZE/SIZE;
      } else
	data--;
    } while ( end != start );
  }
  nTot += BUFFSIZE/SIZE - n;

  /* We can now compute vectorLength, if it wasn't specified. */
  if ( in.vectorLength == 0 ) {
    in.vectorLength = vectorLength = nTot;
    if ( in.vectorLength == 0 ) {
      FREEBUFFERLIST( head );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
    }
#if COMPLEX
    vectorLength /= 2;
    if ( 2*vectorLength < in.vectorLength ) {
      FREEBUFFERLIST( head );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EVLEN, STREAMINPUTH_MSGEVLEN );
    }
    in.vectorLength = vectorLength;
#endif

    /* If length was specified, we can now create the vector sequence,
       and copy the first line into it. */
    if ( length > 0 ) {
      VCREATE ( stat->statusPtr, &(sCopy.data), &in );
      BEGINFAIL( stat ) {
	FREEBUFFERLIST( head );
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      } ENDFAIL( stat );
      here = head;
      data = here->buf.DATACODE;
      sData = sCopy.data->data;
      n = BUFFSIZE/SIZE;
      while ( vectorLength ) {
#if COMPLEX
	n /= 2;
	if ( n > vectorLength )
	  n = vectorLength;
	vectorLength -= n;
	while ( n-- ) {
	  *sData = *(data++);
	  *(sData++) += *(data++) * I;
	}
	here = here->next;
	if ( here )
	  data = here->buf.DATACODE;
	n = BUFFSIZE/SIZE;
#else
	if ( n > vectorLength )
	  n = vectorLength;
	memcpy( sData, here->buf.DATACODE, n*sizeof( TYPE ) );
	sData += n;
	vectorLength -= n;
	here = here->next;
	if ( here )
	  data = here->buf.DATACODE;
	n = BUFFSIZE/SIZE;
#endif
      }
      FREEBUFFERLIST( head );
      n = ( length - 1 )*in.vectorLength;
      numRead = 0;
    }
  }

  /* Read the remaining data with fscanf().  If length was specified,
     then the data sequence has already been created, and n is the
     number of atomic (integer, floating-point, or complex) data still
     to be read. */
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
  if ( length > 0 ) {
#if COMPLEX
    if ( numRead == 1 ) {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      if ( fscanf( stream, "%" FMT, &(zData.x->im) ) != 1 ) {
	TRY( VDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      sData++;
      n--;
    }
    while ( n-- ) {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      if ( fscanf( stream, "%" FMT, &(zData.x->re) ) != 1 ||
	   fscanf( stream, "%" FMT, &(zData.x->im) ) != 1 ) {
	TRY( VDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      sData++;
    }
#else
    while ( n-- )
      if ( fscanf( stream, "%" FMT, sData++ ) != 1 ) {
	TRY( VDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
#endif
  }

  /* Otherwise, if length wasn't specified, we need to continue
     reading into buffers.  In this case n is the space remaining in
     the current buffer. */
  else {
    length = nTot + n;
    while ( numRead == 1 ) {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      length -= n;
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	length += n = BUFFSIZE/SIZE;
      }
    }

    /* Now create the data sequence and copy data into it. */
#if COMPLEX
    length /= 2;
#endif
    in.length = length/in.vectorLength;
    length = in.length*in.vectorLength;
    VCREATE ( stat->statusPtr, &(sCopy.data), &in );
    BEGINFAIL( stat ) {
      FREEBUFFERLIST( head );
    } ENDFAIL( stat );
    here = head;
    data = here->buf.DATACODE;
    sData = sCopy.data->data;
    n = BUFFSIZE/SIZE;
    while ( length ) {
#if COMPLEX
      n /= 2;
      if ( n > length )
	n = length;
      length -= n;
      while ( n-- ) {
	*sData = *(data++);
	*(sData++) += *(data++) * I;
      }
      here = here->next;
      if ( here )
	data = here->buf.DATACODE;
      n = BUFFSIZE/SIZE;
#else
      if ( n > length )
	n = length;
      memcpy( sData, here->buf.DATACODE, n*sizeof( TYPE ) );
      sData += n;
      length -= n;
      here = here->next;
      if ( here )
	data = here->buf.DATACODE;
      n = BUFFSIZE/SIZE;
#endif
    }
    FREEBUFFERLIST( head );
  }

  /* Data have been stored successfully.  So, clean up and exit. */
  *series = sCopy;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


void
AFUNC ( LALStatus *stat, ATYPE *series, FILE *stream )
{ 
  BufferList *head;  /* pointer to head of linked list of buffers */
  BufferList *here;  /* pointer to current position in list */
  CHARVector *line = NULL; /* current line being read */
  CHAR *start, *end; /* start and end of a token on a line */
  UINT4 length = 0;  /* number of sequence elements to be read */
  UINT4 nTot = 0, n; /* number of data read, and countdown index */
  UINT4 arrayDim = 0; /* number of components per element */
  UINT4Vector *dimLength = NULL; /* number of components per index */
  DATA *data;        /* pointer to data in buffers */
  TYPE *sData;       /* pointer to data in output sequence */
  ATYPE sCopy; /* internal copy of series */
  CreateArraySequenceIn in; /* structure to create sequence */
  int numRead; /* number of values read per call of fscanf() */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( series, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !(series->data), stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );
#if COMPLEX
  ASSERT( 2*(BUFFSIZE/SIZE/2) == BUFFSIZE/SIZE, stat,
	  STREAMINPUTH_EBUF, STREAMINPUTH_MSGEBUF );
#endif

  /*******************************************************************
   * PARSE METADATA HEADER                                           *
   *******************************************************************/

  sCopy = *series;

  /* Skip over blank lines; start points to the first non-whitespace
     character (or '\0' if there are none). */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
  start = line->data;
  while ( isspace( *start ) ) {
    start++;
    if ( *start == '\0' ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
      start = line->data;
    }
  }

  /* As long as '#' is the first non-whitespace character of any
     nonblank lines... */
  while ( *start == '#' ) {
    CHAR *startValue; /* start of substring giving metadata value */
    CHAR *endValue;   /* end of substring giving metadata value */

    /* Skip to the start of the metadata field tag. */
    do
      start++;
    while ( isspace( *start ) );

    /* Mark the end of the tag and the start of the metadata value. */
    end = start;
    while ( !isspace( *end ) && *end != '=' && *end != '\0' )
      end++;
    startValue = end;
    while ( isspace( *startValue ) || *startValue == '=' )
      startValue++;
    if ( startValue != end ) {
      *end = '\0';

      /* Parse name field. */
      if ( !strcmp( start, "name" ) ) {
	if ( *startValue != '"' )
	  LALWarning( stat, LALREADSERIESC_HEADER "name" );
	else
	  LALLiteralToString( stat->statusPtr, sCopy.name, startValue,
			      LALNameLength );
      }

      /* Parse epoch field. */
      else if ( !strcmp( start, "epoch" ) ) {
	INT8 sec, nsec;
	LALStringToI8( stat->statusPtr, &sec, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  if ( dimLength ) {
	    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	  }
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "epoch" );
	else {
	  CHAR *temp = endValue;
	  LALStringToI8( stat->statusPtr, &nsec, endValue, &endValue );
	  BEGINFAIL( stat ) {
	    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	    if ( dimLength ) {
	      TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	    }
	  } ENDFAIL( stat );
	  if ( endValue == temp ) {
	    sCopy.epoch.gpsSeconds = (INT4)( sec/1000000000LL );
	    sCopy.epoch.gpsNanoSeconds = (INT4)
	      ( sec - 1000000000LL*sCopy.epoch.gpsSeconds );
	  } else {
	    sCopy.epoch.gpsSeconds = (INT4)( sec );
	    sCopy.epoch.gpsNanoSeconds = (INT4)( nsec );
	  }
	}
      }

      /* Parse deltaT field. */
      else if ( !strcmp( start, "deltaT" ) ) {
	REAL8 deltaT;
	LALStringToD( stat->statusPtr, &deltaT, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  if ( dimLength ) {
	    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	  }
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "deltaT" );
	else
	  sCopy.deltaT = deltaT;
      }

      /* Parse f0 field. */
      else if ( !strcmp( start, "f0" ) ) {
	REAL8 f0;
	LALStringToD( stat->statusPtr, &f0, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  if ( dimLength ) {
	    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	  }
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "f0" );
	else
	  sCopy.f0 = f0;
      }

      /* Parse sampleUnits field. */
      else if ( !strcmp( start, "sampleUnits" ) ) {
	LALUnit unit;
	CHARVector unitString;
	endValue = ++startValue;
	while ( *endValue != '"' && *endValue != '\n' &&
		*endValue != '\0' )
	  endValue++;
	*endValue = '\0';
	unitString.length = strlen( startValue ) + 1;
	unitString.data = startValue;
	LALParseUnitString( stat->statusPtr, &unit, &unitString );
	if ( stat->statusPtr->statusCode == UNITSH_EPARSE ) {
#ifndef NDEBUG
	  if ( lalDebugLevel & LALERROR ) {
	    LALPrintError( "\tCONTINUE: Ignoring preceding error\n" );
	    DETATCHSTATUSPTR( stat );
	    ATTATCHSTATUSPTR( stat );
	  }
#endif
	  LALWarning( stat, LALREADSERIESC_HEADER "sampleUnits" );
	} else {
	  BEGINFAIL( stat ) {
	    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	    if ( dimLength ) {
	      TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	    }
	  } ENDFAIL( stat );
	  sCopy.sampleUnits = unit;
	}
      }

      /* Parse length field. */
      else if ( !strcmp( start, "length" ) ) {
	UINT4 tempLength;
	LALStringToU4( stat->statusPtr, &tempLength, startValue,
		       &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  if ( dimLength ) {
	    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ),
		 stat );
	  }
	} ENDFAIL( stat );
	if ( endValue == startValue || tempLength == 0 )
	  LALWarning( stat, LALREADSERIESC_HEADER "length" );
	else
	  length = tempLength;
      }

      /* Parse arrayDim field. */
      else if ( !strcmp( start, "arrayDim" ) ) {
	UINT4 tempLength;
	LALStringToU4( stat->statusPtr, &tempLength, startValue,
		       &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  if ( dimLength ) {
	    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ),
		 stat );
	  }
	} ENDFAIL( stat );
	if ( endValue == startValue || tempLength == 0 )
	  LALWarning( stat, LALREADSERIESC_HEADER "arrayDim" );
	else
	  arrayDim = tempLength;
      }

      /* Parse dimLength field. */
      else if ( !strcmp( start, "dimLength" ) ) {
	UINT4 *dData;
	here = head = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !here ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  if ( dimLength ) {
	    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ),
		 stat );
	  }
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here->next = NULL;
	dData = here->buf.U4;
	n = BUFFSIZE/4;
	endValue = startValue;
	/* Read components into the buffer list. */
	do {
	  do {
	    LALStringToU4( stat->statusPtr, dData++,
			   startValue = endValue, &endValue );
	    BEGINFAIL( stat ) {
	      FREEBUFFERLIST( head );
	      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	      if ( dimLength ) {
		TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ),
		     stat );
	      }
	    } ENDFAIL( stat );
	  } while ( endValue != startValue && --n );
	  if ( !n ) {
	    here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	    here = here->next;
	    if ( !here ) {
	      FREEBUFFERLIST( head );
	      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	      if ( dimLength ) {
		TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ),
		     stat );
	      }
	      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	    }
	    here->next = NULL;
	    dData = here->buf.U4;
	    nTot += n = BUFFSIZE/4;
	  } else
	    nTot += BUFFSIZE/4 - n;
	} while ( endValue != startValue );
	/* Copy components into dimLength vector. */
	if ( nTot == 0 )
	  LALWarning( stat, LALREADSERIESC_HEADER "dimLength" );
	else {
	  if ( dimLength ) {
	    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ),
		 stat );
	  }
	  dimLength = NULL;
	  LALU4CreateVector( stat->statusPtr, &dimLength, nTot );
	  BEGINFAIL( stat ) {
	    FREEBUFFERLIST( head );
	    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  } ENDFAIL( stat );
	  here = head;
	  dData = dimLength->data;
	  n = BUFFSIZE/4;
	  if ( n > nTot )
	    n = nTot;
	  while ( n ) {
	    memcpy( dData, here->buf.U4, n*sizeof(UINT4) );
	    dData += n;
	    nTot -= n;
	    here = here->next;
	    n = BUFFSIZE/4;
	    if ( n > nTot )
	      n = nTot;
	  }
	}
	/* dimLength is complete, so reset everything else. */
	FREEBUFFERLIST( head );
	nTot = 0;
      }

      /* Parse datatype field. */
#ifndef NDEBUG
      else if ( !strcmp( start, "datatype" ) ) {
	if ( lalDebugLevel & LALWARNING ) {
	  endValue = startValue;
	  while ( !isspace( *endValue ) && *endValue != '\0' )
	    endValue++;
	  *endValue = '\0';
	  if ( strcmp( startValue, "ATYPE" ) ) {
	    LALWarning( stat, "ATYPE data expected" );
	    LALPrintError( "\t%s data being read\n", startValue );
	  }
	}
      }
#endif

      /* No other recognized tags; ignore anything else. */
    }

    /* Read in next line, skipping over blank lines. */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    LALCHARReadVector( stat->statusPtr, &line, stream );
    BEGINFAIL( stat ) {
      if ( dimLength ) {
	TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
      }
    } ENDFAIL( stat );
    start = line->data;
    while ( isspace( *start ) ) {
      start++;
      if ( *start == '\0' ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	LALCHARReadVector( stat->statusPtr, &line, stream );
	BEGINFAIL( stat ) {
	  if ( dimLength ) {
	    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	  }
	} ENDFAIL( stat );
	start = line->data;
      }
    }
  }


  /*******************************************************************
   * PARSE DATA                                                      *
   *******************************************************************/

  end = start;

  /* Make sure dimLength exists, and is consistent with arrayDim (if
     specified). */
  if ( !dimLength ) {
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    ABORT( stat, STREAMINPUTH_EDLEN, STREAMINPUTH_MSGEDLEN );
  } else {
    UINT4 i;       /* index over array indecies */
    UINT4 dim = 1; /* total array outer-product dimension */
    for ( i = 0; i < dimLength->length; i++ )
      dim *= dimLength->data[i];
    if ( dim == 0 || ( arrayDim != 0 && arrayDim != dim ) ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
      ABORT( stat, STREAMINPUTH_EDIM, STREAMINPUTH_MSGEDIM );
    }
    arrayDim = dim;
  }

  /* If length is given, we can simply allocate the required block of
     memory and parse data into it. */
  if ( length > 0 ) {
    in.length = length;
    in.dimLength = dimLength;
    ACREATE ( stat->statusPtr, &(sCopy.data), &in );
    BEGINFAIL( stat ) {
      TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    } ENDFAIL( stat );
    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
    sData = sCopy.data->data;

    /* Begin with line in memory. */
    n = length*arrayDim;
    numRead = 0;
#if COMPLEX
    do {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      STRINGTODATA ( stat->statusPtr, &(zData.x->re), start = end,
		     &end );
      BEGINFAIL( stat ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( ADESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
      } ENDFAIL( stat );
      if ( start != end ) {
	numRead = 1;
	STRINGTODATA ( stat->statusPtr, &(zData.x->im), start = end,
		       &end );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  TRY( ADESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	} ENDFAIL( stat );
	if ( start != end ) {
	  numRead = 0;
	  sData++;
	}
      }
    } while ( end != start && --n );
#else
    do {
      STRINGTODATA ( stat->statusPtr, sData++, start = end, &end );
      BEGINFAIL( stat ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( ADESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
      } ENDFAIL( stat );
    } while ( end != start && --n );
    sData--;
#endif

    /* Read remaining data with fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
#if COMPLEX
    if ( numRead == 1 ) {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      if ( fscanf( stream, "%" FMT, &(zData.x->im) ) != 1 ) {
	TRY( ADESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      sData++;
      n--;
    }
    while ( n-- ) {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      if ( fscanf( stream, "%" FMT, &(zData.x->re) ) != 1 ||
	   fscanf( stream, "%" FMT, &(zData.x->im) ) != 1 ) {
	TRY( ADESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      sData++;
    }
#else
    while ( n-- )
      if ( fscanf( stream, "%" FMT, sData++ ) != 1 ) {
	TRY( ADESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
#endif
  }


  /* If length was not specified, we must read data into buffers until
     we know how many there are. */
  else {
    here = head = (BufferList *)LALMalloc( sizeof(BufferList) );
    if ( !here ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here->next = NULL;
    data = here->buf.DATACODE;
    n = BUFFSIZE/SIZE;

    /* Read from the line currently in memory. */
    do {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	here = here->next;
	if ( !here ) {
	  FREEBUFFERLIST( head );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here->next = NULL;
	data = here->buf.DATACODE;
	length += n = BUFFSIZE/SIZE;
      } else
	data--;
    } while ( end != start );

    /* Read the remaining data using fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    do {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      length += BUFFSIZE/SIZE - n;
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	here = here->next;
	if ( !here ) {
	  FREEBUFFERLIST( head );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here->next = NULL;
	data = here->buf.DATACODE;
	n = BUFFSIZE/SIZE;
      }
    } while ( numRead == 1 );

    /* Create the data sequence, and copy the data into it. */
    length /= arrayDim;
#if COMPLEX
    length /= 2;
#endif
    if ( !length ) {
      FREEBUFFERLIST( head );
      ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
    }
    in.length = length;
    in.dimLength = dimLength;
    ACREATE ( stat->statusPtr, &(sCopy.data), &in );
    BEGINFAIL( stat ) {
      FREEBUFFERLIST( head );
      TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
    } ENDFAIL( stat );
    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
    here = head;
    data = here->buf.DATACODE;
    sData = sCopy.data->data;
    length *= arrayDim;
    n = BUFFSIZE/SIZE;
    while ( length ) {
#if COMPLEX
      n /= 2;
      if ( n > length )
	n = length;
      length -= n;
      while ( n-- ) {
	*sData = *(data++);
	*(sData++) += *(data++) * I;
      }
      here = here->next;
      if ( here )
	data = here->buf.DATACODE;
      n = BUFFSIZE/SIZE;
#else
      if ( n > length )
	n = length;
      memcpy( sData, here->buf.DATACODE, n*sizeof( TYPE ) );
      sData += n;
      length -= n;
      here = here->next;
      if ( here )
	data = here->buf.DATACODE;
      n = BUFFSIZE/SIZE;
#endif
    }
    FREEBUFFERLIST( head );
  }

  /* Data have been stored successfully.  So, clean up and exit. */
  *series = sCopy;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


void
FFUNC ( LALStatus *stat, FTYPE *series, FILE *stream )
{ 
  BufferList *head;  /* pointer to head of linked list of buffers */
  BufferList *here;  /* pointer to current position in list */
  CHARVector *line = NULL; /* current line being read */
  CHAR *start, *end; /* start and end of a token on a line */
  UINT4 length = 0;  /* number of sequence elements to be read */
  UINT4 n;           /* countdown index over elements read */
  DATA *data;        /* pointer to data in buffers */
  TYPE *sData;       /* pointer to data in output sequence */
  FTYPE sCopy; /* internal copy of series */
  int numRead; /* number of values read per call of fscanf() */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( series, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !(series->data), stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );
#if COMPLEX
  ASSERT( 2*(BUFFSIZE/SIZE/2) == BUFFSIZE/SIZE, stat,
	  STREAMINPUTH_EBUF, STREAMINPUTH_MSGEBUF );
#endif

  /*******************************************************************
   * PARSE METADATA HEADER                                           *
   *******************************************************************/

  sCopy = *series;

  /* Skip over blank lines; start points to the first non-whitespace
     character (or '\0' if there are none). */
  TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
  start = line->data;
  while ( isspace( *start ) ) {
    start++;
    if ( *start == '\0' ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
      start = line->data;
    }
  }

  /* As long as '#' is the first non-whitespace character of any
     nonblank lines... */
  while ( *start == '#' ) {
    CHAR *startValue; /* start of substring giving metadata value */
    CHAR *endValue;   /* end of substring giving metadata value */

    /* Skip to the start of the metadata field tag. */
    do
      start++;
    while ( isspace( *start ) );

    /* Mark the end of the tag and the start of the metadata value. */
    end = start;
    while ( !isspace( *end ) && *end != '=' && *end != '\0' )
      end++;
    startValue = end;
    while ( isspace( *startValue ) || *startValue == '=' )
      startValue++;
    if ( startValue != end ) {
      *end = '\0';

      /* Parse name field. */
      if ( !strcmp( start, "name" ) ) {
	if ( *startValue != '"' )
	  LALWarning( stat, LALREADSERIESC_HEADER "name" );
	else
	  LALLiteralToString( stat->statusPtr, sCopy.name, startValue,
			      LALNameLength );
      }

      /* Parse epoch field. */
      else if ( !strcmp( start, "epoch" ) ) {
	INT8 sec, nsec;
	LALStringToI8( stat->statusPtr, &sec, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "epoch" );
	else {
	  CHAR *temp = endValue;
	  LALStringToI8( stat->statusPtr, &nsec, endValue, &endValue );
	  BEGINFAIL( stat ) {
	    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  } ENDFAIL( stat );
	  if ( endValue == temp ) {
	    sCopy.epoch.gpsSeconds = (INT4)( sec/1000000000LL );
	    sCopy.epoch.gpsNanoSeconds = (INT4)
	      ( sec - 1000000000LL*sCopy.epoch.gpsSeconds );
	  } else {
	    sCopy.epoch.gpsSeconds = (INT4)( sec );
	    sCopy.epoch.gpsNanoSeconds = (INT4)( nsec );
	  }
	}
      }

      /* Parse f0 field. */
      else if ( !strcmp( start, "f0" ) ) {
	REAL8 f0;
	LALStringToD( stat->statusPtr, &f0, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "f0" );
	else
	  sCopy.f0 = f0;
      }

      /* Parse deltaF field. */
      else if ( !strcmp( start, "deltaF" ) ) {
	REAL8 deltaF;
	LALStringToD( stat->statusPtr, &deltaF, startValue, &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue )
	  LALWarning( stat, LALREADSERIESC_HEADER "deltaF" );
	else
	  sCopy.deltaF = deltaF;
      }

      /* Parse sampleUnits field. */
      else if ( !strcmp( start, "sampleUnits" ) ) {
	LALUnit unit;
	CHARVector unitString;
	endValue = ++startValue;
	while ( *endValue != '"' && *endValue != '\n' &&
		*endValue != '\0' )
	  endValue++;
	*endValue = '\0';
	unitString.length = strlen( startValue ) + 1;
	unitString.data = startValue;
	LALParseUnitString( stat->statusPtr, &unit, &unitString );
	if ( stat->statusPtr->statusCode == UNITSH_EPARSE ) {
#ifndef NDEBUG
	  if ( lalDebugLevel & LALERROR ) {
	    LALPrintError( "\tCONTINUE: Ignoring preceding error\n" );
	    DETATCHSTATUSPTR( stat );
	    ATTATCHSTATUSPTR( stat );
	  }
#endif
	  LALWarning( stat, LALREADSERIESC_HEADER "sampleUnits" );
	} else {
	  BEGINFAIL( stat ) {
	    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  } ENDFAIL( stat );
	  sCopy.sampleUnits = unit;
	}
      }

      /* Parse length field. */
      else if ( !strcmp( start, "length" ) ) {
	UINT4 tempLength;
	LALStringToU4( stat->statusPtr, &tempLength, startValue,
		       &endValue );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
	if ( endValue == startValue || tempLength == 0 )
	  LALWarning( stat, LALREADSERIESC_HEADER "length" );
	else
	  length = tempLength;
      }

      /* Parse datatype field. */
#ifndef NDEBUG
      else if ( !strcmp( start, "datatype" ) ) {
	if ( lalDebugLevel & LALWARNING ) {
	  endValue = startValue;
	  while ( !isspace( *endValue ) && *endValue != '\0' )
	    endValue++;
	  *endValue = '\0';
	  if ( strcmp( startValue, "FTYPE" ) ) {
	    LALWarning( stat, "FTYPE data expected" );
	    LALPrintError( "\t%s data being read\n", startValue );
	  }
	}
      }
#endif

      /* No other recognized tags; ignore anything else. */
    }

    /* Read in next line, skipping over blank lines. */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
    start = line->data;
    while ( isspace( *start ) ) {
      start++;
      if ( *start == '\0' ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( LALCHARReadVector( stat->statusPtr, &line, stream ), stat );
	start = line->data;
      }
    }
  }


  /*******************************************************************
   * PARSE DATA                                                      *
   *******************************************************************/

  end = start;

  /* If length is given, we can simply allocate the required block of
     memory and parse data into it. */
  if ( length > 0 ) {
    SCREATE ( stat->statusPtr, &(sCopy.data), length );
    BEGINFAIL( stat ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    } ENDFAIL( stat );
    sData = sCopy.data->data;

    /* Begin with line in memory. */
    n = length;
    numRead = 0;
#if COMPLEX
    do {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      STRINGTODATA ( stat->statusPtr, &(zData.x->re), start = end,
		     &end );
      BEGINFAIL( stat ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
      } ENDFAIL( stat );
      if ( start != end ) {
	numRead = 1;
	STRINGTODATA ( stat->statusPtr, &(zData.x->im), start = end,
		       &end );
	BEGINFAIL( stat ) {
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	} ENDFAIL( stat );
	if ( start != end ) {
	  numRead = 0;
	  sData++;
	}
      }
    } while ( end != start && --n );
#else
    do {
      STRINGTODATA ( stat->statusPtr, sData++, start = end, &end );
      BEGINFAIL( stat ) {
	TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
      } ENDFAIL( stat );
    } while ( end != start && --n );
    sData--;
#endif

    /* Read remaining data with fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
#if COMPLEX
    if ( numRead == 1 ) {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      if ( fscanf( stream, "%" FMT, &(zData.x->im) ) != 1 ) {
	TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      sData++;
      n--;
    }
    while ( n-- ) {
      union { TYPE *z; struct { DATA re; DATA im; } *x; } zData;
      zData.z = sData;
      if ( fscanf( stream, "%" FMT, &(zData.x->re) ) != 1 ||
	   fscanf( stream, "%" FMT, &(zData.x->im) ) != 1 ) {
	TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      sData++;
    }
#else
    while ( n-- )
      if ( fscanf( stream, "%" FMT, sData++ ) != 1 ) {
	TRY( SDESTROY ( stat->statusPtr, &(sCopy.data) ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
#endif
  }


  /* If length was not specified, we must read data into buffers until
     we know how many there are. */
  else {
    here = head = (BufferList *)LALMalloc( sizeof(BufferList) );
    if ( !here ) {
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
    }
    here->next = NULL;
    data = here->buf.DATACODE;
    n = BUFFSIZE/SIZE;

    /* Read from the line currently in memory. */
    do {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	here = here->next;
	if ( !here ) {
	  FREEBUFFERLIST( head );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here->next = NULL;
	data = here->buf.DATACODE;
	length += n = BUFFSIZE/SIZE;
      } else
	data--;
    } while ( end != start );

    /* Read the remaining data using fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    do {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      length += BUFFSIZE/SIZE - n;
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	here = here->next;
	if ( !here ) {
	  FREEBUFFERLIST( head );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here->next = NULL;
	data = here->buf.DATACODE;
	n = BUFFSIZE/SIZE;
      }
    } while ( numRead == 1 );

    /* Create the data sequence, and copy the data into it. */
#if COMPLEX
    length /= 2;
#endif
    if ( !length ) {
      FREEBUFFERLIST( head );
      ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
    }
    SCREATE ( stat->statusPtr, &(sCopy.data), length );
    BEGINFAIL( stat )
      FREEBUFFERLIST( head );
    ENDFAIL( stat );

    here = head;
    data = here->buf.DATACODE;
    sData = sCopy.data->data;
    n = BUFFSIZE/SIZE;
    while ( length ) {
#if COMPLEX
      n /= 2;
      if ( n > length )
	n = length;
      length -= n;
      while ( n-- ) {
	*sData = *(data++);
	*(sData++) += *(data++) * I;
      }
      here = here->next;
      if ( here )
	data = here->buf.DATACODE;
      n = BUFFSIZE/SIZE;
#else
      if ( n > length )
	n = length;
      memcpy( sData, here->buf.DATACODE, n*sizeof( TYPE ) );
      sData += n;
      length -= n;
      here = here->next;
      if ( here )
	data = here->buf.DATACODE;
      n = BUFFSIZE/SIZE;
#endif
    }
    FREEBUFFERLIST( head );
  }

  /* Data have been stored successfully.  So, clean up and exit. */
  *series = sCopy;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
