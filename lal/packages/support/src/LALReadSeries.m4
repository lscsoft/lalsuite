/* m4_dnl $Id$
m4_ifelse(TYPECODE,`Z',`m4_define(`TYPE',`COMPLEX16')m4_define(`SIZE',`8')')m4_dnl
m4_ifelse(TYPECODE,`C',`m4_define(`TYPE',`COMPLEX8')m4_define(`SIZE',`4')')m4_dnl
m4_ifelse(TYPECODE,`D',`m4_define(`TYPE',`REAL8')m4_define(`SIZE',`8')')m4_dnl
m4_ifelse(TYPECODE,`S',`m4_define(`TYPE',`REAL4')m4_define(`SIZE',`4')')m4_dnl
m4_ifelse(TYPECODE,`I2',`m4_define(`TYPE',`INT2')m4_define(`SIZE',`2')')m4_dnl
m4_ifelse(TYPECODE,`I4',`m4_define(`TYPE',`INT4')m4_define(`SIZE',`4')')m4_dnl
m4_ifelse(TYPECODE,`I8',`m4_define(`TYPE',`INT8')m4_define(`SIZE',`8')')m4_dnl
m4_ifelse(TYPECODE,`U2',`m4_define(`TYPE',`UINT2')m4_define(`SIZE',`2')')m4_dnl
m4_ifelse(TYPECODE,`U4',`m4_define(`TYPE',`UINT4')m4_define(`SIZE',`4')')m4_dnl
m4_ifelse(TYPECODE,`U8',`m4_define(`TYPE',`UINT8')m4_define(`SIZE',`8')')m4_dnl
m4_define(`DATACODE',TYPECODE)m4_define(`DATA',TYPE)m4_define(`COMPLEX',`0')m4_dnl
m4_ifelse(TYPECODE,`Z',`m4_define(`DATACODE',`D')m4_define(`DATA',`REAL8')m4_define(`COMPLEX',`1')')m4_dnl
m4_ifelse(TYPECODE,`C',`m4_define(`DATACODE',`S')m4_define(`DATA',`REAL4')m4_define(`COMPLEX',`1')')m4_dnl
m4_define(`STYPE',`m4_format(`%sTimeSeries',TYPE)')m4_dnl
m4_define(`VTYPE',`m4_format(`%sTimeVectorSeries',TYPE)')m4_dnl
m4_define(`ATYPE',`m4_format(`%sTimeArraySeries',TYPE)')m4_dnl
m4_define(`FTYPE',`m4_format(`%sFrequencySeries',TYPE)')m4_dnl
m4_define(`SFUNC',`m4_format(`LAL%sReadTSeries',TYPECODE)')m4_dnl
m4_define(`VFUNC',`m4_format(`LAL%sReadTVectorSeries',TYPECODE)')m4_dnl
m4_define(`AFUNC',`m4_format(`LAL%sReadTArraySeries',TYPECODE)')m4_dnl
m4_define(`FFUNC',`m4_format(`LAL%sReadFSeries',TYPECODE)')m4_dnl
m4_define(`SCREATE',`m4_format(`LAL%sCreateVector',TYPECODE)')m4_dnl
m4_define(`VCREATE',`m4_format(`LAL%sCreateVectorSequence',TYPECODE)')m4_dnl
m4_define(`ACREATE',`m4_format(`LAL%sCreateArraySequence',TYPECODE)')m4_dnl
m4_define(`FMT',`m4_format(`LAL_%s_FORMAT',DATA)')m4_dnl
m4_define(`STRINGTODATA',`m4_format(`LALStringTo%s',DATACODE)')m4_dnl
m4_dnl
<lalVerbatim file="StreamSeriesInputCP"> */
void
SFUNC ( LALStatus *stat, STYPE *series, FILE *stream )
{ /* </lalVerbatim> */
  BufferList head;   /* head of linked list of buffers */
  BufferList *here;  /* pointer to current position in list */
  CHARVector *line = NULL; /* current line being read */
  CHAR *start, *end; /* start and end of a token on a line */
  UINT4 length = 0;  /* number of sequence elements to be read */
  DATA *data;        /* pointer to data in buffers */
  TYPE *sData;       /* pointer to data in output sequence */
  STYPE sCopy; /* internal copy of series */
  size_t nTot = 0;   /* total number of values read */
  size_t nMax, n;    /* # of values in buffer, and countdown index */
  int numRead = 1;   /* number of values read per call of fscanf() */

  INITSTATUS( stat, "SFUNC", STREAMSERIESINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( series, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !(series->data), stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  sCopy = *series;

  /*******************************************************************
   * PARSE METADATA HEADER                                           *
   *******************************************************************/

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
	endValue = startValue + 1;
	while ( *endValue != '"' && *endValue != '\0' )
	  endValue++;
	if ( *startValue != '"' || *endValue != '"' ||
	     endValue - ++startValue >= LALNameLength )
	  LALWarning( stat, LALREADSERIESC_HEADER "name" );
	else {
	  if ( endValue - startValue > 0 )
	    memcpy( sCopy.name, startValue,
		    ( endValue - startValue )*sizeof(CHAR) );
	  sCopy.name[endValue-startValue] = '\0';
	}
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

  /* Prepare data buffer. */
  here = &head;
  here->next = NULL;
  data = here->buf.DATACODE;
  end = start;
  start = NULL;
  nMax = BUFFSIZE/SIZE;
`#'if COMPLEX
  length *= 2;
#endif

  /* If length was specified, read data until it is reached.  Start
     with the line currently in memory. */
  if ( length > 0 ) {
    if ( nMax > length - nTot )
      nMax = length - nTot;
    n = nMax;
    while ( end != start && nTot < length ) {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n && nMax + nTot < length ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nTot += nMax;
	nMax = BUFFSIZE/SIZE;
	if ( nMax > length - nTot )
	  nMax = length - nTot;
	n = nMax;
      } else {
	nTot += nMax - n;
	data--;
      }
    }

    /* Read the remaining data using fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    while ( nTot < length ) {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      if ( numRead != 1 ) {
	FREEBUFFERLIST( head.next );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      nTot += nMax - n;
      if ( nTot < length ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nMax = BUFFSIZE/SIZE;
	if ( nMax > length - nTot )
	  nMax = length - nTot;
	n = nMax;
      }
    }
  }

  /* Otherwise, just read to the end of the file.  Start with the line
     currently in memory. */
  else {
    n = nMax;
    while ( end != start ) {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nTot += nMax;
	n = nMax = BUFFSIZE/SIZE;
      } else {
	nTot += nMax - n;
	data--;
      }
    }

    /* Read the remaining data using fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    while ( numRead == 1 ) {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      nTot += nMax - n;
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	nMax = BUFFSIZE/SIZE;
	data = here->buf.DATACODE;
	n = nMax;
      }
    }
  }


  /*******************************************************************
   * STORE DATA                                                      *
   *******************************************************************/

  /* Create the sequence to store the data. */
  length = nTot;
`#'if COMPLEX
  length /= 2;
  nTot = length*2;
#endif
  if ( !length ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
  }
  SCREATE ( stat->statusPtr, &(sCopy.data), length );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Start with head of list. */
  here = &head;
  data = here->buf.DATACODE;
  sData = sCopy.data->data;
  nMax = BUFFSIZE/SIZE;
  if ( nMax > nTot )
    nMax = nTot;
  while ( nTot ) {
`#'if COMPLEX
    /* For complex types, data need to be copied pairwise. */
    n = nMax/2;
    nTot -= 2*n;
    while ( n-- ) {
      sData->re = *(data++);
      (sData++)->im = *(data++);
    }
    if ( 2*(nMax/2) < nMax ) {
      /* Buffer boundary splits real and imaginary parts. */
      sData->re = *data;
      here = here->next;
      data = here->buf.DATACODE;
      (sData++)->im = *(data++);
      nTot -= 2;
      nMax = BUFFSIZE/SIZE - 1;
      if ( nMax > nTot )
	nMax = nTot;
    } else {
      /* Buffer boundary doesn't split real and imaginary parts. */
    here = here->next;
    data = here->buf.DATACODE;
    nMax = BUFFSIZE/SIZE;
    if ( nMax > nTot )
      nMax = nTot;
  }
#else
    /* For base types, data can be copied directly. */
    memcpy( sData, here->buf.DATACODE, nMax*sizeof( TYPE ) );
    sData += nMax;
    nTot -= nMax;
    here = here->next;
    nMax = BUFFSIZE/SIZE;
    if ( nMax > nTot )
      nMax = nTot;
#endif
  }

  /* Data have been stored successfully.  So, clean up and exit. */
  *series = sCopy;
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}

/*<lalVerbatim file="StreamSeriesInputCP"> */
void
VFUNC ( LALStatus *stat, VTYPE *series, FILE *stream )
{ /* </lalVerbatim> */
  BufferList head;         /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHARVector *line = NULL; /* current line being read */
  CHAR *start, *end;       /* start and end of a token on a line */
  UINT4 length = 0;        /* number of sequence elements to be read */
  UINT4 vectorLength = 0;  /* number of components per element */
  DATA *data;              /* pointer to data in buffers */
  TYPE *sData;             /* pointer to data in output sequence */
  VTYPE sCopy; /* internal copy of series */
  CreateVectorSequenceIn in; /* structure to create sequence */
  size_t nTot = 0; /* total number of values read */
  size_t nMax, n;  /* # of values in buffer, and countdown index */
  int numRead = 1; /* number of values read per call of fscanf() */

  INITSTATUS( stat, "VFUNC", STREAMSERIESINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( series, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !(series->data), stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  sCopy = *series;

  /*******************************************************************
   * PARSE METADATA HEADER                                           *
   *******************************************************************/

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
	endValue = startValue + 1;
	while ( *endValue != '"' && *endValue != '\0' )
	  endValue++;
	if ( *startValue != '"' || *endValue != '"' ||
	     endValue - ++startValue >= LALNameLength )
	  LALWarning( stat, LALREADSERIESC_HEADER "name" );
	else {
	  if ( endValue - startValue > 0 )
	    memcpy( sCopy.name, startValue,
		    ( endValue - startValue )*sizeof(CHAR) );
	  sCopy.name[endValue-startValue] = '\0';
	}
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

  /* Prepare data buffer. */
  here = &head;
  here->next = NULL;
  data = here->buf.DATACODE;
  end = start;
  start = NULL;
  nMax = BUFFSIZE/SIZE;
`#'if COMPLEX
  vectorLength *= 2;
#endif

  /* If length and vectorLength were both specified, we know we will
     need length*vectorLength data, so this will limit the number of
     data to be read even on the first line currently in memory. */
  if ( length > 0 && vectorLength > 0 ) {
    length *= vectorLength;
    if ( nMax > length - nTot )
      nMax = length - nTot;
    n = nMax;
    while ( end != start && nTot < length ) {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n && nMax + nTot < length ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nTot += nMax;
	nMax = BUFFSIZE/SIZE;
	if ( nMax > length - nTot )
	  nMax = length - nTot;
	n = nMax;
      } else {
	nTot += nMax - n;
	data--;
      }
    }
  }

  /* Otherwise, we will have to read at least the entire first line
     currently in memory. */
  else {
    n = nMax;
    while ( end != start ) {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nTot += nMax;
	n = nMax = BUFFSIZE/SIZE;
      } else {
	nTot += nMax - n;
	data--;
      }
    }
  }

  /* We can now compute vectorLength, if it wasn't specified. */
  if ( vectorLength == 0 ) {
    vectorLength = nTot + nMax - n;
`#'if COMPLEX
    if ( 2*(vectorLength/2) < vectorLength ) {
      FREEBUFFERLIST( head.next );
      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
      ABORT( stat, STREAMINPUTH_EVLEN, STREAMINPUTH_MSGEVLEN );
    }
#endif
    length *= vectorLength;
    if ( length > 0 && nMax > length - nTot ) {
      n -= nMax - length + nTot;
      nMax = length - nTot;
    }
  }

  /* Now, read the rest of the data using fscanf().  If length was
     specified, read only the required amound. */
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
  if ( length > 0 ) {
    while ( nTot < length ) {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      if ( numRead != 1 ) {
	FREEBUFFERLIST( head.next );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      nTot += nMax - n;
      if ( nTot < length ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nMax = BUFFSIZE/SIZE;
	if ( nMax > length - nTot )
	  nMax = length - nTot;
	n = nMax;
      }
    }
  }

  /* Otherwise, the total amount of data to be read is unknown, so
     just read to the end of the file. */
  else {
    while ( numRead == 1 ) {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      nTot += nMax - n;
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	nMax = BUFFSIZE/SIZE;
	data = here->buf.DATACODE;
	n = nMax;
      }
    }
  }


  /*******************************************************************
   * STORE DATA                                                      *
   *******************************************************************/

  /* Create the sequence to store the data. */
  length = nTot/vectorLength;
  nTot = length*vectorLength;
`#'if COMPLEX
  vectorLength /= 2;
#endif
  if ( !length ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
  }
  in.length = length;
  in.vectorLength = vectorLength;
  VCREATE ( stat->statusPtr, &(sCopy.data), &in );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Start with head of list. */
  here = &head;
  data = here->buf.DATACODE;
  sData = sCopy.data->data;
  nMax = BUFFSIZE/SIZE;
  if ( nMax > nTot )
    nMax = nTot;
  while ( nTot ) {
`#'if COMPLEX
    /* For complex types, data need to be copied pairwise. */
    n = nMax/2;
    nTot -= 2*n;
    while ( n-- ) {
      sData->re = *(data++);
      (sData++)->im = *(data++);
    }
    if ( 2*(nMax/2) < nMax ) {
      /* Buffer boundary splits real and imaginary parts. */
      sData->re = *data;
      here = here->next;
      data = here->buf.DATACODE;
      (sData++)->im = *(data++);
      nTot -= 2;
      nMax = BUFFSIZE/SIZE - 1;
      if ( nMax > nTot )
	nMax = nTot;
    } else {
      /* Buffer boundary doesn't split real and imaginary parts. */
    here = here->next;
    data = here->buf.DATACODE;
    nMax = BUFFSIZE/SIZE;
    if ( nMax > nTot )
      nMax = nTot;
  }
#else
    /* For base types, data can be copied directly. */
    memcpy( sData, here->buf.DATACODE, nMax*sizeof( TYPE ) );
    sData += nMax;
    nTot -= nMax;
    here = here->next;
    nMax = BUFFSIZE/SIZE;
    if ( nMax > nTot )
      nMax = nTot;
#endif
  }

  /* Data have been stored successfully.  So, clean up and exit. */
  *series = sCopy;
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}

/*<lalVerbatim file="StreamSeriesInputCP"> */
void
AFUNC ( LALStatus *stat, ATYPE *series, FILE *stream )
{ /* </lalVerbatim> */
  BufferList head;         /* head of linked list of buffers */
  BufferList *here;        /* pointer to current position in list */
  CHARVector *line = NULL; /* current line being read */
  CHAR *start, *end;       /* start and end of a token on a line */
  UINT4 length = 0;        /* number of sequence elements to be read */
  UINT4 arrayDim = 0;      /* number of components per element */
  UINT4Vector *dimLength = NULL; /* number of components per index */
  DATA *data;              /* pointer to data in buffers */
  TYPE *sData;             /* pointer to data in output sequence */
  ATYPE sCopy; /* internal copy of series */
  CreateArraySequenceIn in; /* structure to create sequence */
  size_t nTot = 0; /* total number of values read */
  size_t nMax, n;  /* # of values in buffer, and countdown index */
  int numRead = 1; /* number of values read per call of fscanf() */

  INITSTATUS( stat, "AFUNC", STREAMSERIESINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( series, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !(series->data), stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  sCopy = *series;

  /*******************************************************************
   * PARSE METADATA HEADER                                           *
   *******************************************************************/

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
	endValue = startValue + 1;
	while ( *endValue != '"' && *endValue != '\0' )
	  endValue++;
	if ( *startValue != '"' || *endValue != '"' ||
	     endValue - ++startValue >= LALNameLength )
	  LALWarning( stat, LALREADSERIESC_HEADER "name" );
	else {
	  if ( endValue - startValue > 0 )
	    memcpy( sCopy.name, startValue,
		    ( endValue - startValue )*sizeof(CHAR) );
	  sCopy.name[endValue-startValue] = '\0';
	}
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
	UINT4 *dData = head.buf.U4;
	here = &head;
	here->next = NULL;
	n = BUFFSIZE/4;
	endValue = startValue;
	startValue = NULL;
	/* Read components into the buffer list. */
	while ( endValue != startValue ) {
	  do {
	    LALStringToU4( stat->statusPtr, dData++,
			   startValue = endValue, &endValue );
	    BEGINFAIL( stat ) {
	      FREEBUFFERLIST( head.next );
	      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	      if ( dimLength ) {
		TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ),
		     stat );
	      }
	    } ENDFAIL( stat );
	  } while ( endValue != startValue && --n );
	  if ( !n ) {
	    here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	    if ( !(here->next) ) {
	      FREEBUFFERLIST( head.next );
	      TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	      if ( dimLength ) {
		TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ),
		     stat );
	      }
	      ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	    }
	    here = here->next;
	    here->next = NULL;
	    dData = here->buf.U4;
	    nTot += BUFFSIZE/4;
	    n = BUFFSIZE/4;
	  } else
	    nTot += BUFFSIZE/4 - n;
	}
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
	    FREEBUFFERLIST( head.next );
	    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  } ENDFAIL( stat );
	  here = &head;
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
	FREEBUFFERLIST( head.next );
	head.next = NULL;
	nTot = 0;
      }

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

  /* Prepare data buffer. */
  here = &head;
  here->next = NULL;
  data = here->buf.DATACODE;
  end = start;
  start = NULL;
  nMax = BUFFSIZE/SIZE;
`#'if COMPLEX
  arrayDim *= 2;
#endif
  length *= arrayDim;

  /* If length was specified, read data until it is reached.  Start
     with the line currently in memory. */
  if ( length > 0 ) {
    if ( nMax > length - nTot )
      nMax = length - nTot;
    n = nMax;
    while ( end != start && nTot < length ) {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n && nMax + nTot < length ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nTot += nMax;
	nMax = BUFFSIZE/SIZE;
	if ( nMax > length - nTot )
	  nMax = length - nTot;
	n = nMax;
      } else {
	nTot += nMax - n;
	data--;
      }
    }

    /* Read the remaining data using fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    while ( nTot < length ) {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      if ( numRead != 1 ) {
	FREEBUFFERLIST( head.next );
	TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      nTot += nMax - n;
      if ( nTot < length ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nMax = BUFFSIZE/SIZE;
	if ( nMax > length - nTot )
	  nMax = length - nTot;
	n = nMax;
      }
    }
  }

  /* Otherwise, just read to the end of the file.  Start with the line
     currently in memory. */
  else {
    n = nMax;
    while ( end != start ) {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nTot += nMax;
	n = nMax = BUFFSIZE/SIZE;
      } else {
	nTot += nMax - n;
	data--;
      }
    }

    /* Read the remaining data using fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    while ( numRead == 1 ) {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      nTot += nMax - n;
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	nMax = BUFFSIZE/SIZE;
	data = here->buf.DATACODE;
	n = nMax;
      }
    }
  }


  /*******************************************************************
   * STORE DATA                                                      *
   *******************************************************************/

  /* Create the sequence to store the data. */
  length = nTot/arrayDim;
  nTot = length*arrayDim;
`#'if COMPLEX
  arrayDim /= 2;
#endif
  if ( !length ) {
    FREEBUFFERLIST( head.next );
    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
    ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
  }
  in.length = length;
  in.dimLength = dimLength;
  ACREATE ( stat->statusPtr, &(sCopy.data), &in );
  BEGINFAIL( stat ) {
    FREEBUFFERLIST( head.next );
    TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );
  } ENDFAIL( stat );
  TRY( LALU4DestroyVector( stat->statusPtr, &dimLength ), stat );

  /* Start with head of list. */
  here = &head;
  data = here->buf.DATACODE;
  sData = sCopy.data->data;
  nMax = BUFFSIZE/SIZE;
  if ( nMax > nTot )
    nMax = nTot;
  while ( nTot ) {
`#'if COMPLEX
    /* For complex types, data need to be copied pairwise. */
    n = nMax/2;
    nTot -= 2*n;
    while ( n-- ) {
      sData->re = *(data++);
      (sData++)->im = *(data++);
    }
    if ( 2*(nMax/2) < nMax ) {
      /* Buffer boundary splits real and imaginary parts. */
      sData->re = *data;
      here = here->next;
      data = here->buf.DATACODE;
      (sData++)->im = *(data++);
      nTot -= 2;
      nMax = BUFFSIZE/SIZE - 1;
      if ( nMax > nTot )
	nMax = nTot;
    } else {
      /* Buffer boundary doesn't split real and imaginary parts. */
    here = here->next;
    data = here->buf.DATACODE;
    nMax = BUFFSIZE/SIZE;
    if ( nMax > nTot )
      nMax = nTot;
  }
#else
    /* For base types, data can be copied directly. */
    memcpy( sData, here->buf.DATACODE, nMax*sizeof( TYPE ) );
    sData += nMax;
    nTot -= nMax;
    here = here->next;
    nMax = BUFFSIZE/SIZE;
    if ( nMax > nTot )
      nMax = nTot;
#endif
  }

  /* Data have been stored successfully.  So, clean up and exit. */
  *series = sCopy;
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}

/* <lalVerbatim file="StreamSeriesInputCP"> */
void
FFUNC ( LALStatus *stat, FTYPE *series, FILE *stream )
{ /* </lalVerbatim> */
  BufferList head;   /* head of linked list of buffers */
  BufferList *here;  /* pointer to current position in list */
  CHARVector *line = NULL; /* current line being read */
  CHAR *start, *end; /* start and end of a token on a line */
  UINT4 length = 0;  /* number of sequence elements to be read */
  DATA *data;        /* pointer to data in buffers */
  TYPE *sData;       /* pointer to data in output sequence */
  FTYPE sCopy; /* internal copy of series */
  size_t nTot = 0;   /* total number of values read */
  size_t nMax, n;    /* # of values in buffer, and countdown index */
  int numRead = 1;   /* number of values read per call of fscanf() */

  INITSTATUS( stat, "FFUNC", STREAMSERIESINPUTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( stream, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( series, stat, STREAMINPUTH_ENUL, STREAMINPUTH_MSGENUL );
  ASSERT( !(series->data), stat, STREAMINPUTH_EOUT, STREAMINPUTH_MSGEOUT );

  sCopy = *series;

  /*******************************************************************
   * PARSE METADATA HEADER                                           *
   *******************************************************************/

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
	endValue = startValue + 1;
	while ( *endValue != '"' && *endValue != '\0' )
	  endValue++;
	if ( *startValue != '"' || *endValue != '"' ||
	     endValue - ++startValue >= LALNameLength )
	  LALWarning( stat, LALREADSERIESC_HEADER "name" );
	else {
	  if ( endValue - startValue > 0 )
	    memcpy( sCopy.name, startValue,
		    ( endValue - startValue )*sizeof(CHAR) );
	  sCopy.name[endValue-startValue] = '\0';
	}
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

  /* Prepare data buffer. */
  here = &head;
  here->next = NULL;
  data = here->buf.DATACODE;
  end = start;
  start = NULL;
  nMax = BUFFSIZE/SIZE;
`#'if COMPLEX
  length *= 2;
#endif

  /* If length was specified, read data until it is reached.  Start
     with the line currently in memory. */
  if ( length > 0 ) {
    if ( nMax > length - nTot )
      nMax = length - nTot;
    n = nMax;
    while ( end != start && nTot < length ) {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n && nMax + nTot < length ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nTot += nMax;
	nMax = BUFFSIZE/SIZE;
	if ( nMax > length - nTot )
	  nMax = length - nTot;
	n = nMax;
      } else {
	nTot += nMax - n;
	data--;
      }
    }

    /* Read the remaining data using fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    while ( nTot < length ) {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      if ( numRead != 1 ) {
	FREEBUFFERLIST( head.next );
	ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
      }
      nTot += nMax - n;
      if ( nTot < length ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nMax = BUFFSIZE/SIZE;
	if ( nMax > length - nTot )
	  nMax = length - nTot;
	n = nMax;
      }
    }
  }

  /* Otherwise, just read to the end of the file.  Start with the line
     currently in memory. */
  else {
    n = nMax;
    while ( end != start ) {
      do {
	STRINGTODATA ( stat->statusPtr, data++, start = end, &end );
	BEGINFAIL( stat ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	} ENDFAIL( stat );
      } while ( end != start && --n );
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	data = here->buf.DATACODE;
	nTot += nMax;
	n = nMax = BUFFSIZE/SIZE;
      } else {
	nTot += nMax - n;
	data--;
      }
    }

    /* Read the remaining data using fscanf(). */
    TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );
    while ( numRead == 1 ) {
      do
	numRead = fscanf( stream, "%" FMT, data++ );
      while ( numRead == 1 && --n );
      nTot += nMax - n;
      if ( !n ) {
	here->next = (BufferList *)LALMalloc( sizeof(BufferList) );
	if ( !(here->next) ) {
	  FREEBUFFERLIST( head.next );
	  ABORT( stat, STREAMINPUTH_EMEM, STREAMINPUTH_MSGEMEM );
	}
	here = here->next;
	here->next = NULL;
	nMax = BUFFSIZE/SIZE;
	data = here->buf.DATACODE;
	n = nMax;
      }
    }
  }


  /*******************************************************************
   * STORE DATA                                                      *
   *******************************************************************/

  /* Create the sequence to store the data. */
  length = nTot;
`#'if COMPLEX
  length /= 2;
  nTot = length*2;
#endif
  if ( !length ) {
    FREEBUFFERLIST( head.next );
    ABORT( stat, STREAMINPUTH_ESLEN, STREAMINPUTH_MSGESLEN );
  }
  SCREATE ( stat->statusPtr, &(sCopy.data), length );
  BEGINFAIL( stat )
    FREEBUFFERLIST( head.next );
  ENDFAIL( stat );

  /* Start with head of list. */
  here = &head;
  data = here->buf.DATACODE;
  sData = sCopy.data->data;
  nMax = BUFFSIZE/SIZE;
  if ( nMax > nTot )
    nMax = nTot;
  while ( nTot ) {
`#'if COMPLEX
    /* For complex types, data need to be copied pairwise. */
    n = nMax/2;
    nTot -= 2*n;
    while ( n-- ) {
      sData->re = *(data++);
      (sData++)->im = *(data++);
    }
    if ( 2*(nMax/2) < nMax ) {
      /* Buffer boundary splits real and imaginary parts. */
      sData->re = *data;
      here = here->next;
      data = here->buf.DATACODE;
      (sData++)->im = *(data++);
      nTot -= 2;
      nMax = BUFFSIZE/SIZE - 1;
      if ( nMax > nTot )
	nMax = nTot;
    } else {
      /* Buffer boundary doesn't split real and imaginary parts. */
    here = here->next;
    data = here->buf.DATACODE;
    nMax = BUFFSIZE/SIZE;
    if ( nMax > nTot )
      nMax = nTot;
  }
#else
    /* For base types, data can be copied directly. */
    memcpy( sData, here->buf.DATACODE, nMax*sizeof( TYPE ) );
    sData += nMax;
    nTot -= nMax;
    here = here->next;
    nMax = BUFFSIZE/SIZE;
    if ( nMax > nTot )
      nMax = nTot;
#endif
  }

  /* Data have been stored successfully.  So, clean up and exit. */
  *series = sCopy;
  FREEBUFFERLIST( head.next );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
