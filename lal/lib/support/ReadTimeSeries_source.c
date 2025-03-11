#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,TimeSeries)
#define VTYPE CONCAT2(TYPE,Sequence)
#define FUNC CONCAT3(LAL,TYPECODE,ReadTimeSeries)

#ifndef BASETYPE
#define BASETYPE TYPE
#endif


/* Maybe for consistent allocation of memory we should include a
companion function called element counter that will notify the calling
routine of how large we should allocate the Series to be before
calling the read frequency series module */


void
FUNC ( LALStatus* status,
                  STYPE *series,
                  const CHAR *filename )

{
  REAL8Vector	*t=NULL;
  REAL8         *tPtr;
  REAL8         *tStopPtr;
  union { TYPE value; BASETYPE array[sizeof(TYPE)/sizeof(BASETYPE)]; } data;
  TYPE	        *outputPtr;
  FILE		*fp;
  CHAR		 line[MaxLineLength];  /*holds data from each line*/
  LALUnit        tempUnit;
  CHARVector    *string=NULL;

  /* need to declare error section here */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT( filename != NULL,status,READFTSERIESH_EFILENOTFOUND,
	  READFTSERIESH_MSGEFILENOTFOUND );

  /* if (filename == NULL) return; */

  fp = LALFopen( filename, "r" );
  if (fp == NULL)
  {
    ABORT( status, READFTSERIESH_EFILENOTFOUND,
	   READFTSERIESH_MSGEFILENOTFOUND );
  }

  /* limited to line of data not exceeding MaxLineLength chars*/

  if (fgets(line,sizeof(line),fp) == NULL)
  {
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }

  if (line[0] != '#' || line[1] != ' ')
  {
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }

  if (snprintf( series->name, sizeof( series->name ), "%s", line + 2) < 0)
  {
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }
  /* Change trailing linefeed to '\0' */
  if ( changeCharToNull(series->name, '\n', series->name + LALNameLength) )
  {
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }

  if (fgets(line,sizeof(line),fp) == NULL)
  {
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }
  if (line[0] != '#' || line[1] != ' ')
  {
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }
  if (line[2] == '\n')
  {
    series->f0 = 0.0;
  }
  else if ( sscanf( line, "# Heterodyned at %lf Hz\n", &(series->f0) ) != 1 )
  {
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }

  if (fgets(line,sizeof(line),fp) == NULL)
  {
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }

  if (line[0] != '#' || line[1] != ' ')
  {
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }
  if (line[2] == '\n')
  {
    series->epoch.gpsSeconds = 0;
    series->epoch.gpsNanoSeconds = 0;
  }
  else if ( sscanf( line, "# Epoch is %d seconds, %d nanoseconds\n",
	            &(series->epoch.gpsSeconds),
                    &(series->epoch.gpsNanoSeconds) )
	    != 2 )
  {
    ABORT( status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE );
  }

  TRY( LALCHARCreateVector(status->statusPtr, &string, MaxLineLength), status );

  if (fgets(line,sizeof(line),fp) == NULL)
  {
    TRY( LALCHARDestroyVector( status->statusPtr, &string ), status );
    ABORT( status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE );
  }

  if (!strcmp(line,"# Units are ()\n"))
  {
    series->sampleUnits = lalDimensionlessUnit;
  }
  else {
    if ( sscanf( line, "# Units are (%[^)]", string->data ) != 1 )
    {
      TRY( LALCHARDestroyVector( status->statusPtr, &string ), status );
      ABORT( status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE );
    }
    if ( XLALParseUnitString(&tempUnit, string->data) == NULL )
    {
      TRY( LALCHARDestroyVector( status->statusPtr, &string ), status );
      ABORT( status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE );
    }
    series->sampleUnits = tempUnit;
  }

  TRY( LALCHARDestroyVector( status->statusPtr, &string ), status );

  TRY( LALDCreateVector(status->statusPtr, &t, series->data->length),
       status );


  tPtr = &(t->data[0]);
  tStopPtr = tPtr + t->length;
  outputPtr = &(series->data->data[0]);

  if(fgets(line,sizeof(line),fp) == NULL)
  {
    TRY( LALDDestroyVector( status->statusPtr, &t ), status );
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }

  while(fgets(line,sizeof(line),fp)!=NULL)
  {
    /*change arg so we dereference pointer */
    if ( sscanf(line, FMT, tPtr, ARG) != 1 + NARGS )
    {
      TRY( LALDDestroyVector( status->statusPtr, &t ), status );
      ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
    }
    *(outputPtr) = data.value;
    tPtr++;
    outputPtr++;

    if (tPtr > tStopPtr)
    {
      TRY( LALDDestroyVector( status->statusPtr, &t ), status );
      ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
    }
  }
  if (tPtr != tStopPtr)
  {
    TRY( LALDDestroyVector( status->statusPtr, &t ), status );
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }

  (series->deltaT) = ( t->data[1] - t->data[0] );

  TRY( LALDDestroyVector( status->statusPtr, &t ), status );

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

#undef BASETYPE
