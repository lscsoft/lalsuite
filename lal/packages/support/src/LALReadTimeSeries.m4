ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`STYPE',`format(`%sTimeSeries',TYPE)')
define(`VTYPE',`format(`%sSequence',TYPE)')
define(`FUNC',`format(`LAL%sReadTimeSeries',TYPECODE)')
define(`FMT',`"%g %i\n"')
ifelse(TYPECODE,`Z',`define(`FMT',`"%lf\t%lf\t%lf\n"')')
ifelse(TYPECODE,`C',`define(`FMT',`"%lf\t%f\t%f\n"')')
ifelse(TYPECODE,`D',`define(`FMT',`"%lf\t%lf\n"')')
ifelse(TYPECODE,`S',`define(`FMT',`"%lf\t%f\n"')')
ifelse(TYPECODE,`',`define(`FMT',`"%lf\t%f\n"')')
define(`ARG',`&data')
define(`NARGS',`1')
ifelse(TYPECODE,`Z',`define(`ARG',`&(data.re),&(data.im)')')
ifelse(TYPECODE,`C',`define(`ARG',`&(data.re),&(data.im)')')
ifelse(TYPECODE,`Z',`define(`NARGS',`2')')
ifelse(TYPECODE,`C',`define(`NARGS',`2')')
/* Maybe for consistent allocation of memory we should include a
companion function called element counter that will notify the calling
routine of how large we should allocate the Series to be before
calling the read frequency series module */

/* <lalVerbatim file="ReadTimeSeriesCP"> */
void 
FUNC ( LALStatus* status,
                  STYPE *series,
                  const CHAR *filename )
 /* </lalVerbatim> */
{
  REAL8Vector	*t=NULL;
  REAL8         *tPtr;
  REAL8         *tStopPtr;
  TYPE           data;
  TYPE	        *outputPtr;
  FILE		*fp;
  CHAR		 line[MaxLineLength];  /*holds data from each line*/
  LALUnit        tempUnit;
  CHARVector    *string=NULL;

  /* need to declare error section here */
  INITSTATUS(status, "ReadTimeSeries", READTIMESERIESC);
  ATTATCHSTATUSPTR(status);  
 
  ASSERT( filename != NULL,status,READFTSERIESH_EFILENOTFOUND,
	  READFTSERIESH_MSGEFILENOTFOUND );

  /* if (filename == NULL) return; */
  
  fp = LALOpenDataFile( filename );
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
    LALParseUnitString(status->statusPtr, &tempUnit, string);
    BEGINFAIL( status ) 
      TRY( LALCHARDestroyVector( status->statusPtr, &string ), status );
    ENDFAIL( status ); 
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
    *(outputPtr) = data;
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
