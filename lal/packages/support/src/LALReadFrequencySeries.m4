ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`STYPE',`format(`%sFrequencySeries',TYPE)')
define(`VTYPE',`format(`%sSequence',TYPE)')
define(`FUNC',`format(`LAL%sReadFrequencySeries',TYPECODE)')
define(`FMT',`"%g %i\n"')
ifelse(TYPECODE,`Z',`define(`FMT',`"%lf\t%lf\t%lf\n"')')
ifelse(TYPECODE,`C',`define(`FMT',`"%lf\t%f\t%f\n"')')
ifelse(TYPECODE,`D',`define(`FMT',`"%lf\t%lf\n"')')
ifelse(TYPECODE,`S',`define(`FMT',`"%lf\t%f\n"')')
ifelse(TYPECODE,`',`define(`FMT',`"%lf\t%f\n"')')
define(`HEADER',`"Freq (Hz)\tValue\n"');
ifelse(TYPECODE,`Z',`define(`HEADER',`"Freq (Hz)\tRe(Value)\tIm(Value)\n"')')
ifelse(TYPECODE,`C',`define(`HEADER',`"Freq (Hz)\tRe(Value)\tIm(Value)\n"')')
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

/* <lalVerbatim file="ReadFrequencySeriesCP"> */
void FUNC ( LALStatus* status, 
                              STYPE *series,
	                      const CHAR *filename )
 /* </lalVerbatim> */
{
  REAL8Vector		*f=NULL;
  REAL8			*fPtr;
  REAL8			*fStopPtr;
  TYPE			data;
  TYPE			*outputPtr;
  FILE			*fp;
  CHAR			line[MaxLineLength];  /*holds data from each line*/
  LALUnit		tempUnit;
  CHARVector            *string=NULL;
  
  /* need to declare error section here */
  INITSTATUS(status, "ReadFrequencySeries", READFREQUENCYSERIESC);
  ATTATCHSTATUSPTR(status);

  ASSERT( filename != NULL, status, READFTSERIESH_EFILENOTFOUND,
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

  if (LALSnprintf( series->name, sizeof( series->name ), "%s", line ) < 0)
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
  if (line[0] == '\n') 
  {
    series->epoch.gpsSeconds = 0;
    series->epoch.gpsNanoSeconds = 0;
  }
  else if ( sscanf( line, "Epoch is %d seconds, %d nanoseconds\n",
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

  if (!strcmp(line,"Units are ()\n"))
  {
    series->sampleUnits = lalDimensionlessUnit; 
  }
  else {
    if ( sscanf( line, "Units are (%[^)]", string->data ) != 1 )
    {
      TRY( LALCHARDestroyVector( status->statusPtr, &string ), status );
      ABORT( status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE );
    }
    LALParseUnitString(status->statusPtr, &tempUnit, string);
    if (status->statusPtr->statusCode) 
    {
      TRY( LALCHARDestroyVector( status->statusPtr, &string ), status );
      ABORT( status, -1, "Recursive error" );
    }
    series->sampleUnits = tempUnit;
  }

  TRY( LALCHARDestroyVector( status->statusPtr, &string ), status );

  TRY( LALDCreateVector(status->statusPtr, &f, series->data->length),
       status );
  
  fPtr = &(f->data[0]);
  fStopPtr = fPtr + f->length;
  outputPtr = &(series->data->data[0]);
  
  if(fgets(line,sizeof(line),fp) == NULL) 
  {
    TRY( LALDDestroyVector( status->statusPtr, &f ), status );
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }
  
  while(fgets(line,sizeof(line),fp)!=NULL)
  {
    /*change arg so we dereference pointer */
    if ( sscanf(line, FMT, fPtr, ARG) != 1 + NARGS ) 
    {
      TRY( LALDDestroyVector( status->statusPtr, &f ), status );
      ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
    }
    *(outputPtr) = data;
    fPtr++;
    outputPtr++;	
    
    if (fPtr > fStopPtr) {
      TRY( LALDDestroyVector( status->statusPtr, &f ), status );
      ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
    }
  }
  if (fPtr != fStopPtr) 
  {
    TRY( LALDDestroyVector( status->statusPtr, &f ), status );
    ABORT(status, READFTSERIESH_EPARSE, READFTSERIESH_MSGEPARSE);
  }
  
  (series->deltaF) = ( f->data[1] - f->data[0] );
  (series->f0) = f->data[0];
  
  TRY( LALDDestroyVector( status->statusPtr, &f ), status );
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

