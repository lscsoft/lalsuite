dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`STYPE',`format(`%sTimeSeries',TYPE)')
define(`VTYPE',`format(`%sSequence',TYPE)')
define(`FUNC',`format(`LAL%sPrintTimeSeries',TYPECODE)')
define(`FMT',`"%g\t%i\n"')
ifelse(TYPECODE,`Z',`define(`FMT',`"%g\t%g\t%g\n"')')
ifelse(TYPECODE,`C',`define(`FMT',`"%g\t%g\t%g\n"')')
ifelse(TYPECODE,`D',`define(`FMT',`"%g\t%g\n"')')
ifelse(TYPECODE,`S',`define(`FMT',`"%g\t%g\n"')')
ifelse(TYPECODE,`I8',`define(`FMT',`"%g\t%0.0f\n"')')
ifelse(TYPECODE,`U8',`define(`FMT',`"%g\t%0.0f\n"')')
ifelse(TYPECODE,`',`define(`FMT',`"%g\t%f\n"')')
define(`HEADER',`"# Seconds since epoch\tValue\n"')
ifelse(TYPECODE,`Z',`define(`HEADER',`"# Seconds since epoch\tRe(Value)\tIm(Value)\n"')')
ifelse(TYPECODE,`C',`define(`HEADER',`"# Seconds since epoch\tRe(Value)\tIm(Value)\n"')')
define(`ARG',`*data')
ifelse(TYPECODE,`Z',`define(`ARG',`data->re,data->im')')
ifelse(TYPECODE,`C',`define(`ARG',`data->re,data->im')')
ifelse(TYPECODE,`I8',`define(`ARG',`(REAL8) *data')')
ifelse(TYPECODE,`U8',`define(`ARG',`(REAL8) *data')')

/* <lalVerbatim file="PrintTimeSeriesCP"> */
void FUNC ( STYPE *series, const CHAR *filename ) 
{ /* </lalVerbatim> */
  REAL8 t;
  TYPE *data;
  FILE *fp;
  UINT4 i;
  static LALStatus status;
  CHARVector *unitString;

  if (series==NULL) return;

  /* *(series->data) is a VTYPE */
  /* series->data->data is a pointer to TYPE */

  /* Make a TYPE pointer which points to the first memory address not
   * belonging to the sequence
   */

  /* open output file */
  fp=LALFopen(filename,"w");
  fprintf(fp,"# %s\n",series->name);
  if (series->f0) {
     fprintf(fp,"# Heterodyned at %g Hz\n",series->f0);
  }
  else {
    fprintf(fp,"# \n");
  }
  fprintf(fp,"# Epoch is %d seconds, %d nanoseconds\n",
          series->epoch.gpsSeconds,series->epoch.gpsNanoSeconds);
  unitString = NULL;
  LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
  LALUnitAsString(&status, unitString, &(series->sampleUnits));
  fprintf(fp,"# Units are (%s)\n",unitString->data);
  fprintf(fp,HEADER);
  LALCHARDestroyVector(&status, &unitString);
  for ( i = 0; i < series->data->length; ++i )
  {
    t = i * series->deltaT;
    data = &(series->data->data[i]);
    fprintf(fp,FMT,t,ARG);
  }	

  LALFclose(fp);

  return;
}
