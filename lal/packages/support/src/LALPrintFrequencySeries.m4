dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')') dnl
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')') dnl
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')') dnl
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')') dnl
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')') dnl
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')') dnl
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')') dnl
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')') dnl
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')') dnl
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')') dnl
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')') dnl
define(`STYPE',`format(`%sFrequencySeries',TYPE)') dnl
define(`VTYPE',`format(`%sSequence',TYPE)') dnl
define(`FUNC',`format(`LAL%sPrintFrequencySeries',TYPECODE)') dnl
define(`FMT',`"%g %i\n"') dnl
ifelse(TYPECODE,`Z',`define(`FMT',`"%g\t%g\t%g\n"')') dnl
ifelse(TYPECODE,`C',`define(`FMT',`"%g\t%g\t%g\n"')') dnl
ifelse(TYPECODE,`D',`define(`FMT',`"%g\t%g\n"')') dnl
ifelse(TYPECODE,`S',`define(`FMT',`"%g\t%g\n"')') dnl
ifelse(TYPECODE,`I8',`define(`FMT',`"%g\t%0.0f\n"')') dnl
ifelse(TYPECODE,`U8',`define(`FMT',`"%g\t%0.0f\n"')') dnl
ifelse(TYPECODE,`',`define(`FMT',`"%g\t%f\n"')') dnl
define(`HEADER',`"# Freq (Hz)\tValue\n"') dnl
ifelse(TYPECODE,`Z',`define(`HEADER',`"# Freq (Hz)\tRe(Value)\tIm(Value)\n"')') dnl
ifelse(TYPECODE,`C',`define(`HEADER',`"# Freq (Hz)\tRe(Value)\tIm(Value)\n"')') dnl
define(`ARG',`*data') dnl
ifelse(TYPECODE,`Z',`define(`ARG',`data->re,data->im')') dnl
ifelse(TYPECODE,`C',`define(`ARG',`data->re,data->im')') dnl
ifelse(TYPECODE,`I8',`define(`ARG',`(REAL8) *data')') dnl
ifelse(TYPECODE,`U8',`define(`ARG',`(REAL8) *data')') dnl

/* <lalVerbatim file="PrintFrequencySeriesCP"> */
void FUNC ( STYPE *series, const CHAR *filename ) 
 /* </lalVerbatim> */
{
  REAL8 f;
  TYPE *data;
  FILE *fp;
  UINT4 i;
  static LALStatus status;
  CHARVector *unitString;

  if (series==NULL) return;

  /* *(series->data) is a VTYPE */
  /* series->data->data is a pointer to TYPE */

  /* open output file */
  fp=LALFopen(filename,"w");
  fprintf(fp,"# %s\n",series->name);
  if (series->epoch.gpsSeconds && series->epoch.gpsNanoSeconds) {
    fprintf(fp,"# Epoch is %d seconds, %d nanoseconds\n",
            series->epoch.gpsSeconds,series->epoch.gpsNanoSeconds);
  }
  else {
    fprintf(fp,"# \n");
  }
  unitString = NULL;
  LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
  LALUnitAsString(&status, unitString, &(series->sampleUnits));
  fprintf(fp,"# Units are (%s)\n",unitString->data);
  fprintf(fp,HEADER);
  LALCHARDestroyVector(&status, &unitString);
  for ( i = 0; i < series->data->length; ++i )
  {
    f = series->f0 + i * series->deltaF;
    data = &(series->data->data[i]);
    fprintf(fp,FMT,f,ARG);
  }

  LALFclose(fp);

  return;
}
