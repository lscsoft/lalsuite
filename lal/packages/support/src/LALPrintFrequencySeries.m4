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
define(`STYPE',`format(`%sFrequencySeries',TYPE)')
define(`VTYPE',`format(`%sSequence',TYPE)')
define(`FUNC',`format(`LAL%sPrintFrequencySeries',TYPECODE)')
define(`FMT',`"%g %i\n"')
ifelse(TYPECODE,`Z',`define(`FMT',`"%g %g %g\n"')')
ifelse(TYPECODE,`C',`define(`FMT',`"%g %g %g\n"')')
ifelse(TYPECODE,`D',`define(`FMT',`"%g %g\n"')')
ifelse(TYPECODE,`S',`define(`FMT',`"%g %g\n"')')
ifelse(TYPECODE,`I8',`define(`FMT',`"%g %0.0f\n"')')
ifelse(TYPECODE,`U8',`define(`FMT',`"%g %0.0f\n"')')
ifelse(TYPECODE,`',`define(`FMT',`"%g %f\n"')')
define(`ARG',`*data')
ifelse(TYPECODE,`Z',`define(`ARG',`data->re,data->im')')
ifelse(TYPECODE,`C',`define(`ARG',`data->re,data->im')')
ifelse(TYPECODE,`I8',`define(`ARG',`(REAL8) *data')')
ifelse(TYPECODE,`U8',`define(`ARG',`(REAL8) *data')')

/* <lalVerbatim file="PrintFrequencySeriesCP"> */
void FUNC ( STYPE *series, const CHAR *filename ) 
{ /* </lalVerbatim> */
  REAL8 f;
  TYPE *data;
  FILE *fp;
  TYPE *endOfSeq;

  if (series==NULL) return;

  /* *(series->data) is a VTYPE */
  /* series->data->data is a pointer to TYPE */

  /* Make a TYPE pointer which points to the first memory address not
   * belonging to the sequence
   */

  endOfSeq = series->data->data + series->data->length;

  /* open output file */
  fp=fopen(filename,"w");
  for ( f = series->f0, data = series->data->data;
        data < endOfSeq;
        f += series->deltaF, ++data )
  {
    fprintf(fp,FMT,f,ARG);
  }	

  fclose(fp);

  return;
}
