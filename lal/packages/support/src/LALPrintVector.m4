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
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`VTYPE',`format(`%sVector',TYPE)')
define(`FUNC',`format(`LAL%sPrintVector',TYPECODE)')
define(`FMT',`"%i %i\n"')
ifelse(TYPECODE,`Z',`define(`FMT',`"%i %g %g\n"')')
ifelse(TYPECODE,`C',`define(`FMT',`"%i %g %g\n"')')
ifelse(TYPECODE,`D',`define(`FMT',`"%i %g\n"')')
ifelse(TYPECODE,`S',`define(`FMT',`"%i %g\n"')')
ifelse(TYPECODE,`I8',`define(`FMT',`"%i %0.0f\n"')')
ifelse(TYPECODE,`U8',`define(`FMT',`"%i %0.0f\n"')')
ifelse(TYPECODE,`',`define(`FMT',`"%i %f\n"')')
ifelse(TYPECODE,`CHAR',`define(`FMT',`"%i %c\n"')')
define(`ARG',`vector->data[i]')
ifelse(TYPECODE,`Z',`define(`ARG',`vector->data[i].re,vector->data[i].im')')
ifelse(TYPECODE,`C',`define(`ARG',`vector->data[i].re,vector->data[i].im')')
ifelse(TYPECODE,`I8',`define(`ARG',`(REAL8) vector->data[i]')')
ifelse(TYPECODE,`U8',`define(`ARG',`(REAL8) vector->data[i]')')

/* <lalVerbatim file="PrintVectorCP"> */
void FUNC ( VTYPE *vector ) 
{ /* </lalVerbatim> */
  int i;
  static int fileno=0;
  FILE *fp;
  char fname[256];


  if (vector==NULL) return;

  /* open output file */
  sprintf(fname,"TYPECODE" "PrintVector.%03d",fileno++);
  fp=fopen(fname,"w");

  for (i=0;i<vector->length;i++)
    fprintf(fp,FMT,i,ARG);

  fclose(fp);

  return;
}
