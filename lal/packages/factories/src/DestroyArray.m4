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
define(`ATYPE',`format(`%sArray',TYPE)')
define(`FUNC',`format(`LAL%sDestroyArray',TYPECODE)')

/* <lalVerbatim file="ArrayFactoriesD"> */
void FUNC ( LALStatus *status, ATYPE **array )
{ /* </lalVerbatim> */
  INITSTATUS (status, "FUNC", ARRAYFACTORIESC);	
  ATTATCHSTATUSPTR (status);
      
  ASSERT (array,          status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
  ASSERT (*array,         status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR);
  ASSERT ((*array)->data, status, AVFACTORIESH_EDPTR, AVFACTORIESH_MSGEDPTR);

  /* Free allocated storage */

  LALU4DestroyVector (status->statusPtr, &(*array)->dimLength);
  CHECKSTATUSPTR (status);

  LALFree ((*array)->data); /* free allocated data */
  LALFree (*array);	    /* free array struct itself */
  *array = NULL;	    /* make sure we don't point to freed struct */

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
