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
define(`FUNC',`format(`LAL%sCreateArray',TYPECODE)')

/* <lalVerbatim file="ArrayFactoriesD"> */
void FUNC ( LALStatus *status, ATYPE **array, UINT4Vector *dimLength ) 
{ /* </lalVerbatim> */
  UINT4 arrayDataSize = 1;
  UINT4 numDims;
  UINT4 dim;

  INITSTATUS (status, "FUNC", ARRAYFACTORIESC);   
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are sane */

  ASSERT (array,             status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
  ASSERT (!*array,           status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR);
  ASSERT (dimLength,         status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
  ASSERT (dimLength->data,   status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
  ASSERT (dimLength->length, status,
          AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH);

  numDims = dimLength->length;

  /* loop over dimensions to compute total size of array data */

  for (dim = 0; dim < numDims; ++dim)
  {
    arrayDataSize *= dimLength->data[dim];
  }

  ASSERT (arrayDataSize, status, AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH);

  /* allocate memory for structure */

  *array = ( ATYPE * ) LALMalloc ( sizeof( ATYPE ) );
  if ( NULL == *array )
  {
    ABORT( status, AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC );
  }

  (*array)->dimLength = NULL;
  (*array)->data      = NULL;

  /* allocate dimLength field and copy information there */

  LALU4CreateVector (status->statusPtr, &(*array)->dimLength, numDims);
  if (status->statusPtr->statusCode)
  {
    LALFree (*array);
    *array = NULL;
    RETURN (status);   /* this returns a recursive error status code (-1) */
  }
  memcpy ((*array)->dimLength->data, dimLength->data, numDims*sizeof(UINT4));
  
  /* allocate storage */

  (*array)->data = ( TYPE * ) LALMalloc( arrayDataSize*sizeof( TYPE ) );

  if (!(*array)->data)
  {
    /* try to free memory                                               */
    /* this should ALWAYS work so we can use the CHECKSTATUSPTR() macro */
    LALU4DestroyVector (status->statusPtr, &(*array)->dimLength);
    CHECKSTATUSPTR (status);  /* if this fails, there is a memory leak  */

    LALFree (*array);
    *array = NULL;
    ABORT (status, AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC);
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

