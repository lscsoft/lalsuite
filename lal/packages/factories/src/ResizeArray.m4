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
define(`ATYPE',`format(`%sArray',TYPE)')
define(`RESIZEARRAY',`format(`LAL%sResizeArray',TYPECODE)')
define(`CREATEARRAY',`format(`LAL%sCreateArray',TYPECODE)')
define(`DESTROYARRAY',`format(`LAL%sDestroyArray',TYPECODE)')

/* <lalVerbatim file="ArrayFactoriesD"> */
void RESIZEARRAY ( LALStatus *status, ATYPE **array, UINT4Vector *dimLength )
{  /* </lalVerbatim> */
  UINT4 arrayDataSize = 1;
  UINT4 numDims;
  UINT4 dim;
  TYPE * p; /* temporary pointer */

  INITSTATUS (status, "RESIZEARRAY", ARRAYFACTORIESC);

  ASSERT ( array != NULL, status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR );

  if ( (*array) == NULL )
    {
      CREATEARRAY ( status, array, dimLength );
    }
  else if ( dimLength == NULL )
    {
      DESTROYARRAY ( status, array );
    }
  else
    {
      ASSERT (dimLength->data,   status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
      ASSERT (dimLength->length, status,
              AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH);

      numDims = dimLength->length;

      /* loop over dimensions to compute total size of array data */
      for (dim = 0; dim < numDims; ++dim)
        arrayDataSize *= dimLength->data[dim];
        
      ASSERT (arrayDataSize > 0, status, AVFACTORIESH_ELENGTH, 
              AVFACTORIESH_MSGELENGTH);

      if (arrayDataSize == 0)
        {
          DESTROYARRAY ( status, array );
        }
      else
        {
          /* copy over the dimension data */
          (*array)->dimLength->length = dimLength->length;
          memcpy((*array)->dimLength->data, dimLength->data, 
                 numDims*sizeof(UINT4));

          /* reallocate memory for array data */
          p = LALRealloc((*array)->data, arrayDataSize*sizeof( TYPE ));

          if (p == NULL)
            {
              LALFree( (*array)->data );
              ABORT( status, AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC );
            }

          (*array)->data = p;
        }
    }
    
  RETURN (status);
}
