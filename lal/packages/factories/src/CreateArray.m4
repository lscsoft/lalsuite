dnl m4 template for defining CreateArray for different types. 
dnl To create the c-source for the typecode x variant, execute
dnl     m4 -DTYPECODE=x template.m4 > xtemplate.c
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
define(`PROG',`format(`%sCreateArray',TYPECODE)')
define(`PROGFILE',`format(`%s.c',PROG)')
/*----------------------------------------------------------------------- 
 * 
 * File Name: PROGFILE
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * PROG
 * 
 * SYNOPSIS 
 * void PROG (Status *, ATYPE **array, UINT4Vector *dimLength);
 * 
 * DESCRIPTION 
 * Create a ATYPE object. 
 * 
 * DIAGNOSTICS 
 * Illegal length, array == NULL, *array != NULL, malloc failure
 *
 * CALLS
 * LALMalloc
 * CreateVector
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _STRING_H
#include <string.h>
#ifndef _STRING_H
#define _STRING_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _AVFACTORIES_H
#include "AVFactories.h"
#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H
#endif
#endif

NRCSID (CREATEARRAYC, "$Id$");

typedef TYPE  dtype;    /* change for different factory */
typedef ATYPE atype;    /* change for different factory */

void PROG (Status       *status, 
           atype       **array,
           UINT4Vector  *dimLength) 
{
  UINT4 arrayDataSize = 0;
  UINT4 numDims;
  UINT4 dim;

  INITSTATUS (status, CREATEARRAYC);   
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are sane */

  ASSERT (array,             status, CREATEARRAY_EVPTR, CREATEARRAY_MSGEVPTR);
  ASSERT (!*array,           status, CREATEARRAY_EUPTR, CREATEARRAY_MSGEUPTR);
  ASSERT (dimLength,         status, CREATEARRAY_EVPTR, CREATEARRAY_MSGEVPTR);
  ASSERT (dimLength->data,   status, CREATEARRAY_EVPTR, CREATEARRAY_MSGEVPTR);
  ASSERT (dimLength->length, status,
          CREATEARRAY_ELENGTH, CREATEARRAY_MSGELENGTH);

  numDims = dimLength->length;

  /* allocate memory for structure */

  *array = (atype *) LALMalloc (sizeof(atype));
  ASSERT (*array, status, CREATEARRAY_EMALLOC, CREATEARRAY_MSGEMALLOC);

  (*array)->dimLength = NULL;
  (*array)->data      = NULL;

  /* allocate dimLength field and copy information there */

  U4CreateVector (status->statusPtr, &(*array)->dimLength, numDims);
  if (status->statusPtr->statusCode)
  {
    LALFree (*array);
    RETURN (status);   /* this returns a recursive error status code (-1) */
  }
  memcpy ((*array)->dimLength->data, dimLength->data, numDims*sizeof(UINT4));

  /* loop over dimensions to compute total size of array data */

  for (dim = 0; dim < numDims; ++dim)
  {
    arrayDataSize *= dimLength->data[dim];
  }

  if (!arrayDataSize)
  {
    LALFree (*array);
    ABORT (status, CREATEARRAY_ELENGTH, CREATEARRAY_MSGELENGTH);
  }
  
  /* allocate storage */

  (*array)->data = (dtype *) LALMalloc (arrayDataSize*sizeof(dtype));

  if (!(*array)->data)
  {
    /* try to free memory                                               */
    /* this should ALWAYS work so we can use the CHECKSTATUSPTR() macro */
    U4DestroyVector (status->statusPtr, &(*array)->dimLength);
    CHECKSTATUSPTR (status);  /* if this fails, there is a memory leak  */

    LALFree (*array);
    ABORT (status, CREATEARRAY_EMALLOC, CREATEARRAY_MSGEMALLOC);
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

