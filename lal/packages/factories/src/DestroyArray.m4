dnl m4 template for defining DestroyArray different types. 
dnl To create the c-source for the typecode x program, execute
dnl     m4 -DTYPECODE=x template.m4
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
define(`PROG',`format(`%sDestroyArray',TYPECODE)')
define(`CPROG',`format(`%sCreateArray',TYPECODE)')
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
 * void PROG ( Status *,  ATYPE **array );
 * 
 * DESCRIPTION 
 * Returns to system storage allocated by CPROG
 * 
 * DIAGNOSTICS 
 * array == NULL, *array == NULL, (*array)->data == NULL, free failure
 *
 * CALLS
 * LALFree
 * DestroyVector
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

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

NRCSID (DESTROYARRAYC, "$Id$");

typedef TYPE  dtype;		/* change for different factory */
typedef ATYPE atype;		/* change for different factory */

void PROG (Status *status, atype **array)
{
  INITSTATUS (status, DESTROYARRAYC);	
  ATTATCHSTATUSPTR (status);
      
  ASSERT (array,          status, DESTROYARRAY_EVPTR, DESTROYARRAY_MSGEVPTR);
  ASSERT (*array,         status, DESTROYARRAY_EUPTR, DESTROYARRAY_MSGEUPTR);
  ASSERT ((*array)->data, status, DESTROYARRAY_EDPTR, DESTROYARRAY_MSGEDPTR);

  /* Free allocated storage */

  U4DestroyVector (status->statusPtr, &(*array)->dimLength);
  CHECKSTATUSPTR (status);

  LALFree ((*array)->data); /* free allocated data */
  LALFree (*array);	    /* free array struct itself */
  *array = NULL;	    /* make sure we don't point to freed struct */

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
