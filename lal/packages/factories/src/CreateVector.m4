dnl m4 template for defining CreateVector for different types. 
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
define(`VTYPE',`format(`%sVector',TYPE)')
define(`PROG',`format(`%sCreateVector',TYPECODE)')
define(`PROGFILE',`format(`%s.c',PROG)')
/*----------------------------------------------------------------------- 
 * 
 * File Name: PROGFILE
 * 
 * Author: Finn, L. S.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * PROG
 * 
 * SYNOPSIS 
 * void PROG (Status *, VTYPE **vector, UINT4 length);
 * 
 * DESCRIPTION 
 * Create a VTYPE object. 
 * 
 * DIAGNOSTICS 
 * Illegal length, vector == NULL, *vector != NULL, malloc failure
 *
 * CALLS
 * LALMalloc
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

NRCSID (CREATEVECTORC, "$Id$");

typedef TYPE dtype;		/* change for different factory */
typedef VTYPE vtype;	/* change for different factory */

void PROG (Status* status, 
	    vtype** vector,
	    UINT4 length) 
{
  /* 
   * Initialize status structure
   */

  INITSTATUS (status, CREATEVECTORC);	
      
  /* Check sequence length: report error if 0 
   * Use of unsigned for length means we can't check if negative
   * length was passed
   */

  ASSERT (length > 0, status, CREATEVECTOR_ELENGTH, CREATEVECTOR_MSGELENGTH);

  /* 
   * Check return structure: If return pointer does not point to a
   *    valid pointer then report an error 
   */

  ASSERT (vector != NULL, status, CREATEVECTOR_EVPTR, CREATEVECTOR_MSGEVPTR);

  ASSERT (*vector == NULL, status, CREATEVECTOR_EUPTR, CREATEVECTOR_MSGEUPTR);

  /*
   * Allocate pointer
   */

  *vector = (vtype *) LALMalloc(sizeof(vtype));
  ASSERT (*vector != NULL, status,
          CREATEVECTOR_EMALLOC, CREATEVECTOR_MSGEMALLOC);

  (*vector)->length = 0;	/* length 0 until storage allocated */
  (*vector)->data   = NULL;	/* NULL data until allocated */

  /* 
   * Allocate storage 
   * Test that storage is properly allocated. Can't handle with ASSERT
   * since we need to de-allocate structure pointer before an error return
   */

  (*vector)->data = (dtype *) LALMalloc(sizeof(dtype)*length);

  if ( NULL == (*vector)->data )
  {
    /* Must free storage pointed to by *vector */
    LALFree((void *) *vector);
    ABORT (status, CREATEVECTOR_EMALLOC, CREATEVECTOR_MSGEMALLOC);
    return;
  }
  (*vector)->length = length;	/* Set length if storage allocated */

  /* We be done: Normal exit */

  RETURN (status);
}
