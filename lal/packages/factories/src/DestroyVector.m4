dnl m4 template for defining DestroyVector for different types. 
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
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`VTYPE',`format(`%sVector',TYPE)')
define(`PROG',`format(`%sDestroyVector',TYPECODE)')
define(`CPROG',`format(`%sCreateVector',TYPECODE)')
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
 * INT4 PROG ( Status *,  VTYPE **vector );
 * 
 * DESCRIPTION 
 * Returns to system storage allocated by CPROG
 * 
 * DIAGNOSTICS 
 * vector == NULL, *vector == NULL, (*vector)->data == NULL, free failure
 *
 * CALLS
 * LALFree
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

NRCSID (DESTROYVECTORC, "$Id$");

typedef TYPE dtype;		/* change for different factory */
typedef VTYPE vtype;		/* change for different factory */

void PROG (Status *status, vtype **vector)
{
  /* 
   * Initialize status
   */

  INITSTATUS (status, DESTROYVECTORC);	
      
  /* 
   * Check vector: is it non-NULL?
   */

  ASSERT (vector != NULL, status, DESTROYVECTOR_EVPTR, DESTROYVECTOR_MSGEVPTR);

  /* 
   * Check vector: does it point to non-NULL?
   */

  ASSERT (*vector != NULL, status, DESTROYVECTOR_EUPTR, DESTROYVECTOR_MSGEUPTR);

  /*
   * Check data in vector: does it point to non-NULL
   */

  ASSERT ((*vector)->data != NULL, status,
          DESTROYVECTOR_EDPTR, DESTROYVECTOR_MSGEDPTR);

  /* Ok, now let's free allocated storage */

  LALFree ( (*vector)->data ); /* free allocated data */
  LALFree ( *vector );	/* free vector struct itself */

  *vector = NULL;		/* make sure we don't point to freed struct */

  RETURN (status);
}
