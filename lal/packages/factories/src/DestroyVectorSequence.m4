dnl m4 template for defining DestroyVector for different types. 
dnl To create the c-source for the typecode x program, execute
dnl     m4 -DTYPECODE=x template.m4
dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`STYPE',`format(`%sVectorSequence',TYPE)')
define(`PROG',`format(`%sDestroyVectorSequence',TYPECODE)')
define(`CPROG',`format(`%sCreateVectorSequence',TYPECODE)')
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
 * void PROG ( Status *,  STYPE **vseq );
 * 
 * DESCRIPTION 
 * Returns to system storage allocated by CPROG
 * 
 * DIAGNOSTICS 
 * vseq == NULL, *vseq == NULL, (*vseq)->data == NULL, free failure
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

#ifndef _SEQFACTORIES_H
#include "SeqFactories.h"
#ifndef _SEQFACTORIES_H
#define _SEQFACTORIES_H
#endif
#endif

NRCSID (DESTROYVECTORSEQUENCEC, "$Id$");

typedef TYPE dtype;		/* change for different factory */
typedef STYPE stype;		/* change for different factory */

void PROG (Status *status, stype **vseq)
{
  /* 
   * Initialize status
   */

  INITSTATUS (status, DESTROYVECTORSEQUENCEC);
      
  /* 
   * Check vseq: is it non-NULL?
   */

  ASSERT (vseq != NULL, status, DESTROYVECSEQ_EVPTR, DESTROYVECSEQ_MSGEVPTR); 

  /* 
   * Check vseq: does it point to non-NULL?
   */

  ASSERT (*vseq != NULL,status, DESTROYVECSEQ_EUPTR, DESTROYVECSEQ_MSGEUPTR);

  /*
   * Check data in vseq: does it point to non-NULL?
   */

  ASSERT ((*vseq)->data != NULL, status,
          DESTROYVECSEQ_EDPTR, DESTROYVECSEQ_MSGEDPTR);

  /* Ok, now let's free allocated storage */

  LALFree ( (*vseq)->data ); /* free allocated data */
  LALFree ( *vseq );	      /* free vseq struct itself */

  *vseq = NULL;		/* make sure we don't point to freed struct */

  RETURN (status);
}
