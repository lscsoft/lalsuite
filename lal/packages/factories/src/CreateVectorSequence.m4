dnl m4 template for defining CreateArraySequence for different types. 
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
define(`STYPE',`format(`%sVectorSequence',TYPE)')
define(`PROG',`format(`%sCreateVectorSequence',TYPECODE)')
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
 * void PROG (Status *, STYPE **vseq, CreateVectorSequenceIn *in);
 * 
 * DESCRIPTION 
 * Create a STYPE object. 
 * 
 * DIAGNOSTICS 
 * Illegal sequence length, illegal vectorLength, vseq == NULL, 
 * *vseq != NULL, malloc failure 
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

#ifndef _SEQFACTORIES_H
#include "SeqFactories.h"
#ifndef _SEQFACTORIES_H
#define _SEQFACTORIES_H
#endif
#endif

NRCSID (CREATEVECTORSEQUENCEC, "$Id$");

typedef TYPE dtype;		/* change for different factory */
typedef STYPE stype;	/* change for different factory */

void PROG (Status *status, 
	    stype **vseq,
	    CreateVectorSequenceIn *in) 
{
  /* 
   * Initialize status
   */
  INITSTATUS (status, CREATEVECTORSEQUENCEC);	

  /* Check input structure: report if NULL */

  ASSERT (in != NULL, status, CREATEVECSEQ_EINPTR, CREATEVECSEQ_MSGEINPTR);
      
  /* Check sequence length: report error if 0 
   * Use of unsigned for length means we can't check if negative
   * length was passed
   */

  ASSERT (in->length > 0, status,
          CREATEVECSEQ_ESLENGTH, CREATEVECSEQ_MSGESLENGTH);

  /* Check vector length: report error if 0 
   * Use of unsigned for length means we can't check if negative
   * length was passed
   */

  ASSERT (in->vectorLength > 0, status,
          CREATEVECSEQ_EVLENGTH, CREATEVECSEQ_MSGEVLENGTH); 

  /* 
   * Check return structure: If return pointer does not point to a
   *    valid pointer then report an error 
   */

  ASSERT (vseq != NULL, status, CREATEVECSEQ_EVPTR, CREATEVECSEQ_MSGEVPTR);
  ASSERT (*vseq == NULL, status, CREATEVECSEQ_EUPTR, CREATEVECSEQ_MSGEUPTR);

  /*
   * Allocate pointer
   */

  *vseq = (stype *) LALMalloc(sizeof(stype));
  ASSERT (*vseq != NULL, status, CREATEVECSEQ_EMALLOC, CREATEVECSEQ_MSGEMALLOC);

  (*vseq)->length = 0;	/* length 0 until storage allocated */
  (*vseq)->vectorLength = 0; /* vector length 0 until storage allocated */
  (*vseq)->data   = NULL;	/* NULL data until allocated */

  /* 
   * Allocate storage 
   */

  {
    size_t tlength;
    tlength = in->vectorLength * in->length * sizeof(dtype);
    (*vseq)->data = (dtype *) LALMalloc (tlength);
  }

  if (NULL == (*vseq)->data)
  {
    /* Must free storage pointed to by *vseq */
    LALFree ((void *) *vseq);
    ABORT (status, CREATEVECSEQ_EMALLOC, CREATEVECSEQ_MSGEMALLOC);
    return;
  }
 
  /* Set length, vectorLength if storage allocated */

  (*vseq)->length = in->length;	
  (*vseq)->vectorLength = in->vectorLength;

  /* We be done: Normal exit */

  RETURN (status);
}
