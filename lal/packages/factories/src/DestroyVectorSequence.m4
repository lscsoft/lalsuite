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
define(`FUNC',`format(`%sDestroyVectorSequence',TYPECODE)')

void FUNC ( Status *status, STYPE **vseq )
{
  /* 
   * Initialize status
   */

  INITSTATUS( status, "FUNC", VECTORSEQUENCEFACTORIESC );
      
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
