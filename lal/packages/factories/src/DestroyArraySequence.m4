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
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`STYPE',`format(`%sArraySequence',TYPE)')
define(`FUNC',`format(`LAL%sDestroyArraySequence',TYPECODE)')

/* <lalVerbatim file="ArraySequenceFactoriesD"> */
void FUNC ( LALStatus *status, STYPE **aseq )
{ /* </lalVerbatim> */
  /*
   * Initialize status
   */

  INITSTATUS( status, "FUNC", ARRAYSEQUENCEFACTORIESC );
  ATTATCHSTATUSPTR( status );

  /*
   * Check aseq: is it non-NULL?
   */

  ASSERT (aseq != NULL, status, SEQFACTORIESH_EVPTR, SEQFACTORIESH_MSGEVPTR);

  /*
   * Check aseq: does it point to non-NULL?
   */

  ASSERT (*aseq != NULL,status, SEQFACTORIESH_EUPTR, SEQFACTORIESH_MSGEUPTR);

  /*
   * Check dimLength in aseq: does it point to non-NULL?
   */

  ASSERT ((*aseq)->dimLength != NULL, status,
          SEQFACTORIESH_EDPTR, SEQFACTORIESH_MSGEDPTR);

  /*
   * Check data in aseq: does it point to non-NULL?
   */

  ASSERT ((*aseq)->data != NULL, status,
          SEQFACTORIESH_EDPTR, SEQFACTORIESH_MSGEDPTR);

  /* Ok, now let's free allocated storage */

  TRY( LALU4DestroyVector( status->statusPtr, &((*aseq)->dimLength) ),
       status );
  LALFree ( (*aseq)->data ); /* free allocated data */
  LALFree ( *aseq );	      /* free aseq struct itself */

  *aseq = NULL;		/* make sure we don't point to freed struct */

  DETATCHSTATUSPTR( status );
  RETURN (status);
}
