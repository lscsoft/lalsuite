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
define(`FUNC',`format(`LAL%sDestroyVectorSequence',TYPECODE)')
ifelse( TYPECODE, `', `define(`XFUNC',`XLALDestroyVectorSequence')', `define(`XFUNC',`format(`XLALDestroy%s',STYPE)')' ) 

void XFUNC ( STYPE *vseq )
{
  if ( ! vseq )
    XLAL_ERROR_VOID( "XFUNC", XLAL_EFAULT );
  if ( ! vseq->data || ! vseq->length || ! vseq->vectorLength )
    XLAL_ERROR_VOID( "XFUNC", XLAL_EINVAL );
  LALFree( vseq->data );
  memset( vseq, 0, sizeof( *vseq ) );
  LALFree( vseq );
  return;
}

/* <lalVerbatim file="VectorSequenceFactoriesD"> */
void FUNC ( LALStatus *status, STYPE **vseq )
{ /* </lalVerbatim> */
  /* 
   * Initialize status
   */

  INITSTATUS( status, "FUNC", VECTORSEQUENCEFACTORIESC );
      
  /* 
   * Check vseq: is it non-NULL?
   */

  ASSERT (vseq != NULL, status, SEQFACTORIESH_EVPTR, SEQFACTORIESH_MSGEVPTR); 

  /* 
   * Check vseq: does it point to non-NULL?
   */

  ASSERT (*vseq != NULL,status, SEQFACTORIESH_EUPTR, SEQFACTORIESH_MSGEUPTR);

  /*
   * Check data in vseq: does it point to non-NULL?
   */

  ASSERT ((*vseq)->data != NULL, status,
          SEQFACTORIESH_EDPTR, SEQFACTORIESH_MSGEDPTR);

  /* Ok, now let's free allocated storage */

  XFUNC ( *vseq );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    if ( code == XLAL_EFAULT )
    {
      ABORT (status, SEQFACTORIESH_EUPTR, SEQFACTORIESH_MSGEUPTR);
    }
    if ( code == XLAL_EINVAL )
    {
      ABORT (status, SEQFACTORIESH_EDPTR, SEQFACTORIESH_MSGEDPTR);
    }
  }

  *vseq = NULL;		/* make sure we don't point to freed struct */

  RETURN (status);
}
