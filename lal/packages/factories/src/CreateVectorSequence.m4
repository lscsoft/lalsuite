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
define(`FUNC',`format(`LAL%sCreateVectorSequence',TYPECODE)')
ifelse( TYPECODE, `', `define(`XFUNC',`XLALCreateVectorSequence')', `define(`XFUNC',`format(`XLALCreate%s',STYPE)')' ) 


STYPE * XFUNC ( UINT4 length, UINT4 veclen )
{
  STYPE *seq;

  if ( ! length || ! veclen )

  seq = LALMalloc( sizeof( *seq ) );
  if ( ! seq )
    XLAL_ERROR_NULL( "XFUNC", XLAL_ENOMEM );
  
  seq->length = length;
  seq->vectorLength = veclen;

  if ( ! length || ! veclen )
    seq->data = NULL;
  else
  {
    seq->data = LALMalloc( length * veclen * sizeof( *seq->data ) );
    if ( ! seq )
    {
      LALFree( seq );
      XLAL_ERROR_NULL( "XFUNC", XLAL_ENOMEM );
    }
  }

  return seq;
}

/* <lalVerbatim file="VectorSequenceFactoriesD"> */
void FUNC ( LALStatus *status, STYPE **vseq, CreateVectorSequenceIn *in ) 
{ /* </lalVerbatim> */
  /* 
   * Initialize status
   */
  INITSTATUS( status, "FUNC", VECTORSEQUENCEFACTORIESC );	

  /* Check input structure: report if NULL */

  ASSERT (in != NULL, status, SEQFACTORIESH_EINPTR, SEQFACTORIESH_MSGEINPTR);
      
  /* Check sequence length: report error if 0 
   * Use of unsigned for length means we can't check if negative
   * length was passed
   */

  ASSERT (in->length > 0, status,
          SEQFACTORIESH_ESLENGTH, SEQFACTORIESH_MSGESLENGTH);

  /* Check vector length: report error if 0 
   * Use of unsigned for length means we can't check if negative
   * length was passed
   */

  ASSERT (in->vectorLength > 0, status,
          SEQFACTORIESH_EVLENGTH, SEQFACTORIESH_MSGEVLENGTH); 

  /* 
   * Check return structure: If return pointer does not point to a
   *    valid pointer then report an error 
   */

  ASSERT (vseq != NULL, status, SEQFACTORIESH_EVPTR, SEQFACTORIESH_MSGEVPTR);
  ASSERT (*vseq == NULL, status, SEQFACTORIESH_EUPTR, SEQFACTORIESH_MSGEUPTR);


  *vseq = XFUNC ( in->length, in->vectorLength );
  if ( ! vseq )
  {
    int code = xlalErrno;
    XLALClearErrno();
    if ( code == XLAL_EBADLEN )
    {
      if ( ! in->length )
      {
        ABORT (status, SEQFACTORIESH_ESLENGTH, SEQFACTORIESH_MSGESLENGTH);
      }
      else
      {
        ABORT (status, SEQFACTORIESH_EVLENGTH, SEQFACTORIESH_MSGEVLENGTH); 
      }
    }
    if ( code == XLAL_ENOMEM )
    {
      ABORT( status, SEQFACTORIESH_EMALLOC, SEQFACTORIESH_MSGEMALLOC );
    }
  }

  /* We be done: Normal exit */

  RETURN (status);
}
