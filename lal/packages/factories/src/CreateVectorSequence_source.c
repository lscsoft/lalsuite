#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,VectorSequence)

#ifdef TYPECODE
#define FUNC CONCAT3(LAL,TYPECODE,CreateVectorSequence)
#define XFUNC CONCAT2(XLALCreate,STYPE)
#else
#define FUNC LALCreateVectorSequence
#define XFUNC XLALCreateVectorSequence
#endif

STYPE * XFUNC ( UINT4 length, UINT4 veclen )
{
  STYPE *seq;

  if ( ! length || ! veclen )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  seq = LALMalloc( sizeof( *seq ) );
  if ( ! seq )
    XLAL_ERROR_NULL( XLAL_ENOMEM );

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
      XLAL_ERROR_NULL( XLAL_ENOMEM );
    }
  }

  return seq;
}


void FUNC ( LALStatus *status, STYPE **vseq, CreateVectorSequenceIn *in )
{
  /*
   * Initialize status
   */
  INITSTATUS(status);

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

#undef STYPE
#undef FUNC
#undef XFUNC
