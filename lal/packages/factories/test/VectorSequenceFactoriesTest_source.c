#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRINGx(a) #a
#define STRING(a) STRINGx(a)

#define VTYPE CONCAT2(TYPE,VectorSequence)

#ifdef TYPECODE
#define CFUNC CONCAT3(LAL,TYPECODE,CreateVectorSequence)
#define DFUNC CONCAT3(LAL,TYPECODE,DestroyVectorSequence)
#define PRINTVEC CONCAT3(LAL,TYPECODE,PrintVectorSequence)
#define FUNC CONCAT2(TYPECODE,VectorSequenceFactoriesTest)
#else
#define CFUNC LALCreateVectorSequence
#define DFUNC LALDestroyVectorSequence
#define PRINTVEC LALPrintVectorSequence
#define FUNC VectorSequenceFactoriesTest
#endif

static void FUNC ( void )
{
  CreateVectorSequenceIn input   = { 2, 8 };
#ifndef LAL_NDEBUG
  CreateVectorSequenceIn badslen = { 0, 8 };
  CreateVectorSequenceIn badvlen = { 2, 0 };
#endif
  static LALStatus  status;
  static VTYPE  *sequence;
  static VTYPE   sstore;


  /*
   *
   * Test ordinary behavior.
   *
   */


  CFUNC ( &status, &sequence, &input );
  TestStatus( &status, CODES( 0 ), 1 );

  memset( sequence->data, 0, input.length*input.vectorLength*sizeof( TYPE ) );

  DFUNC ( &status, &sequence );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCheckMemoryLeaks();


  /*
   *
   * Test error codes.
   *
   */

#ifndef LAL_NDEBUG

  if ( ! lalNoDebug )
  {
    CFUNC ( &status, &sequence, &badslen );
    TestStatus( &status, CODES( SEQFACTORIESH_ESLENGTH ), 1 );

    CFUNC ( &status, &sequence, &badvlen );
    TestStatus( &status, CODES( SEQFACTORIESH_EVLENGTH ), 1 );

    CFUNC ( &status, &sequence, NULL );
    TestStatus( &status, CODES( SEQFACTORIESH_EINPTR ), 1 );

    DFUNC ( &status, NULL );
    TestStatus( &status, CODES( SEQFACTORIESH_EVPTR ), 1 );

    CFUNC ( &status, NULL, &input );
    TestStatus( &status, CODES( SEQFACTORIESH_EVPTR ), 1 );

    DFUNC ( &status, &sequence );
    TestStatus( &status, CODES( SEQFACTORIESH_EUPTR ), 1 );

    sequence = &sstore;
    CFUNC ( &status, &sequence, &input );
    TestStatus( &status, CODES( SEQFACTORIESH_EUPTR ), 1 );

    DFUNC ( &status, &sequence );
    TestStatus( &status, CODES( SEQFACTORIESH_EDPTR ), 1 );
  }

#else
  sequence = &sstore;
#endif

  LALCheckMemoryLeaks();
  printf( "PASS: tests of %s and %s\n", STRING(CFUNC), STRING(DFUNC));

  return;
}

#undef CFUNC
#undef DFUNC
#undef PRINTVEC
#undef FUNC
