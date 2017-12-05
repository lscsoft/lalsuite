#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRINGx(a) #a
#define STRING(a) STRINGx(a)

#define VTYPE CONCAT2(TYPE,ArraySequence)

#ifdef TYPECODE
#define CFUNC CONCAT3(LAL,TYPECODE,CreateArraySequence)
#define DFUNC CONCAT3(LAL,TYPECODE,DestroyArraySequence)
#define FUNC CONCAT2(TYPECODE,ArraySequenceFactoriesTest)
#else
#define CFUNC LALCreateArraySequence
#define DFUNC LALDestroyArraySequence
#define FUNC ArraySequenceFactoriesTest
#endif

static void FUNC ( void )
{
  CreateArraySequenceIn input;
  UINT4Vector dimLength;
  UINT4 data[] = { 2, 3, 2 };
  UINT4 dataBad[] = { 2, 3, 0 };
  static LALStatus  status;
  static VTYPE  *sequence;
  static VTYPE   sstore;


  /*
   *
   * Test ordinary behavior.
   *
   */

  dimLength.length = 3;
  dimLength.data = data;
  input.length = 2;
  input.dimLength = &dimLength;

  CFUNC ( &status, &sequence, &input );
  TestStatus( &status, CODES( 0 ), 1 );

  memset( sequence->data, 0, sequence->length*sequence->arrayDim*sizeof( TYPE ) );

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
    input.length = 0;
    CFUNC ( &status, &sequence, &input );
    TestStatus( &status, CODES( SEQFACTORIESH_ESLENGTH ), 1 );

    input.length = 2;
    input.dimLength->data = dataBad;
    CFUNC ( &status, &sequence, &input );
    TestStatus( &status, CODES( SEQFACTORIESH_EALENGTH ), 1 );

    CFUNC ( &status, &sequence, NULL );
    TestStatus( &status, CODES( SEQFACTORIESH_EINPTR ), 1 );

    DFUNC ( &status, NULL );
    TestStatus( &status, CODES( SEQFACTORIESH_EVPTR ), 1 );

    input.dimLength->data = data;
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
  (void)sstore;
  (void)dataBad;
#endif


  LALCheckMemoryLeaks();
  printf( "PASS: tests of %s and %s\n", STRING(CFUNC), STRING(DFUNC));

  return;
}

#undef CFUNC
#undef RFUNC
#undef DFUNC
#undef FUNC
