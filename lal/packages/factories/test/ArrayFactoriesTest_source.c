#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRINGx(a) #a
#define STRING(a) STRINGx(a)

#define VTYPE CONCAT2(TYPE,Array)

#ifdef TYPECODE
#define CFUNC CONCAT3(LAL,TYPECODE,CreateArray)
#define RFUNC CONCAT3(LAL,TYPECODE,ResizeArray)
#define DFUNC CONCAT3(LAL,TYPECODE,DestroyArray)
#define FUNC CONCAT2(TYPECODE,ArrayFactoriesTest)
#else
#define CFUNC LALCreateArray
#define RFUNC LALResizeArray
#define DFUNC LALDestroyArray
#define FUNC ArrayFactoriesTest
#endif

static void FUNC ( void )
{
  static UINT4   dims[3]    = { 1, 2, 4 };
  UINT4Vector    dimLength  = { 3, dims };
#ifndef LAL_NDEBUG
  static UINT4   dbad[3]    = { 1, 0, 4 };
  UINT4Vector    badLength1 = { 3, NULL };
  UINT4Vector    badLength2 = { 0, dims };
  UINT4Vector    badLength3 = { 3, dbad };
#endif
  static LALStatus  status;
  static VTYPE  *array;
  static VTYPE   astore;
  static TYPE    datum;


  /*
   *
   * Test ordinary behavior.
   *
   */


  CFUNC ( &status, &array, &dimLength );
  TestStatus( &status, CODES( 0 ), 1 );

  memset( array->data, 0, dims[0]*dims[1]*dims[2]*sizeof( TYPE ) );

  /* resize up */
  /*
   * dimLength.data[0] *= 2;
   * dimLength.data[1] *= 3;
   * dimLength.data[2] *= 4;
  */
  dims[0] *= 2;
  dims[1] *= 3;
  dims[2] *= 4;
  RFUNC ( &status, &array, &dimLength );
  TestStatus( &status, CODES( 0 ), 1 );

  memset( array->data, 0, dims[0]*dims[1]*dims[2]*sizeof( TYPE ) );

  /* resize down */
  dims[0] /= 2;
  dims[1] /= 3;
  dims[2] /= 2;
  RFUNC ( &status, &array, &dimLength );
  TestStatus( &status, CODES( 0 ), 1 );

  memset( array->data, 0, dims[0]*dims[1]*dims[2]*sizeof( TYPE ) );

  /* resize down again */
  dims[2] /= 2;
  RFUNC ( &status, &array, &dimLength );
  TestStatus( &status, CODES( 0 ), 1 );

  memset( array->data, 0, dims[0]*dims[1]*dims[2]*sizeof( TYPE ) );

  DFUNC ( &status, &array );
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
    CFUNC ( &status, &array, &badLength1 );
    TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1 );

    RFUNC ( &status, &array, &badLength1 );
    TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1 );

    CFUNC ( &status, &array, &badLength2 );
    TestStatus( &status, CODES( AVFACTORIESH_ELENGTH ), 1 );

    RFUNC ( &status, &array, &badLength2 );
    TestStatus( &status, CODES( AVFACTORIESH_ELENGTH ), 1 );

    CFUNC ( &status, &array, &badLength3 );
    TestStatus( &status, CODES( AVFACTORIESH_ELENGTH ), 1 );

    RFUNC ( &status, &array, &badLength3 );
    TestStatus( &status, CODES( AVFACTORIESH_ELENGTH ), 1 );
    LALCheckMemoryLeaks();

    DFUNC ( &status, NULL );
    TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1 );

    CFUNC ( &status, NULL, &dimLength );
    TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1 );

    RFUNC ( &status, NULL, &badLength1 );
    TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1 );

    DFUNC ( &status, &array );
    TestStatus( &status, CODES( AVFACTORIESH_EUPTR ), 1 );

    array = &astore;
    CFUNC ( &status, &array, &dimLength );
    TestStatus( &status, CODES( AVFACTORIESH_EUPTR ), 1 );

    RFUNC ( &status, &array, &badLength1 );
    TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1);

    RFUNC ( &status, &array, &badLength2 );
    TestStatus( &status, CODES( AVFACTORIESH_ELENGTH ), 1);

    RFUNC ( &status, &array, &badLength3 );
    TestStatus( &status, CODES( AVFACTORIESH_ELENGTH ), 1);

    DFUNC ( &status, &array );
    TestStatus( &status, CODES( AVFACTORIESH_EDPTR ), 1 );

    array->data = &datum;
    DFUNC ( &status, &array );
    TestStatus( &status, CODES( AVFACTORIESH_EDPTR ), 1 );
    ClearStatus( &status );
  }

#else
  array = &astore;
  array->data = &datum;
#endif

  LALCheckMemoryLeaks();
  printf( "PASS: tests of %s, %s, and %s\n", STRING(CFUNC), STRING(RFUNC), STRING(DFUNC));

  return;
}

#undef CFUNC
#undef RFUNC
#undef DFUNC
#undef FUNC
