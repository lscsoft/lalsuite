/* -*- C -*- */

#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRINGx(a) #a
#define STRING(a) STRINGx(a)

#define VTYPE CONCAT2(TYPE,Vector)

#ifdef TYPECODE
#define CFUNC CONCAT3(LAL,TYPECODE,CreateVector)
#define RFUNC CONCAT3(LAL,TYPECODE,ResizeVector)
#define DFUNC CONCAT3(LAL,TYPECODE,DestroyVector)
#define PRINTVEC CONCAT3(LAL,TYPECODE,PrintVector)
#define FUNC CONCAT2(TYPECODE,VectorFactoriesTest)
#else
#define CFUNC LALCreateVector
#define RFUNC LALResizeVector
#define DFUNC LALDestroyVector
#define PRINTVEC LALPrintVector
#define FUNC VectorFactoriesTest
#endif

static void FUNC ( void )
{
  const  UINT4   length = 16;
  static LALStatus  status;
  static VTYPE  *vector = ( VTYPE * )NULL;
  static VTYPE   vstore;


  /*
   *
   * Test ordinary behavior.
   *
   */


  CFUNC ( &status, &vector, length );
  TestStatus( &status, CODES( 0 ), 1 );

  if (verbose)
    {
      printf("VFT line %d:\n", __LINE__);
      printf("  vector->length = %d\n", vector->length);
    }

  memset( vector->data, 0, length*sizeof( TYPE ) );

  if (verbose)
    {
      PRINTVEC ( vector );
    }


  /* Resize up */
  RFUNC ( &status, &vector, length*3 );
  TestStatus( &status, CODES( 0 ), 1 );

  if (verbose)
    {
      printf("VFT line %d:\n", __LINE__);
      printf("  vector->length = %d\n", vector->length);
    }


  memset( vector->data, 0, length*sizeof( TYPE ) );

  if (verbose)
    {
      PRINTVEC ( vector );
    }

  /* Resize down */
  RFUNC ( &status, &vector, length*2 );
  TestStatus( &status, CODES( 0 ), 1 );

  if (verbose)
    {
      printf("VFT line %d:\n", __LINE__);
      printf("  vector->length = %d\n", vector->length);
      PRINTVEC ( vector );
    }

  /* Destroy */
  DFUNC ( &status, &vector );
  TestStatus( &status, CODES( 0 ), 1 );

  if (verbose)
    {
      printf("VFT line %d:\n", __LINE__);
      printf("  vector = %p\n", (void *)vector);
    }

  LALCheckMemoryLeaks();

  /* Check the resize function with create/destroy */
  RFUNC ( &status, &vector, length );
  TestStatus( &status, CODES( 0 ), 1 );

  if (verbose)
    {
      printf("VFT line %d:\n", __LINE__);
      printf("  vector->length = %d\n", vector->length);
    }

  memset( vector->data, 0, length*sizeof( TYPE ) );

  if (verbose)
    {
      PRINTVEC ( vector );
    }

  RFUNC ( &status, &vector, 0 );
  TestStatus( &status, CODES( 0 ), 1);

  if (verbose)
    {
      printf("VFT line %d:\n", __LINE__);
      printf("  vector = %p\n", (void *)vector);
    }

  LALCheckMemoryLeaks();


  /*
   *
   * Test error codes.
   *
   */

#ifndef LAL_NDEBUG

  if ( ! lalNoDebug )
  {
    if (verbose)
      {
        printf("VFT line %d\n", __LINE__);
        printf("  vector = %p\n", (void *)vector);
      }

    CFUNC ( &status, &vector, 0 );
    TestStatus( &status, CODES( AVFACTORIESH_ELENGTH ), 1 );

    RFUNC ( &status, &vector, 0 );
    TestStatus( &status, CODES( AVFACTORIESH_ELENGTH ), 1 );

    DFUNC ( &status, NULL );
    TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1 );

    CFUNC ( &status, NULL, length );
    TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1 );

    DFUNC ( &status, &vector );
    TestStatus( &status, CODES( AVFACTORIESH_EUPTR ), 1 );

    vector = &vstore;
    CFUNC ( &status, &vector, length );
    TestStatus( &status, CODES( AVFACTORIESH_EUPTR ), 1 );

    DFUNC ( &status, &vector );
    TestStatus( &status, CODES( AVFACTORIESH_EDPTR ), 1 );
  }

#else
  vector = &vstore;
#endif

  LALCheckMemoryLeaks();
  printf( "PASS: tests of %s, %s, and %s\n", STRING(CFUNC), STRING(RFUNC), STRING(DFUNC));

  return;
}

#undef CFUNC
#undef RFUNC
#undef DFUNC
#undef PRINTVEC
#undef FUNC
