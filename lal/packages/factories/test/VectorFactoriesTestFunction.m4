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
define(`VTYPE',`format(`%sVector',TYPE)')
define(`CFUNC',`format(`LAL%sCreateVector',TYPECODE)')
define(`DFUNC',`format(`LAL%sDestroyVector',TYPECODE)')
define(`FUNC',`format(`%sVectorFactoriesTest',TYPECODE)')


static void FUNC ( void )
{
  const  UINT4   length = 16;
  static LALStatus  status;
  static VTYPE  *vector;
  static VTYPE   vstore;


  /*
   *
   * Test ordinary behavior.
   *
   */


  CFUNC ( &status, &vector, length );
  TestStatus( &status, CODES( 0 ), 1 );

  memset( vector->data, 0, length*sizeof( TYPE ) );

  DFUNC ( &status, &vector );
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
    CFUNC ( &status, &vector, 0 );
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

#endif

  LALCheckMemoryLeaks();
  printf( "PASS: tests of CFUNC and DFUNC \n" );
          
  return;
}
