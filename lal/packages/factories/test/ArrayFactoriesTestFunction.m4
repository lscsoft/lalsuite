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
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`VTYPE',`format(`%sArray',TYPE)')
define(`CFUNC',`format(`LAL%sCreateArray',TYPECODE)')
define(`DFUNC',`format(`LAL%sDestroyArray',TYPECODE)')
define(`FUNC',`format(`%sArrayFactoriesTest',TYPECODE)')


static void FUNC ( void )
{
  static UINT4   dims[3]    = { 1, 2, 4 };
  static UINT4   dbad[3]    = { 1, 0, 4 };
  UINT4Vector    dimLength  = { 3, dims };
  UINT4Vector    badLength1 = { 3, NULL };
  UINT4Vector    badLength2 = { 0, dims };
  UINT4Vector    badLength3 = { 3, dbad };
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

  DFUNC ( &status, &array );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCheckMemoryLeaks();


  /*
   *
   * Test error codes.
   *
   */


  CFUNC ( &status, &array, &badLength1 );
  TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1 );

  CFUNC ( &status, &array, &badLength2 );
  TestStatus( &status, CODES( AVFACTORIESH_ELENGTH ), 1 );

  CFUNC ( &status, &array, &badLength3 );
  TestStatus( &status, CODES( AVFACTORIESH_ELENGTH ), 1 );
  LALCheckMemoryLeaks();

  DFUNC ( &status, NULL );
  TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1 );

  CFUNC ( &status, NULL, &dimLength );
  TestStatus( &status, CODES( AVFACTORIESH_EVPTR ), 1 );

  DFUNC ( &status, &array );
  TestStatus( &status, CODES( AVFACTORIESH_EUPTR ), 1 );

  array = &astore;
  CFUNC ( &status, &array, &dimLength );
  TestStatus( &status, CODES( AVFACTORIESH_EUPTR ), 1 );

  DFUNC ( &status, &array );
  TestStatus( &status, CODES( AVFACTORIESH_EDPTR ), 1 );

  array->data = &datum;
  DFUNC ( &status, &array );
  TestStatus( &status, CODES( -1 ), 1 );
  ClearStatus( &status );

  LALCheckMemoryLeaks();
  printf( "PASS... tests of CFUNC and DFUNC \n" );
          
  return;
}
