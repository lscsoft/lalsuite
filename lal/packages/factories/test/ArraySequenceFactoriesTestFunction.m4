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
define(`VTYPE',`format(`%sArraySequence',TYPE)')
define(`CFUNC',`format(`LAL%sCreateArraySequence',TYPECODE)')
define(`DFUNC',`format(`LAL%sDestroyArraySequence',TYPECODE)')
define(`FUNC',`format(`%sArraySequenceFactoriesTest',TYPECODE)')


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

#endif

  LALCheckMemoryLeaks();
  printf( "PASS: tests of CFUNC and DFUNC \n" );

  return;
}
