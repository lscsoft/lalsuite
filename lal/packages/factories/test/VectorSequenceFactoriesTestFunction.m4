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
define(`VTYPE',`format(`%sVectorSequence',TYPE)')
define(`CFUNC',`format(`LAL%sCreateVectorSequence',TYPECODE)')
define(`DFUNC',`format(`LAL%sDestroyVectorSequence',TYPECODE)')
define(`FUNC',`format(`%sVectorSequenceFactoriesTest',TYPECODE)')


static void FUNC ( void )
{
  CreateVectorSequenceIn input   = { 2, 8 };
  CreateVectorSequenceIn badslen = { 0, 8 };
  CreateVectorSequenceIn badvlen = { 2, 0 };
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
  badslen.length = 0;
  badvlen.vectorLength = 0;
  sequence = &sstore;
#endif

  LALCheckMemoryLeaks();
  printf( "PASS: tests of CFUNC and DFUNC \n" );

  return;
}
