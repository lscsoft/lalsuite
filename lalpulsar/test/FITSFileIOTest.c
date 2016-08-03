//
// Copyright (C) 2016 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include <lal/FITSFileIO.h>
#include <lal/LALPulsarVCSInfo.h>
#include <lal/StringVector.h>
#include <lal/Date.h>
#include <lal/UserInput.h>
#include <lal/GSLHelpers.h>

#if !defined(HAVE_LIBCFITSIO)

int main( void )
{
  fprintf( stderr, "CFITSIO library is not available; skipping test\n" );
  return 77;
}

#else // defined(HAVE_LIBCFITSIO)

const CHAR *longstring_ref = \
  "This is a long string #1. This is a long string #2. " \
  "This is a long string #3. This is a long string #4. " \
  "This is a long string #5. This is a long string #6. " \
  "This is a long string #7. This is a long string #8. " \
  "This is a long string #9. This is a long string #10." ;

typedef struct {
  INT4 n;
  REAL4 v;
  CHAR desc[4];
} TestSubRecord;

typedef struct {
  INT4 index;
  BOOLEAN flag;
  CHAR name[8];
  LIGOTimeGPS epoch;
  struct {
    REAL4 sky;
    REAL8 freq;
    REAL8 fkdot[5];
  } pos;
  REAL8 values[2];
  COMPLEX8 phasef;
  COMPLEX16 phase;
  const TestSubRecord *sub;
} TestRecord;

const TestSubRecord testsub[3][2] = {
  { { .n=0, .v=1.23, .desc="A1" }, { .n=1, .v=2.34, .desc="X9" } },
  { { .n=0, .v=3.45, .desc="B3" }, { .n=1, .v=4.56, .desc="Y7" } },
  { { .n=0, .v=5.67, .desc="C5" }, { .n=1, .v=6.78, .desc="Z5" } },
};

const TestRecord testtable[3] = {
  { .index=3, .flag=1, .name="CasA", .epoch={123456789, 5}, .pos={.sky=5.4321, .freq=100.999, .fkdot={1e-9, 5e-20}}, .values={13.24, 43.234}, .phasef=crectf( 1.2, 3.4 ), .phase=crect( 4.5, 0.2 ), .sub=testsub[0] },
  { .index=2, .flag=0, .name="Vela", .epoch={452456245, 9}, .pos={.sky=34.454, .freq=1345.34, .fkdot={2e-8, 6e-21}}, .values={14.35, 94.128}, .phasef=crectf( 3.6, 9.3 ), .phase=crect( 8.3, 4.0 ), .sub=testsub[1] },
  { .index=1, .flag=1, .name="Crab", .epoch={467846774, 4}, .pos={.sky=64.244, .freq=15.6463, .fkdot={4e-6,     0}}, .values={153.4, 3.0900}, .phasef=crectf( 6.7, 4.4 ), .phase=crect( 5.6, 6.3 ), .sub=testsub[2] },
};

int main( int argc, char *argv[] )
{

  // Create a dummy user enviroment, for testing XLALFITSFileWriteUVarCmdLine()
  struct uvar_type { INT4 dummy; } uvar_struct = { .dummy = 0 };
  struct uvar_type *const uvar = &uvar_struct;
  XLAL_CHECK_MAIN( XLALRegisterUvarMember( dummy, INT4, 0, OPTIONAL, "Dummy option" ) == XLAL_SUCCESS, XLAL_EFUNC );
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( !should_exit, XLAL_EFAILED );

  // Write an example FITS file
  {
    FITSFile *file = XLALFITSFileOpenWrite( "FITSFileIOTest.fits" );
    XLAL_CHECK_MAIN( file != NULL, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALFITSFileWriteVCSInfo( file, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALFITSFileWriteUVarCmdLine( file ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALFITSFileWriteHistory( file, "%s\n%s", "This is a test history", longstring_ref ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALFITSHeaderWriteComment( file, "%s\n%s", "This is a test comment", longstring_ref ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: opened 'FITSFileIOTest.fits' for writing\n" );

    XLAL_CHECK_MAIN( XLALFITSHeaderWriteBOOLEAN( file, "testbool", 1, "This is a test BOOLEAN" ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: wrote a BOOLEAN\n" );

    XLAL_CHECK_MAIN( XLALFITSHeaderWriteINT4( file, "testint [s]", 2345, "This is a test INT4" ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: wrote a INT4\n" );

    XLAL_CHECK_MAIN( XLALFITSHeaderWriteINT8( file, "testint2", LAL_INT4_MAX + 6789, "This is a test INT8" ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: wrote a INT8\n" );

    XLAL_CHECK_MAIN( XLALFITSHeaderWriteREAL4( file, "testflt", LAL_PI, "This is a test REAL4" ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: wrote a REAL4\n" );

    XLAL_CHECK_MAIN( XLALFITSHeaderWriteREAL8( file, "testdbl [Hz]", LAL_E, "This is a test REAL8" ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: wrote a REAL8\n" );

    XLAL_CHECK_MAIN( XLALFITSHeaderWriteCOMPLEX8( file, "testcmp", crectf( LAL_PI_2, LAL_PI_4 ), "This is a test COMPLEX8" ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: wrote a COMPLEX8\n" );

    XLAL_CHECK_MAIN( XLALFITSHeaderWriteCOMPLEX16( file, "testdblcmp", crect( LAL_LOG2E, LAL_LOG10E ), "This is a test COMPLEX16" ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: wrote a COMPLEX16\n" );

    XLAL_CHECK_MAIN( XLALFITSHeaderWriteString( file, "teststr", "This is a short string", "This is a test string" ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: wrote a string\n" );

    XLAL_CHECK_MAIN( XLALFITSHeaderWriteString( file, "longstring", longstring_ref, "This is a long test string" ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: wrote a long string\n" );

    LALStringVector *testsv = XLALCreateStringVector( "abc", "def", "ghij", NULL );
    XLAL_CHECK_MAIN( testsv != NULL, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALFITSHeaderWriteStringVector( file, "testsv", testsv, "These are test string vector entries" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLALDestroyStringVector( testsv );
    fprintf( stderr, "PASSED: wrote a string vector\n" );

    LIGOTimeGPS ref_time = { 987654321, 123456789 };
    XLAL_CHECK_MAIN( XLALFITSHeaderWriteGPSTime( file, "date-obs", &ref_time, "This is the reference time" ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: wrote a GPS time\n" );

    {
      const size_t dims[4] = {1, 2, 3, 4};
      XLAL_CHECK_MAIN( XLALFITSArrayOpenWrite( file, "array1", XLAL_NUM_ELEM( dims ), dims, "This is a test INT4 array" ) == XLAL_SUCCESS, XLAL_EFUNC );
      size_t idx[4];
      for ( idx[0] = 0; idx[0] < dims[0]; ++idx[0] ) {
        for ( idx[1] = 0; idx[1] < dims[1]; ++idx[1] ) {
          for ( idx[2] = 0; idx[2] < dims[2]; ++idx[2] ) {
            for ( idx[3] = 0; idx[3] < dims[3]; ++idx[3] ) {
              const INT4 value = idx[0] + 2*idx[1] + 3*idx[2] + 4*idx[3];
              XLAL_CHECK_MAIN( XLALFITSArrayWriteINT4( file, idx, value ) == XLAL_SUCCESS, XLAL_EFUNC );
            }
          }
        }
      }
    }
    fprintf( stderr, "PASSED: wrote a INT4 array\n" );

    XLAL_CHECK_MAIN( XLALFITSArrayOpenWrite1( file, "array2", 14, "This is a test REAL4 array" ) == XLAL_SUCCESS, XLAL_EFUNC );
    for ( size_t i = 0; i < 14; ++i ) {
      const size_t idx[] = { i };
      const REAL4 value = 21 + 0.5*i;
      XLAL_CHECK_MAIN( XLALFITSArrayWriteREAL4( file, idx, value ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    fprintf( stderr, "PASSED: wrote a REAL4 array\n" );

    {
      const size_t m = 2, n = 3;
      gsl_matrix *GAMAT( elems, m, n );
      for ( size_t i = 0; i < m; ++i ) {
        for ( size_t j = 0; j < n; ++j ) {
          const double value = 7.0 + 2.5*i + 3.0*j;
          gsl_matrix_set( elems, i, j, value );
        }
      }
      XLAL_CHECK_MAIN( XLALFITSArrayOpenWrite2( file, "array3", m, n, "This is a test REAL8 array" ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLALFITSArrayWriteGSLMatrix( file, NULL, elems ) == XLAL_SUCCESS, XLAL_EFUNC );
      GFMAT( elems );
    }
    fprintf( stderr, "PASSED: wrote a gsl_matrix\n" );

    XLAL_CHECK_MAIN( XLALFITSTableOpenWrite( file, "table1", "This is a test table" ) == XLAL_SUCCESS, XLAL_EFUNC );
    {
      XLAL_FITS_TABLE_COLUMN_BEGIN( TestRecord );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, INT4, index ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, BOOLEAN, flag ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD_ARRAY( file, CHAR, name ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, epoch ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, pos.sky ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, pos.freq ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, pos.fkdot[0], "f1dot [Hz/s]" ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, pos.fkdot[1], "f2dot [Hz/s^2]" ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD_ARRAY( file, REAL8, values ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, COMPLEX8, phasef ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, COMPLEX16, phase ) == XLAL_SUCCESS, XLAL_EFUNC );
      {
        XLAL_FITS_TABLE_COLUMN_PTR_BEGIN( sub, TestSubRecord, 2 );
        XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, 0, INT4, n, "n1" );
        XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, 0, REAL4, v, "v1" );
        XLAL_FITS_TABLE_COLUMN_PTR_ADD_ARRAY_NAMED( file, 0, CHAR, desc, "desc1" );
        XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, 1, INT4, n, "n2" );
        XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, 1, REAL4, v, "v2" );
        XLAL_FITS_TABLE_COLUMN_PTR_ADD_ARRAY_NAMED( file, 1, CHAR, desc, "desc2" );
      }
    }
    for ( size_t i = 0; i < XLAL_NUM_ELEM( testtable ); ++i ) {
      XLAL_CHECK_MAIN( XLALFITSTableWriteRow( file, &testtable[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    fprintf( stderr, "PASSED: wrote a table\n" );

    XLAL_CHECK_MAIN( XLALFITSHeaderWriteComment( file, "%s", "This is another test comment" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLALFITSFileClose( file );
    fprintf( stderr, "PASSED: closed 'FITSFileIOTest.fits'\n" );
  }
  fprintf( stderr, "\n" );
  fflush( stderr );

  // Read the example FITS file and check data is consistent
  {
    FITSFile *file = XLALFITSFileOpenRead( "FITSFileIOTest.fits" );
    XLAL_CHECK_MAIN( file != NULL, XLAL_EFUNC );
    fprintf( stderr, "PASSED: opened 'FITSFileIOTest.fits' for reading\n" );

    {
      BOOLEAN testbool;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadBOOLEAN( file, "testbool", &testbool ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( testbool, XLAL_EFAILED, "testbool is not true" );
    }
    fprintf( stderr, "PASSED: read and verified a BOOLEAN\n" );

    {
      const INT4 testint_ref = 2345;
      INT4 testint;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadINT4( file, "testint", &testint ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( testint == testint_ref, XLAL_EFAILED, "testint = %i != %i", testint, testint_ref );
    }
    fprintf( stderr, "PASSED: read and verified a INT4\n" );

    {
      INT4 testint2_too_small;
      int errnum = 0;
      XLAL_TRY_SILENT( XLALFITSHeaderReadINT4( file, "testint2", &testint2_too_small ), errnum );
      XLAL_CHECK_MAIN( errnum = XLAL_ERANGE, XLAL_EFAILED );
    }
    {
      const INT8 testint2_ref = LAL_INT4_MAX + 6789;
      INT8 testint2;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadINT8( file, "testint2", &testint2 ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( testint2 == testint2_ref, XLAL_EFAILED, "testint = %" LAL_INT8_FORMAT " != %" LAL_INT8_FORMAT, testint2, testint2_ref );
    }
    fprintf( stderr, "PASSED: read and verified a INT8\n" );

    {
      const REAL4 testflt_ref = LAL_PI;
      REAL4 testflt;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadREAL4( file, "testflt", &testflt ) == XLAL_SUCCESS, XLAL_EFUNC );
      const REAL4 err = fabsf( testflt - testflt_ref ), err_tol = 5*FLT_EPSILON;
      XLAL_CHECK_MAIN( err < err_tol, XLAL_EFAILED, "|testflt - testflt_ref| = |%0.*g - %0.*g| = %0.*g >= %0.*g",
                       FLT_DIG, testflt, FLT_DIG, testflt_ref, FLT_DIG, err, FLT_DIG, err_tol );
    }
    fprintf( stderr, "PASSED: read and verified a REAL4\n" );

    {
      const REAL8 testdbl_ref = LAL_E;
      REAL8 testdbl;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadREAL8( file, "testdbl", &testdbl ) == XLAL_SUCCESS, XLAL_EFUNC );
      const REAL8 err = fabs( testdbl - testdbl_ref ), err_tol = 5*DBL_EPSILON;
      XLAL_CHECK_MAIN( err < err_tol, XLAL_EFAILED, "|testdbl - testdbl_ref| = |%0.*g - %0.*g| = %0.*g >= %0.*g",
                       DBL_DIG, testdbl, DBL_DIG, testdbl_ref, DBL_DIG, err, DBL_DIG, err_tol );
    }
    fprintf( stderr, "PASSED: read and verified a REAL8\n" );

    {
      const COMPLEX8 testcmp_ref = crectf( LAL_PI_2, LAL_PI_4 );
      COMPLEX8 testcmp;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadCOMPLEX8( file, "testcmp", &testcmp ) == XLAL_SUCCESS, XLAL_EFUNC );
      const REAL4 err = cabsf( testcmp - testcmp_ref ), err_tol = 5*FLT_EPSILON;
      XLAL_CHECK_MAIN( err < err_tol, XLAL_EFAILED, "|testcmp - testcmp_ref| = |(%0.*g,%0.*g) - (%0.*g,%0.*g)| = %0.*g >= %0.*g",
                       FLT_DIG, crealf( testcmp ), FLT_DIG, cimagf( testcmp ), FLT_DIG, crealf( testcmp_ref ), FLT_DIG, cimagf( testcmp_ref ), FLT_DIG, err, FLT_DIG, err_tol );
    }
    fprintf( stderr, "PASSED: read and verified a COMPLEX8\n" );

    {
      const COMPLEX16 testdblcmp_ref = crect( LAL_LOG2E, LAL_LOG10E );
      COMPLEX16 testdblcmp;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadCOMPLEX16( file, "testdblcmp", &testdblcmp ) == XLAL_SUCCESS, XLAL_EFUNC );
      const REAL8 err = cabs( testdblcmp - testdblcmp_ref ), err_tol = 5*DBL_EPSILON;
      XLAL_CHECK_MAIN( err < err_tol, XLAL_EFAILED, "|testdblcmp - testdblcmp_ref| = |(%0.*g,%0.*g) - (%0.*g,%0.*g)| = %0.*g >= %0.*g",
                       DBL_DIG, creal( testdblcmp ), DBL_DIG, cimag( testdblcmp ), DBL_DIG, creal( testdblcmp_ref ), DBL_DIG, cimag( testdblcmp_ref ), DBL_DIG, err, DBL_DIG, err_tol );
    }
    fprintf( stderr, "PASSED: read and verified a COMPLEX16\n" );

    {
      const CHAR *teststr_ref = "This is a short string";
      CHAR *teststr = NULL;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadString( file, "teststr", &teststr ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( strcmp( teststr, teststr_ref ) == 0, XLAL_EFAILED, "teststr = '%s' != '%s'", teststr, teststr_ref );
      XLALFree( teststr );
    }
    fprintf( stderr, "PASSED: read and verified a string\n" );

    {
      CHAR *longstring = NULL;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadString( file, "longstring", &longstring ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( strcmp( longstring, longstring_ref ) == 0, XLAL_EFAILED, "longstring = '%s' != '%s'", longstring, longstring_ref );
      XLALFree( longstring );
    }
    fprintf( stderr, "PASSED: read and verified a long string\n" );

    {
      LALStringVector *testsv = NULL;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadStringVector( file, "testsv", &testsv ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( testsv->length == 3, XLAL_EFAILED, "testsv->length = %u != 3", testsv->length );
      XLAL_CHECK_MAIN( strcmp( testsv->data[0], "abc" ) == 0, XLAL_EFAILED, "testsv->data[0] = '%s' != 'abc'", testsv->data[0] );
      XLAL_CHECK_MAIN( strcmp( testsv->data[1], "def" ) == 0, XLAL_EFAILED, "testsv->data[1] = '%s' != 'def'", testsv->data[1] );
      XLAL_CHECK_MAIN( strcmp( testsv->data[2], "ghij" ) == 0, XLAL_EFAILED, "testsv->data[2] = '%s' != 'ghij'", testsv->data[2] );
      XLALDestroyStringVector( testsv );
    }
    fprintf( stderr, "PASSED: read and verified a string vector\n" );

    {
      LIGOTimeGPS ref_time = LIGOTIMEGPSZERO;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadGPSTime( file, "date-obs", &ref_time ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( ref_time.gpsSeconds == 987654321 && ref_time.gpsNanoSeconds == 123456789, XLAL_EFAILED,
                       "ref_time = %" LAL_GPS_FORMAT " != 987654321.123456789", LAL_GPS_PRINT( ref_time ) );
    }
    fprintf( stderr, "PASSED: read and verified a GPS time\n" );

    {
      size_t ndim, dims[FFIO_MAX];
      XLAL_CHECK_MAIN( XLALFITSArrayOpenRead( file, "array1", &ndim, dims ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( ndim == 4, XLAL_EFAILED );
      size_t idx[4];
      for ( idx[0] = 0; idx[0] < dims[0]; ++idx[0] ) {
        for ( idx[1] = 0; idx[1] < dims[1]; ++idx[1] ) {
          for ( idx[2] = 0; idx[2] < dims[2]; ++idx[2] ) {
            for ( idx[3] = 0; idx[3] < dims[3]; ++idx[3] ) {
              const INT4 value_ref = idx[0] + 2*idx[1] + 3*idx[2] + 4*idx[3];
              INT4 value = 0;
              XLAL_CHECK_MAIN( XLALFITSArrayReadINT4( file, idx, &value ) == XLAL_SUCCESS, XLAL_EFUNC );
              XLAL_CHECK_MAIN( value == value_ref, XLAL_EFAILED, "value[%zu,%zu,%zu,%zu] = %i != %i", idx[0], idx[1], idx[2], idx[3], value, value_ref );
            }
          }
        }
      }
    }
    fprintf( stderr, "PASSED: read and verified a INT4 array\n" );

    {
      size_t dim = 0;
      XLAL_CHECK_MAIN( XLALFITSArrayOpenRead1( file, "array2", &dim ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( dim == 14, XLAL_EFAILED );
      for ( size_t i = 0; i < dim; ++i ) {
        const size_t idx[] = { i };
        const REAL4 value_ref = 21 + 0.5*i;
        REAL4 value = 0;
        XLAL_CHECK_MAIN( XLALFITSArrayReadREAL4( file, idx, &value ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLAL_CHECK_MAIN( value == value_ref, XLAL_EFAILED, "value[%zu] = %g != %g", i, value, value_ref );
      }
    }
    fprintf( stderr, "PASSED: read and verified a REAL4 array\n" );

    {
      size_t m = 0, n = 0;
      XLAL_CHECK_MAIN( XLALFITSArrayOpenRead2( file, "array3", &m, &n ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( m == 2, XLAL_EFAILED );
      XLAL_CHECK_MAIN( n == 3, XLAL_EFAILED );
      gsl_matrix *elems = NULL;
      XLAL_CHECK_MAIN( XLALFITSArrayReadGSLMatrix( file, NULL, &elems ) == XLAL_SUCCESS, XLAL_EFUNC );
      for ( size_t i = 0; i < m; ++i ) {
        for ( size_t j = 0; j < n; ++j ) {
          const double value_ref = 7.0 + 2.5*i + 3.0*j;
          double value = gsl_matrix_get( elems, i, j );
          XLAL_CHECK_MAIN( value == value_ref, XLAL_EFAILED, "value[%zu,%zu] = %g != %g", i, j, value, value_ref );
        }
      }
      GFMAT( elems );
    }
    fprintf( stderr, "PASSED: read and verified a gsl_matrix\n" );

    UINT8 nrows = LAL_UINT8_MAX;
    XLAL_CHECK_MAIN( XLALFITSTableOpenRead( file, "table1", &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( nrows == XLAL_NUM_ELEM( testtable ), XLAL_EFAILED );
    {
      XLAL_FITS_TABLE_COLUMN_BEGIN( TestRecord );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, INT4, index ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, BOOLEAN, flag ) == XLAL_SUCCESS, XLAL_EFUNC );
      {
        int errnum = 0;
        XLAL_TRY_SILENT( XLAL_FITS_TABLE_COLUMN_ADD( file, INT4, index ), errnum );
        XLAL_CHECK_MAIN( errnum = XLAL_EIO, XLAL_EFAILED );
      }
      {
        int errnum = 0;
        XLAL_TRY_SILENT( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, pos.freq ), errnum );
        XLAL_CHECK_MAIN( errnum = XLAL_EIO, XLAL_EFAILED );
      }
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD_ARRAY( file, CHAR, name ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, epoch ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, pos.sky ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, pos.freq ) == XLAL_SUCCESS, XLAL_EFUNC );
      {
        int errnum = 0;
        XLAL_TRY_SILENT( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, pos.fkdot[0], "wrong" ), errnum );
        XLAL_CHECK_MAIN( errnum = XLAL_EIO, XLAL_EFAILED );
      }
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, pos.fkdot[0], "f1dot [Hz/s]" ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, pos.fkdot[1], "f2dot [Hz/s^2]" ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD_ARRAY( file, REAL8, values ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, COMPLEX8, phasef ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLAL_FITS_TABLE_COLUMN_ADD( file, COMPLEX16, phase ) == XLAL_SUCCESS, XLAL_EFUNC );
      {
        XLAL_FITS_TABLE_COLUMN_PTR_BEGIN( sub, TestSubRecord, 2 );
        for ( size_t s = 0; s < 2; ++s ) {
          char col_name[32];
          snprintf( col_name, sizeof( col_name ), "n%zu", s + 1 );
          XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, INT4, n, col_name );
          snprintf( col_name, sizeof( col_name ), "v%zu", s + 1 );
          XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL4, v, col_name );
          snprintf( col_name, sizeof( col_name ), "desc%zu", s + 1 );
          XLAL_FITS_TABLE_COLUMN_PTR_ADD_ARRAY_NAMED( file, s, CHAR, desc, col_name );
        }
      }
    }
    TestRecord XLAL_INIT_DECL( record );
    TestSubRecord XLAL_INIT_ARRAY_DECL( record_sub, 2 );
    record.sub = record_sub;
    size_t i = 0;
    while ( nrows > 0 ) {
      XLAL_CHECK_MAIN( XLALFITSTableReadRow( file, &record, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( nrows == XLAL_NUM_ELEM( testtable ) - i - 1, XLAL_EFAILED );
      XLAL_CHECK_MAIN( record.index == testtable[i].index, XLAL_EFAILED );
      XLAL_CHECK_MAIN( record.flag == testtable[i].flag, XLAL_EFAILED );
      XLAL_CHECK_MAIN( strcmp( record.name, testtable[i].name ) == 0, XLAL_EFAILED );
      XLAL_CHECK_MAIN( XLALGPSCmp( &record.epoch, &testtable[i].epoch ) == 0, XLAL_EFAILED );
      XLAL_CHECK_MAIN( record.pos.sky == testtable[i].pos.sky, XLAL_EFAILED );
      XLAL_CHECK_MAIN( record.pos.freq == testtable[i].pos.freq, XLAL_EFAILED );
      for ( size_t j = 0; j < XLAL_NUM_ELEM( record.pos.fkdot ); ++j ) {
        XLAL_CHECK_MAIN( record.pos.fkdot[j] == testtable[i].pos.fkdot[j], XLAL_EFAILED );
      }
      for ( size_t j = 0; j < XLAL_NUM_ELEM( record.values ); ++j ) {
        XLAL_CHECK_MAIN( record.values[j] == testtable[i].values[j], XLAL_EFAILED );
      }
      XLAL_CHECK_MAIN( record.phasef == testtable[i].phasef, XLAL_EFAILED );
      XLAL_CHECK_MAIN( record.phase == testtable[i].phase, XLAL_EFAILED );
      for ( size_t s = 0; s < 2; ++s ) {
        XLAL_CHECK_MAIN( record.sub[s].n == testtable[i].sub[s].n, XLAL_EFAILED );
        XLAL_CHECK_MAIN( record.sub[s].v == testtable[i].sub[s].v, XLAL_EFAILED );
        XLAL_CHECK_MAIN( strcmp( record.sub[s].desc, testtable[i].sub[s].desc ) == 0, XLAL_EFAILED );
      }
      ++i;
    }
    XLAL_CHECK_MAIN( i == XLAL_NUM_ELEM( testtable ), XLAL_EFAILED );
    XLAL_CHECK_MAIN( XLALFITSTableReadRow( file, &record, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( nrows == 0, XLAL_EFAILED );
    XLAL_CHECK_MAIN( XLALFITSTableReadRow( file, &record, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
    fprintf( stderr, "PASSED: read and verified a table\n" );

    XLALFITSFileClose( file );
    fprintf( stderr, "PASSED: closed 'FITSFileIOTest.fits'\n" );
  }
  fprintf( stderr, "\n" );
  fflush( stderr );

  // Cleanup
  XLALDestroyUserVars();
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}

#endif // !defined(HAVE_LIBCFITSIO)
