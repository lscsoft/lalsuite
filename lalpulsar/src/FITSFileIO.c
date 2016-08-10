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
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#if defined(HAVE_LIBCFITSIO)

#include <fitsio.h>

// If fffree() is missing, use free() instead
#if !defined(HAVE_FFFREE)
#undef fits_free_memory
#define fits_free_memory(ptr, status) free(ptr)

// If ffree() is present but not declared, declare it
#elif defined(HAVE_DECL_FFFREE) && !HAVE_DECL_FFFREE
int fffree( void *, int * );
#undef fits_free_memory
#define fits_free_memory fffree

#endif // ffree()

#endif // defined(HAVE_LIBCFITSIO)

#include <lal/FITSFileIO.h>
#include <lal/LALString.h>
#include <lal/StringVector.h>
#include <lal/Date.h>
#include <lal/UserInput.h>
#include <lal/GSLHelpers.h>

#if defined(__GNUC__)
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#if defined(HAVE_LIBCFITSIO)

// Call a CFITSIO function, or print error messages on failure
#define CALL_FITS_VAL(errnum, function, ...) \
  do { \
    if (function(__VA_ARGS__, &status) != 0) { \
      CHAR CALL_FITS_buf[FLEN_STATUS + FLEN_ERRMSG]; \
      fits_get_errstatus(status, CALL_FITS_buf); \
      XLAL_PRINT_ERROR("%s() failed: %s", #function, CALL_FITS_buf); \
      while (fits_read_errmsg(CALL_FITS_buf) > 0) { \
        XLAL_PRINT_ERROR("%s() error: %s", #function, CALL_FITS_buf); \
      } \
      XLAL_ERROR_FAIL(XLAL_EIO); \
    } \
  } while(0)
#define CALL_FITS(function, ...) CALL_FITS_VAL(XLAL_EIO, function, __VA_ARGS__)

// Internal representation of a FITS file opened for reading or writing
struct tagFITSFile {
  fitsfile *ff;                         // Pointer to a CFITSIO FITS file representation
  int write;                            // True if the file is open for writing (otherwise reading)
  int hdutype;                          // Type of current HDU
  CHAR hduname[FLEN_VALUE];             // Name of current HDU
  CHAR hducomment[FLEN_COMMENT];        // Comment for name of current HDU
  struct {                              // Parameters of current array
    int naxis;                                  // Number of dimensions of array
    long naxes[FFIO_MAX];                       // Size of dimensions of array
    int bitpix;                                 // Bits per pixel of array
    int datatype;                               // Datatype of array
    size_t size;                                // Number of bytes in element of array
  } array;
  struct {                              // Parameters of current table
    int tfields;                                // Number of columns in table
    size_t noffsets[FFIO_MAX];                  // Number of nested offsets to field in table row record
    size_t offsets[FFIO_MAX][2];                // List of nested offsets to field in table row record
    CHAR ttype[FFIO_MAX][FLEN_VALUE];           // Names of columns in table
    CHAR tform[FFIO_MAX][FLEN_VALUE];           // Format of columns in table
    CHAR tunit[FFIO_MAX][FLEN_VALUE];           // Units of columns in table
    int datatype[FFIO_MAX];                     // Datatype of columns in table
    LONGLONG nelements[FFIO_MAX];               // Number of elements in columns in table
    LONGLONG nrows;                             // Number of rows in table
    LONGLONG irow;                              // Index of current row in table
  } table;
};

///
/// Write a formatted string to a FITS file using the given function
///
static int WriteFormattedString( FITSFile *file, const CHAR *format, va_list ap, int (*fcn)( fitsfile *, const char*, int* ) )
{

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( format != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( fcn != NULL, XLAL_EFAULT );

  // Format the string
  CHAR buf[4096];
  XLAL_CHECK_FAIL( vsnprintf( buf, sizeof( buf ), format, ap ) < (int)sizeof( buf ), XLAL_EERR, "Formatted string is too long" );

  // Split the string by newlines, removing any empty lines
  CHAR *p = buf;
  CHAR *t = XLALStringToken( &p, "\n", 0 );
  while ( t != NULL ) {

    // Replace any nonprintable characters with spaces
    for ( CHAR *u = t; *u != '\0'; ++u ) {
      if ( !isprint( *u ) ) {
        *u = ' ';
      }
    }

    // Write the string
    CALL_FITS( fcn, file->ff, t );

    t = XLALStringToken( &p, "\n", 0 );
  }

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

}

///
/// Extract unit from a keyword or column name
///
static int ExtractUnit( const CHAR *name_unit, CHAR *name, CHAR *unit )
{

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( name_unit != NULL, XLAL_EFAULT );

  // Copy name, possibly with unit
  strcpy( name, name_unit );

  // If name contains a unit, extract it, and remove unit and any trailing whitespace from name
  if ( unit != NULL ) {
    unit[0] = '\0';
    CHAR *unit_start = strchr( name, '[' );
    if ( unit_start != NULL ) {
      CHAR *unit_end = strchr( unit_start, ']' );
      XLAL_CHECK_FAIL( unit_end != NULL, XLAL_EINVAL, "Invalid unit in '%s'", name );
      XLAL_CHECK_FAIL( unit_end - unit_start - 1 < FLEN_VALUE, XLAL_EINVAL, "Unit in '%s' are too long", name );
      *unit_start = *unit_end = '\0';
      strcpy( unit, unit_start + 1 );
      while ( --unit_start > name && isspace( *unit_start ) ) {
        *unit_start = '\0';
      }
    }
  }

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

}

///
/// Format and check a FITS keyword
///
static int CheckFITSKeyword( const CHAR *key, CHAR *keyword, CHAR *unit )
{

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( key != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( strlen( key ) < FLEN_KEYWORD, XLAL_EINVAL, "Key '%s' is too long", key );
  XLAL_CHECK_FAIL( keyword != NULL, XLAL_EFAULT );

  // Extract unit
  XLAL_CHECK_FAIL( ExtractUnit( key, keyword, unit ) == XLAL_SUCCESS, XLAL_EINVAL );

  // Force keyword to upper case
  XLALStringToUpperCase( keyword );

  if ( strlen( keyword ) <= 8 ) {

    // Test for compliant FITS keyword
    CALL_FITS( fits_test_keyword, keyword );

  } else {

    // Format a hierarchical FITS keyword
    XLAL_CHECK_FAIL( 9 + strlen( keyword ) < FLEN_KEYWORD, XLAL_EINVAL, "Key '%s' is too long", keyword );
    CHAR buf[FLEN_KEYWORD];
    strcpy( buf, keyword );
    snprintf( keyword, FLEN_KEYWORD, "HIERARCH %s", buf );

  }

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

}

#endif // defined(HAVE_LIBCFITSIO)

void XLALFITSFileClose( FITSFile UNUSED *file )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR_VOID( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;
  if ( file != NULL ) {
    if ( file->ff != NULL ) {
      fits_close_file( file->ff, &status );
    }
    XLALFree( file );
  }

#endif // !defined(HAVE_LIBCFITSIO)
}

FITSFile *XLALFITSFileOpenWrite( const CHAR UNUSED *file_name )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR_NULL( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;
  FITSFile *file = NULL;
  CHAR *url = NULL;

  // Check input
  XLAL_CHECK_FAIL( file_name != NULL, XLAL_EFAULT );

  // Create FITSFile struct
  file = XLALCalloc( 1, sizeof( *file ) );
  XLAL_CHECK_FAIL( file != NULL, XLAL_ENOMEM );

  // Set FITSFile fields
  file->write = 1;

  // Create FITS file URL which will overwrite any existing file
  url = XLALStringAppendFmt( NULL, "!file://%s", file_name );
  XLAL_CHECK_FAIL( url != NULL, XLAL_EFUNC );

  // Open FITS file for writing
  CALL_FITS_VAL( XLAL_ESYS, fits_create_file, &file->ff, url );

  // By convention, create an empty image for the first HDU,
  // so that the correct FITS header 'SIMPLE = T' is written
  CALL_FITS( fits_create_img, file->ff, SHORT_IMG, 0, NULL );
  file->hdutype = INT_MAX;

  // Write the current system date to the FITS file
  CALL_FITS( fits_write_date, file->ff );

  XLALFree( url );
  return file;

XLAL_FAIL:

  // Delete FITS file and free memory on error
  if ( file != NULL ) {
    if ( file->ff != NULL ) {
      fits_delete_file( file->ff, &status );
    }
    XLALFree( file );
  }

  XLALFree( url );
  return NULL;

#endif // !defined(HAVE_LIBCFITSIO)
}

FITSFile *XLALFITSFileOpenRead( const CHAR UNUSED *file_name )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR_NULL( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;
  FITSFile *file = NULL;
  FILE *f = NULL;

  // Check input
  XLAL_CHECK_FAIL( file_name != NULL, XLAL_EFAULT );

  // Create FITSFile struct
  file = XLALCalloc( 1, sizeof( *file ) );
  XLAL_CHECK_FAIL( file != NULL, XLAL_ENOMEM );

  // Set FITSFile fields
  file->write = 0;

  // Check that file can be opened for reading
  errno = 0;
  f = fopen( file_name, "r" );
  if ( f == NULL ) {
    switch ( errno ) {
    case ENOENT:
      XLAL_ERROR_FAIL( XLAL_ENOENT );
    default:
      XLAL_ERROR_FAIL( XLAL_ESYS );
    }
  }

  // Open FITS file for reading
  CALL_FITS_VAL( XLAL_ESYS, fits_open_diskfile, &file->ff, file_name, READONLY );

  if ( f != NULL ) {
    fclose( f );
  }
  return file;

XLAL_FAIL:

  // Close FITS file and free memory on error
  if ( file != NULL ) {
    if ( file->ff != NULL ) {
      fits_close_file( file->ff, &status );
    }
    XLALFree( file );
  }

  if ( f != NULL ) {
    fclose( f );
  }
  return NULL;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSFileWriteHistory( FITSFile UNUSED *file, const CHAR UNUSED *format, ... )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  XLAL_CHECK_FAIL( format != NULL, XLAL_EFAULT );

  // Check that we are writing to the primary header
  int hdu_num = 0;
  fits_get_hdu_num( file->ff, &hdu_num );
  XLAL_CHECK_FAIL( hdu_num > 0, XLAL_EIO );
  XLAL_CHECK_FAIL( hdu_num == 1, XLAL_EIO, "FITS file is not at primary HDU" );

  // Write history to primary header
  va_list ap;
  va_start( ap, format );
  XLAL_CHECK_FAIL( WriteFormattedString( file, format, ap, fits_write_history ) == XLAL_SUCCESS, XLAL_EFUNC );
  va_end( ap );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSFileWriteVCSInfo( FITSFile UNUSED *file, const LALVCSInfo UNUSED *const vcs_list[] )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  XLAL_CHECK_FAIL( vcs_list != NULL, XLAL_EFAULT );

  // Write VCS information to history
  for ( size_t i = 0; vcs_list[i] != NULL; ++i ) {
    XLAL_CHECK_FAIL( XLALFITSFileWriteHistory( file, "%s version: %s\n%s commit : %s\n%s status : %s",
                                               vcs_list[i]->name, vcs_list[i]->version,
                                               vcs_list[i]->name, vcs_list[i]->vcsId,
                                               vcs_list[i]->name, vcs_list[i]->vcsStatus
                       ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSFileWriteUVarCmdLine( FITSFile UNUSED *file )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;
  CHAR *cmd_line = NULL;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );

  // Write command line to history
  cmd_line = XLALUserVarGetLog( UVAR_LOGFMT_CMDLINE );
  XLAL_CHECK_FAIL( cmd_line != NULL, XLAL_EFUNC );
  XLAL_CHECK_FAIL( XLALFITSFileWriteHistory( file, "Command line: %s", cmd_line ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALFree( cmd_line );
  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  XLALFree( cmd_line );
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteComment( FITSFile UNUSED *file, const CHAR UNUSED *format, ... )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  XLAL_CHECK_FAIL( format != NULL, XLAL_EFAULT );

  // Write comment to current header
  va_list ap;
  va_start( ap, format );
  XLAL_CHECK_FAIL( WriteFormattedString( file, format, ap, fits_write_comment ) == XLAL_SUCCESS, XLAL_EFUNC );
  va_end( ap );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteBOOLEAN( FITSFile UNUSED *file, const CHAR UNUSED *key, const BOOLEAN UNUSED value, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( comment != NULL, XLAL_EFAULT );

  // Write boolean value to current header
  CALL_FITS( fits_write_key_log, file->ff, keyword, value ? 1 : 0, comment );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderReadBOOLEAN( FITSFile UNUSED *file, const CHAR UNUSED *key, BOOLEAN UNUSED *value )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );

  // Read boolean value from current header
  int val = 0;
  CHAR comment[FLEN_COMMENT];
  CALL_FITS( fits_read_key_log, file->ff, keyword, &val, comment );
  *value = val ? 1 : 0;

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteINT4( FITSFile UNUSED *file, const CHAR UNUSED *key, const INT4 UNUSED value, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  CHAR keyword[FLEN_KEYWORD], unit[FLEN_VALUE];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, unit ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( comment != NULL, XLAL_EFAULT );

  // Write 32-bit integer value to current header
  CALL_FITS( fits_write_key_lng, file->ff, keyword, value, comment );
  CALL_FITS( fits_write_key_unit, file->ff, keyword, unit );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderReadINT4( FITSFile UNUSED *file, const CHAR UNUSED *key, INT4 UNUSED *value )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );

  // Read 32-bit integer value from current header
  LONGLONG val = 0;
  CHAR comment[FLEN_COMMENT];
  CALL_FITS( fits_read_key_lnglng, file->ff, keyword, &val, comment );
  XLAL_CHECK_FAIL( INT32_MIN <= val && val <= INT32_MAX, XLAL_ERANGE );
  *value = val;

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteINT8( FITSFile UNUSED *file, const CHAR UNUSED *key, const INT8 UNUSED value, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  CHAR keyword[FLEN_KEYWORD], unit[FLEN_VALUE];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, unit ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( comment != NULL, XLAL_EFAULT );

  // Write 64-bit integer value to current header
  CALL_FITS( fits_write_key_lng, file->ff, keyword, value, comment );
  CALL_FITS( fits_write_key_unit, file->ff, keyword, unit );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderReadINT8( FITSFile UNUSED *file, const CHAR UNUSED *key, INT8 UNUSED *value )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );

  // Read 64-bit integer value from current header
  LONGLONG val = 0;
  CHAR comment[FLEN_COMMENT];
  CALL_FITS( fits_read_key_lnglng, file->ff, keyword, &val, comment );
  XLAL_CHECK_FAIL( INT64_MIN <= val && val <= INT64_MAX, XLAL_ERANGE );
  *value = val;

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteREAL4( FITSFile UNUSED *file, const CHAR UNUSED *key, const REAL4 UNUSED value, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  CHAR keyword[FLEN_KEYWORD], unit[FLEN_VALUE];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, unit ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( comment != NULL, XLAL_EFAULT );

  // Write 32-bit floating-point value to current header
  CALL_FITS( fits_write_key_flt, file->ff, keyword, value, FLT_DIG, comment );
  CALL_FITS( fits_write_key_unit, file->ff, keyword, unit );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderReadREAL4( FITSFile UNUSED *file, const CHAR UNUSED *key, REAL4 UNUSED *value )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );

  // Read 32-bit floating-point value from current header
  CHAR comment[FLEN_COMMENT];
  CALL_FITS( fits_read_key_flt, file->ff, keyword, value, comment );

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteREAL8( FITSFile UNUSED *file, const CHAR UNUSED *key, const REAL8 UNUSED value, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  CHAR keyword[FLEN_KEYWORD], unit[FLEN_VALUE];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, unit ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( comment != NULL, XLAL_EFAULT );

  // Write 64-bit floating-point value to current header
  CALL_FITS( fits_write_key_dbl, file->ff, keyword, value, DBL_DIG, comment );
  CALL_FITS( fits_write_key_unit, file->ff, keyword, unit );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderReadREAL8( FITSFile UNUSED *file, const CHAR UNUSED *key, REAL8 UNUSED *value )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );

  // Read 64-bit floating-point value from current header
  CHAR comment[FLEN_COMMENT];
  CALL_FITS( fits_read_key_dbl, file->ff, keyword, value, comment );

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteCOMPLEX8( FITSFile UNUSED *file, const CHAR UNUSED *key, const COMPLEX8 UNUSED value, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  CHAR keyword[FLEN_KEYWORD], unit[FLEN_VALUE];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, unit ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( comment != NULL, XLAL_EFAULT );

  // Write 64-bit complex floating-point value to current header
  REAL4 val[2] = {crealf( value ), cimagf( value )};
  CALL_FITS( fits_write_key_cmp, file->ff, keyword, val, FLT_DIG, comment );
  CALL_FITS( fits_write_key_unit, file->ff, keyword, unit );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderReadCOMPLEX8( FITSFile UNUSED *file, const CHAR UNUSED *key, COMPLEX8 UNUSED *value )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );

  // Read 64-bit floating-point value from current header
  CHAR comment[FLEN_COMMENT];
  REAL4 val[2] = {0, 0};
  CALL_FITS( fits_read_key_cmp, file->ff, keyword, val, comment );
  *value = crectf( val[0], val[1] );

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteCOMPLEX16( FITSFile UNUSED *file, const CHAR UNUSED *key, const COMPLEX16 UNUSED value, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  CHAR keyword[FLEN_KEYWORD], unit[FLEN_VALUE];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, unit ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( comment != NULL, XLAL_EFAULT );

  // Write 128-bit complex floating-point value to current header
  REAL8 val[2] = {creal( value ), cimag( value )};
  CALL_FITS( fits_write_key_dblcmp, file->ff, keyword, val, DBL_DIG, comment );
  CALL_FITS( fits_write_key_unit, file->ff, keyword, unit );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderReadCOMPLEX16( FITSFile UNUSED *file, const CHAR UNUSED *key, COMPLEX16 UNUSED *value )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );

  // Read 128-bit floating-point value from current header
  CHAR comment[FLEN_COMMENT];
  REAL8 val[2] = {0, 0};
  CALL_FITS( fits_read_key_dblcmp, file->ff, keyword, val, comment );
  *value = crect( val[0], val[1] );

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteString( FITSFile UNUSED *file, const CHAR UNUSED *key, const CHAR UNUSED *value, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( comment != NULL, XLAL_EFAULT );

  // Write warning for FITS long string keyword convention
  CALL_FITS( fits_write_key_longwarn, file->ff );

  // Write string value to current header
  union {
    const CHAR *cc;
    CHAR *c;
  } bad_cast = { .cc = value };
  CALL_FITS( fits_write_key_longstr, file->ff, keyword, bad_cast.c, comment );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderReadString( FITSFile UNUSED *file, const CHAR UNUSED *key, CHAR UNUSED **value )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;
  CHAR *val = NULL;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );

  // Read string value from current header
  CHAR comment[FLEN_COMMENT];
  CALL_FITS( fits_read_key_longstr, file->ff, keyword, &val, comment );
  XLAL_CHECK_FAIL( ( *value = XLALStringDuplicate( val ) ) != NULL, XLAL_EFUNC );

  if ( val != NULL ) {
    fits_free_memory( val, &status );
  }
  return XLAL_SUCCESS;

XLAL_FAIL:
  if ( val != NULL ) {
    fits_free_memory( val, &status );
  }
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteStringVector( FITSFile UNUSED *file, const CHAR UNUSED *key, const LALStringVector UNUSED *values, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( values != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( values->length > 0, XLAL_EFAULT );
  XLAL_CHECK_FAIL( comment != NULL, XLAL_EFAULT );

  // Write string values to current header
  {
    union {
      const CHAR *cc;
      CHAR *c;
    } bad_casts[values->length];
    for ( size_t i = 0; i < values->length; ++i ) {
      bad_casts[i].cc = comment;
    }
    CALL_FITS( fits_write_keys_str, file->ff, keyword, 1, values->length, values->data, ( CHAR ** ) &bad_casts[0] );
  }

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderReadStringVector( FITSFile UNUSED *file, const CHAR UNUSED *key, LALStringVector UNUSED **values )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( values != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( *values == NULL, XLAL_EFAULT );

  // Read string values from current header
  CHAR vals[FFIO_MAX][FLEN_VALUE], *vals_ptr[FFIO_MAX];
  for ( int i = 0; i < FFIO_MAX; ++i ) {
    vals_ptr[i] = vals[i];
  }
  int nfound = 0;
  CALL_FITS( fits_read_keys_str, file->ff, keyword, 1, FFIO_MAX, vals_ptr, &nfound );
  XLAL_CHECK_FAIL( 0 < nfound, XLAL_EIO, "No items to read into string vector '%s'", keyword );
  XLAL_CHECK_FAIL( nfound <= FFIO_MAX, XLAL_EIO, "Too many items to read into string vector '%s'", keyword );
  for ( int i = 0; i < nfound; ++i ) {
    *values = XLALAppendString2Vector( *values, vals[i] );
    XLAL_CHECK_FAIL( *values != NULL, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderWriteGPSTime( FITSFile UNUSED *file, const CHAR UNUSED *key, const LIGOTimeGPS UNUSED *value, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( comment != NULL, XLAL_EFAULT );

  // Ensure GPS time is canonical
  LIGOTimeGPS gps;
  XLAL_CHECK_FAIL( XLALGPSSet( &gps, value->gpsSeconds, value->gpsNanoSeconds ) != NULL, XLAL_EFUNC );

  // Write time in UTC format to current header
  struct tm XLAL_INIT_DECL( utc );
  XLAL_CHECK_FAIL( XLALGPSToUTC( &utc, gps.gpsSeconds ) != NULL, XLAL_EFUNC );
  utc.tm_year += 1900;
  utc.tm_mon += 1;
  CHAR utc_str[FLEN_VALUE];
  CALL_FITS( fits_time2str, utc.tm_year, utc.tm_mon, utc.tm_mday, utc.tm_hour, utc.tm_min, utc.tm_sec + ( gps.gpsNanoSeconds / XLAL_BILLION_REAL8 ), 9, utc_str );
  CALL_FITS( fits_write_key_str, file->ff, keyword, utc_str, comment );

  // Write comment containing time in GPS seconds/nanoseconds format
  CHAR buf[FLEN_COMMENT];
  snprintf( buf, sizeof( buf ), "%s = GPS %" LAL_GPS_FORMAT, keyword, LAL_GPS_PRINT( gps ) );
  CALL_FITS( fits_write_comment, file->ff, buf );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSHeaderReadGPSTime( FITSFile UNUSED *file, const CHAR UNUSED *key, LIGOTimeGPS UNUSED *value )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;
  CHAR *utc_str = NULL;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  CHAR keyword[FLEN_KEYWORD];
  XLAL_CHECK_FAIL( CheckFITSKeyword( key, keyword, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_FAIL( value != NULL, XLAL_EFAULT );

  // Read time in UTC format from current header
  XLAL_CHECK_FAIL( XLALFITSHeaderReadString( file, keyword, &utc_str ) == XLAL_SUCCESS, XLAL_EFUNC );
  struct tm XLAL_INIT_DECL( utc );
  double sec = 0, int_sec = 0, frac_sec = 0;
  CALL_FITS( fits_str2time, utc_str, &utc.tm_year, &utc.tm_mon, &utc.tm_mday, &utc.tm_hour, &utc.tm_min, &sec );
  frac_sec = modf( sec, &int_sec );
  utc.tm_year -= 1900;
  utc.tm_mon -= 1;
  utc.tm_sec = lrint( int_sec );
  INT4 gps_sec = XLALUTCToGPS( &utc );
  XLAL_CHECK_FAIL( xlalErrno == 0, XLAL_EFUNC );
  INT4 gps_nsec = lrint( frac_sec * XLAL_BILLION_REAL8 );
  XLAL_CHECK_FAIL( XLALGPSSet( value, gps_sec, gps_nsec ) != NULL, XLAL_EFUNC );

  XLALFree( utc_str );
  return XLAL_SUCCESS;

XLAL_FAIL:
  XLALFree( utc_str );
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayOpenWrite( FITSFile UNUSED *file, const CHAR UNUSED *name, const size_t UNUSED ndim, const size_t UNUSED dims[], const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  XLAL_CHECK_FAIL( name != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( ndim <= FFIO_MAX, XLAL_ESIZE );
  XLAL_CHECK_FAIL( dims != NULL, XLAL_EFAULT );

  // Set current HDU
  file->hdutype = IMAGE_HDU;
  strncpy( file->hduname, name, sizeof( file->hduname ) - 1 );
  strncpy( file->hducomment, comment, sizeof( file->hducomment ) - 1 );

  // Set current HDU data
  XLAL_INIT_MEM( file->array );
  file->array.bitpix = INT_MAX;
  file->array.datatype = INT_MAX;

  // Save image dimensions
  file->array.naxis = ndim;
  for ( int i = 0; i < file->array.naxis; ++i ) {
    XLAL_CHECK_FAIL( dims[i] <= LONG_MAX, XLAL_ESIZE );
    file->array.naxes[i] = dims[i];
  }

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayOpenRead( FITSFile UNUSED *file, const CHAR UNUSED *name, size_t UNUSED *ndim, size_t UNUSED dims[] )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  XLAL_CHECK_FAIL( name != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( ndim != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( dims != NULL, XLAL_EFAULT );

  // Set current HDU
  file->hdutype = IMAGE_HDU;
  strncpy( file->hduname, name, sizeof( file->hduname ) - 1 );

  // Set current HDU data
  XLAL_INIT_MEM( file->array );
  file->array.bitpix = INT_MAX;
  file->array.datatype = INT_MAX;

  // Go to HDU with given name
  CALL_FITS( fits_movnam_hdu, file->ff, file->hdutype, file->hduname, 0 );

  // Get image dimensions
  CALL_FITS( fits_get_img_dim, file->ff, &file->array.naxis );
  *ndim = file->array.naxis;
  XLAL_CHECK_FAIL( *ndim <= FFIO_MAX, XLAL_ESIZE );
  CALL_FITS( fits_get_img_size, file->ff, file->array.naxis, file->array.naxes );
  for ( int i = 0; i < file->array.naxis; ++i ) {
    XLAL_CHECK_FAIL( 0 < file->array.naxes[i] && ( ( size_t ) file->array.naxes[i] ) <= SIZE_MAX, XLAL_ESIZE );
    dims[i] = file->array.naxes[i];
  }

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayOpenWrite1( FITSFile UNUSED *file, const CHAR UNUSED *name, const size_t UNUSED dim0, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  const size_t dims[1] = { dim0 };
  XLAL_CHECK( XLALFITSArrayOpenWrite( file, name, 1, dims, comment ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayOpenRead1( FITSFile UNUSED *file, const CHAR UNUSED *name, size_t UNUSED *dim0 )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  size_t ndim = 0, dims[FFIO_MAX];
  XLAL_CHECK( XLALFITSArrayOpenRead( file, name, &ndim, dims ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( ndim == 1, XLAL_EIO );
  if ( dim0 != NULL ) {
    *dim0 = dims[0];
  }
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayOpenWrite2( FITSFile UNUSED *file, const CHAR UNUSED *name, const size_t UNUSED dim0, const size_t UNUSED dim1, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  const size_t dims[2] = { dim0, dim1 };
  XLAL_CHECK( XLALFITSArrayOpenWrite( file, name, 2, dims, comment ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayOpenRead2( FITSFile UNUSED *file, const CHAR UNUSED *name, size_t UNUSED *dim0, size_t UNUSED *dim1 )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  size_t ndim = 0, dims[FFIO_MAX];
  XLAL_CHECK( XLALFITSArrayOpenRead( file, name, &ndim, dims ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( ndim == 2, XLAL_EIO );
  if ( dim0 != NULL ) {
    *dim0 = dims[0];
  }
  if ( dim1 != NULL ) {
    *dim1 = dims[1];
  }
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

static int UNUSED XLALFITSArrayWrite( FITSFile UNUSED *file, const size_t UNUSED idx[], const int UNUSED bitpix, const int UNUSED datatype, const void UNUSED *elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  XLAL_CHECK_FAIL( idx != 0, XLAL_EINVAL );
  XLAL_CHECK_FAIL( elem != NULL, XLAL_EFAULT );

  // Check that we are at an array
  XLAL_CHECK_FAIL( file->hdutype == IMAGE_HDU, XLAL_EIO, "Current FITS file HDU is not an array" );

  // Check for valid array indexes
  for ( int i = 0; i < file->array.naxis; ++i ) {
    XLAL_CHECK_FAIL( ( ( long ) idx[i] ) < file->array.naxes[i], XLAL_EDOM, "Array index #%i out of range (%zu >= %li)", i, idx[i], file->array.naxes[i] );
  }

  // Check for bitpix/datatype consistency
  if ( file->array.bitpix == INT_MAX || file->array.datatype == INT_MAX ) {
    file->array.bitpix = bitpix;
    file->array.datatype = datatype;

    // Create a new FITS image to store array
    CALL_FITS( fits_create_img, file->ff, bitpix, file->array.naxis, file->array.naxes );
    CALL_FITS( fits_write_key_str, file->ff, "HDUNAME", file->hduname, file->hducomment );

  } else {
    XLAL_CHECK_FAIL( file->array.bitpix == bitpix && file->array.datatype == datatype, XLAL_EINVAL, "Inconsistent use of %s...() functions", __func__ );
  }

  // Write array element
  long fpixel[FFIO_MAX];
  for ( int i = 0; i < file->array.naxis; ++i ) {
    fpixel[i] = 1 + idx[i];
  }
  union {
    const void *cv;
    CHAR *c;
  } bad_cast = { .cv = elem };
  CALL_FITS( fits_write_pix, file->ff, file->array.datatype, fpixel, 1, bad_cast.c );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

static int UNUSED XLALFITSArrayRead( FITSFile UNUSED *file, const size_t UNUSED idx[], const int UNUSED bitpix, const int UNUSED datatype, void UNUSED *elem, void UNUSED *nulelem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  XLAL_CHECK_FAIL( idx != 0, XLAL_EINVAL );
  XLAL_CHECK_FAIL( elem != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( nulelem != NULL, XLAL_EFAULT );

  // Check that we are at an array
  XLAL_CHECK_FAIL( file->hdutype == IMAGE_HDU, XLAL_EIO, "Current FITS file HDU is not an array" );

  // Check for valid array indexes
  for ( int i = 0; i < file->array.naxis; ++i ) {
    XLAL_CHECK_FAIL( ( ( long ) idx[i] ) < file->array.naxes[i], XLAL_EDOM, "Array index #%i out of range (%zu >= %li)", i, idx[i], file->array.naxes[i] );
  }

  // Check for bitpix/datatype consistency
  if ( file->array.bitpix == INT_MAX || file->array.datatype == INT_MAX ) {
    file->array.bitpix = bitpix;
    file->array.datatype = datatype;
  } else {
    XLAL_CHECK_FAIL( file->array.bitpix == bitpix && file->array.datatype == datatype, XLAL_EINVAL, "Inconsistent use of %s...() functions", __func__ );
  }

  // Read array elememt
  long fpixel[FFIO_MAX];
  for ( int i = 0; i < file->array.naxis; ++i ) {
    fpixel[i] = 1 + idx[i];
  }
  int anynul = 0;
  CALL_FITS( fits_read_pix, file->ff, datatype, fpixel, 1, nulelem, elem, &anynul );

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayWriteUINT2( FITSFile UNUSED *file, const size_t UNUSED idx[], const UINT2 UNUSED elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  const unsigned short e = elem;
  XLAL_CHECK( XLALFITSArrayWrite( file, idx, USHORT_IMG, TUSHORT, &e ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayReadUINT2( FITSFile UNUSED *file, const size_t UNUSED idx[], UINT2 UNUSED *elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  unsigned short e = 0, ne = 0;
  XLAL_CHECK( XLALFITSArrayRead( file, idx, USHORT_IMG, TUSHORT, &e, &ne ) == XLAL_SUCCESS, XLAL_EFUNC );
  *elem = e;
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayWriteUINT4( FITSFile UNUSED *file, const size_t UNUSED idx[], const UINT4 UNUSED elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  const unsigned long e = elem;
  XLAL_CHECK( XLALFITSArrayWrite( file, idx, ULONG_IMG, TULONG, &e ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayReadUINT4( FITSFile UNUSED *file, const size_t UNUSED idx[], UINT4 UNUSED *elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  unsigned long e = 0, ne = 0;
  XLAL_CHECK( XLALFITSArrayRead( file, idx, ULONG_IMG, TULONG, &e, &ne ) == XLAL_SUCCESS, XLAL_EFUNC );
  *elem = e;
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayWriteINT2( FITSFile UNUSED *file, const size_t UNUSED idx[], const INT2 UNUSED elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  const short e = elem;
  XLAL_CHECK( XLALFITSArrayWrite( file, idx, SHORT_IMG, TSHORT, &e ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayReadINT2( FITSFile UNUSED *file, const size_t UNUSED idx[], INT2 UNUSED *elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  short e = 0, ne = 0;
  XLAL_CHECK( XLALFITSArrayRead( file, idx, SHORT_IMG, TSHORT, &e, &ne ) == XLAL_SUCCESS, XLAL_EFUNC );
  *elem = e;
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayWriteINT4( FITSFile UNUSED *file, const size_t UNUSED idx[], const INT4 UNUSED elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  const long e = elem;
  XLAL_CHECK( XLALFITSArrayWrite( file, idx, LONG_IMG, TLONG, &e ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayReadINT4( FITSFile UNUSED *file, const size_t UNUSED idx[], INT4 UNUSED *elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  long e = 0, ne = 0;
  XLAL_CHECK( XLALFITSArrayRead( file, idx, LONG_IMG, TLONG, &e, &ne ) == XLAL_SUCCESS, XLAL_EFUNC );
  *elem = e;
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayWriteREAL4( FITSFile UNUSED *file, const size_t UNUSED idx[], const REAL4 UNUSED elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  const float e = elem;
  XLAL_CHECK( XLALFITSArrayWrite( file, idx, FLOAT_IMG, TFLOAT, &e ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayReadREAL4( FITSFile UNUSED *file, const size_t UNUSED idx[], REAL4 UNUSED *elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  float e = 0, ne = 0;
  XLAL_CHECK( XLALFITSArrayRead( file, idx, FLOAT_IMG, TFLOAT, &e, &ne ) == XLAL_SUCCESS, XLAL_EFUNC );
  *elem = e;
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayWriteREAL8( FITSFile UNUSED *file, const size_t UNUSED idx[], const REAL8 UNUSED elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  const double e = elem;
  XLAL_CHECK( XLALFITSArrayWrite( file, idx, DOUBLE_IMG, TDOUBLE, &e ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayReadREAL8( FITSFile UNUSED *file, const size_t UNUSED idx[], REAL8 UNUSED *elem )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  double e = 0, ne = 0;
  XLAL_CHECK( XLALFITSArrayRead( file, idx, DOUBLE_IMG, TDOUBLE, &e, &ne ) == XLAL_SUCCESS, XLAL_EFUNC );
  *elem = e;
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayWriteGSLMatrix( FITSFile UNUSED *file, const size_t UNUSED idx[], const gsl_matrix UNUSED *elems )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  XLAL_CHECK_FAIL( elems != NULL, XLAL_EFAULT );

  // Check that we are at an array of sufficient size
  XLAL_CHECK_FAIL( file->hdutype == IMAGE_HDU, XLAL_EIO, "Current FITS file HDU is not an array" );
  XLAL_CHECK_FAIL( file->array.naxis >= 2, XLAL_EINVAL, "Array must have at least 2 dimensions" );
  XLAL_CHECK_FAIL( elems->size1 == ( size_t ) file->array.naxes[0], XLAL_EINVAL, "Number of 'elems' rows (%zu) does not match array dimension 0 (%li)", elems->size1, file->array.naxes[0] );
  XLAL_CHECK_FAIL( elems->size2 == ( size_t ) file->array.naxes[1], XLAL_EINVAL, "Number of 'elems' rows (%zu) does not match array dimension 1 (%li)", elems->size2, file->array.naxes[1] );

  // Copy index vector, if given
  size_t XLAL_INIT_ARRAY_DECL( i, FFIO_MAX );
  if ( idx != NULL ) {
    memcpy( i, idx, file->array.naxis * sizeof( i[0] ) );
  }

  // Write GSL matrix elements to first 2 dimensions
  for ( i[0] = 0; i[0] < ( size_t ) file->array.naxes[0]; ++i[0] ) {
    for ( i[1] = 0; i[1] < ( size_t ) file->array.naxes[1]; ++i[1] ) {
      const REAL8 elem = gsl_matrix_get( elems, i[0], i[1] );
      XLAL_CHECK_FAIL( XLALFITSArrayWriteREAL8( file, i, elem ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSArrayReadGSLMatrix( FITSFile UNUSED *file, const size_t UNUSED idx[], gsl_matrix UNUSED **elems )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  XLAL_CHECK_FAIL( elems != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( *elems == NULL, XLAL_EINVAL );

  // Check that we are at an array of sufficient size
  XLAL_CHECK_FAIL( file->hdutype == IMAGE_HDU, XLAL_EIO, "Current FITS file HDU is not an array" );
  XLAL_CHECK_FAIL( file->array.naxis >= 2, XLAL_EINVAL, "Array must have at least 2 dimensions" );

  // Create GSL matrix
  GAMAT( *elems, file->array.naxes[0], file->array.naxes[1] );

  // Copy index vector, if given
  size_t XLAL_INIT_ARRAY_DECL( i, FFIO_MAX );
  if ( idx != NULL ) {
    memcpy( i, idx, file->array.naxis * sizeof( i[0] ) );
  }

  // Read GSL matrix elements from first 2 dimensions
  for ( i[0] = 0; i[0] < ( size_t ) file->array.naxes[0]; ++i[0] ) {
    for ( i[1] = 0; i[1] < ( size_t ) file->array.naxes[1]; ++i[1] ) {
      REAL8 elem = 0;
      XLAL_CHECK_FAIL( XLALFITSArrayReadREAL8( file, i, &elem ) == XLAL_SUCCESS, XLAL_EFUNC );
      gsl_matrix_set( *elems, i[0], i[1], elem );
    }
  }

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableOpenWrite( FITSFile UNUSED *file, const CHAR UNUSED *name, const CHAR UNUSED *comment )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  XLAL_CHECK_FAIL( name != NULL, XLAL_EFAULT );

  // Set current HDU
  file->hdutype = BINARY_TBL;
  strncpy( file->hduname, name, sizeof( file->hduname ) - 1 );
  strncpy( file->hducomment, comment, sizeof( file->hducomment ) - 1 );

  // Set current HDU data
  XLAL_INIT_MEM( file->table );

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableOpenRead( FITSFile UNUSED *file, const CHAR UNUSED *name, UINT8 UNUSED *nrows )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  XLAL_CHECK_FAIL( name != NULL, XLAL_EFAULT );

  // Set current HDU
  file->hdutype = BINARY_TBL;
  strncpy( file->hduname, name, sizeof( file->hduname ) - 1 );

  // Set current HDU data
  XLAL_INIT_MEM( file->table );

  // Go to HDU with given name
  CALL_FITS( fits_movnam_hdu, file->ff, file->hdutype, file->hduname, 0 );

  // Get number of table rows
  CALL_FITS( fits_get_num_rowsll, file->ff, &file->table.nrows );
  if ( nrows != NULL ) {
    *nrows = file->table.nrows;
  }

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

static int UNUSED XLALFITSTableColumnAdd( FITSFile UNUSED *file, const CHAR UNUSED *name, const CHAR UNUSED *unit, const size_t UNUSED noffsets, const size_t UNUSED offsets[2], const void UNUSED *record, const size_t UNUSED record_size, const void UNUSED *field, const size_t UNUSED field_size, const size_t UNUSED elem_size, const int UNUSED typechar, const int UNUSED datatype )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Save previous number of table columns
  const int save_tfields = ( file != NULL ) ? file->table.tfields : 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( name != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( unit != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( 0 < noffsets && noffsets <= 2, XLAL_EINVAL );
  XLAL_CHECK_FAIL( offsets != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( record != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( record_size > 0, XLAL_EINVAL );
  XLAL_CHECK_FAIL( field != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( field_size > 0, XLAL_EINVAL );
  XLAL_CHECK_FAIL( field_size < 128, XLAL_EINVAL, "Record field is too long" );
  XLAL_CHECK_FAIL( ( ( intptr_t ) field ) >= ( ( intptr_t ) record ), XLAL_EINVAL, "Invalid field pointer" );
  XLAL_CHECK_FAIL( ( ( intptr_t ) field ) + field_size <= ( ( intptr_t ) record ) + record_size, XLAL_EINVAL, "Invalid field pointer" );
  XLAL_CHECK_FAIL( elem_size > 0, XLAL_EINVAL );

  // Check that we are at a table
  XLAL_CHECK_FAIL( file->hdutype == BINARY_TBL, XLAL_EIO, "Current FITS file HDU is not a table" );

  // Move to next column
  XLAL_CHECK_FAIL( file->table.tfields <= FFIO_MAX, XLAL_ESIZE );
  const int i = file->table.tfields++;

  // Store field offsets
  file->table.noffsets[i] = noffsets;
  memcpy( file->table.offsets[i], offsets, sizeof( file->table.offsets[i] ) );
  file->table.offsets[i][noffsets - 1] = ( size_t )( ( ( intptr_t ) field ) - ( ( intptr_t ) record ) );

  // Copy column name
  snprintf( file->table.ttype[i], sizeof( file->table.ttype[i] ), "%s", name );

  // Copy column unit
  snprintf( file->table.tunit[i], sizeof( file->table.tunit[i] ), "%s", unit );

  // Create column format
  snprintf( file->table.tform[i], sizeof( file->table.tform[i] ), "%zu%c", field_size / elem_size, typechar );

  // Store column datatype and number of elements
  file->table.datatype[i] = datatype;
  file->table.nelements[i] = ( datatype == TSTRING ) ? 1 : field_size / elem_size;

  if ( !file->write ) {

    // Search for and verify existing column name
    CHAR ttype_chk[FLEN_VALUE];
    int i_chk = 0;
    CALL_FITS( fits_get_colname, file->ff, CASEINSEN, file->table.ttype[i], ttype_chk, &i_chk );
    XLAL_CHECK_FAIL( i_chk > 0, XLAL_EIO, "Column '%s' does not exist in FITS file", file->table.ttype[i] );
    XLAL_CHECK_FAIL( i_chk == 1 + i, XLAL_EIO, "Column '%s' is in a diferent position in FITS file", file->table.ttype[i] );
    XLAL_CHECK_FAIL( strcmp( ttype_chk, file->table.ttype[i] ) == 0, XLAL_EIO, "Inconsistent column #%i name '%s'; should be '%s'", i, ttype_chk, file->table.ttype[i] );

    // Verify existing column format
    int datatype_chk = 0;
    long repeat_chk = 0, width_chk = 0;
    CALL_FITS( fits_get_eqcoltype, file->ff, 1 + i, &datatype_chk, &repeat_chk, &width_chk );
    XLAL_CHECK_FAIL( datatype_chk == datatype, XLAL_EIO, "Inconsistent column #%i datatype %i; should be %i", i, datatype_chk, datatype );
    const size_t elem_size_chk = ( datatype_chk == TSTRING ) ? 1 : width_chk;
    const size_t field_size_chk = repeat_chk * elem_size_chk;
    XLAL_CHECK_FAIL( field_size_chk == field_size, XLAL_EIO, "Inconsistent column #%i field size %zu; should be %zu", i, field_size_chk, field_size );
    XLAL_CHECK_FAIL( elem_size_chk == elem_size, XLAL_EIO, "Inconsistent column #%i element size %zu; should be %zu", i, elem_size_chk, elem_size );

  }

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Restore previous number of table columns
  if ( file != NULL ) {
    file->table.tfields = save_tfields;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableColumnAddBOOLEAN( FITSFile UNUSED *file, const CHAR UNUSED *col_name, const size_t UNUSED noffsets, const size_t UNUSED offsets[2], const void UNUSED *record, const size_t UNUSED record_size, const BOOLEAN UNUSED *field, const size_t UNUSED field_size )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  XLAL_CHECK( col_name != NULL, XLAL_EFAULT );
  XLAL_CHECK( strlen( col_name ) < FLEN_VALUE, XLAL_EINVAL, "Column name '%s' is too long", col_name );
  XLAL_CHECK( XLALFITSTableColumnAdd( file, col_name, "", noffsets, offsets, record, record_size, field, field_size, sizeof( BOOLEAN ), 'L', TLOGICAL ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableColumnAddINT2( FITSFile UNUSED *file, const CHAR UNUSED *col_name, const size_t UNUSED noffsets, const size_t UNUSED offsets[2], const void UNUSED *record, const size_t UNUSED record_size, const INT2 UNUSED *field, const size_t UNUSED field_size )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  XLAL_CHECK( col_name != NULL, XLAL_EFAULT );
  XLAL_CHECK( strlen( col_name ) < FLEN_VALUE, XLAL_EINVAL, "Column name '%s' is too long", col_name );
  CHAR name[FLEN_VALUE], unit[FLEN_VALUE];
  XLAL_CHECK( ExtractUnit( col_name, name, unit ) == XLAL_SUCCESS, XLAL_EINVAL );
  XLAL_CHECK( XLALFITSTableColumnAdd( file, name, unit, noffsets, offsets, record, record_size, field, field_size, sizeof( INT2 ), 'I', TSHORT ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableColumnAddINT4( FITSFile UNUSED *file, const CHAR UNUSED *col_name, const size_t UNUSED noffsets, const size_t UNUSED offsets[2], const void UNUSED *record, const size_t UNUSED record_size, const INT4 UNUSED *field, const size_t UNUSED field_size )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  XLAL_CHECK( col_name != NULL, XLAL_EFAULT );
  XLAL_CHECK( strlen( col_name ) < FLEN_VALUE, XLAL_EINVAL, "Column name '%s' is too long", col_name );
  CHAR name[FLEN_VALUE], unit[FLEN_VALUE];
  XLAL_CHECK( ExtractUnit( col_name, name, unit ) == XLAL_SUCCESS, XLAL_EINVAL );
  XLAL_CHECK( XLALFITSTableColumnAdd( file, name, unit, noffsets, offsets, record, record_size, field, field_size, sizeof( INT4 ), 'J', TINT32BIT ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableColumnAddREAL4( FITSFile UNUSED *file, const CHAR UNUSED *col_name, const size_t UNUSED noffsets, const size_t UNUSED offsets[2], const void UNUSED *record, const size_t UNUSED record_size, const REAL4 UNUSED *field, const size_t UNUSED field_size )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  XLAL_CHECK( col_name != NULL, XLAL_EFAULT );
  XLAL_CHECK( strlen( col_name ) < FLEN_VALUE, XLAL_EINVAL, "Column name '%s' is too long", col_name );
  CHAR name[FLEN_VALUE], unit[FLEN_VALUE];
  XLAL_CHECK( ExtractUnit( col_name, name, unit ) == XLAL_SUCCESS, XLAL_EINVAL );
  XLAL_CHECK( XLALFITSTableColumnAdd( file, name, unit, noffsets, offsets, record, record_size, field, field_size, sizeof( REAL4 ), 'E', TFLOAT ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableColumnAddREAL8( FITSFile UNUSED *file, const CHAR UNUSED *col_name, const size_t UNUSED noffsets, const size_t UNUSED offsets[2], const void UNUSED *record, const size_t UNUSED record_size, const REAL8 UNUSED *field, const size_t UNUSED field_size )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  XLAL_CHECK( col_name != NULL, XLAL_EFAULT );
  XLAL_CHECK( strlen( col_name ) < FLEN_VALUE, XLAL_EINVAL, "Column name '%s' is too long", col_name );
  CHAR name[FLEN_VALUE], unit[FLEN_VALUE];
  XLAL_CHECK( ExtractUnit( col_name, name, unit ) == XLAL_SUCCESS, XLAL_EINVAL );
  XLAL_CHECK( XLALFITSTableColumnAdd( file, name, unit, noffsets, offsets, record, record_size, field, field_size, sizeof( REAL8 ), 'D', TDOUBLE ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableColumnAddCOMPLEX8( FITSFile UNUSED *file, const CHAR UNUSED *col_name, const size_t UNUSED noffsets, const size_t UNUSED offsets[2], const void UNUSED *record, const size_t UNUSED record_size, const COMPLEX8 UNUSED *field, const size_t UNUSED field_size )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  XLAL_CHECK( col_name != NULL, XLAL_EFAULT );
  XLAL_CHECK( strlen( col_name ) < FLEN_VALUE, XLAL_EINVAL, "Column name '%s' is too long", col_name );
  CHAR name[FLEN_VALUE], unit[FLEN_VALUE];
  XLAL_CHECK( ExtractUnit( col_name, name, unit ) == XLAL_SUCCESS, XLAL_EINVAL );
  XLAL_CHECK( XLALFITSTableColumnAdd( file, name, unit, noffsets, offsets, record, record_size, field, field_size, sizeof( COMPLEX8 ), 'C', TCOMPLEX ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableColumnAddCOMPLEX16( FITSFile UNUSED *file, const CHAR UNUSED *col_name, const size_t UNUSED noffsets, const size_t UNUSED offsets[2], const void UNUSED *record, const size_t UNUSED record_size, const COMPLEX16 UNUSED *field, const size_t UNUSED field_size )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  XLAL_CHECK( col_name != NULL, XLAL_EFAULT );
  XLAL_CHECK( strlen( col_name ) < FLEN_VALUE, XLAL_EINVAL, "Column name '%s' is too long", col_name );
  CHAR name[FLEN_VALUE], unit[FLEN_VALUE];
  XLAL_CHECK( ExtractUnit( col_name, name, unit ) == XLAL_SUCCESS, XLAL_EINVAL );
  XLAL_CHECK( XLALFITSTableColumnAdd( file, name, unit, noffsets, offsets, record, record_size, field, field_size, sizeof( COMPLEX16 ), 'M', TDBLCOMPLEX ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableColumnAddCHAR( FITSFile UNUSED *file, const CHAR UNUSED *col_name, const size_t UNUSED noffsets, const size_t UNUSED offsets[2], const void UNUSED *record, const size_t UNUSED record_size, const void UNUSED *field, const size_t UNUSED field_size )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  XLAL_CHECK( col_name != NULL, XLAL_EFAULT );
  XLAL_CHECK( strlen( col_name ) < FLEN_VALUE, XLAL_EINVAL, "Column name '%s' is too long", col_name );
  XLAL_CHECK( XLALFITSTableColumnAdd( file, col_name, "", noffsets, offsets, record, record_size, field, field_size, sizeof( CHAR ), 'A', TSTRING ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableColumnAddGPSTime( FITSFile UNUSED *file, const CHAR UNUSED *col_name, const size_t UNUSED noffsets, const size_t UNUSED offsets[2], const void UNUSED *record, const size_t UNUSED record_size, const LIGOTimeGPS UNUSED *field, const size_t UNUSED field_size )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  XLAL_CHECK( col_name != NULL, XLAL_EFAULT );
  XLAL_CHECK( strlen( col_name ) + 3 < FLEN_VALUE, XLAL_EINVAL, "Column name '%s' is too long", col_name );
  XLAL_CHECK( field_size == sizeof( LIGOTimeGPS ), XLAL_EINVAL, "Array of GPS times is not supported" );
  CHAR name[FLEN_VALUE];
  snprintf( name, sizeof( name ), "%s_s", col_name );
  XLAL_CHECK( XLALFITSTableColumnAdd( file, name, "s", noffsets, offsets, record, record_size, &( field->gpsSeconds ), sizeof( field->gpsSeconds ), sizeof( field->gpsSeconds ), 'J', TINT32BIT ) == XLAL_SUCCESS, XLAL_EFUNC );
  snprintf( name, sizeof( name ), "%s_ns", col_name );
  XLAL_CHECK( XLALFITSTableColumnAdd( file, name, "ns", noffsets, offsets, record, record_size, &( field->gpsNanoSeconds ), sizeof( field->gpsNanoSeconds ), sizeof( field->gpsNanoSeconds ), 'J', TINT32BIT ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableWriteRow( FITSFile UNUSED *file, const void UNUSED *record )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( file->write, XLAL_EINVAL, "FITS file is not open for writing" );
  XLAL_CHECK_FAIL( record != NULL, XLAL_EFAULT );

  // Check that we are at a table
  XLAL_CHECK_FAIL( file->hdutype == BINARY_TBL, XLAL_EIO, "Current FITS file HDU is not a table" );

  // Create new table if required
  if ( file->table.irow == 0 ) {
    CHAR *ttype_ptr[FFIO_MAX], *tform_ptr[FFIO_MAX], *tunit_ptr[FFIO_MAX];
    for ( int i = 0; i < file->table.tfields; ++i ) {
      ttype_ptr[i] = file->table.ttype[i];
      tform_ptr[i] = file->table.tform[i];
      tunit_ptr[i] = file->table.tunit[i];
    }
    CALL_FITS( fits_create_tbl, file->ff, file->hdutype, 0, file->table.tfields, ttype_ptr, tform_ptr, tunit_ptr, NULL );
    CALL_FITS( fits_write_key_str, file->ff, "HDUNAME", file->hduname, file->hducomment );
  }

  // Advance to next row
  ++file->table.irow;

  // Write next table row
  for ( int i = 0; i < file->table.tfields; ++i ) {
    union {
      const void *cv;
      void *v;
    } bad_cast = { .cv = record };
    void *value = bad_cast.v;
    for ( size_t n = 0; n < file->table.noffsets[i]; ++n ) {
      if ( n > 0 ) {
        value = *( ( void ** ) value );
      }
      value = ( void * )( ( ( intptr_t ) value ) + file->table.offsets[i][n] );
    }
    void *pvalue = ( file->table.datatype[i] == TSTRING ) ? ( void * ) &value : value;
    CALL_FITS( fits_write_col, file->ff, file->table.datatype[i], 1 + i, file->table.irow, 1, file->table.nelements[i], pvalue );
  }

  return XLAL_SUCCESS;

XLAL_FAIL:

  // Delete FITS file on error
  if ( file != NULL && file->ff != NULL ) {
    fits_delete_file( file->ff, &status );
    file->ff = NULL;
  }

  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}

int XLALFITSTableReadRow( FITSFile UNUSED *file, void UNUSED *record, UINT8 UNUSED *rem_nrows )
{
#if !defined(HAVE_LIBCFITSIO)
  XLAL_ERROR( XLAL_EFAILED, "CFITSIO is not available" );
#else // defined(HAVE_LIBCFITSIO)

  int UNUSED status = 0;

  // Check input
  XLAL_CHECK_FAIL( file != NULL, XLAL_EFAULT );
  XLAL_CHECK_FAIL( !file->write, XLAL_EINVAL, "FITS file is not open for reading" );
  XLAL_CHECK_FAIL( record != NULL, XLAL_EFAULT );

  // Check that we are at a table
  XLAL_CHECK_FAIL( file->hdutype == BINARY_TBL, XLAL_EIO, "Current FITS file HDU is not a table" );

  // Return if there are no more rows
  if ( file->table.irow == file->table.nrows ) {
    return XLAL_SUCCESS;
  }

  // Advance to next row, and return number of remaining rows
  ++file->table.irow;
  if ( rem_nrows != NULL ) {
    *rem_nrows = file->table.nrows - file->table.irow;
  }

  // Read next table row
  for ( int i = 0; i < file->table.tfields; ++i ) {
    void *value = record;
    for ( size_t n = 0; n < file->table.noffsets[i]; ++n ) {
      if ( n > 0 ) {
        value = *( ( void ** ) value );
      }
      value = ( void * )( ( ( intptr_t ) value ) + file->table.offsets[i][n] );
    }
    void *pvalue = ( file->table.datatype[i] == TSTRING ) ? ( void * ) &value : value;
    CALL_FITS( fits_read_col, file->ff, file->table.datatype[i], 1 + i, file->table.irow, 1, file->table.nelements[i], NULL, pvalue, NULL );
  }

  return XLAL_SUCCESS;

XLAL_FAIL:
  return XLAL_FAILURE;

#endif // !defined(HAVE_LIBCFITSIO)
}
