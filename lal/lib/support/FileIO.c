/*
 *  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Reinhard Prix
 *  Copyright (C) 2012--2015, 2021 Karl Wette
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/* for realpath() for whereami.c on Linux */
#define _GNU_SOURCE

#include <config.h>

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_GLOB_H
#include <glob.h>
#endif

#include <zlib.h>
#define ZLIB_ENABLED

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALString.h>
#include <lal/StringInput.h>
#include <lal/FileIO.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* do not export functions from whereami.c */
#define WAI_FUNCSPEC static UNUSED

#include "whereami.c"

struct tagLALFILE {
  int compression;
  void *fp;
};

LALFILE *lalstdin( void )
{
  static LALFILE _lalstdin;
  _lalstdin.fp = (void*)stdin;
  return &_lalstdin;
}
LALFILE *lalstdout( void )
{
  static LALFILE _lalstdout;
  _lalstdout.fp = (void*)stdout;
  return &_lalstdout;
}
LALFILE *lalstderr( void )
{
  static LALFILE _lalstderr;
  _lalstderr.fp = (void*)stderr;
  return &_lalstderr;
}


/** \cond DONT_DOXYGEN */
/* should be only called from XLAL_FILE_RESOLVE_PATH() */
char *
XLALFileResolvePathLong ( const char *fname, const char *fallbackpath )
{
  XLAL_CHECK_NULL ( fname != NULL, XLAL_EINVAL );

  const UINT4 fname_len = strlen ( fname );

  XLALPrintInfo ( "%s(): trying to resolve fname='%s' ...\n", __func__, fname );

  if ( strchr ( fname, '/' ) != NULL )	// any kind of (absolute or relative) path is given -> use only that
    {
      FILE *tmp;
      XLALPrintInfo ( "%s(): strategy: absolute/relative path\n", __func__ );
      XLALPrintInfo ( "%s(): trying '%s' -> '%s' ...\n", __func__, fname, fname );
      if ( (tmp = LALFopen ( fname, "rb" )) != NULL ) {
        LALFclose ( tmp );
        XLALPrintInfo ( "%s(): success '%s' -> '%s'\n", __func__, fname, fname );
        return XLALStringDuplicate ( fname );
      } // if found
      else {
        XLALPrintInfo ( "%s(): failed '%s'\n", __func__, fname );
        return NULL;
      } // if not found

    } // end: if path given
  else	// if pure filename given: try 1) local directory, then 2) scan LAL_DATA_PATH, then 3) scan 'fallbackpath' (if given)
    {
      FILE *tmp;
      char *resolveFname = NULL;

      // ----- Strategy 1: try local directory explicitly
      XLALPrintInfo ( "%s(): strategy: local directory\n", __func__ );
      XLAL_CHECK_NULL ( (resolveFname = XLALMalloc ( fname_len + 2 + 1 )) != NULL, XLAL_ENOMEM );
      sprintf ( resolveFname, "./%s", fname );
      XLALPrintInfo ( "%s(): trying '%s' -> '%s' ...\n", __func__, fname, resolveFname );
      if ( (tmp = LALFopen ( resolveFname, "rb" )) != NULL ) {
        LALFclose ( tmp );
        XLALPrintInfo ( "%s(): success '%s' -> '%s'\n", __func__, fname, resolveFname );
        return resolveFname;
      } // if found

      // ----- Strategy 2: scan LAL_DATA_PATH
      const char *lal_data_path = getenv( "LAL_DATA_PATH" );
      if ( ! ( lal_data_path && (strlen (lal_data_path ) > 0 ) ) )
        {
          XLALPrintInfo ( "%s(): skip strategy LAL_DATA_PATH: environment variable not set\n", __func__ );
        }
      else
        {
          XLALPrintInfo ( "%s(): strategy: LAL_DATA_PATH='%s'\n", __func__, lal_data_path );
          TokenList *subPaths = NULL;
          XLAL_CHECK_NULL ( XLALCreateTokenList ( &subPaths, lal_data_path, ":" ) == XLAL_SUCCESS, XLAL_EFUNC );
          for ( UINT4 i = 0; i < subPaths->nTokens; i ++ )
            {
              const char *subPath_i = subPaths->tokens[i];
              XLAL_CHECK_NULL ( (resolveFname = XLALRealloc ( resolveFname, strlen(subPath_i) + 1 + fname_len + 1 )) != NULL, XLAL_ENOMEM );
              sprintf ( resolveFname, "%s/%s", subPath_i, fname );
              XLALPrintInfo ( "%s(): trying '%s' -> '%s' ...\n", __func__, fname, resolveFname );
              if ( (tmp = LALFopen ( resolveFname, "rb" )) != NULL ) {
                LALFclose ( tmp );
                XLALDestroyTokenList ( subPaths );
                XLALPrintInfo ( "%s(): success '%s' -> '%s'\n", __func__, fname, resolveFname );
                return resolveFname;
              } // if found

            } // for i < nTokens

          XLALDestroyTokenList ( subPaths );

        } // if LAL_DATA_PATH given

      // ----- Strategy 3: scan 'fallbackpath', if given
      if ( ! fallbackpath ) {
        XLALPrintInfo ( "%s(): skip strategy: fallbackpath=NULL\n", __func__ );
      } else {
        XLALPrintInfo ( "%s(): strategy: fallbackpath='%s'\n", __func__, fallbackpath );

        // get location of liblalsupport.so
        char *module_path = NULL;
        int module_path_length = wai_getModulePath ( NULL, 0, NULL );
        int module_dirname_length = 0;
        if ( module_path_length > 0 ) {
          module_path = XLALMalloc ( module_path_length + 1 );
          XLAL_CHECK_NULL ( module_path != NULL, XLAL_ENOMEM );
          XLAL_CHECK_NULL ( wai_getModulePath ( module_path, module_path_length, &module_dirname_length ) > 0, XLAL_EERR );
          module_path[module_path_length] = '\0';
          XLAL_CHECK_NULL ( module_dirname_length >= 0, XLAL_EERR );
          XLAL_CHECK_NULL ( module_dirname_length < module_path_length, XLAL_EERR );
          module_path[module_dirname_length] = '\0';
        }

        // search 'fallbackpath'; relative paths are considered relative to location of liblalsupport.so (if known)
        TokenList *subPaths = NULL;
        XLAL_CHECK_NULL ( XLALCreateTokenList ( &subPaths, fallbackpath, ":" ) == XLAL_SUCCESS, XLAL_EFUNC );
        for ( UINT4 i = 0; i < subPaths->nTokens; i ++ )
          {
            const char *subPath_i = subPaths->tokens[i];
            const UINT4 subPath_i_len = strlen ( subPath_i );
            if ( subPath_i_len == 0 ) {
              XLALPrintInfo ( "%s(): skip empty fallbackpath element\n", __func__ );
              continue;
            }
            if ( subPath_i[0] == '.' ) {
              if ( ! module_path ) {
                XLALPrintInfo ( "%s(): skip relative fallbackpath element '%s': location of liblalsupport.so unknown\n", __func__, subPath_i );
                continue;
              }
              XLAL_CHECK_NULL ( (resolveFname = XLALRealloc ( resolveFname, module_dirname_length + 1 + subPath_i_len + 1 + fname_len + 1 )) != NULL, XLAL_ENOMEM );
              sprintf ( resolveFname, "%s/%s/%s", module_path, subPath_i, fname );
            } else {
              XLAL_CHECK_NULL ( (resolveFname = XLALRealloc ( resolveFname, subPath_i_len + 1 + fname_len + 1 )) != NULL, XLAL_ENOMEM );
              sprintf ( resolveFname, "%s/%s", subPath_i, fname );
            }
            XLALPrintInfo ( "%s(): trying '%s' -> '%s' ...\n", __func__, fname, resolveFname );
#ifdef HAVE_GLOB_H
            glob_t glob_buf;
            int glob_retn = glob( resolveFname, 0, NULL, &glob_buf );
            if ( glob_retn == 0 ) {
              XLAL_CHECK_NULL ( glob_buf.gl_pathc > 0, XLAL_EERR );
              XLAL_CHECK_NULL ( (resolveFname = XLALRealloc ( resolveFname, strlen( glob_buf.gl_pathv[0] ) + 1 )) != NULL, XLAL_ENOMEM );
              strcpy( resolveFname, glob_buf.gl_pathv[0] );
              globfree( &glob_buf );
              XLALFree ( module_path );
              XLALDestroyTokenList ( subPaths );
              XLALPrintInfo ( "%s(): success '%s' -> '%s'\n", __func__, fname, resolveFname );
              return resolveFname;
            } // if found
            char glob_errmsg[256];
            switch (glob_retn) {
            case GLOB_NOSPACE:
              snprintf( glob_errmsg, sizeof( glob_errmsg ), "out of memory" );
              break;
            case GLOB_ABORTED:
              snprintf( glob_errmsg, sizeof( glob_errmsg ), "read error" );
              break;
            case GLOB_NOMATCH:
              snprintf( glob_errmsg, sizeof( glob_errmsg ), "no matches" );
              break;
            default:
              snprintf( glob_errmsg, sizeof( glob_errmsg ), "unknown error %i", glob_retn );
            }
            XLALPrintInfo ( "%s(): trying '%s' -> '%s' ... glob failed (%s)\n", __func__, fname, resolveFname, glob_errmsg );
            globfree( &glob_buf );
#else
            if ( (tmp = LALFopen ( resolveFname, "rb" )) != NULL ) {
              LALFclose ( tmp );
              XLALFree ( module_path );
              XLALDestroyTokenList ( subPaths );
              XLALPrintInfo ( "%s(): success '%s' -> '%s'\n", __func__, fname, resolveFname );
              return resolveFname;
            } // if found
#endif
          }

        XLALFree ( module_path );
        XLALDestroyTokenList ( subPaths );

      }

      XLALFree ( resolveFname );

    } // end: if pure filename given

  // nothing worked? that's ok: return NULL without failure!
  XLALPrintInfo ( "%s(): failed '%s'\n", __func__, fname);
  return NULL;

} // XLALFileResolvePathLong()

/* Simple wrapper to XLAL_FILE_RESOLVE_PATH(), using hardcoded fallbackdir for LAL
   Use of XLAL_FILE_RESOLVE_PATH() directly is preferred; kept for ABI backward compatibility */
char *
XLALFileResolvePath ( const char *fname )
{
  return XLAL_FILE_RESOLVE_PATH( fname );
} // XLALFileResolvePath()
/** \endcond */ // DONT_DOXYGEN

/** Read a complete data-file into memory as a string
 */
char *
XLALFileLoad ( const char *path		//!< [in] input filepath
               )
{
  XLAL_CHECK_NULL ( path != NULL, XLAL_EINVAL );

  // check that this is actually a regular file, rather than sth else (eg a directory)
  size_t blobLen;
  XLAL_CHECK_NULL ( XLALFileIsRegularAndGetSize ( path, &blobLen ) == 1, XLAL_EINVAL, "Path '%s' does not point to a regular file!\n", path );

  char *dataBuffer = NULL;
  size_t dataBufferLen = 0;
  size_t numReadTotal = 0;
  LALFILE *fp;
  XLAL_CHECK_NULL ( (fp = XLALFileOpenRead (path)) != NULL, XLAL_EFUNC );

  // read file in blobs the size of the (possibly compressed) file
  int blobCounter = 0;
  while ( ! XLALFileEOF( fp ) )
    {
      dataBufferLen += blobLen;
      if ( (dataBuffer = XLALRealloc ( dataBuffer, dataBufferLen + 1)) == NULL ) {
        XLALFileClose(fp);
        XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to XLALRealloc(%zu)\n", dataBufferLen+1 );
      }
      size_t numRead = XLALFileRead ( dataBuffer + blobCounter*blobLen, sizeof(char), blobLen, fp );
      if ( xlalErrno != XLAL_SUCCESS ) {
        XLALFree ( dataBuffer );
        XLALFileClose(fp);
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

      numReadTotal += numRead;

      if ( numRead < blobLen ) {
        break;
      }
      blobCounter ++;

    } // while !eof

  XLALFileClose(fp);

  // adjust buffer to final size of data read
  XLAL_CHECK_NULL ( (dataBuffer = XLALRealloc ( dataBuffer, numReadTotal + 1)) != NULL, XLAL_ENOMEM, "Failed to XLALRealloc(%zu)\n", numReadTotal+1 );
  dataBuffer[numReadTotal] = 0;

  return dataBuffer;

} // XLALFileLoad()

int XLALFileIsCompressed( const char *path )
{
  FILE *fp;
  unsigned char magic[2] = { 0, 0 };
  size_t c;
  if ( ! ( fp = LALFopen( path, "rb" ) ) )
    XLAL_ERROR( XLAL_EIO );
  c = fread( magic, sizeof(*magic), XLAL_NUM_ELEM(magic), fp );
  if (c == 0)
    XLAL_ERROR( XLAL_EIO );
  fclose( fp );
  return magic[0] == 0x1f && magic[1] == 0x8b;
}

LALFILE *XLALFileOpenRead( const char *path )
{
  int compression;
  LALFILE *file;
  if ( 0 > (compression = XLALFileIsCompressed(path) ) )
    XLAL_ERROR_NULL( XLAL_EIO );
#ifndef ZLIB_ENABLED
  if ( compression ) {
    XLALPrintError( "XLAL Error - %s: Cannot read compressed file\n", __func__ );
    XLAL_ERROR_NULL( XLAL_EIO );
  }
#endif
  if ( ! ( file = XLALMalloc( sizeof(*file ) ) ) )
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  file->compression = compression;
#ifdef ZLIB_ENABLED
  file->fp = compression ? (void*)gzopen( path, "rb" ) : (void*)LALFopen( path, "rb" );
#else
  file->fp = (void*)LALFopen( path, "rb" );
#endif
  if ( ! file->fp ) {
    XLALFree( file );
    XLAL_ERROR_NULL( XLAL_EIO );
  }
  return file;
}

LALFILE *XLALFileOpenAppend( const char *path, int compression )
{
  LALFILE *file;
  if ( ! ( file = XLALMalloc( sizeof(*file ) ) ) )
    XLAL_ERROR_NULL( XLAL_ENOMEM );
#ifdef ZLIB_ENABLED
  file->fp = compression ? (void*)gzopen( path, "a+" ) : (void*)LALFopen( path, "a+" );
#else
  if ( compression ) {
    XLALPrintWarning( "XLAL Warning - %s: Compression not supported\n", __func__ );
    compression = 0;
  }
  file->fp = (void*)LALFopen( path, "a+" );
#endif
  file->compression = compression;
  if ( ! file->fp ) {
    XLALFree( file );
    XLAL_ERROR_NULL( XLAL_EIO );
  }
  return file;
}

LALFILE *XLALFileOpenWrite( const char *path, int compression )
{
  LALFILE *file;
  if ( ! ( file = XLALMalloc( sizeof(*file ) ) ) )
    XLAL_ERROR_NULL( XLAL_ENOMEM );
#ifdef ZLIB_ENABLED
  file->fp = compression ? (void*)gzopen( path, "wb" ) : (void*)LALFopen( path, "wb" );
#else
  if ( compression ) {
    XLALPrintWarning( "XLAL Warning - %s: Compression not supported\n", __func__ );
    compression = 0;
  }
  file->fp = (void*)LALFopen( path, "wb" );
#endif
  file->compression = compression;
  if ( ! file->fp ) {
    XLALFree( file );
    XLAL_ERROR_NULL( XLAL_EIO );
  }
  return file;
}

LALFILE *XLALFileOpen( const char *path, const char *mode )
{
  int compression;
  char *ext;
  switch ( *mode ) {
  case 'r':
    return XLALFileOpenRead( path );
  case 'w':
    /* check if filename ends in .gz */
    ext = strrchr( path, '.' );
    compression = ext ? ! strcmp( ext, ".gz" ) : 0;
    return XLALFileOpenWrite( path, compression );
  default:
    break; /* fall-out */
  }
  /* error if code gets here */
  XLAL_ERROR_NULL( XLAL_EINVAL );
}

int XLALFileClose( LALFILE *file )
{
  /* this routine acts as a no-op if the file is NULL */
  /* this behavior is different from BSD fclose */
  if ( file ) {
    int c;
    if ( ! file->fp )
      XLAL_ERROR( XLAL_EINVAL );
#ifdef ZLIB_ENABLED
    c = file->compression ? gzclose(((gzFile)file->fp)) : fclose(((FILE*)file->fp));
#else
    c = fclose(((FILE*)file->fp));
#endif
    if ( c == EOF )
      XLAL_ERROR( XLAL_EIO );
    XLALFree( file );
  }
  return 0;
}

size_t XLALFileRead( void *ptr, size_t size, size_t nobj, LALFILE *file )
{
  size_t c;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  c = file->compression ? (size_t)gzread( ((gzFile)file->fp), ptr, size * nobj ) : fread( ptr, size, nobj, ((FILE*)file->fp) );
#else
  c = fread( ptr, size, nobj, ((FILE*)file->fp) );
#endif
  if ( c == (size_t)(-1) || (file->compression == 0 && ferror((FILE*)(file->fp))) )
    XLAL_ERROR( XLAL_EIO );
  return c;
}

size_t XLALFileWrite( const void *ptr, size_t size, size_t nobj, LALFILE *file )
{
  size_t c;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  c = file->compression ? (size_t)gzwrite( ((gzFile)file->fp), ptr, size * nobj ) : fwrite( ptr, size, nobj, ((FILE*)file->fp) );
#else
  c = fwrite( ptr, size, nobj, (FILE*)(file->fp) );
#endif
  if ( c == 0 || (file->compression == 0 && ferror((FILE*)(file->fp))) )
    XLAL_ERROR( XLAL_EIO );
  return c;
}

int XLALFileGetc( LALFILE *file )
{
  int c;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  c = file->compression ? gzgetc(((gzFile)file->fp)) : fgetc(((FILE*)file->fp));
#else
  c = fgetc(((FILE*)file->fp));
#endif
  return c;
}

int XLALFilePutc( int c, LALFILE *file )
{
  int result;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  result = file->compression ? gzputc(((gzFile)file->fp), c) : fputc(c, ((FILE*)file->fp));
#else
  result = fputc(c, (FILE*)(file->fp));
#endif
  if ( result == -1 )
    XLAL_ERROR( XLAL_EIO );
  return result;
}

char * XLALFileGets( char * s, int size, LALFILE *file )
{
  char *c;
  if ( ! file )
    XLAL_ERROR_NULL( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  c = file->compression ? gzgets( ((gzFile)file->fp), s, size ) : fgets( s, size, ((FILE*)file->fp) );
#else
  c = fgets( s, size, ((FILE*)file->fp) );
#endif
  return c;
}

int XLALFilePuts( const char * s, LALFILE *file )
{
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
  if ( 0 > (int)XLALFileWrite( s, sizeof(*s), strlen(s), file ) )
    XLAL_ERROR( XLAL_EFUNC );
  return 0;
}

int XLALFileVPrintf( LALFILE *file, const char *fmt, va_list ap )
{
  char buf[32768];
  int len;
  int c;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
  len = vsnprintf( buf, sizeof(buf), fmt, ap );
  if ( len < 0 )
    XLAL_ERROR( XLAL_EFAILED );
  if ( len >= (int)sizeof(buf) ) { /* buffer is too small */
    char *s;
    s = XLALMalloc( len + 1 );
    if ( !s )
      XLAL_ERROR( XLAL_ENOMEM );
    len = vsnprintf( s, len + 1, fmt, ap );
    c = XLALFilePuts( s, file );
    XLALFree( s );
  } else {
    c = XLALFilePuts( buf, file );
  }
  if ( c < 0 )
    XLAL_ERROR( XLAL_EFUNC );
  return len;
}

int XLALFilePrintf( LALFILE *file, const char *fmt, ... )
{
  int c;
  va_list ap;
  va_start( ap, fmt );
  c = XLALFileVPrintf( file, fmt, ap );
  va_end( ap );
  if ( c < 0 )
    XLAL_ERROR( XLAL_EFUNC );
  return c;
}


int XLALFileFlush( LALFILE *file )
{
  int c;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  c = file->compression ? gzflush(((gzFile)file->fp), Z_FULL_FLUSH) : fflush(((FILE*)file->fp));
#else
  c = fflush(((FILE*)file->fp));
#endif
  if ( c == -1 )
    XLAL_ERROR( XLAL_EIO );
  return c;
}

int XLALFileSeek( LALFILE *file, long offset, int whence )
{
  int c;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  if ( file->compression && whence == SEEK_END ) {
    XLALPrintError( "XLAL Error - %s: SEEK_END not supported with compressed files\n", __func__ );
    XLAL_ERROR( XLAL_EINVAL );
  }
  c = file->compression ? gzseek(((gzFile)file->fp), offset, whence) : fseek(((FILE*)file->fp), offset, whence);
#else
  c = fseek(((FILE*)file->fp), offset, whence);
#endif
  if ( c == -1 )
    XLAL_ERROR( XLAL_EIO );
  return 0;
}

long XLALFileTell( LALFILE *file )
{
  long c;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  c = file->compression ? (long)gztell(((gzFile)file->fp)) : ftell(((FILE*)file->fp));
#else
  c = ftell(((FILE*)file->fp));
#endif
  if ( c == -1 )
    XLAL_ERROR( XLAL_EIO );
  return c;
}

void XLALFileRewind( LALFILE *file )
{
  if ( ! file )
    XLAL_ERROR_VOID( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  file->compression ? (void)gzrewind(((gzFile)file->fp)) : rewind(((FILE*)file->fp));
#else
  rewind(((FILE*)file->fp));
#endif
  return;
}

/** \brief Set buffering for file I/O
 *
 * For a regular file buffering will be set using \c setvbuf. If \c buf is \c NULL then a buffer of \c size will be
 * automatically allocated. The \c mode can be \c _IONBF, \c _IOLBF or \c _IOFBF for no buffering, line buffering
 * or full buffering respectively.
 *
 * For a compressed file the buffering will be set with \c gzbuffer. The \c buf and \c mode inputs are ignored and a
 * buffer of \c size is set.
 */
int XLALFileSetBuffer( LALFILE *file, char *buf, int mode, size_t size )
{
  int c = 0;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  if ( !file->compression ){
    c = setvbuf(((FILE*)file->fp), buf, mode, size);
  }
  else{
#if defined ZLIB_VER_MAJOR && ZLIB_VER_MAJOR >= 1 && ZLIB_VER_MINOR >= 2 && ZLIB_VER_REVISION >= 4
    c = (int)gzbuffer(((gzFile)file->fp), size);
#else
    XLAL_PRINT_WARNING("Ignored buffering: unsupported in zlib version %s", ZLIB_VERSION);
#endif
  }
#else
  c = setvbuf(((FILE*)file->fp), buf, mode, size);
#endif
  if ( c != 0 )
    XLAL_ERROR( XLAL_EIO );
  return c;
}


int XLALFileEOF( LALFILE *file )
{
  int c;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
#ifdef ZLIB_ENABLED
  c = file->compression ? gzeof(((gzFile)file->fp)) : feof((FILE*)(file->fp));
#else
  c = feof((FILE*)(file->fp));
#endif
  return c;
}

/** Check if path points to a 'regular file', rather than a directory or sth else.
 * and also return the size of the given file.
 *
 * This is simply a wrapper to stat(), and S_ISREG()
 * and could be easily generalized along those lines to test for directories etc.
 *
 * return an error if stat() was unavailable
 * return 1 for true, 0 for false, -1 on error
 */
int
XLALFileIsRegularAndGetSize ( const char *path,	//!< [in] path to file
                              size_t *fileLen	//!< [out] size in bytes of file
                              )
{
#ifndef HAVE_STAT
  XLAL_ERROR ( XLAL_EFAILED, "No implementation of function stat() available, which is required here!\n");
#else
  XLAL_CHECK ( path != NULL, XLAL_EINVAL );

  struct stat sb;
  XLAL_CHECK ( stat( path, &sb) != -1, XLAL_EIO, "stat() failed: perror = '%s'\n", strerror(errno) );

  if ( fileLen ) {	// optional output argument
    (*fileLen) = (size_t)sb.st_size;
  }

  return (S_ISREG(sb.st_mode));
#endif

} // XLALFileIsRegularAndGetSize()


/** Check if given file is 'regular' (rather than a directory or sth else)
 *
 * This is a simple wrapper to XLALFileIsRegularAndGetSize().
 *
 * Return 1 for true, 0 for false, -1 on error
 */
int
XLALFileIsRegular ( const char *path	//!< [in] path to file
                    )
{
  return XLALFileIsRegularAndGetSize ( path, NULL );
} // XLALFileIsRegular()


/** Return the size of given file in bytes.
 *
 * This is simply a wrapper to XLALFileIsRegularAndGetSize()
 *
 * Returns (size_t)-1 on error
 */
size_t
XLALFileSize ( const char *path		//!< [in] path to file
               )
{
  size_t size;

  XLAL_CHECK ( XLALFileIsRegularAndGetSize ( path, &size ) != -1, XLAL_EFUNC );

  return size;

} // XLALFileSize()

/**
 * \brief Use gzip to compress a text file
 * This function will use the gzip compression routines in \c zlib to compress a text file. The compressed file will
 * have the extension ".gz" and the original uncompressed file will be removed.
 *
 * \param filename [in] The input text file.
 */
int XLALGzipTextFile( const char *filename ){
  char *memblock = NULL; /* memory block to read data into */
  long pos;
  char *outname = NULL;

  LALFILE *fp = NULL; /* input file pointer */
  LALFILE *fg = NULL; /* output gzipped file pointer */

  /* open file for reading and work out its size */
  if( (fp = XLALFileOpen( filename, "rb" )) == NULL ){
    XLALPrintError ("%s: Unable to open file %s.\n", __func__, filename );
    XLAL_ERROR( XLAL_EIO );
  }

  XLALFileSeek( fp, 0, SEEK_END ); /* go to end of the file */
  pos = XLALFileTell( fp ); /* get the position of the file pointer */
  XLALFileRewind( fp ); /* return the file pointer to the start of the file */

  /* allocate memory to read in data */
  memblock = XLALMalloc( pos + 1);

  if ( !memblock ){
    XLALPrintError ("%s: Unable to allocate memory for reading file.\n", __func__ );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* read in data */
  if ( XLALFileRead( memblock, pos, 1, fp ) == 0 ){
    XLALPrintError ("%s: Unable to read in file.\n", __func__ );
    XLAL_ERROR( XLAL_EIO );
  }
  memblock[pos] = 0; // properly 0-terminate string

  XLALFileClose( fp );

  /* create output name with .gz appended */
  if ( (outname = XLALCalloc(strlen(filename)+4, 1)) == NULL)      /* we need local copy */{
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  strcpy(outname, filename);
  strcat(outname, ".gz"); /* add the extension */

  if ( !outname ){
    XLALPrintError("%s: Unable to create output filename.\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* open output gzip file */
  if ( (fg = XLALFileOpen( outname, "wb" )) == NULL ){
    XLALPrintError ("%s: Unable to open output file %s.\n", __func__, outname );
    XLAL_ERROR( XLAL_EIO );
  }

  if ( XLALFilePuts( memblock, fg ) != 0 ){
    XLALPrintError ("%s: Unable to output gzipped data.\n", __func__ );
    XLAL_ERROR( XLAL_EIO );
  }

  XLALFileClose( fg );

  /* remove original file */
  if ( remove( filename ) == -1 ){
    XLALPrintError ("%s: Unable to remove original text file.\n", __func__ );
    XLAL_ERROR( XLAL_EFAILED );
  }

  return XLAL_SUCCESS;
}


/**
 * \brief Use gzip to uncompress a compressed text file
 *
 * This function will use the gzip compression routines in \c zlib to uncompress a gzipped text file. The compressed
 * file should have the ".gz" extension otherwise it will be rejected. The output will have the same filename with the
 * ".gz" extension removed. The compressed file will be removed.
 *
 * Note: \c gzopen will check the file's "magic number" to see if it is a gzipped file. If it's not a gzipped file it
 * will still open the file for reading as a standard file.
 *
 * \param filename [in] The input gzipped text file.
 */
int XLALGunzipTextFile( const char *filename ){
  CHAR *memblock = NULL; /* memory block to read data into */
  CHAR *outname = NULL;
  CHAR *gzpos = NULL;

  LALFILE *fp = NULL; /* output file pointer */
  LALFILE *fg = NULL; /* input gzipped file pointer */
  size_t n;

  /* create output file name by striping .gz from the input file name */
  if ( (gzpos = strstr(filename, ".gz")) == NULL ){
    XLALPrintError ("%s: File %s does not contain the .gz extension.\n", __func__, filename );
    XLAL_ERROR ( XLAL_EIO );
  }

  n = gzpos-filename;
  outname = XLALMalloc( sizeof(CHAR)*(n+1) );

  if ( ( outname = strncpy(outname, filename, n) ) == NULL ){
    XLALPrintError ("%s: Unable to strip extension from file string.\n", __func__ );
    XLAL_ERROR ( XLAL_EIO );
  }
  outname[n] = 0;

  /* open gzipped file for reading */
  if ( (fg = XLALFileOpen(filename, "rb")) == NULL ){
    XLALPrintError ("%s: Unable to open gzipped file %s.\n", __func__, filename );
    XLAL_ERROR ( XLAL_EIO );
  }

  /* open output text file */
  if ( (fp = XLALFileOpen(outname, "w")) == NULL ){
    XLALPrintError ("%s: Unable to open output file %s.\n", __func__, outname );
    XLAL_ERROR ( XLAL_EIO );
  }

  /* allocate memory for a single char */
  memblock = XLALMalloc( sizeof(CHAR) );

  /* read in gzipped data one byte at a time and output to text file */
  while( ( memblock = XLALFileGets(memblock, 2, fg) ) != NULL ){ XLALFilePrintf(fp, "%s", memblock); }

  /* close files */
  XLALFileClose( fg );
  XLALFileClose( fp );

  /* remove original file */
  if ( remove( filename ) == -1 ){
    XLALPrintError ("%s: Unable to remove original gzipped file.\n", __func__ );
    XLAL_ERROR( XLAL_EFAILED );
  }

  return XLAL_SUCCESS;
}
