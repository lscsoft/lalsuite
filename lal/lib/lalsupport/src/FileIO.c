/*
 *  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Reinhard Prix
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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/**
 * \addtogroup FileIO_h
 *
 * ### Obsolete LAL Prototypes ###
 *
 * \code
 * FILE *LALFopen( const char *path, const char *mode );
 * int LALFclose( FILE *stream );
 * \endcode
 *
 * ### Description ###
 *
 * The routines <tt>LALFopen()</tt> and <tt>LALFclose()</tt> are macro defined to be
 * the same as the standard C routines <tt>LALFopen()</tt> and <tt>fclose()</tt>.  These
 * should only be used in test programs.
 *
 * The routine <tt>LALOpenDataFile()</tt> is used to open a data file for reading.
 * This routine is also to be used in test programs only.  Unless the data file
 * is specified with an absolute path (beginning with a <tt>/</tt>), or a specific
 * path (beginning with a <tt>./</tt> or a <tt>../</tt>), the directory
 * that the file is in is obtained from the environment variable
 * \c LAL_DATA_PATH, which must be set at run-time.  (If the environment
 * variable is not set, the default path is <tt>.</tt> --- i.e., the current
 * directory.)
 *
 * \c LAL_DATA_PATH should typically be set to
 * <tt>/usr/local/share/lal</tt>, or wherever LAL data is installed in your system
 * (which may be different if you used a <tt>--prefix</tt> argument when
 * configuring LAL), but when the test suite is run with <tt>make check</tt>, the
 * variable \c LAL_DATA_PATH is set to the current source directory.  If the
 * filename (including the directory path) is too long (more than 256
 * characters), <tt>LALOpenDataFile()</tt> returns \c NULL and sets
 * \c errno to \c ENAMETOOLONG.
 *
 * \c LAL_DATA_PATH can be any colon-delimeted list of directories, which
 * are searched in order (just like the \c PATH environment variable).
 * An extra colon inserts the default data directory
 * (\f$\langle\f$prefix\f$\rangle\f$<tt>/share/lal</tt>) into the search path at that
 * point.  E.g., a leading/trailing colon will look for the default data
 * directory at the start/end of the list of directories.
 *
 * It is strongly recommended that <tt>LALOpenDataFile()</tt> be used when writing test code.
 *
 */

#include "config.h"

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <zlib.h>
#define ZLIB_ENABLED

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>

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

#define STR( x ) #x
#define XSTR( x ) STR( x )
#define INFOMSG( msg, file ) ( ( lalDebugLevel & LALINFO ) ? \
                               LALPrintError( "Info: function LALOpenDataFile, file " __FILE__ ", line " \
                                              XSTR( __LINE__ ) ", $Id$\n\t%s %s\n", msg, file ) : 0 )
#define ERRORMSG( file ) ( ( lalDebugLevel & LALERROR ) ? \
                           LALPrintError( "Info: function LALOpenDataFile, file " __FILE__ ", line " \
                                          XSTR( __LINE__ ) ", $Id$\n\tCould not open data file %s\n", file ) : 0 )

#ifndef LAL_PREFIX
#define LAL_PREFIX "/usr/local"
#endif



FILE *
LALOpenDataFile( const char *fname )
{
  FILE *fp;
  const char *path;
  char *datapath;	/* locally allocated copy of env-var LAL_DATA_PATH */
  const char *p0;       /* pointer to current sub-path of datapath*/
  char *p1;		/* pointer to next sub-path */
  char  fdata[265];
  int   n;

  if ( (fname==NULL) || ( strlen(fname)==0) ) {
    return NULL;
  }

  if ( *fname == '/' ) /* absolute path is given */
  {
    fp = LALFopen( fname, "r" );
    if ( ! fp )
      ERRORMSG( fname );
    else
      INFOMSG( "Opening data file", fname );
    return fp;
  }

  n = strlen( fname );
  if ( *fname == '.' && n > 0 && ( fname[1] == '/' || ( n > 1 && fname[1] == '.'
                                                        && fname[2] == '/' ) ) ) /* specific path is given */
  {
    fp = LALFopen( fname, "r" );
    if ( ! fp )
      ERRORMSG( fname );
    else
      INFOMSG( "Opening data file", fname );
    return fp;
  }

  path = getenv( "LAL_DATA_PATH" );

  if ( ! path || ! strlen( path ) ) /* path is NULL or empty */
  {
    fp = LALFopen( fname, "r" );
    if ( ! fp )
      ERRORMSG( fname );
    else
      INFOMSG( "Opening data file", fname );
    return fp;
  }

  /* scan through all directories in colon-delmited list of directories */
  if ( (datapath = LALCalloc (strlen(path)+1, 1)) == NULL)	/* we need local copy */
  {
    ERRORMSG( fname );
    return NULL;
  }
  strcpy (datapath, path);
  p0 = datapath;
  do {
    p1 = strchr( p0, ':' ); /* look for additional directories */
    if ( p1 ) /* there are more things in the list */
      *p1++ = 0; /* NUL-terminate current directory */
    if ( ! strlen( p0 ) ) /* this directory is empty */
      p0 = LAL_PREFIX "/share/lal"; /* default data directory */

    n = snprintf( fdata, sizeof(fdata), "%s/%s", p0 ? p0 : ".", fname );
    if ( n > (int) sizeof( fdata ) ) /* data file name too long */
    {
      errno = ENAMETOOLONG;
      LALFree (datapath);
      return NULL;
    }

    INFOMSG( "Looking for file", fdata );
    fp = LALFopen( fdata, "r" );
    if ( fp ) /* we've found it! */
    {
      INFOMSG( "Opening data file", fdata );
      LALFree (datapath);
      return fp;
    }

    p0 = p1;
  }
  while ( p0 );

  LALFree (datapath);
  ERRORMSG( fname );
  return NULL;
}


int XLALFileIsCompressed( const char *path )
{
  FILE *fp;
  unsigned char magic[2] = { 0, 0 };
  size_t c;
  if ( ! ( fp = LALFopen( path, "rb" ) ) )
    XLAL_ERROR( XLAL_EIO );
  c = fread( magic, sizeof(*magic), sizeof(magic)/sizeof(*magic), fp );
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
  char buf[LAL_PRINTF_BUFSIZE];
  int len;
  int c;
  if ( ! file )
    XLAL_ERROR( XLAL_EFAULT );
  len = vsnprintf( buf, LAL_PRINTF_BUFSIZE, fmt, ap );
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
  memblock = XLALMalloc( pos );

  if ( !memblock ){
    XLALPrintError ("%s: Unable to allocate memory for reading file.\n", __func__ );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* read in data */
  if ( XLALFileRead( memblock, pos, 1, fp ) == 0 ){
    XLALPrintError ("%s: Unable to read in file.\n", __func__ );
    XLAL_ERROR( XLAL_EIO );
  }

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
