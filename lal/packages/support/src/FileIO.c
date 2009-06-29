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

#if 0 /* autodoc block */
<lalVerbatim file="FileIOCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FileIO.c}}

File IO routines for LAL.  These should \emph{not} be used in routines within
the LAL library---only within test programs.

\subsection*{Prototypes}
\input{FileIOCP}
\begin{verbatim}
FILE *LALFopen( const char *path, const char *mode );
int LALFclose( FILE *stream );
\end{verbatim}
\idx{LALOpenDataFile()}
\idx{LALFopen()}
\idx{LALFclose()}

\subsection*{Description}

The routines \verb+LALFopen()+ and \verb+LALFclose()+ are macro defined to be
the same as the standard C routines \verb+fopen()+ and \verb+fclose()+.  These
should only be used in test programs.

The routine \verb+LALOpenDataFile()+ is used to open a data file for reading.
This routine is also to be used in test programs only.  Unless the data file
is specified with an absolute path (beginning with a \verb+/+), or a specific
path (beginning with a \verb+./+ or a \verb+../+), the directory
that the file is in is obtained from the environment variable
\verb+LAL_DATA_PATH+, which must be set at run-time.  (If the environment
variable is not set, the default path is \verb+.+ --- i.e., the current
directory.)

\verb+LAL_DATA_PATH+ should typically be set to
\verb+/usr/local/share/lal+, or wherever LAL data is installed in your system
(which may be different if you used a \verb+--prefix+ argument when
configuring LAL), but when the test suite is run with \verb+make check+, the
variable \verb+LAL_DATA_PATH+ is set to the current source directory.  If the
filename (including the directory path) is too long (more than 256
characters), \verb+LALOpenDataFile()+ returns \verb+NULL+ and sets
\verb+errno+ to \verb+ENAMETOOLONG+.

\verb+LAL_DATA_PATH+ can be any colon-delimeted list of directories, which
are searched in order (just like the \verb+PATH+ environment variable).
An extra colon inserts the default data directory
($\langle$prefix$\rangle$\verb+/share/lal+) into the search path at that
point.  E.g., a leading/trailing colon will look for the default data
directory at the start/end of the list of directories.

It is strongly recommended that
\verb+LALOpenDataFile()+ be used when writing test code.

\vfill{\footnotesize\input{FileIOCV}}

</lalLaTeX>
#endif /* autodoc block */

#include "config.h"

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_ZLIB_H
/* can't actually include zlib.h since older versions have broken prototypes */
/* #include <zlib.h> */
#define Z_FULL_FLUSH    3
typedef void * gzFile;
extern gzFile gzopen (const char *path, const char *mode);
extern gzFile gzdopen (int fd, const char *mode);
extern int gzsetparams (gzFile file, int level, int strategy);
extern int gzread (gzFile file, void * buf, unsigned len);
extern int gzwrite (gzFile file, const void * buf, unsigned len);
extern int gzprintf (gzFile file, const char *format, ...);
extern int gzputs (gzFile file, const char *s);
extern char * gzgets (gzFile file, char *buf, int len);
extern int gzputc (gzFile file, int c);
extern int gzgetc (gzFile file);
extern int gzungetc (int c, gzFile file);
extern int gzflush (gzFile file, int flush);
extern long gzseek (gzFile file, long offset, int whence);
extern int gzrewind (gzFile file);
extern long gztell (gzFile file);
extern int gzeof (gzFile file);
extern int gzdirect (gzFile file);
extern int gzclose (gzFile file);
extern const char * gzerror (gzFile file, int *errnum);
extern void gzclearerr (gzFile file);
#define ZLIB_ENABLED
#endif

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>

#ifdef lalstdin
#undef lalstdin
#endif
#ifdef lalstdout
#undef lalstdout
#endif
#ifdef lalstderr
#undef lalstderr
#endif
LALFILE * lalstdin( void )
{
	static LALFILE _lalstdin;
	_lalstdin.fp = (void*)stdin;
	return &_lalstdin;
}
LALFILE * lalstdout( void )
{
	static LALFILE _lalstdout;
	_lalstdout.fp = (void*)stdout;
	return &_lalstdout;
}
LALFILE * lalstderr( void )
{
	static LALFILE _lalstderr;
	_lalstderr.fp = (void*)stderr;
	return &_lalstderr;
}

/* use boinc_fopen() instead of fopen() for Einstein@Home/BOINC */
#ifdef EAH_BOINC
extern FILE* boinc_fopen(const char* path, const char* mode);
#define fopen boinc_fopen
#endif

NRCSID (FILEIOC,"$Id$");

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


/* <lalVerbatim file="FileIOCP"> */
FILE *
LALOpenDataFile( const char *fname )
{ /* </lalVerbatim> */
  FILE *fp;
  const char *path;
  char *datapath;	/* locally allocated copy of env-var LAL_DATA_PATH */
  const char *p0; 	/* pointer to current sub-path of datapath*/
  char *p1;		/* pointer to next sub-path */
  char  fdata[265];
  int   n;

  if ( ! fname )
    return NULL;

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
	static const char *func = "XLALFileIsCompressed";
	FILE *fp;
	unsigned char magic[2] = { 0, 0 };
	if ( ! ( fp = fopen( path, "rb" ) ) )
		XLAL_ERROR( func, XLAL_EIO );
	fread( magic, sizeof(*magic), sizeof(magic)/sizeof(*magic), fp );
	fclose( fp );
	return magic[0] == 0x1f && magic[1] == 0x8b;
}

LALFILE * XLALFileOpenRead( const char *path )
{
	static const char *func = "XLALFileOpenRead";
	int compression;
	LALFILE *file;
	if ( 0 > (compression = XLALFileIsCompressed(path) ) )
		XLAL_ERROR_NULL( func, XLAL_EIO );
#	ifndef ZLIB_ENABLED
	if ( compression ) {
		XLALPrintError( "XLAL Error - %s: Cannot read compressed file\n", func );
		XLAL_ERROR_NULL( func, XLAL_EIO );
	}
#	endif
	if ( ! ( file = XLALMalloc( sizeof(*file ) ) ) )
		XLAL_ERROR_NULL( func, XLAL_ENOMEM );
	file->compression = compression;
#	ifdef ZLIB_ENABLED
	file->fp = compression ? gzopen( path, "rb" ) : fopen( path, "rb" );
#	else
	file->fp = fopen( path, "rb" );
#	endif
	if ( ! file->fp ) {
		XLALFree( file );
		XLAL_ERROR_NULL( func, XLAL_EIO );
	}
	return file;
}

LALFILE * XLALFileOpenAppend( const char *path, int compression )
{
	static const char *func = "XLALFileOpenAppend";
	LALFILE *file;
	if ( ! ( file = XLALMalloc( sizeof(*file ) ) ) )
		XLAL_ERROR_NULL( func, XLAL_ENOMEM );
#	ifdef ZLIB_ENABLED
	file->fp = compression ? gzopen( path, "a+" ) : fopen( path, "a+" );
#	else
	if ( compression ) {
		XLALPrintWarning( "XLAL Warning - %s: Compression not supported\n", func );
		compression = 0;
	}
	file->fp = fopen( path, "a+" );
#	endif
	file->compression = compression;
	if ( ! file->fp ) {
		XLALFree( file );
		XLAL_ERROR_NULL( func, XLAL_EIO );
	}
	return file;
}

LALFILE * XLALFileOpenWrite( const char *path, int compression )
{
	static const char *func = "XLALFileOpenWrite";
	LALFILE *file;
	if ( ! ( file = XLALMalloc( sizeof(*file ) ) ) )
		XLAL_ERROR_NULL( func, XLAL_ENOMEM );
#	ifdef ZLIB_ENABLED
	file->fp = compression ? gzopen( path, "wb" ) : fopen( path, "wb" );
#	else
	if ( compression ) {
		XLALPrintWarning( "XLAL Warning - %s: Compression not supported\n", func );
		compression = 0;
	}
	file->fp = fopen( path, "wb" );
#	endif
	file->compression = compression;
	if ( ! file->fp ) {
		XLALFree( file );
		XLAL_ERROR_NULL( func, XLAL_EIO );
	}
	return file;
}

LALFILE * XLALFileOpen( const char *path, const char *mode )
{
	static const char *func = "XLALFileOpen";
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
	XLAL_ERROR_NULL( func, XLAL_EINVAL );
}

int XLALFileClose( LALFILE * file )
{
	static const char *func = "XLALFileClose";
	/* this routine acts as a no-op if the file is NULL */
	/* this behavior is different from BSD fclose */
	if ( file ) {
		int c;
		if ( ! file->fp )
			XLAL_ERROR( func, XLAL_EINVAL );
#		ifdef ZLIB_ENABLED
		c = file->compression ? gzclose(file->fp) : fclose(file->fp);
#		else
		c = fclose(file->fp);
#		endif
		if ( c == EOF )
			XLAL_ERROR( func, XLAL_EIO );
		XLALFree( file );
	}
	return 0;
}

size_t XLALFileRead( void *ptr, size_t size, size_t nobj, LALFILE *file )
{
	static const char *func = "XLALFileRead";
	size_t c;
	if ( ! file )
		XLAL_ERROR( func, XLAL_EFAULT );
#	ifdef ZLIB_ENABLED
	c = file->compression ? (size_t)gzread( file->fp, ptr, size * nobj ) : fread( ptr, size, nobj, file->fp );
#	else
	c = fread( ptr, size, nobj, file->fp );
#	endif
	if ( c == (size_t)(-1) || (file->compression == 0 && ferror((FILE*)(file->fp))) )
		XLAL_ERROR( func, XLAL_EIO );
	return c;
}

size_t XLALFileWrite( const void *ptr, size_t size, size_t nobj, LALFILE *file )
{
	static const char *func = "XLALFileRead";
	size_t c;
	if ( ! file )
		XLAL_ERROR( func, XLAL_EFAULT );
#	ifdef ZLIB_ENABLED
	c = file->compression ? (size_t)gzwrite( file->fp, ptr, size * nobj ) : fwrite( ptr, size, nobj, file->fp );
#	else
	c = fwrite( ptr, size, nobj, (FILE*)(file->fp) );
#	endif
	if ( c == 0 || (file->compression == 0 && ferror((FILE*)(file->fp))) )
		XLAL_ERROR( func, XLAL_EIO );
	return c;
}

int XLALFileGetc( LALFILE *file )
{
	static const char *func = "XLALFileGetc";
	int c;
	if ( ! file )
		XLAL_ERROR( func, XLAL_EFAULT );
#	ifdef ZLIB_ENABLED
	c = file->compression ? gzgetc(file->fp) : fgetc(file->fp);
#	else
	c = fgetc(file->fp);
#	endif
	return c;
}

int XLALFilePutc( int c, LALFILE *file )
{
	static const char *func = "XLALFilePutc";
	int result;
	if ( ! file )
		XLAL_ERROR( func, XLAL_EFAULT );
#	ifdef ZLIB_ENABLED
	result = file->compression ? gzputc(file->fp, c) : fputc(c, file->fp);
#	else
	result = fputc(c, (FILE*)(file->fp));
#	endif
	if ( result == -1 )
		XLAL_ERROR( func, XLAL_EIO );
	return result;
}

char * XLALFileGets( char * s, int size, LALFILE *file )
{
	static const char *func = "XLALFileGets";
	char *c;
	if ( ! file )
		XLAL_ERROR_NULL( func, XLAL_EFAULT );
#	ifdef ZLIB_ENABLED
	c = file->compression ? gzgets( file->fp, s, size ) : fgets( s, size, file->fp );
#	else
	c = fgets( s, size, file->fp );
#	endif
	return c;
}

int XLALFilePuts( const char * s, LALFILE *file )
{
	static const char *func = "XLALFilePuts";
	if ( ! file )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( 0 > (int)XLALFileWrite( s, sizeof(*s), strlen(s), file ) )
		XLAL_ERROR( func, XLAL_EFUNC );
	return 0;
}

int XLALFileVPrintf( LALFILE *file, const char *fmt, va_list ap )
{
	static const char *func = "XLALFileVPrintf";
	char buf[LAL_PRINTF_BUFSIZE];
	int len;
	int c;
	if ( ! file )
		XLAL_ERROR( func, XLAL_EFAULT );
	len = vsnprintf( buf, LAL_PRINTF_BUFSIZE, fmt, ap );
	if ( len < 0 )
		XLAL_ERROR( func, XLAL_EFAILED );
	if ( len >= (int)sizeof(buf) ) { /* buffer is too small */
		char *s;
		s = XLALMalloc( len + 1 );
		if ( !s )
			XLAL_ERROR( func, XLAL_ENOMEM );
		len = vsnprintf( s, len + 1, fmt, ap );
		c = XLALFilePuts( s, file );
		XLALFree( s );
	} else {
		c = XLALFilePuts( buf, file );
	}
	if ( c < 0 )
		XLAL_ERROR( func, XLAL_EFUNC );
	return len;
}

int XLALFilePrintf( LALFILE *file, const char *fmt, ... )
{
	static const char *func = "XLALFilePrintf";
	int c;
	va_list ap;
	va_start( ap, fmt );
	c = XLALFileVPrintf( file, fmt, ap );
	va_end( ap );
	if ( c < 0 )
		XLAL_ERROR( func, XLAL_EFUNC );
	return c;
}


int XLALFileFlush( LALFILE *file )
{
	static const char *func = "XLALFileFlush";
	int c;
	if ( ! file )
		XLAL_ERROR( func, XLAL_EFAULT );
#	ifdef ZLIB_ENABLED
	c = file->compression ? gzflush(file->fp, Z_FULL_FLUSH) : fflush(file->fp);
#	else
	c = fflush(file->fp);
#	endif
	if ( c == -1 )
		XLAL_ERROR( func, XLAL_EIO );
	return c;
}

int XLALFileSeek( LALFILE *file, long offset, int whence )
{
	static const char *func = "XLALFileSeek";
	int c;
	if ( ! file )
		XLAL_ERROR( func, XLAL_EFAULT );
#	ifdef ZLIB_ENABLED
	if ( file->compression && whence == SEEK_END ) {
		XLALPrintError( "XLAL Error - %s: SEEK_END not supported with compressed files\n", func );
		XLAL_ERROR( func, XLAL_EINVAL );
	}
	c = file->compression ? gzseek(file->fp, offset, whence) : fseek(file->fp, offset, whence);
#	else
	c = fseek(file->fp, offset, whence);
#	endif
	if ( c == -1 )
		XLAL_ERROR( func, XLAL_EIO );
	return 0;
}

long XLALFileTell( LALFILE *file )
{
	static const char *func = "XLALFileTell";
	long c;
	if ( ! file )
		XLAL_ERROR( func, XLAL_EFAULT );
#	ifdef ZLIB_ENABLED
	c = file->compression ? (long)gztell(file->fp) : ftell(file->fp);
#	else
	c = ftell(file->fp);
#	endif
	if ( c == -1 )
		XLAL_ERROR( func, XLAL_EIO );
	return 0;
}

void XLALFileRewind( LALFILE *file )
{
	static const char *func = "XLALFileRewind";
	if ( ! file )
		XLAL_ERROR_VOID( func, XLAL_EFAULT );
#	ifdef ZLIB_ENABLED
	file->compression ? (void)gzrewind(file->fp) : rewind(file->fp);
#	else
	rewind(file->fp);
#	endif
	return;
}

int XLALFileEOF( LALFILE *file )
{
	static const char *func = "XLALFileEOF";
	int c;
	if ( ! file )
		XLAL_ERROR( func, XLAL_EFAULT );
#	ifdef ZLIB_ENABLED
	c = file->compression ? gzeof(file->fp) : feof((FILE*)(file->fp));
#	else
	c = feof((FILE*)(file->fp));
#	endif
	return c;
}
