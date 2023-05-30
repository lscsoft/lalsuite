/*
 *  Copyright (C) 2007 Jolien Creighton
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

#ifndef _FILEIO_H
#define _FILEIO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdarg.h>
#include <lal/LALStddef.h>
#include <lal/LALStdio.h>

/**
 * \defgroup FileIO_h Header FileIO.h
 * \ingroup lal_support
 *
 * \brief Provides standard LAL support IO functions.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LALStdio.h>
 * #include <lal/FileIO.h>
 * \endcode
 *
 * Only use \ref FileIO.h in test code that links to the \c lalsupport library.
 */
/** @{ */

typedef struct tagLALFILE LALFILE;

LALFILE *lalstdin(void);
LALFILE *lalstdout(void);
LALFILE *lalstderr(void);
#define LALSTDIN  (lalstdin())
#define LALSTDOUT (lalstdout())
#define LALSTDERR (lalstderr())

int XLALFileIsCompressed( const char *path );
LALFILE *XLALFileOpenRead( const char *path );
LALFILE *XLALFileOpenWrite( const char *path, int compression );
LALFILE *XLALFileOpenAppend( const char *path, int compression );
LALFILE *XLALFileOpen( const char *path, const char *mode );
int XLALFileClose( LALFILE *file );
size_t XLALFileRead( void *ptr, size_t size, size_t nobj, LALFILE *file );
size_t XLALFileWrite( const void *ptr, size_t size, size_t nobj, LALFILE *file );
int XLALFileGetc( LALFILE *file );
int XLALFilePutc( int c, LALFILE *file );
char *XLALFileGets( char *s, int size, LALFILE *file );
int XLALFilePuts( const char *s, LALFILE *file );
#ifndef SWIG /* exclude from SWIG interface */
int XLALFileVPrintf( LALFILE *file, const char *fmt, va_list ap );
#endif /* SWIG */
int XLALFilePrintf( LALFILE *file, const char *fmt, ... );
int XLALFileFlush( LALFILE *file );
int XLALFileSeek( LALFILE *file, long offset, int whence );
long XLALFileTell( LALFILE *file );
void XLALFileRewind( LALFILE *file );
int XLALFileSetBuffer( LALFILE *file, char *buf, int mode, size_t size );
int XLALFileEOF( LALFILE *file );

int XLALFileIsRegularAndGetSize ( const char *path, size_t *fileLen );
int XLALFileIsRegular ( const char *path );
size_t XLALFileSize ( const char *path );

/** \cond DONT_DOXYGEN */
char *XLALFileResolvePathLong ( const char *fname, const char *fallbackpath );
char *XLALFileResolvePath ( const char *fname );
/** \endcond */ // DONT_DOXYGEN

/** 'Resolve' a given filename 'fname', returning a file path where the
 *  file can successfully be opened by fopen() using mode='rb'.
 *
 * Return: successful file-path or NULL if failed.
 *
 * Resolving uses the following algorithm: if 'fname' contains a
 * i) (relative or absolute) path: only tries to open that path directly
 * ii) pure filename: try
 *     1) local dir, then
 *     2) search $LAL_DATA_PATH, then
 *     3) search 'fallbackpath'; relative paths are considered relative to location of liblalsupport.so (if known)
 *     return first successful hit
 *
 * Note: it is not an error if the given 'fname' cannot be resolved,
 * this will simply return NULL but xlalErrno will not be set in that case.
 *
 * Note2: successfully resolving an 'fname' doesn't guarantee that the path points
 * to a file, as directories can also be opened in 'rb'.
 *
 * Note3: the returned string is allocated here and must be XLALFree'ed by the caller.
 *
 * Note4: this function should only be used from within LALSuite, as it relies
 * on macros defined by the LALSuite build system. For the public API, functions
 * should be defined for specific types of data files, which then call
 * XLAL_FILE_RESOLVE_PATH() internally.
 */
#define XLAL_FILE_RESOLVE_PATH( fname ) XLALFileResolvePathLong ( fname, LAL_FALLBACK_DATA_PATH )

char *XLALFileLoad ( const char *path );

int XLALGzipTextFile( const char *path );
int XLALGunzipTextFile( const char *filename );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* _FILEIO_H */
