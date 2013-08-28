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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#ifndef _FILEIO_H
#define _FILEIO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdarg.h>
#include <lal/LALStdio.h>

/* maximum string size to print with LAL Printf routines */
#define LAL_PRINTF_BUFSIZE 4096

/**
 * \addtogroup FileIO_h
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
/*@{*/

FILE *LALOpenDataFile( const char* );

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
int XLALFileEOF( LALFILE *file );

int XLALGzipTextFile( const char *path );
int XLALGunzipTextFile( const char *filename );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* _FILEIO_H */
