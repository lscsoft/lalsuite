#ifndef _FILEIO_H
#define _FILEIO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <lal/LALRCSID.h>

NRCSID( FILEIOH, "$Id$" );

#ifndef LALFopen
#define LALFopen fopen
#endif

#ifndef LALFclose
#define LALFclose fclose
#endif

FILE *
LALOpenDataFile( const char * );

#ifdef __cplusplus
}
#endif

#endif /* _FILEIO_H */
