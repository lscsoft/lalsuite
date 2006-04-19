#ifndef _LALSTRING_H
#define _LALSTRING_H

#include <stddef.h>
#include <lal/LALRCSID.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID( LALSTRINGH, "$Id$" );

char * XLALStringAppend( char *s, const char *append );
char * XLALStringDuplicate( const char *s );
size_t XLALStringCopy( char *dst, const char *src, size_t size );
size_t XLALStringConcatenate( char *dst, const char *src, size_t size );

#ifdef __cplusplus
}
#endif

#endif /* _LALSTRING_H */
