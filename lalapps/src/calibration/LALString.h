#include <stddef.h>

char * XLALStringAppend( char *s, const char *append );
char * XLALStringDuplicate( const char *s );
size_t XLALStringCopy( char *dst, const char *src, size_t size );
size_t XLALStringConcatenate( char *dst, const char *src, size_t size );
