#ifndef ERRUTIL_H
#define ERRUTIL_H

#include <lal/LALDatatypes.h>
#include <lal/XLALError.h>

/* global flag to abort on error (rather than just exit) */
extern int abrtflg;

/* global flag to print verbose messages */
extern int vrbflg;

int error( const char *fmt, ... );
int verbose( const char *fmt, ... );

void XLALAbortErrorHandler( const char *func, const char *file, int line,
    int errnum );
void XLALExitErrorHandler( const char *func, const char *file, int line,
    int errnum );

void set_abrt_on_error( void );
void set_exit_on_error( void );

#endif /* ERRUTIL_H */
