#ifndef ERRUTIL_H
#define ERRUTIL_H

#include <lal/LALDatatypes.h>
#include <lal/XLALError.h>

/* global flag to abort on error (rather than just exit) */
extern int abrtflg;

/* global flag to print verbose messages */
extern int vrbflg;

/* print an error message and either exit or abort */
int error( const char *fmt, ... );

/* print an message if verbose messaging is turned on */
int verbose( const char *fmt, ... );

/* XLAL error handler to abort on error */
void XLALAbortErrorHandler( const char *func, const char *file, int line,
    int errnum );

/* XLAL error handler to exit on error */
void XLALExitErrorHandler( const char *func, const char *file, int line,
    int errnum );

/* set handlers to abort on error */
void set_abrt_on_error( void );

/* set handlers to exit on error */
void set_exit_on_error( void );

#endif /* ERRUTIL_H */
