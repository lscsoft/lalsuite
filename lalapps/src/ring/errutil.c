#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>

#include "lalapps.h"
#include "errutil.h"

RCSID( "$Id$" );


/* global flag to abort on error (rather than just exit) */
int abrtflg;

extern int vrbflg;

int error( const char *fmt, ... )
{
  va_list ap;
  va_start( ap, fmt );
  vfprintf( stderr, fmt, ap );
  va_end( ap );
  if ( abrtflg )
    exit( 1 );
  else
    abort();
  return 1;
}

int verbose( const char *fmt, ... )
{
  if ( vrbflg )
  {
    va_list ap;
    va_start( ap, fmt );
    vfprintf( stderr, fmt, ap );
    va_end( ap );
  }
  return 0;
}

void XLALAbortErrorHandler( const char *func, const char *file, int line,
    int errnum )
{
  XLALPerror( func, file, line, errnum );
  abort();
}

void XLALExitErrorHandler( const char *func, const char *file, int line,
    int errnum )
{
  XLALPerror( func, file, line, errnum );
  exit(1);
}

void set_abrt_on_error( void )
{
  XLALSetErrorHandler( XLALAbortErrorHandler );
  lal_errhandler = LAL_ERR_ABRT;
  abrtflg = 1;
  return;
}

void set_exit_on_error( void )
{
  XLALSetErrorHandler( XLALExitErrorHandler );
  lal_errhandler = LAL_ERR_EXIT;
  abrtflg = 0;
  return;
}
