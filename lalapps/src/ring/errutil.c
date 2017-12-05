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

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>

#include "lalapps.h"
#include "errutil.h"

/* global flag to abort on error (rather than just exit) */
int abrtflg;

extern int vrbflg;

/* print an error message and either exit or abort */
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

/* print an verbose message if verbose messaging is turned on */
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

#if 0
/* XLAL error handler to abort on error */
void XLALAbortErrorHandler( const char *func, const char *file, int line,
    int errnum )
{
  XLALPerror( func, file, line, errnum );
  abort();
}

/* XLAL error handler to exit on error */
void XLALExitErrorHandler( const char *func, const char *file, int line,
    int errnum )
{
  XLALPerror( func, file, line, errnum );
  exit(1);
}
#endif

/* set handlers to abort on error */
void set_abrt_on_error( void )
{
  XLALSetErrorHandler( XLALAbortErrorHandler );
  lal_errhandler = LAL_ERR_ABRT;
  abrtflg = 1;
  return;
}

/* set handlers to exit on error */
void set_exit_on_error( void )
{
  XLALSetErrorHandler( XLALExitErrorHandler );
  lal_errhandler = LAL_ERR_EXIT;
  abrtflg = 0;
  return;
}
