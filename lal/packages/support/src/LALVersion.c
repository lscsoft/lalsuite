/*
 *
 * Use LAL's config.h file rather than LALConfig.h which may be from some
 * previous installation.
 *
 */
#include <config.h>

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALVersion.h>

NRCSID( LALVERSIONC, "$Id$" );

#undef LALsnprintf

const char *lalVersion       = LAL_VERSION;
const int   lalVersionMajor  = LAL_VERSION_MAJOR;
const int   lalVersionMinor  = LAL_VERSION_MINOR;
const char *lalConfigureArgs = LAL_CONFIGURE_ARGS;
const char *lalConfigureDate = LAL_CONFIGURE_DATE;

#if defined(NDEBUG) || defined(LAL_NDEBUG)
const int lalNoDebug = 1;
#else
const int lalNoDebug = 0;
#endif


int
LALsnprintf( char *str, size_t size, const char *fmt, ... )
{
  va_list ap;
  int n;
  va_start( ap, fmt );
  n = LALvsnprintf( str, size, fmt, ap );
  va_end( ap );
  return n;
}

void
LALVersion( LALStatus *status, CHAR *message, UINT4 size, INT4 verbose )
{
  INT4 nchar;
  INITSTATUS( status, "LALVersion", LALVERSIONC );

  ASSERT( message,  status, LALVERSIONH_ENULL, LALVERSIONH_MSGENULL );
  ASSERT( size > 0, status, LALVERSIONH_ESIZE, LALVERSIONH_MSGESIZE );

  nchar = verbose ?
    LALsnprintf( message, size, "This is LAL Version %s\nCompiled on %s\n"
        "With arguments %s\n(RCS %s)\n",
        lalVersion, lalConfigureDate, lalConfigureArgs, LALVERSIONC ) :
    LALsnprintf( message, size, "This is LAL Version %s\n", lalVersion ) ;

  if ( nchar < 0 )
  {
    ABORT( status, LALVERSIONH_ESPRN, LALVERSIONH_MSGESPRN );
  }
  if ( nchar > (INT4) size )
  {
    ABORT( status, LALVERSIONH_ESHRT, LALVERSIONH_MSGESHRT );
  }

  RETURN( status );
}
