#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lalapps.h>
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>

#define FAILMSG( stat, func, file, line, id )                                  \
  do {                                                                         \
    if ( lalDebugLevel & LALERROR )                                            \
    {                                                                          \
      LALPrintError( "Error[0]: file %s, line %d, %s\n"                        \
          "\tLAL_CALL: Function call `%s' failed.\n", file, line, id, func );  \
    }                                                                          \
    if ( vrbflg )                                                              \
    {                                                                          \
      fprintf(stderr,"Level 0: %s\n\tFunction call `%s' failed.\n"             \
          "\tfile %s, line %d\n", id, func, file, line );                      \
      REPORTSTATUS( stat );                                                    \
    }                                                                          \
  } while( 0 )

const LALStatus blank_status;
int lalDebugLevel = 0;
int vrbflg = 0;

lal_errhandler_t lal_errhandler = LAL_ERR_DFLT;

int LAL_ERR_EXIT(
    LALStatus  *stat,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    )
{
  if ( stat->statusCode )
  {
    FAILMSG( stat, func, file, line, id );
    exit( 1 );
  }
  return stat->statusCode;
}

int LAL_ERR_ABRT(
    LALStatus  *stat,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    )
{
  if ( stat->statusCode )
  {
    FAILMSG( stat, func, file, line, id );
    abort();
  }
  return 0;
}

int LAL_ERR_RTRN(
    LALStatus  *stat,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    )
{
  if ( stat->statusCode )
  {
    FAILMSG( stat, func, file, line, id );
  }
  return stat->statusCode;
}

int clear_status( LALStatus *stat )
{
  if ( ! stat )
    return 1;
  while ( stat->statusPtr )
  {
    LALStatus *next = stat->statusPtr->statusPtr;
    LALFree( stat->statusPtr );
    stat->statusPtr = next;
  }
  memset( stat, 0, sizeof( *stat ) );
  return 0;
}

int set_debug_level( const char *s )
{
  unsigned level = 0;
  if ( ! s )
  {
    if ( ! ( s = getenv( "LAL_DEBUG_LEVEL" ) ) )
      return lalDebugLevel = 0;
  }

  /* skip whitespace */
  while ( isspace( *s ) )
    ++s;
  
  /* a value is set */
  if ( isdigit( *s ) )
    return lalDebugLevel = atoi( s );

  /* construct the debug level */
  if ( strstr( s, "NDEBUG" ) )
    level |= LALNDEBUG;
  if ( strstr( s, "ERROR" ) )
    level |= LALERROR;
  if ( strstr( s, "WARNING" ) )
    level |= LALWARNING;
  if ( strstr( s, "INFO" ) )
    level |= LALINFO;
  if ( strstr( s, "TRACE" ) )
    level |= LALTRACE;
  if ( strstr( s, "MEMINFO" ) )
    level |= LALMEMINFO;
  if ( strstr( s, "MEMDBG" ) )
    level |= LALMEMDBG;
  if ( strstr( s, "MSGLVL1" ) )
    level |= LALMSGLVL1;
  if ( strstr( s, "MSGLVL2" ) )
    level |= LALMSGLVL2;
  if ( strstr( s, "MSGLVL3" ) )
    level |= LALMSGLVL3;
  if ( strstr( s, "MEMTRACE" ) )
    level |= LALMEMTRACE;
  if ( strstr( s, "ALLDBG" ) )
    level |= LALALLDBG;

  return lalDebugLevel = level;
}
