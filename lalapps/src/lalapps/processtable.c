#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <pwd.h>
#include <unistd.h>
#include <time.h>
#include <lalapps.h>

#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOMetadataTables.h>

#define CVS_REV_STR "$Revision "
#define CVS_SOURCE_STR "$Source "
#define CVS_DATE_STR "$Date "
#define CVS_DATE_FMT "\%Y/\%m/\%d \%T"
#define CVS_DELIM " $"

/* int gethostname(char *name, size_t len); */
char *strptime(const char *s, const char  *format,  struct tm *tm);

NRCSID( PROCESSTABLEC, "$Id$" );

void
populate_process_table (
    LALStatus           *status,
    ProcessTable        *ptable,
    CHAR                *program_name,
    CHAR                *cvs_revision,
    CHAR                *cvs_source,
    CHAR                *cvs_date
    )
{
  const char source_str[] = CVS_SOURCE_STR;
  const char rev_str[] = CVS_REV_STR;
  const char date_str[] = CVS_DATE_STR;
  const char cvs_date_format[] = CVS_DATE_FMT;
  char date_string[256];
  size_t cvsstrlen;
  uid_t userid;
  struct passwd *pwent;
  char *cvsstrstart, *cvsstrend;
  LALDate laldate;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;

  INITSTATUS( status, "populate_process_table", PROCESSTABLEC );
  ATTATCHSTATUSPTR( status );
  ASSERT( ptable, status, 1, "Process table pointer is null" );

  /* program name entry */
  LALSnprintf( ptable->program, LIGOMETA_PROGRAM_MAX, program_name );

  /* cvs version */
  memset( ptable->version, 0, LIGOMETA_VERSION_MAX );
  cvsstrstart = cvs_revision + (size_t) strlen(rev_str) + 1;
  cvsstrend = strstr( cvs_revision, CVS_DELIM );
  cvsstrlen = cvsstrend - cvsstrstart;
  memcpy( ptable->version, cvsstrstart, 
      cvsstrlen < LIGOMETA_VERSION_MAX - 1 ? 
      cvsstrlen : LIGOMETA_VERSION_MAX - 1 );

  /* cvs repository */
  memset( ptable->cvs_repository, 0, LIGOMETA_CVS_REPOSITORY_MAX );
  cvsstrstart = cvs_source + (size_t) strlen(source_str) + 1;
  cvsstrend = strstr( cvs_source, ",v" CVS_DELIM );
  cvsstrlen = cvsstrend - cvsstrstart;
  memcpy( ptable->cvs_repository, cvsstrstart, 
      cvsstrlen < LIGOMETA_CVS_REPOSITORY_MAX - 4 ?
      cvsstrlen : LIGOMETA_CVS_REPOSITORY_MAX - 4 );
  sprintf( ptable->cvs_repository + cvsstrlen, "\\,v" );

  /* cvs check in time */
  memset( date_string, 0, 256 );
  cvsstrstart = cvs_date + (size_t) strlen(date_str) + 1;
  cvsstrend = strstr( cvs_source, " $" );
  cvsstrlen = cvsstrend - cvsstrstart;
  memcpy( date_string, cvsstrstart, 
      cvsstrlen < 255 ? cvsstrlen : 255 );
  if ( ! strptime( date_string, cvs_date_format, &(laldate.unixDate) ) )
  {
    fprintf( stderr, "could not determine cvs checkin date\n" );
    exit( 1 );
  }
  laldate.residualNanoSeconds = 0;
  LALUTCtoGPS( status->statusPtr, &(ptable->cvs_entry_time), &laldate, 
      &accuracy );
  CHECKSTATUSPTR( status );

  /* put blank space in comment */
  LALSnprintf( ptable->comment, LIGOMETA_COMMENT_MAX, " " );

  /* online flag and domain */
  ptable->is_online = 0;
  LALSnprintf( ptable->domain, LIGOMETA_DOMAIN_MAX , "lalapps" );
  
  /* process id, username and host */
  ptable->unix_procid = getpid();
  if ( gethostname( ptable->node, LIGOMETA_NODE_MAX ) < 0 )
  {
    perror( "could not determine host name" );
    exit( 1 );
  }
  userid = geteuid();
  if ( ! (pwent = getpwuid( userid )) )
  {
    LALSnprintf( ptable->username, LIGOMETA_USERNAME_MAX - 1, "%d", userid );
  }
  strncpy( ptable->username, pwent->pw_name, LIGOMETA_USERNAME_MAX - 1);

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
