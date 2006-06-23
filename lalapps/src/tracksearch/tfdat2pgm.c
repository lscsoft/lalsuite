/*
 * Author: Torres C. (Univ of TX at Brownsville)
 */

#include <lalapps.h>
#include "tracksearch.h"
#include "tracksearchAverager.h"


#define TSAPERMS 0666

typedef struct
{
  INT4    argc;
  CHAR**  argv;
}LALInitSearchParams;

/* Code Identifying information */
NRCSID( TRACKSEARCHAVERAGERC, "tracksearchAverager $Id$");
RCSID( "tracksearchAverager $Id$");
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */
#define TRUE     1
#define FALSE    0

/* Usage format string. */
#define USAGE "Still in flux"


/*
 * Begin MAIN
 */
int main (int argc, char *argv[])
{
  CHARVector *file2convert=NULL;
  LALStatus   status=blank_status;
  TSAMap     *map=NULL;

  /*
   *Sleep for Attaching DDD 
   */
  unsigned int doze = 0;
  pid_t myPID;
  myPID = getpid( );
  fprintf( stdout, "pid %d sleeping for %d seconds\n", myPID, doze );
  fflush( stdout );
  sleep( doze );
  fprintf( stdout, "pid %d awake\n", myPID );
  fflush( stdout );

 /* SET LAL DEBUG STUFF */
  set_debug_level("MEMDBG");
  memset(&status, 0, sizeof(status));
  lal_errhandler = LAL_ERR_ABRT;
  /*ls  set_debug_level("ALLDBG");*/


  struct option long_options[] =
    {
      {"filename",            required_argument,  0,    'a'},
      {0,                     0,                  0,      0}
    };
  int  C;
  while (TRUE)
    {
      int option_index=0;
      C = getopt_long_only(argc,
			   argv,
			   "a:b:c:d:e:f:h:i:j:k:l:m:o:p:q:r:s:t:u",
			   long_options, 
			   &option_index);
      /* The end of the arguments is C = -1 */
      if ( C == -1)
	{
	  break;
	}
      switch( C )
	{
	case 0:
	  /* if this option set a flag, do nothing else now */
	  if ( long_options[option_index].flag != 0 )
	    {
	      break;
	    }
	  else
	    {
	      fprintf( stderr, "error parsing option %s with argument %s\n",
		       long_options[option_index].name, optarg );
	      exit( 1 );
	    }
	  break;
	  
	case 'a':
	  /* Setting the GPS start time parameter */
	  {
	    LAL_CALL(LALCHARCreateVector(&status,&file2convert,512),&status);
	    strncpy(file2convert->data,optarg,strlen(optarg)+1);
	  }
	  break;
	};
    };
  /*
   * Open the dat file
   */
  LALappsTSAReadMapFile(&status,
			&map,
			file2convert);
  /*
   * Write the corresponding files
   */
  LALappsTSAWritePGM(&status,
		     map,
		     NULL);
  /*
   * Free RAM
   */
  if (map)
    {
      LALappsTSADestroyMap(&status,&map);
    }
  if (file2convert)
    {
      LAL_CALL(LALCHARDestroyVector(&status,&file2convert),&status);
    }
   LALCheckMemoryLeaks();
  return 0;
}
