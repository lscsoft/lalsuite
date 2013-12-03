 /*
  * Copyright (C) 2004, 2005 Cristina V. Torres
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
  memset(&status, 0, sizeof(status));
  lal_errhandler = LAL_ERR_ABRT;


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
	    file2convert=XLALCreateCHARVector(512);
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
  LALappsTSAWritePGM(
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
      XLALDestroyCHARVector(file2convert);
    }
   LALCheckMemoryLeaks();
  return 0;
}
 
