#include "CheckStatus.h"
#include <stdio.h>
#include <string.h>
extern BOOLEAN optVerbose;

INT4
CheckStatus(LALStatus *status, const INT4 code, const CHAR *message,
	    const INT4 exitcode, const CHAR *error)
{
  
  if (optVerbose)
  {
    REPORTSTATUS (status);
  }
  if (status->statusCode!= code) 
  {
    if (code) printf ("  FAIL: did not recognize \"%s\"\n", message);
    if (optVerbose) printf("Exiting with error: %s\n", error);
    return(exitcode);
  }
  else if (code && strcmp(message, status->statusDescription)) 
  {
    printf("  FAIL: incorrect error message \"%s\" not \"%s\"\n",
	   status->statusDescription, message);
    if (optVerbose) printf("Exiting with error: %s\n", error);
    return(exitcode);
  }
  return 0;
}
