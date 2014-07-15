/*
*  Copyright (C) 2007 Bernd Machenschalk, John Whelan
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

#include <lal/LALStdlib.h>
#include <stdio.h>
#include <string.h>

extern BOOLEAN optVerbose;

static INT4
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
