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
#include <stdlib.h>
#include <lalapps.h>
#include <lal/LALStdlib.h>
#include <lal/LALVersion.h>

const char *rcsid = "$Id$";
extern int vrbflg;

int main( void )
{
  LALStatus status = blank_status;
  char msg[16384];

  vrbflg = 1;
  set_debug_level( "LALMSGLVL3" );
  lal_errhandler = LAL_ERR_EXIT;

  /* print version of this program */
  fprintf( stdout, "LALApps Version:     %s\n", LALAPPS_VERSION );
  fprintf( stdout, "CVS Tag:             %s\n", LALAPPS_CVS_TAG );
  fprintf( stdout, "Configure Date:      %s\n", LALAPPS_CONFIGURE_DATE );
  fprintf( stdout, "Configure Arguments: %s\n", LALAPPS_CONFIGURE_ARGS );
  fprintf( stdout, "(RCS %s)\n\n", rcsid );

  LAL_CALL( LALVersion( &status, msg, sizeof( msg ), vrbflg ), &status );
  puts( msg );

  /* check consistency of lal version */
  if ( strcmp( LAL_VERSION, lalVersion ) ||
       strcmp( LAL_CONFIGURE_ARGS, lalConfigureArgs ) ||
       strcmp( LAL_CONFIGURE_DATE, lalConfigureDate ) )
  {
    fputs( "LAL Version Mismatch!\n\n", stderr );
    fputs( "Header Version:      ",     stderr );
    fputs( LAL_VERSION,                 stderr );
    fputs( "\nConfigure Date:      ",   stderr );
    fputs( LAL_CONFIGURE_DATE,          stderr );
    fputs( "\nConfigure Arguments: ",   stderr );
    fputs( LAL_CONFIGURE_ARGS,          stderr );
    fputs( "\n\n",                      stderr );
    fputs( "Library Version:     ",     stderr );
    fputs( lalVersion,                  stderr );
    fputs( "\nConfigure Date:      ",   stderr );
    fputs( lalConfigureDate,            stderr );
    fputs( "\nConfigure Arguments: ",   stderr );
    fputs( lalConfigureArgs,            stderr );
    fputs( "\n",                        stderr );
    return 1;
  }

  LALCheckMemoryLeaks();
  return 0;
}
