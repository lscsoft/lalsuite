/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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

/**
 * \file
 * \ingroup LALVersion_h
 *
 * \brief Prints the version and configure options of the LAL library being used.
 *
 * ### Usage ###
 *
 * \code
 * LALVersionTest
 * \endcode
 *
 * ### Description ###
 *
 * This program prints the current version of LAL.\@  If the version information
 * in the library differs from the version information in the header file, this
 * program prints the two versions and exits with code 1.  This is useful for
 * determining which version of the LAL library and header files you are linking
 * to.
 *
 * ### Exit codes ###
 *
 * <table>
 * <tr><th>Code</th><th>Explanation</th></tr>
 * <tr><td>   0</td><td>Success, normal exit.</td></tr>
 * <tr><td>   1</td><td>Version info in library disagrees with header file.</td></tr>
 * <tr><td>   2</td><td>Subroutine failed.</td></tr>
 * </table>
 *
 */
/** \cond DONT_DOXYGEN */

#include <stdio.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALVersion.h>


int main( void )
{
  static LALStatus status;
  char msg[16384];
  int verbose = 1;


  if ( strcmp( LAL_VERSION, lalVersion ) ||
       strcmp( LAL_CONFIGURE_ARGS, lalConfigureArgs ) ||
       strcmp( LAL_CONFIGURE_DATE, lalConfigureDate ) )
  {
    fputs( "LAL Version Mismatch!\n\n", stderr );
    fputs( "Header Version ",           stderr );
    fputs( LAL_VERSION,                 stderr );
    fputs( "\nCompiled on ",            stderr );
    fputs( LAL_CONFIGURE_DATE,          stderr );
    fputs( "\nWith arguments ",         stderr );
    fputs( LAL_CONFIGURE_ARGS,          stderr );
    fputs( "\n\n",                      stderr );
    fputs( "Library Version ",          stderr );
    fputs( lalVersion,                  stderr );
    fputs( "\nCompiled on ",            stderr );
    fputs( lalConfigureDate,            stderr );
    fputs( "\nWith arguments ",         stderr );
    fputs( lalConfigureArgs,            stderr );
    fputs( "\n",                        stderr );
    return 1;
  }

  LALVersion( &status, msg, sizeof( msg ), verbose );

  if ( status.statusCode )
  {
    LALStatus *next = &status;
    do
    {
      fputs( next->statusDescription, stderr );
      fputs( "\n", stderr );
      next = next->statusPtr;
    }
    while ( next );
    return 2;
  }

  puts( msg );

  return 0;
}

/** \endcond */
