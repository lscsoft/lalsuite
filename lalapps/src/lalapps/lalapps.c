/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton
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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALString.h>
#include <lalapps.h>

#include <config.h>
#include <LALAppsVCSInfoHeader.h>

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


/**
 * Function that assembles a default VCS info/version string from LAL and LALapps
 * Also checks LAL header<>library version consistency and returns NULL on error.
 *
 * The VCS version string is allocated here and must be freed by caller.
 *
 * \deprecated Use XLALVCSInfoString() instead.
 */
char *
XLALGetVersionString( int level )
{
  XLAL_PRINT_DEPRECATION_WARNING("XLALVCSInfoString");
  char *str = NULL;
  XLAL_CHECK_NULL( ( str = XLALVCSInfoString( lalAppsVCSInfoList, level, "%% " ) ) != NULL, XLAL_EFUNC );
  return str;
} /* XLALGetVersionString() */


/**
 * Simply outputs version information to fp.
 *
 * Returns != XLAL_SUCCESS on error (version-mismatch or writing to fp)
 */
int
XLALOutputVersionString ( FILE *fp, int level )
{
  char *VCSInfoString;

  if (!fp ) {
    XLALPrintError ("%s: invalid NULL input 'fp'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( (VCSInfoString = XLALGetVersionString(level)) == NULL ) {
    XLALPrintError("%s: XLALGetVersionString() failed.\n", __func__);
    XLAL_ERROR ( XLAL_EFUNC );
  }

  if ( fprintf (fp, "%s", VCSInfoString ) < 0 ) {
    XLALPrintError("%s: fprintf failed for given file-pointer 'fp'\n", __func__);
    XLALFree ( VCSInfoString);
    XLAL_ERROR ( XLAL_EIO );
  }

  XLALFree ( VCSInfoString);

  return XLAL_SUCCESS;

} /* XLALOutputVersionString() */
