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

#ifndef LALAPPS_H_
#define LALAPPS_H_

#include <config.h>
#include <stdio.h>
#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

NRCSID( LALAPPSH, "$Id$" );

extern const LALStatus blank_status;

typedef int ( *lal_errhandler_t )(
    LALStatus  *,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    );

#define LAL_ERR_DFLT LAL_ERR_ABRT
extern lal_errhandler_t lal_errhandler;

extern int LAL_ERR_EXIT(
    LALStatus  *,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    );
extern int LAL_ERR_ABRT(
    LALStatus  *,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    );
extern int LAL_ERR_RTRN(
    LALStatus  *,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    );
extern int clear_status( LALStatus * );
extern int set_debug_level( const char *s );

extern char *XLALGetVersionString( int level );
extern int XLALOutputVersionString ( FILE *fp, int level );

#define LAL_CALL( function, statusptr ) \
  ((function),lal_errhandler(statusptr,#function,__FILE__,__LINE__,rcsid))

#define PRINT_VERSION( program ) \
  fprintf(stderr,PACKAGE " %s version " VERSION "\n%s\n",program,rcsid)

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* LALAPPS_H_ */
