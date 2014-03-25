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

/*
 *
 * Use LAL's config.h file rather than LALConfig.h which may be from some
 * previous installation.
 *
 */
#include <config.h>

#include <stdio.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdio.h>
#include <lal/LALVersion.h>
#include <lal/LALVCSInfo.h>
#include <lal/LALBuildInfo.h>

const char *const lalVersion = LAL_VERSION;
const int lalVersionMajor = LAL_VERSION_MAJOR;
const int lalVersionMinor = LAL_VERSION_MINOR;
const int lalVersionMicro = LAL_VERSION_MICRO;
const int lalVersionDevel = LAL_VERSION_DEVEL;

/**
 * \ingroup LALVersion_h
 * Routine that returns the version of LAL.
 * This function writes a version message into the string buffer of specified
 * size (and is truncated if the buffer is too small).  Configuration information
 * is also provided if the verbose flag is set.
 *
 */
void
LALVersion(LALStatus * status, CHAR * message, UINT4 size, INT4 verbose)
{
    INT4 nchar;
    INITSTATUS(status);

    ASSERT(message, status, LALVERSIONH_ENULL, LALVERSIONH_MSGENULL);
    ASSERT(size > 0, status, LALVERSIONH_ESIZE, LALVERSIONH_MSGESIZE);

    nchar = verbose ?
        snprintf(message, size,
                 "LAL Version:         %s\n"
                 "Git ID:              %s\n"
                 "Git Tag:             %s\n"
                 "Build Date:          %s\n"
                 "Configure Date:      %s\n"
                 "Configure Arguments: %s\n",
                 lalVersion, lalHeaderVCSInfo.vcsId,
                 lalHeaderVCSInfo.vcsTag, lalBuildDate, lalConfigureDate,
                 lalConfigureArgs) : snprintf(message, size,
                                              "LAL Version: %s\n",
                                              lalVersion);

    if (nchar < 0) {
        ABORT(status, LALVERSIONH_ESPRN, LALVERSIONH_MSGESPRN);
    }
    if (nchar > (INT4) size) {
        ABORT(status, LALVERSIONH_ESHRT, LALVERSIONH_MSGESHRT);
    }

    RETURN(status);
}
