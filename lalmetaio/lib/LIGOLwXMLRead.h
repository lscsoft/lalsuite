/*
*  Copyright (C) 2007 Alexander Dietz, Duncan Brown, Jolien Creighton, Kipp Cannon, Lisa M. Goggin, Patrick Brady, Robert Adam Mercer, Stephen Fairhurst, Thomas Cokelaer
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

/*-----------------------------------------------------------------------
 *
 * File Name: LIGOLwXMLRead.h
 *
 * Author: Brown, D. A. and Fairhurst, S.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A. and Fairhurst, S.
 * \file
 * \ingroup lalmetaio_general
 * \brief Provides functions for reading LIGO lightweight XML files to LIGO metadata database tables.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LIGOLwXMLRead.h>
 * \endcode
 *
 */

#ifndef _LIGOLWXMLREAD_H
#define _LIGOLWXMLREAD_H

#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#endif

/* Forward declarations of MetaIO types.  Note that metaio.h is not
 * included by this header file, so the metaio API does not become part of
 * the API exported by lalmetaio.  The MetaioParseEnvironment structure is
 * an opaque type, here, and is why the forward declaration is neeed. */
struct MetaioParseEnvironment;

int
XLALLIGOLwFindColumn(
    struct MetaioParseEnvironment *env,
    const char *name,
    unsigned int type,
    int required
);

long long XLALLIGOLwParseIlwdChar(
    const struct MetaioParseEnvironment *env,
    int column_number,
    const char *ilwd_char_table_name,
    const char *ilwd_char_column_name
);

int
XLALLIGOLwHasTable(
    const char *filename,
    const char *table_name
);

ProcessTable *
XLALProcessTableFromLIGOLw (
    const char *filename
);

ProcessParamsTable *
XLALProcessParamsTableFromLIGOLw (
    const char *filename
);

TimeSlide *
XLALTimeSlideTableFromLIGOLw (
    const char *filename
);

/* these functions need to be lalified, but they are in support... */

SearchSummaryTable *
XLALSearchSummaryTableFromLIGOLw (
    const char *fileName
    );

#ifdef  __cplusplus
}
#endif

#endif /* _LIGOLWXMLREAD_H */
