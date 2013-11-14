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
 * \ingroup lalmetaio
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

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**\name Error Codes */ /*@{*/
#define LIGOLWXMLREADH_ENULL 1
#define LIGOLWXMLREADH_ENNUL 2
#define LIGOLWXMLREADH_EALOC 3
#define LIGOLWXMLREADH_EUTAB 4
#define LIGOLWXMLREADH_ENCOL 5
#define LIGOLWXMLREADH_ENTAB 6
#define LIGOLWXMLREADH_EPARS 7
#define LIGOLWXMLREADH_EMTAB 8
#define LIGOLWXMLREADH_EENDT 9
#define LIGOLWXMLREADH_ETMSM 10
#define LIGOLWXMLREADH_ETNOP 11

#define LIGOLWXMLREADH_MSGENULL "Null pointer"
#define LIGOLWXMLREADH_MSGENNUL "Non-null pointer"
#define LIGOLWXMLREADH_MSGEALOC "Memory allocation error"
#define LIGOLWXMLREADH_MSGEUTAB "Unknown metadata table type"
#define LIGOLWXMLREADH_MSGENCOL "Unable to find table column"
#define LIGOLWXMLREADH_MSGENTAB "Requested table not found in file"
#define LIGOLWXMLREADH_MSGEPARS "Error parsing table"
#define LIGOLWXMLREADH_MSGEMTAB "No table type specified"
#define LIGOLWXMLREADH_MSGEENDT "Ending a table without an beginning a table"
#define LIGOLWXMLREADH_MSGETMSM "Table type mismatch"
#define LIGOLWXMLREADH_MSGETNOP "Table not begun for writing"
/*@}*/

/* Forward declarations of MetaIO types */
struct MetaioParseEnvironment;
typedef struct MetaioParseEnvironment* MetaioParseEnv;

/**
 * This structure allows for the association of entries in a MetaDataTable
 * with columns in an xml file.
 * <dl>
 * <dt>name</dt><dd> The name of the column in the XML table.</dd>
 * <dt>pos</dt><dd> The position of this column in the XML table.</dd>
 * <dt>idx</dt><dd> The id number of the column.</dd>
 * </dl>
 *
 */
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagMetaTableDirectory, name));
#endif /* SWIG */
typedef struct
tagMetaTableDirectory
{
  const CHAR *name;
  INT4   pos;
  INT4   idx;
}
MetaTableDirectory;

MetaTableDirectory* XLALCreateMetaTableDir(
    const MetaioParseEnv    env,
    MetadataTableType       table
    );

void
LALCreateMetaTableDir(
    LALStatus              *status,
    MetaTableDirectory    **tableDir,
    const MetaioParseEnv    env,
    MetadataTableType       table
    );

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
    const CHAR          *fileName
    );

#ifdef  __cplusplus
}
#endif

#endif /* _LIGOLWXMLREAD_H */
