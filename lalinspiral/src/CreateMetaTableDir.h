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

/**
 * \author Brown, D. A. and Fairhurst, S.
 * \file
 * \ingroup lalinspiral
 * \brief Provides functions for reading LIGO lightweight XML files to LIGO metadata database tables.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/CreateMetaTableDir.h>
 * \endcode
 *
 */

#ifndef _CREATEMETATABLEDIR_H
#define _CREATEMETATABLEDIR_H

#include <lal/LIGOLwXMLRead.h>

#ifdef  __cplusplus
extern "C" {
#endif

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
    struct MetaioParseEnvironment *const env,
    MetadataTableType       table
    );

#ifdef  __cplusplus
}
#endif

#endif /* _CREATEMETATABLEDIR_H */
