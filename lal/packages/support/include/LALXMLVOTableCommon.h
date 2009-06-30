/*
*  Copyright (C) 2009 Oliver Bock
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
 * \ingroup XML
 * \brief Header file declaring the public common VOTable XML API
 */

/* Double-include protection */
#ifndef _LALXMLVOTABLECOMMON_H
#define _LALXMLVOTABLECOMMON_H

/* C++ protection */
#ifdef __cplusplus
extern "C" {
#endif


#include <libxml/tree.h>


/**
 * \brief List of all supported VOTable data types
 *
 * This enumeration contains all supported VOTable data types.
 * They are used for \c PARAM elements for instance.
 *
 * \sa XLALCreateVOTParamNode
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
typedef enum {
    VOT_BOOL = 1,
    VOT_BIT,
    VOT_CHAR,
    VOT_CHAR_UTF,
    VOT_INT1,
    VOT_INT2,
    VOT_INT4,
    VOT_INT8,
    VOT_REAL4,
    VOT_REAL8,
    VOT_COMPLEX_REAL4,
    VOT_COMPLEX_REAL8,
    VOT_DATATYPE_LAST
} VOTABLE_DATATYPE;

/**
 * \brief List of all supported VOTable \c PARAM element attributes
 *
 * This enumeration contains all supported attributes of the
 * VOTable \c PARAM element
 *
 * \sa XLALGetSingleVOTResourceParamAttribute
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
typedef enum {
    VOT_ID = 1,
    VOT_UNIT,
    VOT_DATATYPE,
    VOT_PRECISION,
    VOT_WIDTH,
    VOT_REF,
    VOT_NAME,
    VOT_UCD,
    VOT_UTYPE,
    VOT_ARRAYSIZE,
    VOT_VALUE,
    VOT_PARAM_ATTRIBUTE_LAST
} VOTABLE_PARAM_ATTRIBUTE;


xmlNodePtr XLALCreateVOTParamNode(const char *name,
                                  const char *unit,
                                  VOTABLE_DATATYPE datatype,
                                  const char *arraysize,
                                  const char *value);

xmlNodePtr XLALCreateVOTResourceNode(const char *type,
                                     const char *identifier,
                                     const xmlNodePtr childNodeList);

xmlDocPtr XLALCreateVOTDocumentFromTree(const xmlNodePtr xmlTree);

INT4 XLALCreateVOTStringFromTree(const xmlNodePtr xmlTree,
                                 xmlChar **xmlStringBuffer,
                                 INT4 *xmlStringBufferSize);

xmlChar *XLALGetSingleVOTResourceParamAttribute(const xmlDocPtr xmlDocument,
                                                const char *resourceType,
                                                const char *resourceName,
                                                const char *paramName,
                                                VOTABLE_PARAM_ATTRIBUTE paramAttribute);


/* C++ protection */
#ifdef  __cplusplus
}
#endif

/* Double-include protection */
#endif
