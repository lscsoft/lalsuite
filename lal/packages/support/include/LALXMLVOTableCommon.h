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

/* ---------- exported includes ---------- */
#include <libxml/tree.h>
#include <libxml/xpath.h>

/* ---------- exported defines and constants ---------- */
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
    VOT_COMPLEX8,
    VOT_COMPLEX16,
    VOT_DATATYPE_LAST
} VOTABLE_DATATYPE;

/**
 * \brief List of all supported VOTable element attributes
 *
 * This enumeration contains all supported attributes of the
 * VOTable \c PARAM and \c FIELD elements
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
    VOT_ATTRIBUTE_LAST
} VOTABLE_ATTRIBUTE;


/** List of supported VOTable "leaf" elements
 *
 */
typedef enum {
  VOT_RESOURCE = 1,
  VOT_TABLE,
  VOT_STREAM,
  VOT_PARAM,
  VOT_FIELD,
  VOT_ELEMENT_LAST
} VOTABLE_ELEMENT;


/**
 * \brief List of all supported VOTable Table Serialization types.
 *
 * Note: this does not yet contain all serialization types allowed
 * within the VOTable standard, such as exteran Fits files etc,
 * but this is left for future extensions, if required.
 *
 * \sa XLALCreateVOTTableNode
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
typedef enum {
  VOT_SERIALIZE_TABLEDATA = 1,	/**< embedded TABLEDATA inside TABLE element */
  VOT_SERIALIZE_BINARY, 	/**< external binary stream, referenced within TABLE element */
  VOT_SERIALIZE_LAST
} VOTABLE_SERIALIZATION_TYPE;


/* ---------- exported API datatypes ---------- */

/** Type holding the attributes of one FIELD node
 * Note: currently this only holds the FIELD attributes that are actually used,
 * but this can further extended as needed. See
 * http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.pdf
 * for a complete list of allowed FIELD attributes.
 *
 */
typedef struct {
  xmlChar *name;		/**< name attribute [required] */
  VOTABLE_DATATYPE datatype;	/**< datatype attribute [required] */
  xmlChar *unit;		/**< unit attribute [optional] */
  xmlChar *arraysize;		/**< arraysize attribute [optional] */
} VOTField;

/** A standard vector of VOTFields
 */
typedef struct {
  UINT4 length;		/**< number of VOTFields */
  VOTField *data;	/**< array of VOTFields */
} VOTFieldVector;


/* ---------- exported API prototypes ---------- */

xmlNodePtr XLALCreateVOTParamNode(const char *name,
                                  const char *unit,
                                  VOTABLE_DATATYPE datatype,
                                  const char *arraysize,
                                  const char *value);

xmlNodePtr
XLALCreateVOTFieldNode ( const char *name,
                         const char *unit,
                         VOTABLE_DATATYPE datatype,
                         const char *arraysize
                         );

xmlNodePtr
XLALCreateVOTTabledataNode ( xmlNode *fieldNodeList, UINT4 numRows, const char *fmt, ... );

xmlNodePtr
XLALCreateVOTTableNode ( const char *name, xmlNode *fieldNodeList, xmlNode *dataContentNode );


xmlNodePtr XLALCreateVOTResourceNode(const char *type,
                                     const char *identifier,
                                     const xmlNodePtr childNodeList);

xmlDoc *XLALCreateVOTDocFromTree(xmlNodePtr xmlTree, BOOLEAN reconcileNamespace );


VOTFieldVector *XLALReadVOTFIELDNodes ( const xmlDocPtr xmlDocument, const CHAR *resourcePath );


void *
XLALReadVOTTabledataSimpleColumn ( const xmlDocPtr xmlDocument,
                                   const CHAR *resourcePath,
                                   UINT4 column,
                                   VOTABLE_DATATYPE datatype,
                                   UINT4 *numRows
                                   );


CHAR *
XLALReadVOTAttributeFromNamedElement ( const xmlDocPtr xmlDocument,
                                       const char *resourcePath,
                                       const char *elementName,
                                       VOTABLE_ELEMENT elementType,
                                       VOTABLE_ATTRIBUTE attrib
                                       );

xmlNodeSet *
XLALFindVOTElementsAtPath ( const xmlDocPtr xmlDocument,
                            const CHAR *extResourcePath
                            );


CHAR *XLALCreateVOTStringFromTree ( xmlNodePtr xmlTree );

VOTFieldVector *XLALCreateVOTFieldVector ( UINT4 numFields );
void XLALDestroyVOTFieldVector ( VOTFieldVector *vect );

const char* XLALVOTDatatype2String ( VOTABLE_DATATYPE datatype );
VOTABLE_DATATYPE XLALVOTString2Datatype ( const CHAR *datatypeString );
const char* XLALVOTElement2String ( VOTABLE_ELEMENT element );
const char* XLALVOTAttribute2String ( VOTABLE_ATTRIBUTE elementAttribute );




/* C++ protection */
#ifdef  __cplusplus
}
#endif

/* Double-include protection */
#endif
