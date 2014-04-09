/*
*  Copyright (C) 2007 Jolien Creighton
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
/* Double-include protection */
#ifndef _LALXML_H
#define _LALXML_H

/* C++ protection */
#ifdef __cplusplus
extern "C" {
#endif

#include <libxml/tree.h>
#include <lal/LALDatatypes.h>

/**
 * \addtogroup pkg_XML
 * \author Oliver Bock, Reinhard Prix, Jolien Creighton
 * \brief Module for reading/writing/manipulating XML as well as data serialization.
 */

/**
 * \defgroup LALXML_h Header LALXML.h
 * \ingroup pkg_XML
 * \brief Header file declaring the public XML API
 */
/*@{*/

/**
 * Cast macro to use instead of libxml2's 'BAD_CAST' macro,
 * with are more descriptive name.
 */
#define CAST_XMLCHAR (xmlChar *)

/**
 * Cast macro to use instead of libxml2's BAD_CAST macro for const pointers.
 */
#define CAST_CONST_XMLCHAR (const xmlChar *)


/**
 * \brief This type represents a XML namespace
 *
 * \sa XML_NAMESPACE_VECTOR
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagXML_NAMESPACE, prefix, url));
#endif /* SWIG */
typedef struct tagXML_NAMESPACE {
    const xmlChar *prefix;
    const xmlChar *url;
} XML_NAMESPACE;


/**
 * \brief This type represents a vector of XML namespaces
 *
 * \sa XML_NAMESPACE
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagXML_NAMESPACE_VECTOR, items, count));
#endif /* SWIG */
typedef struct tagXML_NAMESPACE_VECTOR {
    const XML_NAMESPACE *items;
    const int count;
} XML_NAMESPACE_VECTOR;


int XLALXMLFilePrintElements(const char *fname);

const xmlNode *XLALfindNextNamedNode ( const xmlNode *startNode, const xmlChar *nodeName );
int XLALcountNamedNodes ( const xmlNode *startNode, const xmlChar *nodeName, UINT4 *count );

xmlChar * XLALGetSingleNodeContentByXPath(const xmlDocPtr xmlDocument, const char *xpath, const XML_NAMESPACE_VECTOR *xmlNsVector);
INT4 XLALValidateDocumentByInternalSchema(const xmlDocPtr xmlDocument);
INT4 XLALValidateDocumentByExternalSchema(const xmlDocPtr xmlDocument, const xmlChar *url);
INT4 XLALReconcileDefaultNamespace(const xmlNodePtr xmlRootElement, const xmlNsPtr xmlNamespace);

xmlDocPtr XLALXMLString2Doc ( const char *xmlString );
xmlChar *XLALXMLDoc2String ( xmlDoc *xmlDocument );

/*@}*/

/* C++ protection */
#ifdef __cplusplus
}
#endif

/* Double-include protection */
#endif
