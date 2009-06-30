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
 * \brief Implementation of the common VOTable XML API
 */

/* ---------- includes ---------- */
#include <string.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

#include <lal/LALStdio.h>
#include <lal/XLALError.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXML.h>


/* ---------- defines and macros ---------- */
#define VOTABLE_VERSION     "1.1"
#define VOTABLE_NS_PREFIX   "vot"
#define VOTABLE_NS_URL      "http://www.ivoa.net/xml/VOTable/v"VOTABLE_VERSION
#define VOTABLE_SCHEMA      "http://www.ivoa.net/xml/VOTable/v"VOTABLE_VERSION
#define XPATHSTR_MAXLEN     500
#define TRUE 1
#define FALSE 0

/* ---------- internal prototypes ---------- */
const char* XLALVOTDatatype2String ( VOTABLE_DATATYPE datatype );
const char* XLALVOTParamAttribute2String ( VOTABLE_PARAM_ATTRIBUTE paramAttribute );

/* ---------- function definitions ---------- */

/**
 * \brief Creates a VOTable \c PARAM %node
 *
 * This function creates a VOTable \c PARAM %node with the specified properties.
 *
 * \param name [in] Content of the \c name attribute of the \c PARAM %node (mandatory)
 * \param unit [in] Content of the \c unit attribute of the \c PARAM %node (optional)
 * \param datatype [in] Content of the \c datatype attribute of the \c PARAM %node (mandatory)
 * \param arraysize [in] Content of the \c arraysize attribute of the \c PARAM %node (optional)
 * \param value [in] Content of the \c value attribute of the \c PARAM %node (mandatory, empty value allowed)
 *
 * \return A \c xmlNodePtr that holds the new \c PARAM %node.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * %node isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr XLALCreateVOTParamNode(const char *name,
                                  const char *unit,
                                  VOTABLE_DATATYPE datatype,
                                  const char *arraysize,
                                  const char *value)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALCreateVOTParamNode";
    xmlNodePtr xmlParamNode = NULL;
    static const CHAR *datatypeString;

    /* create node */
    xmlParamNode = xmlNewNode(NULL, CAST_CONST_XMLCHAR("PARAM"));
    if(xmlParamNode == NULL) {
        XLALPrintError("Element instantiation failed: PARAM\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add attributes */
    /* mandatory: name */
    if(!name || strlen(name) <= 0) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Missing mandatory attribute: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!xmlNewProp(xmlParamNode, CAST_CONST_XMLCHAR("name"), CAST_CONST_XMLCHAR(name))) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Attribute instantiation failed: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }
    /* optional: unit */
    if(unit && strlen(unit) > 0) {
        if(!xmlNewProp(xmlParamNode, CAST_CONST_XMLCHAR("unit"), CAST_CONST_XMLCHAR(unit))) {
            /* clean up */
            xmlFreeNode(xmlParamNode);
            XLALPrintError("Attribute instantiation failed: unit\n");
            XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
        }
    }
    /* mandatory: datatype */
    if ( ( datatypeString = XLALVOTDatatype2String ( datatype )) == NULL ) {
      XLALPrintError ("%s: XLALVOTDatatype2String() failed.\n", logReference );
      XLAL_ERROR_NULL ( logReference, XLAL_EFUNC );
    }

    if(!xmlNewProp(xmlParamNode, CAST_CONST_XMLCHAR("datatype"), CAST_CONST_XMLCHAR(datatypeString))) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Attribute instantiation failed: datatype\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }
    /* optional: arraysize */
    if(arraysize && strlen(arraysize) > 0) {
        if(!xmlNewProp(xmlParamNode, CAST_CONST_XMLCHAR("arraysize"), CAST_CONST_XMLCHAR(arraysize))) {
            /* clean up */
            xmlFreeNode(xmlParamNode);
            XLALPrintError("Attribute instantiation failed: arraysize\n");
            XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
        }
    }
    /* mandatory: value (empty value allowed) */
    if(!value) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Missing mandatory attribute: value\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!xmlNewProp(xmlParamNode, CAST_CONST_XMLCHAR("value"), CAST_CONST_XMLCHAR(value))) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Attribute instantiation failed: value\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* return PARAM node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlParamNode;
}


/**
 * \brief Creates a VOTable \c RESOURCE %node
 *
 * This function creates a VOTable \c RESOURCE %node with the specified identifier and assigns
 * the given children to it.
 *
 * \param type [in] Type of the \c RESOURCE %node (typically the \c struct type name)
 * \param identifier [in] Identifier (name) of the \c RESOURCE %node
 * \param children [in] Pointer to an array of \c xmlNodes that are to be assigned as children
 * \param childCount [in] The number of child nodes referenced by \c children
 *
 * \return A \c xmlNodePtr that holds the new \c RESOURCE %node (incl. all children).
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * %node isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr XLALCreateVOTResourceNode(const char *type,
                                     const char *identifier,
                                     const xmlNodePtr childNodeList)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALCreateVOTResourceNode";
    xmlNodePtr xmlResourceNode = NULL;
    xmlNodePtr xmlChildNode = childNodeList;

    /* sanity check */
    if(!type) {
        XLALPrintError("Invalid input parameter: type\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!identifier) {
        XLALPrintError("Invalid input parameter: identifier\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }

    /* create node */
    xmlResourceNode = xmlNewNode(NULL, CAST_CONST_XMLCHAR("RESOURCE"));
    if(xmlResourceNode == NULL) {
        XLALPrintError("Element instantiation failed: RESOURCE\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add attributes */
    if(!xmlNewProp(xmlResourceNode, CAST_CONST_XMLCHAR("utype"), CAST_CONST_XMLCHAR(type))) {
        /* clean up */
        xmlFreeNode(xmlResourceNode);
        XLALPrintError("Attribute instantiation failed: utype\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }
    if(!xmlNewProp(xmlResourceNode, CAST_CONST_XMLCHAR("name"), CAST_CONST_XMLCHAR(identifier))) {
        /* clean up */
        xmlFreeNode(xmlResourceNode);
        XLALPrintError("Attribute instantiation failed: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add children */
    while(xmlChildNode) {
        if(!xmlAddChild(xmlResourceNode, xmlChildNode)) {
            /* clean up */
            xmlFreeNode(xmlResourceNode);
            XLALPrintError("Couldn't add child node to RESOURCE node!\n");
            XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
        }
        /* advance to next sibling in list */
        xmlChildNode = xmlChildNode->next;
    }

    /* return RESOURCE node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlResourceNode;
}


/**
 * \brief Takes a XML fragment (tree) and turns it into a VOTable document
 *
 * This function wraps a given VOTable XML fragment in a \c VOTABLE element to turn it into
 * a valid document. Please make sure that the root element of the given fragment
 * is a valid child of the \c VOTABLE element (VOTable schema 1.1):
 * \li \c DESCRIPTION
 * \li \c COOSYS
 * \li \c PARAM
 * \li \c INFO
 * \li \c RESOURCE
 *
 * \param xmlTree [in] The XML fragment to be turned into a VOTable document
 *
 * \return A pointer to a \c xmlDoc that represents the full VOTable XML document.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * document isn't needed anymore) using \c xmlFreeDoc.
 *
 * \sa XLALCreateVOTStringFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlDocPtr XLALCreateVOTDocumentFromTree(const xmlNodePtr xmlTree)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALCreateVOTDocumentFromTree";
    xmlDocPtr xmlDocument = NULL;
    xmlNodePtr xmlRootElement = NULL;
    xmlNsPtr xmlVOTableNamespace = NULL;
    xmlNsPtr xmlSchemaNamespace = NULL;

    /* make sure that the shared library is the same as the
     * library version the code was compiled against */
    LIBXML_TEST_VERSION

    /* sanity check */
    if(!xmlTree) {
        XLALPrintError("Invalid input parameter: xmlTree\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }

    /* set up XML document */
    xmlDocument = xmlNewDoc(CAST_CONST_XMLCHAR("1.0"));
    if(xmlDocument == NULL) {
        XLALPrintError("VOTable document instantiation failed\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* set up root node */
    xmlRootElement = xmlNewNode(NULL, CAST_CONST_XMLCHAR("VOTABLE"));
    if(xmlRootElement == NULL) {
        /* clean up */
        xmlFreeDoc(xmlDocument);
        XLALPrintError("VOTABLE root element instantiation failed\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add supplemental root node version information */
    if(!xmlNewProp(xmlRootElement, CAST_CONST_XMLCHAR("version"), CAST_CONST_XMLCHAR(VOTABLE_VERSION))) {
        XLALPrintWarning("VOTABLE attribute instantiation failed: version\n");
    }

    /* set up default namespace (required for validation) */
    xmlVOTableNamespace = xmlNewNs(xmlRootElement,
                                   CAST_CONST_XMLCHAR(VOTABLE_NS_URL),
                                   NULL);

    if(xmlVOTableNamespace == NULL) {
        XLALPrintError("VOTABLE namespace instantiation failed\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add supplemental root node schema instance information */
    xmlSchemaNamespace = xmlNewNs(xmlRootElement,
                                  CAST_CONST_XMLCHAR("http://www.w3.org/2001/XMLSchema-instance"),
                                  CAST_CONST_XMLCHAR("xsi"));
    if(!xmlSchemaNamespace) {
        XLALPrintWarning("VOTABLE namespace instantiation failed: xsi\n");
    }
    else if(!xmlNewNsProp(xmlRootElement,
                          xmlSchemaNamespace,
                          CAST_CONST_XMLCHAR("noNamespaceSchemaLocation"),
                          CAST_CONST_XMLCHAR(VOTABLE_SCHEMA)))
    {
        XLALPrintWarning("VOTABLE attribute instantiation failed: xsi:noNamespaceSchemaLocation\n");
    }

    /* append tree to root node */
    if(!xmlAddChild(xmlRootElement, xmlTree)) {
        /* clean up */
        xmlFreeDoc(xmlDocument);
        XLALPrintError("Couldn't append given tree to VOTABLE root element\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* finally, assign root element to document */
    xmlDocSetRootElement(xmlDocument, xmlRootElement);

    /* reconcile default namespace with all document elements */
    if(XLALReconcileDefaultNamespace(xmlRootElement, xmlVOTableNamespace) != XLAL_SUCCESS) {
        /* clean up */
        xmlFreeDoc(xmlDocument);
        XLALPrintError("Default namespace reconciliation failed!\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* return VOTable document (needs to be xmlFreeDoc'd by caller!!!) */
    return xmlDocument;
}


/**
 * \brief Takes a XML fragment (tree) and turns it into a VOTable document string
 *
 * This function takes a VOTable XML fragment and returns a full-fledged VOTable XML string.
 * Please note that all restrictions described for \ref XLALCreateVOTableDocumentFromTree also apply here!
 *
 * \param xmlTree [in] The XML fragment to be turned into a VOTable document
 * \param xmlStringBuffer [out] Pointer to the (uninitialized) buffer that will hold the XML string
 * \param xmlStringBufferSize [out] Pointer to a variable that will hold the size of \c xmlStringBuffer
 *
 * \return \c XLAL_SUCCESS if the specified XML tree could be successfully serialized and dumped into a string.
 * The content will be encoded in UTF-8.\n
 * \b Important: the caller is responsible to free the allocated memory of \c xmlStringBuffer (when the
 * string isn't needed anymore) using \c xmlFree.
 *
 * \sa XLALCreateVOTDocumentFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4 XLALCreateVOTStringFromTree(const xmlNodePtr xmlTree,
                                 xmlChar **xmlStringBuffer,
                                 INT4 *xmlStringBufferSize)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALCreateVOTStringFromTree";
    xmlDocPtr xmlDocument;

    /* sanity check */
    if(!xmlTree) {
        XLALPrintError("Invalid input parameters: xmlTree\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!xmlStringBuffer) {
        XLALPrintError("Invalid input parameters: xmlStringBuffer\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!xmlStringBufferSize) {
        XLALPrintError("Invalid input parameters: xmlStringBufferSize\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }

    /* build VOTable document */
    xmlDocument = XLALCreateVOTDocumentFromTree(xmlTree);
    if(xmlDocument == NULL) {
        XLALPrintError("VOTable document construction failed\n");
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* dump VOTable document to formatted XML string */
    xmlDocDumpFormatMemoryEnc(xmlDocument, xmlStringBuffer, xmlStringBufferSize, "UTF-8", 1);
    if(*xmlStringBufferSize <= 0) {
        /* clean up */
        xmlFreeDoc(xmlDocument);
        XLALPrintError("VOTable document dump failed\n");
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* clean up */
    xmlFreeDoc(xmlDocument);

    return XLAL_SUCCESS;
}


/**
 * \brief Retrieves a specific attribute of a single VOTable \c RESOURCE->PARAM %node relation
 *
 * This function fetches the content of the specified attribute of the specified \c PARAM element,
 * which is a child of the specified \c RESOURCE element, from the given VOTable document.
 *
 * \param xmlDocument [in] The XML document to be searched
 * \param resourceType [in] Value of the \c utype attribute of the \c RESOURCE element to be searched
 * \param resourceName [in] Value of the \c name attribute of the \c RESOURCE element to be searched
 * \param paramName [in] Value of the \c name attribute of the \c PARAM element to be searched
 * \param paramAttribute [in] Attribute of the \c PARAM element to be searched for
 *
 * \return A pointer to a \c xmlChar that holds the content (string) of the specified \c PARAM element
 * attribute. The content will be encoded in UTF-8. In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * string isn't needed anymore) using \c xmlFree.
 *
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlChar * XLALGetSingleVOTResourceParamAttribute(const xmlDocPtr xmlDocument,
                                                 const char *resourceType,
                                                 const char *resourceName,
                                                 const char *paramName,
                                                 VOTABLE_PARAM_ATTRIBUTE paramAttribute)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALGetSingleVOTResourceParamAttribute";
    const CHAR *paramAttributeString = NULL;
    CHAR xpath[XPATHSTR_MAXLEN] = {0};
    static const XML_NAMESPACE xmlVOTableNamespace[1] = {{CAST_CONST_XMLCHAR(VOTABLE_NS_PREFIX), CAST_CONST_XMLCHAR(VOTABLE_NS_URL)}};
    const XML_NAMESPACE_VECTOR xmlNsVector = {xmlVOTableNamespace, 1};

    /* sanity check */
    if(!xmlDocument) {
        XLALPrintError("Invalid input parameters: xmlDocument\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!resourceType) {
        XLALPrintError("Invalid input parameters: resourceType\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!resourceName) {
        XLALPrintError("Invalid input parameters: resourceName\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!paramName) {
        XLALPrintError("Invalid input parameters: paramName\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }


    if ( (paramAttributeString = XLALVOTParamAttribute2String ( paramAttribute )) == NULL ) {
      XLALPrintError ("%s: XLALVOTParamAttribute2String() failed.\n", logReference );
      XLAL_ERROR_NULL ( logReference, XLAL_EFUNC );
    }

    /* prepare XPath search */
    if(snprintf(
            xpath,
            XPATHSTR_MAXLEN,
            "//"VOTABLE_NS_PREFIX":RESOURCE[@utype='%s' and @name='%s']/"VOTABLE_NS_PREFIX":PARAM[@name='%s']/@%s",
            resourceType, resourceName, paramName, paramAttributeString) < 0)
    {
        XLALPrintError("XPath statement construction failed: %s.%s.%s\n", resourceName, paramName, paramAttributeString);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* retrieve specified attribute (content) */
    return (xmlChar *)XLALGetSingleNodeContentByXPath(xmlDocument, xpath, &xmlNsVector);
}



/** Simply returns the string representation of the given VOTABLE_DATATYPE.
 */
const char*
XLALVOTDatatype2String ( VOTABLE_DATATYPE datatype )
{
  static const char *fn = "XLALVOTDatatype2String()";
  const char *datatypeString = NULL;

  switch(datatype)
    {
    case VOT_BOOL:
      datatypeString = "boolean";
      break;
    case VOT_BIT:
      datatypeString = "bit";
      break;
    case VOT_CHAR:
      datatypeString = "char";
      break;
    case VOT_CHAR_UTF:
      datatypeString = "unicodeChar";
      break;
    case VOT_INT1:
      datatypeString = "unsignedByte";
      break;
    case VOT_INT2:
      datatypeString = "short";
      break;
    case VOT_INT4:
      datatypeString = "int";
      break;
    case VOT_INT8:
      datatypeString = "long";
      break;
    case VOT_REAL4:
      datatypeString = "float";
      break;
    case VOT_REAL8:
      datatypeString = "double";
      break;
    case VOT_COMPLEX_REAL4:
      datatypeString = "floatComplex";
      break;
    case VOT_COMPLEX_REAL8:
      datatypeString = "doubleComplex";
      break;
    default:
      XLALPrintError ("%s: invalid datatype passed (%d), has to be within [1, %d].\n", fn, datatype, VOT_DATATYPE_LAST - 1 );
      XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
      break;
    }

  return datatypeString;

} /* XLALVOTDatatype2String() */


/** Simply returns the string representation of the given VOTABLE_PARAM_ATTRIBUTE
 */
const char*
XLALVOTParamAttribute2String ( VOTABLE_PARAM_ATTRIBUTE paramAttribute )
{
  static const char *fn = "XLALVOTParamAttribute2String()";
  const char *paramAttributeString = NULL;

  switch(paramAttribute)
    {
    case VOT_ID:
      paramAttributeString = "ID";
      break;
    case VOT_UNIT:
      paramAttributeString = "unit";
      break;
    case VOT_DATATYPE:
      paramAttributeString = "datatype";
      break;
    case VOT_PRECISION:
      paramAttributeString = "precision";
      break;
    case VOT_WIDTH:
      paramAttributeString = "width";
      break;
    case VOT_REF:
      paramAttributeString = "ref";
      break;
    case VOT_NAME:
      paramAttributeString = "name";
      break;
    case VOT_UCD:
      paramAttributeString = "ucd";
      break;
    case VOT_UTYPE:
      paramAttributeString = "utype";
      break;
    case VOT_ARRAYSIZE:
      paramAttributeString = "arraysize";
      break;
    case VOT_VALUE:
      paramAttributeString = "value";
      break;
    default:
      XLALPrintError ("%s: invalid paramAttribute (%d), must lie within [1, %d].\n", fn, paramAttribute, VOT_PARAM_ATTRIBUTE_LAST - 1 );
      XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
    }

  return paramAttributeString;

} /* XLALVOTParamAttribute2String() */

