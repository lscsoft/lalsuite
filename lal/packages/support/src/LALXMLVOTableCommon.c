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

#include <string.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

#include <lal/XLALError.h>
#include <lal/LALXMLVOTableCommon.h>


#define XPATHSTR_MAXLEN         500


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
xmlNodePtr XLALCreateVOTableParamNode(const char *name,
                                      const char *unit,
                                      VOTABLE_DATATYPE datatype,
                                      const char *arraysize,
                                      const char *value)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALCreateVOTableParamNode";
    xmlNodePtr xmlParamNode = NULL;
    static const CHAR *datatypeString;

    /* create node and add attributes*/
    xmlParamNode = xmlNewNode(NULL, BAD_CAST("PARAM"));
    if(xmlParamNode == NULL) {
        XLALPrintError("Element instantiation failed: PARAM\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }
    /* mandatory: name */
    if(!name || strlen(name) <= 0) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Missing mandatory attribute: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!xmlNewProp(xmlParamNode, BAD_CAST("name"), BAD_CAST(name))) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Attribute instantiation failed: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }
    /* optional: unit */
    if(unit && strlen(unit) > 0) {
        if(!xmlNewProp(xmlParamNode, BAD_CAST("unit"), BAD_CAST(unit))) {
            /* clean up */
            xmlFreeNode(xmlParamNode);
            XLALPrintError("Attribute instantiation failed: unit\n");
            XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
        }
    }
    /* mandatory: datatype */
    if(!datatype) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Missing mandatory attribute: datatype\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    else {
        switch(datatype) {
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
                datatypeString = "UNKNOWN";
        }
    }
    if(!xmlNewProp(xmlParamNode, BAD_CAST("datatype"), BAD_CAST(datatypeString))) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Attribute instantiation failed: datatype\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }
    /* optional: arraysize */
    if(arraysize && strlen(arraysize) > 0) {
        if(!xmlNewProp(xmlParamNode, BAD_CAST("arraysize"), BAD_CAST(arraysize))) {
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
    if(!xmlNewProp(xmlParamNode, BAD_CAST("value"), BAD_CAST(value))) {
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
xmlNodePtr XLALCreateVOTableResourceNode(const char *type,
                                         const char *identifier,
                                         const xmlNodePtr *children,
                                         const INT4 childCount)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALCreateVOTableResourceNode";
    xmlNodePtr xmlResourceNode = NULL;
    INT4 i = 0;

    xmlResourceNode = xmlNewNode(NULL, BAD_CAST("RESOURCE"));
    if(xmlResourceNode == NULL) {
        XLALPrintError("Element instantiation failed: RESOURCE\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }
    if(!xmlNewProp(xmlResourceNode, BAD_CAST("utype"), BAD_CAST(type))) {
        /* clean up */
        xmlFreeNode(xmlResourceNode);

        XLALPrintError("Attribute instantiation failed: utype\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }
    if(!xmlNewProp(xmlResourceNode, BAD_CAST("name"), BAD_CAST(identifier))) {
        /* clean up */
        xmlFreeNode(xmlResourceNode);

        XLALPrintError("Attribute instantiation failed: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add children */
    for(i = 0; i < childCount; ++i) {
        if(!xmlAddChild(xmlResourceNode, children[i])) {
            /* clean up */
            xmlFreeNode(xmlResourceNode);

            XLALPrintError("Couldn't add child node to RESOURCE node: #%i\n", i);
            XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
        }
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
 * \sa XLALCreateVOTableStringFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlDocPtr XLALCreateVOTableXMLFromTree(const xmlNodePtr xmlTree)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALCreateVOTableXMLFromTree";
    xmlDocPtr xmlDocument = NULL;
    xmlNodePtr xmlRootNode = NULL;

    /* make sure that the shared library is the same as the
     * library version the code was compiled against */
    LIBXML_TEST_VERSION

    /* sanity check */
    if(!xmlTree) {
        XLALPrintError("Invalid input parameter: xmlTree\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }

    /* set up XML document */
    xmlDocument = xmlNewDoc(BAD_CAST("1.0"));
    if(xmlDocument == NULL) {
        XLALPrintError("VOTable document instantiation failed\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* set up root node */
    xmlRootNode = xmlNewNode(NULL, BAD_CAST("VOTABLE"));
    if(xmlRootNode == NULL) {
        /* clean up */
        xmlFreeDoc(xmlDocument);

        XLALPrintError("VOTable root element instantiation failed\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    xmlDocSetRootElement(xmlDocument, xmlRootNode);

    /* append tree to root node */
    if(!xmlAddChild(xmlRootNode, xmlTree)) {
        /* clean up */
        xmlFreeNode(xmlRootNode);
        xmlFreeDoc(xmlDocument);

        XLALPrintError("Couldn't append given tree to VOTable root element\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* return VOTable document (needs to be xmlFreeDoc'd by caller!!!) */
    return xmlDocument;
}


/**
 * \brief Takes a XML fragment (tree) and turns it into a VOTable document string
 *
 * This function takes a VOTable XML fragment and returns a full-fledged VOTable XML string.
 * Please note that all restrictions described for \ref XLALCreateVOTableXMLFromTree also apply here!
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
 * \sa XLALCreateVOTableXMLFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4 XLALCreateVOTableStringFromTree(const xmlNodePtr xmlTree,
                                     xmlChar **xmlStringBuffer,
                                     INT4 *xmlStringBufferSize)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALCreateVOTableStringFromTree";
    xmlDocPtr xmlDocument;

    /* build VOTable document */
    xmlDocument = XLALCreateVOTableXMLFromTree(xmlTree);
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
 * \brief Retrieves the content of a single VOTable \c RESOURCE->PARAM %node relation
 *
 * This function fetches the content of the \c value attribute of the specified \c PARAM element,
 * which is a child of the specified \c RESOURCE element, from the given VOTable document.
 *
 * \param xmlDocument [in] The XML document to be searched
 * \param resourceType [in] Value of the \c utype attribute of the \c RESOURCE element to be searched
 * \param resourceName [in] Value of the \c name attribute of the \c RESOURCE element to be searched
 * \param paramName [in] Value of the \c name attribute of the \c PARAM element to be searched
 *
 * \return A pointer to a \c xmlChar that holds the content (string) of the \c PARAM element's
 * value attribute. The content will be encoded in UTF-8. In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * string isn't needed anymore) using \c xmlFree.
 *
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlChar * XLALGetSingleVOTableResourceParamValue(const xmlDocPtr xmlDocument,
                                                 const char *resourceType,
                                                 const char *resourceName,
                                                 const char *paramName)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALGetSingleVOTableResourceParamValue";
    CHAR xpath[XPATHSTR_MAXLEN] = {0};

    /* prepare XPath search */
    if(LALSnprintf(
            xpath,
            XPATHSTR_MAXLEN,
            "//RESOURCE[@utype='%s' and @name='%s']/PARAM[@name='%s']/@value",
            resourceType, resourceName, paramName) < 0)
    {
        XLALPrintError("XPath statement construction failed: %s.%s\n", resourceName, paramName);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* retrieve LIGOTimeGPS.gpsSeconds */
    return (xmlChar *)XLALGetSingleNodeContentByXPath(xmlDocument, xpath);
}
