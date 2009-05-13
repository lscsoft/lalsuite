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

/**
 * \file
 * \ingroup XML
 * \brief Implementation of the XML API
 */

#include <stdio.h>
#include <string.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

#include <lal/XLALError.h>
#include <lal/LALXML.h>


static void print_element_names(xmlNode *node)
{
    xmlNode *cur;
    for (cur = node; cur; cur = cur->next) {
        if (cur->type == XML_ELEMENT_NODE)
            printf("node type: Element, name: %s\n", cur->name);
        print_element_names(cur->children);
    }
    return;
}

int XLALXMLFilePrintElements(const char *fname)
{
    static const char file[] = "XLALXMLFilePrintElements";
    xmlDoc  *doc;

    /* make sure that the shared library is the same as the
     * library version the code was compiled against */
    LIBXML_TEST_VERSION

    if (!(doc = xmlReadFile(fname, NULL, 0)))
        XLAL_ERROR(file, XLAL_EIO);
    print_element_names(xmlDocGetRootElement(doc));
    xmlFreeDoc(doc);
    xmlCleanupParser(); /* free global variables in parser */
    return 0;
}


/**
 * \brief Performs a XPath search on a XML document to retrieve the content of a single %node
 *
 * This function searches the given XML document using the given XPath statement.
 * The XPath statement \b must be specified in such a way that at most a single %node
 * will be found.
 *
 * \param xmlDocument [in] The XML document to be searched
 * \param xpath [in] The XPath statement to be used to search the given XML document
 *
 * \return A pointer to a \c xmlChar that holds the content (string) of the %node
 * specified by the given XPath statement. The content will be encoded in UTF-8.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * string isn't needed anymore) using \c xmlFree.
 *
 * \sa XLALGetSingleVOTableResourceParamValue
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlChar * XLALGetSingleNodeContentByXPath(const xmlDocPtr xmlDocument, const char *xpath)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALGetSingleNodeContentByXPath";
    xmlXPathContextPtr xpathCtx = NULL;
    xmlChar *xpathExpr = NULL;
    xmlXPathObjectPtr xpathObj = NULL;
    xmlNodeSetPtr xmlNodes = NULL;
    xmlChar *nodeContent = NULL;
    INT4 nodeCount;

    /* sanity checks */
    if(!xmlDocument) {
        XLALPrintError("Invalid input parameter: xmlDocument\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!xpath || strlen(xpath) <= 0) {
        XLALPrintError("Invalid input parameter: xpath\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }

    /* prepare xpath context */
    xpathCtx = xmlXPathNewContext(xmlDocument);
    if(xpathCtx == NULL) {
        XLALPrintError("XPATH context instantiation failed\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* prepare xpath expression */
    xpathExpr = xmlCharStrdup(xpath);
    if(xpathExpr == NULL) {
        /* clean up */
        xmlXPathFreeContext(xpathCtx);

        XLALPrintError("XPATH statement preparation failed\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* run xpath query */
    xpathObj = xmlXPathEvalExpression(xpathExpr, xpathCtx);
    if(xpathObj == NULL) {
        /* clean up */
        xmlFree(xpathExpr);
        xmlXPathFreeContext(xpathCtx);

        XLALPrintError("XPATH evaluation failed\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* retrieve node set returned by xpath query */
    xmlNodes = xpathObj->nodesetval;

    /* how many nodes did we find? */
    nodeCount = (xmlNodes) ? xmlNodes->nodeNr : 0;
    if(nodeCount <= 0) {
        /* clean up */
        xmlXPathFreeObject(xpathObj);
        xmlFree(xpathExpr);
        xmlXPathFreeContext(xpathCtx);

        XLALPrintError("XPATH search didn't return any nodes\n");
        XLAL_ERROR_NULL(logReference, XLAL_EDOM);
    }
    else if(nodeCount > 1) {
        /* clean up */
        xmlXPathFreeObject(xpathObj);
        xmlFree(xpathExpr);
        xmlXPathFreeContext(xpathCtx);

        XLALPrintError("XPATH search did return %i nodes where only 1 was expected\n", nodeCount);
        XLAL_ERROR_NULL(logReference, XLAL_EDOM);
    }
    else {
        nodeContent = xmlNodeListGetString(xmlDocument, xmlNodes->nodeTab[0]->xmlChildrenNode, 1);
    }

    /* clean up */
    xmlXPathFreeObject(xpathObj);
    xmlFree(xpathExpr);
    xmlXPathFreeContext(xpathCtx);

    /* return node content (needs to be xmlFree'd by caller!!!) */
    return nodeContent;
}
