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
 * \brief Implementation of the VOTable serializers XML API
 */

#include <string.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

#include <lal/XLALError.h>
#include <lal/LALXML.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXMLVOTableSerializers.h>


#define INT4STR_MAXLEN 15
#define XPATHSTR_MAXLEN 150


/**
 * \brief Serializes a \c LIGOTimeGPS structure into a VOTable XML %node
 *
 * This function takes a \c LIGOTimeGPS structure and serializes it into a VOTable
 * \c RESOURCE %node identified by the given name. The returned \c xmlNode can then be
 * embedded into an existing %node hierarchy or turned into a full VOTable document.
 *
 * \param ltg [in] Pointer to the \c LIGOTimeGPS structure to be serialized
 * \param name [in] Unique identifier of this particular \c LIGOTimeGPS structure instance
 *
 * \return A pointer to a \c xmlNode that holds the VOTable fragment that represents
 * the \c LIGOTimeGPS structure.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * fragment isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \sa XLALCreateVOTableParamNode
 * \sa XLALCreateVOTableResourceNode
 * \sa XLALCreateVOTableXMLFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr XLALLIGOTimeGPS2VOTableNode(const LIGOTimeGPS *const ltg, const char *name)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALLIGOTimeGPS2VOTableNode";
    CHAR gpsSecondsBuffer[INT4STR_MAXLEN] = {0};
    CHAR gpsNanoSecondsBuffer[INT4STR_MAXLEN] = {0};
    xmlNodePtr xmlResourceNode = NULL;
    xmlNodePtr xmlResourceParamNodes[2] = {NULL};

    /* make sure that the shared library is the same as the
     * library version the code was compiled against */
    LIBXML_TEST_VERSION

    /* check and prepare input parameters */
    if(!ltg || LALSnprintf(gpsSecondsBuffer, INT4STR_MAXLEN, "%i", ltg->gpsSeconds) < 0) {
        XLALPrintError("Invalid input parameter: LIGOTimeGPS->gpsSeconds\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!ltg || LALSnprintf(gpsNanoSecondsBuffer, INT4STR_MAXLEN, "%i", ltg->gpsNanoSeconds) < 0) {
        XLALPrintError("Invalid input parameter: LIGOTimeGPS->gpsNanoSeconds\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }

    /* set up RESOURCE node child (first PARAM) */
    xmlResourceParamNodes[0] = XLALCreateVOTableParamNode(GPS_SECONDS, gpsSecondsBuffer);
    if(!xmlResourceParamNodes[0]) {
        XLALPrintError("Couldn't create PARAM node: gpsSeconds\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* set up RESOURCE node child (second PARAM) */
    xmlResourceParamNodes[1] = XLALCreateVOTableParamNode(GPS_NANOSECONDS, gpsNanoSecondsBuffer);
    if(!xmlResourceParamNodes[1]) {
        XLALPrintError("Couldn't create PARAM node: gpsNanoSeconds\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* set up RESOURCE node*/
    xmlResourceNode = XLALCreateVOTableResourceNode("LIGOTimeGPS", name, xmlResourceParamNodes, 2);
    if(!xmlResourceNode) {
        XLALPrintError("Couldn't create RESOURCE node: LIGOTimeGPS\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* return RESOURCE node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlResourceNode;
}


/**
 * \brief Serializes a \c LIGOTimeGPS structure into a VOTable XML string
 *
 * This function takes a \c LIGOTimeGPS structure and serializes it into a full-fledged
 * VOTable XML string containing the serialized structure as the only child element.\n
 * Essentially, this function is just a wrapper for \ref XLALLIGOTimeGPS2VOTableNode and
 * \ref XLALCreateVOTableXMLFromTree followed by a dump of the VOTable document into a
 * string.\n
 *
 * \param ltg [in] Pointer to the \c LIGOTimeGPS structure to be serialized
 * \param name [in] Unique identifier of this particular \c LIGOTimeGPS structure instance
 *
 * \return A pointer to a \c xmlChar (string) that holds the VOTable document containing
 * solely the \c LIGOTimeGPS structure. Please note that the string will be encoded in UTF-8.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * string isn't needed anymore) using \c xmlFree.
 *
 * \sa XLALLIGOTimeGPS2VOTableNode
 * \sa XLALCreateVOTableXMLFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlChar * XLALLIGOTimeGPS2VOTableXML(const LIGOTimeGPS *const ltg, const char *name)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALLIGOTimeGPS2VOTableXML";
    xmlChar *xmlStringBuffer = NULL;
    INT4 xmlStringBufferSize = -1;
    xmlNodePtr xmlTree;

    /* sanity checks */
    if(!ltg) {
        XLALPrintError("Invalid input parameter: ltg\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }

    /* prepare XML serialization */
    xmlThrDefIndentTreeOutput(1);

    /* build VOTable fragment (tree) */
    xmlTree = XLALLIGOTimeGPS2VOTableNode(ltg, name);
    if(xmlTree == NULL) {
        XLALPrintError("VOTable fragment construction failed\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* retrieve full VOTable XML document (string) */
    if(XLALCreateVOTableStringFromTree(xmlTree, &xmlStringBuffer, &xmlStringBufferSize)) {
        /* clean up */
        xmlCleanupParser();

        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* clean up */
    xmlCleanupParser();

    /* return XML string (needs to be xmlFree'd by caller!!!) */
    return xmlStringBuffer;
}


/**
 * \brief Deserializes a \c LIGOTimeGPS structure from a VOTable XML string
 *
 * This function takes a VOTable XML document (string) and deserializes (extracts)
 * the \c LIGOTimeGPS structure identified by the given name.
 *
 * \param xml [in] Pointer to the VOTable XML document (string) containing the structure
 * \param name [in] Unique identifier of the particular \c LIGOTimeGPS structure to be deserialized
 * \param ltg [out] Pointer to an empty \c  LIGOTimeGPS structure to store the deserialized instance
 *
 * \return \c XLAL_SUCCESS if the specified \c LIGOTimeGPS structure could be found and
 * deserialized successfully.
 *
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4 XLALVOTableXML2LIGOTimeGPSByName(const char *xml, const char *name, LIGOTimeGPS *ltg)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALVOTableXML2LIGOTimeGPSByName";
    xmlDocPtr xmlDocument = NULL;
    xmlChar *nodeContent = NULL;
    CHAR xpath[XPATHSTR_MAXLEN] = {0};

    /* sanity checks */
    if(!xml) {
        XLALPrintError("Invalid input parameter: xml\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!ltg) {
        XLALPrintError("Invalid input parameter: ltg\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }

    /* parse XML document */
    xmlDocument = xmlReadMemory(xml, strlen(xml), NULL, "UTF-8", 0);
    if(xmlDocument == NULL) {
        /* clean up */
        xmlCleanupParser();

        XLALPrintError("VOTable document parsing failed\n");
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* prepare XPATH search for LIGOTimeGPS.gpsSeconds */
    if(LALSnprintf(
            xpath,
            XPATHSTR_MAXLEN,
            "//RESOURCE[@utype='LIGOTimeGPS' and @name='%s']/PARAM[@name='gpsSeconds']/@value",
            name) < 0)
    {
        /* clean up */
        xmlFreeDoc(xmlDocument);
        xmlCleanupParser();

        XLALPrintError("XPATH statement construction failed: LIGOTimeGPS.gpsSeconds\n");
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* retrieve LIGOTimeGPS.gpsSeconds */
    nodeContent = (xmlChar *) XLALGetSingleNodeContentByXPath(xmlDocument, xpath);

    /* parse and finally store content */
    if(!nodeContent || sscanf((char*)nodeContent, "%i", &ltg->gpsSeconds) == EOF) {
        /* clean up*/
        xmlFree(nodeContent);
        xmlFreeDoc(xmlDocument);
        xmlCleanupParser();

        XLALPrintError("Invalid node content encountered: gpsSeconds\n");
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* prepare XPATH search for LIGOTimeGPS.gpsNanoSeconds */
    if(LALSnprintf(
            xpath,
            XPATHSTR_MAXLEN,
            "//RESOURCE[@utype='LIGOTimeGPS' and @name='%s']/PARAM[@name='gpsNanoSeconds']/@value",
            name) < 0)
    {
        /* clean up */
        xmlFree(nodeContent);
        xmlFreeDoc(xmlDocument);
        xmlCleanupParser();

        XLALPrintError("XPATH statement construction failed: LIGOTimeGPS.gpsNanoSeconds\n");
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* retrieve LIGOTimeGPS.gpsNanoSeconds */
    nodeContent = (xmlChar *)XLALGetSingleNodeContentByXPath(xmlDocument, xpath);

    /* parse and finally store content */
    if(!nodeContent || sscanf((char*)nodeContent, "%i", &ltg->gpsNanoSeconds) == EOF) {
        /* clean up*/
        xmlFree(nodeContent);
        xmlFreeDoc(xmlDocument);
        xmlCleanupParser();

        XLALPrintError("Invalid node content encountered: gpsNanoSeconds\n");
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* clean up*/
    xmlFree(nodeContent);
    xmlFreeDoc(xmlDocument);
    xmlCleanupParser();

    return XLAL_SUCCESS;
}
