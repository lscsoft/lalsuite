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


#define INT4STR_MAXLEN          15
#define REAL8STR_MAXLEN         25
#define NAMESTR_MAXLEN          50
#define PULSARSPINSTR_MAXLEN    4


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
 * \sa XLALCreateVOTableDocumentFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr XLALLIGOTimeGPS2VOTableNode(const LIGOTimeGPS *const ltg, const char *name)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALLIGOTimeGPS2VOTableNode";
    xmlNodePtr xmlParentNode = NULL;
    xmlNodePtr xmlChildNode = NULL;
    xmlNodePtr xmlChildNodeList = NULL;

    CHAR gpsSecondsBuffer[INT4STR_MAXLEN] = {0};
    CHAR gpsNanoSecondsBuffer[INT4STR_MAXLEN] = {0};

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
    xmlChildNode = XLALCreateVOTableParamNode("gpsSeconds",
                                              "s",
                                              VOT_INT4,
                                              NULL,
                                              gpsSecondsBuffer);
    if(!xmlChildNode) {
        XLALPrintError("Couldn't create PARAM node: , %s.gpsSeconds\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* initialize child node list with first child */
    xmlChildNodeList = xmlChildNode;

    /* set up RESOURCE node child (second PARAM) */
    xmlChildNode = XLALCreateVOTableParamNode("gpsNanoSeconds",
                                              "ns",
                                              VOT_INT4,
                                              NULL,
                                              gpsNanoSecondsBuffer);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.gpsNanoSeconds\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add child as first sibling to child node list */
    xmlChildNodeList->next = xmlChildNode;

    /* set up RESOURCE node*/
    xmlParentNode = XLALCreateVOTableResourceNode("LIGOTimeGPS", name, xmlChildNodeList);
    if(!xmlParentNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: %s\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* return RESOURCE node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlParentNode;
}


/**
 * \brief Deserializes a \c LIGOTimeGPS structure from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts)
 * the \c LIGOTimeGPS structure identified by the given name.
 *
 * \param xml [in] Pointer to the VOTable XML document containing the structure
 * \param name [in] Unique identifier of the particular \c LIGOTimeGPS structure to be deserialized
 * \param ltg [out] Pointer to an empty \c  LIGOTimeGPS structure to store the deserialized instance
 *
 * \return \c XLAL_SUCCESS if the specified \c LIGOTimeGPS structure could be found and
 * deserialized successfully.
 *
 * \sa XLALVOTableXML2LIGOTimeGPSByName
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4 XLALVOTableDoc2LIGOTimeGPSByName(const xmlDocPtr xmlDocument, const char *name, LIGOTimeGPS *ltg)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALVOTableDoc2LIGOTimeGPSByName";
    xmlChar *nodeContent = NULL;

    /* sanity checks */
    if(!xmlDocument) {
        XLALPrintError("Invalid input parameter: xmlDocument\n");
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

    /* retrieve LIGOTimeGPS.gpsSeconds */
    nodeContent = (xmlChar *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, "LIGOTimeGPS", name, "gpsSeconds", VOT_VALUE);

    /* parse and finally store content */
    if(!nodeContent || sscanf((char*)nodeContent, "%i", &ltg->gpsSeconds) == EOF) {
        /* clean up*/
        if(nodeContent) xmlFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.gpsSeconds\n", name);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* retrieve LIGOTimeGPS.gpsNanoSeconds */
    xmlFree(nodeContent);
    nodeContent = (xmlChar *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, "LIGOTimeGPS", name, "gpsNanoSeconds", VOT_VALUE);

    /* parse and finally store content */
    if(!nodeContent || sscanf((char*)nodeContent, "%i", &ltg->gpsNanoSeconds) == EOF) {
        /* clean up*/
        if(nodeContent) xmlFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.gpsNanoSeconds\n", name);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* clean up*/
    xmlFree(nodeContent);

    return XLAL_SUCCESS;
}


/**
 * \brief Serializes a \c BinaryOrbitParams structure into a VOTable XML %node
 *
 * This function takes a \c BinaryOrbitParams structure and serializes it into a VOTable
 * \c RESOURCE %node identified by the given name. The returned \c xmlNode can then be
 * embedded into an existing %node hierarchy or turned into a full VOTable document.
 *
 * \param bop [in] Pointer to the \c BinaryOrbitParams structure to be serialized
 * \param name [in] Unique identifier of this particular \c BinaryOrbitParams structure instance
 *
 * \return A pointer to a \c xmlNode that holds the VOTable fragment that represents
 * the \c BinaryOrbitParams structure.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * fragment isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \sa XLALCreateVOTableParamNode
 * \sa XLALCreateVOTableResourceNode
 * \sa XLALCreateVOTableDocumentFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr XLALBinaryOrbitParams2VOTableNode(const BinaryOrbitParams *const bop, const char *name)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALBinaryOrbitParams2VOTableNode";
    CHAR childNameTp[NAMESTR_MAXLEN] = {0};
    xmlNodePtr xmlParentNode = NULL;
    xmlNodePtr xmlChildNode = NULL;
    xmlNodePtr xmlChildNodeList = NULL;

    CHAR argp[REAL8STR_MAXLEN] = {0};
    CHAR asini[REAL8STR_MAXLEN] = {0};
    CHAR ecc[REAL8STR_MAXLEN] = {0};
    CHAR period[REAL8STR_MAXLEN] = {0};

    /* make sure that the shared library is the same as the
     * library version the code was compiled against */
    LIBXML_TEST_VERSION

    /* check and prepare input parameters */
    if(!bop || LALSnprintf(argp, REAL8STR_MAXLEN, "%g", bop->argp) < 0) {
        XLALPrintError("Invalid input parameter: BinaryOrbitParams->argp\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!bop || LALSnprintf(asini, REAL8STR_MAXLEN, "%g", bop->asini) < 0) {
        XLALPrintError("Invalid input parameter: BinaryOrbitParams->asini\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!bop || LALSnprintf(ecc, REAL8STR_MAXLEN, "%g", bop->ecc) < 0) {
        XLALPrintError("Invalid input parameter: BinaryOrbitParams->ecc\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!bop || LALSnprintf(period, REAL8STR_MAXLEN, "%g", bop->period) < 0) {
        XLALPrintError("Invalid input parameter: BinaryOrbitParams->period\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }

    /* set up PARAM node (argp) */
    xmlChildNode = XLALCreateVOTableParamNode("argp",
                                              "rad",
                                              VOT_REAL8,
                                              NULL,
                                              argp);
    if(!xmlChildNode) {
        XLALPrintError("Couldn't create PARAM node: %s.argp\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* initialize child node list with first node */
    xmlChildNodeList = xmlChildNode;

    /* set up PARAM node (asini) */
    xmlChildNode = XLALCreateVOTableParamNode("asini",
                                              "s",
                                              VOT_REAL8,
                                              NULL,
                                              asini);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.asini\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add child as first sibling to child node list */
    xmlChildNodeList->next = xmlChildNode;

    /* set up PARAM node (ecc) */
    xmlChildNode = XLALCreateVOTableParamNode("ecc",
                                              NULL,
                                              VOT_REAL8,
                                              NULL,
                                              ecc);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.ecc\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add child as second sibling to child node list */
    xmlChildNodeList->next->next = xmlChildNode;

    /* set up PARAM node (period) */
    xmlChildNode = XLALCreateVOTableParamNode("period",
                                              "s",
                                              VOT_REAL8,
                                              NULL,
                                              period);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.period\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add child as third sibling to child node list */
    xmlChildNodeList->next->next->next = xmlChildNode;

    /* compile child name attribute (tp) */
    if(!strncpy(childNameTp, name, NAMESTR_MAXLEN) || !strncat(childNameTp, ".tp", NAMESTR_MAXLEN)) {
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Child attribute preparation failed: BinaryOrbitParams.tp\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* set up RESOURCE node (tp)*/
    xmlChildNode = XLALLIGOTimeGPS2VOTableNode(&bop->tp, childNameTp);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: %s.tp\n", childNameTp);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add child as fourth sibling to child node list */
    xmlChildNodeList->next->next->next->next = xmlChildNode;

    /* set up parent RESOURCE node*/
    xmlParentNode = XLALCreateVOTableResourceNode("BinaryOrbitParams", name, xmlChildNodeList);
    if(!xmlParentNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: %s\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* return RESOURCE node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlParentNode;
}


/**
 * \brief Deserializes a \c BinaryOrbitParams structure from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts)
 * the \c BinaryOrbitParams structure identified by the given name.
 *
 * \param xmlDoc [in] Pointer to the VOTable XML document containing the structure
 * \param name [in] Unique identifier of the particular \c BinaryOrbitParams structure to be deserialized
 * \param bop [out] Pointer to an empty \c  BinaryOrbitParams structure to store the deserialized instance
 *
 * \return \c XLAL_SUCCESS if the specified \c BinaryOrbitParams structure could be found and
 * deserialized successfully.
 *
 * \sa XLALVOTableXML2BinaryOrbitParamsByName
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4 XLALVOTableDoc2BinaryOrbitParamsByName(const xmlDocPtr xmlDocument, const char *name, BinaryOrbitParams *bop)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALVOTableDoc2BinaryOrbitParamsByName";
    CHAR childNameTp[NAMESTR_MAXLEN] = {0};
    xmlChar *nodeContent = NULL;

    /* sanity checks */
    if(!xmlDocument) {
        XLALPrintError("Invalid input parameter: xmlDocument\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!bop) {
        XLALPrintError("Invalid input parameter: bop\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }

    /* compile child name attribute (tp) */
    if(!strncpy(childNameTp, name, NAMESTR_MAXLEN) || !strncat(childNameTp, ".tp", NAMESTR_MAXLEN)) {
        XLALPrintError("Child attribute preparation failed: BinaryOrbitParams.tp\n");
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* retrieve BinaryOrbitParams.tp */
    if(XLALVOTableDoc2LIGOTimeGPSByName(xmlDocument, childNameTp, &bop->tp)) {
        XLALPrintError("Error parsing XML document content: %s.tp\n", childNameTp);
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* retrieve BinaryOrbitParams.argp */
    nodeContent = (xmlChar *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, "BinaryOrbitParams", name, "argp", VOT_VALUE);

    /* parse and finally store content */
    if(!nodeContent || sscanf((char*)nodeContent, "%lf", &bop->argp) == EOF) {
        /* clean up*/
        if(nodeContent) xmlFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.argp\n", name);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* retrieve BinaryOrbitParams.asini */
    xmlFree(nodeContent);
    nodeContent = (xmlChar *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, "BinaryOrbitParams", name, "asini", VOT_VALUE);

    /* parse and finally store content */
    if(!nodeContent || sscanf((char*)nodeContent, "%lf", &bop->asini) == EOF) {
        /* clean up*/
        if(nodeContent) xmlFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.asini\n", name);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* retrieve PulsarDopplerParams.ecc */
    xmlFree(nodeContent);
    nodeContent = (xmlChar *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, "BinaryOrbitParams", name, "ecc", VOT_VALUE);

    /* parse and finally store content */
    if(!nodeContent || sscanf((char*)nodeContent, "%lf", &bop->ecc) == EOF) {
        /* clean up*/
        if(nodeContent) xmlFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.ecc\n", name);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* retrieve PulsarDopplerParams.period */
    xmlFree(nodeContent);
    nodeContent = (xmlChar *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, "BinaryOrbitParams", name, "period", VOT_VALUE);

    /* parse and finally store content */
    if(!nodeContent || sscanf((char*)nodeContent, "%lf", &bop->period) == EOF) {
        /* clean up*/
        if(nodeContent) xmlFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.period\n", name);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* clean up*/
    xmlFree(nodeContent);

    return XLAL_SUCCESS;
}


/**
 * \brief Serializes a \c PulsarSpins array into a VOTable XML %node
 *
 * This function takes a \c PulsarSpins array and serializes it into a VOTable
 * \c PARAM %node identified by the given name. The returned \c xmlNode can then be
 * embedded into an existing %node hierarchy.
 *
 * \param spins [in] Pointer to the \c PulsarSpins array to be serialized
 * \param name [in] Unique identifier of this particular \c PulsarSpins array instance
 *
 * \return A pointer to a \c xmlNode that holds the VOTable fragment that represents
 * the \c PulsarSpins array. In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * fragment isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \sa XLALCreateVOTableParamNode
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr XLALPulsarSpins2VOTableNode(const PulsarSpins *const spins, const char *name)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALPulsarSpins2VOTableNode";
    xmlNodePtr xmlParamNode = NULL;
    int i;

    CHAR spinArraySize[PULSARSPINSTR_MAXLEN] = {0};
    CHAR spinArrayString[PULSAR_MAX_SPINS*REAL8STR_MAXLEN] = {0};
    CHAR spinArrayStringItem[REAL8STR_MAXLEN] = {0};

    /* sanity checks */
    if(!spins || !(*spins)) {
        XLALPrintError("Invalid input parameter: spins\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(sizeof(*spins)/sizeof(REAL8) != PULSAR_MAX_SPINS) {
        XLALPrintError("Invalid input parameter: spins (actual size %i differs from defined size %i)\n", sizeof(*spins)/sizeof(REAL8), PULSAR_MAX_SPINS);
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }

    /* parse input array */
    for(i = 0; i < PULSAR_MAX_SPINS; ++i) {
        if(LALSnprintf(spinArrayStringItem, REAL8STR_MAXLEN, "%g", (*spins)[i]) < 0) {
            XLALPrintError("Invalid input parameter: spins[%i]\n", i);
            XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
        }
        if(!strncat(spinArrayString, spinArrayStringItem, REAL8STR_MAXLEN)) {
            XLALPrintError("Couldn't serialize parameter: spins[%i]\n", i);
            XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
        }
        /* add delimiter (SPACE)*/
        if(i<PULSAR_MAX_SPINS-1 && !strncat(spinArrayString, " ", 1)) {
            XLALPrintError("Couldn't add delimiter to parameter: spins[%i]\n", i);
            XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
        }
    }

    /* set array size attribute */
    if(LALSnprintf(spinArraySize, PULSARSPINSTR_MAXLEN, "%i", PULSAR_MAX_SPINS) < 0) {
        XLALPrintError("Couldn't prepare attribute: arraysize\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* set up PARAM node */
    xmlParamNode = XLALCreateVOTableParamNode(name,
                                              NULL,
                                              VOT_REAL8,
                                              spinArraySize,
                                              spinArrayString);
    if(!xmlParamNode) {
        XLALPrintError("Couldn't create PARAM node: %s\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* return PARAM node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlParamNode;
}


/**
 * \brief Deserializes a \c PulsarSpins array from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts) the \c PulsarSpins array
 * (found as a \c PARAM element in the specified \c RESOURCE element) identified by the given name.
 *
 * \param xmlDocument [in] Pointer to the VOTable XML document containing the array
 * \param resourceType [in] Value of the \c utype attribute of the parent RESOURCE element
 * \param resourceName [in] Unique identifier of the parent RESOURCE element
 * \param paramName [in] Unique identifier of the particular \c PulsarSpins array to be deserialized
 * \param spins [out] Pointer to an empty \c PulsarSpins array to store the deserialized instance
 *
 * \return \c XLAL_SUCCESS if the specified \c PulsarSpins array could be found and
 * deserialized successfully.
 *
 * \sa XLALVOTableXML2PulsarDopplerParamsByName
 * \sa XLALGetSingleVOTableResourceParamAttribute
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4 XLALVOTableDoc2PulsarSpinsByName(const xmlDocPtr xmlDocument,
                                      const char *resourceType,
                                      const char *resourceName,
                                      const char *paramName,
                                      PulsarSpins spins)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALVOTableDoc2PulsarSpinsByName";
    CHAR *nodeContent = NULL;
    CHAR *nodeContentWorker = NULL;
    int arraySize = 0;
    int i;

    /* sanity check */
    if(!xmlDocument) {
        XLALPrintError("Invalid input parameters: xmlDocument\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!resourceType) {
        XLALPrintError("Invalid input parameters: resourceType\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!resourceName) {
        XLALPrintError("Invalid input parameters: resourceName\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!paramName) {
        XLALPrintError("Invalid input parameters: paramName\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!spins) {
        XLALPrintError("Invalid input parameters: spins\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }

    /* retrieve arraysize (number of pulsar spins) */
    nodeContent = (CHAR*)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, resourceType, resourceName, paramName, VOT_ARRAYSIZE);
    if(!nodeContent || sscanf((char*)nodeContent, "%i", &arraySize) == EOF || arraySize == 0) {
        /* clean up*/
        if(nodeContent) xmlFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.%s.arraysize\n", resourceName, paramName);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* retrieve pulsar spin array (string) */
    xmlFree(nodeContent);
    nodeContent = (CHAR *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, resourceType, resourceName, paramName, VOT_VALUE);
    if(!nodeContent) {
        XLALPrintError("Invalid node content encountered: %s.%s\n", resourceName, paramName);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* finally, parse and store individual spin values */
    nodeContentWorker = strtok(nodeContent, " ");
    for(i = 0; i < arraySize; ++i) {
        /* scan current item */
        if(sscanf((char*)nodeContentWorker, "%lf", &spins[i]) == EOF) {
            /* clean up*/
            xmlFree(nodeContent);
            XLALPrintError("Invalid node content encountered: %s.%s[%i]\n", resourceName, paramName, i);
            XLAL_ERROR(logReference, XLAL_EDATA);
        }

        /* advance to next item */
        nodeContentWorker = strtok(NULL, " ");
    }

    /* sanity check */
    if(i != arraySize) {
        /* clean up*/
        xmlFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.%s (found %i of %i items)\n", resourceName, paramName, i, arraySize);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* clean up*/
    xmlFree(nodeContent);

    return XLAL_SUCCESS;
}


/**
 * \brief Serializes a \c PulsarDopplerParams structure into a VOTable XML %node
 *
 * This function takes a \c PulsarDopplerParams structure and serializes it into a VOTable
 * \c RESOURCE %node identified by the given name. The returned \c xmlNode can then be
 * embedded into an existing %node hierarchy or turned into a full VOTable document.
 *
 * \param pdp [in] Pointer to the \c PulsarDopplerParams structure to be serialized
 * \param name [in] Unique identifier of this particular \c PulsarDopplerParams structure instance
 *
 * \return A pointer to a \c xmlNode that holds the VOTable fragment that represents
 * the \c PulsarDopplerParams structure.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * fragment isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \sa XLALCreateVOTableParamNode
 * \sa XLALCreateVOTableResourceNode
 * \sa XLALCreateVOTableDocumentFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr XLALPulsarDopplerParams2VOTableNode(const PulsarDopplerParams *const pdp, const char *name)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALPulsarDopplerParams2VOTableNode";
    CHAR childNameRefTime[NAMESTR_MAXLEN] = {0};
    CHAR childNameOrbit[NAMESTR_MAXLEN] = {0};
    xmlNodePtr xmlParentNode = NULL;
    xmlNodePtr xmlChildNode = NULL;
    xmlNodePtr xmlChildNodeList = NULL;

    CHAR Alpha[REAL8STR_MAXLEN] = {0};
    CHAR Delta[REAL8STR_MAXLEN] = {0};

    /* make sure that the shared library is the same as the
     * library version the code was compiled against */
    LIBXML_TEST_VERSION

    /* check and convert input parameters */
    if(!pdp || LALSnprintf(Alpha, REAL8STR_MAXLEN, "%g", pdp->Alpha) < 0) {
        XLALPrintError("Invalid input parameter: PulsarDopplerParams->Alpha\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!pdp || LALSnprintf(Delta, REAL8STR_MAXLEN, "%g", pdp->Delta) < 0) {
        XLALPrintError("Invalid input parameter: PulsarDopplerParams->Delta\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
    }

    /* set up PARAM node (Alpha) */
    xmlChildNode = XLALCreateVOTableParamNode("Alpha",
                                              "rad",
                                              VOT_REAL8,
                                              NULL,
                                              Alpha);
    if(!xmlChildNode) {
        XLALPrintError("Couldn't create PARAM node: %s.Alpha\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* initialize child node list with first child */
    xmlChildNodeList = xmlChildNode;

    /* set up PARAM node (Delta) */
    xmlChildNode = XLALCreateVOTableParamNode("Delta",
                                              "rad",
                                              VOT_REAL8,
                                              NULL,
                                              Delta);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.Delta\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add child as first sibling to child node list */
    xmlChildNodeList->next = xmlChildNode;

    /* set up PARAM node (fkdot) */
    xmlChildNode = XLALPulsarSpins2VOTableNode(&pdp->fkdot, "fkdot");
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.fkdot\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add child as second sibling to child node list */
    xmlChildNodeList->next->next = xmlChildNode;

    /* compile child name attribute (refTime) */
    if(!strncpy(childNameRefTime, name, NAMESTR_MAXLEN) || !strncat(childNameRefTime, ".refTime", NAMESTR_MAXLEN)) {
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Child attribute preparation failed: PulsarDopplerParams.refTime\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* set up RESOURCE node (refTime)*/
    xmlChildNode = XLALLIGOTimeGPS2VOTableNode(&pdp->refTime, childNameRefTime);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: %s.refTime\n", childNameRefTime);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add child as third sibling to child node list */
    xmlChildNodeList->next->next->next = xmlChildNode;

    /* compile child name attribute (orbit) */
    if(!strncpy(childNameOrbit, name, NAMESTR_MAXLEN) || !strncat(childNameOrbit, ".orbit", NAMESTR_MAXLEN)) {
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Child attribute preparation failed: PulsarDopplerParams.orbit\n");
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* set up RESOURCE node (orbit)*/
    xmlChildNode = XLALBinaryOrbitParams2VOTableNode(pdp->orbit, childNameOrbit);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: %s.orbit\n", childNameOrbit);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* add child as fourth sibling to child node list */
    xmlChildNodeList->next->next->next->next = xmlChildNode;

    /* set up parent RESOURCE node*/
    xmlParentNode = XLALCreateVOTableResourceNode("PulsarDopplerParams", name, xmlChildNodeList);
    if(!xmlParentNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: %s\n", name);
        XLAL_ERROR_NULL(logReference, XLAL_EFAILED);
    }

    /* return RESOURCE node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlParentNode;
}


/**
 * \brief Deserializes a \c PulsarDopplerParams structure from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts)
 * the \c PulsarDopplerParams structure identified by the given name.
 *
 * \param xmlDocument [in] Pointer to the VOTable XML document containing the structure
 * \param name [in] Unique identifier of the particular \c PulsarDopplerParams structure to be deserialized
 * \param pdp [out] Pointer to an empty \c PulsarDopplerParams structure to store the deserialized instance
 *
 * \return \c XLAL_SUCCESS if the specified \c PulsarDopplerParams structure could be found and
 * deserialized successfully.
 *
 * \sa XLALVOTableXML2PulsarDopplerParamsByName
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4 XLALVOTableDoc2PulsarDopplerParamsByName(const xmlDocPtr xmlDocument, const char *name, PulsarDopplerParams *pdp)
{
    /* set up local variables */
    static const CHAR *logReference = "XLALVOTableDoc2PulsarDopplerParamsByName";
    CHAR childNameRefTime[NAMESTR_MAXLEN] = {0};
    CHAR childNameOrbit[NAMESTR_MAXLEN] = {0};
    xmlChar *nodeContent = NULL;

    /* sanity checks */
    if(!xmlDocument) {
        XLALPrintError("Invalid input parameter: xmlDocument\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }
    if(!pdp) {
        XLALPrintError("Invalid input parameter: pdp\n");
        XLAL_ERROR(logReference, XLAL_EINVAL);
    }

    /* compile child name attribute (refTime) */
    if(!strncpy(childNameRefTime, name, NAMESTR_MAXLEN) || !strncat(childNameRefTime, ".refTime", NAMESTR_MAXLEN)) {
        XLALPrintError("Child attribute preparation failed: PulsarDopplerParams.refTime\n");
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* retrieve PulsarDopplerParams.refTime */
    if(XLALVOTableDoc2LIGOTimeGPSByName(xmlDocument, childNameRefTime, &pdp->refTime)) {
        XLALPrintError("Error parsing XML document content: %s.refTime\n", childNameRefTime);
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* retrieve PulsarDopplerParams.Alpha */
    nodeContent = (xmlChar *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, "PulsarDopplerParams", name, "Alpha", VOT_VALUE);

    /* parse and finally store content */
    if(!nodeContent || sscanf((char*)nodeContent, "%lf", &pdp->Alpha) == EOF) {
        /* clean up*/
        if(nodeContent) xmlFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.Alpha\n", name);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* retrieve PulsarDopplerParams.Delta */
    xmlFree(nodeContent);
    nodeContent = (xmlChar *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, "PulsarDopplerParams", name, "Delta", VOT_VALUE);

    /* parse and finally store content */
    if(!nodeContent || sscanf((char*)nodeContent, "%lf", &pdp->Delta) == EOF) {
        /* clean up*/
        if(nodeContent) xmlFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.Delta\n", name);
        XLAL_ERROR(logReference, XLAL_EDATA);
    }

    /* retrieve PulsarDopplerParams.fkdot */
    if(XLALVOTableDoc2PulsarSpinsByName(xmlDocument, "PulsarDopplerParams", name, "fkdot", pdp->fkdot)) {
        /* clean up */
        xmlFree(nodeContent);
        XLALPrintError("Error parsing XML document content: %s.fkdot\n", name);
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* compile child name attribute (orbit) */
    if(!strncpy(childNameOrbit, name, NAMESTR_MAXLEN) || !strncat(childNameOrbit, ".orbit", NAMESTR_MAXLEN)) {
        /* clean up */
        xmlFree(nodeContent);
        XLALPrintError("Child attribute preparation failed: PulsarDopplerParams.orbit\n");
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* retrieve PulsarDopplerParams.orbit */
    if(XLALVOTableDoc2BinaryOrbitParamsByName(xmlDocument, childNameOrbit, pdp->orbit)) {
        /* clean up */
        xmlFree(nodeContent);
        XLALPrintError("Error parsing XML document content: %s.orbit\n", childNameOrbit);
        XLAL_ERROR(logReference, XLAL_EFAILED);
    }

    /* clean up*/
    xmlFree(nodeContent);

    return XLAL_SUCCESS;
}





/**
 * \brief Serializes a \c gsl_vector into a VOTable XML %node
 *
 * This function takes a \c gsl_vector pointer and serializes it into a VOTable
 * \c PARAM %node identified by the given name. The returned \c xmlNode can then be
 * embedded into an existing %node hierarchy.
 *
 * \return A pointer to a \c xmlNode that holds the VOTable fragment that represents
 * the \c gsl_vector. In case of an error, a null-pointer is returned.\n
 *
 * \note All matrix elements are written with maximal precision "%.16g" for a double
 *
 * \b Important: the caller is responsible to free the allocated memory (when the
 * fragment isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \sa XLALCreateVOTableParamNode
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr
XLALgsl_vector2VOTableNode(const gsl_vector *vect,	/**< [in] input gsl_vector to serialize */
			   const CHAR *name,		/**< [in] Unique identifier for this gsl_vector */
			   const CHAR *unitName		/**< [in] optional unit-name (can be NULL) */
			   )
{
  static const CHAR *fn = "XLALgsl_vector2VOTableNode()";

  xmlNodePtr xmlParamNode = NULL;
  CHAR arraySizeStr[INT4STR_MAXLEN] = {0}; 	/* buffer to hold argument to "arraysize" attribute */
  CHAR REAL8Str[REAL8STR_MAXLEN] = {0};		/* buffer to hold string representations of ONE REAL8 element */
  int arrayStrLen;
  CHAR *arrayStr = NULL;			/* holds full list of dim REAL8 strings */
  int i, dim;

  /* sanity checks */
  if ( !vect || !vect->size ) {
    XLALPrintError("%s: Invalid NULL or empty input: vect\n\n", fn);
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }
  if ( !name || strlen(name) <= 0) {
    XLALPrintError("%s: Invalid NULL or empty input parameter: name\n\n", fn);
    XLAL_ERROR_NULL(fn, XLAL_EINVAL);
  }

  /* get input vector dimension */
  dim = vect->size;	/* guaranteed to be >= 1 */

  /* prepare array size attribute */
  if ( LALSnprintf ( arraySizeStr, sizeof(arraySizeStr), "%d", dim ) < 0 ) {
    XLALPrintError("%s: LALSnprintf() failed for arraySizeStr.\n\n", fn );
    XLAL_ERROR_NULL (fn, XLAL_EFAILED);
  }

  /* prepare string to hold array of numbers */
  arrayStrLen = dim * REAL8STR_MAXLEN;
  if ( (arrayStr = XLALCalloc( 1, arrayStrLen )) == NULL ) {
    XLALPrintError ("%s: Failed to XLALCalloc(1, %d).\n\n", fn, arrayStrLen );
    XLAL_ERROR_NULL(fn, XLAL_ENOMEM );
  }

  for ( i = 0; i < dim; i++ )
    {
      /* add vector element [i] to arrayStr */
      if( LALSnprintf ( REAL8Str, sizeof(REAL8Str), "%.16g", gsl_vector_get ( vect, i ) ) < 0 ) {
	XLALFree ( arrayStr );
	XLALPrintError("%s: failed to convert vector element to string: vect[%d]\n", i);
	XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
      }
      strcat ( arrayStr, REAL8Str );

      /* add delimiter (SPACE)*/
      if ( i < dim-1 )
	strcat ( arrayStr, " ");

    } /* for i < dim */

  /* set up PARAM node */
  if ( (xmlParamNode = XLALCreateVOTableParamNode(name, unitName, VOT_REAL8, arraySizeStr, arrayStr )) == NULL ){
    XLALFree ( arrayStr );
    XLALPrintError("%s: XLALCreateVOTableParamNode() failed to create PARAM node: '%s'. xlalErrno = %d\n", fn, name, xlalErrno);
    XLAL_ERROR_NULL (fn, XLAL_EFUNC);
  }

  XLALFree ( arrayStr );

  /* return PARAM node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
  return xmlParamNode;

} /* XLALgsl_vector2VOTableNode() */



/**
 * \brief Deserializes a \c gsl_vector from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts) the \c gsl_vector
 * (found as a \c PARAM element in the specified \c RESOURCE element) identified by the given name.
 *
 * \param xmlDocument [in] Pointer to the VOTable XML document containing the array
 * \param resourceType [in] Value of the \c utype attribute of the parent RESOURCE element
 * \param resourceName [in] Unique identifier of the parent RESOURCE element
 * \param paramName [in] Unique identifier of the particular \c gsl_vector to be deserialized
 * \param unitName [in] Optional unit-name. If != NULL, it MUST agree with the units specified in XML
 *
 * \return \c pointer to deserialized \c gsl_vector, which is allocated by this function.
 *
 * \sa XLALGetSingleVOTableResourceParamAttribute
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
gsl_vector *
XLALVOTableDoc2gsl_vectorByName(const xmlDocPtr xmlDocument,
				const char *resourceType,
				const char *resourceName,
				const char *paramName,
				const CHAR *unitName
				)
{
  static const CHAR *fn = "XLALVOTableDoc2gsl_vectorByName()";
  CHAR *nodeContent = NULL;
  CHAR *nodeContentWorker = NULL;
  int i, dim;
  gsl_vector *vect = NULL;
  double el;

  /* input sanity check */
  if( !xmlDocument ) {
    XLALPrintError ( "%s: Invalid input parameters: xmlDocument\n", fn);
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }
  if( !resourceType ) {
    XLALPrintError ( "%s: Invalid input parameters: resourceType\n", fn);
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }
  if( !resourceName ) {
    XLALPrintError ( "%s; Invalid input parameters: resourceName\n", fn);
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }
  if( !paramName ) {
    XLALPrintError("%s: Invalid input parameters: paramName\n", fn);
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }

  /* retrieve arraysize (number of vector elements) */
  nodeContent = (CHAR*)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, resourceType, resourceName, paramName, VOT_ARRAYSIZE);
  if(!nodeContent || sscanf((char*)nodeContent, "%d", &dim) == EOF || dim == 0) {
    if(nodeContent) xmlFree(nodeContent);
    XLALPrintError("%s: Invalid node content encountered: %s.%s.arraysize\n", fn, resourceName, paramName);
    XLAL_ERROR_NULL (fn, XLAL_EDATA);
  }
  xmlFree(nodeContent);

  /* retrieve and check unitName, if given by caller */
  if ( unitName )
    {
      nodeContent = (CHAR*)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, resourceType, resourceName, paramName, VOT_UNIT );
      if ( !nodeContent ) {
	xmlFree(nodeContent);
	XLALPrintError("%s: failed to find %s.%s.unit, but caller requested '%s'.\n\n", fn, resourceName, paramName, unitName );
	XLAL_ERROR_NULL ( fn, XLAL_EDATA );
      }
      /* check that unit in XML agrees with caller's requirement */
      if ( strcmp ( unitName, nodeContent ) ) {
	XLALPrintError("%s: %s.%s.unit = '%s', but caller requested '%s'. Sorry, I can't do unit conversions.\n",
		       fn, resourceName, paramName, nodeContent, unitName );
	xmlFree(nodeContent);
	XLAL_ERROR_NULL ( fn, XLAL_EDATA );
      }

      xmlFree(nodeContent);
    } /* check unitName */


  /* retrieve pulsar gsl_vector array (string) */
  nodeContent = (CHAR *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, resourceType, resourceName, paramName, VOT_VALUE);

  if(!nodeContent) {
    XLALPrintError("%s: Invalid node content encountered: %s.%s\n", fn, resourceName, paramName);
    XLAL_ERROR_NULL (fn, XLAL_EDATA);
  }

  /* allocate output gsl_vector */
  if ( (vect = gsl_vector_calloc ( dim )) == NULL ) {
    XLALPrintError ("%s: failed to gsl_vector_calloc(%d).\n", fn, dim );
    XLAL_ERROR_NULL (fn, XLAL_ENOMEM);
  }

  /* finally, parse and store individual vector elements */
  if ( (nodeContentWorker = strtok(nodeContent, " ")) == NULL ) {
    gsl_vector_free ( vect );
    xmlFree(nodeContent);
    XLALPrintError ("%s: failed to parse first array element in node: %s.%s.\n",
		    fn, resourceName, paramName );
    XLAL_ERROR_NULL (fn, XLAL_EDATA);
  }
  for(i = 0; i < dim; i++)
    {
      /* scan current item */
      if ( sscanf( nodeContentWorker, "%lf", &el) == EOF) {
	xmlFree(nodeContent);
	gsl_vector_free ( vect );
	XLALPrintError("%s: Invalid node content encountered: %s.%s[%i] = '%s'\n", fn, resourceName, paramName, i, nodeContentWorker );
	XLAL_ERROR_NULL (fn, XLAL_EDATA);
      }

      gsl_vector_set ( vect, i, el );

      /* advance to next item */
      if ( i < dim-1 )
	{
	  if ( (nodeContentWorker = strtok(NULL, " ")) == NULL ) {
	    gsl_vector_free ( vect );
	    xmlFree(nodeContent);
	    XLALPrintError ("%s: failed to parse %d-th array element in node: %s.%s'.\n", fn, i, resourceName, paramName );
	    XLAL_ERROR_NULL (fn, XLAL_EDATA);
	  }
	} /* if not last element yet: advance */

    } /* for i < dim */

  /* clean up*/
  xmlFree(nodeContent);

  return vect;

} /* XLALVOTableDoc2gsl_vectorByName() */





/**
 * \brief Serializes a \c gsl_matrix into a VOTable XML %node
 *
 * This function takes a \c gsl_matrix pointer and serializes it into a VOTable
 * \c PARAM %node identified by the given name. The returned \c xmlNode can then be
 * embedded into an existing %node hierarchy.
 *
 * \return A pointer to a \c xmlNode that holds the VOTable fragment that represents
 * the \c gsl_matrix. In case of an error, a null-pointer is returned.\n
 *
 * \note All matrix elements are written with maximal precision "%.16g" for a double
 *
 * \note We use a VOTable <PARAM> element with arraysize=dim1xdim2. The arraysize-convention
 *  is that the the first index varies fastest, while the last index is the slowest-varying.
 * We therefore write the matrix in *transposed* order, ie. "cols x rows", such that the
 * fastest-varying index is the column-index!
 *
 * \b Important: the caller is responsible to free the allocated memory (when the
 * fragment isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \sa XLALCreateVOTableParamNode
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr
XLALgsl_matrix2VOTableNode(const gsl_matrix *matrix,	/**< [in] input gsl_matrix to serialize */
			   const CHAR *name,		/**< [in] Unique identifier for this gsl_matrix */
			   const CHAR *unitName		/**< [in] optional unit-name (can be NULL) */
			   )
{
  static const CHAR *fn = "XLALgsl_matrix2VOTableNode()";

  xmlNodePtr xmlParamNode = NULL;
  CHAR arraySizeStr[2*INT4STR_MAXLEN+1] = {0}; 	/* buffer to hold argument to "arraysize" attribute */
  CHAR REAL8Str[REAL8STR_MAXLEN] = {0};		/* buffer to hold string representations of ONE REAL8 element */
  int arrayStrLen;
  CHAR *arrayStr = NULL;			/* holds full list of dim1 x dim2 REAL8 strings */

  int row, col, numRows, numCols;

  /* sanity checks */
  if ( !matrix || !matrix->size1 || !matrix->size2 ) {
    XLALPrintError("%s: Invalid NULL or empty/malformed input: matrix\n\n", fn);
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }
  if ( !name || strlen(name) <= 0) {
    XLALPrintError("%s: Invalid NULL or empty input parameter: name\n\n", fn);
    XLAL_ERROR_NULL(fn, XLAL_EINVAL);
  }

  /* get input matrix dimensions */
  numRows = matrix->size1;	/* guaranteed to be >= 1 */
  numCols = matrix->size2;	/* guaranteed to be >= 1 */

  /* prepare array size attribute.
   * As explained in the header, we write the step through
   * the matrix with rows varying slowest, columns fastest, which corresponds to
   * arraysize="cols x rows" !
   */
  if ( LALSnprintf ( arraySizeStr, sizeof(arraySizeStr), "%dx%d", numCols, numRows ) < 0 ) {
    XLALPrintError("%s: LALSnprintf() failed for arraySizeStr (%dx%d).\n\n", fn, numCols, numRows );
    XLAL_ERROR_NULL (fn, XLAL_EFAILED);
  }

  /* prepare string to hold array of numbers */
  arrayStrLen = numRows * numCols * REAL8STR_MAXLEN;
  if ( (arrayStr = XLALCalloc( 1, arrayStrLen )) == NULL ) {
    XLALPrintError ("%s: Failed to XLALCalloc(1, %d).\n\n", fn, arrayStrLen );
    XLAL_ERROR_NULL(fn, XLAL_ENOMEM );
  }

  /* convert all matrix elements to strings and append to arrayStr */
  for ( row = 0; row < numRows; row++ )
    {
      for ( col = 0; col < numCols; col ++ )
	{
	  /* add matrix element [row,col] to arrayStr */
	  if( LALSnprintf ( REAL8Str, sizeof(REAL8Str), "%.16g", gsl_matrix_get ( matrix, row, col ) ) < 0 ) {
	    XLALFree ( arrayStr );
	    XLALPrintError("%s: failed to convert matrix element to string: vect[%d,%d]\n", row, col);
	    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
	  }
	  strcat ( arrayStr, REAL8Str );

	  /* add delimiter (SPACE)*/
	  if ( col < numCols-1 ) {
	    strcat ( arrayStr, " ");
	  }

	} /* for col < numCols */

      /* add row-delimiter (newline!) */
      if ( row < numRows-1 ) {
	strcat ( arrayStr, "; ");
      }

    } /* for rows < numRows */

  /* assemble PARAM node */
  if ( (xmlParamNode = XLALCreateVOTableParamNode(name, unitName, VOT_REAL8, arraySizeStr, arrayStr )) == NULL ){
    XLALFree ( arrayStr );
    XLALPrintError("%s: XLALCreateVOTableParamNode() failed to create PARAM node: '%s'. xlalErrno = %d\n", fn, name, xlalErrno);
    XLAL_ERROR_NULL (fn, XLAL_EFUNC);
  }

  XLALFree ( arrayStr );

  /* return PARAM node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
  return xmlParamNode;

} /* XLALgsl_matrix2VOTableNode() */



/**
 * \brief Deserializes a \c gsl_matrix from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts) the \c gsl_matrix
 * (found as a \c PARAM element in the specified \c RESOURCE element) identified by the given name.
 *
 * \param xmlDocument [in] Pointer to the VOTable XML document containing the array
 * \param resourceType [in] Value of the \c utype attribute of the parent RESOURCE element
 * \param resourceName [in] Unique identifier of the parent RESOURCE element
 * \param paramName [in] Unique identifier of the particular \c gsl_matrix to be deserialized
 * \param unitName [in] Optional unit-name. If != NULL, it MUST agree with the units specified in XML
 *
 * \return \c pointer to deserialized \c gsl_matrix, which is allocated by this function.
 *
 * \sa XLALGetSingleVOTableResourceParamAttribute
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
gsl_matrix *
XLALVOTableDoc2gsl_matrixByName(const xmlDocPtr xmlDocument,
				const char *resourceType,
				const char *resourceName,
				const char *paramName,
				const CHAR *unitName
				)
{
  static const CHAR *fn = "XLALVOTableDoc2gsl_matrixByName()";
  CHAR *nodeContent = NULL;
  CHAR *nodeContentWorker = NULL;
  int row, col, numRows, numCols;
  gsl_matrix *matrix = NULL;
  double el;

  /* input sanity check */
  if( !xmlDocument ) {
    XLALPrintError ( "%s: Invalid input parameters: xmlDocument\n", fn);
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }
  if( !resourceType ) {
    XLALPrintError ( "%s: Invalid input parameters: resourceType\n", fn);
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }
  if( !resourceName ) {
    XLALPrintError ( "%s; Invalid input parameters: resourceName\n", fn);
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }
  if( !paramName ) {
    XLALPrintError("%s: Invalid input parameters: paramName\n", fn);
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }

  /* retrieve arraysize (number of matrixor elements) */
  nodeContent = (CHAR*)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, resourceType, resourceName, paramName, VOT_ARRAYSIZE);
  if(!nodeContent || sscanf((char*)nodeContent, "%dx%d", &numCols, &numRows) == EOF || numRows == 0 || numCols == 0 ) {
    if(nodeContent) xmlFree(nodeContent);
    XLALPrintError("%s: Invalid node content encountered: %s.%s.arraysize\n", fn, resourceName, paramName);
    XLAL_ERROR_NULL (fn, XLAL_EDATA);
  }
  xmlFree(nodeContent);

  /* retrieve and check unitName, if given by caller */
  if ( unitName )
    {
      nodeContent = (CHAR*)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, resourceType, resourceName, paramName, VOT_UNIT );
      if ( !nodeContent ) {
	xmlFree(nodeContent);
	XLALPrintError("%s: failed to find %s.%s.unit, but caller requested '%s'.\n\n", fn, resourceName, paramName, unitName );
	XLAL_ERROR_NULL ( fn, XLAL_EDATA );
      }
      /* check that unit in XML agrees with caller's requirement */
      if ( strcmp ( unitName, nodeContent ) ) {
	XLALPrintError("%s: %s.%s.unit = '%s', but caller requested '%s'. Sorry, I can't do unit conversions.\n",
		       fn, resourceName, paramName, nodeContent, unitName );
	xmlFree(nodeContent);
	XLAL_ERROR_NULL ( fn, XLAL_EDATA );
      }
      xmlFree(nodeContent);
    } /* check unitName */


  /* retrieve pulsar gsl_matrix array (string) */
  nodeContent = (CHAR *)XLALGetSingleVOTableResourceParamAttribute(xmlDocument, resourceType, resourceName, paramName, VOT_VALUE);

  if(!nodeContent) {
    XLALPrintError("%s: Invalid node content encountered: %s.%s\n", fn, resourceName, paramName);
    XLAL_ERROR_NULL (fn, XLAL_EDATA);
  }

  /* allocate output gsl_matrix */
  if ( (matrix = gsl_matrix_calloc ( numRows, numCols )) == NULL ) {
    XLALPrintError ("%s: failed to gsl_matrix_calloc(%d,%d).\n", fn, numRows, numCols );
    XLAL_ERROR_NULL (fn, XLAL_ENOMEM);
  }

  /* finally, parse and store individual matrix elements */
  if ( (nodeContentWorker = strtok(nodeContent, " ;")) == NULL ) {
    gsl_matrix_free ( matrix );
    xmlFree(nodeContent);
    XLALPrintError ("%s: failed to parse first array element in node: %s.%s.\n",
		    fn, resourceName, paramName );
    XLAL_ERROR_NULL (fn, XLAL_EDATA);
  }
  for (row = 0; row < numRows; row++)
    {
      for (col=0; col < numCols; col ++ )
	{
	  /* scan current item */
	  if ( sscanf( nodeContentWorker, "%lf", &el) == EOF) {
	    xmlFree(nodeContent);
	    gsl_matrix_free ( matrix );
	    XLALPrintError("%s: Invalid node content encountered: %s.%s[%d,%d] = '%s'\n", fn, resourceName, paramName, row, col, nodeContentWorker );
	    XLAL_ERROR_NULL (fn, XLAL_EDATA);
	  }

	  gsl_matrix_set ( matrix, row, col, el );

	  /* advance to next item */
	  if ( (row < numRows-1) || (col < numCols-1) )
	    {
	      if ( (nodeContentWorker = strtok(NULL, " ;")) == NULL ) {
		gsl_matrix_free ( matrix );
		xmlFree(nodeContent);
		XLALPrintError ("%s: failed to parse array (%d,%d) element in node: %s.%s'.\n", fn, row, col, resourceName, paramName );
		XLAL_ERROR_NULL (fn, XLAL_EDATA);
	      }
	    } /* not parsed last element yet: advance to next one */

	} /* for col < numCols */

    } /* for row < numRows */

  /* clean up*/
  xmlFree(nodeContent);

  return matrix;

} /* XLALVOTableDoc2gsl_matrixByName() */
