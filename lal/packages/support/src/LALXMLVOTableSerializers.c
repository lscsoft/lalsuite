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
