/*
 *  Copyright (C) 2009 Oliver Bock
 *  Copyright (C) 2009 Reinhard Prix
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

#include <string.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

#include <lal/XLALError.h>
#include <lal/PulsarDataTypes.h>
#include <lal/LALXML.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXMLVOTableSerializers.h>
#include <lal/LALXMLVOTableSerializersPulsar.h>

#define INT4STR_MAXLEN          15
#define REAL8STR_MAXLEN         25
#define PULSARSPINSTR_MAXLEN    4
#define NAMESTR_MAXLEN          256



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
 * \sa XLALCreateVOTParamNode
 * \sa XLALCreateVOTResourceNode
 * \sa XLALCreateVOTDocumentFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr
XLALBinaryOrbitParams2VOTNode ( const BinaryOrbitParams *const bop, const char *name)
{
    /* set up local variables */
    xmlNodePtr xmlParentNode = NULL;
    xmlNodePtr xmlChildNode = NULL;
    xmlNodePtr xmlChildNodeList = NULL;

    CHAR argp[REAL8STR_MAXLEN] = {0};
    CHAR asini[REAL8STR_MAXLEN] = {0};
    CHAR ecc[REAL8STR_MAXLEN] = {0};
    CHAR period[REAL8STR_MAXLEN] = {0};

    /* check and prepare input parameters */
    if(!bop || snprintf(argp, REAL8STR_MAXLEN, "%g", bop->argp) < 0) {
        XLALPrintError("Invalid input parameter: BinaryOrbitParams->argp\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if(!bop || snprintf(asini, REAL8STR_MAXLEN, "%g", bop->asini) < 0) {
        XLALPrintError("Invalid input parameter: BinaryOrbitParams->asini\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if(!bop || snprintf(ecc, REAL8STR_MAXLEN, "%g", bop->ecc) < 0) {
        XLALPrintError("Invalid input parameter: BinaryOrbitParams->ecc\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if(!bop || snprintf(period, REAL8STR_MAXLEN, "%g", bop->period) < 0) {
        XLALPrintError("Invalid input parameter: BinaryOrbitParams->period\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* set up PARAM node (argp) */
    xmlChildNode = XLALCreateVOTParamNode("argp",
                                          "rad",
                                          VOT_REAL8,
                                          NULL,
                                          argp);
    if(!xmlChildNode) {
        XLALPrintError("Couldn't create PARAM node: %s.argp\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* initialize child node list with first node */
    xmlChildNodeList = xmlChildNode;

    /* set up PARAM node (asini) */
    xmlChildNode = XLALCreateVOTParamNode("asini",
                                          "s",
                                          VOT_REAL8,
                                          NULL,
                                          asini);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.asini\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* add child as first sibling to child node list */
    xmlChildNodeList->next = xmlChildNode;

    /* set up PARAM node (ecc) */
    xmlChildNode = XLALCreateVOTParamNode("ecc",
                                          NULL,
                                          VOT_REAL8,
                                          NULL,
                                          ecc);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.ecc\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* add child as second sibling to child node list */
    xmlChildNodeList->next->next = xmlChildNode;

    /* set up PARAM node (period) */
    xmlChildNode = XLALCreateVOTParamNode("period",
                                          "s",
                                          VOT_REAL8,
                                          NULL,
                                          period);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.period\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* add child as third sibling to child node list */
    xmlChildNodeList->next->next->next = xmlChildNode;

     /* set up RESOURCE node (tp)*/
    xmlChildNode = XLALLIGOTimeGPS2VOTNode(&bop->tp, "tp");
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: tp\n");
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* add child as fourth sibling to child node list */
    xmlChildNodeList->next->next->next->next = xmlChildNode;

    /* set up parent RESOURCE node*/
    xmlParentNode = XLALCreateVOTResourceNode("BinaryOrbitParams", name, xmlChildNodeList);
    if(!xmlParentNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: %s\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* return RESOURCE node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlParentNode;

} /* XLALBinaryOrbitParams2VOTNode() */


/**
 * \brief Deserializes a \c BinaryOrbitParams structure from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts)
 * the \c BinaryOrbitParams structure identified by the given name.
 *
 * \param xmlDocument [in] Pointer to the VOTable XML document containing the structure
 * \param name [in] Unique identifier of the particular \c BinaryOrbitParams structure to be deserialized
 * \param bop [out] Pointer to an empty \c  BinaryOrbitParams structure to store the deserialized instance
 *
 * \return \c XLAL_SUCCESS if the specified \c BinaryOrbitParams structure could be found and
 * deserialized successfully.
 *
 * \sa XLALVOTXML2BinaryOrbitParamsByName
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4
XLALVOTDoc2BinaryOrbitParamsByName ( const xmlDocPtr xmlDocument, const char *name, BinaryOrbitParams *bop)
{
    /* set up local variables */
    CHAR childName[NAMESTR_MAXLEN];
    CHAR *nodeContent = NULL;

    /* sanity checks */
    if(!xmlDocument) {
        XLALPrintError("Invalid input parameter: xmlDocument\n");
        XLAL_ERROR(XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR(XLAL_EINVAL);
    }
    if(!bop) {
        XLALPrintError("Invalid input parameter: bop\n");
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* compile child name attribute (tp) */
    if ( snprintf ( childName, NAMESTR_MAXLEN, "%s.%s", name, "tp") >= NAMESTR_MAXLEN ) {
      XLALPrintError("Child attribute preparation failed: BinaryOrbitParams.tp\n");
      XLAL_ERROR(XLAL_EFAILED);
    }

    /* retrieve BinaryOrbitParams.tp */
    if ( XLALVOTDoc2LIGOTimeGPSByName(xmlDocument, childName, &bop->tp)) {
      XLALPrintError("Error parsing XML document content: %s.tp\n", childName);
      XLAL_ERROR(XLAL_EFAILED);
    }

    /* retrieve BinaryOrbitParams.argp content */
    nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, name, "argp", VOT_PARAM, VOT_VALUE );
    if( !nodeContent || sscanf( nodeContent, "%lf", &bop->argp) == EOF) {
      if ( nodeContent ) XLALFree (nodeContent);
      XLALPrintError("Invalid node content encountered: %s.argp\n", name);
      XLAL_ERROR(XLAL_EDATA);
    }
    XLALFree (nodeContent);

    /* retrieve BinaryOrbitParams.asini content */
    nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, name, "asini", VOT_PARAM, VOT_VALUE );
    if( !nodeContent || sscanf( nodeContent, "%lf", &bop->asini) == EOF) {
      if ( nodeContent ) XLALFree (nodeContent);
      XLALPrintError("Invalid node content encountered: %s.asini\n", name);
      XLAL_ERROR(XLAL_EDATA);
    }
    XLALFree (nodeContent);

    /* retrieve PulsarDopplerParams.ecc content */
    nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, name, "ecc", VOT_PARAM, VOT_VALUE );
    if( !nodeContent || sscanf( nodeContent, "%lf", &bop->ecc) == EOF) {
      if ( nodeContent ) XLALFree (nodeContent);
      XLALPrintError("Invalid node content encountered: %s.ecc\n", name);
      XLAL_ERROR(XLAL_EDATA);
    }
    XLALFree (nodeContent);

    /* retrieve PulsarDopplerParams.period content */
    nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, name, "period", VOT_PARAM, VOT_VALUE );
    if(!nodeContent || sscanf( nodeContent, "%lf", &bop->period) == EOF ) {
      if(nodeContent) XLALFree (nodeContent);
      XLALPrintError("Invalid node content encountered: %s.period\n", name);
      XLAL_ERROR(XLAL_EDATA);
    }
    XLALFree(nodeContent);

    return XLAL_SUCCESS;

} /* XLALVOTDoc2BinaryOrbitParamsByName() */


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
 * \sa XLALCreateVOTParamNode
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr XLALPulsarSpins2VOTNode(const PulsarSpins *const spins, const char *name)
{
    /* set up local variables */
    xmlNodePtr xmlParamNode = NULL;
    int i;

    CHAR spinArraySize[PULSARSPINSTR_MAXLEN] = {0};
    CHAR spinArrayString[PULSAR_MAX_SPINS*REAL8STR_MAXLEN] = {0};
    CHAR spinArrayStringItem[REAL8STR_MAXLEN] = {0};

    /* sanity checks */
    if(!spins || !(*spins)) {
        XLALPrintError("Invalid input parameter: spins\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if(sizeof(*spins)/sizeof(REAL8) != PULSAR_MAX_SPINS) {
        XLALPrintError("Invalid input parameter: spins (actual size %i differs from defined size %i)\n", sizeof(*spins)/sizeof(REAL8), PULSAR_MAX_SPINS);
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* parse input array */
    for(i = 0; i < PULSAR_MAX_SPINS; ++i) {
        if(snprintf(spinArrayStringItem, REAL8STR_MAXLEN, "%g", (*spins)[i]) < 0) {
            XLALPrintError("Invalid input parameter: spins[%i]\n", i);
            XLAL_ERROR_NULL(XLAL_EINVAL);
        }
        if(!strncat(spinArrayString, spinArrayStringItem, REAL8STR_MAXLEN)) {
            XLALPrintError("Couldn't serialize parameter: spins[%i]\n", i);
            XLAL_ERROR_NULL(XLAL_EFAILED);
        }
        /* add delimiter (SPACE)*/
        if(i<PULSAR_MAX_SPINS-1 && !strncat(spinArrayString, " ", 1)) {
            XLALPrintError("Couldn't add delimiter to parameter: spins[%i]\n", i);
            XLAL_ERROR_NULL(XLAL_EFAILED);
        }
    }

    /* set array size attribute */
    if(snprintf(spinArraySize, PULSARSPINSTR_MAXLEN, "%i", PULSAR_MAX_SPINS) < 0) {
        XLALPrintError("Couldn't prepare attribute: arraysize\n");
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* set up PARAM node */
    xmlParamNode = XLALCreateVOTParamNode(name,
                                          NULL,
                                          VOT_REAL8,
                                          spinArraySize,
                                          spinArrayString);
    if(!xmlParamNode) {
        XLALPrintError("Couldn't create PARAM node: %s\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* return PARAM node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlParamNode;
}


/**
 * \brief Deserializes a \c PulsarSpins array from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts) the \c PulsarSpins array
 * (found as a \c PARAM element with the specified \c RESOURCE parent-path) identified by the given name.
 *
 * \return \c XLAL_SUCCESS if the specified \c PulsarSpins array could be found and
 * deserialized successfully.
 *
 * \sa XLALVOTXML2PulsarDopplerParamsByName
 * \sa XLALGetSingleVOTResourceParamAttribute
 *
 * \author Oliver Bock, Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4
XLALVOTDoc2PulsarSpinsByName ( const xmlDocPtr xmlDocument,	/**< [in] Pointer to the VOTable XML document containing the array */
                               const char *resourcePath,	/**< [in] (optional) parent RESOURCE path */
                               const char *paramName,		/**< [in] Name of the PulsarSpins PARAM element */
                               PulsarSpins spins		/**< [out] Pointer to an empty \c PulsarSpins array to store the deserialized instance */
                               )
{
    /* set up local variables */
    CHAR *nodeContent = NULL;
    CHAR *nodeContentWorker = NULL;
    int arraySize = 0;
    int i;

    /* sanity check */
    if(!xmlDocument) {
        XLALPrintError("Invalid input parameters: xmlDocument\n");
        XLAL_ERROR(XLAL_EINVAL);
    }
    if(!paramName) {
        XLALPrintError("Invalid input parameters: paramName\n");
        XLAL_ERROR(XLAL_EINVAL);
    }
    if(!spins) {
        XLALPrintError("Invalid input parameters: spins\n");
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* retrieve arraysize (number of pulsar spins) */
    nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, resourcePath, paramName, VOT_PARAM, VOT_ARRAYSIZE );
    if(!nodeContent || sscanf( nodeContent, "%i", &arraySize) == EOF || arraySize == 0) {
      if(nodeContent) XLALFree(nodeContent);
      XLALPrintError("Invalid node content encountered: %s.%s.arraysize\n", resourcePath, paramName);
      XLAL_ERROR(XLAL_EDATA);
    }
    XLALFree(nodeContent);

    /* retrieve pulsar spin array (string) */
    nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, resourcePath, paramName, VOT_PARAM, VOT_VALUE );
    if(!nodeContent) {
      XLALPrintError("Invalid node content encountered: %s.%s\n", resourcePath, paramName);
      XLAL_ERROR(XLAL_EDATA);
    }

    /* finally, parse and store individual spin values */
    nodeContentWorker = strtok(nodeContent, " ");
    for(i = 0; i < arraySize; i++ )
      {
        /* scan current item */
        if(sscanf((char*)nodeContentWorker, "%lf", &spins[i]) == EOF) {
          XLALFree(nodeContent);
          XLALPrintError("Invalid node content encountered: %s.%s[%i]\n", resourcePath, paramName, i);
          XLAL_ERROR(XLAL_EDATA);
        }

        /* advance to next item */
        nodeContentWorker = strtok(NULL, " ");
      } /* for i < arraySize */

    /* sanity check */
    if(i != arraySize) {
      XLALFree(nodeContent);
      XLALPrintError("Invalid node content encountered: %s.%s (found %i of %i items)\n", resourcePath, paramName, i, arraySize);
      XLAL_ERROR(XLAL_EDATA);
    }

    /* clean up*/
    XLALFree(nodeContent);

    return XLAL_SUCCESS;

} /* XLALVOTDoc2PulsarSpinsByName() */


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
 * \sa XLALCreateVOTParamNode
 * \sa XLALCreateVOTResourceNode
 * \sa XLALCreateVOTDocumentFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr XLALPulsarDopplerParams2VOTNode(const PulsarDopplerParams *const pdp, const char *name)
{
    /* set up local variables */
    xmlNodePtr xmlParentNode = NULL;
    xmlNodePtr xmlChildNode = NULL;
    xmlNodePtr xmlChildNodeList = NULL;

    CHAR Alpha[REAL8STR_MAXLEN] = {0};
    CHAR Delta[REAL8STR_MAXLEN] = {0};

    /* check and convert input parameters */
    if(!pdp || snprintf(Alpha, REAL8STR_MAXLEN, "%g", pdp->Alpha) < 0) {
        XLALPrintError("Invalid input parameter: PulsarDopplerParams->Alpha\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if(!pdp || snprintf(Delta, REAL8STR_MAXLEN, "%g", pdp->Delta) < 0) {
        XLALPrintError("Invalid input parameter: PulsarDopplerParams->Delta\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* set up PARAM node (Alpha) */
    xmlChildNode = XLALCreateVOTParamNode("Alpha",
                                          "rad",
                                          VOT_REAL8,
                                          NULL,
                                          Alpha);
    if(!xmlChildNode) {
        XLALPrintError("Couldn't create PARAM node: %s.Alpha\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* initialize child node list with first child */
    xmlChildNodeList = xmlChildNode;

    /* set up PARAM node (Delta) */
    xmlChildNode = XLALCreateVOTParamNode("Delta",
                                          "rad",
                                          VOT_REAL8,
                                          NULL,
                                          Delta);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.Delta\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* add child as first sibling to child node list */
    xmlChildNodeList->next = xmlChildNode;

    /* set up PARAM node (fkdot) */
    xmlChildNode = XLALPulsarSpins2VOTNode(&pdp->fkdot, "fkdot");
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.fkdot\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* add child as second sibling to child node list */
    xmlChildNodeList->next->next = xmlChildNode;

    /* set up RESOURCE node (refTime)*/
    xmlChildNode = XLALLIGOTimeGPS2VOTNode(&pdp->refTime, "refTime" );
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: %s.refTime\n", name );
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* add child as third sibling to child node list */
    xmlChildNodeList->next->next->next = xmlChildNode;


    /* set up RESOURCE node (orbit)*/
    xmlChildNode = XLALBinaryOrbitParams2VOTNode(pdp->orbit, "orbit" );
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: %s.orbit\n", name );
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* add child as fourth sibling to child node list */
    xmlChildNodeList->next->next->next->next = xmlChildNode;

    /* set up parent RESOURCE node*/
    xmlParentNode = XLALCreateVOTResourceNode("PulsarDopplerParams", name, xmlChildNodeList);
    if(!xmlParentNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create RESOURCE node: %s\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
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
 * \sa XLALVOTXML2PulsarDopplerParamsByName
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4
XLALVOTDoc2PulsarDopplerParamsByName ( const xmlDocPtr xmlDocument,
                                       const char *name,
                                       PulsarDopplerParams *pdp
                                       )
{
    /* set up local variables */
    CHAR childNameRefTime[NAMESTR_MAXLEN] = {0};
    CHAR childNameOrbit[NAMESTR_MAXLEN] = {0};
    CHAR *nodeContent = NULL;

    /* sanity checks */
    if(!xmlDocument) {
        XLALPrintError("Invalid input parameter: xmlDocument\n");
        XLAL_ERROR(XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR(XLAL_EINVAL);
    }
    if(!pdp) {
        XLALPrintError("Invalid input parameter: pdp\n");
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* compile child name attribute (refTime) */
    if(!strncpy(childNameRefTime, name, NAMESTR_MAXLEN) || !strncat(childNameRefTime, ".refTime", NAMESTR_MAXLEN)) {
        XLALPrintError("Child attribute preparation failed: PulsarDopplerParams.refTime\n");
        XLAL_ERROR(XLAL_EFAILED);
    }

    /* retrieve PulsarDopplerParams.refTime */
    if(XLALVOTDoc2LIGOTimeGPSByName(xmlDocument, childNameRefTime, &pdp->refTime)) {
        XLALPrintError("Error parsing XML document content: %s.refTime\n", childNameRefTime);
        XLAL_ERROR(XLAL_EFAILED);
    }

    /* retrieve PulsarDopplerParams.Alpha */
    nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, name, "Alpha", VOT_PARAM, VOT_VALUE );
    if(!nodeContent || sscanf( nodeContent, "%lf", &pdp->Alpha) == EOF) {
      if(nodeContent) XLALFree(nodeContent);
      XLALPrintError("Invalid node content encountered: %s.Alpha\n", name);
      XLAL_ERROR(XLAL_EDATA);
    }
    XLALFree ( nodeContent );

    /* retrieve PulsarDopplerParams.Delta */
    nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, name, "Delta", VOT_PARAM, VOT_VALUE );
    if(!nodeContent || sscanf( nodeContent, "%lf", &pdp->Delta) == EOF) {
        if(nodeContent) XLALFree(nodeContent);
        XLALPrintError("Invalid node content encountered: %s.Delta\n", name);
        XLAL_ERROR(XLAL_EDATA);
    }
    XLALFree ( nodeContent );

    /* retrieve PulsarDopplerParams.fkdot */
    if ( XLALVOTDoc2PulsarSpinsByName(xmlDocument, name, "fkdot", pdp->fkdot)) {
      XLALPrintError("Error parsing XML document content: %s.fkdot\n", name);
      XLAL_ERROR(XLAL_EFAILED);
    }

    /* compile child name attribute (orbit) */
    if(!strncpy(childNameOrbit, name, NAMESTR_MAXLEN) || !strncat(childNameOrbit, ".orbit", NAMESTR_MAXLEN)) {
      XLALPrintError("Child attribute preparation failed: PulsarDopplerParams.orbit\n");
      XLAL_ERROR(XLAL_EFAILED);
    }
    /* retrieve PulsarDopplerParams.orbit */
    if(XLALVOTDoc2BinaryOrbitParamsByName(xmlDocument, childNameOrbit, pdp->orbit)) {
        /* clean up */
        XLALFree(nodeContent);
        XLALPrintError("Error parsing XML document content: %s.orbit\n", childNameOrbit);
        XLAL_ERROR(XLAL_EFAILED);
    }

    return XLAL_SUCCESS;

} /* XLALVOTDoc2PulsarDopplerParamsByName() */
