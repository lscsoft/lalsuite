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

#include <stdio.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <lal/LALXML.h>
#include <lal/XLALError.h>
#include <lal/LALDatatypes.h>

#define INT4STR_MAXLEN 15

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

xmlDocPtr XLALCreateVOTableXMLFromTree(const xmlNodePtr xmlTree)
{
	/* set up local variables */
	static const char *logReference = "XLALCreateVOTableXMLFromTree";
	xmlDocPtr xmlDoc = NULL;
    xmlNodePtr xmlRootNode = NULL;

	/* make sure that the shared library is the same as the
	 * library version the code was compiled against */
	LIBXML_TEST_VERSION

	/* TODO: XML exception handling */

	/* set up XML document */
    xmlDoc = xmlNewDoc(BAD_CAST("1.0"));

    /* set up root node */
    xmlRootNode = xmlNewNode(NULL, BAD_CAST("VOTABLE"));
    xmlDocSetRootElement(xmlDoc, xmlRootNode);

    /* append tree to root node */
    xmlAddChild(xmlRootNode, xmlTree);

	/* return VOTable document (needs to be xmlFreeDoc'd by caller!!!) */
	return xmlDoc;
}

xmlNodePtr XLALLIGOTimeGPS2VOTableNode(const LIGOTimeGPS *const ltg, const char *name)
{
	/* set up local variables */
	static const char *logReference = "XLALLIGOTimeGPS2VOTableNode";
	char gpsSecondsBuffer[INT4STR_MAXLEN] = {0};
	char gpsNanoSecondsBuffer[INT4STR_MAXLEN] = {0};
    xmlNodePtr xmlResourceNode = NULL;
    xmlNodePtr xmlParamNodeGpsSeconds = NULL;
    xmlNodePtr xmlParamNodeGpsNanoSeconds = NULL;

	/* make sure that the shared library is the same as the
	 * library version the code was compiled against */
	LIBXML_TEST_VERSION

	/* check and prepare input parameters */
	if(!ltg || snprintf(gpsSecondsBuffer, INT4STR_MAXLEN, "%i", ltg->gpsSeconds) < 0) {
		XLALPrintError("Invalid input parameter: LIGOTimeGPS->gpsSeconds\n");
		XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
	}
	if(!ltg || snprintf(gpsNanoSecondsBuffer, INT4STR_MAXLEN, "%i", ltg->gpsNanoSeconds) < 0) {
		XLALPrintError("Invalid input parameter: LIGOTimeGPS->gpsNanoSeconds\n");
		XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
	}
	if(!name || strlen(name) <= 0) {
		XLALPrintError("Invalid input parameter: name\n");
		XLAL_ERROR_NULL(logReference, XLAL_EINVAL);
	}

	/* TODO: XML exception handling */

    /* set up RESOURCE node*/
    xmlResourceNode = xmlNewNode(NULL, BAD_CAST("RESOURCE"));
    xmlNewProp(xmlResourceNode, BAD_CAST("utype"), BAD_CAST("LIGOTimeGPS"));
    xmlNewProp(xmlResourceNode, BAD_CAST("name"), BAD_CAST(name));

    /* set up RESOURCE node child (first PARAM) */
    xmlParamNodeGpsSeconds = xmlNewChild(xmlResourceNode, NULL, BAD_CAST("PARAM"), NULL);
    xmlNewProp(xmlParamNodeGpsSeconds, BAD_CAST("name"), BAD_CAST("gpsSeconds"));
    xmlNewProp(xmlParamNodeGpsSeconds, BAD_CAST("datatype"), BAD_CAST("int"));
    xmlNewProp(xmlParamNodeGpsSeconds, BAD_CAST("unit"), BAD_CAST("s"));
    xmlNewProp(xmlParamNodeGpsSeconds, BAD_CAST("value"), BAD_CAST(gpsSecondsBuffer));

    /* set up RESOURCE node child (second PARAM) */
    xmlParamNodeGpsNanoSeconds = xmlNewChild(xmlResourceNode, NULL, BAD_CAST("PARAM"), NULL);
    xmlNewProp(xmlParamNodeGpsNanoSeconds, BAD_CAST("name"), BAD_CAST("gpsNanoSeconds"));
    xmlNewProp(xmlParamNodeGpsNanoSeconds, BAD_CAST("datatype"), BAD_CAST("int"));
    xmlNewProp(xmlParamNodeGpsNanoSeconds, BAD_CAST("unit"), BAD_CAST("ns"));
    xmlNewProp(xmlParamNodeGpsNanoSeconds, BAD_CAST("value"), BAD_CAST(gpsNanoSecondsBuffer));

    /* return RESOURCE node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlResourceNode;
}



xmlChar * XLALLIGOTimeGPS2VOTableXML(const LIGOTimeGPS *const ltg, const char *name)
{
	/* set up local variables */
	static const char *logReference = "XLALLIGOTimeGPS2VOTableXML";
	xmlChar *xmlStringBuffer = NULL;
	int xmlStringBufferSize = -1;

	/* TODO: XML exception handling */

	/* prepare XML serialization */
	xmlThrDefIndentTreeOutput(1);

	/* build XML document */
	xmlNodePtr xmlTree = XLALLIGOTimeGPS2VOTableNode(ltg, name);
	xmlDocPtr xmlDoc = XLALCreateVOTableXMLFromTree(xmlTree);

    /* dump XML document */
	xmlDocDumpFormatMemoryEnc(xmlDoc, &xmlStringBuffer, &xmlStringBufferSize, "UTF-8", 1);

	/* clean up */
	xmlFreeDoc(xmlDoc);
	xmlCleanupParser();

	/* return XML string (needs to be xmlFree'd by caller!!!) */
	return xmlStringBuffer;
}
