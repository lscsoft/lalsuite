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

/**
 * \file
 * \ingroup XML
 * \brief Implementation of the common VOTable XML API
 */

/* ---------- includes ---------- */
#include <string.h>
#include <stdarg.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

#include <lal/LALStdio.h>
#include <lal/StringInput.h>
#include <lal/XLALError.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXML.h>
#include <lal/LALMalloc.h>
#include <lal/LALString.h>

/* ---------- defines and macros ---------- */
#define VOTABLE_VERSION     "1.1"
#define VOTABLE_NS_PREFIX   "vot"
#define VOTABLE_NS_URL      "http://www.ivoa.net/xml/VOTable/v"VOTABLE_VERSION
#define VOTABLE_SCHEMA      "http://www.ivoa.net/xml/VOTable/v"VOTABLE_VERSION
#define XPATHSTR_MAXLEN     500
#define TRUE 1
#define FALSE 0

/* ---------- internal prototypes ---------- */

char * XMLCleanVOTTableWhitespace ( const char *xmlString );
const char *XLALgetDefaultFmt4Datatype ( VOTABLE_DATATYPE datatype );
const char *XLALVOTprintfFromArray ( VOTABLE_DATATYPE datatype, const char *fmt, void *dataPtr, UINT4 index );
int XLALVOTsscanfToArray ( VOTABLE_DATATYPE datatype, const xmlChar *content, void *arrayPtr, UINT4 index );

CHAR * XLALgetXPathToVOTElementAttibute ( const CHAR *resourcePath, const CHAR *elementName, VOTABLE_ELEMENT elementType, VOTABLE_ATTRIBUTE attribute );
CHAR *XLALVOTResourcePath2XPath ( const CHAR *resourcePath );
void *XLALVOTAllocateArray ( VOTABLE_DATATYPE datatype, UINT4 rowCounter );

/* ---------- function definitions ---------- */

/**
 * \brief Creates a VOTable \c PARAM %node
 *
 * This function creates a VOTable \c PARAM %node with the specified properties.
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
xmlNodePtr
XLALCreateVOTParamNode ( const char *name,		/**< [in] \c name attribute of the \c PARAM %node (mandatory) */
                         const char *unit,		/**< [in] \c unit attribute of the \c PARAM %node (optional) */
                         VOTABLE_DATATYPE datatype,	/**< [in] \c datatype attribute of the \c PARAM %node (mandatory) */
                         const char *arraysize,		/**< [in] \c arraysize attribute of the \c PARAM %node (optional) */
                         const char *value		/**< [in] \c value attribute of the \c PARAM %node (mandatory, empty value allowed) */
                         )
{
    /* set up local variables */
    static const CHAR *fn = "XLALCreateVOTParamNode";
    xmlNodePtr xmlParamNode = NULL;
    static const CHAR *datatypeString;

    /* create node */
    xmlParamNode = xmlNewNode(NULL, CAST_CONST_XMLCHAR("PARAM"));
    if(xmlParamNode == NULL) {
        XLALPrintError("Element instantiation failed: PARAM\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }

    /* add attributes */
    /* mandatory: name */
    if(!name || strlen(name) <= 0) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Missing mandatory attribute: name\n");
        XLAL_ERROR_NULL(fn, XLAL_EINVAL);
    }
    if(!xmlNewProp(xmlParamNode, CAST_CONST_XMLCHAR("name"), CAST_CONST_XMLCHAR(name))) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Attribute instantiation failed: name\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }
    /* optional: unit */
    if(unit && strlen(unit) > 0) {
        if(!xmlNewProp(xmlParamNode, CAST_CONST_XMLCHAR("unit"), CAST_CONST_XMLCHAR(unit))) {
            /* clean up */
            xmlFreeNode(xmlParamNode);
            XLALPrintError("Attribute instantiation failed: unit\n");
            XLAL_ERROR_NULL(fn, XLAL_EFAILED);
        }
    }
    /* mandatory: datatype */
    if ( ( datatypeString = XLALVOTDatatype2String ( datatype )) == NULL ) {
      XLALPrintError ("%s: XLALVOTDatatype2String() failed.\n", fn );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

    if(!xmlNewProp(xmlParamNode, CAST_CONST_XMLCHAR("datatype"), CAST_CONST_XMLCHAR(datatypeString))) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Attribute instantiation failed: datatype\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }
    /* optional: arraysize */
    if(arraysize && strlen(arraysize) > 0) {
        if(!xmlNewProp(xmlParamNode, CAST_CONST_XMLCHAR("arraysize"), CAST_CONST_XMLCHAR(arraysize))) {
            /* clean up */
            xmlFreeNode(xmlParamNode);
            XLALPrintError("Attribute instantiation failed: arraysize\n");
            XLAL_ERROR_NULL(fn, XLAL_EFAILED);
        }
    }
    /* mandatory: value (empty value allowed) */
    if(!value) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Missing mandatory attribute: value\n");
        XLAL_ERROR_NULL(fn, XLAL_EINVAL);
    }
    if(!xmlNewProp(xmlParamNode, CAST_CONST_XMLCHAR("value"), CAST_CONST_XMLCHAR(value))) {
        /* clean up */
        xmlFreeNode(xmlParamNode);
        XLALPrintError("Attribute instantiation failed: value\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }

    /* return PARAM node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlParamNode;

} /* XLALCreateVOTParamNode() */


/**
 * \brief Creates a VOTable \c FIELD %node
 *
 * This function creates a VOTable \c FIELD %node with the specified properties.
 *
 * \return A \c xmlNodePtr that holds the new \c FIELD %node.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * %node isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr
XLALCreateVOTFieldNode ( const char *name,		/**< [in] \c name attribute of the \c FIELD %node (mandatory) */
                         const char *unit,		/**< [in] \c unit attribute of the \c FIELD %node (optional) */
                         VOTABLE_DATATYPE datatype,	/**< [in] \c datatype attribute of the \c FIELD %node (mandatory) */
                         const char *arraysize		/**< [in] \c arraysize attribute of the \c FIELD %node (optional) */
                         )
{
    /* set up local variables */
    static const CHAR *fn = "XLALCreateVOTFieldNode()";
    xmlNodePtr xmlFieldNode = NULL;
    static const CHAR *datatypeString;

    /* create node */
    xmlFieldNode = xmlNewNode(NULL, CAST_CONST_XMLCHAR("FIELD"));
    if(xmlFieldNode == NULL) {
        XLALPrintError("Element instantiation failed: FIELD\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }

    /* add attributes */
    /* mandatory: name */
    if(!name || strlen(name) <= 0) {
        /* clean up */
        xmlFreeNode(xmlFieldNode);
        XLALPrintError("Missing mandatory attribute: name\n");
        XLAL_ERROR_NULL(fn, XLAL_EINVAL);
    }
    if(!xmlNewProp(xmlFieldNode, CAST_CONST_XMLCHAR("name"), CAST_CONST_XMLCHAR(name))) {
        /* clean up */
        xmlFreeNode(xmlFieldNode);
        XLALPrintError("Attribute instantiation failed: name\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }
    /* optional: unit */
    if(unit && strlen(unit) > 0) {
        if(!xmlNewProp(xmlFieldNode, CAST_CONST_XMLCHAR("unit"), CAST_CONST_XMLCHAR(unit))) {
            /* clean up */
            xmlFreeNode(xmlFieldNode);
            XLALPrintError("Attribute instantiation failed: unit\n");
            XLAL_ERROR_NULL(fn, XLAL_EFAILED);
        }
    }
    /* mandatory: datatype */
    if ( ( datatypeString = XLALVOTDatatype2String ( datatype )) == NULL ) {
      XLALPrintError ("%s: XLALVOTDatatype2String() failed.\n", fn );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

    if(!xmlNewProp(xmlFieldNode, CAST_CONST_XMLCHAR("datatype"), CAST_CONST_XMLCHAR(datatypeString))) {
        /* clean up */
        xmlFreeNode(xmlFieldNode);
        XLALPrintError("Attribute instantiation failed: datatype\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }
    /* optional: arraysize */
    if(arraysize && strlen(arraysize) > 0) {
        if(!xmlNewProp(xmlFieldNode, CAST_CONST_XMLCHAR("arraysize"), CAST_CONST_XMLCHAR(arraysize))) {
            /* clean up */
            xmlFreeNode(xmlFieldNode);
            XLALPrintError("Attribute instantiation failed: arraysize\n");
            XLAL_ERROR_NULL(fn, XLAL_EFAILED);
        }
    }

    /* return FIELD node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlFieldNode;

} /* XLALCreateVOTFieldNode() */


/**
 * \brief Creates a VOTable \c RESOURCE %node
 *
 * This function creates a VOTable \c RESOURCE %node with the specified identifier and assigns
 * the given children to it.
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
xmlNodePtr
XLALCreateVOTResourceNode ( const char *type,			/**< [in] Type of the \c RESOURCE %node (typically the \c struct type name) */
                            const char *identifier,		/**< [in] Identifier (name) of the \c RESOURCE %node */
                            const xmlNodePtr childNodeList	/**< [in] Pointer to an array of \c xmlNodes that are to be assigned as children */
                            )
{
    /* set up local variables */
    static const CHAR *fn = "XLALCreateVOTResourceNode";
    xmlNodePtr xmlResourceNode = NULL;
    xmlNodePtr xmlChildNode = childNodeList;

    /* sanity check */
    if(!type) {
        XLALPrintError("Invalid input parameter: type\n");
        XLAL_ERROR_NULL(fn, XLAL_EINVAL);
    }
    if(!identifier) {
        XLALPrintError("Invalid input parameter: identifier\n");
        XLAL_ERROR_NULL(fn, XLAL_EINVAL);
    }

    /* create node */
    xmlResourceNode = xmlNewNode(NULL, CAST_CONST_XMLCHAR("RESOURCE"));
    if(xmlResourceNode == NULL) {
        XLALPrintError("Element instantiation failed: RESOURCE\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }

    /* add attributes */
    if(!xmlNewProp(xmlResourceNode, CAST_CONST_XMLCHAR("name"), CAST_CONST_XMLCHAR(identifier))) {
        /* clean up */
        xmlFreeNode(xmlResourceNode);
        XLALPrintError("Attribute instantiation failed: name\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }

    if(!xmlNewProp(xmlResourceNode, CAST_CONST_XMLCHAR("utype"), CAST_CONST_XMLCHAR(type))) {
        /* clean up */
        xmlFreeNode(xmlResourceNode);
        XLALPrintError("Attribute instantiation failed: utype\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }

    /* add children */
    while(xmlChildNode) {
        if(!xmlAddChild(xmlResourceNode, xmlChildNode)) {
            /* clean up */
            xmlFreeNode(xmlResourceNode);
            XLALPrintError("Couldn't add child node to RESOURCE node!\n");
            XLAL_ERROR_NULL(fn, XLAL_EFAILED);
        }
        /* advance to next sibling in list */
        xmlChildNode = xmlChildNode->next;
    }

    /* return RESOURCE node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlResourceNode;

} /* XLALCreateVOTResourceNode() */


/**
 * \brief Serialize a list of data-arrays into a VOTable \c TABLEDATA %node.
 *
 * Note: the variable-length argument pointers must be *void pointers*, and match in number,
 * and types with the FIELDs in fieldNodeList
 *
 * \return A \c xmlNodePtr that holds the new \c TABLEDATA %node (incl. all children and data).
 * In case of an error, a null-pointer is returned.\n
 *
 * Note: the input nodes 'fieldNodeList' and 'dataContentNode' become part of the returned table node
 * and are *not* copied. Therefore they should NOT be free'ed individually by the caller!
 *
 * \b Important: the caller is responsible to free the allocated memory (when the
 * %node isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr
XLALCreateVOTTableNode ( const char *name,		/**< [in] optional name attribute to assign to this \c TABLE element (may be NULL) */
                         xmlNode *fieldNodeList, 	/**< [in] linked list of \c xmlNodes that are to be assigned as FIELD children */
                         xmlNode *dataContentNode 	/**< [in] pointer to xmlNode to be inserted under the <DATA> element: TABLEDATA, BINARY or FITS */
                         )
{
  static const char *fn = "XLALCreateVOTTableNode()";

    xmlNodePtr xmlTableNode = NULL;
    xmlNodePtr xmlChildNode = NULL;
    int err;

    /* input sanity check */
    if ( !fieldNodeList ) {
      XLALPrintError ("%s: invalid NULL input 'fieldNodeList'\n", fn );
      XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
    }

    /* create master node */
    if ( (xmlTableNode = xmlNewNode(NULL, CAST_CONST_XMLCHAR("TABLE"))) == NULL ) {
      XLALPrintError("%s: Element instantiation failed: TABLE\n", fn);
      err = XLAL_EFAILED;
      goto failed;
    }

    /* add attributes (if any) */
    if(name && !xmlNewProp(xmlTableNode, CAST_CONST_XMLCHAR("name"), CAST_CONST_XMLCHAR(name))) {
      XLALPrintError("%s: Attribute instantiation failed: name\n", fn);
      err = XLAL_EFAILED;
      goto failed;
    }

    /* add FIELD children */
    xmlChildNode = fieldNodeList;	/* init to first field Node */
    while ( xmlChildNode )
      {
        if ( !xmlAddChild ( xmlTableNode, xmlChildNode ) )
          {
            XLALPrintError("%s: Couldn't add child FIELD node to TABLE node!\n", fn);
            err = XLAL_EFAILED;
            goto failed;
          }

        /* advance to next sibling in list */
        xmlChildNode = xmlChildNode->next;

      } /* while xmlChildNode */

    if ( dataContentNode )
      {
        /* create DATA node */
        xmlNodePtr xmlDataNode = NULL;
        if ( (xmlDataNode = xmlNewChild ( xmlTableNode, NULL, CAST_CONST_XMLCHAR("DATA"), NULL )) == NULL ) {
          XLALPrintError ("%s: xmlNewChild() failed to create child node 'DATA' under 'TABLE' node.\n", fn );
          err = XLAL_ENOMEM;
          goto failed;
        }

        if ( !xmlAddChild ( xmlDataNode, dataContentNode ) ) {
          XLALPrintError("%s: failed to add dataContentNode as child to DATA node!\n", fn);
          xmlFreeNode ( xmlDataNode );
          err = XLAL_EFAILED;
          goto failed;
        }
      }

    return xmlTableNode;

 failed:
    xmlFreeNode ( xmlTableNode );
    XLAL_ERROR_NULL ( fn, err );

} /* XLALCreateVOTTableNode() */


/**
 * \brief Serialize a list of data-arrays into a VOTable \c TABLEDATA %node.
 *
 * Note: the variable-length argument pointers must be *void pointers*, and match in number,
 * and type with the list of comma-separated printf-format strings \c fmt
 *
 * \return A \c xmlNodePtr that holds the new \c TABLEDATA %node (incl. all children and data).
 * In case of an error, a null-pointer is returned.\n
 *
 * \b Important: the caller is responsible to free the allocated memory (when the
 * %node isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr
XLALCreateVOTTabledataNode ( const xmlNode *fieldNodeList, 	/**< [in] linked list of FIELD \c xmlNodes (used for datatype information) */
                             UINT4 numRows,			/**< [in] number of *rows* in the table [*must* be <= than the lenght of the data arrays!] */
                             const char *fmt,			/**< [in] optional CSV list of printf-format strings to use for writing data (may be NULL) */
                             ...				/**< [in] list of void-pointers to field column data: must match FIELD datatype specs! */
                             )
{
    static const char *fn = "XLALCreateVOTTabledataNode()";
    int err;

    va_list ap;	/* pointer to each unnamed argument in turn */

    UINT4 row, col, numFields;
    const xmlNode *xmlChildNode = NULL;
    xmlNodePtr xmlTABLEDATAnode = NULL;
    void **dataColumns = NULL;		/* array of void-pointers to variable-length input arguments */
    VOTABLE_DATATYPE *dataTypes = NULL;	/* array of corresponding datatypes, parsed from fieldNodeList */
    const char **dataFmts = NULL;	/* array of format strings for printing */

    /* input sanity check */
    if ( !fieldNodeList ) {
      XLALPrintError ("%s: invalid NULL input 'fieldNodeList'\n", fn );
      XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
    }

    /* count number of FIELD element */
    xmlChildNode = fieldNodeList;	/* init to first field Node */
    numFields = 0;
    while ( xmlChildNode )
      {
        /* advance to next sibling in list */
        xmlChildNode = xmlChildNode->next;
        numFields ++;
      } /* while xmlChildNode */


    /* ---------- prepare column writing from varargs input list ---------- */
    if ( (dataColumns = XLALCalloc ( numFields, sizeof(*dataColumns) )) == NULL ) {
      XLALPrintError ("%s: XLALCalloc ( %d, %d ) failed.\n", fn, numFields, sizeof(*dataColumns) );
      err = XLAL_ENOMEM;
      goto failed;
    }
    if ( (dataTypes = XLALCalloc ( numFields, sizeof(*dataTypes) )) == NULL ) {
      XLALPrintError ("%s: XLALCalloc ( %d, %d ) failed.\n", fn, numFields, sizeof(*dataTypes) );
      err = XLAL_ENOMEM;
      goto failed;
    }
    if ( (dataFmts = XLALCalloc ( numFields, sizeof(char*) )) == NULL ) {
      XLALPrintError ("%s: XLALCalloc ( %d, %d ) failed.\n", fn, numFields, sizeof(*dataTypes) );
      err = XLAL_ENOMEM;
      goto failed;
    }

    /* parse format-string if given */
    TokenList *fmtList = NULL;
    if ( fmt )
      {
        if ( XLALCreateTokenList(&fmtList, fmt, "," ) != XLAL_SUCCESS ) {
          XLALPrintError ("%s: failed to parse format-string '%s'\n", fn, fmt );
          err = XLAL_EINVAL;
          goto failed;
        }
        /* if given, must be consistent with number of FIELDS */
        if ( fmtList->nTokens != numFields ) {
          XLALPrintError ("%s: inconsistent number of FIELDS (%d) and format strings (%d) in '%s'.\n", fn, numFields, fmtList->nTokens, fmt );
          err = XLAL_EINVAL;
          goto failed;
        }
      } /* if fmt given */

    /* handle variable-length input arguments containing the table data (columns) */
    va_start(ap, fmt);

    xmlChildNode = fieldNodeList;	/* init to first field Node */
    /* ----- in a first pass we just catalog the data-pointers and corresponding data-types into arrays */
    for ( col=0; col < numFields; col ++ )	/* loop over all fields (= columns of table) */
      {
        /* get data-pointers */
        dataColumns[col] = va_arg(ap, void *);	/* assemble a list of data-pointers of all data columns */

        /* get data-types */
        char *datatypeStr = NULL;
        if ( (datatypeStr = (char*)xmlGetProp ( (xmlNodePtr)xmlChildNode, CAST_CONST_XMLCHAR("datatype"))) == NULL ) {
          XLALPrintError ("%s: xmlGetProp() failed to find required attribute 'datatype' in FIELD node Nr %d.\n", fn, col );
          err = XLAL_EINVAL;
          goto failed;
        }
        if ( ( dataTypes[col] = XLALVOTString2Datatype ( datatypeStr ) ) == VOT_DATATYPE_LAST ) {
          XLALPrintError ("%s: invalid data-type attribute encountered '%s' in field node Nr %d.\n", fn, datatypeStr, col );
          xmlFree ( datatypeStr );
          err = XLAL_EINVAL;
          goto failed;
        }
        xmlFree ( datatypeStr );

        /* get format strings for printing */
        if ( fmtList )
          dataFmts[col] = fmtList->tokens[col];
        else
          dataFmts[col] = NULL;

        /* advance to next sibling in list */
        xmlChildNode = xmlChildNode->next;

      } /* for col < numFields */
    va_end(ap);

    /* we're ready for assembling the actual TABLEDATA entries now */

    /* create TABLEDATA node */
    if ( ( xmlTABLEDATAnode = xmlNewNode ( NULL, CAST_CONST_XMLCHAR("TABLEDATA") ))== NULL ) {
      XLALPrintError ("%s: xmlNewNode() failed to create 'TABLEDATA' node.\n", fn );
      err = XLAL_ENOMEM;
      goto failed;
    }

    /* ---------- loop over data-arrays and generate each table-row */
    for ( row = 0; row < numRows; row ++ )
      {
        /* create TR node */
        xmlNodePtr xmlThisRowNode = NULL;
        if ( (xmlThisRowNode = xmlNewNode ( NULL, CAST_CONST_XMLCHAR("TR") )) == NULL ) {
          XLALPrintError ("%s: xmlNewNode() failed to create new 'TR' node.\n", fn );
          err = XLAL_EFAILED;
          goto failed;
        }
        if ( xmlAddChild(xmlTABLEDATAnode, xmlThisRowNode ) == NULL ) {
          XLALPrintError ("%s: failed to insert 'TR' node into 'TABLEDATA' node.\n", fn );
          err = XLAL_EFAILED;
          goto failed;
        }

        /* ----- loop over columns and generate each table element */
        for ( col = 0; col < numFields; col ++ )
          {
            /* create TD node */
            xmlNodePtr xmlThisEntryNode = NULL;
            if ( (xmlThisEntryNode = xmlNewNode ( NULL, CAST_CONST_XMLCHAR("TD") )) == NULL ) {
              XLALPrintError ("%s: xmlNewNode() failed to create new 'TD' node.\n", fn );
              err = XLAL_EFAILED;
              goto failed;
            }
            if ( xmlAddChild(xmlThisRowNode, xmlThisEntryNode ) == NULL ) {
              XLALPrintError ("%s: failed to insert 'TD' node into 'TR' node.\n", fn );
              err = XLAL_EFAILED;
              goto failed;
            }

            const char* tmptxt;
            if ( (tmptxt = XLALVOTprintfFromArray ( dataTypes[col], dataFmts[col], dataColumns[col], row )) == NULL ){
              XLALPrintError ("%s: XLALVOTprintfFromArray() failed for row = %d, col = %d. errno = %d.\n", fn, row, col, xlalErrno );
              err = XLAL_EFUNC;
              goto failed;
            }

            xmlNodePtr xmlTextNode;
            if ( (xmlTextNode = xmlNewText (CAST_CONST_XMLCHAR(tmptxt) )) == NULL ) {
              XLALPrintError("%s: xmlNewText() failed to turn text '%s' into node\n", fn, tmptxt );
              err = XLAL_EFAILED;
              goto failed;
            }
            if ( xmlAddChild(xmlThisEntryNode, xmlTextNode ) == NULL ) {
              XLALPrintError ("%s: failed to insert text-node node into 'TD' node.\n", fn );
              err = XLAL_EFAILED;
              goto failed;
            }

          } /* for col < numFields */

      } /* for row < numRows */

    /* free memory */
    if ( dataTypes ) XLALFree ( dataTypes );
    if ( fmtList ) XLALDestroyTokenList( &fmtList );
    if ( dataFmts ) XLALFree ( dataFmts );
    if ( dataColumns ) XLALFree ( dataColumns );

    /* return complete TABLEDATA node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
    return xmlTABLEDATAnode;

 failed:
    if ( dataTypes ) XLALFree ( dataTypes );
    if ( fmtList ) XLALDestroyTokenList( &fmtList );
    if ( dataFmts ) XLALFree ( dataFmts );
    if ( dataColumns ) XLALFree ( dataColumns );
    if ( xmlTABLEDATAnode ) xmlFreeNode ( xmlTABLEDATAnode );

    XLAL_ERROR_NULL ( fn, err );

} /* XLALCreateVOTTabledataNode () */

/** Function to discover and parse the list of FIELD nodes in a given table,
 * returns a vector holding the attributes of all read-in FIELD elements (in order)
 */
VOTFieldVector *
XLALReadVOTFIELDNodes ( const xmlDocPtr xmlDocument,	/**< [in] The XML document to be searched */
                        const CHAR *resourcePath 	/**< [in] hiearchical RESOURCE path 'res1.res2....resN' to table RESOURCE */
                        )
{
  static const char *fn = "XLALReadVOTFIELDNodes()";
  CHAR ePath[XPATHSTR_MAXLEN];
  static const CHAR *pathExt = ".$TABLE.$FIELD";

  xmlNodeSet *fieldNodeSet;
  UINT4 numFields;
  VOTFieldVector *ret;
  enum XLALErrorValue err;

  /* check input consistency */
  if ( !xmlDocument || !resourcePath || !strlen(resourcePath) ) {
    XLALPrintError ("%s: invalid NULL or empty input.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  /* assemble 'extended resource path' leading to FIELD elements */
  if ( snprintf ( ePath, XPATHSTR_MAXLEN, "%s%s", resourcePath, pathExt ) >= XPATHSTR_MAXLEN ) {
    XLALPrintError ("%s: extended resource-path '%s'+'%s' longer than max allowed (%d).\n", resourcePath, pathExt, XPATHSTR_MAXLEN );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  if ( (fieldNodeSet = XLALFindVOTElementsAtPath ( xmlDocument, ePath )) == NULL ) {
    XLALPrintError ( "%s: failed to find list of FIELD elements in table '%s'.\n", fn, resourcePath );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  numFields = fieldNodeSet->nodeNr;
  if ( (ret = XLALCreateVOTFieldVector ( numFields )) == NULL ) {
    err = XLAL_EFUNC;
    goto failed;
  }

  /* ----- read all FIELD elements and stores in output vector ---------- */
  UINT4 i;
  for ( i = 0 ; i < numFields; i ++ )
    {
      xmlNode *thisNode = fieldNodeSet->nodeTab[i];
      VOTField *thisField = &ret->data[i];
      xmlChar *tmp;

      /* ----- read required attributes, fail if not found */
      /* name */
      if ( (thisField->name = xmlGetProp ( thisNode, CAST_CONST_XMLCHAR("name"))) == NULL ) {
        XLALPrintError ("%s: Required attribute 'name' not found in %d'th FIELD node.\n", fn, i );
        err = XLAL_EDATA;
        goto failed;
      }
      /* datatype */
      if ( (tmp = xmlGetProp ( thisNode, CAST_CONST_XMLCHAR("datatype"))) == NULL ) {
        XLALPrintError ("%s: Required attribute 'datatype' not found in %d'th FIELD node (name='%s').\n", fn, i, (char*)thisField->name );
        err = XLAL_EDATA;
        goto failed;
      }
      if ( (thisField->datatype = XLALVOTString2Datatype ( (char*)tmp )) == VOT_DATATYPE_LAST ) {
        XLALPrintError ("%s: Invalid datatype attribute '%s' in %d'th FIELD node (name='%s').\n", fn, (char*)tmp, i, (char*)thisField->name );
        xmlFree ( tmp );
        err = XLAL_EDATA;
        goto failed;
      }
      xmlFree ( tmp );

      /* ----- read optional attributes, ignore if not found */
      thisField->unit = xmlGetProp ( thisNode, CAST_CONST_XMLCHAR("unit") );
      thisField->arraysize = xmlGetProp ( thisNode, CAST_CONST_XMLCHAR("arraysize") );

      /* check for presence of any other allowed VOTable attributes and warn if found, as currently not used */
      const xmlChar *attrib;
      attrib = CAST_CONST_XMLCHAR("ID");
      if ( ( tmp = xmlGetProp ( thisNode, attrib ) ) != NULL ) {
        XLALPrintWarning ( "%s: WARNING, %d'th FIELD node (name='%s') contains attribute '%s' ('%s'), which is currently not used.\n",
                           fn, i, (char*)thisField->name, attrib, (char*)tmp );
        xmlFree ( tmp );
      }
      attrib = CAST_CONST_XMLCHAR("precision");
      if ( ( tmp = xmlGetProp ( thisNode, attrib ) ) != NULL ) {
        XLALPrintWarning ( "%s: WARNING, %d'th FIELD node (name='%s') contains attribute '%s' ('%s'), which is currently not used.\n",
                           fn, i, (char*)thisField->name, attrib, (char*)tmp );
        xmlFree ( tmp );
      }
      attrib = CAST_CONST_XMLCHAR("width");
      if ( ( tmp = xmlGetProp ( thisNode, attrib ) ) != NULL ) {
        XLALPrintWarning ( "%s: WARNING, %d'th FIELD node (name='%s') contains attribute '%s' ('%s'), which is currently not used.\n",
                           fn, i, (char*)thisField->name, attrib, (char*)tmp );
        xmlFree ( tmp );
      }
      attrib = CAST_CONST_XMLCHAR("ref");
      if ( ( tmp = xmlGetProp ( thisNode, attrib ) ) != NULL ) {
        XLALPrintWarning ( "%s: WARNING, %d'th FIELD node (name='%s') contains attribute '%s' ('%s'), which is currently not used.\n",
                           fn, i, (char*)thisField->name, attrib, (char*)tmp );
        xmlFree ( tmp );
      }
      attrib = CAST_CONST_XMLCHAR("ucd");
      if ( ( tmp = xmlGetProp ( thisNode, attrib ) ) != NULL ) {
        XLALPrintWarning ( "%s: WARNING, %d'th FIELD node (name='%s') contains attribute '%s' ('%s'), which is currently not used.\n",
                           fn, i, (char*)thisField->name, attrib, (char*)tmp );
        xmlFree ( tmp );
      }
      attrib = CAST_CONST_XMLCHAR("utype");
      if ( ( tmp = xmlGetProp ( thisNode, attrib ) ) != NULL ) {
        XLALPrintWarning ( "%s: WARNING, %d'th FIELD node (name='%s') contains attribute '%s' ('%s'), which is currently not used.\n",
                           fn, i, (char*)thisField->name, attrib, (char*)tmp );
        xmlFree ( tmp );
      }
      attrib = CAST_CONST_XMLCHAR("type");
      if ( ( tmp = xmlGetProp ( thisNode, attrib ) ) != NULL ) {
        XLALPrintWarning ( "%s: WARNING, %d'th FIELD node (name='%s') contains attribute '%s' ('%s'), which is currently not used.\n",
                           fn, i, (char*)thisField->name, attrib, (char*)tmp );
        xmlFree ( tmp );
      }

    } /* for i < numFields */

  xmlXPathFreeNodeSet ( fieldNodeSet );

  return ret;


 failed:
  if ( fieldNodeSet ) xmlXPathFreeNodeSet ( fieldNodeSet );
  if ( ret ) XLALDestroyVOTFieldVector ( ret );
  XLAL_ERROR_NULL ( fn, err );

} /* XLALReadVOTFIELDNodes() */

/** Function to read a given table column of elementary dataype from a TABLEDATA.
 *
 * Note: this function does not support general array types in the table column.
 * Strings are always handled by returning an array of pointers, irrespective of the
 * arraysize given in the table FIELD.
 *
 * Returns an array (of type corresponding to datatype) containing the numRows values
 * from the given column of the TABLEDATA table.
 *
 * Note: caller needs to XLALFree() the returned array.
 *
 */
void *
XLALReadVOTTabledataSimpleColumn ( const xmlDocPtr xmlDocument, 	/**< [in] The XML document to be searched */
                                   const CHAR *resourcePath, 		/**< [in] hiearchical RESOURCE path 'res1.res2....resN' to table RESOURCE */
                                   UINT4 column,			/**< [in] index of table-column to read */
                                   VOTABLE_DATATYPE datatype,		/**< [in] data-type of column to be read */
                                   UINT4 *numRows			/**< [out] number of successfully read table rows */
                                   )
{
  static const char *fn = "XLALReadVOTTabledataSimpleColumn()";
  CHAR ePath[XPATHSTR_MAXLEN];
  static const CHAR *pathExt = ".$TABLE.$DATA.$TABLEDATA";
  xmlNodeSet *nodeSet;
  const xmlNode *TR, *TD;	/* pointers to walk around in table elements */
  void *ret=NULL;
  UINT4 i, j, nRows;

  /* check input consistency */
  if ( !xmlDocument || !resourcePath || !strlen(resourcePath) || !numRows ) {
    XLALPrintError ("%s: invalid NULL or empty input.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  /* assemble 'extended resource path' leading to TABLEDATA element */
  if ( snprintf ( ePath, XPATHSTR_MAXLEN, "%s%s", resourcePath, pathExt ) >= XPATHSTR_MAXLEN ) {
    XLALPrintError ("%s: extended resource-path '%s'+'%s' longer than max allowed (%d).\n", resourcePath, pathExt, XPATHSTR_MAXLEN );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  if ( (nodeSet = XLALFindVOTElementsAtPath ( xmlDocument, ePath )) == NULL ) {
    XLALPrintError ( "%s: failed to find TABLEDATA element in table '%s'.\n", fn, resourcePath );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  if ( nodeSet->nodeNr != 1 ) {
    XLALPrintError ("%s: more than one TABLEDATA element found in table '%s' ... something went wrong here.\n", fn, resourcePath );
    xmlXPathFreeNodeSet ( nodeSet );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* initialize TR pointer to first child of TABLEDATA node */
  TR = nodeSet->nodeTab[0]->children;

  /* count number of table rows (doesn't modify TR) */
  if ( XLALcountNamedNodes ( TR, CAST_CONST_XMLCHAR("TR"), &nRows ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: failed to count number of table rows (<TR>) in TABLE '%s'\n", fn, resourcePath );
    xmlXPathFreeNodeSet ( nodeSet );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* allocate output array */
  if ( ( ret = XLALVOTAllocateArray ( datatype, nRows ) ) == NULL ) {
    xmlXPathFreeNodeSet ( nodeSet );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* advance to beginning of first table row, ie first <TR> node element */
  if ( (TR = XLALfindNextNamedNode ( TR, CAST_CONST_XMLCHAR("TR") )) == NULL ) {
    XLALPrintError ("%s: Something is wrong with the TABLEDATA structure: no <TR> element found!\n", fn );
    XLALFree(ret);
    xmlXPathFreeNodeSet ( nodeSet );
    XLAL_ERROR_NULL ( fn, XLAL_EDATA );
  }

  for ( i=0; i < nRows; i ++ )
    {
      /* initialize TD pointer to first child of TR node */
      TD = TR->children;
      /* advance to n'th column, ie n'th <TD> node element */

      for ( j=0; j <= column; j ++ )
        {
          /* find next <TD> node element */
          if ( (TD = XLALfindNextNamedNode ( TD, CAST_CONST_XMLCHAR("TD") )) == NULL ) {
            XLALPrintError ("%s: Something is wrong with the TABLEDATA: does not contain %d columns\n", fn, column );
            XLALFree(ret);
            xmlXPathFreeNodeSet ( nodeSet );
            XLAL_ERROR_NULL ( fn, XLAL_EDATA );
          }
          if ( j < column )
            TD = TD->next;	/* step forward to next column element */

        } /* for j <= columns */

      /* TD now points to the right column <TD> node */

      /* parse data-string into output array according to datatype */
      if ( XLALVOTsscanfToArray ( datatype, TD->children->content, ret, i ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: failed to parse row %d.\n", i );
        XLALFree ( ret );
        xmlXPathFreeNodeSet ( nodeSet );
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }

      /* step forward to next row element */
      TR = TR->next;
      TR = XLALfindNextNamedNode ( TR, CAST_CONST_XMLCHAR("TR") );

    } /* while TR elements found */

  /* clean up */
  xmlXPathFreeNodeSet ( nodeSet );

  /* return results */
  (*numRows) = nRows;
  return ret;


} /* XLALReadVOTTabledataSimpleColumn() */




/** Creator function for a VOTFieldVector of given number of fiels
 * Should be free'ed using XLALDestroyVOTFieldVector.
 */
VOTFieldVector *
XLALCreateVOTFieldVector ( UINT4 numFields )
{
  static const char *fn = "XLALCreateVOTFieldVector()";

  VOTFieldVector *ret;

  if ( (ret = XLALMalloc ( sizeof(*ret) )) == NULL ) {
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  ret->length = numFields;
  if ( ( ret->data = XLALCalloc ( numFields, sizeof(ret->data[0]))) == NULL ) {
    XLALFree ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  return ret;

} /* XLALCreateVOTFieldVector() */


/** Destructor function for a VOTFieldVector.
 * Note: this frees all strings stored in the
 * VOTField elements using xmlFree(), so they must have been alloce'd
 * by an libxml2 function such as xmlGetProp(). Don't alloc them using XLALMalloc()!
 */
void
XLALDestroyVOTFieldVector ( VOTFieldVector *vect )
{
  if ( !vect )
    return;

  if ( vect->data )
    {
      UINT4 numFields, i;
      numFields = vect->length;

      for ( i=0; i < numFields; i ++ )
        {
          VOTField *thisField = &vect->data[i];

          if ( thisField->name ) xmlFree ( thisField->name );
          if ( thisField->unit ) xmlFree ( thisField->unit );
          if ( thisField->arraysize ) xmlFree ( thisField->arraysize );

          /* ... add more here if you extend the VOTField struct ! */

        } /* for i < numFields */

      XLALFree ( vect->data );

    } /* if vect->data */

  XLALFree ( vect );

  return;

} /* XLALDestroyVOTFieldVector() */


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
 * \return A pointer to a \c xmlDoc that represents the full VOTable XML document.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * document isn't needed anymore) using \c xmlFreeDoc.
 *
 * Note: the xmlTree passed as input becomes part of the returned xmlDoc, so
 * be careful not to free the xmlTree in addition to xmlDoc (==> double-free!)!!
 *
 * Note2: if the resulting xmlDoc should pass validation, you need to set reconcileNamespace == true,
 * but this is not required to generate a valid XML string.
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlDoc *
XLALCreateVOTDocFromTree ( xmlNodePtr xmlTree,		/**< [in] The XML fragment to be turned into a VOTable document */
                           BOOLEAN reconcileNamespace	/**< [in] switch to turn on namespace-reconciliation, which is required to validate xmlDoc */
                           )
{
    /* set up local variables */
    static const CHAR *fn = "XLALCreateVOTDocFromTree";
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
        XLAL_ERROR_NULL(fn, XLAL_EINVAL);
    }

    /* set up XML document */
    xmlDocument = xmlNewDoc(CAST_CONST_XMLCHAR("1.0"));
    if(xmlDocument == NULL) {
        XLALPrintError("VOTable document instantiation failed\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }

    /* set up root node */
    xmlRootElement = xmlNewNode(NULL, CAST_CONST_XMLCHAR("VOTABLE"));
    if(xmlRootElement == NULL) {
        /* clean up */
        xmlFreeDoc(xmlDocument);
        XLALPrintError("VOTABLE root element instantiation failed\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
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
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
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
    if(!xmlAddChild(xmlRootElement, (xmlNodePtr)xmlTree)) {
        /* clean up */
        xmlFreeDoc(xmlDocument);
        XLALPrintError("Couldn't append given tree to VOTABLE root element\n");
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }

    /* finally, assign root element to document */
    xmlDocSetRootElement(xmlDocument, xmlRootElement);

    if ( reconcileNamespace )
      {
        /* reconcile default namespace with all document elements */
        if(XLALReconcileDefaultNamespace(xmlRootElement, xmlVOTableNamespace) != XLAL_SUCCESS) {
          /* clean up */
          xmlFreeDoc(xmlDocument);
          XLALPrintError("Default namespace reconciliation failed!\n");
          XLAL_ERROR_NULL(fn, XLAL_EFAILED);
        }
      } /* if reconcileNamespace */

    /* return VOTable document (needs to be xmlFreeDoc'd by caller!!!) */
    return xmlDocument;

} /* XLALCreateVOTDocFromTree() */


/**
 * \brief Retrieves a specific attribute of a hiearchically named VOTable element of the form 'res1.res2...resN.leaf'.
 *
 * This function returns the given VOTable attribute-content for the VOTable element specified
 * by a hiearchical name of the form "res1.res2.res3.leaf", where 'resN' are the names of
 * parent RESOURCE elements, and 'leaf' is the name of the element to be read, which can be any one
 * of RESOURCE, PARAM, FIELD, ..
 *
 * \return a string holding the contents of the given leaf's attribute, if found, NULL on error.
 *
 * \b Note: the caller is responsible to free the allocated memory using \c XLALFree()
 *
 * \author Reinhard Prix, Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
CHAR *
XLALReadVOTAttributeFromNamedElement ( const xmlDocPtr xmlDocument,	/**< [in] The XML document to be searched */
                                       const char *resourcePath,	/**< [in] [optional] hiearchical path 'res1.res2....resN' of parent RESOURCES */
                                       const char *elementName,		/**< [in] name of leaf element to parse */
                                       VOTABLE_ELEMENT elementType,	/**< [in] VOtable element type of 'element' */
                                       VOTABLE_ATTRIBUTE attrib		/**< [in] Attribute to be read */
                                       )
{
  static const CHAR *fn = "XLALReadVOTAttributeFromNamedElement()";

  CHAR *xpath;
  static const XML_NAMESPACE xmlVOTableNamespace[1] = {{CAST_CONST_XMLCHAR(VOTABLE_NS_PREFIX), CAST_CONST_XMLCHAR(VOTABLE_NS_URL)}};
  const XML_NAMESPACE_VECTOR xmlNsVector = {xmlVOTableNamespace, 1};
  xmlChar *xmlContent;
  char *ret;

  /* sanity checks */
  if(!xmlDocument || !elementName ) {
    XLALPrintError("%s: Invalid NULL input\n", fn);
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL);
  }

  /* prepare XPath search */
  if ( (xpath = XLALgetXPathToVOTElementAttibute ( resourcePath, elementName, elementType, attrib )) == NULL ) {
    XLALPrintError ("%s: failed to assemble xpath to named leaf-attribute.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* retrieve specified attribute (content) */
  if ( (xmlContent = XLALGetSingleNodeContentByXPath(xmlDocument, xpath, &xmlNsVector)) == NULL ) {
    XLALPrintError ("%s: XLALGetSingleNodeContentByXPath() failed for xpath='%s'\n", fn, xpath );
    XLALFree ( xpath );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  XLALFree ( xpath );

  /* copy attribute-content into a standard CHAR string, so caller can XLALFree() it */
  UINT4 len = strlen ( (const char*)xmlContent ) + 1;
  if ( (ret = XLALMalloc ( len * sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALMalloc (%d)\n", fn, len * sizeof(*ret) );
    xmlFree ( xmlContent );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  strcpy ( ret, (const char*)xmlContent );
  xmlFree ( xmlContent );

  return ret;

} /* XLALReadVOTAttributeFromNamedElement() */


/** Assemble an xpath-expression to a named elements's attribute, allowing for
 * a hiearachical "path" of parent RESOURCES of the form "resource1.resource2....resourceN",
 */
CHAR *
XLALgetXPathToVOTElementAttibute ( const CHAR *resourcePath,	/**< [in] optional path of parent-RESOURCES "res1.res2...resN" */
                                   const CHAR *elementName,	/**< [in] name of 'leaf' element to parse */
                                   VOTABLE_ELEMENT elementType,	/**< [in] type of element */
                                   VOTABLE_ATTRIBUTE attribute	/**< [in] attribute to read out from element */
                                   )
{
  static const char *fn = "XLALgetXPathToVOTElementAttibute()";

  char *xpath;
  char buf[XPATHSTR_MAXLEN];

  /* check input consistency */
  if ( !elementName ) {
    XLALPrintError ("%s: invalid NULL input.\n", fn);
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  xpath = NULL;
  if ( ( xpath = XLALStringAppend( xpath, "/" )) == NULL ) {
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  /* turn parent hierarchy "resourcePath" into an xpath expression for parent RESOURCE elements */
  if ( resourcePath )
    {
      CHAR *resXPath;
      if ( (resXPath = XLALVOTResourcePath2XPath ( resourcePath )) == NULL ) {
        XLALPrintError ("%s: XLALVOTResourcePath2XPath('%s') failed.\n",fn, resourcePath );
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }
      if ( ( xpath = XLALStringAppend( xpath, resXPath )) == NULL ) {
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }
      XLALFree ( resXPath );

    } /* if resource path */

  /* attach last leaf-element of type 'elementType' */
  const char *elementStr;
  if ( (elementStr = XLALVOTElement2String ( elementType )) == NULL ) {
    XLALFree ( xpath );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  if ( snprintf ( buf, XPATHSTR_MAXLEN, "/%s:%s[@name='%s']", VOTABLE_NS_PREFIX, elementStr, elementName ) < 0 ) {
    XLALPrintError ("%s: snprintf() failed to assemble xpath string.\n", fn );
    XLALFree ( xpath );
    XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
  }
  if ( ( xpath = XLALStringAppend( xpath, buf )) == NULL ) {
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* finally attach search path matching the attribute */
  const char *attribName;
  if ( (attribName = XLALVOTAttribute2String ( attribute )) == NULL ) {
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  if ( snprintf ( buf, XPATHSTR_MAXLEN, "/@%s", attribName ) < 0 ) {
    XLALPrintError ("%s: snprintf() failed to assemble xpath string.\n", fn );
    XLALFree ( xpath );
    XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
  }
  if ( ( xpath = XLALStringAppend( xpath, buf )) == NULL ) {
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  return xpath;

} /* XLALgetXPathToVOTElementAttibute() */


/** Generate an xpath statement reflecting a VOT "extended" resource path.
 *
 * Note: the 'extResourcePath' may contain unnamed element types prefixed by '$', eg
 * 'resource1.resource2....resourceN.$TABLE', in order to parse element inside <TABLE> under
 * the named <RESOURCE name="resourceN">, within a named <RESOURCE name="resourceN-1"> etc.
 *
 * NOTE: caller needs to XLALFree() the returned string.
 */
CHAR *
XLALVOTResourcePath2XPath ( const CHAR *extResourcePath )
{
  static const char *fn = "XLALVOTResourcePath2XPath()";
  CHAR *xpath;
  CHAR buf[XPATHSTR_MAXLEN];
  UINT4 i;

  if ( !extResourcePath ) {
    XLALPrintError ("%s: invalid NULL input 'extResourcePath'\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  TokenList *resParents = NULL;
  if ( XLALCreateTokenList(&resParents, extResourcePath, "." ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: Failed to parse hierachical parent resource path '%s'.\n", fn, extResourcePath );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* assemble path of parent RESOURCE elements and unnamed element-types */
  xpath = NULL;
  for ( i=0; i < resParents->nTokens; i ++ )
    {
      if ( resParents->tokens[i][0] == '$' )
        {
          /* refers to an un-named element TYPE */
          if ( snprintf ( buf, XPATHSTR_MAXLEN, "/%s:%s", VOTABLE_NS_PREFIX, resParents->tokens[i] + 1 ) < 0 ) {
            XLALPrintError ("%s: snprintf() failed to assemble xpath string.\n", fn );
            XLALFree ( xpath );
            XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
          }
          if ( ( xpath = XLALStringAppend( xpath, buf )) == NULL ) {
            XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
          }
        } /* if unnamed element type */
      else
        {
          /* refers to a named RESOURCE element */
          if ( snprintf ( buf, XPATHSTR_MAXLEN, "/%s:RESOURCE[@name='%s']", VOTABLE_NS_PREFIX, resParents->tokens[i] ) < 0 ) {
            XLALPrintError ("%s: snprintf() failed to assemble xpath string.\n", fn );
            XLALFree ( xpath );
            XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
          }
          if ( ( xpath = XLALStringAppend( xpath, buf )) == NULL ) {
            XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
          }
        } /* if named RESOURCE element */

    } /* for i < nTokens */

  XLALDestroyTokenList ( &resParents );

  return(xpath);

} /* XLALVOTResourcePath2XPath() */



/** Convert the given VOTable xmlTree into a complete VOTable XML document string
 *
 * This function takes a VOTable XML fragment and returns a full-fledged VOTable XML string.
 * Please note that all restrictions described for \ref XLALCreateVOTDocFromTree also apply here!
 *
 *  This function should be used for the final string-formatting of a VOTable XML-tree.
 *
 * \author Oliver Bock, Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 *
 */
CHAR *
XLALCreateVOTStringFromTree ( xmlNodePtr xmlTree )
{
  const char *fn = "XLALCreateVOTStringFromTree()";
  xmlChar *xmlString;
  xmlDocPtr xmlDoc;
  CHAR *ret;

  /* check input consistency */
  if ( !xmlTree ) {
    XLALPrintError ("%s: invalid NULL input.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  /* ----- Step 1: convert VOTable tree into full-fledged VOTable document, but skip namespace-reconciliation  */
  if ( (xmlDoc = XLALCreateVOTDocFromTree( xmlTree, 0 )) == NULL ) {
    XLALPrintError ("%s: failed to convert input xmlTree into xmlDoc. errno = %s\n", fn, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* ----- Step 2: convert VOTable document into a string, with standard libxml2 indentation */
  if ( (xmlString = XLALXMLDoc2String ( xmlDoc )) == NULL ) {
    XLALPrintError ("%s: failed to convert xmlDoc into string. errno = %s\n", fn, xlalErrno );
    xmlUnlinkNode ( xmlTree );	/* protect input xmlTree from free'ing */
    xmlSetListDoc ( xmlTree, NULL ); /* unset document info */
    xmlFreeDoc ( xmlDoc );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  xmlUnlinkNode ( xmlTree );	/* protect input xmlTree from free'ing */
  xmlSetListDoc ( xmlTree, NULL ); /* unset document info */
  xmlFreeDoc ( xmlDoc );	/* free the document container */

  /* ----- Step 3: post-process string: clean newlines+white-space from table rows */
  if ( ( ret = XMLCleanVOTTableWhitespace ( (const char*)xmlString )) == NULL ) {
    XLALPrintError ("%s: XMLCleanVOTTableWhitespace() failed to clean table whitespace from string. errno=%d\n", fn, xlalErrno );
    xmlFree ( xmlString );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  xmlFree ( xmlString );

  return ret;

} /* XLALCreateVOTStringFromTree() */



/** Returns the string representation of the given VOTABLE_DATATYPE.
 *
 * Note: the returned string is const and MUST not be freed!
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
    case VOT_COMPLEX8:
      datatypeString = "floatComplex";
      break;
    case VOT_COMPLEX16:
      datatypeString = "doubleComplex";
      break;
    default:
      XLALPrintError ("%s: invalid datatype passed (%d), has to be within [1, %d].\n", fn, datatype, VOT_DATATYPE_LAST - 1 );
      XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
      break;
    }

  return datatypeString;

} /* XLALVOTDatatype2String() */


/** Returns the string representation of the given VOTABLE_DATATYPE.
 *
 * Note: the returned string is const and MUST not be freed!
 */
const char*
XLALVOTElement2String ( VOTABLE_ELEMENT element )
{
  static const char *fn = "XLALVOTElement2String()";
  const char *elString = NULL;

  switch(element)
    {
    case VOT_RESOURCE:
      elString = "RESOURCE";
      break;
    case VOT_TABLE:
      elString = "TABLE";
      break;
    case VOT_STREAM:
      elString = "STREAM";
      break;
    case VOT_PARAM:
      elString = "PARAM";
      break;
    case VOT_FIELD:
      elString = "FIELD";
      break;

      /* not complete yet: add more here as needed ... */

    default:
      XLALPrintError ("%s: invalid element passed (%d), has to be within [1, %d].\n", fn, element, VOT_ELEMENT_LAST - 1 );
      XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
      break;
    } /* switch element */

  return elString;

} /* XLALVOTElement2String() */

/** Simply returns the enum VOTABLE_DATATYPE corresponding to the string representation of 'datatype'
 * returns VOT_DATATYPE_LAST if invalid.
 */
VOTABLE_DATATYPE
XLALVOTString2Datatype ( const char *datatypeString )
{
  static const char *fn = "XLALVOTString2Datatype()";

  if ( !datatypeString ) {
    XLALPrintError("%s: invalid NULL input 'datatypeString'.\n", fn );
    return VOT_DATATYPE_LAST;
  }
  else if ( !strcmp ( datatypeString, "boolean" ) )
    return VOT_BOOL;
  else if ( !strcmp ( datatypeString, "bit" ) )
    return VOT_BIT;
  else if ( !strcmp ( datatypeString, "char" ) )
    return VOT_CHAR;
  else if ( !strcmp ( datatypeString, "unicodeChar" ) )
    return VOT_CHAR_UTF;
  else if ( !strcmp ( datatypeString, "unsignedByte" ) )
    return VOT_INT1;
  else if ( !strcmp ( datatypeString, "short" ) )
    return VOT_INT2;
  else if ( !strcmp ( datatypeString, "int" ) )
    return VOT_INT4;
  else if ( !strcmp ( datatypeString, "long" ) )
    return VOT_INT8;
  else if ( !strcmp ( datatypeString, "float" ) )
    return VOT_REAL4;
  else if ( !strcmp ( datatypeString, "double" ) )
    return VOT_REAL8;
  else if ( !strcmp ( datatypeString, "floatComplex" ) )
    return VOT_COMPLEX8;
  else if ( !strcmp ( datatypeString, "doubleComplex" ) )
    return VOT_COMPLEX16;
  else
    {
      XLALPrintError ("%s: invalid datatype string '%s'\n", fn, datatypeString );
      return VOT_DATATYPE_LAST;
    }

} /* XLALVOTString2Datatype() */



/** Simply returns the string representation of the given VOTABLE_ATTRIBUTE
 */
const char*
XLALVOTAttribute2String ( VOTABLE_ATTRIBUTE elementAttribute )
{
  static const char *fn = "XLALVOTAttribute2String()";
  const char *attributeString = NULL;

  switch(elementAttribute)
    {
    case VOT_ID:
      attributeString = "ID";
      break;
    case VOT_UNIT:
      attributeString = "unit";
      break;
    case VOT_DATATYPE:
      attributeString = "datatype";
      break;
    case VOT_PRECISION:
      attributeString = "precision";
      break;
    case VOT_WIDTH:
      attributeString = "width";
      break;
    case VOT_REF:
      attributeString = "ref";
      break;
    case VOT_NAME:
      attributeString = "name";
      break;
    case VOT_UCD:
      attributeString = "ucd";
      break;
    case VOT_UTYPE:
      attributeString = "utype";
      break;
    case VOT_ARRAYSIZE:
      attributeString = "arraysize";
      break;
    case VOT_VALUE:
      attributeString = "value";
      break;
    default:
      XLALPrintError ("%s: invalid paramAttribute (%d), must lie within [1, %d].\n", fn, elementAttribute, VOT_ATTRIBUTE_LAST - 1 );
      XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
    }

  return attributeString;

} /* XLALVOTAttribute2String() */


/** Cleans whitespace in VOTable string representation from table-rows.
 *
 * Note: this function decides on the final formatting of the XML string:
 * we use standard lalxml2 indentation, except for successive <TD></TD> elements,
 * which all stand on one line. Same applies to <TR><TD> and </TD></TR> elements,
 * so that each table-row takes exactly one line only [for better readabiltiy]
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
CHAR *
XMLCleanVOTTableWhitespace ( const char *xmlString )
{
  static const char *fn = "XMLCleanVOTTableWhitespace()";

  /* check input consistency */
  if ( !xmlString ) {
    XLALPrintError ("%s: invalid NULL input.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  /* we want all successive <TD> elements to be on one line,
   * so we step through the whole xml-string and nuke all whitespace
   * found between </TD> and <TD> elements.
   */
  static const char *TD_OPEN  = "<TD>";
  static const char *TD_CLOSE = "</TD>";
  static const char *TR_OPEN  = "<TR>";
  static const char *TR_CLOSE = "</TR>";
  const size_t TD_OPEN_LEN  = strlen ( TD_OPEN );
  const size_t TD_CLOSE_LEN = strlen ( TD_CLOSE );
  const size_t TR_OPEN_LEN  = strlen ( TR_OPEN );
  const size_t TR_CLOSE_LEN = strlen ( TR_CLOSE );

  const char *tmp, *nextEl;
  char *startEl;

  UINT4 numChars = strlen ( xmlString );
  tmp = xmlString;
  /* ----- first pass: replace all whitespace between <TR>..<TD>, </TD>..<TD>, and </TD>...</TR>  by ZERO's ! */
  while( ( startEl = strchr (tmp, '<' ) ) != NULL )
    {
      /* ---------- case 1: <TR> ... <TD> ---------- */
      if ( !strncmp ( startEl, TR_OPEN, TR_OPEN_LEN ) )
        {
          startEl += TR_OPEN_LEN;
          /* find then next XML element */
          if ( (nextEl = strchr (startEl, '<' )) == NULL ) {	/* find next element */
            XLALPrintError ("%s: invalid XML, no XML elements found after <TR>: '%s'\n", startEl );
            XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
          }
          /* check that next element is <TD>, if not ... something is wrong */
          if (strncmp ( nextEl, TD_OPEN, TD_OPEN_LEN ) ) {
            XLALPrintError ("Malformed XML Table: <TR> is not followed by <TD>: '%s'\n", fn, nextEl );
            XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
          }

        } /* if we found <TR> */

      /* ---------- case 2: </TD> ...[<TD>|</TR>] ---------- */
      else if ( !strncmp ( startEl, TD_CLOSE, TD_CLOSE_LEN ) )
        {
          startEl += TD_CLOSE_LEN;
          /* find the next XML element */
          if ( (nextEl = strchr (startEl, '<' )) == NULL ) {	/* find next element */
            XLALPrintError ("%s: invalid XML, no XML elements found after <TR>: '%s'\n", startEl );
            XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
          }

          /* check that next element is either <TD> or </TR>, if not ... something is wrong */
          if ( strncmp ( nextEl, TD_OPEN, TD_OPEN_LEN ) &&
               strncmp ( nextEl, TR_CLOSE, TR_CLOSE_LEN ) )
            {
              XLALPrintError ("%s: Malformed XML Table: </TD> is not followed by <TD> or </TR>: '%s'\n", fn, nextEl );
              XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
            }

        } /* case 2: </TD> ... <TD>|</TR> */
      else /* not a table element <TR> or </TD>, so continue */
        {
          tmp = startEl + 1;	/* skip '<' */
          continue;
        }

      /* now brutally set everything between startEl and nextEl to ZERO */
      memset ( startEl, 0, (size_t)(nextEl - startEl) * sizeof(char) );

      /* skip all of the treated region */
      tmp = nextEl;

    } /* while more XML elements are found */

  char *ret;
  /* ----- second pass: copy all chars skipping all ZEROS */
  if ( (ret = XLALMalloc ( (numChars + 1) * sizeof(char) )) == NULL ) {
    XLALPrintError ("%s: XLALMalloc ( %d ) failed.\n", fn, numChars * sizeof(char));
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  UINT4 l = 0;	/* count chars converted */
  UINT4 i;
  for ( i=0; i < numChars; i ++ )
    {
      if ( xmlString[i] )
        ret[l++] = xmlString[i];
    } /* for i < numChars */

  /* now reduce converted string to its actual length */
  ret = XLALRealloc ( ret, (l+1) * sizeof(char) );
  /* ZERO-terminate final string */
  ret[l] = 0;

  return ret;

} /* XMLCleanVOTTableWhitespace() */


/** Get (static & constant!) format string suitable for writing 'datatype'.
 * This is used for default printf-format in table writing.
 *
 * Note: the returned string pointer is static const and MUST *NOT* be free'ed!
 */
const char *
XLALgetDefaultFmt4Datatype ( VOTABLE_DATATYPE datatype )
{
  static const char *fn = "XLALgetDefaultFmt4Datatype()";

  static const char *fmtString;

  switch(datatype)
    {
    case VOT_BOOL:
      fmtString = "%s";
      break;
    case VOT_BIT:
      fmtString = "%d";
      break;
    case VOT_CHAR:
      fmtString = "%s";
      break;
    case VOT_CHAR_UTF:
      fmtString = "%s";
      break;
    case VOT_INT1:
      fmtString = "%hhd";
      break;
    case VOT_INT2:
      fmtString = "%hd";
      break;
    case VOT_INT4:
      fmtString = "%d";
      break;
    case VOT_INT8:
      fmtString = "%ld";
      break;
    case VOT_REAL4:
      fmtString = "%.7g";
      break;
    case VOT_REAL8:
      fmtString = "%.16g";
      break;
    case VOT_COMPLEX8:
      fmtString = "%.7g %.7g";
      break;
    case VOT_COMPLEX16:
      fmtString = "%.16g %.16g";
      break;
    default:
      XLALPrintError ("%s: invalid datatype passed (%d), has to be within [1, %d].\n", fn, datatype, VOT_DATATYPE_LAST - 1 );
      XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
      break;
    }

  return fmtString;

} /* XLALgetDefaultFmt4Datatype() */


/** Write a VOTable value of given primitive datatype, using element 'index' from array 'dataPtr',
 *  using the (optional) printf format string.
 *
 * NOTE: the returned string is a static pointer and MUST not be freed.
 * The livetime of the resulting string is only until the next call of this function.
 */
const char *
XLALVOTprintfFromArray ( VOTABLE_DATATYPE datatype,	/**< [in] atomic dataypte of element to write */
                         const char *fmt,		/**< [in] format string: if NULL we use default-fmt for datatype */
                         void *dataPtr,			/**< [in] pointer to array of data values */
                         UINT4 index			/**< [in] index of element to write: dataPtr[index] */
                         )
{
  static const char *fn = "XLALVOTprintfFromArray()";

#define TEXTBUFLEN 1024
  static char textbuf[TEXTBUFLEN];
  BOOLEAN val;
  const char *writeFmt;

  if ( fmt )
    writeFmt = fmt;
  else
    {
      if ( (writeFmt = XLALgetDefaultFmt4Datatype ( datatype )) == NULL ) {
        XLALPrintError ("%s: XLALgetDefaultFmt4Datatype() failed for datatype = %d.\n", fn, datatype );
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }
    }

  switch ( datatype )
    {
    case VOT_BOOL:
      val = ((BOOLEAN*)dataPtr)[index];
      if ( snprintf(textbuf, TEXTBUFLEN, writeFmt, val ? "t" : "f" ) < 0) {
        XLALPrintError("%s: failed to convert BOOLEAN element (index=%d) to string using fmt '%s'.\n", fn, index, writeFmt );
        XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      }
      break;
    case VOT_BIT:
      XLALPrintError ("%s: Sorry, datatype 'VOT_BIT' not currently supported!\n", fn );
      XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      break;
    case VOT_CHAR:
      if ( snprintf(textbuf, TEXTBUFLEN, writeFmt, ((CHAR**)dataPtr)[index] ) < 0) {
        XLALPrintError("%s: failed to convert CHAR element (index=%d) to string using fmt '%s'.\n", fn, index, writeFmt );
        XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      }
      break;
    case VOT_CHAR_UTF:
      XLALPrintError ("%s: Sorry, datatype 'VOT_CHAR_UTF' not currently supported!\n", fn );
      XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      break;
    case VOT_INT1:
      if ( snprintf(textbuf, TEXTBUFLEN, writeFmt, ((CHAR*)dataPtr)[index] ) < 0) {
        XLALPrintError("%s: failed to convert INT1 element (index=%d) to string using fmt '%s'.\n", fn, index, writeFmt );
        XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      }
      break;
    case VOT_INT2:
      if ( snprintf(textbuf, TEXTBUFLEN, writeFmt, ((INT2*)dataPtr)[index] ) < 0) {
        XLALPrintError("%s: failed to convert INT2 element (index=%d) to string using fmt '%s'.\n", fn, index, writeFmt );
        XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      }
      break;
    case VOT_INT4:
      if ( snprintf(textbuf, TEXTBUFLEN, writeFmt, ((INT4*)dataPtr)[index] ) < 0) {
        XLALPrintError("%s: failed to convert INT4 element (index=%d) to string using fmt '%s'.\n", fn, index, writeFmt );
        XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      }
      break;
    case VOT_INT8:
      if ( snprintf(textbuf, TEXTBUFLEN, writeFmt, ((INT8*)dataPtr)[index] ) < 0) {
        XLALPrintError("%s: failed to convert INT8 element (index=%d) to string using fmt '%s'.\n", fn, index, writeFmt );
        XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      }
      break;
    case VOT_REAL4:
      if ( snprintf(textbuf, TEXTBUFLEN, writeFmt, ((REAL4*)dataPtr)[index] ) < 0) {
        XLALPrintError("%s: failed to convert REAL4 element (index=%d) to string using fmt '%s'.\n", fn, index, writeFmt );
        XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      }
      break;
    case VOT_REAL8:
      if ( snprintf(textbuf, TEXTBUFLEN, writeFmt, ((REAL8*)dataPtr)[index] ) < 0) {
        XLALPrintError("%s: failed to convert REAL8 element (index=%d) to string using fmt '%s'.\n", fn, index, writeFmt );
        XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      }
      break;
    case VOT_COMPLEX8:
      if ( snprintf(textbuf, TEXTBUFLEN, writeFmt, ((COMPLEX8*)dataPtr)[index].re, ((COMPLEX8*)dataPtr)[index].im ) < 0) {
        XLALPrintError("%s: failed to convert COMPLEX8 element (index=%d) to string using fmt '%s'.\n", fn, index, writeFmt );
        XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      }
      break;
    case VOT_COMPLEX16:
      if ( snprintf(textbuf, TEXTBUFLEN, writeFmt, ((COMPLEX16*)dataPtr)[index].re, ((COMPLEX16*)dataPtr)[index].im ) < 0) {
        XLALPrintError("%s: failed to convert COMPLEX16 element (index=%d) to string using fmt '%s'.\n", fn, index, writeFmt );
        XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      }
      break;
    default:
      XLALPrintError ("%s: invalid datatype (%d), has to be within [1, %d].\n", fn, datatype, VOT_DATATYPE_LAST - 1 );
      XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      break;
    } /* switch datatype */

  return textbuf;

} /* XLALVOTprintfFromArray() */


/** Parse a string into the n'th entry in a given array, according to the given datatype.
 *
 * This function is the 'scanf' analog to XLALVOTprintfFromArray() and is used to parse
 * table entries into an output array.
 *
 * Note: the caller is responible to make sure arrayPtr points to an array of the right
 * datatype and size! It is recommended to use XLALVOTAllocateArray() to allocate such an array.
 */
int
XLALVOTsscanfToArray ( VOTABLE_DATATYPE datatype,	/**< [in] atomic dataypte of element to write */
                       const xmlChar *content,		/**< [in] content string to be parsed */
                       void *arrayPtr,			/**< [out] pointer to target array to be written into */
                       UINT4 index			/**< [in] index of element to write: arrayPtr[index] */
                       )
{
  static const char *fn = "XLALVOTsscanfToArray()";

  UINT4 len;
  const char *instring = (const char*)content;

  switch ( datatype )
    {
    case VOT_BOOL:
      /* allowed TRUE representations are 't', '1' or 'true', ignoring case */
      if ( !xmlStrcasecmp(content, CAST_CONST_XMLCHAR("t")) ||
           !xmlStrcasecmp(content, CAST_CONST_XMLCHAR("1")) ||
           !xmlStrcasecmp(content, CAST_CONST_XMLCHAR("true")) )
        ((BOOLEAN*)arrayPtr)[index] = TRUE;
      /* allowed FALSE representations are 'f', '0' or 'false', ignoring case */
      else if ( !xmlStrcasecmp(content, CAST_CONST_XMLCHAR("f")) ||
                !xmlStrcasecmp(content, CAST_CONST_XMLCHAR("0")) ||
                !xmlStrcasecmp(content, CAST_CONST_XMLCHAR("false")) )
        ((BOOLEAN*)arrayPtr)[index] = FALSE;
      else
        {
          XLALPrintError ("%s: invalid boolean value encountered '%s' for index=%d\n", fn, instring, index );
          XLAL_ERROR ( fn, XLAL_EDATA );
        }
      break;

    case VOT_BIT:
      XLALPrintError ("%s: Sorry, datatype 'VOT_BIT' not currently supported!\n", fn );
      XLAL_ERROR ( fn, XLAL_EFAILED );
      break;

    case VOT_CHAR:
      len = strlen ( instring );
      if ( (((CHAR**)arrayPtr)[index] = XLALMalloc ( len + 1 )) == NULL ) {
        XLALPrintError ("%s: failed to XLALMalloc(%d).\n", fn, len + 1 );
        XLAL_ERROR ( fn, XLAL_ENOMEM );
      }
      strcpy ( ((CHAR**)arrayPtr)[index], instring );
      break;

    case VOT_CHAR_UTF:
      XLALPrintError ("%s: Sorry, datatype 'VOT_CHAR_UTF' not currently supported!\n", fn );
      XLAL_ERROR ( fn, XLAL_EFAILED );
      break;

    case VOT_INT1:
      if ( sscanf( instring, "%hhd", &( ((signed char*)arrayPtr)[index] ) ) != 1 ) {
        XLALPrintError("%s: failed to parse INT1 element '%s' at index=%d.\n", fn, instring, index );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }
      break;

    case VOT_INT2:
      if ( sscanf( instring, "%" LAL_INT2_FORMAT, &( ((INT2*)arrayPtr)[index] ) ) != 1 ) {
        XLALPrintError("%s: failed to parse INT2 element '%s' at index=%d.\n", fn, instring, index );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }
      break;

    case VOT_INT4:
      if ( sscanf( instring, "%" LAL_INT4_FORMAT, &( ((INT4*)arrayPtr)[index] ) ) != 1 ) {
        XLALPrintError("%s: failed to parse INT4 element '%s' at index=%d.\n", fn, instring, index );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }
      break;

    case VOT_INT8:
      if ( sscanf( instring, "%" LAL_INT8_FORMAT, &( ((INT8*)arrayPtr)[index] ) ) != 1 ) {
        XLALPrintError("%s: failed to parse INT4 element '%s' at index=%d.\n", fn, instring, index );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }
      break;

    case VOT_REAL4:
      if ( sscanf( instring, "%" LAL_REAL4_FORMAT, &( ((REAL4*)arrayPtr)[index] ) ) != 1 ) {
        XLALPrintError("%s: failed to parse REAL4 element '%s' at index=%d.\n", fn, instring, index );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }
      break;

    case VOT_REAL8:
      if ( sscanf( instring, "%" LAL_REAL8_FORMAT, &( ((REAL8*)arrayPtr)[index] ) ) != 1 ) {
        XLALPrintError("%s: failed to parse REAL8 element '%s' at index=%d.\n", fn, instring, index );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }
      break;

    case VOT_COMPLEX8:
      if ( sscanf( instring, "%" LAL_REAL4_FORMAT "%" LAL_REAL4_FORMAT , &( ((COMPLEX8*)arrayPtr)[index].re ), &( ((COMPLEX8*)arrayPtr)[index].im ) ) != 2 ) {
        XLALPrintError("%s: failed to parse COMPLEX8 element '%s' at index=%d.\n", fn, instring, index );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }
      break;

    case VOT_COMPLEX16:
      if ( sscanf( instring, "%" LAL_REAL8_FORMAT "%" LAL_REAL8_FORMAT, &( ((COMPLEX16*)arrayPtr)[index].re ), &( ((COMPLEX16*)arrayPtr)[index].im ) ) != 2 ) {
        XLALPrintError("%s: failed to parse COMPLEX16 element '%s' at index=%d.\n", fn, instring, index );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }
      break;

    default:
      XLALPrintError ("%s: invalid datatype (%d), has to be within [1, %d].\n", fn, datatype, VOT_DATATYPE_LAST - 1 );
      XLAL_ERROR ( fn, XLAL_EFAILED );
      break;
    } /* switch datatype */

  return XLAL_SUCCESS;

} /* XLALVOTsscanfToArray() */






/**
 * \brief Discover the list of elements at a given "extended resourcePath".
 *
 * This function searches the given XML document for VOTable element nodes
 * at a given extended resource-path, prescribing the parent-hierarchy
 * of RESOURCE elements and unnamed element-types.
 *
 * \return A "node set" containing all matching xmlNodes at the given "path".
 *
 * Note: the 'extResourcePath' may contain unnamed element types prefixed by '$', eg
 * 'resource1.resource2....resourceN.$TABLE', in order to parse element inside <TABLE> under
 * the named <RESOURCE name="resourceN">, within a named <RESOURCE name="resourceN-1"> etc.
 *
 * \b Important: the caller is responsible to free the returned node-set using xmlXPathFreeNodeSet()
 *
 * \author Reinhard Prix, Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodeSet *
XLALFindVOTElementsAtPath ( const xmlDocPtr xmlDocument,	/**< [in] xmlDocument to search */
                            const CHAR *extResourcePath		/**< [in] "extended" path 'res1.res2....resN.$TABLE.$DATA...' of parent nodes */
                            )
{
  /* set up local variables */
  static const CHAR *fn = "XLALFindVOTElementsAtPath()";
  xmlXPathContextPtr xpathCtx = NULL;
  CHAR *xpath = NULL;
  xmlXPathObjectPtr xpathObj = NULL;
  xmlNodeSetPtr xmlNodes = NULL;

  UINT4 i;
  xmlChar const *xmlNsPrefix = NULL;
  xmlChar const *xmlNsUrl = NULL;

  static const XML_NAMESPACE xmlVOTableNamespace[1] = {{CAST_CONST_XMLCHAR(VOTABLE_NS_PREFIX), CAST_CONST_XMLCHAR(VOTABLE_NS_URL)}};
  const XML_NAMESPACE_VECTOR xmlNsVector = {xmlVOTableNamespace, 1};

  /* sanity checks */
  if ( !xmlDocument || !extResourcePath ) {
    XLALPrintError("%s: Invalid NULL input parameter\n", fn);
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL);
  }

  /* prepare xpath context */
  if ( (xpathCtx = xmlXPathNewContext ( xmlDocument )) == NULL ) {;
    XLALPrintError("%s: XPath context instantiation xmlXPathNewContext() failed\n", fn);
    XLAL_ERROR_NULL ( fn, XLAL_EFAILED);
  }

  /* register namespaces */
  for(i = 0; i < xmlNsVector.count; ++i)
    {
      xmlNsPrefix = xmlNsVector.items[i].prefix;
      xmlNsUrl= xmlNsVector.items[i].url;
      if ( xmlXPathRegisterNs ( xpathCtx, xmlNsPrefix, xmlNsUrl ) )
        {
          xmlXPathFreeContext(xpathCtx);
          XLALPrintError("%s: XPath namespace registration failed: %s=%s\n", fn, xmlNsPrefix, xmlNsUrl);
          XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
        }
    } /* for i < NSvector->count */

  xpath = NULL;
  if ( ( xpath = XLALStringAppend( xpath, "/" )) == NULL ) {
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  /* turn parent hierarchy "extResourcePath" into an xpath expression for parent elements */
  char *resXPath;
  if ( (resXPath = XLALVOTResourcePath2XPath ( extResourcePath )) == NULL ) {
    XLALPrintError ("%s: XLALVOTResourcePath2XPath('%s') failed.\n",fn, extResourcePath );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  if ( ( xpath = XLALStringAppend( xpath, resXPath )) == NULL ) {
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  XLALFree ( resXPath );

  /* run xpath query */
  if ( (xpathObj = xmlXPathEvalExpression ( (xmlChar*)xpath, xpathCtx)) == NULL ) {
    XLALFree ( xpath );
    xmlXPathFreeContext ( xpathCtx );
    XLALPrintError("%s: XPath evaluation failed\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
  }

  /* retrieve node set returned by xpath query */
  xmlNodes = xpathObj->nodesetval;

  /* clean up */
  xpathObj-> nodesetval = NULL;		/* protect node-set from free'ing */
  xmlXPathFreeObject ( xpathObj );
  xmlXPathFreeContext ( xpathCtx );
  XLALFree ( xpath );

  /* return list of nodes that matched query (to be free'ed by xmlXPathFreeNodeSet() */
  return xmlNodes;

} /* XLALFindVOTElementsAtPath() */


/** Allocate an array of given datatype and number of elements
 */
void *
XLALVOTAllocateArray ( VOTABLE_DATATYPE datatype, UINT4 length )
{
  static const char *fn = "XLALVOTAllocateArray()";
  void *ret;
  size_t typelen;

  switch ( datatype )
    {
    case VOT_BOOL:
      typelen = sizeof(BOOLEAN);
      break;
    case VOT_BIT:
      XLALPrintError ("%s: Sorry, datatype 'VOT_BIT' not currently supported!\n", fn );
      XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      break;
    case VOT_CHAR:
      typelen = sizeof(CHAR*);
      break;
    case VOT_CHAR_UTF:
      XLALPrintError ("%s: Sorry, datatype 'VOT_CHAR_UTF' not currently supported!\n", fn );
      XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      break;
    case VOT_INT1:
      typelen = sizeof(CHAR);
      break;
    case VOT_INT2:
      typelen = sizeof(INT2);
      break;
    case VOT_INT4:
      typelen = sizeof(INT4);
      break;
    case VOT_INT8:
      typelen = sizeof(INT8);
      break;
    case VOT_REAL4:
      typelen = sizeof(REAL4);
      break;
    case VOT_REAL8:
      typelen = sizeof(REAL8);
      break;
    case VOT_COMPLEX8:
      typelen = sizeof(COMPLEX8);
      break;
    case VOT_COMPLEX16:
      typelen = sizeof(COMPLEX16);
      break;
    default:
      XLALPrintError ("%s: invalid datatype (%d), has to be within [1, %d].\n", fn, datatype, VOT_DATATYPE_LAST - 1 );
      XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
      break;
    } /* switch datatype */

  if ( (ret = LALMalloc ( length * typelen )) == NULL ) {
    XLALPrintError ( "%s: LALMalloc ( %s ) failed.\n", fn, length * typelen );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  return ret;

} /* XLALVOTAllocateArray() */
