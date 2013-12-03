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

#include <lal/LALMalloc.h>
#include <lal/XLALError.h>
#include <lal/LALXML.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXMLVOTableSerializers.h>

#define INT4STR_MAXLEN          15
#define REAL8STR_MAXLEN         25
#define PULSARSPINSTR_MAXLEN    4
#define NAMESTR_MAXLEN          256

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
 * \sa XLALCreateVOTParamNode
 * \sa XLALCreateVOTResourceNode
 * \sa XLALCreateVOTDocumentFromTree
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr XLALLIGOTimeGPS2VOTNode(const LIGOTimeGPS *const ltg, const char *name)
{
    /* set up local variables */
    xmlNodePtr xmlParentNode = NULL;
    xmlNodePtr xmlChildNode = NULL;
    xmlNodePtr xmlChildNodeList = NULL;

    CHAR gpsSecondsBuffer[INT4STR_MAXLEN] = {0};
    CHAR gpsNanoSecondsBuffer[INT4STR_MAXLEN] = {0};

    /* check and prepare input parameters */
    if(!ltg || snprintf(gpsSecondsBuffer, INT4STR_MAXLEN, "%i", ltg->gpsSeconds) < 0) {
        XLALPrintError("Invalid input parameter: LIGOTimeGPS->gpsSeconds\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if(!ltg || snprintf(gpsNanoSecondsBuffer, INT4STR_MAXLEN, "%i", ltg->gpsNanoSeconds) < 0) {
        XLALPrintError("Invalid input parameter: LIGOTimeGPS->gpsNanoSeconds\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if(!name || strlen(name) <= 0) {
        XLALPrintError("Invalid input parameter: name\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* set up RESOURCE node child (first PARAM) */
    xmlChildNode = XLALCreateVOTParamNode("gpsSeconds",
                                          "s",
                                          VOT_INT4,
                                          NULL,
                                          gpsSecondsBuffer);
    if(!xmlChildNode) {
        XLALPrintError("Couldn't create PARAM node: , %s.gpsSeconds\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* initialize child node list with first child */
    xmlChildNodeList = xmlChildNode;

    /* set up RESOURCE node child (second PARAM) */
    xmlChildNode = XLALCreateVOTParamNode("gpsNanoSeconds",
                                          "ns",
                                          VOT_INT4,
                                          NULL,
                                          gpsNanoSecondsBuffer);
    if(!xmlChildNode) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.gpsNanoSeconds\n", name);
        XLAL_ERROR_NULL(XLAL_EFAILED);
    }

    /* add child as first sibling to child node list */
    xmlChildNodeList->next = xmlChildNode;

    /* set up RESOURCE node*/
    xmlParentNode = XLALCreateVOTResourceNode("LIGOTimeGPS", name, xmlChildNodeList);
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
 * \brief Deserializes a \c LIGOTimeGPS structure from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts)
 * the \c LIGOTimeGPS structure identified by the given name.
 *
 * \param xmlDocument [in] Pointer to the VOTable XML document containing the structure
 * \param name [in] Unique identifier of the particular \c LIGOTimeGPS structure to be deserialized
 * \param ltg [out] Pointer to an empty \c  LIGOTimeGPS structure to store the deserialized instance
 *
 * \return \c XLAL_SUCCESS if the specified \c LIGOTimeGPS structure could be found and
 * deserialized successfully.
 *
 * \sa XLALVOTXML2LIGOTimeGPSByName
 * \sa XLALGetSingleNodeContentByXPath
 *
 * \author Oliver Bock\n
 * Albert-Einstein-Institute Hannover, Germany
 */
INT4
XLALVOTDoc2LIGOTimeGPSByName ( const xmlDocPtr xmlDocument, const CHAR *name, LIGOTimeGPS *ltg )
{
    /* set up local variables */
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
    if(!ltg) {
        XLALPrintError("Invalid input parameter: ltg\n");
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* retrieve LIGOTimeGPS.gpsSeconds content */
    if ( (nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, name, "gpsSeconds", VOT_PARAM, VOT_VALUE )) == NULL ) {
      XLAL_ERROR ( XLAL_EFUNC );
    }
    /* parse content */
    if( sscanf ( nodeContent, "%i", &ltg->gpsSeconds) == EOF ) {
      XLALFree(nodeContent);
      XLALPrintError("Invalid node content encountered: %s.gpsSeconds\n", name);
      XLAL_ERROR ( XLAL_EDATA );
    }
    XLALFree(nodeContent);

    /* retrieve LIGOTimeGPS.gpsNanoSeconds content */
    nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, name, "gpsNanoSeconds", VOT_PARAM, VOT_VALUE );
    if( !nodeContent || sscanf( nodeContent, "%i", &ltg->gpsNanoSeconds) == EOF) {
      if ( nodeContent ) XLALFree (nodeContent);
      XLALPrintError("Invalid node content encountered: %s.gpsNanoSeconds\n", name);
      XLAL_ERROR(XLAL_EDATA);
    }

    /* clean up*/
    XLALFree ( nodeContent );

    return XLAL_SUCCESS;

} /* XLALVOTDoc2LIGOTimeGPSByName() */


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
 * \sa XLALCreateVOTParamNode
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr
XLALgsl_vector2VOTNode(const gsl_vector *vect,	/**< [in] input gsl_vector to serialize */
                       const CHAR *name,	/**< [in] Unique identifier for this gsl_vector */
                       const CHAR *unitName	/**< [in] optional unit-name (can be NULL) */
                       )
{
  xmlNodePtr xmlParamNode = NULL;
  CHAR arraySizeStr[INT4STR_MAXLEN] = {0}; 	/* buffer to hold argument to "arraysize" attribute */
  CHAR REAL8Str[REAL8STR_MAXLEN] = {0};		/* buffer to hold string representations of ONE REAL8 element */
  int arrayStrLen;
  CHAR *arrayStr = NULL;			/* holds full list of dim REAL8 strings */
  int i, dim;

  /* sanity checks */
  if ( !vect || !vect->size ) {
    XLALPrintError("%s: Invalid NULL or empty input: vect\n\n", __func__);
    XLAL_ERROR_NULL (XLAL_EINVAL);
  }
  if ( !name || strlen(name) <= 0) {
    XLALPrintError("%s: Invalid NULL or empty input parameter: name\n\n", __func__);
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* get input vector dimension */
  dim = vect->size;	/* guaranteed to be >= 1 */

  /* prepare array size attribute */
  if ( snprintf ( arraySizeStr, sizeof(arraySizeStr), "%d", dim ) < 0 ) {
    XLALPrintError("%s: snprintf() failed for arraySizeStr.\n\n", __func__ );
    XLAL_ERROR_NULL (XLAL_EFAILED);
  }

  /* prepare string to hold array of numbers */
  arrayStrLen = dim * REAL8STR_MAXLEN;
  if ( (arrayStr = XLALCalloc( 1, arrayStrLen )) == NULL ) {
    XLALPrintError ("%s: Failed to XLALCalloc(1, %d).\n\n", __func__, arrayStrLen );
    XLAL_ERROR_NULL(XLAL_ENOMEM );
  }

  for ( i = 0; i < dim; i++ )
    {
      /* add vector element [i] to arrayStr */
      if( snprintf ( REAL8Str, sizeof(REAL8Str), "%.16g", gsl_vector_get ( vect, i ) ) < 0 ) {
	XLALFree ( arrayStr );
	XLALPrintError("%s: failed to convert vector element to string: vect[%d]\n", i);
	XLAL_ERROR_NULL ( XLAL_EINVAL );
      }
      strcat ( arrayStr, REAL8Str );

      /* add delimiter (SPACE)*/
      if ( i < dim-1 )
	strcat ( arrayStr, " ");

    } /* for i < dim */

  /* set up PARAM node */
  if ( (xmlParamNode = XLALCreateVOTParamNode(name, unitName, VOT_REAL8, arraySizeStr, arrayStr )) == NULL ){
    XLALFree ( arrayStr );
    XLALPrintError("%s: XLALCreateVOTParamNode() failed to create PARAM node: '%s'. xlalErrno = %d\n", __func__, name, xlalErrno);
    XLAL_ERROR_NULL (XLAL_EFUNC);
  }

  XLALFree ( arrayStr );

  /* return PARAM node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
  return xmlParamNode;

} /* XLALgsl_vector2VOTNode() */



/**
 * \brief Deserializes a \c gsl_vector from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts) the \c gsl_vector
 * (found as a \c PARAM element in the specified \c RESOURCE element) identified by the given name.
 *
 * \return \c pointer to deserialized \c gsl_vector, which is allocated by this function.
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
gsl_vector *
XLALVOTDoc2gsl_vectorByName ( const xmlDocPtr xmlDocument,	/**< [in] Pointer to the VOTable XML document containing the array */
                              const char *resourcePath,		/**< [in] (optional) parent RESOURCE path */
                              const char *paramName,		/**< [in] name of of the particular PARAM \c gsl_vector to be deserialized */
                              const CHAR *unitName		/**< [in] Optional unit-name. If != NULL, it MUST agree with the units specified in XML */
                              )
{
  CHAR *nodeContent = NULL;
  CHAR *nodeContentWorker = NULL;
  int i, dim;
  gsl_vector *vect = NULL;
  double el;

  /* input sanity check */
  if( !xmlDocument ) {
    XLALPrintError ( "%s: Invalid input parameters: xmlDocument\n", __func__);
    XLAL_ERROR_NULL (XLAL_EINVAL);
  }

  if( !paramName ) {
    XLALPrintError("%s: Invalid input parameters: paramName\n", __func__);
    XLAL_ERROR_NULL (XLAL_EINVAL);
  }

  /* retrieve arraysize (number of vector elements) */
  nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, resourcePath, paramName, VOT_PARAM, VOT_ARRAYSIZE);
  if(!nodeContent || sscanf(nodeContent, "%d", &dim) == EOF || dim == 0) {
    if(nodeContent) XLALFree(nodeContent);
    XLALPrintError("%s: Invalid node content encountered: %s.%s.arraysize\n", __func__, resourcePath, paramName);
    XLAL_ERROR_NULL (XLAL_EDATA);
  }
  XLALFree(nodeContent);

  /* retrieve and check unitName, if given by caller */
  if ( unitName )
    {
      if ( (nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, resourcePath, paramName, VOT_PARAM, VOT_UNIT )) == NULL ) {
	XLALPrintError("%s: failed to find %s.%s.unit, but caller requested '%s'.\n\n", __func__, resourcePath, paramName, unitName );
	XLAL_ERROR_NULL ( XLAL_EDATA );
      }
      /* check that unit in XML agrees with caller's requirement */
      if ( strcmp ( unitName, nodeContent ) ) {
	XLALPrintError("%s: %s.%s.unit = '%s', but caller requested '%s'. Sorry, I can't do unit conversions.\n",
		       __func__, resourcePath, paramName, nodeContent, unitName );
	XLALFree (nodeContent);
	XLAL_ERROR_NULL ( XLAL_EDATA );
      }

      XLALFree (nodeContent);
    } /* check unitName */

  /* retrieve pulsar gsl_vector array (string) */
  if ( (nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, resourcePath, paramName, VOT_PARAM, VOT_VALUE )) == NULL ) {
    XLALPrintError("%s: Invalid node content encountered: %s.%s\n", __func__, resourcePath, paramName);
    XLAL_ERROR_NULL (XLAL_EDATA);
  }

  /* allocate output gsl_vector */
  if ( (vect = gsl_vector_calloc ( dim )) == NULL ) {
    XLALPrintError ("%s: failed to gsl_vector_calloc(%d).\n", __func__, dim );
    XLAL_ERROR_NULL (XLAL_ENOMEM);
  }

  /* finally, parse and store individual vector elements */
  if ( (nodeContentWorker = strtok(nodeContent, " ")) == NULL ) {
    gsl_vector_free ( vect );
    XLALFree(nodeContent);
    XLALPrintError ("%s: failed to parse first array element in node: %s.%s.\n", __func__, resourcePath, paramName );
    XLAL_ERROR_NULL (XLAL_EDATA);
  }
  for(i = 0; i < dim; i++)
    {
      /* scan current item */
      if ( sscanf( nodeContentWorker, "%lf", &el) == EOF) {
	XLALFree(nodeContent);
	gsl_vector_free ( vect );
	XLALPrintError("%s: Invalid node content encountered: %s.%s[%i] = '%s'\n", __func__, resourcePath, paramName, i, nodeContentWorker );
	XLAL_ERROR_NULL (XLAL_EDATA);
      }

      gsl_vector_set ( vect, i, el );

      /* advance to next item */
      if ( i < dim-1 )
	{
	  if ( (nodeContentWorker = strtok(NULL, " ")) == NULL ) {
	    gsl_vector_free ( vect );
	    XLALFree(nodeContent);
	    XLALPrintError ("%s: failed to parse %d-th array element in node: %s.%s'.\n", __func__, i, resourcePath, paramName );
	    XLAL_ERROR_NULL (XLAL_EDATA);
	  }
	} /* if not last element yet: advance */

    } /* for i < dim */

  /* clean up*/
  XLALFree(nodeContent);

  return vect;

} /* XLALVOTDoc2gsl_vectorByName() */





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
 * \note We use a VOTable \<PARAM\> element with arraysize=dim1xdim2. The arraysize-convention
 * is that the the first index varies fastest, while the last index is the slowest-varying.
 * We therefore write the matrix in *transposed* order, ie. "cols x rows", such that the
 * fastest-varying index is the column-index!
 *
 * \b Important: the caller is responsible to free the allocated memory (when the
 * fragment isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \sa XLALCreateVOTParamNode
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
xmlNodePtr
XLALgsl_matrix2VOTNode(const gsl_matrix *matrix,	/**< [in] input gsl_matrix to serialize */
                       const CHAR *name,		/**< [in] Unique identifier for this gsl_matrix */
                       const CHAR *unitName		/**< [in] optional unit-name (can be NULL) */
                       )
{
  xmlNodePtr xmlParamNode = NULL;
  CHAR arraySizeStr[2*INT4STR_MAXLEN+1] = {0}; 	/* buffer to hold argument to "arraysize" attribute */
  CHAR REAL8Str[REAL8STR_MAXLEN] = {0};		/* buffer to hold string representations of ONE REAL8 element */
  int arrayStrLen;
  CHAR *arrayStr = NULL;			/* holds full list of dim1 x dim2 REAL8 strings */

  int row, col, numRows, numCols;

  /* sanity checks */
  if ( !matrix || !matrix->size1 || !matrix->size2 ) {
    XLALPrintError("%s: Invalid NULL or empty/malformed input: matrix\n\n", __func__);
    XLAL_ERROR_NULL (XLAL_EINVAL);
  }
  if ( !name || strlen(name) <= 0) {
    XLALPrintError("%s: Invalid NULL or empty input parameter: name\n\n", __func__);
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  /* get input matrix dimensions */
  numRows = matrix->size1;	/* guaranteed to be >= 1 */
  numCols = matrix->size2;	/* guaranteed to be >= 1 */

  /* prepare array size attribute.
   * As explained in the header, we write the step through
   * the matrix with rows varying slowest, columns fastest, which corresponds to
   * arraysize="cols x rows" !
   */
  if ( snprintf ( arraySizeStr, sizeof(arraySizeStr), "%dx%d", numCols, numRows ) < 0 ) {
    XLALPrintError("%s: snprintf() failed for arraySizeStr (%dx%d).\n\n", __func__, numCols, numRows );
    XLAL_ERROR_NULL (XLAL_EFAILED);
  }

  /* prepare string to hold array of numbers */
  arrayStrLen = numRows * numCols * REAL8STR_MAXLEN;
  if ( (arrayStr = XLALCalloc( 1, arrayStrLen )) == NULL ) {
    XLALPrintError ("%s: Failed to XLALCalloc(1, %d).\n\n", __func__, arrayStrLen );
    XLAL_ERROR_NULL(XLAL_ENOMEM );
  }

  /* convert all matrix elements to strings and append to arrayStr */
  for ( row = 0; row < numRows; row++ )
    {
      for ( col = 0; col < numCols; col ++ )
	{
	  /* add matrix element [row,col] to arrayStr */
	  if( snprintf ( REAL8Str, sizeof(REAL8Str), "%.16g", gsl_matrix_get ( matrix, row, col ) ) < 0 ) {
	    XLALFree ( arrayStr );
	    XLALPrintError("%s: failed to convert matrix element to string: vect[%d,%d]\n", row, col);
	    XLAL_ERROR_NULL ( XLAL_EINVAL );
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
  if ( (xmlParamNode = XLALCreateVOTParamNode(name, unitName, VOT_REAL8, arraySizeStr, arrayStr )) == NULL ){
    XLALFree ( arrayStr );
    XLALPrintError("%s: XLALCreateVOTParamNode() failed to create PARAM node: '%s'. xlalErrno = %d\n", __func__, name, xlalErrno);
    XLAL_ERROR_NULL (XLAL_EFUNC);
  }

  XLALFree ( arrayStr );

  /* return PARAM node (needs to be xmlFreeNode'd or xmlFreeDoc'd by caller!!!) */
  return xmlParamNode;

} /* XLALgsl_matrix2VOTNode() */



/**
 * \brief Deserializes a \c gsl_matrix from a VOTable XML document
 *
 * This function takes a VOTable XML document and deserializes (extracts) the \c gsl_matrix
 * (found as a \c PARAM element in the specified \c RESOURCE element) identified by the given name.
 *
 * \return \c pointer to deserialized \c gsl_matrix, which is allocated by this function.
 *
 * \author Reinhard Prix\n
 * Albert-Einstein-Institute Hannover, Germany
 */
gsl_matrix *
XLALVOTDoc2gsl_matrixByName ( const xmlDocPtr xmlDocument,	/**< [in] Pointer to the VOTable XML document containing the array */
                              const char *resourcePath,		/**< [in] (optional) parent RESOURCE path */
                              const char *paramName,		/**< [in] name of of the particular PARAM \c gsl_vector to be deserialized */
                              const CHAR *unitName		/**< [in] Optional unit-name. If != NULL, it MUST agree with the units specified in XML */
                              )
{
  CHAR *nodeContent = NULL;
  CHAR *nodeContentWorker = NULL;
  int row, col, numRows, numCols;
  gsl_matrix *matrix = NULL;
  double el;

  /* input sanity check */
  if( !xmlDocument ) {
    XLALPrintError ( "%s: Invalid input parameters: xmlDocument\n", __func__);
    XLAL_ERROR_NULL (XLAL_EINVAL);
  }
  if( !paramName ) {
    XLALPrintError("%s: Invalid input parameters: paramName\n", __func__);
    XLAL_ERROR_NULL (XLAL_EINVAL);
  }

  /* retrieve arraysize (number of matrixor elements) */
  nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, resourcePath, paramName, VOT_PARAM, VOT_ARRAYSIZE );
  if(!nodeContent || sscanf( nodeContent, "%dx%d", &numCols, &numRows) == EOF || numRows == 0 || numCols == 0 ) {
    if(nodeContent) XLALFree(nodeContent);
    XLALPrintError("%s: Invalid node content encountered: %s.%s.arraysize\n", __func__, resourcePath, paramName);
    XLAL_ERROR_NULL (XLAL_EDATA);
  }
  XLALFree(nodeContent);

  /* retrieve and check unitName, if given by caller */
  if ( unitName )
    {
      if ( (nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, resourcePath, paramName, VOT_PARAM, VOT_UNIT )) == NULL ) {
	XLALPrintError("%s: failed to find %s.%s.unit, but caller requested '%s'.\n\n", __func__, resourcePath, paramName, unitName );
	XLAL_ERROR_NULL ( XLAL_EDATA );
      }
      /* check that unit in XML agrees with caller's requirement */
      if ( strcmp ( unitName, nodeContent ) ) {
	XLALPrintError("%s: %s.%s.unit = '%s', but caller requested '%s'. Sorry, I can't do unit conversions.\n",
		       __func__, resourcePath, paramName, nodeContent, unitName );
	XLALFree(nodeContent);
	XLAL_ERROR_NULL ( XLAL_EDATA );
      }
      XLALFree(nodeContent);
    } /* check unitName */


  /* retrieve pulsar gsl_matrix array (string) */
  if ( (nodeContent = XLALReadVOTAttributeFromNamedElement ( xmlDocument, resourcePath, paramName, VOT_PARAM, VOT_VALUE )) == NULL ) {
    XLALPrintError("%s: Invalid node content encountered: %s.%s\n", __func__, resourcePath, paramName);
    XLAL_ERROR_NULL (XLAL_EDATA);
  }

  /* allocate output gsl_matrix */
  if ( (matrix = gsl_matrix_calloc ( numRows, numCols )) == NULL ) {
    XLALPrintError ("%s: failed to gsl_matrix_calloc(%d,%d).\n", __func__, numRows, numCols );
    XLAL_ERROR_NULL (XLAL_ENOMEM);
  }

  /* finally, parse and store individual matrix elements */
  if ( (nodeContentWorker = strtok(nodeContent, " ;")) == NULL ) {
    gsl_matrix_free ( matrix );
    XLALFree(nodeContent);
    XLALPrintError ("%s: failed to parse first array element in node: %s.%s.\n",
		    __func__, resourcePath, paramName );
    XLAL_ERROR_NULL (XLAL_EDATA);
  }
  for (row = 0; row < numRows; row++)
    {
      for (col=0; col < numCols; col ++ )
	{
	  /* scan current item */
	  if ( sscanf( nodeContentWorker, "%lf", &el) == EOF) {
	    XLALFree(nodeContent);
	    gsl_matrix_free ( matrix );
	    XLALPrintError("%s: Invalid node content encountered: %s.%s[%d,%d] = '%s'\n", __func__, resourcePath, paramName, row, col, nodeContentWorker );
	    XLAL_ERROR_NULL (XLAL_EDATA);
	  }

	  gsl_matrix_set ( matrix, row, col, el );

	  /* advance to next item */
	  if ( (row < numRows-1) || (col < numCols-1) )
	    {
	      if ( (nodeContentWorker = strtok(NULL, " ;")) == NULL ) {
		gsl_matrix_free ( matrix );
		XLALFree(nodeContent);
		XLALPrintError ("%s: failed to parse array (%d,%d) element in node: %s.%s'.\n", __func__, row, col, resourcePath, paramName );
		XLAL_ERROR_NULL (XLAL_EDATA);
	      }
	    } /* not parsed last element yet: advance to next one */

	} /* for col < numCols */

    } /* for row < numRows */

  /* clean up*/
  XLALFree(nodeContent);

  return matrix;

} /* XLALVOTDoc2gsl_matrixByName() */
