/*
 *  Copyright (C) 2011 John Veitch
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
 * \brief Implementation of the VOTable serializers XML LALInference API
 */

#include <string.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

#include <lal/XLALError.h>

#include <lal/LALXML.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXMLVOTableSerializers.h>
#include <lal/LALInferenceXML.h>

#define INT4STR_MAXLEN          15
#define REAL8STR_MAXLEN         25
#define NAMESTR_MAXLEN          256

/**
 * \brief Serializes an array of \c LALInferenceVariables into a VOTable XML %node
 *
 * This function takes a \c LALInferenceVariables structure and serializes it into a VOTable
 * \c RESOURCE %node identified by the given name. The returned \c xmlNode can then be
 * embedded into an existing %node hierarchy or turned into a full VOTable document.
 * A VOTable Table element is returned, with fixed variables as PARAMs and the varying ones as FIELDs.
 * 
 * \param varsArray [in] Pointer to an array of \c LALInferenceVariables structures to be serialized
 * \param N [in] Number of items in the array
 * 
 * \return A pointer to a \c xmlNode that holds the VOTable fragment that represents
 * the \c LALInferenceVariables array.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * fragment isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \sa XLALCreateVOTParamNode
 * \sa XLALCreateVOTResourceNode
 *
 * \author John Veitch\n
 * 
 */



xmlNodePtr XLALInferenceVariablesArray2VOTTable(const LALInferenceVariables **varsArray, UINT4 N)
{
  xmlNode *fieldNodeList=NULL;
  xmlNode *paramNodeList=NULL;
  xmlNode *dataContentNode=NULL;
  xmlNode *VOTtableNode=NULL;
  
  /* Build a list of PARAM and FIELD elements */
  
  /* Build array of DATA */
  
  /* Create TABLEDATA node */
  
  /* Create a TABLE from the FIELDs and TABLEDATA nodes */
  
  VOTtableNode= XLALCreateVOTTableNode (
			 const char *name,		/**< [in] optional name attribute to assign to this \c TABLE element (may be NULL) */
                         xmlNode *fieldNodeList, 	/**< [in] linked list of \c xmlNodes that are to be assigned as FIELD children */
                         xmlNode *dataContentNode 	/**< [in] pointer to xmlNode to be inserted under the \<DATA\> element: TABLEDATA, BINARY or FITS */
                         );
  
  /* Attach PARAMs to TABLE node */
  
  
  return(VOTtableNode);
}

/**
 * \brief Serializes a \c LALInferenceVariables structure into a VOTable XML %node
 *
 * This function takes a \c LALInferenceVariables structure and serializes it into a VOTable
 * \c RESOURCE %node identified by the given name. The returned \c xmlNode can then be
 * embedded into an existing %node hierarchy or turned into a full VOTable document.
 *
 * \param vars [in] Pointer to the \c LALInferenceVariables structure to be serialized
 * \param name [in] Unique identifier of this particular \c LALInferenceVariables structure instance
 *
 * \return A pointer to a \c xmlNode that holds the VOTable fragment that represents
 * the \c LALInferenceVariables structure.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * fragment isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \sa LALInferenceVariableItem2VOTParamNode
 *
 * \author John Veitch\n
 * 
 */
xmlNodePtr XLALInferenceVariables2VOTNode (const LALInferenceVariables *const vars,const char *name)
{
  
  /* set up local variables */
  const char *fn = __func__;
  xmlNodePtr xmlChildNodeList = NULL;
  LALInferenceVariableItem *marker=vars->head;
  xmlNodePtr *xmlChildNodePtr=&xmlChildNodeList;

  
  /* Walk through the LALInferenceVariables adding each one */
  while(marker){
    *xmlChildNodePtr = (LALInferenceVariableItem2VOTParamNode(marker));
        if(!*xmlChildNodePtr) {
        /* clean up */
        xmlFreeNodeList(xmlChildNodeList);
        XLALPrintError("Couldn't create PARAM node: %s.%s\n", name,marker->name);
        XLAL_ERROR_NULL(fn, XLAL_EFAILED);
    }
    marker=marker->next;
    xmlChildNodePtr = &((*xmlChildNodePtr)->next);
  }
  return(xmlChildNodeList);
}

/**
 * \brief Serializes a \c LALInferenceVariableItem structure into a VOTable XML %node
 *
 * This function takes a \c LALInferenceVariableItem structure and serializes it into a VOTable
 * \c RESOURCE %node identified by the given name. The returned \c xmlNode can then be
 * embedded into an existing %node hierarchy or turned into a full VOTable document.
 *
 * \param varitem [in] Pointer to the \c LALInferenceVariables structure to be serialized
 * 
 * \return A pointer to a \c xmlNode that holds the VOTable fragment that represents
 * the \c LALInferenceVariableItem structure.
 * In case of an error, a null-pointer is returned.\n
 * \b Important: the caller is responsible to free the allocated memory (when the
 * fragment isn't needed anymore) using \c xmlFreeNode. Alternatively, \c xmlFreeDoc
 * can be used later on when the returned fragment has been embedded in a XML document.
 *
 * \sa XLALCreateVOTParamNode
 * 
 * \author John Veitch\n
 * 
 */


xmlNodePtr LALInferenceVariableItem2VOTParamNode(LALInferenceVariableItem *varitem)
{
  VOTABLE_DATATYPE vo_type;
  CHAR *unitName={0};
  CHAR valString[VARVALSTRINGSIZE_MAX]="";
  
  /* Special case for matrix */
  if(varitem->type==LALINFERENCE_gslMatrix_t)
    return(XLALgsl_matrix2VOTNode((gsl_matrix *)varitem->value, varitem->name, unitName));

  /* Check the type of the item */
  vo_type=LALInferenceVariableType2VOT(varitem->type);
  
  LALInferencePrintVariableItem(valString, varitem);
  
  return(XLALCreateVOTParamNode(varitem->name,unitName,vo_type,NULL,valString));
  
}

/**
 * \brief Convert a \c LALInferenceVariableType into a VOType
 */
VOTABLE_DATATYPE LALInferenceVariableType2VOT(const LALInferenceVariableType litype){
  
  switch(litype){
    case LALINFERENCE_INT4_t: 		return VOT_INT4;
    case LALINFERENCE_INT8_t: 		return VOT_INT8;
    case LALINFERENCE_UINT4_t: 		return VOT_INT8; /* Need a signed INT8 to store an unsigned UINT4 */
    case LALINFERENCE_REAL4_t:		return VOT_REAL4;
    case LALINFERENCE_REAL8_t:		return VOT_REAL8;
    case LALINFERENCE_COMPLEX8_t: 	return VOT_COMPLEX8;
    case LALINFERENCE_COMPLEX16_t:	return VOT_COMPLEX16;
    case LALINFERENCE_string_t:		return VOT_CHAR;
    default: {XLALPrintError("Unsupported LALInferenceVarableType"); return VOT_DATATYPE_LAST;}
  }
}

  