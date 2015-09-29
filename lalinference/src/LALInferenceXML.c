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

#include <string.h>

#include <lal/XLALError.h>

#include <lal/LALInferenceXML.h>

#define INT4STR_MAXLEN          15
#define REAL8STR_MAXLEN         25
#define NAMESTR_MAXLEN          256


const xmlChar votype_nested_samples[]="lalinference:nestedsampling:samples";
const xmlChar votype_posterior_samples[]="lalinference:posterior:samples";



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
 * \param tablename UNDOCUMENTED
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
 * \author John Veitch
 *
 */
xmlNodePtr XLALInferenceVariablesArray2VOTTable(LALInferenceVariables * const *const varsArray, UINT4 N, const char *tablename)
{
  xmlNodePtr fieldNodeList=NULL;
  xmlNodePtr paramNodeList=NULL;
  xmlNodePtr xmlTABLEDATAnode=NULL;
  xmlNodePtr VOTtableNode=NULL;
  xmlNodePtr tmpNode=NULL;
  xmlNodePtr field_ptr,param_ptr;
  LALInferenceVariableItem *varitem=NULL;
  UINT4 Nfields=0;
  UINT4 bufsize=1024;
  int err;


	/* Sanity check input */
	if(!varsArray) {
		XLALPrintError("Received null varsArray pointer");
		XLAL_ERROR_NULL(XLAL_EFAULT);
	}
	if(N==0) return(NULL);

	field_ptr=fieldNodeList;
	param_ptr=paramNodeList;
    char *field_names[varsArray[0]->dimension];

    /* Build a list of PARAM and FIELD elements */
    for(varitem=varsArray[0]->head;varitem;varitem=varitem->next)
	{
		tmpNode=NULL;
		switch(varitem->vary){
			case LALINFERENCE_PARAM_LINEAR:
			case LALINFERENCE_PARAM_CIRCULAR:
			case LALINFERENCE_PARAM_OUTPUT:
			{
				tmpNode=LALInferenceVariableItem2VOTFieldNode(varitem);
				if(!tmpNode) {
					XLALPrintWarning ("%s: xmlAddNextSibling() failed to add field node for %s.\n", __func__, varitem->name );
					//XLAL_ERROR_NULL(XLAL_EFAILED);
					continue;
				}
				if(field_ptr) field_ptr=xmlAddNextSibling(field_ptr,tmpNode);
				else {field_ptr=tmpNode; fieldNodeList=field_ptr;}
				field_names[Nfields]=varitem->name;
				Nfields++;
				break;
			}
			case LALINFERENCE_PARAM_FIXED:
			{
				tmpNode=LALInferenceVariableItem2VOTParamNode(varitem);
				if(!tmpNode) {
					XLALPrintWarning ("%s: xmlAddNextSibling() failed to add param node for %s.\n", __func__, varitem->name );
					//XLAL_ERROR_NULL(XLAL_EFAILED);
					continue;
				}
				if(param_ptr) param_ptr=xmlAddNextSibling(param_ptr,tmpNode);
				else {param_ptr=tmpNode; paramNodeList=param_ptr;}
				break;
			}
			default:
			{
				XLALPrintWarning("Unknown param vary type");
			}
		}
	}

	if(Nfields>0)
	{
			UINT4 row,col;
			/* create TABLEDATA node */
			if ( ( xmlTABLEDATAnode = xmlNewNode ( NULL, CAST_CONST_XMLCHAR("TABLEDATA") ))== NULL ) {
					XLALPrintError ("%s: xmlNewNode() failed to create 'TABLEDATA' node.\n", __func__ );
					err = XLAL_ENOMEM;
					goto failed;
			}
			/* ---------- loop over data-arrays and generate each table-row */
			for ( row = 0; row < N; row ++ )
			{
					/* create TR node */
					xmlNodePtr xmlThisRowNode = NULL;
					if ( (xmlThisRowNode = xmlNewNode ( NULL, CAST_CONST_XMLCHAR("TR") )) == NULL ) {
							XLALPrintError ("%s: xmlNewNode() failed to create new 'TR' node.\n", __func__ );
							err = XLAL_EFAILED;
							goto failed;
					}
					if ( xmlAddChild(xmlTABLEDATAnode, xmlThisRowNode ) == NULL ) {
							XLALPrintError ("%s: failed to insert 'TR' node into 'TABLEDATA' node.\n", __func__ );
							err = XLAL_EFAILED;
							goto failed;
					}

					/* ----- loop over columns and generate each table element */
					for ( col = 0; col < Nfields; col ++ )
					{
							/* create TD node */
							xmlNodePtr xmlThisEntryNode = NULL;
							if ( (xmlThisEntryNode = xmlNewNode ( NULL, CAST_CONST_XMLCHAR("TD") )) == NULL ) {
									XLALPrintError ("%s: xmlNewNode() failed to create new 'TD' node.\n", __func__ );
									err = XLAL_EFAILED;
									goto failed;
							}
							if ( xmlAddChild(xmlThisRowNode, xmlThisEntryNode ) == NULL ) {
									XLALPrintError ("%s: failed to insert 'TD' node into 'TR' node.\n", __func__ );
									err = XLAL_EFAILED;
									goto failed;
							}

							char *valuestr=XLALCalloc(bufsize,sizeof(char));
							varitem = LALInferenceGetItem(varsArray[row],field_names[col]);
							UINT4 required_size=LALInferencePrintNVariableItem(valuestr,bufsize,varitem);
							if(required_size>bufsize)
							{
								bufsize=required_size;
								valuestr=XLALRealloc(valuestr,required_size*sizeof(char));
								required_size=LALInferencePrintNVariableItem(valuestr,bufsize,varitem);
							}

							xmlNodePtr xmlTextNode= xmlNewText (CAST_CONST_XMLCHAR(valuestr) );
							if ( xmlTextNode  == NULL ) {
									XLALPrintError("%s: xmlNewText() failed to turn text '%s' into node\n", __func__, valuestr );
									err = XLAL_EFAILED;
									XLALFree(valuestr);
									goto failed;
							}
							if ( xmlAddChild(xmlThisEntryNode, xmlTextNode ) == NULL ) {
									XLALPrintError ("%s: failed to insert text-node node into 'TD' node.\n", __func__ );
									err = XLAL_EFAILED;
									XLALFree(valuestr);
									goto failed;
							}
							XLALFree(valuestr);

					} /* for col < numFields */

			} /* for row < numRows */
	}

  /* Create a TABLE from the FIELDs, PARAMs, and TABLEDATA nodes */
  VOTtableNode= XLALCreateVOTTableNode (tablename, fieldNodeList, paramNodeList, xmlTABLEDATAnode );

  return(VOTtableNode);

  failed:
      XLAL_ERROR_NULL ( err );

  return(NULL);

}

xmlNodePtr XLALInferenceStateVariables2VOTResource(LALInferenceRunState *const state, const char *name)
{
	xmlNodePtr algNode=NULL;
	xmlNodePtr priorNode=NULL;
	xmlNodePtr resNode=NULL;
	/* Serialise various params to VOT Table nodes */
	resNode=XLALCreateVOTResourceNode("lalinference:state",name,NULL);
	algNode=XLALInferenceVariablesArray2VOTTable(&(state->algorithmParams),1, "Algorithm Params");
	if(algNode) {
		xmlNewProp(algNode, CAST_CONST_XMLCHAR("utype"), CAST_CONST_XMLCHAR("lalinference:state:algorithmparams"));
		xmlAddChild(resNode,algNode);
	}
	priorNode=XLALInferenceVariablesArray2VOTTable(&(state->priorArgs),1,"Prior Arguments");
    if(priorNode){
		xmlNewProp(priorNode, CAST_CONST_XMLCHAR("utype"), CAST_CONST_XMLCHAR("lalinference:state:priorparams"));
		xmlAddChild(resNode,priorNode);
	}
	return(resNode);

}

/**
 * \brief Serializes a \c LALInferenceVariables structure into a VOTable XML %node
 *
 * This function takes a \c LALInferenceVariables structure and serializes it into a VOTable
 * \c PARAM %node identified by the given name. The returned \c xmlNode can then be
 * embedded into an existing %node hierarchy or turned into a full VOTable document.
 *
 * \param vars [in] Pointer to the \c LALInferenceVariables structure to be serialized
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
 * \author John Veitch
 *
 */
xmlNodePtr XLALInferenceVariables2VOTParamNode (LALInferenceVariables *const vars)
{

  /* set up local variables */
  /*const char *fn = __func__;*/
  xmlNodePtr xmlChildNodeList = NULL;
  xmlNodePtr xmlChild;
  LALInferenceVariableItem *marker=vars->head;

  /* Walk through the LALInferenceVariables adding each one */
  while(marker){
    xmlChild = (LALInferenceVariableItem2VOTParamNode(marker));
    marker=marker->next;
    if(!xmlChild) {
					/* clean up */
					/* if(xmlChildNodeList) xmlFreeNodeList(xmlChildNodeList); */
					XLALPrintWarning("Couldn't create PARAM node for %s\n", marker->name);
					/* XLAL_ERROR_NULL(fn, XLAL_EFAILED); */
					continue;
    }
    if(!xmlChildNodeList) xmlChildNodeList=xmlChild;
    else xmlAddSibling(xmlChildNodeList,xmlChild);

  }
  return(xmlChildNodeList);
}

/**
 * \brief Serializes a \c LALInferenceVariableItem structure into a VOTable XML %node
 *
 * This function takes a \c LALInferenceVariableItem structure and serializes it into a VOTable
 * \c FIELD %node identified by the given name. The returned \c xmlNode can then be
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

xmlNodePtr LALInferenceVariableItem2VOTFieldNode(LALInferenceVariableItem *const varitem)
{
  VOTABLE_DATATYPE vo_type;
  CHAR *unitName={0};

  switch(varitem->type)
  {
		/* Special case for matrix */
		case LALINFERENCE_gslMatrix_t:
				return(XLALgsl_matrix2VOTNode(*(gsl_matrix **)varitem->value, varitem->name, unitName));
		/* Special case for string */
		case LALINFERENCE_string_t:
				return(XLALCreateVOTFieldNode(varitem->name,unitName,VOT_CHAR,"*"));
		case LALINFERENCE_REAL8Vector_t:
		{
				char arraysize[32];
				REAL8Vector *vec=*(REAL8Vector **)varitem->value;
				snprintf(arraysize,32,"%i",vec->length);
				return(XLALCreateVOTFieldNode(varitem->name,unitName,VOT_REAL8,arraysize));
		}
		case LALINFERENCE_INT4Vector_t:
		{
				char arraysize[32];
				INT4Vector *vec=*(INT4Vector **)varitem->value;
				snprintf(arraysize,32,"%i",vec->length);
				return(XLALCreateVOTFieldNode(varitem->name,unitName,VOT_INT4,arraysize));
		}
		default:
		{
				/* Check the type of the item */
				vo_type=LALInferenceVariableType2VOT(varitem->type);
				if(vo_type>=VOT_DATATYPE_LAST){
						XLALPrintError("%s: Unsupported LALInferenceType %i\n",__func__,(int)varitem->type);
						return NULL;
				}
				return(XLALCreateVOTFieldNode(varitem->name,unitName,vo_type,NULL));
		}
  }
}

xmlNodePtr LALInferenceVariableItem2VOTParamNode(LALInferenceVariableItem *const varitem)
{
  VOTABLE_DATATYPE vo_type;
  CHAR *unitName={0};
  UINT4 bufsize=1024;
  UINT4 required_size=0;
  CHAR *valString=XLALCalloc(bufsize,sizeof(CHAR));
  CHAR arraysize[32]="";

  switch(varitem->type)
  {
		  /* Special case for matrix */
		  case LALINFERENCE_gslMatrix_t:
				  return(XLALgsl_matrix2VOTNode(*(gsl_matrix **)varitem->value, varitem->name, unitName));
				  break;

		  /* Special case for string */
		  case LALINFERENCE_string_t:
				  return(XLALCreateVOTParamNode(varitem->name,unitName,VOT_CHAR,"*",*(char **)varitem->value));
				  break;

		  case LALINFERENCE_REAL8Vector_t:
		  {
				  REAL8Vector *vec=*(REAL8Vector **)varitem->value;
				  snprintf(arraysize,32,"%i",vec->length);
				  vo_type=VOT_REAL8;
				  break;
		  }
		  case LALINFERENCE_INT4Vector_t:
		  {
				  INT4Vector *vec=*(INT4Vector **)varitem->value;
				  snprintf(arraysize,32,"%i",vec->length);
				  vo_type=VOT_INT4;
				  break;
		  }
		  default:
		  /* Check the type of the item */
		  vo_type=LALInferenceVariableType2VOT(varitem->type);
  }
  if(vo_type>=VOT_DATATYPE_LAST){
		  XLALPrintError("%s: Unsupported LALInferenceVariableType %i\n",__func__,(int)varitem->type);
		  return NULL;
  }
  required_size=LALInferencePrintNVariableItem(valString, bufsize, varitem);
  if(required_size>bufsize){
		  bufsize=required_size;
		  valString=XLALRealloc(valString,bufsize*sizeof(CHAR));
		  required_size=LALInferencePrintNVariableItem(valString, bufsize, varitem);
  }
  xmlNodePtr node = XLALCreateVOTParamNode(varitem->name,unitName,vo_type,arraysize,valString);
  XLALFree(valString);
  return(node);
}

/**
 * \brief Convert a \c LALInferenceVariableType into a VOType
 */
VOTABLE_DATATYPE LALInferenceVariableType2VOT(LALInferenceVariableType litype){

  switch(litype){
    case LALINFERENCE_INT4_t: 		return VOT_INT4;
    case LALINFERENCE_INT8_t: 		return VOT_INT8;
    case LALINFERENCE_UINT4_t: 		return VOT_INT4; /* Need a signed INT8 to store an unsigned UINT4 */
    case LALINFERENCE_REAL4_t:		return VOT_REAL4;
    case LALINFERENCE_REAL8_t:		return VOT_REAL8;
    case LALINFERENCE_COMPLEX8_t: 	return VOT_COMPLEX8;
    case LALINFERENCE_COMPLEX16_t:	return VOT_COMPLEX16;
    case LALINFERENCE_string_t:		return VOT_CHAR;
    default: {XLALPrintError("%s: Unsupported LALInferenceVarableType %i\n",__func__,(int)litype); return VOT_DATATYPE_LAST;}
  }
}
