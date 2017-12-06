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
/* Double-include protection */
#ifndef _LALXMLVOTABLESERIALIZERLALINFERENCE_H
#define _LALXMLVOTABLESERIALIZERLALINFERENCE_H

/* C++ protection */
#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALInferenceXML_h Header LALInferenceXML.h
 * \ingroup pkg_LALInferenceXML
 * \brief Header file declaring the public VOTable serializers XML LALInference API
 */
/*@{*/

/*#include <libxml/tree.h>
#include <libxml/xpath.h>
*/
#include <lal/LALXML.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXMLVOTableSerializers.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInference.h>

xmlNodePtr XLALInferenceStateVariables2VOTResource(const LALInferenceRunState *state, const char *name);

xmlNodePtr XLALInferenceVariablesArray2VOTTable(const LALInferenceVariables *varsArray, UINT4 N, const char *name);

xmlNodePtr XLALInferenceVariables2VOTParamNode(const LALInferenceVariables *const vars);

/* INT4 XLALVOTDoc2LALInferenceVariablesByName(const xmlDocPtr xmlDocument, const char *name, LALInferenceVariables *bop); */

/* xmlNodePtr LALInferenceVariables2VOTParamNodes(LALInferenceVariables *vars); */

xmlNodePtr LALInferenceVariableItem2VOTParamNode(LALInferenceVariableItem *varitem);

VOTABLE_DATATYPE LALInferenceVariableType2VOT(const LALInferenceVariableType litype);

xmlNodePtr LALInferenceVariableItem2VOTFieldNode(LALInferenceVariableItem *varitem);

/*@}*/

/* C++ protection */
#ifdef __cplusplus
}
#endif

/* Double-include protection */
#endif
