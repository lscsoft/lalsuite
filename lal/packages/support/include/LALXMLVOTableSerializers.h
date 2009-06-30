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
 * \brief Header file declaring the public VOTable serializers XML API
 */

/* Double-include protection */
#ifndef _LALXMLVOTABLESERIALIZERS_H
#define _LALXMLVOTABLESERIALIZERS_H

/* C++ protection */
#ifdef __cplusplus
extern "C" {
#endif


#include <libxml/tree.h>

#include <lal/LALDatatypes.h>
#include <lal/PulsarDataTypes.h>


xmlNodePtr XLALLIGOTimeGPS2VOTNode(const LIGOTimeGPS *const ltg, const char *name);
INT4 XLALVOTDoc2LIGOTimeGPSByName(const xmlDocPtr xmlDocument, const char *name, LIGOTimeGPS *ltg);

xmlNodePtr XLALBinaryOrbitParams2VOTNode(const BinaryOrbitParams *const bop, const char *name);
INT4 XLALVOTDoc2BinaryOrbitParamsByName(const xmlDocPtr xmlDocument, const char *name, BinaryOrbitParams *bop);

xmlNodePtr XLALPulsarSpins2VOTNode(const PulsarSpins *const spins, const char *name);
INT4 XLALVOTDoc2PulsarSpinsByName(const xmlDocPtr xmlDocument,
                                  const char *resourceType,
                                  const char *resourceName,
                                  const char *paramName,
                                  PulsarSpins spins);

xmlNodePtr XLALPulsarDopplerParams2VOTNode(const PulsarDopplerParams *const pdp, const char *name);
INT4 XLALVOTDoc2PulsarDopplerParamsByName(const xmlDocPtr xmlDocument, const char *name, PulsarDopplerParams *pdp);

xmlNodePtr XLALgsl_vector2VOTNode(const gsl_vector *vect, const char *name, const CHAR *unitName );
gsl_vector *XLALVOTDoc2gsl_vectorByName(const xmlDocPtr xmlDocument, const char *resourceType, const char *resourceName, const char *paramName, const CHAR *unitName);

xmlNodePtr XLALgsl_matrix2VOTNode(const gsl_matrix *vect, const char *name, const CHAR *unitName );
gsl_matrix *XLALVOTDoc2gsl_matrixByName(const xmlDocPtr xmlDocument, const char *resourceType, const char *resourceName, const char *paramName, const CHAR *unitName);



/* C++ protection */
#ifdef __cplusplus
}
#endif

/* Double-include protection */
#endif
