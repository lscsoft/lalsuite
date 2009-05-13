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

#include <libxml/tree.h>

#include <lal/LALDatatypes.h>
#include <lal/PulsarDataTypes.h>


xmlNodePtr XLALLIGOTimeGPS2VOTableNode(const LIGOTimeGPS *const ltg, const char *name);
xmlChar * XLALLIGOTimeGPS2VOTableXML(const LIGOTimeGPS *const ltg, const char *name);
INT4 XLALVOTableDoc2LIGOTimeGPSByName(xmlDocPtr xmlDocument, const char *name, LIGOTimeGPS *ltg);
INT4 XLALVOTableXML2LIGOTimeGPSByName(const char *xml, const char *name, LIGOTimeGPS *ltg);

xmlNodePtr XLALBinaryOrbitParams2VOTableNode(const BinaryOrbitParams *const bop, const char *name);
xmlChar * XLALBinaryOrbitParams2VOTableXML(const BinaryOrbitParams *const bop, const char *name);
INT4 XLALVOTableDoc2BinaryOrbitParamsByName(xmlDocPtr xmlDocument, const char *name, BinaryOrbitParams *bop);
INT4 XLALVOTableXML2BinaryOrbitParamsByName(const char *xml, const char *name, BinaryOrbitParams *bop);

xmlNodePtr XLALPulsarDopplerParams2VOTableNode(const PulsarDopplerParams *const pdp, const char *name);
xmlChar * XLALPulsarDopplerParams2VOTableXML(const PulsarDopplerParams *const pdp, const char *name);
INT4 XLALVOTableDoc2PulsarDopplerParamsByName(xmlDocPtr xmlDocument, const char *name, PulsarDopplerParams *pdp);
INT4 XLALVOTableXML2PulsarDopplerParamsByName(const char *xml, const char *name, PulsarDopplerParams *pdp);
