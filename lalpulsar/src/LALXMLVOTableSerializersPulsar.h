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
/* Double-include protection */
#ifndef _LALXMLVOTABLESERIALIZERSPULSAR_H
#define _LALXMLVOTABLESERIALIZERSPULSAR_H

/* C++ protection */
#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALXMLVOTableSerializersPulsar_h Header LALXMLVOTableSerializersPulsar.h
 * \ingroup pkg_pulsarXML
 * \brief Header file declaring the public VOTable serializers XML Pulsar API
 */
/*@{*/

#include <libxml/tree.h>

#include <lal/LALDatatypes.h>
#include <lal/PulsarDataTypes.h>


xmlNodePtr XLALPulsarSpins2VOTNode(const PulsarSpins *const spins, const char *name);
INT4 XLALVOTDoc2PulsarSpinsByName(const xmlDocPtr xmlDocument,
                                  const char *resourcePath,
                                  const char *paramName,
                                  PulsarSpins spins);

xmlNodePtr XLALPulsarDopplerParams2VOTNode(const PulsarDopplerParams *const pdp, const char *name);
INT4 XLALVOTDoc2PulsarDopplerParamsByName(const xmlDocPtr xmlDocument, const char *name, PulsarDopplerParams *pdp);


/*@}*/

/* C++ protection */
#ifdef __cplusplus
}
#endif

/* Double-include protection */
#endif
