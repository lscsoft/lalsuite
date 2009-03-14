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

#include <lal/LALRCSID.h>
#include <libxml/tree.h>
#include <lal/LALDatatypes.h>

NRCSID( LALXMLH, "$Id$" );

int XLALXMLFilePrintElements(const char *fname);
xmlDocPtr XLALCreateVOTableXMLFromTree(const xmlNodePtr xmlTree);
xmlChar * XLALGetSingleNodeContentByXPath(const xmlDocPtr xmlDoc, const char *xpath);
xmlNodePtr XLALLIGOTimeGPS2VOTableNode(const LIGOTimeGPS *const ltg, const char *name);
INT4 XLALVOTableXML2LIGOTimeGPSByName(const char *xml, const char *name, LIGOTimeGPS *ltg);
xmlChar * XLALLIGOTimeGPS2VOTableXML(const LIGOTimeGPS *const ltg, const char *name);
