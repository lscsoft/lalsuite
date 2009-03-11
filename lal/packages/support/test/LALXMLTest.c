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

#include <lal/LALXML.h>

#define LALXMLC_ENOM 0
#define LALXMLC_EFUN 1
#define LALXMLC_MSGENOM "Nominal exit"
#define LALXMLC_MSGEFUN "Subroutine returned error"


int main(void)
{
	LIGOTimeGPS time;

	fprintf(stderr, "Running LALXMLTest...\n");

	fprintf(stderr, "LALXMLTest: testing XLALLIGOTimeGPS2XML...\n");
	time.gpsSeconds = 15;
	time.gpsNanoSeconds = 200;

	xmlChar *xml = (xmlChar *) XLALLIGOTimeGPS2VOTableXML(&time, "test");

	if(!xml) {
		fprintf(stderr, "LALXMLTest: [XLALLIGOTimeGPS2XML: %s]\n", LALXMLC_MSGEFUN);
		return LALXMLC_EFUN;
	}

	fprintf(stderr, "%s", xml);
	xmlFree(xml);

	fprintf(stderr, "LALXMLTest: [XLALLIGOTimeGPS2XML: %s]\n", LALXMLC_MSGENOM);

	/* run XLALXML2LIGOTimeGPS test */

	fprintf(stderr, "LALXMLTest: %s\n", LALXMLC_MSGENOM);
	return LALXMLC_ENOM;
}
