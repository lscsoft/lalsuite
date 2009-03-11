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
#define LALXMLC_EVAL 2
#define LALXMLC_MSGENOM "Nominal exit"
#define LALXMLC_MSGEFUN "Subroutine returned error"
#define LALXMLC_MSGEVAL "Result validation failed"
#define LALXMLC_NAME "LALXMLTest"

int main(void)
{
	/* set up local variables */
	LIGOTimeGPS timeSource;
	LIGOTimeGPS timeDestination;

	/* initialize test data */
	timeSource.gpsSeconds = 15;
	timeSource.gpsNanoSeconds = 200;
	timeDestination.gpsSeconds = 0;
	timeDestination.gpsNanoSeconds = 0;

	fprintf(stderr, "Running LALXMLTest...\n");
	fprintf(stderr, "LALXMLTest: [XLALLIGOTimeGPS2XML: starting...]\n");

	/* invoke and check first function (set) */
	xmlChar *xml = (xmlChar *) XLALLIGOTimeGPS2VOTableXML(&timeSource, LALXMLC_NAME);
	if(!xml) {
		fprintf(stderr, "LALXMLTest: [XLALLIGOTimeGPS2XML: %s]\n", LALXMLC_MSGEFUN);
		return LALXMLC_EFUN;
	}

	/* debug mode only */
	/* fprintf(stderr, xml); */

	fprintf(stderr, "LALXMLTest: [XLALLIGOTimeGPS2XML: %s]\n", LALXMLC_MSGENOM);
	fprintf(stderr, "LALXMLTest: [XLALVOTableXML2LIGOTimeGPSByName: starting...]\n");

	/* invoke and check second function (set) */
	if(XLALVOTableXML2LIGOTimeGPSByName(xml, LALXMLC_NAME, &timeDestination)) {
		fprintf(stderr, "LALXMLTest: [XLALVOTableXML2LIGOTimeGPSByName: %s]\n", LALXMLC_MSGEFUN);
		return LALXMLC_EFUN;
	}

	/* clean up */
	xmlFree(xml);

	/* validate test results */
	if(
		timeSource.gpsSeconds != timeDestination.gpsSeconds ||
		timeSource.gpsNanoSeconds != timeDestination.gpsNanoSeconds)
	{
		fprintf(stderr, "LALXMLTest: [XLALVOTableXML2LIGOTimeGPSByName: %s]\n", LALXMLC_MSGEVAL);
		return LALXMLC_EVAL;
	}

	fprintf(stderr, "LALXMLTest: [XLALVOTableXML2LIGOTimeGPSByName: %s]\n", LALXMLC_MSGENOM);
	fprintf(stderr, "LALXMLTest: %s\n", LALXMLC_MSGENOM);

	return LALXMLC_ENOM;
}
