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

#include <lal/LALXMLVOTableSerializers.h>
#include <LALStatusMacros.h>

#define LALXMLC_ENOM 0
#define LALXMLC_EFUN 1
#define LALXMLC_EVAL 2
#define LALXMLC_MSGENOM "Nominal exit"
#define LALXMLC_MSGEFUN "Subroutine returned error"
#define LALXMLC_MSGEVAL "Result validation failed"

#define LALXMLC_NAMETEST1 "timeGPS"
#define LALXMLC_NAMETEST2 "cand1"


int testLIGOTimeGPS();
int testPulsarDopplerParams();

int main(void)
{
    /* set up local variables */
    int result = LALXMLC_ENOM;

    /* set debug level*/
    lalDebugLevel = LALMSGLVL3 | LALMEMTRACE;

    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "Running LALXMLTest...\n\n");

    result = testLIGOTimeGPS();
    if(result != LALXMLC_ENOM) {
        return result;
    }
    result = testPulsarDopplerParams();
    if(result != LALXMLC_ENOM) {
        return result;
    }

    fprintf(stderr, "======================================================================\n");
    return LALXMLC_ENOM;
}


int testLIGOTimeGPS()
{
    /* set up local variables */
    LIGOTimeGPS timeSource;
    LIGOTimeGPS timeDestination;
    xmlChar *xml;

    /* initialize test data */
    timeSource.gpsSeconds = 15;
    timeSource.gpsNanoSeconds = 200;
    timeDestination.gpsSeconds = 0;
    timeDestination.gpsNanoSeconds = 0;

    fprintf(stderr, "1: Testing LIGOTimeGPS (de)serialization...\n\n");

    fprintf(stderr, "Initial LIGOTimeGPS struct:\n");
    fprintf(stderr, "%s = { %d, %d }\n\n", LALXMLC_NAMETEST1, timeSource.gpsSeconds, timeSource.gpsNanoSeconds );

    /* invoke and check first function (set) */
    xml = (xmlChar *) XLALLIGOTimeGPS2VOTableXML(&timeSource, LALXMLC_NAMETEST1);
    if(!xml) {
        fprintf(stderr, "LALXMLTest: [XLALLIGOTimeGPS2VOTableXML(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }
    fprintf(stderr, "LALXMLTest: [XLALLIGOTimeGPS2VOTableXML(): %s]\n\n", LALXMLC_MSGENOM);

    /* display serialized structure */
    fprintf(stderr, "Serialized VOTable XML:\n");
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, (char*)xml);
    fprintf(stderr, "----------------------------------------------------------------------\n");

    /* invoke and check second function (set) */
    if(XLALVOTableXML2LIGOTimeGPSByName((char*)xml, LALXMLC_NAMETEST1, &timeDestination)) {
        fprintf(stderr, "LALXMLTest: [XLALVOTableXML2LIGOTimeGPSByName(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }
    fprintf(stderr, "\nLALXMLTest: [XLALVOTableXML2LIGOTimeGPSByName(): %s]\n\n", LALXMLC_MSGENOM);

    /* clean up */
    xmlFree(xml);

    fprintf(stderr, "LIGOTimeGPS struct parsed back from VOTable:\n");
    fprintf(stderr, "%s = { %d, %d }\n\n", LALXMLC_NAMETEST1, timeDestination.gpsSeconds, timeDestination.gpsNanoSeconds );

    /* validate test results */
    if(
            timeSource.gpsSeconds != timeDestination.gpsSeconds ||
            timeSource.gpsNanoSeconds != timeDestination.gpsNanoSeconds)
    {
        fprintf(stderr, "LALXMLTest: [XLALVOTableXML2LIGOTimeGPSByName(): %s]\n\n", LALXMLC_MSGEVAL);
        return LALXMLC_EVAL;
    }

    return LALXMLC_ENOM;
}


int testPulsarDopplerParams()
{
    /* set up local variables */
    BinaryOrbitParams bopSource;
    PulsarDopplerParams pdpSource;
    PulsarDopplerParams pdpDestination;
    xmlChar *xml;

    /* initialize test data */
    bopSource.tp.gpsSeconds = 913399939;
    bopSource.tp.gpsNanoSeconds = 15;
    bopSource.argp = 0.5;
    bopSource.asini = 500;
    bopSource.ecc = 0.0167;
    bopSource.period = 31536000;
    pdpSource.refTime.gpsSeconds = 913399939;
    pdpSource.refTime.gpsNanoSeconds = 15;
    pdpSource.Alpha = 3.1452;
    pdpSource.Delta = -0.15;
    pdpSource.fkdot[0] = 100.5;
    pdpSource.fkdot[1] = -1.7e-8;
    pdpSource.fkdot[2] = 0.0;
    pdpSource.fkdot[3] = 0.0;
    pdpSource.orbit = &bopSource;

    fprintf(stderr, "1: Testing PulsarDopplerParams (de)serialization...\n\n");

    fprintf(stderr, "Initial PulsarDopplerParams struct:\n");
    fprintf(stderr, "%s = {\n"
            "\trefTime: {%d, %d}\n"
            "\tAlpha: %g\n"
            "\tDelta: %g\n"
            "\tfkdot: {%g, %g, %g, %g}}\n"
            "\torbit->tp: {%d, %d}\n"
            "\torbit.argp: %g\n"
            "\torbit.asini: %g\n"
            "\torbit.ecc: %g\n"
            "\torbit.period: %g\n}\n\n",
            LALXMLC_NAMETEST2,
            pdpSource.refTime.gpsSeconds, pdpSource.refTime.gpsNanoSeconds,
            pdpSource.Alpha,
            pdpSource.Delta,
            pdpSource.fkdot[0], pdpSource.fkdot[1], pdpSource.fkdot[2], pdpSource.fkdot[3],
            pdpSource.orbit->tp.gpsSeconds, pdpSource.orbit->tp.gpsNanoSeconds,
            bopSource.argp,
            bopSource.asini,
            bopSource.ecc,
            bopSource.period);

    /* invoke and check first function (set) */
    xml = (xmlChar *) XLALPulsarDopplerParams2VOTableXML(&pdpSource, LALXMLC_NAMETEST2);
    if(!xml) {
        fprintf(stderr, "LALXMLTest: [XLALPulsarDopplerParams2VOTableXML(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }
    fprintf(stderr, "LALXMLTest: [XLALPulsarDopplerParams2VOTableXML(): %s]\n\n", LALXMLC_MSGENOM);

    /* display serialized structure */
    fprintf(stderr, "Serialized VOTable XML:\n");
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, (char*)xml);
    fprintf(stderr, "----------------------------------------------------------------------\n");

    return LALXMLC_ENOM;
}
