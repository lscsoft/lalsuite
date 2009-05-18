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

#ifndef LAL_PREFIX
    #define LAL_PREFIX "/usr/local"
#endif

#define PATH_MAXLEN 256


INT4 testLIGOTimeGPS();
INT4 testPulsarDopplerParams();
INT4 validateDocumentString(const xmlChar *xml);
INT4 findFileInLALDataPath(const char *filename, char **validatedPath);


int main(void)
{
    /* set up local variables */
    int result = LALXMLC_ENOM;

    /* set debug level*/
    lalDebugLevel = LALMSGLVL3;

    fprintf(stderr, "**********************************************************************\n");
    fprintf(stderr, "Running LALXMLTest...\n\n");

    result = testLIGOTimeGPS();
    if(result != LALXMLC_ENOM) {
        return result;
    }

    fprintf(stderr, "======================================================================\n\n");

    result = testPulsarDopplerParams();
    if(result != LALXMLC_ENOM) {
        return result;
    }

    fprintf(stderr, "**********************************************************************\n");
    return LALXMLC_ENOM;
}


INT4 testLIGOTimeGPS()
{
    /* set up local variables */
    static LIGOTimeGPS timeSource;
    static LIGOTimeGPS timeDestination;
    xmlChar *xml;
    INT4 result;

    /* initialize test data */
    timeSource.gpsSeconds = 15;
    timeSource.gpsNanoSeconds = 200;

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

    /* validate XML document */
    result = validateDocumentString(xml);
    if(result == XLAL_SUCCESS) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentByInternalSchema(): %s]\n", LALXMLC_MSGENOM);
    }
    else if(result == XLAL_FAILURE) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentByInternalSchema(): %s]\n", LALXMLC_MSGEVAL);
        return LALXMLC_EVAL;
    }
    else {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentByInternalSchema(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }

    /* invoke and check second function (get) */
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


INT4 testPulsarDopplerParams()
{
    /* set up local variables */
    static BinaryOrbitParams bopSource;
    static PulsarDopplerParams pdpSource;
    static BinaryOrbitParams bopDestination;
    static PulsarDopplerParams pdpDestination;
    xmlChar *xml;
    INT4 result;

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

    pdpDestination.orbit = &bopDestination;

    fprintf(stderr, "2: Testing PulsarDopplerParams (de)serialization...\n\n");

    fprintf(stderr, "Initial PulsarDopplerParams struct:\n");
    fprintf(stderr, "%s = {\n"
            "\trefTime: {%d, %d}\n"
            "\tAlpha: %g\n"
            "\tDelta: %g\n"
            "\tfkdot: {%g, %g, %g, %g}\n"
            "\torbit.tp: {%d, %d}\n"
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
            pdpSource.orbit->argp,
            pdpSource.orbit->asini,
            pdpSource.orbit->ecc,
            pdpSource.orbit->period);

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

    /* validate XML document */
    result = validateDocumentString(xml);
    if(result == XLAL_SUCCESS) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentByInternalSchema(): %s]\n", LALXMLC_MSGENOM);
    }
    else if(result == XLAL_FAILURE) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentByInternalSchema(): %s]\n", LALXMLC_MSGEVAL);
        return LALXMLC_EVAL;
    }
    else {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentByInternalSchema(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }

    /* invoke and check second function (get) */
    if(XLALVOTableXML2PulsarDopplerParamsByName((char*)xml, LALXMLC_NAMETEST2, &pdpDestination)) {
        fprintf(stderr, "LALXMLTest: [XLALVOTableXML2PulsarDopplerParamsByName(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }
    fprintf(stderr, "\nLALXMLTest: [XLALVOTableXML2PulsarDopplerParamsByName(): %s]\n\n", LALXMLC_MSGENOM);

    /* clean up */
    xmlFree(xml);

    fprintf(stderr, "PulsarDopplerParams struct parsed back from VOTable:\n");
    fprintf(stderr, "%s = {\n"
            "\trefTime: {%d, %d}\n"
            "\tAlpha: %g\n"
            "\tDelta: %g\n"
            "\tfkdot: {%g, %g, %g, %g}\n"
            "\torbit.tp: {%d, %d}\n"
            "\torbit.argp: %g\n"
            "\torbit.asini: %g\n"
            "\torbit.ecc: %g\n"
            "\torbit.period: %g\n}\n\n",
            LALXMLC_NAMETEST2,
            pdpDestination.refTime.gpsSeconds, pdpDestination.refTime.gpsNanoSeconds,
            pdpDestination.Alpha,
            pdpDestination.Delta,
            pdpDestination.fkdot[0], pdpDestination.fkdot[1], pdpDestination.fkdot[2], pdpDestination.fkdot[3],
            pdpDestination.orbit->tp.gpsSeconds, pdpDestination.orbit->tp.gpsNanoSeconds,
            pdpDestination.orbit->argp,
            pdpDestination.orbit->asini,
            pdpDestination.orbit->ecc,
            pdpDestination.orbit->period);

    /* validate test results */
    if(
            pdpSource.refTime.gpsSeconds != pdpDestination.refTime.gpsSeconds ||
            pdpSource.refTime.gpsNanoSeconds != pdpDestination.refTime.gpsNanoSeconds ||
            pdpSource.Alpha != pdpDestination.Alpha ||
            pdpSource.Delta != pdpDestination.Delta ||
            pdpSource.fkdot[0] != pdpDestination.fkdot[0] ||
            pdpSource.fkdot[1] != pdpDestination.fkdot[1] ||
            pdpSource.fkdot[2] != pdpDestination.fkdot[2] ||
            pdpSource.fkdot[3] != pdpDestination.fkdot[3] ||
            pdpSource.orbit->tp.gpsSeconds != pdpDestination.orbit->tp.gpsSeconds ||
            pdpSource.orbit->tp.gpsNanoSeconds != pdpDestination.orbit->tp.gpsNanoSeconds ||
            pdpSource.orbit->argp != pdpDestination.orbit->argp ||
            pdpSource.orbit->asini != pdpDestination.orbit->asini ||
            pdpSource.orbit->ecc != pdpDestination.orbit->ecc ||
            pdpSource.orbit->period != pdpDestination.orbit->period)
    {
        fprintf(stderr, "LALXMLTest: [XLALVOTableXML2PulsarDopplerParamsByName(): %s]\n\n", LALXMLC_MSGEVAL);
        return LALXMLC_EVAL;
    }

    return LALXMLC_ENOM;
}


INT4 validateDocumentString(const xmlChar *xml)
{
    /* set up local variables */
    xmlDocPtr xmlDocument = NULL;
    char *schemaPath = NULL;
    char schemaUrl[PATH_MAXLEN+10] = "file://";
    INT4 result;

    /* parse XML document */
    xmlDocument = xmlReadMemory(xml, strlen(xml), NULL, "UTF-8", 0);
    if(xmlDocument == NULL) {
        /* clean up */
        xmlCleanupParser();
        fprintf(stderr, "VOTable document parsing failed\n");
        return XLAL_EFAILED;
    }

    /* find schema definition file */
    result = findFileInLALDataPath("VOTable-1.1.xsd", &schemaPath);

    /* validate document */
    if(result == XLAL_SUCCESS) {
        strncat(schemaUrl, schemaPath, PATH_MAXLEN);
        result = XLALValidateDocumentByExternalSchema(xmlDocument, schemaUrl);
        LALFree(schemaPath);
    }
    else {
        fprintf(stderr, "Warning: schema definition file not found! "
                        "Falling back to internal schema definition (online resource)!\n");
        result = XLALValidateDocumentByInternalSchema(xmlDocument);
    }

    /* clean up */
    xmlFreeDoc(xmlDocument);
    xmlCleanupParser();

    return result;
}

INT4 findFileInLALDataPath(const char *filename, char **validatedPath)
{
    /* set up local variables */
    char *absolutePath;
    char workingDir[256] = {0};
    const char *dataPathEnv;
    char *dataPath;
    const char *currentDataPath;
    char *nextDataPath;
    FILE *fileCheck;
    int n;

    /* basic sanity checks */
    if(!filename) {
        fprintf(stderr, "No filename specified!\n");
        return XLAL_EINVAL;
    }
    if(*filename == '/') {
        fprintf(stderr, "Absolute path given!\n");
        return XLAL_EINVAL;
    }
    if(!validatedPath) {
        fprintf(stderr, "No destination buffer specified!\n");
        return XLAL_EINVAL;
    }

    /* allocate buffer for final path */
    if((absolutePath = LALCalloc(PATH_MAXLEN, 1)) == NULL) {
        fprintf(stderr, "Can't allocate memory (%i)!\n", PATH_MAXLEN);
        return XLAL_EFAILED;
    }

    /* get current working directory */
    if(!getcwd(workingDir, PATH_MAXLEN)) {
        fprintf(stderr, "Can't determine current working directory!\n");
        LALFree(absolutePath);
        return XLAL_EFAILED;
    }

    /* get data path (set when using "make check")*/
    dataPathEnv = getenv("LAL_DATA_PATH");

    /* LAL_DATA_PATH unavailable */
    if(!dataPathEnv || !strlen(dataPathEnv)) {
        fprintf(stderr, "Warning: LAL_DATA_PATH not set! Trying working directory...\n");
        fileCheck = LALFopen(filename, "r");
        if(!fileCheck) {
            fprintf(stderr, "Specified file (%s) not found!\n", filename);
            LALFree(absolutePath);
            return XLAL_FAILURE;
        }
        else {
            LALFclose(fileCheck);
        }

        /* build absolute path */
        n = LALSnprintf(absolutePath, PATH_MAXLEN, "%s/./%s", workingDir, filename);
        if(n >= PATH_MAXLEN) {
            /* data file name too long */
            fprintf(stderr, "Absolute path exceeds limit of %i characters!\n", PATH_MAXLEN);
            LALFree(absolutePath);
            return XLAL_EFAILED;
        }

        /* success: return path */
        *validatedPath = absolutePath;
        return XLAL_SUCCESS;
    }

    /* LAL_DATA_PATH available: scan through all directories in colon-delimited list */
    if((dataPath = LALCalloc(strlen(dataPathEnv)+1, 1)) == NULL)
    {
        fprintf(stderr, "Can't allocate memory (%i)!\n", strlen(dataPathEnv)+1);
        return XLAL_EFAILED;
    }

    /* create working copy */
    strcpy(dataPath, dataPathEnv);
    currentDataPath = dataPath;

    do {
        /* look for additional directories */
        nextDataPath = strchr(currentDataPath, ':');
        if(nextDataPath) {
            /* there are more things in the list */
            /* NUL-terminate current directory */
            *nextDataPath++ = 0;
        }
        if(!strlen(currentDataPath)) {
            /* this directory is empty */
            /* default data directory */
            currentDataPath = LAL_PREFIX"/share/lal";
        }

        /* build absolute path (required by "file" URI scheme) */
        n = LALSnprintf(absolutePath,
                        PATH_MAXLEN,
                        "%s/%s/%s",
                        workingDir,
                        currentDataPath ? currentDataPath : ".",
                        filename);

        if(n >= PATH_MAXLEN) {
            /* data file name too long */
            fprintf(stderr, "Absolute path exceeds limit of %i characters!\n", PATH_MAXLEN);
            LALFree(absolutePath);
            LALFree(dataPath);
            return XLAL_EFAILED;
        }

        /* check if file is accessible */
        fileCheck = LALFopen(absolutePath, "r");
        if(fileCheck) {
            /* success: return path */
            *validatedPath = absolutePath;
            LALFclose(fileCheck);
            LALFree(dataPath);
            return XLAL_SUCCESS;
        }

        currentDataPath = nextDataPath;
    }
    while(currentDataPath);

    /* clean up */
    LALFree(dataPath);

    /* return error */
    fprintf(stderr, "Specified file (%s) not found!\n", filename);
    return XLAL_FAILURE;
}
