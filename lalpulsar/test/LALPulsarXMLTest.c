/*
 *  Copyright (C) 2009 Oliver Bock
 *  Copyright (C) 2009 Reinhard Prix
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

#include <config.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALStdio.h>
#include <lal/LogPrintf.h>
#include <lal/LALXML.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXMLVOTableSerializersPulsar.h>
#include <lal/LALStatusMacros.h>


#define LALXMLC_ENOM 0
#define LALXMLC_EFUN 1
#define LALXMLC_EVAL 2
#define LALXMLC_MSGENOM "Nominal exit"
#define LALXMLC_MSGEFUN "Subroutine returned error"
#define LALXMLC_MSGEVAL "Result validation failed"

#define LALXMLC_NAMETEST1 "timeGPS"
#define LALXMLC_NAMETEST2 "cand1"

#define LALXMLC_TYPETEST3 "gsl_vector_wrapper"
#define LALXMLC_NAMETEST3 "myVect"

#define LALXMLC_TYPETEST4 "gsl_matrix_wrapper"
#define LALXMLC_NAMETEST4 "myMatrix"

#define LALXMLC_TYPETEST5 "table_wrapper"
#define LALXMLC_NAMETEST5 "testTable"



#define REAL8TOL 1e-16
#define REAL4TOL 1e-6
#define PATH_MAXLEN 256

/* private test prototypes */
int testPulsarDopplerParams(void);

/* private utility prototypes */
int validateDocument(const xmlDocPtr xmlDocument);
int findFileInLALDataPath(const char *filename, char **validatedPath);

int main(void)
{
    /* set up local variables */
    int result = LALXMLC_ENOM;

    /* set debug level*/
    lalDebugLevel = LALMSGLVL3;

    printf( "\n" );
    printf( "======================================================================\n");
    printf( "1: Test PulsarDopplerParams (de)serialization...\n");
    printf( "======================================================================\n");
    result = testPulsarDopplerParams();
    if(result != LALXMLC_ENOM) {
      LogPrintf (LOG_CRITICAL, "testPulsarDopplerParams() failed. ret = %d\n", result );
      return result;
    }

    LALCheckMemoryLeaks();

    return LALXMLC_ENOM;

} /* main() */

int testPulsarDopplerParams(void)
{
    static BinaryOrbitParams bopSource;
    static PulsarDopplerParams pdpSource;
    static BinaryOrbitParams bopDestination;
    static PulsarDopplerParams pdpDestination;
    xmlNodePtr xmlFragment = NULL;
    xmlDocPtr xmlDocument = NULL;
    char *xmlString = NULL;
    int result;

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

    printf( "--> Initial PulsarDopplerParams struct: ");
    printf( "%s = { \n"
            "\trefTime: {%d, %d}\n"
            "\tAlpha: %g\n"
            "\tDelta: %g\n"
            "\tfkdot: {%g, %g, %g, %g}\n"
            "\torbit.tp: {%d, %d}\n"
            "\torbit.argp: %g\n"
            "\torbit.asini: %g\n"
            "\torbit.ecc: %g\n"
            "\torbit.period: %g}\n",
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

    /* serialize structure into VOTable fragment */
    printf( "--> Serializing into XML string ... ");
    if ( (xmlFragment = XLALPulsarDopplerParams2VOTNode(&pdpSource, LALXMLC_NAMETEST2)) == NULL ) {
      XLALPrintError( "%s: XLALPulsarDopplerParams2VOTNode() failed. err = %d\n", __func__, xlalErrno );
      return LALXMLC_EFUN;
    }
    /* convert VOTable tree into XML string */
    if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
      XLALPrintError( "%s: XLALCreateVOTStringFromTree() failed. err = %d\n", __func__, xlalErrno );
      return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* display serialized structure */
    printf( "----------------------------------------------------------------------\n");
    printf( "%s", xmlString );
    printf( "----------------------------------------------------------------------\n");

    /* convert XML string back into VOTable document */
    printf ("--> Parsing XML string into xmlDoc ... ");
    if( ( xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
      XLALPrintError( "%s: XLALXMLString2Doc() failed. err = %d\n", __func__, xlalErrno );
      return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* validate XML document */
    printf ("--> Validating xmlDoc against VOTable schema ... ");
    if ( (result = validateDocument(xmlDocument)) != XLAL_SUCCESS ) {
      XLALPrintError("%s: validateDocument failed. ret = %d\n", __func__, result );
      if ( result == XLAL_FAILURE )
        return LALXMLC_EVAL;
      else
        return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* parse VOTable document back into C-structure */
    printf ("--> Read out PulsarDopplerParams struct ... ");
    if ( XLALVOTDoc2PulsarDopplerParamsByName ( xmlDocument, LALXMLC_NAMETEST2, &pdpDestination)) {
      XLALPrintError ( "%s: XLALVOTDoc2LIGOTimeGPSByName() failed. errno = %d\n", __func__, xlalErrno );
      return LALXMLC_EFUN;
    }

    printf ( "ok: %s = {\n"
            "\trefTime: {%d, %d}\n"
            "\tAlpha: %g\n"
            "\tDelta: %g\n"
            "\tfkdot: {%g, %g, %g, %g}\n"
            "\torbit.tp: {%d, %d}\n"
            "\torbit.argp: %g\n"
            "\torbit.asini: %g\n"
            "\torbit.ecc: %g\n"
            "\torbit.period: %g}\n",
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
      XLALPrintError ( "%s: PulsarDopplerParams structure differs before and after XML serialization!\n", __func__);
      return LALXMLC_EVAL;
    }

    printf ( "--> Final struct agrees with input struct. Test PASSED.\n");

    /* clean up */
    xmlFreeDoc(xmlDocument);
    xmlFreeNode ( xmlFragment );
    XLALFree ( xmlString );

    return LALXMLC_ENOM;

} /* testPulsarDopplerParams() */


/* -------------------- Helper functions -------------------- */

int validateDocument(const xmlDocPtr xmlDocument)
{
    /* set up local variables */
    char schemaUrl[] = "file://" DATADIR "VOTable-1.1.xsd";
    int result;

    /* validate document */
    printf ("schema-file: %s ... ", schemaUrl );
    result = XLALValidateDocumentByExternalSchema(xmlDocument, BAD_CAST(schemaUrl));

    return result;

} /* validateDocument() */
