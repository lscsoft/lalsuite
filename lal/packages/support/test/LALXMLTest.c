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

#include <config.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LogPrintf.h>
#include <lal/LALXML.h>
#include <lal/LALXMLVOTableCommon.h>
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

#define LALXMLC_TYPETEST3 "gsl_vector_wrapper"
#define LALXMLC_NAMETEST3 "myVect"

#define LALXMLC_TYPETEST4 "gsl_matrix_wrapper"
#define LALXMLC_NAMETEST4 "myMatrix"


#define REAL8TOL 1e-15

#ifndef LAL_PREFIX
    #define LAL_PREFIX "/usr/local"
#endif

#define PATH_MAXLEN 256


/* private test prototypes */
int testLIGOTimeGPS(void);
int testPulsarDopplerParams(void);
int test_gsl_vector(void);
int test_gsl_matrix(void);

/* private utility prototypes */
int validateDocument(const xmlDocPtr xmlDocument);
int findFileInLALDataPath(const char *filename, char **validatedPath);

/* \todo factor out -> generic XML */
int xmlDocument2String(const xmlDocPtr xmlDocument, xmlChar **xmlString);
int xmlString2Document(const xmlChar *xmlString, xmlDocPtr *xmlDocument);


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

    fprintf(stderr, "======================================================================\n\n");

    if ( (result = test_gsl_vector()) != LALXMLC_ENOM ) {
      return result;
    }

    fprintf(stderr, "======================================================================\n\n");

    if ( (result = test_gsl_matrix()) != LALXMLC_ENOM ) {
      return result;
    }

    fprintf(stderr, "**********************************************************************\n");

    return LALXMLC_ENOM;

} /* main() */


int testLIGOTimeGPS(void)
{
    /* set up local variables */
    static LIGOTimeGPS timeSource;
    static LIGOTimeGPS timeDestination;
    xmlNodePtr xmlFragment = NULL;
    xmlDocPtr xmlDocument = NULL;
    xmlChar *xmlString = NULL;
    int result;

    /* initialize test data */
    timeSource.gpsSeconds = 15;
    timeSource.gpsNanoSeconds = 200;

    fprintf(stderr, "1: Testing LIGOTimeGPS (de)serialization...\n\n");

    fprintf(stderr, "Initial LIGOTimeGPS struct:\n");
    fprintf(stderr, "%s = { %d, %d }\n\n", LALXMLC_NAMETEST1, timeSource.gpsSeconds, timeSource.gpsNanoSeconds );

    /* serialize structure into VOTable fragment */
    xmlFragment = XLALLIGOTimeGPS2VOTableNode(&timeSource, LALXMLC_NAMETEST1);
    if(!xmlFragment) {
        fprintf(stderr, "LALXMLTest: [XLALLIGOTimeGPS2VOTableNode(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }
    fprintf(stderr, "LALXMLTest: [XLALLIGOTimeGPS2VOTableNode(): %s]\n", LALXMLC_MSGENOM);

    /* convert VOTable fragment into VOTable document */
    xmlDocument = (xmlDocPtr)XLALCreateVOTableDocumentFromTree((const xmlNodePtr)xmlFragment);
    if(!xmlFragment) {
        xmlFree(xmlFragment);
        fprintf(stderr, "LALXMLTest: [XLALCreateVOTableDocumentFromTree(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }
    fprintf(stderr, "LALXMLTest: [XLALCreateVOTableDocumentFromTree(): %s]\n", LALXMLC_MSGENOM);

    /* convert VOTable document into XML string */
    if(!xmlDocument2String(xmlDocument, &xmlString)) {
        xmlFreeDoc(xmlDocument);
        return LALXMLC_EFUN;
    }

    /* display serialized structure */
    fprintf(stderr, "Serialized VOTable XML:\n");
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, (char*)xmlString);
    fprintf(stderr, "----------------------------------------------------------------------\n");

    /* validate XML document */
    result = validateDocument(xmlDocument);
    if(result == XLAL_SUCCESS) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGENOM);
    }
    else if(result == XLAL_FAILURE) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEVAL);
        xmlFree(xmlString);
        xmlFreeDoc(xmlDocument);
        return LALXMLC_EVAL;
    }
    else {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEFUN);
        xmlFree(xmlString);
        xmlFreeDoc(xmlDocument);
        return LALXMLC_EFUN;
    }

    /* convert XML string to VOTable document (not necessary here, serves as an example!) */
    xmlFreeDoc(xmlDocument);
    if(xmlString2Document(xmlString, &xmlDocument)) {
        xmlFree(xmlString);
        return LALXMLC_EFUN;
    }
    else {
        xmlFree(xmlString);
    }

    /* validate XML document (not necessary here, serves as an example!) */
    result = validateDocument(xmlDocument);
    if(result == XLAL_SUCCESS) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGENOM);
    }
    else if(result == XLAL_FAILURE) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEVAL);
        xmlFreeDoc(xmlDocument);
        return LALXMLC_EVAL;
    }
    else {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEFUN);
        xmlFreeDoc(xmlDocument);
        return LALXMLC_EFUN;
    }

    /* deserialize VOTable document into structure */
    if(XLALVOTableDoc2LIGOTimeGPSByName(xmlDocument, LALXMLC_NAMETEST1, &timeDestination)) {
        fprintf(stderr, "LALXMLTest: [XLALVOTableDoc2LIGOTimeGPSByName(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }
    fprintf(stderr, "LALXMLTest: [XLALVOTableDoc2LIGOTimeGPSByName(): %s]\n", LALXMLC_MSGENOM);

    /* clean up */
    xmlFreeDoc(xmlDocument);

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

} /* testLIGOTimeGPS() */


int testPulsarDopplerParams(void)
{
    /* set up local variables */
    static BinaryOrbitParams bopSource;
    static PulsarDopplerParams pdpSource;
    static BinaryOrbitParams bopDestination;
    static PulsarDopplerParams pdpDestination;
    xmlNodePtr xmlFragment = NULL;
    xmlDocPtr xmlDocument = NULL;
    xmlChar *xmlString = NULL;
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

    /* serialize structure into VOTable fragment */
    xmlFragment = XLALPulsarDopplerParams2VOTableNode(&pdpSource, LALXMLC_NAMETEST2);
    if(!xmlFragment) {
        fprintf(stderr, "LALXMLTest: [XLALPulsarDopplerParams2VOTableNode(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }
    fprintf(stderr, "LALXMLTest: [XLALPulsarDopplerParams2VOTableNode(): %s]\n", LALXMLC_MSGENOM);

    /* convert VOTable fragment into VOTable document */
    xmlDocument = (xmlDocPtr)XLALCreateVOTableDocumentFromTree((const xmlNodePtr)xmlFragment);
    if(!xmlFragment) {
        xmlFree(xmlFragment);
        fprintf(stderr, "LALXMLTest: [XLALCreateVOTableDocumentFromTree(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }
    fprintf(stderr, "LALXMLTest: [XLALCreateVOTableDocumentFromTree(): %s]\n", LALXMLC_MSGENOM);

    /* convert VOTable document into XML string */
    if(!xmlDocument2String(xmlDocument, &xmlString)) {
        xmlFreeDoc(xmlDocument);
        return LALXMLC_EFUN;
    }

    /* display serialized structure */
    fprintf(stderr, "Serialized VOTable XML:\n");
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, (char*)xmlString);
    fprintf(stderr, "----------------------------------------------------------------------\n");

    /* validate XML document */
    result = validateDocument(xmlDocument);
    if(result == XLAL_SUCCESS) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGENOM);
    }
    else if(result == XLAL_FAILURE) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEVAL);
        xmlFree(xmlString);
        xmlFreeDoc(xmlDocument);
        return LALXMLC_EVAL;
    }
    else {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEFUN);
        xmlFree(xmlString);
        xmlFreeDoc(xmlDocument);
        return LALXMLC_EFUN;
    }

    /* convert XML string to VOTable document (not necessary here, serves as an example!) */
    xmlFreeDoc(xmlDocument);
    if(xmlString2Document(xmlString, &xmlDocument)) {
        xmlFree(xmlString);
        return LALXMLC_EFUN;
    }
    else {
        xmlFree(xmlString);
    }

    /* validate XML document (not necessary here, serves as an example!) */
    result = validateDocument(xmlDocument);
    if(result == XLAL_SUCCESS) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGENOM);
    }
    else if(result == XLAL_FAILURE) {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEVAL);
        xmlFreeDoc(xmlDocument);
        return LALXMLC_EVAL;
    }
    else {
        fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEFUN);
        xmlFreeDoc(xmlDocument);
        return LALXMLC_EFUN;
    }

    /* deserialize VOTable document into structure */
    if(XLALVOTableDoc2PulsarDopplerParamsByName(xmlDocument, LALXMLC_NAMETEST2, &pdpDestination)) {
        fprintf(stderr, "LALXMLTest: [XLALVOTableDoc2PulsarDopplerParamsByName(): %s]\n", LALXMLC_MSGEFUN);
        return LALXMLC_EFUN;
    }
    fprintf(stderr, "LALXMLTest: [XLALVOTableDoc2PulsarDopplerParamsByName(): %s]\n", LALXMLC_MSGENOM);

    /* clean up */
    xmlFreeDoc(xmlDocument);

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

} /* testPulsarDopplerParams() */

/* test (de-)serialization of a gsl_vector type */
int
test_gsl_vector(void)
{
  static const char *fn = "test_gsl_vector()";

  int i, dim;
  gsl_vector *in_vect, *out_vect;
  xmlNodePtr xmlChildNodeList = NULL;

  xmlNodePtr xmlFragment = NULL;
  xmlDocPtr xmlDocument = NULL;
  xmlChar *xmlString = NULL;
  int result;
  REAL8 err, maxerr = 0;

  /* ---------- initialize input test data ---------- */

  /* pick a random number of dimensions between [1, 11] */
  srand(time(NULL));	/* pick random seed */
  dim = 1 + (1.0 * rand() / RAND_MAX) * 10;

  if ( (in_vect = gsl_vector_alloc ( dim )) == NULL ) {
    XLALPrintError ("%s: failed to gsl_vector_alloc(%d).\n\n", fn, dim );
    return XLAL_ENOMEM;
  }
  for (i=0; i < dim; i ++ )
    {
      REAL8 val = (2.0 * rand()/RAND_MAX) - 1.0;	/* double between [-1,1] */
      gsl_vector_set ( in_vect, i, val );
    } /* for i < dim */


  fprintf(stderr, "3: Testing gsl_vector (de)serialization...\n\n");

  fprintf(stderr, "Initial gsl_vector:\n");
  fprintf(stderr, "%s = ", LALXMLC_NAMETEST3 );
  XLALfprintfGSLvector ( stderr, "%.16g", in_vect );

  /* ---------- serialize gsl_vector into VOTable fragment */
  if ( (xmlChildNodeList = XLALgsl_vector2VOTableNode(in_vect, "vect", "m")) == NULL ) {
    fprintf(stderr, "LALXMLTest: [XLALgsl_vector2VOTableNode(): %s]\n", LALXMLC_MSGEFUN);
    return LALXMLC_EFUN;
  }
  fprintf(stderr, "LALXMLTest: [XLALgsl_vector2VOTableNode(): %s]\n", LALXMLC_MSGENOM);

  /* wrap into a dummy RESOURCE node*/
  if ( (xmlFragment = XLALCreateVOTableResourceNode(LALXMLC_TYPETEST3, LALXMLC_NAMETEST3, xmlChildNodeList)) == NULL ) {
    xmlFreeNodeList(xmlChildNodeList);
    XLALPrintError("%s: Couldn't create RESOURCE node: %s\n", fn, LALXMLC_NAMETEST3);
    return LALXMLC_EFUN;
  }


  /* convert VOTable fragment into VOTable document */
  if ( (xmlDocument = (xmlDocPtr)XLALCreateVOTableDocumentFromTree((const xmlNodePtr)xmlFragment)) == NULL ) {
    fprintf(stderr, "LALXMLTest: [XLALCreateVOTableDocumentFromTree(): %s]\n", LALXMLC_MSGEFUN);
    return LALXMLC_EFUN;
  }
  fprintf(stderr, "LALXMLTest: [XLALCreateVOTableDocumentFromTree(): %s]\n", LALXMLC_MSGENOM);

  /* convert VOTable document into XML string */
  if(!xmlDocument2String(xmlDocument, &xmlString)) {
    return LALXMLC_EFUN;
  }

  /* ---------- display serialized structure */
  fprintf(stderr, "Serialized VOTable XML:\n");
  fprintf(stderr, "----------------------------------------------------------------------\n");
  fprintf(stderr, (char*)xmlString);
  fprintf(stderr, "----------------------------------------------------------------------\n");

  /* ---------- validate XML document */
  result = validateDocument(xmlDocument);
  if(result == XLAL_SUCCESS) {
    fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGENOM);
  }
  else if(result == XLAL_FAILURE) {
    fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEVAL);
    return LALXMLC_EVAL;
  }
  else {
    fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEFUN);
    return LALXMLC_EFUN;
  }

  /* ---------- deserialize VOTable document into structure */
  if( (out_vect = XLALVOTableDoc2gsl_vectorByName(xmlDocument, LALXMLC_TYPETEST3, LALXMLC_NAMETEST3, "vect", "m")) == NULL ) {
    fprintf(stderr, "LALXMLTest: [XLALVOTableDoc2gsl_vectorByName(): %s]\n", LALXMLC_MSGEFUN);
    return LALXMLC_EFUN;
  }
  fprintf(stderr, "LALXMLTest: [XLALVOTableDoc2gsl_vectorByName(): %s]\n", LALXMLC_MSGENOM);


  fprintf(stderr, "gsl_vector parsed back from VOTable:\n");
  fprintf(stderr, "%s = ", LALXMLC_NAMETEST3 );
  XLALfprintfGSLvector ( stderr, "%.16g", out_vect );

  /* ---------- validate test results */
  for (i=0; i < dim; i ++ )
    {
      REAL8 in, out;
      in  = gsl_vector_get ( in_vect, i );
      out = gsl_vector_get ( out_vect, i );
      err = abs ( (in - out) / (0.5 * (in + out)));
      if ( err > maxerr ) maxerr = err;
      if (  err > REAL8TOL ) {
	XLALPrintError ("%s: element %d in gsl_vector '%s' differs by more (%g) than tolerance %g: in=%.16g, out=%.16g.\n\n",
			fn, i, LALXMLC_NAMETEST3, err, REAL8TOL, in, out );
	fprintf(stderr, "LALXMLTest: [XLALVOTableDoc2gsl_vectorByName(): %s]\n\n", LALXMLC_MSGEVAL);
	return LALXMLC_EVAL;
      }
    } /* for i < dim */
  fprintf (stderr, "%s: maximal relative error %g is within tolerance of %g.\n\n", fn, maxerr, REAL8TOL );

  /* free memory */
  xmlFreeDoc(xmlDocument);
  gsl_vector_free ( in_vect );
  gsl_vector_free ( out_vect );

  return LALXMLC_ENOM;

} /* test_gsl_vector() */



/* test (de-)serialization of a gsl_matrix type */
int
test_gsl_matrix(void)
{
  static const char *fn = "test_gsl_matrix()";

  int row, col, numRows, numCols;
  gsl_matrix *in_matrix, *out_matrix;
  xmlNodePtr xmlChildNodeList = NULL;

  xmlNodePtr xmlFragment = NULL;
  xmlDocPtr xmlDocument = NULL;
  xmlChar *xmlString = NULL;
  int result;
  REAL8 err, maxerr = 0;

  /* ---------- initialize input test data ---------- */

  /* pick a random number of dimensions between [2, 6] */
  srand(time(NULL));	/* pick random seed */
  numRows = 2 + (1.0 * rand() / RAND_MAX) * 4;
  numCols = 2 + (1.0 * rand() / RAND_MAX) * 4;

  if ( (in_matrix = gsl_matrix_alloc ( numRows, numCols )) == NULL ) {
    XLALPrintError ("%s: failed to gsl_matrix_alloc(%d,%d).\n\n", fn, numRows, numCols );
    return XLAL_ENOMEM;
  }
  for (row=0; row < numRows; row ++ )
    {
      for ( col=0; col < numCols; col ++ )
	{
	  REAL8 val = (2.0 * rand()/RAND_MAX) - 1.0;	/* double between [-1,1] */
	  gsl_matrix_set ( in_matrix, row, col, val );
	} /* for col < numCols */
    } /* for row < numRows */

  fprintf(stderr, "4: Testing gsl_matrix (de)serialization...\n\n");

  fprintf(stderr, "Initial gsl_matrix:\n");
  fprintf(stderr, "%s = ", LALXMLC_NAMETEST4 );
  XLALfprintfGSLmatrix ( stderr, "%.16g", in_matrix );

  /* ---------- serialize gsl_matrix into VOTable fragment */
  if ( (xmlChildNodeList = XLALgsl_matrix2VOTableNode(in_matrix, "matrix", "ms")) == NULL ) {
    fprintf(stderr, "LALXMLTest: [XLALgsl_matrix2VOTableNode(): %s]\n", LALXMLC_MSGEFUN);
    return LALXMLC_EFUN;
  }
  fprintf(stderr, "LALXMLTest: [XLALgsl_matrix2VOTableNode(): %s]\n", LALXMLC_MSGENOM);

  /* wrap into a dummy RESOURCE node*/
  if ( (xmlFragment = XLALCreateVOTableResourceNode(LALXMLC_TYPETEST4, LALXMLC_NAMETEST4, xmlChildNodeList)) == NULL ) {
    xmlFreeNodeList(xmlChildNodeList);
    XLALPrintError("%s: Couldn't create RESOURCE node: %s\n", fn, LALXMLC_NAMETEST4);
    return LALXMLC_EFUN;
  }


  /* convert VOTable fragment into VOTable document */
  if ( (xmlDocument = (xmlDocPtr)XLALCreateVOTableDocumentFromTree((const xmlNodePtr)xmlFragment)) == NULL ) {
    fprintf(stderr, "LALXMLTest: [XLALCreateVOTableDocumentFromTree(): %s]\n", LALXMLC_MSGEFUN);
    return LALXMLC_EFUN;
  }
  fprintf(stderr, "LALXMLTest: [XLALCreateVOTableDocumentFromTree(): %s]\n", LALXMLC_MSGENOM);

  /* convert VOTable document into XML string */
  if(!xmlDocument2String(xmlDocument, &xmlString)) {
    return LALXMLC_EFUN;
  }

  /* ---------- display serialized structure */
  fprintf(stderr, "Serialized VOTable XML:\n");
  fprintf(stderr, "----------------------------------------------------------------------\n");
  fprintf(stderr, (char*)xmlString);
  fprintf(stderr, "----------------------------------------------------------------------\n");

  /* ---------- validate XML document */
  result = validateDocument(xmlDocument);
  if(result == XLAL_SUCCESS) {
    fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGENOM);
  }
  else if(result == XLAL_FAILURE) {
    fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEVAL);
    return LALXMLC_EVAL;
  }
  else {
    fprintf(stderr, "LALXMLTest: [XLALValidateDocumentBy[Ex|In]ternalSchema(): %s]\n", LALXMLC_MSGEFUN);
    return LALXMLC_EFUN;
  }

  /* ---------- deserialize VOTable document into structure */
  if( (out_matrix = XLALVOTableDoc2gsl_matrixByName(xmlDocument, LALXMLC_TYPETEST4, LALXMLC_NAMETEST4, "matrix", "ms")) == NULL ) {
    fprintf(stderr, "LALXMLTest: [XLALVOTableDoc2gsl_matrixByName(): %s]\n", LALXMLC_MSGEFUN);
    return LALXMLC_EFUN;
  }
  fprintf(stderr, "LALXMLTest: [XLALVOTableDoc2gsl_matrixByName(): %s]\n", LALXMLC_MSGENOM);


  fprintf(stderr, "gsl_matrix parsed back from VOTable:\n");
  fprintf(stderr, "%s = ", LALXMLC_NAMETEST4 );
  XLALfprintfGSLmatrix ( stderr, "%.16g", out_matrix );

  /* ---------- validate test results */
  for (row=0; row < numRows; row ++ )
    {
      for ( col=0; col < numCols; col ++ )
	{
	  REAL8 in, out;
	  in  = gsl_matrix_get ( in_matrix, row, col );
	  out = gsl_matrix_get ( out_matrix, row, col );
	  err = abs ( (in - out) / (0.5 * (in + out) ) );
	  if ( err > maxerr ) maxerr = err;
	  if ( err > REAL8TOL ) {
	    XLALPrintError ("%s: element (%d,%d) in gsl_matrix '%s' differs by more (%g) than tolerance %g: in=%.16g, out=%.16g.\n\n",
			    fn, row, col, LALXMLC_NAMETEST4, err, REAL8TOL, in, out );
	    fprintf(stderr, "LALXMLTest: [XLALVOTableDoc2gsl_matrixByName(): %s]\n\n", LALXMLC_MSGEVAL);
	    return LALXMLC_EVAL;
	  }
	} /* for col < numCols */
    } /* for row < numRows */
  fprintf (stderr, "%s: maximal relative error %g is within tolerance of %g.\n\n", fn, maxerr, REAL8TOL );

  /* free memory */
  xmlFreeDoc(xmlDocument);
  gsl_matrix_free ( in_matrix );
  gsl_matrix_free ( out_matrix );

  return LALXMLC_ENOM;

} /* test_gsl_matrix() */





/* -------------------- Helper functions -------------------- */

int xmlDocument2String(const xmlDocPtr xmlDocument, xmlChar **xmlString)
{
    /* set up local variables */
    int xmlStringBufferSize = 0;

    /* prepare XML serialization (here: indentation) */
    xmlThrDefIndentTreeOutput(1);

    /* dump document to a string buffer */
    xmlDocDumpFormatMemoryEnc(xmlDocument, xmlString, &xmlStringBufferSize, "UTF-8", 1);
    if(xmlStringBufferSize <= 0) {
        fprintf(stderr, "XML document dump failed!\n");
        return XLAL_EFAILED;
    }

    /* return string size (0 in case of error) */
    return xmlStringBufferSize;

} /* xmlDocument2String() */


int xmlString2Document(const xmlChar *xmlString, xmlDocPtr *xmlDocument)
{
    /* set up local variables */
    int result = XLAL_SUCCESS;

    /* parse XML document */
    *xmlDocument = xmlReadMemory((const char*)xmlString, strlen((const char*)xmlString), NULL, "UTF-8", 0);
    if(*xmlDocument == NULL) {
        fprintf(stderr, "XML document parsing failed!\n");
        result = XLAL_EFAILED;
    }

    /* clean up */
    xmlCleanupParser();

    return result;

} /* xmlString2Document() */


int validateDocument(const xmlDocPtr xmlDocument)
{
    /* set up local variables */
    char *schemaPath = NULL;
    char schemaUrl[PATH_MAXLEN+10] = "file://";
    int result;

    /* find schema definition file */
    result = findFileInLALDataPath("VOTable-1.1.xsd", &schemaPath);

    /* validate document */
    if(result == XLAL_SUCCESS) {
        strncat(schemaUrl, schemaPath, PATH_MAXLEN);
        result = XLALValidateDocumentByExternalSchema(xmlDocument, BAD_CAST(schemaUrl));
        LALFree(schemaPath);
    }
    else {
        fprintf(stderr, "Warning: schema definition file not found! "
                        "Falling back to internal schema definition (online resource)!\n");
        result = XLALValidateDocumentByInternalSchema(xmlDocument);
    }

    return result;

} /* validateDocument() */

int findFileInLALDataPath(const char *filename, char **validatedPath)
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
        n = snprintf(absolutePath, PATH_MAXLEN, "%s/./%s", workingDir, filename);
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
        fprintf(stderr, "Can't allocate memory (%lu)!\n", strlen(dataPathEnv)+1);
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
        n = snprintf(absolutePath,
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

} /* findFileInLALDataPath() */
