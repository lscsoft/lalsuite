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

#define LALXMLC_TYPETEST5 "table_wrapper"
#define LALXMLC_NAMETEST5 "testTable"



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
int testTable ( void );

/* private utility prototypes */
int validateDocument(const xmlDocPtr xmlDocument);
int findFileInLALDataPath(const char *filename, char **validatedPath);

int main(void)
{
    /* set up local variables */
    int result = LALXMLC_ENOM;

    /* set debug level*/
    lalDebugLevel = LALMSGLVL3;

    printf( "======================================================================\n");
    printf( "1: Test LIGOTimeGPS (de)serialization\n");
    printf( "======================================================================\n");
    if ( (result = testLIGOTimeGPS()) != LALXMLC_ENOM ) {
      LogPrintf (LOG_CRITICAL, "testLIGOTimeGPS() failed. ret = %d\n", result );
      return result;
    }

    printf( "\n" );
    printf( "======================================================================\n");
    printf( "2: Test PulsarDopplerParams (de)serialization...\n");
    printf( "======================================================================\n");
    result = testPulsarDopplerParams();
    if(result != LALXMLC_ENOM) {
      LogPrintf (LOG_CRITICAL, "testPulsarDopplerParams() failed. ret = %d\n", result );
      return result;
    }

    printf( "\n" );
    printf( "======================================================================\n");
    printf( "3: Test gsl_vector (de)serialization...\n");
    printf( "======================================================================\n");
    if ( (result = test_gsl_vector()) != LALXMLC_ENOM ) {
      LogPrintf (LOG_CRITICAL, "test_gsl_vector() failed. ret = %d\n", result );
      return result;
    }

    printf( "\n" );
    printf( "======================================================================\n");
    printf( "4: Test gsl_matrix (de)serialization...\n");
    printf( "======================================================================\n");
    if ( (result = test_gsl_matrix()) != LALXMLC_ENOM ) {
      LogPrintf (LOG_CRITICAL, "test_gsl_matrix() failed. ret = %d\n", result );
      return result;
    }

    printf( "\n" );
    printf( "======================================================================\n");
    printf( "5: Test table (de)serialization...\n");
    printf( "======================================================================\n");
    if ( ( result = testTable()) != LALXMLC_ENOM ) {
      LogPrintf (LOG_CRITICAL, "testTable() failed. ret = %d\n", result );
      return result;
    }

    return LALXMLC_ENOM;

} /* main() */


int testLIGOTimeGPS(void)
{
    static const char *fn = "testLIGOTimeGPS()";

    static LIGOTimeGPS timeSource;
    static LIGOTimeGPS timeDestination;
    xmlNodePtr xmlFragment = NULL;
    xmlDocPtr xmlDocument = NULL;
    char *xmlString = NULL;
    int result;

    /* initialize test data */
    timeSource.gpsSeconds = 15;
    timeSource.gpsNanoSeconds = 200;

    printf( "--> Initial LIGOTimeGPS struct: ");
    printf( "%s = { %d, %d }\n", LALXMLC_NAMETEST1, timeSource.gpsSeconds, timeSource.gpsNanoSeconds );

    /* serialize structure into VOTable fragment */
    printf( "--> Serializing into XML string ... ");
    if ( (xmlFragment = XLALLIGOTimeGPS2VOTNode ( &timeSource, LALXMLC_NAMETEST1)) == NULL ) {
      XLALPrintError ( "%s: XLALLIGOTimeGPS2VOTNode() failed.\n", fn );
      return LALXMLC_EFUN;
    }
    /* convert VOTable tree into XML string */
    if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
      XLALPrintError ("%s: XLALCreateVOTStringFromTree() failed.\n", fn);
      return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* display serialized structure */
    printf( "----------------------------------------------------------------------\n");
    printf( xmlString);
    printf( "----------------------------------------------------------------------\n");

    /* ---------- pase XML string back and validate ---------- */
    /* convert VOTable string back into VOTable document */
    printf ("--> Parsing XML string into xmlDoc ... ");
    if ( (xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
      XLALPrintError( "%s: XLALXMLString2Doc() failed.\n", fn);
      return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* validate resulting XML document */
    printf ("--> Validating xmlDoc against VOTable schema ... ");
    if ( (result = validateDocument(xmlDocument)) != XLAL_SUCCESS ) {
      XLALPrintError("%s: validateDocument failed. ret = %d\n", fn, result );
      if ( result == XLAL_FAILURE )
        return LALXMLC_EVAL;
      else
        return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* deserialize VOTable document back into a C-structure */
    printf ("--> Read out LIGOTimeGPSStruct ... ");
    if ( XLALVOTDoc2LIGOTimeGPSByName ( xmlDocument, LALXMLC_NAMETEST1, &timeDestination)) {
      XLALPrintError ( "%s: XLALVOTDoc2LIGOTimeGPSByName() failed. errno = %d\n", fn, xlalErrno );
      return LALXMLC_EFUN;
    }
    /* clean up */
    xmlFreeDoc(xmlDocument);
    printf( "ok: %s = { %d, %d }\n", LALXMLC_NAMETEST1, timeDestination.gpsSeconds, timeDestination.gpsNanoSeconds );

    /* check test results */
    if(
       timeSource.gpsSeconds != timeDestination.gpsSeconds ||
       timeSource.gpsNanoSeconds != timeDestination.gpsNanoSeconds)
      {
        XLALPrintError ( "%s: LIGOTimeGPS structure differs before and after XML serialization!\n", fn );
        return LALXMLC_EVAL;
    }
    else
      printf ( "--> Final struct agrees with input struct. Test PASSED.\n");

    return LALXMLC_ENOM;

} /* testLIGOTimeGPS() */


int testPulsarDopplerParams(void)
{
    static const char *fn = "testPulsarDopplerParams()";

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
      XLALPrintError( "%s: XLALPulsarDopplerParams2VOTNode() failed. err = %d\n", fn, xlalErrno );
      return LALXMLC_EFUN;
    }
    /* convert VOTable tree into XML string */
    if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
      XLALPrintError( "%s: XLALCreateVOTStringFromTree() failed. err = %d\n", fn, xlalErrno );
      return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* display serialized structure */
    printf( "----------------------------------------------------------------------\n");
    printf( xmlString );
    printf( "----------------------------------------------------------------------\n");

    /* convert XML string back into VOTable document */
    printf ("--> Parsing XML string into xmlDoc ... ");
    if( ( xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
      XLALPrintError( "%s: XLALXMLString2Doc() failed. err = %d\n", fn, xlalErrno );
      return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* validate XML document */
    printf ("--> Validating xmlDoc against VOTable schema ... ");
    if ( (result = validateDocument(xmlDocument)) != XLAL_SUCCESS ) {
      XLALPrintError("%s: validateDocument failed. ret = %d\n", fn, result );
      if ( result == XLAL_FAILURE )
        return LALXMLC_EVAL;
      else
        return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* parse VOTable document back into C-structure */
    printf ("--> Read out PulsarDopplerParams struct ... ");
    if ( XLALVOTDoc2PulsarDopplerParamsByName ( xmlDocument, LALXMLC_NAMETEST2, &pdpDestination)) {
      XLALPrintError ( "%s: XLALVOTDoc2LIGOTimeGPSByName() failed. errno = %d\n", fn, xlalErrno );
      return LALXMLC_EFUN;
    }
    /* clean up */
    xmlFreeDoc(xmlDocument);

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
      XLALPrintError ( "%s: PulsarDopplerParams structure differs before and after XML serialization!\n", fn);
      return LALXMLC_EVAL;
    }

    printf ( "--> Final struct agrees with input struct. Test PASSED.\n");
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
  char *xmlString = NULL;
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


  printf("--> Initial gsl_vector: ");
  printf("%s = ", LALXMLC_NAMETEST3 );
  XLALfprintfGSLvector ( stdout, "%.16g", in_vect );

  /* ---------- serialize gsl_vector into VOTable fragment */
  printf( "--> Serializing into XML string ... ");
  if ( (xmlChildNodeList = XLALgsl_vector2VOTNode(in_vect, "vect", "m")) == NULL ) {
    XLALPrintError( "%s: XLALgsl_vector2VOTNode() failed. errno = %d.\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }
  /* wrap into a dummy RESOURCE node*/
  if ( (xmlFragment = XLALCreateVOTResourceNode(LALXMLC_TYPETEST3, LALXMLC_NAMETEST3, xmlChildNodeList)) == NULL ) {
    xmlFreeNodeList(xmlChildNodeList);
    XLALPrintError("%s: Couldn't create RESOURCE node: %s\n", fn, LALXMLC_NAMETEST3);
    return LALXMLC_EFUN;
  }
  /* convert VOTable tree into XML string */
  if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
    return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  /* ---------- display serialized structure */
  printf( "----------------------------------------------------------------------\n");
  printf( xmlString);
  printf( "----------------------------------------------------------------------\n");


  printf ("--> Parsing XML string into xmlDoc ... ");
  if( ( xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
    XLALPrintError( "%s: XLALXMLString2Doc() failed. err = %d\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }
  printf ("ok.\n");


  printf ("--> Validating xmlDoc against VOTable schema ... ");
  if ( (result = validateDocument(xmlDocument)) != XLAL_SUCCESS ) {
    XLALPrintError("%s: validateDocument failed. ret = %d\n", fn, result );
    if ( result == XLAL_FAILURE )
      return LALXMLC_EVAL;
    else
      return LALXMLC_EFUN;
  }
  printf ("ok.\n");


  printf ("--> Read out gsl_vector ... ");
  if( (out_vect = XLALVOTDoc2gsl_vectorByName(xmlDocument, LALXMLC_TYPETEST3, LALXMLC_NAMETEST3, "vect", "m")) == NULL ) {
    return LALXMLC_EFUN;
  }
  printf( "ok : %s = ", LALXMLC_NAMETEST3 );
  XLALfprintfGSLvector ( stdout, "%.16g", out_vect );

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
	return LALXMLC_EVAL;
      }
    } /* for i < dim */
  printf ("--> maximal relative error %g is within tolerance of %g. PASSED.\n", maxerr, REAL8TOL );

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
  char *xmlString = NULL;
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

  printf("--> Initial gsl_matrix: ");
  printf( "%s = ", LALXMLC_NAMETEST4 );
  XLALfprintfGSLmatrix ( stdout, "%.16g", in_matrix );

  /* ---------- serialize gsl_matrix into VOTable fragment */
  printf( "--> Serializing into XML string ... ");
  if ( (xmlChildNodeList = XLALgsl_matrix2VOTNode(in_matrix, "matrix", "ms")) == NULL ) {
    XLALPrintError( "%s: XLALgsl_matrix2VOTNode() failed. errno = %d.\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }
  /* wrap into a dummy RESOURCE node*/
  if ( (xmlFragment = XLALCreateVOTResourceNode(LALXMLC_TYPETEST4, LALXMLC_NAMETEST4, xmlChildNodeList)) == NULL ) {
    xmlFreeNodeList(xmlChildNodeList);
    XLALPrintError("%s: Couldn't create RESOURCE node: %s\n", fn, LALXMLC_NAMETEST4);
    return LALXMLC_EFUN;
  }
  /* convert VOTable tree into XML string */
  if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
    return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  /* ---------- display serialized structure */
  printf( "----------------------------------------------------------------------\n");
  printf( xmlString);
  printf( "----------------------------------------------------------------------\n");


  printf ("--> Parsing XML string into xmlDoc ... ");
  if( ( xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
    XLALPrintError( "%s: XLALXMLString2Doc() failed. err = %d\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  printf ("--> Validating xmlDoc against VOTable schema ... ");
  if ( (result = validateDocument(xmlDocument)) != XLAL_SUCCESS ) {
    XLALPrintError("%s: validateDocument failed. ret = %d\n", fn, result );
    if ( result == XLAL_FAILURE )
      return LALXMLC_EVAL;
    else
      return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  /* ---------- deserialize VOTable document into structure */
  printf ("--> Read out LIGOTimeGPSStruct ... ");
  if( (out_matrix = XLALVOTDoc2gsl_matrixByName ( xmlDocument, LALXMLC_TYPETEST4, LALXMLC_NAMETEST4, "matrix", "ms")) == NULL ) {
    XLALPrintError ("%s: XLALVOTDoc2gsl_matrixByName() failed. errno = %d\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }
  printf( "ok: %s = ", LALXMLC_NAMETEST4 );
  XLALfprintfGSLmatrix ( stdout, "%.16g", out_matrix );

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
	    return LALXMLC_EVAL;
	  }
	} /* for col < numCols */
    } /* for row < numRows */
  printf ("--> maximal relative error %g is within tolerance of %g. PASSED.\n", maxerr, REAL8TOL );

  /* free memory */
  xmlFreeDoc(xmlDocument);
  gsl_matrix_free ( in_matrix );
  gsl_matrix_free ( out_matrix );

  return LALXMLC_ENOM;

} /* test_gsl_matrix() */



int testTable ( void )
{
  const char *fn = "testTable()";

  xmlNodePtr fieldNodeList = NULL, newFieldNode = NULL;

  xmlDocPtr xmlDocument = NULL;
  char *xmlString = NULL;

  if ( (fieldNodeList = XLALCreateVOTFieldNode ( "Freq", "Hz", VOT_REAL8, NULL )) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTFieldNode() failed for 'Freq'. xlalErrno = %d\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }
  if ( (newFieldNode = XLALCreateVOTFieldNode ( "Alpha", "rad", VOT_REAL8, NULL )) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTFieldNode() failed for 'Alpha'. xlalErrno = %d\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }
  if ( xmlAddSibling ( fieldNodeList, newFieldNode ) == NULL ) {
    XLALPrintError ("%s: xmlAddSibling() failed.\n", fn );
    return LALXMLC_EFUN;
  }
  if ( (newFieldNode = XLALCreateVOTFieldNode ( "Delta", "rad", VOT_REAL8, NULL )) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTFieldNode() failed for 'Delta'. xlalErrno = %d\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }
  if ( xmlAddSibling ( fieldNodeList, newFieldNode ) == NULL ) {
    XLALPrintError ("%s: xmlAddSibling() failed.\n", fn );
    return LALXMLC_EFUN;
  }

  xmlNodePtr xmlFragment, xmlTable, xmlTabledataNode;

  REAL8 Freqs[] = { 100.1234, 101.234, 102.345 };
  REAL8 Alphas[] = { 0.1234, 2.123434, 3.2341 };
  REAL8 Deltas[] = { -1.234, -0.5, 1.234 };

  if ( (xmlTabledataNode = XLALCreateVOTTabledataNode ( fieldNodeList, 3, "%.5f,%.1f,%.2f", Freqs, Alphas, Deltas )) == NULL ){
    XLALPrintError("%s: XLALCreateVOTTabledataNode() failed. errno = %d.\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }

  if ( (xmlTable = XLALCreateVOTTableNode ( "testTable", fieldNodeList, xmlTabledataNode )) == NULL ) {
    XLALPrintError("%s: XLALCreateVOTTableNode() failed. errno = %d.\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }
  /* wrap into a dummy RESOURCE node*/
  if ( (xmlFragment = XLALCreateVOTResourceNode(LALXMLC_TYPETEST5, LALXMLC_NAMETEST5, xmlTable)) == NULL ) {
    XLALPrintError("%s: Couldn't create RESOURCE node: %s\n", fn, LALXMLC_NAMETEST5);
    return LALXMLC_EFUN;
  }

  /* convert VOTable tree into XML string */
  if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTStringFromTree() failed. errno = %d.\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }

  /* ---------- display serialized structure */
  printf( "Serialized VOTable XML:\n");
  printf( "----------------------------------------------------------------------\n");
  printf( (char*)xmlString);
  printf( "----------------------------------------------------------------------\n");


  /* convert XML string back into VOTable document */
  printf ("--> Parsing XML string into xmlDoc ... ");
  if( ( xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
    XLALPrintError( "%s: XLALXMLString2Doc() failed. err = %d\n", fn, xlalErrno );
    return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  /* validate XML document */
  int result;
  printf ("--> Validating xmlDoc against VOTable schema ... ");
  if ( (result = validateDocument(xmlDocument)) != XLAL_SUCCESS ) {
    XLALPrintError("%s: validateDocument failed. ret = %d\n", fn, result );
    if ( result == XLAL_FAILURE )
      return LALXMLC_EVAL;
    else
      return LALXMLC_EFUN;
  }
  printf ("ok.\n");




  XLALFree ( xmlString );

  return LALXMLC_ENOM;

} /* testTableCreation() */


/* -------------------- Helper functions -------------------- */

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
            fprintf(stderr, "File (%s) not found!\n", filename);
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
      fprintf(stderr, "Can't allocate memory (%lu)!\n", (long unsigned int)(strlen(dataPathEnv)+1) );
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
            /* this directory is empty, so we skip it */
          currentDataPath = nextDataPath;
          continue;
        }

        /* make sure we got an absolute path (required by "file" URI scheme) */
        if ( currentDataPath[0] == '/' )
          n = snprintf(absolutePath, PATH_MAXLEN, "%s/%s", currentDataPath, filename);
        else
          n = snprintf(absolutePath, PATH_MAXLEN, "%s/%s/%s", workingDir, currentDataPath, filename);

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
    fprintf(stderr, "File (%s) not found!\n", filename);
    return XLAL_FAILURE;

} /* findFileInLALDataPath() */
