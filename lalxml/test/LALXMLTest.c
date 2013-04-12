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

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdio.h>
#include <lal/LogPrintf.h>
#include <lal/LALXML.h>
#include <lal/LALXMLVOTableCommon.h>
#include <lal/LALXMLVOTableSerializers.h>
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
int testLIGOTimeGPS(void);
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
    printf( "2: Test gsl_vector (de)serialization...\n");
    printf( "======================================================================\n");
    if ( (result = test_gsl_vector()) != LALXMLC_ENOM ) {
      LogPrintf (LOG_CRITICAL, "test_gsl_vector() failed. ret = %d\n", result );
      return result;
    }

    printf( "\n" );
    printf( "======================================================================\n");
    printf( "3: Test gsl_matrix (de)serialization...\n");
    printf( "======================================================================\n");
    if ( (result = test_gsl_matrix()) != LALXMLC_ENOM ) {
      LogPrintf (LOG_CRITICAL, "test_gsl_matrix() failed. ret = %d\n", result );
      return result;
    }

    printf( "\n" );
    printf( "======================================================================\n");
    printf( "4: Test table (de)serialization...\n");
    printf( "======================================================================\n");
    if ( ( result = testTable()) != LALXMLC_ENOM ) {
      LogPrintf (LOG_CRITICAL, "testTable() failed. ret = %d\n", result );
      return result;
    }

    LALCheckMemoryLeaks();

    return LALXMLC_ENOM;

} /* main() */


int testLIGOTimeGPS(void)
{
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
      XLALPrintError ( "%s: XLALLIGOTimeGPS2VOTNode() failed.\n", __func__ );
      return LALXMLC_EFUN;
    }

    /* convert VOTable tree into XML string */
    if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
      XLALPrintError ("%s: XLALCreateVOTStringFromTree() failed.\n", __func__);
      return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* display serialized structure */
    printf( "----------------------------------------------------------------------\n");
    printf( "%s", xmlString );
    printf( "----------------------------------------------------------------------\n");

    /* ---------- parse XML string back and validate ---------- */
    /* convert VOTable string back into VOTable document */
    printf ("--> Parsing XML string into xmlDoc ... ");
    if ( (xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
      XLALPrintError( "%s: XLALXMLString2Doc() failed.\n", __func__);
      return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* validate resulting XML document */
    printf ("--> Validating xmlDoc against VOTable schema ... ");
    if ( (result = validateDocument(xmlDocument)) != XLAL_SUCCESS ) {
      XLALPrintError("%s: validateDocument failed. ret = %d\n", __func__, result );
      if ( result == XLAL_FAILURE )
        return LALXMLC_EVAL;
      else
        return LALXMLC_EFUN;
    }
    printf ("ok.\n");

    /* deserialize VOTable document back into a C-structure */
    printf ("--> Read out LIGOTimeGPSStruct ... ");
    if ( XLALVOTDoc2LIGOTimeGPSByName ( xmlDocument, LALXMLC_NAMETEST1, &timeDestination)) {
      XLALPrintError ( "%s: XLALVOTDoc2LIGOTimeGPSByName() failed. errno = %d\n", __func__, xlalErrno );
      return LALXMLC_EFUN;
    }
    printf( "ok: %s = { %d, %d }\n", LALXMLC_NAMETEST1, timeDestination.gpsSeconds, timeDestination.gpsNanoSeconds );

    /* check test results */
    if(
       timeSource.gpsSeconds != timeDestination.gpsSeconds ||
       timeSource.gpsNanoSeconds != timeDestination.gpsNanoSeconds)
      {
        XLALPrintError ( "%s: LIGOTimeGPS structure differs before and after XML serialization!\n", __func__ );
        return LALXMLC_EVAL;
    }
    else
      printf ( "--> Final struct agrees with input struct. Test PASSED.\n");


    /* clean up */
    xmlFreeDoc(xmlDocument);
    xmlFreeNode ( xmlFragment );
    XLALFree ( xmlString );

    return LALXMLC_ENOM;

} /* testLIGOTimeGPS() */


/* test (de-)serialization of a gsl_vector type */
int
test_gsl_vector(void)
{
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
    XLALPrintError ("%s: failed to gsl_vector_alloc(%d).\n\n", __func__, dim );
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
    XLALPrintError( "%s: XLALgsl_vector2VOTNode() failed. errno = %d.\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  /* wrap into a dummy RESOURCE node*/
  if ( (xmlFragment = XLALCreateVOTResourceNode(LALXMLC_TYPETEST3, LALXMLC_NAMETEST3, xmlChildNodeList)) == NULL ) {
    xmlFreeNodeList(xmlChildNodeList);
    XLALPrintError("%s: Couldn't create RESOURCE node: %s\n", __func__, LALXMLC_NAMETEST3);
    return LALXMLC_EFUN;
  }
  /* convert VOTable tree into XML string */
  if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
    return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  /* ---------- display serialized structure */
  printf( "----------------------------------------------------------------------\n");
  printf( "%s", xmlString);
  printf( "----------------------------------------------------------------------\n");


  printf ("--> Parsing XML string into xmlDoc ... ");
  if( ( xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
    XLALPrintError( "%s: XLALXMLString2Doc() failed. err = %d\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  printf ("ok.\n");


  printf ("--> Validating xmlDoc against VOTable schema ... ");
  if ( (result = validateDocument(xmlDocument)) != XLAL_SUCCESS ) {
    XLALPrintError("%s: validateDocument failed. ret = %d\n", __func__, result );
    if ( result == XLAL_FAILURE )
      return LALXMLC_EVAL;
    else
      return LALXMLC_EFUN;
  }
  printf ("ok.\n");


  printf ("--> Read out gsl_vector ... ");
  if( (out_vect = XLALVOTDoc2gsl_vectorByName(xmlDocument, LALXMLC_NAMETEST3, "vect", "m")) == NULL ) {
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
			__func__, i, LALXMLC_NAMETEST3, err, REAL8TOL, in, out );
	return LALXMLC_EVAL;
      }
    } /* for i < dim */
  printf ("--> maximal relative error %g is within tolerance of %g. PASSED.\n", maxerr, REAL8TOL );

  /* free memory */
  xmlFreeDoc(xmlDocument);
  xmlFreeNode ( xmlFragment );
  XLALFree ( xmlString );

  gsl_vector_free ( in_vect );
  gsl_vector_free ( out_vect );

  return LALXMLC_ENOM;

} /* test_gsl_vector() */



/* test (de-)serialization of a gsl_matrix type */
int
test_gsl_matrix(void)
{
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
    XLALPrintError ("%s: failed to gsl_matrix_alloc(%d,%d).\n\n", __func__, numRows, numCols );
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
    XLALPrintError( "%s: XLALgsl_matrix2VOTNode() failed. errno = %d.\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  /* wrap into a dummy RESOURCE node*/
  if ( (xmlFragment = XLALCreateVOTResourceNode(LALXMLC_TYPETEST4, LALXMLC_NAMETEST4, xmlChildNodeList)) == NULL ) {
    xmlFreeNodeList(xmlChildNodeList);
    XLALPrintError("%s: Couldn't create RESOURCE node: %s\n", __func__, LALXMLC_NAMETEST4);
    return LALXMLC_EFUN;
  }
  /* convert VOTable tree into XML string */
  if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
    return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  /* ---------- display serialized structure */
  printf( "----------------------------------------------------------------------\n");
  printf( "%s", xmlString);
  printf( "----------------------------------------------------------------------\n");


  printf ("--> Parsing XML string into xmlDoc ... ");
  if( ( xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
    XLALPrintError( "%s: XLALXMLString2Doc() failed. err = %d\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  printf ("--> Validating xmlDoc against VOTable schema ... ");
  if ( (result = validateDocument(xmlDocument)) != XLAL_SUCCESS ) {
    XLALPrintError("%s: validateDocument failed. ret = %d\n", __func__, result );
    if ( result == XLAL_FAILURE )
      return LALXMLC_EVAL;
    else
      return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  /* ---------- deserialize VOTable document into structure */
  printf ("--> Read out LIGOTimeGPSStruct ... ");
  if( (out_matrix = XLALVOTDoc2gsl_matrixByName ( xmlDocument, LALXMLC_NAMETEST4, "matrix", "ms")) == NULL ) {
    XLALPrintError ("%s: XLALVOTDoc2gsl_matrixByName() failed. errno = %d\n", __func__, xlalErrno );
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
			    __func__, row, col, LALXMLC_NAMETEST4, err, REAL8TOL, in, out );
	    return LALXMLC_EVAL;
	  }
	} /* for col < numCols */
    } /* for row < numRows */
  printf ("--> maximal relative error %g is within tolerance of %g. PASSED.\n", maxerr, REAL8TOL );

  /* free memory */
  xmlFreeDoc(xmlDocument);
  xmlFreeNode ( xmlFragment );
  XLALFree ( xmlString );

  gsl_matrix_free ( in_matrix );
  gsl_matrix_free ( out_matrix );

  return LALXMLC_ENOM;

} /* test_gsl_matrix() */


/* generic test for table- writing and parsing tools */
int testTable ( void )
{
  xmlNodePtr fieldNodeList = NULL, newFieldNode = NULL;

  xmlDocPtr xmlDocument = NULL;
  char *xmlString = NULL;
  UINT4 j;
#define IN_ROWS 3
  REAL8 FreqIn[IN_ROWS]   = { 100.1234, 101.234, 102.345 };
  REAL4 AlphaIn[IN_ROWS]  = { 0.1234, 2.123434, 3.2341 };
  REAL4 DeltaIn[IN_ROWS]  = { -1.234, -0.5, 1.234 };
  const CHAR *NameIn[IN_ROWS]   = { "Pulsar 1", "another pulsar", "PSR J0537-6910" };
  INT4  IndexIn[IN_ROWS]  = { 5, 7, 99 };
  COMPLEX8 FaIn[IN_ROWS]  = { { 1, 2 }, {3, 4}, {5, 6} };

  printf ("--> Input table values ... \n");
  printf ("FreqIn  = ");
  for ( j=0; j < IN_ROWS; j ++ ) printf ("%f%s", FreqIn[j], j<IN_ROWS-1? ", " : "\n" );
  printf ("AlphaIn = ");
  for ( j=0; j < IN_ROWS; j ++ ) printf ("%f%s", AlphaIn[j], j<IN_ROWS-1? ", " : "\n" );
  printf ("DeltaIn = ");
  for ( j=0; j < IN_ROWS; j ++ ) printf ("%f%s", DeltaIn[j], j<IN_ROWS-1? ", " : "\n" );
  printf ("NameIn  = ");
  for ( j=0; j < IN_ROWS; j ++ ) printf ("%s%s", NameIn[j], j<IN_ROWS-1? ", " : "\n" );
  printf ("IndexIn = ");
  for ( j=0; j < IN_ROWS; j ++ ) printf ("%d%s", IndexIn[j], j<IN_ROWS-1? ", " : "\n" );
  printf ("FaIn    = ");
  for ( j=0; j < IN_ROWS; j ++ ) printf ("%f %f%s", crealf(FaIn[j]), cimagf(FaIn[j]), j<IN_ROWS-1? ", " : "\n" );


  /* ---------- create FIELDS */
  /* ----- Freq */
  if ( (fieldNodeList = XLALCreateVOTFieldNode ( "Freq", "Hz", VOT_REAL8, NULL )) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTFieldNode() failed for 'Freq'. xlalErrno = %d\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  /* ----- Alpha */
  if ( (newFieldNode = XLALCreateVOTFieldNode ( "Alpha", "rad", VOT_REAL4, NULL )) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTFieldNode() failed for 'Alpha'. xlalErrno = %d\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  if ( xmlAddSibling ( fieldNodeList, newFieldNode ) == NULL ) {
    XLALPrintError ("%s: xmlAddSibling() failed.\n", __func__ );
    return LALXMLC_EFUN;
  }
  /* ----- Delta */
  if ( (newFieldNode = XLALCreateVOTFieldNode ( "Delta", "rad", VOT_REAL4, NULL )) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTFieldNode() failed for 'Delta'. xlalErrno = %d\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  if ( xmlAddSibling ( fieldNodeList, newFieldNode ) == NULL ) {
    XLALPrintError ("%s: xmlAddSibling() failed.\n", __func__ );
    return LALXMLC_EFUN;
  }
  /* ----- Names */
  if ( (newFieldNode = XLALCreateVOTFieldNode ( "Name", NULL, VOT_CHAR, "*" )) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTFieldNode() failed for 'Names'. xlalErrno = %d\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  if ( xmlAddSibling ( fieldNodeList, newFieldNode ) == NULL ) {
    XLALPrintError ("%s: xmlAddSibling() failed.\n", __func__ );
    return LALXMLC_EFUN;
  }
  /* ----- Fa */
  if ( (newFieldNode = XLALCreateVOTFieldNode ( "Fa", NULL, VOT_COMPLEX8, NULL)) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTFieldNode() failed for 'Fa'. xlalErrno = %d\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  if ( xmlAddSibling ( fieldNodeList, newFieldNode ) == NULL ) {
    XLALPrintError ("%s: xmlAddSibling() failed.\n", __func__ );
    return LALXMLC_EFUN;
  }

  /* ----- Indices */
  if ( (newFieldNode = XLALCreateVOTFieldNode ( "Index", NULL, VOT_INT4, NULL)) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTFieldNode() failed for 'Indices'. xlalErrno = %d\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  if ( xmlAddSibling ( fieldNodeList, newFieldNode ) == NULL ) {
    XLALPrintError ("%s: xmlAddSibling() failed.\n", __func__ );
    return LALXMLC_EFUN;
  }

  xmlNodePtr xmlFragment, xmlTable, xmlTabledataNode;

  if ( (xmlTabledataNode = XLALCreateVOTTabledataNode ( fieldNodeList, 3, "%.16g,%.7g,%.7g,%s,%g %g,%0d", FreqIn, AlphaIn, DeltaIn, NameIn, FaIn, IndexIn  )) == NULL ){
    XLALPrintError("%s: XLALCreateVOTTabledataNode() failed. errno = %d.\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }

  if ( (xmlTable = XLALCreateVOTTableNode ( NULL, fieldNodeList, xmlTabledataNode )) == NULL ) {
    XLALPrintError("%s: XLALCreateVOTTableNode() failed. errno = %d.\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  /* wrap into a dummy RESOURCE node*/
  if ( (xmlFragment = XLALCreateVOTResourceNode(LALXMLC_TYPETEST5, LALXMLC_NAMETEST5, xmlTable)) == NULL ) {
    XLALPrintError("%s: Couldn't create RESOURCE node: %s\n", __func__, LALXMLC_NAMETEST5);
    return LALXMLC_EFUN;
  }

  /* convert VOTable tree into XML string */
  if( (xmlString = XLALCreateVOTStringFromTree ( xmlFragment )) == NULL ) {
    XLALPrintError ("%s: XLALCreateVOTStringFromTree() failed. errno = %d.\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }

  /* ---------- display serialized structure */
  printf( "Serialized VOTable XML:\n");
  printf( "----------------------------------------------------------------------\n");
  printf( "%s", xmlString);
  printf( "----------------------------------------------------------------------\n");


  /* convert XML string back into VOTable document */
  printf ("--> Parsing XML string into xmlDoc ... ");
  if( ( xmlDocument = XLALXMLString2Doc ( xmlString )) == NULL ) {
    XLALPrintError( "%s: XLALXMLString2Doc() failed. err = %d\n", __func__, xlalErrno );
    return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  /* validate XML document */
  int result;
  printf ("--> Validating xmlDoc against VOTable schema ... ");
  if ( (result = validateDocument(xmlDocument)) != XLAL_SUCCESS ) {
    XLALPrintError("%s: validateDocument failed. ret = %d\n", __func__, result );
    if ( result == XLAL_FAILURE )
      return LALXMLC_EVAL;
    else
      return LALXMLC_EFUN;
  }
  printf ("ok.\n");

  /* ---------- de-serialize back into C-struct ---------- */
  UINT4 i, numCols, numRows, nRows;
  VOTFieldVector *fieldVect = NULL;
  REAL8 *FreqOut = NULL;
  REAL4 *AlphaOut = NULL, *DeltaOut = NULL;
  CHAR **NameOut = NULL;
  INT4 *IndexOut = NULL;
  COMPLEX8 *FaOut = NULL;

  if ( ( fieldVect = XLALReadVOTFIELDNodes ( xmlDocument, LALXMLC_NAMETEST5 )) == NULL ) {
    XLALPrintError ("%s: XLALReadVOTFIELDNodes() failed to obtain FIELD elements.\n", __func__ );
    return LALXMLC_EFUN;
  }
  numCols = fieldVect->length;
  for ( i=0; i < numCols; i ++ )
    {
      VOTField *thisField = &fieldVect->data[i];
#if 0
      printf ("i = %d: name='%s', datatype='%s'", i, thisField->name, XLALVOTDatatype2String ( thisField->datatype ) );
      if ( thisField->unit )
        printf ("  unit = '%s'", thisField->unit );
      if ( thisField -> arraysize )
        printf ("  arraysize = '%s' ", thisField->arraysize );
      printf ("\n");
#endif

      /* ----- Freq ----- */
      if ( !strcmp ( (const char*)thisField->name, "Freq" ) ) {
        if ( (FreqOut = XLALReadVOTTabledataSimpleColumn ( xmlDocument, LALXMLC_NAMETEST5, i, VOT_REAL8, &numRows )) == NULL ) {
          XLALPrintError ("%s: XLALReadVOTTableColumn() failed to read column Nr %d: 'Freq'.\n", __func__, i );
          return LALXMLC_EFUN;
        }
      }
      /* ----- Alpha ----- */
      else if ( !strcmp ( (const char*)thisField->name, "Alpha" ) ) {
        if ( (AlphaOut = XLALReadVOTTabledataSimpleColumn ( xmlDocument, LALXMLC_NAMETEST5, i, VOT_REAL4, &nRows )) == NULL ) {
          XLALPrintError ("%s: XLALReadVOTTableColumn() failed to read column Nr %d: 'Alpha'.\n", __func__, i );
          return LALXMLC_EFUN;
        }
      }
      /* ----- Delta ----- */
      else if ( !strcmp ( (const char*)thisField->name, "Delta" ) ) {
        if ( (DeltaOut = XLALReadVOTTabledataSimpleColumn ( xmlDocument, LALXMLC_NAMETEST5, i, VOT_REAL4, &nRows )) == NULL ) {
          XLALPrintError ("%s: XLALReadVOTTableColumn() failed to read column Nr %d: 'Delta'.\n", __func__, i );
          return LALXMLC_EFUN;
        }
      }
      /* ----- Name ----- */
      else if ( !strcmp ( (const char*)thisField->name, "Name" ) ) {
        if ( (NameOut = XLALReadVOTTabledataSimpleColumn ( xmlDocument, LALXMLC_NAMETEST5, i, VOT_CHAR, &nRows )) == NULL ) {
          XLALPrintError ("%s: XLALReadVOTTableColumn() failed to read column Nr %d: 'Name'.\n", __func__, i );
          return LALXMLC_EFUN;
        }
      }
      /* ----- Index ----- */
      else if ( !strcmp ( (const char*)thisField->name, "Index" ) ) {
        if ( (IndexOut = XLALReadVOTTabledataSimpleColumn ( xmlDocument, LALXMLC_NAMETEST5, i, VOT_INT4, &nRows )) == NULL ) {
          XLALPrintError ("%s: XLALReadVOTTableColumn() failed to read column Nr %d: 'Name'.\n", __func__, i );
          return LALXMLC_EFUN;
        }
      }
      /* ----- Fa ----- */
      else if ( !strcmp ( (const char*)thisField->name, "Fa" ) ) {
        if ( (FaOut = XLALReadVOTTabledataSimpleColumn ( xmlDocument, LALXMLC_NAMETEST5, i, VOT_COMPLEX8, &nRows )) == NULL ) {
          XLALPrintError ("%s: XLALReadVOTTableColumn() failed to read column Nr %d: 'Fa'.\n", __func__, i );
          return LALXMLC_EFUN;
        }
      }
      else
        {
          XLALPrintError ( "%s: Unknown table column encountered '%s'\n", __func__, (const char*)thisField->name );
          return LALXMLC_EFUN;
        }

    } /* for i < numCols */

  XLALDestroyVOTFieldVector ( fieldVect );

  if ( numRows != IN_ROWS ) {
    XLALPrintError ("%s: some table rows went missing ... input %d rows, read out %d rows\n", __func__, IN_ROWS, numRows );
    return LALXMLC_EFUN;
  }

  printf ("--> Read back table values ... \n");
  printf ("FreqOut  = ");
  for ( j=0; j < numRows; j ++ ) printf ("%f%s", FreqOut[j], j<numRows-1? ", " : "\n" );
  printf ("AlphaOut = ");
  for ( j=0; j < numRows; j ++ ) printf ("%f%s", AlphaOut[j], j<numRows-1? ", " : "\n" );
  printf ("DeltaOut = ");
  for ( j=0; j < numRows; j ++ ) printf ("%f%s", DeltaOut[j], j<numRows-1? ", " : "\n" );
  printf ("NameOut  = ");
  for ( j=0; j < numRows; j ++ ) printf ("%s%s", NameOut[j], j<numRows-1? ", " : "\n" );
  printf ("IndexOut = ");
  for ( j=0; j < numRows; j ++ ) printf ("%d%s", IndexOut[j], j<numRows-1? ", " : "\n" );
  printf ("FaOut    = ");
  for ( j=0; j < numRows; j ++ ) printf ("%f %f%s", crealf(FaOut[j]), cimagf(FaOut[j]), j<numRows-1? ", " : "\n" );

  /* compare input- and output-values */
  printf ("--> Comparing input values and those parsed back ... ");
  for ( j=0; j < numRows; j++ )
    {
      if ( gsl_fcmp (FreqIn[j], FreqOut[j], REAL8TOL ) ) {
        XLALPrintError ("Input Freq values %.16f and parsed %.16f in row %d differ by more than eps=%g\n", FreqIn[j], FreqOut[j], j, REAL8TOL);
        return LALXMLC_EFUN;
      }
      if ( gsl_fcmp (AlphaIn[j], AlphaOut[j], REAL4TOL ) ) {
        XLALPrintError ("Input Alpha values %.6f and parsed %.6f in row %d differ by more than eps=%g\n", AlphaIn[j], AlphaOut[j], j, REAL4TOL);
        return LALXMLC_EFUN;
      }
      if ( gsl_fcmp (DeltaIn[j], DeltaOut[j], REAL4TOL ) ) {
        XLALPrintError ("Input Delta values %.6f and parsed %.6f in row %d differ by more than eps=%g\n", DeltaIn[j], DeltaOut[j], j, REAL4TOL);
        return LALXMLC_EFUN;
      }
      if ( strcmp (NameIn[j], NameOut[j]) ) {
        XLALPrintError ("Input Name '%s' differs from parsed '%s' in row %d\n", NameIn[j], NameOut[j], j );
        return LALXMLC_EFUN;
      }
      if ( IndexIn[j] != IndexOut[j] ) {
        XLALPrintError ("Input Index '%d' differs from parsed '%d' in row %d\n", IndexIn[j], IndexOut[j], j );
        return LALXMLC_EFUN;
      }
      if ( gsl_fcmp (crealf(FaIn[j]), crealf(FaOut[j]), REAL4TOL ) || gsl_fcmp (cimagf(FaIn[j]), cimagf(FaOut[j]), REAL4TOL ) ) {
        XLALPrintError ("Input Fa {%.6f,%.6f} and parsed {%.6f,%.6f} in row %d differ by more than eps=%g\n",
                        crealf(FaIn[j]), cimagf(FaIn[j]), crealf(FaOut[j]), cimagf(FaOut[j]), j, REAL4TOL);
        return LALXMLC_EFUN;
      }

    } /* for j < numRows */
  printf ("OK. Everything agrees within tolerances.\n");


  XLALFree ( FreqOut );
  XLALFree ( AlphaOut );
  XLALFree ( DeltaOut );
  XLALFree ( IndexOut );
  for ( j=0; j < numRows; j ++ ) XLALFree ( NameOut[j] );
  XLALFree ( NameOut );
  XLALFree ( FaOut );

  /* clean up */
  xmlFreeDoc ( xmlDocument );
  xmlFreeNode ( xmlFragment );
  XLALFree ( xmlString );

  return LALXMLC_ENOM;

} /* testTableCreation() */


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
