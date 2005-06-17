/*
 * Copyright (C) 2005 Reinhard Prix
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

/**
 * \author Reinhard Prix
 * \date 2005
 * \file 
 * \ingroup moduleLatticeCovering
 * \brief test-functions for the module LatticeCovering
 *
 * $Id$
 *
 */


/*---------- INCLUDES ----------*/
#include <math.h>

#include "LatticeCovering.h"

#include <lalapps.h>

RCSID ("$Id$");

/*---------- DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)


/* prototypes */
int main(void);
void testGS(void);
void testCovering(void);

static int writeREAL8VectorList (FILE *fp, const REAL8VectorList *list);
static int print_matrix (const gsl_matrix *m);
static int print_vector (const gsl_vector *v);


/*--------------------------------------------------*/
/* Test function(s) */
/*--------------------------------------------------*/
int main (void)
{

  lalDebugLevel = 1;

  /* for production code: dont' let gsl abort on errors, but return error-codes to caller! */
  /* gsl_set_error_handler_off (); */

  /*  testGS(); */

  testCovering();

  return 0;
}

BOOLEAN testArea1 ( const REAL8Vector *point ); /* example boundary-condition */

void
testCovering (void)
{
  LALStatus lstat = blank_status;
  REAL8VectorList *covering = NULL;
  FILE *fp;
  gsl_matrix *generatorA = NULL;
  gsl_matrix *generator0 = NULL;
  gsl_matrix *generatorB = NULL;
  gsl_matrix_view gij;
  REAL8Vector startPoint;

  LatticeType type = LATTICE_TYPE_ANSTAR;
  REAL8 radius = 0.1;
  double gij_data[] = { 1, 0.1, 0.2,
			0.1, 2, 0.5,
			0.2, 0.5, 3 };
  UINT4 dim = 3;


  /* first method: do steps 'by hand' */
  XLALGetLatticeGenerator (&generator0, dim, type);
  gsl_matrix_scale ( generator0, radius );
  XLALReduceGenerator2FullRank( &generatorA, generator0 );

  /* second: use new function */
  gij  = gsl_matrix_view_array ( gij_data, dim, dim );

  XLALFindCoveringGenerator (&generatorB, type, dim, radius, &(gij.matrix) );
  if ( xlalErrno )
    {
      int code = xlalErrno;
      XLALClearErrno(); 
      LALPrintError ("\nERROR: XLALFindCoveringGenerator() failed (xlalErrno = %d)!\n\n", code);
      return;
    }
  

  printf ("First generator: \n");
  print_matrix ( generatorA );

  printf ("Second generator: \n");
  print_matrix ( generatorB );


  startPoint.length = dim;
  startPoint.data = LALCalloc (dim, sizeof(startPoint.data[0]) ); /* already (0,0) */

  LAL_CALL( LALLatticeFill (&lstat, &covering, generatorB, &startPoint, testArea1), &lstat);

  if ( (fp = fopen ( "test_lattice1.dat", "wb" )) == NULL )
    {
      LALPrintError ("\nFailed to open 'test_lattice1.dat' for writing!\n\n");
      return;
    }
  writeREAL8VectorList (fp, covering);
  fclose(fp);
  REAL8VectorListDestroy ( covering );
  covering = NULL;

  /* do the same again with the new high-level function */
  LAL_CALL( LALLatticeCovering (&lstat, &covering, radius, &(gij.matrix), &startPoint, testArea1), &lstat);
  if ( (fp = fopen ( "test_lattice2.dat", "wb" )) == NULL )
    {
      LALPrintError ("\nFailed to open 'test_lattice2.dat' for writing!\n\n");
      return;
    }
  writeREAL8VectorList (fp, covering);
  fclose(fp);
  REAL8VectorListDestroy ( covering );
  covering = NULL;

  /* free memory */
  gsl_matrix_free ( generator0 );
  gsl_matrix_free ( generatorA );
  gsl_matrix_free ( generatorB );
  LALFree ( startPoint.data );

  /* check all (LAL-)memory for leaks... (unfortunately doesn't cover GSL!) */
  LALCheckMemoryLeaks(); 
  
} /* testCovering() */

/* test boundary-conditions: [-1, 1]^N */
BOOLEAN
testArea1 ( const REAL8Vector *point )
{
  UINT4 i;

  if ( !point || !point->data)
    return FALSE;

  for ( i=0; i < point->length; i++ )
    {
      if ( fabs( point->data[i] ) > 1.0 )
	return FALSE;
    }
  
  return TRUE;
} /* testArea1() */


void
testGS(void)
{
  gsl_matrix_view m1, gij;
  gsl_matrix *orto = NULL;
  gsl_matrix *sGM = NULL;
  gsl_matrix *testM = NULL;

  REAL8 m1data[] = { 1, 	-1, 	0,
		    -2.0/3, 	1.0/3, 	1.0/3};
  
  REAL8 gijdata[] = { 1, 0, 0,
		      0, 1, 0,
		      0, 0, 1  };

  m1 = gsl_matrix_view_array ( m1data, 2, 3 );

  gij = gsl_matrix_view_array ( gijdata, 3, 3 );

  printf ("\nMetric:\n");
  print_matrix (&(gij.matrix));

  printf ("\nInput-matrix:\n");
  print_matrix (&(m1.matrix));

  printf ("\nResulting orthogonal matrix: \n");
  print_matrix (orto);


  XLALReduceGenerator2FullRank( &sGM, &(m1.matrix) );

  printf ("Square generating matrix:\n");
  print_matrix (sGM);

  XLALGetLatticeGenerator (&testM, 4, LATTICE_TYPE_ANSTAR);
  printf ("produced Generating Matrix: \n");
  print_matrix (testM);

} /* testGS() */


int
print_matrix (const gsl_matrix *m)
{
  UINT4 i,j;

  if ( m == NULL ) {
    printf ("\nERROR: print_matrix called with NULL-matrix \n");
    return -1;
  }

  printf ("[");
  for (i=0; i < m->size1; i++)
    {
      for (j=0; j < m->size2; j++)
	{
	  printf ("\t%9.6g", gsl_matrix_get (m, i, j) );
	}
      printf (";\n");
    }

  printf ("]\n");
  return 0;
} /* print_matrix() */


int
print_vector (const gsl_vector *v)
{
  UINT4 i;

  if ( v == NULL ) {
    printf ("\nERROR: print_vector called with NULL-vector \n");
    return -1;
  }

  printf ("[");
  for (i=0; i < v->size; i++)
    {
      printf ("%9.6g ", gsl_vector_get (v, i));
    }
  printf (" ]\n");

  return 0;
} /* print_vector() */


/* write a REAL8VectorList to a file */
int
writeREAL8VectorList (FILE *fp, const REAL8VectorList *list)
{
  UINT4 i;
  const REAL8VectorList *ptr;

  if ( !list || !fp )
    return -1;

  ptr = list;
  do 
    {
      for (i=0; i < ptr->entry.length; i++ )
	fprintf (fp, "%g ", ptr->entry.data[i] );

      fprintf (fp, "\n");

    } while ( (ptr = ptr->next) != NULL );

  return 0;
} /* writeREAL8VectorList() */
