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
 * \ingroup pulsarApps
 * \file
 * \ingroup moduleLatticeCovering
 * \brief test-functions for the module LatticeCovering
 *
 */


/*---------- INCLUDES ----------*/
#include <math.h>
#include <stdlib.h>

#include <lal/UserInput.h>

#include <lalapps.h>

#include <lal/LatticeCovering.h>
#include <lal/DopplerScan.h>

/*---------- DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

#define myrand() (1.0*rand()/RAND_MAX)

/*---------- User variables ---------- */
BOOLEAN uvar_help;
INT4 uvar_dimension;
INT4 uvar_latticeType;
REAL8 uvar_coveringRadius;
INT4 uvar_numMCpoints;
INT4 uvar_metricType;

/*---------- global vars ---------- */
extern int vrbflg;


/* non-unit metrics to use in various dimensions */
#define MAX_DIMENSION	4
REAL8 nonUnitMetric[MAX_DIMENSION][ MAX_DIMENSION * MAX_DIMENSION] = 
  {
    { /* 2D */ 
      1,   0.4,
      0.4, 0.5 
    },
    { /* 3D */
      1.0, 0.1, 0.2,
      0.1, 2.0, 0.5,
      0.2, 0.5, 3.0 
    },
    { /* 4D */
      1.0, 0.4, 0.2, -0.1, 
      0.4, 2.0, 0.3, 0.1,
      0.2, 0.3, 1.5, -0.3,
      -0.1,0.1, -0.3, 1.0
    } 

  };



/* ---------- prototypes ---------- */
int main(int argc, char *argv[]);

void testGS(void);
void testCovering(void);
BOOLEAN testArea1 ( const REAL8Vector *point ); /* example boundary-condition */

int testCovering2(void);
void initUserVars (LALStatus *status);
static int writeREAL8VectorList (FILE *fp, const REAL8VectorList *list);
static int print_matrix (const gsl_matrix *m);

/* static int print_vector (const gsl_vector *v); */

int plot2DCovering (FILE *fp, const REAL8VectorList *list, const REAL8Vector *metric, REAL8 mismatch);

int MC_MeasureCovering(FILE *fp, const REAL8VectorList *covering, 
		       const REAL8Vector *metric, UINT4 numPoints);
int ChooseRandomPoint ( REAL8Vector *point );
REAL8 FindMinDistance (const REAL8Vector *samplePoint, const REAL8VectorList *covering, 
		       const gsl_matrix *metric );

/*--------------------------------------------------*/
/* Test function(s) */
/*--------------------------------------------------*/
int main(int argc, char *argv[])
{
  LALStatus status = blank_status;  

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register user-variables */
  LAL_CALL (initUserVars (&status), &status);	  

  /* read cmdline & cfgfile  */	
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);  

  /* for production code: dont' let gsl abort on errors, but return error-codes to caller! */
  /* gsl_set_error_handler_off (); */
  if ( uvar_help )
    return 0;

  /*  testGS(); */

  /* testCovering(); */

  testCovering2();

  return 0;
}

void
initUserVars (LALStatus *status)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR ( status );

  uvar_dimension = 2;
  uvar_latticeType = LATTICE_TYPE_ANSTAR;
  uvar_coveringRadius = 0.1;
  uvar_numMCpoints = 0;
  uvar_metricType = 0;

  LALregBOOLUserVar(status,	help,		'h', UVAR_HELP, "Print this help/usage message");
  LALregINTUserVar(status,	dimension, 	'd', UVAR_OPTIONAL, "Number of dimensions to cover");
  LALregINTUserVar(status,	latticeType, 	'l', UVAR_OPTIONAL, "Lattice-type: 0=An*, 1=cubic");
  LALregINTUserVar(status,	metricType, 	'M', UVAR_OPTIONAL, "Metric use: 0=unity, 1=non-unity");
  LALregREALUserVar(status,	coveringRadius,	'R', UVAR_OPTIONAL, "Covering-radius (distance)");
  LALregINTUserVar(status,	numMCpoints,      0, UVAR_OPTIONAL, "Nr of points for MC-sampling");

  DETATCHSTATUSPTR ( status );
  RETURN(status);

} /* initUserVars() */

/* generate a covering for given metric and dimensions and write it to disk */
int
testCovering2 (void)
{
  LALStatus status = blank_status;
  REAL8VectorList *covering = NULL;
  REAL8Vector *metric = NULL;
  FILE *fp;
  REAL8Vector startPoint;
  const CHAR *fname = NULL;
  gsl_matrix_view tmp;
  gsl_matrix *gij;
  CHAR fnameRoot[4096];

  if ( uvar_metricType > 0 && uvar_dimension > MAX_DIMENSION ) {
    printf ("No non-unity metric available for dim > %d.\n", MAX_DIMENSION );
    return -1;
  }

  /* prepare root for output-filenames */
  sprintf (fnameRoot, "%s_%dD_%s_R%.2f",
	   (uvar_latticeType==0) ? "Ans" : "cubic", 
	   uvar_dimension, 
	   (uvar_metricType == 0 )? "Unity" : "NonUnity",
	   uvar_coveringRadius
	   );
   
  if ( (gij = gsl_matrix_calloc (uvar_dimension, uvar_dimension)) == NULL )
    return -1;

  if ( uvar_metricType == 0 )
    gsl_matrix_set_identity (gij);
  else
    {
      /* get matrix-view of the non-unity metric */
      tmp  = gsl_matrix_view_array ( nonUnitMetric[uvar_dimension-2], uvar_dimension, uvar_dimension );
      gsl_matrix_memcpy ( gij, &(tmp.matrix) );
    } /* if metricType == 1 */


  if ( NULL == (metric = XLALgsl2LALmetric ( gij )) )
    {
      XLALPrintError ("\nFailed to get proper metric set up \n\n");
      return -1;
    }

  /* setup start-point */
  startPoint.length = uvar_dimension;
  startPoint.data = LALCalloc (uvar_dimension, sizeof(startPoint.data[0]) ); /* already (0,0) */

  /* generate covering */
  LAL_CALL( LALLatticeCovering(&status, &covering, uvar_coveringRadius, gij, &startPoint, testArea1, 
			       uvar_latticeType), &status);


  /* output result into a file */
  {
    CHAR buf[4096];
    sprintf (buf, "%s_lattice.dat", fnameRoot);

    if ( lalDebugLevel ) 
      printf ("Outputting lattice into '%s'\n", buf );

    if ( (fp = fopen ( buf, "wb" )) == NULL ) {
      XLALPrintError ("\nFailed to open '%s' for writing!\n\n", buf );
      return -1;
    }
    if ( writeREAL8VectorList( fp, covering ) < 0 )
      return -1;
    fclose(fp);      

  } /* output lattice */

  
  /* if 2D-lattice AND debugLevel > 0, we output the 2D lattice
   * plus the metric ellipses 
   */
  if ( (uvar_dimension == 2) && (lalDebugLevel > 0) )
    {
      CHAR buf[4096];
      sprintf (buf, "%s_EllipsePlot.dat", fnameRoot);

      printf ("Plotting 2D lattice + ellipses into '%s' ... ", buf);
      if ( (fp = fopen (buf, "wb")) == NULL ) {
	printf ("Failed to open '%s' for writing!\n\n", fname);
	return -1;
      }
      plot2DCovering (fp, covering, metric, uvar_coveringRadius * uvar_coveringRadius );
      fclose(fp);
      printf ("done.\n\n");
    }


  /* produce a histogram of minimal distances of randomly picked points in Rn */
  if ( uvar_numMCpoints )
    {
      CHAR buf[4096];
      sprintf (buf, "%s_MCcovering.dat", fnameRoot);

      if ( lalDebugLevel ) 
	printf ("Outputting MCcovering results into '%s'\n", buf );

      if ( (fp = fopen ( buf, "wb" )) == NULL ) {
	XLALPrintError ("\nFailed to open '%s' for writing!\n\n", buf );
	return -1;
      }
      
      MC_MeasureCovering(fp, covering, metric, uvar_numMCpoints);

      fclose(fp);      
    } /* if uvar_numMCpoints */

  /* free memory */
  gsl_matrix_free (gij);
  XLALREAL8VectorListDestroy ( covering );
  LALFree ( startPoint.data );
  LAL_CALL ( LALDDestroyVector (&status, &metric), &status);

  LAL_CALL ( LALDestroyUserVars(&status), &status);

  /* check all (LAL-)memory for leaks... (unfortunately doesn't cover GSL!) */
  LALCheckMemoryLeaks(); 

  return 0;

} /* testCovering2() */


/* test boundary-conditions: [-1, 1]^N */
BOOLEAN
testArea1 ( const REAL8Vector *point )
{
  UINT4 i;

  if ( !point || !point->data)
    return FALSE;

  for ( i=0; i < point->length; i++ )
    {
      if ( fabs( point->data[i] ) >= 1.0 )
	return FALSE;
    }
  
  return TRUE;
} /* testArea1() */




void
testCovering (void)
{
  LALStatus status = blank_status;
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

  XLALFindCoveringGenerator (&generatorB, type, radius, &(gij.matrix) );
  if ( xlalErrno )
    {
      int code = xlalErrno;
      XLALClearErrno(); 
      XLALPrintError ("\nERROR: XLALFindCoveringGenerator() failed (xlalErrno = %d)!\n\n", code);
      return;
    }

  printf ("First generator: \n");
  print_matrix ( generatorA );

  printf ("Second generator: \n");
  print_matrix ( generatorB );


  startPoint.length = dim;
  startPoint.data = LALCalloc (dim, sizeof(startPoint.data[0]) ); /* already (0,0) */

  LAL_CALL( LALLatticeFill (&status, &covering, generatorB, &startPoint, testArea1), &status);

  if ( (fp = fopen ( "test_lattice1.dat", "wb" )) == NULL )
    {
      XLALPrintError ("\nFailed to open 'test_lattice1.dat' for writing!\n\n");
      return;
    }
  writeREAL8VectorList (fp, covering);
  fclose(fp);
  XLALREAL8VectorListDestroy ( covering );
  covering = NULL;

  /* do the same again with the new high-level function */
  LAL_CALL( LALLatticeCovering (&status, &covering, radius, &(gij.matrix), &startPoint, testArea1,
				LATTICE_TYPE_ANSTAR), &status);

  if ( (fp = fopen ( "test_lattice2.dat", "wb" )) == NULL )
    {
      XLALPrintError ("\nFailed to open 'test_lattice2.dat' for writing!\n\n");
      return;
    }
  writeREAL8VectorList (fp, covering);
  fclose(fp);
  XLALREAL8VectorListDestroy ( covering );
  covering = NULL;

  /* free memory */
  gsl_matrix_free ( generator0 );
  gsl_matrix_free ( generatorA );
  gsl_matrix_free ( generatorB );
  LALFree ( startPoint.data );

  /* check all (LAL-)memory for leaks... (unfortunately doesn't cover GSL!) */
  LALCheckMemoryLeaks(); 
  
} /* testCovering() */

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

#if 0
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
#endif

/* write a REAL8VectorList to a file */
int
writeREAL8VectorList (FILE *fp, const REAL8VectorList *list)
{
  UINT4 i;
  const REAL8VectorList *node;

  if ( !list || !fp )
    return -1;

  for (node = list; node; node=node->next )
    {
      for (i=0; i < node->entry.length; i++ )
	fprintf (fp, "%g ", node->entry.data[i] );
      
      fprintf (fp, "\n");

    } /* for node=node->next */

  return 0;
} /* writeREAL8VectorList() */

#define SPOKES 60  /* spokes for ellipse-plots */
/** write 2D covering-lattice into file, add metric ellipses */
int 
plot2DCovering (FILE *fp, const REAL8VectorList *list, const REAL8Vector *metric, REAL8 mismatch)
{
  REAL8 Alpha, Delta;
  const REAL8VectorList *node;
  MetricEllipse ellipse;
  LALStatus status = blank_status;
  UINT4 i;

  /* ----- first output all grid-points */
  writeREAL8VectorList( fp, list );
  
  /* ----- if metric is given: output all metric ellipses */
  if ( metric ) 
    {
      fprintf (fp, "\n\n");
      getMetricEllipse(&status, &ellipse, mismatch, metric, 0);
      for( node = list; node; node = node->next ) 
	{
	  Alpha =  node->entry.data[0];
	  Delta =  node->entry.data[1];
	  
	  /* Loop around patch ellipse. */
	  for (i=0; i<=SPOKES; i++) {
	    float c, r, b, x, y;
	
	    c = LAL_TWOPI*i/SPOKES;
	    x = ellipse.smajor * cos(c);
	    y = ellipse.sminor * sin(c);
	    r = sqrt( x*x + y*y );
	    b = atan2 ( y, x );
	    fprintf( fp, "%e %e\n", 
		     Alpha + r*cos(ellipse.angle+b), Delta + r*sin(ellipse.angle+b) );
	  }
	  fprintf ( fp, "\n\n");
      
	} /* for node=node->next */

    } /* if metric */

  return 0;

} /* plot2DCovering() */

/**
 * Monte-Carlo run to measure the covering-radius (and distribution) of
 * a given covering.
 * Algorithm: pick a point at random within [-1,1]^n, get the minimal distance
 * to one of the lattice-points, and write this distance into a file, etc.
 * The covering-radius is given as the cutoff-distance of this distribution.
 *
 * Return -1 on error, 0 if ok.
 */
int
MC_MeasureCovering(FILE *fp, 
		   const REAL8VectorList *covering, 
		   const REAL8Vector *metric, 
		   UINT4 numPoints)
{
  UINT4 run;
  gsl_matrix *gij;
  UINT4 dim;
  REAL8Vector *samplePoint = NULL;
  int seed = 1;

  if ( !covering || !metric || !fp )
    return -1;

  if ( ( gij = XLALmetric2gsl (metric)) == NULL )
    return -1;

  dim = gij->size1;
  
  if ( covering->entry.length != dim )
    return -1;

  if ( (samplePoint = XLALCreateREAL8Vector ( dim )) == NULL )
    return -1;

  srand(seed);

  for (run = 0; run < numPoints; run ++)
    {
      REAL8 dist;
      ChooseRandomPoint ( samplePoint );
      
      dist = FindMinDistance (samplePoint, covering, gij );

      fprintf (fp, "%f\n", dist );
      
    } /* for run < numRuns */

  gsl_matrix_free(gij);
  XLALDestroyREAL8Vector (samplePoint);
  
  return 0;

} /* MC_MeasureCovering() */

/* pick a random-point within [-0.8,0.8]^n */
int
ChooseRandomPoint ( REAL8Vector *point )
{
  UINT4 dim, i;

  if (!point)
    return -1;

  dim = point->length;

  for (i = 0; i < dim; i++)
      point->data[i] = (1.6 * myrand() - 0.8);
  
  return 0;

} /* PickRandomPoint() */

REAL8
FindMinDistance (const REAL8Vector *samplePoint, 
		 const REAL8VectorList *covering, 
		 const gsl_matrix *metric )
{
  const REAL8VectorList *node;
  UINT4 dim = metric->size1;
  gsl_vector_view p1, p2;
  gsl_vector *diff;
  REAL8 mindist2 = -1;

  p1  = gsl_vector_view_array ( samplePoint->data, dim );
  diff = gsl_vector_calloc (dim);


  for ( node = covering; node; node = node->next )
    {
      REAL8 dist2;
      p2  = gsl_vector_view_array ( node->entry.data, dim );

      gsl_vector_memcpy (diff, &(p2.vector) );

      gsl_vector_sub (diff, &(p1.vector) );

      dist2 = XLALMetricScalarProduct (diff, diff, metric);

      if ( (mindist2 < 0) || ( dist2 < mindist2 ) )
	mindist2 = dist2;

    } /* for node=node->next */

  gsl_vector_free(diff);

  return ( sqrt(mindist2) );

}/* FindMinDistance() */
