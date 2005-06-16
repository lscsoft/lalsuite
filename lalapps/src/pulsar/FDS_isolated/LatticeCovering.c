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
 * \file LatticeCovering.c
 * \author: Reinhard Prix
 * $Id$
 *
 * Module for covering metric parameter-spaces with a lattice.
 * 
 * These routines should eventually become useable for implementing 
 * a covering that's practical and as "optimal" as possible.
 * Right now this is just a playground...
 *
 */

/*---------- INCLUDES ----------*/
#include <math.h>

#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <lalapps.h>

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>

RCSID ("$Id$");
NRCSID( LATTICECOVERINGC, "$Id$" );

/*---------- DEFINES ----------*/
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

#define TRUE (1==1)
#define FALSE (1==0)

/*----- SWITCHES -----*/
#define HAVE_INLINE	1
/* uncomment the following to turn off range-checking in GSL vector-functions */
/* #define GSL_RANGE_CHECK_OFF 1*/

/*----- Error-codes -----<lalLaTeX>
\subsection*{Error codes}</lalLaTeX><lalErrTable> */

#define LATTICECOVERING_ENULL 		1
#define LATTICECOVERING_ENONULL		2
#define LATTICECOVERING_EMEM		3
#define LATTICECOVERING_EINPUT		4
#define LATTICECOVERING_ELIST		5
#define LATTICECOVERING_EFUNC		6

#define LATTICECOVERING_MSGENULL 	"Arguments contained an unexpected null pointer"
#define LATTICECOVERING_MSGENONULL	"Output pointer is not NULL"
#define LATTICECOVERING_MSGEMEM		"Out of memory"
#define LATTICECOVERING_MSGEINPUT	"Invald input parameter"
#define LATTICECOVERING_MSGELIST	"Error occurred in list-handling ..."
#define LATTICECOVERING_MSGEFUNC	"Sub-routine failed"

/*</lalErrTable> */

/*---------- local types ----------*/

/** enum-type for denoting several types of lattice */
typedef enum
{
  LATTICE_TYPE_CUBIC = 0,	/**< standard cubic grid: Zn */
  LATTICE_TYPE_ANSTAR,		/**< An*: optimal covering grid */
  LATTICE_TYPE_LAST
} LatticeType;

/** doubly linked list of INT4-vectors (lattice-vectors) */
typedef struct tagINT4VectorList
{
  INT4Vector entry;
  struct tagINT4VectorList *next;
  struct tagINT4VectorList *prev;
} INT4VectorList;

/** singly linked list of REAL8-vectors (physical vectors) */
typedef struct tagREAL8VectorList
{
  REAL8Vector entry;
  struct tagREAL8VectorList *next;
  struct tagREAL8VectorList *prev;
} REAL8VectorList;


INT4VectorList empty_INT4VectorList;
REAL8VectorList empty_REAL8VectorList;

/*---------- Global variables ----------*/



/*---------- local prototypes ----------*/
void
LALLatticeFill (LALStatus *lstat,
		REAL8VectorList **fillGrid,
		const gsl_matrix  *generator,
		const REAL8Vector *startPoint,
		BOOLEAN (*isInside)(const REAL8Vector *point)
		);

int
XLALlatticePoint2physicalPoint ( REAL8Vector *physicalPoint, 
				 const INT4Vector *latticePoint, 
				 const gsl_matrix *generator, 
				 const REAL8Vector *startPoint );

REAL8 XLALMetricScalarProduct (const gsl_vector *vector1,
			       const gsl_vector *vector2,	
			       const gsl_matrix *metric);


int
XLALMetricGramSchmidt(gsl_matrix **orth, 		
		      const gsl_matrix *colvect,	
		      const gsl_matrix *metric);
int
XLALCanonicalGenerator ( gsl_matrix **outmatrix,
			 const gsl_matrix *inmatrix,
			 const gsl_matrix *gij );	

int
XLALSquareGeneratingMatrix(gsl_matrix **outmatrix,
			   const gsl_matrix *inmatrix);

int
XLALGetGeneratingMatrix (gsl_matrix **outmatrix,
			 UINT4 dimension,
			 LatticeType type);

/* list-handling functions */
INT4VectorList *INT4VectorListAddEntry (INT4VectorList *head, const INT4Vector *entry);
int INT4VectorListRelinkElement (INT4VectorList *head, INT4VectorList *element);
void INT4VectorListRemoveElement (INT4VectorList *element);
REAL8VectorList *REAL8VectorListAddEntry (REAL8VectorList *head, const REAL8Vector *entry);
BOOLEAN isINT4PointInList ( INT4Vector *point, INT4VectorList *list );
void INT4VectorListDestroy (INT4VectorList *head);
void REAL8VectorListDestroy (REAL8VectorList *head);

/* misc helper functions */
BOOLEAN isSymmetric (const gsl_matrix *Sij);
int writeREAL8VectorList (FILE *fp, const REAL8VectorList *list);
int print_matrix (const gsl_matrix *m);
int print_vector (const gsl_vector *v);

/* test-functions */
int main(void);
void testGS(void);
void testCovering(void);

/*==================== FUNCTION DEFINITIONS ====================*/


/** Fill the given parameter-space by a lattice of specified 
 * generating matrix.
 * 
 * NOTE: the input generating-matrix (generator) needs to be scaled
 * correctly to the required covering radius, also, it needs to be in 
 * canonical _square matrix_ form. 
 *
 * NOTE2: as always in this module, the generating matrix contains
 * the lattice-vectors as _rows_
 */
void
LALLatticeFill (LALStatus *lstat,
		REAL8VectorList **fillGrid,	/**< OUT: final fill-grid (physical points)*/
		const gsl_matrix  *generator,	/**< IN: _SQUARE_ generating matrix for lattice*/
		const REAL8Vector *startPoint, 	/**< IN: physical startpoint for filling */
		BOOLEAN (*isInside)(const REAL8Vector *point) /**< IN: boundary-condition */
		)
{
  UINT4 dim;	/**< dimension of parameter-space to fill */
  UINT4 i;
  INT4VectorList openEnds = empty_INT4VectorList;	/**< list of "open ends" (lattice-points) */
  INT4VectorList gridPoints = empty_INT4VectorList;	/**< resulting grid (lattice-points) */
  REAL8VectorList realPoints = empty_REAL8VectorList;	/**< physical coordinates of grid-points */

  INT4Vector  *latticePoint = NULL;		/* lattice-coordinates (Z^N) */
  REAL8Vector *physicalPoint = NULL;		/* physical coordinates (R^N) */

  INITSTATUS( lstat, "LALLatticeFill", LATTICECOVERINGC );
  ATTATCHSTATUSPTR (lstat); 

  /* This traps coding errors in the calling routine. */
  ASSERT ( fillGrid != NULL, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );  
  ASSERT ( *fillGrid == NULL,lstat, LATTICECOVERING_ENONULL, LATTICECOVERING_MSGENONULL );
  ASSERT ( generator, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( startPoint, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( startPoint->data, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );

  if ( generator->size1 != generator->size2 )	/* need square generator */
    {
      LALPrintError ("\nERROR: LatticeFill() requires a  SQUARE generating matrix!\n\n");
      ABORT (lstat, LATTICECOVERING_EINPUT, LATTICECOVERING_MSGEINPUT);      
    }

  if ( ! (*isInside)(startPoint) )	/* startPoint has to be inside */
    {
      LALPrintError ("\nERROR: startPoint must lie withing the covering-region!\n\n");
      ABORT (lstat, LATTICECOVERING_EINPUT, LATTICECOVERING_MSGEINPUT);
    }

  dim = startPoint->length;	/* dimension of parameter-space to cover */

  if ( (generator->size1 != dim) )
    {
      LALPrintError ("\nERROR: all input-dimensions must agree (dim=%d)\n\n", dim);
      ABORT (lstat, LATTICECOVERING_EINPUT, LATTICECOVERING_MSGEINPUT);
    }

  /* ---------- prepare memory for one grid-coordinate and one physical grid-point */
  TRY ( LALI4CreateVector (lstat->statusPtr, &latticePoint, dim), lstat);
  LALDCreateVector (lstat->statusPtr, &physicalPoint, dim);
  BEGINFAIL (lstat) {
    TRY ( LALI4DestroyVector(lstat->statusPtr, &latticePoint), lstat);
  } ENDFAIL(lstat);
  /* initialize vectors to 0 */
  memset ( latticePoint->data, 0, dim * sizeof (latticePoint->data[0]) );
  memset ( physicalPoint->data, 0, dim * sizeof (physicalPoint->data[0]));
  

  /* ----- start by adding startPoint coordinates (0,0,0,,) to the list of 'open-ends' */
  /* we don't need to set these coordinates, grid-point is already set to (0,0,0,,,)  */
  if ( NULL == INT4VectorListAddEntry (&openEnds, latticePoint)) 
    {	
      /* NOTE: head always stays empty for simplicity! */
      LALPrintError ("\nERROR: INT4VectorListAddEntry () failed!\n\n");
      ABORT (lstat, LATTICECOVERING_ELIST, LATTICECOVERING_MSGELIST);
    }

   
  /* ----- (*) take coordinates of next open-end from hash-list of open-ends */
  while ( openEnds.next )
    {
      INT4Vector *thisPoint = NULL;
      INT4VectorList *thisElement = openEnds.next; 

      /* get lattice-coordinates of this point (pointer to entry) */
      thisPoint = &(thisElement->entry);
      /* find its physical coordinates */
      XLALlatticePoint2physicalPoint ( physicalPoint, thisPoint, generator, startPoint );
      if ( xlalErrno )
	{
	  int code = xlalErrno;
	  XLALClearErrno(); 
	  LALPrintError ("\nERROR: latticePoint2physicalPoint() failed (xlalErrno = %d)!\n\n", code);
	  ABORT (lstat, LATTICECOVERING_EFUNC, LATTICECOVERING_MSGEFUNC);
	}

      /* is it inside the filling-region?: */
      if ( ! (*isInside)(physicalPoint) )
	{ /* NO */
	  /* remove this lattice-point from list of open ends and continue */
	  INT4VectorListRemoveElement ( thisElement );
	  thisPoint = NULL;
	  continue;
	}
      else
	{ /* YES */
	  /* move this point into the list of grid-points (in lattice-coordinates) */
	  INT4VectorListRelinkElement ( &gridPoints, thisElement );

	  /* and store its physical coordinates in a REAL8VectorList as well (avoids
	   * re-calculating physical coordinates again later...) */
	  if ( NULL == REAL8VectorListAddEntry ( &realPoints, physicalPoint) )
	    {	
	      /* NOTE: head always stays empty for simplicity! */
	      LALPrintError ("\nERROR: REAL8VectorListAddEntry () failed!\n\n");
	      ABORT (lstat, LATTICECOVERING_ELIST, LATTICECOVERING_MSGELIST);
	    }

	  /* generate coordinates of this point's 2*dim neighbors, */
	  for (i=0; i < 2 * dim ; i++)
	    {
	      memcpy ( latticePoint->data, thisPoint->data, dim * sizeof(thisPoint->data[0]));
	      if ( i % 2 )
		latticePoint->data[ i / 2 ] ++;
	      else
		latticePoint->data[ i / 2 ] --;

	      /* discard all those already in list of open-ends
	       * and those already in lattice-points (via hash-algorithm..)
	       */
	      if ( isINT4PointInList ( latticePoint, &gridPoints ) )
		continue;
	      if ( isINT4PointInList ( latticePoint, &openEnds ) )
		continue;
	      /* add the other ones to list of open-ends */
	      if ( NULL == INT4VectorListAddEntry ( &openEnds, latticePoint ) )
		{
		  /* NOTE: head always stays empty for simplicity! */
		  LALPrintError ("\nERROR: REAL8VectorListAddEntry () failed!\n\n");
		  ABORT (lstat, LATTICECOVERING_ELIST, LATTICECOVERING_MSGELIST);
		}
	    
	    } /* for i < 2 * dim */

	} /* if point inside filling-region */


      /* start from (1) until no more open ends left */
    } /* while (openEnds->next) */
  
  /* return linked list of physical points */
  (*fillGrid) = realPoints.next;

  /* clean up memory */
  INT4VectorListDestroy ( gridPoints.next );	/* free list of lattice-coordinates */
  TRY ( LALI4DestroyVector(lstat->statusPtr, &latticePoint), lstat);
  TRY ( LALDDestroyVector(lstat->statusPtr, &physicalPoint), lstat);

  DETATCHSTATUSPTR (lstat);

  RETURN( lstat );

} /* LALLatticeFill() */

/** calculate the physical coordinates of a lattice-vector for given
 * generating-matrix and start-point (=0,0,0...0) of the lattice 
 *
 * the algorithm is simply : phys^i = sum_{l=0} lat^l basis_l^i,
 * where the l-th basis-vector is the l-th row of the generating matrix.
 * 
 * NOTE: the memory for physicalPoint needs to be allocated already, and the 
 * dimensions of all vectors and matrices passed to this functions must agree!
 */
int
XLALlatticePoint2physicalPoint ( REAL8Vector *physicalPoint, 	/**< OUT: physical coordinates */
				 const INT4Vector *latticePoint,/**< IN: lattice-coordinates */
				 const gsl_matrix *generator, 	/**< IN: generating vectors as rows */
				 const REAL8Vector *startPoint )/**< IN: phys.coordinates of (0,0,...0)*/
{
  UINT4 dim;
  gsl_matrix *buffer; 
  gsl_vector *res;
  gsl_vector_view tmp;
  UINT4 l;

  /* check validity of input */
  if ( !physicalPoint || !latticePoint || !generator || !startPoint ||
       !physicalPoint->data || !latticePoint->data || !startPoint->data || !generator->data )
    {
      LALPrintError ("\nNULL Input received!\n\n");
      XLAL_ERROR ( "XLALlatticePoint2PhysicalPoint", XLAL_EINVAL);
    }

  dim = physicalPoint->length;
  if ( (latticePoint->length != dim) || (generator->size1 != dim) || (generator->size2 != dim) ) 
    {
      LALPrintError ("\nInconsistent dimensions in input-vectors/matrices!\n\n");
      XLAL_ERROR ( "XLALlatticePoint2PhysicalPoint", XLAL_EINVAL);
    }

  if ( (buffer = gsl_matrix_calloc (dim, dim)) == NULL ) {
    XLAL_ERROR ( "XLALlatticePoint2PhysicalPoint", XLAL_ENOMEM);
  }

  /* create a local copy of the generating-matrix */
  gsl_matrix_memcpy ( buffer, generator );

  /* create a gsl-vector for summing up the lat^l basis_l vectors (->final result) */
  if ( (res = gsl_vector_calloc (dim)) == NULL ) {
    gsl_matrix_free (buffer);
    XLAL_ERROR ( "XLALlatticePoint2PhysicalPoint", XLAL_ENOMEM);
  }

  /* get a vector-view on the startPoint and copy it into res */
  tmp = gsl_vector_view_array (startPoint->data, dim);
  gsl_vector_memcpy ( res, &(tmp.vector) );

  /* now multiply each row-vector base_l with the integer lat^l
   * and add this to the start-vector (already in res) */
  for (l=0; l < dim; l++)
    {
      tmp = gsl_matrix_row ( buffer, l );	/* vector of row(l) of buffer */
      gsl_vector_scale ( &(tmp.vector), latticePoint->data[l]);	/* lat^l * base_l */
      gsl_vector_add ( res, &(tmp.vector) );			/* sum it up */
    }
  
  /* convert final answer back into REAL8Vector */
  /* !!!! NOTE: this assumes that sizeof(REAL8)==sizeof(double) !!!
   * if that's not true, this will fall on its head!!
   * (but that should hopefully not be a worry and 
   * individual element-copying just seems too inefficient
   */
  memcpy (physicalPoint->data, res->data, dim * sizeof (res->data[0]));

  /* free memory */
  gsl_matrix_free (buffer);
  gsl_vector_free (res);
  

  return 0;

} /* XLALlatticePoint2physicalPoint() */

/** Scalar product of two vectors with respect to the given metric.
 *  <v1, v2> = g_ij v1^i v^2j 
 */
REAL8 
XLALMetricScalarProduct (const gsl_vector *v1,
			 const gsl_vector *v2,	
			 const gsl_matrix *gij)
{
  size_t dim;
  UINT4 i, j;
  REAL8 prod;	/* final scalar product */

  /* check that input is non-zero */
  if ( (!v1) || (!v2) || (!gij) ) {
    LALPrintError ("NULL Input received.");
    XLAL_ERROR_REAL8("XLALMetricScalarProduct", XLAL_EINVAL);
  }

  /* check that gij is symmetric */
  if ( !isSymmetric(gij) ) {
    LALPrintError ("Input 'metric' has to be symmetric!");
    XLAL_ERROR_REAL8("XLALMetricScalarProduct", XLAL_EINVAL);
  }

  dim = gij->size1;

  /* check that vectors have correct sizes */
  if ( (v1->size != dim) || (v2->size != dim) ) {
    LALPrintError ("Vectors v1, v2 must have same dimension as metric gij");
    XLAL_ERROR_REAL8("XLALMetricScalarProduct", XLAL_EINVAL);
  }
      
  /* calculate scalar product */
  prod = 0;
  for ( i=0; i < dim; i++)
    for (j=0; j < dim; j++)
      prod += gsl_vector_get (v1, i) * gsl_matrix_get (gij, i, j) * gsl_vector_get (v2, j);

  return prod;

} /* XLALMetricScalarProduct() */


/** Gram-Schmidt orthogonalization of lin. indep. vectors using a given metric.
 * NOTE: this is a straightforward, probably naive implementation of the basic
 * algorithm, completely ignorant about numerically more robust or faster algorithms
 * to do this... [FIXME?].
 *
 * The vectors in input and output are the matrix-rows!
 *
 * NOTE2: the memory for outvects is allocated in here by gsl_matrix_alloc()
 * and needs to be free'ed by the caller via gsl_matrix_free() !
 */
int
XLALMetricGramSchmidt(gsl_matrix **outvects,	/**< OUT: orthonormal row vects */
		      const gsl_matrix *invects,/**< matrix of row-vectors */
		      const gsl_matrix *gij)	/**< metric */
{
  UINT4 numvects, vectdim;
  UINT4 i, j;
  gsl_matrix *orth = NULL;
  gsl_vector_view *ui = NULL;	/* array of orthonormal result-vectors */

  /* check NULL-vectors on input */
  if ( (!invects) || (!gij) || (!outvects) ) {
    LALPrintError ("NULL Input received.");
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EINVAL);
  }
  
  /* check that output 'outvects' points to a NULL-vector! */
  if ( *outvects != NULL ) {
    LALPrintError ("Output-vector not set to NULL");
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EINVAL);
  }

  /* check that gij is symmetric */
  if ( !isSymmetric(gij) ) {
    LALPrintError ("Input 'metric' has to be symmetric!");
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EINVAL);
  }

  numvects = invects->size1;	/* number of rows! */
  vectdim = invects->size2;	/* number of columns */

  /* can't have more vectors than dimensions */
  if ( numvects > vectdim ) {
    LALPrintError ("Input vectors are not linearly independent");
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EINVAL);
  }
  
  /* vector-dimension has to be consistent with metric */
  if ( vectdim != gij->size1 ) {
    LALPrintError ("Dimension of input vectors inconsistent with metric");
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EINVAL);
  }


  /* prepare output-matrix for orthonormalized vectors in rows*/
  orth = gsl_matrix_calloc ( numvects, vectdim);
  if ( orth == NULL ) {
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_ENOMEM);
  }

  /* prepare vector view-arraw on final orthonormal vectors */
  ui = LALCalloc (1, numvects * sizeof(gsl_vector_view) );
  for (i=0; i < numvects; i++)
    ui[i] = gsl_matrix_row ( orth, i );


  /* placeholder for temporary vector */
  gsl_vector *para = gsl_vector_alloc(vectdim);
  if ( para == NULL ) {
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_ENOMEM);
  }

  /*---------- main algorithm ---------- */
  for (i=0; i < numvects; i++)
    {
      REAL8 norm;
      /*---------- Step 1) orthogonalize wrt to previous vectors ----------*/
      
      /* view on input-vector i (convenience) */
      gsl_vector_const_view vi = gsl_matrix_const_row (invects, i);

      printf ("Read row i=%d of in-matrix:\n", i);
      print_vector ( &(vi.vector) );

      gsl_vector_memcpy ( &(ui[i].vector), &(vi.vector) );
      
      for (j=0; (i>0) && (j < i); j++)
	{
	  REAL8 proj;

	  gsl_vector_memcpy ( para, &(ui[j].vector) );

	  proj = XLALMetricScalarProduct(&vi.vector, &(ui[j].vector), gij);
	  if( gsl_vector_scale ( para, proj) ) {	/* para_j = <vi,uj> uj */
	    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EFUNC);
	  }

	  if ( gsl_vector_sub ( &(ui[i].vector), para) ) {	/* ui -= para_j */
	    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EFUNC);
	  }

	} /* for j < i-1 */

      /*---------- Step 2) normalize ---------- */
      
      norm = XLALMetricScalarProduct ( &(ui[i].vector), &(ui[i].vector), gij );
      norm = sqrt(norm);

      if ( gsl_vector_scale ( &(ui[i].vector), 1.0/norm ) ) {
	XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EFUNC);
      }
      
    } /* for i < numvects */

  /* free memory */
  gsl_vector_free (para);
  LALFree (ui);

  /* return result */
  (*outvects) = orth;

  return 0;

} /* XLALMetricGramSchmidt() */

/** Turn a general (non-quadratic) Generating Matrix into a quadratic one.
 * The input matrix must have columns >= rows, the rows reprenting the 
 * lattice vectors. This algorithm simply proceeds by constructing an
 * (Euclidean!) orthonormal basis out of the lattice vectors (using GramSchmidt),
 * and then expressing the lattice-vectors in this new basis.
 * 
 * NOTE: the memory for 'outmatrix' is allocated in here via gsl_matrix_alloc()
 *       and has to be free'ed by the caller via gsl_matrix_free() !
 */
int
XLALSquareGeneratingMatrix(gsl_matrix **outmatrix, 	/**< OUT: square generating matrix */
			   const gsl_matrix *inmatrix)	/**< generating matrix (generally cols >= rows) */
{
  UINT4 rows, cols;
  gsl_matrix *sq = NULL;	/* output: square generating matrix */

  /* check NULL-vectors on input */
  if ( inmatrix == NULL ) 
    {
      LALPrintError ("NULL Input received.");
      XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_EINVAL);
    }
  
  /* check that output 'outmatrix' points to a NULL-vector! */
  if ( *outmatrix != NULL ) 
    {
      LALPrintError ("Output-vector not set to NULL");
      XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_EINVAL);
    }

  rows = inmatrix->size1;
  cols = inmatrix->size2;

  /* rows need to be lattice vectors, and linearly independent */
  if ( rows > cols ) 
    {
      LALPrintError ("ERROR: input-matrix must have full row-rank!\n");
      XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_EINVAL);
    }

  /* allocate output matrix */


  /* if input-matrix is quadratic, we're done */
  if ( rows == cols )
    {
      if ( (sq = gsl_matrix_calloc (rows, rows)) == NULL ) {
	XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_ENOMEM);
      }
    gsl_matrix_memcpy ( sq, inmatrix );
    }
  else
    {
      gsl_matrix *gij = gsl_matrix_alloc( cols, cols );
      if ( !gij ){
	XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_ENOMEM);
      }
      gsl_matrix_set_identity (gij);	/* use Euklidean metric for orthonormalization*/

      /* now express generator in 'canonical coordinates' wrt. Euklidean metric */
      XLALCanonicalGenerator ( &sq, inmatrix, gij );
      if ( xlalErrno )
	{
	  int code = xlalErrno;
	  XLALClearErrno(); 
	  LALPrintError ("\nERROR: XLALCanonicalGenerator() failed (xlalErrno = %d)!\n\n", code);
	  XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_EFUNC);
	}

      gsl_matrix_free (gij);

    } /* if cols > rows  */

  /* we're done: return result */
  (*outmatrix) = sq;

  return 0;
} /* XLALSquareGeneratingMatrix() */

/** express the generating matrix in "canonical coordinates":
 * what we mean be this is: express its coordinates in a basis that is
 *    (a) of identical span as the space spanned by the generating matrix 'inmatrix'
 *    (b) orthonormal with respect to the given metric gij
 *
 * NOTE: the metric gij must obviously be symmetric, and its dimension
 *       MUST be identical to the space in which the generator 'inmatrix'
 *       is expressed, i.e.   dim(gij) == columns(inmatrix)
 *
 * NOTE2: the resulting canonical matrix will _always_ be square
 *        with dim(outmatrix) = rank(inmatrix) = rows(inmatrix)
 *        (this is obvious from the construction, see condition (a) )
 *
 */
int
XLALCanonicalGenerator ( gsl_matrix **outmatrix,	/**< OUT: the resulting (SQUARE) canonical generator */
			 const gsl_matrix *inmatrix,	/**< IN: general generating matrix */
			 const gsl_matrix *gij )	/**< IN: metric */
{    
  UINT4 rows, cols;
  gsl_matrix *basis = NULL;	/* orthonormal basis */  
  gsl_matrix *ret = NULL;	/* resulting canonical square matrix */

  /* ----- check input consistency */
  if ( !inmatrix || !gij || !outmatrix ) {
    LALPrintError ("\nERROR: NULL input\n");
    XLAL_ERROR("XLALCanonicalGenerator", XLAL_EINVAL);
  }

  if ( *outmatrix != NULL ) {
    LALPrintError ("\nERROR: *outmatrix must be NULL!\n\n");
    XLAL_ERROR("XLALCanonicalGenerator", XLAL_EINVAL);
  }

  if (!isSymmetric (gij) ) {
    LALPrintError ("\nERROR: metric gij must be symmetric!\n\n");
    XLAL_ERROR("XLALCanonicalGenerator", XLAL_EINVAL);
  }

  rows = inmatrix->size1;
  cols = inmatrix->size2;

  if ( gij->size1 != cols ) {
    LALPrintError ("\nERROR: dimension of metric gij must equal cols(inmatrix) !\n\n");
    XLAL_ERROR("XLALCanonicalGenerator", XLAL_EINVAL);
  }

  /* ----- find orthonormal basis */
  if ( XLALMetricGramSchmidt (&basis, inmatrix, gij) ) {
    XLAL_ERROR("XLALCanonicalGenerator", XLAL_EFUNC);
  }

  /* prepare memory for resulting-matrix (rows x rows) : */
  if ( (ret = gsl_matrix_calloc (rows, rows)) == NULL ) {
    XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_ENOMEM);
  }
  
  /* ----- express generating matrix in this new Euklidean basis: inmatrix.basis^T */
  
  /* from the gsl-documentation:
   * int gsl_blas_dgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, 
   *                    double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
   * These functions compute the matrix-matrix product and 
   * sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for 
   * TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
   */
  if ( gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, inmatrix, basis, 0.0, ret)) {
    LALPrintError ("ERROR: Call to  gsl_blas_dgemm() failed\n");
    XLAL_ERROR("XLALCanonicalGenerator", XLAL_EFUNC);
  }
  
  /* free memory */
  gsl_matrix_free (basis);
  basis = NULL;
  
  /* return resulting canonical generator */
  (*outmatrix) = ret;

  return 0;

} /* XLALCanonicalGenerator() */





/** Return the (canonical, i.e. not necessarily quadratic) n-dimensional
 * generating matrix for various lattices. 
 * 
 * NOTE: because these lattices are intended for covering, we scale
 *       them so that their covering-radius becomes unity.
 *       This allows the user to later-on scale these easily to any
 *       desired covering-radius without having to know the lattice...
 *       (Remembering that if you scale the generator of a lattice by M' = c M, then
 *       the covering radius scales as R' = c R)
 *
 * NOTE2: the memory for 'outmatrix' is allocated in here via gsl_matrix_alloc()
 *       and has to be free'ed by the caller via gsl_matrix_free() !
 *
 * REFERENCE: for the definition and properties of these lattices, see
 *  	J.H. Conway and N.J.A. Sloane, "Sphere Packings, Lattices and Groups", 
 *  	vol. 290 of Grundlehren der mathematischen Wissenschaften, 
 * 	(Springer, New York, U.S.A., 1999), 3rd edition.
 *
 */
int
XLALGetGeneratingMatrix (gsl_matrix **outmatrix,	/**< OUT: generating matrix */
			 UINT4 dimension,		/**< IN: number of dimensions */
			 LatticeType type)		/**< IN: type of lattice */
{
  gsl_matrix *generator = NULL;	/* output: generating matrix */
  REAL8 coveringRadius;

  /* check that output 'outmatrix' points to a NULL-vector! */
  if ( *outmatrix != NULL ) 
    {
      LALPrintError ("Output-vector not set to NULL");
      XLAL_ERROR("XLALGetGeneratingMatrix", XLAL_EINVAL);
    }

  switch (type)
    {
    case LATTICE_TYPE_CUBIC:
      generator = gsl_matrix_calloc( dimension, dimension );
      gsl_matrix_set_identity (generator);	/* trivial generating matrix */

      /* covering radius is sqrt(n)/2 */
      coveringRadius = sqrt (1.0 * dimension) / 2.0;

      break;

    case LATTICE_TYPE_ANSTAR:
      /* the generating matrix for An* can be written as:
       * 
       * | 1      -1        0       0  0  ....  0      |
       * | 1       0       -1       0  0  ....  0      |
       * | 1       0        0      -1  0  ....  0      | 
       * | ...    ...      ....    ... ....    ...     |
       * |-n/(n+1) 1/(n+1) 1/(n+1) ...        1/(n+1)  |
       * 
       */
      generator = gsl_matrix_calloc( dimension, dimension+1 );
      UINT4 row, col;

      for (row=0; row < dimension; row ++)
	{
	  for (col=0; col < dimension + 1; col++)
	    {
	      /* ---------- find value for that matrix element ----------*/
	      REAL8 val;
	      if ( row < dimension-1 )
		{
		  if ( col == 0 )
		    val = 1.0;
		  else if (col == row + 1)
		    val = -1.0;
		  else
		    continue;
		}
	      else
		{
		  if ( col == 0 )
		    val = - 1.0 * dimension / ( dimension + 1.0);
		  else
		    val =   1.0 / (dimension + 1.0);
		}

	      /* ---------- set matrix element ---------- */
	      gsl_matrix_set (generator, row, col, val);

	    } /* for col < dim + 1 */

	} /* for row < dim */

      /* covering Radius of An* is R = sqrt( n*(n+2) / (12*(n+1)) ), see C&S */
      coveringRadius = sqrt ( 1.0 * dimension * (dimension + 2.0) / (12.0 * (dimension + 1) ));

      break;

    default:
      LALPrintError ("Illegal value for lattice-type (%d)\n", type);
      XLAL_ERROR("XLALGetGeneratingMatrix", XLAL_EINVAL);
      break;
      
    } /* switch(type) */

  /* Now scale the generating matrix by 1/R, such that its covering radius becomes R=1 */
  gsl_matrix_scale ( generator, 1.0 / coveringRadius );

  /* return generating matrix */
  (*outmatrix) = generator;

  return 0;
} /* XLALGetGeneratingMatrix() */


/*----------------------------------------------------------------------
 * list handling tools
 *----------------------------------------------------------------------*/

/** add a new element at the end of the list 'head', _copy_ the given entry there,
 * and return pointer to the new list-entry. 
 * 
 * NOTE: this function is rather permissive in that it takes the vector-dimension 
 * from the new entry 'el', without requiring this to be equal to the dimension of the 
 * other list-entries...
 */
INT4VectorList *
INT4VectorListAddEntry (INT4VectorList *head, const INT4Vector *entry)
{
  UINT4 dim;
  INT4VectorList *ptr = NULL;	/* running list-pointer */
  INT4VectorList *newElement = NULL;	/* new list-element */
  /* check illegal input */
  if ( (head == NULL) || (entry == NULL) )
    return NULL;

  /* find tail of list */
  ptr = head;
  while ( ptr->next )
    ptr = ptr->next;

  /* construct new list-element */
  dim = entry->length;
  if ( (newElement = LALCalloc (1, sizeof (*newElement))) == NULL)
    return NULL;
  if ( (newElement->entry.data = LALCalloc (dim, sizeof(entry->data[0]))) == NULL ) {
    LALFree (newElement);
    return NULL;
  }
  newElement->entry.length = dim;
  memcpy (newElement->entry.data, entry->data, dim * sizeof(entry->data[0]) );

  /* link this to the tail of list */
  ptr->next = newElement;
  newElement->prev = ptr;
  
  return newElement;

} /* INT4VectorListAddEntry() */


/** add a new element at the end of the list 'head', _copy_ the given entry there,
 * and return pointer to the new list-entry. 
 * 
 * NOTE: this function is rather permissive in that it takes the vector-dimension 
 * from the new entry 'el', without requiring this to be equal to the dimension of the 
 * other list-entries...
 */
REAL8VectorList *
REAL8VectorListAddEntry (REAL8VectorList *head, const REAL8Vector *entry)
{
  UINT4 dim;
  REAL8VectorList *ptr = NULL;	/* running list-pointer */
  REAL8VectorList *newElement = NULL;	/* new list-element */
  /* check illegal input */
  if ( (head == NULL) || (entry == NULL) )
    return NULL;

  /* find tail of list */
  ptr = head;
  while ( ptr->next )
    ptr = ptr->next;

  /* construct new list-element */
  dim = entry->length;
  if ( (newElement = LALCalloc (1, sizeof (*newElement))) == NULL)
    return NULL;
  if ( (newElement->entry.data = LALCalloc (dim, sizeof(entry->data[0]))) == NULL ) {
    LALFree (newElement);
    return NULL;
  }
  newElement->entry.length = dim;
  memcpy (newElement->entry.data, entry->data, dim * sizeof(entry->data[0]) );

  /* link this to the tail of list */
  ptr->next = newElement;
  newElement->prev = ptr;
  
  return newElement;

} /* REAL8VectorListAddEntry() */




/** "relink" the given element (from whatever list) to the end of the given list 'head'
 * return 0 if OK, negative number on error...
 */ 
int
INT4VectorListRelinkElement (INT4VectorList *head, INT4VectorList *element)
{
  INT4VectorList *ptr;

  if ( !head || !element )
    return -1;

  /* unlink element from its list */
  if ( element->next )
    (element->next)->prev = element->prev;
  if ( element->prev )
    (element->prev)->next = element->next;

  /* find tail of list */
  ptr = head;
  while ( ptr->next )
    ptr = ptr->next;

  /* link element at the end of given list 'head' */
  ptr->next = element;
  element->next = NULL;
  element->prev = ptr;

  return 0;

} /* INT4VectorListRelinkElement() */


/** remove (and free) the given list-element from the list */
void
INT4VectorListRemoveElement (INT4VectorList *element)
{
  if ( element == NULL )	/* ignore empty input */
    return;

  /* un-link the given element from list */
  if ( element->prev )
    (element->prev)->next = element->next;
  
  if ( element->next )
    (element->next)->prev = element->prev;
  
  /* free the element */
  LALFree (element->entry.data);
  LALFree (element);
  
  return;

} /* INT4VectorListRemoveElement() */

/** 'list-destructor' for INT4VectorList: free a complete list.
 *
 * NOTE: 'head' will be freed too, so make
 * sure not to pass a non-freeable head (like an automatic variabe) 
 */
void
INT4VectorListDestroy (INT4VectorList *head)
{
  INT4VectorList *ptr, *next;

  if ( !head )
    return;

  next = head;

  do 
    {
      /* step to next element */
      ptr = next;
      /* remember pointer to next element */
      next = ptr->next;
      /* free current element */
      LALFree (ptr->entry.data);
      LALFree (ptr);

    } while ( (ptr = next) != NULL );
  
  return;
} /* INT4VectorListDestroy() */

/** 'list-destructor' for REAL8VectorList: free a complete list.
 *
 * NOTE: 'head' will be freed too, so make
 * sure not to pass a non-freeable head (like an automatic variabe) 
 */
void
REAL8VectorListDestroy (REAL8VectorList *head)
{
  REAL8VectorList *ptr, *next;

  if ( !head )
    return;

  next = head;

  do 
    {
      /* step to next element */
      ptr = next;
      /* remember pointer to next element */
      next = ptr->next;
      /* free current element */
      LALFree (ptr->entry.data);
      LALFree (ptr);

    } while ( (ptr = next) != NULL );
  
  return;
} /* REAL8VectorListDestroy() */




/** search given list (or hash-table?) for the given point 
 */
BOOLEAN
isINT4PointInList ( INT4Vector *point, INT4VectorList *list )
{
  UINT4 dim;
  INT4VectorList *ptr;
  
  if ( !point || !list )
    return FALSE;

  dim = point->length;

  /* for now: simple linear search through the list ... */
  ptr = list;
  do 
    {
      if (ptr->entry.length != dim )
	continue;

      if ( memcmp( ptr->entry.data, point->data, dim * sizeof(point->data[0])) == 0 )
	return TRUE;	/* found it! */

    } while ( (ptr = ptr->next) != NULL ); 

  return FALSE;

} /* isINT4PointInList() */

/*----------------------------------------------------------------------
 * misc helper functions
 *----------------------------------------------------------------------*/

/** check if matrix is symmetric */
BOOLEAN
isSymmetric (const gsl_matrix *Sij)
{
  size_t dim; 
  UINT4 i, j;

  if ( !Sij )
    return -1;

  if ( Sij->size1 != Sij->size2 )	/* not quadratic */
    return FALSE;

  dim = Sij->size1;
  for (i=0; i < dim; i++)
    for (j=i; j < dim; j++)
      if ( gsl_matrix_get (Sij, i, j) != gsl_matrix_get (Sij, j, i) )
	return FALSE;

  return TRUE;

} /* isSymmetric() */

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
}


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
}


/*--------------------------------------------------*/
/* Test function(s) */
/*--------------------------------------------------*/
int main (void)
{

  lalDebugLevel = 1;

  /* use this for production code: dont' let gsl abort on errors, but only return error-codes to caller! */
  /*
  gsl_set_error_handler_off ();
  */

  /*  testGS(); */

  testCovering();

  return 0;
}

BOOLEAN testArea1 ( const REAL8Vector *point );

void
testCovering (void)
{
  LALStatus lstat = blank_status;
  REAL8VectorList *covering = NULL;
  UINT4 dim = 2;
  FILE *fp;
  gsl_matrix *generator = NULL;
  gsl_matrix *generator0 = NULL;
  REAL8Vector startPoint;
  
  XLALGetGeneratingMatrix (&generator0, dim, LATTICE_TYPE_ANSTAR);

  gsl_matrix_scale ( generator0, 0.3 );

  XLALSquareGeneratingMatrix( &generator, generator0 );

  startPoint.length = dim;
  startPoint.data = LALCalloc (dim, sizeof(startPoint.data[0]) ); /* already (0,0) */

  LAL_CALL( LALLatticeFill (&lstat, &covering, generator, &startPoint, testArea1), &lstat);

  if ( (fp = fopen ( "test_lattice.dat", "wb" )) == NULL )
    {
      LALPrintError ("\nFailed to open 'test_lattice.dat' for writing!\n\n");
      return;
    }
  writeREAL8VectorList (fp, covering);

  fclose(fp);

  /* free memory */
  gsl_matrix_free ( generator0 );
  gsl_matrix_free ( generator );
  LALFree ( startPoint.data );

  REAL8VectorListDestroy ( covering );

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

  XLALMetricGramSchmidt ( &orto, &(m1.matrix), &(gij.matrix) );

  printf ("\nMetric:\n");
  print_matrix (&(gij.matrix));

  printf ("\nInput-matrix:\n");
  print_matrix (&(m1.matrix));

  printf ("\nResulting orthogonal matrix: \n");
  print_matrix (orto);


  XLALSquareGeneratingMatrix( &sGM, &(m1.matrix) );

  printf ("Square generating matrix:\n");
  print_matrix (sGM);

  XLALGetGeneratingMatrix (&testM, 4, LATTICE_TYPE_ANSTAR);
  printf ("produced Generating Matrix: \n");
  print_matrix (testM);

} /* testGS() */
