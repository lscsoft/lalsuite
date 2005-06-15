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


#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>

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
} REAL8VectorList;


INT4VectorList empty_INT4VectorList;

/*---------- Global variables ----------*/



/*---------- local prototypes ----------*/
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
XLALSquareGeneratingMatrix(gsl_matrix **outmatrix,
			   const gsl_matrix *inmatrix);

int
XLALGetGeneratingMatrix (gsl_matrix **outmatrix,
			 UINT4 dimension,
			 LatticeType type);

/* list-handling functions */
INT4VectorList *append2INT4VectorList (INT4VectorList *head, const INT4Vector *el);

/* misc helper functions */
BOOLEAN isSymmetric (const gsl_matrix *Sij);

/* test-functions */
int main(void);
void testGS(void);
int print_matrix (const gsl_matrix *m);
int print_vector (const gsl_vector *v);

/*==================== FUNCTION DEFINITIONS ====================*/


/** Cover given parameter-space by lattice of specified type.
 * The central function of this module: cover the parameter-space 
 * 
 * NOTE: the input generating-matrix (generator) needs to be scaled
 * correctly to the required covering radius and is supposed to be
 * expressed in the correct coordinates of the parameter-space.
 * Also, it needs to be in canonical *square matrix* form.
 * All this needs to be prepared before this function is called....
 *
 * NOTE2: as always in this module, the generating matrix contains
 * the lattice-vectors as _rows_
 */
void
LALLatticeCovering (LALStatus *lstat,
		    REAL8VectorList **coveringGrid, 	/**< OUT: final covering-grid of physical points */
		    const gsl_matrix  *generator,	/**< IN: _SQUARE_ generating matrix for lattice*/
		    const REAL8Vector *startPoint, 	/**< IN: physical startpoint for covering grid*/
		    BOOLEAN (*isInside)(const REAL8Vector *point) /**< IN: boundary-condition */
		    )
{
  UINT4 dim;	/**< dimension of parameter-space to cover */
  INT4VectorList openEnds = empty_INT4VectorList;	/**< list of "open ends", i.e. unhandled lattice-points */
  INT4VectorList gridPoints = empty_INT4VectorList;	/**< list of lattice-coordinates of resulting grid */

  INT4Vector  *latticePoint = NULL;		/* lattice-coordinates (Z^N) */
  REAL8Vector *physicalPoint = NULL;		/* physical coordinates (R^N) */

  INITSTATUS( lstat, "LALLatticeCovering", LATTICECOVERINGC );
  ATTATCHSTATUSPTR (lstat); 

  /* This traps coding errors in the calling routine. */
  ASSERT ( coveringGrid != NULL, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );  
  ASSERT ( *coveringGrid == NULL,lstat, LATTICECOVERING_ENONULL, LATTICECOVERING_MSGENONULL );
  ASSERT ( generator, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( startPoint, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( startPoint->data, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );

  if ( generator->size1 != generator->size2 )	/* need square generator */
    {
      LALPrintError ("\nERROR: LatticeCovering() required SQUARE generating matrix!\n\n");
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
  if ( NULL == append2INT4VectorList (&openEnds, latticePoint)) {	/* NOTE: head always stays empty for simplicity! */
    LALPrintError ("\nERROR: append2INT4VectorList() failed!\n\n");
    ABORT (lstat, LATTICECOVERING_ELIST, LATTICECOVERING_MSGELIST);
  }

   
  /* ----- (*) take coordinates of next open-end from hash-list of open-ends */
  while ( openEnds.next )
    {
      INT4Vector *p1 = NULL;

      /* get lattice-coordinates of this point (pointer to entry) */
      p1 = &(openEnds.next->entry);
      /* find its physical coordinates */
      XLALlatticePoint2physicalPoint ( physicalPoint, p1, generator, startPoint );
      if ( xlalErrno )
	{
	  int code = xlalErrno;
	  XLALClearErrno(); 
	  LALPrintError ("\nERROR: latticePoint2physicalPoint() failed with xlalErrno = %d!\n\n", code);
	  ABORT (lstat, LATTICECOVERING_EFUNC, LATTICECOVERING_MSGEFUNC);
	}

      /* is it inside the covering-region?: */

      /* NO: ==> move it to (hashed) list of dead-ends, 
	         OR simply remove it from list */

      /* YES: ==> move it to (hashed) list of lattice-coordinates. 
	         generate coordinates of its 2*dim neighbors,
		 discard all those already in list of open-ends
		 and those already in lattice-points (via hash-algorithm..)
		 add the other ones to hash-list of open-ends
      */

      /* start from (1) until no more open ends left */
    } /* while (openEnds->next) */
  
  /* turn linked list of lattice-coordinates into array of lattice-points (in parameter-coordinates) */
  

  /* clean up */
  DETATCHSTATUSPTR (lstat);

  RETURN( lstat );

} /* LALLatticeCovering() */

/** function to calculate the physical coordinates of a lattice-vector with given
 * generating-matrix and start-point of the lattice 
 */
int
XLALlatticePoint2physicalPoint ( REAL8Vector *physicalPoint, 
				 const INT4Vector *latticePoint, 
				 const gsl_matrix *generator, 
				 const REAL8Vector *startPoint )
{
  UINT4 dim;
  
  /* check validity of input */
  if ( !physicalPoint || latticePoint || !generator || !startPoint ||
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
 */
int
XLALSquareGeneratingMatrix(gsl_matrix **outmatrix, 	/**< OUT: square generating matrix */
			   const gsl_matrix *inmatrix)	/**< generating matrix */
{
  UINT4 rows, cols;
  gsl_matrix *sq = NULL;	/* output: square generating matrix */
  gsl_matrix *basis = NULL;	/* orthonormal basis */

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
  sq = gsl_matrix_calloc (rows, rows);
  if ( sq == NULL ) {
    XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_ENOMEM);
  }

  /* if input-matrix is quadratic, we're done */
  if ( rows == cols )
    gsl_matrix_memcpy ( sq, inmatrix );
  else
    {
      gsl_matrix *gij = gsl_matrix_alloc( cols, cols );
      if ( !gij ){
	XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_ENOMEM);
      }
      gsl_matrix_set_identity (gij);	/* use Euklidean metric for orthonormalization*/
      
      /* find orthonormal basis */
      if ( XLALMetricGramSchmidt (&basis, inmatrix, gij) ) {
	XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_EFUNC);
      }
      
      /* free metric again */
      gsl_matrix_free (gij);

      printf ("\nOrthonormal basis:\n");
      print_matrix (basis);

      /* express generating matrix in this new basis: inmatrix.basis^T */
      if ( gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, inmatrix, basis, 0.0, sq)) {
	LALPrintError ("ERROR: Call to  gsl_blas_dgemm() failed\n");
	XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_EFUNC);
      }

    } /* if cols > rows  */

  /* we're done: return result */
  (*outmatrix) = sq;

  return 0;
} /* XLALSquareGeneratingMatrix() */


/** Return the (canonical, i.e. not necessarily quadratic) n-dimensional
 * generating matrix for various lattices. The memory for the matrix
 * is allocated here. 
 * 
 */
int
XLALGetGeneratingMatrix (gsl_matrix **outmatrix,	/**< OUT: generating matrix */
			 UINT4 dimension,		/**< IN: number of dimensions */
			 LatticeType type)		/**< IN: type of lattice */
{
  gsl_matrix *generator = NULL;	/* output: generating matrix */

  /* check that output 'outmatrix' points to a NULL-vector! */
  if ( *outmatrix != NULL ) 
    {
      LALPrintError ("Output-vector not set to NULL");
      XLAL_ERROR("XLALSquareGeneratingMatrix", XLAL_EINVAL);
    }

  switch (type)
    {
    case LATTICE_TYPE_CUBIC:
      generator = gsl_matrix_calloc( dimension, dimension );
      gsl_matrix_set_identity (generator);	/* trivial generating matrix */
      break;

    case LATTICE_TYPE_ANSTAR:
      /* the generating matrix for An* can be written as:
       * [FIXME: check that ]
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

      break;

    default:
      LALPrintError ("Illegal value for lattice-type (%d)\n", type);
      XLAL_ERROR("XLALGetGeneratingMatrix", XLAL_EINVAL);
      break;
      
    } /* switch(type) */

  /* return generating matrix */
  (*outmatrix) = generator;

  return 0;
} /* XLALGetGeneratingMatrix() */


/*----------------------------------------------------------------------
 * list handling tools
 *----------------------------------------------------------------------*/

/** append given INT4Vector element to the list 'head', return pointer to new list-entry.
 * 
 * NOTE: this function is rather permissive in that it takes the dimension from
 * the new element 'el', without requiring this to be equal to the dimension of the 
 * previous list-entries...
 */
INT4VectorList *
append2INT4VectorList (INT4VectorList *head, const INT4Vector *el)
{
  UINT4 dim;
  INT4VectorList *ptr = NULL;	/* running list-pointer */
  INT4VectorList *newElement = NULL;	/* new list-element */
  /* check illegal input */
  if ( (head == NULL) || (el == NULL) )
    return NULL;

  /* find tail of list */
  ptr = head;
  while ( ptr->next )
    ptr = ptr->next;

  /* construct new list-entry */
  dim = el->length;
  if ( (newElement = LALCalloc (1, sizeof (*newElement))) == NULL)
    return NULL;
  if ( (newElement->entry.data = LALCalloc (dim, sizeof(el->data[0]))) == NULL ) {
    LALFree (newElement);
    return NULL;
  }
  newElement->entry.length = dim;
  memcpy (newElement->entry.data, el->data, dim * sizeof(el->data[0]) );

  /* link this to the tail of list */
  ptr->next = newElement;
  newElement->prev = ptr;
  
  return newElement;

} /* append2INT4VectorList() */

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
  
  testGS();

  return 0;
}

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
