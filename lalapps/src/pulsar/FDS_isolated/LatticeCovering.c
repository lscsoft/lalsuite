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
 * \defgroup moduleLatticeCovering Lattice Covering
 *  Covering of parameter-spaces by an optimal lattice
 */

/**
 * \author Reinhard Prix
 * \date 2005
 * \file 
 * \ingroup moduleLatticeCovering
 * \brief Module for covering metric parameter-spaces with a lattice.
 *
 * These routines should eventually become useable for implementing 
 * a covering that's practical and as "optimal" as possible.
 * Right now this is just a playground...
 *
 * $Id$
 *
 */

/** \page References
 * \anchor CS99
 *  	<b>[CS99]</b> J.H. Conway and N.J.A. Sloane, "Sphere Packings, Lattices and Groups", 
 *  	vol. 290 of Grundlehren der mathematischen Wissenschaften, 
 * 	(Springer, New York, U.S.A., 1999), 3rd edition.
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

/** doubly linked list of REAL8-vectors (physical vectors) */
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
void LALLatticeCovering (LALStatus *lstat, REAL8VectorList **covering, REAL8 coveringRadius, const gsl_matrix *metric,
			 const REAL8Vector *startPoint,	BOOLEAN (*isInside)(const REAL8Vector *point) );

void LALLatticeFill (LALStatus *lstat, REAL8VectorList **fillGrid, const gsl_matrix  *generator,
		     const REAL8Vector *startPoint, BOOLEAN (*isInside)(const REAL8Vector *point) );

int XLALlatticePoint2physicalPoint ( REAL8Vector *physicalPoint, const INT4Vector *latticePoint, 
				     const gsl_matrix *generator, const REAL8Vector *startPoint );

/* functions for handling lattice's generating matrix */
int XLALFindCoveringGenerator (gsl_matrix **outmatrix, LatticeType type, UINT4 dimension, 
			       REAL8 coveringRadius, const gsl_matrix *gij);
int XLALReduceGenerator2FullRank(gsl_matrix **outmatrix, const gsl_matrix *inmatrix);
int XLALGetLatticeGenerator (gsl_matrix **outmatrix, UINT4 dimension, LatticeType type);

/* functions to deal with a non-unity metric */
REAL8 XLALMetricScalarProduct (const gsl_vector *vector1, const gsl_vector *vector2,	
			       const gsl_matrix *metric);
int XLALMetricGramSchmidt(gsl_matrix **orth, const gsl_matrix *colvect, const gsl_matrix *metric);


/* list-handling functions */
INT4VectorList* INT4VectorListAddEntry (INT4VectorList *head, const INT4Vector *entry);
int INT4VectorListRelinkElement (INT4VectorList *head, INT4VectorList *element);
void INT4VectorListRemoveElement (INT4VectorList *element);
REAL8VectorList* REAL8VectorListAddEntry (REAL8VectorList *head, const REAL8Vector *entry);
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

/** Central function of this module: produce an optimal covering
 * of the given parameter-space with constant metric.
 *
 * \par Algorithm: 
 * 	\li 1) use XLALFindCoveringGenerator() to get the generator 
 *            of the An* lattice (which is the best known covering-lattice
 *            up to dimension 23, see \ref CS99 "[CS99]"
 *	\li 2) use it to LALLatticeFill() the space.
 *
 */
void
LALLatticeCovering (LALStatus *lstat,
		    REAL8VectorList **covering,		/**< [out] final covering-grid */
		    REAL8 coveringRadius,		/**< [in] covering radius */
		    const gsl_matrix *metric,		/**< [in] constant metric */
		    const REAL8Vector *startPoint,	/**< [in] start-point inside the covering-region */
		    BOOLEAN (*isInside)(const REAL8Vector *point) /**< [in] boundary-condition */
		    )
{
  UINT4 dim;	/* dimension of parameter-space */
  gsl_matrix *generator = NULL;


  INITSTATUS( lstat, "LALLatticeCovering", LATTICECOVERINGC );
  ATTATCHSTATUSPTR (lstat); 


  /* Check validity of input params */
  ASSERT ( covering != NULL, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );  
  ASSERT ( *covering == NULL,lstat, LATTICECOVERING_ENONULL, LATTICECOVERING_MSGENONULL );
  ASSERT ( metric, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( startPoint, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( startPoint->data, lstat, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );

  /* determine dimension of parameter-space from start-Point */
  dim = startPoint->length;

  ASSERT ( (metric->size1 == dim) && ( metric->size2 == dim), lstat, LATTICECOVERING_EINPUT, LATTICECOVERING_MSGEINPUT);


  /* 1) get the generating matrix for a properly scaled An* lattice */
  XLALFindCoveringGenerator (&generator, LATTICE_TYPE_ANSTAR, dim, coveringRadius, metric );
  if ( xlalErrno ) 
    {
      int code = xlalErrno;
      XLALClearErrno(); 
      LALPrintError ("\nERROR: XLALFindCoveringGenerator() failed (xlalErrno = %d)!\n\n", code);
      ABORT (lstat, LATTICECOVERING_EFUNC, LATTICECOVERING_MSGEFUNC);
    }

  /* 2) fill parameter-space with this lattice */
  TRY ( LALLatticeFill(lstat->statusPtr, covering, generator, startPoint, isInside ), lstat );

  /* free memory */
  gsl_matrix_free (generator);

  DETATCHSTATUSPTR (lstat);
  RETURN( lstat );

} /* LALLatticeCovering() */

/** Fill the given parameter-space by a lattice defined by the specified 
 * generating matrix.
 * 
 * \par Note1: 
 * The input generating-matrix (generator) must already be scaled
 * correctly to the required covering radius, also, it needs to be in 
 * canonical full-rank square matrix form.
 *
 * \par Note2: 
 * As always in this module, the generating matrix contains the lattice-vectors as <em>rows</em>
 *
 */
void
LALLatticeFill (LALStatus *lstat,
		REAL8VectorList **fillGrid,	/**< [out] fillGrid final fill-grid (physical points) */
		const gsl_matrix  *generator,	/**< [in] SQUARE generating matrix for lattice*/
		const REAL8Vector *startPoint, 	/**< [in] physical startpoint for filling */
		BOOLEAN (*isInside)(const REAL8Vector *point) /**< [in] boundary-condition */
		)
{
  UINT4 dim;		/* dimension of parameter-space to fill */
  UINT4 i;
  INT4VectorList openEnds = empty_INT4VectorList;	/* list of "open ends" (lattice-points) */
  INT4VectorList gridPoints = empty_INT4VectorList;	/* resulting grid (lattice-points) */
  REAL8VectorList realPoints = empty_REAL8VectorList;	/* physical coordinates of grid-points */

  INT4Vector  *latticePoint = NULL;		/* lattice-coordinates (Z^N) */
  REAL8Vector *physicalPoint = NULL;		/* physical coordinates (R^N) */

  INITSTATUS( lstat, "LALLatticeFill", LATTICECOVERINGC );
  ATTATCHSTATUSPTR (lstat); 

  /* Check input validity */
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

/** Calculate the physical coordinates \f$v^i\f$ of a lattice-vector \f$w^i\f$ for given
 * generating-matrix \f${M_i}^j\f$ and start-point \f$p^i\f$ of the lattice.
 *
 * The algorithm is simply: \f$v^i = p^i + \sum_{l=0} w^l \lambda_{(l)}^i\f$,
 * where \f$\vec{\lambda}_{(l)}\f$ is the l-th lattice-vector, which is stored as
 * the l-th row of the generating matrix, i.e. 
 * \f${M_i}^j = {\lambda_{(i)}}^j\f$
 * 
 * \note The memory for physicalPoint needs to be allocated already, and the 
 * dimensions of all vectors and matrices passed to this functions must agree!
 */
int
XLALlatticePoint2physicalPoint ( REAL8Vector *physicalPoint, 	/**< [out] physical coordinates */
				 const INT4Vector *latticePoint,/**< [in] lattice-coordinates */
				 const gsl_matrix *generator, 	/**< [in] generating vectors as rows */
				 const REAL8Vector *startPoint )/**< [in] phys.coordinates of (0,0,...0)*/
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

/** Scalar product of two vectors with respect to the given metric
 *  \f$\vec{v}_1 \cdot \vec{v}_2 = g_{i j}\, v_1^i \,v_2^j \f$.
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
    LALPrintError ("\nNULL Input received.\n\n");
    XLAL_ERROR_REAL8("XLALMetricScalarProduct", XLAL_EINVAL);
  }

  /* check that gij is symmetric */
  if ( !isSymmetric(gij) ) {
    LALPrintError ("\nInput 'metric' has to be symmetric!\n\n");
    XLAL_ERROR_REAL8("XLALMetricScalarProduct", XLAL_EINVAL);
  }

  dim = gij->size1;

  /* check that vectors have correct sizes */
  if ( (v1->size != dim) || (v2->size != dim) ) {
    LALPrintError ("\nVectors v1, v2 must have same dimension as metric gij\n\n");
    XLAL_ERROR_REAL8("XLALMetricScalarProduct", XLAL_EINVAL);
  }
      
  /* calculate scalar product */
  prod = 0;
  for ( i=0; i < dim; i++)
    for (j=0; j < dim; j++)
      prod += gsl_vector_get (v1, i) * gsl_matrix_get (gij, i, j) * gsl_vector_get (v2, j);

  return prod;

} /* XLALMetricScalarProduct() */


/** Gram-Schmidt orthogonalization of linearly indepenent vectors using a given metric.
 * As usual the vectors in input and output are stored as the matrix-rows!
 * 
 * \par Note1: 
 * this is a straightforward, probably naive implementation of the basic
 * algorithm, completely ignorant about numerically more robust or faster algorithms
 * to do this... [FIXME?].
 *
 * \par Note2: 
 * the memory for outvects is allocated in here by gsl_matrix_alloc()
 * and needs to be free'ed by the caller via gsl_matrix_free() !
 */
int
XLALMetricGramSchmidt(gsl_matrix **outvects,	/**< [out] orthonormal row vects */
		      const gsl_matrix *invects,/**< [in] matrix of row-vectors */
		      const gsl_matrix *gij)	/**< [in] metric */
{
  UINT4 numvects, vectdim;
  UINT4 i, j;
  gsl_matrix *orth = NULL;
  gsl_vector_view *ui = NULL;	/* array of orthonormal result-vectors */

  /* check NULL-vectors on input */
  if ( (!invects) || (!gij) || (!outvects) ) {
    LALPrintError ("\nNULL Input received.\n\n");
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EINVAL);
  }
  
  /* check that output 'outvects' points to a NULL-vector! */
  if ( *outvects != NULL ) {
    LALPrintError ("\nOutput-vector not set to NULL\n\n");
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EINVAL);
  }

  /* check that gij is symmetric */
  if ( !isSymmetric(gij) ) {
    LALPrintError ("\nInput 'metric' has to be symmetric!\n\n");
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EINVAL);
  }

  numvects = invects->size1;	/* number of rows! */
  vectdim = invects->size2;	/* number of columns */

  /* can't have more vectors than dimensions */
  if ( numvects > vectdim ) {
    LALPrintError ("\nInput vectors are not linearly independent\n\n");
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_EINVAL);
  }
  
  /* vector-dimension has to be consistent with metric */
  if ( vectdim != gij->size1 ) {
    LALPrintError ("\nDimension of input vectors inconsistent with metric\n\n");
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

/** Construct the full-rank generating matrix of given type, for covering 
 * a space of given (constant) metric and a given covering Radius.
 * 
 * \par Note1: 
 * the returned generator is a square matrix, with lattice-vectors
 * in the matrix-rows (as always in this module!)
 *
 * \par Note2: 
 * the memory for 'outmatrix' is allocated here with gsl_matrix_alloc()
 * and needs to be free by the caller via gsl_matrix_free()
 *
 * \par Algorithm: 
 * 	\li 1) get the (generally non-square) generating matrix
 *             for the lattice
 *	\li 2) reduce it to a full-rank square generating matrix by expressing the
 *             lattice vectors in the basis of their spanned space
 *	\li 3) translate the basis-vectors to the coordinates cooresponding
 *             to the given metric
 *	\li 4) scale the generating matrix to the given covering radius
 *
 */
int 
XLALFindCoveringGenerator (gsl_matrix **outmatrix, /**< [out] generating matrix for covering lattice */
			  LatticeType type,	/**< [in] type of lattice */
			  UINT4 dimension,	/**< [in] lattice-dimension */
			  REAL8 coveringRadius,	/**< [in] desired covering radius */
			  const gsl_matrix *gij	/**< [in] (constant) metric of covering space */
			  )
{
  gsl_matrix *generator0 = NULL;
  gsl_matrix *generator1 = NULL;
  gsl_matrix *generator2 = NULL;
  gsl_matrix *basis = NULL;

  /* check validity of input */
  if ( !outmatrix || !gij ) {
      LALPrintError ("\nERROR: NULL Input \n\n");
      XLAL_ERROR("XLALFindCoveringGenerator", XLAL_EINVAL);
  }
  if ( *outmatrix != NULL ) {
    LALPrintError ("\nERROR: Output matrix not set to NULL\n\n");
    XLAL_ERROR("XLALFindCoveringGenerator", XLAL_EINVAL);
  }

  /* 1) */
  XLALGetLatticeGenerator (&generator0, dimension, type); /* 'raw' generator (not necessarily full-rank */

  /* 2) */
  XLALReduceGenerator2FullRank( &generator1, generator0 ); /* full-rank, but in unity Cartesian coordinates */
  gsl_matrix_free ( generator0 );
  generator0 = NULL;

  /* 3) now we need to do a coordinate-transformation on the generator
   * to translate it into the coordinates of the parameter-space,
   * which are encoded in the metric gij:
   * for this we construct an orthonormal basis wrt to gij, and then
   * express the generator in this basis.
   */

  /* ----- find orthonormal basis wrt given metric */
  if ( XLALMetricGramSchmidt (&basis, generator1, gij) ) {
    XLAL_ERROR("XLALFindCoveringGenerator", XLAL_EFUNC);
  }
      
  /* ----- express generating matrix in this new Euklidean basis: generator' = generator . basis 
   * where basis is the row-matrix of unit basis-vectors of the new orthonormal basis in 'old' coordinates 
   */
  if ( (generator2 = gsl_matrix_calloc (dimension, dimension)) == NULL ) {
    XLAL_ERROR("XLALFindCoveringGenerator", XLAL_ENOMEM);
  }
  
  /* from the gsl-documentation:
   * int gsl_blas_dgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, 
   *                    double alpha, const gsl_matrix * A, const gsl_matrix * B, 
   *                    double beta, gsl_matrix * C)
   * These functions compute the matrix-matrix product and 
   * sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for 
   * TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
   */
  if ( gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, generator1, basis, 0.0, generator2)) 
    {
      LALPrintError ("\nERROR: Call to  gsl_blas_dgemm() failed\n\n");
      XLAL_ERROR("XLALFindCoveringGenerator", XLAL_EFUNC);
    }
  
  /* 4) ----- finally, scale the generator to the desired covering radius 
   *   (This is simplified by the fact that we scaled the original generating
   *    matrix to R=1)
   */
  gsl_matrix_scale ( generator2, coveringRadius );

  /* free memory */
  gsl_matrix_free (generator1);
  gsl_matrix_free (basis);

  /* return resuling covering-generator */
  (*outmatrix) = generator2;
  
  return 0;

} /* XLALFindCoveringGenerator() */

/** "Reduce" a general (non-quadratic) generating matrix M with rank(M) <= cols(M)
 *  into a quadratic generator of full rank.
 * 
 * The input matrix can have columns >= rows, the rows reprenting the lattice vectors. 
 * This algorithm simply proceeds by constructing an (Euclidean!) orthonormal basis out 
 * of the lattice vectors (using GramSchmidt), and then expressing the lattice-vectors 
 * in this new basis.
 * 
 * \note the memory for 'outmatrix' is allocated in here via gsl_matrix_alloc()
 *       and has to be free'ed by the caller via gsl_matrix_free() !
 */
int
XLALReduceGenerator2FullRank (gsl_matrix **outmatrix, 		/**< [out] full-rank square generating matrix */
			      const gsl_matrix *inmatrix)	/**< [in] generating matrix (generally cols >= rows) */
{
  UINT4 rows, cols;
  gsl_matrix *sq = NULL;	/* output: square generating matrix */
  gsl_matrix *basis = NULL;	/* orthonormal basis */  

  /* check NULL-vectors on input */
  if ( inmatrix == NULL ) 
    {
      LALPrintError ("\nNULL Input received.\n\n");
      XLAL_ERROR("XLALReduceGenerator2FullRank", XLAL_EINVAL);
    }
  
  /* check that output 'outmatrix' points to a NULL-vector! */
  if ( *outmatrix != NULL ) 
    {
      LALPrintError ("\nOutput-vector not set to NULL\n\n");
      XLAL_ERROR("XLALReduceGenerator2FullRank", XLAL_EINVAL);
    }

  rows = inmatrix->size1;
  cols = inmatrix->size2;

  /* rows need to be lattice vectors, and linearly independent */
  if ( rows > cols ) 
    {
      LALPrintError ("\nERROR: input-matrix must have full row-rank!\n\n");
      XLAL_ERROR("XLALReduceGenerator2FullRank", XLAL_EINVAL);
    }

  /* allocate output matrix */
  if ( (sq = gsl_matrix_calloc (rows, rows)) == NULL ) {
    XLAL_ERROR("XLALReduceGenerator2FullRank", XLAL_ENOMEM);
  }

  /* if input-matrix is quadratic, we're done */
  if ( rows == cols ) 
    gsl_matrix_memcpy ( sq, inmatrix );
  else
    {
      gsl_matrix *gij = gsl_matrix_alloc( cols, cols );
      if ( !gij ){
	XLAL_ERROR("XLALReduceGenerator2FullRank", XLAL_ENOMEM);
      }
      gsl_matrix_set_identity (gij);	/* use Euklidean metric for orthonormalization*/

      /* now express generator in 'canonical coordinates' wrt. the Euklidean metric, 
       * i.e. express the coordinates of the lattice-vectors in an orthonormal basis 
       * spanning exactly the space spanned by the lattice. 
       * ==> The resulting generator will therefore be a full-rank square matrix.
       */

      /* ----- find orthonormal basis in the lattice's space*/
      if ( XLALMetricGramSchmidt (&basis, inmatrix, gij) ) {
	XLAL_ERROR("XLALReduceGenerator2FullRank", XLAL_EFUNC);
      }
      
      /* ----- express generating matrix in this new Euklidean basis: inmatrix.basis^T */
  
      /* from the gsl-documentation:
       * int gsl_blas_dgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, 
       *                    double alpha, const gsl_matrix * A, const gsl_matrix * B, 
       *                    double beta, gsl_matrix * C)
       * These functions compute the matrix-matrix product and 
       * sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for 
       * TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
       */
      if ( gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, inmatrix, basis, 0.0, sq)) 
	{
	  LALPrintError ("\nERROR: Call to  gsl_blas_dgemm() failed\n\n");
	  XLAL_ERROR("XLALReduceGenerator2FullRank", XLAL_EFUNC);
	}
  
      /* free memory */
      gsl_matrix_free (basis);
      gsl_matrix_free (gij);

    } /* if cols > rows  */

  /* we're done: return result */
  (*outmatrix) = sq;

  return 0;
} /* XLALReduceGenerator2FullRank() */

/** Return a (not necessarily quadratic) n-dimensional generating matrix 
 *  for one of several possible lattices (currently possible: cubic or \f$A_n^*\f$).
 *  See \ref CS99 "[CS99]" for the definition and properties of these lattices.
 * 
 * \par Note1: 
 * Because these lattices are intended for covering, we scale
 * them so that their covering-radius is unity.
 * This allows the user to later-on scale these easily to any
 * desired covering-radius without having to know anything about the lattice...
 * (Remembering that if you scale the generator \f$M\f$ of a lattice by \f$M' = c M\f$, 
 * then the covering radius R scales as \f$R' = c R\f$)
 *
 * \par Note2: 
 * The memory for 'outmatrix' is allocated in here via gsl_matrix_alloc()
 * and has to be free'ed by the caller via gsl_matrix_free() !
 *
 */
int
XLALGetLatticeGenerator (gsl_matrix **outmatrix,	/**< [out] generating matrix */
			 UINT4 dimension,		/**< [in] number of dimensions */
			 LatticeType type)		/**< [in] type of lattice */
{
  gsl_matrix *generator = NULL;	/* output: generating matrix */
  REAL8 coveringRadius;

  /* check that output 'outmatrix' points to a NULL-vector! */
  if ( *outmatrix != NULL ) 
    {
      LALPrintError ("\nOutput-vector not set to NULL\n\n");
      XLAL_ERROR("XLALGetLatticeGenerator", XLAL_EINVAL);
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

      /* covering Radius of An* is R = sqrt( n*(n+2) / (12*(n+1)) ), see \ref CS99 */
      coveringRadius = sqrt ( 1.0 * dimension * (dimension + 2.0) / (12.0 * (dimension + 1) ));

      break;

    default:
      LALPrintError ("\nIllegal value for lattice-type (%d)\n\n", type);
      XLAL_ERROR("XLALGetLatticeGenerator", XLAL_EINVAL);
      break;
      
    } /* switch(type) */

  /* Now scale the generating matrix by 1/R, such that its covering radius becomes R=1 */
  gsl_matrix_scale ( generator, 1.0 / coveringRadius );

  /* return generating matrix */
  (*outmatrix) = generator;

  return 0;
} /* XLALGetLatticeGenerator() */


/*----------------------------------------------------------------------
 * list handling tools
 *----------------------------------------------------------------------*/

/** Add a new element at the end of the list 'head', _copy_ the given entry there,
 *  and return pointer to the new list-entry. 
 * 
 * \note This function is rather permissive in that it takes the vector-dimension 
 * from the new entry, without requiring this to be equal to the dimension of the 
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


/** Add a new element at the end of the list 'head', _copy_ the given entry there,
 * and return pointer to the new list-entry. 
 * 
 * \note This function is rather permissive in that it takes the vector-dimension 
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


/** Remove (and free) the given list-element from the list */
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

/** 'List-destructor' for INT4VectorList: free a complete list.
 *
 * \note 'head' will be freed too, so make
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

/** 'List-destructor' for REAL8VectorList: free a complete list.
 *
 * \note 'head' will be freed too, so make
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




/** Search the given list (or hash-table?) for the given point.
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

/** Check if matrix is symmetric. */
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


  XLALReduceGenerator2FullRank( &sGM, &(m1.matrix) );

  printf ("Square generating matrix:\n");
  print_matrix (sGM);

  XLALGetLatticeGenerator (&testM, 4, LATTICE_TYPE_ANSTAR);
  printf ("produced Generating Matrix: \n");
  print_matrix (testM);

} /* testGS() */
