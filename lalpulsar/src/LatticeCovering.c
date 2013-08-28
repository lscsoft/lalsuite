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
 * \addtogroup LatticeCovering_h
 * \brief Functions for covering metric parameter-spaces with a lattice.
 *
 * \todo
 * - iterative generation of lattice
 * - optimize for speed
 */

/*---------- INCLUDES ----------*/
#include <math.h>

#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "LatticeCovering.h"

/*---------- DEFINES ----------*/
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

#define TRUE (1==1)
#define FALSE (1==0)

/*----- SWITCHES -----*/
#define HAVE_INLINE	1
/* uncomment the following to turn off range-checking in GSL vector-functions */
/* #define GSL_RANGE_CHECK_OFF 1*/

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */
INT4VectorList empty_INT4VectorList;
REAL8VectorList empty_REAL8VectorList;

/*---------- Global variables ----------*/


/*---------- internal prototypes ----------*/
static int XLALlatticePoint2physicalPoint ( REAL8Vector *physicalPoint, const INT4Vector *latticePoint,
				     const gsl_matrix *generator, const REAL8Vector *startPoint );

/* functions to deal with a non-unity metric */
static int XLALMetricGramSchmidt(gsl_matrix **orth, const gsl_matrix *colvect, const gsl_matrix *metric);

/* list-handling functions */
static INT4VectorList* INT4VectorListAddEntry (INT4VectorList *head, const INT4Vector *entry);
static int INT4VectorListRelinkElement (INT4VectorList *head, INT4VectorList *element);
static void INT4VectorListRemoveElement (INT4VectorList *element);
static BOOLEAN isINT4PointInList ( INT4Vector *point, INT4VectorList *list );
static void INT4VectorListDestroy (INT4VectorList *head);

/* misc helper functions */
static BOOLEAN isSymmetric (const gsl_matrix *Sij);

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * Central function of this module: produce a lattice-covering
 * of given lattice-type for the given parameter-space with constant
 * metric.
 *
 * For optimal covering, use latticeType=0, namely the An* lattice,
 * which is the best known covering-lattice up to dimension 23,
 * see \ref CS99 "[CS99]"
 *
 * ### Algorithm: ###
 *
 * \li 1) use XLALFindCoveringGenerator() to get the generator
 * of the given lattice.
 *
 * \li 2) use it to LALLatticeFill() the space.
 *
 */
void
LALLatticeCovering (LALStatus *status,			/**< pointer to LALStatus structure */
		    REAL8VectorList **covering,		/**< [out] final covering-grid */
		    REAL8 coveringRadius,		/**< [in] covering radius */
		    const gsl_matrix *metric,		/**< [in] constant metric */
		    const REAL8Vector *startPoint,	/**< [in] start-point in the covering-region */
		    BOOLEAN (*isInside)(const REAL8Vector *point), /**< [in] boundary-condition */
		    LatticeType latticeType 		/**< [in] type of lattice to construct */
                    )
{
  gsl_matrix *generator = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* ----- Check validity of input params */
  ASSERT ( covering != NULL, status, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( *covering == NULL,status, LATTICECOVERING_ENONULL, LATTICECOVERING_MSGENONULL );
  ASSERT ( metric, status, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( metric->size1 == metric->size2, status, LATTICECOVERING_EINPUT, LATTICECOVERING_MSGEINPUT );
  ASSERT ( startPoint, status, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( startPoint->data, status, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );

  /* check that startPoint has dimensions consistent with metric */
  ASSERT ( metric->size1 == startPoint->length, status, LATTICECOVERING_EINPUT, LATTICECOVERING_MSGEINPUT);

  /* 1) ----- get the generating matrix for a properly scaled An* lattice */
  if (XLALFindCoveringGenerator (&generator, latticeType, coveringRadius, metric ) < 0)
    {
      int code = xlalErrno;
      XLALClearErrno();
      XLALPrintError ("\nERROR: XLALFindCoveringGenerator() failed (xlalErrno = %d)!\n\n", code);
      ABORT (status, LATTICECOVERING_EFUNC, LATTICECOVERING_MSGEFUNC);
    }

  /* 2) ----- fill parameter-space with this lattice */
  LALLatticeFill(status->statusPtr, covering, generator, startPoint, isInside );
  BEGINFAIL (status) {
    gsl_matrix_free (generator);
  } ENDFAIL(status);

  /* free memory */
  gsl_matrix_free (generator);

  DETATCHSTATUSPTR (status);
  RETURN( status );

} /* LALLatticeCovering() */

/**
 * Fill the given parameter-space by a lattice defined by the specified
 * generating matrix.
 *
 * ### Note1: ###
 *
 * The input generating-matrix (generator) must already be scaled
 * correctly to the required covering radius, also, it needs to be in
 * canonical full-rank square matrix form.
 *
 * ### Note2: ###
 *
 * As always in this module, the generating matrix contains the lattice-vectors as <em>rows</em>
 *
 */
void
LALLatticeFill (LALStatus *status,		/**< pointer to LALStatus structure */
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

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* Check input validity */
  ASSERT ( fillGrid != NULL, status, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( *fillGrid == NULL,status, LATTICECOVERING_ENONULL, LATTICECOVERING_MSGENONULL );
  ASSERT ( generator, status, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( startPoint, status, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );
  ASSERT ( startPoint->data, status, LATTICECOVERING_ENULL, LATTICECOVERING_MSGENULL );

  if ( generator->size1 != generator->size2 )	/* need square generator */
    {
      XLALPrintError ("\nERROR: LatticeFill() requires a  SQUARE generating matrix!\n\n");
      ABORT (status, LATTICECOVERING_EINPUT, LATTICECOVERING_MSGEINPUT);
    }

  if ( ! (*isInside)(startPoint) )	/* startPoint has to be inside */
    {
      XLALPrintError ("\nERROR: startPoint must lie withing the covering-region!\n\n");
      ABORT (status, LATTICECOVERING_EINPUT, LATTICECOVERING_MSGEINPUT);
    }

  dim = startPoint->length;	/* dimension of parameter-space to cover */

  if ( (generator->size1 != dim) )
    {
      XLALPrintError ("\nERROR: all input-dimensions must agree (dim=%d)\n\n", dim);
      ABORT (status, LATTICECOVERING_EINPUT, LATTICECOVERING_MSGEINPUT);
    }

  /* ---------- prepare memory for one grid-coordinate and one physical grid-point */
  TRY ( LALI4CreateVector (status->statusPtr, &latticePoint, dim), status);
  LALDCreateVector (status->statusPtr, &physicalPoint, dim);
  BEGINFAIL (status) {
    TRY ( LALI4DestroyVector(status->statusPtr, &latticePoint), status);
  } ENDFAIL(status);
  /* initialize vectors to 0 */
  memset ( latticePoint->data, 0, dim * sizeof (latticePoint->data[0]) );
  memset ( physicalPoint->data, 0, dim * sizeof (physicalPoint->data[0]));


  /* ----- start by adding startPoint coordinates (0,0,0,,) to the list of 'open-ends' */
  /* we don't need to set these coordinates, grid-point is already set to (0,0,0,,,)  */
  if ( NULL == INT4VectorListAddEntry (&openEnds, latticePoint))
    {
      /* NOTE: head always stays empty for simplicity! */
      XLALPrintError ("\nERROR: INT4VectorListAddEntry () failed!\n\n");
      ABORT (status, LATTICECOVERING_ELIST, LATTICECOVERING_MSGELIST);
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
	  XLALPrintError ("\nERROR: latticePoint2physicalPoint() failed (xlalErrno = %d)!\n\n", code);
	  ABORT (status, LATTICECOVERING_EFUNC, LATTICECOVERING_MSGEFUNC);
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
	  if ( NULL == XLALREAL8VectorListAddEntry ( &realPoints, physicalPoint) )
	    {
	      /* NOTE: head always stays empty for simplicity! */
	      XLALPrintError ("\nERROR: REAL8VectorListAddEntry () failed!\n\n");
	      ABORT (status, LATTICECOVERING_ELIST, LATTICECOVERING_MSGELIST);
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
		  XLALPrintError ("\nERROR: REAL8VectorListAddEntry () failed!\n\n");
		  ABORT (status, LATTICECOVERING_ELIST, LATTICECOVERING_MSGELIST);
		}

	    } /* for i < 2 * dim */

	} /* if point inside filling-region */


      /* start from (1) until no more open ends left */
    } /* while (openEnds->next) */

  /* return linked list of physical points */
  (*fillGrid) = realPoints.next;

  /* clean up memory */
  INT4VectorListDestroy ( gridPoints.next );	/* free list of lattice-coordinates */
  TRY ( LALI4DestroyVector(status->statusPtr, &latticePoint), status);
  TRY ( LALDDestroyVector(status->statusPtr, &physicalPoint), status);

  DETATCHSTATUSPTR (status);

  RETURN( status );

} /* LALLatticeFill() */

/**
 * Calculate the physical coordinates \f$v^i\f$ of a lattice-vector \f$w^i\f$ for given
 * generating-matrix \f${M_i}^j\f$ and start-point \f$p^i\f$ of the lattice.
 *
 * The algorithm is simply: \f$v^i = p^i + \sum_{l=0} w^l \lambda_{(l)}^i\f$,
 * where \f$\vec{\lambda}_{(l)}\f$ is the l-th lattice-vector, which is stored as
 * the l-th row of the generating matrix, i.e.
 * \f${M_i}^j = {\lambda_{(i)}}^j\f$
 *
 * \par Note:
 * The memory for physicalPoint needs to be allocated already, and the
 * dimensions of all vectors and matrices passed to this functions must agree!
 */
int
XLALlatticePoint2physicalPoint ( REAL8Vector *physicalPoint, 	/**< [out] physical coordinates */
				 const INT4Vector *latticePoint,/**< [in] lattice-coordinates */
				 const gsl_matrix *generator, 	/**< [in] generating vectors as rows */
				 const REAL8Vector *startPoint )/**< [in] phys.coordinates of (0,0,..0)*/
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
      XLALPrintError ("\nNULL Input received!\n\n");
      XLAL_ERROR ( XLAL_EINVAL);
    }

  dim = physicalPoint->length;
  if ( (latticePoint->length != dim) || (generator->size1 != dim) || (generator->size2 != dim) )
    {
      XLALPrintError ("\nInconsistent dimensions in input-vectors/matrices!\n\n");
      XLAL_ERROR ( XLAL_EINVAL);
    }

  if ( (buffer = gsl_matrix_calloc (dim, dim)) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM);
  }

  /* create a local copy of the generating-matrix */
  gsl_matrix_memcpy ( buffer, generator );

  /* create a gsl-vector for summing up the lat^l basis_l vectors (->final result) */
  if ( (res = gsl_vector_calloc (dim)) == NULL ) {
    gsl_matrix_free (buffer);
    XLAL_ERROR ( XLAL_ENOMEM);
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

/**
 * Scalar product of two vectors with respect to the given metric
 * \f$\vec{v}_1 \cdot \vec{v}_2 = g_{i j}\, v_1^i \,v_2^j \f$.
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
    XLALPrintError ("\nNULL Input received.\n\n");
    XLAL_ERROR_REAL8(XLAL_EINVAL);
  }

  /* check that gij is symmetric */
  if ( !isSymmetric(gij) ) {
    XLALPrintError ("\nInput 'metric' has to be symmetric!\n\n");
    XLAL_ERROR_REAL8(XLAL_EINVAL);
  }

  dim = gij->size1;

  /* check that vectors have correct sizes */
  if ( (v1->size != dim) || (v2->size != dim) ) {
    XLALPrintError ("\nVectors v1, v2 must have same dimension as metric gij\n\n");
    XLAL_ERROR_REAL8(XLAL_EINVAL);
  }

  /* calculate scalar product */
  prod = 0;
  for ( i=0; i < dim; i++)
    for (j=0; j < dim; j++)
      prod += gsl_vector_get (v1, i) * gsl_matrix_get (gij, i, j) * gsl_vector_get (v2, j);

  return prod;

} /* XLALMetricScalarProduct() */


/**
 * Gram-Schmidt orthogonalization of linearly indepenent vectors using a given metric.
 * As usual the vectors in input and output are stored as the matrix-rows!
 *
 * ### Note1: ###
 *
 * this is a straightforward, probably naive implementation of the basic
 * algorithm, completely ignorant about numerically more robust or faster algorithms
 * to do this... [FIXME?].
 *
 * ### Note2: ###
 *
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
  gsl_vector *para = NULL;

  /* check NULL-vectors on input */
  if ( (!invects) || (!gij) || (!outvects) ) {
    XLALPrintError ("\nNULL Input received.\n\n");
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* check that output 'outvects' points to a NULL-vector! */
  if ( *outvects != NULL ) {
    XLALPrintError ("\nOutput-vector not set to NULL\n\n");
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* check that gij is symmetric */
  if ( !isSymmetric(gij) ) {
    XLALPrintError ("\nInput 'metric' has to be symmetric!\n\n");
    XLAL_ERROR(XLAL_EINVAL);
  }

  numvects = invects->size1;	/* number of rows! */
  vectdim = invects->size2;	/* number of columns */

  /* can't have more vectors than dimensions */
  if ( numvects > vectdim ) {
    XLALPrintError ("\nInput vectors are not linearly independent\n\n");
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* vector-dimension has to be consistent with metric */
  if ( vectdim != gij->size1 ) {
    XLALPrintError ("\nDimension of input vectors inconsistent with metric\n\n");
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* prepare output-matrix for orthonormalized vectors in rows*/
  orth = gsl_matrix_calloc ( numvects, vectdim);
  if ( orth == NULL ) {
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* prepare vector view-arraw on final orthonormal vectors */
  ui = LALCalloc (1, numvects * sizeof(gsl_vector_view) );
  for (i=0; i < numvects; i++)
    ui[i] = gsl_matrix_row ( orth, i );


  /* placeholder for temporary vector */
  para = gsl_vector_alloc(vectdim);
  if ( para == NULL ) {
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /*---------- main algorithm ---------- */
  for (i=0; i < numvects; i++)
    {
      REAL8 norm;
      /*---------- Step 1) orthogonalize wrt to previous vectors ----------*/

      /* view on input-vector i (convenience) */
      gsl_vector_const_view vi = gsl_matrix_const_row (invects, i);

      gsl_vector_memcpy ( &(ui[i].vector), &(vi.vector) );

      for (j=0; (i>0) && (j < i); j++)
	{
	  REAL8 proj;

	  gsl_vector_memcpy ( para, &(ui[j].vector) );

	  proj = XLALMetricScalarProduct(&vi.vector, &(ui[j].vector), gij);
	  if( gsl_vector_scale ( para, proj) ) {	/* para_j = <vi,uj> uj */
	    XLAL_ERROR(XLAL_EFUNC);
	  }

	  if ( gsl_vector_sub ( &(ui[i].vector), para) ) {	/* ui -= para_j */
	    XLAL_ERROR(XLAL_EFUNC);
	  }

	} /* for j < i-1 */

      /*---------- Step 2) normalize ---------- */

      norm = XLALMetricScalarProduct ( &(ui[i].vector), &(ui[i].vector), gij );
      norm = sqrt(norm);

      if ( gsl_vector_scale ( &(ui[i].vector), 1.0/norm ) ) {
	XLAL_ERROR(XLAL_EFUNC);
      }

    } /* for i < numvects */


  /* free memory */
  gsl_vector_free (para);
  LALFree (ui);

  /* return result */
  (*outvects) = orth;

  return 0;

} /* XLALMetricGramSchmidt() */

/**
 * Construct the full-rank generating matrix of given type, for covering
 * a space of given (constant) metric and a given covering Radius.
 *
 * ### Note1: ###
 *
 * the returned generator is a square matrix, with lattice-vectors
 * in the matrix-rows (as always in this module!)
 *
 * ### Note2: ###
 *
 * the memory for 'outmatrix' is allocated here with gsl_matrix_alloc()
 * and needs to be free by the caller via gsl_matrix_free()
 *
 * ### Algorithm: ###
 *
 * \li 1) get the (generally non-square) generating matrix
 * for the lattice
 * \li 2) reduce it to a full-rank square generating matrix by expressing the
 * lattice vectors in the basis of their spanned space
 * \li 3) translate the basis-vectors to the coordinates cooresponding
 * to the given metric
 * \li 4) scale the generating matrix to the given covering radius
 *
 */
int
XLALFindCoveringGenerator (gsl_matrix **outmatrix, /**< [out] generating matrix for covering lattice */
			  LatticeType type,	/**< [in] type of lattice */
			  REAL8 coveringRadius,	/**< [in] desired covering radius */
			  const gsl_matrix *gij	/**< [in] (constant) metric of covering space */
			  )
{
  UINT4 dim;
  gsl_matrix *generator0 = NULL;
  gsl_matrix *generator1 = NULL;
  gsl_matrix *generator2 = NULL;
  gsl_matrix *basis = NULL;

  /* check validity of input */
  if ( !outmatrix || !gij ) {
      XLALPrintError ("\nERROR: NULL Input \n\n");
      XLAL_ERROR(XLAL_EINVAL);
  }
  if ( *outmatrix != NULL ) {
    XLALPrintError ("\nERROR: Output matrix not set to NULL\n\n");
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ( ! isSymmetric (gij) ) {
    XLALPrintError ("\nERROR: metric is not symmetric!!\n\n");
    XLAL_ERROR(XLAL_EINVAL);
  }

  dim = gij->size1;

  /* 1) */
  XLALGetLatticeGenerator (&generator0, dim, type);/* 'raw' generator (not necessarily full-rank)*/

  /* 2) */
  XLALReduceGenerator2FullRank( &generator1, generator0 ); /* full-rank, but in Cartesian coordinates */
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
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* ----- express generating matrix in this new Euklidean basis:
   * generator' = generator . basis
   * where basis is the row-matrix of unit basis-vectors of the new
   * orthonormal basis in 'old' coordinates
   */
  if ( (generator2 = gsl_matrix_calloc (dim, dim)) == NULL ) {
    XLAL_ERROR(XLAL_ENOMEM);
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
      XLALPrintError ("\nERROR: Call to  gsl_blas_dgemm() failed\n\n");
      XLAL_ERROR(XLAL_EFUNC);
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

/**
 * "Reduce" a general (non-quadratic) generating matrix M with rank(M) <= cols(M)
 * into a quadratic generator of full rank.
 *
 * The input matrix can have columns >= rows, the rows reprenting the lattice vectors.
 * This algorithm simply proceeds by constructing an (Euclidean!) orthonormal basis out
 * of the lattice vectors (using GramSchmidt), and then expressing the lattice-vectors
 * in this new basis.
 *
 * \par Note:
 * the memory for 'outmatrix' is allocated in here via gsl_matrix_alloc()
 * and has to be free'ed by the caller via gsl_matrix_free() !
 */
int
XLALReduceGenerator2FullRank (gsl_matrix **outmatrix, 	/**< [out] full-rank square generating matrix */
			      const gsl_matrix *inmatrix)/**< [in] generating matrix (cols >= rows) */
{
  UINT4 rows, cols;
  gsl_matrix *sq = NULL;	/* output: square generating matrix */
  gsl_matrix *basis = NULL;	/* orthonormal basis */

  /* check NULL-vectors on input */
  if ( inmatrix == NULL )
    {
      XLALPrintError ("\nNULL Input received.\n\n");
      XLAL_ERROR(XLAL_EINVAL);
    }

  /* check that output 'outmatrix' points to a NULL-vector! */
  if ( *outmatrix != NULL )
    {
      XLALPrintError ("\nOutput-vector not set to NULL\n\n");
      XLAL_ERROR(XLAL_EINVAL);
    }

  rows = inmatrix->size1;
  cols = inmatrix->size2;

  /* rows need to be lattice vectors, and linearly independent */
  if ( rows > cols )
    {
      XLALPrintError ("\nERROR: input-matrix must have full row-rank!\n\n");
      XLAL_ERROR(XLAL_EINVAL);
    }

  /* allocate output matrix */
  if ( (sq = gsl_matrix_calloc (rows, rows)) == NULL ) {
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* if input-matrix is quadratic, we're done */
  if ( rows == cols )
    gsl_matrix_memcpy ( sq, inmatrix );
  else
    {
      gsl_matrix *gij = gsl_matrix_alloc( cols, cols );
      if ( !gij ){
	XLAL_ERROR(XLAL_ENOMEM);
      }
      gsl_matrix_set_identity (gij);	/* use Euklidean metric for orthonormalization*/

      /* now express generator in 'canonical coordinates' wrt. the Euklidean metric,
       * i.e. express the coordinates of the lattice-vectors in an orthonormal basis
       * spanning exactly the space spanned by the lattice.
       * ==> The resulting generator will therefore be a full-rank square matrix.
       */

      /* ----- find orthonormal basis in the lattice's space*/
      if ( XLALMetricGramSchmidt (&basis, inmatrix, gij) ) {
	XLAL_ERROR(XLAL_EFUNC);
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
	  XLALPrintError ("\nERROR: Call to  gsl_blas_dgemm() failed\n\n");
	  XLAL_ERROR(XLAL_EFUNC);
	}

      /* free memory */
      gsl_matrix_free (basis);
      gsl_matrix_free (gij);

    } /* if cols > rows  */

  /* we're done: return result */
  (*outmatrix) = sq;

  return 0;
} /* XLALReduceGenerator2FullRank() */

/**
 * Return a (not necessarily quadratic) n-dimensional generating matrix
 * for one of several possible lattices (currently possible: cubic or \f$A_n^*\f$).
 * See \ref CS99 "[CS99]" for the definition and properties of these lattices.
 *
 * ### Note1: ###
 *
 * Because these lattices are intended for covering, we scale
 * them so that their covering-radius is unity.
 * This allows the user to later-on scale these easily to any
 * desired covering-radius without having to know anything about the lattice...
 * (Remembering that if you scale the generator \f$M\f$ of a lattice by \f$M' = c M\f$,
 * then the covering radius R scales as \f$R' = c R\f$)
 *
 * ### Note2: ###
 *
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
  UINT4 row, col;

  /* check that output 'outmatrix' points to a NULL-vector! */
  if ( *outmatrix != NULL )
    {
      XLALPrintError ("\nOutput-vector not set to NULL\n\n");
      XLAL_ERROR(XLAL_EINVAL);
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
      XLALPrintError ("\nIllegal value for lattice-type (%d)\n\n", type);
      XLAL_ERROR(XLAL_EINVAL);
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

/**
 * Add a new element at the end of the list 'head', _copy_ the given entry there,
 * and return pointer to the new list-entry.
 *
 * \par Note:
 * This function is rather permissive in that it takes the vector-dimension
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


/**
 * Add a new element at the end of the list 'head', _copy_ the given entry there,
 * and return pointer to the new list-entry.
 *
 * \par Note:
 * This function is rather permissive in that it takes the vector-dimension
 * from the new entry 'el', without requiring this to be equal to the dimension of the
 * other list-entries...
 */
REAL8VectorList *
XLALREAL8VectorListAddEntry (REAL8VectorList *head, const REAL8Vector *entry)
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

} /* XLALREAL8VectorListAddEntry() */


/**
 * "relink" the given element (from whatever list) to the end of the given list 'head'
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

/**
 * 'List-destructor' for INT4VectorList: free a complete list.
 *
 * \par Note:
 * 'head' will be freed too, so make
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

/**
 * 'List-destructor' for REAL8VectorList: free a complete list.
 *
 * \par Note:
 * 'head' will be freed too, so make
 * sure not to pass a non-freeable head (like an automatic variabe)
 */
void
XLALREAL8VectorListDestroy (REAL8VectorList *head)
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
} /* XLALREAL8VectorListDestroy() */




/**
 * Search the given list (or hash-table?) for the given point.
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

/**
 * Translate a symmetric gsl_matrix into a 'LAL-encoded' REAL8Vector, using
 * the index-convention l = a + b*(b+1) if a <= b, see PMETRIC_INDEX(a,b).
 */
REAL8Vector *
XLALgsl2LALmetric (const gsl_matrix *gmetric)
{
  REAL8Vector *metric = NULL;
  UINT4 dim, length, i, j;

  if ( gmetric == NULL ) {
    XLALPrintError ("\nNULL Input received!\n\n");
    XLAL_ERROR_NULL ( XLAL_EINVAL);
  }

  if ( !isSymmetric(gmetric) ) {
    XLALPrintError ("\nInput matrix is not symmetric!\n\n");
    XLAL_ERROR_NULL ( XLAL_EINVAL);
  }

  dim = gmetric->size1;
  length = dim * (dim+1) /2;	/* independent elements in symmetric  matrix */

  metric = XLALCreateREAL8Vector ( length );
  if ( ! metric ) {
    XLAL_ERROR_NULL( XLAL_EFUNC);
  }

  for (i=0; i < dim; i++)
    for (j=i; j < dim; j++)
       metric->data[ PMETRIC_INDEX(i,j) ] = gsl_matrix_get (gmetric, i, j);

  return metric;

} /* XLALgsl2LALmetric() */

/**
 * Convert a LAL-encoded metric (REAL8Vector) into a symmetric gsl_matrix
 */
gsl_matrix *
XLALmetric2gsl (const REAL8Vector *metric)
{
  gsl_matrix *gij;
  INT4 dim;
  INT4 i,j;

  if ( !metric )
    XLAL_ERROR_NULL (XLAL_EINVAL);

  if ( (dim = XLALFindMetricDim ( metric )) <= 0 )
    XLAL_ERROR_NULL (XLAL_EFUNC);

  if ( (gij = gsl_matrix_calloc( dim, dim )) == NULL )
    XLAL_ERROR_NULL (XLAL_ENOMEM);

  for (i=0; i < dim; i++ )
    for (j=0; j < dim; j++ )
      gsl_matrix_set (gij, i, j, metric->data[PMETRIC_INDEX(i,j)] );

  return gij;

} /* XLALmetric2gsl() */
