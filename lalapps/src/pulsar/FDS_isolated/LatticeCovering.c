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


#include <lal/LALStdlib.h>

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
#define LATTICECOVERING_ENONULL	2
#define LATTICECOVERING_EMEM		3
#define LATTICECOVERING_EINPUT		4

#define LATTICECOVERING_MSGENULL 	"Arguments contained an unexpected null pointer"
#define LATTICECOVERING_MSGENONULL	"Output pointer is not NULL"
#define LATTICECOVERING_MSGEMEM	"Out of memory"
#define LATTICECOVERING_MSGEINPUT	"Invald input parameter"

/*</lalErrTable> */

/*---------- local types ----------*/



/*---------- Global variables ----------*/



/*---------- local prototypes ----------*/
REAL8 XLALMetricScalarProduct (const gsl_vector *vector1,
			       const gsl_vector *vector2,	
			       const gsl_matrix *metric);


int
XLALMetricGramSchmidt(gsl_matrix **orth, 		
		      const gsl_matrix *colvect,	
		      const gsl_matrix *metric);

BOOLEAN isSymmetric (const gsl_matrix *Sij);

/*==================== FUNCTION DEFINITIONS ====================*/

/* NOTE: "vector" in the following refers to the covariant compomnents,
 * while the covariant-components would be called the "covector"
 */


/** Scalar product of two vectors with respect to given metric.
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
 * NOTE: this is a straightforward, naive implementation of the basic
 * algorithm, completely ignorant about numerically more robust or faster algorithms
 * to do this... [FIXME?].
 *
 * The set of vectors in input and output consists of the matrix-columns of colvect and orth
 *
 */
int
XLALMetricGramSchmidt(gsl_matrix **outvects,	/**< OUT: orthonormal contravariant vects */
		      const gsl_matrix *invects,/**< matrix of column-vectors */
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

  numvects = invects->size2;	/* number of columns! */
  vectdim = invects->size1;

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


  /* prepare output-matrix for orthonomalized vectors */
  orth = gsl_matrix_calloc ( vectdim, numvects);
  if ( orth == NULL ) {
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_ENOMEM);
  }

  /* prepare vector view-arraw on final orthonormal vectors */
  ui = LALCalloc (1, numvects * sizeof(gsl_vector_view) );
  for (i=0; i < numvects; i++)
    ui[i] = gsl_matrix_column ( orth, i );


  /* placeholder for temporary vector */
  gsl_vector *para = gsl_vector_alloc(vectdim);
  if ( para == NULL ) {
    XLAL_ERROR("XLALMetricGramSchmidt", XLAL_ENOMEM);
  }

  /*---------- main algorithm ---------- */
  for (i=0; i < numvects; i++)
    {

      /*---------- Step 1) orthogonalize wrt to previous vectors ----------*/
      
      /* view on input-vector i (convenience) */
      gsl_vector_const_view vi = gsl_matrix_const_column (invects, i);

      gsl_vector_memcpy ( &(ui[i].vector), &(vi.vector) );
      
      for (j=0; (i>0) && (j < i-1); j++)
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
      
    } /* for i < numvects */

  gsl_vector_free (para);
  LALFree (ui);

  /* return result */
  (*outvects) = orth;

  return 0;

} /* XLALMetricGramSchmidt() */


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
