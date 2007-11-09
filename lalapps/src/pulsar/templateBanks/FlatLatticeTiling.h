/**
 * \author K. Wette
 * \file
 * \brief Flat lattice tiling over multi-dimensioned parameter spaces
 */

#ifndef _FLATLATTICETILING_H
#define _FLATLATTICETILING_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

#include <lal/LALRCSID.h>
#include <lal/LALDatatypes.h>

/**
 * Structure to hold information about the flat lattice tiling
 */
typedef struct tagFlatLatticeTiling {

  /* Dimension of the parameter space */
  INT4 dimension;

  /* Matrix containing the parameter space metric */
  gsl_matrix *metric;

  /* Maximum mismatch of the parameter space templates */
  REAL8 mismatch;

  /* Generator matrix for the lattice tiling, and its normalised thickness */
  gsl_matrix *generator;

  /* Pointer to function to calculate bounds on parameter space */
  void (*bounds)(gsl_vector*, INT4, gsl_vector*, REAL8*, REAL8*);
  gsl_vector *bounds_args;

  /* Increment vectors between lattice points in the parameter space */
  gsl_matrix *increment;

  /* Current point in parameter space and its upper bound */
  gsl_vector *current;
  gsl_vector *upper;

  /* Flags for filling bounds and point to resume from */
  gsl_vector_int *on_upper;
  gsl_vector_int *on_lower;
  gsl_vector *resume;  

  /* Temporary vector */
  gsl_vector *temp;

  /* Number of templates generated */
  UINT8 templates;

} FlatLatticeTiling;

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID(FLATLATTICETILINGH, "$Id$");

  FlatLatticeTiling *XLALCreateFlatLatticeTiling(INT4);

  void XLALDestroyFlatLatticeTiling(FlatLatticeTiling*);

  int XLALOrthonormaliseMatrixWRTMetric(gsl_matrix*, gsl_matrix*);

  gsl_matrix *XLALLowerTriangularLatticeGenerator(gsl_matrix*);

  int XLALNormaliseLatticeGenerator(gsl_matrix*, REAL8);

  gsl_matrix *XLALCubicLatticeGenerator(INT4);
  
  gsl_matrix *XLALAnstarLatticeGenerator(INT4);

  int XLALSquareParameterSpace(FlatLatticeTiling*, ...);

  int XLALSetupFlatLatticeTiling(FlatLatticeTiling*);

  int XLALNextFlatLatticePoint(FlatLatticeTiling*);

  REAL8 XLALCurrentFlatLatticePoint(FlatLatticeTiling*, INT4);

  UINT8 XLALNumberOfFlatLatticeTiles(FlatLatticeTiling*);

#ifdef __cplusplus
}
#endif

#endif
