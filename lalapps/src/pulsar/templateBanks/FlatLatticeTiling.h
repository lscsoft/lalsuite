/*********************************************************************************
 *  \author K. Wette
 *  \file
 *  \ingroup templateBanks
 *  \brief
 *  Flat lattice tiling over multi-dimensioned parameter spaces
 *********************************************************************************/

#ifndef _FLATLATTICETILING_H
#define _FLATLATTICETILING_H

/******** Includes ********/

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

#include <lal/LALRCSID.h>
#include <lal/LALDatatypes.h>

/******** Structures ********/

typedef struct tagFlatLatticeTiling FlatLatticeTiling;
struct tagFlatLatticeTiling {

  /* Dimension of the parameter space, i.e. number of parameters to be searched over */
  UINT4 dimension;

  /* Are some dimensions degenerate */
  gsl_vector_int *is_degenerate;
  int first_degenerate;

  /* Matrix containing the parameter space metric, and a congruent factor */
  gsl_matrix *metric;
  gsl_matrix *metric_congruent_factor;

  /* Maximum mismatch of the parameter space templates */
  REAL8 mismatch;

  /* Generator matrix for the lattice tiling, and its normalised thickness */
  gsl_matrix *generator;
  REAL8 generator_norm_thickness;

  /* Number of bounds on the parameter space */
  UINT4 num_bounds;

  /* Bounds on the parameter space, in real parameter coordinates */
  gsl_matrix *param_bound_normal;
  gsl_matrix *param_bound_origin;

  /* Bounds on the parameter space, in lattice coordinates */
  gsl_matrix *latt_bound_normal;
  gsl_vector *latt_bound_dot;

  /* Lower and upper bounds on the parameter space, in real parameter coordinates */
  gsl_vector *param_lower;
  gsl_vector *param_upper;

  /* Conversion matrices from lattice to parameter space */
  gsl_matrix *latt_to_param;

  /* Lower and upper bounds on the parameter space, in lattice coordinates */
  gsl_vector_long *latt_lower;
  gsl_vector_long *latt_upper;

  /* Index of the dimension along which intersecctions with the parameter space are calculated */
  int line_index;

  /* Order that the parameter space dimensions are iterated over */
  gsl_vector_int *iter_order;

  /* Current lattice point */
  gsl_vector_long *latt_current;
  gsl_vector *current;

  /* Lower and upper bounds on the line_index'th dimension of parameter space */
  long line_latt_lower;
  long line_latt_upper;

  /* Start and direction of the parameter space line and its intersecting point */
  gsl_vector *line_start;
  gsl_vector *line_dir;

  /* Dot products of normal vectors with start and direction of line */
  gsl_vector *line_start_dot;
  gsl_vector *line_dir_dot;

  /* Intersection of normal vectors with line */
  gsl_vector *line_intersect;

  /* Additional function for checking if point is in parameter space */
  BOOLEAN (*in_param_space)(FlatLatticeTiling*);
  double *in_param_space_args;

  /* Number of templates generated */
  UINT8 template_count;

};

/******** Declarations ********/

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID(FLATLATTICETILINGH, "$Id$");

  FlatLatticeTiling *XLALCreateFlatLatticeTiling(UINT4, REAL8);

  void XLALDestroyFlatLatticeTiling(FlatLatticeTiling*);

  int XLALAddParameterSpaceBound(FlatLatticeTiling*, UINT4, UINT4, gsl_vector*, gsl_vector*);

  int XLALSquareParameterSpace(FlatLatticeTiling*, gsl_vector*, gsl_vector*);

  int XLALCubicLatticeGenerator(gsl_matrix*, REAL8*);

  int XLALAnstarLatticeGenerator(gsl_matrix*, REAL8*);

  int XLALSetupFlatLatticeTiling(FlatLatticeTiling*);

  int XLALNextFlatLatticePoint(FlatLatticeTiling*);

#ifdef __cplusplus
}
#endif

#endif
