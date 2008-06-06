/*
 *  Copyright (C) 2007, 2008 Karl Wette
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
 * \author Karl Wette
 * \file
 * \brief Test program for FlatLatticeTiling*
 */

#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_math.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/UserInput.h>
#include <lal/FlatLatticeTiling.h>
#include <lal/FlatLatticeTilingPulsar.h>
#include <lal/FlatLatticeTilingSupport.h>
#include <lal/VeryBasicXMLOutput.h>
#include <lalapps.h>

RCSID("$Id$");

static const BOOLEAN TRUE  = (1==1);
static const BOOLEAN FALSE = (0==1);

int main(int argc, char *argv[]) {
  
  LALStatus status = blank_status;
  BOOLEAN help = FALSE;
  LALStringVector *square = NULL;
  BOOLEAN no_check_linear = FALSE;
  INT4 metric_type = 0;
  REAL8 max_mismatch = 0.25;
  INT4 lattice_type = 0;
  BOOLEAN only_count = FALSE;
  REAL8 Tspan = 12.0 * 86400.0;
  CHAR *output_file = NULL;
  VeryBasicXMLOutput xml = empty_VeryBasicXMLOutput;
  gsl_vector *bounds = NULL;
  INT4 dimension = 0;
  FlatLatticeTiling *tiling = NULL;
  
  /* Initialise LAL error handler, debug level and log level */
  lal_errhandler = LAL_ERR_EXIT;
  LAL_CALL(LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  
  /* Register command line arguments */
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",           'h', UVAR_HELP, "Print this help message", &help), &status);
  LAL_CALL(LALRegisterLISTUserVar  (&status, "square",          0 , UVAR_OPTIONAL, "Square parameter space: start,width,...", &square), &status);
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "no-check-linear", 0 , UVAR_OPTIONAL, "Don't do careful checking of linear bounds", &no_check_linear), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "metric",         'M', UVAR_OPTIONAL, "Metric: 0=spindown, 1=eye", &metric_type), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "max-mismatch",   'X', UVAR_OPTIONAL, "Maximum allowed mismatch between the templates", &max_mismatch), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "lattice",        'L', UVAR_OPTIONAL, "Lattice: 0=Anstar, 1=cubic", &lattice_type), &status);
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "only-count",     'c', UVAR_OPTIONAL, "Only count number of templates", &only_count), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "time-span",      'T', UVAR_OPTIONAL, "Time-span of the data set (in seconds)", &Tspan), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "output",         'o', UVAR_OPTIONAL, "Output file", &output_file), &status);
  
  /* Get command line arguments */
  LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);
  if (help)
    return EXIT_SUCCESS;
  
  /* Open XML output file */
  if (output_file)
    if ((xml.file = fopen(output_file, "w")) == NULL) {
      LALPrintError("Could not open file '%s'\n", output_file);
      return EXIT_FAILURE;
    }
  XLAL_VBXMLO_Header(&xml, 1, 0);
  XLAL_VBXMLO_BeginTag(&xml, "testFlatLatticeTilingPulsar", NULL);
  
  /* Get bounds, and deduce dimension */
  {
    int spaces = 0;

    /* Square parameter space */
    if (square && !spaces) {
      ++spaces;

      /* Get dimension */
      if (GSL_IS_ODD(square->length)) {
	LALPrintError("--square must have an even number of arguments\n");
	return EXIT_FAILURE;
      }
      dimension = square->length/2;
      
      /* Get bounds */
      if ((bounds = XLALGSLVectorFromLALStringVector(square)) == NULL) {
	LALPrintError("XLALGSLVectorFromLALStringVector failed\n");
	return EXIT_FAILURE;
      }

    }

    /* Check only one parameter space was specified */
    if (spaces != 1) {
      LALPrintError("Exactly one of --square must be specified\n");
      return EXIT_FAILURE;
    }
  }
  
  /* Create flat lattice tiling */
  if ((tiling = XLALCreateFlatLatticeTiling(dimension)) == NULL) {
    LALPrintError("XLALCreateFlatLatticeTiling failed\n");
    return EXIT_FAILURE;
  }

  /* Create parameter spaces */
  if (square) {

    gsl_vector *linear_bound = NULL;
    
    /* Allocate memory */
    if ((linear_bound = gsl_vector_alloc(dimension)) == NULL) {
      LALPrintError("Failed to allocate a gsl_vector\n");
      return EXIT_FAILURE;
    }
    {

      int i;
    
      gsl_vector_view start = gsl_vector_subvector_with_stride(bounds, 0, 2, dimension);
      gsl_vector_view width = gsl_vector_subvector_with_stride(bounds, 1, 2, dimension);
      
      /* Set square bounds on each dimension */
      for (i = 0; i < dimension; ++i) {
	
	/* If width is zero */
	if (gsl_vector_get(&width.vector, i) == 0.0) {

	  /* Set singular bound */
	  if (XLAL_SUCCESS != XLALSetFlatLatticeTilingSingularBound(tiling, i, gsl_vector_get(&start.vector, i))) {
	    LALPrintError("XLALSetFlatLatticeTilingSingularBound failed\n");
	    return EXIT_FAILURE;
	  }

	}
	else {
	  
	  /* Set lower bound */
	  gsl_vector_set_zero(linear_bound);
	  gsl_vector_set(linear_bound, i, -1.0);
	  if (XLAL_SUCCESS != XLALAddFlatLatticeTilingLinearBound(tiling, i, linear_bound, -1.0 * gsl_vector_get(&start.vector, i), !no_check_linear)) {
	    LALPrintError("XLALAddFlatLatticeTilingLinearBound failed\n");
	    return EXIT_FAILURE;
	  }

	  /* Set upper bound */
	  gsl_vector_set_zero(linear_bound);
	  gsl_vector_set(linear_bound, i, 1.0);
	  if (XLAL_SUCCESS != XLALAddFlatLatticeTilingLinearBound(tiling, i, linear_bound, 
								  gsl_vector_get(&start.vector, i) + gsl_vector_get(&width.vector, i), !no_check_linear)) {
	    LALPrintError("XLALAddFlatLatticeTilingLinearBound failed\n");
	    return EXIT_FAILURE;
	  }
	  
	}

      }

    }
    
    /* Cleanup */
    gsl_vector_free(linear_bound);

  }
  
  /* Set metric */
  switch (metric_type) {
  case 0:
    if (XLAL_SUCCESS != XLALSetFlatLatticeTilingSpindownFstatMetric(tiling, max_mismatch, Tspan)) {
      LALPrintError("XLALSetFlatLatticeTilingSpindownFstatMetric failed\n");
      return EXIT_FAILURE;
    }
    break;
  case 1:
    {
      gsl_matrix *metric = NULL;
      if ((metric = gsl_matrix_alloc(dimension, dimension)) == NULL) {
	LALPrintError("Failed to allocate a gsl_matrix\n");
	return EXIT_FAILURE;
      }
      gsl_matrix_set_identity(metric);
      if (XLAL_SUCCESS != XLALSetFlatLatticeTilingMetric(tiling, metric, max_mismatch, NULL)) {
 	LALPrintError("XLALSetFlatLatticeTilingMetric failed\n");
	return EXIT_FAILURE;
      }
      gsl_matrix_free(metric);
    }
    break;
  default:
    LALPrintError("Invalid value for --metric\n");
    return EXIT_FAILURE;
  }
  XLAL_VBXMLO_gsl_matrix(&xml, tiling->norm_metric, "norm_metric", NULL);
  XLAL_VBXMLO_Tag(&xml, "max_mismatch", NULL, "%g", tiling->max_mismatch);
  XLAL_VBXMLO_gsl_vector(&xml, tiling->norm_to_real_mul, "norm_to_real_mul", NULL);
  XLAL_VBXMLO_gsl_vector(&xml, tiling->norm_to_real_add, "norm_to_real_add", NULL);
  XLAL_VBXMLO_gsl_matrix(&xml, tiling->linear_bound_A, "linear_bound_A", NULL);
  XLAL_VBXMLO_gsl_vector(&xml, tiling->linear_bound_b, "linear_bound_b", NULL);
  XLAL_VBXMLO_gsl_matrix(&xml, tiling->linear_vertices, "linear_vertices", NULL);

  /* Set lattice */
  switch (lattice_type) {
  case 0:
    if (XLAL_SUCCESS != XLALSetAnstarTilingLattice(tiling)) {
      LALPrintError("XLALSetAnstarTilingLattice failed\n");
      return EXIT_FAILURE;
    }
    break;
  case 1:
    if (XLAL_SUCCESS != XLALSetCubicTilingLattice(tiling)) {
      LALPrintError("XLALSetCubicTilingLattice failed\n");
      return EXIT_FAILURE;
    }
    break;
  }

  /* Generate templates */
  if (!only_count) {
    XLAL_VBXMLO_BeginTag(&xml, "templates", NULL);
    while (XLALNextFlatLatticeTile(tiling) == XLAL_SUCCESS) {
      XLAL_VBXMLO_gsl_vector(&xml, tiling->current_tile, "row", NULL);
      fflush(xml.file);
    }
    XLAL_VBXMLO_EndTag(&xml, "templates");
  }
  XLAL_VBXMLO_Tag(&xml, "template_count", NULL, "%0.0f", XLALTotalFlatLatticeTileCount(tiling));
  
  /* Close XML output file */
  XLAL_VBXMLO_EndTag(&xml, "testFlatLatticeTilingPulsar");
  if (output_file)
    fclose(xml.file);

  /* Cleanup */
  LALDestroyUserVars(&status);
  gsl_vector_free(bounds);
  XLALFreeFlatLatticeTiling(tiling);

  return XLAL_SUCCESS;

}
