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
  REAL8 max_mismatch = 0.25;
  INT4 lattice_type = 0;
  BOOLEAN only_count = FALSE;
  REAL8 Tspan = 12.0 * 86400.0;
  CHAR *output_file = NULL;

  int i;
  VeryBasicXMLOutput xml = empty_VeryBasicXMLOutput;
  int spaces = 0;
  gsl_vector *bounds = NULL;
  FlatLatticeTiling *tiling = NULL;
  
  /* Initialise LAL error handler, debug level and log level */
  lal_errhandler = LAL_ERR_EXIT;
  LAL_CALL(LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  
  /* Register command line arguments */
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",           'h', UVAR_HELP, "Print this help message", &help), &status);
  LAL_CALL(LALRegisterLISTUserVar  (&status, "square",          0 , UVAR_OPTIONAL, "Square parameter space: start,width,...", &square), &status);
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
  XLAL_VBXMLO_BeginTag(&xml, "testFlatLatticeTilingPulsar");
  
  /* Create flat lattice tiling for square parameter space*/
  if (square && ++spaces == 1) {
    
    /* Get bounds */
    if ((bounds = XLALGSLVectorFromLALStringVector(square)) == NULL) {
      LALPrintError("XLALGSLVectorFromLALStringVector failed\n");
      return EXIT_FAILURE;
    }

    /* Create flat lattice tiling */
    if (XLAL_SUCCESS != XLALFlatLatticeTilingSquareParamSpace(&tiling, bounds)) {
      LALPrintError("XLALFlatLatticeTilingSquareParamSpace failed\n");
      return EXIT_FAILURE;
    }

  }

  /* Check only one parameter space was specified */
  if (spaces != 1) {
    LALPrintError("Exactly one of --square must be specified\n");
    return EXIT_FAILURE;
  }

  /* Output tiling dimension */
  XLAL_VBXMLO_Tag(&xml, "dimension", "%i", tiling->dimensions);
  
  /* Set metric */
  if (XLAL_SUCCESS != XLALSetFlatLatticeTilingSpindownFstatMetric(tiling, max_mismatch, Tspan)) {
    LALPrintError("XLALSetFlatLatticeTilingSpindownFstatMetric failed\n");
    return EXIT_FAILURE;
  }

  /* Output metric */
  XLAL_VBXMLO_gsl_matrix(&xml, "norm_metric", "%g", tiling->norm_metric);
  XLAL_VBXMLO_gsl_vector(&xml, "norm_to_real", "%g", tiling->norm_to_real);
  
  /* Set square parameter space bounds */
  if (square) {
    if (XLAL_SUCCESS != XLALFlatLatticeTilingSquareParamSpace(&tiling, bounds)) {
      LALPrintError("XLALFlatLatticeTilingSquareParamSpace failed\n");
      return EXIT_FAILURE;
    }
  }

  /* Output parameter space bounds */
  XLAL_VBXMLO_BeginTag(&xml, "bounds");
  for (i = 0; i < tiling->num_bounds; ++i) {
    const FlatLatticeTilingBound *this = &tiling->bounds[i];
    XLAL_VBXMLO_BeginTag(&xml, "bound");
    XLAL_VBXMLO_Tag(&xml, "dimension", "%i", this->dimension);
    switch (this->type) {
    case FLT_BT_Singular:

      /* Singular bound */
      XLAL_VBXMLO_BeginTag(&xml, "singular");
      XLAL_VBXMLO_Tag(&xml, "value", "%g", this->singular_value);
      XLAL_VBXMLO_EndTag(&xml, "singular");      
      break;
      
    case FLT_BT_Polynomial:

      /* Polynomial bound */
      XLAL_VBXMLO_BeginTag(&xml, "polynomial");
      XLAL_VBXMLO_BeginTag(&xml, "lower");
      XLAL_VBXMLO_gsl_vector(&xml, "const", "%g", this->poly_lower_const);
      if (this->poly_lower_exp)
	XLAL_VBXMLO_gsl_matrix_int(&xml, "exp", "%i", this->poly_lower_exp);
      XLAL_VBXMLO_EndTag(&xml, "lower");
      XLAL_VBXMLO_BeginTag(&xml, "upper");
      XLAL_VBXMLO_gsl_vector(&xml, "const", "%g", this->poly_upper_const);
      if (this->poly_upper_exp)
	XLAL_VBXMLO_gsl_matrix_int(&xml, "exp", "%i", this->poly_upper_exp);
      XLAL_VBXMLO_EndTag(&xml, "upper");
      XLAL_VBXMLO_EndTag(&xml, "polynomial");
      break;

    default:
      break;
    }
    XLAL_VBXMLO_EndTag(&xml, "bound");
  }
  XLAL_VBXMLO_EndTag(&xml, "bounds");

  /* Output dimension maps */
  XLAL_VBXMLO_Tag(&xml, "reduced_dims", "%i", tiling->reduced_dims);
  XLAL_VBXMLO_gsl_vector_int(&xml, "reduced_map", "%i", tiling->reduced_map);
  
  /* Set lattice */
  switch (lattice_type) {
  case 0:
    if (XLAL_SUCCESS != XLALSetFlatTilingAnstarLattice(tiling)) {
      LALPrintError("XLALSetFlatTilingAnstarLattice failed\n");
      return EXIT_FAILURE;
    }
    break;
  case 1:
    if (XLAL_SUCCESS != XLALSetFlatTilingCubicLattice(tiling)) {
      LALPrintError("XLALSetFlatTilingCubicLattice failed\n");
      return EXIT_FAILURE;
    }
    break;
  }

  /* Output lattice transformations */
  XLAL_VBXMLO_gsl_matrix(&xml, "latt_to_norm", "%e", tiling->latt_to_norm);

  /* Generate templates */
  if (!only_count) {
    XLAL_VBXMLO_BeginTag(&xml, "templates");
    while (XLALNextFlatLatticePoint(tiling) == XLAL_SUCCESS) {
      XLAL_VBXMLO_gsl_vector(&xml, "row", "%g", tiling->current);
      fflush(xml.file);
    }
    XLAL_VBXMLO_EndTag(&xml, "templates");
  }
  XLAL_VBXMLO_Tag(&xml, "template_count", "%0.0f", XLALTotalFlatLatticePointCount(tiling));
  
  /* Close XML output file */
  XLAL_VBXMLO_EndTag(&xml, "testFlatLatticeTilingPulsar");
  if (output_file)
    fclose(xml.file);

  /* Cleanup */
  LALDestroyUserVars(&status);
  if (bounds)
    gsl_vector_free(bounds);
  XLALFreeFlatLatticeTiling(tiling);
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
