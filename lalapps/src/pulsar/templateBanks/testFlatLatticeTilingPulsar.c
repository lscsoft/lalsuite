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
#include <lal/GSLSupport.h>
#include <lal/BitField.h>
#include <lal/DopplerFullScan.h>
#include <lal/VeryBasicXMLOutput.h>
#include <lal/FlatLatticeTiling.h>
#include <lal/FlatLatticeTilingPulsar.h>
#include <lalapps.h>

RCSID("$Id$");

#define TRUE  (1==1)
#define FALSE (1==0)

#define LALAPPS_ERROR(message, errno) \
{                                     \
  LALPrintError(message);             \
  return EXIT_FAILURE;                \
}

int main(int argc, char *argv[]) {
  
  LALStatus status = blank_status;
  BOOLEAN help = FALSE;
  REAL8 Tspan = 0.0;
  LALStringVector *square = NULL;
  LALStringVector *age_brake = NULL;
  REAL8 max_mismatch = 0.0;
  INT4 lattice_type = 0;
  INT4 metric_type = 0;
  BOOLEAN only_count = FALSE;
  REAL8 inject_ratio = 0.0;
  INT4 inject_min_mismatch_bins = 100;
  CHAR *output_file = NULL;
  REAL8 scale_padding = 0;

  int i, j;
  UINT4 k;
  VeryBasicXMLOutput xml = empty_VeryBasicXMLOutput;
  int spaces = 0;
  FlatLatticeTiling *tiling = NULL;
  gsl_vector *temp = NULL;
  UINT4 inject_count = 0;
  INT4 inject_seed = 0;
  gsl_vector *inject_point = NULL;
  gsl_vector *inject_min_mismatch = NULL;
  RandomParams *inject_random = NULL;
  REAL8 inject_dist = 0.0;
  gsl_vector_int *inject_min_mismatch_hist = NULL;
  
  /* Initialise LAL error handler, debug level and log level */
  lal_errhandler = LAL_ERR_EXIT;
  LAL_CALL(LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  
  /* Register command line arguments */
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",             'h', UVAR_HELP, "Print this help message", &help), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "time-span",        'T', UVAR_OPTIONAL, "Time-span of the data set (in seconds)", &Tspan), &status);
  LAL_CALL(LALRegisterLISTUserVar  (&status, "square",            0 , UVAR_OPTIONAL, "Square parameter space: start,width,...", &square), &status);
  LAL_CALL(LALRegisterLISTUserVar  (&status, "age-braking",       0 , UVAR_OPTIONAL, "Age/braking index parameter space: "
				                                                     "freq,freqband,alpha,delta,age,minbrake,maxbrake", &age_brake), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "max-mismatch",     'X', UVAR_OPTIONAL, "Maximum allowed mismatch between the templates", &max_mismatch), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "lattice",          'L', UVAR_OPTIONAL, "Lattice: 0=Anstar, 1=cubic", &lattice_type), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "metric",           'M', UVAR_OPTIONAL, "Metric: 0=spindown, 1=eye", &metric_type), &status);
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "only-count",       'c', UVAR_OPTIONAL, "Only count number of templates", &only_count), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "inject-ratio",     'R', UVAR_OPTIONAL, "Inject this ratio of points to templates", &inject_ratio), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "inject-bins",      'B', UVAR_OPTIONAL, "Use this number of bins for the mismatch histogram", &inject_min_mismatch_bins), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "output",           'o', UVAR_OPTIONAL, "Output file", &output_file), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "scale-padding",     0 , UVAR_DEVELOPER, "Scale the padding of the parameter space bounds (for testing)", &scale_padding), &status);
  
  /* Get command line arguments */
  LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);
  if (help)
    return EXIT_SUCCESS;
  
  /* Open XML output file */
  if (output_file)
    if ((xml.file = fopen(output_file, "w")) == NULL)
      LALAPPS_ERROR("Could not open output file\n", 0);
  XLAL_VBXMLO_Header(&xml, 1, 0);
  XLAL_VBXMLO_BeginTag(&xml, "testFlatLatticeTilingPulsar");

  /* Create square parameter space */
  if (LALUserVarWasSet(&square) && ++spaces == 1) {

    gsl_vector *bounds = NULL;
    
    /* Get bounds */
    if ((bounds = XLALGSLVectorFromLALStringVector(square)) == NULL)
      LALAPPS_ERROR("XLALGSLVectorFromLALStringVector failed\n", 0);
    if (GSL_IS_ODD(bounds->size))
      LALAPPS_ERROR("--square must have an even number of arguments", 0);
    
    /* Create flat lattice tiling */
    if (NULL == (tiling = XLALCreateFlatLatticeTiling(bounds->size/2)))
      LALAPPS_ERROR("XLALCreateFlatLatticeTiling failed\n", 0);
    
    /* Create square parameter space */
    for (i = 0; i < (int)(bounds->size/2); ++i) {
      if (XLAL_SUCCESS != XLALAddFlatLatticeTilingConstantBound(tiling, i, gsl_vector_get(bounds, 2*i),
								gsl_vector_get(bounds, 2*i) + gsl_vector_get(bounds, 2*i + 1)))
	LALAPPS_ERROR("XLALAddFlatLatticeTilingConstantBound failed\n", 0);
    }
      
    /* Output parameter space */
    XLAL_VBXMLO_gsl_vector(&xml, "square_param_space", "%0.12g", bounds);

    /* Cleanup */
    FREE_GSL_VECTOR(bounds);

  }

  /* Create age-braking index parameter space */
  if (LALUserVarWasSet(&age_brake) && ++spaces == 1) {

    gsl_vector *values = NULL;

    /* Get bounds */
    if ((values = XLALGSLVectorFromLALStringVector(age_brake)) == NULL)
      LALAPPS_ERROR("XLALGSLVectorFromLALStringVector failed\n", 0);
    if (values->size != 7)
      LALAPPS_ERROR("--age-braking must have exactly 7 arguments", 0);
    {
      const double freq        = gsl_vector_get(values, 0);
      const double freq_band   = gsl_vector_get(values, 1);
      const double alpha       = gsl_vector_get(values, 2);
      const double delta       = gsl_vector_get(values, 3);
      const double age         = gsl_vector_get(values, 4);
      const double min_braking = gsl_vector_get(values, 5);
      const double max_braking = gsl_vector_get(values, 6);
      
      /* Create flat lattice tiling */
      if (NULL == (tiling = XLALCreateFlatLatticeTiling(5)))
	LALAPPS_ERROR("XLALCreateFlatLatticeTiling failed\n", 0);

      /* Add sky position bounds */
      if (XLAL_SUCCESS != XLALAddFlatLatticeTilingConstantBound(tiling, 1, alpha, alpha))
	LALAPPS_ERROR("XLALAddFlatLatticeTilingConstantBound failed\n", 0);
      if (XLAL_SUCCESS != XLALAddFlatLatticeTilingConstantBound(tiling, 2, delta, delta))
	LALAPPS_ERROR("XLALAddFlatLatticeTilingConstantBound failed\n", 0);

      /* Add frequency and spindown bounds */
      if (XLAL_SUCCESS != XLALAddFlatLatticeTilingAgeBrakingIndexBounds(tiling, freq, freq_band, age, min_braking, max_braking, 0, 2))
	LALAPPS_ERROR("XLALAddFlatLatticeTilingAgeBrakingIndexBounds failed\n", 0);
    
    }

    /* Output parameter space */
    XLAL_VBXMLO_gsl_vector(&xml, "age_braking_param_space", "%0.12g", values);

    /* Cleanup */
    FREE_GSL_VECTOR(values);

  }

  /* Check only one parameter space was specified */
  if (spaces != 1)
    LALAPPS_ERROR("Exactly one of --square or --age-braking must be specified\n", 0);
  XLAL_VBXMLO_Tag(&xml, "dimensions", "%i", tiling->dimensions);

  /* Scale padding (testing only) */
  if (LALUserVarWasSet(&scale_padding)) {
    tiling->scale_padding = scale_padding;
  }

  /* Set metric */
  switch (metric_type) {
  case 0:
    if (XLAL_SUCCESS != XLALSetFlatLatticeTilingSpindownFstatMetric(tiling, max_mismatch, Tspan))
      LALAPPS_ERROR("XLALSetFlatLatticeTilingSpindownFstatMetric failed\n", 0);
    break;
  case 1:
    {
      gsl_matrix *identity = NULL;
      ALLOC_GSL_MATRIX(identity, tiling->dimensions, tiling->dimensions, LALAPPS_ERROR);
      gsl_matrix_set_identity(identity);
      if (XLAL_SUCCESS != XLALSetFlatLatticeTilingMetric(tiling, identity, max_mismatch, NULL))
	LALAPPS_ERROR("XLALSetFlatLatticeTilingMetric failed\n", 0);
      gsl_matrix_free(identity);
    }
    break;
  default:
    LALAPPS_ERROR("Invalid --metric\n", 0);
  }
  XLAL_VBXMLO_gsl_matrix(&xml, "metric", "%0.12g", tiling->metric);
  XLAL_VBXMLO_Tag(&xml, "max_mismatch", "%0.12g", tiling->max_mismatch);
  XLAL_VBXMLO_gsl_vector(&xml, "real_scale", "%0.12g", tiling->real_scale);
  XLAL_VBXMLO_gsl_vector(&xml, "real_offset", "%0.12g", tiling->real_offset);

  /* Set lattice */
  switch (lattice_type) {
  case 0:
    if (XLAL_SUCCESS != XLALSetFlatTilingAnstarLattice(tiling))
      LALAPPS_ERROR("XLALSetFlatTilingAnstarLattice failed\n", 0);
    break;
  case 1:
    if (XLAL_SUCCESS != XLALSetFlatTilingCubicLattice(tiling))
      LALAPPS_ERROR("XLALSetFlatTilingCubicLattice failed\n", 0);
    break;
  default:
    LALAPPS_ERROR("Invalid --lattice\n", 0);
  }

  /* Setup injections */
  if (inject_ratio > 0.0) {

    /* Calculate number of injections */
    if (0 == (inject_count = XLALTotalFlatLatticePointCount(tiling)))
      LALAPPS_ERROR("'tiling' did not generate any templates!\n", 0);
    inject_count = (UINT4)ceil(inject_ratio * inject_count);

    /* Create random number seed */
    while (inject_seed == 0)
      inject_seed = time(NULL);

    /* Create random number generator */
    if ((inject_random = XLALCreateRandomParams(inject_seed)) == NULL)
      LALAPPS_ERROR("XLALCreateRandomParams failed", 0);
    
    /* Allocate memory */
    ALLOC_GSL_VECTOR(inject_point, tiling->dimensions, LALAPPS_ERROR);
    ALLOC_GSL_VECTOR(inject_min_mismatch, inject_count, LALAPPS_ERROR);

    /* Initialise minimum mismatch */
    gsl_vector_set_all(inject_min_mismatch, GSL_POSINF);
    
  }

  /* Generate and output templates and injections */
  fflush(xml.file);
  if (!only_count) {
    ALLOC_GSL_VECTOR(temp, tiling->dimensions, LALAPPS_ERROR);
    XLAL_VBXMLO_BeginTag(&xml, "tiling");
  }
  while (XLALNextFlatLatticePoint(tiling) == XLAL_SUCCESS) {
    
    /* Output template */
    if (!only_count) {
      XLAL_VBXMLO_gsl_vector(&xml, "template", "%0.12g", tiling->current);
      fflush(xml.file);
    }

    /* Do injections */
    if (inject_count > 0) {

      /* Reset random number generator */
      XLALResetRandomParams(inject_random, inject_seed);

      /* Generate injections */
      for (k = 0; k < inject_count; ++k) {
	if (XLAL_SUCCESS != XLALRandomPointInFlatLatticeParamSpace(tiling, inject_random, inject_point, tiling->current, &inject_dist))
	  LALAPPS_ERROR("XLALRandomPointInFlatLatticeParamSpace failed\n", 0);
	if (!only_count && tiling->count == 1)
	  XLAL_VBXMLO_gsl_vector(&xml, "injection", "%0.12g", inject_point);	  

	/* Update minimum mismatch */
	if (gsl_vector_get(inject_min_mismatch, k) > inject_dist)
	  gsl_vector_set(inject_min_mismatch, k, inject_dist);

      }
      
    }
    
  }
  if (!only_count)
    XLAL_VBXMLO_EndTag(&xml, "tiling");

  /* Output subspaces */
  XLAL_VBXMLO_BeginTag(&xml, "subspaces");
  for (i = 0; i < tiling->num_subspaces; ++i) {
    XLAL_VBXMLO_BeginTag(&xml, "subspace");
    XLAL_VBXMLO_Tag(&xml, "dimensions", "%i", tiling->subspaces[i]->dimensions);
    XLAL_VBXMLO_BeginTag(&xml, "is_tiled");
    XLAL_VBXMLO_Indent(&xml);
    for (j = 0; j < tiling->dimensions; ++j) {
      XLAL_VBXMLO_Printf(&xml, "%c", GET_BIT(UINT8, tiling->subspaces[i]->is_tiled, j) ? 'Y' : 'N');
    }
    XLAL_VBXMLO_Printf(&xml, "\n");
    XLAL_VBXMLO_EndTag(&xml, "is_tiled");
    XLAL_VBXMLO_gsl_vector(&xml, "padding", "%0.12g", tiling->subspaces[i]->padding);
    XLAL_VBXMLO_gsl_matrix(&xml, "increment", "%0.12g", tiling->subspaces[i]->increment);
    XLAL_VBXMLO_EndTag(&xml, "subspace");
  }
  XLAL_VBXMLO_EndTag(&xml, "subspaces");

  /* Output template count */
  XLAL_VBXMLO_Tag(&xml, "template_count", "%li", XLALTotalFlatLatticePointCount(tiling));

  /* Output injection results */
  if (inject_count > 0) {

    UINT4 inject_unmatched = 0;

    /* Allocate memory */
    ALLOC_GSL_VECTOR_INT(inject_min_mismatch_hist, inject_min_mismatch_bins, LALAPPS_ERROR);

    /* Iterate over injections */
    gsl_vector_int_set_zero(inject_min_mismatch_hist);
    inject_unmatched = 0;
    for (k = 0; k < inject_count; ++k) {

      /* Minimal mistmatch histogram bin */
      const int bin = (int)floor(gsl_vector_get(inject_min_mismatch, k) / tiling->max_mismatch * inject_min_mismatch_bins);

      /* If within range, increase count */
      if (0 <= bin && bin < inject_min_mismatch_bins)
	gsl_vector_int_set(inject_min_mismatch_hist, bin, 1 + gsl_vector_int_get(inject_min_mismatch_hist, bin));

      /* Otherwise increase unmatched injection count */
      else
	++inject_unmatched;
      
    }
    XLAL_VBXMLO_Tag(&xml, "injection_count", "%li", inject_count);
    XLAL_VBXMLO_Tag(&xml, "unmatched_injections", "%li", inject_unmatched);
    XLAL_VBXMLO_gsl_vector_int(&xml, "injection_mismatch_histogram", "%i", inject_min_mismatch_hist);

  }      
  
  /* Close XML output file */
  XLAL_VBXMLO_EndTag(&xml, "testFlatLatticeTilingPulsar");
  if (output_file)
    fclose(xml.file);

  /* Cleanup */
  LALDestroyUserVars(&status);
  XLALFreeFlatLatticeTiling(tiling);
  XLALDestroyRandomParams(inject_random);
  FREE_GSL_VECTOR(temp);
  FREE_GSL_VECTOR(inject_point);
  FREE_GSL_VECTOR(inject_min_mismatch);
  FREE_GSL_VECTOR_INT(inject_min_mismatch_hist);
  LALCheckMemoryLeaks();
  
  return EXIT_SUCCESS;
  
}
