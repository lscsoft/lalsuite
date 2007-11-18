/*
 *  Copyright (C) 2007 Karl Wette
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
 * \author K. Wette
 * \file
 * \brief Flat lattice tiling over multi-dimensioned parameter spaces
 */

#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_math.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/PtoleMetric.h>
#include <lal/FlatLatticeTiling.h>
#include <lalapps.h>

RCSID("$Id$");

int XLALWriteFlatLatticeTilingXMLFile(FlatLatticeTiling*, CHAR *);

int main(int argc, char *argv[]) {

  INT4 i, j, ii, jj;
  INT4 dims = 0;
  LALStatus status = blank_status;
  BOOLEAN is_help = 0;
  LALStringVector *list_bounds = NULL;
  gsl_vector *bounds = NULL;
  double bound = 0.0;
  gsl_vector *start = NULL;
  gsl_vector *width = NULL;
  gsl_matrix *metric = NULL;
  INT4 lattice_type = 0;
  INT4 metric_type = 0;
  REAL8 mismatch = 0.25;
  REAL8 Tspan = 1036800.0;
  CHAR *output_filename = NULL;
  FlatLatticeTiling *tiling = NULL;
  FILE *output_file = NULL;

  /* Initialise LAL error handler, debug level and log level */
  lal_errhandler = LAL_ERR_EXIT;
  LAL_CALL(LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LogSetLevel(lalDebugLevel);

  /* Register command line flags */
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",            'h', UVAR_HELP,     "Print this message", &is_help), &status);
  LAL_CALL(LALRegisterLISTUserVar  (&status, "bounds",          'b', UVAR_REQUIRED, "Start,width,start,... bounds on the parameter space", &list_bounds), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "metric_type",     'm', UVAR_OPTIONAL, "Type of metric (0=spindown, 1=eye)", &metric_type), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "lattice_type",    'L', UVAR_OPTIONAL, "Type of tiling lattice (0=Anstar, 1=cubic)", &lattice_type), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "metric_mismatch", 'u', UVAR_OPTIONAL, "Maximum mismatch of the templates", &mismatch), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "time_span",       'T', UVAR_OPTIONAL, "Time span of the data set (in seconds)", &Tspan), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "output_file",     'o', UVAR_OPTIONAL, "XML output file containing templates and metadata", &output_filename), &status);
  
  /* Read in command line */
  LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);
  if (is_help) {
    return EXIT_SUCCESS;
  }

  /* Parse the parameter space bounds */
  if (list_bounds->length/2 != (list_bounds->length+1)/2) {
    LALPrintError("%s\nERROR: Must be an even number of bounds\n", rcsid);
    return EXIT_FAILURE;
  }
  bounds = gsl_vector_alloc(list_bounds->length);
  for (i = 0; i < (INT4)bounds->size; ++i) {
    if (sscanf(list_bounds->data[i], "%le", &bound) != 1) {
      LALPrintError("%s\nERROR: Bound '%s' must be numberic\n", list_bounds->data[i], rcsid);
      return EXIT_FAILURE;
    }
    gsl_vector_set(bounds, i, bound);
  }

  /* Count the number of dimensions */
  for (i = 1; i < (INT4)bounds->size; i += 2) {
    if (gsl_vector_get(bounds, i) != 0.0) {
      ++dims;
    }
  }
  
  /* Create tiling structure */
  if ((tiling = XLALCreateFlatLatticeTiling(dims)) == NULL) {
    LALPrintError("%s\nERROR: XLALCreateFlatLatticeTiling failed\n", rcsid);
    return EXIT_FAILURE;
  }
  
  /* Setup parameter space */
  start = gsl_vector_alloc(dims);
  width = gsl_vector_alloc(dims);
  for (i = ii = 0; i < dims; ++i, ++ii) {
    
    if (gsl_vector_get(bounds, 2*ii+1) == 0.0) ++ii;

    gsl_vector_set(start, i, gsl_vector_get(bounds, 2*ii  ));
    gsl_vector_set(width, i, gsl_vector_get(bounds, 2*ii+1));

  }
  if (XLALSquareParameterSpace(tiling, start, width) != XLAL_SUCCESS) {
    LALPrintError("%s\nERROR: XLALSquareParameterSpace failed\n", rcsid);
    return EXIT_FAILURE;
  }

  /* Set metric */
  if (metric_type == 0) {

    if ((metric = XLALSpindownMetric(bounds->size/2, Tspan)) == NULL) {
      LALPrintError("%s\nERROR: XLALSpindownOnlyMetric failed\n", rcsid);
      return EXIT_FAILURE;
    }
 
    tiling->metric = gsl_matrix_alloc(dims, dims);
    for (i = ii = 0; i < dims; ++i, ++ii) {
      
      if (gsl_vector_get(bounds, 2*ii+1) == 0.0) ++ii;
      
      for (j = jj = 0; j < dims; ++j, ++jj) {
	
	if (gsl_vector_get(bounds, 2*jj+1) == 0.0) ++jj;

	gsl_matrix_set(tiling->metric, i, j, gsl_matrix_get(metric, ii, jj));
	
      }
      
    }
    
  }
  else if (metric_type == 1) {
    tiling->metric = gsl_matrix_alloc(tiling->dimension, tiling->dimension);
    gsl_matrix_set_identity(tiling->metric);
  }
  else {
    LALPrintError("%s\nERROR: Invalid metric type\n", rcsid);
    return EXIT_FAILURE;
  }

  /* Set mismatch */
  tiling->mismatch = mismatch;

  /* Set lattice generator */
  if (lattice_type == 0) {
    if ((tiling->generator = XLALAnstarLatticeGenerator(tiling->dimension)) == NULL) {
      LALPrintError("%s\nERROR: XLALAnstarLatticeGenerator failed\n", rcsid);
      return EXIT_FAILURE;
    }
  }
  else if (lattice_type == 1) {
    if ((tiling->generator = XLALCubicLatticeGenerator(tiling->dimension)) == NULL) {
      LALPrintError("%s\nERROR: XLALCubicLatticeGenerator failed\n", rcsid);
      return EXIT_FAILURE;
    }
  }
  else {
    LALPrintError("%s\nERROR: Invalid lattice type\n", rcsid);
    return EXIT_FAILURE;
  }

  /* Setup tiling */
  if (XLALSetupFlatLatticeTiling(tiling) != XLAL_SUCCESS) {
    LALPrintError("%s\nERROR: XLALSetupLatticeTiling failed\n", rcsid);
    return EXIT_FAILURE;
  }
  
  /* Write out XML file if requested */
  if (output_filename != NULL) {
    
    /* Open output file */
    if ((output_file = fopen(output_filename, "wt")) == NULL) {
      LALPrintError("%s\nERROR: Could not open output file '%s'\n", rcsid, output_filename);
      XLAL_ERROR("XLALWriteFlatLatticeTilingXMLFile", XLAL_EIO);
    }
    
    /* Write XML header */
    fprintf(output_file, "<?xml version=\"1.0\"?>\n");
    fprintf(output_file, "<flatlatticetiling version=\"%s\">\n", rcsid);
    
    /* Write dimension */
    fprintf(output_file, "  <dimension>%i</dimension>\n", tiling->dimension);

    /* Write indices of non-zero dimensions */
    fprintf(output_file, "  <indices>");
    for (i = 0; i < (INT4)bounds->size/2; ++i) {
      if (gsl_vector_get(bounds, 2*i+1) != 0) {
	fprintf(output_file, "%i ", i);
      }
    }
    fprintf(output_file, "</indices>\n");
    
    /* Write metric */
    fprintf(output_file, "  <metric>\n");
    for (i = 0; i < dims; ++i) {
      fprintf(output_file, "    ");
      for (j = 0; j < dims; ++j) {
	fprintf(output_file, "% 0.16e ", gsl_matrix_get(tiling->metric, i, j));
      }
      fprintf(output_file, ";\n");
    }
    fprintf(output_file, "  </metric>\n");
    
    /* Write mismatch */
    fprintf(output_file, "  <mismatch>%0.16e</mismatch>\n", tiling->mismatch);
    
    /* Write generator */
    fprintf(output_file, "  <generator>\n");
    for (i = 0; i < dims; ++i) {
      fprintf(output_file, "    ");
      for (j = 0; j < dims; ++j) {
	fprintf(output_file, "% 0.16e ", gsl_matrix_get(tiling->generator, i, j));
      }
      fprintf(output_file, ";\n");
    }
    fprintf(output_file, "  </generator>\n");
    
    /* Write XML description of parameter space bounds */
    if (tiling->bounds_xml_desc) {
      fprintf(output_file, "  <bounds>\n    %s\n  </bounds>\n", tiling->bounds_xml_desc);
    }
    else {
      fprintf(output_file, "  <bounds><type>unknown</type></bounds>\n");
    }
    
    /* Write increment vectors */
    fprintf(output_file, "  <increment>\n");
    for (i = 0; i < dims; ++i) {
      fprintf(output_file, "    ");
      for (j = 0; j < dims; ++j) {
	fprintf(output_file, "% 0.16e ", gsl_matrix_get(tiling->increment, i, j));
      }
      fprintf(output_file, ";\n");
    }
    fprintf(output_file, "  </increment>\n");
    
    /* Write padding */
    fprintf(output_file, "  <padding>\n    ");
    for (i = 0; i < dims; ++i) {
      fprintf(output_file, "% 0.16e ", gsl_vector_get(tiling->padding, i));
    }
    fprintf(output_file, "\n  </padding>\n");
    
    /* Generate templates */
    fprintf(output_file, "  <templates>\n");
    while (XLALNextFlatLatticePoint(tiling)) {
      
      /* Write template */
      fprintf(output_file, "    <row>");
      for (j = 0; j < dims; ++j) {
	fprintf(output_file, "% 0.16e ", gsl_vector_get(tiling->current, j));
      }
      fprintf(output_file, "</row>\n");
      
    }
    fprintf(output_file, "  </templates>\n");
    
    /* Write XML footer */
    fprintf(output_file, "</flatlatticetiling>\n");
    
    /* Close file */
    fclose(output_file);
    
  }				      
  
  /* Count the number of templates */
  printf("Number of templates generated: %lli\n", XLALTotalNumberOfFlatLatticePoints(tiling));
  
  /* Cleanup */
  LAL_CALL(LALDestroyUserVars(&status), &status);
  gsl_vector_free(bounds);
  gsl_vector_free(start);
  gsl_vector_free(width);
  gsl_matrix_free(metric);
  XLALDestroyFlatLatticeTiling(tiling);
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
