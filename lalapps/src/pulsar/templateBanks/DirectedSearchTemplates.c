/**
 * \author K. Wette
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
#include <lal/PulsarDataTypes.h>

#include "FlatLatticeTiling.h"
#include "DirectedSearch.h"

#include <lalapps.h>

static const BOOLEAN TRUE  = (1==1);
static const BOOLEAN FALSE = (1==0);

RCSID("$Id$");

int main(int argc, char *argv[]) {

  INT4 i;

  /* LAL status */
  LALStatus status = blank_status;

  /* Print help message? */
  BOOLEAN is_help = FALSE;

  /* Right ascension and declenation */
  REAL8 alpha = 6.12;
  REAL8 delta = 1.02;

  /* Frequency and spindowns */
  REAL8 freq       = 100.0;
  REAL8 freq_band  = 1.0e-2;
  REAL8 f1dot      = 0.0;
  REAL8 f1dot_band = 0.0;
  REAL8 f2dot      = 0.0;
  REAL8 f2dot_band = 0.0;
  REAL8 f3dot      = 0.0;
  REAL8 f3dot_band = 0.0;
  
  /* Type of lattice */
  INT4 lattice_type = 0;

  /* Type of metric */
  INT4 metric_type = 0;

  /* Maximum template mismatch */
  REAL8 mismatch = 0.25;

  /* Time span of the SFTs */
  REAL8 Tspan = 1036800.0;

  /* Output files */
  CHAR *output_filename = NULL;
  FILE *output_file = NULL;

  /* Dimensionality of parameter space */
  INT4 dimension = 0;

  /* Flat lattice tiling structure */
  FlatLatticeTiling *tiling = NULL;

  /* Initialise LAL error handler, debug level and log level */
  lal_errhandler = LAL_ERR_EXIT;
  LAL_CALL(LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LogSetLevel(lalDebugLevel);

  /* Register command line flags */
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",            'h', UVAR_HELP,      "Print this message", &is_help), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "Alpha",           'a', UVAR_OPTIONAL,  "Right ascension of the target object (in radians)", &alpha), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "Delta",           'd', UVAR_OPTIONAL,  "Declination of the target object (in radians)", &delta), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "Freq",            'f', UVAR_OPTIONAL,  "Starting frequency of search band (in Hertz)", &freq), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "FreqBand",        'b', UVAR_OPTIONAL,  "Width of frequency search band (in Hertz)", &freq_band), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f1dot",            0 , UVAR_OPTIONAL,  "First spindown (in Hertz/s)", &f1dot), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f1dotBand",        0 , UVAR_OPTIONAL,  "First spindown band (in Hertz/s)", &f1dot_band), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f2dot",            0 , UVAR_OPTIONAL,  "Second spindown (in Hertz/s)", &f2dot), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f2dotBand",        0 , UVAR_OPTIONAL,  "Second spindown band (in Hertz/s)", &f2dot_band), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f3dot",            0 , UVAR_OPTIONAL,  "Third spindown (in Hertz/s)", &f3dot), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f3dotBand",        0 , UVAR_OPTIONAL,  "Third spindown band (in Hertz/s)", &f3dot_band), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "latticeType",     'L', UVAR_OPTIONAL,  "Type of tiling lattice (0=Anstar, 1=cubic)", &lattice_type), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "metricType",      'M', UVAR_OPTIONAL,  "Type of metric (0=spindown, 1=identity, for testing)", &metric_type), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "metricMismatch",  'X', UVAR_OPTIONAL,  "Maximum mismatch of the search templates", &mismatch), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "timeSpan",        'T', UVAR_OPTIONAL,  "Upper bound on time span of the SFTs (in seconds)", &Tspan), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "outputFile",      'o', UVAR_OPTIONAL,  "Output file containing the generated templates", &output_filename), &status);

  /* Read in command line */
  LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);
  if (is_help) {
    return EXIT_SUCCESS;
  }

  /* Create tiling structure */
  if (LALUserVarWasSet(&freq_band)) {
    ++dimension;
    if (LALUserVarWasSet(&f1dot_band)) {
      ++dimension;
      if (LALUserVarWasSet(&f2dot_band)) {
	++dimension;
	if (LALUserVarWasSet(&f3dot_band)) {
	  ++dimension;
	}
      }
    }
  }
  if ((tiling = XLALCreateFlatLatticeTiling(dimension)) == NULL) {
    LALPrintError("%s\nERROR: XLALCreateFlatLatticeTiling failed\n", rcsid);
    return EXIT_FAILURE;
  }
  
  /* Fill spindown metric and mismatch*/
  if (metric_type == 0) {
    if ((tiling->metric = XLALSpindownMetric(tiling->dimension, Tspan)) == NULL) {
      LALPrintError("%s\nERROR: XLALSpindownOnlyMetric failed\n", rcsid);
      return EXIT_FAILURE;
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
  tiling->mismatch = mismatch;

  /* Fill lattice generator */
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

  /* Setup parameter space */
  if (XLALSquareParameterSpace(tiling, freq, freq + freq_band, f1dot, f1dot + f1dot_band, f2dot, f2dot + f2dot_band, f3dot, f3dot + f3dot_band) != XLAL_SUCCESS) {
    LALPrintError("%s\nERROR: XLALSquareParameterSpace failed\n", rcsid);
    return EXIT_FAILURE;
  }
    
  /* Setup tiling */
  if (XLALSetupFlatLatticeTiling(tiling) != XLAL_SUCCESS) {
    LALPrintError("%s\nERROR: XLALSetupLatticeTiling failed\n", rcsid);
    return EXIT_FAILURE;
  }

  /* Open output file */
  if (output_filename != NULL) {
    if ((output_file = fopen(output_filename, "wt")) == NULL) {
      LALPrintError("%s\nERROR: Could not open output file '%s'\n", rcsid, output_filename);
      return EXIT_FAILURE;
    }
    if (tiling->dimension > 0) {
      fprintf(output_file, "% 0.16e", freq);
      if (tiling->dimension > 1) {
	fprintf(output_file, "% 0.16e", f1dot); 
	if (tiling->dimension > 2) {
	  fprintf(output_file, "% 0.16e", f2dot);
	  if (tiling->dimension > 3) {
	    fprintf(output_file, "% 0.16e", f3dot);
	  }
	}
      }
      fprintf(output_file, "\n");
    }
    if (tiling->dimension > 0) {
      fprintf(output_file, "% 0.16e", freq_band);
      if (tiling->dimension > 1) {
	fprintf(output_file, "% 0.16e", f1dot_band); 
	if (tiling->dimension > 2) {
	  fprintf(output_file, "% 0.16e", f2dot_band);
	  if (tiling->dimension > 3) {
	    fprintf(output_file, "% 0.16e", f3dot_band);
	  }
	}
      }
      fprintf(output_file, "\n");
    }
  }
  
  /* Step through templates */
  while (XLALNextFlatLatticePoint(tiling) == XLAL_SUCCESS) {

    /* If required write template to output file */
    if (output_file != NULL) {

      /* Print frequency, sky position and at least PULSAR_MAX_SPINS spindowns as required by CFS_v2 */
      for (i = 0; i < tiling->dimension; ++i) {
	fprintf(output_file, " % 0.16e", XLALCurrentFlatLatticePoint(tiling, i));
      }
      fprintf(output_file, "\n");
	
    }

  }
  printf("Number of templates generated: %lli\n", XLALNumberOfFlatLatticeTiles(tiling));

  /* Close output file */
  if (output_file != NULL) {
    fclose(output_file);
  }

  /* Cleanup */
  LAL_CALL(LALDestroyUserVars(&status), &status);
  XLALDestroyFlatLatticeTiling(tiling);
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
