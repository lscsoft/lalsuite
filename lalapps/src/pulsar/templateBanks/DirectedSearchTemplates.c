/*********************************************************************************
 *  \author K. Wette
 *  \file
 *  \ingroup templateBanks
 *  \brief
 *  Generate templates for a directed F statistsic search for a single known 
 *  sky position with multiple frequency derivatives, with optimal-ish tiling
 *********************************************************************************/

/******** Includes ********/

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

/******** Constants ********/

static const int TRUE  = (1==1);
static const int FALSE = (1==0);

static const int BUFFER = 1024;

/******** Declarations ********/

RCSID("$Id$");

/******** Main function ********/

int main(int argc, char *argv[]) {

  UINT4 i;

  /* LAL status structure */
  LALStatus status = blank_status;

  /* Print help message? */
  BOOLEAN is_help = FALSE;

  /* Right ascension and declenation */
  REAL8 alpha = 6.12;
  REAL8 delta = 1.02;

  /* Frequency and frequency band */
  REAL8 freq = 100.0;
  REAL8 freq_band = 2e-2;
  
  /* Characteristic age of the target object */
  REAL8 age = 9.5e9;
  
  /* Number of spindowns */
  INT4 spindowns = 2;

  /* Minimum and maximum braking indices */
  INT4 minBrakingIndex = 2;
  INT4 maxBrakingIndex = 6;

  /* Spindowns and spindown bands */
  BOOLEAN useManualSpindowns = FALSE;
  REAL8 f1dot      = -freq / age;
  REAL8 f1dot_band =  freq / age;
  REAL8 f2dot      = 0.0;
  REAL8 f2dot_band = freq / pow(age, 2);
  REAL8 f3dot      = -freq / pow(age, 3);
  REAL8 f3dot_band =  freq / pow(age, 3);

  /* Type of parameter space tiling */
  CHAR *lattice_type = NULL;

  /* Maximum template mismatch */
  REAL8 mismatch = 0.25;

  /* Time span of the SFTs */
  REAL8 Tspan = 777600.0;

  /* Output files */
  CHAR *output_filename = NULL;
  FILE *output_file = NULL;

  /* Format string */
  CHAR *format_string = NULL;

  /* Flat lattice tiling structure */
  FlatLatticeTiling *tiling = NULL;

  /* Manual lower and upper bounds */
  gsl_vector *lower = NULL;
  gsl_vector *upper = NULL;

  /*
   *  Initialise variables
   */

  /* LAL error handler, debug level and log level */
  lal_errhandler = LAL_ERR_EXIT;
  LAL_CALL(LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LogSetLevel(lalDebugLevel);

  /* Lattice type */
  if ((lattice_type = (CHAR*) LALCalloc(10, sizeof(CHAR))) == NULL) {
    LALPrintError("%s\nERROR: LALCalloc failed to allocate lattice_type\n", rcsid);
    return EXIT_FAILURE;
  }
  strcpy(lattice_type, "Anstar");

  /* Format string */
  if ((format_string = (CHAR*) LALCalloc(10, sizeof(CHAR))) == NULL) {
    LALPrintError("%s\nERROR: LALCalloc failed to allocate format_string\n", rcsid);
    return EXIT_FAILURE;
  }
  strcpy(format_string, "%0.16g");

  /*
   *  Read in command line and check consistency
   */

  /* Register command line flags */
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",            'h', UVAR_HELP,      "Print this message", &is_help), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "Alpha",           'a', UVAR_OPTIONAL,  "Right ascension of the target object (in radians)", &alpha), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "Delta",           'd', UVAR_OPTIONAL,  "Declination of the target object (in radians)", &delta), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "Freq",            'f', UVAR_OPTIONAL,  "Frequency (in Hertz)", &freq), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "FreqBand",        'b', UVAR_OPTIONAL,  "Frequency band (in Hertz)", &freq_band), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "numSpindowns",    's', UVAR_OPTIONAL,  "Number of spindown parameters", &spindowns), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "Age",             'A', UVAR_OPTIONAL,  "Characteristic age of the target object (in seconds)", &age), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "minBrakingIndex", 'n', UVAR_OPTIONAL,  "Minimum braking index of the target object", &minBrakingIndex), &status);  
  LAL_CALL(LALRegisterINTUserVar   (&status, "maxBrakingIndex", 'N', UVAR_OPTIONAL,  "Maximum braking index of the target object", &maxBrakingIndex), &status);
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "manualSpindowns", 'S', UVAR_OPTIONAL,  "Use the following manual ranges on the spindown parameters:", &useManualSpindowns), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f1dot",            0 , UVAR_OPTIONAL,  "First spindown (in Hertz/s)", &f1dot), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f1dotBand",        0 , UVAR_OPTIONAL,  "First spindown band (in Hertz/s)", &f1dot_band), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f2dot",            0 , UVAR_OPTIONAL,  "Second spindown (in Hertz/s)", &f2dot), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f2dotBand",        0 , UVAR_OPTIONAL,  "Second spindown band (in Hertz/s)", &f2dot_band), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f3dot",            0 , UVAR_OPTIONAL,  "Third spindown (in Hertz/s)", &f3dot), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "f3dotBand",        0 , UVAR_OPTIONAL,  "Third spindown band (in Hertz/s)", &f3dot_band), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "latticeType",     'L', UVAR_OPTIONAL,  "Type of tiling lattice (cubic, Anstar)", &lattice_type), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "metricMismatch",  'X', UVAR_OPTIONAL,  "Maximum mismatch of the search templates", &mismatch), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "timeSpan",        'T', UVAR_OPTIONAL,  "Upper bound on time span of the SFTs (in seconds)", &Tspan), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "outputFile",      'o', UVAR_OPTIONAL,  "Output file containing the generated templates", &output_filename), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "formatString",    'F', UVAR_DEVELOPER, "Format string for outputting the template values", &format_string), &status);

  /* Read in command line */
  LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);
  if (is_help) {
    return EXIT_SUCCESS;
  }

  /*
   *  Set up lattice tiling over parameter space
   */

  /* Check if spindowns is not larger than hard limit */
  if (1 + spindowns > PULSAR_MAX_SPINS) {
    LALPrintError("%s\nERROR: Number of spindowns exceeds PULSAR_MAX_SPINS\n", rcsid);
    return EXIT_FAILURE;
  }

  /* Create tiling structure */
  if ((tiling = XLALCreateFlatLatticeTiling(1 + spindowns, mismatch)) == NULL) {
    LALPrintError("%s\nERROR: XLALCreateFlatLatticeTiling failed\n", rcsid);
    return EXIT_FAILURE;
  }
  
  /* Fill spindown metric and mismatch*/
  if (XLALSpindownMetric(tiling->metric, tiling->metric_scaling, (REAL8) Tspan) != XLAL_SUCCESS) {
    LALPrintError("%s\nERROR: XLALSpindownOnlyMetric failed\n", rcsid);
    return EXIT_FAILURE;
  }
  
  /* Fill parameter space */
  if (useManualSpindowns) {
    lower = gsl_vector_calloc(tiling->dimension);
    upper = gsl_vector_calloc(tiling->dimension);
    gsl_vector_set(lower, 0, freq);
    gsl_vector_set(upper, 0, freq + freq_band);
    if (tiling->dimension > 1) {
      gsl_vector_set(lower, 1, f1dot);
      gsl_vector_set(upper, 1, f1dot + f1dot_band);
    }
    if (tiling->dimension > 2) {
      gsl_vector_set(lower, 2, f2dot);
      gsl_vector_set(upper, 2, f2dot + f2dot_band);
    }
    if (tiling->dimension > 3) {
      gsl_vector_set(lower, 3, f2dot);
      gsl_vector_set(upper, 3, f2dot + f3dot_band);
    }
    if (XLALSquareParameterSpace(tiling, lower, upper) != XLAL_SUCCESS) {
      LALPrintError("%s\nERROR: XLALSquareParameterSpace failed\n", rcsid);
      return EXIT_FAILURE;
    }
  }
  else {
    if (XLALBrakingIndexParameterSpace(tiling, freq, freq_band, age, minBrakingIndex, maxBrakingIndex) != XLAL_SUCCESS) {
      LALPrintError("%s\nERROR: XLALBrakingIndexParameterSpace failed\n", rcsid);
      return EXIT_FAILURE;
    }
  }
  
  /* Fill lattice generator */
  if (strcmp(lattice_type, "cubic") == 0) {
    if (XLALCubicLatticeGenerator(tiling->generator, &tiling->generator_norm_thickness) != XLAL_SUCCESS) {
      LALPrintError("%s\nERROR: XLALCubicLatticeGenerator failed\n", rcsid);
      return EXIT_FAILURE;
    }
  }    
  else if (strcmp(lattice_type, "Anstar") == 0) {
    if (XLALAnstarLatticeGenerator(tiling->generator, &tiling->generator_norm_thickness) != XLAL_SUCCESS) {
      LALPrintError("%s\nERROR: XLALAnstarLatticeGenerator failed\n", rcsid);
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

  /*
   *  Open output file
   */

  if (output_filename != NULL) {
    if ((output_file = fopen(output_filename, "wt")) == NULL) {
      LALPrintError("%s\nERROR: Could not open output file '%s'\n", rcsid, output_filename);
      return EXIT_FAILURE;
    }
  }

  /*
   *  Step through templates
   */

  while (XLALNextFlatLatticePoint(tiling) == XLAL_SUCCESS) {

    /* If required write template to output file */
    if (output_file != NULL) {

      /* Print frequency, sky position and at least PULSAR_MAX_SPINS spindowns as required by CFS_v2 */
      fprintf(output_file, format_string, gsl_vector_get(tiling->current, 0));
      fprintf(output_file, " ");
      fprintf(output_file, format_string, alpha);
      fprintf(output_file, " ");
      fprintf(output_file, format_string, delta);
      for (i = 1; i < tiling->dimension; ++i) {
	fprintf(output_file, " ");
	fprintf(output_file, format_string, gsl_vector_get(tiling->current, i));
      }
      for (; i < PULSAR_MAX_SPINS; ++i) {
	fprintf(output_file, " 0");
      }
      fprintf(output_file, "\n");
	
    }

  }
  printf("Number of templates generated: %lli\n", tiling->template_count);

  /*
   *  Close output file
   */

  if (output_file != NULL) {
    fclose(output_file);
  }

  /*
   *  Cleanup
   */

  LAL_CALL(LALDestroyUserVars(&status), &status);

  XLALDestroyFlatLatticeTiling(tiling);

  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
