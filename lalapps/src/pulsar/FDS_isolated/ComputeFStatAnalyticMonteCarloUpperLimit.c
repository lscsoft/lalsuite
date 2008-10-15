/*
 *  Copyright (C) 2008 Karl Wette
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
 * \brief Computes an upper limit using Monte Carlo integration
          of the analytic F statistic signal model
 */

#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector_int.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/GSLSupport.h>
#include <lal/LogPrintf.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LALInitBarycenter.h>
#include <lal/DetectorStates.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/ComputeFstat.h>
#include <lal/ConfigFile.h>
#include <lalapps.h>

RCSID("$Id$");

REAL8 pdf_ncx2_4(REAL8, REAL8);
REAL8 d_pdf_ncx2_4(REAL8, REAL8);
REAL8 ran_ncx2_4(const gsl_rng*, REAL8);

#define TRUE  (1==1)
#define FALSE (1==0)

const REAL8 min_cosi = 0;
const REAL8 max_cosi = 1;
const REAL8 min_psi = -LAL_PI_4;
const REAL8 max_psi =  LAL_PI_4;

int main(int argc, char *argv[]) {

  LALStatus status = blank_status;
  BOOLEAN help = FALSE;

  REAL8 alpha = 0.0;
  REAL8 delta = 0.0;
  REAL8 freq = 0.0;
  REAL8 band = 0.0;
  REAL8 twoFs = 0.0;
  CHAR *mism_hist_file = NULL;
  REAL8 max_mismatch = 0.0;
  CHAR *sft_pattern = NULL;
  CHAR *ephem_dir = NULL;
  CHAR *ephem_year = NULL;
  INT4 rng_med_win = 50;
  REAL8 max_rel_err = 1.0e-3;
  REAL8 h0_brake = 0.75;
  INT4 MC_trial_init = 5e5;
  REAL8 MC_trial_fac = 1.5;
  INT4 MC_trial_reset = 5e7;
  REAL8 FDR = 0.05;
  CHAR *output_file = NULL;
  REAL8 twoF_pdf_hist_binw = 1.0;
  CHAR *twoF_pdf_hist_file = NULL;

  FILE *fp = NULL;
  CHAR *cmdline = NULL;
  gsl_matrix *mism_hist = NULL;
  SFTCatalog *catalog = NULL;
  MultiSFTVector *sfts = NULL;
  EphemerisData ephemeris = empty_EphemerisData;
  MultiDetectorStateSeries *detector_states = NULL;
  MultiPSDVector *rng_med = NULL;
  MultiNoiseWeights *noise_weights = NULL;
  MultiAMCoeffs *AM_coeffs = NULL;
  REAL8 A_coeff = 0.0;
  REAL8 B_coeff = 0.0;
  REAL8 C_coeff = 0.0;
  gsl_rng *rng = NULL;
  REAL8 h0 = 0.0;
  REAL8 h0_prev = 0.0;
  REAL8 dh0 = 0.0;
  INT4 h0_iter = 0;
  INT4 MC_trials = 0;
  INT4 MC_iter = 0;
  INT4 twoF_pdf_FD = 0;
  gsl_vector_int *twoF_pdf_hist = NULL;
    
  /* Initialise LAL error handler, debug level and log level */
  lal_errhandler = LAL_ERR_EXIT;
  LAL_CALL(LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LogSetLevel(lalDebugLevel);
  
  /* Register command line arguments */
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",             'h', UVAR_HELP,     "Print this help message", &help), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "alpha",            'a', UVAR_REQUIRED, "Right ascension in radians", &alpha), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "delta",            'd', UVAR_REQUIRED, "Declination in radians", &delta), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "freq",             'f', UVAR_REQUIRED, "Starting frequency", &freq), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "band",             'b', UVAR_REQUIRED, "Frequency band", &band), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "loudest-2F",       'F', UVAR_REQUIRED, "Loudest 2F value in this band", &twoFs), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "mism-hist-file",   'M', UVAR_OPTIONAL, "File containing the mismatch PDF histogram", &mism_hist_file), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "max-mismatch",     'm', UVAR_OPTIONAL, "Maximum mismatch to scale histogram to", &max_mismatch), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "sft-patt",         'D', UVAR_REQUIRED, "File pattern of the input SFTs", &sft_pattern), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "ephem-dir",        'E', UVAR_OPTIONAL, "Directory containing ephemeris files", &ephem_dir), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "ephem-year",       'y', UVAR_REQUIRED, "Year suffix for ephemeris files", &ephem_year), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "rng-med-win",      'k', UVAR_OPTIONAL, "Size of the running median window", &rng_med_win), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "max-rel-err",      'e', UVAR_OPTIONAL, "Maximum error in h0 relative to previous value", &max_rel_err), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "h0-brake",         'b', UVAR_OPTIONAL, "h0 cannot change by more than this fraction of itself", &h0_brake), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "init-MC-tri",      'I', UVAR_OPTIONAL, "Initial number of MC int. trials", &MC_trial_init), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "MC-tri-incr",      'i', UVAR_OPTIONAL, "Multiply number of MC int. trials by this after each step", &MC_trial_fac), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "MC-tri-reset",     'Z', UVAR_OPTIONAL, "Reset if no convergence after this number of MC int. trials", &MC_trial_reset), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "false-dism",       'R', UVAR_OPTIONAL, "Target false dismissal rate", &FDR), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "output-file",      'o', UVAR_OPTIONAL, "Output file for the upper limit and other info (defaults to stdout)", &output_file), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "2F-pdf-hist-binw", 'B', UVAR_OPTIONAL, "Bin width of the histogram of the non-central 2F distribution", &twoF_pdf_hist_binw), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "2F-pdf-hist-file", 'H', UVAR_OPTIONAL, "Output file for the histogram of the non-central 2F distribution", &twoF_pdf_hist_file), &status);

  /* Get command line arguments */
  LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);
  if (help)
    return EXIT_SUCCESS;

  /* Send log output to stdout */
  LogSetFile(stdout);

  /* Load the mismatch PDF histogram */
  if (LALUserVarWasSet(&mism_hist_file)) {

    size_t i;
    LALParsedDataFile *file = NULL;
    double total_count = 0;
    double max_bin = 0;

    /* Check the mismatch */
    if (max_mismatch <= 0.0) {
      LALPrintError("Maximum mismatch must be strictly positive\n");
      return EXIT_FAILURE;
    }
    
    /* Load the file */
    LAL_CALL(LALParseDataFile(&status, &file, mism_hist_file), &status);

    /* Allocate memory */
    if ((mism_hist = gsl_matrix_alloc(file->lines->nTokens, 3)) == NULL) {
      LALPrintError("Couldn't allocate a gsl_matrix\n");
      return EXIT_FAILURE;
    }

    /* Read in histogram */
    for (i = 0; i < mism_hist->size1; ++i) {
      double left_bin, right_bin, count;
      if (sscanf(file->lines->tokens[i], "%lf %lf %lf", &left_bin, &right_bin, &count) != 3) {
	LALPrintError("Couldn't parse line %i of '%s'\n", i, mism_hist_file);
	return EXIT_FAILURE;
      }
      if (left_bin >= right_bin || count < 0.0) {
	LALPrintError("Invalid syntax: line %i of '%s'\n", i, mism_hist_file);
	return EXIT_FAILURE;
      }
      gsl_matrix_set(mism_hist, i, 0, left_bin);
      gsl_matrix_set(mism_hist, i, 1, right_bin);
      gsl_matrix_set(mism_hist, i, 2, count);
      if (max_bin < right_bin)
	max_bin = right_bin;
      total_count += count;
    }

    /* Normalise histogram to unit maximum mismatch and unit area */
    for (i = 0; i < mism_hist->size1; ++i) {
      gsl_matrix_set(mism_hist, i, 0, gsl_matrix_get(mism_hist, i, 0) / max_bin);
      gsl_matrix_set(mism_hist, i, 1, gsl_matrix_get(mism_hist, i, 1) / max_bin);
      gsl_matrix_set(mism_hist, i, 2, gsl_matrix_get(mism_hist, i, 2) / total_count);
    }

    /* Cleanup */
    LAL_CALL(LALDestroyParsedDataFile(&status, &file), &status);

  }

  /* Load the SFTs */
  {
    SFTConstraints constraints = empty_SFTConstraints;
    REAL8 extra = 0.0, fmin = 0.0, fmax = 0.0;
    
    /* Load the catalog */
    LogPrintf(LOG_DEBUG, "Loading SFT catalog ... ");
    LAL_CALL(LALSFTdataFind(&status, &catalog, sft_pattern, &constraints), &status);
    if (!catalog || catalog->length == 0) {
      LALPrintError("Couldn't find SFTs matching '%s'\n", sft_pattern);
      return EXIT_FAILURE;
    }
    LogPrintfVerbatim(LOG_DEBUG, "done: %i SFTs starting at GPS %i\n",
		      catalog->length, catalog->data[0].header.epoch.gpsSeconds);
    
    /* Determine the frequency range */
    extra = catalog->data[0].header.deltaF * (rng_med_win/2 + 1);
    fmin = freq - extra;
    fmax = freq + extra + band;
    
    /* Load the SFTs */
    LogPrintf(LOG_DEBUG, "Loading SFTs (%f to %f) ... ", fmin, fmax);
    LAL_CALL(LALLoadMultiSFTs(&status, &sfts, catalog, fmin, fmax), &status);
    LogPrintfVerbatim(LOG_DEBUG, "done\n");
  }
  
  /* Load the ephemeris data */
  {
    size_t buf = strlen(ephem_year) + 20;
    if (ephem_dir)
      buf += strlen(ephem_dir);
    
    /* Allocate memory */
    if ((ephemeris.ephiles.earthEphemeris = (CHAR*)XLALCalloc(buf, sizeof(CHAR))) == NULL ||
	(ephemeris.ephiles.sunEphemeris = (CHAR*)XLALCalloc(buf, sizeof(CHAR))) == NULL) {
      LALPrintError("Couldn't allocate memory\n");
      return EXIT_FAILURE;
    }
    
    /* Create the file names */	
    if (ephem_dir) {
      LALSnprintf(ephemeris.ephiles.earthEphemeris, buf, "%s/earth%s.dat", ephem_dir, ephem_year);
      LALSnprintf(ephemeris.ephiles.sunEphemeris, buf, "%s/sun%s.dat", ephem_dir, ephem_year);
    }
    else {
      LALSnprintf(ephemeris.ephiles.earthEphemeris, buf, "earth%s.dat", ephem_year);
      LALSnprintf(ephemeris.ephiles.sunEphemeris, buf, "sun%s.dat", ephem_year);
    }
    
    /* Get the leap seconds */
    ephemeris.leap = (INT2)XLALGPSLeapSeconds(catalog->data[0].header.epoch.gpsSeconds);
    if (xlalErrno != XLAL_SUCCESS) {
      LALPrintError("XLALGPSLeapSeconds failed\n");
      return EXIT_FAILURE;
    }
    
    /* Load ephemeris */
    LogPrintf(LOG_DEBUG, "Loading ephemeris ... ");
    LAL_CALL(LALInitBarycenter(&status, &ephemeris), &status);
    LogPrintfVerbatim(LOG_DEBUG, "done\n");
  }
  
  /* Get the detector states */
  LogPrintf(LOG_DEBUG, "Calculating detector states ... ");
  LAL_CALL(LALGetMultiDetectorStates(&status, &detector_states, sfts, &ephemeris), &status);
  LogPrintfVerbatim(LOG_DEBUG, "done\n");
  
  /* Normalise SFTs and compute noise weights */
  {
    LogPrintf(LOG_DEBUG, "Normalising SFTs and computing noise weights ... ");
    LAL_CALL(LALNormalizeMultiSFTVect(&status, &rng_med, sfts, rng_med_win), &status);
    LAL_CALL(LALComputeMultiNoiseWeights(&status, &noise_weights, rng_med, rng_med_win, 0), &status);
    LogPrintfVerbatim(LOG_DEBUG, "done\n");
  }
  
  /* Compute the AM coefficients */
  {
    SkyPosition sky;
    sky.system = COORDINATESYSTEM_EQUATORIAL;
    sky.longitude = alpha;
    sky.latitude = delta;
    LogPrintf(LOG_DEBUG, "Calculating AM coefficients ... ");
    LAL_CALL(LALGetMultiAMCoeffs(&status, &AM_coeffs, detector_states, sky), &status);
    if (XLALWeighMultiAMCoeffs(AM_coeffs, noise_weights) != XLAL_SUCCESS) {
      LALPrintError("XLALWeighMultiAMCoeffs failed\n");
      return EXIT_FAILURE;
    }
    A_coeff = AM_coeffs->Mmunu.Ad * AM_coeffs->Mmunu.Sinv_Tsft;
    B_coeff = AM_coeffs->Mmunu.Bd * AM_coeffs->Mmunu.Sinv_Tsft;
    C_coeff = AM_coeffs->Mmunu.Cd * AM_coeffs->Mmunu.Sinv_Tsft;
    /* DEBUG Single-side PSD correction */
    A_coeff *= 0.5;
    B_coeff *= 0.5;
    C_coeff *= 0.5;
    /* DEBUG Single-side PSD correction */
    LogPrintfVerbatim(LOG_DEBUG, "done: A = %0.4e, B = %0.4e, C = %0.4e\n", A_coeff, B_coeff, C_coeff);
  }      
  
  /* Initialise the random number generator */
  {
    unsigned long seed;
    FILE *fpr = NULL;
    if ((rng = gsl_rng_alloc(gsl_rng_mt19937)) == NULL) {
      LALPrintError("Couldn't allocate a gsl_rng\n");
      return EXIT_FAILURE;
    }
    if ((fpr = fopen("/dev/random", "r")) == NULL) {
      LALPrintError("Couldn't open '/dev/random'\n");
      return EXIT_FAILURE;
    }
    fread(&seed, sizeof(seed), 1, fpr);
    fclose(fpr);
    gsl_rng_set(rng, seed);
  }
  
  /* Begin iterations to find h0 */
  do {

    /* Compute the volume of the Monte Carlo integration volume:
       integrate over twoF but average over cosi,psi (hence not present) */
    const REAL8 MC_int_vol = (twoFs - 0);

    /* Open the output file */
    if (LALUserVarWasSet(&output_file)) {
      if ((fp = fopen(output_file, "wb")) == NULL) {
	LALPrintError("Couldn't open output file '%s'\n", output_file);
	return EXIT_FAILURE;
      }
    }
    else {
      fp = stdout;
    }
    LAL_CALL(LALUserVarGetLog(&status, &cmdline, UVAR_LOGFMT_CMDLINE), &status);
    fprintf(fp, "%%%% %s\n%%%% %s\n", rcsid, cmdline);
    LALFree(cmdline);
    cmdline = NULL;
    fprintf(fp, "freq=%0.4f band=%0.4f\n", freq, band);
  
    /* Compute first guess at h0 */
    h0  = twoFs * pow(A_coeff * B_coeff - C_coeff * C_coeff, -0.25);
    h0 *= GSL_MAX(0.5, 1.0 + gsl_ran_gaussian(rng, 0.1));
  
    /* Begin Newton-Raphson iteration to find h0 */
    h0_prev = GSL_POSINF;
#define H0_ERROR fabs(1.0 - (h0 / h0_prev))
    h0_iter = 0;
    MC_trials = MC_trial_init;
    while (H0_ERROR > max_rel_err && MC_trials < MC_trial_reset) {
      
      /* Integrand and its derivative w.r.t. h0, and total mismatch probability */
      REAL8 J = 0.0;
      REAL8 dJ = 0.0;
      REAL8 total_prob_mism = 0.0;

      /* False dismissal rate from 2F distribution */    
      REAL8 twoF_pdf_FDR = 0.0;
      
      LogPrintf(LOG_DEBUG, "Beginning h0 loop %2i with h0=%0.4e, error=%0.4e, MC_trials=%i\n", h0_iter, h0, H0_ERROR, MC_trials);
      
      fprintf(fp, "MC_trials=%i ", MC_trials);

      /* Destroy any previous histogram */
      if (twoF_pdf_hist) {
	gsl_vector_int_free(twoF_pdf_hist);
	twoF_pdf_hist = NULL;
      }

      /* Begin Monte Carlo integration to find J and dJ */
      MC_iter = MC_trials;
      twoF_pdf_FD = 0;
      while (MC_iter--) {

	REAL8 twoF_mism_factor = 1.0;
	REAL8 prob_mism = 1.0;

	REAL8 cosi = 0.0;
	REAL8 psi = 0.0;
	REAL8 twoF = 0.0;
	REAL8 A1 = 0.0;
	REAL8 A2 = 0.0;
	REAL8 A3 = 0.0;
	REAL8 A4 = 0.0;
	REAL8 rho2 = 0.0;

	REAL8 twoF_from_pdf = 0.0;

	/* Generate random mismatch factor */
	if (LALUserVarWasSet(&mism_hist_file)) {

	  size_t i;
	  REAL8 norm_mism;

	  /* Generate random normalised mismatch */
	  norm_mism = gsl_ran_flat(rng, 0.0, 1.0);

	  /* Find bin in histogram */
	  for (i = 0; i < mism_hist->size1; ++i) {
	    const double left_bin  = gsl_matrix_get(mism_hist, i, 0);
	    const double right_bin = gsl_matrix_get(mism_hist, i, 1);
	    const double prob_bin  = gsl_matrix_get(mism_hist, i, 2);
	    if (left_bin <= norm_mism && norm_mism < right_bin) {

	      /* Create mismatch from centre of bin */
	      const double mismatch = 0.5 * (left_bin + right_bin) * max_mismatch;

	      /* 2F is increase by this factor */
	      twoF_mism_factor = 1.0 / (1.0 - mismatch);

	      /* This factor has this probability */
 	      prob_mism = prob_bin;

	      break;
	    }
	  }

	}	

	/* Generate random cosi, psi, and twoF */
	cosi = gsl_ran_flat(rng, min_cosi, max_cosi);
	psi = gsl_ran_flat(rng, min_psi, max_psi);
	twoF = gsl_ran_flat(rng, 0, twoF_mism_factor * twoFs);
	
	/* Compute the amplitude coefficients vector A */
	A1 =  0.5 * h0 * (1 + cosi * cosi) * cos(2 * psi);
	A2 =  0.5 * h0 * (1 + cosi * cosi) * sin(2 * psi);
	A3 = -1.0 * h0 *             cosi  * sin(2 * psi);
	A4 =  1.0 * h0 *             cosi  * cos(2 * psi);
	
	/* Compute the optimal signal to noise ratio rho^2 */
	rho2 = (
		A1 * A1 * A_coeff +
		A2 * A2 * B_coeff +
		A1 * A2 * C_coeff +
		A2 * A1 * C_coeff +
		A3 * A3 * A_coeff +
		A4 * A4 * B_coeff +
		A3 * A4 * C_coeff +
		A4 * A3 * C_coeff
		);
	
	/* Add to the integrand and its derivative w.r.t. h0, and to total mismatch probability */
	J  += pdf_ncx2_4(rho2, twoF) * prob_mism;
	dJ += 2.0 * rho2 / h0 * d_pdf_ncx2_4(rho2, twoF) * prob_mism;
	total_prob_mism += prob_mism;
	
	/* Break if J and dJ failed */
	if (gsl_isnan(J) || gsl_isnan(dJ))
	  break;
	
	/* Count 2F value if it is below threshold */
	if (twoF_from_pdf < twoFs)
	  ++twoF_pdf_FD;
	
	/* Draw a 2F value from the non-central chi-square with parameter rho^2 */
	twoF_from_pdf = twoF_mism_factor * ran_ncx2_4(rng, rho2);
	
	/* Add 2F to histogram if needed */
	if (LALUserVarWasSet(&twoF_pdf_hist_file)) {
	  
	  /* Compute bin */
	  const size_t bin = twoF_from_pdf / twoF_pdf_hist_binw;
	  
	  /* Resize histogram vector if needed */
	  if (!twoF_pdf_hist || bin >= twoF_pdf_hist->size)
	    if (NULL == (twoF_pdf_hist = XLALResizeGSLVectorInt(twoF_pdf_hist, bin + 1, 0))) {
	      LALPrintError("\nCouldn't (re)allocate 'twoF_pdf_hist'\n");
	      return EXIT_FAILURE;
	    }
	  
	  /* Add to bin */
	  gsl_vector_int_set(twoF_pdf_hist, bin,
			     gsl_vector_int_get(twoF_pdf_hist, bin) + 1);
	  
	}

      }
      
      /* If J and dJ failed, reduce h0 and try again */
      if (gsl_isnan(J) || gsl_isnan(dJ)) {
	h0 /= 2.0;
	LogPrintf(LOG_DEBUG, "Reducing h0 to %0.4e and starting again because either J=%0.4e or dJ=%0.4e\n", h0, J, dJ);
	h0_iter = 0;
	MC_trials = MC_trial_init;
	continue;
      }

      /* Finalise the integrand and derivative */
      J  *= MC_int_vol / total_prob_mism;
      dJ *= MC_int_vol / total_prob_mism;
      
      /* Compute the false dismissal rate from 2F distribution */
      twoF_pdf_FDR = (1.0 * twoF_pdf_FD) / MC_trials;

      /* Compute the increment in h0 from Newton-Raphson */
      dh0 = (FDR - J) / dJ;
      
      /* Limit the increment in h0 to |h0 * h0_brake| */
      dh0 = GSL_SIGN(dh0) * GSL_MIN(fabs(dh0), fabs(h0 * h0_brake));
      
      /* Increment h0 (saving old value) */
      h0_prev = h0;
      h0 += dh0;
      
      fprintf(fp, "h0=%0.4e FDR=%0.4f\n", h0, twoF_pdf_FDR);
      fflush(fp);
      
      LogPrintf(LOG_DEBUG, "Ending    h0 loop %2i with h0=%0.4e, error=%0.4e, J-FDR=% 0.4e, dh0=% 0.4e\n", h0_iter, h0, H0_ERROR, J - FDR, dh0);

      /* Increase iteration count and number of MC trials */
      ++h0_iter;
      MC_trials *= MC_trial_fac;

    }
#undef H0_ERROR
  
    /* Close the output file */
    if (LALUserVarWasSet(&output_file)) {
      fprintf(fp, "%%DONE\n");
      fclose(fp);
    }

    /* If number of MC trials exceeded reset */
    if (MC_trials >= MC_trial_reset)
      LogPrintf(LOG_DEBUG, "Failed to converge after %i iterations (MC_trails=%i): trying again ...\n", h0_iter, MC_trials);

  } while (MC_trials >= MC_trial_reset);
  
  /* Write 2F histogram if needed */
  if (LALUserVarWasSet(&twoF_pdf_hist_file)) {

    size_t i;
    FILE *fpH;

    if ((fpH = fopen(twoF_pdf_hist_file, "wb")) == NULL) {
      LALPrintError("Couldn't open histogram file '%s'\n", twoF_pdf_hist_file);
      return EXIT_FAILURE;
    }
    LAL_CALL(LALUserVarGetLog(&status, &cmdline, UVAR_LOGFMT_CMDLINE), &status);
    fprintf(fpH, "%%%% %s\n%%%% %s\n", rcsid, cmdline);
    LALFree(cmdline);
    cmdline = NULL;
    
    for (i = 0; i < twoF_pdf_hist->size; ++i)
      fprintf(fpH, "%0.3g %0.3g %i\n",
	      twoF_pdf_hist_binw * i,
	      twoF_pdf_hist_binw * (i + 1),
	      gsl_vector_int_get(twoF_pdf_hist, i));
    
    fprintf(fpH, "%%DONE\n");
    fclose(fpH);
    
  }

  /* Cleanup */
  LAL_CALL(LALDestroyUserVars(&status), &status);
  if (mism_hist)
    gsl_matrix_free(mism_hist);
  LAL_CALL(LALDestroySFTCatalog(&status, &catalog), &status);
  LAL_CALL(LALDestroyMultiSFTVector(&status, &sfts), &status);
  XLALFree(ephemeris.ephiles.earthEphemeris);
  XLALFree(ephemeris.ephiles.sunEphemeris);
  LALFree(ephemeris.ephemE);
  LALFree(ephemeris.ephemS);
  XLALDestroyMultiDetectorStateSeries(detector_states);
  LAL_CALL(LALDestroyMultiPSDVector(&status, &rng_med), &status);
  LAL_CALL(LALDestroyMultiNoiseWeights(&status, &noise_weights), &status);
  XLALDestroyMultiAMCoeffs(AM_coeffs);
  gsl_rng_free(rng);
  if (twoF_pdf_hist)
    gsl_vector_int_free(twoF_pdf_hist);
  LALCheckMemoryLeaks();
  
  return EXIT_SUCCESS;
  
}

/* PDF of non-central chi-square with 4 degrees 
   of freedom and non-centrality parameter lambda */
REAL8 pdf_ncx2_4(REAL8 lambda, REAL8 x) {

  const REAL8 z = sqrt(x * lambda);

  gsl_sf_result I1;
  gsl_error_handler_t *h = NULL;

  /* Compute the Bessel functions */
  h = gsl_set_error_handler_off();
  if (gsl_sf_bessel_In_e(1, z, &I1) != GSL_SUCCESS) {
    gsl_set_error_handler(h);
    return GSL_NAN;
  }
  gsl_set_error_handler(h);
  
  /* Compute the PDF */
  return 0.5 * exp(-0.5 * (x + lambda)) * sqrt(x / lambda) * I1.val;

}

/* Derivative of the PDF of non-central chi-square with 4 degrees 
   of freedom and non-centrality parameter lambda w.r.t. lambda */
REAL8 d_pdf_ncx2_4(REAL8 lambda, REAL8 x) {

  const REAL8 z = sqrt(x * lambda);

  int i;
  gsl_sf_result I[3];
  gsl_error_handler_t *h = NULL;

  /* Compute the Bessel functions */
  h = gsl_set_error_handler_off();
  for (i = 0; i < 3; ++i)
    if (gsl_sf_bessel_In_e(i, z, &I[i]) != GSL_SUCCESS) {
      gsl_set_error_handler(h);
      return GSL_NAN;
    }
  gsl_set_error_handler(h);
  
  /* Compute the derivative of the PDF */
  return 0.25 * exp(-0.5 * (x + lambda)) * (0.5 * x / lambda * (I[0].val + I[2].val) - sqrt(x / lambda) * (1 + 1 / lambda) * I[1].val);
  
}

/* Random number drawn from a non-central chi-square with 4 degrees 
   of freedom and non-centrality parameter lambda */
REAL8 ran_ncx2_4(const gsl_rng *rng, REAL8 lambda) {

  const REAL8 a = sqrt(lambda / 4);

  int i;
  REAL8 x = 0.0;

  for (i = 0; i < 4; ++i)
    x += pow(gsl_ran_gaussian(rng, 1.0) + a, 2.0);
  
  return x;

}
