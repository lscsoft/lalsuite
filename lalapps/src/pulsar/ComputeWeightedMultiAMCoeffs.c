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
 * \brief Computes the weighted AM coefficients for given SFTs
 */

#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LogPrintf.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LALInitBarycenter.h>
#include <lal/DetectorStates.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/ComputeFstat.h>
#include <lalapps.h>

RCSID("$Id$");

#define TRUE  (1==1)
#define FALSE (1==0)

int main(int argc, char *argv[]) {

  LALStatus status = blank_status;
  BOOLEAN help = FALSE;

  REAL8 alpha = 0.0;
  REAL8 delta = 0.0;
  REAL8 freq = 0.0;
  REAL8 band = 0.0;
  INT4 num_bands = 1;
  CHAR *sft_pattern = NULL;
  CHAR *ephem_dir = NULL;
  CHAR *ephem_year = NULL;
  INT4 rng_med_win = 50;
  CHAR *output_file = NULL;

  FILE *fp = NULL;
  INT4 i = 0;
  CHAR *cmdline = NULL;
  EphemerisData ephemeris = empty_EphemerisData;
  
  /* Initialise LAL error handler, debug level and log level */
  lal_errhandler = LAL_ERR_EXIT;
  LAL_CALL(LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LogSetLevel(lalDebugLevel);
  
  /* Register command line arguments */
  LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",            'h', UVAR_HELP,     "Print this help message", &help), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "right-ascension", 'a', UVAR_REQUIRED, "Right ascension in radians", &alpha), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "declination",     'd', UVAR_REQUIRED, "Declination in radians", &delta), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "freq",            'f', UVAR_REQUIRED, "Starting frequency", &freq), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "band",            'b', UVAR_REQUIRED, "Frequency band", &band), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "num-bands",       'N', UVAR_OPTIONAL, "Number of bands to compute AM coeffs for", &num_bands), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "sft-patt",        'D', UVAR_REQUIRED, "File pattern of the input SFTs", &sft_pattern), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "ephem-dir",       'E', UVAR_OPTIONAL, "Directory containing ephemeris files", &ephem_dir), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "ephem-year",      'y', UVAR_REQUIRED, "Year suffix for ephemeris files", &ephem_year), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "rng-med-win",     'k', UVAR_OPTIONAL, "Size of the running median window", &rng_med_win), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "output-file",     'o', UVAR_OPTIONAL, "Output file for the AM coefficients (defaults to stdout)", &output_file), &status);

  /* Get command line arguments */
  LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);
  if (help)
    return EXIT_SUCCESS;

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
  fprintf(fp, "%% freq band A B C\n");

  /* Iterate over the bands */
  for (i = 0; i < num_bands; ++i) {
    
    REAL8 freq_i = freq + i*band;

    SFTCatalog *catalog = NULL;
    MultiSFTVector *sfts = NULL;
    MultiDetectorStateSeries *detector_states = NULL;
    MultiPSDVector *rng_med = NULL;
    MultiNoiseWeights *noise_weights = NULL;
    MultiAMCoeffs *AM_coeffs = NULL;
    
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
      fmin = freq_i - extra;
      fmax = freq_i + extra + band;
      
      /* Load the SFTs */
      LogPrintf(LOG_DEBUG, "Loading SFTs (%f to %f) ... ", fmin, fmax);
      LAL_CALL(LALLoadMultiSFTs(&status, &sfts, catalog, fmin, fmax), &status);
      LogPrintfVerbatim(LOG_DEBUG, "done\n");
    }
    
    /* Load the ephemeris data (if not already) */
    if (!(ephemeris.ephemE && ephemeris.ephemS)) {
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
	snprintf(ephemeris.ephiles.earthEphemeris, buf, "%s/earth%s.dat", ephem_dir, ephem_year);
	snprintf(ephemeris.ephiles.sunEphemeris, buf, "%s/sun%s.dat", ephem_dir, ephem_year);
      }
      else {
	snprintf(ephemeris.ephiles.earthEphemeris, buf, "earth%s.dat", ephem_year);
	snprintf(ephemeris.ephiles.sunEphemeris, buf, "sun%s.dat", ephem_year);
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
      LogPrintfVerbatim(LOG_DEBUG, "done\n");
    }      
    
    /* Print the AM coefficients (single-sided PSD in Sinv_Tsft !) */
    {
      const REAL8 A = AM_coeffs->Mmunu.Ad * AM_coeffs->Mmunu.Sinv_Tsft;
      const REAL8 B = AM_coeffs->Mmunu.Bd * AM_coeffs->Mmunu.Sinv_Tsft;
      const REAL8 C = AM_coeffs->Mmunu.Cd * AM_coeffs->Mmunu.Sinv_Tsft;
      fprintf(fp, "%0.3f %0.3f %0.4e %0.4e %0.4e\n", freq_i, band, A, B, C);
    }
    
    /* Cleanup */
    LAL_CALL(LALDestroySFTCatalog(&status, &catalog), &status);
    LAL_CALL(LALDestroyMultiSFTVector(&status, &sfts), &status);
    XLALDestroyMultiDetectorStateSeries(detector_states);
    LAL_CALL(LALDestroyMultiPSDVector(&status, &rng_med), &status);
    LAL_CALL(LALDestroyMultiNoiseWeights(&status, &noise_weights), &status);
    XLALDestroyMultiAMCoeffs(AM_coeffs);
    
  }

  /* Close the output file */
  if (LALUserVarWasSet(&output_file)) {
    fclose(fp);
  }
  
  /* Cleanup */
  LALDestroyUserVars(&status);
  XLALFree(ephemeris.ephiles.earthEphemeris);
  XLALFree(ephemeris.ephiles.sunEphemeris);
  LALFree(ephemeris.ephemE);
  LALFree(ephemeris.ephemS);
  LALCheckMemoryLeaks();
  
  return EXIT_SUCCESS;
  
}
