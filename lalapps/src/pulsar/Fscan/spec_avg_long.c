/*
*  Copyright (C) 2007 Gregory Mendell
*                2020 Evan Goetz
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
* \file
* \ingroup lalapps_pulsar_fscan
*/

#include "config.h"

/*LAL header files*/
#include <LALAppsVCSInfo.h>
#include <lal/SFTutils.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdio.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>


/*---------- DEFINES ----------*/

/*----- Macros ----- */

/*---------- internal types ----------*/


///////////EXTRA FXN DEFS//////////////
int check_and_mark_line(const INT4Vector *line_freq_bin_in_sft, const UINT4 i, const UINT4 j, const UINT4VectorSequence *line_exceeds_thresh_in_sfts);


int main(int argc, char **argv)
{
    FILE *SPECOUT = NULL;
    FILE *WTOUT = NULL;
    FILE *LINEOUT = NULL;

    SFTCatalog *catalog = NULL;
    SFTVector *sft_vect = NULL;
    SFTConstraints XLAL_INIT_DECL(constraints);
    LIGOTimeGPS startTime, endTime;
    REAL8Vector *timeavg = NULL, *timeavgwt = NULL, *sumweight = NULL, *persistency = NULL;
    REAL8 f0 = 0, deltaF = 0;
    CHAR outbase[256], outfile0[512], outfile1[512], outfile2[512];

    CHAR *SFTpatt = NULL, *IFO = NULL, *outputBname = NULL;
    INT4 startGPS = 0, endGPS = 0, blocksRngMean = 21;
    REAL8 f_min = 0.0, f_max = 0.0, timebaseline = 0;

    LALStringVector *line_freq = NULL;
    REAL8Vector *line_freq_array = NULL;

    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&SFTpatt,      "SFTs",         STRING, 'p', REQUIRED, "SFT location/pattern" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&IFO,          "IFO",          STRING, 'I', REQUIRED, "Detector" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&startGPS,     "startGPS",     INT4,   's', REQUIRED, "Starting GPS time (SFT timestamps must be >= this)" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&endGPS,       "endGPS",       INT4,   'e', REQUIRED, "Ending GPS time (SFT timestamps must be < this)" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&f_min,        "fMin",         REAL8,  'f', REQUIRED, "Minimum frequency" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&f_max,        "fMax",         REAL8,  'F', REQUIRED, "Maximum frequency" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&blocksRngMean, "blocksRngMean", INT4,   'w', OPTIONAL, "Running Median window size") == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&outputBname,  "outputBname",  STRING, 'o', OPTIONAL, "Base name of output files" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&timebaseline, "timeBaseline", REAL8,  't', REQUIRED, "The time baseline of sfts") == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&line_freq,    "lineFreq", STRINGVector,  0, OPTIONAL, "CSV list of line frequencies. If set, then an output file with all GPS start times of SFTs with 0s and 1s is given for when line is above threshold (1 indicates above threshold)") == XLAL_SUCCESS, XLAL_EFUNC);

    BOOLEAN should_exit = 0;
    XLAL_CHECK_MAIN(XLALUserVarReadAllInput(&should_exit, argc, argv, lalAppsVCSInfoList) == XLAL_SUCCESS, XLAL_EFUNC);
    if (should_exit)
    {
        return(1);
    }

    printf("Starting spec_avg_long...\n");

    XLAL_CHECK_MAIN( blocksRngMean % 2 == 1, XLAL_EINVAL, "Need to provide an odd value for blocksRngMean");
    INT4 nside = (blocksRngMean - 1) / 2;

    INT4Vector *line_freq_bin_in_sft = NULL;
    if (XLALUserVarWasSet(&line_freq)) {
       XLAL_CHECK_MAIN( (line_freq_array = XLALCreateREAL8Vector(line_freq->length)) != NULL, XLAL_EFUNC );
       XLAL_CHECK_MAIN( (line_freq_bin_in_sft = XLALCreateINT4Vector(line_freq->length)) != NULL, XLAL_EFUNC );
       for (UINT4 n=0; n<line_freq->length; n++) {
	  line_freq_array->data[n] = atof(line_freq->data[n]);
       }
    }

    //Provide the constraints to the catalog
    startTime.gpsSeconds = startGPS;
    startTime.gpsNanoSeconds = 0;
    constraints.minStartTime = &startTime;
    endTime.gpsSeconds = endGPS;
    endTime.gpsNanoSeconds = 0;
    constraints.maxStartTime = &endTime;
    constraints.detector = IFO;

    printf("Calling XLALSFTdataFind with SFTpatt=%s\n", SFTpatt);
    XLAL_CHECK_MAIN( (catalog = XLALSFTdataFind ( SFTpatt, &constraints )) != NULL, XLAL_EFUNC );

    //Ensure that some SFTs were found given the start and end time and IFO constraints
    XLAL_CHECK_MAIN( catalog->length > 0, XLAL_EFAILED, "No SFTs found, please examine start time, end time, frequency range, etc.");

    printf("Now have SFT catalog with %d catalog files\n", catalog->length);

    if (XLALUserVarWasSet(&outputBname)) strcpy(outbase, outputBname);
    else snprintf(outbase, sizeof(outbase), "spec_%.2f_%.2f_%s_%d_%d", f_min, f_max, constraints.detector, startTime.gpsSeconds, endTime.gpsSeconds);
	
    snprintf(outfile0, sizeof(outfile0), "%s.txt", outbase);
    snprintf(outfile1, sizeof(outfile1), "%s_PWA.txt", outbase);
    snprintf(outfile2, sizeof(outfile2), "%s_line_times.csv", outbase);

    SPECOUT = fopen(outfile0, "w");
    WTOUT = fopen(outfile1, "w");

    UINT4VectorSequence *line_exceeds_thresh_in_sfts = NULL;
    if (XLALUserVarWasSet(&line_freq)) {
       LINEOUT = fopen(outfile2, "w");
       XLAL_CHECK_MAIN( (line_exceeds_thresh_in_sfts = XLALCreateUINT4VectorSequence(line_freq_array->length, catalog->length)) != NULL, XLAL_EFUNC );
       memset(line_exceeds_thresh_in_sfts->data, 0, sizeof(UINT4)*line_freq_array->length*catalog->length);
    }

    printf("Looping over SFTs to compute average spectra\n");
    for (UINT4 j = 0; j<catalog->length; j++)
    {
        fprintf(stderr,"Extracting SFT %d...\n", j);

	//Extract one SFT at a time from the catalog
	//we do this by using a catalog timeslice to get just the current SFT
	SFTCatalog XLAL_INIT_DECL(catalogSlice);
	LIGOTimeGPS thisSFTstarttime = catalog->data[j].header.epoch;
	LIGOTimeGPS thisSFTendtime = thisSFTstarttime;
	XLAL_CHECK_MAIN( XLALGPSAdd(&thisSFTendtime, 0.01) != NULL, XLAL_EFUNC );
	XLAL_CHECK_MAIN( XLALSFTCatalogTimeslice(&catalogSlice, catalog, &thisSFTstarttime, &thisSFTendtime) == XLAL_SUCCESS, XLAL_EFUNC );
	XLAL_CHECK_MAIN( (sft_vect = XLALLoadSFTs(&catalogSlice, f_min, f_max)) != NULL, XLAL_EFUNC );
	XLAL_CHECK_MAIN( sft_vect->length == 1, XLAL_EBADLEN, "Oops, got %d SFTs instead of one", sft_vect->length );

	//For the first time through the loop, we allocate some vectors
        if (j == 0)
	{
	    UINT4 numBins = sft_vect->data->data->length;
	    f0 = sft_vect->data->f0;
	    deltaF = sft_vect->data->deltaF;
	    printf("numBins=%d, f0=%f, deltaF=%f\n", numBins, f0, deltaF);

	    XLAL_CHECK_MAIN( (timeavg = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (timeavgwt = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (sumweight = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (persistency = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );

	    if (XLALUserVarWasSet(&line_freq)) {
	       for (UINT4 n=0; n<line_freq_array->length; n++) {
		  line_freq_bin_in_sft->data[n] = (INT4)round((line_freq_array->data[n]-f0)/deltaF);
	       }
	    }
	}

	//Loop over the SFT bins
	for (UINT4 i = 0; i < sft_vect->data->data->length; i++)
	{
	    REAL8 thispower = ((REAL8)crealf(sft_vect->data[0].data->data[i]) *
			       (REAL8)crealf(sft_vect->data[0].data->data[i])) +
		              ((REAL8)cimagf(sft_vect->data[0].data->data[i]) *
			       (REAL8)cimagf(sft_vect->data[0].data->data[i]));
	    REAL8 thisavepower = 0.;
	    UINT4 count = 0;
	    for (INT4 ii = -nside; ii <= nside; ii++)
	    {
	        //Only add to the cumulative average power if the variables are in range
	       if ((INT4)i + ii >= 0 && (INT4)i + ii<(INT4)sft_vect->data->data->length)
		{
		    thisavepower += ((REAL8)crealf(sft_vect->data[0].data->data[i + ii]) *
				     (REAL8)crealf(sft_vect->data[0].data->data[i + ii])) +
				    ((REAL8)cimagf(sft_vect->data[0].data->data[i + ii]) *
				     (REAL8)cimagf(sft_vect->data[0].data->data[i + ii]));
		    count++;
		}
	    }
	    thisavepower /= count;
	    REAL8 weight = 1. / thisavepower;

	    //For the first SFT, just assign the values to the vector, otherwise accumulate
	    if (j == 0)
	    {
	        timeavg->data[i] = thispower;
		timeavgwt->data[i] = thispower*weight;
		sumweight->data[i] = weight;
		// Always check if this bin is above the 3 sigma threshold
		if (thispower >= 3.0*thisavepower) {
		   persistency->data[i] = 1.0;
		   // If user specified a list of line frequencies to monitor
		   if (XLALUserVarWasSet(&line_freq)) {
		      XLAL_CHECK_MAIN( check_and_mark_line(line_freq_bin_in_sft, i, j, line_exceeds_thresh_in_sfts) == XLAL_SUCCESS, XLAL_EFUNC );
		   }
		} else persistency->data[i] = 0.0;
	    }
	    else
	    {
	        timeavg->data[i] += thispower;
		timeavgwt->data[i] += thispower*weight;
		sumweight->data[i] += weight;
		// Always check if this bin is above the 3 sigma threshold
		if (thispower >= 3.0*thisavepower) {
		   persistency->data[i] += 1.0;
		   // If user specified a list of line frequencies to monitor
		   if (XLALUserVarWasSet(&line_freq)) {
		      XLAL_CHECK_MAIN( check_and_mark_line(line_freq_bin_in_sft, i, j, line_exceeds_thresh_in_sfts) == XLAL_SUCCESS, XLAL_EFUNC );
		   }
		}
	    }
	}
	// Destroys current SFT Vector
	XLALDestroySFTVector(sft_vect);
	sft_vect = NULL;
    }
    printf("About to do calculation of averages...\n");
    printf("Sample: timeavg[0]=%g, timeavgwt[0]=%g, sumweight[0]=%g\n", timeavg->data[0], timeavgwt->data[0], sumweight->data[0]);
    /*timeavg records the power of each bin*/
    for (UINT4 i = 0; i < timeavg->length; i++)
    {
        REAL8 f = f0 + ((REAL4)i)*deltaF;
	REAL8 PSD = 2.*timeavg->data[i] / ((REAL4)catalog->length) / timebaseline;
	REAL8 PSDWT = 2.*timeavgwt->data[i] / sumweight->data[i] / timebaseline;
	REAL8 AMPPSD = pow(PSD, 0.5);
	REAL8 AMPPSDWT = pow(PSDWT, 0.5);
	REAL8 persist = persistency->data[i] / ((REAL8)catalog->length);
	fprintf(SPECOUT, "%16.8f %g %g %g %g %g\n", f, PSD, AMPPSD, PSDWT, AMPPSDWT, persist);

	REAL8 PWA_TAVGWT = timeavgwt->data[i];
	REAL8 PWA_SUMWT = sumweight->data[i];
	fprintf(WTOUT, "%16.8f %g %g\n", f, PWA_TAVGWT, PWA_SUMWT);
    }

    // If user specified a list of line frequencies to monitor then print those data to a file
    if (XLALUserVarWasSet(&line_freq)) {
       fprintf(LINEOUT, "# GPS time");
       for (UINT4 n=0; n<line_freq_array->length; n++) {
	  fprintf(LINEOUT,",%g Hz", line_freq_array->data[n]);
       }
       fprintf(LINEOUT, "\n");
       for (UINT4 i = 0; i < catalog->length; i++) {
	  fprintf(LINEOUT, "%u", catalog->data[i].header.epoch.gpsSeconds);
	  for (UINT4 n=0; n<line_freq_array->length; n++) {
	     fprintf(LINEOUT, ",%u", line_exceeds_thresh_in_sfts->data[n*line_exceeds_thresh_in_sfts->vectorLength + i]);
	  }
	  fprintf(LINEOUT, "\n");
       }
    }

	/*------------------------------------------------------------------------------------------------------------------------*/

    fprintf(stderr,"Destroying Variables\n");
    XLALDestroySFTCatalog(catalog);

    XLALDestroyREAL8Vector(timeavg);
    XLALDestroyREAL8Vector(timeavgwt);
    XLALDestroyREAL8Vector(sumweight);
    XLALDestroyREAL8Vector(persistency);

    /*close all the files, spec_avg.c is done, all info written to the files.*/
    fprintf(stderr,"Closing Files\n");
    fclose(SPECOUT);
    fclose(WTOUT);
    if (XLALUserVarWasSet(&line_freq)) {
       fclose(LINEOUT);
       XLALDestroyUINT4VectorSequence(line_exceeds_thresh_in_sfts);
       XLALDestroyREAL8Vector(line_freq_array);
       XLALDestroyINT4Vector(line_freq_bin_in_sft);
    }

    XLALDestroyUserVars();
    fprintf(stderr,"Done Destroying Variables\n");
    fprintf(stderr, "end of spec_avg\n");
    fprintf(stderr, "Spec_avg_done!\n");

    return(0);

}
/* END main */


int check_and_mark_line(const INT4Vector *line_freq_bin_in_sft, const UINT4 i, const UINT4 j, const UINT4VectorSequence *line_exceeds_thresh_in_sfts)
{
    // Loop over the list of lines, checking if the user-specified line is at
    // least equal or larger than the 0th frequency bin, and also if the
    // current bin we are looking at (i) is equal to the frequency specified
    // by the user in units of SFT bins
    // If it is equal, then mark the bin in the j-th SFT
    for (UINT4 n=0; n<line_freq_bin_in_sft->length; n++) {
        if (line_freq_bin_in_sft->data[n]>=0 && i==(UINT4)line_freq_bin_in_sft->data[n]) {
	    line_exceeds_thresh_in_sfts->data[n*line_exceeds_thresh_in_sfts->vectorLength + j] = 1;
	}
    }

    return XLAL_SUCCESS;
}
