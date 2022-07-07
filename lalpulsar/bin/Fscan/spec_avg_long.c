/*
*  Copyright (C) 2007 Gregory Mendell
*                2020-2022 Evan Goetz
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
* \ingroup lalpulsar_bin_fscan
*/

#include "config.h"

#include <lal/SFTutils.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdio.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LALPulsarVCSInfo.h>

/*---------- DEFINES ----------*/
#define POWER(x) (((REAL8)crealf(x)*(REAL8)crealf(x)) + ((REAL8)cimagf(x)*(REAL8)cimagf(x)))

/*----- Macros ----- */

/*---------- internal types ----------*/


///////////EXTRA FXN DEFS//////////////
LIGOTimeGPSVector * setup_epochs(const SFTCatalog *catalog, const INT4 persistAvgOpt, const BOOLEAN persistAvgOptWasSet, const INT4 persistAvgSeconds);
SFTVector * extract_one_sft(const SFTCatalog *full_catalog, const UINT4 sft_index, const REAL8 f_min, const REAL8 f_max);
int set_sft_avg_epoch(struct tm *utc, LIGOTimeGPS *epoch_start, const LIGOTimeGPS first_sft_epoch, const INT4 persistAvgOpt, const BOOLEAN persistAvgOptWasSet);
int validate_line_freq(LALStringVector **line_freq, const REAL8 f0, const REAL8 deltaF, const UINT4 numBins);
REAL8Vector * line_freq_str2dbl(const LALStringVector *line_freq);
int rngmean(const REAL8Vector *input, const REAL8Vector *output, const INT4 blocksRngMean);
int rngstd(const REAL8Vector *input, const REAL8Vector *means, const REAL8Vector *output, const INT4 blocksRngMean);
int select_mean_std_from_vect(REAL8 *mean, REAL8 *std, const REAL8Vector *means, const REAL8Vector *stds, const UINT4 idx, const UINT4 nside);


int main(int argc, char **argv)
{
    FILE *SPECOUT = NULL, *WTOUT = NULL, *LINEOUT = NULL;

    SFTCatalog *catalog = NULL;
    SFTConstraints XLAL_INIT_DECL(constraints);
    LIGOTimeGPS startTime, endTime;
    REAL8Vector *timeavg = NULL, *timeavgwt = NULL, *sumweight = NULL, *persistency = NULL;
    REAL8 f0 = 0, deltaF = 0;
    CHAR outbase[256], outfile0[512], outfile1[512], outfile2[512];

    CHAR *SFTpatt = NULL, *IFO = NULL, *outputBname = NULL;
    INT4 startGPS = 0, endGPS = 0, blocksRngMean = 21;
    REAL8 f_min = 0.0, f_max = 0.0, timebaseline = 0, persistSNRthresh = 3.0, auto_track;

    LALStringVector *line_freq = NULL;
    REAL8Vector *freq_vect = NULL;
    INT4 persistAvgSeconds = 0;
    INT4 persistAvgOpt;
    REAL8Vector *this_epoch_avg = NULL, *this_epoch_avg_wt = NULL, *new_epoch_avg = NULL, *new_epoch_wt = NULL;
    LIGOTimeGPSVector *epoch_gps_times = NULL;
    REAL8VectorSequence *epoch_avg = NULL;

    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&SFTpatt,       "SFTs",          STRING, 'p', REQUIRED, "SFT location/pattern" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&IFO,           "IFO",           STRING, 'I', REQUIRED, "Detector (e.g., H1)" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&startGPS,      "startGPS",      INT4,   's', REQUIRED, "Starting GPS time (SFT timestamps must be >= this time)" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&endGPS,        "endGPS",        INT4,   'e', REQUIRED, "Ending GPS time (SFT timestamps must be < this time)" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&f_min,         "fMin",          REAL8,  'f', REQUIRED, "Minimum frequency" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&f_max,         "fMax",          REAL8,  'F', REQUIRED, "Maximum frequency" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&blocksRngMean, "blocksRngMean", INT4,   'w', OPTIONAL, "Running Median window size") == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&outputBname,   "outputBname",   STRING, 'o', OPTIONAL, "Base name of output files" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&timebaseline,  "timeBaseline",  REAL8,  't', REQUIRED, "The time baseline of sfts") == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&persistAvgSeconds, "persistAvgSeconds", INT4, 'T', OPTIONAL, "Time baseline in seconds for averaging SFTs to measure the persistency, must be >= timeBaseline (cannot also specify --persistAveOption)") == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&persistAvgOpt,     "persistAvgOption",  INT4, 'E', OPTIONAL, "Choose one of 1 = day, 2 = week, or 3 = month averaging for measuring the persistency (cannot also specify --persistAvgSeconds)") == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&persistSNRthresh,  "persistSNRthresh",  REAL8, 'z', OPTIONAL, "SNR of lines for being present in the data") == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&line_freq,     "lineFreq",      STRINGVector,  0, OPTIONAL, "CSV list of line frequencies (e.g., --lineFreq=21.5,22.0). If set, then an output file with all GPS start times of SFTs with float values of number of standard deviations above the mean (>0 indicates above mean). Be careful that the CSV list of values are interpreted as floating point values") == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&auto_track,    "autoTrack",     REAL8,  'a', OPTIONAL, "If specified, also track any frequency whose persistency is >= this threshold within range [0,1]") == XLAL_SUCCESS, XLAL_EFUNC );

    BOOLEAN should_exit = 0;
    XLAL_CHECK_MAIN(XLALUserVarReadAllInput(&should_exit, argc, argv, lalPulsarVCSInfoList) == XLAL_SUCCESS, XLAL_EFUNC);
    if (should_exit)
    {
        return(1);
    }

    printf("Starting spec_avg_long...\n");

    XLAL_CHECK_MAIN( blocksRngMean % 2 == 1, XLAL_EINVAL, "Need to provide an odd value for blocksRngMean");
    XLAL_CHECK_MAIN( blocksRngMean > 2, XLAL_EINVAL, "Need to provide value larger than 2 blocksRngMean");
    INT4 nside = (blocksRngMean - 1) / 2;

    // Make some checks to be sure the persistency options for averaging are set correctly
    XLAL_CHECK_MAIN( !(XLALUserVarWasSet(&persistAvgSeconds) && XLALUserVarWasSet(&persistAvgOpt)), XLAL_EINVAL, "Provide only one of --persistAvgSeconds or --persistAvgOption" );
    if (XLALUserVarWasSet(&persistAvgSeconds)) {
	XLAL_CHECK_MAIN( persistAvgSeconds >= timebaseline, XLAL_EINVAL, "--persistAvgSeconds must be >= --timebaseline" );
    }
    if (XLALUserVarWasSet(&persistAvgOpt)) {
	XLAL_CHECK_MAIN( persistAvgOpt>0 && persistAvgOpt<4, XLAL_EINVAL, "--persistAvgOption can only take a value of 1, 2, or 3" );
    }
    if (!XLALUserVarWasSet(&persistAvgSeconds) && !XLALUserVarWasSet(&persistAvgOpt)) {
	persistAvgSeconds = timebaseline;
    }
    if (XLALUserVarWasSet(&auto_track)) {
	XLAL_CHECK_MAIN( auto_track>=0 && auto_track<=1, XLAL_EINVAL, "--autoTrack must be within range [0,1]" );
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

    // Ensure that some SFTs were found given the start and end time and IFO constraints
    XLAL_CHECK_MAIN( catalog->length > 0, XLAL_EFAILED, "No SFTs found, please examine start time, end time, frequency range, etc.");

    printf("Now have SFT catalog with %d catalog files\n", catalog->length);

    // Count the number of epochs [either the persistAvgSeconds, or the persistAvgOpt (day, week, or month)]
    XLAL_CHECK_MAIN( (epoch_gps_times = setup_epochs(catalog, persistAvgOpt, XLALUserVarWasSet(&persistAvgOpt), persistAvgSeconds)) != NULL, XLAL_EFUNC );

    /* Output files */
    if (XLALUserVarWasSet(&outputBname)) strcpy(outbase, outputBname);
    else snprintf(outbase, sizeof(outbase), "spec_%.2f_%.2f_%s_%d_%d", f_min, f_max, constraints.detector, startTime.gpsSeconds, endTime.gpsSeconds);
	
    snprintf(outfile0, sizeof(outfile0), "%s.txt", outbase);
    snprintf(outfile1, sizeof(outfile1), "%s_PWA.txt", outbase);
    snprintf(outfile2, sizeof(outfile2), "%s_line_times.csv", outbase);

    UINT4 epoch_index = 0;

    printf("Looping over SFTs to compute average spectra\n");
    for (UINT4 j = 0; j<catalog->length; j++)
    {
        fprintf(stderr,"Extracting SFT %d...\n", j);

	//Extract one SFT at a time from the catalog
	//we do this by using a catalog timeslice to get just the current SFT
	SFTVector *sft_vect = NULL;
	XLAL_CHECK_MAIN( (sft_vect = extract_one_sft(catalog, j, f_min, f_max)) != NULL, XLAL_EFUNC );

	//Make sure the SFTs are the same length as what we're expecting from user input
	XLAL_CHECK_MAIN( fabs(timebaseline*sft_vect->data->deltaF - 1.0) <= 10.*LAL_REAL8_EPS, XLAL_EINVAL, "Expected SFTs with length %f but got %f", timebaseline, 1/sft_vect->data->deltaF );

	//For the first time through the loop, we allocate some vectors
        if (j == 0)
	{
	    UINT4 numBins = sft_vect->data->data->length;
	    f0 = sft_vect->data->f0;
	    deltaF = sft_vect->data->deltaF;
	    printf("numBins=%d, f0=%f, deltaF=%f\n", numBins, f0, deltaF);

	    if (line_freq != NULL) {
		XLAL_CHECK_MAIN( validate_line_freq(&line_freq, f0, deltaF, numBins) == XLAL_SUCCESS, XLAL_EFUNC );
		XLAL_CHECK_MAIN( (freq_vect = line_freq_str2dbl(line_freq)) != NULL, XLAL_EFUNC );
	    }

	    XLAL_CHECK_MAIN( (timeavg = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (timeavgwt = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (sumweight = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (persistency = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    memset(persistency->data, 0, sizeof(REAL8)*persistency->length);

	    // Allocate vectors for this epoch average and weights as well as if we need to make a new epoch
	    XLAL_CHECK_MAIN( (this_epoch_avg = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (this_epoch_avg_wt = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (new_epoch_avg = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    memset(new_epoch_avg->data, 0, sizeof(REAL8)*new_epoch_avg->length);
	    XLAL_CHECK_MAIN( (new_epoch_wt = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (epoch_avg = XLALCreateREAL8VectorSequence(epoch_gps_times->length, numBins)) != NULL, XLAL_EFUNC );
	}

	//Loop over the SFT bins
	for (UINT4 i = 0; i < sft_vect->data->data->length; i++)
	{
	    REAL8 thispower = POWER(sft_vect->data[0].data->data[i]);
	    REAL8 thisavepower = 0.;
	    UINT4 count = 0;
	    for (INT4 ii = -nside; ii <= nside; ii++)
	    {
	        //Only add to the cumulative average power if the variables are in range
		if ((INT4)i + ii >= 0 && (INT4)i + ii<(INT4)sft_vect->data->data->length) {
		    thisavepower += POWER(sft_vect->data[0].data->data[i + ii]);
		    count++;
		}
	    }
	    thisavepower /= count;
	    REAL8 weight = 1. / thisavepower;

	    //For the first SFT, just assign the values to the vector, otherwise accumulate
	    if (j == 0) {
	        timeavg->data[i] = thispower;
		timeavgwt->data[i] = thispower*weight;
		sumweight->data[i] = weight;

		this_epoch_avg->data[i] = thispower*weight;
		this_epoch_avg_wt->data[i] = weight;
	    } else {
	        timeavg->data[i] += thispower;
		timeavgwt->data[i] += thispower*weight;
		sumweight->data[i] += weight;

		// Only accumulate in the epoch average if this SFT is within the current epoch
		// Otherwise put the values in the new epoch average and weight vector.
		if ((epoch_index<(epoch_gps_times->length-1) && XLALGPSCmp(&sft_vect->data->epoch, &epoch_gps_times->data[epoch_index])>=0 && XLALGPSCmp(&sft_vect->data->epoch, &epoch_gps_times->data[epoch_index+1])<0) || epoch_index==(epoch_gps_times->length-1)) {
		    this_epoch_avg->data[i] += thispower*weight;
		    this_epoch_avg_wt->data[i] += weight;
		} else {
		    new_epoch_avg->data[i] = thispower*weight;
		    new_epoch_wt->data[i] = weight;
		}
	    }
	} // end loop over this SFT frequency bins

	// If we've started putting new values into the new epoch average, then we should conclude the
	// last epoch and load the values from the new epoch into the current epoch
	if (new_epoch_avg->data[0] != 0.0) {
	    // Compute noise weighted power in this epoch average
	    for (UINT4 i = 0; i<this_epoch_avg->length; i++) {
		this_epoch_avg->data[i] *= 2.0 / this_epoch_avg_wt->data[i] / timebaseline;
	    }

	    // Copy to the vector sequence of epoch averages
	    memcpy(&(epoch_avg->data[epoch_index*this_epoch_avg->length]), this_epoch_avg->data, sizeof(REAL8)*this_epoch_avg->length);
	    epoch_index++;

	    // Copy over the new epoch data into this epoch's average
	    memcpy(this_epoch_avg->data, new_epoch_avg->data, sizeof(REAL8)*new_epoch_avg->length);
	    memcpy(this_epoch_avg_wt->data, new_epoch_wt->data, sizeof(REAL8)*new_epoch_wt->length);

	    // Set the new epoch data to be zero again
	    memset(new_epoch_avg->data, 0, sizeof(REAL8)*new_epoch_avg->length);
	}
	// This just repeats the above if we are in the very last SFT
	if (j == catalog->length-1) {
	    for (UINT4 i = 0; i<this_epoch_avg->length; i++) {
		this_epoch_avg->data[i] = 2.*this_epoch_avg->data[i] / this_epoch_avg_wt->data[i] / timebaseline;
	    }

	    memcpy(&(epoch_avg->data[epoch_index*this_epoch_avg->length]), this_epoch_avg->data, sizeof(REAL8)*this_epoch_avg->length);
	}

	// Destroys current SFT Vector
	XLALDestroySFTVector(sft_vect);
	sft_vect = NULL;
    } // end loop over all SFTs

    XLALDestroyREAL8Vector(this_epoch_avg_wt);
    XLALDestroyREAL8Vector(new_epoch_avg);
    XLALDestroyREAL8Vector(new_epoch_wt);

    // Allocate vectors for running mean and standard deviation of each chunk
    REAL8Vector *means = NULL, *stds = NULL;
    XLAL_CHECK_MAIN( (means = XLALCreateREAL8Vector(epoch_avg->vectorLength - 2*nside)) != NULL, XLAL_EFUNC );
    XLAL_CHECK_MAIN( (stds = XLALCreateREAL8Vector(epoch_avg->vectorLength - 2*nside)) != NULL, XLAL_EFUNC );

    //Loop over epochs for persistency calculation
    for (UINT4 j=0; j<epoch_gps_times->length; j++) {
	// Use the mean and standard deviation of data for this epoch
	memcpy(this_epoch_avg->data, &(epoch_avg->data[j*epoch_avg->vectorLength]), sizeof(REAL8)*epoch_avg->vectorLength);
	XLAL_CHECK_MAIN( rngmean(this_epoch_avg, means, blocksRngMean) == XLAL_SUCCESS, XLAL_EFUNC );
	XLAL_CHECK_MAIN( rngstd(this_epoch_avg, means, stds, blocksRngMean) == XLAL_SUCCESS, XLAL_EFUNC );

	// At the end points of the running mean, we need to re-use the end values
	// This is slightly sub-optimal
	for (UINT4 i=0; i<epoch_avg->vectorLength; i++) {
	    REAL8 mean, std;
	    XLAL_CHECK_MAIN( select_mean_std_from_vect(&mean, &std, means, stds, i, (UINT4)nside) == XLAL_SUCCESS, XLAL_EFUNC );
	    // Compare this SFT frequency data point with the running mean and standard deviation
	    if ((epoch_avg->data[j*epoch_avg->vectorLength + i] - mean)/std > persistSNRthresh) {
	        persistency->data[i] += 1.0;
	    }
	} // end loop over frequencies
    } // end loop over epochs

    // Normalize persistency to be in range 0 - 1 based on the number of epochs
    for (UINT4 i=0; i<persistency->length; i++) {
	persistency->data[i] /= (REAL8)epoch_gps_times->length;

	// if auto tracking, append any frequencies to the list
	if (XLALUserVarWasSet(&auto_track) && (persistency->data[i] >= auto_track)) {
	    if (freq_vect != NULL) {
		XLAL_CHECK_MAIN( (freq_vect = XLALResizeREAL8Vector(freq_vect, freq_vect->length + 1)) != NULL, XLAL_EFUNC );
	    } else {
		XLAL_CHECK_MAIN( (freq_vect = XLALCreateREAL8Vector(1)) != NULL, XLAL_EFUNC );
	    }
	    freq_vect->data[freq_vect->length-1] = f0 + deltaF*i;
	}
    } // end loop over frequency bins

    // allocate arrays for monitoring specific line frequencies
    REAL8VectorSequence *line_excess_in_epochs = NULL;
    if (freq_vect != NULL) {
	INT4Vector *line_freq_bin_in_sft = NULL;
        XLAL_CHECK_MAIN( (line_freq_bin_in_sft = XLALCreateINT4Vector(freq_vect->length)) != NULL, XLAL_EFUNC );
	XLAL_CHECK_MAIN( (line_excess_in_epochs = XLALCreateREAL8VectorSequence(freq_vect->length, epoch_gps_times->length)) != NULL, XLAL_EFUNC );
        memset(line_excess_in_epochs->data, 0, sizeof(REAL8)*line_excess_in_epochs->length*line_excess_in_epochs->vectorLength);

        // Set line_freq_array as the REAL8 values of the line_freq string vector
        // Set line_freq_bin_in_sft as the nearest INT4 values of the line frequencies to track
        for (UINT4 n=0; n<freq_vect->length; n++) {
	    line_freq_bin_in_sft->data[n] = (INT4)round((freq_vect->data[n]-f0)/deltaF);
        }

	// Loop over chunks checking if the line frequencies to track are above threshold in each chunk
	for (UINT4 j=0; j<epoch_gps_times->length; j++) {
	    // Use the standard deviation of data for this epoch
	    memcpy(this_epoch_avg->data, &(epoch_avg->data[j*epoch_avg->vectorLength]), sizeof(REAL8)*epoch_avg->vectorLength);
	    XLAL_CHECK_MAIN( rngmean(this_epoch_avg, means, blocksRngMean) == XLAL_SUCCESS, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( rngstd(this_epoch_avg, means, stds, blocksRngMean) == XLAL_SUCCESS, XLAL_EFUNC );

	    // Check each line frequency
	    for (UINT4 i=0; i<line_freq_bin_in_sft->length; i++) {
		REAL8 mean, std;
		XLAL_CHECK_MAIN( select_mean_std_from_vect(&mean, &std, means, stds, line_freq_bin_in_sft->data[i], (UINT4)nside) == XLAL_SUCCESS, XLAL_EFUNC );

		// assign the value of the number of standard devaiations above the mean for this chunk
		line_excess_in_epochs->data[i*line_excess_in_epochs->vectorLength + j] = (epoch_avg->data[j*epoch_avg->vectorLength + line_freq_bin_in_sft->data[i]] - mean)/std;
	    } // end loop over lines to track
	} // end loop over chunks of SFT averages

	XLALDestroyINT4Vector(line_freq_bin_in_sft);
    } // end if line tracking

    XLALDestroyREAL8Vector(means);
    XLALDestroyREAL8Vector(stds);
    XLALDestroyREAL8Vector(this_epoch_avg);
    XLALDestroyREAL8VectorSequence(epoch_avg);

    /* Print to files */
    SPECOUT = fopen(outfile0, "w");
    WTOUT = fopen(outfile1, "w");
    for (UINT4 i = 0; i < timeavg->length; i++)
    {
        REAL8 f = f0 + ((REAL4)i)*deltaF;
	REAL8 PSD = 2.*timeavg->data[i] / ((REAL4)catalog->length) / timebaseline;
	REAL8 PSDWT = 2.*timeavgwt->data[i] / sumweight->data[i] / timebaseline;
	REAL8 AMPPSD = pow(PSD, 0.5);
	REAL8 AMPPSDWT = pow(PSDWT, 0.5);
	REAL8 persist = persistency->data[i];
	fprintf(SPECOUT, "%16.8f %g %g %g %g %g\n", f, PSD, AMPPSD, PSDWT, AMPPSDWT, persist);

	REAL8 PWA_TAVGWT = timeavgwt->data[i];
	REAL8 PWA_SUMWT = sumweight->data[i];
	fprintf(WTOUT, "%16.8f %g %g\n", f, PWA_TAVGWT, PWA_SUMWT);
    }
    fclose(SPECOUT);
    fclose(WTOUT);

    // If user specified a list of line frequencies to monitor then print those data to a file
    if (freq_vect != NULL) {
        // open the line tracking output file for writing
        LINEOUT = fopen(outfile2, "w");
	if (!XLALUserVarWasSet(&persistAvgOpt)) {
	    fprintf(LINEOUT, "# GPS Epoch (len = %d s)", persistAvgSeconds);
	} else {
	    fprintf(LINEOUT, "# GPS Epoch (option = %d)", persistAvgOpt);
	}
	for (UINT4 n=0; n<freq_vect->length; n++) {
	    fprintf(LINEOUT,",%16.8f Hz", freq_vect->data[n]);
	}
	fprintf(LINEOUT, "\n");
	for (UINT4 i = 0; i < epoch_gps_times->length; i++) {
	    fprintf(LINEOUT, "%d", epoch_gps_times->data[i].gpsSeconds);
	    for (UINT4 n=0; n<freq_vect->length; n++) {
		fprintf(LINEOUT, ",%.4f", line_excess_in_epochs->data[n*line_excess_in_epochs->vectorLength + i]);
	    }
	    fprintf(LINEOUT, "\n");
	}
	fclose(LINEOUT);
	XLALDestroyREAL8VectorSequence(line_excess_in_epochs);
	XLALDestroyREAL8Vector(freq_vect);
    }

	/*------------------------------------------------------------------------------------------------------------------------*/

    fprintf(stderr,"Destroying Variables\n");

    XLALDestroyTimestampVector(epoch_gps_times);
    XLALDestroySFTCatalog(catalog);
    XLALDestroyREAL8Vector(timeavg);
    XLALDestroyREAL8Vector(timeavgwt);
    XLALDestroyREAL8Vector(sumweight);
    XLALDestroyREAL8Vector(persistency);

    XLALDestroyUserVars();
    fprintf(stderr,"Done Destroying Variables\n");
    fprintf(stderr, "end of spec_avg_long\n");

    return(0);

}
/* END main */


/* Set up epochs:
   This counts the number of epochs that would be set;
   allocates a LIGOTimeGPSVector for the number of epochs;
   populates the elements with the GPS times of the epochs */
LIGOTimeGPSVector * setup_epochs(const SFTCatalog *catalog, const INT4 persistAvgOpt, const BOOLEAN persistAvgOptWasSet, const INT4 persistAvgSeconds)
{
    UINT4 epoch_index = 0;
    LIGOTimeGPS XLAL_INIT_DECL(epoch);
    LIGOTimeGPS XLAL_INIT_DECL(test_epoch);
    struct tm epoch_utc;

    // First just figure out how many epochs there are
    for (UINT4 j=0; j<catalog->length; j++) {
	if (j==0) {
	    XLAL_CHECK_NULL( set_sft_avg_epoch(&epoch_utc, &epoch, catalog->data[j].header.epoch, persistAvgOpt, persistAvgOptWasSet) == XLAL_SUCCESS, XLAL_EFUNC );
	    epoch_index++;
	} else {
	    XLAL_CHECK_NULL( set_sft_avg_epoch(&epoch_utc, &test_epoch, catalog->data[j].header.epoch, persistAvgOpt, persistAvgOptWasSet) == XLAL_SUCCESS, XLAL_EFUNC );
	    if ((persistAvgOptWasSet && XLALGPSDiff(&test_epoch, &epoch)!=0) || (!persistAvgOptWasSet && XLALGPSDiff(&test_epoch, &epoch)>=persistAvgSeconds)) {
		epoch = test_epoch;
		epoch_index++;
	    }
	}
    }

    // Now allocate a timestamp vector
    LIGOTimeGPSVector *epoch_gps_times = NULL;
    XLAL_CHECK_NULL( (epoch_gps_times = XLALCreateTimestampVector(epoch_index)) != NULL, XLAL_EFUNC );
    epoch_index = 0;

    // Now assign the GPS times
    for (UINT4 j=0; j<catalog->length; j++) {
	if (j==0) {
	    XLAL_CHECK_NULL( set_sft_avg_epoch(&epoch_utc, &epoch_gps_times->data[epoch_index], catalog->data[j].header.epoch, persistAvgOpt, persistAvgOptWasSet) == XLAL_SUCCESS, XLAL_EFUNC );
	    epoch_index++;
	} else {
	    XLAL_CHECK_NULL( set_sft_avg_epoch(&epoch_utc, &test_epoch, catalog->data[j].header.epoch, persistAvgOpt, persistAvgOptWasSet) == XLAL_SUCCESS, XLAL_EFUNC );
	    if ((persistAvgOptWasSet && XLALGPSDiff(&test_epoch, &epoch_gps_times->data[epoch_index-1])!=0) || (!persistAvgOptWasSet && XLALGPSDiff(&test_epoch, &epoch_gps_times->data[epoch_index-1])>=persistAvgSeconds)) {
		XLAL_CHECK_NULL( set_sft_avg_epoch(&epoch_utc, &epoch_gps_times->data[epoch_index], catalog->data[j].header.epoch, persistAvgOpt, persistAvgOptWasSet) == XLAL_SUCCESS, XLAL_EFUNC );
		epoch_index++;
	    }
	}
    }

    return epoch_gps_times;
}


/* Set the GPS time of each epoch, depending on if the --persistAvgOpt was set.
   If it was set (meaning persistAvgOptWasSet is true and persistAvgOpt=[1,3]),
   then the utc and epoch_start values are determined based on the GPS time of
   first SFT epoch.
   persistAvgOpt = 1: utc and epoch_start are set to the most recent UTC midnight
   persistAvgOpt = 2: utc and epoch_start are set to the most recent UTC midnight Sunday
   persistAvgOpt = 3: utc and epoch_start are set to the most recent UTC midnight first day of the month
   otherwise, set utc and epoch_start to the first_sft_epoch */
int set_sft_avg_epoch(struct tm *utc, LIGOTimeGPS *epoch_start, const LIGOTimeGPS first_sft_epoch, const INT4 persistAvgOpt, const BOOLEAN persistAvgOptWasSet)
{
    //Figure out the UTC time of the SFT epoch
    XLAL_CHECK( (XLALGPSToUTC(utc, first_sft_epoch.gpsSeconds)) != NULL, XLAL_EFUNC );

    // If using the persistAvgOpt method 1, 2, or 3
    if (persistAvgOptWasSet && persistAvgOpt==1) {
	utc->tm_sec = utc->tm_min = utc->tm_hour = 0;
	epoch_start->gpsSeconds = XLALUTCToGPS(utc);
    } else if (persistAvgOptWasSet && persistAvgOpt==2) {
	utc->tm_sec = utc->tm_min = utc->tm_hour = 0;
	epoch_start->gpsSeconds = XLALUTCToGPS(utc) - (utc->tm_wday*24*3600);
    } else if (persistAvgOptWasSet && persistAvgOpt==3) {
	utc->tm_sec = utc->tm_min = utc->tm_hour = 0;
	utc->tm_mday = 1;
	XLAL_CHECK( (XLALFillUTC(utc)) != NULL, XLAL_EFUNC );
	epoch_start->gpsSeconds = XLALUTCToGPS(utc);
    } else {
	//Otherwise, the epoch is just the GPS time
	epoch_start->gpsSeconds = first_sft_epoch.gpsSeconds;
    }

    //We again set the UTC time from the GPS to make sure that the UTC and epoch time return are synchronized
    XLAL_CHECK( (XLALGPSToUTC(utc, epoch_start->gpsSeconds)) != NULL, XLAL_EFUNC );

    return XLAL_SUCCESS;
}


/* Extract a single SFT from an SFTCatalog: the SFT indicated by sft_index with band f_min to f_max*/
SFTVector * extract_one_sft(const SFTCatalog *full_catalog, const UINT4 sft_index, const REAL8 f_min, const REAL8 f_max)
{
    // Initialize an SFTCatalog
    SFTCatalog XLAL_INIT_DECL(catalogSlice);

    //Set start time
    //Set end time just 0.01 seconds after the start time. This is sufficiently small to get just one SFT
    LIGOTimeGPS thisSFTstarttime = full_catalog->data[sft_index].header.epoch;
    LIGOTimeGPS thisSFTendtime = thisSFTstarttime;
    XLAL_CHECK_NULL( XLALGPSAdd(&thisSFTendtime, 0.01) != NULL, XLAL_EFUNC );

    // Get the catalog of the single SFT from the full catalog
    XLAL_CHECK_NULL( XLALSFTCatalogTimeslice(&catalogSlice, full_catalog, &thisSFTstarttime, &thisSFTendtime) == XLAL_SUCCESS, XLAL_EFUNC );

    //Extract the SFT
    SFTVector *sft_vect = NULL;
    XLAL_CHECK_NULL( (sft_vect = XLALLoadSFTs(&catalogSlice, f_min, f_max)) != NULL, XLAL_EFUNC );

    //Check we got only one SFT; no more, no less
    XLAL_CHECK_NULL( sft_vect->length == 1, XLAL_EBADLEN, "Oops, got %d SFTs instead of one", sft_vect->length );

    return sft_vect;
}


/* Validate that the line frequency given is within the range of f_min to f_max.
   If there are any lines outside the range specified, the function will remove those lines from
   the list by making a new vector with the valid lines in range, destroying the old list, and
   returning a pointer to the new list */
int validate_line_freq(LALStringVector **line_freq, const REAL8 f0, const REAL8 deltaF, const UINT4 numBins)
{
    UINT4Vector *valid = NULL;
    XLAL_CHECK( (valid = XLALCreateUINT4Vector((*line_freq)->length)) != NULL, XLAL_EFUNC );
    UINT4 removeLength = 0;
    for (UINT4 i=0; i<(*line_freq)->length; i++) {
	REAL8 freq = atof((*line_freq)->data[i]);
	INT4 bin = (INT4)round((freq-f0)/deltaF);
	if (bin >=0 && bin < (INT4)numBins) {
	    valid->data[i] = 1;
	} else {
	    valid->data[i] = 0;
	    removeLength++;
	}
    }
    if (removeLength > 0) {
	LALStringVector *new_line_freq = NULL;
	for (UINT4 i=0; i<valid->length; i++) {
	    if (valid->data[i] == 1) {
		XLAL_CHECK( (new_line_freq = XLALAppendString2Vector(new_line_freq, (*line_freq)->data[i])) != NULL, XLAL_EFUNC );
	    } else {
		fprintf(stderr, "NOTE: requested frequency to monitor %s Hz is outside of requested band and will not be included in output\n", (*line_freq)->data[i]);
	    }
	}
	XLALDestroyStringVector(*line_freq);
	*line_freq = new_line_freq;
    }

    XLALDestroyUINT4Vector(valid);

    return XLAL_SUCCESS;
}


REAL8Vector * line_freq_str2dbl(const LALStringVector *line_freq)
{
    REAL8Vector *freq_vect = NULL;
    XLAL_CHECK_NULL( (freq_vect = XLALCreateREAL8Vector(line_freq->length)) != NULL, XLAL_EFUNC );
    for (UINT4 i=0; i<line_freq->length; i++) {
	freq_vect->data[i] = atof(line_freq->data[i]);
    }
    return freq_vect;
}


/* Compute the running mean of a REAL8Vector */
int rngmean(const REAL8Vector *input, const REAL8Vector *output, const INT4 blocksRngMean)
{
   REAL8Vector *win = NULL;
   XLAL_CHECK( (win = XLALCreateREAL8Vector(blocksRngMean)) != NULL, XLAL_EFUNC );
   memcpy(win->data, input->data, sizeof(REAL8)*blocksRngMean);
   output->data[0] = 0.0;
   for (UINT4 i = 0; i<win->length; i++) {
      output->data[0] += win->data[i];
   }
   output->data[0] /= (REAL8)blocksRngMean;
   for (UINT4 i = 1; i<output->length; i++) {
      output->data[i] = output->data[i-1] - win->data[0]/(REAL8)blocksRngMean;
      memmove(win->data, &win->data[1], sizeof(REAL8)*(blocksRngMean-1));
      win->data[win->length-1] = input->data[i+(blocksRngMean-1)];
      output->data[i] += win->data[win->length-1]/(REAL8)blocksRngMean;
   }
   XLALDestroyREAL8Vector(win);
   return XLAL_SUCCESS;
}


/* Compute the running standard deviation of a REAL8Vector using a vector of
   pre-computed mean values (use rngmean() above) */
int rngstd(const REAL8Vector *input, const REAL8Vector *means, const REAL8Vector *output, const INT4 blocksRngMean)
{
   REAL8Vector *win = NULL;
   XLAL_CHECK( (win = XLALCreateREAL8Vector(blocksRngMean)) != NULL, XLAL_EFUNC );
   for (UINT4 i = 0; i<output->length; i++) {
      memcpy(win->data, &input->data[i], sizeof(REAL8)*blocksRngMean);
      output->data[i] = 0.0;
      for (UINT4 j = 0; j<win->length; j++) {
	 output->data[i] += (win->data[j] - means->data[i]) * (win->data[j] - means->data[i]);
      }
   }
   for (UINT4 i=0; i<output->length; i++) {
      output->data[i] = sqrt(output->data[i]/(REAL8)(blocksRngMean-1));
   }
   XLALDestroyREAL8Vector(win);
   return XLAL_SUCCESS;
}


/* Select the specific mean and standard deviation from the running means and standard deviations.
   Need to know the size of a "side" of the block size, so that is passed as well as the index. */
int select_mean_std_from_vect(REAL8 *mean, REAL8 *std, const REAL8Vector *means, const REAL8Vector *stds, const UINT4 idx, const UINT4 nside)
{
    // At the end points, we need to re-use the 0-th or end values
    // This is slightly sub-optimal but is simple
    if (idx < nside) {
	*mean = means->data[0];
	*std = stds->data[0];
    } else if (idx >= means->length) {
	*mean = means->data[means->length - 1];
	*std = stds->data[stds->length - 1];
    } else {
	*mean = means->data[idx - nside];
	*std = stds->data[idx - nside];
    }

    return XLAL_SUCCESS;
}
