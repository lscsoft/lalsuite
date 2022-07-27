/*
*  Copyright (C) 2007 Gregory Mendell
*                2021 Evan Goetz
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
 * \ingroup lalpulsar_bin_SFTTools
 */

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LALPulsarVCSInfo.h>


int main(int argc, char **argv)
{
    FILE *fp  = NULL, *fp2 = NULL, *fp3 = NULL, *fp4 = NULL;
    
    SFTCatalog *catalog = NULL;
    SFTVector *sft_vect = NULL;
    UINT4 NumBinsAvg;
    SFTConstraints XLAL_INIT_DECL(constraints);
    LIGOTimeGPS XLAL_INIT_DECL(startTime);
    LIGOTimeGPS XLAL_INIT_DECL(endTime);
    REAL4Vector *timeavg = NULL;
    INT4 timebaseline = 0;
    CHAR outbase[256], outfile[512], outfile2[512], outfile3[512], outfile4[512];
    
    CHAR *SFTpatt = NULL, *IFO = NULL, *outputBname = NULL;
    INT4 startGPS = 0, endGPS = 0;
    REAL8 f_min = 0.0, f_max = 0.0, freqres = 0.0;
    INT4 blocksRngMed = 101, cur_epoch = 0;
    
    /* these varibales are for converting GPS seconds into UTC time and date*/
    struct tm date;
    INT4Vector *timestamps = NULL;

  /*========================================================================================================================*/
    
    
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&SFTpatt,      "SFTs",         STRING, 'p', REQUIRED, "SFT location/pattern" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&IFO,          "IFO",          STRING, 'I', REQUIRED, "Detector" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&startGPS,     "startGPS",     INT4,   's', REQUIRED, "Starting GPS time (SFT timestamps must be >= this)" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&endGPS,       "endGPS",       INT4,   'e', REQUIRED, "Ending GPS time (SFT timestamps must be < this)" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&f_min,        "fMin",         REAL8,  'f', REQUIRED, "Minimum frequency in Hz" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&f_max,        "fMax",         REAL8,  'F', REQUIRED, "Maximum frequency in Hz" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&blocksRngMed, "blocksRngMed", INT4,   'w', OPTIONAL, "Running Median window size") == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&outputBname,  "outputBname",  STRING, 'o', OPTIONAL, "Base name of output files" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&freqres,      "freqRes",      REAL8,  'r', REQUIRED, "Spectrogram freq resolution in Hz" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&timebaseline, "timeBaseline", INT4,  't', REQUIRED, "The time baseline of sfts in seconds") == XLAL_SUCCESS, XLAL_EFUNC);
    
    BOOLEAN should_exit = 0;
    XLAL_CHECK_MAIN( XLALUserVarReadAllInput(&should_exit, argc, argv, lalPulsarVCSInfoList) == XLAL_SUCCESS, XLAL_EFUNC );
    if (should_exit) {
       return(1);
    }

    // Populate the startTime and endTime LIGOTimeGPS variables
    startTime.gpsSeconds = startGPS;
    endTime.gpsSeconds = endGPS;
    startTime.gpsNanoSeconds = endTime.gpsNanoSeconds = 0;

    // Populate the SFT catalog constraints
    constraints.minStartTime = &startTime;
    constraints.maxStartTime = &endTime;
    constraints.detector = IFO;

    // Load SFT catalog
    XLAL_CHECK_MAIN( (catalog = XLALSFTdataFind ( SFTpatt, &constraints )) != NULL, XLAL_EFUNC );

    // Ensure that some SFTs were found given the start and end time and IFO constraints
    XLAL_CHECK_MAIN( catalog->length > 0, XLAL_EFAILED, "No SFTs found, please examine start time, end time, frequency range, etc.");

    // Read in data from the SFT catalog according to the minimum and maximum frequencies
    XLAL_CHECK_MAIN( (sft_vect = XLALLoadSFTs ( catalog, f_min, f_max)) != NULL, XLAL_EFUNC );

    // We no longer need the catalog so it can be destroyed at this point
    XLALDestroySFTCatalog(catalog);

    fprintf(stderr, "SFTs = %d\tSFT bins = %d\tf0 = %f\n", sft_vect->length, sft_vect->data->data->length, sft_vect->data->f0);

    if (XLALUserVarWasSet(&outputBname)) strcpy(outbase, outputBname);
    else snprintf(outbase, sizeof(outbase), "spec_%.2f_%.2f_%s_%d_%d", f_min, f_max, constraints.detector, startTime.gpsSeconds, endTime.gpsSeconds);
    
    snprintf(outfile, sizeof(outfile),  "%s", outbase);/*cg; name of first file to be output*/
    snprintf(outfile2, sizeof(outfile2), "%s_timestamps", outbase);/*cg: name of second file to be output*/
    /* snprintf(outfile3, sizeof(outfile3), "%s.txt", outbase); */ /*cg; name of third file to be output*/
    snprintf(outfile3, sizeof(outfile3), "%s_timeaverage", outbase); /* Use the old outfile3 name from this code; python code will output the .txt version of this file. */
    snprintf(outfile4, sizeof(outfile4), "%s_date", outbase);/*cg;file for outputting the date, which is used in matlab plotting.*/

    fp = fopen(outfile, "w");/*cg;  open all three files for writing, if they don't exist create them, if they do exist overwrite them*/
    fp2 = fopen(outfile2, "w");
    fp3 = fopen(outfile3, "w");
    fp4 = fopen(outfile4, "w");

/*----------------------------------------------------------------------------------------------------------------*/

    // Compute the number of bins in an average; at minimum this should be 1 bin
    // It produces the same frequency resolution as specified in the arguments passed to
    // fscanDriver.py
    REAL8 spectrogram_blocksize = round(freqres * sft_vect->data[0].data->length / (f_max-f_min));
    if (spectrogram_blocksize < 1.0) NumBinsAvg = 1;
    else NumBinsAvg = (UINT4)spectrogram_blocksize;

    // Record timestamps for each SFT or gap
    // This is not known a priori so we end up resizing this vector as we go
    XLAL_CHECK_MAIN( (timestamps = XLALCreateINT4Vector(0)) != NULL, XLAL_EFUNC );

    // Loop over the SFT vector, SFT by SFT
    for (UINT4 j=0; j<sft_vect->length; j++)
    {
        // GPS time of the current SFT and print to file
        cur_epoch = sft_vect->data[j].epoch.gpsSeconds;
        fprintf(fp2, "%d.\t%d\n", timestamps->length, sft_vect->data[j].epoch.gpsSeconds);

	// Get the current UTC time from the GPS seconds of the current SFT and print to file
        XLAL_CHECK_MAIN( XLALGPSToUTC(&date, sft_vect->data[j].epoch.gpsSeconds) != NULL, XLAL_EFUNC );
        fprintf(fp4, "%d\t %i\t %i\t %i\t %i\t %i\t %i\n", timestamps->length, (date.tm_year+1900), date.tm_mon+1, date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec);

	// Here is where the timestamps vector is resized
        XLAL_CHECK_MAIN( XLALResizeINT4Vector(timestamps, timestamps->length + 1) != NULL, XLAL_EFUNC );
        timestamps->data[timestamps->length-1] = sft_vect->data[j].epoch.gpsSeconds;

	// Loop over the number of bins in each SFT, jumping by the average number
	// Start at the highest bin index and average down. This avoids running off the end of the SFT
        for (UINT4 i=NumBinsAvg-1; i<sft_vect->data[j].data->length; i+=NumBinsAvg)
        {
            REAL8 avg = 0.0;
	    for (UINT4 k=0; k<NumBinsAvg; k++) {
              const REAL8 re = (REAL8)crealf(sft_vect->data[j].data->data[i-k]);
              const REAL8 im = (REAL8)cimagf(sft_vect->data[j].data->data[i-k]);
              avg += 2.0*(re*re + im*im)/(REAL8)timebaseline;
	    }
            fprintf(fp, "%e\t", sqrt(avg / (REAL8)NumBinsAvg)); /* 06/15/2017 gam; then take sqrt here. */
        }
        fprintf(fp,"\n");

	// Fill in gaps where there is no SFT data with zeros
        if ( j < (sft_vect->length-1) )/*in all cases except when we are examining the last sft, check that there is no gap to the next sft*/
        {
            /*test to see if the next SFT immediately follows, if not entries in the matrix until there is one*/
            while (cur_epoch+timebaseline < sft_vect->data[j+1].epoch.gpsSeconds)
            {
                cur_epoch += timebaseline;
                fprintf(fp2, "%d.\t%d\n", timestamps->length, cur_epoch );

                XLAL_CHECK_MAIN( XLALGPSToUTC(&date, cur_epoch) != NULL, XLAL_EFUNC );
                fprintf(fp4, "%d\t %i\t %i\t %i\t %i\t %i\t %i\n", timestamps->length, (date.tm_year+1900), date.tm_mon+1, date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec);
                XLAL_CHECK_MAIN( XLALResizeINT4Vector(timestamps, timestamps->length + 1) != NULL, XLAL_EFUNC );
                timestamps->data[timestamps->length-1] = cur_epoch;

		for (UINT4 i=NumBinsAvg-1; i<sft_vect->data[j].data->length; i+=NumBinsAvg)
                {
                    REAL8 avg = 0.0;
                    fprintf(fp, "%e\t", avg);
                }
                fprintf(fp,"\n");
            }
            
        }
    }
    fprintf(stderr,"finished checking for missing sfts, l=%d\n", timestamps->length);

    // Normalize the SFTs in the vector
    XLAL_CHECK_MAIN( XLALNormalizeSFTVect(sft_vect, blocksRngMed, 0.0) == XLAL_SUCCESS, XLAL_EFUNC );

    // Allocate the vector for the normalized data to be averaged together
    XLAL_CHECK_MAIN( (timeavg = XLALCreateREAL4Vector(sft_vect->data->data->length)) != NULL, XLAL_EFUNC );

    // Loop over the SFTs
    for (UINT4 j=0; j<sft_vect->length; j++)
    {
        // Loop over the frequency bins in the SFT
        for (UINT4 i=0; i<sft_vect->data[j].data->length; i++)
        {
            const REAL8 re = (REAL8)crealf(sft_vect->data[j].data->data[i]);
            const REAL8 im = (REAL8)cimagf(sft_vect->data[j].data->data[i]);
            if (j == 0) 
            {
                timeavg->data[i] = re*re + im*im;
            } 
            else 
            {
                timeavg->data[i] += re*re + im*im;
            }
        }
    }
    
    // Write the averaged data to a file
    for (UINT4 i=0; i<sft_vect->data[0].data->length; i++)
    {
        REAL8 f = sft_vect->data->f0 + ((REAL4)i)*sft_vect->data->deltaF;
	REAL4 PWR = timeavg->data[i]/((REAL4)sft_vect->length);
        fprintf(fp3, "%16.6f %16.3f \n", f, PWR);
    }
/*------------------------------------------------------------------------------------------------------------------------*/
/*End of normal spec_avg code.*/

    XLALDestroySFTVector ( sft_vect );
    XLALDestroyREAL4Vector(timeavg);
    XLALDestroyINT4Vector(timestamps);

    XLALDestroyUserVars();

    /*close all the files, spec_avg.c is done, all info written to the files.*/
    fclose(fp);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);

    fprintf(stderr, "end of spec_avg\n");

    return(0);


}
/* END main */
