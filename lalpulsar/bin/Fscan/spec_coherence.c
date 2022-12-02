/*
*  Copyright (C) 2021 Evan Goetz
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

#include <lal/SFTfileIO.h>
#include <lal/LALStdio.h>
#include <lal/UserInput.h>
#include <lal/LALPulsarVCSInfo.h>


int main(int argc, char **argv)
{
    FILE *COHOUT = NULL;

    SFTCatalog *catalog_a = NULL, *catalog_b = NULL;
    SFTVector *sft_vect_a = NULL, *sft_vect_b = NULL;
    SFTConstraints XLAL_INIT_DECL(constraints);
    LIGOTimeGPS startTime, endTime;
    REAL8Vector *psd_a = NULL, *psd_b = NULL;
    COMPLEX16Vector *coh = NULL;
    REAL8 f0 = 0, deltaF = 0;
    CHAR outbase[256], outfile0[512];

    CHAR *SFTpattA = NULL, *SFTpattB = NULL, *outputBname = NULL;
    INT4 startGPS = 0, endGPS = 0;
    REAL8 f_min = 0.0, f_max = 0.0, timebaseline = 0;

    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&SFTpattA,      "ChASFTs",         STRING, 'p', REQUIRED, "SFT location/pattern. Possibilities are:\n"
                                           " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&SFTpattB,      "ChBSFTs",         STRING, 'q', REQUIRED, "SFT location/pattern. Possibilities are:\n"
                                           " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&startGPS,     "startGPS",     INT4,   's', REQUIRED, "Starting GPS time (SFT timestamps must be >= this)" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&endGPS,       "endGPS",       INT4,   'e', REQUIRED, "Ending GPS time (SFT timestamps must be < this)" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&f_min,        "fMin",         REAL8,  'f', REQUIRED, "Minimum frequency" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&f_max,        "fMax",         REAL8,  'F', REQUIRED, "Maximum frequency" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&outputBname,  "outputBname",  STRING, 'o', OPTIONAL, "Base name of output files" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&timebaseline, "timeBaseline", REAL8,  't', REQUIRED, "The time baseline of sfts") == XLAL_SUCCESS, XLAL_EFUNC);

    BOOLEAN should_exit = 0;
    XLAL_CHECK_MAIN(XLALUserVarReadAllInput(&should_exit, argc, argv, lalPulsarVCSInfoList) == XLAL_SUCCESS, XLAL_EFUNC);
    if (should_exit)
    {
        return(1);
    }

    printf("Starting spec_coherence...\n");

    /* Provide the constraints to the catalog */
    startTime.gpsSeconds = startGPS;
    startTime.gpsNanoSeconds = 0;
    constraints.minStartTime = &startTime;
    endTime.gpsSeconds = endGPS;
    endTime.gpsNanoSeconds = 0;
    constraints.maxStartTime = &endTime;

    printf("Calling XLALSFTdataFind with SFTpattA=%s\n", SFTpattA);
    XLAL_CHECK_MAIN( (catalog_a = XLALSFTdataFind ( SFTpattA, &constraints )) != NULL, XLAL_EFUNC );
    printf("Calling XLALSFTdataFind with SFTpattB=%s\n", SFTpattB);
    XLAL_CHECK_MAIN( (catalog_b = XLALSFTdataFind ( SFTpattB, &constraints )) != NULL, XLAL_EFUNC );

    /* Ensure that some SFTs were found given the start and end time and IFO constraints */
    XLAL_CHECK_MAIN( catalog_a->length > 0, XLAL_EFAILED, "No SFTs found for Ch A, please examine start time, end time, frequency range, etc.");
    XLAL_CHECK_MAIN( catalog_b->length > 0, XLAL_EFAILED, "No SFTs found for Ch B, please examine start time, end time, frequency range, etc.");

    printf("Now have Ch A SFT catalog with %d catalog files\n", catalog_a->length);
    printf("Now have Ch B SFT catalog with %d catalog files\n", catalog_b->length);

    if (XLALUserVarWasSet(&outputBname)) strcpy(outbase, outputBname);
    else snprintf(outbase, sizeof(outbase), "spec_%.2f_%.2f_%d_%d_coh", f_min, f_max, startTime.gpsSeconds, endTime.gpsSeconds);
	
    snprintf(outfile0, sizeof(outfile0), "%s.txt", outbase);

    COHOUT = fopen(outfile0, "w");

    UINT4 nAve = 0;
    
    printf("Looping over SFTs to compute coherence\n");
    for (UINT4 j = 0; j<catalog_a->length; j++)
    {
        /* Extract one SFT at a time from the catalog
	   we do this by using a catalog timeslice to get just the current SFT */
	SFTCatalog XLAL_INIT_DECL(catalogSlice_a);
	LIGOTimeGPS thisSFTstarttime = catalog_a->data[j].header.epoch;
	LIGOTimeGPS thisSFTendtime = thisSFTstarttime;
	XLAL_CHECK_MAIN( XLALGPSAdd(&thisSFTendtime, 0.01) != NULL, XLAL_EFUNC );
	XLAL_CHECK_MAIN( XLALSFTCatalogTimeslice(&catalogSlice_a, catalog_a, &thisSFTstarttime, &thisSFTendtime) == XLAL_SUCCESS, XLAL_EFUNC );

	SFTCatalog XLAL_INIT_DECL(catalogSlice_b);
	XLAL_CHECK_MAIN( XLALSFTCatalogTimeslice(&catalogSlice_b, catalog_b, &thisSFTstarttime, &thisSFTendtime) == XLAL_SUCCESS, XLAL_EFUNC );

	/* If no SFT from the B list was found, then just continue with the next SFT in the A list */
	if (catalogSlice_b.length == 0) {
	    continue;
	}
	
	fprintf(stderr,"Extracting SFT %d...\n", j);
	XLAL_CHECK_MAIN( (sft_vect_a = XLALLoadSFTs(&catalogSlice_a, f_min, f_max)) != NULL, XLAL_EFUNC );
	XLAL_CHECK_MAIN( sft_vect_a->length == 1, XLAL_EBADLEN, "Oops, got %d SFTs instead of one", sft_vect_a->length );

	XLAL_CHECK_MAIN( (sft_vect_b = XLALLoadSFTs(&catalogSlice_b, f_min, f_max)) != NULL, XLAL_EFUNC );
	XLAL_CHECK_MAIN( sft_vect_b->length == 1, XLAL_EBADLEN, "Oops, got %d SFTs instead of one", sft_vect_b->length );

	XLAL_CHECK_MAIN( sft_vect_a->data[0].deltaF * timebaseline == 1.0, XLAL_EINVAL, "Time baseline of SFTs and the request do not match" );
	XLAL_CHECK_MAIN( sft_vect_b->data[0].deltaF * timebaseline == 1.0, XLAL_EINVAL, "Time baseline of SFTs and the request do not match" );

	/* For the first time through the loop, we allocate some vectors */
        if (nAve == 0)
	{
	    UINT4 numBins = sft_vect_a->data->data->length;
	    f0 = sft_vect_a->data->f0;
	    deltaF = sft_vect_a->data->deltaF;

	    XLAL_CHECK_MAIN( (coh = XLALCreateCOMPLEX16Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (psd_a = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	    XLAL_CHECK_MAIN( (psd_b = XLALCreateREAL8Vector(numBins)) != NULL, XLAL_EFUNC );
	}

	/* Loop over the SFT bins computing cross spectrum (AB) and auto spectrum (AA and BB) */
	if (nAve == 0) {
	    for (UINT4 i = 0; i < sft_vect_a->data->data->length; i++)
	    {
	        coh->data[i] = sft_vect_a->data[0].data->data[i] * conj(sft_vect_b->data[0].data->data[i]);
		psd_a->data[i] = sft_vect_a->data[0].data->data[i] * conj(sft_vect_a->data[0].data->data[i]);
		psd_b->data[i] = sft_vect_b->data[0].data->data[i] * conj(sft_vect_b->data[0].data->data[i]);
	    }
	}
	else
	{
	    for (UINT4 i = 0; i < sft_vect_a->data->data->length; i++)
	    {
	        coh->data[i] += sft_vect_a->data[0].data->data[i] * conj(sft_vect_b->data[0].data->data[i]);
		psd_a->data[i] += sft_vect_a->data[0].data->data[i] * conj(sft_vect_a->data[0].data->data[i]);
		psd_b->data[i] += sft_vect_b->data[0].data->data[i] * conj(sft_vect_b->data[0].data->data[i]);
	    }
	}
	/* Destroys current SFT Vectors */
	XLALDestroySFTVector(sft_vect_a);
	XLALDestroySFTVector(sft_vect_b);
	sft_vect_a = sft_vect_b = NULL;

	nAve++;
    }

    XLAL_CHECK_MAIN( nAve > 0, XLAL_EFAILED, "No SFTs were found to be matching" );

    /* compute the final coherence (|AB|**2 / (AA * BB)) and print to file */
    for (UINT4 i = 0; i < coh->length; i++)
    {
        REAL8 f = f0 + ((REAL4)i)*deltaF;
	REAL8 COH = coh->data[i] * conj(coh->data[i]) / (psd_a->data[i] * psd_b->data[i]);
	fprintf(COHOUT, "%16.8f %g\n", f, COH);
    }

    fprintf(stderr,"Destroying Variables\n");
    XLALDestroySFTCatalog(catalog_a);
    XLALDestroySFTCatalog(catalog_b);

    XLALDestroyCOMPLEX16Vector(coh);
    XLALDestroyREAL8Vector(psd_a);
    XLALDestroyREAL8Vector(psd_b);

    fprintf(stderr,"Closing Files\n");
    fclose(COHOUT);

    XLALDestroyUserVars();
    fprintf(stderr,"Done Destroying Variables\n");
    fprintf(stderr, "end of spec_coherence\n");
    fprintf(stderr, "Spec_coherence_done!\n");

    return(0);

}
/* END main */


