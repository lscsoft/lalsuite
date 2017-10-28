/*
*  Copyright (C) 2007 Gregory Mendell
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
* \file
* \ingroup lalapps_pulsar_fscan
*/

//#define LAL_USE_OLD_COMPLEX_STRUCTS

/*temporary rubbish bin for headers*/
/*These are included in HeterodyneCrabPulsar files
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/BinaryPulsarTiming.h>*/
/*end of temporary rubbish bin*/

/*LAL header files*/
#include <lalapps.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdio.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LALMalloc.h>

/*normal c header files*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <lal/Date.h>/*cg; needed to use lal routine GPStoUTC, which is used to convert GPS seconds into UTC date*/


/*---------- DEFINES ----------*/
#define NUM 1000 

#define MIN_SFT_VERSION 1
#define MAX_SFT_VERSION 2

#define TRUE    1
#define FALSE   0

#define POLY64 0xd800000000000000ULL
#define TABLELEN 256

/*used for defining structures such as crabOutput*/

/*---------- Global variables ----------*/
static REAL8 fudge_up = 1 + 10 * LAL_REAL8_EPS;	// about ~1 + 2e-15
static REAL8 fudge_down = 1 - 10 * LAL_REAL8_EPS;	// about ~1 - 2e-15

/*----- Macros ----- */

#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define GPSEQUAL(gps1,gps2) (((gps1).gpsSeconds == (gps2).gpsSeconds) && ((gps1).gpsNanoSeconds == (gps2).gpsNanoSeconds))

#define GPSZERO(gps) (((gps).gpsSeconds == 0) && ((gps).gpsNanoSeconds == 0))

/*---------- internal types ----------*/

/* NOTE: the locator is implemented as an OPAQUE type in order to enforce encapsulation
* of the actual physical storage of SFTs and to ease future extensions of the interface.
* DO NOT TRY TO USE THIS TYPE OUTSIDE OF THIS MODULE!!
*/
struct tagSFTLocator
{
	CHAR *fname;		/* name of file containing this SFT */
	long offset;		/* SFT-offset with respect to a merged-SFT */
	UINT4 isft;           /* index of SFT this locator belongs to, used only in XLALLoadSFTs() */
};

typedef struct
{
	REAL8 version;
	INT4 gps_sec;
	INT4 gps_nsec;
	REAL8 tbase;
	INT4 first_frequency_index;
	INT4 nsamples;
} _SFT_header_v1_t;

typedef struct
{
	REAL8 version;
	INT4 gps_sec;
	INT4 gps_nsec;
	REAL8 tbase;
	INT4 first_frequency_index;
	INT4 nsamples;
	UINT8 crc64;
	CHAR detector[2];
	CHAR padding[2];
	INT4 comment_length;
} _SFT_header_v2_t;

/** segments read so far from one SFT */
typedef struct {
	UINT4 first;                     /**< first bin in this segment */
	UINT4 last;                      /**< last bin in this segment */
	LIGOTimeGPS epoch;               /**< timestamp of this SFT */
	struct tagSFTLocator *lastfrom;  /**< last bin read from this locator */
} SFTReadSegment;

///////////EXTRA FXN DEFS//////////////

static BOOLEAN is_valid_detector(const char *channel);
static int compareSFTloc(const void *ptr1, const void *ptr2);
static UINT4 read_sft_bins_from_fp(SFTtype *ret, UINT4 *firstBinRead, UINT4 firstBin2read, UINT4 lastBin2read, FILE *fp);
static int read_sft_header_from_fp(FILE *fp, SFTtype  *header, UINT4 *version, UINT8 *crc64, BOOLEAN *swapEndian, CHAR **SFTcomment, UINT4 *numBins);
static int read_v2_header_from_fp(FILE *fp, SFTtype *header, UINT4 *nsamples, UINT8 *header_crc64, UINT8 *ref_crc64, CHAR **SFTcomment, BOOLEAN swapEndian);
static int read_v1_header_from_fp(FILE *fp, SFTtype *header, UINT4 *nsamples, BOOLEAN swapEndian);
static void endian_swap(CHAR * pdata, size_t dsize, size_t nelements);
static int read_SFTversion_from_fp(UINT4 *version, BOOLEAN *need_swap, FILE *fp);
static UINT8 calc_crc64(const CHAR *data, UINT4 length, UINT8 crc);

SFTVector*
LALExtractSFT(const SFTCatalog *catalog,   /**< The 'catalogue' of SFTs to load */
REAL8 fMin,		   /**< minumum requested frequency (-1 = read from lowest) */
REAL8 fMax,		   /**< maximum requested frequency (-1 = read up to highest) */
REAL8 SFTindex
);


int main(int argc, char **argv)
{
	FILE *fp = NULL;
	FILE *fp2 = NULL;
	FILE *fp3 = NULL;
	FILE *fp4 = NULL;
	FILE *fp5 = NULL;
	LALStatus status = blank_status;

	SFTCatalog *catalog = NULL;
	SFTVector *sft_vect = NULL;
	INT4 i, j, ii, nside, count;
	//    INT4 l,k;
	INT4 numBins = 0, nSFT = 0, nSFTcheck = 0;
	SFTConstraints XLAL_INIT_DECL(constraints);
	LIGOTimeGPS startTime, endTime;
	//    REAL8 avg =0;
	REAL8 *timeavg = NULL, *timeavgwt = NULL, *sumweight = NULL;
	REAL8 PSD, AMPPSD, PSDWT, AMPPSDWT, weight, thispower, thisavepower, scalefactor;
	REAL8 PWA_TAVGWT, PWA_SUMWT;
	REAL8 f = 0, f0 = 0, deltaF = 0;
	CHAR outbase[256], outfile[256], outfile2[256], outfile3[256], outfile4[256] , outfile5[256];
	//    REAL8 NumBinsAvg =0;
	REAL8 timebaseline = 0;

	//BOOLEAN help = 0;
	CHAR *SFTpatt = NULL;
	CHAR *IFO = NULL;
	INT4 startGPS = 0;
	INT4 endGPS = 0;
	REAL8 f_min = 0.0;
	REAL8 f_max = 0.0;
	REAL8 freqres = 0.0;
	INT4 blocksRngMed = 101;
	CHAR *outputBname = NULL;
	//    INT4 cur_epoch = 0, next_epoch = 0;

	/* these varibales are for converting GPS seconds into UTC time and date*/
	//    LALUnixDate       date;
	//    CHARVector        *timestamp = NULL;
	CHARVector	     *year_date = NULL;
	//    REAL8Vector     *timestamps=NULL;

	CHAR *psrInput = NULL;
	CHAR *psrEphemeris = NULL;
	CHAR *earthFile = NULL;
	CHAR *sunFile = NULL;
	/*===============================================================================================================/home/pulsar/public_html/fscan/H2_OneArm/H2_OneArm_SUS/H2_OneArm_SUS/fscans_2012_07_16_23_03_00_PDT_Mon/H2_SUS-ETMY_M0_LOCK_L_IN1_DQ/sfts/tmp=========*/

	printf("Starting spec_avg...\n");

	//lalDebugLevel = 0;
	//LAL_CALL(LALGetDebugLevel(&status, argc, argv, 'v'), &status);


    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&SFTpatt,      "SFTs",         STRING, 'p', REQUIRED, "SFT location/pattern" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&IFO,          "IFO",          STRING, 'I', REQUIRED, "Detector" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&startGPS,     "startGPS",     INT4,   's', REQUIRED, "Starting GPS time" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&endGPS,       "endGPS",       INT4,   'e', REQUIRED, "Ending GPS time" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&f_min,        "fMin",         REAL8,  'f', REQUIRED, "Minimum frequency" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&f_max,        "fMax",         REAL8,  'F', REQUIRED, "Maximum frequency" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&blocksRngMed, "blocksRngMed", INT4,   'w', OPTIONAL, "Running Median window size") == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&outputBname,  "outputBname",  STRING, 'o', OPTIONAL, "Base name of output files" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&freqres,      "freqRes",      REAL8,  'r', REQUIRED, "Spectrogram freq resolution" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&timebaseline, "timeBaseline", REAL8,  't', REQUIRED, "The time baseline of sfts") == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&psrInput,     "psrInput",     STRING, 'P', OPTIONAL, "name of tempo pulsar file" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&psrEphemeris, "psrEphemeris", STRING, 'S', OPTIONAL, "pulsar ephemeris file" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&earthFile,    "earthFile",    STRING, 'y', OPTIONAL, "earth .dat file" ) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN( XLALRegisterNamedUvar(&sunFile,      "sunFile",      STRING, 'z', OPTIONAL, "sun .dat file" ) == XLAL_SUCCESS, XLAL_EFUNC);

    BOOLEAN should_exit = 0;
    XLAL_CHECK_MAIN(XLALUserVarReadAllInput(&should_exit, argc, argv, lalAppsVCSInfoList) == XLAL_SUCCESS, XLAL_EFUNC);
    if (should_exit)
      return(1);



	startTime.gpsSeconds = startGPS;/*cg; startTime is a structure, and gpsSeconds is a member of that structure*/
	startTime.gpsNanoSeconds = 0;/*cg; gps NanoSeconds is also a member of the startTime structure */
	//constraints.startTime = &startTime; /*cg; & operator gets the address of variable, &a is a pointer to a.  This line puts the startTime structure into the structure constraints*/
	constraints.minStartTime = &startTime; /*cg; & operator gets the address of variable, &a is a pointer to a.  This line puts the startTime structure into the structure constraints*/

	endTime.gpsSeconds = endGPS;
	endTime.gpsNanoSeconds = 0;
	//constraints.endTime = &endTime;/*cg; This line puts the end time into the structure constraints*/
	constraints.maxStartTime = &endTime;/*cg; This line puts the end time into the structure constraints*/
	constraints.detector = IFO;/*cg; this adds the interferometer into the contraints structure*/
	printf("Calling LALSFTdataFind with SFTpatt=%s\n", SFTpatt);

        catalog = XLALSFTdataFind ( SFTpatt, &constraints );/*cg; creates SFT catalog, uses the constraints structure*/

	printf("Now have SFT catalog with %d catalog files\n", catalog->length);

	if (catalog == NULL)/*need to check for a NULL pointer, and print info about circumstances if it is null*/
	{
		fprintf(stderr, "SFT catalog pointer is NULL!  There has been an error with LALSFTdataFind\n");
		fprintf(stderr, "LALStatus info.... status code: %d, message: %s, offending function: %s\n", status.statusCode, status.statusDescription, status.function);
		exit(0);
	}
	nSFT = catalog->length;
	if (nSFT == 0)
	{
		fprintf(stderr, "No SFTs found, please examine start time, end time, frequency range etc\n");
		exit(0);
	}

	//    printf("Loading SFTs...\n");
	//    LALLoadSFTs ( &status, &sft_vect, catalog, f_min,f_max);/*cg;reads the SFT data into the structure sft_vect*/
	//   printf("Loaded SFTs\n");

	//    fprintf(stderr, "nSFT = %d\tnumBins = %d\tf0 = %f\n", nSFT, numBins,sft_vect->data->f0);/*print->logs/spectrumAverage_testcg_0.err */
	if (XLALUserVarWasSet(&outputBname))
		strcpy(outbase, outputBname);
	else
		sprintf(outbase, "spec_%.2f_%.2f_%s_%d_%d", f_min, f_max, constraints.detector, startTime.gpsSeconds, endTime.gpsSeconds);/*cg; this is the default name for producing the output files, the different suffixes are just added to this*/
	sprintf(outfile, "%s", outbase);/*cg; name of first file to be output*/
	sprintf(outfile2, "%s_timestamps", outbase);/*cg: name of second file to be output*/
	sprintf(outfile3, "%s.txt", outbase);/*cg; name of third file to be output*/
	sprintf(outfile4, "%s_date", outbase);/*cg;file for outputting the date, which is used in matlab plotting.*/
	// ADDED FOR SPEC_AVG_LONG for use in cumulative plotting
	sprintf(outfile5, "%s_PWA.txt", outbase);

	fp = fopen(outfile, "w");/*cg;  open all three files for writing, if they don't exist create them, if they do exist overwrite them*/
	fp2 = fopen(outfile2, "w");
	fp3 = fopen(outfile3, "w");
	fp4 = fopen(outfile4, "w");
	fp5 = fopen(outfile5, "w");

	LALCHARCreateVector(&status, &year_date, (UINT4)128);


	/*----------------------------------------------------------------------------------------------------------------*/
	/*cg;  Create the third and final file, called   blah_blah_blah.txt.  This file will contain the data used in the matlab plot script to plot the normalised average power vs the frequency.*/

	/* Find time average of normalized SFTs */
	/*    LALNormalizeSFTVect(&status, sft_vect, blocksRngMed);
	LALNormalizeSFTVect(&status, sft_vect, blocksRngMed);   */

	/*    scalefactor = 1.e21; */
	scalefactor = 1.;

	printf("Looping over SFTs to compute average spectra\n");
	for (j = 0; j<nSFT; j++)

	{
		//subcatalog->length = 1;
		//subcatalog->data[0] = catalog->data[j];
		printf("Extracting SFT %d...\n", j);
		fprintf(stderr,"Extracting SFT %d...\n", j);
		//LALLoadSFTs(&status, &sft_vect, subcatalog, f_min, f_max);
		//sft_vect = XLALLoadSFTs(catalog,f_min,f_max,j);
		sft_vect = LALExtractSFT(catalog, f_min, f_max, j);/* reads the SFT data into the structure sft_vect*/
		//	printf("Extracted SFT %d\n",j);

		if (sft_vect == NULL)
		{
			fprintf(stderr, "SFT vector pointer is NULL for file %d!  There has been an error with LALLoadSFTs\n", j);
			exit(0);
		}

		if (j == 0)
		{
			numBins = sft_vect->data->data->length;/*the number of bins in the freq_range*/
			f0 = sft_vect->data->f0;
			deltaF = sft_vect->data->deltaF;
			printf("numBins=%d, f0=%f, deltaF=%f\n", numBins, f0, deltaF);
			//fprintf(stderr, "numBins=%d, f0=%f, deltaF=%f\n", numBins, f0, deltaF);
			timeavg = XLALMalloc(numBins*sizeof(REAL8));
			if (timeavg == NULL) fprintf(stderr, "Timeavg memory not allocated\n");
			timeavgwt = XLALMalloc(numBins*sizeof(REAL8));
			if (timeavgwt == NULL) fprintf(stderr, "Timeavgwt memory not allocated\n");
			sumweight = XLALMalloc(numBins*sizeof(REAL8));
			if (sumweight == NULL) fprintf(stderr, "Sumweight memory not allocated\n");
		}

		nSFTcheck = sft_vect->length;/* the number of sfts.*/
		if (nSFTcheck != 1)
		{
			fprintf(stderr, "Oops, nSFTcheck=%d instead of one\n", nSFTcheck);
			exit(0);
		}

		/*
		for (i = 0; i<numBins; i++)
		{
			sft_vect->data[0].data->data[i].re *= scalefactor;
			sft_vect->data[0].data->data[i].im *= scalefactor;
		}
		*/
		for (i = 0; i < numBins; i++)
		{
			thispower = ((crealf(sft_vect->data[0].data->data[i]) * scalefactor) *
				(crealf(sft_vect->data[0].data->data[i]) * scalefactor)) +
				((cimagf(sft_vect->data[0].data->data[i]) * scalefactor) *
				(cimagf(sft_vect->data[0].data->data[i]) * scalefactor));
			//thispower = sft_vect->data[0].data->data[i].re*sft_vect->data[0].data->data[i].re +
			//	sft_vect->data[0].data->data[i].im*sft_vect->data[0].data->data[i].im;
			thisavepower = 0.;
			nside = 10;
			count = 0;
			for (ii = -nside; ii <= nside; ii++)
			{
				if (i + ii >= 0 && i + ii<numBins) {
					//thisavepower += sft_vect->data[0].data->data[i + ii].re*sft_vect->data[0].data->data[i + ii].re +
					//	sft_vect->data[0].data->data[i + ii].im*sft_vect->data[0].data->data[i + ii].im;
					thisavepower += ((crealf(sft_vect->data[0].data->data[i + ii]) * scalefactor) *
						(crealf(sft_vect->data[0].data->data[i + ii]) * scalefactor)) +
						((cimagf(sft_vect->data[0].data->data[i + ii]) * scalefactor) *
						(cimagf(sft_vect->data[0].data->data[i + ii]) * scalefactor));
					count++;
				}
			}
			thisavepower /= count;
			weight = 1. / thisavepower;
			//	    weight = 1.;
			if (j == 0)
			{
				timeavg[i] = thispower;
				timeavgwt[i] = thispower*weight;
				sumweight[i] = weight;
			}
			else
			{
				timeavg[i] += thispower;
				timeavgwt[i] += thispower*weight;
				sumweight[i] += weight;
			}
		}
		// Destroys current SFT Vector
		XLALDestroySFTVector(sft_vect);
	}
	printf("About to do calculation of averages...\n");
	printf("Sample: timeavg[0]=%g, timeavgwt[0]=%g, sumweight[0]=%g\n", timeavg[0], timeavgwt[0], sumweight[0]);
	/*timeavg records the power of each bin*/
	for (i = 0; i < numBins; i++)
	{
		f = f0 + ((REAL4)i)*deltaF;
		PSD = 2.*timeavg[i] / ((REAL4)nSFT) / scalefactor / scalefactor / timebaseline;
		PSDWT = 2.*timeavgwt[i] / sumweight[i] / scalefactor / scalefactor / timebaseline;
		AMPPSD = pow(PSD, 0.5);
		AMPPSDWT = pow(PSDWT, 0.5);
		/*	SNR=(PWR-1)*(sqrt(((REAL4)nSFT))); */
		//        fprintf(fp3,"%16.8f %g %g %g %g %g %g\n",f, PWR, STRAIN, PWRWT, STRAINWT, timeavgwt[i], sumweight[i]);
		fprintf(fp3, "%16.8f %g %g %g %g\n", f, PSD, AMPPSD, PSDWT, AMPPSDWT);
		PWA_TAVGWT = timeavgwt[i];
		PWA_SUMWT = sumweight[i];
		fprintf(fp5, "%16.8f %g %g\n", f, PWA_TAVGWT, PWA_SUMWT);
		// fprintf(stderr, "%16.8f %g %g\n", f, PWA_TAVGWT, PWA_SUMWT);
	}
	/*------------------------------------------------------------------------------------------------------------------------*/
	/*End of normal spec_avg code, the remaining code is for crab freq calc.*/
	/*================================================================================================================*/

	/*fprintf(stderr,"end of spec_avg 1\n");*/

	/*=======================================================================================================================*/
	/*=======================================================================================================================*/
	fprintf(stderr,"Destroying Variables\n");
	XLALDestroySFTCatalog(catalog);/*cg; desctroys the SFT catalogue*/

	
	/*release a;; the allocaeted memory*/
	LALCHARDestroyVector(&status, &year_date);

	if (timeavg != NULL) XLALFree(timeavg);

    	XLALDestroyUserVars();
	fprintf(stderr,"Done Destroying Variables\n");
/*
	LAL_CALL(LALDestroyUserVars(&status), &status);
*/
	/*fprintf(stderr,"end of spec_avg 4\n");*/
	/*close all the files, spec_avg.c is done, all info written to the files.*/
	fprintf(stderr,"Closing Files\n");
	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);
	fprintf(stderr, "end of spec_avg\n");
	fprintf(stderr, "Spec_avg_done!\n");


	return(0);

}
/* END main */

SFTVector*
LALExtractSFT(const SFTCatalog *catalog,   /**< The 'catalogue' of SFTs to load */
REAL8 fMin,		   /**< minumum requested frequency (-1 = read from lowest) */
REAL8 fMax,		   /**< maximum requested frequency (-1 = read up to highest) */
REAL8 SFTindex
)
{
	UINT4 catPos = SFTindex;                    /**< current file in catalog */
	UINT4 firstbin, lastbin;         /**< the first and last bin we want to read */
	UINT4 minbin, maxbin;            /**< min and max bin of all SFTs in the catalog */
	UINT4 nSFTs = 1;                 /**< number of SFTs, i.e. different GPS timestamps */
	REAL8 deltaF;                    /**< frequency spacing of SFT */
	SFTCatalog locatalog;            /**< local copy of the catalog to be sorted by 'locator' */
	SFTVector* sftVector = NULL;     /**< the vector of SFTs to be returned */
	SFTReadSegment *e_segments = NULL;  /**< array of segments already read of an SFT */
	char empty = '\0';               /**< empty string */
	char* fname = &empty;            /**< name of currently open file, initially "" */
	FILE* fp = NULL;                 /**< open file */
	SFTtype* thisSFT = NULL;         /**< SFT to read from file */

	/* error handler: free memory and return with error */
#define XLALLOADSFTSERROR(eno)	{		\
    if(fp)					\
      fclose(fp);				\
    if(e_segments) 				\
      XLALFree(e_segments);			\
    if(locatalog.data)				\
      XLALFree(locatalog.data);			\
    if(thisSFT)					\
      XLALDestroySFT(thisSFT);			\
    if(sftVector)				\
      XLALDestroySFTVector(sftVector);		\
    XLAL_ERROR_NULL(eno);	                \
		  }

	/* initialize locatalog.data so it doesn't get free()d on early error */
	locatalog.data = NULL;

	/* check function parameters */
	if (!catalog)
		XLALLOADSFTSERROR(XLAL_EINVAL);

	/* determine number of SFTs, i.e. number of different GPS timestamps.
	The catalog should be sorted by GPS time, so just count changes.
	Record the 'index' of GPS time in the 'isft' field of the locator,
	so that we know later in which SFT to put this segment

	while at it, record max and min bin of all SFTs in the catalog */

	//LIGOTimeGPS epoch = catalog->data[catPos].header.epoch;
	catalog->data[catPos].locator->isft = nSFTs - 1;
	deltaF = catalog->data[catPos].header.deltaF; /* Hz/bin */
	minbin = firstbin = lround(catalog->data[catPos].header.f0 / deltaF);
	maxbin = lastbin = firstbin + catalog->data[catPos].numBins - 1;
	//for (catPos = 1; catPos < catalog->length; catPos++) {
	//	firstbin = lround(catalog->data[catPos].header.f0 / deltaF);
	//	lastbin = firstbin + catalog->data[catPos].numBins - 1;
	//	if (firstbin < minbin)
	//		minbin = firstbin;
	//	if (lastbin > maxbin)
	//		maxbin = lastbin;
	//	if (!GPSEQUAL(epoch, catalog->data[catPos].header.epoch)) {
	//		epoch = catalog->data[catPos].header.epoch;
	//		nSFTs++;
	//	}
	//	catalog->data[catPos].locator->isft = nSFTs - 1;
	//}
	XLALPrintInfo("%s: fMin: %f, fMax: %f, deltaF: %f, minbin: %u, maxbin: %u\n", __func__, fMin, fMax, deltaF, minbin, maxbin);
	//fprintf(stderr,"%s: fMin: %f, fMax: %f, deltaF: %f, minbin: %u, maxbin: %u\n", __func__, fMin, fMax, deltaF, minbin, maxbin);

	/* calculate first and last frequency bin to read */
	if (fMin < 0)
		firstbin = minbin;
	else
		firstbin = (UINT4)floor(fMin / deltaF * fudge_up);	// round *down*, but allow for 10*eps 'fudge'
	if (fMax < 0)
		lastbin = maxbin;
	else {
		lastbin = (UINT4)ceil(fMax / deltaF * fudge_down);	// round *up*, but allow for 10*eps fudge
		if ((lastbin == 0) && (fMax != 0)) {
			XLALPrintError("ERROR: last bin to read is 0 (fMax: %f, deltaF: %f)\n", fMax, deltaF);
			XLALLOADSFTSERROR(XLAL_EINVAL);
		}
	}
	XLALPrintInfo("%s: Reading from first bin: %u, last bin: %u\n", __func__, firstbin, lastbin);
	//fprintf(stderr,"%s: Reading from first bin: %u, last bin: %u\n", __func__, firstbin, lastbin);

	/* allocate the SFT vector that will be returned */
	if (!(sftVector = XLALCreateSFTVector(nSFTs, lastbin + 1 - firstbin))) {
		XLALPrintError("ERROR: Couldn't create sftVector\n");
		XLALLOADSFTSERROR(XLAL_ENOMEM);
	}

	/* allocate an additional single SFT where SFTs are read in */
	if (!(thisSFT = XLALCreateSFT(lastbin + 1 - firstbin))) {
		XLALPrintError("ERROR: Couldn't create thisSFT\n");
		XLALLOADSFTSERROR(XLAL_ENOMEM);
	}

	/* make a copy of the catalog that gets sorted by locator.
	Eases maintaing a correctly (epoch-)sorted catalog, particulary in case of errors
	note: only the pointers to the SFTdescriptors are copied & sorted, not the descriptors */
	locatalog.length = catalog->length;
	{
		UINT4 size = catalog->length * sizeof(catalog->data[0]);
		if (!(locatalog.data = XLALMalloc(size))) {
			XLALPrintError("ERROR: Couldn't allocate locatalog.data\n");
			XLALLOADSFTSERROR(XLAL_ENOMEM);
		}
		memcpy(locatalog.data, catalog->data, size);
	}

	/* sort catalog by f0, locator */
	qsort((void*)locatalog.data, locatalog.length, sizeof(locatalog.data[0]), compareSFTloc);

	/* allocate segment vector, one element per final SFT */
	if (!(e_segments = XLALCalloc(nSFTs, sizeof(SFTReadSegment)))) {
		XLALPrintError("ERROR: Couldn't allocate locatalog.data\n");
		XLALLOADSFTSERROR(XLAL_ENOMEM);
	}

	/* loop over all files (actually locators) in the catalog */
	//for (catPos = 0; catPos < catalog->length; catPos++) {

	struct tagSFTLocator* locator = locatalog.data[catPos].locator;
	//tagSFTLocator *locator = locatalog.data[catPos].locator;
	UINT4 isft = locator->isft;
	UINT4 firstBinRead;
	UINT4 lastBinRead;

	if (locatalog.data[catPos].header.data) {
		/* the SFT data has already been read into the catalog,
		copy the relevant part to thisSFT */

		volatile REAL8 tmp = locatalog.data[catPos].header.f0 / deltaF;
		UINT4 firstSFTbin = lround(tmp);
		UINT4 lastSFTbin = firstSFTbin + locatalog.data[catPos].numBins - 1;
		UINT4 firstBin2read = firstbin;
		UINT4 lastBin2read = lastbin;
		UINT4 numBins2read, offsetBins;
                printf("firstSFTbin=%d, lastSFTbin=%d\n",firstSFTbin,lastSFTbin);
                fprintf(stderr,"firstSFTbin=%d, lastSFTbin=%d\n",firstSFTbin,lastSFTbin);
		/* limit the interval to be read to what's actually in the SFT */
		if (firstBin2read < firstSFTbin)
			firstBin2read = firstSFTbin;
		if (lastBin2read > lastSFTbin)
			lastBin2read = lastSFTbin;

		/* check that requested interval is found in SFT */
		if (firstBin2read <= lastBin2read) {

			/* keep a copy of the data pointer */
			COMPLEX8Sequence*data = thisSFT->data;

			firstBinRead = firstBin2read;
			lastBinRead = lastBin2read;
			offsetBins = firstBin2read - firstSFTbin;
			numBins2read = lastBin2read - firstBin2read + 1;

			/* copy the header */
			*thisSFT = locatalog.data[catPos].header;
			/* restore data pointer */
			thisSFT->data = data;
			/* copy the data */
			memcpy(thisSFT->data->data,
				locatalog.data[catPos].header.data + offsetBins,
				numBins2read * sizeof(COMPLEX8));

			/* update the start-frequency entry in the SFT-header to the new value */
			thisSFT->f0 = 1.0 * firstBin2read * thisSFT->deltaF;

		}
		else {
			/* no data was needed from this SFT (segment) */
			firstBinRead = 0;
			lastBinRead = 0;
		}

	}
	else {
		/* SFT data had not yet been read - read it */

		/* open and close a file only when necessary, i.e. reading a different file */
		if (strcmp(fname, locator->fname)) {
			if (fp) {
				fclose(fp);
				fp = NULL;
			}
			fname = locator->fname;
			fp = fopen(fname, "rb");
			XLALPrintInfo("%s: Opening file '%s'\n", __func__, fname);
			if (!fp) {
				XLALPrintError("ERROR: Couldn't open file '%s'\n", fname);
				XLALLOADSFTSERROR(XLAL_EIO);
			}
		}

		/* seek to the position of the SFT in the file (if necessary) */
		if (locator->offset)
			if (fseek(fp, locator->offset, SEEK_SET) == -1) {
				XLALPrintError("ERROR: Couldn't seek to position %ld in file '%s'\n",
					locator->offset, fname);
				XLALLOADSFTSERROR(XLAL_EIO);
			}

		/* read SFT data */
		lastBinRead = read_sft_bins_from_fp(thisSFT, &firstBinRead, firstbin, lastbin, fp);
		XLALPrintInfo("%s: Read data from %s:%lu: %u - %u\n", __func__, locator->fname, locator->offset, firstBinRead, lastBinRead);
		//fprintf(stderr,"%s: Read data from %s:%lu: %u - %u\n", __func__, locator->fname, locator->offset, firstBinRead, lastBinRead);
		//}
		/* SFT data has been read from file or taken from catalog */

		if (lastBinRead) {
			/* data was actually read */

			if (e_segments[isft].last == 0) {

				/* no data was read for this SFT yet: must be first segment */
				if (firstBinRead != firstbin) {
					XLALPrintError("ERROR: data gap or overlap at first bin of SFT#%u (GPS %lf)"
						" expected bin %u, bin %u read from file '%s'\n",
						isft, GPS2REAL8(thisSFT->epoch),
						firstbin, firstBinRead, fname);
					XLALLOADSFTSERROR(XLAL_EIO);
				}
				e_segments[isft].first = firstBinRead;
				e_segments[isft].epoch = thisSFT->epoch;

				/* if not first segment, segment must fit at the end of previous data */
			}
			else if (firstBinRead != e_segments[isft].last + 1) {
				XLALPrintError("ERROR: data gap or overlap in SFT#%u (GPS %lf)"
					" between bin %u read from file '%s' and bin %u read from file '%s'\n",
					isft, GPS2REAL8(thisSFT->epoch),
					e_segments[isft].last, e_segments[isft].lastfrom->fname,
					firstBinRead, fname);
				XLALLOADSFTSERROR(XLAL_EIO);
			}

			/* consistency checks */
			if (deltaF != thisSFT->deltaF) {
				XLALPrintError("ERROR: deltaF mismatch (%f/%f) in SFT read from file '%s'\n",
					thisSFT->deltaF, deltaF, fname);
				XLALLOADSFTSERROR(XLAL_EIO);
			}
			if (!GPSEQUAL(e_segments[isft].epoch, thisSFT->epoch)) {
				XLALPrintError("ERROR: GPS epoch mismatch (%f/%f) in SFT read from file '%s'\n",
					GPS2REAL8(e_segments[isft].epoch), GPS2REAL8(thisSFT->epoch), fname);
				XLALLOADSFTSERROR(XLAL_EIO);
			}

			/* data is ok, add to SFT */
			e_segments[isft].last = lastBinRead;
			e_segments[isft].lastfrom = locator;
			memcpy(sftVector->data[isft].name, locatalog.data[catPos].header.name, sizeof(sftVector->data[isft].name));
			sftVector->data[isft].sampleUnits = locatalog.data[catPos].header.sampleUnits;
			memcpy(sftVector->data[isft].data->data + (firstBinRead - firstbin),
				thisSFT->data->data,
				(lastBinRead - firstBinRead + 1) * sizeof(COMPLEX8));

		}
		else if (!firstBinRead) {
			/* no needed data had been in this segment */
			XLALPrintInfo("%s: No data read from %s:%lu\n", __func__, locator->fname, locator->offset);

			/* set epoch if not yet set, if already set, check it */
			if (GPSZERO(e_segments[isft].epoch))
				e_segments[isft].epoch = thisSFT->epoch;
			else if (!GPSEQUAL(e_segments[isft].epoch, thisSFT->epoch)) {
				XLALPrintError("ERROR: GPS epoch mismatch (%f/%f) in SFT read from file '%s'\n",
					GPS2REAL8(e_segments[isft].epoch), GPS2REAL8(thisSFT->epoch), fname);
				XLALLOADSFTSERROR(XLAL_EIO);
			}

		}
		else {
			/* failed to read data */

			XLALPrintError("ERROR: Error (%u) reading SFT from file '%s'\n", firstBinRead, fname);
			XLALLOADSFTSERROR(XLAL_EIO);
		}
	}

	/* close the last file */
	if (fp) {
		fclose(fp);
		fp = NULL;
	}

	/* check that all SFTs are complete */
	for (UINT4 jsft = 0; jsft < nSFTs; jsft++) {
		if (e_segments[jsft].last == lastbin) {
			sftVector->data[jsft].f0 = 1.0 * firstbin * deltaF;
			sftVector->data[jsft].epoch = e_segments[jsft].epoch;
			sftVector->data[jsft].deltaF = deltaF;
		}
		else {
			if (e_segments[jsft].last)
				XLALPrintError("ERROR: data missing at end of SFT#%u (GPS %lf)"
				" expected bin %u, bin %u read from file '%s'\n",
				isft, GPS2REAL8(e_segments[jsft].epoch),
				lastbin, e_segments[jsft].last,
				e_segments[jsft].lastfrom->fname);


			else
				XLALPrintError("ERROR: no data could be read for SFT#%u (GPS %lf)\n",
				isft, GPS2REAL8(e_segments[jsft].epoch));
			XLALLOADSFTSERROR(XLAL_EIO);
		}
	}

	/* cleanup  */
	XLALFree(e_segments);
	XLALFree(locatalog.data);
	XLALDestroySFT(thisSFT);

	return(sftVector);

} /* LALExtractSFT() */

/* compare two SFT-descriptors by their locator (f0, file, position) */
static int
compareSFTloc(const void *ptr1, const void *ptr2)
{
	const SFTDescriptor *desc1 = ptr1;
	const SFTDescriptor *desc2 = ptr2;
	int s;
	if (desc1->header.f0 < desc2->header.f0)
		return -1;
	else if (desc1->header.f0 > desc2->header.f0)
		return 1;
	s = strcmp(desc1->locator->fname, desc2->locator->fname);
	if (!s) {
		if (desc1->locator->offset < desc2->locator->offset)
			return(-1);
		else if (desc1->locator->offset > desc2->locator->offset)
			return(1);
		else
			return(0);
	}
	return(s);
} /* compareSFTloc() */

/*
This function reads an SFT (segment) from an open file pointer into a buffer.
firstBin2read specifies the first bin to read from the SFT, lastBin2read is the last bin.
If the SFT contains fewer bins than specified, all bins from the SFT are read.
The function returns the last bin actually read, firstBinRead
is set to the first bin actually read. In case of an error, 0 is returned
and firstBinRead is set to a code further decribing the error condition.
*/
static UINT4
read_sft_bins_from_fp(SFTtype *ret, UINT4 *firstBinRead, UINT4 firstBin2read, UINT4 lastBin2read, FILE *fp)
{
	UINT4 version;
	UINT8 crc64;
	BOOLEAN swapEndian;
	UINT4 numBins2read;
	UINT4 firstSFTbin, lastSFTbin, numSFTbins;
	INT4 offsetBins;
	long offsetBytes;
	volatile REAL8 tmp;	/* intermediate results: try to force IEEE-arithmetic */

	if (!firstBinRead)
	{
		XLALPrintError("read_sft_bins_from_fp(): got passed NULL *firstBinRead\n");
		return((UINT4)-1);
	}

	*firstBinRead = 0;

	if ((ret == NULL) ||
		(ret->data == NULL) ||
		(ret->data->data == NULL))
	{
		XLALPrintError("read_sft_bins_from_fp(): got passed NULL SFT*\n");
		*firstBinRead = 1;
		return(0);
	}

	if (!fp)
	{
		XLALPrintError("read_sft_bins_from_fp(): got passed NULL FILE*\n");
		*firstBinRead = 1;
		return(0);
	}

	if (firstBin2read > lastBin2read)
	{
		XLALPrintError("read_sft_bins_from_fp(): Empty frequency-interval requested [%d, %d] bins\n",
			firstBin2read, lastBin2read);
		*firstBinRead = 1;
		return(0);
	}

  {
	  COMPLEX8Sequence*data = ret->data;
	  if (read_sft_header_from_fp(fp, ret, &version, &crc64, &swapEndian, NULL, &numSFTbins) != 0)
	  {
		  XLALPrintError("read_sft_bins_from_fp(): Failed to read SFT-header!\n");
		  *firstBinRead = 2;
		  return(0);
	  }
	  ret->data = data;
  }

  tmp = ret->f0 / ret->deltaF;
  firstSFTbin = lround(tmp);
  lastSFTbin = firstSFTbin + numSFTbins - 1;
  fprintf(stderr,"firstbin=%u, lastSFTbin=%u \n",firstSFTbin, lastSFTbin);
  /* limit the interval to be read to what's actually in the SFT */
  if (firstBin2read < firstSFTbin)
	  firstBin2read = firstSFTbin;
  if (lastBin2read > lastSFTbin)
	  lastBin2read = lastSFTbin;

  /* check that requested interval is found in SFT */
  /* return 0 (no bins read) if this isn't the case */
  if (firstBin2read > lastBin2read) {
	  *firstBinRead = 0;
	  return(0);
  }

  *firstBinRead = firstBin2read;

  offsetBins = firstBin2read - firstSFTbin;
  offsetBytes = offsetBins * 2 * sizeof(REAL4);
  numBins2read = lastBin2read - firstBin2read + 1;

  if (ret->data->length < numBins2read)
  {
	  XLALPrintError("read_sft_bins_from_fp(): passed SFT has not enough bins (%u/%u)\n",
		  ret->data->length, numBins2read);
	  *firstBinRead = 1;
	  return(0);
  }

  /* seek to the desired bins */
  if (fseek(fp, offsetBytes, SEEK_CUR) != 0)
  {
	  XLALPrintError("read_sft_bins_from_fp(): Failed to fseek() to first frequency-bin %d: %s\n",
		  firstBin2read, strerror(errno));
	  *firstBinRead = 3;
	  return(0);
  }

  /* actually read the data */
  if (numBins2read != fread(ret->data->data, sizeof(COMPLEX8), numBins2read, fp))
  {
	  XLALPrintError("read_sft_bins_from_fp(): Failed to read %d bins from SFT!\n", numBins2read);
	  *firstBinRead = 4;
	  return(0);
  }

  /* update the start-frequency entry in the SFT-header to the new value */
  ret->f0 = 1.0 * firstBin2read * ret->deltaF;

  /* take care of normalization and endian-swapping */
  if (version == 1 || swapEndian)
  {
	  UINT4 i;
	  REAL8 band = 1.0 * numSFTbins * ret->deltaF;/* need the TOTAL frequency-band in the SFT-file! */
	  REAL8 fsamp = 2.0 * band;
	  REAL8 dt = 1.0 / fsamp;

	  for (i = 0; i < numBins2read; i++)
	  {
		  REAL4 re = crealf(ret->data->data[i]);
		  REAL4 im = cimagf(ret->data->data[i]);

		  if (swapEndian)
		  {
			  endian_swap((CHAR *)&re, sizeof(re), 1);
			  endian_swap((CHAR *)&im, sizeof(im), 1);
		  }

		  /* if the SFT-file was in v1-Format: need to renormalize the data now by 'Delta t'
		  * in order to follow the correct SFT-normalization
		  * (see LIGO-T040164-01-Z, and LIGO-T010095-00)
		  */
		  if (version == 1)
		  {
			  re *= dt;
			  im *= dt;
		  }

		  ret->data->data[i] = crectf(re, im);
	  } /* for i < numBins2read */
  } /* if SFT-v1 */

  /* return last bin read */
//fprintf(stderr,"lastBin2read = %u\n",lastBin2read);
  return(lastBin2read);

} /* read_sft_bins_from_fp() */

/* Try to read an SFT-header (of ANY VALID SFT-VERSION) at the given FILE-pointer fp,
* and return the SFT-header, SFT-version-number and number of frequency-samples in the SFT.
*
* Sets the filepointer fp at the end of the header if successful, leaves it at
* initial position if not.
*
* RETURN 0 = OK, -1 on ERROR
*
* We do basic checking of compliance with the SFT-spec (<tt>LIGO-T04164-01-Z</tt>)
* as far as a single header is concerned.
*
* NOTE: fatal errors will produce a XLALPrintError() error-message, but
* non-fatal 'SFT-format'-errors will only output error-messages if lalDebugLevel > 0.
* --> this function can therefore be used to check if a given file actually contains SFTs
*
*
*/
static int
read_sft_header_from_fp(FILE *fp, SFTtype *header, UINT4 *version, UINT8 *crc64, BOOLEAN *swapEndian, CHAR **SFTcomment, UINT4 *numBins)
{
	SFTtype XLAL_INIT_DECL(head);
	UINT4 nsamples;
	CHAR *comm = NULL;
	UINT8 ref_crc = 0;
	UINT8 header_crc;

	UINT4 ver;
	BOOLEAN need_swap;
	long save_filepos;

	if (!header || !version || !numBins || !fp)
	{
		XLALPrintError("\nERROR read_sft_header_from_fp(): called with NULL input\n\n");
		return -1;
	}
	if (SFTcomment && ((*SFTcomment) != NULL))
	{
		XLALPrintError("\nERROR: Comment-string passed to read_sft_header_from_fp() is not empty!\n\n");
		return -1;
	}

	/* store fileposition for restoring in case of failure */
	if ((save_filepos = ftell(fp)) == -1)
	{
		XLALPrintError("\nftell() failed: %s\n\n", strerror(errno));
		return -1;
	}

	if (read_SFTversion_from_fp(&ver, &need_swap, fp) != 0)
		return -1;


	/* read this SFT-header with version-specific code */
	XLAL_INIT_MEM(head);

	switch (ver)
	{
	case 1:
		if (read_v1_header_from_fp(fp, &head, &nsamples, need_swap) != 0)
			goto failed;
		break;

	case 2:
		if (read_v2_header_from_fp(fp, &head, &nsamples, &header_crc, &ref_crc, &comm, need_swap) != 0)
			goto failed;
		break;

	default:
		XLALPrintError("\nUnsupported SFT-version %d.\n\n", ver);
		goto failed;
		break;
	} /* switch(ver) */


	/* ----- some general SFT-header consistency-checks */
	if ((head.epoch.gpsSeconds < 0) || (head.epoch.gpsNanoSeconds < 0) || (head.epoch.gpsNanoSeconds >= 1000000000))
	{
		XLALPrintError("\nInvalid GPS-epoch in SFT : [%d, %d]!\n\n",
			head.epoch.gpsSeconds, head.epoch.gpsNanoSeconds);
		goto failed;
	}

	if (head.deltaF <= 0)
	{
		XLALPrintError("\nNegative frequency-spacing in SFT!\n\n");
		goto failed;
	}

	if (head.f0 < 0)
	{
		XLALPrintError("\nNegative start-frequency in SFT!\n\n");
		goto failed;
	}

	/* ok */
	(*header) = head;
	(*version) = ver;

	if (SFTcomment)	  /* return of comment is optional */
		(*SFTcomment) = comm;
	else
		if (comm) LALFree(comm);

	(*swapEndian) = need_swap;
	(*crc64) = ref_crc;
	(*numBins) = nsamples;
	return 0;

	/* ---------- */
failed:
	/* restore filepointer initial position  */
	if (fseek(fp, save_filepos, SEEK_SET) == -1)
		XLALPrintError("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno));

	/* free comment  if we allocated one */
	if (comm)
		LALFree(comm);

	return -1;

} /* read_sft_header_from_fp() */

/* a little endian-swapper needed for SFT reading/writing */
static void
endian_swap(CHAR * pdata, size_t dsize, size_t nelements)
{
	UINT4 i, j, indx;
	CHAR tempbyte;

	if (dsize <= 1) return;

	for (i = 0; i<nelements; i++)
	{
		indx = dsize;
		for (j = 0; j<dsize / 2; j++)
		{
			tempbyte = pdata[j];
			indx = indx - 1;
			pdata[j] = pdata[indx];
			pdata[indx] = tempbyte;
		}

		pdata = pdata + dsize;
	}

	return;

} /* endian swap */

/**
* Read valid SFT version-number at position fp, and determine if we need to
* endian-swap the data.
* Restores filepointer to original position before returning.
*
* RETURN: 0 = OK, -1 = ERROR
*/
static int
read_SFTversion_from_fp(UINT4 *version, BOOLEAN *need_swap, FILE *fp)
{
	long save_filepos;
	REAL8 ver;

	/* store fileposition for restoring in case of failure */
	if ((save_filepos = ftell(fp)) == -1)
	{
		XLALPrintError("\nftell() failed: %s\n\n", strerror(errno));
		return -1;
	}

	/* read version-number */
	if (1 != fread(&ver, sizeof(ver), 1, fp))
	{
		if (lalDebugLevel) XLALPrintError("\nCould not read version-number from file\n\n");
		goto failed;
	}


	/* figure out endian-ness and check version-range */
	for (*version = MAX_SFT_VERSION; *version >= MIN_SFT_VERSION; --(*version))
	{
		REAL8 vertest = *version;
		if (!memcmp(&ver, &vertest, sizeof(ver))) {
			*need_swap = FALSE;
			break;
		}
		endian_swap((char*)(&vertest), sizeof(vertest), 1);
		if (!memcmp(&ver, &vertest, sizeof(ver))) {
			*need_swap = TRUE;
			break;
		}
	}
	if (*version < MIN_SFT_VERSION) {
		if (lalDebugLevel) {
			unsigned char *v = (unsigned char*)(&ver);
			XLALPrintError("\nERROR: illegal SFT-version (%X %X %X %X %X %X %X %X) not within [%.0f, %.0f]\n",
				v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7],
				(float)MIN_SFT_VERSION, (float)MAX_SFT_VERSION);
		}
		goto failed;
	}

	/* restore initial filepointer position */
	if (fseek(fp, save_filepos, SEEK_SET) == -1)
	{
		XLALPrintError("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno));
		goto failed;
	}

	return 0;

failed:
	fseek(fp, save_filepos, SEEK_SET);
	return -1;

} /* read_SFTversion_from_fp() */

/* ----- SFT v1 -specific header-reading function:
*
* return general SFTtype header, place filepointer at the end of the header if it succeeds,
* set fp to initial position if it fails.
* RETURN: 0 = OK, -1 = ERROR
*
* NOTE: fatal errors will produce a XLALPrintError() error-message, but
* non-fatal 'SFT-format'-errors will only output error-messages if lalDebugLevel > 0.
*
*/
static int
read_v1_header_from_fp(FILE *fp, SFTtype *header, UINT4 *nsamples, BOOLEAN swapEndian)
{
	_SFT_header_v1_t rawheader;
	long save_filepos;

	if (!fp || !header || !nsamples)
	{
		XLALPrintError("\nERROR read_v1_header_from_fp(): called with NULL input!\n\n");
		return -1;
	}

	/* store fileposition for restoring in case of failure */
	if ((save_filepos = ftell(fp)) == -1)
	{
		XLALPrintError("\nftell() failed: %s\n\n", strerror(errno));
		return -1;
	}

	/* read the whole header */
	if (fread(&rawheader, sizeof(rawheader), 1, fp) != 1)
	{
		if (lalDebugLevel) XLALPrintError("\nCould not read v1-header. %s\n\n", strerror(errno));
		goto failed;
	}

	if (swapEndian)
	{
		endian_swap((CHAR*)(&rawheader.version), sizeof(rawheader.version), 1);
		endian_swap((CHAR*)(&rawheader.gps_sec), sizeof(rawheader.gps_sec), 1);
		endian_swap((CHAR*)(&rawheader.gps_nsec), sizeof(rawheader.gps_nsec), 1);
		endian_swap((CHAR*)(&rawheader.tbase), sizeof(rawheader.tbase), 1);
		endian_swap((CHAR*)(&rawheader.first_frequency_index), sizeof(rawheader.first_frequency_index), 1);
		endian_swap((CHAR*)(&rawheader.nsamples), sizeof(rawheader.nsamples), 1);
	}

	/* double-check version-number */
	if (rawheader.version != 1)
	{
		XLALPrintError("\nWrong SFT-version in read_v1_header_from_fp()\n\n");
		goto failed;
	}

	if (rawheader.nsamples <= 0)
	{
		XLALPrintError("\nNon-positive number of samples in SFT!\n\n");
		goto failed;
	}


	/* ok: */
	memset(header, 0, sizeof(*header));

	/* NOTE: v1-SFTs don't contain a detector-name, in which case we set it to '??' */
	strcpy(header->name, "??");

	header->epoch.gpsSeconds = rawheader.gps_sec;
	header->epoch.gpsNanoSeconds = rawheader.gps_nsec;
	header->deltaF = 1.0 / rawheader.tbase;
	header->f0 = rawheader.first_frequency_index / rawheader.tbase;

	(*nsamples) = rawheader.nsamples;

	return 0;

failed:
	/* restore filepointer initial position  */
	if (fseek(fp, save_filepos, SEEK_SET) == -1)
		XLALPrintError("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno));

	return -1;

} /* read_v1_header_from_fp() */

/* ----- SFT v2 -specific header-reading function:
*
* return general SFTtype header, place filepointer at the end of the header if it succeeds,
* set fp to initial position if it fails.
* RETURN: 0 = OK, -1 = ERROR
*
* NOTE: fatal errors will produce a XLALPrintError() error-message, but
* non-fatal 'SFT-format'-errors will only output error-messages if lalDebugLevel > 0.
*
*/
static int
read_v2_header_from_fp(FILE *fp, SFTtype *header, UINT4 *nsamples, UINT8 *header_crc64, UINT8 *ref_crc64, CHAR **SFTcomment, BOOLEAN swapEndian)
{
	_SFT_header_v2_t rawheader;
	long save_filepos;
	CHAR *comm = NULL;
	UINT8 crc;


	/* check input-consistency */
	if (!fp || !header || !nsamples || !SFTcomment)
	{
		XLALPrintError("\nERROR read_v2_header_from_fp(): called with NULL input!\n\n");
		return -1;
	}
	if (SFTcomment && (*SFTcomment != NULL))
	{
		XLALPrintError("\nERROR: Comment-string passed to read_v2_header_from_fp() is not NULL!\n\n");
		return -1;
	}

	/* store fileposition for restoring in case of failure */
	if ((save_filepos = ftell(fp)) == -1)
	{
		XLALPrintError("\nERROR: ftell() failed: %s\n\n", strerror(errno));
		return -1;
	}

	/* read the whole header */
	if (fread(&rawheader, sizeof(rawheader), 1, fp) != 1)
	{
		if (lalDebugLevel) XLALPrintError("\nCould not read v2-header. %s\n\n", strerror(errno));
		goto failed;
	}

	/* ----- compute CRC for the header:
	* NOTE: the CRC checksum is computed on the *bytes*, not the numbers,
	* so this must be computed before any endian-swapping.
	*/
  {
	  UINT8 save_crc = rawheader.crc64;
	  rawheader.crc64 = 0;

	  crc = calc_crc64((const CHAR*)&rawheader, sizeof(rawheader), ~(0ULL));

	  rawheader.crc64 = save_crc;
	  /* NOTE: we're not done with crc yet, because we also need to
	  * include the comment's CRC , see below
	  */
  }/* compute crc64 checksum */

  /* ----- swap endian-ness if required ----- */
  if (swapEndian)
  {
	  endian_swap((CHAR*)(&rawheader.version), sizeof(rawheader.version), 1);
	  endian_swap((CHAR*)(&rawheader.gps_sec), sizeof(rawheader.gps_sec), 1);
	  endian_swap((CHAR*)(&rawheader.gps_nsec), sizeof(rawheader.gps_nsec), 1);
	  endian_swap((CHAR*)(&rawheader.tbase), sizeof(rawheader.tbase), 1);
	  endian_swap((CHAR*)(&rawheader.first_frequency_index), sizeof(rawheader.first_frequency_index), 1);
	  endian_swap((CHAR*)(&rawheader.nsamples), sizeof(rawheader.nsamples), 1);

	  /* v2-specific */
	  endian_swap((CHAR*)(&rawheader.crc64), sizeof(rawheader.crc64), 1);
	  endian_swap((CHAR*)(&rawheader.comment_length), sizeof(rawheader.comment_length), 1);
	  /* ----- */

  } /* if endian_swap */

  /* double-check version-number */
  if (rawheader.version != 2)
  {
	  XLALPrintError("\nWrong SFT-version %g in read_v2_header_from_fp()\n\n", rawheader.version);
	  goto failed;
  }

  if (rawheader.nsamples <= 0)
  {
	  XLALPrintError("\nNon-positive number of samples in SFT!\n\n");
	  goto failed;
  }

  /* ----- v2-specific consistency-checks ----- */
  if (rawheader.comment_length < 0)
  {
	  XLALPrintError("\nNegative comment-length in SFT!\n\n");
	  goto failed;
  }

  if (rawheader.comment_length % 8 != 0)
  {
	  XLALPrintError("\nComment-length must be multiple of 8 bytes!\n\n");
	  goto failed;
  }

  if (!is_valid_detector(rawheader.detector))
  {
	  XLALPrintError("\nIllegal detector-name in SFT: '%c%c'\n\n",
		  rawheader.detector[0], rawheader.detector[1]);
	  goto failed;
  }

  /* ----- Now read comment (if any) ----- */
  comm = NULL;
  if (rawheader.comment_length)
  {
	  CHAR *ptr;
	  if ((comm = LALCalloc(1, rawheader.comment_length)) == NULL)
	  {
		  XLALPrintError("\nFATAL: out of memory ...\n\n");
		  goto failed;
	  }
	  if ((size_t)rawheader.comment_length != fread(comm, 1, rawheader.comment_length, fp))
	  {
		  XLALPrintError("\nCould not read %d-bytes comment\n\n", rawheader.comment_length);
		  goto failed;
	  }

	  /* check that comment is 0-terminated */
	  if (comm[rawheader.comment_length - 1] != 0)
	  {
		  XLALPrintError("\nComment is not properly 0-terminated!\n\n");
		  goto failed;
	  }

	  /* check that no NON-NULL bytes after first NULL in comment (->spec) */
	  ptr = strchr(comm, 0);	/* guaranteed to find sth, after previous check */
	  while (ptr < (comm + rawheader.comment_length - 1))
		  if (*ptr++ != 0)
		  {
			  XLALPrintError("\nNon-NULL bytes found after comment-end!\n\n");
			  goto failed;
		  }

	  /* comment length including null terminator to string must be an
	  * integer multiple of eight bytes. comment==NULL means 'no
	  * comment'
	  */
	  if (SFTcomment)
	  {
		  CHAR pad[] = { 0, 0, 0, 0, 0, 0, 0 };	/* for comment-padding */
		  UINT4 comment_len = strlen(comm) + 1;
		  UINT4 pad_len = (8 - (comment_len % 8)) % 8;

		  crc = calc_crc64((const CHAR*)comm, comment_len, crc);
		  crc = calc_crc64((const CHAR*)pad, pad_len, crc);
	  }

  } /* if comment_length > 0 */

  /*  ok: */
  memset(header, 0, sizeof(*header));

  header->name[0] = rawheader.detector[0];
  header->name[1] = rawheader.detector[1];
  header->name[2] = 0;

  header->epoch.gpsSeconds = rawheader.gps_sec;
  header->epoch.gpsNanoSeconds = rawheader.gps_nsec;

  header->f0 = rawheader.first_frequency_index / rawheader.tbase;
  header->deltaF = 1.0 / rawheader.tbase;

  (*nsamples) = rawheader.nsamples;
  (*ref_crc64) = rawheader.crc64;
  (*SFTcomment) = comm;
  (*header_crc64) = crc;


  return 0;

failed:
  /* restore filepointer initial position  */
  if (fseek(fp, save_filepos, SEEK_SET) == -1)
	  XLALPrintError("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno));

  /* free comment  if we allocated one */
  if (comm)
	  LALFree(comm);

  return -1;

} /* read_v2_header_from_fp() */

/* The crc64 checksum of M bytes of data at address data is returned
* by crc64(data, M, ~(0ULL)). Call the function multiple times to
* compute the checksum of data made in contiguous chunks, setting
* final argument to the previously accumulated checksum value. */
static UINT8
calc_crc64(const CHAR *data, UINT4 length, UINT8 crc)
{
	UINT8 CRCTable[TABLELEN];
	UINT4 i;

	/* is there is no data, simply return previous checksum value */
	if (!length || !data)
		return crc;

	/* initialize the CRC table for fast computation.  We could keep
	this table in memory to make the computation faster, but that is
	not re-entrant for multi-threaded code.
	*/
	for (i = 0; i < TABLELEN; i++) {
		UINT4 j;
		UINT8 part = i;
		for (j = 0; j < 8; j++) {
			if (part & 1)
				part = (part >> 1) ^ POLY64;
			else
				part >>= 1;
		}
		CRCTable[i] = part;
	}

	/* compute the CRC-64 code */
	for (i = 0; i<length; i++) {
		UINT8 temp1 = crc >> 8;
		UINT8 temp2 = CRCTable[(crc ^ (UINT8)data[i]) & 0xff];
		crc = temp1 ^ temp2;
	}

	return crc;

} /* calc_crc64() */

/* check that channel-prefix defines a 'known' detector.  The list of
* known detectors implemented here for now follows the list in
* Appendix D of LIGO-T970130-F-E:
*
* returns TRUE if valid, FALSE otherwise */
static BOOLEAN
is_valid_detector(const char *channel)
{
	int i;
	const char *knownDetectors[] =
	{
		"A1",       /* ALLEGRO */
		"B1",       /* NIOBE */
		"E1",       /* EXPLORER */
		"G1",       /* GEO_600 */
		"H1",       /* LHO_4k */
		"H2",       /* LHO_2k */
		"K1",       /* ACIGA */
		"L1",       /* LLO_4k */
		"N1",       /* Nautilus */
		"O1",       /* AURIGA */
		"P1",       /* CIT_40 */
		"T1",       /* TAMA_300 */
		"V1",       /* Virgo_CITF */
		"V2",       /* Virgo (3km) */
		"Z1",	  /* LISA effective IFO 1 */
		"Z2",	  /* LISA effective IFO 2 */
		"Z3",	  /* LISA effective IFO 3 */
		"Z4",	  /* LISA effective IFO 2 minus 3 */
		"Z5",	  /* LISA effective IFO 3 minus 1 */
		"Z6",	  /* LISA effective IFO 1 minus 2 */
		"Z7",	  /* LISA pseudo TDI A */
		"Z8",	  /* LISA pseudo TDI E */
		"Z9",	  /* LISA pseudo TDI T */
		"X1",       /* RXTE PCA */
		"X2",       /* RXTE ASM */
		NULL
	};

	if (!channel)
		return FALSE;

	if (strlen(channel) < 2)
		return FALSE;

	for (i = 0; knownDetectors[i]; i++)
	{
		if ((knownDetectors[i][0] == channel[0]) && (knownDetectors[i][1] == channel[1]))
			return TRUE;
	}

	return FALSE;

} /* is_valid_detector() */
