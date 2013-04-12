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
 * \ingroup pulsarApps
 */

#define LAL_USE_OLD_COMPLEX_STRUCTS

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

#define NUM 1000 
/*used for defining structures such as crabOutput*/

int main(int argc, char **argv)
{
    FILE *fp  = NULL;
    FILE *fp2 = NULL;
    FILE *fp3 = NULL;
    FILE *fp4 = NULL;
    LALStatus status = blank_status;
    
    SFTCatalog *catalog = NULL;
    SFTVector *sft_vect = NULL;
    INT4 i,j,k,l;
    INT4 numBins, nSFT;
    SFTConstraints constraints=empty_SFTConstraints;
    LIGOTimeGPS startTime, endTime; 
    REAL8 avg =0;
    REAL4 *timeavg =NULL;
    REAL4 PWR,SNR;
    REAL8 f =0;
    CHAR outbase[256],outfile[256],outfile2[256],outfile3[256], outfile4[256]; /*, outfile6[256]; */
    REAL8 NumBinsAvg =0;
    REAL8 timebaseline =0;
    
    BOOLEAN help = 0;
    CHAR *SFTpatt = NULL;
    CHAR *IFO = NULL;
    INT4 startGPS = 0;
    INT4 endGPS = 0;
    REAL8 f_min = 0.0;
    REAL8 f_max = 0.0;
    REAL8 freqres =0.0;
    INT4 blocksRngMed = 101;
    CHAR *outputBname = NULL;
    INT4 cur_epoch = 0, next_epoch = 0;
    
    /* these varibales are for converting GPS seconds into UTC time and date*/
    struct tm         date;
    CHARVector        *timestamp = NULL;
    CHARVector	     *year_date = NULL;
    REAL8Vector     *timestamps=NULL;
    
    CHAR *psrInput = NULL;
    CHAR *psrEphemeris = NULL;
    CHAR *earthFile = NULL;
    CHAR *sunFile = NULL;
  /*========================================================================================================================*/
    
    lalDebugLevel = 0;
    LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);
    
    LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",         'h', UVAR_HELP,     "Print this help message",     &help        ), &status);
    LAL_CALL(LALRegisterSTRINGUserVar(&status, "SFTs",         'p', UVAR_REQUIRED, "SFT location/pattern",        &SFTpatt     ), &status);
    LAL_CALL(LALRegisterSTRINGUserVar(&status, "IFO",          'I', UVAR_REQUIRED, "Detector",                    &IFO         ), &status);
    LAL_CALL(LALRegisterINTUserVar   (&status, "startGPS",     's', UVAR_REQUIRED, "Starting GPS time",           &startGPS    ), &status);
    LAL_CALL(LALRegisterINTUserVar   (&status, "endGPS",       'e', UVAR_REQUIRED, "Ending GPS time",             &endGPS      ), &status);
    LAL_CALL(LALRegisterREALUserVar  (&status, "fMin",         'f', UVAR_REQUIRED, "Minimum frequency",           &f_min       ), &status);
    LAL_CALL(LALRegisterREALUserVar  (&status, "fMax",         'F', UVAR_REQUIRED, "Maximum frequency",           &f_max       ), &status);
    LAL_CALL(LALRegisterINTUserVar   (&status, "blocksRngMed", 'w', UVAR_OPTIONAL, "Running Median window size",  &blocksRngMed), &status);
    LAL_CALL(LALRegisterSTRINGUserVar(&status, "outputBname",  'o', UVAR_OPTIONAL, "Base name of output files",   &outputBname ), &status);
    LAL_CALL(LALRegisterREALUserVar  (&status, "freqRes",      'r', UVAR_REQUIRED, "Spectrogram freq resolution", &freqres     ), &status);
    LAL_CALL(LALRegisterREALUserVar  (&status, "timeBaseline", 't', UVAR_REQUIRED, "The time baseline of sfts",   &timebaseline), &status);
    LAL_CALL(LALRegisterSTRINGUserVar(&status, "psrInput",     'P', UVAR_OPTIONAL, "name of tempo pulsar file",   &psrInput ), &status);
    LAL_CALL(LALRegisterSTRINGUserVar(&status, "psrEphemeris", 'S', UVAR_OPTIONAL, "pulsar ephemeris file",   &psrEphemeris ), &status);
    LAL_CALL(LALRegisterSTRINGUserVar(&status, "earthFile",  'y', UVAR_OPTIONAL, "earth .dat file",   &earthFile ), &status);
    LAL_CALL(LALRegisterSTRINGUserVar(&status, "sunFile",   'z', UVAR_OPTIONAL, "sun .dat file",   &sunFile ), &status);
    
    LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);
    if (help)
    return(0);
    
    startTime.gpsSeconds = startGPS;/*cg; startTime is a structure, and gpsSeconds is a member of that structure*/
    startTime.gpsNanoSeconds = 0;/*cg; gps NanoSeconds is also a member of the startTime structure */
    constraints.startTime = &startTime; /*cg; & operator gets the address of variable, &a is a pointer to a.  This line puts the startTime structure into the structure constraints*/
    
    endTime.gpsSeconds = endGPS;
    endTime.gpsNanoSeconds = 0;
    constraints.endTime = &endTime;/*cg; This line puts the end time into the structure constraints*/
    constraints.detector = IFO;/*cg; this adds the interferometer into the contraints structure*/
    LALSFTdataFind ( &status, &catalog, SFTpatt, &constraints );/*cg; creates SFT catalog, uses the constraints structure*/

    if (catalog == NULL)/*need to check for a NULL pointer, and print info about circumstances if it is null*/
    {
        fprintf(stderr, "SFT catalog pointer is NULL!  There has been an error with LALSFTdataFind\n");
        fprintf(stderr, "LALStatus info.... status code: %d, message: %s, offending function: %s\n", status.statusCode, status.statusDescription, status.function);
        exit(0);
    }
    if (catalog->length == 0)
    {
        fprintf(stderr, "No SFTs found, please exmanine start time, end time, frequency range etc\n");
        exit(0);
    }

    LALLoadSFTs ( &status, &sft_vect, catalog, f_min,f_max);/*cg;reads the SFT data into the structure sft_vect*/

    if (sft_vect == NULL)
    {
        fprintf(stderr, "SFT vector pointer is NULL!  There has been an error with LALLoadSFTs\n");
        fprintf(stderr, "LALStatus info.... status code: %d, message: %s, offending function: %s\n", status.statusCode, status.statusDescription, status.function);
        exit(0);
    }

    LALDestroySFTCatalog( &status, &catalog);/*cg; desctroys the SFT catalogue*/
    numBins = sft_vect->data->data->length;/*the number of bins in the freq_range*/
    nSFT = sft_vect->length;/* the number of sfts.*/
    
    fprintf(stderr, "nSFT = %d\tnumBins = %d\tf0 = %f\n", nSFT, numBins,sft_vect->data->f0);/*print->logs/spectrumAverage_testcg_0.err */
    if (LALUserVarWasSet(&outputBname))
    strcpy(outbase, outputBname);
    else
    sprintf(outbase, "spec_%.2f_%.2f_%s_%d_%d", f_min,f_max,constraints.detector,startTime.gpsSeconds,endTime.gpsSeconds);/*cg; this is the default name for producing the output files, the different suffixes are just added to this*/
    sprintf(outfile,  "%s", outbase);/*cg; name of first file to be output*/
    sprintf(outfile2, "%s_timestamps", outbase);/*cg: name of second file to be output*/
    sprintf(outfile3, "%s.txt", outbase);/*cg; name of third file to be output*/
    sprintf(outfile4, "%s_date", outbase);/*cg;file for outputting the date, which is used in matlab plotting.*/

    fp = fopen(outfile, "w");/*cg;  open all three files for writing, if they don't exist create them, if they do exist overwrite them*/
    fp2 = fopen(outfile2, "w");
    fp3 = fopen(outfile3, "w");
    fp4 = fopen(outfile4, "w");

/*----------------------------------------------------------------------------------------------------------------*/
/*cg; Create the first file called spec_blah_blah.  This file outputs the power in each spectrogram bin.  The number of bins depends on the frequency range, freq resolution, and number of SFTs.*/

/*cg;  Create the second file, called    blah_b;ah_blah_timestamps.  This will simply contain the time in GPS seconds of each SFT.*/

    NumBinsAvg = freqres*numBins/(f_max-f_min);/*this calcs the number of bins over which to average the sft data, this is worked out so it produces the same freq resolution as specified in the arguments passed to fscanDriver.py. numBins is the total number of bins in the raw sft data*/

    l=0;/*l is used as a counter to count how many SFTs and fake zero SFTs (used for gaps) are output for specgram.*/
    timestamps = XLALCreateREAL8Vector(l);/*test for getting rid of second loop*/
    LALCHARCreateVector(&status, &year_date, (UINT4)128); 

  /*create output files and check for missing sfts*/
    for (j=0;j<nSFT;j++)/*cg;nSFT is the numnber of SFT files used for the time specified. So process is repeated for each SFT*/
    {
        cur_epoch = sft_vect->data[j].epoch.gpsSeconds;/*finds the gps time of the current sft in the sequence with index j*/
        fprintf(fp2, "%d.\t%d\n", l, cur_epoch);/*cg; this bit writes the second file, i.e. the timestamps*/
    
        XLALResizeREAL8Vector(timestamps, l+1);/*resizes the vector timestamps, so cur_epoch can be added*/
        timestamps->data[l]= cur_epoch;/*number of gaps is not know in advance, hence need for te resizing*/
        XLALGPSToUTC(&date, cur_epoch);/*cg; gets the UTC date in struct tm format from the GPS seconds.*/
        fprintf(fp4, "%d\t %i\t %i\t %i\t %i\t %i\t %i\n", l, (date.tm_year+1900), date.tm_mon+1, date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec);
    
        for ( i=0; i < (numBins-2); i+=NumBinsAvg)/*cg; this loop works out the powers and writes the first file.*/
        {/*cg; each SFT is split up into a number of bins, the number of bins is read in from the SFT file*/
            avg = 0.0;/*cg; the vairable avg is reset each time.*/
            if (i+NumBinsAvg>numBins) {printf("Error\n");return(2);}/*cg; error is detected, to prevent referencing data past the end of sft_vect.*/
            for (k=0;k<NumBinsAvg;k++)/*cg; for each bin, k goes trhough each entry from 0 to 180.*/
                avg += sqrt(2*(crealf(sft_vect->data[j].data->data[i+k])*crealf(sft_vect->data[j].data->data[i+k]) + 
                cimagf(sft_vect->data[j].data->data[i+k])*cimagf(sft_vect->data[j].data->data[i+k]))/timebaseline);/*cg; re amd im are real and imaginary parts of SFT, duh!*/
            fprintf(fp,"%e\t",avg/NumBinsAvg);
        }
        fprintf(fp,"\n");
        /*------------------------------*/
        /*Bit to check if there is a gap in the sfts*/
        if ( j < (nSFT-1) )/*in all cases except when we are examining the last sft, check that there is no gap to the next sft*/
        {
            next_epoch = sft_vect->data[j+1].epoch.gpsSeconds;
            /*test to see if SFT gap is longer than 3/2*timebaseline, if so create another entry in the matrix*/
            while ((cur_epoch+((3/2)*timebaseline)) < next_epoch)
            {
                for ( i=0; i < (numBins-2); i+=NumBinsAvg)
                {
                    avg = 0.0;
                    if (i+NumBinsAvg>numBins) {printf("Error\n");return(2);}
                    fprintf(fp,"%e\t",avg);
                }
                fprintf(fp,"\n");
                l=l+1;
                cur_epoch=cur_epoch+timebaseline;
                fprintf(fp2, "%d.\t%d\n", l, cur_epoch );
                    
                XLALResizeREAL8Vector(timestamps, l+1);
                timestamps->data[l]= cur_epoch;
                XLALGPSToUTC(&date, cur_epoch);/*cg; gets the UTC date in struct tm format from the GPS seconds.*/
                fprintf(fp4, "%d\t %i\t %i\t %i\t %i\t %i\t %i\n", l, (date.tm_year+1900), date.tm_mon+1, date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec);
            }
            
        }
        l=l+1;
    }
    fprintf(stderr,"finished checking for missing sfts, l=%d\n", l);
/*----------------------------------------------------------------------------------------------------------------*/
/*cg;  Create the third and final file, called   blah_blah_blah.txt.  This file will contain the data used in the matlab plot script to plot the normalised average power vs the frequency.*/

    /* Find time average of normalized SFTs */
    LALNormalizeSFTVect(&status, sft_vect, blocksRngMed);   
    LALNormalizeSFTVect(&status, sft_vect, blocksRngMed);   
    timeavg = XLALMalloc(numBins*sizeof(REAL4));
    if (timeavg == NULL) fprintf(stderr,"Timeavg memory not allocated\n");

    for (j=0;j<nSFT;j++)
    { 
        for ( i=0; i < numBins; i++)
        {
            if (j == 0) 
            {
                timeavg[i] = crealf(sft_vect->data[j].data->data[i])*crealf(sft_vect->data[j].data->data[i]) + 
                            cimagf(sft_vect->data[j].data->data[i])*cimagf(sft_vect->data[j].data->data[i]);
            } 
            else 
            {
                timeavg[i] += crealf(sft_vect->data[j].data->data[i])*crealf(sft_vect->data[j].data->data[i]) + 
                            cimagf(sft_vect->data[j].data->data[i])*cimagf(sft_vect->data[j].data->data[i]);
            }
        }
    }
    /*timeavg records the power of each bin*/
    for ( i=0; i < numBins; i++)
    {
        f = sft_vect->data->f0 + ((REAL4)i)*sft_vect->data->deltaF;
	PWR=timeavg[i]/((REAL4)nSFT);
	SNR=(PWR-1)*(sqrt(((REAL4)nSFT)));
        fprintf(fp3,"%16.8f %g %g\n",f, PWR, SNR);
    } 
/*------------------------------------------------------------------------------------------------------------------------*/ 
/*End of normal spec_avg code, the remaining code is for crab freq calc.*/
/*================================================================================================================*/
/*================================================================================================================*/
    /*This next block of code is for the crab specific changes to fscan*/
    #define CRAB 0
    /*change this to CRAB 0 to prevent this section of code from compiling, change to 1 to compile it.*/
    #if CRAB
    /*--------------------------------------------------------------------------------------------------------------*/
    /*some header files for the crab*/
    /*#include "../TDS_isolated/HeterodyneCrabPulsar.h"*/
    /*#include "../TDS_isolated/heterodyne_pulsar.h"*/
    #include<../TDS_isolated/HeterodyneCrabPulsar.h>

    if (psrInput != NULL){

    fprintf(stderr,"--------------------\n\n");
    fprintf(stderr,"start of crab stuff\n");
    LIGOTimeGPS dataEpoch;
    /*below 4 structures are from HeterodyneCrabPulsar.h*/
    GetCrabEphemerisInput input; /*this is needed to get the crab ephemeris*/
    CrabSpindownParamsInput crabEphemerisData;
    CrabSpindownParamsOutput crabOutput;
    ParamsForHeterodyne hetParams;
    
    /*CG; these lines allocate memory for the crab ephemeris and crab output variables...*/
    crabEphemerisData.f1 = NULL;
    LALDCreateVector( &status, &crabEphemerisData.f1, NUM);
    
    crabEphemerisData.f0 = NULL;
    LALDCreateVector( &status, &crabEphemerisData.f0, NUM);
    
    crabEphemerisData.tArr = NULL;
    LALDCreateVector( &status, &crabEphemerisData.tArr, NUM);
    
    crabOutput.tArr = NULL;
    LALDCreateVector( &status, &crabOutput.tArr, NUM);
        
    crabOutput.f0 = NULL;
    LALDCreateVector( &status, &crabOutput.f0, NUM);
    
    crabOutput.f1 = NULL;
    LALDCreateVector( &status, &crabOutput.f1, NUM);
    
    crabOutput.f2 = NULL;
    LALDCreateVector( &status, &crabOutput.f2, NUM);
    
    crabOutput.f3 = NULL;
    LALDCreateVector( &status, &crabOutput.f3, NUM);
    
    crabOutput.f4 = NULL;
    LALDCreateVector( &status, &crabOutput.f4, NUM);
        
        
    /*This next set of variables are to do with the doppler shifts that are then applied to to the crab feq.*/
    REAL8 t2=0., tdt=0.;
    EphemerisData *edat=NULL;
    BarycenterInput baryinput, baryinput2;
    EarthState earth, earth2;
    EmissionTime  emit, emit2;
    REAL8 df=0., freq, finalFreq, max_df, minf, maxf;
    REAL8 dtpos=0.; /* time between position epoch and data timestamp */
    BinaryPulsarParams pulsarParams; /*general pulsar params strcut, despite binary name*/
    /*CHAR *psrInput = NULL;*/ /* pulsar input file containing params f0, f1 etc. */
    LALDetector det;
    CHAR detName[256];
    REAL8 ecliptic_lat, e_tilt=0.409092627;/*the ecliptic latitude of the source, and hte earth's tilt in radians (23.439281 degrees)*/

    char outfile5[256];
    FILE *fp5 = NULL;
    sprintf(outfile5, "%s_crab", outbase);
    fp5 = fopen(outfile5, "w");

    /*----------------------------------------------------------------------------------------------------------------*/
    /*cg; calculating the crab freq, done in three steps. This is done for each sft*/
    /*----------------------------------------------------------------------------------------------------------------*/
    /*    ---1---   */
    /*Find rough guess of crabs freq from ephemeris*/
    /*----------------------------------------------------------------------------------------------------------------*/
    
    /*Get detector position, this is needed for barycentre calcs*/
    det = *XLALGetSiteInfo( IFO );

    /* read in tempo par file for pulsar, This is one of the optional command line arguments*/
    /*fprintf(stderr,"%s\n",psrInput);*/
    fprintf(stderr,"%s\n",psrInput);
    XLALReadTEMPOParFile(&pulsarParams, psrInput);

    /*Make sure that posepoch and pepoch are set*/
    if(pulsarParams.pepoch == 0. && pulsarParams.posepoch != 0.)
	pulsarParams.pepoch = pulsarParams.posepoch;
    else if(pulsarParams.posepoch == 0. && pulsarParams.pepoch != 0.)
	pulsarParams.posepoch = pulsarParams.pepoch;
    fprintf(stderr,"Check on read tempo file, pepoch: %f\n", pulsarParams.pepoch);

    /*input.filename=psrEphemeris;*/ /*/archive/home/colingill/lalsuite/lalapps/src/pulsar/fscan/ /archive/home/colingill/public_html/crab_ephemeris.txt*/
    input.filename = XLALMalloc(sizeof(CHAR)*256);
    if(input.filename == NULL) fprintf(stderr,"input.filename pointer memory not allocated\t");

    strcpy(input.filename,psrEphemeris);
    /*fprintf(stderr,"psrEphemeris:%s\n", input.filename);*/

    /*The first stage is to read in f and fdot from the ephemeris, the crab_ephemeris.txt file is part of lalapps and is maintained by matt*/
    LALGetCrabEphemeris( &status, &crabEphemerisData, &input );
    /*check on the outputs, crabEphemerisData is a struct of type CrabSpindownParamsInput, and has members tArr, f0, f1*/
    fprintf(stderr,"input crab ephemeris present, number of entries: %i\n", crabEphemerisData.numOfData);
    fprintf(stderr,"crabEphemerisData: \ttarr= %f\tf0= %f\tf_dot= %e\n", crabEphemerisData.tArr->data[0], crabEphemerisData.f0->data[0], crabEphemerisData.f1->data[0]);
    
    /*Now I have f and fdot, use function below to compute the higher order derrivatives of the crabs frequency*/
    LALComputeFreqDerivatives( &status, &crabOutput, &crabEphemerisData );
    /*check on this function, crabOutput is datatype CrabSpindownParamsOutput*/
    fprintf(stderr,"crabOutput:\tf0= %f\tf1= %e\tf2= %e\n", crabOutput.f0->data[0], crabOutput.f1->data[0], crabOutput.f2->data[0]);/*need to change the type of printf for F1 and f2 to exponentials.*/
    fprintf(stderr,"--------------------\n");

    /*Allocate memory for edat, no need to do in the loop*/
    edat = XLALMalloc(sizeof(*edat));
    (*edat).ephiles.earthEphemeris = earthFile; /*"/archive/home/colingill/lalsuite/lal/packages/pulsar/test/earth05-09.dat";*/
    (*edat).ephiles.sunEphemeris = sunFile;/*"/archive/home/colingill/lalsuite/lal/packages/pulsar/test/sun05-09.dat";*/

    /*work out the eclitptic lattitude of the source, used in max freq calc.*/
    ecliptic_lat=asin(cos(e_tilt)*sin(pulsarParams.dec)-sin(pulsarParams.ra)*cos(pulsarParams.dec)*sin(e_tilt));
    fprintf(stderr,"eqcliptic_lat: %e\t", ecliptic_lat);

    for (i=0;i<l;i++)
    {
        if (i == (l-1))/*catches the last iteration where there is no i+1 entry in timestamps*/
        {
            cur_epoch = timestamps->data[i]+(timebaseline/2);
        }
        else
        {
            cur_epoch = timestamps->data[i]+((timestamps->data[i+1] - timestamps->data[i])/2);
        }

        fprintf(stderr,"cur_epoch: %d\t", cur_epoch);
        /*The time has to be set so that the nearest entry to that time in the ephemeris can be used*/
        dataEpoch.gpsSeconds = cur_epoch; /*INT8)floor(time->data[j]);*/
        dataEpoch.gpsNanoSeconds = 0;
    
        /*prepare hetParams, which is then used to get the freq derrivatives out and also in the next sub-section for Bary functions*/
        LALSetSpindownParams( &status, &hetParams, &crabOutput, dataEpoch );
        fprintf(stderr,"hetparams epoch: %f\t f0= %f\tf1= %e\n", hetParams.epoch, hetParams.f0, hetParams.f1);
    
        /*----------------------------------------------------------------------------------------------------------------*/
        /*    ---2---   */
        /*Add corrections for timing noise to get a better guess at the freq*/
        /*----------------------------------------------------------------------------------------------------------------*/
    
        /*now I want to use these params to calc the freq at any point in time.  Using the higher order derrivatives is how we adjust for timing noise.*/
    
        /*Get the time difference between the current epoch and the epoch of the ephemeris entry*/
        tdt= cur_epoch - hetParams.epoch;
        fprintf(stderr,"dt: %f,\tf1= %e,\tf2= %e,\tf3= %e,\tf4=%e\n", tdt, hetParams.f1, hetParams.f2, hetParams.f3, hetParams.f4);
    
        freq = 2.0*( hetParams.f0 + ((hetParams.f1)*tdt) + (((hetParams.f2)*tdt*tdt)/2) + (((hetParams.f3)*tdt*tdt*tdt)/6) + (((hetParams.f4)*tdt*tdt*tdt*tdt)/24) );
        fprintf(stderr,"crab fcoarse: %f\t", freq);
    
        /*freq = 2.0*(params->f0 + params->f1*t1 + (params->f2)*t1*t1+ (params->f3)*t1*t1*t1 + (params->f4)*t1*t1*t1*t1);*/  /*cg;line 486 from hetcrabpulsar, works out freq, this is with one order of t removed for each of the derrivatives of f, and also each term is divided by a factorial, 1!, 2!, 3!, but this starts one term along from the oroginal code as we have integrated the orginal code to get freq not phase*/
    
        /*----------------------------------------------------------------------------------------------------------------*/
        /*    ---3---   */
        /*Add doppler shifts for earth's motion back onto the freq to get actual observed freq at detectors.*/
        /*----------------------------------------------------------------------------------------------------------------*/
        
        /*now I have the freq, I need to add the doppler shift for the earths motion around the sun */
        baryinput.dInv = 0.;/*I can always set this to zero, as Matt said so, I must ask him why*/

        LAL_CALL( LALInitBarycenter(&status, edat), &status );/*  */
    
        /*this lines take position of detector which are in xyz coords in meters from earths centre and converts them into seconds (time)*/
        baryinput.site.location[0] = det.location[0]/LAL_C_SI;
        baryinput.site.location[1] = det.location[1]/LAL_C_SI;
        baryinput.site.location[2] = det.location[2]/LAL_C_SI;
    
        /*dtpos should be the time between the entry in the ephemeris and the point in time for which doppler shifts are being calc.ed*/
        dtpos = cur_epoch - pulsarParams.posepoch;
    
        /* set up RA, DEC, and distance variables for LALBarycenter*/
        baryinput.delta = pulsarParams.dec + dtpos*pulsarParams.pmdec;
        baryinput.alpha = pulsarParams.ra + dtpos*pulsarParams.pmra/cos(baryinput.delta);
        
        t2=cur_epoch+1;
    
        baryinput2 = baryinput;
        
        baryinput.tgps.gpsSeconds = (INT4)floor(cur_epoch);
        baryinput.tgps.gpsNanoSeconds = (INT4)floor((fmod(cur_epoch,1.0)*1.e9));
    
        baryinput2.tgps.gpsSeconds = (INT4)floor(t2);
        baryinput2.tgps.gpsNanoSeconds = (INT4)floor((fmod(t2,1.0)*1.e9));
    
        /*the barycentre functions are needed to calc the inputs for the correction to fcoarse, namely emit, earth and baryinput*/
        LAL_CALL( LALBarycenterEarth(&status, &earth, &baryinput.tgps, edat), &status );
        LAL_CALL( LALBarycenter(&status, &emit, &baryinput, &earth), &status );
        
        LAL_CALL( LALBarycenterEarth(&status, &earth2, &baryinput2.tgps, edat), &status );
        LAL_CALL( LALBarycenter(&status, &emit2, &baryinput2, &earth2), &status );
    
        /* I need to calc the correction to the freq for the doppler shifts, the correction is df, from line 1074 heterdyne_pulsar.  deltaT is T_emission in TDB - T_arrrival in GPS + light travel time to SSB.  we are working out delta(delatT) over 1 second, so do not bother with the divide by one bit.*/
        df = freq*(emit2.deltaT - emit.deltaT);
        fprintf(stderr,"df: %f,\t", df);

        /*Calc maximum possible df and then subtract applied df from it, use the remaining df to plot range over which f will wander.*/
        max_df=freq*(( (29.783e3*cos(ecliptic_lat))+(465*cos(pulsarParams.dec)) )/3.0e8);/*max doppler shift, from speed of earth in plane of direction to source over c.*/
        fprintf(stderr,"max df: %f\tbeta: %e\n", max_df, ecliptic_lat);
        maxf=freq+max_df;
        minf=freq-max_df;

        finalFreq=freq+df;
        /*df = fcoarse*(emit2.deltaT - emit.deltaT + binOutput2.deltaT - binOutput.deltaT);*//*use when have binary calcs in here also.*/
        fprintf(fp5,"%f\t%f\t%f\t%f\t%f\n", freq, df, minf, finalFreq, maxf);
        fprintf(stderr,"crab freq calc, i:%d,  minf:%f,  f:%f,  maxf:%f\n", i, minf, finalFreq, maxf);

        /*----------------------------------------------------------------------------------------------------------------*/
    }
    fclose(fp5);
    fprintf(stderr,"end of crab stuff\n");
    
    /*Free up any memory allocated in crab section*/
    XLALFree(edat);
    
    }
    #endif

    /*fprintf(stderr,"end of spec_avg 1\n");*/

    /*=======================================================================================================================*/
    /*=======================================================================================================================*/


    /*release a;; the allocaeted memory*/
    LALCHARDestroyVector(&status, &timestamp);
    LALCHARDestroyVector(&status, &year_date);
    LALDestroySFTVector (&status, &sft_vect );

    /*fprintf(stderr,"end of spec_avg 2\n");*/

    if (timeavg != NULL) XLALFree(timeavg);

    /*fprintf(stderr,"end of spec_avg 3\n");*/

    LAL_CALL(LALDestroyUserVars(&status), &status);

    /*fprintf(stderr,"end of spec_avg 4\n");*/
    /*close all the files, spec_avg.c is done, all info written to the files.*/
    fclose(fp);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);

    fprintf(stderr,"end of spec_avg\n");

    return(0);


}
/* END main */
