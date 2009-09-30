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


/*temporary rubbish bin for headers*/
/*These are included in HeterodyneCrabPulsar files
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/BinaryPulsarTiming.h>*/


/*end of temporary rubbish bin*/



#include <lalapps.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdio.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/NormalizeSFTRngMed.h>

/*LAL malloc header file*/
#include <lal/LALMalloc.h>

/* Not sure if I need these
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SkyCoordinates.h> 
#include <lal/DetectorSite.h>
#include <lal/SFTutils.h>
#include <lal/LALString.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>
#include <lal/LALRCSID.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>*/

/*normal c header files*/
/*test if these two headers are problem*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <lal/Date.h>/*cg; needed to use lal routine GPStoUTC, which is used to convert GPS seconds into UTC date*/

#define NUM 1000 /* used for defining structures such as crabOutput*/

RCSID("$Id");

int main(int argc, char **argv)
{
    FILE *fp  = NULL;
    FILE *fp2 = NULL;
    FILE *fp3 = NULL;
    FILE *fp4 = NULL;
    FILE *fp5 = NULL;
    FILE *fp6 = NULL;
    LALStatus status = blank_status;
    
    SFTCatalog *catalog = NULL;
    SFTVector *sft_vect = NULL;
    INT4 i,j,k,l;
    INT4 numBins, nSFT;
    SFTConstraints constraints=empty_SFTConstraints;
    LIGOTimeGPS startTime, endTime; 
    REAL8 avg =0;
    REAL4 *timeavg =NULL;
    REAL8 f =0;
    CHAR outbase[256],outfile[256],outfile2[256],outfile3[256], outfile4[256], outfile5[256], outfile6[256];
    REAL8 NumBinsAvg =0;
    REAL8 timebaseline =0;
    /*CHAR *missing_sft;*/
    
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
    LALUnixDate       date;
    /*LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;*/
    /*LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;*/
    CHARVector        *timestamp = NULL;
    CHARVector	     *year_date = NULL; /*cg; creates a char vector*/
    REAL8Vector     *timestamps=NULL;
    
    CHAR *psrInput = NULL;
    CHAR *psrEphemeris = NULL;
  /*========================================================================================================================*/



lalDebugLevel = 0;
LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);

LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",         'h', UVAR_HELP,     "Print this help message",     &help        ), &status);
LAL_CALL(LALRegisterSTRINGUserVar(&status, "SFTs",         'I', UVAR_REQUIRED, "SFT location/pattern",        &SFTpatt     ), &status);
LAL_CALL(LALRegisterSTRINGUserVar(&status, "IFO",          'I', UVAR_REQUIRED, "Detector",                    &IFO         ), &status);
LAL_CALL(LALRegisterINTUserVar   (&status, "startGPS",     's', UVAR_REQUIRED, "Starting GPS time",           &startGPS    ), &status);
LAL_CALL(LALRegisterINTUserVar   (&status, "endGPS",       'e', UVAR_REQUIRED, "Ending GPS time",             &endGPS      ), &status);
LAL_CALL(LALRegisterREALUserVar  (&status, "fMin",         'f', UVAR_REQUIRED, "Minimum frequency",           &f_min       ), &status);
LAL_CALL(LALRegisterREALUserVar  (&status, "fMax",         'F', UVAR_REQUIRED, "Maximum frequency",           &f_max       ), &status);
LAL_CALL(LALRegisterINTUserVar   (&status, "blocksRngMed", 'w', UVAR_OPTIONAL, "Running Median window size",  &blocksRngMed), &status);
LAL_CALL(LALRegisterSTRINGUserVar(&status, "outputBname",  'o', UVAR_OPTIONAL, "Base name of output files",   &outputBname ), &status);
LAL_CALL(LALRegisterREALUserVar  (&status, "freqRes",      'r', UVAR_REQUIRED, "Spectrogram freq resolution", &freqres     ), &status);
LAL_CALL(LALRegisterREALUserVar  (&status, "timeBaseline", 't', UVAR_REQUIRED, "The time baseline of sfts",   &timebaseline), &status);
LAL_CALL(LALRegisterSTRINGUserVar(&status, "psrInput",  'S', UVAR_OPTIONAL, "name of tempo pulsar file",   &psrInput ), &status);
LAL_CALL(LALRegisterSTRINGUserVar(&status, "psrEphemeris",  'S', UVAR_OPTIONAL, "pulsar ephemeris file",   &psrEphemeris ), &status);

LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);
if (help)
return(0);

startTime.gpsSeconds = startGPS;/*cg; startTime is a structure, and gpsSeconds is a member of that structure*/
startTime.gpsNanoSeconds = 0;/*cg; gps NanoSeconds is also a member of the startTime structure */
constraints.startTime = &startTime; /*cg; apersand operator gets the address of a variable, &a is a pointer to a.  This line puts the startTime structure into the structure constraints*/

endTime.gpsSeconds = endGPS;
endTime.gpsNanoSeconds = 0;
constraints.endTime = &endTime;/*cg; This line puts the end time into the structure constraints*/
constraints.detector = IFO;/*cg; this adds the interferometer into the contraints structure*/
LALSFTdataFind ( &status, &catalog,SFTpatt, &constraints );/*cg; creates SFT catalog, uses the constraints structure*/
LALLoadSFTs ( &status, &sft_vect, catalog, f_min,f_max);/*cg;reads the SFT data into the structure sft_vect*/
LALDestroySFTCatalog( &status, &catalog);/*cg; obvisouly desctroys the SFT catalogue*/

numBins = sft_vect->data->data->length;/*the number of bins in the freq_range*/
nSFT = sft_vect->length;/* the number of sfts.*/

fprintf(stderr, "nSFT = %d\tnumBins = %d\tf0 = %f\n", nSFT, numBins,sft_vect->data->f0);/*cg; prints out these items into the file called spectrumAverage_testcg_0.err which is written to the logs folder.*/
if (LALUserVarWasSet(&outputBname))
strcpy(outbase, outputBname);
else
sprintf(outbase, "spec_%.2f_%.2f_%s_%d_%d", f_min,f_max,constraints.detector,startTime.gpsSeconds,endTime.gpsSeconds);/*cg; this is the default name for producing the output files, the different suffixes are just added to this*/
sprintf(outfile,  "%s", outbase);/*cg; name of first file to be output*/
sprintf(outfile2, "%s_timestamps", outbase);/*cg: name of second file to be output*/
sprintf(outfile3, "%s_timeaverage", outbase);/*cg; name of third file to be output*/

fp = fopen(outfile, "w");/*cg;  open all three files for writing, if they don't exist create them, if they do exist overwrite them*/
fp2 = fopen(outfile2, "w");
fp3 = fopen(outfile3, "w");


/*----------------------------------------------------------------------------------------------------------------*/
/*cg;  Create the first file, called    blah_blah_blah*/
/*cg;  This file outputs very small numbers, one set of 100 numbers for each 1800 second SFT.  The numbers are the strain values used in the spectrogram plots.  The number of numbers here depends on the frequency range, freq range for this example is from 50Hz to 60Hz, so each Hz is split into 5, giving a freq resolution of 0.2Hz.*/


/*cg;  Create the second file, called    blah_b;ah_blah_timestamps*/
/*cg;  This will simply contain the time in GPS seconds of each SFT.  The time is the name of the SFT, which I am unsure if this refers to the start time, end time, or middle, but I think it is the start time.  So this file will conatin a list of the items, starting with the last SFT first labeled 0.  SFTno up to 2. SFTno for three SFTs, or more for more obvisouly.*/

  NumBinsAvg = freqres*numBins/(f_max-f_min);/*this calcs the number of bins over which to average the sft data, this is worked out so it produces the same frequency resoltuion as specified in the arguments passed to the python script. numBins is the total number of bins in the raw sft data*/
  
  /*startgps+(timebaseline*nSFT) will work out the time of that nSFT would start at IF there were no gaps in the data, where this is the case, replace the sft data with zeros.*/

  l=0;
  /*create output files and check for missing sfts*/
  for (j=0;j<nSFT;j++)/*cg;nSFT is the numnber of SFT files used for the time specified. So process is repeated for each SFT*/
  { 
    /*now I need to check if timebaseline+(n-1)sft_epoch = this sft_epoch, if not I need to create a element in the matrix of sft data that starts at timebaseline+(n-1)sft_epoch*/
    cur_epoch = sft_vect->data[j].epoch.gpsSeconds;/*finds the gps time of the current sft in the sequence with indicie j*/
    fprintf(fp2, "%d.\t%d\n", l, cur_epoch);/*cg; this bit writes the second file, i.e. the timestamps*/
  
    for ( i=0; i < (numBins-2); i+=NumBinsAvg)/*cg; this loop works out the powers and writes the first file. += makes i increment by NumBinsAvg.*/
    {/*cg; each SFT is split up into a number of bins, the number of bins is read in from the SFT file*/
        avg = 0.0;/*cg; the vairable avg is reset each time.*/
        if (i+NumBinsAvg>numBins) {printf("Error\n");return(2);}/*cg; error is detected, to prevent referencing data past the end of sft_vect.*/
        for (k=0;k<NumBinsAvg;k++)/*cg; for each bin, k goes trhough each entry from 0 to 180.*/
            avg += sqrt(sft_vect->data[j].data->data[i+k].re*sft_vect->data[j].data->data[i+k].re + 
            sft_vect->data[j].data->data[i+k].im*sft_vect->data[j].data->data[i+k].im);/*cg; re amd im are real and imaginary parts of SFT, duh!*/
        fprintf(fp,"%e\t",avg/NumBinsAvg);
    }
    fprintf(fp,"\n");
    /*------------------------------*/
    /*Bit to check if there is a gap in the sfts*/
    if ( j < (nSFT-1) )/*in all cases except when we are examining the last sft, check that there is no gap to the next sft*/
    {
      next_epoch = sft_vect->data[j+1].epoch.gpsSeconds;
      
      if (cur_epoch+timebaseline != next_epoch )/*if this returns true then the finishing time of current sft does not match starttime of the next sft.*/
      {
	  l=l+1;
          cur_epoch = cur_epoch+timebaseline;
	  fprintf(fp2, "%d.\t%d\n", l, cur_epoch );
	  
	  /*fill in blanks into the specgram matrix passed to matlab for the gap between the non con-current sfts*/
	  for ( i=0; i < (numBins-2); i+=NumBinsAvg)
	  {
	      avg = 0.0;
	      if (i+NumBinsAvg>numBins) {printf("Error\n");return(2);}
	      fprintf(fp,"%e\t",avg);
	  }
	  fprintf(fp,"\n");
          /*test to see if SFT gap is longer than 3/2*timebaseline, if so create another entry in the matrix*/
          
          while ( (cur_epoch + (0.5*timebaseline)) < next_epoch ) 
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
          }
      }
      
    }
    /*----------------------------*/
  l=l+1;
  }
fprintf(stderr,"finished checking for missing sfts, l=%d\n", l);
/*----------------------------------------------------------------------------------------------------------------*/
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
	    timeavg[i] = sft_vect->data[j].data->data[i].re*sft_vect->data[j].data->data[i].re + 
			sft_vect->data[j].data->data[i].im*sft_vect->data[j].data->data[i].im;
	  } 
	  else 
	  {
	    timeavg[i] += sft_vect->data[j].data->data[i].re*sft_vect->data[j].data->data[i].re + 
			sft_vect->data[j].data->data[i].im*sft_vect->data[j].data->data[i].im;
	  }
      }
  }
/*----------------------------------------------------------------------------------------------------------------*/
/*cg;  Create the third and final file, called   blah_blah_blah_timeaverage.  This file will contain the data used in the matlab plot script to plot the normalised average power vs the frequency.*/

/* Question to the group, the timeaverages are so far output for every freq bin in the raw data, independant of freq res chosen for the spectrogram, do we want this???*/
/*timeavg records the power of each bin*/


for ( i=0; i < numBins; i++)
    {
    f = sft_vect->data->f0 + ((REAL4)i)*sft_vect->data->deltaF;
    fprintf(fp3,"%16.8f %g\n",f,timeavg[i]/((REAL4)nSFT));
    } 

/*------------------------------------------------------------------------------------------------------------------------*/
/* CG; create a file containing the UTC times for each SFT GPS time.*/
sprintf(outfile4, "%s_date", outbase);/*cg; sprintf prints into the char array outfile4, it prints a string followed by _date, the string is outbase.*/
fp4 = fopen(outfile4, "w");

sprintf(outfile6, "%s_numbins", outbase);
fp6 = fopen(outfile6, "w");
fprintf(fp6,"%f\n",NumBinsAvg);
fclose(fp6);

LALCHARCreateVector(&status, &timestamp, (UINT4)256); /*128 is the number of elements in the vector, no need for it to be this big.*/
LALCHARCreateVector(&status, &year_date, (UINT4)128); 
timestamps = XLALCreateREAL8Vector(l+1);
/*times = XLALCreateREAL8Vector( MAXLENGTH );*/

l=0;

/*output the UTC time and date for use by matlab plotting script.*/
for (j=0;j<nSFT;j++)
{
    /*fprintf(stderr,"here 2a, j=%d\n",j);*/
    cur_epoch = sft_vect->data[j].epoch.gpsSeconds;
    /*gps.gpsSeconds = cur_epoch;*//*cg;stores GPS time in seconds into the structure gps, in the field gpsSeconds, for use in LALGPStoUTC function*/
    /*gps.gpsNanoSeconds=0;*//*cg; this sets the gps nanoseconds to 0*/

    /*timestamps->data[l]= gps.gpsSeconds;*/
    timestamps->data[l]= cur_epoch;

    XLALGPSToUTC(&date, cur_epoch);/*cg; gets the UTC date in struct tm format from the GPS seconds.*/
    fprintf(fp4, "%d\t %i\t %i\t %i\t %i\t %i\t %i\n", l, (date.tm_year+1900), (date.tm_mon + 1), date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec);
    /*------------------------------*/
    /*Bit to check if there is a gap in the sfts*/
    if ( j< (nSFT-1) )/*in all cases except when we are examining the last sft, check that there is no gap to the next sft*/
    {
	next_epoch = sft_vect->data[j+1].epoch.gpsSeconds;/*p->x, if p is a pointer then this expression gets the member of that structure. The datat subscript starts at zero and runs to nSFT-1, there is no data[nSFT]!!!  */
	/*cur_epoch = sft_vect->data[j].epoch.gpsSeconds;*/
	if (cur_epoch+timebaseline != next_epoch )/*if this returns true then the finishing time of current sft does not match starttime of the next sft.*/
	{
	    l=l+1;
	    cur_epoch = cur_epoch+timebaseline;/*calcs the timestamp to be used for the gap, this is the end of the last sft before the gap*/

	    timestamps->data[l]= cur_epoch;

	    XLALGPSToUTC(&date, cur_epoch);
	    fprintf(fp4, "%d\t %d\t %d\t %d\t %d\t %d\t %d\n", l, (date.tm_year+1900), (date.tm_mon + 1), date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec);

            while ( (cur_epoch + (0.5*timebaseline)) < next_epoch ) 
            {
                l=l+1;
              
                XLALGPSToUTC(&date, cur_epoch);
                fprintf(fp4, "%d\t %d\t %d\t %d\t %d\t %d\t %d\n", l, (date.tm_year+1900), (date.tm_mon + 1), date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec);
        
                timestamps->data[l]= cur_epoch;
                cur_epoch=cur_epoch+timebaseline;/*increment cur_epoch to test if there is still a long gap that needs spliiting*/
            }
	}
    }
	    /*----------------------------*/
  l=l+1;
  /*fprintf(stderr,"here 2c, j=%d\n",j);*/
  }
  fprintf(stderr,"number of timestamps: %d\n",l);
  /*fprintf(fp4, "%u\n", year_date->length);*/ /*just used this as a test line*/
  /*fprintf(fp4, "%u\n", year_date->length);*/ /*just used this as a test line*/
    
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
    #include<HeterodyneCrabPulsar.h>

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
    REAL8 df=0., freq, finalFreq;
    REAL8 dtpos=0.; /* time between position epoch and data timestamp */
    BinaryPulsarParams pulsarParams; /*general pulsar params strcut, despite binary name*/
    /*CHAR *psrInput = NULL;*/ /* pulsar input file containing params f0, f1 etc. */
    LALDetector det;
    CHAR detName[256];
    
    /*----------------------------------------------------------------------------------------------------------------*/
    /*for debuging */
    sprintf(outfile5, "%s_crab", outbase);
    fp5 = fopen(outfile5, "w");
    
    /*cg; calculating the crab freq, done in three steps*/
    /*this needs to be done for each sft, so need some sort of iteration*/
    
    /*----------------------------------------------------------------------------------------------------------------*/
    /*    ---1---   */
    /*Find rough guess of crabs freq from ephemeris*/
    /*----------------------------------------------------------------------------------------------------------------*/
    
    /*Get detector position, this is needed for barycentre calcs*/
    det = *XLALGetSiteInfo( IFO );

    /* read in tempo par file for pulsar, This is one of the optional command line arguments*/
    /*fprintf(stderr,"%s\n",psrInput);*/
    
    /*psrInput="../B0531+21";*//*../B0531+21" psrInput is now read in from the command line arguments and passed to this code*/
    /*fprintf(stderr,"%s\n",psrInput);*/
    XLALReadTEMPOParFile(&pulsarParams, psrInput);

    /*Make sure that posepoch and pepoch are set*/
    if(pulsarParams.pepoch == 0. && pulsarParams.posepoch != 0.)
	pulsarParams.pepoch = pulsarParams.posepoch;
    else if(pulsarParams.posepoch == 0. && pulsarParams.pepoch != 0.)
	pulsarParams.posepoch = pulsarParams.pepoch;
    /*fprintf(stderr,"Check on read tempo file, pepoch: %f\n", pulsarParams.pepoch);*/

    /*input.filename=psrEphemeris;*/ /*/archive/home/colingill/lalsuite/lalapps/src/pulsar/fscan/ /archive/home/colingill/public_html/crab_ephemeris.txt*/
    input.filename = XLALMalloc(sizeof(CHAR)*256);
    if(input.filename == NULL) fprintf(stderr,"input.filename pointer memory not allocated\t");

    strcpy(input.filename,psrEphemeris);
    /*fprintf(stderr,"psrEphemeris:%s\n", input.filename);*/

    /*The first stage is to read in f and fdot from the ephemeris, the ephemeris file is part of lalapps and is maintained by matt*/
    LALGetCrabEphemeris( &status, &crabEphemerisData, &input );
    /*check on the oputputs, crabEphemerisData is a struct of type CrabSpindownParamsInput, and has members tArr, f0, f1*/
    /*fprintf(stderr,"input crab ephemeris present, number of entries: %i\n", crabEphemerisData.numOfData);
    fprintf(stderr,"crabEphemerisData:\t%f\t%f\t%f\n", crabEphemerisData.tArr->data[0], crabEphemerisData.f0->data[0], crabEphemerisData.f1->data[0]);*/
    
    /*Now I have f and fdot, use function below to compute the higher order derrivatives of the crabs frequency*/
    LALComputeFreqDerivatives( &status, &crabOutput, &crabEphemerisData );
    /*check on this function, crabOutput is datatype CrabSpindownParamsOutput*/
    /*fprintf(stderr,"crabOutput:\t%f\t%f\t%f\n", crabOutput.f0->data[0], crabOutput.f1->data[0], crabOutput.f2->data[0]);*/

    /*Allocate memory for edat, no need to do in the loop*/
    edat = XLALMalloc(sizeof(*edat));
    (*edat).ephiles.earthEphemeris =  "/archive/home/colingill/lalsuite/lalapps/src/pulsar/fscan/earth05-09.dat"; 
    (*edat).ephiles.sunEphemeris = "/archive/home/colingill/lalsuite/lalapps/src/pulsar/fscan/sun05-09.dat";

    
    for (i=0;i<l;i++)
    {
        /*fprintf(stderr,"iteration number:  %d\t", i);*/
        if (i == (l-1))/*catches the last iteration where there is no i+1 entry in timestamps*/
        {
            cur_epoch = timestamps->data[i]+(timebaseline/2);
            /*fprintf(stderr,"EXCEPTION CAUGHT 1\t");*/
        }
        else
        {
            cur_epoch = timestamps->data[i]+((timestamps->data[i+1] - timestamps->data[i])/2);
        }

        /*fprintf(stderr,"cur_epoch:  %d\t", cur_epoch);*/
        /*The time has to be set so that the nearest entry to that time in the ephemeris can be used*/
        dataEpoch.gpsSeconds = cur_epoch; /*INT8)floor(time->data[j]);*/
        dataEpoch.gpsNanoSeconds = 0;
    
        /*prepare hetParams, which is then used to get the freq derrivatives out and also in the next sub-section for Bary functions*/
        LALSetSpindownParams( &status, &hetParams, &crabOutput, dataEpoch );
        /*fprintf(fp5,"hetparams:\t%f\t%f\t%f\n", hetParams.epoch, hetParams.f0, hetParams.f1);*/
    
        /*----------------------------------------------------------------------------------------------------------------*/
        /*    ---2---   */
        /*Add corrections for timing noise to get a better guess at the freq*/
        /*----------------------------------------------------------------------------------------------------------------*/
    
        /*now I want to use these params to calc the freq at any point in time.  Using the higher order derrivatives is how we adjust for timing noise.*/
    
        /*Get the time difference between the current epoch and the epoch of the ephemeris entry*/
        tdt= cur_epoch - hetParams.epoch;
        /*fprintf(fp5,"current epoch: %i\tephemeris epoch: %f\t delta t: %f\n", cur_epoch, hetParams.epoch, tdt);*/
    
        freq = 2.0*( hetParams.f0 + ((hetParams.f1)*tdt) + (((hetParams.f2)*tdt*tdt)/2) + (((hetParams.f3)*tdt*tdt*tdt)/6) + (((hetParams.f4)*tdt*tdt*tdt*tdt)/24) );
        /*fprintf(fp5,"crab fcoarse: %f\n", freq);*/
    
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
        
        /* set leap seconds noting that for all runs prior to S5 that the number
        of leap seconds was 13, 1 leap seconds was added on 31st
        Dec 2005 24:00:00 i.e. GPS 820108813, and another was added on 31st
        Dec 2008 24:00:00 i.e. GPS 91480321 */
        if(cur_epoch <= 820108813)
            (*edat).leap = 13;
        else if(cur_epoch <= 914803214)
            (*edat).leap = 14;
        else
            (*edat).leap = 15;
        
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
    
        /* I need to calc the correction to the freq for the doppler shifts, the correction is df, from line 1074 heterdyne_pulsar.  Note this correction is a gradient.*/
        df = freq*(emit2.deltaT - emit.deltaT);
        finalFreq=freq+df;
        /*df = fcoarse*(emit2.deltaT - emit.deltaT + binOutput2.deltaT - binOutput.deltaT);*//*use when have binary calcs in here also.*/
        fprintf(fp5,"%f\t%f\t%f\n", freq, df, finalFreq);
        /*fprintf(stderr,"crab freq calc, i: %d,  f:  %f\n", i, finalFreq);*/

        /*----------------------------------------------------------------------------------------------------------------*/
    }
    fclose(fp5);
    fprintf(stderr,"end of crab stuff\n");
    
    /*Free up any memory allocated in crab section*/
    XLALFree(edat);
    
    #endif

    fprintf(stderr,"end of spec_avg 1\n");

    /*=======================================================================================================================*/
    /*=======================================================================================================================*/


    /*release a;; the allocaeted memory*/
    LALCHARDestroyVector(&status, &timestamp);
    LALCHARDestroyVector(&status, &year_date);
    LALDestroySFTVector (&status, &sft_vect );

    fprintf(stderr,"end of spec_avg 2\n");

    if (timeavg != NULL) XLALFree(timeavg);

    fprintf(stderr,"end of spec_avg 3\n");

    LAL_CALL(LALDestroyUserVars(&status), &status);

    fprintf(stderr,"end of spec_avg 4\n");
    /*close all the files, spec_avg.c is done, all info written to the files.*/
    fclose(fp);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);

    fprintf(stderr,"end of spec_avg 5\n");
    
    return(0);


}
/* END main */
