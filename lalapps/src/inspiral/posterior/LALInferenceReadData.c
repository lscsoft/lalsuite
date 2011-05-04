/* 
 *  LALInferenceReadData.c:  Bayesian Followup functions
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
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

#include <stdio.h>
#include <stdlib.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>

#include <lal/LALInspiral.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALDetectors.h>
#include <lal/AVFactories.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/VectorOps.h>
#include <lal/Random.h>
#include <lal/LALNoiseModels.h>
#include <lal/XLALError.h>
#include <lal/GenerateInspiral.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLInspiralRead.h>


#include "LALInference.h"

const LALUnit strainPerCount={0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length);


REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length)
{
	LALStatus status;
	memset(&status,0,sizeof(status));
	FrCache *cache = NULL;
	FrStream *stream = NULL;
	REAL8TimeSeries *out = NULL;
	
	cache  = XLALFrImportCache( cachefile );
        int err;
        err = *XLALGetErrnoPtr();
	if(cache==NULL) {fprintf(stderr,"ERROR: Unable to import cache file \"%s\",\n       XLALError: \"%s\".\n",cachefile, XLALErrorString(err)); exit(-1);}
	stream = XLALFrCacheOpen( cache );
	if(stream==NULL) {fprintf(stderr,"ERROR: Unable to open stream from frame cache file\n"); exit(-1);}
	out = XLALFrInputREAL8TimeSeries( stream, channel, &start, length , 0 );
	if(out==NULL) fprintf(stderr,"ERROR: unable to read channel %s from %s at time %i\nCheck the specified data duration is not too long\n",channel,cachefile,start.gpsSeconds);
	LALDestroyFrCache(&status,&cache);
	LALFrClose(&status,&stream);
	return out;
}
#define USAGE \
"Variables needed from command line to read data:\n\
[ --channel [channel1,channel2,channel3,..] ] \n\
--IFO [IFO1,IFO2,IFO3,..] \n\
--cache [cache1,cache2,cache3,..] \n\
  (Use LALLIGO, LAL2kLIGO, LALGEO, LALVirgo, LALAdLIGO to simulate these detectors)\n \
--PSDstart GPSsecs.GPSnanosecs \n\
--PSDlength length \n\
[--srate SampleRate   [4096]] \n\
--seglen segment_length \n\
--trigtime GPSsecs.GPSnanosecs \n\
[--fLow [cutoff1,cutoff2,cutoff3,..] [40Hz]] \n\
[--fHigh [fHigh1,fHigh2,fHigh3,..] [f_Nyquist]]\n\
[--dataseed number]\n"

LALIFOData *readData(ProcessParamsTable *commandLine)
/* Read in the data and store it in a LALIFOData structure */
{
	LALStatus status;
	INT4 dataseed=0;
	memset(&status,0,sizeof(status));
	ProcessParamsTable *procparam=NULL;
	LALIFOData *headIFO=NULL,*IFOdata=NULL;
	REAL8 SampleRate=4096.0,SegmentLength=0;
	if(getProcParamVal(commandLine,"--srate")) SampleRate=atof(getProcParamVal(commandLine,"--srate")->value);
        const REAL8 defaultFLow = 40.0;
        const REAL8 defaultFHigh = SampleRate/2.0;
	int nSegs=0;
	size_t seglen=0;
	REAL8TimeSeries *PSDtimeSeries=NULL;
	REAL8 padding=0.4;//Default was 1.0 second. However for The Event the Common Inputs specify a Tukey parameter of 0.1, so 0.4 second of padding for 8 seconds of data.
	UINT4 Ncache=0,Nifo=0,Nchannel=0,NfLow=0,NfHigh=0;
	UINT4 i,j;
	int FakeFlag=0;
	char strainname[]="LSC-STRAIN";
	
	typedef void (NoiseFunc)(LALStatus *status,REAL8 *psd,REAL8 f);
	NoiseFunc *PSD=NULL;
	REAL8 scalefactor=1;

	RandomParams *datarandparam;

	char *chartmp=NULL;
	char **channels=NULL;
	char **caches=NULL;
	char **IFOnames=NULL;
	char **fLows=NULL,**fHighs=NULL;
	LIGOTimeGPS GPSstart,GPStrig,segStart;
	REAL8 PSDdatalength=0;

	if(!getProcParamVal(commandLine,"--cache")||!getProcParamVal(commandLine,"--IFO")||
	   !getProcParamVal(commandLine,"--PSDstart")||!getProcParamVal(commandLine,"--trigtime")||
	   !getProcParamVal(commandLine,"--PSDlength")||!getProcParamVal(commandLine,"--seglen"))
	{fprintf(stderr,USAGE); return(NULL);}
	
  //TEMPORARY. JUST FOR CHECKING USING SPINSPIRAL PSD
  char **spinspiralPSD=NULL;
  UINT4 NspinspiralPSD = 0;
  if (getProcParamVal(commandLine, "--spinspiralPSD")) {
    parseCharacterOptionString(getProcParamVal(commandLine,"--spinspiralPSD")->value,&spinspiralPSD,&NspinspiralPSD);
  }    
  
	if(getProcParamVal(commandLine,"--channel")){
		parseCharacterOptionString(getProcParamVal(commandLine,"--channel")->value,&channels,&Nchannel);
	}
	parseCharacterOptionString(getProcParamVal(commandLine,"--cache")->value,&caches,&Ncache);
	parseCharacterOptionString(getProcParamVal(commandLine,"--IFO")->value,&IFOnames,&Nifo);
	if(getProcParamVal(commandLine,"--fLow")){
		parseCharacterOptionString(getProcParamVal(commandLine,"--fLow")->value,&fLows,&NfLow);
	}
	if(getProcParamVal(commandLine,"--fHigh")){
		parseCharacterOptionString(getProcParamVal(commandLine,"--fHigh")->value,&fHighs,&NfHigh);
	}
	if(getProcParamVal(commandLine,"--dataseed")){
		procparam=getProcParamVal(commandLine,"--dataseed");
		dataseed=atoi(procparam->value);
	}
								   
	if(Nifo!=Ncache) {fprintf(stderr,"ERROR: Must specify equal number of IFOs and Cache files\n"); exit(1);}
	if(Nchannel!=0 && Nchannel!=Nifo) {fprintf(stderr,"ERROR: Please specify a channel for all caches, or omit to use the defaults\n"); exit(1);}
	
	IFOdata=headIFO=calloc(sizeof(LALIFOData),Nifo);
	
	procparam=getProcParamVal(commandLine,"--PSDstart");
	LALStringToGPS(&status,&GPSstart,procparam->value,&chartmp);
	procparam=getProcParamVal(commandLine,"--trigtime");
	LALStringToGPS(&status,&GPStrig,procparam->value,&chartmp);
	PSDdatalength=atof(getProcParamVal(commandLine,"--PSDlength")->value);
	SegmentLength=atof(getProcParamVal(commandLine,"--seglen")->value);
	seglen=(size_t)(SegmentLength*SampleRate);
	nSegs=(int)floor(PSDdatalength/SegmentLength);
	
	for(i=0;i<Nifo;i++) {
          IFOdata[i].fLow=fLows?atof(fLows[i]):defaultFLow; 
          IFOdata[i].fHigh=fHighs?atof(fHighs[i]):defaultFHigh;
          strncpy(IFOdata[i].name, IFOnames[i], DETNAMELEN);
        }

	/* Only allocate this array if there weren't channels read in from the command line */
	if(!Nchannel) channels=calloc(Nifo,sizeof(char *));
	for(i=0;i<Nifo;i++) {
		if(!Nchannel) channels[i]=malloc(VARNAME_MAX);
		IFOdata[i].detector=calloc(1,sizeof(LALDetector));
		
		if(!strcmp(IFOnames[i],"H1")) {			
			memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
			if(!Nchannel) sprintf((channels[i]),"H1:%s",strainname); continue;}
		if(!strcmp(IFOnames[i],"H2")) {
			memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
			if(!Nchannel) sprintf((channels[i]),"H2:%s",strainname); continue;}
		if(!strcmp(IFOnames[i],"LLO")||!strcmp(IFOnames[i],"L1")) {
			memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexLLODIFF],sizeof(LALDetector));
			if(!Nchannel) sprintf((channels[i]),"L1:%s",strainname); continue;}
		if(!strcmp(IFOnames[i],"V1")||!strcmp(IFOnames[i],"VIRGO")) {
			memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexVIRGODIFF],sizeof(LALDetector));
			if(!Nchannel) sprintf((channels[i]),"V1:h_16384Hz"); continue;}
		if(!strcmp(IFOnames[i],"GEO")||!strcmp(IFOnames[i],"G1")) {
			memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexGEO600DIFF],sizeof(LALDetector));
			if(!Nchannel) sprintf((channels[i]),"G1:DER_DATA_H"); continue;}
		/*		if(!strcmp(IFOnames[i],"TAMA")||!strcmp(IFOnames[i],"T1")) {inputMCMC.detector[i]=&lalCachedDetectors[LALDetectorIndexTAMA300DIFF]; continue;}*/
		fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 H2 L1 V1 GEO\n",IFOnames[i]); exit(-1);
	}
	
	/* Set up FFT structures and window */
	for (i=0;i<Nifo;i++){
		/* Create FFT plans */
		IFOdata[i].timeToFreqFFTPlan = XLALCreateForwardREAL8FFTPlan((UINT4) seglen, 0 );
		IFOdata[i].freqToTimeFFTPlan = XLALCreateReverseREAL8FFTPlan((UINT4) seglen,0);
		
		/* Setup windows */
		IFOdata[i].window=XLALCreateTukeyREAL8Window(seglen,(REAL8)2.0*padding*SampleRate/(REAL8)seglen);
	}
	
	
	/* Trigger time = 2 seconds before end of segment (was 1 second, but Common Inputs for The Events are -6 +2*/
	memcpy(&segStart,&GPStrig,sizeof(LIGOTimeGPS));
	XLALGPSAdd(&segStart,-SegmentLength+2);
	
	
	/* Read the PSD data */
	for(i=0;i<Nifo;i++) {
		memcpy(&(IFOdata[i].epoch),&segStart,sizeof(LIGOTimeGPS));
		/* Check if fake data is requested */
		if(!(strcmp(caches[i],"LALLIGO") && strcmp(caches[i],"LALVirgo") && strcmp(caches[i],"LALGEO") && strcmp(caches[i],"LALEGO")
			 && strcmp(caches[i],"LALAdLIGO")))
		{
			FakeFlag=1;
			datarandparam=XLALCreateRandomParams(dataseed?dataseed+(int)i:dataseed);
			/* Selection of the noise curve */
			if(!strcmp(caches[i],"LALLIGO")) {PSD = &LALLIGOIPsd; scalefactor=9E-46;}
			if(!strcmp(caches[i],"LALVirgo")) {PSD = &LALVIRGOPsd; scalefactor=1.0;}
			if(!strcmp(caches[i],"LALGEO")) {PSD = &LALGEOPsd; scalefactor=1E-46;}
			if(!strcmp(caches[i],"LALEGO")) {PSD = &LALEGOPsd; scalefactor=1.0;}
			if(!strcmp(caches[i],"LALAdLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 10E-49;}
			//if(!strcmp(caches[i],"LAL2kLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 36E-46;}
			if(PSD==NULL) {fprintf(stderr,"Error: unknown simulated PSD: %s\n",caches[i]); exit(-1);}
			IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)
						XLALCreateREAL8FrequencySeries("spectrum",&GPSstart,0.0,
										   (REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
			for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
			{
				PSD(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF);
				IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]*=scalefactor;
			}
			IFOdata[i].freqData = (COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("stilde",&segStart,0.0,IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,&lalDimensionlessUnit,seglen/2 +1);
			/* Create the fake data */
			int j_Lo = (int) IFOdata[i].fLow/IFOdata[i].freqData->deltaF;
			for(j=j_Lo;j<IFOdata[i].freqData->data->length;j++){
				IFOdata[i].freqData->data->data[j].re=XLALNormalDeviate(datarandparam)*(0.5*sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]/IFOdata[i].freqData->deltaF));
				IFOdata[i].freqData->data->data[j].im=XLALNormalDeviate(datarandparam)*(0.5*sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]/IFOdata[i].freqData->deltaF));
			}
			IFOdata[i].freqData->data->data[0].re=0; 			IFOdata[i].freqData->data->data[0].im=0;
			const char timename[]="timeData";
			IFOdata[i].timeData=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries(timename,&segStart,0.0,(REAL8)1.0/SampleRate,&lalDimensionlessUnit,(size_t)seglen);
			XLALREAL8FreqTimeFFT(IFOdata[i].timeData,IFOdata[i].freqData,IFOdata[i].freqToTimeFFTPlan);
			if(*XLALGetErrnoPtr()) printf("XLErr: %s\n",XLALErrorString(*XLALGetErrnoPtr()));
			XLALDestroyRandomParams(datarandparam);

		}
		else{
			fprintf(stderr,"Estimating PSD for %s using %i segments of %i samples (%lfs)\n",IFOnames[i],nSegs,(int)seglen,SegmentLength);
			
			PSDtimeSeries=readTseries(caches[i],channels[i],GPSstart,PSDdatalength);
			if(!PSDtimeSeries) {fprintf(stderr,"Error reading PSD data for %s\n",IFOnames[i]); exit(1);}
			XLALResampleREAL8TimeSeries(PSDtimeSeries,1.0/SampleRate);
			PSDtimeSeries=(REAL8TimeSeries *)XLALShrinkREAL8TimeSeries(PSDtimeSeries,(size_t) 0, (size_t) seglen*nSegs);
			IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("spectrum",&PSDtimeSeries->epoch,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
			if (getProcParamVal(commandLine, "--PSDwelch")) {
        XLALREAL8AverageSpectrumWelch(IFOdata[i].oneSidedNoisePowerSpectrum ,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan);
      }
      else {
        XLALREAL8AverageSpectrumMedian(IFOdata[i].oneSidedNoisePowerSpectrum ,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan);	
			}
        XLALDestroyREAL8TimeSeries(PSDtimeSeries);
			
			/* Read the data segment */
			IFOdata[i].timeData=readTseries(caches[i],channels[i],segStart,SegmentLength);

                        /* FILE *out; */
                        /* char fileName[256]; */
                        /* snprintf(fileName, 256, "readTimeData-%d.dat", i); */
                        /* out = fopen(fileName, "w"); */
                        /* for (j = 0; j < IFOdata[i].timeData->data->length; j++) { */
                        /*   fprintf(out, "%g %g\n", j*IFOdata[i].timeData->deltaT, IFOdata[i].timeData->data->data[j]); */
                        /* } */
                        /* fclose(out); */
                        
			if(!IFOdata[i].timeData) {fprintf(stderr,"Error reading segment data for %s at %i\n",IFOnames[i],segStart.gpsSeconds); exit(1);}
			XLALResampleREAL8TimeSeries(IFOdata[i].timeData,1.0/SampleRate);	 
			if(!IFOdata[i].timeData) {fprintf(stderr,"Error reading segment data for %s\n",IFOnames[i]); exit(1);}
			IFOdata[i].freqData=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("freqData",&(IFOdata[i].timeData->epoch),0.0,1.0/SegmentLength,&lalDimensionlessUnit,seglen/2+1);
			IFOdata[i].windowedTimeData=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("windowed time data",&(IFOdata[i].timeData->epoch),0.0,1.0/SampleRate,&lalDimensionlessUnit,seglen);
			XLALDDVectorMultiply(IFOdata[i].windowedTimeData->data,IFOdata[i].timeData->data,IFOdata[i].window->data);
			XLALREAL8TimeFreqFFT(IFOdata[i].freqData,IFOdata[i].windowedTimeData,IFOdata[i].timeToFreqFFTPlan);
			
			for(j=0;j<IFOdata[i].freqData->data->length;j++){
				IFOdata[i].freqData->data->data[j].re/=sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
				IFOdata[i].freqData->data->data[j].im/=sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
                                IFOdata[i].windowedTimeData->data->data[j] /= sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
			}
			
		}
                /* Now that the PSD is set up, make the TDW. */
                IFOdata[i].timeDomainNoiseWeights = 
                  (REAL8TimeSeries *)XLALCreateREAL8TimeSeries("time domain weights", 
                                                               &(IFOdata[i].oneSidedNoisePowerSpectrum->epoch),
                                                               0.0,
                                                               1.0/SampleRate,
                                                               &lalDimensionlessUnit,
                                                               seglen);
                PSDToTDW(IFOdata[i].timeDomainNoiseWeights, IFOdata[i].oneSidedNoisePowerSpectrum, IFOdata[i].freqToTimeFFTPlan,
                         IFOdata[i].fLow, IFOdata[i].fHigh);

                makeWhiteData(&(IFOdata[i]));
    
    if (getProcParamVal(commandLine, "--spinspiralPSD")) {
      FILE *in;
      //char fileNameIn[256];
      //snprintf(fileNameIn, 256, spinspiralPSD);
      double freq_temp, psd_temp, temp;
      int n=0;
      int k=0;
      int templen=0;
      char buffer[256];
      char * line=buffer;
    
      //in = fopen(fileNameIn, "r");
      in = fopen(spinspiralPSD[i], "r");
      while(fgets(buffer, 256, in)){
        templen++;
      }
    
     // REAL8 *tempPSD = NULL;
     // REAL8 *tempfreq = NULL;
     // tempPSD=calloc(sizeof(REAL8),templen+1);
     // tempfreq=calloc(sizeof(REAL8),templen+1);
    
      rewind(in);
      IFOdata[i].oneSidedNoisePowerSpectrum->data->data[0] = 1.0;
      while(fgets(buffer, 256, in)){
        line=buffer;
      
        sscanf(line, "%lg%n", &freq_temp,&n);
        line+=n;
        sscanf(line, "%lg%n", &psd_temp,&n);
        line+=n;
        sscanf(line, "%lg%n", &temp,&n);
        line+=n;
      
     // tempfreq[k]=freq_temp;
     // tempPSD[k]=psd_temp*psd_temp;
        
        IFOdata[i].oneSidedNoisePowerSpectrum->data->data[k+1]=psd_temp*psd_temp;
        
      k++;
      //fprintf(stdout, "%g %g \n",freq_temp, psd_temp); fflush(stdout);
      }
      fclose(in);
    }
                
                if (getProcParamVal(commandLine, "--data-dump")) {
                  const UINT4 nameLength=256;
                  char filename[nameLength];
                  FILE *out;

                  snprintf(filename, nameLength, "%s-PSD.dat", IFOdata[i].name);
                  out = fopen(filename, "w");
                  for (j = 0; j < IFOdata[i].oneSidedNoisePowerSpectrum->data->length; j++) {
                    REAL8 f = IFOdata[i].oneSidedNoisePowerSpectrum->deltaF*j;
                    REAL8 psd = IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j];

                    fprintf(out, "%g %g\n", f, psd);
                  }
                  fclose(out);

                  snprintf(filename, nameLength, "%s-timeData.dat", IFOdata[i].name);
                  out = fopen(filename, "w");
                  for (j = 0; j < IFOdata[i].timeData->data->length; j++) {
                    REAL8 t = XLALGPSGetREAL8(&(IFOdata[i].timeData->epoch)) + 
                      j * IFOdata[i].timeData->deltaT;
                    REAL8 d = IFOdata[i].timeData->data->data[j];

                    fprintf(out, "%.6f %g\n", t, d);
                  }
                  fclose(out);

                  snprintf(filename, nameLength, "%s-freqData.dat", IFOdata[i].name);
                  out = fopen(filename, "w");
                  for (j = 0; j < IFOdata[i].freqData->data->length; j++) {
                    REAL8 f = IFOdata[i].freqData->deltaF * j;
                    REAL8 dre = IFOdata[i].freqData->data->data[j].re;
                    REAL8 dim = IFOdata[i].freqData->data->data[j].im;

                    fprintf(out, "%g %g %g\n", f, dre, dim);
                  }
                  fclose(out);
                  
                }

	}
	
	for (i=0;i<Nifo-1;i++) IFOdata[i].next=&(IFOdata[i+1]);
	
	for(i=0;i<Nifo;i++) {
		if(channels) if(channels[i]) free(channels[i]);
		if(caches) if(caches[i]) free(caches[i]);
		if(IFOnames) if(IFOnames[i]) free(IFOnames[i]);
		if(fLows) if(fLows[i]) free(fLows[i]);
		if(fHighs) if(fHighs[i]) free(fHighs[i]);
	}
	if(channels) free(channels);
	if(caches) free(caches);
	if(IFOnames) free(IFOnames);
	if(fLows) free(fLows);
	if(fHighs) free(fHighs);
	
	return headIFO;
}

void makeWhiteData(LALIFOData *IFOdata) {
  REAL8 deltaF = IFOdata->freqData->deltaF;
  REAL8 deltaT = IFOdata->timeData->deltaT;

  IFOdata->whiteFreqData = 
    XLALCreateCOMPLEX16FrequencySeries("whitened frequency data", 
                                       &(IFOdata->freqData->epoch),
                                       0.0,
                                       deltaF,
                                       &lalDimensionlessUnit,
                                       IFOdata->freqData->data->length);
  IFOdata->whiteTimeData = 
    XLALCreateREAL8TimeSeries("whitened time data",
                              &(IFOdata->timeData->epoch),
                              0.0,
                              deltaT,
                              &lalDimensionlessUnit,
                              IFOdata->timeData->data->length);


  REAL8 iLow = IFOdata->fLow / deltaF;
  REAL8 iHighDefaultCut = 0.95 * IFOdata->freqData->data->length;
  REAL8 iHighFromFHigh = IFOdata->fHigh / deltaF;
  REAL8 iHigh = (iHighDefaultCut < iHighFromFHigh ? iHighDefaultCut : iHighFromFHigh);
  REAL8 windowSquareSum = 0.0;

  UINT4 i;

  for (i = 0; i < IFOdata->freqData->data->length; i++) {
    IFOdata->whiteFreqData->data->data[i].re = IFOdata->freqData->data->data[i].re / IFOdata->oneSidedNoisePowerSpectrum->data->data[i];
    IFOdata->whiteFreqData->data->data[i].im = IFOdata->freqData->data->data[i].im / IFOdata->oneSidedNoisePowerSpectrum->data->data[i];

    if (i == 0) {
      /* Cut off the average trend in the data. */
      IFOdata->whiteFreqData->data->data[i].re = 0.0;
      IFOdata->whiteFreqData->data->data[i].im = 0.0;
    }
    if (i <= iLow) {
      /* Need to taper to implement the fLow cutoff.  Tukey window
         that starts at zero, and reaches 100% at fLow. */
      REAL8 weight = 0.5*(1.0 + cos(M_PI*(i-iLow)/iLow)); /* Starts at -Pi, runs to zero at iLow. */

      IFOdata->whiteFreqData->data->data[i].re *= weight;
      IFOdata->whiteFreqData->data->data[i].im *= weight;

      windowSquareSum += weight*weight;
    } else if (i >= iHigh) {
      /* Also taper at high freq end, Tukey window that starts at 100%
         at fHigh, then drops to zero at Nyquist.  Except that we
         always taper at least 5% of the data at high freq to avoid a
         sharp edge in freq space there. */
      REAL8 NWind = IFOdata->whiteFreqData->data->length - iHigh;
      REAL8 weight = 0.5*(1.0 + cos(M_PI*(i-iHigh)/NWind)); /* Starts at 0, runs to Pi at i = length */

      IFOdata->whiteFreqData->data->data[i].re *= weight;
      IFOdata->whiteFreqData->data->data[i].im *= weight;

      windowSquareSum += weight*weight;
    } else {
      windowSquareSum += 1.0;
    }
  }

  REAL8 norm = sqrt(IFOdata->whiteFreqData->data->length / windowSquareSum);
  for (i = 0; i < IFOdata->whiteFreqData->data->length; i++) {
    IFOdata->whiteFreqData->data->data[i].re *= norm;
    IFOdata->whiteFreqData->data->data[i].im *= norm;
  }

  XLALREAL8FreqTimeFFT(IFOdata->whiteTimeData, IFOdata->whiteFreqData, IFOdata->freqToTimeFFTPlan);
}

void injectSignal(LALIFOData *IFOdata, ProcessParamsTable *commandLine)
{
	LALStatus status;
	memset(&status,0,sizeof(status));
	SimInspiralTable *injTable=NULL;
	UINT4 Ninj=0;
	UINT4 event=0;
	UINT4 i=0,j=0;
	//CoherentGW InjectGW;
	//PPNParamStruc InjParams;
	LIGOTimeGPS injstart;
	REAL8 SNR=0,NetworkSNR=0;
	DetectorResponse det;
	memset(&injstart,0,sizeof(LIGOTimeGPS));
	//memset(&InjParams,0,sizeof(PPNParamStruc));
	COMPLEX16FrequencySeries *injF=NULL;
	LALIFOData *thisData=IFOdata->next;
	REAL8 minFlow=IFOdata->fLow;
	REAL8 MindeltaT=IFOdata->timeData->deltaT;
	REAL4TimeSeries *injectionBuffer=NULL;
	
	while(thisData){
          minFlow   = minFlow>thisData->fLow ? thisData->fLow : minFlow;
          MindeltaT = MindeltaT>thisData->timeData->deltaT ? thisData->timeData->deltaT : MindeltaT;
          thisData  = thisData->next;
	}
	//InjParams.deltaT = MindeltaT;
	//InjParams.fStartIn=(REAL4)minFlow;
	
	if(!getProcParamVal(commandLine,"--injXML")) {fprintf(stdout,"No injection file specified, not injecting\n"); return;}
	if(getProcParamVal(commandLine,"--event")) event= atoi(getProcParamVal(commandLine,"--event")->value);
	fprintf(stdout,"Injecting event %d\n",event);
	
	Ninj=SimInspiralTableFromLIGOLw(&injTable,getProcParamVal(commandLine,"--injXML")->value,0,0);
	REPORTSTATUS(&status);
	printf("Ninj %d\n", Ninj);
	if(Ninj<event) fprintf(stderr,"Error reading event %d from %s\n",event,getProcParamVal(commandLine,"--injXML")->value);
	while(i<event) {i++; injTable = injTable->next;} /* Select event */

	//memset(&InjectGW,0,sizeof(InjectGW));
	Approximant injapprox;
	LALGetApproximantFromString(&status,injTable->waveform,&injapprox);
    printf("Injecting approximant %s\n", injTable->waveform);
	REPORTSTATUS(&status);
	printf("Approximant %x\n", injapprox);
	//LALGenerateInspiral(&status,&InjectGW,injTable,&InjParams);
	//if(status.statusCode!=0) {fprintf(stderr,"Error generating injection!\n"); REPORTSTATUS(&status); }
		
	/* Begin loop over interferometers */
	while(IFOdata){
		memset(&det,0,sizeof(det));
		det.site=IFOdata->detector;
		COMPLEX8FrequencySeries *resp = XLALCreateCOMPLEX8FrequencySeries("response",&IFOdata->timeData->epoch,
																		  0.0,
																		  IFOdata->freqData->deltaF,
																		  &strainPerCount,
																		  IFOdata->freqData->data->length);
		
		for(i=0;i<resp->data->length;i++) {resp->data->data[i].re=(REAL4)1.0; resp->data->data[i].im=0.0;}
		/* Originally created for injecting into DARM-ERR, so transfer function was needed.  
		But since we are injecting into h(t), the transfer function from h(t) to h(t) is 1.*/

		/* We need a long buffer to inject into so that FindChirpInjectSignals() works properly
		 for low mass systems. Use 100 seconds here */
		REAL8 bufferLength = 100.0;
		UINT4 bufferN = (UINT4) (bufferLength/IFOdata->timeData->deltaT);
		LIGOTimeGPS bufferStart;
		memcpy(&bufferStart,&IFOdata->timeData->epoch,sizeof(LIGOTimeGPS));
		XLALGPSAdd(&bufferStart,(REAL8) IFOdata->timeData->data->length * IFOdata->timeData->deltaT);
		XLALGPSAdd(&bufferStart,-bufferLength);
		injectionBuffer=(REAL4TimeSeries *)XLALCreateREAL4TimeSeries(IFOdata->detector->frDetector.prefix,
																	 &bufferStart, 0.0, IFOdata->timeData->deltaT,
																	 &lalADCCountUnit, bufferN);
		/* This marks the sample in which the real segment starts, within the buffer */
		INT4 realStartSample=(INT4)((IFOdata->timeData->epoch.gpsSeconds - injectionBuffer->epoch.gpsSeconds)/IFOdata->timeData->deltaT);
		realStartSample+=(INT4)((IFOdata->timeData->epoch.gpsNanoSeconds - injectionBuffer->epoch.gpsNanoSeconds)*1e-9/IFOdata->timeData->deltaT);

		/*LALSimulateCoherentGW(&status,injWave,&InjectGW,&det);*/
		LALFindChirpInjectSignals(&status,injectionBuffer,injTable,resp);
		if(status.statusCode) REPORTSTATUS(&status);

		XLALDestroyCOMPLEX8FrequencySeries(resp);
		
		/* Now we cut the injection buffer down to match the time domain wave size */
		injectionBuffer=(REAL4TimeSeries *)XLALCutREAL4TimeSeries(injectionBuffer,realStartSample,IFOdata->timeData->data->length);
		
		if(status.statusCode) REPORTSTATUS(&status);
/*		for(j=0;j<injWave->data->length;j++) printf("%f\n",injWave->data->data[j]);*/
		REAL8TimeSeries *inj8Wave=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("injection8",
																			  &IFOdata->timeData->epoch,
																			  0.0,
																			  IFOdata->timeData->deltaT,
																			  &lalDimensionlessUnit,
																			  IFOdata->timeData->data->length);
		for(i=0;i<injectionBuffer->data->length;i++) inj8Wave->data->data[i]=(REAL8)injectionBuffer->data->data[i];
		XLALDestroyREAL4TimeSeries(injectionBuffer);
		injF=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("injF",
																			&IFOdata->timeData->epoch,
																			0.0,
																			IFOdata->freqData->deltaF,
																			&lalDimensionlessUnit,
																			IFOdata->freqData->data->length);
		/* Window the data */
		REAL4 WinNorm = sqrt(IFOdata->window->sumofsquares/IFOdata->window->data->length);
		for(j=0;j<inj8Wave->data->length;j++) inj8Wave->data->data[j]*=IFOdata->window->data->data[j]/WinNorm;
		XLALREAL8TimeFreqFFT(injF,inj8Wave,IFOdata->timeToFreqFFTPlan);
/*		for(j=0;j<injF->data->length;j++) printf("%lf\n",injF->data->data[j].re);*/
		if(IFOdata->oneSidedNoisePowerSpectrum){
			for(SNR=0.0,j=IFOdata->fLow/injF->deltaF;j<injF->data->length;j++){
				SNR+=pow(injF->data->data[j].re,2.0)/IFOdata->oneSidedNoisePowerSpectrum->data->data[j];
				SNR+=pow(injF->data->data[j].im,2.0)/IFOdata->oneSidedNoisePowerSpectrum->data->data[j];
			}
		}
		NetworkSNR+=SNR;
		
		/* Actually inject the waveform */
		for(j=0;j<inj8Wave->data->length;j++) IFOdata->timeData->data->data[j]+=inj8Wave->data->data[j];

FILE* file=fopen("InjSignal.dat", "w");
//FILE* file2=fopen("Noise.dat", "w");
		for(j=0;j<injF->data->length;j++){
//fprintf(file2, "%lg %lg \t %lg\n", IFOdata->freqData->deltaF*j, IFOdata->freqData->data->data[j].re, IFOdata->freqData->data->data[j].im);

			IFOdata->freqData->data->data[j].re+=injF->data->data[j].re;
			IFOdata->freqData->data->data[j].im+=injF->data->data[j].im;
fprintf(file, "%lg %lg \t %lg\n", IFOdata->freqData->deltaF*j, injF->data->data[j].re, injF->data->data[j].im);
		}
		fprintf(stdout,"Injected SNR in detector %s = %g\n",IFOdata->detector->frDetector.name,sqrt(SNR));
fclose(file);		
//fclose(file2);
		
		
		XLALDestroyREAL8TimeSeries(inj8Wave);
		XLALDestroyCOMPLEX16FrequencySeries(injF);
		IFOdata=IFOdata->next;
	}
	NetworkSNR=sqrt(NetworkSNR);
	REPORTSTATUS(&status);

	fprintf(stdout,"Network SNR of event %d = %g\n",event,NetworkSNR);
	return;
}

/* This function has a Memory Leak!  You cannot free the allocated
   header buffer (of length MAXSIZE).  Don't call it too many times!
   (It's only expected to be called once to initialize the
   differential evolution array, so this should be OK. */
char **getHeaderLine(FILE *inp) {
  const size_t MAXSIZE=1024;
  const char *delimiters = " \n\t";
  char *header = malloc(MAXSIZE*sizeof(char));
  char **colNames = NULL;  /* Will be filled in with the column names,
                              terminated by NULL. */
  size_t colNamesLen=0, colNamesMaxLen=0;
  char *colName = NULL;

  if (!fgets(header, MAXSIZE, inp)) {
    /* Some error.... */
    fprintf(stderr, "Error reading header line from file (in %s, line %d)\n",
            __FILE__, __LINE__);
    exit(1);
  } else if (strlen(header) >= MAXSIZE-1) {
    /* Probably ran out of space before reading the entire line. */
    fprintf(stderr, "Header line too long (more than %ld chars) in %s, line %d.\n",
            MAXSIZE-1, __FILE__, __LINE__);
    exit(1);
  }

  /* Sure hope we read the whole line. */
  colNamesMaxLen=2;
  colNames=(char **)malloc(2*sizeof(char *));

  if (!colNames) {
    fprintf(stderr, "Failed to allocate colNames (in %s, line %d).\n",
            __FILE__, __LINE__);
    exit(1);
  }

  colName=strtok(header, delimiters);
  strcpy(colNames[0],colNameToParamName(colName));
  //colNames[0] = colNameToParamName(colName); /* switched to strcpy() to avoid warning: assignment discards qualifiers from pointer target type */
  colNamesLen=1;
  do {
    colName=strtok(NULL, delimiters);

    strcpy(colNames[colNamesLen],colNameToParamName(colName));
    colNamesLen++;

    /* Expand if necessary. */
    if (colNamesLen >= colNamesMaxLen) {
      colNamesMaxLen *= 2;
      colNames=realloc(colNames, colNamesMaxLen*sizeof(char *));
      if (!colNames) {
        fprintf(stderr, "Failed to realloc colNames (in %s, line %d).\n",
                __FILE__, __LINE__);
        exit(1);
      }
    }

  } while (colName != NULL);

  /* Trim down to size. */
  colNames=realloc(colNames, colNamesLen*sizeof(char *));

  return colNames;
}

const char *colNameToParamName(const char *colName) {
  if (colName == NULL) {
    return NULL;
  }

  if (!strcmp(colName, "dist")) {
    return "distance";
  }

  if (!strcmp(colName, "ra")) {
    return "rightascension";
  }

  if (!strcmp(colName, "iota")) {
    return "inclination";
  }

  if (!strcmp(colName, "psi")) {
    return "polarisation";
  }

  if (!strcmp(colName, "mc")) {
    return "chirpmass";
  }

  if (!strcmp(colName, "phi_orb")) {
    return "phase";
  }

  if (!strcmp(colName, "eta")) {
    return "massratio";
  }

  if (!strcmp(colName, "dec")) {
    return "declination";
  }

  /* Note the 1 <--> 2 swap between the post-proc world and the LI world. */
  if (!strcmp(colName, "phi1")) {
    return "phi_spin2";
  }

  if (!strcmp(colName, "phi2")) {
    return "phi_spin1";
  }

  if (!strcmp(colName, "theta1")) {
    return "theta_spin2";
  }

  if (!strcmp(colName, "theta2")) {
    return "theta_spin1";
  }

  if (!strcmp(colName, "a1")) {
    return "a_spin2";
  }

  if (!strcmp(colName, "a2")) {
    return "a_spin1";
  }

  return colName;
}

int processParamLine(FILE *inp, char **headers, LALVariables *vars) {
  size_t i;

  for (i = 0; headers[i] != NULL; i++) {
    double param;
    int nread;
    
    nread = fscanf(inp, " %lg ", &param);

    if (nread != 1) {
      fprintf(stderr, "Could not read parameter value, the %ld parameter in the row (in %s, line %d)\n",
              i, __FILE__, __LINE__);
      exit(1);
    }

    addVariable(vars, headers[i], &param, REAL8_t, PARAM_FIXED);
  }

  return 0;
}
