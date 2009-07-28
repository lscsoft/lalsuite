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
#include <lal/LALInspiral.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/Units.h>
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

#include "LALInference.h"

REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length)
{
	LALStatus status;
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
[--fHigh [fHigh1,fHigh2,fHigh3,..] [f_Nyquist]]\n"

LALIFOData *readData(ProcessParamsTable *commandLine)
/* Read in the data and store it in a LALIFOData structure */
{
	LALStatus status;
	ProcessParamsTable *procparam=NULL;
	LALIFOData *headIFO=NULL,*IFOdata=NULL;
	REAL8 SampleRate=4096.0,SegmentLength=0;
	int nSegs=0;
	size_t seglen=0;
	REAL8TimeSeries *PSDtimeSeries=NULL,*windowedTimeData=NULL;
	REAL8 padding=1.0;
	int Ncache=0,Nifo=0,Nchannel=0,NfLow=0,NfHigh=0;
	int i,j;
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
	memset(&status,0,sizeof(LALStatus));

	if(!getProcParamVal(commandLine,"--cache")||!getProcParamVal(commandLine,"--IFO")||
	   !getProcParamVal(commandLine,"--PSDstart")||!getProcParamVal(commandLine,"--trigtime")||
	   !getProcParamVal(commandLine,"--PSDlength")||!getProcParamVal(commandLine,"--seglen"))
	{fprintf(stderr,USAGE); return(NULL);}
	
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
	if(Nifo!=Ncache) {fprintf(stderr,"ERROR: Must specify equal number of IFOs and Cache files\n"); exit(1);}
	if(Nchannel!=0 && Nchannel!=Nifo) {fprintf(stderr,"ERROR: Please specify a channel for all caches, or omit to use the defaults\n"); exit(1);}
	
	IFOdata=headIFO=calloc(sizeof(LALIFOData),Nifo);
	
	procparam=getProcParamVal(commandLine,"--PSDstart");
	LALStringToGPS(&status,&GPSstart,procparam->value,&chartmp);
	procparam=getProcParamVal(commandLine,"--trigtime");
	LALStringToGPS(&status,&GPStrig,procparam->value,&chartmp);
	PSDdatalength=atof(getProcParamVal(commandLine,"--PSDlength")->value);
	SegmentLength=atof(getProcParamVal(commandLine,"--seglen")->value);
	if(getProcParamVal(commandLine,"--srate")) SampleRate=atof(getProcParamVal(commandLine,"--srate")->value);
	seglen=(size_t)(SegmentLength*SampleRate);
	nSegs=(int)floor(PSDdatalength/SegmentLength);
	
	for(i=0;i<Nifo;i++) {IFOdata[i].fLow=fLows?atof(fLows[i]):40.0; IFOdata[i].fHigh=fHighs?atof(fHighs[i]):SampleRate;}
	
	if(Nchannel==0)
	{
		channels=calloc(Nifo,sizeof(char *));
		for(i=0;i<Nifo;i++) {
			channels[i]=malloc(VARNAME_MAX);
			IFOdata[i].detector=calloc(1,sizeof(LALDetector));
			if(!strcmp(IFOnames[i],"H1")) {			
				memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
			sprintf((channels[i]),"H1:%s",strainname); continue;}
			if(!strcmp(IFOnames[i],"H2")) {
				memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
			sprintf((channels[i]),"H2:%s",strainname); continue;}
			if(!strcmp(IFOnames[i],"LLO")||!strcmp(IFOnames[i],"L1")) {
				memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexLLODIFF],sizeof(LALDetector));
			sprintf((channels[i]),"L1:%s",strainname); continue;}
			if(!strcmp(IFOnames[i],"V1")||!strcmp(IFOnames[i],"VIRGO")) {
				memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexVIRGODIFF],sizeof(LALDetector));
			sprintf((channels[i]),"V1:h_16384Hz");}
			if(!strcmp(IFOnames[i],"GEO")||!strcmp(IFOnames[i],"G1")) {
				memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexGEO600DIFF],sizeof(LALDetector));
			sprintf((channels[i]),"G1:DER_DATA_H"); continue;}
			/*		if(!strcmp(IFOnames[i],"TAMA")||!strcmp(IFOnames[i],"T1")) {inputMCMC.detector[i]=&lalCachedDetectors[LALDetectorIndexTAMA300DIFF]; continue;}*/
			fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 H2 L1 V1 GEO\n",IFOnames[i]); exit(-1);
		}
	}
	/* Set up FFT structures and window */
	for (i=0;i<Nifo;i++){
		/* Create FFT plans */
		IFOdata[i].timeToFreqFFTPlan = XLALCreateForwardREAL8FFTPlan((UINT4) seglen, 0 );
		IFOdata[i].freqToTimeFFTPlan = XLALCreateReverseREAL8FFTPlan((UINT4) seglen,0);
		
		/* Setup windows */
		IFOdata[i].window=XLALCreateTukeyREAL8Window(seglen,(REAL8)2.0*padding*SampleRate/(REAL8)seglen);
	}
	
	
	/* Trigger time = 1 second before end of segment */
	memcpy(&segStart,&GPStrig,sizeof(LIGOTimeGPS));
	XLALGPSAdd(&segStart,-SegmentLength+1);
	
	
	/* Read the PSD data */
	for(i=0;i<Nifo;i++) {
		/* Check if fake data is requested */
		if(!(strcmp(caches[i],"LALLIGO") && strcmp(caches[i],"LALVirgo") && strcmp(caches[i],"LALGEO") && strcmp(caches[i],"LALEGO")
			 && strcmp(caches[i],"LALAdLIGO")))
		{
			FakeFlag=1;
			datarandparam=XLALCreateRandomParams(0);
			/* Selection of the noise curve */
			if(!strcmp(caches[i],"LALLIGO")) {PSD = &LALLIGOIPsd; scalefactor=9E-46;}
			if(!strcmp(caches[i],"LALVirgo")) {PSD = &LALVIRGOPsd; scalefactor=1.0;}
			if(!strcmp(caches[i],"LALGEO")) {PSD = &LALGEOPsd; scalefactor=1E-46;}
			if(!strcmp(caches[i],"LALEGO")) {PSD = &LALEGOPsd; scalefactor=1.0;}
			if(!strcmp(caches[i],"LALAdLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 10E-49;}
			if(!strcmp(caches[i],"LAL2kLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 36E-46;}
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
			for(j=0;j<IFOdata[i].freqData->data->length;j++){
				IFOdata[i].freqData->data->data[j].re=XLALNormalDeviate(datarandparam)*(0.5*sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]*IFOdata[i].freqData->deltaF));
				IFOdata[i].freqData->data->data[j].im=XLALNormalDeviate(datarandparam)*(0.5*sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]*IFOdata[i].freqData->deltaF));
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
			XLALREAL8AverageSpectrumWelch(IFOdata[i].oneSidedNoisePowerSpectrum ,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan);	
			XLALDestroyREAL8TimeSeries(PSDtimeSeries);
			
			/* Read the data segment */
			IFOdata[i].timeData=readTseries(caches[i],channels[i],segStart,SegmentLength);
			if(!IFOdata[i].timeData) {fprintf(stderr,"Error reading segment data for %s at %i\n",IFOnames[i],segStart.gpsSeconds); exit(1);}
			XLALResampleREAL8TimeSeries(IFOdata[i].timeData,1.0/SampleRate);	 
			if(!IFOdata[i].timeData) {fprintf(stderr,"Error reading segment data for %s\n",IFOnames[i]); exit(1);}
			IFOdata[i].freqData=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("freqData",&(IFOdata[i].timeData->epoch),0.0,1.0/SegmentLength,&lalDimensionlessUnit,seglen/2+1);
			windowedTimeData=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("temp buffer",&(IFOdata[i].timeData->epoch),0.0,1.0/SampleRate,&lalDimensionlessUnit,seglen);
			XLALDDVectorMultiply(windowedTimeData->data,IFOdata[i].timeData->data,IFOdata[i].window->data);
			XLALREAL8TimeFreqFFT(IFOdata[i].freqData,windowedTimeData,IFOdata[i].timeToFreqFFTPlan);
			XLALDestroyREAL8TimeSeries(windowedTimeData);
			
			for(j=0;j<IFOdata[i].freqData->data->length;j++){
				IFOdata[i].freqData->data->data[j].re/=sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
				IFOdata[i].freqData->data->data[j].im/=sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
			}
			
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

void injectSignal(LALIFOData *IFOdata, ProcessParamsTable *commandLine)
{
	LALStatus status;
	SimInspiralTable *injTable=NULL;
	INT4 Ninj=0;
	INT4 event=0;
	int i=0,j=0;
	CoherentGW InjectGW;
	PPNParamStruc InjParams;
	LIGOTimeGPS injstart;
	REAL8 SNR=0,NetworkSNR=0;
	DetectorResponse det;
	memset(&injstart,0,sizeof(LIGOTimeGPS));
	memset(&InjParams,0,sizeof(PPNParamStruc));
	COMPLEX16FrequencySeries *injF=NULL;
	LALIFOData *thisData=IFOdata->next;
	REAL8 minFlow=IFOdata->fLow;
	REAL8 MindeltaT=IFOdata->timeData->deltaT;
	
	while(thisData){
		minFlow=minFlow>thisData->fLow?thisData->fLow:minFlow;
		MindeltaT=MindeltaT>thisData->timeData->deltaT?thisData->timeData->deltaT:MindeltaT;
	}
	InjParams.deltaT = MindeltaT;
	InjParams.fStartIn=(REAL4)minFlow;
	
	if(!getProcParamVal(commandLine,"--injXML")) {fprintf(stdout,"No injection file specified, not injecting\n"); return;}
	if(getProcParamVal(commandLine,"--event")) event=getProcParamVal(commandLine,"--event")->value;
	fprintf(stdout,"Injecting event %d\n",event);
	
	Ninj=SimInspiralTableFromLIGOLw(&injTable,getProcParamVal(commandLine,"--injXML")->value,0,0);
	if(Ninj<event) fprintf(stderr,"Error reading event %d from %s\n",event,getProcParamVal(commandLine,"--injXML")->value);
	while(i<event) {i++; injTable = injTable->next;} /* Select event */

	memset(&InjectGW,0,sizeof(InjectGW));
	Approximant injapprox;
	LALGetApproximantFromString(&status,injTable->waveform,&injapprox);
	LALGenerateInspiral(&status,&InjectGW,injTable,&InjParams);
	if(status.statusCode!=0) {fprintf(stderr,"Error generating injection!!!\n"); REPORTSTATUS(&status); }
	
	/* Begin loop over interferometers */
	while(IFOdata){
		memset(&det,0,sizeof(det));
		det.site=IFOdata->detector;
		REAL4TimeSeries *injWave=(REAL4TimeSeries *)XLALCreateREAL4TimeSeries("injection",
																			  &IFOdata->timeData->epoch,
																			  0.0,
																			  IFOdata->timeData->deltaT,
																			  &lalDimensionlessUnit,
																			  IFOdata->timeData->data->length);
		LALSimulateCoherentGW(&status,injWave,&InjectGW,&det);
		if(status.statusCode) REPORTSTATUS(&status);
		REAL8TimeSeries *inj8Wave=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("injection8",
																			  &IFOdata->timeData->epoch,
																			  0.0,
																			  IFOdata->timeData->deltaT,
																			  &lalDimensionlessUnit,
																			  IFOdata->timeData->data->length);
		for(i=0;i<injWave->data->length;i++) inj8Wave->data->data[i]=(REAL8)injWave->data->data[i];
		XLALDestroyREAL4TimeSeries(injWave);
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
		if(IFOdata->oneSidedNoisePowerSpectrum){
			for(SNR=0.0,j=IFOdata->fLow/injF->deltaF;j<injF->data->length;j++){
				SNR+=pow(injF->data->data[j].re,2.0)*IFOdata->oneSidedNoisePowerSpectrum->data->data[j];
				SNR+=pow(injF->data->data[j].im,2.0)*IFOdata->oneSidedNoisePowerSpectrum->data->data[j];
			}
		}
		NetworkSNR+=SNR;
		
		/* Actually inject the waveform */
		for(j=0;j<inj8Wave->data->length;j++) IFOdata->timeData->data->data[j]+=inj8Wave->data->data[j];
		for(j=0;j<injF->data->length;j++){
			IFOdata->freqData->data->data[j].re+=injF->data->data[j].re;
			IFOdata->freqData->data->data[j].im+=injF->data->data[j].im;
		}
		fprintf(stdout,"Injected SNR in detector %s = %g\n",IFOdata->detector->frDetector.name,SNR);
		XLALDestroyREAL8TimeSeries(inj8Wave);
		XLALDestroyCOMPLEX16FrequencySeries(injF);
		IFOdata=IFOdata->next;
	}
	fprintf(stdout,"Network SNR of event %d = %g\n",event,NetworkSNR);
	return;
}
