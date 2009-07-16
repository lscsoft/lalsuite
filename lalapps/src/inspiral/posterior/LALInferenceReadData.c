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

#include "LALInference.h"

REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length)
{
	LALStatus status;
	FrCache *cache = NULL;
	FrStream *stream = NULL;
	REAL8TimeSeries *out = NULL;

	cache  = XLALFrImportCache( cachefile );
	if(cache==NULL) {fprintf(stderr,"ERROR: Unable to import cache file %s\n",cachefile); exit(-1);}
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
--PSDstart GPSsecs.GPSnanosecs \n\
--PSDlength length \n\
[--srate SampleRate   [4096]] \n\
--seglen segment_length \n\
--trig_time GPSsecs.GPSnanosecs \n\
[--fLow [cutoff1,cutoff2,cutoff3,..] [40Hz]] \n\
[--fHigh [fHigh1,fHigh2,fHigh3,..] [f_Nyquist]]\n"

LALIFOData *ReadData(ProcessParamsTable *commandLine)
/* Read in the data and store it in a LALIFOData structure */
{
 LALStatus status;
 ProcessParamsTable *procparam=NULL;
 LALIFOData *headIFO=NULL,*IFOdata=NULL;
 REAL8 SampleRate=4096.0,SegmentLength=0;
 int nSegs=0;
 int seglen=0;
 REAL8TimeSeries *PSDtimeSeries=NULL,*windowedTimeData=NULL;
 REAL8 padding=1.0;
 int Ncache=0,Nifo=0,Nchannel=0,NfLow=0,NfHigh=0;
 int i;
 char strainname[]="LSC-STRAIN";

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
 seglen=(INT4)(SegmentLength*SampleRate);
 nSegs=(int)floor(PSDdatalength/SegmentLength);

 for(i=0;i<Nifo;i++) {IFOdata[i].fLow=fLows?atof(fLows[i]):40.0; IFOdata[i].fHigh=fHighs?atof(fHighs[i]):1.0/SampleRate;}

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
 /* We now have the number of detectors, let's read the PSD data */

 for(i=0;i<Nifo;i++) {
	fprintf(stderr,"Estimating PSD for %s using %i segments of %i samples (%lfs)\n",IFOnames[i],nSegs,seglen,SegmentLength);
	/* Create FFT plans */
	IFOdata[i].timeToFreqFFTPlan = XLALCreateForwardREAL8FFTPlan((UINT4) seglen, 0 );
	IFOdata[i].freqToTimeFFTPlan = XLALCreateReverseREAL8FFTPlan((UINT4) seglen,0);
	
	/* Setup windows */
	IFOdata[i].window=XLALCreateTukeyREAL8Window(seglen,(REAL8)2.0*padding*SampleRate/(REAL8)seglen);

	PSDtimeSeries=readTseries(caches[i],channels[i],GPSstart,PSDdatalength);
	if(!PSDtimeSeries) {fprintf(stderr,"Error reading PSD data for %s\n",IFOnames[i]); exit(1);}
	XLALResampleREAL8TimeSeries(PSDtimeSeries,1.0/SampleRate);
	PSDtimeSeries=(REAL8TimeSeries *)XLALShrinkREAL8TimeSeries(PSDtimeSeries,(size_t) 0, (size_t) seglen*nSegs);
	IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("inverse spectrum",&PSDtimeSeries->epoch,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
	XLALREAL8AverageSpectrumWelch(IFOdata[i].oneSidedNoisePowerSpectrum ,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan);	
	XLALDestroyREAL8TimeSeries(PSDtimeSeries);
 }

 /* Trigger time = 1 second before end of segment */
 memcpy(&segStart,&GPStrig,sizeof(LIGOTimeGPS));
 XLALGPSAdd(&segStart,-SegmentLength+1);

 /* Read and FFT the data segment */
 for(i=0;i<Nifo;i++){
	IFOdata[i].timeData=readTseries(caches[i],channels[i],segStart,SegmentLength);
	XLALResampleREAL8TimeSeries(IFOdata[i].timeData,1.0/SampleRate);	 
	if(!IFOdata[i].timeData) {fprintf(stderr,"Error reading segment data for %s\n",IFOnames[i]); exit(1);}
	IFOdata[i].freqData=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("freqData",&(IFOdata[i].timeData->epoch),0.0,1.0/SegmentLength,&lalDimensionlessUnit,seglen/2+1);
	windowedTimeData=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("temp buffer",&(IFOdata[i].timeData->epoch),0.0,1.0/SampleRate,&lalDimensionlessUnit,seglen);
	XLALDDVectorMultiply(windowedTimeData->data,IFOdata[i].timeData->data,IFOdata[i].window->data);
	XLALREAL8TimeFreqFFT(IFOdata[i].freqData,windowedTimeData,IFOdata[i].timeToFreqFFTPlan);
	XLALDestroyREAL8TimeSeries(windowedTimeData);
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

