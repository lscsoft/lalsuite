#include <stdio.h>
#include <stdlib.h>
#include <LAL/LALInspiral.h>
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

/* Variables needed from command line *****

[ --channel [channel1,channel2,channel3] ]
--IFO [IFO1,IFO2,IFO3]
--cache [cache1,cache2,cache3]
--PSDstart GPSsecs.GPSnanosecs
--PSDlength length
[--srate SampleRate   [4096]]
--seglen segment_length 
--trig_time GPSsecs.GPSnanosecs
[--fLow lower_cutoff [40Hz]]
[--fHigh high_cutoff [f_Nyquist]]


****************************************************/

LALIFOData *ReadData(ProcessParamsTable *commandLine)
/* Read in the data and store it in a LALIFOData structure */
{
 ProcessParamsTable *procparam;
 LALIFOData *headIFO,*curIFO,*IFOdata;
 REAL8 SampleRate=4096.0,Segment_Length=0,nSegs=0;
 INT4 seglen=0;
 REAL8TimeSeries *PSDtimeSeries;
 REAL8 padding=1.0;
 int Ncache=0,Nifo=0,Nchannel=0;
 int i,j;
 char strainname[]="LSC-STRAIN";
 char parambuffer[512];
 char **chartmp;
 char **channels;
 char **caches;
 char **IFOnames;
 LIGOTimeGPS GPSstart,GPStrig;
 REAL8 PSDdatalength=0;
 if(getProcParamVal(commandLine,"channel")){
	parseCharacterOptionString(getProcParamVal(commandLine,"channel")->value,channels,&Nchannel);
 }
 parseCharacterOptionString(getProcParamVal(commandLine,"cache")->value,caches,&Ncache);
 parseCharacterOptionString(getProcParamVal(commandLine,"IFO")->value,IFOnames,&Nifo);
 if(Nifo!=Ncache) die("ERROR: Must specify equal number of IFOs and Cache files\n");
 if(Nchannels!=0 && Nchannels!=Nifo) die("ERROR: Please specify a channel for all caches, or omit to use the defaults\n");
 IFOdata=headIFO=calloc(sizeof(LALIFOData),Nifo);
 
 procparam=getProcParamVal(commandLine,"PSDstart");
 LALStringToGPS(&status,&GPSstart,procparam->value,chartmp);
 PSDdatalength=atof(getProcParamVal(commandLine,"PSDlength")->value);
 SegmentLength=atof(getProcParamVal(commandLIne,"seglen")->value);
 SampleRate=atof(getProcParamVal(commandLine,"srate")->value);
 seglen=(INT4)(SegmentLength/SampleRate);

 if(Nchannels==0)
 {
	channels=calloc(sizeof(char *),Nifo);
	for(i=0;i<Nifo;i++) {
		channels[i]=malloc(VARNAME_MAX);
		if(!strcmp(IFOnames[i],"H1")) {
			IFOdata[i].detector=&lalCachedDetectors[LALDetectorIndexLHODIFF];
			sprintf((channels[i]),"H1:%s",strainname); continue;}
		if(!strcmp(IFOnames[i],"H2")) {
			IFOdata[i].detector=&lalCachedDetectors[LALDetectorIndexLHODIFF];
			sprintf((channels[i]),"H2:%s",strainname); continue;}
		if(!strcmp(IFOnames[i],"LLO")||!strcmp(IFOnames[i],"L1")) {
			IFOdata[i].detector=&lalCachedDetectors[LALDetectorIndexLLODIFF];
			sprintf((channels[i]),"L1:%s",strainname); continue;}
		if(!strcmp(IFOnames[i],"V1")||!strcmp(IFOnames[i],"VIRGO")) {
			IFOdata[i].detector=&lalCachedDetectors[LALDetectorIndexVIRGODIFF];
			sprintf((channels[i]),"V1:h_16384Hz");
		if(!strcmp(IFOnames[i],"GEO")||!strcmp(IFOnames[i],"G1")) {
			IFOdata[i].detector=&lalCachedDetectors[LALDetectorIndexGEO600DIFF];
			sprintf((channels[i]),"G1:DER_DATA_H"); continue;}
/*		if(!strcmp(IFOnames[i],"TAMA")||!strcmp(IFOnames[i],"T1")) {inputMCMC.detector[i]=&lalCachedDetectors[LALDetectorIndexTAMA300DIFF]; continue;}*/
		fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 H2 L1 V1 GEO\n",IFOnames[i]); exit(-1);
	}
 }
 /* We now have the number of detectors, let's read the PSD data */

 for(i=0;i<Nifo;i++) {
	/* Create FFT plans */
	IFOData[i].timeToFreqFFTPlan = XLALCreateForwardREAL4FFTPlan( seglen, 0 );
	IDOData[i].freqToTimeFFTPlan = XLALCreateReverseREAL4FFTPlan(seglen,0);
	
	/* Setup windows */
	IFOdata[i].window=XLALCreateTukeyREAL8Window(seglen,(REAL8)2.0*padding*SampleRate/(REAL8)seglen);

	PSDTimeSeries=readTseries(caches[i],channels[i],GPSstart,PSDdatalength);
	XLALResampleREAL8TimeSeries(PSDTimeSeries,1.0/SampleRate);
	PSDTimeSeries=(REAL8TimeSeries *)XLALShrinkREAL8TimeSeries(PSDTimeSeries,(size_t) 0, (size_t) seglen*nSegs);
	IFOdata[i].oneSidedNoisePowerSpectrum=XLALCreateREAL8FrequencySeries("inverse spectrum",&PSDTimeSeries->epoch,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
	XLALREAL8AverageSpectrumWelch(IFOdata[i].oneSidedNoisePowerSpectrum ,PSDTimeSeries, seglen, (UINT4)stride, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan);
 }

 /* Read and FFT the data segment */

}