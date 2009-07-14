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


LALIFOData *ReadData(ProcessParamsTable *commandLine)
/* Read in the data and store it in a LALIFOData structure */
{
 LALIFOData *headIFO,*curIFO,*IFOdata; 
 int Ncache=0,Nifo=0,Nchannel=0;
 int i,j;
 char strainname[]="LSC-STRAIN";
 char **channels;
 char **caches;
 char **IFOnames;
 LIGOTimeGPS GPSstart;
 REAL8 PSDdatalength=0;
 if(getProcParamVal(commandLine,"channel")){
	parseCharacterOptionString(getProcParamVal(commandLine,"channel")->value,channels,&Nchannel);
 }
 parseCharacterOptionString(getProcParamVal(commandLine,"cache")->value,caches,&Ncache);
 parseCharacterOptionString(getProcParamVal(commandLine,"IFO")->value,IFOnames,&Nifo);
 if(Nifo!=Ncache) die("ERROR: Must specify equal number of IFOs and Cache files\n");
 if(Nchannels!=0 && Nchannels!=Nifo) die("ERROR: Please specify a channel for all caches, or omit to use the defaults\n");
 IFOdata=headIFO=calloc(sizeof(LALIFOData),Nifo);
 
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
	readTseries(caches[i],channels[i],GPSstart,PSDdatalength);
 }


}