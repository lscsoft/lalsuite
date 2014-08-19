/*
 *  LALInferenceReadData.c:  Bayesian Followup functions
 *
 *  Copyright (C) 2009,2012 Ilya Mandel, Vivien Raymond, Christian
 *  Roever, Marc van der Sluys, John Veitch, Salvatore Vitale, and
 *  Will M. Farr
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
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
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

#include <lal/SeqFactories.h>
#include <lal/DetectorSite.h>
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/Inject.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/LALInspiralBank.h>
#include <lal/FindChirp.h>
#include <lal/LALInspiralBank.h>
#include <lal/GenerateInspiral.h>
#include <lal/NRWaveInject.h>
#include <lal/GenerateInspRing.h>
#include <lal/LALErrno.h>
#include <math.h>
#include <lal/LALInspiral.h>
#include <lal/LALSimulation.h>

#include <lal/LALInference.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALSimNoise.h>
#include <LALInferenceRemoveLines.h>

struct fvec {
	REAL8 f;
	REAL8 x;
};

#define LALINFERENCE_DEFAULT_FLOW "40.0"

char *SNRpath = NULL;

struct fvec *interpFromFile(char *filename);

struct fvec *interpFromFile(char *filename){
	UINT4 fileLength=0;
	UINT4 i=0;
	UINT4 minLength=100; /* size of initial file buffer, and also size of increment */
	FILE *interpfile=NULL;
	struct fvec *interp=NULL;
	interp=XLALCalloc(minLength,sizeof(struct fvec)); /* Initialise array */
	if(!interp) {printf("Unable to allocate memory buffer for reading interpolation file\n");}
	fileLength=minLength;
	REAL8 f=0.0,x=0.0;
	interpfile = fopen(filename,"r");
	if (interpfile==NULL){
		printf("Unable to open file %s\n",filename);
		exit(1);
	}
	while(2==fscanf(interpfile," %lf %lf ", &f, &x )){
		interp[i].f=f; interp[i].x=x*x;
		i++;
		if(i>fileLength-1){ /* Grow the array */
			interp=XLALRealloc(interp,(fileLength+minLength)*sizeof(struct fvec));
			fileLength+=minLength;
		}
	}
	interp[i].f=0; interp[i].x=0;
	fileLength=i+1;
	interp=XLALRealloc(interp,fileLength*sizeof(struct fvec)); /* Resize array */
	fclose(interpfile);
	printf("Read %i records from %s\n",fileLength-1,filename);
	return interp;
}

REAL8 interpolate(struct fvec *fvec, REAL8 f);
REAL8 interpolate(struct fvec *fvec, REAL8 f){
	int i=0;
	REAL8 a=0.0; /* fractional distance between bins */
	REAL8 delta=0.0;
	if(f<fvec[0].f) return(0.0);
	while(fvec[i].f<f && (fvec[i].x!=0.0 )){i++;}; //&& fvec[i].f!=0.0)){i++;};
  //printf("%d\t%lg\t%lg\t%lg\n",i,fvec[i].f,f,fvec[i].x);
	if (fvec[i].f==0.0 && fvec[i].x==0.0) /* Frequency above moximum */
	{
		return (fvec[i-1].x);
	}
  //if(i==0){return (fvec[0].x);}
	a=(fvec[i].f-f)/(fvec[i].f-fvec[i-1].f);
	delta=fvec[i].x-fvec[i-1].x;
	return (fvec[i-1].x + delta*a);
}
void InjectFD(LALInferenceIFOData *IFOdata, SimInspiralTable *inj_table, ProcessParamsTable *commandLine);
void enforce_m1_larger_m2(SimInspiralTable* injEvent);

typedef void (NoiseFunc)(LALStatus *statusPtr,REAL8 *psd,REAL8 f);
void MetaNoiseFunc(LALStatus *status, REAL8 *psd, REAL8 f, struct fvec *interp, NoiseFunc *noisefunc);

void MetaNoiseFunc(LALStatus *status, REAL8 *psd, REAL8 f, struct fvec *interp, NoiseFunc *noisefunc){
	if(interp==NULL&&noisefunc==NULL){
		printf("ERROR: Trying to calculate PSD with NULL inputs\n");
		exit(1);
	}
	if(interp!=NULL && noisefunc!=NULL){
		printf("ERROR: You have specified both an interpolation vector and a function to calculate the PSD\n");
		exit(1);
	}
	if(noisefunc!=NULL){
		noisefunc(status,psd,f);
		return;
	}
	if(interp!=NULL){ /* Use linear interpolation of the interp vector */
		*psd=interpolate(interp,f);
		return;
	}
}

void
LALInferenceLALFindChirpInjectSignals (
                                       LALStatus                  *status,
                                       REAL4TimeSeries            *chan,
                                       SimInspiralTable           *events,
                                       COMPLEX8FrequencySeries    *resp,
                                       LALDetector                *detector
                                       );
static int FindTimeSeriesStartAndEnd (
                                      REAL4Vector *signalvec,
                                      UINT4 *start,
                                      UINT4 *end
                                      );

static const LALUnit strainPerCount={0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

static REAL8TimeSeries *readTseries(LALCache *cache, CHAR *channel, LIGOTimeGPS start, REAL8 length);
static void makeWhiteData(LALInferenceIFOData *IFOdata);
static void PrintSNRsToFile(LALInferenceIFOData *IFOdata , SimInspiralTable *inj_table);


static LALCache *GlobFramesPWD( char *ifo);
static LALCache *GlobFramesPWD(char *ifo)
{
        LALCache *frGlobCache = NULL;

        /* create a frame cache by globbing all *.gwf files in the pwd */
        /* FIXME: This should really open all the files and see if the desired channel is in there */
        char globPattern[8];
        sprintf(globPattern,"%c-*.gwf",ifo[0]);
        frGlobCache = XLALCacheGlob(NULL,globPattern);

        /* check we globbed at least one frame file */
        if ( ! frGlobCache->length )
        {
            fprintf( stderr, "error: no frame file files found\n");
            exit( 1 );
        }
    CHAR ifoRegExPattern[6];
    LALCache *frInCache=NULL;
    /* sieve out the requested data type */
        snprintf( ifoRegExPattern,
                sizeof(ifoRegExPattern) / sizeof(*ifoRegExPattern), ".*%c.*",
                ifo[0] );
    {
        fprintf(stderr,"GlobFramesPWD : Found unseived src files:\n");
        for(UINT4 i=0;i<frGlobCache->length;i++)
            fprintf(stderr,"(%s,%s,%s)\n",frGlobCache->list[i].src,frGlobCache->list[i].dsc,frGlobCache->list[i].url);
    }
    frInCache = XLALCacheDuplicate(frGlobCache);
    XLALCacheSieve(frInCache, 0, 0, ifoRegExPattern, NULL, NULL);
    {
        fprintf(stderr,"GlobFramesPWD : Sieved frames with pattern %s. Found src files:\n",ifoRegExPattern);
        for(UINT4 i=0;i<frInCache->length;i++)
            fprintf(stderr,"(%s,%s,%s)\n",frInCache->list[i].src,frInCache->list[i].dsc,frInCache->list[i].url);
    }
    //XLALDestroyCache(frGlobCache);
    
    return(frGlobCache);
}

static REAL8TimeSeries *readTseries(LALCache *cache, CHAR *channel, LIGOTimeGPS start, REAL8 length)
{
	LALStatus status;
	memset(&status,0,sizeof(status));
	LALFrStream *stream = NULL;
	REAL8TimeSeries *out = NULL;
	if(cache==NULL) fprintf(stderr,"readTseries ERROR: Received NULL pointer for channel %s\n",channel);
	stream = XLALFrStreamCacheOpen( cache );
	if(stream==NULL) {fprintf(stderr,"readTseries ERROR: Unable to open stream from frame cache file\n"); exit(-1);}
	out = XLALFrStreamInputREAL8TimeSeries( stream, channel, &start, length , 0 );
	if(out==NULL) fprintf(stderr,"readTseries ERROR: unable to read channel %s at time %i\nCheck the specified data duration is not too long\n",channel,start.gpsSeconds);
	LALFrClose(&status,&stream);
	return out;
}

/**
 * Parse the command line looking for options of the kind --ifo H1 --H1-channel H1:LDAS_STRAIN --H1-cache H1.cache --H1-flow 40.0 --H1-fhigh 4096.0 --H1-timeslide 100.0 ...
 * It is necessary to use this method instead of the old method for the pipeline to work in DAX mode. Warning: do not mix options between
 * the old and new style.
 */
static INT4 getDataOptionsByDetectors(ProcessParamsTable *commandLine, char ***ifos, char ***caches, char ***channels, char ***flows , char ***fhighs, char ***timeslides, UINT4 *N)
{
    /* Check that the input has no lists with [ifo,ifo] */
    ProcessParamsTable *this=commandLine;
    UINT4 i=0;
    *caches=*ifos=*channels=*flows=*fhighs=*timeslides=NULL;
    *N=0;
    char tmp[128];
    if(!this) {fprintf(stderr,"No command line arguments given!\n"); exit(1);}
    while(this)
    {
        if(!strcmp(this->param,"--ifo") || !strcmp(this->param,"--IFO"))
        for(i=0;this->value[i]!='\0';i++)
            if(this->value[i]=='[' || this->value[i]==']')
            {
                fprintf(stderr,"Found old-style input arguments for %s\n",this->param);
                return(0);
            }
        this=this->next;
    }
    /* Construct a list of IFOs */
    for(this=commandLine;this;this=this->next)
    {
        if(!strcmp(this->param,"--ifo")||!strcmp(this->param,"--IFO"))
        {
            (*N)++;
            *ifos=XLALRealloc(*ifos,*N*sizeof(char *));
            (*ifos)[*N-1]=XLALStringDuplicate(this->value);
        }
    }
    *caches=XLALCalloc(*N,sizeof(char *));
    *channels=XLALCalloc(*N,sizeof(char *));
    *flows=XLALCalloc(*N,sizeof(REAL8));
    *fhighs=XLALCalloc(*N,sizeof(REAL8));
    *timeslides=XLALCalloc(*N,sizeof(REAL8));

    int globFrames=!!LALInferenceGetProcParamVal(commandLine,"--glob-frame-data");

    /* For each IFO, fetch the other options if available */
    for(i=0;i<*N;i++)
    {
        /* Cache */
        if(!globFrames){
            sprintf(tmp,"--%s-cache",(*ifos)[i]);
            this=LALInferenceGetProcParamVal(commandLine,tmp);
            if(!this){fprintf(stderr,"ERROR: Must specify a cache file for %s with --%s-cache\n",(*ifos)[i],(*ifos)[i]); exit(1);}
            (*caches)[i]=XLALStringDuplicate(this->value);
        }
        
        /* Channel */
        sprintf(tmp,"--%s-channel",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*channels)[i]=XLALStringDuplicate(this?this->value:"Unknown channel");

        /* flow */
        sprintf(tmp,"--%s-flow",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*flows)[i]=XLALStringDuplicate(this?this->value:LALINFERENCE_DEFAULT_FLOW);
        
        /* fhigh */
        sprintf(tmp,"--%s-fhigh",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*fhighs)[i]=this?XLALStringDuplicate(this->value):NULL;

        /* timeslides */
        sprintf(tmp,"--%s-timeslide",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*timeslides)[i]=XLALStringDuplicate(this?this->value:"0.0");

    }
    return(1);
}

void LALInferencePrintDataWithInjection(LALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine);
void LALInferencePrintDataWithInjection(LALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine){
  
  UINT4 Nifo=0,i,j;
  LALInferenceIFOData *thisData=IFOdata;
  ProcessParamsTable *ppt=NULL;
  ProcessParamsTable *pptdatadump=NULL;
  while(thisData){
    thisData=thisData->next;
    Nifo++;
  }
  
  
  if (LALInferenceGetProcParamVal(commandLine, "--data-dump")) {
    pptdatadump=LALInferenceGetProcParamVal(commandLine,"--data-dump");
    const UINT4 nameLength=FILENAME_MAX;
    char filename[nameLength];
    FILE *out;
    
    for (i=0;i<Nifo;i++) {
      
      ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");
      if(ppt) {
        snprintf(filename, nameLength, "%s-%s-timeDataWithInjection.dat", ppt->value, IFOdata[i].name);
      }
      else if(strcmp(pptdatadump->value,"")) {
        snprintf(filename, nameLength, "%s/%s-timeDataWithInjection.dat", pptdatadump->value, IFOdata[i].name);
      }
      else
        snprintf(filename, nameLength, "%s-timeDataWithInjection.dat", IFOdata[i].name);
      out = fopen(filename, "w");
      for (j = 0; j < IFOdata[i].timeData->data->length; j++) {
        REAL8 t = XLALGPSGetREAL8(&(IFOdata[i].timeData->epoch)) +
        j * IFOdata[i].timeData->deltaT;
        REAL8 d = IFOdata[i].timeData->data->data[j];
        
        fprintf(out, "%.6f %g\n", t, d);
      }
      fclose(out);
      
      ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");
      if(ppt) {
        snprintf(filename, nameLength, "%s-%s-freqDataWithInjection.dat", ppt->value, IFOdata[i].name);
      }
      else if(strcmp(pptdatadump->value,"")) {
        snprintf(filename, nameLength, "%s/%s-freqDataWithInjection.dat", pptdatadump->value, IFOdata[i].name);
      }
      else
        snprintf(filename, nameLength, "%s-freqDataWithInjection.dat", IFOdata[i].name);
      out = fopen(filename, "w");
      for (j = 0; j < IFOdata[i].freqData->data->length; j++) {
        REAL8 f = IFOdata[i].freqData->deltaF * j;
        REAL8 dre = creal(IFOdata[i].freqData->data->data[j]);
        REAL8 dim = cimag(IFOdata[i].freqData->data->data[j]);
        
        fprintf(out, "%10.10g %10.10g %10.10g\n", f, dre, dim);
      }
      fclose(out);
    }
    
  }
  
}

#define USAGE "\
 --ifo IFO1 [--ifo IFO2 ...]    IFOs can be H1,L1,V1\n\
 --IFO1-cache cache1 [--IFO2-cache2 cache2 ...]    cache files (LALLIGO, LALAdLIGO, LALVirgo to simulate these detectors using lal; LALSimLIGO, LALSimAdLIGO, LALSimVirgo, LALSimAdVirgo to use lalsimuation)\n\
 --psdstart GPStime             GPS start time of PSD estimation data\n\
 --psdlength length             length of PSD estimation data in seconds\n\
 --seglen length                length of segments for PSD estimation and analysis in seconds\n\
(--trigtime GPStime)            GPS time of the trigger to analyse (optional when using --margtime or --margtimephi)\n\
(--segment-start)               GPS time of the start of the segment (optional when --trigtime given, default is seglen-2 s before --trigtime)\n\
(--srate rate)                  Downsample data to rate in Hz (4096.0,)\n\
(--injectionsrate rate)         Downsample injection signal to rate in Hz (--srate)\n\
(--IFO1-flow freq1 [--IFO2-flow freq2 ...])      Specify lower frequency cutoff for overlap integral (40.0)\n\
(--IFO1-fhigh freq1 [--IFO2-fhigh freq2 ...])     Specify higher frequency cutoff for overlap integral (Nyquist freq 0.5*srate)\n\
(--IFO1-channel chan1 [--IFO2-channel chan2 ...])   Specify channel names when reading cache files\n\
(--dataseed number)             Specify random seed to use when generating data\n\
(--lalinspiralinjection)      Enables injections via the LALInspiral package\n\
(--inj-fref)                    Reference frequency for parameters in injection XML\n\
(--inj-lambda1)                 value of lambda1 to be injected, LALSimulation only (0)\n\
(--inj-lambda2)                 value of lambda1 to be injected, LALSimulation only (0)\n\
(--inj-spinOrder PNorder)           Specify twice the PN order (e.g. 5 <==> 2.5PN) of spin effects to use, only for LALSimulation (default: -1 <==> Use all spin effects).\n\
(--inj-tidalOrder PNorder)          Specify twice the PN order (e.g. 10 <==> 5PN) of tidal effects to use, only for LALSimulation (default: -1 <==> Use all tidal effects).\n\
(--snrpath) 			Set a folder where to write a file with the SNRs being injected\n\
(--0noise)                      Sets the noise realisation to be identically zero (for the fake caches above only)\n"


LALInferenceIFOData *LALInferenceReadData(ProcessParamsTable *commandLine)
/* Read in the data and store it in a LALInferenceIFOData structure */
{
    LALStatus status;
    INT4 dataseed=0;
    memset(&status,0,sizeof(status));
    ProcessParamsTable *procparam=NULL,*ppt=NULL;
    ProcessParamsTable *pptdatadump=NULL;
    LALInferenceIFOData *headIFO=NULL,*IFOdata=NULL;
    REAL8 SampleRate=4096.0,SegmentLength=0;
    if(LALInferenceGetProcParamVal(commandLine,"--srate")) SampleRate=atof(LALInferenceGetProcParamVal(commandLine,"--srate")->value);
    REAL8 defaultFLow = atof(LALINFERENCE_DEFAULT_FLOW);
    int nSegs=0;
    size_t seglen=0;
    REAL8TimeSeries *PSDtimeSeries=NULL;
    REAL8 padding=0.4;//Default was 1.0 second. However for The Event the Common Inputs specify a Tukey parameter of 0.1, so 0.4 second of padding for 8 seconds of data.
    UINT4 Ncache=0,Npsd=0,Nifo=0,Nchannel=0,NfLow=0,NfHigh=0;
    UINT4 i,j;
    //int FakeFlag=0; - set but not used
    char strainname[]="LSC-STRAIN";
    UINT4 q=0;	
    //typedef void (NoiseFunc)(LALStatus *statusPtr,REAL8 *psd,REAL8 f);
    NoiseFunc *PSD=NULL;
    REAL8 scalefactor=1;
    SimInspiralTable *injTable=NULL;
    RandomParams *datarandparam;
    UINT4 event=0;
    int globFrames=0; // 0 = no, 1 = will search for frames in PWD
    char *chartmp=NULL;
    char **channels=NULL;
    char **caches=NULL;
    char **psds=NULL;
    char **IFOnames=NULL;
    char **fLows=NULL,**fHighs=NULL;
    char **timeslides=NULL;
    UINT4 Ntimeslides=0;
    LIGOTimeGPS GPSstart,GPStrig,segStart;
    REAL8 PSDdatalength=0;
    REAL8 AIGOang=0.0; //orientation angle for the proposed Australian detector.
    procparam=LALInferenceGetProcParamVal(commandLine,"--aigoang");
    if(!procparam) procparam=LALInferenceGetProcParamVal(commandLine,"--AIGOang");
    if(procparam)
        AIGOang=atof(procparam->value)*LAL_PI/180.0;

    struct fvec *interp;
    int interpFlag=0;

    if(LALInferenceGetProcParamVal(commandLine,"--glob-frame-data")) globFrames=1;

    /* Check if the new style command line arguments are used */
    INT4 dataOpts=getDataOptionsByDetectors(commandLine, &IFOnames, &caches, &channels, &fLows , &fHighs, &timeslides, &Nifo);
    /* Check for options if not given in the new style */
    if(!dataOpts){
        if(!(globFrames||LALInferenceGetProcParamVal(commandLine,"--cache"))||!(LALInferenceGetProcParamVal(commandLine,"--IFO")||LALInferenceGetProcParamVal(commandLine,"--ifo")))
            {fprintf(stderr,USAGE); return(NULL);}
        if(LALInferenceGetProcParamVal(commandLine,"--channel")){
            LALInferenceParseCharacterOptionString(LALInferenceGetProcParamVal(commandLine,"--channel")->value,&channels,&Nchannel);
        }
        LALInferenceParseCharacterOptionString(LALInferenceGetProcParamVal(commandLine,"--cache")->value,&caches,&Ncache);
        ppt=LALInferenceGetProcParamVal(commandLine,"--ifo");
        if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--IFO");
        LALInferenceParseCharacterOptionString(ppt->value,&IFOnames,&Nifo);

        ppt=LALInferenceGetProcParamVal(commandLine,"--flow");
        if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--fLow");
        if(ppt){
            LALInferenceParseCharacterOptionString(ppt->value,&fLows,&NfLow);
        }
        ppt=LALInferenceGetProcParamVal(commandLine,"--fhigh");
        if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--fHigh");
        if(ppt){
            LALInferenceParseCharacterOptionString(ppt->value,&fHighs,&NfHigh);
        }
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--timeslide"))) LALInferenceParseCharacterOptionString(ppt->value,&timeslides,&Ntimeslides);
        if(Nifo!=Ncache) {fprintf(stderr,"ERROR: Must specify equal number of IFOs and Cache files\n"); exit(1);}
        if(Nchannel!=0 && Nchannel!=Nifo) {fprintf(stderr,"ERROR: Please specify a channel for all caches, or omit to use the defaults\n"); exit(1);}
    }
    else
    {
        NfHigh=Ntimeslides=Ncache=Nchannel=NfLow=Nifo;

    }
    /* Check for remaining required options */
	if(!(LALInferenceGetProcParamVal(commandLine,"--PSDstart")||LALInferenceGetProcParamVal(commandLine,"--psdstart")) ||
            !(LALInferenceGetProcParamVal(commandLine,"--PSDlength")||LALInferenceGetProcParamVal(commandLine,"--psdlength")) ||!LALInferenceGetProcParamVal(commandLine,"--seglen"))
    {fprintf(stderr,USAGE); return(NULL);}

    /* ET detectors */
    LALDetector dE1,dE2,dE3;
    /* response of the detectors */
    dE1.type = dE2.type = dE3.type = LALDETECTORTYPE_IFODIFF;
    dE1.location[0] = dE2.location[0] = dE3.location[0] = 4.546374099002599e6;
    dE1.location[1] = dE2.location[1] = dE3.location[1] = 8.42989697626334e5;
    dE1.location[2] = dE2.location[2] = dE3.location[2] = 4.378576962409281e6;
    sprintf(dE1.frDetector.name,"ET-1");
    sprintf(dE1.frDetector.prefix,"E1");
    dE1.response[0][0] = 0.166589852497480;
    dE1.response[1][1] = -0.248382035405337;
    dE1.response[2][2] = 0.081792182907857;
    dE1.response[0][1] = dE1.response[1][0] = -0.218849471235102 ;
    dE1.response[0][2] = dE1.response[2][0] = -0.129963871963915;
    dE1.response[1][2] = dE1.response[2][1] = 0.273214957676611;
    sprintf(dE2.frDetector.name,"ET-2");
    sprintf(dE2.frDetector.prefix,"E2");
    dE2.response[0][0] = -0.199221201378560;
    dE2.response[1][1] = 0.423356724499319;
    dE2.response[2][2]=-0.224135523120759;
    dE2.response[0][1] = dE2.response[1][0] = -0.070223802479191;
    dE2.response[0][2] = dE2.response[2][0] = 0.218900453442919;
    dE2.response[1][2] = dE2.response[2][1] = -0.008534697228688;
    sprintf(dE3.frDetector.name,"ET-3");
    sprintf(dE3.frDetector.prefix,"E3");
    dE3.response[0][0] = 0.032631348881079;
    dE3.response[1][1] = -0.174974689093981 ;
    dE3.response[2][2] = 0.142343340212902;
    dE3.response[0][1] = dE3.response[1][0] = 0.289073273714293;
    dE3.response[0][2] = dE3.response[2][0] = -0.088936581479004;
    dE3.response[1][2] = dE3.response[2][1] = -0.264680260447922;

    //TEMPORARY. JUST FOR CHECKING USING SPINSPIRAL PSD
    char **spinspiralPSD=NULL;
    UINT4 NspinspiralPSD = 0;
    if (LALInferenceGetProcParamVal(commandLine, "--spinspiralPSD")) {
        LALInferenceParseCharacterOptionString(LALInferenceGetProcParamVal(commandLine,"--spinspiralPSD")->value,&spinspiralPSD,&NspinspiralPSD);
    }    

    if(LALInferenceGetProcParamVal(commandLine,"--psd")){
        LALInferenceParseCharacterOptionString(LALInferenceGetProcParamVal(commandLine,"--psd")->value,&psds,&Npsd);
    }

    if(LALInferenceGetProcParamVal(commandLine,"--dataseed")){
        procparam=LALInferenceGetProcParamVal(commandLine,"--dataseed");
        dataseed=atoi(procparam->value);
    }

    IFOdata=headIFO=XLALCalloc(sizeof(LALInferenceIFOData),Nifo);
    if(!IFOdata) XLAL_ERROR_NULL(XLAL_ENOMEM);

    if(LALInferenceGetProcParamVal(commandLine,"--injXML"))
    {
        XLALPrintError("ERROR: --injXML option is deprecated. Use --inj and update your scripts\n");
        exit(1);
    }
    procparam=LALInferenceGetProcParamVal(commandLine,"--inj");
    if(procparam){
        SimInspiralTableFromLIGOLw(&injTable,procparam->value,0,0);
        if(!injTable){
            XLALPrintError("Unable to open injection file(LALInferenceReadData) %s\n",procparam->value);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }
        procparam=LALInferenceGetProcParamVal(commandLine,"--event");
        if(procparam) {
            event=atoi(procparam->value);
            while(q<event) {q++; injTable=injTable->next;}
        }
        else if ((procparam=LALInferenceGetProcParamVal(commandLine,"--event-id")))
        {
            while(injTable)
            {
                if(injTable->event_id->id == (UINT4)atoi(procparam->value)) break;
                else injTable=injTable->next;
            }
            if(!injTable){
                fprintf(stderr,"Error, cannot find simulation id %s in injection file\n",procparam->value);
                exit(1);
            }
        }
    }

    procparam=LALInferenceGetProcParamVal(commandLine,"--psdstart");
    if (!procparam) procparam=LALInferenceGetProcParamVal(commandLine,"--PSDstart");
    LALStringToGPS(&status,&GPSstart,procparam->value,&chartmp);
    if(status.statusCode) REPORTSTATUS(&status);

    if(LALInferenceGetProcParamVal(commandLine,"--trigtime")){
        procparam=LALInferenceGetProcParamVal(commandLine,"--trigtime");
        LALStringToGPS(&status,&GPStrig,procparam->value,&chartmp);
    }
    else{
        if(injTable) memcpy(&GPStrig,&(injTable->geocent_end_time),sizeof(GPStrig));
        else if(!LALInferenceGetProcParamVal(commandLine,"--segment-start")){
            XLALPrintError("Error: No trigger time specifed and no injection given \n");
            XLAL_ERROR_NULL(XLAL_EINVAL);
        }
    }
    if(status.statusCode) REPORTSTATUS(&status);
    ppt=LALInferenceGetProcParamVal(commandLine,"--psdlength");
    if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--PSDlength");
    PSDdatalength=atof(ppt->value);
    SegmentLength=atof(LALInferenceGetProcParamVal(commandLine,"--seglen")->value);
    seglen=(size_t)(SegmentLength*SampleRate);
    nSegs=(int)floor(PSDdatalength/SegmentLength);

    CHAR df_argument_name[128];
    REAL8 dof;

    for(i=0;i<Nifo;i++) {
        IFOdata[i].fLow=fLows?atof(fLows[i]):defaultFLow; 
        if(fHighs) IFOdata[i].fHigh=fHighs[i]?atof(fHighs[i]):(SampleRate/2.0-(1.0/SegmentLength));
        else IFOdata[i].fHigh=(SampleRate/2.0-(1.0/SegmentLength));
        strncpy(IFOdata[i].name, IFOnames[i], DETNAMELEN);

        dof=4.0 / M_PI * nSegs; /* Degrees of freedom parameter */
        sprintf(df_argument_name,"--dof-%s",IFOdata[i].name);
        if((ppt=LALInferenceGetProcParamVal(commandLine,df_argument_name)))
            dof=atof(ppt->value);

        IFOdata[i].STDOF = dof;
        XLALPrintInfo("Detector %s will run with %g DOF if Student's T likelihood used.\n",
                IFOdata[i].name, IFOdata[i].STDOF);
    }

    /* Only allocate this array if there weren't channels read in from the command line */
    if(!dataOpts && !Nchannel) channels=XLALCalloc(Nifo,sizeof(char *));
    for(i=0;i<Nifo;i++) {
        if(!dataOpts && !Nchannel) channels[i]=XLALMalloc(VARNAME_MAX);
        IFOdata[i].detector=XLALCalloc(1,sizeof(LALDetector));

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
        /*		if(!strcmp(IFOnames[i],"TAMA")||!strcmp(IFOnames[i],"T1")) {memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexTAMA300DIFF]); continue;}*/

        if(!strcmp(IFOnames[i],"E1")){
            memcpy(IFOdata[i].detector,&dE1,sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"E1:STRAIN"); continue;}
        if(!strcmp(IFOnames[i],"E2")){
            memcpy(IFOdata[i].detector,&dE2,sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"E2:STRAIN"); continue;}
        if(!strcmp(IFOnames[i],"E3")){
            memcpy(IFOdata[i].detector,&dE3,sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"E3:STRAIN"); continue;}
        if(!strcmp(IFOnames[i],"HM1")){
            /* Note, this is a sqrt(2)*7.5-km 3rd gen detector */
            LALFrDetector ETHomestakeFr;
            sprintf(ETHomestakeFr.name,"ET-HomeStake1");
            sprintf(ETHomestakeFr.prefix,"M1");
            /* Location of Homestake Mine vertex is */
            /* 44d21'23.11" N, 103d45'54.71" W */
            ETHomestakeFr.vertexLatitudeRadians = (44.+ 21./60  + 23.11/3600)*LAL_PI/180.0;
            ETHomestakeFr.vertexLongitudeRadians = - (103. +45./60 + 54.71/3600)*LAL_PI/180.0;
            ETHomestakeFr.vertexElevation=0.0;
            ETHomestakeFr.xArmAltitudeRadians=0.0;
            ETHomestakeFr.xArmAzimuthRadians=LAL_PI/2.0;
            ETHomestakeFr.yArmAltitudeRadians=0.0;
            ETHomestakeFr.yArmAzimuthRadians=0.0;
            ETHomestakeFr.xArmMidpoint = ETHomestakeFr.yArmMidpoint = sqrt(2.0)*7.5/2.0;
            IFOdata[i].detector=XLALCalloc(1,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&ETHomestakeFr,LALDETECTORTYPE_IFODIFF);
            printf("Created Homestake Mine ET detector, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"HM2")){
            /* Note, this is a sqrt(2)*7.5-km 3rd gen detector */
            LALFrDetector ETHomestakeFr;
            sprintf(ETHomestakeFr.name,"ET-HomeStake2");
            sprintf(ETHomestakeFr.prefix,"M2");
            /* Location of Homestake Mine vertex is */
            /* 44d21'23.11" N, 103d45'54.71" W */
            ETHomestakeFr.vertexLatitudeRadians = (44.+ 21./60  + 23.11/3600)*LAL_PI/180.0;
            ETHomestakeFr.vertexLongitudeRadians = - (103. +45./60 + 54.71/3600)*LAL_PI/180.0;
            ETHomestakeFr.vertexElevation=0.0;
            ETHomestakeFr.xArmAltitudeRadians=0.0;
            ETHomestakeFr.xArmAzimuthRadians=3.0*LAL_PI/4.0;
            ETHomestakeFr.yArmAltitudeRadians=0.0;
            ETHomestakeFr.yArmAzimuthRadians=LAL_PI/4.0;
            ETHomestakeFr.xArmMidpoint = ETHomestakeFr.yArmMidpoint = sqrt(2.0)*7500./2.0;
            IFOdata[i].detector=XLALCalloc(1,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&ETHomestakeFr,LALDETECTORTYPE_IFODIFF);
            printf("Created Homestake Mine ET detector, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"EM1")){
            LALFrDetector ETmic1;
            sprintf(ETmic1.name,"ET_Michelson_1");
            sprintf(ETmic1.prefix,"F1");
            ETmic1.vertexLatitudeRadians = (43. + 37./60. + 53.0921/3600)*LAL_PI/180.0;
            ETmic1.vertexLongitudeRadians = (10. + 30./60. + 16.1878/3600.)*LAL_PI/180.0;
            ETmic1.vertexElevation = 0.0;
            ETmic1.xArmAltitudeRadians = ETmic1.yArmAltitudeRadians = 0.0;
            ETmic1.xArmAzimuthRadians = LAL_PI/2.0;
            ETmic1.yArmAzimuthRadians = 0.0;
            ETmic1.xArmMidpoint = ETmic1.yArmMidpoint = sqrt(2.0)*7500./2.;
            IFOdata[i].detector=XLALCalloc(1,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&ETmic1,LALDETECTORTYPE_IFODIFF);
            printf("Created ET L-detector 1 (N/E) arms, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"EM2")){
            LALFrDetector ETmic2;
            sprintf(ETmic2.name,"ET_Michelson_2");
            sprintf(ETmic2.prefix,"F2");
            ETmic2.vertexLatitudeRadians = (43. + 37./60. + 53.0921/3600)*LAL_PI/180.0;
            ETmic2.vertexLongitudeRadians = (10. + 30./60. + 16.1878/3600.)*LAL_PI/180.0;
            ETmic2.vertexElevation = 0.0;
            ETmic2.xArmAltitudeRadians = ETmic2.yArmAltitudeRadians = 0.0;
            ETmic2.xArmAzimuthRadians = 3.0*LAL_PI/4.0;
            ETmic2.yArmAzimuthRadians = LAL_PI/4.0;
            ETmic2.xArmMidpoint = ETmic2.yArmMidpoint = sqrt(2.0)*7500./2.;
            IFOdata[i].detector=XLALCalloc(1,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&ETmic2,LALDETECTORTYPE_IFODIFF);
            printf("Created ET L-detector 2 (NE/SE) arms, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"I1")||!strcmp(IFOnames[i],"LIGOIndia")){
            /* Detector in India with 4k arms */
            LALFrDetector LIGOIndiaFr;
            sprintf(LIGOIndiaFr.name,"LIGO_India");
            sprintf(LIGOIndiaFr.prefix,"I1");
            /* Location of India site is */
            /* 14d14' N 76d26' E */
            LIGOIndiaFr.vertexLatitudeRadians = (14. + 14./60.)*LAL_PI/180.0;
            LIGOIndiaFr.vertexLongitudeRadians = (76. + 26./60.)*LAL_PI/180.0;
            LIGOIndiaFr.vertexElevation = 0.0;
            LIGOIndiaFr.xArmAltitudeRadians = 0.0;
            LIGOIndiaFr.yArmAltitudeRadians = 0.0;
            LIGOIndiaFr.yArmMidpoint = 2000.;
            LIGOIndiaFr.xArmMidpoint = 2000.;
            LIGOIndiaFr.xArmAzimuthRadians = LAL_PI/2.;
            LIGOIndiaFr.yArmAzimuthRadians = 0.;
            IFOdata[i].detector=XLALMalloc(sizeof(LALDetector));
            memset(IFOdata[i].detector,0,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&LIGOIndiaFr,LALDETECTORTYPE_IFODIFF);
            printf("Created LIGO India Detector, location %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("Detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"A1")||!strcmp(IFOnames[i],"LIGOSouth")){
            /* Construct a detector at AIGO with 4k arms */
            LALFrDetector LIGOSouthFr;
            sprintf(LIGOSouthFr.name,"LIGO-South");
            sprintf(LIGOSouthFr.prefix,"A1");
            /* Location of the AIGO detector vertex is */
            /* 31d21'27.56" S, 115d42'50.34"E */
            LIGOSouthFr.vertexLatitudeRadians = - (31. + 21./60. + 27.56/3600.)*LAL_PI/180.0;
            LIGOSouthFr.vertexLongitudeRadians = (115. + 42./60. + 50.34/3600.)*LAL_PI/180.0;
            LIGOSouthFr.vertexElevation=0.0;
            LIGOSouthFr.xArmAltitudeRadians=0.0;
            LIGOSouthFr.xArmAzimuthRadians=AIGOang+LAL_PI/2.;
            LIGOSouthFr.yArmAltitudeRadians=0.0;
            LIGOSouthFr.yArmAzimuthRadians=AIGOang;
            LIGOSouthFr.xArmMidpoint=2000.;
            LIGOSouthFr.yArmMidpoint=2000.;
            IFOdata[i].detector=XLALMalloc(sizeof(LALDetector));
            memset(IFOdata[i].detector,0,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&LIGOSouthFr,LALDETECTORTYPE_IFODIFF);
            printf("Created LIGO South detector, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("Detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        if(!strcmp(IFOnames[i],"J1")||!strcmp(IFOnames[i],"LCGT")){
            /* Construct the LCGT telescope */
            REAL8 LCGTangle=19.0*(LAL_PI/180.0);
            LALFrDetector LCGTFr;
            sprintf(LCGTFr.name,"LCGT");
            sprintf(LCGTFr.prefix,"J1");
            LCGTFr.vertexLatitudeRadians  = 36.25 * LAL_PI/180.0;
            LCGTFr.vertexLongitudeRadians = (137.18 * LAL_PI/180.0);
            LCGTFr.vertexElevation=0.0;
            LCGTFr.xArmAltitudeRadians=0.0;
            LCGTFr.xArmAzimuthRadians=LCGTangle+LAL_PI/2.;
            LCGTFr.yArmAltitudeRadians=0.0;
            LCGTFr.yArmAzimuthRadians=LCGTangle;
            LCGTFr.xArmMidpoint=1500.;
            LCGTFr.yArmMidpoint=1500.;
            IFOdata[i].detector=XLALMalloc(sizeof(LALDetector));
            memset(IFOdata[i].detector,0,sizeof(LALDetector));
            XLALCreateDetector(IFOdata[i].detector,&LCGTFr,LALDETECTORTYPE_IFODIFF);
            printf("Created LCGT telescope, location: %lf, %lf, %lf\n",IFOdata[i].detector->location[0],IFOdata[i].detector->location[1],IFOdata[i].detector->location[2]);
            printf("Detector tensor:\n");
            for(int jdx=0;jdx<3;jdx++){
                for(j=0;j<3;j++) printf("%f ",IFOdata[i].detector->response[jdx][j]);
                printf("\n");
            }
            continue;
        }
        fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 H2 L1 V1 GEO A1 J1 I1 E1 E2 E3 HM1 HM2 EM1 EM2\n",IFOnames[i]); exit(-1);
    }

    /* Set up FFT structures and window */
    for (i=0;i<Nifo;i++){
        /* Create FFT plans (flag 1 to measure performance) */
        IFOdata[i].timeToFreqFFTPlan = XLALCreateForwardREAL8FFTPlan((UINT4) seglen, 1 );
        if(!IFOdata[i].timeToFreqFFTPlan) XLAL_ERROR_NULL(XLAL_EFUNC);
        IFOdata[i].freqToTimeFFTPlan = XLALCreateReverseREAL8FFTPlan((UINT4) seglen, 1 );
        if(!IFOdata[i].freqToTimeFFTPlan) XLAL_ERROR_NULL(XLAL_EFUNC);		
        /* Setup windows */
        IFOdata[i].window=XLALCreateTukeyREAL8Window(seglen,(REAL8)2.0*padding*SampleRate/(REAL8)seglen);
        if(!IFOdata[i].window) XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    if(!(ppt=LALInferenceGetProcParamVal(commandLine,"--segment-start")))
    {
        /* Trigger time = 2 seconds before end of segment (was 1 second, but Common Inputs for The Events are -6 +2*/
        memcpy(&segStart,&GPStrig,sizeof(LIGOTimeGPS));
        XLALGPSAdd(&segStart,-SegmentLength+2);
    }
    else
    {
        /* Segment starts at given time */
        REAL8 segstartR8 = atof(ppt->value);
        XLALGPSSetREAL8(&segStart,segstartR8);
    }


    /* Read the PSD data */
    for(i=0;i<Nifo;i++) {
        memcpy(&(IFOdata[i].epoch),&segStart,sizeof(LIGOTimeGPS));
        /* Check to see if an interpolation file is specified */
        interpFlag=0;
        interp=NULL;
        if( (globFrames)?0:strstr(caches[i],"interp:")==caches[i]){
          /* Extract the file name */
         char *interpfilename=&(caches[i][7]);
         printf("Looking for interpolation file %s\n",interpfilename);
         interpFlag=1;
         interp=interpFromFile(interpfilename);
        }
        /* Check if fake data is requested */
       if( (globFrames)?0:(interpFlag || (!(strcmp(caches[i],"LALLIGO") && strcmp(caches[i],"LALVirgo") && strcmp(caches[i],"LALGEO") && strcmp(caches[i],"LALEGO") && strcmp(caches[i],"LALSimLIGO") && strcmp(caches[i],"LALSimAdLIGO") && strcmp(caches[i],"LALSimVirgo") && strcmp(caches[i],"LALSimAdVirgo") && strcmp(caches[i],"LALAdLIGO")))))
        {
            //FakeFlag=1; - set but not used
            if (!LALInferenceGetProcParamVal(commandLine,"--dataseed")){
                fprintf(stderr,"Error: You need to specify a dataseed when generating data with --dataseed <number>.\n\
                        (--dataseed 0 uses a non-reproducible number from the system clock, and no parallel run is then possible.)\n" );
                exit(-1);
            }
            /* Offset the seed in a way that depends uniquely on the IFO name */
            int ifo_salt=0;
            ifo_salt+=(int)IFOnames[i][0]+(int)IFOnames[i][1];
            datarandparam=XLALCreateRandomParams(dataseed?dataseed+(int)ifo_salt:dataseed);
            if(!datarandparam) XLAL_ERROR_NULL(XLAL_EFUNC);
            IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)
                XLALCreateREAL8FrequencySeries("spectrum",&GPSstart,0.0,
                        (REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
            if(!IFOdata[i].oneSidedNoisePowerSpectrum) XLAL_ERROR_NULL(XLAL_EFUNC);
            
            int LALSimPsd=0;
            /* Selection of the noise curve */
            if(!strcmp(caches[i],"LALLIGO")) {PSD = &LALLIGOIPsd; scalefactor=9E-46;}
            if(!strcmp(caches[i],"LALVirgo")) {PSD = &LALVIRGOPsd; scalefactor=1.0;}
            if(!strcmp(caches[i],"LALGEO")) {PSD = &LALGEOPsd; scalefactor=1E-46;}
            if(!strcmp(caches[i],"LALEGO")) {PSD = &LALEGOPsd; scalefactor=1.0;}
            if(!strcmp(caches[i],"LALAdLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 1E-49;}
            if(!strcmp(caches[i],"LALSimLIGO")) {XLALSimNoisePSD(IFOdata[i].oneSidedNoisePowerSpectrum,IFOdata[i].fLow,XLALSimNoisePSDiLIGOSRD ) ; LALSimPsd=1;}
            if(!strcmp(caches[i],"LALSimVirgo")) {XLALSimNoisePSD(IFOdata[i].oneSidedNoisePowerSpectrum,IFOdata[i].fLow,XLALSimNoisePSDVirgo ); LALSimPsd=1;}
            if(!strcmp(caches[i],"LALSimAdLIGO")) {XLALSimNoisePSD(IFOdata[i].oneSidedNoisePowerSpectrum,IFOdata[i].fLow,XLALSimNoisePSDaLIGOZeroDetHighPower ) ;LALSimPsd=1;}
            if(!strcmp(caches[i],"LALSimAdVirgo")) {XLALSimNoisePSD(IFOdata[i].oneSidedNoisePowerSpectrum,IFOdata[i].fLow,XLALSimNoisePSDAdvVirgo) ;LALSimPsd=1;}
            if(interpFlag) {PSD=NULL; scalefactor=1.0;}
            //if(!strcmp(caches[i],"LAL2kLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 36E-46;}
            if(PSD==NULL && !(interpFlag|| LALSimPsd)) {fprintf(stderr,"Error: unknown simulated PSD: %s\n",caches[i]); exit(-1);}

            if(LALSimPsd==0){
                for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
                {
                    MetaNoiseFunc(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,interp,PSD);
                    //PSD(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF);
                    IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]*=scalefactor;
                }
            }
            
            IFOdata[i].freqData = (COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("stilde",&segStart,0.0,IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,&lalDimensionlessUnit,seglen/2 +1);
            if(!IFOdata[i].freqData) XLAL_ERROR_NULL(XLAL_EFUNC);

            /* Create the fake data */
            int j_Lo = (int) IFOdata[i].fLow/IFOdata[i].freqData->deltaF;
            if(LALInferenceGetProcParamVal(commandLine,"--0noise")){
                for(j=j_Lo;j<IFOdata[i].freqData->data->length;j++){
                    IFOdata[i].freqData->data->data[j] = 0.0;
                }
            } else {
                for(j=j_Lo;j<IFOdata[i].freqData->data->length;j++){
                    IFOdata[i].freqData->data->data[j] = crect(
                      XLALNormalDeviate(datarandparam)*(0.5*sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]/IFOdata[i].freqData->deltaF)),
                      XLALNormalDeviate(datarandparam)*(0.5*sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]/IFOdata[i].freqData->deltaF))
                      );
                }
            }
            IFOdata[i].freqData->data->data[0] = 0;
            const char timename[]="timeData";
            IFOdata[i].timeData=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries(timename,&segStart,0.0,(REAL8)1.0/SampleRate,&lalDimensionlessUnit,(size_t)seglen);
            if(!IFOdata[i].timeData) XLAL_ERROR_NULL(XLAL_EFUNC);
            XLALREAL8FreqTimeFFT(IFOdata[i].timeData,IFOdata[i].freqData,IFOdata[i].freqToTimeFFTPlan);
            if(*XLALGetErrnoPtr()) printf("XLErr: %s\n",XLALErrorString(*XLALGetErrnoPtr()));
            XLALDestroyRandomParams(datarandparam);
        }
        else{ /* Not using fake data, load the data from a cache file */

            LALCache *cache=NULL;
            if(!globFrames)
            {
                cache  = XLALCacheImport( caches[i] );
                int err;
                err = *XLALGetErrnoPtr();
                if(cache==NULL) {fprintf(stderr,"ERROR: Unable to import cache file \"%s\",\n       XLALError: \"%s\".\n",caches[i], XLALErrorString(err)); exit(-1);}
            }
            else
            {
                printf("Looking for frames for %s in PWD\n",IFOnames[i]);
                cache= GlobFramesPWD(IFOnames[i]);

            }
            if(!cache) {fprintf(stderr,"ERROR: Cannot find any frame data!\n"); exit(1);}
            if (LALInferenceGetProcParamVal(commandLine, "--psd")){
                interp=NULL;
                char *interpfilename=&(psds[i][0]);
                fprintf(stderr,"Reading PSD for %s using %s\n",IFOnames[i],interpfilename);
                printf("Looking for psd interpolation file %s\n",interpfilename);
                interp=interpFromFile(interpfilename);
                IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)
                    XLALCreateREAL8FrequencySeries("spectrum",&GPSstart,0.0,
                            (REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
                if(!IFOdata[i].oneSidedNoisePowerSpectrum) XLAL_ERROR_NULL(XLAL_EFUNC);
                for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
                {
                    MetaNoiseFunc(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,interp,NULL);
                }
            }else{
                fprintf(stderr,"Estimating PSD for %s using %i segments of %i samples (%lfs)\n",IFOnames[i],nSegs,(int)seglen,SegmentLength);
                PSDtimeSeries=readTseries(cache,channels[i],GPSstart,PSDdatalength);
                if(!PSDtimeSeries) {XLALPrintError("Error reading PSD data for %s\n",IFOnames[i]); XLAL_ERROR_NULL(XLAL_EFUNC);}
                XLALResampleREAL8TimeSeries(PSDtimeSeries,1.0/SampleRate);
                PSDtimeSeries=(REAL8TimeSeries *)XLALShrinkREAL8TimeSeries(PSDtimeSeries,(size_t) 0, (size_t) seglen*nSegs);
                if(!PSDtimeSeries) {
                    fprintf(stderr,"ERROR while estimating PSD for %s\n",IFOnames[i]);
                    XLAL_ERROR_NULL(XLAL_EFUNC);
                }
                IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("spectrum",&PSDtimeSeries->epoch,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
                if(!IFOdata[i].oneSidedNoisePowerSpectrum) XLAL_ERROR_NULL(XLAL_EFUNC);
                if (LALInferenceGetProcParamVal(commandLine, "--PSDwelch"))
                    XLALREAL8AverageSpectrumWelch(IFOdata[i].oneSidedNoisePowerSpectrum ,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan);
                else
                    XLALREAL8AverageSpectrumMedian(IFOdata[i].oneSidedNoisePowerSpectrum ,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan);	

                if(LALInferenceGetProcParamVal(commandLine, "--binFit")) {

                    LIGOTimeGPS GPStime=segStart;

                    const UINT4 nameLength=256;
                    char filename[nameLength];

                    snprintf(filename, nameLength, "%s-BinFitLines.dat", IFOdata[i].name);

                    printf("Running PSD bin fitting... ");
                    LALInferenceAverageSpectrumBinFit(IFOdata[i].oneSidedNoisePowerSpectrum ,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan,filename,GPStime);
                    printf("completed!\n");
                }

                if (LALInferenceGetProcParamVal(commandLine, "--chisquaredlines")){

                    double deltaF = IFOdata[i].oneSidedNoisePowerSpectrum->deltaF;
                    int lengthF = IFOdata[i].oneSidedNoisePowerSpectrum->data->length;

                    REAL8 *pvalues;
                    pvalues = XLALMalloc( lengthF * sizeof( *pvalues ) );

                    printf("Running chi-squared tests... ");
                    LALInferenceRemoveLinesChiSquared(IFOdata[i].oneSidedNoisePowerSpectrum,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan,pvalues);
                    printf("completed!\n");

                    const UINT4 nameLength=256;
                    char filename[nameLength];
                    FILE *out;

                    double lines_width;
                    ppt = LALInferenceGetProcParamVal(commandLine, "--chisquaredlinesWidth");
                    if(ppt) lines_width = atof(ppt->value);
                    else lines_width = deltaF;

                    double lines_threshold;
                    ppt = LALInferenceGetProcParamVal(commandLine, "--chisquaredlinesThreshold");
                    if(ppt) lines_threshold = atof(ppt->value);
                    else lines_threshold = 2*pow(10.0,-14.0);

                    printf("Using chi squared threshold of %g\n",lines_threshold);

                    snprintf(filename, nameLength, "%s-ChiSquaredLines.dat", IFOdata[i].name);
                    out = fopen(filename, "w");
                    for (int k = 0; k < lengthF; ++k ) {
                        if (pvalues[k] < lines_threshold) {
                            fprintf(out,"%g %g\n",((double) k) * deltaF,lines_width);
                        }
                    }
                    fclose(out);

                    snprintf(filename, nameLength, "%s-ChiSquaredLines-pvalues.dat", IFOdata[i].name);
                    out = fopen(filename, "w");
                    for (int k = 0; k < lengthF; ++k ) {
                        fprintf(out,"%g %g\n",((double) k) * deltaF,pvalues[k]);
                    }
                    fclose(out);

                }

                if (LALInferenceGetProcParamVal(commandLine, "--KSlines")){

                    double deltaF = IFOdata[i].oneSidedNoisePowerSpectrum->deltaF;
                    int lengthF = IFOdata[i].oneSidedNoisePowerSpectrum->data->length;

                    REAL8 *pvalues;
                    pvalues = XLALMalloc( lengthF * sizeof( *pvalues ) );

                    printf("Running KS tests... ");
                    LALInferenceRemoveLinesKS(IFOdata[i].oneSidedNoisePowerSpectrum,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan,pvalues);
                    printf("completed!\n");

                    const UINT4 nameLength=256;
                    char filename[nameLength];
                    FILE *out;

                    double lines_width;
                    ppt = LALInferenceGetProcParamVal(commandLine, "--KSlinesWidth");
                    if(ppt) lines_width = atof(ppt->value);
                    else lines_width = deltaF;

                    double lines_threshold;
                    ppt = LALInferenceGetProcParamVal(commandLine, "--KSlinesThreshold");
                    if(ppt) lines_threshold = atof(ppt->value);
                    else lines_threshold = 0.134558;

                    printf("Using KS threshold of %g\n",lines_threshold);

                    snprintf(filename, nameLength, "%s-KSLines.dat", IFOdata[i].name);
                    out = fopen(filename, "w");
                    for (int k = 0; k < lengthF; ++k ) {
                        if (pvalues[k] < lines_threshold) {
                            fprintf(out,"%g %g\n",((double) k) * deltaF,lines_width);
                        }
                    }
                    fclose(out);

                    snprintf(filename, nameLength, "%s-KSLines-pvalues.dat", IFOdata[i].name);
                    out = fopen(filename, "w");
                    for (int k = 0; k < lengthF; ++k ) {
                        fprintf(out,"%g %g\n",((double) k) * deltaF,pvalues[k]);
                    }
                    fclose(out);

                }

                if (LALInferenceGetProcParamVal(commandLine, "--powerlawlines")){

                    double deltaF = IFOdata[i].oneSidedNoisePowerSpectrum->deltaF;
                    int lengthF = IFOdata[i].oneSidedNoisePowerSpectrum->data->length;

                    REAL8 *pvalues;
                    pvalues = XLALMalloc( lengthF * sizeof( *pvalues ) );

                    printf("Running power law tests... ");
                    LALInferenceRemoveLinesPowerLaw(IFOdata[i].oneSidedNoisePowerSpectrum,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan,pvalues);
                    printf("completed!\n");

                    const UINT4 nameLength=256;
                    char filename[nameLength];
                    FILE *out;

                    double lines_width;
                    ppt = LALInferenceGetProcParamVal(commandLine, "--powerlawlinesWidth");
                    if(ppt) lines_width = atof(ppt->value);
                    else lines_width = deltaF;

                    double lines_threshold;
                    ppt = LALInferenceGetProcParamVal(commandLine, "--powerlawlinesThreshold");
                    if(ppt) lines_threshold = atof(ppt->value);
                    else lines_threshold = 0.7197370;

                    printf("Using power law threshold of %g\n",lines_threshold);

                    snprintf(filename, nameLength, "%s-PowerLawLines.dat", IFOdata[i].name);
                    out = fopen(filename, "w");
                    for (int k = 0; k < lengthF; ++k ) {
                        if (pvalues[k] < lines_threshold) {
                            fprintf(out,"%g %g\n",((double) k) * deltaF,lines_width);
                        }
                    }
                    fclose(out);

                    snprintf(filename, nameLength, "%s-PowerLawLines-pvalues.dat", IFOdata[i].name);
                    out = fopen(filename, "w");
                    for (int k = 0; k < lengthF; ++k ) {
                        fprintf(out,"%g %g\n",((double) k) * deltaF,pvalues[k]);
                    }
                    fclose(out);

                }

                if (LALInferenceGetProcParamVal(commandLine, "--xcorrbands")){

                    //double deltaF = IFOdata[i].oneSidedNoisePowerSpectrum->deltaF;
                    int lengthF = IFOdata[i].oneSidedNoisePowerSpectrum->data->length;

                    REAL8 *pvalues;
                    pvalues = XLALMalloc( lengthF * sizeof( *pvalues ) );

                    const UINT4 nameLength=256;
                    char filename[nameLength];
                    FILE *out;

                    snprintf(filename, nameLength, "%s-XCorrVals.dat", IFOdata[i].name);

                    printf("Running xcorr tests... ");
                    LALInferenceXCorrBands(IFOdata[i].oneSidedNoisePowerSpectrum,PSDtimeSeries, seglen, (UINT4)seglen, IFOdata[i].window, IFOdata[i].timeToFreqFFTPlan,pvalues,filename);
                    printf("completed!\n");

                    snprintf(filename, nameLength, "%s-XCorrBands.dat", IFOdata[i].name);
                    out = fopen(filename, "w");
                    /*
                    for (int k = 0; k < lengthF; ++k ) {
                        if (pvalues[k] < 0.001) {
                            fprintf(out,"%g %g\n",((double) k) * deltaF,lines_width);
                        }
                    }
                    */
                    fprintf(out,"%g %g\n",10.0,75.0);
                    fprintf(out,"%g %g\n",16.0,40.0);
                    fprintf(out,"%g %g\n",40.0,330.0);
                    fclose(out);

                }

                XLALDestroyREAL8TimeSeries(PSDtimeSeries);
            }

            /* Read the data segment */
            LIGOTimeGPS truesegstart=segStart;
            if(Ntimeslides) {
                REAL4 deltaT=-atof(timeslides[i]);
                XLALGPSAdd(&segStart, deltaT);
                fprintf(stderr,"Slid %s by %f s from %10.10lf to %10.10lf\n",IFOnames[i],deltaT,truesegstart.gpsSeconds+1e-9*truesegstart.gpsNanoSeconds,segStart.gpsSeconds+1e-9*segStart.gpsNanoSeconds);
            }
            IFOdata[i].timeData=readTseries(cache,channels[i],segStart,SegmentLength);
            segStart=truesegstart;
            if(Ntimeslides) IFOdata[i].timeData->epoch=truesegstart;
            /* FILE *out; */
            /* char fileName[256]; */
            /* snprintf(fileName, 256, "readTimeData-%d.dat", i); */
            /* out = fopen(fileName, "w"); */
            /* for (j = 0; j < IFOdata[i].timeData->data->length; j++) { */
            /*   fprintf(out, "%g %g\n", j*IFOdata[i].timeData->deltaT, IFOdata[i].timeData->data->data[j]); */
            /* } */
            /* fclose(out); */

            if(!IFOdata[i].timeData) {
                XLALPrintError("Error reading segment data for %s at %i\n",IFOnames[i],segStart.gpsSeconds);
                XLAL_ERROR_NULL(XLAL_EFUNC);
            }
            XLALResampleREAL8TimeSeries(IFOdata[i].timeData,1.0/SampleRate);	 
            if(!IFOdata[i].timeData) {XLALPrintError("Error reading segment data for %s\n",IFOnames[i]); XLAL_ERROR_NULL(XLAL_EFUNC);}
            IFOdata[i].freqData=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("freqData",&(IFOdata[i].timeData->epoch),0.0,1.0/SegmentLength,&lalDimensionlessUnit,seglen/2+1);
            if(!IFOdata[i].freqData) XLAL_ERROR_NULL(XLAL_EFUNC);
            IFOdata[i].windowedTimeData=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("windowed time data",&(IFOdata[i].timeData->epoch),0.0,1.0/SampleRate,&lalDimensionlessUnit,seglen);
            if(!IFOdata[i].windowedTimeData) XLAL_ERROR_NULL(XLAL_EFUNC);
            XLALDDVectorMultiply(IFOdata[i].windowedTimeData->data,IFOdata[i].timeData->data,IFOdata[i].window->data);
            XLALREAL8TimeFreqFFT(IFOdata[i].freqData,IFOdata[i].windowedTimeData,IFOdata[i].timeToFreqFFTPlan);

            for(j=0;j<IFOdata[i].freqData->data->length;j++){
                IFOdata[i].freqData->data->data[j] /= sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
                IFOdata[i].windowedTimeData->data->data[j] /= sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
            }

        XLALDestroyCache(cache); // Clean up cache
        } /* End of data reading process */

        //		/* Now that the PSD is set up, make the TDW. */
        //    IFOdata[i].timeDomainNoiseWeights = 
        //                  (REAL8TimeSeries *)XLALCreateREAL8TimeSeries("time domain weights", 
        //                                                               &(IFOdata[i].oneSidedNoisePowerSpectrum->epoch),
        //                                                               0.0,
        //                                                               1.0/SampleRate,
        //                                                               &lalDimensionlessUnit,
        //                                                               seglen);
        //		if(!IFOdata[i].timeDomainNoiseWeights) XLAL_ERROR_NULL(XLAL_EFUNC);
        //		LALInferencePSDToTDW(IFOdata[i].timeDomainNoiseWeights, IFOdata[i].oneSidedNoisePowerSpectrum, IFOdata[i].freqToTimeFFTPlan,
        //                         IFOdata[i].fLow, IFOdata[i].fHigh);

        makeWhiteData(&(IFOdata[i]));

      /* Store ASD of noise spectrum to whiten glitch model */
      IFOdata[i].noiseASD=(REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("asd",&GPSstart,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
      for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
        IFOdata[i].noiseASD->data->data[j]=sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]);

        if (LALInferenceGetProcParamVal(commandLine, "--spinspiralPSD")) {
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
            // tempPSD=XLALCalloc(sizeof(REAL8),templen+1);
            // tempfreq=XLALCalloc(sizeof(REAL8),templen+1);

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

        if (LALInferenceGetProcParamVal(commandLine, "--data-dump")) {
            pptdatadump=LALInferenceGetProcParamVal(commandLine,"--data-dump");
            const UINT4 nameLength=FILENAME_MAX;
            char filename[nameLength];
            FILE *out;
            ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");
            if(ppt) {
            	snprintf(filename, nameLength, "%s-%s-PSD.dat", ppt->value, IFOdata[i].name);
            }
            else if(strcmp(pptdatadump->value,"")) {
              snprintf(filename, nameLength, "%s/%s-PSD.dat", pptdatadump->value, IFOdata[i].name);
            }
            else
                snprintf(filename, nameLength, "%s-PSD.dat", IFOdata[i].name);
            out = fopen(filename, "w");
            for (j = 0; j < IFOdata[i].oneSidedNoisePowerSpectrum->data->length; j++) {
                REAL8 f = IFOdata[i].oneSidedNoisePowerSpectrum->deltaF*j;
                REAL8 psd = IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j];

                fprintf(out, "%g %g\n", f, psd);
            }
            fclose(out);
          
            if(ppt) {
              snprintf(filename, nameLength, "%s-%s-timeData.dat", ppt->value, IFOdata[i].name);
            }
            else if(strcmp(pptdatadump->value,"")) {
              snprintf(filename, nameLength, "%s/%s-timeData.dat", pptdatadump->value, IFOdata[i].name);
            }
            else
              snprintf(filename, nameLength, "%s-timeData.dat", IFOdata[i].name);
            out = fopen(filename, "w");
            for (j = 0; j < IFOdata[i].timeData->data->length; j++) {
                REAL8 t = XLALGPSGetREAL8(&(IFOdata[i].timeData->epoch)) + 
                    j * IFOdata[i].timeData->deltaT;
                REAL8 d = IFOdata[i].timeData->data->data[j];

                fprintf(out, "%.6f %g\n", t, d);
            }
            fclose(out);
          
            ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");
            if(ppt) {
              snprintf(filename, nameLength, "%s-%s-freqData.dat", ppt->value, IFOdata[i].name);
            }
            else if(strcmp(pptdatadump->value,"")) {
              snprintf(filename, nameLength, "%s/%s-freqData.dat", pptdatadump->value, IFOdata[i].name);
            }
            else
              snprintf(filename, nameLength, "%s-freqData.dat", IFOdata[i].name);
            out = fopen(filename, "w");
            for (j = 0; j < IFOdata[i].freqData->data->length; j++) {
                REAL8 f = IFOdata[i].freqData->deltaF * j;
                REAL8 dre = creal(IFOdata[i].freqData->data->data[j]);
                REAL8 dim = cimag(IFOdata[i].freqData->data->data[j]);

                fprintf(out, "%g %g %g\n", f, dre, dim);
            }
            fclose(out);

        }
    }

    for (i=0;i<Nifo;i++) IFOdata[i].SNR=0.0; //SNR of the injection ONLY IF INJECTION. Set to 0.0 by default.

    for (i=0;i<Nifo-1;i++) IFOdata[i].next=&(IFOdata[i+1]);

    for(i=0;i<Nifo;i++) {
        if(channels) if(channels[i]) XLALFree(channels[i]);
        if(caches) if(caches[i]) XLALFree(caches[i]);
        if(IFOnames) if(IFOnames[i]) XLALFree(IFOnames[i]);
        if(fLows) if(fLows[i]) XLALFree(fLows[i]);
        if(fHighs) if(fHighs[i]) XLALFree(fHighs[i]);
    }
    if(channels) XLALFree(channels);
    if(caches) XLALFree(caches);
    if(IFOnames) XLALFree(IFOnames);
    if(fLows) XLALFree(fLows);
    if(fHighs) XLALFree(fHighs);

    LALSimInspiralWaveformCache *cache=XLALCreateSimInspiralWaveformCache();
    for(i=0;i<Nifo;i++) IFOdata[i].waveformCache=cache;

    return headIFO;
}

static void makeWhiteData(LALInferenceIFOData *IFOdata) {
  REAL8 deltaF = IFOdata->freqData->deltaF;
  REAL8 deltaT = IFOdata->timeData->deltaT;

  IFOdata->whiteFreqData = 
    XLALCreateCOMPLEX16FrequencySeries("whitened frequency data", 
                                       &(IFOdata->freqData->epoch),
                                       0.0,
                                       deltaF,
                                       &lalDimensionlessUnit,
                                       IFOdata->freqData->data->length);
	if(!IFOdata->whiteFreqData) XLAL_ERROR_VOID(XLAL_EFUNC);
  IFOdata->whiteTimeData = 
    XLALCreateREAL8TimeSeries("whitened time data",
                              &(IFOdata->timeData->epoch),
                              0.0,
                              deltaT,
                              &lalDimensionlessUnit,
                              IFOdata->timeData->data->length);
	if(!IFOdata->whiteTimeData) XLAL_ERROR_VOID(XLAL_EFUNC);

  REAL8 iLow = IFOdata->fLow / deltaF;
  REAL8 iHighDefaultCut = 0.95 * IFOdata->freqData->data->length;
  REAL8 iHighFromFHigh = IFOdata->fHigh / deltaF;
  REAL8 iHigh = (iHighDefaultCut < iHighFromFHigh ? iHighDefaultCut : iHighFromFHigh);
  REAL8 windowSquareSum = 0.0;

  UINT4 i;

  for (i = 0; i < IFOdata->freqData->data->length; i++) {
    IFOdata->whiteFreqData->data->data[i] = IFOdata->freqData->data->data[i] / IFOdata->oneSidedNoisePowerSpectrum->data->data[i];
		
    if (i == 0) {
      /* Cut off the average trend in the data. */
      IFOdata->whiteFreqData->data->data[i] = 0.0;
    }
    if (i <= iLow) {
      /* Need to taper to implement the fLow cutoff.  Tukey window
			 that starts at zero, and reaches 100% at fLow. */
      REAL8 weight = 0.5*(1.0 + cos(M_PI*(i-iLow)/iLow)); /* Starts at -Pi, runs to zero at iLow. */
			
      IFOdata->whiteFreqData->data->data[i] *= weight;
			
      windowSquareSum += weight*weight;
    } else if (i >= iHigh) {
      /* Also taper at high freq end, Tukey window that starts at 100%
			 at fHigh, then drops to zero at Nyquist.  Except that we
			 always taper at least 5% of the data at high freq to avoid a
			 sharp edge in freq space there. */
      REAL8 NWind = IFOdata->whiteFreqData->data->length - iHigh;
      REAL8 weight = 0.5*(1.0 + cos(M_PI*(i-iHigh)/NWind)); /* Starts at 0, runs to Pi at i = length */
			
      IFOdata->whiteFreqData->data->data[i] *= weight;
			
      windowSquareSum += weight*weight;
    } else {
      windowSquareSum += 1.0;
    }
  }
	
  REAL8 norm = sqrt(IFOdata->whiteFreqData->data->length / windowSquareSum);
  for (i = 0; i < IFOdata->whiteFreqData->data->length; i++) {
    IFOdata->whiteFreqData->data->data[i] *= norm;
  }
	
  XLALREAL8FreqTimeFFT(IFOdata->whiteTimeData, IFOdata->whiteFreqData, IFOdata->freqToTimeFFTPlan);
}

void LALInferenceInjectInspiralSignal(LALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine)
{
	LALStatus status;
	memset(&status,0,sizeof(status));
	SimInspiralTable *injTable=NULL;
    SimInspiralTable *injEvent=NULL;
	UINT4 Ninj=0;
	UINT4 event=0;
	UINT4 i=0,j=0;
    REAL8 responseScale=1.0;
	//CoherentGW InjectGW;
	//PPNParamStruc InjParams;
	LIGOTimeGPS injstart;
	REAL8 SNR=0,NetworkSNR=0;
	DetectorResponse det;
	memset(&injstart,0,sizeof(LIGOTimeGPS));
	//memset(&InjParams,0,sizeof(PPNParamStruc));
	COMPLEX16FrequencySeries *injF=NULL;
	FILE *rawWaveform=NULL;
	ProcessParamsTable *ppt=NULL;
	REAL8 bufferLength = 2048.0; /* Default length of buffer for injections (seconds) */
	UINT4 bufferN=0;
	LIGOTimeGPS bufferStart;

	
	LALInferenceIFOData *thisData=IFOdata->next;
	REAL8 minFlow=IFOdata->fLow;
	REAL8 MindeltaT=IFOdata->timeData->deltaT;
  REAL8 InjSampleRate=1.0/MindeltaT;
	REAL4TimeSeries *injectionBuffer=NULL;
  REAL8 padding=0.4; //default, set in LALInferenceReadData()
	
  
	while(thisData){
          minFlow   = minFlow>thisData->fLow ? thisData->fLow : minFlow;
          MindeltaT = MindeltaT>thisData->timeData->deltaT ? thisData->timeData->deltaT : MindeltaT;
          thisData  = thisData->next;
	}
	thisData=IFOdata;
	//InjParams.deltaT = MindeltaT;
	//InjParams.fStartIn=(REAL4)minFlow;

	if(!LALInferenceGetProcParamVal(commandLine,"--inj")) {fprintf(stdout,"No injection file specified, not injecting\n"); return;}
	if(LALInferenceGetProcParamVal(commandLine,"--event")){
    event= atoi(LALInferenceGetProcParamVal(commandLine,"--event")->value);
    fprintf(stdout,"Injecting event %d\n",event);
	}
        if(LALInferenceGetProcParamVal(commandLine,"--snrpath")){
                ppt = LALInferenceGetProcParamVal(commandLine,"--snrpath");
		SNRpath = XLALCalloc(strlen(ppt->value)+1,sizeof(char));
		memcpy(SNRpath,ppt->value,strlen(ppt->value)+1);
                fprintf(stdout,"Writing SNRs in %s\n",SNRpath)     ;

	}
	Ninj=SimInspiralTableFromLIGOLw(&injTable,LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
	REPORTSTATUS(&status);
	printf("Ninj %d\n", Ninj);
	if(Ninj<=event){
          fprintf(stderr,"Error reading event %d from %s\n",event,LALInferenceGetProcParamVal(commandLine,"--inj")->value);
          exit(1);
        }
	while(i<event) {i++; injTable = injTable->next;} /* Select event */
	injEvent = injTable;
	injEvent->next = NULL;
  enforce_m1_larger_m2(injEvent);
	//memset(&InjectGW,0,sizeof(InjectGW));
	Approximant injapprox;
	injapprox = XLALGetApproximantFromString(injTable->waveform);
        if( (int) injapprox == XLAL_FAILURE)
          ABORTXLAL(&status);
	printf("Injecting approximant %i: %s\n", injapprox, injTable->waveform);
	REPORTSTATUS(&status);
	//LALGenerateInspiral(&status,&InjectGW,injTable,&InjParams);
	//if(status.statusCode!=0) {fprintf(stderr,"Error generating injection!\n"); REPORTSTATUS(&status); }

	/* Check for frequency domain injection. All aproximants supported by XLALSimInspiralImplementedFDApproximants will work.
   * CAVEAT: FD spinning approximants will refer the spin to the lower frequency as given in the xml table. Templates instead will refer it to the lower cutoff of the likelihood integral. This means *different* values of spin will be recovered if one doesn't pay attention! */
	if(XLALSimInspiralImplementedFDApproximants(XLALGetApproximantFromString(injEvent->waveform)))
	{
	 InjectFD(IFOdata, injTable, commandLine);
	 LALInferencePrintDataWithInjection(IFOdata,commandLine);
	 return;
	}
	/* Begin loop over interferometers */
	while(thisData){
		Approximant       approximant;        /* Get approximant value      */
		approximant = XLALGetApproximantFromString(injEvent->waveform);
		if( (int) approximant == XLAL_FAILURE)
			ABORTXLAL(&status);

		InjSampleRate=1.0/thisData->timeData->deltaT;
		if(LALInferenceGetProcParamVal(commandLine,"--injectionsrate")) InjSampleRate=atof(LALInferenceGetProcParamVal(commandLine,"--injectionsrate")->value);
		if(approximant == NumRelNinja2 && InjSampleRate != 16384) {
			fprintf(stderr, "WARNING: NINJA2 injections only work with 16384 Hz sampling rates.  Generating injection in %s at this rate, then downsample to the run's sampling rate.\n", thisData->name);
			InjSampleRate = 16384;
		}

		memset(&det,0,sizeof(det));
		det.site=thisData->detector;
		COMPLEX8FrequencySeries *resp = XLALCreateCOMPLEX8FrequencySeries("response",&thisData->timeData->epoch,
																		  0.0,
																		  thisData->freqData->deltaF,
																		  &strainPerCount,
																		  thisData->freqData->data->length);
		
		for(i=0;i<resp->data->length;i++) {resp->data->data[i] = 1.0;}
		/* Originally created for injecting into DARM-ERR, so transfer function was needed.  
		But since we are injecting into h(t), the transfer function from h(t) to h(t) is 1.*/

		/* We need a long buffer to inject into so that FindChirpInjectSignals() works properly
		 for low mass systems. Use 100 seconds here */
		bufferN = (UINT4) (bufferLength*InjSampleRate);// /thisData->timeData->deltaT);
		memcpy(&bufferStart,&thisData->timeData->epoch,sizeof(LIGOTimeGPS));
		XLALGPSAdd(&bufferStart,(REAL8) thisData->timeData->data->length * thisData->timeData->deltaT);
		XLALGPSAdd(&bufferStart,-bufferLength);
		injectionBuffer=(REAL4TimeSeries *)XLALCreateREAL4TimeSeries(thisData->detector->frDetector.prefix,
																	 &bufferStart, 0.0, 1.0/InjSampleRate,//thisData->timeData->deltaT,
																	 &lalADCCountUnit, bufferN);
		REAL8TimeSeries *inj8Wave=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("injection8",
                                                                           &thisData->timeData->epoch,
                                                                           0.0,
                                                                           thisData->timeData->deltaT,
                                                                           //&lalDimensionlessUnit,
                                                                           &lalStrainUnit,
                                                                           thisData->timeData->data->length);
		if(!inj8Wave) XLAL_ERROR_VOID(XLAL_EFUNC);
		/* This marks the sample in which the real segment starts, within the buffer */
		for(i=0;i<injectionBuffer->data->length;i++) injectionBuffer->data->data[i]=0.0;
		for(i=0;i<inj8Wave->data->length;i++) inj8Wave->data->data[i]=0.0;
		INT4 realStartSample=(INT4)((thisData->timeData->epoch.gpsSeconds - injectionBuffer->epoch.gpsSeconds)/thisData->timeData->deltaT);
		realStartSample+=(INT4)((thisData->timeData->epoch.gpsNanoSeconds - injectionBuffer->epoch.gpsNanoSeconds)*1e-9/thisData->timeData->deltaT);

		/*LALSimulateCoherentGW(&status,injWave,&InjectGW,&det);*/
    //LALFindChirpInjectSignals(&status,injectionBuffer,injEvent,resp);

    if(LALInferenceGetProcParamVal(commandLine,"--lalinspiralinjection")){
      if ( approximant == NumRelNinja2) {
        XLALSimInjectNinjaSignals(injectionBuffer, thisData->name, 1./responseScale, injEvent);
      } else {
        /* Use this custom version for extra sites - not currently maintained */
        // LALInferenceLALFindChirpInjectSignals (&status,injectionBuffer,injEvent,resp,det.site);
	      /* Normal find chirp simulation cannot handle the extra sites */
	      LALFindChirpInjectSignals (&status,injectionBuffer,injEvent,resp);
      }
      printf("Using LALInspiral for injection\n");
      XLALResampleREAL4TimeSeries(injectionBuffer,thisData->timeData->deltaT); //downsample to analysis sampling rate.
      if(status.statusCode) REPORTSTATUS(&status);
      XLALDestroyCOMPLEX8FrequencySeries(resp);
      
      if ( approximant != NumRelNinja2 ) {
        /* Checking the lenght of the injection waveform with respect of thisData->timeData->data->length */
        CoherentGW            waveform;
        PPNParamStruc         ppnParams;
        memset( &waveform, 0, sizeof(CoherentGW) );
        memset( &ppnParams, 0, sizeof(PPNParamStruc) );
        ppnParams.deltaT   = 1.0/InjSampleRate;//thisData->timeData->deltaT;
        ppnParams.lengthIn = 0;
        ppnParams.ppn      = NULL;
        unsigned lengthTest = 0;
        
        LALGenerateInspiral(&status, &waveform, injEvent, &ppnParams ); //Recompute the waveform just to get access to ppnParams.tc and waveform.h->data->length or waveform.phi->data->length
        if(status.statusCode) REPORTSTATUS(&status);
        
        if(waveform.h){
          lengthTest = waveform.h->data->length*(thisData->timeData->deltaT*InjSampleRate);
        }
        if(waveform.phi){
          XLALResampleREAL8TimeSeries(waveform.phi,thisData->timeData->deltaT);
          lengthTest = waveform.phi->data->length;
        }
        
        
        if(lengthTest>thisData->timeData->data->length-(UINT4)ceil((2.0*padding+2.0)/thisData->timeData->deltaT)){
          fprintf(stderr, "WARNING: waveform length = %u is longer than thisData->timeData->data->length = %d minus the window width = %d and the 2.0 seconds after tc (total of %d points available).\n", lengthTest, thisData->timeData->data->length, (INT4)ceil((2.0*padding)/thisData->timeData->deltaT) , thisData->timeData->data->length-(INT4)ceil((2.0*padding+2.0)/thisData->timeData->deltaT));
          fprintf(stderr, "The waveform injected is %f seconds long. Consider increasing the %f seconds segment length (--seglen) to be greater than %f. (in %s, line %d)\n",ppnParams.tc , thisData->timeData->data->length * thisData->timeData->deltaT, ppnParams.tc + 2.0*padding + 2.0, __FILE__, __LINE__);
        }
        if(ppnParams.tc>bufferLength){
          fprintf(stderr, "ERROR: The waveform injected is %f seconds long and the buffer for FindChirpInjectSignal is %f seconds long. The end of the waveform will be cut ! (in %s, line %d)\n",ppnParams.tc , bufferLength, __FILE__, __LINE__);
          exit(1);
        }
      }
      
      /* Now we cut the injection buffer down to match the time domain wave size */
      injectionBuffer=(REAL4TimeSeries *)XLALCutREAL4TimeSeries(injectionBuffer,realStartSample,thisData->timeData->data->length);
      if (!injectionBuffer) XLAL_ERROR_VOID(XLAL_EFUNC);
      if(status.statusCode) REPORTSTATUS(&status);
      /*		for(j=0;j<injWave->data->length;j++) printf("%f\n",injWave->data->data[j]);*/
      for(i=0;i<injectionBuffer->data->length;i++) inj8Wave->data->data[i]=(REAL8)injectionBuffer->data->data[i];
    }else{
      printf("Using LALSimulation for injection\n");
      REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform */
      REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform */
      REAL8TimeSeries       *signalvecREAL8=NULL;
      LALPNOrder        order;              /* Order of the model             */
      INT4              amporder=0;         /* Amplitude order of the model   */
      
      order = XLALGetOrderFromString(injEvent->waveform);
      if ( (int) order == XLAL_FAILURE)
        ABORTXLAL(&status);
      amporder = injEvent->amp_order;
      //if(amporder<0) amporder=0;
      /* FIXME - tidal lambda's and interactionFlag are just set to command line values here.
       * They should be added to injEvent and set to appropriate values
       */
      REAL8 lambda1 = 0.;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")) {
        lambda1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")->value);
        fprintf(stdout,"Injection lambda1 set to %f\n",lambda1);
      }
      REAL8 lambda2 = 0.;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")) {
        lambda2= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")->value);
        fprintf(stdout,"Injection lambda2 set to %f\n",lambda2);
      }
      REAL8 lambdaT = 0.;
      REAL8 dLambdaT = 0.;
      REAL8 m1=injEvent->mass1;
      REAL8 m2=injEvent->mass2;
      REAL8 Mt=m1+m2;
      REAL8 eta=m1*m2/(Mt*Mt);
      if(LALInferenceGetProcParamVal(commandLine,"--inj-lambdaT")&&LALInferenceGetProcParamVal(commandLine,"--inj-dLambdaT")) {
        lambdaT= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambdaT")->value);
        dLambdaT= atof(LALInferenceGetProcParamVal(commandLine,"--inj-dLambdaT")->value);
        LALInferenceLambdaTsEta2Lambdas(lambdaT,dLambdaT,eta,&lambda1,&lambda2);
        fprintf(stdout,"Injection lambdaT set to %f\n",lambdaT);
        fprintf(stdout,"Injection dLambdaT set to %f\n",dLambdaT);
        fprintf(stdout,"lambda1 set to %f\n",lambda1);
        fprintf(stdout,"lambda2 set to %f\n",lambda2);
      }

      REAL8 fref = 100.;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-fref")) {
        fref = atoi(LALInferenceGetProcParamVal(commandLine,"--inj-fref")->value);
      }

      LALSimInspiralWaveformFlags *waveFlags = XLALSimInspiralCreateWaveformFlags();
      LALSimInspiralSpinOrder spinO = -1;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-spinOrder")) {
        spinO = atoi(LALInferenceGetProcParamVal(commandLine,"--inj-spinOrder")->value);
        XLALSimInspiralSetSpinOrder(waveFlags, spinO);
      }
      LALSimInspiralTidalOrder tideO = -1;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-tidalOrder")) {
        tideO = atoi(LALInferenceGetProcParamVal(commandLine,"--inj-tidalOrder")->value);
        XLALSimInspiralSetTidalOrder(waveFlags, tideO);
      }
      LALSimInspiralTestGRParam *nonGRparams = NULL;
      /* Print a line with information about approximant, amporder, phaseorder, tide order and spin order */
      fprintf(stdout,"Injection will run using Approximant %i (%s), phase order %i, amp order %i, spin order %i, tidal order %i, in the time domain with a reference frequency of %f.\n",approximant,XLALGetStringFromApproximant(approximant),order,amporder,(int) spinO, (int) tideO, (float) fref);

      /* ChooseWaveform starts the (2,2) mode of the waveform at the given minimum frequency.  We want the highest order contribution to start at the f_lower of the injection file */
      REAL8 f_min = fLow2fStart(injEvent->f_lower, amporder, approximant);
      printf("Injecting with f_min = %f.\n", f_min);

      XLALSimInspiralChooseTDWaveform(&hplus, &hcross, injEvent->coa_phase, 1.0/InjSampleRate,
                                      injEvent->mass1*LAL_MSUN_SI, injEvent->mass2*LAL_MSUN_SI, injEvent->spin1x,
                                      injEvent->spin1y, injEvent->spin1z, injEvent->spin2x, injEvent->spin2y,
                                      injEvent->spin2z, f_min, fref, injEvent->distance*LAL_PC_SI * 1.0e6,
                                      injEvent->inclination, lambda1, lambda2, waveFlags,
                                      nonGRparams, amporder, order, approximant);
      if(!hplus || !hcross) {
        fprintf(stderr,"Error: XLALSimInspiralChooseWaveform() failed to produce waveform.\n");
        exit(-1);
      }
      XLALSimInspiralDestroyWaveformFlags(waveFlags);
      XLALSimInspiralDestroyTestGRParam(nonGRparams);
      XLALResampleREAL8TimeSeries(hplus,thisData->timeData->deltaT);
      XLALResampleREAL8TimeSeries(hcross,thisData->timeData->deltaT);
      /* XLALSimInspiralChooseTDWaveform always ends the waveform at t=0 */
      /* So we can adjust the epoch so that the end time is as desired */
      XLALGPSAddGPS(&(hplus->epoch), &(injEvent->geocent_end_time));
      XLALGPSAddGPS(&(hcross->epoch), &(injEvent->geocent_end_time));
      //XLALGPSAdd(&(hplus->epoch), -(REAL8)hplus->data->length*hplus->deltaT);
      //XLALGPSAdd(&(hcross->epoch), -(REAL8)hcross->data->length*hplus->deltaT);
      
      signalvecREAL8=XLALSimDetectorStrainREAL8TimeSeries(hplus, hcross, injEvent->longitude, injEvent->latitude, injEvent->polarization, det.site);
      if (!signalvecREAL8) XLAL_ERROR_VOID(XLAL_EFUNC);
      
      for(i=0;i<signalvecREAL8->data->length;i++){
        if(isnan(signalvecREAL8->data->data[i])) signalvecREAL8->data->data[i]=0.0;
      }

      if(signalvecREAL8->data->length > thisData->timeData->data->length-(UINT4)ceil((2.0*padding+2.0)/thisData->timeData->deltaT)){
        fprintf(stderr, "WARNING: waveform length = %u is longer than thisData->timeData->data->length = %d minus the window width = %d and the 2.0 seconds after tc (total of %d points available).\n", signalvecREAL8->data->length, thisData->timeData->data->length, (INT4)ceil((2.0*padding)/thisData->timeData->deltaT) , thisData->timeData->data->length-(INT4)ceil((2.0*padding+2.0)/thisData->timeData->deltaT));
        fprintf(stderr, "The waveform injected is %f seconds long. Consider increasing the %f seconds segment length (--seglen) to be greater than %f. (in %s, line %d)\n",signalvecREAL8->data->length * thisData->timeData->deltaT , thisData->timeData->data->length * thisData->timeData->deltaT, signalvecREAL8->data->length * thisData->timeData->deltaT + 2.0*padding + 2.0, __FILE__, __LINE__);
      }
      
      XLALSimAddInjectionREAL8TimeSeries(inj8Wave, signalvecREAL8, NULL);
      
      if ( hplus ) XLALDestroyREAL8TimeSeries(hplus);
      if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
      
    }
    XLALDestroyREAL4TimeSeries(injectionBuffer);
    injF=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("injF",
										&thisData->timeData->epoch,
										0.0,
										thisData->freqData->deltaF,
										&lalDimensionlessUnit,
										thisData->freqData->data->length);
    if(!injF) {
      XLALPrintError("Unable to allocate memory for injection buffer\n");
      XLAL_ERROR_VOID(XLAL_EFUNC);
    }
    /* Window the data */
    REAL4 WinNorm = sqrt(thisData->window->sumofsquares/thisData->window->data->length);
        for(j=0;j<inj8Wave->data->length;j++) inj8Wave->data->data[j]*=thisData->window->data->data[j]; /* /WinNorm; */ /* Window normalisation applied only in freq domain */
    XLALREAL8TimeFreqFFT(injF,inj8Wave,thisData->timeToFreqFFTPlan);
    /*for(j=0;j<injF->data->length;j++) printf("%lf\n",injF->data->data[j].re);*/
    if(thisData->oneSidedNoisePowerSpectrum){
      for(SNR=0.0,j=thisData->fLow/injF->deltaF;j<thisData->fHigh/injF->deltaF;j++){
        SNR += pow(creal(injF->data->data[j]), 2.0)/thisData->oneSidedNoisePowerSpectrum->data->data[j];
        SNR += pow(cimag(injF->data->data[j]), 2.0)/thisData->oneSidedNoisePowerSpectrum->data->data[j];
      }
      SNR*=4.0*injF->deltaF;
    }
    thisData->SNR=sqrt(SNR);
    NetworkSNR+=SNR;

    /* Actually inject the waveform */
    for(j=0;j<inj8Wave->data->length;j++) thisData->timeData->data->data[j]+=inj8Wave->data->data[j];
      fprintf(stdout,"Injected SNR in detector %s = %g\n",thisData->name,thisData->SNR);
      char filename[256];
      sprintf(filename,"%s_timeInjection.dat",thisData->name);
      FILE* file=fopen(filename, "w");
      for(j=0;j<inj8Wave->data->length;j++){   
	fprintf(file, "%.6f\t%lg\n", XLALGPSGetREAL8(&thisData->timeData->epoch) + thisData->timeData->deltaT*j, inj8Wave->data->data[j]);
      }
      fclose(file);
      sprintf(filename,"%s_freqInjection.dat",thisData->name);
      file=fopen(filename, "w");
      for(j=0;j<injF->data->length;j++){   
	thisData->freqData->data->data[j] += injF->data->data[j] / WinNorm;
	fprintf(file, "%lg %lg \t %lg\n", thisData->freqData->deltaF*j, creal(injF->data->data[j]), cimag(injF->data->data[j]));
      }
      fclose(file);
    
      XLALDestroyREAL8TimeSeries(inj8Wave);
      XLALDestroyCOMPLEX16FrequencySeries(injF);
      thisData=thisData->next;
    }
     if (!(SNRpath==NULL)){ /* If the user provided a path with --snrpath store a file with injected SNRs */
      PrintSNRsToFile(IFOdata , injTable);
    }

    NetworkSNR=sqrt(NetworkSNR);
    fprintf(stdout,"Network SNR of event %d = %g\n",event,NetworkSNR);

    /* Output waveform raw h-plus mode */
    if( (ppt=LALInferenceGetProcParamVal(commandLine,"--rawwaveform")) )
    {
        rawWaveform=fopen(ppt->value,"w");
        bufferN = (UINT4) (bufferLength/IFOdata->timeData->deltaT);
        memcpy(&bufferStart,&IFOdata->timeData->epoch,sizeof(LIGOTimeGPS));
        XLALGPSAdd(&bufferStart,(REAL8) IFOdata->timeData->data->length * IFOdata->timeData->deltaT);
        XLALGPSAdd(&bufferStart,-bufferLength);
        COMPLEX8FrequencySeries *resp = XLALCreateCOMPLEX8FrequencySeries("response",&IFOdata->timeData->epoch,0.0,IFOdata->freqData->deltaF,&strainPerCount,IFOdata->freqData->data->length);
        if(!resp) XLAL_ERROR_VOID(XLAL_EFUNC);
        injectionBuffer=(REAL4TimeSeries *)XLALCreateREAL4TimeSeries("None",&bufferStart, 0.0, IFOdata->timeData->deltaT,&lalADCCountUnit, bufferN);
        if(!injectionBuffer) XLAL_ERROR_VOID(XLAL_EFUNC);
        /* This marks the sample in which the real segment starts, within the buffer */
        INT4 realStartSample=(INT4)((IFOdata->timeData->epoch.gpsSeconds - injectionBuffer->epoch.gpsSeconds)/IFOdata->timeData->deltaT);
        realStartSample+=(INT4)((IFOdata->timeData->epoch.gpsNanoSeconds - injectionBuffer->epoch.gpsNanoSeconds)*1e-9/IFOdata->timeData->deltaT);
        LALFindChirpInjectSignals(&status,injectionBuffer,injEvent,resp);
        if(status.statusCode) REPORTSTATUS(&status);
        XLALDestroyCOMPLEX8FrequencySeries(resp);
        injectionBuffer=(REAL4TimeSeries *)XLALCutREAL4TimeSeries(injectionBuffer,realStartSample,IFOdata->timeData->data->length);
        for(j=0;j<injectionBuffer->data->length;j++) fprintf(rawWaveform,"%.6f\t%g\n", XLALGPSGetREAL8(&IFOdata->timeData->epoch) + IFOdata->timeData->deltaT*j, injectionBuffer->data->data[j]);
        fclose(rawWaveform);
        XLALDestroyREAL4TimeSeries(injectionBuffer);
    }
  
    LALInferencePrintDataWithInjection(IFOdata,commandLine);
  
    return;
}

//temporary? replacement function for FindChirpInjectSignals in order to accept any detector.site and not only the ones in lalCachedDetectors.
void
LALInferenceLALFindChirpInjectSignals (
    LALStatus                  *status,
    REAL4TimeSeries            *chan,
    SimInspiralTable           *events,
    COMPLEX8FrequencySeries    *resp,
    LALDetector                *LALInference_detector
    )

{
  UINT4                 k;
  DetectorResponse      detector;
  SimInspiralTable     *thisEvent = NULL;
  PPNParamStruc         ppnParams;
  CoherentGW            waveform;
  INT8                  waveformStartTime;
  REAL4TimeSeries       signalvec;
  COMPLEX8Vector       *unity = NULL;
  CHAR                  warnMsg[512];
  CHAR                  ifo[LIGOMETA_IFO_MAX];
  REAL8                 timeDelay;
  UINT4                  i; 
  REAL8TimeSeries       *hplus=NULL;
  REAL8TimeSeries       *hcross=NULL;
  REAL8TimeSeries       *signalvecREAL8=NULL;
 // REAL4TimeSeries       *signalvecREAL4=NULL;
  
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( chan, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( chan->data, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( chan->data->data, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  ASSERT( events, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  ASSERT( resp, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( resp->data, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( resp->data->data, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );


  /*
   *
   * set up structures and parameters needed
   *
   */


  /* fixed waveform injection parameters */
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  ppnParams.deltaT   = chan->deltaT;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;


  /*
   *
   * compute the transfer function from the given response function
   *
   */


  /* allocate memory and copy the parameters describing the freq series */
  memset( &detector, 0, sizeof( DetectorResponse ) );
  detector.transfer = (COMPLEX8FrequencySeries *)
    LALCalloc( 1, sizeof(COMPLEX8FrequencySeries) );
  if ( ! detector.transfer )
  {
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }
  memcpy( &(detector.transfer->epoch), &(resp->epoch),
      sizeof(LIGOTimeGPS) );
  detector.transfer->f0 = resp->f0;
  detector.transfer->deltaF = resp->deltaF;

  detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );
  /* set the detector site */
  
  detector.site = LALInference_detector;
  strcpy(ifo, LALInference_detector->frDetector.prefix);
  printf("computing waveform for %s\n",LALInference_detector->frDetector.name);

  /* set up units for the transfer function */
  if (XLALUnitDivide( &(detector.transfer->sampleUnits),
                      &lalADCCountUnit, &lalStrainUnit ) == NULL) {
    ABORTXLAL(status);
  }

  /* invert the response function to get the transfer function */
  LALCCreateVector( status->statusPtr, &( detector.transfer->data ),
      resp->data->length );
  CHECKSTATUSPTR( status );

  LALCCreateVector( status->statusPtr, &unity, resp->data->length );
  CHECKSTATUSPTR( status );
  for ( k = 0; k < resp->data->length; ++k )
  {
    unity->data[k] = 1.0;
  }

  LALCCVectorDivide( status->statusPtr, detector.transfer->data, unity,
      resp->data );
  CHECKSTATUSPTR( status );

  LALCDestroyVector( status->statusPtr, &unity );
  CHECKSTATUSPTR( status );


  /*
   *
   * loop over the signals and inject them into the time series
   *
   */


  for ( thisEvent = events; thisEvent; thisEvent = thisEvent->next )
  {
    /*
     *
     * generate waveform and inject it into the data
     *
     */


    /* clear the waveform structure */
    memset( &waveform, 0, sizeof(CoherentGW) );
    
    LALGenerateInspiral(status->statusPtr, &waveform, thisEvent, &ppnParams );
    CHECKSTATUSPTR( status );
    
    LALInfo( status, ppnParams.termDescription );

    if ( strstr( thisEvent->waveform, "KludgeIMR") ||
         strstr( thisEvent->waveform, "KludgeRingOnly") )
     {
       CoherentGW *wfm;
       SimRingdownTable *ringEvent;
       int injectSignalType = LALRINGDOWN_IMR_INJECT;


       ringEvent = (SimRingdownTable *)
         LALCalloc( 1, sizeof(SimRingdownTable) );
       wfm = XLALGenerateInspRing( &waveform, thisEvent, ringEvent,
           injectSignalType);
       LALFree(ringEvent);

       if ( !wfm )
       {
         LALInfo( status, "Unable to generate merger/ringdown, "
             "injecting inspiral only");
         ABORT( status, FINDCHIRPH_EIMRW, FINDCHIRPH_MSGEIMRW );
       }
       waveform = *wfm;
     }


    if ( thisEvent->geocent_end_time.gpsSeconds )
    {
      /* get the gps start time of the signal to inject */
      waveformStartTime = XLALGPSToINT8NS( &(thisEvent->geocent_end_time) );
      waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );
    }
    else
    {
      LALInfo( status, "Waveform start time is zero: injecting waveform "
          "into center of data segment" );

      /* center the waveform in the data segment */
      waveformStartTime = XLALGPSToINT8NS( &(chan->epoch) );

      waveformStartTime += (INT8) ( 1000000000.0 *
          ((REAL8) (chan->data->length - ppnParams.length) / 2.0) * chan->deltaT
          );
    }

    snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
        "Injected waveform timing:\n"
        "thisEvent->geocent_end_time.gpsSeconds = %d\n"
        "thisEvent->geocent_end_time.gpsNanoSeconds = %d\n"
        "ppnParams.tc = %e\n"
        "waveformStartTime = %" LAL_INT8_FORMAT "\n",
        thisEvent->geocent_end_time.gpsSeconds,
        thisEvent->geocent_end_time.gpsNanoSeconds,
        ppnParams.tc,
        waveformStartTime );
    LALInfo( status, warnMsg );

      /* clear the signal structure */
      memset( &signalvec, 0, sizeof(REAL4TimeSeries) );

      /* set the start time of the signal vector to the appropriate start time of the injection */
      if ( detector.site )
      {
        timeDelay = XLALTimeDelayFromEarthCenter( detector.site->location, thisEvent->longitude,
          thisEvent->latitude, &(thisEvent->geocent_end_time) );
        if ( XLAL_IS_REAL8_FAIL_NAN( timeDelay ) )
        {
          ABORTXLAL( status );
        }
      }
      else
      {
        timeDelay = 0.0;
      }
      /* Give a little more breathing space to aid band-passing */
      XLALGPSSetREAL8( &(signalvec.epoch), (waveformStartTime * 1.0e-9) - 0.25 + timeDelay );
      /* set the parameters for the signal time series */
      signalvec.deltaT = chan->deltaT;
      if ( ( signalvec.f0 = chan->f0 ) != 0 )
      {
        ABORT( status, FINDCHIRPH_EHETR, FINDCHIRPH_MSGEHETR );
      }
      signalvec.sampleUnits = lalADCCountUnit;
      
      if(waveform.h == NULL){
      /* set the start times for injection */
      XLALINT8NSToGPS( &(waveform.a->epoch), waveformStartTime );
      /* put a rug on a polished floor? */
      waveform.f->epoch = waveform.a->epoch;
      waveform.phi->epoch = waveform.a->epoch;
      /* you might as well set a man trap */
      if ( waveform.shift )
      {
        waveform.shift->epoch = waveform.a->epoch;
      }
      /* and to think he'd just come from the hospital */
      }else{
        /* set the start times for injection */
        XLALINT8NSToGPS( &(waveform.h->epoch), waveformStartTime );  
      }
      /* simulate the detectors response to the inspiral */
      LALSCreateVector( status->statusPtr, &(signalvec.data), chan->data->length );
      CHECKSTATUSPTR( status );

      if(waveform.h == NULL){ //LALSimulateCoherentGW only for waveform generators filling CoherentGW.a and CoherentGW.phi
        LALSimulateCoherentGW( status->statusPtr, &signalvec, &waveform, &detector );
      }else{
      hplus=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("hplus",
                                                                &(waveform.h->epoch),
                                                                0.0,
                                                                waveform.h->deltaT,
                                                                &lalDimensionlessUnit,
                                                                waveform.h->data->length);

      hcross=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("hcross",
                                                                  &(waveform.h->epoch),
                                                                  0.0,
                                                                  waveform.h->deltaT,
                                                                  &lalDimensionlessUnit,
                                                                  waveform.h->data->length);
      for( i = 0; i < waveform.h->data->length; i++)
      {
        hplus->data->data[i] = waveform.h->data->data[2*i];
        hcross->data->data[i] = waveform.h->data->data[(2*i)+1];
      }

      signalvecREAL8=XLALSimDetectorStrainREAL8TimeSeries(hplus, 
                                                          hcross,
                                                          thisEvent->longitude,
                                                          thisEvent->latitude,
                                                          thisEvent->polarization,
                                                          LALInference_detector);
        
      INT8 offset = ( signalvecREAL8->epoch.gpsSeconds - signalvec.epoch.gpsSeconds ) / signalvec.deltaT;
      offset += ( signalvecREAL8->epoch.gpsNanoSeconds - signalvec.epoch.gpsNanoSeconds ) * 1.0e-9 / signalvec.deltaT;

      
      int Nnans=0;
      for (i=0; i<signalvec.data->length; i++){
        if(i<offset || i>=signalvecREAL8->data->length+offset || isnan(signalvecREAL8->data->data[i-offset])) signalvec.data->data[i]=0.0; //The isnan() condition should not be necessary. To be investigated.
	else signalvec.data->data[i]=(REAL4) signalvecREAL8->data->data[i-offset];
	if((i>=offset)&&(i<signalvecREAL8->data->length+offset) && isnan(signalvecREAL8->data->data[i-offset])) Nnans++;
      }
      if(Nnans>0) fprintf(stderr,"Trimmed %i NaNs from the injection waveform!\n",Nnans);
      }
      CHECKSTATUSPTR( status );
      
      /* Taper the signal */
      {

          if ( ! strcmp( "TAPER_START", thisEvent->taper ) )
          {
              XLALSimInspiralREAL4WaveTaper( signalvec.data, LAL_SIM_INSPIRAL_TAPER_START );
          }
          else if (  ! strcmp( "TAPER_END", thisEvent->taper ) )
          {
              XLALSimInspiralREAL4WaveTaper( signalvec.data, LAL_SIM_INSPIRAL_TAPER_END );
          }
          else if (  ! strcmp( "TAPER_STARTEND", thisEvent->taper ) )
          {
              XLALSimInspiralREAL4WaveTaper( signalvec.data, LAL_SIM_INSPIRAL_TAPER_STARTEND );
          }
          else if ( strcmp( "TAPER_NONE", thisEvent->taper ) )
          {
              XLALPrintError( "Invalid injection tapering option specified: %s\n",
                 thisEvent->taper );
              ABORT( status, LAL_BADPARM_ERR, LAL_BADPARM_MSG );
          }
      }
      
      /* Band pass the signal */
      if ( thisEvent->bandpass )
      {
          UINT4 safeToBandPass = 0;
          UINT4 start=0, end=0;
          REAL4Vector *bandpassVec = NULL;

          safeToBandPass = FindTimeSeriesStartAndEnd (
                  signalvec.data, &start, &end );

          if ( safeToBandPass )
          {
              /* Check if we can grab some padding at the extremeties.
               * This will make the bandpassing better
               */

              if (((INT4)start - (int)(0.25/chan->deltaT)) > 0 )
                    start -= (int)(0.25/chan->deltaT);
              else
                    start = 0;

              if ((end + (int)(0.25/chan->deltaT)) < signalvec.data->length )
                    end += (int)(0.25/chan->deltaT);
              else
                    end = signalvec.data->length - 1;

              bandpassVec = (REAL4Vector *)
                      LALCalloc(1, sizeof(REAL4Vector) );

              bandpassVec->length = (end - start + 1);
              bandpassVec->data = signalvec.data->data + start;

              if ( XLALBandPassInspiralTemplate( bandpassVec,
                          1.1*thisEvent->f_lower,
                          1.05*thisEvent->f_final,
                          1./chan->deltaT) != XLAL_SUCCESS )
              {
                  LALError( status, "Failed to Bandpass signal" );
                  ABORT (status, LALINSPIRALH_EBPERR, LALINSPIRALH_MSGEBPERR);
              };

              LALFree( bandpassVec );
          }
      }
      /* inject the signal into the data channel */
      LALSSInjectTimeSeries( status->statusPtr, chan, &signalvec );

      CHECKSTATUSPTR( status );


    if ( waveform.shift )
    {
      LALSDestroyVector( status->statusPtr, &(waveform.shift->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.shift );
    }

    if( waveform.h )
    {
      LALSDestroyVectorSequence( status->statusPtr, &(waveform.h->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.h );
    }
    if( waveform.a )
    {
      LALSDestroyVectorSequence( status->statusPtr, &(waveform.a->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.a );
      /*
       * destroy the signal only if waveform.h is NULL as otherwise it won't
       * be created
       * */
      if ( waveform.h == NULL )
      {
	LALSDestroyVector( status->statusPtr, &(signalvec.data) );
        CHECKSTATUSPTR( status );
      }
    }
    if( waveform.f )
    {
      LALSDestroyVector( status->statusPtr, &(waveform.f->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.f );
    }
    if( waveform.phi )
    {
      LALDDestroyVector( status->statusPtr, &(waveform.phi->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.phi );
    }
  }

  
  if(hplus) XLALDestroyREAL8TimeSeries(hplus);
  if(hcross) XLALDestroyREAL8TimeSeries(hcross);
  if(signalvecREAL8) XLALDestroyREAL8TimeSeries(signalvecREAL8);
  
  LALCDestroyVector( status->statusPtr, &( detector.transfer->data ) );
  CHECKSTATUSPTR( status );

//  if ( detector.site ) LALFree( detector.site );
  LALFree( detector.transfer );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

static int FindTimeSeriesStartAndEnd (
                                      REAL4Vector *signalvec,
                                      UINT4 *start,
                                      UINT4 *end
                                      )
{
  UINT4 i; /* mid, n; indices */
  UINT4 flag, safe = 1;
  UINT4 length;
  
#ifndef LAL_NDEBUG
  if ( !signalvec )
    XLAL_ERROR( XLAL_EFAULT );
  
  if ( !signalvec->data )
    XLAL_ERROR( XLAL_EFAULT );
#endif
  
  length = signalvec->length;
  
  /* Search for start and end of signal */
  flag = 0;
  i = 0;
  while(flag == 0 && i < length )
  {
    if( signalvec->data[i] != 0.)
    {
      *start = i;
      flag = 1;
    }
    i++;
  }
  if ( flag == 0 )
  {
    return flag;
  }
  
  flag = 0;
  i = length - 1;
  while(flag == 0)
  {
    if( signalvec->data[i] != 0.)
    {
      *end = i;
      flag = 1;
    }
    i--;
  }
  
  /* Check we have more than 2 data points */
  if(((*end) - (*start)) <= 1)
  {
    XLALPrintWarning( "Data less than 3 points in this signal!\n" );
    safe = 0;
  }
  
  return safe;
  
}

void InjectFD(LALInferenceIFOData *IFOdata, SimInspiralTable *inj_table, ProcessParamsTable *commandLine)
///*-------------- Inject in Frequency domain -----------------*/
{
  /* Inject a gravitational wave into the data in the frequency domain */
  LALStatus status;
  memset(&status,0,sizeof(LALStatus));
  INT4 errnum;

  Approximant approximant = XLALGetApproximantFromString(inj_table->waveform);
  if( (int) approximant == XLAL_FAILURE)
      ABORTXLAL(&status);

  LALPNOrder phase_order = XLALGetOrderFromString(inj_table->waveform);
  if ( (int) phase_order == XLAL_FAILURE)
      ABORTXLAL(&status);

  LALPNOrder amp_order = (LALPNOrder) inj_table->amp_order;

  enforce_m1_larger_m2(inj_table);

  REAL8 injtime=0.0;
  injtime=(REAL8) inj_table->geocent_end_time.gpsSeconds + (REAL8) inj_table->geocent_end_time.gpsNanoSeconds*1.0e-9;

  REAL8 lambda1 = 0.;
  if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")) {
    lambda1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")->value);
    fprintf(stdout,"Injection lambda1 set to %f\n",lambda1);
  }

  REAL8 lambda2 = 0.;
  if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")) {
    lambda2= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")->value);
    fprintf(stdout,"Injection lambda2 set to %f\n",lambda2);
  }

  REAL8 lambdaT = 0.;
  REAL8 dLambdaT = 0.;

  if(LALInferenceGetProcParamVal(commandLine,"--inj-lambdaT")&&LALInferenceGetProcParamVal(commandLine,"--inj-dLambdaT")) {
    lambdaT= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambdaT")->value);
    dLambdaT= atof(LALInferenceGetProcParamVal(commandLine,"--inj-dLambdaT")->value);
    LALInferenceLambdaTsEta2Lambdas(lambdaT, dLambdaT, inj_table->eta, &lambda1, &lambda2);
    fprintf(stdout,"Injection lambdaT set to %f\n",lambdaT);
    fprintf(stdout,"Injection dLambdaT set to %f\n",dLambdaT);
    fprintf(stdout,"lambda1 set to %f\n",lambda1);
    fprintf(stdout,"lambda2 set to %f\n",lambda2);
  }

  /* Set up wave flags */
  LALSimInspiralWaveformFlags *waveFlags = XLALSimInspiralCreateWaveformFlags();

  LALSimInspiralSpinOrder spinO = LAL_SIM_INSPIRAL_SPIN_ORDER_ALL;
  if(LALInferenceGetProcParamVal(commandLine, "--inj-spinOrder")) {
    spinO = atoi(LALInferenceGetProcParamVal(commandLine, "--inj-spinOrder")->value);
    XLALSimInspiralSetSpinOrder(waveFlags, spinO);
  }

  LALSimInspiralTidalOrder tideO = LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL;
  if(LALInferenceGetProcParamVal(commandLine, "--inj-tidalOrder")) {
    tideO = atoi(LALInferenceGetProcParamVal(commandLine, "--inj-tidalOrder")->value);
    XLALSimInspiralSetTidalOrder(waveFlags, tideO);
  }

  REAL8 deltaT = IFOdata->timeData->deltaT;
  REAL8 deltaF = IFOdata->freqData->deltaF;

  REAL8 f_min = fLow2fStart(inj_table->f_lower, amp_order, approximant);
  REAL8 f_max = 0.0;

  REAL8 fref = 100.;
  if(LALInferenceGetProcParamVal(commandLine,"--inj-fref")) {
    fref = atoi(LALInferenceGetProcParamVal(commandLine,"--inj-fref")->value);
  }

  LALSimInspiralTestGRParam *nonGRparams = NULL;

 /* Print a line with information about approximant, amp_order, phaseorder, tide order and spin order */
  fprintf(stdout,"\n\n---\t\t ---\n");
 fprintf(stdout,"Injection will run using Approximant %i (%s), phase order %i, amp order %i, spin order %i, tidal order %i, in the frequency domain.\n",approximant,XLALGetStringFromApproximant(approximant),phase_order,amp_order,(int) spinO,(int) tideO);
   fprintf(stdout,"---\t\t ---\n\n");

  COMPLEX16FrequencySeries *hptilde=NULL, *hctilde=NULL;

  XLALSimInspiralChooseFDWaveform(&hptilde, &hctilde, inj_table->coa_phase, deltaF,
                                  inj_table->mass1*LAL_MSUN_SI, inj_table->mass2*LAL_MSUN_SI, inj_table->spin1x,
                                  inj_table->spin1y, inj_table->spin1z, inj_table->spin2x, inj_table->spin2y,
                                  inj_table->spin2z, f_min, f_max, fref, inj_table->distance*LAL_PC_SI * 1.0e6,
                                  inj_table->inclination, lambda1, lambda2, waveFlags,
                                  nonGRparams, amp_order, phase_order, approximant);

  /* Fail if injection waveform generation was not successful */
  errnum = *XLALGetErrnoPtr();
  if (errnum != XLAL_SUCCESS) {
    XLALPrintError(" ERROR in InjectFD(): error encountered when injecting waveform. errnum=%d\n",errnum);
    exit(1);
  }

  XLALSimInspiralDestroyWaveformFlags(waveFlags);
  XLALSimInspiralDestroyTestGRParam(nonGRparams);

  LALInferenceIFOData *dataPtr;
  REAL8 Fplus, Fcross;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  LIGOTimeGPS GPSlal;
  REAL8 gmst;
  REAL8 chisquared;
  REAL8 timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  REAL8 timeshift;  /* time shift (not necessarily same as above)                   */
  REAL8 twopit, f, re, im, dre, dim, newRe, newIm;
  INT4 i, lower, upper;

  REAL8 temp=0.0;
  REAL8 NetSNR=0.0;

  /* figure out GMST: */
  XLALGPSSetREAL8(&GPSlal, injtime);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);

  /* loop over data (different interferometers): */
  dataPtr = IFOdata;

  while (dataPtr != NULL) {
    /*-- WF to inject is now in hptilde and hctilde. --*/
    /* determine beam pattern response (Fplus and Fcross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                                (const REAL4(*)[3])dataPtr->detector->response,
                                inj_table->longitude, inj_table->latitude,
                                inj_table->polarization, gmst);

    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                                inj_table->longitude, inj_table->latitude,
                                                &GPSlal);

    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */
    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    REAL8 instant = dataPtr->timeData->epoch.gpsSeconds + 1e-9*dataPtr->timeData->epoch.gpsNanoSeconds;

    timeshift = (injtime - instant) + timedelay;
    twopit    = LAL_TWOPI * (timeshift);

    dataPtr->fPlus = Fplus;
    dataPtr->fCross = Fcross;
    dataPtr->timeshift = timeshift;

    char InjFileName[50];
    sprintf(InjFileName,"injection_%s.dat",dataPtr->name);
    FILE *outInj=fopen(InjFileName,"w");

     /* determine frequency range & loop over frequency bins: */
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
    chisquared = 0.0;

    re = cos(twopit * deltaF * lower);
    im = -sin(twopit * deltaF * lower);
    for (i=lower; i<=upper; ++i){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      if (i < hptilde->data->length) {
          plainTemplateReal = Fplus * creal(hptilde->data->data[i])
                              +  Fcross * creal(hctilde->data->data[i]);
          plainTemplateImag = Fplus * cimag(hptilde->data->data[i])
                              +  Fcross * cimag(hctilde->data->data[i]);
      } else {
          plainTemplateReal = 0.0;
          plainTemplateImag = 0.0;
      }

      /* do time-shifting...             */
      /* (also un-do 1/deltaT scaling): */
      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      templateReal = (plainTemplateReal*re - plainTemplateImag*im);
      templateImag = (plainTemplateReal*im + plainTemplateImag*re);

      /* Incremental values, using cos(theta) - 1 = -2*sin(theta/2)^2 */
      dim = -sin(twopit*deltaF);
      dre = -2.0*sin(0.5*twopit*deltaF)*sin(0.5*twopit*deltaF);
      newRe = re + re*dre - im * dim;
      newIm = im + re*dim + im*dre;
      re = newRe;
      im = newIm;

      fprintf(outInj,"%lf %e %e %e\n",i*deltaF ,templateReal,templateImag,1.0/dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      dataPtr->freqData->data->data[i] += crect( templateReal, templateImag );

      temp = ((2.0/( deltaT*(double) dataPtr->timeData->data->length) * (templateReal*templateReal+templateImag*templateImag)) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      chisquared  += temp;
    }
    printf("injected SNR %.1f in IFO %s\n",sqrt(2.0*chisquared),dataPtr->name);
    NetSNR+=2.0*chisquared;
    dataPtr->SNR=sqrt(2.0*chisquared);
    dataPtr = dataPtr->next;

    fclose(outInj);
  }
  printf("injected Network SNR %.1f \n",sqrt(NetSNR));

  if (!(SNRpath==NULL)){ /* If the user provided a path with --snrpath store a file with injected SNRs */
    PrintSNRsToFile(IFOdata , inj_table);
  }

  XLALDestroyCOMPLEX16FrequencySeries(hctilde);
  XLALDestroyCOMPLEX16FrequencySeries(hptilde);
}


static void PrintSNRsToFile(LALInferenceIFOData *IFOdata , SimInspiralTable *inj_table){
  char SnrName[200];
  char ListOfIFOs[10]="";
  REAL8 NetSNR=0.0;
  LALInferenceIFOData *thisData=IFOdata;
  int nIFO=0;

  while(thisData){
       sprintf(ListOfIFOs,"%s%s",ListOfIFOs,thisData->name);
       thisData=thisData->next;
  nIFO++;
      }

  (void) ListOfIFOs;
  (void) inj_table;
  sprintf(SnrName,"%s/snr_IMR.dat",SNRpath);
  FILE * snrout = fopen(SnrName,"a");
  if(!snrout){
    fprintf(stderr,"Unable to open the path %s for writing SNR files\n",SNRpath);
    exit(1);
  }

  thisData=IFOdata; // restart from the first IFO
  while(thisData){
      NetSNR+=(thisData->SNR*thisData->SNR);
      thisData=thisData->next;
  }

  fprintf(snrout,"%4.2f\n",sqrt(NetSNR));
  fclose(snrout);
}

/**
* Fill the variables passed in vars with the parameters of the injection passed in event
* will over-write and destroy any existing parameters. Param vary type will be fixed
*/
void LALInferenceInjectionToVariables(SimInspiralTable *theEventTable, LALInferenceVariables *vars)
{
  UINT4 spinCheck=0;
  if(!vars) {
  XLALPrintError("Encountered NULL variables pointer");
  XLAL_ERROR_VOID(XLAL_EINVAL);
  }
  enforce_m1_larger_m2(theEventTable);
  REAL8 q = theEventTable->mass2 / theEventTable->mass1;
  if (q > 1.0) q = 1.0/q;

  REAL8 sx = theEventTable->spin1x;
  REAL8 sy = theEventTable->spin1y;
  REAL8 s1z = theEventTable->spin1z;

  REAL8 a_spin1 = sqrt(sx*sx + sy*sy + s1z*s1z);

  REAL8 theta_spin1, phi_spin1;
  if (a_spin1 == 0.0) {
    theta_spin1 = 0.0;
    phi_spin1 = 0.0;
  } else {
    theta_spin1 = acos(s1z / a_spin1);
    phi_spin1 = atan2(sy, sx);
    if (phi_spin1 < 0.0) phi_spin1 += 2.0*M_PI;
  }

  sx = theEventTable->spin2x;
  sy = theEventTable->spin2y;
  REAL8 s2z = theEventTable->spin2z;

  REAL8 a_spin2 = sqrt(sx*sx + sy*sy + s2z*s2z), theta_spin2, phi_spin2;
  if (a_spin2 == 0.0) {
    theta_spin2 = 0.0;
    phi_spin2 = 0.0;
  } else {
    theta_spin2 = acos(s2z / a_spin2);
    phi_spin2 = atan2(sy, sx);
    if (phi_spin2 < 0.0) phi_spin2 += 2.0*M_PI;
  }

  /* Check for presence of spin in the injection */
  if(a_spin1!=0.0 || a_spin2!=0.0) spinCheck=1;

  REAL8 psi = theEventTable->polarization;
  if (psi>=M_PI) psi -= M_PI;

  REAL8 injGPSTime = XLALGPSGetREAL8(&(theEventTable->geocent_end_time));

  REAL8 dist = theEventTable->distance;
  REAL8 inclination = theEventTable->inclination;
  REAL8 phase = theEventTable->coa_phase;
  REAL8 dec = theEventTable->latitude;
  REAL8 ra = theEventTable->longitude;

  Approximant injapprox = XLALGetApproximantFromString(theEventTable->waveform);
  LALPNOrder order = XLALGetOrderFromString(theEventTable->waveform);

  REAL8 m1=theEventTable->mass1;
  REAL8 m2=theEventTable->mass2;
  REAL8 chirpmass = theEventTable->mchirp;
  LALInferenceAddVariable(vars, "mass1", &m1, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "mass2", &m2, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "chirpmass", &chirpmass, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "asym_massratio", &q, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "time", &injGPSTime, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "distance", &dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "inclination", &inclination, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "theta_JN", &inclination, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "polarisation", &(psi), LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "phase", &phase, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "declination", &dec, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "rightascension", &ra, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "LAL_APPROXIMANT", &injapprox, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "LAL_PNORDER",&order,LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "LAL_AMPORDER", &(theEventTable->amp_order), LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  if(spinCheck){
      if (theEventTable->spin1x==0 && theEventTable->spin1y==0.0 && theEventTable->spin2x==0.0 && theEventTable->spin2y==0.0){
        LALInferenceAddVariable(vars, "spin1", &s1z, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "spin2", &s2z, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      }
      else{
        LALInferenceAddVariable(vars, "a_spin1", &a_spin1, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "a_spin2", &a_spin2, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "theta_spin1", &theta_spin1, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "theta_spin2", &theta_spin2, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "phi_spin1", &phi_spin1, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "phi_spin2", &phi_spin2, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      }
  }

}

void LALInferencePrintInjectionSample(LALInferenceRunState *runState)
{
    ProcessParamsTable *ppt=LALInferenceGetProcParamVal(runState->commandLine,"--inj");
    LALInferenceVariables backup;
    LALInferenceVariables injparams;
    memset(&injparams,0,sizeof(LALInferenceVariables));
    memset(&backup,0,sizeof(LALInferenceVariables));
    char *fname=NULL;
    char defaultname[]="injection_params.dat";
    FILE *outfile=NULL;
    if(!ppt) return;
    SimInspiralTable *injTable=NULL,*theEventTable=NULL;
    SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);

    ppt=LALInferenceGetProcParamVal(runState->commandLine,"--outfile");
    if(ppt) {
      fname = XLALCalloc((strlen(ppt->value)+255)*sizeof(char),1);
      sprintf(fname,"%s.injection",ppt->value);
    }
    else fname=defaultname;

    ppt=LALInferenceGetProcParamVal(runState->commandLine,"--event");
    if (ppt) {
      UINT4 event = atoi(ppt->value);
      UINT4 i;
      theEventTable = injTable;
      for (i = 0; i < event; i++) {
        theEventTable = theEventTable->next;
      }
      theEventTable->next = NULL;
    } else {
      theEventTable=injTable;
      theEventTable->next = NULL;
    }

    /* Save old variables */
    LALInferenceCopyVariables(runState->currentParams,&backup);
    //LALInferenceClearVariables(runState->currentParams);
    LALPNOrder *order=LALInferenceGetVariable(&backup,"LAL_PNORDER");
    Approximant *approx=LALInferenceGetVariable(&backup,"LAL_APPROXIMANT");
    /* Fill named variables */
    LALInferenceInjectionToVariables(theEventTable,runState->currentParams);

    /* If the time prior is stored in currentParams for the margtime likelihood, this will copy its range over */
    if(LALInferenceCheckVariable(&backup,"time_min"))
    {
            REAL8 time_min=LALInferenceGetREAL8Variable(&backup,"time_min");
            LALInferenceAddVariable(runState->currentParams,"time_min",&time_min,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
    }
    if(LALInferenceCheckVariable(&backup,"time_max"))
    {
            REAL8 time_max=LALInferenceGetREAL8Variable(&backup,"time_max");
            LALInferenceAddVariable(runState->currentParams,"time_max",&time_max,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
    }

    if(order && approx){
      /* Set the waveform to the one used in the analysis */
      LALInferenceRemoveVariable(runState->currentParams,"LAL_APPROXIMANT");
      LALInferenceRemoveVariable(runState->currentParams,"LAL_PNORDER");
      LALInferenceAddVariable(runState->currentParams,"LAL_PNORDER",order,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
      LALInferenceAddVariable(runState->currentParams,"LAL_APPROXIMANT",approx,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
    }
    REAL8 injPrior = runState->prior(runState,runState->currentParams);
    LALInferenceAddVariable(runState->currentParams,"logPrior",&injPrior,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    int errnum=0;
    REAL8 injL=0.;
    if ( (int) *approx == XLALGetApproximantFromString(theEventTable->waveform)){
      XLAL_TRY(injL = runState->likelihood(runState->currentParams, runState->data, runState->templt), errnum);
      if(errnum){
          fprintf(stderr,"ERROR: Cannot print injection sample. Received error code %s\n",XLALErrorString(errnum));
      }
    }
    LALInferenceAddVariable(runState->currentParams,"logL",(void *)&injL,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    if(LALInferenceCheckVariable(runState->algorithmParams,"logZnoise")){
        REAL8 tmp=injL-*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logZnoise");
        LALInferenceAddVariable(runState->currentParams,"deltalogL",(void *)&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    }
    LALInferenceIFOData *data=runState->data;
    while(data)
    {
        char tmpName[50];
        REAL8 tmp=data->loglikelihood - data->nullloglikelihood;
        sprintf(tmpName,"deltalogl%s",data->name);
        LALInferenceAddVariable(runState->currentParams,tmpName,&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
        data=data->next;
    }
    /* Save to file */
    outfile=fopen(fname,"w");
    if(!outfile) {fprintf(stderr,"ERROR: Unable to open file %s for injection saving\n",fname); exit(1);}
    LALInferenceSortVariablesByName(runState->currentParams);
    for(LALInferenceVariableItem *this=runState->currentParams->head; this; this=this->next)
        fprintf(outfile,"%s\t",this->name);
    fprintf(outfile,"\n");
    LALInferencePrintSample(outfile,runState->currentParams);
    fclose(outfile);
    
    /* Set things back the way they were */    
    //LALInferenceCopyVariables(&backup,runState->currentParams);
    //if(runState->currentParams && runState->currentParams->head) runState->likelihood(runState->currentParams,runState->data,runState->templt);
    return;
}

void enforce_m1_larger_m2(SimInspiralTable* injEvent){	
    /* Template generator assumes m1>=m2 thus we must enfore the same convention while injecting, otherwise spin2 will be assigned to mass1
    *        We also shift the phase by pi to be sure the same WF in injected 
    */
    REAL8 m1,m2,tmp;
    m1=injEvent->mass1;
    m2=injEvent->mass2;
   
    if (m1>=m2) return;
    else{
        fprintf(stdout, "Injtable has m1<m2. Flipping masses and spins in injection. Shifting phase by pi. \n");
        tmp=m1;
        injEvent->mass1=injEvent->mass2;
        injEvent->mass2=tmp;
        tmp=injEvent->spin1x;
        injEvent->spin1x=injEvent->spin2x;
        injEvent->spin2x=tmp;
        tmp=injEvent->spin1y;
        injEvent->spin1y=injEvent->spin2y;
        injEvent->spin2y=tmp;
        tmp=injEvent->spin1z;
        injEvent->spin1z=injEvent->spin2z;
        injEvent->spin2z=tmp;
	injEvent->coa_phase=injEvent->coa_phase+LAL_PI;
        }
    return ;
}

void LALInferenceSetupROQ(LALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine){
  
  LALStatus status;
  memset(&status,0,sizeof(status));
  UINT4 Nifo=0;
  LALInferenceIFOData *thisData=IFOdata;
  UINT4 i;
  UINT4 q=0;
  UINT4 event=0;
  char *chartmp=NULL;
  ProcessParamsTable *procparam=NULL,*ppt=NULL;
  SimInspiralTable *injTable=NULL;
  FILE *tempfp;
  unsigned int n_basis,n_samples,time_steps;
  n_basis = 965;//TODO: have it read from file or from command line.
  n_samples = 31489;
  REAL8 delta_tc = 0.0001;
  REAL8 dt=0.1;
  //REAL8 tc=0;
  LIGOTimeGPS GPStrig;
  REAL8 endtime=0.0;
  REAL8 timeMin=0.0,timeMax=0.0;
  const UINT4 nameLength=FILENAME_MAX;
  char filename[nameLength];
  FILE *out;
	char tmp[128];

  
  while(thisData){
    thisData=thisData->next;
    Nifo++;
  }
  
  procparam=LALInferenceGetProcParamVal(commandLine,"--inj");
  if(procparam){
    SimInspiralTableFromLIGOLw(&injTable,procparam->value,0,0);
    if(!injTable){
      fprintf(stderr,"Unable to open injection file(LALInferenceReadData) %s\n",procparam->value);
      exit(1);
    }
    procparam=LALInferenceGetProcParamVal(commandLine,"--event");
    if(procparam) {
      event=atoi(procparam->value);
      while(q<event) {q++; injTable=injTable->next;}
    }
    else if ((procparam=LALInferenceGetProcParamVal(commandLine,"--event-id")))
    {
      while(injTable)
      {
        if(injTable->event_id->id == (UINT4)atoi(procparam->value)) break;
        else injTable=injTable->next;
      }
      if(!injTable){
        fprintf(stderr,"Error, cannot find simulation id %s in injection file\n",procparam->value);
        exit(1);
      }
    }
  }
  
  if(LALInferenceGetProcParamVal(commandLine,"--trigtime")){
    procparam=LALInferenceGetProcParamVal(commandLine,"--trigtime");
    LALStringToGPS(&status,&GPStrig,procparam->value,&chartmp);
  }
  else{
    if(injTable) memcpy(&GPStrig,&(injTable->geocent_end_time),sizeof(GPStrig));
    else {
      fprintf(stderr,"Error: No trigger time specifed and no injection given \n");
      exit(1);
    }
  }
  
  endtime=XLALGPSGetREAL8(&GPStrig);
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
  if(ppt){
    dt=atof(ppt->value);
  }
  ppt=LALInferenceGetProcParamVal(commandLine,"--delta_tc");
  if(ppt){
    delta_tc=atof(ppt->value);
  }
  
  timeMin=endtime-dt-0.022; timeMax=endtime+dt+0.022;
  
  timeMin -= XLALGPSGetREAL8(&IFOdata[0].whiteFreqData->epoch);
  timeMax -= XLALGPSGetREAL8(&IFOdata[0].whiteFreqData->epoch);
  
  time_steps = (unsigned int)((timeMax-timeMin)/delta_tc)+1;
	
	if(LALInferenceGetProcParamVal(commandLine,"--roqtime_steps")){
		ppt=LALInferenceGetProcParamVal(commandLine,"--roqtime_steps");
		tempfp = fopen (ppt->value,"r");
		fscanf (tempfp, "%u", &time_steps);
		fscanf (tempfp, "%u", &n_basis);
		fscanf (tempfp, "%u", &n_samples);
	}

  printf("endtime = %f, timeMin = %f, timeMax = %f\n", endtime, timeMin, timeMax);
  printf("time steps = %d\n", time_steps);
  
  double deltaF = IFOdata[0].oneSidedNoisePowerSpectrum->deltaF; //assumes same deltaF for all IFOs
  /*
  gsl_matrix_complex *vandermonde_matrix=NULL;
  gsl_matrix_complex *rb_matrix=NULL;
  
  gsl_matrix_complex *whitened_data_matrix = gsl_matrix_complex_calloc(n_samples, time_steps); //IFO data / PSD
  gsl_vector_complex *whitened_data = gsl_vector_complex_calloc(n_samples);
  gsl_matrix_complex *E_matrix = gsl_matrix_complex_calloc(n_basis, time_steps);
  gsl_vector_complex *exp_2pi_i_f_tc = gsl_vector_complex_calloc(n_samples);
  gsl_complex alpha;
  gsl_complex beta;
  gsl_complex wd; // element of whiteFreqData
  gsl_complex shifted_data_element;
  */
  
  /*
  if(LALInferenceGetProcParamVal(commandLine,"--roqvandermonde")){
    ppt=LALInferenceGetProcParamVal(commandLine,"--roqvandermonde");
    
    vandermonde_matrix = gsl_matrix_complex_calloc(n_basis, n_basis);
    tempfp = fopen(ppt->value, "rb");
    gsl_matrix_complex_fread(tempfp, vandermonde_matrix);
  }
  
  if(LALInferenceGetProcParamVal(commandLine,"--roqrb")){
    ppt=LALInferenceGetProcParamVal(commandLine,"--roqrb");
    
    rb_matrix = gsl_matrix_complex_calloc(n_samples, n_basis);
    tempfp = fopen(ppt->value, "rb");
    gsl_matrix_complex_fread(tempfp, rb_matrix);
  }
  */
  
  //if(LALInferenceGetProcParamVal(commandLine,"--roqnodes") && LALInferenceGetProcParamVal(commandLine,"--roqvandermonde") && LALInferenceGetProcParamVal(commandLine,"--roqrb")){
  if(LALInferenceGetProcParamVal(commandLine,"--roqnodes")){
    
    ppt=LALInferenceGetProcParamVal(commandLine,"--roqnodes");
    for (i=0;i<Nifo;i++) {
			IFOdata[i].roqData = XLALCalloc(1, sizeof(LALInferenceROQData));
      IFOdata[i].roqData->frequencyNodes = gsl_vector_calloc(n_basis);
      tempfp = fopen(ppt->value, "rb");
      gsl_vector_fread(tempfp, IFOdata[i].roqData->frequencyNodes);
    }
    
    //printf("n_basis=%d\ttime_steps=%d\n", n_basis, time_steps);
		
    //ppt=LALInferenceGetProcParamVal(commandLine,"--roqweights");
    for (i=0;i<Nifo;i++) {
			sprintf(tmp,"--%s-roqweights",IFOdata[i].name);
			ppt=LALInferenceGetProcParamVal(commandLine,tmp);
      IFOdata[i].roqData->weights = gsl_matrix_complex_calloc(n_basis, time_steps);
      tempfp = fopen(ppt->value, "rb");
      gsl_matrix_complex_fread(tempfp, IFOdata[i].roqData->weights);
    }
    
    //printf("time steps = %d\n", (int)IFOdata[0].roqData->weights->size2);
  
    for (i=0;i<Nifo;i++) {
      
      IFOdata[i].roqData->trigtime = endtime;
      IFOdata[i].roqData->time_weights_width = timeMax-timeMin;
      
      /*** compute the weights ***/
      /*
      GSL_SET_COMPLEX(&alpha, 4.*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF, 0);
      GSL_SET_COMPLEX(&beta, 0, 0);
      
      

      snprintf(filename, nameLength, "%s-ROQFreqData.dat", IFOdata[i].name);
      out = fopen(filename, "w");
      for(unsigned int k = 0; k < n_samples; k++){
            fprintf(out,"%g %g %g\n", IFOdata[i].freqData->deltaF*(k + (unsigned int)(IFOdata[i].fLow/deltaF)),
              creal(IFOdata[i].freqData->data->data[k + (unsigned int)(IFOdata[i].fLow/deltaF)]),
              cimag(IFOdata[i].freqData->data->data[k + (unsigned int)(IFOdata[i].fLow/deltaF)]));
      }
      fclose(out);
      //REAL8 deltaT = IFOdata[i].timeData->deltaT;
      
      for(unsigned int jj = 0; jj < time_steps; jj++){
        
        tc = timeMin + 2.0*(dt+0.022)*jj / time_steps;
        //printf("%f = %f + 2*%f*%d / %d;\n",tc,timeMin,dt,jj,time_steps);
        
        for(unsigned int ii=0; ii < n_samples; ii++){
          gsl_vector_complex_set(exp_2pi_i_f_tc, ii, gsl_complex_polar (1, (IFOdata[i].fLow/deltaF+ii)*deltaF*tc*2*LAL_PI));
        }
        
        for(unsigned int k = 0; k < n_samples; k++){
          
          if(!isnan(creal(IFOdata[i].freqData->data->data[k + (unsigned int)(IFOdata[i].fLow/deltaF)])) &&
             IFOdata[i].oneSidedNoisePowerSpectrum->data->data[k + (unsigned int)(IFOdata[i].fLow/deltaF)] != 0.0 ){
            //GSL_SET_COMPLEX(&wd, creal(IFOdata[i].whiteFreqData->data->data[k + (unsigned int)(IFOdata[i].fLow/deltaF)]), cimag(IFOdata[i].whiteFreqData->data->data[k + (unsigned int)(IFOdata[i].fLow/deltaF)]));
            GSL_SET_COMPLEX(&wd, creal(IFOdata[i].freqData->data->data[k + (unsigned int)(IFOdata[i].fLow/deltaF)])/
                                       IFOdata[i].oneSidedNoisePowerSpectrum->data->data[k + (unsigned int)(IFOdata[i].fLow/deltaF)],
                            cimag(IFOdata[i].freqData->data->data[k + (unsigned int)(IFOdata[i].fLow/deltaF)])/
                                  IFOdata[i].oneSidedNoisePowerSpectrum->data->data[k + (unsigned int)(IFOdata[i].fLow/deltaF)]);
          }else{
            GSL_SET_COMPLEX(&wd, 0.0, 0.0);
          }
          
          shifted_data_element = gsl_complex_mul (wd, gsl_vector_complex_get(exp_2pi_i_f_tc, k));
          
          //take conjugate
          shifted_data_element = gsl_complex_conjugate (shifted_data_element);
          
          gsl_vector_complex_set (whitened_data, k, shifted_data_element);
          
        }
        
        gsl_matrix_complex_set_col (whitened_data_matrix, jj, whitened_data);
        
      }
      
      printf("Computing weights for %s\n",IFOdata[i].name);
      
      gsl_blas_zgemm (CblasTrans, CblasNoTrans, alpha, rb_matrix, whitened_data_matrix, beta, E_matrix);
      
      GSL_SET_COMPLEX (&alpha, 1, 0);
      
      gsl_blas_zgemm (CblasTrans, CblasNoTrans, alpha, vandermonde_matrix, E_matrix, beta, IFOdata[i].roqData->weights);
      
      printf("Weights have been computed for %s\n",IFOdata[i].name);
      */

      if (LALInferenceGetProcParamVal(commandLine, "--data-dump")) {
        snprintf(filename, nameLength, "%s-ROQWeights.dat", IFOdata[i].name);
        out = fopen(filename, "w");
        for(unsigned int size2 = 0; size2 < IFOdata[i].roqData->weights->size2; size2++){
          for(unsigned int size1 = 0; size1 < IFOdata[i].roqData->weights->size1; size1++){
            fprintf(out,"(%g+%gj)\t",GSL_REAL(gsl_matrix_complex_get(IFOdata[i].roqData->weights,size1,size2)),GSL_IMAG(gsl_matrix_complex_get(IFOdata[i].roqData->weights,size1,size2)));
          }
          fprintf(out,"\n");
        }
      fclose(out);
      }
      //printf("IFOdata[%d].whiteFreqData->epoch=%e\n",i,XLALGPSGetREAL8(&IFOdata[i].whiteFreqData->epoch));
      //printf("timeMin=%e\tendtime=%e\ttimeMax=%e\n",timeMin,endtime,timeMax);
      //printf("---------\n");
      
      // compute int_f_7_over_3
      IFOdata[i].roqData->int_f_7_over_3 = 0;
      for(unsigned int kk = 0; kk < n_samples; kk++){
        if(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[kk + (unsigned int)(IFOdata[i].fLow/deltaF)] != 0.0){
          IFOdata[i].roqData->int_f_7_over_3 += 4.*deltaF*pow((IFOdata[i].fLow + kk*deltaF), -7./3.) / IFOdata[i].oneSidedNoisePowerSpectrum->data->data[kk + (unsigned int)(IFOdata[i].fLow/deltaF)];
        }
      }
    }
    
  }
  
  //if(vandermonde_matrix) gsl_matrix_complex_free(vandermonde_matrix);
  //if(rb_matrix) gsl_matrix_complex_free(rb_matrix);
  //if(rb_matrix){
  //  gsl_matrix_complex_free(whitened_data_matrix);
  //  gsl_matrix_complex_free(E_matrix);
  //}
  
}
