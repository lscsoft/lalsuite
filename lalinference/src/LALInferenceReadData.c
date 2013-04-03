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

#define LAL_USE_OLD_COMPLEX_STRUCTS
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
	interp=calloc(minLength,sizeof(struct fvec)); /* Initialise array */
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
			interp=realloc(interp,(fileLength+minLength)*sizeof(struct fvec));
			fileLength+=minLength;
		}
	}
	interp[i].f=0; interp[i].x=0;
	fileLength=i+1;
	interp=realloc(interp,fileLength*sizeof(struct fvec)); /* Resize array */
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
void InjectTaylorF2(LALInferenceIFOData *IFOdata, SimInspiralTable *inj_table, ProcessParamsTable *commandLine);

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

static REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length);
static void makeWhiteData(LALInferenceIFOData *IFOdata);
static void PrintSNRsToFile(LALInferenceIFOData *IFOdata , SimInspiralTable *inj_table);

static REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length)
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

/** Parse the command line looking for options of the kind --ifo H1 --H1-channel H1:LDAS_STRAIN --H1-cache H1.cache --H1-flow 40.0 --H1-fhigh 4096.0 --H1-timeslide 100.0 ...
    It is necessary to use this method instead of the old method for the pipeline to work in DAX mode. Warning: do not mix options between
    the old and new style.
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
            *ifos=realloc(*ifos,*N*sizeof(char *));
            (*ifos)[*N-1]=strdup(this->value);
        }
    }
    *caches=calloc(*N,sizeof(char *));
    *channels=calloc(*N,sizeof(char *));
    *flows=calloc(*N,sizeof(REAL8));
    *fhighs=calloc(*N,sizeof(REAL8));
    *timeslides=calloc(*N,sizeof(REAL8));
    /* For each IFO, fetch the other options if available */
    for(i=0;i<*N;i++)
    {
        /* Cache */
        sprintf(tmp,"--%s-cache",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        if(!this){fprintf(stderr,"ERROR: Must specify a cache file for %s with --%s-cache\n",(*ifos)[i],(*ifos)[i]); exit(1);}
        (*caches)[i]=strdup(this->value);
        
        /* Channel */
        sprintf(tmp,"--%s-channel",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*channels)[i]=strdup(this?this->value:"Unknown channel");

        /* flow */
        sprintf(tmp,"--%s-flow",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*flows)[i]=strdup(this?this->value:LALINFERENCE_DEFAULT_FLOW);
        
        /* fhigh */
        sprintf(tmp,"--%s-fhigh",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*fhighs)[i]=this?strdup(this->value):NULL;

        /* timeslides */
        sprintf(tmp,"--%s-timeslide",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*timeslides)[i]=strdup(this?this->value:"0.0");
    }
    return(1);
}

#define USAGE "\
 --ifo IFO1 [--ifo IFO2 ...]    IFOs can be H1,L1,V1\n\
 --IFO1-cache cache1 [--IFO2-cache2 cache2 ...]    LAL cache files (LALLIGO, LALAdLIGO, LALVirgo to simulate these detectors)\n\
 --psdstart GPStime             GPS start time of PSD estimation data\n\
 --psdlength length             length of PSD estimation data in seconds\n\
 --seglen length                length of segments for PSD estimation and analysis in seconds\n\
 --trigtime GPStime             GPS time of the trigger to analyse\n\
(--srate rate)                  Downsample data to rate in Hz (4096.0,)\n\
(--injectionsrate rate)         Downsample injection signal to rate in Hz (--srate)\n\
(--IFO1-flow freq1 [--IFO2-flow freq2 ...])      Specify lower frequency cutoff for overlap integral (40.0)\n\
(--IFO1-fhigh freq1 [--IFO2-fhigh freq2 ...])     Specify higher frequency cutoff for overlap integral (2048.0)\n\
(--IFO1-channel chan1 [--IFO2-channel chan2 ...])   Specify channel names when reading cache files\n\
(--dataseed number)             Specify random seed to use when generating data\n\
(--lalinspiralinjection)      Enables injections via the LALInspiral package\n\
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

    /* Check if the new style command line arguments are used */
    INT4 dataOpts=getDataOptionsByDetectors(commandLine, &IFOnames, &caches, &channels, &fLows , &fHighs, &timeslides, &Nifo);
    /* Check for options if not given in the new style */
    if(!dataOpts){
        if(!LALInferenceGetProcParamVal(commandLine,"--cache")||!(LALInferenceGetProcParamVal(commandLine,"--IFO")||LALInferenceGetProcParamVal(commandLine,"--ifo")))
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
    dE1.location[0] = dE2.location[0] = dE3.location[0] = 4.5464e6;
    dE1.location[1] = dE2.location[1] = dE3.location[1] = 8.4299e5;
    dE1.location[2] = dE2.location[2] = dE3.location[2] = 4.3786e6;
    sprintf(dE1.frDetector.name,"ET-1");
    sprintf(dE1.frDetector.prefix,"E1");
    dE1.response[0][0] = 0.1666;
    dE1.response[1][1] = -0.2484;
    dE1.response[2][2] = 0.0818;
    dE1.response[0][1] = dE1.response[1][0] = -0.2188;
    dE1.response[0][2] = dE1.response[2][0] = -0.1300;
    dE1.response[1][2] = dE1.response[2][1] = 0.2732;
    sprintf(dE2.frDetector.name,"ET-2");
    sprintf(dE2.frDetector.prefix,"E2");
    dE2.response[0][0] = -0.1992;
    dE2.response[1][1] = 0.4234;
    dE2.response[2][2] = 0.0818;
    dE2.response[0][1] = dE2.response[1][0] = -0.0702;
    dE2.response[0][2] = dE2.response[2][0] = 0.2189;
    dE2.response[1][2] = dE2.response[2][1] = -0.0085;
    sprintf(dE3.frDetector.name,"ET-3");
    sprintf(dE3.frDetector.prefix,"E3");
    dE3.response[0][0] = 0.0326;
    dE3.response[1][1] = -0.1750;
    dE3.response[2][2] = 0.1423;
    dE3.response[0][1] = dE3.response[1][0] = 0.2891;
    dE3.response[0][2] = dE3.response[2][0] = -0.0889;
    dE3.response[1][2] = dE3.response[2][1] = -0.2647;  

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

    IFOdata=headIFO=calloc(sizeof(LALInferenceIFOData),Nifo);
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
        else {
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

    for(i=0;i<Nifo;i++) {
        IFOdata[i].fLow=fLows?atof(fLows[i]):defaultFLow; 
        if(fHighs) IFOdata[i].fHigh=fHighs[i]?atof(fHighs[i]):(SampleRate/2.0-(1.0/SegmentLength));
        else IFOdata[i].fHigh=(SampleRate/2.0-(1.0/SegmentLength));
        strncpy(IFOdata[i].name, IFOnames[i], DETNAMELEN);
        IFOdata[i].STDOF = 4.0 / M_PI * nSegs;
        fprintf(stderr, "Detector %s will run with %g DOF if Student's T likelihood used.\n",
                IFOdata[i].name, IFOdata[i].STDOF);
    }

    /* Only allocate this array if there weren't channels read in from the command line */
    if(!dataOpts && !Nchannel) channels=calloc(Nifo,sizeof(char *));
    for(i=0;i<Nifo;i++) {
        if(!dataOpts && !Nchannel) channels[i]=malloc(VARNAME_MAX);
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
            IFOdata[i].detector=calloc(1,sizeof(LALDetector));
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
            IFOdata[i].detector=calloc(1,sizeof(LALDetector));
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
            IFOdata[i].detector=calloc(1,sizeof(LALDetector));
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
            IFOdata[i].detector=calloc(1,sizeof(LALDetector));
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
            IFOdata[i].detector=malloc(sizeof(LALDetector));
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
            IFOdata[i].detector=malloc(sizeof(LALDetector));
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
            IFOdata[i].detector=malloc(sizeof(LALDetector));
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
        /* Create FFT plans */
        IFOdata[i].timeToFreqFFTPlan = XLALCreateForwardREAL8FFTPlan((UINT4) seglen, 0 );
        if(!IFOdata[i].timeToFreqFFTPlan) XLAL_ERROR_NULL(XLAL_EFUNC);
        IFOdata[i].freqToTimeFFTPlan = XLALCreateReverseREAL8FFTPlan((UINT4) seglen,0);
        if(!IFOdata[i].freqToTimeFFTPlan) XLAL_ERROR_NULL(XLAL_EFUNC);		
        /* Setup windows */
        IFOdata[i].window=XLALCreateTukeyREAL8Window(seglen,(REAL8)2.0*padding*SampleRate/(REAL8)seglen);
        if(!IFOdata[i].window) XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    /* Trigger time = 2 seconds before end of segment (was 1 second, but Common Inputs for The Events are -6 +2*/
    memcpy(&segStart,&GPStrig,sizeof(LIGOTimeGPS));
    XLALGPSAdd(&segStart,-SegmentLength+2);


    /* Read the PSD data */
    for(i=0;i<Nifo;i++) {
        memcpy(&(IFOdata[i].epoch),&segStart,sizeof(LIGOTimeGPS));
        /* Check to see if an interpolation file is specified */
        interpFlag=0;
        interp=NULL;
        if(strstr(caches[i],"interp:")==caches[i]){
            /* Extract the file name */
            char *interpfilename=&(caches[i][7]);
            printf("Looking for interpolation file %s\n",interpfilename);
            interpFlag=1;
            interp=interpFromFile(interpfilename);
        }    
        /* Check if fake data is requested */
        if(interpFlag || (!(strcmp(caches[i],"LALLIGO") && strcmp(caches[i],"LALVirgo") && strcmp(caches[i],"LALGEO") && strcmp(caches[i],"LALEGO")
                        && strcmp(caches[i],"LALAdLIGO"))))
        {
            //FakeFlag=1; - set but not used
            if (!LALInferenceGetProcParamVal(commandLine,"--dataseed")){
                fprintf(stderr,"Error: You need to specify a dataseed when generating data with --dataseed <number>.\n\
                        (--dataseed 0 uses a non-reproducible number from the system clock, and no parallel run is then possible.)\n" );
                exit(-1);
            }
            datarandparam=XLALCreateRandomParams(dataseed?dataseed+(int)i:dataseed);
            if(!datarandparam) XLAL_ERROR_NULL(XLAL_EFUNC);
            /* Selection of the noise curve */
            if(!strcmp(caches[i],"LALLIGO")) {PSD = &LALLIGOIPsd; scalefactor=9E-46;}
            if(!strcmp(caches[i],"LALVirgo")) {PSD = &LALVIRGOPsd; scalefactor=1.0;}
            if(!strcmp(caches[i],"LALGEO")) {PSD = &LALGEOPsd; scalefactor=1E-46;}
            if(!strcmp(caches[i],"LALEGO")) {PSD = &LALEGOPsd; scalefactor=1.0;}
            if(!strcmp(caches[i],"LALAdLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 1E-49;}
            if(interpFlag) {PSD=NULL; scalefactor=1.0;}
            //if(!strcmp(caches[i],"LAL2kLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 36E-46;}
            if(PSD==NULL && !interpFlag) {fprintf(stderr,"Error: unknown simulated PSD: %s\n",caches[i]); exit(-1);}


            IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)
                XLALCreateREAL8FrequencySeries("spectrum",&GPSstart,0.0,
                        (REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
            if(!IFOdata[i].oneSidedNoisePowerSpectrum) XLAL_ERROR_NULL(XLAL_EFUNC);
            for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
            {
                MetaNoiseFunc(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,interp,PSD);
                //PSD(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF);
                IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]*=scalefactor;
            }
            IFOdata[i].freqData = (COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("stilde",&segStart,0.0,IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,&lalDimensionlessUnit,seglen/2 +1);
            if(!IFOdata[i].freqData) XLAL_ERROR_NULL(XLAL_EFUNC);

            /* Create the fake data */
            int j_Lo = (int) IFOdata[i].fLow/IFOdata[i].freqData->deltaF;
            if(LALInferenceGetProcParamVal(commandLine,"--0noise")){
                for(j=j_Lo;j<IFOdata[i].freqData->data->length;j++){
                    IFOdata[i].freqData->data->data[j].re=IFOdata[i].freqData->data->data[j].im=0.0;
                }
            } else {
                for(j=j_Lo;j<IFOdata[i].freqData->data->length;j++){
                    IFOdata[i].freqData->data->data[j].re=XLALNormalDeviate(datarandparam)*(0.5*sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]/IFOdata[i].freqData->deltaF));
                    IFOdata[i].freqData->data->data[j].im=XLALNormalDeviate(datarandparam)*(0.5*sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]/IFOdata[i].freqData->deltaF));
                }
            }
            IFOdata[i].freqData->data->data[0].re=0; 			IFOdata[i].freqData->data->data[0].im=0;
            const char timename[]="timeData";
            IFOdata[i].timeData=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries(timename,&segStart,0.0,(REAL8)1.0/SampleRate,&lalDimensionlessUnit,(size_t)seglen);
            if(!IFOdata[i].timeData) XLAL_ERROR_NULL(XLAL_EFUNC);
            XLALREAL8FreqTimeFFT(IFOdata[i].timeData,IFOdata[i].freqData,IFOdata[i].freqToTimeFFTPlan);
            if(*XLALGetErrnoPtr()) printf("XLErr: %s\n",XLALErrorString(*XLALGetErrnoPtr()));
            XLALDestroyRandomParams(datarandparam);
        }
        else{ /* Not using fake data, load the data from a cache file */
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
                PSDtimeSeries=readTseries(caches[i],channels[i],GPSstart,PSDdatalength);
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

                XLALDestroyREAL8TimeSeries(PSDtimeSeries);
            }

            /* Read the data segment */
            LIGOTimeGPS truesegstart=segStart;
            if(Ntimeslides) {
                REAL4 deltaT=-atof(timeslides[i]);
                XLALGPSAdd(&segStart, deltaT);
                fprintf(stderr,"Slid %s by %f s from %10.10lf to %10.10lf\n",IFOnames[i],deltaT,truesegstart.gpsSeconds+1e-9*truesegstart.gpsNanoSeconds,segStart.gpsSeconds+1e-9*segStart.gpsNanoSeconds);
            }
            IFOdata[i].timeData=readTseries(caches[i],channels[i],segStart,SegmentLength);
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
                IFOdata[i].freqData->data->data[j].re/=sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
                IFOdata[i].freqData->data->data[j].im/=sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
                IFOdata[i].windowedTimeData->data->data[j] /= sqrt(IFOdata[i].window->sumofsquares / IFOdata[i].window->data->length);
            }
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

        if (LALInferenceGetProcParamVal(commandLine, "--data-dump")) {
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

    for (i=0;i<Nifo;i++) IFOdata[i].SNR=0.0; //SNR of the injection ONLY IF INJECTION. Set to 0.0 by default.

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
		SNRpath = calloc(strlen(ppt->value)+1,sizeof(char));
		memcpy(SNRpath,ppt->value,strlen(ppt->value)+1);
                fprintf(stdout,"Writing SNRs in %s\n",SNRpath)     ;

	}
	Ninj=SimInspiralTableFromLIGOLw(&injTable,LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
	REPORTSTATUS(&status);
	printf("Ninj %d\n", Ninj);
	if(Ninj<event) fprintf(stderr,"Error reading event %d from %s\n",event,LALInferenceGetProcParamVal(commandLine,"--inj")->value);
	while(i<event) {i++; injTable = injTable->next;} /* Select event */
	injEvent = injTable;
	injEvent->next = NULL;
	
	//memset(&InjectGW,0,sizeof(InjectGW));
	Approximant injapprox;
	injapprox = XLALGetApproximantFromString(injTable->waveform);
        if( (int) injapprox == XLAL_FAILURE)
          ABORTXLAL(&status);
	printf("Injecting approximant %i: %s\n", injapprox, injTable->waveform);
	REPORTSTATUS(&status);
	//LALGenerateInspiral(&status,&InjectGW,injTable,&InjParams);
	//if(status.statusCode!=0) {fprintf(stderr,"Error generating injection!\n"); REPORTSTATUS(&status); }
	/* Check for frequency domain injection (TF2 only at present) */
	if(strstr(injTable->waveform,"TaylorF2"))
	{ printf("Injecting TaylorF2 in the frequency domain...\n");
	 InjectTaylorF2(IFOdata, injTable, commandLine);
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
		
		for(i=0;i<resp->data->length;i++) {resp->data->data[i].re=(REAL4)1.0; resp->data->data[i].im=0.0;}
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
      LALSimInspiralWaveformFlags *waveFlags = XLALSimInspiralCreateWaveformFlags();
      LALSimInspiralSpinOrder spinO = -1;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-spinOrder")) {
        spinO = atoi(LALInferenceGetProcParamVal(commandLine,"--inj-spinOrder")->value);
        XLALSimInspiralSetSpinOrder(waveFlags, spinO);
        fprintf(stdout,"Injection (twice) PN spin order set to %i\n",spinO);
      }
      LALSimInspiralTidalOrder tideO = -1;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-tidalOrder")) {
        tideO = atoi(LALInferenceGetProcParamVal(commandLine,"--inj-tidalOrder")->value);
        XLALSimInspiralSetTidalOrder(waveFlags, tideO);
        fprintf(stdout,"Injection (twice) PN tidal order set to %i\n",tideO);
      }
      LALSimInspiralTestGRParam *nonGRparams = NULL;
      
      XLALSimInspiralChooseTDWaveform(&hplus, &hcross, injEvent->coa_phase, 1.0/InjSampleRate,
                                      injEvent->mass1*LAL_MSUN_SI, injEvent->mass2*LAL_MSUN_SI, injEvent->spin1x,
                                      injEvent->spin1y, injEvent->spin1z, injEvent->spin2x, injEvent->spin2y,
                                      injEvent->spin2z, injEvent->f_lower, 0., injEvent->distance*LAL_PC_SI * 1.0e6,
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
	for(SNR=0.0,j=thisData->fLow/injF->deltaF;j<injF->data->length;j++){
	  SNR+=pow(injF->data->data[j].re,2.0)/thisData->oneSidedNoisePowerSpectrum->data->data[j];
	  SNR+=pow(injF->data->data[j].im,2.0)/thisData->oneSidedNoisePowerSpectrum->data->data[j];
	}
        SNR*=4.0*injF->deltaF;
    }
    thisData->SNR=sqrt(SNR);
    NetworkSNR+=SNR;

    if (!(SNRpath==NULL)){ /* If the user provided a path with --snrpath store a file with injected SNRs */
      PrintSNRsToFile(IFOdata , injTable);
    }
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
	thisData->freqData->data->data[j].re+=injF->data->data[j].re/WinNorm;
	thisData->freqData->data->data[j].im+=injF->data->data[j].im/WinNorm;
	fprintf(file, "%lg %lg \t %lg\n", thisData->freqData->deltaF*j, injF->data->data[j].re, injF->data->data[j].im);
      }
      fclose(file);
    
      XLALDestroyREAL8TimeSeries(inj8Wave);
      XLALDestroyCOMPLEX16FrequencySeries(injF);
      thisData=thisData->next;
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
  {
    RAT4 negOne = { -1, 0 };
    LALUnit unit;
    LALUnitPair pair;
    pair.unitOne = &lalADCCountUnit;
    pair.unitTwo = &lalStrainUnit;
    LALUnitRaise( status->statusPtr, &unit, pair.unitTwo, &negOne );
    CHECKSTATUSPTR( status );
    pair.unitTwo = &unit;
    LALUnitMultiply( status->statusPtr, &(detector.transfer->sampleUnits),
        &pair );
    CHECKSTATUSPTR( status );
  }

  /* invert the response function to get the transfer function */
  LALCCreateVector( status->statusPtr, &( detector.transfer->data ),
      resp->data->length );
  CHECKSTATUSPTR( status );

  LALCCreateVector( status->statusPtr, &unity, resp->data->length );
  CHECKSTATUSPTR( status );
  for ( k = 0; k < resp->data->length; ++k )
  {
    unity->data[k].re = 1.0;
    unity->data[k].im = 0.0;
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

void InjectTaylorF2(LALInferenceIFOData *IFOdata, SimInspiralTable *inj_table, ProcessParamsTable *commandLine)
///*-------------- Inject in Frequency domain -----------------*/
{
    /* Inject a gravitational wave into the data in the frequency domain */ 
    LALStatus status;
    memset(&status,0,sizeof(LALStatus));
    REAL8 mc=0.0;
    Approximant injapprox;
    LALPNOrder phase_order=-1;
    LALPNOrder amp_order=-1;
    injapprox = XLALGetApproximantFromString(inj_table->waveform);
    if( (int) injapprox == XLAL_FAILURE)
        ABORTXLAL(&status);
    phase_order = XLALGetOrderFromString(inj_table->waveform);
    if ( (int) phase_order == XLAL_FAILURE)
        ABORTXLAL(&status);
    LALInferenceVariables *modelParams=NULL;
    LALInferenceIFOData * tmpdata=IFOdata;
    REAL8 eta =0.0;
    REAL8 startPhase = 0.0;
    REAL8 inclination = 0.0;
    REAL8 distance=0.0;
    REAL8 longitude=0.0;
    REAL8 latitude=0.0;
    REAL8 polarization=0.0;
    REAL8 injtime=0.0;
   
    
    tmpdata->modelParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	modelParams=tmpdata->modelParams;
    memset(modelParams,0,sizeof(LALInferenceVariables));

    eta = inj_table->eta;
    mc=inj_table->mchirp;
    startPhase = inj_table->coa_phase;
    inclination = inj_table->inclination;
    distance=inj_table->distance;
    longitude=inj_table->longitude;
    latitude=inj_table->latitude;
    polarization=inj_table->polarization;
    injtime=(REAL8) inj_table->geocent_end_time.gpsSeconds + (REAL8) inj_table->geocent_end_time.gpsNanoSeconds*1.0e-9;
    amp_order=(LALPNOrder) inj_table->amp_order;

    LALInferenceAddVariable(tmpdata->modelParams, "chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(tmpdata->modelParams, "phase",&startPhase,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);  
    LALInferenceAddVariable(tmpdata->modelParams, "rightascension",&longitude,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);  
    LALInferenceAddVariable(tmpdata->modelParams, "declination",&latitude,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);  
    LALInferenceAddVariable(tmpdata->modelParams, "polarisation",&polarization,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);  
    LALInferenceAddVariable(tmpdata->modelParams, "time",&injtime,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(tmpdata->modelParams, "inclination",&inclination,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(tmpdata->modelParams, "massratio",&eta,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(tmpdata->modelParams, "distance",&distance,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(tmpdata->modelParams, "LAL_APPROXIMANT",&injapprox,LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(tmpdata->modelParams, "LAL_PNORDER",&phase_order,LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(tmpdata->modelParams, "LAL_AMPORDER",&amp_order,LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);

      REAL8 lambda1 = 0.;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")) {
        lambda1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")->value);
        fprintf(stdout,"Injection lambda1 set to %f\n",lambda1);
        LALInferenceAddVariable(tmpdata->modelParams, "lambda1",&lambda1,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
      }
      REAL8 lambda2 = 0.;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")) {
        lambda2= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")->value);
        fprintf(stdout,"Injection lambda2 set to %f\n",lambda2);
        LALInferenceAddVariable(tmpdata->modelParams, "lambda2",&lambda2,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
      }
      
    LALSimInspiralSpinOrder spinO = LAL_SIM_INSPIRAL_SPIN_ORDER_ALL;

    if(LALInferenceGetProcParamVal(commandLine, "--inj-spinOrder")) {
        spinO = atoi(LALInferenceGetProcParamVal(commandLine, "--inj-spinOrder")->value);
        LALInferenceAddVariable(tmpdata->modelParams, "spinO", &spinO,   LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    }
    else
        fprintf(stdout,"No --inj-spinOrder option given. Injecting the highest spin order for this waveform!\n");
    LALSimInspiralTidalOrder tideO = LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL;

    if(LALInferenceGetProcParamVal(commandLine, "--inj-tidalOrder")) {
        tideO = atoi(LALInferenceGetProcParamVal(commandLine, "--inj-tidalOrder")->value);
        LALInferenceAddVariable(tmpdata->modelParams, "tideO", &tideO,   LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    }
    else
        fprintf(stdout,"No --inj-tidalOrder option given. Injecting the highest tidal order for this waveform!\n");
    fprintf(stdout,"injectTaylorF2 will run using Approximant %i (%s), phase order %i, amp order %i, spinOrder %i TidalOrder %i in the Frequency domain.\n",injapprox,XLALGetStringFromApproximant(injapprox),phase_order,amp_order,(int) spinO,(int) tideO);
    
    COMPLEX16FrequencySeries *freqModelhCross=NULL;
   freqModelhCross=XLALCreateCOMPLEX16FrequencySeries("freqDatahC",&(tmpdata->timeData->epoch),0.0,tmpdata->freqData->deltaF,&lalDimensionlessUnit,tmpdata->freqData->data->length);
    COMPLEX16FrequencySeries *freqModelhPlus=NULL;
    freqModelhPlus=XLALCreateCOMPLEX16FrequencySeries("freqDatahP",&(tmpdata->timeData->epoch),0.0,tmpdata->freqData->deltaF,&lalDimensionlessUnit,tmpdata->freqData->data->length);
    tmpdata->freqModelhPlus=freqModelhPlus;
    tmpdata->freqModelhCross=freqModelhCross;
    if(LALInferenceGetProcParamVal(commandLine,"--lalinspiralinjection")){
      LALInferenceTemplateLAL(tmpdata);
    }else{
      tmpdata->modelDomain = LAL_SIM_DOMAIN_FREQUENCY;
      LALInferenceTemplateXLALSimInspiralChooseWaveform(tmpdata);
    }

     
    LALInferenceVariables *currentParams=IFOdata->modelParams;
       
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  int i, lower, upper;
  LALInferenceIFOData *dataPtr;
  double ra, dec, psi, distMpc, gmst;
  LIGOTimeGPS GPSlal;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift;  /* time shift (not necessarily same as above)                   */
  double deltaT, deltaF, twopit, f, re, im;
 
  REAL8 temp=0.0;
	UINT4 logDistFlag=0;
    REAL8 NetSNR=0.0;
  LALInferenceVariables intrinsicParams;

  logDistFlag=LALInferenceCheckVariable(currentParams, "logdistance");

  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
	if(logDistFlag)
		 distMpc = exp(*(REAL8*)LALInferenceGetVariable(currentParams,"logdistance"));
	else
		 distMpc   = *(REAL8*) LALInferenceGetVariable(currentParams, "distance");       /* Mpc         */

  /* figure out GMST: */
  //XLALGPSSetREAL8(&GPSlal, GPSdouble); //This is what used in the likelihood. It seems off by two seconds (should not make a big difference as the antenna patterns would not change much in such a short interval)
  XLALGPSSetREAL8(&GPSlal, injtime);
  //UandA.units    = MST_RAD;
  //UandA.accuracy = LALLEAPSEC_LOOSE;
  //LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
  intrinsicParams.head      = NULL;
  intrinsicParams.dimension = 0;
  LALInferenceCopyVariables(currentParams, &intrinsicParams);
  LALInferenceRemoveVariable(&intrinsicParams, "rightascension");
  LALInferenceRemoveVariable(&intrinsicParams, "declination");
  LALInferenceRemoveVariable(&intrinsicParams, "polarisation");
  LALInferenceRemoveVariable(&intrinsicParams, "time");
	if(logDistFlag)
			LALInferenceRemoveVariable(&intrinsicParams, "logdistance");
	else
			LALInferenceRemoveVariable(&intrinsicParams, "distance");
  // TODO: add pointer to template function here.
  // (otherwise same parameters but different template will lead to no re-computation!!)

  
  /* loop over data (different interferometers): */
  dataPtr = IFOdata;
  
  while (dataPtr != NULL) {
     
      if (IFOdata->modelDomain == LAL_SIM_DOMAIN_TIME) {
	  printf("There is a problem. You seem to be using a time domain model into the frequency domain injection function!. Exiting....\n"); 
      exit(1);
    }
      
    /*-- WF to inject is now in dataPtr->freqModelhPlus and dataPtr->freqModelhCross. --*/
    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                             dataPtr->detector->response,
			     ra, dec, psi, gmst);
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */
    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (injtime - (*(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * (timeshift);
    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    dataPtr->fPlus = FplusScaled;
    dataPtr->fCross = FcrossScaled;
    dataPtr->timeshift = timeshift;

  //char InjFileName[50];
   //       sprintf(InjFileName,"injection_%s.dat",dataPtr->name);
   //       FILE *outInj=fopen(InjFileName,"w");
 
     /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
     chisquared = 0.0;
    for (i=lower; i<=upper; ++i){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      plainTemplateReal = FplusScaled * IFOdata->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * IFOdata->freqModelhCross->data->data[i].re;
      plainTemplateImag = FplusScaled * IFOdata->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * IFOdata->freqModelhCross->data->data[i].im;

      /* do time-shifting...             */
      /* (also un-do 1/deltaT scaling): */
      f = ((double) i) * deltaF;
      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      re = cos(twopit * f);
      im = - sin(twopit * f);
      templateReal = (plainTemplateReal*re - plainTemplateImag*im);
      templateImag = (plainTemplateReal*im + plainTemplateImag*re);
  
  
       //  fprintf(outInj,"%lf %e %e %e %e %e\n",i*deltaF ,dataPtr->freqData->data->data[i].re,dataPtr->freqData->data->data[i].im,templateReal,templateImag,1.0/dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      dataPtr->freqData->data->data[i].re+=templateReal;
      dataPtr->freqData->data->data[i].im+=templateImag;
   
      temp = ((2.0/( deltaT*(double) dataPtr->timeData->data->length) * (templateReal*templateReal+templateImag*templateImag)) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      chisquared  += temp;
    }
    printf("injected SNR %.1f in IFO %s\n",sqrt(2.0*chisquared),dataPtr->name);
    NetSNR+=2.0*chisquared;
    dataPtr->SNR=sqrt(2.0*chisquared);
    dataPtr = dataPtr->next;
    
// fclose(outInj);
  }

    LALInferenceClearVariables(&intrinsicParams);
    printf("injected Network SNR %.1f \n",sqrt(NetSNR)); 
    
    if (!(SNRpath==NULL)){ /* If the user provided a path with --snrpath store a file with injected SNRs */
	PrintSNRsToFile(IFOdata , inj_table);
	}
	XLALDestroyCOMPLEX16FrequencySeries(freqModelhCross);
    XLALDestroyCOMPLEX16FrequencySeries(freqModelhPlus);
}


static void PrintSNRsToFile(LALInferenceIFOData *IFOdata , SimInspiralTable *inj_table){
    char SnrName[200];
    char ListOfIFOs[10]="";
    REAL8 NetSNR=0.0;
    // sprintf(ListOfIFOs,"");
    LALInferenceIFOData *thisData=IFOdata;
    int nIFO=0;

    while(thisData){
         sprintf(ListOfIFOs,"%s%s",ListOfIFOs,thisData->name);
         thisData=thisData->next;
	nIFO++;
        }
    
    sprintf(SnrName,"%s/snr_%s_%10.1f.dat",SNRpath,ListOfIFOs,(REAL8) inj_table->geocent_end_time.gpsSeconds+ (REAL8) inj_table->geocent_end_time.gpsNanoSeconds*1.0e-9);
    FILE * snrout = fopen(SnrName,"w");
    if(!snrout){
	fprintf(stderr,"Unable to open the path %s for writing SNR files\n",SNRpath);
	exit(1);
    }
    
    thisData=IFOdata; // restart from the first IFO
    while(thisData){
        fprintf(snrout,"%s:\t %4.2f\n",thisData->name,thisData->SNR);
        NetSNR+=(thisData->SNR*thisData->SNR);
        thisData=thisData->next;
    }		
    if (nIFO>1){  fprintf(snrout,"Network:\t");
    fprintf(snrout,"%4.2f\n",sqrt(NetSNR));}
    fclose(snrout);
}

/** Fill the variables passed in vars with the parameters of the injection passed in event
    will over-write and destroy any existing parameters. Param vary type will be fixed */
void LALInferenceInjectionToVariables(SimInspiralTable *theEventTable, LALInferenceVariables *vars)
{
    UINT4 spinCheck=0;
    if(!vars) {
	XLALPrintError("Encountered NULL variables pointer");
   	XLAL_ERROR_VOID(XLAL_EINVAL);
	}
    /* Destroy existing parameters */
    if(vars->head!=NULL) LALInferenceClearVariables(vars);
    REAL8 q = theEventTable->mass2 / theEventTable->mass1;
    if (q > 1.0) q = 1.0/q;

    REAL8 sx = theEventTable->spin1x;
    REAL8 sy = theEventTable->spin1y;
    REAL8 sz = theEventTable->spin1z;

    REAL8 a_spin1 = sqrt(sx*sx + sy*sy + sz*sz);

    REAL8 theta_spin1, phi_spin1;
    if (a_spin1 == 0.0) {
      theta_spin1 = 0.0;
      phi_spin1 = 0.0;
    } else {
      theta_spin1 = acos(sz / a_spin1);
      phi_spin1 = atan2(sy, sx);
      if (phi_spin1 < 0.0) phi_spin1 += 2.0*M_PI;
    }

    sx = theEventTable->spin2x;
    sy = theEventTable->spin2y;
    sz = theEventTable->spin2z;

    REAL8 a_spin2 = sqrt(sx*sx + sy*sy + sz*sz), theta_spin2, phi_spin2;
    if (a_spin2 == 0.0) {
      theta_spin2 = 0.0;
      phi_spin2 = 0.0;
    } else {
      theta_spin2 = acos(sz / a_spin2);
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
    LALInferenceAddVariable(vars, "polarisation", &(psi), LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "phase", &phase, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "declination", &dec, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "rightascension", &ra, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "LAL_APPROXIMANT", &injapprox, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "LAL_PNORDER",&order,LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(vars, "LAL_AMPORDER", &(theEventTable->amp_order), LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    if(spinCheck){
        LALInferenceAddVariable(vars, "a_spin1", &a_spin1, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "a_spin2", &a_spin2, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "theta_spin1", &theta_spin1, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "theta_spin2", &theta_spin2, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "phi_spin1", &phi_spin1, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(vars, "phi_spin2", &phi_spin2, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
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
    LALPNOrder *order=LALInferenceGetVariable(&backup,"LAL_PNORDER");
    Approximant *approx=LALInferenceGetVariable(&backup,"LAL_APPROXIMANT");
    /* Fill named variables */
    LALInferenceInjectionToVariables(theEventTable,runState->currentParams);
    if(order && approx){
      /* Set the waveform to the one used in the analysis */
      LALInferenceRemoveVariable(runState->currentParams,"LAL_APPROXIMANT");
      LALInferenceRemoveVariable(runState->currentParams,"LAL_PNORDER");
      LALInferenceAddVariable(runState->currentParams,"LAL_PNORDER",order,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
      LALInferenceAddVariable(runState->currentParams,"LAL_APPROXIMANT",approx,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
    }
    REAL8 injPrior = runState->prior(runState,runState->currentParams);
    LALInferenceAddVariable(runState->currentParams,"logPrior",&injPrior,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    REAL8 injL = runState->likelihood(runState->currentParams, runState->data, runState->templt);
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

