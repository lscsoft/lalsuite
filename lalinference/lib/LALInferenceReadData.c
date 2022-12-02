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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>

#include <lal/LALInspiral.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALDetectors.h>
#include <lal/AVFactories.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/Sequence.h>
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
#include <math.h>
#include <lal/LALInspiral.h>
#include <lal/LALSimulation.h>

#include <lal/LALInference.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceInit.h>
#include <lal/LALSimNoise.h>
#include <LALInferenceRemoveLines.h>
/* LIB deps */
#include <lal/LALInferenceBurstRoutines.h>
#include <lal/LIGOLwXMLBurstRead.h>

struct fvec {
  REAL8 f;
  REAL8 x;
};

#define LALINFERENCE_DEFAULT_FLOW "20.0"

static void LALInferenceSetGPSTrigtime(LIGOTimeGPS *GPStrig, ProcessParamsTable *commandLine);
struct fvec *interpFromFile(char *filename, REAL8 squareinput);

struct fvec *interpFromFile(char *filename, REAL8 squareinput){
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
		interp[i].f=f;
		if (squareinput) {
			interp[i].x=x*x;
		}
		else {
			interp[i].x=x;
		}
		i++;
		if(i>fileLength-1){ /* Grow the array */
			interp=XLALRealloc(interp,(fileLength+minLength)*sizeof(struct fvec));
			fileLength+=minLength;
		}
	}
	interp[i].f=0.0; interp[i].x=0.0;
	fileLength=i+1;
	interp=XLALRealloc(interp,fileLength*sizeof(struct fvec)); /* Resize array */
	fclose(interpfile);
	printf("Read %i records from %s\n",fileLength-1,filename);
    if(fileLength-1 <1)
    {
            XLALPrintError("Error: read no records from %s\n",filename);
            exit(1);
    }
    return interp;
}

REAL8 interpolate(struct fvec *fvec, REAL8 f);
REAL8 interpolate(struct fvec *fvec, REAL8 f){
	int i=0;
	REAL8 a=0.0; /* fractional distance between bins */
	if(f<fvec[1].f) return(INFINITY); /* Frequency below minimum */
	while((i==0) || (fvec[i].f<f && (fvec[i].x!=0.0 ))){i++;}; //&& fvec[i].f!=0.0)){i++;};
	if (fvec[i].f==0.0 && fvec[i].x==0.0) /* Frequency above maximum */
	{
		return (INFINITY);
	}
	a=(fvec[i].f-f)/(fvec[i].f-fvec[i-1].f);
	return (fvec[i-1].x*a + fvec[i].x*(1.0-a));
}

void InjectFD(LALInferenceIFOData *IFOdata, SimInspiralTable *inj_table, ProcessParamsTable *commandLine);
int enforce_m1_larger_m2(SimInspiralTable* injEvent);

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


static const LALUnit strainPerCount={0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

static REAL8TimeSeries *readTseries(LALCache *cache, CHAR *channel, LIGOTimeGPS start, REAL8 length);
static void makeWhiteData(LALInferenceIFOData *IFOdata);
static void PrintSNRsToFile(LALInferenceIFOData *IFOdata , char SNRpath[] );


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
                        XLAL_NUM_ELEM(ifoRegExPattern), ".*%c.*",
                        ifo[0] );
        fprintf(stderr,"GlobFramesPWD : Found unseived src files:\n");
        for(UINT4 i=0;i<frGlobCache->length;i++)
                        fprintf(stderr,"(%s,%s,%s)\n",frGlobCache->list[i].src,frGlobCache->list[i].dsc,frGlobCache->list[i].url);
        frInCache = XLALCacheDuplicate(frGlobCache);
        XLALCacheSieve(frInCache, 0, 0, ifoRegExPattern, NULL, NULL);
        if ( ! frGlobCache->length )
        {
            fprintf( stderr, "error: no frame file files found after sieving\n");
            exit( 1 );
        }
        else
        {
                fprintf(stderr,"GlobFramesPWD : Sieved frames with pattern %s. Found src files:\n",ifoRegExPattern);
                for(UINT4 i=0;i<frInCache->length;i++)
                        fprintf(stderr,"(%s,%s,%s)\n",frInCache->list[i].src,frInCache->list[i].dsc,frInCache->list[i].url);
        }

        return(frGlobCache);
}

static REAL8TimeSeries *readTseries(LALCache *cache, CHAR *channel, LIGOTimeGPS start, REAL8 length)
{
    /* This function attempts to read the data from the given cache file. If it failes to open the data
	 * it waits for a period and tries again, up to 5 times. This is an attempt to get around overloading
	 * of file servers when many jobs are run at the time same */
	LALStatus status;
	UINT4 max_tries=7,tries=0,delay=5;
	memset(&status,0,sizeof(status));
	LALFrStream *stream = NULL;
	REAL8TimeSeries *out = NULL;
	if(cache==NULL) fprintf(stderr,"readTseries ERROR: Received NULL pointer for channel %s\n",channel);
	for (tries=0,delay=5;tries<max_tries;tries++,delay*=2) {
			stream = XLALFrStreamCacheOpen( cache );
			if(stream) break;
			sleep(delay);
	}
	if(stream==NULL) {fprintf(stderr,"readTseries ERROR: Unable to open stream from frame cache file\n"); exit(-1);}

	for (tries=0,delay=5;tries<max_tries;tries++,delay*=2) {
			out = XLALFrStreamInputREAL8TimeSeries( stream, channel, &start, length , 0 );
			if(out) break;
			sleep(delay);
	}
	if(out==NULL) fprintf(stderr,"readTseries ERROR: unable to read channel %s at times %i - %f\nCheck the specified data duration is not too long\n",channel,start.gpsSeconds,start.gpsSeconds+length);
	XLALFrStreamClose(stream);
	return out;
}

/**
 * Parse the command line looking for options of the kind --ifo H1 --H1-channel H1:LDAS_STRAIN --H1-cache H1.cache --H1-flow 20.0 --H1-fhigh 4096.0 --H1-timeslide 100.0 --H1-asd asd_ascii.txt --H1-psd psd_ascii.txt ...
 * It is necessary to use this method instead of the old method for the pipeline to work in DAX mode. Warning: do not mix options between
 * the old and new style.
 */
static INT4 getDataOptionsByDetectors(ProcessParamsTable *commandLine, char ***ifos, char ***caches, char ***channels, char ***flows , char ***fhighs, char ***timeslides, char ***asds, char ***psds, UINT4 *N)
{
    /* Check that the input has no lists with [ifo,ifo] */
    ProcessParamsTable *this=commandLine;
    UINT4 i=0;
    *caches=*ifos=*channels=*flows=*fhighs=*timeslides=*asds=*psds=NULL;
    *N=0;
    char tmp[128];
    if(!this) {fprintf(stderr,"No command line arguments given!\n"); exit(1);}
    /* Construct a list of IFOs */
    for(this=commandLine;this;this=this->next)
    {
        if(!strcmp(this->param,"--ifo"))
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
    *asds=XLALCalloc(*N,sizeof(char *));
    *psds=XLALCalloc(*N,sizeof(char *));

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

        /* ASD */
        sprintf(tmp,"--%s-asd",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*asds)[i]=this?XLALStringDuplicate(this->value):NULL;

        /* PSD */
        sprintf(tmp,"--%s-psd",(*ifos)[i]);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*psds)[i]=this?XLALStringDuplicate(this->value):NULL;
    }
    return(1);
}

/**
 * Parse the command line looking for options of the kind ---IFO-name value
 * Unlike the function above, this one does not have a preset list of names to lookup, but instead uses the option "name"
 * It is necessary to use this method instead of the old method for the pipeline to work in DAX mode. Warning: do not mix options between
 * the old and new style.
 * Return 0 if the number of options --IFO-name doesn't much the number of ifos, 1 otherwise. Fills in the pointer out with the values that were found.
 */
static INT4 getNamedDataOptionsByDetectors(ProcessParamsTable *commandLine, char ***ifos, char ***out, const char *name, UINT4 *N)
{
    /* Check that the input has no lists with [ifo,ifo] */
    ProcessParamsTable *this=commandLine;
    UINT4 i=0;
    *out=*ifos=NULL;
    *N=0;
    char tmp[128];
    if(!this) {fprintf(stderr,"No command line arguments given!\n"); exit(1);}
    /* Construct a list of IFOs */
    for(this=commandLine;this;this=this->next)
    {
        if(!strcmp(this->param,"--ifo"))
        {
            (*N)++;
            *ifos=XLALRealloc(*ifos,*N*sizeof(char *));
            (*ifos)[*N-1]=XLALStringDuplicate(this->value);
        }
    }
    *out=XLALCalloc(*N,sizeof(char *));

    UINT4 found=0;
    /* For each IFO, fetch the other options if available */
    for(i=0;i<*N;i++)
    {
        /* Channel */
        sprintf(tmp,"--%s-%s",(*ifos)[i],name);
        this=LALInferenceGetProcParamVal(commandLine,tmp);
        (*out)[i]=this?XLALStringDuplicate(this->value):NULL;
	if (this) found++;

    }
    if (found==*N)
      return(1);
    else
      return 0;
}

static void LALInferencePrintDataWithInjection(LALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine){

  LALStatus status;
  memset(&status,0,sizeof(status));
  UINT4 Nifo=0,i,j;
  LALInferenceIFOData *thisData=IFOdata;
  UINT4 q=0;
  UINT4 event=0;
  ProcessParamsTable *procparam=NULL,*ppt=NULL;
  SimInspiralTable *injTable=NULL;
  //ProcessParamsTable *pptdatadump=NULL;
  LIGOTimeGPS GPStrig;
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
        if(injTable->simulation_id == atol(procparam->value)) break;
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
    XLALStrToGPS(&GPStrig,procparam->value,NULL);
  }
  else{
    if(injTable) memcpy(&GPStrig,&(injTable->geocent_end_time),sizeof(GPStrig));
    else {
      fprintf(stderr,"+++ Error: No trigger time specifed and no injection given \n");
      exit(1);
    }
  }

  if (LALInferenceGetProcParamVal(commandLine, "--data-dump")) {
    //pptdatadump=LALInferenceGetProcParamVal(commandLine,"--data-dump");
    const UINT4 nameLength=FILENAME_MAX+50;
    char filename[nameLength];
    FILE *out;

    for (i=0;i<Nifo;i++) {

      ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");
      if(ppt) {
        if((int)nameLength<=snprintf(filename, nameLength, "%s%s-timeDataWithInjection.dat", ppt->value, IFOdata[i].name))
            XLAL_ERROR_VOID(XLAL_EINVAL, "Output filename too long!");
      }
      //else if(strcmp(pptdatadump->value,"")) {
      //  snprintf(filename, nameLength, "%s/%s-timeDataWithInjection.dat", pptdatadump->value, IFOdata[i].name);
      //}
      else
        snprintf(filename, nameLength, "%.3f_%s-timeDataWithInjection.dat", GPStrig.gpsSeconds+1e-9*GPStrig.gpsNanoSeconds, IFOdata[i].name);
      out = fopen(filename, "w");
      if(!out){
        fprintf(stderr,"Unable to open the path %s for writing time data with injection files\n",filename);
        exit(1);
      }
      for (j = 0; j < IFOdata[i].timeData->data->length; j++) {
        REAL8 t = XLALGPSGetREAL8(&(IFOdata[i].timeData->epoch)) +
        j * IFOdata[i].timeData->deltaT;
        REAL8 d = IFOdata[i].timeData->data->data[j];

        fprintf(out, "%.6f %g\n", t, d);
      }
      fclose(out);

      ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");
      if(ppt) {
        snprintf(filename, nameLength, "%s%s-freqDataWithInjection.dat", ppt->value, IFOdata[i].name)        ;
      }
      //else if(strcmp(pptdatadump->value,"")) {
      //  snprintf(filename, nameLength, "%s/%s-freqDataWithInjection.dat", pptdatadump->value, IFOdata[i].name);
      //}
      else
        snprintf(filename, nameLength, "%.3f_%s-freqDataWithInjection.dat", GPStrig.gpsSeconds+1e-9*GPStrig.gpsNanoSeconds, IFOdata[i].name);
      out = fopen(filename, "w");
      if(!out){
        fprintf(stderr,"Unable to open the path %s for writing freq data with injection files\n",filename);
        exit(1);
      }
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
    ----------------------------------------------\n\
    --- Data Parameters --------------------------\n\
    ----------------------------------------------\n\
    Options for reading/generating data. User should specify which interferometers\n\
    to use and their data source using the following options. The data source\n\
    can be either a LFS cache file (generated by gw_data_find) with channel name\n\
    (e.g. --H1-cache H1.cache --H1-channel H1:DCS-CALIB-STRAIN_C02 )\n\
    or internally-generated gaussian noise with a given detector PSD model\n\
    (e.g. --H1-cache LALSimAdLIGO --dataseed 1234)\n\
    or by searching for frame files in the local directory\n\
    (e.g. --glob-frame-data --H1-channel H1:DCS-CALIB-STRAIN_C02)\n\
    \n\
    Additional noise curves for simulated data can be specified by providing\n\
    their PSD or ASD as a text file. See detailed options below.\n\
    \n\
    --ifo IFO1 [--ifo IFO2 ...] IFOs can be H1,L1,V1\n\
    --IFO1-cache cache1         Cache files \n\
    [--IFO2-cache2 cache2 ...]      lal PSDs: LAL{Ad}LIGO, LALVirgo\n\
                                    lalsimuation PSDs: LALSim{Ad}LIGO, LALSim{Ad}Virgo\n\
                                    interpolate from file: interp:asd_file.txt\n\
    --psdstart GPStime          GPS start time of PSD estimation data\n\
    --psdlength length          Length of PSD estimation data in seconds\n\
    --seglen length             Length of segments for PSD estimation and analysis in seconds\n\
    (--glob-frame-data)         Will search for frame files containing data in the PWD.\n\
     				Filenames must begin with the IFO's 1-letter code, e.g. H-*.gwf\n\
    (--dont-dump-extras)        If given, won't save PSD and SNR files\n\
    (--dump-geocenter-pols)     If given, print out the TD/FD h_plus and h_cross polarisations\n\
    (--trigtime GPStime)        GPS time of the trigger to analyse\n\
                                    (optional when using --margtime or --margtimephi)\n\
    (--segment-start)           GPS time of the start of the segment\n\
                                     (optional with --trigtime,\n\
                                      default: seglen-2 s before --trigtime)\n\
    (--srate rate)              Downsample data to rate in Hz (4096.0,)\n\
    (--padding PAD [sec]        Override default 0.4 seconds padding\n\
    (--injectionsrate rate)     Downsample injection signal to rate in Hz (--srate)\n\
    (--IFO1-flow freq1          Specify lower frequency cutoff for overlap integral (20.0)\n\
     [--IFO2-flow freq2 ...])\n\
    (--IFO1-fhigh freq1         Specify higher frequency cutoff for overlap integral (Nyquist\n\
     [--IFO2-fhigh freq2 ...])      freq 0.5*srate)\n\
    (--IFO1-channel chan1       Specify channel names when reading cache files\n\
     [--IFO2-channel chan2 ...])\n\
         (--IFO1-asd asd1-ascii.txt        Read in ASD from ascii file. This is not equivalent \n\
     [--IFO2-asd asd2-ascii.txt ...])     to using --IFO1-cache interp:asd_file.txt since the former\n\
                                          won't use the ascii ASD to generate fake noise. \n\
    (--IFO1-psd psd1-ascii.txt        Read in PSD from ascii file. This is not equivalent \n\
     [--IFO2-psd psd2-ascii.txt ...])     to using --IFO1-cache interp:asd_file.txt since the former\n\
                                          won't use the ascii PSD to generate fake noise. \n\
    (--dataseed number)         Specify random seed to use when generating data\n\
    (--lalinspiralinjection)    Enables injections via the LALInspiral package\n\
    (--inj-fref)                Reference frequency of parameters in injection XML (default 100Hz)\n\
    (--inj-lambda1)             value of lambda1 to be injected, LALSimulation only (0)\n\
    (--inj-lambda2)             value of lambda2 to be injected, LALSimulation only (0)\n\
    (--inj-lambdaT              value of lambdaT to be injected (0)\n\
    (--inj-dlambdaT             value of dlambdaT to be injected (0)\n\
    (--inj-logp1)               value of logp1 to be injected\n\
    (--inj-gamma1)              value of gamma1 to be injected\n\
    (--inj-gamma2)              value of gamma2 to be injected\n\
    (--inj-gamma3)              value of gamma3 to be injected\n\
    (--inj-SDgamma0)            value of SDgamma0 to be injected (0)\n\
    (--inj-SDgamma1)            value of SDgamma1 to be injected (0)\n\
    (--inj-SDgamma2)            value of SDgamma2 to be injected (0)\n\
    (--inj-SDgamma3)            value of SDgamma3 to be injected (0)\n\
    (--inj-spinOrder PNorder)   Specify twice the injection PN order (e.g. 5 <==> 2.5PN)\n\
                                    of spin effects effects to use, only for LALSimulation\n\
                                    (default: -1 <==> Use all spin effects).\n\
    (--inj-tidalOrder PNorder)  Specify twice the injection PN order (e.g. 10 <==> 5PN)\n\
                                    of tidal effects to use, only for LALSimulation\n\
                                    (default: -1 <==> Use all tidal effects).\n\
    (--inj-spin-frame FRAME     Specify injection spin frame: choice of total-j, orbital-l, view.\n\
                                    (Default = OrbitalL).\n\
    (--inj-numreldata FileName) Location of NR data file for the injection of NR waveforms (with NR_hdf5 in injection XML file).\n\
    (--0noise)                  Sets the noise realisation to be identically zero\n\
                                    (for the fake caches above only)\n\
    \n"

LALInferenceIFOData *LALInferenceReadData(ProcessParamsTable *commandLine)
/* Read in the data and store it in a LALInferenceIFOData structure */
{
    LALStatus status;
    INT4 dataseed=0;
    memset(&status,0,sizeof(status));
    ProcessParamsTable *procparam=NULL,*ppt=NULL;
    //ProcessParamsTable *pptdatadump=NULL;
    LALInferenceIFOData *headIFO=NULL,*IFOdata=NULL;
    REAL8 SampleRate=4096.0,SegmentLength=0;
    if(LALInferenceGetProcParamVal(commandLine,"--srate")) SampleRate=atof(LALInferenceGetProcParamVal(commandLine,"--srate")->value);
    REAL8 defaultFLow = atof(LALINFERENCE_DEFAULT_FLOW);
    int nSegs=0;
    size_t seglen=0;
    REAL8TimeSeries *PSDtimeSeries=NULL;
    REAL8 padding=0.4;//Default was 1.0 second. However for The Event the Common Inputs specify a Tukey parameter of 0.1, so 0.4 second of padding for 8 seconds of data.
    UINT4 Ncache=0,Nifo=0,Nchannel=0,NfLow=0,NfHigh=0;
    UINT4 i,j;
    //int FakeFlag=0; - set but not used
    char strainname[]="LSC-STRAIN";
    //typedef void (NoiseFunc)(LALStatus *statusPtr,REAL8 *psd,REAL8 f);
    NoiseFunc *PSD=NULL;
    REAL8 scalefactor=1;
    RandomParams *datarandparam;
    int globFrames=0; // 0 = no, 1 = will search for frames in PWD
    char **channels=NULL;
    char **caches=NULL;
    char **asds=NULL;
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
    REAL8 asdFlag=0;

    if(LALInferenceGetProcParamVal(commandLine,"--glob-frame-data")) globFrames=1;

    /* Check if the new style command line arguments are used */
    INT4 dataOpts=getDataOptionsByDetectors(commandLine, &IFOnames, &caches, &channels, &fLows, &fHighs, &timeslides, &asds, &psds, &Nifo);
    /* Check for options if not given in the new style */
    if(!dataOpts){
        if(!(globFrames||LALInferenceGetProcParamVal(commandLine,"--cache"))||!(LALInferenceGetProcParamVal(commandLine,"--IFO")||LALInferenceGetProcParamVal(commandLine,"--ifo")))
            {fprintf(stderr,USAGE); return(NULL);}
        if(LALInferenceGetProcParamVal(commandLine,"--channel")){
            LALInferenceParseCharacterOptionString(LALInferenceGetProcParamVal(commandLine,"--channel")->value,&channels,&Nchannel);
        }
        LALInferenceParseCharacterOptionString(LALInferenceGetProcParamVal(commandLine,"--cache")->value,&caches,&Ncache);
        ppt=LALInferenceGetProcParamVal(commandLine,"--ifo");
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
	if(!LALInferenceGetProcParamVal(commandLine,"--seglen"))
    {fprintf(stderr,USAGE); return(NULL);}

    if(LALInferenceGetProcParamVal(commandLine,"--dataseed")){
        procparam=LALInferenceGetProcParamVal(commandLine,"--dataseed");
        dataseed=atoi(procparam->value);
    }

    IFOdata=headIFO=XLALCalloc(sizeof(LALInferenceIFOData),Nifo);
    if(!IFOdata) XLAL_ERROR_NULL(XLAL_ENOMEM);

    procparam=LALInferenceGetProcParamVal(commandLine,"--psdstart");
    if (procparam) {
        XLALStrToGPS(&GPSstart,procparam->value,NULL);
        if(status.statusCode) REPORTSTATUS(&status);
    } else
        XLALINT8NSToGPS(&GPSstart, 0);

    /*Set trigtime in GPStrig using either inj file or --trigtime*/
    LALInferenceSetGPSTrigtime(&GPStrig,commandLine);

    if(status.statusCode) REPORTSTATUS(&status);

    SegmentLength=atof(LALInferenceGetProcParamVal(commandLine,"--seglen")->value);
    seglen=(size_t)(SegmentLength*SampleRate);

    ppt=LALInferenceGetProcParamVal(commandLine,"--psdlength");
    if(ppt) {
        PSDdatalength=atof(ppt->value);
        nSegs=(int)floor(PSDdatalength/SegmentLength);
    }

    CHAR df_argument_name[262];
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

        if(!strcmp(IFOnames[i],"E1")){
            memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexE1DIFF],sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"E1:STRAIN"); continue;}
        if(!strcmp(IFOnames[i],"E2")){
            memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexE2DIFF],sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"E2:STRAIN"); continue;}
        if(!strcmp(IFOnames[i],"E3")){
            memcpy(IFOdata[i].detector,&lalCachedDetectors[LALDetectorIndexE3DIFF],sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"E3:STRAIN"); continue;}
        if(!strcmp(IFOnames[i],"K1")){
            memcpy(IFOdata[i].detector, &lalCachedDetectors[LALDetectorIndexKAGRADIFF],sizeof(LALDetector));
            if(!Nchannel) sprintf((channels[i]),"K1:STRAIN"); continue;}
	    if(!strcmp(IFOnames[i],"I1")){
	        memcpy(IFOdata[i].detector, &lalCachedDetectors[LALDetectorIndexLIODIFF],sizeof(LALDetector));
	        if(!Nchannel) sprintf((channels[i]),"I1:STRAIN"); continue;}
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
        fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 H2 L1 V1 GEO A1 K1 I1 E1 E2 E3 HM1 HM2 EM1 EM2\n",IFOnames[i]); exit(-1);
    }

    /* Set up FFT structures and window */
    for (i=0;i<Nifo;i++){
        /* Create FFT plans (flag 1 to measure performance) */
        IFOdata[i].timeToFreqFFTPlan = XLALCreateForwardREAL8FFTPlan((UINT4) seglen, 1 );
        if(!IFOdata[i].timeToFreqFFTPlan) XLAL_ERROR_NULL(XLAL_ENOMEM);
        IFOdata[i].freqToTimeFFTPlan = XLALCreateReverseREAL8FFTPlan((UINT4) seglen, 1 );
        if(!IFOdata[i].freqToTimeFFTPlan) XLAL_ERROR_NULL(XLAL_ENOMEM);
        IFOdata[i].margFFTPlan = XLALCreateReverseREAL8FFTPlan((UINT4) seglen, 1);
        if(!IFOdata[i].margFFTPlan) XLAL_ERROR_NULL(XLAL_ENOMEM);
        /* Setup windows */
        ppt=LALInferenceGetProcParamVal(commandLine,"--padding");
        if (ppt){
            padding=atof(ppt->value);
            fprintf(stdout,"Using %lf seconds of padding for IFO %s \n",padding, IFOdata[i].name);
        }
        if ((REAL8)2.0*padding*SampleRate/(REAL8)seglen <0.0 ||(REAL8)2.0*padding*SampleRate/(REAL8)seglen >1 ){
            fprintf(stderr,"Padding is negative or 2*padding is bigger than the whole segment. Consider reducing it using --padding or increase --seglen. Exiting\n");
            exit(1);
        }
        IFOdata[i].padding=padding;
        IFOdata[i].window=XLALCreateTukeyREAL8Window(seglen,(REAL8)2.0*padding*SampleRate/(REAL8)seglen);
        if(!IFOdata[i].window) XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    if(!(ppt=LALInferenceGetProcParamVal(commandLine,"--segment-start")))
    {
        /* Trigger time = 2 seconds before end of segment (was 1 second, but Common Inputs for The Events are -6 +2*/
        memcpy(&segStart,&GPStrig,sizeof(LIGOTimeGPS));
        REAL8 offset=SegmentLength-2.;
        /* If we are using a burst approximant, put at the center */
        if ((ppt=LALInferenceGetProcParamVal(commandLine,"--approx"))){
          if (XLALCheckBurstApproximantFromString(ppt->value)) offset=SegmentLength/2.;
        }
        XLALGPSAdd(&segStart,-offset);
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
         printf("Looking for ASD interpolation file %s\n",interpfilename);
         interpFlag=1;
         asdFlag=1;
         interp=interpFromFile(interpfilename, asdFlag);
        }
        /* Check if fake data is requested */
       if( (globFrames)?0:(interpFlag || (!(strcmp(caches[i],"LALLIGO") && strcmp(caches[i],"LALVirgo") && strcmp(caches[i],"LALGEO") && strcmp(caches[i],"LALEGO") && strcmp(caches[i],"LALSimLIGO") && strcmp(caches[i],"LALSimAdLIGO") && strcmp(caches[i],"LALSimVirgo") && strcmp(caches[i],"LALSimAdVirgo") && strcmp(caches[i],"LALAdLIGO")))))
        {
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
            if(!strcmp(caches[i],"LALSimAdLIGO")) {XLALSimNoisePSDaLIGODesignSensitivityT1800044(IFOdata[i].oneSidedNoisePowerSpectrum,IFOdata[i].fLow) ;LALSimPsd=1;}
            if(!strcmp(caches[i],"LALSimAdVirgo")) {XLALSimNoisePSD(IFOdata[i].oneSidedNoisePowerSpectrum,IFOdata[i].fLow,XLALSimNoisePSDAdvVirgo) ;LALSimPsd=1;}
            if(interpFlag) {PSD=NULL; scalefactor=1.0;}
            if(PSD==NULL && !(interpFlag|| LALSimPsd)) {fprintf(stderr,"Error: unknown simulated PSD: %s\n",caches[i]); exit(-1);}

            if(LALSimPsd==0){
                for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
                {
                    MetaNoiseFunc(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,interp,PSD);
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
            if ((!((psds[i])==NULL)) && (!((asds[i])==NULL))) {fprintf(stderr,"ERROR: Cannot provide both ASD and PSD file from command line!\n"); exit(1);}
            if (!((asds)==NULL || (asds[i])==NULL)){
                interp=NULL;
                asdFlag=1;
                char *interpfilename=&(asds[i][0]);
                fprintf(stderr,"Reading ASD for %s using %s\n",IFOnames[i],interpfilename);
                printf("Looking for ASD file %s for PSD interpolation\n",interpfilename);
                interp=interpFromFile(interpfilename, asdFlag);
                IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)
                    XLALCreateREAL8FrequencySeries("spectrum",&GPSstart,0.0,
                            (REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
                if(!IFOdata[i].oneSidedNoisePowerSpectrum) XLAL_ERROR_NULL(XLAL_EFUNC);
                for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
                {
                    MetaNoiseFunc(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,interp,NULL);
                    //fprintf(stdout,"%lf\n",IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]);
                }
             }else if(!((psds)==NULL || (psds[i])==NULL)){
                interp=NULL;
                asdFlag=0;
                char *interpfilename=&(psds[i][0]);
                fprintf(stderr,"Reading PSD for %s using %s\n",IFOnames[i],interpfilename);
                printf("Looking for PSD file %s for PSD interpolation\n",interpfilename);
                interp=interpFromFile(interpfilename, asdFlag);
                IFOdata[i].oneSidedNoisePowerSpectrum=(REAL8FrequencySeries *)
                    XLALCreateREAL8FrequencySeries("spectrum",&GPSstart,0.0,
                            (REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
                if(!IFOdata[i].oneSidedNoisePowerSpectrum) XLAL_ERROR_NULL(XLAL_EFUNC);
                for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
                {
                    MetaNoiseFunc(&status,&(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]),j*IFOdata[i].oneSidedNoisePowerSpectrum->deltaF,interp,NULL);
                    //fprintf(stdout,"%lf\n",IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]);
                }
            }else{
                /* Make sure required PSD arguments were provided */
                if(!LALInferenceGetProcParamVal(commandLine,"--psdstart") ||
                        !LALInferenceGetProcParamVal(commandLine,"--psdlength"))
                {fprintf(stderr,USAGE); return(NULL);}

                fprintf(stderr,"Estimating PSD for %s using %i segments of %i samples (%lfs)\n",IFOnames[i],nSegs,(int)seglen,SegmentLength);
                LIGOTimeGPS trueGPSstart=GPSstart;
                if(Ntimeslides) {
                  REAL4 deltaT=-atof(timeslides[i]);
                  XLALGPSAdd(&GPSstart, deltaT);
                  fprintf(stderr,"Slid PSD estimation of %s by %f s from %10.10lf to %10.10lf\n",IFOnames[i],deltaT,trueGPSstart.gpsSeconds+1e-9*trueGPSstart.gpsNanoSeconds,GPSstart.gpsSeconds+1e-9*GPSstart.gpsNanoSeconds);
                }
                PSDtimeSeries=readTseries(cache,channels[i],GPSstart,PSDdatalength);
                GPSstart=trueGPSstart;
                if(!PSDtimeSeries) {XLALPrintError("Error reading PSD data for %s\n",IFOnames[i]); exit(1);}
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

        makeWhiteData(&(IFOdata[i]));

      /* Store ASD of noise spectrum to whiten glitch model */
      IFOdata[i].noiseASD=(REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("asd",&GPSstart,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
      for(j=0;j<IFOdata[i].oneSidedNoisePowerSpectrum->data->length;j++)
        IFOdata[i].noiseASD->data->data[j]=sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]);
        /* Save to file the PSDs so that they can be used in the PP pages */
        const UINT4 nameLength=FILENAME_MAX+100;
        char filename[nameLength];
        FILE *out;
        ppt=LALInferenceGetProcParamVal(commandLine,"--dont-dump-extras");
        if (!ppt){
          ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");
          if(ppt) {
            snprintf(filename, nameLength, "%s%s-PSD.dat", ppt->value, IFOdata[i].name);
          }
          else
            snprintf(filename, nameLength, "%.3f_%s-PSD.dat",GPStrig.gpsSeconds+1e-9*GPStrig.gpsNanoSeconds, IFOdata[i].name);

          out = fopen(filename, "w");
          if(!out){
            fprintf(stderr,"Unable to open the path %s for writing PSD files\n",filename);
            exit(1);
          }
          for (j = 0; j < IFOdata[i].oneSidedNoisePowerSpectrum->data->length; j++) {
              REAL8 f = IFOdata[i].oneSidedNoisePowerSpectrum->deltaF*j;
              REAL8 psd = IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j];

              fprintf(out, "%10.10g %10.10g\n", f, psd);
          }
          fclose(out);
        }
        if (LALInferenceGetProcParamVal(commandLine, "--data-dump")) {
            //pptdatadump=LALInferenceGetProcParamVal(commandLine,"--data-dump");

            ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");

            if(ppt) {
              snprintf(filename, nameLength, "%s%s-timeData.dat", ppt->value, IFOdata[i].name);
            }
            //else if(strcmp(pptdatadump->value,"")) {
            //  snprintf(filename, nameLength, "%s/%s-timeData.dat", pptdatadump->value, IFOdata[i].name);
            //}
            else
              snprintf(filename, nameLength, "%.3f_%s-timeData.dat",GPStrig.gpsSeconds+1e-9*GPStrig.gpsNanoSeconds, IFOdata[i].name);
            out = fopen(filename, "w");
            if(!out){
                fprintf(stderr,"Unable to open the path %s for writing time data files\n",filename);
                exit(1);
            }
            for (j = 0; j < IFOdata[i].timeData->data->length; j++) {
                REAL8 t = XLALGPSGetREAL8(&(IFOdata[i].timeData->epoch)) +
                    j * IFOdata[i].timeData->deltaT;
                REAL8 d = IFOdata[i].timeData->data->data[j];

                fprintf(out, "%.6f %g\n", t, d);
            }
            fclose(out);

            ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");
            if(ppt) {
              snprintf(filename, nameLength, "%s%s-freqData.dat", ppt->value, IFOdata[i].name);
            }
            //else if(strcmp(pptdatadump->value,"")) {
            //  snprintf(filename, nameLength, "%s/%s-freqData.dat", pptdatadump->value, IFOdata[i].name);
            //}
            else
              snprintf(filename, nameLength, "%.3f_%s-freqData.dat",GPStrig.gpsSeconds+1e-9*GPStrig.gpsNanoSeconds, IFOdata[i].name);
            out = fopen(filename, "w");
            if(!out){
                fprintf(stderr,"Unable to open the path %s for writing freq data files\n",filename);
                exit(1);
            }
            for (j = 0; j < IFOdata[i].freqData->data->length; j++) {
                REAL8 f = IFOdata[i].freqData->deltaF * j;
                REAL8 dre = creal(IFOdata[i].freqData->data->data[j]);
                REAL8 dim = cimag(IFOdata[i].freqData->data->data[j]);

                fprintf(out, "%10.10g %10.10g %10.10g\n", f, dre, dim);
            }
            fclose(out);

            ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");
            if(ppt) {
              snprintf(filename, nameLength, "%s%s-ASD.dat", ppt->value, IFOdata[i].name);
            }
            else
              snprintf(filename, nameLength, "%.3f_%s-ASD.dat",GPStrig.gpsSeconds+1e-9*GPStrig.gpsNanoSeconds, IFOdata[i].name);
            out = fopen(filename, "w");
            if(!out){
              fprintf(stderr,"Unable to open the path %s for writing freq ASD files\n",filename);
              exit(1);
            }
            for (j = 0; j < IFOdata[i].oneSidedNoisePowerSpectrum->data->length; j++) {
              REAL8 f = IFOdata[i].oneSidedNoisePowerSpectrum->deltaF*j;
              REAL8 asd = sqrt(IFOdata[i].oneSidedNoisePowerSpectrum->data->data[j]);

              fprintf(out, "%10.10g %10.10g\n", f, asd);
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

    if (LALInferenceGetProcParamVal(commandLine, "--roqtime_steps")){

        LALInferenceSetupROQdata(IFOdata, commandLine);
        fprintf(stderr, "done LALInferenceSetupROQdata\n");

     }


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
	LIGOTimeGPS injstart;
	REAL8 SNR=0,NetworkSNR=0;
	DetectorResponse det;
	memset(&injstart,0,sizeof(LIGOTimeGPS));
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
    char SNRpath[FILENAME_MAX+50]="";
    int flipped_masses=0;

	while(thisData){
          minFlow   = minFlow>thisData->fLow ? thisData->fLow : minFlow;
          MindeltaT = MindeltaT>thisData->timeData->deltaT ? thisData->timeData->deltaT : MindeltaT;
          thisData  = thisData->next;
	}
	thisData=IFOdata;

	if(!LALInferenceGetProcParamVal(commandLine,"--inj")) {fprintf(stdout,"No injection file specified, not injecting\n"); return;}
	if(LALInferenceGetProcParamVal(commandLine,"--event")){
    event= atoi(LALInferenceGetProcParamVal(commandLine,"--event")->value);
    fprintf(stdout,"Injecting event %d\n",event);
	}

	ppt = LALInferenceGetProcParamVal(commandLine,"--outfile");
	if (ppt)
	    sprintf(SNRpath, "%s_snr.txt", ppt->value);
	else
		sprintf(SNRpath, "snr.txt");

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
	Approximant injapprox;
	injapprox = XLALGetApproximantFromString(injTable->waveform);
        if( (int) injapprox == XLAL_FAILURE)
          ABORTXLAL(&status);
	printf("Injecting approximant %i: %s\n", injapprox, injTable->waveform);
	REPORTSTATUS(&status);

	/* Check for frequency domain injection. All aproximants supported by XLALSimInspiralImplementedFDApproximants will work.
   * CAVEAT: FD spinning approximants will refer the spin to the lower frequency as given in the xml table. Templates instead will refer it to the lower cutoff of the likelihood integral. This means *different* values of spin will be recovered if one doesn't pay attention! */
	if(XLALSimInspiralImplementedFDApproximants(XLALGetApproximantFromString(injEvent->waveform)))
	{
	 InjectFD(IFOdata, injTable, commandLine);
	 LALInferencePrintDataWithInjection(IFOdata,commandLine);
	 return;
	}
        flipped_masses = enforce_m1_larger_m2(injEvent);
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

    if(LALInferenceGetProcParamVal(commandLine,"--lalinspiralinjection")){
      if ( approximant == NumRelNinja2) {
        XLALSimInjectNinjaSignals(injectionBuffer, thisData->name, 1./responseScale, injEvent);
      } else {
        /* Use this custom version for extra sites - not currently maintained */
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

      // Inject (lambda1,lambda2)
      REAL8 lambda1 = 0.;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")) {
        lambda1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")->value);
      }
      REAL8 lambda2 = 0.;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")) {
        lambda2= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")->value);
      }
	  if(flipped_masses)
	  {
			  /* Having previously flipped the masses so m1>m2, also flip lambdas */
			  REAL8 lambda_tmp = lambda1;
			  lambda1 = lambda2;
			  lambda2 = lambda_tmp;
			  fprintf(stdout,"Flipping lambdas since masses are flipped\n");
	  }
      if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")) {
        fprintf(stdout,"Injection lambda1 set to %f\n",lambda1);
      }
      if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")) {
        fprintf(stdout,"Injection lambda1 set to %f\n",lambda1);
      }

      // Inject (lambdaT,dLambdaT)
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

      // Inject 4-piece polytrope eos
      REAL8 logp1=0.0;
      REAL8 gamma1=0.0;
      REAL8 gamma2=0.0;
      REAL8 gamma3=0.0;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-logp1")&&LALInferenceGetProcParamVal(commandLine,"--inj-gamma1")&&LALInferenceGetProcParamVal(commandLine,"--inj-gamma2")&&LALInferenceGetProcParamVal(commandLine,"--inj-gamma3")){
        logp1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-logp1")->value);
        gamma1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-gamma1")->value);
        gamma2= atof(LALInferenceGetProcParamVal(commandLine,"--inj-gamma2")->value);
        gamma3= atof(LALInferenceGetProcParamVal(commandLine,"--inj-gamma3")->value);
        // Find lambda1,2(m1,2|eos)
        LALInferenceLogp1GammasMasses2Lambdas(logp1, gamma1, gamma2, gamma3, m1, m2, &lambda1, &lambda2);
        fprintf(stdout,"Injection logp1 set to %f\n",logp1);
        fprintf(stdout,"Injection gamma1 set to %f\n",gamma1);
        fprintf(stdout,"Injection gamma2 set to %f\n",gamma2);
        fprintf(stdout,"Injection gamma3 set to %f\n",gamma3);
        fprintf(stdout,"lambda1 set to %f\n",lambda1);
        fprintf(stdout,"lambda2 set to %f\n",lambda2);

        /*
        For when tagSimInspiralTable is updated
        to include EOS params

        injEvent->logp1= logp1;
        injEvent->gamma1= gamma1;
        injEvent->gamma2= gamma2;
        injEvent->gamma3= gamma3;
        */
      }

      // Inject 4-coef. spectral eos
      REAL8 SDgamma0=0.0;
      REAL8 SDgamma1=0.0;
      REAL8 SDgamma2=0.0;
      REAL8 SDgamma3=0.0;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma0") && LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma1") && LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma2") && LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma3")){
        SDgamma0= atof(LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma0")->value);
        SDgamma1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma1")->value);
       SDgamma2= atof(LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma2")->value);
        SDgamma3= atof(LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma3")->value);
        REAL8 gamma[]={SDgamma0,SDgamma1,SDgamma2,SDgamma3};
        LALInferenceSDGammasMasses2Lambdas(gamma,m1,m2,&lambda1,&lambda2,4);
        fprintf(stdout,"Injection SDgamma0 set to %lf\n",SDgamma0);
        fprintf(stdout,"Injection SDgamma1 set to %lf\n",SDgamma1);
        fprintf(stdout,"Injection SDgamma2 set to %lf\n",SDgamma2);
        fprintf(stdout,"Injection SDgamma3 set to %lf\n",SDgamma3);
        fprintf(stdout,"lambda1 set to %f\n",lambda1);
        fprintf(stdout,"lambda2 set to %f\n",lambda2);
      }


      REAL8 fref = 100.;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-fref")) {
        fref = atoi(LALInferenceGetProcParamVal(commandLine,"--inj-fref")->value);
      }

      LALDict *LALpars=XLALCreateDict();

      /* Set the spin-frame convention */

      LALSimInspiralSpinOrder spinO = -1;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-spinOrder")) {
        spinO = atoi(LALInferenceGetProcParamVal(commandLine,"--inj-spinOrder")->value);
        XLALSimInspiralWaveformParamsInsertPNSpinOrder(LALpars, spinO);
      }
      LALSimInspiralTidalOrder tideO = -1;
      if(LALInferenceGetProcParamVal(commandLine,"--inj-tidalOrder")) {
        tideO = atoi(LALInferenceGetProcParamVal(commandLine,"--inj-tidalOrder")->value);
        XLALSimInspiralWaveformParamsInsertPNTidalOrder(LALpars, tideO);
      }

      LALSimInspiralFrameAxis frameAxis = LAL_SIM_INSPIRAL_FRAME_AXIS_DEFAULT;
      if((ppt=LALInferenceGetProcParamVal(commandLine,"--inj-spin-frame"))) {
        frameAxis = XLALSimInspiralGetFrameAxisFromString(ppt->value);
      }
      XLALSimInspiralWaveformParamsInsertFrameAxis(LALpars,(int) frameAxis);

      if((ppt=LALInferenceGetProcParamVal(commandLine,"--inj-numreldata"))) {
	XLALSimInspiralWaveformParamsInsertNumRelData(LALpars, ppt->value);
	fprintf(stdout,"Injection will use %s.\n",ppt->value);
      }
      else if (strlen(injEvent->numrel_data) > 0) {
        XLALSimInspiralWaveformParamsInsertNumRelData(LALpars, injEvent->numrel_data);
      }

      /* Print a line with information about approximant, amporder, phaseorder, tide order and spin order */
      fprintf(stdout,"Injection will run using Approximant %i (%s), phase order %i, amp order %i, spin order %i, tidal order %i, in the time domain with a reference frequency of %f.\n",approximant,XLALSimInspiralGetStringFromApproximant(approximant),order,amporder,(int) spinO, (int) tideO, (float) fref);

      /* ChooseWaveform starts the (2,2) mode of the waveform at the given minimum frequency.  We want the highest order contribution to start at the f_lower of the injection file */
      REAL8 f_min = XLALSimInspiralfLow2fStart(injEvent->f_lower, amporder, approximant);
      printf("Injecting with f_min = %f.\n", f_min);

      XLALSimInspiralWaveformParamsInsertTidalLambda1(LALpars,lambda1);
      XLALSimInspiralWaveformParamsInsertTidalLambda2(LALpars,lambda2);
      XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(LALpars,amporder);
      XLALSimInspiralWaveformParamsInsertPNPhaseOrder(LALpars,order);

      XLALSimInspiralChooseTDWaveform(&hplus, &hcross, injEvent->mass1*LAL_MSUN_SI, injEvent->mass2*LAL_MSUN_SI,
				      injEvent->spin1x, injEvent->spin1y, injEvent->spin1z,
				      injEvent->spin2x, injEvent->spin2y, injEvent->spin2z,
				      injEvent->distance*LAL_PC_SI * 1.0e6, injEvent->inclination,
				      injEvent->coa_phase, 0., 0., 0.,
				      1.0/InjSampleRate, f_min, fref,
				      LALpars,  approximant);
      XLALDestroyDict(LALpars);

      if(!hplus || !hcross) {
        fprintf(stderr,"Error: XLALSimInspiralChooseWaveform() failed to produce waveform.\n");
        exit(-1);
      }
      XLALResampleREAL8TimeSeries(hplus,thisData->timeData->deltaT);
      XLALResampleREAL8TimeSeries(hcross,thisData->timeData->deltaT);
      /* XLALSimInspiralChooseTDWaveform always ends the waveform at t=0 */
      /* So we can adjust the epoch so that the end time is as desired */
      XLALGPSAddGPS(&(hplus->epoch), &(injEvent->geocent_end_time));
      XLALGPSAddGPS(&(hcross->epoch), &(injEvent->geocent_end_time));

      if((ppt=LALInferenceGetProcParamVal(commandLine,"--dump-geocenter-pols"))) {

        fprintf(stdout,"Dump injected TimeDomain h_plus and h_cross at geocenter (for IFO %s)\n", thisData->name);

        char filename[320];
        sprintf(filename,"%s_TD_geocenter_pols.dat",thisData->name);

        FILE* file=fopen(filename, "w");
        if(!file){
          fprintf(stderr,"Unable to open the path %s for writing injected TimeDomain h_plus and h_cross at geocenter\n",filename);
          exit(1);
        }

        for(j=0; j<hplus->data->length; j++){
          fprintf(file, "%.6f %g %g \n", XLALGPSGetREAL8(&hplus->epoch) + hplus->deltaT*j, hplus->data->data[j], hcross->data->data[j]);
        }
        fclose(file);
      }

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
      char filename[320];
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
    ppt=LALInferenceGetProcParamVal(commandLine,"--dont-dump-extras");
    if (!ppt){
      PrintSNRsToFile(IFOdata , SNRpath);
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


void InjectFD(LALInferenceIFOData *IFOdata, SimInspiralTable *inj_table, ProcessParamsTable *commandLine)
///*-------------- Inject in Frequency domain -----------------*/
{
  /* Inject a gravitational wave into the data in the frequency domain */
  LALStatus status;
  memset(&status,0,sizeof(LALStatus));
  INT4 errnum;
  char SNRpath[FILENAME_MAX+50];
  ProcessParamsTable *ppt=NULL;
  int flipped_masses=0;
  LALInferenceIFOData *dataPtr;

  ppt = LALInferenceGetProcParamVal(commandLine,"--outfile");
  if (ppt)
    sprintf(SNRpath, "%s_snr.txt", ppt->value);
  else
    sprintf(SNRpath, "snr.txt");

  Approximant approximant = XLALGetApproximantFromString(inj_table->waveform);
  if( (int) approximant == XLAL_FAILURE)
      ABORTXLAL(&status);

  LALPNOrder phase_order = XLALGetOrderFromString(inj_table->waveform);
  if ( (int) phase_order == XLAL_FAILURE)
      ABORTXLAL(&status);

  LALPNOrder amp_order = (LALPNOrder) inj_table->amp_order;

  flipped_masses = enforce_m1_larger_m2(inj_table);

  REAL8 injtime=0.0;
  injtime=(REAL8) inj_table->geocent_end_time.gpsSeconds + (REAL8) inj_table->geocent_end_time.gpsNanoSeconds*1.0e-9;

  // Inject (lambda1,lambda2)
  REAL8 lambda1 = 0.;
  if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")) {
    lambda1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")->value);
  }
  REAL8 lambda2 = 0.;
  if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")) {
    lambda2= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")->value);
  }
  if(flipped_masses) /* If flipped masses, also flip lambda */
  {
		REAL8 lambda_tmp=lambda1;
		lambda1=lambda2;
		lambda2=lambda_tmp;
		fprintf(stdout,"Flipping lambdas since masses are flipped\n");
  }
  if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda1")) {
    fprintf(stdout,"Injection lambda1 set to %f\n",lambda1);
  }
  if(LALInferenceGetProcParamVal(commandLine,"--inj-lambda2")) {
    fprintf(stdout,"Injection lambda2 set to %f\n",lambda2);
  }


  // Inject (lambdaT,dLambdaT)
  REAL8 lambdaT = 0.;
  REAL8 dLambdaT = 0.;
  if(LALInferenceGetProcParamVal(commandLine,"--inj-lambdaT")&&LALInferenceGetProcParamVal(commandLine,"--inj-dLambdaT")) {
    lambdaT= atof(LALInferenceGetProcParamVal(commandLine,"--inj-lambdaT")->value);
    dLambdaT= atof(LALInferenceGetProcParamVal(commandLine,"--inj-dLambdaT")->value);
    // Find lambda1,2(LambdaT,dLambdaT,eta)
    LALInferenceLambdaTsEta2Lambdas(lambdaT, dLambdaT, inj_table->eta, &lambda1, &lambda2);
    fprintf(stdout,"Injection lambdaT set to %f\n",lambdaT);
    fprintf(stdout,"Injection dLambdaT set to %f\n",dLambdaT);
    fprintf(stdout,"lambda1 set to %f\n",lambda1);
    fprintf(stdout,"lambda2 set to %f\n",lambda2);
  }

   // Inject 4-piece polytrope eos
   REAL8 logp1=0.0;
   REAL8 gamma1=0.0;
   REAL8 gamma2=0.0;
   REAL8 gamma3=0.0;
   if(LALInferenceGetProcParamVal(commandLine,"--inj-logp1")&&LALInferenceGetProcParamVal(commandLine,"--inj-gamma1")&&LALInferenceGetProcParamVal(commandLine,"--inj-gamma2")&&LALInferenceGetProcParamVal(commandLine,"--inj-gamma3")){
     logp1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-logp1")->value);
     gamma1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-gamma1")->value);
     gamma2= atof(LALInferenceGetProcParamVal(commandLine,"--inj-gamma2")->value);
     gamma3= atof(LALInferenceGetProcParamVal(commandLine,"--inj-gamma3")->value);
     // Find lambda1,2(m1,2|eos)
     LALInferenceLogp1GammasMasses2Lambdas(logp1, gamma1, gamma2, gamma3, inj_table->mass1, inj_table->mass2, &lambda1, &lambda2);
     fprintf(stdout,"Injection logp1 set to %f\n",logp1);
     fprintf(stdout,"Injection gamma1 set to %f\n",gamma1);
     fprintf(stdout,"Injection gamma2 set to %f\n",gamma2);
     fprintf(stdout,"Injection gamma3 set to %f\n",gamma3);
     fprintf(stdout,"lambda1 set to %f\n",lambda1);
     fprintf(stdout,"lambda2 set to %f\n",lambda2);

     /*
     For when tagSimInspiralTable is updated
     to include EOS params

     injEvent->logp1= logp1;
     injEvent->gamma1= gamma1;
     injEvent->gamma2= gamma2;
     injEvent->gamma3= gamma3;
     */
   }

   // Inject 4-coef. spectral eos
   REAL8 SDgamma0=0.0;
   REAL8 SDgamma1=0.0;
   REAL8 SDgamma2=0.0;
   REAL8 SDgamma3=0.0;
   if(LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma0") && LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma1") && LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma2") && LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma3")){
     SDgamma0= atof(LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma0")->value);
     SDgamma1= atof(LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma1")->value);
     SDgamma2= atof(LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma2")->value);
     SDgamma3= atof(LALInferenceGetProcParamVal(commandLine,"--inj-SDgamma3")->value);
     REAL8 gamma[]={SDgamma0,SDgamma1,SDgamma2,SDgamma3};
     LALInferenceSDGammasMasses2Lambdas(gamma,inj_table->mass1,inj_table->mass2,&lambda1,&lambda2,4);
     fprintf(stdout,"Injection SDgamma0 set to %lf\n",SDgamma0);
     fprintf(stdout,"Injection SDgamma1 set to %lf\n",SDgamma1);
     fprintf(stdout,"Injection SDgamma2 set to %lf\n",SDgamma2);
     fprintf(stdout,"Injection SDgamma3 set to %lf\n",SDgamma3);
     fprintf(stdout,"lambda1 set to %f\n",lambda1);
     fprintf(stdout,"lambda2 set to %f\n",lambda2);
   }


	/* FIXME: One also need to code the same manipulations done to f_max in LALInferenceTemplateXLALSimInspiralChooseWaveform line 856 circa*/

  /* Set up LAL dictionary */
  LALDict* LALpars=XLALCreateDict();

  /* Set the spin-frame convention */
  ppt = LALInferenceGetProcParamVal(commandLine,"--inj-spin-frame");
  if(ppt) {
      if (!strcmp(ppt->value, "view"))
          XLALSimInspiralWaveformParamsInsertFrameAxis(LALpars, LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW);
      else if (!strcmp(ppt->value, "orbital-l"))
          XLALSimInspiralWaveformParamsInsertFrameAxis(LALpars, LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L);
      else if (!strcmp(ppt->value, "total-j"))
          XLALSimInspiralWaveformParamsInsertFrameAxis(LALpars, LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J);
  }

  LALSimInspiralSpinOrder spinO = LAL_SIM_INSPIRAL_SPIN_ORDER_ALL;
  if(LALInferenceGetProcParamVal(commandLine, "--inj-spinOrder")) {
    spinO = atoi(LALInferenceGetProcParamVal(commandLine, "--inj-spinOrder")->value);
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(LALpars, spinO);
  }

  LALSimInspiralTidalOrder tideO = LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL;
  if(LALInferenceGetProcParamVal(commandLine, "--inj-tidalOrder")) {
    tideO = atoi(LALInferenceGetProcParamVal(commandLine, "--inj-tidalOrder")->value);
    XLALSimInspiralWaveformParamsInsertPNTidalOrder(LALpars, tideO);
  }

  REAL8 deltaT = IFOdata->timeData->deltaT;
  REAL8 deltaF = IFOdata->freqData->deltaF;

  REAL8 f_min = XLALSimInspiralfLow2fStart(inj_table->f_lower, amp_order, approximant);
  REAL8 f_max = 0.0;

  REAL8 fref = 100.;
  if(LALInferenceGetProcParamVal(commandLine,"--inj-fref")) {
    fref = atoi(LALInferenceGetProcParamVal(commandLine,"--inj-fref")->value);
  }
  
  if (approximant == TaylorF2)
      f_max = 0.0; /* this will stop at ISCO */
  else{
	  /* get the max f_max as done in Template. Has to do this since I cannot  access LALInferenceModel here*/
	  dataPtr=IFOdata;
	  while(dataPtr){
      	if(dataPtr->fHigh>f_max) f_max=dataPtr->fHigh;
		dataPtr=dataPtr->next;
	  }
  }
 /* Print a line with information about approximant, amp_order, phaseorder, tide order and spin order */
  fprintf(stdout,"\n\n---\t\t ---\n");
 fprintf(stdout,"Injection will run using Approximant %i (%s), phase order %i, amp order %i, spin order %i, tidal order %i, in the frequency domain.\n",approximant,XLALSimInspiralGetStringFromApproximant(approximant),phase_order,amp_order,(int) spinO,(int) tideO);
   fprintf(stdout,"---\t\t ---\n\n");

  COMPLEX16FrequencySeries *hptilde=NULL, *hctilde=NULL;

  XLALSimInspiralWaveformParamsInsertTidalLambda1(LALpars,lambda1);
  XLALSimInspiralWaveformParamsInsertTidalLambda2(LALpars,lambda2);
  XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(LALpars,amp_order);
  XLALSimInspiralWaveformParamsInsertPNPhaseOrder(LALpars,phase_order);

  XLALSimInspiralChooseFDWaveform(&hptilde, &hctilde, inj_table->mass1*LAL_MSUN_SI, inj_table->mass2*LAL_MSUN_SI,
				  inj_table->spin1x, inj_table->spin1y, inj_table->spin1z,
				  inj_table->spin2x, inj_table->spin2y, inj_table->spin2z,
				  inj_table->distance*LAL_PC_SI * 1.0e6, inj_table->inclination,
				  inj_table->coa_phase, 0., 0., 0., deltaF, f_min, f_max, fref,
				  LALpars, approximant);
  XLALDestroyDict(LALpars);

  /* Fail if injection waveform generation was not successful */
  errnum = *XLALGetErrnoPtr();
  if (errnum != XLAL_SUCCESS) {
    XLALPrintError(" ERROR in InjectFD(): error encountered when injecting waveform. errnum=%d\n",errnum);
    exit(1);
  }

  if((ppt=LALInferenceGetProcParamVal(commandLine,"--dump-geocenter-pols"))) {
    fprintf(stdout,"Dump injected FreqDomain h_plus and h_cross at geocenter (for IFO %s)\n", IFOdata->name);
    char filename[320];
    sprintf(filename,"%s_FD_geocenter_pols.dat",IFOdata->name);
    FILE* file=fopen(filename, "w");
    if(!file){
      fprintf(stderr,"Unable to open the path %s for writing injected FreqDomain h_plus and h_cross\n",filename);
    exit(1);
    }
    UINT4 j;
    for(j=0; j<hptilde->data->length; j++){     
      fprintf(file, "%10.10g %10.10g %10.10g %10.10g %10.10g\n", deltaF*j, 
        creal(hptilde->data->data[j]), cimag(hptilde->data->data[j]), 
        creal(hctilde->data->data[j]), cimag(hctilde->data->data[j]));
    }
    fclose(file);
  }

  REAL8 Fplus, Fcross;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  LIGOTimeGPS GPSlal;
  REAL8 gmst;
  REAL8 chisquared;
  REAL8 timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  REAL8 timeshift;  /* time shift (not necessarily same as above)                   */
  REAL8 twopit, re, im, dre, dim, newRe, newIm;
  UINT4 i, lower, upper;

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

    char InjFileName[320];
    sprintf(InjFileName,"injection_%s.dat",dataPtr->name);
    FILE *outInj=fopen(InjFileName,"w");

     /* determine frequency range & loop over frequency bins: */
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
    chisquared = 0.0;

    re = cos(twopit * deltaF * lower);
    im = -sin(twopit * deltaF * lower);

    /* When analysing TD data, FD templates need to be amplified by the window power loss factor */
    double windowFactor;
    windowFactor=sqrt(dataPtr->window->data->length/dataPtr->window->sumofsquares);

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

      /* apply windowing factor */
      templateReal *= ((REAL8) windowFactor);
      templateImag *= ((REAL8) windowFactor);

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
  ppt=LALInferenceGetProcParamVal(commandLine,"--dont-dump-extras");
  if (!ppt){
    PrintSNRsToFile(IFOdata , SNRpath);
  }
  XLALDestroyCOMPLEX16FrequencySeries(hctilde);
  XLALDestroyCOMPLEX16FrequencySeries(hptilde);
}

static void PrintSNRsToFile(LALInferenceIFOData *IFOdata , char SNRpath[] ){
  REAL8 NetSNR=0.0;
  LALInferenceIFOData *thisData=IFOdata;
  int nIFO=0;

  while(thisData){
    thisData=thisData->next;
    nIFO++;
  }
  FILE * snrout = fopen(SNRpath,"w");
  if(!snrout){
    fprintf(stderr,"Unable to open the path %s for writing SNR files\n",SNRpath);
    fprintf(stderr,"Error code %i: %s\n",errno,strerror(errno));
    exit(errno);
  }
  thisData=IFOdata;
  while(thisData){
    fprintf(snrout,"%s:\t %4.2f\n",thisData->name,thisData->SNR);
    nIFO++;
    NetSNR+=(thisData->SNR*thisData->SNR);
    thisData=thisData->next;
  }
  if (nIFO>1){
    fprintf(snrout,"Network:\t");
    fprintf(snrout,"%4.2f\n",sqrt(NetSNR));
  }
  fclose(snrout);
}

/**
* Fill the variables passed in vars with the parameters of the injection passed in event
* will over-write and destroy any existing parameters. Param vary type will be fixed
*/
void LALInferenceInjectionToVariables(SimInspiralTable *theEventTable, LALInferenceVariables *vars)
{
  //UINT4 spinCheck=0;
  if(!vars) {
  XLALPrintError("Encountered NULL variables pointer");
  XLAL_ERROR_VOID(XLAL_EINVAL);
  }
  enforce_m1_larger_m2(theEventTable);
  REAL8 q = theEventTable->mass2 / theEventTable->mass1;
  if (q > 1.0) q = 1.0/q;

  REAL8 psi = theEventTable->polarization;
  if (psi>=M_PI) psi -= M_PI;

  REAL8 injGPSTime = XLALGPSGetREAL8(&(theEventTable->geocent_end_time));

  REAL8 dist = theEventTable->distance;
  //REAL8 cosinclination = cos(theEventTable->inclination);
  REAL8 phase = theEventTable->coa_phase;
  REAL8 dec = theEventTable->latitude;
  REAL8 ra = theEventTable->longitude;

  Approximant injapprox = XLALGetApproximantFromString(theEventTable->waveform);
  LALPNOrder order = XLALGetOrderFromString(theEventTable->waveform);

  REAL8 m1=theEventTable->mass1;
  REAL8 m2=theEventTable->mass2;
  REAL8 chirpmass = theEventTable->mchirp;

  LALInferenceAddVariable(vars, "chirpmass", &chirpmass, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "q", &q, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  if  (LALInferenceCheckVariable(vars,"distance"))
  	LALInferenceAddVariable(vars, "distance", &dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  else if  (LALInferenceCheckVariable(vars,"logdistance")){
	  REAL8 logdistance=log(dist);
  	LALInferenceAddVariable(vars, "logdistance", &logdistance, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  }
  LALInferenceAddVariable(vars, "polarisation", &(psi), LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "phase", &phase, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

  /* Those will work even if the user is working with the detector-frame variables because SKY_FRAME is set 
  to zero while calculating the injected logL in LALInferencePrintInjectionSample */  
  LALInferenceAddVariable(vars, "declination", &dec, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "rightascension", &ra, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "time", &injGPSTime, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);


  LALInferenceAddVariable(vars, "LAL_APPROXIMANT", &injapprox, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "LAL_PNORDER",&order,LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(vars, "LAL_AMPORDER", &(theEventTable->amp_order), LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);

	REAL8 thetaJN,phiJL,theta1,theta2,phi12,chi1,chi2;
	/* Convert cartesian spin coordinates to system-frame variables*/
	
	/* This fref is --inj-fref, which has been set previously by LALInferencePrintInjectionSample
	I don't call it inj_fref since LALInferenceTemplate looks for fref, and that is what will be called to calculate 
	the logL at injval
	*/
	REAL8 fref=100.0;
	if (LALInferenceCheckVariable(vars,"f_ref"))
		fref= *(REAL8*)  LALInferenceGetVariable(vars,"f_ref");
	
	XLALSimInspiralTransformPrecessingWvf2PE(&thetaJN,&phiJL,&theta1,&theta2,&phi12,&chi1,&chi2,theEventTable->inclination,theEventTable->spin1x,theEventTable->spin1y,theEventTable->spin1z,  theEventTable->spin2x, theEventTable->spin2y, theEventTable->spin2z,m1,m2,fref,phase);
	
	if (LALInferenceCheckVariable(vars,"a_spin1"))
		LALInferenceAddVariable(vars,"a_spin1", &chi1, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	if (LALInferenceCheckVariable(vars,"a_spin2"))
		LALInferenceAddVariable(vars,"a_spin2", &chi2, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	if (LALInferenceCheckVariable(vars,"tilt_spin1"))
		LALInferenceAddVariable(vars,"tilt_spin1", &theta1, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	if (LALInferenceCheckVariable(vars,"tilt_spin2"))
		LALInferenceAddVariable(vars,"tilt_spin2", &theta2, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	if (LALInferenceCheckVariable(vars,"phi_jl"))
		LALInferenceAddVariable(vars,"phi_jl", &phiJL, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	if (LALInferenceCheckVariable(vars,"phi12"))
		LALInferenceAddVariable(vars,"phi12", &phi12, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	REAL8 costhetajn=cos(thetaJN);
	LALInferenceAddVariable(vars, "costheta_jn", &costhetajn, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

}

LALInferenceVariables *LALInferencePrintInjectionSample(LALInferenceRunState *runState) {
    int errnum=0;
    char *fname=NULL;
    char defaultname[]="injection_params.dat";
    FILE *outfile=NULL;

    /* check if the --inj argument has been passed */
    ProcessParamsTable *ppt = LALInferenceGetProcParamVal(runState->commandLine,"--inj");
    if (!ppt)
        return(NULL);

    SimInspiralTable *injTable=NULL, *theEventTable=NULL;
    LALInferenceModel *model = LALInferenceInitCBCModel(runState);
    if (LALInferenceGetProcParamVal(runState->commandLine, "--roqtime_steps")){
      LALInferenceSetupROQmodel(model, runState->commandLine);
      fprintf(stderr, "done LALInferenceSetupROQmodel\n");
    } else {
      model->roq_flag=0;
    }
    LALInferenceVariables *injparams = XLALCalloc(1, sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(model->params, injparams);

    SimInspiralTableFromLIGOLw(&injTable, ppt->value, 0, 0);

    ppt = LALInferenceGetProcParamVal(runState->commandLine, "--outfile");
    if (ppt) {
        fname = XLALCalloc((strlen(ppt->value)+255)*sizeof(char),1);
        sprintf(fname,"%s.injection",ppt->value);
    }
    else
        fname = defaultname;

    ppt = LALInferenceGetProcParamVal(runState->commandLine, "--event");
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

    LALPNOrder *order = LALInferenceGetVariable(injparams, "LAL_PNORDER");
    Approximant *approx = LALInferenceGetVariable(injparams, "LAL_APPROXIMANT");

    if (!(approx && order)){
        fprintf(stdout,"Unable to print injection sample: No approximant/PN order set\n");
        return(NULL);
    }

    REAL8 fref = 100.;
    if(LALInferenceGetProcParamVal(runState->commandLine,"--inj-fref")) {
      fref = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--inj-fref")->value);
    }
	LALInferenceAddVariable(injparams,"f_ref",(void *)&fref,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

	UINT4 azero=0;
	LALInferenceAddVariable(injparams,"SKY_FRAME",(void *)&azero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	/* remove eventual SKY FRAME vars since they will contain garbage*/
	if (LALInferenceCheckVariable(injparams,"t0"))
		LALInferenceRemoveVariable(injparams,"t0");
	if (LALInferenceCheckVariable(injparams,"cosalpha"))
		LALInferenceRemoveVariable(injparams,"cosalpha");
	if (LALInferenceCheckVariable(injparams,"azimuth"))
		LALInferenceRemoveVariable(injparams,"azimuth");
	
    /* Fill named variables */
    LALInferenceInjectionToVariables(theEventTable, injparams);

    REAL8 injPrior = runState->prior(runState, injparams, model);
    LALInferenceAddVariable(injparams, "logPrior", &injPrior, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
    REAL8 injL=0.;
    if ( (int) *approx == XLALGetApproximantFromString(theEventTable->waveform)){
        XLAL_TRY(injL = runState->likelihood(injparams, runState->data, model), errnum);
        if(errnum){
            fprintf(stderr,"ERROR: Cannot print injection sample. Received error code %s\n",XLALErrorString(errnum));
        }
    }
    LALInferenceAddVariable(injparams, "logL", (void *)&injL,LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
	REAL8 logZnoise=LALInferenceNullLogLikelihood(runState->data);
    REAL8 tmp2=injL-logZnoise;
    LALInferenceAddVariable(injparams,"deltalogL",(void *)&tmp2,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    
    LALInferenceIFOData *data=runState->data;
    while(data) {
        char tmpName[320];
        REAL8 tmp=model->loglikelihood - data->nullloglikelihood;
        sprintf(tmpName,"deltalogl%s",data->name);
        LALInferenceAddVariable(injparams, tmpName, &tmp, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
        data=data->next;
    }
	
    /* Save to file */
    outfile=fopen(fname,"w");
    if(!outfile) {fprintf(stderr,"ERROR: Unable to open file %s for injection saving\n",fname); exit(1);}
    LALInferenceSortVariablesByName(injparams);
    LALInferenceFprintParameterHeaders(outfile, injparams);
    fprintf(outfile,"\n");
    LALInferencePrintSample(outfile, injparams);

    fclose(outfile);
    //LALInferenceClearVariables(injparams);
    return(injparams);
}

int enforce_m1_larger_m2(SimInspiralTable* injEvent){
    /* Template generator assumes m1>=m2 thus we must enfore the same convention while injecting, otherwise spin2 will be assigned to mass1
    *        We also shift the phase by pi to be sure the same WF in injected
	*        Returns: 1 if masses were flipped
    */
    REAL8 m1,m2,tmp;
    m1=injEvent->mass1;
    m2=injEvent->mass2;

    if (m1>=m2) return(0);
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
        return(1);
        }
}

void LALInferenceSetupROQmodel(LALInferenceModel *model, ProcessParamsTable *commandLine){

  LALStatus status;
  memset(&status,0,sizeof(status));
  UINT4 q=0;
  UINT4 event=0;
  ProcessParamsTable *procparam=NULL,*ppt=NULL;
  SimInspiralTable *injTable=NULL;
  FILE *tempfp;
  unsigned int n_basis_linear=0, n_basis_quadratic=0, n_samples=0, time_steps=0;
  int errsave;

  LIGOTimeGPS GPStrig;
  REAL8 endtime=0.0;

	  model->roq = XLALMalloc(sizeof(LALInferenceROQModel));
	  model->roq_flag = 1;
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
		if(injTable->simulation_id == atol(procparam->value)) break;
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
	    XLALStrToGPS(&GPStrig,procparam->value,NULL);
	  }
	  else{
	    if(injTable) memcpy(&GPStrig,&(injTable->geocent_end_time),sizeof(GPStrig));
	    else {
	      fprintf(stderr,">>> Error: No trigger time specifed and no injection given \n");
	      exit(1);
	    }
	  }

	  endtime=XLALGPSGetREAL8(&GPStrig);

	  if(LALInferenceGetProcParamVal(commandLine,"--roqtime_steps")){
	    ppt=LALInferenceGetProcParamVal(commandLine,"--roqtime_steps");
	    tempfp = fopen (ppt->value,"r");
	    if (tempfp == NULL){
	        errsave = errno;
		fprintf(stderr, "Error: cannot find file %s \n", ppt->value);
		fprintf(stderr, "Error code %i: %s\n", errsave, strerror(errsave));
		exit(errsave);
	    }
	    fscanf(tempfp, "%u", &time_steps);
	    fscanf(tempfp, "%u", &n_basis_linear);
	    fscanf(tempfp, "%u", &n_basis_quadratic);
	    fscanf(tempfp, "%u", &n_samples);
	    fclose(tempfp);
	    fprintf(stderr, "loaded --roqtime_steps\n");
	  }


	  model->roq->frequencyNodesLinear = XLALCreateREAL8Sequence(n_basis_linear);
	  model->roq->frequencyNodesQuadratic = XLALCreateREAL8Sequence(n_basis_quadratic);

	  model->roq->trigtime = endtime;

          model->roq->frequencyNodesLinear = XLALCreateREAL8Sequence(n_basis_linear);
          model->roq->frequencyNodesQuadratic = XLALCreateREAL8Sequence(n_basis_quadratic);
          model->roq->calFactorLinear = XLALCreateCOMPLEX16Sequence(model->roq->frequencyNodesLinear->length);
          model->roq->calFactorQuadratic = XLALCreateCOMPLEX16Sequence(model->roq->frequencyNodesQuadratic->length);

	  if(LALInferenceGetProcParamVal(commandLine,"--roqnodesLinear")){
	    ppt=LALInferenceGetProcParamVal(commandLine,"--roqnodesLinear");

	    model->roq->nodesFileLinear = fopen(ppt->value, "rb");
	    if (!(model->roq->nodesFileLinear)) {
	        errsave = errno;
		fprintf(stderr, "Error: cannot find file %s \n", ppt->value);
		fprintf(stderr, "Error code %i: %s\n", errsave, strerror(errsave));
		exit(errsave);
	    } // check file exists
	    fprintf(stderr, "read model->roq->frequencyNodesLinear");

	    for(unsigned int linsize = 0; linsize < n_basis_linear; linsize++){
	      fread(&(model->roq->frequencyNodesLinear->data[linsize]), sizeof(REAL8), 1, model->roq->nodesFileLinear);
	    }
	    fclose(model->roq->nodesFileLinear);
	    model->roq->nodesFileLinear = NULL;
	    fprintf(stderr, "loaded --roqnodesLinear\n");
	  }

	  if(LALInferenceGetProcParamVal(commandLine,"--roqnodesQuadratic")){
	    ppt=LALInferenceGetProcParamVal(commandLine,"--roqnodesQuadratic");
	    model->roq->nodesFileQuadratic = fopen(ppt->value, "rb");
	    if (!(model->roq->nodesFileQuadratic)) {
	      errsave = errno;
	      fprintf(stderr, "Error: cannot find file %s \n", ppt->value);
	      fprintf(stderr, "Error code %i: %s\n", errsave, strerror(errsave));
	      exit(errsave);
	    } // check file exists

	    for(unsigned int quadsize = 0; quadsize < n_basis_quadratic; quadsize++){
	      fread(&(model->roq->frequencyNodesQuadratic->data[quadsize]), sizeof(REAL8), 1, model->roq->nodesFileQuadratic);
	    }
	    fclose(model->roq->nodesFileQuadratic);
	    model->roq->nodesFileQuadratic = NULL;
	    fprintf(stderr, "loaded --roqnodesQuadratic\n");



	  }
}

void LALInferenceSetupROQdata(LALInferenceIFOData *IFOdata, ProcessParamsTable *commandLine){
  LALStatus status;
  memset(&status,0,sizeof(status));
  LALInferenceIFOData *thisData=IFOdata;
  UINT4 q=0;
  UINT4 event=0;
  ProcessParamsTable *procparam=NULL,*ppt=NULL;
  SimInspiralTable *injTable=NULL;
  unsigned int n_basis_linear, n_basis_quadratic, n_samples, time_steps;
  float dt=0.1;
  //REAL8 timeMin=0.0,timeMax=0.0;
  FILE *tempfp;
  char tmp[320];
  int errsave;

  procparam=LALInferenceGetProcParamVal(commandLine,"--inj");
  if(procparam) {
    SimInspiralTableFromLIGOLw(&injTable,procparam->value,0,0);
    if (!injTable) {
      fprintf(stderr,"Unable to open injection file(LALInferenceReadData) %s\n",procparam->value);
      exit(1);
    }
    procparam=LALInferenceGetProcParamVal(commandLine,"--event");
    if (procparam) {
      event=atoi(procparam->value);
      while(q<event) {
	q++;
	injTable=injTable->next;
      }
    } else if ((procparam=LALInferenceGetProcParamVal(commandLine,"--event-id"))) {
      while (injTable) {
	if (injTable->simulation_id == atol(procparam->value)) {
	  break;
	} else {
	  injTable=injTable->next;
	}
      }
      if (!injTable) {
	fprintf(stderr, "Error, cannot find simulation id %s in injection file\n", procparam->value);
	exit(1);
      }
    }
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
  if(ppt){
    dt=atof(ppt->value);
  }

  if (LALInferenceGetProcParamVal(commandLine,"--roqtime_steps")) {
    ppt = LALInferenceGetProcParamVal(commandLine,"--roqtime_steps");
    tempfp = fopen (ppt->value,"r");
    if (tempfp == NULL){
      errsave = errno;
      fprintf(stderr, "Error: cannot find file %s \n", ppt->value);
      fprintf(stderr, "Error code %i: %s\n", errsave, strerror(errsave));
      exit(errsave);
    }
    fscanf(tempfp, "%u", &time_steps);
    fscanf(tempfp, "%u", &n_basis_linear);
    fscanf(tempfp, "%u", &n_basis_quadratic);
    fscanf(tempfp, "%u", &n_samples);
    fclose(tempfp);
    fprintf(stderr, "loaded --roqtime_steps\n");
  }


  thisData=IFOdata;
    while (thisData) {
      thisData->roq = XLALMalloc(sizeof(LALInferenceROQData));

      thisData->roq->weights_linear = XLALMalloc(n_basis_linear*sizeof(LALInferenceROQSplineWeights));

      sprintf(tmp, "--%s-roqweightsLinear", thisData->name);
      ppt = LALInferenceGetProcParamVal(commandLine,tmp);

      thisData->roq->weightsFileLinear = fopen(ppt->value, "rb");
      if (thisData->roq->weightsFileLinear == NULL){
	errsave = errno;
	fprintf(stderr, "Error: cannot find file %s \n", ppt->value);
	fprintf(stderr, "Error code %i: %s\n", errsave, strerror(errsave));
	exit(errsave);
      }
      thisData->roq->weightsLinear = (double complex*)malloc(n_basis_linear*time_steps*(sizeof(double complex)));

      //0.045 comes from the diameter of the earth in light seconds: the maximum time-delay between earth-based observatories
      thisData->roq->time_weights_width = 2*dt + 2*0.045;
      thisData->roq->time_step_size = thisData->roq->time_weights_width/time_steps;
      thisData->roq->n_time_steps = time_steps;


      fprintf(stderr, "basis_size = %d\n", n_basis_linear);
      fprintf(stderr, "time steps = %d\n", time_steps);

      double *tmp_real_weight = malloc(time_steps*(sizeof(double)));
      double *tmp_imag_weight = malloc(time_steps*(sizeof(double)));

      double *tmp_tcs = malloc(time_steps*(sizeof(double)));

      sprintf(tmp, "--roq-times");
      ppt = LALInferenceGetProcParamVal(commandLine,tmp);

      FILE *tcFile = fopen(ppt->value, "rb");
      if (tcFile == NULL) {
	errsave = errno;
	fprintf(stderr, "Error: cannot find file %s \n", ppt->value);
	fprintf(stderr, "Error code %i: %s\n", errsave, strerror(errsave));
	exit(errsave);
      }

      for(unsigned int gg=0;gg < time_steps; gg++){
	fread(&(tmp_tcs[gg]), sizeof(double), 1, tcFile);
      }

      for(unsigned int ii=0; ii<n_basis_linear;ii++){
	for(unsigned int jj=0; jj<time_steps;jj++){
	  fread(&(thisData->roq->weightsLinear[ii*time_steps + jj]), sizeof(double complex), 1, thisData->roq->weightsFileLinear);
	  tmp_real_weight[jj] = creal(thisData->roq->weightsLinear[ii*time_steps + jj]);
	  tmp_imag_weight[jj] = cimag(thisData->roq->weightsLinear[ii*time_steps + jj]);
	}
  //gsl_interp_accel is not thread-safe, and each OpenMP thread will need its
  //own gsl_interp_accel object.
	thisData->roq->weights_linear[ii].acc_real_weight_linear = NULL;
 	thisData->roq->weights_linear[ii].acc_imag_weight_linear = NULL;

	thisData->roq->weights_linear[ii].spline_real_weight_linear = gsl_spline_alloc (gsl_interp_cspline, time_steps);
        gsl_spline_init(thisData->roq->weights_linear[ii].spline_real_weight_linear, tmp_tcs, tmp_real_weight, time_steps);

        thisData->roq->weights_linear[ii].spline_imag_weight_linear = gsl_spline_alloc (gsl_interp_cspline, time_steps);
        gsl_spline_init(thisData->roq->weights_linear[ii].spline_imag_weight_linear, tmp_tcs, tmp_imag_weight, time_steps);
      }
      fclose(thisData->roq->weightsFileLinear);
      thisData->roq->weightsFileLinear = NULL;
      fclose(tcFile);

      sprintf(tmp, "--%s-roqweightsQuadratic", thisData->name);
      ppt = LALInferenceGetProcParamVal(commandLine,tmp);
      thisData->roq->weightsQuadratic = (double*)malloc(n_basis_quadratic*sizeof(double));
      thisData->roq->weightsFileQuadratic = fopen(ppt->value, "rb");
      if (thisData->roq->weightsFileQuadratic == NULL){
	errsave = errno;
	fprintf(stderr, "Error: cannot find file %s \n", ppt->value);
	fprintf(stderr, "Error code %i: %s\n", errsave, strerror(errsave));
	exit(errsave);
      }
      for(unsigned int ii=0; ii<n_basis_quadratic;ii++){
	fread(&(thisData->roq->weightsQuadratic[ii]), sizeof(double), 1, thisData->roq->weightsFileQuadratic);
      }
      fclose(thisData->roq->weightsFileQuadratic);
      thisData->roq->weightsFileQuadratic = NULL;
      fprintf(stderr, "loaded %s ROQ weights\n", thisData->name);
      thisData = thisData->next;
    }
}

static void LALInferenceSetGPSTrigtime(LIGOTimeGPS *GPStrig, ProcessParamsTable *commandLine){

    ProcessParamsTable *procparam;
    SimInspiralTable *inspiralTable=NULL;
    SimBurst *burstTable=NULL;
    UINT4 event=0;
    UINT4 q=0;
    LALStatus status;
    memset(&status,0,sizeof(LALStatus));

    /* First check if trigtime has been given as an option */
    if(LALInferenceGetProcParamVal(commandLine,"--trigtime")){
        procparam=LALInferenceGetProcParamVal(commandLine,"--trigtime");
        XLALStrToGPS(GPStrig,procparam->value,NULL);
        fprintf(stdout,"Set trigtime to %.10f\n",GPStrig->gpsSeconds+1.0e-9 * GPStrig->gpsNanoSeconds);
        return;

    }
    else{
        /* If not check if we have an injtable passed with --inj */

        if(LALInferenceGetProcParamVal(commandLine,"--injXML"))
        {
            XLALPrintError("ERROR: --injXML option is deprecated. Use --inj and update your scripts\n");
            exit(1);
        }
        if((procparam=LALInferenceGetProcParamVal(commandLine,"--inj"))){
            fprintf(stdout,"Checking if the xml table is an inspiral table... \n");
            /* Check if it is a SimInspiralTable */
            SimInspiralTableFromLIGOLw(&inspiralTable,procparam->value,0,0);

            if (inspiralTable){
                procparam=LALInferenceGetProcParamVal(commandLine,"--event");
                if(procparam) {
                event=atoi(procparam->value);
                while(q<event) {q++; inspiralTable=inspiralTable->next;}
                }
                else if ((procparam=LALInferenceGetProcParamVal(commandLine,"--event-id")))
                {
                while(inspiralTable)
                {
                if(inspiralTable->simulation_id == atol(procparam->value)) break;
                else inspiralTable=inspiralTable->next;
                }
                if(!inspiralTable){
                fprintf(stderr,"Error, cannot find simulation id %s in injection file\n",procparam->value);
                exit(1);
                }
                }
                else
                fprintf(stdout,"You did not provide an event number with the injtable. Using event 0 which may not be what you want!!!!!\n");
                memcpy(GPStrig,&(inspiralTable->geocent_end_time),sizeof(LIGOTimeGPS));
                printf("Set inspiral injtime %.10f\n",inspiralTable->geocent_end_time.gpsSeconds+1.0e-9* inspiralTable->geocent_end_time.gpsNanoSeconds);
                return;
            }
        }
        else if((procparam=LALInferenceGetProcParamVal(commandLine,"--binj"))){
            /* Check if it is a SimBurst table */
            fprintf(stdout,"Checking if the xml table is a burst table... \n");
            burstTable=XLALSimBurstTableFromLIGOLw(procparam->value,0,0);
            if(burstTable){
                procparam=LALInferenceGetProcParamVal(commandLine,"--event");
                if(procparam) {
                    event=atoi(procparam->value);
                    while(q<event) {q++; burstTable=burstTable->next;}
                }
                else if ((procparam=LALInferenceGetProcParamVal(commandLine,"--event-id")))
                {
                    fprintf(stderr,"Error, SimBurst tables do not currently support event_id tags \n");
                    exit(1);
                }
                else
                    fprintf(stdout,"You did not provide an event number with the injtable. Using event 0 which may not be what you want!!!!!\n");
                memcpy(GPStrig,&(burstTable->time_geocent_gps),sizeof(LIGOTimeGPS));
                fprintf(stdout,"Set trigtime from burstable to %.10f\n",GPStrig->gpsSeconds+1.0e-9 * GPStrig->gpsNanoSeconds);
                return;
            }
        }
        else if(!LALInferenceGetProcParamVal(commandLine,"--segment-start")){
            XLALPrintError("Error: No trigger time specifed and no injection given \n");
            //XLAL_ERROR_NULL(XLAL_EINVAL);
            exit(1);
        }

    }
}

void LALInferenceInjectFromMDC(ProcessParamsTable *commandLine, LALInferenceIFOData *IFOdata){

    /* Read time domain WF present in an mdc frame file, FFT it and inject into the frequency domain stream */

    char mdcname[]="GW";
    char **mdc_caches=NULL;
    char **mdc_channels=NULL;
    ProcessParamsTable * ppt=commandLine;

    UINT4 nIFO=0;
    int i=0;
    UINT4 j=0;
    LALInferenceIFOData *data=IFOdata;
    REAL8 prefactor =1.0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--mdc-prefactor");
    if (ppt){

        prefactor=atof(ppt->value);
        fprintf(stdout,"Using prefactor=%f to scale the MDC injection\n",prefactor);
    }

    ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
    if (ppt){

        fprintf(stderr,"You cannot use both injfile (--inj) and MDCs (--inject_from_mdc) Exiting... \n");
        exit(1);

    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--binj");
    if (ppt){

        fprintf(stderr,"You cannot use both injfile (--binj) and MDCs (--inject_from_mdc) Exiting... \n");
        exit(1);

    }

    REAL8 tmp=0.0;
    REAL8 net_snr=0.0;
    while (data) {nIFO++; data=data->next;}
    UINT4 Nmdc=0,Nchannel=0;

    char mdc_caches_name[] = "injcache";
    char mdc_channels_name[] = "injchannel";
    char **IFOnames=NULL;
    INT4 rlceops= getNamedDataOptionsByDetectors(commandLine, &IFOnames,&mdc_caches ,mdc_caches_name, &Nmdc);
    if (!rlceops){
      fprintf(stderr,"Must provide a --IFO-injcache option for each IFO if --inject_from_mdc is given\n");
      exit(1);
    }

    rlceops= getNamedDataOptionsByDetectors(commandLine, &IFOnames,&mdc_channels ,mdc_channels_name, &Nchannel);
    if (!rlceops){
        fprintf(stdout,"WARNING: You did not provide the name(s) of channel(s) to use with the injection mdc. Using the default which may not be what you want!\n");
        mdc_channels=  malloc((nIFO+1)*sizeof(char*));
        data=IFOdata;
        i=0;
        while (data){
           mdc_channels[i] =  malloc(512*sizeof(char));
            if(!strcmp(data->name,"H1")) {
               sprintf(mdc_channels[i],"H1:%s-H",mdcname);}
            else if(!strcmp(data->name,"L1")) {
                 sprintf(mdc_channels[i],"L1:%s-H",mdcname); }
            else if(!strcmp(data->name,"V1")) {
                 sprintf(mdc_channels[i],"V1:%s-16K",mdcname);}
            data=data->next;
            i++;

            }
    }

    LIGOTimeGPS epoch=IFOdata->timeData->epoch;
    REAL8 deltaT=IFOdata->timeData->deltaT ;
    int seglen=IFOdata->timeData->data->length;
    REAL8 SampleRate=4096.0,SegmentLength=0.0;
    if(LALInferenceGetProcParamVal(commandLine,"--srate")) SampleRate=atof(LALInferenceGetProcParamVal(commandLine,"--srate")->value);
    SegmentLength=(REAL8) seglen/SampleRate;

    REAL8TimeSeries * timeData=NULL;
    REAL8TimeSeries * windTimeData=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("WindMDCdata",&epoch,0.0,deltaT,&lalDimensionlessUnit,(size_t)seglen);
    COMPLEX16FrequencySeries* injF=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("injF",&IFOdata->timeData->epoch,0.0,IFOdata->freqData->deltaF,&lalDimensionlessUnit,	IFOdata->freqData->data->length);

    if(!injF) {
      XLALPrintError("Unable to allocate memory for injection buffer\n");
      XLAL_ERROR_VOID(XLAL_EFUNC);
    }

    REAL4 WinNorm = sqrt(IFOdata->window->sumofsquares/IFOdata->window->data->length);

    data=IFOdata;
    i=0;
    UINT4 lower = (UINT4)ceil(data->fLow / injF->deltaF);
    UINT4 upper = (UINT4)floor(data->fHigh /injF-> deltaF);
    //FIXME CHECK WNORM
    /* Inject into FD data stream and calculate optimal SNR */
    while(data){
      tmp=0.0;
        LALCache *mdc_cache=NULL;
        mdc_cache  = XLALCacheImport(mdc_caches[i] );

        /* Read MDC frame */
        timeData=readTseries(mdc_cache,mdc_channels[i],epoch,SegmentLength);
        /* downsample */
        XLALResampleREAL8TimeSeries(timeData,1.0/SampleRate);
        /* window timeData and store it in windTimeData */
        XLALDDVectorMultiply(windTimeData->data,timeData->data,IFOdata->window->data);

        /*for(j=0;j< timeData->data->length;j++)
            fprintf(out,"%lf %10.10e %10.10e %10.10e \n",epoch.gpsSeconds + j*deltaT,data->timeData->data->data[j],data->timeData->data->data[j]+timeData->data->data[j],timeData->data->data[j]);
        fclose(out);
        */

        /* set the whole seq to 0 */
        for(j=0;j<injF->data->length;j++) injF->data->data[j]=0.0;

        /* FFT */
        XLALREAL8TimeFreqFFT(injF,windTimeData,IFOdata->timeToFreqFFTPlan);


        for(j=lower;j<upper;j++){
                windTimeData->data->data[j] /= sqrt(data->window->sumofsquares / data->window->data->length);
                /* Add data in freq stream */
                data->freqData->data->data[j]+=crect(prefactor *creal(injF->data->data[j])/WinNorm,prefactor *cimag(injF->data->data[j])/WinNorm);
                tmp+= prefactor*prefactor*(creal(injF ->data->data[j])*creal(injF ->data->data[j])+cimag(injF ->data->data[j])*cimag(injF ->data->data[j]))/data->oneSidedNoisePowerSpectrum->data->data[j];
        }

        tmp*=2.*injF->deltaF;
        printf("Injected SNR %.3f in IFO %s from MDC \n",sqrt(2*tmp),data->name);
        data->SNR=sqrt(2*tmp);
        net_snr+=2*tmp;
        i++;
        data=data->next;
    }
    printf("Injected network SNR %.3f from MDC\n",sqrt(net_snr));

    char SNRpath[FILENAME_MAX+100];
    ppt=LALInferenceGetProcParamVal(commandLine,"--outfile");
    if(!ppt){
      fprintf(stderr,"Must specify --outfile <filename.dat>\n");
      exit(1);
    }
    char *outfile=ppt->value;
    snprintf(SNRpath,sizeof(SNRpath),"%s_snr.txt",outfile);
    ppt=LALInferenceGetProcParamVal(commandLine,"--dont-dump-extras");
    if (!ppt){
      PrintSNRsToFile(IFOdata , SNRpath);
    }
    return ;

}
