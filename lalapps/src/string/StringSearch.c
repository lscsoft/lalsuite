/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon, Patrick Brady, Xavier Siemens
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

/*********************************************************************************/
/*                            Cosmic string search code                          */
/*                                                                               */
/*                           X. Siemens and J. Creighton                         */
/*                                                                               */
/*                                 UWM - July 2004                               */
/*********************************************************************************/

#include <config.h>
#if !defined HAVE_LIBGSL || !defined HAVE_LIBLALFRAME
#include <stdio.h>
int main(void) {fputs("disabled, no gsl or no lal frame library support.\n", stderr);return 1;}
#else

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <glob.h>
#include <errno.h>
#include <getopt.h>
#include <stdarg.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/AVFactories.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/Window.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/AVFactories.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/PrintFTSeries.h>
#include <lal/Random.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/lalGitID.h>

#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>

#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>

#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/GenerateBurst.h>


#include <lalapps.h>
#include <processtable.h>
#include <lalappsGitID.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

#define ADD_PROCESS_PARAM(process, type) \
	do { paramaddpoint = add_process_param(paramaddpoint, process, type, long_options[option_index].name, optarg); } while(0)

#define SCALE 1e20
#define MAXTEMPLATES 1000

NRCSID( STRINGSEARCHC, "StringSearch $Id$");
RCSID( "StringSearch $Id$");

#define PROGRAM_NAME "StringSearch"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"


/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL4 flow;                 /* High pass filtering frequency */
  REAL8 samplerate;           /* desired sample rate */
  REAL4 fbankstart;           /* lowest frequency of templates */
  REAL4 fbankhighfcutofflow;  /* lowest high frequency cut-off */
  REAL4 fmismatchmax;         /* maximal mismatch allowed from 1 template to the next */
  char *FrCacheFile;          /* Frame cache file */
  char *InjectionFile;        /* LIGO xml injection file */
  char *ChannelName;          /* Name of channel to be read in from frames */
  char *outputFileName;       /* Name of xml output filename */
  INT4 GPSStart;              /* GPS start time of segment to be analysed */
  INT4 GPSEnd;                /* GPS end time of segment to be analysed */
  INT4 ShortSegDuration;      /* Number of fixed length sub-segments between GPSStart and GPSEnd */
  REAL4 TruncSecs;            /* Half the number of seconds truncated at beginning and end of a chunk */
  REAL4 power;                /* Kink (-5/3) or cusp (-4/3) frequency power law */
  REAL4 threshold;            /* event SNR threshold */
  INT4 fakenoiseflag;         /* =0 if real noise =1 if fake gaussian noise */
  INT4 whitespectrumflag;     /* =0 if spectrum is to be computed =1 for white spectrum */
  INT4 trigstarttime;         /* start-time of allowed triggers */
  REAL4 cluster;              /* =0.0 if events are not to be clustered = clustering time otherwise */
  INT4 pad;                   /* seconds of padding */
  INT4 printspectrumflag;     /* flag set to 1 if user wants to print the spectrum */
  INT4 printfilterflag;       /* flag set to 1 if user wants to print the filter in the frequency domain */
  INT4 printfirflag;          /* flag set to 1 if user wants to print the filter in the time domain */
  INT4 printsnrflag;          /* flag set to 1 if user wants to print the snr */
  INT4 printdataflag;         /* flag set to 1 if user wants to print the data */  
  INT4 printinjectionflag;    /* flag set to 1 if user wants to print the injection(s) */  
} CommandLineArgs;

typedef 
struct GlobalVariablesTag {
  INT4 duration;              /* duration of entire segment to be analysed */
  LIGOTimeGPS gpsepoch;       /* GPS epoch of start of entire segment to be analysed */ 
  REAL8TimeSeries *ht;        /* raw input data (LIGO data) */
  REAL4TimeSeries *ht_V;      /* raw input data (Virgo data) */
  REAL4TimeSeries *ht_proc;   /* processed (band-pass filtered and down-sampled) input data */
  REAL4FrequencySeries Spec;  /* average spectrum */
  RealFFTPlan *fplan;         /* fft plans */
  RealFFTPlan *rplan;         /* fft plans */
  INT4 seg_length;
} GlobalVariables;

typedef
struct StringTemplateTag {
  INT4 findex;                /* Template frequency index */
  REAL4 f;                    /* Template frequency */
  REAL4 norm;                 /* Template normalisation */
  REAL4 mismatch;             /* Template mismatch relative to last one*/
  REAL4FrequencySeries StringFilter; /* Frequency domain filter corresponding to this template */
} StringTemplate;

/***************************************************************************/

/* GLOBAL VARIABLES */
INT4 lalDebugLevel=3;
FrCache *framecache;          /* frame reading variables */
FrStream *framestream=NULL;

GlobalVariables GV;           /* A bunch of stuff is stored in here; mainly to protect it from accidents */

StringTemplate strtemplate[MAXTEMPLATES];

int NTemplates;

SnglBurst *events=NULL;
MetadataTable  procTable;
MetadataTable  procparams;
MetadataTable  searchsumm;

CHAR outfilename[256];
CHAR ifo[4]; 

REAL4 SAMPLERATE;

int Nevents=0;

PassBandParamStruc highpassParams;



/***************************************************************************/

/* FUNCTION PROTOTYPES */

/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Reads raw data (or puts in fake gaussian noise with a sigma=10^-20) */
int ReadData(struct CommandLineArgsTag CLA);

/* Adds injections if an xml injection file is given */
int AddInjections(struct CommandLineArgsTag CLA);

/* windows the data with a Tukey window */
int WindowData(void);

/* High pass filters and casts data to REAL4 */
int ProcessData(void);

/* DownSamples data */
int DownSample(struct CommandLineArgsTag CLA);

/* Computes the average spectrum  */
int AvgSpectrum(struct CommandLineArgsTag CLA);

/* Creates the template bank based on the spectrum  */
int CreateTemplateBank(struct CommandLineArgsTag CLA);

/* Creates the frequency domain string cusp or kink filters  */
int CreateStringFilters(struct CommandLineArgsTag CLA);

/* Filters the data through the template banks  */
int FindStringBurst(struct CommandLineArgsTag CLA);

/* Finds events above SNR threshold specified  */
int FindEvents(struct CommandLineArgsTag CLA, REAL4Vector *vector, 
	       INT4 i, INT4 m, SnglBurst **thisEvent);

/* Writes out the xml file with the events it found  */
int OutputEvents(struct CommandLineArgsTag CLA);

/* Frees the memory */
int FreeMem(void);                                        

/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  /****** ReadCommandLine ******/
  printf("ReadCommandLine()\n");
  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
  
  highpassParams.nMax =  4;
  highpassParams.f1   = -1;
  highpassParams.a1   = -1;
  highpassParams.f2   = CommandLineArgs.flow;
  highpassParams.a2   = 0.9; /* this means 90% of amplitude at f2 */
  
  /****** ReadData ******/
  printf("ReadData()\n");
  if (ReadData(CommandLineArgs)) return 2;
  
  /****** AddInjections ******/
  if (CommandLineArgs.InjectionFile != NULL) {
    printf("AddInjections()\n");
    if (AddInjections(CommandLineArgs)) return 3;
    /* at this stage, ht_proc contains only the injection */
    if ( CommandLineArgs.printinjectionflag ) LALSPrintTimeSeries( GV.ht_proc, "injection.txt" );
  }
  
  /****** WindowData ******/
  printf("WindowData()\n");
  if (WindowData()) return 4;
  
  /****** ProcessData ******/
  printf("ProcessData()\n");
  if (ProcessData()) return 5;
  
  if ( CommandLineArgs.printdataflag ){
    int p;
    for (p=0; p<(int)GV.ht_proc->data->length; p++)
      fprintf(stdout,"%1.15e\n",GV.ht_proc->data->data[p]);
    return 0;
  }
  
  /****** DownSample ******/
  printf("DownSample()\n");
  if (DownSample(CommandLineArgs)) return 6;
  
  /****** XLALResizeREAL4TimeSeries ******/
  printf("XLALResizeREAL4TimeSeries()\n");	
  /* re-size the time series to remove the pad */
  GV.ht_proc = XLALResizeREAL4TimeSeries(GV.ht_proc, 
					 (int)(CommandLineArgs.pad/GV.ht_proc->deltaT+0.5),
					 GV.ht_proc->data->length-2*(UINT4)(CommandLineArgs.pad/GV.ht_proc->deltaT+0.5));
    
  /* reduce duration of segment appropriately */
  GV.duration -= 2*CommandLineArgs.pad; 
  
  /****** AvgSpectrum ******/
  printf("AvgSpectrum()\n");
  if (AvgSpectrum(CommandLineArgs)) return 7;  
  if (CommandLineArgs.printspectrumflag) LALSPrintFrequencySeries( &(GV.Spec), "Spectrum.txt" );
  
  /****** CreateTemplateBank ******/
  printf("CreateTemplateBank()\n");
  if (CreateTemplateBank(CommandLineArgs)) return 8;
  
  /****** CreateStringFilters ******/
  printf("CreateStringFilters()\n");
  if (CreateStringFilters(CommandLineArgs)) return 9;
  
  /****** FindStringBurst ******/
  printf("FindStringBurst()\n");
  if (FindStringBurst(CommandLineArgs)) return 10;
  
  /****** XLALClusterSnglBurstTable ******/
  printf("XLALClusterSnglBurstTable()\n");
  if (CommandLineArgs.cluster != 0.0 && events){
    XLALClusterSnglBurstTable(&events, XLALCompareStringBurstByTime, XLALCompareStringBurstByTime, XLALStringBurstCluster);
    XLALSortSnglBurst(&events, XLALCompareSnglBurstByPeakTimeAndSNR);
  }
  
  /****** XLALSnglBurstAssignIDs ******/
  printf("XLALSnglBurstAssignIDs()\n");
  XLALSnglBurstAssignIDs(events, procTable.processTable->process_id, 0);
  
  /****** OutputEvents ******/
  printf("OutputEvents()\n");
  if (OutputEvents(CommandLineArgs)) return 12;
  
  /****** FreeMem ******/
  printf("FreeMem()\n");
  if (FreeMem()) return 13;
  
  return 0;
}


/************************************* MAIN PROGRAM ENDS *************************************/


/*******************************************************************************/

int WindowData(void){
  REAL8Window *window;
  REAL8 r = 0.001;
  int k;
  int N=GV.ht->data->length;
  
  window = XLALCreateTukeyREAL8Window(N, r);
    
  for (k = 0; k < N; k++)
    GV.ht->data->data[k] *= window->data->data[k];
      
  XLALDestroyREAL8Window(window);
  
  return 0;
}

/*******************************************************************************/

int AddInjections(struct CommandLineArgsTag CLA){
  REAL8TimeSeries *injections;
  LIGOTimeGPS startTime = GV.ht->epoch;
  LIGOTimeGPS stopTime = GV.ht->epoch;
  SimBurst *sim_burst;
  unsigned p;
  
  XLALGPSAdd(&startTime, CLA.ShortSegDuration / 4 + CLA.pad);
  XLALGPSAdd(&stopTime, GV.ht->data->length * GV.ht->deltaT - CLA.ShortSegDuration / 4 - CLA.pad);

  /* Get info from injection file */
  sim_burst = XLALSimBurstTableFromLIGOLw(CLA.InjectionFile, &startTime, &stopTime);
  
  /* new injection code is double precision, so we need to create a
   * buffer to put the injections in and then quantize to single precision
   * for the string code */
  injections = XLALCreateREAL8TimeSeries(GV.ht_proc->name, &GV.ht_proc->epoch, GV.ht_proc->f0, GV.ht_proc->deltaT, &GV.ht_proc->sampleUnits, GV.ht_proc->data->length);
  memset(injections->data->data, 0, injections->data->length * sizeof(*injections->data->data));

  /* Inject the signals into ht_proc -> for printing
     Inject the signals into ht */
  if(XLALBurstInjectSignals(injections, sim_burst, NULL)) return 1;
  for(p = 0; p < GV.ht_proc->data->length; p++){
    GV.ht_proc->data->data[p] += (REAL4)injections->data->data[p];
    GV.ht->data->data[p] += injections->data->data[p];
  }

  /* injection time series is no more needed */
  XLALDestroyREAL8TimeSeries(injections);

  /* free the injection table */
  while(sim_burst) {
    SimBurst *next = sim_burst->next;
    XLALDestroySimBurst(sim_burst);
    sim_burst = next;
  }

  return 0;
}

/*******************************************************************************/

static ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param,
					      const ProcessTable *process,
					      const char *type, const char *param, const char *value){
  *proc_param = XLALCreateProcessParamsTableRow(process);
  snprintf((*proc_param)->program, LIGOMETA_PROGRAM_MAX, PROGRAM_NAME);
  snprintf((*proc_param)->type, LIGOMETA_TYPE_MAX, type);
  snprintf((*proc_param)->param, LIGOMETA_PARAM_MAX, "--%s", param);
  snprintf((*proc_param)->value, LIGOMETA_VALUE_MAX, value);
  
  return(&(*proc_param)->next);
}

/*******************************************************************************/

int OutputEvents(struct CommandLineArgsTag CLA){  
  LIGOLwXMLStream *xml;
  
  if (!CLA.outputFileName){
    snprintf(outfilename, sizeof(outfilename)-1, "%s-STRINGSEARCH-%d-%d.xml", ifo,
	     searchsumm.searchSummaryTable->in_start_time.gpsSeconds,
	     searchsumm.searchSummaryTable->in_end_time.gpsSeconds - 
	     searchsumm.searchSummaryTable->in_start_time.gpsSeconds);
    outfilename[sizeof(outfilename)-1] = '\0';
  }
  else{
    snprintf(outfilename, sizeof(outfilename)-1, "%s", CLA.outputFileName);
    outfilename[sizeof(outfilename)-1] = '\0';
  }

  xml = XLALOpenLIGOLwXMLFile(outfilename);

  /* process table */
  snprintf(procTable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo);
  XLALGPSTimeNow(&(procTable.processTable->start_time));
  
  if(XLALWriteLIGOLwXMLProcessTable(xml, procTable.processTable)) return -1;
  
  /* process params table */
  if(XLALWriteLIGOLwXMLProcessParamsTable(xml, procparams.processParamsTable)) return -1;
  
  /* search summary table */
  snprintf(searchsumm.searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo);
  searchsumm.searchSummaryTable->nevents = Nevents;
  
  if(XLALWriteLIGOLwXMLSearchSummaryTable(xml, searchsumm.searchSummaryTable)) return -1;

  /* burst table */
  if(XLALWriteLIGOLwXMLSnglBurstTable(xml, events)) return -1;
  
  XLALCloseLIGOLwXMLFile(xml);
  
  /* free event list, process table, search summary and process params */

  XLALDestroySnglBurstTable(events);
  XLALDestroyProcessTable(procTable.processTable);
  XLALDestroySearchSummaryTable(searchsumm.searchSummaryTable);
  XLALDestroyProcessParamsTable(procparams.processParamsTable);

  return 0;
}

/*******************************************************************************/

int FindEvents(struct CommandLineArgsTag CLA, REAL4Vector *vector, INT4 i, INT4 m, SnglBurst **thisEvent){
  int p;
  REAL4 maximum;
  REAL8 duration;
  INT4 pmax, pend, pstart;
  INT8  peaktime, starttime;
  INT8  timeNS;

  /* print the snr to stdout */
  if (CLA.printsnrflag)
    for ( p = (int)vector->length/4 ; p < (int)(3*vector->length/4); p++ )
      fprintf(stdout,"%d %e\n",m, vector->data[p]);
	

  /* Now find thisEvent in the inner half */
  for ( p = (int)vector->length/4 ; p < (int)(3*vector->length/4); p++ ){
    maximum = 0.0;
    pmax=p;
    timeNS  = (INT8)( 1000000000 ) * (INT8)(GV.ht_proc->epoch.gpsSeconds+GV.seg_length*i/2*GV.ht_proc->deltaT);
    timeNS  +=   (INT8)( 1e9 * GV.ht_proc->deltaT * p );

    /* Do we have the start of a cluster? */
    if ( (fabs(vector->data[p]) > CLA.threshold) && ( (double)(1e-9*timeNS) > (double)CLA.trigstarttime)){
      pend=p; pstart=p;
      
      timeNS  = (INT8)( 1000000000 ) * (INT8)(GV.ht_proc->epoch.gpsSeconds+GV.seg_length*i/2*GV.ht_proc->deltaT);
      
      if ( *thisEvent ){ /* create a new event */
	(*thisEvent)->next = XLALCreateSnglBurst();
	*thisEvent = (*thisEvent)->next;
      }
      else /* create the list */
	*thisEvent = events = XLALCreateSnglBurst();
            
      if ( ! *thisEvent ){ /* allocation error */
	fprintf(stderr,"Could not allocate memory for event. Memory allocation error. Exiting. \n");
	return 1;
      }

      /* Clustering in time: While we are above threshold, or within clustering time of the last point above threshold... */
      while( ((fabs(vector->data[p]) > CLA.threshold) || ((p-pend)* GV.ht_proc->deltaT < (float)(CLA.cluster)) ) 
	     && p<(int)(3*vector->length/4)){
	
	/* This keeps track of the largest SNR point of the cluster */
	if(fabs(vector->data[p]) > maximum){
	  maximum=fabs(vector->data[p]);
	  pmax=p;
	}
	/* pend is the last point above threshold */
	if ( (fabs(vector->data[p]) > CLA.threshold))
	  pend =  p; 
	
	p++;
      }

      peaktime = timeNS + (INT8)( 1e9 * GV.ht_proc->deltaT * pmax );
      duration = GV.ht_proc->deltaT * ( pend - pstart );

      starttime = timeNS + (INT8)( 1e9 * GV.ht_proc->deltaT * pstart );

      /* Now copy stuff into event */
      strncpy( (*thisEvent)->ifo, CLA.ChannelName, sizeof(ifo)-2 );
      strncpy( (*thisEvent)->search, "StringCusp", sizeof( (*thisEvent)->search ) );
      strncpy( (*thisEvent)->channel, CLA.ChannelName, sizeof( (*thisEvent)->channel ) );
      
      /* give trigger a 1 sample fuzz on either side */
      starttime -= GV.ht_proc->deltaT *1e9;
      duration += 2*GV.ht_proc->deltaT;
      
      (*thisEvent)->start_time.gpsSeconds     = starttime / 1000000000;
      (*thisEvent)->start_time.gpsNanoSeconds = starttime % 1000000000;
      (*thisEvent)->peak_time.gpsSeconds      = peaktime / 1000000000;
      (*thisEvent)->peak_time.gpsNanoSeconds  = peaktime % 1000000000;
      (*thisEvent)->duration     = duration;
      (*thisEvent)->central_freq = (strtemplate[m].f+CLA.fbankstart)/2.0;	   
      (*thisEvent)->bandwidth    = strtemplate[m].f-CLA.fbankstart;				     
      (*thisEvent)->snr          = maximum;
      (*thisEvent)->amplitude   = vector->data[pmax]/strtemplate[m].norm;
      (*thisEvent)->confidence   = -fabs((*thisEvent)->amplitude); /* FIXME */
      (*thisEvent)->string_cluster_t = CLA.cluster;
    }
  }
  
  return 0;
}

/*******************************************************************************/

int FindStringBurst(struct CommandLineArgsTag CLA){
  int i,p,m; 
  REAL4Vector *vector = NULL;
  COMPLEX8Vector *vtilde = NULL;
  SnglBurst *thisEvent = NULL;

  /* create vector that will hold the data for each overlapping chunk */ 
  vector = XLALCreateREAL4Vector( GV.seg_length);

  /* create vector that will hold FFT of data*/
  vtilde = XLALCreateCOMPLEX8Vector( GV.seg_length / 2 + 1 );

  /* loop over templates  */
  for (m = 0; m < NTemplates; m++){
    /* loop over overlapping chunks */ 
    for(i=0; i < 2*GV.duration/CLA.ShortSegDuration - 1 ;i++){
      /* populate vector that will hold the data for each overlapping chunk */
      memcpy( vector->data, GV.ht_proc->data->data + i*GV.seg_length/2,vector->length*sizeof( *vector->data ) );
	  
      /* fft it */
      if(XLALREAL4ForwardFFT( vtilde, vector, GV.fplan )) return 1;
      
      /* multiply FT of data and String Filter and deltaT (latter not included in LALForwardRealFFT) */
      for ( p = 0 ; p < (int) vtilde->length; p++ ){
	vtilde->data[p].re *= strtemplate[m].StringFilter.data->data[p]*GV.ht_proc->deltaT;
	vtilde->data[p].im *= strtemplate[m].StringFilter.data->data[p]*GV.ht_proc->deltaT;
      }
      
      if(XLALREAL4ReverseFFT( vector, vtilde, GV.rplan )) return 1;
      
      /* normalise the result by template normalisation and multiply by 
	 df (not inluded in LALReverseRealFFT)  factor of 2 is from 
	 match-filter definition */
      
      for ( p = 0 ; p < (int)vector->length; p++ )
	vector->data[p] *= 2.0 * GV.Spec.deltaF / strtemplate[m].norm;
      	
      if(FindEvents(CLA, vector, i, m, &thisEvent)) return 1;
    }
  }

  /* sort events in time; if there are any */
  if (events) /* first sort list in increasing GPS peak time */
    XLALSortSnglBurst(&events, XLALCompareSnglBurstByPeakTimeAndSNR);
  
  XLALDestroyCOMPLEX8Vector( vtilde );
  XLALDestroyREAL4Vector( vector );

  return 0;
}


/*******************************************************************************/

int CreateStringFilters(struct CommandLineArgsTag CLA){

  int p, m, f_low_cutoff_index, f_high_cutoff_index; 
  COMPLEX8Vector *vtilde; /* frequency-domain vector workspace */
  REAL4Vector    *vector; /* time-domain vector workspace */
  REAL4 f, re, im;
  REAL4TimeSeries series;
  CHAR filterfilename[256];

  vector = XLALCreateREAL4Vector( GV.seg_length);
  vtilde = XLALCreateCOMPLEX8Vector( GV.seg_length / 2 + 1 );
 
  f_low_cutoff_index = (int) (CLA.fbankstart/ GV.Spec.deltaF+0.5);

  for (m = 0; m < NTemplates; m++){
    f_high_cutoff_index = (int) (strtemplate[m].f/ GV.Spec.deltaF+0.5);

    /* create the space for the filter */
    strtemplate[m].StringFilter.deltaF=GV.Spec.deltaF;
    strtemplate[m].StringFilter.data = XLALCreateREAL4Vector(GV.Spec.data->length);
            
    /* populate vtilde with the template divided by the noise */
    for ( p = f_low_cutoff_index; p < (int) vtilde->length; p++ ){
      f=p*GV.Spec.deltaF;
	  
      if(f<=strtemplate[m].f) vtilde->data[p].re = sqrt(pow(f,CLA.power)/(GV.Spec.data->data[p]));
      else vtilde->data[p].re = sqrt(pow(f,CLA.power)*exp(1-f/strtemplate[m].f)/(GV.Spec.data->data[p]));
      vtilde->data[p].im = 0;
    }
      
    /* set all frequencies below the low freq cutoff to zero */
    memset( vtilde->data, 0, f_low_cutoff_index  * sizeof( *vtilde->data ) );
    
    /* set DC and Nyquist to zero anyway */
    vtilde->data[0].re = vtilde->data[vtilde->length - 1].re = 0;
    vtilde->data[0].im = vtilde->data[vtilde->length - 1].im = 0;

    /* reverse FFT vtilde into vector */
    if(XLALREAL4ReverseFFT( vector, vtilde, GV.rplan )) return 1;
             
    /* multiply times df to make sure units are correct */
    for ( p = 0 ; p < (int)vector->length; p++ )
      vector->data[p] *= GV.Spec.deltaF;

    /* perform the truncation; the truncation is CLA.TruncSecs/2 because 
       we are dealing with the sqrt of the filter at the moment*/
    if(CLA.TruncSecs != 0.0)
      memset( vector->data + (INT4)(CLA.TruncSecs/2/GV.ht_proc->deltaT +0.5), 0,
	      ( vector->length -  2 * (INT4)(CLA.TruncSecs/2/GV.ht_proc->deltaT +0.5)) 
	      * sizeof( *vector->data ) );
    
    
    /* forward fft the truncated vector into vtilde */
    if(XLALREAL4ForwardFFT( vtilde, vector, GV.fplan )) return 1;
    
    for ( p = 0 ; p < (int)vtilde->length-1; p++ ){
      re = vtilde->data[p].re * GV.ht_proc->deltaT;
      im = vtilde->data[p].im * GV.ht_proc->deltaT;
      strtemplate[m].StringFilter.data->data[p] = (re * re + im * im);
    }

    /* set DC and Nyquist to 0*/
    strtemplate[m].StringFilter.data->data[0] =
      strtemplate[m].StringFilter.data->data[vtilde->length-1] = 0;
        
    /* print out the frequency domain filter */
    if (CLA.printfilterflag){
      snprintf(filterfilename, sizeof(filterfilename)-1, "Filter-%d.txt", m);
      filterfilename[sizeof(outfilename)-1] = '\0';
      LALSPrintFrequencySeries( &(strtemplate[m].StringFilter), filterfilename );
    }

    /* print out the time domain FIR filter */
    if (CLA.printfirflag){
      series.deltaT=GV.ht_proc->deltaT;
      series.f0 = 0.0;
      strncpy(series.name, "fir filter", LALNameLength);
      series.epoch=GV.ht_proc->epoch;
      series.sampleUnits=GV.ht_proc->sampleUnits;
      
      for ( p = 0 ; p < (int)vtilde->length-1; p++ ){
	re = vtilde->data[p].re * GV.ht_proc->deltaT;
	im = vtilde->data[p].im * GV.ht_proc->deltaT;
	
	vtilde->data[p].re = (re * re + im * im);
	vtilde->data[p].im = 0.0;
      }
      vtilde->data[0].re = vtilde->data[0].im = 0.0;
      vtilde->data[vtilde->length-1].re = vtilde->data[vtilde->length-1].im =0;
      
      if(XLALREAL4ReverseFFT( vector, vtilde, GV.rplan )) return 1;
      
      series.data = vector;
      
      snprintf(filterfilename, sizeof(filterfilename)-1, "FIRFilter-%d.txt", m);
      filterfilename[sizeof(outfilename)-1] = '\0';
      LALSPrintTimeSeries( &series, filterfilename );
    }
  }
  
  XLALDestroyCOMPLEX8Vector( vtilde );
  XLALDestroyREAL4Vector( vector );
  
  return 0;
}

/*******************************************************************************/

int CreateTemplateBank(struct CommandLineArgsTag CLA){
  REAL8 fmax, f_cut, f, t1t1, t2t2, t1t2, epsilon, previous_epsilon;
  int p, pcut, f_min_index, f_max_index, f_cut_index, k;
  REAL4Vector *integral;

  fmax = (1.0/GV.ht_proc->deltaT) / 2.0;
  f_min_index = CLA.fbankstart / GV.Spec.deltaF;
  f_max_index = fmax / GV.Spec.deltaF;
  integral = XLALCreateREAL4Vector(f_max_index-f_min_index);

  /* first template : f_cutoff = fbankhighfcutofflow */
  f_cut = (int)(CLA.fbankhighfcutofflow/GV.Spec.deltaF)*GV.Spec.deltaF;
  f_cut_index = CLA.fbankhighfcutofflow / GV.Spec.deltaF;

  /* compute (t1|t1) */
  t1t1=0.0;
  integral->data[0]=4*pow( pow(CLA.fbankstart,CLA.power),2)/GV.Spec.data->data[f_min_index]*GV.Spec.deltaF;
  for( p = f_min_index ; p < f_max_index; p++ ){
    f = p*GV.Spec.deltaF;
    
    if(f<=f_cut) t1t1 += 4*pow(pow(f,CLA.power),2)/GV.Spec.data->data[p]*GV.Spec.deltaF;
    else t1t1 += 4*pow( pow(f,CLA.power)*exp(1-f/f_cut) ,2)/GV.Spec.data->data[p]*GV.Spec.deltaF;

    if(p>0) /* keep the integral in memory (to run faster) */
      integral->data[p-f_min_index] = integral->data[p-f_min_index-1]+4*pow(pow(f,CLA.power),2)/GV.Spec.data->data[p]*GV.Spec.deltaF;
  }
  
  strtemplate[0].findex=f_cut_index;
  strtemplate[0].f=f_cut;
  strtemplate[0].mismatch=0.0;
  strtemplate[0].norm=sqrt(t1t1);
  k=1;
  fprintf(stdout,"%% Templ. frequency      sigma      mismatch\n");  
  fprintf(stdout,"%% %d      %1.3e    %1.3e    %1.3e\n",k-1,strtemplate[0].f,strtemplate[0].norm, strtemplate[0].mismatch);
  
  /* find the next cutoffs given the maximal mismatch */
  for(pcut=f_cut_index+1; pcut<f_max_index; pcut++){
    f_cut = pcut*GV.Spec.deltaF;
   
    t2t2=integral->data[strtemplate[k-1].findex-f_min_index];
    t1t2=integral->data[strtemplate[k-1].findex-f_min_index];
    
    /* compute (t2|t2) and (t1|t2) */
    for( p = strtemplate[k-1].findex+1 ; p < f_max_index; p++ ){
      f = p*GV.Spec.deltaF;
      
      /* (t2|t2) */
      if(f<=f_cut)
	t2t2 += 4*pow(pow(f,CLA.power),2)/GV.Spec.data->data[p]*GV.Spec.deltaF;
      else 
	t2t2 += 4*pow( pow(f,CLA.power)*exp(1-f/f_cut) ,2)/GV.Spec.data->data[p]*GV.Spec.deltaF;

      /* (t1|t2) */
      if(f<=f_cut)
	t1t2 += 4*pow(pow(f,CLA.power),2)*exp(1-f/strtemplate[k-1].f) /GV.Spec.data->data[p]*GV.Spec.deltaF;
      else 
	t1t2 += 4*pow( pow(f,CLA.power),2)*exp(1-f/strtemplate[k-1].f)*exp(1-f/f_cut) /GV.Spec.data->data[p]*GV.Spec.deltaF;
    }
        
    previous_epsilon = epsilon;
    epsilon=1-t1t2/sqrt(t1t1*t2t2);

    /*if(pcut%50==0) printf("%d %f %f\n",pcut, f_cut, epsilon);*/

    if(epsilon >= CLA.fmismatchmax || pcut==f_max_index-1){
      strtemplate[k].findex=pcut;
      strtemplate[k].f=f_cut;
      strtemplate[k].norm=sqrt(t2t2);
      strtemplate[k].mismatch=epsilon;
      k++;
      fprintf(stdout,"%% %d      %1.3e    %1.3e    %1.3e\n",k-1,strtemplate[k-1].f,strtemplate[k-1].norm, strtemplate[k-1].mismatch);
      t1t1=t2t2;
      if(k == MAXTEMPLATES){
	fprintf(stderr,"Too many templates for code... Exiting\n");
	return 1;
      }
    }
    
    /* to get faster (not so smart though) */
    if(pcut<f_max_index-16 && (epsilon-previous_epsilon)<0.005) 
      pcut+=15;

  }
 
  NTemplates=k;
  XLALDestroyREAL4Vector( integral );

  return 0;
}


/*******************************************************************************/
int AvgSpectrum(struct CommandLineArgsTag CLA){
  
  int p;
  int segmentLength;
  int segmentStride;
  REAL4Window  *window4;

  GV.seg_length = (int)(CLA.ShortSegDuration/GV.ht_proc->deltaT + 0.5);
  GV.Spec.data = XLALCreateREAL4Vector(GV.seg_length / 2 + 1);
  GV.fplan = XLALCreateForwardREAL4FFTPlan( GV.seg_length, 0 );
  GV.rplan = XLALCreateReverseREAL4FFTPlan( GV.seg_length, 0 );

  if (CLA.fakenoiseflag && CLA.whitespectrumflag){
    for ( p = 0 ; p < (int)GV.Spec.data->length; p++ )
      GV.Spec.data->data[p]=2/SAMPLERATE;
    GV.Spec.deltaF=1/(GV.seg_length*GV.ht_proc->deltaT);
  }
  else{
    segmentLength = GV.seg_length;
    segmentStride = GV.seg_length/2;
    window4  = NULL;

    window4 = XLALCreateHannREAL4Window( segmentLength );
    if(XLALREAL4AverageSpectrumMedianMean( &GV.Spec, GV.ht_proc, segmentLength,
					   segmentStride, window4, GV.fplan ))
	return 1;

      XLALDestroyREAL4Window( window4 );
    }

  return 0;
}

/*******************************************************************************/

int DownSample(struct CommandLineArgsTag CLA){
  if(XLALResampleREAL4TimeSeries(GV.ht_proc, 1.0/CLA.samplerate)) return 1;
  return 0;
}

/*******************************************************************************/

int ProcessData(void){
  int p;

  if(XLALButterworthREAL8TimeSeries(GV.ht, &highpassParams)) return 1;
  for (p=0; p<(int)GV.ht->data->length; p++)  
    GV.ht_proc->data->data[p]=(REAL4)GV.ht->data->data[p]; 
    
  /* ht is no more needed -> free memory */
  XLALDestroyREAL8TimeSeries(GV.ht);

  return 0;
}

/*******************************************************************************/

int ReadData(struct CommandLineArgsTag CLA){
  int p;
  
  /* create Frame cache, open frame stream and delete frame cache */
  framecache = XLALFrImportCache(CLA.FrCacheFile);
  framestream = XLALFrCacheOpen(framecache);
  XLALFrDestroyCache(framecache);
  
  GV.duration                = CLA.GPSEnd-CLA.GPSStart;
  GV.gpsepoch.gpsSeconds     = CLA.GPSStart;
  GV.gpsepoch.gpsNanoSeconds = 0;
  
  /* Double vs. simple precision data for LIGO vs. Virgo */
  if(CLA.ChannelName[0]=='V'){

    /* create and initialize _simple_ precision time series */
    GV.ht_V  = XLALCreateREAL4TimeSeries(CLA.ChannelName, &GV.gpsepoch, 0, 0, &lalStrainUnit, 1);

    /* get the meta data */
    XLALFrGetREAL4TimeSeriesMetadata(GV.ht_V,framestream);

    /* resize ht to the correct number of samples */
    XLALResizeREAL4TimeSeries(GV.ht_V, 0, (UINT4)(GV.duration/GV.ht_V->deltaT +0.5));

  }
  else{

    /* create and initialize _double_ precision time series */
    GV.ht  = XLALCreateREAL8TimeSeries(CLA.ChannelName, &GV.gpsepoch, 0, 0, &lalStrainUnit, 1);

    /* get the meta data */
    XLALFrGetREAL8TimeSeriesMetadata(GV.ht,framestream);

    /* resize ht to the correct number of samples */
    XLALResizeREAL8TimeSeries(GV.ht, 0, (UINT4)(GV.duration/GV.ht->deltaT +0.5));
    
  }  


  /* If we are reading real noise then read it*/
  if(!CLA.fakenoiseflag){
    /* seek to and read data */
    XLALFrSeek( framestream, &GV.gpsepoch );
    
    if(CLA.ChannelName[0]=='V'){
      XLALFrGetREAL4TimeSeries(GV.ht_V,framestream);
      
      /* Allocate space for REAL8 data */
      GV.ht  = XLALCreateREAL8TimeSeries(GV.ht_V->name, 
					 &GV.ht_V->epoch, 
					 GV.ht_V->f0, 
					 GV.ht_V->deltaT, 
					 &lalStrainUnit, 
					 (UINT4)(GV.duration/GV.ht_V->deltaT +0.5));
	
      /* Fill REAL8 data vector */
      for (p=0; p<(int)GV.ht_V->data->length; p++)
	GV.ht->data->data[p] = (REAL8)GV.ht_V->data->data[p];
    }
    else XLALFrGetREAL8TimeSeries(GV.ht,framestream);

    /* Scale data to avoid single float precision problems */
    for (p=0; p<(int)GV.ht->data->length; p++)
      GV.ht->data->data[p] *= SCALE;
  }
  /* otherwise create random data set */
  else{
    FILE *devrandom;
    RandomParams   *randpar=NULL;
    REAL4Vector    *v1=NULL;
    int seed, errorcode;
    
    if(!(devrandom=fopen("/dev/urandom","r"))){
      fprintf(stderr,"Unable to open device /dev/urandom\n");
      return 1;
    }
    errorcode=fread((void*)&seed,sizeof(INT4),1,devrandom);
    if (errorcode!=1){
      fprintf( stderr,"Error reading /dev/urandom file!\n");
      return 1;
    }
    fclose(devrandom);
    
    v1 = XLALCreateREAL4Vector(GV.ht->data->length);
    
    randpar = XLALCreateRandomParams (seed);
    if(XLALNormalDeviates(v1, randpar)) return 1;;
    XLALDestroyRandomParams (randpar);
    
    for (p=0; p<(int)GV.ht->data->length; p++)
      GV.ht->data->data[p] = v1->data[p];
    
    XLALDestroyREAL4Vector(v1);
  }
  
  SAMPLERATE=1.0/GV.ht->deltaT;

  /* Allocate space for processed data */
  GV.ht_proc  = XLALCreateREAL4TimeSeries(GV.ht->name, 
					  &GV.ht->epoch, 
					  GV.ht->f0, 
					  GV.ht->deltaT, 
					  &lalStrainUnit, 
					  (UINT4)(GV.duration/GV.ht->deltaT +0.5));
  /* zero out processed data */
  for (p=0; p<(int)GV.ht_proc->data->length; p++) GV.ht_proc->data->data[p] = 0.0;
  
  XLALFrClose(framestream);
  if(CLA.ChannelName[0]=='V') XLALDestroyREAL4TimeSeries(GV.ht_V);
  
  return 0;
}




/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA){
  INT4 errflg = 0;
  ProcessParamsTable **paramaddpoint = &procparams.processParamsTable;
  struct option long_options[] = {
    {"bw-flow",                     required_argument, NULL,           'f'},
    {"bank-freq-start",             required_argument, NULL,           'L'}, 
    {"bank-lowest-hifreq-cutoff",   required_argument, NULL,           'H'},
    {"max-mismatch",                required_argument, NULL,           'M'},
    {"threshold",                   required_argument, NULL,           't'},
    {"frame-cache",                 required_argument, NULL,           'F'},
    {"channel-name",                required_argument, NULL,           'C'},
    {"outfile",                     required_argument, NULL,           'o'},
    {"gps-end-time",                required_argument, NULL,           'E'},
    {"gps-start-time",              required_argument, NULL,           'S'},
    {"injection-file",              required_argument, NULL,           'i'},
    {"short-segment-duration",      required_argument, NULL,           'd'},
    {"settling-time",               required_argument, NULL,           'T'},
    {"sample-rate",                 required_argument, NULL,           's'},
    {"trig-start-time",             required_argument, NULL,           'g'},
    {"pad",                         required_argument, NULL,           'p'},
    {"cusp-search",                 no_argument, NULL,          'c' },
    {"kink-search",                 no_argument, NULL,          'k' },
    {"test-gaussian-data",          no_argument, NULL,          'n' },
    {"test-white-spectrum",         no_argument, NULL,          'w' },
    {"cluster-events",              required_argument, NULL,          'l' },
    {"print-spectrum",              no_argument, NULL,          'a' },
    {"print-fd-filter",             no_argument, NULL,          'b' },    
    {"print-snr",                   no_argument, NULL,          'r' },        
    {"print-td-filter",             no_argument, NULL,          'x' },        
    {"print-data",                  no_argument, NULL,          'y' },        
    {"print-injection",             no_argument, NULL,          'z' },        
    {"help",                        no_argument, NULL,          'h' },
    {0, 0, 0, 0}
  };
  char args[] = "hnckwabrxyzl:f:L:M:H:t:F:C:E:S:i:d:T:s:g:o:p:";

  optarg = NULL;
  /* set up xml output stuff */
  /* create the process and process params tables */
  procTable.processTable = XLALCreateProcessTableRow();
  XLALGPSTimeNow(&(procTable.processTable->start_time));
  if (strcmp(CVS_REVISION, "$Revi" "sion$"))
    {
      if(XLALPopulateProcessTable(procTable.processTable, PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE, 0))
	exit(1);
    }
  else
    {
      if(XLALPopulateProcessTable(procTable.processTable, PROGRAM_NAME, lalappsGitCommitID, lalappsGitGitStatus, lalappsGitCommitDate, 0))
	exit(1);
    }
  procparams.processParamsTable = NULL;
  /* create the search summary table */
  searchsumm.searchSummaryTable = XLALCreateSearchSummaryTableRow(procTable.processTable);
  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;

  /* Initialize default values */
  CLA->flow=0.0;
  CLA->fbankstart=0.0;
  CLA->fbankhighfcutofflow=0.0;
  CLA->fmismatchmax=0.05;
  CLA->FrCacheFile=NULL;
  CLA->InjectionFile=NULL;
  CLA->ChannelName=NULL;
  CLA->outputFileName=NULL;
  CLA->GPSStart=0;
  CLA->GPSEnd=0;
  CLA->ShortSegDuration=0;
  CLA->TruncSecs=0;
  CLA->power=0.0;
  CLA->threshold=0.0;
  CLA->fakenoiseflag=0;
  CLA->whitespectrumflag=0;
  CLA->samplerate=4096.0;
  CLA->trigstarttime=0;
  CLA->cluster=0.0;
  CLA->pad=0;
  CLA->printfilterflag=0;
  CLA->printspectrumflag=0;
  CLA->printsnrflag=0;
  CLA->printfirflag=0;
  CLA->printdataflag=0;
  CLA->printinjectionflag=0;
  
  /* initialise ifo string */
  memset(ifo, 0, sizeof(ifo));

  /* Scan through list of command line arguments */
  while ( 1 )
  {
    int option_index = 0; /* getopt_long stores long option here */
    int c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch ( c )
    {

    case 'f':
      /* low frequency cutoff */
      CLA->flow=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "float");
      break;
    case 's':
      /* low frequency cutoff */
      CLA->samplerate=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "float");
      break;
    case 'H':
      /* lowest high frequency cutoff */
      CLA->fbankhighfcutofflow=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "float");
      break;
    case 'M':
      /* Maximal mismatch */
      CLA->fmismatchmax=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "float");
      break;
    case 'L':
      /* low frequency cutoff */
      CLA->fbankstart=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "float");
      break;
    case 't':
      /* low frequency cutoff */
      CLA->threshold=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "float");
      break;
    case 'F':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'C':
      /* name channel */
      CLA->ChannelName=optarg;
      memcpy(ifo, optarg, sizeof(ifo) - 2);
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'i':
      /* name of xml injection file */
      CLA->InjectionFile=optarg;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'o':
      /* name of xml injection file */
      CLA->outputFileName=optarg;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'S':
      /* GPS start time of search */
       CLA->GPSStart=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "int");
      break;
    case 'E':
       /* GPS end time time of search */
      CLA->GPSEnd=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "int");
      break;
    case 'd':
       /* Number of segment to break-up search into */
      CLA->ShortSegDuration=atoi(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "int");
      break;
    case 'T':
      /* Half the number of seconds that are trown out at the start and at the end of a short chunk */
      CLA->TruncSecs=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "int");
      break;
    case 'g':
      /* start time of allowed triggers */
      CLA->trigstarttime=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "int");
      break;
    case 'p':
      /* start time of allowed triggers */
      CLA->pad=atoi(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "int");
      break;
    case 'c':
      /* cusp power law */
      CLA->power=-4.0/3.0;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'k':
      /* kink power law */
      CLA->power=-5.0/3.0;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'n':
      /* fake gaussian noise flag */
      CLA->fakenoiseflag=1;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'w':
      /* fake gaussian noise flag */
      CLA->whitespectrumflag=1;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'l':
      /* fake gaussian noise flag */
      CLA->cluster=atof(optarg);
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'a':
      /* fake gaussian noise flag */
      CLA->printspectrumflag=1;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'b':
      /* fake gaussian noise flag */
      CLA->printfilterflag=1;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'r':
      /* fake gaussian noise flag */
      CLA->printsnrflag=1;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'x':
      /* fake gaussian noise flag */
      CLA->printfirflag=1;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'y':
      /* fake gaussian noise flag */
      CLA->printdataflag=1;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'z':
      /* fake gaussian noise flag */
      CLA->printinjectionflag=1;
      ADD_PROCESS_PARAM(procTable.processTable, "string");
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"All arguments are required except -n, -h, -w, -g, -o, -x, -y, -z, -b, -r -a, -l, -p, -T  and -i. One of -k or -c must be specified. They are:\n");
      fprintf(stdout,"\t--bw-flow (-f)\t\tFLOAT\t Low frequency cut-off.\n");
      fprintf(stdout,"\t--sample-rate (-s)\t\tFLOAT\t Desired sample rate (Hz).\n");
      fprintf(stdout,"\t--bank-lowest-hifreq-cutoff (-H)\tFLOAT\t Template bank lowest high frequency cut-off.\n");
      fprintf(stdout,"\t--max-mismatch (-M)\tFLOAT\t Maximal mismatch allowed from 1 template to the next.\n");
      fprintf(stdout,"\t--bank-freq-start (-L)\tFLOAT\t Template bank low frequency cut-off.\n");
      fprintf(stdout,"\t--threshold (-t)\t\tFLOAT\t SNR threshold.\n");
      fprintf(stdout,"\t--frame-cache (-F)\t\tSTRING\t Name of frame cache file.\n");
      fprintf(stdout,"\t--channel-name (-C)\t\tSTRING\t Name of channel.\n");
      fprintf(stdout,"\t--injection-file (-i)\t\tSTRING\t Name of xml injection file.\n");
      fprintf(stdout,"\t--outfile (-o)\t\tSTRING\t Name of xml output file.\n");
      fprintf(stdout,"\t--gps-start-time (-S)\t\tINTEGER\t GPS start time.\n");
      fprintf(stdout,"\t--gps-end-time (-E)\t\tINTEGER\t GPS end time.\n");
      fprintf(stdout,"\t--settling-time (-T)\t\tINTEGER\t Number of seconds to truncate filter.\n");
      fprintf(stdout,"\t--trig-start-time (-g)\t\tINTEGER\t GPS start time of triggers to consider.\n");
      fprintf(stdout,"\t--pad (-p)\t\tINTEGER\t Pad the data with these many seconds at beginning and end.\n");
      fprintf(stdout,"\t--short-segment-duration (-d)\t\tINTEGER\t Duration of shor segments. They will overlap by 50%s. \n","%");
      fprintf(stdout,"\t--kink-search (-k)\t\tFLAG\t Specifies a search for string kinks.\n");
      fprintf(stdout,"\t--cusp-search (-c)\t\tFLAG\t Specifies a search for string cusps.\n");
      fprintf(stdout,"\t--test-gaussian-data (-n)\tFLAG\t Use unit variance fake gaussian noise.\n");
      fprintf(stdout,"\t--test-white-spectrum (-w)\tFLAG\t Use constant white noise (used only in combination with fake gaussian noise; otherwise ignored).\n");
      fprintf(stdout,"\t--cluster-events (-l)\tREAL4\t Cluster events with input timescale.\n");
      fprintf(stdout,"\t--print-spectrum (-a)\tFLAG\t Prints the spectrum to Spectrum.txt.\n");
      fprintf(stdout,"\t--print-fd-filter (-b)\tFLAG\t Prints the frequency domain filter to Filter-<template no>.txt.\n");      
      fprintf(stdout,"\t--print-td-filter (-r)\tFLAG\t Prints the time domain filter to FIRFilter-<template no>.txt.\n");      
      fprintf(stdout,"\t--print-snr (-x)\tFLAG\t Prints the snr to stdout.\n");      
      fprintf(stdout,"\t--print-data (-y)\tFLAG\t Prints the post-processed (HP filtered, downsampled, padding removed, with injections) data to data.txt.\n");
      fprintf(stdout,"\t--print-injection (-z)\tFLAG\t Prints the injeciton data to injection.txt.\n");      
      fprintf(stdout,"\t--help (-h)\t\t\tFLAG\t Print this message.\n");
      fprintf(stdout,"eg %s  --sample-rate 4096 --bw-flow 39 --bank-freq-start 30 --bank-lowest-hifreq-cutoff 200 --settling-time 0.1 --short-segment-duration 4 --cusp-search --cluster-events 0.1 --pad 4 --threshold 4 --outfile ladida.xml --frame-cache cache/H-H1_RDS_C01_LX-795169179-795171015.cache --channel-name H1:LSC-STRAIN --gps-start-time 795170318 --gps-end-time 795170396\n", argv[0]);
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
  }

  if(CLA->flow == 0.0)
    {
      fprintf(stderr,"No low cutoff frequency specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->fbankstart == 0.0)
    {
      fprintf(stderr,"No low frequency for frequency bank specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->fbankhighfcutofflow == 0.0)
    {
      fprintf(stderr,"No template bank lowest high frequency cutoff specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->fmismatchmax == 0.0){
    fprintf(stderr,"No maximal mismatch specified.\n");
    fprintf(stderr,"Try %s -h \n",argv[0]);
    return 1;
  }      
  if(CLA->threshold == 0.0)
    {
      fprintf(stderr,"No SNR threshold specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->power == 0.0)
    {
      fprintf(stderr,"Cusp or kink search not specified. \n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->FrCacheFile == NULL)
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->ChannelName == NULL)
    {
      fprintf(stderr,"No channel name specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->GPSStart == 0)
    {
      fprintf(stderr,"No GPS start time specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->GPSEnd == 0)
    {
      fprintf(stderr,"No GPS end time specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->ShortSegDuration == 0)
    {
      fprintf(stderr,"Short segment duration not specified (they overlap by 50%s).\n","%");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      

  /* Some consistency checking */
  {
    int big_seg_length=CLA->GPSEnd-CLA->GPSStart-2*CLA->pad;
    int small_seg_length=CLA->ShortSegDuration;

    REAL4 x=((float)big_seg_length/(float)small_seg_length)-0.5;

    if((int)x != x){
      fprintf(stderr,"The total duration of the segment T and the short segment duration\n");
      fprintf(stderr,"Should obey the following rule: T/t - 0.5 shold be an odd integer.\n");
      return 1;
    } 
    if(((int)x)%2 != 1){
      fprintf(stderr,"The total duration of the segment T and the short segment duration\n");
      fprintf(stderr,"Should obey the following rule: T/t - 0.5 shold be an odd integer.\n");
      return 1;
    }     

    if( CLA->ShortSegDuration/4.0  < CLA->TruncSecs){
      fprintf(stderr,"Short segment length t=%d is too small to accomodate truncation time requested.\n", small_seg_length);
	fprintf(stderr,"Need short segment t(=%d) to be >= 4 x Truncation length (%f).\n",CLA->ShortSegDuration,CLA->TruncSecs);
	return 1;
    }    
  }

  /* check mismatch */
  {
    if(CLA->fmismatchmax < 0.0 || CLA->fmismatchmax > 1.0){
      fprintf(stderr,"ERROR : the maximal mismatch is not authorized.\n");
      return 1;
    }      
  }
  /* check frequencies */
  {
    REAL4 f99=CLA->flow*pow((1/0.9-1)/(1/0.99-1),0.25);
    if(CLA->fbankstart < f99)
      fprintf(stderr,"WARNING: Template starting frequency and BW high pass frequency are close. f99=%e, fbw=%e\n",f99, CLA->flow);
  }

  /* store the input start and end times */
  /* set the start and end time for the search summary */
  {
    int small_seg_length=CLA->ShortSegDuration;

    searchsumm.searchSummaryTable->in_start_time.gpsSeconds = CLA->GPSStart;
    searchsumm.searchSummaryTable->in_start_time.gpsNanoSeconds =0;
    searchsumm.searchSummaryTable->in_end_time.gpsSeconds = CLA->GPSEnd;
    searchsumm.searchSummaryTable->in_end_time.gpsNanoSeconds =0;

    if (CLA->trigstarttime > 0)
      searchsumm.searchSummaryTable->out_start_time.gpsSeconds = CLA->trigstarttime;
    else
      searchsumm.searchSummaryTable->out_start_time.gpsSeconds = CLA->GPSStart+small_seg_length/4+CLA->pad;
      
    searchsumm.searchSummaryTable->out_start_time.gpsNanoSeconds =0;
    searchsumm.searchSummaryTable->out_end_time.gpsSeconds = CLA->GPSEnd-small_seg_length/4-CLA->pad;
    searchsumm.searchSummaryTable->out_end_time.gpsNanoSeconds =0;
  }

  return errflg;
}

/*******************************************************************************/

int FreeMem(void){
  int m;
  
  XLALDestroyREAL4TimeSeries(GV.ht_proc);
  XLALDestroyREAL4Vector(GV.Spec.data);
  
  for (m=0; m < NTemplates; m++)
    XLALDestroyREAL4Vector(strtemplate[m].StringFilter.data);
  
  XLALDestroyREAL4FFTPlan( GV.fplan );
  XLALDestroyREAL4FFTPlan( GV.rplan );
  
  LALCheckMemoryLeaks();
  
  return 0;
}

/*******************************************************************************/
#endif
