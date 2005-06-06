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

#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>

#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>

#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/GenerateBurst.h>


#include <lal/BurstUtils.h>

#include <lalapps.h>
#include <processtable.h>

extern char *optarg;
extern int optind, opterr, optopt;

int snprintf(char *str, size_t size, const char *format, ...);
int vsnprintf(char *str, size_t size, const char *format, va_list ap);

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

#define ADD_PROCESS_PARAM(type) \
	do { paramaddpoint = add_process_param(paramaddpoint, type, long_options[option_index].name, optarg); } while(0)

#define SCALE 1e20
#define MAXTEMPLATES 10000

NRCSID( STRINGSEARCHC, "StringSearch $Id$");
RCSID( "StringSearch $Id$");

#define PROGRAM_NAME "StringSearch"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"


/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL4 flow;                 /* Low frequency cut-off */
  REAL8 samplerate;           /* desired sample rate */
  REAL4 fbanklow;             /* Template bank low frequency cut-off */
  char *FrCacheFile;          /* Frame cache file */
  char *InjectionFile;        /* LIGO xml injection file */
  char *ChannelName;          /* Name of channel to be read in from frames */
  char *outputFileName;         /* Name of xml output filename */
  INT4 GPSStart;              /* GPS start time of segment to be analysed */
  INT4 GPSEnd;                /* GPS end time of segment to be analysed */
  INT4 ShortSegDuration;      /* Number of fixed length sub-segments between GPSStart and GPSEnd */
  REAL4 TruncSecs;             /* Half the number of seconds truncated at beginning and end of a chunk */
  REAL4 power;                /* Kink (-5/3) or cusp (-4/3) frequency power law */
  REAL4 threshold;            /* event SNR threshold */
  INT4 fakenoiseflag;         /* =0 if real noise =1 if fake gaussian noise */
  INT4 whitespectrumflag;     /* =0 if spectrum is to be computed =1 for white spectrum */
  INT4 trigstarttime;         /* start-time of allowed triggers */
  INT4 cluster;               /* =0 if events are not to be clustered =1 otherwise */
  INT4 pad;                   /* seconds of padding */
} CommandLineArgs;

typedef 
struct GlobalVariablesTag {
  INT4 duration;                /* duration of entire segment to be analysed */
  LIGOTimeGPS gpsepoch;         /* GPS epoch of start of entire segment to be analysed */ 
  REAL8TimeSeries ht;           /* raw input data */
  REAL4TimeSeries ht_proc;      /* processed (band-pass filtered and down-sampled) input data */
  REAL4FrequencySeries Spec;    /* average spectrum */
  REAL4FrequencySeries StringFilter;    /* inverse truncated average spectrum x f^(-4/3) or f^(-5/3) */
  RealFFTPlan *fplan;           /* fft plans */
  RealFFTPlan *rplan;           /* fft plans */
  INT4 seg_length;
} GlobalVariables;

typedef
struct StringTemplateTag {
  REAL4 f;                      /* Template frequency */
  REAL4 norm;                   /* Template normalisation */
  REAL4 mismatch;               /* Template mismatch relative to last one*/
} StringTemplate;

/***************************************************************************/

/* GLOBAL VARIABLES */

static LALStatus status;
INT4 lalDebugLevel=3;
FrCache *framecache;                                           /* frame reading variables */
FrStream *framestream=NULL;

REAL4 Templateflow;

GlobalVariables GV;   /* A bunch of stuff is stored in here; mainly to protect it from accidents */

StringTemplate strtemplate[MAXTEMPLATES];

int NTemplates;

LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;

SnglBurstTable *events=NULL;
MetadataTable  procTable;
MetadataTable  procparams;
MetadataTable  searchsumm;

CHAR outfilename[256];
CHAR ifo[3]; 

REAL4 SAMPLERATE;

int Nevents=0;

PassBandParamStruc highpassParams;

REAL8Vector *window;

/***************************************************************************/

/* FUNCTION PROTOTYPES */

/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Reads raw data (or puts in fake gaussian noise with a sigma=10^-20) */
int ReadData(struct CommandLineArgsTag CLA);

/* windows the data with a Tukey window */
int WindowData(void);

/* High pass filters and casts data to REAL4 */
int ProcessData(void);

/* Adds injections if an xml injection file is given */
int AddInjections(struct CommandLineArgsTag CLA);

/* High pass filters data again if an injection has been made */
int ProcessData2(void);

/* DownSamples data */
int DownSample(struct CommandLineArgsTag CLA);

/* Computes the average spectrum  */
int AvgSpectrum(struct CommandLineArgsTag CLA);

/* Creates the frequency domain string cusp or kink filter  */
int CreateStringFilter(struct CommandLineArgsTag CLA);

/* Creates the template bank based on the spectrum  */
int CreateTemplateBank(struct CommandLineArgsTag CLA);

/* Filters the data through the template banks  */
int FindStringBurst(struct CommandLineArgsTag CLA);

/* Finds events above SNR threshold specified  */
int FindEvents(struct CommandLineArgsTag CLA, REAL4Vector *vector, 
	       INT4 i, INT4 m, SnglBurstTable **thisEvent);

/* Writes out the xml file with the events it found  */
int OutputEvents(struct CommandLineArgsTag CLA);

/* Frees the memory */
int FreeMem(void);                                        

/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{

 if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
 
 /* set lowest frequency of template bank; 
    should be higher than BW high pass frequency */ 
 Templateflow=CommandLineArgs.flow+10.0;

 highpassParams.nMax =  4;
 highpassParams.f1   = -1;
 highpassParams.a1   = -1;
 highpassParams.f2   = CommandLineArgs.flow;
 highpassParams.a2   = 0.5; /* this means 50% attenuation at f2 */

 if (ReadData(CommandLineArgs)) return 2;

 if (WindowData()) return 3;

 if (ProcessData()) return 4;

 /* after this there's no more GV.ht, only GV.ht_proc */

 if (CommandLineArgs.InjectionFile != NULL) 
   {
     if (AddInjections(CommandLineArgs)) return 3;
     /* high pass filter data again with added injection */ 
     if (ProcessData2()) return 3;
   }
  	 
 if (DownSample(CommandLineArgs)) return 5;
	
 /* Here I need to re-size the time series to remove the pad */

 LALResizeREAL4TimeSeries(&status, &(GV.ht_proc), (int)(CommandLineArgs.pad/GV.ht_proc.deltaT+0.5),
			  GV.ht_proc.data->length-2*(UINT4)(CommandLineArgs.pad/GV.ht_proc.deltaT+0.5));
 TESTSTATUS( &status );
 /* reduce duration of segment appropriately */
 GV.duration -= 2*CommandLineArgs.pad; 

 if (AvgSpectrum(CommandLineArgs)) return 6;
 
 if (CreateStringFilter(CommandLineArgs)) return 7;
 
 if (CreateTemplateBank(CommandLineArgs)) return 8;

 if (FindStringBurst(CommandLineArgs)) return 9;

 if (CommandLineArgs.cluster == 1) 
   {
     XLALClusterStringBurstTable(&events, NULL, XLALCompareStringBurstByTime);
   }

 if (OutputEvents(CommandLineArgs)) return 11;

 if (FreeMem()) return 12;

 return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/


/*******************************************************************************/


int WindowData(void)
{
  REAL8 r = 0.001;
  int k;
  int N=GV.ht.data->length;
  int kl=r/2*(N-1)+1;
  int kh=N-r/2*(N-1)+1;

  LALDCreateVector(&status,&window,GV.ht.data->length);
  TESTSTATUS( &status );

  for(k = 1; k < kl; k++) 
    {
      window->data[k-1]=0.5*( 1 + cos( LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI));
    }

  for(k = kl; k < kh; k++) 
    {
      window->data[k-1]=1.0;
    }

  for(k = kh; k <= N; k++) 
    {
      window->data[k-1]=0.5*( 1 + cos( LAL_TWOPI/r - LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI));
    }

  for (k = 0; k < N; k++)
    {
      GV.ht.data->data[k] *= window->data[k];
    }
  
  LALDDestroyVector(&status,&window);
  TESTSTATUS( &status );

  return 0;
}

/*******************************************************************************/

int ProcessData2(void)
{

  LALDButterworthREAL4TimeSeries( &status, &GV.ht_proc, &highpassParams ); 
  TESTSTATUS( &status ); 

  return 0;
}

/*******************************************************************************/

int AddInjections(struct CommandLineArgsTag CLA)
{

  INT4 startTime = GV.ht.epoch.gpsSeconds;
  INT4 stopTime = startTime + GV.ht_proc.data->length * GV.ht_proc.deltaT;
  SimBurstTable *injections = NULL;
  int i;
  INT4 calType=0;

  COMPLEX8 one = {1.0, 0.0};
  COMPLEX8FrequencySeries *response = NULL;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

  /* Get info from injection file */
  LALSimBurstTableFromLIGOLw(&status, &injections, CLA.InjectionFile, startTime, stopTime);

  /* make the response function */
  LALCreateCOMPLEX8FrequencySeries(&status, &response, CLA.ChannelName, GV.ht_proc.epoch, 0.0, 
				   1.0 / (GV.ht_proc.data->length * GV.ht_proc.deltaT), 
				   strainPerCount, GV.ht_proc.data->length / 2 + 1);
  
  for(i = 0; i < (int)response->data->length; i++)
			response->data->data[i] = one;

  /* Inject the signals into the data */
  LALBurstInjectSignals(&status, &GV.ht_proc, injections, response, calType); 

  /* free the injection table */
  while(injections) {
    SimBurstTable *thisEvent = injections;
    injections = injections->next;
    LALFree(thisEvent);
  }

  LALDestroyCOMPLEX8FrequencySeries(&status, response);

  return 0;
}

/*******************************************************************************/

static ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param, 
					      const char *type, const char *param, const char *value)
{
  *proc_param = LALCalloc(1, sizeof(**proc_param));
  (*proc_param)->next = NULL;
  snprintf((*proc_param)->program, LIGOMETA_PROGRAM_MAX, PROGRAM_NAME);
  snprintf((*proc_param)->type, LIGOMETA_TYPE_MAX, type);
  snprintf((*proc_param)->param, LIGOMETA_PARAM_MAX, "--%s", param);
  snprintf((*proc_param)->value, LIGOMETA_VALUE_MAX, value);
  
  return(&(*proc_param)->next);
}

/*******************************************************************************/

int OutputEvents(struct CommandLineArgsTag CLA)
{  
  LIGOLwXMLStream xml;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;
  MetadataTable myTable;

  if (!CLA.outputFileName)
    {
      snprintf(outfilename, sizeof(outfilename)-1, "%s-STRINGSEARCH-%d-%d.xml", ifo,
	       searchsumm.searchSummaryTable->in_start_time.gpsSeconds,
	       searchsumm.searchSummaryTable->in_end_time.gpsSeconds - 
	       searchsumm.searchSummaryTable->in_start_time.gpsSeconds);
      outfilename[sizeof(outfilename)-1] = '\0';
    }else{
      snprintf(outfilename, sizeof(outfilename)-1, "%s", CLA.outputFileName);
      outfilename[sizeof(outfilename)-1] = '\0';
    }

  memset(&xml, 0, sizeof(LIGOLwXMLStream));
  LALOpenLIGOLwXMLFile(&status, &xml, outfilename);

  /* process table */
  snprintf(procTable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo);
  LALGPSTimeNow(&status, &(procTable.processTable->start_time), &accuracy);
  LALBeginLIGOLwXMLTable(&status, &xml, process_table);
  LALWriteLIGOLwXMLTable(&status, &xml, procTable, process_table);
  LALEndLIGOLwXMLTable(&status, &xml);

  /* process params table */
  LALBeginLIGOLwXMLTable(&status, &xml, process_params_table);
  LALWriteLIGOLwXMLTable(&status, &xml, procparams, process_params_table);
  LALEndLIGOLwXMLTable(&status, &xml);

  /* search summary table */
  snprintf(searchsumm.searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo);
  searchsumm.searchSummaryTable->nevents = Nevents;

  LALBeginLIGOLwXMLTable(&status, &xml, search_summary_table);
  LALWriteLIGOLwXMLTable(&status, &xml, searchsumm, search_summary_table);
  LALEndLIGOLwXMLTable(&status, &xml);

  /* burst table */
  LALBeginLIGOLwXMLTable(&status, &xml, sngl_burst_table);

  myTable.snglBurstTable = events;      
    
  LALWriteLIGOLwXMLTable(&status, &xml, myTable, sngl_burst_table);
  LALEndLIGOLwXMLTable(&status, &xml);

  LALCloseLIGOLwXMLFile(&status, &xml);
  
  /* free event list, process table, search summary and process params */

  while ( events )
  {
    SnglBurstTable *next = events->next;
    LALFree( events );
    events = next;
  }


  LALFree(procTable.processTable);
  LALFree(searchsumm.searchSummaryTable);

  while(procparams.processParamsTable) {
    ProcessParamsTable *table = procparams.processParamsTable;
    procparams.processParamsTable = table->next;
    LALFree(table);
  }

  return 0;
}

/*******************************************************************************/

int FindEvents(struct CommandLineArgsTag CLA, REAL4Vector *vector, INT4 i, INT4 m, SnglBurstTable **thisEvent)
{
  int p,pend=0,pstart=0;

  /* Now find thisEvent in the inner half */
  for ( p = (int)vector->length/4 ; p < (int)(3*vector->length/4); p++ )
    {
      REAL4 maximum = 0.0;
      INT4 pmax=p;
      INT8 timeNS  = (INT8)( 1000000000 ) * (INT8)(GV.ht_proc.epoch.gpsSeconds+GV.seg_length*i/2*GV.ht_proc.deltaT);
      timeNS += (INT8)( 1e9 * GV.ht_proc.deltaT * p );

      if ( (fabs(vector->data[p]) > CLA.threshold) && ( (double)(1e-9*timeNS) > (double)CLA.trigstarttime))
	{
          INT8  timeNS, peaktime;
	  REAL8 duration;
	  pstart=p;
	  timeNS  = (INT8)( 1000000000 ) * (INT8)(GV.ht_proc.epoch.gpsSeconds+GV.seg_length*i/2*GV.ht_proc.deltaT);
          timeNS += (INT8)( 1e9 * GV.ht_proc.deltaT * p );

	  if ( *thisEvent ) /* create a new event */
            {
              (*thisEvent)->next = LALCalloc( 1, sizeof( *(*thisEvent)->next ) );
              *thisEvent = (*thisEvent)->next;
            }
	  else /* create the list */
            {
              *thisEvent = events = LALCalloc( 1, sizeof( *events ) );
            }

	  if ( ! *thisEvent ) /* allocation error */
	    {
	      fprintf(stderr,"Could not allocate memory for event. Memory allocation error. Exiting. \n");
	      return 1;
	    }
	  while( ((fabs(vector->data[p]) > CLA.threshold) || ((p-pstart)* GV.ht_proc.deltaT < (float)CLA.TruncSecs)) 
		 && p<(int)(3*vector->length/4))
	    {
	      if(fabs(vector->data[p]) > maximum) 
		{
		  maximum=fabs(vector->data[p]);
		  pmax=p;
		}
	      p++;
	    }
	  pend=p-1;

	  peaktime = timeNS + (INT8)( 1e9 * GV.ht_proc.deltaT * (pmax-pstart) );
	  duration = GV.ht_proc.deltaT * ( (p-1) - pstart );

	  /* Now copy stuff into event */
	  strncpy( (*thisEvent)->ifo, CLA.ChannelName, sizeof(ifo)-1 );
	  strncpy( (*thisEvent)->search, "StringCusp", sizeof( (*thisEvent)->search ) );
	  strncpy( (*thisEvent)->channel, CLA.ChannelName, sizeof( (*thisEvent)->channel ) );

	  /* give trigger a 1 sample fuzz on either side */
	  timeNS -= GV.ht_proc.deltaT *1e9;
	  duration += 2*GV.ht_proc.deltaT;

	  (*thisEvent)->start_time.gpsSeconds     = timeNS / 1000000000;
	  (*thisEvent)->start_time.gpsNanoSeconds = timeNS % 1000000000;
	  (*thisEvent)->peak_time.gpsSeconds      = peaktime / 1000000000;
	  (*thisEvent)->peak_time.gpsNanoSeconds  = peaktime % 1000000000;
	  (*thisEvent)->duration     = duration;
	  (*thisEvent)->central_freq = (strtemplate[m].f+Templateflow)/2.0;	   
	  (*thisEvent)->bandwidth    = strtemplate[m].f-Templateflow;				     
	  (*thisEvent)->snr          = maximum;
	  (*thisEvent)->amplitude   = vector->data[pmax]/strtemplate[m].norm;
	  (*thisEvent)->confidence   = -maximum; /* FIXME */

	}
    }

  return 0;
}

/*******************************************************************************/

int FindStringBurst(struct CommandLineArgsTag CLA)
{
  int i,p,m; 
  REAL4Vector *vector = NULL;
  COMPLEX8Vector *vtilde = NULL;
  SnglBurstTable *thisEvent = NULL;

  int f_low_index = Templateflow / GV.StringFilter.deltaF;
  
  /* create vector that will hold the data for each overlapping chunk */ 
  LALSCreateVector( &status, &vector, GV.seg_length);
  TESTSTATUS( &status );
  /* create vector that will hold FFT of data*/
  LALCCreateVector( &status, &vtilde, GV.seg_length / 2 + 1 );
  TESTSTATUS( &status );

  /* loop over overlapping chunks */ 
  for(i=0; i < 2*GV.duration/CLA.ShortSegDuration - 1 ;i++)
    {
      /* populate vector that will hold the data for each overlapping chunk */
      memcpy( vector->data, GV.ht_proc.data->data + i*GV.seg_length/2,vector->length*sizeof( *vector->data ) );

      /* fft it */
      LALForwardRealFFT( &status, vtilde, vector, GV.fplan );
      TESTSTATUS( &status );

      /* this sets to zero all data below the low frequency cutoff */ 
      memset( vtilde->data, 0, f_low_index  * sizeof( *vtilde->data ) );
   
      /* multiply FT of data and String Filter and deltaT (latter not included in LALForwardRealFFT) */
      for ( p = (int)f_low_index ; p < (int)GV.StringFilter.data->length; p++ )
	{
	  vtilde->data[p].re *= GV.StringFilter.data->data[p]*GV.ht_proc.deltaT;
	  vtilde->data[p].im *= GV.StringFilter.data->data[p]*GV.ht_proc.deltaT;
	}

      /* loop over templates  */
      for (m = 0; m < NTemplates; m++)
	{
	  int f_high_index = strtemplate[m].f/ GV.StringFilter.deltaF;
	  
	  /* set to zero all values greater than the high frequency cutoff of the template */
	  memset( vtilde->data+f_high_index, 0, (vtilde->length-f_high_index) * sizeof( *vtilde->data ) );

	  LALReverseRealFFT( &status, vector, vtilde,  GV.rplan);
	  TESTSTATUS( &status );

	  /* normalise the result by template normalisation and multiply by 
	     df (not inluded in LALReverseRealFFT)  factor of 2 is from 
	     match-filter definition */

	  for ( p = 0 ; p < (int)vector->length; p++ )
	    {
	      vector->data[p] *= 2.0 * GV.StringFilter.deltaF / strtemplate[m].norm;
	    }

	  if(FindEvents(CLA, vector, i, m, &thisEvent)) return 1;
	}
    }

 /* sort events in time; if there are any */
 if (events) 
   {
     /* first sort list in increasing GPS peak time */
     LALSortSnglBurst(&status, &events, XLALCompareSnglBurstByPeakTimeAndSNR);
     TESTSTATUS( &status );
   }
  
  LALSDestroyVector( &status, &vector );
  TESTSTATUS( &status );
  LALCDestroyVector( &status, &vtilde );
  TESTSTATUS( &status );

  return 0;
}


/*******************************************************************************/

int CreateTemplateBank(struct CommandLineArgsTag CLA)
{
  REAL8 fmax,f,t1t1,t2t2, epsilon;
  int p,f_low_index,f_high_index,k;

  fmax = (1.0/GV.ht_proc.deltaT) / 2.0;

  f_low_index = Templateflow / GV.StringFilter.deltaF;
  f_high_index = fmax / GV.StringFilter.deltaF;
  
  t1t1=0.0;
  for ( p = f_low_index ; p < f_high_index; p++ )
    {
      f= p*GV.StringFilter.deltaF;
      t1t1 += 4*pow(f,CLA.power)*GV.StringFilter.data->data[p]*GV.StringFilter.deltaF;
    }

  /* This is the first template, all others will be slightly smaller */ 
  strtemplate[0].f=fmax;
  strtemplate[0].norm=sqrt(t1t1);
  strtemplate[0].mismatch=0.0;
  
  t2t2=t1t1;
  k=1;

  f_low_index = CLA.fbanklow / GV.StringFilter.deltaF;
  /* now we loop through and take away from the integral one point at a time */
  for ( p = f_high_index-2 ; p >= f_low_index; p-- )
    {
      f= p*GV.StringFilter.deltaF;
      t2t2 -= 4*pow(f,CLA.power)*GV.StringFilter.data->data[p]*GV.StringFilter.deltaF;
      
      epsilon=1-sqrt(t2t2/t1t1);
      
      /* FIXME: epsilon should be a command line argument */
      if(epsilon >= 0.02)
	{
	  t1t1=t2t2;
	  strtemplate[k].f=f;
	  strtemplate[k].norm=sqrt(t1t1);
	  strtemplate[k].mismatch=epsilon;
	  k++;
	}
      if(k == MAXTEMPLATES)
	{
	  fprintf(stderr,"Too many templates for code... Exiting\n");
	  return 1;
	}
    }

  NTemplates=k-1;

  return 0;
}


/*******************************************************************************/
int CreateStringFilter(struct CommandLineArgsTag CLA)
{
  COMPLEX8Vector *vtilde = NULL;
  REAL4Vector    *vector = NULL;

  int p,f_cutoff_index; 
 
  /* Create vector that will hold the string filter */
  GV.StringFilter.deltaF=GV.Spec.deltaF;
  LALSCreateVector( &status, &GV.StringFilter.data, GV.Spec.data->length );
  TESTSTATUS( &status );

  /* We need to truncate the filter (in the time domain) to ensure it has a 
     finite (and small) response time */

  LALCCreateVector( &status, &vtilde, GV.Spec.data->length );
  TESTSTATUS( &status );

  f_cutoff_index = CLA.flow/ GV.StringFilter.deltaF;
  
  memset( vtilde->data, 0, f_cutoff_index  * sizeof( *vtilde->data ) );
 
  for ( p = (int)f_cutoff_index ; p < (int)vtilde->length - 1; p++ )
    {
      REAL4 f=p*GV.StringFilter.deltaF;
      
      vtilde->data[p].re = sqrt(pow(f,CLA.power)/(GV.Spec.data->data[p]));
      vtilde->data[p].im = 0;
    }
  vtilde->data[vtilde->length - 1].re = 0;
  vtilde->data[vtilde->length - 1].im = 0;

  LALSCreateVector( &status, &vector, GV.seg_length );
  TESTSTATUS( &status );

  LALReverseRealFFT( &status, vector, vtilde,  GV.rplan);
  TESTSTATUS( &status );
  
  /* multiply times df to make sure units are correct */
  for ( p = 0 ; p < (int)vector->length; p++ )
      vector->data[p] *= GV.StringFilter.deltaF;
  
  if(CLA.TruncSecs != 0.0) 
    {
      memset( vector->data + (INT4)(CLA.TruncSecs/GV.ht_proc.deltaT +0.5), 0,
	      ( vector->length -  2 * (INT4)(CLA.TruncSecs/GV.ht_proc.deltaT +0.5)) 
	      * sizeof( *vector->data ) );
    }
  
  LALForwardRealFFT( &status, vtilde, vector,  GV.fplan);
  TESTSTATUS( &status );

  for ( p = 0 ; p < (int)vtilde->length - 1; p++ )
    {
      REAL4 re = vtilde->data[p].re * GV.ht_proc.deltaT;
      REAL4 im = vtilde->data[p].im * GV.ht_proc.deltaT;

      GV.StringFilter.data->data[p] = (re * re + im * im);
    }

  /* set all values below the cutoff frequency to 0 */
  memset( GV.StringFilter.data->data, 0, f_cutoff_index  * 
	  sizeof( *GV.StringFilter.data->data ) );

  LALCDestroyVector( &status, &vtilde );
  TESTSTATUS( &status );

  LALSDestroyVector( &status, &vector );
  TESTSTATUS( &status );

  return 0;
}

/*******************************************************************************/
int AvgSpectrum(struct CommandLineArgsTag CLA)
{
 
  GV.seg_length = (int)(CLA.ShortSegDuration/GV.ht_proc.deltaT + 0.5);

  LALSCreateVector( &status, &GV.Spec.data, GV.seg_length / 2 + 1 );
  TESTSTATUS( &status );

  LALCreateForwardRealFFTPlan( &status, &GV.fplan, GV.seg_length, 0 );
  TESTSTATUS( &status );
  
  LALCreateReverseRealFFTPlan( &status, &GV.rplan, GV.seg_length, 0 );
  TESTSTATUS( &status );

  if (CLA.fakenoiseflag && CLA.whitespectrumflag)
    {
      int p;
      for ( p = 0 ; p < (int)GV.Spec.data->length; p++ )
	{
	  GV.Spec.data->data[p]=2/SAMPLERATE;
	}
      GV.Spec.deltaF=1/(GV.seg_length*GV.ht_proc.deltaT);
    }
  else
    {
      int segmentLength = GV.seg_length;
      int segmentStride = GV.seg_length/2;
      REAL4Window  *window  = NULL;

      window = XLALCreateWelchREAL4Window( segmentLength );
      XLALREAL4AverageSpectrumMedianMean( &GV.Spec, &GV.ht_proc, segmentLength,
					  segmentStride, window, GV.fplan );
      XLALDestroyREAL4Window( window );
    }

  return 0;
}

/*******************************************************************************/

int DownSample(struct CommandLineArgsTag CLA)
{
  ResampleTSParams resamplepar;

  memset( &resamplepar, 0, sizeof( resamplepar ) );
  resamplepar.deltaT     = 1.0/CLA.samplerate;
  resamplepar.filterType = defaultButterworth;

  LALResampleREAL4TimeSeries( &status, &GV.ht_proc, &resamplepar );
  TESTSTATUS( &status );

  return 0;
}

/*******************************************************************************/

int ProcessData(void)
{
  int p;

  LALButterworthREAL8TimeSeries( &status, &GV.ht, &highpassParams ); 
  TESTSTATUS( &status ); 
  
  TESTSTATUS( &status ); 

  for (p=0; p<(int)GV.ht.data->length; p++)  
    {
      GV.ht_proc.data->data[p]=GV.ht.data->data[p]; 
    } 

  strncpy(GV.ht_proc.name, GV.ht.name, LALNameLength);
  GV.ht_proc.deltaT=GV.ht.deltaT;
  GV.ht_proc.epoch=GV.ht.epoch;
  GV.ht_proc.sampleUnits=GV.ht.sampleUnits;

  /* destroy double precision vector */
  LALDDestroyVector(&status,&GV.ht.data); 
  TESTSTATUS( &status ); 

  return 0;
}

/*******************************************************************************/

int ReadData(struct CommandLineArgsTag CLA)
{
  static FrChanIn chanin_ht;
  int p;

  /* create Frame cache, open frame stream and delete frame cache */
  LALFrCacheImport(&status,&framecache,CLA.FrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&framestream,framecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&framecache);
  TESTSTATUS( &status );
  
  /* Define channel */
  chanin_ht.type  = ProcDataChannel;
  chanin_ht.name  = CLA.ChannelName;
  
  GV.duration=CLA.GPSEnd-CLA.GPSStart;
  GV.gpsepoch.gpsSeconds=CLA.GPSStart; /* Set global variable epoch */
  GV.gpsepoch.gpsNanoSeconds=0;

  LALFrSeek(&status,&GV.gpsepoch,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL8TimeSeries(&status,&GV.ht,&chanin_ht,framestream);
  TESTSTATUS( &status );

  /* Allocate space for data vectors */
  LALDCreateVector(&status,&GV.ht.data,(UINT4)(GV.duration/GV.ht.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&GV.ht_proc.data,(UINT4)(GV.duration/GV.ht.deltaT +0.5));
  TESTSTATUS( &status );

  SAMPLERATE=1.0/GV.ht.deltaT;

  /* If we are reading real noise then read it*/
  if (!CLA.fakenoiseflag)
    {
      /* Read data */
      LALFrSeek(&status,&GV.gpsepoch,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL8TimeSeries(&status,&GV.ht,&chanin_ht,framestream);
      TESTSTATUS( &status );

      /* Scale data to avoid single float precision problems */
      for (p=0; p<(int)GV.ht_proc.data->length; p++)
	{
	  GV.ht.data->data[p] *= SCALE;
	}
    }
  /* otherwise create random data set */
  else
    {
      FILE *devrandom;
      RandomParams   *randpar=NULL;
      REAL4Vector    *v1=NULL;
      int seed, errorcode;

      if (!(devrandom=fopen("/dev/urandom","r")))
	{
	  fprintf(stderr,"Unable to open device /dev/urandom\n");
	  return 1;
	}
      errorcode=fread((void*)&seed,sizeof(INT4),1,devrandom);
      if (errorcode!=1)
	{
	  fprintf( stderr,"Error reading /dev/urandom file!\n");
	  return 1;
	}
      fclose(devrandom);

/*       seed = 8030655; */

      LALSCreateVector (&status, &v1, GV.ht.data->length);
      TESTSTATUS( &status );
                                                                                                                                    
      LALCreateRandomParams (&status, &randpar, seed);                                                                                                        
      TESTSTATUS( &status );
      LALNormalDeviates(&status, v1, randpar);
      TESTSTATUS( &status );
      LALDestroyRandomParams (&status, &randpar);
      TESTSTATUS( &status );
     
      for (p=0; p<(int)GV.ht.data->length; p++)
	{
	  GV.ht.data->data[p] = v1->data[p];
	}

      /* GV.ht.data->data[0]=1e7;  */
      LALSDestroyVector (&status, &v1);
      TESTSTATUS( &status );   
    }

  LALFrClose(&status,&framestream);
  TESTSTATUS( &status );

  return 0;
}




/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{
  INT4 errflg = 0;
  optarg = NULL;
  ProcessParamsTable **paramaddpoint = &procparams.processParamsTable;
  struct option long_options[] = {
    {"low-freq-cutoff",     required_argument, NULL,           'f'},
    {"bank-low-freq-cutoff",        required_argument, NULL,   'b'},
    {"threshold",           required_argument, NULL,           't'},
    {"frame-cache",         required_argument, NULL,           'F'},
    {"channel-name",        required_argument, NULL,           'C'},
    {"outfile",             required_argument, NULL,           'o'},
    {"gps-end-time",        required_argument, NULL,           'E'},
    {"gps-start-time",      required_argument, NULL,           'S'},
    {"injection-file",      required_argument, NULL,           'i'},
    {"short-segment-duration",      required_argument, NULL,   'd'},
    {"settling-time",       required_argument, NULL,           'T'},
    {"sample-rate",         required_argument, NULL,           's'},
    {"trig-start-time",     required_argument, NULL,           'g'},
    {"pad",                 required_argument, NULL,           'p'},
    {"cusp-search",                no_argument, NULL,          'c' },
    {"kink-search",                no_argument, NULL,          'k' },
    {"test-gaussian-data",         no_argument, NULL,          'n' },
    {"test-white-spectrum",        no_argument, NULL,          'w' },
    {"cluster-events",             no_argument, NULL,          'l' },
    {"help",        no_argument, NULL,         'h' },
    {0, 0, 0, 0}
  };
  char args[] = "hnckwlf:b:t:F:C:E:S:i:d:T:s:g:o:p:";

  /* set up xml output stuff */
  /* create the process and process params tables */
  procTable.processTable = LALCalloc(1, sizeof(ProcessTable));
  LALGPSTimeNow(&status, &(procTable.processTable->start_time), &accuracy);
  populate_process_table(&status, procTable.processTable, PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE);
  procparams.processParamsTable = NULL;
  /* create the search summary table */
  searchsumm.searchSummaryTable = LALCalloc(1, sizeof(SearchSummaryTable));
  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;

  /* Initialize default values */
  CLA->flow=0.0;
  CLA->FrCacheFile=NULL;
  CLA->InjectionFile=NULL;
  CLA->ChannelName=NULL;
  CLA->outputFileName=NULL;
  CLA->GPSStart=0;
  CLA->GPSEnd=0;
  CLA->ShortSegDuration=0;
  CLA->TruncSecs=0;
  CLA->power=0.0;
  CLA->fbanklow=0.0;
  CLA->threshold=0.0;
  CLA->fakenoiseflag=0;
  CLA->whitespectrumflag=0;
  CLA->whitespectrumflag=0;
  CLA->samplerate=4096.0;
  CLA->trigstarttime=0;
  CLA->cluster=0;
  CLA->pad=0;
 
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
      ADD_PROCESS_PARAM("float");
      break;
    case 's':
      /* low frequency cutoff */
      CLA->samplerate=atof(optarg);
      ADD_PROCESS_PARAM("float");
      break;
    case 'b':
      /* low frequency cutoff */
      CLA->fbanklow=atof(optarg);
      ADD_PROCESS_PARAM("float");
      break;
    case 't':
      /* low frequency cutoff */
      CLA->threshold=atof(optarg);
      ADD_PROCESS_PARAM("float");
      break;
    case 'F':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      ADD_PROCESS_PARAM("string");
      break;
    case 'C':
      /* name channel */
      CLA->ChannelName=optarg;
      memcpy(ifo, optarg, sizeof(ifo) - 1);
      ADD_PROCESS_PARAM("string");
      break;
    case 'i':
      /* name of xml injection file */
      CLA->InjectionFile=optarg;
      ADD_PROCESS_PARAM("string");
      break;
    case 'o':
      /* name of xml injection file */
      CLA->outputFileName=optarg;
      ADD_PROCESS_PARAM("string");
      break;
    case 'S':
      /* GPS start time of search */
       CLA->GPSStart=atof(optarg);
      ADD_PROCESS_PARAM("int");
      break;
    case 'E':
       /* GPS end time time of search */
      CLA->GPSEnd=atof(optarg);
      ADD_PROCESS_PARAM("int");
      break;
    case 'd':
       /* Number of segment to break-up search into */
      CLA->ShortSegDuration=atoi(optarg);
      ADD_PROCESS_PARAM("int");
      break;
    case 'T':
      /* Half the number of seconds that are trown out at the start and at the end of a short chunk */
      CLA->TruncSecs=atof(optarg);
      ADD_PROCESS_PARAM("int");
      break;
    case 'g':
      /* start time of allowed triggers */
      CLA->trigstarttime=atof(optarg);
      ADD_PROCESS_PARAM("int");
      break;
    case 'p':
      /* start time of allowed triggers */
      CLA->pad=atoi(optarg);
      ADD_PROCESS_PARAM("int");
      break;
    case 'c':
      /* cusp power law */
      CLA->power=-4.0/3.0;
      ADD_PROCESS_PARAM("string");
      break;
    case 'k':
      /* kink power law */
      CLA->power=-5.0/3.0;
      ADD_PROCESS_PARAM("string");
      break;
    case 'n':
      /* fake gaussian noise flag */
      CLA->fakenoiseflag=1;
      ADD_PROCESS_PARAM("string");
      break;
    case 'w':
      /* fake gaussian noise flag */
      CLA->whitespectrumflag=1;
      ADD_PROCESS_PARAM("string");
      break;
    case 'l':
      /* fake gaussian noise flag */
      CLA->cluster=1;
      ADD_PROCESS_PARAM("string");
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"All arguments are required except -n, -h, -w, -g, -o  and -i. One of -k or -c must be specified. They are:\n");
      fprintf(stdout,"\t--low-freq-cutoff (-f)\t\tFLOAT\t Low frequency cut-off.\n");
      fprintf(stdout,"\t--sample-rate (-s)\t\tFLOAT\t Desired sample rate (Hz).\n");
      fprintf(stdout,"\t--bank-low-freq-cutoff (-b)\tFLOAT\t Template bank low frequency cut-off.\n");
      fprintf(stdout,"\t--threshold (-t)\t\tFLOAT\t SNR threshold.\n");
      fprintf(stdout,"\t--frame-cache (-F)\t\tSTRING\t Name of frame cache file.\n");
      fprintf(stdout,"\t--channel-name (-C)\t\tSTRING\t Name of channel.\n");
      fprintf(stdout,"\t--injection-file (-i)\t\tSTRING\t Name of xml injection file.\n");
      fprintf(stdout,"\t--outfile (-o)\t\tSTRING\t Name of xml output file.\n");
      fprintf(stdout,"\t--gps-start-time (-S)\t\tINTEGER\t GPS start time.\n");
      fprintf(stdout,"\t--gps-end-time (-E)\t\tINTEGER\t GPS end time.\n");
      fprintf(stdout,"\t--settling-time (-T)\t\tINTEGER\t Number of seconds to truncate inverse square root of power spectrum.\n");
      fprintf(stdout,"\t--trig-start-time (-g)\t\tINTEGER\t GPS start time of triggers to consider.\n");
      fprintf(stdout,"\t--pad (-p)\t\tINTEGER\t Pad the data with these many seconds at beginning and end.\n");
      fprintf(stdout,"\t--short-segment-duration (-d)\t\tINTEGER\t Duration of shor segments. They will overlap by 50%s. \n","%");
      fprintf(stdout,"\t--kink-search (-k)\t\tFLAG\t Specifies a search for string kinks.\n");
      fprintf(stdout,"\t--cusp-search (-c)\t\tFLAG\t Specifies a search for string cusps.\n");
      fprintf(stdout,"\t--test-gaussian-data (-n)\tFLAG\t Use unit variance fake gaussian noise.\n");
      fprintf(stdout,"\t--test-white-spectrum (-w)\tFLAG\t Use constant white noise (used only in combination with fake gaussian noise; otherwise ignored).\n");
      fprintf(stdout,"\t--cluster-events (-l)\tFLAG\t Cluster events.\n");
      fprintf(stdout,"\t--help (-h)\t\t\tFLAG\t Print this message.\n");
      fprintf(stdout,"eg ./lalapps_StringSearch --low-freq-cutoff 40.0 --bank-low-freq-cutoff 40.0 --threshold 80.0 --frame-cache ht_local_cache --channel-name H1:Calibrated-Strain --gps-start-time 732847600 --gps-end-time 732849648 --no-of-segments 32 --settling-time 8 --cusp-search\n");
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
      fprintf(stderr,"Try ./StringSearch -h \n");
      return 1;
    }      
  if(CLA->fbanklow == 0.0)
    {
      fprintf(stderr,"No template bank low cutoff frequency specified.\n");
      fprintf(stderr,"Try ./StringSearch -h \n");
      return 1;
    }      
  if(CLA->threshold == 0.0)
    {
      fprintf(stderr,"No SNR threshold specified.\n");
      fprintf(stderr,"Try ./StringSearch -h \n");
      return 1;
    }      
  if(CLA->power == 0.0)
    {
      fprintf(stderr,"Cusp or kink search not specified. \n");
      fprintf(stderr,"Try ./StringSearch -h \n");
      return 1;
    }      
  if(CLA->FrCacheFile == NULL)
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try ./StringSearch -h \n");
      return 1;
    }      
  if(CLA->ChannelName == NULL)
    {
      fprintf(stderr,"No channel name specified.\n");
      fprintf(stderr,"Try ./StringSearch -h \n");
      return 1;
    }      
  if(CLA->GPSStart == 0)
    {
      fprintf(stderr,"No GPS start time specified.\n");
      fprintf(stderr,"Try ./StringSearch -h \n");
      return 1;
    }      
  if(CLA->GPSEnd == 0)
    {
      fprintf(stderr,"No GPS end time specified.\n");
      fprintf(stderr,"Try ./StringSearch -h \n");
      return 1;
    }      
  if(CLA->ShortSegDuration == 0)
    {
      fprintf(stderr,"Short segment duration not specified (they overlap by 50%s).\n","%");
      fprintf(stderr,"Try ./StringSearch -h \n");
      return 1;
    }      

  /* Some consistency checking */
  {
    int big_seg_length=CLA->GPSEnd-CLA->GPSStart-2*CLA->pad;
    int small_seg_length=CLA->ShortSegDuration;

    REAL4 x=((float)big_seg_length/(float)small_seg_length)-0.5;

    if ((int)x != x)
      {
        fprintf(stderr,"The total duration of the segment T and the short segment duration\n");
        fprintf(stderr,"Should obey the following rule: T/t - 0.5 shold be an odd integer.\n");
	return 1;
      } 
    if (((int)x)%2 != 1)
      {
        fprintf(stderr,"The total duration of the segment T and the short segment duration\n");
        fprintf(stderr,"Should obey the following rule: T/t - 0.5 shold be an odd integer.\n");
	return 1;
      }     
    if( CLA->ShortSegDuration/4  < CLA->TruncSecs)
      {
        fprintf(stderr,"Short segment length t=%d is too small to accomodate truncation time requested.\n",
		small_seg_length);
	fprintf(stderr,"Need short segment t(=%d) to be >= 4 x Truncation length (%f).\n",CLA->ShortSegDuration,CLA->TruncSecs);
	return 1;
      }    
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
      {
	searchsumm.searchSummaryTable->out_start_time.gpsSeconds = CLA->trigstarttime;
      }else{
	searchsumm.searchSummaryTable->out_start_time.gpsSeconds = CLA->GPSStart+small_seg_length/4+CLA->pad;
      }
    searchsumm.searchSummaryTable->out_start_time.gpsNanoSeconds =0;
    searchsumm.searchSummaryTable->out_end_time.gpsSeconds = CLA->GPSEnd-small_seg_length/4-CLA->pad;
    searchsumm.searchSummaryTable->out_end_time.gpsNanoSeconds =0;
  }

  return errflg;
}

/*******************************************************************************/

int FreeMem(void)
{


  LALSDestroyVector(&status,&GV.ht_proc.data);
  TESTSTATUS( &status );

  LALSDestroyVector(&status,&GV.Spec.data);
  TESTSTATUS( &status );

  LALSDestroyVector(&status,&GV.StringFilter.data);
  TESTSTATUS( &status );

  LALDestroyRealFFTPlan( &status, &GV.fplan );
  TESTSTATUS( &status );

  LALDestroyRealFFTPlan( &status, &GV.rplan );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
 
  return 0;
}

/*******************************************************************************/
#endif
