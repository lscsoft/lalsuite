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
/*            X. Siemens, J. Creighton, F. Robinet, and K. Cannon                */
/*                                                                               */
/*                             UWM/LAL - December 2009                           */
/*********************************************************************************/

#include <config.h>

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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
#include <lal/Date.h>
#include <lal/Units.h>

#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataBurstUtils.h>

#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLBurstRead.h>

#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/GenerateBurst.h>

#include <lalapps.h>
#include <processtable.h>
#include <LALAppsVCSInfo.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

#define SCALE 1e20
#define MAXTEMPLATES 50

NRCSID( STRINGSEARCHC, "StringSearch $Id$");
RCSID( "StringSearch $Id$");

/* FIXME:  should be "lalapps_StringSearch" to match the executable */
/* requires post-processing codes to be updated */
#define PROGRAM_NAME "StringSearch"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"


/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 samplerate;           /* desired sample rate */
  REAL8 fbankstart;           /* lowest frequency of templates */
  REAL8 fbankhighfcutofflow;  /* lowest high frequency cut-off */
  REAL8 fmismatchmax;         /* maximal mismatch allowed from 1 template to the next */
  char *FrCacheFile;          /* Frame cache file */
  char *InjectionFile;        /* LIGO/Virgo xml injection file */
  char *ChannelName;          /* Name of channel to be read in from frames */
  char *outputFileName;       /* Name of xml output filename */
  LIGOTimeGPS GPSStart;       /* GPS start time of segment to be analysed */
  LIGOTimeGPS GPSEnd;         /* GPS end time of segment to be analysed */
  INT4 ShortSegDuration;      /* Number of fixed length sub-segments between GPSStart and GPSEnd */
  REAL8 TruncSecs;            /* Half the number of seconds truncated at beginning and end of a chunk */
  REAL8 power;                /* Kink (-5/3) or cusp (-4/3) frequency power law */
  REAL8 threshold;            /* event SNR threshold */
  INT4 fakenoiseflag;         /* =0 if real noise =1 if fake gaussian noise */
  INT4 whitespectrumflag;     /* =0 if spectrum is to be computed =1 for white spectrum */
  LIGOTimeGPS trigstarttime;  /* GPS start time of allowed triggers */
  REAL8 cluster;              /* =0.0 if events are not to be clustered = clustering time otherwise */
  INT4 pad;                   /* seconds of padding */
  double chi2cut[3];          /* chi2 cut parameters */
  INT4 printspectrumflag;     /* flag set to 1 if user wants to print the spectrum */
  INT4 printfilterflag;       /* flag set to 1 if user wants to print the filter in the frequency domain */
  INT4 printfirflag;          /* flag set to 1 if user wants to print the filter in the time domain */
  INT4 printsnrflag;          /* flag set to 1 if user wants to print the snr */
  INT4 printdataflag;         /* flag set to 1 if user wants to print the data */  
  INT4 printinjectionflag;    /* flag set to 1 if user wants to print the injection(s) */  
  char *comment;              /* for "comment" columns in some tables */
} CommandLineArgs;

typedef 
struct GlobalVariablesTag {
  REAL8FrequencySeries *Spec; /* average spectrum */
  REAL8FFTPlan *fplan;        /* fft plan */
  REAL8FFTPlan *rplan;        /* fft plan */
  unsigned seg_length;
} GlobalVariables;

typedef
struct StringTemplateTag {
  INT4 findex;                /* Template frequency index */
  REAL8 f;                    /* Template frequency */
  REAL8 norm;                 /* Template normalisation */
  REAL8 mismatch;             /* Template mismatch relative to last one */
  REAL8FrequencySeries *StringFilter; /* Frequency domain filter corresponding to this template */
  REAL8Vector *waveform_t;    /* Template waveform - time-domain */
  COMPLEX16Vector *waveform_f; /* Template waveform - frequency domain */
  REAL8Vector *auto_cor;      /* Auto-correlation vector */
  INT4 chi2_index;            /* index to compute chi2 */
} StringTemplate;

/***************************************************************************/

/* GLOBAL VARIABLES */
GlobalVariables GV;           /* A bunch of stuff is stored in here; mainly to protect it from accidents */


/***************************************************************************/

/* FUNCTION PROTOTYPES */

/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA, const ProcessTable *process, ProcessParamsTable **paramaddpoint);

/* Reads raw data (or puts in fake gaussian noise with a sigma=10^-20) */
REAL8TimeSeries *ReadData(struct CommandLineArgsTag CLA);

/* Adds injections if an xml injection file is given */
int AddInjections(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht);

/* DownSamples data */
int DownSample(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht);

/* Computes the average spectrum  */
int AvgSpectrum(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht);

/* Creates the template bank based on the spectrum  */
int CreateTemplateBank(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, StringTemplate *strtemplate, int *NTemplates);

/* Creates the frequency domain string cusp or kink filters  */
int CreateStringFilters(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, StringTemplate *strtemplate, int NTemplates);

/* Filters the data through the template banks  */
int FindStringBurst(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, const StringTemplate *strtemplate, int NTemplates, SnglBurst **head);

/* Finds events above SNR threshold specified  */
int FindEvents(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, const StringTemplate *strtemplate,
               const REAL8Vector *vector, INT4 i, SnglBurst **head);

/* Writes out the xml file with the events it found  */
int OutputEvents(const struct CommandLineArgsTag *CLA, ProcessTable *proctable, ProcessParamsTable *procparamtable, SnglBurst *events);

/* Frees the memory */
int FreeMem(StringTemplate *strtemplate, int NTemplates);                                        

/* Clustering comparison function */
static int XLALCompareStringBurstByTime(const SnglBurst * const *, const SnglBurst * const *);

/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  StringTemplate strtemplate[MAXTEMPLATES];
  int NTemplates;
  SnglBurst *events=NULL;
  MetadataTable  process;
  MetadataTable  procparams;
  REAL8TimeSeries *ht;        /* strain data */

  lalDebugLevel = LALINFO | LALWARNING | LALERROR | LALNMEMDBG | LALNMEMPAD | LALNMEMTRK;

  /* create the process and process params tables */
  procparams.processParamsTable = NULL;
  process.processTable = XLALCreateProcessTableRow();
  XLALGPSTimeNow(&(process.processTable->start_time));
  if(XLALPopulateProcessTable(process.processTable, PROGRAM_NAME, LALAPPS_VCS_IDENT_ID, LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 0))
    exit(1);

  /****** ReadCommandLine ******/
  XLALPrintInfo("ReadCommandLine()\n");
  if (ReadCommandLine(argc,argv,&CommandLineArgs, process.processTable, &procparams.processParamsTable)) return 1;
  XLALPrintInfo("\t%c%c detector\n",CommandLineArgs.ChannelName[0],CommandLineArgs.ChannelName[1]);

  /****** ReadData ******/
  XLALPrintInfo("ReadData()\n");
  ht = ReadData(CommandLineArgs);
  if (!ht) return 4;

  /****** AddInjections ******/
  if (CommandLineArgs.InjectionFile != NULL) {
    XLALPrintInfo("AddInjections()\n");
    if (AddInjections(CommandLineArgs, ht)) return 5;
    if ( CommandLineArgs.printinjectionflag ) LALDPrintTimeSeries( ht, "injection.txt" );
  }

  if ( CommandLineArgs.printdataflag ){
    unsigned p;
    for (p=0; p<ht->data->length; p++)
      fprintf(stdout,"%1.15e\n",ht->data->data[p]);
    return 0;
  }

  /****** DownSample ******/
  XLALPrintInfo("DownSample()\n");
  if (DownSample(CommandLineArgs, ht)) return 8;

  /****** XLALResizeREAL8TimeSeries ******/
  XLALPrintInfo("XLALResizeREAL8TimeSeries()\n");	
  /* re-size the time series to remove the pad */
  ht = XLALResizeREAL8TimeSeries(ht, 
					 (int)round(CommandLineArgs.pad/ht->deltaT),
					 ht->data->length-2*(int)round(CommandLineArgs.pad/ht->deltaT));

  /****** AvgSpectrum ******/
  XLALPrintInfo("AvgSpectrum()\n");
  if (AvgSpectrum(CommandLineArgs, ht)) return 9;  
  if (CommandLineArgs.printspectrumflag) LALDPrintFrequencySeries( GV.Spec, "Spectrum.txt" );

  /****** CreateTemplateBank ******/
  XLALPrintInfo("CreateTemplateBank()\n");
  if (CreateTemplateBank(CommandLineArgs, ht, strtemplate, &NTemplates)) return 10;

  /****** CreateStringFilters ******/
  XLALPrintInfo("CreateStringFilters()\n");
  if (CreateStringFilters(CommandLineArgs, ht, strtemplate, NTemplates)) return 11;

  /****** FindStringBurst ******/
  XLALPrintInfo("FindStringBurst()\n");
  if (FindStringBurst(CommandLineArgs, ht, strtemplate, NTemplates, &events)) return 12;
  if(!XLALSortSnglBurst(&events, XLALCompareSnglBurstByExactPeakTime)) return 12;
  XLALDestroyREAL8TimeSeries(ht);
  ht = NULL;

  /****** XLALClusterSnglBurstTable ******/
  XLALPrintInfo("XLALClusterSnglBurstTable()\n");
  if (CommandLineArgs.cluster != 0.0 && events){
    XLALClusterSnglBurstTable(&events, XLALCompareStringBurstByTime, XLALCompareStringBurstByTime, XLALStringBurstCluster);
    XLALSortSnglBurst(&events, XLALCompareSnglBurstByPeakTimeAndSNR);
  }

  /****** XLALSnglBurstAssignIDs ******/
  XLALPrintInfo("XLALSnglBurstAssignIDs()\n");
  XLALSnglBurstAssignIDs(events, process.processTable->process_id, 0);

  /****** OutputEvents ******/
  XLALPrintInfo("OutputEvents()\n");
  if (OutputEvents(&CommandLineArgs, process.processTable, procparams.processParamsTable, events)) return 13;

  /****** FreeMem ******/
  XLALPrintInfo("FreeMem()\n");
  XLALDestroyProcessParamsTable(procparams.processParamsTable);
  XLALDestroySnglBurstTable(events);
  XLALDestroyProcessTable(process.processTable);
  if (FreeMem(strtemplate, NTemplates)) return 14;

  XLALPrintInfo("StringJob is done\n");

  return 0;
}


/**************************** MAIN PROGRAM ENDS ********************************/

/*******************************************************************************/

/*
 * Check if two string events overlap in time. The peak times are uncertain
 * to whatever the high frequency cutoff is
 */

/* <lalVerbatim file="SnglBurstUtilsCP"> */
static int XLALCompareStringBurstByTime(
  const SnglBurst * const *a,
  const SnglBurst * const *b
)
/* </lalVerbatim> */
{
  double delta_t = XLALGPSDiff(&(*a)->peak_time, &(*b)->peak_time);
  /* FIXME:  global variables = BAD BAD BAD! (my fault -- Kipp) */
  double epsilon = CommandLineArgs.cluster;

  if(delta_t > epsilon)
    return(1);
  if(delta_t < -epsilon)
    return(-1);
  return(0);
}

/*******************************************************************************/

int AddInjections(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht){
  LIGOTimeGPS startTime = ht->epoch;
  LIGOTimeGPS stopTime = ht->epoch;
  SimBurst *sim_burst;

  XLALGPSAdd(&stopTime, ht->data->length * ht->deltaT);

  /* Get info from injection file */
  sim_burst = XLALSimBurstTableFromLIGOLw(CLA.InjectionFile, &startTime, &stopTime);

  /* Inject the signals into ht */
  if(XLALBurstInjectSignals(ht, sim_burst, NULL))
    return 1;

  /* free the injection table */
  XLALDestroySimBurstTable(sim_burst);

  return 0;
}

/*******************************************************************************/

static ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param,
					      const ProcessTable *process,
					      const char *type, const char *param, const char *value){
  *proc_param = XLALCreateProcessParamsTableRow(process);
  snprintf((*proc_param)->program, LIGOMETA_PROGRAM_MAX, PROGRAM_NAME);
  snprintf((*proc_param)->type, LIGOMETA_TYPE_MAX, "%s", type);
  snprintf((*proc_param)->param, LIGOMETA_PARAM_MAX, "--%s", param);
  snprintf((*proc_param)->value, LIGOMETA_VALUE_MAX, "%s", value);
  
  return(&(*proc_param)->next);
}

#define ADD_PROCESS_PARAM(process, type) \
	do { paramaddpoint = add_process_param(paramaddpoint, process, type, long_options[option_index].name, optarg); } while(0)

/*******************************************************************************/

int OutputEvents(const struct CommandLineArgsTag *CLA, ProcessTable *proctable, ProcessParamsTable *procparamtable, SnglBurst *events){
  LIGOLwXMLStream *xml;
  MetadataTable  searchsumm;
  char ifo[3];

  strncpy( ifo, CLA->ChannelName, 2 );
  ifo[sizeof(ifo) - 1] = '\0';

  if (!CLA->outputFileName){
    CHAR outfilename[256];
    snprintf(outfilename, sizeof(outfilename)-1, "%s-STRINGSEARCH-%d-%d.xml", ifo,
	     searchsumm.searchSummaryTable->in_start_time.gpsSeconds,
	     searchsumm.searchSummaryTable->in_end_time.gpsSeconds - 
	     searchsumm.searchSummaryTable->in_start_time.gpsSeconds);
    outfilename[sizeof(outfilename)-1] = '\0';
    xml = XLALOpenLIGOLwXMLFile(outfilename);
  }
  else
    xml = XLALOpenLIGOLwXMLFile(CLA->outputFileName);

  /* process table */
  snprintf(proctable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo);
  XLALGPSTimeNow(&(proctable->end_time));
  
  if(XLALWriteLIGOLwXMLProcessTable(xml, proctable)) return -1;
  
  /* process params table */
  if(XLALWriteLIGOLwXMLProcessParamsTable(xml, procparamtable)) return -1;
  
  /* search summary table */
  /* create the search summary table */
  searchsumm.searchSummaryTable = XLALCreateSearchSummaryTableRow(proctable);
  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;
  /* store the input and output start and end times */
  searchsumm.searchSummaryTable->in_start_time = CLA->GPSStart;
  searchsumm.searchSummaryTable->in_end_time = CLA->GPSEnd;
  if (XLALGPSToINT8NS(&CLA->trigstarttime) > 0)
    searchsumm.searchSummaryTable->out_start_time = CLA->trigstarttime;
  else {
    searchsumm.searchSummaryTable->out_start_time = CLA->GPSStart;
    XLALGPSAdd(&searchsumm.searchSummaryTable->out_start_time, CLA->ShortSegDuration/4+CLA->pad);
  }
  searchsumm.searchSummaryTable->out_end_time = CLA->GPSEnd;
  XLALGPSAdd(&searchsumm.searchSummaryTable->out_end_time, -CLA->ShortSegDuration/4-CLA->pad);
  snprintf(searchsumm.searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo);
  searchsumm.searchSummaryTable->nevents = XLALSnglBurstTableLength(events);
  if(XLALWriteLIGOLwXMLSearchSummaryTable(xml, searchsumm.searchSummaryTable)) return -1;
  XLALDestroySearchSummaryTable(searchsumm.searchSummaryTable);

  /* burst table */
  if(XLALWriteLIGOLwXMLSnglBurstTable(xml, events)) return -1;
  
  XLALCloseLIGOLwXMLFile(xml);

  return 0;
}

/*******************************************************************************/

int FindEvents(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, const StringTemplate *strtemplate, const REAL8Vector *vector, INT4 i, SnglBurst **head){
  unsigned p;
  REAL8 maximum, chi2, ndof;
  REAL8 duration;
  INT4 pmax, pend, pstart;

  /* print the snr to stdout */
  if (CLA.printsnrflag)
    for ( p = vector->length/4 ; p < 3*vector->length/4; p++ )
      fprintf(stdout,"%p %e\n", strtemplate, vector->data[p]);
  
  /* Now find event in the inner half */
  for ( p = vector->length/4 ; p < 3*vector->length/4; p++ ){
    SnglBurst *new;
    LIGOTimeGPS peaktime, starttime;
    LIGOTimeGPS t;

    maximum = 0.0;
    pmax=p;
    t = ht->epoch;
    XLALGPSAdd(&t, (GV.seg_length*i/2 + p) * ht->deltaT);

    /* Do we have the start of a cluster? */
    if ( (fabs(vector->data[p]) > CLA.threshold) && (XLALGPSCmp(&t, &CLA.trigstarttime) >= 0)){
      int pp;
      pend=p; pstart=p;

      t = ht->epoch;
      XLALGPSAdd(&t, GV.seg_length*i/2*ht->deltaT);

      /* Clustering in time: While we are above threshold, or within clustering time of the last point above threshold... */
      while( ((fabs(vector->data[p]) > CLA.threshold) || ((p-pend)* ht->deltaT < (float)(CLA.cluster)) ) 
	     && p<3*vector->length/4){
	
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

      starttime = peaktime = t;
      XLALGPSAdd(&peaktime, ht->deltaT * pmax);
      XLALGPSAdd(&starttime, ht->deltaT * pstart);
      duration = ht->deltaT * ( pend - pstart );

      /* compute \chi^{2} */
      chi2=0, ndof=0;
      for(pp=-strtemplate->chi2_index; pp<strtemplate->chi2_index; pp++){
        chi2 += (vector->data[pmax+pp]-vector->data[pmax]*strtemplate->auto_cor->data[GV.seg_length/2+pp])*(vector->data[pmax+pp]-vector->data[pmax]*strtemplate->auto_cor->data[GV.seg_length/2+pp]);
        ndof += (1-strtemplate->auto_cor->data[GV.seg_length/2+pp]*strtemplate->auto_cor->data[GV.seg_length/2+pp]);
      }

 
      /* Apply the \chi^{2} cut */
      if( CLA.chi2cut[0]    > -9999
	  && CLA.chi2cut[1] > -9999
	  && CLA.chi2cut[2] > -9999 )
	if(log10(chi2/ndof)>CLA.chi2cut[0]
	   && log10(chi2/ndof)> CLA.chi2cut[1]*log10(fabs(maximum))+CLA.chi2cut[2]) continue;

      /* prepend a new event to the linked list */
      new = XLALCreateSnglBurst();
      if ( ! new ){ /* allocation error */
	XLALPrintError("Could not allocate memory for event. Memory allocation error. Exiting.\n");
	return 1;
      }
      new->next = *head;
      *head = new;

      /* Now copy stuff into event */
      strncpy( new->ifo, CLA.ChannelName, 2 );
      new->ifo[3] = 0;
      strncpy( new->search, "StringCusp", sizeof( new->search ) );
      strncpy( new->channel, CLA.ChannelName, sizeof( new->channel ) );
      
      /* give trigger a 1 sample fuzz on either side */
      XLALGPSAdd(&starttime, -ht->deltaT);
      duration += 2 * ht->deltaT;

      new->start_time = starttime;
      new->peak_time = peaktime;
      new->duration     = duration;
      new->central_freq = (strtemplate->f+CLA.fbankstart)/2.0;	   
      new->bandwidth    = strtemplate->f-CLA.fbankstart;				     
      new->snr          = maximum;
      new->amplitude   = vector->data[pmax]/strtemplate->norm;
      new->chisq = chi2;
      new->chisq_dof = ndof;
    }
  }
    
  return 0;
}

/*******************************************************************************/

int FindStringBurst(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, const StringTemplate *strtemplate, int NTemplates, SnglBurst **head){
  int i,m;
  unsigned p;
  REAL8Vector *vector = NULL;
  COMPLEX16Vector *vtilde = NULL;

  /* create vector that will hold the data for each overlapping chunk */ 
  vector = XLALCreateREAL8Vector( GV.seg_length);
  
  /* create vector that will hold FFT of data*/
  vtilde = XLALCreateCOMPLEX16Vector( GV.seg_length / 2 + 1 );
  
  /* loop over templates  */
  for (m = 0; m < NTemplates; m++){
    /* loop over overlapping chunks */ 
    for(i=0; i < 2*(ht->data->length*ht->deltaT)/CLA.ShortSegDuration - 1 ;i++){
      /* populate vector that will hold the data for each overlapping chunk */
      memcpy( vector->data, ht->data->data + i*GV.seg_length/2,vector->length*sizeof( *vector->data ) );
	  
      /* fft it */
      if(XLALREAL8ForwardFFT( vtilde, vector, GV.fplan )) return 1;
      
      /* multiply FT of data and String Filter and deltaT */
      for ( p = 0 ; p < vtilde->length; p++ ){
	vtilde->data[p].re *= strtemplate[m].StringFilter->data->data[p]*ht->deltaT;
	vtilde->data[p].im *= strtemplate[m].StringFilter->data->data[p]*ht->deltaT;
      }
      
      if(XLALREAL8ReverseFFT( vector, vtilde, GV.rplan )) return 1;

      /* normalise the result by template normalisation and multiply by 
	 df (not inluded in LALReverseRealFFT)  factor of 2 is from 
	 match-filter definition */

      for ( p = 0 ; p < vector->length; p++ )
	vector->data[p] *= 2.0 * GV.Spec->deltaF / strtemplate[m].norm;

      if(FindEvents(CLA, ht, &strtemplate[m], vector, i, head)) return 1;
    }
  }

  XLALDestroyCOMPLEX16Vector( vtilde );
  XLALDestroyREAL8Vector( vector );

  return 0;
}


/*******************************************************************************/

int CreateStringFilters(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, StringTemplate *strtemplate, int NTemplates){
  int m;
  unsigned p;
  COMPLEX16Vector *vtilde; /* frequency-domain vector workspace */
  REAL8Vector    *vector; /* time-domain vector workspace */
  REAL8 re, im;

  vector = XLALCreateREAL8Vector( GV.seg_length);
  vtilde = XLALCreateCOMPLEX16Vector( GV.seg_length / 2 + 1 );
 
  for (m = 0; m < NTemplates; m++){

    /* Initialize the filter */
    strtemplate[m].StringFilter = XLALCreateREAL8FrequencySeries(CLA.ChannelName, &CLA.GPSStart, 0, 0, &lalStrainUnit, GV.Spec->data->length);
    strtemplate[m].StringFilter->deltaF=GV.Spec->deltaF;

    /* populate vtilde with the template divided by the noise */
    for ( p = 0; p < vtilde->length; p++ ){
      vtilde->data[p].re = sqrt(strtemplate[m].waveform_f->data[p].re/(GV.Spec->data->data[p]));
      vtilde->data[p].im = sqrt(strtemplate[m].waveform_f->data[p].im/(GV.Spec->data->data[p]));
    }

    /* reverse FFT vtilde into vector */
    if(XLALREAL8ReverseFFT( vector, vtilde, GV.rplan )) return 1;

    /* multiply times df to make sure units are correct */
    for ( p = 0 ; p < vector->length; p++ )
      vector->data[p] *= GV.Spec->deltaF;

    /* perform the truncation; the truncation is CLA.TruncSecs/2 because 
       we are dealing with the sqrt of the filter at the moment*/
    if(CLA.TruncSecs != 0.0)
      memset( vector->data + (INT4)(CLA.TruncSecs/2/ht->deltaT +0.5), 0,
	      ( vector->length -  2 * (INT4)(CLA.TruncSecs/2/ht->deltaT +0.5)) 
	      * sizeof( *vector->data ) );

    /* forward fft the truncated vector into vtilde */
    if(XLALREAL8ForwardFFT( vtilde, vector, GV.fplan )) return 1;

    for ( p = 0 ; p < vtilde->length-1; p++ ){
      re = vtilde->data[p].re * ht->deltaT;
      im = vtilde->data[p].im * ht->deltaT;
      strtemplate[m].StringFilter->data->data[p] = (re * re + im * im);
    }

    /* set DC and Nyquist to 0*/
    strtemplate[m].StringFilter->data->data[0] =
      strtemplate[m].StringFilter->data->data[vtilde->length-1] = 0;

    /* print out the frequency domain filter */
    if (CLA.printfilterflag){
      CHAR filterfilename[256];
      snprintf(filterfilename, sizeof(filterfilename)-1, "Filter-%d.txt", m);
      filterfilename[sizeof(filterfilename)-1] = '\0';
      LALDPrintFrequencySeries( strtemplate[m].StringFilter, filterfilename );
    }

    /* print out the time domain FIR filter */
    if (CLA.printfirflag){
      REAL8TimeSeries series;
      CHAR filterfilename[256];
      series.deltaT=ht->deltaT;
      series.f0 = 0.0;
      strncpy(series.name, "fir filter", LALNameLength);
      series.epoch=ht->epoch;
      series.sampleUnits=ht->sampleUnits;

      for ( p = 0 ; p < vtilde->length-1; p++ ){
	re = vtilde->data[p].re * ht->deltaT;
	im = vtilde->data[p].im * ht->deltaT;

	vtilde->data[p].re = (re * re + im * im);
	vtilde->data[p].im = 0.0;
      }
      vtilde->data[0].re = vtilde->data[0].im = 0.0;
      vtilde->data[vtilde->length-1].re = vtilde->data[vtilde->length-1].im =0;

      if(XLALREAL8ReverseFFT( vector, vtilde, GV.rplan )) return 1;

      series.data = vector;

      snprintf(filterfilename, sizeof(filterfilename)-1, "FIRFilter-%d.txt", m);
      filterfilename[sizeof(filterfilename)-1] = '\0';
      LALDPrintTimeSeries( &series, filterfilename );
    }
  }

  XLALDestroyCOMPLEX16Vector( vtilde );
  XLALDestroyREAL8Vector( vector );

  return 0;
}

/*******************************************************************************/

int CreateTemplateBank(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, StringTemplate *strtemplate, int *NTemplates){
  REAL8 fNyq, f_cut, f, t1t1, t2t2, t1t2, epsilon, previous_epsilon, norm, slope0, slope1;
  int m, f_cut_index, f_low_cutoff_index, extr_ctr;
  unsigned p, pcut, f_min_index, f_max_index;
  REAL8Vector *integral;
  REAL8Vector    *vector; /* time-domain vector workspace */
  COMPLEX16Vector *vtilde; /* frequency-domain vector workspace */

  fNyq = (1.0/ht->deltaT) / 2.0;
  f_min_index = round(CLA.fbankstart / GV.Spec->deltaF);
  f_max_index = round(fNyq / GV.Spec->deltaF);
  integral = XLALCreateREAL8Vector(f_max_index-f_min_index);
  epsilon=0;

  /* first template : f_cutoff = fbankhighfcutofflow */
  f_cut_index = round(CLA.fbankhighfcutofflow / GV.Spec->deltaF);
  f_cut = f_cut_index * GV.Spec->deltaF;

  /* compute (t1|t1) */
  t1t1=0.0;
  integral->data[0]=4*pow( pow(CLA.fbankstart,CLA.power),2)/GV.Spec->data->data[f_min_index]*GV.Spec->deltaF;
  for( p = f_min_index ; p < f_max_index; p++ ){
    f = p*GV.Spec->deltaF;

    if(f<=f_cut) t1t1 += 4*pow(pow(f,CLA.power),2)/GV.Spec->data->data[p]*GV.Spec->deltaF;
    else t1t1 += 4*pow( pow(f,CLA.power)*exp(1-f/f_cut) ,2)/GV.Spec->data->data[p]*GV.Spec->deltaF;

    if(p>f_min_index) /* keep the integral in memory (to run faster) */
      integral->data[p-f_min_index] = integral->data[p-f_min_index-1]+4*pow(pow(f,CLA.power),2)/GV.Spec->data->data[p]*GV.Spec->deltaF;
  }

  strtemplate[0].findex=f_cut_index;
  strtemplate[0].f=f_cut;
  strtemplate[0].mismatch=0.0;
  strtemplate[0].norm=sqrt(t1t1);
  *NTemplates=1;
  XLALPrintInfo("%% Templ. frequency      sigma      mismatch\n");  
  XLALPrintInfo("%% %d      %1.3e    %1.3e    %1.3e\n",*NTemplates-1,strtemplate[0].f,strtemplate[0].norm, strtemplate[0].mismatch);
  
  /* find the next cutoffs given the maximal mismatch */
  for(pcut=f_cut_index+1; pcut<f_max_index; pcut++){
    f_cut = pcut*GV.Spec->deltaF;
   
    t2t2=integral->data[strtemplate[*NTemplates-1].findex-f_min_index];
    t1t2=integral->data[strtemplate[*NTemplates-1].findex-f_min_index];
    
    /* compute (t2|t2) and (t1|t2) */
    for( p = strtemplate[*NTemplates-1].findex+1 ; p < f_max_index; p++ ){
      f = p*GV.Spec->deltaF;
      
      /* (t2|t2) */
      if(f<=f_cut)
	t2t2 += 4*pow(pow(f,CLA.power),2)/GV.Spec->data->data[p]*GV.Spec->deltaF;
      else 
	t2t2 += 4*pow( pow(f,CLA.power)*exp(1-f/f_cut) ,2)/GV.Spec->data->data[p]*GV.Spec->deltaF;

      /* (t1|t2) */
      if(f<=f_cut)
	t1t2 += 4*pow(pow(f,CLA.power),2)*exp(1-f/strtemplate[*NTemplates-1].f) /GV.Spec->data->data[p]*GV.Spec->deltaF;
      else 
	t1t2 += 4*pow( pow(f,CLA.power),2)*exp(1-f/strtemplate[*NTemplates-1].f)*exp(1-f/f_cut) /GV.Spec->data->data[p]*GV.Spec->deltaF;
    }
        
    previous_epsilon = epsilon;
    epsilon=1-t1t2/sqrt(t1t1*t2t2);

    /*if(pcut%50==0) XLALPrintInfo("%d %f %f\n",pcut, f_cut, epsilon);*/

    if(epsilon >= CLA.fmismatchmax || pcut==f_max_index-1){
      strtemplate[*NTemplates].findex=pcut;
      strtemplate[*NTemplates].f=f_cut;
      strtemplate[*NTemplates].norm=sqrt(t2t2);
      strtemplate[*NTemplates].mismatch=epsilon;
      (*NTemplates)++;
      XLALPrintInfo("%% %d      %1.3e    %1.3e    %1.3e\n",*NTemplates-1,strtemplate[*NTemplates-1].f,strtemplate[*NTemplates-1].norm, strtemplate[*NTemplates-1].mismatch);
      t1t1=t2t2;
      if(*NTemplates == MAXTEMPLATES){
	XLALPrintError("Too many templates for code... Exiting\n");
	return 1;
      }
    }
    
    /* to get faster (not so smart though) */
    if(pcut<f_max_index-16 && (epsilon-previous_epsilon)<0.005) 
      pcut+=15;

  }

  XLALDestroyREAL8Vector( integral );

  /* Now, the point is to store the template waveform vector */
  vector = XLALCreateREAL8Vector( GV.seg_length);
  vtilde = XLALCreateCOMPLEX16Vector( GV.seg_length / 2 + 1 );
  f_low_cutoff_index = (int) (CLA.fbankstart/ GV.Spec->deltaF+0.5);
  for (m = 0; m < *NTemplates; m++){
    
    /* create the space for the waveform vectors */
    strtemplate[m].waveform_f = XLALCreateCOMPLEX16Vector( GV.seg_length / 2 + 1 );
    strtemplate[m].waveform_t = XLALCreateREAL8Vector( GV.seg_length);
    strtemplate[m].auto_cor   = XLALCreateREAL8Vector( GV.seg_length);
    
    /* populate with the template waveform */
    for ( p = f_low_cutoff_index; p < strtemplate[m].waveform_f->length; p++ ){
      f=p*GV.Spec->deltaF;
      if(f<=strtemplate[m].f) 
	strtemplate[m].waveform_f->data[p].re = pow(f,CLA.power);
      else 
	strtemplate[m].waveform_f->data[p].re = pow(f,CLA.power)*exp(1-f/strtemplate[m].f);
      strtemplate[m].waveform_f->data[p].im = 0;
    }
    
    /* set all frequencies below the low freq cutoff to zero */
    memset(strtemplate[m].waveform_f->data, 0, f_low_cutoff_index*sizeof(*strtemplate[m].waveform_f->data));
    
    /* set DC and Nyquist to zero anyway */
    strtemplate[m].waveform_f->data[0].re = strtemplate[m].waveform_f->data[strtemplate[m].waveform_f->length - 1].re = 0;
    strtemplate[m].waveform_f->data[0].im = strtemplate[m].waveform_f->data[strtemplate[m].waveform_f->length - 1].im = 0;
    
    for (p=0 ; p< vtilde->length; p++){
      vtilde->data[p].re = strtemplate[m].waveform_f->data[p].re*strtemplate[m].waveform_f->data[p].re/GV.Spec->data->data[p];
      vtilde->data[p].im = 0;
    }
    
    /* reverse FFT */
    if(XLALREAL8ReverseFFT(vector, strtemplate[m].waveform_f, GV.rplan)) return 1;
    if(XLALREAL8ReverseFFT(strtemplate[m].auto_cor, vtilde, GV.rplan)) return 1;

    /* The vector is reshuffled in the right order */
    for ( p = 0 ; p < GV.seg_length/2; p++ ){
      strtemplate[m].waveform_t->data[p] = vector->data[GV.seg_length/2+p]*GV.Spec->deltaF;;
      strtemplate[m].waveform_t->data[GV.seg_length/2+p] = vector->data[p]*GV.Spec->deltaF;;
    }

    /* Normalize the autocorrelation by the central value */
    norm=strtemplate[m].auto_cor->data[0];
    for ( p = 0 ; p < strtemplate[m].auto_cor->length; p++ ){
      strtemplate[m].auto_cor->data[p] /= norm;
      vector->data[p]=strtemplate[m].auto_cor->data[p];
    }

    /* The vector is reshuffled in the right order */
    for ( p = 0 ; p < GV.seg_length/2; p++ ){
      strtemplate[m].auto_cor->data[p] = vector->data[GV.seg_length/2+p];
      strtemplate[m].auto_cor->data[GV.seg_length/2+p] = vector->data[p];
    }

    /* search for the index of the 3rd extremum */
    extr_ctr=0;
    strtemplate[m].chi2_index=0;
    for ( p = GV.seg_length/2+1; p< GV.seg_length-1; p++ ){

      slope1 = strtemplate[m].waveform_t->data[p+1]-strtemplate[m].waveform_t->data[p];
      slope0 = strtemplate[m].waveform_t->data[p]-strtemplate[m].waveform_t->data[p-1];
      strtemplate[m].chi2_index++;
      if(slope0*slope1<0){
	extr_ctr++;
	if(extr_ctr==2) break;
      }

    }
      

  }
  
  XLALDestroyREAL8Vector( vector );
  XLALDestroyCOMPLEX16Vector( vtilde );

  return 0;
}


/*******************************************************************************/
int AvgSpectrum(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht){
  GV.seg_length = (int)(CLA.ShortSegDuration/ht->deltaT + 0.5);
  GV.Spec  = XLALCreateREAL8FrequencySeries(CLA.ChannelName, &CLA.GPSStart, 0, 0, &lalStrainUnit, GV.seg_length / 2 + 1);
  GV.fplan = XLALCreateForwardREAL8FFTPlan( GV.seg_length, 0 );
  GV.rplan = XLALCreateReverseREAL8FFTPlan( GV.seg_length, 0 );
  
  if (CLA.fakenoiseflag && CLA.whitespectrumflag){
    unsigned p;
    for ( p = 0 ; p < GV.Spec->data->length; p++ )
      /* FIXME:  shouldn't this be 2 * \Delta f */
      GV.Spec->data->data[p]=2/(1.0/ht->deltaT);
    GV.Spec->deltaF=1/(GV.seg_length*ht->deltaT);
  } else{
    int segmentLength = GV.seg_length;
    int segmentStride = GV.seg_length/2;
    REAL8Window *window = XLALCreateHannREAL8Window( segmentLength );
    
    if(XLALREAL8AverageSpectrumMedianMean( GV.Spec, ht, segmentLength,
					   segmentStride, window, GV.fplan ))
      return 1;

    XLALDestroyREAL8Window( window );
  }

  return 0;
}

/*******************************************************************************/

int DownSample(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht){
  if(XLALResampleREAL8TimeSeries(ht, 1.0/CLA.samplerate)) return 1;
  return 0;
}

/*******************************************************************************/

REAL8TimeSeries *ReadData(struct CommandLineArgsTag CLA){
  unsigned p;
  REAL8TimeSeries *ht = NULL;

  if(CLA.fakenoiseflag) {
    /* create random data set */
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    FILE *devrandom;
    int errorcode;
    unsigned long seed;

    if(!(devrandom=fopen("/dev/urandom","r"))){
      XLALPrintError("Unable to open device /dev/urandom\n");
      return NULL;
    }
    errorcode=fread(&seed, sizeof(seed), 1, devrandom);
    if(errorcode!=1){
      XLALPrintError("Error reading /dev/urandom file!\n");
      return NULL;
    }
    fclose(devrandom);
    gsl_rng_set(rng, seed);

    /* hard-code sample rate of simulated noise to 16384 Hz */
    ht = XLALCreateREAL8TimeSeries("white noise", &CLA.GPSStart, 0.0, 1.0 / 16384, &lalDimensionlessUnit, XLALGPSDiff(&CLA.GPSEnd, &CLA.GPSStart) * 16384);
    for(p = 0; p < ht->data->length; p++)
      ht->data->data[p] = gsl_ran_gaussian(rng, 1.0);
    gsl_rng_free(rng);
  } else {
    FrCache *cache;
    FrStream *stream;
    LALTYPECODE series_type;

    /* create Frame cache, open frame stream and delete frame cache */
    cache = XLALFrImportCache(CLA.FrCacheFile);
    if(!cache)
      XLAL_ERROR_NULL(__func__, XLAL_EFUNC);
    stream = XLALFrCacheOpen(cache);
    XLALFrDestroyCache(cache);
    if(!stream)
      XLAL_ERROR_NULL(__func__, XLAL_EFUNC);

    /* turn on checking for missing data */
    XLALFrSetMode(stream, LAL_FR_VERBOSE_MODE);

    /* get the data type */
    series_type = XLALFrGetTimeSeriesType(CLA.ChannelName, stream);
    if((int) series_type < 0) {
      XLALFrClose(stream);
      XLAL_ERROR_NULL(__func__, XLAL_EFUNC);
    }

    /* read data */
    switch(series_type) {
    case LAL_S_TYPE_CODE: {
      /* read single-precision data */
      REAL4TimeSeries *ht_V = XLALFrReadREAL4TimeSeries(stream, CLA.ChannelName, &CLA.GPSStart, XLALGPSDiff(&CLA.GPSEnd, &CLA.GPSStart), 0);
      if(!ht_V) {
        XLALFrClose(stream);
        XLAL_ERROR_NULL(__func__, XLAL_EFUNC);
      }

      /* cast to double precision */
      ht = XLALCreateREAL8TimeSeries(ht_V->name, &ht_V->epoch, ht_V->f0, ht_V->deltaT, &ht_V->sampleUnits, ht_V->data->length);
      for(p = 0; p < ht_V->data->length; p++)
        ht->data->data[p] = ht_V->data->data[p];

      /* clean up */
      XLALDestroyREAL4TimeSeries(ht_V);
      break;
    }

    case LAL_D_TYPE_CODE:
      /* read double-precision data */
      ht = XLALFrReadREAL8TimeSeries(stream, CLA.ChannelName, &CLA.GPSStart, XLALGPSDiff(&CLA.GPSEnd, &CLA.GPSStart), 0);
      if(!ht) {
        XLALFrClose(stream);
        XLAL_ERROR_NULL(__func__, XLAL_EFUNC);
      }
      break;

    default:
      XLAL_ERROR_NULL(__func__, XLAL_EINVAL);
    }

    /* close */
    XLALFrClose(stream);

    /* Scale data to avoid single float precision problems */
    for (p=0; p<ht->data->length; p++)
      ht->data->data[p] *= SCALE;

    /* FIXME:  ARGH!!!  frame files cannot be trusted to provide units for
     * their contents! */
    ht->sampleUnits = lalStrainUnit;
  }

  return ht;
}

/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA, const ProcessTable *process, ProcessParamsTable **paramaddpoint){
  static char default_comment[] = "";
  INT4 errflg = 0;
  struct option long_options[] = {
    {"bank-freq-start",           required_argument,	NULL,	'L'},
    {"bank-lowest-hifreq-cutoff", required_argument,	NULL,	'H'},
    {"max-mismatch",              required_argument,	NULL,	'M'},
    {"threshold",                 required_argument,	NULL,	't'},
    {"frame-cache",               required_argument,	NULL,	'F'},
    {"channel",                   required_argument,	NULL,	'C'},
    {"output",                    required_argument,	NULL,	'o'},
    {"gps-end-time",              required_argument,	NULL,	'E'},
    {"gps-start-time",            required_argument,	NULL,	'S'},
    {"injection-file",            required_argument,	NULL,	'i'},
    {"short-segment-duration",    required_argument,	NULL,	'd'},
    {"settling-time",             required_argument,	NULL,	'T'},
    {"sample-rate",               required_argument,	NULL,	's'},
    {"trig-start-time",           required_argument,	NULL,	'g'},
    {"pad",                       required_argument,	NULL,	'p'},
    {"chi2par0",                  required_argument,	NULL,	'A'},
    {"chi2par1",                  required_argument,	NULL,	'B'},
    {"chi2par2",                  required_argument,	NULL,	'G'},
    {"cusp-search",               no_argument,	NULL,	'c'},
    {"kink-search",               no_argument,	NULL,	'k'},
    {"test-gaussian-data",        no_argument,	NULL,	'n'},
    {"test-white-spectrum",       no_argument,	NULL,	'w'},
    {"cluster-events",            required_argument,	NULL,	'l'},
    {"print-spectrum",            no_argument,	NULL,	'a'},
    {"print-fd-filter",           no_argument,	NULL,	'b'},
    {"print-snr",                 no_argument,	NULL,	'r'},
    {"print-td-filter",           no_argument,	NULL,	'x'},
    {"print-data",                no_argument,	NULL,	'y'},
    {"print-injection",           no_argument,	NULL,	'z'},
    {"user-tag",                  required_argument,	NULL,	'j'},
    {"help",                      no_argument,	NULL,	'h'},
    {0, 0, 0, 0}
  };
  char args[] = "hnckwabrxyzlj:f:L:M:D:H:t:F:C:E:S:i:v:d:T:s:g:o:p:A:B:G:";

  optarg = NULL;

  /* Initialize default values */
  CLA->fbankstart=0.0;
  CLA->fbankhighfcutofflow=0.0;
  CLA->fmismatchmax=0.05;
  CLA->FrCacheFile=NULL;
  CLA->InjectionFile=NULL;
  CLA->ChannelName=NULL;
  CLA->outputFileName=NULL;
  XLALINT8NSToGPS(&CLA->GPSStart, 0);
  XLALINT8NSToGPS(&CLA->GPSEnd, 0);
  CLA->ShortSegDuration=0;
  CLA->TruncSecs=0;
  CLA->power=0.0;
  CLA->threshold=0.0;
  CLA->fakenoiseflag=0;
  CLA->whitespectrumflag=0;
  CLA->samplerate=4096.0;
  XLALINT8NSToGPS(&CLA->trigstarttime, 0);
  CLA->cluster=0.0;
  CLA->pad=0;
  CLA->printfilterflag=0;
  CLA->printspectrumflag=0;
  CLA->printsnrflag=0;
  CLA->printfirflag=0;
  CLA->printdataflag=0;
  CLA->printinjectionflag=0;
  CLA->comment=default_comment;

  /* initialise chi2cut */
  memset(CLA->chi2cut, 0, sizeof(CLA->chi2cut));
  CLA->chi2cut[0]=CLA->chi2cut[1]=CLA->chi2cut[2]=-9999.1;

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
    case 's':
      /* resample to this sample rate */
      CLA->samplerate=atof(optarg);
      ADD_PROCESS_PARAM(process, "float");
      break;
    case 'H':
      /* lowest high frequency cutoff */
      CLA->fbankhighfcutofflow=atof(optarg);
      ADD_PROCESS_PARAM(process, "float");
      break;
    case 'M':
      /* Maximal mismatch */
      CLA->fmismatchmax=atof(optarg);
      ADD_PROCESS_PARAM(process, "float");
      break;
    case 'L':
      /* low frequency cutoff */
      CLA->fbankstart=atof(optarg);
      ADD_PROCESS_PARAM(process, "float");
      break;
    case 't':
      /* low frequency cutoff */
      CLA->threshold=atof(optarg);
      ADD_PROCESS_PARAM(process, "float");
      break;
    case 'F':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'C':
      /* name channel */
      CLA->ChannelName=optarg;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'i':
      /* name of xml injection file */
      CLA->InjectionFile=optarg;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'o':
      /* name of xml injection file */
      CLA->outputFileName=optarg;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'S':
      /* GPS start time of search */
      if(XLALStrToGPS(&CLA->GPSStart, optarg, NULL)) {
        fprintf(stderr,"range error parsing \"%s\"", optarg);
        return 1;
      }
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'E':
      /* GPS end time time of search */
      if(XLALStrToGPS(&CLA->GPSEnd, optarg, NULL)) {
        fprintf(stderr,"range error parsing \"%s\"", optarg);
        return 1;
      }
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'd':
       /* Number of segment to break-up search into */
      CLA->ShortSegDuration=atoi(optarg);
      ADD_PROCESS_PARAM(process, "int");
      break;
    case 'T':
      /* Half the number of seconds that are trown out at the start and at the end of a short chunk */
      CLA->TruncSecs=atof(optarg);
      ADD_PROCESS_PARAM(process, "int");
      break;
    case 'g':
      /* start time of allowed triggers */
      if(XLALStrToGPS(&CLA->trigstarttime, optarg, NULL)) {
        fprintf(stderr,"range error parsing \"%s\"", optarg);
        return 1;
      }
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'p':
      /* start time of allowed triggers */
      CLA->pad=atoi(optarg);
      ADD_PROCESS_PARAM(process, "int");
      break;
    case 'A':
      /* chi2 cut parameter 0 */
      CLA->chi2cut[0]=atof(optarg);
      ADD_PROCESS_PARAM(process, "float");
      break;
    case 'B':
      /* chi2 cut parameter 1 */
      CLA->chi2cut[1]=atof(optarg);
      ADD_PROCESS_PARAM(process, "float");
      break;
    case 'G':
      /* chi2 cut parameter 2 */
      CLA->chi2cut[2]=atof(optarg);
      ADD_PROCESS_PARAM(process, "float");
      break;
    case 'c':
      /* cusp power law */
      CLA->power=-4.0/3.0;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'k':
      /* kink power law */
      CLA->power=-5.0/3.0;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'n':
      /* fake gaussian noise flag */
      CLA->fakenoiseflag=1;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'w':
      /* fake gaussian noise flag */
      CLA->whitespectrumflag=1;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'l':
      /* fake gaussian noise flag */
      CLA->cluster=atof(optarg);
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'a':
      /* fake gaussian noise flag */
      CLA->printspectrumflag=1;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'b':
      /* fake gaussian noise flag */
      CLA->printfilterflag=1;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'j':
      /* --user-tag */
      CLA->comment = optarg;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'r':
      /* fake gaussian noise flag */
      CLA->printsnrflag=1;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'x':
      /* fake gaussian noise flag */
      CLA->printfirflag=1;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'y':
      /* fake gaussian noise flag */
      CLA->printdataflag=1;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'z':
      /* fake gaussian noise flag */
      CLA->printinjectionflag=1;
      ADD_PROCESS_PARAM(process, "string");
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"All arguments are required except -n, -h, -w, -g, -o, -x, -y, -z, -b, -r -a, -l, -p, -T  and -i. One of -k or -c must be specified. They are:\n");
      fprintf(stdout,"\t--sample-rate (-s)\t\tFLOAT\t Desired sample rate (Hz).\n");
      fprintf(stdout,"\t--bank-lowest-hifreq-cutoff (-H)\tFLOAT\t Template bank lowest high frequency cut-off.\n");
      fprintf(stdout,"\t--max-mismatch (-M)\tFLOAT\t Maximal mismatch allowed from 1 template to the next.\n");
      fprintf(stdout,"\t--bank-freq-start (-L)\tFLOAT\t Template bank low frequency cut-off.\n");
      fprintf(stdout,"\t--threshold (-t)\t\tFLOAT\t SNR threshold.\n");
      fprintf(stdout,"\t--frame-cache (-F)\t\tSTRING\t Name of frame cache file.\n");
      fprintf(stdout,"\t--channel (-C)\t\tSTRING\t Name of channel.\n");
      fprintf(stdout,"\t--injection-file (-i)\t\tSTRING\t Name of xml injection file.\n");
      fprintf(stdout,"\t--output (-o)\t\tSTRING\t Name of xml output file.\n");
      fprintf(stdout,"\t--gps-start-time (-S)\t\tINTEGER\t GPS start time.\n");
      fprintf(stdout,"\t--gps-end-time (-E)\t\tINTEGER\t GPS end time.\n");
      fprintf(stdout,"\t--settling-time (-T)\t\tINTEGER\t Number of seconds to truncate filter.\n");
      fprintf(stdout,"\t--trig-start-time (-g)\t\tINTEGER\t GPS start time of triggers to consider.\n");
      fprintf(stdout,"\t--pad (-p)\t\tINTEGER\t Pad the data with these many seconds at beginning and end.\n");
      fprintf(stdout,"\t--chi2par0 (-A)\t\tFLOAT\t parameter[0] for the chi2 selection.\n");
      fprintf(stdout,"\t--chi2par1 (-B)\t\tFLOAT\t parameter[1] for the chi2 selection.\n");
      fprintf(stdout,"\t--chi2par2 (-G)\t\tFLOAT\t parameter[2] for the chi2 selection.\n");
      fprintf(stdout,"\t--short-segment-duration (-d)\t\tINTEGER\t Duration of short segments. They will overlap by 50%s. \n","%");
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
      fprintf(stdout,"eg %s  --sample-rate 4096 --bank-freq-start 30 --bank-lowest-hifreq-cutoff 200 --settling-time 0.1 --short-segment-duration 4 --cusp-search --cluster-events 0.1 --pad 4 --threshold 4 --output ladida.xml --frame-cache cache/H-H1_RDS_C01_LX-795169179-795171015.cache --channel H1:LSC-STRAIN --gps-start-time 795170318 --gps-end-time 795170396\n", argv[0]);
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
  if(!(CLA->ChannelName[0] == 'V' || CLA->ChannelName[0] == 'H' || CLA->ChannelName[0] == 'L'))
    {
      fprintf(stderr,"The channel name is  not well specified\n");
      fprintf(stderr,"It should start with H1, H2, L1 or V1\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(XLALGPSToINT8NS(&CLA->GPSStart) == 0)
    {
      fprintf(stderr,"No GPS start time specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->GPSStart.gpsNanoSeconds)
    {
      fprintf(stderr,"Only integer values allowed for --gps-start-time.\n");
      return 1;
    }
  if(XLALGPSToINT8NS(&CLA->GPSEnd) == 0)
    {
      fprintf(stderr,"No GPS end time specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->GPSEnd.gpsNanoSeconds)
    {
      fprintf(stderr,"Only integer values allowed for --gps-end-time.\n");
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
    int big_seg_length=XLALGPSDiff(&CLA->GPSEnd, &CLA->GPSStart)-2*CLA->pad;

    REAL4 x=((float)big_seg_length/(float)CLA->ShortSegDuration)-0.5;

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
      fprintf(stderr,"Short segment length t=%d is too small to accomodate truncation time requested.\n", CLA->ShortSegDuration);
	fprintf(stderr,"Need short segment t(=%d) to be >= 4 x Truncation length (%f).\n",CLA->ShortSegDuration,CLA->TruncSecs);
	return 1;
    }    
  }

  return errflg;
}

/*******************************************************************************/

int FreeMem(StringTemplate *strtemplate, int NTemplates){
  int m;
  
  XLALDestroyREAL8FrequencySeries(GV.Spec);
  
  for (m=0; m < NTemplates; m++)
    XLALDestroyREAL8FrequencySeries(strtemplate[m].StringFilter);
    
  XLALDestroyREAL8FFTPlan( GV.fplan );
  XLALDestroyREAL8FFTPlan( GV.rplan );
  
  LALCheckMemoryLeaks();
  
  return 0;
}

/*******************************************************************************/
