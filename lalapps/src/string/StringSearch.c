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

#include <complex.h>
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
#include <gsl/gsl_roots.h>

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

#define MAXTEMPLATES 50


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
  char *TemplateFile;         /* File with the list of fcutoff for a static template bank */
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
};

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

/* time window for trigger clustering */
/* FIXME:  global variables = BAD BAD BAD! (my fault -- Kipp) */
static double cluster_window;	/* seconds */

/***************************************************************************/

/* FUNCTION PROTOTYPES */

/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA, const ProcessTable *process, ProcessParamsTable **paramaddpoint);

/* Reads the template bank file */
int ReadTemplateFile(struct CommandLineArgsTag CLA, int *NTemplates_fix, REAL8 *fcutoff_fix);

/* Reads raw data (or puts in fake gaussian noise with a sigma=10^-20) */
REAL8TimeSeries *ReadData(struct CommandLineArgsTag CLA);

/* Adds injections if an xml injection file is given */
int AddInjections(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht);

/* DownSamples data */
int DownSample(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht);

/* Computes the average spectrum  */
REAL8FrequencySeries *AvgSpectrum(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, unsigned seg_length, REAL8FFTPlan *fplan);

/* Creates the template bank based on the spectrum  */
int CreateTemplateBank(struct CommandLineArgsTag CLA, unsigned seg_length, REAL8FrequencySeries *Spec, StringTemplate *strtemplate, int *NTemplates, REAL8 *fcutoff_fix, int NTemplates_fix, REAL8FFTPlan *rplan);

/* Creates the frequency domain string cusp or kink filters  */
int CreateStringFilters(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, unsigned seg_length, REAL8FrequencySeries *Spec, StringTemplate *strtemplate, int NTemplates, REAL8FFTPlan *fplan, REAL8FFTPlan *rplan);

/* Filters the data through the template banks  */
int FindStringBurst(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, unsigned seg_length, const StringTemplate *strtemplate, int NTemplates, REAL8FFTPlan *fplan, REAL8FFTPlan *rplan, SnglBurst **head);

/* Finds events above SNR threshold specified  */
int FindEvents(struct CommandLineArgsTag CLA, const StringTemplate *strtemplate,
               const REAL8TimeSeries *vector, SnglBurst **head);

/* Writes out the xml file with the events it found  */
int OutputEvents(const struct CommandLineArgsTag *CLA, ProcessTable *proctable, ProcessParamsTable *procparamtable, SnglBurst *events);

/* Frees the memory */
int FreeMem(StringTemplate *strtemplate, int NTemplates);

/* Clustering comparison function */
static int XLALCompareStringBurstByTime(const SnglBurst * const *, const SnglBurst * const *);

/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  struct CommandLineArgsTag CommandLineArgs;
  unsigned seg_length;
  StringTemplate strtemplate[MAXTEMPLATES];
  int NTemplates;
  int NTemplates_fix; /* number of template given by the template bank file */
  REAL8 fcutoff_fix[MAXTEMPLATES]; /* high frequency cutoffs given by the template bank file */
  SnglBurst *events=NULL;
  MetadataTable  process;
  MetadataTable  procparams;
  REAL8TimeSeries *ht;        /* strain data */
  REAL8FFTPlan *fplan;        /* fft plan */
  REAL8FFTPlan *rplan;        /* fft plan */
  REAL8FrequencySeries *Spec; /* average spectrum */


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
  /* set the trigger cluster window global variable */
  cluster_window = CommandLineArgs.cluster;

  /****** ReadTemplatefile ******/
  if (CommandLineArgs.TemplateFile != NULL) {
    XLALPrintInfo("ReadTemplateFile()\n");
    if (ReadTemplateFile(CommandLineArgs,&NTemplates_fix,fcutoff_fix)) return 3;
  }

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
  ht = XLALResizeREAL8TimeSeries(ht, (int)round(CommandLineArgs.pad/ht->deltaT),
				 ht->data->length-2*(int)round(CommandLineArgs.pad/ht->deltaT));

  /****** FFT plans ******/
  seg_length = round(CommandLineArgs.ShortSegDuration / ht->deltaT);
  fplan = XLALCreateForwardREAL8FFTPlan( seg_length, 1 );
  rplan = XLALCreateReverseREAL8FFTPlan( seg_length, 1 );
  if(!fplan || !rplan) {
    XLALDestroyREAL8FFTPlan(fplan);
    XLALDestroyREAL8FFTPlan(rplan);
    return 9;
  }

  /****** AvgSpectrum ******/
  XLALPrintInfo("AvgSpectrum()\n");
  Spec = AvgSpectrum(CommandLineArgs, ht, seg_length, fplan);
  if (!Spec) return 9;
  if (CommandLineArgs.printspectrumflag) LALDPrintFrequencySeries( Spec, "Spectrum.txt" );

  /****** CreateTemplateBank ******/
  XLALPrintInfo("CreateTemplateBank()\n");
  if (CreateTemplateBank(CommandLineArgs, seg_length, Spec, strtemplate, &NTemplates, fcutoff_fix, NTemplates_fix, rplan)) return 10;

  /****** CreateStringFilters ******/
  XLALPrintInfo("CreateStringFilters()\n");
  if (CreateStringFilters(CommandLineArgs, ht, seg_length, Spec, strtemplate, NTemplates, fplan, rplan)) return 11;
  XLALDestroyREAL8FrequencySeries(Spec);
  Spec = NULL;

  /****** FindStringBurst ******/
  XLALPrintInfo("FindStringBurst()\n");
  if (FindStringBurst(CommandLineArgs, ht, seg_length, strtemplate, NTemplates, fplan, rplan, &events)) return 12;
  if(!XLALSortSnglBurst(&events, XLALCompareSnglBurstByExactPeakTime)) return 12;
  XLALDestroyREAL8TimeSeries(ht);
  XLALDestroyREAL8FFTPlan(fplan);
  XLALDestroyREAL8FFTPlan(rplan);

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


static int XLALCompareStringBurstByTime(
  const SnglBurst * const *a,
  const SnglBurst * const *b
)

{
  double delta_t = XLALGPSDiff(&(*a)->peak_time, &(*b)->peak_time);

  if(delta_t > cluster_window)
    return(1);
  if(delta_t < -cluster_window)
    return(-1);
  return(0);
}

/*******************************************************************************/

int AddInjections(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht){
  TimeSlide *time_slide_table_head;
  SimBurst *sim_burst_table_head;
  COMPLEX16FrequencySeries *response = NULL;

  /* Get info from injection file */
  time_slide_table_head = XLALTimeSlideTableFromLIGOLw(CLA.InjectionFile);
  sim_burst_table_head = XLALSimBurstTableFromLIGOLw(CLA.InjectionFile, NULL, NULL);
  if(!time_slide_table_head || !sim_burst_table_head)
    return 1;

  /* Construct response function for null stream */
  if(0) {	/* FIXME:  put proper test for null stream here */
    /* reduce injection amplitude by 10x for null stream.  injection will
     * be Fourier transformed, and the transform divided by this function
     * bin-by-bin, rounding to the closest available bin */
    response = XLALCreateCOMPLEX16FrequencySeries("", &ht->epoch, 0.0, 1.0, &lalDimensionlessUnit, 1);
    if(!response)
      return 1;
    response->data->data[0] = 10;
  }

  /* Inject the signals into ht */
  if(XLALBurstInjectSignals(ht, sim_burst_table_head, time_slide_table_head, response)) return 1;
  XLALDestroyCOMPLEX16FrequencySeries(response);

  /* free the injection table */
  XLALDestroyTimeSlideTable(time_slide_table_head);
  XLALDestroySimBurstTable(sim_burst_table_head);

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
    snprintf(outfilename, sizeof(outfilename)-1, "%s-STRINGSEARCH-%d-%d.xml", ifo, CLA->GPSStart.gpsSeconds, CLA->GPSEnd.gpsSeconds - CLA->GPSEnd.gpsSeconds);
    outfilename[sizeof(outfilename)-1] = '\0';
    xml = XLALOpenLIGOLwXMLFile(outfilename);
  } else
    xml = XLALOpenLIGOLwXMLFile(CLA->outputFileName);

  /* finish populating process table */
  snprintf(proctable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo);
  XLALGPSTimeNow(&(proctable->end_time));

  /* write process table */
  if(XLALWriteLIGOLwXMLProcessTable(xml, proctable)) return -1;

  /* write process params table */
  if(XLALWriteLIGOLwXMLProcessParamsTable(xml, procparamtable)) return -1;

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

  /* write search_summary table */
  if(XLALWriteLIGOLwXMLSearchSummaryTable(xml, searchsumm.searchSummaryTable)) return -1;
  XLALDestroySearchSummaryTable(searchsumm.searchSummaryTable);

  /* burst table */
  if(XLALWriteLIGOLwXMLSnglBurstTable(xml, events)) return -1;

  XLALCloseLIGOLwXMLFile(xml);

  return 0;
}

/*******************************************************************************/

int FindEvents(struct CommandLineArgsTag CLA, const StringTemplate *strtemplate, const REAL8TimeSeries *vector, SnglBurst **head){
  unsigned p;
  INT4 pmax, pend, pstart;

  /* print the snr to stdout */
  if (CLA.printsnrflag)
    for ( p = vector->data->length/4 ; p < 3*vector->data->length/4; p++ )
      fprintf(stdout,"%p %e\n", strtemplate, vector->data->data[p]);

  /* Now find event in the inner half */
  for ( p = vector->data->length/4 ; p < 3*vector->data->length/4; p++ ){
    REAL8 maximum_snr = 0.0;
    pmax=p;

    /* Do we have the start of a cluster? */
    if ( (fabs(vector->data->data[p]) > CLA.threshold) && (XLALGPSDiff(&vector->epoch, &CLA.trigstarttime) + p * vector->deltaT >= 0)){
      SnglBurst *new;
      REAL8 chi2, ndof;
      int pp;
      pend=p; pstart=p;

      /* Clustering in time: While we are above threshold, or within clustering time of the last point above threshold... */
      while( ((fabs(vector->data->data[p]) > CLA.threshold) || ((p-pend)* vector->deltaT < (float)(CLA.cluster)) )
	     && p<3*vector->data->length/4){

	/* This keeps track of the largest SNR point of the cluster */
	if(fabs(vector->data->data[p]) > maximum_snr){
	  maximum_snr=fabs(vector->data->data[p]);
	  pmax=p;
	}
	/* pend is the last point above threshold */
	if ( (fabs(vector->data->data[p]) > CLA.threshold))
	  pend =  p;

	p++;
      }

      /* compute \chi^{2} */
      chi2=0, ndof=0;
      for(pp=-strtemplate->chi2_index; pp<strtemplate->chi2_index; pp++){
        chi2 += (vector->data->data[pmax+pp]-vector->data->data[pmax]*strtemplate->auto_cor->data[vector->data->length/2+pp])*(vector->data->data[pmax+pp]-vector->data->data[pmax]*strtemplate->auto_cor->data[vector->data->length/2+pp]);
        ndof += (1-strtemplate->auto_cor->data[vector->data->length/2+pp]*strtemplate->auto_cor->data[vector->data->length/2+pp]);
      }

      /* Apply the \chi^{2} cut */
      if( CLA.chi2cut[0]    > -9999
	  && CLA.chi2cut[1] > -9999
	  && CLA.chi2cut[2] > -9999 )
	if(log10(chi2/ndof)>CLA.chi2cut[0]
	   && log10(chi2/ndof)> CLA.chi2cut[1]*log10(fabs(maximum_snr))+CLA.chi2cut[2]) continue;

      /* prepend a new event to the linked list */
      new = XLALCreateSnglBurst();
      if ( ! new )
        XLAL_ERROR(XLAL_EFUNC);
      new->next = *head;
      *head = new;

      /* Now copy stuff into event */
      strncpy( new->ifo, CLA.ChannelName, 2 );
      new->ifo[2] = 0;
      strncpy( new->search, "StringCusp", sizeof( new->search ) );
      strncpy( new->channel, CLA.ChannelName, sizeof( new->channel ) );

      /* compute start and peak time and duration, give 1 sample of fuzz on
       * both sides */
      new->start_time = new->peak_time = vector->epoch;
      XLALGPSAdd(&new->peak_time, pmax * vector->deltaT);
      XLALGPSAdd(&new->start_time, (pstart - 1) * vector->deltaT);
      new->duration = vector->deltaT * ( pend - pstart + 2 );

      new->central_freq = (strtemplate->f+CLA.fbankstart)/2.0;
      new->bandwidth    = strtemplate->f-CLA.fbankstart;
      new->snr          = maximum_snr;
      new->amplitude    = vector->data->data[pmax]/strtemplate->norm;
      new->chisq = chi2;
      new->chisq_dof = ndof;
    }
  }

  return 0;
}

/*******************************************************************************/

int FindStringBurst(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, unsigned seg_length, const StringTemplate *strtemplate, int NTemplates, REAL8FFTPlan *fplan, REAL8FFTPlan *rplan, SnglBurst **head){
  int i,m;
  unsigned p;
  COMPLEX16FrequencySeries *vtilde;

  /* create vector that will hold FFT of data;  metadata will be populated
   * by FFT function */
  vtilde = XLALCreateCOMPLEX16FrequencySeries( ht->name, &ht->epoch, ht->f0, 0.0, &lalDimensionlessUnit, seg_length / 2 + 1 );

  /* loop over templates  */
  for (m = 0; m < NTemplates; m++){
    /* loop over overlapping chunks */
    for(i=0; i < 2*(ht->data->length*ht->deltaT)/CLA.ShortSegDuration - 1 ;i++){
      /* extract overlapping chunk of data */
      REAL8TimeSeries *vector = XLALCutREAL8TimeSeries(ht, i * seg_length / 2, seg_length);

      /* FFT it */
      if(XLALREAL8TimeFreqFFT( vtilde, vector, fplan )) return 1;

      /* multiply FT of data and String Filter */
      for ( p = 0 ; p < vtilde->data->length; p++ )
        vtilde->data->data[p] *= strtemplate[m].StringFilter->data->data[p];

      /* reverse FFT it */
      if(XLALREAL8FreqTimeFFT( vector, vtilde, rplan )) return 1;
      vector->deltaT = ht->deltaT;	/* gets mucked up by round-off */

      /* normalise the result by template normalisation
	 factor of 2 is from match-filter definition */
      for ( p = 0 ; p < vector->data->length; p++ )
	vector->data->data[p] *= 2.0 / strtemplate[m].norm;

      /* find triggers */
      if(FindEvents(CLA, &strtemplate[m], vector, head)) return 1;

      /*
      if(m==10){
      	for ( int ii = 0 ; ii < vector->data->length; ii++ )
	  printf("%.3e\n",vector->data->data[ii]);
      }
      */
      /* free chunk */
      XLALDestroyREAL8TimeSeries( vector );
    }
  }

  XLALDestroyCOMPLEX16FrequencySeries( vtilde );

  return 0;
}


/*******************************************************************************/

int CreateStringFilters(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, unsigned seg_length, REAL8FrequencySeries *Spec, StringTemplate *strtemplate, int NTemplates, REAL8FFTPlan *fplan, REAL8FFTPlan *rplan){
  int m;
  unsigned p;
  COMPLEX16FrequencySeries *vtilde; /* frequency-domain vector workspace */
  REAL8TimeSeries *vector; /* time-domain vector workspace */

  vector = XLALCreateREAL8TimeSeries( ht->name, &ht->epoch, ht->f0, ht->deltaT, &ht->sampleUnits, seg_length );
  vtilde = XLALCreateCOMPLEX16FrequencySeries( ht->name, &ht->epoch, ht->f0, 1.0 / (vector->data->length * vector->deltaT), &lalDimensionlessUnit, vector->data->length / 2 + 1 );

  for (m = 0; m < NTemplates; m++){
    /* Initialize the filter */
    strtemplate[m].StringFilter = XLALCreateREAL8FrequencySeries(CLA.ChannelName, &CLA.GPSStart, 0, Spec->deltaF, &lalStrainUnit, Spec->data->length);

    /* populate vtilde with the template divided by the noise */
    for ( p = 0; p < vtilde->data->length; p++ ){
      vtilde->data->data[p] = sqrt(strtemplate[m].waveform_f->data[p]/Spec->data->data[p]);
    }

    /* reverse FFT vtilde into vector */
    if(XLALREAL8FreqTimeFFT( vector, vtilde, rplan )) return 1;
    /* this gets mucked up by round-off each time through the loop */
    vector->deltaT = ht->deltaT;

    /* perform the truncation; the truncation is CLA.TruncSecs/2 because
       we are dealing with the sqrt of the filter at the moment*/
    if(CLA.TruncSecs != 0.0)
      memset( vector->data->data + (int)round(CLA.TruncSecs/2/vector->deltaT), 0,
	      ( vector->data->length -  2 * (int)round(CLA.TruncSecs/2/vector->deltaT))
	      * sizeof( *vector->data->data ) );

    /* forward fft the truncated vector into vtilde */
    if(XLALREAL8TimeFreqFFT( vtilde, vector, fplan )) return 1;
    /* this gets mucked up by round-off each time through the loop */
    vtilde->deltaF = 1.0 / (vector->data->length * vector->deltaT);

    /* store the square magnitude in the filter */
    for ( p = 0 ; p < vtilde->data->length; p++ )
      strtemplate[m].StringFilter->data->data[p] = creal(vtilde->data->data[p]) * creal(vtilde->data->data[p]) + cimag(vtilde->data->data[p]) * cimag(vtilde->data->data[p]);

    /* set DC and Nyquist to 0 */
    strtemplate[m].StringFilter->data->data[0] = strtemplate[m].StringFilter->data->data[strtemplate[m].StringFilter->data->length-1] = 0;

    /* print out the frequency domain filter */
    if (CLA.printfilterflag){
      CHAR filterfilename[256];
      snprintf(filterfilename, sizeof(filterfilename)-1, "Filter-%d.txt", m);
      filterfilename[sizeof(filterfilename)-1] = '\0';
      LALDPrintFrequencySeries( strtemplate[m].StringFilter, filterfilename );
    }

    /* print out the time domain FIR filter */
    if (CLA.printfirflag){
      char filterfilename[256];

      strncpy(vector->name, "fir filter", LALNameLength);
      for ( p = 0 ; p < vtilde->data->length; p++ ) {
        vtilde->data->data[p] = creal(vtilde->data->data[p]) * creal(vtilde->data->data[p]) + cimag(vtilde->data->data[p]) * cimag(vtilde->data->data[p]);
      }
      vtilde->data->data[0] = vtilde->data->data[vtilde->data->length - 1] = 0.0;
      XLALREAL8FreqTimeFFT( vector, vtilde, rplan );

      snprintf(filterfilename, sizeof(filterfilename)-1, "FIRFilter-%d.txt", m);
      filterfilename[sizeof(filterfilename)-1] = '\0';
      LALDPrintTimeSeries( vector, filterfilename );
    }
  }

  XLALDestroyCOMPLEX16FrequencySeries( vtilde );
  XLALDestroyREAL8TimeSeries( vector );

  return 0;
}

/*******************************************************************************/

/* compute (t2|t2) and (t1|t2) */
static void compute_t2t2_and_t1t2(double power, const REAL8FrequencySeries *Spec, const REAL8Vector *integral, double last_templates_f_cut, double f_cut, double *t2t2, double *t1t2)
{
  unsigned i = round((last_templates_f_cut - Spec->f0) / Spec->deltaF);

  *t2t2 = *t1t2 = integral->data[i];

  for(i++; i < Spec->data->length; i++) {
    double f = Spec->f0 + i * Spec->deltaF;

    if(f < f_cut) {
      *t2t2 += 4 * pow(pow(f, power), 2) / Spec->data->data[i] * Spec->deltaF;
      *t1t2 += 4 * pow(pow(f, power), 2) * exp(1 - f / last_templates_f_cut) / Spec->data->data[i] * Spec->deltaF;
    } else {
      *t2t2 += 4 * pow(pow(f, power) * exp(1 - f / f_cut), 2) / Spec->data->data[i] * Spec->deltaF;
      *t1t2 += 4 * pow(pow(f, power), 2) * exp(1 - f / last_templates_f_cut) * exp(1 - f / f_cut) / Spec->data->data[i] * Spec->deltaF;
    }
  }
}

struct compute_epsilon_minus_desired_params {
  double desired_epsilon;
  double string_spectrum_power;
  const REAL8FrequencySeries *Spec;
  const REAL8Vector *integral;
  double last_templates_f_cut;
  double last_templates_norm;
};

static double compute_epsilon_minus_desired(double f_cut, void *params)
{
  struct compute_epsilon_minus_desired_params *p = params;
  double epsilon;
  double t1t1 = pow(p->last_templates_norm, 2);
  double t2t2, t1t2;

  compute_t2t2_and_t1t2(p->string_spectrum_power, p->Spec, p->integral, p->last_templates_f_cut, f_cut, &t2t2, &t1t2);

  /* "epsilon" is the mismatch between two templates.  in this case we've
   * computed it between a template with a cut-off frequency placed at
   * f_cut a template with a cut-off frequency placed at
   * last_templates_f_cut */
  epsilon = 1 - t1t2 / sqrt(t1t1 * t2t2);

  /* the "desired epsilon" is the template bank mismatch provided on the
   * command line, by writing a function that returns the difference
   * between the measured epsilon and the value requested by the user we
   * can use a root-solver to find the f_cut that gives the desired epsilon
   * with repsect to the previous template */
  return epsilon - p->desired_epsilon;
}

static double next_f_cut(double desired_epsilon, double string_spectrum_power, const REAL8FrequencySeries *Spec, const REAL8Vector *integral, double last_templates_f_cut, double last_templates_norm)
{
  struct compute_epsilon_minus_desired_params params = {
    .desired_epsilon = desired_epsilon,
    .string_spectrum_power = string_spectrum_power,
    .Spec = Spec,
    .integral = integral,
    .last_templates_f_cut = last_templates_f_cut,
    .last_templates_norm = last_templates_norm
  };
  gsl_function F = {
    .function = compute_epsilon_minus_desired,
    .params = &params
  };
  gsl_root_fsolver *solver;
  /* we seek an f_cut bracketed by these values */
  double flo = last_templates_f_cut;
  double fhi = Spec->f0 + (Spec->data->length - 1) * Spec->deltaF;

  /* if there isn't enough mismatch to place another template between the
   * previous one and fhi, return fhi.  note that we must ensure that the
   * last template has exactly this frequency to cause the template
   * construction loop to terminate */
  if(compute_epsilon_minus_desired(fhi, &params) <= 0)
    return fhi;

  /* iterate the bisection algorithm until fhi and flo are less than 1
   * frequency bin apart */
  solver = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
  gsl_root_fsolver_set(solver, &F, flo, fhi);
  do {
    gsl_root_fsolver_iterate(solver);
    flo = gsl_root_fsolver_x_lower(solver);
    fhi = gsl_root_fsolver_x_upper(solver);
  } while(fhi - flo >= Spec->deltaF);
  gsl_root_fsolver_free(solver);

  /* return the mid-point as f_cut */
  return (fhi + flo) / 2;
}

int CreateTemplateBank(struct CommandLineArgsTag CLA, unsigned seg_length, REAL8FrequencySeries *Spec, StringTemplate *strtemplate, int *NTemplates, REAL8 *fcutoff_fix, int NTemplates_fix, REAL8FFTPlan *rplan){
  REAL8 f_cut, t1t1, t2t2, t1t2, norm, slope0, slope1;
  int m, f_low_cutoff_index, extr_ctr;
  unsigned p;
  REAL8Vector *integral;
  REAL8Vector *vector; /* time-domain vector workspace */
  COMPLEX16Vector *vtilde; /* frequency-domain vector workspace */

  *NTemplates = 0;

  /* populate integral */
  integral = XLALCreateREAL8Vector(Spec->data->length);
  memset(integral->data, 0, integral->length * sizeof(*integral->data));
  for( p = round((CLA.fbankstart - Spec->f0) / Spec->deltaF) ; p < integral->length; p++ ) {
    integral->data[p] = 4 * pow(pow(Spec->f0 + p * Spec->deltaF, CLA.power), 2) / Spec->data->data[p] * Spec->deltaF;
    if(p > 0)
      integral->data[p] += integral->data[p - 1];
  }

  /* Use static template bank or...*/
  if(CLA.TemplateFile){
    compute_t2t2_and_t1t2(CLA.power, Spec, integral, CLA.fbankstart, fcutoff_fix[0], &t1t1, &t1t2);

    strtemplate[0].findex = round((fcutoff_fix[0] - Spec->f0) / Spec->deltaF);
    strtemplate[0].f = fcutoff_fix[0];
    strtemplate[0].mismatch = 0.0;
    strtemplate[0].norm = sqrt(t1t1);
    XLALPrintInfo("%% Templ. frequency      sigma      mismatch\n");
    XLALPrintInfo("%% %d      %1.3e    %1.3e    %1.3e\n",0,strtemplate[0].f,strtemplate[0].norm, strtemplate[0].mismatch);
    *NTemplates = NTemplates_fix;

    for (m = 1; m < NTemplates_fix; m++){
      compute_t2t2_and_t1t2(CLA.power, Spec, integral, fcutoff_fix[m-1], fcutoff_fix[m], &t2t2, &t1t2);

      strtemplate[m].findex = round((fcutoff_fix[m] - Spec->f0) / Spec->deltaF);
      strtemplate[m].f = fcutoff_fix[m];
      strtemplate[m].norm = sqrt(t2t2);
      strtemplate[m].mismatch = 1 - t1t2 / sqrt(t1t1 * t2t2);
      XLALPrintInfo("%% %d      %1.3e    %1.3e    %1.3e\n", m, strtemplate[m].f, strtemplate[m].norm, strtemplate[m].mismatch);
      t1t1 = t2t2;
    }
  }

  /* ... or compute an "adapted" bank */
  else{
    f_cut = CLA.fbankhighfcutofflow;

    /* compute (t1|t1) for fist template.  we can do this by re-using the
     * (t2|t2),(t1|t2) function with the correct inputs.  t1t2 result is
     * meaningless and not used */
    compute_t2t2_and_t1t2(CLA.power, Spec, integral, CLA.fbankstart, f_cut, &t1t1, &t1t2);

    strtemplate[0].findex = round((f_cut - Spec->f0) / Spec->deltaF);
    strtemplate[0].f = f_cut;
    strtemplate[0].mismatch = 0.0;
    strtemplate[0].norm = sqrt(t1t1);
    XLALPrintInfo("%% Templ. frequency      sigma      mismatch\n");
    XLALPrintInfo("%% %d      %1.3e    %1.3e    %1.3e\n",*NTemplates,strtemplate[0].f,strtemplate[0].norm, strtemplate[0].mismatch);
    *NTemplates = 1;

    /* find the next cutoffs given the maximal mismatch, until we hit the
     * highest frequency.  note that the algorithm will hit that frequency
     * bin by construction */
    while(strtemplate[*NTemplates - 1].findex < (int) Spec->data->length - 1) {
      f_cut = next_f_cut(CLA.fmismatchmax, CLA.power, Spec, integral, strtemplate[*NTemplates-1].f, strtemplate[*NTemplates-1].norm);

      compute_t2t2_and_t1t2(CLA.power, Spec, integral, strtemplate[*NTemplates-1].f, f_cut, &t2t2, &t1t2);

      strtemplate[*NTemplates].findex = round((f_cut - Spec->f0) / Spec->deltaF);
      strtemplate[*NTemplates].f = f_cut;
      strtemplate[*NTemplates].norm = sqrt(t2t2);
      strtemplate[*NTemplates].mismatch = 1 - t1t2 / sqrt(t1t1 * t2t2);
      XLALPrintInfo("%% %d      %1.3e    %1.3e    %1.3e\n", *NTemplates, strtemplate[*NTemplates].f, strtemplate[*NTemplates].norm, strtemplate[*NTemplates].mismatch);
      (*NTemplates)++;
      if(*NTemplates == MAXTEMPLATES){
	XLALPrintError("Too many templates for code... Exiting\n");
	return 1;
      }
      t1t1 = t2t2;
    }
  }

  XLALDestroyREAL8Vector( integral );

  /* Now, the point is to store the template waveform vector */
  vector = XLALCreateREAL8Vector( seg_length);
  vtilde = XLALCreateCOMPLEX16Vector( vector->length / 2 + 1 );
  f_low_cutoff_index = round(CLA.fbankstart/ Spec->deltaF);
  for (m = 0; m < *NTemplates; m++){
    /* create the space for the waveform vectors */
    strtemplate[m].waveform_f = XLALCreateCOMPLEX16Vector( vtilde->length );
    strtemplate[m].waveform_t = XLALCreateREAL8Vector( vector->length );
    strtemplate[m].auto_cor   = XLALCreateREAL8Vector( vector->length );

    /* set all frequencies below the low freq cutoff to zero */
    memset(strtemplate[m].waveform_f->data, 0, f_low_cutoff_index*sizeof(*strtemplate[m].waveform_f->data));

    /* populate the rest with the template waveform */
    for ( p = f_low_cutoff_index; p < strtemplate[m].waveform_f->length; p++ ){
      double f = Spec->f0 + p * Spec->deltaF;
      if(f<=strtemplate[m].f)
	strtemplate[m].waveform_f->data[p] = pow(f, CLA.power);
      else
	strtemplate[m].waveform_f->data[p] = pow(f, CLA.power)*exp(1-f/strtemplate[m].f);
    }

    /* set DC and Nyquist to zero */
    strtemplate[m].waveform_f->data[0] = strtemplate[m].waveform_f->data[strtemplate[m].waveform_f->length - 1] = 0.0;

    /* whiten and convolve the template with itself, store in vtilde.
     * template is assumed to be real-valued */
    for (p=0 ; p< vtilde->length; p++)
      vtilde->data[p] = pow(creal(strtemplate[m].waveform_f->data[p]), 2) / Spec->data->data[p];

    /* reverse FFT */
    if(XLALREAL8ReverseFFT(vector, strtemplate[m].waveform_f, rplan)) return 1;
    if(XLALREAL8ReverseFFT(strtemplate[m].auto_cor, vtilde, rplan)) return 1;

    /* The vector is reshuffled in the right order */
    for ( p = 0 ; p < vector->length / 2; p++ ){
      strtemplate[m].waveform_t->data[p] = vector->data[vector->length/2+p]*Spec->deltaF;
      strtemplate[m].waveform_t->data[vector->length/2+p] = vector->data[p]*Spec->deltaF;
    }

    /* Normalize the autocorrelation by the central value */
    norm=strtemplate[m].auto_cor->data[0];
    for ( p = 0 ; p < strtemplate[m].auto_cor->length; p++ )
      strtemplate[m].auto_cor->data[p] /= norm;

    /* The vector is reshuffled in the right order */
    memcpy(vector->data, strtemplate[m].auto_cor->data, vector->length * sizeof(*vector->data));
    for ( p = 0 ; p < vector->length/2; p++ ){
      strtemplate[m].auto_cor->data[p] = vector->data[vector->length/2+p];
      strtemplate[m].auto_cor->data[vector->length/2+p] = vector->data[p];
    }

    /* search for the index of the 3rd extremum */
    extr_ctr=0;
    strtemplate[m].chi2_index=0;
    for ( p = strtemplate[m].waveform_t->length/2+1; p< strtemplate[m].waveform_t->length-1; p++ ){
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
REAL8FrequencySeries *AvgSpectrum(struct CommandLineArgsTag CLA, REAL8TimeSeries *ht, unsigned seg_length, REAL8FFTPlan *fplan){
  REAL8FrequencySeries *Spec;
  Spec  = XLALCreateREAL8FrequencySeries(CLA.ChannelName, &CLA.GPSStart, 0, 0, &lalStrainUnit, seg_length / 2 + 1);
  if(!Spec)
    XLAL_ERROR_NULL(XLAL_EFUNC);

  if (CLA.fakenoiseflag && CLA.whitespectrumflag){
    unsigned p;
    for ( p = 0 ; p < Spec->data->length; p++ )
      /* FIXME:  shouldn't this be 2 * \Delta f */
      Spec->data->data[p]=2/(1.0/ht->deltaT);
    Spec->deltaF=1/(seg_length*ht->deltaT);
  } else{
    unsigned segmentStride = seg_length/2;
    REAL8Window *window = XLALCreateHannREAL8Window( seg_length );
    if(!window)
      XLAL_ERROR_NULL(XLAL_EFUNC);

    if(XLALREAL8AverageSpectrumMedianMean( Spec, ht, seg_length, segmentStride, window, fplan ))
      XLAL_ERROR_NULL(XLAL_EFUNC);

    XLALDestroyREAL8Window( window );
  }

  return Spec;
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
      XLAL_ERROR_NULL(XLAL_EFUNC);
    stream = XLALFrCacheOpen(cache);
    XLALFrDestroyCache(cache);
    if(!stream)
      XLAL_ERROR_NULL(XLAL_EFUNC);

    /* turn on checking for missing data */
    XLALFrSetMode(stream, LAL_FR_VERBOSE_MODE);

    /* get the data type */
    series_type = XLALFrGetTimeSeriesType(CLA.ChannelName, stream);
    if((int) series_type < 0) {
      XLALFrClose(stream);
      XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    /* read data */
    switch(series_type) {
    case LAL_S_TYPE_CODE: {
      /* read single-precision data */
      REAL4TimeSeries *ht_V = XLALFrReadREAL4TimeSeries(stream, CLA.ChannelName, &CLA.GPSStart, XLALGPSDiff(&CLA.GPSEnd, &CLA.GPSStart), 0);
      if(!ht_V) {
        XLALFrClose(stream);
        XLAL_ERROR_NULL(XLAL_EFUNC);
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
        XLAL_ERROR_NULL(XLAL_EFUNC);
      }
      break;

    default:
      XLALFrClose(stream);
      XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* close */
    XLALFrClose(stream);

    /* FIXME:  ARGH!!!  frame files cannot be trusted to provide units for
     * their contents! */
    ht->sampleUnits = lalStrainUnit;
  }

  return ht;
}

/*******************************************************************************/

int ReadTemplateFile(struct CommandLineArgsTag CLA, int *NTemplates_fix, REAL8 *fcutoff_fix){

  SnglBurst *templates=NULL, *templates_root=NULL;
  int i;

  /* Initialize */
  *NTemplates_fix=0;

  /* Get templates from burst table */
  templates_root = XLALSnglBurstTableFromLIGOLw(CLA.TemplateFile);

  for(templates = templates_root; templates != NULL; templates = templates->next){
    fcutoff_fix[*NTemplates_fix]=templates->central_freq + (templates->bandwidth)/2;
    *NTemplates_fix=*NTemplates_fix+1;
    if(*NTemplates_fix==MAXTEMPLATES){
      XLALPrintError("Too many templates (> %d)\n",MAXTEMPLATES);
      return 1;
    }
  }

  /* Check that the bank has at least one template */
  if(*NTemplates_fix<=0){
    XLALPrintError("Empty template bank\n");
    return 1;
  }

  /* Check that the frequencies are well ordered */
  for(i=0; i<*NTemplates_fix-1; i++){
    if(fcutoff_fix[i]>fcutoff_fix[i+1]){
      XLALPrintError("The templates frequencies are not sorted by frequencies\n");
      return 1;
    }
  }

  /* check that  the highest frequency is below the Nyquist frequency */
  if(fcutoff_fix[*NTemplates_fix-1]>CLA.samplerate/2){
    XLALPrintError("The templates frequencies go beyond the Nyquist frequency\n");
    return 1;
  }

  return 0;
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
    {"template-bank",             required_argument,	NULL,	'K'},
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
  char args[] = "hnckwabrxyzlj:f:L:M:D:H:t:F:C:E:S:i:v:d:T:s:g:o:p:A:B:G:K:";

  optarg = NULL;

  /* Initialize default values */
  CLA->fbankstart=0.0;
  CLA->fbankhighfcutofflow=0.0;
  CLA->fmismatchmax=0.05;
  CLA->FrCacheFile=NULL;
  CLA->InjectionFile=NULL;
  CLA->TemplateFile=NULL;
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
    case 'K':
      /* name of txt template file */
      CLA->TemplateFile=optarg;
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
      fprintf(stdout,"\t--template-bank (-K)\t\tSTRING\t Name of txt template file.\n");
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
  if(! CLA->TemplateFile && CLA->fbankhighfcutofflow == 0.0)
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
      fprintf(stderr,"The channel name is not well specified\n");
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

  for (m=0; m < NTemplates; m++)
    XLALDestroyREAL8FrequencySeries(strtemplate[m].StringFilter);

  LALCheckMemoryLeaks();

  return 0;
}

/*******************************************************************************/
