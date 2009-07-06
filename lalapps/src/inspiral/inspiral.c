/*
*  Copyright (C) 2007 Andres C. Rodriguez, Stas Babak, Badri Krishnan, Sukanta Bose, Chad Hanna, Darren Woods, Alexander Dietz, Drew Keppel, Duncan Brown, Eirini Messaritaki, Gareth Jones, Jolien Creighton, Patrick Brady, Anand Sengupta, Stephen Fairhurst, Craig Robinson , Sean Seader, Thomas Cokelaer
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

/*-----------------------------------------------------------------------
 *
 * File Name: inspiral.c
 *
 * Author: Brown, D. A.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <math.h>

#include <FrameL.h>

#include <lalapps.h>
#include <series.h>
#include <processtable.h>
#include <lalappsfrutils.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/FrameStream.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/Calibration.h>
#include <lal/FrameCalibration.h>
#include <lal/Window.h>
#include <lal/TimeFreqFFT.h>
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpBCVSpin.h>
#include <lal/FindChirpPTF.h>
#include <lal/FindChirpChisq.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALTrigScanCluster.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>

#include "inspiral.h"

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "inspiral"

/* define the parameters for a 1.4,1.4 sloar mass standard candle with snr 8 */
#define CANDLE_MASS1 1.4
#define CANDLE_MASS2 1.4
#define CANDLE_RHOSQ 64.0


#define ADD_SUMM_VALUE( sv_name, sv_comment, val, intval ) \
if ( this_summ_value ) \
{ \
  this_summ_value = this_summ_value->next = (SummValueTable *) \
  LALCalloc( 1, sizeof(SummValueTable) );\
} \
else \
{ \
  summvalue.summValueTable = this_summ_value = (SummValueTable *) \
  LALCalloc( 1, sizeof(SummValueTable) );\
} \
LALSnprintf( this_summ_value->program, LIGOMETA_PROGRAM_MAX, "%s", \
  PROGRAM_NAME );\
this_summ_value->version = 0;\
this_summ_value->start_time = searchsumm.searchSummaryTable->out_start_time;\
this_summ_value->end_time = searchsumm.searchSummaryTable->out_end_time;\
this_summ_value->value = (REAL4) val;\
this_summ_value->intvalue = (INT4) intval;\
LALSnprintf( this_summ_value->name, LIGOMETA_SUMMVALUE_NAME_MAX, "%s", \
    sv_name );\
LALSnprintf( this_summ_value->ifo, LIGOMETA_IFO_MAX, "%s", ifo );\
LALSnprintf( this_summ_value->comment, LIGOMETA_SUMMVALUE_COMM_MAX, \
    "%s", sv_comment );\

double rint(double x);
int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

#ifdef LALAPPS_CONDOR
extern int condor_compress_ckpt;
void init_image_with_file_name( char *ckpt_file_name );
void ckpt_and_exit( void );
#endif

typedef struct
tagFrameHNode
{
  FrameH *frHeader;
  struct tagFrameHNode *next;
}
FrameHNode;



/*
 *
 * variables that control program behaviour
 *
 */


/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* checkpointing */
INT4  dataCheckpoint = 0;               /* condor checkpoint after data */
CHAR  ckptPath[FILENAME_MAX];           /* input and ckpt file path     */
CHAR  outputPath[FILENAME_MAX];         /* output data file path        */

/* input data parameters */
INT8  gpsStartTimeNS   = 0;             /* input data GPS start time ns */
LIGOTimeGPS gpsStartTime;               /* input data GPS start time    */
INT8  gpsEndTimeNS     = 0;             /* input data GPS end time ns   */
LIGOTimeGPS gpsEndTime;                 /* input data GPS end time      */
INT4  padData = 0;                      /* saftety margin on input data */
CHAR  *fqChanName       = NULL;         /* name of data channel         */
INT4  globFrameData     = 0;            /* glob *.gwf to get frame data */
CHAR  *frInCacheName    = NULL;         /* cache file containing frames */
CHAR  *frInType         = NULL;         /* type of data frames          */
INT4  numPoints         = -1;           /* points in a segment          */
INT4  numSegments       = -1;           /* number of segments           */
INT4  ovrlap            = -1;           /* overlap between segments     */
CHAR  ifo[3];                           /* two character ifo code       */
CHAR *channelName = NULL;               /* channel string               */
UINT4 inputDataLength = 0;              /* number of points in input    */
REAL4 strainHighPassFreq = -1;          /* h(t) high pass frequency     */
INT4  strainHighPassOrder = -1;         /* h(t) high pass filter order  */
REAL4 strainHighPassAtten = -1;         /* h(t) high pass attenuation   */
enum { undefined, real_4, real_8 } calData = undefined; /* cal data type*/

/* data conditioning parameters */
LIGOTimeGPS slideData   = {0,0};        /* slide data for time shifting */
INT4   resampFiltType   = -1;           /* low pass filter used for res */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
INT4   highPass         = -1;           /* enable high pass on raw data */
REAL4  highPassFreq     = 0;            /* high pass frequency          */
INT4   highPassOrder    = -1;           /* order of the td iir filter   */
REAL4  highPassAtten    = -1;           /* attenuation of the td filter */

REAL4  fLow             = -1;           /* low frequency cutoff         */

/* which type of PSD to use */
enum
{
  specType_mean,
  specType_median,
  specType_gaussian,
  specType_LIGO,
  specType_AdvLIGO,
  specType_undefined
} specType = specType_undefined;

INT4   badMeanPsd       = 0;            /* use a mean with no overlap   */
INT4   invSpecTrunc     = -1;           /* length of inverse spec (s)   */
REAL4  dynRangeExponent = -1;           /* exponent of dynamic range    */
CHAR  *calCacheName     = NULL;         /* location of calibration data */
INT4   globCalData      = 0;            /* glob for calibration frames  */
INT4   pointCal         = 0;            /* don't average cal over chunk */
CHAR  *injectionFile    = NULL;         /* name of file containing injs */
CHAR **tdFollowUpFiles  = NULL;         /* name of file containing td f */
INT4   numTDFiles       = 0;            /* Number of files to follow up */
int    injectOverhead   = 0;            /* inject h+ into detector      */
REAL4  mmFast           = -1.0;         /* match for the --fast option  */
INT4   hardwareInjection= 0;            /* hardware injection flag      */
Approximant injApproximant;             /* injected waveform approximant*/

/* matched filter parameters */
CHAR *bankFileName      = NULL;         /* name of input template bank  */
INT4  startTemplate     = -1;           /* index of first template      */
INT4  stopTemplate      = -1;           /* index of last template       */
INT4  numChisqBins      = -1;           /* number of chisq bins         */
REAL4 chisqDelta        = -1;           /* set chisq delta param        */
REAL4 snrThresh         = -1;           /* signal to noise thresholds   */
REAL4 chisqThresh       = -1;           /* chisq veto thresholds        */
FindChirpClustering clusterMethod;      /* chosen clustering algorithm  */
REAL4 clusterWindow     = -1;           /* cluster over time window     */
Approximant approximant;                /* waveform approximant         */
CHAR *approximantName   = NULL;         /* waveform approximant name    */
LALPNOrder order;                       /* pN order of waveform         */
CHAR *orderName = NULL;                 /* pN order of the waveform     */
INT4 bcvConstraint      = 0;            /* constraint BCV filter        */
INT4 flagFilterInjOnly  = -1;           /* flag for filtering inj. only */

/* rsq veto params */
INT4 enableRsqVeto      = -1;           /* enable the r^2 veto          */
REAL4 rsqVetoWindow     = -1;           /* r^2 veto time window         */
REAL4 rsqVetoThresh     = -1;           /* r^2 veto threshold           */
INT4 doRsqVeto          = 0;            /* do the r^2 veto              */
REAL4 rsqVetoTimeThresh = -1;           /* r^2 veto time threshold      */
REAL4 rsqVetoMaxSNR     = -1;           /* r^2 veto maximum snr         */
REAL4 rsqVetoCoeff      = -1;           /* r^2 veto coefficient         */
REAL4 rsqVetoPow        = -1;           /* r^2 veto power               */

/* generic simulation parameters */
enum { unset, urandom, user } randSeedType = unset;    /* sim seed type */
INT4  randomSeed        = 0;            /* value of sim rand seed       */
REAL4 gaussVar          = 64.0;         /* variance of Gaussian noise   */
INT4  whiteGaussian     = 0;            /* make input data Gaussian     */
INT4  unitResponse      = 0;            /* set the response to unity    */
INT4  colorSpec         = 0;            /* set the spectrum for colored */
                                        /* Gaussian noise               */
INT4  coloredGaussian   = 0;            /* generate colored Gaussian    */
                                        /* noise                        */
/* template bank simulation params */
INT4  bankSim           = 0;            /* number of template bank sims */
CHAR *bankSimFileName   = NULL;         /* file contining sim_inspiral  */
FindChirpBankSimParams bankSimParams = { 0, 0, -1, -1, NULL, -1, NULL, NULL, -1 };
                                        /* template bank sim params     */

/* reverse chirp bank option */
INT4 reverseChirpBank      = 0;         /* enable the reverse chirp     */
                                        /* template bank option         */

/* taper/band pass template options */
INT4  bandPassTmplt           = 0;
InspiralApplyTaper taperTmplt = INSPIRAL_TAPER_NONE;

/* template bank veto options */
UINT4 subBankSize          = 0;         /* num templates in a subbank   */
UINT4 ccFlag = 0;
/* output parameters */
CHAR  *userTag          = NULL;         /* string the user can tag with */
CHAR  *ifoTag           = NULL;         /* string to tag parent IFOs    */
CHAR   fileName[FILENAME_MAX];          /* name of output files         */
INT4   outCompress      = 0;
INT4   maximizationInterval = 0;        /* Max over template in this    */
                                        /* maximizationInterval Nanosec */
                                        /* interval                     */
trigScanType trigScanMethod = trigScanNone;
                                        /* Switch for clustering        */
                                        /* triggers in template         */
                                        /* parameters and end time      */
REAL8  trigScanDeltaEndTime = -1.0;      /* Use this interval (msec)     */
                                        /* over trigger end time while  */
                                        /* using trigScanCluster        */
REAL8  trigScanMetricScalingFac = -1.0;
/* Use this scaling factor for the volume spanned by a trigger in the   */
/* parameter space. When set to x, the volume is taken to be that of the*/
/* ambiguity ellipsoid at a 'minimal match' of (1.0-x).                 */
                                        /* original bank entered at the */
                                        /* command line                 */
INT2  trigScanAppendStragglers = -1;    /* Switch to append cluster     */
                                        /* out-liers (stragglers)       */

INT8   trigStartTimeNS  = 0;            /* write triggers only after    */
INT8   trigEndTimeNS    = 0;            /* write triggers only before   */
INT8   outTimeNS        = 0;            /* search summ out time         */
int    enableOutput     = -1;           /* write out inspiral events    */
MetadataTableType outputMask = sngl_inspiral_table; /* default to all   */
int    writeRawData     = 0;            /* write the raw data to a file */
int    writeFilterData  = 0;            /* write post injection data    */
int    writeResponse    = 0;            /* write response function used */
int    writeSpectrum    = 0;            /* write computed psd to file   */
int    writeRhosq       = 0;            /* write rhosq time series      */
int    writeChisq       = 0;            /* write chisq time series      */
int    writeCData       = 0;            /* write complex time series c  */
int    writeTemplate    = 0;            /* write the template time series */

/* other command line args */
CHAR comment[LIGOMETA_COMMENT_MAX];     /* process param comment        */

int main( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  /* frame input data */
  FrCache      *frInCache = NULL;
  FrCache      *frGlobCache = NULL;
  FrCache      *calCache = NULL;
  FrStream     *frStream = NULL;
  FrChanIn      frChan;
  FrCacheSieve  sieve;
  const size_t  calGlobLen = FILENAME_MAX;
  CHAR         *calGlobPattern;

  /* frame output data */
  struct FrFile *frOutFile  = NULL;
  struct FrameH *outFrame   = NULL;
  FrameHNode *coherentFrames = NULL;
  FrameHNode *thisCoherentFrame = NULL;
  UINT4          nRhosqFr = 0;
  UINT4          nChisqFr = 0;
  REAL4TimeSeries 		templateTimeSeries;
  REAL4Vector                  *templateTimeSeriesVector = NULL;
  /* raw input data storage */
  REAL4TimeSeries               chan;
  REAL8TimeSeries               strainChan;
  REAL4FrequencySeries          spec;
  COMPLEX8FrequencySeries       resp;
  DataSegmentVector            *dataSegVec = NULL;
  COMPLEX8TimeSeries           *coherentInputData = NULL;

  /* structures for preconditioning */
  COMPLEX8FrequencySeries       injResp;
  COMPLEX8FrequencySeries      *injRespPtr;
  ResampleTSParams              resampleParams;
  LALWindowParams               wpars;
  AverageSpectrumParams         avgSpecParams;

  /* findchirp data structures */
  FindChirpInitParams          *fcInitParams   = NULL;
  FindChirpSegmentVector       *fcSegVec       = NULL;
  FindChirpDataParams          *fcDataParams   = NULL;
  FindChirpTmpltParams         *fcTmpltParams  = NULL;
  FindChirpFilterParams        *fcFilterParams = NULL;
  FindChirpFilterInput         *fcFilterInput  = NULL;
  FindChirpStandardCandle       candle;

  /* inspiral template structures */
  INT4                          numTmplts    = 0;
  InspiralTemplate             *bankHead     = NULL;
  InspiralTemplate             *bankCurrent  = NULL;

  /* template bank veto structures */
  UINT4                         maxSubBankSize = 0;
  UINT4                         subBankIndex = 0;
  FindChirpSubBank             *subBankHead = NULL;
  FindChirpSubBank             *subBankCurrent = NULL;
  FindChirpBankVetoData         bankVetoData;

  /* inspiral events */
  INT4                          numEvents   = 0;
  SnglInspiralTable            *event       = NULL;
  SnglInspiralTable            *eventList   = NULL;
  SnglInspiralTable            *tempTmplt   = NULL;
  MetadataTable                 savedEvents;

  /* output data */
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsumm;
  MetadataTable         searchsummvars;
  MetadataTable         siminspiral;
  MetadataTable         siminstparams;
  MetadataTable         summvalue;
  MetadataTable         filtertable;
  SearchSummvarsTable  *this_search_summvar = NULL;
  SummValueTable       *this_summ_value = NULL;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       results;

  /* counters and other variables */
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  UINT4 i, j, k;
  INT4  cDataForFrame = 0;
  CHAR  fname[FILENAME_MAX];
  CHAR  cdataStr[LALNameLength];
  REAL8 inputLengthNS;
  UINT4 numInputPoints;
  const REAL8 epsilon = 1.0e-8;
  UINT4 resampleChan = 0;
  REAL8 tsLength;
  INT8  durationNS = 0;
  CalibrationUpdateParams calfacts, inj_calfacts;
  REAL4 alpha = 0;
  REAL4 alphabeta = 0;
  REAL4 inj_alpha = 0;
  REAL4 inj_alphabeta = 0;
  CHAR tmpChName[LALNameLength];
  REAL8 inputDeltaT;
  REAL8 dynRange = 1.0;
  REAL8 trigTime = 0.0;
  REAL8 lowerBound = 0.0;
  REAL8 upperBound = 0.0;

  /* random number generator parameters */
  RandomParams *randParams = NULL;

  /* template bank simulation variables */
  UINT4 bankSimCutLowIndex = 0;
  INT4  bankSimCount = 0;
  REAL4 matchNorm = 0;
  SnglInspiralTable  *loudestEventHead = NULL;
  SnglInspiralTable  *thisLoudestEvent = NULL;
  SimInspiralTable   *thisSimInspiral = NULL;
  SimInstParamsTable *thisSimInstParams = NULL;
  SimInspiralTable   *bankSimHead = NULL;
  SimInspiralTable   *thisBankSim = NULL;
  FindChirpBankSimParams *bankSimParamsPtr = NULL;

  /* injection information */
  int                  numInjections = 0;
  SimInspiralTable    *injections = NULL;
  SimInspiralTable    *thisInj = NULL;

  /* Time domain follow-up events */
  SnglInspiralTable    *tdFollowUpEvents = NULL;
  SnglInspiralTable    *thisFollowUpEvent  = NULL;

  /* --fast option related variables */
  UINT4  *analyseThisTmplt = NULL;
  INT4    thisTemplateIndex = 0;
  UINT4   analyseTag;
  

  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  if (strcmp(CVS_REVISION,"$Revi" "sion$"))
    {
      LAL_CALL( populate_process_table( &status, proctable.processTable,
					PROGRAM_NAME, CVS_REVISION,
					CVS_SOURCE, CVS_DATE ), &status );
    }
  else
    {
      LAL_CALL( populate_process_table( &status, proctable.processTable,
					PROGRAM_NAME, lalappsGitCommitID,
					lalappsGitGitStatus,
					lalappsGitCommitDate ), &status );
    }
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *)
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* create the search summary and zero out the summvars table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );
  searchsummvars.searchSummvarsTable = NULL;

  /* create the filter table */
  filtertable.filterTable = (FilterTable *)
    calloc( 1, sizeof(FilterTable) );

  /* zero out the checkpoint and output paths */
  memset( ckptPath, 0, FILENAME_MAX * sizeof(CHAR) );
  memset( outputPath, 0, FILENAME_MAX * sizeof(CHAR) );

  /* call the argument parse and check function */
  arg_parse_check( argc, argv, procparams );

  /* wind to the end of the process params table */
  for ( this_proc_param = procparams.processParamsTable; this_proc_param->next;
      this_proc_param = this_proc_param->next );

  /* can use LALMalloc() and LALCalloc() from here onwards */

  /* populate the filter table */
  LALSnprintf( filtertable.filterTable->program, LIGOMETA_PROGRAM_MAX, "%s",
      PROGRAM_NAME );
  filtertable.filterTable->start_time = gpsStartTime.gpsSeconds;
  LALSnprintf( filtertable.filterTable->filter_name, LIGOMETA_COMMENT_MAX,
      "%s%s", approximantName, orderName );

  /* fill the comment, if a user has specified on, or leave it blank */
  if ( ! *comment )
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
    LALSnprintf( filtertable.filterTable->comment, LIGOMETA_COMMENT_MAX, " " );
  }
  else
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
    LALSnprintf( filtertable.filterTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  /* put the name of the search in the search_summary comment */
  LALSnprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
       "%s%s", approximantName, orderName );

  /* set the name of the output file */
  if ( userTag && ifoTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INSPIRAL_%s_%s-%d-%d", ifo,
        ifoTag, userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if ( userTag && !ifoTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INSPIRAL_%s-%d-%d", ifo,
        userTag,  gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if ( !userTag && ifoTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INSPIRAL_%s-%d-%d", ifo,
        ifoTag,  gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INSPIRAL-%d-%d", ifo,
        gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }

  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;

  /* fill the ifos field of the search summary table */
  LALSnprintf( searchsumm.searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, ifo );

  /* make sure all the output table pointers are null */
  savedEvents.snglInspiralTable = NULL;
  siminspiral.simInspiralTable = NULL;
  siminstparams.simInstParamsTable = NULL;

  /* create the standard candle and database table */
  memset( &candle, 0, sizeof(FindChirpStandardCandle) );
  strncpy( candle.ifo, ifo, 2 * sizeof(CHAR) );
  candle.tmplt.mass1 = CANDLE_MASS1;
  candle.tmplt.mass2 = CANDLE_MASS2;
  candle.rhosq       = CANDLE_RHOSQ;
  candle.tmplt.totalMass = candle.tmplt.mass1 + candle.tmplt.mass2;
  candle.tmplt.mu = candle.tmplt.mass1 * candle.tmplt.mass2 /
    candle.tmplt.totalMass;
  candle.tmplt.eta = candle.tmplt.mu / candle.tmplt.totalMass;

  /* create the dynamic range exponent */
  dynRange = pow( 2.0, dynRangeExponent );

  /*
   *
   * create a template bank
   *
   */


  /* read in the template bank from a ligo lw xml file */
  numTmplts = InspiralTmpltBankFromLIGOLw( &bankHead, bankFileName,
      startTemplate, stopTemplate );
  if ( numTmplts < 0 )
  {
    fprintf( stderr, "error: unable to read templates from %s\n",
        bankFileName );
    exit( 1 );
  }
  else if ( numTmplts == 0 )
  {
    /* if there are no tmplts, store the time we would have analyzed and exit */
    fprintf( stdout, "no templates found in template bank file: %s\n"
        "exiting without searching for events.\n", bankFileName );

    searchsumm.searchSummaryTable->out_start_time.gpsSeconds =
      gpsStartTime.gpsSeconds + (numPoints / (4 * sampleRate));

    outTimeNS = XLALGPSToINT8NS( &(searchsumm.searchSummaryTable->out_start_time) );

    if ( ! bankSim && ( trigStartTimeNS && (trigStartTimeNS > outTimeNS) ) )
    {
      XLALINT8NSToGPS( &(searchsumm.searchSummaryTable->out_start_time),
            trigStartTimeNS );
    }

    searchsumm.searchSummaryTable->out_end_time.gpsSeconds =
      gpsEndTime.gpsSeconds - (numPoints / (4 * sampleRate));

    outTimeNS = XLALGPSToINT8NS( &(searchsumm.searchSummaryTable->out_end_time) );

    if ( ! bankSim && ( trigEndTimeNS && (trigEndTimeNS < outTimeNS) ) )
    {
      XLALINT8NSToGPS( &(searchsumm.searchSummaryTable->out_end_time),
            trigEndTimeNS );
    }
  }

  if ( vrbflg ) fprintf( stdout, "parsed %d templates from %s\n",
      numTmplts, bankFileName );


  /*
   *
   * read in the input data channel
   *
   */


  /* set the time series parameters of the input data and resample params */
  memset( &resampleParams, 0, sizeof(ResampleTSParams) );
  resampleParams.deltaT = 1.0 / (REAL8) sampleRate;

  /* set the params of the input data time series */
  memset( &chan, 0, sizeof(REAL4TimeSeries) );
  memset( &strainChan, 0, sizeof(REAL8TimeSeries) );
  chan.epoch = gpsStartTime;
  chan.epoch.gpsSeconds -= padData;   /* subtract pad seconds from start */
  /* subtract slide from start */
  chan.epoch.gpsSeconds -= slideData.gpsSeconds;
  chan.epoch.gpsNanoSeconds -= slideData.gpsNanoSeconds;
  /* copy the start time into the REAL8 h(t) time series */
  strainChan.epoch = chan.epoch;

  if ( globFrameData )
  {
    CHAR ifoRegExPattern[6];

    if ( vrbflg ) fprintf( stdout, "globbing for *.gwf frame files from %c "
        "of type %s in current directory\n", fqChanName[0], frInType );

    frGlobCache = NULL;

    /* create a frame cache by globbing all *.gwf files in the pwd */
    LAL_CALL( LALFrCacheGenerate( &status, &frGlobCache, NULL, NULL ),
        &status );

    /* check we globbed at least one frame file */
    if ( ! frGlobCache->numFrameFiles )
    {
      fprintf( stderr, "error: no frame file files of type %s found\n",
          frInType );
      exit( 1 );
    }

    /* sieve out the requested data type */
    memset( &sieve, 0, sizeof(FrCacheSieve) );
    LALSnprintf( ifoRegExPattern,
        sizeof(ifoRegExPattern) / sizeof(*ifoRegExPattern), ".*%c.*",
        fqChanName[0] );
    sieve.srcRegEx = ifoRegExPattern;
    sieve.dscRegEx = frInType;
    LAL_CALL( LALFrCacheSieve( &status, &frInCache, frGlobCache, &sieve ),
        &status );

    /* check we got at least one frame file back after the sieve */
    if ( ! frInCache->numFrameFiles )
    {
      fprintf( stderr, "error: no frame files of type %s globbed as input\n",
          frInType );
      exit( 1 );
    }

    LAL_CALL( LALDestroyFrCache( &status, &frGlobCache ), &status );
  }
  else
  {
    if ( vrbflg ) fprintf( stdout,
        "reading frame file locations from cache file: %s\n", frInCacheName );

    /* read a frame cache from the specified file */
    LAL_CALL( LALFrCacheImport( &status, &frInCache, frInCacheName), &status );
  }

  /* open the input data frame stream from the frame cache */
  LAL_CALL( LALFrCacheOpen( &status, &frStream, frInCache ), &status );

  /* set the mode of the frame stream to fail on gaps or time errors */
  frStream->mode = LAL_FR_VERBOSE_MODE;

  /* enable frame-file checksum checking */
  XLALFrSetMode( frStream, frStream->mode | LAL_FR_CHECKSUM_MODE );

  /* seek to required epoch and set chan name */
  LAL_CALL( LALFrSeek( &status, &(chan.epoch), frStream ), &status );
  frChan.name = fqChanName;

  if ( calData == real_8 )
  {
    /* determine the sample rate of the raw data */
    LAL_CALL( LALFrGetREAL8TimeSeries( &status, &strainChan, &frChan,
          frStream ), &status );

    /* copy the data paramaters from the h(t) channel to input data channel */
    LALSnprintf( chan.name, LALNameLength * sizeof(CHAR), "%s",
        strainChan.name );
    chan.epoch          = strainChan.epoch;
    chan.deltaT         = strainChan.deltaT;
    chan.f0             = strainChan.f0;
    chan.sampleUnits    = strainChan.sampleUnits;
  }
  else
  {
    /* determine the sample rate of the raw data */
    LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
        &status );
  }

  /* store the input sample rate */
  this_search_summvar = searchsummvars.searchSummvarsTable =
    (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
  LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
      "raw data sample rate" );
  this_search_summvar->value = inputDeltaT = chan.deltaT;

  /* determine if we need to resample the channel */
  if ( vrbflg )
  {
    fprintf( stdout, "resampleParams.deltaT = %e\n", resampleParams.deltaT );
    fprintf( stdout, "chan.deltaT = %e\n", chan.deltaT );
  }
  if ( ! ( fabs( resampleParams.deltaT - chan.deltaT ) < epsilon ) )
  {
    resampleChan = 1;
    if ( vrbflg )
      fprintf( stdout, "input channel will be resampled\n" );

    if ( resampFiltType == 0 )
    {
      resampleParams.filterType = LDASfirLP;
    }
    else if ( resampFiltType == 1 )
    {
      resampleParams.filterType = defaultButterworth;
    }
  }

  /* determine the number of points to get and create storage for the data */
  inputLengthNS =
    (REAL8) ( gpsEndTimeNS - gpsStartTimeNS + 2000000000LL * padData );
  numInputPoints = (UINT4) floor( inputLengthNS / (chan.deltaT * 1.0e9) + 0.5 );
  if ( calData == real_8 )
  {
    /* create storage for the REAL8 h(t) input data */
    LAL_CALL( LALDCreateVector( &status, &(strainChan.data), numInputPoints ),
        &status );
  }
  LAL_CALL( LALSCreateVector( &status, &(chan.data), numInputPoints ),
      &status );

  if ( vrbflg ) fprintf( stdout, "input channel %s has sample interval "
      "(deltaT) = %e\nreading %d points from frame stream\n", fqChanName,
      chan.deltaT, numInputPoints );

  if ( calData == real_8 )
  {
    /* read in the REAL8 h(t) data here */
    PassBandParamStruc strainHighpassParam;

    /* read the REAL8 h(t) data from the time series into strainChan      */
    /* which already has the correct amount of memory allocated */
    if ( vrbflg ) fprintf( stdout, "reading REAL8 h(t) data from frames... " );

    LAL_CALL( LALFrGetREAL8TimeSeries( &status, &strainChan, &frChan,
          frStream ), &status);

    if ( vrbflg ) fprintf( stdout, "done\n" );

    /* high pass the h(t) data using the parameters specified on the cmd line*/
    strainHighpassParam.nMax = strainHighPassOrder;
    strainHighpassParam.f1 = -1.0;
    strainHighpassParam.f2 = (REAL8) strainHighPassFreq;
    strainHighpassParam.a1 = -1.0;
    strainHighpassParam.a2 = (REAL8)(1.0 - strainHighPassAtten);
    if ( vrbflg ) fprintf( stdout,
        "applying %d order high pass to REAL8 h(t) data: "
        "%3.2f of signal passes at %4.2f Hz\n",
        strainHighpassParam.nMax, strainHighpassParam.a2,
        strainHighpassParam.f2 );

    LAL_CALL( LALButterworthREAL8TimeSeries( &status, &strainChan,
          &strainHighpassParam ), &status );

    /* cast the REAL8 h(t) data to REAL4 in the chan time series       */
    /* which already has the correct amount of memory allocated */
    for ( j = 0 ; j < numInputPoints ; ++j )
    {
      chan.data->data[j] = (REAL4) ( strainChan.data->data[j] * dynRange );
    }

    /* re-copy the data paramaters from the h(t) channel to input data channel*/
    LALSnprintf( chan.name, LALNameLength * sizeof(CHAR), "%s",
        strainChan.name );
    chan.epoch          = strainChan.epoch;
    chan.deltaT         = strainChan.deltaT;
    chan.f0             = strainChan.f0;
    chan.sampleUnits    = strainChan.sampleUnits;

    /* free the REAL8 h(t) input data */
    LAL_CALL( LALDDestroyVector( &status, &(strainChan.data) ), &status );
    strainChan.data = NULL;
  }
  else
  {
    /* read the data channel time series from frames */
    LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
        &status );
    if ( calData == real_4 )
    {
      /* multiply the input data by dynRange */
      for ( j = 0 ; j < numInputPoints ; ++j )
      {
        chan.data->data[j] *= dynRange;
      }
    }
  }
  memcpy( &(chan.sampleUnits), &lalADCCountUnit, sizeof(LALUnit) );

  /* store the start and end time of the raw channel in the search summary */
  searchsumm.searchSummaryTable->in_start_time = chan.epoch;
  LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan.epoch) ),
      &status );
  tsLength += chan.deltaT * (REAL8) chan.data->length;
  LAL_CALL( LALFloatToGPS( &status,
        &(searchsumm.searchSummaryTable->in_end_time), &tsLength ), &status );

  /* close the frame file stream and destroy the cache */
  LAL_CALL( LALFrClose( &status, &frStream ), &status );
  LAL_CALL( LALDestroyFrCache( &status, &frInCache ), &status );

  /* write the raw channel data as read in from the frame files */
  if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
      &chan, "ct", "RAW" );

  if ( vrbflg ) fprintf( stdout, "read channel %s from frame stream\n"
      "got %d points with deltaT %e\nstarting at GPS time %d sec %d ns\n",
      chan.name, chan.data->length, chan.deltaT,
      chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );


  /*
   *
   * create the random seed if it will be needed
   *
   */


  if ( randSeedType != unset )
  {
    /* store the seed in the search summvars table */
    this_search_summvar = this_search_summvar->next =
      (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
    LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
        "template bank simulation seed" );

    if ( randSeedType == urandom )
    {
      FILE   *fpRand = NULL;
      INT4    randByte;

      if ( vrbflg )
        fprintf( stdout, "obtaining random seed from /dev/urandom: " );

      randomSeed = 0;
      fpRand = fopen( "/dev/urandom", "r" );
      if ( fpRand )
      {
        for ( randByte = 0; randByte < 4 ; ++randByte )
        {
          INT4 tmpSeed = (INT4) fgetc( fpRand );
          randomSeed += tmpSeed << ( randByte * 8 );
        }
        fclose( fpRand );
      }
      else
      {
        perror( "error obtaining random seed from /dev/urandom" );
        exit( 1 );
      }
    }
    else if ( randSeedType == user )
    {
      if ( vrbflg )
        fprintf( stdout, "using user specified random seed: " );
    }
    else
    {
      /* should never get here */
      fprintf( stderr, "error obtaining random seed\n" );
      exit( 1 );
    }

    this_search_summvar->value = randomSeed;
    if ( vrbflg ) fprintf( stdout, "%d\n", randomSeed );

    /* create the tmplt bank random parameter structure */
    LAL_CALL( LALCreateRandomParams( &status, &randParams, randomSeed ),
        &status );
  }

  /* replace the input data with Gaussian noise if necessary */
  if ( whiteGaussian )
  {
    if ( vrbflg ) fprintf( stdout,
        "setting input data to white Gaussian noise with variance %e... ",
        gaussVar );
    memset( chan.data->data, 0, chan.data->length * sizeof(REAL4) );
    LAL_CALL( LALNormalDeviates( &status, chan.data, randParams ), &status );
    for ( j = 0; j < chan.data->length; ++j )
    {
      chan.data->data[j] *= gaussVar;
    }
    if ( vrbflg ) fprintf( stdout, "done\n" );

    /* write the raw channel data as read in from the frame files */
    if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
        &chan, "ct", "RAW_GAUSSIAN" );
  }

  /* replace data with colored Gaussian noise if requested */
  if ( coloredGaussian )
  {
    COMPLEX8Sequence *ntilde         = NULL;
    REAL4Sequence    *ntilde_re      = NULL;
    REAL4Sequence    *ntilde_im      = NULL;
    REAL8Sequence    *spectrum       = NULL;
    REAL4FFTPlan     *invPlan        = NULL;
    INT4              length         = chan.data->length;
    REAL8             deltaT         = chan.deltaT;
    REAL8             deltaF         = 1.0 / (deltaT * (REAL8) length);
    INT4              kmin           = strainHighPassFreq / deltaF > 1 ?
                                       strainHighPassFreq / deltaF : 1 ;

    if ( vrbflg ) fprintf( stdout, "replacing data with colored Gaussian noise...." );

    /* Generate white Gaussian noise with unit variance */

    LAL_CALL( LALSCreateVector( &status, &ntilde_re, length / 2 + 1 ),
        &status );
    LAL_CALL( LALNormalDeviates( &status, ntilde_re, randParams ), &status );
    LAL_CALL( LALSCreateVector( &status, &ntilde_im, length / 2 + 1 ),
        &status );
    LAL_CALL( LALNormalDeviates( &status, ntilde_im, randParams ), &status );

    /* create storage for the frequency domain noise and psd*/
    LAL_CALL( LALCCreateVector( &status, &ntilde, length / 2 + 1 ),
        &status );
    LAL_CALL( LALDCreateVector( &status, &spectrum, length / 2 + 1 ),
        &status );
    if (colorSpec == 3 )
    {
      /* set the spectrum to the Initial LIGO design noise curve */
      REAL8 psd_value;
      LALLIGOIPsd( NULL, &psd_value, strainHighPassFreq );
      for ( k = 0; k < (UINT4)kmin ; ++k )
      {
        spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
      }
      for ( k = kmin; k < spectrum->length ; ++k )
      {
        REAL8 psd_freq = (REAL8) k * deltaF;
        LALLIGOIPsd( NULL, &psd_value, psd_freq );
        spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
      }
    }
    else if( colorSpec == 4)
    {
      /* set the spectrum to the Advanced LIGO design noise curve */
      REAL8 psd_value;
      LALAdvLIGOPsd( NULL, &psd_value, strainHighPassFreq );
      for ( k = 0; k < (UINT4)kmin ; ++k )
      {
        spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
      }
      for ( k = kmin; k < spectrum->length ; ++k )
      {
        REAL8 psd_freq = (REAL8) k * deltaF;
        LALAdvLIGOPsd( NULL, &psd_value, psd_freq );
        spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
      }
    }

    /* Color white noise with given psd */
    for ( k=0; k < ntilde->length; k++ )
    {
      ntilde->data[k].re = ntilde_re->data[k] * sqrt(( (REAL4) length * 0.25 /
                           (REAL4) deltaT ) * (REAL4) spectrum->data[k] );
      ntilde->data[k].im = ntilde_im->data[k] * sqrt(( (REAL4) length * 0.25 /
                           (REAL4) deltaT ) * (REAL4) spectrum->data[k] );
    }
    /* setting d.c. and Nyquist to zero */
    ntilde->data[0].im = ntilde->data[length / 2].im = 0.0;

    /* Fourier transform back in the time domain */
    LAL_CALL( LALCreateReverseRealFFTPlan( &status, &invPlan, length, 0 ),
        &status );
    LAL_CALL( LALReverseRealFFT( &status, chan.data, ntilde, invPlan ), &status);
    /* normalize the noise */
    for ( j = 0; j < (UINT4)length; ++j )
    {
      chan.data->data[j] /= (REAL4) length ;
    }

    if ( vrbflg ) fprintf( stdout, "done\n" );

    /* write the simulated noise to frame file*/
    if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
        &chan, "ct", "RAW_COLORED_GAUSSIAN" );

    LALSDestroyVector( &status, &ntilde_re );
    LALSDestroyVector( &status, &ntilde_im );
    LALCDestroyVector( &status, &ntilde );
    LALDDestroyVector( &status, &spectrum);
    LALDestroyRealFFTPlan( &status, &invPlan);
  }


  /*
   *
   * generate the response function for the requested time
   *
   */

  /* create storage for the response function */

  memset( &resp, 0, sizeof(COMPLEX8FrequencySeries) );
  LAL_CALL( LALCCreateVector( &status, &(resp.data), numPoints / 2 + 1 ),
      &status );

  /* set the parameters of the response to match the data */
  resp.epoch.gpsSeconds = chan.epoch.gpsSeconds + padData;
  resp.epoch.gpsNanoSeconds = chan.epoch.gpsNanoSeconds;
  resp.deltaF = (REAL8) sampleRate / (REAL8) numPoints;
  resp.sampleUnits = strainPerCount;
  strcpy( resp.name, chan.name );

  /* generate the response function for the current time */
  if ( vrbflg ) fprintf( stdout, "generating response at time %d sec %d ns\n",
      resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds );

  /* initialize the calfacts */
  memset( &calfacts, 0, sizeof(CalibrationUpdateParams) );
  calfacts.ifo = ifo;

  /* determine length of chunk */
  if ( pointCal )
  {
    calfacts.duration.gpsSeconds = 1;
    calfacts.duration.gpsNanoSeconds = 0;
  }
  else
  {
    durationNS = gpsEndTimeNS - gpsStartTimeNS;
    XLALINT8NSToGPS( &(calfacts.duration), durationNS );
  }

  if ( calData )
  {
    /* if we are using calibrated data set the response to unity */
    for( k = 0; k < resp.data->length; ++k )
    {
      resp.data->data[k].re = (REAL4) (1.0 / dynRange);
      resp.data->data[k].im = 0.0;
    }
    if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries(
        outFrame, &resp, "strain/ct", "RESPONSE_h(t)" );
  }
  else
  {
    /* create the lal calibration frame cache */
    if ( globCalData )
    {
      calGlobPattern = (CHAR *) LALCalloc( calGlobLen, sizeof(CHAR) );
      LALSnprintf( calGlobPattern, calGlobLen * sizeof(CHAR),
          "*CAL*%s*.gwf", ifo );
      if ( vrbflg ) fprintf( stdout, "globbing for %s calibration frame files "
          "in current directory\n", calGlobPattern );
    }
    else
    {
      calGlobPattern = NULL;
      if ( vrbflg ) fprintf( stdout,
          "reading calibration data from cache: %s\n", calCacheName );
    }

    LAL_CALL( LALCreateCalibFrCache( &status, &calCache, calCacheName,
          NULL, calGlobPattern ), &status );

    if ( calGlobPattern ) LALFree( calGlobPattern );

    /* store the name of the calibration files used */
    for ( i = 0; i < calCache->numFrameFiles; ++i )
    {
      this_search_summvar = this_search_summvar->next =
        (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
      LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
          "calibration frame %d", i );
      LALSnprintf( this_search_summvar->string,
          LIGOMETA_STRING_MAX * sizeof(CHAR), "%s",
          calCache->frameFiles[i].url );
    }

    /* get the response from the frame data */
    LAL_CALL( LALExtractFrameResponse( &status, &resp, calCache,
          &calfacts), &status );
    LAL_CALL( LALDestroyFrCache( &status, &calCache), &status );
    alpha = (REAL4) calfacts.alpha.re;
    alphabeta = (REAL4) calfacts.alphabeta.re;
    if ( vrbflg ) fprintf( stdout,
        "for calibration of data, alpha = %f and alphabeta = %f\n",
        alpha, alphabeta);

    if ( writeResponse )
      outFrame = fr_add_proc_COMPLEX8FrequencySeries( outFrame, &resp,
          "strain/ct", "RESPONSE" );
  }

  if ( unitResponse )
  {
    /* replace the response function with unity if */
    /* we are filtering white Gaussian noise             */
    if ( vrbflg ) fprintf( stdout, "setting response to unity... " );
    for ( k = 0; k < resp.data->length; ++k )
    {
      resp.data->data[k].re = 1.0;
      resp.data->data[k].im = 0;
    }
    if ( vrbflg ) fprintf( stdout, "done\n" );

    if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries(
        outFrame, &resp, "strain/ct", "RESPONSE_UNITY" );
  }
  if ( ( specType == specType_LIGO ) || ( specType == specType_AdvLIGO ) )
  {
    /* replace the response function with 1/dynRange if we are using the */
    /* design LIGO or AdvLIGO psd                                        */
    if ( vrbflg ) fprintf( stdout, "setting response to inverse dynRange... " );
    for( k = 0; k < resp.data->length; ++k )
    {
      resp.data->data[k].re = (REAL4) (1.0 / dynRange);
      resp.data->data[k].im = 0.0;
    }
    if ( vrbflg ) fprintf( stdout, "done\n" );
    if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries(
        outFrame, &resp, "strain/ct", "RESPONSE_UNITY" );
  }

  /* slide the channel back to the fake time for background studies */
  chan.epoch.gpsSeconds += slideData.gpsSeconds;
  chan.epoch.gpsNanoSeconds += slideData.gpsNanoSeconds;


  if ( tdFollowUpFiles )
  {
    INT4 injSafety = 10;
    LIGOTimeGPS startKeep = gpsStartTime;
    LIGOTimeGPS endKeep   = gpsEndTime;

    startKeep.gpsSeconds -= injSafety;
    endKeep.gpsSeconds   += injSafety;

    /* read in the time domain follow-up data from XML */
    for (i = 0; i < (UINT4)numTDFiles; i++ )
    {
      INT4 thisTDNum = 0;
      if ( !tdFollowUpEvents )
      {
        thisTDNum = LALSnglInspiralTableFromLIGOLw(&tdFollowUpEvents,
             tdFollowUpFiles[i], 0, -1);
        thisFollowUpEvent = tdFollowUpEvents;
      }
      else
      {
        thisTDNum = LALSnglInspiralTableFromLIGOLw(&(thisFollowUpEvent->next),
             tdFollowUpFiles[i], 0, -1);
      }
      if ( thisTDNum < 0 )
      {
        fprintf( stderr, "Error reading file %s\n", tdFollowUpFiles[i] );
        exit( 1 );
      }

      if ( thisTDNum )
      {
        while ( thisFollowUpEvent->next )
        {
          thisFollowUpEvent = thisFollowUpEvent->next;
        }
      }
    }

    tdFollowUpEvents = XLALTimeCutSingleInspiral( tdFollowUpEvents,
        &startKeep, &endKeep);
  }

  /*
   *
   * inject signals into the raw, unresampled data
   *
   */


  if ( injectionFile )
  {
    /* get injections within 500 seconds of either end of the segment.   */
    /* a 0.4,0.4 MACHO starting at 30.0 Hz has length 435.374683 seconds */
    /* so this should be plenty of safety. better to waste cpu than miss */
    /* injected signals...                                               */
    INT4 injSafety = 500;

    /* read in the injection data from XML */
    numInjections = SimInspiralTableFromLIGOLw( &injections, injectionFile,
        gpsStartTime.gpsSeconds - injSafety,
        gpsEndTime.gpsSeconds + injSafety );

    if ( numInjections < 0 )
    {
      fprintf( stderr, "error: cannot read injection file" );
      exit( 1 );
    }
    else if ( numInjections && !hardwareInjection )
    {
      /* see if we need a higher resolution response to do the injections */
      if ( resampleChan )
      {
        /* we need a different resolution of response function for injections */
        UINT4 rateRatio = floor( resampleParams.deltaT / chan.deltaT + 0.5 );
        UINT4 rawNumPoints = rateRatio * numPoints;

        if ( vrbflg ) fprintf( stdout, "rateRatio = %d\nrawNumPoints = %d\n"
            "chan.deltaT = %e\n", rateRatio, rawNumPoints, chan.deltaT );

        memset( &injResp, 0, sizeof(COMPLEX8FrequencySeries) );
        LAL_CALL( LALCCreateVector( &status, &(injResp.data),
              rawNumPoints / 2 + 1 ), &status );
        injResp.epoch = resp.epoch;
        injResp.deltaF = 1.0 / ( rawNumPoints * chan.deltaT );
        injResp.sampleUnits = strainPerCount;
        strcpy( injResp.name, chan.name );

        if ( calData || coloredGaussian )
        {
          /* if we are using calibrated data set the response to unity */
          if ( vrbflg ) fprintf( stdout,
              "setting injection response to inverse dynRange..." );
          for ( k = 0; k < injResp.data->length; ++k )
          {
            injResp.data->data[k].re = (REAL4)(1.0/dynRange);
            injResp.data->data[k].im = 0.0;
          }
          injRespPtr = &injResp;
          if ( vrbflg ) fprintf( stdout, "done\n" );
          if ( writeResponse )
            outFrame = fr_add_proc_COMPLEX8FrequencySeries( outFrame,
                &injResp, "strain/ct", "RESPONSE_INJ_CAL" );
        }
        else
        {
          /* generate the response function for the current time */
          if ( vrbflg ) fprintf( stdout,
              "generating high resolution response at time %d sec %d ns\n"
              "length = %d points, deltaF = %e Hz\n",
              resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds,
              injResp.data->length, injResp.deltaF );

          /* initialize the inj_calfacts */
          memset( &inj_calfacts, 0, sizeof(CalibrationUpdateParams) );
          inj_calfacts.ifo = ifo;
          XLALINT8NSToGPS( &(inj_calfacts.duration), durationNS );

          /* create the lal calibration frame cache */
          if ( globCalData )
          {
            calGlobPattern = (CHAR *) LALCalloc( calGlobLen, sizeof(CHAR) );
            LALSnprintf( calGlobPattern, calGlobLen * sizeof(CHAR),
                "*CAL*%s*.gwf", ifo );
            if ( vrbflg ) fprintf( stdout,
                "globbing for %s calibration frame files "
                "in current directory\n", calGlobPattern );
          }
          else
          {
            calGlobPattern = NULL;
            if ( vrbflg ) fprintf( stdout,
                "reading calibration data from cache: %s\n", calCacheName );
          }

          LAL_CALL( LALCreateCalibFrCache( &status, &calCache, calCacheName,
                NULL, calGlobPattern ), &status );

          if ( calGlobPattern ) LALFree( calGlobPattern );

          /* store the name of the calibration files used */
          for ( i = 0; i < calCache->numFrameFiles; ++i )
          {
            this_search_summvar = this_search_summvar->next =
              (SearchSummvarsTable *)
              LALCalloc( 1, sizeof(SearchSummvarsTable) );
            LALSnprintf( this_search_summvar->name,
                LIGOMETA_NAME_MAX * sizeof(CHAR),
                "injection calibration frame %d", i );
            LALSnprintf( this_search_summvar->string,
                LIGOMETA_STRING_MAX * sizeof(CHAR), "%s",
                calCache->frameFiles[i].url );
          }

          /* extract the calibration from frames */
          LAL_CALL( LALExtractFrameResponse( &status, &injResp, calCache,
                &inj_calfacts ), &status );
          LAL_CALL( LALDestroyFrCache( &status, &calCache), &status );

          inj_alpha = (REAL4) calfacts.alpha.re;
          inj_alphabeta = (REAL4) calfacts.alphabeta.re;
          if ( vrbflg ) fprintf( stdout,
              "for injections, alpha = %f and alphabeta = %f\n",
              inj_alpha, inj_alphabeta);

          injRespPtr = &injResp;

          if ( writeResponse )
            outFrame = fr_add_proc_COMPLEX8FrequencySeries(
                outFrame, &injResp, "strain/ct", "RESPONSE_INJ" );
        }

        if ( whiteGaussian )
        {
          /* replace the response function with unity if */
          /* we are filtering white Gaussian noise             */
          if ( vrbflg ) fprintf( stdout, "setting response to unity... " );
          for ( k = 0; k < injResp.data->length; ++k )
          {
            injResp.data->data[k].re = 1.0;
            injResp.data->data[k].im = 0;
          }
          if ( vrbflg ) fprintf( stdout, "done\n" );

          if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries(
              outFrame, &injResp, "strain/ct", "RESPONSE_INJ_GAUSSIAN" );
        }
      }
      else
      {
        /* the data is already at the correct sample rate, just do injections */
        injRespPtr = &resp;
        memset( &injResp, 0, sizeof(COMPLEX8FrequencySeries) );
      }

      /* inject the signals, preserving the channel name (Tev mangles it) */
      LALSnprintf( tmpChName, LALNameLength * sizeof(CHAR), "%s", chan.name );

      /* if injectOverhead option, then set chan.name to "ZENITH".
       * This causes no detector site to be found in the injection code so
       * that the injection is done directly overhead (i.e. with a response
       * function of F+ = 1; Fx = 0) */
      if ( injectOverhead )
      {
        LALSnprintf( chan.name, LALNameLength * sizeof(CHAR), "ZENITH" );
      }

      /* read the event waveform approximant to see if we've been asked to
       perform NumRel injections */
      LAL_CALL( LALGetApproximantFromString( &status, injections->waveform,
				  &injApproximant ), &status);

      if (injApproximant == NumRel)
      {
        REAL4TimeVectorSeries *tempStrain=NULL;

	/* loop over injections */
        for ( thisInj = injections; thisInj; thisInj = thisInj->next )
	  {
            LAL_CALL( AddNumRelStrainModes( &status, &tempStrain, thisInj),
                      &status);

            LAL_CALL( LALInjectStrainGW( &status, &chan, tempStrain, thisInj,
                      ifo, dynRange), &status);

            XLALDestroyREAL4VectorSequence ( tempStrain->data);
            tempStrain->data = NULL;
            LALFree(tempStrain);
            tempStrain = NULL;

          } /* loop over injectionsj */
      }
      else
      {
	LAL_CALL( LALFindChirpInjectSignals( &status, &chan, injections,
					       injRespPtr ), &status );
      }
      LALSnprintf( chan.name,  LALNameLength * sizeof(CHAR), "%s", tmpChName );

      if ( vrbflg ) fprintf( stdout, "injected %d signals from %s into %s\n",
          numInjections, injectionFile, chan.name );

      /* write the raw channel data plus injections to the output frame file */
      if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
          &chan, "ct", "RAW_INJ" );

      if ( injResp.data )
        LAL_CALL( LALCDestroyVector( &status, &(injResp.data) ), &status );
    }
    else
    {
      if ( vrbflg && !hardwareInjection ) fprintf( stdout,
              "no injections in this chunk\n" );
    }
  }


  /*
   *
   * resample the data to the requested rate
   *
   */


  if ( resampleChan )
  {
    if (vrbflg) fprintf( stdout, "resampling input data from %e to %e\n",
        chan.deltaT, resampleParams.deltaT );

    LAL_CALL( LALResampleREAL4TimeSeries( &status, &chan, &resampleParams ),
        &status );

    if ( vrbflg ) fprintf( stdout, "channel %s resampled:\n"
        "%d points with deltaT %e\nstarting at GPS time %d sec %d ns\n",
        chan.name, chan.data->length, chan.deltaT,
        chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

    /* write the resampled channel data as read in from the frame files */
    if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
        &chan, "ct", "RAW_RESAMP" );
  }

  /* store the filter data sample rate */
  this_search_summvar = this_search_summvar->next =
    (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
  LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
      "filter data sample rate" );
  this_search_summvar->value = chan.deltaT;


  /*
   *
   * we have all the input data that we need, so checkpoint if requested
   *
   */


  if ( dataCheckpoint )
  {
#ifdef LALAPPS_CONDOR
    condor_compress_ckpt = 1;
    if ( ckptPath[0] )
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s/%s.ckpt",
          ckptPath, fileName );
    }
    else
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s.ckpt", fileName );
    }
    if ( vrbflg ) fprintf( stdout, "checkpointing to file %s\n", fname );
    init_image_with_file_name( fname );
    ckpt_and_exit();
#else
    fprintf( stderr, "--data-checkpoint cannot be used unless "
        "lalapps is condor compiled\n" );
    exit( 1 );
#endif
  }


  /*
   *
   * high pass the data, removed pad from time series and check length of data
   *
   */


  /* iir filter to remove low frequencies from data channel */
  if ( highPass )
  {
    PassBandParamStruc highpassParam;
    highpassParam.nMax = highPassOrder;
    highpassParam.f1 = -1.0;
    highpassParam.f2 = (REAL8) highPassFreq;
    highpassParam.a1 = -1.0;
    highpassParam.a2 = (REAL8)(1.0 - highPassAtten); /* a2 is not attenuation */

    if ( vrbflg ) fprintf( stdout, "applying %d order high pass: "
        "%3.2f of signal passes at %4.2f Hz\n",
        highpassParam.nMax, highpassParam.a2, highpassParam.f2 );

    LAL_CALL( LALDButterworthREAL4TimeSeries( &status, &chan, &highpassParam ),
        &status );
  }

  /* remove pad from requested data from start and end of time series */
  memmove( chan.data->data, chan.data->data + padData * sampleRate,
      (chan.data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
  LALRealloc( chan.data->data,
      (chan.data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
  chan.data->length -= 2 * padData * sampleRate;
  chan.epoch.gpsSeconds += padData;

  if ( vrbflg ) fprintf( stdout, "after removal of %d second padding at "
      "start and end:\ndata channel sample interval (deltaT) = %e\n"
      "data channel length = %d\nstarting at %d sec %d ns\n",
      padData , chan.deltaT , chan.data->length,
      chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

  if ( writeFilterData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
      &chan, "ct", "FILTER" );

  /* check data length */
  if ( chan.data->length != inputDataLength )
  {
    fprintf( stderr, "error: computed channel length and requested\n"
        "input data length do not match:\nchan.data->length = %d\n"
        "inputDataLength = %d\nyou have found a bug in the code.\n"
        "please report this to <duncan@gravity.phys.uwm.edu>\n",
        chan.data->length, inputDataLength );
    exit( 1 );
  }

  /* store the start and end time of the filter channel in the search summ */
  /* noting that we don't look for events in the first and last quarter    */
  /* of each findchirp segment of the input data                           */
  LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan.epoch) ),
      &status );
  tsLength += (REAL8) (numPoints / 4) * chan.deltaT;
  LAL_CALL( LALFloatToGPS( &status,
        &(searchsumm.searchSummaryTable->out_start_time), &tsLength ),
      &status );

  outTimeNS = XLALGPSToINT8NS( &(searchsumm.searchSummaryTable->out_start_time) );

  if ( ! bankSim && ( trigStartTimeNS && (trigStartTimeNS > outTimeNS) ) )
  {
    /* override with trigger start time */
    XLALINT8NSToGPS( &(searchsumm.searchSummaryTable->out_start_time),
          trigStartTimeNS );
  }

  LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan.epoch) ),
      &status );
  tsLength += chan.deltaT * ((REAL8) chan.data->length - (REAL8) (numPoints/4));
  LAL_CALL( LALFloatToGPS( &status,
        &(searchsumm.searchSummaryTable->out_end_time), &tsLength ),
      &status );

  outTimeNS = XLALGPSToINT8NS( &(searchsumm.searchSummaryTable->out_end_time) );

  if ( ! bankSim && ( trigEndTimeNS && (trigEndTimeNS < outTimeNS) ) )
  {
    /* override with trigger end time */
    XLALINT8NSToGPS( &(searchsumm.searchSummaryTable->out_end_time),
          trigEndTimeNS );
  }

  /*
   *
   * create and populate findchip initialization structure
   *
   */


  if ( ! ( fcInitParams = (FindChirpInitParams *)
        LALCalloc( 1, sizeof(FindChirpInitParams) ) ) )
  {
    fprintf( stderr, "could not allocate memory for findchirp init params\n" );
    exit( 1 );
  }

  fcInitParams->numPoints      = numPoints;
  fcInitParams->numSegments    = numSegments;
  fcInitParams->numChisqBins   = numChisqBins;
  fcInitParams->createRhosqVec = writeRhosq;
  fcInitParams->ovrlap         = ovrlap;
  fcInitParams->approximant    = approximant;
  fcInitParams->order          = order;
  fcInitParams->createCVec     = writeCData;


  /*
   *
   * power spectrum estimation
   *
   */


  /* initialize findchirp data conditioning routine */
  LAL_CALL( LALFindChirpDataInit( &status, &fcDataParams, fcInitParams ),
      &status );
  fcDataParams->invSpecTrunc = invSpecTrunc * sampleRate;
  fcDataParams->fLow = fLow;

  /* create storage for the power spectral estimate */
  memset( &spec, 0, sizeof(REAL4FrequencySeries) );
  LAL_CALL( LALSCreateVector( &status, &(spec.data), numPoints / 2 + 1 ),
      &status );

  /* compute the windowed power spectrum for the data channel */
  avgSpecParams.window = NULL;
  switch ( specType )
  {
    case specType_mean:
      avgSpecParams.method = useMean;
      if ( vrbflg ) fprintf( stdout, "computing mean psd" );
      break;
    case specType_median:
      avgSpecParams.method = useMedian;
      if ( vrbflg ) fprintf( stdout, "computing median psd" );
      break;
    case specType_gaussian:
    case specType_LIGO:
    case specType_AdvLIGO:
      avgSpecParams.method = useUnity;
      if ( vrbflg ) fprintf( stdout, "computing constant psd with unit value" );
      break;
    default:
      fprintf( stderr, "Error: unknown PSD estimation type: %d\n", specType );
      exit( 1 );
  }

  /* use the fft plan created by findchirp */
  avgSpecParams.plan = fcDataParams->fwdPlan;

  wpars.type = Hann;
  wpars.length = numPoints;
  if ( badMeanPsd )
  {
    avgSpecParams.overlap = 0;
    if ( vrbflg ) fprintf( stdout, " without overlap\n" );
  }
  else
  {
    avgSpecParams.overlap = numPoints / 2;
    if ( vrbflg )
      fprintf( stdout, " with overlap %d\n", avgSpecParams.overlap );
  }

  LAL_CALL( LALCreateREAL4Window( &status, &(avgSpecParams.window),
        &wpars ), &status );
  LAL_CALL( LALREAL4AverageSpectrum( &status, &spec, &chan, &avgSpecParams ),
      &status );
  XLALDestroyREAL4Window( avgSpecParams.window );
  strcpy( spec.name, chan.name );

  if ( specType == specType_gaussian )
  {
    /* multiply the unit power spectrum by the variance to get the white psd */
    REAL4 gaussVarSq = gaussVar * gaussVar;
    for ( k = 0; k < spec.data->length; ++k )
    {
      spec.data->data[k] *= 2.0 * gaussVarSq * (REAL4) inputDeltaT;
    }
    if ( vrbflg )
      fprintf( stdout, "set psd to constant value = %e\n", spec.data->data[0] );
  }
  else if ( specType == specType_LIGO )
  {
    /* replace the spectrum with the Initial LIGO design noise curve */
    INT4  k_min          = strainHighPassFreq / spec.deltaF > 1 ?
                           strainHighPassFreq / spec.deltaF : 1 ;
    REAL8 sim_psd_value, sim_psd_freq;
    LALLIGOIPsd( NULL, &sim_psd_value, strainHighPassFreq );

    for ( k = 0; k < (UINT4)k_min; ++k )
    {
      spec.data->data[k] = (REAL4) (9.0e-46 * sim_psd_value * dynRange *
                                    dynRange);
    }
    for ( k = k_min; k < spec.data->length; ++k )
    {
      sim_psd_freq = (REAL8) k * spec.deltaF;
      LALLIGOIPsd( NULL, &sim_psd_value, sim_psd_freq );
      spec.data->data[k] = (REAL4) (9.0e-46 * sim_psd_value * dynRange *
                                    dynRange);
    }
    if ( vrbflg ) fprintf( stdout, "set psd to Initial LIGO design\n" );
  }
  else if ( specType == specType_AdvLIGO )
  {
    /* replace the spectrum with the Advanced LIGO design noise curve */
    INT4  k_min          = strainHighPassFreq / spec.deltaF > 1 ?
                           strainHighPassFreq / spec.deltaF : 1 ;
    REAL8 sim_psd_value, sim_psd_freq;
    LALAdvLIGOPsd( NULL, &sim_psd_value, strainHighPassFreq );

    for ( k = 0; k < (UINT4)k_min; ++k )
    {
      spec.data->data[k] = (REAL4) (9.0e-46 * sim_psd_value * dynRange *
                                    dynRange);
    }
    for ( k = k_min; k < spec.data->length; ++k )
    {
      sim_psd_freq = (REAL8) k * spec.deltaF;
      LALAdvLIGOPsd( NULL, &sim_psd_value, sim_psd_freq );
      spec.data->data[k] = (REAL4) (9.0e-46 * sim_psd_value * dynRange *
                                    dynRange);
    }
    if ( vrbflg ) fprintf( stdout, "set psd to Advanced LIGO design\n" );
  }

  /* write the spectrum data to a file */
  if ( writeSpectrum ) outFrame = fr_add_proc_REAL4FrequencySeries( outFrame,
      &spec, "ct^2/Hz", "PSD" );


  /*
   *
   * initialize the template bank simulation
   *
   */


  if ( bankSim )
  {
    if ( vrbflg ) fprintf( stdout, "initializing template bank simulation\n" );

    /* use the created random parameters */
    bankSimParams.randParams = randParams;

    /* override the initialization parameters for a bank simulation */
    fcInitParams->numSegments  = 1;
    fcInitParams->ovrlap       = 0;
    fcInitParams->numChisqBins = 0;
    snrThresh                  = 0;
    chisqThresh                = LAL_REAL4_MAX;
    chisqDelta                 = LAL_REAL4_MAX;

    if ( bankSimFileName )
    {
      /* read in the simuation parameters from a sim_inspiral_table */
      bankSim = SimInspiralTableFromLIGOLw( &bankSimHead, bankSimFileName,
          0, 0 );
      for ( thisBankSim = bankSimHead; thisBankSim;
          thisBankSim = thisBankSim->next )
      {
        /* set the time of the injection to zero so it is injected in  */
        /* the middle of the data segment: the other times are ignored */
        thisBankSim->geocent_end_time.gpsSeconds = 0;
        thisBankSim->geocent_end_time.gpsNanoSeconds = 0;
      }
      thisBankSim = bankSimHead;
      if ( vrbflg ) fprintf( stdout, "Read %d bank sim parameters from %s\n",
            bankSim, bankSimFileName );
    }
    else
    {
      /* use the bank sim params to generate the simulated signals */
      bankSimParamsPtr = &bankSimParams;
      if ( vrbflg )
        fprintf( stdout, "internally generating bank sim siganls\n" );
    }

    bankSimCutLowIndex = XLALFindChirpBankSimInitialize( &spec, &resp, fLow );

    if ( vrbflg )
      fprintf( stdout, "psd low frequency cutoff index = %d\n", bankSimCutLowIndex );
    if ( writeSpectrum ) outFrame = fr_add_proc_REAL4FrequencySeries( outFrame,
        &spec, "strain^2/Hz", "PSD_SIM" );
    if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries(
        outFrame, &resp, "strain/ct", "RESPONSE_SIM" );
  }


  /*
   *
   * create the data structures needed for findchirp
   *
   */


  if ( vrbflg ) fprintf( stdout, "initializing findchirp... " );

  /* create the data segment vector from the input data */
  LAL_CALL( LALInitializeDataSegmentVector( &status, &dataSegVec,
        &chan, &spec, &resp, fcInitParams ), &status );

  /* set the analyze flag according to what we are doing */
  if ( tdFollowUpFiles || injectionFile )
  {
    if (tdFollowUpFiles)
    {
      /* Only analyze segments containing coincident BCV triggers */
      XLALFindChirpSetFollowUpSegment (dataSegVec, &tdFollowUpEvents);
    }
    else
    {
      /* set the analyzeSegment flag only on segments with injections */
      XLALFindChirpSetAnalyzeSegment (dataSegVec, injections);
    }
  }

  /* create the findchirp data storage */
  LAL_CALL( LALCreateFindChirpSegmentVector( &status, &fcSegVec,
        fcInitParams ), &status );

  /* initialize the template functions */
  LAL_CALL( LALFindChirpTemplateInit( &status, &fcTmpltParams,
        fcInitParams ), &status );

  fcDataParams->dynRange = fcTmpltParams->dynRange = dynRange;
  fcTmpltParams->deltaT = chan.deltaT;
  fcTmpltParams->fLow = fLow;
  fcTmpltParams->reverseChirpBank = reverseChirpBank;
  fcTmpltParams->order = order;
  fcTmpltParams->bandPassTmplt = bandPassTmplt;
  fcTmpltParams->taperTmplt = taperTmplt;

  /* initialize findchirp filter functions */
  LAL_CALL( LALFindChirpFilterInit( &status, &fcFilterParams, fcInitParams ),
      &status );
  fcFilterParams->deltaT = chan.deltaT;
  fcFilterParams->chisqParams->approximant = approximant;

  /* set up parameters for the filter output veto */
  if ( enableRsqVeto )
  {
    fcFilterParams->filterOutputVetoParams = (FindChirpFilterOutputVetoParams *)
      LALCalloc( 1, sizeof(FindChirpFilterOutputVetoParams) );
    fcFilterParams->filterOutputVetoParams->rsqvetoWindow = rsqVetoWindow;
    fcFilterParams->filterOutputVetoParams->rsqvetoThresh = rsqVetoThresh;
    if ( doRsqVeto )
    {
      fcFilterParams->filterOutputVetoParams->rsqvetoTimeThresh = rsqVetoTimeThresh;
      fcFilterParams->filterOutputVetoParams->rsqvetoMaxSNR = rsqVetoMaxSNR;
      fcFilterParams->filterOutputVetoParams->rsqvetoCoeff = rsqVetoCoeff;
      fcFilterParams->filterOutputVetoParams->rsqvetoPow = rsqVetoPow;
    }
  }
  else
  {
    fcFilterParams->filterOutputVetoParams = NULL;
  }

  LAL_CALL( LALFindChirpChisqVetoInit( &status, fcFilterParams->chisqParams,
        fcInitParams->numChisqBins, fcInitParams->numPoints ),
      &status );

  /* initialize findchirp clustering params */
  fcFilterParams->clusterMethod = clusterMethod;
  fcFilterParams->clusterWindow = clusterWindow;

  /* parse the thresholds */
  fcFilterParams->rhosqThresh = snrThresh * snrThresh;
  fcFilterParams->chisqThresh = chisqThresh;
  fcFilterParams->chisqDelta  = chisqDelta;

  if ( vrbflg ) fprintf( stdout, "done\n" );


  /*
   *
   * matched filtering engine
   *
   */


  /* begin loop over number of template bank simulation */
  /* if we are not doing a template bank simulation,    */
  /* this exectues the main part of the filtering code  */
  /* just one (which is what we want to do)             */
  bankSimCount = 0;
  do
  {
    if ( bankSim )
    {
      if ( vrbflg )
        fprintf( stdout, "bank simulation %d/%d\n", bankSimCount, bankSim );

      /* zero out the input data segment and the injection params */
      if ( vrbflg ) fprintf( stdout,
          "zeroing data stream and adding random injection for bank sim... " );
      memset( chan.data->data, 0, chan.data->length * sizeof(REAL4) );

      /* inject a random signal if we are doing a template bank simulation */
      if ( ! siminspiral.simInspiralTable )
        thisSimInspiral = siminspiral.simInspiralTable =
          XLALFindChirpBankSimInjectSignal( dataSegVec, &resp, thisBankSim,
              bankSimParamsPtr );
      else
        thisSimInspiral = thisSimInspiral->next =
          XLALFindChirpBankSimInjectSignal( dataSegVec, &resp, thisBankSim,
              bankSimParamsPtr );

      /* write the channel data plus injection to the output frame file */
      if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
          &chan, "ct", "SIM" );
      if ( vrbflg ) fprintf( stdout, "done\n" );
    }


    /*
     *
     * condition data segments for filtering
     *
     */


    switch ( approximant )
    {
      case TaylorT1:
      case TaylorT2:
      case TaylorT3:
      case GeneratePPN:
      case PadeT1:
      case EOB:
      case EOBNR:
      case FindChirpPTF:
        if ( vrbflg )
          fprintf( stdout, "findchirp conditioning data for TD or PTF\n" );
        LAL_CALL( LALFindChirpTDData( &status, fcSegVec, dataSegVec,
              fcDataParams ), &status );
        break;

      case FindChirpSP:
        if ( vrbflg )
          fprintf( stdout, "findchirp conditioning data for SP\n" );
        LAL_CALL( LALFindChirpSPData( &status, fcSegVec, dataSegVec,
              fcDataParams ), &status );
        break;

      case BCV:
        if ( vrbflg )
          fprintf( stdout, "findchirp conditioning data for BCV\n" );
        LAL_CALL( LALFindChirpBCVData( &status, fcSegVec, dataSegVec,
              fcDataParams ), &status );
        break;

      case BCVSpin:
        if ( vrbflg )
          fprintf( stdout, "findchirp conditioning data for BCVSpin\n" );
        LAL_CALL( LALFindChirpTDData( &status, fcSegVec, dataSegVec,
              fcDataParams ), &status );
        break;
      default:
        fprintf( stderr, "error: unknown waveform approximant for data\n" );
        exit( 1 );
        break;
    }

    if ( bankSim )
    {
      if ( vrbflg ) fprintf( stdout,
          "computing minimal match normalization... " );

      matchNorm = XLALFindChirpBankSimSignalNorm( fcDataParams, fcSegVec,
          bankSimCutLowIndex );

      if ( vrbflg ) fprintf( stdout, "%e\n", matchNorm );
    }
    else
    {
      if ( approximant == FindChirpSP )
      {
        /* compute the standard candle */
        REAL4 cannonDist = 1.0; /* Mpc */
        REAL4 m  = (REAL4) candle.tmplt.totalMass;
        REAL4 mu = (REAL4) candle.tmplt.mu;
        REAL4 distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1e6 * LAL_PC_SI);
        REAL4 candleTmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
          pow( m / (LAL_PI*LAL_PI) , 1.0/3.0 ) *
          pow( LAL_MTSUN_SI / (REAL4) chan.deltaT, -1.0/6.0 );

        distNorm *= fcTmpltParams->dynRange;
        candleTmpltNorm *= candleTmpltNorm;
        candleTmpltNorm *= distNorm * distNorm;

        candle.sigmasq = 4.0 * ( (REAL4) chan.deltaT / (REAL4) numPoints );
        candle.sigmasq *= candleTmpltNorm *
          fcSegVec->data->segNorm->data[fcSegVec->data->segNorm->length-1];

        candle.distance = sqrt( candle.sigmasq / candle.rhosq );

        if ( vrbflg )
        {
          fprintf( stdout, "candle m = %e\ncandle mu = %e\n"
              "candle.rhosq = %e\nchan.deltaT = %e\nnumPoints = %d\n"
              "fcSegVec->data->segNorm->data[fcSegVec->data->segNorm->length-1]"
              " = %e\ncandleTmpltNorm = %e\ncandle.distance = %e Mpc\n"
              "candle.sigmasq = %e\n",
              m, mu, candle.rhosq, chan.deltaT, numPoints,
              fcSegVec->data->segNorm->data[fcSegVec->data->segNorm->length-1],
              candleTmpltNorm, candle.distance, candle.sigmasq );
          fflush( stdout );
        }
      }
      else if ( approximant == BCV )
      {
        /* compute the standard candle for a  5-5 Msun inspiral*/
        REAL4 cannonDist = 1.0; /* Mpc */
        REAL4 m  = 10.0;
        REAL4 mu = 2.5;
        REAL4 k1 = numPoints * chan.deltaT ;
        UINT4 kmax = 432 * k1 ; /* 432 = fISCO for 5-5 */
        REAL4 distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1e6 * LAL_PC_SI);
        REAL4 candleTmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
          pow( m/(LAL_PI*LAL_PI) , 1.0/3.0 ) *
          pow(LAL_MTSUN_SI / (REAL4) chan.deltaT, -1.0/6.0);
        distNorm *= fcTmpltParams->dynRange;
        candleTmpltNorm *= candleTmpltNorm;
        candleTmpltNorm *= distNorm * distNorm;
        candle.sigmasq = 4.0 * ( (REAL4) chan.deltaT / (REAL4) numPoints );
        candle.sigmasq *= candleTmpltNorm *
          fcSegVec->data->segNorm->data[kmax];
        candle.distance = sqrt( candle.sigmasq / candle.rhosq );

        /* for 5-5 Msun... */
        candle.tmplt.mass1 = 5.0;
        candle.tmplt.mass2 = 5.0;
        candle.tmplt.totalMass = 10.0;
        candle.tmplt.mu = 2.5;
        candle.tmplt.eta = 0.25;

        if ( vrbflg )
        {
          fprintf( stdout, "candle m = %e\ncandle mu = %e\n"
              "candle.rhosq = %e\nchan.deltaT = %e\nnumPoints = %d\n"
              "fcSegVec->data->segNorm->data[kmax] = %e\n"
              "kmax = %d\ncandleTmpltNorm = %e\ncandle.distance = %e Mpc \n"
              "candle.sigmasq=%e\n",
              m,mu,candle.rhosq,chan.deltaT,
              numPoints,fcSegVec->data->segNorm->data[kmax],kmax,
              candleTmpltNorm,candle.distance,candle.sigmasq);
          fflush(stdout);
        }
      }
      else
      {
        if ( vrbflg )
        {
          fprintf( stdout, "standard candle not calculated;\n"
              "chan.deltaT = %e\nnumPoints = %d\n",
              chan.deltaT, numPoints );
          fflush( stdout );
        }
      }
    }


    /*================================================================
     * This is a dense bit of logical constructs and basically
     * implements the following truth table
     *
     *   Injection   Follow-Up    --fast      Action
     *   ---------   ---------    --------  ----------
     *      T           F           F        SetAnalyseTemplate()
     *      T           F           T              "
     *      T           T           F              "
     *      T           T           T        TagTemplateAndSegment()
     *      F           T           T              "
     *      F           T           F        Do Nothing
     *      F           F           F        Do Nothing
     *      F           F           T        Irrelevant (do nothing)
     *================================================================*/
    if ( numTmplts > 0 &&
        ( injectionFile  || (tdFollowUpFiles && mmFast >= 0.0) ) )
    {
      /* Make space for analyseThisTmplt */
      analyseThisTmplt = (UINT4 *) LALCalloc (numTmplts, sizeof(UINT4));

      /* If we are performing a follow up, we only want to follow up templates
       * near coincident events, regardless of the injection file used
       */
      if ( injectionFile && (!tdFollowUpFiles || mmFast < 0) )
      {

      /* set the analyseThisTmplt flag on templates     */
      /* that are 'close' to the injections             */
        LAL_CALL( LALFindChirpSetAnalyseTemplate( &status, analyseThisTmplt,
            mmFast, fcSegVec->data[0].data->deltaF, sampleRate, fcDataParams,
            numTmplts, bankHead, numInjections, injections ), &status );

      }
      else
      {
        XLALFindChirpTagTemplateAndSegment(dataSegVec, bankHead,
            &tdFollowUpEvents, ifo, mmFast, analyseThisTmplt );
      }
    }

    /*
     *
     * split the template bank into subbanks for the bank veto
     *
     */

    /* only sort the bank if doing bank veto */

    if (!bankSimCount && numTmplts > 0) /*just doing this once is fine*/
    {
      if (subBankSize > 1)
         bankHead = XLALFindChirpSortTemplates( bankHead, numTmplts, subBankSize);

      if ( vrbflg ) fprintf( stdout,
        "splitting bank in to subbanks of size ~ %d\n", subBankSize );

      subBankHead = XLALFindChirpCreateSubBanks( &maxSubBankSize,
        subBankSize, numTmplts, bankHead );

      if ( vrbflg ) fprintf( stdout, "maximum subbank size = %d\n",
        maxSubBankSize );

      bankVetoData.length = maxSubBankSize;

      /* free the existing qVec structure so that we can use the bankveto one */
      XLALDestroyCOMPLEX8Vector( fcFilterParams->qtildeVec );
      fcFilterParams->qtildeVec = NULL;
      XLALDestroyCOMPLEX8Vector( fcFilterParams->qVec );
      fcFilterParams->qVec = NULL;


      /* create an array of fcFilterInput structs */
      bankVetoData.qVecArray = (COMPLEX8Vector **)
        LALCalloc( bankVetoData.length, sizeof(COMPLEX8Vector*) );
      bankVetoData.qtildeVecArray = (COMPLEX8Vector **)
        LALCalloc( bankVetoData.length, sizeof(COMPLEX8Vector*) );
      bankVetoData.fcInputArray = (FindChirpFilterInput **)
        LALCalloc( bankVetoData.length, sizeof(FindChirpFilterInput*) );
      /* create ccMat for bank veto */
      bankVetoData.ccMat =
        XLALCreateVector( bankVetoData.length * bankVetoData.length );
      bankVetoData.normMat =
        XLALCreateVector( bankVetoData.length * bankVetoData.length );
      /* point to response and spectrum */
      bankVetoData.spec = spec.data;
      bankVetoData.resp = resp.data;

      for ( i = 0; i < bankVetoData.length; ++i )
      {
        bankVetoData.qVecArray[i] =
          XLALCreateCOMPLEX8Vector( fcInitParams->numPoints );
        bankVetoData.qtildeVecArray[i] =
          XLALCreateCOMPLEX8Vector( fcInitParams->numPoints );
        LAL_CALL( LALCreateFindChirpInput( &status,
              &(bankVetoData.fcInputArray[i]), fcInitParams ), &status );
      }
    }

    /*
     *
     * search engine
     *
     */

    for ( subBankCurrent = subBankHead, thisTemplateIndex = 0;
        subBankCurrent;
        subBankCurrent = subBankCurrent->next, thisTemplateIndex++ )
    {

      /* If we are injecting or in td-follow-up mode and the    */
      /* analyseThisTmpltFlag is down look no further:          */
      /* simply continue to the next template                   */
      if ( injectionFile )
      {
        if ( ! analyseThisTmplt[thisTemplateIndex] )
          continue;
        else
          if ( vrbflg ) fprintf( stdout,
              "\n\n === Template %d going through \n\n", thisTemplateIndex );
      }

      for ( bankCurrent = subBankCurrent->bankHead, subBankIndex = 0;
          bankCurrent; bankCurrent = bankCurrent->next, ++subBankIndex )
      {
        /* set fcFilterInput to the correct one */
        fcFilterInput = bankVetoData.fcInputArray[subBankIndex];
        if ( vrbflg ) fprintf( stdout,
            "Creating template in fcInputArray[%d] at %p\n", subBankIndex,
            fcFilterInput );

        /*  generate template */
        switch ( approximant )
        {
          case TaylorT1:
          case TaylorT2:
          case TaylorT3:
          case GeneratePPN:
          case PadeT1:
          case EOB:
          case EOBNR:
            LAL_CALL( LALFindChirpTDTemplate( &status, fcFilterInput->fcTmplt,
                  bankCurrent, fcTmpltParams ), &status );
            break;

          case FindChirpSP:
            LAL_CALL( LALFindChirpSPTemplate( &status, fcFilterInput->fcTmplt,
                  bankCurrent, fcTmpltParams ), &status );
            break;

          case BCV:
            LAL_CALL( LALFindChirpBCVTemplate( &status, fcFilterInput->fcTmplt,
                  bankCurrent, fcTmpltParams ), &status );
            break;

          case BCVSpin:
            LAL_CALL( LALFindChirpBCVSpinTemplate( &status,
                  fcFilterInput->fcTmplt, bankCurrent,
                  fcTmpltParams, fcDataParams ), &status );
            break;

          case FindChirpPTF:
            LAL_CALL( LALFindChirpPTFTemplate( &status, fcFilterInput->fcTmplt,
                  bankCurrent, fcTmpltParams ), &status );
            break;

          default:
            fprintf( stderr,
                "error: unknown waveform template approximant \n" );
            exit( 1 );
            break;
        }
      }

      /* Write out the template time series if requested */
      if (writeTemplate )
      {
        /* Locally scoped variables */
        CHAR snrsqStr[LALNameLength];
        COMPLEX8Vector *templateFFTDataVector = NULL;
        REAL4FFTPlan *plan = NULL;
        REAL8 deltaF;
        INT4 kmax, numPoints, nb2;
        LALSnprintf( snrsqStr, LALNameLength*sizeof(CHAR),
                  "TEMPLATE");
        memcpy(&templateTimeSeries, &chan, sizeof(REAL4TimeSeries));
        strcpy( templateTimeSeries.name, chan.name );
        /* Things are complicated if the template is in the Frequency domain */
        if ( approximant==FindChirpSP )
          {
          nb2 = fcTmpltParams->xfacVec->length;
          numPoints = (int) floor( (nb2-1) * 2.0);
          templateTimeSeriesVector = XLALCreateREAL4Vector(numPoints);
          templateFFTDataVector = XLALCreateCOMPLEX8Vector(nb2);
          numPoints = templateTimeSeriesVector->length;

          deltaF = 1.0 / ( fcTmpltParams->deltaT * (REAL8) numPoints / 2.0);
          kmax = fcFilterInput->fcTmplt->tmplt.fFinal / deltaF < numPoints/2 ?
            fcFilterInput->fcTmplt->tmplt.fFinal / deltaF : numPoints/2;

          for (i = 0; i < nb2; i++)
            {
            if (1 || (i * deltaF) > fLow && i < kmax)
              {
              templateFFTDataVector->data[i].re = fcFilterInput->fcTmplt->data->data[i].re * fcTmpltParams->xfacVec->data[i];
              templateFFTDataVector->data[i].im = fcFilterInput->fcTmplt->data->data[i].im * fcTmpltParams->xfacVec->data[i];
              }
            else
              {
              templateFFTDataVector->data[i].re = 0;
              templateFFTDataVector->data[i].im = 0;
              }
            }
          plan = XLALCreateReverseREAL4FFTPlan( numPoints, 0);
          if ( XLALREAL4ReverseFFT( templateTimeSeriesVector, templateFFTDataVector, plan) )  fprintf(stderr, "\n\nFFT FAILED\n\n");
          templateTimeSeries.data = templateTimeSeriesVector;
          }
        /* Much easier for a TD template */
        else
          templateTimeSeries.data = fcTmpltParams->xfacVec;

        if (templateTimeSeries.data)
               outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
                  &templateTimeSeries, "none", snrsqStr );
        else fprintf(stderr,"fcTmpltParams->xfacVec is NULL\n");
        /* if we had to make intermediate data products, destroy them */
        /* The template time series persists and will be destroyed later */
        if (templateFFTDataVector) XLALDestroyCOMPLEX8Vector(templateFFTDataVector);
        if (plan) XLALDestroyREAL4FFTPlan( plan );
      }

      ccFlag = 1;
      /* loop over data segments */
      for ( i = 0; i < fcSegVec->length ; ++i )
      {
        INT8 fcSegStartTimeNS;
        INT8 fcSegEndTimeNS;

        fcSegStartTimeNS = XLALGPSToINT8NS( &(fcSegVec->data[i].data->epoch) );
        fcSegEndTimeNS = fcSegStartTimeNS + (INT8)
          ( (REAL8) numPoints * 1e9 * fcSegVec->data[i].deltaT );

        /* skip segment if it is not contained in the trig start or end times */
        if ( ! bankSim &&
            ( (trigStartTimeNS && (trigStartTimeNS > fcSegEndTimeNS)) ||
            (trigEndTimeNS && (trigEndTimeNS < fcSegStartTimeNS)) ) )
        {
          if ( vrbflg ) fprintf( stdout,
              "skipping segment %d/%d [%ld-%ld] (outside trig time)\n",
              fcSegVec->data[i].number, fcSegVec->length,
              fcSegStartTimeNS, fcSegEndTimeNS );

          continue;
        }

        /* By default assume that we are going to analyse the   */
        /* segment with the template (analyseTag = 1)           */
        analyseTag = 1;

        /* If injections are being done or if in td follow-up mode */
        /* check if for any reason the analyseTag flag needs to be */
        /* brought down                                            */
        if ( tdFollowUpFiles || ( injectionFile  && flagFilterInjOnly ))
        {
          if ( tdFollowUpFiles )
          {
            if ( mmFast >= 0.0 )
            {
                /* Compare tmplt flag with segment number using bitwise and */
                if ( !(analyseThisTmplt[thisTemplateIndex] & ( 1u << i )))
                {
                    analyseTag = 0;
                }
            }
            else if ( !fcSegVec->data[i].analyzeSegment )
            {
              analyseTag = 0;
            }

          }
          else
          {
            if (subBankSize > 1)
               analyseTag = fcSegVec->data[i].analyzeSegment;
            else
               analyseTag = XLALCmprSgmntTmpltFlags( numInjections,
                  analyseThisTmplt[thisTemplateIndex],
                  fcSegVec->data[i].analyzeSegment );
          }
        }

        for ( bankCurrent = subBankCurrent->bankHead, subBankIndex = 0;
            bankCurrent; bankCurrent = bankCurrent->next, ++subBankIndex )
        {
          /* filter data segment */
          if ( analyseTag )
          {
            /* set fcFilterInput and qVec to the correct ones */
            fcFilterInput = bankVetoData.fcInputArray[subBankIndex];
            if ( vrbflg ) fprintf( stdout,
                "Using template in fcInputArray[%d] at %p\n", subBankIndex,
                fcFilterInput );
            fcFilterParams->qVec=bankVetoData.qVecArray[subBankIndex];
            fcFilterParams->qtildeVec=bankVetoData.qtildeVecArray[subBankIndex];
            if ( vrbflg ) fprintf( stdout,
                "Using qVec in qVecArray[%d] at %p\n", subBankIndex,
                fcFilterParams->qVec );
            if ( vrbflg ) fprintf( stdout,
                "filtering segment %d/%d [%ld-%ld] "
                "against template %d/%d (%e,%e)\n",
                fcSegVec->data[i].number, fcSegVec->length,
                fcSegStartTimeNS, fcSegEndTimeNS,
                bankCurrent->number, numTmplts,
                fcFilterInput->fcTmplt->tmplt.mass1,
                fcFilterInput->fcTmplt->tmplt.mass2 );

            fcFilterInput->segment = fcSegVec->data + i;


            /* decide which filtering routine to use */
            switch ( approximant )
            {
              case TaylorT1:
              case TaylorT2:
              case TaylorT3:
              case GeneratePPN:
              case PadeT1:
              case EOB:
              case EOBNR:
                /* construct normalization for time domain templates... */
                LAL_CALL( LALFindChirpTDNormalize( &status,
                      fcFilterInput->fcTmplt, fcFilterInput->segment,
                      fcDataParams ), &status );

                /* ...and fall through to FindChirpFilterSegment() */
              case FindChirpSP:
                LAL_CALL( LALFindChirpFilterSegment( &status,
                      &eventList, fcFilterInput, fcFilterParams ), &status );
                break;

              case BCV:
                if ( ! bcvConstraint )
                {
                  LAL_CALL( LALFindChirpBCVFilterSegment( &status,
                        &eventList, fcFilterInput, fcFilterParams ), &status );
                }
                else
                {
                  LAL_CALL( LALFindChirpBCVCFilterSegment( &status,
                        &eventList, fcFilterInput, fcFilterParams ), &status );
                }
                break;

              case BCVSpin:
                LAL_CALL( LALFindChirpBCVSpinFilterSegment( &status,
                      &eventList, fcFilterInput, fcFilterParams, fcDataParams
                      ), &status );
                break;

              case FindChirpPTF:
                LAL_CALL( LALFindChirpPTFNormalize( &status,
                      fcFilterInput->fcTmplt, fcFilterInput->segment,
                      fcDataParams ), &status );
                LAL_CALL( LALFindChirpPTFFilterSegment( &status,
                      &eventList, fcFilterInput, fcFilterParams ), &status );
                break;

              default:
                fprintf( stderr,
                    "error: unknown waveform approximant for filter\n" );
                exit( 1 );
                break;
            }

            if ( writeRhosq )
            {
              CHAR snrsqStr[LALNameLength];
              LALSnprintf( snrsqStr, LALNameLength*sizeof(CHAR),
                  "SNRSQ_%d", nRhosqFr++ );
              strcpy( fcFilterParams->rhosqVec->name, chan.name );
              outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
                  fcFilterParams->rhosqVec, "none", snrsqStr );
            }

            if ( vrbflg )
              fprintf( stdout, "epoch = %d\n",
                  fcFilterInput->segment->data->epoch.gpsSeconds );
          }
          else /* not analyzeTag */
          {
            if ( vrbflg ) fprintf( stdout,
                "skipping segment %d/%d [%ld-%ld]\n",
                fcSegVec->data[i].number, fcSegVec->length,
                fcSegStartTimeNS, fcSegEndTimeNS );
          }


          /*  test if filter returned any events */
          if ( eventList )
          {
            /* this can only happen for BCV and BCVSpin */
            if ( vrbflg ) fprintf( stdout,
                "segment %d rang template [m (%e,%e)] [psi (%e,%e)]\n",
                fcSegVec->data[i].number,
                fcFilterInput->fcTmplt->tmplt.mass1,
                fcFilterInput->fcTmplt->tmplt.mass2,
                fcFilterInput->fcTmplt->tmplt.psi0,
                fcFilterInput->fcTmplt->tmplt.psi3 );

            if ( vrbflg ) fprintf( stdout, "***>  dumping events  <***\n" );

            if ( ! savedEvents.snglInspiralTable )
            {
              savedEvents.snglInspiralTable = eventList;
            }
            else
            {
              event->next = eventList;
            }

            /* save a ptr to the last event in the list and count the events */
            ++numEvents;
            while ( eventList->next )
            {
              eventList = eventList->next;
              ++numEvents;
            }
            event = eventList;
            eventList = NULL;

          } /* end if ( events ) for BCV and BCVSpin */
        } /* end of loop over templates in subbank */

        /* If doing bank veto compute CC Matrix */
        /* I removed the ccFlag dependence - this is being computed
           for each segment now!!! */
        if (ccFlag && (subBankCurrent->subBankSize > 1) && analyseTag)
        {
          if (vrbflg) fprintf(stderr, "doing ccmat\n");
          XLALBankVetoCCMat( &bankVetoData, subBankCurrent, fcDataParams,
          dynRange, fLow, spec.deltaF, chan.deltaT);
          ccFlag = 0;
          /*char filename[10];
          sprintf(filename, "ccmat%d.dat",i);
          FILE *FP = NULL;
          FP = fopen(filename,"w");
          for(j = 0; j < bankVetoData.ccMat->length; ++j)
          {
            fprintf(FP, "%e\n",bankVetoData.ccMat->data[j]);
          }*/
        }
        /* now look through the filter outputs of the subbank for events */
        for ( bankCurrent = subBankCurrent->bankHead, subBankIndex = 0;
            bankCurrent; bankCurrent = bankCurrent->next, ++subBankIndex )
        {
          if ( analyseTag )
          {
            /* set fcFilterInput and qVec to the correct ones */
            fcFilterInput = bankVetoData.fcInputArray[subBankIndex];
            if ( vrbflg ) fprintf( stdout,
                "Finding Events - Using template in fcInputArray[%d] at %p\n",
                subBankIndex, fcFilterInput );
            fcFilterParams->qtildeVec=bankVetoData.qtildeVecArray[subBankIndex];
            fcFilterParams->qVec = bankVetoData.qVecArray[subBankIndex];
            if ( vrbflg ) fprintf( stdout,
                "Finding Events - Using qVec in qVecArray[%d] at %p\n",
                subBankIndex, fcFilterParams->qVec );

            /* determine if FindChirpFilterSegment returned any events */
            switch ( approximant )
            {
              case TaylorT1:
              case TaylorT2:
              case TaylorT3:
              case GeneratePPN:
              case PadeT1:
              case EOB:
              case EOBNR:
                /* recompute the template norm since it has been over written */
                /* ( When doing the chisq test and the bank veto for          */
                /* time domain searches that don't use the                    */
                /* the template power vec that goes as f^-(7/3) )             */
                if (subBankCurrent->subBankSize > 1)
                {
                LAL_CALL( LALFindChirpTDNormalize( &status,
                      fcFilterInput->fcTmplt, fcFilterInput->segment,
                      fcDataParams ), &status );
                }
              case FindChirpSP:
                /* find any events in the time series of snr and chisq */
                LAL_CALL( LALFindChirpClusterEvents( &status,
                      &eventList, fcFilterInput, fcFilterParams,
                      &bankVetoData, subBankIndex ), &status );

                if ( writeChisq )
                {
                  CHAR chisqStr[LALNameLength];
                  REAL4TimeSeries chisqts;
                  LALSnprintf( chisqStr, LALNameLength*sizeof(CHAR),
                      "CHISQ_%d", nChisqFr++ );
                  chisqts.epoch = fcFilterInput->segment->data->epoch;
                  memcpy( &(chisqts.name), fcFilterInput->segment->data->name,
                      LALNameLength * sizeof(CHAR) );
                  chisqts.deltaT = fcFilterInput->segment->deltaT;
                  chisqts.data = fcFilterParams->chisqVec;
                  outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
                      &chisqts, "none", chisqStr );
                }

                if ( writeCData && ! strcmp(ifo, bankCurrent->ifo) )
                {
                  trigTime = bankCurrent->end_time.gpsSeconds + 1e-9 *
                    bankCurrent->end_time.gpsNanoSeconds;
                  lowerBound =
                    gpsStartTime.gpsSeconds + numPoints/(4 * sampleRate );
                  upperBound =
                    gpsEndTime.gpsSeconds - numPoints/(4 * sampleRate );

                  if ( trigTime >= lowerBound && trigTime <= upperBound )
                  {
                    tempTmplt = (SnglInspiralTable *)
                      LALCalloc(1, sizeof(SnglInspiralTable) );
                    tempTmplt->event_id = (EventIDColumn *)
                      LALCalloc(1, sizeof(EventIDColumn) );
                    tempTmplt->mass1 = bankCurrent->mass1;
                    tempTmplt->end_time.gpsSeconds =
                      bankCurrent->end_time.gpsSeconds;
                    tempTmplt->end_time.gpsNanoSeconds =
                      bankCurrent->end_time.gpsNanoSeconds;
                    tempTmplt->event_id->id = bankCurrent->event_id->id;

                    tempTmplt->sigmasq = eventList->sigmasq;

                    LAL_CALL( LALFindChirpCreateCoherentInput( &status,
                          &coherentInputData, fcFilterParams->cVec,
                          tempTmplt, 0.5, numPoints / 4 ), &status );

                    LALFree( tempTmplt->event_id );
                    LALFree( tempTmplt );

                    if ( coherentInputData )
                    {
                      cDataForFrame = 1;
                      LALSnprintf( cdataStr, LALNameLength*sizeof(CHAR),
                                   "%Ld", bankCurrent->event_id->id );
                      LALSnprintf( coherentInputData->name,
                                   LALNameLength*sizeof(CHAR),
                                   "%s:CBC-CData", ifo );
                      if ( ! coherentFrames )
                      {
                        thisCoherentFrame = coherentFrames = (FrameHNode *)
                          LALCalloc( 1, sizeof(FrameHNode) );
                      }
                      else
                      {
                        thisCoherentFrame = thisCoherentFrame->next =
                          (FrameHNode *) LALCalloc( 1, sizeof(FrameHNode) );
                      }
                      thisCoherentFrame->frHeader = fr_add_proc_COMPLEX8TimeSeries(
                          outFrame, coherentInputData, "none", cdataStr );
                      LAL_CALL( LALCDestroyVector( &status,
                            &(coherentInputData->data) ), &status );
                      LALFree( coherentInputData );
                      coherentInputData = NULL;
                    }
                  }
                }

                /* apply the rsq veto to any surviving events */
                if ( fcFilterParams->filterOutputVetoParams )
                {
                  LAL_CALL( LALFindChirpFilterOutputVeto( &status,
                        &eventList, fcFilterInput, fcFilterParams ), &status );
                }
                break;
              default:
                break;
            }

            /*  test if filter returned any events */
            if ( eventList )
            {
              /* this can only happen for FindChirpFilter */

              if ( vrbflg )
              {
                if( approximant == FindChirpPTF)
                {
                  fprintf( stdout,
                      "segment %d rang template [m (%e,%e)], chi = %e,"
                      "kappa = %e)]\n",
                      fcSegVec->data[i].number,
                      fcFilterInput->fcTmplt->tmplt.mass1,
                      fcFilterInput->fcTmplt->tmplt.mass2,
                      fcFilterInput->fcTmplt->tmplt.chi,
                      fcFilterInput->fcTmplt->tmplt.kappa );
                }
                else
                  fprintf( stdout,
                      "segment %d rang template [m (%e,%e)] [psi (%e,%e)]\n",
                      fcSegVec->data[i].number,
                      fcFilterInput->fcTmplt->tmplt.mass1,
                      fcFilterInput->fcTmplt->tmplt.mass2,
                      fcFilterInput->fcTmplt->tmplt.psi0,
                      fcFilterInput->fcTmplt->tmplt.psi3 );
              }

              if ( vrbflg ) fprintf( stdout, "***>  dumping events  <***\n" );

              if ( ! savedEvents.snglInspiralTable )
              {
                savedEvents.snglInspiralTable = eventList;
              }
              else
              {
                event->next = eventList;
              }

              /* save a ptr to the last event in the list , count the events */
              ++numEvents;
              while ( eventList->next )
              {
                eventList = eventList->next;
                ++numEvents;
              }
              event = eventList;
              eventList = NULL;

            } /* end if ( events ) for FindChirpFilter */

          } /* end if (analyseTag) */

        } /* end loop over template in subbank */

      } /* end loop over data segments */

      /* delete the chisq bins for templates */
      switch ( approximant )
      {
        case TaylorT1:
        case TaylorT2:
        case TaylorT3:
        case GeneratePPN:
        case PadeT1:
        case EOB:
        case EOBNR:
        case FindChirpSP:
          /* the chisq bins need to be re-computed for the next template */
          for ( i = 0; i < fcSegVec->length ; ++i )
          {
            if ( fcSegVec->data[i].chisqBinVec->data )
            {
              LALFree( fcSegVec->data[i].chisqBinVec->data );
              fcSegVec->data[i].chisqBinVec->data = NULL;
            }
          }
        default:
          break;
      }

    } /* end loop over linked list */

    if ( bankSim  )
    {
      if ( bankSimParams.maxMatch )
      {
        /* compute the best match over the template bank */
        if ( ! siminstparams.simInstParamsTable )
          thisSimInstParams = siminstparams.simInstParamsTable =
            XLALFindChirpBankSimMaxMatch( &(savedEvents.snglInspiralTable),
                matchNorm );
        else
          thisSimInstParams = thisSimInstParams->next =
            XLALFindChirpBankSimMaxMatch( &(savedEvents.snglInspiralTable),
                matchNorm );

        /* save only the template that gives the loudest event */
        if ( ! loudestEventHead )
          thisLoudestEvent = loudestEventHead =
            savedEvents.snglInspiralTable;
        else
          thisLoudestEvent = thisLoudestEvent->next =
            savedEvents.snglInspiralTable;
      }
      else
      {
        /* compute the match for all the templates in the bank */
        for ( event = savedEvents.snglInspiralTable; event;
            event = event->next )
        {
          if ( ! siminstparams.simInstParamsTable )
            thisSimInstParams = siminstparams.simInstParamsTable =
              XLALFindChirpBankSimComputeMatch( event, matchNorm );
          else
            thisSimInstParams = thisSimInstParams->next =
              XLALFindChirpBankSimComputeMatch( event, matchNorm );
        }

        /* append all the inspiral templates to the loudest event list */
        if ( ! loudestEventHead )
          thisLoudestEvent = loudestEventHead = savedEvents.snglInspiralTable;
        else
          thisLoudestEvent = thisLoudestEvent->next =
            savedEvents.snglInspiralTable;

        while ( thisLoudestEvent && thisLoudestEvent->next )
          thisLoudestEvent = thisLoudestEvent->next;
      }

      /* null out the list of inspiral events */
      savedEvents.snglInspiralTable = NULL;
    }

    /* increment the bank sim sim_inspiral table if necessary */
    if ( bankSimHead )
    {
      thisBankSim = thisBankSim->next;
    }

  } while ( ++bankSimCount < bankSim ); /* end loop over bank simulations */

  if ( bankSim )
  {
    /* save the number of bank simulation in the search summary table */
    searchsumm.searchSummaryTable->nevents = bankSim;

    /* point the saved events to the linked list of loudest events */
    savedEvents.snglInspiralTable = loudestEventHead;
  }
  else
  {
    /* save the number of events in the search summary table */
    searchsumm.searchSummaryTable->nevents = numEvents;
  }


  /*
   *
   * free memory used by filtering code
   *
   */


  if ( vrbflg ) fprintf( stdout, "freeing memory\n" );

  /* free memory used by findchirp */
  LAL_CALL( LALFindChirpChisqVetoFinalize( &status,
        fcFilterParams->chisqParams, fcInitParams->numChisqBins ),
      &status );

  if ( numTmplts > 0 )
  {
    for ( i = 0; i < bankVetoData.length; ++i )
    {
      XLALDestroyCOMPLEX8Vector( bankVetoData.qVecArray[i] );
      XLALDestroyCOMPLEX8Vector( bankVetoData.qtildeVecArray[i] );
      bankVetoData.qVecArray[i] = NULL;
      bankVetoData.qtildeVecArray[i] = NULL;
      LAL_CALL( LALDestroyFindChirpInput( &status,
            &(bankVetoData.fcInputArray[i]) ), &status );
      bankVetoData.fcInputArray[i] = NULL;
    }
    LALFree( bankVetoData.qVecArray );
    LALFree( bankVetoData.qtildeVecArray );
    LALFree( bankVetoData.fcInputArray );
    XLALDestroyVector( bankVetoData.ccMat );
    XLALDestroyVector( bankVetoData.normMat );
    /* XLALDestroyVector( bankVetoData.normMat ); */
    fcFilterParams->qVec = NULL;
    fcFilterParams->qtildeVec = NULL;
  }


  if ( fcFilterParams->filterOutputVetoParams )
  {
    LALFree( fcFilterParams->filterOutputVetoParams );
  }
  LAL_CALL( LALFindChirpFilterFinalize( &status, &fcFilterParams ),
      &status );
  LAL_CALL( LALFindChirpTemplateFinalize( &status, &fcTmpltParams ),
      &status );
  LAL_CALL( LALFindChirpDataFinalize( &status, &fcDataParams ),
      &status );
  LAL_CALL( LALDestroyFindChirpSegmentVector( &status, &fcSegVec ),
      &status );
  LALFree( fcInitParams );

  /* free the template bank */
  if ( subBankHead )
  {
    while ( subBankHead )
    {
      subBankCurrent = subBankHead;
      while ( subBankCurrent->bankHead )
      {
        bankCurrent = subBankCurrent->bankHead;
        subBankCurrent->bankHead = subBankCurrent->bankHead->next;
        if ( bankCurrent->event_id )
        {
          LALFree( bankCurrent->event_id );
        }
        LALFree( bankCurrent );
        bankCurrent = NULL;
      }
      subBankHead = subBankHead->next;
      LALFree( subBankCurrent );
      subBankCurrent = NULL;
    }
  }

  /* free any bank sim parameters */
  if ( bankSimHead )
  {
    while ( bankSimHead )
    {
      thisBankSim = bankSimHead;
      bankSimHead = bankSimHead->next;
      LALFree( thisBankSim );
    }
  }

  /* free the data storage */
  LAL_CALL( LALFinalizeDataSegmentVector( &status, &dataSegVec ), &status );
  LAL_CALL( LALSDestroyVector( &status, &(chan.data) ), &status );
  LAL_CALL( LALSDestroyVector( &status, &(spec.data) ), &status );
  LAL_CALL( LALCDestroyVector( &status, &(resp.data) ), &status );

  /* free the random parameters structure */
  if ( randSeedType != unset )
  {
    LAL_CALL( LALDestroyRandomParams( &status, &randParams ), &status );
  }


  /*
   *
   * write the result results to disk
   *
   */


  /* write the output frame */
  if ( writeRawData || writeFilterData || writeResponse || writeSpectrum ||
       writeRhosq || writeChisq || writeCData )
  {
    if ( outputPath[0] )
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s/%s.gwf",
          outputPath, fileName );
    }
    else
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s.gwf", fileName );
    }
    if ( vrbflg ) fprintf( stdout, "writing frame data to %s... ", fname );
    frOutFile = FrFileONew( fname, 0 );
    if ( writeRawData || writeFilterData || writeResponse || writeSpectrum ||
        writeRhosq || writeChisq )
    {
      FrameWrite( outFrame, frOutFile );
    }
    while( coherentFrames )
    {
      thisCoherentFrame = coherentFrames;
      FrameWrite( coherentFrames->frHeader, frOutFile );
      coherentFrames = coherentFrames->next;
      LALFree( thisCoherentFrame );
      thisCoherentFrame = NULL;
    }
    FrFileOEnd( frOutFile );
    if ( vrbflg ) fprintf( stdout, "done\n" );
  }

  /* open the output xml file */
  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  if ( outputPath[0] )
  {
    if ( outCompress )
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s/%s.xml.gz",
          outputPath, fileName );
    }
    else
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s/%s.xml",
          outputPath, fileName );
    }
  }
  else
  {
    if ( outCompress )
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s.xml.gz", fileName );
    }
     else
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s.xml", fileName );
    }
  }
  if ( vrbflg ) fprintf( stdout, "writing XML data to %s...\n", fname );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, fname ), &status );

  /* write the process table */
  if ( vrbflg ) fprintf( stdout, "  process table...\n" );
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable,
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  free( proctable.processTable );

  /* write the process params table */
  if ( vrbflg ) fprintf( stdout, "  process_params table...\n" );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams,
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  /* write the search summary table */
  if ( vrbflg ) fprintf( stdout, "  search_summary table...\n" );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results,
        search_summary_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsumm,
        search_summary_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

  /* write the filter table */
  if ( vrbflg ) fprintf( stdout, "  filter table...\n" );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results,
        filter_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, filtertable,
        filter_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  free( filtertable.filterTable );

  /* write the search summvars table */
  if ( numTmplts )
  {
    if ( vrbflg ) fprintf( stdout, "  search_summvars table...\n" );
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results,
          search_summvars_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsummvars,
          search_summvars_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }
  while( searchsummvars.searchSummvarsTable )
  {
    this_search_summvar = searchsummvars.searchSummvarsTable;
    searchsummvars.searchSummvarsTable = this_search_summvar->next;
    LALFree( this_search_summvar );
  }

  /* write the summ_value table with the standard candle distance */
  if ( ! bankSim )
  {
    if ( approximant == FindChirpSP )
    {
      if ( vrbflg ) fprintf( stdout, "  summ_value table...\n" );
      ADD_SUMM_VALUE( "inspiral_effective_distance", "1.4_1.4_8",
          candle.distance, 0);
    }
    else if ( approximant == BCV )
    {
      if ( vrbflg ) fprintf( stdout, "  summ_value table...\n" );
      ADD_SUMM_VALUE( "inspiral_effective_distance", "5.0_5.0_8",
          candle.distance, 0);
    }
  }

  /* store calibration information */
  ADD_SUMM_VALUE( "calibration alpha", "analysis", alpha, 0 );
  ADD_SUMM_VALUE( "calibration alphabeta", "analysis", alphabeta, 0 );
  if (injectionFile)
  {
    ADD_SUMM_VALUE( "calibration alpha", "injection", inj_alpha, 0 );
    ADD_SUMM_VALUE( "calibration alphabeta", "injection", inj_alphabeta, 0 );
  }
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, summ_value_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, summvalue,
        summ_value_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

  while ( summvalue.summValueTable )
  {
    this_summ_value = summvalue.summValueTable;
    summvalue.summValueTable = summvalue.summValueTable->next;
    LALFree( this_summ_value );
  }

  /* free the search summary table after the summ_value table is written */
  free( searchsumm.searchSummaryTable );

  /* cut triggers based on start/end times and do trig_scan clustering */
  if ( savedEvents.snglInspiralTable )
  {
    SnglInspiralTable *tmpEventHead = NULL;
    SnglInspiralTable *lastEvent = NULL;

    /* sort the inspiral events by time */
    if ( vrbflg ) fprintf( stdout, "  sorting events by time... " );
    LAL_CALL( LALSortSnglInspiral( &status, &(savedEvents.snglInspiralTable),
          LALCompareSnglInspiralByTime), &status );
    if ( vrbflg ) fprintf( stdout, "done\n" );

    /* discard any triggers outside the trig start/end time window */
    event = savedEvents.snglInspiralTable;
    if ( ! bankSim && ( trigStartTimeNS || trigEndTimeNS ) )
    {
      if ( vrbflg ) fprintf( stdout,
          "  discarding triggers outside trig start/end time... " );

      while ( event )
      {
        INT8 trigTimeNS;
        trigTimeNS = XLALGPSToINT8NS( &(event->end_time) );

        if ( trigTimeNS &&
            ((trigStartTimeNS && (trigTimeNS < trigStartTimeNS)) ||
             (trigEndTimeNS && (trigTimeNS >= trigEndTimeNS))) )
        {
          /* throw this trigger away */
          SnglInspiralTable *tmpEvent = event;

          if ( lastEvent )
          {
            lastEvent->next = event->next;
          }

          /* increment the linked list by one and free the event */
          event = event->next;
          LALFree( tmpEvent );
        }
        else
        {
          /* store the first event as the head of the new linked list */
          if ( ! tmpEventHead ) tmpEventHead = event;

          /* save the last event and increment the linked list by one */
          lastEvent = event;
          event = event->next;
        }
      }

      savedEvents.snglInspiralTable = tmpEventHead;

      if ( vrbflg ) fprintf( stdout, "done\n" );
    }

    if (maximizationInterval)
    {
      XLALMaxSnglInspiralOverIntervals( &(savedEvents.snglInspiralTable),
          maximizationInterval);
    }

    /* trigScanClustering */ 
    if ( trigScanMethod ) 
    { 
        if ( savedEvents.snglInspiralTable) 
        { 
            
           /* Call the clustering routine */ 
           if (XLALTrigScanClusterTriggers( &(savedEvents.snglInspiralTable),
                                       trigScanMethod,
                                       trigScanMetricScalingFac,
                                       trigScanAppendStragglers ) == XLAL_FAILURE )
           {
             fprintf( stderr, "New trig scan has failed!!\n" );
             exit(1);
           }

        }
        else
        {
            if ( vrbflg )
                  fprintf (stderr,
                          "The event head appears to be null containing %d triggers \n",
                          XLALCountSnglInspiral ( (savedEvents.snglInspiralTable) ));
        }
    }
  }

  /* write sngl_inspiral table */
  if ( vrbflg ) fprintf( stdout, "  sngl_inspiral table...\n" );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status,
        &results, outputMask ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, savedEvents,
        outputMask ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

  while ( savedEvents.snglInspiralTable )
  {
    event = savedEvents.snglInspiralTable;
    savedEvents.snglInspiralTable = savedEvents.snglInspiralTable->next;
    LALFree( event );
  }

  /* write the signals injected in a template bank simulation */
  if ( bankSim && siminspiral.simInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sim_inspiral_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, siminspiral,
          sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

    while ( siminspiral.simInspiralTable )
    {
      thisSimInspiral = siminspiral.simInspiralTable;
      siminspiral.simInspiralTable = siminspiral.simInspiralTable->next;
      LALFree( thisSimInspiral );
    }
  }

  /* write the template bank simulation results to the xml */
  if ( bankSim && siminstparams.simInstParamsTable )
  {
    if ( vrbflg ) fprintf( stdout, "  sim_inst table...\n" );
    LAL_CALL( LALBeginLIGOLwXMLTable( &status,
          &results, sim_inst_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, siminstparams,
          sim_inst_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

    while ( siminstparams.simInstParamsTable )
    {
      thisSimInstParams = siminstparams.simInstParamsTable;
      siminstparams.simInstParamsTable = siminstparams.simInstParamsTable->next;
      LALFree( thisSimInstParams );
    }
  }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );
  if ( vrbflg ) fprintf( stdout, "done. XML file closed\n" );

  /* free the rest of the memory, check for memory leaks and exit */
  if ( tdFollowUpFiles )
  {
    free ( tdFollowUpFiles );
    while ( tdFollowUpEvents )
    {
      thisFollowUpEvent = tdFollowUpEvents;
      tdFollowUpEvents = tdFollowUpEvents->next;
      XLALFreeSnglInspiral( &thisFollowUpEvent );
    }
  }

  if ( injectionFile )
  {
    free ( injectionFile );
    while ( injections )
    {
      thisInj = injections;
      injections = injections->next;
      LALFree( thisInj );
    }
  }
  if ( approximantName) free( approximantName );
  if ( orderName )     free( orderName );
  if ( calCacheName )  free( calCacheName );
  if ( frInCacheName ) free( frInCacheName );
  if ( frInType )      free( frInType );
  if ( bankFileName )  free( bankFileName );
  if ( channelName )   free( channelName );
  if ( fqChanName )    free( fqChanName );

  /* Free the analyseThisTmplt vector */
  if ( analyseThisTmplt) LALFree (analyseThisTmplt);

  /* Free template time series vector */
  if (templateTimeSeriesVector) XLALDestroyREAL4Vector(templateTimeSeriesVector);

  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();

  /* print a success message to stdout for parsing by exitcode */
  fprintf( stdout, "%s: EXITCODE0\n", argv[0] );
  exit( 0 );
}

/* ------------------------------------------------------------------------- */

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) );\
LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME );\
LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name );\
LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype );\
LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define USAGE(a) \
fprintf( a,   "lalapps_inspiral [options]\n\n");\
fprintf( a, "  --help                       display this message\n");\
fprintf( a, "  --verbose                    print progress information\n");\
fprintf( a, "  --version                    print version information and exit\n");\
fprintf( a, "  --debug-level LEVEL          set the LAL debug level to LEVEL\n");\
fprintf( a, "  --user-tag STRING            set the process_params usertag to STRING\n");\
fprintf( a, "  --ifo-tag STRING             set the ifotag to STRING - for file naming\n");\
fprintf( a, "  --comment STRING             set the process table comment to STRING\n");\
fprintf( a, "\n");\
fprintf( a, "  --gps-start-time SEC         GPS second of data start time\n");\
fprintf( a, "  --gps-start-time-ns NS       GPS nanosecond of data start time\n");\
fprintf( a, "  --gps-end-time SEC           GPS second of data end time\n");\
fprintf( a, "  --gps-end-time-ns NS         GPS nanosecond of data end time\n");\
fprintf( a, "  --pad-data T                 pad the data start and end time by T seconds\n");\
fprintf( a, "  --slide-time T               slide data start epoch by T seconds\n");\
fprintf( a, "  --slide-time-ns T            slide data start epoch by T nanoseconds\n");\
fprintf( a, "\n");\
fprintf( a, "  --glob-frame-data            glob *.gwf files in the pwd to obtain frame data\n");\
fprintf( a, "  --frame-type TAG             input data is contained in frames of type TAG\n");\
fprintf( a, "  --frame-cache                obtain frame data from LAL frame cache FILE\n");\
fprintf( a, "  --calibration-cache FILE     obtain calibration from LAL frame cache FILE\n");\
fprintf( a, "  --glob-calibration-data      obtain calibration by globbing in working dir\n");\
fprintf( a, "\n");\
fprintf( a, "  --channel-name CHAN          read data from interferometer channel CHAN\n");\
fprintf( a, "  --calibrated-data TYPE       calibrated data of TYPE real_4 or real_8\n");\
fprintf( a, "  --strain-high-pass-freq F    high pass REAL8 h(t) data above F Hz\n");\
fprintf( a, "  --strain-high-pass-order O   set the order of the h(t) high pass filter to O\n");\
fprintf( a, "  --strain-high-pass-atten A   set the attenuation of the high pass filter to A\n");\
fprintf( a, "  --point-calibration          use the first point in the chunk to calibrate\n");\
fprintf( a, "\n");\
fprintf( a, "  --injection-file FILE        inject simulated inspiral signals from FILE\n");\
fprintf( a, "  --fast F                     analyse injections templates within a match F\n");\
fprintf( a, "  --inject-overhead            inject signals from overhead detector\n");\
fprintf( a, "  --enable-filter-inj-only     filter only segments with injections\n");\
fprintf( a, "  --disable-filter-inj-only    filter all segments when doing injections\n");\
fprintf( a, "  --hardware-injection         Injections are in the frame files don't make them again!\n");\
fprintf( a, "                               All segments are filtered.\n");\
fprintf( a, "\n");\
fprintf( a, "  --td-follow-up FILE          Follow up coincident BCV events in FILE\n");\
fprintf( a, "  --bank-file FILE             read template bank parameters from FILE\n");\
fprintf( a, "  --start-template N           start filtering at template number N in bank\n");\
fprintf( a, "  --stop-template N            stop filtering at template number N in bank\n");\
fprintf( a, "  --reverse-chirp-bank         filters data using a reverse chirp template bank\n");\
fprintf( a, "\n");\
fprintf( a, "  --sample-rate F              filter data at F Hz, downsampling if necessary\n");\
fprintf( a, "  --resample-filter TYPE       set resample filter to TYPE (ldas|butterworth) \n");\
fprintf( a, "\n");\
fprintf( a, "  --disable-high-pass          turn off the IIR highpass filter\n");\
fprintf( a, "  --enable-high-pass F         high pass data above F Hz using an IIR filter\n");\
fprintf( a, "  --high-pass-order O          set the order of the high pass filter to O\n");\
fprintf( a, "  --high-pass-attenuation A    set the attenuation of the high pass filter to A\n");\
fprintf( a, "  --spectrum-type TYPE         use PSD estimator TYPE\n");\
fprintf( a, "                                 (mean|median|gaussian|white|LIGO|AdvLIGO) \n");\
fprintf( a, "                               TYPE 'gaussian' is equivalent to 'white' - deprecated option\n");\
fprintf( a, "\n");\
fprintf( a, "  --segment-length N           set data segment length to N points\n");\
fprintf( a, "  --number-of-segments N       set number of data segments to N\n");\
fprintf( a, "  --segment-overlap N          overlap data segments by N points\n");\
fprintf( a, "\n");\
fprintf( a, "  --low-frequency-cutoff F     do not filter below F Hz\n");\
fprintf( a, "  --inverse-spec-length T      set length of inverse spectrum to T seconds\n");\
fprintf( a, "  --dynamic-range-exponent X   set dynamic range scaling to 2^X\n");\
fprintf( a, "\n");\
fprintf( a, "  --approximant APPROX         set approximant of the waveform to APPROX\n");\
fprintf( a, "                               (FindChirpSP|BCV|BCVC|BCVSpin|TaylorT1|TaylorT2|\n");\
fprintf( a, "                                  TaylorT3|PadeT1|EOB|GeneratePPN|FindChirpPTF) \n");\
fprintf( a, "  --order ORDER                set the pN order of the waveform to ORDER\n");\
fprintf( a, "                               (twoPN|twoPointFivePN|threePN|threePointFivePN|\n");\
fprintf( a, "                                  pseudoFourPN) \n");\
fprintf( a, "  --snr-threshold RHO          set signal-to-noise threshold to RHO\n");\
fprintf( a, "  --chisq-bins P               set number of chisq veto bins to P\n");\
fprintf( a, "  --chisq-delta DELTA          set chisq delta parameter to DELTA\n");\
fprintf( a, "  --chisq-threshold X          threshold on chi^2 < X * ( p + DELTA *rho^2 ) \n");\
fprintf( a, "  --cluster-method MTHD        max over chirp MTHD (tmplt|window|tmpltwindow|none) \n");\
fprintf( a, "  --cluster-window SEC         set length of clustering time window if required\n");\
fprintf( a, "\n");\
fprintf( a, "  --enable-rsq-veto            enable the r^2 veto test\n");\
fprintf( a, "  --disable-rsq-veto           disable the r^2 veto test\n");\
fprintf( a, "  --rsq-veto-window SEC        set the r^2 veto window to SEC\n");\
fprintf( a, "  --rsq-veto-threshold RSQ     set r^2 veto threshold to RSQ\n");\
fprintf( a, "\n");\
fprintf( a, "  --do-rsq-veto                do the r^2 veto\n");\
fprintf( a, "  --rsq-veto-time-thresh SEC   set the r^2 veto window to SEC\n");\
fprintf( a, "  --rsq-veto-max-snr MAXSNR    set the r^2 veto maximum snr to MAXSNR\n");\
fprintf( a, "  --rsq-veto-coeff COEFF       set the r^2 veto coefficient to COEFF\n");\
fprintf( a, "  --rsq-veto-pow POW           set the r^2 veto power to POW\n");\
fprintf( a, "\n");\
fprintf( a, "  --bank-veto-subbank-size N   set the number of tmplts in a subbank to N\n");\
fprintf( a, "\n");\
fprintf( a, "  --maximization-interval MSEC set length of interval (in ms) for\n");\
fprintf( a, "                                 maximization of triggers over the template bank.\n");\
fprintf( a, "                                 Cannot be used with --ts-cluster.\n");\
fprintf( a, "\n");\
fprintf( a, "  --ts-cluster   MTHD          max over template and end time MTHD \n");\
fprintf( a, "                                 (T0T3Tc|T0T3TcAS|Psi0Psi3Tc|Psi0Psi3TcAS) \n");\
fprintf( a, "                                 Cannot be used with --maximization-interval.\n");\
fprintf( a, "  --ts-endtime-interval msec   set end-time interval for TrigScan clustering\n");\
fprintf( a, "  --ts-metric-scaling fac      scale the metric which defines the ellipsoids for TrigScan\n");\
fprintf( a, "                               Scaling must be > 0\n");\
fprintf( a, "\n");\
fprintf( a, "\n");\
fprintf( a, "  --band-pass-template         Band-pass filter the time-domain inspiral template\n");\
fprintf( a, "  --taper-template OPT         Taper the inspiral template using option OPT\n");\
fprintf( a, "                                 (start|end|startend) \n");\
fprintf( a, "\n");\
fprintf( a, "  --enable-output              write the results to a LIGO LW XML file\n");\
fprintf( a, "  --output-mask MASK           write the output sngl_inspiral table\n");\
fprintf( a, "                                 with optional MASK (bns|bcv) \n");\
fprintf( a, "  --write-compress             write a compressed xml file\n");\
fprintf( a, "  --disable-output             do not write LIGO LW XML output file\n");\
fprintf( a, "  --trig-start-time SEC        only output triggers after GPS time SEC\n");\
fprintf( a, "  --trig-end-time SEC          only output triggers before GPS time SEC\n");\
fprintf( a, "\n");\
fprintf( a, "  --white-gaussian VAR         replace data with white gaussian noise of variance VAR\n");\
fprintf( a, "  --gaussian-noise VAR         same as --white-gaussian - deprecated option\n");\
fprintf( a, "  --colored-gaussian PSD       replace data with colored gaussian noise with psd PSD\n");\
fprintf( a, "                               (LIGO|AdvLIGO) \n");\
fprintf( a, "\n");\
fprintf( a, "  --random-seed SEED           set random number seed for injections to SEED\n");\
fprintf( a, "                                 (urandom|integer) \n");\
fprintf( a, "\n");\
fprintf( a, "  --bank-simulation N          perform N injections to test the template bank\n");\
fprintf( a, "                                (sim_inspiral.xml|integer) \n");\
fprintf( a, "  --enable-bank-sim-max        compute the maximum match over the bank\n");\
fprintf( a, "  --disable-bank-sim-max       do not maximize the match over the bank\n");\
fprintf( a, "  --sim-approximant APX        set approximant of the injected waveform to APX\n");\
fprintf( a, "                                 (TaylorT1|TaylorT2|TaylorT3|PadeT1|EOB|\n");\
fprintf( a, "                                  GeneratePPN|FrameFile) \n");\
fprintf( a, "  --sim-frame-file F           read the bank sim waveform from frame named F\n");\
fprintf( a, "  --sim-frame-channel C        read the bank sim waveform from frame channel C\n");\
fprintf( a, "  --sim-minimum-mass M         set minimum mass of bank injected signal to M\n");\
fprintf( a, "  --sim-maximum-mass M         set maximum mass of bank injected signal to M\n");\
fprintf( a, "  --bank-sim-flower F          set low frequency of signal to F\n");\
fprintf( a, "\n");\
fprintf( a, "  --data-checkpoint            checkpoint and exit after data is read in\n");\
fprintf( a, "  --checkpoint-path PATH       write checkpoint file under PATH\n");\
fprintf( a, "  --output-path PATH           write output data to PATH\n");\
fprintf( a, "\n");\
fprintf( a, "  --write-raw-data             write raw data to a frame file\n");\
fprintf( a, "  --write-filter-data          write data that is passed to filter to a frame\n");\
fprintf( a, "  --write-response             write the computed response function to a frame\n");\
fprintf( a, "  --write-spectrum             write the uncalibrated psd to a frame\n");\
fprintf( a, "  --write-snrsq                write the snr time series for each data segment\n");\
fprintf( a, "  --write-chisq                write the r^2 time series for each data segment\n");\
fprintf( a, "  --write-cdata                write the complex filter output\n");\
fprintf( a, "  --write-template                write the template time series\n");\
fprintf( a, "\n");

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"enable-output",           no_argument,       &enableOutput,     1 },
    {"write-compress",          no_argument,       &outCompress,      1 },
    {"disable-output",          no_argument,       &enableOutput,     0 },
    {"disable-high-pass",       no_argument,       &highPass,         0 },
    {"inject-overhead",         no_argument,       &injectOverhead,   1 },
    {"data-checkpoint",         no_argument,       &dataCheckpoint,   1 },
    {"glob-frame-data",         no_argument,       &globFrameData,    1 },
    {"glob-calibration-data",   no_argument,       &globCalData,      1 },
    {"point-calibration",       no_argument,       &pointCal,         1 },
    {"enable-rsq-veto",         no_argument,       &enableRsqVeto,    1 },
    {"disable-rsq-veto",        no_argument,       &enableRsqVeto,    0 },
    {"enable-filter-inj-only",  no_argument,       &flagFilterInjOnly,1 },
    {"disable-filter-inj-only", no_argument,       &flagFilterInjOnly,0 },
    {"hardware-injection",      no_argument,       &hardwareInjection,1 },
    {"reverse-chirp-bank",      no_argument,       &reverseChirpBank, 1 },
    {"do-rsq-veto",             no_argument,       &doRsqVeto,        1 },
    /* these options don't set a flag */
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-start-time-ns",       required_argument, 0,                'A'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"gps-end-time-ns",         required_argument, 0,                'B'},
    {"channel-name",            required_argument, 0,                'c'},
    {"segment-length",          required_argument, 0,                'd'},
    {"trig-start-time",         required_argument, 0,                'C'},
    {"trig-end-time",           required_argument, 0,                'D'},
    {"number-of-segments",      required_argument, 0,                'e'},
    {"segment-overlap",         required_argument, 0,                'f'},
    {"sample-rate",             required_argument, 0,                'g'},
    {"calibrated-data",         required_argument, 0,                'y'},
    {"strain-high-pass-freq",   required_argument, 0,                'E'},
    {"strain-high-pass-order",  required_argument, 0,                'P'},
    {"strain-high-pass-atten",  required_argument, 0,                'Q'},
    {"help",                    no_argument,       0,                'h'},
    {"low-frequency-cutoff",    required_argument, 0,                'i'},
    {"spectrum-type",           required_argument, 0,                'j'},
    {"inverse-spec-length",     required_argument, 0,                'k'},
    {"dynamic-range-exponent",  required_argument, 0,                'l'},
    {"start-template",          required_argument, 0,                'm'},
    {"chisq-delta",             required_argument, 0,                'M'},
    {"stop-template",           required_argument, 0,                'n'},
    {"chisq-bins",              required_argument, 0,                'o'},
    {"calibration-cache",       required_argument, 0,                'p'},
    {"approximant",             required_argument, 0,                'F'},
    {"order",                   required_argument, 0,                '^'},
    {"snr-threshold",           required_argument, 0,                'q'},
    {"chisq-threshold",         required_argument, 0,                'r'},
    {"resample-filter",         required_argument, 0,                'R'},
    {"comment",                 required_argument, 0,                's'},
    {"enable-high-pass",        required_argument, 0,                't'},
    {"high-pass-order",         required_argument, 0,                'H'},
    {"high-pass-attenuation",   required_argument, 0,                'T'},
    {"frame-cache",             required_argument, 0,                'u'},
    {"frame-type",              required_argument, 0,                'S'},
    {"bank-file",               required_argument, 0,                'v'},
    {"injection-file",          required_argument, 0,                'w'},
    {"pad-data",                required_argument, 0,                'x'},
    {"slide-time",              required_argument, 0,                'X'},
    {"slide-time-ns",           required_argument, 0,                'Y'},
    {"bank-simulation",         required_argument, 0,                'K'},
    {"sim-approximant",         required_argument, 0,                'L'},
    {"random-seed",             required_argument, 0,                'J'},
    {"sim-minimum-mass",        required_argument, 0,                'U'},
    {"sim-maximum-mass",        required_argument, 0,                'W'},
    {"bank-sim-flower",         required_argument, 0,                '?'},
    {"white-gaussian",          required_argument, 0,                'G'},
    {"gaussian-noise",          required_argument, 0,                'G'},
    {"colored-gaussian",        required_argument, 0,                '.'},
    {"checkpoint-path",         required_argument, 0,                'N'},
    {"output-path",             required_argument, 0,                'O'},
    {"debug-level",             required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"ifo-tag",                 required_argument, 0,                'I'},
    {"version",                 no_argument,       0,                'V'},
    {"maximization-interval",   required_argument, 0,                '@'},
    {"cluster-window",          required_argument, 0,                '#'},
    {"cluster-method",          required_argument, 0,                '0'},
    {"output-mask",             required_argument, 0,                '1'},
    {"fast",                    required_argument, 0,                '2'},
    {"rsq-veto-window",         required_argument, 0,                '3'},
    {"rsq-veto-threshold",      required_argument, 0,                '4'},
    {"enable-bank-sim-max",     no_argument,       0,                '5'},
    {"disable-bank-sim-max",    no_argument,       0,                '6'},
    {"sim-frame-file",          required_argument, 0,                '7'},
    {"sim-frame-channel",       required_argument, 0,                '8'},
    {"td-follow-up",            required_argument, 0,                '9'},
    {"ts-cluster",              required_argument, 0,                '*'},
    {"ts-endtime-interval",     required_argument, 0,                '<'},
    {"ts-metric-scaling",       required_argument, 0,                '>'},
    {"rsq-veto-time-thresh",    required_argument, 0,                '('},
    {"rsq-veto-max-snr",        required_argument, 0,                ')'},
    {"rsq-veto-coeff",          required_argument, 0,                '['},
    {"rsq-veto-pow",            required_argument, 0,                ']'},
    {"bank-veto-subbank-size",  required_argument, 0,                ','},
    {"band-pass-template",      no_argument,       0,                '}'},
    {"taper-template",          required_argument, 0,                '{'},
    /* frame writing options */
    {"write-raw-data",          no_argument,       &writeRawData,     1 },
    {"write-filter-data",       no_argument,       &writeFilterData,  1 },
    {"write-response",          no_argument,       &writeResponse,    1 },
    {"write-spectrum",          no_argument,       &writeSpectrum,    1 },
    {"write-snrsq",             no_argument,       &writeRhosq,       1 },
    {"write-chisq",             no_argument,       &writeChisq,       1 },
    {"write-cdata",             no_argument,       &writeCData,       1 },
    {"write-template",          no_argument,       &writeTemplate,       1 },
    {0, 0, 0, 0}
  };
  int c;
  INT4 haveDynRange = 0;
  INT4 haveApprox = 0;
  INT4 haveOrder = 0;
  INT4 haveClusterMethod = 0;
  INT4 haveBankSimApprox = 0;
  ProcessParamsTable *this_proc_param = procparams.processParamsTable;


  /*
   *
   * parse command line arguments
   *
   */


  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv,
        "-A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:VW:?:X:Y:Z:"
        "a:b:c:d:e:f:g:hi:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:"
        "0:1::2:3:4:567:8:9:*:>:<:(:):[:],:{:}:+:=:^:.:",
        long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
    {
      break;
    }

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'a':
        {
          long int gstartt = atol( optarg );
          if ( gstartt < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time is prior to "
                "Jan 01, 1994  00:00:00 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gstartt );
            exit( 1 );
          }
          if ( gstartt > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time is after "
                "Sep 14, 2011  01:46:26 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gstartt );
            exit( 1 );
          }
          gpsStartTimeNS += (INT8) gstartt * 1000000000LL;
          ADD_PROCESS_PARAM( "int", "%ld", gstartt );
        }
        break;

      case 'A':
        {
          long int gstarttns = atol( optarg );
          if ( gstarttns < 0 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time nanoseconds is negative\n",
                long_options[option_index].name );
            exit( 1 );
          }
          if ( gstarttns > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time nanoseconds is greater than unity:\n"
                "Must be <= 999999999 (%ld specified)\n",
                long_options[option_index].name, gstarttns );
            exit( 1 );
          }
          gpsStartTimeNS += (INT8) gstarttns;
          ADD_PROCESS_PARAM( "int", "%ld", gstarttns );
        }
        break;

      case 'b':
        {
          long int gendt = atol( optarg );
          if ( gendt > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time is after "
                "Sep 14, 2011  01:46:26 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gendt );
            exit( 1 );
          }
          else if ( gendt < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time is prior to "
                "Jan 01, 1994  00:00:00 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gendt );
            exit( 1 );
          }
          gpsEndTimeNS += (INT8) gendt * 1000000000LL;
          ADD_PROCESS_PARAM( "int", "%ld", gendt );
        }
        break;

      case 'B':
        {
          long int gendtns = atol( optarg );
          if ( gendtns < 0 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time nanoseconds is negative\n",
                long_options[option_index].name );
            exit( 1 );
          }
          else if ( gendtns > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time nanoseconds is greater than unity:\n"
                "Must be <= 999999999:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gendtns );
            exit( 1 );
          }
          gpsEndTimeNS += (INT8) gendtns;
          ADD_PROCESS_PARAM( "int", "%ld", gendtns );
        }
        break;

      case 'c':
        {
          /* create storage for the channel name and copy it */
          char *channamptr = NULL;
          optarg_len = strlen( optarg ) + 1;
          fqChanName = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
          memcpy( fqChanName, optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );

          /* check that we have a proper channel name */
          if ( ! (channamptr = strstr( fqChanName, ":" ) ) )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "channel name must be a full LIGO channel name "
                "e.g. L1:LSC-AS_Q\n(%s specified)\n",
                long_options[option_index].name, optarg );
            exit( 1 );
          }
          optarg_len = strlen( ++channamptr ) + 1;
          channelName = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
          memcpy( channelName, channamptr, optarg_len );

          /* copy the first two characters to the ifo name */
          memset( ifo, 0, sizeof(ifo) );
          memcpy( ifo, optarg, sizeof(ifo) - 1 );
        }
        break;

      case 'd':
        numPoints = (INT4) atoi( optarg );
        if ( numPoints < 2 || numPoints % 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of points must be a non-zero power of 2: "
              "(%d specified) \n",
              long_options[option_index].name, numPoints );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numPoints );
        break;

      case 'C':
        {
          long int gstartt = atol( optarg );
          /* ignore a value of zero */
          if ( gstartt )
          {
            if ( gstartt < 441417609 )
            {
              fprintf( stderr, "invalid argument to --%s:\n"
                  "GPS start time is prior to "
                  "Jan 01, 1994  00:00:00 UTC:\n"
                  "(%ld specified)\n",
                  long_options[option_index].name, gstartt );
              exit( 1 );
            }
            if ( gstartt > 999999999 )
            {
              fprintf( stderr, "invalid argument to --%s:\n"
                  "GPS start time is after "
                  "Sep 14, 2011  01:46:26 UTC:\n"
                  "(%ld specified)\n",
                  long_options[option_index].name, gstartt );
              exit( 1 );
            }
            trigStartTimeNS = (INT8) gstartt * 1000000000LL;
          }
          ADD_PROCESS_PARAM( "int", "%ld", gstartt );
        }
        break;

      case 'D':
        {
          long int gendt = atol( optarg );
          /* ignore a value of zero */
          if ( gendt )
          {
            if ( gendt > 999999999 )
            {
              fprintf( stderr, "invalid argument to --%s:\n"
                  "GPS end time is after "
                  "Sep 14, 2011  01:46:26 UTC:\n"
                  "(%ld specified)\n",
                  long_options[option_index].name, gendt );
              exit( 1 );
            }
            else if ( gendt < 441417609 )
            {
              fprintf( stderr, "invalid argument to --%s:\n"
                  "GPS end time is prior to "
                  "Jan 01, 1994  00:00:00 UTC:\n"
                  "(%ld specified)\n",
                  long_options[option_index].name, gendt );
              exit( 1 );
            }
            trigEndTimeNS = (INT8) gendt * 1000000000LL;
          }
          ADD_PROCESS_PARAM( "int", "%ld", gendt );
        }
        break;

      case 'e':
        numSegments = (INT4) atoi( optarg );
        if ( numSegments < 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of data segment must be greater than 0: "
              "(%d specified)\n",
              long_options[option_index].name, numSegments );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numSegments );
        break;

      case 'f':
        ovrlap = (INT4) atoi( optarg );
        if ( ovrlap < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "data segment overlap must be positive: "
              "(%d specified)\n",
              long_options[option_index].name, ovrlap );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", ovrlap );
        break;

      case 'g':
        sampleRate = (INT4) atoi( optarg );
        if ( sampleRate < 2 || sampleRate > 16384 || sampleRate % 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "rate must be power of 2 between 2 and 16384 inclusive: "
              "(%d specified)\n",
              long_options[option_index].name, sampleRate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", sampleRate );
        break;

      case 'y':
        /* specify which type of calibrated data */
        {
          if ( ! strcmp( "real_4", optarg ) )
          {
            calData = real_4;
          }
          else if ( ! strcmp( "real_8", optarg ) )
          {
            calData = real_8;
          }
          else
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "unknown data type specified;\n"
                "%s (must be one of: real_4, real_8)\n",
                long_options[option_index].name, optarg);
          }
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'E':
        strainHighPassFreq = (REAL4) atof( optarg );
        if ( strainHighPassFreq <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "REAL8 h(t) high pass filter frequency must be greater than 0 Hz:"
              "(%f Hz specified)\n",
              long_options[option_index].name, strainHighPassFreq );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", strainHighPassFreq );
        break;

      case 'P':
        strainHighPassOrder = (INT4) atoi( optarg );
        if ( strainHighPassOrder <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "REAL8 h(t) high pass filter order must be greater than 0: "
              "(%d specified)\n",
              long_options[option_index].name, strainHighPassOrder );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", strainHighPassOrder );
        break;

      case 'Q':
        strainHighPassAtten = (REAL4) atof( optarg );
        if ( strainHighPassAtten < 0.0 || strainHighPassAtten > 1.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "REAL8 h(t) high pass attenuation must be in the range [0:1]: "
              "(%f specified)\n",
              long_options[option_index].name, strainHighPassAtten );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", strainHighPassAtten );
        break;

      case 'h':
        USAGE( stdout );
        exit( 0 );
        break;

      case 'i':
        fLow = (REAL4) atof( optarg );
        if ( fLow < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "low frequency cutoff is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, fLow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fLow );
        break;

      case 'j':
        if ( ! strcmp( "mean", optarg ) )
        {
          specType = specType_mean;
        }
        else if ( ! strcmp( "median", optarg ) )
        {
          specType = specType_median;
        }
        else if ( ! strcmp( "bad-mean", optarg ) )
        {
          specType = specType_mean;
          badMeanPsd = 1;
        }
        else if ( (! strcmp( "gaussian", optarg )) || (! strcmp( "white", optarg )) )
        {
          specType = specType_gaussian;
          fprintf( stderr,
              "WARNING: replacing psd with white Gaussian spectrum\n" );
        }
        else if ( ! strcmp( "LIGO", optarg ) )
        {
          specType = specType_LIGO;
          fprintf( stderr,
              "WARNING: replacing psd with Initial LIGO design spectrum\n"
              "WARNING: replacing response function with a constant\n" );
        }
        else if ( ! strcmp( "AdvLIGO", optarg ) )
        {
          specType = specType_AdvLIGO;
          fprintf( stderr,
              "WARNING: replacing psd with Advanced LIGO design spectrum\n"
              "WARNING: replacing response function with a constant\n" );
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown power spectrum type: "
              "%s (must be mean, median, gaussian, white, LIGO or AdvLIGO)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'k':
        invSpecTrunc = (INT4) atoi( optarg );
        if ( invSpecTrunc < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "inverse spectrum length must be positive or zero: "
              "(%d specified)\n",
              long_options[option_index].name, invSpecTrunc );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", invSpecTrunc );
        break;

      case 'l':
        dynRangeExponent = (REAL4) atof( optarg );
        haveDynRange = 1;
        ADD_PROCESS_PARAM( "float", "%e", dynRangeExponent );
        break;

      case 'm':
        startTemplate = (INT4) atoi( optarg );
        if ( startTemplate < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "template bank start index must be positive: "
              "(%d specified)\n",
              long_options[option_index].name, startTemplate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", startTemplate );
        break;

      case 'M':
        chisqDelta = (REAL4) atof( optarg );
        if ( chisqDelta < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "chi squared delta parameter must be positive: "
              "(%e specified)\n",
              long_options[option_index].name, chisqDelta );
        }
        ADD_PROCESS_PARAM( "float", "%e", chisqDelta );
        break;

      case 'n':
        stopTemplate = (INT4) atoi( optarg );
        if ( stopTemplate < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "template bank stop index must be positive: "
              "(%d specified)\n",
              long_options[option_index].name, stopTemplate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", stopTemplate );
        break;

      case 'o':
        numChisqBins = (INT4) atoi( optarg );
        if ( numChisqBins < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of chisq veto bins must be positive: "
              "(%d specified)\n",
              long_options[option_index].name, numChisqBins );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numChisqBins );
        break;

      case 'p':
        /* create storage for the calibration frame cache name */
        optarg_len = strlen( optarg ) + 1;
        calCacheName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( calCacheName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'F':
        /* create storage for the approximant name */
        optarg_len = strlen( optarg ) + 1;
        approximantName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( approximantName, optarg, optarg_len );
        if ( ! strcmp( "TaylorT1", optarg ) )
        {
          approximant = TaylorT1;
        }
        else if ( ! strcmp( "TaylorT2", optarg ) )
        {
          approximant = TaylorT2;
        }
        else if ( ! strcmp( "TaylorT3", optarg ) )
        {
          approximant = TaylorT3;
        }
        else if ( ! strcmp( "GeneratePPN", optarg ) )
        {
          approximant = GeneratePPN;
        }
        else if ( ! strcmp( "PadeT1", optarg ) )
        {
          approximant = PadeT1;
        }
        else if ( ! strcmp( "EOB", optarg ) )
        {
          approximant = EOB;
        }
        else if ( ! strcmp( "EOBNR", optarg ) )
        {
          approximant = EOBNR;
        }
        else if ( ! strcmp( "FindChirpSP", optarg ) )
        {
          approximant = FindChirpSP;
        }
        else if ( ! strcmp( "BCV", optarg ) )
        {
          approximant = BCV;
        }
        else if ( ! strcmp( "BCVC", optarg ) )
        {
          approximant = BCV;
          bcvConstraint = 1;
        }
        else if ( ! strcmp( "BCVSpin", optarg ) )
        {
          approximant = BCVSpin;
        }
        else if ( ! strcmp( "FindChirpPTF", optarg ) )
        {
          approximant = FindChirpPTF;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be either FindChirpSP, BCV, BCVC, BCVSpin, FindChirpPTF\n"
              "TaylorT1, TaylorT2, TaylorT3, GeneratePPN, PadeT1 or EOB)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveApprox = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '^':
        /* create storage for the order name */
        optarg_len = strlen( optarg ) + 1;
        orderName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( orderName, optarg, optarg_len );
        if ( ! strcmp( "twoPN", optarg ) )
        {
          order = LAL_PNORDER_TWO;
        }
        else if ( ! strcmp( "twoPointFivePN", optarg ) )
        {
          order = LAL_PNORDER_TWO_POINT_FIVE;
        }
        else if ( ! strcmp( "threePN", optarg ) )
        {
          order = LAL_PNORDER_THREE;
        }
        else if ( ! strcmp( "threePointFivePN", optarg ) )
        {
          order = LAL_PNORDER_THREE_POINT_FIVE;
        }
        else if ( ! strcmp( "pseudoFourPN", optarg ) )
        {
          order = LAL_PNORDER_PSEUDO_FOUR;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be one of twoPN, twoPointFivePN, threePN, threePointFivePN, pseudoFourPN)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveOrder = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'q':
        snrThresh = atof( optarg );
        if ( snrThresh < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "signal to noise threshold must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, snrThresh );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'r':
        chisqThresh = atof( optarg );
        if ( chisqThresh < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "chi squared threshold must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, chisqThresh );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'R':
        if ( ! strcmp( "ldas", optarg ) )
        {
          resampFiltType = 0;
        }
        else if ( ! strcmp( "butterworth", optarg ) )
        {
          resampFiltType = 1;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown resampling filter type: "
              "%s (must be ldas or butterworth)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 's':
        if ( strlen( optarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          LALSnprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 't':
        highPass = 1;
        highPassFreq = (REAL4) atof( optarg );
        if ( highPassFreq <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "high pass filter frequency must be greater than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, highPassFreq );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", highPassFreq );
        break;

      case 'H':
        highPassOrder = (INT4) atoi( optarg );
        if ( highPassOrder <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "high pass filter order must be greater than 0: "
              "(%d specified)\n",
              long_options[option_index].name, highPassOrder );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", highPassOrder );
        break;

      case 'T':
        highPassAtten = (REAL4) atof( optarg );
        if ( highPassAtten < 0.0 || highPassAtten > 1.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "high pass attenuation must be in the range [0:1]: "
              "(%f specified)\n",
              long_options[option_index].name, highPassAtten );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", highPassAtten );
        break;

      case 'u':
        /* create storage for the input frame cache name */
        optarg_len = strlen( optarg ) + 1;
        frInCacheName = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( frInCacheName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'S':
        optarg_len = strlen( optarg ) + 1;
        frInType = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( frInType, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'v':
        /* create storage for the calibration frame cache name */
        optarg_len = strlen( optarg ) + 1;
        bankFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( bankFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'w':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectionFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectionFile, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'x':
        padData = (INT4) atoi( optarg );
        if ( padData < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of seconds to pad from input data"
              "must be greater than 0: (%d specified)\n",
              long_options[option_index].name, padData );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", padData );
        break;

      case 'X':
        slideData.gpsSeconds = (INT4) atoi( optarg );
        ADD_PROCESS_PARAM( "int", "%d", slideData.gpsSeconds );
        break;

      case 'Y':
        slideData.gpsNanoSeconds = (INT4) atoi( optarg );
        ADD_PROCESS_PARAM( "int", "%d", slideData.gpsNanoSeconds );
        break;

      case 'K':
        if ( strstr( optarg, "xml" ) )
        {
          /* we have been passed an xml file name */
          optarg_len = strlen( optarg ) + 1;
          bankSimFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
          memcpy( bankSimFileName, optarg, optarg_len );
          bankSim = 1;
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        else
        {
          /* parse the agument as an integer numer of simulations */
          bankSim = (INT4) atoi( optarg );
          if ( bankSim < 1 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "number of template bank simulations"
                "must be greater than 1: (%d specified)\n",
                long_options[option_index].name, bankSim );
            exit( 1 );
          }
          ADD_PROCESS_PARAM( "int", "%d", bankSim );
        }
        break;

      case 'L':
        if ( ! strcmp( "TaylorT1", optarg ) )
        {
          bankSimParams.approx = TaylorT1;
        }
        else if ( ! strcmp( "TaylorT2", optarg ) )
        {
          bankSimParams.approx = TaylorT2;
        }
        else if ( ! strcmp( "TaylorT3", optarg ) )
        {
          bankSimParams.approx = TaylorT3;
        }
        else if ( ! strcmp( "PadeT1", optarg ) )
        {
          bankSimParams.approx = PadeT1;
        }
        else if ( ! strcmp( "EOB", optarg ) )
        {
          bankSimParams.approx = EOB;
        }
        else if ( ! strcmp( "GeneratePPN", optarg ) )
        {
          bankSimParams.approx = GeneratePPN;
        }
        else if ( ! strcmp( "FrameFile", optarg ) )
        {
          bankSimParams.approx = FrameFile;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: %s\n(must be one of TaylorT1, "
              "TaylorT2, TaylorT3, PadeT1, EOB, GeneratePPN, FrameFile)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveBankSimApprox = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'J':
        if ( ! strcmp( "urandom", optarg ) )
        {
          randSeedType = urandom;
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        else
        {
          randSeedType = user;
          randomSeed = (INT4) atoi( optarg );
          ADD_PROCESS_PARAM( "int", "%d", randomSeed );
        }
        break;

      case 'U':
        bankSimParams.minMass = (REAL4) atof( optarg );
        if ( bankSimParams.minMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "miniumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, bankSimParams.minMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", bankSimParams.minMass );
        break;

      case 'W':
        bankSimParams.maxMass = (REAL4) atof( optarg );
        if ( bankSimParams.maxMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maxiumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, bankSimParams.maxMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", bankSimParams.maxMass );
        break;

      case 'G':
        gaussVar = (REAL4) atof( optarg );
        if ( gaussVar < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "variance of white Gaussian noise must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, gaussVar );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", gaussVar );
        whiteGaussian = 1;
        unitResponse = 1;
        fprintf( stderr,
            "WARNING: replacing input data with white Gaussian noise\n"
            "WARNING: replacing response function with unity\n" );
        break;

      case '.':
        colorSpec = (INT4) atoi( optarg );
        if ( ! strcmp( "LIGO", optarg ) )
        {
          colorSpec = 3;
          fprintf( stderr,
              "WARNING: replacing input data with colored Gaussian noise: "
              "psd = Initial LIGO\n");
        }
        else if ( ! strcmp( "AdvLIGO", optarg ) )
        {
          colorSpec = 4;
          fprintf( stderr,
              "WARNING: replacing input data with colored Gaussian noise: "
              "psd = Advanced LIGO\n");
        }
        else
        {
          fprintf(stderr,"invalid power spectrum for colored Gaussian noise;"
              "colorSpec must be either LIGO or advLIGO");
        }
        coloredGaussian = 1;
        break;

      case 'N':
        if ( LALSnprintf( ckptPath, FILENAME_MAX * sizeof(CHAR),
              "%s", optarg ) < 0 )
        {
          fprintf( stderr, "invalid argument to --%s\n"
              "local path %s too long: string truncated\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );

      case 'O':
        if ( LALSnprintf( outputPath, FILENAME_MAX * sizeof(CHAR),
              "%s", optarg ) < 0 )
        {
          fprintf( stderr, "invalid argument to --%s\n"
              "output path %s too long: string truncated\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );

      case 'z':
        set_debug_level( optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
        break;

      case 'I':
        /* create storaged for the ifo-tag */
        optarg_len = strlen( optarg ) + 1;
        ifoTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( ifoTag, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "LIGO/LSC Standalone Inspiral Search Engine\n"
            "Duncan Brown <duncan@gravity.phys.uwm.edu>\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
	fprintf( stdout, lalappsGitID );
        exit( 0 );
        break;

      case '#':
        clusterWindow = (REAL4) atof( optarg );
        if ( clusterWindow <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "clustering time window is less than or equal to 0 secs: "
              "(%f secs specified)\n",
              long_options[option_index].name, clusterWindow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", clusterWindow );
        break;

      case '@':
        {
          long int maxms = atol( optarg );
          if ( maxms < 0 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "maximization interval is less than 0 msecs: "
                "(%ld msecs specified)\n",
                long_options[option_index].name, maxms );
            exit( 1 );
          }
	  /* internally we require maximizationInterval to be in nano seconds */
 	  /* This will be passed as an argument in the call to                */
          /* XLALMaxSnglInspiralOverIntervals (). Therefore multiply by       */
          /* 1000000 to convert msec to nano seconds                          */
          maximizationInterval = (INT4) maxms * 1000000;
          ADD_PROCESS_PARAM( "int", "%ld", maxms );
        }
        break;

      case '0':
        if ( ! strcmp( "none", optarg ) )
        {
          clusterMethod = FindChirpClustering_none;
        }
        else if ( ! strcmp( "template", optarg ) )
        {
          clusterMethod = FindChirpClustering_tmplt;
        }
        else if ( ! strcmp( "window", optarg ) )
        {
          clusterMethod = FindChirpClustering_window;
        }
        else if ( ! strcmp( "tmpltwindow", optarg ) )
        {
          clusterMethod = FindChirpClustering_tmpltwindow;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown clustering method: "
              "%s (must be 'none', 'template', 'window' or 'tmpltwindow')\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveClusterMethod = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '1':
        if ( ! strcmp( "bns", optarg ) )
        {
          outputMask = sngl_inspiral_table_bns;
        }
        else if ( ! strcmp( "bcv", optarg ) )
        {
          outputMask = sngl_inspiral_table_bcv;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "must be either bns or bcv: (%s specified)\n",
              long_options[option_index].name, optarg);
          exit (1);
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '2':
        {
          float  mmFastArg = atof (optarg);
          if ( mmFastArg < (float)(0.0) || mmFastArg > (float)(1.0) )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "should be between 0.0 and 1.0: (%f specified)\n",
                long_options[option_index].name, mmFastArg);
            exit (1);
          }
          mmFast = (REAL4) mmFastArg;
          ADD_PROCESS_PARAM( "float", "%e", mmFast );
        }
        break;

      case '3':
        rsqVetoWindow = atol( optarg );
        if ( rsqVetoWindow < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "r^2 veto time window is less than or equal to 0 msec: "
              "(%f msecs specified)\n",
              long_options[option_index].name, rsqVetoWindow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '4':
        rsqVetoThresh = atof( optarg );
        if ( rsqVetoThresh < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              " r^2 veto threshold must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, rsqVetoThresh );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '5':
        bankSimParams.maxMatch = 1;
        ADD_PROCESS_PARAM( "string", "%s", " " );
        break;

      case '6':
        bankSimParams.maxMatch = 0;
        ADD_PROCESS_PARAM( "string", "%s", " " );
        break;

      case '7':
        /* create storage for the bank sim frame file name */
        optarg_len = strlen( optarg ) + 1;
        bankSimParams.frameName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( bankSimParams.frameName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '8':
        /* create storage for the bank sim frame file channel */
        optarg_len = strlen( optarg ) + 1;
        bankSimParams.frameChan = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( bankSimParams.frameChan, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '9':
        numTDFiles = 1;
        /* Count the number of thinca files to follow up */
        while ( !strstr( argv[optind], "--" ) )
        {
          numTDFiles++;
          optind++;
        }
        optind = optind - numTDFiles;

        /* Set pointers to the relevant filenames */
        tdFollowUpFiles = (CHAR **) calloc( numTDFiles, sizeof(CHAR *));

        numTDFiles = 0;
        while ( !strstr( argv[optind], "--" ) )
        {
          tdFollowUpFiles[numTDFiles++] = argv[optind];
          ADD_PROCESS_PARAM( "string", "%s", argv[optind] );
          optind++;

        }
        break;

      case '*':
        /* store trigSanClustering method */
        if ( ! strcmp( "T0T3Tc", optarg ) )
        {
            trigScanMethod = T0T3Tc;
            trigScanAppendStragglers = 0;
        }
        else if ( ! strcmp( "T0T3TcAS", optarg ) )
        {
            trigScanMethod = T0T3Tc;
            trigScanAppendStragglers = 1;
        }
        else if ( ! strcmp( "Psi0Psi3Tc", optarg ) )
        {
            trigScanMethod = Psi0Psi3Tc;
            trigScanAppendStragglers = 0;
        }
        else if ( ! strcmp( "Psi0Psi3TcAS", optarg ) )
        {
            trigScanMethod = Psi0Psi3Tc;
            trigScanAppendStragglers = 1;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown scan method specified: %s\n"
              "(Must be one of T0T3Tc, T0T3TcAS, Psi0Psi3Tc, Psi0Psi3TcAS)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '<':
        /* TrigScan Delta End Time */
        trigScanDeltaEndTime = atof( optarg );
        if ( trigScanDeltaEndTime < 0.0L )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "ts-endtime-interval must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, trigScanDeltaEndTime );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '>':
        /* TrigScan Template Metric Scaling Factor */
        trigScanMetricScalingFac = atof( optarg );
        if ( trigScanMetricScalingFac <= 0.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "ts-volume-safety must be > 0.0 : "
              "(%f specified)\n",
              long_options[option_index].name, trigScanMetricScalingFac );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '?':
 	bankSimParams.f_lower = (REAL4) atof( optarg );
 	if (bankSimParams.f_lower <= 0 )
 	{
 		fprintf( stderr, "invalid argument to --%s:\n"
 				"bank-sim-flower must be > 0.0 : "
 				"(%f specified)\n",
 				long_options[option_index].name, bankSimParams.f_lower );
 		exit( 1 );
 	}
 	ADD_PROCESS_PARAM( "float", "%e", bankSimParams.f_lower );
        break;

      case '(':
        rsqVetoTimeThresh = atof( optarg );
        if ( rsqVetoTimeThresh < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              " r^2 veto time threshold must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, rsqVetoTimeThresh );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case ')':
        rsqVetoMaxSNR = atof( optarg );
        if ( rsqVetoMaxSNR < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              " r^2 veto maximum snr must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, rsqVetoMaxSNR );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '[':
        rsqVetoCoeff = atof( optarg );
        if ( rsqVetoCoeff < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              " r^2 veto coefficient must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, rsqVetoCoeff );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case ']':
        rsqVetoPow = atof( optarg );
        if ( rsqVetoPow < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              " r^2 veto power must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, rsqVetoPow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case ',':
        {
          int subBankSizeArg = atoi( optarg );

          if ( subBankSizeArg < 1 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                " number of templates in a subbank must be greater than zero: "
                "(%d specified)\n",
                long_options[option_index].name, subBankSizeArg );
            exit( 1 );
          }
          ADD_PROCESS_PARAM( "int", "%s", optarg );

          subBankSize = (UINT4) subBankSizeArg;
        }
        break;
      case '}':
        bandPassTmplt = 1;
        ADD_PROCESS_PARAM( "int", "%s", " " );
        break;

      case '{':
        if ( !strcmp( "start", optarg ) )
        {
          taperTmplt = INSPIRAL_TAPER_START;
        }
        else if ( !strcmp( "end", optarg ) )
        {
          taperTmplt = INSPIRAL_TAPER_END;
        }
        else if ( !strcmp( "startend", optarg ) )
        {
          taperTmplt = INSPIRAL_TAPER_STARTEND;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "Taper must be set to start, end or startend:(%s specified)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;



      default:
        fprintf( stderr, "unknown error while parsing options (%d)\n", c );
        exit( 1 );
    }
  }

  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc )
    {
      fprintf ( stderr, "%s\n", argv[optind++] );
    }
    exit( 1 );
  }

  /* enable output is stored in the first process param row */
  if ( enableOutput == 1 )
  {
    LALSnprintf( procparams.processParamsTable->program,
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    LALSnprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--enable-output" );
    LALSnprintf( procparams.processParamsTable->type,
        LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( procparams.processParamsTable->value,
        LIGOMETA_VALUE_MAX, " " );
  }
  else if ( enableOutput == 0 )
  {
    LALSnprintf( procparams.processParamsTable->program,
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    LALSnprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--disable-output" );
    LALSnprintf( procparams.processParamsTable->type,
        LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( procparams.processParamsTable->value,
        LIGOMETA_VALUE_MAX, " " );
  }
  else
  {
    fprintf( stderr, "--enable-output or --disable-output "
        "argument must be specified\n" );
    exit( 1 );
  }

  /* check inject-overhead option */
  if ( injectOverhead )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--inject-overhead" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  /* store point calibration option */
  if ( pointCal )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--point-calibration" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  if ( outCompress)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--write-compress" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  /*
   *
   * check validity of arguments
   *
   *
   */


  /* check validity of input data time */
  if ( ! gpsStartTimeNS )
  {
    fprintf( stderr, "--gps-start-time must be specified\n" );
    exit( 1 );
  }
  XLALINT8NSToGPS( &gpsStartTime, gpsStartTimeNS );
  if ( ! gpsEndTimeNS )
  {
    fprintf( stderr, "--gps-end-time must be specified\n" );
    exit( 1 );
  }
  XLALINT8NSToGPS( &gpsEndTime, gpsEndTimeNS );
  if ( gpsEndTimeNS <= gpsStartTimeNS )
  {
    fprintf( stderr, "invalid gps time range: "
        "start time: %d, end time %d\n",
        gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds );
    exit( 1 );
  }


      /* check trigger generation time is within input time */
  if ( trigStartTimeNS )
  {
    if ( trigStartTimeNS < gpsStartTimeNS )
    {
      fprintf( stderr,
          "trigStartTimeNS = %ld\nis less than gpsStartTimeNS = %ld",
          trigStartTimeNS, gpsStartTimeNS );
    }
  }
  if ( trigEndTimeNS )
  {
    if ( trigEndTimeNS > gpsEndTimeNS )
    {
      fprintf( stderr,
          "trigEndTimeNS = %ld\nis greater than gpsEndTimeNS = %ld",
          trigEndTimeNS, gpsEndTimeNS );
    }
  }

  /* check validity of data length parameters */
  if ( numPoints < 0 )
  {
    fprintf( stderr, "--segment-length must be specified\n" );
    exit( 1 );
  }
  if ( numSegments < 0 )
  {
    fprintf( stderr, "--number-of-segments must be specified\n" );
    exit( 1 );
  }
  if ( ovrlap < 0 )
  {
    fprintf( stderr, "--segment-overlap must be specified\n" );
    exit( 1 );
  }

  /* check sample rate has been given */
  if ( sampleRate < 0 )
  {
    fprintf( stderr, "--sample-rate must be specified\n" );
    exit( 1 );
  }

  /* check high pass option has been given */
  if ( highPass < 0 )
  {
    fprintf( stderr, "--disable-high-pass or --enable-high-pass (freq)"
        " must be specified\n" );
    exit( 1 );
  }
  else if ( ! highPass )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--disable-high-pass" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }
  else
  {
    /* check that all the high pass parameters have been specified */
    if ( highPassOrder < 0 )
    {
      fprintf( stderr, "--high-pass-order must be specified\n" );
      exit( 1 );
    }
    if ( highPassAtten < 0 )
    {
      fprintf( stderr, "--high-pass-attenuation must be specified\n" );
      exit( 1 );
    }
  }

  if ( calData == real_8 )
  {
    /* check that strain high pass parameters have been specified */
    if ( strainHighPassFreq < 0 )
    {
      fprintf( stderr,
          "--strain-high-pass-freq must be specified for REAL8 h(t) data\n" );
      exit( 1 );
    }
    if ( strainHighPassOrder < 0 )
    {
      fprintf( stderr,
          "--strain-high-pass-order must be specified for REAL8 h(t) data\n" );
      exit( 1 );
    }
    if ( strainHighPassAtten < 0 )
    {
      fprintf( stderr,
          "--strain-high-pass-atten must be specified for REAL8 h(t) data\n" );
      exit( 1 );
    }
  }

  /* check validity of input data length */
  inputDataLength = numPoints * numSegments - ( numSegments - 1 ) * ovrlap;
  {
    INT8 gpsChanIntervalNS = gpsEndTimeNS - gpsStartTimeNS;
    INT8 inputDataLengthNS = (INT8) inputDataLength * 1000000000LL /
      (INT8) sampleRate;

    if ( inputDataLengthNS != gpsChanIntervalNS )
    {
      fprintf( stderr, "length of input data and data chunk do not match\n" );
      fprintf( stderr, "start time: %ld, end time %ld\n",
          gpsStartTimeNS / 1000000000LL, gpsEndTimeNS / 1000000000LL );
      fprintf( stderr, "gps channel time interval: %ld ns\n"
          "computed input data length: %ld ns\n",
          gpsChanIntervalNS, inputDataLengthNS );
      exit( 1 );
    }
  }

  if ( fLow < bankSimParams.f_lower )
  {
	  fprintf( stderr, "--low-frequency-cutoff must be greater than bank sim injection starting frequency\n" );
	  exit( 1 );
  }

  /* check filter parameters have been specified */
  if ( numChisqBins < 0 )
  {
    fprintf( stderr, "--chisq-bins must be specified\n" );
    exit( 1 );
  }
  if ( fLow < 0 )
  {
    fprintf( stderr, "--low-frequency-cutoff must be specified\n" );
    exit( 1 );
  }
  if ( resampFiltType < 0 )
  {
    fprintf( stderr, "--resample-filter must be specified\n" );
    exit( 1 );
  }
  if ( specType == specType_undefined )
  {
    fprintf( stderr, "--spectrum-type must be specified\n" );
    exit( 1 );
  }
  if ( clusterWindow * sampleRate > numPoints )
  {
    fprintf( stderr, "--cluster-window must be less than "
        "--segment-length\n" );
    exit( 1 );
  }
  if ( ! haveClusterMethod )
  {
    fprintf( stderr, "--cluster-method must be specified\n" );
    exit( 1 );
  }
  if ( ( clusterMethod == FindChirpClustering_window || clusterMethod == FindChirpClustering_tmpltwindow ) && clusterWindow == -1 )
  {
    fprintf( stderr, "--cluster-window must be specified "
        "if --clustering method 'window' chosen\n" );
    exit( 1 );
  }
  if ( clusterMethod != FindChirpClustering_window && clusterMethod != FindChirpClustering_tmpltwindow && clusterWindow != -1 )
  {
    fprintf( stderr, "--cluster-window specified "
        "but --clustering method 'window' or 'tmpltwindow' not chosen\n" );
    exit( 1 );
  }
  if ( invSpecTrunc < 0 )
  {
    fprintf( stderr, "--inverse-spec-length must be specified\n" );
    exit( 1 );
  }
  else if ( invSpecTrunc * sampleRate > numPoints )
  {
    fprintf( stderr, "--inverse-spec-length must be less than "
        "--segment-length\n" );
    exit( 1 );
  }

  if ( ! haveDynRange )
  {
    fprintf( stderr, "--dynamic-range-exponent must be specified\n" );
    exit( 1 );
  }
  if ( ! haveApprox )
  {
    fprintf( stderr, "--approximant must be specified\n" );
    exit( 1 );
  }

  /* check that the pN order is specified */
  if ( ! haveOrder )
  {
    fprintf( stderr, "--order must be specified\n" );
    exit( 1 );
  }


  /* check that a channel has been requested and fill the ifo */
  if ( ! fqChanName )
  {
    fprintf( stderr, "--channel-name must be specified\n" );
    exit( 1 );
  }

  if ( ! bankSim )
  {
    /* check that the thresholds have been specified */
    if ( snrThresh < 0 )
    {
      fprintf( stderr, "--snr-threshold must be specified\n" );
      exit( 1 );
    }
    if ( chisqThresh < 0 )
    {
      fprintf( stderr, "--chisq-threshold must be specified\n" );
      exit( 1 );
    }
    if ( chisqDelta < 0 )
    {
      fprintf( stderr, "--chisq-delta must be specified\n" );
      exit( 1 );
    }
  }

  /* check that we can correctly obtain the input frame data */
  if ( globFrameData )
  {
    if ( frInCacheName )
    {
      fprintf( stderr,
          "--frame-cache must not be specified when globbing frame data\n" );
      exit( 1 );
    }

    if ( ! frInType )
    {
      fprintf( stderr,
          "--frame-type must be specified when globbing frame data\n" );
      exit( 1 );
    }
  }
  else
  {
    if ( ! frInCacheName )
    {
      fprintf( stderr,
          "--frame-cache must be specified when not globbing frame data\n" );
      exit( 1 );
    }

    if ( frInType )
    {
      fprintf( stderr, "--frame-type must not be specified when obtaining "
          "frame data from a cache file\n" );
      exit( 1 );
    }
  }

  /* record the glob frame data option in the process params */
  if ( globFrameData )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--glob-frame-data" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  /* check we can calibrate the data if it's not h(t) */
  if ( ! calData )
  {
    if ( ! ( calCacheName || globCalData ) )
    {
      fprintf( stderr, "either --calibration-cache or "
          "--glob-calibration-data must be specified\n" );
      exit( 1 );
    }
    else if ( calCacheName && globCalData )
    {
      fprintf( stderr, "only one of --calibration-cache or "
          "--glob-calibration-data can be specified\n" );
      exit( 1 );
    }
  }
  else
  {
    if ( calCacheName || globCalData )
    {
      fprintf( stderr, "neither --calibration-cache nor "
          "--glob-calibration-data\nshould be given for calibrated data\n" );
      exit( 1 );
    }
  }

  /* record the glob calibration data option in the process params */
  if ( globCalData )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--glob-calibration-data" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  /* check that a template bank has been specified */
  if ( ! bankFileName )
  {
    fprintf( stderr, "--bank-file must be specified\n" );
    exit( 1 );
  }

  /* check FindChirpSP is used if reverse chirp bank is specified */
  if ( reverseChirpBank )
  {
    if ( reverseChirpBank && ! ( approximant == FindChirpSP ) )
    {
      fprintf( stderr, "--approximant must be FindChirpSP if "
          "--reverse-chirp-bank is specified\n" );
      exit( 1 );
    }
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--reverse-chirp-bank" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  /* check that a random seed for gaussian noise generation has been given */
  if ( (whiteGaussian || coloredGaussian ) && randSeedType == unset )
  {
    fprintf( stderr, "--random-seed must be specified if "
        "--white-gaussian (--gaussian-noise) or --colored-gaussian are given\n" );
    exit( 1 );
  }


  /* check that specType is either LIGO or advLIGO for colored Gaussian noise
   */
  if (coloredGaussian)
  {
    /* check that --white-gaussian has not been specified */
    if (whiteGaussian)
    {
      fprintf( stderr, "cannot specify options --white-gaussian (or "
          " --gaussian-noise) and --colored-gaussian at the same time;"
          " spectrum must be either white or colored \n");
      exit( 1 );
    }

    /* check that specType is either LIGO or advLIGO for colored Gaussian noise
     */
    if (specType == specType_gaussian)
    {
      fprintf( stderr, "--spectrum-type cannot be white for "
                        "colored Gaussian noise: specify either LIGO, "
                        "AdvLIGO, mean or median\n");
      exit( 1 );
    }
    if (specType == specType_LIGO && colorSpec == 4)
    {
      fprintf(stderr,"Error: if "
        "--colored-gaussian is AdvLIGO --spectrum-type cannot be LIGO\n");
      exit( 1 );
    }
    if (specType == specType_AdvLIGO && colorSpec == 3)
    {
      fprintf(stderr,"Error: if "
        "--colored-gaussian is LIGO --spectrum-type cannot be AdvLIGO\n");
      exit( 1 );
    }
     /* check that strain high pass parameters have been specified */
    if ( strainHighPassFreq < 0 )
    {
      fprintf( stderr,
          "--strain-high-pass-freq must be specified for simulated colored "
          "Gaussian noise\n" );
      exit( 1 );
    }
    if ( strainHighPassOrder < 0 )
    {
      fprintf( stderr,
          "--strain-high-pass-order must be specified for simulated colored "
          "Gaussian noise \n" );
      exit( 1 );
    }
    if ( strainHighPassAtten < 0 )
    {
      fprintf( stderr,
          "--strain-high-pass-atten must be specified for simulated colored "
          "Gaussian noise\n" );
      exit( 1 );
    }
  }

  /* check the bank simulation parameters */
  if ( bankSim )
  {
    if ( bankSimParams.maxMatch < 0 )
    {
      fprintf( stderr, "--enable-bank-sim-max or --disable-bank-sim-max "
          "must be specified\n" );
      exit( 1 );
    }

    if ( ! bankSimFileName )
    {
      if ( ! haveBankSimApprox )
      {
        fprintf( stderr, "--sim-approximant must be specified if "
            "--bank-simulation is not a sim_inspiral file\n" );
        exit( 1 );
      }

      if ( bankSimParams.approx == FrameFile )
      {
        if ( ! bankSimParams.frameName || ! bankSimParams.frameChan )
        {
          fprintf( stderr, "--sim-frame-file and --sim-frame-channel\n"
              "must be specified if --sim-approximant is FrameFile\n" );
          exit( 1 );
        }
      }
      else
      {
        if ( randSeedType == unset )
        {
          fprintf( stderr, "--random-seed must be specified if\n"
              "--bank-simulation is given and approx is not FrameFile\n" );
          exit( 1 );
        }
        if ( bankSimParams.minMass < 0 )
        {
          fprintf( stderr, "--sim-minimum-mass must be specified if\n"
              "--bank-simulation is given and approx is not FrameFile\n" );
          exit( 1 );
        }
        if ( bankSimParams.maxMass < 0 )
        {
          fprintf( stderr, "--sim-maximum-mass must be specified if\n"
              "--bank-simulation is given and approx is not FrameFile\n" );
          exit( 1 );
        }
        if ( bankSimParams.f_lower < 0 )
        {
          fprintf( stderr, "--bank-sim-flower must be specified if\n"
              "--bank-simulation is given and approx is not FrameFile\n" );
          exit( 1 );
        }
      }
    }
  }

  /* check rsq veto flags */
  if ( enableRsqVeto == 0 )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--disable-rsq-veto" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }
  else if ( enableRsqVeto == 1 )
  {
    if ( (rsqVetoWindow < 0) || (rsqVetoThresh < 0) || (numChisqBins < 2) )
    {
      fprintf( stderr, "\n rsq-veto-window, rsq-veto-thresh and chisq-bins "
          " > 1\n must "
          "be specified if the --enable-rsq-veto argument is given\n" );
      exit( 1 );
    }

    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--enable-rsq-veto" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }
  else
  {
    fprintf( stderr,
        "one of --enable-rsq-veto or --disable-rsq-veto must be specified\n" );
    exit( 1 );
  }

  if ( doRsqVeto == 1 )
  {
    if ( enableRsqVeto == 0 )
    {
      fprintf( stderr, "--enable-rsq-veto must be specified \n"
          "if the --do-rsq-veto argument is specified\n" );
      exit( 1 );
    }
    else if ( ( rsqVetoTimeThresh < 0 ) || ( rsqVetoMaxSNR < 0 ) )
    {
      fprintf( stderr, "--rsq-veto-time-thresh and "
          "--rsq-veto-max-snr must be \n"
          "specified if the --do-rsq-veto argument is specified\n" );
      exit( 1 );
    }
    else if ( ( rsqVetoCoeff > 0 ) &&
        ( ( rsqVetoTimeThresh < 0 ) ||
          ( rsqVetoMaxSNR < 0 ) || ( rsqVetoPow < 0 ) ) )
    {
      fprintf( stderr, "--rsq-veto-time-thresh --rsq-veto-max-snr "
          "and --rsq-veto-pow \n"
          "must be specified if --rsq-veto-coeff is specified\n" );
      exit( 1 );
    }
    else if ( ( rsqVetoPow > 0 ) &&
        ( ( rsqVetoTimeThresh < 0 ) || ( rsqVetoMaxSNR < 0 ) ||
          ( rsqVetoCoeff < 0 ) ) )
    {
      fprintf( stderr, "--rsq-veto-time-thresh "
          "--rsq-veto-max-snr and --rsq-veto-coeff \n"
          "must be specified if --rsq-veto-pow is specified\n" );
      exit( 1 );
    }

    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--do-rsq-veto" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  /* check to filter injection segments only */
  if ( flagFilterInjOnly == 1 )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--enable-filter-inj-only" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }
  else if ( flagFilterInjOnly == 0 )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--disable-filter-inj-only" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }
  else if ( flagFilterInjOnly == -1 && injectionFile )
  {
    fprintf( stderr, "one of --enable-filter-inj-only or "
        "--disable-filter-inj-only must be specified\n" );
    exit( 1 );
  }

  /* Check the trigScan input parameters */
  if ( trigScanMethod )
  {
    if ( trigScanMetricScalingFac <= 0.0 )
    {
      fprintf ( stderr, "You must specify --ts-metric-scaling\n" );
      exit(1);
    }

    if ( maximizationInterval )
    {
      fprintf ( stderr, "Cannot specify both --maximization-interval"
          " and --ts-cluster \nChoose any one of the two methods"
          " for clustering raw inspiral triggers.\n" );
      exit(1);
    }
  }
  /* If the trigScan parameters are reasonable, set trigScanDeltaEndTime*/
  /* Here we change trigScanDeltaEndTime to msec                        */
  trigScanDeltaEndTime /= 1000.L;

  /* check the bank veto input parameters */
  if ( ! subBankSize )
  {
    fprintf( stderr, "--bank-veto-subbank-size must be specifed.\n"
        "Use --bank-veto-subbank-size 1 to disable template bank veto\n" );
    exit( 1 );
  }
  else
  {
    /*if ( numChisqBins > 0 && subBankSize != 1 )
    {
      fprintf( stderr,
          "Template bank veto is incompatible with --num-chisq-bins %d\n"
          "Use --bank-veto-subbank-size 1 to disable template bank veto\n",
          numChisqBins );
      exit( 1 );
    }*/
    if ( mmFast >= 0 )
    {
      fprintf( stderr, "Template bank veto is incompatible --fast %f\n"
          "Use --bank-veto-subbank-size 1 to disable template bank veto\n",
          mmFast );
      exit( 1 );
    }
    if ( numTDFiles > 0 )
    {
      fprintf( stderr, "Template bank veto is incompatible --td-follow-up\n"
          "Use --bank-veto-subbank-size 1 to disable template bank veto\n" );
      exit( 1 );
    }
  }

  return 0;
}


#undef ADD_PROCESS_PARAM
