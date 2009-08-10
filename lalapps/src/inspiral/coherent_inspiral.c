/*
*  Copyright (C) 2007 Sukanta Bose, Duncan Brown, Robert Adam Mercer, Stephen Fairhurst, Sean Seader
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
 * File Name: coherent_inspiral.c
 *
 * Author: Bose, S. and Seader, S. and Rogan, A.
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

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <lal/LALFrameIO.h>
#include <lal/TimeSeries.h>

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

#include <lal/LALRCSID.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/FrameStream.h>
#include <lal/DataBuffer.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpBCVSpin.h>
#include <lal/FindChirpChisq.h>
#include <lal/DetectorSite.h>
#include <lal/Random.h>
#include <lal/LALInspiral.h>
#include <lal/CoherentInspiral.h>
#include <lal/LALStatusMacros.h>
#include <lal/SkyCoordinates.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "coherent_inspiral"
#define CVS_NAME_STRING "$Name$"

#define ADD_SUMM_VALUE( sv_name, sv_comment, val, intval ) \
if ( this_summ_value ) \
{ \
  this_summ_value = this_summ_value->next = (SummValueTable *) \
  LALCalloc( 1, sizeof(SummValueTable) ); \
} \
else \
{ \
  summvalue.summValueTable = this_summ_value = (SummValueTable *) \
  LALCalloc( 1, sizeof(SummValueTable) ); \
} \
snprintf( this_summ_value->program, LIGOMETA_PROGRAM_MAX, "%s", \
  PROGRAM_NAME ); \
this_summ_value->version = 0; \
this_summ_value->start_time = searchsumm.searchSummaryTable->out_start_time; \
this_summ_value->end_time = searchsumm.searchSummaryTable->out_end_time; \
this_summ_value->value = (REAL4) val; \
this_summ_value->intvalue = (INT4) intval; \
snprintf( this_summ_value->name, LIGOMETA_SUMMVALUE_NAME_MAX, "%s", \
    sv_name ); \
snprintf( this_summ_value->ifo, LIGOMETA_IFO_MAX, "%s", ifo ); \
snprintf( this_summ_value->comment, LIGOMETA_SUMMVALUE_COMM_MAX, \
    "%s", sv_comment ); \

#define rint(x) (floor((x)+0.5))

int arg_parse_check( int argc, char *argv[], MetadataTable procparams );


/*
 *
 * variables that control program behaviour
 *
 */

/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* input data parameters */

INT4   sampleRate           = -1;       /* sample rate of filter data   */
INT4   numPointsSeg         = -1;/* set to segment-length used in inspiral.c */
REAL4  fLow                 = -1;       /* low frequency cutoff         */
REAL4  dynRangeExponent     = -1;/* set to same value used in inspiral.c */

/*Coherent code specific inputs*/

char  *ifoframefile[6] = {NULL,NULL,NULL,NULL,NULL,NULL};

INT4 H1file = 0;
INT4 H2file = 0;
INT4 L1file = 0;
INT4 G1file = 0;
INT4 T1file = 0;
INT4 V1file = 0;
                                
/* input time-slide parameters */
REAL8  slideStep[6]     = {0.0,0.0,0.0,0.0,0.0,0.0};
int    bankDuration     = 0;
CHAR   cohbankFileName[FILENAME_MAX]; /* name of input template bank */
/* CHAR  *cohbankFileName = NULL; name of input template bank  */
UINT4  cohSNROut            = 0;    /* default is not to write frame */
UINT4  cohH1H2SNROut      = 0;    /* default is not to write frame */
UINT4  nullStatOut       = 0;    /* default is not to write frame */
UINT4  nullStatH1H2Out       = 0;    /* default is not to write frame */
UINT4  eventsOut            = 0;    /* default is not to write events */
REAL4  cohSNRThresh         = -1;
UINT4  maximizeOverChirp    = 0;    /* default is no clustering */
INT4   verbose              = 0;
CHAR   outputPath[FILENAME_MAX];
UINT4  outCompress = 0;
INT4   numCohTrigs      = 400; /*optional argument that can be reset */
INT8  gpsStartTimeNS   = 0;         /* input data GPS start time ns */
LIGOTimeGPS gpsStartTime;           /* input data GPS start time    */
INT8  gpsEndTimeNS     = 0;         /* input data GPS end time ns   */
LIGOTimeGPS gpsEndTime;             /* input data GPS end time      */

double raStep = 6.0;
double decStep = 3.0;
int  gpsStartTimeTemp   = 0;         /* input data GPS start time ns */
int  gpsEndTimeTemp   = 0;         /* input data GPS start time ns */
INT8   outTimeNS        = 0;            /* search summ out time    */

LALStatus             status;
LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

CHAR  *userTag          = NULL;         /* string the user can tag with */
CHAR  *ifos           = NULL;         /* string to tag parent IFOs    */

/* Params to convert from geocentric to equatorial */

ConvertSkyParams             convertParams;
SkyPosition                  tempSky;
MultiInspiralTable           *thisEventTemp = NULL;

int main( int argc, char *argv[] )
{
  /* output */

  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsumm;
  MetadataTable         searchsummvars;
  MetadataTable         summvalue;
  ProcessParamsTable   *this_proc_param = NULL;
  SearchSummvarsTable  *this_search_summvar = NULL;
  SummValueTable       *this_summ_value = NULL;

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummaryTable   *searchSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;


  LIGOLwXMLStream       results;
  REAL4 alpha = 0;
  REAL4 alphabeta = 0;
  REAL4 inj_alpha = 0;
  REAL4 inj_alphabeta = 0;

  CHAR   fileName[FILENAME_MAX];
  CHAR   framename[FILENAME_MAX];
  CHAR   xmlname[FILENAME_MAX];
  CHAR   cohdataStr[LALNameLength];
  CHAR   caseIDChars[6][LIGOMETA_IFOS_MAX] = {"0","0","0","0","0","0"};
  /* type of channels, e.g., LSC-STRAIN, in c-data frames          */
  CHAR   channelNameArray[6][LALNameLength] = {"0","0","0","0","0","0"}; 
  CHAR  *ifo = NULL;

  REAL4   cohSegLength     = 1.0; /* This should match hardcoded value of 0.125s in inspiral.c */
  INT4   numPoints        = 0;
  UINT4  numSegments      = 1; /* CHECK: number of segments; hardwired to 1 */
  UINT4  numBeamPoints    = 0; /* CHECK: number of sky position templates; now redundant */
  UINT8  eventID          = 0;
  REAL8  cohSeriesLength;

  REAL8  m1               = 0.0;
  REAL8  m2               = 0.0;
  REAL8  dynRange         = 0.0;

  /* variables for initializing tempTime to account for time-slides */
  UINT8  triggerNumber    = 0;
  UINT8  slideNumber      = 0;
  UINT8  slideSign        = 0;

  /* counters and other variables */
  INT4   j,k,l, cohtrigs, cohtriglimit;
  INT4   kmax = 0;
  UINT4  numDetectors     = 0;
  REAL8  tempTime[6]      = {0.0,0.0,0.0,0.0,0.0,0.0};
  INT4   timeptDiff[5]    = {0,0,0,0,0};
  UINT2  caseID[6]        = {0,0,0,0,0,0}; /* H1 L V G T H2 */
  INT4   numTriggers      = 0;
  INT4   numCoincs        = 0;
  UINT4  numCohFiles      = 1;
  UINT4  cohFileID        = 1;

  REAL4 totMass = 0.0;
  REAL8 muMass = 0.0;
  REAL8 deltaT= 0.0;
  REAL4 distNorm = 0.0;
  REAL4 templateNorm = 0.0;
  REAL8 c0= 0.0;
  REAL8 c2= 0.0;
  REAL4 c3= 0.0;
  REAL8 c4= 0.0;
  REAL4 x1= 0.0;
  REAL4 x2= 0.0;
  REAL4 x3= 0.0;
  REAL4 x4= 0.0;
  REAL4 x8= 0.0;
  
  SnglInspiralTable    *currentTrigger = NULL;
  SnglInspiralTable    *cohbankEventList=NULL;

  CoincInspiralTable   *coincHead = NULL;
  CoincInspiralTable   *thisCoinc = NULL;

  CoherentInspiralInitParams   *cohInspInitParams = NULL;
  CoherentInspiralFilterParams *cohInspFilterParams = NULL;
  CoherentInspiralFilterInput  *cohInspFilterInput = NULL;
  CoherentInspiralCVector      *cohInspCVec = NULL;
  CoherentInspiralBeamVector   *cohInspBeamVec = NULL;
  MultiInspiralTable           *thisEvent = NULL;
  MultiInspiralTable           *tempTable = NULL;
  MetadataTable                 savedEvents;
  COMPLEX8TimeSeries            tempSnippet;

  char nameArrayCData[6][256] = {"0","0","0","0","0","0"}; /* cData chan names */

  /* set default debugging level */
  XLALSetErrorHandler( XLALAbortErrorHandler );
  set_debug_level( "1" ); /* change with parse option */

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );

  if (strcmp(CVS_REVISION, "$Revi" "sion$"))
  {
    XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME,
        CVS_REVISION, CVS_SOURCE, CVS_DATE, 0);
  }
  else
  {
    XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME,
        lalappsGitCommitID, lalappsGitGitStatus, lalappsGitCommitDate, 0);
  }
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *)
    calloc( 1, sizeof(ProcessParamsTable) );

  /* create the search summary and zero out the summvars table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );
  searchsummvars.searchSummvarsTable = NULL;

  arg_parse_check( argc, argv, procparams );
  if (vrbflg)  fprintf(stdout, "called parse options..\n");

  /* wind to the end of the process params table */
  for ( this_proc_param = procparams.processParamsTable; this_proc_param->next;
        this_proc_param = this_proc_param->next );

  /* Set other variables */
  numPoints = rint(sampleRate * cohSegLength);
  savedEvents.multiInspiralTable = NULL;
  k = 0; 

  /* store the input sample rate */
  this_search_summvar = searchsummvars.searchSummvarsTable = 
    (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
  snprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
	       "data sample rate" );
  this_search_summvar->value = (REAL8) sampleRate;
  
  
  /* Set the dynamic range; needed for distNorm, templateNorm calculation */
  dynRange = pow( 2.0, dynRangeExponent );
  
  /* read in the cohbank trigger ligo lw xml file */
  numTriggers = XLALReadInspiralTriggerFile( &cohbankEventList,
					     &currentTrigger, &searchSummList, 
					     &inputFiles, cohbankFileName );
  
  fprintf(stdout,"Reading templates from %s\n",cohbankFileName);
  
  if ( numTriggers < 0 )
    {
      fprintf(stderr, "Error reading triggers from file %s", cohbankFileName);
      exit( 1 );
    }
  else if ( numTriggers == 0 )
    { 
      if( vrbflg )
	{
	  fprintf( stdout, 
		   "%s contains no triggers - the coherent bank will be empty\n",
		   cohbankFileName );
	}
    }
  else {
    if( vrbflg )
      {
	fprintf( stdout, 
		 "Read in %d triggers from the file %s\n", numTriggers, 
		 cohbankFileName );
      }
    
    /* reconstruct the coincs */
    numCoincs = XLALRecreateCoincFromSngls( &coincHead, &cohbankEventList );

    /* save the number of coincidences in the search summary table */
    searchsumm.searchSummaryTable->nevents = numCoincs;

    if( numCoincs < 0 )
      {
	fprintf(stderr, "Unable to reconstruct coincs from single ifo triggers");
	exit( 1 );
      }
    else
      {
	fprintf( stdout,
		 "Recreated %d coincs from the %d triggers\n", numCoincs, 
		 numTriggers );
	if ( numCoincs % numCohTrigs == 0) {
	  
	  numCohFiles = floor( numCoincs / numCohTrigs );
	}
	else {
	  numCohFiles = floor( numCoincs / numCohTrigs ) + 1;
	}
      }
    
    thisCoinc=coincHead;
    
    numCohFiles++;

    /* loop over coincident triggers to compute cohSNR */
    
    for( cohFileID = 1 ; cohFileID < numCohFiles ; cohFileID++ ) {
      
      /* frame output data */
      struct FrFile *frOutFile  = NULL;
      struct FrameH *outFrameCoh   = NULL;
      struct FrameH *outFrameCohH1H2SNR   = NULL;
      struct FrameH *outFrameNullStat   = NULL;
      struct FrameH *outFrameNullStatH1H2   = NULL;      
      
      if( (cohFileID == (numCohFiles-1) && !((numCoincs % numCohTrigs) == 0)) )
	{
	  cohtriglimit = (numCoincs % numCohTrigs);
	}
      else 
	{
	  cohtriglimit = numCohTrigs; 
	}
      
      for ( cohtrigs = 0 ; cohtrigs < cohtriglimit ; cohtrigs++) {
	
	numDetectors = thisCoinc->numIfos;
	
	/* l is another detector index, which can have a max value of 4 */
	l=0;
	
	/* Note the participating ifos and the eventID
	   for this coincident trigger */
	for( k=0 ; k<LAL_NUM_IFO ; k++)
	  {
	    if( thisCoinc->snglInspiral[k] )
	      {
		kmax = k; /* final trigger's k value */
		caseID[k] = 1;
		memcpy( caseIDChars[k], &thisCoinc->snglInspiral[k]->ifo, sizeof(caseIDChars[k] - 1) );
		snprintf( channelNameArray[k], LALNameLength*sizeof(CHAR), "%s", &thisCoinc->snglInspiral[k]->channel );
		eventID = thisCoinc->snglInspiral[k]->event_id->id;
		if( vrbflg ) fprintf(stdout,"eventID = %Ld\n",eventID );
		
		/* Parse eventID to get the slide number */
		triggerNumber = eventID % 100000;
		slideNumber = ((eventID % 100000000) - triggerNumber)/100000;
		slideSign = (eventID % 1000000000) - slideNumber*100000 - triggerNumber;
		
		/* Store CData frame name now for reading its frame-file 
		   later, within thisCoinc-ident loop
		*/
		snprintf( nameArrayCData[k], LALNameLength*sizeof(CHAR), "%s:CBC-CData_%Ld", caseIDChars[k], eventID );
	      }
	  }/* Closes loop for k; finished noting the participating ifos 
	      and the eventID for this coincident trigger*/
	
	/* Initialize tempTime to account for time-slides (wrt H1, hence j<5)*/
	for( j=0; j<(LAL_NUM_IFO - 1) ; j++)
	  {
	    /* slideSign=0 is the same as a positive time slide */
	    if(slideSign != 0)
	      {
		tempTime[j] = slideStep[j]*slideNumber;
	      }
	    else
	      {
		tempTime[j] -= slideStep[j]*slideNumber;
	      }
	  }
	
	l=0;
	if( G1file ) l++;
	if( H1file ) l++;
	if( H2file ) l++;
	if( L1file ) l++;
	if( T1file ) l++;
	if( V1file ) l++;
	
	if( (INT4)numDetectors != l )
	  {
	    fprintf( stderr, "You have events for %d detectors, but specified frame files for %d detectors\n",numDetectors,l);
	    if( (INT4)numDetectors > l )
	      {
		fprintf( stderr, "You must specify more frame files. Exiting...\n");
		exit(1);
	      }
	    else
	      {
		if( vrbflg ) fprintf( stdout, "One or more of the frame files specified will not be used for this event since the number of detectors is less than the number of frame files you specified.\n");
	      }
	  }
	
	l = 0;
	
	if( vrbflg ) fprintf(stdout,"numDetectors = %d\n", numDetectors);
	if( vrbflg ) fprintf(stdout,"caseID = %d %d %d %d %d %d (G1,H1,H2,L1,T1,V1)\n", caseID[0], caseID[1], caseID[2], caseID[3], caseID[4], caseID[5]);
	
	
	/* Initialize the necessary structures for thisCoinc-ident trigger*/
	
	if( !(cohInspInitParams = (CoherentInspiralInitParams *) calloc(1,sizeof(CoherentInspiralInitParams)) ))
	  {
	    fprintf( stdout, "could not allocate memory for coherentInspiral init params\n" );
	    exit(1);
	  }
	
	/* Initialize the coherent param structure for thisCoinc trigger */      
	cohInspInitParams->numDetectors            = numDetectors;
	cohInspInitParams->numSegments             = numSegments;
	cohInspInitParams->numPoints               = numPoints;
	cohInspInitParams->numBeamPoints           = numBeamPoints;
	cohInspInitParams->cohSNROut               = cohSNROut;
	/* In addition to the network cohSNR, output the cohH1H2SNR if
	   the user wants it and the network has the ifos H1 and H2; since
	   in a 2D network this will make one of cohSNR and cohH1H2SNR
	   redundant, do not output the latter for < 3D networks */
	if( cohH1H2SNROut && caseID[1] && caseID[2] && (numDetectors > 2) ) {
	  cohInspInitParams->cohH1H2SNROut       = 1;
	}
	else if( cohSNROut && caseID[1] && caseID[2] && (numDetectors == 2)){
	  cohInspInitParams->cohH1H2SNROut       = 1;
	}
	else {
	  if ( vrbflg && cohH1H2SNROut ) fprintf( stdout, "Not outputting cohH1H2SNR because either numDetectors < 3 or at least one of H1 and H2 is missing ...\n " );
	  cohInspInitParams->cohH1H2SNROut       = 0;
	}
	
	/* Determine if the H1-H2 null-statistic should be computed */
	if( (nullStatH1H2Out && ( (caseID[1] && caseID[2])
				  || (numDetectors > 2) ) ) ) {
	  cohInspInitParams->nullStatH1H2Out         = 1;
	}
	else {
	  if ( vrbflg && nullStatH1H2Out ) fprintf( stdout, "Not outputting nullStatH1H2Out because either this is a two-detector network or at least one of H1 and H2 is missing ...\n " );
	  cohInspInitParams->nullStatH1H2Out         = 0;
	}
	
	if( (nullStatOut && ( ( numDetectors >3 ) || 
			      ( (numDetectors == 3) && !(caseID[1] && caseID[2]) ) ) ) ) {
	  cohInspInitParams->nullStatOut         = 1;
	}
	else {
	  if ( vrbflg && nullStatOut ) fprintf( stdout, "Not outputting nullStatOut because either numDetectors < 3 or at least one of H1 and H2 is missing ...\n " );
	  cohInspInitParams->nullStatOut         = 0;
	}
	
	/* create the data structures needed for coherentInspiral */
	
	if ( vrbflg ) fprintf( stdout, "initializing coherentInspiral...\n " );
	
	/* initialize coherentInspiral filter functions */
	LAL_CALL( LALCoherentInspiralFilterInputInit (&status, &cohInspFilterInput,
						      cohInspInitParams), &status );
	
	cohInspCVec = cohInspFilterInput->multiCData;
	
	m1 = thisCoinc->snglInspiral[kmax]->mass1;
	m2 = thisCoinc->snglInspiral[kmax]->mass2;
	muMass =  m1 * m2 / (m1 + m2);
	
	cohInspFilterInput->tmplt = (InspiralTemplate *)
	  LALCalloc(1,sizeof(InspiralTemplate) );
	cohInspFilterInput->tmplt->mass1 = m1;
	cohInspFilterInput->tmplt->mass2 = m2;
	cohInspFilterInput->tmplt->totalMass = m1 + m2;
	cohInspFilterInput->tmplt->mu = m1 * m2 / (m1 + m2);
	cohInspFilterInput->tmplt->eta = (m1 * m2) / ((m1 + m2) * (m1 + m2 ));
	
	if (vrbflg)  fprintf( stdout, "m1:%f m2:%f totalmass:%f mu:%f eta:%f\n", cohInspFilterInput->tmplt->mass1, cohInspFilterInput->tmplt->mass2,cohInspFilterInput->tmplt->totalMass,cohInspFilterInput->tmplt->mu,cohInspFilterInput->tmplt->eta);
	
	LAL_CALL( LALCoherentInspiralFilterParamsInit (&status, &cohInspFilterParams,
						       cohInspInitParams),&status );
	
	/* Initialize the filter param structure for thisCoinc trigger */ 
	
	cohInspFilterParams->deltaT                  = 1.0/((REAL8) sampleRate);
	cohInspFilterParams->cohSNRThresh            = cohSNRThresh;  
	cohInspFilterParams->cohSNROut               = cohInspInitParams->cohSNROut;
	cohInspFilterParams->cohH1H2SNROut           = cohInspInitParams->cohH1H2SNROut;
	cohInspFilterParams->nullStatH1H2Out         = cohInspInitParams->nullStatH1H2Out;
	cohInspFilterParams->nullStatOut             = cohInspInitParams->nullStatOut;
	cohInspFilterParams->numTmplts               = 1;
	cohInspFilterParams->fLow                    = fLow;
	cohInspFilterParams->maximizeOverChirp       = maximizeOverChirp;
	cohInspFilterParams->numDetectors            = cohInspInitParams->numDetectors;
	cohInspFilterParams->raStep                  = raStep;
	cohInspFilterParams->decStep                 = decStep;
	/* initParams not needed anymore */
	free( cohInspInitParams );
	cohInspInitParams = NULL;
	
	if (vrbflg)  fprintf( stdout, "deltaT:%f cohSNRThresh:%f numTmplts:%d\n", cohInspFilterParams->deltaT,cohInspFilterParams->cohSNRThresh,cohInspFilterParams->numTmplts);
	
	for( j=0; j<LAL_NUM_IFO; j++ ) 
	  {
	    cohInspFilterParams->detIDVec->data[j] = caseID[j];
	    cohInspFilterParams->sigmasqVec->data[j] = 1.0;
	  }
	
	/* Read in the snippets associated with thisCoinc trigger */
	l = 0; /* A counter to step through the cohInspCVec-tors */
	for( j=0; j<LAL_NUM_IFO; j++ )
	  {
	    if( caseID[j] )		{
	      FrFile *frfileIn[j];
	      FrameH *ifoFrame = NULL;
	      FrProcData *proc = NULL;
	      LIGOTimeGPS tmpEpoch;
	      
	      if( vrbflg ) fprintf(stdout, "getting the COMPLEX8TimeSeries %s \n", nameArrayCData[j] );
	      
	      if ( !(frfileIn[j]= XLALFrOpenURL( ifoframefile[j] )) ) {
		XLALPrintError( "XLAL Error: could not open frame file %s - exiting...\n", ifoframefile[j] );
		exit(1);
	      }		
	      while( ! proc && (ifoFrame = FrameRead( frfileIn[j] )) ) {
		proc = ifoFrame->procData;
                if( vrbflg ) fprintf(stdout, "Proc name is %s \n", proc->name );
		/*while ( proc && ( strcmp( nameArrayCData[j], proc->name) == 0 ) ) {*/
                while ( proc && strcmp( nameArrayCData[j], proc->name ) ) {
		  proc = proc->next;
		}
		if ( ! proc )
		  FrameFree( ifoFrame ); 
	      }
	     
	      if ( ! proc ) {
		XLALPrintError( "XLAL Error: could not find channel %s in file %s - exiting...\n", nameArrayCData[j], ifoframefile[j] );
		FrFileIEnd( frfileIn[j] );
		exit(1);
	      }
	      
	      XLALGPSSet( &tmpEpoch, ifoFrame->GTimeS, ifoFrame->GTimeN );
	      XLALGPSAdd( &tmpEpoch, proc->timeOffset );
	      XLALGPSAdd( &tmpEpoch, proc->data->startX[0] );
	      
	      /* Scale up fShift for all ifos */
	      if ( proc->fShift ) {
		cohInspCVec->cData[l] = XLALCreateCOMPLEX8TimeSeries( proc->data->name, &tmpEpoch, (REAL8) proc->fShift, proc->data->dx[0], &lalDimensionlessUnit, proc->data->nData );
	      }
	      else {
		cohInspCVec->cData[l] = XLALCreateCOMPLEX8TimeSeries( proc->data->name, &tmpEpoch, 0.0, proc->data->dx[0], &lalDimensionlessUnit, proc->data->nData );
	      }
	      if ( ! cohInspCVec->cData[l] ) {
		FrFileIEnd( frfileIn[j] );
		XLALPrintError( "XLAL Error: could not create cData from channel %s in file %s - exiting...\n", nameArrayCData[j], ifoframefile[j] );
		exit(1);
	      }
	      
	      memcpy( cohInspCVec->cData[l]->data->data, proc->data->data, cohInspCVec->cData[l]->data->length * sizeof( *(cohInspCVec->cData[l])->data->data ) );
	      FrFileIEnd( frfileIn[j] );
	      
	      /*CHECK: Note that inspiral.c has a temporary fix for using 
		the f0 member of the COMPLEX8TimeSeries structure for 
		communicating sigmasq.
	      */
	      if ( !(cohInspCVec->cData[l]->f0 == 0.0) ) {
		cohInspFilterParams->sigmasqVec->data[j] = cohInspCVec->cData[l]->f0;
	      }
	      else {
		cohInspFilterParams->sigmasqVec->data[j] = 1.0;
	      }		
	      
	      if (vrbflg)  fprintf( stdout, "sigmasq:%f\n",cohInspFilterParams->sigmasqVec->data[j]);
	      
	      l++;
	    }/* Closes "if( caseID[j] )" */
	  }/* Closes "for( j=0; j<LAL_NUM_IFO; j++ )" */
	
	for ( l=0 ; l<(INT4)numDetectors ; l++)
	  { 
	    tempTime[l] += cohInspCVec->cData[l]->epoch.gpsSeconds + cohInspCVec->cData[l]->epoch.gpsNanoSeconds * 1e-9;
	  }
	
	/* store the start and end time of the first ifo cdata in the search summary */
	/* FIXME:  loss of precision;  consider
	searchsumm.searchSummaryTable->in_start_time = searchsumm.searchSummaryTable->in_end_time = cohInspCVec->cData[0]->epoch;
	XLALGPSAdd(&searchsumm.searchSummaryTable->in_end_time, (REAL8) cohSegLength / sampleRate);
	*/
	searchsumm.searchSummaryTable->in_start_time = cohInspCVec->cData[0]->epoch;
	cohSeriesLength = XLALGPSGetREAL8( &(cohInspCVec->cData[0]->epoch) );
	cohSeriesLength += (REAL8) cohSegLength / sampleRate ;
	XLALGPSSetREAL8( &(searchsumm.searchSummaryTable->in_end_time), cohSeriesLength );

	/* If we can estimate distance then compute templateNorm */
	/* At present, this is only good for frequency domain tmplts */
	/* Since each detector's data has been filtered with templates */
	/* that have the same mass pair, templateNorm is the same for */
	/* every detector and needs to be computed only once.         */
	
	totMass  = (REAL4) cohInspFilterInput->tmplt->totalMass;
	deltaT = (REAL4) cohInspFilterParams->deltaT;
	distNorm = 2.0 * LAL_MRSUN_SI / (1.0 * 1e6 * LAL_PC_SI);
	templateNorm = sqrt( (5.0*((REAL4)muMass)) / 96.0 ) *  pow( totMass / (LAL_PI*LAL_PI) , 1.0/3.0 ) * pow( LAL_MTSUN_SI / deltaT, -1.0/6.0 );
	distNorm *= dynRange;
	templateNorm *= templateNorm;
	templateNorm *= distNorm * distNorm;
	cohInspFilterParams->templateNorm = templateNorm;
	cohInspFilterParams->segmentLength = numPointsSeg;
	
	/* calculate the length of the chirp for clustering over chirp-length*/
	c0 = 5.0*(cohInspFilterInput->tmplt->totalMass)*LAL_MTSUN_SI/(256.0*(cohInspFilterInput->tmplt->eta));
	c2 = 743.0/252.0 + (cohInspFilterInput->tmplt->eta)*11.0/3.0;
	c3 = -32*LAL_PI/3;
	c4 = 3058673.0/508032.0 + (cohInspFilterInput->tmplt->eta)*(5429.0/504.0 + (cohInspFilterInput->tmplt->eta)*617.0/72.0);
	x1 = pow(LAL_PI*((REAL4)(cohInspFilterInput->tmplt->totalMass))*LAL_MTSUN_SI*cohInspFilterParams->fLow , 1.0/3.0);
	x2 = x1*x1;
	x3 = x1*x2;
	x4 = x2*x2;
	x8 = x4*x4;
	
	cohInspFilterParams->chirpTime = ((REAL4)c0)*(1 + ((REAL4)c2)*x2 + c3*x3 + ((REAL4)c4)*x4)/x8;
	
	if (vrbflg) fprintf(stdout,"filtering the data..\n");
	if ( maximizeOverChirp && vrbflg )
	  {
	    fprintf(stdout,"clustering events\n");
	  }                    
	
	/* Before the data gets filtered, I need to make the c-data snippets commensurate */
	for(j=0;j<(INT4)numDetectors - 1;j++)
	  {
	    timeptDiff[j] = rint((tempTime[0] - tempTime[j+1]) * sampleRate);
	    if( vrbflg ) fprintf(stdout,"timeptDiff = %d\n",timeptDiff[j]);
	  }
	
	/* Now allocate memory for a temporary storage vector */
	memset( &tempSnippet, 0, sizeof(COMPLEX8TimeSeries) );
	LAL_CALL( LALCCreateVector( &status, &(tempSnippet.data), numPoints ), &status );
	
	/* If cohSNR is being output, then copy epoch */
	if( cohInspFilterParams->cohSNROut ) {
	  cohInspFilterParams->cohSNRVec->epoch = cohInspCVec->cData[0]->epoch;
	  cohInspFilterParams->cohSNRVec->deltaT = cohInspCVec->cData[0]->deltaT;
	}
	
	if( cohInspFilterParams->cohH1H2SNROut ) {
	  cohInspFilterParams->cohH1H2SNRVec->epoch = cohInspCVec->cData[0]->epoch;
	  cohInspFilterParams->cohH1H2SNRVec->deltaT = cohInspCVec->cData[0]->deltaT;
	}
	
	if( cohInspFilterParams->nullStatH1H2Out ) {
	  cohInspFilterParams->nullStatH1H2Vec->epoch = cohInspCVec->cData[0]->epoch;
	  cohInspFilterParams->nullStatH1H2Vec->deltaT = cohInspCVec->cData[0]->deltaT;
	}
	
	if( cohInspFilterParams->nullStatOut ) {
	  cohInspFilterParams->nullStatVec->epoch = cohInspCVec->cData[0]->epoch;
	  cohInspFilterParams->nullStatVec->deltaT = cohInspCVec->cData[0]->deltaT;
	}
	
	/* Now that the time series are commensurate, do the filtering... */
	LALCoherentInspiralFilterSegment (&status, &thisEvent, cohInspFilterInput, cohInspFilterParams);
	
	/* Save event id in multi_inspiral table */
	thisEventTemp =thisEvent;
	while( thisEventTemp )
	  {
	    thisEventTemp->event_id = (EventIDColumn *) 
	      LALCalloc(1, sizeof(EventIDColumn) );
	    thisEventTemp->event_id->id=eventID;
	    thisEventTemp = thisEventTemp->next;
	  }
	
	if ( cohInspFilterParams->cohSNROut )
	  {
	    snprintf( cohdataStr, LALNameLength*sizeof(CHAR),
			 "SNR_%Ld", eventID );
	    strcpy( cohInspFilterParams->cohSNRVec->name, "Coherent");
	    outFrameCoh = fr_add_proc_REAL4TimeSeries( outFrameCoh, cohInspFilterParams->cohSNRVec, "none", cohdataStr );
	  }
	
	/* save the coherent-snr of the H1-H2 pair */
	if ( cohInspFilterParams->cohH1H2SNROut )
	  {
	    snprintf( cohdataStr, LALNameLength*sizeof(CHAR),
			 "H1H2SNR_%Ld", eventID );
	    strcpy( cohInspFilterParams->cohH1H2SNRVec->name, "Coherent");
	    outFrameCohH1H2SNR = fr_add_proc_REAL4TimeSeries( outFrameCohH1H2SNR, cohInspFilterParams->cohH1H2SNRVec, "none", cohdataStr );
	  }
	
	/* save H1-H2 null-stream statistic in frames */
	if ( cohInspFilterParams->nullStatH1H2Out )
	  {
	    snprintf( cohdataStr, LALNameLength*sizeof(CHAR),
			 "H1H2_NullStat_%Ld", eventID );
	    strcpy( cohInspFilterParams->nullStatH1H2Vec->name, "Coherent");
	    outFrameNullStatH1H2 = fr_add_proc_REAL4TimeSeries( outFrameNullStatH1H2, cohInspFilterParams->nullStatH1H2Vec, "none", cohdataStr );
	  }
	
	/* save network null-stream statistic in frames */
	if ( cohInspFilterParams->nullStatOut )
	  {
	    snprintf( cohdataStr, LALNameLength*sizeof(CHAR),
			 "NullStat_%Ld", eventID );
	    strcpy( cohInspFilterParams->nullStatVec->name, "Coherent");
	    outFrameNullStat = fr_add_proc_REAL4TimeSeries( outFrameNullStat, cohInspFilterParams->nullStatVec, "none", cohdataStr );
	  }
	
	if ( !eventsOut )
	  {
	    while( thisEvent )
	      {
		MultiInspiralTable *tempEvent = thisEvent;
		thisEvent = thisEvent->next;
		LALFree( tempEvent->event_id );
		LALFree( tempEvent );
	      }
	  }
	
	if( thisEvent )
	  {
	    if( vrbflg ) fprintf( stdout,"******> Dumping Events <******\n");
	    if( !savedEvents.multiInspiralTable )
	      {
		MultiInspiralTable *tempEvent = thisEvent;
		tempTable = (MultiInspiralTable *) LALCalloc( 1, sizeof(MultiInspiralTable) );
		memcpy(tempTable,thisEvent,sizeof(MultiInspiralTable) );
		savedEvents.multiInspiralTable = tempTable;
		thisEvent = thisEvent->next;
		LALFree( tempEvent );
		tempEvent = NULL;
		if( thisEvent )
		  {
		    while( thisEvent )
		      {
			MultiInspiralTable *tempEvent = thisEvent;
			tempTable->next = (MultiInspiralTable *) LALCalloc( 1, sizeof(MultiInspiralTable) );
			tempTable = tempTable->next;
			memcpy(tempTable, thisEvent, sizeof(MultiInspiralTable) );
			thisEvent = thisEvent->next;
			
			LALFree( tempEvent );
			tempEvent = NULL;
		      }
		  }
	      }
	    else
	      {
		while( thisEvent )
		  {
		    MultiInspiralTable *tempEvent = thisEvent;
		    tempTable->next = (MultiInspiralTable *) LALCalloc( 1, sizeof(MultiInspiralTable) );
		    tempTable = tempTable->next;
		    memcpy(tempTable, thisEvent, sizeof(MultiInspiralTable) );
		    thisEvent = thisEvent->next;
		    LALFree( tempEvent );
		    tempEvent = NULL;
		  }
	      }
	    
	    while( thisEvent )
	      {
		MultiInspiralTable *tempEvent = thisEvent;
		thisEvent = thisEvent->next;
		LALFree( tempEvent );
		tempEvent = NULL;
	      }
	    
	  }/*end if( thisEvent ) */
	
	/* Finalize */
	LAL_CALL( LALCDestroyVector( &status, &(tempSnippet.data) ), &status );
	LALFree( cohInspFilterInput->tmplt );
	cohInspFilterInput->tmplt = NULL;
	
	for ( l = 0 ; l<(INT4)numDetectors ; l++ ) {
	  XLALDestroyCOMPLEX8TimeSeries( cohInspCVec->cData[l] );
	}
	
	/* Destroy params structure for coherent filter code */
	LAL_CALL( LALCoherentInspiralFilterParamsFinalize (&status,&cohInspFilterParams), &status );
	cohInspFilterParams = NULL;
	
	/* Destroy input structure for coherent filter code */
	LAL_CALL( LALCoherentInspiralFilterInputFinalize (&status, &cohInspFilterInput), &status);
	cohInspFilterInput = NULL;
	cohInspBeamVec = NULL;
	cohInspCVec = NULL;
	
	for(j=0;j<5;j++)
	  {
	    tempTime[j] = 0.0;
	    timeptDiff[j] = 0;
	  }
	tempTime[5] = 0.0;
	
	thisCoinc=thisCoinc->next;
	
      }/* Close loop over cohtrigs*/
      
      /* Write the summary information */
      if ( userTag )	{
	snprintf( fileName, FILENAME_MAX, "H1H2-CHIA_COHSNR_%d_%s-%d-%d", cohFileID,
		     userTag, gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
      }
      else	  {
	snprintf( fileName, FILENAME_MAX, "H1H2-CHIA_COHSNR_%d-%d-%d",
		     cohFileID, gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
      }
      
      if( outFrameCohH1H2SNR )
	{
	  if ( outputPath[0] )
	    {
	      snprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s/%s.gwf", outputPath, fileName);
	    }
	  else 
	    {
	      snprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s.gwf", fileName );
	    }
	  
	  if ( vrbflg ) fprintf( stdout, "writing H1-H2 coherent-snr frame data to %s....", framename );
	  frOutFile = FrFileONew( framename, 0);
	  FrameWrite( outFrameCohH1H2SNR, frOutFile);
	  FrFileOEnd( frOutFile );
	  if ( vrbflg ) fprintf(stdout, "done\n");
	  
	}
      
      if ( userTag )	{
	snprintf( fileName, FILENAME_MAX, "H1H2-CHIA_NULL_STAT_%d_%s-%d-%d", cohFileID,
		     userTag, gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
      }
      else	{
	snprintf( fileName, FILENAME_MAX, "H1H2-CHIA_NULL_STAT_%d-%d-%d",
		     cohFileID, gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
      }
      
      if( outFrameNullStatH1H2 )
	{
	  if ( outputPath[0] )
	    {
	      snprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s/%s.gwf", outputPath, fileName);
	    }
	  else 
	    {
	      snprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s.gwf", fileName );
	    }
	  
	  if ( vrbflg ) fprintf( stdout, "writing null statistic frame data to %s....", framename );
	  frOutFile = FrFileONew( framename, 0);
	  FrameWrite( outFrameNullStatH1H2, frOutFile);
	  FrFileOEnd( frOutFile );
	  if ( vrbflg ) fprintf(stdout, "done\n");
	  
	}
      
      if ( userTag )    {
	snprintf( fileName, FILENAME_MAX, "%s-CHIA_NULL_STAT_%d_%s-%d-%d", ifos, cohFileID,
		     userTag, gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
      }
      else      {
	snprintf( fileName, FILENAME_MAX, "%s-CHIA_NULL_STAT_%d-%d-%d", ifos,
		     cohFileID, gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
      }
      
      if( outFrameNullStat )
        {
          if ( outputPath[0] )
            {
              snprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s/%s.gwf", outputPath, fileName);
            }
          else
            {
              snprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s.gwf", fileName );
            }
	  
          if ( vrbflg ) fprintf( stdout, "writing null statistic frame data to %s....", framename );
          frOutFile = FrFileONew( framename, 0);
          FrameWrite( outFrameNullStat, frOutFile);
          FrFileOEnd( frOutFile );
          if ( vrbflg ) fprintf(stdout, "done\n");
	  
	}
      
      if ( userTag )	  {
	snprintf( fileName, FILENAME_MAX, "%s-CHIA_%d_%s-%d-%d", ifos, cohFileID,
		     userTag, gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
      }
      else	{
	snprintf( fileName, FILENAME_MAX, "%s-CHIA_%d-%d-%d", ifos,
		     cohFileID, gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
      }
      
      if( outFrameCoh )
	{
	  if ( outputPath[0] )
	    {
	      snprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s/%s.gwf", outputPath, fileName);
	    }
	  else 
	    {
	      snprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s.gwf", fileName );
	    }
	  
	  if ( vrbflg ) fprintf( stdout, "writing coherent frame data to %s....", framename );
	  frOutFile = FrFileONew( framename, 0);
	  FrameWrite( outFrameCoh, frOutFile);
	  FrFileOEnd( frOutFile );
	  if ( vrbflg ) fprintf(stdout, "done\n");
	  
	}
      
      if (eventsOut )
	{ 
	  memset( &results, 0, sizeof(LIGOLwXMLStream) );
	  if ( outputPath[0] )
	    {
	      if ( outCompress )
		{
		  snprintf( xmlname, FILENAME_MAX * sizeof(CHAR), "%s/%s.xml.gz", outputPath, fileName);
		}
	      else
		{
		  snprintf( xmlname, FILENAME_MAX * sizeof(CHAR), "%s/%s.xml", outputPath, fileName);
		}
	    }
	  else 
	    {
	      if ( outCompress )
		{
		  snprintf( xmlname, FILENAME_MAX * sizeof(CHAR), "%s.xml.gz", fileName );		
		}
	      else 
		{
		  snprintf( xmlname, FILENAME_MAX * sizeof(CHAR), "%s.xml", fileName );
		}	    
	    }
	  if ( vrbflg ) fprintf( stdout, "writing XML data to %s...\n", xmlname );
	  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, xmlname), &status );
	  
	  /* write the process table */
	  j=0;
	  while( !(caseID[j]) ) {
	    ifo = caseIDChars[j];
	    j++;
	  }
	  
	  if ( vrbflg ) fprintf( stdout, "  process table...\n" );
	  /*      snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", caseID );*/
	  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time), &accuracy ), &status );
	  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ), &status );
	  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable, process_table ), &status );
	  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
	  
	  /* write the process params table */
	  if ( vrbflg ) fprintf( stdout, "  process_params table...\n" );
	  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ), &status );
	  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams, process_params_table ), &status );
	  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
	  
	  /* write the search summary table */
	  if ( vrbflg ) fprintf( stdout, "  search_summary table...\n" );
	  searchsumm.searchSummaryTable->out_start_time.gpsSeconds = 
	    gpsStartTime.gpsSeconds + (numPointsSeg / (4 * sampleRate));
	  searchsumm.searchSummaryTable->out_end_time.gpsSeconds = 
	    gpsEndTime.gpsSeconds - (numPoints / (4 * sampleRate));
	  
	  /* the number of nodes for a standalone job is always 1 */
	  searchsumm.searchSummaryTable->nnodes = 1; 
	  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
					    search_summary_table ), &status );
	  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsumm, 
					    search_summary_table ), &status );
	  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

	  if ( numTriggers )
	    {
	      if ( vrbflg ) fprintf( stdout, "  search_summvars table...\n" );
	      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
						search_summvars_table ), &status );
	      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsummvars, 
						search_summvars_table ), &status );
	      LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
	    }
	  
	  /* store calibration information */
	  if ( cohFileID == 1) {
	    ADD_SUMM_VALUE( "calibration alpha", "analysis", alpha, 0 );
	    ADD_SUMM_VALUE( "calibration alphabeta", "analysis", alphabeta, 0 );
	    ADD_SUMM_VALUE( "calibration alpha", "injection", inj_alpha, 0 );
	    ADD_SUMM_VALUE( "calibration alphabeta", "injection", inj_alphabeta, 0 );
	  }
	  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, summ_value_table ), 
		    &status );
	  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, summvalue, 
					    summ_value_table ), &status );
	  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

	  /* write multi_inspiral table */
	  if ( savedEvents.multiInspiralTable ) { 
	    if( vrbflg ) fprintf(stdout,"  event params table\n ");
	    
	    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, multi_inspiral_table ), &status );
	    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, savedEvents, multi_inspiral_table ), &status );
	    LAL_CALL( LALEndLIGOLwXMLTable( &status, &results), &status );
	  }

	  while( savedEvents.multiInspiralTable )
	    {  
	      MultiInspiralTable *tempEvent2 = savedEvents.multiInspiralTable;
	      savedEvents.multiInspiralTable = savedEvents.multiInspiralTable->next;
	      LALFree( tempEvent2->event_id );
	      LALFree( tempEvent2 );
	    }
	  
	  /* close the output xml file */
	  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );
	  if ( vrbflg ) fprintf( stdout, "done. XML file closed\n" );
	  
	}/*end "if ( eventsOut )" */
      
    }/* close "for( cohFileID...)" */
    
  }/* closes "if ( numTriggers < 0 )" */
  
  free( proctable.processTable ); 
  while( procparams.processParamsTable )
    {
      this_proc_param = procparams.processParamsTable;
      procparams.processParamsTable = this_proc_param->next;
      free( this_proc_param );
    }
  
  while( searchsummvars.searchSummvarsTable )
    {
      this_search_summvar = searchsummvars.searchSummvarsTable;
      searchsummvars.searchSummvarsTable = this_search_summvar->next;
      LALFree( this_search_summvar );
    }
  free( searchsummvars.searchSummaryTable );
  
  while ( summvalue.summValueTable )
  {
    this_summ_value = summvalue.summValueTable;
    summvalue.summValueTable = summvalue.summValueTable->next;
    LALFree( this_summ_value );
  }
  /* free the search summary table after the summ_value table is written */
  free( searchsumm.searchSummaryTable );

  while ( cohbankEventList )  {
    currentTrigger = cohbankEventList;
    cohbankEventList = cohbankEventList->next;
    LAL_CALL( LALFreeSnglInspiral( &status, &currentTrigger ), &status );
  }
  
  while ( inputFiles ) {
    thisInputFile = inputFiles;
    inputFiles = thisInputFile->next;
    LALFree( thisInputFile );
  }
  
  while ( searchSummList )  {
    thisSearchSumm = searchSummList;
    searchSummList = searchSummList->next;
    LALFree( thisSearchSumm );
  }

  while ( coincHead )	{
    thisCoinc = coincHead;
    coincHead = coincHead->next;
    LALFree( thisCoinc );
  }

  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
 
  return 0;
 
}/* main function end */


/* ------------------------------------------------------------------------- */


#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
      snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
          long_options[option_index].name ); \
          snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
          snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );


#define USAGE1 \
"lalapps_inspiral [options]\n\n"\
"  --help                       display this message\n"\
"  --verbose                    print progress information\n"\
"  --version                    print version information and exit\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"  --low-frequency-cutoff F     low f cutoff of previously filtered data\n"\
"  --ifo-tag STRING             set STRING to whatever the ifo-tag of \n"\
                                "the bank file(needed for file naming) \n"\
"  --user-tag STRING            set STRING to tag the file names\n"\
"\n"
#define USAGE2 \
"  --bank-file FILE             read template bank parameters from FILE\n"\
"  --sample-rate N              set data sample rate to N\n"\
"  --segment-length N           set N to same value used in inspiral.c\n"\
"  --dynamic-range-exponent N   set N to same value used in inspiral.c\n"\
"  [--g1-slide]      g1_slide    Slide G1 data by multiples of g1_slide\n"\
"  [--h1-slide]      h1_slide    Slide H1 data by multiples of h1_slide\n"\
"  [--h2-slide]      h2_slide    Slide H2 data by multiples of h2_slide\n"\
"  [--l1-slide]      l1_slide    Slide L1 data by multiples of l1_slide\n"\
"  [--t1-slide]      t1_slide    Slide T1 data by multiples of t1_slide\n"\
"  [--v1-slide]      v1_slide    Slide V1 data by multiples of v1_slide\n"\
"  [--numCohTrigs]   numCohTrigs number of coherent snr frames per output frame-file\n"\
"  --cohsnr-threshold RHO       set signal-to-noise threshold to RHO\n"\
"  --maximize-over-chirp        do clustering\n"\
"  --gps-start-time SEC         GPS second of data start time (needed if globbing)\n"\
"  --gps-end-time SEC           GPS second of data end time (needed if globbing)\n"\
"  --ra-step         raStep     right-ascension step-size (in degrees)\n"\
"  --dec-step        decStep    declination step-size (in degrees)\n"\
"\n"
#define USAGE3 \
"  --write-events               write events\n"\
"  --write-cohsnr               write cohsnr\n"\
"  --write-cohnullstat          write coherent network null statistic \n"\
"  --write-h1h2nullstat         write H1-H2 null statistic \n"\
"  --write-cohh1h2snr           write H1-H2 coherent-snr when data from both ifos are present\n"\
"  --output-path                write files here\n"\
"  --write-compress             write compressed xml files\n"\
"  --H1-framefile               frame data for H1\n"\
"  --H2-framefile               frame data for H2\n"\
"  --L1-framefile                frame data for L\n"\
"  --V1-framefile                frame data for V\n"\
"  --G1-framefile                frame data for G\n"\
"  --T1-framefile                frame data for T\n"\
"\n"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
   struct option long_options[] = 
     {
     {"verbose",                  no_argument,       &vrbflg,            1 },
     {"write-compress",           no_argument,       &outCompress,       1 },
     {"help",                     no_argument,       0,                 'h'},
     {"version",                  no_argument,       0,                 'v'},
     {"debug-level",              required_argument, 0,                 'd'},
     {"ifo-tag",                  required_argument, 0,                 'I'},
     {"user-tag",                 required_argument, 0,                 'B'},
     {"low-frequency-cutoff",     required_argument, 0,                 'f'},
     {"bank-file",                required_argument, 0,                 'u'},
     {"sample-rate",              required_argument, 0,                 'r'},
     {"segment-length",           required_argument, 0,                 'l'},
     {"dynamic-range-exponent",   required_argument, 0,                 'e'}, 
     {"g1-slide",                 required_argument, 0,                 'g'},
     {"h1-slide",                 required_argument, 0,                 'W'},
     {"h2-slide",                 required_argument, 0,                 'X'},
     {"l1-slide",                 required_argument, 0,                 'Y'},
     {"t1-slide",                 required_argument, 0,                 't'},
     {"v1-slide",                 required_argument, 0,                 'w'},
     {"numCohTrigs",              required_argument, 0,                 'n'},
     {"cohsnr-threshold",         required_argument, 0,                 'p'},
     {"maximize-over-chirp",      no_argument,       &maximizeOverChirp, 1 },
     {"write-events",             no_argument,       &eventsOut,         1 },
     {"write-cohsnr",             no_argument,       &cohSNROut,         1 },
     {"write-cohnullstat",        no_argument,       &nullStatOut,       1 },
     {"write-h1h2nullstat",       no_argument,       &nullStatH1H2Out,   1 },
     {"write-cohh1h2snr",         no_argument,       &cohH1H2SNROut,     1 },
     {"gps-start-time",           required_argument, 0,                 'a'},
     {"gps-end-time",             required_argument, 0,                 'b'},
     {"ra-step",                  required_argument, 0,                 'R'},
     {"dec-step",                 required_argument, 0,                 'D'},
     {"output-path",              required_argument, 0,                 'P'},
     {"H1-framefile",             required_argument, 0,                 'A'},
     {"H2-framefile",             required_argument, 0,                 'Z'},
     {"L1-framefile",             required_argument, 0,                 'L'},
     {"V1-framefile",             required_argument, 0,                 'V'},
     {"G1-framefile",             required_argument, 0,                 'G'},
     {"T1-framefile",             required_argument, 0,                 'T'},
     {0, 0, 0, 0}
   };

   int c;
   ProcessParamsTable *this_proc_param = procparams.processParamsTable;

   while (1)
     {
       /* getopt_long stores long options here */
       int option_index = 0;
       size_t optarg_len;

       c = getopt_long_only( argc, argv,
           "A:B:a:b:D:G:I:L:l:e:g:W:X:Y:t:w:n:P:R:T:V:Z:d:f:h:p:r:u:v:",
           long_options, &option_index );

       if ( c == -1 )
         {
           break;
         }

       switch ( c )
         {
         case 0:
        /* if this option set a flag, do nothing else now */
           if ( long_options[option_index].flag != 0 )
             {
               break;        }
           else 
             {
               fprintf( stderr, "error parsing option %s with argument %s\n",
                        long_options[option_index].name, optarg );
               exit( 1 );
             }
           break;

         case 'A':
	   optarg_len = strlen( optarg ) + 1;
	   ifoframefile[1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	   memcpy( ifoframefile[1], optarg, optarg_len );
           H1file = 1;
	   ADD_PROCESS_PARAM( "string", "%s", optarg );
	   break;

         case 'D': /* set right-ascension step-size (in degrees) */
           decStep = atof (optarg);
           ADD_PROCESS_PARAM( "float", "%e", decStep );
           break;

         case 'R': /* set right-ascension step-size (in degrees) */
           raStep = atof (optarg);
           ADD_PROCESS_PARAM( "float", "%e", raStep );
           break;

         case 'G':
	   optarg_len = strlen( optarg ) + 1;
	   ifoframefile[0] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	   memcpy( ifoframefile[0], optarg, optarg_len );
           G1file = 1;
	   ADD_PROCESS_PARAM( "string", "%s", optarg );
           break;

         case 'L':
	   optarg_len = strlen( optarg ) + 1;
	   ifoframefile[3] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	   memcpy( ifoframefile[3], optarg, optarg_len );
           L1file = 1;
	   ADD_PROCESS_PARAM( "string", "%s", optarg );
           break;

         case 'T':
	   optarg_len = strlen( optarg ) + 1;
	   ifoframefile[4] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	   memcpy( ifoframefile[4], optarg, optarg_len );
           T1file = 1;
	   ADD_PROCESS_PARAM( "string", "%s", optarg );
           break;

         case 'V':
	   optarg_len = strlen( optarg ) + 1;
	   ifoframefile[5] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	   memcpy( ifoframefile[5], optarg, optarg_len );
           V1file = 1;
	   ADD_PROCESS_PARAM( "string", "%s", optarg );
           break;

         case 'Z':
	   optarg_len = strlen( optarg ) + 1;
	   ifoframefile[2] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	   memcpy( ifoframefile[2], optarg, optarg_len );
           H2file = 1;
	   ADD_PROCESS_PARAM( "string", "%s", optarg );
           break;

         case 'B':
           /* create storage for the user-tag */
           optarg_len = strlen( optarg ) + 1;
           userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
           memcpy( userTag, optarg, optarg_len );
           ADD_PROCESS_PARAM( "string", "%s", optarg );
           break;

         case 'I':
           /* create storage for the ifo-tag */
           optarg_len = strlen( optarg ) + 1;
           ifos = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
           memcpy( ifos, optarg, optarg_len );
           ADD_PROCESS_PARAM( "string", "%s", optarg );
           break;

         case 'P':
           memset( outputPath, 0, FILENAME_MAX * sizeof(CHAR) );
           snprintf( outputPath, FILENAME_MAX * sizeof(CHAR),"%s", optarg );
           ADD_PROCESS_PARAM( "string", "%s", outputPath );
           break;
	   
         case 'd': /* set debuglevel */
           set_debug_level( optarg );
           ADD_PROCESS_PARAM( "string", "%s", optarg );
           break;
	   
         case 'f': /* set fLow */
           fLow = (REAL4) atof (optarg);
           ADD_PROCESS_PARAM( "float", "%e", fLow );
           break;

         case 'h':
           fprintf( stdout, USAGE1 );
           fprintf( stdout, USAGE2 );
           fprintf( stdout, USAGE3 );
           exit( 0 );
           break;

         case 'p': /* set coherent SNR threshold */
           cohSNRThresh = atof (optarg);
           ADD_PROCESS_PARAM( "float", "%e", cohSNRThresh );
           break;

         case 'l':
           numPointsSeg = (INT4) atoi(optarg);
           ADD_PROCESS_PARAM("int", "%d", numPointsSeg );
           break;

         case 'e':
           dynRangeExponent = (REAL4) atof(optarg);
           ADD_PROCESS_PARAM("float", "%e", dynRangeExponent );
           break;

         case 'r': 
           sampleRate = (INT4) atoi(optarg);
           ADD_PROCESS_PARAM("int", "%d", sampleRate );
           break;

         case 'u':
           /* create storage for the bank filename */
           /*optarg_len = strlen( optarg ) + 1;
           bankFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
           memcpy( bankFileName, optarg, optarg_len );*/
	   strcpy(cohbankFileName, optarg);
	   char tempName[256];
	   char *duration =NULL;
	   strcpy(tempName, cohbankFileName);
	   duration = strtok(tempName,"-");
	   duration = strtok(NULL,"-");
	   duration = strtok(NULL,"-");
	   duration = strtok(NULL,".");
	   bankDuration=atoi(duration);        
           ADD_PROCESS_PARAM( "string", "%s", cohbankFileName );
	   duration=NULL;
       break;

           /* Read in time-slide steps for all detectors */
           /* Read in time-slide step for G1 */
         case 'g': 
           slideStep[0] = (REAL8) atof(optarg);
           ADD_PROCESS_PARAM("float", "%e", slideStep[0] );
           break;

           /* Read in time-slide step for H1 */
         case 'W': 
           slideStep[1] = (REAL8) atof(optarg);
           ADD_PROCESS_PARAM("float", "%e", slideStep[1]);
           break;

           /* Read in time-slide step for H2 */
         case 'X': 
           slideStep[2] = (REAL8) atof(optarg);
           ADD_PROCESS_PARAM("float", "%e", slideStep[2]);
           break;

           /* Read in time-slide step for L1 */
         case 'Y': 
           slideStep[3] = (REAL8) atof(optarg);
           ADD_PROCESS_PARAM("float", "%e", slideStep[3]);
           break;

           /* Read in time-slide step for T1 */
         case 't': 
           slideStep[4] = (REAL8) atof(optarg);
           ADD_PROCESS_PARAM("float", "%e", slideStep[4]);
           break;

           /* Read in time-slide step for V1 */
         case 'w': 
           slideStep[5] = (REAL8) atof(optarg);
           ADD_PROCESS_PARAM("float", "%e", slideStep[5]);
           break;

           /* Read in time-slide step for V1 */
         case 'n': 
           numCohTrigs = (INT4) atol(optarg);
           ADD_PROCESS_PARAM("int", "%d", numCohTrigs );
           break;

         case 'v':
           /* print version information and exit */
           fprintf( stdout, "LIGO/LSC Multi-Detector Search Code\n" 
                 "Bose/Seader <sukanta@wsu.edu>\n"
                 "CVS Version: " CVS_ID_STRING "\n"
                 "CVS Tag: " CVS_NAME_STRING "\n" );
           fprintf( stdout, lalappsGitID );
           exit( 0 );
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
	     gpsStartTimeTemp=gstartt;
             ADD_PROCESS_PARAM( "int", "%ld", gstartt );
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
	     gpsEndTimeTemp=gendt;
             ADD_PROCESS_PARAM( "int", "%ld", gendt );
           }
           break;

         case '?':
           exit( 1 );
           break;

         default:
           fprintf( stderr, "unknown error while parsing options\n" );
           exit( 1 );

         }

     }

   if (optind < argc)
     {
       fprintf( stderr, "extraneous command line arguments:\n" );
       while ( optind < argc )
         {
           fprintf ( stderr, "%s\n", argv[optind++] );
         }
       exit( 1 );      
     }

   /*  Store optional arguments in the process param table */
   if ( eventsOut == 1 )
     {
       snprintf( procparams.processParamsTable->program, 
		    LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
       snprintf( procparams.processParamsTable->param,
		    LIGOMETA_PARAM_MAX, "--write-events" );
       snprintf( procparams.processParamsTable->type, 
		    LIGOMETA_TYPE_MAX, "string" );
       snprintf( procparams.processParamsTable->value, 
		    LIGOMETA_VALUE_MAX, " " );
     }
   
   if ( cohSNROut == 1 )
     {
       this_proc_param = this_proc_param->next = (ProcessParamsTable *)
	 calloc( 1, sizeof(ProcessParamsTable) );
       snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
		    "%s", PROGRAM_NAME );
       snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
		    "--write-cohsnr" );
       snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
       snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
     }
   
   if ( cohH1H2SNROut == 1 )
     {    
       this_proc_param = this_proc_param->next = (ProcessParamsTable *)
	 calloc( 1, sizeof(ProcessParamsTable) );
       snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
		    "%s", PROGRAM_NAME );
       snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
		    "--write-cohh1h2snr" );
       snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
       snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
     }
   
   if ( nullStatOut == 1 )
     {
       this_proc_param = this_proc_param->next = (ProcessParamsTable *)
	 calloc( 1, sizeof(ProcessParamsTable) );
       snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
		    "%s", PROGRAM_NAME );
       snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
		    "--write-cohnullstat" );
       snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
       snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
     }

   if ( nullStatH1H2Out == 1 )
     {
       this_proc_param = this_proc_param->next = (ProcessParamsTable *)
         calloc( 1, sizeof(ProcessParamsTable) );
       snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
                    "%s", PROGRAM_NAME );
       snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
                    "--write-h1h2nullstat" );
       snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
       snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
     }

   /* check validity of input data time if globbing */
   /* the times should be that spanned by the bank(trigger) file */
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

   /* check sample rate has been given */
   if ( sampleRate < 0 )
     {
       fprintf( stderr, "--sample-rate must be specified\n" );
       exit( 1 );
     }

   if ( numPointsSeg < 0 )
     {
       fprintf( stderr, "--segment-length must be specified.\n" );
       fprintf( stderr,"It must be set to the same value as was used in inspiral.c when the C-data was generated.\n");
       exit( 1 );
     }

   if ( dynRangeExponent < 0 )
     {
       fprintf( stderr, "--dynamic-range-exponent must be specified.\n" );
       fprintf( stderr,"It must be set to the same value as was used in inspiral.c when the C-data was generated.\n");
       exit( 1 );
     }

   if ( ! cohbankFileName )
     {
       fprintf( stderr, "--bank-file must be specified\n" );
       exit( 1 );
     }

   if ( fLow < 0 )
     {
       fprintf( stderr, "--low-frequency-cutoff must be specified\n" );
       exit( 1 );
     }

   if ( cohSNRThresh < 0 )
     {
       fprintf( stderr, "--cohsnr-threshold must be specified\n" );
       exit( 1 );
     }

   if( !ifos )
     {
       fprintf(stderr, "--ifo-tag must be specified for file naming\n" );
       exit( 1 );
     }

   return 0;
}

#undef ADD_PROCESS_PARAM
