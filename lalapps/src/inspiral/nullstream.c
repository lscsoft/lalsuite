/*
*  Copyright (C) 2007 Eirini Messaritaki
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
 * File Name: nullstream.c
 *
 * Author: Messaritaki, E.
 *
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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <math.h>

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
#include <lal/StochasticCrossCorrelation.h>
#include <lal/DetectorSite.h>
#include <lal/Random.h>
#include <lal/LALInspiral.h>
#include <lal/CoherentInspiral.h>
#include <lal/NullStatistic.h>
#include <lal/LALStatusMacros.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALFrameL.h>

#include <lal/lalGitID.h>
#include <lalappsGitID.h>

#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "nullstream"
#define CVS_NAME_STRING "$Name$"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* input data parameters */

INT4   sampleRate           = -1;  /* sample rate of filter data   */
INT4   numPointsSeg         = -1;  /* set to segment-length from inspiral.c */
REAL4  dynRangeExponent     = -1;  /* set to same value used in inspiral.c */

/*null stream specific inputs*/

char   ifoframefile[LAL_NUM_IFO][256];

INT4 H1file = 0;
INT4 H2file = 0;
INT4 L1file = 0;
INT4 G1file = 0;
INT4 T1file = 0;
INT4 V1file = 0;

/* input time-slide parameters */
REAL8  slideStep[LAL_NUM_IFO]        = {0.0,0.0,0.0,0.0,0.0,0.0};
int    bankDuration        = 0;
CHAR  *cohbankFileName     = NULL;   /* name of input template bank  */
CHAR  *xmlFileNameH1       = NULL;   /* name of H1 inspiral trigger file */
CHAR  *xmlFileNameH2       = NULL;   /* name of H2 inspiral trigger file */
int    nullStatOut         = 0;      /* default is not to write frame */
int    eventsOut           = 0;      /* default is not to write events */
REAL4  nullStatThresh      = -1;
INT4   maximizeOverChirp   = 0;      /* default is no clustering */
CHAR   outputPath[FILENAME_MAX];
CHAR  *frInType            = NULL;   /* type of data frames */

INT8          gpsStartTimeNS     = 0;   /* input data GPS start time ns */
LIGOTimeGPS   gpsStartTime;             /* input data GPS start time    */
INT8          gpsEndTimeNS       = 0;   /* input data GPS end time ns   */
LIGOTimeGPS   gpsEndTime;               /* input data GPS end time      */
int           gpsStartTimeTemp   = 0;   /* input data GPS start time ns */
int           gpsEndTimeTemp     = 0;   /* input data GPS start time ns */

LALStatus             status;

CHAR  *userTag          = NULL;         /* string the user can tag with */
CHAR  *ifoTag           = NULL;         /* string to tag IFOs    */


int main( int argc, char *argv[] )
{
  /* FrChanIn      frChan; */

  /* frame data */
  FrStream       *frStream   = NULL;
  struct FrFile  *frOutFile  = NULL;
  struct FrameH  *outFrame   = NULL;

  /* output */
  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;

  /* FILE *filePtr[4]; */

  CHAR  fileName[FILENAME_MAX];
  CHAR  framename[FILENAME_MAX];
  CHAR  xmlname[FILENAME_MAX];
  CHAR  nullStatStr[LALNameLength];

  INT4   segLength      = 4;  /* should match hardcoded value in inspiral.c */
  INT4   numPoints      = 0;
  UINT4  numSegments    = 1;  /* number of segments */
  UINT4  numNullStatFr  = 0;  
  UINT8  eventID        = 0;

  REAL4  m1             = 0.0;
  REAL4  m2             = 0.0;
  REAL4  dynRange       = 0.0;

  /* variables for initializing tempTime to account for time-slides */
  /*UINT8  triggerNumber  = 0;
  UINT8  slideNumber    = 0;
  UINT8  slideSign      = 0; */

  /* counters and other variables */
  INT4   j, k, l;
  INT4   kidx            = 0;
  UINT4  numDetectors    = 0;
  REAL8  tempTime[LAL_NUM_IFO]     = {0.0,0.0,0.0,0.0,0.0,0.0}; 
  INT4   numTriggers     = 0;
  INT4   numTrigsH1      = 0;
  INT4   numTrigsH2      = 0;
  INT4   numCoincs       = 0;
  INT4   numEvents       = 0;

  FrCache              *frInCache        = NULL;

  SnglInspiralTable    *currentTrigger   = NULL;
  SnglInspiralTable    *currentTrigH1    = NULL;
  SnglInspiralTable    *currentTrigH2    = NULL;
  SnglInspiralTable    *cohbankEventList = NULL;
  SnglInspiralTable    *trigListH1       = NULL;
  SnglInspiralTable    *trigListH2       = NULL;

  CoincInspiralTable   *coincHead = NULL;
  CoincInspiralTable   *thisCoinc = NULL;

  SearchSummvarsTable  *inputFiles       = NULL;
  SearchSummvarsTable  *inputFilesH1     = NULL;
  SearchSummvarsTable  *inputFilesH2     = NULL;
  SearchSummaryTable   *searchSummList   = NULL;
  SearchSummaryTable   *searchSummListH1 = NULL;
  SearchSummaryTable   *searchSummListH2 = NULL;


  NullStatInitParams      *nullStatInitParams   = NULL;
  NullStatParams          *nullStatParams       = NULL;
  NullStatInputParams     *nullStatInputParams  = NULL;
  CVector                 *CVec                 = NULL;
  MultiInspiralTable      *event                = NULL;
  MultiInspiralTable      *thisEvent            = NULL;
  MetadataTable            savedEvents;

  /* cData channel names */
  ChanNames  *cDataChanNames ; 

  if ( vrbflg ) fprintf( stdout, "%d.\n", LAL_NUM_IFO );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) LALCalloc(1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  if (strcmp(CVS_REVISION,"$Revi" "sion$"))
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
    LALCalloc( 1, sizeof(ProcessParamsTable) );

  arg_parse_check( argc, argv, procparams );
  if ( vrbflg )  fprintf(stdout, "Called parse options.\n");

  /* wind to the end of the process params table */
  for ( this_proc_param = procparams.processParamsTable; this_proc_param->next;
        this_proc_param = this_proc_param->next );


  /* set the dynamic range */
  dynRange = pow( 2.0, dynRangeExponent );

  /* set other variables */
  numPoints = sampleRate * segLength;
  savedEvents.multiInspiralTable = NULL;
  k = 0;

  /* read in the frame files: why is this here? it shows up later on too */
#if 0
  if ( vrbflg ) fprintf(stdout, "Reading in the frame files.\n");
  for ( k=0; k<LAL_NUM_IFO ; k++)
  {
    if ( (k == LAL_IFO_H1) || (k == LAL_IFO_H2) )
    {
      if ( ifoframefile[k] )
      { 
        frStream = XLALFrOpen( NULL, ifoframefile[k]); 
        if (!frStream)
        {
          fprintf(stdout,"The file %s does not exist - exiting.\n", 
                  ifoframefile[k]);
          goto cleanexit;
        }
      }   
    }
    else
    {
      fprintf( stdout, 
        "Frame file not needed and not read for interferometer %d.\n",k);
    }
  }
#endif


  /* read in the cohbank trigger ligo lw xml file */
  /* cohbankEventList is the list of the events in the COHBANK xml file */
  /* currentTrigger is the last trigger */
  numTriggers = XLALReadInspiralTriggerFile( &cohbankEventList, 
       &currentTrigger, &searchSummList, &inputFiles, cohbankFileName );
  fprintf(stdout,"Reading templates from %s.\n",cohbankFileName);

  numTrigsH1 = XLALReadInspiralTriggerFile( &trigListH1, &currentTrigH1,
        &searchSummListH1, &inputFilesH1, xmlFileNameH1);
  fprintf(stdout,"Read %d H1 triggers from %s.\n", numTrigsH1, xmlFileNameH1);

  numTrigsH2 = XLALReadInspiralTriggerFile( &trigListH2, &currentTrigH2,
        &searchSummListH2, &inputFilesH2, xmlFileNameH2);
  fprintf(stdout,"Read %d H2 triggers from %s.\n", numTrigsH2, xmlFileNameH2);

  if ( numTriggers < 0 )  /* no triggers found */
  {
    fprintf(stderr, "Error reading triggers from file %s.\n", cohbankFileName);
    exit( 1 );
  }
  else if ( numTriggers == 0 )  /* no triggers found */
  {
    if ( vrbflg )
    {
      fprintf( stdout, "%s contains no triggers - the cohbank will be empty.\n",
         cohbankFileName );
    }
  }
  else  /* triggers do exist */
  {
    if ( vrbflg )
    {
      fprintf( stdout, "Read in %d triggers from the file %s.\n", numTriggers,
               cohbankFileName );
    }

    /* pair up the coincidences */
    /* coincHead points to a CoincInspiralTable that was generated
       by the triggers in cohbankEventList */
    numCoincs = XLALRecreateCoincFromSngls( &coincHead, cohbankEventList );
    if ( numCoincs < 0 )
    {
      fprintf(stderr, "Unable to reconstruct coincs from single ifo triggers.");
      exit( 1 );
    }
    else if ( vrbflg )
    {
      fprintf( stdout, "Recreated %d coincs from the %d triggers.\n", 
               numCoincs, numTriggers );
    }

    cDataChanNames = (ChanNames *) LALMalloc( sizeof(ChanNames) );

    if ( !cDataChanNames )
    {
        fprintf( stdout, "Could not allocate memory for channel names.\n" );
        goto cleanexit;
    }




    /* loop over coincident triggers to compute the null statistic */
    for ( thisCoinc=coincHead; thisCoinc; thisCoinc=thisCoinc->next)
    {
      /* numDetectors = thisCoinc->numIfos;  detector number for this coinc */
      numDetectors = 2; /* hardcoded for now, change to previous line later */

      /* l is another detector index */
      l=0;

      /* Note the participating ifos and the eventID for this coincidence  */
      /* channel naming is trivial in the H1-H2 case; to be revisited when */
      /* more ifos are added; in that case use LAL_NUM_IFO and k instead of*/
      /* LAL_IFO_H1.                                                       */

      if ( thisCoinc->snglInspiral[LAL_IFO_H1] && 
           thisCoinc->snglInspiral[LAL_IFO_H2] )
      {
        /* record the eventid for this coincidence */
        eventID = thisCoinc->snglInspiral[LAL_IFO_H1]->event_id->id; 
        if ( vrbflg ) fprintf( stdout,"eventID = %Ld.\n",eventID );

        /* Parse eventID to get the slide number */
#if 0
        triggerNumber = eventID % 100000;
        slideNumber = ((eventID % 100000000) - triggerNumber)/100000;
        slideSign = (eventID % 1000000000)-(slideNumber*100000)-triggerNumber;
#endif

        if ( vrbflg ) fprintf( stdout, "Input Frame = %s.\n", frInType);

        /* Store CData frame name  */

        snprintf( cDataChanNames->chanNameH1,
          LALNameLength, "%s:%s_CData_%Ld",
          &thisCoinc->snglInspiral[LAL_IFO_H1]->ifo, frInType, eventID );
        if (vrbflg) fprintf( stdout, "H1 channel: %s\n",
           cDataChanNames->chanNameH1 );

        snprintf( cDataChanNames->chanNameH2,      
          LALNameLength, "%s:%s_CData_%Ld",             
          (*(thisCoinc->snglInspiral[LAL_IFO_H2])).ifo, frInType, eventID );
        if (vrbflg) fprintf( stdout, "H2 channel: %s\n",
           cDataChanNames->chanNameH2 );

        kidx =  LAL_IFO_H1;
      }

#if 0
      /* Initialize tempTime to account for time-slides */
      for( j=0; j<5; j++)
      {
        /* slideSign=0 is the same as a positive time slide */
        if(slideSign != 0)
        {
          tempTime[j] = slideStep[j]*slideNumber*slideSign;
        }
        else
        {
          tempTime[j] -= slideStep[j]*slideNumber*slideSign;
        }
      }
#endif

      l=0;
      if ( G1file ) l++; 
      if ( H1file ) l++;
      if ( H2file ) l++;
      if ( L1file ) l++;
      if ( T1file ) l++;
      if ( V1file ) l++;


      if ( (INT4)numDetectors != l )
      {
        fprintf( stderr, 
         "You have events for %d detectors but frame files for %d detectors.\n",
         numDetectors,l);
        if ( (INT4)numDetectors > l )
        {
          fprintf( stderr, "Too few frame files specified. Exiting.\n");
          exit(1);
        }
        else
        {
          fprintf( stderr, "Too many frame files specified. Exiting.\n");
          exit(1);
        }
      }

      l = 0;

      if ( vrbflg ) fprintf(stdout,"numDetectors = %d\n", numDetectors);

      /* Initialize the necessary structures for thisCoinc-ident trigger*/

      if ( !(nullStatInitParams = (NullStatInitParams *) 
          LALCalloc(1,sizeof(NullStatInitParams)) ))
      {
        fprintf( stdout, 
         "Could not allocate memory for nullStat init params.\n" );
        goto cleanexit;
      }

      /* Initialize the null param structure for thisCoinc trigger */ 
      nullStatInitParams->numDetectors    = numDetectors;
      nullStatInitParams->numSegments     = numSegments;
      nullStatInitParams->numPoints       = numPoints;
      nullStatInitParams->nullStatOut     = nullStatOut;

      /* create the data structures needed */

      if ( vrbflg ) fprintf( stdout, "Initializing.\n " );

      /* initialize null statistic functions */
      XLALNullStatisticInputInit(&nullStatInputParams, nullStatInitParams);

      m1 = thisCoinc->snglInspiral[kidx]->mass1;
      m2 = thisCoinc->snglInspiral[kidx]->mass2;

      nullStatInputParams->tmplt = (InspiralTemplate *) 
         LALCalloc(1,sizeof(InspiralTemplate) );
      nullStatInputParams->tmplt->mass1 = m1;
      nullStatInputParams->tmplt->mass2 = m2;
      nullStatInputParams->tmplt->totalMass = m1 + m2;
      nullStatInputParams->tmplt->mu = m1 * m2 / (m1 + m2);
      nullStatInputParams->tmplt->eta = (m1 * m2) / ((m1 + m2) * (m1 + m2 ));

      if ( vrbflg ) fprintf(stdout,"m1:%f, m2:%f, Mtotal:%f, mu:%f, eta:%f.\n", 
         m1, m2, nullStatInputParams->tmplt->totalMass,
         nullStatInputParams->tmplt->mu, nullStatInputParams->tmplt->eta);

      XLALNullStatisticParamsInit(&nullStatParams, nullStatInitParams);
      if ( vrbflg ) fprintf( stdout, "Initializing.\n " );

      nullStatParams->numTmplts         = 1;
      nullStatParams->maxOverChirp      = maximizeOverChirp;
      nullStatParams->nullStatThresh    = nullStatThresh;
      nullStatParams->nullStatOut       = nullStatOut;

      /*
       * we only need H1 and H2 at the moment, but assign
       * all, for future use
       */
#if 0
      nullStatParams->detVector->detector[LAL_IFO_G1] = 
        lalCachedDetectors[LALDetectorIndexGEO600DIFF];
      nullStatParams->detVector->detector[LAL_IFO_H1] =
        lalCachedDetectors[0];
      nullStatParams->detVector->detector[LAL_IFO_H2] =
        lalCachedDetectors[0];
      nullStatParams->detVector->detector[LAL_IFO_L1] =
        lalCachedDetectors[1];
      nullStatParams->detVector->detector[LAL_IFO_T1] =
        lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
      nullStatParams->detVector->detector[LAL_IFO_V1] =
        lalCachedDetectors[LALDetectorIndexVIRGODIFF];
#endif

      CVec = nullStatInputParams->CData;

      /* Read in the snippets associated with thisCoinc trigger */
      for ( j=1; j<LAL_NUM_IFO; j++ )
      {
        if ( (j == LAL_IFO_H1) || (j == LAL_IFO_H2) ) 
        { 
          if (vrbflg) fprintf(stdout, " j = %d \n", j );
          frStream = XLALFrOpen( NULL, ifoframefile[j] ); 
          if ( vrbflg ) fprintf( stdout, 
                 "Getting the c-data time series for %s.\n",
                 thisCoinc->snglInspiral[j]->ifo );

          if (!frStream)
          {  
            fprintf(stdout,
              "The file %s does not exist - exiting.\n", ifoframefile[j] );
            goto cleanexit;
          }

          if ( vrbflg ) fprintf( stdout, "error");

          /* assign the channel names */
          if ( j == LAL_IFO_H1 ) 
          {
            strcpy( CVec->cData[j]->name, &(cDataChanNames->chanNameH1) ); 
            if ( vrbflg ) fprintf( stdout, "H1 channel name: %s \n",
                CVec->cData[j]->name );
            nullStatParams->sigmasq[j] = trigListH1->sigmasq;
            if ( vrbflg ) fprintf( stdout, "sigma-sq for H1: %f \n",
                 nullStatParams->sigmasq[j]);
          }


          if ( j == LAL_IFO_H2 )               
          {
            strcpy( CVec->cData[j]->name, &(cDataChanNames->chanNameH2) );
            if ( vrbflg ) fprintf( stdout, "H2 channel name: %s \n",
                CVec->cData[j]->name );
            nullStatParams->sigmasq[j] = trigListH2->sigmasq;
            if ( vrbflg ) fprintf( stdout, "sigma-sq for H2: %f \n",
                 nullStatParams->sigmasq[j]);
          }

          if ( vrbflg ) fprintf( stdout, "error");
          XLALFrGetCOMPLEX8TimeSeries( CVec->cData[j], frStream );
          if ( vrbflg ) fprintf( stdout, "error");

          /* Need to worry about WRAPPING of time-slides             */
          /* tempTime is the start time of cData plus - (time slide) */
          tempTime[j] += CVec->cData[j]->epoch.gpsSeconds + 
                         CVec->cData[j]->epoch.gpsNanoSeconds * 1e-9;
          if ( vrbflg ) fprintf( stdout,"tempTime = %f\n", tempTime[j] );
 
          XLALFrClose( frStream );
 
          if (j == 2) j = LAL_NUM_IFO+1;
        }
        else
        {
          if ( vrbflg ) fprintf( stdout, "No data needed for %s.\n", 
                 thisCoinc->snglInspiral[j]->ifo );
        }
      }      /* closes for( j=0; j<LAL_NUM_IFO; j++ ) */


      /*
       * skipping the commensuration of the c-data snippets. This
       * is necessary in principle, but the H1-H2 snippets are always
       * commensurate so we postpone writing that part of the code.
       */

      /* calculation of the null statistic for this coincident event */
      XLALComputeNullStatistic(&thisEvent, nullStatInputParams, nullStatParams);
      if (vrbflg) fprintf( stdout, "error \n" );

      if ( nullStatOut )
      {
        snprintf( nullStatStr, LALNameLength, "NULL_STAT_%d", 
                     numNullStatFr++ );
        strcpy( nullStatParams->nullStatVec->name, "NullStatistic");
        outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
                     nullStatParams->nullStatVec, "none", nullStatStr);
      }

      if ( !eventsOut )
      {
        while( thisEvent )
        {
          MultiInspiralTable *tempEvent = thisEvent;
          thisEvent = thisEvent->next;
          LALFree( tempEvent );
        }
      }
        

      /*  test if any events got returned                             */
      /*  not necessary now, but will be when a threshold is applied  */
      if ( thisEvent )
      {
        if ( vrbflg ) fprintf( stdout, "***>  dumping events  <***\n" );

        if ( ! savedEvents.multiInspiralTable )
        {
          savedEvents.multiInspiralTable = thisEvent;
        }
        else
        {
          thisEvent->next = thisEvent;
        }

        /* save a ptr to the last event in the list , count the events */
        ++numEvents;
        while ( thisEvent->next )
        {
          thisEvent = thisEvent->next;
          ++numEvents;
        }
        event = thisEvent;
        thisEvent = NULL;

      } /* close  if ( thisEvent ) */

      XLALNullStatisticParamsFinal( &nullStatParams );

      XLALNullStatisticInputFinal( &nullStatInputParams );

    } /* close for (thisCoinc=coincHead; thisCoinc; thisCoinc=thisCoinc->next */

  } /* close  else (if triggers do exist) */

  
  /* write the output xml file */
  if ( userTag )
  {
    snprintf( fileName, FILENAME_MAX, "%s-NULLSTAT_%s-%d-%d", ifoTag, 
                 userTag, gpsStartTime.gpsSeconds, 
                 gpsEndTime.gpsSeconds-gpsStartTime.gpsSeconds );
  }
  else
  {
    snprintf( fileName, FILENAME_MAX, "%s-NULLSTAT-%d-%d", ifoTag,
       gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds-gpsStartTime.gpsSeconds );
  }

  if ( nullStatOut )
  {
    if ( outputPath[0] )
    {
      snprintf( framename, FILENAME_MAX, "%s/%s.gwf",
                   outputPath, fileName );
    }
    else
    {
      snprintf( framename, FILENAME_MAX, "%s.gwf", fileName);
    }

    if ( vrbflg ) fprintf( stdout, "Writing null statistic time series to %s.",
                            framename );
    frOutFile = FrFileONew( framename, 0 );
    FrameWrite( outFrame, frOutFile );
    FrFileOEnd( frOutFile );
  }

  if ( eventsOut )
  {
    memset( &results, 0, sizeof(LIGOLwXMLStream) );
    if ( outputPath[0] )
    {
      snprintf( xmlname, FILENAME_MAX, "%s/%s.xml", 
                   outputPath, fileName );
    }
    else
    {
      snprintf( xmlname, FILENAME_MAX, "%s.xml", fileName );
    }
    if ( vrbflg ) fprintf( stdout, "Writing xml data to %s.", xmlname );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, xmlname ), &status );

    /* process table */
    if ( vrbflg ) fprintf( stdout, "Writing the process table..." );
    XLALGPSTimeNow(&(proctable.processTable->end_time));
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ), 
                                      &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable, 
                                      process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
    free( proctable.processTable );
    if ( vrbflg ) fprintf( stdout, " done.\n" );

    /* process params table */
    if ( vrbflg ) fprintf( stdout, "Writing the process_params table..." );
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
                                      process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams, 
                                      process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
    while( procparams.processParamsTable )
    {
      this_proc_param = procparams.processParamsTable;
      procparams.processParamsTable = this_proc_param->next;
      free( this_proc_param );
    }
    if ( vrbflg ) fprintf( stdout, " done.\n" );

    /* multi inspiral table */
    if( vrbflg ) fprintf( stdout,"Writing the multi-inspiral table...");
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
                                      multi_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, savedEvents, 
                                      multi_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &results), &status );

    while( savedEvents.multiInspiralTable )
    {
      MultiInspiralTable * tempEvent = savedEvents.multiInspiralTable;
      savedEvents.multiInspiralTable = savedEvents.multiInspiralTable->next;
      LALFree( tempEvent->event_id );
      LALFree( tempEvent );
      tempEvent = NULL;
    }
    if ( vrbflg ) fprintf( stdout, " done.\n" );

    /* close the xml file */
    LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );
    if ( vrbflg ) fprintf( stdout, "Closed the xml file.\n" );

  }

  goto cleanexit;

  cleanexit:

  /* free the frame cache */
  if( frInCache ) LAL_CALL( LALDestroyFrCache( &status, &frInCache ), &status );
  if ( frInType ) free( frInType );

  if ( vrbflg ) fprintf( stdout, "Checking memory leaks and exiting.\n" );
  LALCheckMemoryLeaks();

  exit(0);

}/* main function end */ 
/* -------------------------------------------------------------------------- */


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
"lalapps_nullstream [options]\n"\
"  --help                      display this message\n"\
"  --verbose                   print progress information\n"\
"  --version                   print version information and exit\n"\
"  --ifo-tag STRING            set STRING to the ifo-tag of the bank file\n"\
"  --user-tag STRING           set STRING to tag the file names\n"\
"\n"
#define USAGE2 \
"  --cohbank-file COHBANKFILE  read template bank parameters from COHBANKFILE\n"\
"  --sample-rate RATE          set data sample rate to RATE\n"\
"  --segment-length LEN        set LEN to same value used in inspiral.c\n"\
"  --dynamic-range-exponent N  set N to same value used in inspiral.c\n"\
"  [--h1-slide H1_SLIDE        Slide H1 data by multiples of H1_SLIDE]\n"\
"  [--h2-slide H2_SLIDE        Slide H2 data by multiples of H2_SLIDE]\n"\
"  --frame-type TAG            input data is contained in frames of type TAG\n"\
"\n"
#define USAGE3 \
"  --write-events              write events\n"\
"  --write-nullstat-series     write null statistic time series\n"\
"  --output-path               write files here\n"\
"  --h1-framefile              frame data for H1\n"\
"  --h2-framefile              frame data for H2\n"\
"  --h1-xml-file               H1 inspiral triggers\n"\
"  --h2-xml-file               H2 inspiral triggers\n"\
"\n"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
   struct option long_options[] =
   {
     {"verbose",                no_argument,       &vrbflg,            1 },
     {"help",                   no_argument,       0,                 'h'},
     {"version",                no_argument,       0,                 'v'},
     {"ifo-tag",                required_argument, 0,                 'I'},
     {"user-tag",               required_argument, 0,                 'B'},
     {"cohbank-file",           required_argument, 0,                 'u'},
     {"sample-rate",            required_argument, 0,                 'r'},
     {"segment-length",         required_argument, 0,                 'l'},
     {"dynamic-range-exponent", required_argument, 0,                 'e'},
     {"h1-slide",               required_argument, 0,                 'W'},
     {"h2-slide",               required_argument, 0,                 'X'},
     {"write-events",           no_argument,       &eventsOut,         1 },
     {"write-nullstat-series",  no_argument,       &nullStatOut,       1 },
     {"output-path",            required_argument, 0,                 'P'},
     {"h1-framefile",           required_argument, 0,                 'A'},
     {"h2-framefile",           required_argument, 0,                 'Z'},
     {"h1-xml-file",            required_argument, 0,                 'a'},
     {"h2-xml-file",            required_argument, 0,                 'b'},
     {"frame-type",             required_argument, 0,                 'S'},
     {0, 0, 0, 0}
   };

   int c;
   ProcessParamsTable *this_proc_param = procparams.processParamsTable;

   while (1)
   {
     /* getopt_long stores long options here */
     int option_index = 0;
     size_t optarg_len;

     c = getopt_long_only( argc, argv, "A:B:S:I:l:e:W:X:P:Z:h:r:u:v:a:b:",
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
           break;        
         }
         else
         {
           fprintf( stderr, "error parsing option %s with argument %s\n",
                     long_options[option_index].name, optarg );
           exit( 1 );
          }
          break;

       case 'A':
         strcpy(ifoframefile[1],optarg);
         H1file = 1;
         ADD_PROCESS_PARAM( "string", "%s", ifoframefile[1] );
         break;

       case 'Z':
         strcpy(ifoframefile[2],optarg);
         H2file = 1;
         ADD_PROCESS_PARAM( "string", "%s", ifoframefile[2] );
         break;

       case 'B':
         /* create storage for the user tag */
         optarg_len = strlen( optarg ) + 1;
         userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
         memcpy( userTag, optarg, optarg_len );
         ADD_PROCESS_PARAM( "string", "%s", optarg );
         break;

       case 'I':
         /* create storaged for the ifo-tag */
         optarg_len = strlen( optarg ) + 1;
         ifoTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
         memcpy( ifoTag, optarg, optarg_len );
         ADD_PROCESS_PARAM( "string", "%s", optarg );
         break;

       case 'P':
         memset( outputPath, 0, FILENAME_MAX * sizeof(CHAR) );
         snprintf( outputPath, FILENAME_MAX,"%s", optarg );
         ADD_PROCESS_PARAM( "string", "%s", outputPath );
         break;

       case 'S':
         optarg_len = strlen( optarg ) + 1;
         frInType = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
         memcpy( frInType, optarg, optarg_len );
         ADD_PROCESS_PARAM( "string", "%s", optarg );
         break;

       case 'h':
         fprintf( stdout, USAGE1 );
         fprintf( stdout, USAGE2 );
         fprintf( stdout, USAGE3 );
         exit( 0 );
         break;

       case 'l':
         numPointsSeg = atof(optarg);
         ADD_PROCESS_PARAM("float", "%e", numPointsSeg );
         break;

       case 'e':
         dynRangeExponent = atof(optarg);
         ADD_PROCESS_PARAM("float", "%e", dynRangeExponent );
         break;

       case 'r':
         sampleRate = atof(optarg);
         ADD_PROCESS_PARAM("float", "%e", sampleRate );
         break;

       case 'u':
         /* create storage for the cohbank filename */
         optarg_len = strlen( optarg ) + 1;
         cohbankFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
         memcpy( cohbankFileName, optarg, optarg_len );
         /*strcpy(cohbankFileName, optarg);*/
         ADD_PROCESS_PARAM( "string", "%s", cohbankFileName );
         char tempName[256];
         char *duration =NULL;
         strcpy(tempName, cohbankFileName);
         duration = strtok(tempName,"-");
         duration = strtok(NULL,"-");
         duration = strtok(NULL,"-");
         duration = strtok(NULL,".");
         bankDuration = atoi(duration);
         duration = NULL;
         break;

       case 'a':
         /*create storage for the H1 inspiral xml file name */
         optarg_len = strlen( optarg ) + 1;
         xmlFileNameH1 = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
         memcpy( xmlFileNameH1, optarg, optarg_len );
         ADD_PROCESS_PARAM( "string", "%s", xmlFileNameH1 );
         break;

       case 'b':
         /*create storage for the H2 inspiral xml file name */
         optarg_len = strlen( optarg ) + 1;
         xmlFileNameH2 = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
         memcpy( xmlFileNameH2, optarg, optarg_len );
         ADD_PROCESS_PARAM( "string", "%s", xmlFileNameH2 );
         break;
         

       case 'W':
         /* Read in time-slide step for H1 */
         slideStep[1] = atof(optarg);
         ADD_PROCESS_PARAM("float", "%e", slideStep[1]);
         break;

       /* Read in time-slide step for H2 */
       case 'X':
         slideStep[2] = atof(optarg);
         ADD_PROCESS_PARAM("float", "%e", slideStep[2]);
         break;

       case 'v':
         /* print version information and exit */
         fprintf( stdout, "Null stream code\n"
               "Messaritaki <emess@caltech.ed>\n"
               "CVS Version: " CVS_ID_STRING "\n"
               "CVS Tag: " CVS_NAME_STRING "\n" );
         fprintf( stdout, lalappsGitID );
         exit( 0 );
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


   /* check sample rate has been given */
   if ( sampleRate < 0 )
   {
     fprintf( stderr, "--sample-rate must be specified\n" );
     exit( 1 );
   }

   if ( numPointsSeg < 0 )
   {
     fprintf( stderr, 
       "--segment-length must be specified and set to the same value used when the C-data was generated.\n");
     exit( 1 );
   }

   if ( dynRangeExponent < 0 )
   {
     fprintf( stderr, 
       "--dynamic-range-exponent must be specified and set to the same value as was used when the C-data was generated.\n");
     exit( 1 );
   }

   if ( ! cohbankFileName )
   {
     fprintf( stderr, "--cohbank-file must be specified\n" );
     exit( 1 );
   }

   /* check that a channel has been requested and fill the ifo */
   if ( ! frInType )
   {
     fprintf( stderr, "--channel-name must be specified\n" );
     exit( 1 );
   }

   if( !ifoTag )
   {
     fprintf(stderr, "--ifo-tag must be specified for file naming\n" );
     exit( 1 );
   }

   return 0;
}

#undef ADD_PROCESS_PARAM

