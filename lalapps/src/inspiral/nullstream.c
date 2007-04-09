/*-----------------------------------------------------------------------
 *
 * File Name: nullstream.c
 *
 * Author: Messaritaki, E.
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
#include <lal/StochasticCrossCorrelation.h>
#include <lal/DetectorSite.h>
#include <lal/Random.h>
#include <lal/LALInspiral.h>
#include <lal/CoherentInspiral.h>
#include <lal/NullStatistic.h>
#include <lal/LALStatusMacros.h>
#include <lal/SkyCoordinates.h>

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "nullstream"
#define CVS_NAME_STRING "$Name$"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

#define rint(x) (floor((x)+0.5))


/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* input data parameters */

INT4   sampleRate           = -1;  /* sample rate of filter data   */
INT4   numPointsSeg         = -1;  /* set to segment-length from inspiral.c */
REAL4  fLow                 = -1;  /* low frequency cutoff         */
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
REAL8  slideStep[LAL_NUM_IFO]     = {0.0,0.0,0.0,0.0,0.0,0.0};
int    bankDuration     = 0;

CHAR  *cohbankFileName     = NULL;   /* name of input template bank  */
UINT4  nullStatOut         = 0;      /* default is not to write frame */
UINT4  eventsOut           = 0;      /* default is not to write events */
REAL4  nullStatThresh      = -1;
INT4   maximizeOverChirp   = 0;      /* default is no clustering */
INT4   verbose             = 0;
CHAR   outputPath[FILENAME_MAX];
CHAR  *frInType            = NULL;   /* type of data frames */

INT8          gpsStartTimeNS     = 0;   /* input data GPS start time ns */
LIGOTimeGPS   gpsStartTime;             /* input data GPS start time    */
INT8          gpsEndTimeNS       = 0;   /* input data GPS end time ns   */
LIGOTimeGPS   gpsEndTime;               /* input data GPS end time      */
int           gpsStartTimeTemp   = 0;   /* input data GPS start time ns */
int           gpsEndTimeTemp     = 0;   /* input data GPS start time ns */

LALStatus             status;
LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

CHAR  *userTag          = NULL;         /* string the user can tag with */
CHAR  *ifoTag           = NULL;         /* string to tag IFOs    */


int main( int argc, char *argv[] )
{
  FrChanIn      frChan;

  /* frame output data */
  struct FrFile *frOutFile  = NULL;
  struct FrameH *outFrame   = NULL;
  FrStream      *frStream   = NULL;

  /* output */
  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;

  FILE *filePtr[4];

  CHAR  fileName[FILENAME_MAX];
  CHAR  framename[FILENAME_MAX];
  CHAR  xmlname[FILENAME_MAX];
  CHAR  cohdataStr[LALNameLength];

  INT4   segLength      = 4;  /* should match hardcoded value in inspiral.c */
  INT4   numPoints      = 0;
  UINT4  numSegments    = 1;  /* number of segments */
  UINT4  nCohDataFr     = 0;  /* what is this? */
  UINT8  eventID        = 0;

  REAL4  m1             = 0.0;
  REAL4  m2             = 0.0;
  REAL4  dynRange       = 0.0;

  /* variables for initializing tempTime to account for time-slides */
  UINT8  triggerNumber  = 0;
  UINT8  slideNumber    = 0;
  UINT8  slideSign      = 0;

  /* counters and other variables */
  INT4   j, k, l, w, kidx;
  UINT4  numDetectors            = 0;
  REAL8  tempTime[LAL_NUM_IFO]   = {0.0,0.0,0.0,0.0,0.0,0.0}; 
  INT4   timeptDiff[5]           = {0,0,0,0,0};
  INT4   numTriggers             = 0;
  INT4   numCoincs               = 0;

  FrCache              *frInCache        = NULL;

  SnglInspiralTable    *currentTrigger   = NULL;
  SnglInspiralTable    *cohbankEventList = NULL;

  CoincInspiralTable   *coincHead = NULL;
  CoincInspiralTable   *thisCoinc = NULL;

  SearchSummvarsTable  *inputFiles     = NULL;
  SearchSummaryTable   *searchSummList = NULL;


  NullStatInitParams      *nullStatInitParams   = NULL;
  NullStatParams          *nullStatParams       = NULL;
  NullStatInputParams     *nullStatInputParams  = NULL;
  CVector                 *CVec                 = NULL;
  MultiInspiralTable      *thisEvent            = NULL;
  MultiInspiralTable      *tempTable            = NULL;
  MetadataTable           savedEvents;
  COMPLEX8TimeSeries      tempSnippet;

  /* cData channel names */
  char nameArrayCData[LAL_NUM_IFO][256] = {"0","0","0","0","0","0"}; 

  set_debug_level( "1" ); /* change with parse option */

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable,
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *)
    calloc( 1, sizeof(ProcessParamsTable) );

  arg_parse_check( argc, argv, procparams );
  if (verbose)  fprintf(stdout, "Called parse options.\n");

  /* wind to the end of the process params table */
  for ( this_proc_param = procparams.processParamsTable; this_proc_param->next;
        this_proc_param = this_proc_param->next );


  /* set the dynamic range */
  dynRange = pow( 2.0, dynRangeExponent );

  /* set other variables */
  numPoints = sampleRate * segLength;
  savedEvents.multiInspiralTable = NULL;
  k = 0;

  /* read in the frame files */
  if ( verbose ) fprintf(stdout, "Reading in the frame files.\n");
  for ( k=0; k<LAL_NUM_IFO ; k++)
  {
    if ( ifoframefile[k] )
    {
      LAL_CALL( LALFrOpen(&status, &frStream, NULL, ifoframefile[k]), &status); 
      if (!frStream)
      {
        fprintf(stdout,"The file %s does not exist - exiting.\n", 
                ifoframefile[k]);
        goto cleanexit;
      }
    }
  }


  /* read in the cohbank trigger ligo lw xml file */
  numTriggers = XLALReadInspiralTriggerFile( &cohbankEventList, 
                  &currentTrigger, &searchSummList, &inputFiles, 
                  cohbankFileName );

  fprintf(stdout,"Reading templates from %s.\n",cohbankFileName);

  if ( numTriggers < 0 )  /* no triggers found */
  {
    fprintf(stderr, "Error reading triggers from file %s.\n", cohbankFileName);
    exit( 1 );
  }
  else if ( numTriggers == 0 )  /* no triggers found */
  {
    if ( vrbflg )
    {
      fprintf( stdout,
               "%s contains no triggers - the coherent bank will be empty.\n",
               cohbankFileName );
    }
  }
  else  /* triggers do exist */
  {
    if ( vrbflg )
    {
      fprintf( stdout,
               "Read in %d triggers from the file %s.\n", numTriggers,
               cohbankFileName );
    }

    /* pair up the coincidences */
    numCoincs = XLALRecreateCoincFromSngls( &coincHead, cohbankEventList );
    if ( numCoincs < 0 )
    {
      fprintf(stderr, "Unable to reconstruct coincs from single ifo triggers.");
      exit( 1 );
    }
    else if ( vrbflg )
    {
      fprintf( stdout,
               "Recreated %d coincs from the %d triggers.\n", numCoincs,
               numTriggers );
    }

    /* loop over coincident triggers to compute the null statistic */
    for ( thisCoinc=coincHead; thisCoinc; thisCoinc=thisCoinc->next)
    {
      numDetectors = thisCoinc->numIfos;

      /* l is another detector index */
      l=0;

      /* Note the participating ifos and the eventID for this coincidenice */
      for ( k=0 ; k<LAL_NUM_IFO ; k++)
      {
        if ( thisCoinc->snglInspiral[k] )
        {
          eventID = thisCoinc->snglInspiral[k]->event_id;
          if ( verbose ) fprintf(stdout,"eventID = %Ld.\n",eventID );

          /* Parse eventID to get the slide number */
          triggerNumber = eventID % 100000;
          slideNumber = ((eventID % 100000000) - triggerNumber)/100000;
          slideSign = (eventID % 1000000000)-(slideNumber*100000)-triggerNumber;

          /* Store CData frame name  */
          LALSnprintf( nameArrayCData[k], LALNameLength*sizeof(CHAR), 
            "%s:%s_CData_%d", &thisCoinc->snglInspiral[k]->ifo, 
            frInType, eventID );
          kidx = k;
        }
      } 

      /* Initialize tempTime to account for time-slides */
      for( j=0; j<(LAL_NUM_IFO-1); j++)
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

      l=0;
      if ( G1file ) l++;
      if ( H1file ) l++;
      if ( H2file ) l++;
      if ( L1file ) l++;
      if ( T1file ) l++;
      if ( V1file ) l++;


      if ( (INT4)numDetectors != l )
      {
        fprintf( stderr, "You have events for %d detectors, but specified 
                          frame files for %d detectors.\n",numDetectors,l);
        if ( (INT4)numDetectors > l )
        {
          fprintf( stderr, "Too few frame files specified. Exiting.\n");
          exit(1);
        }
        else
        {
          fprintf( stderr, "Too many frame files specified. Exiting.\n");
          exit(1)
        }
      }

      l = 0;

      if ( verbose ) fprintf(stdout,"numDetectors = %d\n", numDetectors);


      /* Initialize the necessary structures for thisCoinc-ident trigger*/

      if ( !(nullStatInitParams = (NullStatInitParams *) 
          calloc(1,sizeof(NullStatInitParams)) ))
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

      if ( verbose ) fprintf( stdout, "Initializing.\n " );

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

      if (verbose) fprintf(stdout,"m1:%f, m2:%f, Mtotal:%f, mu:%f, eta:%f.\n", 
         m1, m2, nullStatInputParams->tmplt->totalMass,
         nullStatInputParams->tmplt->mu, nullStatInputParams->tmplt->eta);

       XLALNullStatisticParamsInit(&nullStatParams, nullStatInitParams);

       nullStatParams->numTmplts         = 1;
       nullStatParams->maxOverChirp      = maximizeOverChirp;
       nullStatParams->fLow              = fLow;
       nullStatParams->deltaT            = 1/((REAL4) sampleRate);  
       nullStatParams->nullStatThresh    = nullStatThresh;
       nullStatParams->nullStatOut       = nullStatOut;

       /*
        * we only need H1 and H2 at the moment, but assign
        * all, for future use
        */
       nullStatParams->detVector->detector[LAL_IFO_G1] = 
         lalCachedDetectors[LALDetectorIndexGEO600DIFF];
       nullStatParams->detVector->detector[LAL_IFO_H1] =
         lalCachedDetectors[LALDetectorIndexLHODIFF];
       nullStatParams->detVector->detector[LAL_IFO_H2] =
         lalCachedDetectors[LALDetectorIndexLHODIFF];
       nullStatParams->detVector->detector[LAL_IFO_L1] =
         lalCachedDetectors[LALDetectorIndexLLODIFF];
       nullStatParams->detVector->detector[LAL_IFO_T1] =
         lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
       nullStatParams->detVector->detector[LAL_IFO_V1] =
         lalCachedDetectors[LALDetectorIndexVIRGODIFF];


       /* Read in the snippets associated with thisCoinc trigger */
       for ( j=0; j<LAL_NUM_IFO; j++ )
       {
         if ( (j == LAL_IFO_H1) || (j == LAL_IFO_H2) ) 
         { 
           if ( verbose ) fprintf(stdout, "Getting the COMPLEX8TimeSeries.\n");
           LAL_CALL( LALFrOpen( &status, &frStream, NULL, ifoframefile[j]), 
              &status);
           if (!frStream)
           {  
             fprintf(stdout,
               "The file %s does not exist - exiting.\n", ifoframefile[j] );
             goto cleanexit;
           }

           frChan.name = nameArrayCData[j];

           LAL_CALL( LALFrGetCOMPLEX8TimeSeries( &status, 
              &(CVec->cData[j]), &frChan, frStream), &status);

           /* Need to worry about WRAPPING of time-slides             */
           /* tempTime is the start time of cData plus - (time slide) */
           tempTime[j] += CVec->cData[j].epoch.gpsSeconds + 
                          CVec->cData[j].epoch.gpsNanoSeconds * 1e-9;
           if ( verbose ) fprintf( stdout, "tempTime = %f\n", tempTime[j]);
  
           LAL_CALL( LALFrClose( &status, &frStream ), &status );
  
           /* set sigma-squared */
           nullStatParams->sigmasq[j] = thisCoinc->snglInspiral[j]->sigmasq;
           j++;
         }
         else
         {
           if (verbose) fprintf(stdout,"No data to read for G1,L1,T1 or V1.\n");
         }
       }      /* closes for( j=0; j<LAL_NUM_IFO; j++ ) */


       /*
        * skipping the commensuration of the c-data snippets. This
        * is necessary in principle, but the H1-H2 snippets are always
        * commensurate so we postpone writing that part of the code.
        */

        /* calculation of the null statistic */
        
