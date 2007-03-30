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


/*
 *
 * variables that control program behaviour
 *
 */

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

/*
 *
 * more initialization parameters to be inserted here, 
 * once we know what we need.
 *
 */


int main( int argc, char *argv[] )
{
  FrChanIn      frChan;

  /* frame output data */
  struct FrFile *frOutFile  = NULL;
  struct FrameH *outFrame   = NULL;
  FrStream     *frStream = NULL;

  /* output */
  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;

  FILE *filePtr[4];

  CHAR   fileName[FILENAME_MAX];
  CHAR   framename[FILENAME_MAX];
  CHAR   xmlname[FILENAME_MAX];
  CHAR   cohdataStr[LALNameLength];
  CHAR   caseIDChars[6][LIGOMETA_IFOS_MAX] = {"0","0","0","0","0","0"};

  INT4   cohSegLength     = 4; /* should match hardcoded value in inspiral.c */
  INT4   numPoints        = 0;
  UINT4  numSegments      = 1;      /* number of segments */
  UINT4  numBeamPoints    = 3721;   /* number of sky position templates */
  UINT4  nCohDataFr       = 0;
  UINT8  eventID          = 0;

  REAL4  m1               = 0.0;
  REAL4  m2               = 0.0;
  REAL4  dynRange         = 0.0;

  /* variables for initializing tempTime to account for time-slides */
  UINT8  triggerNumber    = 0;
  UINT8  slideNumber      = 0;
  UINT8  slideSign        = 0;

  /* counters and other variables */
  INT4   j,k,l,w,kmax;
  REAL4  theta,phi,vPlus,vMinus;
  UINT4  numDetectors     = 0;
  REAL8  tempTime[6]      = {0.0,0.0,0.0,0.0,0.0,0.0};
  INT4   timeptDiff[5]    = {0,0,0,0,0};
  UINT2  caseID[6]        = {0,0,0,0,0,0}; /* H1 L V G T H2 */
  INT4   numTriggers      = 0;
  INT4   numCoincs        = 0;

  FrCache      *frInCache   = NULL;

  SnglInspiralTable    *currentTrigger = NULL;
  SnglInspiralTable    *cohbankEventList=NULL;

  CoincInspiralTable   *coincHead = NULL;
  CoincInspiralTable   *thisCoinc = NULL;

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummaryTable   *searchSummList = NULL;


  NullStatInitParams   *nullInspInitParams = NULL;
  NullStatParams       *nullInspFilterParams = NULL;
  NullStatInputInput   *nullInspFilterInput = NULL;
  NullStatCVector      *nullInspCVec = NULL;
  MultiInspiralTable   *thisEvent = NULL;
  MultiInspiralTable   *tempTable = NULL;
  MetadataTable        savedEvents;
  COMPLEX8TimeSeries   tempSnippet;

  char nameArrayBeam[6][256] = {"0","0","0","0","0","0"}; /*beam files*/
  char nameArrayCData[6][256] = {"0","0","0","0","0","0"};/*cData chan names*/

  set_debug_level( "1" ); /* change with parse option */

  /* set the dynamic range; needed for distNorm, templateNorm calculation */
  dynRange = pow( 2.0, dynRangeExponent );

  /* set other variables */
  numPoints = sampleRate * cohSegLength;
  savedEvents.multiInspiralTable = NULL;
  k = 0;

  /* read in the frame files */

  if( verbose ) fprintf(stdout, "reading in the frame files\n");
  for ( k=0; k<LAL_NUM_IFO ; k++)
  {
    if( ifoframefile[k] )
    {
      LAL_CALL( LALFrOpen(&status,&frStream,NULL,ifoframefile[k]), &status);   
      if(!frStream)
      {
        fprintf(stdout,"The file %s does not exist - exiting...\n", 
                ifoframefile[k]);
        goto cleanexit;
      }
    }
  }


  /* read in the cohbank trigger ligo lw xml file */
  numTriggers = XLALReadInspiralTriggerFile( &cohbankEventList, 
                   &currentTrigger, &searchSummList, &inputFiles, 
                   cohbankFileName );

  fprintf(stdout,"Reading templates from %s\n",cohbankFileName);

  if ( numTriggers < 0 )  /* no triggers found */
  {
    fprintf(stderr, "Error reading triggers from file %s", cohbankFileName);
    exit( 1 );
  }
  else if ( numTriggers == 0 )  /* no triggers found */
  {
    if( vrbflg )
    {
      fprintf( stdout,
               "%s contains no triggers - the coherent bank will be empty\n",
               cohbankFileName );
    }
  }
  else  /* triggers do exist */
  {
    if( vrbflg )
    {
      fprintf( stdout,
               "Read in %d triggers from the file %s\n", numTriggers,
               cohbankFileName );
    }

    /* pair up the coincidences */
    numCoincs = XLALRecreateCoincFromSngls( &coincHead, cohbankEventList );
    if( numCoincs < 0 )
    {
      fprintf(stderr, "Unable to reconstruct coincs from single ifo triggers");
      exit( 1 );
    }
    else if ( vrbflg )
    {
      fprintf( stdout,
               "Recreated %d coincs from the %d triggers\n", numCoincs,
               numTriggers );
    }

    /* loop over coincident triggers to compute the null statistic */
    for( thisCoinc=coincHead; thisCoinc; thisCoinc=thisCoinc->next)
    {
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
          memcpy( caseIDChars[k], &thisCoinc->snglInspiral[k]->ifo, 
                  sizeof(caseIDChars[k] - 1) );
          eventID = thisCoinc->snglInspiral[k]->event_id;
          if( verbose ) fprintf(stdout,"eventID = %Ld\n",eventID );

          /* Parse eventID to get the slide number */
          triggerNumber = eventID % 100000;
          slideNumber = ((eventID % 100000000) - triggerNumber)/100000;
          slideSign = (eventID % 1000000000) - slideNumber*100000 - 
                      triggerNumber;

          /* Store CData frame name now for reading its frame-file */
          /* later within thisCoinc-ident loop*/
          LALSnprintf( nameArrayCData[k], LALNameLength*sizeof(CHAR), 
                      "%s:%s_CData_%d", caseIDChars[k],frInType,eventID );
        }
      }/* Closes loop for k; finished noting the participating ifos
                and the eventID for this coincident trigger*/
      
      /* Initialize tempTime to account for time-slides (wrt H1, hence j<5)*/
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

      l=0;
      if( G1file ) l++;
      if( H1file ) l++;
      if( H2file ) l++;
      if( L1file ) l++;
      if( T1file ) l++;
      if( V1file ) l++;

      if( (INT4)numDetectors != l )
      {
        fprintf( stderr, "You have events for %d detectors, 
                but specified frame files for %d detectors\n",numDetectors,l);
        if( (INT4)numDetectors > l )
        {
          fprintf( stderr, "You must specify more frame files. Exiting...\n");
          exit(1);
        }
        else
        {
          if( verbose ) fprintf( stderr, "the number of detectors is smaller 
             than the number of frame files you specified. Exiting...\n");
          exit(1)
        }
      }

      l = 0;

      if( verbose ) fprintf(stdout,"numDetectors = %d\n", numDetectors);
      if( verbose ) fprintf(stdout,"caseID = %d %d %d %d %d %d 
         (G1,H1,H2,L1,V1,T1)\n", caseID[0], caseID[1], caseID[2], caseID[3], 
         caseID[4], caseID[5]);


      /* Initialize the necessary structures for thisCoinc-ident trigger*/

      if( !(nullStatInitParams = (NullStatisticInitParams *) 
          calloc(1,sizeof(NullStatisticInitParams)) ))
      {
        fprintf( stdout, 
         "could not allocate memory for nullInspiral init params\n" );
        goto cleanexit;
      }

      /* Initialize the null param structure for thisCoinc trigger */ 
      nullStatInitParams->numDetectors    = numDetectors;
      nullStatInitParams->numSegments     = numSegments;
      nullStatInitParams->numPoints       = numPoints;
      nullStatInitParams->numBeamPoints   = numBeamPoints;
      nullStatInitParams->nullStatOut     = nullStatOut;

      /* create the data structures needed */

      if ( verbose ) fprintf( stdout, "initializing...\n " );

      /* initialize null statistic functions */
      LAL_CALL( LALNullStatisticInputInit (&status, &nullStatInput,
                nullStatInitParams), &status );

      m1 = thisCoinc->snglInspiral[kmax]->mass1;
      m2 = thisCoinc->snglInspiral[kmax]->mass2;

      nullStatInput->tmplt = (InspiralTemplate *) 
         LALCalloc(1,sizeof(InspiralTemplate) );
      nullStatInput->tmplt->mass1 = m1;
      nullStatInput->tmplt->mass2 = m2;
      nullStatInput->tmplt->totalMass = m1 + m2;
      nullStatInput->tmplt->mu = m1 * m2 / (m1 + m2);
      nullStatInput->tmplt->eta = (m1 * m2) / ((m1 + m2) * (m1 + m2 ));

      if (verbose)  fprintf( stdout, "m1:%f m2:%f totalmass:%f mu:%f eta%f\n", 
         m1, m2,nullStatInput->tmplt->totalMass,nullStatInput->tmplt->mu,
         nullStatInput->tmplt->eta);

       LAL_CALL( LALNullStatisticParamsInit (&status, &nullStatParams,
                 nullStatInitParams),&status );

       nullStatParams->deltaT                  = 1/((REAL4) sampleRate);  
       nullStatParams->cohSNRThresh            = nullStatThresh;
       nullStatParams->cohSNROut               = nullStatOut;
       nullStatParams->numTmplts               = 1;
       nullStatParams->fLow                    = fLow;
       nullStatParams->maximizeOverChirp       = maximizeOverChirp;

       for( j=0; j<LAL_NUM_IFO; j++ )  /* what does this do??? */
       {
         nullStatParams->detIDVec->data[j] = caseID[j];
       }


       /* Read in the snippets associated with thisCoinc trigger */
       l = 0;
       for( j=0; j<LAL_NUM_IFO; j++ )
       {
         if( caseID[j] )
         {
           if( verbose ) fprintf(stdout, "getting the COMPLEX8TimeSeries\n");
           LAL_CALL(LALFrOpen(&status,&frStream,NULL,ifoframefile[j]), &status);
           if(!frStream)
           {
             fprintf(stdout,
               "The file %s does not exist - exiting...\n", ifoframefile[j] );
             goto cleanexit;
           }

           frChan.name = nameArrayCData[j];
           /*CHECK: Need to replace index l with an index that assigns these
            CData to the correct ifo*/
           LAL_CALL( LALFrGetCOMPLEX8TimeSeries( &status, 
              &(cohInspCVec->cData[l]), &frChan, frStream), &status);
           /*CHECK: Need to worry about WRAPPING of time-slides here */
           /* tempTime is the start time of cData plus -(time slide)*/
           tempTime[l] += cohInspCVec->cData[l].epoch.gpsSeconds + 
                          cohInspCVec->cData[l].epoch.gpsNanoSeconds * 1e-9;
           if( verbose ) fprintf(stdout,"tempTime = %f\n",tempTime[l]);
           LAL_CALL( LALFrClose( &status, &frStream ), &status );

           /* CHECK: delete this after updating the cohinspfilterparams defn
             cohInspFilterParams->segNorm[l] = thisCoinc->snglInspiral[k]
           */
           nullStatParams->sigmasq[l] = thisCoinc->snglInspiral[j]->sigmasq;
           l++;
         }/* closes if( caseID[j] ) */
       } /* closes for( j=0; j<LAL_NUM_IFO; j++ ) */


       /*
        * lines that do the calculation of templateNorm have been skipped
        * since that is probably not necessary due to the fact that the
        * sigmasq is read.
        *
        * also skipping the commensuration of the c-data snippets. This
        * is necessary in principle, but the H1-H2 snippets are always
        * commensurate so we postpone writing that part of the code.
        */

        /* calculation of the null statistic */
        









