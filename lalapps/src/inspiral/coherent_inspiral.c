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
#include <lal/LALStatusMacros.h>
#include <lal/SkyCoordinates.h>

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "coherent_inspiral"
#define CVS_NAME_STRING "$Name$"



static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

static void
ClearStatus (LALStatus *status);

int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

#define rint(x) (floor((x)+0.5))


/*
 *
 * variables that control program behaviour
 *
 */


/* input data parameters */

INT4   sampleRate           = -1;       /* sample rate of filter data   */
INT4   numPointsSeg         = -1;/* set to segment-length used in inspiral.c */
REAL4  fLow                 = -1;       /* low frequency cutoff         */
REAL4  dynRangeExponent     = -1;/* set to same value used in inspiral.c */

/*Coherent code specific inputs*/

char   H1filename[256];        
char   Lfilename[256];         
char   GEOfilename[256];       
char   VIRGOfilename[256];      
char   TAMAfilename[256];       
char   H2filename[256];        

INT4 H1file = 0;
INT4 H2file = 0;
INT4 Lfile = 0;
INT4 GEOfile = 0;
INT4 TAMAfile = 0;
INT4 VIRGOfile = 0;
                                
/* input time-slide parameters */
REAL8  slideStep[6]     = {0.0,0.0,0.0,0.0,0.0,0.0};
int    bankDuration     = 0;
CHAR   bankFileName[FILENAME_MAX]; /* name of input template bank  */
UINT4  cohSNROut            = 0;    /* default is not to write frame */
UINT4  eventsOut            = 0;    /* default is not to write events */
REAL4  cohSNRThresh         = -1;
INT4   maximizeOverChirp    = 0;    /* default is no clustering */
INT4   verbose              = 0;
CHAR   outputPath[FILENAME_MAX];
CHAR  *frInType         = NULL;         /* type of data frames          */

INT8  gpsStartTimeNS   = 0;         /* input data GPS start time ns */
LIGOTimeGPS gpsStartTime;           /* input data GPS start time    */
INT8  gpsEndTimeNS     = 0;         /* input data GPS end time ns   */
LIGOTimeGPS gpsEndTime;             /* input data GPS end time      */

int  gpsStartTimeTemp   = 0;         /* input data GPS start time ns */
int  gpsEndTimeTemp   = 0;         /* input data GPS start time ns */

LALStatus             status;
LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

CHAR  *userTag          = NULL;         /* string the user can tag with */
CHAR  *ifoTag           = NULL;         /* string to tag parent IFOs    */
INT4  globFrameData     = 0;            /* glob to get frame data */

/* Params to convert from geocentric to equatorial */

ConvertSkyParams             convertParams;
SkyPosition		     tempSky;
MultiInspiralTable           *thisEventTemp = NULL;

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
  CHAR   caseIDChars[4][LIGOMETA_IFOS_MAX] = {"0","0","0","0"};
  CHAR   tempBankName[FILENAME_MAX]; /* name of input template bank  */
  CHAR   nextURL[256];
  char   networkName[256] = "INSPIRAL_";
  
  int    frameIndex         = 0;
  int    incrementH1        = 1;
  int    incrementH2        = 1;
  int    incrementL1        = 1;
  int    incrementV1        = 1;
  int    incrementG1        = 1;
  int    incrementT1        = 1;
  int    firstPass        = 0;
  
  INT4   numTmplts        = 0; /* number of templates */
  INT4   cohSegLength     = 4; /* This should match hardcoded value in inspiral.c */
  INT4   numPoints        = 0;
  INT4   startTemplate    = -1;           
  INT4   stopTemplate     = -1;
  INT4   numChannels      = 0;
  INT4   frEnd            = 0;
  UINT4  numSegments      = 1;       /* number of segments */
  UINT4  numBeamPoints    = 3721;       /* number of sky position templates */
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
  INT4   i,j,k,l,w,p;
  REAL4  theta,phi,vPlus,vMinus;
  UINT4  numDetectors     = 0;
  UINT4  totalNumFrames   = 0;
  UINT4  framesPerDetector= 0;
  REAL8  tempTime[6]      = {0.0,0.0,0.0,0.0,0.0,0.0};
  INT4   timeptDiff[5]    = {0,0,0,0,0};
  UINT2  caseID[6]        = {0,0,0,0,0,0}; /* H1 L V G T H2 */
  INT4   h1ChanNum        = 0;
  INT4   h2ChanNum        = 0;
  INT4   lChanNum         = 0;
  INT4   geoChanNum       = 0;
  INT4   virgoChanNum     = 0;
  INT4   tamaChanNum      = 0;
  INT4   h1ChanNumTemp    = 0;
  INT4   h2ChanNumTemp    = 0;
  INT4   lChanNumTemp     = 0;
  INT4   geoChanNumTemp   = 0;
  INT4   virgoChanNumTemp = 0;
  INT4   tamaChanNumTemp  = 0;
  INT4   chanNumArray[6]  = {-1, -1, -1, -1, -1, -1};

  FrCache      *frGlobCache = NULL;
  FrCache      *frInCache   = NULL;
  FrCache      *tempCache   = NULL;
  FrCache      *tempCache2  = NULL;
  FrStat       *tempFrame   = NULL;
  FrCacheSieve  sieve;
  FrCacheSieve  tempSieve;

  CoherentInspiralInitParams   *cohInspInitParams = NULL;
  CoherentInspiralFilterParams *cohInspFilterParams = NULL;
  CoherentInspiralFilterInput  *cohInspFilterInput = NULL;
  CoherentInspiralBeamVector   *cohInspBeamVec = NULL;
  CoherentInspiralCVector      *cohInspCVec = NULL;
  MultiInspiralTable           *thisEvent = NULL;
  MultiInspiralTable           *tempTable = NULL;
  MetadataTable                savedEvents;
  InspiralTemplate             *bankHead = NULL;
  InspiralTemplate             *bankTemp = NULL;
  InspiralTemplate             *bankTemp2 = NULL;
  COMPLEX8TimeSeries            tempSnippet;
  REAL4FrequencySeries         *segNormVector = NULL;
  
  char namearray[6][256]  = {"0","0","0","0","0","0"}; /* input frame files */
  char namearray2[6][256] = {"0","0","0","0","0","0"}; /* beam files */
  char namearray3[6][256] = {"0","0","0","0","0","0"}; /* cData chan names */
  char namearray4[6][256] = {"0","0","0","0","0","0"}; /* segNorm chan names */


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
  if (verbose)  fprintf(stdout, "called parse options..\n");

  /* wind to the end of the process params table */
  for ( this_proc_param = procparams.processParamsTable; this_proc_param->next;
      this_proc_param = this_proc_param->next );

  fprintf(stdout,"Reading templates from %s\n",bankFileName);

  /* read in the template bank from a ligo lw xml file */
    numTmplts = InspiralTmpltBankFromLIGOLw( &bankHead, bankFileName,
      startTemplate, stopTemplate );
  if ( numTmplts < 0 )
    {
    fprintf( stderr, "error: unable to read templates from %s\n", 
        bankFileName );
    goto cleanexit;
    }
  else if ( numTmplts == 0 )
    {
    /* if there are no tmplts, exit */
    fprintf( stderr, "no templates found in template bank file: %s\n"
        "exiting without searching for events.\n", bankFileName ); 
    goto cleanexit;
    }
  
  if ( verbose ) 
      fprintf( stdout, "parsed %d templates from %s\n", numTmplts, bankFileName );

  /* Now glob for frame data based on the gps start and duration in the */
  /* thinca input file name */

  if( globFrameData )
      {
      
      if ( verbose ) 
          fprintf( stdout, "globbing for *.gwf frame files in current directory\n");
      LAL_CALL( LALFrCacheGenerate( &status, &frGlobCache, NULL, NULL ), 
		&status );
      /* check we globbed at least one frame file */
      if ( ! frGlobCache->numFrameFiles )
         {
         fprintf( stderr, "error: no frame file files found\n");
         exit( 1 );
         }
      totalNumFrames=frGlobCache->numFrameFiles;
      /* sieve out the requested data type */
      memset( &sieve, 0, sizeof(FrCacheSieve) );
      sieve.srcRegEx = NULL;
      sieve.dscRegEx = NULL;
      sieve.urlRegEx = NULL;
      sieve.earliestTime = gpsStartTime.gpsSeconds;
      sieve.latestTime = gpsEndTime.gpsSeconds;
      LAL_CALL( LALFrCacheSieve( &status, &frInCache, frGlobCache, &sieve ), 
		&status );
      if( verbose ) fprintf(stdout,"num files after sieve: %d\n",frInCache->numFrameFiles);
      LAL_CALL( LALDestroyFrCache( &status, &frGlobCache ), &status );
      /* open the input data frame stream from the frame cache */
      /*LAL_CALL( LALFrCacheOpen( &status, &frStream, frInCache ), &status );*/

      /* set the mode of the frame stream to fail on gaps or time errors */
      /*frStream->mode = LAL_FR_DEFAULT_MODE;*/
      }

  /* Set the dynamic range */
  dynRange = pow( 2.0, dynRangeExponent );
  
  /* Loop over templates (or event id's) */
  
  numPoints = sampleRate * cohSegLength;
  bankTemp = bankHead;
  savedEvents.multiInspiralTable = NULL;
  k = 0; 

  if( !globFrameData )
      {
      h1ChanNum        = -1;
      h2ChanNum        = -1;
      lChanNum         = -1;
      geoChanNum       = -1;
      virgoChanNum     = -1;
      tamaChanNum      = -1;
      }
  
  framesPerDetector=totalNumFrames;
  
 /* while(framesPerDetector>=1)
  {*/
  
  for( i=0; i<numTmplts; i++)
    {
      memcpy( caseIDChars[k], &bankTemp->ifo, sizeof(caseIDChars[k] - 1) );
      if( verbose ) fprintf(stdout,"caseIDChars = %s %s %s %s\n",caseIDChars[0],caseIDChars[1],caseIDChars[2],caseIDChars[3] );
      eventID = bankTemp->event_id->id;
      if( verbose ) fprintf(stdout,"eventID = %Ld\n",eventID );

      /* Parse eventID to get the slide number */
      triggerNumber = eventID % 100000;
      slideNumber = ((eventID % 100000000) - triggerNumber)/100000;
      slideSign = (eventID % 1000000000) - slideNumber*100000 - triggerNumber;

      /* Initialize tempTime to account for time-slides */
      for( l=0; l<5; l++)
        {
        /* slideSign=0 is the same as a positive time slide */
        if(slideSign != 0)
            {
	        tempTime[l] = slideStep[l]*slideNumber*slideSign;
            }
        else
            {
            tempTime[l] -= slideStep[l]*slideNumber*slideSign;
            }
        }
      if( i != numTmplts - 1 )
          {
          bankTemp2 = bankTemp;
          bankTemp = bankTemp->next;
          }
      if( i != 0 && ((bankTemp->event_id->id != eventID) || (i == numTmplts - 1)))
          {
          if( k > 0 )
            {
            /* Here we should combine data for the last event */
            if( verbose ) fprintf(stdout,"combining data for the last event\n");
            for(j=0;j<4;j++)
              {
              if( !strcmp( caseIDChars[j],"H1" ) )
                {
                caseID[0] = 1;
                if( !globFrameData ) h1ChanNum++;
                }
              else if( !strcmp( caseIDChars[j], "L1" ) )
    		    {
    		    caseID[1] = 1;
		        if( !globFrameData ) lChanNum++;
		        }
		      else if( !strcmp( caseIDChars[j], "V1" ) )
		        {
		        caseID[2] = 1;
		        if( !globFrameData ) virgoChanNum++;
		        }
		      else if( !strcmp( caseIDChars[j], "G1" ) )
		        {
	            caseID[3] = 1;
		        if( !globFrameData ) geoChanNum++;
		        }
		      else if( !strcmp( caseIDChars[j], "T1" ) )
		        {
		        caseID[4] = 1;
		        if( !globFrameData ) tamaChanNum++;
		        }
		      else if( !strcmp( caseIDChars[j], "H2" ) )
		        {
		        caseID[5] = 1;
		        if( !globFrameData ) h2ChanNum++;
		        }
		      }


	        /* Now get the number of detectors */
	        l = 0;
	        for(j=0;j<6;j++)
	          {
	          if(caseID[j])
	            {
	            l++;
		        }
		      }
	        numDetectors = l;

	        /* Now check that the number of detectors matches the number of frame files provided */

	        l=0;
	        if( H1file ) l++;
	        if( H2file ) l++;
	        if( Lfile ) l++;
	        if( GEOfile ) l++;
	        if( TAMAfile ) l++;
	        if( VIRGOfile ) l++;
	      
	        if( ! globFrameData )
		      {
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
			      if( verbose ) fprintf( stdout, "One or more of the frame files specified will not be used for this event since the number of detectors is less than the number of frame files you specified.\n");
			      }
		        }
		      }

	        l = 0;

	        if( verbose ) fprintf(stdout,"numDetectors = %d\n", numDetectors);
	        if( verbose ) fprintf(stdout,"caseID = %d %d %d %d %d %d (H1,L1,V1,G1,T1,H2)\n", caseID[0], caseID[1], caseID[2], caseID[3], caseID[4], caseID[5]);


	        /* Initialize the necessary structures */

	        if( !(cohInspInitParams = (CoherentInspiralInitParams *) calloc(1,sizeof(CoherentInspiralInitParams)) ))
		      {
		      fprintf( stdout, "could not allocate memory for coherentInspiral init params\n" );
		      goto cleanexit;
		      }

	        cohInspInitParams->numDetectors            = numDetectors;
	        cohInspInitParams->numSegments             = numSegments;
	        cohInspInitParams->numPoints               = numPoints;
	        cohInspInitParams->numBeamPoints           = numBeamPoints;
	        cohInspInitParams->cohSNROut               = cohSNROut;

             
	        /* create the data structures needed for coherentInspiral */

	        if ( verbose ) fprintf( stdout, "initializing coherentInspiral...\n " );

	        /* initialize coherentInspiral filter functions */
            LAL_CALL( LALCoherentInspiralFilterInputInit (&status, &cohInspFilterInput,
						cohInspInitParams), &status );

            m1 = bankTemp2->mass1;
	        m2 = bankTemp2->mass2;
 
	        cohInspFilterInput->tmplt = (InspiralTemplate *)
		   LALCalloc(1,sizeof(InspiralTemplate) );
	        cohInspFilterInput->tmplt->mass1 = m1;
	        cohInspFilterInput->tmplt->mass2 = m2;
	        cohInspFilterInput->tmplt->totalMass = m1 + m2;
	        cohInspFilterInput->tmplt->mu = m1 * m2 / (m1 + m2);
	        cohInspFilterInput->tmplt->eta = (m1 * m2) / ((m1 + m2) * (m1 + m2 ));

	        if (verbose)  fprintf( stdout, "m1:%f m2:%f totalmass:%f mu:%f eta%f\n", m1, m2,cohInspFilterInput->tmplt->totalMass,cohInspFilterInput->tmplt->mu,cohInspFilterInput->tmplt->eta);

	        LAL_CALL( LALCoherentInspiralFilterParamsInit (&status, &cohInspFilterParams,
							     cohInspInitParams),&status );

	        cohInspFilterParams->deltaT                  = 1/((REAL4) sampleRate);
	        cohInspFilterParams->cohSNRThresh            = cohSNRThresh;  
	        cohInspFilterParams->cohSNROut               = cohSNROut;
	        cohInspFilterParams->numTmplts               = 1;
	        cohInspFilterParams->fLow                    = fLow;
	        cohInspFilterParams->maximizeOverChirp       = maximizeOverChirp;
  
	        for( j=0; j<6; j++ ) 
		      {
		      cohInspFilterParams->detIDVec->data[j] = caseID[j];
	          }

	        w=0;
	        for ( j=0; j<6; j++ ) 
		      {
		      if ( caseID[j] ) 
		        { 
		        cohInspFilterParams->detectorVec->detector[w++] = lalCachedDetectors[j];
		        }
		      }

	        if (caseID[5]) 
		      {
		      cohInspFilterParams->detectorVec->detector[numDetectors-1] = lalCachedDetectors[0];
		      }


	      /* Get the data that corresponds to this event from the glob cache */
              /* or from the specified frame files. */
              /* First, the file names and channel names must be set correctly */
            /*NEW CHANGES ALL ARE HERE*/
            if( globFrameData ) framesPerDetector=floor(totalNumFrames/numDetectors);
            else framesPerDetector=1; fprintf(stdout,"frames per detector is %i total number of frames are %i\n",framesPerDetector,totalNumFrames );
           
            
            if( globFrameData )
              {
              if(firstPass==0)
                {
                if(caseID[3]==1) strcat(networkName, "G1");
                if(caseID[0]==1) strcat(networkName, "H1");
                if(caseID[5]==1) strcat(networkName, "H2");
                if(caseID[1]==1) strcat(networkName, "L1");
                if(caseID[4]==1) strcat(networkName, "T1");
                if(caseID[2]==1) strcat(networkName, "V1");
                strcat(networkName, "C");
                firstPass=1;
                }
              
              /* fprintf(stdout, "the network is %s\n", networkName);*/
              if( verbose ) fprintf(stdout,"caseID = %d %d %d %d %d %d (H1,L1,V1,G1,T1,H2)\n", caseID[0], caseID[1], caseID[2], caseID[3], caseID[4], caseID[5]);
              if( caseID[0] )
                {
                memset( &tempSieve, 0, sizeof(FrCacheSieve) );
                /*memset( &tempCache, 0, sizeof(FrCache) );*/
                tempSieve.srcRegEx = "H1";
                tempSieve.dscRegEx = networkName;
                tempSieve.urlRegEx = NULL;
                tempSieve.earliestTime = bankTemp2->end_time.gpsSeconds+incrementH1;
                tempSieve.latestTime = bankTemp2->end_time.gpsSeconds + incrementH1;
                LAL_CALL( LALFrCacheSieve( &status, &tempCache, frInCache, &tempSieve ), &status );
                if( !tempCache->numFrameFiles )
                    {
                    
                    fprintf(stderr,"CacheSieve should only have returned 1 frame file for H1. Exiting..\n");
                    goto cleanexit;
                    }
                fprintf(stdout,"The first URL for H1 is %s and the number of frames in the cache is %i\n",tempCache->frameFiles->url,tempCache->numFrameFiles);
                if(tempCache->numFrameFiles>1)
                    {
                    tempFrame=tempCache->frameFiles+1;
                    strcpy(nextURL,tempFrame->url);
                    fprintf(stdout,"The second URL for H1 is %s\n",nextURL);
                    }
                else framesPerDetector=1;
                    
                    
    		    /* Need to strip the file name out of the url */
                char *tempName = NULL;
                char tempName2[256];
                char tempName3[256];
                char *nextStartTime =NULL;
                tempName = strtok(tempCache->frameFiles->url,"//");
                tempName = strtok(NULL,"/");
                while( tempName != NULL)
                    {
                    strcpy(tempName2,tempName);
                    tempName = strtok(NULL,"/");
                    }
                tempName = NULL;
                strcpy(tempName3, tempName2);
                nextStartTime=strtok(tempName3, "-");
                nextStartTime=strtok(NULL, "-");
                nextStartTime=strtok(NULL, "-");
                if(i>=numTmplts-1&&bankDuration>2048)
                  {
                  incrementH1=2049-(bankTemp2->end_time.gpsSeconds-atoi(nextStartTime));
                  }
                nextStartTime=NULL;
                fprintf(stdout, "The incrementH1 is %i\n", incrementH1);
                if( !strcmp( tempName2, namearray[0] ) )
                    {
                    /* if this file has been used, increment the chan # */
                    h1ChanNum++;
                        
                    }
                else
                    {
                    h1ChanNum = 0;/*start at zero unless debugging*/
                    strcpy( namearray[0], tempName2 );
                    }
                strcpy(namearray2[0],"HanfordBeamCoeff.dat");
                LALSnprintf( namearray3[0], LALNameLength*sizeof(CHAR), "H1:%s_CData_%d",frInType, h1ChanNum );
                LALSnprintf( namearray4[0], LALNameLength*sizeof(CHAR), "_SegNorm_%d", h1ChanNum );
                LAL_CALL( LALDestroyFrCache(&status, &tempCache), &status );
                tempCache = NULL;
                tempFrame = NULL;
                tempCache2 = NULL;
                }
              

		    if( caseID[1] )
                {
                memset( &tempSieve, 0, sizeof(FrCacheSieve) );
                tempSieve.srcRegEx = "L1";
                tempSieve.dscRegEx = networkName;;
                tempSieve.urlRegEx = NULL;
                tempSieve.earliestTime = bankTemp2->end_time.gpsSeconds+incrementL1;
                tempSieve.latestTime = bankTemp2->end_time.gpsSeconds + incrementL1;
                LAL_CALL( LALFrCacheSieve( &status, &tempCache, 
				       frInCache, &tempSieve ), &status );
    		    if( !tempCache->numFrameFiles )
                    {
                    fprintf(stderr,"The total number of templates returned was %i\n",tempCache->numFrameFiles);
                    fprintf(stderr,"CacheSieve should only have returned 1 frame file for L1. Exiting..\n");
                    goto cleanexit;
                    }
                 if(tempCache->numFrameFiles==1) framesPerDetector=1;
                 /* Need to strip the file name out of the url */
                char *tempName = NULL;
                char tempName2[256];
                char tempName3[256];
                char *nextStartTime =NULL;
                tempName = strtok(tempCache->frameFiles->url,"//");
                tempName = strtok(NULL,"/");
                while( tempName != NULL)
                    {
                    strcpy(tempName2,tempName);
                    tempName = strtok(NULL,"/");
                    }
                tempName = NULL;
                strcpy(tempName3, tempName2);
                nextStartTime=strtok(tempName3, "-");
                nextStartTime=strtok(NULL, "-");
                nextStartTime=strtok(NULL, "-");
                if(i>=numTmplts-1&&bankDuration>2048)
                  {
                  incrementL1=2049-(bankTemp2->end_time.gpsSeconds-atoi(nextStartTime));
                  }
                nextStartTime=NULL;
                fprintf(stdout, "The incrementL1 is %i\n", incrementL1);
                if( !strcmp( tempName2, namearray[1] ) )
                    {
                    /* if this file has been used, increment the chan # */
                    lChanNum++;
                    }
                else
                    {
                    lChanNum = 0;
                    strcpy( namearray[1], tempName2 );
                    }
                strcpy(namearray2[1],"LivingstonBeamCoeff.dat");
                LALSnprintf( namearray3[1], LALNameLength*sizeof(CHAR), "L1:%s_CData_%d",frInType ,lChanNum );
                LALSnprintf( namearray4[1], LALNameLength*sizeof(CHAR), "_SegNorm_%d", lChanNum );
                LAL_CALL( LALDestroyFrCache(&status, &tempCache), &status );
                tempCache = NULL;
                }

            if( caseID[2] )
                {
                memset( &tempSieve, 0, sizeof(FrCacheSieve) );
                tempSieve.srcRegEx = "V1";
                tempSieve.dscRegEx = networkName;;
                tempSieve.urlRegEx = NULL;
                tempSieve.earliestTime = bankTemp2->end_time.gpsSeconds+incrementV1;
                tempSieve.latestTime = bankTemp2->end_time.gpsSeconds + incrementV1;
                LAL_CALL( LALFrCacheSieve( &status, &tempCache, 
				       frInCache, &tempSieve ), &status );
                if( !tempCache->numFrameFiles )
                    {
                    fprintf(stderr,"The total number of templates returned was %i\n",tempCache->numFrameFiles);
                    fprintf(stderr,"CacheSieve should only have returned 1 frame file for V1. Exiting..\n");
                    goto cleanexit;
                    }
                if(tempCache->numFrameFiles==1) framesPerDetector=1;
                 /* Need to strip the file name out of the url */
                char *tempName = NULL;
                char tempName2[256];
                char tempName3[256];
                char *nextStartTime =NULL;
                tempName = strtok(tempCache->frameFiles->url,"//");
                tempName = strtok(NULL,"/");
                while( tempName != NULL)
                    {
                    strcpy(tempName2,tempName);
                    tempName = strtok(NULL,"/");
                    }
                tempName = NULL;
                strcpy(tempName3, tempName2);
                nextStartTime=strtok(tempName3, "-");
                nextStartTime=strtok(NULL, "-");
                nextStartTime=strtok(NULL, "-");
                if(i>=numTmplts-1&&bankDuration>2048)
                  {
                  incrementV1=2049-(bankTemp2->end_time.gpsSeconds-atoi(nextStartTime));
                  }
                nextStartTime=NULL;
                fprintf(stdout, "The incrementV1 is %i\n", incrementV1);
                if( !strcmp( tempName2, namearray[2] ) )
                    {
                    /* if this file has been used, increment the chan # */
                    virgoChanNum++;
                    }
                else
                    {
                    virgoChanNum = 0;
                    strcpy( namearray[2], tempName2 );
                    }
                strcpy(namearray2[2],"VirgoBeamCoeff.dat");
                LALSnprintf( namearray3[2], LALNameLength*sizeof(CHAR), "V1:%s_CData_%d",frInType, virgoChanNum );
                LALSnprintf( namearray4[2], LALNameLength*sizeof(CHAR), "_SegNorm_%d", virgoChanNum );
                LAL_CALL( LALDestroyFrCache(&status, &tempCache), &status );
                tempCache = NULL;
    		    }

		    if( caseID[3] )
                {
                memset( &tempSieve, 0, sizeof(FrCacheSieve) );
                tempSieve.srcRegEx = "G1";
                tempSieve.dscRegEx = networkName;;
                tempSieve.urlRegEx = NULL;
                tempSieve.earliestTime = bankTemp2->end_time.gpsSeconds+incrementG1;
                tempSieve.latestTime = bankTemp2->end_time.gpsSeconds + incrementG1;
                LAL_CALL( LALFrCacheSieve( &status, &tempCache, 
				       frInCache, &tempSieve ), &status );
                if( !tempCache->numFrameFiles )
                    {
                    fprintf(stderr,"The total number of templates returned was %i\n",tempCache->numFrameFiles);
                    fprintf(stderr,"CacheSieve should only have returned 1 frame file for G1. Exiting..\n");
                    goto cleanexit;
                    }
                if(tempCache->numFrameFiles==1) framesPerDetector=1;
                /* Need to strip the file name out of the url */
                char *tempName = NULL;
                char tempName2[256];
                char tempName3[256];
                char *nextStartTime =NULL;
                tempName = strtok(tempCache->frameFiles->url,"//");
                tempName = strtok(NULL,"/");
                while( tempName != NULL)
                    {
                    strcpy(tempName2,tempName);
                    tempName = strtok(NULL,"/");
                    }
                tempName = NULL;
                strcpy(tempName3, tempName2);
                nextStartTime=strtok(tempName3, "-");
                nextStartTime=strtok(NULL, "-");
                nextStartTime=strtok(NULL, "-");
                if(i>=numTmplts-1&&bankDuration>2048)
                  {
                  incrementG1=2049-(bankTemp2->end_time.gpsSeconds-atoi(nextStartTime));
                  }
                nextStartTime=NULL;
                fprintf(stdout, "The incrementG1 is %i\n", incrementG1);
                if( !strcmp(tempName2, namearray[3]) )
                    {
                    /* if this file has been used, increment the chan # */
                    geoChanNum++;
                    }
                else
                    {
                    geoChanNum = 0;
                    strcpy( namearray[3], tempName2 );
                    }
                strcpy(namearray2[3],"GeoBeamCoeff.dat");
                LALSnprintf( namearray3[3], LALNameLength*sizeof(CHAR), "G1:DER_DATA_H_CData_%d", geoChanNum );
                LALSnprintf( namearray4[3], LALNameLength*sizeof(CHAR), "_SegNorm_%d", geoChanNum );
                LAL_CALL( LALDestroyFrCache(&status, &tempCache), &status );
                tempCache = NULL;
                }

            if( caseID[4] )
                {
                memset( &tempSieve, 0, sizeof(FrCacheSieve) );
                tempSieve.srcRegEx = "T1";
                tempSieve.dscRegEx = networkName;;
                tempSieve.urlRegEx = NULL;
                tempSieve.earliestTime = bankTemp2->end_time.gpsSeconds+incrementT1;
                tempSieve.latestTime = bankTemp2->end_time.gpsSeconds + incrementT1;
                LAL_CALL( LALFrCacheSieve( &status, &tempCache, 
    				       frInCache, &tempSieve ), &status );
                if( !tempCache->numFrameFiles )
                    {
                    fprintf(stderr,"The total number of templates returned was %i\n",tempCache->numFrameFiles);
                    fprintf(stderr,"CacheSieve should only have returned 1 frame file for T1. Exiting..\n");
                    goto cleanexit;
                    }
                if(tempCache->numFrameFiles==1) framesPerDetector=1;
                 /* Need to strip the file name out of the url */
                char *tempName = NULL;
                char tempName2[256];
                char tempName3[256];
                char *nextStartTime =NULL;
                tempName = strtok(tempCache->frameFiles->url,"//");
                tempName = strtok(NULL,"/");
                while( tempName != NULL)
                    {
                    strcpy(tempName2,tempName);
                    tempName = strtok(NULL,"/");
                    }
                tempName = NULL;
                strcpy(tempName3, tempName2);
                nextStartTime=strtok(tempName3, "-");
                nextStartTime=strtok(NULL, "-");
                nextStartTime=strtok(NULL, "-");
                if(i>=numTmplts-1&&bankDuration>2048)
                  {
                  incrementT1=2049-(bankTemp2->end_time.gpsSeconds-atoi(nextStartTime));
                  }
                nextStartTime=NULL;
                fprintf(stdout, "The incrementT1 is %i\n", incrementT1);
                if( !strcmp(tempName2, namearray[4]) )
                    {
                    /* if this file has been used, increment the chan # */
                    tamaChanNum++;
                    }
                else
                    {
                    tamaChanNum = 0;
                    strcpy( namearray[4], tempName2 );
                    }
                strcpy(namearray2[4],"TamaBeamCoeff.dat");
                LALSnprintf( namearray3[4], LALNameLength*sizeof(CHAR), "T1:%s_CData_%d", frInType,tamaChanNum );
                LALSnprintf( namearray4[4], LALNameLength*sizeof(CHAR), "_SegNorm_%d", tamaChanNum );
                LAL_CALL( LALDestroyFrCache(&status, &tempCache), &status );
                tempCache = NULL;
                }

		    if( caseID[5] )
                {
                memset( &tempSieve, 0, sizeof(FrCacheSieve) );
                tempSieve.srcRegEx = "H2";
                tempSieve.dscRegEx = networkName;;
                tempSieve.urlRegEx = NULL;
                tempSieve.earliestTime = bankTemp2->end_time.gpsSeconds+incrementH2;
                tempSieve.latestTime = bankTemp2->end_time.gpsSeconds + incrementH2;
                LAL_CALL( LALFrCacheSieve( &status, &tempCache, 
    				       frInCache, &tempSieve ), &status );
  		        if( !tempCache->numFrameFiles )
                    {
                    fprintf(stderr,"The total number of templates returned was %i\n",tempCache->numFrameFiles);
                    fprintf(stderr,"CacheSieve should only have returned 1 frame file for H2. Exiting..\n");
                    goto cleanexit;
                    }
                if(tempCache->numFrameFiles==1) framesPerDetector=1;
                fprintf(stdout,"The first URL for H2 is %s and the number of frames in the cache is %i\n",tempCache->frameFiles->url,tempCache->numFrameFiles);
                if(tempCache->numFrameFiles>1)
                    {
                    tempFrame=tempCache->frameFiles+1;
                    strcpy(nextURL,tempFrame->url);
                    fprintf(stdout,"The second URL for H2 is %s\n",nextURL);
                    }
                    
                 /* Need to strip the file name out of the url */
                char *tempName = NULL;
                char tempName2[256];
                char tempName3[256];
                char *nextStartTime =NULL;
                tempName = strtok(tempCache->frameFiles->url,"//");
                tempName = strtok(NULL,"/");
                while( tempName != NULL)
                    {
                    strcpy(tempName2,tempName);
                    tempName = strtok(NULL,"/");
                    }
                tempName = NULL;
                strcpy(tempName3, tempName2);
                nextStartTime=strtok(tempName3, "-");
                nextStartTime=strtok(NULL, "-");
                nextStartTime=strtok(NULL, "-");
                if(i>=numTmplts-1&&bankDuration>2048)
                  {
                  incrementH2=2049-(bankTemp2->end_time.gpsSeconds-atoi(nextStartTime));
                  }
                nextStartTime=NULL;
                fprintf(stdout, "The incrementH2 is %i\n", incrementH2);
                if( !strcmp(tempName2, namearray[5]) )
                    {
                    /* if this file has been used, increment the chan # */
                    h2ChanNum++;
                    }
                else
                    {
                    h2ChanNum = 0;/*change back to 0*/
                    strcpy( namearray[5], tempName2 );
                    }
                strcpy(namearray2[5],"HanfordBeamCoeff.dat");
                LALSnprintf( namearray3[5], LALNameLength*sizeof(CHAR), "H2:%s_CData_%d",frInType, h2ChanNum );
                LALSnprintf( namearray4[5], LALNameLength*sizeof(CHAR), "_SegNorm_%d", h2ChanNum );
                LAL_CALL( LALDestroyFrCache(&status, &tempCache), &status );
                tempCache = NULL;
                tempFrame = NULL;
                }

		    chanNumArray[0] = h1ChanNum;
		    chanNumArray[1] = lChanNum;
	        chanNumArray[2] = virgoChanNum;
	        chanNumArray[3] = geoChanNum;
	        chanNumArray[4] = tamaChanNum;
	        chanNumArray[5] = h2ChanNum;

              }
	        else
		      {
              /* If we arent globbing, names need to be set differently*/

              if(caseID[0])
                {
                strcpy(namearray[0],H1filename);
                strcpy(namearray2[0],"HanfordBeamCoeff.dat");
                LALSnprintf( namearray3[0], LALNameLength*sizeof(CHAR), "H1:%s_CData_%d", frInType,h1ChanNum );
                LALSnprintf( namearray4[0], LALNameLength*sizeof(CHAR), "_SegNorm_%d", h1ChanNum );
                }

              if(caseID[1])
                {
                strcpy(namearray[1],Lfilename);
                strcpy(namearray2[1],"LivingstonBeamCoeff.dat");
                LALSnprintf( namearray3[1], LALNameLength*sizeof(CHAR), "L1:%s_CData_%d", frInType,lChanNum );
		        LALSnprintf( namearray4[1], LALNameLength*sizeof(CHAR), "_SegNorm_%d", lChanNum );
                }

              if(caseID[2])
	            {
		        strcpy(namearray[2],VIRGOfilename);
                strcpy(namearray2[2],"VirgoBeamCoeff.dat"); 
                LALSnprintf(namearray3[2], LALNameLength*sizeof(CHAR), "V1:%s_CData_%d", frInType,virgoChanNum );
                LALSnprintf( namearray4[2], LALNameLength*sizeof(CHAR), "_SegNorm_%d", virgoChanNum );
                }

	          if(caseID[3])
                {
                strcpy(namearray[3],GEOfilename);
                strcpy(namearray2[3],"GeoBeamCoeff.dat");
                LALSnprintf(namearray3[3], LALNameLength*sizeof(CHAR), "G1:DER_DATA_H_CData_%d",geoChanNum );
                LALSnprintf( namearray4[3], LALNameLength*sizeof(CHAR), "_SegNorm_%d", geoChanNum );
                }

              if(caseID[4])
                { 
                strcpy(namearray[4],TAMAfilename);
                strcpy(namearray2[4],"TamaBeamCoeff.dat");
                LALSnprintf(namearray3[4], LALNameLength*sizeof(CHAR), "T1:%s_CData_%d", frInType,tamaChanNum );
                LALSnprintf( namearray4[4], LALNameLength*sizeof(CHAR), "_SegNorm_%d", tamaChanNum );
                }

              if(caseID[5])
                {
                strcpy(namearray[5],H2filename);
                strcpy(namearray2[5],"HanfordBeamCoeff.dat");
                LALSnprintf(namearray3[5], LALNameLength*sizeof(CHAR), "H2:%s_CData_%d", frInType,h2ChanNum );
                LALSnprintf( namearray4[5], LALNameLength*sizeof(CHAR), "_SegNorm_%d", h2ChanNum );
                }

		      chanNumArray[0] = h1ChanNum;
		      chanNumArray[1] = lChanNum;
	          chanNumArray[2] = virgoChanNum;
	          chanNumArray[3] = geoChanNum;
	          chanNumArray[4] = tamaChanNum;
	          chanNumArray[5] = h2ChanNum;

		      }/*end else line 786*/

            if(verbose) fprintf(stdout,"built namearrays\n");

            if(verbose) 
              {
              fprintf(stdout,"frame files: %s\n %s\n %s\n %s\n %s\n %s\n",namearray[0],namearray[1],namearray[2],namearray[3],namearray[4],namearray[5]);
              fprintf(stdout,"channels: %s\n %s\n %s\n %s\n %s\n %s\n",namearray3[0],namearray3[1],namearray3[2],namearray3[3],namearray3[4],namearray3[5]);
		      }

            /* get beam pattern coefficients if necessary */
	        if ( (numDetectors == 3 && !( caseID[0] && caseID[5])) || numDetectors == 4 )
                {
                cohInspBeamVec = cohInspFilterInput->beamVec;
                w=0;
                if( verbose ) fprintf(stdout, "This network requires beam-pattern coefficients - reading them in...\n");
                for ( j=0; j<6; j++ ) {
                    if ( caseID[j] ) 
                        { 
                        filePtr[w] = fopen(namearray2[j], "r");
                        if(!filePtr[w])
                            {
                            fprintf(stderr,"The file %s containing the coefficients could not be found - exiting...\n",namearray2[j]);
                            goto cleanexit;
                            }
                        for ( l=0 ; l < (INT4) numBeamPoints ; l++)
                            { 
            			    fscanf(filePtr[w],"%f %f %f %f",&theta,&phi,&vPlus,&vMinus);
                            cohInspBeamVec->detBeamArray[w].thetaPhiVs[l].data->data[0] = theta;
                            cohInspBeamVec->detBeamArray[w].thetaPhiVs[l].data->data[1] = phi;
                            cohInspBeamVec->detBeamArray[w].thetaPhiVs[l].data->data[2] = vPlus;
                            cohInspBeamVec->detBeamArray[w].thetaPhiVs[l].data->data[3] = vMinus;
                            }
                        fclose(filePtr[w++]);
                        }
                    } 
                }
	        cohInspCVec = cohInspFilterInput->multiCData;

            /* Read in the snippets associated with this event */
	        l = 0;
	        for( j=0; j<6; j++ )
                {
                if( caseID[j] )
                    {
                    if( verbose ) fprintf(stdout, "getting the COMPLEX8TimeSeries\n");
                      LAL_CALL( LALFrOpen(&status,&frStream,NULL,namearray[j]), &status);
                    if(!frStream)
                        {
                        fprintf(stdout,"The file %s does not exist - exiting...\n", namearray[j]);
                        goto cleanexit;
                        }
                    /* Since there is a segnorm frame written at the front */
                    /* of the frame file for each c-data frame, I need to  */
                    /* determine the number of frames contained in the     */
                    /* frame file so that I can advance the stream past    */
                    /* the segNorm vectors                                 */

                    p = 0;
                    frEnd = 0;
                    while( !frEnd )
                        {
                        LAL_CALL( LALFrEnd( &status, &frEnd, frStream), &status);
                        if( frEnd == 0 )
                            {
                            LAL_CALL( LALFrNext( &status,frStream ), &status);
                            p++;
                            }
                        LAL_CALL( LALFrEnd( &status, &frEnd, frStream), &status);
                        }
                    numChannels = p;
                    LAL_CALL( LALFrRewind( &status, frStream ), &status );

                    /* get segnorm */
			
                    for( w=0; w<chanNumArray[j]; w++)
                        {
                        LAL_CALL( LALFrNext( &status,frStream ), &status );
                        }
                    frChan.name = namearray4[j];
                    segNormVector = (REAL4FrequencySeries *) 
                    LALCalloc( 1, sizeof(REAL4FrequencySeries) );
                    LAL_CALL( LALFrGetREAL4FrequencySeries( &status, 
                            segNormVector, &frChan, frStream), &status);

                    REAL4 fFinal = 1.0 / (6.0 * sqrt(6.0) * LAL_PI * cohInspFilterInput->tmplt->totalMass * LAL_MTSUN_SI);
                    REAL4 deltaF = 1.0 / ( cohInspFilterParams->deltaT * (REAL4) numPointsSeg );
                    INT4 kmax = fFinal / deltaF < numPointsSeg/2 ? 
                                fFinal / deltaF : numPointsSeg/2;
                    cohInspFilterParams->segNorm[l] = segNormVector->data->data[kmax];
                    LAL_CALL( LALDestroyVector( &status, &(segNormVector->data) ), &status );
                    LALFree( segNormVector );
                    segNormVector = NULL;
                    LAL_CALL( LALFrRewind( &status, frStream ), &status );
			

                    /* Advance the stream to the appropriate frame */
                    for( w=0; w<chanNumArray[j]+numChannels/2; w++)
                        {
                        LAL_CALL( LALFrNext( &status,frStream ), &status );
                        }
                    frChan.name = namearray3[j];
                    LAL_CALL( LALFrGetCOMPLEX8TimeSeries( &status, &(cohInspCVec->cData[l]), &frChan, frStream), &status);
                    /* tempTime is the start time of cData plus -(time slide)*/
                    tempTime[l] += cohInspCVec->cData[l].epoch.gpsSeconds + cohInspCVec->cData[l].epoch.gpsNanoSeconds * 1e-9;
                    if( verbose ) fprintf(stdout,"tempTime = %f\n",tempTime[l]);
                    LAL_CALL( LALFrClose( &status, &frStream ), &status );
                    l++;
                    }/*end line 888*/
                }/*end line 886*/

	        /* If we can estimate distance then compute templateNorm */
            /* At present, this is only good for frequency domain tmplts */
              /* Since each detectors data has been filtered with a templates*/
              /* that have the same mass pair, templateNorm is the same for */
              /* every detector and needs to be computed only once.         */

            REAL4 cannonDist = 1.0; /* Mpc */
		    REAL4 m  = cohInspFilterInput->tmplt->totalMass;
		    REAL4 mu = cohInspFilterInput->tmplt->mu;
		    REAL4 deltaT = cohInspFilterParams->deltaT;
		    REAL4 distNorm = 2.0 *  LAL_MRSUN_SI / (cannonDist * 1e6 * LAL_PC_SI);
            REAL4 templateNorm = sqrt( (5.0*mu) / 96.0 ) *  pow( m / (LAL_PI*LAL_PI) , 1.0/3.0 ) * pow( LAL_MTSUN_SI / deltaT, -1.0/6.0 );
		    distNorm *= dynRange;
		    templateNorm *= templateNorm;
		    templateNorm *= distNorm * distNorm;
		    cohInspFilterParams->templateNorm = templateNorm;
		    cohInspFilterParams->segmentLength = numPointsSeg;
		
	        if (verbose) fprintf(stdout,"filtering the data..\n");
	        if ( maximizeOverChirp && verbose )
                {
                fprintf(stdout,"clustering events\n");
                }		     

	        /* Before the data gets filtered, I need to make the c-data snippets commensurate */
	        for(j=0;j<(INT4)numDetectors - 1;j++)
                {
                timeptDiff[j] = rint((tempTime[0] - tempTime[j+1]) * sampleRate);
                if( verbose ) fprintf(stdout,"timeptDiff = %d\n",timeptDiff[j]);
                }

	        /* Now allocate memory for a temporary storage vector */
	        memset( &tempSnippet, 0, sizeof(COMPLEX8TimeSeries) );
	        LAL_CALL( LALCCreateVector( &status, &(tempSnippet.data), numPoints ), &status );

	        /* The following switch statement accomplishes the commensuration of the time series */
	        switch( numDetectors )
                {
                case 2:
                    if( timeptDiff[0] < 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[1].data->data, numPoints * sizeof(COMPLEX8) );
            			if( verbose ) 
                            fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[0].re,tempSnippet.data->data[0].im,tempSnippet.data->data[1].re,tempSnippet.data->data[1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[1].data->data[j].re = 0;
                            cohInspCVec->cData[1].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[1].data->data + abs(timeptDiff[0]), tempSnippet.data->data, (numPoints - abs(timeptDiff[0])) * sizeof(COMPLEX8) );
                        if( verbose ) 
                            fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[1].data->data[abs(timeptDiff[0])].re, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])].im, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])+1].re, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])+1].im);
                        cohInspCVec->cData[1].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[1].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else if( timeptDiff[0] > 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[1].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[abs(timeptDiff[0])].re,tempSnippet.data->data[abs(timeptDiff[0])].im,tempSnippet.data->data[abs(timeptDiff[0])+1].re,tempSnippet.data->data[abs(timeptDiff[0])+1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[1].data->data[j].re = 0;
                            cohInspCVec->cData[1].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[1].data->data, tempSnippet.data->data + abs( timeptDiff[0] ), (numPoints - abs(timeptDiff[0])) * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[1].data->data[0].re, cohInspCVec->cData[1].data->data[0].im, cohInspCVec->cData[1].data->data[1].re, cohInspCVec->cData[1].data->data[1].im);
                        cohInspCVec->cData[1].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[1].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else
                        {
                        if( verbose ) fprintf(stdout,"Time series are commensurate...directly combining them.\n");
                        }
                    break;
		        case 3:
                    if( timeptDiff[0] < 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[1].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[0].re,tempSnippet.data->data[0].im,tempSnippet.data->data[1].re,tempSnippet.data->data[1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[1].data->data[j].re = 0;
                            cohInspCVec->cData[1].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[1].data->data + abs(timeptDiff[0]), tempSnippet.data->data, (numPoints - abs(timeptDiff[0])) * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[1].data->data[abs(timeptDiff[0])].re, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])].im, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])+1].re, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])+1].im);
                        cohInspCVec->cData[1].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[1].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else if( timeptDiff[0] > 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[1].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[0].re,tempSnippet.data->data[0].im,tempSnippet.data->data[1].re,tempSnippet.data->data[1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[1].data->data[j].re = 0;
                            cohInspCVec->cData[1].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[1].data->data, tempSnippet.data->data + abs( timeptDiff[0] ), (numPoints - abs(timeptDiff[0])) * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[1].data->data[abs(timeptDiff[0])].re, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])].im, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])+1].re, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])+1].im);
                        cohInspCVec->cData[1].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[1].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else
                        {
                        if( verbose ) fprintf(stdout,"Time series are commensurate...directly combining them.\n");
                        }

                    if( timeptDiff[1] < 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[2].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[0].re,tempSnippet.data->data[0].im,tempSnippet.data->data[1].re,tempSnippet.data->data[1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[2].data->data[j].re = 0;
                            cohInspCVec->cData[2].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[2].data->data + abs(timeptDiff[1]), tempSnippet.data->data, (numPoints - abs(timeptDiff[1])) * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[2].data->data[abs(timeptDiff[1])].re, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])].im, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])+1].re, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])+1].im);
                        cohInspCVec->cData[2].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[2].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else if( timeptDiff[1] > 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[2].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[0].re,tempSnippet.data->data[0].im,tempSnippet.data->data[1].re,tempSnippet.data->data[1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[2].data->data[j].re = 0;
                            cohInspCVec->cData[2].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[2].data->data, tempSnippet.data->data + abs( timeptDiff[1] ), (numPoints - abs(timeptDiff[1])) * sizeof(COMPLEX8) );
                        if( verbose )fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[2].data->data[abs(timeptDiff[1])].re, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])].im, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])+1].re, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])+1].im);
    
                        cohInspCVec->cData[2].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[2].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else
                        {
                        if( verbose ) fprintf(stdout,"Time series are commensurate...directly combining them.\n");
                        }		   
                    break;
		        case 4:
                    if( timeptDiff[0] < 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[1].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[0].re,tempSnippet.data->data[0].im,tempSnippet.data->data[1].re,tempSnippet.data->data[1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[1].data->data[j].re = 0;
                            cohInspCVec->cData[1].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[1].data->data + abs(timeptDiff[0]), tempSnippet.data->data, (numPoints - abs(timeptDiff[0])) * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[1].data->data[abs(timeptDiff[0])].re, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])].im, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])+1].re, cohInspCVec->cData[1].data->data[abs(timeptDiff[0])+1].im);
                        cohInspCVec->cData[1].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[1].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else if( timeptDiff[0] > 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[1].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[abs(timeptDiff[0])].re,tempSnippet.data->data[abs(timeptDiff[0])].im,tempSnippet.data->data[abs(timeptDiff[0])+1].re,tempSnippet.data->data[abs(timeptDiff[0])+1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[1].data->data[j].re = 0;
                            cohInspCVec->cData[1].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[1].data->data, tempSnippet.data->data + abs( timeptDiff[0] ), (numPoints - abs(timeptDiff[0])) * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[1].data->data[0].re, cohInspCVec->cData[1].data->data[0].im, cohInspCVec->cData[1].data->data[1].re, cohInspCVec->cData[1].data->data[1].im);
                        cohInspCVec->cData[1].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[1].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else
                        {
                        if( verbose ) fprintf(stdout,"Time series are commensurate...directly combining them.\n");
                        }


                    if( timeptDiff[1] < 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[2].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[0].re,tempSnippet.data->data[0].im,tempSnippet.data->data[1].re,tempSnippet.data->data[1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[2].data->data[j].re = 0;
                            cohInspCVec->cData[2].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[2].data->data + abs(timeptDiff[1]), tempSnippet.data->data, (numPoints - abs(timeptDiff[1])) * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[2].data->data[abs(timeptDiff[1])].re, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])].im, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])+1].re, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])+1].im);
            			cohInspCVec->cData[2].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[2].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else if( timeptDiff[1] > 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[2].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[0].re,tempSnippet.data->data[0].im,tempSnippet.data->data[1].re,tempSnippet.data->data[1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[2].data->data[j].re = 0;
                            cohInspCVec->cData[2].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[2].data->data, tempSnippet.data->data + abs( timeptDiff[1] ), (numPoints - abs(timeptDiff[1])) * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[2].data->data[abs(timeptDiff[1])].re, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])].im, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])+1].re, cohInspCVec->cData[2].data->data[abs(timeptDiff[1])+1].im);
                        cohInspCVec->cData[2].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[2].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else
                        {
                        if( verbose ) fprintf(stdout,"Time series are commensurate...directly combining them.\n");
                        }

                    if( timeptDiff[2] < 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[3].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[0].re,tempSnippet.data->data[0].im,tempSnippet.data->data[1].re,tempSnippet.data->data[1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[3].data->data[j].re = 0;
                            cohInspCVec->cData[3].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[3].data->data + abs(timeptDiff[2]), tempSnippet.data->data, (numPoints - abs(timeptDiff[2])) * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[3].data->data[abs(timeptDiff[2])].re, cohInspCVec->cData[3].data->data[abs(timeptDiff[2])].im, cohInspCVec->cData[3].data->data[abs(timeptDiff[2])+1].re, cohInspCVec->cData[3].data->data[abs(timeptDiff[2])+1].im);
                        cohInspCVec->cData[3].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
                        cohInspCVec->cData[3].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else if( timeptDiff[2] > 0 )
                        {
                        memcpy( tempSnippet.data->data, cohInspCVec->cData[3].data->data, numPoints * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"tempSnippet data: %f %f %f %f\n",tempSnippet.data->data[0].re,tempSnippet.data->data[0].im,tempSnippet.data->data[1].re,tempSnippet.data->data[1].im);
                        for(j=0;j<(INT4)numPoints;j++)
                            {
                            cohInspCVec->cData[3].data->data[j].re = 0;
                            cohInspCVec->cData[3].data->data[j].im = 0;
                            }
                        memcpy( cohInspCVec->cData[3].data->data, tempSnippet.data->data + abs( timeptDiff[2] ), (numPoints - abs(timeptDiff[2])) * sizeof(COMPLEX8) );
                        if( verbose ) fprintf(stdout,"Some frame data: %f %f %f %f\n",cohInspCVec->cData[3].data->data[abs(timeptDiff[2])].re, cohInspCVec->cData[3].data->data[abs(timeptDiff[2])].im, cohInspCVec->cData[3].data->data[abs(timeptDiff[2])+1].re, cohInspCVec->cData[3].data->data[abs(timeptDiff[2])+1].im);
            
            			cohInspCVec->cData[3].epoch.gpsSeconds = cohInspCVec->cData[0].epoch.gpsSeconds;
            			cohInspCVec->cData[3].epoch.gpsNanoSeconds = cohInspCVec->cData[0].epoch.gpsNanoSeconds;
                        }
                    else
                        {
                        if( verbose ) fprintf(stdout,"Time series are commensurate...directly combining them.\n");
                        }		   
                    break;
                }/*end switch line 992*/

	        /* Now that the time series are commensurate, do the filtering... */
            cohInspFilterParams->cohSNRVec->epoch = cohInspCVec->cData[0].epoch;  
	        LALCoherentInspiralFilterSegment (&status, &thisEvent, cohInspFilterInput, cohInspFilterParams);

            /* Now change from geocentric to equatorial */
	        thisEventTemp =thisEvent;
	        if ( (numDetectors == 3 && !( caseID[0] && caseID[5])) || numDetectors == 4 )
                {
                while( thisEventTemp )
                    {
                    convertParams.system=COORDINATESYSTEM_EQUATORIAL;
                    convertParams.zenith= NULL;
                    convertParams.gpsTime=&(thisEventTemp->end_time);
                    /*The next few lines gets the parameters ready for conversion. Need to check to make sure that thisEvent->ligo_axis_dec is in proper notation, i.e. */
    
                    tempSky.system=COORDINATESYSTEM_GEOGRAPHIC;
                    tempSky.latitude=thisEventTemp->ligo_axis_dec;
                    tempSky.longitude=thisEventTemp->ligo_axis_ra;
                    LALConvertSkyCoordinates(&status,&tempSky, &tempSky,&convertParams);	     
                    thisEventTemp->ligo_axis_dec=LAL_PI*0.5-tempSky.latitude;
                    thisEventTemp->ligo_axis_ra=tempSky.longitude;
                    thisEventTemp = thisEventTemp->next;
                    }
                }
	        thisEventTemp =thisEvent;
	        while( thisEventTemp )
                {
                thisEventTemp->event_id= (EventIDColumn *)LALCalloc(1, sizeof(EventIDColumn) ); 
                thisEventTemp->event_id->id=eventID;
                thisEventTemp = thisEventTemp->next;
                }

	        if ( cohSNROut )
                {
                LALSnprintf( cohdataStr, LALNameLength*sizeof(CHAR),
			       "SNR_%d", nCohDataFr++ );
                strcpy( cohInspFilterParams->cohSNRVec->name, "Coherent");
                outFrame = fr_add_proc_REAL4TimeSeries( outFrame, cohInspFilterParams->cohSNRVec, "none", cohdataStr );
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


	        if( thisEvent )
                {
                if( verbose ) fprintf( stdout,"******> Dumping Events <******\n");
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

                }/*end if( thisEvent ) around line 1253*/

            /* initialize */
            LAL_CALL( LALCDestroyVector( &status, &(tempSnippet.data) ), &status );
	        LALFree( cohInspFilterInput->tmplt );
	        cohInspFilterInput->tmplt = NULL;
	        free( cohInspInitParams );

	        /* Destroy params structure for coherent filter code */
	        LAL_CALL( LALCoherentInspiralFilterParamsFinalize (&status,&cohInspFilterParams), &status );
            TestStatus (&status, "0", 1);
	        ClearStatus (&status);

            /* Destroy input structure for coherent filter code */
	        LAL_CALL( LALCoherentInspiralFilterInputFinalize (&status, &cohInspFilterInput), &status);
            TestStatus (&status, "0", 1);
	        ClearStatus (&status);
 
    	    k = -1;
	        for(j=0;j<4;j++)
              {
              memcpy(caseIDChars[j], "0",3);
              }
	        for(j=0;j<6;j++)
              {
              caseID[j] = 0;
              if( !globFrameData )
    		    {
                chanNumArray[j] = -1; 
                strcpy( namearray[j], "0" );
                }
              strcpy( namearray2[j], "0" );
              strcpy( namearray3[j], "0" );
              }
	        for(j=0;j<5;j++)
	          {
              tempTime[j] = 0.0;
              timeptDiff[j] = 0;
              }
	        tempTime[5] = 0.0;
            }/*end if( k > 0 ) around line 361*/
	  
          else
	        {
	        fprintf(stderr,"Error - there is only one detector associated with event %Ld, aborting...\n", eventID );
	        exit( 1 );
            }
	    }/*end if( i != 0 && ((bankTemp->event_id->id != eventID) || (i == numTmplts - 1))) around line 359*/
        k++;    
    }/*end for loop over number of templates*/
 
    framesPerDetector--;
    /*if(framesPerDetector<1) 
      {
      incrementL1=1;
      incrementT1=1;
      incrementG1=1;
      incrementV1=1;
      incrementH2=1;
      incrementH1=1;
      }*/
               
    /*}*/
   /*end while(framesPerDetector>=1)*/
  
 /* set the name of the output file(s) (without extension) */

  if ( userTag )
    {
    LALSnprintf( fileName, FILENAME_MAX, "%s-COHERENT_%s-%d-%d", ifoTag, userTag, 
        gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds ); 
    }
  else
    {
    LALSnprintf( fileName, FILENAME_MAX, "%s-COHERENT-%d-%d", ifoTag, 
        gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds ); 
    }


  if( cohSNROut )
    {
    if ( outputPath[0] )
        {
        LALSnprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s/%s.gwf", outputPath, fileName);
        }
    else 
        {
        LALSnprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s.gwf", fileName );
        }

    if ( verbose ) fprintf( stdout, "writing frame data to %s....", framename );
    frOutFile = FrFileONew( framename, 0);
    FrameWrite( outFrame, frOutFile);
    FrFileOEnd( frOutFile );
    if ( verbose ) fprintf(stdout, "done\n");
  
    }
    

  
  if ( eventsOut )
    { 
    memset( &results, 0, sizeof(LIGOLwXMLStream) );
    if ( outputPath[0] )
        {
        LALSnprintf( xmlname, FILENAME_MAX * sizeof(CHAR), "%s/%s.xml", outputPath, fileName);
        }
    else 
        {
        LALSnprintf( xmlname, FILENAME_MAX * sizeof(CHAR), "%s.xml", fileName );
        }
    if ( verbose ) fprintf( stdout, "writing XML data to %s...\n", xmlname );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, xmlname), &status );

    /* write the process table */
    if ( verbose ) fprintf( stdout, "  process table...\n" );
    /*      LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", caseID );*/
    LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time), &accuracy ), &status );
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable, process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
    free( proctable.processTable );

    /* write the process params table */
    if ( verbose ) fprintf( stdout, "  process_params table...\n" );
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams, process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
    while( procparams.processParamsTable )
        {
        this_proc_param = procparams.processParamsTable;
        procparams.processParamsTable = this_proc_param->next;
        free( this_proc_param );
        }

    if( verbose ) fprintf(stdout,"  event params table\n ");

    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, multi_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, savedEvents, multi_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &results), &status );

    while( savedEvents.multiInspiralTable )
	    {  
	    MultiInspiralTable *tempEvent = savedEvents.multiInspiralTable;
	    savedEvents.multiInspiralTable = savedEvents.multiInspiralTable->next;
	    LALFree( tempEvent->event_id );
	    LALFree( tempEvent );
	    tempEvent = NULL;
	    }
    /* close the output xml file */
    LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );
    if ( verbose ) fprintf( stdout, "done. XML file closed\n" );

    }/*end if ( eventsOut ) around line 1389*/

  goto cleanexit;

  cleanexit:

  /* Free the template bank */
  while( bankHead )
    {
      InspiralTemplate *tempTmplt = bankHead;
      bankHead = bankHead->next;
      LALFree( tempTmplt->event_id );
      LALFree( tempTmplt );
      tempTmplt = NULL;
    }

  /* free the frame cache */
  if( frInCache ) LAL_CALL( LALDestroyFrCache( &status, &frInCache ), &status );

  if ( frInType )    free( frInType );

  if ( verbose ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
  exit(0);

}/* main function end */


/* ------------------------------------------------------------------------- */



static void
TestStatus (
    LALStatus  *status, 
    const char *ignored, 
    int         exitcode
           )
{
  char  str[64];
  char *tok;

  if (verbose)
  {
    REPORTSTATUS (status);
  }

  if (strncpy (str, ignored, sizeof (str)))
  {
    if ((tok = strtok (str, " ")))
    {
      do
      {
        if (status->statusCode == atoi (tok))
        {
          return;
        }
      }
      while ((tok = strtok (NULL, " ")));
    }
    else
    {
      if (status->statusCode == atoi (tok))
      {
        return;
      }
    }
  }

  fprintf (stderr, "\nExiting to system with code %d\n", exitcode);
  exit (exitcode);
}


/*
 *
 * ClearStatus ()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
static void
ClearStatus (
    LALStatus   *status
            )
{
  if (status->statusPtr)
  {
    ClearStatus      (status->statusPtr);
    DETATCHSTATUSPTR (status);
  }
}

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
      LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
          long_options[option_index].name ); \
          LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
          LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );


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
"  --cohsnr-threshold RHO       set signal-to-noise threshold to RHO\n"\
"  --cohsnr-threshold RHO       set signal-to-noise threshold to RHO\n"\
"  --maximize-over-chirp        do clustering\n"\
"  --glob-frame-data            glob files in the pwd to obtain frame data\n"\
"  --frame-type TAG             input data is contained in frames of type TAG\n"\
"  --gps-start-time SEC         GPS second of data start time (needed if globbing)\n"\
"  --gps-end-time SEC           GPS second of data end time (needed if globbing)\n"\
"\n"
#define USAGE3 \
"  --write-events               write events\n"\
"  --write-cohsnr               write cohsnr\n"\
"  --output-path                write files here\n"\
"  --H1-framefile               frame data for H1\n"\
"  --H2-framefile               frame data for H2\n"\
"  --L-framefile                frame data for L\n"\
"  --V-framefile                frame data for V\n"\
"  --G-framefile                frame data for G\n"\
"  --T-framefile                frame data for T\n"\
"\n"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
   struct option long_options[] = 
   {
     {"verbose",                  no_argument,       &verbose,           1 },
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
     {"cohsnr-threshold",         required_argument, 0,                 'p'},
     {"maximize-over-chirp",      no_argument,       &maximizeOverChirp, 1 },
     {"glob-frame-data",          no_argument,       &globFrameData,     1 },
     {"write-events",             no_argument,       &eventsOut,         1 },
     {"write-cohsnr",             no_argument,       &cohSNROut,         1 },
     {"gps-start-time",           required_argument, 0,                 'a'},
     {"gps-end-time",             required_argument, 0,                 'b'},
     {"output-path",              required_argument, 0,                 'P'},
     {"H1-framefile",             required_argument, 0,                 'A'},
     {"H2-framefile",             required_argument, 0,                 'Z'},
     {"L-framefile",              required_argument, 0,                 'L'},
     {"V-framefile",              required_argument, 0,                 'V'},
     {"G-framefile",              required_argument, 0,                 'G'},
     {"frame-type",               required_argument, 0,                 'S'},
     {"T-framefile",              required_argument, 0,                 'T'},
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
	   "A:B:a:b:G:S:I:L:l:e:g:W:X:Y:t:w:P:T:V:Z:d:f:h:p:r:u:v:",
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
	       break;	     }
	   else 
	     {
	       fprintf( stderr, "error parsing option %s with argument %s\n",
			long_options[option_index].name, optarg );
	       exit( 1 );
	     }
	   break;

	 case 'A':
	   strcpy(H1filename,optarg);
	   H1file = 1;
	   ADD_PROCESS_PARAM( "string", "%s", H1filename );
	   break;

	 case 'G':
	   strcpy(GEOfilename,optarg);
	   GEOfile = 1;
	   ADD_PROCESS_PARAM( "string", "%s", GEOfilename );
	   break;

	 case 'L':
	   strcpy(Lfilename,optarg);
	   Lfile = 1;
	   ADD_PROCESS_PARAM( "string", "%s", Lfilename );
	   break;

	 case 'T':
	   strcpy(TAMAfilename,optarg);
	   TAMAfile = 1;
	   ADD_PROCESS_PARAM( "string", "%s", TAMAfilename );
	   break;

	 case 'V':
	   strcpy(VIRGOfilename,optarg);
	   VIRGOfile = 1; 
	   ADD_PROCESS_PARAM( "string", "%s", VIRGOfilename );
	   break;

	 case 'Z':
	   strcpy(H2filename,optarg);
	   H2file = 1;
	   ADD_PROCESS_PARAM( "string", "%s", H2filename );
	   break;

         case 'B':
           /* create storaged for the ifo-tag */
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
	   LALSnprintf( outputPath, FILENAME_MAX * sizeof(CHAR),"%s", optarg );
	   ADD_PROCESS_PARAM( "string", "%s", outputPath );
	   break;
     case 'S':
        optarg_len = strlen( optarg ) + 1;
        frInType = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( frInType, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
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
           /* create storage for the bank filename */
	   /*optarg_len = strlen( optarg ) + 1;
	   bankFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	   memcpy( bankFileName, optarg, optarg_len );*/
       strcpy(bankFileName, optarg);
       char tempName[256];
       char *duration =NULL;
       strcpy(tempName, bankFileName);
       duration = strtok(tempName,"-");
       duration = strtok(NULL,"-");
       duration = strtok(NULL,"-");
       duration = strtok(NULL,".");
       bankDuration=atoi(duration);        
	   ADD_PROCESS_PARAM( "string", "%s", bankFileName );
       duration=NULL;
       break;

	   /* Read in time-slide steps for all detectors */
	   /* Read in time-slide step for G1 */
	 case 'g': 
	   slideStep[0] = atof(optarg);
	   ADD_PROCESS_PARAM("float", "%e", slideStep[0] );
	   break;

	   /* Read in time-slide step for H1 */
	 case 'W': 
	   slideStep[1] = atof(optarg);
	   ADD_PROCESS_PARAM("float", "%e", slideStep[1]);
	   break;

	   /* Read in time-slide step for H2 */
	 case 'X': 
	   slideStep[2] = atof(optarg);
	   ADD_PROCESS_PARAM("float", "%e", slideStep[2]);
	   break;

	   /* Read in time-slide step for L1 */
	 case 'Y': 
	   slideStep[3] = atof(optarg);
	   ADD_PROCESS_PARAM("float", "%e", slideStep[3]);
	   break;

	   /* Read in time-slide step for T1 */
	 case 't': 
	   slideStep[4] = atof(optarg);
	   ADD_PROCESS_PARAM("float", "%e", slideStep[4]);
	   break;

	   /* Read in time-slide step for V1 */
	 case 'w': 
	   slideStep[5] = atof(optarg);
	   ADD_PROCESS_PARAM("float", "%e", slideStep[5]);
	   break;

	 case 'v':
	   /* print version information and exit */
           fprintf( stdout, "LIGO/LSC Multi-Detecter Search Code\n" 
                 "Bose/Seader <sukanta@wsu.edu> <sseader@wsu.edu>\n"
                 "CVS Version: " CVS_ID_STRING "\n"
		 "CVS Tag: " CVS_NAME_STRING "\n" );
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

   if ( globFrameData && (H1file || H2file || Lfile || GEOfile || VIRGOfile || TAMAfile) )
     {
       fprintf( stderr, "Specify frame files or to glob for frames - not both\n" );
       exit( 1 );
     }

   /* check validity of input data time if globbing */
   /* the times should be that spanned by the bank(trigger) file */
   if ( ! gpsStartTimeNS )
     {
       fprintf( stderr, "--gps-start-time must be specified\n" );
       exit( 1 );
     }
   LAL_CALL( LALINT8toGPS( &status, &gpsStartTime, &gpsStartTimeNS ),
	     &status );
   if ( ! gpsEndTimeNS )
     {
       fprintf( stderr, "--gps-end-time must be specified\n" );
       exit( 1 );
     }
   LAL_CALL( LALINT8toGPS( &status, &gpsEndTime, &gpsEndTimeNS ), 
	     &status );
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

   if ( ! bankFileName )
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
       LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
     }
  
   return 0;
}

#undef ADD_PROCESS_PARAM


