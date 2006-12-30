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

INT4   sampleRate           = -1;       /* sample rate of filter data   */
INT4   numPointsSeg         = -1;/* set to segment-length used in inspiral.c */
REAL4  fLow                 = -1;       /* low frequency cutoff         */
REAL4  dynRangeExponent     = -1;/* set to same value used in inspiral.c */

/*Coherent code specific inputs*/

/* CHECK: DELETE:
   char   H1filename[256];        
   char   Lfilename[256];         
   char   GEOfilename[256];       
   char   VIRGOfilename[256];      
   char   TAMAfilename[256];       
   char   H2filename[256];        
*/

char   ifoframefile[LAL_NUM_IFO][256];

INT4 H1file = 0;
INT4 H2file = 0;
INT4 L1file = 0;
INT4 G1file = 0;
INT4 T1file = 0;
INT4 V1file = 0;
                                
/* input time-slide parameters */
REAL8  slideStep[6]     = {0.0,0.0,0.0,0.0,0.0,0.0};
int    bankDuration     = 0;
/* CHECK: CHAR   bankFileName[FILENAME_MAX]; */
CHAR  *cohbankFileName = NULL; /* name of input template bank  */
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

/* Params to convert from geocentric to equatorial */

ConvertSkyParams             convertParams;
SkyPosition                  tempSky;
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
  CHAR   caseIDChars[6][LIGOMETA_IFOS_MAX] = {"0","0","0","0","0","0"};
  
  INT4   cohSegLength     = 4; /* This should match hardcoded value in inspiral.c */
  INT4   numPoints        = 0;
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


  CoherentInspiralInitParams   *cohInspInitParams = NULL;
  CoherentInspiralFilterParams *cohInspFilterParams = NULL;
  CoherentInspiralFilterInput  *cohInspFilterInput = NULL;
  CoherentInspiralBeamVector   *cohInspBeamVec = NULL;
  CoherentInspiralCVector      *cohInspCVec = NULL;
  MultiInspiralTable           *thisEvent = NULL;
  MultiInspiralTable           *tempTable = NULL;
  MetadataTable                savedEvents;
  COMPLEX8TimeSeries            tempSnippet;
  
  char nameArrayBeam[6][256] = {"0","0","0","0","0","0"}; /* beam files */
  char nameArrayCData[6][256] = {"0","0","0","0","0","0"}; /* cData chan names */

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
  
  /* Set the dynamic range; needed for distNorm, templateNorm calculation */
  dynRange = pow( 2.0, dynRangeExponent );
  
  /* Set other variables */
  numPoints = sampleRate * cohSegLength;
  savedEvents.multiInspiralTable = NULL;
  k = 0; 
  
  /* Read in the frame files */
  
  if( verbose ) fprintf(stdout, "reading in the frame files\n");
  for ( k=0; k<LAL_NUM_IFO ; k++)
    {
      if( ifoframefile[k] )
	{
	  LAL_CALL( LALFrOpen(&status,&frStream,NULL,ifoframefile[k]), &status);
	  if(!frStream)
	    {
	      fprintf(stdout,"The file %s does not exist - exiting...\n", ifoframefile[k]);
	      goto cleanexit;
	    }
	}
    }
  
  /* read in the cohbank trigger ligo lw xml file */
  numTriggers = XLALReadInspiralTriggerFile( &cohbankEventList,
					     &currentTrigger, &searchSummList, &inputFiles, cohbankFileName );
  
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
  else
    {
      if( vrbflg )
	{
	  fprintf( stdout, 
		   "Read in %d triggers from the file %s\n", numTriggers, 
		   cohbankFileName );
	}
      
      /* reconstruct the coincs */
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
      
      /* loop over coincident triggers to compute cohSNR */
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
		  memcpy( caseIDChars[k], &thisCoinc->snglInspiral[k]->ifo, sizeof(caseIDChars[k] - 1) );
		  eventID = thisCoinc->snglInspiral[k]->event_id;
		  if( verbose ) fprintf(stdout,"eventID = %Ld\n",eventID );
		  
		  /* Parse eventID to get the slide number */
		  triggerNumber = eventID % 100000;
		  slideNumber = ((eventID % 100000000) - triggerNumber)/100000;
		  slideSign = (eventID % 1000000000) - slideNumber*100000 - triggerNumber;
		  
		  /* Store CData frame name now for reading it's frame-file later 
		     within thisCoinc-ident loop*/	  
		  LALSnprintf( nameArrayCData[k], LALNameLength*sizeof(CHAR), "%s:%s_CData_%d", caseIDChars[k],frInType,eventID );
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
	  
	  l = 0;
	  
	  if( verbose ) fprintf(stdout,"numDetectors = %d\n", numDetectors);
	  if( verbose ) fprintf(stdout,"caseID = %d %d %d %d %d %d (G1,H1,H2,L1,V1,T1)\n", caseID[0], caseID[1], caseID[2], caseID[3], caseID[4], caseID[5]);
	  
	  
	  /* Initialize the necessary structures for thisCoinc-ident trigger*/
	  
	  if( !(cohInspInitParams = (CoherentInspiralInitParams *) calloc(1,sizeof(CoherentInspiralInitParams)) ))
	    {
	      fprintf( stdout, "could not allocate memory for coherentInspiral init params\n" );
	      goto cleanexit;
	    }
	  
	  /* Initialize the coherent param structure for thisCoinc trigger */      
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
	  
	  m1 = thisCoinc->snglInspiral[kmax]->mass1;
	  m2 = thisCoinc->snglInspiral[kmax]->mass2;
	  
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
	  
	  for( j=0; j<LAL_NUM_IFO; j++ ) 
	    {
	      cohInspFilterParams->detIDVec->data[j] = caseID[j];
	    }
	  
	  /* CHECK: Steve says that this can be replaced with ~XLALReadIfo functions*/
	  /* Note that the ifo orders in InterferometerNumber and 
	     lalCachedDetectors are different:
	     caseID[0,..,5]=(G1,H1,H2,L1,V1,T1), whereas
	     lalCachedDetectors[0,...]=(LHO,LLO,V1,G1,T1,CIT,...)*/
	  w=0;
	  if( caseID[0] ) cohInspFilterParams->detectorVec->detector[w++] = lalCachedDetectors[3];
	  if( caseID[1] ) cohInspFilterParams->detectorVec->detector[w++] = lalCachedDetectors[0];
	  if( caseID[2] ) cohInspFilterParams->detectorVec->detector[w++] = lalCachedDetectors[0];
	  if( caseID[3] ) cohInspFilterParams->detectorVec->detector[w++] = lalCachedDetectors[1];
	  if( caseID[4] ) cohInspFilterParams->detectorVec->detector[w++] = lalCachedDetectors[4];
	  if( caseID[5] ) cohInspFilterParams->detectorVec->detector[w++] = lalCachedDetectors[2];
	  
	  
	  /* Get beam pattern coefficients if necessary */
	  if ( (numDetectors == 3 && !( caseID[1] && caseID[2])) || numDetectors == 4 )
	    {
	      cohInspBeamVec = cohInspFilterInput->beamVec;
	      w=0;
	      if( verbose ) fprintf(stdout, "This network requires beam-pattern coefficients - reading them in...\n");
	      for ( j=0; j<LAL_NUM_IFO; j++ ) {
		if ( caseID[j] ) 
		  {
		    /*CHECK:*/
		    LALSnprintf( nameArrayBeam[j], LALNameLength*sizeof(CHAR), "%sBeamCoeff.dat", caseIDChars[j]);
		    filePtr[w] = fopen(nameArrayBeam[j], "r");
		    if(!filePtr[w])
		      {
			fprintf(stderr,"The file %s containing the coefficients could not be found - exiting...\n",nameArrayBeam[j]);
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
	  
	  /* Read in the snippets associated with thisCoinc trigger */
	  l = 0;
	  for( j=0; j<LAL_NUM_IFO; j++ )
	    {
	      if( caseID[j] )
		{
		  if( verbose ) fprintf(stdout, "getting the COMPLEX8TimeSeries\n");
		  LAL_CALL( LALFrOpen(&status,&frStream,NULL,ifoframefile[j]), &status);
		  if(!frStream)
		    {
		      fprintf(stdout,"The file %s does not exist - exiting...\n", ifoframefile[j] );
		      goto cleanexit;
		    }
		  
		  frChan.name = nameArrayCData[j];
		  /*CHECK: Need to replace index l with an index that assigns these
		    CData to the correct ifo*/
		  LAL_CALL( LALFrGetCOMPLEX8TimeSeries( &status, &(cohInspCVec->cData[l]), &frChan, frStream), &status);
		  /*CHECK: Need to worry about WRAPPING of time-slides here */
		  /* tempTime is the start time of cData plus -(time slide)*/
		  tempTime[l] += cohInspCVec->cData[l].epoch.gpsSeconds + cohInspCVec->cData[l].epoch.gpsNanoSeconds * 1e-9;
		  if( verbose ) fprintf(stdout,"tempTime = %f\n",tempTime[l]);
		  LAL_CALL( LALFrClose( &status, &frStream ), &status );
		  
		  /* CHECK: delete this after updating the cohinspfilterparams defn
		     cohInspFilterParams->segNorm[l] = thisCoinc->snglInspiral[k]
		  */
		  cohInspFilterParams->sigmasq[l] = thisCoinc->snglInspiral[j]->sigmasq;
		  l++;
		}/* Closes if( caseID[j] ) */
	    }
	  
	  /* If we can estimate distance then compute templateNorm */
	  /* At present, this is only good for frequency domain tmplts */
	  /* Since each detector's data has been filtered with templates */
	  /* that have the same mass pair, templateNorm is the same for */
	  /* every detector and needs to be computed only once.         */
	  
	  /* CHECK: with the sigmasq being read in now, 
	     the lines below may be redundant */
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
	  if ( (numDetectors == 3 && !( caseID[1] && caseID[2])) || numDetectors == 4 )
	    {
	      while( thisEventTemp )
		{
		  convertParams.system=COORDINATESYSTEM_EQUATORIAL;
		  convertParams.zenith= NULL;
		  convertParams.gpsTime=&(thisEventTemp->end_time);
		  /*The next few lines get the parameters ready for conversion. 
		    Need to check to make sure that thisEvent->ligo_axis_dec is 
		    in proper notation, i.e. */
		  
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
	      
	    }/*end if( thisEvent ) */
	  
	  /* Finalize */
	  LAL_CALL( LALCDestroyVector( &status, &(tempSnippet.data) ), &status );
	  LALFree( cohInspFilterInput->tmplt );
	  cohInspFilterInput->tmplt = NULL;
	  free( cohInspInitParams );
	  
	  /* Destroy params structure for coherent filter code */
	  LAL_CALL( LALCoherentInspiralFilterParamsFinalize (&status,&cohInspFilterParams), &status );
	  
	  /* Destroy input structure for coherent filter code */
	  LAL_CALL( LALCoherentInspiralFilterInputFinalize (&status, &cohInspFilterInput), &status);
	  
	  for(j=0;j<LAL_NUM_IFO;j++)
	    {
	      memcpy(caseIDChars[j], "0",3);
	      caseID[j] = 0;
	      strcpy( nameArrayBeam[j], "0" );
	      strcpy( nameArrayCData[j], "0" );
	    }
	  for(j=0;j<5;j++)
	    {
	      tempTime[j] = 0.0;
	      timeptDiff[j] = 0;
	    }
	  tempTime[5] = 0.0;
	}/* Close loop over thisCoinc*/
    } /* End the condition "if ( numTriggers < 0 )" */
  
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
  
  /* free the frame cache */
  if( frInCache ) LAL_CALL( LALDestroyFrCache( &status, &frInCache ), &status );
  
  if ( frInType )    free( frInType );
  
  if ( verbose ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
  exit(0);
  
}/* main function end */


/* ------------------------------------------------------------------------- */


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
"  --maximize-over-chirp        do clustering\n"\
"  --frame-type TAG             input data is contained in frames of type TAG\n"\
"\n"
#define USAGE3 \
"  --write-events               write events\n"\
"  --write-cohsnr               write cohsnr\n"\
"  --output-path                write files here\n"\
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
     {"write-events",             no_argument,       &eventsOut,         1 },
     {"write-cohsnr",             no_argument,       &cohSNROut,         1 },
     {"output-path",              required_argument, 0,                 'P'},
     {"H1-framefile",             required_argument, 0,                 'A'},
     {"H2-framefile",             required_argument, 0,                 'Z'},
     {"L1-framefile",             required_argument, 0,                 'L'},
     {"V1-framefile",             required_argument, 0,                 'V'},
     {"G1-framefile",             required_argument, 0,                 'G'},
     {"frame-type",               required_argument, 0,                 'S'},
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
           "A:B:G:S:I:L:l:e:g:W:X:Y:t:w:P:T:V:Z:d:f:h:p:r:u:v:",
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
           strcpy(ifoframefile[1],optarg);
           H1file = 1;
           ADD_PROCESS_PARAM( "string", "%s", ifoframefile[1] );
           break;

         case 'G':
           strcpy(ifoframefile[0],optarg);
           G1file = 1;
           ADD_PROCESS_PARAM( "string", "%s", ifoframefile[0] );
           break;

         case 'L':
           strcpy(ifoframefile[3],optarg);
           L1file = 1;
           ADD_PROCESS_PARAM( "string", "%s", ifoframefile[3] );
           break;

         case 'T':
           strcpy(ifoframefile[4],optarg);
           T1file = 1;
           ADD_PROCESS_PARAM( "string", "%s", ifoframefile[4] );
           break;

         case 'V':
           strcpy(ifoframefile[5],optarg);
           V1file = 1; 
           ADD_PROCESS_PARAM( "string", "%s", ifoframefile[5] );
           break;

         case 'Z':
           strcpy(ifoframefile[2],optarg);
           H2file = 1;
           ADD_PROCESS_PARAM( "string", "%s", ifoframefile[2] );
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
