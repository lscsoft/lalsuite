/*----------------------------------------------------------------------- 
 * 
 * File Name: coherent_inspiral.c
 *
 * Author: Bose, S. and Seader, S.
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

INT4   sampleRate           = -1;            /* sample rate of filter data   */
REAL4  fLow                 = -1;            /* low frequency cutoff         */

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
                                
CHAR   bankFileName[FILENAME_MAX]; /* name of input template bank  */
UINT4  cohSNROut            = 0;    /* default is not to write frame */
UINT4  eventsOut            = 0;    /* default is not to write events */
REAL4  cohSNRThresh         = -1;
INT4   maximizeOverChirp    = 0;    /* default is no clustering */
INT4   verbose              = 0;
CHAR   outputPath[FILENAME_MAX];

LALStatus             status;
LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

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

  CHAR   framename[FILENAME_MAX];
  CHAR   xmlname[FILENAME_MAX];
  CHAR   cohdataStr[LALNameLength];
  CHAR   caseIDChars[4][LIGOMETA_IFOS_MAX] = {"0","0","0","0"};
  CHAR tempStr[FILENAME_MAX];

  INT4   numTmplts        = 0; /* number of templates */
  INT4   cohSegLength     = 4; /* This should match hardcoded value in inspiral.c */
  INT4   numPoints        = 0;
  INT4   startTemplate    = -1;           
  INT4   stopTemplate     = -1;
  INT4   bankFileNameLength = 0;
  UINT4  numSegments      = 1;       /* number of segments */
  UINT4  numBeamPoints    = 0;       /* number of sky position templates */
  UINT4  nCohDataFr       = 0;
  UINT8  eventID          = 0;

  REAL4  m1               = 0;
  REAL4  m2               = 0;

  /* counters and other variables */
  INT4   i,j,k,l,w;
  REAL4  theta,phi,vPlus,vMinus;
  UINT4  numDetectors     = 0;
  REAL8  tempTime[6]      = {0.0,0.0,0.0,0.0,0.0,0.0};
  INT4   timeptDiff[5]    = {0,0,0,0,0};
  UINT2  caseID[6]        = {0,0,0,0,0,0}; /* H1 L V G T H2 */
  INT4   h1ChanNum        = -1;
  INT4   h2ChanNum        = -1;
  INT4   lChanNum         = -1;
  INT4   geoChanNum       = -1;
  INT4   virgoChanNum     = -1;
  INT4   tamaChanNum      = -1;
  INT4   chanNumArray[6]  = {-1, -1, -1, -1, -1, -1};

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

  char namearray[6][256]  = {"0","0","0","0","0","0"}; /* input frame files */
  char namearray2[6][256] = {"0","0","0","0","0","0"}; /* beam files */
  char namearray3[6][256] = {"0","0","0","0","0","0"}; /* chan names */

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
  
  if ( verbose ) fprintf( stdout, "parsed %d templates from %s\n", 
      numTmplts, bankFileName );


  /* Loop over templates (or event id's) */
  
  numPoints = sampleRate * cohSegLength;
  bankTemp = bankHead;
  savedEvents.multiInspiralTable = NULL;
  k = 0;   

  for( i=0; i<numTmplts; i++)
    {
      memcpy( caseIDChars[k], &bankTemp->ifo, sizeof(caseIDChars[k] - 1) );
      if( verbose ) fprintf(stdout,"caseIDChars = %s %s %s %s\n",caseIDChars[0],caseIDChars[1],caseIDChars[2],caseIDChars[3] );
      eventID = bankTemp->event_id->id;
      if( verbose ) fprintf(stdout,"eventID = %Ld\n",eventID );
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
		      h1ChanNum++;
		    }
		  else if( !strcmp( caseIDChars[j], "L1" ) )
		    {
		      caseID[1] = 1;
		      lChanNum++;
		    }
		  else if( !strcmp( caseIDChars[j], "V1" ) )
		    {
		      caseID[2] = 1;
		      virgoChanNum++;
		    }
		  else if( !strcmp( caseIDChars[j], "G1" ) )
		    {
	              caseID[3] = 1;
		      geoChanNum++;
		    }
		  else if( !strcmp( caseIDChars[j], "T1" ) )
		    {
		      caseID[4] = 1;
		      tamaChanNum++;
		    }
		  else if( !strcmp( caseIDChars[j], "H2" ) )
		    {
		      caseID[5] = 1;
		      h2ChanNum++;
		    }
		}

	      chanNumArray[0] = h1ChanNum;
	      chanNumArray[1] = lChanNum;
	      chanNumArray[2] = virgoChanNum;
	      chanNumArray[3] = geoChanNum;
	      chanNumArray[4] = tamaChanNum;
	      chanNumArray[5] = h2ChanNum;

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

	      /* Now check that the number of detectors matches the number of cache files provided */

	      l=0;
	      if( H1file ) l++;
	      if( H2file ) l++;
	      if( Lfile ) l++;
	      if( GEOfile ) l++;
	      if( TAMAfile ) l++;
	      if( VIRGOfile ) l++;
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
	      if( verbose ) fprintf(stdout,"caseID = %d %d %d %d %d %d (H1,L1,V1,G1,T1,H2)\n", caseID[0], caseID[1], caseID[2], caseID[3], caseID[4], caseID[5]);

	      if(caseID[0])
		{
		  strcpy(namearray[0],H1filename);
		  strcpy(namearray2[0],"HBeam.dat");
		  LALSnprintf( namearray3[0], LALNameLength*sizeof(CHAR), "H1:LSC-AS_Q_CData_%d", h1ChanNum );
		}

              if(caseID[1])
                {
                  strcpy(namearray[1],Lfilename);
                  strcpy(namearray2[1],"LBeam.dat");
                  LALSnprintf( namearray3[1], LALNameLength*sizeof(CHAR), "L1:LSC-AS_Q_CData_%d", lChanNum );
                }

              if(caseID[2])
                {
		  /* for the moment, give virgo H2 data and beam pattern functions since */
		  /* virgo data is not available */
		  strcpy(namearray[2],VIRGOfilename);
		  /*strcpy(namearray2[2],"VIRGOBeam.dat");*/ 
		  strcpy(namearray2[2],"VIRGOBeam.dat");
		  LALSnprintf(namearray3[2], LALNameLength*sizeof(CHAR), "H2:LSC-AS_Q_CData_%d", virgoChanNum );
		}

	      if(caseID[3])
		{
		  strcpy(namearray[3],GEOfilename);
		  strcpy(namearray2[3],"GEOBeam.dat");
		  LALSnprintf(namearray3[3], LALNameLength*sizeof(CHAR), "G1:LSC-AS_Q_CData_%d", geoChanNum );
		}

              if(caseID[4])
                { 
                  strcpy(namearray[4],TAMAfilename);
                  strcpy(namearray2[4],"TAMABeam.dat");
                  LALSnprintf(namearray3[4], LALNameLength*sizeof(CHAR), "T1:LSC-AS_Q_CData_%d", tamaChanNum );
                }

              if(caseID[5])
                {
                  strcpy(namearray[5],H2filename);
                  strcpy(namearray2[5],"HBeam.dat");
                  LALSnprintf(namearray3[5], LALNameLength*sizeof(CHAR), "H2:LSC-AS_Q_CData_%d", h2ChanNum );
                }
		  
              if(verbose)  fprintf(stdout,"built namearrays\n");

              if(verbose) {
                  fprintf(stdout,"cache files: %s\n %s\n %s\n %s\n %s\n %s\n",namearray[0],namearray[1],namearray[2],namearray[3],namearray[4],namearray[5]);
                  fprintf(stdout,"channels: %s\n %s\n %s\n %s\n %s\n %s\n",namearray3[0],namearray3[1],namearray3[2],namearray3[3],namearray3[4],namearray3[5]);
	      }

	      /* Now get the data for the involved detectors */

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

              /* get beam pattern coefficients if necessary */
	      if ( (numDetectors == 3 && !( caseID[0] && caseID[5])) || numDetectors == 4 )
		{
		  cohInspBeamVec = cohInspFilterInput->beamVec;
		  w=0;
		  if( verbose ) fprintf(stdout, "This network requires beam-pattern coefficients - reading them in...\n");
		  for ( j=0; j<6; j++ ) {
		    if ( caseID[j] ) { 
		      filePtr[w] = fopen(namearray2[j], "r");
		      if(!filePtr[w])
			{
			  fprintf(stderr,"The file %s containing the coefficients could not be found - exiting...\n",namearray2[j]);
			  goto cleanexit;
			}
		      for ( l=0 ; l < (INT4) numBeamPoints ; l++)
			{ 
			  if( verbose ) fprintf(stdout,"scanning a beam file...");
			  fscanf(filePtr[w],"%f %f %f %f",&theta,&phi,&vPlus,&vMinus);
			  if( verbose ) fprintf(stdout,"done\n");
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
		      /* Advance the stream to the appropriate frame */
		      for( w=0; w<chanNumArray[j]; w++)
			{
			  LAL_CALL( LALFrNext( &status,frStream ), &status );
			}
		      frChan.name = namearray3[j];
		      LAL_CALL( LALFrGetCOMPLEX8TimeSeries( &status, &(cohInspCVec->cData[l]), &frChan, frStream), &status);
		      tempTime[l] = cohInspCVec->cData[l].epoch.gpsSeconds + cohInspCVec->cData[l].epoch.gpsNanoSeconds * 1e-9;
		      if( verbose ) fprintf(stdout,"tempTime = %f\n",tempTime[l]);
		      LAL_CALL( LALFrClose( &status, &frStream ), &status );
		      l++;
		    }
		}

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
		}

	      /* Now that the time series are commensurate, do the filtering... */
                
	      LALCoherentInspiralFilterSegment (&status, &thisEvent, cohInspFilterInput, cohInspFilterParams);
	      cohInspFilterParams->cohSNRVec->epoch = cohInspCVec->cData[0].epoch;

	      if ( cohSNROut )
		{
		  LALSnprintf( cohdataStr, LALNameLength*sizeof(CHAR),
			       "SNR_%d", nCohDataFr++ );
		  strcpy( cohInspFilterParams->cohSNRVec->name, "Coherent");
		  outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
			       cohInspFilterParams->cohSNRVec, "none", cohdataStr );
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

		}

	      /* initialize */
	      LAL_CALL( LALCDestroyVector( &status, &(tempSnippet.data) ), &status );
	      LALFree( cohInspFilterInput->tmplt );
	      cohInspFilterInput->tmplt = NULL;
	      free( cohInspInitParams );

	      /* Destroy params structure for coherent filter code */
	      LAL_CALL( LALCoherentInspiralFilterParamsFinalize (&status, 
								   &cohInspFilterParams), &status );
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
		  chanNumArray[j] = -1; 
		  strcpy( namearray[j], "0" );
		  strcpy( namearray2[j], "0" );
		  strcpy( namearray3[j], "0" );
		}
	      for(j=0;j<5;j++)
	        {
		  tempTime[j] = 0.0;
		  timeptDiff[j] = 0;
		}
	      tempTime[5] = 0.0;
	  
	    }
	  else
	    {
	      fprintf(stderr,"Error - there is only one detector associated with event %Ld, aborting...\n", eventID );
	      exit( 1 );
	    }
	}
      k++;    
    }

  bankFileNameLength = strlen(bankFileName);
  memcpy(tempStr, bankFileName, bankFileNameLength - 4);


  if( cohSNROut )
    {
      if ( outputPath[0] )
	{
	  LALSnprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s/%s-COHERENT.gwf", outputPath, tempStr);
	}
      else 
	{
	  LALSnprintf( framename, FILENAME_MAX * sizeof(CHAR), "%s-COHERENT.gwf", tempStr );
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
	  LALSnprintf( xmlname, FILENAME_MAX * sizeof(CHAR), "%s/%s-COHERENT.xml", outputPath, tempStr);
	}
      else 
	{
	  LALSnprintf( xmlname, FILENAME_MAX * sizeof(CHAR), "%s-COHERENT.xml", tempStr );
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
	      LALFree( tempEvent );
	      tempEvent = NULL;
	    }
      /* close the output xml file */
      LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );
      if ( verbose ) fprintf( stdout, "done. XML file closed\n" );

    }

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
"\n"
#define USAGE2 \
"  --bank-file FILE             read template bank parameters from FILE\n"\
"  --sample-rate N              set data sample rate to N\n"\
"  --cohsnr-threshold RHO          set signal-to-noise threshold to RHO\n"\
"  --maximize-over-chirp        do clustering\n"\
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
     {"low-frequency-cutoff",     required_argument, 0,                 'f'},
     {"bank-file",                required_argument, 0,                 'u'},
     {"sample-rate",              required_argument, 0,                 'r'},
     {"cohsnr-threshold",         required_argument, 0,                 'p'},
     {"maximize-over-chirp",      no_argument,       &maximizeOverChirp, 1 },
     {"write-events",             no_argument,       &eventsOut,         1 },
     {"write-cohsnr",             no_argument,       &cohSNROut,         1 },
     {"output-path",              required_argument, 0,                 'P'},
     {"H1-framefile",             required_argument, 0,                 'A'},
     {"H2-framefile",             required_argument, 0,                 'Z'},
     {"L-framefile",              required_argument, 0,                 'L'},
     {"V-framefile",              required_argument, 0,                 'V'},
     {"G-framefile",              required_argument, 0,                 'G'},
     {"T-framefile",              required_argument, 0,                 'T'},
     {0, 0, 0, 0}
   };

   int c;
   ProcessParamsTable *this_proc_param = procparams.processParamsTable;

   while (1)
     {
       /* getopt_long stores long options here */
       int option_index = 0;

       c = getopt_long_only( argc, argv,
	   "A:G:L:P:T:V:Z:d:f:h:p:r:u:v:",
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

	 case 'P':
	   memset( outputPath, 0, FILENAME_MAX * sizeof(CHAR) );
	   LALSnprintf( outputPath, FILENAME_MAX * sizeof(CHAR),"%s", optarg );
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
	   ADD_PROCESS_PARAM( "string", "%s", bankFileName );
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
  
   return 0;
}

#undef ADD_PROCESS_PARAM


