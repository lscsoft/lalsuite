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

double rint(double x);
int arg_parse_check( int argc, char *argv[], MetadataTable procparams );


/*
 *
 * variables that control program behaviour
 *
 */


/* input data parameters */

CHAR  *chanNumber       = NULL;         /* name of data channel         */
INT4  numPoints         = -1;           /* points in a segment          */
INT4  numSegments       = 1;            /* number of segments           */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
REAL4  fLow             = -1;           /* low frequency cutoff         */

/* matched filter parameters */

CHAR *bankFileName      = NULL;         /* name of input template bank  */
INT4  startTemplate     = 0;            /* index of first template      */
INT4  stopTemplate      = 1;            /* index of last template       */
INT4  eventCluster      = -1;           /* perform chirplen clustering  */

/*Coherent code specific inputs*/

char H1filename[256];
char Lfilename[256];
char GEOfilename[256];
char VIRGOfilename[256];
char TAMAfilename[256];
char H2filename[256];

/*hard code the beam file names*/
char H1Beam[256];
char H2Beam[256];
char L1Beam[256];
char GEOBeam[256];
char VIRGOBeam[256];
char TAMABeam[256];

UINT4 siteID[6] = {0,1,2,3,4,0}; /*H1 L1 V G T H2*/
UINT2 caseID[6] = {0,0,0,0,0,0};

static BOOLEAN verbose = 0;

static BOOLEAN   cohSNROut            = 0;
CHAR             *cohSNROutFrame      = NULL;
CHAR             outputPath[FILENAME_MAX];
INT4             numTmplts            = -1;
UINT4            numBeamPoints        = 100;
REAL4            cohSNRThresh         = -1;
UINT4            maximiseOverChirp    = 0; 

FILE *fp[4];


int main( int argc, char *argv[] )
{
  /* lal function variables */
  static LALStatus             status;

  /* frame input data */
  FrStream     *frStream = NULL;
  FrChanIn      frChan;

  /* frame output data */

  struct FrFile *frOutFile  = NULL;
  struct FrameH *outFrame   = NULL;

/* inspiral events */
  INT4                          numEvents   = 0;
  MultiInspiralTable           *event       = NULL;
  MultiInspiralTable           *eventList   = NULL;
  MetadataTable                 savedEvents;

  /* output */

  MetadataTable procparams;
  CHAR fname[FILENAME_MAX];

 /* counters and other variables */
  INT4                          i,j,k,l;
  UINT4                         numDetectors=0;
  REAL4                         theta,phi,vPlus,vMinus;

  UINT2Vector                  *detIDVec = NULL;
  DetectorVector               *detectorVec = NULL;

  CoherentInspiralInitParams   *cohInspInitParams = NULL;
  CoherentInspiralFilterParams *cohInspFilterParams = NULL;
  CoherentInspiralFilterInput  *cohInspFilterInput = NULL;
  CoherentInspiralBeamVector   *cohInspBeamVec = NULL;
  CoherentInspiralCVector      *cohInspCVec = NULL;
  CoherentInspiralEvent        *cohInspEvent = NULL;

  InspiralTemplate             *bankHead = NULL;

  char namearray[6][256]  = {"0","0","0","0","0","0"}; /* input data frame files */
  char namearray2[6][256] = {"0","0","0","0","0","0"}; /* beam files */
  char namearray3[6][256] = {"0","0","0","0","0","0"}; /* chan names */

  set_debug_level( "1" ); /* change with parse option */
  memset( outputPath, 0, FILENAME_MAX * sizeof(CHAR) );

if (verbose)  fprintf(stdout, "calling parse options\n");
  arg_parse_check( argc, argv, procparams );
if (verbose)  fprintf(stdout, "parse options called\n");


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
    /* if there are no tmplts, store the time we would have analyzed and exit */
    fprintf( stdout, "no templates found in template bank file: %s\n"
        "exiting without searching for events.\n", bankFileName ); 

    goto cleanexit;
  }
  else if ( numTmplts != 1 )
  {
    fprintf( stderr, "must only have 1 template in the bank file %s\n", 
        bankFileName );
    goto cleanexit;
  }

  if ( verbose ) fprintf( stdout, "parsed %d templates from %s\n", 
      numTmplts, bankFileName );

/* read in the network detectors */
  for (l=0;l<6;l++)
    {
      if(caseID[l] == 1)
	numDetectors++;
    }

if (verbose)  fprintf(stdout, "You have specified %2d  detector(s).\n",numDetectors);
if (verbose)  fprintf(stdout, "The caseID is: %d %d %d %d %d %d (H1,L1,VIRGO,GEO,TAMA,H2) \n",caseID[0],caseID[1],caseID[2],caseID[3],caseID[4],caseID[5]);

  if (numDetectors > 4)
    {
    fprintf(stdout, "Too many detectors specified - exiting. \n");
    goto cleanexit;
    }

  if (numDetectors == 0)
    {
      fprintf(stdout, "You must specify data filename(s) for 1 to 4 detectors - exiting. \n");
      goto cleanexit;
    }


  /*
   *
   * read in the input data channels
   *
   */


 /* 
   *
   * create and populate coherentInspiral initialization structure 
   *
   */
  cohInspInitParams = (CoherentInspiralInitParams *) LALMalloc (sizeof(CoherentInspiralInitParams));

  if ( ! ( cohInspInitParams = (CoherentInspiralInitParams *) 
	   LALCalloc (1,sizeof(CoherentInspiralInitParams)) ))
  {
    fprintf( stderr, "could not allocate memory for coherentInspiral init params\n" );
    goto cleanexit;
  }
  
  cohInspInitParams->numDetectors            = numDetectors;
  cohInspInitParams->numSegments             = numSegments;
  cohInspInitParams->numPoints               = numPoints;
  cohInspInitParams->numBeamPoints           = numBeamPoints;
  cohInspInitParams->cohSNROut               = cohSNROut;


  /*
   *
   * create the data structures needed for coherentInspiral
   *
   */

  if ( verbose ) fprintf( stdout, "initializing coherentInspiral...\n " );

  /* initialize coherentInspiral filter functions */
  LAL_CALL( LALCoherentInspiralFilterInputInit (&status, &cohInspFilterInput,
						cohInspInitParams), &status );

  TestStatus (&status, "0", 1);
  ClearStatus (&status);

   /*
   * information for calculating chirp time
   */

 
  cohInspFilterInput->tmplt = (InspiralTemplate *)
    LALMalloc (sizeof(InspiralTemplate) );
  memset( cohInspFilterInput->tmplt, 0, sizeof(InspiralTemplate) );
  cohInspFilterInput->tmplt->mass1 = bankHead->mass1;
  cohInspFilterInput->tmplt->mass2 = bankHead->mass2;
  REAL8 m1;
  REAL8 m2;
  m1 = cohInspFilterInput->tmplt->mass1;
  m2 = cohInspFilterInput->tmplt->mass2;
  cohInspFilterInput->tmplt->totalMass = m1 + m2;
  cohInspFilterInput->tmplt->mu = m1 * m2 / (m1 + m2);
  cohInspFilterInput->tmplt->eta = (m1 * m2) / ((m1 + m2) * (m1 + m2 ));

if (verbose)  fprintf( stdout, "m1:%f m2:%f totalmass:%f mu:%f eta%f\n", m1, m2,cohInspFilterInput->tmplt->totalMass,cohInspFilterInput->tmplt->mu,cohInspFilterInput->tmplt->eta);
  

  /* fill the params structure */
  LAL_CALL( LALCoherentInspiralFilterParamsInit (&status, &cohInspFilterParams,
						 cohInspInitParams),&status );
  TestStatus (&status, "0", 1);
  ClearStatus (&status);
  
  cohInspFilterParams->numDetectors            = numDetectors;
  cohInspFilterParams->numSegments             = numSegments;
  cohInspFilterParams->numPoints               = numPoints;
  cohInspFilterParams->numBeamPoints           = numBeamPoints;
  cohInspFilterParams->deltaT                  = 1/((REAL4) sampleRate);
  cohInspFilterParams->cohSNROut               = cohSNROut;
  cohInspFilterParams->cohSNRThresh            = cohSNRThresh;
  cohInspFilterParams->numTmplts               = numTmplts;
  cohInspFilterParams->fLow                    = fLow;
  cohInspFilterParams->maximiseOverChirp       = maximiseOverChirp;
  
  detIDVec = cohInspFilterParams->detIDVec;


/*assign detIDs to the coincident detectors in the network */
  for ( i=0 ; i < 6 ; i++) {
    detIDVec->data[i] = caseID[i];
  }
  
 /* create and fill the DetectorVector structure of detector IDs*/
  
  detectorVec = cohInspFilterParams->detectorVec;
  
  i=0;
  for ( j=0 ; j < 6 ; j++ ) {
    if ( caseID[j] ) { 
      detectorVec->detector[i++] = lalCachedDetectors[j];
    }
  }

  if (caseID[5]) {
    detectorVec->detector[numDetectors-1] = lalCachedDetectors[0];
  }



  /* Note that some of the namearray3 elements will need to be changed */
  /* depending what channels are in the corresponding frame files */
  if(caseID[0])
    {
    strcpy(namearray[0],H1filename);
    strcpy(namearray2[0],"HBeam.dat");
    LALSnprintf( namearray3[0], LALNameLength*sizeof(CHAR), "H1:LSC-AS_Q_CData_%s", chanNumber );
    }

  if(caseID[1])
    {
    strcpy(namearray[1],Lfilename);
    strcpy(namearray2[1],"LBeam.dat");
    LALSnprintf( namearray3[1], LALNameLength*sizeof(CHAR), "L1:LSC-AS_Q_CData_%s", chanNumber );
    }

  if(caseID[2])
    {
      /* for the moment, give virgo H2 data and beam pattern functions since */
      /* virgo data is not available */
    strcpy(namearray[2],VIRGOfilename);
    /*strcpy(namearray2[2],"VIRGOBeam.dat");*/ 
    strcpy(namearray2[2],"HBeam.dat");
    LALSnprintf(namearray3[2], LALNameLength*sizeof(CHAR), "H2:LSC-AS_Q_CData_%s", chanNumber );
    }

  if(caseID[3])
    {
    strcpy(namearray[3],GEOfilename);
    strcpy(namearray2[3],"GEOBeam.dat");
    LALSnprintf(namearray3[3], LALNameLength*sizeof(CHAR), "G1:LSC-AS_Q_CData_%s", chanNumber );
    }

  if(caseID[4])
    {
    strcpy(namearray[4],TAMAfilename);
    strcpy(namearray2[4],"TAMABeam.dat");
    LALSnprintf(namearray3[4], LALNameLength*sizeof(CHAR), "T1:LSC-AS_Q_CData_%s", chanNumber );
    }

  if (caseID[5])
    {
    strcpy(namearray[5],H2filename);
    strcpy(namearray2[5],"HBeam.dat");
    LALSnprintf(namearray3[5], LALNameLength*sizeof(CHAR), "H2:LSC-AS_Q_CData_%s", chanNumber );
    }

if (verbose)  fprintf(stdout,"built namearrays\n");

  /* Now read in the beam coefficients if necessary */

  cohInspBeamVec = cohInspFilterInput->beamVec;

 if ( numDetectors > 2 && !( caseID[0] && caseID[5]) ) {
  l=0;
  for ( j=0 ; j < 6 ; j++ ) {
    if ( caseID[j] ) { 
      fp[l] = fopen(namearray2[j], "r");
      if(!fp[l])
	{
	  fprintf(stdout,"The file %s containing the coefficients could not be found - exiting...\n",namearray2[j]);
	  goto cleanexit;
	}
      for ( k=0 ; k< (INT4) numBeamPoints ; k++)
	{
	  fscanf(fp[l],"%f %f %f %f",&theta,&phi,&vPlus,&vMinus);
	  cohInspBeamVec->detBeamArray[l].thetaPhiVs[k].data->data[0] = theta;
	  cohInspBeamVec->detBeamArray[l].thetaPhiVs[k].data->data[1] = phi;
	  cohInspBeamVec->detBeamArray[l].thetaPhiVs[k].data->data[2] = vPlus;
	  cohInspBeamVec->detBeamArray[l].thetaPhiVs[k].data->data[3] = vMinus;
	}
      fclose(fp[l++]);
    }
  } 
 }

 cohInspCVec = cohInspFilterInput->multiCData;

 TestStatus (&status, "0", 1);
 ClearStatus (&status);

 if (verbose) {
   fprintf(stdout,"reading data from frames\n");
   fprintf(stdout,"namearray: %s\n %s\n %s\n %s\n %s\n %s\n",namearray[0],namearray[1],namearray[2],namearray[3],namearray[4],namearray[5]);
   fprintf(stdout,"namearray3: %s\n %s\n %s\n %s\n %s\n %s\n",namearray3[0],namearray3[1],namearray3[2],namearray3[3],namearray3[4],namearray3[5]);
 }

   l=0;
  for ( j=0 ; j < 6 ; j++ ) {
    if ( caseID[j] ) {
      if (verbose) fprintf(stdout,"opening C framefile %s\n",namearray[j]);
      LAL_CALL( LALFrOpen(&status,&frStream,0,namearray[j]), &status);
      if(!frStream)
	{
	  fprintf(stdout,"The file %s does not exist - exiting...\n",
		  namearray[j]);
	  goto cleanexit;
	}
      if (verbose) fprintf(stdout,"getting the complex8timeseries\n");
      frChan.name = namearray3[j];
      LAL_CALL( LALFrGetCOMPLEX8TimeSeries( &status, &(cohInspCVec->cData[l]), &frChan, frStream), &status);
      if (verbose) fprintf(stdout,"closing frame stream\n");
      LAL_CALL( LALFrClose( &status, &frStream), &status);
      l++;
      TestStatus (&status, "0", 1);
      ClearStatus (&status);

    }
  }

  /*Do the filtering and output events */
  cohInspEvent = NULL;
  if (verbose) fprintf(stdout,"filtering the data\n");
  LALCoherentInspiralFilterSegment (&status, &cohInspEvent, cohInspFilterInput, cohInspFilterParams); 

  if ( cohSNROut )
    {
      strcpy( cohInspFilterParams->cohSNRVec->name, "coherent");
      outFrame = fr_add_proc_REAL4TimeSeries( outFrame,
	  cohInspFilterParams->cohSNRVec, "none", "SNR" );
    }
  if (outputPath[0] )
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s/%s", outputPath, cohSNROutFrame);
    }
  else
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s", cohSNROutFrame );
    }
  if ( verbose ) fprintf( stdout, "writing frame data to %s....", fname );
  frOutFile = FrFileONew( fname, 0);
  FrameWrite( outFrame, frOutFile);
  FrFileOEnd( frOutFile );
  if ( verbose ) fprintf(stdout, "done\n");



 cleanexit:

  /* Destroy params structure for coherent filter code */
  LALCoherentInspiralFilterParamsFinalize (&status, &cohInspFilterParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  /* Destroy input structure for coherent filter code */
  LALCoherentInspiralFilterInputFinalize (&status, &cohInspFilterInput);
    TestStatus (&status, "0", 1);
    ClearStatus (&status);

 /* free the rest of the memory, check for memory leaks and exit */
  if ( bankFileName ) free( bankFileName );
  if ( chanNumber ) free( chanNumber );

  if ( verbose ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();

  return 0;

}


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
void
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

#define USAGE \
"lalapps_inspiral [options]\n\n"\
"  --help                       display this message\n"\
"  --verbose                    print progress information\n"\
"  --version                    print version information and exit\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"  --low-frequency-cutoff F     low f cutoff of previously filtered data\n"\
"\n"\
"  --bank-file FILE             read template bank parameters from FILE\n"\
"  --channel-number             the channel (segment) number\n"\
"  --sample-rate F              filter data at F Hz, downsampling if necessary\n"\
"  --segment-length N           set data segment length to N points\n"\
"  --cohsnr-threshold RHO          set signal-to-noise threshold to RHO\n"\
"  --maximize-over-chirp        do clustering\n"\
"\n"\
"  --write-cohsnr <filename(.gwf)>    write cohsnr to specified frame file\n"\
"  --output-path                 the path to where the user wants to the output file written\n"\
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
     {"verbose",                  no_argument,       0,                 'w'},
     {"help",                     no_argument,       0,                 'h'},
     {"version",                  no_argument,       0,                 'v'},
     {"debug-level",              no_argument,       0,                 'd'},
     {"low-frequency-cutoff",     required_argument, 0,                 'f'},
     {"bank-file",                required_argument, 0,                 'u'},
     {"sample-rate",              required_argument, 0,                 'r'},
     {"segment-length",           required_argument, 0,                 'n'},
     {"cohsnr-threshold",         required_argument, 0,                 'p'},
     {"maximize-over-chirp",      required_argument, 0,                 'c'},
     {"write-cohsnr",             required_argument, 0,                 'o'},
     {"output-path",              required_argument, 0,                 'P'},
     {"H1-framefile",             no_argument,       0,                 'A'},
     {"H2-framefile",             no_argument,       0,                 'Z'},
     {"L-framefile",              no_argument,       0,                 'L'},
     {"V-framefile",              no_argument,       0,                 'V'},
     {"G-framefile",              no_argument,       0,                 'G'},
     {"T-framefile",              no_argument,       0,                 'T'},
     {"channel-number",           required_argument, 0,                 'k'},
     {0, 0, 0, 0,}
   };

   int c;

   while (1)
     {
       /* getopt_long stores long options here */
       int option_index = 0;
       size_t optarg_len;

       c = getopt_long_only( argc, argv,
	   "A:G:L:P:T:V:Z:c:d:f:h:k:n:o:p:r:u:v:w:",
	   long_options, &option_index );

       if ( c == -1 )
	 {
	   break;
	 }

       switch ( c )
	 {
	 case 'A':
	   strcpy(H1filename,optarg);
	   caseID[0] = 1;
	   break;

	 case 'G':
	   strcpy(GEOfilename,optarg);
	   caseID[3] = 1;  
	   break;

	 case 'L':
	   strcpy(Lfilename,optarg);
	   caseID[1] = 1;
	   break;

	 case 'T':
	   strcpy(TAMAfilename,optarg);
	   caseID[4] = 1;      
	   break;

	 case 'V':
	   strcpy(VIRGOfilename,optarg);
	   caseID[2] = 1;	   
	   break;

	 case 'Z':
	   strcpy(H2filename,optarg);
	   caseID[5] = 1;
	   break;

	 case 'P':
	   memset( outputPath, 0, FILENAME_MAX * sizeof(CHAR) );
	   LALSnprintf( outputPath, FILENAME_MAX * sizeof(CHAR),
			"%s", optarg );
	   break;

	 case 'c': /* set maximiseOverChirp */
	   maximiseOverChirp  = 1;
	   break;

	 case 'd': /* set debuglevel */
	   lalDebugLevel = atoi (optarg);
	   break;

	 case 'f': /* set fLow */
	   fLow = (REAL4) atof (optarg);
	   break;

	 case 'h':
	   fprintf( stdout, USAGE );
	   exit( 0 );
	   break;

	 case 'n': /* set number of points in a segment */
	   numPoints = atoi (optarg);
	   break;

	 case 'o': /* sets flag to write coherent SNR */
	   cohSNROut = 1;
	   optarg_len = strlen( optarg ) +1;
	   cohSNROutFrame = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	   memcpy( cohSNROutFrame, optarg, optarg_len );
	   break;

	 case 'p': /* set coherent SNR threshold */
	   cohSNRThresh = atof (optarg);
	   break;

	 case 'r': /* set sampling rate */
	   sampleRate = atoi (optarg);
	   break;

	 case 'u':
           /* create storage for the calibration frame cache name */
	   optarg_len = strlen( optarg ) + 1;
	   bankFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	   memcpy( bankFileName, optarg, optarg_len );
	   /* ADD_PROCESS_PARAM( "string", "%s", optarg );*/
	   break;

	 case 'v':
         /* print version information and exit */
           fprintf( stdout, "LIGO/LSC Multi-Detecter Search Code\n" 
                 "Bose/Seader <sukanta@wsu.edu> <sseader@wsu.edu>\n"
                 "CVS Version: " CVS_ID_STRING "\n"
		 "CVS Tag: " CVS_NAME_STRING "\n" );
	   exit( 0 );
	   break;

	 case 'w': /* set verbosity */
	   verbose = 1;
	   break;

	 case 'k':
	   optarg_len = strlen( optarg ) +1;
	   chanNumber = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
	   memcpy( chanNumber, optarg, optarg_len );
	   break;

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
