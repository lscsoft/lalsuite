/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes  
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

/**
 * \file DriveHough_v3.c
 * \author Badri Krishnan, Alicia Sintes 
 * \brief Driver code for performing Hough transform search on non-demodulated
   data using SFTs from possible multiple IFOs

   Revision: $Id$
 
   History:   Created by Sintes and Krishnan July 04, 2003
              Modifications for S4 January 2006

   \par Description
   
   This is the main driver for the Hough transform routines. It takes as input 
   a set of SFTs from possibly more than one IFO and outputs the number counts 
   using the Hough transform.  For a single IFO, this should be essentially equivalent 
   to DriveHough_v3.  

   \par User input

   The user inputs the following parameters:

   - Search frequency range

   - A file containing list of skypatches to search over.  For each skypatch, 
      the information is:
      - RA and dec of skypatch center.
      - Size in RA and dec.

   - Location of Directory containing the SFTs (must be v2 SFTs).

   - Interferometer (optional)

   - List of linefiles containing information about known spectral disturbances

   - Location of output directory and basename of output files.

   - Block size of running median for estimating psd.

   - The parameter nfSizeCylinder which determines the range of spindown parameters
      to be searched over.

   - Boolean variable for deciding if the SFTs should be inverse noise weighed.

   - Boolean variable for deciding whether amplitude modulation weights should be used.

   - Boolean variables for deciding whether the Hough maps, the statistics, list of 
      events above a threshold, and logfile should be written

   /par Output

   The output is written in several sub-directories of the specified output directory.  The
   first two items are default while the rest are written according to the user input:

   - A directory called logfiles records the user input, contents of the skypatch file 
      and cvs tags contained in the executable (if the user has required logging)

   - A directory called nstarfiles containing the loudest event for each search frequency 
      bin maximised over all sky-locations and spindown parameters.  An event is said to be
      the loudest if it has the maximum significance defined as: (number count - mean)/sigma.

   - A directory for each skypatch containing the number count statistics, the histogram, 
      the list of events, and the Hough maps 
*/


#include "./DriveHoughColor.h"


RCSID( "$Id$");



/* globals, constants and defaults */


extern int lalDebugLevel;

/* boolean global variables for controlling output */
BOOLEAN uvar_printEvents, uvar_printTemplates, uvar_printMaps, uvar_printStats, uvar_printSigma;

/* #define EARTHEPHEMERIS "./earth05-09.dat" */
/* #define SUNEPHEMERIS "./sun05-09.dat"    */

#define EARTHEPHEMERIS "./earth00-04.dat"
#define SUNEPHEMERIS "./sun00-04.dat"

#define MAXFILENAMELENGTH 512 /* maximum # of characters  of a filename */

#define DIROUT "./outMulti"   /* output directory */
#define BASENAMEOUT "HM"    /* prefix file output */

#define THRESHOLD 1.6 /* thresold for peak selection, with respect to the
                              the averaged power in the search band */
#define FALSEALARM 1.0e-9 /* Hough false alarm for candidate selection */
#define SKYFILE "./sky1"      
#define F0 505.0   /*  frequency to build the LUT and start search */
#define FBAND 0.05   /* search frequency band  (in Hz) */
#define NFSIZE  21   /* n-freq. span of the cylinder, to account for spin-down search */
#define BLOCKSRNGMED 101 /* Running median window size */

#define TRUE (1==1)
#define FALSE (1==0)

/* local function prototype */
void PrintLogFile (LALStatus *status, CHAR *dir, CHAR *basename, CHAR *skyfile, LALStringVector *linefiles, CHAR *executable );



/******************************************/

int main(int argc, char *argv[]){

  /* LALStatus pointer */
  static LALStatus  status;  

  /* LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;
  
  /* standard pulsar sft types */ 
  MultiSFTVector *inputSFTs = NULL;
  UINT4 binsSFT;
  
  /* information about all the ifos */
  MultiDetectorStateSeries *mdetStates = NULL;
  UINT4 numifo;

  /* vector of weights */
  REAL8Vector weightsV, weightsNoise;


  /* miscellaneous */
  INT4   houghThreshold, iHmap, nSpin1Max;
  UINT4  mObsCoh;
  INT8   f0Bin, fLastBin, fBin;
  REAL8  alpha, delta, timeBase, deltaF, f1jump;
  REAL8  patchSizeX, patchSizeY;
  UINT2  xSide, ySide;
  UINT2  maxNBins, maxNBorders;

  /* user input variables */
  BOOLEAN  uvar_help, uvar_weighAM, uvar_weighNoise, uvar_printLog, uvar_printWeights;
  INT4     uvar_blocksRngMed, uvar_nfSizeCylinder, uvar_maxBinsClean;
  REAL8    uvar_f0, uvar_peakThreshold, uvar_houghFalseAlarm, uvar_fSearchBand;
  CHAR     *uvar_earthEphemeris=NULL;
  CHAR     *uvar_sunEphemeris=NULL;
  CHAR     *uvar_sftDir=NULL;
  CHAR     *uvar_dirnameOut=NULL;
  CHAR     *uvar_fbasenameOut=NULL;
  CHAR     *uvar_skyfile=NULL;
  LALStringVector *uvar_linefiles=NULL;


  /* Set up the default parameters */
  
  lalDebugLevel = 0;  /* LALDebugLevel must be called before anything else */
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);
  
  uvar_help = FALSE;
  uvar_weighAM = TRUE;
  uvar_weighNoise = TRUE;
  uvar_printLog = FALSE;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_nfSizeCylinder = NFSIZE;
  uvar_f0 = F0;
  uvar_fSearchBand = FBAND;
  uvar_peakThreshold = THRESHOLD;
  uvar_houghFalseAlarm = FALSEALARM;
  uvar_printEvents = FALSE;
  uvar_printTemplates = FALSE;
  uvar_printMaps = FALSE;
  uvar_printStats = FALSE;
  uvar_printSigma = FALSE;
  uvar_maxBinsClean = 100;
  uvar_printWeights = FALSE;

  uvar_earthEphemeris = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);

  uvar_dirnameOut = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_dirnameOut,DIROUT);

  uvar_fbasenameOut = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_fbasenameOut,BASENAMEOUT);

  uvar_skyfile = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_skyfile,SKYFILE);

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",                    &uvar_help),            &status);  
  LAL_CALL( LALRegisterREALUserVar(   &status, "f0",              'f', UVAR_OPTIONAL, "Start search frequency",                &uvar_f0),              &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fSearchBand",     'b', UVAR_OPTIONAL, "Search frequency band",                 &uvar_fSearchBand),     &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyfile",          0,  UVAR_OPTIONAL, "Input skypatch file",                   &uvar_skyfile),         &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThreshold",    0,  UVAR_OPTIONAL, "Peak selection threshold",              &uvar_peakThreshold),   &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "weighAM",          0,  UVAR_OPTIONAL, "Use amplitude modulation weights",      &uvar_weighAM),         &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "weighNoise",       0,  UVAR_OPTIONAL, "Use SFT noise weights",                 &uvar_weighNoise),      &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printLog",         0,  UVAR_OPTIONAL, "Print Log file",                        &uvar_printLog),        &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_OPTIONAL, "Earth Ephemeris file",                  &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_OPTIONAL, "Sun Ephemeris file",                    &uvar_sunEphemeris),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",          'D', UVAR_REQUIRED, "SFT filename pattern",                  &uvar_sftDir),          &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "dirnameOut",      'o', UVAR_OPTIONAL, "Output directory",                      &uvar_dirnameOut),      &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fbasenameOut",     0,  UVAR_OPTIONAL, "Output file basename",                  &uvar_fbasenameOut),    &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printMaps",        0,  UVAR_OPTIONAL, "Print Hough maps",                      &uvar_printMaps),       &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printTemplates",   0,  UVAR_OPTIONAL, "Print templates file",                  &uvar_printTemplates),  &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "houghFalseAlarm",  0,  UVAR_OPTIONAL, "Hough false alarm to set threshold",    &uvar_houghFalseAlarm), &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printEvents",      0,  UVAR_OPTIONAL, "Print events above threshold",          &uvar_printEvents),     &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printStats",       0,  UVAR_OPTIONAL, "Print Hough statistics",                &uvar_printStats),      &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printSigma",       0,  UVAR_OPTIONAL, "Print expected number count stdev.",    &uvar_printSigma),      &status);
  LAL_CALL( LALRegisterLISTUserVar(   &status, "linefiles",        0,  UVAR_OPTIONAL, "Comma separated List of linefiles (filenames must contain IFO name)",
				                                                                                               &uvar_linefiles),       &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nfSizeCylinder",   0,  UVAR_OPTIONAL, "Size of cylinder of PHMDs",             &uvar_nfSizeCylinder),  &status);

  /* developer input variables */
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",     0, UVAR_DEVELOPER, "Running Median block size",             &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "maxBinsClean",     0, UVAR_DEVELOPER, "Maximum number of bins in cleaning",    &uvar_maxBinsClean),    &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printWeights",     0, UVAR_DEVELOPER, "Print relative noise weights of ifos",  &uvar_printWeights),    &status);  

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 


  /* write log file with command line arguments, cvs tags, and contents of skypatch file */
  if ( uvar_printLog ) {
    LAL_CALL( PrintLogFile( &status, uvar_dirnameOut, uvar_fbasenameOut, uvar_skyfile, uvar_linefiles, argv[0]), &status);
  }

  /***** start main calculations *****/



  /* read sft Files and set up weights and nstar vector */
  {
    /* new SFT I/O data types */
    SFTCatalog *catalog = NULL;
    static SFTConstraints constraints;

    REAL8 doppWings, fmin, fmax;
    INT4 length;

    /* set detector constraint */
    constraints.detector = NULL;

    /* get sft catalog */
    LAL_CALL( LALSFTdataFind( &status, &catalog, uvar_sftDir, &constraints), &status);
    if ( (catalog == NULL) || (catalog->length == 0) ) {
      fprintf (stderr,"Unable to match any SFTs with pattern '%s'\n", uvar_sftDir );
      exit(1);
    }

    /* get some sft parameters */
    mObsCoh = catalog->length; /* number of sfts */
    deltaF = catalog->data->header.deltaF;  /* frequency resolution */
    timeBase= 1.0/deltaF; /* coherent integration time */
    f0Bin = floor( uvar_f0 * timeBase + 0.5); /* initial search frequency */
    length =  uvar_fSearchBand * timeBase; /* total number of search bins - 1 */
    fLastBin = f0Bin + length;   /* final frequency bin to be analyzed */
  
    /* add wings for Doppler modulation and running median block size*/
    doppWings = (uvar_f0 + uvar_fSearchBand) * VTOT;    
    fmin = uvar_f0 - doppWings - (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;
    fmax = uvar_f0 + uvar_fSearchBand + doppWings + (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;

    /* read sft files making sure to add extra bins for running median */
    /* read the sfts */
    LAL_CALL( LALLoadMultiSFTs ( &status, &inputSFTs, catalog, fmin, fmax), &status);


    /* clean sfts if required */
    if ( LALUserVarWasSet( &uvar_linefiles ) )
      {
	RandomParams *randPar=NULL;
	FILE *fpRand=NULL;
	INT4 seed, ranCount;  

	if ( (fpRand = fopen("/dev/urandom", "r")) == NULL ) {
	  fprintf(stderr,"Error in opening /dev/urandom" ); 
	  exit(1);
	} 

	if ( (ranCount = fread(&seed, sizeof(seed), 1, fpRand)) != 1 ) {
	  fprintf(stderr,"Error in getting random seed" );
	  exit(1);
	}

	LAL_CALL ( LALCreateRandomParams (&status, &randPar, seed), &status );

	LAL_CALL( LALRemoveKnownLinesInMultiSFTVector ( &status, inputSFTs, uvar_maxBinsClean, uvar_blocksRngMed, uvar_linefiles, randPar), &status);

	LAL_CALL ( LALDestroyRandomParams (&status, &randPar), &status);
	fclose(fpRand);
      } /* end cleaning */


    /* SFT info -- assume all SFTs have same length */
    numifo = inputSFTs->length;
    binsSFT = inputSFTs->data[0]->data->data->length;

    LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);  	

    /* set up weights -- this should be done before normalizing the sfts */
    weightsNoise.length = mObsCoh;
    weightsNoise.data = (REAL8 *)LALCalloc(mObsCoh, sizeof(REAL8));

    /* initialize all weights to unity */
    LAL_CALL( LALHOUGHInitializeWeights( &status, &weightsNoise), &status);


  } /* end of sft reading block */



  /* get weights vector */
  { 
    MultiNoiseWeights *multweight = NULL;    
    MultiPSDVector *multPSD = NULL;  
    REAL8 dmpNormalization;
    INT4 tmpLeap;
    UINT4 iIFO, iSFT, numsft, j;


    /* normalize sfts */
    LAL_CALL( LALNormalizeMultiSFTVect (&status, &multPSD, inputSFTs, uvar_blocksRngMed), &status);
   
    /* compute multi noise weights */
    if ( uvar_weighNoise ) {
      LAL_CALL ( LALComputeMultiNoiseWeights ( &status, &multweight, &dmpNormalization, multPSD, uvar_blocksRngMed, 0), &status);
    }
    
    /* we are now done with the psd */
    LAL_CALL ( LALDestroyMultiPSDVector  ( &status, &multPSD), &status);



    /* copy the timestamps, weights, and velocity vector */
    for (j = 0, iIFO = 0; iIFO < numifo; iIFO++ ) {

      numsft = multweight->data[iIFO]->length;
      
      for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {


	if ( uvar_weighNoise )
	  weightsNoise.data[j] = multweight->data[iIFO]->data[iSFT];


      } /* loop over SFTs */

    } /* loop over IFOs */

    if ( uvar_weighNoise ) {
      LAL_CALL( LALHOUGHNormalizeWeights( &status, &weightsNoise), &status);
    }

    if ( uvar_weighNoise ) {    
      LAL_CALL ( LALDestroyMultiNoiseWeights ( &status, &multweight), &status);
    }

  } /* end block for weights, velocity and time */
  

  /* print relative weights of ifos to stdout */  
  if ( uvar_printWeights )
    {
      UINT4 iIFO, iSFT, numsft, j;
      REAL8 *sumweights=NULL;
    
      sumweights = (REAL8 *)LALCalloc(1, numifo*sizeof(REAL8));
      for (j = 0, iIFO = 0; iIFO < numifo; iIFO++ ) {
	
	numsft = mdetStates->data[iIFO]->length;
	
	for ( iSFT = 0; iSFT < numsft; iSFT++, j++) 	  
	  sumweights[iIFO] += weightsNoise.data[j];	
      
      } /* loop over IFOs */
      
      for ( iIFO = 0; iIFO < numifo; iIFO++ )
	fprintf(stdout, "%d  %f\n", iIFO, sumweights[iIFO]);
       
     LALFree(sumweights);
      
    } /* end debugging */
  

  LALFree(weightsNoise.data);  
 

  LAL_CALL (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  if ( lalDebugLevel )
    REPORTSTATUS ( &status);

  return status.statusCode;
}



  

void PrintLogFile (LALStatus       *status,
		   CHAR            *dir,
		   CHAR            *basename,
		   CHAR            *skyfile,
		   LALStringVector *linefiles,
		   CHAR            *executable )
{
  CHAR *fnameLog=NULL; 
  FILE *fpLog=NULL;
  CHAR *logstr=NULL; 
  UINT4 k;

  INITSTATUS (status, "PrintLogFile", rcsid);
  ATTATCHSTATUSPTR (status);
  
  /* open log file for writing */
  fnameLog = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(fnameLog,dir);
  strcat(fnameLog, "/logfiles/");
  /* now create directory fdirOut/logfiles using mkdir */
  errno = 0;
  {
    /* check whether file can be created or if it exists already 
       if not then exit */
    INT4 mkdir_result;
    mkdir_result = mkdir(fnameLog, S_IRWXU | S_IRWXG | S_IRWXO);
    if ( (mkdir_result == -1) && (errno != EEXIST) )
      {
	fprintf(stderr, "unable to create logfiles directory %s\n", fnameLog);
        LALFree(fnameLog);
	exit(1);  /* stop the program */
      }
  }

  /* create the logfilename in the logdirectory */
  strcat(fnameLog, basename);
  strcat(fnameLog,".log");
  /* open the log file for writing */
  if ((fpLog = fopen(fnameLog, "w")) == NULL) {
    fprintf(stderr, "Unable to open file %s for writing\n", fnameLog);
    LALFree(fnameLog);
    exit(1);
  }
  
  /* get the log string */
  TRY( LALUserVarGetLog(status->statusPtr, &logstr, UVAR_LOGFMT_CFGFILE), status);  

  fprintf( fpLog, "## LOG FILE FOR Hough Driver\n\n");
  fprintf( fpLog, "# User Input:\n");
  fprintf( fpLog, "#-------------------------------------------\n");
  fprintf( fpLog, logstr);
  LALFree(logstr);

  /* copy contents of skypatch file into logfile */
  fprintf(fpLog, "\n\n# Contents of skypatch file:\n");
  fclose(fpLog);
  {
    CHAR command[1024] = "";
    sprintf(command, "cat %s >> %s", skyfile, fnameLog);
    system(command);
  }


  /* copy contents of linefile if necessary */
  if ( linefiles ) {

    for ( k = 0; k < linefiles->length; k++) {
      
      if ((fpLog = fopen(fnameLog, "a")) != NULL) {
	CHAR command[1024] = "";
	fprintf (fpLog, "\n\n# Contents of linefile %s :\n", linefiles->data[k]);
	fprintf (fpLog, "# -----------------------------------------\n");
	fclose (fpLog);
	sprintf(command, "cat %s >> %s", linefiles->data[k], fnameLog);      
	system (command);	 
      } 
    } 
  }

  /* append an ident-string defining the exact CVS-version of the code used */
  if ((fpLog = fopen(fnameLog, "a")) != NULL) 
    {
      CHAR command[1024] = "";
      fprintf (fpLog, "\n\n# CVS-versions of executable:\n");
      fprintf (fpLog, "# -----------------------------------------\n");
      fclose (fpLog);
      
      sprintf (command, "ident %s | sort -u >> %s", executable, fnameLog);
      system (command);	/* we don't check this. If it fails, we assume that */
    			/* one of the system-commands was not available, and */
    			/* therefore the CVS-versions will not be logged */ 
    }

  LALFree(fnameLog); 
  	 
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}    


