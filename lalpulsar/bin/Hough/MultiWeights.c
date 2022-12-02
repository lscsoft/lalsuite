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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
 *  MA  02110-1301  USA
 */

/**
 * \file MultiWeights.c
 * \ingroup lalpulsar_bin_Hough
 * \author Badri Krishnan, Alicia Sintes
 * \brief Utility code for calculating some properties of the noise weights
 * for the multi-IFO Hough search
 *
 * History:   Created by Sintes and Krishnan July 04, 2003
 * Modifications for S4 January 2006
 *
 */

#include "DriveHoughColor.h"

/* globals, constants and defaults */

/* boolean global variables for controlling output */
BOOLEAN uvar_printEvents, uvar_printTemplates, uvar_printMaps, uvar_printStats, uvar_printSigma;

#define EARTHEPHEMERIS "./earth05-09.dat" 
#define SUNEPHEMERIS "./sun05-09.dat"    

/* #define EARTHEPHEMERIS "./earth00-04.dat" */
/* #define SUNEPHEMERIS "./sun00-04.dat" */

#define MAXFILENAMELENGTH 512 /* maximum # of characters  of a filename */

#define DIROUT "./outMulti"   /* output directory */
#define BASENAMEOUT "HM"    /* prefix file output */
#define TEMPOUT "./tempout"    /* output file */

#define THRESHOLD 1.6 /* thresold for peak selection, with respect to the
                              the averaged power in the search band */
#define FALSEALARM 1.0e-9 /* Hough false alarm for candidate selection */
#define SKYFILE "./sky1"      
#define F0 505.0   /*  frequency to build the LUT and start search */
#define FBAND 0.05   /* search frequency band  (in Hz) */
#define NFSIZE  51   /* n-freq. span of the cylinder, to account for spin-down search */
#define BLOCKSRNGMED 101 /* Running median window size */

#define TRUE (1==1)
#define FALSE (1==0)

/* local function prototype */
void PrintLogFile (LALStatus *status, LALStringVector *linefiles, CHAR *executable );



/*******************************************************/

int main(int argc, char *argv[]){

  /* LALStatus pointer */
  static LALStatus  status;  

  REAL8 scalepow=1.0e40;
  
  /* ephemeris */
  EphemerisData    *edat=NULL;
 
  /* standard pulsar sft types */ 
  MultiSFTVector *inputSFTs = NULL;
 
   /* sft constraint variables */
  LIGOTimeGPS startTimeGPS, endTimeGPS;
  LIGOTimeGPSVector *inputTimeStampsVector=NULL;
 
 
 /* information about all the ifos */
  MultiDetectorStateSeries *mdetStates = NULL;
  UINT4 numifo;

 /* vector of weights */
  REAL8Vector weightsV;
     
  UINT4 numsft;
  INT4 k;
 
  /* miscellaneous */
  UINT4  mObsCoh;
  REAL8  deltaF;

  /* user input variables */
  BOOLEAN  uvar_printLog;
  BOOLEAN  uvar_weighAM, uvar_weighNoise;
  BOOLEAN  uvar_dumpAllW, uvar_dumpRelW, uvar_dumpNoise;
  CHAR     *uvar_sftDir=NULL;
  REAL8    uvar_startTime, uvar_endTime;
  CHAR     *uvar_timeStampsFile=NULL;
  LALStringVector *uvar_linefiles=NULL;
  REAL8    uvar_AlphaWeight, uvar_DeltaWeight;
  CHAR     *uvar_earthEphemeris=NULL;
  CHAR     *uvar_sunEphemeris=NULL;
  INT4     uvar_blocksRngMed, uvar_nfSizeCylinder, uvar_maxBinsClean;
  REAL8    uvar_f0, uvar_fSearchBand;
  CHAR     *uvar_outfile=NULL;

  /* Set up the default parameters */  

  /* LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;
 
  uvar_printLog = FALSE;
  uvar_weighAM = FALSE;
  uvar_weighNoise = FALSE; 
  uvar_dumpAllW = FALSE;
  uvar_dumpRelW = FALSE;
  uvar_dumpNoise= TRUE;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_f0 = F0;
  uvar_fSearchBand = FBAND;
  uvar_maxBinsClean = 100;
  uvar_startTime= 0;
  uvar_endTime = LAL_INT4_MAX;

  /* nfsizecylinder is really not required here but retained
     just to make sure that the number of bins read in the sft
     is the same as in the driver */
  uvar_nfSizeCylinder = NFSIZE;  
  uvar_AlphaWeight = 1.0;
  uvar_DeltaWeight = 1.0;

  uvar_outfile = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_outfile,TEMPOUT);

  uvar_earthEphemeris = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);
  
  /* register user input variables */
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_f0,             "f0",                  REAL8,        'f', OPTIONAL,  "Start search frequency") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_fSearchBand,    "fSearchBand",         REAL8,        'b', OPTIONAL,  "Search frequency band") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_printLog,       "printLog",            BOOLEAN,      0,   OPTIONAL,  "Print Log file") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_sftDir,         "sftDir",              STRING,       'D', REQUIRED,  "SFT filename pattern. Possibilities are:\n"
                                          " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_linefiles,      "linefiles",           STRINGVector, 0,   OPTIONAL,  "Comma separated List of linefiles (filenames must contain IFO name)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_startTime,      "startTime",           REAL8,        0,   OPTIONAL,  "GPS start time of observation (SFT timestamps must be >= this)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_endTime,        "endTime",             REAL8,        0,   OPTIONAL,  "GPS end time of observation (SFT timestamps must be < this)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_timeStampsFile, "timeStampsFile",      STRING,       0,   OPTIONAL,  "Input time-stamps file") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_weighAM,        "weightAM",            BOOLEAN,      0,   OPTIONAL,  "Use amplitude modulation weights") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_weighNoise,     "weightNoise",         BOOLEAN,      0,   OPTIONAL,  "Use SFT noise weights") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_dumpAllW,       "dumpAllWeights",      BOOLEAN,      0,   OPTIONAL,  "Dump all weights") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_dumpRelW,       "dumpRelativeWeights", BOOLEAN,      0,   OPTIONAL,  "Dump IFO relative weights") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_dumpNoise,      "dumpNoise",           BOOLEAN,      0,   OPTIONAL,  "Dump Noise estimate") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_earthEphemeris, "earthEphemeris",      STRING,       'E', OPTIONAL,  "Earth Ephemeris file") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_sunEphemeris,   "sunEphemeris",        STRING,       'S', OPTIONAL,  "Sun Ephemeris file") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_AlphaWeight,    "AlphaWeight",         REAL8,        0,   OPTIONAL,  "sky Alpha for weight calculation") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_DeltaWeight,    "DeltaWeight",         REAL8,        0,   OPTIONAL,  "sky Delta for weight calculation") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_outfile,        "outfile",             STRING,       0,   OPTIONAL,  "output file name") == XLAL_SUCCESS, XLAL_EFUNC);
 /* developer input variables */
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_blocksRngMed,   "blocksRngMed",        INT4,         0,   DEVELOPER, "Running Median block size") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_maxBinsClean,   "maxBinsClean",        INT4,         0,   DEVELOPER, "Maximum number of bins in cleaning") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_nfSizeCylinder, "nfSizeCylinder",      INT4,         0,   DEVELOPER, "Size of cylinder of PHMDs") == XLAL_SUCCESS, XLAL_EFUNC);

  /* read all command line variables */
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput(&should_exit, argc, argv, lalPulsarVCSInfoList) == XLAL_SUCCESS, XLAL_EFUNC);
  if (should_exit)
    exit(1);

 /* very basic consistency checks on user input */
  if ( uvar_f0 < 0 ) {
    fprintf(stderr, "start frequency must be positive\n");
    exit(1);
  }

  if ( uvar_fSearchBand < 0 ) {
    fprintf(stderr, "search frequency band must be positive\n");
    exit(1);
  }
 
   /* write log file with command line arguments, cvs tags, and contents of skypatch file */
  if ( uvar_printLog ) {
    LAL_CALL( PrintLogFile( &status, uvar_linefiles, argv[0]), &status);
  }
    
  /* read sfts and clean them */
  {
    /* new SFT I/O data types */
    SFTCatalog *catalog = NULL;
    static SFTConstraints constraints;
    
    REAL8 doppWings, f_min, f_max;

    /* set detector constraint */
    constraints.detector = NULL;

    if ( XLALUserVarWasSet( &uvar_startTime ) ) {
      XLALGPSSetREAL8(&startTimeGPS, uvar_startTime);
      constraints.minStartTime = &startTimeGPS;
    }

    if ( XLALUserVarWasSet( &uvar_endTime ) ) {
      XLALGPSSetREAL8(&endTimeGPS, uvar_endTime);
      constraints.maxStartTime = &endTimeGPS;
    }

    if ( XLALUserVarWasSet( &uvar_timeStampsFile ) ) {
      XLAL_CHECK_MAIN ( ( inputTimeStampsVector = XLALReadTimestampsFile ( uvar_timeStampsFile) ) != NULL, XLAL_EFUNC);
      constraints.timestamps = inputTimeStampsVector;
    }

    /* get sft catalog */
    XLAL_CHECK_MAIN( ( catalog = XLALSFTdataFind( uvar_sftDir, &constraints) ) != NULL, XLAL_EFUNC);
    if ( (catalog == NULL) || (catalog->length == 0) ) {
      fprintf (stderr,"Unable to match any SFTs with pattern '%s'\n", uvar_sftDir );
      exit(1);
    }
    
    /* now we can free the inputTimeStampsVector */
    if ( XLALUserVarWasSet( &uvar_timeStampsFile ) ) {
      XLALDestroyTimestampVector(inputTimeStampsVector );
    }
    
    /* catalog is ordered in time so we can get start, end time and tObs*/
    mObsCoh = catalog->length; /* not always correct number of sfts */ 
    deltaF = catalog->data->header.deltaF;  /* frequency resolution */
  
    /* add wings for Doppler modulation and running median block size*/
    doppWings = (uvar_f0 + uvar_fSearchBand) * VTOT;    
    f_min = uvar_f0 - doppWings - (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;
    f_max = uvar_f0 + uvar_fSearchBand + doppWings + (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;

    /* read the sfts */
    XLAL_CHECK_MAIN( ( inputSFTs = XLALLoadMultiSFTs ( catalog, f_min, f_max) ) != NULL, XLAL_EFUNC);
    numifo = inputSFTs->length;

    /* find number of sfts */     
    /* mObsCoh = catalog->length; not always correct number of sfts */   
    /* loop over ifos and calculate number of sfts */
    /* note that we can't use the catalog to determine the number of SFTs
       because SFTs might be segmented in frequency */
    mObsCoh = 0; /* initialization */
    for (k = 0; k < (INT4)numifo; k++ ) {
      mObsCoh += inputSFTs->data[k]->length;
    } 

    /* clean sfts if required */
    if ( XLALUserVarWasSet( &uvar_linefiles ) )
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

    XLALDestroySFTCatalog(catalog );  	
    
  } /* end of sft reading block */
    
    /*******************************************************/

  /* get detector velocities, positions, weights vector */ 
  { 
    MultiNoiseWeights *multweight = NULL;    
    MultiPSDVector *multPSD = NULL;  
    UINT4 iIFO, iSFT, j;

    /*  get ephemeris  */
    XLAL_CHECK_MAIN( ( edat = XLALInitBarycenter( uvar_earthEphemeris, uvar_sunEphemeris ) ) != NULL, XLAL_EFUNC);

   /* set up weights */
    weightsV.length = mObsCoh;
    weightsV.data = (REAL8 *)LALCalloc(1, mObsCoh * sizeof(REAL8));

    /* initialize all weights to unity */
    LAL_CALL( LALHOUGHInitializeWeights( &status, &weightsV), &status);
 

    /* get information about all detectors including velocity and timestamps */
    /* note that this function returns the velocity at the 
       mid-time of the SFTs -- should not make any difference */
    const REAL8 tOffset = 0.5 / inputSFTs->data[0]->data[0].deltaF;
    XLAL_CHECK_MAIN ( ( mdetStates = XLALGetMultiDetectorStatesFromMultiSFTs ( inputSFTs, edat, tOffset ) ) != NULL, XLAL_EFUNC);

    /* normalize sfts and get power running-median rngmed[ |data|^2] from SFTs */
    XLAL_CHECK_MAIN( ( multPSD = XLALNormalizeMultiSFTVect(  inputSFTs, uvar_blocksRngMed, NULL ) ) != NULL, XLAL_EFUNC); 

    if ( uvar_weighNoise ) {      
      /* compute multi noise weights if required */ 
      XLAL_CHECK_MAIN ( ( multweight = XLALComputeMultiNoiseWeights ( multPSD, uvar_blocksRngMed, 0) ) != NULL, XLAL_EFUNC);
      /* copy the weights */
      for (j = 0, iIFO = 0; iIFO < numifo; iIFO++ ) {
        numsft = mdetStates->data[iIFO]->length;     
        for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {
	  weightsV.data[j] = multweight->data[iIFO]->data[iSFT];
        } /* loop over SFTs */
      } /* loop over IFOs */

      LAL_CALL( LALHOUGHNormalizeWeights( &status, &weightsV), &status);     
      XLALDestroyMultiNoiseWeights (multweight);
    }
 
    if (uvar_dumpNoise)  { 
    /* calculate (1/ (sum of 1/Sn^2) )^(1/4) where Sn is avg psd of each sft */ 
      REAL8 sumSn, normPSD,sumSnInv=0.0;
      INT8 binsSFT;

      /* nomalization factor to get proper single-sided PSD: Sn=(2/Tsft) rngmed[ |data|^2]*/
      normPSD=2.0*deltaF;
	
      for ( iIFO = 0; iIFO < numifo; iIFO++ ) {
        numsft = multPSD->data[iIFO]->length;
      
        for ( iSFT = 0; iSFT < numsft; iSFT++) {	
	  binsSFT = multPSD->data[iIFO]->data[iSFT].data->length;
	  sumSn = 0.0;
	  /* use a scale to make numbers closer to unity */
	  for ( j = 0; j < binsSFT; j++)
	    sumSn += scalepow * multPSD->data[iIFO]->data[iSFT].data->data[j];
	  
	  sumSnInv += binsSFT * binsSFT / (sumSn * sumSn);
        } /* end loop over sfts */     
      } /* end loop over IFOs */
       
      sumSnInv = sqrt(sumSnInv);
      sumSnInv *= (scalepow/normPSD);
      sumSnInv = sqrt(sumSnInv);  
      fprintf(stdout, "%f  %g\n", uvar_f0, 1.0/sumSnInv);
    } /* end block for Sn-weights calculation */

    /* we are now done with the psd */
    XLALDestroyMultiPSDVector  ( multPSD);

  } /* end block for noise weights, velocity and time */
  

     /* calculate amplitude modulation weights if required */
  if (uvar_weighAM) 
    {
      MultiAMCoeffs *multiAMcoef = NULL;
      UINT4 iIFO, iSFT;
      SkyPosition skypos;      

      /* get the amplitude modulation coefficients */
      skypos.longitude = uvar_AlphaWeight;
      skypos.latitude = uvar_DeltaWeight;
      skypos.system = COORDINATESYSTEM_EQUATORIAL;
      XLAL_CHECK_MAIN ( ( multiAMcoef = XLALComputeMultiAMCoeffs ( mdetStates, NULL, skypos) ) != NULL, XLAL_EFUNC);
      
      /* loop over the weights and multiply them by the appropriate
	 AM coefficients */
      for ( k = 0, iIFO = 0; iIFO < numifo; iIFO++) {	
	numsft = mdetStates->data[iIFO]->length;	

	for ( iSFT = 0; iSFT < numsft; iSFT++, k++) {	  	  
	  REAL8 a, b;
	  
	  a = multiAMcoef->data[iIFO]->a->data[iSFT];
	  b = multiAMcoef->data[iIFO]->b->data[iSFT];    
	  weightsV.data[k] *= (a*a + b*b);
	} /* loop over SFTs */
      } /* loop over IFOs */
      
      LAL_CALL( LALHOUGHNormalizeWeights( &status, &weightsV), &status);
      
      XLALDestroyMultiAMCoeffs ( multiAMcoef );
    } /* end AM weights calculation */

  /*******************************************************/
  /* dump weight vector if required */
  if (uvar_dumpAllW) 
    {
      FILE *fp=NULL;
       
      fp = fopen(uvar_outfile   , "w");
      setvbuf(fp, (char *)NULL, _IOLBF, 0);
      
      for (k=0; k < (INT4)mObsCoh; k++){
        fprintf(fp, "%g  \n", weightsV.data[k]);
      }
      fclose(fp);
    }
 
  /*******************************************************/

   /* print relative weights of ifos to stdout if required */ 
   if (uvar_dumpRelW) 
    {
      REAL8 *sumweights=NULL;    
      sumweights = (REAL8 *)LALCalloc(1, numifo*sizeof(REAL8));
      UINT4 j, iIFO, iSFT;
      for (j=0, iIFO = 0; iIFO < numifo; iIFO++ ) {      
        numsft = mdetStates->data[iIFO]->length;
      
        for ( iSFT = 0; iSFT < numsft; iSFT++, j++) 	  
	  sumweights[iIFO] += weightsV.data[j];
      
      } /* end loop over IFOs */

      /*print relative sum of weights */
      for ( iIFO = 0; iIFO < numifo; iIFO++ )
        fprintf(stdout, "%s  %f\n", inputSFTs->data[iIFO]->data[0].name, sumweights[iIFO]/mObsCoh);
    
      LALFree(sumweights);
      
    } /* end printing of relative weights */

  /*******************************************************/

  /* free memory and exit */
  LALFree(weightsV.data);
  XLALDestroyMultiDetectorStateSeries ( mdetStates );

  XLALDestroyEphemerisData(edat);
    
  XLALDestroyMultiSFTVector( inputSFTs);

  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  if ( lalDebugLevel )
    REPORTSTATUS ( &status);

  return status.statusCode;
}

/*******************************************************/  

void PrintLogFile (LALStatus       *status,
		   LALStringVector *linefiles,
		   CHAR            *executable )
{
  CHAR *fnameLog=NULL; 
  FILE *fpLog=NULL;
  CHAR *logstr=NULL; 
  UINT4 k;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /* open log file for writing */
  fnameLog = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(fnameLog,"./MultiWeight.log");

  /* open the log file for writing */
  if ((fpLog = fopen(fnameLog, "w")) == NULL) {
    fprintf(stderr, "Unable to open file %s for writing\n", fnameLog);
    LALFree(fnameLog);
    exit(1);
  }
  
  /* get the log string */
  XLAL_CHECK_LAL( status, ( logstr = XLALUserVarGetLog(UVAR_LOGFMT_CFGFILE) ) != NULL, XLAL_EFUNC);  

  fprintf( fpLog, "## LOG FILE FOR Hough Driver\n\n");
  fprintf( fpLog, "# User Input:\n");
  fprintf( fpLog, "#-------------------------------------------\n");
  fprintf( fpLog, "%s", logstr);
  LALFree(logstr);

  fclose(fpLog);


  /* copy contents of linefile if necessary */
  if ( linefiles ) {

    for ( k = 0; k < linefiles->length; k++) {
      
      if ((fpLog = fopen(fnameLog, "a")) != NULL) {
	CHAR command[1024] = "";
	fprintf (fpLog, "\n\n# Contents of linefile %s :\n", linefiles->data[k]);
	fprintf (fpLog, "# -----------------------------------------\n");
	fclose (fpLog);
	sprintf(command, "cat %s >> %s", linefiles->data[k], fnameLog);      
        if ( system(command) ) fprintf (stderr, "\nsystem('%s') returned non-zero status!\n\n", command );
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
      /* we don't check this. If it fails, we assume that */
      /* one of the system-commands was not available, and */
      /* therefore the CVS-versions will not be logged */
      if ( system(command) ) fprintf (stderr, "\nsystem('%s') returned non-zero status!\n\n", command );
    }

  LALFree(fnameLog); 
  	 
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}    


