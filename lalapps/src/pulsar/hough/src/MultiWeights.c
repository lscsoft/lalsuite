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
 * \file MultiWeights.c
 * \author Badri Krishnan, Alicia Sintes 
 * \brief Utility code for calculating some properties of the noise weights 
   for the multi-IFO Hough search

   Revision: $Id$
 
   History:   Created by Sintes and Krishnan July 04, 2003
              Modifications for S4 January 2006

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
void PrintLogFile (LALStatus *status, LALStringVector *linefiles, CHAR *executable );



/******************************************/

int main(int argc, char *argv[]){

  /* LALStatus pointer */
  static LALStatus  status;  

  /* LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;
  
  /* standard pulsar sft types */ 
  MultiSFTVector *inputSFTs = NULL;
  UINT4 numifo;


  MultiNoiseWeights *multweight = NULL;    
  MultiPSDVector *multPSD = NULL;  
  REAL8 dmpNormalization;
  UINT4 iIFO, iSFT, numsft;
  

  /* miscellaneous */
  UINT4  mObsCoh;
  REAL8  deltaF;

  /* user input variables */
  BOOLEAN  uvar_help, uvar_printLog;
  INT4     uvar_blocksRngMed, uvar_nfSizeCylinder, uvar_maxBinsClean;
  REAL8    uvar_f0, uvar_fSearchBand;
  CHAR     *uvar_sftDir=NULL;
  LALStringVector *uvar_linefiles=NULL;


  /* Set up the default parameters */  
  lalDebugLevel = 0;  /* LALDebugLevel must be called before anything else */
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);  
  uvar_help = FALSE;
  uvar_printLog = FALSE;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_f0 = F0;
  uvar_fSearchBand = FBAND;
  uvar_maxBinsClean = 100;

  /* nfsizecylinder is really not required here but retained
     just to make sure that the number of bins read in the sft
     is the same as in the driver */
  uvar_nfSizeCylinder = NFSIZE; 

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",                    &uvar_help),            &status);  
  LAL_CALL( LALRegisterREALUserVar(   &status, "f0",              'f', UVAR_OPTIONAL, "Start search frequency",                &uvar_f0),              &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fSearchBand",     'b', UVAR_OPTIONAL, "Search frequency band",                 &uvar_fSearchBand),     &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printLog",         0,  UVAR_OPTIONAL, "Print Log file",                        &uvar_printLog),        &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",          'D', UVAR_REQUIRED, "SFT filename pattern",                  &uvar_sftDir),          &status);
  LAL_CALL( LALRegisterLISTUserVar(   &status, "linefiles",        0,  UVAR_OPTIONAL, "Comma separated List of linefiles (filenames must contain IFO name)",
				      &uvar_linefiles),       &status);		      
  /* developer input variables */
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",     0, UVAR_DEVELOPER, "Running Median block size",             &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "maxBinsClean",     0, UVAR_DEVELOPER, "Maximum number of bins in cleaning",    &uvar_maxBinsClean),    &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nfSizeCylinder",   0, UVAR_DEVELOPER, "Size of cylinder of PHMDs",             &uvar_nfSizeCylinder),  &status);

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 


  /* write log file with command line arguments, cvs tags, and contents of skypatch file */
  if ( uvar_printLog ) {
    LAL_CALL( PrintLogFile( &status, uvar_linefiles, argv[0]), &status);
  }


  /* read sfts and clean them */
  {
    /* new SFT I/O data types */
    SFTCatalog *catalog = NULL;
    static SFTConstraints constraints;
    REAL8 doppWings, fmin, fmax;

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

    LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);  	
    
  } /* end of sft reading block */
    
  /* normalize sfts and compute psds */
  LAL_CALL( LALNormalizeMultiSFTVect (&status, &multPSD, inputSFTs, uvar_blocksRngMed), &status);

  /* calculate sum of 1/Sn where Sn is avg psd of each sft */
  {
    REAL8 sumSn, sumSnInv=0.0;
    INT8 binsSFT;
    INT4 j;

    for ( iIFO = 0; iIFO < numifo; iIFO++ ) {
      numsft = multPSD->data[iIFO]->length;
      
      for ( iSFT = 0; iSFT < numsft; iSFT++) {
	
	binsSFT = multPSD->data[iIFO]->data[iSFT].data->length;
	sumSn = 0.0;
	for ( j = 0; j < binsSFT; j++)
	  sumSn += multPSD->data[iIFO]->data[iSFT].data->data[j];

	sumSnInv += binsSFT/sumSn;
      } /* end loop over sfts */
      
    } /* end loop over IFOs */
    

    fprintf(stdout, "%g\n", sumSnInv);
  } /* end block for 1/Sn calculation */


  /* compute multi noise weights */
  LAL_CALL ( LALComputeMultiNoiseWeights ( &status, &multweight, &dmpNormalization, multPSD, uvar_blocksRngMed, 0), &status);
    
  /* we are now done with the psd */
  LAL_CALL ( LALDestroyMultiPSDVector  ( &status, &multPSD), &status);
  
  LAL_CALL ( LALDestroyMultiNoiseWeights ( &status, &multweight), &status);
  
  
  
  {
    /* print relative weights of ifos to stdout */ 
    REAL8 *sumweights=NULL;
    
    sumweights = (REAL8 *)LALCalloc(1, numifo*sizeof(REAL8));
    for ( iIFO = 0; iIFO < numifo; iIFO++ ) {
      
      numsft = multweight->data[iIFO]->length;
      
      for ( iSFT = 0; iSFT < numsft; iSFT++) 	  
	sumweights[iIFO] += multweight->data[iIFO]->data[iSFT];
      
    } /* end loop over IFOs */

    /* print relative sum of weights */    
    for ( iIFO = 0; iIFO < numifo; iIFO++ )
      fprintf(stdout, "%d  %f\n", iIFO, sumweights[iIFO]);
    
    LALFree(sumweights);
      
  } /* end debugging */
				      
 

  LAL_CALL (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  if ( lalDebugLevel )
    REPORTSTATUS ( &status);

  return status.statusCode;
}



  

void PrintLogFile (LALStatus       *status,
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
  strcpy(fnameLog,"./MultiWeight.log");

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


