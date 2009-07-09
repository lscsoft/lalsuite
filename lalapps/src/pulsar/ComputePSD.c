/*
*  Copyright (C) 2007 Badri Krishnan, Iraj Gholami, Reinhard Prix, Karl Wette, Alicia Sintes
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
 * File Name: ComputePSD.c
 * Authors:  Krishnan, B.  ; Sintes, A. M.
 *
 * Revision: $Id$
 *
 * History:   Created by Krishnan 
 *            Modified by Sintes 
 *
 *-----------------------------------------------------------------------
 */
#include <glob.h> 
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <gsl/gsl_sort.h>
#include <gsl/gsl_math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SFTfileIO.h>
#include <lal/Random.h>
#include <lal/PulsarDataTypes.h>
#include <lal/UserInput.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SFTClean.h>

#include <lal/LogPrintf.h>

#include <lalapps.h>

RCSID( "$Id$");

/* ---------- Error codes and messages ---------- */
#define COMPUTEPSDC_ENORM 0
#define COMPUTEPSDC_ESUB  1
#define COMPUTEPSDC_EARG  2
#define COMPUTEPSDC_EBAD  3
#define COMPUTEPSDC_EFILE 4
#define COMPUTEPSDC_ENULL 5
#define COMPUTEPSDC_EMEM  6

#define COMPUTEPSDC_MSGENORM "Normal exit"
#define COMPUTEPSDC_MSGESUB  "Subroutine failed"
#define COMPUTEPSDC_MSGEARG  "Error parsing arguments"
#define COMPUTEPSDC_MSGEBAD  "Bad argument values"
#define COMPUTEPSDC_MSGEFILE "Could not create output file"
#define COMPUTEPSDC_MSGENULL "Null Pointer"
#define COMPUTEPSDC_MSGEMEM  "Out of memory"

/*---------- local defines ---------- */

#define TRUE (1==1)
#define FALSE (1==0)

/* ----- Macros ----- */

/* ---------- local types ---------- */

/** types of mathematical operations */
enum tagMATH_OP_TYPE {
  MATH_OP_ARITHMETIC_SUM = 0,   /** sum(x)     */
  MATH_OP_ARITHMETIC_MEAN,      /** sum(x) / n */
  MATH_OP_ARITHMETIC_MEDIAN,    /** x_1 <= ... x_{n/2} <= .. <= x_n */
  MATH_OP_HARMONIC_SUM,         /** 1 / sum(1/x) */
  MATH_OP_HARMONIC_MEAN,        /** n / sum(1/x) */
  MATH_OP_POWERMINUS2_SUM,        /** 1 / sqrt( sum(1/x/x) )*/
  MATH_OP_POWERMINUS2_MEAN,        /** 1 / sqrt( sum(1/x/x) /n )*/
  MATH_OP_MINIMUM,              /** x_1 <= ... */
  MATH_OP_MAXIMUM,              /** ... <= x_n */
  MATH_OP_LAST
};    

/** user input variables */
typedef struct 
{ 
  BOOLEAN help;

  CHAR *inputData;    /* directory for unclean sfts */
  CHAR *outputPSD;   /* directory for cleaned sfts */
  CHAR *outputSpectBname;
  REAL8 fStart; 
  REAL8 fBand;
  REAL8 startTime; 
  REAL8 endTime;
  CHAR *IFO;
  CHAR  *timeStampsFile;
  LALStringVector *linefiles;
  INT4 blocksRngMed;
  INT4 maxBinsClean;

  INT4 mthopOverSFTs;        /* type of math. operation over SFTs */
  INT4 mthopOverIFOs;        /* type of math. operation over IFOs */
  REAL8 finalPSDBinSize;     /* bin final PSD size */
  INT4 finalPSDBinMthOp;     /* type of math. operation over bins */
  REAL8 finalPSDBinStep;     /* final PSD bin step size */

  CHAR *outputFILE;   /* OBSOLETE: don't use. Superceded by outputPSD */
} UserVariables_t;

/*---------- empty structs for initializations ----------*/
UserVariables_t empty_UserVariables;
/* ---------- global variables ----------*/

extern int vrbflg;
extern INT4 lalDebugLevel;


/* ---------- local prototypes ---------- */
void initUserVars (LALStatus *status, int argc, char *argv[], UserVariables_t *uvar);
void ReadTimeStampsFile (LALStatus *status, LIGOTimeGPSVector  *ts, CHAR *filename);
void LALfwriteSpectrograms ( LALStatus *status, const CHAR *bname, const MultiPSDVector *multiPSD );
static REAL8 math_op(REAL8*, size_t, INT4);

/*============================================================
 * FUNCTION definitions
 *============================================================*/

int 
main(int argc, char *argv[]){ 

  static LALStatus       status;  /* LALStatus pointer */ 
  UserVariables_t uvar = empty_UserVariables;
  
  UINT4 k, numBins, numIFOs, maxNumSFTs, X, alpha;
  REAL8 Freq0, dFreq, normPSD;
  UINT4 finalBinSize, finalBinStep, finalNumBins;

  FILE *fpOut = NULL;

  REAL8Vector *overSFTs = NULL; /* one frequency bin over SFTs */
  REAL8Vector *overIFOs = NULL; /* one frequency bin over IFOs */
  REAL8Vector *finalPSD = NULL; /* math. operation PSD over SFTs and IFOs */

  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_SFTConstraints;

  MultiSFTVector *inputSFTs = NULL;
  MultiPSDVector *multiPSD = NULL;  

  LIGOTimeGPS startTimeGPS, endTimeGPS;
  LIGOTimeGPSVector inputTimeStampsVector;
  
  /* LALDebugLevel must be called before anything else */
  lalDebugLevel = 0;
  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  /* set log-level */
  LogSetLevel ( lalDebugLevel );

  LAL_CALL (initUserVars (&status, argc, argv, &uvar), &status);	  

  /* exit if help was required */
  if (uvar.help)
    exit(0); 

  /* set detector constraint */
  if (LALUserVarWasSet(&uvar.IFO))
    constraints.detector = uvar.IFO;
  else
    constraints.detector = NULL;
  
  if ( LALUserVarWasSet( &uvar.startTime ) ) {
    XLALGPSSetREAL8(&startTimeGPS, uvar.startTime);
    constraints.startTime = &startTimeGPS;
  }
  
  if ( LALUserVarWasSet( &uvar.endTime ) ) {
    XLALGPSSetREAL8(&endTimeGPS, uvar.endTime);
    constraints.endTime = &endTimeGPS;
  }
  
  if ( LALUserVarWasSet( &uvar.timeStampsFile ) ) {
    LAL_CALL ( ReadTimeStampsFile ( &status, &inputTimeStampsVector, uvar.timeStampsFile), &status);
    constraints.timestamps = &inputTimeStampsVector;
  }
  
  /* get sft catalog */
  LogPrintf ( LOG_DEBUG, "Finding all SFTs to load ... ");
  LAL_CALL( LALSFTdataFind( &status, &catalog, uvar.inputData, &constraints), &status);
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");
  if ( (catalog == NULL) || (catalog->length == 0) ) {
    fprintf (stderr,"Unable to match any SFTs with pattern '%s'\n", uvar.inputData );
    exit(1);
  }
  
  /* now we can free the inputTimeStampsVector */
  if ( LALUserVarWasSet( &uvar.timeStampsFile ) ) {
    LALFree( inputTimeStampsVector.data );
  }

  /* read the sfts */  
  LogPrintf (LOG_DEBUG, "Loading all SFTs ... ");
  LAL_CALL( LALLoadMultiSFTs ( &status, &inputSFTs, catalog, uvar.fStart, uvar.fStart + uvar.fBand), &status);
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

  LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);  	

  /* clean sfts if required */
  if ( LALUserVarWasSet( &uvar.linefiles ) )
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
      
      LAL_CALL( LALRemoveKnownLinesInMultiSFTVector ( &status, inputSFTs, uvar.maxBinsClean, uvar.blocksRngMed, uvar.linefiles, randPar), &status);
      LAL_CALL ( LALDestroyRandomParams (&status, &randPar), &status);
      fclose(fpRand);
    } /* end cleaning */
  
  LogPrintf (LOG_DEBUG, "Computing spectrogram and PSD ... ");
  /* get power running-median rngmed[ |data|^2 ] from SFTs */
  LAL_CALL( LALNormalizeMultiSFTVect (&status, &multiPSD, inputSFTs, uvar.blocksRngMed), &status);
  /* Throw away SFTs, not needed any more */
  LAL_CALL (LALDestroyMultiSFTVector(&status, &inputSFTs), &status );

  /* start frequency and frequency spacing */
  Freq0 = multiPSD->data[0]->data[0].f0;
  dFreq = multiPSD->data[0]->data[0].deltaF;

  /* number of raw bins in final PSD */
  numBins = multiPSD->data[0]->data[0].data->length;
  if ( (finalPSD = XLALCreateREAL8Vector ( numBins )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Out of memory!\n");
    exit (-1);
  }

  /* number of IFOs */
  numIFOs = multiPSD->length;
  if ( (overIFOs = XLALCreateREAL8Vector ( numIFOs )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Out of memory!\n");
    exit (-1);
  }

  /* maximum number of SFTs */
  maxNumSFTs = 0;
  for (X = 0; X < numIFOs; ++X) {
    maxNumSFTs = GSL_MAX(maxNumSFTs, multiPSD->data[X]->length);    
  }
  if ( (overSFTs = XLALCreateREAL8Vector ( maxNumSFTs )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Out of memory!\n");
    exit (-1);
  }

  /* normalize rngmd(power) to get proper *single-sided* PSD: Sn = (2/Tsft) rngmed[|data|^2]] */
  normPSD = 2.0 * dFreq;

  /* loop over frequency bins in final PSD */
  for (k = 0; k < numBins; ++k) {
 
    /* loop over IFOs */
    for (X = 0; X < numIFOs; ++X) {

      /* number of SFTs for this IFO */
      UINT4 numSFTs = multiPSD->data[X]->length;

      /* copy PSD frequency bins and normalise multiPSD for later use */
      for (alpha = 0; alpha < numSFTs; ++alpha) {
	multiPSD->data[X]->data[alpha].data->data[k] *= normPSD;
	overSFTs->data[alpha] = multiPSD->data[X]->data[alpha].data->data[k];
      }

      /* compute math. operation over SFTs for this IFO */
      overIFOs->data[X] = math_op(overSFTs->data, numSFTs, uvar.mthopOverSFTs);

    } /* over IFOs */

    /* compute math. operation over IFOs for this frequency */
    finalPSD->data[k] = math_op(overIFOs->data, numIFOs, uvar.mthopOverIFOs);

  } /* over freq bins */
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

  /* output spectrograms */
  if ( uvar.outputSpectBname ) {
    LAL_CALL ( LALfwriteSpectrograms ( &status, uvar.outputSpectBname, multiPSD ), &status );
  }

  /* open output file */
  if (  (fpOut = fopen(uvar.outputPSD, "wb")) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Unable to open output file %s for writing...exiting \n", uvar.outputPSD );
    exit(1);
  }

  /* work out bin size */
  if (uvar.finalPSDBinSize <= 0.0) {
    finalBinSize = 1;
  }
  else {
    finalBinSize = (UINT4)floor(uvar.finalPSDBinSize / dFreq + 0.5); /* round to nearest bin */
  }

  /* work out bin step */
  if (uvar.finalPSDBinStep <= 0.0) {
    finalBinStep = finalBinSize;
  }
  else {
    finalBinStep = (UINT4)floor(uvar.finalPSDBinStep / dFreq + 0.5); /* round to nearest bin */
  }

  /* work out total number of bins */
  finalNumBins = (UINT4)floor((numBins - finalBinSize) / finalBinStep) + 1;

  /* write final PSD to file */
  LogPrintf(LOG_DEBUG, "Printing PSD to file ... ");
  for (k = 0; k < finalNumBins; ++k) {
    UINT4 b = k * finalBinStep;
    REAL8 f = Freq0 + b * dFreq;
    REAL8 p = math_op(&(finalPSD->data[b]), finalBinSize, uvar.finalPSDBinMthOp);
    fprintf(fpOut, "%f   %e\n", f, p);
  }
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

  fclose(fpOut);

  /* we are now done with the psd */
  LAL_CALL ( LALDestroyMultiPSDVector  ( &status, &multiPSD), &status);
    
  LAL_CALL (LALDestroyUserVars(&status), &status);

  XLALDestroyREAL8Vector ( overSFTs );
  XLALDestroyREAL8Vector ( overIFOs );
  XLALDestroyREAL8Vector ( finalPSD );

  LALCheckMemoryLeaks(); 

  return COMPUTEPSDC_ENORM;

} /* main() */



/* read timestamps file */
void ReadTimeStampsFile (LALStatus          *status,
			 LIGOTimeGPSVector  *ts,
			 CHAR               *filename)
{

  FILE  *fp = NULL;
  INT4  numTimeStamps, r;
  UINT4 j;
  REAL8 temp1, temp2;

  INITSTATUS (status, "ReadTimeStampsFile", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT(ts, status, COMPUTEPSDC_ENULL,COMPUTEPSDC_MSGENULL); 
  ASSERT(ts->data == NULL, status, COMPUTEPSDC_ENULL,COMPUTEPSDC_MSGENULL); 
  ASSERT(ts->length == 0, status, COMPUTEPSDC_ENULL,COMPUTEPSDC_MSGENULL); 
  ASSERT(filename, status, COMPUTEPSDC_ENULL,COMPUTEPSDC_MSGENULL); 

  if ( (fp = fopen(filename, "r")) == NULL) {
    ABORT( status, COMPUTEPSDC_EFILE, COMPUTEPSDC_MSGEFILE);
  }

  /* count number of timestamps */
  numTimeStamps = 0;     

  do {
    r = fscanf(fp,"%lf%lf\n", &temp1, &temp2);
    /* make sure the line has the right number of entries or is EOF */
    if (r==2) numTimeStamps++;
  } while ( r != EOF);
  rewind(fp);

  ts->length = numTimeStamps;
  ts->data = LALCalloc (1, numTimeStamps * sizeof(LIGOTimeGPS));;
  if ( ts->data == NULL ) {
    fclose(fp);
    ABORT( status, COMPUTEPSDC_ENULL, COMPUTEPSDC_MSGENULL);
  }
  
  for (j = 0; j < ts->length; j++)
    {
      r = fscanf(fp,"%lf%lf\n", &temp1, &temp2);
      ts->data[j].gpsSeconds = (INT4)temp1;
      ts->data[j].gpsNanoSeconds = (INT4)temp2;
    }
  
  fclose(fp);
  	 
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}    

/** compute the various kinds of math. operation */
REAL8 math_op(REAL8* data, size_t length, INT4 type) {

  UINT4 i;
  REAL8 res = 0.0;

  switch (type) {

  case MATH_OP_ARITHMETIC_SUM: /* sum(data) */

    for (i = 0; i < length; ++i) res += *(data++);

    break;

  case MATH_OP_ARITHMETIC_MEAN: /* sum(data)/length  */

    for (i = 0; i < length; ++i) res += *(data++);
    res /= (REAL8)length;

    break;

  case MATH_OP_ARITHMETIC_MEDIAN: /* middle element of sort(data) */

    gsl_sort(data, 1, length);
    if (length/2 == (length+1)/2) /* length is even */ {
      res = (data[length/2] + data[length/2+1])/2;
    }
    else /* length is odd */ {
      res = data[length/2+1];
    }

    break;

  case MATH_OP_HARMONIC_SUM: /* 1 / sum(1 / data) */

    for (i = 0; i < length; ++i) res += 1.0 / *(data++);
    res = 1.0 / res;

    break;

  case MATH_OP_HARMONIC_MEAN: /* length / sum(1 / data) */

    for (i = 0; i < length; ++i) res += 1.0 / *(data++);
    res = (REAL8)length / res;

    break;

  case MATH_OP_POWERMINUS2_SUM: /*   1 / sqrt ( sum(1 / data/data) )*/

    for (i = 0; i < length; ++i) res += 1.0 / (data[i]*data[i]);
    res = 1.0 / sqrt(res);

    break;
 
   case MATH_OP_POWERMINUS2_MEAN: /*   1 / sqrt ( sum(1/data/data) / length )*/

    for (i = 0; i < length; ++i) res += 1.0 / (data[i]*data[i]);
    res = 1.0 / sqrt(res / (REAL8)length);

    break;
   
  case MATH_OP_MINIMUM: /* first element of sort(data) */

    gsl_sort(data, 1, length);
    res = data[0];
    break;

  case MATH_OP_MAXIMUM: /* first element of sort(data) */

    gsl_sort(data, 1, length);
    res = data[length-1];
    break;

  default:

    LALPrintError("'%s' is not a valid math. operation", type);
    exit (-1);

  }

  return res;

}


/** register all "user-variables" */
void
initUserVars (LALStatus *status, int argc, char *argv[], UserVariables_t *uvar)
{
  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
  uvar->help = FALSE;

  uvar->maxBinsClean = 100;
  uvar->blocksRngMed = 101;

  uvar->startTime = 0.0;
  uvar->endTime = 0.0;

  uvar->inputData = NULL;

  uvar->IFO = NULL;

  /* default: read all SFT bins */
  uvar->fStart = -1;
  uvar->fBand = 0;

#define OUTPUTPSDFILE "./psd"
  uvar->outputPSD = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar->outputPSD, OUTPUTPSDFILE);

  uvar->mthopOverSFTs = MATH_OP_HARMONIC_MEAN;
  uvar->mthopOverIFOs = MATH_OP_HARMONIC_SUM;
  uvar->finalPSDBinSize = -1.0;
  uvar->finalPSDBinMthOp = MATH_OP_ARITHMETIC_MEDIAN;
  uvar->finalPSDBinStep = -1.0;

  /* register user input variables */
  LALregBOOLUserStruct ( status, 	help, 		'h', UVAR_HELP,    	"Print this message" );
  LALregSTRINGUserStruct ( status, 	inputData, 	'i', UVAR_REQUIRED, 	"Input SFT pattern");
  LALregSTRINGUserStruct ( status, 	outputPSD, 	'o', UVAR_OPTIONAL, 	"Output PSD into this file");
  LALregSTRINGUserStruct ( status,     outputSpectBname, 0,  UVAR_OPTIONAL, 	"Filename-base for (binary) spectrograms (one per IFO)");

  LALregREALUserStruct ( status, 	fStart, 	'f', UVAR_OPTIONAL, 	"Frequency to start from (-1 = all freqs)");
  LALregREALUserStruct ( status, 	fBand, 		'b', UVAR_OPTIONAL, 	"Frequency Band");
  LALregREALUserStruct ( status, 	startTime, 	's', UVAR_OPTIONAL, 	"GPS start time");
  LALregREALUserStruct ( status, 	endTime,  	'e', UVAR_OPTIONAL, 	"GPS end time");
  LALregSTRINGUserStruct ( status, 	timeStampsFile, 't', UVAR_OPTIONAL, 	"Time-stamps file");
  LALregSTRINGUserStruct ( status,      IFO,             0 , UVAR_OPTIONAL,     "Detector filter");

  LALregINTUserStruct ( status, 	blocksRngMed,  	'w', UVAR_OPTIONAL, 	"Running Median window size");
  LALregINTUserStruct ( status,         mthopOverSFTs,  'S', UVAR_OPTIONAL,     "Type of math. operation over SFTs:  0=arith-sum, 1=arith-mean, 2=arith-median, 3=harm-sum, 4=harm-mean, 5=power-2-sum, 6=power-2-mean, 7=min, 8=max");
  LALregINTUserStruct ( status,         mthopOverIFOs,  'I', UVAR_OPTIONAL,     "Type of math. operation over IFOs: as for mthopOverSFTs");
  LALregREALUserStruct( status,         finalPSDBinSize,'B', UVAR_OPTIONAL,     "Bin the final PSD into bins of size (default: no binning)");
  LALregINTUserStruct ( status,        finalPSDBinMthOp,'A', UVAR_OPTIONAL,     "If binning, type of math. operation over frequency bins: as for mthopOverSFTs");
  LALregREALUserStruct( status,         finalPSDBinStep,'D', UVAR_OPTIONAL,     "If binning, step size to move bin along (default: size of bin)");

  LALregINTUserStruct ( status, 	maxBinsClean, 	'm', UVAR_OPTIONAL, 	"Maximum Cleaning Bins");
  LALregLISTUserStruct ( status, 	linefiles, 	 0, UVAR_OPTIONAL, 	"Comma separated list of linefiles (names must contain IFO name)");

  LALregSTRINGUserStruct ( status, 	outputFILE, 	 0, UVAR_DEVELOPER, 	"Output PSD file [OBSOLETE: use '--outputPSD' instead]");

  /* read all command line variables */
  TRY( LALUserVarReadAllInput(status->statusPtr, argc, argv), status);

  /* check user-input consistency */
  if ( LALUserVarWasSet ( &uvar->outputFILE ) )
    {
      LogPrintf (LOG_NORMAL, "Warning: --outputFILE is obsolete, use --outputPSD instead!\n");
      if ( LALUserVarWasSet ( &uvar->outputPSD ) ) {
	LogPrintf ( LOG_CRITICAL, "Both --outputFILE and --outputPSD was specified!\n" );
	ABORT ( status, COMPUTEPSDC_EBAD, COMPUTEPSDC_MSGEBAD );
      }
      if ( uvar->outputPSD ) LALFree ( uvar->outputPSD );
      if ( (uvar->outputPSD = LALMalloc ( strlen(uvar->outputFILE ) + 1)) == NULL ) {
	ABORT ( status, COMPUTEPSDC_EMEM, COMPUTEPSDC_MSGEMEM );
      }
      strcpy ( uvar->outputPSD, uvar->outputFILE );
    }

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */


/** Write a multi-PSD into spectrograms for each IFO.
 * Using gnuplot 'binary' matrix format
 * The filename for each <IFO> is generated as 'bname-<IFO>'
 */
void
LALfwriteSpectrograms ( LALStatus *status, const CHAR* bname, const MultiPSDVector *multiPSD )
{
  UINT4 X;
  CHAR *fname;
  float num, *row_data;		/* cast to float for writing (gnuplot binary format) */
  FILE *fp;
      
  INITSTATUS( status, "LALfwriteSpectrograms", rcsid );
  ATTATCHSTATUSPTR (status);

  if ( !bname || !multiPSD || multiPSD->length == 0 ) {
    ABORT ( status, COMPUTEPSDC_ENULL, COMPUTEPSDC_MSGENULL );
  }

  /* loop over IFOs */
  for ( X = 0; X < multiPSD->length ; X ++ )
    {
      UINT4 len = strlen ( bname ) + 4;	/* append '-XN' to get IFO-specific filename */
      UINT4 numSFTs, numBins;
      UINT4 j, k;
      const CHAR *tmp;
      REAL8 f0, df;

      numSFTs = multiPSD->data[X]->length;
      numBins = multiPSD->data[X]->data[0].data->length;

      /* allocate memory for data row-vector */
      if ( ( row_data = LALMalloc ( numBins * sizeof(float) )) == NULL ) {
	ABORT ( status, COMPUTEPSDC_EMEM, COMPUTEPSDC_MSGEMEM );
      }

      if ( ( fname = LALMalloc ( len * sizeof(CHAR) )) == NULL ) {
	LALFree ( row_data );
	ABORT ( status, COMPUTEPSDC_EMEM, COMPUTEPSDC_MSGEMEM );
      }
      tmp = multiPSD->data[X]->data[0].name;
      sprintf ( fname, "%s-%c%c", bname, tmp[0], tmp[1] );

      if ( ( fp = fopen( fname, "wb" ))  == NULL ) {
	LogPrintf (LOG_CRITICAL, "Failed to open spectrogram file '%s' for writing!\n", fname );
	goto failed;

      }

      /* write number of columns: i.e. frequency-bins */
      num = (float)numBins;
      if ((fwrite((char *) &num, sizeof(float), 1, fp)) != 1) {
	LogPrintf (LOG_CRITICAL, "Failed to fwrite() to spectrogram file '%s'\n", fname );
	goto failed;
      }

      /* write frequencies as column-titles */
      f0 = multiPSD->data[X]->data[0].f0;
      df = multiPSD->data[X]->data[0].deltaF;
      for ( k=0; k < numBins; k ++ )
	row_data[k] = (float) ( f0 + 1.0 * k * df );
      if ( fwrite((char *) row_data, sizeof(float), numBins, fp) != numBins ) {
	LogPrintf (LOG_CRITICAL, "Failed to fwrite() to spectrogram file '%s'\n", fname );
	goto failed;
      }

      /* write PSDs of successive SFTs in rows, first column is GPS-time in seconds */
      for ( j = 0; j < numSFTs ; j ++ )
	{
	  num = (float) multiPSD->data[X]->data[j].epoch.gpsSeconds;
	  for ( k = 0; k < numBins; k ++ )
	    row_data[k] = (float) sqrt ( multiPSD->data[X]->data[j].data->data[k] );

	  if ( ( fwrite((char *) &num, sizeof(float), 1, fp) != 1 ) ||
	       ( fwrite((char *) row_data, sizeof(float), numBins, fp) != numBins ) ) {
	    LogPrintf (LOG_CRITICAL, "Failed to fwrite() to spectrogram file '%s'\n", fname );
	    goto failed;
	  }
	  
	} /* for j < numSFTs */

      fclose ( fp );
      LALFree ( fname );
      LALFree ( row_data );

    } /* for X < numIFOs */

  DETATCHSTATUSPTR (status);
  RETURN (status);
  
  /* cleanup and exit on write-error */
 failed:
  if ( fname ) LALFree ( fname );
  if ( row_data ) LALFree ( row_data );
  if ( fp ) fclose ( fp );
  ABORT ( status, COMPUTEPSDC_EFILE, COMPUTEPSDC_MSGEFILE );

} /* LALfwriteSpectrograms() */
