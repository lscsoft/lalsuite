/*
*  Copyright (C) 2007, 2010 Karl Wette
*  Copyright (C) 2007 Badri Krishnan, Iraj Gholami, Reinhard Prix, Alicia Sintes
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
 * \file
 * \ingroup pulsarApps
 * \author Badri Krishnan, Iraj Gholami, Reinhard Prix, Alicia Sintes, Karl Wette
 * \brief Compute power spectral densities
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
#include <lal/Segments.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>

#include <lal/LogPrintf.h>

#include <lalapps.h>

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
  MATH_OP_ARITHMETIC_SUM = 0,   /**< sum(x)     */
  MATH_OP_ARITHMETIC_MEAN,      /**< sum(x) / n */
  MATH_OP_ARITHMETIC_MEDIAN,    /**< x_1 <= ... x_{n/2} <= .. <= x_n */
  MATH_OP_HARMONIC_SUM,         /**< 1 / sum(1/x) */
  MATH_OP_HARMONIC_MEAN,        /**< n / sum(1/x) */
  MATH_OP_POWERMINUS2_SUM,      /**< 1 / sqrt( sum(1/x/x) )*/
  MATH_OP_POWERMINUS2_MEAN,     /**< 1 / sqrt( sum(1/x/x) /n )*/
  MATH_OP_MINIMUM,              /**< x_1 <= ... */
  MATH_OP_MAXIMUM,              /**< ... <= x_n */
  MATH_OP_LAST
};

/** user input variables */
typedef struct
{
  BOOLEAN help;

  CHAR *inputData;    	/**< directory for input sfts */
  CHAR *outputPSD;    	/**< directory for output sfts */
  CHAR *outputSpectBname;

  REAL8 Freq;		/**< *physical* start frequency to compute PSD for (excluding rngmed wings) */
  REAL8 FreqBand;	/**< *physical* frequency band to compute PSD for (excluding rngmed wings) */

  REAL8 startTime;	/**< earliest SFT-timestamp to include */
  REAL8 endTime;	/**< last SFT-timestamp to include */
  CHAR *IFO;
  CHAR  *timeStampsFile;
  LALStringVector *linefiles;
  INT4 blocksRngMed;	/**< number of running-median bins to use */
  INT4 maxBinsClean;

  INT4 PSDmthopSFTs;     /**< for PSD, type of math. operation over SFTs */
  INT4 PSDmthopIFOs;     /**< for PSD, type of math. operation over IFOs */
  BOOLEAN outputNormSFT; /**< output normalised SFT power? */
  INT4 nSFTmthopSFTs;    /**< for norm. SFT, type of math. operation over SFTs */
  INT4 nSFTmthopIFOs;    /**< for norm. SFT, type of math. operation over IFOs */

  REAL8 binSizeHz;       /**< output PSD bin size in Hz */
  INT4  binSize;         /**< output PSD bin size in no. of bins */
  INT4  PSDmthopBins;    /**< for PSD, type of math. operation over bins */
  INT4  nSFTmthopBins;   /**< for norm. SFT, type of math. operation over bins */
  REAL8 binStepHz;       /**< output PSD bin step in Hz */
  INT4  binStep;         /**< output PSD bin step in no. of bins */
  BOOLEAN outFreqBinEnd; /**< output the end frequency of each bin? */

  BOOLEAN dumpMultiPSDVector; /**< output multi-PSD vector over IFOs, timestamps, and frequencies into file(s) */
  CHAR *outputQ;	/**< output the 'data-quality factor' Q(f) into this file */

  REAL8 fStart;		/**< Start Frequency to load from SFT and compute PSD, including wings (it is RECOMMENDED to use --Freq instead) */
  REAL8 fBand;		/**< Frequency Band to load from SFT and compute PSD, including wings (it is RECOMMENDED to use --FreqBand instead) */

  BOOLEAN version;	/**< Output version information */

} UserVariables_t;

/**
 * Config variables 'derived' from user-input
 */
typedef struct
{
  UINT4 firstBin;		/**< first PSD bin for output */
  UINT4 lastBin;		/**< last PSD bin for output */
  LALSeg dataSegment;		/**< the data-segment for which PSD was computed */
} ConfigVariables_t;


/*---------- empty structs for initializations ----------*/
UserVariables_t empty_UserVariables;
ConfigVariables_t empty_ConfigVariables;
/* ---------- global variables ----------*/

extern int vrbflg;


/* ---------- local prototypes ---------- */
int initUserVars (int argc, char *argv[], UserVariables_t *uvar);
void LALfwriteSpectrograms ( LALStatus *status, const CHAR *bname, const MultiPSDVector *multiPSD );
static REAL8 math_op(REAL8*, size_t, INT4);
int XLALDumpMultiPSDVector ( const CHAR *outbname, const MultiPSDVector *multiPSDVect );
MultiSFTVector *XLALReadSFTs ( ConfigVariables_t *cfg, const UserVariables_t *uvar );

int XLALCropMultiPSDandSFTVectors ( MultiPSDVector *multiPSDVect, MultiSFTVector *multiSFTVect, UINT4 firstBin, UINT4 lastBin );
REAL8FrequencySeries *XLALComputeSegmentDataQ ( const MultiPSDVector *multiPSDVect, LALSeg segment );
int XLALWriteREAL8FrequencySeries_to_file ( const REAL8FrequencySeries *series, const char *fname );

/*============================================================
 * FUNCTION definitions
 *============================================================*/
int
main(int argc, char *argv[])
{
  static LALStatus       status;  /* LALStatus pointer */
  UserVariables_t uvar = empty_UserVariables;
  ConfigVariables_t cfg = empty_ConfigVariables;

  UINT4 k, numBins, numIFOs, maxNumSFTs, X, alpha;
  REAL8 Freq0, dFreq, normPSD;
  UINT4 finalBinSize, finalBinStep, finalNumBins;

  REAL8Vector *overSFTs = NULL; /* one frequency bin over SFTs */
  REAL8Vector *overIFOs = NULL; /* one frequency bin over IFOs */
  REAL8Vector *finalPSD = NULL; /* math. operation PSD over SFTs and IFOs */
  REAL8Vector *finalNormSFT = NULL; /* normalised SFT power */

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* set log-level */
  LogSetLevel ( lalDebugLevel );

  /* register and read user variables */
  if (initUserVars(argc, argv, &uvar) != XLAL_SUCCESS)
    return EXIT_FAILURE;

  /* assemble version string */
  CHAR *VCSInfoString;
  if ( (VCSInfoString = XLALGetVersionString(0)) == NULL ) {
    XLALPrintError("XLALGetVersionString(0) failed with xlalErrno = %d\n", xlalErrno );
    return EXIT_FAILURE;
  }

  if ( uvar.version )
    {
      printf ("%s\n", VCSInfoString );
      return (0);
    }

  /* exit if help was required */
  if (uvar.help)
    return EXIT_SUCCESS;

  MultiSFTVector *inputSFTs = NULL;
  if ( ( inputSFTs = XLALReadSFTs ( &cfg, &uvar ) ) == NULL )
    {
      XLALPrintError ("Call to XLALReadSFTs() failed with xlalErrno = %d\n", xlalErrno );
      return EXIT_FAILURE;
    }

  /* clean sfts if required */
  if ( XLALUserVarWasSet( &uvar.linefiles ) )
    {
      RandomParams *randPar=NULL;
      FILE *fpRand=NULL;
      INT4 seed, ranCount;

      if ( (fpRand = fopen("/dev/urandom", "r")) == NULL ) {
	fprintf(stderr,"Error in opening /dev/urandom" );
	return EXIT_FAILURE;
      }

      if ( (ranCount = fread(&seed, sizeof(seed), 1, fpRand)) != 1 ) {
	fprintf(stderr,"Error in getting random seed" );
	return EXIT_FAILURE;
      }

      LAL_CALL ( LALCreateRandomParams (&status, &randPar, seed), &status );

      LAL_CALL( LALRemoveKnownLinesInMultiSFTVector ( &status, inputSFTs, uvar.maxBinsClean, uvar.blocksRngMed, uvar.linefiles, randPar), &status);
      LAL_CALL ( LALDestroyRandomParams (&status, &randPar), &status);
      fclose(fpRand);
    } /* end cleaning */

  LogPrintf (LOG_DEBUG, "Computing spectrogram and PSD ... ");

  /* get power running-median rngmed[ |data|^2 ] from SFTs */
  MultiPSDVector *multiPSD = NULL;
  LAL_CALL( LALNormalizeMultiSFTVect (&status, &multiPSD, inputSFTs, uvar.blocksRngMed), &status);
  /* restrict this PSD to just the "physical" band if requested using {--Freq, --FreqBand} */
  if ( ( XLALCropMultiPSDandSFTVectors ( multiPSD, inputSFTs, cfg.firstBin, cfg.lastBin )) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALCropMultiPSDandSFTVectors (inputPSD, inputSFTs, %d, %d) failed with xlalErrno = %d\n", __func__, cfg.firstBin, cfg.lastBin, xlalErrno );
    return EXIT_FAILURE;
  }

  /* start frequency and frequency spacing */
  Freq0 = multiPSD->data[0]->data[0].f0;
  dFreq = multiPSD->data[0]->data[0].deltaF;

  /* number of raw bins in final PSD */
  numBins = multiPSD->data[0]->data[0].data->length;
  if ( (finalPSD = XLALCreateREAL8Vector ( numBins )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Out of memory!\n");
    return EXIT_FAILURE;
  }

  /* number of IFOs */
  numIFOs = multiPSD->length;
  if ( (overIFOs = XLALCreateREAL8Vector ( numIFOs )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Out of memory!\n");
    return EXIT_FAILURE;
  }

  /* maximum number of SFTs */
  maxNumSFTs = 0;
  for (X = 0; X < numIFOs; ++X) {
    maxNumSFTs = GSL_MAX(maxNumSFTs, multiPSD->data[X]->length);
  }
  if ( (overSFTs = XLALCreateREAL8Vector ( maxNumSFTs )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Out of memory!\n");
    return EXIT_FAILURE;
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
      overIFOs->data[X] = math_op(overSFTs->data, numSFTs, uvar.PSDmthopSFTs);
      if ( isnan( overIFOs->data[X] ) )
        XLAL_ERROR ( EXIT_FAILURE, "Found Not-A-Number in overIFOs->data[X=%d] = NAN ... exiting\n", X );

    } /* for IFOs X */

    /* compute math. operation over IFOs for this frequency */
    finalPSD->data[k] = math_op(overIFOs->data, numIFOs, uvar.PSDmthopIFOs);
    if ( isnan ( finalPSD->data[k] ) )
      XLAL_ERROR ( EXIT_FAILURE, "Found Not-A-Number in finalPSD->data[k=%d] = NAN ... exiting\n", k );

  } /* for freq bins k */
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

  /* compute normalised SFT power */
  if (uvar.outputNormSFT) {
    LogPrintf (LOG_DEBUG, "Computing normalised SFT power ... ");

    if ( (finalNormSFT = XLALCreateREAL8Vector ( numBins )) == NULL ) {
      LogPrintf (LOG_CRITICAL, "Out of memory!\n");
      return EXIT_FAILURE;
    }

    /* loop over frequency bins in SFTs */
    for (k = 0; k < numBins; ++k) {

      /* loop over IFOs */
      for (X = 0; X < numIFOs; ++X) {

	/* number of SFTs for this IFO */
	UINT4 numSFTs = inputSFTs->data[X]->length;

	/* compute SFT power */
	for (alpha = 0; alpha < numSFTs; ++alpha) {
	  COMPLEX8 bin = inputSFTs->data[X]->data[alpha].data->data[k];
	  overSFTs->data[alpha] = crealf(bin)*crealf(bin) + cimagf(bin)*cimagf(bin);
	}

	/* compute math. operation over SFTs for this IFO */
	overIFOs->data[X] = math_op(overSFTs->data, numSFTs, uvar.nSFTmthopSFTs);
	if ( isnan ( overIFOs->data[X] ))
          XLAL_ERROR ( EXIT_FAILURE, "Found Not-A-Number in overIFOs->data[X=%d] = NAN ... exiting\n", X );

      } /* over IFOs */

      /* compute math. operation over IFOs for this frequency */
      finalNormSFT->data[k] = math_op(overIFOs->data, numIFOs, uvar.nSFTmthopIFOs);
      if ( isnan( finalNormSFT->data[k] ) )
        XLAL_ERROR ( EXIT_FAILURE, "Found Not-A-Number in bin finalNormSFT->data[k=%d] = NAN ... exiting\n", k );

    } /* over freq bins */
    LogPrintfVerbatim ( LOG_DEBUG, "done.\n");
  }

  /* output spectrograms */
  if ( uvar.outputSpectBname ) {
    LAL_CALL ( LALfwriteSpectrograms ( &status, uvar.outputSpectBname, multiPSD ), &status );
  }

  /* ---------- if user requested it, output complete MultiPSDVector over IFOs X, timestamps and freq-bins into ASCI file(s) */
  if ( uvar.dumpMultiPSDVector ) {
    if ( XLALDumpMultiPSDVector ( uvar.outputPSD, multiPSD ) != XLAL_SUCCESS ) {
      XLALPrintError ("%s: XLALDumpMultiPSDVector() failed, xlalErrnor = %d\n", __func__, xlalErrno );
      return EXIT_FAILURE;
    }
  } /* if uvar.dumpMultiPSDVector */

  /* ----- if requested, compute data-quality factor 'Q' -------------------- */
  if ( uvar.outputQ )
    {
      REAL8FrequencySeries *Q;
      if ( (Q = XLALComputeSegmentDataQ ( multiPSD, cfg.dataSegment )) == NULL ) {
        XLALPrintError ("%s: XLALComputeSegmentDataQ() failed with xlalErrno = %d\n", __func__, xlalErrno );
        return EXIT_FAILURE;
      }
      if ( XLAL_SUCCESS != XLALWriteREAL8FrequencySeries_to_file ( Q, uvar.outputQ ) ) {
        return EXIT_FAILURE;
      }
      XLALDestroyREAL8FrequencySeries ( Q );
    } /* if outputQ */

  /* ---------- BINNING if requested ---------- */
  /* work out bin size */
  if (XLALUserVarWasSet(&uvar.binSize)) {
    finalBinSize = uvar.binSize;
  }
  else if (XLALUserVarWasSet(&uvar.binSizeHz)) {
    finalBinSize = (UINT4)floor(uvar.binSizeHz / dFreq + 0.5); /* round to nearest bin */
  }
  else {
    finalBinSize = 1;
  }

  /* work out bin step */
  if (XLALUserVarWasSet(&uvar.binStep)) {
    finalBinStep = uvar.binStep;
  }
  else if (XLALUserVarWasSet(&uvar.binStepHz)) {
    finalBinStep = (UINT4)floor(uvar.binStepHz / dFreq + 0.5); /* round to nearest bin */
  }
  else {
    finalBinStep = finalBinSize;
  }

  /* work out total number of bins */
  finalNumBins = (UINT4)floor((numBins - finalBinSize) / finalBinStep) + 1;

  /* write final PSD to file */
  if (XLALUserVarWasSet(&uvar.outputPSD)) {

    FILE *fpOut = NULL;

    if ((fpOut = fopen(uvar.outputPSD, "wb")) == NULL) {
      LogPrintf ( LOG_CRITICAL, "Unable to open output file %s for writing...exiting \n", uvar.outputPSD );
      return EXIT_FAILURE;
    }

    /* write header info in comments */
    if ( XLAL_SUCCESS != XLALOutputVersionString ( fpOut, 0 ) )
      XLAL_ERROR ( XLAL_EFUNC );

    /* write the command-line */
    for (int a = 0; a < argc; a++)
      fprintf(fpOut,"%%%% argv[%d]: '%s'\n", a, argv[a]);

    /* write column headings */
    fprintf(fpOut,"%%%% columns:\n%%%% FreqBinStart");
    if (uvar.outFreqBinEnd)
      fprintf(fpOut," FreqBinEnd");
    fprintf(fpOut," PSD");
    if (uvar.outputNormSFT)
      fprintf(fpOut," normSFTpower");
    fprintf(fpOut,"\n");

    LogPrintf(LOG_DEBUG, "Printing PSD to file ... ");
    for (k = 0; k < finalNumBins; ++k) {
      UINT4 b = k * finalBinStep;

      REAL8 f0 = Freq0 + b * dFreq;
      REAL8 f1 = f0 + finalBinStep * dFreq;
      fprintf(fpOut, "%f", f0);
      if (uvar.outFreqBinEnd)
	fprintf(fpOut, "   %f", f1);

      REAL8 psd = math_op(&(finalPSD->data[b]), finalBinSize, uvar.PSDmthopBins);
      if ( isnan ( psd ))
        XLAL_ERROR ( EXIT_FAILURE, "Found Not-A-Number in psd[k=%d] = NAN ... exiting\n", k );

      fprintf(fpOut, "   %e", psd);

      if (uvar.outputNormSFT) {
	REAL8 nsft = math_op(&(finalNormSFT->data[b]), finalBinSize, uvar.nSFTmthopBins);
	if ( isnan ( nsft ))
          XLAL_ERROR ( EXIT_FAILURE, "Found Not-A-Number in nsft[k=%d] = NAN ... exiting\n", k );

	fprintf(fpOut, "   %f", nsft);
      }

      fprintf(fpOut, "\n");
    } // k < finalNumBins
    LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

    fclose(fpOut);

  }

  /* we are now done with the psd */
  LAL_CALL ( LALDestroyMultiPSDVector  ( &status, &multiPSD), &status);
  LAL_CALL ( LALDestroyMultiSFTVector  (&status, &inputSFTs), &status);

  LAL_CALL (LALDestroyUserVars(&status), &status);

  XLALDestroyREAL8Vector ( overSFTs );
  XLALDestroyREAL8Vector ( overIFOs );
  XLALDestroyREAL8Vector ( finalPSD );
  XLALDestroyREAL8Vector ( finalNormSFT );
  XLALFree ( VCSInfoString );

  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

} /* main() */

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
      res = data[length/2];
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

    XLALPrintError("'%i' is not a valid math. operation", type);
    XLAL_ERROR_REAL8(XLAL_EINVAL);

  }

  return res;

}


/** register all "user-variables" */
int
initUserVars (int argc, char *argv[], UserVariables_t *uvar)
{

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

  uvar->outputPSD = NULL;
  uvar->outputNormSFT = FALSE;
  uvar->outFreqBinEnd = FALSE;

  uvar->PSDmthopSFTs = MATH_OP_HARMONIC_MEAN;
  uvar->PSDmthopIFOs = MATH_OP_HARMONIC_SUM;

  uvar->nSFTmthopSFTs = MATH_OP_ARITHMETIC_MEAN;
  uvar->nSFTmthopIFOs = MATH_OP_MAXIMUM;
  uvar->dumpMultiPSDVector = FALSE;

  uvar->binSizeHz = 0.0;
  uvar->binSize   = 1;
  uvar->PSDmthopBins  = MATH_OP_ARITHMETIC_MEDIAN;
  uvar->nSFTmthopBins = MATH_OP_MAXIMUM;
  uvar->binStep   = 0.0;
  uvar->binStep   = 1;

  uvar->version = FALSE;

  /* register user input variables */
  XLALregBOOLUserStruct  (help,             'h', UVAR_HELP,     "Print this message" );
  XLALregSTRINGUserStruct(inputData,        'i', UVAR_REQUIRED, "Input SFT pattern");
  XLALregSTRINGUserStruct(outputPSD,        'o', UVAR_OPTIONAL, "Output PSD into this file");
  XLALregSTRINGUserStruct(outputQ,	     0,  UVAR_OPTIONAL, "Output the 'data-quality factor' Q(f) into this file");
  XLALregSTRINGUserStruct(outputSpectBname,  0 , UVAR_OPTIONAL, "Filename-base for (binary) spectrograms (one per IFO)");

  XLALregREALUserStruct  (Freq,              0,  UVAR_OPTIONAL, "physical start frequency to compute PSD for (excluding rngmed wings)");
  XLALregREALUserStruct  (FreqBand,          0,  UVAR_OPTIONAL, "physical frequency band to compute PSD for (excluding rngmed wings)");

  XLALregREALUserStruct  (startTime,        's', UVAR_OPTIONAL, "GPS timestamp of earliest SFT to include");
  XLALregREALUserStruct  (endTime,          'e', UVAR_OPTIONAL, "GPS timestamp of last SFT to include (NOTE: this refers to the SFT start-time!)");
  XLALregSTRINGUserStruct(timeStampsFile,   't', UVAR_OPTIONAL, "Time-stamps file");
  XLALregSTRINGUserStruct(IFO,               0 , UVAR_OPTIONAL, "Detector filter");

  XLALregINTUserStruct   (blocksRngMed,     'w', UVAR_OPTIONAL, "Running Median window size");

  XLALregINTUserStruct   (PSDmthopSFTs,     'S', UVAR_OPTIONAL, "For PSD, type of math. operation over SFTs: "
                                                                "0=arith-sum, 1=arith-mean, 2=arith-median, "
                                                                "3=harm-sum, 4=harm-mean, "
                                                                "5=power-2-sum, 6=power-2-mean, "
                                                                "7=min, 8=max");
  XLALregINTUserStruct   (PSDmthopIFOs,     'I', UVAR_OPTIONAL, "For PSD, type of math. op. over IFOs: "
                                                                "see --PSDmthopSFTs");
  XLALregBOOLUserStruct  (outputNormSFT,    'n', UVAR_OPTIONAL, "Output normalised SFT power to PSD file");
  XLALregINTUserStruct   (nSFTmthopSFTs,    'N', UVAR_OPTIONAL, "For norm. SFT, type of math. op. over SFTs: "
                                                                "see --PSDmthopSFTs");
  XLALregINTUserStruct   (nSFTmthopIFOs,    'J', UVAR_OPTIONAL, "For norm. SFT, type of math. op. over IFOs: "
                                                                "see --PSDmthopSFTs");

  XLALregINTUserStruct   (binSize,          'z', UVAR_OPTIONAL, "Bin the output into bins of size (in number of bins)");
  XLALregREALUserStruct  (binSizeHz,        'Z', UVAR_OPTIONAL, "Bin the output into bins of size (in Hz)");
  XLALregINTUserStruct   (PSDmthopBins,     'A', UVAR_OPTIONAL, "If binning, for PSD type of math. op. over bins: "
                                                                "see --PSDmthopSFTs");
  XLALregINTUserStruct   (nSFTmthopBins,    'B', UVAR_OPTIONAL, "If binning, for norm. SFT type of math. op. over bins: "
                                                                "see --PSDmthopSFTs");
  XLALregINTUserStruct   (binStep,          'p', UVAR_OPTIONAL, "If binning, step size to move bin along "
                                                                "(in number of bins, default is bin size)");
  XLALregREALUserStruct  (binStepHz,        'P', UVAR_OPTIONAL, "If binning, step size to move bin along "
                                                                "(in Hz, default is bin size)");
  XLALregBOOLUserStruct  (outFreqBinEnd,    'E', UVAR_OPTIONAL, "Output the end frequency of each bin");

  XLALregINTUserStruct   (maxBinsClean,     'm', UVAR_OPTIONAL, "Maximum Cleaning Bins");
  XLALregLISTUserStruct  (linefiles,         0 , UVAR_OPTIONAL, "Comma separated list of linefiles "
								"(names must contain IFO name)");

  XLALregBOOLUserStruct  (dumpMultiPSDVector,'d',UVAR_OPTIONAL, "Output multi-PSD vector over IFOs, timestamps, and frequencies into file(s) '<outputPSD>-IFO'");

  /* ----- developer options ---------- */
  XLALregREALUserStruct  (fStart,           'f', UVAR_DEVELOPER, "Start Frequency to load from SFT and compute PSD, including rngmed wings (BETTER: use --Freq instead)");
  XLALregREALUserStruct  (fBand,            'b', UVAR_DEVELOPER, "Frequency Band to load from SFT and compute PSD, including rngmed wings (BETTER: use --FreqBand instead)");

  XLALregBOOLUserStruct  (version,          'V', UVAR_SPECIAL, "Output version information");


  /* read all command line variables */
  if (XLALUserVarReadAllInput(argc, argv) != XLAL_SUCCESS)
    return XLAL_FAILURE;

  /* check user-input consistency */
  if (XLALUserVarWasSet(&(uvar->PSDmthopSFTs)) && !(0 <= uvar->PSDmthopSFTs && uvar->PSDmthopSFTs < MATH_OP_LAST)) {
    XLALPrintError("ERROR: --PSDmthopSFTs(-S) must be between 0 and %i", MATH_OP_LAST - 1);
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->PSDmthopIFOs)) && !(0 <= uvar->PSDmthopIFOs && uvar->PSDmthopIFOs < MATH_OP_LAST)) {
    XLALPrintError("ERROR: --PSDmthopIFOs(-I) must be between 0 and %i", MATH_OP_LAST - 1);
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->nSFTmthopSFTs)) && !(0 <= uvar->nSFTmthopSFTs && uvar->nSFTmthopSFTs < MATH_OP_LAST)) {
    XLALPrintError("ERROR: --nSFTmthopSFTs(-N) must be between 0 and %i", MATH_OP_LAST - 1);
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->nSFTmthopIFOs)) && !(0 <= uvar->nSFTmthopIFOs && uvar->nSFTmthopIFOs < MATH_OP_LAST)) {
    XLALPrintError("ERROR: --nSFTmthopIFOs(-J) must be between 0 and %i", MATH_OP_LAST - 1);
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->PSDmthopBins)) && !(0 <= uvar->PSDmthopBins && uvar->PSDmthopBins < MATH_OP_LAST)) {
    XLALPrintError("ERROR: --PSDmthopBins(-A) must be between 0 and %i", MATH_OP_LAST - 1);
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->nSFTmthopBins)) && !(0 <= uvar->nSFTmthopBins && uvar->nSFTmthopBins < MATH_OP_LAST)) {
    XLALPrintError("ERROR: --nSFTmthopBins(-B) must be between 0 and %i", MATH_OP_LAST - 1);
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->binSize)) && XLALUserVarWasSet(&(uvar->binSizeHz))) {
    XLALPrintError("ERROR: --binSize(-z) and --binSizeHz(-Z) are mutually exclusive");
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->binSize)) && uvar->binSize <= 0) {
    XLALPrintError("ERROR: --binSize(-z) must be strictly positive");
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->binSizeHz)) && uvar->binSizeHz <= 0.0) {
    XLALPrintError("ERROR: --binSizeHz(-Z) must be strictly positive");
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->binStep)) && XLALUserVarWasSet(&(uvar->binStepHz))) {
    XLALPrintError("ERROR: --binStep(-p) and --binStepHz(-P) are mutually exclusive");
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->binStep)) && uvar->binStep <= 0) {
    XLALPrintError("ERROR: --binStep(-p) must be strictly positive");
    return XLAL_FAILURE;
  }
  if (XLALUserVarWasSet(&(uvar->binStepHz)) && uvar->binStepHz <= 0.0) {
    XLALPrintError("ERROR: --binStepHz(-P) must be strictly positive");
    return XLAL_FAILURE;
  }

  return XLAL_SUCCESS;

} /* initUserVars() */


/**
 * Write a multi-PSD into spectrograms for each IFO.
 * Using gnuplot 'binary' matrix format
 * The filename for each IFO is generated as 'bname-IFO'
 */
void
LALfwriteSpectrograms ( LALStatus *status, const CHAR* bname, const MultiPSDVector *multiPSD )
{
  UINT4 X;
  CHAR *fname;
  float num, *row_data;		/* cast to float for writing (gnuplot binary format) */
  FILE *fp;

  INITSTATUS(status);
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

/**
 * Dump complete multi-PSDVector over IFOs, timestamps and frequency-bins into
 * per-IFO ASCII output-files 'outbname-IFO'
 *
 */
int
XLALDumpMultiPSDVector ( const CHAR *outbname,			/**< output basename 'outbname' */
                         const MultiPSDVector *multiPSDVect	/**< multi-psd vector to output */
                  )
{
  /* check input consistency */
  if ( outbname == NULL ) {
    XLALPrintError ("%s: NULL input 'outbname'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( multiPSDVect == NULL ) {
    XLALPrintError ("%s: NULL input 'multiPSDVect'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( multiPSDVect->length == 0 || multiPSDVect->data==0 ) {
    XLALPrintError ("%s: invalid multiPSDVect input (length=0 or data=NULL)\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  CHAR *fname;
  FILE *fp;

  UINT4 len = strlen ( outbname ) + 4;
  if ( ( fname = XLALMalloc ( len * sizeof(CHAR) )) == NULL ) {
    XLALPrintError ("%s: XLALMalloc(%d) failed.\n", __func__, len );
    XLAL_ERROR ( XLAL_ENOMEM);
  }

  UINT4 numIFOs = multiPSDVect->length;
  UINT4 X;
  for ( X = 0; X < numIFOs; X ++ )
    {
      PSDVector *thisPSDVect = multiPSDVect->data[X];
      char buf[100];

      sprintf ( fname, "%s-%c%c", outbname, thisPSDVect->data[0].name[0], thisPSDVect->data[0].name[1] );

      if ( ( fp = fopen( fname, "wb" ))  == NULL ) {
        XLALPrintError ("%s: Failed to open PSDperSFT file '%s' for writing!\n", __func__, fname );
        XLALFree ( fname );
        XLAL_ERROR ( XLAL_ESYS );
      }

      REAL8 f0       = thisPSDVect->data[0].f0;
      REAL8 dFreq    = thisPSDVect->data[0].deltaF;
      UINT4 numFreqs = thisPSDVect->data[0].data->length;
      UINT4 iFreq;

      /* write comment header line into this output file */
      /* FIXME: output code-version/cmdline/history info */
      fprintf(fp,"%%%% first line holds frequencies [Hz] of PSD-Columns\n");
      /* loop over frequency and output comnment-header markers */
      fprintf(fp,"%%%%%-17s", "dummy");
      for ( iFreq = 0; iFreq < numFreqs; iFreq ++ )
        {
          sprintf (buf, "f%d [Hz]", iFreq + 1 );
          fprintf (fp, " %-21s", buf );
        }
      fprintf (fp, "\n");

      /* write parseable header-line giving bin frequencies for PSDs */
      fprintf (fp, "%-19d", -1 );
      for (iFreq = 0; iFreq < numFreqs; iFreq++ )
        fprintf (fp, " %-21.16g", f0 + iFreq * dFreq );
      fprintf (fp, "\n\n\n");

      /* output another header line describing the following format "ti[GPS] PSD(f1) ... " */
      fprintf(fp,"%%%%%-17s", "ti[GPS]");
      for ( iFreq = 0; iFreq < numFreqs; iFreq ++ )
        {
          sprintf (buf, "PSD(f%d)", iFreq + 1 );
          fprintf (fp, " %-21s", buf );
        }
      fprintf (fp, "\n");

      /* loop over timestamps: dump all PSDs over frequencies into one line */
      UINT4 numTS = thisPSDVect->length;
      UINT4 iTS;
      for ( iTS = 0; iTS < numTS; iTS++ )
        {
          REAL8FrequencySeries *thisPSD = &thisPSDVect->data[iTS];

          /* first output timestamp GPS time for this line */
          REAL8 tGPS = XLALGPSGetREAL8( &thisPSD->epoch );
          fprintf (fp, "%-19.18g", tGPS );

          /* some internal consistency/paranoia checks */
          if ( ( f0 != thisPSD->f0) || ( dFreq != thisPSD->deltaF ) || (numFreqs != thisPSD->data->length ) ) {
            XLALPrintError ("%s: %d-th timestamp %f: inconsistent PSDVector: f0 = %g : %g,  dFreq = %g : %g, numFreqs = %d : %d \n",
                            __func__, iTS, tGPS, f0, thisPSD->f0, dFreq, thisPSD->deltaF, numFreqs, thisPSD->data->length );
            XLALFree ( fname );
            fclose ( fp );
            XLAL_ERROR ( XLAL_EDOM );
          }

          /* loop over all frequencies and dump PSD-value */
          for ( iFreq = 0; iFreq < numFreqs; iFreq ++ )
            fprintf (fp, " %-21.16g", thisPSD->data->data[iFreq] );

          fprintf (fp, "\n");

        } /* for iTS < numTS */

      fclose ( fp );

    } /* for X < numIFOs */

  XLALFree ( fname );

  return XLAL_SUCCESS;

} /* XLALDumpMultiPSDVector() */


/**
 * Load all SFTs according to user-input, returns multi-SFT vector.
 * \return cfg:
 * Returns 'effective' range of SFT-bins [firstBin, lastBin], which which the PSD will be estimated:
 * - if the user input {fStart, fBand} then these are loaded from SFTs and directly translated into bins
 * - if user input {Freq, FreqBand}, we load a wider frequency-band ADDING running-median/2 on either side
 * from the SFTs, and firstBind, lastBin correspond to {Freq,FreqBand} (rounded to closest bins)
 * Also returns the 'data-segment' for which SFTs were loaded
 *
 */
MultiSFTVector *
XLALReadSFTs ( ConfigVariables_t *cfg,		/**< [out] return derived configuration info (firstBin, lastBin, segment) */
               const UserVariables_t *uvar	/**< [in] complete user-input */
               )
{
  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_SFTConstraints;
  LIGOTimeGPS startTimeGPS = {0,0}, endTimeGPS = {0,0};
  LIGOTimeGPSVector *inputTimeStampsVector = NULL;

  /* check input */
  if ( !uvar || !uvar->inputData ) {
    XLALPrintError ("%s: invalid NULL input 'uvar' or 'uvar->inputData'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( !cfg ) {
    XLALPrintError ("%s: invalid NULL input 'cfg'", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* set detector constraint */
  if ( XLALUserVarWasSet ( &uvar->IFO ) )
    constraints.detector = uvar->IFO;
  else
    constraints.detector = NULL;

  if ( XLALUserVarWasSet( &uvar->startTime ) ) {
    XLALGPSSetREAL8 ( &startTimeGPS, uvar->startTime);
    constraints.minStartTime = &startTimeGPS;
  }

  if ( XLALUserVarWasSet( &uvar->endTime ) ) {
    XLALGPSSetREAL8 ( &endTimeGPS, uvar->endTime);
    constraints.maxStartTime = &endTimeGPS;
  }

  if ( XLALUserVarWasSet( &uvar->timeStampsFile ) ) {
    if ( (inputTimeStampsVector = XLALReadTimestampsFile ( uvar->timeStampsFile )) == NULL )
      XLAL_ERROR_NULL ( XLAL_EFUNC );

    constraints.timestamps = inputTimeStampsVector;
  }

  /* get sft catalog */
  LogPrintf ( LOG_DEBUG, "Finding all SFTs to load ... ");
  LALStatus status = blank_status;
  LALSFTdataFind ( &status, &catalog, uvar->inputData, &constraints);
  if ( status.statusCode != 0 ) {
    XLALPrintError ("%s: LALSFTdataFind() failed with statusCode = %d\n", __func__, status.statusCode );
    XLAL_ERROR_NULL ( XLAL_EFAILED );
  }
  if ( (catalog == NULL) || (catalog->length == 0) ) {
    XLALPrintError ("%s: Unable to match any SFTs with pattern '%s'\n", __func__, uvar->inputData );
    XLAL_ERROR_NULL ( XLAL_EFAILED );
  }
  LogPrintfVerbatim ( LOG_DEBUG, "done (found %i SFTs).\n", catalog->length);

  /* now we can free the inputTimeStampsVector */
  if ( inputTimeStampsVector )
    XLALDestroyTimestampVector ( inputTimeStampsVector );

  /* ----- some user-input consistency checks */
  BOOLEAN have_fStart   = XLALUserVarWasSet ( &uvar->fStart );
  BOOLEAN have_Freq     = XLALUserVarWasSet ( &uvar->Freq );
  BOOLEAN have_fBand    = XLALUserVarWasSet ( &uvar->fBand );
  BOOLEAN have_FreqBand = XLALUserVarWasSet ( &uvar->FreqBand );
  if ( have_fStart && have_Freq ) {
    XLALPrintError ("%s: use only one of --fStart OR --Freq (see --help)\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( have_fBand && have_FreqBand ) {
    XLALPrintError ("%s: use only one of --fBand OR --FreqBand (see --help)\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( ( have_fStart && have_FreqBand ) || ( have_Freq && have_fBand ) ) {
    XLALPrintError ("%s: don't mix {--fStart,--fBand} with {--Freq,--FreqBand} inputs (see --help)\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  /* ---------- figure out the right frequency-band to read from the SFTs, depending on user-input ----- */
  REAL8 fMin, fMax;
  UINT4 binsOffset; /* rngmed bin offset from start and end */
  UINT4 binsBand=0; /* width of physical FreqBand in bins */
  if ( have_Freq )
    {
      REAL8 dFreq = catalog->data[0].header.deltaF;
      binsOffset = uvar->blocksRngMed / 2 + 1;	/* truncates down plus add one bin extra safety! */
      binsBand   = ceil ( (uvar->FreqBand - 1e-9) / dFreq ) + 1; /* round up ! */

      REAL8 rngmedSideBand = binsOffset * dFreq;

      fMin = uvar->Freq - rngmedSideBand;
      fMax = uvar->Freq + uvar->FreqBand + rngmedSideBand;
    }
  else	/* NOTE: if no user-input on freq-band, we fall back to defaults on {fStart, fBand} */
    {
      fMin = uvar->fStart;
      fMax = uvar->fStart + uvar->fBand;
      binsOffset = 0;	/* no truncation of rngmed sidebands */
    }

  /* ----- figure out the data-segment span from the user-input and SFT-catalog ----- */
  /* if used passed these, then 'startTimeGPS' and 'endTimeGPS' are already set */
  if ( startTimeGPS.gpsSeconds == 0 )
    startTimeGPS = catalog->data[0].header.epoch;
  if ( endTimeGPS.gpsSeconds == 0 )
    endTimeGPS = catalog->data[catalog->length-1].header.epoch;
  /* SFT 'constraints' only refer to SFT *start-times*, for segment we need the end-time */
  REAL8 Tsft = 1.0 / catalog->data[0].header.deltaF;
  XLALGPSAdd ( &endTimeGPS, Tsft );

  /* ---------- read the sfts ---------- */
  LogPrintf (LOG_DEBUG, "Loading all SFTs ... ");
  MultiSFTVector *multi_sfts;
  if ( ( multi_sfts = XLALLoadMultiSFTs ( catalog, fMin, fMax ) ) == NULL ) {
    XLALPrintError ("%s: XLALLoadMultiSFTs( %f, %f ) failed with xlalErrno = %d\n", __func__, fMin, fMax, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  XLALDestroySFTCatalog ( catalog );
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");
  /* ---------- end loading SFTs ---------- */

  /* figure out effective PSD bin-boundaries for user */
  UINT4 numBins = multi_sfts->data[0]->data[0].data->length;
  INT4 bin0, bin1;
  if ( have_Freq )
    {
      bin0 = 0 + binsOffset;
      bin1 = bin0 + binsBand - 1;
    }
  else	/* output all bins loaded from SFTs (includes rngmed-sidebands) */
    {
      bin0 = 0;
      bin1 = numBins - 1;
    }

  /* return results */
  cfg->firstBin = (UINT4) bin0;
  cfg->lastBin = (UINT4) bin1;
  cfg->dataSegment.start = startTimeGPS;
  cfg->dataSegment.end   = endTimeGPS;

  XLALPrintInfo ("%s: loaded SFTs have %d bins, effective PSD output band is [%d, %d]\n", __func__, numBins, bin0, bin1 );

  return multi_sfts;

} /* XLALReadSFTs() */

/**
 * Function that *truncates the PSD in place* to the requested frequency-bin interval [firstBin, lastBin] for the given multiPSDVector.
 * Now also truncates the original SFT vector, as necessary for correct computation of normalized SFT power.
 */
int
XLALCropMultiPSDandSFTVectors ( MultiPSDVector *multiPSDVect,
                         MultiSFTVector *multiSFTVect,
                         UINT4 firstBin,
                         UINT4 lastBin
                         )
{
  /* check user input */
  if ( !multiPSDVect ) {
    XLALPrintError ("%s: invalid NULL input 'multiPSDVect'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( !multiSFTVect ) {
    XLALPrintError ("%s: invalid NULL input 'multiSFTVect'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( lastBin < firstBin ) {
    XLALPrintError ("%s: empty bin interval requested [%d, %d]\n", __func__, firstBin, lastBin );
    XLAL_ERROR ( XLAL_EDOM );
  }

  UINT4 numIFOs = multiPSDVect->length;
  UINT4 numBins = multiPSDVect->data[0]->data[0].data->length;

  if ( numIFOs != multiSFTVect->length ) {
    XLALPrintError ("%s: inconsistent number of IFOs between PSD (%d) and SFT (%d) vectors.\n", __func__, numIFOs, multiSFTVect->length );
    XLAL_ERROR ( XLAL_EDOM );
  }

  if ( numBins != multiSFTVect->data[0]->data[0].data->length ) {
    XLALPrintError ("%s: inconsistent number of bins between PSD (%d bins) and SFT (%d bins) vectors.\n", __func__, numBins, multiSFTVect->data[0]->data[0].data->length );
    XLAL_ERROR ( XLAL_EDOM );
  }

  if ( (firstBin >= numBins) || (lastBin >= numBins ) ) {
    XLALPrintError ("%s: requested bin-interval [%d, %d] outside of PSD bins [0, %d]\n", __func__, firstBin, lastBin, 0, numBins - 1 );
    XLAL_ERROR ( XLAL_EDOM );
  }

  /* ----- check if there's anything to do at all? ----- */
  if ( (firstBin == 0)  && (lastBin == numBins - 1) )
    return XLAL_SUCCESS;

  /* ----- loop over detectors, timestamps, then crop each PSD ----- */
  UINT4 X;
  for ( X=0; X < numIFOs; X ++ )
    {
      PSDVector *thisPSDVect = multiPSDVect->data[X];
      SFTVector *thisSFTVect = multiSFTVect->data[X];
      UINT4 numTS   = thisPSDVect->length;

      UINT4 iTS;
      for ( iTS = 0; iTS < numTS; iTS ++ )
        {
          REAL8FrequencySeries *thisPSD = &thisPSDVect->data[iTS];
          COMPLEX8FrequencySeries *thisSFT = &thisSFTVect->data[iTS];

          if ( numBins != thisPSD->data->length ) {
            XLALPrintError ("%s: inconsistent number of frequency-bins across multiPSDVector: X=%d, iTS=%d: numBins = %d != %d\n",
                            __func__, X, iTS, numBins, thisPSD->data->length );
            XLAL_ERROR ( XLAL_EDOM );
          }

          if ( numBins != thisSFT->data->length ) {
            XLALPrintError ("%s: inconsistent number of frequency-bins across multiSFTVector: X=%d, iTS=%d: numBins = %d != %d\n",
                            __func__, X, iTS, numBins, thisSFT->data->length );
            XLAL_ERROR ( XLAL_EDOM );
          }

          UINT4 numNewBins = lastBin - firstBin + 1;

          /* crop PSD */
          XLAL_CHECK(XLALResizeREAL8FrequencySeries(thisPSD, firstBin, numNewBins) != NULL, XLAL_EFUNC);

          /* crop SFT */
          XLAL_CHECK(XLALResizeCOMPLEX8FrequencySeries(thisSFT, firstBin, numNewBins) != NULL, XLAL_EFUNC);

        } /* for iTS < numTS */

    } /* for X < numIFOs */

  /* that should be all ... */
  return XLAL_SUCCESS;

} /* XLALCropMultiPSDandSFTVectors() */

/**
 * Compute the "data-quality factor" \f$\mathcal{Q}(f) = \sum_X \frac{\epsilon_X}{\mathcal{S}_X(f)}\f$ over the given SFTs.
 * The input \a multiPSD is a pre-computed PSD map \f$\mathcal{S}_{X,i}(f)\f$, over IFOs \f$X\f$, SFTs \f$i\f$
 * and frequency \f$f\f$.
 *
 * \return the output is a vector \f$\mathcal{Q}(f)\f$.
 *
 */
REAL8FrequencySeries *
XLALComputeSegmentDataQ ( const MultiPSDVector *multiPSDVect, 	/**< input PSD map over IFOs, SFTs, and frequencies */
                          LALSeg segment		  	/**< segment to compute Q for */
                          )
{
  /* check input consistency */
  if ( multiPSDVect == NULL ) {
    XLALPrintError ("%s: NULL input 'multiPSDVect'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( multiPSDVect->length == 0 || multiPSDVect->data==0 ) {
    XLALPrintError ("%s: invalid multiPSDVect input (length=0 or data=NULL)\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  REAL8 Tseg = XLALGPSDiff ( &segment.end, &segment.start );
  if ( Tseg <= 0 ) {
    XLALPrintError ("%s: negative segment-duration '%g'\n", __func__, Tseg );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  REAL8 Tsft 	 = 1.0 / multiPSDVect->data[0]->data[0].deltaF;
  REAL8 f0       = multiPSDVect->data[0]->data[0].f0;
  REAL8 dFreq    = multiPSDVect->data[0]->data[0].deltaF;
  UINT4 numFreqs = multiPSDVect->data[0]->data[0].data->length;

  REAL8FrequencySeries *Q, *SXinv;
  if ( (Q = XLALCreateREAL8FrequencySeries ( "Qfactor", &segment.start, f0, dFreq, &lalHertzUnit, numFreqs )) == NULL ) {
    XLALPrintError ("%s: Q = XLALCreateREAL8FrequencySeries(numFreqs=%d) failed with xlalErrno = %d\n", __func__, numFreqs, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  if ( (SXinv = XLALCreateREAL8FrequencySeries ( "SXinv", &segment.start, f0, dFreq, &lalHertzUnit, numFreqs )) == NULL ) {
    XLALPrintError ("%s: SXinv = XLALCreateREAL8FrequencySeries(numFreqs=%d) failed with xlalErrno = %d\n", __func__, numFreqs, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  /* initialize Q-array to zero, as we'll keep adding to it */
  memset ( Q->data->data, 0, Q->data->length * sizeof(Q->data->data[0]) );

  /* ----- loop over IFOs ----- */
  UINT4 numIFOs = multiPSDVect->length;
  UINT4 X;
  for ( X = 0; X < numIFOs; X ++ )
    {
      PSDVector *thisPSDVect = multiPSDVect->data[X];

      /* initialize SXinv-array to zero, as we'll keep adding to it */
      memset ( SXinv->data->data, 0, SXinv->data->length * sizeof(SXinv->data->data[0]) );
      UINT4 numSFTsInSeg = 0;	/* reset counter of SFTs within this segment */

      /* ----- loop over all timestamps ----- */
      /* find SFTs inside segment, count them and combine their PSDs */
      UINT4 numTS = thisPSDVect->length;
      UINT4 iTS;
      for ( iTS = 0; iTS < numTS; iTS++ )
        {
          REAL8FrequencySeries *thisPSD = &thisPSDVect->data[iTS];

          /* some internal consistency/paranoia checks */
          if ( ( f0 != thisPSD->f0) || ( dFreq != thisPSD->deltaF ) || (numFreqs != thisPSD->data->length ) ) {
            XLALPrintError ("%s: %d-th timestamp %f: inconsistent PSDVector: f0 = %g : %g,  dFreq = %g : %g, numFreqs = %d : %d \n",
                            __func__, iTS, XLALGPSGetREAL8( &thisPSD->epoch ), f0, thisPSD->f0, dFreq, thisPSD->deltaF, numFreqs, thisPSD->data->length );
            XLAL_ERROR_NULL ( XLAL_EDOM );
          }

          int cmp = XLALCWGPSinRange( thisPSD->epoch, &segment.start, &segment.end );

          if ( cmp < 0 )	/* SFT-end before segment => advance to the next one */
            continue;
          if ( cmp > 0 )	/* SFT-start past end of segment: ==> terminate loop */
            break;

          if ( cmp == 0 )	/* this SFT is inside segment */
            {
              numSFTsInSeg ++;
              /* add SXinv(f) += 1/SX_i(f) over all frequencies */
              UINT4 iFreq;
              for ( iFreq = 0; iFreq < numFreqs; iFreq++ )
                SXinv->data->data[iFreq] += 1.0 / thisPSD->data->data[iFreq] ;

            }	/* if SFT inside segment */

        } /* for iTS < numTS */

      /* compute duty-cycle eps_X = nX * Tsft / Tseg for this IFO */
      REAL8 duty_X = numSFTsInSeg * Tsft / Tseg;
      /* sanity check: eps in [0, 1]*/
      if ( (duty_X < 0) && (duty_X > 1 ) ) {
        XLALPrintError ("%s: something is WRONG: duty-cyle = %g not within [0,1]!\n", __func__, duty_X );
        XLAL_ERROR_NULL ( XLAL_EFAILED );
      }

      /* add duty_X-weighted SXinv to Q */
      UINT4 iFreq;
      for ( iFreq = 0; iFreq < numFreqs; iFreq ++ )
        Q->data->data[iFreq] += duty_X * SXinv->data->data[iFreq] / numSFTsInSeg;

    } /* for X < numIFOs */

  /* clean up, free memory */
  XLALDestroyREAL8FrequencySeries ( SXinv );


  return Q;

} /* XLALComputeSegmentDataQ() */

/**
 * Write given REAL8FrequencySeries into file
 */
int
XLALWriteREAL8FrequencySeries_to_file ( const REAL8FrequencySeries *series,	/**< [in] frequency-series to write to file */
                                        const char *fname			/**< [in] filename to write into */
                                        )
{
  /* check input consistency */
  if ( !series || !fname ) {
    XLALPrintError ("%s: invalid NULL input.\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  FILE *fp;
  if ( ( fp = fopen ( fname, "wb" )) == NULL ) {
    XLALPrintError ("%s: failed to open file '%s' for writing.\n", __func__, fname );
    XLAL_ERROR ( XLAL_ESYS );
  }

  /* write header info in comments */
  if ( XLAL_SUCCESS != XLALOutputVersionString ( fp, 0 ) )
    XLAL_ERROR ( XLAL_EFUNC );

  fprintf ( fp, "%%%% name = '%s'\n", series->name );
  fprintf ( fp, "%%%% epoch = {%d, %d}\n", series->epoch.gpsSeconds, series->epoch.gpsNanoSeconds );
  fprintf ( fp, "%%%% f0 = %f Hz\n", series->f0 );
  fprintf ( fp, "%%%% deltaF = %g Hz\n", series->deltaF );

CHAR unitStr[1024];
 if ( XLALUnitAsString( &unitStr[0], sizeof(unitStr)-1, &series->sampleUnits ) == NULL ) {
   XLALPrintError ("%s: XLALUnitAsString() failed with xlalErrno = %d.\n", __func__, xlalErrno );
   XLAL_ERROR ( XLAL_EFUNC );
 }
 fprintf ( fp, "%%%% Units = %s\n", unitStr );

 fprintf ( fp, "%%%% Freq [Hz]           Data(Freq)\n");
 UINT4 numBins = series->data->length;
 UINT4 iFreq;
 for ( iFreq = 0; iFreq < numBins; iFreq ++ )
   {
     REAL8 thisFreq = series->f0 + iFreq * series->deltaF;
     fprintf (fp, "%20.16f  %20.16g\n", thisFreq, series->data->data[iFreq] );
   }

  fclose ( fp );

  return XLAL_SUCCESS;

} /* XLALWriteREAL8FrequencySeries_to_file() */
