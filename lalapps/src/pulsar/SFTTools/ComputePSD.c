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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \file
 * \ingroup lalapps_pulsar_SFTTools
 * \author Badri Krishnan, Iraj Gholami, Reinhard Prix, Alicia Sintes, Karl Wette
 * \brief Compute power spectral densities
 */

#include "config.h"

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
#include <lal/SFTutils.h>

#include <lal/LogPrintf.h>

#include <LALAppsVCSInfo.h>

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

/** user input variables */
typedef struct
{
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
  BOOLEAN normalizeByTotalNumSFTs; /**< apply normalization factor from total number of SFTs over all IFOs */

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

} UserVariables_t;

/**
 * Config variables 'derived' from user-input
 */
typedef struct
{
  REAL8 FreqMin;		/**< frequency of first PSD bin for output */
  REAL8 FreqBand;		/**< width of frequency band for output */
  LALSeg dataSegment;		/**< the data-segment for which PSD was computed */
} ConfigVariables_t;

/* ---------- global variables ----------*/
extern int vrbflg;

/* ---------- local prototypes ---------- */
int initUserVars (int argc, char *argv[], UserVariables_t *uvar);
void LALfwriteSpectrograms ( LALStatus *status, const CHAR *bname, const MultiPSDVector *multiPSD );
MultiSFTVector *XLALReadSFTs ( ConfigVariables_t *cfg, const UserVariables_t *uvar );

int XLALWriteREAL8FrequencySeries_to_file ( const REAL8FrequencySeries *series, const char *fname );

/*============================================================
 * FUNCTION definitions
 *============================================================*/
int
main(int argc, char *argv[])
{
  static LALStatus       status;  /* LALStatus pointer */
  UserVariables_t XLAL_INIT_DECL(uvar);
  ConfigVariables_t XLAL_INIT_DECL(cfg);

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register and read user variables */
  if (initUserVars(argc, argv, &uvar) != XLAL_SUCCESS)
    return EXIT_FAILURE;

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

  /* call the main loop function; output vectors will be allocated inside
   * NOTE: inputSFTs will be normalized in place for efficiency reasons
   */
  REAL8Vector *finalPSD = NULL;
  MultiPSDVector *multiPSDVector = NULL;
  REAL8Vector *normSFT = NULL;
  BOOLEAN returnMultiPSDVector = ( uvar.outputSpectBname || uvar.dumpMultiPSDVector || uvar.outputQ );
  XLAL_CHECK_MAIN ( XLALComputePSDandNormSFTPower ( &finalPSD,
                      &multiPSDVector,
                      &normSFT,
                      inputSFTs,
                      returnMultiPSDVector,
                      uvar.outputNormSFT,
                      uvar.blocksRngMed,
                      uvar.PSDmthopSFTs,
                      uvar.PSDmthopIFOs,
                      uvar.nSFTmthopSFTs,
                      uvar.nSFTmthopIFOs,
                      uvar.normalizeByTotalNumSFTs,
                      cfg.FreqMin,
                      cfg.FreqBand,
                      TRUE // normalizeSFTsInPlace
                  ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* output spectrograms */
  if ( uvar.outputSpectBname ) {
    LAL_CALL ( LALfwriteSpectrograms ( &status, uvar.outputSpectBname, multiPSDVector ), &status );
  }

  /* ---------- if user requested it, output complete MultiPSDVector over IFOs X, timestamps and freq-bins into ASCI file(s) */
  if ( uvar.dumpMultiPSDVector ) {
    if ( XLALDumpMultiPSDVector ( uvar.outputPSD, multiPSDVector ) != XLAL_SUCCESS ) {
      XLALPrintError ("%s: XLALDumpMultiPSDVector() failed, xlalErrnor = %d\n", __func__, xlalErrno );
      return EXIT_FAILURE;
    }
  } /* if uvar.dumpMultiPSDVector */

  /* ----- if requested, compute data-quality factor 'Q' -------------------- */
  if ( uvar.outputQ )
    {
      REAL8FrequencySeries *Q;
      if ( (Q = XLALComputeSegmentDataQ ( multiPSDVector, cfg.dataSegment )) == NULL ) {
        XLALPrintError ("%s: XLALComputeSegmentDataQ() failed with xlalErrno = %d\n", __func__, xlalErrno );
        return EXIT_FAILURE;
      }
      if ( XLAL_SUCCESS != XLALWriteREAL8FrequencySeries_to_file ( Q, uvar.outputQ ) ) {
        return EXIT_FAILURE;
      }
      XLALDestroyREAL8FrequencySeries ( Q );
    } /* if outputQ */

  /* ---------- BINNING if requested ---------- */
  /* start frequency and frequency spacing */
  REAL8 Freq0 = inputSFTs->data[0]->data[0].f0;
  REAL8 dFreq = inputSFTs->data[0]->data[0].deltaF;
  /* work out bin size */
  UINT4 finalBinSize;
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
  UINT4 finalBinStep;
  if (XLALUserVarWasSet(&uvar.binStep)) {
    finalBinStep = uvar.binStep;
  }
  else if (XLALUserVarWasSet(&uvar.binStepHz)) {
    finalBinStep = (UINT4)floor(uvar.binStepHz / dFreq + 0.5); /* round to nearest bin */
  }
  else {
    finalBinStep = finalBinSize;
  }

  /* write final PSD to file */
  if (XLALUserVarWasSet(&uvar.outputPSD)) {
    LogPrintf(LOG_DEBUG, "Printing PSD to file ...\n");
    FILE *fpOut = NULL;
    XLAL_CHECK_MAIN ( (fpOut = fopen(uvar.outputPSD, "wb")) != NULL, XLAL_EIO, "Unable to open output file %s for writing", uvar.outputPSD );
    XLAL_CHECK_MAIN ( XLALOutputVCSInfo(fpOut, lalAppsVCSInfoList, 0, "%% ") == XLAL_SUCCESS, XLAL_EFUNC );
    for (int a = 0; a < argc; a++) { /* write the command-line */
      fprintf(fpOut,"%%%% argv[%d]: '%s'\n", a, argv[a]);
    }
    XLAL_CHECK_MAIN ( XLALWritePSDtoFilePointer ( fpOut, finalPSD, normSFT, uvar.outputNormSFT, uvar.outFreqBinEnd, uvar.PSDmthopBins, uvar.nSFTmthopBins, finalBinSize, finalBinStep, Freq0, dFreq ) == XLAL_SUCCESS, XLAL_EFUNC );
    LogPrintfVerbatim ( LOG_DEBUG, "done.\n");
    fclose(fpOut);
  }

  /* we are now done with the psd */
  if ( returnMultiPSDVector ) {
    XLALDestroyMultiPSDVector ( multiPSDVector);
  }
  XLALDestroyMultiSFTVector ( inputSFTs);
  XLALDestroyREAL8Vector ( finalPSD );
  XLALDestroyREAL8Vector ( normSFT );

  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

} /* main() */


/** register all "user-variables" */
int
initUserVars (int argc, char *argv[], UserVariables_t *uvar)
{

  /* set a few defaults */
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

  uvar->normalizeByTotalNumSFTs = FALSE;

  uvar->dumpMultiPSDVector = FALSE;

  uvar->binSizeHz = 0.0;
  uvar->binSize   = 1;
  uvar->PSDmthopBins  = MATH_OP_ARITHMETIC_MEDIAN;
  uvar->nSFTmthopBins = MATH_OP_MAXIMUM;
  uvar->binStep   = 0.0;
  uvar->binStep   = 1;

  /* register user input variables */
  XLALRegisterUvarMember(inputData,        STRING, 'i', REQUIRED, "Input SFT pattern");
  XLALRegisterUvarMember(outputPSD,        STRING, 'o', OPTIONAL, "Output PSD into this file");
  XLALRegisterUvarMember(outputQ,	     STRING, 0,  OPTIONAL, "Output the 'data-quality factor' Q(f) into this file");
  XLALRegisterUvarMember(outputSpectBname,  STRING, 0 , OPTIONAL, "Filename-base for (binary) spectrograms (one per IFO)");

  XLALRegisterUvarMember(Freq,              REAL8, 0,  OPTIONAL, "physical start frequency to compute PSD for (excluding rngmed wings)");
  XLALRegisterUvarMember(FreqBand,          REAL8, 0,  OPTIONAL, "physical frequency band to compute PSD for (excluding rngmed wings)");

  XLALRegisterUvarMember(startTime,        REAL8, 's', OPTIONAL, "SFT timestamps must be >= this GPS timestamp");
  XLALRegisterUvarMember(endTime,          REAL8, 'e', OPTIONAL, "SFT timestamps must be < this GPS timestamp");
  XLALRegisterUvarMember(timeStampsFile,   STRING, 't', OPTIONAL, "Time-stamps file");
  XLALRegisterUvarMember(IFO,               STRING, 0 , OPTIONAL, "Detector filter");

  XLALRegisterUvarMember(blocksRngMed,     INT4, 'w', OPTIONAL, "Running Median window size");

  XLALRegisterUvarAuxDataMember(PSDmthopSFTs,     UserEnum, &MathOpTypeChoices, 'S', OPTIONAL, "For PSD, type of math. operation over SFTs, can be given by string names (preferred) or legacy numbers: \n"
                                                                "arithsum (0), arithmean (1), arithmedian (2), "
                                                                "harmsum (3), harmmean (4), "
                                                                "powerminus2sum (5), powerminus2mean (6), "
                                                                "min (7), max (8)");
  XLALRegisterUvarAuxDataMember(PSDmthopIFOs,     UserEnum, &MathOpTypeChoices, 'I', OPTIONAL, "For PSD, type of math. op. over IFOs: "
                                                                "see --PSDmthopSFTs");
  XLALRegisterUvarMember(outputNormSFT,    BOOLEAN, 'n', OPTIONAL, "Output normalised SFT power to PSD file");
  XLALRegisterUvarAuxDataMember(nSFTmthopSFTs,    UserEnum, &MathOpTypeChoices, 'N', OPTIONAL, "For norm. SFT, type of math. op. over SFTs: "
                                                                "see --PSDmthopSFTs");
  XLALRegisterUvarAuxDataMember(nSFTmthopIFOs,    UserEnum, &MathOpTypeChoices, 'J', OPTIONAL, "For norm. SFT, type of math. op. over IFOs: "
                                                                "see --PSDmthopSFTs");

  XLALRegisterUvarMember(normalizeByTotalNumSFTs,    BOOLEAN, 0, OPTIONAL, "For harmsum/powerminus2sum, apply normalization factor from total number of SFTs over all IFOs (mimics harmmean/powerminus2mean over a combined set of all SFTs)");

  XLALRegisterUvarMember(binSize,          INT4, 'z', OPTIONAL, "Bin the output into bins of size (in number of bins)");
  XLALRegisterUvarMember(binSizeHz,        REAL8, 'Z', OPTIONAL, "Bin the output into bins of size (in Hz)");
  XLALRegisterUvarAuxDataMember(PSDmthopBins,     UserEnum, &MathOpTypeChoices, 'A', OPTIONAL, "If binning, for PSD type of math. op. over bins: "
                                                                "see --PSDmthopSFTs");
  XLALRegisterUvarAuxDataMember(nSFTmthopBins,    UserEnum, &MathOpTypeChoices, 'B', OPTIONAL, "If binning, for norm. SFT type of math. op. over bins: "
                                                                "see --PSDmthopSFTs");
  XLALRegisterUvarMember(binStep,          INT4, 'p', OPTIONAL, "If binning, step size to move bin along "
                                                                "(in number of bins, default is bin size)");
  XLALRegisterUvarMember(binStepHz,        REAL8, 'P', OPTIONAL, "If binning, step size to move bin along "
                                                                "(in Hz, default is bin size)");
  XLALRegisterUvarMember(outFreqBinEnd,    BOOLEAN, 'E', OPTIONAL, "Output the end frequency of each bin");

  XLALRegisterUvarMember(maxBinsClean,     INT4, 'm', OPTIONAL, "Maximum Cleaning Bins");
  XLALRegisterUvarMember(linefiles,         STRINGVector, 0 , OPTIONAL, "Comma separated list of linefiles "
								"(names must contain IFO name)");

  XLALRegisterUvarMember(dumpMultiPSDVector,BOOLEAN, 'd',OPTIONAL, "Output multi-PSD vector over IFOs, timestamps, and frequencies into file(s) '<outputPSD>-IFO'");

  /* ----- developer options ---------- */
  XLALRegisterUvarMember(fStart,           REAL8, 'f', DEVELOPER, "Start Frequency to load from SFT and compute PSD, including rngmed wings (BETTER: use --Freq instead)");
  XLALRegisterUvarMember(fBand,            REAL8, 'b', DEVELOPER, "Frequency Band to load from SFT and compute PSD, including rngmed wings (BETTER: use --FreqBand instead)");


  /* read all command line variables */
  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv, lalAppsVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit(1);
  }

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
  BOOLEAN have_fStart   = XLALUserVarWasSet ( &uvar->fStart );
  BOOLEAN have_Freq     = XLALUserVarWasSet ( &uvar->Freq );
  BOOLEAN have_fBand    = XLALUserVarWasSet ( &uvar->fBand );
  BOOLEAN have_FreqBand = XLALUserVarWasSet ( &uvar->FreqBand );
  XLAL_CHECK ( !(have_fStart && have_Freq), XLAL_EINVAL, "use only one of --fStart OR --Freq (see --help)" );
  XLAL_CHECK ( !(have_fBand && have_FreqBand), XLAL_EINVAL, "use only one of --fBand OR --FreqBand (see --help)" );
  XLAL_CHECK ( ! (( have_fStart && have_FreqBand ) || ( have_Freq && have_fBand )), XLAL_EINVAL, "don't mix {--fStart,--fBand} with {--Freq,--FreqBand} inputs (see --help)");

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
  CHAR *fname;
  float num, *row_data;		/* cast to float for writing (gnuplot binary format) */
  FILE *fp;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  if ( !bname || !multiPSD || multiPSD->length == 0 ) {
    ABORT ( status, COMPUTEPSDC_ENULL, COMPUTEPSDC_MSGENULL );
  }

  /* loop over IFOs */
  for ( UINT4 X = 0; X < multiPSD->length ; X ++ )
    {
      UINT4 len = strlen ( bname ) + 4;	/* append '-XN' to get IFO-specific filename */
      const CHAR *tmp;
      REAL8 f0, df;

      UINT4 numSFTs = multiPSD->data[X]->length;
      UINT4 numBins = multiPSD->data[X]->data[0].data->length;

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
      for ( UINT4 k=0; k < numBins; k ++ )
        row_data[k] = (float) ( f0 + 1.0 * k * df );
      if ( fwrite((char *) row_data, sizeof(float), numBins, fp) != numBins ) {
        LogPrintf (LOG_CRITICAL, "Failed to fwrite() to spectrogram file '%s'\n", fname );
        goto failed;
      }

      /* write PSDs of successive SFTs in rows, first column is GPS-time in seconds */
      for ( UINT4 j = 0; j < numSFTs ; j ++ ) {
        num = (float) multiPSD->data[X]->data[j].epoch.gpsSeconds;
        for ( UINT4 k = 0; k < numBins; k ++ )
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
  SFTConstraints XLAL_INIT_DECL(constraints);
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
  if ( ( catalog = XLALSFTdataFind ( uvar->inputData, &constraints) ) == NULL ) {
    XLALPrintError ("%s: XLALSFTdataFind() failed with xlalErrno = %d\n", __func__, xlalErrno );
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

  /* ---------- figure out the right frequency-band to read from the SFTs, depending on user-input ----- */
  REAL8 fMin, fMax;
  if ( XLALUserVarWasSet ( &uvar->Freq ) )
    {
      REAL8 dFreq = catalog->data[0].header.deltaF;
      /* rngmed bin offset from start and end */
      UINT4 rngmedSideBandBins = uvar->blocksRngMed / 2 + 1; /* truncates down plus add one bin extra safety! */
      REAL8 rngmedSideBand = rngmedSideBandBins * dFreq;
      fMin = uvar->Freq - rngmedSideBand;
      fMax = uvar->Freq + uvar->FreqBand + rngmedSideBand;
      cfg->FreqMin  = uvar->Freq;
      cfg->FreqBand = uvar->FreqBand;
    }
  else
    {
      /* if no user-input on freq-band, we fall back to defaults on {fStart, fBand} */
      /* (no truncation of rngmed sidebands) */
      fMin = uvar->fStart;
      fMax = uvar->fStart + uvar->fBand;
      cfg->FreqMin  = uvar->fStart;
      cfg->FreqBand = uvar->fBand;
    }

  /* ----- figure out the data-segment span from the user-input and SFT-catalog ----- */
  /* if used passed these, then 'startTimeGPS' and 'endTimeGPS' are already set */
  if ( startTimeGPS.gpsSeconds == 0 )
    startTimeGPS = catalog->data[0].header.epoch;
  if ( endTimeGPS.gpsSeconds == 0 )
    endTimeGPS = catalog->data[catalog->length-1].header.epoch;
  /* SFT 'constraints' only refer to SFT *start-times*, for segment we need the end-time */
  REAL8 deltaF = catalog->data[0].header.deltaF;
  REAL8 Tsft = 1.0 / deltaF;
  XLALGPSAdd ( &endTimeGPS, Tsft );

  /* ---------- read the sfts ---------- */
  LogPrintf (LOG_DEBUG, "Loading all SFTs over frequency band [%f,%f]...\n", fMin, fMax);
  MultiSFTVector *multi_sfts;
  if ( ( multi_sfts = XLALLoadMultiSFTs ( catalog, fMin, fMax ) ) == NULL ) {
    XLALPrintError ("%s: XLALLoadMultiSFTs( %f, %f ) failed with xlalErrno = %d\n", __func__, fMin, fMax, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  XLALDestroySFTCatalog ( catalog );
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");
  /* ---------- end loading SFTs ---------- */

  /* return results */
  cfg->dataSegment.start = startTimeGPS;
  cfg->dataSegment.end   = endTimeGPS;

  UINT4 numBins = multi_sfts->data[0]->data[0].data->length;
  REAL8 f0sfts = multi_sfts->data[0]->data[0].f0;
  LogPrintf (LOG_DEBUG, "Loaded SFTs have %d bins, sampled at %fHz, covering frequency band [%f, %f]\n", numBins, deltaF, f0sfts, f0sfts+numBins*deltaF );

  return multi_sfts;

} /* XLALReadSFTs() */


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
  if ( XLAL_SUCCESS != XLALOutputVCSInfo(fp, lalAppsVCSInfoList, 0, "%% ") )
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
 for ( UINT4 iFreq = 0; iFreq < numBins; iFreq ++ )
   {
     REAL8 thisFreq = series->f0 + iFreq * series->deltaF;
     fprintf (fp, "%20.16f  %20.16g\n", thisFreq, series->data->data[iFreq] );
   }

  fclose ( fp );

  return XLAL_SUCCESS;

} /* XLALWriteREAL8FrequencySeries_to_file() */
