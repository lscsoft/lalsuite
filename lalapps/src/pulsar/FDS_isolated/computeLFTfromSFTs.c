/*
 * Copyright (C) 2009 Reinhard Prix
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

/*********************************************************************************/
/** \author R. Prix
 * \file
 * \brief
 * Calculate the Fourier transform over the total timespan from a set of SFTs
 *
 *********************************************************************************/
#include "config.h"

/* System includes */
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif


/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LogPrintf.h>
#include <lal/TimeSeries.h>
#include <lal/ComplexFFT.h>

#include <lal/ComputeFstat.h>

#include <lal/lalGitID.h>
#include <lalappsGitID.h>

#include <lalapps.h>

/* local includes */

RCSID( "$Id$");

/*---------- DEFINES ----------*/
#define DTERMS 32

/** convert GPS-time to REAL8 */
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define SQ(x) ((x)*(x))

#define NhalfPosDC(N) ((UINT4)(ceil ( ((N)/2.0 - 1e-6 ))))	/* round up */
#define NhalfNeg(N) ((UINT4)( (N) - NhalfPosDC(N) ))		/* round down (making sure N+ + N- = (N-1) */


#define LAL_INT4_MAX    LAL_UINT8_C(2147483647)

/*----- Error-codes -----*/
#define LFTFROMSFTS_ENULL 	1
#define LFTFROMSFTS_ESYS     	2
#define LFTFROMSFTS_EINPUT   	3
#define LFTFROMSFTS_EMEM   	4
#define LFTFROMSFTS_ENONULL 	5
#define LFTFROMSFTS_EXLAL	6

#define LFTFROMSFTS_MSGENULL 	"Arguments contained an unexpected null pointer"
#define LFTFROMSFTS_MSGESYS	"System call failed (probably file IO)"
#define LFTFROMSFTS_MSGEINPUT  "Invalid input"
#define LFTFROMSFTS_MSGEMEM   	"Out of memory. Bad."
#define LFTFROMSFTS_MSGENONULL "Output pointer is non-NULL"
#define LFTFROMSFTS_MSGEXLAL	"XLALFunction-call failed"

/** Structure containing input SFTs plus useful meta-data about those SFTs.
 */
typedef struct {
  MultiSFTVector *multiSFTs;	/**< input SFT vector */
  CHAR *dataSummary;            /**< descriptive string describing the data */
  UINT4 numDet;			/**< number of detectors in multiSFT vector */
  REAL8 Tsft;			/**< duration of each SFT */
  LIGOTimeGPS startTime;	/**< start-time of SFTs */
  LIGOTimeGPS endTime;
  REAL8 fmin;			/**< smallest frequency contained in input SFTs */
  UINT4 numBins;		/**< number of frequency bins in input SFTs */
} InputSFTData;


/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/* ----- User-variables: can be set from config-file or command-line */
typedef struct {
  BOOLEAN help;		/**< trigger output of help string */

  CHAR *inputSFTs;	/**< SFT input-files to use */
  CHAR *outputLFT;	/**< output ('Long Fourier Transform') file to write total-time Fourier transform into */
  INT4 minStartTime;	/**< limit start-time of input SFTs to use */
  INT4 maxEndTime;	/**< limit end-time of input SFTs to use */

  REAL8 fmin;		/**< minimal frequency to include */
  REAL8 fmax;		/**< maximal frequency to include */

  REAL8 upsampling;	/**< factor by which to upsample the frequency resolution by */

  BOOLEAN version;	/**< output version-info */
} UserInput_t;

static UserInput_t empty_UserInput;
static LALUnit empty_LALUnit;

/* ---------- local prototypes ---------- */
int main(int argc, char *argv[]);

void LALInitUserVars ( LALStatus *status, UserInput_t *uvar);
void LoadInputSFTs ( LALStatus *status, InputSFTData *data, const UserInput_t *uvar );
int XLALWeightedSumOverSFTVector ( SFTVector **outSFT, const SFTVector *inVect, const COMPLEX8Vector *weights );

SFTtype *XLALSFTVectorToLFT ( const SFTVector *sfts, REAL8 upsampling );
int XLALReorderFFTWtoSFT (COMPLEX8Vector *X);
int XLALReorderSFTtoFFTW (COMPLEX8Vector *X);

int XLALTimeShiftSFT ( SFTtype *sft, REAL8 shift );

/*---------- empty initializers ---------- */

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function
 */
int main(int argc, char *argv[])
{
  UserInput_t uvar = empty_UserInput;
  LALStatus status = blank_status;	/* initialize status */
  InputSFTData inputData;		/**< container for input-data and its meta-info */
  SFTtype *outputLFT = NULL;		/**< output 'Long Fourier Transform */
  MultiSFTVector *SSBmultiSFTs = NULL;	/**< SFT vector transferred to the SSB */
  UINT4 X;

  lalDebugLevel = 0;
  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  LAL_CALL ( LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LAL_CALL ( LALInitUserVars(&status, &uvar), &status);

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput(&status, argc, argv), &status);

  if (uvar.help)	/* if help was requested, we're done here */
    exit (0);

  if ( uvar.version ) {
    printf ( "%s\n", lalGitID );
    printf ( "%s\n", lalappsGitID );
    return 0;
  }

  /* ----- Load SFTs */
  LAL_CALL ( LoadInputSFTs(&status, &inputData, &uvar ), &status);

  if ( inputData.numDet != 1 )
    {
      LogPrintf ( LOG_CRITICAL, "Sorry, can only deal with SFTs from single IFO at the moment!\n");
      return -1;
    }

  /* ----- allocate container for SSB-demodulated multi-SFTs */
  {
    UINT4Vector *numSFTs = NULL;

    if ( (numSFTs = XLALCreateUINT4Vector ( inputData.numDet )) == NULL )
      {
	LogPrintf ( LOG_CRITICAL, "Failed to XLALCreateUINT4Vector(%d)!\n", inputData.numDet );
	return -1;
      }
    for ( X=0; X < inputData.numDet; X ++ ) {
      numSFTs->data[X] = inputData.multiSFTs->data[X]->length;		/* number of sfts for IFO X */
    }

    LAL_CALL ( LALCreateMultiSFTVector ( &status, &SSBmultiSFTs, inputData.numBins, numSFTs ), &status );

    XLALDestroyUINT4Vector ( numSFTs );

  } /* allocate SSBmultiSFTs */

  /* ----- Central Demodulation step: time-shift each SFT into SSB */
  /* we do this in the frequency-domain, ie multiply each bin by a complex phase-factor */

  /* loop over detectors X */
  for ( X = 0; X < inputData.numDet; X ++ )
    {
      UINT4 n;
      SFTVector *SSBthisVect = SSBmultiSFTs->data[X];

      /* loop over SFTs n */
      for ( n = 0; n < SSBthisVect->length; n ++ )
	{
	  SFTtype *inputSFT = &(inputData.multiSFTs->data[X]->data[n]);
	  SFTtype *SSBthisSFT = &SSBthisVect->data[n];
	  COMPLEX8Vector *SSBsftData;

	  /* prepare new SFT copy */
	  SSBsftData = SSBthisSFT->data;	/* backup copy of data pointer  */
	  (*SSBthisSFT) = (*inputSFT);				/* struct copy (kills data-pointer) */
	  SSBthisSFT->data = SSBsftData;			/* restore data-pointer */

	  /* copy SFT data */
	  memcpy ( SSBthisSFT->data->data, inputSFT->data->data, SSBthisSFT->data->length * sizeof(SSBthisSFT->data->data[0]) );

	} /* for n < numSFTs */

    } /* for X < numDet */


  /* ----- turn SFT vectors into long Fourier-transforms */
  for ( X=0; X < SSBmultiSFTs->length; X ++ )
    {
      if ( (outputLFT = XLALSFTVectorToLFT ( SSBmultiSFTs->data[X], uvar.upsampling )) == NULL )
	{
	  LogPrintf ( LOG_CRITICAL, "%s: call to XLALSFTVectorToLFT() failed (upsample=%f) ! errno = %d\n",
		      argv[0], uvar.upsampling, xlalErrno );
	  return -1;
	}
    } /* for X < numDet */

  /* write output LFT */
  if ( uvar.outputLFT ) {
    LAL_CALL ( LALWriteSFT2file ( &status, outputLFT, uvar.outputLFT, inputData.dataSummary), &status);
  }

  /* Free config-Variables and userInput stuff */
  LAL_CALL (LALDestroyUserVars (&status), &status);

  LALFree ( inputData.dataSummary );
  LAL_CALL ( LALDestroyMultiSFTVector (&status, &inputData.multiSFTs), &status );
  LAL_CALL ( LALDestroySFTtype ( &status, &outputLFT ), &status );
  LAL_CALL ( LALDestroyMultiSFTVector (&status, &SSBmultiSFTs), &status);

  /* did we forget anything ? */
  LALCheckMemoryLeaks();

  return 0;

} /* main() */

/**
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
void
LALInitUserVars ( LALStatus *status, UserInput_t *uvar )
{
  static const CHAR *fn = "XLALInitUserVars()";

  INITSTATUS( status, fn, rcsid );
  ATTATCHSTATUSPTR (status);

  uvar->help = 0;

  uvar->minStartTime = 0;
  uvar->maxEndTime = LAL_INT4_MAX;

  uvar->upsampling = 1;

  /* register all our user-variables */
  LALregBOOLUserStruct(status,	help, 		'h', UVAR_HELP,     "Print this message");

  LALregSTRINGUserStruct(status,inputSFTs, 	'D', UVAR_OPTIONAL, "File-pattern specifying input SFT-files");
  LALregSTRINGUserStruct(status,outputLFT,      'o', UVAR_OPTIONAL, "Output 'Long Fourier Transform' (LFT) file" );

  LALregINTUserStruct ( status,	minStartTime, 	 0,  UVAR_OPTIONAL, "Earliest SFT-timestamp to include");
  LALregINTUserStruct ( status,	maxEndTime, 	 0,  UVAR_OPTIONAL, "Latest SFT-timestamps to include");

  LALregBOOLUserStruct(status,	version,        'V', UVAR_SPECIAL,   "Output code version");

  LALregREALUserStruct(status,   fmin,		'f', UVAR_OPTIONAL, "Lowest frequency to extract from SFTs. [Default: lowest in inputSFTs]");
  LALregREALUserStruct(status,   fmax,		'F', UVAR_OPTIONAL, "Highest frequency to extract from SFTs. [Default: highest in inputSFTs]");

  LALregREALUserStruct(status, upsampling,	'u', UVAR_OPTIONAL, "Factor to upsample the frequency resolution by");


  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */

/** Initialize code: handle user-input and set everything up. */
void
LoadInputSFTs ( LALStatus *status, InputSFTData *sftData, const UserInput_t *uvar )
{
  static const CHAR *fn = "LoadInputSFTs()";
  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_SFTConstraints;

  LIGOTimeGPS minStartTimeGPS, maxEndTimeGPS;
  MultiSFTVector *multiSFTs = NULL;	    	/* multi-IFO SFT-vectors */
  UINT4 numSFTs;
  REAL8 Tspan, Tdata;

  INITSTATUS( status, fn, rcsid );
  ATTATCHSTATUSPTR (status);

  minStartTimeGPS.gpsSeconds = uvar->minStartTime;
  minStartTimeGPS.gpsNanoSeconds = 0;
  maxEndTimeGPS.gpsSeconds = uvar->maxEndTime;
  maxEndTimeGPS.gpsNanoSeconds = 0;
  constraints.startTime = &minStartTimeGPS;
  constraints.endTime = &maxEndTimeGPS;

  /* ----- get full SFT-catalog of all matching (multi-IFO) SFTs */
  LogPrintf (LOG_DEBUG, "Finding all SFTs to load ... ");

  TRY ( LALSFTdataFind ( status->statusPtr, &catalog, uvar->inputSFTs, &constraints ), status);
  LogPrintfVerbatim (LOG_DEBUG, "done. (found %d SFTs)\n", catalog->length);

  if ( catalog->length == 0 )
    {
      LogPrintf (LOG_CRITICAL, "No matching SFTs for pattern '%s'!\n", uvar->inputSFTs );
      ABORT ( status,  LFTFROMSFTS_EINPUT,  LFTFROMSFTS_MSGEINPUT);
    }

  /* ----- deduce start- and end-time of the observation spanned by the data */
  {
    numSFTs = catalog->length;	/* total number of SFTs */
    sftData->Tsft = 1.0 / catalog->data[0].header.deltaF;
    sftData->startTime = catalog->data[0].header.epoch;
    sftData->endTime   = catalog->data[numSFTs-1].header.epoch;
    LALAddFloatToGPS(status->statusPtr, &sftData->endTime, &sftData->endTime, sftData->Tsft );	/* can't fail */
    Tspan = XLALGPSGetREAL8(&sftData->endTime) - XLALGPSGetREAL8(&sftData->startTime);
    Tdata = numSFTs * sftData->Tsft;
  }

  {/* ----- load the multi-IFO SFT-vectors ----- */
    REAL8 fMin = -1, fMax = -1; /* default: all */

    if ( LALUserVarWasSet ( &uvar->fmin ) )
      fMin = uvar->fmin;
    if ( LALUserVarWasSet ( &uvar->fmax ) )
      fMax = uvar->fmax;

    LogPrintf (LOG_DEBUG, "Loading SFTs ... ");
    TRY ( LALLoadMultiSFTs ( status->statusPtr, &multiSFTs, catalog, fMin, fMax ), status );
    LogPrintfVerbatim (LOG_DEBUG, "done.\n");
    TRY ( LALDestroySFTCatalog ( status->statusPtr, &catalog ), status );

    sftData->fmin = multiSFTs->data[0]->data[0].f0;
    sftData->numBins = multiSFTs->data[0]->data[0].data->length;

    sftData->numDet = multiSFTs->length;
  }

  /* ----- produce a log-string describing the data-specific setup ----- */
  {
    struct tm utc;
    time_t tp;
    CHAR dateStr[512], line[512], summary[1024];
    CHAR *cmdline = NULL;
    UINT4 i;
    tp = time(NULL);
    sprintf (summary, "%%%% Date: %s", asctime( gmtime( &tp ) ) );
    strcat (summary, "%% Loaded SFTs: [ " );
    for ( i=0; i < sftData->numDet; i ++ ) {
      sprintf (line, "%s:%d%s",  multiSFTs->data[i]->data->name, multiSFTs->data[i]->length,
	       (i < sftData->numDet - 1)?", ":" ]\n");
      strcat ( summary, line );
    }
    utc = *XLALGPSToUTC( &utc, (INT4)XLALGPSGetREAL8(&sftData->startTime) );
    strcpy ( dateStr, asctime(&utc) );
    dateStr[ strlen(dateStr) - 1 ] = 0;
    sprintf (line, "%%%% Start GPS time tStart = %12.3f    (%s GMT)\n", XLALGPSGetREAL8(&sftData->startTime), dateStr);
    strcat ( summary, line );
    sprintf (line, "%%%% Total time spanned    = %12.3f s  (%.1f hours)\n", Tspan, Tspan/3600 );
    strcat ( summary, line );

    TRY ( LALUserVarGetLog (status->statusPtr, &cmdline,  UVAR_LOGFMT_CMDLINE ), status);

    if ( (sftData->dataSummary = LALCalloc(1, strlen(summary) + strlen(cmdline) + 20)) == NULL ) {
      ABORT (status, LFTFROMSFTS_EMEM, LFTFROMSFTS_MSGEMEM);
    }
    sprintf( sftData->dataSummary, "\nCommandline: %s\n", cmdline);
    strcat ( sftData->dataSummary, summary );

    LogPrintfVerbatim( LOG_DEBUG, sftData->dataSummary );

    XLALFree ( cmdline );
  } /* write dataSummary string */

  sftData->multiSFTs = multiSFTs;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LoadInputSFTs() */

/** Turn the given multi-IFO SFTvectors into one long Fourier transform (LFT) over the total observation time
 */
SFTtype *
XLALSFTVectorToLFT ( const SFTVector *sfts,	/**< input SFT vector */
		     REAL8 upsampling )		/**< upsampling factor >= 1 */
{
  static const CHAR *fn = "XLALSFTVectorToLFT()";

  COMPLEX8FFTPlan *SFTplan, *LFTplan;

  COMPLEX8TimeSeries *lTS = NULL;		/* long time-series corresponding to full set of SFTs */
  COMPLEX8TimeSeries *sTS = NULL; 		/* short time-series corresponding to a single SFT */

  REAL8 fHet;				/* heterodyning frequency */
  LIGOTimeGPS epoch = {0,0};

  /* constant quantities for all SFTs */
  SFTtype *firstSFT;
  SFTtype *outputLFT;
  UINT4 numBinsSFT;
  REAL8 dfSFT;
  REAL8 Tsft;

  REAL8 deltaT;

  UINT4 numSFTs;
  UINT4 n;
  REAL8 Tspan;
  UINT4 numTimeSamples;
  UINT4 NnegSFT, NnegLFT;
  REAL8 f0SFT, f0LFT;
  UINT4 numSFTsFit;

  if ( !sfts || (sfts->length == 0) )
    {
      XLALPrintError ("%s: empty SFT input!\n", fn );
      XLAL_ERROR_NULL (fn, XLAL_EINVAL);
    }
  if ( upsampling < 1 )
    {
      XLALPrintError ("%s: upsampling factor (%f) must be >= 1 \n", fn, upsampling );
      XLAL_ERROR_NULL (fn, XLAL_EINVAL);
    }

  /* some useful shorthands */
  numSFTs = sfts->length;
  firstSFT = &(sfts->data[0]);
  numBinsSFT = firstSFT->data->length;
  dfSFT = firstSFT->deltaF;
  Tsft = 1.0 / dfSFT;
  deltaT = 1.0 / ( numBinsSFT * dfSFT );

  f0SFT = firstSFT->f0;

  /* ---------- determine time-span of the final long time-series */
  Tspan = XLALGPSDiff ( &sfts->data[numSFTs-1].epoch, &sfts->data[0].epoch ) + Tsft;
  Tspan *= upsampling;

  /* NOTE: Tspan MUST be an integer multiple of Tsft,
   * in order for the frequency bins of the final FFT
   * to be commensurate with the SFT bins.
   * This is required so that fHet is an exact
   * frequency-bin in both cases
   */
  numSFTsFit = (UINT4)floor(Tspan / Tsft + 0.5);	/* round */
  Tspan = numSFTsFit * Tsft;
  numTimeSamples = (UINT4)floor(Tspan / deltaT + 0.5);	/* round */

  /* determine the heterodyning frequency */
  /* fHet = DC of our internal DFTs */
  NnegSFT = NhalfNeg ( numBinsSFT );
  fHet = f0SFT + 1.0 * NnegSFT * dfSFT;

  /* translate this back into fmin for the LFT (counting down from DC==fHet) */
  NnegLFT = NhalfNeg ( numTimeSamples );
  f0LFT = fHet - NnegLFT / Tspan;

  printf ("NSFT = %d, NLFT = %d, fminSFT = %.9f, fHet = %.9f, fminLFT = %.9f\n",
	  numBinsSFT, numTimeSamples, f0SFT, fHet, f0LFT );


  /* ----- Prepare invFFT: compute plan for FFTW */
  if ( (SFTplan = XLALCreateReverseCOMPLEX8FFTPlan( numBinsSFT, 0 )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateReverseCOMPLEX8FFTPlan(%d, ESTIMATE ) failed! errno = %d!\n", fn, numBinsSFT, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

  /* ----- prepare long TimeSeries container ---------- */
  if ( (lTS = XLALCreateCOMPLEX8TimeSeries ( firstSFT->name, &firstSFT->epoch, 0, deltaT, &empty_LALUnit, numTimeSamples )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numTimeSamples, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }
  memset ( lTS->data->data, 0, numTimeSamples * sizeof(*lTS->data->data)); /* set all time-samples to zero */

  /* ----- Prepare short time-series holding ONE invFFT of a single SFT */
  if ( (sTS = XLALCreateCOMPLEX8TimeSeries ( "short timeseries", &epoch, 0, deltaT, &empty_LALUnit, numBinsSFT )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numBinsSFT, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

  /* ----- prepare output LFT ---------- */
  if ( (outputLFT = LALCalloc ( 1, sizeof(*outputLFT) )) == NULL )
    {
      XLALPrintError ( "%s: LALCalloc ( 1, %d ) failed!\n", fn, sizeof(*outputLFT) );
      XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
    }

  /* prepare LFT header */
  strcpy ( outputLFT->name, firstSFT->name );
  strcat ( outputLFT->name, ":long Fourier transform");
  outputLFT->epoch  = firstSFT->epoch;
  outputLFT->f0     = f0LFT;
  outputLFT->deltaF = 1.0 / Tspan;
  outputLFT->sampleUnits = firstSFT->sampleUnits;

  if ( (outputLFT->data = XLALCreateCOMPLEX8Vector ( numTimeSamples )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateCOMPLEX8Vector(%d) failed!\n", fn, numTimeSamples );
      XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
    }

  /* ---------- loop over all SFTs and inverse-FFT them ---------- */
  for ( n = 0; n < numSFTs; n ++ )
    {
      SFTtype *thisSFT = &(sfts->data[n]);
      UINT4 bin0_n;
      REAL8 offset_n, nudge_n;
      UINT4 copyLen, binsLeft;

      UINT4 k;
      REAL8 offset0, offsetEff, hetCycles;
      REAL4 hetCorr_re, hetCorr_im;
      REAL4 fact_re, fact_im;
      REAL4 norm = 1.0 / numBinsSFT;


      /* find bin in long timeseries corresponding to starttime of *this* SFT */
      offset_n = XLALGPSDiff ( &thisSFT->epoch, &firstSFT->epoch );
      bin0_n = (UINT4) ( offset_n / deltaT + 0.5 );	/* round to closest bin */

      nudge_n = bin0_n * deltaT - offset_n;		/* rounding error */
      nudge_n = 1e-9 * (floor)(nudge_n * 1e9 + 0.5);	/* round to closest nanosecond */
      {
	REAL8 t0 = XLALGPSGetREAL8 ( &firstSFT->epoch );
	printf ("n = %d: t0_n = %f, sft_tn =(%d,%d), bin-offset = %g s, corresponding to %g timesteps\n",
		n, t0 + bin0_n * deltaT, sfts->data[n].epoch.gpsSeconds,  sfts->data[n].epoch.gpsNanoSeconds, nudge_n, nudge_n/deltaT );
      }

      if ( (nudge_n != 0) && (XLALTimeShiftSFT ( thisSFT, nudge_n ) != XLAL_SUCCESS) )
	{
	  XLALPrintError ( "%s: XLALTimeShiftSFT(sft-%d, %g) failed! errno = %d!\n", fn, n, nudge_n, xlalErrno );
	  XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
	}


      /* determine heterodyning phase-correction for this SFT */
      offset0 = XLALGPSDiff ( &thisSFT->epoch, &firstSFT->epoch );
      /* we know that fHet * Tsft = integer, so we only need the remainder of offset_t0 % Tsft */
      offsetEff = fmod ( offset0, Tsft );
      offsetEff = 1e-9 * (floor)( offsetEff * 1e9 + 0.5 );	/* round to closest integer multiple of nanoseconds */

      hetCycles = fmod ( fHet * offsetEff, 1);	/* required heterodyning phase-correction for this SFT */

      sin_cos_2PI_LUT (&hetCorr_im, &hetCorr_re, -hetCycles );

      /* apply all phase- and normalizing factors */
      fact_re = norm * hetCorr_re;
      fact_im = norm * hetCorr_im;

      printf ("SFT n = %d: (tn - t0) = %g s EQUIV %g s, hetCycles = %g, ==> fact = %g + i %g\n", n, offset0, offsetEff, hetCycles, fact_re, fact_im );

      for ( k = 0; k < numBinsSFT; k ++ )
	{
	  REAL8 binReal, binImag;

	  binReal = fact_re * thisSFT->data->data[k].re - fact_im * thisSFT->data->data[k].im;
	  binImag = fact_re * thisSFT->data->data[k].im + fact_im * thisSFT->data->data[k].re;

	  thisSFT->data->data[k].re = binReal;
	  thisSFT->data->data[k].im = binImag;
	} /* k < numBins */

      if ( XLALReorderSFTtoFFTW (thisSFT->data) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALReorderSFTtoFFTW() failed! errno = %d!\n", fn, xlalErrno );
	  XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
	}

      if ( XLALCOMPLEX8VectorFFT( sTS->data, thisSFT->data, SFTplan ) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALCOMPLEX8VectorFFT() failed! errno = %d!\n", fn, xlalErrno );
	  XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
	}

      /* copy short (shifted) timeseries into correct location within long timeseries */
      binsLeft = numTimeSamples - bin0_n - 1;
      copyLen = MYMIN ( numBinsSFT, binsLeft );		/* make sure not to write past the end of the long TS */
      memcpy ( &lTS->data->data[bin0_n], sTS->data->data, copyLen * sizeof(lTS->data->data[0]) );

    } /* for n < numSFTs */

  /* ---------- now FFT the long timeseries ---------- */

  /* ----- compute plan for FFTW */
  if ( (LFTplan = XLALCreateForwardCOMPLEX8FFTPlan( numTimeSamples, 0 )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateForwardCOMPLEX8FFTPlan(%d, ESTIMATE ) failed! errno = %d!\n", fn, numTimeSamples, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

  if ( XLALCOMPLEX8VectorFFT( outputLFT->data, lTS->data, LFTplan ) != XLAL_SUCCESS )
    {
      XLALPrintError ( "%s: XLALCOMPLEX8VectorFFT() failed! errno = %d!\n", fn, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

  if ( XLALReorderFFTWtoSFT (outputLFT->data) != XLAL_SUCCESS )
    {
      XLALPrintError ( "%s: XLALReorderFFTWtoSFT() failed! errno = %d!\n", fn, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }


  /* ---------- debug-output timeseries */
  {
    REAL8 t0 = XLALGPSGetREAL8 ( &lTS->epoch );
    REAL8 dt = lTS->deltaT;
    UINT4 i;
    for ( i = 0; i < lTS->data->length; i ++ )
      fprintf ( stderr, "%f   %f %f\n", t0 + i * dt, lTS->data->data[i].re, lTS->data->data[i].im);
  }


  /* cleanup memory */
  XLALDestroyCOMPLEX8TimeSeries ( sTS );
  XLALDestroyCOMPLEX8TimeSeries ( lTS );
  XLALDestroyCOMPLEX8FFTPlan ( SFTplan );
  XLALDestroyCOMPLEX8FFTPlan ( LFTplan );

  return outputLFT;

} /* XLALSFTVectorToLFT() */


/** Change frequency-bin order from fftw-convention to a 'SFT'
 * ie. from FFTW: f[0], f[1],...f[N/2], f[-(N-1)/2], ... f[-2], f[-1]
 * to: f[-(N-1)/2], ... f[-1], f[0], f[1], .... f[N/2]
 */
int
XLALReorderFFTWtoSFT (COMPLEX8Vector *X)
{
  static const CHAR *fn = "XLALReorderFFTWtoSFT()";
  UINT4 N, Npos_and_DC, Nneg;

  /* temporary storage for data */
  COMPLEX8 *tmp;


  if ( !X || (X->length==0) )
    {
      XLALPrintError ("%s: empty input vector 'X'!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }

  N = X -> length;
  Npos_and_DC = NhalfPosDC ( N );
  Nneg = NhalfNeg ( N );

  /* allocate temporary storage for swap */
  if ( (tmp = XLALMalloc ( N * sizeof(*tmp) )) == NULL )
    {
      XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, N * sizeof(*tmp));
      XLAL_ERROR (fn, XLAL_ENOMEM);
    }

  memcpy ( tmp, X->data, N * sizeof(*tmp) );

  /* Copy second half from FFTW: 'negative' frequencies */
  memcpy ( X->data, tmp + Npos_and_DC, Nneg * sizeof(*tmp) );

  /* Copy first half from FFTW: 'DC + positive' frequencies */
  memcpy ( X->data + Nneg, tmp, Npos_and_DC * sizeof(*tmp) );

  XLALFree(tmp);

  return XLAL_SUCCESS;

}  /* XLALReorderFFTWtoSFT() */

/** Change frequency-bin order from 'SFT' to fftw-convention
 * ie. from f[-(N-1)/2], ... f[-1], f[0], f[1], .... f[N/2]
 * to FFTW: f[0], f[1],...f[N/2], f[-(N-1)/2], ... f[-2], f[-1]
 */
int
XLALReorderSFTtoFFTW (COMPLEX8Vector *X)
{
  static const CHAR *fn = "XLALReorderSFTtoFFTW()";
  UINT4 N, Npos_and_DC, Nneg;

  /* temporary storage for data */
  COMPLEX8 *tmp;


  if ( !X || (X->length==0) )
    {
      XLALPrintError ("%s: empty input vector 'X'!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }

  N = X->length;
  Npos_and_DC = NhalfPosDC ( N );
  Nneg = NhalfNeg ( N );

  /* allocate temporary storage for swap */
  if ( (tmp = XLALMalloc ( N * sizeof(*tmp) )) == NULL )
    {
      XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, N * sizeof(*tmp));
      XLAL_ERROR (fn, XLAL_ENOMEM);
    }

  memcpy ( tmp, X->data, N * sizeof(*tmp) );

  /* Copy second half of FFTW: 'negative' frequencies */
  memcpy ( X->data + Npos_and_DC, tmp, Nneg * sizeof(*tmp) );

  /* Copy first half of FFTW: 'DC + positive' frequencies */
  memcpy ( X->data, tmp + Nneg, Npos_and_DC * sizeof(*tmp) );

  XLALFree(tmp);

  return XLAL_SUCCESS;

}  /* XLALReorderSFTtoFFTW() */


/** Time-shift the given SFT by an amount of 'shift' seconds,
 * using the frequency-domain expression y(f) = x(f) * e^(-i 2pi f tau),
 * which shifts x(t) into y(t) = x(t - tau)
 */
int
XLALTimeShiftSFT ( SFTtype *sft, REAL8 shift )
{
  const CHAR *fn = "XLALTimeShiftSFT()";
  UINT4 k;

  if ( !sft || !sft->data )
    {
      XLALPrintError ("%s: empty input SFT!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }

  for ( k=0; k < sft->data->length; k++)
    {
      REAL8 fk = sft->f0 + k * sft->deltaF;	/* frequency of k-th bin */
      REAL8 shiftCyles = shift * fk;
      REAL4 fact_re, fact_im;			/* complex phase-shift factor e^(-2pi f tau) */
      REAL4 yRe, yIm;

      sin_cos_2PI_LUT ( &fact_im, &fact_re, shiftCyles );

      yRe = fact_re * sft->data->data[k].re - fact_im * sft->data->data[k].im;
      yIm = fact_re * sft->data->data[k].im + fact_im * sft->data->data[k].re;

      sft->data->data[k].re = yRe;
      sft->data->data[k].im = yIm;

    } /* for k < numBins */

  /* adjust SFTs epoch to the shift */
  XLALGPSAdd( &sft->epoch, shift );

  return XLAL_SUCCESS;

} /* XLALTimeShiftSFT() */
