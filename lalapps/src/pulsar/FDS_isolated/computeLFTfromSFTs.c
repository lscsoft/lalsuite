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

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  MultiSFTVector *multiSFTs;	/**< input SFT vector */
  CHAR *dataSummary;            /**< descriptive string describing the data */
  REAL8 Tsft;			/**< duration of each SFT */
  LIGOTimeGPS startTime;	/**< start-time of SFTs */
  LIGOTimeGPS endTime;
  REAL8 fmin;			/**< smallest frequency contained in input SFTs */
  UINT4 numBins;		/**< number of frequency bins in input SFTs */
} ConfigVariables;


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

  BOOLEAN version;	/**< output version-info */
} UserInput_t;

static UserInput_t empty_UserInput;
static LALUnit empty_LALUnit;

/* ---------- local prototypes ---------- */
int main(int argc, char *argv[]);

void LALInitUserVars ( LALStatus *status, UserInput_t *uvar);
void LoadInputSFTs ( LALStatus *status, ConfigVariables *cfg, const UserInput_t *uvar );
int XLALWeightedSumOverSFTVector ( SFTVector **outSFT, const SFTVector *inVect, const COMPLEX8Vector *weights );

SFTtype *XLALSFTVectorToLFT ( const SFTVector *sfts );

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
  ConfigVariables GV;			/**< container for various derived configuration settings */
  SFTtype *outputLFT = NULL;		/**< output 'Long Fourier Transform */
  UINT4 X, numDet;

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

  /* Load SFTs */
  LAL_CALL ( LoadInputSFTs(&status, &GV, &uvar ), &status);

  numDet = GV.multiSFTs->length;

  if ( numDet != 1 )
    {
      LogPrintf ( LOG_CRITICAL, "Sorry, can only deal with SFTs from single IFO at the moment!\n");
      return -1;
    }

  for ( X=0; X < numDet; X ++ )
    {

      if ( (outputLFT = XLALSFTVectorToLFT ( GV.multiSFTs->data[X] )) == NULL )
	{
	  LogPrintf ( LOG_CRITICAL, "%s: call to XLALSFTVectorToLFT() failed! errno = %d\n", argv[0], xlalErrno );
	  return -1;
	}
    } /* for X < numDet */

  /* write output LFT */
  if ( uvar.outputLFT ) {
    LAL_CALL ( LALWriteSFT2file ( &status, outputLFT, uvar.outputLFT, GV.dataSummary), &status);
  }

  /* Free config-Variables and userInput stuff */
  LAL_CALL (LALDestroyUserVars (&status), &status);

  LALFree ( GV.dataSummary );
  LAL_CALL ( LALDestroyMultiSFTVector (&status, &GV.multiSFTs), &status );
  LAL_CALL ( LALDestroySFTtype ( &status, &outputLFT ), &status );

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
  const CHAR *fn = "XLALInitUserVars()";

  INITSTATUS( status, fn, rcsid );
  ATTATCHSTATUSPTR (status);

  uvar->help = 0;

  uvar->minStartTime = 0;
  uvar->maxEndTime = LAL_INT4_MAX;


  /* register all our user-variables */
  LALregBOOLUserStruct(status,	help, 		'h', UVAR_HELP,     "Print this message");

  LALregSTRINGUserStruct(status,inputSFTs, 	'D', UVAR_OPTIONAL, "File-pattern specifying input SFT-files");
  LALregSTRINGUserStruct(status,outputLFT,      'o', UVAR_OPTIONAL, "Output 'Long Fourier Transform' (LFT) file" );

  LALregINTUserStruct ( status,	minStartTime, 	 0,  UVAR_OPTIONAL, "Earliest SFT-timestamp to include");
  LALregINTUserStruct ( status,	maxEndTime, 	 0,  UVAR_OPTIONAL, "Latest SFT-timestamps to include");

  LALregBOOLUserStruct(status,	version,        'V', UVAR_SPECIAL,   "Output code version");

  LALregREALUserStruct(status,   fmin,		'f', UVAR_OPTIONAL, "Lowest frequency to extract from SFTs. [Default: lowest in inputSFTs]");
  LALregREALUserStruct(status,   fmax,		'F', UVAR_OPTIONAL, "Highest frequency to extract from SFTs. [Default: highest in inputSFTs]");

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */

/** Initialize code: handle user-input and set everything up. */
void
LoadInputSFTs ( LALStatus *status, ConfigVariables *cfg, const UserInput_t *uvar )
{
  const CHAR *fn = "LoadInputSFTs()";
  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_SFTConstraints;

  LIGOTimeGPS minStartTimeGPS, maxEndTimeGPS;
  MultiSFTVector *multiSFTs = NULL;	    	/* multi-IFO SFT-vectors */
  UINT4 numDet, numSFTs;
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
    cfg->Tsft = 1.0 / catalog->data[0].header.deltaF;
    cfg->startTime = catalog->data[0].header.epoch;
    cfg->endTime   = catalog->data[numSFTs-1].header.epoch;
    LALAddFloatToGPS(status->statusPtr, &cfg->endTime, &cfg->endTime, cfg->Tsft );	/* can't fail */
    Tspan = XLALGPSGetREAL8(&cfg->endTime) - XLALGPSGetREAL8(&cfg->startTime);
    Tdata = numSFTs * cfg->Tsft;
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

    cfg->fmin = multiSFTs->data[0]->data[0].f0;
    cfg->numBins = multiSFTs->data[0]->data[0].data->length;

    numDet = multiSFTs->length;
  }

  /* ----- produce a log-string describing the data-specific setup ----- */
  {
    struct tm utc;
    time_t tp;
    CHAR dateStr[512], line[512], summary[1024];
    UINT4 i;
    tp = time(NULL);
    sprintf (summary, "%%%% Date: %s", asctime( gmtime( &tp ) ) );
    strcat (summary, "%% Loaded SFTs: [ " );
    for ( i=0; i < numDet; i ++ ) {
      sprintf (line, "%s:%d%s",  multiSFTs->data[i]->data->name, multiSFTs->data[i]->length,
	       (i < numDet - 1)?", ":" ]\n");
      strcat ( summary, line );
    }
    utc = *XLALGPSToUTC( &utc, (INT4)XLALGPSGetREAL8(&cfg->startTime) );
    strcpy ( dateStr, asctime(&utc) );
    dateStr[ strlen(dateStr) - 1 ] = 0;
    sprintf (line, "%%%% Start GPS time tStart = %12.3f    (%s GMT)\n", XLALGPSGetREAL8(&cfg->startTime), dateStr);
    strcat ( summary, line );
    sprintf (line, "%%%% Total time spanned    = %12.3f s  (%.1f hours)\n", Tspan, Tspan/3600 );
    strcat ( summary, line );

    if ( (cfg->dataSummary = LALCalloc(1, strlen(summary) + 1 )) == NULL ) {
      ABORT (status, LFTFROMSFTS_EMEM, LFTFROMSFTS_MSGEMEM);
    }
    strcpy ( cfg->dataSummary, summary );

    LogPrintfVerbatim( LOG_DEBUG, cfg->dataSummary );
  } /* write dataSummary string */

  cfg->multiSFTs = multiSFTs;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LoadInputSFTs() */

/** Turn the given multi-IFO SFTvectors into COMPLEX8 (heterodyned) timeseries, one per IFO
 */
SFTtype *
XLALSFTVectorToLFT ( const SFTVector *sfts )
{
  const CHAR *fn = "XLALSFTVectorToTimeseries()";

  COMPLEX8FFTPlan *SFTplan, *LFTplan;
  LIGOTimeGPS epoch = {0,0};

  COMPLEX8TimeSeries *lTS = NULL;		/* long time-series corresponding to full set of SFTs */
  COMPLEX8TimeSeries *sTS = NULL; 		/* short time-series corresponding to a single SFT */

  /* constant quantities for all SFTs */
  SFTtype *firstHead;
  SFTtype *outputLFT;
  UINT4 numBins;
  REAL8 deltaF;
  REAL8 Tsft;

  REAL8 f0;
  REAL8 deltaT;

  UINT4 numSFTs;
  REAL8 startTime;
  REAL8 endTime;
  UINT4 n;
  REAL8 Tspan;
  REAL8 numTimeSamples;

  if ( !sfts || (sfts->length == 0) )
    {
      XLALPrintError ("%s: empty SFT input!\n", fn );
      XLAL_ERROR_NULL (fn, XLAL_EINVAL);
    }

  /* obtain quantities that are constant for all SFTs */
  numSFTs = sfts->length;
  firstHead = &(sfts->data[0]);
  numBins = firstHead->data->length;
  deltaF = firstHead->deltaF;
  Tsft = 1.0 / deltaF;

  f0 = firstHead->f0;
  deltaT = 1.0 / ( numBins * deltaF );

  startTime = XLALGPSGetREAL8 ( &sfts->data[0].epoch );
  endTime   = XLALGPSGetREAL8 ( &sfts->data[numSFTs-1].epoch ) + Tsft;

  Tspan = endTime - startTime;
  numTimeSamples = (UINT4)(Tspan / deltaT + 0.5);	/* round */
  Tspan = numTimeSamples * deltaT;			/* correct Tspan to actual time-samples */

  /* ----- Prepare invFFT: compute plan for FFTW */
  if ( (SFTplan = XLALCreateReverseCOMPLEX8FFTPlan( numBins, 0 )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateReverseCOMPLEX8FFTPlan(%d, ESTIMATE ) failed! errno = %d!\n", fn, numBins, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

  /* ----- prepare long TimeSeries container ---------- */
  if ( (lTS = XLALCreateCOMPLEX8TimeSeries ( firstHead->name, &firstHead->epoch, f0, deltaT, &empty_LALUnit, numTimeSamples )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numTimeSamples, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }
  memset ( lTS->data->data, 0, numTimeSamples * sizeof(*lTS->data->data)); /* set all time-samples to zero */

  /* ----- Prepare short time-series holding ONE invFFT of a single SFT */
  if ( (sTS = XLALCreateCOMPLEX8TimeSeries ( "short timeseries", &epoch, f0, deltaT, &empty_LALUnit, numBins )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numBins, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

  /* ----- prepare output LFT ---------- */
  if ( (outputLFT = LALCalloc ( 1, sizeof(*outputLFT) )) == NULL )
    {
      XLALPrintError ( "%s: LALCalloc ( 1, %d ) failed!\n", fn, sizeof(*outputLFT) );
      XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
    }

  strcpy ( outputLFT->name, firstHead->name );
  strcat ( outputLFT->name, ":long Fourier transform");
  outputLFT->epoch  = firstHead->epoch;
  outputLFT->f0     = firstHead->f0;
  outputLFT->deltaF = 1.0 / Tspan;
  outputLFT->sampleUnits = firstHead->sampleUnits;

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
      REAL8 startTime_n;

      if ( XLALCOMPLEX8VectorFFT( sTS->data, thisSFT->data, SFTplan ) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALCOMPLEX8VectorFFT() failed! errno = %d!\n", fn, xlalErrno );
	  XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
	}

      /* find bin in long timeseries corresponding to starttime of *this* SFT */
      startTime_n = XLALGPSGetREAL8 ( &sfts->data[n].epoch );
      bin0_n = (UINT4) ( ( startTime_n - startTime ) / deltaT + 0.5 );	/* round to closest bin */

      /* copy short timeseries into correct location within long timeseries */
      memcpy ( &lTS->data->data[bin0_n], sTS->data->data, numBins * sizeof(lTS->data->data[bin0_n]) );

    } /* for n < numSFTs */


  /* ---------- now FFT the complete timeseries ---------- */

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

  /* cleanup memory */
  XLALDestroyCOMPLEX8TimeSeries ( sTS );
  XLALDestroyCOMPLEX8TimeSeries ( lTS );
  XLALDestroyCOMPLEX8FFTPlan ( SFTplan );
  XLALDestroyCOMPLEX8FFTPlan ( LFTplan );

  return outputLFT;

} /* XLALMultiSFTsToTimeseries() */
