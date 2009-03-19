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

#include <fftw3.h>

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LogPrintf.h>
#include <lal/TimeSeries.h>
#include <lal/RealFFT.h>

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

/* ---------- local prototypes ---------- */
int main(int argc, char *argv[]);

void LALInitUserVars ( LALStatus *status, UserInput_t *uvar);
void LoadInputSFTs ( LALStatus *status, ConfigVariables *cfg, const UserInput_t *uvar );
int XLALWeightedSumOverSFTVector ( SFTVector **outSFT, const SFTVector *inVect, const COMPLEX8Vector *weights );
MultiREAL4TimeSeries *XLALMultiSFTsToTimeseries ( const MultiSFTVector *multiSFTs );

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
  SFTVector *outputLFT = NULL;		/**< output 'Long Fourier Transform */
  UINT4 j, LFTnumBins;			/**< number of frequency bins in output LFT */
  MultiREAL4TimeSeries *multiTS = NULL;

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


  if ( (multiTS = XLALMultiSFTsToTimeseries ( GV.multiSFTs )) == NULL )
    {
      LogPrintf ( LOG_CRITICAL, "%s: call to XLALMultiSFTsToTimeseries() failed! errno = %d\n", argv[0], xlalErrno );
      return -1;
    }


  /* write output LFT */
  if ( uvar.outputLFT ) {
    LAL_CALL ( LALWriteSFT2file ( &status, &outputLFT->data[0], uvar.outputLFT, GV.dataSummary), &status);
  }

  /* Free config-Variables and userInput stuff */
  LAL_CALL (LALDestroyUserVars (&status), &status);

  LALFree ( GV.dataSummary );
  LAL_CALL ( LALDestroyMultiSFTVector (&status, &GV.multiSFTs), &status);


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

/** Turn the given multi-IFO SFTvectors into REAL4 timeseries, one per IFO
 */
MultiREAL4TimeSeries *
XLALMultiSFTsToTimeseries ( const MultiSFTVector *multiSFTs )
{
  const CHAR *fn = "XLALMultiSFTsToTimeseries()";

  MultiREAL4TimeSeries *mts = NULL;
  UINT4 X, numDet;
  REAL4FFTPlan *plan;

  UINT4 numBinsSFT;	/* identical for all SFTs ! */
  REAL8 f0;
  REAL8 deltaF;
  LIGOTimeGPS epoch = {0,0};
  REAL8 deltaT;

  REAL4TimeSeries *sTS = NULL; 		/* short time-series corresponding to a single SFT */

  if ( !multiSFTs || (multiSFTs->length == 0) )
    {
      XLALPrintError ("%s: empty SFT input!\n", fn );
      XLAL_ERROR_NULL (fn, XLAL_EINVAL);
    }

  numDet = multiSFTs->length;


  /* ---------- first create the multi-IFO TimeSeries vector ---------- */
  if ( (mts = XLALMalloc ( sizeof(*mts) )) == NULL )
    {
      XLALPrintError ("%s: XLALMalloc (%s) failed!\n", fn, sizeof(*mts) );
      XLAL_ERROR_NULL (fn, XLAL_ENOMEM);
    }

  mts->length = numDet;
  for ( X=0; X < numDet; X ++)
    {
      SFTtype *thisHead = &(multiSFTs->data[X]->data[0]);
      UINT4 numSFTs = multiSFTs->data[X]->length;

      REAL8 FreqBand = f0 + 2.0 * numBinsSFT * thisHead->deltaF;	/* count positive and negative frequencies */
      REAL8 Tsft = 1.0 / deltaF;
      REAL8 startTime = XLALGPSGetREAL8 ( &multiSFTs->data[X]->data[0].epoch );
      REAL8 endTime   = XLALGPSGetREAL8 ( &multiSFTs->data[X]->data[numSFTs-1].epoch ) + Tsft;
      REAL8 Tspan = endTime - startTime;
      REAL8 numBinsTS = floor ( Tspan / deltaT );

      numBinsSFT = thisHead->data->length;
      f0 = thisHead->f0;
      deltaF = thisHead->deltaF;
      deltaT  = 0.5 / FreqBand;	/* sampling-frequency = 2*Band */

      mts->data[X] = XLALCreateREAL4TimeSeries ( thisHead->name, &thisHead->epoch, f0, deltaT, NULL, numBinsTS );
      if ( mts->data[X] == NULL )
	{
	  XLALPrintError ("%s: XLALCreateREAL4TimeSeries() failed to create timseries with %d bins!\n", fn, numBinsTS );
	  XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
	}

      /* set timeseries to zeros */
      memset( mts->data[X]->data->data, 0, numBinsTS * sizeof(REAL4) );

    } /* for X < numDet */


  /* Prepare invFFT: compute plan for FFTW */
  if ( (plan = XLALCreateReverseREAL4FFTPlan( numBinsSFT, FFTW_ESTIMATE )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateReverseREAL4FFTPlan(%d, ESTIMATE ) failed! errno = %d!\n", fn, numBinsSFT, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

  /* Prepare short time-series holding the invFFT of a single SFT */
  if ( (sTS = XLALCreateREAL4TimeSeries ( "short timeseries", &epoch, f0, deltaT, NULL, numBinsSFT )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateREAL4TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numBinsSFT, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }


  /* cleanup memory */
  XLALDestroyREAL4TimeSeries ( sTS );
  XLALDestroyREAL4FFTPlan ( plan );

  return mts;

} /* XLALMultiSFTsToTimeseries() */
