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
 * \file CFSv3
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
#include "Fstat_v3.h"


#include <lal/lalGitID.h>
#include <lalappsGitID.h>

#include <lalapps.h>

/* local includes */

/*---------- DEFINES ----------*/
#define LAL_INT4_MAX    LAL_UINT8_C(2147483647)

#define CFSV3_ERROR   1
#define CFSV3_ERROR_MSG  "CFSv3 failed"

/* ---------- local types ---------- */
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

/* ---------- local prototypes ---------- */
int main(int argc, char *argv[]);

void LALInitUserVars ( LALStatus *status, UserInput_t *uvar);
void LoadInputSFTs ( LALStatus *status, InputSFTData *data, const UserInput_t *uvar );

/*---------- empty initializers ---------- */
RCSID( "$Id$");

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
	  (*SSBthisSFT) = (*inputSFT);		/* struct copy (kills data-pointer) */
	  SSBthisSFT->data = SSBsftData;	/* restore data-pointer */

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
      ABORT ( status,  CFSV3_ERROR, CFSV3_ERROR_MSG );
    }

  /* ----- deduce start- and end-time of the observation spanned by the data */
  {
    numSFTs = catalog->length;	/* total number of SFTs */
    sftData->Tsft = 1.0 / catalog->data[0].header.deltaF;
    sftData->startTime = catalog->data[0].header.epoch;
    sftData->endTime   = catalog->data[numSFTs-1].header.epoch;
    XLALGPSAdd( &sftData->endTime, sftData->Tsft );
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
      ABORT ( status,  CFSV3_ERROR, CFSV3_ERROR_MSG );
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
