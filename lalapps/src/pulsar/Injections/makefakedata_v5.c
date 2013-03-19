/*
*  Copyright (C) 2013 Reinhard Prix
*  Copyright (C) 2008, 2010 Karl Wette
*  Copyright (C) 2008 Chris Messenger
*  Copyright (C) 2007 Badri Krishnan, Reinhard Prix
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
 * \author R. Prix, M.A. Papa, X. Siemens, B. Allen, C. Messenger
 */

/*-----------------------------------------------------------------------
 *
 * File Name: makefakedata_v5.c
 *
 * Authors: R. Prix, M.A. Papa, X. Siemens, B. Allen, C. Messenger
 *
 * This code is a descendant of an earlier implementation 'makefakedata_v4.c',
 * which itself descended from 'makefakedata_v2.c'
 * by Badri Krishnan, Bruce Allen, Maria Alessandra Papa, Reinhard Prix, Xavier Siemens, Yousuke Itoh
 *
 *-----------------------------------------------------------------------
 */

/* ---------- includes ---------- */
#include <sys/stat.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lalapps.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FrequencySeries.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Random.h>
#include <gsl/gsl_math.h>

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/SimulatePulsarSignal.h>
#include <lal/TimeSeries.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/Window.h>

#include <lal/TransientCW_utils.h>
#include <lal/CWMakeFakeData.h>

#include <lalapps.h>

/***************************************************/

/*----------------------------------------------------------------------*/


/** configuration-variables derived from user-variables */
typedef struct
{
  PulsarParamsVector *injectionSources;		/**< list of injection source parameters */

  EphemerisData *edat;		/**< ephemeris-data */
  MultiLIGOTimeGPSVector *multiTimestamps;/**< a vector of timestamps to generate time-series/SFTs for */
  MultiDetectorInfo detInfo;	//!< detectors and noise-floors (for Gaussian noise) to generate data for

  SFTCatalog *noiseCatalog; 			/**< catalog of noise-SFTs */
  MultiSFTCatalogView *multiNoiseCatalogView; 	/**< multi-IFO 'view' of noise-SFT catalogs */

  transientWindow_t transientWindow;	/**< properties of transient-signal window */
  CHAR *VCSInfoString;          /**< LAL + LALapps Git version string */
} ConfigVars_t;

// ----- User variables
typedef struct
{
  BOOLEAN help;		/**< Print this help/usage message */

  /* output */
  CHAR *outSFTdir;		/**< Output directory for SFTs */
  CHAR *outMiscField;		/**< 'misc' field entry in the SFT-filename (part of 'description' field) */

  BOOLEAN outSingleSFT;	        /**< use to output a single concatenated SFT */

  CHAR *TDDfile;		/**< Filename for ASCII output time-series */
  CHAR *logfile;		/**< name of logfile */

  /* specify start + duration */
  LALStringVector *timestampsFiles;        /**< Names of numDet timestamps files */
  INT4 startTime;		/**< Start-time of requested signal in detector-frame (GPS seconds) */
  INT4 duration;		/**< Duration of requested signal in seconds */

  /* time-series sampling + heterodyning frequencies */
  REAL8 fmin;		/**< Lowest frequency in output SFT (= heterodyning frequency) */
  REAL8 Band;		/**< bandwidth of output SFT in Hz (= 1/2 sampling frequency) */

  /* SFT params */
  REAL8 Tsft;		        /**< SFT time baseline Tsft */
  REAL8 SFToverlap;		/**< overlap SFTs by this many seconds */

  LALStringVector* IFOs;	/**< list of detector-names "H1,H2,L1,.." or single detector*/

  /* noise to add [OPTIONAL] */
  LALStringVector* sqrtSX; 	/**< Add Gaussian noise: list of respective detectors' noise-floors sqrt{Sn}" */

  CHAR *noiseSFTs;		/**< Glob-like pattern specifying noise-SFTs to be added to signal */

  /* Window function [OPTIONAL] */
  CHAR *SFTWindowType;		/**< Windowing function to apply to the SFT time series */
  REAL8 SFTWindowBeta;         	/**< 'beta' parameter required for certain window-types */

  CHAR *ephemDir;		/**< Directory path for ephemeris files (optional), use LAL_DATA_PATH if unset. */
  CHAR *ephemYear;		/**< Year (or range of years) of ephemeris files to be used */

  /* pulsar parameters [REQUIRED] */
  CHAR *injectionSources;	///< either a file-specification ("@file-pattern") or a config-string defining the sources to inject

  BOOLEAN version;		/**< output version information */

  INT4 randSeed;		/**< allow user to specify random-number seed for reproducible noise-realizations */

} UserVariables_t;


// ---------- exportable API types ----------

// ----- global variables ----------

// ----- empty structs for initializations
static const UserVariables_t empty_UserVariables;
static const ConfigVars_t empty_GV;
static const LALUnit empty_LALUnit;

int XLALWriteREAL4TimeSeries2fp ( FILE *fp, const REAL4TimeSeries *TS );

// ---------- local prototypes ----------
int XLALInitUserVars ( UserVariables_t *uvar, int argc, char *argv[] );
int XLALInitMakefakedata ( ConfigVars_t *cfg, UserVariables_t *uvar );

int XLALWriteMFDlog ( const char *logfile, const ConfigVars_t *cfg );
int XLALFreeMem ( ConfigVars_t *cfg );

BOOLEAN is_directory ( const CHAR *fname );

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  int len;
  ConfigVars_t GV = empty_GV;
  UserVariables_t uvar = empty_UserVariables;

  /* ------------------------------
   * read user-input and set up shop
   *------------------------------*/
  XLAL_CHECK ( XLALInitUserVars ( &uvar, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALInitMakefakedata ( &GV, &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  MultiSFTVector *mSFTs = NULL;
  MultiREAL4TimeSeries *mTseries = NULL;

  CWMFDataParams DataParams   = empty_CWMFDataParams;
  DataParams.fMin               = uvar.fmin;
  DataParams.Band               = uvar.Band;
  DataParams.detInfo            = GV.detInfo;
  DataParams.multiTimestamps 	= (*GV.multiTimestamps);
  DataParams.randSeed           = uvar.randSeed;
  DataParams.SFTWindowType      = uvar.SFTWindowType;
  DataParams.SFTWindowBeta      = uvar.SFTWindowBeta;

  XLAL_CHECK ( XLALCWMakeFakeMultiData ( &mSFTs, &mTseries, GV.injectionSources, &DataParams, GV.edat ) == XLAL_SUCCESS, XLAL_EFUNC );

  // if noiseSFTs specified, load them and add them to the resulting SFT-vector
  if ( GV.multiNoiseCatalogView )
    {
      SFTtype *sft0 = &(mSFTs->data[0]->data[0]);
      /* load effective frequency-band from noise-SFTs */
      UINT4 numBins = sft0->data->length;
      REAL8 dFreq   = sft0->deltaF;
      REAL8 fMin    = sft0->f0;
      REAL8 fMax    = fMin + ( numBins - 1 ) * dFreq;
      MultiSFTVector *mNoiseSFTs;
      XLAL_CHECK ( (mNoiseSFTs = XLALLoadMultiSFTsFromView ( GV.multiNoiseCatalogView, fMin, fMax )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( XLALMultiSFTVectorAdd ( mSFTs, mNoiseSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroyMultiSFTVector ( mNoiseSFTs );
    }

  if (uvar.outSFTdir)
    {
      XLAL_CHECK ( is_directory ( uvar.outSFTdir ), XLAL_EINVAL );

      /* generate comment string */
      CHAR *logstr;
      XLAL_CHECK ( (logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
      char *comment = XLALCalloc ( 1, len = strlen ( logstr ) + strlen(GV.VCSInfoString) + 512 );
      XLAL_CHECK ( comment != NULL, XLAL_ENOMEM, "XLALCalloc(1,%d) failed.\n", len );
      sprintf ( comment, "Generated by:\n%s\n%s\n", logstr, GV.VCSInfoString );

      for ( UINT4 X=0; X < mSFTs->length; X ++ )
        {
          SFTVector *sfts = mSFTs->data[X];
          /* either write whole SFT-vector to single concatenated file */
          if ( uvar.outSingleSFT ) {
            XLAL_CHECK ( XLALWriteSFTVector2File( sfts, uvar.outSFTdir, comment, uvar.outMiscField ) == XLAL_SUCCESS, XLAL_EFUNC );
          } else {	// or as individual SFT-files
            XLAL_CHECK ( XLALWriteSFTVector2Dir( sfts, uvar.outSFTdir, comment, uvar.outMiscField ) == XLAL_SUCCESS, XLAL_EFUNC );
          }
        } // for X < numIFOs

      XLALFree ( logstr );
      XLALFree ( comment );

    } /* if outSFTdir */


   /* output ASCII time-series if requested */
  if ( uvar.TDDfile )
    {
      FILE *fp;
      CHAR *fname = XLALCalloc (1, len = strlen(uvar.TDDfile) + 10 );
      XLAL_CHECK ( fname != NULL, XLAL_ENOMEM, "XLALCalloc(1,%d) failed\n", len );

      for ( UINT4 X=0; X < mTseries->length; X ++ )
        {
          const REAL4TimeSeries *TS = mTseries->data[X];
          sprintf (fname, "%c%c-%s", TS->name[0], TS->name[1], uvar.TDDfile );

          XLAL_CHECK ( (fp = fopen (fname, "wb")) != NULL, XLAL_EIO, "Failed to fopen TDDfile = '%s' for writing\n", fname );

          XLAL_CHECK ( XLALWriteREAL4TimeSeries2fp ( fp, TS ) == XLAL_SUCCESS, XLAL_EFUNC );

	  fclose(fp);
        } // for X < numDet

      XLALFree (fname);
    } /* if outputting ASCII time-series */

  /* ---------- free memory ---------- */
  XLALDestroyMultiREAL4TimeSeries ( mTseries );
  XLALDestroyMultiSFTVector ( mSFTs );

  XLALFreeMem ( &GV );	/* free the config-struct */

  LALCheckMemoryLeaks();

  return 0;
} /* main */

/**
 *  Handle user-input and set up shop accordingly, and do all
 * consistency-checks on user-input.
 */
int
XLALInitMakefakedata ( ConfigVars_t *cfg, UserVariables_t *uvar )
{
  int len;
  XLAL_CHECK ( cfg != NULL, XLAL_EINVAL, "Invalid NULL input 'cfg'\n" );
  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL, "Invalid NULL input 'uvar'\n");

  cfg->VCSInfoString = XLALGetVersionString(0);
  XLAL_CHECK ( cfg->VCSInfoString != NULL, XLAL_EFUNC, "XLALGetVersionString(0) failed.\n" );

  // version info was requested: output then exit
  if ( uvar->version )
    {
      printf ("%s\n", cfg->VCSInfoString );
      exit (0);
    }

  /* if requested, log all user-input and code-versions */
  if ( uvar->logfile ) {
    XLAL_CHECK ( XLALWriteMFDlog ( uvar->logfile, cfg ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALWriteMFDlog() failed with xlalErrno = %d\n", xlalErrno );
  }

   /* check for negative fMin and Band, which would break the fMin_eff, fBand_eff calculation below */
  XLAL_CHECK ( uvar->fmin >= 0, XLAL_EDOM, "Invalid negative frequency fMin=%f!\n\n", uvar->fmin );
  XLAL_CHECK ( uvar->Band >= 0, XLAL_EDOM, "Invalid negative frequency band Band=%f!\n\n", uvar->Band );

  XLAL_CHECK ( XLALParseMultiDetectorInfo ( &(cfg->detInfo), uvar->IFOs, uvar->sqrtSX ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* ---------- determine timestamps to produce signal for  ---------- */
  {
    /* check input consistency: *uvar->timestampsFiles, uvar->startTime, uvar->duration */
    BOOLEAN haveStart = XLALUserVarWasSet ( &uvar->startTime );
    BOOLEAN haveDuration = XLALUserVarWasSet ( &uvar->duration );
    BOOLEAN haveOverlap = ( uvar->SFToverlap > 0 );
    BOOLEAN haveTimestampsFiles = ( uvar->timestampsFiles != NULL );

    if ( ( haveDuration && !haveStart) || ( !haveDuration && haveStart ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "Need BOTH --startTime AND --duration if you give one of them !\n\n");
    }

    /* don't allow using --SFToverlap with anything other than pure (--startTime,--duration) */
    if ( haveOverlap && ( uvar->noiseSFTs || haveTimestampsFiles ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "I can't combine --SFToverlap with --noiseSFTs or --timestampsFiles, only use with (--startTime, --duration)!\n\n");
    }

    /* ----- load timestamps from file if given  */
    if ( haveTimestampsFiles )
      {
	if ( haveStart || haveDuration || haveOverlap ) {
          XLAL_ERROR ( XLAL_EINVAL, "Using --timestampsFiles is incompatible with either of --startTime, --duration or --SFToverlap\n\n");
        }
	XLAL_CHECK ( (cfg->multiTimestamps = XLALReadMultiTimestampsFiles ( uvar->timestampsFiles )) != NULL, XLAL_EFUNC );

	if ( ( cfg->multiTimestamps->length == 0 ) || ( cfg->multiTimestamps->data == NULL ) ) {
          XLAL_ERROR ( XLAL_EINVAL, "Got empty timestamps-list from file '%s'\n", uvar->timestampsFiles );
        }
        for ( UINT4 X=0; X < cfg->multiTimestamps->length; X ++ ) {
          cfg->multiTimestamps->data[X]->deltaT = uvar->Tsft;
        }

        XLAL_CHECK ( uvar->noiseSFTs == NULL, XLAL_EINVAL, "--timestampsFiles is incompatible with --noiseSFTs\n" );

      } /* if haveTimestampsFile */

    /* ----- if real noise-SFTs given: load them now using optional (start,start+duration) as constraints,
     * Require window option to be given to ensure consistency with noise-SFT's windowing
     */
    if ( uvar->noiseSFTs )
      {
        XLAL_CHECK ( XLALUserVarWasSet ( &uvar->SFTWindowType ), XLAL_EINVAL, "--SFTWindowType required when given noiseSFTs.\n" );
        XLAL_CHECK ( !XLALUserVarWasSet ( &uvar->IFOs ), XLAL_EINVAL, "Specifying --IFOs not allowed with --noiseSFTs.\n" );

	/* use additional time-constraints from user input in noiseSFT loading */
        SFTConstraints constraints = empty_SFTConstraints;
	if ( haveStart && haveDuration )
	  {
            LIGOTimeGPS minStartTime, maxEndTime;
	    XLALGPSSetREAL8 ( &minStartTime, uvar->startTime );
	    XLALGPSSetREAL8 ( &maxEndTime, uvar->startTime + uvar->duration );

	    constraints.startTime = &minStartTime;
	    constraints.endTime   = &maxEndTime;

            XLALPrintWarning ( "Only noise-SFTs between GPS [%d, %d] will be used!\n", uvar->startTime, uvar->startTime + uvar->duration );
	  } /* if start+duration given */

        XLAL_CHECK ( (cfg->noiseCatalog = XLALSFTdataFind ( uvar->noiseSFTs, &constraints )) != NULL, XLAL_EFUNC );
	XLAL_CHECK (  cfg->noiseCatalog->length > 0, XLAL_EINVAL, "No noise-SFTs matching (start+duration, timestamps) were found!\n" );
	XLAL_CHECK ( (cfg->multiNoiseCatalogView = XLALMultiSFTCatalogView ( cfg->noiseCatalog )) != NULL, XLAL_EFUNC );

	/* extract multi-timestamps from the multi-SFT-catalog view */
        XLAL_CHECK ( (cfg->multiTimestamps = XLALTimestampsFromMultiSFTCatalogView ( cfg->multiNoiseCatalogView )) != NULL, XLAL_EFUNC );

      } /* if uvar->noiseSFTs */

    /* have we got our timestamps yet?: If not, we must get them from (start, duration) user-input */
    if ( cfg->multiTimestamps == NULL )
      {
	if ( !haveStart || !haveDuration ) {
          XLAL_ERROR ( XLAL_EINVAL, "Need to have either --timestampsFiles OR (--startTime,--duration) OR --noiseSFTs\n\n");
        }

	/* internally always use timestamps, so we generate them  */
	LIGOTimeGPS tStart;
	XLALGPSSetREAL8 ( &tStart, uvar->startTime );
        XLAL_CHECK ( ( cfg->multiTimestamps = XLALMakeMultiTimestamps ( tStart, uvar->duration, uvar->Tsft, uvar->SFToverlap, cfg->detInfo.length )) != NULL, XLAL_EFUNC );

      } /* if !cfg->multiTimestamps */

  } /* END: setup signal start + duration */


  /* -------------------- Prepare quantities for barycentering -------------------- */
  {
    CHAR *earthdata, *sundata;

    len = strlen(uvar->ephemYear) + 20;
    if ( uvar->ephemDir ) {
      len += strlen ( uvar->ephemDir );
    }
    XLAL_CHECK ( (earthdata = XLALCalloc(1, len)) != NULL, XLAL_ENOMEM );
    XLAL_CHECK ( (sundata   = XLALCalloc(1, len)) != NULL, XLAL_ENOMEM );
    const char *sep = uvar->ephemDir ? "/" : "";
    const char *ephemDir = uvar->ephemDir ? uvar->ephemDir : "";
    sprintf ( earthdata, "%s%searth%s.dat", ephemDir, sep, uvar->ephemYear);
    sprintf ( sundata,   "%s%ssun%s.dat",   ephemDir, sep, uvar->ephemYear);

    /* Init ephemerides */
    XLAL_CHECK ( ( cfg->edat = XLALInitBarycenter ( earthdata, sundata ) ) != NULL, XLAL_EFUNC );
    XLALFree(earthdata);
    XLALFree(sundata);

  } /* END: prepare barycentering routines */

  // --------------------------------------------------------------------------------
  // handle signal input parameters
  // --------------------------------------------------------------------------------
  XLAL_CHECK ( ( cfg->injectionSources = XLALPulsarParamsFromUserInput ( uvar->injectionSources ) ) != NULL, XLAL_EFUNC );

  return XLAL_SUCCESS;

} /* XLALInitMakefakedata() */


/**
 * Register all "user-variables", and initialized them from command-line and config-files
 */
int
XLALInitUserVars ( UserVariables_t *uvar, int argc, char *argv[] )
{
  int ret, len;

  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL, "Invalid NULL input 'uvar'\n");
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL, "Invalid NULL input 'argv'\n");

  // ---------- set a few defaults ----------
#define EPHEM_YEARS  "00-19-DE405"
  uvar->ephemYear = XLALCalloc ( 1, len = strlen(EPHEM_YEARS)+1 );
  XLAL_CHECK ( uvar->ephemYear != NULL, XLAL_ENOMEM, "XLALCalloc ( 1, %d ) failed.\n", len );
  strcpy ( uvar->ephemYear, EPHEM_YEARS );

  uvar->ephemDir = NULL;
  uvar->Tsft = 1800;
  uvar->fmin = 0;	/* no heterodyning by default */
  uvar->Band = 8192;	/* 1/2 LIGO sampling rate by default */

#define MISC_DEFAULT "mfdv5"
  XLAL_CHECK ( (uvar->outMiscField = XLALCalloc ( 1, strlen(MISC_DEFAULT)+1 ) ) != NULL, XLAL_ENOMEM );
  strcpy ( uvar->outMiscField, MISC_DEFAULT );

  // ---------- register all our user-variable ----------
  XLALregBOOLUserStruct (  help,                'h', UVAR_HELP    , "Print this help/usage message");

  /* output options */
  XLALregBOOLUserStruct (   outSingleSFT,       's', UVAR_OPTIONAL, "Write a single concatenated SFT file instead of individual files" );
  XLALregSTRINGUserStruct ( outSFTdir,          'n', UVAR_OPTIONAL, "Output SFTs:  Output directory for SFTs");
  XLALregSTRINGUserStruct(  outMiscField,	  0, UVAR_OPTIONAL, "'misc' field entry in the SFT-filename (part of 'description' field)");

  XLALregSTRINGUserStruct ( TDDfile,            't', UVAR_OPTIONAL, "Filename to output time-series into");

  XLALregSTRINGUserStruct ( logfile,            'l', UVAR_OPTIONAL, "Filename for log-output");

  /* detectors and respective noise-floors */
  XLALregLISTUserStruct ( IFOs,			'I', UVAR_OPTIONAL, "CSV list of detectors, eg. \"H1,H2,L1,G1, ...\" ");
  XLALregLISTUserStruct ( sqrtSX,	 	 0,  UVAR_OPTIONAL, "Add Gaussian Noise: CSV list of detectors' noise-floors sqrt{Sn}");

  XLALregSTRINGUserStruct ( ephemDir,           'E', UVAR_OPTIONAL, "Directory path for ephemeris files (use LAL_DATA_PATH if unspecified)");
  XLALregSTRINGUserStruct ( ephemYear,          'y', UVAR_OPTIONAL, "Year-range string of ephemeris files to be used");

  /* start + duration of timeseries */
  XLALregINTUserStruct (  startTime,            'G', UVAR_OPTIONAL, "Start-time of requested signal in detector-frame (GPS seconds)");
  XLALregINTUserStruct (  duration,              0,  UVAR_OPTIONAL, "Duration of requested signal in seconds");
  XLALregLISTUserStruct ( timestampsFiles,       0,  UVAR_OPTIONAL, "ALTERNATIVE: File to read timestamps from (file-format: lines with <seconds> <nanoseconds>)");

  /* sampling and heterodyning frequencies */
  XLALregREALUserStruct (  fmin,                 0, UVAR_OPTIONAL, "Lowest frequency in output SFT (= heterodyning frequency)");
  XLALregREALUserStruct (  Band,                 0, UVAR_OPTIONAL, "Bandwidth of output SFT in Hz (= 1/2 sampling frequency)");

  /* SFT properties */
  XLALregREALUserStruct (  Tsft,                 0, UVAR_OPTIONAL, "Time baseline of one SFT in seconds");
  XLALregREALUserStruct (  SFToverlap,           0, UVAR_OPTIONAL, "Overlap between successive SFTs in seconds (conflicts with --noiseSFTs or --timestampsFile)");
  XLALregSTRINGUserStruct( SFTWindowType,        0, UVAR_OPTIONAL, "Window function to be applied to the SFTs (required when using --noiseSFTs)");
  XLALregREALUserStruct (  SFTWindowBeta,        0, UVAR_OPTIONAL, "Window 'beta' parameter required for a few window-types (eg. 'tukey')");

  /* pulsar params */
  XLALregSTRINGUserStruct( injectionSources,     0, UVAR_REQUIRED, "Either a file-specification (\"@file-pattern\") or a config-string defining the sources to inject" );

  /* noise */
  XLALregSTRINGUserStruct ( noiseSFTs,          'D', UVAR_OPTIONAL, "Noise-SFTs to be added to signal (Used also to set IFOs and timestamps)");

  XLALregBOOLUserStruct (  version,             'V', UVAR_SPECIAL, "Output version information");

  /* ----- 'expert-user/developer' options ----- */
  XLALregINTUserStruct (   randSeed,             0, UVAR_DEVELOPER, "Specify random-number seed for reproducible noise (use /dev/urandom otherwise).");

  /* read cmdline & cfgfile  */
  ret = XLALUserVarReadAllInput ( argc, argv );
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "Failed to parse user-input\n");

  if ( uvar->help ) {	/* if help was requested, we're done */
    exit (0);
  }

  return XLAL_SUCCESS;

} /* XLALInitUserVars() */


/**
 * This routine frees up all the memory.
 */
int
XLALFreeMem ( ConfigVars_t *cfg )
{
  XLAL_CHECK ( cfg != NULL, XLAL_EINVAL );

  /* Free config-Variables and userInput stuff */
  XLALDestroyUserVars();

  /* free timestamps if any */
  XLALDestroyMultiTimestamps ( cfg->multiTimestamps );

  XLALDestroyPulsarParamsVector ( cfg->injectionSources );

  XLALFree ( cfg->VCSInfoString );

  // free noise-SFT catalog
  XLALDestroySFTCatalog ( cfg->noiseCatalog );
  XLALDestroyMultiSFTCatalogView ( cfg->multiNoiseCatalogView );

  /* Clean up earth/sun Ephemeris tables */
  XLALDestroyEphemerisData ( cfg->edat );

  return XLAL_SUCCESS;

} /* XLALFreeMem() */

/** Log the all relevant parameters of this run into a log-file.
 * The name of the log-file used is uvar_logfile
 * <em>NOTE:</em> Currently this function only logs the user-input and code-versions.
 */
int
XLALWriteMFDlog ( const char *logfile, const ConfigVars_t *cfg )
{
  XLAL_CHECK ( logfile != NULL, XLAL_EINVAL, "Invalid NULL input 'logfile'\n" );
  XLAL_CHECK ( cfg != NULL, XLAL_EINVAL, "Invalid NULL input 'cfg'\n");

  FILE *fplog = fopen ( logfile, "wb" );
  XLAL_CHECK ( fplog != NULL, XLAL_EIO, "Failed to fopen log-file '%s' for writing ('wb').\n", logfile );

  /* write out a log describing the complete user-input (in cfg-file format) */
  CHAR *logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CFGFILE );
  XLAL_CHECK ( logstr != NULL, XLAL_EFUNC, "XLALUserVarGetLog(UVAR_LOGFMT_CFGFILE) failed.\n" );

  fprintf ( fplog, "## LOG-FILE of Makefakedata run\n\n");
  fprintf ( fplog, "# User-input: [formatted as config-file]\n");
  fprintf ( fplog, "# ----------------------------------------------------------------------\n\n");
  fprintf ( fplog, "%s", logstr);
  XLALFree ( logstr );

  /* write out a log describing the complete user-input (in commandline format) */
  logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE );
  XLAL_CHECK ( logstr != NULL, XLAL_EFUNC, "XLALUserVarGetLog(UVAR_LOGFMT_CMDLINE) failed.\n" );

  fprintf ( fplog, "\n\n# User-input: [formatted as commandline]\n");
  fprintf ( fplog, "# ----------------------------------------------------------------------\n\n");
  fprintf ( fplog, "%s", logstr );
  XLALFree ( logstr );

  /* append an VCS-version string of the code used */
  fprintf ( fplog, "\n\n# VCS-versions of executable:\n");
  fprintf ( fplog, "# ----------------------------------------------------------------------\n");
  fprintf ( fplog, "\n%s\n", cfg->VCSInfoString );
  fclose ( fplog );

  return XLAL_SUCCESS;

} /* XLALWriteMFDLog() */

/* determine if the given filepath is an existing directory or not */
BOOLEAN
is_directory ( const CHAR *fname )
{
  struct stat stat_buf;

  if ( stat (fname, &stat_buf) )
    return 0;

  if ( ! S_ISDIR (stat_buf.st_mode) )
    return 0;
  else
    return 1;

} /* is_directory() */

/**
 * Write a REAL4TimeSeries in a 2-column ASCII format (lines of "GPS_i  data_i")
 * into an open file 'fp'
 */
int
XLALWriteREAL4TimeSeries2fp ( FILE *fp, const REAL4TimeSeries *TS )
{
  XLAL_CHECK ( fp != NULL, XLAL_EINVAL );
  XLAL_CHECK ( TS != NULL, XLAL_EINVAL );
  XLAL_CHECK ( (TS->data != NULL) && (TS->data->length >0) && (TS->data->data != NULL), XLAL_EINVAL );

  REAL8 t0 = XLALGPSGetREAL8( &(TS->epoch) );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFAULT, "XLALGPSGetREAL8(%d,%d) failed\n", TS->epoch.gpsSeconds, TS->epoch.gpsNanoSeconds );

  for( UINT4 i = 0; i < TS->data->length; i++)
  {
    REAL8 ti = t0 + i * TS->deltaT;
    XLAL_CHECK ( fprintf( fp, "%16.9f %e\n", ti, TS->data->data[i] ) > 0, XLAL_EIO );
  }

  return XLAL_SUCCESS;

} /* XLALWriteREAL4TimeSeries2fp() */
