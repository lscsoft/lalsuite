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

#include <lalapps.h>

/***************************************************/
#define SQ(x) ( (x) * (x) )

/*----------------------------------------------------------------------*/

/**
 * Lists of CW signal params and "line-feature" sinusoid params to generate
 * data for
 */
typedef struct tagInjectionSources
{
  UINT4 numPulsars;		/**< number of pulsar-signals to inject */
  PulsarParams *pulsarParams;	/**< array of pulsar-signal parameters to inject */

  /* ... */
  /* may be extended to carry further injection sources, eg 'lines' here */
} InjectionSources;


/** configuration-variables derived from user-variables */
typedef struct
{
  InjectionSources *injectionSources;		/**< list of injection source parameters */

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
  REAL8 refTime;		/**< Pulsar reference time tRef in SSB ('0' means: use startTime converted to SSB) */
  REAL8 refTimeMJD;             /**< Pulsar reference time tRef in MJD ('0' means: use startTime converted to SSB) */

  REAL8 h0;			/**< overall signal amplitude h0 */
  REAL8 cosi;		/**< cos(iota) of inclination angle iota */
  REAL8 aPlus;		/**< ALTERNATIVE to {h0,cosi}: Plus polarization amplitude aPlus */
  REAL8 aCross;		/**< ALTERNATIVE to {h0,cosi}: Cross polarization amplitude aCross */
  REAL8 psi;			/**< Polarization angle psi */
  REAL8 phi0;		/**< Initial phase phi */

  REAL8 Alpha;		/**< Right ascension [radians] alpha of pulsar */
  REAL8 Delta;		/**< Declination [radians] delta of pulsar */
  CHAR *RA;		/**< Right ascension [hh:mm:ss.ssss] alpha of pulsar */
  CHAR *Dec;	        /**< Declination [dd:mm:ss.ssss] delta of pulsar */

  REAL8 Freq;

  REAL8 f1dot;		/**< First spindown parameter f' */
  REAL8 f2dot;		/**< Second spindown parameter f'' */
  REAL8 f3dot;		/**< Third spindown parameter f''' */

  /* orbital parameters [OPTIONAL] */

  REAL8 orbitasini;	        /**< Projected orbital semi-major axis in seconds (a/c) */
  REAL8 orbitEcc;	        /**< Orbital eccentricity */
  INT4  orbitTpSSBsec;	        /**< 'observed' (SSB) time of periapsis passage. Seconds. */
  INT4  orbitTpSSBnan;	        /**< 'observed' (SSB) time of periapsis passage. Nanoseconds. */
  REAL8 orbitTpSSBMJD;          /**< 'observed' (SSB) time of periapsis passage. MJD. */
  REAL8 orbitPeriod;		/**< Orbital period (seconds) */
  REAL8 orbitArgp;	        /**< Argument of periapsis (radians) */

  BOOLEAN lineFeature;	        /**< generate a monochromatic line instead of a pulsar-signal */

  BOOLEAN version;		/**< output version information */

  INT4 randSeed;		/**< allow user to specify random-number seed for reproducible noise-realizations */

  CHAR *transientWindowType;	/**< name of transient window ('rect', 'exp',...) */
  REAL8 transientStartTime;	/**< GPS start-time of transient window */
  REAL8 transientTauDays;	/**< time-scale in days of transient window */
} UserVariables_t;


// ---------- exportable API types ----------
/**
 * Struct controlling all the aspects of the fake data (time-series + SFTs)
 * to be produced by XLALMakeFakeCWData()
 */
typedef struct tagCWMultiDataParams
{
  REAL8 fMin;					//!< smallest frequency guaranteed to be generated ["effective" fMin can be smaller]
  REAL8 Band;					//!< smallest frequency band guaranteed to be generated ["effective" band can be larger]
  MultiDetectorInfo detInfo;			//!< detectors and noise-floors (for Gaussian noise) to generate data for
  MultiLIGOTimeGPSVector *multiTimestamps;	//!< timestamps to generate SFTs for
  const char *SFTWindowType;			//!< window to apply to the SFT timeseries
  REAL8 SFTWindowBeta;				//!< 'beta' parameter required for *some* windows [otherwise must be 0]
  UINT4 randSeed;				//!< seed value for random-number generator
} CWMultiDataParams;

typedef struct tagDataParams
{
  REAL8 fMin;					//!< smallest frequency guaranteed to be generated ["effective" fMin can be smaller]
  REAL8 Band;					//!< smallest frequency band guaranteed to be generated ["effective" band can be larger]
  LALDetector site;
  REAL8 sqrtSn;
  LIGOTimeGPSVector *timestamps;		//!< timestamps to generate SFTs for
  const char *SFTWindowType;			//!< window type ('rectangular', 'hann', ...) to apply to the SFT timeseries
  REAL8 SFTWindowBeta;				//!< 'beta' parameter required for *some* windows [otherwise must be 0]
  UINT4 randSeed;				//!< seed value for random-number generator
} CWDataParams;

// ----- global variables ----------

// ----- empty structs for initializations
static const UserVariables_t empty_UserVariables;
static const ConfigVars_t empty_GV;
static const LALUnit empty_LALUnit;
static const CWDataParams empty_CWDataParams;
static const CWMultiDataParams empty_CWMultiDataParams;

// ---------- exportable API prototypes ----------
int XLALFindSmallestValidSamplingRate ( UINT4 *n1, UINT4 n0, const LIGOTimeGPSVector *timestamps );
int
XLALCWMakeFakeMultiData ( MultiSFTVector **multiSFTs, MultiREAL4TimeSeries **multiTseries,
                          const InjectionSources *injectionSources, const CWMultiDataParams *multiDataParams, EphemerisData *edat );
int
XLALCWMakeFakeData ( SFTVector **SFTVect, REAL4TimeSeries **Tseries,
                     const InjectionSources *injectionSources, const CWDataParams *dataParams, EphemerisData *edat );
REAL4TimeSeries *
XLALGenerateCWSignalTS ( const PulsarParams *pulsarParams, const LALDetector *site, LIGOTimeGPS startTime, REAL8 duration, REAL8 fSamp, REAL8 fHet, EphemerisData *edat );


void XLALDestroyInjectionSources ( InjectionSources *sources );
int XLALWriteREAL4TimeSeries2fp ( FILE *fp, const REAL4TimeSeries *TS );

// ---------- local prototypes ----------
int XLALInitUserVars ( UserVariables_t *uvar, int argc, char *argv[] );
int XLALInitMakefakedata ( ConfigVars_t *cfg, UserVariables_t *uvar );

int XLALWriteMFDlog ( const char *logfile, const ConfigVars_t *cfg );
int XLALFreeMem ( ConfigVars_t *cfg );

BOOLEAN is_directory ( const CHAR *fname );
static UINT4 gcd (UINT4 numer, UINT4 denom);

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

  CWMultiDataParams DataParams   = empty_CWMultiDataParams;
  DataParams.fMin               = uvar.fmin;
  DataParams.Band               = uvar.Band;
  DataParams.detInfo            = GV.detInfo;
  DataParams.multiTimestamps 	= GV.multiTimestamps;
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

      const char *misc = "mfdv5";
      for ( UINT4 X=0; X < mSFTs->length; X ++ )
        {
          SFTVector *sfts = mSFTs->data[X];
          /* either write whole SFT-vector to single concatenated file */
          if ( uvar.outSingleSFT ) {
            XLAL_CHECK ( XLALWriteSFTVector2File( sfts, uvar.outSFTdir, comment, misc ) == XLAL_SUCCESS, XLAL_EFUNC );
          } else {	// or as individual SFT-files
            XLAL_CHECK ( XLALWriteSFTVector2Dir( sfts, uvar.outSFTdir, comment, misc ) == XLAL_SUCCESS, XLAL_EFUNC );
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

  // FIXME: hardcoded to 1 CW signal injection right now
  UINT4 numPulsars = 1;
  XLAL_CHECK ( ( cfg->injectionSources = XLALCalloc ( 1, sizeof(*cfg->injectionSources) )) != NULL, XLAL_ENOMEM );
  cfg->injectionSources->numPulsars = numPulsars;
  XLAL_CHECK ( ( cfg->injectionSources->pulsarParams = XLALCalloc ( numPulsars, sizeof(*cfg->injectionSources->pulsarParams) )) != NULL, XLAL_ENOMEM );
  PulsarParams *pulsarParams = &( cfg->injectionSources->pulsarParams[0] );

  { /* ========== translate user-input into 'PulsarParams' struct ========== */
    BOOLEAN have_h0     = XLALUserVarWasSet ( &uvar->h0 );
    BOOLEAN have_cosi   = XLALUserVarWasSet ( &uvar->cosi );
    BOOLEAN have_aPlus  = XLALUserVarWasSet ( &uvar->aPlus );
    BOOLEAN have_aCross = XLALUserVarWasSet ( &uvar->aCross );
    BOOLEAN have_Freq   = XLALUserVarWasSet ( &uvar->Freq );
    BOOLEAN have_Alpha  = XLALUserVarWasSet ( &uvar->Alpha );
    BOOLEAN have_Delta  = XLALUserVarWasSet ( &uvar->Delta );
    BOOLEAN have_RA = XLALUserVarWasSet ( &uvar->RA );
    BOOLEAN have_Dec = XLALUserVarWasSet ( &uvar->Dec );

    /* ----- {h0,cosi} or {aPlus,aCross} ----- */
    if ( (have_aPlus || have_aCross) && ( have_h0 || have_cosi ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "Need to specify EITHER {h0,cosi} OR {aPlus, aCross}!\n\n");
    }
    if ( (have_h0 ^ have_cosi) ) {
      XLAL_ERROR ( XLAL_EINVAL, "Need BOTH --h0 and --cosi!\n\n");
    }
    if ( (have_aPlus ^ have_aCross) ) {
      XLAL_ERROR ( XLAL_EINVAL, "Need BOTH --aPlus and --aCross !\n\n");
    }

    if ( have_h0 && have_cosi )
      {
	pulsarParams->Amp.h0 = uvar->h0;
	pulsarParams->Amp.cosi = uvar->cosi;
      }
    else if ( have_aPlus && have_aCross )
      {  /* translate A_{+,x} into {h_0, cosi} */
	REAL8 disc;
	if ( fabs(uvar->aCross) > uvar->aPlus ) {
          XLAL_ERROR ( XLAL_EINVAL, "Invalid input parameters: |aCross| = %g must be <= than aPlus = %g.\n", fabs(uvar->aCross), uvar->aPlus );
        }
	disc = sqrt ( SQ(uvar->aPlus) - SQ(uvar->aCross) );
	pulsarParams->Amp.h0   = uvar->aPlus + disc;
        if ( pulsarParams->Amp.h0 > 0 )
          pulsarParams->Amp.cosi = uvar->aCross / pulsarParams->Amp.h0;	// avoid division by 0!
        else
          pulsarParams->Amp.cosi = 0;
      }
    else {
      pulsarParams->Amp.h0 = 0.0;
      pulsarParams->Amp.cosi = 0.0;
    }
    pulsarParams->Amp.phi0 = uvar->phi0;
    pulsarParams->Amp.psi  = uvar->psi;

    /* ----- signal Frequency ----- */
    if ( have_Freq )
      pulsarParams->Doppler.fkdot[0] = uvar->Freq;
    else
      pulsarParams->Doppler.fkdot[0] = 0.0;

    /* ----- skypos ----- */
    if ( (have_Alpha || have_Delta) && (have_RA || have_Dec) ) {
      XLAL_ERROR ( XLAL_EINVAL, "Use EITHER {Alpha, Delta} OR {RA, Dec}\n\n");
    }
    if ( (have_Alpha && !have_Delta) || ( !have_Alpha && have_Delta ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "\nSpecify skyposition: need BOTH --Alpha and --Delta!\n\n");
    }
    if ( (have_RA && !have_Dec) || ( !have_RA && have_Dec ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "\nSpecify skyposition: need BOTH --RA and --Dec!\n\n");
    }
    if ( have_Alpha )
      {
	pulsarParams->Doppler.Alpha = uvar->Alpha;
	pulsarParams->Doppler.Delta = uvar->Delta;
      }
    else if ( have_RA )
      {
	pulsarParams->Doppler.Alpha = LALDegsToRads(uvar->RA,"alpha");
	pulsarParams->Doppler.Delta = LALDegsToRads(uvar->Dec,"delta");
      }

  } /* Pulsar signal parameters */

  /* ---------- prepare vector of spindown parameters ---------- */
  {
    pulsarParams->Doppler.fkdot[1] = uvar->f1dot;
    pulsarParams->Doppler.fkdot[2] = uvar->f2dot;
    pulsarParams->Doppler.fkdot[3] = uvar->f3dot;

  } /* END: prepare spindown parameters */


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

    if (XLALUserVarWasSet(&uvar->ephemDir) )
      len += strlen (uvar->ephemDir);

    if ( (earthdata = XLALCalloc(1, len)) == NULL) {
      XLAL_ERROR ( XLAL_ENOMEM, "earthdata = XLALCalloc(1, %d) failed.\n", len );
    }
    if ( (sundata = XLALCalloc(1, len)) == NULL) {
      XLAL_ERROR ( XLAL_ENOMEM, "sundata = XLALCalloc(1, %d) failed.\n", len );
    }

    if (XLALUserVarWasSet(&uvar->ephemDir) )
      {
	sprintf ( earthdata, "%s/earth%s.dat", uvar->ephemDir, uvar->ephemYear);
	sprintf ( sundata, "%s/sun%s.dat", uvar->ephemDir, uvar->ephemYear);
      }
    else
      {
	sprintf ( earthdata, "earth%s.dat", uvar->ephemYear);
	sprintf ( sundata, "sun%s.dat",  uvar->ephemYear);
      }

    /* Init ephemerides */
    cfg->edat = XLALInitBarycenter ( earthdata, sundata );
    XLAL_CHECK ( cfg->edat != NULL, XLAL_EFUNC, "XLALInitBarycenter() failed.\n" );
    XLALFree(earthdata);
    XLALFree(sundata);

  } /* END: prepare barycentering routines */

  /* -------------------- handle binary orbital params if given -------------------- */

  /* Consistency check: if any orbital parameters specified, we need all of them (except for nano-seconds)! */
  {
    BOOLEAN set1 = XLALUserVarWasSet(&uvar->orbitasini);
    BOOLEAN set2 = XLALUserVarWasSet(&uvar->orbitEcc);
    BOOLEAN set3 = XLALUserVarWasSet(&uvar->orbitPeriod);
    BOOLEAN set4 = XLALUserVarWasSet(&uvar->orbitArgp);
    BOOLEAN set5 = XLALUserVarWasSet(&uvar->orbitTpSSBsec);
    BOOLEAN set6 = XLALUserVarWasSet(&uvar->orbitTpSSBnan);
    BOOLEAN set7 = XLALUserVarWasSet(&uvar->orbitTpSSBMJD);
    BinaryOrbitParams *orbit = NULL;

    if (set1 || set2 || set3 || set4 || set5 || set6 || set7)
    {
      if ( (uvar->orbitasini > 0) && !(set1 && set2 && set3 && set4 && (set5 || set7)) ) {
        XLAL_ERROR ( XLAL_EINVAL, "\nPlease either specify  ALL orbital parameters or NONE!\n\n");
      }
      if ( (uvar->orbitEcc < 0) || (uvar->orbitEcc > 1) ) {
        XLAL_ERROR ( XLAL_EINVAL, "\nEccentricity = %g has to lie within [0, 1]\n\n", uvar->orbitEcc );
      }
      orbit = XLALCalloc ( 1, len = sizeof(BinaryOrbitParams) );
      XLAL_CHECK ( orbit != NULL, XLAL_ENOMEM, "XLALCalloc (1, %d) failed.\n", len );

      if ( set7 && (!set5 && !set6) )
	{
	  /* convert MJD peripase to GPS using Matt Pitkins code found at lal/packages/pulsar/src/BinaryPulsarTimeing.c */
	  REAL8 GPSfloat;
	  GPSfloat = XLALTTMJDtoGPS(uvar->orbitTpSSBMJD);
	  XLALGPSSetREAL8(&(orbit->tp),GPSfloat);
	}
      else if ((set5 && set6) && !set7)
	{
	  orbit->tp.gpsSeconds = uvar->orbitTpSSBsec;
	  orbit->tp.gpsNanoSeconds = uvar->orbitTpSSBnan;
	}
      else if ((set7 && set5) || (set7 && set6))
	{
	  XLAL_ERROR ( XLAL_EINVAL, "\nPlease either specify time of periapse in GPS OR MJD, not both!\n\n");
	}

      /* fill in orbital parameter structure */
      orbit->period = uvar->orbitPeriod;
      orbit->asini = uvar->orbitasini;
      orbit->argp = uvar->orbitArgp;
      orbit->ecc = uvar->orbitEcc;

      pulsarParams->Doppler.orbit = orbit;     /* struct copy */
    } /* if one or more orbital parameters were set */
    else
      pulsarParams->Doppler.orbit = NULL;
  } /* END: binary orbital params */

  /* ----- set "pulsar reference time", i.e. SSB-time at which pulsar params are defined ---------- */
  if (XLALUserVarWasSet(&uvar->refTime) && XLALUserVarWasSet(&uvar->refTimeMJD))
    {
      XLAL_ERROR ( XLAL_EINVAL, "\nUse only one of '--refTime' and '--refTimeMJD' to specify SSB reference time!\n\n");
    }
  else if (XLALUserVarWasSet(&uvar->refTime))
    {
      XLALGPSSetREAL8(&(pulsarParams->Doppler.refTime), uvar->refTime);
    }
  else if (XLALUserVarWasSet(&uvar->refTimeMJD))
    {

      /* convert MJD to GPS using Matt Pitkins code found at lal/packages/pulsar/src/BinaryPulsarTimeing.c */
      REAL8 GPSfloat;
      GPSfloat = XLALTTMJDtoGPS(uvar->refTimeMJD);
      XLALGPSSetREAL8(&(cfg->pulsar.Doppler.refTime),GPSfloat);
    }
  else {
    XLAL_ERROR ( XLAL_EINVAL, "Must specify explicit reference time (--refTime or --refTimeMJD)!\n");
  }

  /* ----- handle transient-signal window if given ----- */
  if ( !XLALUserVarWasSet ( &uvar->transientWindowType ) || !strcmp ( uvar->transientWindowType, "none") )
    cfg->transientWindow.type = TRANSIENT_NONE;                /* default: no transient signal window */
  else
    {
      if ( !strcmp ( uvar->transientWindowType, "rect" ) )
       {
         cfg->transientWindow.type = TRANSIENT_RECTANGULAR;              /* rectangular window [t0, t0+tau] */
       }
      else if ( !strcmp ( uvar->transientWindowType, "exp" ) )
        {
          cfg->transientWindow.type = TRANSIENT_EXPONENTIAL;            /* exponential decay window e^[-(t-t0)/tau for t>t0, 0 otherwise */
        }
      else
        {
          XLAL_ERROR ( XLAL_EINVAL, "Illegal transient window '%s' specified: valid are 'none', 'rect' or 'exp'\n", uvar->transientWindowType );
        }

      cfg->transientWindow.t0   = uvar->transientStartTime;
      cfg->transientWindow.tau  = uvar->transientTauDays * LAL_DAYSID_SI;

    } /* if transient window != none */

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

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar->ephemDir = XLALCalloc ( 1, len = strlen(DEFAULT_EPHEMDIR)+1 );
  XLAL_CHECK ( uvar->ephemDir != NULL, XLAL_ENOMEM, "XLALCalloc ( 1, %d ) failed.\n", len );
  strcpy (uvar->ephemDir, DEFAULT_EPHEMDIR );

  uvar->Tsft = 1800;
  uvar->fmin = 0;	/* no heterodyning by default */
  uvar->Band = 8192;	/* 1/2 LIGO sampling rate by default */

#define DEFAULT_TRANSIENT "none"
  uvar->transientWindowType = XLALCalloc ( 1, len = strlen(DEFAULT_TRANSIENT)+1 );
  XLAL_CHECK ( uvar->transientWindowType != NULL, XLAL_ENOMEM, "XLALCalloc ( 1, %d ) failed.\n", len );
  strcpy ( uvar->transientWindowType, DEFAULT_TRANSIENT );

  // ---------- register all our user-variable ----------
  XLALregBOOLUserStruct (  help,                'h', UVAR_HELP    , "Print this help/usage message");

  /* output options */
  XLALregBOOLUserStruct (   outSingleSFT,       's', UVAR_OPTIONAL, "Write a single concatenated SFT file instead of individual files" );
  XLALregSTRINGUserStruct ( outSFTdir,          'n', UVAR_OPTIONAL, "Output SFTs:  Output directory for SFTs");

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
  XLALregREALUserStruct (  refTime,             'S', UVAR_OPTIONAL, "Pulsar SSB reference time in GPS seconds (if 0: use startTime)");
  XLALregREALUserStruct (  refTimeMJD,           0 , UVAR_OPTIONAL, "ALTERNATIVE: Pulsar SSB reference time in MJD (if 0: use startTime)");

  XLALregREALUserStruct (  Alpha,                0, UVAR_OPTIONAL, "Right-ascension/longitude of pulsar in radians");
  XLALregSTRINGUserStruct (RA,                   0, UVAR_OPTIONAL, "ALTERNATIVE: Righ-ascension/longitude of pulsar in HMS 'hh:mm:ss.ssss'");

  XLALregREALUserStruct (  Delta,                0, UVAR_OPTIONAL, "Declination/latitude of pulsar in radians");
  XLALregSTRINGUserStruct (Dec,                  0, UVAR_OPTIONAL, "ALTERNATIVE: Declination/latitude of pulsar in DMS 'dd:mm:ss.ssss'");

  XLALregREALUserStruct (  h0,                   0, UVAR_OPTIONAL, "Overall signal-amplitude h0");
  XLALregREALUserStruct (  cosi,                 0, UVAR_OPTIONAL, "cos(iota) of inclination-angle iota");
  XLALregREALUserStruct (  aPlus,                0, UVAR_OPTIONAL, "ALTERNATIVE to {--h0,--cosi}: A_+ amplitude");
  XLALregREALUserStruct (  aCross,               0, UVAR_OPTIONAL, "ALTERNATIVE to {--h0,--cosi}: A_x amplitude");

  XLALregREALUserStruct (  psi,                  0, UVAR_OPTIONAL, "Polarization angle psi");
  XLALregREALUserStruct (  phi0,                 0, UVAR_OPTIONAL, "Initial GW phase phi");
  XLALregREALUserStruct (  Freq,                 0, UVAR_OPTIONAL, "Intrinsic GW-frequency at refTime");

  XLALregREALUserStruct (  f1dot,                0, UVAR_OPTIONAL, "First spindown parameter f' at refTime");
  XLALregREALUserStruct (  f2dot,                0, UVAR_OPTIONAL, "Second spindown parameter f'' at refTime");
  XLALregREALUserStruct (  f3dot,                0, UVAR_OPTIONAL, "Third spindown parameter f''' at refTime");

  /* binary-system orbital parameters */
  XLALregREALUserStruct (  orbitasini,           0, UVAR_OPTIONAL, "Projected orbital semi-major axis in seconds (a/c)");
  XLALregREALUserStruct (  orbitEcc,             0, UVAR_OPTIONAL, "Orbital eccentricity");
  XLALregINTUserStruct (   orbitTpSSBsec,        0, UVAR_OPTIONAL, "'true' (SSB) time of periapsis passage. Seconds.");
  XLALregINTUserStruct (   orbitTpSSBnan,        0, UVAR_OPTIONAL, "'true' (SSB) time of periapsis passage. Nanoseconds.");
  XLALregREALUserStruct (  orbitTpSSBMJD,        0, UVAR_OPTIONAL, "'true' (SSB) time of periapsis passage. MJD.");
  XLALregREALUserStruct (  orbitPeriod,          0, UVAR_OPTIONAL, "Orbital period (seconds)");
  XLALregREALUserStruct (  orbitArgp,            0, UVAR_OPTIONAL, "Argument of periapsis (radians)");

  /* noise */
  XLALregSTRINGUserStruct ( noiseSFTs,          'D', UVAR_OPTIONAL, "Noise-SFTs to be added to signal (Used also to set IFOs and timestamps)");

  XLALregBOOLUserStruct (  lineFeature,          0, UVAR_OPTIONAL, "Generate a monochromatic 'line' of amplitude h0 and frequency 'Freq'}");

  XLALregBOOLUserStruct (  version,             'V', UVAR_SPECIAL, "Output version information");

  /* transient signal window properties (name, start, duration) */
  XLALregSTRINGUserStruct (transientWindowType,  0, UVAR_OPTIONAL, "Type of transient signal window to use. ('none', 'rect', 'exp').");
  XLALregREALUserStruct (  transientStartTime,   0, UVAR_OPTIONAL, "GPS start-time 't0' of transient signal window.");
  XLALregREALUserStruct (  transientTauDays,     0, UVAR_OPTIONAL, "Timescale 'tau' of transient signal window in days.");

  /* ----- 'expert-user/developer' options ----- */
  XLALregINTUserStruct (   randSeed,             0, UVAR_DEVELOPER, "Specify random-number seed for reproducible noise (use /dev/urandom otherwise).");

  /* read cmdline & cfgfile  */
  ret = XLALUserVarReadAllInput ( argc, argv );
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "Failed to parse user-input\n");

  if ( uvar->help ) 	/* if help was requested, we're done */
    exit (0);

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

  XLALDestroyInjectionSources ( cfg->injectionSources );

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


// ----------------------------------------------------------------------------------------------------
// ---------- from here: will be exported into lalpulsar module --------------------
// ----------------------------------------------------------------------------------------------------

/**
 * Generate fake 'CW' data, returned either as SFTVector or REAL4Timeseries or both,
 * for given CW-signal ("pulsar") parameters and output parameters (frequency band etc)
 */
int
XLALCWMakeFakeMultiData ( MultiSFTVector **multiSFTs,		//< [out] pointer to optional SFT-vector for output [FIXME! Multi-]
                          MultiREAL4TimeSeries **multiTseries,	//< [out] pointer to optional timeseries-vector for output [FIXME! Multi-]
                          const InjectionSources *injectionSources,	//< [in] array of sources inject
                          const CWMultiDataParams *dataParams,		//< [in] parameters specifying the type of data to generate
                          EphemerisData *edat			//< [in] ephemeris data
                          )
{
  XLAL_CHECK ( (multiSFTs == NULL) || ((*multiSFTs) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (multiTseries == NULL) || ((*multiTseries) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (multiSFTs != NULL) || (multiTseries != NULL), XLAL_EINVAL );

  XLAL_CHECK ( injectionSources != NULL, XLAL_EINVAL );
  XLAL_CHECK ( dataParams != NULL, XLAL_EINVAL );
  XLAL_CHECK ( edat != NULL, XLAL_EINVAL );

  XLAL_CHECK ( dataParams->multiTimestamps != NULL, XLAL_EINVAL );
  MultiLIGOTimeGPSVector *multiTimestamps = dataParams->multiTimestamps;

  // check multi-detector input
  XLAL_CHECK ( dataParams->detInfo.length >= 1, XLAL_EINVAL );
  UINT4 numDet = dataParams->detInfo.length;
  XLAL_CHECK ( multiTimestamps->length == numDet, XLAL_EINVAL, "Inconsistent number of IFOs: detInfo says '%d', multiTimestamps says '%d'\n", numDet, multiTimestamps->length );

  // check Tsft, consistent over detectors
  REAL8 Tsft = multiTimestamps->data[0]->deltaT;
  XLAL_CHECK ( Tsft > 0, XLAL_EINVAL, "Got invalid Tsft = %g must be > 0\n", Tsft );
  for ( UINT4 X=0; X < numDet; X ++ ) {
    XLAL_CHECK ( multiTimestamps->data[X]->deltaT == Tsft, XLAL_EINVAL, "Inconsistent Tsft, for Tsft[X=0]=%g, while Tsft[X=%d]=%g\n", Tsft, X, multiTimestamps->data[X]->deltaT );
  }

  // ----- prepare output containers, as required
  MultiSFTVector *outMSFTs = NULL;
  MultiREAL4TimeSeries *outMTS = NULL;
  if ( multiSFTs != NULL )
    {
      XLAL_CHECK ( (outMSFTs = XLALCalloc ( 1, sizeof(*outMSFTs) )) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (outMSFTs->data = XLALCalloc ( numDet, sizeof(outMSFTs->data[0]))) != NULL, XLAL_ENOMEM );
      outMSFTs->length = numDet;
    } // if multiSFTs
  if ( multiTseries != NULL )
    {
      XLAL_CHECK ( (outMTS = XLALCalloc ( 1, sizeof(*outMTS) )) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (outMTS->data = XLALCalloc ( numDet, sizeof(outMTS->data[0]))) != NULL, XLAL_ENOMEM );
      outMTS->length = numDet;
    } // if multiTseries

  for ( UINT4 X=0; X < numDet; X ++ )
    {
      /* detector params */
      CWDataParams dataParamsX = empty_CWDataParams;
      dataParamsX.fMin = dataParams->fMin;
      dataParamsX.Band = dataParams->Band;
      dataParamsX.site = dataParams->detInfo.sites[X];
      dataParamsX.sqrtSn=dataParams->detInfo.sqrtSn[X];
      dataParamsX.timestamps = dataParams->multiTimestamps->data[X];
      dataParamsX.SFTWindowType = dataParams->SFTWindowType;
      dataParamsX.SFTWindowBeta = dataParams->SFTWindowBeta;
      dataParamsX.randSeed = dataParams->randSeed + X;	// increase seed in deterministic way: allows comparison w mfd_v4 !!

      SFTVector **svp = NULL;
      REAL4TimeSeries **tsp = NULL;
      if ( outMSFTs != NULL ) {
        svp = &(outMSFTs->data[X]);
      }
      if ( outMTS != NULL ) {
        tsp = &(outMTS->data[X]);
      }
      XLAL_CHECK ( XLALCWMakeFakeData ( svp, tsp, injectionSources, &dataParamsX, edat ) == XLAL_SUCCESS, XLAL_EFUNC );

    } // for X < numDet

  // return multi-Timeseries if requested
  if ( multiTseries ) {
    (*multiTseries) = outMTS;
  }
  // return multi-SFTs if requested
  if ( multiSFTs ) {
    (*multiSFTs) = outMSFTs;
  }

  return XLAL_SUCCESS;

} // XLALCWMakeFakeMultiData()

int
XLALCWMakeFakeData ( SFTVector **SFTvect,
                     REAL4TimeSeries **Tseries,
                     const InjectionSources *injectionSources,
                     const CWDataParams *dataParams,
                     EphemerisData *edat
                     )
{
  XLAL_CHECK ( (SFTvect == NULL) || ((*SFTvect) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (Tseries == NULL) || ((*Tseries) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (SFTvect != NULL) || (Tseries != NULL), XLAL_EINVAL );

  XLAL_CHECK ( injectionSources != NULL, XLAL_EINVAL );
  XLAL_CHECK ( dataParams != NULL, XLAL_EINVAL );
  XLAL_CHECK ( edat != NULL, XLAL_EINVAL );
  XLAL_CHECK ( dataParams->timestamps != NULL, XLAL_EINVAL );

  // initial default values fMin, sampling rate from caller input
  REAL8 fMin   = dataParams->fMin;
  REAL8 fBand = dataParams->Band;
  REAL8 fSamp = 2.0 * fBand;

  REAL8 Tsft = dataParams->timestamps->deltaT;

  // if SFT output requested: need *effective* fMin and Band consistent with SFT bins
  if ( SFTvect )
    {
      UINT4 firstBinEff, numBinsEff;
      XLAL_CHECK ( XLALFindCoveringSFTBins ( &firstBinEff, &numBinsEff, dataParams->fMin, dataParams->Band, Tsft ) == XLAL_SUCCESS, XLAL_EFUNC );

      REAL8 fBand_eff = (numBinsEff - 1.0) / Tsft;
      REAL8 fMin_eff  = firstBinEff / Tsft;
      REAL8 fMax = fMin + dataParams->Band;
      REAL8 fMax_eff = fMin_eff + fBand_eff;
      if ( (fMin_eff != fMin) || (fBand_eff != dataParams->Band ) ) {
        XLALPrintWarning("Caller asked for Band [%.16g, %.16g] Hz, effective SFT-Band produced is [%.16g, %.16g] Hz\n",
                         fMin, fMax, fMin_eff, fMax_eff );
      }
      fMin = fMin_eff;		// (potentially) lower minimal frequency to fit SFT bins
      fBand = fBand_eff;
      fSamp = 2.0 * fBand_eff;	// (potentially) higher sampling rate required to fit SFT bins
    } // if SFT-output

  /* characterize the output time-series */
  UINT4 n0_fSamp = (UINT4) round ( Tsft * fSamp );

  // by construction, fSamp * Tsft = integer, but if there are gaps between SFTs,
  // then we might have to sample at higher rates in order for all SFT-boundaries to
  // fall on exact timesteps of the timeseries.
  // ie we start from fsamp0 = n0_fSamp/Tsft, and then try to find the smallest
  // n1_fSamp >= n0_fSamp, such that for fsamp1 = n1_fSamp/Tsft, for all gaps i: Dt_i * fsamp1 = int
  UINT4 n1_fSamp;
  XLAL_CHECK ( XLALFindSmallestValidSamplingRate ( &n1_fSamp, n0_fSamp, dataParams->timestamps ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( n1_fSamp != n0_fSamp )
    {
      REAL8 fSamp1 = n1_fSamp / Tsft;	// increased sampling rate to fit all gaps
      XLALPrintWarning ( "GAPS: Initial SFT sampling frequency fSamp0= %d/%.0f = %g had to be increased to fSamp1 = %d/%.0f = %g\n",
                         n0_fSamp, Tsft, fSamp, n1_fSamp, Tsft, fSamp1 );
      fSamp = fSamp1;
    } // if higher effective sampling rate required

   /* ----- start-time and duration ----- */
  LIGOTimeGPS firstGPS = dataParams->timestamps->data[0];
  LIGOTimeGPS lastGPS  = dataParams->timestamps->data [ dataParams->timestamps->length - 1 ];
  REAL8 duration = XLALGPSDiff ( &lastGPS, &firstGPS ) + Tsft;
  XLAL_CHECK ( duration >= Tsft, XLAL_EINVAL, "Requested duration=%.0f sec is less than Tsft =%.0f sec.\n\n", duration, Tsft);

  REAL4TimeSeries *Tseries_sum = NULL;
  for ( UINT4 iInj = 0; iInj < injectionSources->numPulsars; iInj ++ )
    {
      const PulsarParams *sourceParams = &( injectionSources->pulsarParams[iInj] );

      REAL4TimeSeries *Tseries_i = NULL;
      XLAL_CHECK ( (Tseries_i = XLALGenerateCWSignalTS ( sourceParams, &(dataParams->site), firstGPS, duration, fSamp, fMin, edat )) != NULL, XLAL_EFUNC );

      if ( Tseries_sum == NULL )
        {
          Tseries_sum = Tseries_i;
        }
      else
        {
          XLAL_CHECK ( (Tseries_sum = XLALAddREAL4TimeSeries ( Tseries_sum, Tseries_i )) != NULL, XLAL_EFUNC );
          XLALDestroyREAL4TimeSeries ( Tseries_i );
        }

    } // for iInj < numSources

  /* add Gaussian noise if requested */
  REAL8 sqrtSn = dataParams->sqrtSn;
  if ( sqrtSn > 0)
    {
      REAL8 noiseSigma = sqrtSn * sqrt ( fBand );
      XLAL_CHECK ( XLALAddGaussianNoise ( Tseries_sum, noiseSigma, dataParams->randSeed ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  /*---------------------------------------------
   * turn this timeseries into SFTs, if requested
   *---------------------------------------------*/
  if ( SFTvect )
    {
      // Prepare windowing of time series for SFTs
      REAL4Window *window = NULL;
      if ( dataParams->SFTWindowType )
        {
          REAL8 dt = Tseries_sum->deltaT;
          UINT4 numTimesteps = round ( Tsft / dt );	/* number of time-samples in an Tsft (should be exact integer!) */

          XLAL_CHECK ( (window = XLALCreateNamedREAL4Window ( dataParams->SFTWindowType, dataParams->SFTWindowBeta, numTimesteps )) != NULL, XLAL_EFUNC );
        } // if uvar->SFTwindowName

      SFTParams sftParams = empty_SFTParams;
      sftParams.Tsft = Tsft;
      sftParams.timestamps = dataParams->timestamps;
      sftParams.noiseSFTs = NULL;	// not used here any more!
      sftParams.window = window;

      /* compute SFTs from timeseries */
      SFTVector *sftVect;
      XLAL_CHECK ( (sftVect = XLALSignalToSFTs (Tseries_sum, &sftParams)) != NULL, XLAL_EFUNC );

      // extract effective band from this, if neccessary (ie if faster-sampled output SFTs)
      if ( n1_fSamp != n0_fSamp )
        {
          XLAL_CHECK ( ((*SFTvect) = XLALExtractBandFromSFTVector ( sftVect, fMin, fBand )) != NULL, XLAL_EFUNC );
          XLALDestroySFTVector ( sftVect );
        }
      else
        {
          (*SFTvect) = sftVect;
        }

      XLALDestroyREAL4Window ( window );

    } // if SFTvect

  // return timeseries if requested
  if ( Tseries )
    {
      (*Tseries) = Tseries_sum;
    }
  else
    {
      XLALDestroyREAL4TimeSeries ( Tseries_sum );
    }

  return XLAL_SUCCESS;

} // XLALCWMakeFakeData()


/**
 * Generate a (heterodyned) REAL4 timeseries of a CW signal for given pulsarParams,
 * site, start-time, duration, and sampling-rate
 *
 * NOTE: this is mostly an API-wrapper to the more 'old-style' function
 * XLALGeneratePulsarSignal() [which will become deprecated in the future],
 * extended for the option to generate transient-CW signals
 */
REAL4TimeSeries *
XLALGenerateCWSignalTS ( const PulsarParams *pulsarParams,	///< input CW pulsar-signal parameters
                         const LALDetector *site,		///< detector
                         LIGOTimeGPS startTime,			///< time-series start-time GPS
                         REAL8 duration,
                         REAL8 fSamp,
                         REAL8 fHet,
                         EphemerisData *edat
                         )
{
  XLAL_CHECK_NULL ( pulsarParams != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( site != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( duration > 0, XLAL_EDOM );
  XLAL_CHECK_NULL ( fSamp > 0, XLAL_EDOM );
  XLAL_CHECK_NULL ( fHet >= 0, XLAL_EDOM );

  // translate amplitude params
  REAL8 h0     = pulsarParams->Amp.h0;
  REAL8 cosi   = pulsarParams->Amp.cosi;
  REAL8 aPlus  = 0.5 * h0 * ( 1.0 + SQ(cosi) );
  REAL8 aCross = h0 * cosi;
  // translate 'modern' fkdot into 'old-style' spindown-vector
  UINT4 s_max;
  for ( s_max = PULSAR_MAX_SPINS-1; s_max > 0; s_max -- )
    {
      if ( pulsarParams->Doppler.fkdot[s_max] != 0 )
        break;
    } // for s_max = max ... 0
  REAL8Vector *spindown = NULL;
  if ( s_max > 0 )
    {
      XLAL_CHECK_NULL ( (spindown = XLALCreateREAL8Vector ( s_max )) != NULL, XLAL_EFUNC );
      for ( UINT4 s = 0; s < s_max; s ++ ) {
        spindown->data[s] = pulsarParams->Doppler.fkdot[s+1];
      }
    }

  /*----------------------------------------
   * fill old-style PulsarSignalParams struct
   *----------------------------------------*/
  PulsarSignalParams params = empty_PulsarSignalParams;
  params.pulsar.refTime            = pulsarParams->Doppler.refTime;
  params.pulsar.position.system    = COORDINATESYSTEM_EQUATORIAL;
  params.pulsar.position.longitude = pulsarParams->Doppler.Alpha;
  params.pulsar.position.latitude  = pulsarParams->Doppler.Delta;
  params.pulsar.aPlus              = aPlus;
  params.pulsar.aCross             = aCross;
  params.pulsar.phi0               = pulsarParams->Amp.phi0;
  params.pulsar.psi                = pulsarParams->Amp.psi;
  params.pulsar.f0                 = pulsarParams->Doppler.fkdot[0];
  params.pulsar.spindown           = spindown;
  params.orbit                     = pulsarParams->Doppler.orbit;
  params.transfer                  = NULL;
  params.ephemerides               = edat;
  params.fHeterodyne               = fHet;

  // detector-specific settings
  params.startTimeGPS              = startTime;
  params.duration                  = ceil ( duration );
  params.samplingRate              = fSamp;
  params.site                      = site;

  /*----------------------------------------
   * generate the signal time-series
   *----------------------------------------*/
  REAL4TimeSeries *Tseries;
  XLAL_CHECK_NULL ( (Tseries = XLALGeneratePulsarSignal ( &params )) != NULL, XLAL_EFUNC );
  // ----- free internal memory
  XLALDestroyREAL8Vector ( spindown );

  // ----- apply transient-CW window
  XLAL_CHECK_NULL ( XLALApplyTransientWindow ( Tseries, pulsarParams->Transient ) == XLAL_SUCCESS, XLAL_EFUNC );

  return Tseries;

} // XLALGenerateCWSignalTS()


/**
 * Find the smallest sampling rate of the form fsamp = n / Tsft, with n>=n0,
 * such that all gap sizes Dt_i between SFTs of the given timestamps are also
 * exactly resolved, ie. that Dt_i * fsamp = integer, for all i
 *
 * The smallest allowed sampling rate is the user-specified fsamp0 = n0 / Tsft,
 * which guarantees by construction that fSamp0 * Tsft = n0 = integer
 * This sampling rate would be valid if there are no gaps between SFTs,
 * so it's only in cases of gaps that are non-integer multiples of Tsft that we'll
 * (potentially) have to increase the sampling rate.
 *
 * NOTE: This approach replaces the old mfdv4 habit of 'nudging' noise SFTs start-times
 * to fall on integer timesteps of the fsamp0 timeseries. The purpose of this function
 * is to avoid this behaviour, by appropriately increasing the sampling rate
 * as required.
 *
 * NOTE2: we only allow integer-second gaps, everything else will be
 * rejected with an error-message.
 *
 * NOTE3: Tsft=timestamps->deltaT must be integer seconds, everything else is rejected
 * with an error as well
 */
int
XLALFindSmallestValidSamplingRate ( UINT4 *n1,				//< [out] minimal valid sampling rate n1/Tsft
                                    UINT4 n0, 				//< [in] minimal sampling rate n0/Tsft
                                    const LIGOTimeGPSVector *timestamps	//< [in] start-timestamps and length Tsft of SFTs
                                    )
{
  XLAL_CHECK ( n1 != NULL, XLAL_EINVAL );
  XLAL_CHECK ( n0 > 0, XLAL_EINVAL );
  XLAL_CHECK ( timestamps && (timestamps->length > 0), XLAL_EINVAL );
  REAL8 TsftREAL = timestamps->deltaT;
  XLAL_CHECK ( TsftREAL == round(TsftREAL), XLAL_EDOM, "Only exact integer-second Tsft allowed, got Tsft = %g s\n", TsftREAL );
  UINT4 Tsft = (UINT4)TsftREAL;
  XLAL_CHECK ( Tsft > 0, XLAL_EINVAL );

  // NOTE: all 'valid' sampling rates are of the form  fSamp = n / Tsft, where n >= n0,
  // therefore we label sampling rates by their index 'n'. The purpose is to find the
  // smallest 'valid' n, which is the max of the required 'n' over all gaps
  UINT4 nCur = n0;

  // ----- We just step through the vector and figure out for each gap if we need to
  // decrease the stepsize Tsft/n from the previous value
  UINT4 numSFTs = timestamps->length;

  for ( UINT4 i = 1; i < numSFTs; i ++ )
    {
      LIGOTimeGPS *t1 = &(timestamps->data[i]);
      LIGOTimeGPS *t0 = &(timestamps->data[i-1]);

      INT4 nsdiff = t1->gpsNanoSeconds - t0->gpsNanoSeconds;
      XLAL_CHECK ( nsdiff == 0, XLAL_EDOM, "Only integer-second gaps allowed, found %d ns excess in gap between i=%d and i-1\n", nsdiff, i );

      INT4 gap_i0 = t1->gpsSeconds - t0->gpsSeconds;
      XLAL_CHECK ( gap_i0 > 0, XLAL_EDOM, "Timestamps must be sorted in increasing order, found negative gap %d s between i=%d and i-1\n", gap_i0, i );

      // now reduce gap to remainder wrt Tsft
      INT4 gap_i = gap_i0 % Tsft;

      XLALPrintInfo ("gap_i = %d s, remainder wrt Tsft=%d s = %d s: ", gap_i0, Tsft, gap_i );
      if ( (gap_i * nCur) % Tsft == 0 ) {
        XLALPrintInfo ("Fits exactly with fSamp = %d / Tsft = %g\n", nCur, nCur / TsftREAL );
        continue;
      }

      // otherwise:
      // solve for required new smaller step-size 'dt = Tsft/nNew' to fit integer cycles
      // both into Tsft (by construction) and into (gap_i mod Tsft), with nNew > nCur
      //
      // gap[i] == (t[i+1] - t[i]) mod Tsft
      // such that 0 <= gap[i] < Tsft
      // we want integers nNew, m , such that nNew > nCur, and m < nNew, with
      // nNew * dtNew = Tsft   AND  m * dtNew = gap[i]
      // ==> m / nNew = gap[i] / Tsft
      // This could be solved easily by rounding nCur to next highest
      // multiple of Tsft: nNew' = ceil(nCur/Tsft) * Tsft > nCur
      // but this can be wasteful if the fraction simplifies: so we compute
      // the greatest common divisor 'g'
      // g = gcd(gap[i], Tsft), and then use
      // nNew = ceil ( nCur  * g  / Tsft ) * Tsft / g > nCur

      UINT4 g = gcd ( gap_i, Tsft );
      REAL8 Tg = TsftREAL / g;
      UINT4 nNew = (UINT4) ceil ( nCur / Tg ) * Tg;

      XLAL_CHECK ( nNew > nCur, XLAL_ETOL, "This seems wrong: nNew = %d !> nCur = %d, but should be greater!\n", nNew, nCur );
      XLALPrintInfo ("Need to increase to fSamp = %d / Tsft = %g\n", nNew, nNew / TsftREAL );

      nCur = nNew;

    } // for i < numSFTs


  // our final minimal valid sampling rate is therefore n1/Tsft
  (*n1) = nCur;

  return XLAL_SUCCESS;

} // XLALFindSmallestValidSamplingRate()


/* Find greatest common divisor between two numbers,
 * where numer <= denom
 * this is an implementation of the Euclidean Algorithm,
 * taken from John's UnitNormalize.c, and extended to UINT4's
 *
 * For reference, see also
 * https://en.wikipedia.org/wiki/Euclidean_algorithm#Implementations
 */
static UINT4
gcd (UINT4 numer, UINT4 denom)
{
  UINT4 next_numer, next_denom, remainder;

  next_numer = numer;
  next_denom = denom;
  while ( next_denom != 0 )
    {
      remainder = next_numer % next_denom;
      next_numer = next_denom;
      next_denom = remainder;
    }
  return next_numer;
} // gcd

/**
 * Simple destructor for an 'InjectionSources' object.
 */
void
XLALDestroyInjectionSources ( InjectionSources *sources )
{
  if ( sources == NULL ) {
    return;
  }
  if ( sources->pulsarParams != NULL )
    {
      for ( UINT4 i = 0 ; i < sources->numPulsars; i ++ )
        {
          BinaryOrbitParams *orbit = sources->pulsarParams[i].Doppler.orbit;
          if ( orbit != NULL ) {
            XLALFree ( orbit );
          }
        } // for i < numPulsars
      XLALFree ( sources->pulsarParams );
    }

  XLALFree ( sources );

  return;

} /* XLALDestroyInjectionSources() */

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

