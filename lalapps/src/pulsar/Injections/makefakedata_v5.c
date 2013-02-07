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
 * File Name: makefakedata_v4.c
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
/** configuration-variables derived from user-variables */
typedef struct
{
  PulsarParams pulsar;		/**< pulsar signal-parameters (amplitude + doppler */

  EphemerisData *edat;		/**< ephemeris-data */
  LALDetector site;  		/**< detector-site info */

  LIGOTimeGPS startTimeGPS;	/**< start-time of observation */
  UINT4 duration;		/**< total duration of observation in seconds */

  LIGOTimeGPSVector *timestamps;/**< a vector of timestamps to generate time-series/SFTs for */

  SFTCatalog *noiseCatalog;/**< catalog of noise-SFTs to be added to signal */

  COMPLEX8FrequencySeries *transfer;  /**< detector's transfer function for use in hardware-injection */

  INT4 randSeed;		/**< random-number seed: either taken from user or /dev/urandom */

  transientWindow_t transientWindow;	/**< properties of transient-signal window */
  CHAR *VCSInfoString;          /**< LAL + LALapps Git version string */
} ConfigVars_t;

// ----- User variables
typedef struct
{
  BOOLEAN help;		/**< Print this help/usage message */

  /* output */
  CHAR *outSFTbname;		/**< Path and basefilename of output SFT files */
  BOOLEAN outSingleSFT;	/**< use to output a single concatenated SFT */

  CHAR *TDDfile;		/**< Filename for ASCII output time-series */
  BOOLEAN hardwareTDD;	/**< Binary output timeseries in chunks of Tsft for hardware injections. */

  CHAR *logfile;		/**< name of logfile */

  /* specify start + duration */
  CHAR *timestampsFile;	/**< Timestamps file */
  INT4 startTime;		/**< Start-time of requested signal in detector-frame (GPS seconds) */
  INT4 duration;		/**< Duration of requested signal in seconds */

  /* time-series sampling + heterodyning frequencies */
  REAL8 fmin;		/**< Lowest frequency in output SFT (= heterodyning frequency) */
  REAL8 Band;		/**< bandwidth of output SFT in Hz (= 1/2 sampling frequency) */

  /* SFT params */
  REAL8 Tsft;		/**< SFT time baseline Tsft */
  REAL8 SFToverlap;		/**< overlap SFTs by this many seconds */

  /* noise to add [OPTIONAL] */
  CHAR *noiseSFTs;		/**< Glob-like pattern specifying noise-SFTs to be added to signal */
  REAL8 noiseSqrtSh;		/**< ALTERNATIVE: single-sided sqrt(Sh) for Gaussian noise */

  /* Window function [OPTIONAL] */
  CHAR *windowType;		/**< Windowing function for the time series */
  REAL8 windowBeta;          	/**< 'beta' parameter required for certain window-types */

  /* Detector and ephemeris */
  CHAR *IFO;			/**< Detector: H1, L1, H2, V1, ... */

  CHAR *actuation;		/**< filename containg detector actuation function */
  REAL8 actuationScale;	/**< Scale-factor to be applied to actuation-function */
  CHAR *ephemDir;		/**< Directory path for ephemeris files (optional), use LAL_DATA_PATH if unset. */
  CHAR *ephemYear;		/**< Year (or range of years) of ephemeris files to be used */

  /* pulsar parameters [REQUIRED] */
  REAL8 refTime;		/**< Pulsar reference time tRef in SSB ('0' means: use startTime converted to SSB) */
  REAL8 refTimeMJD;          /**< Pulsar reference time tRef in MJD ('0' means: use startTime converted to SSB) */

  REAL8 h0;			/**< overall signal amplitude h0 */
  REAL8 cosi;		/**< cos(iota) of inclination angle iota */
  REAL8 aPlus;		/**< ALTERNATIVE to {h0,cosi}: Plus polarization amplitude aPlus */
  REAL8 aCross;		/**< ALTERNATIVE to {h0,cosi}: Cross polarization amplitude aCross */
  REAL8 psi;			/**< Polarization angle psi */
  REAL8 phi0;		/**< Initial phase phi */

  REAL8 Alpha;		/**< Right ascension [radians] alpha of pulsar */
  REAL8 Delta;		/**< Declination [radians] delta of pulsar */
  CHAR *RA;		        /**< Right ascension [hh:mm:ss.ssss] alpha of pulsar */
  CHAR *Dec;	         	/**< Declination [dd:mm:ss.ssss] delta of pulsar */

  REAL8 Freq;

  REAL8 f1dot;		/**< First spindown parameter f' */
  REAL8 f2dot;		/**< Second spindown parameter f'' */
  REAL8 f3dot;		/**< Third spindown parameter f''' */

  /* orbital parameters [OPTIONAL] */

  REAL8 orbitasini;	        /**< Projected orbital semi-major axis in seconds (a/c) */
  REAL8 orbitEcc;	        /**< Orbital eccentricity */
  INT4  orbitTpSSBsec;	/**< 'observed' (SSB) time of periapsis passage. Seconds. */
  INT4  orbitTpSSBnan;	/**< 'observed' (SSB) time of periapsis passage. Nanoseconds. */
  REAL8 orbitTpSSBMJD;       /**< 'observed' (SSB) time of periapsis passage. MJD. */
  REAL8 orbitPeriod;		/**< Orbital period (seconds) */
  REAL8 orbitArgp;	        /**< Argument of periapsis (radians) */

  BOOLEAN lineFeature;	/**< generate a monochromatic line instead of a pulsar-signal */

  BOOLEAN version;		/**< output version information */

  INT4 randSeed;		/**< allow user to specify random-number seed for reproducible noise-realizations */

  CHAR *parfile;             /** option .par file path */
  CHAR *transientWindowType;	/**< name of transient window ('rect', 'exp',...) */
  REAL8 transientStartTime;	/**< GPS start-time of transient window */
  REAL8 transientTauDays;	/**< time-scale in days of transient window */
} UserVariables_t;

// ----- global variables ----------

// ----- empty structs for initializations
static const UserVariables_t empty_UserVariables;
static const ConfigVars_t empty_GV;
static const LALUnit empty_LALUnit;

// ---------- exportable API prototypes ----------
int XLALFindSmallestValidSamplingRate ( REAL8 *fSamp, UINT4 n0, const LIGOTimeGPSVector *timestamps );
int XLALMakeFakeCWData ( SFTVector **SFTs, REAL4TimeSeries **Tseries, const UserVariables_t *uvar, const ConfigVars_t *cfg );

// ---------- local prototypes ----------
int XLALInitUserVars ( UserVariables_t *uvar, int argc, char *argv[] );
int XLALInitMakefakedata ( ConfigVars_t *cfg, UserVariables_t *uvar );

int XLALWriteMFDlog ( const char *logfile, const ConfigVars_t *cfg );
COMPLEX8FrequencySeries *XLALLoadTransferFunctionFromActuation ( REAL8 actuationScale, const CHAR *fname );
int XLALFreeMem ( ConfigVars_t *cfg );

extern void write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series);
extern void write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series);
BOOLEAN is_directory ( const CHAR *fname );
static UINT4 gcd (UINT4 numer, UINT4 denom);

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  int ret, len;
  ConfigVars_t GV = empty_GV;
  FILE *fpSingleSFT = NULL;

  UserVariables_t uvar = empty_UserVariables;

  /* ------------------------------
   * read user-input and set up shop
   *------------------------------*/
  XLAL_CHECK ( XLALInitUserVars ( &uvar, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALInitMakefakedata ( &GV, &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  SFTVector *SFTs = NULL;
  REAL4TimeSeries *Tseries = NULL;

  XLAL_CHECK ( XLALMakeFakeCWData ( &SFTs, &Tseries, &uvar, &GV ) == XLAL_SUCCESS, XLAL_EFUNC );


  // if noiseSFTs specified, load them and add them to the resulting SFT-vector
  if ( GV.noiseCatalog )
    {
      /* load effective frequency-band from noise-SFTs */
      UINT4 numBins = SFTs->data[0].data->length;
      REAL8 dFreq   = SFTs->data[0].deltaF;
      REAL8 fMin    = SFTs->data[0].f0;
      REAL8 fMax    = fMin + ( numBins - 1 ) * dFreq;
      SFTVector *noiseSFTs;
      XLAL_CHECK ( (noiseSFTs = XLALLoadSFTs ( GV.noiseCatalog, fMin, fMax )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( XLALSFTVectorAdd ( SFTs, noiseSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroySFTVector ( noiseSFTs );
    }

  /* if user requesting single concatenated SFT */
  if ( uvar.outSingleSFT )
    {
      /* check that user isn't giving a directory */
      if ( is_directory ( uvar.outSFTbname ) ) {
        XLAL_ERROR ( XLAL_ETYPE, "'%s' is a directory, but --outSingleSFT expects a filename!\n", uvar.outSFTbname);
      }

      /* open concatenated SFT file for writing */
      if ( (fpSingleSFT = LALFopen ( uvar.outSFTbname, "wb" )) == NULL ) {
        XLAL_ERROR ( XLAL_EIO, "Failed to open file '%s' for writing: %s\n\n", uvar.outSFTbname, strerror(errno));
      }
    } // if outSingleSFT


  /* output ASCII time-series if requested */
  if ( uvar.TDDfile )
    {
      FILE *fp;
      if ( (fp = fopen (uvar.TDDfile, "w")) == NULL)
        {
          perror ("Error opening outTDDfile for writing");
          XLAL_ERROR ( XLAL_EIO, "Failed to fopen TDDfile = '%s' for writing\n", uvar.TDDfile );
        }

      write_timeSeriesR4(fp, Tseries);
      fclose(fp);
    } /* if outputting ASCII time-series */


  /* for HARDWARE-INJECTION:
   * before the first chunk we send magic number and chunk-length to stdout
   */
  if ( uvar.hardwareTDD )
    {
      REAL4 magic = 1234.5;
      UINT4 length = Tseries->data->length;
      if ( (1 != fwrite ( &magic, sizeof(magic), 1, stdout )) || (1 != fwrite(&length, sizeof(INT4), 1, stdout)) )
        {
          perror ("Failed to write to stdout");
          XLAL_ERROR ( XLAL_EIO, "fwrite() failed\n");
        }

      REAL4 *datap = Tseries->data->data;

      if ( length != fwrite (datap, sizeof(datap[0]), length, stdout) )
        {
          perror( "Fatal error in writing binary data to stdout\n");
          XLAL_ERROR ( XLAL_EIO, "fwrite() failed\n");
        }
      fflush (stdout);

    } /* if hardware injections */


  if (uvar.outSFTbname)
    {
      /* generate comment string */
      CHAR *logstr;
      XLAL_CHECK ( (logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
      char *comment = XLALCalloc ( 1, len = strlen ( logstr ) + strlen(GV.VCSInfoString) + 512 );
      XLAL_CHECK ( comment != NULL, XLAL_ENOMEM, "XLALCalloc(1,%d) failed.\n", len );
      sprintf ( comment, "Generated by:\n%s\n%s\n", logstr, GV.VCSInfoString );

      /* if user requesting single concatenated SFT */
      if ( uvar.outSingleSFT )
        {
          /* write all SFTs to concatenated file */
          for ( UINT4 k = 0; k < SFTs->length; k++ )
            {
              ret = XLALWriteSFT2fp ( &(SFTs->data[k]), fpSingleSFT, comment );
              XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALWriteSFT2fp() failed to write SFT %d / %d to '%s'!\n", k+1, SFTs->length, uvar.outSFTbname );
            } // for k < numSFTs
        } // if outSingleSFT
      else
        {
          /* check that user didn't follow v1-habits and tried to give us a base filename */
          if ( ! is_directory ( uvar.outSFTbname ) )
            {
              fprintf (stderr, "\nERROR: the --outSFTbname '%s' is not a directory!\n", uvar.outSFTbname);
              fprintf (stderr, "------------------------------------------------------------\n");
              fprintf (stderr, "NOTE: v2-SFTs are written using the SFTv2 filename-convention!\n");
              fprintf (stderr, "==> therefore, only an (existing) directory-name may be given to --outSFTbname, but no basename!\n");
              fprintf (stderr, "------------------------------------------------------------\n\n");
              XLAL_ERROR ( XLAL_EINVAL );
            }
          XLAL_CHECK ( XLALWriteSFTVector2Dir( SFTs, uvar.outSFTbname, comment, "mfdv4" ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
      XLALFree ( logstr );
      XLALFree ( comment );

    } /* if outSFTbname */

  /* free memory */
  if (Tseries)
    {
      XLALDestroyREAL4TimeSeries ( Tseries );
      Tseries = NULL;
    }

  if (SFTs)
    {
      XLALDestroySFTVector ( SFTs );
      SFTs = NULL;
    }

  /* if user requesting single concatenated SFT */
  if ( uvar.outSingleSFT ) {

    /* close concatenated SFT */
    LALFclose( fpSingleSFT );

  }

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

  BOOLEAN have_parfile = XLALUserVarWasSet (&uvar->parfile);
  BinaryPulsarParams pulparams;

  /* read in par file parameters if given */
   if (have_parfile)
     {
       XLALReadTEMPOParFile( &pulparams, uvar->parfile);
       XLAL_CHECK ( xlalErrno == XLAL_SUCCESS, XLAL_EFUNC, "XLALReadTEMPOParFile() failed for parfile = '%s', xlalErrno = %d\n", uvar->parfile, xlalErrno );
       XLAL_CHECK ( pulparams.f0 > 0, XLAL_EINVAL, "Invalid .par file values, need f0 > 0!\n" );
       XLAL_CHECK ( (pulparams.pepoch > 0) || (pulparams.posepoch > 0), XLAL_EINVAL, "Invalid .par file values, need PEPOCH or POSEPOCH!\n");
     }

  /* if requested, log all user-input and code-versions */
  if ( uvar->logfile ) {
    XLAL_CHECK ( XLALWriteMFDlog ( uvar->logfile, cfg ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALWriteMFDlog() failed with xlalErrno = %d\n", xlalErrno );
  }

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

    /*check .par file for gw parameters*/
    if (have_parfile){
      if (pulparams.h0 != 0){
	uvar->h0 = pulparams.h0;
	uvar->cosi = pulparams.cosiota;
	uvar->phi0 = pulparams.phi0;
	uvar->psi = pulparams.psi;
	have_h0 = 1; /*Set to TRUE as uvar->h0 not declared on command line -- same for rest*/
	have_cosi = 1;
      }
      else{
	uvar->aPlus = pulparams.Aplus;
	uvar->aCross = pulparams.Across;
	uvar->phi0 = pulparams.phi0;
	uvar->psi = pulparams.psi;
	have_aPlus = 1;
	have_aCross = 1;
      }
      uvar->Freq = 2.*pulparams.f0;
      uvar->Alpha = pulparams.ra;
      uvar->Delta = pulparams.dec;
      have_Freq = 1;
      have_Alpha = 1;
      have_Delta = 1;
    }

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
	cfg->pulsar.Amp.h0 = uvar->h0;
	cfg->pulsar.Amp.cosi = uvar->cosi;
      }
    else if ( have_aPlus && have_aCross )
      {  /* translate A_{+,x} into {h_0, cosi} */
	REAL8 disc;
	if ( fabs(uvar->aCross) > uvar->aPlus ) {
          XLAL_ERROR ( XLAL_EINVAL, "Invalid input parameters: |aCross| = %g must be <= than aPlus = %g.\n", fabs(uvar->aCross), uvar->aPlus );
        }
	disc = sqrt ( SQ(uvar->aPlus) - SQ(uvar->aCross) );
	cfg->pulsar.Amp.h0   = uvar->aPlus + disc;
        if ( cfg->pulsar.Amp.h0 > 0 )
          cfg->pulsar.Amp.cosi = uvar->aCross / cfg->pulsar.Amp.h0;	// avoid division by 0!
        else
          cfg->pulsar.Amp.cosi = 0;
      }
    else {
      cfg->pulsar.Amp.h0 = 0.0;
      cfg->pulsar.Amp.cosi = 0.0;
    }
    cfg->pulsar.Amp.phi0 = uvar->phi0;
    cfg->pulsar.Amp.psi  = uvar->psi;

    /* ----- signal Frequency ----- */
    if ( have_Freq )
      cfg->pulsar.Doppler.fkdot[0] = uvar->Freq;
    else
      cfg->pulsar.Doppler.fkdot[0] = 0.0;

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
	cfg->pulsar.Doppler.Alpha = uvar->Alpha;
	cfg->pulsar.Doppler.Delta = uvar->Delta;
      }
    else if ( have_RA )
      {
	cfg->pulsar.Doppler.Alpha = LALDegsToRads(uvar->RA,"alpha");
	cfg->pulsar.Doppler.Delta = LALDegsToRads(uvar->Dec,"delta");
      }

  } /* Pulsar signal parameters */

  /* ---------- prepare vector of spindown parameters ---------- */
  {
    if ( have_parfile )
      {
	uvar->f1dot = 2.*pulparams.f1;
	uvar->f2dot = 2.*pulparams.f2;
	uvar->f3dot = 2.*pulparams.f3;
      }
    cfg->pulsar.Doppler.fkdot[1] = uvar->f1dot;
    cfg->pulsar.Doppler.fkdot[2] = uvar->f2dot;
    cfg->pulsar.Doppler.fkdot[3] = uvar->f3dot;

  } /* END: prepare spindown parameters */


   /* check for negative fmin and Band, which would break the fmin_eff, fBand_eff calculation below */
  XLAL_CHECK ( uvar->fmin >= 0, XLAL_EDOM, "Invalid negative frequency fmin=%f!\n\n", uvar->fmin );
  XLAL_CHECK ( uvar->Band >= 0, XLAL_EDOM, "Invalid negative frequency band Band=%f!\n\n", uvar->Band );

  /* ---------- determine timestamps to produce signal for  ---------- */
  {
    /* check input consistency: *uvar->timestampsFile, uvar->startTime, uvar->duration */
    BOOLEAN haveStart, haveDuration, haveTimestampsFile, haveOverlap;

    haveStart = XLALUserVarWasSet(&uvar->startTime);
    haveDuration = XLALUserVarWasSet(&uvar->duration);
    haveTimestampsFile = ( uvar->timestampsFile != NULL );
    haveOverlap = ( uvar->SFToverlap > 0 );

    if ( ( haveDuration && !haveStart) || ( !haveDuration && haveStart ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "Need BOTH --startTime AND --duration if you give one of them !\n\n");
    }

    /* don't allow using --SFToverlap with anything other than pure (--startTime,--duration) */
    if ( haveOverlap && ( uvar->noiseSFTs || haveTimestampsFile ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "I can't combine --SFToverlap with --noiseSFTs or --timestampsFile, only use with (--startTime, --duration)!\n\n");
    }

    /*-------------------- check special case: Hardware injection ---------- */
    /* don't allow timestamps-file, noise-SFTs or SFT-output */
    if ( uvar->hardwareTDD )
      {
	if (haveTimestampsFile || uvar->noiseSFTs ) {
          XLAL_ERROR ( XLAL_EINVAL, "\nHardware injection: don't specify --timestampsFile or --noiseSFTs\n\n");
        }
	if ( !haveStart || !haveDuration ) {
          XLAL_ERROR ( XLAL_EINVAL, "Hardware injection: need to specify --startTime and --duration!\n\n");
        }
	if ( XLALUserVarWasSet( &uvar->outSFTbname ) ) {
          XLAL_ERROR ( XLAL_EINVAL, "Hardware injection mode is incompatible with producing SFTs\n\n");
        }
      } /* if hardware-injection */

    /* ----- load timestamps from file if given  */
    if ( haveTimestampsFile )
      {
	if ( haveStart || haveDuration || haveOverlap ) {
          XLAL_ERROR ( XLAL_EINVAL, "Using --timestampsFile is incompatible with either of --startTime, --duration or --SFToverlap\n\n");
        }
	if ( ( cfg->timestamps = XLALReadTimestampsFile ( uvar->timestampsFile )) == NULL ) {
          XLAL_ERROR ( XLAL_EFUNC, "XLALReadTimestampsFile() failed to read timestamps file '%s'\n", uvar->timestampsFile );
        }
	if ( ( cfg->timestamps->length == 0 ) || ( cfg->timestamps->data == NULL ) ) {
          XLAL_ERROR ( XLAL_EINVAL, "Got empty timestamps-list from file '%s'\n", uvar->timestampsFile );
        }
        cfg->timestamps->deltaT = uvar->Tsft;

        XLAL_CHECK ( uvar->noiseSFTs == NULL, XLAL_EINVAL, "--timestampsFile is incompatible with --noiseSFTs\n" );

      } /* if haveTimestampsFile */

    /* ----- if real noise-SFTs given: load them now using EITHER (start,start+duration) OR timestamps
     * as constraints if given, otherwise load all of them. Also require window option to be given
     */
    if ( uvar->noiseSFTs )
      {
	SFTConstraints constraints = empty_SFTConstraints;
	LIGOTimeGPS minStartTime, maxEndTime;
        BOOLEAN have_window = XLALUserVarWasSet ( &uvar->windowType );

        /* user must specify the window function used for the noiseSFTs */
        if ( !have_window ) {
          XLAL_ERROR ( XLAL_EINVAL, "Require window option to be given when specifying noiseSFTs.\n" );
        }

        XLALPrintWarning ( "\nWARNING: only SFTs corresponding to the noiseSFT-timestamps will be produced!\n" );

	/* use all additional constraints user might have given */
	if ( haveStart && haveDuration )
	  {
	    XLALGPSSetREAL8 ( &minStartTime, uvar->startTime );
	    constraints.startTime = &minStartTime;
	    XLALGPSSetREAL8 ( &maxEndTime, uvar->startTime + uvar->duration );
	    constraints.endTime = &maxEndTime;
            XLALPrintWarning ( "\nWARNING: only noise-SFTs between GPS [%d, %d] will be used!\n", uvar->startTime, uvar->startTime + uvar->duration );
	  } /* if start+duration given */

        CHAR *channelName = NULL;
        XLAL_CHECK ( (channelName = XLALGetChannelPrefix ( uvar->IFO )) != NULL, XLAL_EFUNC, "XLALGetChannelPrefix('%s') failed.\n", uvar->IFO );
	/* use detector-constraint */
	constraints.detector = channelName ;

	XLAL_CHECK ( (cfg->noiseCatalog = XLALSFTdataFind ( uvar->noiseSFTs, &constraints )) != NULL, XLAL_EFUNC );
	XLAL_CHECK ( cfg->noiseCatalog->length > 0, XLAL_EINVAL, "No noise-SFTs matching the constraints (IFO, start+duration, timestamps) were found!\n" );
        XLALFree ( channelName );

	/* extract timestamps from the SFT-catalog */
        XLAL_CHECK ( (cfg->timestamps = XLALTimestampsFromSFTCatalog ( cfg->noiseCatalog )) != NULL, XLAL_EFUNC );

      } /* if uvar->noiseSFTs */

    /* have we got our timestamps yet?: If not, we must get them from (start, duration) user-input */
    if ( ! cfg->timestamps )
      {
	REAL8 tStep = uvar->Tsft - uvar->SFToverlap;
	LIGOTimeGPS tStart;
	REAL8 t0, tLast;
	if ( !haveStart || !haveDuration ) {
          XLAL_ERROR ( XLAL_EINVAL, "Need to have either --timestampsFile OR (--startTime,--duration) OR --noiseSFTs\n\n");
        }
	if ( uvar->SFToverlap > uvar->Tsft ) {
          XLAL_ERROR ( XLAL_EINVAL, "--SFToverlap cannot be larger than --Tsft!\n\n");
        }

	/* internally always use timestamps, so we generate them  */
	XLALGPSSetREAL8 ( &tStart, uvar->startTime );

        cfg->timestamps = XLALMakeTimestamps ( tStart, uvar->duration, tStep );
        XLAL_CHECK ( cfg->timestamps != NULL, XLAL_EFUNC );

	/* "prune" last timestamp(s) if the one before-last also covers the end
	 * (this happens for overlapping SFTs with LALMakeTimestamps() as used above
	 */
	t0 = XLALGPSGetREAL8( &(cfg->timestamps->data[0]) );
	tLast = XLALGPSGetREAL8 ( &(cfg->timestamps->data[cfg->timestamps->length - 1 ]) );
	while ( tLast - t0  + uvar->Tsft > uvar->duration + 1e-6)
	  {
	    cfg->timestamps->length --;
	    tLast = XLALGPSGetREAL8 ( &(cfg->timestamps->data[cfg->timestamps->length - 1 ]) );
	  }

      } /* if !cfg->timestamps */

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
    if (have_parfile){
      if (pulparams.model != NULL) {
	uvar->orbitasini = pulparams.x;
	uvar->orbitPeriod = pulparams.Pb*86400;
	if (strstr(pulparams.model,"ELL1") != NULL) {
	  REAL8 w,e,eps1,eps2;
	  eps1 = pulparams.eps1;
	  eps2 = pulparams.eps2;
	  w = atan2(eps1,eps2);
	  e = sqrt(eps1*eps1+eps2*eps2);
	  uvar->orbitArgp = w;
	  uvar->orbitEcc = e;
	}
	else {
	  uvar->orbitArgp = pulparams.w0;
	  uvar->orbitEcc = pulparams.e;
	}
	if (strstr(pulparams.model,"ELL1") != NULL) {
	  REAL8 fe, uasc,Dt;
	  fe = sqrt((1.0-uvar->orbitEcc)/(1.0+uvar->orbitEcc));
	  uasc = 2.0*atan(fe*tan(uvar->orbitArgp/2.0));
	  Dt = (uvar->orbitPeriod/LAL_TWOPI)*(uasc-uvar->orbitEcc*sin(uasc));
	  pulparams.T0 = pulparams.Tasc + Dt;
	}
	uvar->orbitTpSSBsec = (UINT4)floor(pulparams.T0);
	uvar->orbitTpSSBnan = (UINT4)floor((pulparams.T0 - uvar->orbitTpSSBsec)*1e9);
      }
    }
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
	  GPSfloat = LALTTMJDtoGPS(uvar->orbitTpSSBMJD);
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

      cfg->pulsar.Doppler.orbit = orbit;     /* struct copy */
    } /* if one or more orbital parameters were set */
    else
      cfg->pulsar.Doppler.orbit = NULL;
  } /* END: binary orbital params */

  /* ----- set "pulsar reference time", i.e. SSB-time at which pulsar params are defined ---------- */
  if (XLALUserVarWasSet (&uvar->parfile)) {
    uvar->refTime = pulparams.pepoch; /*XLALReadTEMPOParFile already converted pepoch to GPS*/
    XLALGPSSetREAL8(&(cfg->pulsar.Doppler.refTime),uvar->refTime);
  }
  else if (XLALUserVarWasSet(&uvar->refTime) && XLALUserVarWasSet(&uvar->refTimeMJD))
    {
      XLAL_ERROR ( XLAL_EINVAL, "\nUse only one of '--refTime' and '--refTimeMJD' to specify SSB reference time!\n\n");
    }
  else if (XLALUserVarWasSet(&uvar->refTime))
    {
      XLALGPSSetREAL8(&(cfg->pulsar.Doppler.refTime), uvar->refTime);
    }
  else if (XLALUserVarWasSet(&uvar->refTimeMJD))
    {

      /* convert MJD to GPS using Matt Pitkins code found at lal/packages/pulsar/src/BinaryPulsarTimeing.c */
      REAL8 GPSfloat;
      GPSfloat = LALTTMJDtoGPS(uvar->refTimeMJD);
      XLALGPSSetREAL8(&(cfg->pulsar.Doppler.refTime),GPSfloat);
    }
  else
    cfg->pulsar.Doppler.refTime = cfg->timestamps->data[0];	/* internal startTime always found in here*/


  /* ---------- has the user specified an actuation-function file ? ---------- */
  if ( uvar->actuation )
    {
      /* currently we only allow using a transfer-function for hardware-injections */
      if (!uvar->hardwareTDD )
	{
	  XLAL_ERROR ( XLAL_EINVAL, "Error: use of an actuation/transfer function restricted to hardare-injections\n\n");
	}
      else
        {
          cfg->transfer = XLALLoadTransferFunctionFromActuation( uvar->actuationScale, uvar->actuation );
          XLAL_CHECK ( cfg->transfer != NULL, XLAL_EFUNC );
        }

    } /* if uvar->actuation */

  if ( !uvar->actuation && XLALUserVarWasSet(&uvar->actuationScale) ) {
    XLAL_ERROR ( XLAL_EINVAL, "Actuation-scale was specified without actuation-function file!\n\n");
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

  uvar->actuationScale = - 1.0;	/* seems to be the LIGO-default */

#define DEFAULT_TRANSIENT "none"
  uvar->transientWindowType = XLALCalloc ( 1, len = strlen(DEFAULT_TRANSIENT)+1 );
  XLAL_CHECK ( uvar->transientWindowType != NULL, XLAL_ENOMEM, "XLALCalloc ( 1, %d ) failed.\n", len );
  strcpy ( uvar->transientWindowType, DEFAULT_TRANSIENT );

  // ---------- register all our user-variable ----------
  XLALregBOOLUserStruct (  help,                'h', UVAR_HELP    , "Print this help/usage message");

  /* output options */
  XLALregBOOLUserStruct (   outSingleSFT,       's', UVAR_OPTIONAL, "Write a single concatenated SFT (name given by --outSFTbname)" );
  XLALregSTRINGUserStruct ( outSFTbname,        'n', UVAR_OPTIONAL, "Output SFTs: target Directory (if --outSingleSFT=false) or filename (if --outSingleSFT=true)");

  XLALregSTRINGUserStruct ( TDDfile,            't', UVAR_OPTIONAL, "Filename to output time-series into");

  XLALregSTRINGUserStruct ( logfile,            'l', UVAR_OPTIONAL, "Filename for log-output");

  /* detector and ephemeris */
  XLALregSTRINGUserStruct ( IFO,                'I', UVAR_OPTIONAL, "Detector: one of 'G1','L1','H1,'H2','V1', ...");

  XLALregSTRINGUserStruct ( ephemDir,           'E', UVAR_OPTIONAL, "Directory path for ephemeris files (use LAL_DATA_PATH if unspecified)");
  XLALregSTRINGUserStruct ( ephemYear,          'y', UVAR_OPTIONAL, "Year-range string of ephemeris files to be used");

  /* start + duration of timeseries */
  XLALregINTUserStruct (  startTime,            'G', UVAR_OPTIONAL, "Start-time of requested signal in detector-frame (GPS seconds)");
  XLALregINTUserStruct (  duration,              0,  UVAR_OPTIONAL, "Duration of requested signal in seconds");
  XLALregSTRINGUserStruct ( timestampsFile,      0,  UVAR_OPTIONAL, "ALTERNATIVE: File to read timestamps from (file-format: lines with <seconds> <nanoseconds>)");

  /* sampling and heterodyning frequencies */
  XLALregREALUserStruct (  fmin,                 0, UVAR_OPTIONAL, "Lowest frequency in output SFT (= heterodyning frequency)");
  XLALregREALUserStruct (  Band,                 0, UVAR_OPTIONAL, "Bandwidth of output SFT in Hz (= 1/2 sampling frequency)");

  /* SFT properties */
  XLALregREALUserStruct (  Tsft,                 0, UVAR_OPTIONAL, "Time baseline of one SFT in seconds");
  XLALregREALUserStruct (  SFToverlap,           0, UVAR_OPTIONAL, "Overlap between successive SFTs in seconds (conflicts with --noiseSFTs or --timestampsFile)");
  XLALregSTRINGUserStruct( windowType,           0, UVAR_OPTIONAL, "Window function to be applied to the SFTs (required when using --noiseSFTs)");
  XLALregREALUserStruct (  windowBeta,           0, UVAR_OPTIONAL, "Window 'beta' parameter required for a few window-types (eg. 'tukey')");

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
  XLALregSTRINGUserStruct ( noiseSFTs,          'D', UVAR_OPTIONAL, "Noise-SFTs to be added to signal (Uses ONLY SFTs falling within (--startTime,--duration) or the set given in --timstampsFile)");
  XLALregREALUserStruct (  noiseSqrtSh,          0, UVAR_OPTIONAL,  "ALTERNATIVE: Gaussian noise with single-sided PSD sqrt(Sh)");

  XLALregBOOLUserStruct (  lineFeature,          0, UVAR_OPTIONAL, "Generate a monochromatic 'line' of amplitude h0 and frequency 'Freq'}");

  XLALregBOOLUserStruct (  version,             'V', UVAR_SPECIAL, "Output version information");

  XLALregSTRINGUserStruct (parfile,             'p', UVAR_OPTIONAL, "Directory path for optional .par files");            /*registers .par file in mfd*/

  /* transient signal window properties (name, start, duration) */
  XLALregSTRINGUserStruct (transientWindowType,  0, UVAR_OPTIONAL, "Type of transient signal window to use. ('none', 'rect', 'exp').");
  XLALregREALUserStruct (  transientStartTime,   0, UVAR_OPTIONAL, "GPS start-time 't0' of transient signal window.");
  XLALregREALUserStruct (  transientTauDays,     0, UVAR_OPTIONAL, "Timescale 'tau' of transient signal window in days.");

  /* ----- 'expert-user/developer' options ----- */
  XLALregBOOLUserStruct (  hardwareTDD,         'b', UVAR_DEVELOPER, "Hardware injection: output TDD in binary format");
  XLALregSTRINGUserStruct (actuation,            0,  UVAR_DEVELOPER, "Filname containing actuation function of this detector");
  XLALregREALUserStruct (  actuationScale,       0,  UVAR_DEVELOPER,  "(Signed) scale-factor to apply to the actuation-function.");
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
  XLALDestroyTimestampVector ( cfg->timestamps );

  XLALFree ( cfg->pulsar.Doppler.orbit );

  /* free transfer-function if we have one.. */
  XLALDestroyCOMPLEX8FrequencySeries ( cfg->transfer );

  XLALFree ( cfg->VCSInfoString );

  // free noise-SFT catalog
  XLALDestroySFTCatalog ( cfg->noiseCatalog );

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

/**
 * Reads an actuation-function in format (r,phi) from file 'fname',
 * and returns the associated transfer-function as a COMPLEX8FrequencySeries (Re,Im)
 * The transfer-function T is simply the inverse of the actuation A, so T=A^-1.
 */
COMPLEX8FrequencySeries *
XLALLoadTransferFunctionFromActuation ( REAL8 actuationScale, /**< overall scale-factor to actuation */
                                        const CHAR *fname     /**< file containing actuation-function */
                                        )
{
  int len;

  XLAL_CHECK_NULL ( fname != NULL, XLAL_EINVAL, "Invalid NULL input 'fname'\n" );

  LALParsedDataFile *fileContents = NULL;
  XLAL_CHECK_NULL ( XLALParseDataFile ( &fileContents, fname ) == XLAL_SUCCESS, XLAL_EFUNC );


  /* skip first line if containing NaN's ... */
  UINT4 startline = 0;
  if ( strstr ( fileContents->lines->tokens[startline], "NaN" ) != NULL ) {
    startline ++;
  }

  UINT4 numLines = fileContents->lines->nTokens - startline;
  COMPLEX8Vector *data = XLALCreateCOMPLEX8Vector ( numLines );
  XLAL_CHECK_NULL ( data != NULL, XLAL_EFUNC, "XLALCreateCOMPLEX8Vector(%d) failed\n", numLines );

  COMPLEX8FrequencySeries *ret = XLALCalloc (1, len = sizeof(*ret) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc (1, %d)\n", len );

  snprintf ( ret->name, LALNameLength-1, "Transfer-function from: %s", fname );
  ret->name[LALNameLength-1]=0;

  /* initialize loop */
  REAL8 f0 = 0;
  REAL8 f1 = 0;

  const CHAR *readfmt = "%" LAL_REAL8_FORMAT "%" LAL_REAL8_FORMAT "%" LAL_REAL8_FORMAT;

  for ( UINT4 i = startline; i < fileContents->lines->nTokens; i++ )
    {
      CHAR *thisline = fileContents->lines->tokens[i];
      REAL8 amp, phi;

      f0 = f1;
      if ( 3 != sscanf (thisline, readfmt, &f1, &amp, &phi) ) {
        XLAL_ERROR_NULL ( XLAL_EIO, "Failed to read 3 floats from line %d of file %s\n\n", i, fname );
      }

      if ( !gsl_finite ( amp ) || !gsl_finite ( phi ) ) {
        XLAL_ERROR_NULL ( XLAL_EINVAL, "ERROR: non-finite number encountered in actuation-function at f!=0. Line=%d\n\n", i);
      }

      /* first line: set f0 */
      if ( i == startline )
	ret->f0 = f1;

      /* second line: set deltaF */
      if ( i == startline + 1 )
	ret->deltaF = f1 - ret->f0;

      /* check constancy of frequency-step */
      if ( (f1 - f0 != ret->deltaF) && (i > startline) ) {
        XLAL_ERROR_NULL ( XLAL_EINVAL, "Illegal frequency-step %f != %f in line %d of file %s\n\n", (f1-f0), ret->deltaF, i, fname );
      }

      /* now convert into transfer-function and (Re,Im): T = A^-1 */
      data->data[i-startline].re =  cos(phi) / ( amp * actuationScale );
      data->data[i-startline].im = -sin(phi) / ( amp * actuationScale );

    } /* for i < numlines */

  XLALDestroyParsedDataFile ( fileContents);

  ret->data = data;

  return ret;

} /* XLALLoadTransferFunctionFromActuation() */

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
 * Generate fake 'CW' data, returned either as SFTVector or REAL4Timeseries or both,
 * for given CW-signal ("pulsar") parameters and output parameters (frequency band etc)
 */
int
XLALMakeFakeCWData ( SFTVector **outSFTs,		//< [out] pointer to optional SFT-vector for output [FIXME! Multi-]
                     REAL4TimeSeries **outTseries,	//< [out] pointer to optional timeseries-vector for output [FIXME! Multi-]
                     const UserVariables_t *uvar,	//< [in] user-input vars [FIXME!]
                     const ConfigVars_t *cfg		//< [in] config-vars	[FIXME!]
                     )
{
  XLAL_CHECK ( (outSFTs == NULL) || ((*outSFTs) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (outTseries == NULL) || ((*outTseries) == NULL ), XLAL_EINVAL );

  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );
  XLAL_CHECK ( cfg != NULL, XLAL_EINVAL );


  // ---------- for SFT output: calculate effective fmin and Band consistent with SFT bins

  /* lowest and highest frequency boundaries have to fall on exact SFT bins, ie they must be
   * exact integer multiples of 1/Tsft:
   * ==> calculate "effective" fmin by rounding down from uvar->fmin to closest fmin_eff,
   * such that fmin_eff * Tsft = integer
   * and idem for fmax_eff, the effective Band is then fBand_eff = fmax_eff - fmin_eff
   *
   * 'exact' here means being within 10*LAL_REAL8_EPS ~2e-15 relative deviation, to
   * avoid unneccessary increases in the created frequency-Band due to numerical noise alone
   * (this is consistent with what XLALExtractBandFromSFTVector() does, for example)
   */
  volatile REAL8 dFreq = 1.0 / uvar->Tsft;
  volatile REAL8 tmp;
  // NOTE: don't "simplify" the above: we try to make sure
  // the result of this will be guaranteed to be IEEE-compliant,
  // and identical to other locations, such as in SFT-IO

  REAL8 eps = 10 * LAL_REAL8_EPS;	// about ~2e-15
  REAL8 fudge_up   = 1 + eps;
  REAL8 fudge_down = 1 - eps;

  REAL8 fMin = uvar->fmin;
  tmp = fMin / dFreq;
  UINT4 imin = (UINT4) floor( tmp * fudge_up );	// round *down*, allowing for eps 'fudge'
  REAL8 fmin_eff = imin * dFreq;

  REAL8 fMax = uvar->fmin + uvar->Band;
  tmp = fMax / dFreq;
  UINT4 imax = (UINT4) ceil ( tmp * fudge_down );  // round *up*, allowing for eps fudge
  UINT4 numBins = (UINT4) (imax - imin + 1);
  REAL8 fmax_eff = imax * dFreq;

  REAL8 fBand_eff = (numBins - 1) * dFreq;

  if ( fBand_eff != uvar->Band ) {
    XLALPrintWarning("Asked for Band [%.16g, %.16g] Hz, effective Band produced is [%.16g, %.16g] Hz (numSFTBins=%d)\n", fMin, fMax, fmin_eff, fmax_eff, numBins);
  }

  /*----------------------------------------
   * fill the PulsarSignalParams struct
   *----------------------------------------*/
  PulsarSignalParams params = empty_PulsarSignalParams;
  /* pulsar params */
  params.pulsar.refTime            = cfg->pulsar.Doppler.refTime;
  params.pulsar.position.system    = COORDINATESYSTEM_EQUATORIAL;
  params.pulsar.position.longitude = cfg->pulsar.Doppler.Alpha;
  params.pulsar.position.latitude  = cfg->pulsar.Doppler.Delta;
  {
    REAL8 h0   = cfg->pulsar.Amp.h0;
    REAL8 cosi = cfg->pulsar.Amp.cosi;
    params.pulsar.aPlus		   = 0.5 * h0 * ( 1.0 + SQ(cosi) );
    params.pulsar.aCross	   = h0 * cosi;
  }
  params.pulsar.phi0		   = cfg->pulsar.Amp.phi0;
  params.pulsar.psi 		   = cfg->pulsar.Amp.psi;

  // translate 'modern' fkdot into 'old-style' spindown-vector
  UINT4 maxSpindownOrder = 0;
  for ( UINT4 s = PULSAR_MAX_SPINS-1; s > 0; s -- ) {
    if ( cfg->pulsar.Doppler.fkdot[s] != 0 )
      {
        maxSpindownOrder = s;
        break;
      }
  } // for s = sMax ... 1
  REAL8Vector *spindown = NULL;
  XLAL_CHECK ( (spindown = XLALCreateREAL8Vector ( maxSpindownOrder )) != NULL, XLAL_EFUNC );
  for ( UINT4 s = 0; s < maxSpindownOrder; s ++ )
    {
      spindown->data[s] = cfg->pulsar.Doppler.fkdot[s+1];
    }
  params.pulsar.f0		   = cfg->pulsar.Doppler.fkdot[0];
  params.pulsar.spindown           = spindown;
  params.orbit                     = cfg->pulsar.Doppler.orbit;

  /* detector params */

  /* ---------- prepare detector ---------- */
  LALDetector *site;
  XLAL_CHECK ( XLALUserVarWasSet ( &uvar->IFO ), XLAL_EINVAL, "Need detector input --IFO!\n\n");
  XLAL_CHECK ( (site = XLALGetSiteInfo ( uvar->IFO )) != NULL, XLAL_EFUNC , "XLALGetSiteInfo('%s') failed\n", uvar->IFO );

  params.transfer = cfg->transfer;	/* detector transfer function (NULL if not used) */
  params.site = site;
  params.ephemerides = cfg->edat;

  /* characterize the output time-series */
  REAL8 min_fSamp = 2.0 * fBand_eff;	/* minimal sampling rate, via Nyquist theorem */

  // by construction, Tsft / min_fSamp = integer, but if there are gaps between SFTs,
  // then we might have to sample at higher rates in order for all SFT-boundaries to
  // lie on an exact timestep of the timeseries.
  // ie we start from min_fsamp = Tsft / n0, and then try to find the smallest
  // n >= n0, such that for fsamp = Tsft / n, for all gaps i: Dt_i * fsamp = int
  UINT4 n0 = (UINT4) round ( uvar->Tsft * min_fSamp );
  REAL8 fSamp;
  XLAL_CHECK ( XLALFindSmallestValidSamplingRate ( &fSamp, n0, cfg->timestamps ) == XLAL_SUCCESS, XLAL_EFUNC );

  params.samplingRate 	= fSamp;
  params.fHeterodyne  	= fmin_eff;		/* heterodyning frequency for output time-series */

  /* ----- figure out start-time and duration ----- */
  LIGOTimeGPS firstGPS = cfg->timestamps->data[0];
  LIGOTimeGPS  lastGPS  = cfg->timestamps->data[cfg->timestamps->length - 1 ];
  REAL8 duration = XLALGPSDiff ( &lastGPS, &firstGPS ) + uvar->Tsft;

  XLAL_CHECK ( duration >= uvar->Tsft, XLAL_EINVAL, "Requested duration=%d sec is less than Tsft =%.0f sec.\n\n", duration, uvar->Tsft);

  params.duration     = (UINT4) ceil ( duration );
  params.startTimeGPS = cfg->timestamps->data[0];

  /*----------------------------------------
   * generate the signal time-series
   *----------------------------------------*/
  REAL4TimeSeries *Tseries = NULL;
  if ( uvar->lineFeature )
    {
      XLAL_CHECK ( (Tseries = XLALGenerateLineFeature (  &params )) != NULL, XLAL_EFUNC );
    }
  else
    {
      XLAL_CHECK ( (Tseries = XLALGeneratePulsarSignal ( &params )) != NULL, XLAL_EFUNC );
    }

  // ----- free temporary memory
  XLALFree ( site );
  XLALDestroyREAL8Vector ( spindown );

  XLAL_CHECK ( XLALApplyTransientWindow ( Tseries, cfg->transientWindow ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* add Gaussian noise if requested */
  if ( uvar->noiseSqrtSh > 0)
    {
      REAL8 noiseSigma = uvar->noiseSqrtSh * sqrt ( fBand_eff );
      XLAL_CHECK ( XLALAddGaussianNoise ( Tseries, noiseSigma, uvar->randSeed ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  /*----------------------------------------
   * last step: turn this timeseries into SFTs
   *----------------------------------------*/
  SFTParams sftParams = empty_SFTParams;

  sftParams.Tsft = uvar->Tsft;

  sftParams.timestamps = cfg->timestamps;
  sftParams.noiseSFTs = NULL;	// not used here any more!

  /*--------------------- Prepare windowing of time series ---------------------*/

  REAL4Window *window = NULL;
  if ( uvar->windowType )
    {
      /* Determine SFT window */
      REAL8 dt = Tseries->deltaT;
      REAL8 REALnumTimesteps = uvar->Tsft / dt;
      UINT4 numTimesteps = round ( REALnumTimesteps );	/* number of time-samples in an Tsft (should be exact integer!) */

      XLAL_CHECK ( (window = XLALCreateNamedREAL4Window ( uvar->windowType, uvar->windowBeta, numTimesteps )) != NULL, XLAL_EFUNC );

    } // if uvar->windowType

  sftParams.window = window;

  /* get SFTs from timeseries */
  SFTVector *SFTs;
  XLAL_CHECK ( (SFTs = XLALSignalToSFTs (Tseries, &sftParams)) != NULL, XLAL_EFUNC );

  XLALDestroyREAL4Window ( window );

  // return Timeseries if requested
  if ( outTseries )
    {
      (*outTseries) = Tseries;
    }
  else
    {
      XLALDestroyREAL4TimeSeries ( Tseries );
    }
  // return SFTs if requested
  if ( outSFTs )
    {
      (*outSFTs) = SFTs;
    }
  else
    {
      XLALDestroySFTVector ( SFTs );
    }

  return XLAL_SUCCESS;

} // XLALMakeFakeCWData()


  /**
   * Find the smallest sampling rate of the form fsamp = n / Tsft, with n>=n0,
   * such that all gap sizes Dt_i between SFTs of the given timestamps are also
   * exactly resolved, ie. that Dt_i * fsamp = integer, for all i
   *
   * The smallest possible sampling rate is fsamp0 = n0 / Tsft, which is specified
   * by the user. This sampling rate would be valid if there are no gaps between SFTs,
   * so it's only in cases of gaps (which are not integer-multiples of Tsft) that we'll
   * have to increase the sampling rate.
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
XLALFindSmallestValidSamplingRate ( REAL8 *fSamp,			//< [out] minimal valid sampling rate
                                    UINT4 n0, 				//< [in] minimal sampling rate n0/Tsft
                                    const LIGOTimeGPSVector *timestamps	//< [in] start-timestamps and length Tsft of SFTs
                                    )
{
  XLAL_CHECK ( fSamp != NULL, XLAL_EINVAL );
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


  // our final minimal valid sampling rate is therefore
  (*fSamp) = nCur / TsftREAL;

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
