/*
*  Copyright (C) 2008, 2010, 2022 Karl Wette
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \file
 * \ingroup lalpulsar_bin_Tools
 * \author R. Prix, M.A. Papa, X. Siemens, B. Allen, C. Messenger
 */

/*-----------------------------------------------------------------------
 *
 * File Name: makefakedata_v4.c
 *
 * Authors: R. Prix, M.A. Papa, X. Siemens, B. Allen, C. Messenger
 *
 * This code is a descendant of an earlier implementation 'makefakedata_v2.[ch]'
 * by Badri Krishnan, Bruce Allen, Maria Alessandra Papa, Reinhard Prix, Xavier Siemens, Yousuke Itoh
 *
 *-----------------------------------------------------------------------
 */

/* ---------- includes ---------- */
#include "config.h"

#include <sys/stat.h>

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FrequencySeries.h>
#include <lal/LALInitBarycenter.h>
#include <gsl/gsl_math.h>

#include <lal/LALString.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/SimulatePulsarSignal.h>
#include <lal/TimeSeries.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/Window.h>
#include <lal/TranslateAngles.h>
#include <lal/TranslateMJD.h>
#include <lal/TransientCW_utils.h>
#include <lal/LALPulsarVCSInfo.h>

#ifdef HAVE_LIBLALFRAME
#include <lal/LALFrameIO.h>
#endif

/***************************************************/
#define SQ(x) ( (x) * (x) )

#define TRUE    1
#define FALSE   0

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

  REAL8 fmin_eff;		/**< 'effective' fmin: round down such that fmin*Tsft = integer! */
  REAL8 fBand_eff;		/**< 'effective' frequency-band such that samples/SFT is integer */
  REAL8Vector *spindown;	/**< vector of frequency-derivatives of GW signal */

  SFTVector *noiseSFTs;		/**< vector of noise-SFTs to be added to signal */
  REAL8 noiseSigma;		/**< sigma for Gaussian noise to be added */

  REAL4Window *window;		/**< window function for the time series */
  CHAR *window_type;            /**< name of window function */
  REAL8 window_param;           /**< parameter of window function */

  COMPLEX8FrequencySeries *transfer;  /**< detector's transfer function for use in hardware-injection */

  transientWindow_t transientWindow;	/**< properties of transient-signal window */
  CHAR *VCSInfoString;          /**< Git version string */
} ConfigVars_t;

typedef enum
  {
    GENERATE_ALL_AT_ONCE = 0,	/**< generate whole time-series at once before turning into SFTs */
    GENERATE_PER_SFT,		/**< generate timeseries individually for each SFT */
    GENERATE_LAST		/**< end-marker */
  } GenerationMode;

// ----- User variables
typedef struct
{
  /* output */
  CHAR *outSFTbname;		/**< Path and basefilename of output SFT files */
  BOOLEAN outSingleSFT;	/**< use to output a single concatenated SFT */

  CHAR *TDDfile;		/**< Filename for ASCII output time-series */
  CHAR *TDDframedir;		/**< directory for frame file output time-series */
  CHAR *frameDesc;           	/**< description field entry in the frame filename */

  BOOLEAN hardwareTDD;	/**< Binary output timeseries in chunks of Tsft for hardware injections. */

  CHAR *logfile;		/**< name of logfile */

  /* specify start + duration */
  CHAR *timestampsFile;	/**< Timestamps file */
  LIGOTimeGPS startTime;		/**< Start-time of requested signal in detector-frame (GPS seconds) */
  INT4 duration;		/**< Duration of requested signal in seconds */

  /* generation mode of timer-series: all-at-once or per-sft */
  INT4 generationMode;	/**< How to generate the timeseries: all-at-once or per-sft */

  /* time-series sampling + heterodyning frequencies */
  REAL8 fmin;		/**< Lowest frequency in output SFT (= heterodyning frequency) */
  REAL8 Band;		/**< bandwidth of output SFT in Hz (= 1/2 sampling frequency) */

  /* SFT params */
  REAL8 Tsft;		/**< SFT time baseline Tsft */
  REAL8 SFToverlap;		/**< overlap SFTs by this many seconds */

  /* noise to add [OPTIONAL] */
  CHAR *noiseSFTs;		/**< Glob-like pattern specifying noise-SFTs to be added to signal */
  REAL8 noiseSqrtSh;		/**< single-sided sqrt(Sh) for Gaussian noise */

  /* Window function [OPTIONAL] */
  CHAR *window;               /**< Windowing function for the time series */
  REAL8 windowParam;          /**< Hann fraction of Tukey window (0.0=rect; 1,0=han; 0.5=default */

  /* Detector and ephemeris */
  CHAR *IFO;			/**< Detector: H1, L1, H2, V1, ... */

  CHAR *actuation;		/**< filename containg detector actuation function */
  REAL8 actuationScale;	/**< Scale-factor to be applied to actuation-function */

  CHAR *ephemEarth;		/**< Earth ephemeris file to use */
  CHAR *ephemSun;		/**< Sun ephemeris file to use */

  /* pulsar parameters [REQUIRED] */
  LIGOTimeGPS refTime;		/**< Pulsar reference epoch tRef in SSB ('0' means: use startTime converted to SSB) */

  REAL8 h0;			/**< overall signal amplitude h0 */
  REAL8 cosi;		/**< cos(iota) of inclination angle iota */
  REAL8 aPlus;		/**< ALTERNATIVE to {h0,cosi}: Plus polarization amplitude aPlus */
  REAL8 aCross;		/**< ALTERNATIVE to {h0,cosi}: Cross polarization amplitude aCross */
  REAL8 psi;			/**< Polarization angle psi */
  REAL8 phi0;		/**< Initial phase phi */

  REAL8 Alpha;		/**< Right ascension [radians] alpha of pulsar */
  REAL8 Delta;		/**< Declination [radians] delta of pulsar */

  REAL8 Freq;
  REAL8 f1dot;		/**< First spindown parameter f' */
  REAL8 f2dot;		/**< Second spindown parameter f'' */
  REAL8 f3dot;		/**< Third spindown parameter f''' */

  /* orbital parameters [OPTIONAL] */

  REAL8 orbitasini;	        /**< Projected orbital semi-major axis in seconds (a/c) */
  REAL8 orbitEcc;	        /**< Orbital eccentricity */
  LIGOTimeGPS  orbitTp;		/**< 'true' epoch of periapsis passage */
  REAL8 orbitPeriod;		/**< Orbital period (seconds) */
  REAL8 orbitArgp;	        /**< Argument of periapsis (radians) */

  /* precision-level of signal-generation */
  BOOLEAN exactSignal;	/**< generate signal timeseries as exactly as possible (slow) */
  BOOLEAN lineFeature;	/**< generate a monochromatic line instead of a pulsar-signal */

  INT4 randSeed;		/**< allow user to specify random-number seed for reproducible noise-realizations */

  CHAR *parfile;             /** option .par file path */
  CHAR *transientWindowType;	/**< name of transient window ('rect', 'exp',...) */
  REAL8 transientStartTime;	/**< GPS start-time of transient window */
  REAL8 transientTauDays;	/**< time-scale in days of transient window */

  REAL8 sourceDeltaT;   /**< source-frame sampling period. '0' implies previous internal defaults */

} UserVariables_t;

// ----- global variables ----------

// ---------- local prototypes ----------
int XLALInitUserVars ( UserVariables_t *uvar, int argc, char *argv[] );
int XLALInitMakefakedata ( ConfigVars_t *cfg, UserVariables_t *uvar );
int XLALWriteMFDlog ( const char *logfile, const ConfigVars_t *cfg );
COMPLEX8FrequencySeries *XLALLoadTransferFunctionFromActuation ( REAL8 actuationScale, const CHAR *fname );
int XLALFreeMem ( ConfigVars_t *cfg );

BOOLEAN is_directory ( const CHAR *fname );
int XLALIsValidDescriptionField ( const char *desc );

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  ConfigVars_t XLAL_INIT_DECL(GV);
  PulsarSignalParams XLAL_INIT_DECL(params);
  REAL4TimeSeries *Tseries = NULL;
  UINT4 i_chunk, numchunks;
  FILE *fpSingleSFT = NULL;
  size_t len;
  UserVariables_t XLAL_INIT_DECL(uvar);


  /* ------------------------------
   * read user-input and set up shop
   *------------------------------*/
  XLAL_CHECK ( XLALInitUserVars ( &uvar, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALInitMakefakedata ( &GV, &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  /*----------------------------------------
   * fill the PulsarSignalParams struct
   *----------------------------------------*/
  /* pulsar params */
  params.pulsar.refTime            = GV.pulsar.Doppler.refTime;
  params.pulsar.position.system    = COORDINATESYSTEM_EQUATORIAL;
  params.pulsar.position.longitude = GV.pulsar.Doppler.Alpha;
  params.pulsar.position.latitude  = GV.pulsar.Doppler.Delta;
  params.pulsar.aPlus		   = GV.pulsar.Amp.aPlus;
  params.pulsar.aCross         = GV.pulsar.Amp.aCross;
  params.pulsar.phi0		   = GV.pulsar.Amp.phi0;
  params.pulsar.psi 		   = GV.pulsar.Amp.psi;

  params.pulsar.f0		   = GV.pulsar.Doppler.fkdot[0];
  params.pulsar.spindown           = GV.spindown;
  params.orbit.tp                  = GV.pulsar.Doppler.tp;
  params.orbit.argp                = GV.pulsar.Doppler.argp;
  params.orbit.asini               = GV.pulsar.Doppler.asini;
  params.orbit.ecc                 = GV.pulsar.Doppler.ecc;
  params.orbit.period              = GV.pulsar.Doppler.period;

  params.sourceDeltaT              = uvar.sourceDeltaT;

  /* detector params */
  params.transfer = GV.transfer;	/* detector transfer function (NULL if not used) */
  params.site = &(GV.site);
  params.ephemerides = GV.edat;

  /* characterize the output time-series */
  if ( ! uvar.exactSignal )	/* GeneratePulsarSignal() uses 'idealized heterodyning' */
    {
      params.samplingRate 	= 2.0 * GV.fBand_eff;	/* sampling rate of time-series (=2*frequency-Band) */
      params.fHeterodyne  	= GV.fmin_eff;		/* heterodyning frequency for output time-series */
    }
  else	/* in the exact-signal case: don't do heterodyning, sample at least twice highest frequency */
    {
      params.samplingRate 	= fmax ( 2.0 * (params.pulsar.f0 + 2 ), 2*(GV.fmin_eff + GV.fBand_eff ) );
      params.fHeterodyne 	= 0;
    }

  /* set-up main-loop according to 'generation-mode' (all-at-once' or 'per-sft') */
  switch ( uvar.generationMode )
    {
    case GENERATE_ALL_AT_ONCE:
      params.duration     = GV.duration;
      numchunks = 1;
      break;
    case GENERATE_PER_SFT:
      params.duration = (UINT4) ceil(uvar.Tsft);
      numchunks = GV.timestamps->length;
      break;
    default:
      XLAL_ERROR ( XLAL_EINVAL, "Illegal value for generationMode %d\n\n", uvar.generationMode );
      break;
    } /* switch generationMode */

  /* if user requesting single concatenated SFT */
  if ( uvar.outSingleSFT )
    {
      /* check that user isn't giving a directory */
      if ( is_directory ( uvar.outSFTbname ) ) {
        XLAL_ERROR ( XLAL_ETYPE, "'%s' is a directory, but --outSingleSFT expects a filename!\n", uvar.outSFTbname);
      }

      /* open concatenated SFT file for writing */
      if ( (fpSingleSFT = fopen ( uvar.outSFTbname, "wb" )) == NULL ) {
        XLAL_ERROR ( XLAL_EIO, "Failed to open file '%s' for writing: %s\n\n", uvar.outSFTbname, strerror(errno));
      }
    } // if outSingleSFT

  /* ----------
   * Main loop: produce time-series and turn it into SFTs,
   * either all-at-once or per-sft
   * ----------*/
  for ( i_chunk = 0; i_chunk < numchunks; i_chunk++ )
    {
      params.startTimeGPS = GV.timestamps->data[i_chunk];

      /*----------------------------------------
       * generate the signal time-series
       *----------------------------------------*/
      if ( uvar.exactSignal )
	{
          XLAL_CHECK ( (Tseries = XLALSimulateExactPulsarSignal ( &params )) != NULL, XLAL_EFUNC );
	}
      else if ( uvar.lineFeature )
	{
          XLAL_CHECK ( (Tseries = XLALGenerateLineFeature (  &params )) != NULL, XLAL_EFUNC );
	}
      else
	{
	  XLAL_CHECK ( (Tseries = XLALGeneratePulsarSignal ( &params )) != NULL, XLAL_EFUNC );
	}

      XLAL_CHECK ( XLALApplyTransientWindow ( Tseries, GV.transientWindow ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* for HARDWARE-INJECTION:
       * before the first chunk we send magic number and chunk-length to stdout
       */
      if ( uvar.hardwareTDD && (i_chunk == 0) )
	{
	  REAL4 magic = 1234.5;
	  UINT4 length = Tseries->data->length;
	  if ( (1 != fwrite ( &magic, sizeof(magic), 1, stdout )) || (1 != fwrite(&length, sizeof(INT4), 1, stdout)) )
	    {
	      perror ("Failed to write to stdout");
	      XLAL_ERROR ( XLAL_EIO );
	    }
	} /* if hardware-injection and doing first chunk */


      /* add Gaussian noise if requested */
      if ( GV.noiseSigma > 0) {
        // NOTE: seed=0 means randomize seed from /dev/urandom, otherwise we'll have to increment it for each chunk here
	XLAL_CHECK ( XLALAddGaussianNoise ( Tseries, GV.noiseSigma, (uvar.randSeed == 0) ? 0 : (uvar.randSeed + i_chunk) ) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      /* output ASCII time-series if requested */
      if ( uvar.TDDfile )
	{
	  CHAR *fname = XLALCalloc (1, len = strlen(uvar.TDDfile) + 10 );
          XLAL_CHECK ( fname != NULL, XLAL_ENOMEM, "XLALCalloc(1,%zu) failed\n", len );
	  sprintf (fname, "%s.%02d", uvar.TDDfile, i_chunk);
	  XLAL_CHECK ( XLALdumpREAL4TimeSeries ( fname, Tseries ) == XLAL_SUCCESS, XLAL_EFUNC );
	  XLALFree (fname);
	} /* if outputting ASCII time-series */

      /* output time-series to frames if requested */
      if ( uvar.TDDframedir )
	{
#ifndef HAVE_LIBLALFRAME
          XLAL_ERROR ( XLAL_EINVAL, "--TDDframedir option not supported, code has to be compiled with lalframe\n" );
#else
	  /* use standard frame output filename format */
          XLAL_CHECK ( XLALCheckValidDescriptionField ( uvar.frameDesc ) == XLAL_SUCCESS, XLAL_EFUNC );
          len = strlen(uvar.TDDframedir) + strlen(uvar.frameDesc) + 100;
	  char *fname;
          char IFO[2] = { Tseries->name[0], Tseries->name[1] };
          XLAL_CHECK ( (fname = LALCalloc (1, len )) != NULL, XLAL_ENOMEM );
          size_t written = snprintf ( fname, len, "%s/%c-%c%c_%s-%d-%d.gwf",
                                      uvar.TDDframedir, IFO[0], IFO[0], IFO[1], uvar.frameDesc, params.startTimeGPS.gpsSeconds, (int)params.duration );
          XLAL_CHECK ( written < len, XLAL_ESIZE, "Frame-filename exceeds expected maximal length (%zu): '%s'\n", len, fname );

	  /* define the output frame */
	  LALFrameH *outFrame;
	  XLAL_CHECK ( (outFrame = XLALFrameNew( &(params.startTimeGPS), params.duration, uvar.frameDesc, 1, 0, 0 )) != NULL, XLAL_EFUNC );

	  /* add timeseries to the frame - make sure to change the timeseries name since this is used as the channel name */
          char buffer[LALNameLength];
	  written = snprintf ( buffer, LALNameLength, "%s:%s", Tseries->name, uvar.frameDesc );
          XLAL_CHECK ( written < LALNameLength, XLAL_ESIZE, "Updated frame name exceeds max length (%d): '%s'\n", LALNameLength, buffer );
          strcpy ( Tseries->name, buffer );

	  XLAL_CHECK ( (XLALFrameAddREAL4TimeSeriesProcData ( outFrame, Tseries ) == XLAL_SUCCESS ) , XLAL_EFUNC );

	  /* Here's where we add extra information into the frame - first we add the command line args used to generate it */
	  char *hist = XLALUserVarGetLog (UVAR_LOGFMT_CMDLINE);
          XLALFrameAddFrHistory ( outFrame, __FILE__, hist );

	  /* then we add the version string */
	  XLALFrameAddFrHistory ( outFrame, __FILE__, GV.VCSInfoString );

	  /* output the frame to file - compression level 1 (higher values make no difference) */
	  XLAL_CHECK ( (XLALFrameWrite(outFrame, fname) == 0) , XLAL_EFUNC );

	  /* free the frame, frame file name and history memory */
	  XLALFrameFree ( outFrame );
	  LALFree ( fname );
          LALFree ( hist );
#endif
	} /* if outputting time-series to frames */


      /* if hardware injection: send timeseries in binary-format to stdout */
      if ( uvar.hardwareTDD )
	{
	  UINT4  length = Tseries->data->length;
	  REAL4 *datap = Tseries->data->data;

	  if ( length != fwrite (datap, sizeof(datap[0]), length, stdout) )
	    {
	      perror( "Fatal error in writing binary data to stdout\n");
	      XLAL_ERROR ( XLAL_EIO, "fwrite() failed\n");
	    }
	  fflush (stdout);

	} /* if hardware injections */

      /*----------------------------------------
       * last step: turn this timeseries into SFTs
       * and output them to disk
       *----------------------------------------*/
      SFTVector *SFTs = NULL;
      if (uvar.outSFTbname)
	{
	  SFTParams XLAL_INIT_DECL(sftParams);
	  LIGOTimeGPSVector ts;
	  SFTVector noise;

	  sftParams.Tsft = uvar.Tsft;

	  switch (uvar.generationMode)
	    {
	    case GENERATE_ALL_AT_ONCE:
	      sftParams.timestamps = GV.timestamps;
	      sftParams.noiseSFTs = GV.noiseSFTs;
	      break;
	    case GENERATE_PER_SFT:
	      ts.length = 1;
	      ts.data = &(GV.timestamps->data[i_chunk]);
	      sftParams.timestamps = &(ts);
	      sftParams.noiseSFTs = NULL;

	      if ( GV.noiseSFTs )
		{
		  noise.length = 1;
		  noise.data = &(GV.noiseSFTs->data[i_chunk]);
		  sftParams.noiseSFTs = &(noise);
		}

	      break;

	    default:
	      XLAL_ERROR ( XLAL_EINVAL, "Invalid Value --generationMode=%d\n", uvar.generationMode );
	      break;
	    }

	  /* Enter the window function into the SFTparams struct */
	  sftParams.window = GV.window;

	  /* get SFTs from timeseries */
          XLAL_CHECK ( (SFTs = XLALSignalToSFTs (Tseries, &sftParams)) != NULL, XLAL_EFUNC );

          /* extract requested band */
          {
            SFTVector *outSFTs;
            XLAL_CHECK ( (outSFTs = XLALExtractStrictBandFromSFTVector ( SFTs, uvar.fmin, uvar.Band )) != NULL, XLAL_EFUNC );
            XLALDestroySFTVector ( SFTs );
            SFTs = outSFTs;
          }

          /* generate comment string */
          CHAR *logstr;
          XLAL_CHECK ( (logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
          char *comment = XLALCalloc ( 1, len = strlen ( logstr ) + strlen(GV.VCSInfoString) + 512 );
          XLAL_CHECK ( comment != NULL, XLAL_ENOMEM, "XLALCalloc(1,%zu) failed.\n", len );
          sprintf ( comment, "Generated by:\n%s\n%s\n", logstr, GV.VCSInfoString );

          /* if user requesting single concatenated SFT */
          SFTFilenameSpec XLAL_INIT_DECL(spec);
          XLAL_CHECK( XLALFillSFTFilenameSpecStrings( &spec, uvar.outSFTbname, NULL, NULL, GV.window_type, "mfdv4", NULL, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
          spec.window_param = GV.window_param;
          if ( uvar.outSingleSFT )
            {
              /* write all SFTs to concatenated file */
              for ( UINT4 k = 0; k < SFTs->length; k++ )
                {
                  int ret = XLALWriteSFT2FilePointer ( &(SFTs->data[k]), fpSingleSFT, spec.window_type, spec.window_param, comment );
                  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALWriteSFT2FilePointer() failed to write SFT %d / %d to '%s'!\n", k+1, SFTs->length, uvar.outSFTbname );
                } // for k < numSFTs
            } // if outSingleSFT
          else
            {
              XLAL_CHECK( is_directory ( uvar.outSFTbname ), XLAL_EINVAL, "ERROR: the --outSFTbname '%s' is not a directory!\n", uvar.outSFTbname );
              XLAL_CHECK ( XLALWriteSFTVector2StandardFile( SFTs, &spec, comment, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );
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

    } /* for i_chunk < numchunks */

  /* if user requesting single concatenated SFT */
  if ( uvar.outSingleSFT ) {

    /* close concatenated SFT */
    fclose( fpSingleSFT );

  }

  XLALFreeMem ( &GV );	/* free the config-struct */

  LALCheckMemoryLeaks();

  return 0;
} /* main */

/**
 * Handle user-input and set up shop accordingly, and do all
 * consistency-checks on user-input.
 */
int
XLALInitMakefakedata ( ConfigVars_t *cfg, UserVariables_t *uvar )
{
  XLAL_CHECK ( cfg != NULL, XLAL_EINVAL, "Invalid NULL input 'cfg'\n" );
  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL, "Invalid NULL input 'uvar'\n");

  cfg->VCSInfoString = XLALVCSInfoString(lalPulsarVCSInfoList, 0, "%% ");
  XLAL_CHECK ( cfg->VCSInfoString != NULL, XLAL_EFUNC, "XLALVCSInfoString failed.\n" );

  BOOLEAN have_parfile = XLALUserVarWasSet (&uvar->parfile);
  BinaryPulsarParams pulparams;

  char* window_type_from_noiseSFTs = NULL;
  REAL8 window_param_from_noiseSFTs = 0;

  /* read in par file parameters if given */
   if (have_parfile)
     {
       XLALReadTEMPOParFileOrig( &pulparams, uvar->parfile);
       XLAL_CHECK ( xlalErrno == XLAL_SUCCESS, XLAL_EFUNC, "XLALReadTEMPOParFileOrig() failed for parfile = '%s', xlalErrno = %d\n", uvar->parfile, xlalErrno );
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
    BOOLEAN have_Alpha  = XLALUserVarWasSet ( &uvar->Alpha );
    BOOLEAN have_Delta  = XLALUserVarWasSet ( &uvar->Delta );

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
      {   /* translate {h_0, cosi} into A_{+,x} */
          /* assume at 2f */
          REAL8 h0 = uvar->h0;
          REAL8 cosi = uvar->cosi;
          cfg->pulsar.Amp.aPlus = 0.5 * h0 * (1.0 + SQ(cosi));
          cfg->pulsar.Amp.aCross = h0 * cosi;
      }
    else if ( have_aPlus && have_aCross )
      {
          cfg->pulsar.Amp.aPlus = uvar->aPlus;
          cfg->pulsar.Amp.aCross = uvar->aCross;
      }
    else {
        cfg->pulsar.Amp.aPlus = 0.0;
        cfg->pulsar.Amp.aCross = 0.0;
    }
    cfg->pulsar.Amp.phi0 = uvar->phi0;
    cfg->pulsar.Amp.psi  = uvar->psi;

    /* ----- signal Frequency ----- */
    cfg->pulsar.Doppler.fkdot[0] = uvar->Freq;

    /* ----- skypos ----- */
    if ( (have_Alpha && !have_Delta) || ( !have_Alpha && have_Delta ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "\nSpecify skyposition: need BOTH --Alpha and --Delta!\n\n");
    }
    cfg->pulsar.Doppler.Alpha = uvar->Alpha;
    cfg->pulsar.Doppler.Delta = uvar->Delta;

  } /* Pulsar signal parameters */

  /* ---------- prepare vector of spindown parameters ---------- */
  {
    UINT4 msp = 0;	/* number of spindown-parameters */
    if ( have_parfile )
      {
	uvar->f1dot = 2.*pulparams.f1;
	uvar->f2dot = 2.*pulparams.f2;
	uvar->f3dot = 2.*pulparams.f3;
      }
    cfg->pulsar.Doppler.fkdot[1] = uvar->f1dot;
    cfg->pulsar.Doppler.fkdot[2] = uvar->f2dot;
    cfg->pulsar.Doppler.fkdot[3] = uvar->f3dot;

    if (uvar->f3dot != 0) 	msp = 3;	/* counter number of spindowns */
    else if (uvar->f2dot != 0)	msp = 2;
    else if (uvar->f1dot != 0)	msp = 1;
    else 			msp = 0;
    if (msp)
      {
        /* memory not initialized, but ALL alloc'ed entries will be set below! */
        cfg->spindown = XLALCreateREAL8Vector ( msp );
        XLAL_CHECK ( cfg->spindown != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %d ) failed.\n", msp );
      }
    switch (msp)
      {
      case 3:
	cfg->spindown->data[2] = uvar->f3dot;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
	__attribute__ ((fallthrough));
#endif
      case 2:
	cfg->spindown->data[1] = uvar->f2dot;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
	__attribute__ ((fallthrough));
#endif
      case 1:
	cfg->spindown->data[0] = uvar->f1dot;
	break;
      case 0:
	break;
      default:
	XLAL_ERROR ( XLAL_EINVAL, "\nmsp = %d makes no sense to me...\n\n", msp);
      } /* switch(msp) */

  } /* END: prepare spindown parameters */

  CHAR *channelName = NULL;
  /* ---------- prepare detector ---------- */
  {
    const LALDetector *site;

    site = XLALGetSiteInfo ( uvar->IFO );
    XLAL_CHECK ( site != NULL, XLAL_EFUNC, "XLALGetSiteInfo('%s') failed\n", uvar->IFO );
    channelName = XLALGetChannelPrefix ( uvar->IFO );
    XLAL_CHECK ( channelName != NULL, XLAL_EFUNC, "XLALGetChannelPrefix('%s') failed.\n", uvar->IFO );

    cfg->site = (*site);
  }

   /* check for negative fmin and Band, which would break the fmin_eff, fBand_eff calculation below */
  XLAL_CHECK ( uvar->fmin >= 0, XLAL_EDOM, "Invalid negative frequency fmin=%f!\n\n", uvar->fmin );
  XLAL_CHECK ( uvar->Band >= 0, XLAL_EDOM, "Invalid negative frequency band Band=%f!\n\n", uvar->Band );

  /* ---------- for SFT output: calculate effective fmin and Band ---------- */
  // Note: this band is only used for internal data operations; ultimately SFTs covering
  // the half-open interval uvar->[fmin,fmin+Band) are returned to the user using
  // XLALExtractStrictBandFromSFTVector()
  if ( XLALUserVarWasSet( &uvar->outSFTbname ) )
    {
      UINT4 firstBin, numBins;
      /* calculate "effective" fmin from uvar->fmin:
       * make sure that fmin_eff * Tsft = integer, such that freqBinIndex corresponds
       * to a frequency-index of the non-heterodyned signal.
       */

      XLAL_CHECK ( XLALFindCoveringSFTBins ( &firstBin, &numBins, uvar->fmin, uvar->Band, uvar->Tsft ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* Adjust Band correspondingly */
      REAL8 dFreq = 1.0 / uvar->Tsft;
      cfg->fmin_eff  = firstBin * dFreq;
      cfg->fBand_eff = (numBins-1) * dFreq;

      if ( lalDebugLevel )
	{
	  if ( fabs(cfg->fmin_eff - uvar->fmin)> LAL_REAL8_EPS
	       || fabs(cfg->fBand_eff - uvar->Band) > LAL_REAL8_EPS )
	    printf("\nWARNING: for SFT-creation we had to adjust (fmin,Band) to"
		   " fmin_eff=%.15g and Band_eff=%.15g\n\n", cfg->fmin_eff, cfg->fBand_eff);
	}

    } /* END: SFT-specific corrections to fmin and Band */
  else
    { /* producing pure time-series output ==> no adjustments necessary */
      cfg->fmin_eff = uvar->fmin;
      cfg->fBand_eff = uvar->Band;
    }


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

    /* don't allow --SFToverlap with generationMode==1 (one timeseries per SFT), as this would
     * result in independent noise realizations in the overlapping regions
     */
    if ( haveOverlap && (uvar->generationMode != 0) ) {
      XLAL_ERROR ( XLAL_EINVAL, "--SFToverlap can only be used with --generationMode=0, otherwise we'll get overlapping independent noise!\n\n");
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

      } /* if haveTimestampsFile */

    /* ----- if real noise-SFTs given: load them now using EITHER (start,start+duration) OR timestamps
     * as constraints if given, otherwise load all of them. Also require window option to be given
     */
    if ( uvar->noiseSFTs )
      {
	REAL8 fMin, fMax;
	SFTConstraints XLAL_INIT_DECL(constraints);
	LIGOTimeGPS minStartTime, maxStartTime;

        XLALPrintWarning ( "\nWARNING: only SFTs corresponding to the noiseSFT-timestamps will be produced!\n" );

	/* use all additional constraints user might have given */
	if ( haveStart && haveDuration )
	  {
	    minStartTime = maxStartTime = uvar->startTime;
            XLALGPSAdd ( &maxStartTime, uvar->duration );
            constraints.minStartTime = &minStartTime;
	    constraints.maxStartTime = &maxStartTime;
            char bufGPS1[32], bufGPS2[32];
            XLALPrintWarning ( "\nWARNING: only noise-SFTs between GPS [%s, %s] will be used!\n", XLALGPSToStr (bufGPS1, &minStartTime), XLALGPSToStr (bufGPS2, &maxStartTime) );
	  } /* if start+duration given */
	if ( cfg->timestamps )
	  {
	    constraints.timestamps = cfg->timestamps;
            XLALPrintWarning ( "\nWARNING: only noise-SFTs corresponding to given timestamps '%s' will be used!\n", uvar->timestampsFile );
	  } /* if we have timestamps already */

	/* use detector-constraint */
	constraints.detector = channelName ;

	SFTCatalog *catalog = NULL;
	XLAL_CHECK ( (catalog = XLALSFTdataFind ( uvar->noiseSFTs, &constraints )) != NULL, XLAL_EFUNC );

	/* check if anything matched */
	if ( catalog->length == 0 ) {
          XLAL_ERROR ( XLAL_EFAILED, "No noise-SFTs matching the constraints (IFO, start+duration, timestamps) were found!\n" );
        }

        /* extract SFT window function */
        window_type_from_noiseSFTs = XLALStringDuplicate(catalog->data[0].window_type);
        XLAL_CHECK ( window_type_from_noiseSFTs != NULL, XLAL_ENOMEM );
        window_param_from_noiseSFTs = catalog->data[0].window_param;

	/* load effective frequency-band from noise-SFTs */
	fMin = cfg->fmin_eff;
	fMax = fMin + cfg->fBand_eff;

        cfg->noiseSFTs = XLALLoadSFTs ( catalog, fMin, fMax );
        XLAL_CHECK ( cfg->noiseSFTs != NULL, XLAL_EFUNC, "XLALLoadSFTs() failed\n" );

        XLALDestroySFTCatalog ( catalog );

	/* get timestamps from the loaded noise SFTs */
        if ( cfg->timestamps ) {
          XLALDestroyTimestampVector ( cfg->timestamps );
        }
        cfg->timestamps = XLALExtractTimestampsFromSFTs ( cfg->noiseSFTs );
	XLAL_CHECK ( cfg->timestamps != NULL, XLAL_EFUNC, "XLALExtractTimestampsFromSFTs() failed\n" );

      } /* if uvar->noiseSFTs */

    /* have we got our timestamps yet?: If not, we must get them from (start, duration) user-input */
    if ( ! cfg->timestamps )
      {
	if ( !haveStart || !haveDuration ) {
          XLAL_ERROR ( XLAL_EINVAL, "Need to have either --timestampsFile OR (--startTime,--duration) OR --noiseSFTs\n\n");
        }
	if ( uvar->SFToverlap > uvar->Tsft ) {
          XLAL_ERROR ( XLAL_EINVAL, "--SFToverlap cannot be larger than --Tsft!\n\n");
        }

	/* internally always use timestamps, so we generate them  */
        XLAL_CHECK ( ( cfg->timestamps = XLALMakeTimestamps ( uvar->startTime, uvar->duration, uvar->Tsft, uvar->SFToverlap )) != NULL, XLAL_EFUNC );

      } /* if !cfg->timestamps */

    /* ----- figure out start-time and duration ----- */
    {
      LIGOTimeGPS t1, t0;
      REAL8 duration;

      t0 = cfg->timestamps->data[0];
      t1 = cfg->timestamps->data[cfg->timestamps->length - 1 ];

      duration = XLALGPSDiff(&t1, &t0);
      duration += uvar->Tsft;

      cfg->startTimeGPS = cfg->timestamps->data[0];
      cfg->duration = (UINT4)ceil ( duration );	/* round up to seconds */
    }

    if ( cfg->duration < uvar->Tsft ) {
      XLAL_ERROR ( XLAL_EINVAL, "Requested duration of %d sec is less than minimal chunk-size of Tsft =%.0f sec.\n\n", uvar->duration, uvar->Tsft);
    }

  } /* END: setup signal start + duration */

  /*----------------------------------------------------------------------*/
  /* currently there are only two modes: [uvar->generationMode]
   * 1) produce the whole timeseries at once [default],
   *       [GENERATE_ALL_AT_ONCE]
   * 2) produce timeseries for each SFT,
   *       [GENERATE_PER_SFT]
   *
   * Mode 2 which is useful if 1) is too memory-intensive (e.g. for hardware-injections)
   *
   * intermediate modes would require a bit more work because the
   * case with specified timestamps (allowing for gaps) would be
   * a bit more tricky ==> this is left for future extensions if found useful
   */
  if ( (uvar->generationMode < 0) || (uvar->generationMode >= GENERATE_LAST) ) {
    XLAL_ERROR ( XLAL_EINVAL, "Illegal input for 'generationMode': must lie within [0, %d]\n\n", GENERATE_LAST -1);
  }

  /* this is a violation of the UserInput-paradigm that no user-variables
   * should by modified by the code, but this is the easiest way to do
   * this here, and the log-output of user-variables will not be changed
   * by this [it's already been written], so in this case it should be safe..
   */
  if ( uvar->hardwareTDD )
    uvar->generationMode = GENERATE_PER_SFT;

  /*--------------------- Prepare windowing of time series ---------------------*/
  cfg->window = NULL;
  cfg->window_type = NULL;
  cfg->window_param = 0;

  {
    const BOOLEAN have_window = XLALUserVarWasSet ( &uvar->window );
    const BOOLEAN have_param = XLALUserVarWasSet( &uvar->windowParam );
    const CHAR* window_type_from_uvar = uvar->window;
    const REAL8 window_param_from_uvar = have_param ? uvar->windowParam : 0;
    if ( have_window )
      {
        XLAL_CHECK ( XLALCheckNamedWindow ( window_type_from_uvar, have_param ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    if ( uvar->noiseSFTs )
      {
        const BOOLEAN have_window_type_from_noiseSFTs = ( XLALStringCaseCompare( window_type_from_noiseSFTs, "unknown" ) != 0 );
        if ( have_window_type_from_noiseSFTs ^ have_window )
          {
            cfg->window_type = XLALStringDuplicate( have_window_type_from_noiseSFTs ? window_type_from_noiseSFTs : window_type_from_uvar );
            XLAL_CHECK ( cfg->window_type != NULL, XLAL_ENOMEM );
            cfg->window_param = have_window_type_from_noiseSFTs ? window_param_from_noiseSFTs : window_param_from_uvar;
          }
        else
          {
            /* EITHER noise SFTs must have a known window OR user must specify the window function */
            XLAL_ERROR ( XLAL_EINVAL, "When --noiseSFTs is given, --window is required ONLY if noise SFTs have unknown window.\n" );
          }
      }
    else if ( have_window )
      {
        cfg->window_type = XLALStringDuplicate( window_type_from_uvar );
        XLAL_CHECK ( cfg->window_type != NULL, XLAL_ENOMEM );
        cfg->window_param = window_param_from_uvar;
      }
    if ( cfg->window_type )
      {
        /* NOTE: a timeseries of length N*dT has no timestep at N*dT !! (convention) */
        UINT4 lengthOfTimeSeries = (UINT4)round(uvar->Tsft * 2 * cfg->fBand_eff);

        cfg->window = XLALCreateNamedREAL4Window( cfg->window_type, cfg->window_param, lengthOfTimeSeries );
        XLAL_CHECK ( cfg->window != NULL, XLAL_EFUNC, "XLALCreateNamedREAL4Window('%s', %g, %d) failed\n", cfg->window_type, cfg->window_param, lengthOfTimeSeries );
      }
  }

  /* Init ephemerides */
  XLAL_CHECK ( (cfg->edat = XLALInitBarycenter ( uvar->ephemEarth, uvar->ephemSun )) != NULL, XLAL_EFUNC );

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
	uvar->orbitTp.gpsSeconds = (UINT4)floor(pulparams.T0);
	uvar->orbitTp.gpsNanoSeconds = (UINT4)floor((pulparams.T0 - uvar->orbitTp.gpsSeconds)*1e9);
      }
    }
    BOOLEAN set1 = XLALUserVarWasSet(&uvar->orbitasini);
    BOOLEAN set2 = XLALUserVarWasSet(&uvar->orbitEcc);
    BOOLEAN set3 = XLALUserVarWasSet(&uvar->orbitPeriod);
    BOOLEAN set4 = XLALUserVarWasSet(&uvar->orbitArgp);
    BOOLEAN set5 = XLALUserVarWasSet(&uvar->orbitTp);

    if (set1 || set2 || set3 || set4 || set5 )
    {
      if ( (uvar->orbitasini > 0) && !(set1 && set2 && set3 && set4 && set5 ) ) {
        XLAL_ERROR ( XLAL_EINVAL, "\nPlease either specify  ALL orbital parameters or NONE!\n\n");
      }
      if ( (uvar->orbitEcc < 0) || (uvar->orbitEcc > 1) ) {
        XLAL_ERROR ( XLAL_EINVAL, "\nEccentricity = %g has to lie within [0, 1]\n\n", uvar->orbitEcc );
      }

      cfg->pulsar.Doppler.tp = uvar->orbitTp;

      /* fill in orbital parameter structure */
      cfg->pulsar.Doppler.period = uvar->orbitPeriod;
      cfg->pulsar.Doppler.asini = uvar->orbitasini;
      cfg->pulsar.Doppler.argp = uvar->orbitArgp;
      cfg->pulsar.Doppler.ecc = uvar->orbitEcc;

    } /* if one or more orbital parameters were set */
    else
      cfg->pulsar.Doppler.asini = 0 /* isolated pulsar */;
  } /* END: binary orbital params */


  /* -------------------- handle NOISE params -------------------- */
  cfg->noiseSigma = uvar->noiseSqrtSh * sqrt ( cfg->fBand_eff );	/* convert Sh -> sigma */


  /* ----- set "pulsar reference time", i.e. SSB-time at which pulsar params are defined ---------- */
  if (XLALUserVarWasSet (&uvar->parfile)) {
    XLALGPSSetREAL8( &(uvar->refTime), pulparams.pepoch ); /*XLALReadTEMPOParFileOrig converted pepoch to REAL8 */
    XLALGPSSetREAL8( &(cfg->pulsar.Doppler.refTime), pulparams.pepoch);
  }
  else if (XLALUserVarWasSet(&uvar->refTime))
    {
      cfg->pulsar.Doppler.refTime = uvar->refTime;
    }
  else
    {
      cfg->pulsar.Doppler.refTime = cfg->timestamps->data[0];	/* internal startTime always found in here*/
    }


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

  XLALFree ( channelName );

  /* ----- handle transient-signal window if given ----- */
  int twtype;
  XLAL_CHECK ( (twtype = XLALParseTransientWindowName ( uvar->transientWindowType )) >= 0, XLAL_EFUNC );
  cfg->transientWindow.type = twtype;

  cfg->transientWindow.t0   = uvar->transientStartTime;
  cfg->transientWindow.tau  = uvar->transientTauDays * 86400;

  /* free memory */
  XLALFree ( window_type_from_noiseSFTs );

  return XLAL_SUCCESS;

} /* XLALInitMakefakedata() */


/**
 * Register all "user-variables", and initialized them from command-line and config-files
 */
int
XLALInitUserVars ( UserVariables_t *uvar, int argc, char *argv[] )
{
  int len;

  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL, "Invalid NULL input 'uvar'\n");
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL, "Invalid NULL input 'argv'\n");

  // ---------- set a few defaults ----------
  uvar->ephemEarth = XLALStringDuplicate("earth00-40-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-40-DE405.dat.gz");

  uvar->Tsft = 1800;

  // per default we now generate a timeseries per SFT: slower, but avoids potential confusion about sft-"nudging"
  uvar->generationMode = GENERATE_PER_SFT;

  uvar->actuationScale = + 1.0;

#define DEFAULT_TRANSIENT "none"
  uvar->transientWindowType = XLALCalloc ( 1, len = strlen(DEFAULT_TRANSIENT)+1 );
  XLAL_CHECK ( uvar->transientWindowType != NULL, XLAL_ENOMEM, "XLALCalloc ( 1, %d ) failed.\n", len );
  strcpy ( uvar->transientWindowType, DEFAULT_TRANSIENT );

  // ---------- register all our user-variable ----------
  /* output options */
  XLALRegisterUvarMember(   outSingleSFT,       BOOLEAN, 's', OPTIONAL, "Write a single concatenated SFT (name given by --outSFTbname)" );
  XLALRegisterUvarMember( outSFTbname,        STRING, 'n', OPTIONAL, "Output SFTs: target Directory (if --outSingleSFT=false) or filename (if --outSingleSFT=true)");

  XLALRegisterUvarMember( TDDfile,		STRING, 't', OPTIONAL, "Basename for output of ASCII time-series");
  XLALRegisterUvarMember( TDDframedir,		STRING, 'F', OPTIONAL, "Directory to output frame time-series into");
  XLALRegisterUvarMember( frameDesc,	 	 STRING, 0,  OPTIONAL,  "Description-field entry in frame filename");

  XLALRegisterUvarMember( logfile,            STRING, 'l', OPTIONAL, "Filename for log-output");

  /* detector and ephemeris */
  XLALRegisterUvarMember( IFO,                STRING, 'I', REQUIRED, "Detector: one of 'G1','L1','H1,'H2','V1', ...");

  XLALRegisterUvarMember( ephemEarth, 	 	STRING, 0,  OPTIONAL, "Earth ephemeris file to use");
  XLALRegisterUvarMember( ephemSun, 	 	STRING, 0,  OPTIONAL, "Sun ephemeris file to use");

  /* start + duration of timeseries */
  XLALRegisterUvarMember( startTime,           EPOCH, 'G', OPTIONAL, "Start-time of requested signal in detector-frame (format 'xx.yy[GPS|MJD]')");
  XLALRegisterUvarMember(  duration,              INT4, 0,  OPTIONAL, "Duration of requested signal in seconds");
  XLALRegisterUvarMember( timestampsFile,      STRING, 0,  OPTIONAL, "ALTERNATIVE: File to read timestamps from (file-format: lines with <seconds> <nanoseconds>)");

  /* sampling and heterodyning frequencies */
  XLALRegisterUvarMember(  fmin,                 REAL8, 0, REQUIRED, "Lowest frequency in output SFT (= heterodyning frequency)");
  XLALRegisterUvarMember(  Band,                 REAL8, 0, REQUIRED, "Bandwidth of output SFT in Hz (= 1/2 sampling frequency)");

  /* SFT properties */
  XLALRegisterUvarMember(  Tsft,                 REAL8, 0, OPTIONAL, "Time baseline of one SFT in seconds");
  XLALRegisterUvarMember(  SFToverlap,           REAL8, 0, OPTIONAL, "Overlap between successive SFTs in seconds (conflicts with --noiseSFTs or --timestampsFile)");
  XLALRegisterUvarMember( window,               STRING, 0, OPTIONAL, "Window function to apply to the SFTs ('rectangular', 'hann', 'tukey', etc.); when --noiseSFTs is given, required ONLY if noise SFTs have unknown window");
  XLALRegisterUvarMember( windowParam,           REAL8, 0, OPTIONAL, "Window parameter required for a few window-types (eg. 'tukey')");
  XLALRegisterNamedUvar( NULL, "tukeyBeta",      REAL8, 0, DEFUNCT,  "Use " UVAR_STR( windowParam ) " instead");

  /* pulsar params */
  XLALRegisterUvarMember( refTime,             EPOCH, 'S', OPTIONAL, "Pulsar SSB reference epoch: format 'xx.yy[GPS|MJD]' [default: startTime]");

  XLALRegisterUvarMember(  Alpha,		 RAJ, 0, OPTIONAL, "Sky: equatorial J2000 right ascension (in radians or hours:minutes:seconds)");
  XLALRegisterUvarMember(  Delta,		 DECJ, 0, OPTIONAL, "Sky: equatorial J2000 declination (in radians or degrees:minutes:seconds)");

  XLALRegisterUvarMember(  h0,                   REAL8, 0, OPTIONAL, "Overall signal-amplitude h0 (for emission at twice spin frequency only)");
  XLALRegisterUvarMember(  cosi,                 REAL8, 0, OPTIONAL, "cos(iota) of inclination-angle iota (for emission at twice spin frequency only)");
  XLALRegisterUvarMember(  aPlus,                REAL8, 0, OPTIONAL, "ALTERNATIVE to {--h0,--cosi}: A_+ amplitude");
  XLALRegisterUvarMember(  aCross,               REAL8, 0, OPTIONAL, "ALTERNATIVE to {--h0,--cosi}: A_x amplitude");

  XLALRegisterUvarMember(  psi,                  REAL8, 0, OPTIONAL, "Polarization angle psi");
  XLALRegisterUvarMember(  phi0,                 REAL8, 0, OPTIONAL, "Initial GW phase phi");
  XLALRegisterUvarMember(  Freq,                 REAL8, 0, OPTIONAL, "Intrinsic GW-frequency at refTime");

  XLALRegisterUvarMember(  f1dot,                REAL8, 0, OPTIONAL, "First spindown parameter f' at refTime");
  XLALRegisterUvarMember(  f2dot,                REAL8, 0, OPTIONAL, "Second spindown parameter f'' at refTime");
  XLALRegisterUvarMember(  f3dot,                REAL8, 0, OPTIONAL, "Third spindown parameter f''' at refTime");

  /* binary-system orbital parameters */
  XLALRegisterUvarMember(  orbitasini,           REAL8, 0, OPTIONAL, "Projected orbital semi-major axis in seconds (a/c)");
  XLALRegisterUvarMember(  orbitEcc,             REAL8, 0, OPTIONAL, "Orbital eccentricity");
  XLALRegisterUvarMember( orbitTp,        	 EPOCH, 0, OPTIONAL, "True epoch of periapsis passage: format 'xx.yy[GPS|MJD]'");
  XLALRegisterUvarMember(  orbitPeriod,          REAL8, 0, OPTIONAL, "Orbital period (seconds)");
  XLALRegisterUvarMember(  orbitArgp,            REAL8, 0, OPTIONAL, "Argument of periapsis (radians)");

  /* noise */
  XLALRegisterUvarMember( noiseSFTs,          STRING, 'D', OPTIONAL, "Noise-SFTs to be added to signal. Possibilities are:\n"
                          " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'\n"
                          "(Uses ONLY SFTs falling within (--startTime,--duration) or the set given in --timstampsFile, and ONLY within (--fmin,--Band).)");
  XLALRegisterUvarMember(  noiseSqrtSh,          REAL8, 0, OPTIONAL,  "Gaussian noise with single-sided PSD sqrt(Sh)");

  XLALRegisterUvarMember(  lineFeature,          BOOLEAN, 0, OPTIONAL, "Generate a monochromatic 'line' of amplitude h0 and frequency 'Freq'}");

  XLALRegisterUvarMember(parfile,             STRING, 'p', OPTIONAL, "Directory path for optional .par files");            /*registers .par file in mfd*/

  /* transient signal window properties (name, start, duration) */
  XLALRegisterUvarMember(transientWindowType,  STRING, 0, OPTIONAL, "Type of transient signal window to use. ('none', 'rect', 'exp').");
  XLALRegisterUvarMember(  transientStartTime,   REAL8, 0, OPTIONAL, "GPS start-time 't0' of transient signal window.");
  XLALRegisterUvarMember(  transientTauDays,     REAL8, 0, OPTIONAL, "Timescale 'tau' of transient signal window in days.");

  /* ----- 'expert-user/developer' options ----- */
  XLALRegisterUvarMember(  sourceDeltaT,        REAL8,  0, DEVELOPER, "Source-frame sampling period. '0' implies previous internal defaults" );
  XLALRegisterUvarMember(   generationMode,       INT4, 0,  DEVELOPER, "How to generate timeseries: 0=all-at-once (faster), 1=per-sft (slower)");

  XLALRegisterUvarMember(  hardwareTDD,         BOOLEAN, 'b', DEVELOPER, "Hardware injection: output TDD in binary format (implies generationMode=1)");
  XLALRegisterUvarMember(actuation,            STRING, 0,  DEVELOPER, "Filname containing actuation function of this detector");
  XLALRegisterUvarMember(  actuationScale,       REAL8, 0,  DEVELOPER,  "(Signed) scale-factor to apply to the actuation-function.");

  XLALRegisterUvarMember(  exactSignal,          BOOLEAN, 0, DEVELOPER, "Generate signal time-series as exactly as possible (slow).");
  XLALRegisterUvarMember(   randSeed,             INT4, 0, DEVELOPER, "Specify random-number seed for reproducible noise (0 means use /dev/urandom for seeding).");

  /* read cmdline & cfgfile  */
  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit )
    exit (1);

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

  /* free window if any */
  XLALDestroyREAL4Window ( cfg->window );
  XLALFree ( cfg->window_type );

  /* free spindown-vector (REAL8) */
  XLALDestroyREAL8Vector ( cfg->spindown );

  /* free noise-SFTs */
  XLALDestroySFTVector( cfg->noiseSFTs );

  /* free transfer-function if we have one.. */
  XLALDestroyCOMPLEX8FrequencySeries ( cfg->transfer );

  XLALFree ( cfg->VCSInfoString );

  /* Clean up earth/sun Ephemeris tables */
  XLALDestroyEphemerisData ( cfg->edat );

  return XLAL_SUCCESS;

} /* XLALFreeMem() */

/**
 * Log the all relevant parameters of this run into a log-file.
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
      data->data[i-startline] = crectf( cos(phi) / ( amp * actuationScale ), -sin(phi) / ( amp * actuationScale ) );

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
