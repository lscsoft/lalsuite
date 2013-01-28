/*
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
 * This code is a descendant of an earlier implementation 'makefakedata_v2.[ch]'
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

#include <lal/LALString.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/SimulatePulsarSignal.h>
#include <lal/TimeSeries.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/Window.h>

#ifdef HAVE_LIBLALFRAME
#include <lal/LALFrameIO.h>
#endif

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

  REAL8 fmin_eff;		/**< 'effective' fmin: round down such that fmin*Tsft = integer! */
  REAL8 fBand_eff;		/**< 'effective' frequency-band such that samples/SFT is integer */
  REAL8Vector *spindown;	/**< vector of frequency-derivatives of GW signal */

  SFTVector *noiseSFTs;		/**< vector of noise-SFTs to be added to signal */
  REAL8 noiseSigma;		/**< sigma for Gaussian noise to be added */

  REAL4Window *window;		/**< window function for the time series */

  COMPLEX8FrequencySeries *transfer;  /**< detector's transfer function for use in hardware-injection */

  INT4 randSeed;		/**< random-number seed: either taken from user or /dev/urandom */

  transientWindow_t transientWindow;	/**< properties of transient-signal window */
  CHAR *VCSInfoString;          /**< LAL + LALapps Git version string */
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
  BOOLEAN help;		/**< Print this help/usage message */

  /* output */
  CHAR *outSFTbname;		/**< Path and basefilename of output SFT files */
  BOOLEAN outSFTv1;		/**< use v1-spec for output-SFTs */
  BOOLEAN outSingleSFT;	/**< use to output a single concatenated SFT */

  CHAR *TDDfile;		/**< Filename for ASCII output time-series */
  CHAR *TDDframedir;		/**< directory for frame file output time-series */
  CHAR *frameDesc;           	/**< description field entry in the frame filename */

  BOOLEAN hardwareTDD;	/**< Binary output timeseries in chunks of Tsft for hardware injections. */

  CHAR *logfile;		/**< name of logfile */

  /* specify start + duration */
  CHAR *timestampsFile;	/**< Timestamps file */
  INT4 startTime;		/**< Start-time of requested signal in detector-frame (GPS seconds) */
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
  REAL8 noiseSigma;		/**< Gaussian noise with standard-deviation sigma */
  REAL8 noiseSqrtSh;		/**< ALTERNATIVE: single-sided sqrt(Sh) for Gaussian noise */

  /* Window function [OPTIONAL] */
  CHAR *window;		/**< Windowing function for the time series */
  REAL8 tukeyBeta;          /**< Hann fraction of Tukey window (0.0=rect; 1,0=han; 0.5=default */

  /* Detector and ephemeris */
  CHAR *IFO;			/**< Detector: H1, L1, H2, V1, ... */
  CHAR *detector;		/**< [DEPRECATED] Detector: LHO, LLO, VIRGO, GEO, TAMA, CIT, ROME */

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
  REAL8 longitude;		/**< [DEPRECATED] Right ascension [radians] alpha of pulsar */
  REAL8 latitude;		/**< [DEPRECATED] Declination [radians] delta of pulsar */

  REAL8 Freq;
  REAL8 f0;			/**< [DEPRECATED] Gravitational wave-frequency f0 at refTime */

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

  /* precision-level of signal-generation */
  BOOLEAN exactSignal;	/**< generate signal timeseries as exactly as possible (slow) */
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

// ---------- local prototypes ----------
int XLALInitUserVars ( UserVariables_t *uvar, int argc, char *argv[] );
int XLALInitMakefakedata ( ConfigVars_t *cfg, UserVariables_t *uvar );
int XLALWriteMFDlog ( const char *logfile, const ConfigVars_t *cfg );
COMPLEX8FrequencySeries *XLALLoadTransferFunctionFromActuation ( REAL8 actuationScale, const CHAR *fname );
SFTVector *XLALExtractSFTBand ( const SFTVector *inSFTs, REAL8 f_min, REAL8 Band );
int XLALAddGaussianNoise ( REAL4TimeSeries *inSeries, REAL4 sigma, INT4 seed );
int XLALFreeMem ( ConfigVars_t *cfg );
REAL4TimeSeries *XLALGenerateLineFeature ( const PulsarSignalParams *params );

extern void write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series);
extern void write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series);
BOOLEAN is_directory ( const CHAR *fname );
int XLALIsValidDescriptionField ( const char *desc );

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  ConfigVars_t GV = empty_GV;
  PulsarSignalParams params = empty_PulsarSignalParams;
  REAL4TimeSeries *Tseries = NULL;
  UINT4 i_chunk, numchunks;
  FILE *fpSingleSFT = NULL;
  size_t len;
  UserVariables_t uvar = empty_UserVariables;

  lalDebugLevel = 0;	/* default value */

  /* ------------------------------
   * read user-input and set up shop
   *------------------------------*/
  XLAL_CHECK ( XLALGetDebugLevel(argc, argv, 'v') == XLAL_SUCCESS, XLAL_EFUNC );

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
  {
    REAL8 h0   = GV.pulsar.Amp.h0;
    REAL8 cosi = GV.pulsar.Amp.cosi;
    params.pulsar.aPlus		   = 0.5 * h0 * ( 1.0 + SQ(cosi) );
    params.pulsar.aCross	   = h0 * cosi;
  }
  params.pulsar.phi0		   = GV.pulsar.Amp.phi0;
  params.pulsar.psi 		   = GV.pulsar.Amp.psi;

  params.pulsar.f0		   = GV.pulsar.Doppler.fkdot[0];
  params.pulsar.spindown           = GV.spindown;
  params.orbit                     = GV.pulsar.Doppler.orbit;

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
      if ( (fpSingleSFT = LALFopen ( uvar.outSFTbname, "wb" )) == NULL ) {
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
	XLAL_CHECK ( XLALAddGaussianNoise ( Tseries, GV.noiseSigma, GV.randSeed + i_chunk ) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      /* output ASCII time-series if requested */
      if ( uvar.TDDfile )
	{
	  FILE *fp;
	  CHAR *fname = XLALCalloc (1, len = strlen(uvar.TDDfile) + 10 );
          XLAL_CHECK ( fname != NULL, XLAL_ENOMEM, "XLALCalloc(1,%d) failed\n", len );
	  sprintf (fname, "%s.%02d", uvar.TDDfile, i_chunk);

	  if ( (fp = fopen (fname, "w")) == NULL)
	    {
	      perror ("Error opening outTDDfile for writing");
              XLAL_ERROR ( XLAL_EIO, "Failed to fopen TDDfile = '%s' for writing\n", fname );
	    }

	  write_timeSeriesR4(fp, Tseries);
	  fclose(fp);
	  XLALFree (fname);
	} /* if outputting ASCII time-series */

      /* output time-series to frames if requested */
      if ( uvar.TDDframedir )
	{
#ifndef HAVE_LIBLALFRAME
          XLAL_ERROR ( XLAL_EINVAL, "--TDDframedir option not supported, code has to be compiled with lalframe\n" );
#else
	  /* use standard frame output filename format */
          XLAL_CHECK ( XLALIsValidDescriptionField ( uvar.frameDesc ) == XLAL_SUCCESS, XLAL_EFUNC );
          len = strlen(uvar.TDDframedir) + strlen(uvar.frameDesc) + 100;
	  char *fname;
          char IFO[2] = { Tseries->name[0], Tseries->name[1] };
          XLAL_CHECK ( (fname = LALCalloc (1, len )) != NULL, XLAL_ENOMEM );
          size_t written = snprintf ( fname, len, "%s/%c-%c%c_%s-%d-%d.gwf",
                                      uvar.TDDframedir, IFO[0], IFO[0], IFO[1], uvar.frameDesc, params.startTimeGPS.gpsSeconds, (int)params.duration );
          XLAL_CHECK ( written < len, XLAL_ESIZE, "Frame-filename exceeds expected maximal length (%d): '%s'\n", len, fname );

	  /* define the output frame */
	  struct FrameH *outFrame;
	  XLAL_CHECK ( (outFrame = XLALFrameNew( &(params.startTimeGPS), params.duration, uvar.frameDesc, 1, 0, 0 )) != NULL, XLAL_EFUNC );

	  /* add timeseries to the frame - make sure to change the timeseries name since this is used as the channel name */
          char buffer[LALNameLength];
	  written = snprintf ( buffer, LALNameLength, "%s:%s", Tseries->name, uvar.frameDesc );
          XLAL_CHECK ( written < LALNameLength, XLAL_ESIZE, "Updated frame name exceeds max length (%d): '%s'\n", LALNameLength, buffer );
          strcpy ( Tseries->name, buffer );

	  XLAL_CHECK ( (XLALFrameAddREAL4TimeSeriesProcData ( outFrame, Tseries ) == XLAL_SUCCESS ) , XLAL_EFUNC );

	  /* Here's where we add extra information into the frame - first we add the command line args used to generate it */
	  char *hist = XLALUserVarGetLog (UVAR_LOGFMT_CMDLINE);
          FrHistoryAdd ( outFrame, hist );

	  /* then we add the version string */
	  FrHistoryAdd ( outFrame, GV.VCSInfoString );

	  /* output the frame to file - compression level 1 (higher values make no difference) */
	  XLAL_CHECK ( (XLALFrameWrite(outFrame, fname,1) == 0) , XLAL_EFUNC );

	  /* free the frame, frame file name and history memory */
	  FrameFree ( outFrame );
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
	  SFTParams sftParams = empty_SFTParams;
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

	  /* extract requested band if necessary (eg in the exact-signal case) */
	  if ( uvar.exactSignal )
	    {
	      SFTVector *outSFTs;
              XLAL_CHECK ( (outSFTs = XLALExtractSFTBand ( SFTs, GV.fmin_eff, GV.fBand_eff )) != NULL, XLAL_EFUNC );
	      XLALDestroySFTVector ( SFTs );
	      SFTs = outSFTs;
	    }

	  if ( uvar.outSFTv1 ) 		/* write output-SFTs using the SFT-v1 format */
	    {
	      CHAR *fname = XLALCalloc (1, len = strlen (uvar.outSFTbname) + 10 );
              XLAL_CHECK ( fname != NULL, XLAL_ENOMEM, "XLALCalloc(1,%d) failed.\n", len );

              LALStatus status = blank_status;
	      for (UINT4 i=0; i < SFTs->length; i++)
		{
		  sprintf (fname, "%s.%05d", uvar.outSFTbname, i_chunk*SFTs->length + i);
		  LALWrite_v2SFT_to_v1file ( &status, &(SFTs->data[i]), fname );
                  XLAL_CHECK ( status.statusCode == 0, XLAL_EFAILED, "LALWrite_v2SFT_to_v1file('%s') failed with status=%d : '%s'\n",
                               fname, status.statusCode, status.statusDescription );
		}
	      XLALFree (fname);
	    } /* if outSFTv1 */
	  else
	    {	/* write standard v2-SFTs */

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
                      int ret = XLALWriteSFT2fp ( &(SFTs->data[k]), fpSingleSFT, comment );
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

            } // if v2-SFTs
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

  BOOLEAN have_parfile = LALUserVarWasSet (&uvar->parfile);
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
    BOOLEAN have_h0     = LALUserVarWasSet ( &uvar->h0 );
    BOOLEAN have_cosi   = LALUserVarWasSet ( &uvar->cosi );
    BOOLEAN have_aPlus  = LALUserVarWasSet ( &uvar->aPlus );
    BOOLEAN have_aCross = LALUserVarWasSet ( &uvar->aCross );
    BOOLEAN have_Freq   = LALUserVarWasSet ( &uvar->Freq );
    BOOLEAN have_f0     = LALUserVarWasSet ( &uvar->f0 );
    BOOLEAN have_longitude = LALUserVarWasSet ( &uvar->longitude );
    BOOLEAN have_latitude  = LALUserVarWasSet ( &uvar->latitude );
    BOOLEAN have_Alpha  = LALUserVarWasSet ( &uvar->Alpha );
    BOOLEAN have_Delta  = LALUserVarWasSet ( &uvar->Delta );
    BOOLEAN have_RA = LALUserVarWasSet ( &uvar->RA );
    BOOLEAN have_Dec = LALUserVarWasSet ( &uvar->Dec );

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
    if ( have_f0 && have_Freq ) {
      XLAL_ERROR ( XLAL_EINVAL, "Specify signal-frequency using EITHER --Freq [preferred] OR --f0 [deprecated]!\n\n");
    }
    if ( have_Freq )
      cfg->pulsar.Doppler.fkdot[0] = uvar->Freq;
    else if ( have_f0 )
      cfg->pulsar.Doppler.fkdot[0] = uvar->f0;
    else
      cfg->pulsar.Doppler.fkdot[0] = 0.0;

    /* ----- skypos ----- */
    if ( (have_Alpha || have_Delta) && ( have_longitude || have_latitude ) && (have_RA || have_Dec) ) {
      XLAL_ERROR ( XLAL_EINVAL, "Use EITHER {Alpha, Delta} [preferred] OR {RA, Dec} [preferred] OR {longitude,latitude} [deprecated]\n\n");
    }
    if ( (have_Alpha && !have_Delta) || ( !have_Alpha && have_Delta ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "\nSpecify skyposition: need BOTH --Alpha and --Delta!\n\n");
    }
    if ( (have_longitude && !have_latitude) || ( !have_longitude && have_latitude ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "\nSpecify skyposition: need BOTH --longitude and --latitude!\n\n");
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
	cfg->pulsar.Doppler.Alpha = XLALhmsToRads(uvar->RA);
	cfg->pulsar.Doppler.Delta = XLALdmsToRads(uvar->Dec);
      }
    else
      {
	cfg->pulsar.Doppler.Alpha = uvar->longitude;
	cfg->pulsar.Doppler.Delta = uvar->latitude;
      }

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
      case 2:
	cfg->spindown->data[1] = uvar->f2dot;
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
    LALDetector *site;
    BOOLEAN have_detector = LALUserVarWasSet ( &uvar->detector );
    BOOLEAN have_IFO = LALUserVarWasSet ( &uvar->IFO );
    if ( !have_detector && !have_IFO ) {
      XLAL_ERROR ( XLAL_EINVAL, "Need detector input --IFO!\n\n");
    }
    if ( have_detector && have_IFO ) {
      XLAL_ERROR ( XLAL_EINVAL, "\nUse EITHER --IFO [preferred] OR --detector [deprecated]!\n\n");
    }

    if ( have_detector )
      {
	site = XLALGetSiteInfo ( uvar->detector );
	channelName = XLALGetChannelPrefix ( uvar->detector );
      }
    else
      {
	site = XLALGetSiteInfo ( uvar->IFO );
        XLAL_CHECK ( site != NULL, XLAL_EFUNC, "XLALGetSiteInfo('%s') failed\n", uvar->IFO );
	channelName = XLALGetChannelPrefix ( uvar->IFO );
        XLAL_CHECK ( channelName != NULL, XLAL_EFUNC, "XLALGetChannelPrefix('%s') failed.\n", uvar->IFO );
      }

    cfg->site = (*site);
    XLALFree ( site );
  }

   /* check for negative fmin and Band, which would break the fmin_eff, fBand_eff calculation below */
  XLAL_CHECK ( uvar->fmin >= 0, XLAL_EDOM, "Invalid negative frequency fmin=%f!\n\n", uvar->fmin );
  XLAL_CHECK ( uvar->Band >= 0, XLAL_EDOM, "Invalid negative frequency band Band=%f!\n\n", uvar->Band );

  /* ---------- for SFT output: calculate effective fmin and Band ---------- */
  if ( XLALUserVarWasSet( &uvar->outSFTbname ) )
    {
      UINT4 imin, imax;
      volatile REAL8 dFreq = 1.0 / uvar->Tsft;
      volatile REAL8 tmp;
      REAL8 fMax, fMin_eff;

      /* calculate "effective" fmin from uvar->fmin: following makefakedata_v2, we
       * make sure that fmin_eff * Tsft = integer, such that freqBinIndex corresponds
       * to a frequency-index of the non-heterodyned signal.
       */
      tmp = uvar->fmin / dFreq;	/* NOTE: don't "simplify" this: we try to make sure
				 * the result of this will be guaranteed to be IEEE-compliant,
				 * and identical to other locations, such as in SFT-IO
				 */
      imin = (UINT4) floor( tmp );
      fMin_eff = (REAL8)imin * dFreq;

      fMax = uvar->fmin + uvar->Band;
      tmp = fMax / dFreq;
      imax = (UINT4) ceil (tmp);

      /* Increase Band correspondingly. */
      cfg->fmin_eff = fMin_eff;
      cfg->fBand_eff = 1.0 * (imax - imin) * dFreq;

      if ( lalDebugLevel )
	{
	  if ( abs(cfg->fmin_eff - uvar->fmin)> LAL_REAL8_EPS
	       || abs(cfg->fBand_eff - uvar->Band) > LAL_REAL8_EPS )
	    printf("\nWARNING: for SFT-creation we had to adjust (fmin,Band) to"
		   " fmin_eff=%.20g and Band_eff=%.20g\n\n", cfg->fmin_eff, cfg->fBand_eff);
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

    haveStart = LALUserVarWasSet(&uvar->startTime);
    haveDuration = LALUserVarWasSet(&uvar->duration);
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
	if ( LALUserVarWasSet( &uvar->outSFTbname ) ) {
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
	SFTConstraints constraints = empty_SFTConstraints;
	LIGOTimeGPS minStartTime, maxEndTime;
        BOOLEAN have_window = LALUserVarWasSet ( &uvar->window );

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

  if ( uvar->window )
    {
      XLALStringToLowerCase ( uvar->window );	// get rid of case

      if ( LALUserVarWasSet( &uvar->tukeyBeta ) && strcmp ( uvar->window, "tukey" ) ) {
        XLAL_ERROR ( XLAL_EINVAL, "Tukey beta value '%f' was specified with window %s; only allowed for Tukey windowing.\n\n", uvar->tukeyBeta, uvar->window );
      }

      /* NOTE: a timeseries of length N*dT has no timestep at N*dT !! (convention) */
      UINT4 lengthOfTimeSeries = (UINT4)round(uvar->Tsft * 2 * cfg->fBand_eff);

      if ( !strcmp ( uvar->window, "hann" ) || !strcmp ( uvar->window, "hanning" ) )
        {
          REAL4Window *win = XLALCreateHannREAL4Window( lengthOfTimeSeries );
          cfg->window = win;
        }
      else if ( !strcmp ( uvar->window, "tukey" ) )
	{
	  if ( !LALUserVarWasSet( &uvar->tukeyBeta ) )
	    {
	      uvar->tukeyBeta = 0.5;   /* If Tukey window specified, default transition fraction is 1/2 */
	    }
	  else if ( uvar->tukeyBeta < 0.0 || uvar->tukeyBeta > 1.0 ) {
            XLAL_ERROR ( XLAL_EINVAL, "Tukey beta value '%f' was specified; must be between 0 and 1.\n\n", uvar->tukeyBeta );
          }
          cfg->window = XLALCreateTukeyREAL4Window( lengthOfTimeSeries, uvar->tukeyBeta );
          XLAL_CHECK ( cfg->window != NULL, XLAL_EFUNC, "XLALCreateTukeyREAL4Window(%d, %g) failed\n", lengthOfTimeSeries, uvar->tukeyBeta );
	}
      else if ( !strcmp ( uvar->window, "none" ) || !strcmp ( uvar->window, "rectangular" ) || !strcmp ( uvar->window, "boxcar" ) || !strcmp ( uvar->window, "tophat" ) ) {
        cfg->window = NULL;
      }
      else
        {
          XLAL_ERROR ( XLAL_EINVAL, "Invalid window function '%s', allowed are ['None', 'Hann', or 'Tukey'].\n\n", uvar->window );
        }
    }
  else
    {
      if ( LALUserVarWasSet( &uvar->tukeyBeta ) ) {
        XLAL_ERROR ( XLAL_EINVAL, "Tukey beta value '%f' was specified; only relevant if Tukey windowing specified.\n\n", uvar->tukeyBeta );
      }
    } /* if uvar->window */

  /* -------------------- Prepare quantities for barycentering -------------------- */
  {
    CHAR *earthdata, *sundata;

    len = strlen(uvar->ephemYear) + 20;

    if (LALUserVarWasSet(&uvar->ephemDir) )
      len += strlen (uvar->ephemDir);

    if ( (earthdata = XLALCalloc(1, len)) == NULL) {
      XLAL_ERROR ( XLAL_ENOMEM, "earthdata = XLALCalloc(1, %d) failed.\n", len );
    }
    if ( (sundata = XLALCalloc(1, len)) == NULL) {
      XLAL_ERROR ( XLAL_ENOMEM, "sundata = XLALCalloc(1, %d) failed.\n", len );
    }

    if (LALUserVarWasSet(&uvar->ephemDir) )
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
    BOOLEAN set1 = LALUserVarWasSet(&uvar->orbitasini);
    BOOLEAN set2 = LALUserVarWasSet(&uvar->orbitEcc);
    BOOLEAN set3 = LALUserVarWasSet(&uvar->orbitPeriod);
    BOOLEAN set4 = LALUserVarWasSet(&uvar->orbitArgp);
    BOOLEAN set5 = LALUserVarWasSet(&uvar->orbitTpSSBsec);
    BOOLEAN set6 = LALUserVarWasSet(&uvar->orbitTpSSBnan);
    BOOLEAN set7 = LALUserVarWasSet(&uvar->orbitTpSSBMJD);
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

      cfg->pulsar.Doppler.orbit = orbit;     /* struct copy */
    } /* if one or more orbital parameters were set */
    else
      cfg->pulsar.Doppler.orbit = NULL;
  } /* END: binary orbital params */


  /* -------------------- handle NOISE params -------------------- */
  {
    if ( LALUserVarWasSet ( &uvar->noiseSigma ) && LALUserVarWasSet ( &uvar->noiseSqrtSh ) )
      {
	XLAL_ERROR ( XLAL_EINVAL, "Use only one of '--noiseSigma' and '--noiseSqrtSh' to specify Gaussian noise!\n\n");
      }
    if ( LALUserVarWasSet ( &uvar->noiseSigma ) )
      cfg->noiseSigma = uvar->noiseSigma;
    else if ( LALUserVarWasSet ( &uvar->noiseSqrtSh ) )	/* convert Sh -> sigma */
      cfg->noiseSigma = uvar->noiseSqrtSh * sqrt ( cfg->fBand_eff );
    else
      cfg->noiseSigma = 0;

    /* set random-number generator seed: either taken from user or from /dev/urandom */
    if ( LALUserVarWasSet ( &uvar->randSeed ) )
      {
        if ( uvar->randSeed == 0 ) {
          XLALPrintError ("WARNING: setting randSeed==0 results in the system clock being used as a random seed!\n");
        }
        cfg->randSeed = uvar->randSeed;
      }
    else
      {
	/*
	 * Modified so as to not create random number parameters with seed
	 * drawn from clock.  Seconds don't change fast enough and sft's
	 * look alike.  We open /dev/urandom and read a 4 byte integer from
	 * it and use that as our seed.  Note: /dev/random is slow after the
	 * first, few accesses.
	 */
	FILE *devrandom = fopen ( "/dev/urandom", "r" );
        XLAL_CHECK ( devrandom != NULL, XLAL_EIO, "fopen() failed to open '/dev/urandom' for reading\n\n");

	if ( fread( (void*)&(cfg->randSeed), sizeof(INT4), 1, devrandom) != 1 )
	  {
	    fclose ( devrandom );
	    XLAL_ERROR ( XLAL_EIO, "Failed to read from '/dev/urandom'\n\n");
	  }
	fclose ( devrandom );
      }

  } /* END: Noise params */


  /* ----- set "pulsar reference time", i.e. SSB-time at which pulsar params are defined ---------- */
  if (LALUserVarWasSet (&uvar->parfile)) {
    uvar->refTime = pulparams.pepoch; /*XLALReadTEMPOParFile already converted pepoch to GPS*/
    XLALGPSSetREAL8(&(cfg->pulsar.Doppler.refTime),uvar->refTime);
  }
  else if (LALUserVarWasSet(&uvar->refTime) && LALUserVarWasSet(&uvar->refTimeMJD))
    {
      XLAL_ERROR ( XLAL_EINVAL, "\nUse only one of '--refTime' and '--refTimeMJD' to specify SSB reference time!\n\n");
    }
  else if (LALUserVarWasSet(&uvar->refTime))
    {
      XLALGPSSetREAL8(&(cfg->pulsar.Doppler.refTime), uvar->refTime);
    }
  else if (LALUserVarWasSet(&uvar->refTimeMJD))
    {

      /* convert MJD to GPS using Matt Pitkins code found at lal/packages/pulsar/src/BinaryPulsarTimeing.c */
      REAL8 GPSfloat;
      GPSfloat = XLALTTMJDtoGPS(uvar->refTimeMJD);
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

  if ( !uvar->actuation && LALUserVarWasSet(&uvar->actuationScale) ) {
    XLAL_ERROR ( XLAL_EINVAL, "Actuation-scale was specified without actuation-function file!\n\n");
  }

  XLALFree ( channelName );

  /* ----- handle transient-signal window if given ----- */
  if ( !LALUserVarWasSet ( &uvar->transientWindowType ) || !strcmp ( uvar->transientWindowType, "none") )
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

  // per default we now generate a timeseries per SFT: slower, but avoids potential confusion about sft-"nudging"
  uvar->generationMode = GENERATE_PER_SFT;

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

  XLALregSTRINGUserStruct( TDDfile,		't', UVAR_OPTIONAL, "Basename for output of ASCII time-series");
  XLALregSTRINGUserStruct( TDDframedir,		'F', UVAR_OPTIONAL, "Directory to output frame time-series into");
  XLALregSTRINGUserStruct( frameDesc,	 	 0,  UVAR_OPTIONAL,  "Description-field entry in frame filename");

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
  XLALregSTRINGUserStruct (window,               0, UVAR_OPTIONAL, "Window function to be applied to the SFTs; required when using --noiseSFTs ('None', 'Hann', or 'Tukey'; when --noiseSFTs is not given, default is 'None')");
  XLALregREALUserStruct (  tukeyBeta,            0, UVAR_OPTIONAL, "Fraction of Tukey window which is transition (0.0=rect, 1.0=Hann)");

  /* pulsar params */
  XLALregREALUserStruct (  refTime,             'S', UVAR_OPTIONAL, "Pulsar SSB reference time in GPS seconds (default: use startTime)");
  XLALregREALUserStruct (  refTimeMJD,           0 , UVAR_OPTIONAL, "ALTERNATIVE: Pulsar SSB reference time in MJD (default: use startTime)");

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
  XLALregREALUserStruct (  noiseSigma,           0, UVAR_OPTIONAL,  "Gaussian noise with standard-deviation sigma");
  XLALregREALUserStruct (  noiseSqrtSh,          0, UVAR_OPTIONAL,  "ALTERNATIVE: Gaussian noise with single-sided PSD sqrt(Sh)");

  XLALregBOOLUserStruct (  lineFeature,          0, UVAR_OPTIONAL, "Generate a monochromatic 'line' of amplitude h0 and frequency 'Freq'}");

  XLALregBOOLUserStruct (  version,             'V', UVAR_SPECIAL, "Output version information");

  XLALregSTRINGUserStruct (parfile,             'p', UVAR_OPTIONAL, "Directory path for optional .par files");            /*registers .par file in mfd*/

  /* transient signal window properties (name, start, duration) */
  XLALregSTRINGUserStruct (transientWindowType,  0, UVAR_OPTIONAL, "Type of transient signal window to use. ('none', 'rect', 'exp').");
  XLALregREALUserStruct (  transientStartTime,   0, UVAR_OPTIONAL, "GPS start-time 't0' of transient signal window.");
  XLALregREALUserStruct (  transientTauDays,     0, UVAR_OPTIONAL, "Timescale 'tau' of transient signal window in days.");

  /* ----- 'expert-user/developer' and deprecated options ----- */
  XLALregINTUserStruct (   generationMode,       0,  UVAR_DEVELOPER, "How to generate timeseries: 0=all-at-once (faster), 1=per-sft (slower)");

  XLALregBOOLUserStruct (  hardwareTDD,         'b', UVAR_DEVELOPER, "Hardware injection: output TDD in binary format (implies generationMode=1)");
  XLALregSTRINGUserStruct (actuation,            0,  UVAR_DEVELOPER, "Filname containing actuation function of this detector");
  XLALregREALUserStruct (  actuationScale,       0,  UVAR_DEVELOPER,  "(Signed) scale-factor to apply to the actuation-function.");

  XLALregBOOLUserStruct (  exactSignal,          0, UVAR_DEVELOPER, "Generate signal time-series as exactly as possible (slow).");
  XLALregINTUserStruct (   randSeed,             0, UVAR_DEVELOPER, "Specify random-number seed for reproducible noise (use /dev/urandom otherwise).");

  XLALregREALUserStruct (  longitude,            0, UVAR_DEVELOPER, "[DEPRECATED] Use --Alpha instead!");
  XLALregREALUserStruct (  latitude,             0, UVAR_DEVELOPER, "[DEPRECATED] Use --Delta instead!");
  XLALregREALUserStruct (  f0,                   0, UVAR_DEVELOPER, "[DEPRECATED] Use --Freq instead!");
  XLALregSTRINGUserStruct (detector,             0,  UVAR_DEVELOPER, "[DEPRECATED] Detector: use --IFO instead!.");
  XLALregBOOLUserStruct( outSFTv1,	 	 0, UVAR_DEVELOPER, "[DEPRECATED]Write output-SFTs in obsolete SFT-v1 format." );

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

  /* free window if any */
  XLALDestroyREAL4Window ( cfg->window );

  /* free spindown-vector (REAL8) */
  XLALDestroyREAL8Vector ( cfg->spindown );

  XLALFree ( cfg->pulsar.Doppler.orbit );

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
 * Generate Gaussian noise with standard-deviation sigma, add it to inSeries.
 *
 * NOTE2: if seed==0, then time(NULL) is used as random-seed!
 *
 */
int
XLALAddGaussianNoise ( REAL4TimeSeries *inSeries, REAL4 sigma, INT4 seed )
{
  XLAL_CHECK ( inSeries != NULL, XLAL_EINVAL );

  UINT4 numPoints = inSeries->data->length;

  REAL4Vector *v1;
  XLAL_CHECK ( (v1 = XLALCreateREAL4Vector ( numPoints )) != NULL, XLAL_EFUNC );

  RandomParams *randpar;
  XLAL_CHECK ( (randpar = XLALCreateRandomParams ( seed )) != NULL, XLAL_EFUNC );

  XLAL_CHECK ( XLALNormalDeviates ( v1, randpar) == XLAL_SUCCESS, XLAL_EFUNC );

  for (UINT4 i = 0; i < numPoints; i++ ) {
    inSeries->data->data[i] += sigma * v1->data[i];
  }

  /* destroy randpar*/
  XLALDestroyRandomParams ( randpar );

  /*   destroy v1 */
  XLALDestroyREAL4Vector ( v1 );

  return XLAL_SUCCESS;

} /* XLALAddGaussianNoise() */

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


/** Return a vector of SFTs containg only the bins in [fmin, fmin+Band].
 */
SFTVector *
XLALExtractSFTBand ( const SFTVector *inSFTs, REAL8 f_min, REAL8 Band )
{
  XLAL_CHECK_NULL ( inSFTs != NULL, XLAL_EINVAL, "Invalid NULL input SFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( inSFTs->length > 0, XLAL_EINVAL, "Invalid zero-length input SFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( f_min >= 0, XLAL_EDOM, "Invalid negative frequency f_min = %g\n", f_min );
  XLAL_CHECK_NULL ( Band > 0, XLAL_EDOM, "Invalid non-positive Band = %g\n", Band );

  UINT4 numSFTs = inSFTs->length;
  REAL8 SFTf0   = inSFTs->data[0].f0;
  REAL8 df      = inSFTs->data[0].deltaF;
  REAL8 SFTBand = df * inSFTs->data[0].data->length;

  XLAL_CHECK_NULL ( (f_min >= SFTf0) && ( f_min + Band <= SFTf0 + SFTBand ), XLAL_EINVAL,
                    "Requested frequency-band [%f,%f] Hz not contained SFTs [%f, %f] Hz.\n", f_min, f_min + Band, SFTf0, SFTf0 + SFTBand );

  UINT4 firstBin = round ( f_min / df );
  UINT4 numBins =  round ( Band / df ) + 1;

  SFTVector *ret = XLALCreateSFTVector ( numSFTs, numBins );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_EFUNC, "XLALCreateSFTVector ( %d, %d ) failed.\n", numSFTs, numBins );

  for ( UINT4 i = 0; i < numSFTs; i ++ )
    {
      SFTtype *dest = &(ret->data[i]);
      SFTtype *src =  &(inSFTs->data[i]);
      COMPLEX8Vector *ptr = dest->data;

      /* copy complete header first */
      memcpy ( dest, src, sizeof(*dest) );
      /* restore data-pointer */
      dest->data = ptr;
      /* set correct f_min */
      dest->f0 = firstBin * df ;

      /* copy the relevant part of the data */
      memcpy ( dest->data->data, src->data->data + firstBin, numBins * sizeof( dest->data->data[0] ) );

    } /* for i < numSFTs */

  /* return final SFT-vector */
  return ret;

} /* XLALExtractSFTBand() */


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
 * Generate a REAL4TimeSeries containing a sinusoid with
 * amplitude 'h0', frequency 'Freq-fHeterodyne' and initial phase 'phi0'.
 */
REAL4TimeSeries *
XLALGenerateLineFeature ( const PulsarSignalParams *params )
{
  XLAL_CHECK_NULL ( params != NULL, XLAL_EINVAL );

  /* set 'name'-field of timeseries to contain the right "channel prefix" for the detector */
  char *name;
  XLAL_CHECK_NULL ( (name = XLALGetChannelPrefix ( params->site->frDetector.name )) != NULL, XLAL_EFUNC );

  /* NOTE: a timeseries of length N*dT has no timestep at N*dT !! (convention) */
  UINT4 length = (UINT4) ceil( params->samplingRate * params->duration);
  REAL8 deltaT = 1.0 / params->samplingRate;
  REAL8 tStart = XLALGPSGetREAL8 ( &params->startTimeGPS );

  LALUnit units = empty_LALUnit;
  REAL4TimeSeries *ret = XLALCreateREAL4TimeSeries (name, &(params->startTimeGPS), params->fHeterodyne, deltaT, &units, length);
  XLAL_CHECK_NULL ( ret != NULL, XLAL_EFUNC, "XLALCreateREAL4TimeSeries() failed.\n");

  REAL8 h0 = params->pulsar.aPlus + sqrt ( pow(params->pulsar.aPlus,2) - pow(params->pulsar.aCross,2) );
  REAL8 omH = LAL_TWOPI * ( params->pulsar.f0 - params->fHeterodyne );

  for ( UINT4 i = 0; i < length; i++ )
    {
      REAL8 ti = tStart + i * deltaT;
      ret->data->data[i] = h0 * sin( omH * ti  + params->pulsar.phi0 );
    }

  /* return final timeseries */
  return ret;

} /* XLALGenerateLineFeature() */


/**
 * Check whether given string qualifies as a valid 'description' field of a FRAME (or SFT)
 * filename, according to  LIGO-T010150-00-E "Naming Convention for Frame Files which are to be Processed by LDAS",
 * LIGO-T040164-01 at https://dcc.ligo.org/LIGO-T040164-x0/public
 *
 * NOTE: this function will be moved into SFTutils.c later
 */
#include <ctype.h>
int
XLALIsValidDescriptionField ( const char *desc )
{
  XLAL_CHECK ( desc != NULL, XLAL_EINVAL );

  size_t len = strlen ( desc );

  if ( len == 1 && isupper(desc[0]) ) {
    XLAL_ERROR ( XLAL_EINVAL, "Single uppercase description reserved for class-1 raw frames!\n" );
  }

  for ( UINT4 i=0; i < len; i ++ )
    {
      int c = desc[i];
      if ( !isalnum(c) && (c!='_') && (c!='+') && (c!='#') ) {	// all the valid characters allowed
        XLAL_ERROR ( XLAL_EINVAL, "Invalid chacter '%c' found, only alphanumeric and ['_', '+', '#'] are allowed\n", c );
      }
    } // for i < len

  return XLAL_SUCCESS;

} // XLALIsValidDescriptionField()
