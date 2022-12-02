/*
*  Copyright (C) 2013, 2015 Reinhard Prix
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
#include "config.h"

#include <sys/stat.h>

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
#include <lal/LALString.h>
#include <lal/LALCache.h>
#include <lal/Units.h>
#include <lal/TransientCW_utils.h>
#include <lal/CWMakeFakeData.h>
#include <lal/LALPulsarVCSInfo.h>

#ifdef HAVE_LIBLALFRAME
#include <lal/LALFrameIO.h>
#include <lal/LALFrStream.h>
#endif

/***************************************************/

/*----------------------------------------------------------------------*/


/** configuration-variables derived from user-variables */
typedef struct
{
  EphemerisData *edat;				/**< ephemeris-data */
  MultiLIGOTimeGPSVector *multiTimestamps;	/**< a vector of timestamps to generate time-series/SFTs for */
  MultiLALDetector multiIFO;			//!< detectors to generate data for
  MultiNoiseFloor multiNoiseFloor;		//!< ... and corresponding noise-floors to generate Gaussian white noise for

  SFTCatalog *noiseCatalog; 			/**< catalog of noise-SFTs */
  MultiSFTCatalogView *multiNoiseCatalogView; 	/**< multi-IFO 'view' of noise-SFT catalogs */
  MultiREAL8TimeSeries *inputMultiTS;	/**< 'input' time-series to add other stuff to, and output as frames or SFTs */
  REAL8 fminOut;				/**< Lowest frequency in output SFT (= heterodyning frequency) */
  REAL8 BandOut;				/**< bandwidth of output SFT in Hz (= 1/2 sampling frequency) */

  const CHAR *window_type;			/**< name of window function */
  REAL8 window_param;				/**< parameter of window function */

  transientWindow_t transientWindow;	/**< properties of transient-signal window */
  CHAR *VCSInfoString;			/**< Git version string */
  CHAR *outFrameDir;			/**< output frame directory */
} ConfigVars_t;

// ----- User variables
typedef struct
{
  /* SFT output */
  CHAR *outSFTdir;		/**< Output directory for SFTs */
  CHAR *outLabel;		/**< 'misc' entry in SFT-filenames, and description entry of output frame filenames */
  UINT4 outPubObsRun;           /**< if >0, names SFTs using the public filename convention with this observing run number */
  UINT4 outPubRevision;         /**< if outPubObsRun>0, names SFTs using the public filename convention with this revision number */
  BOOLEAN outSingleSFT;	        /**< use to output a single concatenated SFT */

  CHAR *TDDfile;		/**< Filename for ASCII output time-series */
  CHAR *logfile;		/**< name of logfile */

  /* specify start + duration */
  LALStringVector *timestampsFiles;        /**< Names of numDet timestamps files */
  LIGOTimeGPS startTime;	/**< Start-time of requested signal in detector-frame (GPS seconds) */
  INT4 duration;		/**< Duration of requested signal in seconds */

  /* time-series sampling + heterodyning frequencies */
  REAL8 fmin;		/**< Lowest frequency in output SFT (= heterodyning frequency) */
  REAL8 Band;		/**< bandwidth of output SFT in Hz (= 1/2 sampling frequency) */
  REAL8 sourceDeltaT;   /**< source-frame sampling period. '0' implies defaults set in XLALGeneratePulsarSignal() */

  /* SFT params */
  REAL8 Tsft;		        /**< SFT time baseline Tsft */
  REAL8 SFToverlap;		/**< overlap SFTs by this many seconds */

  LALStringVector* IFOs;	/**< list of detector-names "H1,H2,L1,.." or single detector*/

  /* noise to add [OPTIONAL] */
  LALStringVector* sqrtSX; 	/**< Add Gaussian noise: list of respective detectors' noise-floors sqrt{Sn}" */

  CHAR *noiseSFTs;		/**< Glob-like pattern specifying noise-SFTs to be added to signal */

  /* Frame input and output */
  CHAR *outFrameDir;		/**< directory for writing output timeseries in frame files */
  LALStringVector *inFrames;	/**< frame glob-patterns or cache files of time-series data to be added to output (one per IFO) */
  LALStringVector *inFrChannels;/**< list of frame channels to read time-series data from (one per IFO) */
  LALStringVector *outFrChannels;/**< list of output frame channel names to write (one per IFO) */

  /* Window function [OPTIONAL] */
  CHAR *SFTWindowType;		/**< Windowing function to apply to the SFT time series */
  REAL8 SFTWindowParam;        	/**< parameter required for certain window-types */

  CHAR *ephemEarth;		/**< Earth ephemeris file to use */
  CHAR *ephemSun;		/**< Sun ephemeris file to use */

  /* pulsar parameters */
  LALStringVector *injectionSources;	///< Source parameters to inject: comma-separated list of file-patterns and/or direct config-strings ('{...}')

  INT4 randSeed;		/**< allow user to specify random-number seed for reproducible noise-realizations */

} UserVariables_t;


// ---------- exportable API types ----------

// ----- global variables ----------

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
  size_t len;
  ConfigVars_t XLAL_INIT_DECL(GV);
  UserVariables_t XLAL_INIT_DECL(uvar);

  /* ------------------------------
   * read user-input and set up shop
   *------------------------------*/
  XLAL_CHECK ( XLALInitUserVars ( &uvar, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALInitMakefakedata ( &GV, &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  MultiSFTVector *mSFTs = NULL;
  MultiREAL8TimeSeries *mTseries = NULL;

  PulsarParamsVector *injectionSources = NULL;
  if ( uvar.injectionSources ) {
    XLAL_CHECK ( (injectionSources = XLALPulsarParamsFromUserInput ( uvar.injectionSources, NULL ) ) != NULL, XLAL_EFUNC );
  }

  CWMFDataParams XLAL_INIT_DECL(DataParams);
  DataParams.multiIFO           = GV.multiIFO;
  DataParams.multiNoiseFloor    = GV.multiNoiseFloor;
  DataParams.multiTimestamps 	= GV.multiTimestamps;
  DataParams.randSeed           = uvar.randSeed;
  DataParams.SFTWindowType      = GV.window_type;
  DataParams.SFTWindowParam     = GV.window_param;
  DataParams.sourceDeltaT       = uvar.sourceDeltaT;
  DataParams.inputMultiTS       = GV.inputMultiTS;
  DataParams.fMin               = GV.fminOut;
  DataParams.Band               = GV.BandOut;

  XLAL_CHECK ( XLALCWMakeFakeMultiData ( &mSFTs, &mTseries, injectionSources, &DataParams, GV.edat ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALDestroyPulsarParamsVector ( injectionSources );
  injectionSources = NULL;

  // if noiseSFTs specified, load them and add them to the resulting SFT-vector
  if ( GV.multiNoiseCatalogView )
    {
      SFTtype *sft0 = &(mSFTs->data[0]->data[0]);
      /* load effective frequency-band from noise-SFTs */
      UINT4 numBins = sft0->data->length;
      REAL8 dFreq   = sft0->deltaF;
      REAL8 fMin    = sft0->f0;
      REAL8 fMax    = fMin + numBins * dFreq;
      MultiSFTVector *mNoiseSFTs;
      XLAL_CHECK ( (mNoiseSFTs = XLALLoadMultiSFTsFromView ( GV.multiNoiseCatalogView, fMin, fMax )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( XLALMultiSFTVectorAdd ( mSFTs, mNoiseSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroyMultiSFTVector ( mNoiseSFTs );
    }

  // determine output channel names for frames and public SFT filenames
#ifdef HAVE_LIBLALFRAME
  if ( XLALUserVarWasSet ( &uvar.outFrChannels ) ) {
    XLAL_CHECK ( uvar.outFrChannels->length == mTseries->length, XLAL_EINVAL, "--outFrChannels: number of channel names (%d) must agree with number of IFOs (%d)\n",
                 uvar.outFrChannels->length, mTseries->length );
  }
#endif
  LALStringVector *outChannelNames = XLALCreateEmptyStringVector( mTseries->length );
  XLAL_CHECK ( outChannelNames != NULL, XLAL_EFUNC );
  for ( UINT4 X=0; X < mTseries->length; X ++ )
    {
      REAL8TimeSeries *Tseries = mTseries->data[X];
      char buffer[LALNameLength];

      size_t written = 0;
      if ( 0 ) { // needed so that the 'else if' blocks can be removed if HAVE_LIBLALFRAME is not defined
#ifdef HAVE_LIBLALFRAME
      } else if ( XLALUserVarWasSet ( &uvar.outFrChannels ) ) { // if output frame channel names given, use those
        written = snprintf ( buffer, sizeof(buffer), "%s", uvar.outFrChannels->data[X] );
        if ( buffer[2] == ':' ) { // check we got correct IFO association
          XLAL_CHECK ( (buffer[0] == Tseries->name[0]) && (buffer[1] == Tseries->name[1]), XLAL_EINVAL,
                       "Possible IFO mismatch: outFrChannel[%d] = '%s', IFO = '%c%c': be careful about --outFrChannel ordering\n", X, buffer, Tseries->name[0], Tseries->name[1] );
        } // if buffer[2]==':'
      } else if ( XLALUserVarWasSet ( &uvar.inFrChannels ) ) { // otherwise: if input frame channel names given, use them for output, append "-<outLabel>"
        written = snprintf ( buffer, sizeof(buffer), "%s-%s", uvar.inFrChannels->data[X], uvar.outLabel );
#endif
      } else { // otherwise: fall back to <IFO>:<outLabel> channel name
        written = snprintf ( buffer, sizeof(buffer), "%c%c:%s", Tseries->name[0], Tseries->name[1], uvar.outLabel );
      }
      XLAL_CHECK ( written < LALNameLength, XLAL_ESIZE, "Output frame name exceeded max length (%d): '%s'\n", LALNameLength, buffer );

      outChannelNames->data[X] = XLALStringDuplicate( buffer );
    }

  if (uvar.outSFTdir)
    {
      XLAL_CHECK ( is_directory ( uvar.outSFTdir ), XLAL_EINVAL );

      /* generate comment string */
      CHAR *logstr;
      XLAL_CHECK ( (logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
      char *comment = XLALCalloc ( 1, len = strlen ( logstr ) + strlen(GV.VCSInfoString) + 512 );
      XLAL_CHECK ( comment != NULL, XLAL_ENOMEM, "XLALCalloc(1,%zu) failed.\n", len );
      sprintf ( comment, "Generated by:\n%s\n%s\n", logstr, GV.VCSInfoString );

      /* default window is rectangular */
      const char *window_type = ( GV.window_type != NULL ) ? GV.window_type : "rectangular";
      const REAL8 window_param = ( GV.window_type != NULL ) ? GV.window_param : 0;

      for ( UINT4 X=0; X < mSFTs->length; X ++ )
        {
          SFTVector *sfts = mSFTs->data[X];

          /* set up SFT filename */
          SFTFilenameSpec XLAL_INIT_DECL(spec);
          XLAL_CHECK ( XLALFillSFTFilenameSpecStrings( &spec, uvar.outSFTdir, NULL, NULL, window_type, NULL, NULL, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
          spec.window_param = window_param;

          if ( uvar.outPubObsRun > 0 ) { // public SFT filename

            /* get channel name and check for consistency */
            const char* channelName = outChannelNames->data[X];
            const SFTtype *sft0 = &sfts->data[0];
            XLAL_CHECK ( (channelName[0] == sft0->name[0]) && (channelName[1] == sft0->name[1]), XLAL_EINVAL,
                         "Possible IFO mismatch: output channel[%d] = '%s', IFO = '%c%c': be careful about output channel ordering\n", X, channelName, sft0->name[0], sft0->name[1] );

            /* set public observing run, revision, channel name */
            XLAL_CHECK ( XLALFillSFTFilenameSpecStrings( &spec, NULL, NULL, NULL, NULL, NULL, "SIM", channelName ) == XLAL_SUCCESS, XLAL_EFUNC );
            spec.pubObsRun = uvar.outPubObsRun;
            spec.pubRevision = uvar.outPubRevision;

          } else { // private SFT filename; see private description from --outLabel
            XLAL_CHECK ( XLALFillSFTFilenameSpecStrings( &spec, NULL, NULL, NULL, NULL, uvar.outLabel, NULL, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
          }

          /* either write whole SFT-vector to single concatenated file */
          XLAL_CHECK ( XLALWriteSFTVector2StandardFile( sfts, &spec, comment, uvar.outSingleSFT ) == XLAL_SUCCESS, XLAL_EFUNC );
        } // for X < numIFOs

      XLALFree ( logstr );
      XLALFree ( comment );

    } /* if outSFTdir */


   /* output ASCII time-series if requested */
  if ( uvar.TDDfile )
    {
      CHAR *fname = XLALCalloc (1, len = strlen(uvar.TDDfile) + 10 );
      XLAL_CHECK ( fname != NULL, XLAL_ENOMEM, "XLALCalloc(1,%zu) failed\n", len );

      for ( UINT4 X=0; X < mTseries->length; X ++ )
        {
          const REAL8TimeSeries *TS = mTseries->data[X];
          sprintf (fname, "%c%c-%s", TS->name[0], TS->name[1], uvar.TDDfile );
          XLAL_CHECK ( XLALdumpREAL8TimeSeries ( fname, TS ) == XLAL_SUCCESS, XLAL_EFUNC );

        } // for X < numDet

      XLALFree (fname);
    } /* if outputting ASCII time-series */

  /* output time-series to frames if requested */
#ifdef HAVE_LIBLALFRAME
  if ( GV.outFrameDir != NULL )
    {
      XLAL_CHECK ( is_directory ( uvar.outFrameDir ), XLAL_EINVAL );
      XLAL_CHECK ( XLALCheckValidDescriptionField ( uvar.outLabel ) == XLAL_SUCCESS, XLAL_EFUNC );
      len = strlen(GV.outFrameDir) + strlen(uvar.outLabel) + 100;
      char *fname;

      char *hist = XLALUserVarGetLog (UVAR_LOGFMT_CMDLINE);

      for ( UINT4 X=0; X < mTseries->length; X ++ )
        {
          REAL8TimeSeries *Tseries = mTseries->data[X];

          /* use standard frame output filename format */
          char IFO[2] = { Tseries->name[0], Tseries->name[1] };
          LIGOTimeGPS startTimeGPS = Tseries->epoch;
          REAL8 duration = Tseries->data->length * Tseries->deltaT;
          XLAL_CHECK ( (fname = LALCalloc (1, len )) != NULL, XLAL_ENOMEM );
          size_t written = snprintf ( fname, len, "%s/%c-%c%c_%s-%d-%d.gwf",
                                      GV.outFrameDir, IFO[0], IFO[0], IFO[1], uvar.outLabel, startTimeGPS.gpsSeconds, (int)duration );
          XLAL_CHECK ( written < len, XLAL_ESIZE, "Frame-filename exceeds expected maximal length (%zu): '%s'\n", len, fname );

          /* define the output frame */
          LALFrameH *outFrame;
          XLAL_CHECK ( (outFrame = XLALFrameNew ( &startTimeGPS, duration, uvar.outLabel, 1, 0, 0 )) != NULL, XLAL_EFUNC );

          /* add timeseries to the frame - make sure to change the timeseries name since this is used as the channel name */
          strcpy ( Tseries->name, outChannelNames->data[X] );
          XLAL_CHECK ( (XLALFrameAddREAL8TimeSeriesProcData ( outFrame, Tseries ) == XLAL_SUCCESS ) , XLAL_EFUNC );

          /* Here's where we add extra information into the frame - first we add the command line args used to generate it */
          XLALFrameAddFrHistory ( outFrame, __FILE__, hist );

          /* then we add the version string */
          XLALFrameAddFrHistory ( outFrame, __FILE__, GV.VCSInfoString );

          /* output the frame to file - compression level 1 (higher values make no difference) */
          XLAL_CHECK ( XLALFrameWrite ( outFrame, fname ) == XLAL_SUCCESS , XLAL_EFUNC );

          /* free the frame, frame file name and history memory */
          XLALFrameFree ( outFrame );
          XLALFree ( fname );

        } // for X < numDetectors

      XLALFree ( hist );

    } /* if GV.outFrameDir: outputting time-series to frames */
#endif // HAVE_LIBLALFRAME

  /* ---------- free memory ---------- */
  XLALDestroyMultiREAL8TimeSeries ( mTseries );
  XLALDestroyMultiSFTVector ( mSFTs );

  XLALDestroyStringVector( outChannelNames );

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

  /* if requested, log all user-input and code-versions */
  if ( uvar->logfile ) {
    XLAL_CHECK ( XLALWriteMFDlog ( uvar->logfile, cfg ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALWriteMFDlog() failed with xlalErrno = %d\n", xlalErrno );
  }

  /* Init ephemerides */
  XLAL_CHECK ( (cfg->edat = XLALInitBarycenter ( uvar->ephemEarth, uvar->ephemSun )) != NULL, XLAL_EFUNC );

  // ---------- check user-input consistency ----------

  // ----- check if frames + frame channels given
  BOOLEAN have_frames  = (uvar->inFrames != NULL);
  BOOLEAN have_channels= (uvar->inFrChannels != NULL);
  XLAL_CHECK ( !(have_frames || have_channels) || (have_frames && have_channels), XLAL_EINVAL, "Need both --inFrames and --inFrChannels, or NONE\n");

  // ----- IFOs : only from one of {--IFOs, --noiseSFTs, --inFrChannels}: mutually exclusive
  BOOLEAN have_IFOs      = (uvar->IFOs != NULL);
  BOOLEAN have_noiseSFTs = (uvar->noiseSFTs != NULL);
  XLAL_CHECK ( have_frames || have_IFOs || have_noiseSFTs, XLAL_EINVAL, "Need one of --IFOs, --noiseSFTs or --inFrChannels to determine detectors\n");

  if ( have_frames ) {
    XLAL_CHECK ( !have_IFOs && !have_noiseSFTs, XLAL_EINVAL, "If --inFrames given, cannot handle --IFOs or --noiseSFTs input\n");
    XLAL_CHECK ( XLALParseMultiLALDetector ( &(cfg->multiIFO), uvar->inFrChannels ) == XLAL_SUCCESS, XLAL_EFUNC );
  } else { // !have_frames
    XLAL_CHECK ( !(have_IFOs && have_noiseSFTs), XLAL_EINVAL, "Cannot handle both --IFOs and --noiseSFTs input\n");
  }
  if ( have_IFOs ) {
    XLAL_CHECK ( XLALParseMultiLALDetector ( &(cfg->multiIFO), uvar->IFOs ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // ----- TIMESTAMPS: either from --timestampsFiles, --startTime+duration, or --noiseSFTs
  BOOLEAN have_startTime = XLALUserVarWasSet ( &uvar->startTime );
  BOOLEAN have_duration = XLALUserVarWasSet ( &uvar->duration );
  BOOLEAN have_timestampsFiles = ( uvar->timestampsFiles != NULL );
  // need BOTH startTime+duration or none
  XLAL_CHECK ( ( have_duration && have_startTime) || !( have_duration || have_startTime ), XLAL_EINVAL, "Need BOTH {--startTime,--duration} or NONE\n");
  // at least one of {startTime,timestamps,noiseSFTs,inFrames} required
  XLAL_CHECK ( have_timestampsFiles || have_startTime || have_noiseSFTs || have_frames, XLAL_EINVAL, "Need at least one of {--timestampsFiles, --startTime+duration, --noiseSFTs, --inFrames}\n" );
  // don't allow timestamps + {startTime+duration OR noiseSFTs}
  XLAL_CHECK ( !have_timestampsFiles || !(have_startTime||have_noiseSFTs), XLAL_EINVAL, "--timestampsFiles incompatible with {--noiseSFTs or --startTime+duration}\n");
  // note, however, that we DO allow --noiseSFTs and --startTime+duration, which will act as a constraint
  // on the noise-SFTs to load in

  // don't allow --SFToverlap with either --noiseSFTs OR --timestampsFiles
  XLAL_CHECK ( uvar->SFToverlap >= 0, XLAL_EDOM );
  BOOLEAN haveOverlap = ( uvar->SFToverlap > 0 );
  XLAL_CHECK ( !haveOverlap || !( have_noiseSFTs || have_timestampsFiles ), XLAL_EINVAL, "--SFToverlap incompatible with {--noiseSFTs or --timestampsFiles}\n" );


  // frequency-band either taken here from user-input, or later from noiseSFTs (no defaults)
  UINT4 nSetfminBand = UVAR_SET2(fmin,Band);
  XLAL_CHECK ( (nSetfminBand==2) || ( (nSetfminBand==0) && have_noiseSFTs ), XLAL_EINVAL, "Need either 'noiseSFTs' or BOTH of 'fmin' and 'Band'!\n");
  if ( nSetfminBand==2 )
    {
      /* check for negative fmin and Band, which would break the fMin_eff, fBand_eff calculation below */
      XLAL_CHECK ( uvar->fmin >= 0, XLAL_EDOM, "Invalid negative frequency fmin=%f!\n\n", uvar->fmin );
      XLAL_CHECK ( uvar->Band > 0, XLAL_EDOM, "Invalid non-positive frequency band Band=%f!\n\n", uvar->Band );
      cfg->fminOut = uvar->fmin;
      cfg->BandOut = uvar->Band;
    }

  // now handle the 3 mutually-exclusive cases: have_noiseSFTs || have_timestampsFiles || have_startTime (only)
  if ( have_noiseSFTs )
    {
      SFTConstraints XLAL_INIT_DECL(constraints);
      if ( have_startTime && have_duration )	 // use optional (startTime+duration) as constraints,
        {
          LIGOTimeGPS minStartTime, maxStartTime;
          minStartTime = uvar->startTime;
          maxStartTime = uvar->startTime;
          XLALGPSAdd ( &maxStartTime, uvar->duration );
          constraints.minStartTime = &minStartTime;
          constraints.maxStartTime = &maxStartTime;
          char bufGPS1[32], bufGPS2[32];
          XLALPrintWarning ( "Only noise-SFTs between GPS [%s, %s] will be used!\n", XLALGPSToStr(bufGPS1, &minStartTime), XLALGPSToStr(bufGPS2, &maxStartTime) );
        } /* if start+duration given */
      XLAL_CHECK ( (cfg->noiseCatalog = XLALSFTdataFind ( uvar->noiseSFTs, &constraints )) != NULL, XLAL_EFUNC );
      XLAL_CHECK (  cfg->noiseCatalog->length > 0, XLAL_EINVAL, "No noise-SFTs matching (start+duration, timestamps) were found!\n" );
      XLAL_CHECK ( (cfg->multiNoiseCatalogView = XLALGetMultiSFTCatalogView ( cfg->noiseCatalog )) != NULL, XLAL_EFUNC );

      // extract multi-timestamps from the multi-SFT-catalog view
      XLAL_CHECK ( (cfg->multiTimestamps = XLALTimestampsFromMultiSFTCatalogView ( cfg->multiNoiseCatalogView )) != NULL, XLAL_EFUNC );
      // extract IFOs from multi-SFT catalog
      XLAL_CHECK ( XLALMultiLALDetectorFromMultiSFTCatalogView ( &(cfg->multiIFO), cfg->multiNoiseCatalogView ) == XLAL_SUCCESS, XLAL_EFUNC );

      // if not given by user: extract frequency start + band from noise SFT catalog
      if ( !UVAR_ANYSET2(fmin,Band) )
        {
          const SFTDescriptor *desc = &(cfg->multiNoiseCatalogView->data[0].data[0]);
          REAL8 noise_fmin    = desc->header.f0;
          REAL8 noise_dFreq   = desc->header.deltaF;
          UINT4 noise_numBins = desc->numBins;
          REAL8 noise_band    = noise_numBins * noise_dFreq;
          cfg->fminOut = noise_fmin;
          cfg->BandOut = noise_band;
        }

    } // endif have_noiseSFTs
  else if ( have_timestampsFiles )
    {
      XLAL_CHECK ( (cfg->multiTimestamps = XLALReadMultiTimestampsFiles ( uvar->timestampsFiles )) != NULL, XLAL_EFUNC );

      XLAL_CHECK ( (cfg->multiTimestamps->length > 0) && (cfg->multiTimestamps->data != NULL), XLAL_EINVAL, "Got empty timestamps-list from XLALReadMultiTimestampsFiles()\n" );

      for ( UINT4 X=0; X < cfg->multiTimestamps->length; X ++ ) {
        XLAL_CHECK ( (cfg->multiTimestamps->data[X]->length > 0) && (cfg->multiTimestamps->data[X]->data != NULL), XLAL_EINVAL, "Got empty timestamps-list for detector %s", cfg->multiIFO.sites[X].frDetector.prefix );
        cfg->multiTimestamps->data[X]->deltaT = uvar->Tsft;	// Tsft information not given by timestamps-file
      }
    } // endif have_timestampsFiles
  else if ( have_startTime && have_duration )
    {
      XLAL_CHECK ( ( cfg->multiTimestamps = XLALMakeMultiTimestamps ( uvar->startTime, uvar->duration, uvar->Tsft, uvar->SFToverlap, cfg->multiIFO.length )) != NULL, XLAL_EFUNC );
    } // endif have_startTime

  // check if the user asked for Gaussian white noise to be produced (sqrtSn[X]!=0), otherwise leave noise-floors at 0
  if ( uvar->sqrtSX != NULL ) {
    XLAL_CHECK ( XLALParseMultiNoiseFloor ( &(cfg->multiNoiseFloor), uvar->sqrtSX, cfg->multiIFO.length ) == XLAL_SUCCESS, XLAL_EFUNC );
  } else {
    cfg->multiNoiseFloor.length = cfg->multiIFO.length;
    // values remain at their default sqrtSn[X] = 0;
  }

  // determine SFT window to use
  {
    const BOOLEAN have_window = XLALUserVarWasSet( &uvar->SFTWindowType );
    const BOOLEAN have_param = XLALUserVarWasSet( &uvar->SFTWindowParam );
    const CHAR* window_type_from_uvar = uvar->SFTWindowType;
    const REAL8 window_param_from_uvar = have_param ? uvar->SFTWindowParam : 0;
    if ( have_window )
      {
        XLAL_CHECK ( XLALCheckNamedWindow ( window_type_from_uvar, have_param ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    if ( have_noiseSFTs )
      {
        const CHAR* window_type_from_noiseSFTs = cfg->noiseCatalog->data[0].window_type;
        const REAL8 window_param_from_noiseSFTs = cfg->noiseCatalog->data[0].window_param;
        const BOOLEAN have_window_type_from_noiseSFTs = ( XLALStringCaseCompare( window_type_from_noiseSFTs, "unknown" ) != 0 );
        if ( have_window_type_from_noiseSFTs ^ have_window )
          {
            cfg->window_type = have_window_type_from_noiseSFTs ? window_type_from_noiseSFTs : window_type_from_uvar;
            cfg->window_param = have_window_type_from_noiseSFTs ? window_param_from_noiseSFTs : window_param_from_uvar;
          }
        else
          {
            /* EITHER noise SFTs must have a known window OR user must specify the window function */
            XLAL_ERROR ( XLAL_EINVAL, "When --noiseSFTs is given, --SFTWindowType is required ONLY if noise SFTs have unknown window.\n" );
          }
      }
    else if ( have_window )
      {
        cfg->window_type = window_type_from_uvar;
        cfg->window_param = window_param_from_uvar;
      }
  }

#ifdef HAVE_LIBLALFRAME
  // if user requested time-series data from frames to be added: try to read the frames now
  if ( have_frames )
    {
      UINT4 numDetectors = uvar->inFrChannels->length;
      XLAL_CHECK ( uvar->inFrames->length == numDetectors, XLAL_EINVAL, "Need equal number of channel names (%d) as frame specifications (%d)\n", uvar->inFrChannels->length, numDetectors );

      XLAL_CHECK ( (cfg->inputMultiTS = XLALCalloc ( 1, sizeof(*cfg->inputMultiTS))) != NULL, XLAL_ENOMEM );
      cfg->inputMultiTS->length = numDetectors;
      XLAL_CHECK ( (cfg->inputMultiTS->data = XLALCalloc ( numDetectors, sizeof(cfg->inputMultiTS->data[0]) )) != NULL, XLAL_ENOMEM );
      if ( cfg->multiTimestamps == NULL )
        {
          XLAL_CHECK ( (cfg->multiTimestamps = XLALCalloc ( 1, sizeof(*cfg->multiTimestamps) )) != NULL, XLAL_ENOMEM );
          XLAL_CHECK ( (cfg->multiTimestamps->data = XLALCalloc ( numDetectors, sizeof(cfg->multiTimestamps->data[0]))) != NULL, XLAL_ENOMEM );
          cfg->multiTimestamps->length = numDetectors;
        }
      for ( UINT4 X = 0; X < numDetectors; X ++ )
        {
          LALCache *cache;
          XLAL_CHECK ( (cache = XLALCacheImport ( uvar->inFrames->data[X] )) != NULL, XLAL_EFUNC, "Failed to import cache file '%s'\n", uvar->inFrames->data[X] );
          // ----- open frame stream ----------
          LALFrStream *stream;
          XLAL_CHECK ( (stream = XLALFrStreamCacheOpen ( cache )) != NULL, XLAL_EFUNC, "Failed to open stream from cache file '%s'\n", uvar->inFrames->data[X] );
          XLALDestroyCache ( cache );

          // ----- determine time span of frame data from the stream ----------
          // determine time spanned by the frames in this cache
          LIGOTimeGPS frames_startGPS, frames_endGPS;
          XLAL_CHECK ( XLALFrStreamSeekO ( stream, 0, SEEK_SET ) == 0, XLAL_EFUNC, "Failed to move to start of input frame-stream\n");
          XLAL_CHECK ( XLALFrStreamTell ( &frames_startGPS, stream ) == XLAL_SUCCESS, XLAL_EFUNC );
          XLAL_CHECK ( XLALFrStreamSeekO ( stream, 0, SEEK_END ) == 0, XLAL_EFUNC, "Failed to move to end of input frame-stream\n");
          XLAL_CHECK ( XLALFrStreamTell ( &frames_endGPS, stream ) == XLAL_SUCCESS, XLAL_EFUNC );
          REAL8 frames_start = XLALGPSGetREAL8 ( &frames_startGPS );
          REAL8 frames_end   = XLALGPSGetREAL8 ( &frames_endGPS );
          REAL8 frames_span = frames_end - frames_start;

          LIGOTimeGPS ts_startGPS;
          REAL8 ts_duration;
          // check that it's consistent with timestamps, if given, otherwise create timestamps from this
          if ( cfg->multiTimestamps->data[X] != NULL )	// FIXME: implicitly assumes timestamps are sorted, which is not guaranteed by timestamps-reading from file
            {
              const LIGOTimeGPSVector *timestampsX = cfg->multiTimestamps->data[X];
              REAL8 tStart = XLALGPSGetREAL8( &timestampsX->data[0] );
              REAL8 tEnd   = XLALGPSGetREAL8( &timestampsX->data[timestampsX->length-1]) + timestampsX->deltaT;
              XLAL_CHECK ( tStart >= frames_start && tEnd <= frames_end, XLAL_EINVAL, "Detector X=%d: Requested timestamps-range [%.0f, %.0f]s outside of cache range [%.0f,%.0f]s\n",
                           X, tStart, tEnd, frames_start, frames_end );
              XLALGPSSetREAL8 ( &ts_startGPS, tStart );
              ts_duration = (tEnd - tStart);
            } // if have_timestamps
          else
            {
              ts_startGPS = frames_startGPS;
              ts_duration = frames_span;
              XLAL_CHECK ( (cfg->multiTimestamps->data[X] = XLALMakeTimestamps ( ts_startGPS, ts_duration, uvar->Tsft, uvar->SFToverlap ) ) != NULL, XLAL_EFUNC );
              // in this case we can't allow timestamps extending beyond the end of the frame timeseries, so we drop the last timestamp if necessary
              LIGOTimeGPSVector *timestampsX = cfg->multiTimestamps->data[X];
              REAL8 tEnd;
              while ( (tEnd = XLALGPSGetREAL8 ( &timestampsX->data[timestampsX->length-1] ) + timestampsX->deltaT) > frames_end ) {
                timestampsX->length --;
              }
            } // if have no timestamps

          const char *channel = uvar->inFrChannels->data[X];
          size_t limit = 0;	// unlimited read
          REAL8TimeSeries *ts;
          XLAL_CHECK ( (ts = XLALFrStreamInputREAL8TimeSeries ( stream, channel, &ts_startGPS, ts_duration, limit )) != NULL,
                       XLAL_EFUNC, "Frame reading failed for stream created for '%s': ts_start = {%d,%d}, duration=%.0f\n",
                       uvar->inFrames->data[X], ts_startGPS.gpsSeconds, ts_startGPS.gpsNanoSeconds, ts_duration );

          if ( XLALUnitCompare ( &(ts->sampleUnits), &lalStrainUnit ) != 0 )
            {
              if ( XLALUnitCompare ( &(ts->sampleUnits), &lalADCCountUnit ) == 0 )
                {
                  XLALPrintWarning ("Note: Input timeseries has dimensions 'ADC count' [count] instead of strain.\n");
                }
              else
                {
                  char unitName[64];
                  XLAL_CHECK ( XLALUnitAsString( unitName, sizeof(unitName), &(ts->sampleUnits) ) != NULL, XLAL_EUNIT,
                               "Failed to express time-series units as a string\n");
                  XLAL_CHECK ( XLALUnitIsDimensionless ( &(ts->sampleUnits) ), XLAL_EUNIT,
                               "Timeseries input data has invalid units '%s': must be either strain or dimensionless!", unitName );
                  XLALPrintWarning ("Note: Input timeseries is dimensionless instead of strain units.\n");
                }
              ts->sampleUnits = lalStrainUnit;	// needed to make XLALCWMakeFakeData() work
            }

          cfg->inputMultiTS->data[X] = ts;

          XLAL_CHECK ( XLALFrStreamClose ( stream ) == XLAL_SUCCESS, XLAL_EFUNC, "Stream closing failed for cache file '%s'\n", uvar->inFrames->data[X] );
        } // for X < numDetectors
    } // if inFrames

  if ( uvar->outFrameDir ) {
    cfg->outFrameDir = uvar->outFrameDir;
  }
#endif

  return XLAL_SUCCESS;

} /* XLALInitMakefakedata() */


/**
 * Register all "user-variables", and initialized them from command-line and config-files
 */
int
XLALInitUserVars ( UserVariables_t *uvar, int argc, char *argv[] )
{

  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL, "Invalid NULL input 'uvar'\n");
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL, "Invalid NULL input 'argv'\n");

  // ---------- set a few defaults ----------
  uvar->ephemEarth = XLALStringDuplicate("earth00-40-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-40-DE405.dat.gz");

  uvar->Tsft = 1800;
  uvar->outSingleSFT = 1; /* write our a single SFT file by default */

#define MISC_DEFAULT "mfdv5"
  XLAL_CHECK ( (uvar->outLabel = XLALStringDuplicate ( MISC_DEFAULT ))  != NULL, XLAL_EFUNC );

  // ---------- register all our user-variable ----------
  /* output options */
  XLALRegisterUvarMember(   outSingleSFT,       BOOLEAN, 's', OPTIONAL, "Write a single concatenated SFT file instead of individual files" );
  XLALRegisterUvarMember( outSFTdir,          STRING, 'n', OPTIONAL, "Output SFTs:  directory for output SFTs");
  XLALRegisterUvarMember(  outLabel,	         STRING, 0, OPTIONAL, "'misc' entry in SFT-filenames or 'description' entry of frame filenames" );
  XLALRegisterUvarMember( outPubObsRun,       UINT4, 'O', OPTIONAL, "if >0, names SFTs using the public filename convention with this observing run number" );
  XLALRegisterUvarMember( outPubRevision,     UINT4, 'R', OPTIONAL, "if " UVAR_STR(outPubObsRun) ">0, names SFTs using the public filename convention with this revision number" );
  XLALRegisterUvarMember( TDDfile,            STRING, 't', OPTIONAL, "Filename to output time-series into");

  XLALRegisterUvarMember( logfile,            STRING, 'l', OPTIONAL, "Filename for log-output");

  /* detectors and respective noise-floors */
  XLALRegisterUvarMember( IFOs,			STRINGVector, 'I', OPTIONAL, "CSV list of detectors, eg. \"H1,H2,L1,G1, ...\" ");
  XLALRegisterUvarMember( sqrtSX,	 	 STRINGVector, 0,  OPTIONAL, "Add Gaussian Noise: CSV list of detectors' noise-floors sqrt{Sn}");

  XLALRegisterUvarMember( ephemEarth, 	 	STRING, 0,  OPTIONAL, "Earth ephemeris file to use");
  XLALRegisterUvarMember( ephemSun, 	 	STRING, 0,  OPTIONAL, "Sun ephemeris file to use");

  /* start + duration of timeseries */
  XLALRegisterUvarMember(startTime,            EPOCH, 'G', OPTIONAL, "Start-time of requested signal in detector-frame (format 'xx.yy[GPS|MJD]')");
  XLALRegisterUvarMember(  duration,              INT4, 0,  OPTIONAL, "Duration of requested signal in seconds");
  XLALRegisterUvarMember( timestampsFiles,       STRINGVector, 0,  OPTIONAL, "ALTERNATIVE: File to read timestamps from (file-format: lines with <seconds> <nanoseconds>)");

  /* sampling and heterodyning frequencies */
  XLALRegisterUvarMember(  fmin,                 REAL8, 0, NODEFAULT, "Lowest frequency (Hz) of output SFT, and heterodyning frequency of time series written to frames; REQUIRED unless --noiseSFTs given.\n\nNote that, since the output time series written to frames are heterodyned at this frequency, signals at frequency 'f' will be shifted to frequency 'f - fmin'.");
  XLALRegisterUvarMember(  Band,                 REAL8, 0, NODEFAULT, "Bandwidth (Hz) of output SFT, and of time series written to frames (i.e. half the time series sampling frequency); REQUIRED unless --noiseSFTs given.");

  /* SFT properties */
  XLALRegisterUvarMember(  Tsft,                 REAL8, 0, OPTIONAL, "Time baseline of one SFT in seconds");
  XLALRegisterUvarMember(  SFToverlap,           REAL8, 0, OPTIONAL, "Overlap between successive SFTs in seconds (conflicts with --noiseSFTs or --timestampsFiles)");
  XLALRegisterUvarMember( SFTWindowType,        STRING, 0, OPTIONAL, "Window function to apply to the SFTs ('rectangular', 'hann', 'tukey', etc.); when --noiseSFTs is given, required ONLY if noise SFTs have unknown window");
  XLALRegisterUvarMember( SFTWindowParam,        REAL8, 0, OPTIONAL, "Window parameter required for a few window-types (eg. 'tukey')");

  /* pulsar params */
  XLALRegisterUvarMember( injectionSources,     STRINGVector, 0, OPTIONAL, "%s", InjectionSourcesHelpString );
  XLALRegisterUvarMember( sourceDeltaT,         REAL8,  0, OPTIONAL, "Source-frame sampling period. '0' implies implies defaults set in XLALGeneratePulsarSignal()." );

  /* noise */
  XLALRegisterUvarMember( noiseSFTs,          STRING, 'D', OPTIONAL, "Noise-SFTs to be added to signal. Possibilities are:\n"
                          " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'\n"
                          "(Used also to set IFOs and timestamps, and frequency range unless separately specified.)");
  XLALRegisterUvarMember( randSeed,           INT4, 0, OPTIONAL, "Specify random-number seed for reproducible noise (0 means use /dev/urandom for seeding).");

  /* frame input/output options */
#ifdef HAVE_LIBLALFRAME
  XLALRegisterUvarMember ( outFrameDir,	  STRING,      'F', OPTIONAL,  "Output Frames: directory for output timeseries frame files");
  XLALRegisterUvarMember ( inFrames, 	  STRINGVector,'C', OPTIONAL,  "CSV list (one per IFO) of input frame cache files");
  XLALRegisterUvarMember ( inFrChannels,  STRINGVector,'N', OPTIONAL,  "CSV list (one per IFO) of frame channels to read timeseries from");
  XLALRegisterUvarMember ( outFrChannels, STRINGVector, 0,  NODEFAULT, "CSV list (one per IFO) of output frame channel names [default: <inFrChannels>-<outLabel> or <IFO>:<outLabel>]");
#else
  XLALRegisterUvarMember ( outFrameDir,	 STRING,       'F', DEPRECATED, "Need to compile with lalframe support for this option to work");
  XLALRegisterUvarMember ( inFrames,     STRINGVector, 'C', DEPRECATED, "Need to compile with lalframe support for this option to work");
  XLALRegisterUvarMember ( inFrChannels, STRINGVector, 'N', DEPRECATED, "Need to compile with lalframe support for this option to work");
  XLALRegisterUvarMember ( outFrChannels, STRINGVector, 0,  DEPRECATED, "Need to compile with lalframe support for this option to work");
#endif

  // ----- 'expert-user/developer' options ----- (only shown in help at lalDebugLevel >= warning)

  /* read cmdline & cfgfile  */
  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
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

  XLALFree ( cfg->VCSInfoString );

  // free noise-SFT catalog
  XLALDestroySFTCatalog ( cfg->noiseCatalog );
  XLALDestroyMultiSFTCatalogView ( cfg->multiNoiseCatalogView );

  // free noise time-series data
  XLALDestroyMultiREAL8TimeSeries ( cfg->inputMultiTS );

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
