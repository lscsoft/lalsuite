/*
 * Copyright (C) 2010 Reinhard Prix, 2013 David Keitel
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
 * \author David Keitel, Reinhard Prix
 * \brief
 * Simple standalone code to provide ASCII output for antenna-pattern function and matrix
 * for given detector(s) and sky-location
 * based on PrintDetectorState
 *
 */

/* ---------- includes ---------- */
#include <math.h>
#include <errno.h>
#include <string.h>

#include <lal/UserInput.h>
#include <lal/LALInitBarycenter.h>
#include <lal/ComputeFstat.h>
#include <lal/LALString.h>
#include <lal/StringVector.h>

#include <lalapps.h>

/* ----- compile switches ----- */

/*---------- local defines ---------- */

/* ----- Macros ----- */
#define SQ(x) ( (x) * (x) )

/* ---------- local types ---------- */

typedef struct
{
  EphemerisData *edat;			/**< ephemeris data (from LALInitBarycenter()) */
  UINT4 numDetectors;			/**< number of detectors */
  MultiDetectorStateSeries *multiDetStates;	/**< detector state time series */
  MultiNoiseWeights *multiNoiseWeights;		/**< per-detector noise weights */
  MultiLIGOTimeGPSVector *multiTimestamps;	/**< timestamps vector (LIGOtimeGPS format) */
  UINT4 numTimeStamps;			/**< number of timestamps (SFTs) for all detectors */
  UINT4Vector *numTimeStampsX;		/**< number of timestamps (SFTs) per detector */
  REAL8Vector *Alpha;			/**< skyposition Alpha: radians, equatorial coords */
  REAL8Vector *Delta;			/**< skyposition Delta: radians, equatorial coords */
  UINT4 numSkyPoints;			/**< common length of Alpha and Delta vectors */
} ConfigVariables;


typedef struct
{
  BOOLEAN help;

  LALStringVector* IFOs; /**< list of detector-names "H1,H2,L1,.." or single detector*/

  REAL8 Alpha;		/**< a single skyposition Alpha: radians, equatorial coords. */
  REAL8 Delta;		/**< a single skyposition Delta: radians, equatorial coords. */
  CHAR *skyGridFile;	/**< alternative: matrix of (Alpha,Delta) pairs from a file */

  CHAR *ephemEarth;	/**< Earth ephemeris file to use */
  CHAR *ephemSun;	/**< Sun ephemeris file to use */

  LALStringVector* timeGPS;	/**< GPS timestamps to compute detector state for (REAL8 format) */
  CHAR *timeStampsFile;		/**< alternative: read in timestamps from a single file (deprecated) */
  LALStringVector *timeStampsFiles;		/**< alternative: read in timestamps from file(s) */
  INT4 Tsft;			/**< assumed length of SFTs, needed for offset to timestamps when comparing to CFS_v2, PFS etc */

  LALStringVector* noiseSqrtShX; /**< per-detector noise PSD sqrt(SX) */

  CHAR *outab; 			/**< output file for antenna pattern functions a(t), b(t) at each timestamp */
  CHAR *outABCD; 		/**< output file for antenna pattern matrix elements A, B, C, D averaged over timestamps */

  BOOLEAN version;	/**< output code versions */

} UserVariables_t;

/* ---------- global variables ----------*/
extern int vrbflg;

/* ---------- local prototypes ---------- */
int XLALInitUserVars ( UserVariables_t *uvar );
int XLALInitCode ( ConfigVariables *cfg, const UserVariables_t *uvar, const char *app_name);
int XLALDestroyConfig ( ConfigVariables *cfg );

/*============================================================
 * FUNCTION definitions
 *============================================================*/

int
main(int argc, char *argv[])
{

  ConfigVariables XLAL_INIT_DECL(config);
  UserVariables_t XLAL_INIT_DECL(uvar);

  /* register user-variables */

  XLAL_CHECK ( XLALInitUserVars ( &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* read cmdline & cfgfile  */
  XLAL_CHECK ( XLALUserVarReadAllInput ( argc,argv ) == XLAL_SUCCESS, XLAL_EFUNC );

  if (uvar.help) { 	/* help requested: we're done */
    exit(0);
  }

  if ( uvar.version )
    {
      XLALOutputVersionString ( stdout, lalDebugLevel );
      exit(0);
    }

  /* basic setup and initializations */
  XLAL_CHECK ( XLALInitCode( &config, &uvar, argv[0] ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* prepare output files */
  FILE *fpOutab = NULL;
  if ( uvar.outab ) {

      XLAL_CHECK ( (fpOutab = fopen (uvar.outab, "wb")) != NULL, XLAL_EIO, "Error opening file '%s' for writing...", uvar.outab );

      /* write header info in comments */
      XLAL_CHECK ( XLAL_SUCCESS == XLALOutputVersionString ( fpOutab, 0 ), XLAL_EFUNC );

      /* write the command-line */
      for (int a = 0; a < argc; a++)
        {
          fprintf(fpOutab,"%%%% argv[%d]: '%s'\n", a, argv[a]);
        }

      /* write column headings */
      fprintf(fpOutab, "%%%% columns:\n%%%% Alpha   Delta       tGPS ");
      if ( config.numDetectors == 1 ) {
        fprintf(fpOutab, "      a(t)         b(t)");
      }
      else {
        for ( UINT4 X=0; X < config.numDetectors; X++ ) {
          fprintf(fpOutab, "      a[%d](t)      b[%d](t)", X, X);
        }
      }
      fprintf(fpOutab, "\n");

  }

  FILE *fpOutABCD = NULL;
  if ( uvar.outABCD ) {

      XLAL_CHECK ( (fpOutABCD = fopen (uvar.outABCD, "wb")) != NULL, XLAL_EIO, "Error opening file '%s' for writing...", uvar.outABCD );

      /* write header info in comments */
      XLAL_CHECK ( XLAL_SUCCESS == XLALOutputVersionString ( fpOutABCD, 0 ), XLAL_EFUNC );

      /* write the command-line */
      for (int a = 0; a < argc; a++)
        {
          fprintf(fpOutABCD,"%%%% argv[%d]: '%s'\n", a, argv[a]);
        }

      /* write column headings */
      fprintf(fpOutABCD, "%%%% columns:\n%%%% Alpha   Delta");
      fprintf(fpOutABCD, "        A            B            C            D");
      if ( config.numDetectors > 1 ) {
        fprintf(fpOutABCD, "   ");
        for ( UINT4 X=0; X < config.numDetectors; X++ ) {
          fprintf(fpOutABCD, "         A[%d]         B[%d]         C[%d]         D[%d]", X, X, X, X);
        }
      }
      fprintf(fpOutABCD, "\n");

  }

  /* loop over sky positions (outer loop, could allow for buffering if necessary) */
  for (UINT4 n = 0; n < config.numSkyPoints; n++) {
    SkyPosition skypos;
    skypos.system = COORDINATESYSTEM_EQUATORIAL;
    skypos.longitude = config.Alpha->data[n];
    skypos.latitude  = config.Delta->data[n];

    /* do the actual computation of the antenna pattern functions */
    MultiAMCoeffs *multiAM;
    XLAL_CHECK ( ( multiAM = XLALComputeMultiAMCoeffs ( config.multiDetStates, config.multiNoiseWeights, skypos ) ) != NULL, XLAL_EFUNC, "XLALComputeAMCoeffs() failed." );

    /* for multi-IFO run with weights, do it again, without weights, to get single-IFO quantities consistent with single-IFO runs
     * FIXME: remove this temporary hack when MultiAmCoeffs have been changed to include non-weighted single-IFO quantities
     */
    MultiAMCoeffs *multiAMforSingle = NULL;
    MultiAMCoeffs *multiAMunweighted = NULL;
    if ( ( config.numDetectors > 1 ) && ( config.multiNoiseWeights != NULL ) ) {
      XLAL_CHECK ( ( multiAMunweighted = XLALComputeMultiAMCoeffs ( config.multiDetStates, NULL, skypos ) ) != NULL, XLAL_EFUNC, "XLALComputeAMCoeffs() failed." );
      multiAMforSingle = multiAMunweighted;
    }
    else {
      multiAMforSingle = multiAM;
    }

    /* write out the data for this sky point */
    if ( uvar.outab ) { // output a(t), b(t) at each timestamp
      for (UINT4 t = 0; t < config.numTimeStampsX->data[0]; t++) { // FIXME: does not work for different multi-IFO numTimeStampsX
         fprintf (fpOutab, "%.7f  %.7f  %d", config.Alpha->data[n], config.Delta->data[n], config.multiTimestamps->data[0]->data[t].gpsSeconds );
         for ( UINT4 X=0; X < config.numDetectors; X++ ) {
           fprintf(fpOutab, " %12.8f %12.8f", multiAMforSingle->data[X]->a->data[t], multiAMforSingle->data[X]->b->data[t]);
         } // for ( UINT4 X=0; X < config.numDetectors; X++ )
         fprintf(fpOutab, "\n");
       } // for (UINT4 t = 0; t < config.numTimeStamps; t++)
    } // if ( uvar.outab )

    if ( uvar.outABCD ) { // output ABCD averaged over all timestamps
      // FIXME: stop doing average manually when AMCoeffs is changed to contain averaged values
      REAL8 A = multiAM->Mmunu.Ad/config.numTimeStamps;
      REAL8 B = multiAM->Mmunu.Bd/config.numTimeStamps;
      REAL8 C = multiAM->Mmunu.Cd/config.numTimeStamps;
      REAL8 D = A*B-SQ(C);
      fprintf (fpOutABCD, "%.7f  %.7f %12.8f %12.8f %12.8f %12.8f", config.Alpha->data[n], config.Delta->data[n], A, B, C, D );
      if ( config.numDetectors > 1 ) {
        for ( UINT4 X=0; X < config.numDetectors; X++ ) {
          REAL4 AX = multiAMforSingle->data[X]->A/config.numTimeStampsX->data[X];
          REAL4 BX = multiAMforSingle->data[X]->B/config.numTimeStampsX->data[X];
          REAL4 CX = multiAMforSingle->data[X]->C/config.numTimeStampsX->data[X];
          REAL4 DX = AX*BX-SQ(CX);
          fprintf(fpOutABCD, " %12.8f %12.8f %12.8f %12.8f", AX, BX, CX, DX);
        }
      }
      fprintf(fpOutABCD, "\n");
    } // if ( uvar.outABCD )

    XLALDestroyMultiAMCoeffs ( multiAM );
    if ( multiAMunweighted ) {
      XLALDestroyMultiAMCoeffs ( multiAMunweighted );
    }

  } // for (UINT4 n = 0; n < config.numSkyPoints; n++)

  /* ----- close output files ----- */
  if ( fpOutab ) {
    fprintf (fpOutab, "\n");
    fclose ( fpOutab );
  }
  if ( fpOutABCD ) {
    fprintf (fpOutABCD, "\n");
    fclose ( fpOutABCD );
  }

  /* ----- done: free all memory */
  XLAL_CHECK ( XLALDestroyConfig( &config ) == XLAL_SUCCESS, XLAL_EFUNC );

  LALCheckMemoryLeaks();

  return 0;
} /* main */


/** register all "user-variables" */
int
XLALInitUserVars ( UserVariables_t *uvar )
{
  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

  /* set a few defaults */
  uvar->help = 0;

  XLAL_CHECK ( (uvar->IFOs = XLALCreateStringVector ( "H1", NULL )) != NULL, XLAL_ENOMEM, "Call to XLALCreateStringVector() failed." );

  uvar->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");

  uvar->Alpha     = 0.0;
  uvar->Delta     = 0.0;
  uvar->skyGridFile = NULL;

  uvar->timeGPS = NULL;
  uvar->timeStampsFile = NULL;
  uvar->outab = 0;
  uvar->outABCD = 0;
  uvar->Tsft = 1800;

  uvar->noiseSqrtShX = NULL;

  /* register all user-variables */
  XLALregBOOLUserStruct(	help,		'h', UVAR_HELP,		"Print this help/usage message");
  XLALregLISTUserStruct( IFOs,                  'I', UVAR_OPTIONAL, "Comma-separated list of detectors, eg. \"H1,H2,L1,G1, ...\" [only 1 detector supported at the moment] ");

  XLALregREALUserStruct(	Alpha,		'a', UVAR_OPTIONAL,	"single skyposition Alpha in radians, equatorial coords.");
  XLALregREALUserStruct(	Delta, 		'd', UVAR_OPTIONAL,	"single skyposition Delta in radians, equatorial coords.");

  XLALregSTRINGUserStruct( skyGridFile,		's', UVAR_OPTIONAL,	"Alternatively: sky-grid file");

  XLALregLISTUserStruct( 	timeGPS,        't', UVAR_OPTIONAL, 	"GPS time at which to compute detector states (separate multiple timestamps by commata)");
  XLALregLISTUserStruct(	timeStampsFiles, 'T', UVAR_OPTIONAL,	"Alternative: time-stamps file(s) (comma-separated list per IFO, or one for all)");
  XLALregINTUserStruct(		Tsft,		 0, UVAR_OPTIONAL,	"Assumed length of one SFT in seconds; needed for timestamps offset consistency with F-stat based codes");

  XLALregLISTUserStruct ( noiseSqrtShX,		 0, UVAR_OPTIONAL, "Per-detector noise PSD sqrt(SX). Only ratios relevant to compute noise weights. Defaults to 1,1,...");

  XLALregSTRINGUserStruct (	ephemEarth,	 0,  UVAR_OPTIONAL,	"Earth ephemeris file to use");
  XLALregSTRINGUserStruct (	ephemSun,	 0,  UVAR_OPTIONAL,	"Sun ephemeris file to use");

  XLALregSTRINGUserStruct(	outab,		'o', UVAR_OPTIONAL,	"output file for antenna pattern functions a(t), b(t) at each timestamp");
  XLALregSTRINGUserStruct(	outABCD,	'O', UVAR_OPTIONAL,	"output file for antenna pattern matrix elements A, B, C, D averaged over timestamps");

  XLALregBOOLUserStruct(	version,        'V', UVAR_SPECIAL,      "Output code version");

  /* developer user variables */
  XLALregSTRINGUserStruct(	timeStampsFile,	  0, UVAR_OPTIONAL,	"Alternative: single time-stamps file (deprecated, use --timeStampsFiles instead");

  return XLAL_SUCCESS;

} /* XLALInitUserVars() */


/**
 * basic initializations: deal with user input and return standardized 'ConfigVariables'
 */
int
XLALInitCode ( ConfigVariables *cfg, const UserVariables_t *uvar, const char *app_name)
{
  XLAL_CHECK ( cfg && uvar && app_name, XLAL_EINVAL, "Illegal NULL pointer input." );

  /* init ephemeris data */
  XLAL_CHECK ( ( cfg->edat = XLALInitBarycenter( uvar->ephemEarth, uvar->ephemSun ) ) != NULL, XLAL_EFUNC, "XLALInitBarycenter failed: could not load Earth ephemeris '%s' and Sun ephemeris '%s.", uvar->ephemEarth, uvar->ephemSun);

  cfg->numDetectors = uvar->IFOs->length;

  cfg->numTimeStamps = 0;
  XLAL_CHECK ( (cfg->numTimeStampsX = XLALCreateUINT4Vector ( cfg->numDetectors )) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(%d) failed.", cfg->numDetectors );

  BOOLEAN haveTimeGPS = XLALUserVarWasSet( &uvar->timeGPS );
  BOOLEAN haveTimeStampsFile = XLALUserVarWasSet( &uvar->timeStampsFile );
  BOOLEAN haveTimeStampsFiles = XLALUserVarWasSet( &uvar->timeStampsFiles );

  XLAL_CHECK ( !(haveTimeStampsFiles && haveTimeStampsFile), XLAL_EINVAL, "Can't handle both timeStampsFiles and (deprecated) haveTimeStampsFiles input options." );
  XLAL_CHECK ( !(haveTimeGPS && haveTimeStampsFile), XLAL_EINVAL, "Can't handle both (deprecated) timeStampsFile and timeGPS input options." );
  XLAL_CHECK ( !(haveTimeGPS && haveTimeStampsFiles), XLAL_EINVAL, "Can't handle both timeStampsFiles and timeGPS input options." );
  XLAL_CHECK ( haveTimeGPS || haveTimeStampsFiles || haveTimeStampsFile, XLAL_EINVAL, "Need either timeStampsFiles or timeGPS input option." );
  if ( haveTimeStampsFiles ) {
    XLAL_CHECK ( (uvar->timeStampsFiles->length == 1 ) || ( uvar->timeStampsFiles->length == cfg->numDetectors ), XLAL_EINVAL, "Length of timeStampsFiles list is neither 1 (one file for all detectors) nor does it match the number of detectors. (%d != %d)", uvar->timeStampsFiles->length, cfg->numDetectors );
    XLAL_CHECK ( (uvar->timeStampsFiles->length == 1 ) || !uvar->outab, XLAL_EINVAL, "At the moment, can't produce a(t), b(t) output (--outab) when given per-IFO --timeStampsFiles.");
  }

  if ( haveTimeStampsFiles && ( uvar->timeStampsFiles->length == cfg->numDetectors ) ) {

    XLAL_CHECK ( ( cfg->multiTimestamps = XLALReadMultiTimestampsFiles ( uvar->timeStampsFiles ) ) != NULL, XLAL_EFUNC );

    XLAL_CHECK ( (cfg->multiTimestamps->length > 0) && (cfg->multiTimestamps->data != NULL), XLAL_EINVAL, "Got empty timestamps-list from '%s'.", uvar->timeStampsFiles );

  }

  else {

    /* prepare multiTimestamps structure */
    UINT4 nTS = 0;
    XLAL_CHECK ( ( cfg->multiTimestamps = XLALCalloc ( 1, sizeof(*cfg->multiTimestamps))) != NULL, XLAL_ENOMEM, "Allocating multiTimestamps failed." );
    XLAL_CHECK ( ( cfg->multiTimestamps->data = XLALCalloc ( cfg->numDetectors, sizeof(cfg->multiTimestamps->data) )) != NULL, XLAL_ENOMEM, "Allocating multiTimestamps->data failed." );
    cfg->multiTimestamps->length = cfg->numDetectors;

    if ( haveTimeGPS ) { /* set up timestamps vector from timeGPS, use same for all IFOs */

      nTS = uvar->timeGPS->length;
      XLAL_CHECK ( (cfg->multiTimestamps->data[0] = XLALCreateTimestampVector ( nTS ) ) != NULL, XLAL_EFUNC, "XLALCreateTimestampVector( %d ) failed.",  nTS );

      /* convert input REAL8 times into LIGOTimeGPS for first detector */
      for (UINT4 t = 0; t < nTS; t++) {
        REAL8 temp_real8_timestamp = 0;
        XLAL_CHECK ( 1 == sscanf ( uvar->timeGPS->data[t], "%" LAL_REAL8_FORMAT, &temp_real8_timestamp ), XLAL_EINVAL, "Illegal REAL8 commandline argument to --timeGPS[%d]: '%s'", t, uvar->timeGPS->data[t] );
        XLAL_CHECK ( XLALGPSSetREAL8( &cfg->multiTimestamps->data[0]->data[t], temp_real8_timestamp ) != NULL, XLAL_EFUNC, "Failed to convert input GPS %g into LIGOTimeGPS", temp_real8_timestamp );
       } // for (UINT4 t = 0; t < nTS; t++)

    } // if ( haveTimeGPS )

    else { // haveTimeStampsFiles || haveTimeStampsFile

     CHAR *singleTimeStampsFile = NULL;
     if ( haveTimeStampsFiles ) {
      singleTimeStampsFile = uvar->timeStampsFiles->data[0];
     }
     else if ( haveTimeStampsFile ) {
      singleTimeStampsFile = uvar->timeStampsFile;
     }

     XLAL_CHECK ( ( cfg->multiTimestamps->data[0] = XLALReadTimestampsFile ( singleTimeStampsFile ) ) != NULL, XLAL_EFUNC );
     nTS = cfg->multiTimestamps->data[0]->length;

    } // else: haveTimeStampsFiles || haveTimeStampsFile

    /* copy timestamps from first detector to all others */
    if ( cfg->numDetectors > 1 ) {
      for ( UINT4 X=1; X < cfg->numDetectors; X++ ) {
        XLAL_CHECK ( (cfg->multiTimestamps->data[X] = XLALCreateTimestampVector ( nTS ) ) != NULL, XLAL_EFUNC, "XLALCreateTimestampVector( %d ) failed.", nTS );
        for (UINT4 t = 0; t < nTS; t++) {
          cfg->multiTimestamps->data[X]->data[t].gpsSeconds = cfg->multiTimestamps->data[0]->data[t].gpsSeconds;
          cfg->multiTimestamps->data[X]->data[t].gpsNanoSeconds = cfg->multiTimestamps->data[0]->data[t].gpsNanoSeconds;
        } // for (UINT4 t = 0; t < nTS; t++)
      } // for ( UINT4 X=1; X < cfg->numDetectors X++ )
    } // if ( cfg->numDetectors > 1 )

  } // if !( haveTimeStampsFiles && ( uvar->timeStampsFiles->length == cfg->numDetectors ) )

  for ( UINT4 X=0; X < cfg->numDetectors; X++ ) {
    cfg->numTimeStampsX->data[X] = cfg->multiTimestamps->data[X]->length;
    cfg->numTimeStamps += cfg->numTimeStampsX->data[X];
  }

  /* convert detector names into site-info */
  MultiLALDetector multiDet;
  XLAL_CHECK ( XLALParseMultiLALDetector ( &multiDet, uvar->IFOs ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* get detector states */
  XLAL_CHECK ( (cfg->multiDetStates = XLALGetMultiDetectorStates ( cfg->multiTimestamps, &multiDet, cfg->edat, 0.5 * uvar->Tsft )) != NULL, XLAL_EFUNC, "XLALGetDetectorStates() failed." );

  BOOLEAN haveAlphaDelta = ( XLALUserVarWasSet(&uvar->Alpha) && XLALUserVarWasSet(&uvar->Delta) );
  BOOLEAN haveSkyGrid = XLALUserVarWasSet( &uvar->skyGridFile );

  XLAL_CHECK ( !(haveAlphaDelta && haveSkyGrid), XLAL_EINVAL, "Can't handle both Alpha/Delta and skyGridFile input options." );
  XLAL_CHECK ( haveAlphaDelta || haveSkyGrid, XLAL_EINVAL, "Need either Alpha/Delta or skyGridFile input option." );

  if (haveAlphaDelta) { /* parse this into one-element Alpha, Delta vectors */
    XLAL_CHECK ( (cfg->Alpha = XLALCreateREAL8Vector ( 1 )) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(1) failed." );
    cfg->Alpha->data[0] = uvar->Alpha;
    XLAL_CHECK ( (cfg->Delta = XLALCreateREAL8Vector ( 1 )) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(1) failed." );
    cfg->Delta->data[0] = uvar->Delta;
    cfg->numSkyPoints = 1;
  } // if (haveAlphaDelta)

  else if ( haveSkyGrid ) {
    LALParsedDataFile *data = NULL;
    XLAL_CHECK ( XLALParseDataFile (&data, uvar->skyGridFile) == XLAL_SUCCESS, XLAL_EFUNC, "Failed to parse data file '%s'.", uvar->skyGridFile );
    cfg->numSkyPoints = data->lines->nTokens;
    XLAL_CHECK ( (cfg->Alpha = XLALCreateREAL8Vector ( cfg->numSkyPoints )) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector( %d ) failed.", cfg->numSkyPoints  );
    XLAL_CHECK ( (cfg->Delta = XLALCreateREAL8Vector ( cfg->numSkyPoints )) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector( %d ) failed.", cfg->numSkyPoints  );
    for (UINT4 n=0; n < cfg->numSkyPoints; n++) {
      XLAL_CHECK ( 2 == sscanf( data->lines->tokens[n], "%" LAL_REAL8_FORMAT "%" LAL_REAL8_FORMAT, &cfg->Alpha->data[n], &cfg->Delta->data[n] ), XLAL_EDATA, "Could not parse 2 numbers from line %d in candidate-file '%s':\n'%s'", n, uvar->skyGridFile, data->lines->tokens[n] );
    } // for (UINT4 n=0; n < cfg->numSkyPoints; n++)
    XLALDestroyParsedDataFile ( data );
  } // else if ( haveSkyGrid )

  if ( uvar->noiseSqrtShX ) { /* translate user-input PSD sqrt(SX) to noise-weights (this actually does not care whether they were normalized or not) */

    if (  uvar->noiseSqrtShX->length != cfg->numDetectors ) {
      fprintf(stderr, "Length of noiseSqrtShX vector does not match number of detectors! (%d != %d)\n", uvar->noiseSqrtShX->length, cfg->numDetectors);
      XLAL_ERROR ( XLAL_EINVAL );
    }
    REAL8Vector *noiseSqrtShX = NULL;
    if ( (noiseSqrtShX = XLALCreateREAL8Vector ( cfg->numDetectors )) == NULL ) {
      fprintf(stderr, "Failed call to XLALCreateREAL8Vector( %d )\n", cfg->numDetectors );
      XLAL_ERROR ( XLAL_EFUNC );
    }

    REAL8 psd_normalization = 0;

    for (UINT4 X = 0; X < cfg->numDetectors; X++) {

      if ( 1 != sscanf ( uvar->noiseSqrtShX->data[X], "%" LAL_REAL8_FORMAT, &noiseSqrtShX->data[X] ) ) {
        fprintf(stderr, "Illegal REAL8 commandline argument to --noiseSqrtShX[%d]: '%s'\n", X, uvar->noiseSqrtShX->data[X]);
        XLAL_ERROR ( XLAL_EINVAL );
      }

      if ( noiseSqrtShX->data[X] <= 0.0 ) {
        fprintf(stderr, "Non-positive input PSD ratio for detector X=%d: noiseSqrtShX[X]=%f\n", X, noiseSqrtShX->data[X] );
        XLAL_ERROR ( XLAL_EINVAL );
      }

      psd_normalization += 1.0/SQ(noiseSqrtShX->data[X]);

    } /* for X < cfg->numDetectors */

    psd_normalization = (REAL8)cfg->numDetectors/psd_normalization; /* S = NSFT / sum S_Xalpha^-1, no per-SFT variation here -> S = Ndet / sum S_X^-1 */

    /* create multi noise weights */
    if ( (cfg->multiNoiseWeights = XLALCalloc(1, sizeof(*cfg->multiNoiseWeights))) == NULL ) {
     XLALPrintError ("%s: failed to XLALCalloc ( 1, %d )\n", __func__, sizeof(*cfg->multiNoiseWeights) );
     XLAL_ERROR ( XLAL_ENOMEM );
    }
    if ( (cfg->multiNoiseWeights->data = XLALCalloc(cfg->numDetectors, sizeof(*cfg->multiNoiseWeights->data))) == NULL ) {
     XLALPrintError ("%s: failed to XLALCalloc ( %d, %d )\n", __func__, cfg->numDetectors, sizeof(*cfg->multiNoiseWeights->data) );
     XLAL_ERROR ( XLAL_ENOMEM );
    }
    cfg->multiNoiseWeights->length = cfg->numDetectors;

    for (UINT4 X = 0; X < cfg->numDetectors; X++) {

      REAL8 noise_weight_X = psd_normalization/SQ(noiseSqrtShX->data[X]); /* w_Xalpha = S_Xalpha^-1/S^-1 = S / S_Xalpha */

      /* create k^th weights vector */
      if( ( cfg->multiNoiseWeights->data[X] = XLALCreateREAL8Vector ( cfg->numTimeStampsX->data[X] ) ) == NULL )
        {
          /* free weights vectors created previously in loop */
          XLALDestroyMultiNoiseWeights ( cfg->multiNoiseWeights );
          XLAL_ERROR ( XLAL_EFUNC, "Failed to allocate noiseweights for IFO X = %d\n", X );
        } /* if XLALCreateREAL8Vector() failed */

      /* loop over rngmeds and calculate weights -- one for each sft */
      for ( UINT4 alpha = 0; alpha < cfg->numTimeStampsX->data[X]; alpha++) {
        cfg->multiNoiseWeights->data[X]->data[alpha] = noise_weight_X;
      }

    } /* for X < cfg->numDetectors */

    XLALDestroyREAL8Vector ( noiseSqrtShX );

  } /* if ( uvar->noiseSqrtShX ) */

  else {
    cfg->multiNoiseWeights =  NULL;
  }

  return XLAL_SUCCESS;

} /* XLALInitCode() */


/** Destructor for internal configuration struct */
int
XLALDestroyConfig ( ConfigVariables *cfg )
{
  XLAL_CHECK ( cfg != NULL, XLAL_EINVAL );

  XLALDestroyUserVars ();

  XLALDestroyREAL8Vector ( cfg->Alpha );
  XLALDestroyREAL8Vector ( cfg->Delta );

  XLALDestroyMultiTimestamps ( cfg->multiTimestamps );
  XLALDestroyUINT4Vector ( cfg->numTimeStampsX );

  XLALDestroyEphemerisData ( cfg->edat );

  XLALDestroyMultiDetectorStateSeries ( cfg->multiDetStates );
  XLALDestroyMultiNoiseWeights ( cfg->multiNoiseWeights );

  return XLAL_SUCCESS;

} /* XLALDestroyConfig() */
