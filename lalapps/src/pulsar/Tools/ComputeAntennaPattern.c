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
  UINT4 numTimeStamps;			/**< number of timestamps (SFTs) PER DETECTOR */
  UINT4 numSFTstotal;                   /**< number of timestamps (SFTs) for ALL detectors together */
  BOOLEAN averageABCD;			/**< output antenna pattern matrix elements averaged over timestamps, suppress a(t), b(t) */
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
  CHAR  *timeStampsFile;	/**< alternative: read in timestamps from a file (expect same format) */
  BOOLEAN averageABCD; 		/**< output antenna pattern matrix elements averaged over timestamps, suppress a(t), b(t) */
  INT4 Tsft;			/**< assumed length of SFTs, needed for offset to timestamps when comparing to CFS_v2, PFS etc */

  LALStringVector* noiseSqrtShX; /**< per-detector noise PSD sqrt(SX) */

  CHAR *outputFile;	/**< output file to write antenna pattern functions into */

  BOOLEAN version;	/**< output code versions */

} UserVariables_t;


/*---------- empty structs for initializations ----------*/
static ConfigVariables empty_ConfigVariables;
static UserVariables_t empty_UserVariables;

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

  ConfigVariables config = empty_ConfigVariables;
  UserVariables_t uvar = empty_UserVariables;


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

  /* prepare output file */
  FILE *fpOut = NULL;
  if ( uvar.outputFile ) {

      XLAL_CHECK ( (fpOut = fopen (uvar.outputFile, "wb")) != NULL, XLAL_EIO, "Error opening file '%s' for writing...", uvar.outputFile );

      /* write header info in comments */
      XLAL_CHECK ( XLAL_SUCCESS == XLALOutputVersionString ( fpOut, 0 ), XLAL_EFUNC );

      /* write the command-line */
      for (int a = 0; a < argc; a++)
        {
          fprintf(fpOut,"%%%% argv[%d]: '%s'\n", a, argv[a]);
        }

      /* write column headings */
      fprintf(fpOut, "%%%% columns:\n%%%% Alpha  Delta");
      if ( !config.averageABCD ) {
        fprintf(fpOut, "     tGPS       a(t)         b(t)");
      }
      fprintf(fpOut, "         A            B            C            D");
      if ( config.numDetectors > 1 ) {
        fprintf(fpOut, "   ");
        for ( UINT4 X=0; X < config.numDetectors; X++ ) {
          fprintf(fpOut, "         A[%d]         B[%d]         C[%d]         D[%d]", X, X, X, X);
        }
      }
      fprintf(fpOut, "\n");

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

    if ( uvar.outputFile ) {

      /* write out the data for this sky point*/
     if ( config.averageABCD ) { // output only ABCD averaged over all timestamps, suppress a(t), b(t) (no meaningful use known for their averages)
       // FIXME: stop doing average manually when AMCoeffs is changed to contain averaged values
       REAL8 A = multiAM->Mmunu.Ad/config.numSFTstotal;
       REAL8 B = multiAM->Mmunu.Bd/config.numSFTstotal;
       REAL8 C = multiAM->Mmunu.Cd/config.numSFTstotal;
       REAL8 D = A*B-SQ(C);
       fprintf (fpOut, "%.7f %.7f %12.8f %12.8f %12.8f %12.8f", config.Alpha->data[n], config.Delta->data[n], A, B, C, D );
       if ( config.numDetectors > 1 ) {
         fprintf(fpOut, "   ");
         for ( UINT4 X=0; X < config.numDetectors; X++ ) {
           REAL4 AX = multiAM->data[X]->A/config.numTimeStamps;
           REAL4 BX = multiAM->data[X]->B/config.numTimeStamps;
           REAL4 CX = multiAM->data[X]->C/config.numTimeStamps;
           REAL4 DX = AX*BX-SQ(CX);
           fprintf(fpOut, " %12.8f %12.8f %12.8f %12.8f", AX, BX, CX, DX);
         }
       }
       fprintf(fpOut, "\n");
     } // if ( config.averageABCD )

     else { // output all values at each timestamp
       for (UINT4 t = 0; t < config.numTimeStamps; t++) {
         REAL4 A = SQ(multiAM->data[0]->a->data[t]);
         REAL4 B = SQ(multiAM->data[0]->b->data[t]);
         REAL4 C = multiAM->data[0]->a->data[t]*multiAM->data[0]->b->data[t];
         REAL4 D = A*B-SQ(C);
         fprintf (fpOut, "%.7f %.7f %d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f", config.Alpha->data[n], config.Delta->data[n], config.multiTimestamps->data[0]->data[t].gpsSeconds, multiAM->data[0]->a->data[t], multiAM->data[0]->b->data[t], A, B, C, D );
         if ( config.numDetectors > 1 ) {
           fprintf(fpOut, "   ");
           for ( UINT4 X=0; X < config.numDetectors; X++ ) {
             REAL4 AX = multiAM->data[X]->A/config.numTimeStamps;
             REAL4 BX = multiAM->data[X]->B/config.numTimeStamps;
             REAL4 CX = multiAM->data[X]->C/config.numTimeStamps;
             REAL4 DX = AX*BX-SQ(CX);
             fprintf(fpOut, " %12.8f %12.8f %12.8f %12.8f", AX, BX, CX, DX);
           }
         }
         fprintf(fpOut, "\n");
       }
     } // if ( !config.averageABCD )

    } /* if outputFile */

  XLALDestroyMultiAMCoeffs ( multiAM );

  } // for (UINT4 n = 0; n < config.numSkyPoints; n++)

  /* ----- close output file ----- */
  fprintf (fpOut, "\n");
  if ( fpOut ) fclose ( fpOut );

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
  uvar->averageABCD = 0;
  uvar->Tsft = 1800;

  uvar->noiseSqrtShX = NULL;

  uvar->outputFile = NULL;

  /* register all user-variables */
  XLALregBOOLUserStruct(	help,		'h', UVAR_HELP,		"Print this help/usage message");
  XLALregLISTUserStruct( IFOs,                  'I', UVAR_OPTIONAL, "Comma-separated list of detectors, eg. \"H1,H2,L1,G1, ...\" [only 1 detector supported at the moment] ");

  XLALregREALUserStruct(	Alpha,		'a', UVAR_OPTIONAL,	"single skyposition Alpha in radians, equatorial coords.");
  XLALregREALUserStruct(	Delta, 		'd', UVAR_OPTIONAL,	"single skyposition Delta in radians, equatorial coords.");

  XLALregSTRINGUserStruct( skyGridFile,		's', UVAR_OPTIONAL,	"Alternatively: sky-grid file");

  XLALregLISTUserStruct( 	timeGPS,        't', UVAR_OPTIONAL, 	"GPS time at which to compute detector states (separate multiple timestamps by commata)");
  XLALregSTRINGUserStruct(	timeStampsFile, 'T', UVAR_OPTIONAL,	"Alternative: time-stamps file");
  XLALregBOOLUserStruct(	averageABCD,	0, UVAR_OPTIONAL,	"output only time-averaged antenna pattern matrix elements");
  XLALregINTUserStruct(		Tsft,		0, UVAR_OPTIONAL,	"Assumed length of one SFT in seconds; needed for timestamps offset consistency with F-stat based codes");

  XLALregLISTUserStruct ( noiseSqrtShX,		 0, UVAR_OPTIONAL, "Per-detector noise PSD sqrt(SX). Only ratios relevant to compute noise weights. Defaults to 1,1,...");

  XLALregSTRINGUserStruct (	ephemEarth,	 0,  UVAR_OPTIONAL,	"Earth ephemeris file to use");
  XLALregSTRINGUserStruct (	ephemSun,	 0,  UVAR_OPTIONAL,	"Sun ephemeris file to use");

  XLALregSTRINGUserStruct(	outputFile, 	'o', UVAR_OPTIONAL,     "Output file for antenna pattern functions");

  XLALregBOOLUserStruct(	version,        'V', UVAR_SPECIAL,      "Output code version");

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

  BOOLEAN haveTimeGPS = XLALUserVarWasSet( &uvar->timeGPS );
  BOOLEAN haveTimeStampsFile = XLALUserVarWasSet( &uvar->timeStampsFile );

  XLAL_CHECK ( !(haveTimeGPS && haveTimeStampsFile), XLAL_EINVAL, "Can't handle both timeStampsFile and timeGPS input options." );
  XLAL_CHECK ( haveTimeGPS || haveTimeStampsFile, XLAL_EINVAL, "Need either timeStampsFile or timeGPS input option." );

  /* FIXME: only using identical timestamps for all detectors */
  XLAL_CHECK ( ( cfg->multiTimestamps = XLALCalloc ( 1, sizeof(*cfg->multiTimestamps))) != NULL, XLAL_ENOMEM, "Allocating multiTimestamps failed." );
  XLAL_CHECK ( ( cfg->multiTimestamps->data = XLALCalloc ( cfg->numDetectors, sizeof(cfg->multiTimestamps->data) )) != NULL, XLAL_ENOMEM, "Allocating multiTimestamps->data failed." );
  cfg->multiTimestamps->length = cfg->numDetectors;

  if ( haveTimeStampsFile ) { // read in timestamps vector from file
    XLAL_CHECK ( (cfg->multiTimestamps->data[0] = XLALReadTimestampsFile ( uvar->timeStampsFile )) != NULL, XLAL_EFUNC, "illegal NULL pointer returned when reading timestampsfile '%s'.", uvar->timeStampsFile );
    cfg->numTimeStamps = cfg->multiTimestamps->data[0]->length;
  } // timestamps from timeStampsFile

  else if ( haveTimeGPS ) { // set up timestamps vector from timeGPS

    cfg->numTimeStamps = uvar->timeGPS->length;

    XLAL_CHECK ( (cfg->multiTimestamps->data[0] = XLALCreateTimestampVector ( cfg->numTimeStamps ) ) != NULL, XLAL_EFUNC, "XLALCreateTimestampVector( %d ) failed.",  cfg->numTimeStamps );

    /* convert input REAL8 time into LIGOTimeGPS */
    for (UINT4 t = 0; t < cfg->numTimeStamps; t++) {
      REAL8 temp_real8_timestamp = 0;
      XLAL_CHECK ( 1 == sscanf ( uvar->timeGPS->data[t], "%" LAL_REAL8_FORMAT, &temp_real8_timestamp ), XLAL_EINVAL, "Illegal REAL8 commandline argument to --timeGPS[%d]: '%s'", t, uvar->timeGPS->data[t] );
      XLAL_CHECK ( XLALGPSSetREAL8( &cfg->multiTimestamps->data[0]->data[t], temp_real8_timestamp ) != NULL, XLAL_EFUNC, "Failed to convert input GPS %g into LIGOTimeGPS", temp_real8_timestamp );
    } // for (UINT4 t = 0; t < cfg->numTimeStamps; t++)

  } // timestamps from timeGPS

  /* in either case, copy timestamps from first detector to all others (for now) */
  if ( cfg->numDetectors > 1 ) {
    for ( UINT4 X=1; X < cfg->numDetectors; X++ ) {
      XLAL_CHECK ( (cfg->multiTimestamps->data[X] = XLALCreateTimestampVector ( cfg->numTimeStamps ) ) != NULL, XLAL_EFUNC, "XLALCreateTimestampVector( %d ) failed.", cfg->numTimeStamps );
      for (UINT4 t = 0; t < cfg->numTimeStamps; t++) {
       cfg->multiTimestamps->data[X]->data[t].gpsSeconds = cfg->multiTimestamps->data[0]->data[t].gpsSeconds;
       cfg->multiTimestamps->data[X]->data[t].gpsNanoSeconds = cfg->multiTimestamps->data[0]->data[t].gpsNanoSeconds;
      } // for (UINT4 t = 0; t < cfg->numTimeStamps; t++)
    } // for ( UINT4 X=1; X < cfg->numDetectors; X++ )
  } // if ( cfg->numDetectors > 1 )

  cfg->numTimeStamps = cfg->multiTimestamps->data[0]->length;

  cfg->numSFTstotal = cfg->numDetectors*cfg->numTimeStamps;

  cfg->averageABCD = uvar->averageABCD;

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
      if( ( cfg->multiNoiseWeights->data[X] = XLALCreateREAL8Vector ( cfg->numTimeStamps ) ) == NULL )
        {
          /* free weights vectors created previously in loop */
          XLALDestroyMultiNoiseWeights ( cfg->multiNoiseWeights );
          XLAL_ERROR ( XLAL_EFUNC, "Failed to allocate noiseweights for IFO X = %d\n", X );
        } /* if XLALCreateREAL8Vector() failed */

      /* loop over rngmeds and calculate weights -- one for each sft */
      for ( UINT4 alpha = 0; alpha < cfg->numTimeStamps; alpha++) {
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

  XLALDestroyEphemerisData ( cfg->edat );

  XLALDestroyMultiDetectorStateSeries ( cfg->multiDetStates );
  XLALDestroyMultiNoiseWeights ( cfg->multiNoiseWeights );

  return XLAL_SUCCESS;

} /* XLALDestroyConfig() */
