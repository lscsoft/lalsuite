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
  MultiDetectorStateSeries *multiDetStates;	/**< detector state time series */
  MultiNoiseWeights *multiWeights;		/**< per-detector noise weights */
  MultiLIGOTimeGPSVector *multiTimestamps;	/**< timestamps vector (LIGOtimeGPS format) */
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
      fprintf(fpOut, "         A            B            C            D\n");

    }

  UINT4 numTimeStamps = config.multiTimestamps->data[0]->length;

  /* loop over sky positions (outer loop, could allow for buffering if necessary) */
  for (UINT4 n = 0; n < config.numSkyPoints; n++) {
    SkyPosition skypos;
    skypos.system = COORDINATESYSTEM_EQUATORIAL;
    skypos.longitude = config.Alpha->data[n];
    skypos.latitude  = config.Delta->data[n];

    /* do the actual computation of the antenna pattern functions */
    MultiAMCoeffs *multiAM;
    XLAL_CHECK ( ( multiAM = XLALComputeMultiAMCoeffs ( config.multiDetStates, config.multiWeights, skypos ) ) != NULL, XLAL_EFUNC, "XLALComputeAMCoeffs() failed." );

    if ( uvar.outputFile ) {

      /* write out the data for this sky point*/
     if ( config.averageABCD ) { // output only ABCD averaged over all timestamps, suppress a(t), b(t) (no meaningful use known for their averages)
       // FIXME: stop doing average manually when AMCoeffs is changed to contain averaged values
       REAL4 A = multiAM->data[0]->A/numTimeStamps;
       REAL4 B = multiAM->data[0]->B/numTimeStamps;
       REAL4 C = multiAM->data[0]->C/numTimeStamps;
       REAL4 D = A*B-SQ(C);
       fprintf (fpOut, "%.7f %.7f %12.8f %12.8f %12.8f %12.8f\n", config.Alpha->data[n], config.Delta->data[n], A, B, C, D );
     } // if ( config.averageABCD )

     else { // output all values at each timestamp
       for (UINT4 t = 0; t < numTimeStamps; t++) {
         REAL4 A = SQ(multiAM->data[0]->a->data[t]);
         REAL4 B = SQ(multiAM->data[0]->b->data[t]);
         REAL4 C = multiAM->data[0]->a->data[t]*multiAM->data[0]->b->data[t];
         REAL4 D = A*B-SQ(C);
         fprintf (fpOut, "%.7f %.7f %d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", config.Alpha->data[n], config.Delta->data[n], config.multiTimestamps->data[0]->data[t].gpsSeconds, multiAM->data[0]->a->data[t], multiAM->data[0]->b->data[t], A, B, C, D );
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

  UINT4 numDetectors = uvar->IFOs->length;
  XLAL_CHECK ( numDetectors == 1, XLAL_EINVAL, "Can't handle more than one IFO at the moment." ); // FIXME

  BOOLEAN haveTimeGPS = XLALUserVarWasSet( &uvar->timeGPS );
  BOOLEAN haveTimeStampsFile = XLALUserVarWasSet( &uvar->timeStampsFile );

  XLAL_CHECK ( !(haveTimeGPS && haveTimeStampsFile), XLAL_EINVAL, "Can't handle both timeStampsFile and timeGPS input options." );
  XLAL_CHECK ( haveTimeGPS || haveTimeStampsFile, XLAL_EINVAL, "Need either timeStampsFile or timeGPS input option." );

  /* FIXME: only using identical timestamps for all detectors */
  XLAL_CHECK ( ( cfg->multiTimestamps = XLALCalloc ( 1, sizeof(*cfg->multiTimestamps))) != NULL, XLAL_ENOMEM, "Allocating multiTimestamps failed." );
  XLAL_CHECK ( ( cfg->multiTimestamps->data = XLALCalloc ( numDetectors, sizeof(cfg->multiTimestamps->data[0]) )) != NULL, XLAL_ENOMEM, "Allocating multiTimestamps->data failed." );
  cfg->multiTimestamps->length = numDetectors;

  if ( haveTimeStampsFile ) { // read in timestamps vector from file
    XLAL_CHECK ( (cfg->multiTimestamps->data[0] = XLALReadTimestampsFile ( uvar->timeStampsFile )) != NULL, XLAL_EFUNC, "illegal NULL pointer returned when reading timestampsfile '%s'.", uvar->timeStampsFile );
  } // timestamps from timeStampsFile

  else if ( haveTimeGPS ) { // set up timestamps vector from timeGPS

    XLAL_CHECK ( (cfg->multiTimestamps->data[0] = XLALCreateTimestampVector ( uvar->timeGPS->length ) ) != NULL, XLAL_EFUNC, "XLALCreateTimestampVector( %d ) failed.", uvar->timeGPS->length );

    /* convert input REAL8 time into LIGOTimeGPS */
    for (UINT4 t = 0; t < uvar->timeGPS->length; t++) {
      REAL8 temp_real8_timestamp = 0;
      XLAL_CHECK ( 1 == sscanf ( uvar->timeGPS->data[t], "%" LAL_REAL8_FORMAT, &temp_real8_timestamp ), XLAL_EINVAL, "Illegal REAL8 commandline argument to --timeGPS[%d]: '%s'", t, uvar->timeGPS->data[t] );
      XLAL_CHECK ( XLALGPSSetREAL8( &cfg->multiTimestamps->data[0]->data[t], temp_real8_timestamp ) != NULL, XLAL_EFUNC, "Failed to convert input GPS %g into LIGOTimeGPS", temp_real8_timestamp );
    } // for (UINT4 t = 0; t < uvar->timeGPS->length; t++)

  } // timestamps from timeGPS

  /* in either case, copy timestamps from first detector to all others (for now) */
  if ( numDetectors > 1 ) {
    for ( UINT4 X=1; X < numDetectors; X++ ) {
      XLAL_CHECK ( (cfg->multiTimestamps->data[X] = XLALCreateTimestampVector ( uvar->timeGPS->length ) ) != NULL, XLAL_EFUNC, "XLALCreateTimestampVector( %d ) failed.", uvar->timeGPS->length );
      for (UINT4 t = 0; t < uvar->timeGPS->length; t++) {
       cfg->multiTimestamps->data[X]->data[t].gpsSeconds = cfg->multiTimestamps->data[0]->data[t].gpsNanoSeconds;
       cfg->multiTimestamps->data[X]->data[t].gpsSeconds = cfg->multiTimestamps->data[0]->data[t].gpsNanoSeconds;
      } // for (UINT4 t = 0; t < uvar->timeGPS->length; t++)
    } // for ( UINT4 X=1; X < numDetectors; X++ )
  } // if ( numDetectors > 1 )

  cfg->averageABCD = uvar->averageABCD;

  /* convert detector names into site-info */
  MultiLALDetector *multiDet;
  XLAL_CHECK ( (multiDet = XLALCreateMultiLALDetector ( numDetectors )) != NULL, XLAL_EFUNC, "XLALCreateMultiLALDetector(1) failed." );
  LALDetector *site = NULL;
  for ( UINT4 X=0; X < numDetectors; X++ ) {
    XLAL_CHECK ( (site = XLALGetSiteInfo ( uvar->IFOs->data[X] )) != NULL, XLAL_EFUNC, "Failed to get site-info for detector '%s'.", uvar->IFOs->data[X] );
    multiDet->data[X] = (*site); 	/* copy! */
    XLALFree ( site );
  }

  /* get detector states */
  XLAL_CHECK ( (cfg->multiDetStates = XLALGetMultiDetectorStates ( cfg->multiTimestamps, multiDet, cfg->edat, 0.5 * uvar->Tsft )) != NULL, XLAL_EFUNC, "XLALGetDetectorStates() failed." );

  XLALDestroyMultiLALDetector ( multiDet );

  cfg->multiWeights =  NULL; // FIXME: noise weights not supported (yet)

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

  return XLAL_SUCCESS;

} /* XLALDestroyConfig() */
