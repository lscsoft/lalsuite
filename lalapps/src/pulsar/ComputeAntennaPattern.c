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
 * for given detector and sky-location
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

#include <lalapps.h>

/* ----- compile switches ----- */

/*---------- local defines ---------- */

/* ----- Macros ----- */
#define SQ(x) ( (x) * (x) )

/* ---------- local types ---------- */

/** types of mathematical operations */
enum tagMATH_OP_TYPE {
  MATH_OP_ARITHMETIC_SINGLE = 0, /**< individual x values */
  MATH_OP_ARITHMETIC_SUM,        /**< sum(x)     */
  MATH_OP_ARITHMETIC_MEAN,       /**< sum(x) / n */
};


typedef struct
{
  EphemerisData *edat;			/**< ephemeris data (from LALInitBarycenter()) */
  LALDetector *det;			/**< LIGODetector struct holding detector info */
  LIGOTimeGPSVector *timestamps;	/**< timestamps vector (LIGOtimeGPS format) */
  UINT4 mthopTimeStamps;		/**< math operation over timestamps */
  REAL8 sinzeta;			/**< detector-arm angle correction (needed for GEO) */
  REAL8Vector *Alpha;			/**< skyposition Alpha: radians, equatorial coords */
  REAL8Vector *Delta;			/**< skyposition Delta: radians, equatorial coords */
  UINT4 numTimeStamps;			/**< length of timestamps vector */
  UINT4 numSkyPoints;			/**< common length of Alpha and Delta vectors */
  UINT4 numData;			/**< total number of data points, numTimeStamps*numSkyPoints */
} ConfigVariables;


typedef struct
{
  BOOLEAN help;

  CHAR *detector;	/**< detector name */

  REAL8 Alpha;		/**< a single skyposition Alpha: radians, equatorial coords. */
  REAL8 Delta;		/**< a single skyposition Delta: radians, equatorial coords. */
  CHAR *skyGridFile;	/**< alternative: matrix of (Alpha,Delta) pairs from a file */

  CHAR *ephemDir;	/**< directory to look for ephemeris files */
  CHAR *ephemYear;	/**< date-range string on ephemeris-files to use */

  LALStringVector* timeGPS;	/**< GPS timestamps to compute detector state for (REAL8 format) */
  CHAR  *timeStampsFile;	/**< alternative: read in timestamps from a file (expect same format) */
  INT4 mthopTimeStamps; 	/**< math operation over timestamps */
  INT4 timeStampOffset;		/**< offset to timestamps, needed when comparing to CFS_v2, PFS etc */

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

EphemerisData *InitEphemeris (const CHAR *ephemDir, const CHAR *ephemYear );

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

  REAL8Vector *a0, *b0;
  if ( (a0 = XLALCreateREAL8Vector ( config.numData )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateREAL8Vector( %d )\n", __func__, config.numData );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  if ( (b0 = XLALCreateREAL8Vector ( config.numData )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateREAL8Vector( %d )\n", __func__, config.numData );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* loop over sky positions (outer loop, could allow for buffering if necessary) */
  for (UINT4 n = 0; n < config.numSkyPoints; n++) {
    SkyPosition skypos;
    skypos.system = COORDINATESYSTEM_EQUATORIAL;
    skypos.longitude = config.Alpha->data[n];
    skypos.latitude  = config.Delta->data[n];

    /* do the actual computation of the antenna pattern functions */
    for (UINT4 t = 0; t < config.numTimeStamps; t++) {
      if ( XLALComputeAntennaPatternCoeffs ( &a0->data[n*config.numTimeStamps+t], &b0->data[n*config.numTimeStamps+t], &skypos, &config.timestamps->data[t], config.det, config.edat ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: XLALComputeAntennaPatternCoeffs() failed.\n", __func__ );
        XLAL_ERROR ( XLAL_EFUNC );
      }
      a0->data[n*config.numTimeStamps+t] /= config.sinzeta;
      b0->data[n*config.numTimeStamps+t] /= config.sinzeta;
    } // for (UINT4 t = 0; t < config.numTimeStamps; t++)

  } // for (UINT4 n = 0; n < config.numSkyPoints; n++)

  /* output the results */
  FILE *fpOut = NULL;
  if ( uvar.outputFile )
    {

      if ( (fpOut = fopen (uvar.outputFile, "wb")) == NULL)
        {
          XLALPrintError ("%s: Error opening file '%s' for writing...\n", __func__, uvar.outputFile );
          XLAL_ERROR ( XLAL_EIO );
        }

      /* write header info in comments */
      if ( XLAL_SUCCESS != XLALOutputVersionString ( fpOut, 0 ) )
        {
          XLAL_ERROR ( XLAL_EFUNC );
        }

      /* write the command-line */
      for (int a = 0; a < argc; a++)
        {
          fprintf(fpOut,"%%%% argv[%d]: '%s'\n", a, argv[a]);
        }

      /* write column headings */
      fprintf(fpOut, "%%%% columns:\n%%%% Alpha  Delta");
      if ( config.mthopTimeStamps == MATH_OP_ARITHMETIC_SINGLE ) {
        fprintf(fpOut, "     tGPS");
      }
      fprintf(fpOut, "       a(t)         b(t)         A            B            C            D\n");

      /* write the actual data */
      for (UINT4 n = 0; n < config.numSkyPoints; n++) {

       if ( config.mthopTimeStamps == MATH_OP_ARITHMETIC_SINGLE ) { // in this case, just output all antenna pattern function values at each timestamp
         for (UINT4 t = 0; t < config.numTimeStamps; t++) {
           REAL8 A = SQ(a0->data[n*config.numTimeStamps+t]);
           REAL8 B = SQ(b0->data[n*config.numTimeStamps+t]);
           REAL8 C = a0->data[n*config.numTimeStamps+t]*b0->data[n*config.numTimeStamps+t];
           REAL8 D = A*B-SQ(C);
           fprintf (fpOut, "%.7f %.7f %d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", config.Alpha->data[n], config.Delta->data[n], config.timestamps->data[t].gpsSeconds, a0->data[n*config.numTimeStamps+t], b0->data[n*config.numTimeStamps+t], A, B, C, D );
         }
       } // if ( config.mthopTimeStamps == MATH_OP_ARITHMETIC_SINGLE )

       else { // in all other cases, output a single antenna pattern function value per sky point computed over all timestamps
         REAL8 atotn = 0;
         REAL8 btotn = 0;
         REAL8 A = 0;
         REAL8 B = 0;
         REAL8 C = 0;
         if ( ( config.mthopTimeStamps == MATH_OP_ARITHMETIC_SUM ) || ( config.mthopTimeStamps == MATH_OP_ARITHMETIC_MEAN ) ) {
           for (UINT4 t = 0; t < config.numTimeStamps; t++) {
             atotn += a0->data[n*config.numTimeStamps+t];
             btotn += b0->data[n*config.numTimeStamps+t];
             A += SQ(a0->data[n*config.numTimeStamps+t]);
             B += SQ(b0->data[n*config.numTimeStamps+t]);
             C += a0->data[n*config.numTimeStamps+t]*b0->data[n*config.numTimeStamps+t];
           }
         }
         if ( config.mthopTimeStamps == MATH_OP_ARITHMETIC_MEAN ) {
           atotn /= config.numTimeStamps;
           btotn /= config.numTimeStamps;
           A /= config.numTimeStamps;
           B /= config.numTimeStamps;
           C /= config.numTimeStamps;
         }
         REAL8 D = A*B-SQ(C);
         fprintf (fpOut, "%.7f %.7f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", config.Alpha->data[n], config.Delta->data[n], atotn, btotn, A, B, C, D );
       } // !if ( mthopTimeStamps == MATH_OP_ARITHMETIC_SINGLE )

      } // for (UINT4 n = 0; n < config.numSkyPoints; n++)
      fprintf (fpOut, "\n");

    } /* if outputFile */

  /* ----- close output file ----- */
  if ( fpOut ) fclose ( fpOut );

  /* ----- done: free all memory */
  XLALDestroyREAL8Vector ( a0 );
  XLALDestroyREAL8Vector ( b0 );
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

#define EPHEM_YEAR  "00-19-DE405"
  XLAL_CHECK ( (uvar->ephemYear = XLALCalloc (1, strlen(EPHEM_YEAR)+1)) != NULL, XLAL_ENOMEM );
  strcpy (uvar->ephemYear, EPHEM_YEAR);

  uvar->Alpha     = 0.0;
  uvar->Delta     = 0.0;
  uvar->skyGridFile = NULL;

  uvar->timeGPS = NULL;
  uvar->timeStampsFile = NULL;
  uvar->mthopTimeStamps = MATH_OP_ARITHMETIC_SINGLE;
  uvar->timeStampOffset = 0;

  uvar->outputFile = NULL;

  /* register all user-variables */
  XLALregBOOLUserStruct(	help,		'h', UVAR_HELP,		"Print this help/usage message");
  XLALregSTRINGUserStruct( 	detector,	'I', UVAR_REQUIRED, 	"Detector name (eg. H1,H2,L1,G1,etc).");

  XLALregREALUserStruct(	Alpha,		'a', UVAR_OPTIONAL,	"single skyposition Alpha in radians, equatorial coords.");
  XLALregREALUserStruct(	Delta, 		'd', UVAR_OPTIONAL,	"single skyposition Delta in radians, equatorial coords.");

  XLALregSTRINGUserStruct( skyGridFile,		's', UVAR_OPTIONAL,	"Alternatively: sky-grid file");

  XLALregLISTUserStruct( 	timeGPS,        't', UVAR_OPTIONAL, 	"GPS time at which to compute detector state (separate multiple timestamps by commata)");
  XLALregSTRINGUserStruct(	timeStampsFile, 'T', UVAR_OPTIONAL,	"Alternative: time-stamps file");
  XLALregINTUserStruct(		mthopTimeStamps,'m', UVAR_OPTIONAL,	"type of math. operation over timestamps: 0=individual values, 1=arith-sum, 2=arith-mean");
  XLALregINTUserStruct(		timeStampOffset,'O', UVAR_OPTIONAL,	"Compute at offset from given timestamps - usually Tsft/2 when comparing to F-stat based codes");

  XLALregSTRINGUserStruct(	ephemDir, 	'E', UVAR_OPTIONAL,     "Directory where Ephemeris files are located");
  XLALregSTRINGUserStruct(	ephemYear, 	'y', UVAR_OPTIONAL,     "Year (or range of years) of ephemeris files to be used");

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
  if ( !cfg || !uvar || !app_name ) {
    XLALPrintError ("%s: illegal NULL pointer input.\n\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* init ephemeris data */
  if ( (cfg->edat = InitEphemeris ( uvar->ephemDir, uvar->ephemYear)) == NULL ) {
    XLALPrintError ("%s: InitEphemeris() Failed to initialize ephemeris data!\n\n", __func__);
    XLAL_ERROR ( XLAL_EINVAL );
  }

  BOOLEAN haveTimeGPS = XLALUserVarWasSet( &uvar->timeGPS );
  BOOLEAN haveTimeStampsFile = XLALUserVarWasSet( &uvar->timeStampsFile );

  if ( haveTimeGPS && haveTimeStampsFile ) {
    XLALPrintError ("%s: Can't handle both timeStampsFile and timeGPS input options.\n\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  if ( !haveTimeGPS && !haveTimeStampsFile ) {
    XLALPrintError ("%s: Need either timeStampsFile or timeGPS input options.\n\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  if ( haveTimeStampsFile ) { // read in timestamps vector from file
    if ( (cfg->timestamps = XLALReadTimestampsFile ( uvar->timeStampsFile )) == NULL ) {
      XLALPrintError ("%s: illegal NULL pointer returned when reading timestampsfile '%s'.\n\n", __func__, uvar->timeStampsFile );
      XLAL_ERROR ( XLAL_EFUNC );
    }
  } // timestamps from timeStampsFile

  else if ( haveTimeGPS ) { // set up timestamps vector from timeGPS

    if ( (cfg->timestamps = XLALCreateTimestampVector( uvar->timeGPS->length ) ) == NULL ) {
      XLALPrintError ("%s: XLALCreateTimestampVector( %d ) failed.", __func__, uvar->timeGPS->length );
      XLAL_ERROR ( XLAL_EFUNC );
    }

    /* convert input REAL8 time into LIGOTimeGPS */
    for (UINT4 t = 0; t < uvar->timeGPS->length; t++) {
      REAL8 temp_real8_timestamp = 0;
      if ( 1 != sscanf ( uvar->timeGPS->data[t], "%" LAL_REAL8_FORMAT, &temp_real8_timestamp ) ) {
        fprintf(stderr, "Illegal REAL8 commandline argument to --timeGPS[%d]: '%s'\n", t, uvar->timeGPS->data[t]);
        XLAL_ERROR ( XLAL_EINVAL );
      }
      if ( XLALGPSSetREAL8( &cfg->timestamps->data[t], temp_real8_timestamp ) == NULL ) {
        XLALPrintError ("%s: failed to convert input GPS %g into LIGOTimeGPS\n", __func__, temp_real8_timestamp );
        XLAL_ERROR ( XLAL_EFUNC );
      }
    } // for (UINT4 t = 0; t < uvar->timeGPS->length; t++)

  } // timestamps from timeGPS

  cfg->numTimeStamps = cfg->timestamps->length;
  cfg->mthopTimeStamps = uvar->mthopTimeStamps;

  if ( uvar->timeStampOffset > 0 ) { // apply an offset, e.g. Tsft/2, to all timestamps
    for (UINT4 t = 0; t < cfg->numTimeStamps; t++) {
      XLALGPSAdd ( &cfg->timestamps->data[t], uvar->timeStampOffset );
    }
  }

  /* convert detector name into site-info */
  if ( ( cfg->det = XLALGetSiteInfo ( uvar->detector )) == NULL )
    {
      XLALPrintError ("%s: XLALGetSiteInfo('%s') failed.\n", __func__, uvar->detector );
      XLAL_ERROR ( XLAL_EFUNC );
    }

  /* NOTE: contrary to ComputeAM() and LALGetAMCoffs(), the new function LALNewGetAMCoeffs()
   * computes 'a * sinzeta' and 'b * sinzeta': for the comparison we therefore need to correct
   * for GEO's opening-angle of 94.33degrees [JKS98]: */
  if ( ! strcmp ( cfg->det->frDetector.name, "GEO_600" ) )
    cfg->sinzeta = 0.997146;
  else
    cfg->sinzeta = 1;


  BOOLEAN haveAlphaDelta = ( XLALUserVarWasSet(&uvar->Alpha) && XLALUserVarWasSet(&uvar->Delta) );
  BOOLEAN haveSkyGrid = XLALUserVarWasSet( &uvar->skyGridFile );

  if ( haveAlphaDelta && haveSkyGrid ) {
    XLALPrintError ("%s: Can't handle both Alpha/Delta and skyGridFile input options.\n\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  if ( !haveAlphaDelta && !haveSkyGrid ) {
    XLALPrintError ("%s: Need either Alpha/Delta or skyGridFile input options.\n\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  if (haveAlphaDelta) { /* parse this into one-element Alpha, Delta vectors */
    if ( (cfg->Alpha = XLALCreateREAL8Vector ( 1 )) == NULL ) {
      XLALPrintError ("%s: failed to XLALCreateREAL8Vector( %d )\n", __func__, 1 );
      XLAL_ERROR ( XLAL_EFUNC );
    }
    cfg->Alpha->data[0] = uvar->Alpha;
    if ( (cfg->Delta = XLALCreateREAL8Vector ( 1 )) == NULL ) {
      XLALPrintError ("%s: failed to XLALCreateREAL8Vector( %d )\n", __func__, 1 );
      XLAL_ERROR ( XLAL_EFUNC );
    }
    cfg->Delta->data[0] = uvar->Delta;
    cfg->numSkyPoints = 1;
  } // if (haveAlphaDelta)

  else if ( haveSkyGrid ) {
    LALParsedDataFile *data = NULL;
    if ( XLALParseDataFile (&data, uvar->skyGridFile) != XLAL_SUCCESS ) {
      XLALPrintError ("%s: failed to parse data file '%s'\n", __func__, uvar->skyGridFile );
      XLAL_ERROR ( XLAL_EFUNC );
    }
    cfg->numSkyPoints = data->lines->nTokens;
    if ( (cfg->Alpha = XLALCreateREAL8Vector ( cfg->numSkyPoints )) == NULL ) {
      XLALPrintError ("%s: failed to XLALCreateREAL8Vector( %d )\n", __func__, cfg->numSkyPoints );
      XLAL_ERROR ( XLAL_EFUNC );
    }
    if ( (cfg->Delta = XLALCreateREAL8Vector ( cfg->numSkyPoints )) == NULL ) {
      XLALPrintError ("%s: failed to XLALCreateREAL8Vector( %d )\n", __func__, cfg->numSkyPoints );
      XLAL_ERROR ( XLAL_EFUNC );
    }
    for (UINT4 n=0; n < cfg->numSkyPoints; n++) {
      if ( 2 != sscanf( data->lines->tokens[n], "%" LAL_REAL8_FORMAT "%" LAL_REAL8_FORMAT, &cfg->Alpha->data[n], &cfg->Delta->data[n] ))
        {
          XLALPrintError ( "%s: could not parse 2 numbers from line %d in candidate-file '%s':\n", __func__, n, uvar->skyGridFile);
          XLALPrintError ("'%s'\n", data->lines->tokens[n] );
          XLAL_ERROR (   XLAL_EDATA );
        }
    } // for (UINT4 n=0; n < cfg->numSkyPoints; n++)
    XLALDestroyParsedDataFile ( data );
  } // else if ( haveSkyGrid )

  cfg->numData = cfg->timestamps->length*cfg->numSkyPoints;

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

  XLALDestroyTimestampVector ( cfg->timestamps );

  XLALDestroyEphemerisData ( cfg->edat );

  XLALFree ( cfg->det );

  return XLAL_SUCCESS;

} /* XLALDestroyConfig() */



/** Load Ephemeris from ephemeris data-files  */
EphemerisData *
InitEphemeris (const CHAR *ephemDir,	/**< directory containing ephems */
	       const CHAR *ephemYear	/**< which years do we need? */
	       )
{
#define FNAME_LENGTH 1024
  EphemerisData *edat;
  CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */

  XLAL_CHECK_NULL ( ephemYear != NULL, XLAL_EINVAL );

  if ( ephemDir )
    {
      snprintf(EphemEarth, FNAME_LENGTH, "%s/earth%s.dat", ephemDir, ephemYear);
      snprintf(EphemSun, FNAME_LENGTH, "%s/sun%s.dat", ephemDir, ephemYear);
    }
  else
    {
      snprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", ephemYear);
      snprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  ephemYear);
    }

  EphemEarth[FNAME_LENGTH-1] = 0;
  EphemSun[FNAME_LENGTH-1] = 0;

  XLAL_CHECK_NULL ( ( edat = XLALInitBarycenter ( EphemEarth, EphemSun )) != NULL, XLAL_EFUNC );

  return edat;

} /* InitEphemeris() */
