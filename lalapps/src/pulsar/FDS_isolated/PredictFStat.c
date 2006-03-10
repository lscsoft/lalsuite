/*
 * Copyright (C) 2006 Iraj Gholami, Reinhard Prix
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

/*********************************************************************************/
/** \author I. Gholami, R. Prix
 * \file 
 * \brief
 * Calculate the F-statistic for a given parameter-space of pulsar GW signals.
 * Implements the so-called "F-statistic" as introduced in \ref JKS98.
 *                                                                          
 *********************************************************************************/
#include "config.h"

/* System includes */
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* LAL-includes */
#include <lal/AVFactories.h>

#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/ExtrapolatePulsarSpins.h>

#include <lal/NormalizeSFTRngMed.h>
#include <lal/ComputeFstat.h>
#include <lal/LALHough.h>

#include <lalapps.h>

/* local includes */
#include "DopplerScan.h"

RCSID( "$Id$");

/*---------- DEFINES ----------*/

#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

#define EPHEM_YEARS  "00-04"	/**< default range: override with --ephemYear */

#define TRUE (1==1)
#define FALSE (1==0)

/*----- SWITCHES -----*/
#define NUM_SPINS 4		/* number of spin-values to consider: {f, fdot, f2dot, ... } */ 


/*----- Error-codes -----*/
#define COMPUTEFSTATISTIC_ENULL 	1
#define COMPUTEFSTATISTIC_ESYS     	2
#define COMPUTEFSTATISTIC_EINPUT   	3
#define COMPUTEFSTATISTIC_EMEM   	4
#define COMPUTEFSTATISTIC_ENONULL 	5
#define COMPUTEFSTATISTIC_EXLAL		6

#define COMPUTEFSTATISTIC_MSGENULL 	"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATISTIC_MSGESYS	"System call failed (probably file IO)"
#define COMPUTEFSTATISTIC_MSGEINPUT   	"Invalid input"
#define COMPUTEFSTATISTIC_MSGEMEM   	"Out of memory. Bad."
#define COMPUTEFSTATISTIC_MSGENONULL 	"Output pointer is non-NULL"
#define COMPUTEFSTATISTIC_MSGEXLAL	"XLALFunction-call failed"

/*----- Macros -----*/

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

/*---------- internal types ----------*/

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  LIGOTimeGPS startTime;		    /**< start time of observation */
  REAL8 duration;			    /**< total time-span of the data (all streams) in seconds */
  LIGOTimeGPS refTime;			    /**< reference-time for pulsar-parameters in SBB frame */
  EphemerisData *edat;			    /**< ephemeris data (from LALInitBarycenter()) */
  DopplerRegion searchRegion;		    /**< parameter-space region to search over (FIXME) */
  UINT4 numDetectors;			    /**< number of detectors */
  MultiSFTVector *multiSFTs;		    /**< multi-IFO SFT-vectors */
  MultiDetectorStateSeries *multiDetStates; /**< pos, vel and LMSTs for detector at times t_i */
  MultiNoiseWeights *multiNoiseWeights;	    /**< normalized noise-weights of those SFTs */
  ComputeFParams CFparams;		    /**< parameters for the computation of Fstat (e.g Dterms, SSB-precision,...) */
} ConfigVariables;

#define SQ(x) ((x)*(x))

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

ConfigVariables GV;		/**< global container for various derived configuration settings */

/* ----- User-variables: can be set from config-file or command-line */
REAL8 uvar_aPlus;
REAL8 uvar_aCross;
INT4 uvar_gpsStart;
REAL8 uvar_phi;
REAL8 uvar_psi;
REAL8 uvar_h0;
REAL8 uvar_cosiota;
REAL8 uvar_sqrtSh;

REAL8 uvar_fmin;
REAL8 uvar_fmax;
CHAR *uvar_IFO;
BOOLEAN uvar_SignalOnly;
REAL8 uvar_Freq;
REAL8 uvar_FreqBand;
REAL8 uvar_Alpha;
REAL8 uvar_dAlpha;
REAL8 uvar_AlphaBand;
REAL8 uvar_Delta;
REAL8 uvar_dDelta;
REAL8 uvar_DeltaBand;
/* --- */
REAL8 uvar_Fthreshold;
CHAR *uvar_ephemDir;
CHAR *uvar_ephemYear;
CHAR *uvar_skyRegion;
CHAR *uvar_DataFiles;
BOOLEAN uvar_help;
CHAR *uvar_outputLabel;
CHAR *uvar_outputFstat;
CHAR *uvar_outputLoudest;
CHAR *uvar_skyGridFile;
CHAR *uvar_outputSkyGrid;
CHAR *uvar_workingDir;
INT4 uvar_RngMedWindow;
REAL8 uvar_refTime;
INT4 uvar_SSBprecision;
REAL8 uvar_dopplermax;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);
void initUserVars (LALStatus *);
void InitFStat ( LALStatus *, ConfigVariables *cfg );

void Freemem(LALStatus *,  ConfigVariables *cfg);

void WriteFStatLog (LALStatus *, CHAR *argv[]);
void checkUserInputConsistency (LALStatus *);
int outputBeamTS( const CHAR *fname, const AMCoeffs *amcoe, const DetectorStateSeries *detStates );
void InitEphemeris (LALStatus *, EphemerisData *edat, const CHAR *ephemDir, const CHAR *ephemYear, LIGOTimeGPS epoch);
void getUnitWeights ( LALStatus *, MultiNoiseWeights **multiWeights, const MultiSFTVector *multiSFTs );

const char *va(const char *format, ...);	/* little var-arg string helper function */

/*---------- empty initializers ---------- */
static const PulsarTimesParamStruc empty_PulsarTimesParamStruc;
static const BarycenterInput empty_BarycenterInput;
static const SFTConstraints empty_SFTConstraints;
static const ComputeFBuffer empty_ComputeFBuffer;

/*----------------------------------------------------------------------*/
/* Function definitions start here */
/*----------------------------------------------------------------------*/

/** 
 * MAIN function of ComputeFStatistic code.
 * Calculate the F-statistic over a given portion of the parameter-space
 * and write a list of 'candidates' into a file(default: 'Fstats').
 */
int main(int argc,char *argv[]) 
{
  LALStatus status = blank_status;	/* initialize status */

  SkyPosition thisPoint;
  
  FILE *fpFstat = NULL;
  CWParamSpacePoint psPoint;		/* parameter-space point to compute Fstat for */
   

  lalDebugLevel = 0;  
  vrbflg = 1;	/* verbose error-messages */
  
  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LAL_CALL (initUserVars(&status), &status); 	

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput(&status, argc,argv), &status);	

  if (uvar_help)	/* if help was requested, we're done here */
    exit (0);

  /* keep a log-file recording all relevant parameters of this search-run */
  LAL_CALL (WriteFStatLog (&status, argv), &status);

  /* do some sanity checks on the user-input before we proceed */
  LAL_CALL ( checkUserInputConsistency(&status), &status);

  /* Initialization the common variables of the code, */
  /* like ephemeries data and template grids: */
  LAL_CALL ( InitFStat(&status, &GV), &status);

  
  /* if a complete output of the F-statistic file was requested,
   * we open and prepare the output-file here */
  if (uvar_outputFstat)
    {
      if ( (fpFstat = fopen (uvar_outputFstat, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputFstat);
	  return (COMPUTEFSTATISTIC_ESYS);
	}
    } /* if outputFstat */
  
  if (lalDebugLevel) printf ("\nStarting main search-loop.. \n");
  
  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
   */
  
   /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
  thisPoint.longitude = uvar_Alpha;
  thisPoint.latitude = uvar_Delta;
  thisPoint.system = COORDINATESYSTEM_EQUATORIAL;
  LAL_CALL (LALNormalizeSkyPosition(&status, &thisPoint, &thisPoint), &status);
  
  /* set parameter-space point: sky-position */
  psPoint.skypos = thisPoint;
  
  { /* Calculating the F-Statistic */
    
    REAL8 A = 0.0, B = 0.0 ,C = 0.0, At = 0.0 ,Bt = 0.0 ,Ct = 0.0 ,Dt = 0.0; 
    REAL8 A1, A2, A3, A4, h0, cosi, F = 0.0;
    REAL8 aPlus, aCross;
    REAL8 twopsi, twophi;
    REAL8 Tsft;
    UINT4 X;  

    MultiAMCoeffs *multiAMcoef = NULL;
    MultiSSBtimes *multiSSB = NULL;

    LAL_CALL ( LALGetMultiSSBtimes ( &status, &multiSSB, GV.multiDetStates, psPoint.skypos, psPoint.refTime, GV.CFparams.SSBprec ), &status );
    LAL_CALL ( LALGetMultiAMCoeffs ( &status, &multiAMcoef, GV.multiDetStates, psPoint.skypos ), &status);

    XLALDestroyMultiSSBtimes ( multiSSB );

    /* noise-weight Antenna-patterns and compute A,B,C */
    for ( X=0; X < GV.numDetectors; X ++)
      {
	UINT4 alpha;
	UINT4 numSFTsX = GV.multiSFTs->data[X]->length;
	AMCoeffs *amcoeX = multiAMcoef->data[X];
	REAL8Vector *weightsX = GV.multiNoiseWeights->data[X];

	Tsft = 1 / (GV.multiSFTs->data[X]->data[0].deltaF);	

	for(alpha = 0; alpha < numSFTsX; alpha++)
	  {
	    REAL8 SqWN = sqrt ( weightsX->data[alpha] );
	    REAL8 ahat = SqWN * amcoeX->a->data[alpha] ;
	    REAL8 bhat = SqWN * amcoeX->b->data[alpha] ;
	    
	    /* sum A, B, C on the fly */
	    A += ahat * ahat;
	    B += bhat * bhat;
	    C += ahat * bhat;
	    
	  } /* for alpha < numSFTsX */

	At += (Tsft ) * A;
	Bt += (Tsft ) * B;
	Ct += (Tsft ) * C;
	Dt += A * B - C * C; 
	
      } /* for X < numDetectors */
        
    twophi = 2.0 * uvar_phi;
    twopsi = 2.0 * uvar_psi;
    
    h0 = uvar_h0;
    cosi = uvar_cosiota;
    
    if ( h0 != 0 ) 
      {
	aPlus = h0 * 0.5 * (1.0 + cosi*cosi );
	aCross= h0 * cosi;
      } 
    else   /* alternative way to specify amplitude (compatible with mfd_v4) */
      {
	aPlus = uvar_aPlus;
	aCross = uvar_aCross;
      }
    
    A1 = aPlus * cos(twopsi) * cos(twophi) - aCross * sin(twopsi) * sin(twophi);
    A2 = aPlus * sin(twopsi) * cos(twophi) + aCross * cos(twopsi) * sin(twophi);
    A3 =-aPlus * cos(twopsi) * sin(twophi) - aCross * sin(twopsi) * cos(twophi);
    A4 =-aPlus * sin(twopsi) * sin(twophi) + aCross * cos(twopsi) * cos(twophi);
    
    F += At * ( SQ(A1) + SQ(A3) ) + Bt * ( SQ(A2) + SQ(A4) )+ 2.0 * Ct * (A1 * A2 + A3 * A4 );
    
    /* Note: the expectation-value of 2F is 4 + lambda ==> add 2 to Fstat*/
    F += 2.0;
    
    fprintf( stdout, "%g\n", F );
    
  }/* Compute FStat */ 
  
  
  if ( fpFstat )
    {
      fprintf (fpFstat, "%%DONE\n");
      fclose (fpFstat);
      fpFstat = NULL;
    }
  
  if (lalDebugLevel) printf ("\nSearch finished.\n");
  
  LAL_CALL ( Freemem(&status, &GV), &status);
  
  /* did we forget anything ? */
  LALCheckMemoryLeaks();
  
  return 0;
  
} /* main() */


/** 
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
void
initUserVars (LALStatus *status)
{
INITSTATUS( status, "initUserVars", rcsid );
ATTATCHSTATUSPTR (status);

/* set a few defaults */
 uvar_phi = 0.0;
 uvar_psi = 0.0;
 uvar_h0 = 0.0;
 uvar_cosiota = 0.0;
 uvar_sqrtSh = 1.0;
 
 uvar_FreqBand = 0.0;
 uvar_Alpha 	= 0.0;
 uvar_Delta 	= 0.0;
 uvar_AlphaBand = 0;
 uvar_DeltaBand = 0;
 uvar_dAlpha 	= 0.001;
 uvar_dDelta 	= 0.001;
 uvar_Freq       = 100.0;
 uvar_fmin       = 99.0;
 uvar_fmax       = 101.0;
 uvar_skyRegion = NULL;
 
 uvar_ephemYear = LALCalloc (1, strlen(EPHEM_YEARS)+1);
 strcpy (uvar_ephemYear, EPHEM_YEARS);
 
#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
 uvar_ephemDir = LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
 strcpy (uvar_ephemDir, DEFAULT_EPHEMDIR);
 
 uvar_SignalOnly = FALSE;
 
 uvar_Fthreshold = 10.0;
 
 uvar_help = FALSE;
 uvar_outputLabel = NULL;
 
 uvar_outputFstat = NULL;
  
 uvar_skyGridFile = NULL;
 
 uvar_workingDir = LALMalloc(512);
 strcpy(uvar_workingDir, ".");
 
 uvar_RngMedWindow = 50;	/* for running-median */
 
 uvar_SSBprecision = SSBPREC_RELATIVISTIC;
 
 /* register all our user-variables */
 LALregBOOLUserVar(status, 	help, 		'h', UVAR_HELP,     "Print this message"); 
 
 LALregREALUserVar(status,      phi,            'Q', UVAR_OPTIONAL, "Phi_0: Initial phase in radians");
 LALregREALUserVar(status,      psi,            'Y', UVAR_OPTIONAL, "Polarisation in radians");
 LALregREALUserVar(status,      cosiota,        'i', UVAR_OPTIONAL, "Cos(iota)");
 LALregREALUserVar(status,      h0,             's', UVAR_OPTIONAL, "Strain amplitude h_0");
 LALregREALUserVar(status,      sqrtSh,         'N', UVAR_OPTIONAL, "Noise floor: one-sided sqrt(Sh) in 1/sqrt(Hz)");
 LALregREALUserVar(status,      aPlus,            0, UVAR_OPTIONAL, "Strain amplitude h_0");
 LALregREALUserVar(status,      aCross,           0, UVAR_OPTIONAL, "Noise floor: one-sided sqrt(Sh) in 1/sqrt(Hz)");
	   
 LALregREALUserVar(status, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
 LALregREALUserVar(status, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
 LALregREALUserVar(status, 	AlphaBand, 	'z', UVAR_OPTIONAL, "Band in alpha (equatorial coordinates) in radians");
 LALregREALUserVar(status, 	DeltaBand, 	'c', UVAR_OPTIONAL, "Band in delta (equatorial coordinates) in radians");
 LALregREALUserVar(status, 	dAlpha, 	'l', UVAR_OPTIONAL, "Resolution in alpha (equatorial coordinates) in radians");
 LALregREALUserVar(status, 	dDelta, 	'g', UVAR_OPTIONAL, "Resolution in delta (equatorial coordinates) in radians");
 LALregREALUserVar(status,      Freq, 	        'F', UVAR_OPTIONAL, "Search Frequency");
 LALregREALUserVar(status,      FreqBand, 	'B', UVAR_OPTIONAL, "Search Frequency Band");
 LALregREALUserVar(status,      fmin,   	  0, UVAR_OPTIONAL, "Minimum Search Frequency");
 LALregREALUserVar(status,      fmax,   	  0, UVAR_OPTIONAL, "Maximum Search Frequency Band");
 
 LALregSTRINGUserVar(status,	skyRegion, 	'R', UVAR_OPTIONAL, "ALTERNATIVE: Specify sky-region by polygon (or use 'allsky')");
 LALregSTRINGUserVar(status,	DataFiles, 	'D', UVAR_REQUIRED, "File-pattern specifying (multi-IFO) input SFT-files"); 
 LALregSTRINGUserVar(status, 	IFO, 		'I', UVAR_OPTIONAL, 
		     "Detector-constraint: 'G1', 'L1', 'H1', 'H2' ...(useful for single-IFO v1-SFTs only!)");
 LALregSTRINGUserVar(status,	ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
 LALregSTRINGUserVar(status,	ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
 LALregBOOLUserVar(status, 	SignalOnly, 	'S', UVAR_OPTIONAL, "Signal only flag");
 LALregREALUserVar(status, 	Fthreshold,	'F', UVAR_OPTIONAL, "Signal Set the threshold for selection of 2F");
 LALregSTRINGUserVar(status,	outputLabel,	'o', UVAR_OPTIONAL, "Label to be appended to all output file-names");
 LALregSTRINGUserVar(status,	skyGridFile,	 0,  UVAR_OPTIONAL, "Load sky-grid from this file.");
 LALregREALUserVar(status,	refTime,	 0,  UVAR_OPTIONAL, "SSB reference time for pulsar-paramters");
 LALregSTRINGUserVar(status,	outputFstat,	 0,  UVAR_OPTIONAL, "Output-file for F-statistic field over the parameter-space");
  
 /* more experimental and unofficial stuff follows here */
 LALregINTUserVar (status, 	SSBprecision,	 0,  UVAR_DEVELOPER, "Precision to use for time-transformation to SSB: 0=Newtonian 1=relativistic");
 LALregINTUserVar(status, 	RngMedWindow,	'k', UVAR_DEVELOPER, "Running-Median window size");
 LALregSTRINGUserVar(status,    workingDir,     'w', UVAR_DEVELOPER, "Directory to be made the working directory, . is default");
 LALregSTRINGUserVar(status,	outputSkyGrid,	 0,  UVAR_DEVELOPER, "Write sky-grid into this file.");
 LALregSTRINGUserVar(status,    outputLoudest,	 0,  UVAR_DEVELOPER, 
		     "Output-file for the loudest F-statistic candidate in this search");
 
 
 DETATCHSTATUSPTR (status);
 RETURN (status);
} /* initUserVars() */

/** Load Ephemeris from ephemeris data-files  */
void
InitEphemeris (LALStatus * status,   
  EphemerisData *edat,	/**< [out] the ephemeris-data */
	       const CHAR *ephemDir,	/**< directory containing ephems */
	       const CHAR *ephemYear,	/**< which years do we need? */
	       LIGOTimeGPS epoch	/**< epoch of observation */
	       )
{
#define FNAME_LENGTH 1024
  CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4 leap;

  INITSTATUS( status, "InitEphemeris", rcsid );
  ATTATCHSTATUSPTR (status);
  
  ASSERT ( edat, status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL );
  ASSERT ( ephemYear, status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL );
  
  if ( ephemDir )
    {
      LALSnprintf(EphemEarth, FNAME_LENGTH, "%s/earth%s.dat", ephemDir, ephemYear);
      LALSnprintf(EphemSun, FNAME_LENGTH, "%s/sun%s.dat", ephemDir, ephemYear);
    }
  else
    {
      LALSnprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", ephemYear);
      LALSnprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  ephemYear);
    }
  EphemEarth[FNAME_LENGTH-1]=0;
  EphemSun[FNAME_LENGTH-1]=0;
  
  /* NOTE: the 'ephiles' are ONLY ever used in LALInitBarycenter, which is
   * why we can use local variables (EphemEarth, EphemSun) to initialize them.
   */
  edat->ephiles.earthEphemeris = EphemEarth;
  edat->ephiles.sunEphemeris = EphemSun;
  
  TRY (LALLeapSecs (status->statusPtr, &leap, &epoch, &formatAndAcc), status);
  edat->leap = (INT2) leap;
  
  TRY (LALInitBarycenter(status->statusPtr, edat), status);
  
  DETATCHSTATUSPTR ( status );
  RETURN ( status );
  
} /* InitEphemeris() */



/** Initialized Fstat-code: handle user-input and set everything up.
 * NOTE: the logical *order* of things in here is very important, so be careful
 */
void
InitFStat ( LALStatus *status, ConfigVariables *cfg )
{
  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_SFTConstraints;

  UINT4 numSFTs;
  REAL8 Tsft;
  LIGOTimeGPS endTime;

  INITSTATUS (status, "InitFStat", rcsid);
  ATTATCHSTATUSPTR (status);

  /* use IFO-contraint if one given by the user */
  if ( LALUserVarWasSet ( &uvar_IFO ) )
    if ( (constraints.detector = XLALGetChannelPrefix ( uvar_IFO )) == NULL ) {
      ABORT ( status,  COMPUTEFSTATISTIC_EINPUT,  COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* get full SFT-catalog of all matching (multi-IFO) SFTs */
  TRY ( LALSFTdataFind ( status->statusPtr, &catalog, uvar_DataFiles, &constraints ), status);    
  if ( constraints.detector ) 
    LALFree ( constraints.detector );

  if ( catalog->length == 0 ) 
    {
      LALPrintError ("\nSorry, didn't find any matching SFTs with pattern '%s'!\n\n", uvar_DataFiles );
      ABORT ( status,  COMPUTEFSTATISTIC_EINPUT,  COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* deduce start- and end-time of the observation spanned by the data */
  numSFTs = catalog->length;
  Tsft = 1.0 / catalog->data[0].header.deltaF;
  cfg->startTime = catalog->data[0].header.epoch;
  endTime   = catalog->data[numSFTs-1].header.epoch;
  LALAddFloatToGPS(status->statusPtr, &endTime, &endTime, Tsft );	/* can't fail */
  cfg->duration = GPS2REAL8(endTime) - GPS2REAL8 (cfg->startTime);

  /* ----- get reference-time (from user if given, use startTime otherwise): ----- */
  if ( LALUserVarWasSet(&uvar_refTime)) {
    TRY ( LALFloatToGPS (status->statusPtr, &(cfg->refTime), &uvar_refTime), status);
  } else
    cfg->refTime = cfg->startTime;


  { /* ----- get sky-region to search ----- */
    BOOLEAN haveAlphaDelta = LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta);
    if (uvar_skyRegion)
      {
	cfg->searchRegion.skyRegionString = (CHAR*)LALCalloc(1, strlen(uvar_skyRegion)+1);
	if ( cfg->searchRegion.skyRegionString == NULL ) {
	  ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
	}
	strcpy (cfg->searchRegion.skyRegionString, uvar_skyRegion);
      }
    else if (haveAlphaDelta)    /* parse this into a sky-region */
      {
	REAL8 eps = 1e-9;	/* hack for backwards compatbility */
	TRY ( SkySquare2String( status->statusPtr, &(cfg->searchRegion.skyRegionString),
				uvar_Alpha, uvar_Delta,
				uvar_AlphaBand + eps, uvar_DeltaBand + eps), status);
      }
  } /* get sky-region */

  {/* ----- load the multi-IFO SFT-vectors ----- */
 
    REAL8 fMax = uvar_fmax;
    REAL8 fMin = uvar_fmin;

    TRY ( LALLoadMultiSFTs ( status->statusPtr, &(cfg->multiSFTs), catalog, fMin, fMax ), status );
    TRY ( LALDestroySFTCatalog ( status->statusPtr, &catalog ), status );
  }

  cfg->numDetectors = cfg->multiSFTs->length;
  
  /* ----- normalize SFTs and calculate noise-weights ----- */
  if ( uvar_SignalOnly )
    {
      TRY ( getUnitWeights ( status->statusPtr, &(cfg->multiNoiseWeights), cfg->multiSFTs ), status );
    }
  else
    {
      MultiPSDVector *psds = NULL;
      REAL8 S_hat;

       TRY ( LALNormalizeMultiSFTVect (status->statusPtr, &psds, cfg->multiSFTs, uvar_RngMedWindow ), status ); 
      /* note: the normalization S_hat would be required to compute the ML-estimator for A^\mu */
       TRY ( LALComputeMultiNoiseWeights  (status->statusPtr, &(cfg->multiNoiseWeights), &S_hat, psds, uvar_RngMedWindow ), status );
   
      TRY ( LALDestroyMultiPSDVector (status->statusPtr, &psds ), status );
      
    } /* if ! SignalOnly */

  { /* ----- load ephemeris-data ----- */
    CHAR *ephemDir;

    cfg->edat = LALCalloc(1, sizeof(EphemerisData));
    if ( LALUserVarWasSet ( &uvar_ephemDir ) )
      ephemDir = uvar_ephemDir;
    else
      ephemDir = NULL;
    TRY( InitEphemeris (status->statusPtr, cfg->edat, ephemDir, uvar_ephemYear, cfg->startTime ),status);
  }
  
  /* ----- obtain the (multi-IFO) 'detector-state series' for all SFTs ----- */
  TRY ( LALGetMultiDetectorStates ( status->statusPtr, &(cfg->multiDetStates), cfg->multiSFTs, cfg->edat ), status );

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* InitFStat() */


/***********************************************************************/
/** Log the all relevant parameters of the present search-run to a log-file.
 * The name of the log-file is "Fstats{uvar_outputLabel}.log".
 * <em>NOTE:</em> Currently this function only logs the user-input and code-versions.
 */
void
WriteFStatLog (LALStatus *status, char *argv[])
{
  CHAR *logstr = NULL;
  const CHAR *head = "Fstats";
  CHAR command[512] = "";
  UINT4 len;
  CHAR *fname = NULL;
  FILE *fplog;

  INITSTATUS (status, "WriteFStatLog", rcsid);
  ATTATCHSTATUSPTR (status);

  /* prepare log-file for writing */
  len = strlen(head) + strlen(".log") +10;
  if (uvar_outputLabel)
    len += strlen(uvar_outputLabel);

  if ( (fname=LALCalloc(len,1)) == NULL) {
    ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM);
  }
  strcpy (fname, head);
  if (uvar_outputLabel)
    strcat (fname, uvar_outputLabel);
  strcat (fname, ".log");

  if ( (fplog = fopen(fname, "w" )) == NULL) {
    LALPrintError ("\nFailed to open log-file '%f' for writing.\n\n", fname);
    LALFree (fname);
    ABORT (status, COMPUTEFSTATISTIC_ESYS, COMPUTEFSTATISTIC_MSGESYS);
  }

  /* write out a log describing the complete user-input (in cfg-file format) */
  TRY (LALUserVarGetLog (status->statusPtr, &logstr,  UVAR_LOGFMT_CFGFILE), status);

  fprintf (fplog, "## LOG-FILE of ComputeFStatistic run\n\n");
  fprintf (fplog, "# User-input:\n");
  fprintf (fplog, "# ----------------------------------------------------------------------\n\n");

  fprintf (fplog, logstr);
  LALFree (logstr);

  /* append an ident-string defining the exact CVS-version of the code used */
  fprintf (fplog, "\n\n# CVS-versions of executable:\n");
  fprintf (fplog, "# ----------------------------------------------------------------------\n");
  fclose (fplog);
    
  sprintf (command, "ident %s 2> /dev/null | sort -u >> %s", argv[0], fname);
  system (command);	/* we don't check this. If it fails, we assume that */
    			/* one of the system-commands was not available, and */
    			/* therefore the CVS-versions will not be logged */

  LALFree (fname);

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* WriteFStatLog() */


/** Free all globally allocated memory. */
void
Freemem(LALStatus *status,  ConfigVariables *cfg) 
{
  INITSTATUS (status, "Freemem", rcsid);
  ATTATCHSTATUSPTR (status);


  /* Free SFT data */
  TRY ( LALDestroyMultiSFTVector (status->statusPtr, &(cfg->multiSFTs) ), status );
  /* and corresponding noise-weights */
  TRY ( LALDestroyMultiNoiseWeights (status->statusPtr, &(cfg->multiNoiseWeights) ), status );

  /* destroy DetectorStateSeries */
  XLALDestroyMultiDetectorStateSeries ( cfg->multiDetStates );


  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (status->statusPtr), status);

  if ( GV.searchRegion.skyRegionString )
    LALFree ( GV.searchRegion.skyRegionString );
  
  /* Free ephemeris data */
  LALFree(cfg->edat->ephemE);
  LALFree(cfg->edat->ephemS);
  LALFree(cfg->edat);

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* Freemem() */


/*----------------------------------------------------------------------*/
/** Some general consistency-checks on user-input.
 * Throws an error plus prints error-message if problems are found.
 */
void
checkUserInputConsistency (LALStatus *status)
{

  INITSTATUS (status, "checkUserInputConsistency", rcsid);  
  
  if (uvar_ephemYear == NULL)
    {
      LALPrintError ("\nNo ephemeris year specified (option 'ephemYear')\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }      

  /* don't allow negative bands (for safty in griding-routines) */
  if ( (uvar_AlphaBand < 0) ||  (uvar_DeltaBand < 0) )
    {
      LALPrintError ("\nNegative value of sky-bands not allowed (alpha or delta)!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  /* check for negative stepsizes in Freq, Alpha, Delta */
  if ( LALUserVarWasSet(&uvar_dAlpha) && (uvar_dAlpha < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dAlpha not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar_dDelta) && (uvar_dDelta < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dDelta not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* grid-related checks */
  {
    BOOLEAN haveAlphaBand = LALUserVarWasSet( &uvar_AlphaBand );
    BOOLEAN haveDeltaBand = LALUserVarWasSet( &uvar_DeltaBand );
    BOOLEAN haveSkyRegion, haveAlphaDelta;

    haveSkyRegion  = (uvar_skyRegion != NULL);
    haveAlphaDelta = (LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta) );

    if ( (haveAlphaBand && !haveDeltaBand) || (haveDeltaBand && !haveAlphaBand) )
      {
	LALPrintError ("\nERROR: Need either BOTH (AlphaBand, DeltaBand) or NONE.\n\n"); 
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

    if ( haveSkyRegion && haveAlphaDelta )
      {
        LALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta)"
		       " OR skyRegion!\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
  } /* Grid-related checks */

  RETURN (status);
} /* checkUserInputConsistency() */

/* debug-output a(t) and b(t) into given file.
 * return 0 = OK, -1 on error
 */
int
outputBeamTS( const CHAR *fname, const AMCoeffs *amcoe, const DetectorStateSeries *detStates )
{
  FILE *fp;
  UINT4 i, len;

  if ( !fname || !amcoe || !amcoe->a || !amcoe->b || !detStates)
    return -1;
  
  len = amcoe->a->length;
  if ( (len != amcoe->b->length) || ( len != detStates->length ) )
    return -1;

  if ( (fp = fopen(fname, "wb")) == NULL )
    return -1;

  for (i=0; i < len; i ++ )
    {
      INT4 ret;
      ret = fprintf (fp, "%9d %f %f %f \n", 
		     detStates->data[i].tGPS.gpsSeconds, detStates->data[i].LMST, amcoe->a->data[i], amcoe->b->data[i] );
      if ( ret < 0 )
	{
	  fprintf (fp, "ERROR\n");
	  fclose(fp);
	  return -1;
	}
    }

  fclose(fp);
  return 0;
} /* outputBeamTS() */

/** Helper-function to generate unity noise-weights for given multiSFT-vector */
void
getUnitWeights ( LALStatus *status, MultiNoiseWeights **multiWeights, const MultiSFTVector *multiSFTs )
{
  MultiNoiseWeights *ret = NULL;
  UINT4 X, numDet;

  INITSTATUS (status, "getUnitWeights", rcsid);
  ATTATCHSTATUSPTR (status); 

  ASSERT ( multiSFTs, status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL);
  ASSERT ( multiSFTs->data, status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL);
  ASSERT ( multiSFTs->length, status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);

  ASSERT ( multiWeights, status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL);
  ASSERT ( (*multiWeights) == NULL, status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL);

  numDet = multiSFTs->length;

  if ( (ret = LALCalloc(1, sizeof(MultiNoiseWeights))) == NULL ){
    ABORT (status,  COMPUTEFSTATISTIC_EMEM,  COMPUTEFSTATISTIC_MSGEMEM);
  }
  ret->length = numDet;

  if ( (ret->data = LALCalloc( numDet, sizeof(REAL8Vector *))) == NULL) {
    ABORT (status,  COMPUTEFSTATISTIC_EMEM,  COMPUTEFSTATISTIC_MSGEMEM);      
  }

  for ( X = 0; X < numDet; X ++) 
    {
      UINT4 numsfts = multiSFTs->data[X]->length;
      UINT4 alpha;

      /* create k^th weights vector */
      if ( (ret->data[X] = XLALCreateREAL8Vector ( numsfts )) == NULL ) {
	ABORT (status,  COMPUTEFSTATISTIC_EMEM,  COMPUTEFSTATISTIC_MSGEMEM);      
      }

      /* set all weights to unity */
      for ( alpha = 0; alpha < numsfts; alpha ++ ) 
	ret->data[X]->data[alpha] = 1.0;
      
    } /* for X < numDet */


  (*multiWeights) = ret;
  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* getUnitWeights() */


/*
  ============
  va ['stolen' from Quake2 (GPL'ed)]

  does a varargs printf into a temp buffer, so I don't need to have
  varargs versions of all text functions.
  FIXME: make this buffer size safe someday
  ============
*/
const char *va(const char *format, ...)
{
  va_list         argptr;
  static char     string[1024];

  va_start (argptr, format);
  vsprintf (string, format,argptr);
  va_end (argptr);

  return string;
}
