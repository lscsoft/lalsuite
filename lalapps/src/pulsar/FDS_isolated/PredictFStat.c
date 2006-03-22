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
 * Calculate the F-statistic Semi-Analytically of pulsar GW signals.
 * Implements the so-called "F-statistic" as introduced in \ref JKS98 and Cutler-Schutz 2005.
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

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define SQ(x) ((x)*(x))

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  LIGOTimeGPS startTime;		    /**< start time of observation */
  REAL8 duration;			    /**< total time-span of the data (all streams) in seconds */
  EphemerisData *edat;			    /**< ephemeris data (from LALInitBarycenter()) */
  DopplerRegion searchRegion;		    /**< parameter-space region to search over (FIXME) */
  UINT4 numDetectors;			    /**< number of detectors */
  MultiSFTVector *multiSFTs;		    /**< multi-IFO SFT-vectors */
  MultiDetectorStateSeries *multiDetStates; /**< pos, vel and LMSTs for detector at times t_i */
} ConfigVariables;

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

ConfigVariables GV;		/**< global container for various derived configuration settings */

/* ----- User-variables: can be set from config-file or command-line */
BOOLEAN uvar_help;

INT4 uvar_gpsStart;
INT4 uvar_RngMedWindow;

REAL8 uvar_aPlus;
REAL8 uvar_aCross;
REAL8 uvar_phi;
REAL8 uvar_psi;
REAL8 uvar_h0;
REAL8 uvar_cosiota;
REAL8 uvar_Freq;
REAL8 uvar_FreqBand;
REAL8 uvar_Alpha;
REAL8 uvar_Delta;

CHAR *uvar_IFO;
CHAR *uvar_ephemDir;
CHAR *uvar_ephemYear;
CHAR *uvar_DataFiles;
CHAR *uvar_outputLabel;
CHAR *uvar_outputFstat;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

void initUserVars (LALStatus *);
void InitFStat ( LALStatus *, ConfigVariables *cfg );
void Freemem(LALStatus *,  ConfigVariables *cfg);
void WriteFStatLog (LALStatus *, CHAR *argv[]);
void checkUserInputConsistency (LALStatus *);
void InitEphemeris (LALStatus *, EphemerisData *edat, const CHAR *ephemDir, const CHAR *ephemYear, LIGOTimeGPS epoch);

/*---------- empty initializers ---------- */
static const SFTConstraints empty_SFTConstraints;

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/** 
 * MAIN function of PredictFStat code.
 * Calculates the F-statistic for a given position in the sky and detector 
 * semi-analytically and outputs the final 2F value.
 */
int main(int argc,char *argv[]) 
{
  LALStatus status = blank_status;	/* initialize status */

  FILE *fpFstat = NULL;

  CWParamSpacePoint psPoint;		/* parameter-space point to compute Fstat for */
  SkyPosition thisPoint;   

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
  /* like ephemeries data*/
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
  
   /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
  thisPoint.longitude = uvar_Alpha;
  thisPoint.latitude = uvar_Delta;
  thisPoint.system = COORDINATESYSTEM_EQUATORIAL;
  LAL_CALL (LALNormalizeSkyPosition(&status, &thisPoint, &thisPoint), &status);
  
  /* set parameter-space point: sky-position */
  psPoint.skypos = thisPoint;

  { /* Calculating the F-Statistic */
    
    REAL8 F = 0.0;
    REAL8 At = 0.0 ,Bt = 0.0 ,Ct = 0.0; 
    REAL8 A1, A2, A3, A4;
    REAL8 h0, cosi; 
    REAL8 Sh;
    REAL8 aPlus, aCross;
    REAL8 twopsi, phi0;
    REAL8 Tsft;
    UINT4 X;  

    MultiAMCoeffs *multiAMcoef = NULL;
    MultiPSDVector *multiPSDs = NULL;

    LAL_CALL ( LALGetMultiAMCoeffs ( &status, &multiAMcoef, GV.multiDetStates, psPoint.skypos ), &status);
    LAL_CALL ( LALNormalizeMultiSFTVect ( &status, &multiPSDs, GV.multiSFTs, uvar_RngMedWindow ), &status );

    /* Antenna-patterns and compute A,B,C */
    for ( X=0; X < GV.numDetectors; X ++)
      {
	REAL8 A = 0.0, B = 0.0 ,C = 0.0;
	UINT4 alpha;
	UINT4 numSFTsX = GV.multiSFTs->data[X]->length;
	AMCoeffs *amcoeX = multiAMcoef->data[X];
	PSDVector *psdsX = multiPSDs->data[X]; 

	Tsft = 1 / (GV.multiSFTs->data[X]->data[0].deltaF);	

	for(alpha = 0; alpha < numSFTsX; alpha++)
	  {
	    UINT4 lengthPSD, i;
	    REAL8 sumPSD = 0.0, meanPSD = 0.0, PSD;
	    REAL8 ahat, bhat;

	    lengthPSD = psdsX->data[0].data->length;

	    for(i = 0; i < lengthPSD; i++)
	      {
		PSD =  psdsX->data[alpha].data->data[i];
		sumPSD += PSD;  
	      }
	    meanPSD = sumPSD / lengthPSD;

	    Sh = 2 * meanPSD / Tsft;

	    ahat = (amcoeX->a->data[alpha]);
	    bhat = (amcoeX->b->data[alpha]);
	    
	    /* sum A, B, C on the fly */
	    A += ahat * ahat / Sh;
	    B += bhat * bhat / Sh;
	    C += ahat * bhat / Sh;
	    
	  } /* for alpha < numSFTsX */

	At += 2 * Tsft * A;
	Bt += 2 * Tsft * B;
	Ct += 2 * Tsft * C;

      } /* for X < numDetectors */

    /* Free AM Coefficients */
    XLALDestroyMultiAMCoeffs ( multiAMcoef );
    /* Free MultiPSDVector  */
    LAL_CALL ( LALDestroyMultiPSDVector (&status, &multiPSDs ), &status );
      
    phi0 = uvar_phi;
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
    
    A1 = aPlus * cos(twopsi) * cos(phi0) - aCross * sin(twopsi) * sin(phi0);
    A2 = aPlus * sin(twopsi) * cos(phi0) + aCross * cos(twopsi) * sin(phi0);
    A3 =-aPlus * cos(twopsi) * sin(phi0) - aCross * sin(twopsi) * cos(phi0);
    A4 =-aPlus * sin(twopsi) * sin(phi0) + aCross * cos(twopsi) * cos(phi0);
    
    F += (At * ( SQ(A1) + SQ(A3) ) + Bt * ( SQ(A2) + SQ(A4) )+ 2.0 * Ct * (A1 * A2 + A3 * A4 )) / 4;
    
    /* Note: the expectation-value of 2F is 4 + lambda ==> add 2 to Fstat*/
    F += 2.0;
    
/*     fprintf(stdout, "\n2F = %g,   sqrtSh =  %g\n\n", 2 * F , sqrt(Sh)); */
    fprintf(stdout, "\n%g\n", 2 * F);

    if ( fpFstat )
      {
	fprintf (fpFstat, "%g\n%%DONE\n", 2 * F);
	fclose (fpFstat);
	fpFstat = NULL;
      }
    
  }/* Compute FStat */ 
  
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
 uvar_phi       = 0.0;
 uvar_psi       = 0.0;
 uvar_h0        = 0.0;
 uvar_cosiota   = 0.0;
  uvar_Freq      = 100.0; 
 uvar_FreqBand  = 1.0;
 uvar_Alpha 	= 0.0;
 uvar_Delta 	= 0.0;
 uvar_RngMedWindow = 50;	/* for running-median */
 
 uvar_ephemYear = LALCalloc (1, strlen(EPHEM_YEARS)+1);
 strcpy (uvar_ephemYear, EPHEM_YEARS);
 
#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
 uvar_ephemDir = LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
 strcpy (uvar_ephemDir, DEFAULT_EPHEMDIR);

 uvar_help = FALSE;
 uvar_outputLabel = NULL;
 uvar_outputFstat = NULL;
  
 /* register all our user-variables */
 LALregBOOLUserVar(status, 	help, 		'h', UVAR_HELP,     "Print this message"); 

 LALregINTUserVar(status, 	RngMedWindow,	'k', UVAR_DEVELOPER, "Running-Median window size");

 LALregREALUserVar(status,      phi,            'Q', UVAR_OPTIONAL, "Phi_0: Initial phase in radians");
 LALregREALUserVar(status,      psi,            'Y', UVAR_OPTIONAL, "Polarisation in radians");
 LALregREALUserVar(status,      cosiota,        'i', UVAR_OPTIONAL, "Cos(iota)");
 LALregREALUserVar(status,      h0,             's', UVAR_OPTIONAL, "Strain amplitude h_0");
  LALregREALUserVar(status,      aPlus,            0, UVAR_OPTIONAL, "Strain amplitude h_0");
 LALregREALUserVar(status,      aCross,           0, UVAR_OPTIONAL, "Noise floor: one-sided sqrt(Sh) in 1/sqrt(Hz)");
 LALregREALUserVar(status, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
 LALregREALUserVar(status, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
 LALregREALUserVar(status,      Freq, 	        'F', UVAR_OPTIONAL, "Search Frequency");
 LALregREALUserVar(status,      FreqBand, 	'B', UVAR_OPTIONAL, "Search Frequency Band");
 
 LALregSTRINGUserVar(status,	DataFiles, 	'D', UVAR_REQUIRED, "File-pattern specifying (multi-IFO) input SFT-files"); 
 LALregSTRINGUserVar(status, 	IFO, 		'I', UVAR_OPTIONAL, "Detector-constraint: 'G1', 'L1', 'H1', 'H2' ...(useful for single-IFO v1-SFTs only!)");
 LALregSTRINGUserVar(status,	ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
 LALregSTRINGUserVar(status,	ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
 LALregSTRINGUserVar(status,	outputLabel,	'o', UVAR_OPTIONAL, "Label to be appended to all output file-names");
 LALregSTRINGUserVar(status,	outputFstat,	 0,  UVAR_OPTIONAL, "Output-file for F-statistic field over the parameter-space");

 DETATCHSTATUSPTR (status);
 RETURN (status);
} /* initUserVars() */

/** Load Ephemeris from ephemeris data-files  */
void
InitEphemeris (LALStatus * status,   
  EphemerisData *edat,	                /**< [out] the ephemeris-data */
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

/** Initialized Fstat-code: handle user-input and set everything up. */
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

  {/* ----- load the multi-IFO SFT-vectors ----- */
 
    REAL8 fMin = MYMIN ( uvar_Freq, uvar_Freq + uvar_FreqBand );
    REAL8 fMax = MYMAX ( uvar_Freq, uvar_Freq + uvar_FreqBand );

    TRY ( LALLoadMultiSFTs ( status->statusPtr, &(cfg->multiSFTs), catalog, fMin, fMax ), status );
    TRY ( LALDestroySFTCatalog ( status->statusPtr, &catalog ), status );
  }

  cfg->numDetectors = cfg->multiSFTs->length;

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

  /* destroy DetectorStateSeries */
  XLALDestroyMultiDetectorStateSeries ( cfg->multiDetStates );
  
  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (status->statusPtr), status);

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

   RETURN (status);
} /* checkUserInputConsistency() */
