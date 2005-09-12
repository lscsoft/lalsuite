/*
 * Copyright (C) 2005 Reinhard Prix, Iraj Gholami
 * Copyright (C) 2004 Reinhard Prix
 * Copyright (C) 2002, 2003, 2004 M.A. Papa, X. Siemens, Y. Itoh
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
/** \author R. Prix, Y. Ioth, Papa, X. Siemens 
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
 

/*----- Error-codes -----*/
#define COMPUTEFSTATISTICC_ENULL 		1
#define COMPUTEFSTATISTICC_ESYS     		2
#define COMPUTEFSTATISTICC_EINPUT   		3
#define COMPUTEFSTATISTICC_EMEM   		4
#define COMPUTEFSTATISTICC_ENONULL 		5
#define COMPUTEFSTATISTICC_EXLAL		6

#define COMPUTEFSTATISTICC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATISTICC_MSGESYS		"System call failed (probably file IO)"
#define COMPUTEFSTATISTICC_MSGEINPUT   		"Invalid input"
#define COMPUTEFSTATISTICC_MSGEMEM   		"Out of memory. Bad."
#define COMPUTEFSTATISTICC_MSGENONULL 		"Output pointer is non-NULL"
#define COMPUTEFSTATISTICC_MSGEXLAL		"XLALFunction-call failed"

/*----- Macros -----*/

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

/*---------- internal types ----------*/

/** Detectors Vector; contains all quantities that are 'detector-specific', i.e for which we need an entry per detector */
typedef struct {
  UINT4 length;
  LALDetector *Detectors;         
  SFTVector **sftVects;
  DetectorStateSeries **DetectorStates;	/**< pos, vel and LMSTs for detector at times t_i */
  SSBtimes **tSSB;			/**< SSB-times DeltaT_alpha and Tdot_alpha */
  AMCoeffs **amcoe;         		/**< Amplitude Modulation coefficients */
} IFOspecifics;

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  REAL8 dFreq;			/**< frequency resolution */
  REAL8 df1dot;			/**< spindown resolution (f1 = df/dt!!) */
  LIGOTimeGPS startTime;	/**< start time of observation */
  LIGOTimeGPS refTime;		/**< reference-time for pulsar-parameters in SBB frame */
  REAL8Vector *fkdot;		/**< vector of frequency + derivatives (spindowns)*/
  DopplerRegion searchRegion;	/**< parameter-space region to search over */
  EphemerisData *edat;		/**< ephemeris data (from LALInitBarycenter()) */
  CHAR *skyRegionString;	/**< sky-region to search (polygon defined by list of points) */
  IFOspecifics ifos;		/**< IFO-specific configuration data  */
} ConfigVariables;

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

ConfigVariables GV;		/**< global container for various derived configuration settings */

/* ----- User-variables: can be set from config-file or command-line */
INT4 uvar_Dterms;
CHAR *uvar_IFO;
CHAR *uvar_IFO2;
BOOLEAN uvar_SignalOnly;
REAL8 uvar_Freq;
REAL8 uvar_FreqBand;
REAL8 uvar_dFreq;
REAL8 uvar_Alpha;
REAL8 uvar_dAlpha;
REAL8 uvar_AlphaBand;
REAL8 uvar_Delta;
REAL8 uvar_dDelta;
REAL8 uvar_DeltaBand;
REAL8 uvar_f1dot;
REAL8 uvar_df1dot;
REAL8 uvar_f1dotBand;
REAL8 uvar_Fthreshold;
CHAR *uvar_ephemDir;
CHAR *uvar_ephemYear;
INT4  uvar_gridType;
INT4  uvar_metricType;
REAL8 uvar_metricMismatch;
CHAR *uvar_skyRegion;
CHAR *uvar_DataFiles;
CHAR *uvar_DataFiles2;
BOOLEAN uvar_help;
CHAR *uvar_outputLabel;
CHAR *uvar_outputFstat;
CHAR *uvar_outputLoudest;
CHAR *uvar_skyGridFile;
CHAR *uvar_outputSkyGrid;
CHAR *uvar_workingDir;
REAL8 uvar_dopplermax;
INT4 uvar_RngMedWindow;
REAL8 uvar_refTime;
INT4 uvar_SSBprecision;

/* ---------- local prototypes ---------- */

int main(int argc,char *argv[]);
void initUserVars (LALStatus *);

void InitFStat (LALStatus *, ConfigVariables *cfg);
void InitFStatDetector (LALStatus *, ConfigVariables *cfg, UINT4 nD);

void CreateNautilusDetector (LALStatus *, LALDetector *Detector);
void Freemem(LALStatus *,  ConfigVariables *cfg);

void WriteFStatLog (LALStatus *, CHAR *argv[]);
void checkUserInputConsistency (LALStatus *);

const char *va(const char *format, ...);	/* little var-arg string helper function */

/*---------- empty initializers ---------- */
static const PulsarTimesParamStruc empty_PulsarTimesParamStruc;
static const BarycenterInput empty_BarycenterInput;

/*----------------------------------------------------------------------*/
/* CODE starts here */
/*----------------------------------------------------------------------*/

/** 
 * MAIN function of ComputeFStatistic code.
 * Calculate the F-statistic over a given portion of the parameter-space
 * and write a list of 'candidates' into a file(default: 'Fstats').
 */
int main(int argc,char *argv[]) 
{
  LALStatus status = blank_status;	/* initialize status */

  DopplerScanInit scanInit;		/* init-structure for DopperScanner */
  DopplerScanState thisScan = empty_DopplerScanState; /* current state of the Doppler-scan */
  DopplerPosition dopplerpos;		/* current search-parameters */
  SkyPosition thisPoint;
  FILE *fpOut=NULL;
  UINT4 loopcounter;
  CHAR loudestEntry[512];
  CHAR buf[512];
  REAL8 loudestF = 0;

  UINT4 nD;         /** index over number of Detectors**/
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


  /* prepare initialization of DopplerScanner to step through paramter space */
  scanInit.dAlpha = uvar_dAlpha;
  scanInit.dDelta = uvar_dDelta;
  scanInit.gridType = uvar_gridType;
  scanInit.metricType = uvar_metricType;
  scanInit.metricMismatch = uvar_metricMismatch;

  /*----- figure out total observation time */
  {
    LIGOTimeGPS t0, t1;
    UINT4 numSFTs = GV.ifos.sftVects[0]->length;
    REAL8 tObs;

    t0 = GV.ifos.sftVects[0]->data[0].epoch;
    t1 = GV.ifos.sftVects[0]->data[numSFTs-1].epoch;
    LAL_CALL (LALDeltaFloatGPS (&status, &tObs, &t1, &t0), &status);	/* t1 - t0 */
    tObs += 1.0 / (GV.ifos.sftVects[0]->data[0].deltaF );		/* +tSFT */
    GV.startTime = t0;
    scanInit.obsDuration = tObs;
  }

  scanInit.obsBegin = GV.startTime;
  scanInit.Detector = &(GV.ifos.Detectors[0]);
  scanInit.ephemeris = GV.edat;		/* used by Ephemeris-based metric */
  scanInit.skyGridFile = uvar_skyGridFile;

  scanInit.searchRegion = GV.searchRegion;
  
  if (lalDebugLevel) printf ("\nSetting up template grid ...");
  
  LAL_CALL ( InitDopplerScan ( &status, &thisScan, &scanInit), &status); 
  
  /* ---------- set Frequency- and spindown-resolution if not input by user ----------*/
  if ( !LALUserVarWasSet( &uvar_dFreq ) )
    GV.dFreq = thisScan.dFreq;
  else
    GV.dFreq = uvar_dFreq;
  
  if( !LALUserVarWasSet( &uvar_df1dot) ) 
    GV.df1dot = thisScan.df1dot;
  else
    GV.df1dot = uvar_df1dot;
  
  if ( lalDebugLevel ) {
    printf ("\nDEBUG: actual grid-spacings: dFreq = %g, df1dot = %g\n\n",
	    GV.dFreq, GV.df1dot);
  }
  /*----------------------------------------------------------------------*/
  if (lalDebugLevel) printf ("done.\n");
  if ( uvar_outputSkyGrid ) {
    printf ("\nNow writing sky-grid into file '%s' ...", uvar_outputSkyGrid);
    LAL_CALL (writeSkyGridFile( &status, thisScan.grid, uvar_outputSkyGrid, &scanInit), &status);
    printf (" done.\n\n");
  }
  
  /* if a complete output of the F-statistic file was requested,
   * we open and prepare the output-file here */
  if (uvar_outputFstat)
    {
      if ( (fpOut = fopen (uvar_outputFstat, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputFstat);
	  exit(-1);
	}
    } /* if outputFstat */
  
  if (lalDebugLevel) printf ("\nStarting main search-loop.. \n");
  
  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
   */
  loopcounter = 0;
  while (1)
    {
      UINT4 nFreq, nf1dot;	/* number of frequency- and f1dot-bins */
      UINT4 iFreq, if1dot;  	/* counters over freq- and f1dot- bins */

      LAL_CALL (NextDopplerPos( &status, &dopplerpos, &thisScan ), &status);
      if (thisScan.state == STATE_FINISHED) /* scanned all DopplerPositions yet? */
	break;
      
      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
      thisPoint.longitude = dopplerpos.Alpha;
      thisPoint.latitude = dopplerpos.Delta;
      thisPoint.system = COORDINATESYSTEM_EQUATORIAL;
      LAL_CALL (LALNormalizeSkyPosition(&status, &thisPoint, &thisPoint), &status);


      nFreq =  (UINT4)(uvar_FreqBand  / GV.dFreq + 0.5) + 1;  
      nf1dot = (UINT4)(uvar_f1dotBand / GV.df1dot+ 0.5) + 1; 
      /*----- loop over first-order spindown values */
      for (if1dot = 0; if1dot < nf1dot; if1dot ++)
	{
	  GV.fkdot->data[1] = uvar_f1dot + if1dot * GV.df1dot;
	  
	  /* Loop over frequencies to be demodulated */
	  for ( iFreq = 0 ; iFreq < nFreq ; iFreq ++ )
	    {
	      GV.fkdot->data[0] = uvar_Freq + iFreq * GV.dFreq;	
	      
	      for(nD=0; nD < GV.ifos.length; nD++)
		{
		  
		  /*----- calculate SSB-times DeltaT_alpha and Tdot_alpha for this skyposition */
		  LAL_CALL ( LALGetSSBtimes (&status, 
					     GV.ifos.tSSB[nD],
					     GV.ifos.DetectorStates[nD], 
					     thisPoint, 
					     GV.refTime,
					     uvar_SSBprecision), &status);
		  
		  /*----- calculate skypos-specific coefficients a_i, b_i, A, B, C, D */
		  LAL_CALL ( LALGetAMCoeffs (&status, GV.ifos.amcoe[nD], GV.ifos.DetectorStates[nD], 
					     thisPoint), &status);
		  
		  /** Caculate F-statistic using XLALComputeFaFb() */
		  {
		    Fcomponents FaFb;
		    REAL4 fact;
		    REAL4 At, Bt, Ct;
		    REAL4 FaRe, FaIm, FbRe, FbIm;
		    REAL8 Fstat;

		    /* prepare quantities to calculate Fstat from Fa and Fb */
		    fact = 4.0f / (1.0f * GV.ifos.sftVects[nD]->length * GV.ifos.amcoe[nD]->D);

		    /* In the signal-only case (only for testing using fake data),
		     * we did not apply any normalization to the data, and we need
		     * the use the correct factor now (taken from CFS_v1)
		     * [see Xavie's notes on LALDemod() for details]
		     */
		    if ( uvar_SignalOnly )
		      {
			REAL8 Ns = 2.0 * GV.ifos.sftVects[0]->data[0].data->length;
			REAL8 Tsft = 1.0 / (GV.ifos.sftVects[0]->data[0].deltaF );
			REAL8 norm = Tsft / ( Ns * Ns ); 

			fact *= norm;
		      }
		    else
		      {
			/* NOTE: we normalized the data by a double-sided PSD, (LALNormalizeSFTVect),
			 * therefore we apply another factor of 1/2 now with respect to the 
			 * equations in JKS, which are based on the single-sided PSD:
			 */
			fact *= 0.5f;
		      }

		    At = GV.ifos.amcoe[nD]->A;
		    Bt = GV.ifos.amcoe[nD]->B;
		    Ct = GV.ifos.amcoe[nD]->C;

		    if ( XLALComputeFaFb (&FaFb, GV.ifos.sftVects[nD], GV.fkdot, GV.ifos.tSSB[nD], 
					  GV.ifos.amcoe[nD], uvar_Dterms) != 0)
		      {
			LALPrintError ("\nXALNewLALDemod() failed\n");
			XLAL_ERROR ("XLALcomputeFStat", XLAL_EFUNC);
		      }
		    
		    FaRe = FaFb.Fa.re;
		    FaIm = FaFb.Fa.im;
		    
		    FbRe = FaFb.Fb.re;
		    FbIm = FaFb.Fb.im;
		    
		    /* calculate F-statistic from  Fa and Fb */
		    Fstat = fact * (Bt * (FaRe*FaRe + FaIm*FaIm) 
				    + At * (FbRe*FbRe + FbIm*FbIm) 
				    - 2.0f * Ct *(FaRe*FbRe + FaIm*FbIm) );

		    /* now, if user requested it, we output ALL F-statistic results above threshold */
		    if ( uvar_outputFstat || uvar_outputLoudest )
		      {
			REAL8 freq = GV.fkdot->data[0];

			if ( Fstat > uvar_Fthreshold )
			  {
			    LALSnprintf (buf, 511, "%8.7f %8.7f %16.12f %.17g %10.6g\n", 
				      dopplerpos.Alpha, dopplerpos.Delta, 
				      freq, GV.fkdot->data[1], 2.0 * Fstat);
			    buf[511] = 0;
			    if ( fpOut )
			      fprintf (fpOut, buf );
			  } /* if F > threshold */
		      }  

		    if ( Fstat > loudestF )
		      {
			loudestF = Fstat;
			strcpy ( loudestEntry,  buf );
		      }
			

	      
		  } /* Calculate F-statistic */
		  
		} /* End of loop over detectors */
	      
	    } /* for i < nBins: loop over frequency-bins */
	  
	} /* For GV.spinImax: loop over spindowns */
      
      loopcounter ++;
      if (lalDebugLevel) 
	printf ("\
Search progress: %5.1f%%", (100.0* loopcounter / thisScan.numGridPoints));
      
    } /*  while SkyPos : loop over skypositions */
  
  if (uvar_outputFstat && fpOut)
    fclose (fpOut);

  /* now write loudest canidate into separate file ".loudest" */
  if ( uvar_outputLoudest )
    {
      FILE *fpLoudest;
      if ( (fpLoudest = fopen (uvar_outputLoudest, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputLoudest);
	  return COMPUTEFSTATISTICC_ESYS;
	}
      fprintf (fpLoudest, "%s", loudestEntry );
      fclose(fpLoudest);
    } /* write loudest candidate to file */


  
  if (lalDebugLevel) printf ("\nSearch finished.\n");
  
  /* Free memory */
  LAL_CALL ( FreeDopplerScan(&status, &thisScan), &status);
  
  if ( GV.searchRegion.skyRegionString )
    LALFree ( GV.searchRegion.skyRegionString );

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
  uvar_Dterms 	= 16;
  uvar_FreqBand = 0.0;
  uvar_dFreq 	= 0.0;
  uvar_Alpha 	= 0.0;
  uvar_Delta 	= 0.0;
  uvar_AlphaBand = 0;
  uvar_DeltaBand = 0;
  uvar_dAlpha 	= 0.001;
  uvar_dDelta 	= 0.001;
  uvar_skyRegion = NULL;

  uvar_ephemYear = LALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar_ephemYear, EPHEM_YEARS);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar_ephemDir = LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar_ephemDir, DEFAULT_EPHEMDIR);

  uvar_SignalOnly = FALSE;
 
  uvar_f1dot = 0.0;
  uvar_df1dot 	= 0.0;
  uvar_f1dotBand = 0.0;
  
  uvar_Fthreshold = 10.0;
  uvar_metricType =  LAL_PMETRIC_NONE;
  uvar_gridType = GRID_FLAT;

  uvar_metricMismatch = 0.02;

  uvar_help = FALSE;
  uvar_outputLabel = NULL;

  uvar_outputFstat = NULL;

  uvar_skyGridFile = NULL;

  uvar_workingDir = LALMalloc(512);
  strcpy(uvar_workingDir, ".");

  uvar_dopplermax =  1.05e-4;
  uvar_RngMedWindow = 50;	/* for running-median */

  uvar_SSBprecision = SSBPREC_RELATIVISTIC;

  /* register all our user-variables */
  LALregBOOLUserVar(status, 	help, 		'h', UVAR_HELP,     "Print this message"); 
  LALregREALUserVar(status, 	Freq, 		'f', UVAR_REQUIRED, "Starting search frequency in Hz");
  LALregREALUserVar(status, 	FreqBand, 	'b', UVAR_OPTIONAL, "Search frequency band in Hz");
  LALregREALUserVar(status,     dFreq,          'r', UVAR_OPTIONAL, "Frequency resolution in Hz (default: 1/(2*Tsft*Nsft)");
  LALregREALUserVar(status, 	f1dot, 		's', UVAR_OPTIONAL, "First spindown parameter f1dot");
  LALregREALUserVar(status, 	f1dotBand, 	'm', UVAR_OPTIONAL, "Search-band for f1dot");
  LALregREALUserVar(status, 	df1dot, 	'e', UVAR_OPTIONAL, "Resolution for f1dot (default 1/(2*Tobs*tSFT*Nsft)");
  LALregREALUserVar(status, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	AlphaBand, 	'z', UVAR_OPTIONAL, "Band in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	dAlpha, 	'l', UVAR_OPTIONAL, "Resolution in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	DeltaBand, 	'c', UVAR_OPTIONAL, "Band in delta (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	dDelta, 	'g', UVAR_OPTIONAL, "Resolution in delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(status,	skyRegion, 	'R', UVAR_OPTIONAL, "ALTERNATIVE: Specify sky-region by polygon (or use 'allsky')");
  LALregSTRINGUserVar(status,	DataFiles, 	'D', UVAR_REQUIRED, "File-pattern specifying (first) set of data SFT-files"); 
  LALregSTRINGUserVar(status, 	IFO, 		'I', UVAR_REQUIRED, "Detector: GEO(0), LLO(1), LHO(2), NAUTILUS(3), VIRGO(4), TAMA(5), CIT(6)");
  LALregSTRINGUserVar(status,	DataFiles2, 	 0,  UVAR_OPTIONAL, "File-pattern specifying second set of data SFT-files");
  LALregSTRINGUserVar(status, 	IFO2, 		 0,  UVAR_OPTIONAL, "Detector corresponding to second data-set (--DataFiles2)"); 
  LALregSTRINGUserVar(status,	ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(status,	ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregBOOLUserVar(status, 	SignalOnly, 	'S', UVAR_OPTIONAL, "Signal only flag");
  LALregREALUserVar(status, 	Fthreshold,	'F', UVAR_OPTIONAL, "Signal Set the threshold for selection of 2F");
  LALregINTUserVar(status, 	gridType,	 0 , UVAR_OPTIONAL, "Template grid: 0=flat, 1=isotropic, 2=metric, 3=file");
  LALregINTUserVar(status, 	metricType,	'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserVar(status, 	metricMismatch,	'X', UVAR_OPTIONAL, "Maximal allowed mismatch for metric tiling");
  LALregSTRINGUserVar(status,	outputLabel,	'o', UVAR_OPTIONAL, "Label to be appended to all output file-names");
  LALregSTRINGUserVar(status,	skyGridFile,	 0,  UVAR_OPTIONAL, "Load sky-grid from this file.");
  LALregREALUserVar(status,	refTime,	 0,  UVAR_OPTIONAL, "SSB reference time for pulsar-paramters");
  LALregREALUserVar(status, 	dopplermax, 	'q', UVAR_OPTIONAL, "Maximum doppler shift expected");  
  LALregSTRINGUserVar(status,	outputFstat,	 0,  UVAR_OPTIONAL, "Output-file for F-statistic field over the parameter-space");

  /* more experimental and unofficial stuff follows here */
  LALregINTUserVar (status, 	SSBprecision,	 0,  UVAR_DEVELOPER, "Precision to use for time-transformation to SSB: 0=Newtonian 1=relativistic");
  LALregINTUserVar(status, 	RngMedWindow,	'k', UVAR_DEVELOPER, "Running-Median window size");
  LALregINTUserVar(status,	Dterms,		't', UVAR_DEVELOPER, "Number of terms to keep in Dirichlet kernel sum");
  LALregSTRINGUserVar(status,   workingDir,     'w', UVAR_DEVELOPER, "Directory to be made the working directory, . is default");
  LALregSTRINGUserVar(status,	outputSkyGrid,	 0,  UVAR_DEVELOPER, "Write sky-grid into this file.");
  LALregSTRINGUserVar(status,   outputLoudest,	 0,  UVAR_DEVELOPER, 
		      "Output-file for the loudest F-statistic candidate in this search");


  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* initUserVars() */


/** Do some basic initializations of the F-statistic code before starting the main-loop.
 * Things we do in this function: 
 * \li check consistency of user-input
 * \li prepare ephemeris-data and determine SFT input-files to be loaded
 * \li set some defaults + allocate memory 
 * \li Return 'derived' configuration settings in the struct \em ConfigVariables
 * 
 */
void
InitFStat (LALStatus *status, ConfigVariables *cfg)
{
  UINT4 nDet;
  UINT4 nD;
  INITSTATUS (status, "InitFStat", rcsid);
  ATTATCHSTATUSPTR (status);

  /*---------- Initialize Ephemeris-data ---------- */
  {
#define FNAME_LENGTH 1024
    CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
    CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */
    LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
    INT4 leap;

    if (LALUserVarWasSet (&uvar_ephemDir) )
      {
	LALSnprintf(EphemEarth, FNAME_LENGTH, "%s/earth%s.dat", uvar_ephemDir, uvar_ephemYear);
	LALSnprintf(EphemSun, FNAME_LENGTH, "%s/sun%s.dat", uvar_ephemDir, uvar_ephemYear);
      }
    else
      {
	LALSnprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", uvar_ephemYear);
	LALSnprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  uvar_ephemYear);
      }
    EphemEarth[FNAME_LENGTH-1]=0;
    EphemSun[FNAME_LENGTH-1]=0;
    
    cfg->edat = LALCalloc(1, sizeof(EphemerisData));
    /* NOTE: the 'ephiles' are ONLY ever used in LALInitBarycenter, which is
     * why we could use local variables (EphemEarth, EphemSun) to initialize them.
     */
    cfg->edat->ephiles.earthEphemeris = EphemEarth;
    cfg->edat->ephiles.sunEphemeris = EphemSun;
    
    TRY (LALLeapSecs (status->statusPtr, &leap, &(cfg->startTime), &formatAndAcc), status);
    cfg->edat->leap = (INT2) leap;

    TRY (LALInitBarycenter(status->statusPtr, cfg->edat), status);               

  } /* end: init ephemeris data */


  /* ---------- initialize+check  template-grid related parameters ---------- */
  {
    BOOLEAN haveSkyRegion, haveAlphaDelta;

    haveSkyRegion  = (uvar_skyRegion != NULL);
    haveAlphaDelta = (LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta) );

    /* pre-process template-related input */
    if (haveSkyRegion)
      {
	cfg->skyRegionString = LALCalloc(1, strlen(uvar_skyRegion)+1);
	strcpy (cfg->skyRegionString, uvar_skyRegion);
      }
    else if (haveAlphaDelta)	/* parse this into a sky-region */
      {
	TRY ( SkySquare2String( status->statusPtr, &(cfg->skyRegionString), 
				uvar_Alpha, uvar_Delta, 
				uvar_AlphaBand, uvar_DeltaBand), status);
      }
    
  } /* end: template-grid stuff */

  /* ----- count number of detectors */
  if ( LALUserVarWasSet(&uvar_DataFiles2) )
    nDet = 2;
  else
    nDet = 1;

  cfg->ifos.length = nDet;

  /* ----- set up detector-specific data */
  cfg->ifos.Detectors = LALCalloc ( nDet,  sizeof( *(cfg->ifos.Detectors) ) );
  cfg->ifos.sftVects =  LALCalloc ( nDet,  sizeof( *(cfg->ifos.sftVects) ) );
  /* cfg->ifos.DetStates =  LALCalloc ( nDet,  sizeof( *(cfg->ifos.DetStates) ) );*/
  cfg->ifos.DetectorStates =  LALCalloc ( nDet,  sizeof( *(cfg->ifos.DetectorStates) ) );
  cfg->ifos.tSSB = LALCalloc( nDet, sizeof( *(cfg->ifos.tSSB) ) );
  cfg->ifos.amcoe = LALCalloc ( nDet,  sizeof( *(cfg->ifos.amcoe) ) );
  TRY ( LALDCreateVector (status->statusPtr, &(cfg->fkdot), 2), status);




  /**---------------------------------------------------------**/
  /** Starting the Loop for different Detectors **/
  /** At this moment we are trying to match it for single detector **/

  for(nD=0; nD < cfg->ifos.length; nD++)
    {
      REAL8 tSFT;
      LIGOTimeGPSVector *timestamps = NULL;
      UINT4 i;

      /* main initialization of the code: */
      TRY ( InitFStatDetector(status->statusPtr, cfg, nD), status);
      
      /* figure out duration of an SFT */
      tSFT = 1.0 / cfg->ifos.sftVects[nD]->data[0].deltaF ;

      /* extract vector of timestamps from SFT-vector (for call to GetDetectorStates()) */
      TRY (LALCreateTimestampVector (status->statusPtr, &timestamps, cfg->ifos.sftVects[nD]->length ), 
	   status);
      for (i=0; i < timestamps->length; i++)
	timestamps->data[i] = cfg->ifos.sftVects[nD]->data[i].epoch;	
      
      /* obtain detector positions and velocities, together with LMSTs for the SFT midpoints 
       * (i.e. shifted by tSFT/2) */
      TRY (LALGetDetectorStates(status->statusPtr, &(GV.ifos.DetectorStates[nD]), timestamps, 
				&(GV.ifos.Detectors[nD]), cfg->edat, tSFT / 2.0 ), status);

      /* free timstamps-vector */
      TRY ( LALDestroyTimestampVector (status->statusPtr, &timestamps), status );
      timestamps = NULL;

    } /* end of loop over different detectors */
  
  /* determine start-time from first set of SFTs */
  cfg->startTime = cfg->ifos.sftVects[0]->data[0].epoch;
  

  /*---------- Set reference-time: ----------*/
  if ( LALUserVarWasSet(&uvar_refTime)) {
    TRY ( LALFloatToGPS (status->statusPtr, &(cfg->refTime), &uvar_refTime), status);
  } else
    cfg->refTime = cfg->startTime;


  DETATCHSTATUSPTR (status);
  RETURN (status);


} /* InitFStat() */


void
InitFStatDetector (LALStatus *status, ConfigVariables *cfg, UINT4 nD)
{
  CHAR *IFO;

  INITSTATUS (status, "InitFStatDetector", rcsid);
  ATTATCHSTATUSPTR (status);

  /*----------------------------------------------------------------------
   * load SFT data-files
   */
  /* min and max physical frequency to read */
  {
    REAL8 f_min, f_max, f1dot_min, f1dot_max;
    REAL8 deltaFreq_min, deltaFreq_max;	/* frequency-shift due to spindown */

    UINT4 wings;
    CHAR *fname = NULL;
    SFTVector *headers = NULL;
    LIGOTimeGPS startTime, endTime;
    REAL8 refTime;
    REAL8 duration;

    REAL8 refShift;		/* time between reference-time and first-SFT timestamp */

    if ( nD == 0 ) fname = uvar_DataFiles;
    if ( nD == 1 ) fname = uvar_DataFiles2;
    if ( nD > 1 ) {
      LALPrintError ("\nERROR: currently maximally two detectors are supported!\n\n");
      ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
    }

    /* frequency/spindown ranges (at reference-time tRef) */
    f_min = MYMIN(uvar_Freq, uvar_Freq + uvar_FreqBand);
    f_max = MYMAX(uvar_Freq, uvar_Freq + uvar_FreqBand);

    f1dot_min = MYMIN(uvar_f1dot, uvar_f1dot + uvar_f1dotBand );
    f1dot_max = MYMAX(uvar_f1dot, uvar_f1dot + uvar_f1dotBand );

    /* read in all SFT-header in order to get data time-span */
    TRY ( LALGetSFTheaders (status->statusPtr, &headers, fname, NULL, NULL ), status );

    /* figure out observation time-span */
    startTime = headers->data[0].epoch;
    endTime = headers->data[ headers->length - 1 ].epoch;
    duration = GPS2REAL8(endTime) - GPS2REAL8(startTime);


    TRY ( LALDestroySFTVector (status->statusPtr, &headers), status );

    if ( LALUserVarWasSet ( &uvar_refTime ) )
      refTime = uvar_refTime;
    else
      refTime = GPS2REAL8 (startTime);

    /* time from refTime to startTime */
    refShift = GPS2REAL8(startTime) - refTime;

    /* first propagate frequency-interval from tRef to startTime: */
    deltaFreq_max = f1dot_max * refShift;
    deltaFreq_min = f1dot_min * refShift;


    /* now propagate from startTime to endTime, keeping the 'covering' Band */
    if ( f1dot_max > 0 )
      deltaFreq_max += f1dot_max * duration;
    if ( f1dot_min < 0 )
      deltaFreq_min += f1dot_min * duration;

    /* final 'propagated' physical frequency-interval */
    f_min += deltaFreq_min;
    f_max += deltaFreq_max;

    /* ----- correct for maximal dopper-shift due to earth's motion */
    f_min *= (1.0 - uvar_dopplermax);
    f_max *= (1.0 + uvar_dopplermax);
    
    /* ----- load the SFT-vector ----- */
    wings = MYMAX(uvar_Dterms, uvar_RngMedWindow/2 +1);

    TRY ( LALReadSFTfiles(status->statusPtr, &(cfg->ifos.sftVects[nD]), f_min, f_max, wings, fname), 
	  status);

    /* this is normalized by 1/sqrt(Sh), where Sh is the median of |X|^2  
     * NOTE: this corresponds to a double-sided PSD, therefore we need to 
     * divide by another factor of 2 with respect to the JKS formulae.
     */
    if ( ! uvar_SignalOnly ) {
      TRY ( LALNormalizeSFTVect (status->statusPtr, cfg->ifos.sftVects[nD], uvar_RngMedWindow, 0 ), 
	    status );
    }
    

    /* ----- find search region ----- */
    {
      BOOLEAN haveAlphaDelta = LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta);

      cfg->searchRegion.Freq = f_min;
      cfg->searchRegion.FreqBand = f_max - f_min;

      cfg->searchRegion.f1dot = f1dot_min;
      cfg->searchRegion.f1dotBand = f1dot_max - f1dot_min;

      /* get sky-region */
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
    } /* find search-region */


  } /* SFT-loading */

  if(nD == 0)
    IFO=uvar_IFO;
  else 
    IFO=uvar_IFO2;







  /*----------------------------------------------------------------------
   * initialize detector 
   */
  if ( !strcmp (IFO, "GEO") || !strcmp (IFO, "0") ) 
    cfg->ifos.Detectors[nD] = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  else if ( !strcmp (IFO, "LLO") || ! strcmp (IFO, "1") ) 
    cfg->ifos.Detectors[nD] = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if ( !strcmp (IFO, "LHO") || !strcmp (IFO, "2") )
    cfg->ifos.Detectors[nD] = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else if ( !strcmp (IFO, "NAUTILUS") || !strcmp (IFO, "3"))
    {
      TRY (CreateNautilusDetector (status->statusPtr, &(cfg->ifos.Detectors[0])), status);
    }
  else if ( !strcmp (IFO, "VIRGO") || !strcmp (IFO, "4") )
    cfg->ifos.Detectors[nD] = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  else if ( !strcmp (IFO, "TAMA") || !strcmp (IFO, "5") )
    cfg->ifos.Detectors[nD] = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  else if ( !strcmp (IFO, "CIT") || !strcmp (IFO, "6") )
    cfg->ifos.Detectors[nD] = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
  else
    {
      LALPrintError ("\nUnknown detector. Currently allowed are \
'GEO', 'LLO', 'LHO', 'NAUTILUS', 'VIRGO', 'TAMA', 'CIT' or '0'-'6'\n\n");
      ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
    }


  /* ----------------------------------------------------------------------
   * initialize + allocate space for AM-coefficients and SSB-times 
   */
  {
    AMCoeffs *amc = NULL;
    SSBtimes *tSSB = NULL;

    /* Allocate space for AMCoeffs */
    if ( (amc = LALCalloc(1, sizeof(AMCoeffs))) == NULL) {
      ABORT (status, COMPUTEFSTATISTICC_EMEM, COMPUTEFSTATISTICC_MSGEMEM);
    }
    TRY (LALSCreateVector(status->statusPtr, &(amc->a), (UINT4) cfg->ifos.sftVects[nD]->length), status);
    TRY (LALSCreateVector(status->statusPtr, &(amc->b), (UINT4) cfg->ifos.sftVects[nD]->length), status);
    
    cfg->ifos.amcoe[nD] = amc;
    
    /* allocate memory of the SSB-times: DeltaT_alpha and Tdot_alpha */
    if ( (tSSB = LALCalloc(1, sizeof(SSBtimes))) == NULL) {
      ABORT (status, COMPUTEFSTATISTICC_EMEM, COMPUTEFSTATISTICC_MSGEMEM);
    }
    TRY ( LALDCreateVector(status->statusPtr, &(tSSB->DeltaT), cfg->ifos.sftVects[nD]->length), status );
    TRY ( LALDCreateVector(status->statusPtr, &(tSSB->Tdot),   cfg->ifos.sftVects[nD]->length), status );

    cfg->ifos.tSSB[nD] = tSSB;

  } /* end: init AM- and SSBtimes */

  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* InitFStatDetector() */

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
      ABORT (status, COMPUTEFSTATISTICC_EMEM, COMPUTEFSTATISTICC_MSGEMEM);
    }
    strcpy (fname, head);
    if (uvar_outputLabel)
      strcat (fname, uvar_outputLabel);
    strcat (fname, ".log");

    if ( (fplog = fopen(fname, "w" )) == NULL) {
      LALPrintError ("\nFailed to open log-file '%f' for writing.\n\n", fname);
      LALFree (fname);
      ABORT (status, COMPUTEFSTATISTICC_ESYS, COMPUTEFSTATISTICC_MSGESYS);
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


/*******************************************************************************/
/** Set up the \em LALDetector struct representing the NAUTILUS detector */
void
CreateNautilusDetector (LALStatus *status, LALDetector *Detector)
{
  /*   LALDetector Detector;  */
  LALFrDetector detector_params;
  LALDetectorType bar;
  LALDetector Detector1;

  INITSTATUS (status, "CreateNautilusDetector", rcsid);
  ATTATCHSTATUSPTR (status);

/*   detector_params=(LALFrDetector )LALMalloc(sizeof(LALFrDetector)); */
 
  bar=LALDETECTORTYPE_CYLBAR;
  strcpy(detector_params.name, "NAUTILUS");
  detector_params.vertexLongitudeRadians=12.67*LAL_PI/180.0;
  detector_params.vertexLatitudeRadians=41.82*LAL_PI/180.0;
  detector_params.vertexElevation=300.0;
  detector_params.xArmAltitudeRadians=0.0;
  detector_params.xArmAzimuthRadians=44.0*LAL_PI/180.0;

  TRY (LALCreateDetector(status->statusPtr, &Detector1, &detector_params, bar), status);
  
  *Detector=Detector1;

  DETATCHSTATUSPTR (status);
  RETURN (status);
  
} /* CreateNautilusDetector() */

/*******************************************************************************/
/** Free all globally allocated memory. */
void
Freemem(LALStatus *status,  ConfigVariables *cfg) 
{
  UINT4 i;
  INITSTATUS (status, "Freemem", rcsid);
  ATTATCHSTATUSPTR (status);

  for(i=0; i< cfg->ifos.length; i++)
    {
      
      /* Free SFT data */
      TRY (LALDestroySFTVector (status->statusPtr, &(cfg->ifos.sftVects[i]) ), status);	 /* the new way*/
                  
      /* Free AM-coefficients */
      TRY (LALSDestroyVector(status->statusPtr, &(cfg->ifos.amcoe[i]->a)), status);
      TRY (LALSDestroyVector(status->statusPtr, &(cfg->ifos.amcoe[i]->b)), status);
      LALFree ( cfg->ifos.amcoe[i]);
      /* Free SSB-times */
      TRY (LALDDestroyVector(status->statusPtr, &(cfg->ifos.tSSB[i]->DeltaT)), status);
      TRY (LALDDestroyVector(status->statusPtr, &(cfg->ifos.tSSB[i]->Tdot)), status);
      LALFree ( cfg->ifos.tSSB[i] );

      /* destroy DetectorStateSeries */
      TRY ( LALDestroyDetectorStateSeries (status->statusPtr, &(cfg->ifos.DetectorStates[i]) ), status);
    }
  
  LALFree ( cfg->ifos.Detectors );
  LALFree ( cfg->ifos.sftVects );
  LALFree ( cfg->ifos.DetectorStates );
  LALFree ( cfg->ifos.tSSB );
  LALFree ( cfg->ifos.amcoe );
  
  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (status->statusPtr), status);
  
  LALFree ( cfg->skyRegionString );
  
  /* Free ephemeris data */
  LALFree(cfg->edat->ephemE);
  LALFree(cfg->edat->ephemS);
  LALFree(cfg->edat);
    
  TRY (LALDDestroyVector (status->statusPtr, &(cfg->fkdot)), status);
    
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
      ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
    }      
  /* don't allow negative frequency-band for safety */
  if ( uvar_FreqBand < 0)
    {
      LALPrintError ("\nNegative value of frequency-band not allowed !\n\n");
      ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
    }

  /* don't allow negative bands (for safty in griding-routines) */
  if ( (uvar_AlphaBand < 0) ||  (uvar_DeltaBand < 0) )
    {
      LALPrintError ("\nNegative value of sky-bands not allowed (alpha or delta)!\n\n");
      ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
    }
  /* check for negative stepsizes in Freq, Alpha, Delta */
  if ( LALUserVarWasSet(&uvar_dAlpha) && (uvar_dAlpha < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dAlpha not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar_dDelta) && (uvar_dDelta < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dDelta not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar_dFreq) && (uvar_dFreq < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dFreq not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
    }

  /* grid-related checks */
  {
    BOOLEAN haveAlphaBand = LALUserVarWasSet( &uvar_AlphaBand );
    BOOLEAN haveDeltaBand = LALUserVarWasSet( &uvar_DeltaBand );
    BOOLEAN haveSkyRegion, haveAlphaDelta, haveGridFile, useGridFile, haveMetric, useMetric;

    haveSkyRegion  = (uvar_skyRegion != NULL);
    haveAlphaDelta = (LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta) );
    haveGridFile   = (uvar_skyGridFile != NULL);
    useGridFile   = (uvar_gridType == GRID_FILE);
    haveMetric     = (uvar_metricType > LAL_PMETRIC_NONE);
    useMetric     = (uvar_gridType == GRID_METRIC);


    if ( (haveAlphaBand && !haveDeltaBand) || (haveDeltaBand && !haveAlphaBand) )
      {
	LALPrintError ("\nERROR: Need either BOTH (AlphaBand, DeltaBand) or NONE.\n\n"); 
        ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
      }

    if ( !useGridFile && !(haveSkyRegion || haveAlphaDelta) )
      {
        LALPrintError ("\nNeed sky-region: either use (Alpha,Delta) or skyRegion!\n\n");
        ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
      }
    if ( haveSkyRegion && haveAlphaDelta )
      {
        LALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta)"
		       " OR skyRegion!\n\n");
        ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);
      }
    if ( useGridFile && !haveGridFile )
      {
        LALPrintError ("\nERROR: gridType=FILE, but no skyGridFile specified!\n\n");
        ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);  
      }
    if ( !useGridFile && haveGridFile )
      {
        LALWarning (status, "\nWARNING: skyGridFile was specified but not needed ..."
		    " will be ignored\n");
      }
    if ( useGridFile && (haveSkyRegion || haveAlphaDelta) )
      {
        LALWarning (status, "\nWARNING: We are using skyGridFile, but sky-region was"
		    " also specified ... will be ignored!\n");
      }
    if ( !useMetric && haveMetric) 
      {
        LALWarning (status, "\nWARNING: Metric was specified for non-metric grid..."
		    " will be ignored!\n");
      }
    if ( useMetric && !haveMetric) 
      {
        LALPrintError ("\nERROR: metric grid-type selected, but no metricType selected\n\n");
        ABORT (status, COMPUTEFSTATISTICC_EINPUT, COMPUTEFSTATISTICC_MSGEINPUT);      
      }

  } /* grid-related checks */

  RETURN (status);
} /* checkUserInputConsistency() */


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

