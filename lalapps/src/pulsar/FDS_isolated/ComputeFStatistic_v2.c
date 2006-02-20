/*
 * Copyright (C) 2005, 2006 Reinhard Prix, Iraj Gholami
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
/** \author R. Prix, Y. Ioth, Papa, X. Siemens, I. Gholami
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
#define NUM_SPINS 2		/* number of spin-values to consider: {f, fdot, f2dot, ... } */ 

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

/** struct holding all quantities that are 'detector-specific' */
typedef struct {
  SSBtimes *tSSB;			/**< SSB-times DeltaT_alpha and Tdot_alpha */
  AMCoeffs *amcoe;         		/**< Amplitude Modulation coefficients */
  REAL8Vector *weightsNoise;             /**< vector of weights */
} IFOspecifics;

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  LIGOTimeGPS startTime;		/**< start time of observation */
  REAL8 duration;			/**< total time-span of the data (all streams) in seconds */
  LIGOTimeGPS refTime;			/**< reference-time for pulsar-parameters in SBB frame */
  LALPulsarSpinRange *spinRangeRef; 	/**< pulsar spin-range at reference-time 'refTime' */
  LALPulsarSpinRange *spinRangeStart; 	/**< pulsar spin-range at start of observation 'startTime; */
  DopplerRegion searchRegion;		/**< parameter-space region to search over (FIXME) */
  EphemerisData *edat;			/**< ephemeris data (from LALInitBarycenter()) */
  UINT4 numDetectors;			/**< number of detectors */
  IFOspecifics *ifos;			/**< IFO-specific configuration data  */
  MultiSFTVector *multiSFTs;		/**< multi-IFO SFT-vectors */
  MultiDetectorStateSeries *multiDetStates; /**< pos, vel and LMSTs for detector at times t_i */
} ConfigVariables;

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

ConfigVariables GV;		/**< global container for various derived configuration settings */

/* ----- User-variables: can be set from config-file or command-line */
INT4 uvar_Dterms;
CHAR *uvar_IFO;
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
BOOLEAN uvar_help;
CHAR *uvar_outputLabel;
CHAR *uvar_outputFstat;
CHAR *uvar_outputBstat;
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
void InitFStat ( LALStatus *, ConfigVariables *cfg );

void Freemem(LALStatus *,  ConfigVariables *cfg);

void WriteFStatLog (LALStatus *, CHAR *argv[]);
void checkUserInputConsistency (LALStatus *);
int outputBeamTS( const CHAR *fname, const AMCoeffs *amcoe, const DetectorStateSeries *detStates );
void InitEphemeris (LALStatus *, EphemerisData *edat, const CHAR *ephemDir, const CHAR *ephemYear, LIGOTimeGPS epoch);

const char *va(const char *format, ...);	/* little var-arg string helper function */

/*---------- empty initializers ---------- */
static const PulsarTimesParamStruc empty_PulsarTimesParamStruc;
static const BarycenterInput empty_BarycenterInput;
static const SFTConstraints empty_SFTConstraints;

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

  DopplerScanInit scanInit;		/* init-structure for DopperScanner */
  DopplerScanState thisScan = empty_DopplerScanState; /* current state of the Doppler-scan */
  DopplerPosition dopplerpos;		/* current search-parameters */
  SkyPosition thisPoint;
  FILE *fpFstat = NULL;
  FILE *fpBstat = NULL;
  UINT4 loopcounter;
  CHAR loudestEntry[512];
  CHAR buf[512];
  REAL8 loudestF = 0;
  REAL8Vector *fkdot = NULL;
  
  
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
  scanInit.projectMetric = TRUE;
  scanInit.obsDuration = GV.duration;
  scanInit.obsBegin = GV.startTime;
  scanInit.Detector = &(GV.multiDetStates->data[0]->detector);	/* FIXME: need multi-IFO metric */
  scanInit.ephemeris = GV.edat;		/* used by Ephemeris-based metric */
  scanInit.skyGridFile = uvar_skyGridFile;

  scanInit.searchRegion = GV.searchRegion;
  
  if (lalDebugLevel) printf ("\nSetting up template grid ...");
  
  LAL_CALL ( InitDopplerScan ( &status, &thisScan, &scanInit), &status); 
  
  /* ---------- set Frequency- and spindown-resolution if not input by user ----------*/
  if ( LALUserVarWasSet( &uvar_dFreq ) )
    thisScan.dFreq = uvar_dFreq;
  
  if( LALUserVarWasSet( &uvar_df1dot) ) 
    thisScan.df1dot = uvar_df1dot;
  
  if ( lalDebugLevel ) {
    printf ("\nDEBUG: actual grid-spacings: dFreq = %g, df1dot = %g\n\n", 
	    thisScan.dFreq, thisScan.df1dot );
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
      if ( (fpFstat = fopen (uvar_outputFstat, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputFstat);
	  return (COMPUTEFSTATISTICC_ESYS);
	}
    } /* if outputFstat */
  
  if (uvar_outputBstat)
    {
      if ( (fpBstat = fopen (uvar_outputBstat, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputBstat);
	  return (COMPUTEFSTATISTICC_ESYS);
	}
    } /* if outputFstat */

  if (uvar_outputBstat)
    {
      if ( (fpBstat = fopen (uvar_outputBstat, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputBstat);
	  return (COMPUTEFSTATISTICC_ESYS);
	}
    } /* if outputBstat */


  if (lalDebugLevel) printf ("\nStarting main search-loop.. \n");
  
  if ( ( fkdot = XLALCreateREAL8Vector ( NUM_SPINS ) ) == NULL ) {
    return COMPUTEFSTATISTICC_EMEM;
  }
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


      nFreq =  (UINT4)(GV.spinRangeStart->fkdotBand->data[0] / thisScan.dFreq  + 0.5) + 1;  
      nf1dot = (UINT4)(GV.spinRangeStart->fkdotBand->data[1] / thisScan.df1dot + 0.5) + 1; 

      /*----- loop over first-order spindown values */
      for (if1dot = 0; if1dot < nf1dot; if1dot ++)
	{
	  fkdot->data[1] = GV.spinRangeStart->fkdot->data[1] + if1dot * thisScan.df1dot;
	  
	  /* Loop over frequencies to be demodulated */
	  for ( iFreq = 0 ; iFreq < nFreq ; iFreq ++ )
	    {
	      Fcomponents FaFb;
	      REAL4 fact = 0;
	      REAL4 At = 0.0, Bt = 0.0, Ct = 0.0, Dt = 0.0;
	      REAL4 FaRe = 0.0, FaIm = 0.0, FbRe = 0.0, FbIm = 0.0;
	      REAL8 Fstat;
	      REAL8 Bstat;
	      UINT4 M, alpha;
	      	 
	      fkdot->data[0] = GV.spinRangeStart->fkdot->data[0] + iFreq * thisScan.dFreq;
	      
	      for(nD=0; nD < GV.numDetectors; nD++)
		{
		  /*----- calculate SSB-times DeltaT_alpha and Tdot_alpha for this skyposition */
		  LAL_CALL ( LALGetSSBtimes (&status, 
					     GV.ifos[nD].tSSB,
					     GV.multiDetStates->data[nD], 
					     thisPoint, 
					     GV.startTime,
					     uvar_SSBprecision), &status);
		  
		  /*----- calculate skypos-specific coefficients a_i, b_i, A, B, C, D */
		  LAL_CALL ( LALGetAMCoeffs (&status, 
					     GV.ifos[nD].amcoe, 
					     GV.multiDetStates->data[nD],
					     thisPoint), &status);

		  /* initilizing the A, B and C to be set at zero */
		  GV.ifos[nD].amcoe->A = 0.0;
		  GV.ifos[nD].amcoe->B = 0.0;
		  GV.ifos[nD].amcoe->C = 0.0;
		  
		  /* Calculate the Coefficients a(t) and b(t) */
		  
		  M = 1.0f * GV.multiSFTs->data[nD]->length;
		  
		  for(alpha = 0; alpha < M; alpha++)
		    {
		      REAL4 ahat;
		      REAL4 bhat;
		      
		      ahat = (GV.ifos[nD].amcoe->a->data[alpha]) * sqrt(GV.ifos[nD].weightsNoise->data[alpha]);
		      bhat = (GV.ifos[nD].amcoe->b->data[alpha]) * sqrt(GV.ifos[nD].weightsNoise->data[alpha]);
		      
		      GV.ifos[nD].amcoe->a->data[alpha] = ahat;
		      GV.ifos[nD].amcoe->b->data[alpha] = bhat;
		      
		      /* sum A, B, C on the fly */
		      GV.ifos[nD].amcoe->A += ahat * ahat;
		      GV.ifos[nD].amcoe->B += bhat * bhat;
		      GV.ifos[nD].amcoe->C += ahat * bhat;
		    }
		  
		  At += GV.ifos[nD].amcoe->A;
		  Bt += GV.ifos[nD].amcoe->B;
		  Ct += GV.ifos[nD].amcoe->C;
		  
		  /** Caculate F-statistic using XLALComputeFaFb() */
		  /* prepare quantities to calculate Fstat from Fa and Fb */
		  
		  if ( XLALComputeFaFb (&FaFb, GV.multiSFTs->data[nD], fkdot, GV.ifos[nD].tSSB, 
					GV.ifos[nD].amcoe, uvar_Dterms) != 0)
		    {
		      LALPrintError ("\nXALNewLALDemod() failed\n");
		      XLAL_ERROR ("XLALcomputeFStat", XLAL_EFUNC);
		    }
		  
		  FaRe += FaFb.Fa.re;
		  FaIm += FaFb.Fa.im;
		  
		  FbRe += FaFb.Fb.re;
		  FbIm += FaFb.Fb.im;

		  
		}/* End of loop over detectors */

	      Dt = At * Bt - Ct * Ct;
	      fact = 2.0f / Dt; /* by simplifying the formula, the factor M vanishes, look at the Iraj's note */

	      /* In the signal-only case (only for testing using fake data),
	       * we did not apply any normalization to the data, and we need
	       * the use the correct factor now (taken from CFS_v1)
	       * [see Xavie's notes on LALDemod() for details]
	       */
	      if ( uvar_SignalOnly )
		{
		  REAL8 Tsft = 1.0 / (GV.multiSFTs->data[0]->data[0].deltaF );
 		      
		  fact /= Tsft;
		}
	      else
		{
		  /* NOTE: we normalized the data by a double-sided PSD, (LALNormalizeSFTVect),
		   * therefore we apply another factor of 1/2 now with respect to the 
		   * equations in JKS, which are based on the single-sided PSD:
		   */
		  fact *= 0.5f;
		}
		       
	      /* calculate F-statistic from  Fa and Fb */
	      Fstat = fact * (Bt * (FaRe*FaRe + FaIm*FaIm) 
			      + At * (FbRe*FbRe + FbIm*FbIm) 
			      - 2.0f * Ct *(FaRe*FbRe + FaIm*FbIm) );


	      /* calculate the baysian-marginalized 'B-statistic' */
	      if ( fpBstat )
		{
		  Bstat = exp( Fstat ) / Dt ;
		  fprintf (fpBstat, "%16.12f %8.7f %8.7f %.17g %10.6g\n", 
			   fkdot->data[0], dopplerpos.Alpha, dopplerpos.Delta, fkdot->data[1], 
			   Bstat );
		}
	      
	      /* now, if user requested it, we output ALL F-statistic results above threshold */
	      if ( uvar_outputFstat || uvar_outputLoudest )
		{
		  REAL8 freq = fkdot->data[0];
		  
		  if ( Fstat > uvar_Fthreshold )
		    {
		      LALSnprintf (buf, 511, "%16.12f %8.7f %8.7f %.17g %10.6g\n", 
				   freq, dopplerpos.Alpha, dopplerpos.Delta, 
				   fkdot->data[1], 2.0 * Fstat);
		      buf[511] = 0;
		      if ( fpFstat )
			fprintf (fpFstat, buf );
		    } /* if F > threshold */
		}  
	      
	      if ( Fstat > loudestF )
		{
		  loudestF = Fstat;
		  strcpy ( loudestEntry,  buf );
		}
	      
	    } /* for i < nBins: loop over frequency-bins */

	} /* For GV.spinImax: loop over spindowns */
      

      loopcounter ++;
      if (lalDebugLevel) 
	printf ("\
Search progress: %5.1f%%", (100.0* loopcounter / thisScan.numGridPoints));
      
    } /*  while SkyPos : loop over skypositions */
  
  if ( fpFstat )
    {
      fprintf (fpFstat, "%%DONE\n");
      fclose (fpFstat);
      fpFstat = NULL;
    }
  if ( fpBstat )
    {
      fprintf (fpBstat, "%%DONE\n");
      fclose (fpBstat);
      fpBstat = NULL;
    }

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
  
  XLALDestroyREAL8Vector ( fkdot );

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
  uvar_outputBstat = NULL;

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
  LALregSTRINGUserVar(status,	DataFiles, 	'D', UVAR_REQUIRED, "File-pattern specifying (multi-IFO) input SFT-files"); 
  LALregSTRINGUserVar(status, 	IFO, 		'I', UVAR_OPTIONAL, 
		      "Detector-constraint: 'G1', 'L1', 'H1', 'H2' ...(useful for single-IFO v1-SFTs only!)");
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
  LALregSTRINGUserVar(status,	outputBstat,	 0,  UVAR_OPTIONAL, "Output-file for 'B-statistic' field over the parameter-space");

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

  ASSERT ( edat, status, COMPUTEFSTATISTICC_ENULL, COMPUTEFSTATISTICC_MSGENULL );
  ASSERT ( ephemYear, status, COMPUTEFSTATISTICC_ENULL, COMPUTEFSTATISTICC_MSGENULL );

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
  REAL8 fCoverMin, fCoverMax;	/* covering frequency-band to read from SFTs */
  UINT4 X;		/* index over different detectors */
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
      ABORT ( status,  COMPUTEFSTATISTICC_EINPUT,  COMPUTEFSTATISTICC_MSGEINPUT);
    }

  /* get full SFT-catalog of all matching (multi-IFO) SFTs */
  TRY ( LALSFTdataFind ( status->statusPtr, &catalog, uvar_DataFiles, &constraints ), status);    
  if ( constraints.detector ) 
    LALFree ( constraints.detector );

  if ( catalog->length == 0 ) 
    {
      LALPrintError ("\nSorry, didn't find any matching SFTs with pattern '%s'!\n\n", uvar_DataFiles );
      ABORT ( status,  COMPUTEFSTATISTICC_EINPUT,  COMPUTEFSTATISTICC_MSGEINPUT);
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


  { /* ----- prepare spin-range at refTime (in 'canonical format': Bands >= 0) ----- */
    REAL8 fMin = MYMIN ( uvar_Freq, uvar_Freq + uvar_FreqBand );
    REAL8 fMax = MYMAX ( uvar_Freq, uvar_Freq + uvar_FreqBand );
    REAL8 f1dotMin = MYMIN ( uvar_f1dot, uvar_f1dot + uvar_f1dotBand );
    REAL8 f1dotMax = MYMAX ( uvar_f1dot, uvar_f1dot + uvar_f1dotBand );
    
    if ( ( cfg->spinRangeRef = XLALCreatePulsarSpinRange ( NUM_SPINS )) == NULL ) {
      ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }
    cfg->spinRangeRef->epoch = cfg->refTime;
    cfg->spinRangeRef->fkdot->data[0] = fMin;
    cfg->spinRangeRef->fkdotBand->data[0] = fMax - fMin;
    cfg->spinRangeRef->fkdot->data[1] = f1dotMin;
    cfg->spinRangeRef->fkdotBand->data[1] = f1dotMax - f1dotMin;
  }

  { /* ----- propage spin-range from refTime to startTime and endTime of observation ----- */
    LALPulsarSpinRange *spinRangeEnd;	/* temporary only */
    REAL8 fmaxStart, fmaxEnd, fminStart, fminEnd;

    if ( ( cfg->spinRangeStart = XLALCreatePulsarSpinRange ( NUM_SPINS )) == NULL ) {
      ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }
    if ( ( spinRangeEnd = XLALCreatePulsarSpinRange ( NUM_SPINS )) == NULL ) {
      ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }

    /* compute spin-range at startTime of observation */
    TRY ( LALExtrapolatePulsarSpinRange (status->statusPtr, 
					 cfg->spinRangeStart, cfg->startTime, cfg->spinRangeRef ), status );
    /* compute spin-range at endTime of these SFTs */
    TRY ( LALExtrapolatePulsarSpinRange (status->statusPtr, 
					 spinRangeEnd, endTime, cfg->spinRangeRef ), status );

    fminStart = cfg->spinRangeStart->fkdot->data[0];
    /* ranges are in canonical format! */
    fmaxStart = fminStart + cfg->spinRangeStart->fkdotBand->data[0];  
    fminEnd   = spinRangeEnd->fkdot->data[0];
    fmaxEnd   = fminEnd + spinRangeEnd->fkdotBand->data[0];

    XLALDestroyPulsarSpinRange ( spinRangeEnd );
    /*  get covering frequency-band  */
    fCoverMax = MYMAX ( fmaxStart, fmaxEnd );
    fCoverMin = MYMIN ( fminStart, fminEnd );

  } /* extrapolate spin-range */

  /* ----- correct for maximal doppler-shift due to earth's motion */
 
  
  {/* ----- load the multi-IFO SFT-vectors ----- */
    UINT4 wings = MYMAX(uvar_Dterms, uvar_RngMedWindow/2 +1);	/* extra frequency-bins needed for rngmed, and Dterms */
    REAL8 fMax = (1.0 + uvar_dopplermax) * fCoverMax + wings / Tsft;	/* correct for doppler-shift and wings */
    REAL8 fMin = (1.0 - uvar_dopplermax) * fCoverMin - wings / Tsft;

    TRY ( LALLoadMultiSFTs ( status->statusPtr, &(cfg->multiSFTs), catalog, fMin, fMax ), status );
  }
  TRY ( LALDestroySFTCatalog ( status->statusPtr, &catalog ), status );

  cfg->numDetectors = cfg->multiSFTs->length;

  /*============================== FIXME: the rest is just adapted to keep things working for now: */

  /* ----- allocate array of ifo-specific configuration data ----- */
  if ( ( cfg->ifos = LALCalloc ( cfg->numDetectors, sizeof( *(cfg->ifos) ) )) == NULL ) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  


  for ( X=0; X < cfg->numDetectors; X ++ )
    {
      if ( (cfg->ifos[X].weightsNoise  = XLALCreateREAL8Vector ( cfg->multiSFTs->data[X]->length )) == NULL ) {
	ABORT (status, COMPUTEFSTATISTICC_EMEM, COMPUTEFSTATISTICC_MSGEMEM);
      }
      TRY( LALHOUGHInitializeWeights( status->statusPtr, (cfg->ifos[X].weightsNoise) ), status);
      
      /* Normalize this by 1/sqrt(Tsft * Sh), where Tsft*Sh is the median of |X|^2  
       * NOTE: this corresponds to a double-sided PSD, therefore we need to 
       * divide by another factor of 2 with respect to the JKS formulae.
       */
      if ( ! uvar_SignalOnly ) 
	{
	  TRY( LALComputeNoiseWeights( status->statusPtr, (cfg->ifos[X].weightsNoise), cfg->multiSFTs->data[X], uvar_RngMedWindow), status); 
	  TRY( LALNormalizeSFTVect (status->statusPtr, cfg->multiSFTs->data[X], uvar_RngMedWindow), status );
	}
    } /* for X < numDetectors */



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


  /* Normalyzing the Noise Weights. Sum of weighwNoise over the all detectors is going to be 1.0 */
  {
    UINT4 alpha;
    REAL8 norm = 0.0;
    
    for( X = 0; X < cfg->numDetectors; X ++ )
      {
	UINT4 M;
	REAL8 W = 0.0;
	
	M = cfg->ifos[X].weightsNoise->length;
	
	for (alpha = 0; alpha < M; alpha++)
	  {
	    W += cfg->ifos[X].weightsNoise->data[alpha]; 
	  }
	norm += W / M;
      } /* for X < numDetector */

    for ( X = 0; X < cfg->numDetectors; X++ )
      for (alpha = 0; alpha < cfg->ifos[X].weightsNoise->length; alpha++ )
	cfg->ifos[X].weightsNoise->data[alpha] /= norm;
    
  } /* normalize noise-weights */

  /* ----- initialize + allocate space for AM-coefficients and SSB-times ----- */
  {
    AMCoeffs *amc = NULL;
    SSBtimes *tSSB = NULL;

    for ( X = 0; X < cfg->numDetectors; X ++ )
      {
	/* Allocate space for AMCoeffs */
	if ( (amc = LALCalloc(1, sizeof(AMCoeffs))) == NULL) {
	  ABORT (status, COMPUTEFSTATISTICC_EMEM, COMPUTEFSTATISTICC_MSGEMEM);
	}
	TRY (LALSCreateVector(status->statusPtr, &(amc->a), cfg->multiSFTs->data[X]->length ), status);
	TRY (LALSCreateVector(status->statusPtr, &(amc->b), cfg->multiSFTs->data[X]->length ), status);
    
	cfg->ifos[X].amcoe = amc;
	amc = NULL;
    
	/* allocate memory of the SSB-times: DeltaT_alpha and Tdot_alpha */
	if ( (tSSB = LALCalloc(1, sizeof(SSBtimes))) == NULL) {
	  ABORT (status, COMPUTEFSTATISTICC_EMEM, COMPUTEFSTATISTICC_MSGEMEM);
	}
	TRY(LALDCreateVector(status->statusPtr, &(tSSB->DeltaT), cfg->multiSFTs->data[X]->length),status);
	TRY(LALDCreateVector(status->statusPtr, &(tSSB->Tdot),  cfg->multiSFTs->data[X]->length), status );

	cfg->ifos[X].tSSB = tSSB;
	tSSB = NULL;
      } /* for X < numDetectors */

  } /* init AM- and SSBtimes */

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


/** Free all globally allocated memory. */
void
Freemem(LALStatus *status,  ConfigVariables *cfg) 
{
  UINT4 X;
  INITSTATUS (status, "Freemem", rcsid);
  ATTATCHSTATUSPTR (status);


  /* Free SFT data */
  TRY (LALDestroyMultiSFTVector (status->statusPtr, &(cfg->multiSFTs) ), status );

  for ( X=0; X < cfg->numDetectors; X++)
    {
      /* Free Noise Weight */
      TRY (LALDDestroyVector(status->statusPtr, &(cfg->ifos[X].weightsNoise)), status); 
                  
      /* Free AM-coefficients */
      TRY (LALSDestroyVector(status->statusPtr, &(cfg->ifos[X].amcoe->a)), status);
      TRY (LALSDestroyVector(status->statusPtr, &(cfg->ifos[X].amcoe->b)), status);
      LALFree ( cfg->ifos[X].amcoe);
      /* Free SSB-times */
      TRY (LALDDestroyVector(status->statusPtr, &(cfg->ifos[X].tSSB->DeltaT)), status);
      TRY (LALDDestroyVector(status->statusPtr, &(cfg->ifos[X].tSSB->Tdot)), status);
      LALFree ( cfg->ifos[X].tSSB );

      /* destroy DetectorStateSeries */
    }  

  XLALDestroyMultiDetectorStateSeries ( cfg->multiDetStates );

  LALFree ( cfg->ifos );

  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (status->statusPtr), status);

  if ( GV.searchRegion.skyRegionString )
    LALFree ( GV.searchRegion.skyRegionString );
  
  /* Free ephemeris data */
  LALFree(cfg->edat->ephemE);
  LALFree(cfg->edat->ephemS);
  LALFree(cfg->edat);

  XLALDestroyPulsarSpinRange ( cfg->spinRangeStart );
  XLALDestroyPulsarSpinRange ( cfg->spinRangeRef );

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
