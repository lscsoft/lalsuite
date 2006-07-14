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
/** \author R. Prix, I. Gholami, Y. Ioth, Papa, X. Siemens
 * \file 
 * \brief
 * Calculate the F-statistic for a given parameter-space of pulsar GW signals.
 * Implements the so-called "F-statistic" as introduced in \ref JKS98.
 *                                                                          
 *********************************************************************************/
#include "config.h"

/* System includes */
#include <math.h>
#include <stdio.h>
#include <time.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif


int finite(double);

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
#include "LogPrintf.h"

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
#define SQ(x) ( (x) * (x) )

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define LAL_INT4_MAX 2147483647

/*---------- internal types ----------*/

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  LIGOTimeGPS startTime;		    /**< start time of observation */
  REAL8 Tsft;                               /**< length of one SFT in seconds */
  REAL8 duration;			    /**< total time-span of the data (all streams) in seconds */
  LIGOTimeGPS refTime;			    /**< reference-time for pulsar-parameters in SBB frame */
  LALPulsarSpinRange *spinRangeRef; 	    /**< pulsar spin-range at reference-time 'refTime' */
  LALPulsarSpinRange *spinRangeStart; 	    /**< pulsar spin-range at start of observation 'startTime; */
  DopplerRegion searchRegion;		    /**< parameter-space region to search over (FIXME) */
  EphemerisData *edat;			    /**< ephemeris data (from LALInitBarycenter()) */
  MultiSFTVector *multiSFTs;		    /**< multi-IFO SFT-vectors */
  MultiDetectorStateSeries *multiDetStates; /**< pos, vel and LMSTs for detector at times t_i */
  MultiNoiseWeights *multiNoiseWeights;	    /**< normalized noise-weights of those SFTs */
  REAL8 S_hat;                              /**< Sum over the 1/Sn */ 
  ComputeFParams CFparams;		    /**< parameters for Fstat (e.g Dterms, SSB-prec,...) */
  CHAR *logstring;                          /**< log containing max-info on this search setup */
} ConfigVariables;


/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

ConfigVariables GV;		/**< global container for various derived configuration settings */

/* ----- User-variables: can be set from config-file or command-line */
INT4 uvar_Dterms;
CHAR *uvar_IFO;
BOOLEAN uvar_SignalOnly;
BOOLEAN uvar_noHeader;
BOOLEAN uvar_addOutput;
REAL8 uvar_Freq;
REAL8 uvar_FreqBand;
REAL8 uvar_dFreq;
REAL8 uvar_Alpha;
REAL8 uvar_dAlpha;
REAL8 uvar_AlphaBand;
REAL8 uvar_Delta;
REAL8 uvar_dDelta;
REAL8 uvar_DeltaBand;
/* 1st spindown */
REAL8 uvar_f1dot;
REAL8 uvar_df1dot;
REAL8 uvar_f1dotBand;
/* 2nd spindown */
REAL8 uvar_f2dot;
REAL8 uvar_df2dot;
REAL8 uvar_f2dotBand;
/* 3rd spindown */
REAL8 uvar_f3dot;
REAL8 uvar_df3dot;
REAL8 uvar_f3dotBand;
/* --- */
REAL8 uvar_TwoFthreshold;
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
REAL8 uvar_dopplermax;
INT4 uvar_RngMedWindow;
REAL8 uvar_refTime;
INT4 uvar_SSBprecision;

INT4 uvar_minStartTime;
INT4 uvar_maxEndTime;
CHAR *uvar_workingDir;
REAL8 uvar_timerCount;

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

int XLALwriteCandidate2file ( FILE *fp,  const PulsarCandidate *cand, const Fcomponents *Fstat );


const char *va(const char *format, ...);	/* little var-arg string helper function */

/*---------- empty initializers ---------- */
static const PulsarTimesParamStruc empty_PulsarTimesParamStruc;
static const BarycenterInput empty_BarycenterInput;
static const SFTConstraints empty_SFTConstraints;
static const ComputeFBuffer empty_ComputeFBuffer;
static const LIGOTimeGPS empty_LIGOTimeGPS;
static const PulsarCandidate empty_PulsarCandidate;
static const PulsarDopplerParams empty_DopplerParams;
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
  SkyPosition skypos;
  FILE *fpFstat = NULL;
  FILE *fpBstat = NULL;
  CHAR buf[512];
  REAL8Vector *fkdotRef = NULL;
  REAL8Vector *fkdotTmp = NULL;
  ComputeFBuffer cfBuffer = empty_ComputeFBuffer;
  UINT4 nFreq, nf1dot, nf2dot, nf3dot;	/* number of frequency- and f1dot-bins */
  UINT4 iFreq, if1dot, if2dot, if3dot;  /* counters over freq- and f1dot- bins */
  REAL8 numTemplates, templateCounter;
  REAL8 tickCounter;
  time_t clock0;
  Fcomponents Fstat, loudestFstat;
  PulsarDopplerParams doppler = empty_DopplerParams;
  PulsarDopplerParams loudestDoppler = empty_DopplerParams;


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

  /* set log-level */
  LogSetLevel ( lalDebugLevel );

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
  /* this is a bit obsolete, but DopplerScan still uses it..: */
  GV.searchRegion.Freq = GV.spinRangeStart->fkdot->data[0];
  GV.searchRegion.FreqBand = GV.spinRangeStart->fkdotBand->data[0];
  GV.searchRegion.f1dot = GV.spinRangeStart->fkdot->data[1];
  GV.searchRegion.f1dotBand = GV.spinRangeStart->fkdotBand->data[1];
  scanInit.searchRegion = GV.searchRegion;
  
  LogPrintf (LOG_DEBUG, "Setting up template grid ... ");
  LAL_CALL ( InitDopplerScan ( &status, &thisScan, &scanInit), &status); 
  LogPrintfVerbatim (LOG_DEBUG, "done.\n");
  
  /* ---------- set Frequency- and spindown-resolution if not input by user ----------*/
  if ( LALUserVarWasSet( &uvar_dFreq ) )
    thisScan.dFreq = uvar_dFreq;
  
  if( LALUserVarWasSet( &uvar_df1dot) ) 
    thisScan.df1dot = uvar_df1dot;

  LogPrintf (LOG_DEBUG, "Actual grid-spacings: dFreq = %g, df1dot = %g\n\n", 
	     thisScan.dFreq, thisScan.df1dot );
  
  /*----------------------------------------------------------------------*/
  if ( uvar_outputSkyGrid ) {
    LogPrintf (LOG_NORMAL, "Now writing sky-grid into file '%s' ...", uvar_outputSkyGrid);
    LAL_CALL (writeSkyGridFile( &status, thisScan.grid, uvar_outputSkyGrid, &scanInit), &status);
    LogPrintfVerbatim (LOG_NORMAL, " done.\n");
  }
  
  /* if a complete output of the F-statistic file was requested,
   * we open and prepare the output-file here */
  if (uvar_outputFstat)
    {
      if(uvar_addOutput)
	{
	  if ( (fpFstat = fopen (uvar_outputFstat, "ab")) == NULL)
	    {
	      LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputFstat);
	      return (COMPUTEFSTATISTIC_ESYS);
	    }
	}
      else
	{
	  if ( (fpFstat = fopen (uvar_outputFstat, "wb")) == NULL)
	    {
	      LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputFstat);
	      return (COMPUTEFSTATISTIC_ESYS);
	    }
	}
      
      if(!uvar_noHeader)
	fprintf (fpFstat, "%s", GV.logstring );
    } /* if outputFstat */
  
  if (uvar_outputBstat)
    {
      if ( (fpBstat = fopen (uvar_outputBstat, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputBstat);
	  return (COMPUTEFSTATISTIC_ESYS);
	}
    } /* if outputFstat */
  
  if (uvar_outputBstat)
    {
      if ( (fpBstat = fopen (uvar_outputBstat, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputBstat);
	  return (COMPUTEFSTATISTIC_ESYS);
	}
    } /* if outputBstat */

  if ( ( fkdotRef = XLALCreateREAL8Vector ( NUM_SPINS ) ) == NULL ) {
    return COMPUTEFSTATISTIC_EMEM;
  }
  if ( ( fkdotTmp = XLALCreateREAL8Vector ( NUM_SPINS ) ) == NULL ) {
    return COMPUTEFSTATISTIC_EMEM;
  }

  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
   */
  
  doppler.refTime = GV.startTime;	/* we compute at startTime, not refTime right now */  
  doppler.orbit = NULL;			/* binary pulsars not implemented yet */
  /* loop-counters for spin-loops fkdot */
  nFreq =  (UINT4)(GV.spinRangeStart->fkdotBand->data[0] / thisScan.dFreq  + 0.5) + 1;  
  nf1dot = (UINT4)(GV.spinRangeStart->fkdotBand->data[1] / thisScan.df1dot + 0.5) + 1; 
  
  /* the 2nd and 3rd spindown stepsizes are not controlled by DopplerScan (and the metric) yet */
  nf2dot = (UINT4)(GV.spinRangeStart->fkdotBand->data[2] / uvar_df2dot + 0.5) + 1; 
  nf3dot = (UINT4)(GV.spinRangeStart->fkdotBand->data[3] / uvar_df3dot + 0.5) + 1; 

  numTemplates = 1.0 * thisScan.numGridPoints * nFreq * nf1dot * nf2dot * nf3dot;
  
  LogPrintf (LOG_DEBUG, "N = Sky x Freq x f1dot x f2dot x f3dot = %d x %d x %d x %d x %d = %g\n",
	     thisScan.numGridPoints, nFreq, nf1dot, nf2dot, nf3dot, 
	     numTemplates);

  LogPrintf (LOG_DEBUG, "Progress: 0/%g = 0 %% done, Estimated time left: ?? s\n", numTemplates );

  templateCounter = 0.0; 
  tickCounter = 0;
  clock0 = time(NULL);
  while (1)
    {
      LAL_CALL (NextDopplerPos( &status, &dopplerpos, &thisScan ), &status);
      if (thisScan.state == STATE_FINISHED) /* scanned all DopplerPositions yet? */
	break;
      
      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
      skypos.longitude = dopplerpos.Alpha;
      skypos.latitude = dopplerpos.Delta;
      skypos.system = COORDINATESYSTEM_EQUATORIAL;
      LAL_CALL (LALNormalizeSkyPosition(&status, &skypos, &skypos), &status);
      /* set parameter-space point: sky-position */
      doppler.Alpha = skypos.longitude;
      doppler.Delta = skypos.latitude;

      /*----- loop over spindown values */
      for ( if3dot = 0; if3dot < nf3dot; if3dot ++ )
	{
	  doppler.fkdot[3] = GV.spinRangeStart->fkdot->data[3] + if3dot * uvar_df3dot;
	  
	  for ( if2dot = 0; if2dot < nf2dot; if2dot ++ )
	    {
	      doppler.fkdot[2] = GV.spinRangeStart->fkdot->data[2] + if2dot * uvar_df2dot;

	      for (if1dot = 0; if1dot < nf1dot; if1dot ++)
		{
		  doppler.fkdot[1] = GV.spinRangeStart->fkdot->data[1] + if1dot * thisScan.df1dot;
	  
		  /* Loop over frequencies to be demodulated */
		  for ( iFreq = 0 ; iFreq < nFreq ; iFreq ++ )
		    {
		      doppler.fkdot[0] = GV.spinRangeStart->fkdot->data[0] + iFreq * thisScan.dFreq;
		      
		      LAL_CALL( ComputeFStat(&status, &Fstat, &doppler, GV.multiSFTs, GV.multiNoiseWeights, 
					     GV.multiDetStates, &GV.CFparams, &cfBuffer ), &status );

		      templateCounter += 1.0;
		      tickCounter += 1.0;
		      if ( lalDebugLevel && ( tickCounter > uvar_timerCount) )
			{
			  REAL8 diffSec = time(NULL) - clock0 ;  /* seconds since start of loop*/
			  REAL8 taup = diffSec / templateCounter ;
			  REAL8 timeLeft = (numTemplates - templateCounter) *  taup;
			  tickCounter = 0.0;
			  LogPrintf (LOG_DEBUG, 
				     "Progress: %g/%g = %.2f %% done, Estimated time left: %.0f s\n",
				     templateCounter, numTemplates, 
				     templateCounter/numTemplates * 100.0, timeLeft);
			}

		      if ( !finite(Fstat.F) )
			{
			  LogPrintf(LOG_CRITICAL, 
				    "non-finite F = %.16g, Fa=(%.16g,%.16g), Fb=(%.16g,%.16g)\n", 
				    Fstat.F, Fstat.Fa.re, Fstat.Fa.im, Fstat.Fb.re, Fstat.Fb.im );
			  LogPrintf (LOG_CRITICAL, 
				     "[Alpha,Delta] = [%.16g,%.16g],\n"
				     "fkdot=[%.16g,%.16g,%.16g,%16.g]\n",
				     dopplerpos.Alpha, dopplerpos.Delta, fkdotRef->data[0], 
				     fkdotRef->data[1], fkdotRef->data[2], fkdotRef->data[3] );
			  return -1;
			}

		      /* correct results in --SignalOnly case:
		       * this means we didn't normalize data by 1/sqrt(Tsft * 0.5 * Sh) in terms of 
		       * the single-sided PSD Sh: the SignalOnly case is characterized by
		         
		       * setting Sh->1, so we need to divide Fa,Fb by sqrt(0.5*Tsft)
		       * and F by (0.5*Tsft)
		       */
		      if ( uvar_SignalOnly )
			{
			  REAL8 norm = 1.0 / sqrt( 0.5 * GV.Tsft );
			  Fstat.Fa.re *= norm;
			  Fstat.Fa.im *= norm;
			  Fstat.Fb.re *= norm;
			  Fstat.Fb.im *= norm;
			  Fstat.F *= norm * norm;
			} /* if SignalOnly */
		      
		      /* propagate fkdot back to reference-time for outputting results */
		      fkdotTmp->data[0] = doppler.fkdot[0];
		      fkdotTmp->data[1] = doppler.fkdot[1];
		      fkdotTmp->data[2] = doppler.fkdot[2];
		      fkdotTmp->data[3] = doppler.fkdot[3];

		      LAL_CALL ( LALExtrapolatePulsarSpins(&status, fkdotRef, GV.refTime, fkdotTmp, 
							   GV.startTime ), &status );

		      /* calculate the baysian-marginalized 'B-statistic' */
		      if ( fpBstat )
			{
			  fprintf (fpBstat, "%16.12f %8.7f %8.7f %.17g %10.6g\n", 
				   fkdotRef->data[0], dopplerpos.Alpha, dopplerpos.Delta, fkdotRef->data[1], 
				   Fstat.Bstat );
			}
		      
		      /* now, if user requested it, we output ALL F-statistic results above threshold */
		      if ( uvar_outputFstat && ( 2.0 * Fstat.F >= uvar_TwoFthreshold ) )
			{
 			  LALSnprintf (buf, 511, "%.16g %.16g %.16g %.6g %.5g %.5g %.9g\n",
				       fkdotRef->data[0], 
				       dopplerpos.Alpha, dopplerpos.Delta, 
				       fkdotRef->data[1], fkdotRef->data[2], fkdotRef->data[3],
				       2.0 * Fstat.F );
			  buf[511] = 0;
			  if ( fpFstat )
			    fprintf (fpFstat, buf );
			} /* if F > threshold */

		      /* keep track of loudest candidate */
		      if ( Fstat.F > loudestFstat.F )
			{
			  loudestFstat = Fstat;
			  loudestDoppler = doppler;
			  /* correct spins to reference-time */
			  loudestDoppler.fkdot[0] = fkdotRef->data[0];
			  loudestDoppler.fkdot[1] = fkdotRef->data[1];
			  loudestDoppler.fkdot[2] = fkdotRef->data[2];
			  loudestDoppler.fkdot[3] = fkdotRef->data[3];
			}
		    } /* for i < nBins: loop over frequency-bins */
		} /* For GV.spinImax: loop over 1st spindowns */
	    } /* for if2dot < nf2dot */
	} /* for if3dot < nf3dot */
    } /*  while SkyPos : loop over skypositions */

 
  if ( fpFstat )
    {
      if(!uvar_noHeader)
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

  /* do full parameter-estimation for loudest canidate and output into separate file */
  if ( uvar_outputLoudest )
    {
      MultiAMCoeffs *multiAMcoef = NULL;
      REAL8 norm = GV.Tsft * GV.S_hat;
      FILE *fpLoudest;
      PulsarAmplitudeParams Amp, dAmp;
      PulsarCandidate cand = empty_PulsarCandidate;

      skypos.longitude = loudestDoppler.Alpha;
      skypos.latitude  = loudestDoppler.Delta;
      skypos.system = COORDINATESYSTEM_EQUATORIAL;
      LAL_CALL ( LALGetMultiAMCoeffs ( &status, &multiAMcoef, GV.multiDetStates, skypos ), &status);
      /* noise-weigh Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiAMCoeffs ( multiAMcoef, GV.multiNoiseWeights ) != XLAL_SUCCESS ) {
	LALPrintError("\nXLALWeighMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	return COMPUTEFSTATC_EXLAL;
      }

      LAL_CALL(LALEstimatePulsarAmplitudeParams (&status, &Amp, &dAmp, &loudestFstat, cfBuffer.multiAMcoef, norm), 
	       &status);
      XLALDestroyMultiAMCoeffs ( multiAMcoef );
      
      /* propagate initial-phase from internal reference-time 'startTime' to refTime of Doppler-params */
      fkdotTmp->data[0] = loudestDoppler.fkdot[0];
      fkdotTmp->data[1] = loudestDoppler.fkdot[1];
      fkdotTmp->data[2] = loudestDoppler.fkdot[2];
      fkdotTmp->data[3] = loudestDoppler.fkdot[3];
      LAL_CALL(LALExtrapolatePulsarPhase (&status, &Amp.phi0, fkdotTmp, GV.refTime, Amp.phi0, GV.startTime),&status);
      if ( Amp.phi0 < 0 )	      /* make sure phi0 in [0, 2*pi] */
	Amp.phi0 += LAL_TWOPI;
      Amp.phi0 = fmod ( Amp.phi0, LAL_TWOPI );

      if(uvar_addOutput)
	{
	  if ( (fpLoudest = fopen (uvar_outputLoudest, "ab")) == NULL)
	    {
	      LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputLoudest);
	      return COMPUTEFSTATISTIC_ESYS;
	    }
	}
      else
	{
	  if ( (fpLoudest = fopen (uvar_outputLoudest, "wb")) == NULL)
	    {
	      LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputLoudest);
	      return COMPUTEFSTATISTIC_ESYS;
	    }
	}
      /* write header with run-info */
      if(!uvar_noHeader)
	fprintf (fpLoudest, "%s", GV.logstring );

      /* assemble 'candidate' structure */
      cand.Amp = Amp;
      cand.dAmp = dAmp;
      cand.Doppler = loudestDoppler;
      if ( XLALwriteCandidate2file ( fpLoudest,  &cand, &loudestFstat ) != XLAL_SUCCESS )
	{
	  LALPrintError("\nXLALwriteCandidate2file() failed with error = %d\n\n", xlalErrno );
	  return COMPUTEFSTATC_EXLAL;
	}
      fclose (fpLoudest);
    } /* write loudest candidate to file */
  
  LogPrintf (LOG_DEBUG, "Search finished.\n");
  
  /* Free memory */
  LAL_CALL ( FreeDopplerScan(&status, &thisScan), &status);

  XLALEmptyComputeFBuffer ( cfBuffer );

  XLALDestroyREAL8Vector ( fkdotRef );
  XLALDestroyREAL8Vector ( fkdotTmp );

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
  uvar_noHeader = FALSE;
  uvar_addOutput =FALSE;

  uvar_f1dot     = 0.0;
  uvar_df1dot    = 1.0;
  uvar_df2dot    = 1.0;
  uvar_df3dot    = 1.0;
  uvar_f1dotBand = 0.0;
  
  uvar_TwoFthreshold = 10.0;
  uvar_metricType =  LAL_PMETRIC_NONE;
  uvar_gridType = GRID_FLAT;

  uvar_metricMismatch = 0.02;

  uvar_help = FALSE;
  uvar_outputLabel = NULL;

  uvar_outputFstat = NULL;
  uvar_outputBstat = NULL;

  uvar_skyGridFile = NULL;

  uvar_dopplermax =  1.05e-4;
  uvar_RngMedWindow = 50;	/* for running-median */

  uvar_SSBprecision = SSBPREC_RELATIVISTIC;

  uvar_minStartTime = 0;
  uvar_maxEndTime = LAL_INT4_MAX;

  uvar_workingDir = (CHAR*)LALMalloc(512);
  strcpy(uvar_workingDir, ".");

  uvar_timerCount = 1e4;	/* output a timer/progress count every N templates */

  /* ---------- register all user-variables ---------- */
  LALregBOOLUserVar(status, 	help, 		'h', UVAR_HELP,     "Print this message"); 

  LALregREALUserVar(status, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	Freq, 		'f', UVAR_REQUIRED, "Starting search frequency in Hz");
  LALregREALUserVar(status, 	f1dot, 		's', UVAR_OPTIONAL, "First spindown parameter  dFreq/dt");
  LALregREALUserVar(status, 	f2dot, 		 0 , UVAR_OPTIONAL, "Second spindown parameter d^2Freq/dt^2");
  LALregREALUserVar(status, 	f3dot, 		 0 , UVAR_OPTIONAL, "Third spindown parameter  d^3Freq/dt^2");

  LALregREALUserVar(status, 	AlphaBand, 	'z', UVAR_OPTIONAL, "Band in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	DeltaBand, 	'c', UVAR_OPTIONAL, "Band in delta (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	FreqBand, 	'b', UVAR_OPTIONAL, "Search frequency band in Hz");
  LALregREALUserVar(status, 	f1dotBand, 	'm', UVAR_OPTIONAL, "Search-band for f1dot");
  LALregREALUserVar(status, 	f2dotBand, 	 0 , UVAR_OPTIONAL, "Search-band for f2dot");
  LALregREALUserVar(status, 	f3dotBand, 	 0 , UVAR_OPTIONAL, "Search-band for f3dot");

  LALregREALUserVar(status, 	dAlpha, 	'l', UVAR_OPTIONAL, "Resolution in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	dDelta, 	'g', UVAR_OPTIONAL, "Resolution in delta (equatorial coordinates) in radians");
  LALregREALUserVar(status,     dFreq,          'r', UVAR_OPTIONAL, "Frequency resolution in Hz (default: 1/( 2 * Tobs )");
  LALregREALUserVar(status, 	df1dot, 	'e', UVAR_OPTIONAL, "Stepsize for f1dot (default 1/( 2 * Tobs^2 )");
  LALregREALUserVar(status, 	df2dot, 	 0 , UVAR_OPTIONAL, "Stepsize for f2dot");
  LALregREALUserVar(status, 	df3dot, 	 0 , UVAR_OPTIONAL, "Stepsize for f3dot");

  LALregSTRINGUserVar(status,	skyRegion, 	'R', UVAR_OPTIONAL, "ALTERNATIVE: Specify sky-region by polygon (or use 'allsky')");
  LALregSTRINGUserVar(status,	DataFiles, 	'D', UVAR_REQUIRED, "File-pattern specifying (multi-IFO) input SFT-files"); 
  LALregSTRINGUserVar(status, 	IFO, 		'I', UVAR_OPTIONAL, 
		      "Detector-constraint: 'G1', 'L1', 'H1', 'H2' ...(useful for single-IFO v1-SFTs only!)");
  LALregSTRINGUserVar(status,	ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(status,	ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregBOOLUserVar(status, 	SignalOnly, 	'S', UVAR_OPTIONAL, "Signal only flag");
  LALregREALUserVar(status, 	TwoFthreshold,	'F', UVAR_OPTIONAL, "Set the threshold for selection of 2F");
  LALregINTUserVar(status, 	gridType,	 0 , UVAR_OPTIONAL, "Template grid: 0=flat, 1=isotropic, 2=metric, 3=file");
  LALregINTUserVar(status, 	metricType,	'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserVar(status, 	metricMismatch,	'X', UVAR_OPTIONAL, "Maximal allowed mismatch for metric tiling");
  LALregSTRINGUserVar(status,	outputLabel,	'o', UVAR_OPTIONAL, "Label to be appended to all output file-names");
  LALregSTRINGUserVar(status,	skyGridFile,	 0,  UVAR_OPTIONAL, "Load sky-grid from this file.");
  LALregREALUserVar(status,	refTime,	 0,  UVAR_OPTIONAL, "SSB reference time for pulsar-paramters");
  LALregREALUserVar(status, 	dopplermax, 	'q', UVAR_OPTIONAL, "Maximum doppler shift expected");  
  LALregSTRINGUserVar(status,	outputFstat,	 0,  UVAR_OPTIONAL, "Output-file for F-statistic field over the parameter-space");
  LALregSTRINGUserVar(status,	outputBstat,	 0,  UVAR_OPTIONAL, "Output-file for 'B-statistic' field over the parameter-space");

  LALregINTUserVar ( status, 	minStartTime, 	 0,  UVAR_OPTIONAL, "Earliest SFT-timestamp to include");
  LALregINTUserVar ( status, 	maxEndTime, 	 0,  UVAR_OPTIONAL, "Latest SFT-timestamps to include");

  /* ----- more experimental/expert options follow here ----- */
  LALregINTUserVar (status, 	SSBprecision,	 0,  UVAR_DEVELOPER, "Precision to use for time-transformation to SSB: 0=Newtonian 1=relativistic");
  LALregINTUserVar(status, 	RngMedWindow,	'k', UVAR_DEVELOPER, "Running-Median window size");
  LALregINTUserVar(status,	Dterms,		't', UVAR_DEVELOPER, "Number of terms to keep in Dirichlet kernel sum");
  LALregSTRINGUserVar(status,	outputSkyGrid,	 0,  UVAR_DEVELOPER, "Write sky-grid into this file.");
  LALregSTRINGUserVar(status,   outputLoudest,	 0,  UVAR_DEVELOPER, 
		      "Output-file for the loudest F-statistic candidate in this search");

  LALregSTRINGUserVar(status,   workingDir,     'w', UVAR_DEVELOPER, "Directory to use as work directory.");
  LALregREALUserVar(status, 	timerCount, 	 0,  UVAR_DEVELOPER, "N: Output progress/timer info every N templates");  
  LALregBOOLUserVar(status, 	noHeader, 	'H', UVAR_DEVELOPER, "Do not print the header if this option is used");
  LALregBOOLUserVar(status, 	addOutput, 	'A', UVAR_DEVELOPER, "Add output result to the previous run");


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
  REAL8 fCoverMin, fCoverMax;	/* covering frequency-band to read from SFTs */
  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_SFTConstraints;
  LIGOTimeGPS minStartTimeGPS = empty_LIGOTimeGPS;
  LIGOTimeGPS maxEndTimeGPS = empty_LIGOTimeGPS;

  UINT4 numSFTs;
  LIGOTimeGPS endTime;

  INITSTATUS (status, "InitFStat", rcsid);
  ATTATCHSTATUSPTR (status);


  /* set the current working directory */
  if(chdir(uvar_workingDir) != 0)
    {
      LogPrintf (LOG_CRITICAL,  "Unable to change directory to workinDir '%s'\n", uvar_workingDir);
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }

  /* use IFO-contraint if one given by the user */
  if ( LALUserVarWasSet ( &uvar_IFO ) )
    if ( (constraints.detector = XLALGetChannelPrefix ( uvar_IFO )) == NULL ) {
      ABORT ( status,  COMPUTEFSTATISTIC_EINPUT,  COMPUTEFSTATISTIC_MSGEINPUT);
    }
  minStartTimeGPS.gpsSeconds = uvar_minStartTime;
  maxEndTimeGPS.gpsSeconds = uvar_maxEndTime;
  constraints.startTime = &minStartTimeGPS;
  constraints.endTime = &maxEndTimeGPS;
  
  /* get full SFT-catalog of all matching (multi-IFO) SFTs */
  LogPrintf (LOG_DEBUG, "Finding all SFTs to load ... ");
  TRY ( LALSFTdataFind ( status->statusPtr, &catalog, uvar_DataFiles, &constraints ), status);    
  LogPrintfVerbatim (LOG_DEBUG, "done. (found %d SFTs)\n", catalog->length);

  if ( constraints.detector ) 
    LALFree ( constraints.detector );

  if ( !catalog || catalog->length == 0 ) 
    {
      LALPrintError ("\nSorry, didn't find any matching SFTs with pattern '%s'!\n\n", 
		     uvar_DataFiles );
      ABORT ( status,  COMPUTEFSTATISTIC_EINPUT,  COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* deduce start- and end-time of the observation spanned by the data */
  numSFTs = catalog->length;
  cfg->Tsft = 1.0 / catalog->data[0].header.deltaF;
  cfg->startTime = catalog->data[0].header.epoch;
  endTime   = catalog->data[numSFTs-1].header.epoch;
  LALAddFloatToGPS(status->statusPtr, &endTime, &endTime, cfg->Tsft );	/* can't fail */
  cfg->duration = GPS2REAL8(endTime) - GPS2REAL8 (cfg->startTime);

  /* ----- get reference-time (from user if given, use startTime otherwise): ----- */
  if ( LALUserVarWasSet(&uvar_refTime)) {
    TRY ( LALFloatToGPS (status->statusPtr, &(cfg->refTime), &uvar_refTime), status);
  } else
    cfg->refTime = cfg->startTime;

  { /* ----- prepare spin-range at refTime (in *canonical format*, ie all Bands >= 0) ----- */
    REAL8 fMin = MYMIN ( uvar_Freq, uvar_Freq + uvar_FreqBand );
    REAL8 fMax = MYMAX ( uvar_Freq, uvar_Freq + uvar_FreqBand );

    REAL8 f1dotMin = MYMIN ( uvar_f1dot, uvar_f1dot + uvar_f1dotBand );
    REAL8 f1dotMax = MYMAX ( uvar_f1dot, uvar_f1dot + uvar_f1dotBand );

    REAL8 f2dotMin = MYMIN ( uvar_f2dot, uvar_f2dot + uvar_f2dotBand );
    REAL8 f2dotMax = MYMAX ( uvar_f2dot, uvar_f2dot + uvar_f2dotBand );

    REAL8 f3dotMin = MYMIN ( uvar_f3dot, uvar_f3dot + uvar_f3dotBand );
    REAL8 f3dotMax = MYMAX ( uvar_f3dot, uvar_f3dot + uvar_f3dotBand );
    
    if ( ( cfg->spinRangeRef = XLALCreatePulsarSpinRange ( NUM_SPINS )) == NULL ) {
      ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }
    cfg->spinRangeRef->epoch = cfg->refTime;
    cfg->spinRangeRef->fkdot->data[0] = fMin;
    cfg->spinRangeRef->fkdotBand->data[0] = fMax - fMin;
    cfg->spinRangeRef->fkdot->data[1] = f1dotMin;
    cfg->spinRangeRef->fkdotBand->data[1] = f1dotMax - f1dotMin;
    cfg->spinRangeRef->fkdot->data[2] = f2dotMin;
    cfg->spinRangeRef->fkdotBand->data[2] = f2dotMax - f2dotMin;
    cfg->spinRangeRef->fkdot->data[3] = f3dotMin;
    cfg->spinRangeRef->fkdotBand->data[3] = f3dotMax - f3dotMin;
  } /* spin-range at refTime */

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

  { /* ----- propagate spin-range from refTime to startTime and endTime of observation ----- */
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
    REAL8 fMax = (1.0 + uvar_dopplermax) * fCoverMax + wings / cfg->Tsft; /* correct for doppler-shift and wings */
    REAL8 fMin = (1.0 - uvar_dopplermax) * fCoverMin - wings / cfg->Tsft;
    
    LogPrintf (LOG_DEBUG, "Loading SFTs ... ");
    TRY ( LALLoadMultiSFTs ( status->statusPtr, &(cfg->multiSFTs), catalog, fMin, fMax ), status );
    LogPrintfVerbatim (LOG_DEBUG, "done.\n");
    TRY ( LALDestroySFTCatalog ( status->statusPtr, &catalog ), status );
  }

  /* ----- normalize SFTs and calculate noise-weights ----- */
  if ( uvar_SignalOnly )
    {
      cfg->multiNoiseWeights = NULL; 
      GV.S_hat = 2;
    }
  else
    {
      MultiPSDVector *psds = NULL;

      TRY ( LALNormalizeMultiSFTVect (status->statusPtr, &psds, cfg->multiSFTs, uvar_RngMedWindow ), status );
      /* note: the normalization S_hat would be required to compute the ML-estimator for A^\mu */
      TRY ( LALComputeMultiNoiseWeights  (status->statusPtr, &(cfg->multiNoiseWeights), &GV.S_hat, psds, uvar_RngMedWindow, 0 ), status );

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

  /* ----- set computational parameters for F-statistic from User-input ----- */
  cfg->CFparams.Dterms = uvar_Dterms;
  cfg->CFparams.SSBprec = uvar_SSBprecision;


  /* ----- produce a log-string describing the data-specific setup ----- */
  {
    struct tm utc;
    time_t tp;
    CHAR dateStr[512], line[512], summary[4096];
    CHAR *cmdline = NULL;
    UINT4 i, numDet, numSpins;
    const CHAR *codeID = "$Id$";

    /* first get full commandline describing search*/
    TRY ( LALUserVarGetLog (status->statusPtr, &cmdline,  UVAR_LOGFMT_CMDLINE ), status );
    sprintf (summary, "%%%% %s\n%%%% %s\n", codeID, cmdline );
    LALFree ( cmdline );

    numDet = cfg->multiSFTs->length;
    tp = time(NULL);
    sprintf (line, "%%%% Started search: %s", asctime( gmtime( &tp ) ) );
    strcat ( summary, line );
    strcat (summary, "%% Loaded SFTs: [ " );
    for ( i=0; i < numDet; i ++ ) 
      {
	sprintf (line, "%s:%d%s",  cfg->multiSFTs->data[i]->data->name, 
		 cfg->multiSFTs->data[i]->length,
		 (i < numDet - 1)?", ":" ]\n");
	strcat ( summary, line );
      }
    utc = *XLALGPSToUTC( &utc, (INT4)GPS2REAL8(cfg->startTime) );
    strcpy ( dateStr, asctime(&utc) );
    dateStr[ strlen(dateStr) - 1 ] = 0;
    sprintf (line, "%%%% Start GPS time tStart = %12.3f    (%s GMT)\n", 
	     GPS2REAL8(cfg->startTime), dateStr);
    strcat ( summary, line );
    sprintf (line, "%%%% Total time spanned    = %12.3f s  (%.1f hours)\n", 
	     cfg->duration, cfg->duration/3600 );
    strcat ( summary, line );
    sprintf (line, "%%%% Effective spin-range at tStart: " );
    strcat ( summary, line );
    numSpins = cfg->spinRangeStart->fkdot->length;
    strcat (summary, "fkdot = [ " );
    for (i=0; i < numSpins; i ++ ) 
      {
	sprintf (line, "%.16g:%.16g%s", 
		 cfg->spinRangeStart->fkdot->data[i],
		 cfg->spinRangeStart->fkdot->data[i] + cfg->spinRangeStart->fkdotBand->data[i],
		 (i < numSpins - 1)?", ":" ]\n");
	strcat ( summary, line );
      }

    if ( (cfg->logstring = LALCalloc(1, strlen(summary) + 1 )) == NULL ) {
      ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM);
    }
    strcpy ( cfg->logstring, summary );
  } /* write dataSummary string */

  LogPrintfVerbatim( LOG_DEBUG, cfg->logstring );


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

    fprintf (fplog, "%%%% LOG-FILE of ComputeFStatistic run\n\n");
    fprintf (fplog, "%% User-input:\n");
    fprintf (fplog, "%%----------------------------------------------------------------------\n\n");

    fprintf (fplog, logstr);
    LALFree (logstr);

    /* append an ident-string defining the exact CVS-version of the code used */
    fprintf (fplog, "\n\n%% CVS-versions of executable:\n");
    fprintf (fplog, "%% ----------------------------------------------------------------------\n");
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

  XLALDestroyPulsarSpinRange ( cfg->spinRangeStart );
  XLALDestroyPulsarSpinRange ( cfg->spinRangeRef );

  if ( cfg->logstring ) 
    LALFree ( cfg->logstring );

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
  if ( LALUserVarWasSet(&uvar_dFreq) && (uvar_dFreq < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dFreq not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
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
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

    if ( !useGridFile && !(haveSkyRegion || haveAlphaDelta) )
      {
        LALPrintError ("\nNeed sky-region: either use (Alpha,Delta) or skyRegion!\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
    if ( haveSkyRegion && haveAlphaDelta )
      {
        LALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta)"
		       " OR skyRegion!\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
    if ( useGridFile && !haveGridFile )
      {
        LALPrintError ("\nERROR: gridType=FILE, but no skyGridFile specified!\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);  
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


int
XLALwriteCandidate2file ( FILE *fp,  const PulsarCandidate *cand, const Fcomponents *Fstat )
{
  fprintf (fp, "\n");

  fprintf (fp, "refTime  = % 9d;\n", cand->Doppler.refTime.gpsSeconds );   /* forget about ns... */

  fprintf (fp, "\n");

  /* Amplitude parameters with error-estimates */
  fprintf (fp, "h0       = % .6g;\n", cand->Amp.h0 );
  fprintf (fp, "dh0      = % .6g;\n", cand->dAmp.h0 );
  fprintf (fp, "cosi     = % .6g;\n", cand->Amp.cosi );
  fprintf (fp, "dcosi    = % .6g;\n", cand->dAmp.cosi );
  fprintf (fp, "phi0     = % .6g;\n", cand->Amp.phi0 );
  fprintf (fp, "dphi0    = % .6g;\n", cand->dAmp.phi0 );
  fprintf (fp, "psi      = % .6g;\n", cand->Amp.psi );
  fprintf (fp, "dpsi     = % .6g;\n", cand->dAmp.psi );

  fprintf (fp, "\n");

  /* Doppler parameters */
  fprintf (fp, "Alpha    = % .16g;\n", cand->Doppler.Alpha );
  fprintf (fp, "Delta    = % .16g;\n", cand->Doppler.Delta );
  fprintf (fp, "Freq     = % .16g;\n", cand->Doppler.fkdot[0] );
  fprintf (fp, "f1ddot   = % .16g;\n", cand->Doppler.fkdot[1] );
  fprintf (fp, "f2ddot   = % .16g;\n", cand->Doppler.fkdot[2] );
  fprintf (fp, "f3ddot   = % .16g;\n", cand->Doppler.fkdot[3] );

  fprintf (fp, "\n");

  /* Fstat-values */
  fprintf (fp, "Fa       = % .6g  %+.6gi;\n", Fstat->Fa.re, Fstat->Fa.im );
  fprintf (fp, "Fb       = % .6g  %+.6gi;\n", Fstat->Fb.re, Fstat->Fb.im );
  fprintf (fp, "twoF     = % .6g;\n", 2.0 * Fstat->F );


  return XLAL_SUCCESS;

} /* XLALwriteCandidate2file() */

