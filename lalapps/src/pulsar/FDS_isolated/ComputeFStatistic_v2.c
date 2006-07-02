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

/* GSL includes */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>



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
#define MYSIGN(x) ( ((x) < 0) ? (-1.0):(+1.0) )

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
  ComputeFParams CFparams;		    /**< parameters for the computation of Fstat (e.g Dterms, SSB-precision,...) */
  CHAR *dataSummary;                        /**< descriptive string describing the data (e.g. #SFTs, startTime etc. ..*/
} ConfigVariables;

/** Container to hold all relevant parameters of a 'candidate' 
 */
typedef struct {
  LIGOTimeGPS refTime;		/**< SSB reference GPS-time at which spins are defined */
  LIGOTimeGPS startTime;	/**< internal reference time used to compute Fa,Fb */
  REAL8Vector *fkdotRef;	/**< spin-vector {f, f1dot, f2dot, ... } @ refTime */
  SkyPosition skypos;
  Fcomponents Fstat;		/**< Fstat-value Fa,Fb, using internal reference-time 'startTime' */
  REAL8 aPlus, daPlus;		/**< amplitude-parameters with (Cramer-Rao) error estimators */
  REAL8 aCross, daCross;
  REAL8 phi0, dphi0;		/**< signal-phase @ reference-epoch */
  REAL8 psi, dpsi;
  REAL8 h0, dh0;
  REAL8 cosi, dcosi;
} candidate_t;

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
void EstimateSigParams (LALStatus *, candidate_t *cand, const MultiAMCoeffs *multiAMcoef, REAL8 TsftShat);
void Freemem(LALStatus *,  ConfigVariables *cfg);

void WriteFStatLog (LALStatus *, CHAR *argv[]);
void checkUserInputConsistency (LALStatus *);
int outputBeamTS( const CHAR *fname, const AMCoeffs *amcoe, const DetectorStateSeries *detStates );
void InitEphemeris (LALStatus *, EphemerisData *edat, const CHAR *ephemDir, const CHAR *ephemYear, LIGOTimeGPS epoch);
void getUnitWeights ( LALStatus *, MultiNoiseWeights **multiWeights, const MultiSFTVector *multiSFTs );
int XLALwriteCandidate2file ( FILE *fp,  const candidate_t *cand );
int printGSLmatrix4 ( FILE *fp, const CHAR *prefix, const gsl_matrix *gij );

const char *va(const char *format, ...);	/* little var-arg string helper function */

/*---------- empty initializers ---------- */
static const PulsarTimesParamStruc empty_PulsarTimesParamStruc;
static const BarycenterInput empty_BarycenterInput;
static const SFTConstraints empty_SFTConstraints;
static const ComputeFBuffer empty_ComputeFBuffer;
static const LIGOTimeGPS empty_LIGOTimeGPS;
static const candidate_t empty_candidate;
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
  CHAR buf[512];
  REAL8Vector *fkdotStart = NULL;	/* temporary storage for fkdots */
  REAL8Vector *fkdotRef = NULL;
  ComputeFBuffer cfBuffer = empty_ComputeFBuffer;
  CWParamSpacePoint psPoint;		/* parameter-space point to compute Fstat for */
  UINT4 nFreq, nf1dot, nf2dot, nf3dot;	/* number of frequency- and f1dot-bins */
  UINT4 iFreq, if1dot, if2dot, if3dot;  /* counters over freq- and f1dot- bins */
  REAL8 numTemplates, templateCounter;
  REAL8 tickCounter;
  REAL8 clock0;
  Fcomponents Fstat;
  candidate_t loudestCandidate = empty_candidate;

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

  /* prepare container for loudest candidate */
  loudestCandidate.fkdotRef = XLALCreateREAL8Vector ( GV.spinRangeRef->fkdot->length );
  if ( loudestCandidate.fkdotRef == NULL )
    return COMPUTEFSTATISTIC_EMEM;
  loudestCandidate.refTime = GV.refTime;
  loudestCandidate.startTime = GV.startTime;
  
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
      CHAR *logstr = NULL;
      if ( (fpFstat = fopen (uvar_outputFstat, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputFstat);
	  return (COMPUTEFSTATISTIC_ESYS);
	}
      /* log search-footprint at head of output-file */
      LAL_CALL( LALUserVarGetLog (&status, &logstr,  UVAR_LOGFMT_CMDLINE ), &status );
      fprintf(fpFstat, "## %s\n## %s\n",
	      "$Id$",
	      logstr );
      LALFree ( logstr );
      /* append 'dataSummary' */
      fprintf (fpFstat, "%s", GV.dataSummary );
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

  if ( ( fkdotStart = XLALCreateREAL8Vector ( NUM_SPINS ) ) == NULL ) {
    return COMPUTEFSTATISTIC_EMEM;
  }
  if ( ( fkdotRef = XLALCreateREAL8Vector ( NUM_SPINS ) ) == NULL ) {
    return COMPUTEFSTATISTIC_EMEM;
  }
  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
   */
  
  psPoint.refTime = GV.startTime;	/* we compute at startTime, not refTime right now */  
  psPoint.binary = NULL;	/* binary pulsars not implemented yet */
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

  templateCounter = 0.0; 
  tickCounter = uvar_timerCount - 100;	/* do 100 iterations before first progress-report */
  clock0 = (REAL8)clock() / CLOCKS_PER_SEC;

  while (1)
    {
      LAL_CALL (NextDopplerPos( &status, &dopplerpos, &thisScan ), &status);
      if (thisScan.state == STATE_FINISHED) /* scanned all DopplerPositions yet? */
	break;
      
      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
      thisPoint.longitude = dopplerpos.Alpha;
      thisPoint.latitude = dopplerpos.Delta;
      thisPoint.system = COORDINATESYSTEM_EQUATORIAL;
      LAL_CALL (LALNormalizeSkyPosition(&status, &thisPoint, &thisPoint), &status);

      /* set parameter-space point: sky-position */
      psPoint.skypos = thisPoint;

      /*----- loop over spindown values */
      for ( if3dot = 0; if3dot < nf3dot; if3dot ++ )
	{
	  fkdotStart->data[3] = GV.spinRangeStart->fkdot->data[3] + if3dot * uvar_df3dot;
	  
	  for ( if2dot = 0; if2dot < nf2dot; if2dot ++ )
	    {
	      fkdotStart->data[2] = GV.spinRangeStart->fkdot->data[2] + if2dot * uvar_df2dot;

	      for (if1dot = 0; if1dot < nf1dot; if1dot ++)
		{
		  fkdotStart->data[1] = GV.spinRangeStart->fkdot->data[1] + if1dot * thisScan.df1dot;
	  
		  /* Loop over frequencies to be demodulated */
		  for ( iFreq = 0 ; iFreq < nFreq ; iFreq ++ )
		    {
		      fkdotStart->data[0] = GV.spinRangeStart->fkdot->data[0] + iFreq * thisScan.dFreq;

		      /* set parameter-space point: spin-vector fkdot */
		      psPoint.fkdot = fkdotStart;
		      
		      LAL_CALL( ComputeFStat(&status, &Fstat, &psPoint, GV.multiSFTs, GV.multiNoiseWeights, 
					     GV.multiDetStates, &GV.CFparams, &cfBuffer ), &status );

		      templateCounter += 1.0;
		      tickCounter += 1.0;
		      if ( lalDebugLevel && ( tickCounter > uvar_timerCount) )
			{
			  REAL8 diffSec = (REAL8)(clock()) / CLOCKS_PER_SEC - clock0;
			  REAL8 taup = diffSec / templateCounter ;
			  REAL8 timeLeft = (numTemplates - templateCounter) *  taup;
			  tickCounter = 0.0;
			  LogPrintf (LOG_DEBUG, "Progres: %g/%g = %.2f %% done, Estimated time left: %.0f s\n",
				     templateCounter, numTemplates, templateCounter/numTemplates * 100.0, 
				     timeLeft);
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
		      LAL_CALL ( LALExtrapolatePulsarSpins(&status, fkdotRef, GV.refTime, fkdotStart, 
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
		      if ( Fstat.F > loudestCandidate.Fstat.F )
			{
			  UINT4 len = fkdotRef->length * sizeof( fkdotRef->data[0] );
			  memcpy ( loudestCandidate.fkdotRef->data, fkdotRef->data, len );
			  loudestCandidate.skypos = thisPoint;
			  loudestCandidate.Fstat = Fstat;
			}
		      
		    } /* for i < nBins: loop over frequency-bins */
		} /* For GV.spinImax: loop over 1st spindowns */
	    } /* for if2dot < nf2dot */
	} /* for if3dot < nf3dot */
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

  /* do full parameter-estimation for loudest canidate and output into separate file */
  if ( uvar_outputLoudest )
    {
      MultiAMCoeffs *multiAMcoef = NULL;
      REAL8 norm = GV.Tsft * GV.S_hat;
      FILE *fpLoudest;
      
      LAL_CALL ( LALGetMultiAMCoeffs ( &status, &multiAMcoef, GV.multiDetStates, loudestCandidate.skypos ), 
		 &status);
      /* noise-weigh Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiAMCoeffs ( multiAMcoef, GV.multiNoiseWeights ) != XLAL_SUCCESS ) {
	LALPrintError("\nXLALWeighMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	return COMPUTEFSTATC_EXLAL;
      }

      LAL_CALL ( EstimateSigParams(&status, &loudestCandidate, cfBuffer.multiAMcoef, norm),  &status);

      XLALDestroyMultiAMCoeffs ( multiAMcoef );

      if ( (fpLoudest = fopen (uvar_outputLoudest, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputLoudest);
	  return COMPUTEFSTATISTIC_ESYS;
	}
      /* write header with run-info */
      fprintf (fpLoudest, GV.dataSummary );
      if ( XLALwriteCandidate2file ( fpLoudest,  &loudestCandidate ) != XLAL_SUCCESS )
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

  XLALDestroyREAL8Vector ( fkdotStart );
  XLALDestroyREAL8Vector ( fkdotRef );
  XLALDestroyREAL8Vector ( loudestCandidate.fkdotRef );

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
      LALPrintError ("\nSorry, didn't find any matching SFTs with pattern '%s'!\n\n", uvar_DataFiles );
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
    CHAR dateStr[512], line[512], summary[1024];
    UINT4 i, numDet, numSpins;
    numDet = cfg->multiSFTs->length;
    tp = time(NULL);
    sprintf (summary, "## Started search: %s", asctime( gmtime( &tp ) ) );
    strcat (summary, "## Loaded SFTs: [ " );
    for ( i=0; i < numDet; i ++ ) {
      sprintf (line, "%s:%d%s",  cfg->multiSFTs->data[i]->data->name, 
	       cfg->multiSFTs->data[i]->length,
	       (i < numDet - 1)?", ":" ]\n");
      strcat ( summary, line );
    }
    utc = *XLALGPSToUTC( &utc, (INT4)GPS2REAL8(cfg->startTime) );
    strcpy ( dateStr, asctime(&utc) );
    dateStr[ strlen(dateStr) - 1 ] = 0;
    sprintf (line, "## Start GPS time tStart = %12.3f    (%s GMT)\n", 
	     GPS2REAL8(cfg->startTime), dateStr);
    strcat ( summary, line );
    sprintf (line, "## Total time spanned    = %12.3f s  (%.1f hours)\n", 
	     cfg->duration, cfg->duration/3600 );
    strcat ( summary, line );
    sprintf (line, "## Effective spin-range at tStart: " );
    strcat ( summary, line );
    numSpins = cfg->spinRangeStart->fkdot->length;
    strcat (summary, "fkdot = [ " );
    for (i=0; i < numSpins; i ++ ) {
      sprintf (line, "%.16g:%.16g%s", 
	       cfg->spinRangeStart->fkdot->data[i],
	       cfg->spinRangeStart->fkdot->data[i] + cfg->spinRangeStart->fkdotBand->data[i],
	       (i < numSpins - 1)?", ":" ]\n");
      strcat ( summary, line );
    }
    if ( (cfg->dataSummary = LALCalloc(1, strlen(summary) + 1 )) == NULL ) {
      ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM);
    }
    strcpy ( cfg->dataSummary, summary );
  } /* write dataSummary string */

  LogPrintfVerbatim( LOG_DEBUG, cfg->dataSummary );


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

  XLALDestroyPulsarSpinRange ( cfg->spinRangeStart );
  XLALDestroyPulsarSpinRange ( cfg->spinRangeRef );

  if ( cfg->dataSummary ) 
    LALFree ( cfg->dataSummary );

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


/** Parameter-estimation: based on large parts on Yousuke's notes and implemention (in CFSv1)
 */
void
EstimateSigParams (LALStatus *status, 
		   candidate_t *cand, 
		   const MultiAMCoeffs *multiAMcoef, 
		   REAL8 TsftShat)
{
  REAL8 A1h, A2h, A3h, A4h;
  REAL8 Ad, Bd, Cd, Dd;
  REAL8 normAmu;
  REAL8 A1check, A2check, A3check, A4check;

  REAL8 Asq, Da, disc;
  REAL8 aPlus, aCross;
  REAL8 Ap2, Ac2;
  REAL8 beta;
  REAL8 phi0, psi;
  REAL8 b1, b2, b3;
  REAL8 h0, cosi;
  
  REAL8 cosphi0, sinphi0, cos2psi, sin2psi;

  REAL8 tolerance = LAL_REAL4_EPS;

  gsl_vector *x_mu, *A_Mu;
  gsl_matrix *M_Mu_Nu;
  gsl_matrix *J_Mu_nu, *Jh_Mu_nu;
  gsl_permutation *perm, *permh;
  gsl_matrix *tmp, *tmp2;
  int signum;

  INITSTATUS (status, "EstimateSigParams", rcsid);
  ATTATCHSTATUSPTR (status);

  Ad = multiAMcoef->A;
  Bd = multiAMcoef->B;
  Cd = multiAMcoef->C;
  Dd = multiAMcoef->D;

  normAmu = 2.0 / sqrt(TsftShat);	/* generally *very* small!! */

  /* ----- GSL memory allocation ----- */
  if ( ( x_mu = gsl_vector_calloc (4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }
  if ( ( A_Mu = gsl_vector_calloc (4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }
  if ( ( M_Mu_Nu = gsl_matrix_calloc (4, 4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }
  if ( ( J_Mu_nu = gsl_matrix_calloc (4, 4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }
  if ( ( Jh_Mu_nu = gsl_matrix_calloc (4, 4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }
  if ( ( perm = gsl_permutation_calloc ( 4 )) == NULL ) {
    ABORT ( status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }
  if ( ( permh = gsl_permutation_calloc ( 4 )) == NULL ) {
    ABORT ( status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }
  if ( ( tmp = gsl_matrix_calloc (4, 4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }
  if ( ( tmp2 = gsl_matrix_calloc (4, 4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }


  /* ----- fill vector x_mu */
  gsl_vector_set (x_mu, 0,   cand->Fstat.Fa.re );	/* x_1 */
  gsl_vector_set (x_mu, 1,   cand->Fstat.Fb.re ); 	/* x_2 */
  gsl_vector_set (x_mu, 2, - cand->Fstat.Fa.im );	/* x_3 */
  gsl_vector_set (x_mu, 3, - cand->Fstat.Fb.im );	/* x_4 */

  /* ----- fill matrix M^{mu,nu} [symmetric: use UPPER HALF ONLY!!]*/
  gsl_matrix_set (M_Mu_Nu, 0, 0,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 1, 1,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 0, 1, - Cd / Dd );  

  gsl_matrix_set (M_Mu_Nu, 2, 2,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 3, 3,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 2, 3, - Cd / Dd );  

  /* get (un-normalized) MLE's for amplitudes A^mu  = M^{mu,nu} x_nu */

  /* GSL-doc: int gsl_blas_dsymv (CBLAS_UPLO_t Uplo, double alpha, const gsl_matrix * A, 
   *                              const gsl_vector * x, double beta, gsl_vector * y )
   * 
   * compute the matrix-vector product and sum: y = alpha A x + beta y 
   * for the symmetric matrix A. Since the matrix A is symmetric only its 
   * upper half or lower half need to be stored. When Uplo is CblasUpper 
   * then the upper triangle and diagonal of A are used, and when Uplo 
   * is CblasLower then the lower triangle and diagonal of A are used. 
   */
  gsl_blas_dsymv (CblasUpper, 1.0, M_Mu_Nu, x_mu, 0.0, A_Mu);

  A1h = gsl_vector_get ( A_Mu, 0 );
  A2h = gsl_vector_get ( A_Mu, 1 );
  A3h = gsl_vector_get ( A_Mu, 2 );
  A4h = gsl_vector_get ( A_Mu, 3 );

  LogPrintf (LOG_NORMAL, "norm= %g; A1 = %g, A2 = %g, A3 = %g, A4 = %g\n", 
	     normAmu, A1h, A2h, A3h, A4h );

  Asq = SQ(A1h) + SQ(A2h) + SQ(A3h) + SQ(A4h);
  Da = A1h * A4h - A2h * A3h;
  disc = sqrt ( SQ(Asq) - 4.0 * SQ(Da) );

  Ap2  = 0.5 * ( Asq + disc );
  aPlus = sqrt(Ap2);		/* not yet normalized */

  Ac2 = 0.5 * ( Asq - disc );
  aCross = sqrt( Ac2 );	
  aCross *= MYSIGN ( Da ); 	/* not yet normalized */

  beta = aCross / aPlus;
  
  b1 =   A4h - beta * A1h;
  b2 =   A3h + beta * A2h;
  b3 = - A1h + beta * A4h ;

  psi  = 0.5 * atan2 ( b1, b2 );
  phi0 =       atan2 ( b2, b3 );

  /* Fix remaining sign-ambiguity by checking sign of reconstructed A1 */
  A1check = aPlus * cos(phi0) * cos(2.0*psi) - aCross * sin(phi0) * sin(2*psi);
  if ( A1check * A1h <  0 )
    phi0 += LAL_PI;

  A1check =   aPlus * cos(phi0) * cos(2*psi) - aCross * sin(phi0) * sin(2*psi);  
  A2check =   aPlus * cos(phi0) * sin(2*psi) + aCross * sin(phi0) * cos(2*psi);  
  A3check = - aPlus * sin(phi0) * cos(2*psi) - aCross * cos(phi0) * sin(2*psi);  
  A4check = - aPlus * sin(phi0) * sin(2*psi) + aCross * cos(phi0) * cos(2*psi);  

  LogPrintf (LOG_NORMAL, "reconstructed:    A1 = %g, A2 = %g, A3 = %g, A4 = %g\n", 
	     A1check, A2check, A3check, A4check );

  if ( ( fabs( (A1check - A1h)/A1h ) > tolerance ) ||
       ( fabs( (A2check - A2h)/A2h ) > tolerance ) ||
       ( fabs( (A3check - A3h)/A3h ) > tolerance ) ||
       ( fabs( (A4check - A4h)/A4h ) > tolerance ) )
    {
      LogPrintf (LOG_CRITICAL, 
		 "Difference between estimated and reconstructed Amu exceeds tolerance of %g\n", 
		 tolerance );
    }




  /* propagate phase from internal reference-time 'startTime' to refTime */
  TRY ( LALExtrapolatePulsarPhase (status->statusPtr, &phi0, cand->fkdotRef, cand->refTime, 
				   phi0, cand->startTime ), status);

  /* use gauge-freedom to fix gauge to [0,pi] for both phi0, psi */
  if ( phi0 < 0 )
    phi0 += LAL_TWOPI;

  if ( phi0 > LAL_PI )
    {
      phi0 -= LAL_PI;
      psi  -= LAL_PI_2;
    }
  if ( psi > LAL_PI  )
    psi -= LAL_PI;
  else if ( psi < 0 )
    psi += LAL_PI;

  /* translate A_{+,x} into {h_0, cosi} */
  h0 = aPlus + sqrt ( disc );  /* not yet normalized ! */
  cosi = aCross / h0;

  /* check numerical consistency of estimated Amu and reconstructed */
  cosphi0 = cos(phi0);
  sinphi0 = sin(phi0);
  cos2psi = cos(2*psi);
  sin2psi = sin(2*psi);

  /* fill candidate-struct with the obtained signal-parameters and error-estimations */
  cand->aPlus  = normAmu * aPlus;
  cand->aCross = normAmu * aCross;
  cand->phi0   = phi0;
  cand->psi    = psi;
  cand->h0     = normAmu * h0;
  cand->cosi   = cosi;

  /* ========== Estimate the errors ========== */
  
  /* ----- compute derivatives \partial A^\mu / \partial B^\nu, where
   * we consider the output-variables B^\nu = (A_+, A_x, phi0, psi) 
   */

  /* ----- A1 =   aPlus * cosphi0 * cos2psi - aCross * sinphi0 * sin2psi; ----- */
  gsl_matrix_set (J_Mu_nu, 0, 0,   cosphi0 * cos2psi );	/* dA1/daPlus */
  gsl_matrix_set (J_Mu_nu, 0, 1, - sinphi0 * sin2psi );	/* dA1/daCross */
  gsl_matrix_set (J_Mu_nu, 0, 2,   A3h );		/* dA1/dphi0 */
  gsl_matrix_set (J_Mu_nu, 0, 3, - 2.0 * A2h );		/* dA1/dpsi */
  
  /* ----- A2 =   aPlus * cosphi0 * sin2psi + aCross * sinphi0 * cos2psi; ----- */ 
  gsl_matrix_set (J_Mu_nu, 1, 0,   cosphi0 * sin2psi );	/* dA2/daPlus */
  gsl_matrix_set (J_Mu_nu, 1, 1,   sinphi0 * cos2psi );	/* dA2/daCross */
  gsl_matrix_set (J_Mu_nu, 1, 2,   A4h );		/* dA2/dphi0 */
  gsl_matrix_set (J_Mu_nu, 1, 3,   2.0 * A1h );		/* dA2/dpsi */

  /* ----- A3 = - aPlus * sinphi0 * cos2psi - aCross * cosphi0 * sin2psi; ----- */ 
  gsl_matrix_set (J_Mu_nu, 2, 0, - sinphi0 * cos2psi );	/* dA3/daPlus */ 
  gsl_matrix_set (J_Mu_nu, 2, 1, - cosphi0 * sin2psi );	/* dA3/daCross */
  gsl_matrix_set (J_Mu_nu, 2, 2, - A1h );		/* dA3/dphi0 */  
  gsl_matrix_set (J_Mu_nu, 2, 3, - 2.0 * A4h );		/* dA3/dpsi */

  /* ----- A4 = - aPlus * sinphi0 * sin2psi + aCross * cosphi0 * cos2psi; ----- */ 
  gsl_matrix_set (J_Mu_nu, 3, 0, - sinphi0 * sin2psi );	/* dA4/daPlus */
  gsl_matrix_set (J_Mu_nu, 3, 1,   cosphi0 * cos2psi );	/* dA4/daCross */
  gsl_matrix_set (J_Mu_nu, 3, 2, - A2h );		/* dA4/dphi0 */
  gsl_matrix_set (J_Mu_nu, 3, 3,   2.0 * A3h );		/* dA4/dpsi */

  /* ----- compute derivatives \partial A^\mu / \partial Bh^\nu, where
   * we consider the output-variables Bh^\nu = (h0, cosi, phi0, psi) 
   * where aPlus = 0.5 * h0 * (1 + cosi^2)  and aCross = h0 * cosi
   */
  { /* Ahat^mu is just A^mu with the replacements: A_+ --> A_x, and A_x --> h0 */
    REAL8 A1hat =   aCross * cosphi0 * cos2psi - h0 * sinphi0 * sin2psi;  
    REAL8 A2hat =   aCross * cosphi0 * sin2psi + h0 * sinphi0 * cos2psi;  
    REAL8 A3hat = - aCross * sinphi0 * cos2psi - h0 * cosphi0 * sin2psi;  
    REAL8 A4hat = - aCross * sinphi0 * sin2psi + h0 * cosphi0 * cos2psi;  

    /* ----- A1 =   aPlus * cosphi0 * cos2psi - aCross * sinphi0 * sin2psi; ----- */
    gsl_matrix_set (Jh_Mu_nu, 0, 0,   A1h / h0 );	/* dA1/h0 */
    gsl_matrix_set (Jh_Mu_nu, 0, 1,   A1hat ); 		/* dA1/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 0, 2,   A3h );		/* dA1/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 0, 3, - 2.0 * A2h );	/* dA1/dpsi */
  
    /* ----- A2 =   aPlus * cosphi0 * sin2psi + aCross * sinphi0 * cos2psi; ----- */ 
    gsl_matrix_set (Jh_Mu_nu, 1, 0,   A2h / h0 );	/* dA2/h0 */
    gsl_matrix_set (Jh_Mu_nu, 1, 1,   A2hat ); 		/* dA2/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 1, 2,   A4h );		/* dA2/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 1, 3,   2.0 * A1h );	/* dA2/dpsi */

    /* ----- A3 = - aPlus * sinphi0 * cos2psi - aCross * cosphi0 * sin2psi; ----- */ 
    gsl_matrix_set (Jh_Mu_nu, 2, 0,   A3h / h0 );	/* dA3/h0 */ 
    gsl_matrix_set (Jh_Mu_nu, 2, 1,   A3hat ); 		/* dA3/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 2, 2, - A1h );		/* dA3/dphi0 */  
    gsl_matrix_set (Jh_Mu_nu, 2, 3, - 2.0 * A4h );	/* dA3/dpsi */

    /* ----- A4 = - aPlus * sinphi0 * sin2psi + aCross * cosphi0 * cos2psi; ----- */ 
    gsl_matrix_set (Jh_Mu_nu, 3, 0,   A4h / h0 );	/* dA4/h0 */
    gsl_matrix_set (Jh_Mu_nu, 3, 1,   A4hat ); 		/* dA4/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 3, 2, - A2h );		/* dA4/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 3, 3,   2.0 * A3h );	/* dA4/dpsi */
  }

  
  /* ----- compute inverse matrices J^{-1} by LU-decomposition ----- */
  gsl_linalg_LU_decomp (J_Mu_nu,  perm, &signum );
  gsl_linalg_LU_decomp (Jh_Mu_nu, permh, &signum );

  /* inverse matrix */
  gsl_linalg_LU_invert (J_Mu_nu,  perm,  tmp );	
  gsl_matrix_memcpy ( J_Mu_nu, tmp );

  gsl_linalg_LU_invert (Jh_Mu_nu, permh, tmp );
  gsl_matrix_memcpy ( Jh_Mu_nu, tmp );

  /* ----- compute J^-1 . Minv . (J^-1)^T ----- */
  
  /* GSL-doc: gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, 
   *                          const gsl_matrix *A, const gsl_matrix *B, double beta, gsl_matrix *C)
   * These functions compute the matrix-matrix product and sum 
   * C = \alpha op(A) op(B) + \beta C 
   * where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans 
   * and similarly for the parameter TransB.
   */

  /* first tmp = Minv . (J^-1)^T */
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, M_Mu_Nu, J_Mu_nu, 0.0, tmp );
  /* then J^-1 . tmp , store result in tmp2 */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, J_Mu_nu, tmp, 0.0, tmp2 );
  gsl_matrix_memcpy ( J_Mu_nu, tmp2 );

  /* same for Jh^-1 */
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, M_Mu_Nu, Jh_Mu_nu, 0.0, tmp );
  /* then J^-1 . tmp , store result in J_Mu_nu */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Jh_Mu_nu, tmp, 0.0, tmp2 );
  gsl_matrix_memcpy ( Jh_Mu_nu, tmp2 );
  
  /* ===== debug-output resulting matrices ===== */
  if ( lalDebugLevel )
    {
      printGSLmatrix4 ( stdout, "var(dB^mu, dB^nu) = \n", J_Mu_nu );
      printGSLmatrix4 ( stdout, "var(dBh^mu, dBh^nu) = \n", Jh_Mu_nu );
    }
  
  /* read out principal estimation-errors from diagonal elements */
  cand->daPlus  = normAmu * sqrt( gsl_matrix_get (J_Mu_nu, 0, 0 ) );
  cand->daCross = normAmu * sqrt( gsl_matrix_get (J_Mu_nu, 1, 1 ) );
  cand->dphi0   = sqrt( gsl_matrix_get (J_Mu_nu, 2, 2 ) );
  cand->dpsi    = sqrt( gsl_matrix_get (J_Mu_nu, 3, 3 ) );
  cand->dh0     = normAmu * sqrt( gsl_matrix_get (Jh_Mu_nu, 0, 0 ) );
  cand->dcosi   = sqrt( gsl_matrix_get (Jh_Mu_nu, 1, 1 ) );

  LogPrintf (LOG_NORMAL, "aPlus  = %g +- %g\n", cand->aPlus, cand->daPlus );
  LogPrintf (LOG_NORMAL, "aCross = %g +- %g\n", cand->aCross, cand->daCross );
  LogPrintf (LOG_NORMAL, "h0     = %g +- %g\n", cand->h0, cand->dh0 );
  LogPrintf (LOG_NORMAL, "cosi   = %g +- %g\n", cand->cosi, cand->dcosi );
  LogPrintf (LOG_NORMAL, "phi0   = %g +- %g\n", cand->phi0, cand->dphi0 );
  LogPrintf (LOG_NORMAL, "psi    = %g +- %g\n", cand->psi,  cand->dpsi );

  /* ----- free GSL memory ----- */
  gsl_vector_free ( x_mu );
  gsl_vector_free ( A_Mu );
  gsl_matrix_free ( M_Mu_Nu );
  gsl_matrix_free ( J_Mu_nu );
  gsl_matrix_free ( Jh_Mu_nu );
  gsl_permutation_free ( perm );
  gsl_permutation_free ( permh );
  gsl_matrix_free ( tmp );
  gsl_matrix_free ( tmp2 );

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /*EstimateSigParams()*/


int
XLALwriteCandidate2file ( FILE *fp,  const candidate_t *cand )
{
  UINT4 i, s;

  fprintf (fp, "refTime = %9d;\n", cand->refTime.gpsSeconds );   /* forget about ns... */

  fprintf (fp, "Freq = %.16g;\n", cand->fkdotRef->data[0] );
  s = cand->fkdotRef->length;
  for ( i=1; i < s; i ++ )
    fprintf (fp, "f%ddot  = %.16g;\n", i, cand->fkdotRef->data[i] );

  fprintf (fp, "Alpha = %.16g;\n", cand->skypos.longitude );
  fprintf (fp, "Delta = %.16g;\n", cand->skypos.latitude );

  fprintf (fp, "Fa  = %.6g  %+.6gi;\n", cand->Fstat.Fa.re, cand->Fstat.Fa.im );
  fprintf (fp, "Fb  = %.6g  %+.6gi;\n", cand->Fstat.Fb.re, cand->Fstat.Fb.im );
  fprintf (fp, "F   = %.6g;\n", cand->Fstat.F );

  fprintf (fp, "aPlus  = %.6g;\n", cand->aPlus );
  fprintf (fp, "daPlus   = %.6g;\n", cand->daPlus );
  fprintf (fp, "\n");

  fprintf (fp, "aCross = %.6g;\n", cand->aCross );
  fprintf (fp, "daCross  = %.6g;\n", cand->daCross );
  fprintf (fp, "\n");

  fprintf (fp, "phi0   = %.6g;\n", cand->phi0 );
  fprintf (fp, "dphi0    = %.6g;\n", cand->dphi0 );
  fprintf (fp, "\n");

  fprintf (fp, "psi    = %.6g;\n", cand->psi );
  fprintf (fp, "dpsi     = %.6g;\n", cand->dpsi );
  fprintf (fp, "\n");

  fprintf (fp, "h0     = %.6g;\n", cand->h0 );
  fprintf (fp, "dh0      = %.6g;\n", cand->dh0 );
  fprintf (fp, "\n");

  fprintf (fp, "cosiota= %.6g;\n", cand->cosi );
  fprintf (fp, "dcosiota = %.6g;\n", cand->dcosi );

  return XLAL_SUCCESS;

} /* XLALwriteCandidate2file() */


/** Output 4x4 gsl_matrix in octave-format.
 * return -1 on error, 0 if OK.
 */
int
printGSLmatrix4 ( FILE *fp, const CHAR *prefix, const gsl_matrix *gij )
{
  if ( !gij )
    return -1;
  if ( !fp )
    return -1;
  if ( (gij->size1 != 4) || (gij->size2 != 4 ) )
    return -1;
  if ( !prefix ) 
    return -1;

  fprintf (fp, prefix );
  fprintf (fp, " [ %.9f, %.9f, %.9f, %.9f;\n",
	   gsl_matrix_get ( gij, 0, 0 ), gsl_matrix_get ( gij, 0, 1 ), 
	   gsl_matrix_get ( gij, 0, 2 ), gsl_matrix_get ( gij, 0, 3 ) );
  fprintf (fp, "   %.9f, %.9f, %.9f, %.9f;\n",
	   gsl_matrix_get ( gij, 1, 0 ), gsl_matrix_get ( gij, 1, 1 ),
	   gsl_matrix_get ( gij, 1, 2 ), gsl_matrix_get ( gij, 1, 3 ) );
  fprintf (fp, "   %.9f, %.9f, %.9f, %.9f;\n",
	   gsl_matrix_get ( gij, 2, 0 ), gsl_matrix_get ( gij, 2, 1 ),
	   gsl_matrix_get ( gij, 2, 2 ), gsl_matrix_get ( gij, 2, 3 ) );
  fprintf (fp, "   %.9f, %.9f, %.9f, %.9f ];\n",
	   gsl_matrix_get ( gij, 3, 0 ), gsl_matrix_get ( gij, 3, 1 ),
	   gsl_matrix_get ( gij, 3, 2 ), gsl_matrix_get ( gij, 3, 3 ) );
  
  return 0;

} /* printGSLmatrix4() */
