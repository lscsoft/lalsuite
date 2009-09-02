/*
 * Copyright (C) 2008 Pinkesh Patel, Xavier Siemens
 * Copyright (C) 2007 Reinhard Prix, Iraj Gholami , Pinkesh Patel, Xavier Siemens
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
/** \author P.Patel, X.Siemens, R. Prix, I. Gholami, Y. Ioth,M. Papa
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
#include <lal/FrequencySeries.h>
#include <lal/GSLSupport.h>

#include <lal/NormalizeSFTRngMed.h>
#include <lal/ComputeFstat.h>
#include <lal/LALHough.h>

#include <lal/LogPrintf.h>
#include <lal/DopplerFullScan.h>
#include <lal/ComplexFFT.h>
#include <lal/LALBarycenter.h>

#include <lalapps.h>
#include <lal/Window.h>
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

/* local includes */

#include "../HeapToplist.h"

RCSID("$Id: ComputeFStatistic_resamp.c,v 1.46 2009/03/10 08:40:27 ppatel Exp $");
NRCSID(TEMPORARY,"$Blah$");

/*---------- DEFINES ----------*/

#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

#define EPHEM_YEARS  "00-04"	/**< default range: override with --ephemYear */

#define TRUE (1==1)
#define FALSE (1==0)

/*----- SWITCHES -----*/
#define NUM_SPINS 3		/* number of spin-values to consider: {f, fdot, f2dot, ... } */

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

/* convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )
#define SQ(x) ( (x) * (x) )

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define LAL_INT4_MAX 2147483647
UINT4 FactorialLookup[7] = {1,1,2,6,24,120,720,5040};

/*---------- internal types ----------*/

/** What info do we want to store in our toplist? */
typedef struct {
  PulsarDopplerParams doppler;		/**< Doppler params of this 'candidate' */
  Fcomponents  Fstat;			/**< the Fstat-value (plus Fa,Fb) for this candidate */
  CmplxAntennaPatternMatrix Mmunu;		/**< antenna-pattern matrix Mmunu = 0.5* Sinv*Tsft * [ Ad, Cd; Cd; Bd ] */
} FstatCandidate;


/** moving 'Scanline window' of candidates on the scan-line,
 * which is used to find local 1D maxima.
 */
typedef struct
{
  UINT4 length;
  FstatCandidate *window; 		/**< array holding candidates */
  FstatCandidate *center;		/**< pointer to middle candidate in window */
} scanlineWindow_t;

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  REAL8 Alpha;                              /**< sky position alpha in radians */
  REAL8 Delta;                              /**< sky position delta in radians */
  REAL8 Tsft;                               /**< length of one SFT in seconds */
  LIGOTimeGPS refTime;			    /**< reference-time for pulsar-parameters in SBB frame */
  /* -------------------- Resampling -------------------- */
  REAL8 FFTFreqBand;                        /**< treated outside of DopplerScan for resampling-technique */
  /* ---------------------------------------------------- */
  DopplerRegion searchRegion;		    /**< parameter-space region (at *internalRefTime*) to search over */
  DopplerFullScanState *scanState;          /**< current state of the Doppler-scan */
  PulsarDopplerParams stepSizes;	    /**< user-preferences on Doppler-param step-sizes */
  EphemerisData *ephemeris;		    /**< ephemeris data (from LALInitBarycenter()) */
  MultiSFTVector *multiSFTs;		    /**< multi-IFO SFT-vectors */
  MultiDetectorStateSeries *multiDetStates; /**< pos, vel and LMSTs for detector at times t_i */
  MultiNoiseWeights *multiNoiseWeights;	    /**< normalized noise-weights of those SFTs */
  ComputeFParams CFparams;		    /**< parameters for Fstat (e.g Dterms, SSB-prec,...) */
  CHAR *logstring;                          /**< log containing max-info on this search setup */
  toplist_t* FstatToplist;		    /**< sorted 'toplist' of the NumCandidatesToKeep loudest candidates */
  scanlineWindow_t *scanlineWindow;         /**< moving window of candidates on scanline to find local maxima */
} ConfigVariables;

LALUnit empty_Unit;

/*---------- Global variables ----------*/
extern UINT4 vrbflg;		/**< defined in lalapps.c */

/* ----- User-variables: can be set from config-file or command-line */
INT4 uvar_Dterms;
CHAR *uvar_IFO;
BOOLEAN uvar_SignalOnly;
BOOLEAN uvar_UseNoiseWeights;
REAL8 uvar_Freq;
REAL8 uvar_FreqBand;
REAL8 uvar_dFreq;
REAL8 uvar_Alpha;
CHAR* uvar_RA;
REAL8 uvar_dAlpha;
REAL8 uvar_AlphaBand;
REAL8 uvar_Delta;
CHAR* uvar_Dec;
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
BOOLEAN uvar_projectMetric;
REAL8 uvar_metricMismatch;
CHAR *uvar_skyRegion;
CHAR *uvar_DataFiles;
BOOLEAN uvar_help;
CHAR *uvar_outputLogfile;
CHAR *uvar_outputFstat;
CHAR *uvar_outputLoudest;
CHAR *uvar_outputTimeSeries;
BOOLEAN uvar_countTemplates;

INT4 uvar_NumCandidatesToKeep;
INT4 uvar_clusterOnScanline;

CHAR *uvar_outputFstatHist;        /**< output discrete histogram of all Fstatistic values */
REAL8 uvar_FstatHistBin;           /**< width of an Fstatistic histogram bin */ 

CHAR *uvar_gridFile;
REAL8 uvar_dopplermax;
INT4 uvar_RngMedWindow;
REAL8 uvar_refTime;
REAL8 uvar_internalRefTime;
INT4 uvar_SSBprecision;

INT4 uvar_minStartTime;
INT4 uvar_maxEndTime;
CHAR *uvar_workingDir;
REAL8 uvar_timerCount;
INT4 uvar_upsampleSFTs;

REAL8 BaryRefTime;

/** Defining a multi-IFO complex time series. Keeping the Real and Imaginary parts as seperate vectors, since it is necessary to interpolate them seperately */

typedef struct
{
  UINT4 length;                   /**< Number of IFOs */
  REAL8Sequence** Real;           /**< Real part of the time series */
  REAL8Sequence** Imag;           /**< Imaginary part of the time series */
  REAL8Sequence** Times;          /**< Time Stamps for each point */
  REAL8 f_het;                    /**< The heterodyne frequency */
  REAL8 deltaT;                   /**< The spacing between points if only 1 SFT was used */
  LIGOTimeGPS epoch;              /**< StartTime of the Analysis */
  REAL8 Tspan;                    /**< Span of the Analysis */
  REAL8Sequence *Tdata;           /**< Amount of data time in the Analysis, Tspan - time of Gaps */ 
}MultiCOMPLEX8TimeSeries;

/** A container for resampling variables */
typedef struct
{
  REAL8 dF;            /**< Old dF of fstatvector  */
  REAL8 dF_closest;    /**< New dF of fstatvector */
  UINT4 length;        /**< Old length of fstatvector */
  UINT4 new_length;    /**< New length of fstatvector */
}Resamp_Variables;

/* A container for fftw_complex vector data */
typedef struct
{
  UINT4 length;
  fftw_complex* data;
}FFTWCOMPLEXSeries;

typedef struct
{
  UINT4 length;
  FFTWCOMPLEXSeries** data;
}MultiFFTWCOMPLEXSeries;

/** MultiREAL8Sequence is a Vector of REAL8Sequences */
typedef struct
{
  UINT4 length; /**< Number of IFO's */
  REAL8Sequence** data; /**< REAL8Sequences */
}MultiREAL8Sequence;

/* A buffer for resampling */
typedef struct
{
  const MultiDetectorStateSeries *multiDetStates;/**< buffer for each detStates (store pointer) and skypos */
  REAL8 Alpha, Delta;				/**< skyposition of candidate */
  MultiSSBtimes *multiSSB;
  MultiSSBtimes *multiBinary;
  MultiAMCoeffs *multiAMcoef;
  MultiFFTWCOMPLEXSeries *Saved_a,*Saved_b;    
  MultiREAL8Sequence *MultiCorrDetTimes;             /**< This stores the times in the detector frame that correspond to a linear spacing in the barycentric frame */
  REAL8FrequencySeries *fstatVector;       
  REAL8 StartTimeinBaryCenter;
}ReSampBuffer;


/* A contiguity structure required by the preprocessing function in order to store the information pertaining to the contiguity of SFTs and the gaps between them */

typedef struct
{
  UINT4  length;                    /* Number of Contiguous blocks */
  UINT4* NumContinuous;             /* Number of Contiguous SFTs in each block */
  REAL8* Gap;                       /* Gap between two Contiguous blocks in seconds */
  REAL8* StartTime;                 /* StartTime of each block */
  REAL8* dt;                        /* dt of each block */
  UINT4* N;                         /* Number of Points in each Block including gaps*/
  UINT4* StartIndex;                /* StartIndex of each block */
  UINT4* Ndata;                     /* Number of Points corresponding to data */
}Contiguity;


/* ---------- local prototypes ---------- */
/* Resampling prototypes (start) */

LIGOTimeGPS REAL82GPS(REAL8 Time);
void ComputeFStat_resamp (LALStatus *, const PulsarDopplerParams *doppler, const MultiSFTVector *multiSFTs, const MultiNoiseWeights *multiWeights, const MultiDetectorStateSeries *multiDetStates,const ComputeFParams *params, ReSampBuffer *Buffer, MultiCOMPLEX8TimeSeries *TSeries,Resamp_Variables* Vars);
MultiCOMPLEX8TimeSeries* CalcTimeSeries(MultiSFTVector *multiSFTs,FILE *Out,Resamp_Variables* Vars);
MultiCOMPLEX8TimeSeries* XLALCreateMultiCOMPLEX8TimeSeries(UINT4 i);
void XLALDestroyMultiCOMPLEX8TimeSeries(MultiCOMPLEX8TimeSeries* T);
void XLALDestroyReSampBuffer ( ReSampBuffer *cfb);
void XLALDestroyMultiFFTWCOMPLEXSeries(MultiFFTWCOMPLEXSeries *X);
void XLALDestroyFFTWCOMPLEXSeries(FFTWCOMPLEXSeries *X);
void XLALDestroyREAL8Sequence(REAL8Sequence *X);
void XLALDestroyMultiCmplxAMCoeffs(MultiCmplxAMCoeffs *X);
REAL8Sequence* XLALCreateREAL8Sequence(UINT4 length);
MultiFFTWCOMPLEXSeries *XLALCreateMultiFFTWCOMPLEXSeries(UINT4 length);
FFTWCOMPLEXSeries *XLALCreateFFTWCOMPLEXSeries(UINT4 length);
REAL8 magsquare(fftw_complex f);
INT4 CombineSFTs(COMPLEX16Vector *L,SFTVector *sft_vect,REAL8 FMIN,INT4 number,INT4 startindex);
void ApplyWindow(REAL8Window *Win, COMPLEX16Vector *X);
void Reshuffle(COMPLEX16Vector *X);
void PrintL(FFTWCOMPLEXSeries* L,REAL8 min,REAL8 step);
void ApplyHetCorrection(REAL8Sequence *BaryTimes, REAL8Sequence *DetectorTimes,  const REAL8Sequence *Real, const REAL8Sequence *Imag, REAL8Sequence *Times, MultiCOMPLEX8TimeSeries *TSeries, REAL8Sequence* Real_Corrected, REAL8Sequence* Imag_Corrected);
void ApplySpinDowns(const PulsarSpins *SpinDowns, REAL8 dt, const FFTWCOMPLEXSeries *FaIn, const FFTWCOMPLEXSeries *FbIn, REAL8 BaryStartTime,REAL8Sequence *CorrTimes, REAL8 RefTime, FFTWCOMPLEXSeries *FaInSpinCorrected, FFTWCOMPLEXSeries *FbInSpinCorrected);
void ApplyAandB(REAL8Sequence *FineBaryTimes,REAL8Sequence *BaryTimes,REAL8Sequence *a,REAL8Sequence *b,REAL8Sequence *Real,REAL8Sequence *Imag,FFTWCOMPLEXSeries *FaIn, FFTWCOMPLEXSeries *FbIn, REAL8 TSFT);
double sinc(double t);
void retband(REAL8 t0, REAL8 dt, REAL8* t,REAL8* x, REAL8* y,UINT4 n,UINT4 size, UINT4 terms);
REAL8 strob(REAL8* Xdata, REAL8* Ydata, REAL8 X, UINT4 N1);
REAL8Sequence* ResampleSeries(REAL8Sequence *X_Real,REAL8Sequence *X_Imag,REAL8Sequence *Y_Real,REAL8Sequence *Y_Imag,REAL8 dt,REAL8Vector *BaryTimes, REAL8Sequence *DetectorTimes, REAL8Sequence *Times,REAL8 StartTimeinBaryCenter);
void Heterodyne(REAL8 f_het,REAL8 dt,REAL8 StartTime,REAL8Sequence *Real,REAL8Sequence *Imag);
MultiREAL8Sequence* XLALCreateMultiREAL8Sequence(UINT4 length);
void XLALDestroyMultiREAL8Sequence(MultiREAL8Sequence *X);

/* Resampling prototypes (end) */

int main(int argc,char *argv[]);
void initUserVars (LALStatus *);
void InitFStat ( LALStatus *, ConfigVariables *cfg );
void Freemem(LALStatus *,  ConfigVariables *cfg);

void WriteFStatLog (LALStatus *, CHAR *argv[], const CHAR *log_fname);
void checkUserInputConsistency (LALStatus *);
int outputBeamTS( const CHAR *fname, const AMCoeffs *amcoe, const DetectorStateSeries *detStates );
void InitEphemeris (LALStatus *, EphemerisData *edat, const CHAR *ephemDir, const CHAR *ephemYear, LIGOTimeGPS epoch, BOOLEAN isLISA);
void getUnitWeights ( LALStatus *, MultiNoiseWeights **multiWeights, const MultiSFTVector *multiSFTs );

int write_FstatCandidate_to_fp ( FILE *fp, const FstatCandidate *thisFCand );
int write_PulsarCandidate_to_fp ( FILE *fp,  const PulsarCandidate *pulsarParams, const FstatCandidate *Fcand );

int compareFstatCandidates ( const void *candA, const void *candB );
void getLogString ( LALStatus *status, CHAR **logstr, const ConfigVariables *cfg );

const char *va(const char *format, ...);	/* little var-arg string helper function */

/* ---------- scanline window functions ---------- */
scanlineWindow_t *XLALCreateScanlineWindow ( UINT4 windowWings );
void XLALDestroyScanlineWindow ( scanlineWindow_t *scanlineWindow );
int XLALAdvanceScanlineWindow ( const FstatCandidate *nextCand, scanlineWindow_t *scanWindow );
BOOLEAN XLALCenterIsLocalMax ( const scanlineWindow_t *scanWindow );


/*---------- empty initializers ---------- */
static const ConfigVariables empty_ConfigVariables;
static const FstatCandidate empty_FstatCandidate;
static const ReSampBuffer empty_ReSampBuffer;

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
  MultiCOMPLEX8TimeSeries* TSeries;

  FILE *fpFstat = NULL;
  FILE *fpTSeries = NULL;
  ReSampBuffer Buffer = empty_ReSampBuffer;
  REAL8 numTemplates, templateCounter;
  REAL8 tickCounter;
  time_t clock0;
  PulsarDopplerParams dopplerpos = empty_PulsarDopplerParams;		/* current search-parameters */
  FstatCandidate loudestFCand = empty_FstatCandidate, thisFCand = empty_FstatCandidate;
  UINT4 k;
  ConfigVariables GV = empty_ConfigVariables;		/**< global container for various derived configuration settings */
  REAL8FrequencySeries *fstatVector = Buffer.fstatVector;
  Resamp_Variables Vars;
  gsl_vector_int *Fstat_histogram = NULL;
  Buffer.fstatVector = NULL;

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
  if ( uvar_outputLogfile ) {
    LAL_CALL (WriteFStatLog ( &status, argv, uvar_outputLogfile ), &status );
  }

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

      fprintf (fpFstat, "%s", GV.logstring );
    } /* if outputFstat */

  /* start Fstatistic histogram with a single empty bin */
  if (uvar_outputFstatHist) {
    if ((Fstat_histogram = gsl_vector_int_alloc(1)) == NULL) {
      LALPrintError("\nCouldn't allocate 'Fstat_histogram'\n");
      return COMPUTEFSTATISTIC_EMEM;
    }
    gsl_vector_int_set_zero(Fstat_histogram);
  }

  /* if a complete output of the Time Series was requested,
   * we open and prepare the output-file here */
  if (uvar_outputTimeSeries)
    {
      if ( (fpTSeries = fopen (uvar_outputTimeSeries, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputTimeSeries);
	  return (COMPUTEFSTATISTIC_ESYS);
	}
    } /* if outputTimeSeries */


  /* count number of templates */
  numTemplates = XLALNumDopplerTemplates ( GV.scanState );
  if(uvar_countTemplates)
    {
      printf("%%%% Number of templates: %0.0f\n",numTemplates);
      exit(0);
    }

  /*Call the CalcTimeSeries Function Here*/
  LogPrintf (LOG_DEBUG, "Calculating Time Series.\n");
  TSeries = CalcTimeSeries(GV.multiSFTs,fpTSeries,&Vars);

  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
  */
  templateCounter = 0.0;
  tickCounter = 0;
  clock0 = time(NULL);

  LogPrintf (LOG_DEBUG, "Done Calculating Time Series.\n");

  LogPrintf (LOG_DEBUG, "Starting Main Resampling Loop.\n");
  while ((XLALNextDopplerPos( &dopplerpos, GV.scanState ) == 0) )
    {
      /* main function call: compute F-statistic over frequency-band  */ 
      LAL_CALL( ComputeFStat_resamp ( &status, &dopplerpos, GV.multiSFTs, GV.multiNoiseWeights,GV.multiDetStates, &GV.CFparams, &Buffer, TSeries,&Vars), &status );

      fstatVector = Buffer.fstatVector;
      /* Progress meter */
      templateCounter += 1.0;
      if ( lalDebugLevel && ( ++tickCounter > uvar_timerCount) )
	{
	  REAL8 diffSec = time(NULL) - clock0 ;  /* seconds since start of loop*/
	  REAL8 taup = diffSec / templateCounter ;
	  REAL8 timeLeft = (numTemplates - templateCounter) *  taup;
	  tickCounter = 0.0;
	  LogPrintf (LOG_DEBUG, "Progress: %g/%g = %.2f %% done, Estimated time left: %.0f s\n",
		     templateCounter, numTemplates, templateCounter/numTemplates * 100.0, timeLeft);
	}

      for ( k=0; k < fstatVector->data->length; k++)
	{
	  REAL8 thisF = fstatVector->data->data[k];
	  REAL8 thisFreq = fstatVector->f0 + k * fstatVector->deltaF;
	  /* sanity check on the result */
	  if ( !finite ( thisF ) )
	    {
	      LogPrintf(LOG_CRITICAL, "non-finite F = %.16g\n", thisF );
	      LogPrintf (LOG_CRITICAL, "[Alpha,Delta] = [%.16g,%.16g],\nfkdot=[%.16g,%.16g,%.16g,%16.g]\n",
			 dopplerpos.Alpha, dopplerpos.Delta,
			 thisFreq, dopplerpos.fkdot[1], dopplerpos.fkdot[2], dopplerpos.fkdot[3] );
	  return -1;
	    }

	  /* propagate fkdot from internalRefTime back to refTime for outputting results */
	  /* FIXE: only do this for candidates we're going to write out */
	  dopplerpos.fkdot[0] = thisFreq;
	  /*LAL_CALL ( LALExtrapolatePulsarSpins ( &status, dopplerpos.fkdot, GV.refTime, dopplerpos.fkdot, GV.searchRegion.refTime ), &status );*/
	  dopplerpos.refTime = GV.refTime;

	  /* correct normalization in --SignalOnly case:
	   * we didn't normalize data by 1/sqrt(Tsft * 0.5 * Sh) in terms of
	   * the single-sided PSD Sh: the SignalOnly case is characterized by
	   * setting Sh->1, so we need to divide Fa,Fb by sqrt(0.5*Tsft) and F by (0.5*Tsft)
	   */
	  if ( uvar_SignalOnly )
	    {
	      REAL8 norm = 1.0 / sqrt( 0.5 * GV.Tsft );
	      thisF *= norm * norm;
	      thisF += 2;		/* compute E[2F]:= 4 + SNR^2 */
	    }
	  thisFCand.Fstat.F = thisF;
	  thisFCand.doppler = dopplerpos;

	  /* push new value onto scan-line buffer */
	  XLALAdvanceScanlineWindow ( &thisFCand, GV.scanlineWindow );

	  /* two types of threshold: fixed (TwoFThreshold) and dynamic (NumCandidatesToKeep) */
	  if ( XLALCenterIsLocalMax ( GV.scanlineWindow ) 					/* must be 1D local maximum */
	       && (2.0 * GV.scanlineWindow->center->Fstat.F >= uvar_TwoFthreshold) )	/* fixed threshold */
	    {
	      FstatCandidate *writeCand = GV.scanlineWindow->center;

	      /* insert this into toplist if requested */
	      if ( GV.FstatToplist  )			/* dynamic threshold */
		{
		  if ( insert_into_toplist(GV.FstatToplist, (void*)writeCand ) )
		    LogPrintf ( LOG_DETAIL, "Added new candidate into toplist: 2F = %f\n", 2.0 * writeCand->Fstat.F );
		  else
		    LogPrintf ( LOG_DETAIL, "NOT added the candidate into toplist: 2F = %f\n", 2 * writeCand->Fstat.F );
		}
	      else if ( fpFstat ) 				/* no toplist :write out immediately */
		{
		  if ( write_FstatCandidate_to_fp ( fpFstat, writeCand ) != 0 )
		    {
		      LogPrintf (LOG_CRITICAL, "Failed to write candidate to file.\n");
		      return -1;
		    }
		} /* if outputFstat */

	    } /* if 2F > threshold */

	  /* separately keep track of loudest candidate (for --outputLoudest) */
	  if ( thisFCand.Fstat.F > loudestFCand.Fstat.F )
	    loudestFCand = thisFCand;

	  /* add Fstatistic to histogram if needed */
      if (uvar_outputFstatHist) 
	{
	  
	  /* compute bin */
	  const size_t bin = 2.0 * thisFCand.Fstat.F / uvar_FstatHistBin;

	  /* resize histogram vector if needed */
	  if (!Fstat_histogram || bin >= Fstat_histogram->size)
	    if (NULL == (Fstat_histogram = XLALResizeGSLVectorInt(Fstat_histogram, bin + 1, 0))) {
	      LALPrintError("\nCouldn't (re)allocate 'Fstat_histogram'\n");
	      return COMPUTEFSTATISTIC_EMEM;
	    }
	  
	  /* add to bin */
	  gsl_vector_int_set(Fstat_histogram, bin,
			     gsl_vector_int_get(Fstat_histogram, bin) + 1);
	  
	}
	} /* inner loop about frequency-bins from resampling frequ-band */

    } /* while more Doppler positions to scan */

  /* ----- if using toplist: sort and write it out to file now ----- */
  if ( fpFstat && GV.FstatToplist )
    {
      UINT4 el;

      /* sort toplist */
      LogPrintf ( LOG_DEBUG, "Sorting toplist ... ");
      qsort_toplist ( GV.FstatToplist, compareFstatCandidates );
      LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

      for ( el=0; el < GV.FstatToplist->elems; el ++ )
	{
	  const FstatCandidate *candi;
	  if ( ( candi = (const FstatCandidate *) toplist_elem ( GV.FstatToplist, el )) == NULL ) {
	    LogPrintf ( LOG_CRITICAL, "Internal consistency problems with toplist: contains fewer elements than expected!\n");
	    return -1;
	  }
	  if ( write_FstatCandidate_to_fp ( fpFstat, candi ) != 0 )
	    {
	      LogPrintf (LOG_CRITICAL, "Failed to write candidate to file.\n");
	      return -1;
	    }
	} /* for el < elems in toplist */

    } /* if fpFstat && toplist */

  XLALDestroyMultiCOMPLEX8TimeSeries(TSeries);

  if ( fpFstat )
    {
      fprintf (fpFstat, "%%DONE\n");
      fclose (fpFstat);
      fpFstat = NULL;
    }

  /* ----- estimate amplitude-parameters for the loudest canidate and output into separate file ----- */
  if ( uvar_outputLoudest )
    {
      FILE *fpLoudest;
      PulsarCandidate pulsarParams = empty_PulsarCandidate;
      pulsarParams.Doppler = loudestFCand.doppler;

      LAL_CALL(LALEstimatePulsarAmplitudeParams (&status, &pulsarParams, &loudestFCand.Fstat, &GV.searchRegion.refTime, &loudestFCand.Mmunu ),
	       &status );

      if ( (fpLoudest = fopen (uvar_outputLoudest, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputLoudest);
	  return COMPUTEFSTATISTIC_ESYS;
	}

      /* write header with run-info */
      fprintf (fpLoudest, "%s", GV.logstring );

      /* write this 'candidate' to disc */
      if ( write_PulsarCandidate_to_fp ( fpLoudest,  &pulsarParams, &loudestFCand) != XLAL_SUCCESS )
	{
	  LogPrintf(LOG_CRITICAL, "call to write_PulsarCandidate_to_fp() failed!\n");
	  return COMPUTEFSTATISTIC_ESYS;
	}

      fclose (fpLoudest);

      gsl_matrix_free ( pulsarParams.AmpFisherMatrix );

    } /* write loudest candidate to file */

  LogPrintf (LOG_DEBUG, "Search finished.\n");

  /* write out the Fstatistic histogram */
  if (uvar_outputFstatHist) 
    {
      
      size_t i = 0;
      FILE *fpFstatHist = fopen(uvar_outputFstatHist, "wb");
      
      if (fpFstatHist == NULL) {
	LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputFstat);
	return (COMPUTEFSTATISTIC_ESYS);
      }
      fprintf(fpFstatHist, "%s", GV.logstring);
      
      for (i = 0; i < Fstat_histogram->size; ++i)
	fprintf(fpFstatHist, "%0.3g %0.3g %i\n",
		uvar_FstatHistBin * i,
		uvar_FstatHistBin * (i + 1),
		gsl_vector_int_get(Fstat_histogram, i));
      
      fprintf(fpFstatHist, "%%DONE\n");
      fclose(fpFstatHist);
      
    }
  

  /* Free memory */
  LogPrintf (LOG_DEBUG, "Freeing Doppler grid ... ");
  LAL_CALL ( FreeDopplerFullScan(&status, &GV.scanState), &status);
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

  XLALDestroyReSampBuffer ( &Buffer );

  LAL_CALL ( Freemem(&status, &GV), &status);

  if (Fstat_histogram)
    gsl_vector_int_free(Fstat_histogram);
  
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
  uvar_upsampleSFTs = 1;
  uvar_Dterms 	= 16;
  uvar_FreqBand = 0.0;
  uvar_Alpha 	= 0.0;
  uvar_Delta 	= 0.0;
  uvar_AlphaBand = 0;
  uvar_DeltaBand = 0;
  uvar_skyRegion = NULL;
  uvar_RA = NULL;
  uvar_Dec = NULL;

  uvar_ephemYear = LALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar_ephemYear, EPHEM_YEARS);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar_ephemDir = LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar_ephemDir, DEFAULT_EPHEMDIR);

  uvar_SignalOnly = FALSE;
  uvar_UseNoiseWeights = TRUE;

  uvar_f1dot     = 0.0;
  uvar_f1dotBand = 0.0;

  /* default step-sizes for GRID_FLAT */
  uvar_dAlpha 	= 0.001;
  uvar_dDelta 	= 0.001;
  uvar_dFreq 	 = 0.0; 
  uvar_df1dot    = 0.0;
  uvar_df2dot    = 0.0;
  uvar_df3dot    = 0.0;

  uvar_countTemplates = FALSE;

  uvar_TwoFthreshold = 10.0;
  uvar_NumCandidatesToKeep = 0;
  uvar_clusterOnScanline = 0;

  uvar_metricType =  LAL_PMETRIC_NONE;
  uvar_projectMetric = TRUE;
  uvar_gridType = GRID_FLAT;

  uvar_metricMismatch = 0.02;

  uvar_help = FALSE;
  uvar_outputLogfile = NULL;

  uvar_outputFstat = NULL;
  uvar_outputLoudest = NULL;
  
  uvar_outputFstatHist = NULL;
  uvar_FstatHistBin = 0.1;

  uvar_gridFile = NULL;

  uvar_dopplermax =  1.05e-4;
  uvar_RngMedWindow = 50;	/* for running-median */

  uvar_SSBprecision = SSBPREC_RELATIVISTIC;

  uvar_minStartTime = 0;
  uvar_maxEndTime = LAL_INT4_MAX;

  uvar_workingDir = (CHAR*)LALMalloc(512);
  strcpy(uvar_workingDir, ".");

  uvar_timerCount = 1e5;	/* output a timer/progress count every N templates */

  /* ---------- register all user-variables ---------- */

  /* Pinkesh's Additions */
  
  LALregSTRINGUserVar(status,	outputTimeSeries, 0,  UVAR_OPTIONAL, "Output the Time Series File ");

  /* Original Stuff */

  LALregBOOLUserVar(status, 	help, 		'h', UVAR_HELP,     "Print this message");

  LALregREALUserVar(status, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(status,RA, 		 0 , UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in format hh:mm:ss.sss");
  LALregSTRINGUserVar(status,Dec, 		 0 , UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in format dd:mm:ss.sss");
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
  LALregREALUserVar(status,     dFreq,          'r', UVAR_OPTIONAL, "Frequency resolution in Hz [Default: 1/(2T)]");
  LALregREALUserVar(status, 	df1dot, 	'e', UVAR_OPTIONAL, "Stepsize for f1dot [Default: 1/(2T^2)");
  LALregREALUserVar(status, 	df2dot, 	 0 , UVAR_OPTIONAL, "Stepsize for f2dot [Default: 1/(2T^3)");
  LALregREALUserVar(status, 	df3dot, 	 0 , UVAR_OPTIONAL, "Stepsize for f3dot [Default: 1/(2T^4)");

  LALregSTRINGUserVar(status,	skyRegion, 	'R', UVAR_OPTIONAL, "ALTERNATIVE: Specify sky-region by polygon (or use 'allsky')");
  LALregSTRINGUserVar(status,	DataFiles, 	'D', UVAR_REQUIRED, "File-pattern specifying (multi-IFO) input SFT-files");
  LALregSTRINGUserVar(status, 	IFO, 		'I', UVAR_OPTIONAL, "Detector: 'G1', 'L1', 'H1', 'H2' ...(useful for single-IFO v1-SFTs only!)");
  LALregSTRINGUserVar(status,	ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(status,	ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregBOOLUserVar(status, 	SignalOnly, 	'S', UVAR_OPTIONAL, "Signal only flag");
  LALregBOOLUserVar(status, 	UseNoiseWeights,'W', UVAR_OPTIONAL, "Use SFT-specific noise weights");

  LALregREALUserVar(status, 	TwoFthreshold,	'F', UVAR_OPTIONAL, "Set the threshold for selection of 2F");
  LALregINTUserVar(status, 	gridType,	 0 , UVAR_OPTIONAL, "Grid: 0=flat, 1=isotropic, 2=metric, 3=skygrid-file, 6=grid-file, 7=An*lattice");
  LALregINTUserVar(status, 	metricType,	'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserVar(status, 	metricMismatch,	'X', UVAR_OPTIONAL, "Maximal allowed mismatch for metric tiling");
  LALregSTRINGUserVar(status,	outputLogfile,	 0,  UVAR_OPTIONAL, "Name of log-file identifying the code + search performed");
  LALregSTRINGUserVar(status,	gridFile,	 0,  UVAR_OPTIONAL, "Load grid from this file: sky-grid or full-grid depending on --gridType.");
  LALregREALUserVar(status,	refTime,	 0,  UVAR_OPTIONAL, "SSB reference time for pulsar-paramters [Default: startTime]");
  LALregREALUserVar(status, 	dopplermax, 	'q', UVAR_OPTIONAL, "Maximum doppler shift expected");

  LALregSTRINGUserVar(status,	outputFstat,	 0,  UVAR_OPTIONAL, "Output-file for F-statistic field over the parameter-space");
  LALregSTRINGUserVar(status,   outputLoudest,	 0,  UVAR_OPTIONAL, "Loudest F-statistic candidate + estimated MLE amplitudes");

 LALregSTRINGUserVar(status,outputFstatHist, 0,  UVAR_OPTIONAL, "Output-file for a discrete histogram of all Fstatistic values");
 LALregREALUserVar(status,  FstatHistBin,    0,  UVAR_OPTIONAL, "Width of an Fstatistic histogram bin");

  LALregINTUserVar(status,      NumCandidatesToKeep,0, UVAR_OPTIONAL, "Number of Fstat 'candidates' to keep. (0 = All)");
  LALregINTUserVar(status,      clusterOnScanline, 0, UVAR_OPTIONAL, "Neighbors on each side for finding 1D local maxima on scanline");


  LALregINTUserVar ( status, 	minStartTime, 	 0,  UVAR_OPTIONAL, "Earliest SFT-timestamp to include");
  LALregINTUserVar ( status, 	maxEndTime, 	 0,  UVAR_OPTIONAL, "Latest SFT-timestamps to include");

  /* ----- more experimental/expert options ----- */
  LALregINTUserVar (status, 	SSBprecision,	 0,  UVAR_DEVELOPER, "Precision to use for time-transformation to SSB: 0=Newtonian 1=relativistic");
  LALregINTUserVar(status, 	RngMedWindow,	'k', UVAR_DEVELOPER, "Running-Median window size");
  LALregINTUserVar(status,	Dterms,		't', UVAR_DEVELOPER, "Number of terms to keep in Dirichlet kernel sum");

  LALregSTRINGUserVar(status,   workingDir,     'w', UVAR_DEVELOPER, "Directory to use as work directory.");
  LALregREALUserVar(status, 	timerCount, 	 0,  UVAR_DEVELOPER, "N: Output progress/timer info every N templates");
  LALregREALUserVar(status,	internalRefTime, 0,  UVAR_DEVELOPER, "internal reference time to use for Fstat-computation [Default: startTime]");

  LALregINTUserVar(status,	upsampleSFTs,	 0,  UVAR_DEVELOPER, "(integer) Factor to up-sample SFTs by");
  LALregBOOLUserVar(status, 	projectMetric, 	 0,  UVAR_DEVELOPER, "Use projected metric on Freq=const subspact");
  LALregBOOLUserVar(status, 	countTemplates,  0,  UVAR_DEVELOPER, "Count number of templates (if supported) instead of search");

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* initUserVars() */

/** Load Ephemeris from ephemeris data-files  */
void
InitEphemeris (LALStatus * status,
	       EphemerisData *edat,	/**< [out] the ephemeris-data */
	       const CHAR *ephemDir,	/**< directory containing ephems */
	       const CHAR *ephemYear,	/**< which years do we need? */
	       LIGOTimeGPS epoch,	/**< epoch of observation */
	       BOOLEAN isLISA		/**< hack this function for LISA ephemeris */
	       )
{
#define FNAME_LENGTH 1024
  CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_LOOSE};
  INT4 leap;

  INITSTATUS( status, "InitEphemeris", rcsid );
  ATTATCHSTATUSPTR (status);

  ASSERT ( edat, status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL );
  ASSERT ( ephemYear, status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL );

  if ( ephemDir )
    {
      if ( isLISA )
	snprintf(EphemEarth, FNAME_LENGTH, "%s/ephemMLDC.dat", ephemDir);
      else
	snprintf(EphemEarth, FNAME_LENGTH, "%s/earth%s.dat", ephemDir, ephemYear);

      snprintf(EphemSun, FNAME_LENGTH, "%s/sun%s.dat", ephemDir, ephemYear);
    }
  else
    {
      if ( isLISA )
	snprintf(EphemEarth, FNAME_LENGTH, "ephemMLDC.dat");
      else
	snprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", ephemYear);
      snprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  ephemYear);
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
  PulsarSpinRange spinRangeRef = empty_PulsarSpinRange;

  UINT4 numSFTs;
  LIGOTimeGPS startTime, endTime;

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
  startTime = catalog->data[0].header.epoch;
  endTime   = catalog->data[numSFTs-1].header.epoch;
  XLALGPSAdd(&endTime, cfg->Tsft);	/* add on Tsft to last SFT start-time */

  /* ----- get reference-times (from user if given, use startTime otherwise): ----- */
  if ( LALUserVarWasSet(&uvar_refTime)) {
    XLALGPSSetREAL8(&(cfg->refTime), uvar_refTime);
  }
  else
    cfg->refTime = startTime;

  { /* ----- prepare spin-range at refTime (in *canonical format*, ie all Bands >= 0) ----- */
    REAL8 fMin = MYMIN ( uvar_Freq - uvar_FreqBand/2.0 , uvar_Freq + 3.0*uvar_FreqBand/2.0 );
    REAL8 fMax = MYMAX ( uvar_Freq - uvar_FreqBand/2.0 , uvar_Freq + 3.0*uvar_FreqBand/2.0 );

    REAL8 f1dotMin = MYMIN ( uvar_f1dot, uvar_f1dot + uvar_f1dotBand );
    REAL8 f1dotMax = MYMAX ( uvar_f1dot, uvar_f1dot + uvar_f1dotBand );

    REAL8 f2dotMin = MYMIN ( uvar_f2dot, uvar_f2dot + uvar_f2dotBand );
    REAL8 f2dotMax = MYMAX ( uvar_f2dot, uvar_f2dot + uvar_f2dotBand );

    REAL8 f3dotMin = MYMIN ( uvar_f3dot, uvar_f3dot + uvar_f3dotBand );
    REAL8 f3dotMax = MYMAX ( uvar_f3dot, uvar_f3dot + uvar_f3dotBand );

    spinRangeRef.refTime = cfg->refTime;
    spinRangeRef.fkdot[0] = fMin;
    spinRangeRef.fkdot[1] = f1dotMin;
    spinRangeRef.fkdot[2] = f2dotMin;
    spinRangeRef.fkdot[3] = f3dotMin;

    spinRangeRef.fkdotBand[0] = fMax - fMin;
    spinRangeRef.fkdotBand[1] = f1dotMax - f1dotMin;
    spinRangeRef.fkdotBand[2] = f2dotMax - f2dotMin;
    spinRangeRef.fkdotBand[3] = f3dotMax - f3dotMin;
  } /* spin-range at refTime */

  { /* ----- What frequency-band do we need to read from the SFTs?
     * propagate spin-range from refTime to startTime and endTime of observation
     */
    PulsarSpinRange spinRangeStart, spinRangeEnd;	/* temporary only */
    REAL8 fmaxStart, fmaxEnd, fminStart, fminEnd;

    /* compute spin-range at startTime of observation */
    TRY ( LALExtrapolatePulsarSpinRange (status->statusPtr, &spinRangeStart, startTime, &spinRangeRef ), status );
    /* compute spin-range at endTime of these SFTs */
    TRY ( LALExtrapolatePulsarSpinRange (status->statusPtr, &spinRangeEnd, endTime, &spinRangeStart ), status );

    fminStart = spinRangeStart.fkdot[0];
    /* ranges are in canonical format! */
    fmaxStart = fminStart + spinRangeStart.fkdotBand[0];
    fminEnd   = spinRangeEnd.fkdot[0];
    fmaxEnd   = fminEnd + spinRangeEnd.fkdotBand[0];

    /*  get covering frequency-band  */
    fCoverMax = MYMAX ( fmaxStart, fmaxEnd );
    fCoverMin = MYMIN ( fminStart, fminEnd );
  } /* extrapolate spin-range */

  /* Now Calculate the Covering frequencies with refTime = StartTime */
  spinRangeRef.refTime = startTime;
  { /* ----- What frequency-band do we need to read from the SFTs?
     * propagate spin-range from refTime to startTime and endTime of observation
     */
    PulsarSpinRange spinRangeStart, spinRangeEnd;	/* temporary only */
    REAL8 fmaxStart, fmaxEnd, fminStart, fminEnd;

    /* compute spin-range at startTime of observation */
    TRY ( LALExtrapolatePulsarSpinRange (status->statusPtr, &spinRangeStart, startTime, &spinRangeRef ), status );
    /* compute spin-range at endTime of these SFTs */
    TRY ( LALExtrapolatePulsarSpinRange (status->statusPtr, &spinRangeEnd, endTime, &spinRangeStart ), status );

    fminStart = spinRangeStart.fkdot[0];
    /* ranges are in canonical format! */
    fmaxStart = fminStart + spinRangeStart.fkdotBand[0];
    fminEnd   = spinRangeEnd.fkdot[0];
    fmaxEnd   = fminEnd + spinRangeEnd.fkdotBand[0];

  } /* extrapolate spin-range */

  spinRangeRef.refTime = cfg->refTime;

  {/* ----- load the multi-IFO SFT-vectors ----- */
    UINT4 wings = MYMAX(uvar_Dterms, uvar_RngMedWindow/2 +1);	/* extra frequency-bins needed for rngmed, and Dterms */
    REAL8 fMax = (1.0 + uvar_dopplermax) * fCoverMax + wings / cfg->Tsft; /* correct for doppler-shift and wings */
    REAL8 fMin = (1.0 - uvar_dopplermax) * fCoverMin - wings / cfg->Tsft;

    LogPrintf (LOG_DEBUG, "Loading SFTs ... ");
    TRY ( LALLoadMultiSFTs ( status->statusPtr, &(cfg->multiSFTs), catalog, fMin, fMax ), status );
    LogPrintfVerbatim (LOG_DEBUG, "done.\n");
    TRY ( LALDestroySFTCatalog ( status->statusPtr, &catalog ), status );
  }
  { /* ----- load ephemeris-data ----- */
    CHAR *ephemDir;
    BOOLEAN isLISA = FALSE;

    cfg->ephemeris = LALCalloc(1, sizeof(EphemerisData));
    if ( LALUserVarWasSet ( &uvar_ephemDir ) )
      ephemDir = uvar_ephemDir;
    else
      ephemDir = NULL;

    /* hack: if first detector is LISA, we load MLDC-ephemeris instead of 'earth' files */
    if ( cfg->multiSFTs->data[0]->data[0].name[0] == 'Z' )
      isLISA = TRUE;

    TRY( InitEphemeris (status->statusPtr, cfg->ephemeris, ephemDir, uvar_ephemYear, startTime, isLISA ), status);
  }

  /* ----- obtain the (multi-IFO) 'detector-state series' for all SFTs ----- */
  TRY ( LALGetMultiDetectorStates ( status->statusPtr, &(cfg->multiDetStates), cfg->multiSFTs, cfg->ephemeris ), status );

  /* ----- normalize SFTs and calculate noise-weights ----- */
  if ( uvar_SignalOnly )
      cfg->multiNoiseWeights = NULL;   /* noiseWeights == NULL is equivalent to unit noise-weights in ComputeFstat() */
  else
    {
      UINT4 X, alpha;
      MultiPSDVector *rngmed = NULL;
      cfg->multiNoiseWeights = NULL;
      TRY ( LALNormalizeMultiSFTVect (status->statusPtr, &rngmed, cfg->multiSFTs, uvar_RngMedWindow ), status );
      TRY ( LALComputeMultiNoiseWeights  (status->statusPtr, &(cfg->multiNoiseWeights), rngmed, uvar_RngMedWindow, 0 ), status );

      TRY ( LALDestroyMultiPSDVector (status->statusPtr, &rngmed ), status );
      if ( !uvar_UseNoiseWeights )	/* in that case simply set weights to 1.0 */
	for ( X = 0; X < cfg->multiNoiseWeights->length; X ++ )
	  for ( alpha = 0; alpha < cfg->multiNoiseWeights->data[X]->length; alpha ++ )
	    cfg->multiNoiseWeights->data[X]->data[alpha] = 1.0;
    } /* if ! SignalOnly */

  /* ----- upsample SFTs ----- */
  if ( (lalDebugLevel >= 2) && (uvar_upsampleSFTs > 1) )
  {
    UINT4 X, numDet = cfg->multiSFTs->length;
    LogPrintf (LOG_DEBUG, "Writing original SFTs for debugging ... ");
    for (X=0; X < numDet ; X ++ )
      {
	TRY ( LALWriteSFTVector2Dir ( status->statusPtr, cfg->multiSFTs->data[X], "./", "original", "orig"), status );
      }
    LogPrintfVerbatim ( LOG_DEBUG, "done.\n");
  }

  LogPrintf (LOG_DEBUG, "Upsampling SFTs by factor %d ... ", uvar_upsampleSFTs );
  TRY ( upsampleMultiSFTVector ( status->statusPtr, cfg->multiSFTs, uvar_upsampleSFTs, 16 ), status );
  LogPrintfVerbatim (LOG_DEBUG, "done.\n");

  if ( lalDebugLevel >= 2 && (uvar_upsampleSFTs > 1) )
  {
    UINT4 X, numDet = cfg->multiSFTs->length;
    CHAR tag[60];
    sprintf (tag, "upsampled%02d", uvar_upsampleSFTs );
    LogPrintf (LOG_DEBUG, "Writing upsampled SFTs for debugging ... ");
    for (X=0; X < numDet ; X ++ )
      {
	TRY ( LALWriteSFTVector2Dir ( status->statusPtr, cfg->multiSFTs->data[X], "./", tag, tag), status );
      }
    LogPrintfVerbatim ( LOG_DEBUG, "done.\n");
  }
 
  { /* ----- set up Doppler region (at internalRefTime) to scan ----- */
    LIGOTimeGPS internalRefTime = empty_LIGOTimeGPS;
    PulsarSpinRange spinRangeInt = empty_PulsarSpinRange;
    BOOLEAN haveAlphaDelta = (LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta)) || (LALUserVarWasSet(&uvar_RA) && LALUserVarWasSet(&uvar_Dec));

    /* define sky position variables from user input */
    if (LALUserVarWasSet(&uvar_RA)) 
      {
	/* use Matt Pitkins conversion code found in lal/packages/pulsar/src/BinaryPulsarTiming.c */
	cfg->Alpha = LALDegsToRads(uvar_RA, "alpha");
      }
    else cfg->Alpha = uvar_Alpha;
    if (LALUserVarWasSet(&uvar_Dec)) 
      {
	/* use Matt Pitkins conversion code found in lal/packages/pulsar/src/BinaryPulsarTiming.c */
	cfg->Delta = LALDegsToRads(uvar_Dec, "delta");
      }
    else cfg->Delta = uvar_Delta;
    

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
	TRY ( SkySquare2String( status->statusPtr, &(cfg->searchRegion.skyRegionString),
				cfg->Alpha, cfg->Delta,	uvar_AlphaBand, uvar_DeltaBand), status);
      }

    if ( LALUserVarWasSet ( &uvar_internalRefTime ) ) {
      XLALGPSSetREAL8(&(internalRefTime), uvar_internalRefTime);
    }
    else
      internalRefTime = startTime;

    /* spin searchRegion defined by spin-range at *internal* reference-time */
    TRY ( LALExtrapolatePulsarSpinRange (status->statusPtr, &spinRangeInt, internalRefTime, &spinRangeRef ), status );
    cfg->searchRegion.refTime = spinRangeInt.refTime;
    memcpy ( &cfg->searchRegion.fkdot, &spinRangeInt.fkdot, sizeof(spinRangeInt.fkdot) );
    memcpy ( &cfg->searchRegion.fkdotBand, &spinRangeInt.fkdotBand, sizeof(spinRangeInt.fkdotBand) );

    /* special treatment of frequency band: take out of Doppler search-region for resampling technique */
    cfg->FFTFreqBand = cfg->searchRegion.fkdotBand[0];
    cfg->searchRegion.fkdotBand[0] = 0;		/* Doppler region contains no frequency-band */

  } /* get DopplerRegion */

  /* ----- set computational parameters for F-statistic from User-input ----- */
  cfg->CFparams.Dterms = uvar_Dterms;
  cfg->CFparams.SSBprec = uvar_SSBprecision;
  cfg->CFparams.upsampling = 1.0 * uvar_upsampleSFTs;

  /* ----- set fixed grid step-sizes from user-input for GRID_FLAT ----- */
  cfg->stepSizes.Alpha = uvar_dAlpha;
  cfg->stepSizes.Delta = uvar_dDelta;
  cfg->stepSizes.fkdot[0] = 0;  	/* set default stepsize to FFT spacing: 1/Tspan */
  cfg->stepSizes.fkdot[1] = uvar_df1dot;
  cfg->stepSizes.fkdot[2] = uvar_df2dot;
  cfg->stepSizes.fkdot[3] = uvar_df3dot;
  cfg->stepSizes.orbit = NULL;

  /* ----- set up toplist if requested ----- */
  if ( uvar_NumCandidatesToKeep > 0 )
    if ( create_toplist( &(cfg->FstatToplist), uvar_NumCandidatesToKeep, sizeof(FstatCandidate), compareFstatCandidates) != 0 ) {
      ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
    }

  /* ----- set up scanline-window if requested for 1D local-maximum clustering on scanline ----- */
  if ( (cfg->scanlineWindow = XLALCreateScanlineWindow ( uvar_clusterOnScanline )) == NULL ) {
    ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }

  /* initialize full multi-dimensional Doppler-scanner */
  {
    DopplerFullScanInit scanInit;			/* init-structure for DopperScanner */

    scanInit.searchRegion = cfg->searchRegion;
    scanInit.gridType = uvar_gridType;
    scanInit.gridFile = uvar_gridFile;
    scanInit.metricType = uvar_metricType;
    scanInit.projectMetric = uvar_projectMetric;
    scanInit.metricMismatch = uvar_metricMismatch;
    scanInit.stepSizes = cfg->stepSizes;
    scanInit.ephemeris = cfg->ephemeris;		/* used by Ephemeris-based metric */
    scanInit.startTime = cfg->multiDetStates->startTime;
    scanInit.Tspan     = cfg->multiDetStates->Tspan;
    scanInit.Detector  = &(cfg->multiDetStates->data[0]->detector);	/* just use first IFO for metric */

    LogPrintf (LOG_DEBUG, "Setting up template grid ... ");
    TRY ( InitDopplerFullScan ( status->statusPtr, &cfg->scanState, &scanInit), status);
    LogPrintf (LOG_DEBUG, "template grid ready: %.0f templates.\n", XLALNumDopplerTemplates ( cfg->scanState ) );
  }

  /* ----- produce a log-string describing the data-specific setup ----- */
  TRY ( getLogString ( status->statusPtr, &(cfg->logstring), cfg ), status );
  LogPrintfVerbatim( LOG_DEBUG, cfg->logstring );


  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* InitFStat() */

/** Produce a log-string describing the present run-setup
 */
void
getLogString ( LALStatus *status, CHAR **logstr, const ConfigVariables *cfg )
{
  struct tm utc;
  time_t tp;
  CHAR dateStr[512], line[512], summary[4096];
  CHAR *cmdline = NULL;
  UINT4 i, numDet, numSpins = PULSAR_MAX_SPINS;
  const CHAR *codeID = "$Id: ComputeFStatistic_resamp.c,v 1.46 2009/03/10 08:40:27 ppatel Exp $";
  CHAR *ret = NULL;

  INITSTATUS( status, "getLogString", rcsid );
  ATTATCHSTATUSPTR (status);

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
  utc = *XLALGPSToUTC( &utc, (INT4)GPS2REAL8(cfg->multiDetStates->startTime) );
  strcpy ( dateStr, asctime(&utc) );
  dateStr[ strlen(dateStr) - 1 ] = 0;
  sprintf (line, "%%%% Start GPS time tStart = %12.3f    (%s GMT)\n",
	   GPS2REAL8(cfg->multiDetStates->startTime), dateStr);
  strcat ( summary, line );
  sprintf (line, "%%%% Total time spanned    = %12.3f s  (%.1f hours)\n",
	   cfg->multiDetStates->Tspan, cfg->multiDetStates->Tspan/3600 );
  strcat ( summary, line );
  sprintf (line, "%%%% Pulsar-params refTime = %12.3f \n", GPS2REAL8(cfg->refTime) );
  strcat ( summary, line );
  sprintf (line, "%%%% InternalRefTime       = %12.3f \n", GPS2REAL8(cfg->searchRegion.refTime) );
  strcat ( summary, line );
  sprintf (line, "%%%% Spin-range at internalRefTime: " );
  strcat ( summary, line );

  strcat (summary, "fkdot = [ " );
  for (i=0; i < numSpins; i ++ )
    {
      sprintf (line, "%.16g:%.16g%s",
	       cfg->searchRegion.fkdot[i],
	       cfg->searchRegion.fkdot[i] + cfg->searchRegion.fkdotBand[i],
	       (i < numSpins - 1)?", ":" ]\n");
      strcat ( summary, line );
    }

  if ( (ret = LALCalloc(1, strlen(summary) + 1 )) == NULL ) {
    ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM);
  }

  strcpy ( ret, summary );

  /* return result */
  (*logstr) = ret;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* getLogString() */



/***********************************************************************/
/** Log the all relevant parameters of the present search-run to a log-file.
 * The name of the log-file is log_fname
 * <em>NOTE:</em> Currently this function only logs the user-input and code-versions.
 */
void
WriteFStatLog (LALStatus *status, char *argv[], const CHAR *log_fname )
{
  CHAR *logstr = NULL;
  CHAR command[512] = "";
  FILE *fplog;

  INITSTATUS (status, "WriteFStatLog", rcsid);
  ATTATCHSTATUSPTR (status);

  if ( !log_fname )	/* no logfile given */
    return;

  /* prepare log-file for writing */
  if ( (fplog = fopen(log_fname, "wb" )) == NULL) {
    LogPrintf ( LOG_CRITICAL , "Failed to open log-file '%s' for writing.\n\n", log_fname );
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

  sprintf (command, "ident %s 2> /dev/null | sort -u >> %s", argv[0], log_fname);
  system (command);	/* we don't check this. If it fails, we assume that */
    			/* one of the system-commands was not available, and */
    			/* therefore the CVS-versions will not be logged */

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

  /* destroy FstatToplist if any */
  if ( cfg->FstatToplist )
    free_toplist( &(cfg->FstatToplist) );

  if ( cfg->scanlineWindow )
    XLALDestroyScanlineWindow ( cfg->scanlineWindow );

  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (status->statusPtr), status);

  if ( cfg->searchRegion.skyRegionString )
    LALFree ( cfg->searchRegion.skyRegionString );

  /* Free ephemeris data */
  LALFree(cfg->ephemeris->ephemE);
  LALFree(cfg->ephemeris->ephemS);
  LALFree(cfg->ephemeris);

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
    BOOLEAN haveSkyRegion, haveAlphaDelta, haveGridFile;
    BOOLEAN useSkyGridFile, useFullGridFile, haveMetric, useMetric;

    haveSkyRegion  	= (uvar_skyRegion != NULL);
    haveAlphaDelta 	= (LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta) );
    haveGridFile      	= (uvar_gridFile != NULL);
    useSkyGridFile   	= (uvar_gridType == GRID_FILE_SKYGRID);
    useFullGridFile	= (uvar_gridType == GRID_FILE_FULLGRID);
    haveMetric     	= (uvar_metricType > LAL_PMETRIC_NONE);
    useMetric     	= (uvar_gridType == GRID_METRIC);

    if ( !useFullGridFile && !useSkyGridFile && haveGridFile )
      {
        LALWarning (status, "\nWARNING: gridFile was specified but not needed ... will be ignored\n\n");
      }
    if ( useSkyGridFile && !haveGridFile )
      {
        LALPrintError ("\nERROR: gridType=SKY-FILE, but no --gridFile specified!\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
    if ( useFullGridFile && !haveGridFile )
      {
	LALPrintError ("\nERROR: gridType=GRID-FILE, but no --gridFile specified!\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

    if ( (haveAlphaBand && !haveDeltaBand) || (haveDeltaBand && !haveAlphaBand) )
      {
	LALPrintError ("\nERROR: Need either BOTH (AlphaBand, DeltaBand) or NONE.\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

    if ( haveSkyRegion && haveAlphaDelta )
      {
        LALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta) OR skyRegion!\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
    if ( !useMetric && haveMetric)
      {
        LALWarning (status, "\nWARNING: Metric was specified for non-metric grid... will be ignored!\n");
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

/** write full 'PulsarCandidate' (i.e. Doppler params + Amplitude params + error-bars + Fa,Fb, F, + A,B,C,D
 * RETURN 0 = OK, -1 = ERROR
 */
int
write_PulsarCandidate_to_fp ( FILE *fp,  const PulsarCandidate *pulsarParams, const FstatCandidate *Fcand )
{
  if ( !fp || !pulsarParams || !Fcand  )
    return -1;

  fprintf (fp, "\n");

  fprintf (fp, "refTime  = % 9d;\n", pulsarParams->Doppler.refTime.gpsSeconds );   /* forget about ns... */

  fprintf (fp, "\n");

  /* Amplitude parameters with error-estimates */
  fprintf (fp, "h0       = % .6g;\n", pulsarParams->Amp.h0 );
  fprintf (fp, "dh0      = % .6g;\n", pulsarParams->dAmp.h0 );
  fprintf (fp, "cosi     = % .6g;\n", pulsarParams->Amp.cosi );
  fprintf (fp, "dcosi    = % .6g;\n", pulsarParams->dAmp.cosi );
  fprintf (fp, "phi0     = % .6g;\n", pulsarParams->Amp.phi0 );
  fprintf (fp, "dphi0    = % .6g;\n", pulsarParams->dAmp.phi0 );
  fprintf (fp, "psi      = % .6g;\n", pulsarParams->Amp.psi );
  fprintf (fp, "dpsi     = % .6g;\n", pulsarParams->dAmp.psi );

  fprintf (fp, "\n");

  /* Doppler parameters */
  fprintf (fp, "Alpha    = % .16g;\n", pulsarParams->Doppler.Alpha );
  fprintf (fp, "Delta    = % .16g;\n", pulsarParams->Doppler.Delta );
  fprintf (fp, "Freq     = % .16g;\n", pulsarParams->Doppler.fkdot[0] );
  fprintf (fp, "f1dot    = % .16g;\n", pulsarParams->Doppler.fkdot[1] );
  fprintf (fp, "f2dot    = % .16g;\n", pulsarParams->Doppler.fkdot[2] );
  fprintf (fp, "f3dot    = % .16g;\n", pulsarParams->Doppler.fkdot[3] );

  fprintf (fp, "\n");

  /* Amplitude Modulation Coefficients */
  fprintf (fp, "Ad       = % .6g;\n", Fcand->Mmunu.Ad );
  fprintf (fp, "Bd       = % .6g;\n", Fcand->Mmunu.Bd );
  fprintf (fp, "Cd       = % .6g;\n", Fcand->Mmunu.Cd );
  fprintf (fp, "Sinv_Tsft= % .6g;\n", Fcand->Mmunu.Sinv_Tsft );
  fprintf (fp, "\n");

  /* Fstat-values */
  fprintf (fp, "Fa       = % .6g  %+.6gi;\n", Fcand->Fstat.Fa.re, Fcand->Fstat.Fa.im );
  fprintf (fp, "Fb       = % .6g  %+.6gi;\n", Fcand->Fstat.Fb.re, Fcand->Fstat.Fb.im );
  fprintf (fp, "twoF     = % .6g;\n", 2.0 * Fcand->Fstat.F );

  fprintf (fp, "\nAmpFisher = \\\n" );
  XLALfprintfGSLmatrix ( fp, "%.9g",pulsarParams->AmpFisherMatrix );

  return 0;

} /* write_PulsarCandidate_to_fp() */

/** comparison function for our candidates toplist */
int
compareFstatCandidates ( const void *candA, const void *candB )
{
  if ( ((const FstatCandidate *)candA)->Fstat.F < ((const FstatCandidate *)candB)->Fstat.F )
    return 1;
  else
    return -1;

} /* compareFstatCandidates() */

/** write one 'FstatCandidate' (i.e. only Doppler-params + Fstat) into file 'fp'.
 * Return: 0 = OK, -1 = ERROR
 */
int
write_FstatCandidate_to_fp ( FILE *fp, const FstatCandidate *thisFCand )
{

  if ( !fp || !thisFCand )
    return -1;

  fprintf (fp, "%.16g %.16g %.16g %.6g %.5g %.5g %.9g\n",
	   thisFCand->doppler.fkdot[0], thisFCand->doppler.Alpha, thisFCand->doppler.Delta,
	   thisFCand->doppler.fkdot[1], thisFCand->doppler.fkdot[2], thisFCand->doppler.fkdot[3],
	   2.0 * thisFCand->Fstat.F );

  return 0;

} /* write_candidate_to_fp() */

/* --------------------------------------------------------------------------------
 * Scanline window functions
 * FIXME: should go into a separate file once implementation is settled down ...
 *
 * --------------------------------------------------------------------------------*/

/** Create a scanline window, with given windowWings >= 0.
 * Note: the actual window-size is 1 + 2 * windowWings
 */
scanlineWindow_t *
XLALCreateScanlineWindow ( UINT4 windowWings ) /**< number of neighbors on each side in scanlineWindow */
{
  const CHAR *fn = "XLALCreateScanlineWindow()";
  scanlineWindow_t *ret = NULL;
  UINT4 windowLen = 1 + 2 * windowWings;

  if ( ( ret = LALCalloc ( 1, sizeof(*ret)) ) == NULL ) {
    XLAL_ERROR_NULL( fn, COMPUTEFSTATISTIC_EMEM );
  }

  ret->length = windowLen;

  if ( (ret->window = LALCalloc ( windowLen, sizeof( ret->window[0] ) )) == NULL ) {
    LALFree ( ret );
    XLAL_ERROR_NULL( fn, COMPUTEFSTATISTIC_EMEM );
  }

  ret->center = &(ret->window[ windowWings ]);	/* points to central bin */

  return ret;

} /* XLALCreateScanlineWindow() */

void
XLALDestroyScanlineWindow ( scanlineWindow_t *scanlineWindow )
{
  if ( !scanlineWindow )
    return;

  if ( scanlineWindow->window )
    LALFree ( scanlineWindow->window );

  LALFree ( scanlineWindow );

  return;

} /* XLALDestroyScanlineWindow() */

/** Advance by pushing a new candidate into the scanline-window
 */
int
XLALAdvanceScanlineWindow ( const FstatCandidate *nextCand, scanlineWindow_t *scanWindow )
{
  const CHAR *fn = "XLALAdvanceScanlineWindow()";
  UINT4 i;

  if ( !nextCand || !scanWindow || !scanWindow->window ) {
    XLAL_ERROR (fn, XLAL_EINVAL );
  }

  for ( i=1; i < scanWindow->length; i ++ )
    scanWindow->window[i - 1] = scanWindow->window[i];

  scanWindow->window[ scanWindow->length - 1 ] = *nextCand;	/* struct-copy */

  return XLAL_SUCCESS;

} /* XLALAdvanceScanlineWindow() */

/** check wether central candidate in Scanline-window is a local maximum
 */
BOOLEAN
XLALCenterIsLocalMax ( const scanlineWindow_t *scanWindow )
{
  UINT4 i;
  REAL8 F0;

  if ( !scanWindow || !scanWindow->center )
    return FALSE;

  F0 = scanWindow->center->Fstat.F;

  for ( i=0; i < scanWindow->length; i ++ )
    if ( scanWindow->window[i].Fstat.F > F0 )
      return FALSE;

  return TRUE;

} /* XLALCenterIsLocalMax() */


/*******************************************************************/

/* Assign memory to a MultiREAL8TimeSeries */
MultiCOMPLEX8TimeSeries* XLALCreateMultiCOMPLEX8TimeSeries(UINT4 length)
{

  /*static const char func[] = "XLALCreateMultiCOMPLEX8TimeSeries";*/
  MultiCOMPLEX8TimeSeries *new;
  REAL8Sequence *Temp;
  REAL8Sequence **Real;
  REAL8Sequence **Imag;
  REAL8Sequence **Times;
  REAL8Sequence *Tdata;
  new = XLALMalloc(sizeof(*new));
  Real = XLALMalloc(sizeof(Temp)*length);
  Imag = XLALMalloc(sizeof(Temp)*length);
  Times = XLALMalloc(sizeof(Temp)*length);
  Tdata = XLALCreateREAL8Sequence(length);
  new->Tdata = Tdata;
  new->Times = Times;
  new->Real = Real;
  new->Imag = Imag;
  new->length = length;
  new->f_het = 0;
  new->deltaT = 0;
  new->epoch = REAL82GPS(0);
  return new;
}

MultiREAL8Sequence* XLALCreateMultiREAL8Sequence(UINT4 length)
{
  MultiREAL8Sequence *new;
  REAL8Sequence *Temp;
  REAL8Sequence **data;
  new = XLALMalloc(sizeof(*new));
  data = XLALMalloc(sizeof(Temp)*length);
  new->data = data;
  new->length = length;
  return new;
}

MultiFFTWCOMPLEXSeries *XLALCreateMultiFFTWCOMPLEXSeries(UINT4 length)
{
  MultiFFTWCOMPLEXSeries *new;
  FFTWCOMPLEXSeries *Temp;
  FFTWCOMPLEXSeries **Temp2;
  new = XLALMalloc(sizeof(*new));
  Temp2 = XLALMalloc(sizeof(Temp)*length);
  new->data = Temp2;
  new->length = length;
  return new;
}

void XLALDestroyMultiREAL8Sequence(MultiREAL8Sequence *X)
{
  UINT4 i;
  for(i=0;i<X->length;i++)
    XLALDestroyREAL8Sequence(X->data[i]);
  XLALFree(X->data);
  XLALFree(X);
}

void XLALDestroyMultiFFTWCOMPLEXSeries(MultiFFTWCOMPLEXSeries *X)
{
  UINT4 i;
  for(i=0;i<X->length;i++)
    XLALDestroyFFTWCOMPLEXSeries(X->data[i]);
  XLALFree(X->data);
  XLALFree(X);
}
  

/* Destroy a MultiCOMPLEX8TimeSeries */
void XLALDestroyMultiCOMPLEX8TimeSeries(MultiCOMPLEX8TimeSeries *T)
{
  UINT4 i;
  for(i=0;i<T->length;i++)
    {
      XLALDestroyREAL8Sequence(T->Real[i]);
      XLALDestroyREAL8Sequence(T->Imag[i]);
      XLALDestroyREAL8Sequence(T->Times[i]);
    }
  XLALDestroyREAL8Sequence(T->Tdata);
  XLALFree(T->Real);
  XLALFree(T->Imag);
  XLALFree(T->Times);
  XLALFree(T);
}

/* Create an FFTWCOMPLEXSeries */
FFTWCOMPLEXSeries *XLALCreateFFTWCOMPLEXSeries(UINT4 length)
{
  FFTWCOMPLEXSeries *new;
  fftw_complex *data;
  new = XLALMalloc(sizeof(*new));
  new->length = length;
  data = fftw_malloc(sizeof(fftw_complex)*length);
  new->data = data;
  return(new);
}

/* Destroy an FFTWCOMPLEXSeries */
void XLALDestroyFFTWCOMPLEXSeries(FFTWCOMPLEXSeries *X)
{
  if(X)
    fftw_free(X->data);
  XLALFree(X);
}


/** Destruction of a ReSampBuffer *contents*,
 * i.e. the multiSSB and multiAMcoeff, while the
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...) */
void XLALDestroyReSampBuffer ( ReSampBuffer *cfb)
{
  XLALDestroyMultiSSBtimes ( cfb->multiSSB );
  cfb->multiSSB = NULL;
  XLALDestroyMultiSSBtimes ( cfb->multiBinary );
  cfb->multiBinary = NULL;
  XLALDestroyMultiAMCoeffs ( cfb->multiAMcoef );
  cfb->multiAMcoef = NULL;
  XLALDestroyMultiFFTWCOMPLEXSeries(cfb->Saved_a);
  cfb->Saved_a = NULL;
  XLALDestroyMultiFFTWCOMPLEXSeries(cfb->Saved_b);
  cfb->Saved_b = NULL; 
  XLALDestroyMultiREAL8Sequence(cfb->MultiCorrDetTimes);
  cfb->MultiCorrDetTimes = NULL;
  XLALDestroyREAL8FrequencySeries(cfb->fstatVector);
  cfb->fstatVector = NULL;
  return;
} /* XLALDestroyReSampBuffer() */



/* Returns the magnitude square of a complex number */
REAL8 magsquare(fftw_complex f)
{
  return(f[0]*f[0]+f[1]*f[1]);
}

/* Convert REAL8 to GPS */
LIGOTimeGPS REAL82GPS(REAL8 Time)
{
  LIGOTimeGPS GPS;
  GPS.gpsSeconds = floor((UINT4)(Time));
  GPS.gpsNanoSeconds = (UINT4)((Time - GPS.gpsSeconds)*1e9);
  return(GPS);
}


/* CombineSFTs combines a series of contiguous SFTs into a coherent longer time baseline one SFT */
/* Written by Xavier, Modified by Pinkesh (Xavier please document further if you think it is necessary) */

INT4 CombineSFTs(COMPLEX16Vector *L,SFTVector *sft_vect,REAL8 FMIN,INT4 number,INT4 startindex)
{
  REAL8* sinVal;
  REAL8* cosVal;
  UINT4 index;
  INT4 k = 0;
  INT4 res=64;
  REAL8 STimeBaseLine = 0;
  REAL8 LTimeBaseLine = 0;
  REAL8 deltaF = 0;
  INT4  nDeltaF = 0;            /* Number of Frequency Bins per SFT band */
  INT4 alpha,m;                 /* loop indices */
  REAL8	xTemp;	                /* temp variable for phase model */
  INT4 k1;	                /* defining the sum over which is calculated */
  REAL8	x;		        /* local variable for holding x */
  REAL8	realXP, imagXP; 	/* temp variables used in computation of */
  REAL8	realP, imagP;	        /* real and imaginary parts of P, see CVS */
  INT4	sftIndex;	        /* more temp variables */
  REAL8	y;		        /* local variable for holding y */
  REAL8 realQ, imagQ;

  COMPLEX16 llSFT;

  REAL8 f,if0,ifmin;
  
  sinVal=(REAL8 *)XLALMalloc((res+1)*sizeof(REAL8));
  cosVal=(REAL8 *)XLALMalloc((res+1)*sizeof(REAL8)); 
  for (k=0; k<=res; k++)
    {
      sinVal[k]=sin((LAL_TWOPI*k)/res);
      cosVal[k]=cos((LAL_TWOPI*k)/res);
    }
  

  /* Variable redefinitions for code readability */
  deltaF  = sft_vect->data->deltaF;
  nDeltaF = sft_vect->data->data->length;
  STimeBaseLine = 1.0/deltaF;
  LTimeBaseLine = number*STimeBaseLine;

  if0 = floor(FMIN*STimeBaseLine);

  ifmin = floor(FMIN*STimeBaseLine)-uvar_Dterms;

  /* Loop over frequencies to be demodulated */
  /*for(m = -number ; m < (number)*(if1-if0-1) ; m++ )*/
  for(m = -number; m < ((INT4)L->length)-number; m++)
  {
    llSFT.re =0.0;
    llSFT.im =0.0;

    f=if0*deltaF+m*deltaF/number;

    /* Loop over SFTs that contribute to F-stat for a given frequency */
    for(alpha=0;alpha<number;alpha++)
      {
	REAL8 tsin, tcos, tempFreq;
	COMPLEX8 *Xalpha = sft_vect->data[alpha+startindex].data->data;
	xTemp = (REAL8)if0+(REAL8)m/(REAL8)number;
	realXP = 0.0;
	imagXP = 0.0;
	/* find correct index into LUT -- pick closest point */
	tempFreq = xTemp-(INT4)xTemp;
	index=(INT4)(tempFreq*64 + 0.5); /*just like res above */
	      
	{
	  REAL8 d=LAL_TWOPI*(tempFreq-(REAL8)index/64.0);/*just like res above */
	  REAL8 d2=0.5*d*d;
	  REAL8 ts=sinVal[index];
	  REAL8 tc=cosVal[index];
		
	  tsin=ts+d*tc-d2*ts;
	  tcos=tc-d*ts-d2*tc-1.0;
	}

        tempFreq=LAL_TWOPI*(tempFreq+uvar_Dterms-1);
        k1=(INT4)xTemp-uvar_Dterms+1;
        /* Loop over terms in dirichlet Kernel */
        for(k=0;k<2*uvar_Dterms;k++)
	  {
	    COMPLEX8 Xalpha_k;
	    x = tempFreq-LAL_TWOPI*(REAL8)k;
	    realP = tsin/x;
	    imagP = tcos/x;

	    /* If x is small we need correct x->0 limit of Dirichlet kernel */
	    if(fabs(x) < 0.000001) 
	      {
		realP = 1.0;
		imagP = 0.0;
	      }	 
 
	    sftIndex=k1+k-ifmin+1;

	   
	    /* these four lines compute P*xtilde */
	    Xalpha_k = Xalpha[sftIndex];

	    realXP += Xalpha_k.re*realP;
	    realXP -= Xalpha_k.im*imagP;
	    imagXP += Xalpha_k.re*imagP;
	    imagXP += Xalpha_k.im*realP;
	  }
	y = -LAL_TWOPI*alpha*(if0+(REAL8)m/(REAL8)number);

	realQ = cos(y);
	imagQ = sin(y);

	/* implementation of amplitude demodulation */
	{
	  REAL8 realQXP = realXP*realQ-imagXP*imagQ;
	  REAL8 imagQXP = realXP*imagQ+imagXP*realQ;
	  llSFT.re += realQXP;
	  llSFT.im += imagQXP;
	}
      }      

    L->data[m+number].re = llSFT.re; 
    L->data[m+number].im = llSFT.im; 
    
  }
  XLALFree(sinVal);
  XLALFree(cosVal);
  return 0;

}/*CombineSFTs()*/
   
/* Apply Window applies a window to a complex time series */
void ApplyWindow(REAL8Window *Win, COMPLEX16Vector *X)
{
  UINT4 i = 0;

  /* Check if the length of the Window and the Complex Series match */
  if(Win->data->length != X->length)
    {
      fprintf(stderr,"Window length = %d != Vector length = %d\n",Win->data->length,X->length);
      exit(0);
    }
  
  /* Multiply it to both the Real and Imaginary Parts */
  for(i=0;i<Win->data->length;i++)
    {
      X->data[i].re = Win->data->data[i] * X->data[i].re; /* Real */
      X->data[i].im = Win->data->data[i] * X->data[i].im; /* Imag */
    }

}/*ApplyWindow*/



/* Reshuffle, reshuffles the frequency bins in a format which is compatible with fftw3 */
/* In fftw3, the first half of the frequency series encodes 0 to Nyquist frequency and second half goes from -Ny to 0 */

void Reshuffle(COMPLEX16Vector *X)
{
  UINT4 length = X->length;
  UINT4 N = floor((length/2.0) + 0.51);
  UINT4 M = length - N;
  /* book-keeping variables */
  UINT4 k = 0;
  UINT4 i = 0;

  /* Create a Copy */
  COMPLEX8 *Temp;
  Temp = (COMPLEX8*)XLALMalloc(sizeof(COMPLEX8)*length);
  for(i=0;i<length;i++)
    {
      Temp[i].re = X->data[i].re; /* Real */
      Temp[i].im = X->data[i].im; /* Imag */
    }
  
  /* Copy first half */
  for(i=M;i<length;i++)
    {
      X->data[k].re = Temp[i].re;
      X->data[k].im = Temp[i].im;
      k++;
    }

  /* Copy Second half */
  for(i=0;i<M;i++)
    {
      X->data[k].re = Temp[i].re;
      X->data[k].im = Temp[i].im;
      k++;
    }

  XLALFree(Temp);
}/*Reshuffle*/


/* Prints a Complex Series, with minimum and steps specifying the abscissa of the plot (Used for debugging purposes)*/
void PrintL(FFTWCOMPLEXSeries* L,REAL8 min,REAL8 step)
{
  UINT4 i=0;
  REAL8 x,y;
  fprintf(stderr,"%f %f\n",min,step);
  for(i=0;i<L->length;i++)
    {
      x = L->data[i][0];
      y = L->data[i][1];
      printf("%f %6.12f %6.12f %6.12f\n",min+step*i,x,y,sqrt(x*x+y*y));
    
    }
}


/* CalcTimeSeries calculates a heterodyned downsampled time series.
   It heterodynes the middle of the band to zero and downsamples
   appropriately. The resulting time series is complex and is stored
   in the MultiComplex8TimesSeries structure.
*/
MultiCOMPLEX8TimeSeries* CalcTimeSeries(MultiSFTVector *multiSFTs,FILE *Out,Resamp_Variables* Vars)
{
  
  MultiCOMPLEX8TimeSeries* TSeries;
  Contiguity C;
  UINT4 i,j,k,p,m;         /* Counters */
  
  /* This is the starting time of the time series */ 
  REAL8 StartTime = 0;

  /* This is the end time of the time series */
  REAL8 EndTime = 0;
      
  /* Minimum Frequency used to calculate the TimeSeries */
  REAL8 Fmin = 0;
  
  /* Maximum Frequency used to calculate the TimeSeries */
  REAL8 Fmax = 0;
  
  /* The Time Base line of the SFTs */ 
  REAL8 SFTTimeBaseline = 0;

  /* The Frequency spacing of the SFTs */
  REAL8 deltaF = 0; 
  
  /* Closest Frequency to the minimum frequency asked by the user */
  REAL8 UserFmin_Closest = 0;
  
  /* Difference between user asked minimum frequency and the one that is an exact multiple of dF_closest from TSeries->f_het */
  REAL8 UserFmin_Diff = 0;

  /* Allocate memory for the Time Series */
  TSeries = XLALCreateMultiCOMPLEX8TimeSeries(multiSFTs->length);

  /* Loop over IFOs*/
  for(i=0;i<multiSFTs->length;i++)
    {
      SFTVector *SFT_Vect = multiSFTs->data[i]; /* Copy local SFTVect */
      UINT4 NumofSFTs = SFT_Vect->length;       /* Number of SFTs */
      
      /* Time_Baseline = 1/deltaF*/
      deltaF = SFT_Vect->data[0].deltaF;
      SFTTimeBaseline = floor(1.0/deltaF + 0.5);

      /* Set StartTime and EndTime for minimization/maximization respectively. */
      /* Also calculate the Fmin and Fmax */
      if(i == 0)
	{
	  REAL8 f0 = SFT_Vect->data[0].f0;
	  REAL8 DtermsWings = (REAL8)(uvar_Dterms)*deltaF;
	  INT4 lengthofBand = (SFT_Vect->data[0].data->length); 

	  /* The length of the Band used has to be odd because we want to preserve the bin used as the heterodyne frequency band for each continuous chunk in the data. If the data is odd, then the resulting number of bins from patched SFTs will still be odd and the middle bin will be the heterodyne frequency bin always. If the data is even, then the resulting patched SFTs bins will be odd or even depending on the number of SFTs in that chunk. Thus the bin which will be used as DC or the heterodyne frequency will change all the time. In order to avoid this problem, we must ensure that then length is odd to begin with */
	  /* Ensure that the length is odd */
	  if((lengthofBand % 2) == 0)
	    {
	      lengthofBand--;
	    }

	  StartTime = GPS2REAL8(SFT_Vect->data[0].epoch);
	  EndTime = GPS2REAL8(SFT_Vect->data[NumofSFTs-1].epoch)+SFTTimeBaseline;
	  /* Keep everything except for the Dirichlet Terms */
	  Fmin = f0+DtermsWings;
	  Fmax = f0+((lengthofBand-1)*deltaF)-DtermsWings;

	  /* Middle of Band */
	  TSeries->f_het = Fmin + floor(lengthofBand/2.0-uvar_Dterms)*deltaF;
	  /* Store deltaF */
	  TSeries->deltaT = 1.0/(Fmax-Fmin+deltaF);	  
	}
      
      else
	{
	  if(StartTime > GPS2REAL8(SFT_Vect->data[0].epoch))
	    StartTime = GPS2REAL8(SFT_Vect->data[0].epoch);
	  if(EndTime < (GPS2REAL8(SFT_Vect->data[NumofSFTs-1].epoch) + SFTTimeBaseline))
	    EndTime = (GPS2REAL8(SFT_Vect->data[NumofSFTs-1].epoch) + SFTTimeBaseline);
	} 
 
    }/* Loop over IFOs */
  
  /* The nominal frequency spacing */
  Vars->dF = 1.0/(EndTime-StartTime);
  
  /* Closest dF to Requested frequency spacing */
  Vars->dF_closest = Vars->dF;
  
  /* Nominal length of fstatvector */
  Vars->length = floor((EndTime-StartTime)/TSeries->deltaT + 0.5);
  Vars->new_length = Vars->length;
  
  /* Pick lengths and dF's */
  if(uvar_dFreq <= Vars->dF && (uvar_dFreq > 0))
    {
      Vars->new_length = floor(Vars->length/uvar_dFreq*Vars->dF + 0.5);
      Vars->dF_closest = Vars->length*Vars->dF/Vars->new_length;
    }
  
  if(TSeries->f_het > uvar_Freq)
    UserFmin_Closest = TSeries->f_het - floor((TSeries->f_het-uvar_Freq)/Vars->dF_closest+0.5)*Vars->dF_closest;
  else
    UserFmin_Closest = TSeries->f_het + floor((uvar_Freq-TSeries->f_het)/Vars->dF_closest+0.5)*Vars->dF_closest;
  
  UserFmin_Diff = uvar_Freq - UserFmin_Closest;
  TSeries->f_het += UserFmin_Diff;
  
  /* Store the Starting time */
  TSeries->epoch = REAL82GPS(StartTime); 
  TSeries->Tspan = EndTime - StartTime;

  /*Loop over IFOs*/
  for(i=0;i<multiSFTs->length;i++)
    {
      SFTVector *SFT_Vect = multiSFTs->data[i]; /* Copy local  SFTVect */
      BOOLEAN IsFirst = TRUE;                   /* Bookkeeping Variable */
      UINT4 NumCount = 1;                       /* Number of SFTs in each block */
      /* The number of SFTs is stored for repeated use */
      UINT4 NumofSFTs = SFT_Vect->length;

      /* An Index to keep track of how many SFTs have been processesd so far */
      UINT4 StartIndex = 0; 

      UINT4 PointsinTimeSeries;

      /* Another Bookkeeping variable */
      UINT4 NumofBlocks = 0;

      REAL8Sequence *Times;

      REAL8 CurrentTime = 0;

      UINT4 err = 0;

      /* Initialize C, length = 0 to begin with. But we need to assign memory to Gap and NumContinuous. The maximum number of continuous blocks is the total number of SFTs, therefore it is appropriate to assign that much memory */
      C.length = 0;
      C.Gap = (REAL8*)XLALMalloc(sizeof(REAL8)*NumofSFTs);
      C.StartTime = (REAL8*)XLALMalloc(sizeof(REAL8)*NumofSFTs);
      C.NumContinuous = (UINT4*)XLALMalloc(sizeof(UINT4)*NumofSFTs);
      C.StartIndex = (UINT4*)XLALMalloc(sizeof(UINT4)*NumofSFTs);
      C.dt = (REAL8*)XLALMalloc(sizeof(REAL8)*NumofSFTs);
      C.N = (UINT4*)XLALMalloc(sizeof(UINT4)*NumofSFTs);
      C.Ndata = (UINT4*)XLALMalloc(sizeof(UINT4)*NumofSFTs);
      
      /* Loop over all SFTs in this SFTVector */
      for(j=0;j<NumofSFTs;j++)
	{
	  /* Stores difference in times between two consecutive SFTs */
	  REAL8 TimeDiff;             

	  if(IsFirst)
	    {
	      IsFirst = FALSE;        /* No Longer First */
	      NumCount = 1;           /* Minimum SFTs in each block is 1 */
	      C.StartTime[NumofBlocks] = GPS2REAL8(SFT_Vect->data[j].epoch);
	    }
	  else
	    {
	        /* Calculate the difference in start times between this SFT and the one before it, since this one isnt the first */
	      TimeDiff = GPS2REAL8(SFT_Vect->data[j].epoch)-GPS2REAL8(SFT_Vect->data[j-1].epoch);                   

	      /* If true, it means that these two SFTs are next to each other in time and hence add 1 to the Number continuous */
	      if(TimeDiff == SFTTimeBaseline) 
		NumCount++;           
	      
	      /* Now we are dealing with a new block */
	      else                    
		{
		  IsFirst = TRUE;

		  /* Restart Cycle with this SFT being first */
		  j--;      

		  /* Record the Gap between these two blocks */
		  C.Gap[NumofBlocks] = TimeDiff-SFTTimeBaseline;

		  /* Also Record how many SFTs in this block */
		  C.NumContinuous[NumofBlocks] = NumCount;
		  
		  /* One more in this Block */
		  NumofBlocks += 1;
		}
	    }/*Top most else() */
	}/* Loop over SFTs */

      /* Record information for the last block */
      C.NumContinuous[NumofBlocks] = NumCount;
      C.Gap[NumofBlocks] = EndTime - C.StartTime[NumofBlocks] - C.NumContinuous[NumofBlocks]*SFTTimeBaseline;
      C.length = NumofBlocks + 1; 

      for(k=0;k<C.length;k++)
	{
	  UINT4 NforBlock = floor(SFTTimeBaseline/TSeries->deltaT + 0.5)*C.NumContinuous[k] - C.NumContinuous[k] + 1;
	  C.Ndata[k] = NforBlock;
	  C.dt[k] = SFTTimeBaseline*C.NumContinuous[k]/(REAL8)NforBlock;
	  C.N[k] = NforBlock + floor(C.Gap[k]/C.dt[k] + 0.5);
	  if(k == 0)
	    C.StartIndex[k] = floor((C.StartTime[k]-StartTime)/C.dt[k] + 0.5);
	  else
	    C.StartIndex[k] = C.StartIndex[k-1] + C.N[k-1];
	}


      PointsinTimeSeries = C.StartIndex[0];
      for(k=0;k<C.length;k++)
	{
	  PointsinTimeSeries += C.N[k];
	}

      /* Allocate Memory */
      TSeries->Real[i]  = (REAL8Sequence*)XLALCreateREAL8Sequence(PointsinTimeSeries);
      TSeries->Imag[i]  = (REAL8Sequence*)XLALCreateREAL8Sequence(PointsinTimeSeries);
      TSeries->Times[i] = (REAL8Sequence*)XLALCreateREAL8Sequence(PointsinTimeSeries);
      Times = TSeries->Times[i];

     
      /* Initialize Tdata */
      TSeries->Tdata->data[i] = 0;
      
      /* Set the TSeries to Zeros to deal with gaps. */
      for(p=0;p<TSeries->Real[i]->length;p++)
	{
	  TSeries->Real[i]->data[p] = 0;
	  TSeries->Imag[i]->data[p] = 0;
	}      
      
	
      /* Create a time series for each contiguous block and fit it in the right time period. */

      /* Loop over all contiguous blocks */
      for(k=0;k<C.length;k++)
	{
	  /* We need a small time series variable and a small large SFT variable */
	  COMPLEX16Vector *L = NULL;
	  COMPLEX16Vector *SmallT = NULL;
	  COMPLEX16FFTPlan *plan;

	  /* Number of data points in this contiguous block */
	  UINT4 N = C.Ndata[k];

	   /* Since the data in the Frequency and Time domain both have the same length, we can use one Tukey window for it all */
	  REAL8Window *Win;

	  /* The window will use the Dterms as a rise and a fall, So have to appropriately calculate the fraction */

	  /* Beta is a parameter used to create a Tukey window and signifies what fraction of the total length will constitute a rise and a fall */
	  REAL8 Win_beta = (2.0*uvar_Dterms)/(REAL8)(N);

	  if(Win_beta > 0.01)
	    Win_beta = 0.01;

	  /* Create the Window */
	  Win = XLALCreateTukeyREAL8Window(N,Win_beta);

	  /* Assign some memory */
	  L = XLALCreateCOMPLEX16Vector(N);
	  SmallT = XLALCreateCOMPLEX16Vector(N);

	  TSeries->Tdata->data[i] += SFTTimeBaseline*C.NumContinuous[k];

	  /* Call the CombineSFTs function only if C.NumContinuous  > 1 */
	  if(C.NumContinuous[k] > 1)
	    {
	      err =  CombineSFTs(L,SFT_Vect,Fmin,C.NumContinuous[k],StartIndex);
	    }

	  /* Else just assign the lone SFT to L */
	  else
	    {
	      for(p=0;p<N;p++)
		{
		  L->data[p].re = SFT_Vect->data[StartIndex].data->data[p+uvar_Dterms].re;
		  L->data[p].im = SFT_Vect->data[StartIndex].data->data[p+uvar_Dterms].im;
		}
	    }

	  /*fprintf(stderr," N = %d\n",N);
	  fprintf(stderr,"Fmax = %f , Fmin = %f \n",Fmax,Fmin);
	  PrintL(L,Fmin,deltaF/C.NumContinuous[k]);
	  exit(0);*/
	 	  
	  plan = XLALCreateReverseCOMPLEX16FFTPlan(N,0);
	  
	  /* Now apply a window to L */
	  ApplyWindow(Win,L);
	  
	  /* Reshuffle L to be in a format fftw3 uses */
	  Reshuffle(L);
	  
	  /* Execute plan */
	  p = XLALCOMPLEX16VectorFFT(SmallT,L,plan);

	  ApplyWindow(Win,SmallT);

	  /* Correct for Phase Change due to different starting times of each chunk */
	  for(p=0;p<N;p++)
	    {
	      REAL8 cosphis = cos(LAL_TWOPI*(TSeries->f_het)*(C.StartTime[k]-StartTime));
	      REAL8 sinphis = -sin(LAL_TWOPI*(TSeries->f_het)*(C.StartTime[k]-StartTime));
	      REAL8 Realpart = SmallT->data[p].re;
	      REAL8 Imagpart = SmallT->data[p].im;
	      SmallT->data[p].re = Realpart*cosphis - Imagpart*sinphis;
	      SmallT->data[p].im = Realpart*sinphis + Imagpart*cosphis;
	    }
	  
	  /* Add into appropriate chunk */
	  for(p=0;p<N;p++)
	    {
	      TSeries->Real[i]->data[C.StartIndex[k]+p] = SmallT->data[p].re*deltaF/C.NumContinuous[k];
	      TSeries->Imag[i]->data[C.StartIndex[k]+p] = SmallT->data[p].im*deltaF/C.NumContinuous[k];
	      
	    }
 
	  StartIndex += C.NumContinuous[k];
	  XLALDestroyCOMPLEX16Vector(L);
	  XLALDestroyCOMPLEX16Vector(SmallT);
	  XLALDestroyCOMPLEX16FFTPlan(plan);
	  XLALDestroyREAL8Window(Win);

	}/* Loop over Contiguous Blocks (k) */

      /* Create TimeStamps vector Times */
      CurrentTime = C.StartTime[0] - StartTime;
      m = 0;
      if(C.StartTime[0] > StartTime)
	{
	  for(p=0;p<floor((C.StartTime[0]-StartTime)/C.dt[0] + 0.5);p++)
	    Times->data[m++] = p*C.dt[0];
	}
      for(k=0;k<C.length;k++)
	{
	  CurrentTime = C.StartTime[k] - StartTime;
	  for(p=0;p<C.N[k];p++)
	    {
	      Times->data[C.StartIndex[k] + p] = CurrentTime + p*C.dt[k];
	    }
	}

      XLALFree(C.Gap);
      XLALFree(C.NumContinuous);
      XLALFree(C.StartTime);
      XLALFree(C.StartIndex);
      XLALFree(C.dt);
      XLALFree(C.N);
      XLALFree(C.Ndata);
      
      /* Heterodyne to ensure that the minimum Freq is the one requested by the user */
      for(p=0;p<TSeries->Real[i]->length;p++)
	{
	  REAL8 cosphis = cos(LAL_TWOPI*(UserFmin_Diff)*TSeries->Times[i]->data[p]);
	  REAL8 sinphis = -sin(LAL_TWOPI*(UserFmin_Diff)*TSeries->Times[i]->data[p]);
	  REAL8 Realpart = TSeries->Real[i]->data[p];
	  REAL8 Imagpart = TSeries->Imag[i]->data[p];
	  TSeries->Real[i]->data[p] = Realpart*cosphis - Imagpart*sinphis;
	  TSeries->Imag[i]->data[p] = Realpart*sinphis + Imagpart*cosphis;
	}
	  
      if(Out)
	{
	  fprintf(Out,"# Time Series Output \n");
	  fprintf(Out,"$ %f\n",TSeries->f_het);
	  for(p=0;p<TSeries->Real[i]->length;p++)
	    {
	      fprintf(Out,"%f %f %f %f %f\n",TSeries->Times[i]->data[p],TSeries->Real[i]->data[p], TSeries->Imag[i]->data[p],TSeries->Real[i]->data[p]*TSeries->Real[i]->data[p] + TSeries->Imag[i]->data[p]*TSeries->Imag[i]->data[p],atan2(TSeries->Imag[i]->data[p],TSeries->Real[i]->data[p]));
	    }
	}	
    }/*Loop over Multi-IFOs */
  return(TSeries);
}/*CalcTimeSeries()*/


/****** Returns the sinc function -> sin(pi*x)/x/pi ******/
double sinc(double t)
{
  if(t == 0)
    return 1.0;
  else
    return sin(LAL_PI*t)/t/LAL_PI;
}

/* Function to perform Stroboscopic Resampling (Nearest Point)*/
REAL8 strob(REAL8* Xdata, REAL8* Ydata, REAL8 X, UINT4 N1)
{
  REAL8 fraction;
  INT4 index;
  INT4 N = N1;
  REAL8 X0 = Xdata[0];
  REAL8 dX = Xdata[1] - Xdata[0];
  fraction = (X-X0)/dX;
  index = floor(fraction);
  if(fraction-floor(fraction) > 0.5)
    index++;
  if(index < 0)
    return(Ydata[0]);
  if(index > N-1)
    return(Ydata[N-1]);
  return(Ydata[index]);
}
  

/****** This function returns the interpolated valuegiven the two arrays x and t and other parameters like number of terms and starting time and dt ******/
void retband(REAL8 t0, REAL8 dt, REAL8* t,REAL8* x, REAL8* y,UINT4 n,UINT4 size, UINT4 terms)
{
  UINT4 i,j;
  REAL8 f_s = 1/dt;
  REAL8 Value = 0;
  UINT4 lterms = terms;
  UINT4 rterms = terms;
  for(i=0;i<size;i++)
    {
      UINT4 index;
      Value = 0;
      index = (UINT4)((t[i]-t0)/dt);
      if(index < terms)
	lterms = index;
      if(index > (n-terms))
	rterms = (n-index);
      for(j=0;j<lterms;j++)
	Value += sinc(f_s*(t[i]-t0-(index-j)*dt))*x[index-j];
      for(j=1;j<rterms;j++)
	Value += sinc(f_s*(t[i]-t0-(index+j)*dt))*x[index+j];
      y[i] = Value;
    }
}

/* Resamples the Time Series and returns a timestamps vector, which corresponds to detector times linearly sampled in the barycentric frame */
REAL8Sequence* ResampleSeries(REAL8Sequence *X_Real,REAL8Sequence *X_Imag,REAL8Sequence *Y_Real,REAL8Sequence *Y_Imag,REAL8 dt,REAL8Vector *BaryTimes, REAL8Sequence *DetectorTimes, REAL8Sequence *Times,REAL8 StartTimeinBaryCenter)
{
  UINT4 length = Y_Real->length; /* length of data */
  
  UINT4 i;
  REAL8 x,y;
  REAL8 X0;
  REAL8Sequence *CorrespondingDetTimes; /* Timestamps vector to be returned with detector times corresponding to linear sampling in the barycentric frame */
  gsl_interp_accel *accl;
  gsl_spline *splineinter;
  REAL8 MinTimeStamp,MaxTimeStamp;

  /* GSL's memory allocated for BaryTimes->length, which is the lenght of the data set */
  gsl_interp *lininter = gsl_interp_alloc(gsl_interp_linear,BaryTimes->length);

  CorrespondingDetTimes = (REAL8Sequence*)XLALCreateREAL8Sequence(length);
  MinTimeStamp = Times->data[0];
  MaxTimeStamp = Times->data[Times->length-1];

  /* Initialize GSL */
  accl = gsl_interp_accel_alloc();

  /* Xdata = BaryTimes, Ydata = DetectorTimes */
  gsl_interp_init(lininter,BaryTimes->data,DetectorTimes->data,BaryTimes->length);  
  /* Starting point  in Barycentric frame, not exact but close enough*/
  X0 = StartTimeinBaryCenter;

  for(i=0;i<length;i++)
    {
      x = X0 + i*dt; /* Linear Sampling in Barycentric Frame */
      y = gsl_interp_eval(lininter,BaryTimes->data,DetectorTimes->data,x,accl); /* Calculate Detector Time */
      CorrespondingDetTimes->data[i] = y; /* Store it */
    }

  /* Free GSL stuff */
  gsl_interp_free(lininter);
  gsl_interp_accel_free(accl);

  
  /* Reallocate GSL stuff */
  accl = gsl_interp_accel_alloc();
  splineinter = gsl_spline_alloc(gsl_interp_cspline,X_Real->length);

  /* Xdata is TimeStamps and Ydata is Real part of Time Series */
  gsl_spline_init(splineinter,Times->data,X_Real->data,X_Real->length);
  
  for(i=0;i<length;i++)
    {
      x = CorrespondingDetTimes->data[i]; /* New time to calculate at */
      if(x > MinTimeStamp && x < MaxTimeStamp)
	y = gsl_spline_eval(splineinter,x,accl); /* Interpolate and calculate */
      else
	y = 0;
      Y_Real->data[i] = y; /* Store */
    }

  gsl_spline_free(splineinter);
  gsl_interp_accel_free(accl);

  /* Repeat for Imaginary Part */

  accl = gsl_interp_accel_alloc();
  splineinter = gsl_spline_alloc(gsl_interp_cspline,X_Imag->length);
  gsl_spline_init(splineinter,Times->data,X_Imag->data,X_Imag->length);
    
  for(i=0;i<length;i++)
    {
      x = CorrespondingDetTimes->data[i];
      if(x > MinTimeStamp && x < MaxTimeStamp)
	y = gsl_spline_eval(splineinter,x,accl);
      else
	y = 0;
      Y_Imag->data[i] = y;
    }
  
  gsl_spline_free(splineinter);
  gsl_interp_accel_free(accl);

  /* Return TimeStamps Vector */
  return(CorrespondingDetTimes);

}

void ApplySpinDowns(const PulsarSpins *SpinDowns, REAL8 dt, const FFTWCOMPLEXSeries *FaIn, const FFTWCOMPLEXSeries *FbIn, REAL8 BaryStartTime,REAL8Sequence *CorrTimes, REAL8 RefTime, FFTWCOMPLEXSeries *FaInSpinCorrected, FFTWCOMPLEXSeries *FbInSpinCorrected)
{
  UINT4 i;
  UINT4 j;
  REAL8 Phi;
  REAL8 sinphi, cosphi;
  REAL8 Fareal,Faimag;
  REAL8 Fbreal,Fbimag;
  REAL8 DT;
  REAL8 Phi_M;
  for(i=0;i<CorrTimes->length;i++)
    {
      DT = CorrTimes->data[i] - RefTime;
      Phi_M = i*dt + BaryStartTime - CorrTimes->data[i];

      /* Phi is the sum of all terms */
      Phi = 0;

      for(j=1;j<NUM_SPINS;j++)
	{
	  Phi += 2.0*LAL_PI* (*SpinDowns)[j] *pow(DT,j+1)/FactorialLookup[j+1] + 2.0*LAL_PI*Phi_M* (*SpinDowns)[j] * pow(DT,j)/FactorialLookup[j];
	}
      
      sinphi = sin(Phi);
      cosphi = cos(Phi);
      Fareal = FaIn->data[i][0];
      Faimag = FaIn->data[i][1];
      Fbreal = FbIn->data[i][0];
      Fbimag = FbIn->data[i][1];

      FaInSpinCorrected->data[i][0] = Fareal*cosphi + Faimag*sinphi;
      FaInSpinCorrected->data[i][1] = Faimag*cosphi - Fareal*sinphi;
      FbInSpinCorrected->data[i][0] = Fbreal*cosphi + Fbimag*sinphi;
      FbInSpinCorrected->data[i][1] = Fbimag*cosphi - Fbreal*sinphi;
    }
  
}

/* Applies the extra factor termed by PP as heterodyne correction */
void ApplyHetCorrection(REAL8Sequence *BaryTimes, REAL8Sequence *DetectorTimes, const REAL8Sequence *Real, const REAL8Sequence *Imag, REAL8Sequence *Times, MultiCOMPLEX8TimeSeries *TSeries, REAL8Sequence* Real_Corrected, REAL8Sequence* Imag_Corrected)
{
  UINT4 i; /* Counter */
  REAL8 Phi,retemp,imtemp; /* Temporary variables */ 
  REAL8 shift,cosshift,sinshift; 
  REAL8 x;
  REAL8 y;
  REAL8 Het = TSeries->f_het; /* Heterodyne Frequency */

  /* GSL's interpolation accelerator initialized */
  gsl_interp_accel *accl = gsl_interp_accel_alloc();
  /* GSL's memory allocated for DetectorTimes->length, which is the lenght of the data set */
  gsl_interp *lininter = gsl_interp_alloc(gsl_interp_linear,DetectorTimes->length);
  /* Xdata = DetectorTimes, Ydata = BaryTimes */
  gsl_interp_init(lininter,DetectorTimes->data,BaryTimes->data,DetectorTimes->length);
 
  for(i=0;i<Real->length;i++)
    {
      x = Times->data[i]; /* Select Detector Time */
      y = gsl_interp_eval(lininter,DetectorTimes->data,BaryTimes->data,x,accl); /* Calculate Barycenter time by interpolating */
      Phi = y-x; /* Phi_m as termed in JKS */
      shift = -2.0*LAL_PI*Phi*Het; /* Calculated Shift */
      cosshift = cos(shift);
      sinshift = sin(shift);
      retemp = Real->data[i];
      imtemp = Imag->data[i];
      Real_Corrected->data[i] = retemp*cosshift - imtemp*sinshift;
      Imag_Corrected->data[i] = retemp*sinshift + imtemp*cosshift; /* Apply it */
    }
  
  /* Free GSL stuff */
  gsl_interp_free(lininter);
  gsl_interp_accel_free(accl);

}

/* Multiply the Time Series with the Antenna Patterns */
void ApplyAandB(REAL8Sequence *CorrTimes,REAL8Sequence *DetTimes,REAL8Sequence *a,REAL8Sequence *b,REAL8Sequence *Real,REAL8Sequence *Imag,FFTWCOMPLEXSeries *FaIn, FFTWCOMPLEXSeries *FbIn,REAL8 TSFT)
{
  UINT4 i;  /* Counter */
  REAL8 y;  /* Temporary Storage*/
  UINT4 DetTimesIndex = 0; /* Index of the SFT, whose antenna pattern is to be used */

  for(i=0;i<CorrTimes->length;i++)
    {
      /* Check if the time currently is smaller than the SFT being used */
      if(CorrTimes->data[i] > (DetTimes->data[DetTimesIndex] + TSFT/2.0))
	if(DetTimesIndex < (DetTimes->length - 1))
	  DetTimesIndex++; /* If not, increment, i.e. use next SFT */
	
      y = a->data[DetTimesIndex]; /* Store the a */
      FaIn->data[i][0] = Real->data[i] * y; 
      FaIn->data[i][1] = Imag->data[i] * y;
    }
  /* Reset index to 0 */
  DetTimesIndex = 0;
  /* Apply same logic as in a for b */
  for(i=0;i<CorrTimes->length;i++)
    {
      if(CorrTimes->data[i] > (DetTimes->data[DetTimesIndex] + TSFT/2.0))
	if(DetTimesIndex < (DetTimes->length - 1))
	  DetTimesIndex++;
	
      y = b->data[DetTimesIndex];
      FbIn->data[i][0] = Real->data[i] * y;
      FbIn->data[i][1] = Imag->data[i] * y;
    }
}

/* Heterodynes Time Series */
void Heterodyne(REAL8 f_het,REAL8 dt,REAL8 StartTime,REAL8Sequence *Real,REAL8Sequence *Imag)
{
  UINT4 p;
  REAL8 phase,sinphi,cosphi;
  REAL8 RealPart,ImagPart;
  for(p=0;p<Real->length;p++)
    {
      RealPart = Real->data[p];
      ImagPart = Imag->data[p];
      phase = -2*LAL_PI*f_het*(StartTime+p*dt);
      cosphi = cos(phase);
      sinphi = sin(phase);
      Real->data[p] = RealPart*cosphi-ImagPart*sinphi;
      Imag->data[p] = RealPart*sinphi+ImagPart*cosphi;
    }
}

void ComputeFStat_resamp(LALStatus *status, const PulsarDopplerParams *doppler, const MultiSFTVector *multiSFTs, const MultiNoiseWeights *multiWeights, const MultiDetectorStateSeries *multiDetStates,const ComputeFParams *params,ReSampBuffer *Buffer, MultiCOMPLEX8TimeSeries *TSeries,Resamp_Variables* Vars)
{

  UINT4 numDetectors;
  REAL8 SFTTimeBaseline = 0;
  UINT4 i,p;
  MultiSSBtimes *multiSSB = NULL;
  MultiSSBtimes *multiBinary = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  REAL8 Ad, Bd, Cd, Dd_inv;
  SkyPosition skypos;
  MultiFFTWCOMPLEXSeries *Saved_a;
  MultiFFTWCOMPLEXSeries *Saved_b;
  BOOLEAN SAMESKYPOSITION = FALSE;

  /* Starttime of the analysis */
  REAL8 StartTime = GPS2REAL8(TSeries->epoch);

  /* Span of the Analysis */
  REAL8 Tspan = TSeries->Tspan;

  /* End Time of the Analysis */
  REAL8 EndTime = StartTime + Tspan;

   /* Time spacing */
  REAL8 dt = TSeries->deltaT;

  /* This stores the times in the detector frame that correspond to a linear spacing in the barycentric frame */
  MultiREAL8Sequence *MultiCorrDetTimes = NULL; 

  /* Length of the Data set (Delta Freq) * (Time Span) */
  UINT4 length = Vars->length;

  /* Lenght of actual FFT, if nominal dF is used, then it is same as length */
  UINT4 new_length = Vars->new_length;

  /* Store F-Statistic in a temporary variable before reshuffling */
  REAL8Sequence *Fstat_temp; 

  /* Store Fa and Fb for the whole search */
  REAL8Sequence *Fa_Real;
  REAL8Sequence *Fa_Imag;
  REAL8Sequence *Fb_Real;
  REAL8Sequence *Fb_Imag;

  /* The nominal frequency spacing */
  REAL8 dF = Vars->dF;

  /* Closest dF to Requested frequency spacing */
  REAL8 dF_closest = Vars->dF_closest;

  /* The output fstatVector */
  REAL8FrequencySeries *fstatVector = NULL;

  /* Allocate Memory to some common variables */
  Fa_Real = XLALCreateREAL8Sequence(new_length);
  Fb_Real = XLALCreateREAL8Sequence(new_length);
  Fa_Imag = XLALCreateREAL8Sequence(new_length);
  Fb_Imag = XLALCreateREAL8Sequence(new_length);
  Fstat_temp = XLALCreateREAL8Sequence(new_length);

  for(p=0;p<new_length;p++)
    {
      Fa_Real->data[p] = 0;
      Fb_Real->data[p] = 0;
      Fa_Imag->data[p] = 0;
      Fb_Imag->data[p] = 0;
      Fstat_temp->data[p] = 0;
    }
 

  /* Check if fstatVector exists */
  if(Buffer && Buffer->fstatVector)
    {
      fstatVector = Buffer->fstatVector;
    }
  else
    {
      /* prepare Fstat-vector over frequencies to hold output results */
      /* Number of Frequency Bins */
      UINT4 numFreqBins = floor(uvar_FreqBand/dF_closest + 0.5);
      if(Buffer->fstatVector)
	XLALDestroyREAL8FrequencySeries(Buffer->fstatVector);
      fstatVector = XLALCreateREAL8FrequencySeries ("Fstat vector", &doppler->refTime, uvar_Freq, dF_closest, &empty_Unit, numFreqBins );
      Buffer->fstatVector = fstatVector;
    }
      

  
  /* Check if it the previous SkyPosition */
  if ( Buffer
       && ( Buffer->multiDetStates == multiDetStates )
       && ( Buffer->Alpha == doppler->Alpha )
       && ( Buffer->Delta == doppler->Delta )
       )
    SAMESKYPOSITION = TRUE;

  ATTATCHSTATUSPTR (status);

  ASSERT ( multiSFTs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  numDetectors = multiSFTs->length;

  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );

  /* check input */
  ASSERT ( multiSFTs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  numDetectors = multiSFTs->length;
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  if ( multiWeights ) {
    ASSERT ( multiWeights->length == numDetectors , status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  }

  /* check if that skyposition SSB+AMcoef were already buffered */
  if ( Buffer
       && ( Buffer->multiDetStates == multiDetStates )
       && ( Buffer->Alpha == doppler->Alpha )
       && ( Buffer->Delta == doppler->Delta )
       && Buffer->multiSSB )
    { /* yes ==> reuse */
      multiSSB = Buffer->multiSSB;
 

      /* re-use (LWL) AM coefficients whenever available */
      if ( Buffer->multiAMcoef )
	multiAMcoef = Buffer->multiAMcoef;

    } /* if have buffered stuff to reuse */
  else
    {
      skypos.system =   COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;
      TRY ( LALGetMultiSSBtimes ( status->statusPtr, &multiSSB, multiDetStates, skypos, doppler->refTime, params->SSBprec ), status );
      if ( Buffer )
	{
	  XLALDestroyMultiSSBtimes ( Buffer->multiSSB );
	  Buffer->multiSSB = multiSSB;
	  Buffer->Alpha = doppler->Alpha;
	  Buffer->Delta = doppler->Delta;
	  Buffer->multiDetStates = multiDetStates;
	} /* buffer new SSB times */

    } /* could not reuse previously buffered quantites */

    /* new orbital parameter corrections if not already buffered */
  if ( doppler->orbit )
    {
      /* if already buffered */
      if ( Buffer && Buffer->multiBinary )
	{ /* yes ==> reuse */
	  multiBinary = Buffer->multiBinary;
	}
      else
	{
	  /* compute binary time corrections to the SSB time delays and SSB time derivitive */
	  TRY ( LALGetMultiBinarytimes ( status->statusPtr, &multiBinary, multiSSB, multiDetStates, doppler->orbit, doppler->refTime ), status );

	  /* store these in buffer if available */
	  if ( Buffer )
	    {
	      XLALDestroyMultiSSBtimes ( Buffer->multiBinary );
	      Buffer->multiBinary = multiBinary;
	    } /* if Buffer */
 	}
    }
  else multiBinary = multiSSB;

  if ( !multiAMcoef )
    {
      /* compute new AM-coefficients */
      LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, skypos );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      /* noise-weight Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	LALPrintError("\nXLALWeighMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
      }

      /* store these in buffer if available */
      if ( Buffer )
	{
	  XLALDestroyMultiAMCoeffs ( Buffer->multiAMcoef );
	  Buffer->multiAMcoef = multiAMcoef;
	} /* if Buffer */

    } /* if LWL AM coefficient need to be computed */

  if ( multiAMcoef )
    {
      Ad = multiAMcoef->Mmunu.Ad;
      Bd = multiAMcoef->Mmunu.Bd;
      Cd = multiAMcoef->Mmunu.Cd;
      Dd_inv = 1.0 / multiAMcoef->Mmunu.Dd;
    }
  else
    {
      LALPrintError ( "Programming error: 'multiAMcoef' not available!\n");
      ABORT ( status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
    }


  /* Store the SFT Time Baseline */
  SFTTimeBaseline = floor(1.0/multiSFTs->data[0]->data->deltaF + 0.5);
	 
  /* Saved_a and Saved_b are the local copies of F_a and F_b before applying SpinDowns. The Buffered versions are stored in the Buffer */
  if(SAMESKYPOSITION)
    {
      Saved_a = Buffer->Saved_a;
      Saved_b = Buffer->Saved_b;
      MultiCorrDetTimes = Buffer->MultiCorrDetTimes;
    }

  /* If not same sky position, create new ones for the buffer */
  if(!SAMESKYPOSITION)
    {
      if(Buffer->Saved_a)
	XLALDestroyMultiFFTWCOMPLEXSeries(Buffer->Saved_a);
      Saved_a = XLALCreateMultiFFTWCOMPLEXSeries(numDetectors);
      Buffer->Saved_a = Saved_a;

      if(Buffer->Saved_b)
	XLALDestroyMultiFFTWCOMPLEXSeries(Buffer->Saved_b);
      Saved_b = XLALCreateMultiFFTWCOMPLEXSeries(numDetectors);
      Buffer->Saved_b = Saved_b;

      if(Buffer->MultiCorrDetTimes)
	XLALDestroyMultiREAL8Sequence(Buffer->MultiCorrDetTimes);
      MultiCorrDetTimes = XLALCreateMultiREAL8Sequence(numDetectors);
      Buffer->MultiCorrDetTimes = MultiCorrDetTimes;
    }

  /* Store the earliest BaryStartTime */
  if(!SAMESKYPOSITION)
    {
      REAL8 BaryTime = multiSSB->data[0]->DeltaT->data[0] + GPS2REAL8(doppler->refTime) - StartTime - SFTTimeBaseline/2.0*multiSSB->data[0]->Tdot->data[0];
      Buffer->StartTimeinBaryCenter = BaryTime;
      for(i=1;i<numDetectors;i++)
	{
	  BaryTime = multiSSB->data[i]->DeltaT->data[0] + GPS2REAL8(doppler->refTime) - StartTime - SFTTimeBaseline/2.0*multiSSB->data[i]->Tdot->data[0];
	  if(BaryTime < Buffer->StartTimeinBaryCenter)
	    Buffer->StartTimeinBaryCenter = BaryTime;
	} 
    }

  /* Loop over Detectors */
  for(i=0;i<numDetectors;i++)
    {
      /* Complex Series to store  F_a and F_b integrands in */
      FFTWCOMPLEXSeries *FaIn, *FbIn;
      /* Complex Series to store  F_a and F-b integrands after correcting for spin */
      FFTWCOMPLEXSeries *FaInSpinCorrected,*FbInSpinCorrected;
      /* The plans to fft FaIn and FbIn */
      fftw_plan plan_a,plan_b;
      /* The integrals once calculated are stored in FaOut and FbOut */
      FFTWCOMPLEXSeries *FaOut, *FbOut;

      /* Go through the whole process if not the same sky position */
      if(!SAMESKYPOSITION)
	{
	  /* Seperate the Real and Imaginary Parts */
	  REAL8Sequence* ResampledReal;
	  REAL8Sequence* ResampledImag;
   
	  /* BaryTimes is a sequence containing the times at the barycenter at the centers of each SFT */
	  REAL8Vector* BaryTimes;

	  /* DetectorTimes, like BaryTimes stores the times at the Detector at the centers of each SFT */
	  REAL8Sequence* DetectorTimes;

	  /* Times at which a and b are sampled */
	  REAL8Sequence* AandBTimes;

	  /* Real and Imag are for readability and will be the raw time series from now on */
	  REAL8Sequence* Real = TSeries->Real[i];
	  REAL8Sequence* Imag = TSeries->Imag[i];

	  /* Real_Corrected and Imag_Corrected will store the heterodyne corrected time Series */
	  REAL8Sequence* Real_Corrected = (REAL8Sequence*)XLALCreateREAL8Sequence(TSeries->Real[i]->length);
	  REAL8Sequence* Imag_Corrected = (REAL8Sequence*)XLALCreateREAL8Sequence(TSeries->Imag[i]->length);

	  /* Antenna Patterns */
	  REAL8Sequence* a_at_DetectorTimes = (REAL8Sequence*)XLALCreateREAL8Sequence(multiAMcoef->data[i]->a->length);
	  REAL8Sequence* b_at_DetectorTimes = (REAL8Sequence*)XLALCreateREAL8Sequence(multiAMcoef->data[i]->b->length);

	   /* Detector Times Corresponding to a linear sampling in the Barycentric frame */
	  REAL8Sequence* CorrDetTimes;

	  /* Store the Antenna Patterns */
	  for(p=0;p<multiAMcoef->data[i]->a->length;p++)
	    {
	      a_at_DetectorTimes->data[p] = (REAL8)multiAMcoef->data[i]->a->data[p];
	      b_at_DetectorTimes->data[p] = (REAL8)multiAMcoef->data[i]->b->data[p];
	    }

	  /* Allocate and Store DetectorTimes */
	  DetectorTimes = (REAL8Sequence*)XLALCreateREAL8Sequence(multiSSB->data[i]->DeltaT->length*2);
	  /* Allocate and Store AandBTimes */
	  AandBTimes = (REAL8Sequence*)XLALCreateREAL8Sequence(multiSSB->data[i]->DeltaT->length);
      
	  /* Store the StartTime and EndTime of Each SFT. In order to avoid discontinuities, store an EndTime that is off by 1e-3 seconds */
	  for(p=0;p<DetectorTimes->length;p+=2)
	    {
	      DetectorTimes->data[p] = GPS2REAL8(multiDetStates->data[i]->data[p/2].tGPS) - SFTTimeBaseline/2.0 - StartTime;
	      DetectorTimes->data[p+1] = GPS2REAL8(multiDetStates->data[i]->data[p/2].tGPS) + SFTTimeBaseline/2.0 - StartTime - 1e-3;
	    }
	  
	  /* Store AandBTimes */
	  for(p=0;p<AandBTimes->length;p++)
	    {
	      AandBTimes->data[p] = GPS2REAL8(multiDetStates->data[i]->data[p].tGPS) - StartTime;
	    }
	  
	  /* These will store the Resampled Real and Imaginary parts */
	  ResampledReal = (REAL8Sequence*)XLALCreateREAL8Sequence(length);
	  ResampledImag = (REAL8Sequence*)XLALCreateREAL8Sequence(length);
      
	  /* Allocate and Store BaryTimes */
	  BaryTimes = (REAL8Sequence*)XLALCreateREAL8Sequence(multiSSB->data[i]->DeltaT->length*2);


	  /* doppler->refTime is something called the internal refTime != real refTime */
	  /* Store corresponding BaryTimes */
	  for(p=0;p<BaryTimes->length;p+=2)
	    {
	    BaryTimes->data[p] = multiSSB->data[i]->DeltaT->data[p/2] + GPS2REAL8(doppler->refTime) - StartTime - SFTTimeBaseline/2.0*multiSSB->data[i]->Tdot->data[p/2];
	    BaryTimes->data[p+1] = multiSSB->data[i]->DeltaT->data[p/2] + GPS2REAL8(doppler->refTime) - StartTime + (SFTTimeBaseline/2.0-1e-3)*multiSSB->data[i]->Tdot->data[p/2];
	    
	    }

	  /* Apply the correction term for the heterodyning done to the data */
	  ApplyHetCorrection(BaryTimes,DetectorTimes,Real,Imag,TSeries->Times[i],TSeries, Real_Corrected, Imag_Corrected);

	  /* Resample Real and Image and store in ResampledReal and ResampledImag , also return the calculated CorrDetTimes */
	  CorrDetTimes = ResampleSeries(Real_Corrected,Imag_Corrected,ResampledReal,ResampledImag,dt,BaryTimes,DetectorTimes,TSeries->Times[i],Buffer->StartTimeinBaryCenter); 

	   /* Store CorrDetTimes in MultiCorrDetTimes */
	  MultiCorrDetTimes->data[i] = CorrDetTimes;

	  /* FaIn and FbIn get allocated */
	  FaIn = XLALCreateFFTWCOMPLEXSeries(new_length);
	  FbIn = XLALCreateFFTWCOMPLEXSeries(new_length);
	  for(p=length;p<new_length;p++)
	    {
	      FaIn->data[p][0] = 0;
	      FaIn->data[p][1] = 0;
	      FbIn->data[p][0] = 0;
	      FbIn->data[p][1] = 0;
	    }
	  
	  /* Store in Buffer */
	  Saved_a->data[i] = FaIn;
	  Saved_b->data[i] = FbIn;

	  /* Multiply with A and B */
	  ApplyAandB(CorrDetTimes,AandBTimes,a_at_DetectorTimes,b_at_DetectorTimes,ResampledReal,ResampledImag,FaIn,FbIn,SFTTimeBaseline);

	  XLALDestroyREAL8Sequence(ResampledReal);
	  XLALDestroyREAL8Sequence(ResampledImag);
	  XLALDestroyREAL8Vector(BaryTimes);
	  XLALDestroyREAL8Sequence(DetectorTimes);
	  XLALDestroyREAL8Sequence(AandBTimes);
	  XLALDestroyREAL8Sequence(a_at_DetectorTimes);
	  XLALDestroyREAL8Sequence(b_at_DetectorTimes);
	  XLALDestroyREAL8Sequence(Real_Corrected);
	  XLALDestroyREAL8Sequence(Imag_Corrected);

	}
      /* If Same Sky Postion, Skip all the steps and reuse FaIn and FbIn*/

      else
	{
	  FaIn = Saved_a->data[i];
	  FbIn = Saved_b->data[i];
	}
    
      FaInSpinCorrected = XLALCreateFFTWCOMPLEXSeries(new_length);
      FbInSpinCorrected = XLALCreateFFTWCOMPLEXSeries(new_length);

      for(p=0;p<new_length;p++)
	{
	  FaInSpinCorrected->data[p][0] = 0;
	  FaInSpinCorrected->data[p][1] = 0;
	  FbInSpinCorrected->data[p][0] = 0;
	  FbInSpinCorrected->data[p][1] = 0;
	} 

      /* Allocate Memory for FaOut and FbOut*/
      FaOut = XLALCreateFFTWCOMPLEXSeries(new_length);
      FbOut = XLALCreateFFTWCOMPLEXSeries(new_length);
      
      /* Make Plans */
      plan_a = fftw_plan_dft_1d(FaInSpinCorrected->length,FaInSpinCorrected->data,FaOut->data,FFTW_FORWARD,FFTW_ESTIMATE);
      plan_b = fftw_plan_dft_1d(FbInSpinCorrected->length,FbInSpinCorrected->data,FbOut->data,FFTW_FORWARD,FFTW_ESTIMATE);

      ApplySpinDowns(&(doppler->fkdot),dt,FaIn,FbIn,Buffer->StartTimeinBaryCenter,MultiCorrDetTimes->data[i],uvar_refTime-StartTime,FaInSpinCorrected,FbInSpinCorrected);

      /* FFT!! */
      fftw_execute(plan_a);
      fftw_execute(plan_b);
       
      for(p=0;p<new_length;p++)
	{
	  Fa_Real->data[p] += FaOut->data[p][0]*dt;
	  Fa_Imag->data[p] += FaOut->data[p][1]*dt;
	  Fb_Real->data[p] += FbOut->data[p][0]*dt;
	  Fb_Imag->data[p] += FbOut->data[p][1]*dt;
	}

      
      XLALDestroyFFTWCOMPLEXSeries(FaOut);
      XLALDestroyFFTWCOMPLEXSeries(FbOut);
      XLALDestroyFFTWCOMPLEXSeries(FaInSpinCorrected);
      XLALDestroyFFTWCOMPLEXSeries(FbInSpinCorrected);
           
    }

  /* Calculate the F-Statistic and store it in a temporary variable */
  for(p=0;p<new_length;p++)
    {
      REAL8 Fa_magsquare = pow(Fa_Real->data[p],2)+pow(Fa_Imag->data[p],2);
      REAL8 Fb_magsquare = pow(Fb_Real->data[p],2)+pow(Fb_Imag->data[p],2);
      REAL8 CrossTerm = Fa_Real->data[p]*Fb_Real->data[p]+Fa_Imag->data[p]*Fb_Imag->data[p];
      Fstat_temp->data[p] = Dd_inv*((Bd*Fa_magsquare)+(Ad*Fb_magsquare)-(2*Cd*CrossTerm));
    }

  /* Store the Values in the appropriate spot on the fstatvector */
  {
    UINT4 fmin_index = 0;
    UINT4 fstatVectorlength = fstatVector->data->length;
    UINT4 q,r;
    if(fstatVectorlength > new_length)
      {
	fprintf(stderr," fstatVector's length is greater than total number of bins calculated. Something went wrong allocating fstatVector \n");
	exit(1);
      }
    if(TSeries->f_het > fstatVector->f0 && TSeries->f_het < fstatVector->f0 + uvar_FreqBand)
      {
	fmin_index = floor((TSeries->f_het-fstatVector->f0)/dF_closest + 0.5);
	q = 0;
	for(r=new_length-fmin_index;r<new_length;r++)
	  fstatVector->data->data[q++] = Fstat_temp->data[r];
	r = 0;
	while(q<fstatVectorlength)
	  fstatVector->data->data[q++] = Fstat_temp->data[r++];
      }
    else if(TSeries->f_het < fstatVector->f0)
      {
	fmin_index = floor((fstatVector->f0-TSeries->f_het)/dF_closest + 0.5);
	q = 0;
	for(r = fmin_index;r<fstatVector->data->length+fmin_index;r++)
	  fstatVector->data->data[q++] = Fstat_temp->data[r];
      }
    else
      {
	fmin_index = floor((TSeries->f_het-fstatVector->f0)/dF_closest + 0.5);
	q = 0;
	for(r=new_length-fmin_index;r<new_length-fmin_index+fstatVectorlength;r++)
	  fstatVector->data->data[q++] = Fstat_temp->data[r];
	
      }
  }
  
  XLALDestroyREAL8Sequence(Fa_Real);
  XLALDestroyREAL8Sequence(Fb_Real);
  XLALDestroyREAL8Sequence(Fa_Imag);
  XLALDestroyREAL8Sequence(Fb_Imag);
  XLALDestroyREAL8Sequence(Fstat_temp);

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

