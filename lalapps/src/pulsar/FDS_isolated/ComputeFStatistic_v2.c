/*
 * Copyright (C) 2008 Karl Wette
 * Copyright (C) 2007 Chris Messenger
 * Copyright (C) 2007 Reinhard Prix
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
/** \author R. Prix, I. Gholami, Y. Ioth, Papa, X. Siemens, C. Messenger, K. Wette
 * \file
 * \brief
 * Calculate the F-statistic for a given parameter-space of pulsar GW signals.
 * Implements the so-called "F-statistic" as introduced in \ref JKS98.
 *
 * This code is a descendant of an earlier implementation 'ComputeFStatistic.[ch]'
 * by Bruce Allen, Bernd Machenschalk, David Hammer, Jolien Creighton, Maria Alessandra Papa,
 * Reinhard Prix, Xavier Siemens, Scott Koranda, Yousuke Itoh
 *
 *
 *********************************************************************************/
#include "config.h"

/* System includes */
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <strings.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif


int finite(double);

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/GSLSupport.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/ExtrapolatePulsarSpins.h>

#include <lal/NormalizeSFTRngMed.h>
#include <lal/ComputeFstat.h>
#include <lal/LALHough.h>

#include <lal/LogPrintf.h>
#include <lal/DopplerFullScan.h>
#include <lal/BinaryPulsarTiming.h>

#include <lalapps.h>

#include <lal/lalGitID.h>
#include <lalappsGitID.h>

/* local includes */
#include "HeapToplist.h"

#include "ComputeFstatGPU.h"

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

/** What info do we want to store in our toplist? */
typedef struct {
  PulsarDopplerParams doppler;		/**< Doppler params of this 'candidate' */
  Fcomponents  Fstat;			/**< the Fstat-value (plus Fa,Fb) for this candidate */
  CmplxAntennaPatternMatrix Mmunu;		/**< antenna-pattern matrix Mmunu = Sinv*Tsft * [ Ad, Cd; Cd; Bd ] */
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


/* ----- User-variables: can be set from config-file or command-line */
typedef struct {
  INT4 Dterms;			/**< number of terms in LALDemod Dirichlet kernel is Dterms+1 */
  CHAR *IFO;			/**< IFO name: only required if using v1 SFTs */
  BOOLEAN SignalOnly;		/**< FALSE: estimate noise-floor from data, TRUE: assume Sh=1 */
  BOOLEAN UseNoiseWeights;	/**< use SFT-specific noise-weights for each segment in Fstat-computation */

  REAL8 Freq;			/**< start-frequency of search */
  REAL8 FreqBand;		/**< Frequency-band to search over */
  REAL8 dFreq;			/**< user-specifyable Frequency stepsize */

  REAL8 Alpha;			/**< equatorial right-ascension in rad */
  CHAR *RA;
  REAL8 dAlpha;
  REAL8 AlphaBand;

  REAL8 Delta;			/**< equatorial declination in rad */
  CHAR *Dec;
  REAL8 dDelta;
  REAL8 DeltaBand;

  REAL8 f1dot;			/**< 1st spindown dFreq/dt */
  REAL8 df1dot;
  REAL8 f1dotBand;

  REAL8 f2dot;			/**< 2nd spindown d^2Freq/dt^2 */
  REAL8 df2dot;
  REAL8 f2dotBand;

  REAL8 f3dot;			/**< 3rd spindown d^3Freq/dt^3 */
  REAL8 df3dot;
  REAL8 f3dotBand;

  /* orbital parameters */
  REAL8 orbitPeriod;		/**< binary-system orbital period in s */
  REAL8 orbitasini;		/**< amplitude of radial motion */
  INT4 orbitTpSSBsec;		/**< time of periapse passage */
  INT4 orbitTpSSBnan;
  REAL8 orbitTpSSBMJD;		/**< in MJD format */
  REAL8 orbitArgp;		/**< angle of periapse */
  REAL8 orbitEcc;		/**< orbital eccentricity */

  /* extra parameters for --gridType=GRID_SPINDOWN_AGEBRK parameter space */
  REAL8 spindownAge;            /**< spindown age of the object */
  REAL8 minBraking;             /**< minimum braking index */
  REAL8 maxBraking;             /**< maximum braking index */

  /* --- */
  REAL8 TwoFthreshold;		/**< output threshold on 2F */
  CHAR *ephemDir;		/**< directory to look for ephemeris files */
  CHAR *ephemYear;		/**< date-range string on ephemeris-files to use */

  INT4  gridType;		/**< type of template grid in parameter space */
  INT4  metricType;		/**< type of metric to use in metric template grids */
  BOOLEAN projectMetric;	/**< should we use frequency-projected metric? */
  REAL8 metricMismatch;		/**< maximal *nominal* metric mismatch *per* dimension */
  CHAR *skyRegion;		/**< list of skypositions defining a search-polygon */
  CHAR *DataFiles;		/**< glob-pattern for SFT data-files to use */

  BOOLEAN help;			/**< output help-string */
  CHAR *outputLogfile;		/**< write a log-file */
  CHAR *outputFstat;		/**< filename to output Fstatistic in */
  CHAR *outputLoudest;		/**< filename for loudest F-candidate plus parameter estimation */
  CHAR *outputLogPrintf;        /**< send output from LogPrintf statements to this file */

  CHAR *outputFstatHist;        /**< output discrete histogram of all Fstatistic values */
  REAL8 FstatHistBin;           /**< width of an Fstatistic histogram bin */

  BOOLEAN countTemplates;       /**< just count templates (if supported) instead of search */

  INT4 NumCandidatesToKeep;	/**< maximal number of toplist candidates to output */
  REAL8 FracCandidatesToKeep;	/**< fractional number of candidates to output in toplist */
  INT4 clusterOnScanline;	/**< number of points on "scanline" to use for 1-D local maxima finding */

  CHAR *gridFile;		/**< read template grid from this file */
  REAL8 dopplermax;		/**< maximal possible doppler-effect */

  INT4 RngMedWindow;		/**< running-median window for noise floor estimation */
  REAL8 refTime;		/**< reference-time for definition of pulsar-parameters [GPS] */
  REAL8 refTimeMJD;		/**< the same in MJD */

  REAL8 internalRefTime;	/**< which reference time to use internally for template-grid */
  INT4 SSBprecision;		/**< full relativistic timing or Newtonian */

  BOOLEAN useRAA;               /**< use rigid adiabatic response instead of long-wavelength */
  BOOLEAN bufferedRAA;		/**< approximate RAA by using only middle frequency */

  INT4 minStartTime;		/**< earliest start-time to use data from */
  INT4 maxEndTime;		/**< latest end-time to use data from */
  CHAR *workingDir;		/**< directory to use for output files */
  REAL8 timerCount;		/**< output progress-meter every <timerCount> templates */

  INT4 upsampleSFTs;		/**< use SFT-upsampling by this factor */

  BOOLEAN GPUready;		/**< use special single-precision  'GPU-ready' version */

  BOOLEAN version;		/**< output version information */
} UserInput_t;

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/* empty initializers */
static UserInput_t empty_UserInput;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);
void initUserVars (LALStatus *, UserInput_t *uvar);
void InitFStat ( LALStatus *, ConfigVariables *cfg, const UserInput_t *uvar );
void Freemem(LALStatus *,  ConfigVariables *cfg);

void checkUserInputConsistency (LALStatus *, const UserInput_t *uvar);
int outputBeamTS( const CHAR *fname, const AMCoeffs *amcoe, const DetectorStateSeries *detStates );
void InitEphemeris (LALStatus *, EphemerisData *edat, const CHAR *ephemDir, const CHAR *ephemYear, LIGOTimeGPS epoch, BOOLEAN isLISA);
void getUnitWeights ( LALStatus *, MultiNoiseWeights **multiWeights, const MultiSFTVector *multiSFTs );

int write_FstatCandidate_to_fp ( FILE *fp, const FstatCandidate *thisFCand );
int write_PulsarCandidate_to_fp ( FILE *fp,  const PulsarCandidate *pulsarParams, const FstatCandidate *Fcand );

int compareFstatCandidates ( const void *candA, const void *candB );

void WriteFStatLog ( LALStatus *status, const CHAR *log_fname, const CHAR *logstr );
void getLogString ( LALStatus *status, CHAR **logstr, const ConfigVariables *cfg );
void OutputVersion ( void );

const char *va(const char *format, ...);	/* little var-arg string helper function */
CHAR *append_string ( CHAR *str1, const CHAR *append );

/* ---------- scanline window functions ---------- */
scanlineWindow_t *XLALCreateScanlineWindow ( UINT4 windowWings );
void XLALDestroyScanlineWindow ( scanlineWindow_t *scanlineWindow );
int XLALAdvanceScanlineWindow ( const FstatCandidate *nextCand, scanlineWindow_t *scanWindow );
BOOLEAN XLALCenterIsLocalMax ( const scanlineWindow_t *scanWindow );


/*---------- empty initializers ---------- */
static const ConfigVariables empty_ConfigVariables;
static const FstatCandidate empty_FstatCandidate;

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
  static const char *fn = "main()";
  LALStatus status = blank_status;	/* initialize status */

  FILE *fpFstat = NULL;
  ComputeFBuffer cfBuffer = empty_ComputeFBuffer;
  ComputeFBuffer_REAL4 cfBuffer4 = empty_ComputeFBuffer_REAL4;
  REAL8 numTemplates, templateCounter;
  REAL8 tickCounter;
  time_t clock0;
  Fcomponents Fstat = empty_Fcomponents;
  PulsarDopplerParams dopplerpos = empty_PulsarDopplerParams;		/* current search-parameters */
  FstatCandidate loudestFCand = empty_FstatCandidate, thisFCand = empty_FstatCandidate;
  BinaryOrbitParams *orbitalParams = NULL;
  FILE *fpLogPrintf = NULL;
  gsl_vector_int *Fstat_histogram = NULL;

  UserInput_t uvar = empty_UserInput;
  ConfigVariables GV = empty_ConfigVariables;		/**< global container for various derived configuration settings */

  lalDebugLevel = 0;
  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LAL_CALL (initUserVars(&status, &uvar), &status); 	

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput(&status, argc, argv), &status);	

  if (uvar.help)	/* if help was requested, we're done here */
    exit (0);

  if ( uvar.version )
    {
      OutputVersion();
      exit (0);
    }

  /* set log-level and open log-file */
  LogSetLevel ( lalDebugLevel );
  if (LALUserVarWasSet(&uvar.outputLogPrintf)) {
    if ((fpLogPrintf = fopen(uvar.outputLogPrintf, "wb")) == NULL) {
      LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputLogPrintf);
      return (COMPUTEFSTATISTIC_ESYS);
    }
    LogSetFile(fpLogPrintf);
  }

  /* do some sanity checks on the user-input before we proceed */
  LAL_CALL ( checkUserInputConsistency(&status, &uvar), &status);

  /* Initialization the common variables of the code, */
  /* like ephemeries data and template grids: */
  LAL_CALL ( InitFStat(&status, &GV, &uvar), &status);

  /* ----- produce a log-string describing the specific run setup ----- */
  LAL_CALL ( getLogString ( &status, &(GV.logstring), &GV ), &status );
  LogPrintfVerbatim( LOG_DEBUG, "%s", GV.logstring );

  /* keep a log-file recording all relevant parameters of this search-run */
  if ( uvar.outputLogfile ) {
    LAL_CALL (WriteFStatLog ( &status, uvar.outputLogfile, GV.logstring ), &status );
  }

  /* if a complete output of the F-statistic file was requested,
   * we open and prepare the output-file here */
  if (uvar.outputFstat)
    {
      if ( (fpFstat = fopen (uvar.outputFstat, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputFstat);
	  return (COMPUTEFSTATISTIC_ESYS);
	}

      fprintf (fpFstat, "%s", GV.logstring );
    } /* if outputFstat */

  /* start Fstatistic histogram with a single empty bin */
  if (uvar.outputFstatHist) {
    if ((Fstat_histogram = gsl_vector_int_alloc(1)) == NULL) {
      LALPrintError("\nCouldn't allocate 'Fstat_histogram'\n");
      return COMPUTEFSTATISTIC_EMEM;
    }
    gsl_vector_int_set_zero(Fstat_histogram);
  }

  /* setup binary parameters */
  orbitalParams = NULL;
  if ( LALUserVarWasSet(&uvar.orbitasini) && (uvar.orbitasini > 0) )
    {
      orbitalParams = (BinaryOrbitParams *)LALCalloc(1,sizeof(BinaryOrbitParams));
      orbitalParams->tp.gpsSeconds = uvar.orbitTpSSBsec;
      orbitalParams->tp.gpsNanoSeconds = uvar.orbitTpSSBnan;
      orbitalParams->argp = uvar.orbitArgp;
      orbitalParams->asini = uvar.orbitasini;
      orbitalParams->ecc = uvar.orbitEcc;
      orbitalParams->period = uvar.orbitPeriod;
      if (LALUserVarWasSet(&uvar.orbitTpSSBMJD))
	{
	  /* convert MJD peripase to GPS using Matt Pitkins code found at lal/packages/pulsar/src/BinaryPulsarTimeing.c */
	  REAL8 GPSfloat;
	  GPSfloat = LALTDBMJDtoGPS(uvar.orbitTpSSBMJD);
	  XLALGPSSetREAL8(&(orbitalParams->tp),GPSfloat);
	}
      else
	{
	  orbitalParams->tp.gpsSeconds = uvar.orbitTpSSBsec;
	  orbitalParams->tp.gpsNanoSeconds = uvar.orbitTpSSBnan;
	}
    }

  /* count number of templates */
  numTemplates = XLALNumDopplerTemplates ( GV.scanState );
  if (uvar.countTemplates)
    printf("%%%% Number of templates: %0.0f\n", numTemplates);    

  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
   */
  templateCounter = 0.0; 
  tickCounter = 0;
  clock0 = time(NULL);

  /* skip search if user supplied --countTemplates */
  while ( !uvar.countTemplates && (XLALNextDopplerPos( &dopplerpos, GV.scanState ) == 0) )
    {
      dopplerpos.orbit = orbitalParams;		/* temporary solution until binary-gridding exists */
      
      /* main function call: compute F-statistic for this template */
      if ( ! uvar.GPUready )
        {
          LAL_CALL( ComputeFStat(&status, &Fstat, &dopplerpos, GV.multiSFTs, GV.multiNoiseWeights,
                                       GV.multiDetStates, &GV.CFparams, &cfBuffer ), &status );
        }
      else
        {
          REAL4 F;

          XLALDriverFstatGPU ( &F, &dopplerpos, GV.multiSFTs, GV.multiNoiseWeights, GV.multiDetStates, GV.CFparams.Dterms, &cfBuffer4 );
          if ( xlalErrno ) {
            XLALPrintError ("%s: XLALDriverFstatGPU() failed with errno=%d\n", fn, xlalErrno );
            return xlalErrno;
          }
          /* this function only returns F, not Fa, Fb */
          Fstat = empty_Fcomponents;
          Fstat.F = F;

        } /* if GPUready==true */

      /* Progress meter */
      templateCounter += 1.0;
      if ( lalDebugLevel && ( ++tickCounter > uvar.timerCount) )
	{
	  REAL8 diffSec = time(NULL) - clock0 ;  /* seconds since start of loop*/
	  REAL8 taup = diffSec / templateCounter ;
	  REAL8 timeLeft = (numTemplates - templateCounter) *  taup;
	  tickCounter = 0.0;
	  LogPrintf (LOG_DEBUG, "Progress: %g/%g = %.2f %% done, Estimated time left: %.0f s\n",
		     templateCounter, numTemplates, templateCounter/numTemplates * 100.0, timeLeft);
	}
      
      /* sanity check on the result */
      if ( !finite(Fstat.F) )
	{
	  LogPrintf(LOG_CRITICAL, "non-finite F = %.16g, Fa=(%.16g,%.16g), Fb=(%.16g,%.16g)\n", 
		    Fstat.F, Fstat.Fa.re, Fstat.Fa.im, Fstat.Fb.re, Fstat.Fb.im );
	  LogPrintf (LOG_CRITICAL, "[Alpha,Delta] = [%.16g,%.16g],\nfkdot=[%.16g,%.16g,%.16g,%16.g]\n",
		     dopplerpos.Alpha, dopplerpos.Delta, 
		     dopplerpos.fkdot[0], dopplerpos.fkdot[1], dopplerpos.fkdot[2], dopplerpos.fkdot[3] );
	  if (dopplerpos.orbit != NULL) 
	    {
	      LogPrintf (LOG_CRITICAL, "tp = {%d s, %d ns}, argp = %f, asini = %f, ecc = %f, period = %f\n", 
			 dopplerpos.orbit->tp.gpsSeconds, dopplerpos.orbit->tp.gpsNanoSeconds, 
			 dopplerpos.orbit->argp, dopplerpos.orbit->asini,
			 dopplerpos.orbit->ecc, dopplerpos.orbit->period);
	    }
	  return -1;
	}
      
      /* propagate fkdot from internalRefTime back to refTime for outputting results */
      /* FIXE: only do this for candidates we're going to write out */
      LAL_CALL ( LALExtrapolatePulsarSpins ( &status, dopplerpos.fkdot, GV.refTime, dopplerpos.fkdot, GV.searchRegion.refTime ), &status );
      dopplerpos.refTime = GV.refTime;

      /* collect data on current 'Fstat-candidate' */
      thisFCand.doppler = dopplerpos;
      if ( !uvar.GPUready )
        {
          if ( cfBuffer.multiCmplxAMcoef ) {
            thisFCand.Mmunu = cfBuffer.multiCmplxAMcoef->Mmunu;
          } else {
            thisFCand.Mmunu.Ad = cfBuffer.multiAMcoef->Mmunu.Ad;
            thisFCand.Mmunu.Bd = cfBuffer.multiAMcoef->Mmunu.Bd;
            thisFCand.Mmunu.Cd = cfBuffer.multiAMcoef->Mmunu.Cd;
            thisFCand.Mmunu.Sinv_Tsft = cfBuffer.multiAMcoef->Mmunu.Sinv_Tsft;
            thisFCand.Mmunu.Ed = 0.0;
          }
        }
      else
        {
          thisFCand.Mmunu.Ad = cfBuffer4.multiAMcoef->Mmunu.Ad;
          thisFCand.Mmunu.Bd = cfBuffer4.multiAMcoef->Mmunu.Bd;
          thisFCand.Mmunu.Cd = cfBuffer4.multiAMcoef->Mmunu.Cd;
          thisFCand.Mmunu.Sinv_Tsft = cfBuffer4.multiAMcoef->Mmunu.Sinv_Tsft;
          thisFCand.Mmunu.Ed = 0.0;
        }

      /* correct normalization in --SignalOnly case:
       * we didn't normalize data by 1/sqrt(Tsft * 0.5 * Sh) in terms of 
       * the single-sided PSD Sh: the SignalOnly case is characterized by
       * setting Sh->1, so we need to divide Fa,Fb by sqrt(0.5*Tsft) and F by (0.5*Tsft)
       */
      if ( uvar.SignalOnly )
	{
	  REAL8 norm = 1.0 / sqrt( 0.5 * GV.Tsft );
	  Fstat.Fa.re *= norm;  Fstat.Fa.im *= norm;
	  Fstat.Fb.re *= norm;  Fstat.Fb.im *= norm;
	  Fstat.F *= norm * norm;
	  Fstat.F += 2;		/* compute E[2F]:= 4 + SNR^2 */
	  thisFCand.Mmunu.Sinv_Tsft = GV.Tsft; 
	} 
      thisFCand.Fstat   = Fstat;

      /* push new value onto scan-line buffer */
      XLALAdvanceScanlineWindow ( &thisFCand, GV.scanlineWindow );

      /* two types of threshold: fixed (TwoFThreshold) and dynamic (NumCandidatesToKeep) */
      if ( XLALCenterIsLocalMax ( GV.scanlineWindow ) 					/* must be 1D local maximum */
	   && (2.0 * GV.scanlineWindow->center->Fstat.F >= uvar.TwoFthreshold) )	/* fixed threshold */
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
      if (uvar.outputFstatHist) {

	/* compute bin */
	const size_t bin = 2.0 * thisFCand.Fstat.F / uvar.FstatHistBin;

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

  if ( fpFstat )
    {
      fprintf (fpFstat, "%%DONE\n");
      fclose (fpFstat);
      fpFstat = NULL;
    }

  /* ----- estimate amplitude-parameters for the loudest canidate and output into separate file ----- */
  if ( uvar.outputLoudest )
    {
      FILE *fpLoudest;
      PulsarCandidate pulsarParams = empty_PulsarCandidate;
      pulsarParams.Doppler = loudestFCand.doppler;

      LAL_CALL(LALEstimatePulsarAmplitudeParams (&status, &pulsarParams, &loudestFCand.Fstat, &GV.searchRegion.refTime, &loudestFCand.Mmunu ), 
	       &status );

      if ( (fpLoudest = fopen (uvar.outputLoudest, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputLoudest);
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
  if (uvar.outputFstatHist) {

    size_t i = 0;
    FILE *fpFstatHist = fopen(uvar.outputFstatHist, "wb");
   
    if (fpFstatHist == NULL) {
      LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputFstat);
      return (COMPUTEFSTATISTIC_ESYS);
    }
    fprintf(fpFstatHist, "%s", GV.logstring);
    
    for (i = 0; i < Fstat_histogram->size; ++i)
      fprintf(fpFstatHist, "%0.3g %0.3g %i\n",
	      uvar.FstatHistBin * i,
	      uvar.FstatHistBin * (i + 1),
	      gsl_vector_int_get(Fstat_histogram, i));
    
    fprintf(fpFstatHist, "%%DONE\n");
    fclose(fpFstatHist);
    
  }

  /* Free memory */
  LogPrintf (LOG_DEBUG, "Freeing Doppler grid ... ");
  LAL_CALL ( FreeDopplerFullScan(&status, &GV.scanState), &status);
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

  XLALEmptyComputeFBuffer ( &cfBuffer );
  XLALEmptyComputeFBuffer_REAL4 ( &cfBuffer4 );

  LAL_CALL ( Freemem(&status, &GV), &status);

  if (Fstat_histogram)
    gsl_vector_int_free(Fstat_histogram);

  /* close log-file */
  if (fpLogPrintf) {
    fclose(fpLogPrintf);
    LogSetFile(fpLogPrintf = NULL);
  }
  
  /* did we forget anything ? */
  LALCheckMemoryLeaks();
    
  return 0;
  
} /* main() */


/** 
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
void
initUserVars (LALStatus *status, UserInput_t *uvar)
{
  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
  uvar->upsampleSFTs = 1;
  uvar->Dterms 	= 16;
  uvar->FreqBand = 0.0;
  uvar->Alpha 	= 0.0;
  uvar->Delta 	= 0.0;
  uvar->RA       = NULL;
  uvar->Dec      = NULL;
  uvar->AlphaBand = 0;
  uvar->DeltaBand = 0;
  uvar->skyRegion = NULL;

  uvar->ephemYear = LALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar->ephemYear, EPHEM_YEARS);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar->ephemDir = LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar->ephemDir, DEFAULT_EPHEMDIR);

  uvar->SignalOnly = FALSE;
  uvar->UseNoiseWeights = TRUE;

  uvar->f1dot     = 0.0;
  uvar->f1dotBand = 0.0;

  /* default step-sizes for GRID_FLAT */
  uvar->dAlpha 	= 0.001;
  uvar->dDelta 	= 0.001;
  uvar->dFreq 	 = 0.0;		/* zero indicates 'not set by user==> i.e. use canonical default */
  uvar->df1dot    = 0.0;
  uvar->df2dot    = 0.0;
  uvar->df3dot    = 0.0;

  /* define default orbital semi-major axis */
  uvar->orbitasini = 0.0;

  uvar->TwoFthreshold = 10.0;
  uvar->NumCandidatesToKeep = 0;
  uvar->FracCandidatesToKeep = 0.0;
  uvar->clusterOnScanline = 0;

  uvar->metricType =  LAL_PMETRIC_NONE;
  uvar->projectMetric = TRUE;
  uvar->gridType = GRID_FLAT;

  uvar->metricMismatch = 0.02;

  uvar->help = FALSE;
  uvar->version = FALSE;

  uvar->outputLogfile = NULL;
  uvar->outputFstat = NULL;
  uvar->outputLoudest = NULL;
  uvar->outputLogPrintf = NULL;

  uvar->outputFstatHist = NULL;
  uvar->FstatHistBin = 0.1;

  uvar->countTemplates = FALSE;

  uvar->gridFile = NULL;

  uvar->dopplermax =  1.05e-4;
  uvar->RngMedWindow = 50;	/* for running-median */

  uvar->SSBprecision = SSBPREC_RELATIVISTIC;

  uvar->useRAA = FALSE;
  uvar->bufferedRAA = FALSE;

  uvar->minStartTime = 0;
  uvar->maxEndTime = LAL_INT4_MAX;

  uvar->workingDir = (CHAR*)LALMalloc(512);
  strcpy(uvar->workingDir, ".");

  uvar->timerCount = 1e5;	/* output a timer/progress count every N templates */

  uvar->spindownAge = 0.0;
  uvar->minBraking = 0.0;
  uvar->maxBraking = 0.0;

  uvar->GPUready = 0;

  /* ---------- register all user-variables ---------- */
  LALregBOOLUserStruct(status, 	help, 		'h', UVAR_HELP,     "Print this message"); 

  LALregREALUserStruct(status, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserStruct(status, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
  LALregSTRINGUserStruct(status,RA, 		 0 , UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in format hh:mm:ss.sss");
  LALregSTRINGUserStruct(status,Dec, 		 0 , UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in format dd:mm:ss.sss");
  LALregREALUserStruct(status, 	Freq, 		'f', UVAR_REQUIRED, "Starting search frequency in Hz");
  LALregREALUserStruct(status, 	f1dot, 		's', UVAR_OPTIONAL, "First spindown parameter  dFreq/dt");
  LALregREALUserStruct(status, 	f2dot, 		 0 , UVAR_OPTIONAL, "Second spindown parameter d^2Freq/dt^2");
  LALregREALUserStruct(status, 	f3dot, 		 0 , UVAR_OPTIONAL, "Third spindown parameter  d^3Freq/dt^2");

  LALregREALUserStruct(status, 	AlphaBand, 	'z', UVAR_OPTIONAL, "Band in alpha (equatorial coordinates) in radians");
  LALregREALUserStruct(status, 	DeltaBand, 	'c', UVAR_OPTIONAL, "Band in delta (equatorial coordinates) in radians");
  LALregREALUserStruct(status, 	FreqBand, 	'b', UVAR_OPTIONAL, "Search frequency band in Hz");
  LALregREALUserStruct(status, 	f1dotBand, 	'm', UVAR_OPTIONAL, "Search-band for f1dot");
  LALregREALUserStruct(status, 	f2dotBand, 	 0 , UVAR_OPTIONAL, "Search-band for f2dot");
  LALregREALUserStruct(status, 	f3dotBand, 	 0 , UVAR_OPTIONAL, "Search-band for f3dot");

  LALregREALUserStruct(status, 	dAlpha, 	'l', UVAR_OPTIONAL, "Resolution in alpha (equatorial coordinates) in radians");
  LALregREALUserStruct(status, 	dDelta, 	'g', UVAR_OPTIONAL, "Resolution in delta (equatorial coordinates) in radians");
  LALregREALUserStruct(status,  dFreq,          'r', UVAR_OPTIONAL, "Frequency resolution in Hz [Default: 1/(2T)]");
  LALregREALUserStruct(status, 	df1dot, 	'e', UVAR_OPTIONAL, "Stepsize for f1dot [Default: 1/(2T^2)");
  LALregREALUserStruct(status, 	df2dot, 	 0 , UVAR_OPTIONAL, "Stepsize for f2dot [Default: 1/(2T^3)");
  LALregREALUserStruct(status, 	df3dot, 	 0 , UVAR_OPTIONAL, "Stepsize for f3dot [Default: 1/(2T^4)");

  LALregREALUserStruct(status, 	orbitPeriod, 	 0,  UVAR_OPTIONAL, "Orbital period in seconds");
  LALregREALUserStruct(status, 	orbitasini, 	 0,  UVAR_OPTIONAL, "Orbital projected semi-major axis (normalised by the speed of light) in seconds [Default: 0.0]");
  LALregINTUserStruct(status, 	orbitTpSSBsec, 	 0,  UVAR_OPTIONAL, "The true time of periapsis in the SSB frame (seconds part) in GPS seconds");
  LALregINTUserStruct(status, 	orbitTpSSBnan, 	 0,  UVAR_OPTIONAL, "The true time of periapsis in the SSB frame (nanoseconds part)");
  LALregREALUserStruct(status, 	orbitTpSSBMJD, 	 0,  UVAR_OPTIONAL, "The true time of periapsis in the SSB frame (in MJD)");
  LALregREALUserStruct(status, 	orbitArgp, 	 0,  UVAR_OPTIONAL, "The orbital argument of periapse in radians");
  LALregREALUserStruct(status, 	orbitEcc, 	 0,  UVAR_OPTIONAL, "The orbital eccentricity");
  
  LALregSTRINGUserStruct(status,skyRegion, 	'R', UVAR_OPTIONAL, "ALTERNATIVE: Specify sky-region by polygon (or use 'allsky')");
  LALregSTRINGUserStruct(status,DataFiles, 	'D', UVAR_REQUIRED, "File-pattern specifying (multi-IFO) input SFT-files"); 
  LALregSTRINGUserStruct(status,IFO, 		'I', UVAR_OPTIONAL, "Detector: 'G1', 'L1', 'H1', 'H2' ...(useful for single-IFO v1-SFTs only!)");
  LALregSTRINGUserStruct(status,ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserStruct(status,ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregBOOLUserStruct(status, 	SignalOnly, 	'S', UVAR_OPTIONAL, "Signal only flag");
  LALregBOOLUserStruct(status, 	UseNoiseWeights,'W', UVAR_OPTIONAL, "Use SFT-specific noise weights");

  LALregREALUserStruct(status, 	TwoFthreshold,	'F', UVAR_OPTIONAL, "Set the threshold for selection of 2F");
  LALregINTUserStruct(status, 	gridType,	 0 , UVAR_OPTIONAL, "Grid: 0=flat, 1=isotropic, 2=metric, 3=skygrid-file, 6=grid-file, 7=An*lattice, 8=spin-square, 9=spin-age-brk");
  LALregINTUserStruct(status, 	metricType,	'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserStruct(status, 	metricMismatch,	'X', UVAR_OPTIONAL, "Maximal allowed mismatch for metric tiling");
  LALregSTRINGUserStruct(status,outputLogfile,	 0,  UVAR_OPTIONAL, "Name of log-file identifying the code + search performed");
  LALregSTRINGUserStruct(status,gridFile,	 0,  UVAR_OPTIONAL, "Load grid from this file: sky-grid or full-grid depending on --gridType.");
  LALregREALUserStruct(status,	refTime,	 0,  UVAR_OPTIONAL, "SSB reference time for pulsar-parameters [Default: startTime]");
  LALregREALUserStruct(status,	refTimeMJD,	 0,  UVAR_OPTIONAL, "SSB reference time for pulsar-parameters in MJD [Default: startTime]");
  LALregREALUserStruct(status, 	dopplermax, 	'q', UVAR_OPTIONAL, "Maximum doppler shift expected");  

  LALregSTRINGUserStruct(status,outputFstat,	 0,  UVAR_OPTIONAL, "Output-file for F-statistic field over the parameter-space");
  LALregSTRINGUserStruct(status,outputLoudest,	 0,  UVAR_OPTIONAL, "Loudest F-statistic candidate + estimated MLE amplitudes");

  LALregSTRINGUserStruct(status,outputFstatHist, 0,  UVAR_OPTIONAL, "Output-file for a discrete histogram of all Fstatistic values");
  LALregREALUserStruct(status,  FstatHistBin,    0,  UVAR_OPTIONAL, "Width of an Fstatistic histogram bin");

  LALregINTUserStruct(status,  NumCandidatesToKeep,0, UVAR_OPTIONAL, "Number of Fstat 'candidates' to keep. (0 = All)");
  LALregREALUserStruct(status,FracCandidatesToKeep,0, UVAR_OPTIONAL, "Fraction of Fstat 'candidates' to keep.");
  LALregINTUserStruct(status,   clusterOnScanline, 0, UVAR_OPTIONAL, "Neighbors on each side for finding 1D local maxima on scanline");


  LALregINTUserStruct ( status, minStartTime, 	 0,  UVAR_OPTIONAL, "Earliest SFT-timestamp to include");
  LALregINTUserStruct ( status, maxEndTime, 	 0,  UVAR_OPTIONAL, "Latest SFT-timestamps to include");

  LALregBOOLUserStruct( status, version,	'V', UVAR_SPECIAL,  "Output version information");

  /* ----- more experimental/expert options ----- */
  LALregINTUserStruct (status, 	SSBprecision,	 0,  UVAR_DEVELOPER, "Precision to use for time-transformation to SSB: 0=Newtonian 1=relativistic");
  
  LALregBOOLUserStruct(status, 	useRAA, 	 0,  UVAR_DEVELOPER, "Use rigid adiabatic approximation (RAA) for detector response");
  LALregBOOLUserStruct(status, 	bufferedRAA, 	 0,  UVAR_DEVELOPER, "Approximate RAA by using only middle-frequency");

  LALregINTUserStruct(status, 	RngMedWindow,	'k', UVAR_DEVELOPER, "Running-Median window size");
  LALregINTUserStruct(status,	Dterms,		't', UVAR_DEVELOPER, "Number of terms to keep in Dirichlet kernel sum");

  LALregSTRINGUserStruct(status,workingDir,     'w', UVAR_DEVELOPER, "Directory to use as work directory.");
  LALregREALUserStruct(status, 	timerCount, 	 0,  UVAR_DEVELOPER, "N: Output progress/timer info every N templates");  
  LALregREALUserStruct(status,	internalRefTime, 0,  UVAR_DEVELOPER, "internal reference time to use for Fstat-computation [Default: startTime]");

  LALregINTUserStruct(status,	upsampleSFTs,	 0,  UVAR_DEVELOPER, "(integer) Factor to up-sample SFTs by");
  LALregBOOLUserStruct(status, 	projectMetric, 	 0,  UVAR_DEVELOPER, "Use projected metric on Freq=const subspact");

  LALregSTRINGUserStruct(status,outputLogPrintf, 0,  UVAR_DEVELOPER, "Send all output from LogPrintf statements to this file");

  LALregBOOLUserStruct(status, 	countTemplates,  0,  UVAR_DEVELOPER, "Count number of templates (if supported) instead of search"); 

  LALregREALUserStruct(status,  spindownAge,     0,  UVAR_DEVELOPER, "Spindown age for --gridType=9");
  LALregREALUserStruct(status,  minBraking,      0,  UVAR_DEVELOPER, "Minimum braking index for --gridType=9");
  LALregREALUserStruct(status,  maxBraking,      0,  UVAR_DEVELOPER, "Maximum braking index for --gridType=9");

  LALregBOOLUserStruct(status,  GPUready,        0,  UVAR_OPTIONAL,  "Use single-precision 'GPU-ready' core routines");

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

  INITSTATUS( status, "InitEphemeris", rcsid );
  ATTATCHSTATUSPTR (status);

  ASSERT ( edat, status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL );
  ASSERT ( ephemYear, status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL );

  if ( ephemDir )
    {
      if ( isLISA )
	LALSnprintf(EphemEarth, FNAME_LENGTH, "%s/ephemMLDC.dat", ephemDir);
      else
	LALSnprintf(EphemEarth, FNAME_LENGTH, "%s/earth%s.dat", ephemDir, ephemYear);
      
      LALSnprintf(EphemSun, FNAME_LENGTH, "%s/sun%s.dat", ephemDir, ephemYear);
    }
  else
    {
      if ( isLISA )
	LALSnprintf(EphemEarth, FNAME_LENGTH, "ephemMLDC.dat");
      else
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

  edat->leap = XLALGPSLeapSeconds( epoch.gpsSeconds );
  {
    INT4 err = xlalErrno;
    if ( err != XLAL_SUCCESS ) {
      ABORT ( status, err, "XLALLeapSeconds() failed!\n");
    }
  }

  TRY (LALInitBarycenter(status->statusPtr, edat), status);

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* InitEphemeris() */



/** Initialized Fstat-code: handle user-input and set everything up.
 * NOTE: the logical *order* of things in here is very important, so be careful
 */
void
InitFStat ( LALStatus *status, ConfigVariables *cfg, const UserInput_t *uvar )
{
  REAL8 fCoverMin, fCoverMax;	/* covering frequency-band to read from SFTs */
  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_SFTConstraints;
  LIGOTimeGPS minStartTimeGPS = empty_LIGOTimeGPS;
  LIGOTimeGPS maxEndTimeGPS = empty_LIGOTimeGPS;
  PulsarSpinRange spinRangeRef = empty_PulsarSpinRange;

  UINT4 numSFTs;
  LIGOTimeGPS startTime, endTime;
  size_t toplist_length = uvar->NumCandidatesToKeep;

  INITSTATUS (status, "InitFStat", rcsid);
  ATTATCHSTATUSPTR (status);

  /* set the current working directory */
  if(chdir(uvar->workingDir) != 0)
    {
      LogPrintf (LOG_CRITICAL,  "Unable to change directory to workinDir '%s'\n", uvar->workingDir);
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }

  /* use IFO-contraint if one given by the user */
  if ( LALUserVarWasSet ( &uvar->IFO ) )
    if ( (constraints.detector = XLALGetChannelPrefix ( uvar->IFO )) == NULL ) {
      ABORT ( status,  COMPUTEFSTATISTIC_EINPUT,  COMPUTEFSTATISTIC_MSGEINPUT);
    }
  minStartTimeGPS.gpsSeconds = uvar->minStartTime;
  maxEndTimeGPS.gpsSeconds = uvar->maxEndTime;
  constraints.startTime = &minStartTimeGPS;
  constraints.endTime = &maxEndTimeGPS;
  
  /* get full SFT-catalog of all matching (multi-IFO) SFTs */
  LogPrintf (LOG_DEBUG, "Finding all SFTs to load ... ");
  TRY ( LALSFTdataFind ( status->statusPtr, &catalog, uvar->DataFiles, &constraints ), status);    
  LogPrintfVerbatim (LOG_DEBUG, "done. (found %d SFTs)\n", catalog->length);

  if ( constraints.detector ) 
    LALFree ( constraints.detector );

  if ( !catalog || catalog->length == 0 ) 
    {
      LALPrintError ("\nSorry, didn't find any matching SFTs with pattern '%s'!\n\n", uvar->DataFiles );
      ABORT ( status,  COMPUTEFSTATISTIC_EINPUT,  COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* deduce start- and end-time of the observation spanned by the data */
  numSFTs = catalog->length;
  cfg->Tsft = 1.0 / catalog->data[0].header.deltaF;
  startTime = catalog->data[0].header.epoch;
  endTime   = catalog->data[numSFTs-1].header.epoch;
  XLALGPSAdd(&endTime, cfg->Tsft);	/* add on Tsft to last SFT start-time */

  /* ----- get reference-times (from user if given, use startTime otherwise): ----- */
  if ( LALUserVarWasSet(&uvar->refTime)) 
    {
      TRY ( LALFloatToGPS (status->statusPtr, &(cfg->refTime), &uvar->refTime), status);
    } 
  else if (LALUserVarWasSet(&uvar->refTimeMJD)) 
    {
      /* convert MJD peripase to GPS using Matt Pitkins code found at lal/packages/pulsar/src/BinaryPulsarTimeing.c */
      REAL8 GPSfloat;
      GPSfloat = LALTDBMJDtoGPS(uvar->refTimeMJD);
      XLALGPSSetREAL8(&(cfg->refTime),GPSfloat);
    }
  else
    cfg->refTime = startTime;

  { /* ----- prepare spin-range at refTime (in *canonical format*, ie all Bands >= 0) ----- */

    REAL8 fMin, fMax, f1dotMin, f1dotMax, f2dotMin, f2dotMax, f3dotMin, f3dotMax;

    fMin = MYMIN ( uvar->Freq, uvar->Freq + uvar->FreqBand );
    fMax = MYMAX ( uvar->Freq, uvar->Freq + uvar->FreqBand );

    /* Specific to --gridType=GRID_SPINDOWN_AGEBRK parameter space */
    if (uvar->gridType == GRID_SPINDOWN_AGEBRK) {

      /* Set the first and second spindown ranges to the maximum that will be
       * encountered by the age-braking index parameter space. These will *ONLY*
       * be used by older code to, e.g. load the correct band of SFTs, and
       * will *NOT* be used by the tiling code itself.
       * The formulas used here are, with age=a, min braking=n, max braking=N:
       *
       * -f0/((n-1)*a) <= f1 <= -f0/((N-1)*a)
       * n*(f1^2)/f0 <= f2 <= N*(f1^2)/f0
       *
       * f0/f1 are taken as the maximum/minimum value appropriate for getting the
       * maximum/minimum values of f1/f2 (note that f1 is strictly negative here).
       */
      
      f1dotMin = -1.0 * fMin / ((uvar->minBraking - 1.0) * uvar->spindownAge);
      f1dotMax = -1.0 * fMax / ((uvar->maxBraking - 1.0) * uvar->spindownAge);

      f2dotMin = uvar->minBraking * (f1dotMax * f1dotMax) / fMax;
      f2dotMax = uvar->maxBraking * (f1dotMin * f1dotMin) / fMin;
      
    }

    /* Used for all other --gridTypes */
    else {

      f1dotMin = MYMIN ( uvar->f1dot, uvar->f1dot + uvar->f1dotBand );
      f1dotMax = MYMAX ( uvar->f1dot, uvar->f1dot + uvar->f1dotBand );
      
      f2dotMin = MYMIN ( uvar->f2dot, uvar->f2dot + uvar->f2dotBand );
      f2dotMax = MYMAX ( uvar->f2dot, uvar->f2dot + uvar->f2dotBand );
      
    }

    f3dotMin = MYMIN ( uvar->f3dot, uvar->f3dot + uvar->f3dotBand );
    f3dotMax = MYMAX ( uvar->f3dot, uvar->f3dot + uvar->f3dotBand );
    
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

  {/* ----- load the multi-IFO SFT-vectors ----- */
    UINT4 wings = MYMAX(uvar->Dterms, uvar->RngMedWindow/2 +1);	/* extra frequency-bins needed for rngmed, and Dterms */
    REAL8 fMax = (1.0 + uvar->dopplermax) * fCoverMax + wings / cfg->Tsft; /* correct for doppler-shift and wings */
    REAL8 fMin = (1.0 - uvar->dopplermax) * fCoverMin - wings / cfg->Tsft;
    
    LogPrintf (LOG_DEBUG, "Loading SFTs ... ");
    TRY ( LALLoadMultiSFTs ( status->statusPtr, &(cfg->multiSFTs), catalog, fMin, fMax ), status );
    LogPrintfVerbatim (LOG_DEBUG, "done.\n");
    TRY ( LALDestroySFTCatalog ( status->statusPtr, &catalog ), status );
  }
  { /* ----- load ephemeris-data ----- */
    CHAR *ephemDir;
    BOOLEAN isLISA = FALSE;

    cfg->ephemeris = LALCalloc(1, sizeof(EphemerisData));
    if ( LALUserVarWasSet ( &uvar->ephemDir ) )
      ephemDir = uvar->ephemDir;
    else
      ephemDir = NULL;

    /* hack: if first detector is LISA, we load MLDC-ephemeris instead of 'earth' files */
    if ( cfg->multiSFTs->data[0]->data[0].name[0] == 'Z' )
      isLISA = TRUE;

    TRY( InitEphemeris (status->statusPtr, cfg->ephemeris, ephemDir, uvar->ephemYear, startTime, isLISA ), status);
  }

  /* ----- obtain the (multi-IFO) 'detector-state series' for all SFTs ----- */
  TRY ( LALGetMultiDetectorStates ( status->statusPtr, &(cfg->multiDetStates), cfg->multiSFTs, cfg->ephemeris ), status );

  /* ----- normalize SFTs and calculate noise-weights ----- */
  if ( uvar->SignalOnly )
      cfg->multiNoiseWeights = NULL;   /* noiseWeights == NULL is equivalent to unit noise-weights in ComputeFstat() */
  else 
    {
      UINT4 X, alpha;
      MultiPSDVector *rngmed = NULL;
      cfg->multiNoiseWeights = NULL; 
      TRY ( LALNormalizeMultiSFTVect (status->statusPtr, &rngmed, cfg->multiSFTs, uvar->RngMedWindow ), status );
      TRY ( LALComputeMultiNoiseWeights  (status->statusPtr, &(cfg->multiNoiseWeights), rngmed, uvar->RngMedWindow, 0 ), status );
      TRY ( LALDestroyMultiPSDVector (status->statusPtr, &rngmed ), status );
      if ( !uvar->UseNoiseWeights )	/* in that case simply set weights to 1.0 */
	for ( X = 0; X < cfg->multiNoiseWeights->length; X ++ )
	  for ( alpha = 0; alpha < cfg->multiNoiseWeights->data[X]->length; alpha ++ )
	    cfg->multiNoiseWeights->data[X]->data[alpha] = 1.0;
    } /* if ! SignalOnly */

  /* ----- upsample SFTs ----- */
  if ( (lalDebugLevel >= 2) && (uvar->upsampleSFTs > 1) )
  {
    UINT4 X, numDet = cfg->multiSFTs->length;
    LogPrintf (LOG_DEBUG, "Writing original SFTs for debugging ... ");
    for (X=0; X < numDet ; X ++ ) 
      {
	TRY ( LALWriteSFTVector2Dir ( status->statusPtr, cfg->multiSFTs->data[X], "./", "original", "orig"), status );
      }
    LogPrintfVerbatim ( LOG_DEBUG, "done.\n");
  }

  LogPrintf (LOG_DEBUG, "Upsampling SFTs by factor %d ... ", uvar->upsampleSFTs );
  TRY ( upsampleMultiSFTVector ( status->statusPtr, cfg->multiSFTs, uvar->upsampleSFTs, 16 ), status );
  LogPrintfVerbatim (LOG_DEBUG, "done.\n");

  if ( lalDebugLevel >= 2 && (uvar->upsampleSFTs > 1) )
  {
    UINT4 X, numDet = cfg->multiSFTs->length;
    CHAR tag[60];
    sprintf (tag, "upsampled%02d", uvar->upsampleSFTs );
    LogPrintf (LOG_DEBUG, "Writing upsampled SFTs for debugging ... ");
    for (X=0; X < numDet ; X ++ ) 
      {
	TRY ( LALWriteSFTVector2Dir ( status->statusPtr, cfg->multiSFTs->data[X], "./", tag, tag), status );
      }
    LogPrintfVerbatim ( LOG_DEBUG, "done.\n");
  }

  /* define sky position variables from user input */
  if (LALUserVarWasSet(&uvar->RA)) 
    {
      /* use Matt Pitkins conversion code found in lal/packages/pulsar/src/BinaryPulsarTiming.c */
      cfg->Alpha = LALDegsToRads(uvar->RA, "alpha");
    }
  else cfg->Alpha = uvar->Alpha;
  if (LALUserVarWasSet(&uvar->Dec)) 
    {
      /* use Matt Pitkins conversion code found in lal/packages/pulsar/src/BinaryPulsarTiming.c */
      cfg->Delta = LALDegsToRads(uvar->Dec, "delta");
    }
  else cfg->Delta = uvar->Delta;

  { /* ----- set up Doppler region (at internalRefTime) to scan ----- */
    LIGOTimeGPS internalRefTime = empty_LIGOTimeGPS;
    PulsarSpinRange spinRangeInt = empty_PulsarSpinRange;
    BOOLEAN haveAlphaDelta = (LALUserVarWasSet(&uvar->Alpha) && LALUserVarWasSet(&uvar->Delta)) || (LALUserVarWasSet(&uvar->RA) && LALUserVarWasSet(&uvar->Dec));
    
    if (uvar->skyRegion)
      {
	cfg->searchRegion.skyRegionString = (CHAR*)LALCalloc(1, strlen(uvar->skyRegion)+1);
	if ( cfg->searchRegion.skyRegionString == NULL ) {
	  ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
	}
	strcpy (cfg->searchRegion.skyRegionString, uvar->skyRegion);
      }
    else if (haveAlphaDelta)    /* parse this into a sky-region */
      {
	TRY ( SkySquare2String( status->statusPtr, &(cfg->searchRegion.skyRegionString),
				cfg->Alpha, cfg->Delta,	uvar->AlphaBand, uvar->DeltaBand), status);
      }

    if ( LALUserVarWasSet ( &uvar->internalRefTime ) ) {
      TRY ( LALFloatToGPS (status->statusPtr, &(internalRefTime), &uvar->internalRefTime), status);
    }
    else
      internalRefTime = startTime;

    /* spin searchRegion defined by spin-range at *internal* reference-time */
    TRY ( LALExtrapolatePulsarSpinRange (status->statusPtr, &spinRangeInt, internalRefTime, &spinRangeRef ), status );
    cfg->searchRegion.refTime = spinRangeInt.refTime;
    memcpy ( &cfg->searchRegion.fkdot, &spinRangeInt.fkdot, sizeof(spinRangeInt.fkdot) );
    memcpy ( &cfg->searchRegion.fkdotBand, &spinRangeInt.fkdotBand, sizeof(spinRangeInt.fkdotBand) );

  } /* get DopplerRegion */

  /* ----- set computational parameters for F-statistic from User-input ----- */
  cfg->CFparams.Dterms = uvar->Dterms;
  cfg->CFparams.SSBprec = uvar->SSBprecision;
  cfg->CFparams.useRAA = uvar->useRAA;
  cfg->CFparams.bufferedRAA = uvar->bufferedRAA;
  cfg->CFparams.upsampling = 1.0 * uvar->upsampleSFTs; 

  /* ----- set fixed grid step-sizes from user-input for GRID_FLAT ----- */
  cfg->stepSizes.Alpha = uvar->dAlpha;
  cfg->stepSizes.Delta = uvar->dDelta;
  cfg->stepSizes.fkdot[0] = uvar->dFreq;
  cfg->stepSizes.fkdot[1] = uvar->df1dot;
  cfg->stepSizes.fkdot[2] = uvar->df2dot;
  cfg->stepSizes.fkdot[3] = uvar->df3dot;
  cfg->stepSizes.orbit = NULL;

  /* ----- set up scanline-window if requested for 1D local-maximum clustering on scanline ----- */
  if ( (cfg->scanlineWindow = XLALCreateScanlineWindow ( uvar->clusterOnScanline )) == NULL ) {
    ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM ); 
  }

  /* initialize full multi-dimensional Doppler-scanner */
  {
    DopplerFullScanInit scanInit;			/* init-structure for DopperScanner */
    
    scanInit.searchRegion = cfg->searchRegion;
    scanInit.gridType = uvar->gridType;
    scanInit.gridFile = uvar->gridFile;
    scanInit.metricType = uvar->metricType;
    scanInit.projectMetric = uvar->projectMetric;
    scanInit.metricMismatch = uvar->metricMismatch;
    scanInit.stepSizes = cfg->stepSizes;
    scanInit.ephemeris = cfg->ephemeris;		/* used by Ephemeris-based metric */
    scanInit.startTime = cfg->multiDetStates->startTime;
    scanInit.Tspan     = cfg->multiDetStates->Tspan;
    scanInit.Detector  = &(cfg->multiDetStates->data[0]->detector);	/* just use first IFO for metric */

    /* Specific to --gridType=GRID_SPINDOWN_AGEBRK parameter space */
    if (uvar->gridType == GRID_SPINDOWN_AGEBRK) {
      scanInit.extraArgs[0] = uvar->spindownAge;
      scanInit.extraArgs[1] = uvar->minBraking;
      scanInit.extraArgs[2] = uvar->maxBraking;
    }

    LogPrintf (LOG_DEBUG, "Setting up template grid ... ");
    TRY ( InitDopplerFullScan ( status->statusPtr, &cfg->scanState, &scanInit), status); 
    LogPrintfVerbatim (LOG_DEBUG, "template grid ready.\n");
    XLALNumDopplerTemplates ( cfg->scanState );
  }

  /* set number of toplist candidates from fraction if asked to */
  if (0.0 < uvar->FracCandidatesToKeep && uvar->FracCandidatesToKeep <= 1.0) {
    if (XLALNumDopplerTemplates(cfg->scanState) <= 0.0) {
      LogPrintf(LOG_CRITICAL, "Cannot use FracCandidatesToKeep because number of templates was counted to be zero!\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }
    toplist_length = ceil(XLALNumDopplerTemplates(cfg->scanState) * uvar->FracCandidatesToKeep);
  }

  /* ----- set up toplist if requested ----- */
  if ( toplist_length > 0 )
    if ( create_toplist( &(cfg->FstatToplist), toplist_length, sizeof(FstatCandidate), compareFstatCandidates) != 0 ) {
      ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM ); 
    }

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
#define BUFLEN 1024
  CHAR dateStr[BUFLEN], line[BUFLEN];
  CHAR *cmdline = NULL;
  UINT4 i, numDet, numSpins = PULSAR_MAX_SPINS;
  CHAR *ret = NULL;
  CHAR *id1, *id2;

  INITSTATUS( status, "getLogString", rcsid );
  ATTATCHSTATUSPTR (status);

  /* get full commandline describing search*/
  TRY ( LALUserVarGetLog (status->statusPtr, &cmdline,  UVAR_LOGFMT_CMDLINE ), status );
  sprintf (line, "%%%% cmdline: %s\n", cmdline );
  LALFree ( cmdline );
  ret = append_string ( ret, line );

  /* add code version ID (only useful for git-derived versions) */
  id1 = XLALClearLinebreaks ( lalGitID );
  id2 = XLALClearLinebreaks ( lalappsGitID );
  LALSnprintf (line, BUFLEN, "%%%% %s\n%%%% %s\n", id1, id2 );
  LALFree ( id1 );
  LALFree ( id2 );
  ret = append_string ( ret, line );


  numDet = cfg->multiSFTs->length;
  tp = time(NULL);
  sprintf (line, "%%%% Started search: %s", asctime( gmtime( &tp ) ) );
  ret = append_string ( ret, line );
  ret = append_string ( ret, "%% Loaded SFTs: [ " );
  for ( i=0; i < numDet; i ++ )
    {
      sprintf (line, "%s:%d%s",  cfg->multiSFTs->data[i]->data->name, 
	       cfg->multiSFTs->data[i]->length,
	       (i < numDet - 1)?", ":" ]\n");
      ret = append_string ( ret, line );
    }
  utc = *XLALGPSToUTC( &utc, (INT4)GPS2REAL8(cfg->multiDetStates->startTime) );
  strcpy ( dateStr, asctime(&utc) );
  dateStr[ strlen(dateStr) - 1 ] = 0;
  sprintf (line, "%%%% Start GPS time tStart = %12.3f    (%s GMT)\n", 
	   GPS2REAL8(cfg->multiDetStates->startTime), dateStr);
  ret = append_string ( ret, line );
  sprintf (line, "%%%% Total time spanned    = %12.3f s  (%.1f hours)\n", 
	   cfg->multiDetStates->Tspan, cfg->multiDetStates->Tspan/3600 );
  ret = append_string ( ret, line );
  sprintf (line, "%%%% Pulsar-params refTime = %12.3f \n", GPS2REAL8(cfg->refTime) );
  ret = append_string ( ret, line );
  sprintf (line, "%%%% InternalRefTime       = %12.3f \n", GPS2REAL8(cfg->searchRegion.refTime) );
  ret = append_string ( ret, line );
  sprintf (line, "%%%% Spin-range at internalRefTime: " );
  ret = append_string ( ret, line );

  ret = append_string (ret, "fkdot = [ " );
  for (i=0; i < numSpins; i ++ ) 
    {
      sprintf (line, "%.16g:%.16g%s", 
	       cfg->searchRegion.fkdot[i], 
	       cfg->searchRegion.fkdot[i] + cfg->searchRegion.fkdotBand[i], 
	       (i < numSpins - 1)?", ":" ]\n");
      ret = append_string ( ret, line );
    }

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
WriteFStatLog ( LALStatus *status, const CHAR *log_fname, const CHAR *log_string )
{
  FILE *fplog;

  INITSTATUS (status, "WriteFStatLog", rcsid);
  ATTATCHSTATUSPTR (status);

  if ( !log_fname || !log_string ) {	/* no logfile given */
    ABORT (status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL);
  }

  /* prepare log-file for writing */
  if ( (fplog = fopen(log_fname, "wb" )) == NULL) {
    LogPrintf ( LOG_CRITICAL , "Failed to open log-file '%s' for writing.\n\n", log_fname );
    ABORT (status, COMPUTEFSTATISTIC_ESYS, COMPUTEFSTATISTIC_MSGESYS);
  }

  fprintf (fplog, "%%%% LOG-FILE for ComputeFStatistic run\n\n");
  fprintf (fplog, log_string);
  fclose (fplog);


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
checkUserInputConsistency (LALStatus *status, const UserInput_t *uvar)
{

  INITSTATUS (status, "checkUserInputConsistency", rcsid);  
  
  if (uvar->ephemYear == NULL)
    {
      LALPrintError ("\nNo ephemeris year specified (option 'ephemYear')\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }      

  /* check that only alpha OR RA has been set */
  if ( LALUserVarWasSet(&uvar->Alpha) && (LALUserVarWasSet(&uvar->RA)) )
    {
      LALPrintError ("\nInput either Alpha OR RA, not both!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  /* check that only delta OR Dec has been set */
  if ( LALUserVarWasSet(&uvar->Delta) && (LALUserVarWasSet(&uvar->Dec)) )
    {
      LALPrintError ("\nInput either Delta OR Dec, not both!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* check for negative stepsizes in Freq, Alpha, Delta */
  if ( LALUserVarWasSet(&uvar->dAlpha) && (uvar->dAlpha < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dAlpha not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->dDelta) && (uvar->dDelta < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dDelta not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->dFreq) && (uvar->dFreq < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dFreq not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* check that reference time has not been set twice */
  if ( LALUserVarWasSet(&uvar->refTime) && LALUserVarWasSet(&uvar->refTimeMJD) )
    {
      LALPrintError ("\nSet only uvar->refTime OR uvar->refTimeMJD OR leave empty to use SSB start time as Tref!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* binary parameter checks */
  if ( LALUserVarWasSet(&uvar->orbitPeriod) && (uvar->orbitPeriod <= 0) )
    {
      LALPrintError ("\nNegative or zero value of orbital period not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->orbitasini) && (uvar->orbitasini < 0) )
    {
      LALPrintError ("\nNegative value of projected orbital semi-major axis not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
   if ( LALUserVarWasSet(&uvar->orbitTpSSBMJD) && (LALUserVarWasSet(&uvar->orbitTpSSBsec) || LALUserVarWasSet(&uvar->orbitTpSSBnan)))
    {
      LALPrintError ("\nSet only uvar->orbitTpSSBMJD OR uvar->orbitTpSSBsec/nan to specify periapse passage time!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
   if ( LALUserVarWasSet(&uvar->orbitTpSSBMJD) && (uvar->orbitTpSSBMJD < 0) )
    {
      LALPrintError ("\nNegative value of the true time of orbital periapsis not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->orbitTpSSBsec) && (uvar->orbitTpSSBsec < 0) )
    {
      LALPrintError ("\nNegative value of seconds part of the true time of orbital periapsis not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->orbitTpSSBnan) && ((uvar->orbitTpSSBnan < 0) || (uvar->orbitTpSSBnan >= 1e9)) )
    {
      LALPrintError ("\nTime of nanoseconds part the true time of orbital periapsis must lie in range (0, 1e9]!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->orbitArgp) && ((uvar->orbitArgp < 0) || (uvar->orbitArgp >= LAL_TWOPI)) )
    {
      LALPrintError ("\nOrbital argument of periapse must lie in range [0 2*PI)!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->orbitEcc) && (uvar->orbitEcc < 0) )
    {
      LALPrintError ("\nNegative value of orbital eccentricity not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* grid-related checks */
  {
    BOOLEAN haveAlphaBand = LALUserVarWasSet( &uvar->AlphaBand );
    BOOLEAN haveDeltaBand = LALUserVarWasSet( &uvar->DeltaBand );
    BOOLEAN haveSkyRegion, haveAlphaDelta, haveGridFile; 
    BOOLEAN useSkyGridFile, useFullGridFile, haveMetric, useMetric;

    haveSkyRegion  	= (uvar->skyRegion != NULL);
    haveAlphaDelta 	= (LALUserVarWasSet(&uvar->Alpha) && LALUserVarWasSet(&uvar->Delta) ) || (LALUserVarWasSet(&uvar->RA) && LALUserVarWasSet(&uvar->Dec) );
    haveGridFile      	= (uvar->gridFile != NULL);
    useSkyGridFile   	= (uvar->gridType == GRID_FILE_SKYGRID);
    useFullGridFile	= (uvar->gridType == GRID_FILE_FULLGRID);
    haveMetric     	= (uvar->metricType > LAL_PMETRIC_NONE);
    useMetric     	= (uvar->gridType == GRID_METRIC);

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

    /* Specific checks for --gridType=GRID_SPINDOWN_{SQUARE,AGEBRK} parameter spaces */
    if (uvar->gridType == GRID_SPINDOWN_SQUARE || uvar->gridType == GRID_SPINDOWN_AGEBRK) {

      /* Check that no third spindown range were given */
      if (uvar->f3dot != 0.0 || uvar->f3dotBand != 0.0) {
        LALPrintError ("\nERROR: f3dot and f3dotBand cannot be used with gridType={8,9}\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);      
      }

      /* Check that no grid spacings were given */
      if (uvar->df1dot != 0.0 || uvar->df2dot != 0.0 || uvar->df3dot != 0.0) {
        LALPrintError ("\nERROR: df{1,2,3}dot cannot be used with gridType={8,9}\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);      
      }

    }

    /* Specific checks for --gridType=GRID_SPINDOWN_AGEBRK parameter space */
    if (uvar->gridType == GRID_SPINDOWN_AGEBRK) {

      /* Check age and braking indices */
      if (uvar->spindownAge <= 0.0) {
        LALPrintError ("\nERROR: spindownAge must be strictly positive with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);      
      }
      if (uvar->minBraking <= 0.0) {
        LALPrintError ("\nERROR: minBraking must be strictly positive with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);      
      }
      if (uvar->maxBraking <= 0.0) {
        LALPrintError ("\nERROR: minBraking must be strictly positive with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);      
      }
      if (uvar->minBraking >= uvar->maxBraking) {
        LALPrintError ("\nERROR: minBraking must be strictly less than maxBraking with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);      
      }

      /* Check that no first and second spindown ranges were given */
      if (uvar->f1dot != 0.0 || uvar->f1dotBand != 0.0) {
        LALPrintError ("\nERROR: f1dot and f1dotBand cannot be used with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);      
      }
      if (uvar->f2dot != 0.0 || uvar->f2dotBand != 0.0) {
        LALPrintError ("\nERROR: f2dot and f2dotBand cannot be used with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);      
      }

    }

  } /* Grid-related checks */

  /* check NumCandidatesToKeep and FracCandidatesToKeep */
  if (LALUserVarWasSet(&uvar->NumCandidatesToKeep) && LALUserVarWasSet(&uvar->FracCandidatesToKeep)) {
    LALPrintError ("\nERROR: NumCandidatesToKeep and FracCandidatesToKeep are mutually exclusive\n\n");
    ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);      
  }
  if (LALUserVarWasSet(&uvar->FracCandidatesToKeep) && (uvar->FracCandidatesToKeep <= 0.0 || 1.0 < uvar->FracCandidatesToKeep)) {
    LALPrintError ("\nERROR: FracCandidatesToKeep must be greater than 0.0 and less than or equal to 1.0\n\n");
    ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);      
  }

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

  /* Binary parameters */
  if (pulsarParams->Doppler.orbit) 
    {
      fprintf (fp, "orbitPeriod       = % .16g;\n", pulsarParams->Doppler.orbit->period );
      fprintf (fp, "orbitasini        = % .16g;\n", pulsarParams->Doppler.orbit->asini );
      fprintf (fp, "orbitTpSSBsec     = % .8d;\n", pulsarParams->Doppler.orbit->tp.gpsSeconds );
      fprintf (fp, "orbitTpSSBnan     = % .8d;\n", pulsarParams->Doppler.orbit->tp.gpsNanoSeconds );
      fprintf (fp, "orbitArgp         = % .16g;\n", pulsarParams->Doppler.orbit->argp );
      fprintf (fp, "orbitEcc          = % .16g;\n", pulsarParams->Doppler.orbit->ecc );
    }

  /* Amplitude Modulation Coefficients */
  fprintf (fp, "Ad       = % .6g;\n", Fcand->Mmunu.Ad );
  fprintf (fp, "Bd       = % .6g;\n", Fcand->Mmunu.Bd );
  fprintf (fp, "Cd       = % .6g;\n", Fcand->Mmunu.Cd );
  fprintf (fp, "Ed       = % .6g;\n", Fcand->Mmunu.Ed );
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

/** Simply output version information to stdout */
void
OutputVersion ( void )
{
  printf ( "%s\n", lalGitID );
  printf ( "%s\n", lalappsGitID );

  return;

} /* OutputVersion() */

/** Mini helper-function: append string 'str2' to string 'str1',
 * returns pointer to new concatenated string
 */
CHAR *append_string ( CHAR *str1, const CHAR *str2 )
{
  UINT4 len1 = 0, len2 = 0;
  CHAR *outstr;

  if ( str1 )
    len1 = strlen(str1);
  if ( str2 )
    len2 = strlen(str2);

  if ( ( outstr = LALRealloc ( str1, len1 + len2 + 1 ) ) == NULL )
    {
      XLALPrintError ("Seems like we're out of memory!\n");
      return NULL;
    }

  if ( len1 == 0 )
    outstr[0] = 0;

  if ( str2 )
    outstr = strcat ( outstr, str2 );

  return outstr;

} /* append_string() */

