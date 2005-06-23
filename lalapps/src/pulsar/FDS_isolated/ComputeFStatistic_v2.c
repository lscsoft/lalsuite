/*
 * 
 * Copyright (C) 2002, 2003, 2004 M.A. Papa, X. Siemens, Y. Itoh
 * Copyright (C) 2004, 2005 Reinhard Prix
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
/** \author Y. Ioth, C. Messenger, M.A. Papa, R.Prix, X. Siemens 
 * \file 
 * \brief
 * Calculate the F-statistic for a given parameter-space of pulsar GW signals.
 * Implements the so-called "F-statistic" as introduced in JKS98.
 *                                                                          
 *                                                                          
 *                 Albert Einstein Institute/UWM - started September 2002   
 *********************************************************************************/
#include "config.h"

/* System includes */
#include <stdio.h>
#define __USE_ISOC99 1
#include <math.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/RngMedBias.h>
#include <lal/LALComputeAM.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/PulsarDataTypes.h>
#include <lal/SFTfileIO.h>

#include <lalapps.h>

/* local includes */
#include "clusters.h"
#include "DopplerScan.h"


/* this is defined in C99 and *should* be in math.h.  Long term
   protect this with a HAVE_FINITE */
int finite(double);


RCSID( "$Id$");


#define LD_SMALL        (1.0e-9 / LAL_TWOPI)
#define OOTWOPI         (1.0 / LAL_TWOPI)

#define TWOPI_FLOAT     6.28318530717958f  /* 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/* 1 / (2pi) */ 


#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

/** local parameter-type for computeFStat() */
typedef struct {
  REAL8Vector	*fkdot;		/**< vector of frequency + derivatives (spindowns) */
  REAL8Vector	*DeltaT;	/**< vector of DeltaT_alpha's (depend on skyposition)*/
  REAL8Vector	*Tdot;		/**< vector of Tdot_alpha's (depend on skyposition)*/ 
  AMCoeffs      *amcoe;         /**< Amplitude Modulation coefficients */
  INT4          Dterms;         /**< Terms used in the computation of the dirichlet kernel*/
  REAL8 	tSFT;		/**< length of an SFT in seconds */
} computeFStatPar;


/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  CHAR EphemEarth[MAXFILENAMELENGTH];	/**< filename of earth-ephemeris data */
  CHAR EphemSun[MAXFILENAMELENGTH];	/**< filename of sun-ephemeris data */
  UINT4 FreqImax;  		/**< number of frequency-bins to compute F-stat for */
  UINT4 SpinImax;		/**< number of spindown-bins to compute F for */
  UINT4 spdnOrder;		/**< highest spindown orders (0=none) */ 
  REAL8 dFreq;			/**< frequency resolution */
  REAL8 df1dot;			/**< spindown resolution (f1 = df/dt!!) */
  LIGOTimeGPS Tstart;		/**< start time of observation */
  REAL8 tSFT;			/**< length of an SFT in seconds */
  REAL8 tObs;			/**< total observation time in seconds */
  UINT4 SFTno;			/**< number of SFTs in input */
  UINT4 nsamples;		/**< number of frequency-bins in an SFT */
  LALDetector Detector;         /**< Our detector*/
  EphemerisData *edat;		/**< ephemeris data (from LALInitBarycenter()) */
  CHAR *skyRegionString;	/**< sky-region to search (polygon defined by list of points) */
  LIGOTimeGPSVector timestamps; /**< SFT timestamps */
  LIGOTimeGPSVector midTS;	/**< midpoints of SFT's */
  computeFStatPar *CFSparams;  	/**< Demodulation parameters for computeFStat() */
} ConfigVariables;

/** FIXME: OBSOLETE: used to hold result from LALDemod(). 
    Only kept for the moment to make things work (FIXME)*/
typedef struct {
  const REAL8         *F;            /* Array of value of the F statistic */
  const COMPLEX16     *Fa;           /* Results of match filter with a(t) */
  const COMPLEX16     *Fb;           /* Results of match filter with b(t) */
} LALFstat;


/** Local type for storing F-statistic output from NewLALDemod().
 * Note that length has to be the same for all vectors, anything else
 * is considered an error.
 */
typedef struct {
  UINT4 length;		/**< number of frequency-bins */
  REAL8 f0;		/**< lowest frequency in the band */
  REAL8 fBand;		/**< user-requested frequency-band */
  REAL8 df;		/**< frequency-resolution */
  REAL8Vector *F;	/**< Values of the F-statistic proper (F) over frequency-bins */
  COMPLEX16Vector *Fa;	/**< Values Fa over frequency-bins */
  COMPLEX16Vector *Fb;	/**< Values Fb */
} FStatisticVector;


typedef struct {
  COMPLEX16 Fa;
  COMPLEX16 Fb;
} Fcomponents;


/*----------------------------------------------------------------------*/
/* conditional compilation-switches */

/*----------------------------------------------------------------------*/
/* Error-codes */

#define COMPUTEFSTATC_ENULL 		1
#define COMPUTEFSTATC_ESYS     		2
#define COMPUTEFSTATC_EINPUT   		3
#define COMPUTEFSTATC_EMEM   		4
#define COMPUTEFSTATC_ENONULL 		5
#define COMPUTEFSTATC_EXLAL		6

#define COMPUTEFSTATC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATC_MSGESYS		"System call failed (probably file IO)"
#define COMPUTEFSTATC_MSGEINPUT   	"Invalid input"
#define COMPUTEFSTATC_MSGEMEM   	"Out of memory. Bad."
#define COMPUTEFSTATC_MSGENONULL 	"Output pointer is non-NULL"
#define COMPUTEFSTATC_MSGEXLAL		"XLALFunction-call failed"

/*----------------------------------------------------------------------
 * User-variables: can be set from config-file or command-line */

INT4 uvar_Dterms;
CHAR* uvar_IFO;
/*================*/
CHAR* uvar_IFO2;
/*================*/
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
CHAR *uvar_DataDir;
/*================*/
CHAR *uvar_DataDir2;
/*================*/
CHAR *uvar_mergedSFTFile;
CHAR *uvar_BaseName;
/*================*/
CHAR *uvar_BaseName2;
/*================*/
BOOLEAN uvar_help;
CHAR *uvar_outputLabel;
CHAR *uvar_outputFstat;
REAL8 uvar_FstatMin;
CHAR *uvar_skyGridFile;
CHAR *uvar_outputSkyGrid;
CHAR *uvar_workingDir;
REAL8 uvar_dopplermax;
INT4 uvar_windowsize;

/* interpolation stuff */
BOOLEAN uvar_useEphemeris;

/*----------------------------------------------------------------------*/
/* some other global variables (FIXME) */
SFTVector *SFTvect = NULL;	/**< holds the SFT-data to analyze */
LALFstat Fstat;			/**< output from LALDemod(): F-statistic and amplitudes Fa and Fb */
REAL8 Alpha,Delta;		/**< sky-position currently searched (equatorial coords, radians) */
Clusters HFLines;		/**< stores information about outliers/clusters in F-statistic */
Clusters HPLines;		/**< stores information about outliers/clusters in SFT-power spectrum */

Clusters HFLines, HPLines;

Clusters *highSpLines=&HPLines, *highFLines=&HFLines;

REAL8 medianbias=1.0;		/**< bias in running-median depending on window-size 
				 * (set in NormaliseSFTDataRngMdn()) */

FILE *fpstat=NULL;		/**< output-file: F-statistic candidates and cluster-information */

ConfigVariables GV;		/**< global container for various derived configuration settings */


/* ---------- local prototypes ---------- */

int main(int argc,char *argv[]);
void initUserVars (LALStatus *);
void InitFStat (LALStatus *, ConfigVariables *cfg);
void CreateDemodParams (LALStatus *, computeFStatPar *CFSparams, const SkyPosition *pos,
			const ConfigVariables *cfg);

void CreateNautilusDetector (LALStatus *, LALDetector *Detector);
void Freemem(LALStatus *,  ConfigVariables *cfg);

void EstimateFLines(LALStatus *, const FStatisticVector *FVect);
void NormaliseSFTDataRngMdn (LALStatus *);
INT4 writeFLines(INT4 *maxIndex, REAL8 f0, REAL8 df);
int compare(const void *ip, const void *jp);
INT4 writeFaFb(INT4 *maxIndex);
void WriteFStatLog (LALStatus *, CHAR *argv[]);


void writeFVect(LALStatus *, const FStatisticVector *FVect, const CHAR *fname);

const char *va(const char *format, ...);	/* little var-arg string helper function */

int
XLALNewLALDemod(Fcomponents *FaFb,
		const SFTVector *sfts, 
		const computeFStatPar *params);

int 
XLALcomputeFStat (REAL8 *Fval, const SFTVector *sfts, const computeFStatPar *params);

void
LALGetSSBtimes (LALStatus *status, 
		REAL8Vector *DeltaT,
		REAL8Vector *Tdot,
		LIGOTimeGPS refTime,
		const LIGOTimeGPSVector *GPStimes,
		const SkyPosition *pos,
		const LALDetector *site,
		const EphemerisData *ephem);


/*----------------------------------------------------------------------*/
/* some local defines */

#define EPHEM_YEARS  "00-04"

#define TRUE (1==1)
#define FALSE (1==0)

#ifndef max
#define max(a,b) ( (a) > (b) ? (a) : (b) )
#endif
#ifndef min
#define min(a,b) ( (a) < (b) ? (a) : (b) )
#endif

extern int vrbflg;

#define FILE_FSTATS 1

/*----------------------------------------------------------------------*/
/* empty structs for initialziations */

static const DopplerScanState empty_DopplerScanState;
static const PulsarTimesParamStruc empty_PulsarTimesParamStruc;

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

  INT4 *maxIndex=NULL; 			/*  array that contains indexes of maximum of each cluster */
  CHAR Fstatsfilename[256]; 		/* Fstats file name*/
  UINT4 spdwn;				/* counter over spindown-params */
  DopplerScanInit scanInit;		/* init-structure for DopperScanner */
  DopplerScanState thisScan = empty_DopplerScanState; /* current state of the Doppler-scan */
  DopplerPosition dopplerpos;		/* current search-parameters */
  SkyPosition thisPoint;
  FILE *fpOut=NULL;
  UINT4 loopcounter;
  FStatisticVector *FVect = NULL;   /* new type to store F-statistic results in a frequency-band */
  REAL8Vector *DeltaT = NULL;	/* SSB times at SFT half-points: DeltaT_alpha */
  REAL8Vector *Tdot = NULL;	/* derivatives of SSB-times at SFT half-points: Tdot_alpha */

  UINT4 nBins; 			/* number of frequency-bins */
  UINT4 i;			/* loop index over frequency-bins */


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

  /* main initialization of the code: */
  LAL_CALL ( InitFStat(&status, &GV), &status);

  /* allocate memory of the SSB-times: DeltaT_alpha and Tdot_alpha */
  LAL_CALL ( LALDCreateVector (&status, &DeltaT, GV.SFTno), &status );
  LAL_CALL ( LALDCreateVector (&status, &Tdot, GV.SFTno), &status );


  /* normalize SFTs by running median */
  LAL_CALL (NormaliseSFTDataRngMdn(&status), &status);

  /*      open file */
  strcpy(Fstatsfilename,"Fstats");
  if ( uvar_outputLabel )
    strcat(Fstatsfilename, uvar_outputLabel);

  if ( FILE_FSTATS && !(fpstat=fopen(Fstatsfilename,"w")))
    {
      fprintf(stderr,"in Main: unable to open Fstats file\n");
      return 2;
    }

  /* prepare initialization of DopplerScanner to step through paramter space */
  scanInit.dAlpha = uvar_dAlpha;
  scanInit.dDelta = uvar_dDelta;
  scanInit.gridType = uvar_gridType;
  scanInit.metricType = uvar_metricType;
  scanInit.metricMismatch = uvar_metricMismatch;
  scanInit.obsBegin = GV.Tstart;
  scanInit.obsDuration = GV.tObs;
  scanInit.fmax  = uvar_Freq + uvar_FreqBand;
  scanInit.Detector = &GV.Detector;
  scanInit.ephemeris = GV.edat;		/* used by Ephemeris-based metric */
  scanInit.searchRegion.skyRegionString = GV.skyRegionString;
  scanInit.skyGridFile = uvar_skyGridFile;
  
  if (lalDebugLevel) LALPrintError ("\nSetting up template grid ...");
  
  LAL_CALL ( InitDopplerScan ( &status, &thisScan, &scanInit), &status); 

  /*----------------------------------------------------------------------*/
  if (lalDebugLevel) LALPrintError ("done.\n");
  if ( uvar_outputSkyGrid ) {
    LALPrintError ("\nNow writing sky-grid into file '%s' ...", uvar_outputSkyGrid);
    LAL_CALL (writeSkyGridFile( &status, thisScan.grid, uvar_outputSkyGrid, &scanInit), &status);
    LALPrintError (" done.\n\n");
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

  /* prepare memory to hold F-statistic array over frequency (for cluster-stuff) */
  nBins = GV.FreqImax;
  if ( (FVect = (FStatisticVector*)LALCalloc(1, sizeof(FStatisticVector))) == NULL) {
    LALPrintError ("\nOut of memory..!\n\n");
    return (COMPUTEFSTATC_EMEM);
  }
  LAL_CALL ( LALDCreateVector (&status, &(FVect->F), nBins), &status);

  FVect->f0 = uvar_Freq;
  FVect->df = GV.dFreq;
  FVect->fBand = (nBins - 1) * FVect->df;
  FVect->length = nBins;


  if (lalDebugLevel) LALPrintError ("\nStarting main search-loop.. \n");

  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
   */
  loopcounter = 0;
  while (1)
    {
      EphemerisData *ephem = NULL;
      LAL_CALL (NextDopplerPos( &status, &dopplerpos, &thisScan ), &status);
      /* Have we scanned all DopplerPositions yet? */
      if (thisScan.state == STATE_FINISHED)
	break;

      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
      thisPoint.longitude = dopplerpos.Alpha;
      thisPoint.latitude = dopplerpos.Delta;
      thisPoint.system = COORDINATESYSTEM_EQUATORIAL;
      LAL_CALL (LALNormalizeSkyPosition(&status, &thisPoint, &thisPoint), &status);

      /*----- calculate SSB-times DeltaT_alpha and Tdot_alpha for this skyposition */
      if ( uvar_useEphemeris )
	ephem = GV.edat;
      else
	ephem = NULL;
      LAL_CALL ( LALGetSSBtimes (&status, DeltaT, Tdot, GV.Tstart, &(GV.midTS), 
				 &thisPoint, &(GV.Detector), ephem ), &status);
      GV.CFSparams->DeltaT = DeltaT;
      GV.CFSparams->Tdot = Tdot;

      if ( lalDebugLevel ) 
	{
	  UINT4 N, j;
	  N = GV.timestamps.length;
	  printf ( "\nDEBUG: result from LALGetSSBtimes(): \n");
	  for ( j=0; j < N; j++ )
	    {
	      printf ("i=%d: t_i = %d s, deltaTau_i = %22.20g, dTau_i = %.20g\n",
		      j, GV.timestamps.data[j].gpsSeconds,
		      DeltaT->data[j], Tdot->data[j] );
	    }
	  
	} /* if lalDebugLevel */

      LAL_CALL (CreateDemodParams(&status, GV.CFSparams, &thisPoint, &GV), &status);
      
      /* loop over first-order spindown params */
      for (spdwn=0; spdwn <= GV.SpinImax; spdwn++)
	{
	  if ( GV.spdnOrder > 0 )
	    GV.CFSparams->fkdot->data[1] = uvar_f1dot + spdwn * GV.df1dot;

	  /* Loop over frequencies to be demodulated */
	  for(i=0 ; i< nBins ; i++ )
	    {
	      GV.CFSparams->fkdot->data[0] = uvar_Freq + i * GV.dFreq;	

	      if ( XLALcomputeFStat(&(FVect->F->data[i]), SFTvect, GV.CFSparams) )
		{
		  int code = xlalErrno;
		  XLALClearErrno(); 
		  LALPrintError ("\nERROR: XLALcomputeFStat() failed (xlalErrno = %d)!\n\n", code);
		  return (COMPUTEFSTATC_EXLAL);
		}

	    } /* for i < nBins: loop over frequency-bins */

	  
	  /* now, if user requested it, we output ALL F-statistic results above threshold */
	  if ( fpOut )
	    {
	      UINT4 j;
	      for(j=0; j < FVect->length; j++)
		{
		  REAL8 FStat = 2.0 * medianbias * FVect->F->data[j];
		  REAL8 freq = uvar_Freq + j * GV.dFreq;
		  REAL8 f1dot;
		  
		  if ( GV.spdnOrder > 0 )
		    f1dot = GV.CFSparams->fkdot->data[1];
		  else
		    f1dot = 0;
		  
		  if ( FStat > uvar_Fthreshold )
		    fprintf (fpOut, "%16.12f %8.7f %8.7f %.17g %10.6g\n", 
			     freq, dopplerpos.Alpha, dopplerpos.Delta, 
			     f1dot, FStat);
		}
	      
	    } /* if outputFstat */





	  /* FIXME: to keep cluster-stuff working, we provide the "translation" 
	   * from FVect back into old Fstats-struct */
	  Fstat.F  = FVect->F->data;
	  
	  /*  This fills-in highFLines */
	  if (GV.FreqImax > 5) {
	    LAL_CALL (EstimateFLines(&status, FVect), &status);
	  }
	  
	  /*  This fills-in highFLines  */
	  if (highFLines != NULL && highFLines->Nclusters > 0)
	    {
	      maxIndex=(INT4 *)LALMalloc(highFLines->Nclusters*sizeof(INT4));
	      
	      /*  for every cluster writes the information about it in file Fstats */
	      if (writeFLines(maxIndex, FVect->f0, FVect->df)){
		fprintf(stderr, "%s: trouble making file Fstats\n", argv[0]);
		return 6;
	      }
	  
	      LALFree(maxIndex);
	    } /* if highFLines found */

	  
	  /* Set the number of the clusters detected to 0 at each iteration 
	   * of the sky-direction and the spin down */
	  highFLines->Nclusters=0;

	} /* For GV.spinImax */

      loopcounter ++;
      if (lalDebugLevel) 
	LALPrintError ("\
Search progress: %5.1f%%", (100.0* loopcounter / thisScan.numGridPoints));

    } /*  while SkyPos */

  if (uvar_outputFstat && fpOut)
    fclose (fpOut);

  if (lalDebugLevel) LALPrintError ("\nSearch finished.\n");

  /* properly terminate Fstats-file by 'DONE' marker: */ 
  if (fpstat) fprintf(fpstat, "%%DONE\n");
  if (fpstat) fclose(fpstat);

  /* Free memory */
  LAL_CALL ( FreeDopplerScan(&status, &thisScan), &status);

  LAL_CALL ( Freemem(&status, &GV), &status);
  
  LAL_CALL ( LALDDestroyVector (&status, &DeltaT), &status);
  LAL_CALL ( LALDDestroyVector (&status, &Tdot), &status);

  LAL_CALL (LALDDestroyVector (&status, &(FVect->F)), &status);
  LALFree (FVect);

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
  uvar_FstatMin = 0.0;

  uvar_skyGridFile = NULL;

  uvar_workingDir = LALMalloc(512);
  strcpy(uvar_workingDir, ".");

  uvar_dopplermax =  1.05e-4;
  uvar_windowsize = 50;	/* for running-median */

  uvar_useEphemeris = TRUE;

  /* register all our user-variables */
  LALregBOOLUserVar(status, 	help, 		'h', UVAR_HELP,     "Print this message"); 
  LALregINTUserVar(status,	Dterms,		't', UVAR_OPTIONAL, "Number of terms to keep in Dirichlet kernel sum");
  LALregREALUserVar(status, 	Freq, 		'f', UVAR_REQUIRED, "Starting search frequency in Hz");
  LALregREALUserVar(status, 	FreqBand, 	'b', UVAR_OPTIONAL, "Search frequency band in Hz");
  LALregREALUserVar(status,     dFreq,          'r', UVAR_OPTIONAL, "Frequency resolution in Hz (default: 1/(2*Tsft*Nsft)");
  LALregREALUserVar(status, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	AlphaBand, 	'z', UVAR_OPTIONAL, "Band in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	DeltaBand, 	'c', UVAR_OPTIONAL, "Band in delta (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	dAlpha, 	'l', UVAR_OPTIONAL, "Resolution in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status, 	dDelta, 	'g', UVAR_OPTIONAL, "Resolution in delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(status,	DataDir, 	'D', UVAR_OPTIONAL, "Directory where SFT's are located");
  /*================*/
  LALregSTRINGUserVar(status,	DataDir2, 	 0, UVAR_OPTIONAL, "Directory where SFT's are located");
  /*================*/
  LALregSTRINGUserVar(status,	mergedSFTFile, 	'B', UVAR_OPTIONAL, "Merged SFT's file to be used"); 
  LALregSTRINGUserVar(status,	BaseName, 	'i', UVAR_OPTIONAL, "The base name of the input  file you want to read");
  /*================*/
  LALregSTRINGUserVar(status,	BaseName2, 	 0, UVAR_OPTIONAL, "The base name of the input  file you want to read");  
  /*================*/
  LALregSTRINGUserVar(status,	ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(status,	ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregSTRINGUserVar(status, 	IFO, 		'I', UVAR_REQUIRED, "Detector: GEO(0), LLO(1), LHO(2), NAUTILUS(3), VIRGO(4), TAMA(5), CIT(6)");
  /*================*/
  LALregSTRINGUserVar(status, 	IFO2, 		 0,  UVAR_OPTIONAL, "Detector: GEO(0), LLO(1), LHO(2), NAUTILUS(3), VIRGO(4), TAMA(5), CIT(6)"); 
  /*================*/
  LALregBOOLUserVar(status, 	SignalOnly, 	'S', UVAR_OPTIONAL, "Signal only flag");
  LALregREALUserVar(status, 	dopplermax, 	'q', UVAR_OPTIONAL, "Maximum doppler shift expected");  
  LALregREALUserVar(status, 	f1dot, 		's', UVAR_OPTIONAL, "First spindown parameter f1dot");
  LALregREALUserVar(status, 	f1dotBand, 	'm', UVAR_OPTIONAL, "Search-band for f1dot");
  LALregREALUserVar(status, 	df1dot, 	'e', UVAR_OPTIONAL, "Resolution for f1dot (default 1/(2*Tobs*tSFT*Nsft)");
  LALregREALUserVar(status, 	Fthreshold,	'F', UVAR_OPTIONAL, "Signal Set the threshold for selection of 2F");
  LALregINTUserVar(status, 	windowsize,	'k', UVAR_OPTIONAL, "Running-Median window size");
  LALregINTUserVar(status, 	gridType,	 0 , UVAR_OPTIONAL, "Template grid: 0=flat, 1=isotropic, 2=metric, 3=file");
  LALregINTUserVar(status, 	metricType,	'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserVar(status, 	metricMismatch,	'X', UVAR_OPTIONAL, "Maximal mismatch for metric tiling");
  LALregSTRINGUserVar(status,	skyRegion, 	'R', UVAR_OPTIONAL, "Specify sky-region by polygon");
  LALregSTRINGUserVar(status,	outputLabel,	'o', UVAR_OPTIONAL, "Label to be appended to all output file-names");
  LALregSTRINGUserVar(status,	skyGridFile,	 0,  UVAR_OPTIONAL, "Load sky-grid from this file.");
  LALregSTRINGUserVar(status,	outputSkyGrid,	 0,  UVAR_OPTIONAL, "Write sky-grid into this file.");
  LALregSTRINGUserVar(status,     workingDir,     'w', UVAR_OPTIONAL, "Directory to be made the working directory, . is default");
  /* more experimental and unofficial stuff follows here */
  LALregSTRINGUserVar(status,	outputFstat,	 0,  UVAR_OPTIONAL, "Output the F-statistic field over the parameter-space");
  LALregREALUserVar(status, 	FstatMin,	 0,  UVAR_OPTIONAL, "Minimum F-Stat value to written into outputFstat-file");
  LALregBOOLUserVar (status,	useEphemeris,    0,  UVAR_DEVELOPER, "Use ephemeris or Ptolemaic model");

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* initUserVars() */



void 
CreateDemodParams (LALStatus *status, 
		   computeFStatPar *CFSparams,
		   const SkyPosition *pos,
		   const ConfigVariables *cfg)
{
  EarthState earth;
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  AMCoeffsParams amParams;
  LALDetAndSource das;
  LALSource pSource;

  INITSTATUS (status, "CreateDemodParams", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT ( CFSparams, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( pos, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( pos->system == COORDINATESYSTEM_EQUATORIAL, status, 
	   COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( cfg, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  
  /* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.site.location[0] = cfg->Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1] = cfg->Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2] = cfg->Detector.location[2]/LAL_C_SI;
  baryinput.alpha = pos->longitude;
  baryinput.delta = pos->latitude;
  baryinput.dInv=0.e0;

  /* Fill up AMCoeffsParams structure */
  amParams.das = &das;
  amParams.das->pSource = &pSource;
  amParams.baryinput = &baryinput;
  amParams.earth = &earth; 
  amParams.edat = cfg->edat;
  amParams.das->pDetector = &(cfg->Detector); 
  amParams.das->pSource->equatorialCoords = *pos;
  amParams.das->pSource->orientation = 0.0;

  amParams.polAngle = amParams.das->pSource->orientation ; /* These two have to be the same!!!!!!!!!*/
  amParams.leapAcc = LALLEAPSEC_STRICT;

  TRY ( LALComputeAM(status->statusPtr, CFSparams->amcoe, cfg->midTS.data, &amParams), status); 

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* CreateDemodParams() */

/*******************************************************************************/
/*  for every cluster writes the information about it in file Fstats */
/*  precisely it writes: */
/*  fr_max alpha delta N_points_of_cluster mean std max (of 2F) */
int writeFLines(INT4 *maxIndex, REAL8 f0, REAL8 df){

  INT4 i,j,j1,j2,k,N;
  REAL8 max,logof2,mean,var,std,R,fr;
  INT4 imax;
  INT4 err = 0;

  logof2=medianbias;
 
  j1=0;
  j2=0;

  for (i=0;i<highFLines->Nclusters;i++){
    N=highFLines->NclustPoints[i];
    
    /*  determine maximum of j-th cluster */
    /*  and compute mean */
    max=0.0;
    imax=0;
    mean=0.0;
    std=0.0;
    for (j=0;j<N;j++){
      R=2.0*logof2*highFLines->clusters[j1];
      k=highFLines->Iclust[j1];
      j1=j1+1;
      mean=mean+R;
      if( R > max){
	max=R;
	imax=k;
      }
    }/*  end j loop over points of i-th cluster  */
    /*  and start again to compute variance */
    maxIndex[i]=imax;
    mean=mean/N;
    var=0.0;
    for (j=0;j<N;j++){
      R=2.0*logof2*highFLines->clusters[j2];
      j2=j2+1;
      var=var+(R-mean)*(R-mean);
    }/*  end j loop over points of i-th cluster  */
    var=var/N;
    std=sqrt(var);
    fr = f0 + imax * df;
    /*    print the output */
    if (fpstat)
      err=fprintf(fpstat,"%16.12f %10.8f %10.8f    %d %10.5f %10.5f %10.5f\n",fr,
		Alpha, Delta, N, mean, std, max);

    if (err<=0) {
    fprintf(stderr,"writeFLines couldn't print to Fstas!\n");
    return 4;
  }

  }/*  end i loop over different clusters */

  return 0;
}

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
  UINT4 i;

  INITSTATUS (status, "InitFStat", rcsid);
  ATTATCHSTATUSPTR (status);

  /* ----------------------------------------------------------------------
   * do some sanity checks on the user-input before we proceed 
   */
  if(!uvar_DataDir && !uvar_mergedSFTFile)
    {
      LALPrintError ( "\nMust specify 'DataDir' OR 'mergedSFTFile'\n"
		      "No SFT directory specified; input directory with -D option.\n"
		      "No merged SFT file specified; input file with -B option.\n"
		      "Try ./ComputeFStatistic_v2 -h \n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);      
    }

  if(uvar_DataDir && uvar_mergedSFTFile)
    {
      LALPrintError ( "\nCannot specify both 'DataDir' and 'mergedSFTfile'.\n"
		      "Try ./ComputeFStatistic_v2 -h \n\n" );
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }      

   if (uvar_ephemYear == NULL)
    {
      LALPrintError ("\nNo ephemeris year specified (option 'ephemYear')\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }      
  /* don't allow negative frequency-band for safety */
  if ( uvar_FreqBand < 0)
    {
      LALPrintError ("\nNegative value of frequency-band not allowed !\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }

  /* don't allow negative bands (for safty in griding-routines) */
  if ( (uvar_AlphaBand < 0) ||  (uvar_DeltaBand < 0) )
    {
      LALPrintError ("\nNegative value of sky-bands not allowed (alpha or delta)!\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }

  /* set the current working directory */
  if(chdir(uvar_workingDir) != 0)
    {
      fprintf(stderr, "in Main: unable to change directory to %s\n", uvar_workingDir);
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }

  /*----------------------------------------------------------------------
   * load SFT data-files
   */
  /* min and max frequency index to read */
  {
    REAL8 f0, f1;
    CHAR *fpattern = NULL;
    UINT4 len = 0;

    f0 = (1.0 - uvar_dopplermax) * uvar_Freq;	/* lower physical frequency bound */
    f1 = (1.0 + uvar_dopplermax) * (uvar_Freq + uvar_FreqBand); /* upper physical bound */

    if (uvar_DataDir) 
      len = strlen(uvar_DataDir);
    else
      len = 1;	/* '.' */
 
    if (uvar_BaseName) 
      len += strlen(uvar_BaseName);

    len += 10;	/* to allow for '/' and '*' to be inserted */

    if ( (fpattern = LALCalloc(1, len)) == NULL)
      {
	LALPrintError ("\nOut of memory !\n");
	ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
      }

    if (uvar_DataDir)
      strcpy (fpattern, uvar_DataDir);
    else
      strcpy (fpattern, ".");

    strcat (fpattern, "/*");
    
    if (uvar_BaseName)
      {
	strcat (fpattern, uvar_BaseName);
	strcat (fpattern, "*");
      }
      
    TRY ( LALReadSFTfiles(status->statusPtr, &SFTvect, f0, f1, uvar_Dterms, fpattern), status);
    LALFree (fpattern);

  } /* SFT-loading */

  /* deduce search-parameters determined by input-data */
  cfg->tSFT = 1.0 / SFTvect->data[0].deltaF;
  cfg->SFTno = SFTvect->length;
  cfg->nsamples = SFTvect->data[0].data->length;

  /*----------------------------------------------------------------------*/
  /* extract timestamps-vector from SFT-data */
  /* and prepare SFT-data in FFT-struct to feed it into LALDemod() (FIXME)*/
  cfg->timestamps.length = cfg->SFTno;
  if ( (cfg->timestamps.data =  (LIGOTimeGPS *) LALCalloc(cfg->SFTno, sizeof(LIGOTimeGPS))) == NULL) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  cfg->midTS.length = cfg->SFTno;
  if ( (cfg->midTS.data = (LIGOTimeGPS *) LALCalloc(cfg->SFTno, sizeof(LIGOTimeGPS))) == NULL) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  for (i=0; i < cfg->SFTno; i++)
    {
      cfg->timestamps.data[i] = SFTvect->data[i].epoch;
      /* to get midpoints, simply add Tsft/2 to each timestamp */
      TRY (LALAddFloatToGPS (status->statusPtr, &(cfg->midTS.data[i]), &(cfg->timestamps.data[i]), 
			     0.5*cfg->tSFT), status);
    }

  /* figure out total observation time */
  {
    LIGOTimeGPS t0, t1;
    t0 = SFTvect->data[0].epoch;
    t1 = SFTvect->data[cfg->SFTno-1].epoch;
    TRY (LALDeltaFloatGPS (status->statusPtr, &(cfg->tObs), &t1, &t0), status);
    cfg->tObs += cfg->tSFT;
    cfg->Tstart = t0;
  }

  /*----------------------------------------------------------------------
   * set up and check ephemeris-file locations
   */

  if (LALUserVarWasSet (&uvar_ephemDir) )
    {
      sprintf(cfg->EphemEarth, "%s/earth%s.dat", uvar_ephemDir, uvar_ephemYear);
      sprintf(cfg->EphemSun, "%s/sun%s.dat", uvar_ephemDir, uvar_ephemYear);
    }
  else
    {
      sprintf(cfg->EphemEarth, "earth%s.dat", uvar_ephemYear);
      sprintf(cfg->EphemSun, "sun%s.dat",  uvar_ephemYear);
    }

  /*----------------------------------------------------------------------
   * initialize detector 
   */
  if ( !strcmp (uvar_IFO, "GEO") || !strcmp (uvar_IFO, "0") ) 
    cfg->Detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  else if ( !strcmp (uvar_IFO, "LLO") || ! strcmp (uvar_IFO, "1") ) 
    cfg->Detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if ( !strcmp (uvar_IFO, "LHO") || !strcmp (uvar_IFO, "2") )
    cfg->Detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else if ( !strcmp (uvar_IFO, "NAUTILUS") || !strcmp (uvar_IFO, "3"))
    {
      TRY (CreateNautilusDetector (status->statusPtr, &(cfg->Detector)), status);
    }
  else if ( !strcmp (uvar_IFO, "VIRGO") || !strcmp (uvar_IFO, "4") )
    cfg->Detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  else if ( !strcmp (uvar_IFO, "TAMA") || !strcmp (uvar_IFO, "5") )
    cfg->Detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  else if ( !strcmp (uvar_IFO, "CIT") || !strcmp (uvar_IFO, "6") )
    cfg->Detector = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
  else
    {
      LALPrintError ("\nUnknown detector. Currently allowed are \
'GEO', 'LLO', 'LHO', 'NAUTILUS', 'VIRGO', 'TAMA', 'CIT' or '0'-'6'\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }
    
  /*----------------------------------------------------------------------
   * set some defaults
   */
  /* 'natural' FFT frequency-resolution */
  cfg->dFreq = 1.0 / (2.0 * cfg->tObs);
  
  /*Number of frequency values to calculate F for */
  cfg->FreqImax = (INT4)(uvar_FreqBand/cfg->dFreq + 0.5) + 1;  


  if (LALUserVarWasSet (&uvar_f1dotBand) && (uvar_f1dotBand != 0) )
    {
      cfg->SpinImax=(int)(uvar_f1dotBand/uvar_df1dot+.5) + 1; 
      cfg->spdnOrder = 1;
    }
  else
    {
      cfg->SpinImax = 0;
      cfg->spdnOrder = 0;
    }

  /*----------------------------------------------------------------------
   * initialize+check  template-grid related parameters 
   */
  {
    BOOLEAN haveSkyRegion, haveAlphaDelta, haveGridFile, useGridFile, haveMetric, useMetric;

    haveSkyRegion  = (uvar_skyRegion != NULL);
    haveAlphaDelta = (LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta) );
    haveGridFile   = (uvar_skyGridFile != NULL);
    useGridFile   = (uvar_gridType == GRID_FILE);
    haveMetric     = (uvar_metricType > LAL_PMETRIC_NONE);
    useMetric     = (uvar_gridType == GRID_METRIC);

    /* some consistency checks on input to help catch errors */
    if ( !useGridFile && !(haveSkyRegion || haveAlphaDelta) )
      {
	LALPrintError ("\nNeed sky-region: either use (Alpha,Delta) or skyRegion!\n\n");
	ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
      }
    if ( haveSkyRegion && haveAlphaDelta )
      {
	LALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta) OR skyRegion!\n\n");
	ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
      }
    if ( useGridFile && !haveGridFile )
      {
	LALPrintError ("\nERROR: gridType=FILE, but no skyGridFile specified!\n\n");
	ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);	
      }
    if ( !useGridFile && haveGridFile )
      {
	LALWarning (status, "\nWARNING: skyGridFile was specified but not needed ... will be ignored\n");
      }
    if ( useGridFile && (haveSkyRegion || haveAlphaDelta) )
      {
	LALWarning (status, "\nWARNING: Using skyGridFile, but sky-region was specified ... ignored!\n");
      }
    if ( !useMetric && haveMetric) 
      {
	LALWarning (status, "\nWARNING: Metric was specified for non-metric grid... will be ignored!\n");
      }
    if ( useMetric && !haveMetric) 
      {
	LALPrintError ("\nERROR: metric grid-type selected, but no metricType selected\n\n");
	ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);      
      }

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

  /* ----------------------------------------------------------------------*/
  /*
   * initialize Ephemeris-data 
   */
  {
    LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
    INT4 leap;

    cfg->edat = LALCalloc(1, sizeof(EphemerisData));
    cfg->edat->ephiles.earthEphemeris = cfg->EphemEarth;
    cfg->edat->ephiles.sunEphemeris = cfg->EphemSun;

    TRY (LALLeapSecs (status->statusPtr, &leap, &(cfg->Tstart), &formatAndAcc), status);
    cfg->edat->leap = leap;

    TRY (LALInitBarycenter(status->statusPtr, cfg->edat), status);               

  } /* end: init ephemeris data */

  /* ----------------------------------------------------------------------
   * initialize + allocate space for AM-coefficients and Demod-params
   */
  {
    AMCoeffs *amc;
    /* Allocate DemodParams structure */
    if ( (cfg->CFSparams = (computeFStatPar*) LALCalloc(1, sizeof(computeFStatPar))) == NULL) {
      ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }
    /* Allocate space for AMCoeffs */
    if ( (amc = LALCalloc(1, sizeof(AMCoeffs))) == NULL) {
      ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }
    TRY (LALSCreateVector(status->statusPtr, &(amc->a), (UINT4) cfg->SFTno), status);
    TRY (LALSCreateVector(status->statusPtr, &(amc->b), (UINT4) cfg->SFTno), status);
    
    cfg->CFSparams->amcoe = amc;
    cfg->CFSparams->tSFT = cfg->tSFT;
    cfg->CFSparams->Dterms = uvar_Dterms;
    /* prepare memory for fkdot - vector : (f, f1dot, f2dot, ..) */
    TRY ( LALDCreateVector (status->statusPtr, &(cfg->CFSparams->fkdot), cfg->spdnOrder + 1), status );

  } /* end: init AM- and demod-params */

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
      ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }
    strcpy (fname, head);
    if (uvar_outputLabel)
      strcat (fname, uvar_outputLabel);
    strcat (fname, ".log");

    if ( (fplog = fopen(fname, "w" )) == NULL) {
      LALPrintError ("\nFailed to open log-file '%f' for writing.\n\n", fname);
      LALFree (fname);
      ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
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

  INITSTATUS (status, "Freemem", rcsid);
  ATTATCHSTATUSPTR (status);

  /* Free SFT data */
  TRY (LALDestroySFTVector (status->statusPtr, &SFTvect), status);	 /* the new way*/

  /* Free timestamps */
  LALFree (cfg->timestamps.data);
  LALFree (cfg->midTS.data);

  /* FIXME: quite a few missing here */

  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (status->statusPtr), status);

  LALFree ( cfg->skyRegionString );

  /* this comes from clusters.c */
  if (highFLines->clusters) LALFree(highFLines->clusters);
  if (highFLines->Iclust) LALFree(highFLines->Iclust);
  if (highFLines->NclustPoints) LALFree(highFLines->NclustPoints);


  /* Free ephemeris data */
  LALFree(cfg->edat->ephemE);
  LALFree(cfg->edat->ephemS);
  LALFree(cfg->edat);


  TRY ( LALDDestroyVector (status->statusPtr, &(cfg->CFSparams->fkdot)), status);

  TRY (LALSDestroyVector(status->statusPtr, &(cfg->CFSparams->amcoe->a)), status);
  TRY (LALSDestroyVector(status->statusPtr, &(cfg->CFSparams->amcoe->b)), status);
  LALFree ( cfg->CFSparams->amcoe);

  LALFree ( cfg->CFSparams);

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* Freemem() */

/*******************************************************************************/
/** Sorting function to sort into DECREASING order. Used in PrintTopValues(). */
int compare(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=Fstat.F[*(const int *)ip];
  dj=Fstat.F[*(const int *)jp];

  if (di<dj)
    return 1;
  
  if (di==dj)
    return 0;

  return -1;
}

/** Find outliers and then clusters in the F-statistic array over frequency. 
 * These clusters get written in the global highFLines. 
 */
void
EstimateFLines(LALStatus *status, const FStatisticVector *FVect)
{
  UINT4 i,j,Ntot;   
  UINT4 nbins;                	/**< Number of bins in F */
  REAL8Vector *F1=NULL; 
  REAL8Vector *FloorF1=NULL;             /* Square of SFT */
  REAL4 THR=10.0;
  REAL8 dFreq, f0;
  
  OutliersInput  *outliersInput;
  OutliersParams *outliersParams;
  Outliers       *outliers;
  ClustersInput  *clustersInput;
  ClustersParams *SpClParams;
  Clusters       *SpLines=highFLines;
    
  INT2 smallBlock=1;
  INT4 wings;

  INITSTATUS( status, "EstimateFLines", rcsid );
  ATTATCHSTATUSPTR (status);

  nbins = FVect->length;
  dFreq = FVect->df;
  f0 = FVect->f0;

  THR=uvar_Fthreshold;

/* 0.0002 is the max expected width of the F status curve for signal */
/* with ~ 10 h observation time */

  wings = (UINT4) (0.5 + 0.0002 / dFreq );

  TRY ( LALDCreateVector(status->statusPtr, &F1, nbins), status);
  TRY ( LALDCreateVector(status->statusPtr, &FloorF1, nbins), status);

  /* loop over SFT data to estimate noise */
  for (j=0;j<nbins;j++)
    {
      F1->data[j] = FVect->F->data[j];
      FloorF1->data[j] = 1.0;
    }
  
  F1->length = nbins;
  FloorF1->length = nbins;

  if ( (outliers = (Outliers *)LALCalloc(1, sizeof(Outliers))) == NULL) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  outliers->Noutliers=0;

  if ( (outliersParams = (OutliersParams *)LALCalloc(1,sizeof(OutliersParams))) == NULL) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  if ( (outliersInput = (OutliersInput *)LALCalloc(1,sizeof(OutliersInput))) == NULL) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  
  outliersParams->Thr = THR/(2.0*medianbias);
  outliersParams->Floor = FloorF1;
  outliersParams->wings = wings; /*these must be the same as ClustersParams->wings */
  outliersInput->ifmin = (INT4) ((f0 / dFreq) + 0.5);
  outliersInput->data = F1;

  /*find values of F above THR and populate outliers with them */
  ComputeOutliers(outliersInput, outliersParams, outliers);


  /*if no outliers were found clean and exit */
   if (outliers->Noutliers == 0){

     LALFree(outliers->ratio);
     LALFree(outliers);
     LALFree(outliersParams);
     LALFree(outliersInput);
     TRY ( LALDDestroyVector(status->statusPtr, &F1), status);
     TRY ( LALDDestroyVector(status->statusPtr, &FloorF1), status);

     /*      fprintf(stderr,"Nclusters zero \n"); */
     /*      fflush(stderr); */

     goto finished;

   } /* if Noutliers == 0 */
  

   /* if outliers are found get ready to identify clusters of outliers*/
   if ( (SpClParams = (ClustersParams*) LALCalloc(1,sizeof(ClustersParams))) == NULL) {
     ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
   }
   
   if ( (clustersInput = (ClustersInput *) LALCalloc(1,sizeof(ClustersInput))) == NULL) {
     ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
   }
      
   SpClParams->wings = wings;
   SpClParams->smallBlock = smallBlock;
   
   clustersInput->outliersInput = outliersInput;
   clustersInput->outliersParams= outliersParams;
   clustersInput->outliers      = outliers;     
   
   /* clusters of outliers in F get written in SpLines which is the global highFLines*/
   TRY (DetectClusters(status->statusPtr, clustersInput, SpClParams, SpLines), status);
   
   /*  sum of points in all lines */
   Ntot=0;
   for (i=0; i < (UINT4)SpLines->Nclusters; i++){ 
     Ntot = Ntot + SpLines->NclustPoints[i];
   }

   TRY ( LALDDestroyVector(status->statusPtr, &F1), status);
   TRY ( LALDDestroyVector(status->statusPtr, &FloorF1), status);

   LALFree(outliers->ratio);
   LALFree(outliers->outlierIndexes);
   LALFree(outliers);
   LALFree(outliersParams);
   LALFree(outliersInput);
   LALFree(SpClParams);
   LALFree(clustersInput);

 finished:
   DETATCHSTATUSPTR(status);
   RETURN(status);

} /* EstimateFLines() */


/** Normalise the SFT-array \em SFTData by the running median.
 * The running median windowSize in this routine determines 
 * the sample bias which, instead of log(2.0), must be 
 * multiplied by F statistics.
 */
void 
NormaliseSFTDataRngMdn(LALStatus *status)
{
  INT4 m, il;                         /* loop indices */
  UINT4 i, j, lpc;
  UINT4 Ntot;
  REAL8Vector *Sp=NULL, *RngMdnSp=NULL;   /* |SFT|^2 and its rngmdn  */
  REAL8 B;                          /* SFT Bandwidth */
  REAL8 deltaT,norm,*N, *Sp1;
  INT2 windowSize=uvar_windowsize;                  /* Running Median Window Size*/
  REAL4 xre,xim,xreNorm,ximNorm;

  INITSTATUS( status, "NormaliseSFTDataRngMdn", rcsid );
  ATTATCHSTATUSPTR (status);

  if ( !uvar_SignalOnly ) {
    TRY( LALRngMedBias (status->statusPtr, &medianbias, windowSize), status);
  }

  TRY ( LALDCreateVector(status->statusPtr, &Sp, GV.nsamples), status);
  TRY ( LALDCreateVector(status->statusPtr, &RngMdnSp, GV.nsamples), status);

  if( (N = (REAL8 *) LALCalloc(GV.nsamples,sizeof(REAL8))) == NULL) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  if( (Sp1 = (REAL8 *) LALCalloc(GV.nsamples,sizeof(REAL8))) == NULL) { 
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  /* loop over each SFTs */
  for (i=0;i<GV.SFTno;i++)         
    {
      /* Set to zero the values */
      for (j=0;j<GV.nsamples;j++){
	RngMdnSp->data[j] = 0.0;
	Sp->data[j]       = 0.0;
      }
      
      /* loop over SFT data to estimate noise */
      for (j=0;j<GV.nsamples;j++){
	xre=SFTvect->data[i].data->data[j].re;
	xim=SFTvect->data[i].data->data[j].im;
	Sp->data[j]=(REAL8)(xre*xre+xim*xim);
      }
      
      /* Compute running median */
      TRY ( EstimateFloor(status->statusPtr, Sp, windowSize, RngMdnSp), status);

      /* compute how many cluster points in all */
      /* substitute the line profiles value in RngMdnSp */
      Ntot=0;
      if (highSpLines != NULL){
	for (il=0;il<highSpLines->Nclusters;il++){
	  Ntot=Ntot+highSpLines->NclustPoints[il];
	}
	for (j=0;j<Ntot;j++){
	  m=highSpLines->Iclust[j];
	  RngMdnSp->data[m]=RngMdnSp->data[m]*highSpLines->clusters[j];	
	}
      }

      /*Compute Normalization factor*/
      /* for signal only case as well */  
      for (lpc=0;lpc<GV.nsamples;lpc++){
	N[lpc]=1.0/sqrt(2.0*RngMdnSp->data[lpc]);
      }
      
      if(uvar_SignalOnly == 1){
	B=(1.0*GV.nsamples)/(1.0*GV.tSFT);
	deltaT=1.0/(2.0*B);
	norm=deltaT/sqrt(GV.tSFT);
	for (lpc=0;lpc<GV.nsamples;lpc++){
	  N[lpc]=norm;
	}
      }
      
      /*  loop over SFT data to normalise it (with N) */
      /*  also compute Sp1, average normalized PSD */
      /*  and the sum of the PSD in the band, SpSum */
      for (j=0;j<GV.nsamples;j++){
	xre=SFTvect->data[i].data->data[j].re;
	xim=SFTvect->data[i].data->data[j].im;
	xreNorm=N[j]*xre; 
	ximNorm=N[j]*xim; 
	SFTvect->data[i].data->data[j].re = xreNorm;    
	SFTvect->data[i].data->data[j].im = ximNorm;
	Sp1[j]=Sp1[j]+xreNorm*xreNorm+ximNorm*ximNorm;
      }
      
    } /* end loop over SFTs*/

  LALFree(N);
  LALFree(Sp1);

  TRY ( LALDDestroyVector(status->statusPtr, &RngMdnSp), status);
  TRY ( LALDDestroyVector(status->statusPtr, &Sp), status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
  
} /* NormaliseSFTDataRngMed() */


/** v2-specific version of LALDemod() (based on TestLALDemod() in CFS)
 */
#define OOTWOPI         (1.0 / LAL_TWOPI)
#define LUT_RES         64      /* resolution of lookup-table */

#define TWOPI_FLOAT     6.28318530717958f  /* 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/* 1 / (2pi) */ 

int
XLALNewLALDemod(Fcomponents *FaFb,
		const SFTVector *sfts, 
		const computeFStatPar *params) 
{ 
  UINT4 alpha;                 	/* loop index over SFTs */
  UINT4 spdnOrder;		/* maximal spindown-orders */
  UINT4 numSFTs;		/* number of SFTs (M in the Notes) */
  UINT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  COMPLEX16 Fa, Fb;
#define LUT_RES         64      /* resolution of lookup-table */
  static REAL8 sinVal[LUT_RES+1], cosVal[LUT_RES+1];/*LUT values computed by the routine do_trig_lut*/
  static BOOLEAN firstCall = 1;
  UINT4 index; 
  REAL8 f;

  /* ----- check validity of input */
  if ( !FaFb ) {
    LALPrintError ("\nOutput-pointer is NULL !\n\n");
    XLAL_ERROR ( "XLALNewLALDemod", XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    LALPrintError ("\nInput SFTs are NULL!\n\n");
    XLAL_ERROR ( "XLALNewLALDemod", XLAL_EINVAL);
  }
  
  if ( !params || !params->fkdot ) {
    LALPrintError ("\nIllegal NULL in input !\n\n");
    XLAL_ERROR ( "XLALNewLALDemod", XLAL_EINVAL);
  }

  /* This size LUT gives errors ~ 10^-7 with a three-term Taylor series */
  if ( firstCall )
    {
      UINT4 k;
      for (k=0; k <= LUT_RES; k++)
        {
          sinVal[k] = sin( (LAL_TWOPI*k)/LUT_RES );
          cosVal[k] = cos( (LAL_TWOPI*k)/LUT_RES );
        }
      firstCall = 0;
    }

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  freqIndex0 = (UINT4) ( sfts->data[0].f0 / sfts->data[0].deltaF + 0.5); /* lowest freqency-index */

  spdnOrder = params->fkdot->length - 1;

  f = params->fkdot->data[0];

  Fa.re = 0.0;
  Fa.im = 0.0;
  Fb.re = 0.0;
  Fb.im = 0.0;


  /* Loop over all SFTs  */
  for ( alpha = 0; alpha < numSFTs; alpha++ )
    {
      REAL4 a = params->amcoe->a->data[alpha];
      REAL4 b = params->amcoe->b->data[alpha];

      REAL8 xhat_alpha, y_alpha;	/* xhat(alpha), x(alpha,k) and y(alpha) */
      REAL4 x0;
      REAL8 rem; 		/* remainder of x_alpha0: how close to an int? */

      UINT4 k;			/* loop index over frequency-bins */
      UINT4 kstar;		/* central frequency-bin k* = round(xhat_alpha) */

      COMPLEX8 *Xalpha = sfts->data[alpha].data->data; /* pointer to current SFT-data */
      COMPLEX8 *Xalpha_k; 	/* pointer to frequency-bin k in current SFT */
      REAL4 sinx, cosxm1;	/* sin(x_alpha) and (cos(x_alpha)-1) */
      REAL4 realXP, imagXP;	/* the sum_k X_alpha_k P_alpha_k */
      REAL8 realQ, imagQ;	/* Re and Im of Q = e^{-i y} */
      REAL8 realQXP, imagQXP;	/* Re/Im of Q_alpha XP_alpha */
      UINT4 k0;
      /* ----- calculate x(alpha,0) and y(alpha) */
      {
	UINT4 s; 		/* loop-index over spindown-order */
	REAL8 Tas; 		/* temporary variable to calculate (DeltaT_alpha)^2 */
	UINT4 sfact = 1;	/* store for s! */

	Tas = params->DeltaT->data[alpha]; 	/* DeltaT_alpha = T^1 */

	/* Step 1: s = 0 */
	xhat_alpha = f;		/* f^{0) T^0 / 0! */
	y_alpha = f * Tas;	/* f^{0} T^1 / 1! */

	/* Step 2: sum s >= 1 */
	for (s=1; s <= spdnOrder; s++)
	  {
	    REAL8 fsdot = params->fkdot->data[s];
	    xhat_alpha += fsdot * Tas / sfact; 	/* Tas = T^s here, sfact=s! */
	    Tas *= Tas; 		/* T^(s+1) */
	    sfact *= s;		/* (s+1)! */	  
	    y_alpha += fsdot * Tas / sfact; 
	  } /* for s <= spdnOrder */

	/* Step 3: apply global factors and complete y_alpha */
	xhat_alpha *= params->tSFT * params->Tdot->data[alpha];	/* guaranteed > 0 ! */
	y_alpha -= 0.5 * xhat_alpha;
	
	/* still missing: prefactor 2 pi in y_alpha */
	y_alpha *= LAL_TWOPI;

	/* real- and imaginary part of e^{-i y } */
	realQ = cos(y_alpha);
	imagQ = - sin(y_alpha);
      }
      /* ---------------------------------------- */

      /* xhat_alpha determines the 'central' frequency-bin k* in the sum */
      kstar = (UINT4) (xhat_alpha + 0.5);	/* k* = round(xhat_alpha) */

      /* Trick: sin[ 2pi (xhat - k) ] = sin [ 2pi xhat ], therefore
       * the trig-functions need to be calculated only once!
       * We choose the value sin[ 2pi(xhat - kstar) ] because it is the 
       * smallest and will pose no numerical difficulties !
       */

      /*-------------------- calculate sin(x), cos(x) */
      /* (1) don't use LUT */
      /*
      rem = xhat_alpha - kstar; 	
      x_alpha_k = LAL_TWOPI * rem;	
      sinx = sinf(x_alpha_k);
      cosxm1 = cosf(x_alpha_k) - 1.0f;
      */
      /*  (2) use LUT */
      rem = xhat_alpha - (UINT4)xhat_alpha;	/* positive ! */
      index = (UINT4)( rem * LUT_RES + 0.5 );   /* positive! */
      {
	REAL8 d=LAL_TWOPI*(rem - (REAL8)index/(REAL8)LUT_RES);
	REAL8 d2=0.5*d*d;
	REAL8 ts=sinVal[index];
	REAL8 tc=cosVal[index];
                
	sinx = ts+d*tc-d2*ts;
	cosxm1 = tc-d*ts-d2*tc-1.0;
      }
      /* -------------------- */

      realXP = 0;
      imagXP = 0;

      k0 = kstar - params->Dterms;
      if ( k0 < freqIndex0 ) {
	LALPrintError ("\nLowest frequency-index k0=%d outside of SFT-interval (%d)\n\n",
		       k0, freqIndex0 );
	XLAL_ERROR("XLALNewLALDemod", XLAL_EDOM);
      }

      /* ---------- calculate the (truncated to Dterms) sum over k ---------- */

      /* ---------- ATTENTION: this the "hot-loop", which will be 
       * executed many millions of times, so anything in here 
       * has a HUGE impact on the whole performance of the code.
       * 
       * DON'T touch *anything* in here unless you really know 
       * what you're doing !!
       *------------------------------------------------------------
       */

      Xalpha_k = Xalpha + k0 - freqIndex0;  /* first frequency-bin in sum */
      x0 = xhat_alpha - k0;	/* first xhat-value in the loop */
      /* count down 2*Dterms values */
      for ( k = 2 * params->Dterms; k != 0;  k -- )
	{
	  REAL4 realP, imagP;	/* real and imaginary parts of Dirichlet-kernel P_alpha_k */
	  REAL4 Xre, Xim;

	  REAL4 xinv = OOTWOPI_FLOAT / x0;

	  /* calculate P_alpha_k */
	  realP = sinx * xinv;
	  imagP = cosxm1 * xinv;

	  Xre = Xalpha_k->re;
	  Xim = Xalpha_k->im;

	  /* calculate P_alpha_k * X_alpha_k */
	  realXP += realP * Xre - imagP * Xim;
	  imagXP += imagP * Xre + realP * Xim;

	  Xalpha_k ++;	/* point to next frequency-bin */
	  x0 -- ;	/* x0-value for next iteration */

	} /* for k=kstar-Dterms to kstar+Dterms */
      
      realQXP = realQ * realXP - imagQ * imagXP;
      imagQXP = realQ * imagXP + imagQ * realXP;

      /* we're done: ==> combine these into Fa and Fb */
      Fa.re += a * realQXP;
      Fa.im += a * imagQXP;
      
      Fb.re += b * realQXP;
      Fb.im += b * imagQXP;

    } /* for alpha < numSFTs */
      
  /* return result */
  FaFb->Fa = Fa;
  FaFb->Fb = Fb;

  return 0;

} /* XLALNewLALDemod() */



/** Caculate F-statistic using XLALNewLALDemod() 
 */
int
XLALcomputeFStat (REAL8 *Fval,			/**< [out] the resulting F-statistic value */
		  const SFTVector *sfts, 	/**< input: SFT-vector */
		  const computeFStatPar *params /**< antenna-pattern coefficients a,b,..*/
		  )
{
  Fcomponents FaFb;
  REAL4 fact;
  REAL4 At, Bt, Ct;
  REAL8 FaRe, FaIm, FbRe, FbIm;

  if (!sfts || !params || !Fval) 
    {
      LALPrintError ("\nInput contains illegal NULL-pointer !\n\n");
      XLAL_ERROR ("XLALcomputeFStat", XLAL_EINVAL);
    }

  /* prepare quantities to calculate Fstat from Fa and Fb */
  fact = 4.0f / (1.0f * sfts->length * params->amcoe->D);
  At = params->amcoe->A;
  Bt = params->amcoe->B;
  Ct = params->amcoe->C;
  
  if ( XLALNewLALDemod (&FaFb, sfts, params) != 0) 
    {
      LALPrintError ("\nXALNewLALDemod() failed\n");
      XLAL_ERROR ("XLALcomputeFStat", XLAL_EFUNC);
    }


  FaRe = FaFb.Fa.re;
  FaIm = FaFb.Fa.im;
      
  FbRe = FaFb.Fb.re;
  FbIm = FaFb.Fb.im;

  /* calculate F-statistic from  Fa and Fb */
  (*Fval) = fact * (Bt * (FaRe*FaRe + FaIm*FaIm) 
		     + At * (FbRe*FbRe + FbIm*FbIm) 
		     - 2.0 * Ct *(FaRe*FbRe + FaIm*FbIm) );
  
  return 0;

} /* XLALcomputeFStat() */


/** write out the F-statistic over the searched frequency-band.
 */
void
writeFVect(LALStatus *status, const FStatisticVector *FVect, const CHAR *fname)
{
  FILE *fp;
  UINT4 i;
  REAL8 fi;

  INITSTATUS( status, "LALDemod", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (FVect, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (fname, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  if ( (fp = fopen(fname, "wb")) == NULL) 
    {
      LALPrintError ("Failed to open file '%f' for writing.\n", fname);
      ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
    }

  for (i=0; i < FVect->length; i++) {
    fi = FVect->f0 + 1.0*i*(FVect->df);

    fprintf (fp, "%20.17f %20.17f\n", fi, FVect->F->data[i]);
  } /* for i < FVect->length */

  fclose (fp);

  DETATCHSTATUSPTR (status);
  RETURN( status );

} /* writeCOMPLEX16Vector() */


/** Convert a whole vector of detector-times 'GPStimes' into SSB-frame 'SSBtimes'.
 * Using ephemeris-timing if ephem!=NULL, or using the Ptolemaic approximation otherwise.
 */
void
LALGetSSBtimes (LALStatus *status, 
		REAL8Vector *DeltaT,		/**< [out] DeltaT_alpha = T(t_alpha) - T(t0)*/
		REAL8Vector *Tdot,		/**< [out] Tdot(t_alpha) */
		LIGOTimeGPS refTime,		/**< reference-time t0 */
		const LIGOTimeGPSVector *GPStimes,/**< input detector times t_i */
		const SkyPosition *pos,		/**< source sky-location */
		const LALDetector *site,	/**< detector location and orientation */  
		const EphemerisData *ephem	/**< ephemeris-data, NULL for Ptolemaic timing*/
		)
{
  PulsarTimesParamStruc baryParams = empty_PulsarTimesParamStruc;
  REAL8Vector *var = NULL;	/* input-params: (t, alpha, delta) */
  REAL8Vector *dt = NULL;	/* hold answer (tau, dtau/dt, dtau/dalpha, dtau/ddelta ) */
  UINT4 N, i;
  REAL8 t0;

  INITSTATUS( status, "LALConvertGPS2SSBVector", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (DeltaT, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (Tdot, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (GPStimes, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (DeltaT->length == GPStimes->length, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT (Tdot->length == GPStimes->length, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT (pos, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (site, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  /*NOTE: ephem is allowed to be NULL ==> use Ptolemaic approximation */
  
  ASSERT ( pos->system == COORDINATESYSTEM_EQUATORIAL, status,
	   COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);

  N = GPStimes->length;	/* number of timestamps */

  /* Set up constant parameters for barycentre transformation. */
  baryParams.epoch = refTime;
  baryParams.t0 = 0;

  baryParams.latitude = site->frDetector.vertexLatitudeRadians;	
  baryParams.longitude = site->frDetector.vertexLongitudeRadians;

  baryParams.site = site;
  baryParams.ephemeris = ephem;

  /* set some time-constants depending on epoch */
  TRY ( LALGetEarthTimes( status->statusPtr, &baryParams ), status);

  /* create the 'variables'-vector : (t, alpha, delta) */
  TRY ( LALDCreateVector ( status->statusPtr, &var, 3 ), status);
  TRY ( LALDCreateVector ( status->statusPtr, &dt, 4 ), status);

  var->data[1] = pos->longitude;
  var->data[2] = pos->latitude;

  /*---------- first get reference-time t0 */
  var->data[0] = 0.0; 
  if ( ephem ) {
    TRY ( LALDTEphemeris( status->statusPtr, dt, var, &baryParams ), status );
  } else {
    TRY ( LALDTBaryPtolemaic( status->statusPtr, dt, var, &baryParams ), status );
  }
  t0 = dt->data[0];
  
  /*---------- Now loop over timestamps */
  for (i=0; i < N; i++ )
    {
      var->data[0] = XLALDeltaFloatGPS( &(GPStimes->data[i]), &refTime );

      /* use timing-function for earth-motion: either ptolemaic or ephemeris */
      if ( ephem ) {
	TRY ( LALDTEphemeris( status->statusPtr, dt, var, &baryParams ), status );
      } else {
	TRY ( LALDTBaryPtolemaic( status->statusPtr, dt, var, &baryParams ), status );
      }

      DeltaT->data[i] = dt->data[0] - t0;	/* DeltaT_alpha */
      Tdot->data[i] = dt->data[1];		/* Tdot_alpha */

    } /* for i <= N */

  /* free memory */
  TRY ( LALDDestroyVector ( status->statusPtr, &var ), status);
  TRY ( LALDDestroyVector ( status->statusPtr, &dt ), status);


  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetSSBtimes() */

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


