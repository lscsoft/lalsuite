/*
 * 
 * Copyright (C) 2004, 2005 Reinhard Prix
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
#include <lal/SFTfileIO.h>
#include <lal/ExtrapolatePulsarSpins.h>

#include <lalapps.h>

/* local includes */
#include "clusters.h"
#include "DopplerScan.h"

RCSID( "$Id$");

/*---------- DEFINES ----------*/

#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

#define EPHEM_YEARS  "00-04"	/**< default range: override with --ephemYear */

#define TRUE (1==1)
#define FALSE (1==0)

/*----- SWITCHES -----*/
#define FILE_FSTATS 1		/**< write out an 'Fstats' file containing cluster-output */


/*----- Error-codes -----*/
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

/*----- Macros -----*/
/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

/*---------- internal types ----------*/

/* ----- Output types for LALGetDetectorStates() */
/** State-info about position, velocity and LMST of a detector together 
 * with corresponding EarthState.
 */
typedef struct
{
  LIGOTimeGPS tGPS;		/**< GPS timestamps corresponding to this entry */
  REAL8 rDetector[3];		/**< Cartesian coords of detector position in ICRS J2000. Units=sec */
  REAL8 vDetector[3];		/**< Cart. coords. of detector velocity, in dimensionless units (v/c)*/
  REAL8 LMST;			/**< local mean sidereal time at the detector-location in radians */
  EarthState earthState;	/**< pointer to EarthState information */
} DetectorState;

/** Timeseries of DetectorState's, representing the detector-info at different timestamps.
 * In addition to the standard 'vector'-fields we also store the detector-info in here.
 */
typedef struct
{
  UINT4 length;			/**< total number of entries */
  DetectorState *data;		/**< array of DetectorState entries */
  LALDetector detector;		/**< detector-info corresponding to this timeseries */
} DetectorStateSeries;



typedef struct {
  REAL8Vector	*fkdot;		/**< vector of frequency + derivatives (spindowns) */
  REAL8Vector	*DeltaT;	/**< vector of DeltaT_alpha's (depend on skyposition)*/
  REAL8Vector	*Tdot;		/**< vector of Tdot_alpha's (depend on skyposition)*/ 
  AMCoeffs      *amcoe;         /**< Amplitude Modulation coefficients */
  INT4          Dterms;         /**< Terms used in the computation of the dirichlet kernel*/
} computeFStatPar;


/** Detectors Vector; specify's the number of detectors and SFTs.*/
typedef struct {
  UINT4 length;
  LALDetector *Detectors;         
  SFTVector **sftVects;
  DetectorStateSeries **detStates;
  LIGOTimeGPSVector **timestamps;	/**< SFT timestamps */
  LIGOTimeGPSVector **midTS;		/**< GPS midpoints of SFT's */
  DetectorStateSeries **DetectorStates;	/**< pos, vel and LMSTs for detector at times t_i */
} IFOspecifics;

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  CHAR EphemEarth[MAXFILENAMELENGTH];	/**< filename of earth-ephemeris data */
  CHAR EphemSun[MAXFILENAMELENGTH];	/**< filename of sun-ephemeris data */
  REAL8 dFreq;			/**< frequency resolution */
  REAL8 df1dot;			/**< spindown resolution (f1 = df/dt!!) */
  LIGOTimeGPS startTime;	/**< start time of observation */
  LIGOTimeGPS refTime;		/**< reference-time for pulsar-parameters in SBB frame */
  REAL8 refTime0;		/**< *internal* SSB reference time: e.g. start of observation */
  REAL8Vector *fkdot0;		/**< start frequency- and spindowns- at internal reference-time */
  EphemerisData *edat;		/**< ephemeris data (from LALInitBarycenter()) */
  CHAR *skyRegionString;	/**< sky-region to search (polygon defined by list of points) */
  computeFStatPar CFSparams;  	/**< Demodulation parameters for computeFStat() */
  IFOspecifics ifos;		/**< IFO-specific configuration parameters */
} ConfigVariables;

/** FIXME: OBSOLETE: used to hold result from LALDemod(). 
    Only kept for the moment to make things work (FIXME)*/
typedef struct {
  const REAL8         *F;       /* Array of value of the F statistic */
  const COMPLEX16     *Fa;      /* Results of match filter with a(t) */
  const COMPLEX16     *Fb;      /* Results of match filter with b(t) */
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

/** The precision in calculating the barycentric transformation */
typedef enum {
  SSBPREC_NEWTONIAN,		/**< simple Newtonian: \f$\tau = t + \vec{r}\cdot\vec{n}/c\f$ */
  SSBPREC_RELATIVISTIC,		/**< detailed relativistic: \f$\tau=\tau(t; \vec{n}, \vec{r})\f$ */
  SSBPREC_LAST			/**< end marker */
} SSBprecision;




/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/* *ifos.sftVects[0] = NULL;*/	/**< holds the SFT-data to analyze */
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
REAL8 uvar_FstatMin;
CHAR *uvar_skyGridFile;
CHAR *uvar_outputSkyGrid;
CHAR *uvar_workingDir;
REAL8 uvar_dopplermax;
INT4 uvar_windowsize;
REAL8 uvar_refTime;
INT4 uvar_SSBprecision;

/* ---------- local prototypes ---------- */

int main(int argc,char *argv[]);
void initUserVars (LALStatus *);

void InitFStatCommon (LALStatus *, ConfigVariables *cfg);
void InitFStat (LALStatus *, ConfigVariables *cfg);

void CreateNautilusDetector (LALStatus *, LALDetector *Detector);
void Freemem(LALStatus *,  ConfigVariables *cfg);


void NormaliseSFTDataRngMdn (LALStatus *);
void WriteFStatLog (LALStatus *, CHAR *argv[]);
void writeFVect(LALStatus *, const FStatisticVector *FVect, const CHAR *fname);
void checkUserInputConsistency (LALStatus *lstat);

int
XLALNewLALDemod(Fcomponents *FaFb,
		const SFTVector *sfts, 
		const computeFStatPar *params);

int 
XLALcomputeFStat (REAL8 *Fval, const SFTVector *sfts, const computeFStatPar *params);

void
LALGetSSBtimes (LALStatus *, 
		REAL8Vector *DeltaT, REAL8Vector *Tdot, 
		const DetectorStateSeries *DetStates, 
		SkyPosition pos,
		REAL8 refTime,
		SSBprecision precision);

void
LALGetAMCoeffs(LALStatus *status,
	       AMCoeffs *coeffs, 
	       const DetectorStateSeries *DetStates,
	       SkyPosition skypos);


const char *va(const char *format, ...);	/* little var-arg string helper function */
int sin_cos_LUT (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x); /* LUT-calculation of sin/cos */


void LALCreateDetectorStateSeries (LALStatus *, DetectorStateSeries **vect, UINT4 length );
void LALDestroyDetectorStateSeries(LALStatus *, DetectorStateSeries **vect );

void
LALGetDetectorStates (LALStatus *, DetectorStateSeries **DetStates,
		      const LIGOTimeGPSVector *timestamps,
		      const LALDetector *detector,
		      const EphemerisData *edat);

/****** black list ******/
void EstimateFLines(LALStatus *, const FStatisticVector *FVect);
INT4 writeFLines(INT4 *maxIndex, REAL8 f0, REAL8 df);
int compare(const void *ip, const void *jp);
INT4 writeFaFb(INT4 *maxIndex);

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

  INT4 *maxIndex=NULL; 			/*  array that contains indexes of maximum of each cluster */
  CHAR Fstatsfilename[256]; 		/* Fstats file name*/
  DopplerScanInit scanInit;		/* init-structure for DopperScanner */
  DopplerScanState thisScan = empty_DopplerScanState; /* current state of the Doppler-scan */
  DopplerPosition dopplerpos;		/* current search-parameters */
  SkyPosition thisPoint;
  FILE *fpOut=NULL;
  UINT4 loopcounter;
  FStatisticVector *FVect = NULL;   /* new type to store F-statistic results in a frequency-band */

  UINT4 nBins; 			/* number of frequency-bins */
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
  LAL_CALL ( InitFStatCommon(&status, &GV), &status);

  /**---------------------------------------------------------**/
  /** Starting the Loop for different Detectors **/
  /** At this moment we are trying to match it for singel detector **/


  for(nD=0; nD < GV.ifos.length ; nD++)
    {
      /* main initialization of the code: */
      LAL_CALL ( InitFStat(&status, &GV), &status);
      
      /* normalize SFTs by running median */
      LAL_CALL (NormaliseSFTDataRngMdn(&status), &status);
      
    }
  


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

  /*----- figure out total observation time */
  {
    LIGOTimeGPS t0, t1;
    UINT4 numSFTs = GV.ifos.sftVects[0]->length;
    REAL8 tObs;

    t0 = GV.ifos.sftVects[0]->data[0].epoch;
    t1 = GV.ifos.sftVects[0]->data[numSFTs-1].epoch;
    LAL_CALL (LALDeltaFloatGPS (&status, &tObs, &t1, &t0), &status);	/* t1 - t0 */
    tObs += 1.0 / (GV.ifos.sftVects[0]->data[0].deltaF );			/* +tSFT */
    GV.startTime = t0;

    scanInit.obsDuration = tObs;
  }

  scanInit.obsBegin = GV.startTime;
  /* scanInit.fmax  = uvar_Freq + uvar_FreqBand;*/
  scanInit.Detector = &(GV.ifos.Detectors[0]);
  scanInit.ephemeris = GV.edat;		/* used by Ephemeris-based metric */
  scanInit.searchRegion.skyRegionString = GV.skyRegionString;
  scanInit.skyGridFile = uvar_skyGridFile;
  
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
  
      /* prepare memory to hold F-statistic array over frequency (for cluster-stuff) [FIXME]*/
  nBins =  (UINT4)(uvar_FreqBand/GV.dFreq + 0.5) + 1;
  if ( (FVect = (FStatisticVector*)LALCalloc(1, sizeof(FStatisticVector))) == NULL) {
    LALPrintError ("\nOut of memory..!\n\n");
    return (COMPUTEFSTATC_EMEM);
  }
  LAL_CALL ( LALDCreateVector (&status, &(FVect->F), nBins), &status);
  
  FVect->f0 = GV.fkdot0->data[0];
  FVect->df = GV.dFreq;
  FVect->fBand = (nBins - 1) * FVect->df;
  FVect->length = nBins;
  
  /* obtain detector positions and velocities, together with LMSTs for the 
   * SFT midpoints 
   */
  LAL_CALL(LALGetDetectorStates(&status, &(GV.ifos.DetectorStates[0]), GV.ifos.midTS[0], &(GV.ifos.Detectors[0]), GV.edat), &status);
  
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
      
      nFreq =  (UINT4)(uvar_FreqBand  / GV.dFreq + 0.5) + 1;  
      nf1dot = (UINT4)(uvar_f1dotBand / GV.df1dot+ 0.5) + 1; 
      
      LAL_CALL (NextDopplerPos( &status, &dopplerpos, &thisScan ), &status);
      if (thisScan.state == STATE_FINISHED) /* scanned all DopplerPositions yet? */
	break;
      
      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
      thisPoint.longitude = dopplerpos.Alpha;
      thisPoint.latitude = dopplerpos.Delta;
      thisPoint.system = COORDINATESYSTEM_EQUATORIAL;
      LAL_CALL (LALNormalizeSkyPosition(&status, &thisPoint, &thisPoint), &status);
      
      
      /*----- loop over first-order spindown values */
      for (if1dot = 0; if1dot < nf1dot; if1dot ++)
	{
	  GV.CFSparams.fkdot->data[1] = GV.fkdot0->data[1] + if1dot * GV.df1dot;
	  
	  /* Loop over frequencies to be demodulated */
	  for ( iFreq = 0 ; iFreq < nFreq ; iFreq ++ )
	    {
	      GV.CFSparams.fkdot->data[0] = GV.fkdot0->data[0] + iFreq * GV.dFreq;	
	      
	      for(nD=0; nD < GV.ifos.length; nD++)
		{
		  
		  /*----- calculate SSB-times DeltaT_alpha and Tdot_alpha for this skyposition */
		  LAL_CALL ( LALGetSSBtimes (&status, GV.CFSparams.DeltaT, GV.CFSparams.Tdot, 
					     GV.ifos.DetectorStates[nD], thisPoint, GV.refTime0,
					     uvar_SSBprecision), &status);
		  
		  /*----- calculate skypos-specific coefficients a_i, b_i, A, B, C, D */
		  LAL_CALL ( LALGetAMCoeffs (&status, GV.CFSparams.amcoe, GV.ifos.DetectorStates[nD], thisPoint), &status);
		  
		  /** Caculate F-statistic using XLALNewLALDemod() */
		  
		  {
		    Fcomponents FaFb;
		    REAL4 fact;
		    REAL4 At, Bt, Ct;
		    REAL4 FaRe, FaIm, FbRe, FbIm;
		    
		    /* prepare quantities to calculate Fstat from Fa and Fb */
		    fact = 4.0f / (1.0f * GV.ifos.sftVects[0]->length * GV.CFSparams.amcoe->D);
		    At = GV.CFSparams.amcoe->A;
		    Bt = GV.CFSparams.amcoe->B;
		    Ct = GV.CFSparams.amcoe->C;
		    
		    if ( XLALNewLALDemod (&FaFb, GV.ifos.sftVects[nD], &(GV.CFSparams)) != 0) 
		      {
			LALPrintError ("\nXALNewLALDemod() failed\n");
			XLAL_ERROR ("XLALcomputeFStat", XLAL_EFUNC);
		      }
		    
		    
		    FaRe = FaFb.Fa.re;
		    FaIm = FaFb.Fa.im;
		    
		    FbRe = FaFb.Fb.re;
		    FbIm = FaFb.Fb.im;
		    
		    /* calculate F-statistic from  Fa and Fb */
		    FVect->F->data[iFreq] = fact * (Bt * (FaRe*FaRe + FaIm*FaIm) 
						    + At * (FbRe*FbRe + FbIm*FbIm) 
						    - 2.0f * Ct *(FaRe*FbRe + FaIm*FbIm) );
		    
		  } /* XLALcomputeFStat() */
		  
		} /* End of loop over detectors */
	      
	    } /* for i < nBins: loop over frequency-bins */
	  
	  
	      /* now, if user requested it, we output ALL F-statistic results above threshold */
	  if ( fpOut )
	    {
	      UINT4 j;
	      for(j=0; j < FVect->length; j++)
		{
		  REAL8 FStat = 2.0 * medianbias * FVect->F->data[j];
		  REAL8 freq = uvar_Freq + j * GV.dFreq;
		  
		  if ( FStat > uvar_Fthreshold )
		    fprintf (fpOut, "%16.12f %8.7f %8.7f %.17g %10.6g\n", 
			     freq, dopplerpos.Alpha, dopplerpos.Delta, 
			     GV.CFSparams.fkdot->data[1], FStat);
		}
	      
	    } /* if outputFstat */
	  
	  
	      /* FIXME: to keep cluster-stuff working, we provide the "translation" 
	       * from FVect back into old Fstats-struct */
	  Fstat.F  = FVect->F->data;
	  
	  /*  This fills-in highFLines */
	  if (nFreq > 5) {
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
	printf ("\
Search progress: %5.1f%%", (100.0* loopcounter / thisScan.numGridPoints));
      
    } /*  while SkyPos */
  
  if (uvar_outputFstat && fpOut)
    fclose (fpOut);
  
  if (lalDebugLevel) printf ("\nSearch finished.\n");
  
  /* properly terminate Fstats-file by 'DONE' marker: */ 
  if (fpstat) fprintf(fpstat, "%%DONE\n");
  if (fpstat) fclose(fpstat);
  
  /* Free memory */
  LAL_CALL ( FreeDopplerScan(&status, &thisScan), &status);
  
  LAL_CALL ( Freemem(&status, &GV), &status);
  
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

  uvar_SSBprecision = SSBPREC_RELATIVISTIC;

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
  LALregSTRINGUserVar(status,	DataFiles, 	'D', UVAR_REQUIRED, "Directory where SFT's are located"); 
  LALregSTRINGUserVar(status,	DataFiles2, 	 0,  UVAR_OPTIONAL, "Directory where SFT's are located");
  LALregSTRINGUserVar(status,	ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(status,	ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregSTRINGUserVar(status, 	IFO, 		'I', UVAR_REQUIRED, "Detector: GEO(0), LLO(1), LHO(2), NAUTILUS(3), VIRGO(4), TAMA(5), CIT(6)");
  LALregSTRINGUserVar(status, 	IFO2, 		 0,  UVAR_OPTIONAL, "Detector: GEO(0), LLO(1), LHO(2), NAUTILUS(3), VIRGO(4), TAMA(5), CIT(6)"); 
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
  LALregSTRINGUserVar(status,   workingDir,     'w', UVAR_OPTIONAL, "Directory to be made the working directory, . is default");
  LALregREALUserVar(status,	refTime,	 0,  UVAR_OPTIONAL, "SSB reference time for pulsar-paramters");

  /* more experimental and unofficial stuff follows here */
  LALregSTRINGUserVar(status,	outputFstat,	 0,  UVAR_OPTIONAL, "Output the F-statistic field over the parameter-space");
  LALregREALUserVar(status, 	FstatMin,	 0,  UVAR_OPTIONAL, "Minimum F-Stat value to written into outputFstat-file");
  LALregINTUserVar (status, 	SSBprecision,	 0,  UVAR_DEVELOPER, "Precision to use for time-transformation to SSB: 0=Newtonian 1=relativistic");

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* initUserVars() */

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
InitFStatCommon (LALStatus *status, ConfigVariables *cfg)
{
  UINT4 nDet;

  INITSTATUS (status, "InitFStat", rcsid);
  ATTATCHSTATUSPTR (status);

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
   * initialize+check  template-grid related parameters 
   */
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

    TRY (LALLeapSecs (status->statusPtr, &leap, &(cfg->startTime), &formatAndAcc), status);
    cfg->edat->leap = leap;

    TRY (LALInitBarycenter(status->statusPtr, cfg->edat), status);               

  } /* end: init ephemeris data */


  /* ----- count number of detectors */
  if ( LALUserVarWasSet(&uvar_DataFiles2) )
    nDet = 2;
  else
    nDet = 1;

  cfg->ifos.length = nDet;

  /* ----- set up detector-specific data */
  cfg->ifos.Detectors = LALCalloc ( nDet,  sizeof( *(cfg->ifos.Detectors) ) );
  cfg->ifos.sftVects =  LALCalloc ( nDet,  sizeof( *(cfg->ifos.sftVects) ) );
  cfg->ifos.detStates =  LALCalloc ( nDet,  sizeof( *(cfg->ifos.detStates) ) );
  cfg->ifos.timestamps =  LALCalloc ( nDet,  sizeof( *(cfg->ifos.timestamps) ) );
  cfg->ifos.midTS =  LALCalloc ( nDet,  sizeof( *(cfg->ifos.midTS) ) );
  cfg->ifos.DetectorStates =  LALCalloc ( nDet,  sizeof( *(cfg->ifos.DetectorStates) ) );


  DETATCHSTATUSPTR (status);
  RETURN (status);


} /* InitFStatCommon() */


void
InitFStat (LALStatus *status, ConfigVariables *cfg)
{
  UINT4 i;

  INITSTATUS (status, "InitFStat", rcsid);
  ATTATCHSTATUSPTR (status);

  /*----------------------------------------------------------------------
   * load SFT data-files
   */
  /* min and max physical frequency to read */
  {
    REAL8 f_min, f_max; 
    f_min = uvar_Freq;			/* lower physical frequency requested by user */
    f_max = f_min + uvar_FreqBand; 	/* upper physical frequency requested by user */

    /* NOTE: the following correction for the frequency-drift due to spindown
     * requires to know the total observation time, which right now we don't
     * have access to (haven't read the SFTs yet..).
     * ==> we either need a new function in SFTIO-lib, or we load the 
     * SFTs twice if the initial guess turns out to be to small .. (->easy)
     * but it's not a high priority right now, as spindowns are not yet the 
     * focus of v2-development [FIXME]
     */
    /* correct for spindown-shift of frequency: extend the frequency-band */
#if 0
      {
	REAL8 f1dot_1 = uvar_f1dot;
	REAL8 f1dot_2 = f1dot_1 + uvar_f1dotBand;
	REAL8 f1dot_max = fabs(f1dot_1) > fabs(f1dot_2) ? f1dot_1 : f1dot_2;
	REAL8 df = f1dot_max * (GV.tObs);	/* don't use!! undefined tObs */
	if ( df < 0)
	  f_min += df;
	else
	  f_max += df;
#endif
    
    /* ----- correct for maximal dopper-shift due to earth's motion */
    f_min *= (1.0 - uvar_dopplermax);
    f_max *= (1.0 + uvar_dopplermax);
    
    /* ----- contruct file-patterns and load the SFTs */
    
    if (!uvar_DataFiles)
      strcpy (uvar_DataFiles, ".");
    
    TRY ( LALReadSFTfiles(status->statusPtr, &(cfg->ifos.sftVects[0]), f_min, f_max, uvar_Dterms, uvar_DataFiles), status);

  } /* SFT-loading */

  /*----------  prepare vectors of timestamps ---------- */
  {
    TRY ( LALCreateTimestampVector (status->statusPtr, &(cfg->ifos.timestamps[0]), cfg->ifos.length ), status);
    TRY ( LALCreateTimestampVector (status->statusPtr, &(cfg->ifos.midTS[0]), cfg->ifos.length ), status);

    for (i=0; i < cfg->ifos.length; i++)
      {
	cfg->ifos.timestamps[0]->data[i] = cfg->ifos.sftVects[0]->data[i].epoch;	/* SFT start-timestamps */
	/* SFT midpoints */
	TRY (LALAddFloatToGPS(status->statusPtr, 
			      &(cfg->ifos.midTS[0]->data[i]), &(cfg->ifos.timestamps[0]->data[i]),0.5*1.0 / (GV.ifos.sftVects[0]->data[0].deltaF )),status);
      }/* for i < ifos.length */

  } 

  /*---------- Standardise reference-time: ----------*/
  /* translate spindown-paramters {f, fdot, fdotdot..} from the user-specified 
   * reference-time uvar_refTime to the internal reference-time, which 
   * we chose as the start-time of the first SFT (*verbatim*, i.e. not translated to SSB! )
   */
  {
    UINT4 spdnOrder = 1;	/* hard-coded default FIXME. DON'T change without fixing main() */

    REAL8Vector *fkdotRef = NULL;
    LIGOTimeGPS refTime0;	/* internal reference-time */

    if ( LALUserVarWasSet(&uvar_refTime)) {
      TRY ( LALFloatToGPS (status->statusPtr, &(cfg->refTime), &uvar_refTime), status);
    } else
      cfg->refTime = cfg->startTime;

    TRY ( LALDCreateVector (status->statusPtr, &fkdotRef, 1 + spdnOrder), status);
    TRY ( LALDCreateVector (status->statusPtr, &(cfg->fkdot0), 1 + spdnOrder), status);
    fkdotRef->data[0] = uvar_Freq;
    if ( spdnOrder > 0 )
      fkdotRef->data[1] = uvar_f1dot;	    /* currently not more spindowns implemented... */

    /* currently we use the observation GPS start-time as internal SSB reference-time: */
    refTime0 = cfg->startTime;
    cfg->refTime0 = GPS2REAL8 (refTime0);

    /*----- now translate spin-params to internal reference-time */
    if ( XLALExtrapolatePulsarSpins ( cfg->fkdot0, refTime0, fkdotRef, cfg->refTime) ) 
      {
	int code = xlalErrno;
	XLALClearErrno(); 
	LALPrintError ("\nERROR: XLALExtrapolatePulsarSpins() failed (xlalErrno = %d)!\n\n", code);
	ABORT (status,  COMPUTEFSTATC_EXLAL,  COMPUTEFSTATC_MSGEXLAL);
      }

    TRY ( LALDDestroyVector (status->statusPtr, &fkdotRef), status);
  }



  /*----------------------------------------------------------------------
   * initialize detector 
   */
  if ( !strcmp (uvar_IFO, "GEO") || !strcmp (uvar_IFO, "0") ) 
    cfg->ifos.Detectors[0] = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  else if ( !strcmp (uvar_IFO, "LLO") || ! strcmp (uvar_IFO, "1") ) 
    cfg->ifos.Detectors[0] = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if ( !strcmp (uvar_IFO, "LHO") || !strcmp (uvar_IFO, "2") )
    cfg->ifos.Detectors[0] = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else if ( !strcmp (uvar_IFO, "NAUTILUS") || !strcmp (uvar_IFO, "3"))
    {
      TRY (CreateNautilusDetector (status->statusPtr, &(cfg->ifos.Detectors[0])), status);
    }
  else if ( !strcmp (uvar_IFO, "VIRGO") || !strcmp (uvar_IFO, "4") )
    cfg->ifos.Detectors[0] = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  else if ( !strcmp (uvar_IFO, "TAMA") || !strcmp (uvar_IFO, "5") )
    cfg->ifos.Detectors[0] = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  else if ( !strcmp (uvar_IFO, "CIT") || !strcmp (uvar_IFO, "6") )
    cfg->ifos.Detectors[0] = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
  else
    {
      LALPrintError ("\nUnknown detector. Currently allowed are \
'GEO', 'LLO', 'LHO', 'NAUTILUS', 'VIRGO', 'TAMA', 'CIT' or '0'-'6'\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }


  /* ----------------------------------------------------------------------
   * initialize + allocate space for AM-coefficients and Demod-params
   */
  {
    AMCoeffs *amc;
    REAL8Vector *DeltaT = NULL;
    REAL8Vector *Tdot = NULL;

    /* Allocate space for AMCoeffs */
    if ( (amc = LALCalloc(1, sizeof(AMCoeffs))) == NULL) {
      ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }
    TRY (LALSCreateVector(status->statusPtr, &(amc->a), (UINT4) cfg->ifos.length), status);
    TRY (LALSCreateVector(status->statusPtr, &(amc->b), (UINT4) cfg->ifos.length), status);

    cfg->CFSparams.amcoe = amc;
    
    /* allocate memory of the SSB-times: DeltaT_alpha and Tdot_alpha */
    TRY ( LALDCreateVector(status->statusPtr, &DeltaT, cfg->ifos.length), status );
    TRY ( LALDCreateVector(status->statusPtr, &Tdot,   cfg->ifos.length), status );

    cfg->CFSparams.DeltaT = DeltaT;
    cfg->CFSparams.Tdot = Tdot;

    cfg->CFSparams.Dterms = uvar_Dterms;

    /* prepare memory for fkdot - vector : (f, f1dot, f2dot, ..) */
    TRY(LALDCreateVector(status->statusPtr, &(cfg->CFSparams.fkdot), cfg->fkdot0->length), status);
    
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
  TRY (LALDestroySFTVector (status->statusPtr, &(cfg->ifos.sftVects[0]) ), status);	 /* the new way*/

  /* Free timestamps */
  TRY ( LALDestroyTimestampVector (status->statusPtr, &(cfg->ifos.timestamps[0]) ), status );
  TRY ( LALDestroyTimestampVector (status->statusPtr, &(cfg->ifos.midTS[0]) ), status );

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

  TRY (LALDDestroyVector (status->statusPtr, &(cfg->fkdot0)), status);

  TRY (LALDDestroyVector (status->statusPtr, &(cfg->CFSparams.fkdot)), status);
  TRY (LALSDestroyVector(status->statusPtr, &(cfg->CFSparams.amcoe->a)), status);
  TRY (LALSDestroyVector(status->statusPtr, &(cfg->CFSparams.amcoe->b)), status);
  LALFree ( cfg->CFSparams.amcoe);
  TRY (LALDDestroyVector(status->statusPtr, &(cfg->CFSparams.DeltaT)), status);
  TRY (LALDDestroyVector(status->statusPtr, &(cfg->CFSparams.Tdot)), status);


  /* destroy DetectorStateSeries */
  TRY ( LALDestroyDetectorStateSeries (status->statusPtr, &(GV.ifos.DetectorStates[0]) ), status);



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

  TRY ( LALDCreateVector(status->statusPtr, &Sp, GV.ifos.sftVects[0]->data[0].data->length), status);
  TRY ( LALDCreateVector(status->statusPtr, &RngMdnSp, GV.ifos.sftVects[0]->data[0].data->length), status);

  if( (N = (REAL8 *) LALCalloc(GV.ifos.sftVects[0]->data[0].data->length,sizeof(REAL8))) == NULL) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  if( (Sp1 = (REAL8 *) LALCalloc(GV.ifos.sftVects[0]->data[0].data->length,sizeof(REAL8))) == NULL) { 
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  /* loop over each SFTs */
  for (i=0;i<GV.ifos.length;i++)         
    {
      /* Set to zero the values */
      for (j=0;j<GV.ifos.sftVects[0]->data[0].data->length;j++){
	RngMdnSp->data[j] = 0.0;
	Sp->data[j]       = 0.0;
      }
      
      /* loop over SFT data to estimate noise */
      for (j=0;j<GV.ifos.sftVects[0]->data[0].data->length;j++){
	xre=GV.ifos.sftVects[0]->data[i].data->data[j].re;
	xim=GV.ifos.sftVects[0]->data[i].data->data[j].im;
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
      for (lpc=0;lpc<GV.ifos.sftVects[0]->data[0].data->length;lpc++){
	N[lpc]=1.0/sqrt(2.0*RngMdnSp->data[lpc]);
      }
      
      if(uvar_SignalOnly == 1){
	B=(1.0*GV.ifos.sftVects[0]->data[0].data->length)/(1.0 / (GV.ifos.sftVects[0]->data[0].deltaF ));
	deltaT=1.0/(2.0*B);
	norm=deltaT/sqrt(1.0 / (GV.ifos.sftVects[0]->data[0].deltaF ));
	for (lpc=0;lpc<GV.ifos.sftVects[0]->data[0].data->length;lpc++){
	  N[lpc]=norm;
	}
      }
      
      /*  loop over SFT data to normalise it (with N) */
      /*  also compute Sp1, average normalized PSD */
      /*  and the sum of the PSD in the band, SpSum */
      for (j=0;j<GV.ifos.sftVects[0]->data[0].data->length;j++){
	xre=GV.ifos.sftVects[0]->data[i].data->data[j].re;
	xim=GV.ifos.sftVects[0]->data[i].data->data[j].im;
	xreNorm=N[j]*xre; 
	ximNorm=N[j]*xim; 
	GV.ifos.sftVects[0]->data[i].data->data[j].re = xreNorm;    
	GV.ifos.sftVects[0]->data[i].data->data[j].im = ximNorm;
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


#define LD_SMALL        (1.0e-9 / LAL_TWOPI)	/**< "small" number */
#define OOTWOPI         (1.0 / LAL_TWOPI)	/**< 1/2pi */

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */ 

/** v2-specific version of LALDemod() (based on TestLALDemod() in CFS)
 */
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
  REAL8 f;		/* !! MUST be REAL8, or precision breaks down !! */
  REAL8 Tsft; 			/* length of SFTs in seconds */

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

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;

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

      REAL8 xhat_alpha, y_alpha;	/* xhat(alpha), y(alpha): need to be REAL8 !! */
      REAL4 x0;
      UINT4 k;			/* loop index over frequency-bins */
      UINT4 kstar;		/* central frequency-bin k* = round(xhat_alpha) */

      COMPLEX8 *Xalpha = sfts->data[alpha].data->data; /* pointer to current SFT-data */
      COMPLEX8 *Xalpha_k; 	/* pointer to frequency-bin k in current SFT */
      REAL4 sinx, cosxm1;	/* sin(x_alpha) and (cos(x_alpha)-1) */
      REAL4 realXP, imagXP;	/* the sum_k X_alpha_k P_alpha_k */
      REAL4 realQ, imagQ;	/* Re and Im of Q = e^{-i y} */
      REAL4 realQXP, imagQXP;	/* Re/Im of Q_alpha XP_alpha */
      UINT4 k0;
      /* ----- calculate x(alpha,0) and y(alpha) */
      {
	UINT4 s; 		/* loop-index over spindown-order */
	REAL8 Tas; 		/* temporary variable to calculate (DeltaT_alpha)^2 */
	UINT4 sfact = 1;	/* store for s! */
	REAL8 DeltaTalpha = params->DeltaT->data[alpha];
	Tas = 1.0; 	/* DeltaT_alpha = T^1 */

	/* Step 1: s = 0 */
	xhat_alpha = f * Tas;	/* f^{0) T^0 / 0! */
	Tas *= DeltaTalpha;
	y_alpha = f * Tas;	/* f^{0} T^1 / 1! */

	/* Step 2: sum s >= 1 */
	for (s=1; s <= spdnOrder; s++)
	  {
	    REAL8 fsdot = params->fkdot->data[s];
	    xhat_alpha += fsdot * Tas / sfact; 	/* Tas = T^s here, sfact=s! */
	    Tas *= DeltaTalpha; 		/* T^(s+1) */
	    sfact *= (s+1);			/* (s+1)! */	  
	    y_alpha += fsdot * Tas / sfact; 
	  } /* for s <= spdnOrder */

	/* Step 3: apply global factors and complete y_alpha */
	xhat_alpha *= Tsft * params->Tdot->data[alpha];	/* guaranteed > 0 ! */
	y_alpha -= 0.5 * xhat_alpha;
	
	/* real- and imaginary part of e^{-i 2 pi y } */
	if ( sin_cos_LUT ( &imagQ, &realQ, y_alpha ) ) {
	  XLAL_ERROR ( "XLALNewLALDemod", XLAL_EFUNC);
	}
	imagQ = -imagQ;
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
      sin_cos_LUT ( &sinx, &cosxm1, xhat_alpha );
      cosxm1 -= 1.0f; 
      /*-------------------- */

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
      x0 = (REAL4)(xhat_alpha - (REAL8)k0);	/* first xhat-value in the loop */
      /* count down 2*Dterms values */
      for ( k = 2 * params->Dterms; k != 0;  k -- )
	{
	  REAL4 realP, imagP;	/* real and imaginary parts of Dirichlet-kernel P_alpha_k */
	  COMPLEX8 Xa = *Xalpha_k;
	  REAL4 xinv = OOTWOPI_FLOAT / x0;

	  /* calculate P_alpha_k */
	  realP = sinx * xinv;
	  imagP = cosxm1 * xinv;

	  /* calculate P_alpha_k * X_alpha_k */
	  realXP += realP * Xa.re - imagP * Xa.im;
	  imagXP += imagP * Xa.re + realP * Xa.im;

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

  return XLAL_SUCCESS;

} /* XLALNewLALDemod() */


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
      LALPrintError ("\nFailed to open file '%f' for writing.\n\n", fname);
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
  /* check for negative stepsizes in Freq, Alpha, Delta */
  if ( LALUserVarWasSet(&uvar_dAlpha) && (uvar_dAlpha < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dAlpha not allowed!\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar_dDelta) && (uvar_dDelta < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dDelta not allowed!\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar_dFreq) && (uvar_dFreq < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dFreq not allowed!\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
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
        ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
      }

    if ( !useGridFile && !(haveSkyRegion || haveAlphaDelta) )
      {
        LALPrintError ("\nNeed sky-region: either use (Alpha,Delta) or skyRegion!\n\n");
        ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
      }
    if ( haveSkyRegion && haveAlphaDelta )
      {
        LALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta)"
		       " OR skyRegion!\n\n");
        ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
      }
    if ( useGridFile && !haveGridFile )
      {
        LALPrintError ("\nERROR: gridType=FILE, but no skyGridFile specified!\n\n");
        ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);  
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
        ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);      
      }

  } /* grid-related checks */

  RETURN (status);
} /* checkUserInputConsistency() */

/** Calculate sin(2 pi x) and cos(2 pi x) to roughly 1e-7 error using 
 * a lookup-table and Tayler-expansion.
 * This is meant to be fast, so we don't even check the input-pointers...
 *
 * However, for numerical sanity&safty, we *DO* check if the resulting
 * index is within bounds, which can fail in case the argument x is too large..
 *
 * return = 0: OK, nonzero=ERROR
 */
int
sin_cos_LUT (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x)
{
#define LUT_RES         64      /* resolution of lookup-table */
  UINT4 ind; 
  REAL8 rem;
  static BOOLEAN firstCall = TRUE;
  static REAL4 sinVal[LUT_RES+1], cosVal[LUT_RES+1];


  if ( firstCall )
    {
      UINT4 k;
      for (k=0; k <= LUT_RES; k++)
        {
          sinVal[k] = sin( (LAL_TWOPI*k)/LUT_RES );
          cosVal[k] = cos( (LAL_TWOPI*k)/LUT_RES );
        }
      firstCall = FALSE;
    }

  rem = x - (INT8)x;	/* rem in (-1, 1) */
  if ( rem < 0 )
    rem += 1.0;		/* rem in [0, 1) */

  /* security check if we didn't overstretch the numerics here (can happen for x too large) */
  if ( (rem < 0) || (rem > 1) )
    {
      LALPrintError ("\nLUT-index out of bounds. Input argument was probably too large!\n\n");
      XLAL_ERROR ( "sin_cos_LUT", XLAL_EDOM);
    }
  
			   
  ind = (UINT4)( rem * LUT_RES + 0.5 );   /* closest LUT-entry */
  {
    REAL8 d = LAL_TWOPI *(rem - (REAL8)ind/(REAL8)LUT_RES);
    REAL8 d2 = 0.5 * d * d;
    REAL8 ts = sinVal[ind];
    REAL8 tc = cosVal[ind];
                
    (*sin2pix) = ts + d * tc - d2 * ts;
    (*cos2pix) = tc - d * ts - d2 * tc;
  }
  
  return XLAL_SUCCESS;
} /* sin_cos_LUT() */


/** Compute the 'amplitude coefficients' \f$a(t), b(t)\f$ as defined in 
 * \ref JKS98 for a series of timestamps.
 * 
 * The input consists of the DetectorState-timeseries, which contains
 * the detector-info and the LMST's corresponding to the different times.
 * 
 * In order to allow re-using the output-structure AMCoeffs for subsequent
 * calls, we require the REAL4Vectors a and b to be allocated already and 
 * to have the same length as the DetectoStates-timeseries.
 *
 * \note This is an alternative implementation to LALComputeAM() with 
 * the aim to be both simpler and faster.
 * The difference being that we don't implicitly re-derive the final expression
 * here but simply try to implement the final expressions (12), (13) in \ref JKS98
 * in the most economical way possible.
 */ 
void
LALGetAMCoeffs(LALStatus *status,
	       AMCoeffs *coeffs,			/**< [out] amplitude-coeffs {a(t_i), b(t_i)} */
	       const DetectorStateSeries *DetStates,	/**< timeseries of detector states */
	       SkyPosition skypos			/**< {alpha,delta} of the source */
	       )
{
  REAL4 ah1, ah2, ah3, ah4, ah5;
  REAL4 a1, a2, a3, a4, a5;
  
  REAL4 bh1, bh2, bh3, bh4;
  REAL4 b1, b2, b3, b4;

  REAL4 delta, alpha;
  REAL4 sin1delta, cos1delta, sin2delta, cos2delta;

  REAL4 gamma, lambda;
  REAL4 norm;
  UINT4 i, numSteps;

  INITSTATUS (status, "LALGetAMCoeffs", rcsid);

  /*---------- check input ---------- */
  ASSERT ( DetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  numSteps = DetStates->length;

  /* require the coeffients-vectors to be allocated and consistent with timestamps */
  ASSERT ( coeffs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( coeffs->a && coeffs->b, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( (coeffs->a->length == numSteps) && (coeffs->b->length == numSteps), status,
	   COMPUTEFSTATC_EINPUT,  COMPUTEFSTATC_MSGEINPUT);

  /* require sky-pos to be in equatorial coordinates */
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, 
	   SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );

  /*---------- detector paramters: lambda, L, gamma */
  {
    /* FIXME: put into DetectorStateSeries */
    /* orientation of detector arms */
    REAL8 xAzi = DetStates->detector.frDetector.xArmAzimuthRadians;
    REAL8 yAzi = DetStates->detector.frDetector.yArmAzimuthRadians;

    /* get detector orientation gamma */
    gamma = LAL_PI_2 - 0.5 * (xAzi + yAzi);
    /* get detector position latitude (lambda) */
    lambda = DetStates->detector.frDetector.vertexLatitudeRadians;
  }

  /*---------- coefficient ahN, bhN dependent ONLY on detector-position  ---------- */
  /* FIXME: put these coefficients into DetectorStateSeries */
  {
    REAL4 sin2gamma = sinf ( 2.0f * gamma );
    REAL4 cos2gamma = cosf ( 2.0f * gamma );
    REAL4 sin1lambda = sinf ( lambda );
    REAL4 cos1lambda = cosf ( lambda );
    REAL4 sin2lambda = 2.0f * sin1lambda * cos1lambda;
    REAL4 cos2lambda = cos1lambda * cos1lambda - sin1lambda * sin1lambda;

    /* coefficients for a(t) */
    ah1 = 0.0625f * sin2gamma * (3.0f - cos2lambda);	/* 1/16 = 0.0625 */
    ah2 = - 0.25f * cos2gamma * sin1lambda;
    ah3 =   0.25f * sin2gamma * sin2lambda;
    ah4 =  -0.5f  * cos2gamma * cos1lambda;
    ah5 =  0.75f  * sin2gamma * cos1lambda * cos1lambda;

    /* coefficients for b(t) */
    bh1 =           cos2gamma * sin1lambda;
    bh2 =   0.25f * sin2gamma * (3.0f - cos2lambda);
    bh3 =           cos2gamma * cos1lambda;
    bh4 =   0.5f  * sin2gamma * sin2lambda;
  }

  /*---------- coefficients aN, bN dependent ONLY on {ahN, bhN} and source-latitude delta */
  alpha = skypos.longitude;
  delta = skypos.latitude;
  sin1delta = sinf (delta);
  cos1delta = cosf (delta);
  sin2delta = 2.0f * sin1delta * cos1delta;
  cos2delta = cos1delta * cos1delta - sin1delta * sin1delta;

  /* coefficients for a(t) */
  a1 = ah1 * ( 3.0f - cos2delta );
  a2 = ah2 * ( 3.0f - cos2delta );
  a3 = ah3 * sin2delta;
  a4 = ah4 * sin2delta;
  a5 = ah5 * cos1delta * cos1delta;

  /* coefficients for b(t) */
  b1 = bh1 * sin1delta;
  b2 = bh2 * sin1delta;
  b3 = bh3 * cos1delta;
  b4 = bh4 * cos1delta;


  /*---------- Compute the a(t_i) and b(t_i) ---------- */
  coeffs->A = 0;
  coeffs->B = 0;
  coeffs->C = 0;
  coeffs->D = 0;
  for ( i=0; i < numSteps; i++ )
    {
      REAL4 ah;
      REAL4 cos1ah, sin1ah, cos2ah, sin2ah;
      REAL4 ai, bi;

      ah = alpha - DetStates->data[i].LMST;

      sin1ah = sinf ( ah );
      cos1ah = cosf ( ah );
      sin2ah = 2.0f * sin1ah * cos1ah;
      cos2ah = cos1ah * cos1ah - sin1ah * sin1ah;

      ai = a1 * cos2ah + a2 * sin2ah + a3 * cos1ah + a4 * sin1ah + a5;
      bi = b1 * cos2ah + b2 * sin2ah + b3 * cos1ah + b4 * sin1ah;
      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;

      /* sum A, B, C on the fly */
      coeffs->A += ai * ai;
      coeffs->B += bi * bi;
      coeffs->C += ai * bi;

    } /* for i < numSteps */

  /* finish calculation of A,B,C, D */
  norm = 2.0f / numSteps;
  coeffs->A *= norm;
  coeffs->B *= norm;
  coeffs->C *= norm;

  coeffs->D = coeffs->A * coeffs->B - coeffs->C * coeffs->C;

  RETURN(status);

} /* LALGetAMCoeffs() */


/** For a given vector of GPS-times, calculate the time-differences
 *  \f$\Delta T_\alpha\equiv T(t_\alpha) - T_0\f$, and their
 *  derivatives \f$Tdot_\alpha \equiv d T / d t (t_\alpha)\f$.
 * 
 *  \note The return-vectors \a DeltaT and \a Tdot must be allocated already
 *  and have the same length as the input time-series \a DetStates.
 *
 */
void
LALGetSSBtimes (LALStatus *status, 
		REAL8Vector *DeltaT,		/**< [out] DeltaT_alpha = T(t_alpha) - T_0*/
		REAL8Vector *Tdot,		/**< [out] Tdot(t_alpha) */
		const DetectorStateSeries *DetStates,/**< [in] detector-states at timestamps t_i */
		SkyPosition pos,		/**< source sky-location */
		REAL8 refTime,			/**< SSB reference-time T_0 of pulsar-parameters */
		SSBprecision precision		/**< relativistic or Newtonian SSB transformation? */
		)
{
  UINT4 numSteps, i;
  REAL8 vn[3];		/* unit-vector pointing to source in Cart. coord. */
  REAL8 alpha, delta;	/* source position */

  INITSTATUS( status, "LALGetSSBtimes", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (DetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  numSteps = DetStates->length;		/* number of timestamps */

  ASSERT (DeltaT, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (Tdot, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  ASSERT (DeltaT->length == numSteps, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT (Tdot->length == numSteps, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);

  ASSERT (precision < SSBPREC_LAST, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT ( pos.system == COORDINATESYSTEM_EQUATORIAL, status,
	   COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  
  /*----- get the cartesian source unit-vector */
  alpha = pos.longitude;
  delta = pos.latitude;
  vn[0] = cos(alpha) * cos(delta);
  vn[1] = sin(alpha) * cos(delta);
  vn[2] = sin(delta);

  /*----- now calculate the SSB transformation in the precision required */
  switch (precision)
    {
    case SSBPREC_NEWTONIAN:	/* use simple vr.vn to calculate time-delay */
      /*----- first figure out reference-time in SSB */

      for (i=0; i < numSteps; i++ )
	{
	  LIGOTimeGPS *ti = &(DetStates->data[i].tGPS);
	  /* DeltaT_alpha */
	  DeltaT->data[i]  = GPS2REAL8 ( (*ti) );
	  DeltaT->data[i] += SCALAR(vn, DetStates->data[i].rDetector);
	  DeltaT->data[i] -= refTime;

	  /* Tdot_alpha */
	  Tdot->data[i] = 1.0 + SCALAR(vn, DetStates->data[i].vDetector);
	  
	} /* for i < numSteps */

      break;

    case SSBPREC_RELATIVISTIC:	/* use LALBarycenter() to get SSB-times and derivative */
      for (i=0; i < numSteps; i++ )
	{
	  BarycenterInput baryinput = empty_BarycenterInput;
	  EmissionTime emit;
	  DetectorState *state = &(DetStates->data[i]);

	  baryinput.tgps = state->tGPS;
	  baryinput.site = DetStates->detector;
	  /* ARGHHH!!! */
	  baryinput.site.location[0] /= LAL_C_SI;
	  baryinput.site.location[1] /= LAL_C_SI;
	  baryinput.site.location[2] /= LAL_C_SI;

	  baryinput.alpha = alpha;
	  baryinput.delta = delta;
	  baryinput.dInv = 0;

	  TRY ( LALBarycenter(status->statusPtr, &emit, &baryinput, &(state->earthState)), status);

	  DeltaT->data[i] = GPS2REAL8 ( emit.te ) - refTime;
	  Tdot->data[i] = emit.tDot;

	} /* for i < numSteps */

      break;
    default:
      LALPrintError ("\n?? Something went wrong.. this should never be called!\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
      break;
    } /* switch precision */


  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetSSBtimes() */


/** Get the 'detector state' (ie position, velocity, etc) for the given
 * vector of timestamps.
 *
 * This function just calls LALBarycenterEarth() and LALBarycenter() for the
 * given vector of timestamps and returns the positions, velocities and LMSTs
 * of the detector, stored in a DetectorStateSeries. There is also an entry
 * containing the EarthState at each timestamp, which can be used as input for
 * subsequent calls to LALBarycenter().
 *
 * \note the DetectorStateSeries is allocated here and should be free'ed with
 * LALDestroyDetectorStateSeries().
 *
 */
void
LALGetDetectorStates (LALStatus *status,
		      DetectorStateSeries **DetStates,		/**< [out] series of DetectorStates */
		      const LIGOTimeGPSVector *timestamps,	/**< array of GPS timestamps t_i */
		      const LALDetector *detector,		/**< detector info */
		      const EphemerisData *edat			/**< ephemeris-files */		
		      )
{
  UINT4 i, j, numSteps;
  DetectorStateSeries *ret = NULL;

  INITSTATUS( status, "LALGetDetectorStates", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT ( DetStates, status,  COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *DetStates == NULL, status,  COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);

  ASSERT ( timestamps, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( detector, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( edat, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  numSteps = timestamps->length;

  TRY ( LALCreateDetectorStateSeries (status->statusPtr, &ret, numSteps), status);

  /* enter detector-info into the head of the state-vector */
  ret->detector = (*detector);
  
  /* now fill all the vector-entries corresponding to different timestamps */
  for ( i=0; i < numSteps; i++ )
    {
      BarycenterInput baryinput;
      EmissionTime emit;
      DetectorState *state = &(ret->data[i]);
      EarthState *earth = &(state->earthState);
      LIGOTimeGPS *tgps = &(timestamps->data[i]);

      /*----- first get earth-state */
      LALBarycenterEarth (status->statusPtr, earth, tgps, edat );
      BEGINFAIL(status){
	TRY ( LALDestroyDetectorStateSeries(status->statusPtr, &ret), status);
      }ENDFAIL(status);
      /*----- then get detector-specific info */
      baryinput.tgps = (*tgps);			/* irrelevant here! */
      baryinput.site = (*detector);
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;
      baryinput.alpha = baryinput.delta = 0;	/* irrelevant */
      baryinput.dInv = 0;

      LALBarycenter (status->statusPtr, &emit, &baryinput, earth);
      BEGINFAIL(status) {
	TRY ( LALDestroyDetectorStateSeries(status->statusPtr, &ret), status);
      }ENDFAIL(status);

      /*----- extract the output-data from this */
      for (j=0; j < 3; j++)	/* copy detector's position and velocity */
	{
	  state->rDetector[j] = emit.rDetector[j];
	  state->vDetector[j] = emit.vDetector[j];
	} /* for j < 3 */

      /* local mean sidereal time = GMST + longitude */
      state->LMST = earth->gmstRad + detector->frDetector.vertexLongitudeRadians;
      state->LMST = fmod (state->LMST, LAL_TWOPI );	/* normalize */

      /* enter timestamps */
      state->tGPS = (*tgps);

    } /* for i < numSteps */

  /* return result */
  (*DetStates) = ret;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALGetDetectorStates() */
 

/** Create a DetectorStateSeries */
void
LALCreateDetectorStateSeries (LALStatus *status, 
			      DetectorStateSeries **vect,	/**< output vector */
			      UINT4 length )			/**< number of entries */
{
  DetectorStateSeries *ret = NULL;

  INITSTATUS (status, "LALCreateDetectorStateSeries", rcsid );

  ASSERT ( vect, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *vect == NULL, status, COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);

  if ( (ret = LALCalloc(1, sizeof(DetectorStateSeries) )) == NULL ) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  if ( (ret->data = LALCalloc (length, sizeof(DetectorState) )) == NULL ) {
    LALFree (ret);
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  ret->length = length;

  /* return result */
  (*vect) = ret;

  RETURN (status);

} /* LALCreateDetectorStateSeries() */


/** Destroy a DetectorStateSeries (and set it to NULL) */
void
LALDestroyDetectorStateSeries (LALStatus *status, 
			       DetectorStateSeries **vect ) /**< pointer to vector to be destroyed */
{
  INITSTATUS (status, "LALDestroyDetectorStateSeries", rcsid );

  ASSERT ( vect, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  if ( *vect != NULL ) 
    {
      LALFree ( (*vect)->data );
      LALFree ( *vect );
      *vect = NULL;
    }

  RETURN (status);
} /* LALDestroyDetectorStateSeries() */
