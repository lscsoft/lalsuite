/*********************************************************************************/
/** \file ComputeFStatistic.c
 * Calculate the F-statistic for a given parameter-space of pulsar GW signals.
 * Implements the so-called "F-statistic" as introduced in JKS98.
 *                                                                          
 * \author Y. Ioth, C. Messenger, M.A. Papa, R.Prix, X. Siemens 
 *                                                                          
 *                 Albert Einstein Institute/UWM - started September 2002   
 *********************************************************************************/
#include "config.h"

/* System includes */
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/RngMedBias.h>
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/ComputeSkyBinary.h>
#include <lal/PulsarDataTypes.h>
#include <lal/SFTfileIO.h>

#include <lalapps.h>

/* local includes */
#include "clusters.h"
#include "DopplerScan.h"


RCSID( "$Id$");

#define BUFFERSIZE 1024                                                                   

/* Maximum fractional doppler shift */
#define DOPPLERMAX 1.e-4

/* (Half the ) number of terms to keep in the Dirichlet kernel sum */
#define NTERMS 32

#define MAXFILES 40000         /* Maximum # of files in a directory  */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  CHAR EphemEarth[MAXFILENAMELENGTH];	/**< filename of earth-ephemeris data */
  CHAR EphemSun[MAXFILENAMELENGTH];	/**< filename of sun-ephemeris data */
  UINT4 FreqImax;  		/**< number of frequency-bins to compute F-stat for */
  UINT4 SpinImax;		/**< number of spindown-bins to compute F for */
  REAL8 dFreq;			/**< frequency resolution */
  LIGOTimeGPS Tstart;		/**< start time of observation */
  REAL8 tSFT;			/**< length of an SFT in seconds */
  REAL8 tObs;			/**< total observation time in seconds */
  UINT4 SFTno;			/**< number of SFTs in input */
  UINT4 nsamples;		/**< number of frequency-bins in an SFT */
  LALDetector Detector;         /**< Our detector*/
  EphemerisData *edat;		/**< ephemeris data (from LALInitBarycenter()) */
  CHAR *skyRegion;		/**< sky-region to search (polygon defined by list of points) */
} ConfigVariables;

/** Local type for storing F-statistic output from NewLALDemod().
 * Note that length has to be the same for all vectors, anything else
 * has to be considered an error.
 */
typedef struct {
  UINT4 length;		/**< number of frequency-bins in the frequency-band */
  REAL8Vector *freqs;	/**< frequency-values of frequency-bins in Band */
  REAL8Vector *F;	/**< Values of the F-statistic proper (F) over frequency-bins */
  COMPLEX16Vector *Fa;	/**< Values Fa over frequency-bins */
  COMPLEX16Vector *Fb;	/**< Values Fb */
} FStatisticBand;

/** FIXME: OBSOLETE: holds a single FFT. Only kept for the moment to make things work (FIXME)*/
typedef struct FFTTag 
{
  COMPLEX8FrequencySeries *fft;   
} FFT;
/** FIXME: OBSOLETE: used to hold result from LALDemod(). Only kept for the moment to make things work (FIXME)*/
typedef struct {
  const REAL8         *F;            /* Array of value of the F statistic */
  const COMPLEX16     *Fa;           /* Results of match filter with a(t) */
  const COMPLEX16     *Fb;           /* Results of match filter with b(t) */
} LALFstat;


/** local parameter-type for local NewLALDemod() */
typedef struct {
  INT4		spinDwnOrder;	/* Maximum order of spdwn parameter */
  REAL8		*skyConst;	/* Constants computed in ComputeSky.c */
  REAL8		*spinDwn;	/* Spindown parameter set */
  AMCoeffs      *amcoe;         /*Amplitude Modulation coefficients */
  REAL8         f0;            	/*Starting Frequency to be demodulated*/
  REAL8         df;            	/*Frequency index resolution*/
  INT4          SFTno;          /* No. of SFTs*/
  INT4          Dterms;         /*Terms used in the computation of the dirichlet kernel*/
  INT4          ifmin;          /*smallest frequency index in SFTs */
  INT4          imax;           /*maximum # of values of F to calculate */
} NewDemodPar;


/* BINARY-MOD - structure to store a single binary signal template */  
typedef struct BinaryTemplatetag {             
    REAL8       ProjSMaxis;
    REAL8       Period;
    LIGOTimeGPS TperiSSB;
    REAL8       Eccentricity;
    REAL8       ArgPeri;
} BinaryTemplate;

typedef struct BinaryTemplateBanktag {
    REAL8          BTBfmax;
    REAL8          BTBTspan;
    LIGOTimeGPS    BTBTobsStart;
    INT4           BTBNFilters;
    REAL8          BTBMismatch;
    REAL8          BTBProjSMaxisMIN;
    REAL8          BTBProjSMaxisMAX;
    LIGOTimeGPS    BTBTperiSSBMIN;
    LIGOTimeGPS    BTBTperiSSBMAX;
    REAL8          BTBEccMIN;
    REAL8          BTBEccMAX;
    REAL8          BTBArgPeriMIN;
    REAL8          BTBArgPeriMAX;
    REAL8          BTBPeriodMIN;
    REAL8          BTBPeriodMAX;
    CHAR           BTBversion[256];
    BinaryTemplate *BTB;       
} BinaryTemplateBank;

/*----------------------------------------------------------------------*/
/* conditional compilation-switches */

/*----------------------------------------------------------------------*/
/* Error-codes */

#define COMPUTEFSTATC_ENULL 		1
#define COMPUTEFSTATC_ESYS     		2
#define COMPUTEFSTATC_EINPUT   		3
#define COMPUTEFSTATC_EMEM   		4
#define COMPUTEFSTATC_ENONULL 		5

#define COMPUTEFSTATC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATC_MSGESYS		"System call failed (probably file IO)"
#define COMPUTEFSTATC_MSGEINPUT   	"Invalid input"
#define COMPUTEFSTATC_MSGEMEM   	"Out of memory. Bad."
#define COMPUTEFSTATC_MSGENONULL 	"Output pointer is non-NULL"
/*----------------------------------------------------------------------
 * User-variables: can be set from config-file or command-line */

INT4 uvar_Dterms;
CHAR* uvar_IFO;
BOOLEAN uvar_SignalOnly;
BOOLEAN uvar_EstimSigParam;
BOOLEAN uvar_binary;  /* BINARY-MOD - added a flag to indicate binary search */
REAL8 uvar_dopplermax;  /* BINARY-MOD - added this because different binary systems will require different wings */
CHAR *uvar_binarytemplatefile;  /* BINARY-MOD - added for binary template file */
INT4 uvar_windowsize;          /* BINARY-MOD - added because of shorter SFT's */
REAL8 uvar_Freq;
REAL8 uvar_overSampling;
BOOLEAN uvar_useInterpolation;
REAL8 uvar_FreqBand;
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
CHAR *uvar_EphemDir;
CHAR *uvar_EphemYear;
INT4  uvar_gridType;
INT4  uvar_metricType;
REAL8 uvar_metricMismatch;
CHAR *uvar_skyRegion;
CHAR *uvar_DataDir;
CHAR *uvar_mergedSFTFile;
CHAR *uvar_BaseName;
BOOLEAN uvar_help;
CHAR *uvar_outputLabel;
CHAR *uvar_outputFstat;
CHAR *uvar_skyGridFile;
CHAR *uvar_outputSkyGrid;
CHAR *uvar_workingDir;
BOOLEAN uvar_searchNeighbors;
BOOLEAN uvar_openDX;

/*----------------------------------------------------------------------*/
/* some other global variables (FIXME) */
FFT **SFTData=NULL; 		/**< SFT Data for LALDemod (FIXME: to be removed) */
SFTVector *SFTvect = NULL;	/**< holds the SFT-data to analyze */
NewDemodPar *DemodParams  = NULL;	/**< Demodulation parameters for LALDemod */
LIGOTimeGPS *timestamps=NULL;	/**< Time stamps from SFT data */
LALFstat Fstat;			/**< output from LALDemod(): F-statistic and amplitudes Fa and Fb */
AMCoeffs amc;			/**< amplitude-modulation coefficients (and derived quantities) */
REAL8 Alpha,Delta;		/**< sky-position currently searched (equatorial coords, radians) */
Clusters HFLines;		/**< stores information about outliers/clusters in F-statistic */
Clusters HPLines;		/**< stores information about outliers/clusters in SFT-power spectrum */

BinaryTemplateBank *BinaryBank=NULL;    /* BINARY-MOD - pointer to structure to store binary template bank */
BinaryTemplate thisBinaryTemplate; /* BINARY-MOD - structure to store current binary template */
REAL8 bin_SMaxis;                 /* BINARY-MOD - Orbital semi-major axis of source system */
REAL8 bin_Period;                 /* BINARY-MOD - Orbital period of source system */
REAL8 bin_Eccentricity;           /* BINARY-MOD - Orbital eccentricity of source system */
REAL8 bin_ArgPeri;                /* BINARY-MOD - Argument of periapse of source system */
LIGOTimeGPS bin_TperiSSB;         /* BINARY-MOD - Time of periapse passage as defined in the SSB */
Clusters HFLines, HPLines;

Clusters *highSpLines=&HPLines, *highFLines=&HFLines;

REAL8 medianbias=1.0;		/**< bias in running-median depending on window-size (set in NormaliseSFTDataRngMdn()) */

FILE *fpmax;			/**< output-file: maximum of F-statistic over frequency-range */
FILE *fpstat;			/**< output-file: F-statistic candidates and cluster-information */

ConfigVariables GV;		/**< global container for various derived configuration settings */

/*----------------------------------------------------------------------*/
/* local prototypes */

int main(int argc,char *argv[]);
void initUserVars (LALStatus *stat);
INT4 ReadSFTData (void);
void InitFStat (LALStatus *status, ConfigVariables *cfg);
void CreateDemodParams (LALStatus *status);
void CreateBinaryDemodParams (LALStatus *status);
void CreateNautilusDetector (LALStatus *status, LALDetector *Detector);
void Freemem (LALStatus *status);
void EstimateFLines(LALStatus *status);
void NormaliseSFTDataRngMdn (LALStatus *status);
INT4 EstimateSignalParameters(INT4 * maxIndex);
INT4 writeFLines(INT4 *maxIndex);
int compare(const void *ip, const void *jp);
INT4 writeFaFb(INT4 *maxIndex);
void WriteFStatLog (LALStatus *stat, CHAR *argv[]);
int ReadBinaryTemplateBank(LALStatus *status);

void NewLALDemod(LALStatus *stat, FStatisticBand **FBand, FFT **input, NewDemodPar *params); 
void writeCOMPLEX16Vector(LALStatus *stat, COMPLEX16Vector *vect, const CHAR *fname);
const char *va(const char *format, ...);	/* little var-arg string helper function */

/*----------------------------------------------------------------------*/
/* some local defines */

#define EPHEM_YEARS  "00-04"

#define TRUE (1==1)
#define FALSE (1==0)

extern int vrbflg;

/*----------------------------------------------------------------------*/
/* empty structs for initialziations */

DopplerScanState emptyScan;


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
  CHAR Fmaxfilename[256]; 		/* Fmax file name*/
  UINT4 spdwn;				/* counter over spindown-params */
  DopplerScanInit scanInit;		/* init-structure for DopperScanner */
  DopplerScanState thisScan = emptyScan; /* current state of the Doppler-scan */
  DopplerPosition dopplerpos;		/* current search-parameters */
  SkyPosition thisPoint;
  FILE *fpOut=NULL;
  UINT4 loopcounter;
  FStatisticBand *FBand = NULL;		/* new type to store F-statistic results in a frequency-band */

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

  /* normalize SFTs by running median */
  LAL_CALL (NormaliseSFTDataRngMdn(&status), &status);

  /*   open file */
  strcpy(Fmaxfilename,"Fmax");
  if (uvar_outputLabel)
    strcat(Fmaxfilename,uvar_outputLabel);

  if (!(fpmax=fopen(Fmaxfilename,"w"))){
    fprintf(stderr,"in Main: unable to open Fmax file %s\n", Fmaxfilename);
    return 2;
  }

  /*      open file */
  strcpy(Fstatsfilename,"Fstats");
  if ( uvar_outputLabel )
    strcat(Fstatsfilename, uvar_outputLabel);

  if (!(fpstat=fopen(Fstatsfilename,"w"))){
    fprintf(stderr,"in Main: unable to open Fstats file\n");
    return 2;
  }

  if (!uvar_binary) {
    if (lalDebugLevel) LALPrintError ("\nSetting up template grid ...");
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
    scanInit.skyRegion = GV.skyRegion;
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
      if ( (fpOut = fopen (uvar_outputFstat, "w")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputFstat);
	  exit(-1);
	}
      if ( uvar_openDX )	/* prepend openDX header */
	{
	  UINT4 nFreq, nAlpha, nDelta;
	  nFreq = GV.FreqImax;
	  nAlpha = (UINT4)(uvar_AlphaBand / uvar_dAlpha) + 1;
	  nDelta = (UINT4)(uvar_DeltaBand / uvar_dDelta) + 1;

	  /* regular grid or not? */
	  if ( uvar_gridType == GRID_FLAT )
	    fprintf (fpOut, "grid = %d x %d x %d \n", nFreq, nAlpha, nDelta);
	  else
	    fprintf (fpOut, "points = %d \n", nFreq * thisScan.numGridPoints);

	  fprintf (fpOut, 
		   "format = ascii\n"
		   "field = locations, Fstat\n"
		   "structure = 3-vector, scalar\n"
		   "dependency = positions, positions\n"
		   "interleaving = field\n"
		   "end\n" );
	}
    } /* if outputFstat */

  } /* if (!uvar_binary) */

  /* BINARY-MOD - Call function to read in Binary template bank */
  if (uvar_binary) {
    LAL_CALL (ReadBinaryTemplateBank(&status), &status);  
    Alpha=uvar_Alpha;   /* BINARY-MOD - Also set sky location */
    Delta=uvar_Delta;
  }

  if (uvar_outputFstat)   
    {
      if ( (fpOut = fopen (uvar_outputFstat, "w")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputFstat);
	  exit(-1);
	}      
    } 

  if (lalDebugLevel) LALPrintError ("\nStarting main search-loop.. \n");

  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
   */
  loopcounter = 0;
  while (1)
    {
      /* BINARY-MOD - option to select next skyposition */
      if (!uvar_binary) 
	{
	  LAL_CALL (NextDopplerPos( &status, &dopplerpos, &thisScan ), &status);
	  /* Have we scanned all DopplerPositions yet? */
	  if (thisScan.state == STATE_FINISHED)
	    break;
	  LALNormalizeSkyPosition (&status, &thisPoint, &(dopplerpos.skypos) );
	  
	  Alpha = thisPoint.longitude;
	  Delta = thisPoint.latitude;
	  
	  LAL_CALL (CreateDemodParams(&status), &status);
	  
	} /* if (!uvar_binary) */

      /* BINARY-MOD - option to select next Binary template */
      if (uvar_binary) 
	{
	  
	  /* Have we scanned all DopplerPositions yet? */
	  if (((INT4)loopcounter) >= BinaryBank->BTBNFilters)
	    break;

	  printf("about to fill current template\n");
	  thisBinaryTemplate=(BinaryBank->BTB[loopcounter]); 
	  /*thisBinaryTemplate.ProjSMaxis=2.0;*/
	  /* thisBinaryTemplate.ProjSMaxis=(BinaryBank->BTB[counter]).ProjSMaxis; */
	  /*LAL_CALL (NextBinaryTemplate(&status, &BinaryBank ,&thisBinaryTemplate, &BTNumber ) , &status); */
	  /*printf("filled current template\n");
	  printf("the bank has smaxis as %le\n",BinaryBank->BTB[counter].ProjSMaxis);
	  printf("current smaxis is %le\n",thisBinaryTemplate.ProjSMaxis);
	  exit(0); */
	  printf("Current template is\n[%le,%le,%d,%d,%le,%le]\n",thisBinaryTemplate.ProjSMaxis, \
		 thisBinaryTemplate.Period,thisBinaryTemplate.TperiSSB.gpsSeconds, \
		 thisBinaryTemplate.TperiSSB.gpsNanoSeconds,thisBinaryTemplate.Eccentricity,  \
		 thisBinaryTemplate.ArgPeri);

	  /* printf("counter is %d Nfilters is %d\n",counter,BinaryBank->BTBNFilters); */
	  
	  printf("doing createbinarydemodparams\n");
	  LAL_CALL (CreateBinaryDemodParams(&status), &status);
	  printf("doing createbinarydemodparams\n");

	}  /* if (uvar_binary) */
      
      /* loop over spin params */
      for (spdwn=0; spdwn <= GV.SpinImax; spdwn++)
	{
	  /* BINARY-MOD - only deal with spin down for isolated at this stage */
	  if (!uvar_binary) 
	    DemodParams->spinDwn[0]=uvar_f1dot + spdwn*uvar_df1dot;

	  LAL_CALL (NewLALDemod (&status, &FBand, SFTData , DemodParams), &status);

	  /* FIXME: TMP very specific testing Fa/Fb refining by interpolation */
	  if (loopcounter == 0)
	    {
	      COMPLEX16Vector *FaRefined = NULL;
	      COMPLEX16Vector *FbRefined = NULL;
	      REAL8Vector *FstatRefined = NULL;
	      UINT4 nbins, ni; 
	      REAL4 A, B, C, D;
	      REAL8 FaRe, FaIm, FbRe, FbIm;
	      UINT4 M; 

	      if ( uvar_useInterpolation && (uvar_overSampling != 1) )
		{
		  nbins =  uvar_overSampling * FBand->Fa->length;
		  LAL_CALL (refineCOMPLEX16Vector(&status, &FaRefined, FBand->Fa, nbins, uvar_Dterms), &status);
		  LAL_CALL (refineCOMPLEX16Vector(&status, &FbRefined, FBand->Fb, nbins, uvar_Dterms), &status);

		  /* calculate F-statistic by combining Fa and Fb properly */
		  LAL_CALL (LALDCreateVector (&status, &FstatRefined, nbins), &status);

		  A = DemodParams->amcoe->A;
		  B = DemodParams->amcoe->B;
		  C = DemodParams->amcoe->C;
		  D = DemodParams->amcoe->D;
		  
		  M = DemodParams->SFTno;

		  for (ni = 0; ni < nbins; ni++)
		    {
		      FaRe = FaRefined->data[ni].re;
		      FaIm = FaRefined->data[ni].im;
		      FbRe = FbRefined->data[ni].re;
		      FbIm = FbRefined->data[ni].im;

		      FstatRefined->data[ni] = B * (FaRe*FaRe + FaIm*FaIm) + A * (FbRe*FbRe + FbIm*FbIm) - 2.0*C*(FaRe*FbRe + FaIm*FbIm);
		      FstatRefined->data[ni] *= (4.0/(M*D));

		    } /* for ni < nbins */

		  if (lalDebugLevel) 
		    {
		      LAL_CALL (writeCOMPLEX16Vector(&status, FaRefined, va("Fa_interpolated_%f.dat", uvar_overSampling) ), &status);
		      LAL_CALL (writeCOMPLEX16Vector(&status, FBand->Fa, "Fa_natural.dat"), &status);
		    }

		  LAL_CALL (LALZDestroyVector (&status, &FaRefined), &status);
		  LAL_CALL (LALZDestroyVector (&status, &FbRefined), &status);
		  LAL_CALL (LALDDestroyVector (&status, &FstatRefined), &status);

		} /* if useInterpolation */

	      if (lalDebugLevel && !uvar_useInterpolation) {
		LAL_CALL (writeCOMPLEX16Vector(&status, FBand->Fa, va("Fa_oversampled_%f.dat", uvar_overSampling) ), &status);
	      }

	    } /* loopcounter == 0 (interpolation-testing) */

	  /* FIXME: to keep cluster-stuff working, we provide the "translation" from FBand back into old Fstats-struct */
	  Fstat.F = FBand->F->data;
	  Fstat.Fa = FBand->Fa->data;
	  Fstat.Fb = FBand->Fb->data;

	  /*  This fills-in highFLines */
	  if (GV.FreqImax > 5) {
	    LAL_CALL (EstimateFLines(&status), &status);
	  }
	  
	  /* now, if user requested it, we output ALL F-statistic results */
	  if ((fpOut)&&(!uvar_binary)) 
	    {
	      UINT4 i;
	      printf("outputting isolated Fstat\n");
	      for(i=0;i < GV.FreqImax ;i++)
		{
		  fprintf (fpOut, "%20.17f %20.17f %20.17f %20.17f\n", 
			   uvar_Freq + i*GV.dFreq, Alpha, Delta, 2.0*medianbias*Fstat.F[i]);
		}

	    } /* if outputFstat and not binary */

	  /* BINARY-MOD - output binary search F-stat with binary params */
	  if ((fpOut)&&(uvar_binary)) 
	    {
	      UINT4 i;
	      printf("printing binary output\n");
	      for(i=0;i < GV.FreqImax ;i++)
		{
		  fprintf (fpOut, "%6.12f %6.12f %12.12f %d %d %6.12f %6.12f %12.12f\n", 
			   uvar_Freq + i*GV.dFreq, thisBinaryTemplate.ProjSMaxis, \
			   thisBinaryTemplate.Period, thisBinaryTemplate.TperiSSB.gpsSeconds, \
			   thisBinaryTemplate.TperiSSB.gpsNanoSeconds, \
			   thisBinaryTemplate.Eccentricity, thisBinaryTemplate.ArgPeri, \
			   2.0*medianbias*Fstat.F[i]);
		}

	    } /* if outputFstat and binary */

	  
	  /*  This fills-in highFLines  */
	  if (highFLines != NULL && highFLines->Nclusters > 0)
	    {
	      maxIndex=(INT4 *)LALMalloc(highFLines->Nclusters*sizeof(INT4));
	    
	      /*  for every cluster writes the information about it in file Fstats */
	      if (writeFLines(maxIndex)){
		fprintf(stderr, "%s: trouble making file Fstats\n", argv[0]);
		return 6;
	      }
	   	  
	      if (uvar_EstimSigParam)
		{
		  if(writeFaFb(maxIndex)) return 255;
		  if (EstimateSignalParameters(maxIndex)) return 7;
		} /* if signal-estimation */
	  
	  
	      LALFree(maxIndex);
	    } /* if highFLines found */

	  
	  /* Set the number of the clusters detected to 0 at each iteration 
	     of the sky-direction and the spin down */
	  highFLines->Nclusters=0;

	} /* For GV.spinImax */


      /* free Fstats-results from this skypos */
      LAL_CALL (LALDDestroyVector (&status, &(FBand->F)), &status);
      LAL_CALL (LALDDestroyVector (&status, &(FBand->freqs)), &status);
      LAL_CALL (LALZDestroyVector (&status, &(FBand->Fa)), &status);
      LAL_CALL (LALZDestroyVector (&status, &(FBand->Fb)), &status);
      LALFree (FBand);
      FBand = NULL;

      loopcounter ++;
      if (lalDebugLevel) LALPrintError ("Search progress: %5.1f%%", 
					(100.0* loopcounter / thisScan.numGridPoints));
    } /*  while SkyPos */

  if (uvar_outputFstat && fpOut)
    fclose (fpOut);

  if (lalDebugLevel) LALPrintError ("\nSearch finished.\n");

  fclose(fpmax);

  /* properly terminate Fstats-file by '.*DONE.*' marker: 
   * we use "%DONE" for matlab/octave compatibility */
  fprintf(fpstat, "%%DONE");
  fclose(fpstat);

  /* Free DopplerScan-stuff (grid) */
  if (!uvar_binary)  {
      LAL_CALL ( FreeDopplerScan(&status, &thisScan), &status);
  }


  LAL_CALL ( Freemem(&status), &status);

  return 0;

} /* main() */


/** 
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
void
initUserVars (LALStatus *stat)
{
  INITSTATUS( stat, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (stat);

  /* set a few defaults */
  uvar_Dterms 	= 16;
  uvar_FreqBand = 0.0;
  uvar_overSampling = 2.0;
  uvar_useInterpolation = FALSE;	/* use LALDemod() for oversampling by default */

  uvar_binary=FALSE;       /* BINARY-MOD - Set binary flag to FALSE for safety */

  uvar_Alpha 	= 0.0;
  uvar_Delta 	= 0.0;
  uvar_AlphaBand = 0;
  uvar_DeltaBand = 0;
  uvar_dAlpha 	= 0.001;
  uvar_dDelta 	= 0.001;
  uvar_skyRegion = NULL;

  uvar_dopplermax=1e-4;       /* BINARY-MOD - Need to have this as an option (important for binary) */
  uvar_windowsize=50;      /* BINARY-MOD - Need to have this as an option (important for shorter SFT's ?) */

  uvar_EphemYear = LALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar_EphemYear, EPHEM_YEARS);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar_EphemDir = LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar_EphemDir, DEFAULT_EPHEMDIR);

  uvar_SignalOnly = FALSE;
  uvar_EstimSigParam = FALSE;
 
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
  uvar_openDX = FALSE;

  uvar_skyGridFile = NULL;

  uvar_workingDir = LALMalloc(512);
  strcpy(uvar_workingDir, ".");

  uvar_searchNeighbors = FALSE;



  /* register all our user-variables */
  LALregBOOLUserVar(stat, 	help, 		'h', UVAR_HELP,     "Print this message"); 
  LALregINTUserVar(stat,	Dterms,		't', UVAR_OPTIONAL, "Number of terms to keep in Dirichlet kernel sum");
  LALregREALUserVar(stat, 	Freq, 		'f', UVAR_REQUIRED, "Starting search frequency in Hz");
  LALregREALUserVar(stat, 	FreqBand, 	'b', UVAR_OPTIONAL, "Search frequency band in Hz");
  LALregREALUserVar(stat, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserVar(stat, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
  LALregREALUserVar(stat, 	AlphaBand, 	'z', UVAR_OPTIONAL, "Band in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(stat, 	DeltaBand, 	'c', UVAR_OPTIONAL, "Band in delta (equatorial coordinates) in radians");
  LALregREALUserVar(stat, 	dAlpha, 	'l', UVAR_OPTIONAL, "Resolution in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(stat, 	dDelta, 	'g', UVAR_OPTIONAL, "Resolution in delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(stat,	DataDir, 	'D', UVAR_OPTIONAL, "Directory where SFT's are located");
  LALregSTRINGUserVar(stat,	mergedSFTFile, 	'B', UVAR_OPTIONAL, "Merged SFT's file to be used"); 
  LALregSTRINGUserVar(stat,	BaseName, 	'i', UVAR_OPTIONAL, "The base name of the input  file you want to read");
  LALregSTRINGUserVar(stat,	EphemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(stat,	EphemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregSTRINGUserVar(stat, 	IFO, 		'I', UVAR_REQUIRED, "Detector: GEO(0), LLO(1), LHO(2), NAUTILUS(3), VIRGO(4), TAMA(5), CIT(6)");
  LALregBOOLUserVar(stat, 	SignalOnly, 	'S', UVAR_OPTIONAL, "Signal only flag");
  LALregBOOLUserVar(stat, 	binary, 	'u', UVAR_OPTIONAL, "Binary search flag");
  LALregREALUserVar(stat, 	dopplermax, 	'q', UVAR_OPTIONAL, "Maximum doppler shift expected");  
  LALregREALUserVar(stat, 	f1dot, 		's', UVAR_OPTIONAL, "First spindown parameter f1dot");
  LALregREALUserVar(stat, 	f1dotBand, 	'm', UVAR_OPTIONAL, "Search-band for f1dot");
  LALregREALUserVar(stat, 	df1dot, 	'e', UVAR_OPTIONAL, "Resolution for f1dot (default 1/(2*Tobs*tSFT*Nsft)");
  LALregBOOLUserVar(stat, 	EstimSigParam, 	'p', UVAR_OPTIONAL, "Do Signal Parameter Estimation");
  LALregREALUserVar(stat, 	Fthreshold,	'F', UVAR_OPTIONAL, "Signal Set the threshold for selection of 2F");
  LALregINTUserVar(stat, 	windowsize,	'k', UVAR_OPTIONAL, "Running-Median window size");
  LALregINTUserVar(stat, 	gridType,	 0 , UVAR_OPTIONAL, "Template grid: 0=flat, 1=isotropic, 2=metric, 3=file");
  LALregINTUserVar(stat, 	metricType,	'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserVar(stat, 	metricMismatch,	'X', UVAR_OPTIONAL, "Maximal mismatch for metric tiling");
  LALregSTRINGUserVar(stat,	skyRegion, 	'R', UVAR_OPTIONAL, "Specify sky-region by polygon");
  LALregSTRINGUserVar(stat,	outputLabel,	'o', UVAR_OPTIONAL, "Label to be appended to all output file-names");
  LALregSTRINGUserVar(stat,	outputFstat,	 0,  UVAR_OPTIONAL, "Output-file for the F-statistic field over the parameter-space");
  LALregSTRINGUserVar(stat,	skyGridFile,	 0,  UVAR_OPTIONAL, "Load sky-grid from this file.");
  LALregSTRINGUserVar(stat,	outputSkyGrid,	 0,  UVAR_OPTIONAL, "Write sky-grid into this file.");
  LALregSTRINGUserVar(stat,	binarytemplatefile,	 0,  UVAR_OPTIONAL, "Read binary templates from this file.");

  LALregBOOLUserVar(stat,	openDX,	 	 0,  UVAR_OPTIONAL, "Make output-files openDX-readable (adds proper header)");
  LALregSTRINGUserVar(stat,     workingDir,     'w', UVAR_OPTIONAL, "Directory to be made the working directory, . is default");

  LALregREALUserVar(stat, 	overSampling,	'r', UVAR_OPTIONAL, "Oversampling factor for frequency resolution.");
  LALregBOOLUserVar(stat,	useInterpolation, 0, UVAR_OPTIONAL, "Use Interpolation instead of LALDemod() for the oversampling.");


  DETATCHSTATUSPTR (stat);
  RETURN (stat);
} /* initUserVars() */



/*******************************************************************************/
/*  Note that there is a degeneracy where the shifts taken at the same time,  */
/*  psi -> psi+Pi/2 and Phi0 -> Phi0 + Pi/2,  */
/*  give the same A1,A2,A3,A4.  */
/*  */
/*******************************************************************************/


int EstimateSignalParameters(INT4 * maxIndex)
{
  INT4 irec,jrec;
  REAL8 A1,A2,A3,A4,A=amc.A,B=amc.B,C=amc.C,D=amc.D;
  REAL8 beta,Asq,detA,ampratio,A1test;
  REAL8 psi_mle,Phi0_mle,mu_mle;
  REAL8 h0mle,h0mleSq;
  REAL8 error_tol=1.0/pow(10,14);
  REAL8 norm;
  FILE * fpMLEParam;
  CHAR Paramfilename[256];

  strcpy(Paramfilename,"ParamMLE");
  if (uvar_outputLabel)
    strcat(Paramfilename,uvar_outputLabel);
  
  if(!(fpMLEParam=fopen(Paramfilename,"w")))
    fprintf(stderr,"Error in EstimateSignalParameters: unable to open the file");


  norm=2.0*sqrt(GV.tSFT)/(GV.tSFT*GV.SFTno);


  for(jrec=0;jrec < highFLines->Nclusters ;jrec++)
    {

      irec=maxIndex[jrec];

      A1 =  2.0*( B * Fstat.Fa[irec].re - C * Fstat.Fb[irec].re) / D;
      A2 =  2.0*( A * Fstat.Fb[irec].re - C * Fstat.Fa[irec].re) / D;
      A3 = - 2.0*( B * Fstat.Fa[irec].im - C * Fstat.Fb[irec].im) / D;
      A4 = - 2.0*( A * Fstat.Fb[irec].im - C * Fstat.Fa[irec].im) / D;



      Asq = A1*A1 + A2*A2 + A3*A3 + A4*A4;
      detA = A1*A4-A2*A3;

      h0mle = 0.5*pow(
		      pow(((A1-A4)*(A1-A4)+(A2+A3)*(A2+A3)),0.25)+
		      pow(((A1+A4)*(A1+A4)+(A2-A3)*(A2-A3)),0.25)
		      ,2);

      h0mleSq = pow(h0mle,2.0);
      ampratio=Asq/h0mleSq;


      if(ampratio<0.25-error_tol||ampratio>2.0+error_tol) 
	{
	  fprintf(stderr,"Imaginary Cos[iota]; cannot compute parameters");
	  fprintf(stderr,"in the EstimateSignalParameters routine");
	  fprintf(stderr,"in ComputeFStatistic code");
	  fprintf(stderr,"Now exitting...");
	  /* 	  break; */
	  exit(1);
	}

      if(fabs(ampratio-0.25)<error_tol) {
	mu_mle =0.0;
      } else if(fabs(ampratio-2.0)<error_tol) {
	mu_mle = 1.0;
      } else {
	mu_mle = sqrt(-3.0+2.0*sqrt(2.0+ampratio));
      }

      if(detA<0) 
	mu_mle = - 1.0*mu_mle;


      if(Asq*Asq < 4.0*detA*detA)
	{
	  fprintf(stderr,"Imaginary beta; cannot compute parameters");
	  break;
	}

      /* Compute MLEs of psi and Phi0 up to sign of Cos[2*Phi0] */
      /* Make psi and Phi0 always in -Pi/2 to Pi/2 */ 
      beta  = (Asq + sqrt(Asq*Asq - 4.0*detA*detA))/(2.0*detA);
      psi_mle  = atan( (beta*A4-A1)/(beta*A3+A2) )/2.0;
      Phi0_mle  = atan( (A1-beta*A4)/(A3+beta*A2) )/2.0;


      A1test=h0mle*(0.5*(1+mu_mle*mu_mle)*cos(2.0*psi_mle)*cos(2.0*Phi0_mle)
		    -mu_mle*sin(2.0*psi_mle)*sin(2.0*Phi0_mle));

      /* Determine the sign of Cos[2*Phi0] */
      if(A1*A1test<0) {
	if(Phi0_mle>0) {
	  Phi0_mle=Phi0_mle - LAL_PI/2.0;
	} else {
	  Phi0_mle=Phi0_mle + LAL_PI/2.0;
	}
      }

      /* normalization */
      h0mle=h0mle*norm;


      /* For the real data, we need to multiply long(2.0) */
      /* Because we use running median to estimate the S_h. */
      /* if(GV.SignalOnly!=1) 
	h0mle=h0mle*sqrt(medianbias);
      */
      /* medianbias is 1 when GV.SignalOnly==1 */
      h0mle=h0mle*sqrt(medianbias);

      /* Note that we print out MLE of 2.0*Phi0_JKS */
      /* because Phi0_PULGROUPDOC=2.0*Phi0_JKS */
      /* and Phi0_PULGROUPDOC is the one used in In.data. */
 
      /* medianbias is 1 if GV.SignalOnly==1 */
      fprintf(fpMLEParam,"%16.8lf %22E", uvar_Freq + irec*GV.dFreq, 2.0*medianbias*Fstat.F[irec]);


      fprintf(fpMLEParam,"  %10.6f",(1.0+mu_mle*mu_mle)*h0mle/2.0);
      fprintf(fpMLEParam,"  %10.6f",mu_mle*h0mle);
      fprintf(fpMLEParam,"  %10.6f",psi_mle);
      fprintf(fpMLEParam,"  %10.6f",2.0*Phi0_mle);
      fprintf(fpMLEParam,"\n");
    }

  fclose(fpMLEParam);

  return 0;

} /* EstimateSignalParameters() */


/*******************************************************************************/

/* Write the Fa and Fb for the later use of Fstatistic Shape test */
/* the explicit format specifier like %22.12f looks ugly. */
int writeFaFb(INT4 *maxIndex)
{
  INT4 irec,jrec;
  INT4 index,krec=0;
  CHAR filename[256];         /* Base of the output file name */
  CHAR noiseswitch[16];
  CHAR clusterno[16];
  INT4 N;
  FILE * fp=NULL;
  REAL8 bias=1.0;
  CHAR FaFbfilename[256];
  
  strcpy(FaFbfilename,"FaFb");
  if (uvar_outputLabel)
    strcat(FaFbfilename,uvar_outputLabel);
  sprintf(noiseswitch,"%02d", uvar_SignalOnly);
  strcat(FaFbfilename,noiseswitch);

  bias=sqrt(medianbias);


  for (irec=0;irec<highFLines->Nclusters;irec++){
    sprintf(clusterno,".%03d",irec+1);
    strcpy(filename,FaFbfilename);
    strcat(filename,clusterno);

    if((fp=fopen(filename,"w"))==NULL) {
      fprintf(stderr,"Unable to open a file %s\n",filename);
      return 1;
    }

    N=highFLines->NclustPoints[irec];


    /* The header contains */
    /* N the number of points in the cluster */
    /* the frequency where the maximum amplitude is */
    /* A,B,C coefficients */
    index=highFLines->Iclust[krec];

    fprintf(fp,"%10d\n",N);
    fprintf(fp,"%22.12f %22.12f\n",
	    uvar_Freq+maxIndex[irec]*GV.dFreq,
	    Fstat.F[maxIndex[irec]]*bias*bias);
    fprintf(fp,"%22.12f %22.12f\n",uvar_Freq + index * GV.dFreq, GV.dFreq);
    fprintf(fp,"%22.12f %22.12f %22.12f\n",amc.A,amc.B,amc.C);



    /* WARNING 1: 
       Here we make it sure that Fa and Fb are already 
       normaized by M (the number of sfts). 
       See the pulsargroup document or LALDemod document.

       This means that we use in the FstatShapeTest code  
       F = (4.0/D)*(B*FaSq + A*FbSq - 2.0*C*FaFb); 
       instead of 
       F = (4.0/(M*D))*(B*FaSq + A*FbSq - 2.0*C*FaFb); 
    */

    /* WARNING 2: 
       Here we assume that the one sided noise power spectrum density 
       Sh is properly normalized. Namely, the Fa and Fb must be 
       multiplied by B*log(2.0) if one uses running median for an 
       estimate of Sh. B is the sample bias of the sample median 
       estimate.
    */


    /* WARNING 3:
       The information stored in the FaFb file will be used by 
       FstatShapeTest code. The format of the file must not be 
       changed without the appropriate change in FstatShapeTest 
       code. 
     */ 

    /* BINARY-MOD - need to check with FstatShapeTest to get input format */

    for(jrec=0;jrec<N;jrec++) {
      index=highFLines->Iclust[krec];
      krec++;

      /* Freqency, Re[Fa],Im[Fa],Re[Fb],Im[Fb], F */
      fprintf(fp,"%22.16f %22.12f %22.12f %22.12f %22.12f %22.12f\n",
	      uvar_Freq+index*GV.dFreq,
	      Fstat.Fa[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fa[index].im/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].im/sqrt(GV.SFTno)*bias,
	      Fstat.F[index]*bias*bias);

    }
    fclose(fp);
  }

  return 0;
}



/*******************************************************************************/

void CreateDemodParams (LALStatus *status)
{
  CSParams *csParams  = NULL;        /* ComputeSky parameters */
  AMCoeffsParams *amParams;
  EarthState earth;
  EmissionTime emit;
  LIGOTimeGPS *midTS=NULL;           /* Time stamps for amplitude modulation coefficients */
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  UINT4 k;

  INITSTATUS (status, "CreateDemodParams", rcsid);
  ATTATCHSTATUSPTR (status);
  
  /* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.site.location[0]=GV.Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=GV.Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=GV.Detector.location[2]/LAL_C_SI;
  baryinput.alpha=Alpha;
  baryinput.delta=Delta;
  baryinput.dInv=0.e0;

/* amParams structure to compute a(t) and b(t) */

/* Allocate space for amParams stucture */
/* Here, amParams->das is the Detector and Source info */
  amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
/* Fill up AMCoeffsParams structure */
  amParams->baryinput = &baryinput;
  amParams->earth = &earth; 
  amParams->edat = GV.edat;
  amParams->das->pDetector = &GV.Detector; 
  amParams->das->pSource->equatorialCoords.latitude = Delta;
  amParams->das->pSource->equatorialCoords.longitude = Alpha;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!!!!!!!!*/
  amParams->leapAcc = LALLEAPSEC_STRICT;

 /* Mid point of each SFT */
   midTS = (LIGOTimeGPS *)LALCalloc(GV.SFTno, sizeof(LIGOTimeGPS));
   for(k=0; k<GV.SFTno; k++)
     { 
       REAL8 teemp=0.0;
       TRY (LALGPStoFloat(status->statusPtr, &teemp, &(timestamps[k])), status);
       teemp += 0.5*GV.tSFT;
       TRY (LALFloatToGPS(status->statusPtr, &(midTS[k]), &teemp), status);
     }
   
   TRY (LALComputeAM(status->statusPtr, &amc, midTS, amParams), status); 

/* ComputeSky stuff*/
  csParams=(CSParams *)LALMalloc(sizeof(CSParams));
  csParams->tGPS=timestamps;  
  csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
  csParams->mObsSFT=GV.SFTno;     /* Changed this from GV.mobssft !!!!!! */
  csParams->tSFT=GV.tSFT;
  csParams->edat = GV.edat;
  csParams->baryinput=&baryinput;
  csParams->spinDwnOrder=1;
  csParams->skyPos[0]=Alpha;
  csParams->skyPos[1]=Delta;
  csParams->earth = &earth;
  csParams->emit = &emit;

/* Finally, DemodParams */
  DemodParams->amcoe=&amc;
  DemodParams->spinDwnOrder=1;
  DemodParams->SFTno=GV.SFTno;

  DemodParams->f0=uvar_Freq;
  DemodParams->imax=GV.FreqImax;
  DemodParams->df=GV.dFreq;

  DemodParams->Dterms=uvar_Dterms;
  DemodParams->ifmin=(INT4) (SFTvect->data[0].f0 / SFTvect->data[0].deltaF + 0.5);

  /* compute the "sky-constants" A and B */
  TRY ( ComputeSky (status->statusPtr, DemodParams->skyConst, 0, csParams), status);  
  LALFree(midTS);

  LALFree(csParams->skyPos);
  LALFree(csParams);

  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* CreateDemodParams() */

/*******************************************************************************/

/* BINARY-MOD - I decided to add this function rather than change the existing one */
void CreateBinaryDemodParams (LALStatus *status)  
{
  CSBParams *csbParams  = NULL;        /* ComputeSkyBinary parameters */
  AMCoeffsParams *amParams;
  EarthState earth;
  EmissionTime emit;
  LIGOTimeGPS *midTS=NULL;           /* Time stamps for amplitude modulation coefficients */
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  UINT4 k;

  INITSTATUS (status, "CreateBinaryDemodParams", rcsid);
  ATTATCHSTATUSPTR (status);
  
  printf("alpha is %le delta is %le\n",Alpha,Delta);

  /* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.site.location[0]=GV.Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=GV.Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=GV.Detector.location[2]/LAL_C_SI;
  baryinput.alpha=Alpha;
  baryinput.delta=Delta;
  baryinput.dInv=0.e0;

/* amParams structure to compute a(t) and b(t) */

/* Allocate space for amParams stucture */
/* Here, amParams->das is the Detector and Source info */
  amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
/* Fill up AMCoeffsParams structure */
  amParams->baryinput = &baryinput;
  amParams->earth = &earth; 
  amParams->edat = GV.edat;
  amParams->das->pDetector = &GV.Detector; 
  amParams->das->pSource->equatorialCoords.latitude = Delta;
  amParams->das->pSource->equatorialCoords.longitude = Alpha;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!!!!!!!!*/
  amParams->leapAcc = LALLEAPSEC_STRICT;

 /* Mid point of each SFT */
   midTS = (LIGOTimeGPS *)LALCalloc(GV.SFTno, sizeof(LIGOTimeGPS));
   for(k=0; k<GV.SFTno; k++)
     { 
       REAL8 teemp=0.0;
       TRY (LALGPStoFloat(status->statusPtr, &teemp, &(timestamps[k])), status);
       teemp += 0.5*GV.tSFT;
       TRY (LALFloatToGPS(status->statusPtr, &(midTS[k]), &teemp), status);
     }
   
   TRY (LALComputeAM(status->statusPtr, &amc, midTS, amParams), status); 

/* ComputeSkyBinary stuff*/
  csbParams=(CSBParams *)LALMalloc(sizeof(CSBParams));
  csbParams->tGPS=timestamps;  
  csbParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
  csbParams->mObsSFT=GV.SFTno;     /* Changed this from GV.mobssft !!!!!! */
  csbParams->tSFT=GV.tSFT;
  csbParams->edat = GV.edat;
  csbParams->baryinput=&baryinput;
  csbParams->spinDwnOrder=0;    /* need to set up an input variable for this */
  csbParams->skyPos[0]=Alpha;
  csbParams->skyPos[1]=Delta;
  csbParams->earth = &earth;
  csbParams->emit = &emit;
  csbParams->SemiMajorAxis=thisBinaryTemplate.ProjSMaxis;     /* Added the orbital parameter inputs */
  csbParams->OrbitalPeriod=thisBinaryTemplate.Period;
  csbParams->TperiapseSSB.gpsSeconds=thisBinaryTemplate.TperiSSB.gpsSeconds;
  csbParams->TperiapseSSB.gpsNanoSeconds=thisBinaryTemplate.TperiSSB.gpsNanoSeconds;
  csbParams->ArgPeriapse=thisBinaryTemplate.ArgPeri;
  csbParams->OrbitalEccentricity=thisBinaryTemplate.Eccentricity;

/* Finally, DemodParams */
  DemodParams->amcoe=&amc;
  DemodParams->spinDwnOrder=0;
  DemodParams->SFTno=GV.SFTno;

  DemodParams->f0=uvar_Freq;
  DemodParams->imax=GV.FreqImax;
  DemodParams->df=GV.dFreq;

  DemodParams->Dterms=uvar_Dterms;
  
  DemodParams->ifmin=(INT4) (SFTvect->data[0].f0 / SFTvect->data[0].deltaF + 0.5);

  /* compute the "sky-constants" A and B */
  TRY ( ComputeSkyBinary (status->statusPtr, DemodParams->skyConst, 0, csbParams), status);  
  LALFree(midTS);

  LALFree(csbParams->skyPos);
  LALFree(csbParams);

  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* CreateBinaryDemodParams() */

/*******************************************************************************/

/*  for every cluster writes the information about it in file Fstats */
/*  precisely it writes: */
/*  fr_max alpha delta N_points_of_cluster mean std max (of 2F) */
int writeFLines(INT4 *maxIndex){

  INT4 i,j,j1,j2,k,N;
  REAL8 max,log2,mean,var,std,R,fr;
  INT4 imax;
  INT4 err;

  log2=medianbias;
 
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
      R=2.0*log2*highFLines->clusters[j1];
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
      R=2.0*log2*highFLines->clusters[j2];
      j2=j2+1;
      var=var+(R-mean)*(R-mean);
    }/*  end j loop over points of i-th cluster  */
    var=var/N;
    std=sqrt(var);
    fr=uvar_Freq + imax*GV.dFreq;
/*    print the output */
    if (!uvar_binary)
      {
	err=fprintf(fpstat,"%16.12f %10.8f %10.8f    %d %10.5f %10.5f %10.5f\n",fr,
		    Alpha, Delta, N, mean, std, max);
      }
    if (uvar_binary)  /* BINARY-MOD - changed to also output biniary search parameters */
      {
	err=fprintf(fpstat,"%16.12f %10.8f %10.8f %10.8f %10.8f %d %d %10.8f %10.8f %d %10.5f %10.5f %10.5f\n",
		    fr, Alpha, Delta, thisBinaryTemplate.ProjSMaxis, thisBinaryTemplate.Period, 
		    thisBinaryTemplate.TperiSSB.gpsSeconds, thisBinaryTemplate.TperiSSB.gpsNanoSeconds,
		    thisBinaryTemplate.Eccentricity, thisBinaryTemplate.ArgPeri, 
		    N, mean, std, max);
      }
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
		      "Try ./ComputeFStatistic -h \n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);      
    }

  if(uvar_DataDir && uvar_mergedSFTFile)
    {
      LALPrintError ( "\nCannot specify both 'DataDir' and 'mergedSFTfile'.\n"
		      "Try ./ComputeFStatistic -h \n\n" );
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }      

  if (uvar_EphemYear == NULL)
    {
      LALPrintError ("\nNo ephemeris year specified (option 'EphemYear')\n\n");
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

  /* extract timestamps-vector from SFT-data */
  /* and prepare SFT-data in FFT-struct to feed it into LALDemod() (FIXME)*/
  if ( (timestamps = LALCalloc(1, cfg->SFTno * sizeof(LIGOTimeGPS))) == NULL) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  if ( (SFTData = LALCalloc(1, cfg->SFTno * sizeof(FFT*))) == NULL) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }


  for (i=0; i < cfg->SFTno; i++)
    {
      timestamps[i] = SFTvect->data[i].epoch;
      SFTData[i] = LALCalloc(1, sizeof(FFT));
      SFTData[i]->fft = &(SFTvect->data[i]);
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

  if (LALUserVarWasSet (&uvar_EphemDir) )
    {
      sprintf(cfg->EphemEarth, "%s/earth%s.dat", uvar_EphemDir, uvar_EphemYear);
      sprintf(cfg->EphemSun, "%s/sun%s.dat", uvar_EphemDir, uvar_EphemYear);
    }
  else
    {
      sprintf(cfg->EphemEarth, "earth%s.dat", uvar_EphemYear);
      sprintf(cfg->EphemSun, "sun%s.dat",  uvar_EphemYear);
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
      LALPrintError ("\nUnknown detector. Currently allowed are 'GEO', 'LLO', 'LHO', 'NAUTILUS', 'VIRGO', 'TAMA', 'CIT' or '0'-'6'\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }

  /*----------------------------------------------------------------------
   * set some defaults
   */
 if (!uvar_binary)
   {
     cfg->dFreq = 1.0 / cfg->tObs;
     if (! uvar_useInterpolation)
       cfg->dFreq /= uvar_overSampling;

     cfg->FreqImax = (INT4)(uvar_FreqBand/cfg->dFreq + 0.5) + 1;  /*Number of frequency values to calculate F for */
   }

 if (LALUserVarWasSet (&uvar_f1dotBand) && (uvar_f1dotBand != 0) )
    cfg->SpinImax=(int)(uvar_f1dotBand/uvar_df1dot+.5)+1;  /*Number of spindown values to calculate F for */
  else
    cfg->SpinImax = 0;

  /*----------------------------------------------------------------------
   * initialize+check  template-grid related parameters 
   */
  if (!uvar_binary)
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
	LALWarning (status, "\nWARNING: We are using skyGridFile, but sky-region was also specified ... will be ignored!\n");
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
	cfg->skyRegion = LALCalloc(1, strlen(uvar_skyRegion)+1);
	strcpy (cfg->skyRegion, uvar_skyRegion);
      }
    else if (haveAlphaDelta)	/* parse this into a sky-region */
      {
	REAL8 eps = 1.0e-9;
	REAL8 a, d, Da, Dd;
	a = uvar_Alpha;
	d = uvar_Delta;
	Da = uvar_AlphaBand + eps;	/* slightly push outwards to make sure boundary-points are included */
	Dd = uvar_DeltaBand + eps;
	/* consistency check either one point given or a 2D region! */
	ASSERT ( (Da && Dd)  || ((Da == 0) && (Dd == 0.0)), status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
      
	cfg->skyRegion = LALMalloc (512); /* should be enough for max 4 points... */
	if ( (Da == 0) || (Dd == 0) ) 	/* only one point */
	  sprintf (cfg->skyRegion, "(%.16f, %.16f)", a, d);
	else				/* or a rectangle */
	  sprintf (cfg->skyRegion, "(%.16f, %.16f), (%.16f, %.16f), (%.16f, %.16f), (%.16f, %.16f)", 
		   a, d, 
		   a + Da, d, 
		   a + Da, d + Dd,
		   a, d + Dd );
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
    INT4 k;

    /* Allocate space for AMCoeffs */
    amc.a = NULL;
    amc.b = NULL;
    TRY (LALSCreateVector(status->statusPtr, &(amc.a), (UINT4) cfg->SFTno), status);
    TRY (LALSCreateVector(status->statusPtr, &(amc.b), (UINT4) cfg->SFTno), status);

    /* Allocate DemodParams structure */
    DemodParams = (NewDemodPar *)LALCalloc(1, sizeof(NewDemodPar));
  
    /* space for sky constants */
    /* Based on maximum index for array of as and bs sky constants as from ComputeSky.c */
    k = 4*(cfg->SFTno-1) + 4; 
    DemodParams->skyConst = (REAL8 *)LALMalloc(k*sizeof(REAL8));
    /* space for spin down params */
    DemodParams->spinDwn = (REAL8 *)LALMalloc(sizeof(REAL8));
  
  }/* end: init AM- and demod-params */

    DETATCHSTATUSPTR (status);
    RETURN (status);

} /* InitFStat() */

/***********************************************************************/
/** Log the all relevant parameters of the present search-run to a log-file.
 * The name of the log-file is "Fstats{uvar_outputLabel}.log".
 * <em>NOTE:</em> Currently this function only logs the user-input and code-versions.
 */
void
WriteFStatLog (LALStatus *stat, char *argv[])
{
    CHAR *logstr = NULL;
    const CHAR *head = "Fstats";
    CHAR command[512] = "";
    UINT4 len;
    CHAR *fname = NULL;
    FILE *fplog;

    INITSTATUS (stat, "WriteFStatLog", rcsid);
    ATTATCHSTATUSPTR (stat);

    /* prepare log-file for writing */
    len = strlen(head) + strlen(".log") +10;
    if (uvar_outputLabel)
      len += strlen(uvar_outputLabel);

    if ( (fname=LALCalloc(len,1)) == NULL) {
      ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }
    strcpy (fname, head);
    if (uvar_outputLabel)
      strcat (fname, uvar_outputLabel);
    strcat (fname, ".log");

    if ( (fplog = fopen(fname, "w" )) == NULL) {
      LALPrintError ("\nFailed to open log-file '%f' for writing.\n\n", fname);
      LALFree (fname);
      ABORT (stat, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
    }

    /* write out a log describing the complete user-input (in cfg-file format) */
    TRY (LALUserVarGetLog (stat->statusPtr, &logstr,  UVAR_LOGFMT_CFGFILE), stat);

    fprintf (fplog, "## LOG-FILE of ComputeFStatistic run\n\n");
    fprintf (fplog, "# User-input:\n");
    fprintf (fplog, "# ----------------------------------------------------------------------\n\n");

    fprintf (fplog, logstr);
    LALFree (logstr);

    /* append an ident-string defining the exact CVS-version of the code used */
    fprintf (fplog, "\n\n# CVS-versions of executable:\n");
    fprintf (fplog, "# ----------------------------------------------------------------------\n");
    fclose (fplog);
    
    sprintf (command, "ident %s | sort -u >> %s", argv[0], fname);
    system (command);	/* we don't check this. If it fails, we assume that */
    			/* one of the system-commands was not available, and */
    			/* therefore the CVS-versions will not be logged */

    LALFree (fname);

    DETATCHSTATUSPTR (stat);
    RETURN(stat);

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
void Freemem(LALStatus *status) 
{
  UINT4 i;

  INITSTATUS (status, "Freemem", rcsid);
  ATTATCHSTATUSPTR (status);

  /* Free SFT data */
  TRY (LALDestroySFTVector (status->statusPtr, &SFTvect), status);	 /* the new way*/
  
  for (i=0; i < GV.SFTno; i++)	/* the old way */
    LALFree (SFTData[i]);
  LALFree (SFTData);


  /* Free timestamps */
  LALFree(timestamps);

  /* Free DemodParams */
  if (DemodParams->amcoe)
    {
      TRY (LALSDestroyVector(status->statusPtr, &(DemodParams->amcoe->a)), status);
      TRY (LALSDestroyVector(status->statusPtr, &(DemodParams->amcoe->b)), status);
    }
  if (DemodParams)
    {
      LALFree(DemodParams->skyConst);
      LALFree(DemodParams->spinDwn);
      LALFree(DemodParams);
    }

     

  if (uvar_binary) {
    LALFree(BinaryBank->BTB);
    LALFree(BinaryBank);
  }

  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (status->statusPtr), status);

  LALFree ( GV.skyRegion );

  /* this comes from clusters.c */
  if (highFLines->clusters) LALFree(highFLines->clusters);
  if (highFLines->Iclust) LALFree(highFLines->Iclust);
  if (highFLines->NclustPoints) LALFree(highFLines->NclustPoints);


  /* Free ephemeris data */
  LALFree(GV.edat->ephemE);
  LALFree(GV.edat->ephemS);
  LALFree(GV.edat);

  DETATCHSTATUSPTR (status);

  /* did we forget anything ? */
  LALCheckMemoryLeaks();

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
EstimateFLines(LALStatus *stat)
{
#ifdef FILE_FTXT  
  FILE *outfile;
#endif
#ifdef FILE_FLINES  
  FILE *outfile1;
#endif
  INT4 i,j,Ntot;   
  INT4 nbins=GV.FreqImax;                /* Number of bins in F */
  REAL8Vector *F1=NULL; 
  REAL8Vector *FloorF1=NULL;             /* Square of SFT */
  REAL4 THR=10.0;
  
  REAL8 dmp;

  OutliersInput  *outliersInput;
  OutliersParams *outliersParams;
  Outliers       *outliers;
  ClustersInput  *clustersInput;
  ClustersParams *SpClParams;
  Clusters       *SpLines=highFLines;
    
  INT2 smallBlock=1;
  INT4 wings;

  INITSTATUS( stat, "EstimateFLines", rcsid );
  ATTATCHSTATUSPTR (stat);

  nbins=(UINT4)nbins;

  THR=uvar_Fthreshold;

/* 0.0002 is the max expected width of the F stat curve for signal */
/* with ~ 10 h observation time */

  dmp=0.5+0.0002/GV.dFreq;
  wings = (INT4) dmp;


#ifdef FILE_FTXT
  /*  file contains freq, PSD, noise floor */
  if(!(outfile=fopen("F.txt","wb"))){
    fprintf(stderr, "Cannot open F.txt file\n");
    return 1;
  }
#endif
  /*  file contanis freq, PSD, noise floor,lines */
#ifdef FILE_FLINES  
  if(!(outfile1=fopen("FLines.txt","wb"))){
    fprintf(stderr, "Cannot open FLines.txt file\n");
    return 1;
  }
#endif


  TRY ( LALDCreateVector(stat->statusPtr, &F1, nbins), stat);
  TRY ( LALDCreateVector(stat->statusPtr, &FloorF1, nbins), stat);
    
  /* loop over SFT data to estimate noise */
  for (j=0;j<nbins;j++){
    F1->data[j]=Fstat.F[j];
    FloorF1->data[j]=1.0;
  }
  
  F1->length=nbins;
  FloorF1->length=nbins;

  /*   j=EstimateFloor(F1, windowSize, FloorF1); */
 
  if ( (outliers = (Outliers *)LALCalloc(1, sizeof(Outliers))) == NULL) {
    ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  outliers->Noutliers=0;

  if ( (outliersParams = (OutliersParams *)LALCalloc(1,sizeof(OutliersParams))) == NULL) {
    ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  if ( (outliersInput = (OutliersInput *)LALCalloc(1,sizeof(OutliersInput))) == NULL) {
    ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  
  outliersParams->Thr=THR/(2.0*medianbias);
  outliersParams->Floor = FloorF1;
  outliersParams->wings=wings; /*these must be the same as ClustersParams->wings */
  outliersInput->ifmin = (INT4) ((uvar_Freq/GV.dFreq)+0.5);
  outliersInput->data = F1;

  /*find values of F above THR and populate outliers with them */
  ComputeOutliers(outliersInput, outliersParams, outliers);


  /*if no outliers were found clean and exit */
   if (outliers->Noutliers == 0){

#ifdef FILE_FTXT
     /*  F.txt file contains freq, F, noise floor of F   */
     for (i=0;i<nbins;i++){ 
       REAL4 freq;
       REAL8 r0,r1;
       freq=uvar_Freq + i*GV.dFreq;
       r0=F1->data[i]*2.0*medianbias;
       r1=FloorF1->data[i]*2.0*medianbias;
       fprintf(outfile,"%f %E %E\n",freq,r0,r1);
     }
#endif     

#ifdef FILE_FTXT
     fclose(outfile);
#endif
#ifdef FILE_FLINES  
     fclose(outfile1);
#endif 
     LALFree(outliers->ratio);
     LALFree(outliers);
     LALFree(outliersParams);
     LALFree(outliersInput);
     TRY ( LALDDestroyVector(stat->statusPtr, &F1), stat);
     TRY ( LALDDestroyVector(stat->statusPtr, &FloorF1), stat);

     /*      fprintf(stderr,"Nclusters zero \n"); */
     /*      fflush(stderr); */

     goto finished;

   } /* if Noutliers == 0 */
  

   /* if outliers are found get ready to identify clusters of outliers*/
   if ( (SpClParams = (ClustersParams*)LALCalloc(1,sizeof(ClustersParams))) == NULL) {
     ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
   }
   
   if ( (clustersInput = (ClustersInput *)LALCalloc(1,sizeof(ClustersInput))) == NULL) {
     ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
   }
      
   SpClParams->wings=wings;
   SpClParams->smallBlock=smallBlock;
   
   clustersInput->outliersInput = outliersInput;
   clustersInput->outliersParams= outliersParams;
   clustersInput->outliers      = outliers;     
   
   /* clusters of outliers in F get written in SpLines which is the global highFLines*/
   TRY (DetectClusters(stat->statusPtr, clustersInput, SpClParams, SpLines), stat);
   
   /*  sum of points in all lines */
   Ntot=0;
   for (i=0;i<SpLines->Nclusters;i++){ 
     Ntot=Ntot+SpLines->NclustPoints[i];
   }
   


#ifdef FILE_FLINES  
   /*  FLines file contains: F, noise floor and lines. */
   for (i=0;i<Ntot;i++){ 
     REAL8 freq;
     REAL8 r0,r1,r2;
     j=SpLines->Iclust[i];
     freq=(uvar_Freq+SpLines->Iclust[i]*GV.dFreq);
     r0=F1->data[j];
     r1=FloorF1->data[j]*2.0*medianbias;
     r2=SpLines->clusters[i]*FloorF1->data[j]*2.0*medianbias;
     fprintf(outfile1,"%20.17f %E %E %E\n",freq,r0,r1,r2);
   }
#endif
#ifdef FILE_FTXT   
   /*  PSD.txt file contains freq, PSD, noise floor   */
   for (i=0;i<nbins;i++){ 
     REAL8 freq;
     REAL8 r0,r1;
     freq=uvar_Freq + i*GV.dFreq;
     r0=F1->data[i]*2.0*medianbias;
     r1=FloorF1->data[i]*2.0*medianbias;
     fprintf(outfile,"%20.17f %E %E\n",freq,r0,r1);
   }
#endif   

#ifdef FILE_FTXT
   fclose(outfile);
#endif
#ifdef FILE_FLINES  
   fclose(outfile1);
#endif   

   TRY ( LALDDestroyVector(stat->statusPtr, &F1), stat);
   TRY ( LALDDestroyVector(stat->statusPtr, &FloorF1), stat);

   LALFree(outliers->ratio);
   LALFree(outliers->outlierIndexes);
   LALFree(outliers);
   LALFree(outliersParams);
   LALFree(outliersInput);
   LALFree(SpClParams);
   LALFree(clustersInput);

 finished:
   DETATCHSTATUSPTR(stat);
   RETURN(stat);

} /* EstimateFLines() */


/** Normalise the SFT-array \em SFTData by the running median.
 * The running median windowSize in this routine determines 
 * the sample bias which, instead of log(2.0), must be 
 * multiplied by F statistics.
 */
void 
NormaliseSFTDataRngMdn(LALStatus *stat)
{
  INT4 i,j,m,lpc,il;                         /* loop indices */
  INT4 Ntot;
  REAL8Vector *Sp=NULL, *RngMdnSp=NULL;   /* |SFT|^2 and its rngmdn  */
  REAL8 B;                          /* SFT Bandwidth */
  REAL8 deltaT,norm,*N, *Sp1;
  INT2 windowSize=uvar_windowsize;                  /* Running Median Window Size*/
  REAL4 xre,xim,xreNorm,ximNorm;

  INITSTATUS( stat, "NormaliseSFTDataRngMdn", rcsid );
  ATTATCHSTATUSPTR (stat);

  if ( !uvar_SignalOnly ) {
    TRY( LALRngMedBias (stat->statusPtr, &medianbias, windowSize), stat);
  }

  TRY ( LALDCreateVector(stat->statusPtr, &Sp, GV.nsamples), stat);
  TRY ( LALDCreateVector(stat->statusPtr, &RngMdnSp, GV.nsamples), stat);

  if( (N = (REAL8 *) LALCalloc(GV.nsamples,sizeof(REAL8))) == NULL) {
    ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  if( (Sp1 = (REAL8 *) LALCalloc(GV.nsamples,sizeof(REAL8))) == NULL) { 
    ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
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
      TRY ( EstimateFloor(stat->statusPtr, Sp, windowSize, RngMdnSp), stat);

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

  TRY ( LALDDestroyVector(stat->statusPtr, &RngMdnSp), stat);
  TRY ( LALDDestroyVector(stat->statusPtr, &Sp), stat);

  DETATCHSTATUSPTR(stat);
  RETURN(stat);
  
} /* NormaliseSFTDataRngMed() */

/**************************************************************************************/

int ReadBinaryTemplateBank(LALStatus *status)
{

  FILE *BTBfp;
  char dmp[256];
  char filename[256];
  INT4 i;
  REAL8 temp1,temp2;

  
  strcpy(filename,uvar_binarytemplatefile);
  /*  something that reads the binary input file to memory */
  if(!(BTBfp=fopen(filename,"r"))){ 
    printf("Cannot open BinaryTemplate file %s\n",filename); 
    return 1;
  } 

  /* allocate memory for header info */
  BinaryBank=(BinaryTemplateBank *)LALMalloc(1*sizeof(BinaryTemplateBank));

  /*printf("opened binary file %s\n",filename);
  fscanf(BTBfp,"%s %le\n",dmp,&(BinaryBank->BTBfmax));
  printf("temp1 read as %le\n",BinaryBank->BTBfmax);
  exit(0);*/
  /* read header information and fill in TemplateBank fields */
  fscanf(BTBfp,"%s%le\n %s%le\n %s%d\n %s%d\n %s%d\n %s%lf\n %s%le\n %s%le\n" \
	 "%s%d\n %s%d\n %s%d\n %s%d\n %s%le\n %s%le\n %s%le\n %s%le\n %s%le\n %s%le\n %s%s\n\n", \
	 dmp,&(BinaryBank->BTBfmax), \
	 dmp,&(BinaryBank->BTBTspan), \
	 dmp,&(BinaryBank->BTBTobsStart).gpsSeconds, \
	 dmp,&(BinaryBank->BTBTobsStart).gpsNanoSeconds, \
	 dmp,&(BinaryBank->BTBNFilters), \
	 dmp,&(BinaryBank->BTBMismatch), \
	 dmp,&(BinaryBank->BTBProjSMaxisMIN), \
	 dmp,&(BinaryBank->BTBProjSMaxisMAX), \
	 dmp,&(BinaryBank->BTBTperiSSBMIN).gpsSeconds, \
	 dmp,&(BinaryBank->BTBTperiSSBMIN).gpsNanoSeconds, \
	 dmp,&(BinaryBank->BTBTperiSSBMAX).gpsSeconds, \
	 dmp,&(BinaryBank->BTBTperiSSBMAX).gpsNanoSeconds, \
	 dmp,&(BinaryBank->BTBEccMIN), \
	 dmp,&(BinaryBank->BTBEccMAX), \
	 dmp,&(BinaryBank->BTBArgPeriMIN), \
	 dmp,&(BinaryBank->BTBArgPeriMAX), \
	 dmp,&(BinaryBank->BTBPeriodMIN), \
	 dmp,&(BinaryBank->BTBPeriodMAX), \
	 dmp,(BinaryBank->BTBversion));

  printf("read header info\n");
  /* Do initial validation of header information */
  if (BinaryBank->BTBfmax<0.0) {
    printf("In BinaryTemplate file %s : header value of fmax < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBTspan<0.0) {
    printf("In BinaryTemplate file %s : header value of Tspan < 0.0\n",filename);
    return 1;
  }
    if ((BinaryBank->BTBTobsStart.gpsSeconds!=GV.Tstart.gpsSeconds)||(BinaryBank->BTBTobsStart.gpsNanoSeconds!=0)) {
      printf("In BinaryTemplate file %s : header value of Tstart != Tstart of data\n",filename);
      return 1;
    }
  if (BinaryBank->BTBNFilters<1) {
    printf("In BinaryTemplate file %s : header value of NFilters < 1\n",filename);
    return 1;
  }
  if (BinaryBank->BTBMismatch<0.0) {
    printf("In BinaryTemplate file %s : header value of Mismatch < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBProjSMaxisMIN<0.0) {
    printf("In BinaryTemplate file %s : header value of Minimum Projected semi-major axis < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBProjSMaxisMAX<0.0) {
    printf("In BinaryTemplate file %s : header value of Maximum Projected semi-major axis < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBProjSMaxisMAX<BinaryBank->BTBProjSMaxisMIN) {
    printf("In BinaryTemplate file %s : header value of Maximum Projected semi-major axis < Minimum\n",filename);
    return 1;
  }
  if (BinaryBank->BTBTperiSSBMIN.gpsSeconds<0) {
    printf("In BinaryTemplate file %s : header value of Minimum SSB time of periapse passage < 0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBTperiSSBMAX.gpsSeconds<0) {
    printf("In BinaryTemplate file %s : header value of Maximum SSB time of periapse passage < 0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBTperiSSBMIN.gpsNanoSeconds<0) {
    printf("In BinaryTemplate file %s : header value of Minimum SSB time of periapse passage (nanoseconds) < 0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBTperiSSBMAX.gpsNanoSeconds<0) {
    printf("In BinaryTemplate file %s : header value of Maximum SSB time of periapse passage (nanoseconds) < 0\n",filename);
    return 1;
  }
  temp1=BinaryBank->BTBTperiSSBMIN.gpsSeconds+1e-9*BinaryBank->BTBTperiSSBMIN.gpsNanoSeconds;
  temp2=BinaryBank->BTBTperiSSBMAX.gpsSeconds+1e-9*BinaryBank->BTBTperiSSBMAX.gpsNanoSeconds;
  if (temp2<temp1) {
    printf("In BinaryTemplate file %s : header value of Maximum Projected semi-major axis < Minimum\n",filename);
    return 1;
  }
  if (BinaryBank->BTBEccMIN<0.0) {
    printf("In BinaryTemplate file %s : header value of Minimum eccentricity < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBEccMAX<0.0) {
    printf("In BinaryTemplate file %s : header value of Maximum eccentricity < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBEccMAX<BinaryBank->BTBEccMIN) {
    printf("In BinaryTemplate file %s : header value of Maximum eccentricity < Minimum\n",filename);
    return 1;
  }
  if ((BinaryBank->BTBArgPeriMIN<0.0)||(BinaryBank->BTBArgPeriMIN>LAL_TWOPI)) {
    printf("In BinaryTemplate file %s : header value of Minimum argument of periapse not in range (0 - 2*PI)\n",filename);
    return 1;
  }
  if ((BinaryBank->BTBArgPeriMAX<0.0)||(BinaryBank->BTBArgPeriMAX>LAL_TWOPI)) {
    printf("In BinaryTemplate file %s : header value of Maximum argument of periapse not in range (0 - 2*PI)\n",filename);
    return 1;
  }
  if (BinaryBank->BTBArgPeriMAX<BinaryBank->BTBArgPeriMIN) {
    printf("In BinaryTemplate file %s : header value of Maximum argument of periapse < Minimum\n",filename);
    return 1;
  }
  if ((BinaryBank->BTBArgPeriMIN<0.0)||(BinaryBank->BTBArgPeriMIN>LAL_TWOPI)) {
    printf("In BinaryTemplate file %s : header value of Minimum argument of periapse not in range (0 - 2*PI)\n",filename);
    return 1;
  }
  if ((BinaryBank->BTBArgPeriMAX<0.0)||(BinaryBank->BTBArgPeriMAX>LAL_TWOPI)) {
    printf("In BinaryTemplate file %s : header value of Maximum argument of periapse not in range (0 - 2*PI)\n",filename);
    return 1;
  }
  if (BinaryBank->BTBArgPeriMAX<BinaryBank->BTBArgPeriMIN) {
    printf("In BinaryTemplate file %s : header value of Maximum argument of periapse < Minimum\n",filename);
    return 1;
  }  
  if (BinaryBank->BTBPeriodMIN<0.0) {
    printf("In BinaryTemplate file %s : header value of Minimum period < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBPeriodMAX<0.0) {
    printf("In BinaryTemplate file %s : header value of Maximum period < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BTBPeriodMAX<BinaryBank->BTBPeriodMIN) {
    printf("In BinaryTemplate file %s : header value of Maximum period < Minimum\n",filename);
    return 1;
  }  
  
  /* allocate memory for templates */
  (BinaryBank->BTB)=(BinaryTemplate *)LALMalloc(BinaryBank->BTBNFilters*sizeof(BinaryTemplate));
  
  /* Now read in all templates into memory */
  i=0;
  while (i<BinaryBank->BTBNFilters) {
    fscanf(BTBfp,"%le%le%d%d%le%le\n",
		 &(BinaryBank->BTB[i]).ProjSMaxis,
		 &(BinaryBank->BTB[i]).Period,
		 &(BinaryBank->BTB[i]).TperiSSB.gpsSeconds,
		 &(BinaryBank->BTB[i]).TperiSSB.gpsNanoSeconds,
		 &(BinaryBank->BTB[i]).Eccentricity,
		 &(BinaryBank->BTB[i]).ArgPeri); 
  
    /* printf("ProjSMaxis = %le Period = %le TperiSSB.sec = %d TperiSSB.nano = %d Eccentricity = %le ArgPeri = %le\n", \
	   BinaryBank->BTB[i].ProjSMaxis,BinaryBank->BTB[i].Period,BinaryBank->BTB[i].TperiSSB.gpsSeconds, \
	   BinaryBank->BTB[i].TperiSSB.gpsNanoSeconds,BinaryBank->BTB[i].Eccentricity,BinaryBank->BTB[i].ArgPeri); */
    
    /* check each template before being used in the search */
    /* may also add a check that the templates are within the specified boundaries */
    if ((BinaryBank->BTB[i]).ProjSMaxis<0.0) {
      printf("In BinaryTemplate file %s : template #[%d] has negative projected semi-major axis\n",filename,i);
      return 1;
    }
    if ((BinaryBank->BTB[i]).Period<0.0) {
      printf("In BinaryTemplate file %s : template #[%d] has negative orbital period\n",filename,i);
      return 1;
    }  
    if ((BinaryBank->BTB[i]).TperiSSB.gpsSeconds<0) {
      printf("In BinaryTemplate file %s : template #[%d] has negative time of periapse passage (sec)\n",filename,i);
      return 1;
    }  
    if ((BinaryBank->BTB[i]).TperiSSB.gpsNanoSeconds<0) {
      printf("In BinaryTemplate file %s : template #[%d] has negative time of periapse passage (nanosec)\n",filename,i);
      return 1;
    }  
    if ((BinaryBank->BTB[i]).Eccentricity<0.0) {
      printf("In BinaryTemplate file %s : template #[%d] has negative eccentricity\n",filename,i);
      return 1;
    }
    if (((BinaryBank->BTB[i]).ArgPeri<0.0)||((BinaryBank->BTB[i]).ArgPeri>LAL_TWOPI)) {
      printf("In BinaryTemplate file %s : template #[%d] has argument of periapse out of range (0 - 2*PI)\n",filename,i);
      return 1;
    }  
    i++;
  }

  

  return 1;
  
}

/**************************************************************************************/

/** Local version of LALDemod() for local testing and general fooling around...
 */
#define SMALL	0.000000001
void 
NewLALDemod(LALStatus *stat, FStatisticBand **FBand, FFT **input, NewDemodPar *params) 
{ 

  INT4 alpha,i;                 /* loop indices */
  REAL8	*xSum=NULL, *ySum=NULL;	/* temp variables for computation of fs*as and fs*bs */
  INT4 s;		        /* local variable for spinDwn calcs. */
  REAL8	xTemp;	                /* temp variable for phase model */
  REAL8	deltaF;	                /* width of SFT band */
  INT4	k, k1;	                /* defining the sum over which is calculated */
  REAL8 *skyConst;	        /* vector of sky constants data */
  REAL8 *spinDwn;	        /* vector of spinDwn parameters (maybe a structure? */
  INT4	spOrder;	        /* maximum spinDwn order */
  REAL8	x;		        /* local variable for holding x */
  REAL8	realXP, imagXP; 	/* temp variables used in computation of */
  REAL8	realP, imagP;	        /* real and imaginary parts of P, see CVS */
  INT4	nDeltaF;	        /* number of frequency bins per SFT band */
  INT4	sftIndex;	        /* more temp variables */
  REAL8	y;		        /* local variable for holding y */
  REAL8 realQ, imagQ;
  REAL8 *sinVal,*cosVal;        /*LUT values computed by the routine do_trig_lut*/
  INT4  res;                    /*resolution of the argument of the trig functions: 2pi/res.*/
  INT4 *tempInt1;
  INT4 index;
  REAL8 FaSq;
  REAL8 FbSq;
  REAL8 FaFb;
  COMPLEX16 Fa, Fb;
  REAL8 f;
  REAL8 A=params->amcoe->A,B=params->amcoe->B,C=params->amcoe->C,D=params->amcoe->D;
  INT4 M=params->SFTno;

  FStatisticBand *ret = NULL;	/* return-struct */
  UINT4 nBins;			/* number of frequency-bins */


  INITSTATUS( stat, "LALDemod", rcsid);
  ATTATCHSTATUSPTR (stat);

  ASSERT ( FBand, stat, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *FBand == NULL, stat, COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);

  nBins = params->imax;	/* number of frequency-bins to calculate Fstatistic for */

  /* allocate memory for FBand return-struct.
   * FIXME: for now we're a bit sloppy with leaks in  case of errors 
  */
  if ( (ret = (FStatisticBand*)LALCalloc(1, sizeof(FStatisticBand))) == NULL) {
    ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  TRY ( LALDCreateVector (stat->statusPtr, &(ret->freqs), nBins), stat);
  TRY ( LALDCreateVector (stat->statusPtr, &(ret->F), nBins), stat);
  TRY ( LALZCreateVector (stat->statusPtr, &(ret->Fa), nBins), stat);
  TRY ( LALZCreateVector (stat->statusPtr, &(ret->Fb), nBins), stat);

  /* variable redefinitions for code readability */
  spOrder=params->spinDwnOrder;
  spinDwn=params->spinDwn;
  skyConst=params->skyConst;
  deltaF=(*input)->fft->deltaF;
  nDeltaF=(*input)->fft->data->length;

  /* res=10*(params->mCohSFT); */
  /* This size LUT gives errors ~ 10^-7 with a three-term Taylor series */
   res=64;
   sinVal=(REAL8 *)LALMalloc((res+1)*sizeof(REAL8));
   cosVal=(REAL8 *)LALMalloc((res+1)*sizeof(REAL8)); 
   for (k=0; k<=res; k++){
     sinVal[k]=sin((LAL_TWOPI*k)/res);
     cosVal[k]=cos((LAL_TWOPI*k)/res);
   }

  /* this loop computes the values of the phase model */
  xSum=(REAL8 *)LALMalloc(params->SFTno*sizeof(REAL8));
  ySum=(REAL8 *)LALMalloc(params->SFTno*sizeof(REAL8));
  tempInt1=(INT4 *)LALMalloc(params->SFTno*sizeof(INT4));
  for(alpha=0;alpha<params->SFTno;alpha++){
    tempInt1[alpha]=2*alpha*(spOrder+1)+1;
    xSum[alpha]=0.0;
    ySum[alpha]=0.0;
    for(s=0; s<spOrder;s++) {
      xSum[alpha] += spinDwn[s] * skyConst[tempInt1[alpha]+2+2*s]; 	
      ySum[alpha] += spinDwn[s] * skyConst[tempInt1[alpha]+1+2*s];
    }
  }


  /* Loop over frequencies to be demodulated */
  for(i=0 ; i< params->imax  ; i++ )
  {
    Fa.re =0.0;
    Fa.im =0.0;
    Fb.re =0.0;
    Fb.im =0.0;

    f=params->f0+i*params->df;

    /* Loop over SFTs that contribute to F-stat for a given frequency */
    for(alpha=0;alpha<params->SFTno;alpha++)
      {
	REAL8 tsin, tcos, tempFreq;
	COMPLEX8 *Xalpha=input[alpha]->fft->data->data;
	REAL4 a = params->amcoe->a->data[alpha];
	REAL4 b = params->amcoe->b->data[alpha];

	xTemp=f*skyConst[tempInt1[alpha]]+xSum[alpha];

	realXP=0.0;
	imagXP=0.0;
	      
	/* find correct index into LUT -- pick closest point */
	tempFreq=xTemp-(INT4)xTemp;
	index=(INT4)(tempFreq*res+0.5);
	      
	{
	  REAL8 d=LAL_TWOPI*(tempFreq-(REAL8)index/(REAL8)res);
	  REAL8 d2=0.5*d*d;
	  REAL8 ts=sinVal[index];
	  REAL8 tc=cosVal[index];
		
	  tsin=ts+d*tc-d2*ts;
	  tcos=tc-d*ts-d2*tc-1.0;
	}
		     
	tempFreq=LAL_TWOPI*(tempFreq+params->Dterms-1);
	k1=(INT4)xTemp-params->Dterms+1;
	/* Loop over terms in dirichlet Kernel */
	for(k=0;k<2*params->Dterms;k++)
	  {
	    COMPLEX8 Xalpha_k;
	    x=tempFreq-LAL_TWOPI*(REAL8)k;
	    realP=tsin/x;
	    imagP=tcos/x;

	    /* If x is small we need correct x->0 limit of Dirichlet kernel */
	    if(fabs(x) < SMALL) 
	      {
		realP=1.0;
		imagP=0.0;
	      }	 
 
	    sftIndex=k1+k-params->ifmin;

	    /* these four lines compute P*xtilde */
	    Xalpha_k=Xalpha[sftIndex];
	    realXP += Xalpha_k.re*realP;
	    realXP -= Xalpha_k.im*imagP;
	    imagXP += Xalpha_k.re*imagP;
	    imagXP += Xalpha_k.im*realP;
	  }
      
	y=-LAL_TWOPI*(f*skyConst[tempInt1[alpha]-1]+ySum[alpha]);

	realQ = cos(y);
	imagQ = sin(y);

	/* implementation of amplitude demodulation */
	{
	  REAL8 realQXP = realXP*realQ-imagXP*imagQ;
	  REAL8 imagQXP = realXP*imagQ+imagXP*realQ;
	  Fa.re += a*realQXP;
	  Fa.im += a*imagQXP;
	  Fb.re += b*realQXP;
	  Fb.im += b*imagQXP;
	}
    }      

    FaSq = Fa.re*Fa.re+Fa.im*Fa.im;
    FbSq = Fb.re*Fb.re+Fb.im*Fb.im;
    FaFb = Fa.re*Fb.re+Fa.im*Fb.im;
	  	  	
    ret->F->data[i] = (4.0/(M*D))*(B*FaSq + A*FbSq - 2.0*C*FaFb);
    ret->Fa->data[i] = Fa;
    ret->Fb->data[i] = Fb;
    ret->freqs->data[i] = f;

  } /* for i < imax */

  /* Clean up */
  LALFree(tempInt1);
  LALFree(xSum);
  LALFree(ySum);
  
  LALFree(sinVal);
  LALFree(cosVal);

  /* return result */
  *FBand = ret;

  DETATCHSTATUSPTR (stat);
  RETURN( stat );

} /* NewLALDemod() */


void
writeCOMPLEX16Vector(LALStatus *stat, COMPLEX16Vector *vect, const CHAR *fname)
{
  FILE *fp;
  UINT4 i;

  INITSTATUS( stat, "LALDemod", rcsid);
  ATTATCHSTATUSPTR (stat);

  ASSERT (vect, stat, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (fname, stat, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  if ( (fp = fopen(fname, "wb")) == NULL) 
    {
      LALPrintError ("Failed to open file '%f' for writing.\n", fname);
      ABORT (stat, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
    }

  for (i=0; i < vect->length; i++)
    fprintf (fp, "%d %g %g \n", i, vect->data[i].re, vect->data[i].im );

  fclose (fp);

  DETATCHSTATUSPTR (stat);
  RETURN( stat );

} /* writeCOMPLEX16Vector() */
