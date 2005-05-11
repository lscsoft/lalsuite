/*********************************************************************************/
/*                    F-statistic generation code for known pulsars              */
/*                                                                               */
/*			      Y. Ioth, M.A. Papa, X. Siemens, R.Prix             */
/*                                                                               */
/*                 Albert Einstein Institute/UWM - started September 2002        */
/*********************************************************************************/
/*                                                                               */
/* Version of ComputeFStatistic that has an input option to deal with Binary     */
/* systems.  With this option switched off the code runs as designed for         */
/* isolated pulsar searches.                                                     */
/*                                                                               */
/* Binary modifications added by Chris Messenger (University of Birmingham UK)   */
/*********************************************************************************/
#include <lal/UserInput.h>
#include <lal/LALDemod.h>
#include <lal/RngMedBias.h>
#include <lalapps.h>


#include "ComputeFStatisticBinary_v2.h"    /* BINARY-MOD - Using modified header file */
#include "GenerateBinaryMesh_v1.h"
#include "ReadSourceFile_v1.h"
#include "clusters.h"
#include "../FDS_isolated/DopplerScan.h"

RCSID( "$Id$");

/* BOINC should be set to 1 to be run under BOINC */
#ifndef USE_BOINC
#define USE_BOINC 0
#endif
#if USE_BOINC
#define USE_BOINC_DEBUG 0
/* for getpid() */
#include <sys/types.h>
#include <unistd.h>

typedef int bool;
extern int boinc_init(bool standalone);
extern int boinc_finish(int);
extern int boinc_resolve_filename(const char*, char*, int len);
extern int boinc_init_graphics();
extern int boinc_finish_graphics();
void use_boinc_filename1(char** orig_name);
void use_boinc_filename0(char* orig_name);
#endif /* USE_BOINC */


/*
#define NEARESTGRIDPOINTS_ON


#define DEBG_FAFB                
#define DEBG_ESTSIGPAR
#define DEBG_SGV 
*/

/* If FILE_FILENAME is defined, then print the corresponding file  */
#define FILE_FSTATS /* outputs the basic statistics of each Fstat cluster found for each template */ 
#define FILE_FMAX  /* outputs the loudest Fstat for each template */
/* #define FILE_FLINES */
/* #define FILE_FTXT */
/* #define FILE_PSD */    /* outputs frequency + <UNnormalised PSD> + <UNnormalised PSD RngMed> */
/* #define FILE_PSDLINES */  /* outputs nothing at present */
/* #define FILE_SPRNG */  /* outputs frequency + <RngMed> */


/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
#define COMPUTEFSTATC_ENULL 		1
#define COMPUTEFSTATC_ESYS     		2
#define COMPUTEFSTATC_EINPUT   		3

#define COMPUTEFSTATC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATC_MSGESYS		"System call failed (probably file IO)"
#define COMPUTEFSTATC_MSGEINPUT   	"Invalid input"
/*************************************************** </lalErrTable> */


/*----------------------------------------------------------------------
 * User-variables: provided either by default, config-file or command-line */
INT4 uvar_dterms;
CHAR* uvar_IFO;
BOOLEAN uvar_SignalOnly;
BOOLEAN uvar_EstimSigParam;
BOOLEAN uvar_binary;  /* BINARY-MOD - added a flag to indicate binary search */
REAL8 uvar_dopplermax;  /* BINARY-MOD - added this because different binary systems will require different wings */
CHAR *uvar_binarytemplatefile;  /* BINARY-MOD - added for binary template file */
INT4 uvar_windowsize;          /* BINARY-MOD - added because of shorter SFT's */
REAL8 uvar_Freq;
REAL8 uvar_dFreq;
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
CHAR *uvar_sourcefile;
CHAR *uvar_source;
REAL8 uvar_overres;

/* try this */
BOOLEAN uvar_openDX;
/*----------------------------------------------------------------------*/

FFT **SFTData=NULL;                 /* SFT Data for LALDemod */
DemodPar *DemodParams  = NULL;      /* Demodulation parameters for LALDemod */
LIGOTimeGPS *timestamps=NULL;       /* Time stamps from SFT data */
LALFstat Fstat;
AMCoeffs amc;
REAL8 MeanOneOverSh=0.0;
REAL8 Alpha,Delta;
BinaryTemplateBank *BinaryBank=NULL;    /* BINARY-MOD - pointer to structure to store binary template bank */
BinaryTemplate thisBinaryTemplate; /* BINARY-MOD - structure to store current binary template */
REAL8 bin_SMaxis;                 /* BINARY-MOD - Orbital semi-major axis of source system */
REAL8 bin_Period;                 /* BINARY-MOD - Orbital period of source system */
REAL8 bin_Eccentricity;           /* BINARY-MOD - Orbital eccentricity of source system */
REAL8 bin_ArgPeri;                /* BINARY-MOD - Argument of periapse of source system */
LIGOTimeGPS bin_TperiSSB;         /* BINARY-MOD - Time of periapse passage as defined in the SSB */
Clusters HFLines, HPLines;
Clusters *highSpLines=&HPLines, *highFLines=&HFLines;
/* #ifdef FILE_FMAX     */
FILE *fpmax;
/* #endif */
/* #ifdef FILE_STATS */
FILE *fpstat;
/* #endif     */
REAL8 medianbias=1.0;

/* we use this to avoid closing then re-opening the same file */
FILE *fp_mergedSFT=NULL;

DopplerScanState thisScan;
ConfigVariables GV;

/* local prototypes */
void CreateDemodParams (LALStatus *status);
void CreateBinaryDemodParams (LALStatus *status);    /* BINARY-MOD - declaring a slightly modified version of CreateDemodParams */
void AllocateMem (LALStatus *status);
void SetGlobalVariables (LALStatus *status, ConfigVariables *cfg);
void CreateNautilusDetector (LALStatus *status, LALDetector *Detector);
void Freemem (LALStatus *status);

INT4 EstimatePSDLines(LALStatus *status);
INT4 EstimateFLines(LALStatus *status);
INT4 NormaliseSFTDataRngMdn(LALStatus *status);

INT4 NormaliseSFTData(void);
INT4 ReadSFTData (void);
INT4 EstimateSignalParameters(INT4 * maxIndex);
INT4 writeFLines(INT4 *maxIndex);
INT4 PrintTopValues(REAL8 TwoFthr, INT4 ReturnMaxN);
INT4 EstimateFloor(REAL8Vector *Sp, INT2 windowSize, REAL8Vector *SpFloor);
int compare(const void *ip, const void *jp);
INT4 writeFaFb(INT4 *maxIndex);
INT4 NormaliseSFTData(void);

void initUserVars (LALStatus *Stat);

/* BINARY-MOD - declaration of binary template reading function */
INT4 ReadBinaryTemplateBank(void);

/*----------------------------------------------------------------------
 * Helper function (Yousuke): 
 * Refine the skyRegion to search only at neighboring grid points of the 
 * center of the original skyRegion. 
 *----------------------------------------------------------------------*/
void InitDopplerScanOnRefinedGrid ( LALStatus *status, DopplerScanState *theScan, DopplerScanInit *scanInit);

#define EPHEM_YEARS  "00-04"
#define SFT_BNAME  ""

#define TRUE (1==1)
#define FALSE (1==0)

extern int vrbflg;


/*----------------------------------------------------------------------
 * MAIN
 *----------------------------------------------------------------------*/
int main(int argc,char *argv[]) 
{
  INT4 *maxIndex=NULL; /*  array that contains indexes of maximum of each cluster */
  DopplerPosition dopplerpos;
  SkyPosition thisPoint; 
  CHAR Fstatsfilename[256];         /* Fstats file name*/
  CHAR Fmaxfilename[256];           /* Fmax file name*/
  INT4 s;
  DopplerScanInit scanInit;
  LIGOTimeGPS t0, t1;
  REAL8 duration;
  FILE *fpOut=NULL;
  UINT4 counter;
  LALStatus status = blank_status;	/* initialize status */
  binarysource sourceparams;
  LIGOTimeGPS *dummyGPS=NULL;

  lalDebugLevel = 0;  
  vrbflg = 1;	/* verbose error-messages */
  
#if USE_BOINC
  /* boinc_init() needs to be run before any boinc_api functions are used */
  boinc_init(FALSE);
  boinc_init_graphics();
#if USE_BOINC_DEBUG
  {
    char commandstring[256];
    /* char *cmd_name = argv[0]; */
    pid_t process_id=getpid();
    sprintf(commandstring,"ddd %s %d &","../../projects/ein*/einstein*" ,process_id);
    system(commandstring);
    sleep(20);
  }
#endif /*USE_BOINC_DEBUG*/
#endif /*USE_BOINC*/
  
  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  LAL_CALL (LALGetDebugLevel (&status, argc, argv, 'v'), &status);
  LAL_CALL (initUserVars (&status), &status); 	

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);	
  
  if (uvar_help)
    exit (0);  
  
  if ( !uvar_binary ) {
    LALPrintError ("\nSorry, this code is not functional in the NON-binary case\n\n");
    exit (COMPUTEFSTATC_EINPUT);
  }

  LAL_CALL ( SetGlobalVariables (&status, &GV), &status);
  
  LAL_CALL ( AllocateMem(&status), &status);
  

  if (ReadSFTData()) return 4;
  
  /*  This fills-in highSpLines that are then used by NormaliseSFTRngMdn */
#if 0
  if (uvar_SignalOnly!=1){
    if (EstimatePSDLines(&status)) return 6;
  }
#endif

  LAL_CALL (NormaliseSFTDataRngMdn(&status), &status);


#ifdef FILE_FMAX  
  /*   open file */
  strcpy(Fmaxfilename,"Fmax");
  if (uvar_outputLabel)
    strcat(Fmaxfilename,uvar_outputLabel);
#if USE_BOINC
  use_boinc_filename0(Fmaxfilename);
#endif /* USE_BOINC */
  if (!(fpmax=fopen(Fmaxfilename,"w"))){
    fprintf(stderr,"in Main: unable to open Fmax file %s\n", Fmaxfilename);
    return 2;
  }
#endif
#ifdef FILE_FSTATS  
  /*      open file */
  strcpy(Fstatsfilename,"Fstats");
  if ( LALUserVarWasSet(&uvar_outputLabel) )
    strcat(Fstatsfilename,uvar_outputLabel);
#if USE_BOINC
  use_boinc_filename0(Fstatsfilename);
#endif /* USE_BOINC */
  if (!(fpstat=fopen(Fstatsfilename,"w"))){
    fprintf(stderr,"in Main: unable to open Fstats file\n");
    return 2;
  }
#endif



  if (!uvar_binary) {

    /* prepare initialization of DopplerScanner to step through paramter space */
    scanInit.dAlpha = uvar_dAlpha;
    scanInit.dDelta = uvar_dDelta;
    scanInit.gridType = uvar_gridType;
    scanInit.metricType = uvar_metricType;
    scanInit.metricMismatch = uvar_metricMismatch;
    t0 = SFTData[0]->fft->epoch;
    t1 = SFTData[GV.SFTno-1]->fft->epoch;
    scanInit.obsBegin = t0;
    LAL_CALL ( LALDeltaFloatGPS ( &status, &duration, &t1, &t0), &status);
    scanInit.obsDuration = duration + GV.tsft;
    scanInit.fmax  = uvar_Freq;
    if (uvar_FreqBand > 0) scanInit.fmax += uvar_FreqBand;
    scanInit.Detector = &GV.Detector;
    scanInit.ephemeris = GV.edat;		/* used by Ephemeris-based metric */
    scanInit.skyGridFile = uvar_skyGridFile;
  
  


    if (lalDebugLevel) LALPrintError ("\nSetting up template grid ...");

    LAL_CALL ( InitDopplerScan ( &status, &thisScan, &scanInit), &status); 

    /*----------------------------------------------------------------------*/
    if (lalDebugLevel) LALPrintError ("done.\n");
    if ( uvar_outputSkyGrid ) {
      LALPrintError ("\nNow writing sky-grid into file '%s' ...", uvar_outputSkyGrid);
      LAL_CALL (writeSkyGridFile ( &status, thisScan.grid, uvar_outputSkyGrid, &scanInit), &status);
      LALPrintError (" done.\n\n");
    }

  
  
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
    LAL_CALL (ReadBinaryTemplateBank(), &status);  
    
    /* read in source file paramters (only care about RA, dec) */
    if (uvar_sourcefile!=NULL) {
      
      if (ReadSource(uvar_sourcefile,uvar_source,dummyGPS,&sourceparams)) return 1;;
      Alpha=sourceparams.skypos.ra;
      Delta=sourceparams.skypos.dec;
      
    }
    else{
      Alpha=uvar_Alpha;   /* BINARY-MOD - Also set sky location */
      Delta=uvar_Delta;
    }
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

  counter = 0;
  while (1)
    {
      /* BINARY-MOD - option to select next skyposition */
      if (!uvar_binary) 
	{

	  LAL_CALL (NextDopplerPos( &status, &dopplerpos, &thisScan ), &status);
	  /* Have we scanned all DopplerPositions yet? */
	  if (thisScan.state == STATE_FINISHED)
	    break;
	  
	  Alpha = thisPoint.longitude;
	  Delta = thisPoint.latitude;
	  
	  LAL_CALL (CreateDemodParams(&status), &status);
	  
	} /* if (!uvar_binary) */

      /* BINARY-MOD - option to select next Binary template */
      if (uvar_binary) 
	{
	  
	  /* Have we scanned all DopplerPositions yet? */
	  if (((UINT4)counter) >= BinaryBank->BMFheader.Nfilters)
	    break;


	  thisBinaryTemplate=(BinaryBank->BTB[counter]); 
	  /*thisBinaryTemplate.ProjSMaxis=2.0;*/
	  /* thisBinaryTemplate.ProjSMaxis=(BinaryBank->BTB[counter]).ProjSMaxis; */
	  /*LAL_CALL (NextBinaryTemplate(&status, &BinaryBank ,&thisBinaryTemplate, &BTNumber ) , &status); */
	  /*printf("filled current template\n");
	  printf("the bank has smaxis as %le\n",BinaryBank->BTB[counter].ProjSMaxis);
	  printf("current smaxis is %le\n",thisBinaryTemplate.ProjSMaxis);
	  exit(0); */
	  /*printf("Current template is\n[%le,%le,%d,%d,%le,%le]\n",thisBinaryTemplate.ProjSMaxis, \
		 thisBinaryTemplate.Period,thisBinaryTemplate.TperiSSB.gpsSeconds, \
		 thisBinaryTemplate.TperiSSB.gpsNanoSeconds,thisBinaryTemplate.Eccentricity,  \
		 thisBinaryTemplate.ArgPeri);*/

	  /* printf("counter is %d Nfilters is %d\n",counter,BinaryBank->BTBNFilters); */
	  
	  
	  LAL_CALL (CreateBinaryDemodParams(&status), &status);
	  

	}  /* if (uvar_binary) */

      /* loop over spin params */
      for(s=0;s<=GV.SpinImax;s++)
	{
	  /* BINARY-MOD - only deal with spin down for isolated at this stage */
	  if (!uvar_binary) DemodParams->spinDwn[0]=uvar_f1dot + s*uvar_df1dot;
	  LAL_CALL (LALDemod (&status, &Fstat, SFTData, DemodParams), &status);

	  /*  This fills-in highFLines that are then used by discardFLines */
	  if (GV.FreqImax > 5) {
	  
	    LAL_CALL (EstimateFLines(&status), &status);
	  }
	  
	  /* now, if user requested it, we output ALL F-statistic results */
	  if ((fpOut)&&(!uvar_binary)) 
	    {
	      INT4 i;
	    
	      for(i=0;i < GV.FreqImax ;i++)
		{
		  fprintf (fpOut, "%20.17f %20.17f %20.17f %20.17f\n", 
			   uvar_Freq + i*GV.dFreq, Alpha, Delta, 2.0*medianbias*Fstat.F[i]);
		}

	    } /* if outputFstat and not binary */

	  /* BINARY-MOD - output binary search F-stat with binary params */
	  if ((fpOut)&&(uvar_binary)) 
	    {
	      INT4 i;
	    
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
	  if (highFLines != NULL && highFLines->Nclusters > 0){
	    
	    maxIndex=(INT4 *)LALMalloc(highFLines->Nclusters*sizeof(INT4));
	    
	    /*  for every cluster writes the information about it in file Fstats */

	    if (writeFLines(maxIndex)){
	      fprintf(stderr, "%s: trouble making file Fstats\n", argv[0]);
	      return 6;
	    }
	  }
	  
	  if( uvar_EstimSigParam &&(highFLines !=NULL) && (highFLines->Nclusters >0))
	    if(writeFaFb(maxIndex)) return 255;
	  
	  
	  if( uvar_EstimSigParam &&(highFLines !=NULL) &&(highFLines->Nclusters >0)) {
	    if (EstimateSignalParameters(maxIndex)) return 7;
	  }
	 
	    if (PrintTopValues(/* thresh */ 0.0, /* max returned */ 1))
	    LALPrintError ("%s: trouble making files Fmax and/or Fstats\n", argv[0]);
	 

	  if (highFLines != NULL && highFLines->Nclusters > 0){
	    LALFree(maxIndex);
	  }

	  /* Set the number of the clusters detected to 0 at each iteration 
	     of the sky-direction and the spin down */
	  highFLines->Nclusters=0;

	} /* For GV.spinImax */

      counter ++;
      if (lalDebugLevel) LALPrintError ("Search progress: %5.1f%%", 
					(100.0* counter / thisScan.numGridPoints));
      
    } /*  while SkyPos */

  if (uvar_outputFstat && fpOut)
    fclose (fpOut);

  if (lalDebugLevel) LALPrintError ("\nSearch finished.\n");

#ifdef FILE_FMAX  
  fclose(fpmax);
#endif
#ifdef FILE_FSTATS  
  fclose(fpstat);
#endif
  LAL_CALL (Freemem(&status), &status);

#if USE_BOINC
  boinc_finish_graphics();
  boinc_finish(0);
#endif

  return 0;

} /* main() */


/* register all our "user-variables", which can be read from cmd-line and config-file */
void
initUserVars (LALStatus *Stat)
{
  INITSTATUS( Stat, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (Stat);

  /* set a few defaults */
  uvar_dterms 	= 16;
  uvar_FreqBand = 0.0;
  uvar_dFreq = 0;

  uvar_binary=FALSE;       /* BINARY-MOD - Set binary flag to FALSE for safety */
  uvar_sourcefile = LALMalloc(512);
  uvar_source = LALMalloc(512);
  uvar_sourcefile=NULL;
  uvar_source=NULL;
  uvar_overres=0.0;

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

  uvar_BaseName	= LALCalloc (1, strlen(SFT_BNAME)+1);
  strcpy (uvar_BaseName, SFT_BNAME);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar_EphemDir = LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar_EphemDir, DEFAULT_EPHEMDIR);

  uvar_SignalOnly = FALSE;
  uvar_EstimSigParam = FALSE;
 
  uvar_f1dot = 0.0;
  uvar_df1dot 	= 0.0;
  uvar_f1dotBand = 0.0;
  
  uvar_Fthreshold = 10.0;
  uvar_metricType =  LAL_PMETRIC_COH_PTOLE_ANALYTIC;
  uvar_gridType = GRID_FLAT;

  uvar_metricMismatch = 0.02;

  uvar_help = FALSE;
  uvar_outputLabel = NULL;

  uvar_outputFstat = NULL;
  uvar_openDX = FALSE;	/* write openDX-compatible output-file */

  uvar_skyGridFile = NULL;

  uvar_workingDir = LALMalloc(512);
  strcpy(uvar_workingDir, ".");


  /* register all our user-variables */
 
  /* BINARY-MOD - added binary flag, dopplermax variable, running median window size variable */ 

  LALregINTUserVar(Stat,	dterms,		't', UVAR_OPTIONAL, "Number of terms to keep in Dirichlet kernel sum");
  LALregREALUserVar(Stat, 	Freq, 		'f', UVAR_REQUIRED, "Starting search frequency in Hz");
  LALregREALUserVar(Stat, 	FreqBand, 	'b', UVAR_OPTIONAL, "Search frequency band in Hz");
  LALregREALUserVar(Stat, 	dFreq, 		'r', UVAR_OPTIONAL, "Frequency resolution in Hz (default: 1/(2*Tsft*Nsft)");
  LALregREALUserVar(Stat, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserVar(Stat, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
  LALregREALUserVar(Stat, 	AlphaBand, 	'z', UVAR_OPTIONAL, "Band in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(Stat, 	DeltaBand, 	'c', UVAR_OPTIONAL, "Band in delta (equatorial coordinates) in radians");
  LALregREALUserVar(Stat, 	dAlpha, 	'l', UVAR_OPTIONAL, "Resolution in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(Stat, 	dDelta, 	'g', UVAR_OPTIONAL, "Resolution in delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(Stat,	DataDir, 	'D', UVAR_OPTIONAL, "Directory where SFT's are located");
  LALregSTRINGUserVar(Stat,	mergedSFTFile, 	'B', UVAR_OPTIONAL, "Merged SFT's file to be used"); 
  LALregSTRINGUserVar(Stat,	BaseName, 	'i', UVAR_OPTIONAL, "The base name of the input  file you want to read");
  LALregSTRINGUserVar(Stat,	EphemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(Stat,	EphemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregSTRINGUserVar(Stat, 	IFO, 		'I', UVAR_REQUIRED, "Detector: GEO(0), LLO(1), LHO(2), NAUTILUS(3), VIRGO(4), TAMA(5), CIT(6)");
  LALregBOOLUserVar(Stat, 	SignalOnly, 	'S', UVAR_OPTIONAL, "Signal only flag");
  LALregBOOLUserVar(Stat, 	binary, 	'u', UVAR_OPTIONAL, "Binary search flag");
  LALregREALUserVar(Stat, 	dopplermax, 	'q', UVAR_OPTIONAL, "Maximum doppler shift expected");  
  LALregREALUserVar(Stat, 	f1dot, 		's', UVAR_OPTIONAL, "First spindown parameter f1dot");
  LALregREALUserVar(Stat, 	f1dotBand, 	'm', UVAR_OPTIONAL, "Search-band for f1dot");
  LALregREALUserVar(Stat, 	df1dot, 	'e', UVAR_OPTIONAL, "Resolution for f1dot (default 1/(2*Tobs*Tsft*Nsft)");
  LALregBOOLUserVar(Stat, 	EstimSigParam, 	'p', UVAR_OPTIONAL, "Do Signal Parameter Estimation");
  LALregREALUserVar(Stat, 	Fthreshold,	'F', UVAR_OPTIONAL, "Signal Set the threshold for selection of 2F");
  LALregINTUserVar(Stat, 	windowsize,	'k', UVAR_OPTIONAL, "Running-Median window size");
  LALregINTUserVar(Stat, 	gridType,	 0 , UVAR_OPTIONAL, "Template grid: 0=flat, 1=isotropic, 2=metric, 3=file");
  LALregINTUserVar(Stat, 	metricType,	'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserVar(Stat, 	metricMismatch,	'X', UVAR_OPTIONAL, "Maximal mismatch for metric tiling");
  LALregBOOLUserVar(Stat, 	help, 		'h', UVAR_HELP,     "Print this message");
  LALregSTRINGUserVar(Stat,	skyRegion, 	'R', UVAR_OPTIONAL, "Specify sky-region by polygon");
  LALregSTRINGUserVar(Stat,	outputLabel,	'o', UVAR_OPTIONAL, "Label to be appended to all output file-names");
  LALregSTRINGUserVar(Stat,	outputFstat,	 0,  UVAR_OPTIONAL, "Output-file for the F-statistic field over the parameter-space");
  LALregSTRINGUserVar(Stat,	skyGridFile,	 0,  UVAR_OPTIONAL, "Load sky-grid from this file.");
  LALregSTRINGUserVar(Stat,	outputSkyGrid,	 0,  UVAR_OPTIONAL, "Write sky-grid into this file.");
  LALregSTRINGUserVar(Stat,	binarytemplatefile,	 0,  UVAR_OPTIONAL, "Read binary templates from this file.");
  LALregSTRINGUserVar(Stat,	sourcefile,	 0,  UVAR_OPTIONAL, "Read in source parameters from this file.");
  LALregSTRINGUserVar(Stat,	source,	 0,          UVAR_OPTIONAL, "Name of the binary source to be analysed .");
  LALregREALUserVar(Stat,	overres,	 0,          UVAR_OPTIONAL, "The frequency over resolution factor.");

  LALregBOOLUserVar(Stat,	openDX,	 	 0,  UVAR_OPTIONAL, "Make output-files openDX-readable (adds proper header)");
  LALregSTRINGUserVar(Stat,     workingDir,     'w', UVAR_OPTIONAL, "Directory to be made the working directory, . is default");


  DETATCHSTATUSPTR (Stat);
  RETURN (Stat);
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
  

#ifdef DEBG_ESTSIGPAR
  REAL8 Ftest;
  REAL8 A2test,A3test,A4test;
#endif


  strcpy(Paramfilename,"ParamMLE");
  if (uvar_outputLabel)
    strcat(Paramfilename,uvar_outputLabel);
  
  if(!(fpMLEParam=fopen(Paramfilename,"w")))
    fprintf(stderr,"Error in EstimateSignalParameters: unable to open the file");


  norm=2.0*sqrt(GV.tsft)/(GV.tsft*GV.SFTno);


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


#ifdef DEBG_ESTSIGPAR
      /* Reconstruct A1,A2,A3,A4. Compare them with the original values. */

      A1test=h0mle*(0.5*(1+mu_mle*mu_mle)*cos(2.0*psi_mle)*cos(2.0*Phi0_mle)
		    -mu_mle*sin(2.0*psi_mle)*sin(2.0*Phi0_mle));
      A2test=h0mle*(0.5*(1+mu_mle*mu_mle)*sin(2.0*psi_mle)*cos(2.0*Phi0_mle)
		    +mu_mle*cos(2.0*psi_mle)*sin(2.0*Phi0_mle));
      A3test=h0mle*(-0.5*(1+mu_mle*mu_mle)*cos(2.0*psi_mle)*sin(2.0*Phi0_mle)
		    -mu_mle*sin(2.0*psi_mle)*cos(2.0*Phi0_mle));
      A4test=h0mle*(-0.5*(1+mu_mle*mu_mle)*sin(2.0*psi_mle)*sin(2.0*Phi0_mle)
		    +mu_mle*cos(2.0*psi_mle)*cos(2.0*Phi0_mle));


      fprintf(stderr,"LALDemod_Estimate output: "
              "A1=%20.15f A2=%20.15f A3=%20.15f A4=%20.15f\n"
	      ,A1,A2,A3,A4);
      fprintf(stderr,"Reconstructed from MLE: "
              "A1=%20.15f A2=%20.15f A3=%20.15f A4=%20.15f !!!!\n\n",
	      A1test,A2test,A3test,A4test);
      fflush(stderr);


      if(fabs(A1-A1test)>fabs(A1)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A1\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A1=%f,A1test=%f\n",
		irec,uvar_Freq+irec*GV.dFreq,A1,A1test);
	fprintf(stderr,"relative error Abs((A1-A1test)/A1)=%lf\n",
		fabs(A1-A1test)/fabs(A1));
	exit(1);
      }
      if(fabs(A2-A2test)>fabs(A2)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A2\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A2=%f,A2test=%f\n",
		irec,uvar_Freq+irec*GV.dFreq,A2,A2test);
	fprintf(stderr,"relative error Abs((A2-A2test)/A2)=%lf\n",
		fabs(A2-A2test)/fabs(A2));
	exit(1);
      }
      if(fabs(A3-A3test)>fabs(A3)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A3\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A3=%f,A3test=%f\n",
		irec,uvar_Freq+irec*GV.dFreq,A3,A3test);
	fprintf(stderr,"relative error Abs((A3-A3test)/A3)=%lf\n",
		fabs(A3-A3test)/fabs(A3));
	exit(1);
      }
      if(fabs(A4-A4test)>fabs(A4)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A4\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A4=%f,A4test=%f\n",
		irec,uvar_Freq+irec*GV.dFreq,A1,A1test);
	fprintf(stderr,"relative error Abs((A4-A4test)/A4)=%lf\n",
		fabs(A4-A4test)/fabs(A4));
	exit(1);
      }

      
      /* Reconstruct F. Compare it with the original value. */
      Ftest=(A*(A1*A1+A3*A3)+B*(A2*A2+A4*A4)+2.0*C*(A1*A2+A3*A4))/
	4.0*4.0/GV.SFTno;


      if(fabs(Fstat.F[irec] - Ftest)> fabs(Ftest)/10e5){ 
	fprintf(stderr,"Something is wrong with Estimate in F\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),F=%f,Ftest=%f\n",
		irec,uvar_Freq+irec*GV.dFreq,Fstat.F[irec],Ftest);
	fprintf(stderr,"relative error Abs((F-Ftest)/Ftest)=%lf\n",
		fabs(Fstat.F[irec]-Ftest)/fabs(Ftest));
	exit(1);
      }
#endif


      /* normalization */
      h0mle=h0mle*norm;


      /* For the real data, we need to multiply long(2.0) */
      /* Because we use running median to estimate the S_h. */
      /* if(GV.SignalOnly!=1) 
	h0mle=h0mle*sqrt(medianbias);
      */
      /* medianbias is 1 when GV.SignalOnly==1 */
      h0mle=h0mle*sqrt(medianbias);

#ifdef DEBG_ESTSIGPAR
      {double hp,hc,ds;
      hp=(1.0+mu_mle*mu_mle)*h0mle/2.0;
      hc=mu_mle*h0mle;
      ds=GV.SFTno*GV.tsft/2.0*(hp*hp*((A+B)/2.0+(A-B)/2.0*cos(4.0*psi_mle)
				   +C*sin(4.0*psi_mle))
			    +hc*hc*((A+B)/2.0-(A-B)/2.0*cos(4.0*psi_mle)
				    -C*sin(4.0*psi_mle)));
      fprintf(stderr,"A=%f,B=%f,C=%f,f=%f,h0=%f,F=%f\n",
	      A,B,C,uvar_Freq+irec*GV.dFreq,h0mle,Fstat.F[irec]*medianbias);
      }
#endif

      /* Note that we print out MLE of 2.0*Phi0_JKS */
      /* because Phi0_PULGROUPDOC=2.0*Phi0_JKS */
      /* and Phi0_PULGROUPDOC is the one used in In.data. */
 
      /* medianbias is 1 if GV.SignalOnly==1 */
      fprintf(fpMLEParam,"%16.8f %22E", uvar_Freq + irec*GV.dFreq, 2.0*medianbias*Fstat.F[irec]);


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

#ifdef DEBG_FAFB
      fprintf(fp,"%22.16f %22.16f "
                 "%E %20.17f %20.17f "
                 "%22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n",
	      uvar_Freq+index*GV.dFreq,Fstat.F[index]*bias*bias,
	      DemodParams->spinDwn[0], Alpha, Delta,
	      Fstat.Fa[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fa[index].im/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].im/sqrt(GV.SFTno)*bias,
	      amc.A,amc.B,amc.C);
#else
      /* Freqency, Re[Fa],Im[Fa],Re[Fb],Im[Fb], F */
      fprintf(fp,"%22.16f %22.12f %22.12f %22.12f %22.12f %22.12f\n",
	      uvar_Freq+index*GV.dFreq,
	      Fstat.Fa[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fa[index].im/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].im/sqrt(GV.SFTno)*bias,
	      Fstat.F[index]*bias*bias);
#endif


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
  INT4 k;

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
       teemp += 0.5*GV.tsft;
       TRY (LALFloatToGPS(status->statusPtr, &(midTS[k]), &teemp), status);
     }
   
   TRY (LALComputeAM(status->statusPtr, &amc, midTS, amParams), status); 

/* ComputeSky stuff*/
  csParams=(CSParams *)LALMalloc(sizeof(CSParams));
  csParams->tGPS=timestamps;  
  csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
  csParams->mObsSFT=GV.SFTno;     /* Changed this from GV.mobssft !!!!!! */
  csParams->tSFT=GV.tsft;
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

  DemodParams->Dterms=uvar_dterms;
  DemodParams->ifmin=GV.ifmin;

  DemodParams->returnFaFb = uvar_EstimSigParam;

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
  INT4 k;

  INITSTATUS (status, "CreateBinaryDemodParams", rcsid);
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
       teemp += 0.5*GV.tsft;
       TRY (LALFloatToGPS(status->statusPtr, &(midTS[k]), &teemp), status);
     }
   
   TRY (LALComputeAM(status->statusPtr, &amc, midTS, amParams), status); 

/* ComputeSkyBinary stuff*/
  csbParams=(CSBParams *)LALMalloc(sizeof(CSBParams));
  csbParams->tGPS=timestamps;  
  csbParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
  csbParams->mObsSFT=GV.SFTno;     /* Changed this from GV.mobssft !!!!!! */
  csbParams->tSFT=GV.tsft;
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

  DemodParams->Dterms=uvar_dterms;
  DemodParams->ifmin=GV.ifmin;

  DemodParams->returnFaFb = uvar_EstimSigParam;

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
void AllocateMem(LALStatus *status)
{
  INT4 k;

  INITSTATUS (status, "AllocateMem", rcsid);
  ATTATCHSTATUSPTR (status);

  /* Allocate space for AMCoeffs */
  amc.a = NULL;
  amc.b = NULL;
  TRY (LALSCreateVector(status->statusPtr, &(amc.a), (UINT4) GV.SFTno), status);
  TRY (LALSCreateVector(status->statusPtr, &(amc.b), (UINT4) GV.SFTno), status);

  /* Allocate DemodParams structure */
  DemodParams=(DemodPar *)LALCalloc(1, sizeof(DemodPar));
  
  /* space for sky constants */
  /* Based on maximum index for array of as and bs sky constants as from ComputeSky.c */
  k=4*(GV.SFTno-1)+4; 
  DemodParams->skyConst = (REAL8 *)LALMalloc(k*sizeof(REAL8));

  /* space for spin down params */
  if (!uvar_binary) DemodParams->spinDwn = (REAL8 *)LALMalloc(sizeof(REAL8));
  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* AllocateMem() */

/*******************************************************************************/

/*  for every cluster writes the information about it in file Fstats */
/*  precisely it writes: */
/*  fr_max alpha delta N_points_of_cluster mean std max (of 2F) */
int writeFLines(INT4 *maxIndex){

  INT4 i,j,j1,j2,k,N;
  REAL8 max,log2,mean,var,std,R,fr;
  INT4 imax;
  INT4 err=0;

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
#ifdef FILE_FSTATS  
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
#endif

  }/*  end i loop over different clusters */

  return 0;
}

/*******************************************************************************/

int NormaliseSFTData(void)
{
  INT4 k,j;                         /* loop indices */
  INT4 nbins=GV.ifmax-GV.ifmin+1;   /* Number of points in SFT's */
  REAL8 SFTsqav;                  /* Average of Square of SFT */
  REAL8 B;                        /* SFT Bandwidth */
  REAL8 deltaT,N;


  /* loop over each SFTs */
  for (k=0;k<GV.SFTno;k++)         
    {
      
      SFTsqav=0.0;
      /* loop over SFT data to estimate noise */
      for (j=0;j<nbins;j++)               
	  {
	    SFTsqav=SFTsqav+
	      SFTData[k]->fft->data->data[j].re * SFTData[k]->fft->data->data[j].re+
	      SFTData[k]->fft->data->data[j].im * SFTData[k]->fft->data->data[j].im;
	  }
      SFTsqav=SFTsqav/(1.0*nbins);              /* Actual average of Square of SFT */
      MeanOneOverSh=2.0*GV.nsamples*GV.nsamples/(SFTsqav*GV.tsft)+MeanOneOverSh;      

      N=1.0/sqrt(2.0*(REAL8)SFTsqav);

      /* signal only case */  
      if(uvar_SignalOnly == 1)
	{
	  B=(1.0*GV.nsamples)/(1.0*GV.tsft);
	  deltaT=1.0/(2.0*B);
	  N=deltaT/sqrt(GV.tsft);
	}

      /* loop over SFT data to Normalise it*/
      for (j=0;j<nbins;j++)               
	{
	  SFTData[k]->fft->data->data[j].re = N*SFTData[k]->fft->data->data[j].re; 
	  SFTData[k]->fft->data->data[j].im = N*SFTData[k]->fft->data->data[j].im;
	}
    } 
	
   MeanOneOverSh=MeanOneOverSh/(1.0*GV.SFTno); 
  return 0;
}

/*******************************************************************************/

int ReadSFTData(void)
{
  INT4 fileno=0,offset;
  FILE *fp=NULL;
  size_t errorcode;
  UINT4 ndeltaf;
  INT4 k=0;

  SFTData=(FFT **)LALMalloc(GV.SFTno*sizeof(FFT *));
  timestamps=(LIGOTimeGPS *)LALMalloc(GV.SFTno*sizeof(LIGOTimeGPS));
  
  /* if using a merged SFT file, it's already open */
  if (uvar_mergedSFTFile)
    fp=fp_mergedSFT;

  for (fileno=0;fileno<GV.SFTno;fileno++)
    {
      /* seek to fileno'th SFT k bytes in */
      if (uvar_mergedSFTFile){
	if (fseek(fp,k,SEEK_SET)) {
	  fprintf(stderr,"Unable to seek to the start of %s !\n",GV.filelist[fileno]);
	  return 1;
	}
      }
      else {
	if (!(fp=fopen(GV.filelist[fileno],"rb"))) {
	  fprintf(stderr,"Weird... %s doesn't exist!\n",GV.filelist[fileno]);
	  return 1;
	}
      }

      /* Read in the header from the file */
      errorcode=fread((void*)&header,sizeof(header),1,fp);
      if (errorcode!=1) 
	{
	  fprintf(stderr,"No header in data file %s\n",GV.filelist[fileno]);
	  return 1;
	}

      /* Check that data is correct endian order */
      if (header.endian!=1.0)
	{
	  fprintf(stderr,"First object in file %s is not (double)1.0!\n",GV.filelist[fileno]);
	  fprintf(stderr,"It could be a file format error (big/little\n");
	  fprintf(stderr,"endian) or the file might be corrupted\n\n");
	  return 2;
	}
    
      /* Check that the time base is positive */
      if (header.tbase<=0.0)
	{
	  fprintf(stderr,"Timebase %f from data file %s non-positive!\n",
		  header.tbase,GV.filelist[fileno]);
	  return 3;
	}
        
      /* Check that are frequency bins needed are in data set */
      if (GV.ifmin<header.firstfreqindex || 
	  GV.ifmax>header.firstfreqindex+header.nsamples) 
	{
	  fprintf(stderr,"Freq index range %d->%d not in %d to %d (file %s)\n",
		GV.ifmin,GV.ifmax,header.firstfreqindex,
		  header.firstfreqindex+header.nsamples,GV.filelist[fileno]);
	  return 4;
	}
      /* Put time stamps from file into array */
      timestamps[fileno].gpsSeconds = header.gps_sec;
      timestamps[fileno].gpsNanoSeconds = header.gps_nsec;

      /* Move forward in file */
      offset=(GV.ifmin-header.firstfreqindex)*2*sizeof(REAL4);
      errorcode=fseek(fp,offset,SEEK_CUR);
      if (errorcode) 
	{
	  perror(GV.filelist[fileno]);
	  fprintf(stderr,"Can't get to offset %d in file %s\n",offset,GV.filelist[fileno]);
	  return 5;
	}

      /* Make data structures */
      ndeltaf=GV.ifmax-GV.ifmin+1;
      SFTData[fileno]=(FFT *)LALMalloc(sizeof(FFT));
      SFTData[fileno]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
      SFTData[fileno]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
      SFTData[fileno]->fft->data->data=(COMPLEX8 *)LALMalloc(ndeltaf*sizeof(COMPLEX8));

      /* Fill in actual SFT data, and housekeeping */
      errorcode=fread((void*)(SFTData[fileno]->fft->data->data), sizeof(COMPLEX8), ndeltaf, fp);
      if (errorcode!=ndeltaf){
	perror(GV.filelist[fileno]);
	fprintf(stderr, "The SFT data was truncated.  Only read %d not %d complex floats\n", errorcode, ndeltaf);
	return 6;
      }
      SFTData[fileno]->fft->epoch=timestamps[fileno];
      SFTData[fileno]->fft->f0 = GV.ifmin / GV.tsft;
      SFTData[fileno]->fft->deltaF = 1.0 / GV.tsft;
      SFTData[fileno]->fft->data->length = ndeltaf;

         if (uvar_mergedSFTFile)
	   k+=sizeof(header)+header.nsamples*8;
	 else
	   fclose(fp);     /* close file */
    
    }
  return 0;  
}

/*******************************************************************************/

void
SetGlobalVariables(LALStatus *status, ConfigVariables *cfg)
{

  CHAR command[512];
  FILE *fp;
  size_t errorcode;
  INT4 fileno=0;   
#ifndef NOGLOB
  glob_t globbuf;
#endif
  LIGOTimeGPS starttime;

  INITSTATUS (status, "SetGlobalVariables", rcsid);
  ATTATCHSTATUSPTR (status);

  /* do some sanity checks on the user-input before we proceed */
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

  /* don't allow negative bands (for safty in griding-routines) */
  if ( (uvar_AlphaBand < 0) ||  (uvar_DeltaBand < 0) )
    {
      LALPrintError ("\nNegative value of sky-bands not allowed (alpha or delta)!\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }

  /* set the current working directory */
#ifndef _MSC_VER
  if(chdir(uvar_workingDir) != 0)
#else
  if(_chdir(uvar_workingDir) != 0)
#endif
    {
      fprintf(stderr, "in Main: unable to change directory to %s\n", uvar_workingDir);
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }

#if USE_BOINC
  strcat(cfg->EphemEarth,"earth");
  strcat(cfg->EphemSun,"sun");
  use_boinc_filename0(cfg->EphemEarth);
  use_boinc_filename0(cfg->EphemSun);
#else 
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
#endif /* !USE_BOINC */

  if (uvar_mergedSFTFile) {
    long k=0;
#if USE_BOINC
    use_boinc_filename1(&(uvar_mergedSFTFile));
#endif
    if (!(fp=fp_mergedSFT=fopen(uvar_mergedSFTFile,"rb"))){
      fprintf(stderr,"Unable to open SFT file %s\n", uvar_mergedSFTFile);
      ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
    }
    
    while (fread((void*)&header,sizeof(header),1,fp) == 1) {
      char tmp[256];
      /* check that we've still got space for more data */
      if (fileno >= MAXFILES) {
	fprintf(stderr,"Too many SFT's in merged file! Exiting... \n");
	ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
      }
      
      /* store a "file name" composed of merged name + block number */
      sprintf(tmp, "%s (block %d)", uvar_mergedSFTFile, fileno+1);
      strcpy(cfg->filelist[fileno],tmp);
      
      /* check that data is correct endian order */
      if (header.endian!=1.0) {
	fprintf(stderr,"First object in file %s is not (double)1.0!\n",cfg->filelist[fileno]);
	fprintf(stderr,"It could be a file format error (big/little\n");
	fprintf(stderr,"endian) or the file might be corrupted\n\n");
	ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
      }
      
      /* if this is the first SFT, save initial time */
      if (fileno==0) {
	cfg->Ti=header.gps_sec;
	starttime.gpsSeconds = header.gps_sec;
	starttime.gpsNanoSeconds = header.gps_nsec;
      }
      /* increment file no and pointer to the start of the next header */
      fileno++;
      k=header.nsamples*8;
      fseek(fp,k,SEEK_CUR);
    }
    /* save final time and time baseline */
    cfg->Tf=header.gps_sec+header.tbase;  /* FINAL TIME */
    cfg->tsft=header.tbase;  /* Time baseline of SFTs */
    
    /* NOTE: we do NOT close fp here.  If we are using merged SFT file
       for data, we keep it open because we'll need it again in
       ReadSFTData().
    */
    
  } 
  else {

    strcpy(command, uvar_DataDir);
    strcat(command,"/*");
    
    strcat(command, uvar_BaseName);
    strcat(command,"*");
    
#ifndef NOGLOB
    globbuf.gl_offs = 1;
    glob(command, GLOB_ERR, NULL, &globbuf);
    
    /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */
    
    if(globbuf.gl_pathc==0)
      {
	LALPrintError ("\nNo SFTs in directory %s ... Exiting.\n\n", uvar_DataDir);
	ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
      }
    
    while ((UINT4)fileno < globbuf.gl_pathc) 
      {
	strcpy(cfg->filelist[fileno],globbuf.gl_pathv[fileno]);
	fileno++;
	if (fileno > MAXFILES)
	  {
	    LALPrintError ("\nToo many files in directory! Exiting... \n\n");
	    ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
	  }
      }
    globfree(&globbuf);
#endif
  }
  cfg->SFTno=fileno; /* remember this is 1 more than the index value */

  /* initialize detector */
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

  if (!uvar_mergedSFTFile) {
    /* open FIRST file and get info from it*/
    fp=fopen(cfg->filelist[0],"rb");
    /* read in the header from the file */
    errorcode=fread((void*)&header,sizeof(header),1,fp);
    if (errorcode!=1) 
      {
	LALPrintError ("\nNo header in data file %s\n\n", cfg->filelist[0]);
	ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
      }
    
    /* check that data is correct endian order */
    if (header.endian!=1.0)
      {
	fprintf(stderr,"First object in file %s is not (double)1.0!\n",cfg->filelist[0]);
	fprintf(stderr,"It could be a file format error (big/little\n");
	fprintf(stderr,"endian) or the file might be corrupted\n\n");
	ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
      }
    fclose(fp);
    
    /* INITIAL TIME */
    starttime.gpsSeconds = header.gps_sec;
    starttime.gpsNanoSeconds = header.gps_nsec;
    cfg->Ti = header.gps_sec; 
    
    /* open LAST file and get info from it*/
    fp=fopen(cfg->filelist[fileno-1],"rb");
    /* read in the header from the file */
    errorcode=fread((void*)&header,sizeof(header),1,fp);
    if (errorcode!=1) 
      {
	fprintf(stderr,"No header in data file %s\n",cfg->filelist[fileno-1]);
	ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
      }
    /* check that data is correct endian order */
    if (header.endian!=1.0)
      {
	fprintf(stderr,"First object in file %s is not (double)1.0!\n",cfg->filelist[fileno-1]);
	fprintf(stderr,"It could be a file format error (big/little\n");
	fprintf(stderr,"endian) or the file might be corrupted\n\n");
	ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
      }
    fclose(fp);

    cfg->Tf=header.gps_sec+header.tbase;  /* FINAL TIME */

    cfg->tsft=header.tbase;  /* Time baseline of SFTs */
  }
    
  if (!uvar_binary)
    {
      /* if user has not input demodulation frequency resolution; set to 1/2*Tobs */
      if( !LALUserVarWasSet (&uvar_dFreq) ) 
	cfg->dFreq=1.0/(2.0*header.tbase*cfg->SFTno);
      else
	cfg->dFreq = uvar_dFreq;
      
      cfg->FreqImax=(INT4)(uvar_FreqBand/cfg->dFreq+.5)+1;  /*Number of frequency values to calculate F for */
      
      /* if user has not input spin down increment then set it to some value */
      if( !LALUserVarWasSet (&uvar_df1dot) ) 
	uvar_df1dot=1.0/(2.0*header.tbase*cfg->SFTno*(cfg->Tf - cfg->Ti));
    }
  if (uvar_binary)  /* BINARY-MOD - Just a safety thing, setting default df=1/5T and df1dot=0.0 */
    {
      /* if user has not input demodulation frequency resolution; set to 1/5*Tspan */
      if(( !LALUserVarWasSet (&uvar_dFreq) )&(!LALUserVarWasSet(&uvar_overres))) 
	cfg->dFreq=1.0/(5.0*header.tbase*(cfg->Tf-cfg->Ti));  /* note here we use tspan !!!! */
      else if (LALUserVarWasSet(&uvar_overres))
	cfg->dFreq=1.0/(uvar_overres*header.tbase*(cfg->Tf-cfg->Ti));
      else
	cfg->dFreq = uvar_dFreq;
      
      cfg->FreqImax=(INT4)(uvar_FreqBand/cfg->dFreq+.5)+1;  /*Number of frequency values to calculate F for */
      
      /* if user has not input spin down increment then set it to zero (default safety value ) */
      if( !LALUserVarWasSet (&uvar_df1dot) ) 
	uvar_df1dot=0.0;

    }

  if (LALUserVarWasSet (&uvar_f1dotBand) && (uvar_f1dotBand != 0) )
    cfg->SpinImax=(int)(uvar_f1dotBand/uvar_df1dot+.5)+1;  /*Number of spindown values to calculate F for */
  else
    cfg->SpinImax = 0;

  cfg->nsamples=header.nsamples;    /* # of freq. bins */

  cfg->ifmax=ceil((1.0+uvar_dopplermax)*(uvar_Freq+uvar_FreqBand)*cfg->tsft)+uvar_dterms;
  cfg->ifmin=floor((1.0-uvar_dopplermax)*uvar_Freq*cfg->tsft)-uvar_dterms;

  /* allocate F-statistic arrays */
  Fstat.F =(REAL8*)LALMalloc(cfg->FreqImax*sizeof(REAL8));
  if(uvar_EstimSigParam) 
    {
      Fstat.Fa =(COMPLEX16*)LALMalloc(cfg->FreqImax*sizeof(COMPLEX16));
      Fstat.Fb =(COMPLEX16*)LALMalloc(cfg->FreqImax*sizeof(COMPLEX16));
    } else {
      Fstat.Fa = NULL;
      Fstat.Fb = NULL;
    }

  /*----------------------------------------------------------------------
   * initialize+check  template-grid related parameters 
   */
  if (!uvar_binary)
    {
      UINT2 haveSkyRegion, haveAlphaDelta, haveGridFile, needGridFile, haveMetric, needMetric, setMetric;

      haveSkyRegion  = (uvar_skyRegion != NULL);
      haveAlphaDelta = (LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta) );
      haveGridFile   = (uvar_skyGridFile != NULL);
      needGridFile   = (uvar_gridType == GRID_FILE);
      haveMetric     = (uvar_metricType > LAL_PMETRIC_NONE);
      needMetric     = (uvar_gridType == GRID_METRIC);
      setMetric      = LALUserVarWasSet (&uvar_metricType);
      
      /* some consistency checks on input to help catch errors */
      if ( !needGridFile && !(haveSkyRegion || haveAlphaDelta) )
	{
	  LALPrintError ("\nNeed sky-region: either use (Alpha,Delta) or skyRegion!\n\n");
	ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
	}
      if ( !needGridFile && haveSkyRegion && haveAlphaDelta )
	{
	  LALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta) OR skyRegion!\n\n");
	  ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
	}
      if ( needGridFile && !haveGridFile )
	{
	  LALPrintError ("\nERROR: gridType=FILE, but no skyGridFile specified!\n\n");
	  ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);	
	}
      if ( !needGridFile && haveGridFile )
	{
	  LALWarning (status, "\nWARNING: skyGridFile was specified but not needed ... will be ignored\n");
	}
      if ( needGridFile && (haveSkyRegion || haveAlphaDelta) )
	{
	  LALWarning (status, "\nWARNING: We are using skyGridFile, but sky-region was also specified ... will be ignored!\n");
	}
      if ( setMetric && haveMetric && !needMetric) 
	{
	  LALWarning (status, "\nWARNING: Metric was specified for non-metric grid... will be ignored!\n");
	}
      if (needMetric && !haveMetric) 
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
      
      
    } /* end if not binary: template-grid stuff */


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

    TRY (LALLeapSecs (status->statusPtr, &leap, &starttime, &formatAndAcc), status);
    cfg->edat->leap = leap;

    TRY (LALInitBarycenter(status->statusPtr, cfg->edat), status);               
  }
  /*------------------------------------------------------------*/

  /* Tell the user what we have arrived at */
#ifdef DEBG_SGV
    fprintf(stdout,"\n");
    fprintf(stdout,"# SFT time baseline:                  %f min\n",header.tbase/60.0);
    fprintf(stdout,"# SFT freq resolution:                %f Hz\n",df);
    fprintf(stdout,"# Starting search frequency:          %f Hz\n",uvar_Freq);
    fprintf(stdout,"# Demodulation frequency band:        %f Hz\n",uvar_FreqBand);
    fprintf(stdout,"# no of SFT in a DeFT:                %f\n",ceil((1.0*(cfg->Tf - cfg->Ti))/header.tbase));
    fprintf(stdout,"# Actual # of SFTs:                   %d\n",cfg->SFTno);
    fprintf(stdout,"# ==> DeFT baseline:                  %f hours\n",(cfg->Tf - cfg->Ti)/3600.0);
#endif

    DETATCHSTATUSPTR (status);
    RETURN (status);

} /* SetGlobalVariables() */

/*******************************************************************************/



/*******************************************************************************/

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

void Freemem(LALStatus *status) 
{

  INT4 k;

  INITSTATUS (status, "Freemem", rcsid);
  ATTATCHSTATUSPTR (status);

  /* Free SFTData */
  for (k=0;k<GV.SFTno;k++)
    {
      LALFree(SFTData[k]->fft->data->data);
      LALFree(SFTData[k]->fft->data);
      LALFree(SFTData[k]->fft);
      LALFree(SFTData[k]);
    }
  LALFree(SFTData);

  /* Free timestamps */
  LALFree(timestamps);

  LALFree(Fstat.F);
  if(uvar_EstimSigParam) 
    {
      LALFree(Fstat.Fa);
      LALFree(Fstat.Fb);
    }

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
    /*LALFree(BinaryBank->BMFheader);*/
    LALFree(BinaryBank);
  }

  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (status->statusPtr), status);


  if (!uvar_binary) 
    {
      /* Free DopplerScan-stuff (grid) */
      TRY (FreeDopplerScan (status->statusPtr, &thisScan), status);
      if (GV.skyRegion)
	LALFree ( GV.skyRegion );
    }
      
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

/* Sorting function to sort into DECREASING order */
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

/*******************************************************************************/

/* This routine prints the values of (f,FF) above a certain threshold */
/* in 2F, called 2Fthr.  If there are more than ReturnMaxN of these, */
/* then it simply returns the top ReturnMaxN of them. If there are */
/* none, then it returns none.  It also returns some basic statisical */
/* information about the distribution of 2F: the mean and standard */
/* deviation. */
/* Returns zero if all is well, else nonzero if a problem was encountered. */
/*    Basic strategy: sort the array by values of F, then look at the */
/*    top ones. Then search for the points above threshold. */


INT4 PrintTopValues(REAL8 TwoFthr, INT4 ReturnMaxN)
{

  INT4 *indexes,i,j,iF,ntop,err,N;
  REAL8 mean=0.0, std=0.0 ,log2 /*=log(2.0)*/;

  log2=medianbias;



  for (i=0;i<highFLines->Nclusters;i++){
    N=highFLines->NclustPoints[i];
    for (j=0;j<N;j++){
      iF=highFLines->Iclust[j];
      Fstat.F[iF]=0.0;
    }/*  end j loop over points of i-th cluster  */
  }/*  end i loop over different clusters */




/*    check that ReturnMaxN is sensible */
 
  if (ReturnMaxN>GV.FreqImax){
    fprintf(stderr,"PrintTopValues() WARNING: resetting ReturnMaxN=%d to %d\n",
	    ReturnMaxN, (INT4)GV.FreqImax);
    ReturnMaxN=GV.FreqImax;
  }

/*    create an array of indexes */
  if (!(indexes=(INT4 *)LALMalloc(sizeof(INT4)*GV.FreqImax))){
    fprintf(stderr,"Unable to allocate index array in PrintTopValues()\n");
    return 1;
  }

/*    populate it */
  for (i=0;i<GV.FreqImax;i++)
    indexes[i]=i;

/*   sort array of indexes */
  qsort((void *)indexes, (size_t)GV.FreqImax, sizeof(int), compare);


/*    Normalize */
  TwoFthr*=0.5/log2;

#ifdef FILE_FMAX  
#if 1
  if (!uvar_binary) 
    {
      /*    print out the top ones */
      for (ntop=0; ntop<ReturnMaxN; ntop++)
	if (Fstat.F[indexes[ntop]]>TwoFthr){
	  err=fprintf(fpmax, "%20.10f %10.8f %10.8f %20.15f\n",
		      uvar_Freq+indexes[ntop]*GV.dFreq,
		      Alpha, Delta, 2.0*log2*Fstat.F[indexes[ntop]]);
	  if (err<=0) {
	    fprintf(stderr,"PrintTopValues couldn't print to Fmax!\n");
	    LALFree(indexes);
	    return 3;
	  }
	}
	else
	  /*        Since array sorted, as soon as value too small we can break */
	  break;
    }
  if (uvar_binary) 
    {
      /*    print out the top ones */
      for (ntop=0; ntop<ReturnMaxN; ntop++)
	if (Fstat.F[indexes[ntop]]>TwoFthr){
	  err=fprintf(fpmax, "%20.10f %10.8f %10.8f %10.8f %10.8f %d %d %10.8f %10.8f %20.15f\n",
		      uvar_Freq+indexes[ntop]*GV.dFreq,
		      Alpha, Delta,thisBinaryTemplate.ProjSMaxis,thisBinaryTemplate.Period,
		      thisBinaryTemplate.TperiSSB.gpsSeconds,thisBinaryTemplate.TperiSSB.gpsNanoSeconds,
		      thisBinaryTemplate.Eccentricity,thisBinaryTemplate.ArgPeri,2.0*log2*Fstat.F[indexes[ntop]]);
	  if (err<=0) {
	    fprintf(stderr,"PrintTopValues couldn't print to Fmax!\n");
	    LALFree(indexes);
	    return 3;
	  }
	}
	else
	  /*        Since array sorted, as soon as value too small we can break */
	  break;
    }
#endif
#endif

  /*  find out how many points have been set to zero (N) */
  N=0;
  if (highFLines) {
    for (i=0;i<highFLines->Nclusters; i++){
      N=N+highFLines->NclustPoints[i];
    }
  }
/*    get mean of F[] */
  for (i=0;i<GV.FreqImax; i++)
    mean+=Fstat.F[i];
  mean/=(GV.FreqImax-N);

/*    get sigma for F[] */
  for (i=0; i<GV.FreqImax; i++){
    REAL8 diff=Fstat.F[i]-mean;
    std+=diff*diff;
  }
  std/=(GV.FreqImax-1-N);

/*    normalize by appropriate factors */
  mean*=(2.0*log2);
  std=2.0*log2*sqrt(std);


/*    Find number of values above threshold (could be made O(log N) */
/*    with binary search! */
  for (i=0; i<GV.FreqImax; i++)
    if (Fstat.F[indexes[i]]<=TwoFthr)
      break;

#ifdef FILE_FMAX_DEBG    
/*    print the output */
  if (!uvar_binary) { 
  err=fprintf(fpmax,"%10.5f %10.8f %10.8f    %d %10.5f %10.5f %10.5f\n",uvar_Freq,
	      Alpha, Delta, GV.FreqImax-N, mean, std, 2.0*log2*Fstat.F[indexes[0]]);
  }
  if (uvar_binary) {
    err=fprintf(fpmax,"%10.5f %10.8f %10.8f %10.8f %10.8f %d %d %10.8f %10.8f %d %10.5f %10.5f %10.5f\n",
		uvar_Freq,Alpha, Delta,thisBinaryTemplate.SMaxis,thisBinaryTemplate.Period,
		thisBinaryTemplate.TperiSSB.gpsSeconds,thisBinaryTemplate.TperiSSB.gpsNanoSeconds,
		thisBinaryTemplate.Eccentricity,thisBinaryTemplate.ArgPeri, 
		GV.FreqImax-N, mean, std, 2.0*log2*Fstat.F[indexes[0]]);
  }
#endif
  LALFree(indexes);
#ifdef FILE_FMAX_DEBG    
  if (err<=0) {
    fprintf(stderr,"PrintTopValues couldn't print to Fmax!\n");
    return 4;
  }
#endif

  return 0;
}

/*******************************************************************************/

INT4 EstimatePSDLines(LALStatus *status)
{
#ifdef FILE_PSD
  FILE *outfile;
#endif
#ifdef FILE_PSDLINES
  FILE *outfile1;
#endif
  INT4 i,j,Ntot;                         /* loop indices */
  INT4 nbins=GV.ifmax-GV.ifmin+1;   /* Number of points in SFT's */
  REAL8Vector *Sp=NULL; 
  REAL8Vector *FloorSp=NULL;                        /* Square of SFT */
  INT2 windowSize=uvar_windowsize;                  /* Running Median Window Size*/
  REAL4 THR=10000.0;
  
  REAL4 xre,xim;

  OutliersInput  *outliersInput;
  OutliersParams *outliersParams;
  Outliers       *outliers;
  ClustersInput  *clustersInput;
  ClustersParams *SpClParams;
  Clusters       *SpLines=highSpLines;
    
  INT2 smallBlock=3;
  INT4 wings;

  nbins=(UINT4)nbins;
  wings=windowSize/2;

#ifdef FILE_PSD
  /*  file contains freq, PSD, noise floor */
  if(!(outfile=fopen("PSD.txt","w"))){
    printf("Cannot open PSD.txt file");
    return 1;
  } 
#endif 
#ifdef FILE_PSDLINES
  /*  file contains freq, PSD, noise floor,lines */
  if(!(outfile1=fopen("PSDLines.txt","w"))){
    printf("Cannot open PSD.txt file");
    return 1;
  }
#endif

  /* Allocate memory for input & output */
  /* if (!(Sp = (double *) calloc(nbins,sizeof(double)))){ */
  /*   printf("Memory allocation failure"); */
  /*   return 0; */
  /* } */
  
  LALDCreateVector(status, &Sp, nbins);
  LALDCreateVector(status, &FloorSp, nbins);
  
  
  /* loop over each SFTs */
  for (i=0;i<GV.SFTno;i++){
    
    /* loop over SFT data to estimate noise */
    for (j=0;j<nbins;j++){
      xre=SFTData[i]->fft->data->data[j].re;
      xim=SFTData[i]->fft->data->data[j].im;
      Sp->data[j]=Sp->data[j]+(REAL8)(xre*xre+xim*xim);
    }
  }/*end loop over SFTs*/
  
  /*Average Sp*/
  for (j=0;j<nbins;j++){
    Sp->data[j]=Sp->data[j]/GV.SFTno;
  }
  Sp->length=nbins;
  FloorSp->length=nbins;

  
  j=EstimateFloor(Sp, windowSize, FloorSp);
 
  if (!(outliers=(Outliers *)LALMalloc(sizeof(Outliers)))){
    fprintf(stderr,"Memory allocation failure for SpOutliers\n");
    return 1;
  }
  outliers->Noutliers=0;

  if (!(outliersParams=(OutliersParams *)LALMalloc(sizeof(OutliersParams)))){
    fprintf(stderr,"Memory allocation failure for OutliersParams\n");
    return 1;
  }
  if (!(outliersInput=(OutliersInput *)LALMalloc(sizeof(OutliersInput)))){
    fprintf(stderr,"Memory allocation failure for OutliersParams\n");
    return 1;
  }
  
  outliersParams->Thr=THR;
  outliersParams->Floor = FloorSp;
  outliersParams->wings=wings; /*these must be the same as ClustersParams->wings */
  outliersInput->ifmin=GV.ifmin;
  outliersInput->data = Sp;

  ComputeOutliers(outliersInput, outliersParams, outliers);

   if (outliers->Noutliers == 0){

#ifdef FILE_PSD
     /*  PSD.txt file contains freq, PSD, noise floor   */
     for (i=0;i<nbins;i++){ 
       REAL4 freq;
       REAL8 r0,r1;
       freq=(GV.ifmin+i)/GV.tsft;
       r0=Sp->data[i];
       r1=FloorSp->data[i];
       fprintf(outfile,"%f %E %E\n",freq,r0,r1);
     }
#endif

#ifdef FILE_PSD     
     fclose(outfile);
#endif
#ifdef FILE_PSDLINES
     fclose(outfile1);
#endif

     LALFree(outliers->ratio);
     LALFree(outliers);
     LALFree(outliersParams);
     LALFree(outliersInput);
     LALDDestroyVector(status,&Sp);
     LALDDestroyVector(status,&FloorSp);

     return 0;

   }
  


   if (!(SpClParams=(ClustersParams *)LALMalloc(sizeof(ClustersParams)))){ 
     printf("Memory allocation failure for SpClusterParams");
     return 1;
   }

   if (!(clustersInput=(ClustersInput *)LALMalloc(sizeof(ClustersInput)))){ 
     printf("Memory allocation failure for SpClusters");
     return 1;
   }
      
   SpClParams->wings=wings;
   SpClParams->smallBlock=smallBlock;
   
   clustersInput->outliersInput = outliersInput;
   clustersInput->outliersParams= outliersParams;
   clustersInput->outliers      = outliers;     
   
   j=DetectClusters(clustersInput, SpClParams, SpLines);
   if (j!=0){
     printf("DetectClusters problem");
     return 1;
   }
      
   /*  sum of points in all lines */
   Ntot=0;
   for (i=0;i<SpLines->Nclusters;i++){ 
     Ntot=Ntot+SpLines->NclustPoints[i];
   }
   
#ifdef FILE_PSDLINES
   /*  PSDLines file contains: PSD, noise floor and lines. */
   for (i=0;i<Ntot;i++){ 
     REAL4 freq;
     REAL8 r0,r1,r2;
     j=SpLines->Iclust[i];
     freq=(GV.ifmin+SpLines->Iclust[i])/GV.tsft;
     r0=Sp->data[j];
     r1=FloorSp->data[j];
     r2=SpLines->clusters[i]*FloorSp->data[j];
     fprintf(outfile1,"%f %E %E %E\n",freq,r0,r1,r2);
   }
#endif

#ifdef FILE_PSD   
   /*  PSD.txt file contains freq, PSD, noise floor   */
   for (i=0;i<nbins;i++){ 
     REAL4 freq;
     REAL8 r0,r1;
     freq=(GV.ifmin+i)/GV.tsft;
     r0=Sp->data[i];
     r1=FloorSp->data[i];
     fprintf(outfile,"%f %E %E\n",freq,r0,r1);
   }
#endif

#ifdef FILE_PSD   
   fclose(outfile);
#endif
#ifdef FILE_PSDLINES
   fclose(outfile1);
#endif

   LALFree(outliers->ratio);
   LALFree(outliers->outlierIndexes);
   LALFree(outliers);
   LALFree(outliersParams);
   LALFree(outliersInput);
   LALDDestroyVector(status,&Sp);
   LALDDestroyVector(status,&FloorSp);
   LALFree(SpClParams);
   LALFree(clustersInput);

   return 0;
}

/*******************************************************************************/

INT4 EstimateFLines(LALStatus *status) /* BINARY-MOD - may only need to change wings size */
{
#ifdef FILE_FTXT  
  FILE *outfile;
#endif
#ifdef FILE_FLINES  
  FILE *outfile1;
#endif
  INT4 i,j,Ntot;                         /* loop indices */
  INT4 nbins=GV.FreqImax;                /* Number of points in F */
  REAL8Vector *F1=NULL; 
  REAL8Vector *FloorF1=NULL;                        /* Square of SFT */
  /* INT2 windowSize=(0.01/GV.dFreq);               0.1 is 1E-4*1000 */
  INT2 windowSize=100;    /* this window size is fine for Fstat as we have high resolution */
  REAL4 THR;
  
  REAL8 dmp;

  OutliersInput  *outliersInput;
  OutliersParams *outliersParams;
  Outliers       *outliers;
  ClustersInput  *clustersInput;
  ClustersParams *SpClParams;
  Clusters       *SpLines=highFLines;
    
  INT2 smallBlock=1;
  
  INT4 wings;

  nbins=(UINT4)nbins;

  THR=uvar_Fthreshold;


  /* wings=windowSize/2; */
  /*  0.0002 is the max expected width of the F stat curve for signal */
  /*  with ~ 10 h observation time */
  /*  0.0001 = 0.0002/2 */
  /*  let me put 0.005 */
 
  /*dmp=0.5+0.0002/GV.dFreq; */  /* ISOLATED value */
  dmp=20;  /* BINARY-MOD - value chosen based on MC simulation.  No signals were found to have width FWHM > 20 search bins  */
  wings=dmp;


  if (windowSize > nbins){
    windowSize = nbins/2.0;
    /* printf("Had to change windowSize for running median in F floor estimate\n"); */
  }

#ifdef FILE_FTXT
  /*  file contains freq, PSD, noise floor */
  if(!(outfile=fopen("F.txt","w"))){
    printf("Cannot open F.txt file\n");
    return 1;
  }
#endif
#ifdef FILE_FLINES  
  /*  file contains freq, PSD, noise floor,lines */
  if(!(outfile1=fopen("FLines.txt","w"))){
    printf("Cannot open FLines.txt file\n");
    return 1;
  }
#endif


  LALDCreateVector(status, &F1, nbins);
  LALDCreateVector(status, &FloorF1, nbins);
    
  /* loop over SFT data to estimate noise */
  for (j=0;j<nbins;j++){
    F1->data[j]=Fstat.F[j];
    FloorF1->data[j]=1.0;
  }
  
  F1->length=nbins;
  FloorF1->length=nbins;

  /*   j=EstimateFloor(F1, windowSize, FloorF1); */
 
  if (!(outliers=(Outliers *)LALMalloc(sizeof(Outliers)))){
    fprintf(stderr,"Memory allocation failure for SpOutliers\n");
    return 1;
  }
  outliers->Noutliers=0;

  if (!(outliersParams=(OutliersParams *)LALMalloc(sizeof(OutliersParams)))){
    fprintf(stderr,"Memory allocation failure for OutliersParams\n");
    return 1;
  }
  if (!(outliersInput=(OutliersInput *)LALMalloc(sizeof(OutliersInput)))){
    fprintf(stderr,"Memory allocation failure for OutliersParams\n");
    return 1;
  }
  
  outliersParams->Thr=THR/(2.0*medianbias);
  outliersParams->Floor = FloorF1;
  outliersParams->wings=wings; /*these must be the same as ClustersParams->wings */
  outliersInput->ifmin=((uvar_Freq/GV.dFreq)+0.5);
  outliersInput->data = F1;

  ComputeOutliers(outliersInput, outliersParams, outliers);

   if (outliers->Noutliers == 0){

#ifdef FILE_FTXT
     if (!uvar_binary) {
       /*  F.txt file contains freq, F, noise floor of F   */
       for (i=0;i<nbins;i++){ 
	 REAL4 freq;
	 REAL8 r0,r1;
	 freq=uvar_Freq + i*GV.dFreq;
	 r0=F1->data[i];
	 r1=FloorF1->data[i];
	 fprintf(outfile,"%12.17f %E %E\n",freq,r0,r1);
       }
     }
     if (uvar_binary) {
       /*  F.txt file contains freq, F, noise floor of F   */
       for (i=0;i<nbins;i++){ 
	 REAL4 freq;
	 REAL8 r0,r1;
	 freq=uvar_Freq + i*GV.dFreq;
	 r0=F1->data[i];
	 r1=FloorF1->data[i];
	 fprintf(outfile1,"%20.17f %12.12f %12.12f %d %d %6.12f %6.12f %E %E\n",freq, \
		 thisBinaryTemplate.ProjSMaxis, \
		 thisBinaryTemplate.Period, \
		 thisBinaryTemplate.TperiSSB.gpsSeconds, \
		 thisBinaryTemplate.TperiSSB.gpsNanoSeconds, \
		 thisBinaryTemplate.Eccentricity, \
		 thisBinaryTemplate.ArgPeri,r0,r1);
       }
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
     LALDDestroyVector(status,&F1);
     LALDDestroyVector(status,&FloorF1);

     /*      fprintf(stderr,"Nclusters zero \n"); */
     /*      fflush(stderr); */

     return 0;

   }
  
   /*      fprintf(stderr,"Nclusters non zero \n"); */
   /*      fflush(stderr); */

   if (!(SpClParams=(ClustersParams *)LALMalloc(sizeof(ClustersParams)))){ 
     printf("Memory allocation failure for SpClusterParams");
     return 1;
   }

   
   if (!(clustersInput=(ClustersInput *)LALMalloc(sizeof(ClustersInput)))){ 
     printf("Memory allocation failure for SpClusters");
     return 1;
   }
      
   SpClParams->wings=wings;
   SpClParams->smallBlock=smallBlock;
   
   clustersInput->outliersInput = outliersInput;
   clustersInput->outliersParams= outliersParams;
   clustersInput->outliers      = outliers;     
   
   j=DetectClusters(clustersInput, SpClParams, SpLines);
   if (j!=0){
     printf("DetectClusters problem");
     return 1;
   }
   
   
   /*  sum of points in all lines */
   Ntot=0;
   for (i=0;i<SpLines->Nclusters;i++){ 
     Ntot=Ntot+SpLines->NclustPoints[i];
   }
   


#ifdef FILE_FLINES
   if (!uvar_binary)   /* BINARY-MOD - Need to select isolated or binary output */
     {
       /*  FLines file contains: F, noise floor and lines. */
       for (i=0;i<Ntot;i++){ 
	 REAL4 freq;
	 REAL8 r0,r1,r2;
	 j=SpLines->Iclust[i];
	 freq=(uvar_Freq+SpLines->Iclust[i]*GV.dFreq);
	 r0=F1->data[j];
	 r1=FloorF1->data[j];
	 r2=SpLines->clusters[i]*FloorF1->data[j];
	 fprintf(outfile1,"%f %12.12f %12.12f %E %E %E\n",freq,Alpha,Delta,r0,r1,r2);
       }
     }
   if (uvar_binary)
     {
       /*  FLines file contains: F, noise floor and lines. */
       for (i=0;i<Ntot;i++){ 
	 REAL4 freq;
	 REAL8 r0,r1,r2;
	 j=SpLines->Iclust[i];
	 freq=(uvar_Freq+SpLines->Iclust[i]*GV.dFreq);
	 r0=F1->data[j];
	 r1=FloorF1->data[j];
	 r2=SpLines->clusters[i]*FloorF1->data[j];
	 fprintf(outfile1,"%20.17f %12.12f %12.12f %d %d %6.12f %6.12f %E %E %E\n",freq, \
		 thisBinaryTemplate.ProjSMaxis, \
		 thisBinaryTemplate.Period, \
		 thisBinaryTemplate.TperiSSB.gpsSeconds, \
		 thisBinaryTemplate.TperiSSB.gpsNanoSeconds, \
		 thisBinaryTemplate.Eccentricity, \
		 thisBinaryTemplate.ArgPeri,r0,r1,r2);
       }
     }
#endif
#ifdef FILE_FTXT   
   if (!uvar_binary) /* BINARY-MOD - Need to select isolated or binary output */
     {
       /*  PSD.txt file contains freq, PSD, noise floor   */
       for (i=0;i<nbins;i++){ 
	 REAL8 freq;
	 REAL8 r0,r1;
	 freq=uvar_Freq + (REAL8)i*GV.dFreq;
	 r0=F1->data[i];
	 r1=FloorF1->data[i];
	 fprintf(outfile,"%20.17f %E %E\n",freq,r0,r1);
       }
     }
   if (uvar_binary) 
     {
       /*  PSD.txt file contains freq, PSD, noise floor   */
       for (i=0;i<nbins;i++){ 
	 REAL8 freq;
	 REAL8 r0,r1;
	 freq=uvar_Freq + (REAL8)i*GV.dFreq;
	 r0=F1->data[i];
	 r1=FloorF1->data[i];
	 fprintf(outfile,"%20.17f %12.12f %12.12f %d %d %6.12f %6.12f %E %E\n",freq, \
		 thisBinaryTemplate.ProjSMaxis, \
		 thisBinaryTemplate.Period, \
		 thisBinaryTemplate.TperiSSB.gpsSeconds, \
		 thisBinaryTemplate.TperiSSB.gpsNanoSeconds, \
		 thisBinaryTemplate.Eccentricity, \
		 thisBinaryTemplate.ArgPeri,r0,r1);
       }
     }
#endif   

#ifdef FILE_FTXT
   fclose(outfile);
#endif
#ifdef FILE_FLINES  
   fclose(outfile1);
#endif   

   LALFree(outliers->ratio);
   LALFree(outliers->outlierIndexes);
   LALFree(outliers);
   LALFree(outliersParams);
   LALFree(outliersInput);
   LALDDestroyVector(status,&F1);
   LALDDestroyVector(status,&FloorF1);
   LALFree(SpClParams);
   LALFree(clustersInput);

   return 0;


}

/*******************************************************************************/
/*******************************************************************************/

INT4 NormaliseSFTDataRngMdn(LALStatus *status)
{
#ifdef FILE_SPRNG  
  FILE *outfile;
#endif
  INT4 i,j,m,lpc,il;                         /* loop indices */
  INT4 Ntot,nbins=GV.ifmax-GV.ifmin+1;   /* Number of points in SFT's */
  REAL8Vector *Sp=NULL, *RngMdnSp=NULL;   /* |SFT|^2 and its rngmdn  */
  REAL8 B;                          /* SFT Bandwidth */
  REAL8 deltaT,norm,*N, *Sp1;
  INT2 windowSize=uvar_windowsize;                  /* Running Median Window Size*/
  REAL4 xre,xim,xreNorm,ximNorm;


  /* The running median windowSize in this routine determines 
     the sample bias which, instead of log(2.0), must be 
     multiplied by F statistics.
  */

  if (uvar_SignalOnly != 1)
    LALRngMedBias (status, &medianbias, windowSize);


  LALDCreateVector(status, &Sp, (UINT4)nbins);
  LALDCreateVector(status, &RngMdnSp, (UINT4)nbins);

  nbins=(INT2)nbins;

  if(!(N= (REAL8 *) LALCalloc(nbins,sizeof(REAL8)))){ 
    printf("Memory allocation failure of N");
    return 0;
  }
   if(!(Sp1= (REAL8 *) LALCalloc(nbins,sizeof(REAL8)))){ 
    printf("Memory allocation failure of Sp1");
    return 0;
  }

   if( nbins < windowSize ) {
     fprintf( stderr, "The frequency band has too small bins (= %d) compared to window size (= %d) used in EstimateFloor().\n", nbins,windowSize );
     exit(1);
   }
   

  /* loop over each SFTs */
  for (i=0;i<GV.SFTno;i++)         
    {
      /* Set to zero the values */
      for (j=0;j<nbins;j++){
	RngMdnSp->data[j] = 0.0;
	Sp->data[j]       = 0.0;
      }
      
      /* loop over SFT data to estimate noise */
      for (j=0;j<nbins;j++){
	xre=SFTData[i]->fft->data->data[j].re;
	xim=SFTData[i]->fft->data->data[j].im;
	Sp->data[j]=(REAL8)(xre*xre+xim*xim);
      }

      /* Compute running median */
      EstimateFloor(Sp, windowSize, RngMdnSp);

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
      for (lpc=0;lpc<nbins;lpc++){
	N[lpc]=1.0/sqrt(2.0*RngMdnSp->data[lpc]);
      }
      
      if(uvar_SignalOnly == 1){
	B=(1.0*GV.nsamples)/(1.0*GV.tsft);
	deltaT=1.0/(2.0*B);
	norm=deltaT/sqrt(GV.tsft);
	for (lpc=0;lpc<nbins;lpc++){
	  N[lpc]=norm;
	}
      }
      
      /*  loop over SFT data to normalise it (with N) */
      /*  also compute Sp1, average normalized PSD */
      /*  and the sum of the PSD in the band, SpSum */
      for (j=0;j<nbins;j++){
	xre=SFTData[i]->fft->data->data[j].re;
	xim=SFTData[i]->fft->data->data[j].im;
	xreNorm=N[j]*xre; 
	ximNorm=N[j]*xim; 
	SFTData[i]->fft->data->data[j].re = xreNorm;    
	SFTData[i]->fft->data->data[j].im = ximNorm;
	Sp1[j]=Sp1[j]+xreNorm*xreNorm+ximNorm*ximNorm;
      }
      
    } /* end loop over SFTs*/

#ifdef FILE_SPRNG  
  if(!(outfile=fopen("SpRng.txt","w"))){ 
    printf("Cannot open output file"); 
    return 1;
  } 


  for (j=0;j<nbins;j++){
    Sp1[j]=2.0*Sp1[j]/(1.0*GV.SFTno);
    fprintf(outfile,"%f %E \n",(GV.ifmin+j)/GV.tsft,Sp1[j]); 
  }
  
  fclose(outfile);
#endif

  
  LALFree(N);
  LALFree(Sp1);
  LALDDestroyVector(status, &RngMdnSp);
  LALDDestroyVector(status, &Sp);
  
  return 0;
}

/*******************************************************************************/

#if USE_BOINC
void use_boinc_filename0(char *orig_name ) {
  char resolved_name[512];
  if (boinc_resolve_filename(orig_name, resolved_name, sizeof(resolved_name))) {
    fprintf(stderr, "Can't resolve file %s\n", orig_name);
    boinc_finish(2);
  }
  strcpy(orig_name, resolved_name);
  return;
}

void use_boinc_filename1(char **orig_name ) {
  char resolved_name[512];
  if (boinc_resolve_filename(*orig_name, resolved_name, sizeof(resolved_name))) {
    fprintf(stderr, "Can't resolve file %s\n", *orig_name);
    boinc_finish(2);
  }
  *orig_name = calloc(strlen(resolved_name)+1,1);
  strcpy(*orig_name, resolved_name);
  return;
}

#endif /*USE_BOINC*/

/**************************************************************************************/

int ReadBinaryTemplateBank(void)
{

  FILE *BTBfp;
  char filename[256];
  UINT4 i;
  REAL8 temp1,temp2;

  
  strcpy(filename,uvar_binarytemplatefile);
  /*  something that reads the binary input file to memory */
  if(!(BTBfp=fopen(filename,"r"))){ 
    printf("Cannot open BinaryTemplate file %s\n",filename); 
    return 1;
  } 
 
  /* allocate memory for header and template info */
  BinaryBank=(BinaryTemplateBank *)LALMalloc(1*sizeof(BinaryTemplateBank));
  /*BinaryBank->BMFheader=(BinaryMeshFileHeader *)(sizeof(BinaryMeshFileHeader));*/
 
  /* read the header information in */
  if (ReadMeshFileHeader(BTBfp,&(BinaryBank->BMFheader))) return 1;

 
  /* Do initial validation of header information */
  if (BinaryBank->BMFheader.f_max<0.0) {
    printf("In BinaryTemplate file %s : header value of fmax < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.tspan<0.0) {
    printf("In BinaryTemplate file %s : header value of Tspan < 0.0\n",filename);
    return 1;
  }
    if ((BinaryBank->BMFheader.tstart.gpsSeconds!=GV.Ti)||(BinaryBank->BMFheader.tstart.gpsNanoSeconds!=0)) {
      printf("In BinaryTemplate file %s : header value of Tstart != Tstart of data\n",filename);
      return 1;
    }
  if (BinaryBank->BMFheader.Nfilters<1) {
    printf("In BinaryTemplate file %s : header value of NFilters < 1\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.mismatch<0.0) {
    printf("In BinaryTemplate file %s : header value of Mismatch < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.sma_MIN<0.0) {
    printf("In BinaryTemplate file %s : header value of Minimum Projected semi-major axis < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.sma_MAX<0.0) {
    printf("In BinaryTemplate file %s : header value of Maximum Projected semi-major axis < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.sma_MAX<BinaryBank->BMFheader.sma_MIN) {
    printf("In BinaryTemplate file %s : header value of Maximum Projected semi-major axis < Minimum\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.tperi_MIN.gpsSeconds<0) {
    printf("In BinaryTemplate file %s : header value of Minimum SSB time of periapse passage < 0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.tperi_MAX.gpsSeconds<0) {
    printf("In BinaryTemplate file %s : header value of Maximum SSB time of periapse passage < 0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.tperi_MIN.gpsNanoSeconds<0) {
    printf("In BinaryTemplate file %s : header value of Minimum SSB time of periapse passage (nanoseconds) < 0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.tperi_MAX.gpsNanoSeconds<0) {
    printf("In BinaryTemplate file %s : header value of Maximum SSB time of periapse passage (nanoseconds) < 0\n",filename);
    return 1;
  }
  temp1=BinaryBank->BMFheader.tperi_MIN.gpsSeconds+1e-9*BinaryBank->BMFheader.tperi_MIN.gpsNanoSeconds;
  temp2=BinaryBank->BMFheader.tperi_MAX.gpsSeconds+1e-9*BinaryBank->BMFheader.tperi_MAX.gpsNanoSeconds;
  if (temp2<temp1) {
    printf("In BinaryTemplate file %s : header value of Maximum Projected semi-major axis < Minimum\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.ecc_MIN<0.0) {
    printf("In BinaryTemplate file %s : header value of Minimum eccentricity < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.ecc_MAX<0.0) {
    printf("In BinaryTemplate file %s : header value of Maximum eccentricity < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.ecc_MAX<BinaryBank->BMFheader.ecc_MIN) {
    printf("In BinaryTemplate file %s : header value of Maximum eccentricity < Minimum\n",filename);
    return 1;
  }
  if ((BinaryBank->BMFheader.argp_MIN<0.0)||(BinaryBank->BMFheader.argp_MIN>LAL_TWOPI)) {
    printf("In BinaryTemplate file %s : header value of Minimum argument of periapse not in range (0 - 2*PI)\n",filename);
    return 1;
  }
  if ((BinaryBank->BMFheader.argp_MAX<0.0)||(BinaryBank->BMFheader.argp_MAX>LAL_TWOPI)) {
    printf("In BinaryTemplate file %s : header value of Maximum argument of periapse not in range (0 - 2*PI)\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.argp_MAX<BinaryBank->BMFheader.argp_MIN) {
    printf("In BinaryTemplate file %s : header value of Maximum argument of periapse < Minimum\n",filename);
    return 1;
  }
  if ((BinaryBank->BMFheader.argp_MIN<0.0)||(BinaryBank->BMFheader.argp_MIN>LAL_TWOPI)) {
    printf("In BinaryTemplate file %s : header value of Minimum argument of periapse not in range (0 - 2*PI)\n",filename);
    return 1;
  }
  if ((BinaryBank->BMFheader.argp_MAX<0.0)||(BinaryBank->BMFheader.argp_MAX>LAL_TWOPI)) {
    printf("In BinaryTemplate file %s : header value of Maximum argument of periapse not in range (0 - 2*PI)\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.argp_MAX<BinaryBank->BMFheader.argp_MIN) {
    printf("In BinaryTemplate file %s : header value of Maximum argument of periapse < Minimum\n",filename);
    return 1;
  }  
  if (BinaryBank->BMFheader.period_MIN<0.0) {
    printf("In BinaryTemplate file %s : header value of Minimum period < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.period_MAX<0.0) {
    printf("In BinaryTemplate file %s : header value of Maximum period < 0.0\n",filename);
    return 1;
  }
  if (BinaryBank->BMFheader.period_MAX<BinaryBank->BMFheader.period_MIN) {
    printf("In BinaryTemplate file %s : header value of Maximum period < Minimum\n",filename);
    return 1;
  }  
  
  /* allocate memory for templates */
  (BinaryBank->BTB)=(BinaryTemplate *)LALMalloc(BinaryBank->BMFheader.Nfilters*sizeof(BinaryTemplate));
  
  /* Now read in all templates into memory */
  i=0;
  while (i<BinaryBank->BMFheader.Nfilters) {
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

