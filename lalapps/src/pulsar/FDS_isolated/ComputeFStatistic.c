/*********************************************************************************
 ** \file ComputeFStatistic.c
 * Calculate the F-statistic for a given parameter-space of pulsar GW signals.
 * Implements the so-called "F-statistic" as introduced in JKS98.
 *                                                                          
 * \author		   Y. Ioth, M.A. Papa, X. Siemens, R.Prix, B. Machenschalk
 *                                                                          
 *                 Albert Einstein Institute/UWM - started September 2002   
 *********************************************************************************/
#include "config.h"

/* System includes */
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_GLOB_H
#include <glob.h>
#endif

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/RngMedBias.h>
#include <lal/LALDemod.h>
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>

#include <lalapps.h>

/* local includes */
#include "ComputeFStatistic.h"
#include "clusters.h"
#include "DopplerScan.h"

RCSID( "$Id$");

/*----------------------------------------------------------------------*/
/* conditional compilation-switches */

/*
#define NEARESTGRIDPOINTS_ON
#define DEBG_FAFB                
#define DEBG_ESTSIGPAR
#define DEBG_SGV 
*/

/* If FILE_FILENAME is defined, then print the corresponding file  */

BOOLEAN FILE_FSTATS = 1;
/*
#define FILE_FMAX
#define FILE_FLINES
#define FILE_FTXT 
#define FILE_SPRNG
#define FILE_AMCOEFFS
*/

/*----------------------------------------------------------------------*/
/* USE_BOINC should be set to 1 to be run under BOINC */
/* NO_BOINC_GRAPHICS should be unset or 0 to have BOINC application graphics */

#ifndef USE_BOINC
#define USE_BOINC 0
#endif

#if USE_BOINC

/* this includes patches for chdir() and sleep() */
#ifdef _WIN32
#include "win_lib.h"
#endif

#ifdef NO_BOINC_GRAPHICS
#define BOINC_GRAPHICS 0
#endif
#ifndef BOINC_GRAPHICS
#define BOINC_GRAPHICS 0
#endif

/* compress earth and sun files and output Fstats file, using zip.
   This is BOINC-only code */
#ifndef BOINC_COMPRESS
#define BOINC_COMPRESS 0
#endif

/* Boinc diag constants */
#define BOINC_DIAG_DUMPCALLSTACKENABLED     0x00000001L
#define BOINC_DIAG_HEAPCHECKENABLED         0x00000002L
#define BOINC_DIAG_MEMORYLEAKCHECKENABLED   0x00000004L
#define BOINC_DIAG_ARCHIVESTDERR            0x00000008L
#define BOINC_DIAG_ARCHIVESTDOUT            0x00000010L
#define BOINC_DIAG_REDIRECTSTDERR           0x00000020L
#define BOINC_DIAG_REDIRECTSTDOUT           0x00000040L
#define BOINC_DIAG_REDIRECTSTDERROVERWRITE  0x00000080L
#define BOINC_DIAG_REDIRECTSTDOUTOVERWRITE  0x00000100L
#define BOINC_DIAG_TRACETOSTDERR            0x00000200L
#define BOINC_DIAG_TRACETOSTDOUT            0x00000400L

#include <signal.h>

#ifndef _WIN32
#include <pthread.h>
#endif

#define USE_BOINC_DEBUG 0
/* for getpid() */
#include <sys/types.h>
#include <unistd.h>

#include "boinc_api.h"
#include "filesys.h"
#include "diagnostics.h"
#if BOINC_COMPRESS
#include "boinc_zip.h"
#endif

#ifdef HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#if BOINC_GRAPHICS
#include "graphics_api.h"
#include "graphics_lib.h"
#endif

#define fopen boinc_fopen

char *fstatbuff = NULL;
int boincmain(int argc, char *argv[]);
void worker();

/* hooks for communication with the graphics thread */
void (*set_search_pos_hook)(float,float) = NULL;
int (*boinc_init_graphics_hook)(void (*worker)()) = NULL;
double *fraction_done_hook = NULL;

void use_boinc_filename1(char** orig_name);
void use_boinc_filename0(char* orig_name);

#ifdef __cplusplus
extern "C" {
#endif

#if (BOINC_GRAPHICS == 1)
extern double fraction_done;
/* FIXME: include proper header for this! */
extern void set_search_pos(float RAdeg, float DEdeg);
extern int boinc_init_graphics(void (*worker)());
#endif

void sighandler(int sig);

/* polka prototype */
#ifdef RUN_POLKA
int polka(int argc,char *argv[]);
#endif

#ifdef __cplusplus
}
#endif

#endif /* USE_BOINC */

/*----------------------------------------------------------------------*/
/* Error-codes */

#define COMPUTEFSTATC_ENULL 		1
#define COMPUTEFSTATC_ESYS     		2
#define COMPUTEFSTATC_EINPUT   		3
#define COMPUTEFSTATC_EMEM   		4
#define COMPUTEFSTATC_ECHECKPOINT	5
#define COMPUTEFSTATC_ECLUSTER		6

#define COMPUTEFSTATC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATC_MSGESYS		"System call failed (probably file IO)"
#define COMPUTEFSTATC_MSGEINPUT   	"Invalid input"
#define COMPUTEFSTATC_MSGEMEM   	"Out of memory. Bad."
#define COMPUTEFSTATC_MSGECHECKPOINT	"Illegal checkpoint"
#define COMPUTEFSTATC_MSGECLUSTER	"Unspecified error in cluster-related routine"

/*----------------------------------------------------------------------*/
/* Exit values */
#define COMPUTEFSTAT_EXIT_OK              0  /* normal exit */
#define COMPUTEFSTAT_EXIT_USAGE           7  /* user requested help */
#define COMPUTEFSTAT_EXIT_READSFTFAIL     8  /* ReadSFT failed */
#define COMPUTEFSTAT_EXIT_OPENFMAX        9  /* error opening Fmax file */
#define COMPUTEFSTAT_EXIT_OPENFSTAT      10  /* error opening FStats file */
#define COMPUTEFSTAT_EXIT_OPENFSTAT2     11  /* error opening FStats file for append (checkpointing) */
#define COMPUTEFSTAT_EXIT_WRITEFSTAT     12  /* error opening FStats file */
#define COMPUTEFSTAT_EXIT_WRITEFAFB      13  /* writeFaFb failed */
#define COMPUTEFSTAT_EXIT_ESTSIGPAR      14  /* EstimateSignalParameters failed */
#define COMPUTEFSTAT_EXIT_NOMEM          15  /* out of memory */
#define COMPUTEFSTAT_EXIT_CANTZIP        16  /* unable to zip Fstats file */
#define COMPUTEFSTAT_EXIT_CANTUNZIP      17  /* unable to zip Fstats file */
#define COMPUTEFSTAT_EXIT_CANTRENAME     18  /* unable to zip Fstats file */
#define COMPUTEFSTAT_EXIT_NOPOLKADEL     19  /* no // found in command line */
#define COMPUTEFSTAT_EXIT_LALCALLERROR  100  /* this is added to the LAL status to get BOINC exit value */

/*----------------------------------------------------------------------
 * User-variables: can be set from config-file or command-line */

INT4 uvar_Dterms;
CHAR* uvar_IFO;
BOOLEAN uvar_SignalOnly;
BOOLEAN uvar_EstimSigParam;
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
CHAR *uvar_ephemDir;
CHAR *uvar_ephemYear;
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
BOOLEAN uvar_doCheckpointing;
BOOLEAN uvar_expLALDemod;
REAL8 uvar_startTime;
REAL8 uvar_endTime;
#if BOINC_COMPRESS
BOOLEAN uvar_useCompression;
#endif
/*----------------------------------------------------------------------*/
/* some other global variables */

FFT **SFTData=NULL; 		/**< SFT Data for LALDmod */
DemodPar *DemodParams  = NULL;	/**< Demodulation parameters for LALDemod */
LIGOTimeGPS *timestamps=NULL;	/**< Time stamps from SFT data */
LALFstat Fstat;			/**< output from LALDemod(): F-statistic and amplitudes Fa and Fb */
AMCoeffs amc;			/**< amplitude-modulation coefficients (and derived quantities) */
REAL8 Alpha,Delta;		/**< sky-position currently searched (equatorial coords, radians) */
Clusters HFLines;		/**< stores information about outliers/clusters in F-statistic */
Clusters HPLines;		/**< stores information about outliers/clusters in SFT-power spectrum */
Clusters *highSpLines=&HPLines, *highFLines=&HFLines;

REAL8 medianbias=1.0;		/**< bias in running-median depending on window-size (set in NormaliseSFTDataRngMdn()) */

FILE *fp_mergedSFT;		/**< input-file containing merged SFTs */
FILE *fpmax;			/**< output-file: maximum of F-statistic over frequency-range */
FILE *fpstat=NULL;		/**< output-file: F-statistic candidates and cluster-information */

ConfigVariables GV;		/**< global container for various derived configuration settings */
int reverse_endian=-1;          /**< endian order of SFT data.  -1: unknown, 0: native, 1: reversed */
CHAR Fstatsfilename[256]; 		/* Fstats file name*/
CHAR ckp_fname[260];			/* filename of checkpoint-file, global for polka */
CHAR *Outputfilename; 		/* Name of output file, either Fstats- or Polka file name*/
INT4 cfsRunNo = 0; 		    /* the CFS run no if ran multiple times 0=run only once, 1=first run, 2=second run */


/*----------------------------------------------------------------------*/
/* local prototypes */
#ifdef __cplusplus
extern "C" {
#endif
  int main(int argc,char *argv[]);
  void initUserVars (LALStatus *stat);
  INT4 ReadSFTData (void);
  void InitFStat (LALStatus *status, ConfigVariables *cfg);
  void CreateDemodParams (LALStatus *status);
  void CreateNautilusDetector (LALStatus *status, LALDetector *Detector);
  void Freemem (LALStatus *status);
  void EstimateFLines(LALStatus *status);
  void NormaliseSFTDataRngMdn (LALStatus *status);
  INT4 EstimateSignalParameters(INT4 * maxIndex);
  int writeFLines(INT4 *maxIndex, int *bytes_written, UINT4 *checksum);
  INT4 PrintTopValues(REAL8 TwoFthr, INT4 ReturnMaxN);
  int compare(const void *ip, const void *jp);
  INT4 writeFaFb(INT4 *maxIndex);
  void InitDopplerScanOnRefinedGrid ( LALStatus *status, DopplerScanState *theScan, DopplerScanInit *scanInit);

  void WriteFStatLog (LALStatus *stat, CHAR *argv[]);
  static void swap2(char *location);
  static void swap4(char *location);
  static void swap8(char *location);
  void swapheader(struct headertag *thisheader);
  void getCheckpointCounters(LALStatus *stat, UINT4 *loopcounter, UINT4 *checksum, long *bytecounter, const CHAR *fstat_fname, const CHAR *ckpfn);
  void TestLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params);
  void OrigLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params);
#ifdef FILE_AMCOEFFS
  void PrintAMCoeffs (REAL8 Alpha, REAL8 Delta, AMCoeffs* amc);
#endif
#ifdef __cplusplus
}
#endif


/*----------------------------------------------------------------------*/
/* some local defines */

#define EPHEM_YEARS  "00-04"
#define SFT_BNAME  "SFT"

#ifndef TRUE 
#define TRUE (1==1)
#endif
#ifndef FALSE
#define FALSE (1==0)
#endif

extern int vrbflg;

/*----------------------------------------------------------------------*/
/* empty structs for initialziations */

DopplerScanState emptyScan;

/* initialize LAL-status.  
 * Note: for this to be thread safe, you MUST not use *this* status structure 
 * for LAL functions that are called from the graphics thread.  Of course there 
 * is probably no need to call LAL functions from the graphics thread, so this 
 * should not be an issue.
 *
 * NOTE2: This is the global *head* of main's status structure.
 *        The ONLY reason for this to be global is for the signal-handler to 
 *        be able to access it. 
 *        You should NEVER use this struct directly ANYWHERE in the code !!
 */
LALStatus global_status;

/*----------------------------------------------------------------------*/
/* CODE starts here */
/*----------------------------------------------------------------------*/


#if USE_BOINC
int BOINC_ERR_EXIT(LALStatus  *stat, const char *func, const char *file, const int line, volatile const char *id) {
  if (stat->statusCode) {
    fprintf(stderr,
	    "Level 0: %s\n"
	    "\tFunction call `%s' failed.\n"
	    "\tfile %s, line %d\n",
	    id, func, file, line );
    REPORTSTATUS(stat);
    boinc_finish( COMPUTEFSTAT_EXIT_LALCALLERROR+stat->statusCode );
  }
  /* should this call boinc_finish too?? */
  return 0;
}
#endif

/** 
 * MAIN function of ComputeFStatistic code.
 * Calculate the F-statistic over a given portion of the parameter-space
 * and write a list of 'candidates' into a file(default: 'Fstats').
 */
 
#if USE_BOINC
int boincmain(int argc, char *argv[])
#else
int main(int argc,char *argv[]) 
#endif
{
  LALStatus *stat = &global_status;

  INT4 *maxIndex=NULL; 			/*  array that contains indexes of maximum of each cluster */
  INT4 spdwn;				/* counter over spindown-params */
  DopplerScanInit scanInit;		/* init-structure for DopperScanner */
  DopplerScanState thisScan = emptyScan; /* current state of the Doppler-scan */
  DopplerPosition dopplerpos;		/* current search-parameters */
  SkyPosition thisPoint;
  LIGOTimeGPS t0, t1;
  REAL8 duration;
  FILE *fpOut = NULL;

  UINT4 loopcounter;		/* Checkpoint-counters for restarting checkpointed search */
  long fstat_bytecounter;
  UINT4 checksum=0;             /* Checksum of fstats file contents */

  /* set LAL error-handler */
#if USE_BOINC
  CHAR resfname[256];	/* buffer for boinc-resolving config-file name */

  lal_errhandler = BOINC_ERR_EXIT;
#else
  lal_errhandler = LAL_ERR_EXIT;
#endif

  lalDebugLevel = 0 ;  
  vrbflg = 1;	/* verbose error-messages */

  /* register all user-variable */
  LAL_CALL (LALGetDebugLevel(stat, argc, argv, 'v'), stat);
  LAL_CALL (initUserVars(stat), stat); 	

#if USE_BOINC
  /* handle config file request with boinc_resolve_filename */
  /* NOTE: @configfile must be at the beginning of the command line! */
  if ( argv[1] && argv[1][0] == '@') 
    {
      resfname[0] = '@';
      if (boinc_resolve_filename(argv[1]+1,resfname+1,sizeof(resfname)))
	fprintf(stderr,"WARNING: Can't boinc-resolve config file \"%s\"\n", argv[1]+1);
      else
	{
	  /* hack the command-line: replace config-file by boinc-resolved path */
	  argv[1] = resfname;
	  /* note: we don't free the previous argv[1] in order to avoid possible problems..*/
	}
    } /* if config-file given as first argument */
#endif

  LAL_CALL (LALUserVarReadAllInput(stat,argc,argv),stat);	


  if (uvar_help)	/* if help was requested, we're done here */
    return COMPUTEFSTAT_EXIT_USAGE;

  /* This is dangerous for BOINC since it calls system() and makes
     assumptions that might not be true */
#if !USE_BOINC
  /* keep a log-file recording all relevant parameters of this search-run */
  LAL_CALL (WriteFStatLog (stat, argv), stat);
#endif /* !USE_BOINC */

  /* main initialization of the code: */
  LAL_CALL ( InitFStat(stat, &GV), stat);

  /* read in SFT-data */
  if (ReadSFTData()) return COMPUTEFSTAT_EXIT_READSFTFAIL;

  /* normalize SFTs by running median */
  LAL_CALL (NormaliseSFTDataRngMdn(stat), stat);

#ifdef FILE_FMAX
  {
    CHAR Fmaxfilename[256]; 		/* Fmax file name*/

    /*   open file */
    strcpy(Fmaxfilename,"Fmax");
    if (uvar_outputLabel)
      strcat(Fmaxfilename,uvar_outputLabel);
#if USE_BOINC
    use_boinc_filename0(Fmaxfilename);
#endif /* USE_BOINC */
    if (!(fpmax=fopen(Fmaxfilename,"wb"))){
      fprintf(stderr,"in Main: unable to open Fmax file %s\n", Fmaxfilename);
      return COMPUTEFSTAT_EXIT_OPENFMAX;
    }
  }
#endif
  
  /*----------------------------------------------------------------------
   * prepare initialization of the DopplerScanner to step through paramter space 
   * In input-structure for this initialization is DopplerScanInit 
   */
  if (lalDebugLevel) LALPrintError ("\nSetting up template grid ...");
  scanInit.dAlpha = uvar_dAlpha;
  scanInit.dDelta = uvar_dDelta;
  scanInit.gridType = (DopplerGridType) uvar_gridType;
  scanInit.metricType = (LALPulsarMetricType) uvar_metricType;
  scanInit.metricMismatch = uvar_metricMismatch;
  t0 = SFTData[0]->fft->epoch;
  t1 = SFTData[GV.SFTno-1]->fft->epoch;
  scanInit.obsBegin = t0;
  LAL_CALL ( LALDeltaFloatGPS( stat, &duration, &t1, &t0), stat);
  scanInit.obsDuration = duration + GV.tsft;
  scanInit.fmax  = uvar_Freq;
  scanInit.fmax += uvar_FreqBand;
  scanInit.Detector = &GV.Detector;
  scanInit.ephemeris = GV.edat;		/* used by Ephemeris-based metric */
  scanInit.skyRegion = GV.skyRegion;
  scanInit.skyGridFile = uvar_skyGridFile;	/* if applicable */

  if ( uvar_searchNeighbors ) {
    LAL_CALL ( InitDopplerScanOnRefinedGrid( stat, &thisScan, &scanInit ), stat );
  } else {
    LAL_CALL ( InitDopplerScan( stat, &thisScan, &scanInit), stat); 
  }
  if (lalDebugLevel) LALPrintError ("done.\n");

  /* should we write the sky-grid to disk? */
  if ( uvar_outputSkyGrid ) 
    {
      LALPrintError ("\nNow writing sky-grid into file '%s' ...", uvar_outputSkyGrid);
      LAL_CALL (writeSkyGridFile( stat, thisScan.grid, uvar_outputSkyGrid, &scanInit), stat);
      LALPrintError (" done.\n\n");
    }

  /* if a complete output of the F-statistic file was requested,
   * we open and prepare the output-file here */
  if (uvar_outputFstat) 
    {
      if ( (fpOut = fopen (uvar_outputFstat, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar_outputFstat);
	  return COMPUTEFSTAT_EXIT_OPENFSTAT;
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

  /* prepare main output-file "Fstats": get actual filename */
  strcpy(Fstatsfilename,"Fstats");
  if ( uvar_outputLabel )
    strcat(Fstatsfilename,uvar_outputLabel);
  /* prepare checkpointing file */
  strcpy(ckp_fname, Fstatsfilename);
  strcat(ckp_fname, ".ckp");

#if USE_BOINC
  /* only boinc_resolve the filename if we run CFS once */
  if (cfsRunNo == 0)
	  use_boinc_filename0(Fstatsfilename);
  /* use_boinc_filename0(ckp_fname); */
#endif /* USE_BOINC */

  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
   */

  /* we allow for checkpointing here, so we might not actually start this
   * loop at 0, but at some previously stored loop-counter.
   * In addition we retrieve the number of bytes that SHOULD have
   * been written to fstats-file so far, in order to catch file-corruptions.
   */

  loopcounter = 0;
  fstat_bytecounter = 0;
  checksum =0;
  
  /* Checkpointed information is retrieved from checkpoint-file (if found) */
  if (uvar_doCheckpointing)
    LAL_CALL (getCheckpointCounters( stat, &loopcounter, &checksum, &fstat_bytecounter, Fstatsfilename, ckp_fname ), stat); 


  /* allow for checkpointing: 
   * open fstats file for writing or appending, depending on loopcounter. 
   */
  if ( FILE_FSTATS && ( (fpstat=fopen( Fstatsfilename, fstat_bytecounter>0 ? "rb+" : "wb")) == NULL) )
    {
      fprintf(stderr,"in Main: unable to open Fstats file\n");
      return COMPUTEFSTAT_EXIT_OPENFSTAT2;
    }

#if USE_BOINC
  { /* set a buffer large enough that no output is written to disk
       unless we fflush().  Policy is fully buffered. Note: the man
       page says "The setvbuf function may only be used after opening
       a stream and BEFORE ANY OTHER OPERATIONS HAVE BEEN PERFORMED ON
       IT."
    */
    const int bsize=2*1024*1024;
    if ((fstatbuff=(char *)LALCalloc(bsize, 1)))
      setvbuf(fpstat, fstatbuff, _IOFBF, bsize);
    else {
      LALPrintError ("Failed to allocate memory for Fstats file buffering. Exiting.\n");
      return COMPUTEFSTAT_EXIT_NOMEM;
    }
  }
#endif

  /* if we checkpointed in the loop we need to "spool forward" accordingly in our sky-position list */
  if ( loopcounter > 0 ) {
    UINT4 i;
    for (i=0; i < loopcounter; i++) {
      LAL_CALL (NextDopplerPos( stat, &dopplerpos, &thisScan ), stat);
      if (thisScan.state == STATE_FINISHED) {
	LALPrintError ("Error: checkpointed loopcounter already at the end of main-loop\n");
	return COMPUTEFSTATC_ECHECKPOINT; 
      }
    }
    /* seek to right point of fstats file (truncate what's left over) */
    if ( 0 != fseek( fpstat, fstat_bytecounter, SEEK_SET) ) {	/* something gone wrong seeking .. */
      if (lalDebugLevel) LALPrintError ("broken fstats-file.\nStarting main-loop from beginning.\n");
      return COMPUTEFSTATC_ECHECKPOINT;;
    }
  } /* if loopcounter > 0 */
    
  while (1)
    {
      /* flush fstats-file and write checkpoint-file */
      if (uvar_doCheckpointing && fpstat)
	{
#if USE_BOINC
	  if (boinc_time_to_checkpoint())
	    {
#endif
	      FILE *fp;
	      fflush (fpstat);
	      if ( (fp = fopen(ckp_fname, "wb")) == NULL) {
		LALPrintError ("Failed to open checkpoint-file for writing. Exiting.\n");
		return COMPUTEFSTATC_ECHECKPOINT;
	      }
	      if ( fprintf (fp, "%" LAL_UINT4_FORMAT " %" LAL_UINT4_FORMAT " %ld\nDONE\n", loopcounter, checksum, fstat_bytecounter) < 0) {
		LALPrintError ("Error writing to checkpoint-file. Exiting.\n");
		return COMPUTEFSTATC_ECHECKPOINT;
	      }
	      fclose (fp);
#if USE_BOINC
	      boinc_checkpoint_completed();
	    } /* if boinc_time_to_checkpoint() */
#endif
	} /* if doCheckpointing && fpstat */
      

      /* Show some progress */
#if USE_BOINC
      {
	double local_fraction_done;
	if (cfsRunNo == 1)
	  local_fraction_done=(((double)loopcounter)/((double)thisScan.numGridPoints))/2;
	else if (cfsRunNo == 2)
	  local_fraction_done=(((double)loopcounter)/((double)thisScan.numGridPoints))/2+0.5;
	else
	  local_fraction_done=((double)loopcounter)/((double)thisScan.numGridPoints);
	if (local_fraction_done<0.0)
	  local_fraction_done=0.0;
	if (local_fraction_done>1.0)
	  local_fraction_done=1.0;
	boinc_fraction_done(local_fraction_done);
 	    /* pass variable externally to graphics routines */
	    if (fraction_done_hook != NULL)
	      *fraction_done_hook=local_fraction_done;
      }
#endif
      if (lalDebugLevel) LALPrintError ("Search progress: %5.1f%%", 
					(100.0* loopcounter / thisScan.numGridPoints));
      
      LAL_CALL (NextDopplerPos( stat, &dopplerpos, &thisScan ), stat);

      /* Have we scanned all DopplerPositions yet? */
      if (thisScan.state == STATE_FINISHED)
	break;

      LAL_CALL (LALNormalizeSkyPosition(stat, &thisPoint, &(dopplerpos.skypos) ), stat);

      Alpha = thisPoint.longitude;
      Delta = thisPoint.latitude;
#if USE_BOINC
      /* pass current search position, for use with starsphere.C
	 revision 4.6 or greater. Need to convert radians to
	 degrees. */
	  if (set_search_pos_hook != NULL){
		float pAlpha=(180.0*Alpha/LAL_PI), pDelta=(180.0*Delta/LAL_PI);
        set_search_pos_hook(pAlpha,pDelta);
	  }
#endif
	  
      LAL_CALL (CreateDemodParams(stat), stat);
#ifdef FILE_AMCOEFFS
	  PrintAMCoeffs(Alpha, Delta, DemodParams->amcoe);
#endif
	  /* loop over spin params */
      for (spdwn=0; spdwn <= GV.SpinImax; spdwn++)
	{
	  DemodParams->spinDwn[0] = uvar_f1dot + spdwn * uvar_df1dot;

	  if ( uvar_expLALDemod )
	    LAL_CALL ( TestLALDemod(stat, &Fstat, SFTData, DemodParams), stat);
	  else
	    LAL_CALL ( OrigLALDemod(stat, &Fstat, SFTData, DemodParams), stat);

	  /*  This fills-in highFLines that contains the outliers of F*/
	  if (GV.FreqImax > 5) {
	    LAL_CALL (EstimateFLines(stat), stat);
	  }
	  
	  /* if user requested it, we output ALL F-statistic results */
	  if (fpOut) 
	    {
	      INT4 i;
	      for(i=0;i < GV.FreqImax ;i++)
		{
		  fprintf (fpOut, "%20.17f %20.17f %20.17f %20.17f\n", 
			   uvar_Freq + i*GV.dFreq, Alpha, Delta, 2.0*medianbias*Fstat.F[i]);
		}

	    } /* if outputFstat */

	  
	  /*  This fills-in highFLines  */
	  if (highFLines != NULL && highFLines->Nclusters > 0)
	    {
	      int bytesWritten = 0;

	      maxIndex=(INT4 *)LALMalloc(highFLines->Nclusters*sizeof(INT4));
	    
	      /*  for every cluster writes the information about it in file Fstats */
	      if (writeFLines(maxIndex, &bytesWritten, &checksum))
		{
		  fprintf(stderr, "%s: trouble making file Fstats\n", argv[0]);
		  return COMPUTEFSTAT_EXIT_WRITEFSTAT;
		}
	      fstat_bytecounter += (long)bytesWritten;	/* bookkeeping of nominal length of Fstats-file */
	   	  
	      if (uvar_EstimSigParam)
		{
		  if(writeFaFb(maxIndex)) return COMPUTEFSTAT_EXIT_WRITEFAFB;
		  if (EstimateSignalParameters(maxIndex)) return COMPUTEFSTAT_EXIT_ESTSIGPAR;
		} /* if signal-estimation */
	  
	  
	      LALFree(maxIndex);
	    } /* if highFLines found */

	  if (PrintTopValues(/* thresh */ 0.0, /* max returned */ 1))
	    LALPrintError ("%s: trouble making files Fmax and/or Fstats\n", argv[0]);

	  
	  /* Set the number of the clusters detected to 0 at each iteration 
	     of the sky-direction and the spin down */
	  highFLines->Nclusters=0;

	} /* For GV.spinImax */

      loopcounter ++;		/* number of *completed* loops */

    } /*  while SkyPos */
  
  if (uvar_outputFstat && fpOut)
    fclose (fpOut);

  if (lalDebugLevel) LALPrintError ("\nSearch finished.\n");
  
#ifdef FILE_FMAX  
  fclose(fpmax);
#endif

  if (fpstat) {
    /* this is our marker indicating 'finished'.  What appears in the file is:
    %DONE
    on the final line */
    fprintf(fpstat, "%%DONE\n");
    fclose(fpstat);
  }

#ifndef RUN_POLKA
  /* remove checkpoint-file */
  remove (ckp_fname);
#endif

  /* Free DopplerScan-stuff (grid) */
  LAL_CALL ( FreeDopplerScan(stat, &thisScan), stat);
  
  LAL_CALL ( Freemem(stat), stat);
  
  return COMPUTEFSTAT_EXIT_OK;
  
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
  uvar_dFreq = 0;

  uvar_Alpha 	= 0.0;
  uvar_Delta 	= 0.0;
  uvar_AlphaBand = 0;
  uvar_DeltaBand = 0;
  uvar_dAlpha 	= 0.001;
  uvar_dDelta 	= 0.001;
  uvar_skyRegion = NULL;

  uvar_ephemYear = (CHAR*)LALCalloc (1, strlen(EPHEM_YEARS)+1);


#if USE_BOINC
  strcpy (uvar_ephemYear, "");		  /* the default year-string under BOINC is empty! */
#else
  strcpy (uvar_ephemYear, EPHEM_YEARS);
#endif

  uvar_BaseName	= (CHAR*)LALCalloc (1, strlen(SFT_BNAME)+1);
  strcpy (uvar_BaseName, SFT_BNAME);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar_ephemDir = (CHAR*)LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar_ephemDir, DEFAULT_EPHEMDIR);

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

  uvar_workingDir = (CHAR*)LALMalloc(512);
  strcpy(uvar_workingDir, ".");

  uvar_searchNeighbors = FALSE;
  
  /* if user does not set start/end times, use all SFTs */
  uvar_startTime = -1.0e308;
  uvar_endTime   =  1.0e308;

/* some BOINC-defaults differ: checkpointing is ON, and we use experimental LALDemod() */
#if USE_BOINC
  uvar_doCheckpointing = TRUE;
  uvar_expLALDemod = TRUE;
#else
  uvar_doCheckpointing = FALSE;
  uvar_expLALDemod = FALSE;
#endif

#if BOINC_COMPRESS
  uvar_useCompression = TRUE;
#endif

  /* register all our user-variables */
  LALregBOOLUserVar(stat, 	help, 		'h', UVAR_HELP,     "Print this message"); 
  LALregSTRINGUserVar(stat, 	IFO, 		'I', UVAR_REQUIRED, "Detector: GEO(0),LLO(1),LHO(2),NAUTILUS(3),VIRGO(4),TAMA(5),CIT(6)");
  LALregREALUserVar(stat, 	Freq, 		'f', UVAR_REQUIRED, "Starting search frequency in Hz");
  LALregREALUserVar(stat, 	FreqBand, 	'b', UVAR_OPTIONAL, "Search frequency band in Hz");
  LALregREALUserVar(stat, 	dFreq, 		'r', UVAR_OPTIONAL, "Frequency resolution in Hz (default: 1/(2*Tsft*Nsft)");
  LALregREALUserVar(stat, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserVar(stat, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
  LALregREALUserVar(stat, 	AlphaBand, 	'z', UVAR_OPTIONAL, "Band in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(stat, 	DeltaBand, 	'c', UVAR_OPTIONAL, "Band in delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(stat,	skyRegion, 	'R', UVAR_OPTIONAL, "ALTERNATIVE: specify sky-region by polygon");
  LALregINTUserVar(stat, 	gridType,	 0 , UVAR_OPTIONAL, "Template grid: 0=flat, 1=isotropic, 2=metric, 3=file");
  LALregINTUserVar(stat, 	metricType,	'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserVar(stat, 	metricMismatch,	'X', UVAR_OPTIONAL, "Maximal mismatch for metric tiling");
  LALregREALUserVar(stat, 	dAlpha, 	'l', UVAR_OPTIONAL, "Resolution in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(stat, 	dDelta, 	'g', UVAR_OPTIONAL, "Resolution in delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(stat,	skyGridFile,	 0,  UVAR_OPTIONAL, "Load sky-grid from this file.");
  LALregSTRINGUserVar(stat,	outputSkyGrid,	 0,  UVAR_OPTIONAL, "Write sky-grid into this file.");
  LALregREALUserVar(stat, 	f1dot, 		's', UVAR_OPTIONAL, "First spindown parameter f1dot");
  LALregREALUserVar(stat, 	f1dotBand, 	'm', UVAR_OPTIONAL, "Search-band for f1dot");
  LALregREALUserVar(stat, 	df1dot, 	'e', UVAR_OPTIONAL, "Resolution for f1dot (default 1/(2*Tobs*Tsft*Nsft)");
  LALregSTRINGUserVar(stat,	DataDir, 	'D', UVAR_OPTIONAL, "Directory where SFT's are located");
  LALregSTRINGUserVar(stat,	BaseName, 	'i', UVAR_OPTIONAL, "The base name of the input  file you want to read");
  LALregSTRINGUserVar(stat,	mergedSFTFile, 	'B', UVAR_OPTIONAL, "Merged SFT's file to be used"); 
  LALregSTRINGUserVar(stat,	ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(stat,	ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregBOOLUserVar(stat, 	SignalOnly, 	'S', UVAR_OPTIONAL, "Signal only flag");
  LALregBOOLUserVar(stat, 	EstimSigParam, 	'p', UVAR_OPTIONAL, "Do Signal Parameter Estimation");
  LALregREALUserVar(stat, 	Fthreshold,	'F', UVAR_OPTIONAL, "Signal Set the threshold for selection of 2F");
  LALregINTUserVar(stat,	Dterms,		't', UVAR_OPTIONAL, "Number of terms to keep in Dirichlet kernel sum");
  LALregSTRINGUserVar(stat,	outputLabel,	'o', UVAR_OPTIONAL, "Label to be appended to all output file-names");
  LALregSTRINGUserVar(stat,	outputFstat,	 0,  UVAR_OPTIONAL, "Output-file for the F-statistic field over the parameter-space");
  LALregBOOLUserVar(stat,	openDX,	 	 0,  UVAR_OPTIONAL, "Make output-files openDX-readable (adds proper header)");
  LALregSTRINGUserVar(stat,     workingDir,     'w', UVAR_OPTIONAL, "Directory to be made the working directory.");
  LALregBOOLUserVar(stat,	searchNeighbors, 0,  UVAR_OPTIONAL, "Refine skyregion to neighboring points of original center.");
  LALregBOOLUserVar(stat,	doCheckpointing, 0,  UVAR_OPTIONAL, "Do checkpointing and resume for previously checkpointed state.");
  LALregBOOLUserVar(stat,	expLALDemod, 	 0,  UVAR_OPTIONAL, "Use new experimental LALDemod-version.");

  LALregREALUserVar(stat, 	startTime, 	 0,  UVAR_OPTIONAL, "Ignore SFTs with GPS_time <  this value. Default:");
  LALregREALUserVar(stat, 	endTime, 	 0,  UVAR_OPTIONAL, "Ignore SFTs with GPS_time >= this value. Default:");
#if BOINC_COMPRESS
  LALregBOOLUserVar(stat,	useCompression,  0,  UVAR_OPTIONAL, "BOINC: use compression for download/uploading data");
#endif

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
  

#ifdef DEBG_ESTSIGPAR
  REAL8 Ftest;
  REAL8 A2test,A3test,A4test;
#endif


  strcpy(Paramfilename,"ParamMLE");
  if (uvar_outputLabel)
    strcat(Paramfilename,uvar_outputLabel);
  
  if(!(fpMLEParam=fopen(Paramfilename,"wb"))) {
    fprintf(stderr,"Error in EstimateSignalParameters: unable to open the file");
    return 1;
  }

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


      if(ampratio<0.25-error_tol||ampratio>2.0+error_tol) {
	fprintf(stderr,
		"Imaginary Cos[iota]; cannot compute parameters\n"
		"in the EstimateSignalParameters routine\n"
		"in ComputeFStatistic code\n"
		"Now exiting...\n");
	/* 	  break; */
	return 1;
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
	return 1;
      }
      if(fabs(A2-A2test)>fabs(A2)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A2\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A2=%f,A2test=%f\n",
		irec,uvar_Freq+irec*GV.dFreq,A2,A2test);
	fprintf(stderr,"relative error Abs((A2-A2test)/A2)=%lf\n",
		fabs(A2-A2test)/fabs(A2));
	return 1;
      }
      if(fabs(A3-A3test)>fabs(A3)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A3\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A3=%f,A3test=%f\n",
		irec,uvar_Freq+irec*GV.dFreq,A3,A3test);
	fprintf(stderr,"relative error Abs((A3-A3test)/A3)=%lf\n",
		fabs(A3-A3test)/fabs(A3));
	return 1;
      }
      if(fabs(A4-A4test)>fabs(A4)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A4\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A4=%f,A4test=%f\n",
		irec,uvar_Freq+irec*GV.dFreq,A1,A1test);
	fprintf(stderr,"relative error Abs((A4-A4test)/A4)=%lf\n",
		fabs(A4-A4test)/fabs(A4));
	return 1;
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
	return 1;
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

    if((fp=fopen(filename,"wb"))==NULL) {
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

  DemodParams->Dterms=uvar_Dterms;
  DemodParams->ifmin=GV.ifmin;

  DemodParams->returnFaFb = uvar_EstimSigParam;

  /* compute the "sky-constants" A and B */
  TRY ( ComputeSky(status->statusPtr, DemodParams->skyConst, 0, csParams), status);  
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

/*  for every cluster writes the information about it in file Fstats */
/*  precisely it writes: */
/*  fr_max alpha delta N_points_of_cluster mean std max (of 2F) */
int writeFLines(INT4 *maxIndex, int *bytes_written, UINT4 *checksum)
{
  INT4 i,j,j1,j2,k,N;
  REAL8 max,log2,mean,var,std,R,fr;
  INT4 imax;
  int numBytes = 0;
  unsigned char tmpline[1024];
  UINT4 localchecksum=*checksum;

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
    if (fpstat) {
      int l;
      int howmany2=0;
      int howmany=sprintf((char *)tmpline, "%16.12f %10.8f %10.8f    %d %10.5f %10.5f %20.17f\n",
			 fr, Alpha, Delta, N, mean, std, max);
      
      if (howmany <= 0) {
	fprintf(stderr,"writeFLines couldn't print to Fstats-file placeholder!\n");
	return 4;
      }
      
      howmany2=fwrite(tmpline, 1, howmany, fpstat);
      
      if (howmany2!=howmany) {
	fprintf(stderr,"writeFLines couldn't print to Fstats-file!\n");
	return 4;
      }
      
      /* update bytecount */
      numBytes+=howmany;
      
      /* update checksum sum */
      for (l=0; l<howmany; l++)
	localchecksum+=(int)tmpline[l];
      
    } /* if fpstat */
    
  }/*  end i loop over different clusters */
  
  /* return total number of bytes written into fstats-file */
  *bytes_written = numBytes;
  *checksum=localchecksum;

  return 0;

} /* writeFLines() */


/** Reads in data from SFT-files.
 *
 * This function reads in the SFTs from the list of files in \em ConfigVariables GV.filelist 
 * or from merged SFTs in uvar_mergedSFTFile.  If user has specified --startTime or --endTime
 * The read SFT-data is stored in the global array \em SFTData and the timestamps 
 * of the SFTs are stored in the global array \em timestamps (both are allocated here).
 *
 * NOTE: this function is obsolete and should be replaced by the use of the SFT-IO lib in LAL.
 */
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

  for (fileno=0;fileno<GV.SFTno; /* INCREMENT IN LOOP */ ) 
    {
      REAL8 thisSFTtime;

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

      if (reverse_endian)
	swapheader(&header);
      
      if (header.endian!=1.0)
	{
	  fprintf(stderr,"First object in file %s is not (double)1.0!\n",GV.filelist[fileno]);
	  fprintf(stderr,"The file might be corrupted\n\n");
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

      /* Move forward in file */
      offset=(GV.ifmin-header.firstfreqindex)*2*sizeof(REAL4);
      errorcode=fseek(fp,offset,SEEK_CUR);
      if (errorcode) 
	{
	  perror(GV.filelist[fileno]);
	  fprintf(stderr,"Can't get to offset %d in file %s\n",offset,GV.filelist[fileno]);
	  return 5;
	}


      /* determine if THIS SFT is in the range of times of those which
	 we need to use.  If so, read the data into arrays, else
	 ignore it.  For my first CVS commit I will not indent this
	 correctly so that the differences are obvious. A later commit
	 will just clean up the indentation but make no changes to the
	 actual non-whitespace code. */
      thisSFTtime=(REAL8)header.gps_sec+(1.e-9)*(REAL8)header.gps_nsec;
      if (uvar_startTime<=thisSFTtime && thisSFTtime<uvar_endTime) {

      /* Make data structures and put time stamps from file into array */
      timestamps[fileno].gpsSeconds = header.gps_sec;
      timestamps[fileno].gpsNanoSeconds = header.gps_nsec;
      ndeltaf=GV.ifmax-GV.ifmin+1;
      SFTData[fileno]=(FFT *)LALMalloc(sizeof(FFT));
      SFTData[fileno]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
      SFTData[fileno]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
      SFTData[fileno]->fft->data->data=(COMPLEX8 *)LALMalloc(ndeltaf*sizeof(COMPLEX8));

      /* Fill in actual SFT data, and housekeeping */
      errorcode=fread((void*)(SFTData[fileno]->fft->data->data), sizeof(COMPLEX8), ndeltaf, fp);
      if (errorcode!=ndeltaf){
	perror(GV.filelist[fileno]);
	fprintf(stderr, "The SFT data was truncated.  Only read %d not %d complex floats\n", (int)errorcode, ndeltaf);
	return 6;
      }
      /* reverse byte order if needed */
      if (reverse_endian) {
	unsigned int cnt;
	for (cnt=0; cnt<ndeltaf; cnt++) {
	  swap4((char *)&(SFTData[fileno]->fft->data->data[cnt].re));
	  swap4((char *)&(SFTData[fileno]->fft->data->data[cnt].im));
	}
      }
      
      SFTData[fileno]->fft->epoch=timestamps[fileno];
      SFTData[fileno]->fft->f0 = GV.ifmin / GV.tsft;
      SFTData[fileno]->fft->deltaF = 1.0 / GV.tsft;
      SFTData[fileno]->fft->data->length = ndeltaf;
      fileno++;
      }

      if (uvar_mergedSFTFile)
	k+=sizeof(header)+header.nsamples*sizeof(COMPLEX8);
      else
	fclose(fp);     /* close file */
      
    }
  return 0;  

} /* ReadSFTData() */

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

  CHAR command[512];
  FILE *fp;
  size_t errorcode;
  INT4 fileno=0;   
#ifdef HAVE_GLOB_H
  glob_t globbuf;
#endif
  LIGOTimeGPS starttime;
  INT4 last_time_used=0;
  REAL8 thisSFTtime;
  INT4 blockno=0;

  INITSTATUS (status, "InitFStat", rcsid);
  ATTATCHSTATUSPTR (status);

  /* ----------------------------------------------------------------------
   * do some sanity checks on the user-input before we proceed 
   */
  if (!uvar_mergedSFTFile && (uvar_startTime>-1.0e308 || uvar_endTime<1.0e308)) {
    LALPrintError ( "\nThe --startTime and --endTIme options may ONLY be used\n"
		    "with merged SFT files, specified with the -B option.\n"
		    "Try ./ComputeFStatistic -h \n\n");
    ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);      
  }
  
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
  /*----------------------------------------------------------------------*/

  /* set the current working directory */
  if(chdir(uvar_workingDir) != 0)
    {
      fprintf(stderr, "in Main: unable to change directory to %s\n", uvar_workingDir);
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
    }


  /*----------------------------------------------------------------------
   * set up and check ephemeris-file locations, and SFT input data
   */
#if USE_BOINC
#define EPHEM_EXT ""
#else
#define EPHEM_EXT ".dat"
#endif  
  if (LALUserVarWasSet (&uvar_ephemDir) )
    {
      sprintf(cfg->EphemEarth, "%s/earth%s" EPHEM_EXT, uvar_ephemDir, uvar_ephemYear);
      sprintf(cfg->EphemSun, "%s/sun%s" EPHEM_EXT, uvar_ephemDir, uvar_ephemYear);
    }
  else
    {
      sprintf(cfg->EphemEarth, "earth%s" EPHEM_EXT, uvar_ephemYear);
      sprintf(cfg->EphemSun, "sun%s" EPHEM_EXT,  uvar_ephemYear);
    }
  
#if BOINC_COMPRESS
  /* logic: look for files 'earth.zip' and 'sun,zip'.  If found, use
     boinc_resolve() to get the 'real' file name and then unzip.  If
     not found, look for files named 'earth' and 'sun', use
     boinc_resolve() to get the 'real' file names, and use those
     instead.
  */
  if (uvar_useCompression) {
    char zippedname[256];
    int boinczipret;
    
    /* see if there is a softlink to earth.zip */
    if (!boinc_resolve_filename("earth.zip", zippedname,sizeof(zippedname))) {
      /* if there is, unzip it into the current directory */
      if ((boinczipret=boinc_zip(UNZIP_IT, zippedname, "./"))) {
	fprintf(stderr, "Error in unzipping file %s to earth.  Return value %d\n", zippedname, boinczipret);
	boinc_finish(COMPUTEFSTAT_EXIT_CANTUNZIP);
      }
    }
    /* see if there is a softlink to sun.zip */
    if (!boinc_resolve_filename("sun.zip", zippedname, sizeof(zippedname))) {
      /* if there is, unzip it into the current directory */
      if ((boinczipret=boinc_zip(UNZIP_IT, zippedname, "./"))) {
	fprintf(stderr, "Error in unzipping file %s to sun.  Return value %d\n", zippedname, boinczipret);
	boinc_finish(COMPUTEFSTAT_EXIT_CANTUNZIP);
      }
    }
  }
#endif

#if USE_BOINC
  /* resolve the name of the ephemeris-files.  Does nothing if files are NOT softlinks */
  use_boinc_filename0(cfg->EphemEarth);
  use_boinc_filename0(cfg->EphemSun);
  
#endif  /* not USE_BOINC */
  
  if (uvar_mergedSFTFile) {
    long k=0;
#if USE_BOINC
    use_boinc_filename1(&(uvar_mergedSFTFile));
#endif
    if (!(fp=fp_mergedSFT=fopen(uvar_mergedSFTFile,"rb"))){
      fprintf(stderr,"Unable to open SFT file %s\n", uvar_mergedSFTFile);
      ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
    }

    fileno = 0;
    while (fread((void*)&header,sizeof(header),1,fp) == 1) {
      char tmp[256];

      /* prepare memory for another filename */
      if ( (cfg->filelist = (CHAR**)LALRealloc( cfg->filelist, (fileno+1) * sizeof(CHAR*))) == NULL) {
	ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
      }
      /* store a "file name" composed of merged name + block number */
      sprintf(tmp, "%s (block %d)", uvar_mergedSFTFile, ++blockno);
      if ( (cfg->filelist[fileno] = (CHAR*)LALCalloc (1, strlen(tmp)+1)) == NULL) {
	ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
      }
      strcpy(cfg->filelist[fileno],tmp);
      
      /* check that data is correct endian order and swap if needed */
      if (-1 == reverse_endian) {
	if (header.endian==1.0)
	  reverse_endian=0;
	else
	  reverse_endian=1;
      }
      
      if (reverse_endian)
	swapheader(&header);

      if (header.endian!=1.0) {
	fprintf(stderr,"First object in file %s is not (double)1.0!\n",cfg->filelist[fileno]);
	fprintf(stderr,"The file might be corrupted\n\n");
	ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
      }
      
      /* if this is the first SFT, save initial time */
      if (fileno==0) {
	cfg->Ti=header.gps_sec;
	starttime.gpsSeconds = header.gps_sec;
	starttime.gpsNanoSeconds = header.gps_nsec;
      }
      /* increment file no and pointer to the start of the next header */
      thisSFTtime=(REAL8)header.gps_sec+(1.e-9)*(REAL8)header.gps_nsec;
      if (uvar_startTime<=thisSFTtime && thisSFTtime<uvar_endTime) {
	fileno++;
	last_time_used=header.gps_sec;
      }
      else
	LALFree(cfg->filelist[fileno]);

      k=header.nsamples*8;
      fseek(fp,k,SEEK_CUR);
    }
    /* save final time and time baseline */
    cfg->Tf = (INT4)(last_time_used + header.tbase);  /* FINAL TIME */
    cfg->tsft=header.tbase;  /* Time baseline of SFTs */
    
    /* NOTE: we do NOT close fp here.  If we are using merged SFT file
       for data, we keep it open because we'll need it again in
       ReadSFTData().
    */
    
  }  /* if mergedSFTFile */
  else {

    strcpy(command, uvar_DataDir);
    strcat(command,"/*");
    
    strcat(command, uvar_BaseName);
    strcat(command,"*");
    
#ifdef HAVE_GLOB_H
    globbuf.gl_offs = 1;
    glob(command, GLOB_ERR, NULL, &globbuf);
    
    /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */
    
    if(globbuf.gl_pathc==0)
      {
	LALPrintError ("\nNo SFTs in directory %s ... Exiting.\n\n", uvar_DataDir);
	ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
      }
    /* prepare memory for all filenames */
    if ( (cfg->filelist = (CHAR**)LALCalloc(globbuf.gl_pathc, sizeof(CHAR*))) == NULL) {
      ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }
    while ((UINT4)fileno < (UINT4)globbuf.gl_pathc) 
      {
	if ( (cfg->filelist[fileno] = (CHAR*)LALCalloc(1, strlen(globbuf.gl_pathv[fileno])+1)) == NULL) {
	  ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
	}
	strcpy(cfg->filelist[fileno],globbuf.gl_pathv[fileno]);
	fileno++;
      }
    globfree(&globbuf);
#endif
  }
  cfg->SFTno=fileno; /* remember this is 1 more than the index value */

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

    /* check that data is correct endian order and swap if needed */
    if (-1 == reverse_endian) {
      if (header.endian==1.0)
	reverse_endian=0;
      else
	reverse_endian=1;
    }
    
    if (reverse_endian)
      swapheader(&header);
    
    /* check that data is correct endian order */
    if (header.endian!=1.0)
      {
	fprintf(stderr,"First object in file %s is not (double)1.0!\n",cfg->filelist[0]);
	fprintf(stderr,"The file might be corrupted\n\n");
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
    if (reverse_endian)
      swapheader(&header);

    /* check that data is correct endian order */
    if (header.endian!=1.0)
      {
	fprintf(stderr,"First object in file %s is not (double)1.0!\n",cfg->filelist[fileno-1]);
	fprintf(stderr,"The file might be corrupted\n\n");
	ABORT (status, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
      }
    fclose(fp);

    /* FINAL TIME */
    cfg->Tf = (INT4) (header.gps_sec+header.tbase);  
    /* Time baseline of SFTs */
    cfg->tsft=header.tbase;  
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

  /* if user has not input demodulation frequency resolution; set to 1/2*Tobs */
  if( !LALUserVarWasSet (&uvar_dFreq) ) 
    cfg->dFreq=1.0/(2.0*header.tbase*cfg->SFTno);
  else
    cfg->dFreq = uvar_dFreq;

  cfg->FreqImax=(INT4)(uvar_FreqBand/cfg->dFreq+.5)+1;  /*Number of frequency values to calculate F for */
    
  /* if user has not input demodulation frequency resolution; set to 1/Tobs */
  if( !LALUserVarWasSet (&uvar_df1dot) ) 
    uvar_df1dot=1.0/(2.0*header.tbase*cfg->SFTno*(cfg->Tf - cfg->Ti));

  if (LALUserVarWasSet (&uvar_f1dotBand) && (uvar_f1dotBand != 0) )
    cfg->SpinImax=(int)(uvar_f1dotBand/uvar_df1dot+.5)+1;  /*Number of spindown values to calculate F for */
  else
    cfg->SpinImax = 0;

  cfg->nsamples=header.nsamples;    /* # of freq. bins */

  cfg->ifmax = (INT4) ceil((1.0+DOPPLERMAX)*(uvar_Freq+uvar_FreqBand)*cfg->tsft)+uvar_Dterms;
  cfg->ifmin = (INT4) floor((1.0-DOPPLERMAX)*uvar_Freq*cfg->tsft)-uvar_Dterms;

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
  {
    BOOLEAN haveSkyRegion, haveAlphaDelta, haveGridFile, useGridFile, haveMetric, useMetric;

    haveSkyRegion  = (uvar_skyRegion != NULL);
    haveAlphaDelta = (LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta) );
    haveGridFile   = (uvar_skyGridFile != NULL);
    useGridFile   = (uvar_gridType == GRID_FILE);
    haveMetric     = (uvar_metricType > LAL_PMETRIC_NONE);
    useMetric     = (uvar_gridType == GRID_METRIC);

#if USE_BOINC
    if (haveGridFile)
      use_boinc_filename1(&(uvar_skyGridFile));
    /*
      if storage is allocated dynamically, use this instead!
          use_boinc_filename0(uvar_skyGridFile);
    */
#endif

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
	cfg->skyRegion = (CHAR*)LALCalloc(1, strlen(uvar_skyRegion)+1);
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
      
	cfg->skyRegion = (CHAR*)LALMalloc (512); /* should be enough for max 4 points... */
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

    cfg->edat = (EphemerisData*)LALCalloc(1, sizeof(EphemerisData));
    cfg->edat->ephiles.earthEphemeris = cfg->EphemEarth;
    cfg->edat->ephiles.sunEphemeris = cfg->EphemSun;

    TRY (LALLeapSecs(status->statusPtr, &leap, &starttime, &formatAndAcc), status);
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
    DemodParams = (DemodPar *)LALCalloc(1, sizeof(DemodPar));
  
    /* space for sky constants */
    /* Based on maximum index for array of as and bs sky constants as from ComputeSky.c */
    k = 4*(cfg->SFTno-1) + 4; 
    DemodParams->skyConst = (REAL8 *)LALMalloc(k*sizeof(REAL8));
    /* space for spin down params */
    DemodParams->spinDwn = (REAL8 *)LALMalloc(sizeof(REAL8));
  
  }/* end: init AM- and demod-params */

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

    if ( (fname = (CHAR*)LALCalloc(len,1)) == NULL) {
      ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
    }
    strcpy (fname, head);
    if (uvar_outputLabel)
      strcat (fname, uvar_outputLabel);
    strcat (fname, ".log");

    if ( (fplog = fopen(fname, "wb" )) == NULL) {
      LALPrintError ("\nFailed to open log-file '%f' for writing.\n\n", fname);
      LALFree (fname);
      ABORT (stat, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
    }

    /* write out a log describing the complete user-input (in cfg-file format) */
    TRY (LALUserVarGetLog(stat->statusPtr, &logstr,  UVAR_LOGFMT_CFGFILE), stat);

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
    system (command);	/* we currently don't check this. If it fails, we assume that */
    			/* one of the system-commands was not available, and */
    			/* therefore the CVS-versions will simply not be logged */

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

  INT4 k;

  INITSTATUS (status, "Freemem", rcsid);
  ATTATCHSTATUSPTR (status);

  /* Free SFTData and filenames */
  for (k=0;k<GV.SFTno;k++)
    {
      /* data */
      LALFree(SFTData[k]->fft->data->data);
      LALFree(SFTData[k]->fft->data);
      LALFree(SFTData[k]->fft);
      LALFree(SFTData[k]);

      /* filenames */
      LALFree (GV.filelist[k]);
    } /* for k < SFTno */
  LALFree(SFTData);
  SFTData = NULL;

  LALFree(GV.filelist);
  GV.filelist = NULL;

  /* Free timestamps */
  LALFree(timestamps);
  timestamps = NULL;

  LALFree(Fstat.F);
  Fstat.F = NULL;

  if(uvar_EstimSigParam) 
    {
      LALFree(Fstat.Fa);
      Fstat.Fa = NULL;
      LALFree(Fstat.Fb);
      Fstat.Fb = NULL;
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
      DemodParams = NULL;
    }
     

  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars(status->statusPtr), status);

  if (GV.skyRegion)
    {
      LALFree ( GV.skyRegion );
      GV.skyRegion = NULL;
    }

  /* this comes from clusters.c */
  if (highFLines->clusters) { 
    LALFree(highFLines->clusters);
    highFLines->clusters = NULL;
  }
  if (highFLines->Iclust) {
    LALFree(highFLines->Iclust);
    highFLines->Iclust = NULL;
  }
  if (highFLines->NclustPoints) {
    LALFree(highFLines->NclustPoints);
    highFLines->NclustPoints = NULL;
  }


  /* Free ephemeris data */
  LALFree(GV.edat->ephemE);
  LALFree(GV.edat->ephemS);
  LALFree(GV.edat);
  GV.edat = NULL;

#if USE_BOINC
  /* free buffer used for fstat.  Its safe to do this because we already did fclose(fpstat) earlier */
  if (fstatbuff) {
    LALFree(fstatbuff);
    fstatbuff = NULL;
  }
#endif
  
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

/*******************************************************************************/

/** Print the values of (f,FF) above a certain threshold 
 * in 2F, called 2Fthr.  If there are more than ReturnMaxN of these, 
 * then it simply returns the top ReturnMaxN of them. If there are 
 * none, then it returns none.  It also returns some basic statisical
 * information about the distribution of 2F: the mean and standard 
 * deviation.
 * Returns zero if all is well, else nonzero if a problem was encountered. 
 * Basic strategy: sort the array by values of F, then look at the 
 * top ones. Then search for the points above threshold. 
 */
INT4 PrintTopValues(REAL8 TwoFthr, INT4 ReturnMaxN)
{

  INT4 *indexes,i,j,iF,N;
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
	    ReturnMaxN, GV.FreqImax);
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
  {
    INT4 ntop,err;
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
  err=fprintf(fpmax,"%10.5f %10.8f %10.8f    %d %10.5f %10.5f %10.5f\n",uvar_Freq,
	      Alpha, Delta, GV.FreqImax-N, mean, std, 2.0*log2*Fstat.F[indexes[0]]);
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
  if( (outfile = fopen("F.txt","wb")) == NULL) {
    ABORT (stat, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
  }
#endif
  /*  file contanis freq, PSD, noise floor,lines */
#ifdef FILE_FLINES  
  if( (outfile1 = fopen("FLines.txt","wb")) == NULL) {
    ABORT (stat, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
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
#ifdef FILE_SPRNG  
  FILE *outfile;
#endif
  INT4 i,j,m,lpc,il;                         /* loop indices */
  INT4 Ntot,nbins=GV.ifmax-GV.ifmin+1;   /* Number of points in SFT's */
  REAL8Vector *Sp=NULL, *RngMdnSp=NULL;   /* |SFT|^2 and its rngmdn  */
  REAL8 B;                          /* SFT Bandwidth */
  REAL8 deltaT,norm,*N, *Sp1;
  INT2 windowSize=50;                  /* Running Median Window Size*/
  REAL4 xre,xim,xreNorm,ximNorm;

  INITSTATUS( stat, "NormaliseSFTDataRngMdn", rcsid );
  ATTATCHSTATUSPTR (stat);


  /* The running median windowSize in this routine determines 
     the sample bias which, instead of log(2.0), must be 
     multiplied by F statistics.
  */

  if ( !uvar_SignalOnly ) {
    TRY ( LALRngMedBias (stat->statusPtr, &medianbias, windowSize), stat);
  }

  TRY ( LALDCreateVector(stat->statusPtr, &Sp, (UINT4)nbins), stat);
  TRY ( LALDCreateVector(stat->statusPtr, &RngMdnSp, (UINT4)nbins), stat);

  nbins=(INT2)nbins;	/* RP: What the heck is this????? (FIXME) */

  if( (N = (REAL8 *) LALCalloc(nbins,sizeof(REAL8))) == NULL) {
    ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  if( (Sp1 = (REAL8 *) LALCalloc(nbins,sizeof(REAL8))) == NULL) {
    ABORT (stat, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
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
	Sp->data[j]=((REAL8)xre)*((REAL8)xre)+((REAL8)xim)*((REAL8)xim);
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
  if( (outfile = fopen("SpRng.txt","wb")) == NULL) {
    ABORT (stat, COMPUTEFSTATC_ESYS, COMPUTEFSTATC_MSGESYS);
  } 


  for (j=0;j<nbins;j++){
    Sp1[j]=2.0*Sp1[j]/(1.0*GV.SFTno);
    fprintf(outfile,"%f %E \n",(GV.ifmin+j)/GV.tsft,Sp1[j]); 
  }
  
  fclose(outfile);
#endif

  
  LALFree(N);
  LALFree(Sp1);
  TRY ( LALDDestroyVector(stat->statusPtr, &RngMdnSp), stat);
  TRY ( LALDDestroyVector(stat->statusPtr, &Sp), stat);

  DETATCHSTATUSPTR(stat);
  RETURN(stat);
  
} /* NormaliseSFTDataRngMed() */


/* swap 2, 4 or 8 bytes.  Point to low address */
static void swap2(char *location){
  char tmp=*location;
  *location=*(location+1);
  *(location+1)=tmp;
  return;
}
  
static void swap4(char *location){
  char tmp=*location;
  *location=*(location+3);
  *(location+3)=tmp;
  swap2(location+1);
  return;
}
  
static void swap8(char *location){
  char tmp=*location;
  *location=*(location+7);
  *(location+7)=tmp;
  tmp=*(location+1);
  *(location+1)=*(location+6);
  *(location+6)=tmp;
  swap4(location+2);
  return;
}

void swapheader(struct headertag *thisheader) {
  swap8((char *)&(thisheader->endian));
  swap8((char *)&(thisheader->tbase));
  swap4((char *)&(thisheader->gps_sec));
  swap4((char *)&(thisheader->gps_nsec));
  swap4((char *)&(thisheader->firstfreqindex));
  swap4((char *)&(thisheader->nsamples));
  return;
}

/*******************************************************************************/
/* BOINC-specific functions follow here */
#if USE_BOINC
void use_boinc_filename0(char *orig_name ) {
  char resolved_name[512];
  if (boinc_resolve_filename(orig_name, resolved_name, sizeof(resolved_name))) {
    fprintf(stderr, 
	    "Can't resolve file \"%s\"\n"
	    "If running a non-BOINC test, create [INPUT] or touch [OUTPUT] file\n",
	    orig_name);

    boinc_finish(2);
  }
  strcpy(orig_name, resolved_name);
  return;
}

void use_boinc_filename1(char **orig_name ) {
  char resolved_name[512];
  if (boinc_resolve_filename(*orig_name, resolved_name, sizeof(resolved_name))) {
    fprintf(stderr, 
	    "Can't resolve file \"%s\"\n"
	    "If running a non-BOINC test, create [INPUT] or touch [OUTPUT] file\n",
	    *orig_name);
    boinc_finish(2);
  }
  LALFree(*orig_name);
  *orig_name = (CHAR*) LALCalloc(strlen(resolved_name)+1,1);
  strcpy(*orig_name, resolved_name);
  return;
}

int globargc=0;
char **globargv=NULL;

void worker() {

#if (BOINC_GRAPHICS == 2) 
  if (graphics_lib_handle) {
    if (!(set_search_pos_hook = dlsym(graphics_lib_handle,"set_search_pos"))) {
      fprintf(stderr, "unable to resolve set_search_pos(): %s\n", dlerror());
      boinc_finish(1);
    }
    if (!(fraction_done_hook = dlsym(graphics_lib_handle,"fraction_done"))) {
      fprintf(stderr, "unable to resolve fraction_done(): %s\n", dlerror());
      boinc_finish(1);
    }
  }
  else
    fprintf(stderr,"graphics_lib_handle NULL: running without graphics\n");
#endif

#ifndef RUN_POLKA
  int retval=boincmain(globargc,globargv);
  Outputfilename=Fstatsfilename;
#else
  int a1,a2,retval;
  CHAR ckptfname1[260];
  CHAR Fstatsfilename1[260];

  /* find first // delimiter */ 
  for(a1=0;(a1<globargc)&&(strncmp(globargv[a1],"//",3));a1++);
  if(a1==globargc)
    cfsRunNo = 0;
  else
    cfsRunNo = 1;
  retval=boincmain(a1,globargv);
  /* if there was no //, globargc==a1 and this is old-style command line */
  if(a1<globargc) {
	/* remember first file names */
    strncpy(Fstatsfilename1,Fstatsfilename,260);
    strncpy(ckptfname1,ckp_fname,260);
    if (!retval){
      /* find second // delimiter */ 
      for(a2=a1+1;(a2<globargc)&&(strncmp(globargv[a2],"//",3));a2++);
      cfsRunNo = 2;
      if(a2==globargc)
        retval=COMPUTEFSTAT_EXIT_NOPOLKADEL;
      else
        retval=boincmain(a2-a1,&(globargv[a1]));
      if (!retval)
        retval=polka(globargc-a2, &(globargv[a2]));
    }
    /* remove checkpoint-files */
    remove (ckp_fname);
    remove (ckptfname1);
    /* keep Fstats files while testing - should be deleted as
       temp files with the BOINC slots directory anyway
      remove (Fstatsfilename);
      remove (Fstatsfilename1);
    */
  } else {
    Outputfilename=Fstatsfilename;
    remove (ckp_fname);
  }
#endif
#if BOINC_COMPRESS
  /* compress the file if it exists */
  if (uvar_useCompression && !retval) {
    int boinczipret;
    boinczipret=boinc_zip(ZIP_IT, "temp.zip" , Outputfilename);
    if (boinczipret) {
      fprintf(stderr, "Error in zipping file %s to temp.zip.  Return value %d\n", Outputfilename, boinczipret);
      boinc_finish(COMPUTEFSTAT_EXIT_CANTZIP);
    }
    if ((boinczipret=boinc_rename("temp.zip", Outputfilename))) {
      fprintf(stderr, "Error in renaming file temp.zip to %s.  rename() returned %d\n", Outputfilename, boinczipret);
      boinc_finish(COMPUTEFSTAT_EXIT_CANTRENAME);
    }
  } /* if useCompression && ok */
#endif
  boinc_finish(retval);
  return;
}

int main(int argc, char *argv[]){

  int skipsighandler=0;
  
#if defined(__GNUC__)
  /* see if user has created a DEBUG_CFS file at the top level... */
  FILE *fp_debug=NULL;
  
  if ((fp_debug=fopen("../../DEBUG_CFS", "r")) || (fp_debug=fopen("./DEBUG_CFS", "r"))) {
    
    char commandstring[256];
    char resolved_name[256];
    char *ptr;
    pid_t process_id=getpid();
    
    fclose(fp_debug);
    fprintf(stderr, "Found ../../DEBUG_CFS file, so trying real-time debugging\n");
    
    /* see if the path is absolute or has slashes.  If it has
       slashes, take tail name */
    if ((ptr = strrchr(argv[0], '/'))) {
      ptr++;
    } else {
      ptr = argv[0];
    }
    
    /* if file name is an XML soft link, resolve it */
    if (boinc_resolve_filename(ptr, resolved_name, 256))
      fprintf(stderr, "Unable to boinc_resolve_filename(%s), so no debugging\n", ptr);
    else {
      skipsighandler=1;
      sprintf(commandstring,"ddd %s %d &", resolved_name ,process_id);
      system(commandstring);
      sleep(20);
    }
  }
#endif

  globargc=argc;
  globargv=argv;
    
  /* install signal handler (for ALL threads) for catching
     Segmentation violations, floating point exceptions, Bus
     violations and Illegal instructions */
  if (!skipsighandler) {
    if (signal(SIGSEGV, sighandler)==SIG_IGN)
      signal(SIGSEGV, SIG_IGN);
    if (signal(SIGFPE, sighandler)==SIG_IGN)
      signal(SIGFPE, SIG_IGN);
#ifndef _MSC_VER
    if (signal(SIGBUS, sighandler)==SIG_IGN)
      signal(SIGBUS, SIG_IGN);
#endif
    if (signal(SIGILL, sighandler)==SIG_IGN)
      signal(SIGILL, SIG_IGN);
  }
  
#if (BOINC_GRAPHICS == 1)
  set_search_pos_hook = set_search_pos;
  fraction_done_hook = &fraction_done;
#endif

#ifdef _WIN32
  boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED |
			 BOINC_DIAG_HEAPCHECKENABLED |
			 BOINC_DIAG_ARCHIVESTDERR |
			 BOINC_DIAG_REDIRECTSTDERR |
			 BOINC_DIAG_TRACETOSTDERR);
#endif

  /* boinc_init() or boinc_init_graphics() needs to be run before any
   * boinc_api functions are used */
  
#if BOINC_GRAPHICS>0
  {
    int retval;
#if BOINC_GRAPHICS==2
    /* Try loading screensaver-graphics as a dynamic library.  If this
       succeeds then extern void* graphics_lib_handle is set, and can
       be used with dlsym() to resolve symbols from that library as
       needed.
    */
    retval=boinc_init_graphics_lib(worker, argv[0]);
#endif /* BOINC_GRAPHICS==2 */
#if BOINC_GRAPHICS==1
    /* no dynamic library, just call boinc_init_graphics() */
    retval = boinc_init_graphics(worker);
#endif /* BOINC_GRAPHICS==1 */
    if (retval)
      fprintf(stderr,"boinc_init_graphics[_lib]() returned %d: unable to create worker thread\n", retval);
    boinc_finish(1234+retval);
  }
#endif /*  BOINC_GRAPHICS>0 */
    
    boinc_init();
    worker();
    boinc_finish(222);
    /* we never get here!! */
    return 222;
  }
  
  
#ifdef __GLIBC__
  /* needed to define backtrace() which is glibc specific*/
#include <execinfo.h>
#endif /* __GLIBC__ */

/* signal handlers */
void sighandler(int sig){
  void *array[64];
  size_t size;
  
#ifndef _WIN32
  sigset_t signalset;
  sigemptyset(&signalset);
  sigaddset(&signalset, sig);
  pthread_sigmask(SIG_BLOCK, &signalset, NULL);
#endif

  LALStatus *mystat = &global_status;	  /* the only place in the universe where we're
					   * allowed to use the global_status struct !! */
  /* lets start by ignoring ANY further occurences of this signal
     (hopefully just in THIS thread, if truly implementing POSIX threads */
  
  fprintf(stderr, "Application caught signal %d\n", sig);
  if (mystat)
    fprintf(stderr, "Stack trace of LAL functions in worker thread:\n");
  while (mystat) {
    fprintf(stderr, "%s at line %d of file %s\n", mystat->function, mystat->line, mystat->file);
    if (!(mystat->statusPtr)) {
      const char *p=mystat->statusDescription;
      fprintf(stderr, "At lowest level status code = %d, description: %s\n", mystat->statusCode, p?p:"NO LAL ERROR REGISTERED");
    }
    mystat=mystat->statusPtr;
  }
  
#ifdef __GLIBC__
  /* now get TRUE stacktrace */
  size = backtrace (array, 64);
  fprintf(stderr, "Obtained %zd stack frames for this thread.\n", size);
  fprintf(stderr, "Use gdb command: 'info line *0xADDRESS' to print corresponding line numbers.\n");
  backtrace_symbols_fd(array, size, fileno(stderr));
#endif /* __GLIBC__ */
  /* sleep a few seconds to let the OTHER thread(s) catch the signal too... */
  sleep(5);
  boinc_finish(4321);
  return;
} /* sighandler */
#endif /*USE_BOINC*/


/** Check presence and consistency of checkpoint-file and use to set loopcounter if valid.

 *  The name of the checkpoint-file is <fname>.ckp
 *  @param[OUT] loopcounter	number of completed loops (refers to main-loop in main())
 *  @param[OUT] checksum        checksum of file (up the bytecounter bytes)
 *  @param[OUT] bytecounter	bytes nominally written to fstats file (for consistency-check)
 *  @param[IN]  fstat_fname	Name of Fstats-file. 
 */
void
getCheckpointCounters(LALStatus *stat, UINT4 *loopcounter, UINT4 *checksum, long *bytecounter, const CHAR *fstat_fname, const CHAR *ckpfn)
{
  FILE *fp;
  UINT4 lcount; 	/* loopcounter */
  long bcount;		/* and bytecounter read from ckp-file */
  long flen;
  char lastnewline='\0';
  UINT4 cksum;           /* as read for .ckp file */
  UINT4 computecksum;    /* as computed from contents of Fstats file */
  int i;
#if 0
  int savelaldebuglevel=lalDebugLevel;
  lalDebugLevel=1;
#endif
 
  INITSTATUS( stat, "getChkptCounters", rcsid );
  ASSERT ( fstat_fname, stat, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( ckpfn, stat, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  /* if anything goes wrong in here: start main-loop from beginning  */
  *loopcounter = 0;	
  *bytecounter = 0;
  
  /* try opening checkpoint-file read-only */
  if (lalDebugLevel) LALPrintError("Checking presence of checkpoint-file ... ");
  if (!(fp = fopen(ckpfn, "rb"))) {
    if (lalDebugLevel) LALPrintError ("none found. \nStarting main-loop from beginning.\n");
    RETURN(stat);
  }
  
  /* try reading checkpoint-counters: three INT's loopcounter, checksum, and fstat_bytecounter */
  if (lalDebugLevel) LALPrintError ("found! \nTrying to read checkpoint counters from it...");
  if ( 4 != fscanf (fp, "%" LAL_UINT4_FORMAT " %" LAL_UINT4_FORMAT " %ld\nDONE%c", &lcount, &cksum, &bcount, &lastnewline) || lastnewline!='\n') {
    if (lalDebugLevel) LALPrintError ("failed! \nStarting main-loop from beginning.\n");
    goto exit;
  }
  fclose( fp );
  
  /* checkpoint-file read successfully: check consistency with fstats-file */
  if (lalDebugLevel) LALPrintError ("ok.\nChecking if fstats-file is ok ...");
  if (!(fp = fopen(fstat_fname, "rb"))) {
    if (lalDebugLevel) LALPrintError ("none found.\nStarting main-loop from beginning.\n");
    RETURN(stat);
  }


  /* seek to end of fstats file */
  if (fseek( fp, 0, SEEK_END)) {	/* something gone wrong seeking .. */
    if (lalDebugLevel) LALPrintError ("broken fstats-file.\nStarting main-loop from beginning.\n");
    goto exit;
  }
  
  flen = ftell(fp);

#ifdef USE_BOINC
  fprintf(stderr,"Resuming computation at %d/%ld/%ld\n", lcount, bcount, flen);
#endif

  /* is bytecounter consistent with length of this file? */
  if ( bcount > flen) {
    if (lalDebugLevel) 
      LALPrintError ("seems corrupted: has %ld bytes instead of %ld.\nStarting main-loop from beginning.\n", flen, bcount);
    goto exit;
  }
  
  /* compute checksum */
  computecksum=0;
  if (fseek( fp, 0, SEEK_SET)) {	/* something gone wrong seeking .. */
    if (lalDebugLevel) LALPrintError ("broken fstats-file.\nStarting main-loop from beginning.\n");
    goto exit;
  }
  for (i=0; i<bcount; i++) {
    int onechar=getc(fp);
    if (onechar==EOF)
      goto exit;
    else
      computecksum+=onechar;
  }
  if (computecksum!=cksum) {
    if (lalDebugLevel) 
      LALPrintError ("fstats file seems corrupted: has incorrect checksum.\nStarting main-loop from beginning.\n");
    goto exit;
  }

  if (lalDebugLevel) LALPrintError ("seems ok.\nWill resume from loopcounter = %ld\n", lcount);

  *loopcounter = lcount;
  *bytecounter = bcount;
  *checksum=cksum;

#if 0
  lalDebugLevel=savelaldebuglevel;
#endif

 exit:
  fclose( fp );
  RETURN(stat);
  
} /* getChkptCounters() */ 


#ifdef FILE_AMCOEFFS
void PrintAMCoeffs (REAL8 Alpha, REAL8 Delta, AMCoeffs* amc) {
	UINT4 i;
	FILE*fp;
	fp=fopen("AMCoeffs.dump","a");
	if (fp==NULL) return;
	fprintf(fp,"Alpha=%20.17f, Delta=%20.17f\n", Alpha, Delta);
	fprintf(fp,"a:\n");
	for(i=0;i<amc->a->length;i++)
		fprintf(fp,"%20.17f\n",amc->a->data[i]);
	fprintf(fp,"b:\n");
	for(i=0;i<amc->b->length;i++)
		fprintf(fp,"%20.17f\n",amc->b->data[i]);
	fprintf(fp,"A=%20.17f, B=%20.17f, C=%20.17f, D=%20.17f\n\n",
		amc->A, amc->B, amc->C, amc->D);
	fclose(fp);
}
#endif

/* <lalVerbatim file="LALDemodCP"> */
void OrigLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params) 
/* </lalVerbatim> */
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

  INITSTATUS( status, "LALDemod", rcsid );

  /* catch some obvious programming errors */
  ASSERT ( (Fs != NULL)&&(Fs->F != NULL), status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  if (params->returnFaFb)
    {
      ASSERT ( (Fs->Fa != NULL)&&(Fs->Fb != NULL), status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
    }

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
	  	  	
    Fs->F[i] = (4.0/(M*D))*(B*FaSq + A*FbSq - 2.0*C*FaFb);
    if (params->returnFaFb)
      {
	Fs->Fa[i] = Fa;
	Fs->Fb[i] = Fb;
      }


  }
  /* Clean up */
  LALFree(tempInt1);
  LALFree(xSum);
  LALFree(ySum);
  
  LALFree(sinVal);
  LALFree(cosVal);

  RETURN( status );

} /* OrigLALDemod() */


#define LD_SMALL	(1.0e-9 / LAL_TWOPI)
#define OOTWOPI		(1.0 / LAL_TWOPI)
#define LUT_RES		64	/* resolution of lookup-table */

/* <lalVerbatim file="LALDemodCP"> */
void TestLALDemod(LALStatus *status, LALFstat *Fs, FFT **input, DemodPar *params) 
/* </lalVerbatim> */
{ 

  INT4 alpha,i;                 /* loop indices */
  REAL8	*xSum=NULL, *ySum=NULL;	/* temp variables for computation of fs*as and fs*bs */
  INT4 s;		        /* local variable for spinDwn calcs. */
  REAL8	xTemp;	                /* temp variable for phase model */
  REAL8	deltaF;	                /* width of SFT band */
  INT4	k1;	                /* defining the sum over which is calculated */
  UINT4 k;
  REAL8 *skyConst;	        /* vector of sky constants data */
  REAL8 *spinDwn;	        /* vector of spinDwn parameters (maybe a structure? */
  INT4	spOrder;	        /* maximum spinDwn order */
  REAL8	realXP, imagXP; 	/* temp variables used in computation of */
  INT4	nDeltaF;	        /* number of frequency bins per SFT band */
  INT4	sftIndex;	        /* more temp variables */
  REAL8 realQ, imagQ;
  INT4 *tempInt1;
  UINT4 index;
  REAL8 FaSq;
  REAL8 FbSq;
  REAL8 FaFb;
  COMPLEX16 Fa, Fb;
  UINT4 klim = 2*params->Dterms;
  REAL8 f;
  static REAL8 sinVal[LUT_RES+1], cosVal[LUT_RES+1];        /*LUT values computed by the routine do_trig_lut*/
  static BOOLEAN firstCall = 1;

  REAL8 A=params->amcoe->A;
  REAL8 B=params->amcoe->B;
  REAL8 C=params->amcoe->C;
  REAL8 D=params->amcoe->D;

  UINT4 M=params->SFTno;

  INITSTATUS( status, "LALDemod", rcsid );

  /* catch some obvious programming errors */
  ASSERT ( (Fs != NULL)&&(Fs->F != NULL), status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  if (params->returnFaFb)
    {
      ASSERT ( (Fs->Fa != NULL)&&(Fs->Fb != NULL), status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
    }

  /* variable redefinitions for code readability */
  spOrder=params->spinDwnOrder;
  spinDwn=params->spinDwn;
  skyConst=params->skyConst;
  deltaF=(*input)->fft->deltaF;
  nDeltaF=(*input)->fft->data->length;

  /* res=10*(params->mCohSFT); */
  /* This size LUT gives errors ~ 10^-7 with a three-term Taylor series */
  if ( firstCall )
    {
      for (k=0; k <= LUT_RES; k++)
	{
	  sinVal[k] = sin( (LAL_TWOPI*k)/LUT_RES );
	  cosVal[k] = cos( (LAL_TWOPI*k)/LUT_RES );
	}
      firstCall = 0;
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
	REAL8 tempFreq0, tempFreq1;
	REAL4 tsin, tcos;
	COMPLEX8 *Xalpha=input[alpha]->fft->data->data;
	REAL4 a = params->amcoe->a->data[alpha];
	REAL4 b = params->amcoe->b->data[alpha];
	REAL8 x,y;
	REAL4 realP, imagP;	        /* real and imaginary parts of P, see CVS */

	/* NOTE: sky-constants are always positive!!
	 * this can be seen from there definition (-> documentation)
	 * we will use this fact in the following! 
	 */
	xTemp= f * skyConst[ tempInt1[ alpha ] ] + xSum[ alpha ];	/* >= 0 !! */
	/* this will now be assumed positive, but we double-check this
	 * for now to be sure: 
	 */
	if ( xTemp < 0 )
	  {
	    printf ( "RP is an idiot: the assumption xTemp >= 0 is not generally true!!\n");
	    exit (-1);
	  }

	/* find correct index into LUT -- pick closest point */
	tempFreq0 = xTemp - (UINT4)xTemp;  /* lies in [0, +1) by definition */

	index = (UINT4)( tempFreq0 * LUT_RES + 0.5 );	/* positive! */
	{
	  REAL8 d=LAL_TWOPI*(tempFreq0 - (REAL8)index/(REAL8)LUT_RES);
	  REAL8 d2=0.5*d*d;
	  REAL8 ts=sinVal[index];
	  REAL8 tc=cosVal[index];
		
	  tsin = ts+d*tc-d2*ts;
	  tcos = tc-d*ts-d2*tc-1.0;
	}

	y = - LAL_TWOPI * ( f * skyConst[ tempInt1[ alpha ]-1 ] + ySum[ alpha ] );
	realQ = cos(y);
	imagQ = sin(y);

	/*
	REAL8 yTemp = f * skyConst[ tempInt1[ alpha ]-1 ] + ySum[ alpha ];
	REAL8 yRem = yTemp - (UINT4)yTemp;

	index = (UINT4)( yRem * LUT_RES + 0.5 );
	{
	  REAL8 d = LAL_TWOPI*(yRem - (REAL8)index/(REAL8)LUT_RES);
	  REAL8 d2=0.5*d*d;
	  REAL8 ts = sinVal[index];
	  REAL8 tc = cosVal[index];
		
	  imagQ = ts + d * tc - d2 * ts;
	  imagQ = -imagQ;
	  realQ = tc - d * ts - d2 * tc;
	}
	*/

	k1 = (UINT4)xTemp - params->Dterms + 1;

	sftIndex = k1 - params->ifmin;

	tempFreq1 = tempFreq0 + params->Dterms - 1; 	/* positive if Dterms > 1 (trivial) */

	x = LAL_TWOPI * tempFreq1;	/* positive! */

	/* we branch now (instead of inside the central loop)
	 * depending on wether x can ever become SMALL in the loop or not, 
	 * because it requires special treatment in the Dirichlet kernel
	 */
	if ( tempFreq0 < LD_SMALL ) 
	  {

	    realXP=0.0;
	    imagXP=0.0;

	    /* Loop over terms in Dirichlet Kernel */
	    for(k=0; k < klim ; k++)
	      {
		COMPLEX8 Xalpha_k = Xalpha[sftIndex];
		sftIndex ++;
		/* If x is small we need correct x->0 limit of Dirichlet kernel */
		if( fabs(x) <  SMALL) 
		  {
		    realXP += Xalpha_k.re;
		    imagXP += Xalpha_k.im;
		  }	 
		else
		  {
		    realP = tsin / x;
		    imagP = tcos / x;
		    /* these four lines compute P*xtilde */
		    realXP += Xalpha_k.re * realP;
		    realXP -= Xalpha_k.im * imagP;
		    imagXP += Xalpha_k.re * imagP;
		    imagXP += Xalpha_k.im * realP;
		  }
		
		tempFreq1 --;
		x = LAL_TWOPI * tempFreq1;
		
	      } /* for k < klim */

	  } /* if x could become close to 0 */
	else
	  {
	    COMPLEX8 *Xalpha_k = Xalpha + sftIndex;

	    realXP=0.0;
	    imagXP=0.0;

	    /* Loop over terms in dirichlet Kernel */
	    for(k=0; k < klim ; k++)
	      {
		REAL4 xinv = (REAL4)OOTWOPI / (REAL4)tempFreq1;
		COMPLEX8 Xa = *Xalpha_k;
		Xalpha_k ++;
		tempFreq1 --;
		
		realP = tsin * xinv;
		imagP = tcos * xinv;
		/* these four lines compute P*xtilde */
		realXP += Xa.re * realP - Xa.im * imagP;
		imagXP += Xa.re * imagP + Xa.im * realP;

	      } /* for k < klim */

	  } /* if x cannot be close to 0 */
	

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
	  	  	
    Fs->F[i] = (4.0/(M*D))*(B*FaSq + A*FbSq - 2.0*C*FaFb);
    if (params->returnFaFb)
      {
	Fs->Fa[i] = Fa;
	Fs->Fb[i] = Fb;
      }


  }
  /* Clean up */
  LALFree(tempInt1);
  LALFree(xSum);
  LALFree(ySum);
  
  RETURN( status );

} /* TestLALDemod() */
