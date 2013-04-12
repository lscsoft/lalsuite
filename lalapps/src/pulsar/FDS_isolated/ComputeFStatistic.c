/*
*  Copyright (C) 2007 Bruce Allen, Bernd Machenschalk, David Hammer, Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Xavier Siemens, Scott Koranda, Yousuke Itoh
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

/*********************************************************************************
 * \file
 * \ingroup pulsarApps
 * \author  B. Allen, Y. Ioth, B. Machenschalk, M.A. Papa, R.Prix, X. Siemens
 * \brief
 * Calculate the F-statistic for a given parameter-space of pulsar GW signals.
 * Implements the so-called "F-statistic" as introduced in JKS98.
 *                                                                          
 *          Albert Einstein Institute/UWM - started September 2002   
 *********************************************************************************/

#include "config.h"


/* System includes */
#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_GLOB_H
#include <glob.h>
#endif

#include <string.h>
#include <errno.h>

#define	 __USE_ISOC99 1
#include <math.h>

/* LAL-includes */
#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/AVFactories.h>
#include <lal/RngMedBias.h>
#include <lal/LALDemod.h>
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/LogPrintf.h>

#include <lalapps.h>

/* local includes */
#include "ComputeFStatistic.h"
#include "clusters.h"

#include "FstatToplist.h"



/* this is defined in C99 and *should* be in math.h.  Long term
   protect this with a HAVE_FINITE */
int finite(double);

/*----------------------------------------------------------------------*/
/* conditional compilation-switches */
/*
#define DEBG_FAFB                
#define DEBG_ESTSIGPAR
*/

/* Uncomment the following if you want to activate output f1dot in 'Fstats' */
/* #define OUTPUT_F1DOT  1 */

/*
#define FINDFWHM
#define FILE_FMAX
#define FILE_FLINES
#define FILE_FTXT 
#define FILE_SPRNG
#define FILE_AMCOEFFS
*/

/*----------------------------------------------------------------------*/
/* BOINC stuff for running Einstein@Home */
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

/* for getpid() */
#include <sys/types.h>

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
#define remove boinc_delete_file

int boincmain(int argc, char *argv[]);
void worker(void);

/* hooks for communication with the graphics thread */
void (*set_search_pos_hook)(float,float) = NULL;
int (*boinc_init_graphics_hook)(void (*worker)(void)) = NULL;
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
extern int boinc_init_graphics(void (*worker)(void));
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

struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
} header;
  
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
CHAR *uvar_DataFiles;
BOOLEAN uvar_help;
CHAR *uvar_outputLabel;
CHAR *uvar_outputFstat;
CHAR *uvar_outputLoudest;
CHAR *uvar_skyGridFile;
CHAR *uvar_outputSkyGrid;
CHAR *uvar_workingDir;
BOOLEAN uvar_doCheckpointing;
INT4 uvar_expLALDemod;
REAL8 uvar_startTime;
REAL8 uvar_endTime;
REAL8 uvar_refTime;
CHAR *uvar_outputClusters;
INT4 uvar_RngMedWindow;
#if BOINC_COMPRESS
BOOLEAN uvar_useCompression;
#endif

INT4 uvar_OutputBufferKB;
INT4 uvar_MaxFileSizeKB;
INT4 uvar_NumCandidatesToKeep;

BOOLEAN uvar_projectMetric;
REAL8 uvar_WUfpops;

/*----------------------------------------------------------------------*/
/* some other global variables */

FFT **SFTData=NULL;	/**< SFT Data for LALDemod */
FFT **UpSFTData=NULL;	/**< Upsampled SFT Data for LALDemodFAST */
DemodPar *DemodParams;	/**< Demodulation parameters for LALDemod */
LIGOTimeGPS *timestamps;/**< Time stamps from SFT data */
LALFstat Fstat;		/**< output from LALDemod(): F-statistic and amplitudes Fa and Fb */
AMCoeffs amc; 		/**< amplitude-modulation coefficients (and derived quantities) */
Clusters HFLines;	/**< stores information about outliers/clusters in F-statistic */
Clusters HPLines;	/**< stores information about outliers/clusters in SFT-power spectrum */
Clusters *highSpLines=&HPLines, *highFLines=&HFLines;

REAL8 medianbias=1.0;	/**< median-bias depending on window-size */

FILE *fp_mergedSFT;	/**< input-file containing merged SFTs */
FILE *fpmax;		/**< output-file: maximum of F-statistic over frequency-range */
FILE *fpClusters;	/**< output-file pointer to clustered Fstat output */
FILE *fpFstat;		/**< output-file pointer to *unclustered* Fstat output */

ConfigVariables GV;	/**< global container for various derived configuration settings */
int reverse_endian=-1;	/**< endian order of SFT data.  -1: unknown, 0: native, 1: reversed */
CHAR CFstatFilename[MAXFILENAMELENGTH]; /**< clustered Fstats file name*/
CHAR FstatFilename[MAXFILENAMELENGTH]; 	/**< (unclustered) Fstats file name*/
CHAR ckp_fname[MAXFILENAMELENGTH+4];    /**< filename of checkpoint-file, global for polka */
CHAR *Outputfilename;	/**< Name of output file, either Fstats- or Polka file name*/
INT4 cfsRunNo = 0;	/**< CFS run-number: 0=run only once, 1=first run, 2=second run */

toplist_t* toplist = NULL;
char *fstatbuff = NULL;

FstatOutputEntry empty_FstatOutputEntry;

INT4 maxSFTindex = 0;

/*----------------------------------------------------------------------*/
/* local prototypes */
#ifdef __cplusplus
extern "C" {
#endif
  int main(int argc,char *argv[]);
  void initUserVars (LALStatus *);
  INT4 ReadSFTData (void);
  INT4 UpsampleSFTData (void);
  void InitFStat (LALStatus *, ConfigVariables *cfg);
  void CreateDemodParams (LALStatus *, PulsarDopplerParams dopplerpos);
  void CreateNautilusDetector (LALStatus *, LALDetector *Detector);
  void Freemem (LALStatus *);
  void EstimateFLines(LALStatus *);
  void NormaliseSFTDataRngMdn (LALStatus *, INT4 windowSize);
  INT4 EstimateSignalParameters(INT4 * maxIndex);
  int writeFLines(INT4 *maxIndex, PulsarDopplerParams searchpos, FILE *fpOut);
  int writeFLinesCS(INT4 *maxIndex, PulsarDopplerParams searchpos, FILE *fpOut, long*bytecount, UINT4*checksum);
  INT4 PrintTopValues(REAL8 TwoFthr, INT4 ReturnMaxN, PulsarDopplerParams searchpos);
  int compare(const void *ip, const void *jp);
  INT4 writeFaFb(INT4 *maxIndex, PulsarDopplerParams searchpos);
  void checkUserInputConsistency (LALStatus *);

  void InitSearchGrid ( LALStatus *, DopplerSkyScanState *scan, ConfigVariables *cfg);

  void WriteFStatLog (LALStatus *, CHAR *argv[]);
  static void swap2(char *location);
  static void swap4(char *location);
  static void swap8(char *location);
  void swapheader(struct headertag *thisheader);
  void getCheckpointCounters(LALStatus *, UINT4 *loopcounter, UINT4 *checksum, 
			     long *bytecounter, const CHAR *fstat_fname, const CHAR *ckpfn);
  
  int debug_dump_commandline (int argc,  char *argv[]);
  int is_zip_file ( const char *fname );
  void ourREPORTSTATUS( LALStatus * );

#ifdef FILE_AMCOEFFS
  void PrintAMCoeffs (REAL8 Alpha, REAL8 Delta, AMCoeffs* amc);
#endif
#ifdef __cplusplus
}
#endif


/*----------------------------------------------------------------------*/
/* some local defines */

#define EPHEM_YEARS  "00-19-DE405"
#define SFT_BNAME  "SFT"

#ifndef TRUE 
#define TRUE (1==1)
#endif
#ifndef FALSE
#define FALSE (1==0)
#endif

#define MYMAX(x,y) ((x) > (y) ? (x) : (y) )

extern int vrbflg;

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


/* freaking LAL's REPORTSTATUS just won't work with any of NDEBUG or 
 * LAL_NDEBUG set, so it's time to use our own...
 */
void
ourREPORTSTATUS( LALStatus *status )
{ 
  LALStatus *ptr;                                                    
  for ( ptr = status; ptr ; ptr = ptr->statusPtr )                  
  {                                                           
    LogPrintf ( LOG_NORMAL, "\nLevel %i: %s\n", ptr->level, ptr->Id );    
    if ( ptr->statusCode )                                    
    {                                                      
      LogPrintf ( LOG_NORMAL, "\tStatus code %i: %s\n", ptr->statusCode,      
                     ptr->statusDescription );                      
    }                                                            
    else                                                        
    {                                                          
      LogPrintf ( LOG_NORMAL, "\tStatus code 0: Nominal\n" );            
    }                                                        
    LogPrintf ( LOG_NORMAL, "\tfunction %s, file %s, line %i\n",      
                   ptr->function, ptr->file, ptr->line );      
  }                                                       
  return;

} /* ourREPORTSTATUS() */



#if USE_BOINC
int BOINC_ERR_EXIT(LALStatus  *stat, const char *func, const char *file, const int line, volatile const char *id) {
  if (stat->statusCode) {
    fprintf(stderr,
            "Level 0: %s\n"
            "\tFunction call `%s' failed.\n"
            "\tfile %s, line %d\n",
            id, func, file, line );
    ourREPORTSTATUS(stat);
    LogPrintf (LOG_CRITICAL, "BOINC_ERR_EXIT(): now calling boinc_finish()\n");
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
  LALStatus *status = &global_status;

  INT4 *maxIndex=NULL;           /* array that contains indexes of maximum of each cluster */
  INT4 spdwn;                    /* counter over spindown-params */
  DopplerSkyScanState thisScan = empty_DopplerSkyScanState; /* current state of the Doppler-scan */
  PulsarDopplerParams dopplerpos;    /* current search-parameters */
  SkyPosition thisPoint;
  UINT4 loopcounter;             /* Checkpoint-counters for restarting checkpointed search */
  long fstat_bytecounter;
  UINT4 fstat_checksum = 0;      /* Checksum of fstats file contents */
  FstatOutputEntry loudest = empty_FstatOutputEntry; /* loudest canidate in search-region */

#ifdef RUN_POLKA
  INT4 fstats_completed = FALSE; /* did we find a completed fstats file? */
#endif

#if USE_BOINC
  CHAR resfname[MAXFILENAMELENGTH]; /* buffer for boinc-resolving config-file name */

  /* to avoid unwanted filename resolution in LALOpenDataFile() */
#ifdef _MSC_VER /* for MSC */
  _putenv("LAL_DATA_PATH=");
#elif HAVE_UNSETENV /* should have been checked in configure */
  unsetenv("LAL_DATA_PATH");
#else
  /* hope that this does the right thing on Systems that don't have unsetenv... */
  putenv("LAL_DATA_PATH");
#endif

  /* set LAL error-handler */
  lal_errhandler = BOINC_ERR_EXIT;
#else
  lal_errhandler = LAL_ERR_EXIT;
#endif

  vrbflg = 1;   /* verbose error-messages */

  /* register all user-variable */
  LAL_CALL (LALGetDebugLevel(status, argc, argv, 'v'), status);
  LAL_CALL (initUserVars(status), status);  

  /* we use lalDebugLevel for our logging */
#if USE_BOINC
  LogSetLevel ( lalDebugLevel );
#else
  LogSetLevel ( lalDebugLevel - 1 );
#endif

  LogPrintf (LOG_NORMAL, "Started search at lalDebugLevel = %d\n", lalDebugLevel);

  debug_dump_commandline (argc, argv);

#if USE_BOINC
  /* handle config file request with boinc_resolve_filename */
  /* NOTE: @configfile must be at the beginning of the command line! */
  if ( argv[1] && argv[1][0] == '@') 
    {
      resfname[0] = '@';
      if (boinc_resolve_filename(argv[1]+1,resfname+1,sizeof(resfname)))
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve config file '%s'\n", argv[1]+1);
      else
        {
          /* hack the command-line: replace config-file by boinc-resolved path */
          argv[1] = resfname;
          /* note: we don't free the previous argv[1] in order to avoid possible problems..*/
        }
    } /* if config-file given as first argument */
#endif

  LAL_CALL (LALUserVarReadAllInput(status,argc,argv),status);       


  if (uvar_help)        /* if help was requested, we're done here */
    return COMPUTEFSTAT_EXIT_USAGE;

  if ( XLALUserVarWasSet ( &uvar_refTime ) )
    XLAL_ERROR (XLAL_EFAILED, "Sorry, using --refTime in CFSv1 is deprecated, the reference time is forced to be = startTime!\n");

  /* This is dangerous for BOINC since it calls system() and makes
     assumptions that might not be true */
#if !USE_BOINC
  /* keep a log-file recording all relevant parameters of this search-run */
  LAL_CALL (WriteFStatLog (status, argv), status);
#endif /* !USE_BOINC */

  /* do some sanity checks on the user-input before we proceed */
  LAL_CALL ( checkUserInputConsistency(status), status);

  /* main initialization of the code: */
  LAL_CALL ( InitFStat(status, &GV), status);
  
  /* ----- initialization of the DopplerSkyScanner to step through paramter space ----- */
  LAL_CALL ( InitSearchGrid(status, &thisScan, &GV), status);

  /* ---------- overload Frequency- and spindown-resolution if input by user ----------*/
  if ( LALUserVarWasSet( &uvar_dFreq ) )
    thisScan.dfkdot[0] = uvar_dFreq;

  if( LALUserVarWasSet( &uvar_df1dot) ) 
    thisScan.dfkdot[1] = uvar_df1dot;

  DemodParams->df   = thisScan.dfkdot[0];

  /* Number of Freq- and spindown values to calculate F for */
  GV.FreqImax = (INT4)(GV.spinRange.fkdotBand[0] / thisScan.dfkdot[0] + 1e-6) + 1;  
  GV.SpinImax = (INT4)(GV.spinRange.fkdotBand[1] / thisScan.dfkdot[1] + 1e-6) + 1;  

  /* debug output about search-parameters */
  {
    const CHAR *str;
    LogPrintf (LOG_DEBUG, "Total number of templates (sky x f x fdot) = %d x %d x %d = %g\n",
	       thisScan.numSkyGridPoints, GV.FreqImax, GV.SpinImax, 
	       1.0 * thisScan.numSkyGridPoints * GV.FreqImax * GV.SpinImax );

    str = GV.skyRegionString;
    LogPrintf (LOG_DETAIL, "skyRegion = '%s'\n", str ? str : "NULL" );
    LogPrintf (LOG_DETAIL, "Frequency-range [%.16g, %.16g]\n", 
	       GV.spinRange.fkdot[0], GV.spinRange.fkdot[0] + GV.spinRange.fkdotBand[0]);
    LogPrintf (LOG_DETAIL, "Spindown-range [%.16g, %.16g]\n", 
	       GV.spinRange.fkdot[1], GV.spinRange.fkdot[1] + GV.spinRange.fkdotBand[1]);
    LogPrintf (LOG_DETAIL, "Grid-spacings: dFreq = %g, df1dot = %g\n\n",
	       thisScan.dfkdot[0], thisScan.dfkdot[1]);
  }

  /* determine smallest required band of frequency-bins for the search-parameters */
  {
    REAL8 f_min, f_max;
    REAL8 f1dot_1, f1dot_2, f1dot_max, df;
    f_min = GV.spinRange.fkdot[0];
    f_max = f_min + GV.spinRange.fkdotBand[0];

    /* correct for spindown-shift of frequency: extend frequency-band */
    f1dot_1 = GV.spinRange.fkdot[1];
    f1dot_2 = f1dot_1 + GV.spinRange.fkdotBand[1];
    f1dot_max = fabs(f1dot_1) > fabs(f1dot_2) ? f1dot_1 : f1dot_2;
    df = f1dot_max * (GV.Tf - GV.Ti);
    if ( df < 0)
      f_min += df;
    else
      f_max += df;

    GV.ifmin = (INT4) floor( (1.0-DOPPLERMAX)* f_min * GV.tsft) - MYMAX(uvar_Dterms, uvar_RngMedWindow/2 + 1 );
    GV.ifmax = (INT4) ceil( (1.0+DOPPLERMAX) * f_max * GV.tsft) + MYMAX(uvar_Dterms, uvar_RngMedWindow/2 + 1);
    maxSFTindex = GV.ifmax - GV.ifmin;
  }

  /* allocate F-statistic arrays */
  Fstat.F =(REAL8*)LALMalloc(GV.FreqImax*sizeof(REAL8));
  if(uvar_EstimSigParam) 
    {
      Fstat.Fa =(COMPLEX16*)LALMalloc(GV.FreqImax*sizeof(COMPLEX16));
      Fstat.Fb =(COMPLEX16*)LALMalloc(GV.FreqImax*sizeof(COMPLEX16));
    } else {
      Fstat.Fa = NULL;
      Fstat.Fb = NULL;
    }
  
  /* ----------------------------------------------------------------------*/

  /* read in SFT-data */
  if (ReadSFTData()) return COMPUTEFSTAT_EXIT_READSFTFAIL;
  
  /* normalize SFTs by running median */
  LAL_CALL (NormaliseSFTDataRngMdn(status, uvar_RngMedWindow), status);

  /* Here's where I need to upsample the SFTs */
  if ( uvar_expLALDemod == 3)
    {
      if (UpsampleSFTData()) return COMPUTEFSTAT_EXIT_UPSAMPLESFTFAIL;
      fprintf(stdout,"Finished upsampling data by factor of %f\n", 1.0 / ( DemodParams->df * GV.tsft) );
    }

#ifdef FILE_FMAX
  {
    CHAR Fmaxfilename[MAXFILENAMELENGTH]; /* Fmax file name*/

    /*   open file */
    strcpy(Fmaxfilename,"Fmax");
    if (uvar_outputLabel)
      strcat(Fmaxfilename,uvar_outputLabel);
#if USE_BOINC
    use_boinc_filename0(Fmaxfilename);
#endif /* USE_BOINC */
    if (!(fpmax=fopen(Fmaxfilename,"wb"))){
      LogPrintf (LAL_ERROR, "in Main: unable to open Fmax file %s\n", Fmaxfilename);
      return COMPUTEFSTAT_EXIT_OPENFMAX;
    }
  }
#endif

  /* ----- prepare cluster-output filename if given and append outputLabel */
  if ( uvar_outputClusters && (strlen(uvar_outputClusters) > 0) )
    {
      strncpy ( CFstatFilename, uvar_outputClusters, sizeof(CFstatFilename) );
      if ( uvar_outputLabel )
	strncat ( CFstatFilename, uvar_outputLabel, sizeof(CFstatFilename) - strlen(CFstatFilename)  - 1 );

      if ( (fpClusters = fopen (CFstatFilename, "wb")) == NULL ) {
	LogPrintf (LOG_CRITICAL, "Failed to open Clusters-file '%s' for writing!\n", CFstatFilename );
	return ( COMPUTEFSTAT_ESYS );
      }
    }
  else {
    strcpy (CFstatFilename,"");
    fpClusters = NULL;
  }

  /* ----- prepare (unclustered) Fstat-output filename if given and apprend outputLabel */ 
  if ( uvar_outputFstat )
    {
      UINT4 len = strlen (uvar_outputFstat );
      if ( uvar_outputLabel )
	len += strlen ( uvar_outputLabel );

      strncpy ( FstatFilename, uvar_outputFstat, sizeof(FstatFilename) );
      if ( uvar_outputLabel )
	strncat ( FstatFilename, uvar_outputLabel, sizeof(FstatFilename) - strlen(FstatFilename) - 1 );
    }
  else
    strncpy ( FstatFilename, "", sizeof(FstatFilename) );


  if ( uvar_NumCandidatesToKeep > 0 )
    create_fstat_toplist(&toplist, uvar_NumCandidatesToKeep);


  /* prepare checkpointing file */
#ifdef CLUSTERED_OUTPUT
  if ( strlen(CFstatFilename) ) 
    strncpy(ckp_fname, CFstatFilename, sizeof(ckp_fname));
#else
  if ( strlen(FstatFilename) ) 
    strncpy(ckp_fname, FstatFilename, sizeof(ckp_fname));
#endif
  else
    strcpy(ckp_fname, "Fstats");
  strncat(ckp_fname, ".ckp", sizeof(ckp_fname) - strlen(ckp_fname) - 1 );

#if USE_BOINC
  /* only boinc_resolve the filename if we run CFS once */
  if (cfsRunNo == 0)
#ifdef CLUSTERED_OUTPUT
    if(strlen(CFstatFilename))
      use_boinc_filename0(CFstatFilename);
#else
    if(strlen(FstatFilename))
      use_boinc_filename0(FstatFilename);
#endif
  /* use_boinc_filename0(ckp_fname); */
#endif 

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
  fstat_checksum =0;

#ifdef RUN_POLKA
  /* look for a Fstat file ending in "%DONE" and quit (boinc)main if found */
#ifdef CLUSTERED_OUTPUT
  fpFstat = fopen( CFstatFilename,"rb");
#else
  fpFstat = fopen( FstatFilename,"rb");
#endif
  if(fpFstat){
    char done[6];
    done[0] = '\0';
    if(!fseek(fpFstat,-6,SEEK_END))
      if(fread(done,6,1,fpFstat)==1)
        if(strncmp(done,"%DONE",5)==0){
          LogPrintf (LOG_NORMAL, "detected finished Fstat file - skipping Fstat run %d\n", cfsRunNo);
	  fstats_completed = TRUE;
        }
    fclose(fpFstat);
    fpFstat = NULL;
  }

  if (!fstats_completed) {
#endif /* RUN_POLKA */

#ifdef CLUSTERED_OUTPUT
  /* Checkpointed information is retrieved from checkpoint-file (if found) */
  if (uvar_doCheckpointing)
    LAL_CALL (getCheckpointCounters( status, &loopcounter, &fstat_checksum, &fstat_bytecounter, 
				     CFstatFilename, ckp_fname ), status); 
  /* allow for checkpointing: 
   * open Fstat file for writing or appending, depending on loopcounter. 
   */
  if ( strlen(CFstatFilename)
       && ( (fpClusters = fopen( CFstatFilename, fstat_bytecounter>0 ? "rb+" : "wb")) == NULL) )
    {
      LogPrintf (LOG_CRITICAL, "in Main: unable to open Fstats file '%s'\n", FstatFilename);
      return COMPUTEFSTAT_EXIT_OPENFSTAT2;
    }

  /* set a buffer large enough that no output is written to disk
   *  unless we fflush().  Policy is fully buffered. Note: the man
   * page says "The setvbuf function may only be used after opening
   * a stream and BEFORE ANY OTHER OPERATIONS HAVE BEEN PERFORMED ON IT."
   */
  if ( fpClusters && uvar_OutputBufferKB )
    {
      if ( (fstatbuff=(char *)LALCalloc(uvar_OutputBufferKB * 1024, 1)))
	setvbuf(fpClusters, fstatbuff, _IOFBF, uvar_OutputBufferKB * 1024);
      else {
	LogPrintf (LOG_CRITICAL, "Failed to allocate memory for Fstats file buffering. Exiting.\n");
	return COMPUTEFSTAT_EXIT_NOMEM;
      }
    }
#else /* checkpointing clustered output */
  /* Checkpointed information is retrieved from checkpoint-file (if found) */
  if (uvar_doCheckpointing)
    LAL_CALL (getCheckpointCounters( status, &loopcounter, &fstat_checksum, &fstat_bytecounter, 
				     FstatFilename, ckp_fname ), status); 

  /* allow for checkpointing: 
   * open Fstat file for writing or appending, depending on loopcounter. 
   */
  if ( strlen(FstatFilename)
       && ( (fpFstat = fopen( FstatFilename, fstat_bytecounter>0 ? "rb+" : "wb")) == NULL) )
    {
      LogPrintf (LOG_CRITICAL, "in Main: unable to open Fstats file '%s'\n", FstatFilename);
      return COMPUTEFSTAT_EXIT_OPENFSTAT2;
    }

  /* set a buffer large enough that no output is written to disk
   *  unless we fflush().  Policy is fully buffered. Note: the man
   * page says "The setvbuf function may only be used after opening
   * a stream and BEFORE ANY OTHER OPERATIONS HAVE BEEN PERFORMED ON IT."
   */
  if ( fpFstat && uvar_OutputBufferKB )
    {
      if ( (fstatbuff=(char *)LALCalloc(uvar_OutputBufferKB * 1024, 1)))
	setvbuf(fpFstat, fstatbuff, _IOFBF, uvar_OutputBufferKB * 1024);
      else {
	LogPrintf (LOG_CRITICAL, "Failed to allocate memory for Fstats file buffering. Exiting.\n");
	return COMPUTEFSTAT_EXIT_NOMEM;
      }
    }
#endif

  /* log-output of checkpoint status */
  if ( uvar_doCheckpointing )
    {
      if ( loopcounter && fstat_bytecounter )
	LogPrintf (LOG_NORMAL, "Resuming computation at (%d/%d/%ld).\n", 
		   loopcounter, fstat_checksum, fstat_bytecounter);
      else
	LogPrintf (LOG_NORMAL, "No usable checkpoint found, starting from beginning.\n");
    }


  /* if we checkpointed we need to "spool forward" to the right entry in the sky-position list */
  if ( loopcounter > 0 ) 
    {
      UINT4 i;
      for (i=0; i < loopcounter; i++) {
	XLALNextDopplerSkyPos( &dopplerpos, &thisScan );
	if (thisScan.state == STATE_FINISHED) {
	  LogPrintf (LOG_CRITICAL, "Checkpointed loopcounter already at the end of main-loop\n");
	  return COMPUTEFSTAT_ECHECKPOINT; 
	}
      }
      /* seek to right point of fstats file (truncate what's left over) */
      if ( fpFstat )
	{
#ifdef CLUSTERED_OUTPUT
	  if ( 0 != fseek( fpClusters, fstat_bytecounter, SEEK_SET) ) 
#else
	  if ( 0 != fseek( fpFstat, fstat_bytecounter, SEEK_SET) ) 
#endif
	    {   /* something gone wrong seeking .. */
	      LogPrintf (LOG_CRITICAL, "fseek(SEEK_SET) failed Fstat-file '%s'.Giving up...\n", 
			 FstatFilename);
	      return COMPUTEFSTAT_ECHECKPOINT;;
	    }
	} /* if fpFstat */

    } /* if loopcounter > 0 */
  
  while (1)
    {
#if USE_BOINC
      BOOLEAN need2checkpoint = FALSE;
#endif
      /* If our Fstat-file reaches a MaxFileSizeKB-threshold, we 'compatify'
       * the Fstat-file back to the contents of toplist with NumCandidatesToKeep
       * 
       * NOTE: if we use the toplist AND checkpointing, we must make sure that after
       * this compatification, we also checkpoint the state, 
       * because otherwise the *true* Fstat filestate will differ from the checkpointed one!
       */
      if ((uvar_NumCandidatesToKeep>0) && (fstat_bytecounter > uvar_MaxFileSizeKB * 1024) )
	{
	  INT4 howmany;
	  
	  LogPrintf ( LOG_NORMAL, "Fstat file reached MaxFileSizeKB ==> compactifying ...");

	  fclose(fpFstat);
	  howmany = atomic_write_fstat_toplist_to_file(toplist, FstatFilename, &fstat_checksum);
	  if (howmany < 0) {
	    LogPrintf (LOG_CRITICAL, "Couldn't write compacted toplist to '%s'\n", FstatFilename);
	    return (COMPUTEFSTAT_EXIT_OPENFSTAT);
	  }
	  if (howmany >= uvar_MaxFileSizeKB * 1024) {
	    LogPrintf (LOG_CRITICAL, "Size of compacted list exceeds MaxFileSize - reduce candidates to keep\n", FstatFilename);
	    return (COMPUTEFSTAT_EINPUT);
	  }
	  fstat_bytecounter = howmany;
	  
	  if ( (fpFstat = fopen(FstatFilename, "ab")) == NULL )
	    {
	      LogPrintf (LOG_CRITICAL, "Couldn't open compacted toplist for appending\n");
	      return (COMPUTEFSTAT_EXIT_OPENFSTAT2);
	    }
	  if ( fstatbuff )
	    setvbuf(fpFstat, fstatbuff, _IOFBF, uvar_OutputBufferKB * 1024);
#if USE_BOINC
	  need2checkpoint = TRUE;
#endif
	  LogPrintfVerbatim ( LOG_NORMAL, " done.\n");

	} /* if maxFileSizeKB atteined => re-compactify output file by toplist */


      /* flush fstats-file and write checkpoint-file */
#ifdef CLUSTERED_OUTPUT
      if ( uvar_doCheckpointing && fpClusters )
#else
      if ( uvar_doCheckpointing && fpFstat )
#endif
        {
#if USE_BOINC
          /* make sure the last checkpoint is written even if is not time_to_checkpoint */
	  if ( boinc_is_standalone() || need2checkpoint || boinc_time_to_checkpoint() || (loopcounter >= (thisScan.numSkyGridPoints-1)) )
            {
#endif
              FILE *fp;
#ifdef CLUSTERED_OUTPUT
              fflush (fpClusters);
#else
              fflush (fpFstat);
#endif
              if ( (fp = fopen(ckp_fname, "wb")) == NULL) {
                LogPrintf (LOG_CRITICAL, "Failed to open checkpoint-file '%s' for writing. Exiting.\n", 
			   ckp_fname);
                return COMPUTEFSTAT_ECHECKPOINT;
              }
              if ( fprintf (fp, "%" LAL_UINT4_FORMAT " %" LAL_UINT4_FORMAT " %ld\nDONE\n", 
			    loopcounter, fstat_checksum, fstat_bytecounter) < 0) {
                LogPrintf (LOG_CRITICAL, "Error writing to checkpoint-file '%s'. Exiting.\n", ckp_fname);
                return COMPUTEFSTAT_ECHECKPOINT;
              }
              fclose (fp);
	      LogPrintf (LOG_DETAIL, "Checkpointed state (%d/%d/%ld).\n", 
			 loopcounter, fstat_checksum, fstat_bytecounter );


#if USE_BOINC
              boinc_checkpoint_completed();
	      
	      need2checkpoint = FALSE;

            } /* if boinc_time_to_checkpoint() */
#endif
        } /* if doCheckpointing && fpstat */
      
      
      /* Show some progress */
#if USE_BOINC
      {
        /* update progress, the last % is reserved for polka in polka commandline runs */
        double local_fraction_done;
        if (cfsRunNo == 1)
          local_fraction_done=(((double)loopcounter)/((double)thisScan.numSkyGridPoints))*0.99/2;
        else if (cfsRunNo == 2)
          local_fraction_done=(((double)loopcounter)/((double)thisScan.numSkyGridPoints))*0.99/2+0.495;
        else
          local_fraction_done=((double)loopcounter)/((double)thisScan.numSkyGridPoints);
        if (local_fraction_done<0.0)
          local_fraction_done=0.0;
        if (local_fraction_done>1.0)
          local_fraction_done=1.0;
        boinc_fraction_done(local_fraction_done);
	/* pass variable externally to graphics routines */
	if (fraction_done_hook != NULL)
	  *fraction_done_hook=local_fraction_done;

	/* pass currently estimated fpops to client: */
	if ( LALUserVarWasSet ( &uvar_WUfpops ) ) /* ignore IOPS  */
	  boinc_ops_cumulative( local_fraction_done * uvar_WUfpops, 0); 
      }
#endif

#if !USE_BOINC
      LogPrintfVerbatim (LOG_DETAIL, "\r %5.1f%% ", 
			 (100.0* loopcounter / thisScan.numSkyGridPoints));
#endif
      
      XLALNextDopplerSkyPos( &dopplerpos, &thisScan );
      
      /* Have we scanned all PulsarDopplerParamss yet? */
      if (thisScan.state == STATE_FINISHED)
	{
	  LogPrintfVerbatim(LOG_DETAIL, "\n");
	  break;
	}

      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
      thisPoint.longitude = dopplerpos.Alpha;
      thisPoint.latitude = dopplerpos.Delta;
      LAL_CALL (LALNormalizeSkyPosition(status, &thisPoint, &thisPoint), status);
      dopplerpos.Alpha = thisPoint.longitude;
      dopplerpos.Delta = thisPoint.latitude;
      
#if USE_BOINC
      /* pass current search position, for use with starsphere.C
         revision 4.6 or greater. Need to convert radians to
         degrees. */
      if (set_search_pos_hook != NULL)
	{
	  REAL4 pAlpha=(180.0*dopplerpos.Alpha / LAL_PI); 
	  REAL4 pDelta=(180.0*dopplerpos.Delta / LAL_PI);
	  set_search_pos_hook(pAlpha,pDelta);
	}
#endif
      
      LAL_CALL (CreateDemodParams(status, dopplerpos), status);
#ifdef FILE_AMCOEFFS
      PrintAMCoeffs(dopplerpos.Alpha, dopplerpos.Delta, DemodParams->amcoe);
#endif

      /* loop over spin params */
      for (spdwn=0; spdwn < GV.SpinImax; spdwn++)
        {
          DemodParams->spinDwn[0] = GV.spinRange.fkdot[1] + spdwn * thisScan.dfkdot[1];
	  
          switch ( uvar_expLALDemod )
	    {
	    case 0: 
	      LAL_CALL ( LALDemod(status, &Fstat, SFTData, DemodParams), status);
	      break;
	    case 1:
	      LAL_CALL ( TestLALDemod(status, &Fstat, SFTData, DemodParams), status);
	      break;
	    case 3:
	      LAL_CALL ( LALDemodFAST(status, &Fstat, UpSFTData, DemodParams), status);
	      break;
	    default:
	      LogPrintf (LOG_CRITICAL, "Invalid expLALDemod value %d\n", uvar_expLALDemod);
	      return COMPUTEFSTAT_EINPUT;
	    }
	  
          /*  CLUSTER-OUTPUT: This fills-in highFLines that contains the outliers of F*/
          if ( fpClusters && (GV.FreqImax > 5) ) {
	    LAL_CALL (EstimateFLines(status), status);
	  }
          
          /* output top unclustered F-statistic results or the loudest */
          if (uvar_outputFstat || uvar_outputLoudest) 
            {
              INT4 i;
	      FstatOutputEntry outputLine = empty_FstatOutputEntry;

              for(i=0;i < GV.FreqImax ;i++)
                {
		  REAL8 Fval = 2.0*medianbias*Fstat.F[i];
		  /* correct frequency back to reference-time: assume maximaly 1 spindown */
		  REAL8 freq = GV.spinRange.fkdot[0] + i * thisScan.dfkdot[0] + GV.DeltaFreqRef; 

		  if ( uvar_SignalOnly )
		    Fval += 4;

		  if ( Fval > uvar_Fthreshold )
		    {
		      /* candidate found: insert into Top-candidate list and output */
		      outputLine.Freq  = freq;
		      outputLine.f1dot = DemodParams->spinDwn[0];
		      outputLine.Alpha = dopplerpos.Alpha;
		      outputLine.Delta = dopplerpos.Delta;
		      outputLine.Fstat = Fval;

		      /* we append this candidate to the fpFstat file if either
		       * 1) we don't use a toplist of NumCandidatesToKeep at all, or if
		       * 2) we use a toplist and the candidate was inserted.
		       */
		      if ( uvar_outputFstat && 
			   ((uvar_NumCandidatesToKeep<=0)|| insert_into_fstat_toplist(toplist,outputLine)))
			{	
			  INT4 howmany = write_fstat_toplist_item_to_fp(outputLine, fpFstat, &fstat_checksum );
			  if (howmany < 0 ) 
			    {
			      fclose(fpFstat);
			      return (COMPUTEFSTAT_EXIT_WRITEFSTAT);
			    } else
			      fstat_bytecounter += howmany;
			} /* if candidate is to be written to fpFstat file */

		    } /* if Fval > Ftheshold */

		  /* keep track of  loudest candidate */
		  if ( Fval > loudest.Fstat )
		    loudest = outputLine;

                } /* for i < FreqImax */

            } /* if outputFstat || outputLoudest */

        
          /* CLUSTER-OUTPUT: This fills-in highFLines  */
          if (fpClusters)
	    {
	      if ( (highFLines != NULL) && (highFLines->Nclusters > 0) )
		{
		  maxIndex=(INT4 *)LALMalloc(highFLines->Nclusters*sizeof(INT4));
		  
		  /*  for every cluster writes the information about it in file fpClusters */
#ifdef CLUSTERED_OUTPUT
		  if (writeFLinesCS(maxIndex, dopplerpos, fpClusters, &fstat_bytecounter, &fstat_checksum))
#else
		  if (writeFLines(maxIndex, dopplerpos, fpClusters))
#endif
		    {
		      LogPrintf (LOG_CRITICAL, "Error in writeFLines()\n\n" );
		      return COMPUTEFSTAT_EXIT_WRITEFSTAT;
		    }
		  
		  if (uvar_EstimSigParam)
		    {
		      if(writeFaFb(maxIndex, dopplerpos))
			return COMPUTEFSTAT_EXIT_WRITEFAFB;
		      if (EstimateSignalParameters(maxIndex))
			return COMPUTEFSTAT_EXIT_ESTSIGPAR;
		    } /* if signal-estimation */
          
		  LALFree(maxIndex);
		} /* if highFLines found */

	      if (PrintTopValues(/* thresh */ 0.0, /* max returned */ 1, dopplerpos))
		LogPrintf (LOG_CRITICAL, "%s: trouble making files Fmax and/or Fstats\n", argv[0]);

	      /* Set the number of the clusters detected to 0 at each iteration 
		 of the sky-direction and the spin down */
	      highFLines->Nclusters=0;

	    } /* if fpClusters */
	      
        } /* For spdwn < GV.spinImax */

      loopcounter ++;           /* number of *completed* loops */

    } /*  while SkyPos */

  if ( fpClusters )
    {
      fprintf(fpClusters, "%%DONE\n");
      fclose(fpClusters);
    }

  if (fpFstat) 
    {
      fclose(fpFstat);
    
      /* final compactification */
      if ( uvar_NumCandidatesToKeep > 0 )
	{
	  if( final_write_fstat_toplist_to_file ( toplist, FstatFilename, &fstat_checksum ) < 0 ) 
	    {
	      LogPrintf (LOG_CRITICAL, "Couldn't write compacted toplist in '%s'\n", FstatFilename);
	      return (COMPUTEFSTAT_EXIT_OPENFSTAT);
	    }
	}
      /* this is our marker indicating 'finished'.  What appears in the file is:
	 %DONE
	 on the final line */
      if ( (fpFstat = fopen (FstatFilename, "ab")) == NULL ) 
	{
	  LogPrintf ( LOG_CRITICAL, "Failed to open Fstat-file '%s' for final '%%DONE' marker!\n\n", 
		      FstatFilename);
	  return (COMPUTEFSTAT_EXIT_OPENFSTAT);
	}
      else
	{
	  fprintf(fpFstat, "%%DONE\n");
	  fclose(fpFstat);
	}

    } /* if fpFstat */

  if ( toplist ) 
    free_fstat_toplist(&toplist);

#ifdef RUN_POLKA
  } /* if (!fstats_completed) */
#endif  

  /* now write loudest canidate into separate file ".loudest" */
  if ( uvar_outputLoudest )
    {
      FILE *fpLoudest;
      if ( (fpLoudest = fopen (uvar_outputLoudest, "wb")) == NULL)
	{
	  LogPrintf (LOG_CRITICAL, "Error opening file '%s' for writing..\n\n", uvar_outputLoudest);
	  return COMPUTEFSTAT_EXIT_OPENFSTAT;
	}
      fprintf (fpLoudest, "%10f %10f %10f %10g %10f\n", 
	       loudest.Freq, loudest.Alpha, loudest.Delta, loudest.f1dot, loudest.Fstat );

      fclose(fpLoudest);
    } /* write loudest candidate to file */

  LogPrintf (LOG_NORMAL, "Search finished successfully.\n");
  
#ifdef FILE_FMAX  
  fclose(fpmax);
#endif

#ifndef RUN_POLKA
  /* remove checkpoint-file */
  remove (ckp_fname);
#endif
  /* Free DopplerSkyScan-stuff (grid) */
  LAL_CALL ( FreeDopplerSkyScan(status, &thisScan), status);
  
  LAL_CALL ( Freemem(status), status);
  
  return COMPUTEFSTAT_EXIT_OK;
  
} /* main() */


/** 
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
void
initUserVars (LALStatus *status)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
  uvar_Dterms   = 16;
  uvar_FreqBand = 0.0;
  uvar_dFreq = 0;

  uvar_Alpha    = 0.0;
  uvar_Delta    = 0.0;
  uvar_AlphaBand = 0;
  uvar_DeltaBand = 0;
  uvar_dAlpha   = 0.001;
  uvar_dDelta   = 0.001;
  uvar_skyRegion = NULL;

  uvar_ephemYear = (CHAR*)LALCalloc (1, strlen(EPHEM_YEARS)+1);

#if USE_BOINC
  strcpy (uvar_ephemYear, "");            /* the default year-string UNDER BOINC is empty! */
#else
  strcpy (uvar_ephemYear, EPHEM_YEARS);
#endif

  uvar_BaseName = (CHAR*)LALCalloc (1, strlen(SFT_BNAME)+1);
  strcpy (uvar_BaseName, SFT_BNAME);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar_ephemDir = (CHAR*)LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar_ephemDir, DEFAULT_EPHEMDIR);

  uvar_SignalOnly = FALSE;
  uvar_EstimSigParam = FALSE;
 
  uvar_f1dot = 0.0;
  uvar_df1dot   = 0.0;
  uvar_f1dotBand = 0.0;
  
  uvar_Fthreshold = 10.0;
  uvar_metricType =  LAL_PMETRIC_NONE;
  uvar_gridType = GRID_FLAT;

  uvar_metricMismatch = 0.02;

  uvar_help = FALSE;
  uvar_outputLabel = NULL;

  uvar_outputFstat = NULL;
  uvar_outputLoudest = NULL;

  uvar_skyGridFile = NULL;

  uvar_workingDir = (CHAR*)LALMalloc(512);
  strcpy(uvar_workingDir, ".");

  uvar_projectMetric = TRUE;
  
  /* if user does not set start/end times, use all SFTs */
  uvar_startTime = -1.0e308;
  uvar_endTime   =  1.0e308;

  uvar_OutputBufferKB = 2048; /* to keep the previous behavior so far */

/* some BOINC-defaults differ: checkpointing is ON, and we use experimental LALDemod() */
#if USE_BOINC
  uvar_doCheckpointing = TRUE;
  uvar_expLALDemod = 1;
#ifdef CLUSTERED_OUTPUT
  uvar_outputClusters = "Fstats";
#else
  uvar_outputClusters = NULL;	/* by default: no more cluster-output */
#endif
#else
  uvar_doCheckpointing = FALSE;
  uvar_expLALDemod = 0;
#define CLUSTERED_FNAME "Fstats"	/* provide backwards-compatible default for now */
  uvar_outputClusters = LALCalloc(1, strlen(CLUSTERED_FNAME) + 1);
  strcpy ( uvar_outputClusters, CLUSTERED_FNAME );
#endif

#if BOINC_COMPRESS
  uvar_useCompression = TRUE;
#endif

  uvar_MaxFileSizeKB = 5*1024;
  uvar_NumCandidatesToKeep = 0;	/* default: 0 = don't use toplist at all*/

  uvar_RngMedWindow = 50;

  /* register all our user-variables */
  LALregBOOLUserVar(status,       help,           'h', UVAR_HELP,     "Print this message"); 
  LALregSTRINGUserVar(status,     IFO,            'I', UVAR_REQUIRED, "Detector: GEO(0),LLO(1),LHO(2),NAUTILUS(3),VIRGO(4),TAMA(5),CIT(6)");
  LALregREALUserVar(status,       Freq,           'f', UVAR_REQUIRED, "Starting search frequency in Hz");
  LALregREALUserVar(status,       FreqBand,       'b', UVAR_OPTIONAL, "Search frequency band in Hz");
  LALregREALUserVar(status,       dFreq,          'r', UVAR_OPTIONAL, "Frequency resolution in Hz (default: 1/(2*Tsft*Nsft)");
  LALregREALUserVar(status,       Alpha,          'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status,       Delta,          'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
  LALregREALUserVar(status,       AlphaBand,      'z', UVAR_OPTIONAL, "Band in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status,       DeltaBand,      'c', UVAR_OPTIONAL, "Band in delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(status,     skyRegion,      'R', UVAR_OPTIONAL, "ALTERNATIVE: specify sky-region by polygon (or use 'allsky')");
  LALregINTUserVar(status,        gridType,        0 , UVAR_OPTIONAL, "Grid: 0=flat, 1=isotropic, 2=metric, 3=sky-file, 4=SkyGridFILE+metric");
  LALregINTUserVar(status,        metricType,     'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserVar(status,       metricMismatch, 'X', UVAR_OPTIONAL, "Maximal mismatch per dimension for metric grid");
  LALregREALUserVar(status,       dAlpha,         'l', UVAR_OPTIONAL, "Resolution in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(status,       dDelta,         'g', UVAR_OPTIONAL, "Resolution in delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(status,     skyGridFile,     0,  UVAR_OPTIONAL, "Load sky-grid from this file.");
  LALregSTRINGUserVar(status,     outputSkyGrid,   0,  UVAR_OPTIONAL, "Write sky-grid into this file.");
  LALregREALUserVar(status,       f1dot,          's', UVAR_OPTIONAL, "First spindown parameter f1dot");
  LALregREALUserVar(status,       f1dotBand,      'm', UVAR_OPTIONAL, "Search-band for f1dot");
  LALregREALUserVar(status,       df1dot,         'e', UVAR_OPTIONAL, "Resolution for f1dot (default: use metric or 1/(2*T^2))");
  LALregSTRINGUserVar(status,     DataDir,        'D', UVAR_OPTIONAL, "Directory where SFT's are located");
  LALregSTRINGUserVar(status,     BaseName,       'i', UVAR_OPTIONAL, "The base name of the input  file you want to read");
  LALregSTRINGUserVar(status,     DataFiles,	   0 , UVAR_OPTIONAL, "ALTERNATIVE: path+file-pattern specifying data SFT-files");
  LALregSTRINGUserVar(status,     ephemDir,       'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(status,     ephemYear,      'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregBOOLUserVar(status,       SignalOnly,     'S', UVAR_OPTIONAL, "Signal only flag");
  LALregBOOLUserVar(status,       EstimSigParam,  'p', UVAR_OPTIONAL, "Do Signal Parameter Estimation");
  LALregREALUserVar(status,       Fthreshold,     'F', UVAR_OPTIONAL, "Output-threshold on 2F");
  LALregSTRINGUserVar(status,     outputLabel,    'o', UVAR_OPTIONAL, "Label to be appended to all output file-names");
  LALregREALUserVar(status,       startTime,       0,  UVAR_OPTIONAL, "Ignore SFTs with GPS_time <  this value. Default:");
  LALregREALUserVar(status,       endTime,         0,  UVAR_OPTIONAL, "Ignore SFTs with GPS_time >= this value. Default:");
  LALregREALUserVar(status,	  refTime,	   0,  UVAR_OPTIONAL, "SSB reference time for pulsar-parameters");

  LALregSTRINGUserVar(status,     outputFstat,	   0,  UVAR_OPTIONAL,
		      "Output-file for the (unclustered) F-statistic field over the parameter-space");
  LALregSTRINGUserVar(status,     outputClusters,	0,  UVAR_OPTIONAL,
		      "Output-file for the *clustered* F-statistic field over the parameter-space");
  LALregINTUserVar(status,        NumCandidatesToKeep,0, UVAR_OPTIONAL, "Number of Fstat 'candidates' to keep. (0 = All)");

  /* the following are 'developer'-options */
  LALregSTRINGUserVar(status,     workingDir,     'w', UVAR_DEVELOPER, "Directory to be made the working directory.");
  LALregBOOLUserVar(status,       doCheckpointing, 0,  UVAR_DEVELOPER, "Do checkpointing and resume for previously checkpointed statuse.");
  LALregINTUserVar(status,        expLALDemod,     0,  UVAR_DEVELOPER, "Type of LALDemod to use. 0=standard, 1=exp1, 2=REAL4, 3=LALDemodFAST");
  LALregINTUserVar(status,        Dterms,         't', UVAR_DEVELOPER, "Number of terms to keep in Dirichlet kernel sum");
  LALregSTRINGUserVar(status,     mergedSFTFile,  'B', UVAR_DEVELOPER, "Merged SFT's file to be used"); 
#if BOINC_COMPRESS
  LALregBOOLUserVar(status,       useCompression,  0,  UVAR_DEVELOPER, "BOINC: use compression for download/uploading data");
#endif

  LALregBOOLUserVar(status,	  projectMetric,	0,   UVAR_DEVELOPER, 
		    "Use projected metric for skygrid");
  LALregSTRINGUserVar(status,     outputLoudest,	0,  UVAR_DEVELOPER, 
		      "Output-file for the loudest F-statistic candidate in this search");

  LALregINTUserVar(status,        OutputBufferKB, 0, UVAR_DEVELOPER, 
		   "Size of the output buffer in kB");

  LALregINTUserVar(status,        MaxFileSizeKB,  0, UVAR_DEVELOPER, 
		   "Size to which the Fstat-output file can grow until re-compactified (in kB)");

  LALregINTUserVar(status,        RngMedWindow,     0, UVAR_DEVELOPER, "Window-size to use in running-median normalization of data");

  LALregREALUserVar(status,       WUfpops, 	0, UVAR_DEVELOPER, "Einstein@Home: estimated fpops for this WU");

  DETATCHSTATUSPTR (status);
  RETURN (status);
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
  CHAR Paramfilename[MAXFILENAMELENGTH];
  

#ifdef DEBG_ESTSIGPAR
  REAL8 Ftest;
  REAL8 A2test,A3test,A4test;
#endif


  strcpy(Paramfilename,"ParamMLE");
  if (uvar_outputLabel)
    strcat(Paramfilename,uvar_outputLabel);
  
  if(!(fpMLEParam=fopen(Paramfilename,"wb"))) {
    LogPrintf (LOG_CRITICAL, "Error in EstimateSignalParameters: unable to open the file");
    return 1;
  }

  norm=2.0*sqrt(GV.tsft)/(GV.tsft*GV.SFTno);


  for(jrec=0;jrec < highFLines->Nclusters ;jrec++)
    {

      irec=maxIndex[jrec];

      A1 =  2.0*( B * creal(Fstat.Fa[irec]) - C * creal(Fstat.Fb[irec])) / D;
      A2 =  2.0*( A * creal(Fstat.Fb[irec]) - C * creal(Fstat.Fa[irec])) / D;
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
        XLALPrintError ( 
                "Imaginary Cos[iota]; cannot compute parameters\n"
                "in the EstimateSignalParameters routine\n"
                "in ComputeFStatistic code\n"
                "Now exiting...\n");
        /*        break; */
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
          XLALPrintError ( "Imaginary beta; cannot compute parameters");
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


      XLALPrintError ( "LALDemod_Estimate output: "
              "A1=%20.15f A2=%20.15f A3=%20.15f A4=%20.15f\n"
              ,A1,A2,A3,A4);
      XLALPrintError ( "Reconstructed from MLE: "
              "A1=%20.15f A2=%20.15f A3=%20.15f A4=%20.15f !!!!\n\n",
              A1test,A2test,A3test,A4test);
      fflush(stderr);


      if(fabs(A1-A1test)>fabs(A1)/(10e5)){ 
        XLALPrintError ( "Something is wrong with Estimate A1\n");
        XLALPrintError ( "Frequency index %d, %lf (Hz),A1=%f,A1test=%f\n",
                irec,GV.fkdot0->data[0]+irec*DemodParams->df,A1,A1test);
        XLALPrintError ( "relative error Abs((A1-A1test)/A1)=%lf\n",
                fabs(A1-A1test)/fabs(A1));
        return 1;
      }
      if(fabs(A2-A2test)>fabs(A2)/(10e5)){ 
        XLALPrintError ( "Something is wrong with Estimate A2\n");
        XLALPrintError ( "Frequency index %d, %lf (Hz),A2=%f,A2test=%f\n",
                irec,GV.fkdot0->data[0]+irec*DemodParams->df,A2,A2test);
        XLALPrintError ( "relative error Abs((A2-A2test)/A2)=%lf\n",
                fabs(A2-A2test)/fabs(A2));
        return 1;
      }
      if(fabs(A3-A3test)>fabs(A3)/(10e5)){ 
        XLALPrintError ( "Something is wrong with Estimate A3\n");
        XLALPrintError ( "Frequency index %d, %lf (Hz),A3=%f,A3test=%f\n",
                irec,GV.fkdot0->data[0]+irec*DemodParams->df,A3,A3test);
        XLALPrintError ( "relative error Abs((A3-A3test)/A3)=%lf\n",
                fabs(A3-A3test)/fabs(A3));
        return 1;
      }
      if(fabs(A4-A4test)>fabs(A4)/(10e5)){ 
        XLALPrintError ( "Something is wrong with Estimate A4\n");
        XLALPrintError ( "Frequency index %d, %lf (Hz),A4=%f,A4test=%f\n",
                irec,GV.fkdot0->data[0]+irec*DemodParams->df,A1,A1test);
        XLALPrintError ( "relative error Abs((A4-A4test)/A4)=%lf\n",
                fabs(A4-A4test)/fabs(A4));
        return 1;
      }

      
      /* Reconstruct F. Compare it with the original value. */
      Ftest=(A*(A1*A1+A3*A3)+B*(A2*A2+A4*A4)+2.0*C*(A1*A2+A3*A4))/
        4.0*4.0/GV.SFTno;


      if(fabs(Fstat.F[irec] - Ftest)> fabs(Ftest)/10e5){ 
        XLALPrintError ( "Something is wrong with Estimate in F\n");
        XLALPrintError ( "Frequency index %d, %lf (Hz),F=%f,Ftest=%f\n",
                irec,GV.fkdot0->data[0]+irec*DemodParams->df,Fstat.F[irec],Ftest);
        XLALPrintError ( "relative error Abs((F-Ftest)/Ftest)=%lf\n",
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
      XLALPrintError ( "A=%f,B=%f,C=%f,f=%f,h0=%f,F=%f\n",
              A,B,C,GV.spinRange.fkdot[0]+irec* DemodParams->df,h0mle,Fstat.F[irec]*medianbias);
      }
#endif

      /* Note that we print out MLE of 2.0*Phi0_JKS */
      /* because Phi0_PULGROUPDOC=2.0*Phi0_JKS */
      /* and Phi0_PULGROUPDOC is the one used in In.data. */
 
      /* medianbias is 1 if GV.SignalOnly==1 */
      fprintf(fpMLEParam,"%16.8f %22E", GV.spinRange.fkdot[0] + irec*DemodParams->df, 
	      2.0*medianbias*Fstat.F[irec]);


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
int writeFaFb(INT4 *maxIndex, PulsarDopplerParams searchpos)
{
  INT4 irec,jrec;
  INT4 ind,krec=0;
  CHAR filename[MAXFILENAMELENGTH]; /* Base of the output file name */
  CHAR noiseswitch[16];
  CHAR clusterno[16];
  INT4 N;
  FILE * fp=NULL;
  REAL8 bias=1.0;
  CHAR FaFbfilename[MAXFILENAMELENGTH];

  if ( lalDebugLevel > 10)	/* dummy: avoid warnings */
    XLALPrintError (  "%f, %f", searchpos.Alpha, searchpos.Delta);
  
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
      LogPrintf (LOG_CRITICAL,  "Unable to open a file %s\n",filename);
      return 1;
    }

    N=highFLines->NclustPoints[irec];


    /* The header contains */
    /* N the number of points in the cluster */
    /* the Frequency where the maximum amplitude is */
    /* A,B,C coefficients */
    ind=highFLines->Iclust[krec];

    fprintf(fp,"%10d\n",N);
    fprintf(fp,"%22.12f %22.12f\n",
            GV.spinRange.fkdot[0] + maxIndex[irec]* DemodParams->df,
            Fstat.F[maxIndex[irec]]*bias*bias);
    fprintf(fp,"%22.12f %22.12f\n",GV.spinRange.fkdot[0] + ind * DemodParams->df, DemodParams->df);
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
      ind=highFLines->Iclust[krec];
      krec++;

#ifdef DEBG_FAFB
      fprintf(fp,"%22.16f %22.16f "
                 "%E %20.17f %20.17f "
                 "%22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n",
              GV.spinRange.fkdot[0] + ind* DemodParams->df,Fstat.F[ind]*bias*bias,
              DemodParams->spinDwn[0], searchpos.Alpha, searchpos.Delta,
              Fstat.Fa[ind].re/sqrt(GV.SFTno)*bias,
              Fstat.Fa[ind].im/sqrt(GV.SFTno)*bias,
              Fstat.Fb[ind].re/sqrt(GV.SFTno)*bias,
              Fstat.Fb[ind].im/sqrt(GV.SFTno)*bias,
              amc.A,amc.B,amc.C);
#else
      /* Freqency, Re[Fa],Im[Fa],Re[Fb],Im[Fb], F */
      fprintf(fp,"%22.16f %22.12f %22.12f %22.12f %22.12f %22.12f\n",
              GV.spinRange.fkdot[0] + ind* DemodParams->df,
              creal(Fstat.Fa[ind])/sqrt(GV.SFTno)*bias,
              Fstat.Fa[ind].im/sqrt(GV.SFTno)*bias,
              creal(Fstat.Fb[ind])/sqrt(GV.SFTno)*bias,
              Fstat.Fb[ind].im/sqrt(GV.SFTno)*bias,
              Fstat.F[ind]*bias*bias);
#endif


    }
    fclose(fp);
  }

  return 0;

} /* writeFaFb() */



/*******************************************************************************/

void CreateDemodParams (LALStatus *status, PulsarDopplerParams searchpos)
{
  CSParams *csParams  = NULL;        /* ComputeSky parameters */
  AMCoeffsParams *amParams;
  EarthState earth;
  EmissionTime emit;
  LIGOTimeGPS *midTS=NULL;           /* Time stamps for amplitude modulation coefficients */
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  INT4 k;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.site.location[0] = GV.Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1] = GV.Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2] = GV.Detector.location[2]/LAL_C_SI;
  baryinput.alpha = searchpos.Alpha;
  baryinput.delta = searchpos.Delta;
  baryinput.dInv = 0.e0;

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
  amParams->das->pSource->equatorialCoords.latitude = searchpos.Delta;
  amParams->das->pSource->equatorialCoords.longitude = searchpos.Alpha;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!*/

 /* Mid point of each SFT */
   midTS = (LIGOTimeGPS *)LALCalloc(GV.SFTno, sizeof(LIGOTimeGPS));
   for(k=0; k<GV.SFTno; k++)
     { 
       /* FIXME: loss of precision; consider
       midTS[k] = timestamps[k];
       XLALGPSAdd(&midTS[k], 0.5*GV.tsft);
       */
       REAL8 teemp=0.0;
       teemp = XLALGPSGetREAL8(&(timestamps[k]));
       teemp += 0.5*GV.tsft;
       XLALGPSSetREAL8(&(midTS[k]), teemp);
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
  csParams->skyPos[0] = searchpos.Alpha;
  csParams->skyPos[1] = searchpos.Delta;
  csParams->earth = &earth;
  csParams->emit = &emit;

/* Finally, DemodParams */
  DemodParams->amcoe = &amc;
  DemodParams->spinDwnOrder = 1;
  DemodParams->SFTno = GV.SFTno;

  DemodParams->f0   = GV.spinRange.fkdot[0];
  DemodParams->imax = GV.FreqImax;

  DemodParams->Dterms = uvar_Dterms;
  DemodParams->ifmin = GV.ifmin;

  DemodParams->returnFaFb = uvar_EstimSigParam;

  /* compute the "sky-constants" A and B */
  TRY ( LALComputeSky(status->statusPtr, DemodParams->skyConst, 0, csParams), status);  
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
int 
writeFLines(INT4 *maxIndex, PulsarDopplerParams searchpos, FILE *fpOut)
{
  INT4 i,j,j_1,j_2,k,N;
  REAL8 max,log2val,mean,var,std,R;
  REAL8 freq;
  INT4 imax;

  log2val=medianbias;
 
  j_1=0;
  j_2=0;

  for (i=0;i<highFLines->Nclusters;i++){
    N=highFLines->NclustPoints[i];
    
    /*  determine maximum of j-th cluster */
    /*  and compute mean */
    max=0.0;
    imax=0;
    mean=0.0;
    std=0.0;
    for (j=0;j<N;j++){
      R=2.0*log2val*highFLines->clusters[j_1];
      k=highFLines->Iclust[j_1];
      j_1=j_1+1;
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
      R=2.0*log2val*highFLines->clusters[j_2];
      j_2=j_2+1;
      var=var+(R-mean)*(R-mean);
    }/*  end j loop over points of i-th cluster  */
    var=var/N;
    std=sqrt(var);
    /* correct frequency back to reference-time: assume maximally 1 spindown */
    freq = GV.spinRange.fkdot[0] + imax * DemodParams->df + GV.DeltaFreqRef;

    /* print the output */
    if (fpOut) 
      fprintf( fpOut, "%16.12f %10.8f %10.8f    %d %10.5f %10.5f %20.17f\n",
	       freq, searchpos.Alpha, searchpos.Delta, N, mean, std, max );
    
  } /*  end i loop over different clusters */
  
  return 0;

} /* writeFLines() */


/* checksumming version of WriteFLines for checkpointing the clustered output */
int 
writeFLinesCS(INT4 *maxIndex, PulsarDopplerParams searchpos, FILE *fpOut, long*bytecount, UINT4*checksum)
{
  INT4 i,j,j_1,j_2,k,N;
  REAL8 max,log2val,mean,var,std,R;
  REAL8 freq;
  INT4 imax;
  CHAR buf[256];
  INT4 len,chr;

  log2val=medianbias;
 
  j_1=0;
  j_2=0;

  for (i=0;i<highFLines->Nclusters;i++){
    N=highFLines->NclustPoints[i];
    
    /*  determine maximum of j-th cluster */
    /*  and compute mean */
    max=0.0;
    imax=0;
    mean=0.0;
    std=0.0;
    for (j=0;j<N;j++){
      R=2.0*log2val*highFLines->clusters[j_1];
      k=highFLines->Iclust[j_1];
      j_1=j_1+1;
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
      R=2.0*log2val*highFLines->clusters[j_2];
      j_2=j_2+1;
      var=var+(R-mean)*(R-mean);
    }/*  end j loop over points of i-th cluster  */
    var=var/N;
    std=sqrt(var);
    /* correct frequency back to reference-time: assume maximally 1 spindown */
    freq = GV.spinRange.fkdot[0] + imax * DemodParams->df + GV.DeltaFreqRef;

    /* print the output */
    if (fpOut) {
	len =
	    snprintf( buf, sizeof(buf), "%16.12f %10.8f %10.8f    %d %10.5f %10.5f %20.17f\n",
		    freq, searchpos.Alpha, searchpos.Delta, N, mean, std, max );
	if ((UINT4)len > sizeof(buf))
	    return(-1);
	*bytecount += len;
	for(chr=0;chr<len;chr++)
	    *checksum += buf[chr];
	if (fprintf(fpOut,"%s",buf) <0)
	    return(-1);
    }
    
  } /*  end i loop over different clusters */
  
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
  INT4 filenum=0,offset;
  FILE *fp=NULL;
  size_t errorcode;
  UINT4 ndeltaf;
  INT4 k=0;

  SFTData=(FFT **)LALMalloc(GV.SFTno*sizeof(FFT *));
  timestamps=(LIGOTimeGPS *)LALMalloc(GV.SFTno*sizeof(LIGOTimeGPS));
  
  /* if using a merged SFT file, it's already open */
  if (uvar_mergedSFTFile)
    fp=fp_mergedSFT;

  for (filenum=0;filenum<GV.SFTno; /* INCREMENT IN LOOP */ ) 
    {
      REAL8 thisSFTtime;

      /* seek to filenum'th SFT k bytes in */
      if (uvar_mergedSFTFile){
        if (fseek(fp,k,SEEK_SET)) {
          LogPrintf (LOG_CRITICAL, "Unable to seek to the start of %s !\n",GV.filelist[filenum]);
          return 1;
        }
      }
      else {
        if (!(fp=fopen(GV.filelist[filenum],"rb"))) {
          LogPrintf (LOG_CRITICAL, "Weird... %s doesn't exist!\n", GV.filelist[filenum]);
          return 1;
        }
      }

      /* Read in the header from the file */
      errorcode=fread((void*)&header,sizeof(header),1,fp);
      if (errorcode!=1) 
        {
          LogPrintf( LOG_CRITICAL, "No header in data file %s\n", GV.filelist[filenum]);
          return 1;
        }

      if (reverse_endian)
        swapheader(&header);
      
      if (header.endian!=1.0)
        {
          LogPrintf (LOG_CRITICAL,  "First object in file %s is not (double)1.0!\n",GV.filelist[filenum]);
          LogPrintf (LOG_CRITICAL,  "The file might be corrupted\n\n");
          return 2;
        }
    
      /* Check that the time base is positive */
      if (header.tbase<=0.0)
        {
          LogPrintf ( LOG_CRITICAL, "Timebase %f from data file %s non-positive!\n",
                  header.tbase,GV.filelist[filenum]);
          return 3;
        }
        
      /* Check that are Frequency bins needed are in data set */
      if (GV.ifmin<header.firstfreqindex || 
          GV.ifmax>header.firstfreqindex+header.nsamples) 
        {
          LogPrintf (LOG_CRITICAL, "Freq index range %d->%d not in %d to %d (file %s)\n",
		     GV.ifmin,GV.ifmax,header.firstfreqindex,
		     header.firstfreqindex+header.nsamples,GV.filelist[filenum]);
          return 4;
        }

      /* Move forward in file */
      offset=(GV.ifmin-header.firstfreqindex)*2*sizeof(REAL4);
      errorcode=fseek(fp,offset,SEEK_CUR);
      if (errorcode) 
        {
          perror(GV.filelist[filenum]);
          LogPrintf (LOG_CRITICAL, "Can't get to offset %d in file %s\n",offset,GV.filelist[filenum]);
          return 5;
        }


      /* determine if THIS SFT is in the range of times of those which
       * we need to use.  If so, read the data into arrays, else ignore it.
       */

      thisSFTtime=(REAL8)header.gps_sec+(1.e-9)*(REAL8)header.gps_nsec;
      if (uvar_startTime<=thisSFTtime && thisSFTtime<uvar_endTime) {

      /* Make data structures and put time stamps from file into array */
      timestamps[filenum].gpsSeconds = header.gps_sec;
      timestamps[filenum].gpsNanoSeconds = header.gps_nsec;
      ndeltaf=GV.ifmax-GV.ifmin+1;
      SFTData[filenum]=(FFT *)LALMalloc(sizeof(FFT));
      SFTData[filenum]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
      SFTData[filenum]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
      SFTData[filenum]->fft->data->data=(COMPLEX8 *)LALMalloc(ndeltaf*sizeof(COMPLEX8));

      /* Fill in actual SFT data, and housekeeping */
      errorcode=fread((void*)(SFTData[filenum]->fft->data->data), sizeof(COMPLEX8), ndeltaf, fp);
      if (errorcode!=ndeltaf){
        perror(GV.filelist[filenum]);
        LogPrintf (LOG_CRITICAL,  "The SFT data was truncated.  Only read %d not %d complex floats\n", 
		(int)errorcode, ndeltaf);
        return 6;
      }
      /* reverse byte order if needed */
      if (reverse_endian) {
        unsigned int cnt;
        for (cnt=0; cnt<ndeltaf; cnt++) {
          swap4((char *)&(crealf(SFTData[filenum]->fft->data->data[cnt])));
          swap4((char *)&(cimagf(SFTData[filenum]->fft->data->data[cnt])));
        }
      }
      
      SFTData[filenum]->fft->epoch=timestamps[filenum];
      SFTData[filenum]->fft->f0 = GV.ifmin / GV.tsft;
      SFTData[filenum]->fft->deltaF = 1.0 / GV.tsft;
      SFTData[filenum]->fft->data->length = ndeltaf;
      filenum++;
      }

      if (uvar_mergedSFTFile)
        k+=sizeof(header)+header.nsamples*sizeof(COMPLEX8);
      else
        fclose(fp);     /* close file */
      
    }
  return 0;  

} /* ReadSFTData() */

/** Upsamples the SFT data
 *
 * This function upsamples the SFTs using the Dirichlet kernel
 *
 */
int UpsampleSFTData(void)
{
  INT4 filenum=0;
  INT4 k=0,i;
  /* sorry - variable array size is a gcc extension.
     As this is not changed anywhere in this function, I put this as a macro.
     BM */
#define res 64                    /*resolution of the argument of the trig functions: 2pi/res.*/
  REAL8 sinVal[res+1], cosVal[res+1];
  INT4 ndeltaf=GV.ifmax-GV.ifmin+1, imax;
  REAL8 df=DemodParams->df; /* frequency resolution of upsampled SFTs */

  /* number of fine ferquency bins in upsampled SFTs */
  imax = (INT4) ( ((REAL4)(GV.ifmax - GV.ifmin))/(REAL4)GV.tsft /df );
  
  /* This size LUT gives errors ~ 10^-7 with a three-term Taylor series */
   for (k=0; k<=res; k++){
     sinVal[k]=sin((LAL_TWOPI*k)/res);
     cosVal[k]=cos((LAL_TWOPI*k)/res);
   }

  UpSFTData=(FFT **)LALMalloc(GV.SFTno*sizeof(FFT *));

  for (filenum=0;filenum<GV.SFTno; filenum++ ) 
    {
  
      UpSFTData[filenum]=(FFT *)LALMalloc(sizeof(FFT));

      UpSFTData[filenum]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
      UpSFTData[filenum]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
      UpSFTData[filenum]->fft->data->data=(COMPLEX8 *)LALMalloc(imax*sizeof(COMPLEX8));

      /* Loop over frequencies to be upsampled to */
      for(i=0 ; i <  imax ; i++ )
	{
	  REAL8 f =  GV.ifmin/GV.tsft+ i*df;
	  REAL8 xTemp, realXP, imagXP, tempFreq, tsin, tcos;
	  INT4 ind, k1;

	  xTemp=f*GV.tsft;
	  realXP=0.0;
	  imagXP=0.0;
	  /* find correct index into LUT -- pick closest point */
	  tempFreq=xTemp-(INT4)xTemp;
	  ind=(INT4)(tempFreq*res+0.5); 

	  {
	    REAL8 d=LAL_TWOPI*(tempFreq-(REAL8)ind/64.0);/*just like res above*/
	    REAL8 d2=0.5*d*d;
	    REAL8 ts=sinVal[ind];
	    REAL8 tc=cosVal[ind];
	    
	    tsin=ts+d*tc-d2*ts;
	    tcos=tc-d*ts-d2*tc-1.0;
	  }

	  tempFreq=LAL_TWOPI*( tempFreq + uvar_Dterms-1);
	  k1=(INT4)xTemp-uvar_Dterms+1;

	  /* Upsample each SFT using the Dirichlet kernel */
	  for( k = 0; k < 2*uvar_Dterms; k++)
	    {
	      COMPLEX8 Xalpha_k;
	      REAL8 x, realP, imagP;
	      INT4 SFTIndex;

	      x=tempFreq-LAL_TWOPI*(REAL8)k;
	      realP=tsin/x;
	      imagP=tcos/x;
	      
	      /* If x is small we need correct x->0 limit of Dirichlet kernel */
	      if(fabs(x) < 0.000001)
		{
		  realP=1.0;
		  imagP=0.0;
		}
	      
	      SFTIndex=k1+k-GV.ifmin;

/* 	      fprintf(stdout,"%d %d %d %d %d\n",k1,k,GV.ifmin,k1+k-GV.ifmin,SFTIndex); */

	      if (SFTIndex < 0 || SFTIndex > ndeltaf-1)
		{
		  Xalpha_k.realf_FIXME= Xalpha_k.imagf_FIXME=0.0;
		}else{
		Xalpha_k=SFTData[filenum]->fft->data->data[SFTIndex];
	      }

	      /* these four lines compute P*xtilde */
	      realXP += crealf(Xalpha_k)*realP;
	      realXP -= cimagf(Xalpha_k)*imagP;
	      imagXP += crealf(Xalpha_k)*imagP;
	      imagXP += cimagf(Xalpha_k)*realP;
	    }
     
	  /* fill in the data here */
	  UpSFTData[filenum]->fft->data->data[i].realf_FIXME= realXP;
	  UpSFTData[filenum]->fft->data->data[i].imagf_FIXME= imagXP;
	  
/*  	  fprintf(stdout,"%d %d %e %e\n", i, filenum, realXP, imagXP); */

	}

      /* Fill in actual SFT data, and housekeeping */
      UpSFTData[filenum]->fft->epoch=timestamps[filenum];
      UpSFTData[filenum]->fft->f0 = GV.ifmin / GV.tsft;
      /* The line below is a total hack. I need to have this be the wrong resolution to pass it to 
	 LALDemodFAST; should move this to DemodParams */
      UpSFTData[filenum]->fft->deltaF = 1.0 / GV.tsft;
      UpSFTData[filenum]->fft->data->length = imax;

    }

/*   for (filenum=0;filenum<GV.SFTno; filenum++ )  */
/*     { */
/*       for(i=0 ; i <  imax ; i++ ) */
/* 	{ */
/* 	 fprintf(stdout,"%d %d %e %e\n", i, filenum, */
/* 		 UpSFTData[filenum]->fft->data->data[i].re, */
/* 		 UpSFTData[filenum]->fft->data->data[i].im); */
/* 	} */
/*     } */


  return 0;  

} /* UpsampleSFTData() */


/*----------------------------------------------------------------------*/
/** Do some basic initializations of the F-statistic code before starting the main-loop.
 * Things we do in this function: 
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
  INT4 filenum=0;   
#ifdef HAVE_GLOB_H
  glob_t globbuf;
#endif
  LIGOTimeGPS starttime = {0,0};
  INT4 last_time_used=0;
  REAL8 thisSFTtime;
  INT4 blockno=0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set the current working directory */
  if(chdir(uvar_workingDir) != 0)
    {
      LogPrintf (LOG_CRITICAL,  "in Main: unable to change directory to %s\n", uvar_workingDir);
      ABORT (status, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
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
    char zippedname[MAXFILENAMELENGTH];
    int boinczipret;
    
    /* see if there is a softlink to earth.zip */
    if (!boinc_resolve_filename("earth.zip", zippedname,sizeof(zippedname))) {
      /* if there is, unzip it into the current directory */
      if ((boinczipret=boinc_zip(UNZIP_IT, zippedname, "./"))) {
        LogPrintf (LOG_CRITICAL,  "Error in unzipping file %s to earth.  Return value %d\n", 
		   zippedname, boinczipret);
        boinc_finish(COMPUTEFSTAT_EXIT_CANTUNZIP);
      }
    }
    /* see if there is a softlink to sun.zip */
    if (!boinc_resolve_filename("sun.zip", zippedname, sizeof(zippedname))) {
      /* if there is, unzip it into the current directory */
      if ((boinczipret=boinc_zip(UNZIP_IT, zippedname, "./"))) {
        LogPrintf (LOG_CRITICAL,  "Error in unzipping file %s to sun.  Return value %d\n", 
		zippedname, boinczipret);
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
      LogPrintf (LOG_CRITICAL, "Unable to open SFT file %s\n", uvar_mergedSFTFile);
      ABORT (status, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
    }

    filenum = 0;
    while (fread((void*)&header,sizeof(header),1,fp) == 1) {
      char tmp[MAXFILENAMELENGTH];

      /* prepare memory for another filename */
      if ( (cfg->filelist=(CHAR**)LALRealloc(cfg->filelist, (filenum+1)*sizeof(CHAR*)))==NULL) {
        ABORT (status, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
      }
      /* store a "file name" composed of merged name + block number */
      sprintf(tmp, "%s (block %d)", uvar_mergedSFTFile, ++blockno);
      if ( (cfg->filelist[filenum] = (CHAR*)LALCalloc (1, strlen(tmp)+1)) == NULL) {
        ABORT (status, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
      }
      strcpy(cfg->filelist[filenum],tmp);
      
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
        LogPrintf (LOG_CRITICAL, "First object in file %s is not (double)1.0!\n",cfg->filelist[filenum]);
        LogPrintf (LOG_CRITICAL, "The file might be corrupted\n\n");
        ABORT (status, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
      }
      
      /* if this is the first SFT, save initial time */
      if (filenum==0) {
        cfg->Ti=header.gps_sec;
        starttime.gpsSeconds = header.gps_sec;
        starttime.gpsNanoSeconds = header.gps_nsec;
      }
      /* increment file no and pointer to the start of the next header */
      thisSFTtime=(REAL8)header.gps_sec+(1.e-9)*(REAL8)header.gps_nsec;
      if (uvar_startTime<=thisSFTtime && thisSFTtime<uvar_endTime) {
        filenum++;
        last_time_used=header.gps_sec;
      }
      else
        LALFree(cfg->filelist[filenum]);

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
  else if ( uvar_DataDir || uvar_DataFiles )
    {
      if ( uvar_DataDir )
	{
	  strcpy(command, uvar_DataDir);
	  strcat(command,"/*");
	  strcat(command, uvar_BaseName);
	  strcat(command,"*");
	}
      else
	{
	  strcpy (command, uvar_DataFiles);
	}
    
#ifdef HAVE_GLOB_H
    globbuf.gl_offs = 1;
    glob(command, GLOB_ERR, NULL, &globbuf);
    
    /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */
    
    if(globbuf.gl_pathc==0)
      {
        LogPrintf (LOG_CRITICAL, "No SFTs in directory %s ... Exiting.\n\n", uvar_DataDir);
        ABORT (status, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
      }
    /* prepare memory for all filenames */
    if ( (cfg->filelist = (CHAR**)LALCalloc(globbuf.gl_pathc, sizeof(CHAR*))) == NULL) {
      ABORT (status, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
    }
    while ((UINT4)filenum < (UINT4)globbuf.gl_pathc) 
      {
        if ((cfg->filelist[filenum]=(CHAR*)LALCalloc(1,strlen(globbuf.gl_pathv[filenum])+1))== NULL)
	  {
	    ABORT (status, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
	  }
        strcpy(cfg->filelist[filenum],globbuf.gl_pathv[filenum]);
        filenum++;
      }
    globfree(&globbuf);
#endif
  }
  cfg->SFTno = filenum; /* remember this is 1 more than the index value */

  /* check that we found any suitable SFTs at all!! */
  if ( cfg->SFTno  == 0 )
    {
      LogPrintf (LOG_CRITICAL, "No suitable SFTs found within given time-range [%f, %g]\n\n", 
		 uvar_startTime, uvar_endTime);
      ABORT (status, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
    }


  if (!uvar_mergedSFTFile) {
    /* open FIRST file and get info from it*/
    fp=fopen(cfg->filelist[0],"rb");
    /* read in the header from the file */
    errorcode=fread((void*)&header,sizeof(header),1,fp);

    if (errorcode!=1) 
      {
        LogPrintf (LOG_CRITICAL, "No header in data file %s\n\n", cfg->filelist[0]);
        ABORT (status, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
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
        LogPrintf (LOG_CRITICAL, "First object in file %s is not (double)1.0!\n",cfg->filelist[0]);
        LogPrintf (LOG_CRITICAL, "The file might be corrupted\n\n");
        ABORT (status, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
      }
    fclose(fp);
    
    /* INITIAL TIME */
    starttime.gpsSeconds = header.gps_sec;
    starttime.gpsNanoSeconds = header.gps_nsec;
    cfg->Ti = header.gps_sec; 
    
    /* open LAST file and get info from it*/
    fp=fopen(cfg->filelist[filenum-1],"rb");
    /* read in the header from the file */
    errorcode=fread((void*)&header,sizeof(header),1,fp);
    if (errorcode!=1) 
      {
        LogPrintf (LOG_CRITICAL, "No header in data file %s\n", cfg->filelist[filenum-1]);
        ABORT (status, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
      }
    if (reverse_endian)
      swapheader(&header);

    /* check that data is correct endian order */
    if (header.endian!=1.0)
      {
        LogPrintf (LOG_CRITICAL, "First object in file %s is not (double)1.0!\n",cfg->filelist[filenum-1]);
        LogPrintf (LOG_CRITICAL, "The file might be corrupted\n\n");
        ABORT (status, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
      }
    fclose(fp);

    /* FINAL TIME */
    cfg->Tf = (INT4) (header.gps_sec+header.tbase);  
    /* Time baseline of SFTs */
    cfg->tsft=header.tbase;  
    cfg->nsamples=header.nsamples;    /* # of Freq. bins in SFT*/
  }


  /* ----- determine search-region in parameter-space from user-input ----- */
  {
    LIGOTimeGPS refTime;
    PulsarSpins fkdotRef, fkdot0;

    BOOLEAN haveAlphaDelta = LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta);

    /* Standardise reference-time:
     * translate spindown-paramters {f, fdot, fdotdot..} from the user-specified 
     * reference-time uvar_refTime to the internal reference-time, which 
     * we chose as the start-time of the first SFT (*verbatim*, i.e. not translated to SSB! )
     */
    if ( LALUserVarWasSet(&uvar_refTime)) {
      XLALGPSSetREAL8(&refTime, uvar_refTime);
    } else {
      refTime.gpsSeconds = cfg->Ti; refTime.gpsNanoSeconds = 0;
    }
    
    memset ( &fkdotRef, 0, sizeof(fkdotRef) );
    memset ( &fkdot0, 0, sizeof(fkdotRef) );
    fkdotRef[0] = uvar_Freq;
    fkdotRef[1] = uvar_f1dot;	    /* currently not more spindowns implemented... */
    
    /* now translate spin-params to internal reference-time (ie. startTime) */
    TRY ( LALExtrapolatePulsarSpins ( status->statusPtr, fkdot0, starttime, fkdotRef, refTime), status );

    /* we assume for now that we only have 1 spindown, which is therefore
     * constant over time, and we only need to correct frequencies by a constant
     * accounting for the difference between reference-time and data start-time,
     * i.e. DeltaFreqRef =  f(tRef) - f(tStart)
     */
    cfg->DeltaFreqRef = fkdotRef[0] - fkdot0[0];


    cfg->spinRange.fkdot[0] = fkdot0[0];
    cfg->spinRange.fkdot[1] = fkdot0[1];

    cfg->spinRange.fkdotBand[0] = uvar_FreqBand;
    cfg->spinRange.fkdotBand[1] = uvar_f1dotBand;
    /* 'normalize' if negative bands where given */
    if ( cfg->spinRange.fkdotBand[0] < 0 )
      {
	cfg->spinRange.fkdotBand[0] *= -1.0;
	cfg->spinRange.fkdot[0] -= cfg->spinRange.fkdotBand[0];
      }
    if ( cfg->spinRange.fkdotBand[1] < 0 )
      {
	cfg->spinRange.fkdotBand[1] *= -1.0;
	cfg->spinRange.fkdot[1] -= cfg->spinRange.fkdotBand[1];
      }

    /* get sky-region */
    if (uvar_skyRegion)
      {
	cfg->skyRegionString = (CHAR*)LALCalloc(1, strlen(uvar_skyRegion)+1);
	if ( cfg->skyRegionString == NULL ) {
	  ABORT (status, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
	}
	strcpy (cfg->skyRegionString, uvar_skyRegion);
      }
    else if (haveAlphaDelta)    /* parse this into a sky-region */
      {
	REAL8 eps = 1e-9;	/* hack for backwards compatbility */
	TRY ( SkySquare2String( status->statusPtr, &(cfg->skyRegionString),
				uvar_Alpha, uvar_Delta,
				uvar_AlphaBand + eps, uvar_DeltaBand + eps), status);
      }
  } /* find search-region */


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
      LogPrintf (LOG_CRITICAL, "Unknown detector. Currently allowed are 'GEO', 'LLO', 'LHO',"
		     " 'NAUTILUS', 'VIRGO', 'TAMA', 'CIT' or '0'-'6'\n\n");
      ABORT (status, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
    }

  /* ----- handle skygrid-file ----- */
  cfg->skyGridFile = NULL;
  if (uvar_skyGridFile)  
    {
      UINT4 len = strlen (uvar_skyGridFile);

      /* +3 to allow prepending of './' if no path was given */
      if ( ( cfg->skyGridFile = LALCalloc (1, len + 3)) == NULL ) {
	ABORT (status, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
      }

      /* does it contain a path? */
      if ( strchr(uvar_skyGridFile, '/') == NULL )
	sprintf ( cfg->skyGridFile, "./%s", uvar_skyGridFile );
      else
	strcpy ( cfg->skyGridFile, uvar_skyGridFile );
    }

#if USE_BOINC
  if ( uvar_skyGridFile )
      {
	char resolved_name[MAXFILENAMELENGTH];
	int is_zipped;
	/* get physical file-name for sky-grid file */
	LogPrintf (LOG_DEBUG, "boinc-resolving sky-grid file '%s' ... ", cfg->skyGridFile );
	if ( boinc_resolve_filename ( cfg->skyGridFile, resolved_name, sizeof(resolved_name)) )
	  {
	    LogPrintf (LOG_CRITICAL,  "Can't resolve file '%s'\n", cfg->skyGridFile );
	    ABORT (status, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
	  }
	LogPrintfVerbatim ( LOG_DEBUG, "resolved to '%s'.\n", resolved_name );
	
	is_zipped = is_zip_file ( resolved_name );
	if ( is_zipped == -1 ) {
	  LogPrintf (LOG_CRITICAL, "Error determining whether '%s' is a zip-file.\n", resolved_name );
	  ABORT (status, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
	}
	
	/* if sky-Grid is not zipped ==> use the resolved name */
	if ( !is_zipped )
	  {
	    LogPrintf ( LOG_DEBUG, "Skygrid file is not zip-compressed. Using boinc-resolved file.\n");
	    LALFree ( cfg->skyGridFile );
	    if ( (cfg->skyGridFile = LALCalloc ( 1, strlen ( resolved_name ) + 1 )) == NULL ) {
	      ABORT (status, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
	    }
	    strcpy ( cfg->skyGridFile, resolved_name );	/* done */
	  }
	/* ELSE, if physical skyGrid-file is zipped ==> unzip it *locally*:
	 * NOTE: a-priori we don't know the output-file(s) of the unzipping-process.
	 * we rely here on the unzipped skygrid-files ALWAYS having the same name
	 * as the original skyGridFile (i.e. BEFORE boinc-resolving!)
	 **/
	else 	
	  {
	    int ret;
	    char *outdir = NULL;
	    char *ptr;
	    char *zipfile;
	    char *tmpfile = "./__tmp.zip";

	    LogPrintf (LOG_DEBUG, "Sky-grid file is zip-compressed.\n");

	    /* unzip the resolved-name locally: first make sure to remove the original skyGridFile,
	     * as the unzipped file is expected to replace it
	     */
	    if ( strcmp ( cfg->skyGridFile, resolved_name )  ) /* only if different files! */
	      {
		LogPrintf (LOG_DEBUG, "Removing boinc-link file '%s' ... ", cfg->skyGridFile );
		if ( (ret = boinc_delete_file ( cfg->skyGridFile ) ) != 0 )
		  {
		    LogPrintf (LOG_CRITICAL,  "Error removing '%s'. Return value: %d\n", cfg->skyGridFile, ret);
		    LogPrintf (LOG_CRITICAL, "System says: %s\n", strerror (errno) );
		    ABORT (status, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
		  }
		zipfile = resolved_name;
		LogPrintfVerbatim ( LOG_DEBUG, "done.\n");
	      } /* if boinc-link is different from physical file */
	    else
	      {	/* same file: move to temporary file to allow replacing original */
		LogPrintf (LOG_DEBUG, "Now rename'ing '%s' to '%s' ... ", resolved_name, tmpfile );
		if ( (ret = boinc_rename ( resolved_name, tmpfile)) != 0 )
		  {
		    LogPrintf (LOG_CRITICAL, "Error in renaming '%s' to '%s'. boinc_rename() returned %d.\n",
			       resolved_name, tmpfile, ret);
		    LogPrintf (LOG_CRITICAL, "System says: %s\n", strerror (errno) );
		  }
		else
		  LogPrintfVerbatim (LOG_DEBUG, " done.\n");

		zipfile = tmpfile;

	      } /* if boinc-link == physical file */

	    /* find directory of boinc-link file */
	    if ( (outdir = LALCalloc (1, strlen ( cfg->skyGridFile ) + 1 )) == NULL ) {
 	      ABORT (status, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
 	    }
 	    strcpy ( outdir, cfg->skyGridFile ); 
 	    /* simply truncate at last '/' */
 	    if ( (ptr = strrchr ( outdir, '/' )) == NULL ) /* no  path? -> local dir */
 	      strcpy ( outdir, "./" );
 	    else
 	      *(ptr + 1) = 0;	/* truncate outdir to basename-path */
 
	    LogPrintf (LOG_DEBUG, "Now boinc_unzipping '%s' to '%s' ... ", zipfile, outdir );
	    if ( (ret = boinc_zip ( UNZIP_IT, zipfile, outdir ) ) != 0 )
	      {
		LogPrintf (LOG_CRITICAL,  "Error in unzipping file '%s'. Return value: %d\n", 
			   zipfile, ret);
		LogPrintf (LOG_CRITICAL, "System says: %s\n", strerror (errno) );
		ABORT (status, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
	      }

	    LALFree(outdir);

	    /* the original 'cfg->skyGridfile' should now point to the correct 
	     * locally unzipped skygrid-file */
	    LogPrintfVerbatim ( LOG_DEBUG, " done.\n");
	    
	  } /* if is_zipped */

      } /* if uvar_skyGridFile */
    
#endif

  /* ----------------------------------------------------------------------*/
  /*
   * initialize Ephemeris-data 
   */
  {
    cfg->edat = (EphemerisData*)LALCalloc(1, sizeof(EphemerisData));
    cfg->edat->ephiles.earthEphemeris = cfg->EphemEarth;
    cfg->edat->ephiles.sunEphemeris = cfg->EphemSun;

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
  LogPrintf ( LOG_DETAIL, "# SFT time baseline:                  %f min\n",header.tbase/60.0);
  LogPrintf ( LOG_DETAIL, "# Starting search frequency:          %f Hz\n", cfg->spinRange.fkdot[0]);
  LogPrintf ( LOG_DETAIL, "# Demodulation frequency band:        %f Hz\n",cfg->spinRange.fkdotBand[0]);
  LogPrintf ( LOG_DETAIL, "# Actual # of SFTs:                   %d\n", cfg->SFTno);
  LogPrintf ( LOG_DETAIL, "# total observation time:             %f hours\n",1.0*(cfg->Tf-cfg->Ti)/3600.0);

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* InitFStat() */


/*----------------------------------------------------------------------*/
/** Some general consistency-checks on user-input.
 * Throws an error plus prints error-message if problems are found.
 */
void
checkUserInputConsistency (LALStatus *lstat)
{

  INITSTATUS(lstat);  

  if (!uvar_mergedSFTFile && (uvar_startTime>-1.0e308 || uvar_endTime<1.0e308)) {
    LogPrintf (LOG_CRITICAL,  "The --startTime and --endTIme options may ONLY be used\n"
                    "with merged SFT files, specified with the -B option.\n"
                    "Try ./ComputeFStatistic -h \n\n");
    ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);      
  }
  
  if(!uvar_DataDir && !uvar_mergedSFTFile && !uvar_DataFiles)
    {
      LogPrintf (LOG_CRITICAL,  "Must specify 'DataDir' OR 'mergedSFTFile' OR 'DataFiles'\n"
                      "No SFT directory specified; input directory with -D option.\n"
                      "No merged SFT file specified; input file with -B option.\n"
                      "Try ./ComputeFStatistic -h \n\n");
      ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);      
    }

  if(uvar_DataDir && uvar_mergedSFTFile)
    {
      LogPrintf (LOG_CRITICAL,  "Cannot specify both 'DataDir' and 'mergedSFTfile'.\n"
                      "Try ./ComputeFStatistic -h \n\n" );
      ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
    }      

  if (uvar_ephemYear == NULL)
    {
      LogPrintf (LOG_CRITICAL, "No ephemeris year specified (option 'ephemYear')\n\n");
      ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
    }      

  /* don't allow negative bands (for safty in griding-routines) */
  if ( (uvar_AlphaBand < 0) ||  (uvar_DeltaBand < 0) )
    {
      LogPrintf (LOG_CRITICAL, "Negative value of sky-bands not allowed (alpha or delta)!\n\n");
      ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
    }
  /* check for negative stepsizes in Freq, Alpha, Delta */
  if ( LALUserVarWasSet(&uvar_dAlpha) && (uvar_dAlpha < 0) )
    {
      LogPrintf (LOG_CRITICAL, "Negative value of stepsize dAlpha not allowed!\n\n");
      ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar_dDelta) && (uvar_dDelta < 0) )
    {
      LogPrintf (LOG_CRITICAL, "Negative value of stepsize dDelta not allowed!\n\n");
      ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar_dFreq) && (uvar_dFreq < 0) )
    {
      LogPrintf (LOG_CRITICAL, "Negative value of stepsize dFreq not allowed!\n\n");
      ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
    }

  /* grid-related checks */
  {
    BOOLEAN haveAlphaBand = LALUserVarWasSet( &uvar_AlphaBand );
    BOOLEAN haveDeltaBand = LALUserVarWasSet( &uvar_DeltaBand );
    BOOLEAN haveSkyRegion, haveAlphaDelta, haveSkyGridFile, useSkyGridFile, haveMetric, useMetric;

    haveSkyRegion  = (uvar_skyRegion != NULL);
    haveAlphaDelta = (LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta) );
    haveSkyGridFile   = (uvar_skyGridFile != NULL);
    useSkyGridFile   = (uvar_gridType == GRID_FILE_SKYGRID) || (uvar_gridType == GRID_METRIC_SKYFILE);
    haveMetric     = (uvar_metricType > LAL_PMETRIC_NONE);
    useMetric     = (uvar_gridType == GRID_METRIC) || (uvar_gridType == GRID_METRIC_SKYFILE);

    if ( (haveAlphaBand && !haveDeltaBand) || (haveDeltaBand && !haveAlphaBand) )
      {
	LogPrintf (LOG_CRITICAL, "ERROR: Need either BOTH (AlphaBand, DeltaBand) or NONE.\n\n"); 
        ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
      }

    if ( !useSkyGridFile && !(haveSkyRegion || haveAlphaDelta) )
      {
        LogPrintf (LOG_CRITICAL, "Need sky-region: either use (Alpha,Delta) or skyRegion!\n\n");
        ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
      }
    if ( haveSkyRegion && haveAlphaDelta )
      {
        LogPrintf (LOG_CRITICAL, "Overdetermined sky-region: only use EITHER (Alpha,Delta)"
		       " OR skyRegion!\n\n");
        ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);
      }
    if ( useSkyGridFile && !haveSkyGridFile )
      {
        LogPrintf (LOG_CRITICAL, "ERROR: gridType=(FILE|METRIC_SKYFILE), but no skyGridFile specified!\n\n");
        ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);  
      }
    if ( !useSkyGridFile && haveSkyGridFile )
      {
        LALWarning (lstat, "WARNING: skyGridFile was specified but not needed ..."
		    " will be ignored\n");
      }
    if ( useSkyGridFile && (haveSkyRegion || haveAlphaDelta) )
      {
        LALWarning (lstat, "WARNING: We are using skyGridFile, but sky-region was"
		    " also specified ... will be ignored!\n");
      }
    if ( !useMetric && haveMetric) 
      {
        LALWarning (lstat, "WARNING: Metric was specified for non-metric grid..."
		    " will be ignored!\n");
      }
    if ( useMetric && !haveMetric) 
      {
        LogPrintf (LOG_CRITICAL, "Metric grid-type selected, but no metricType selected\n\n");
        ABORT (lstat, COMPUTEFSTAT_EINPUT, COMPUTEFSTAT_MSGEINPUT);      
      }

  } /* grid-related checks */

  RETURN (lstat);
} /* checkUserInputConsistency() */


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

    INITSTATUS(stat);
    ATTATCHSTATUSPTR (stat);

    /* prepare log-file for writing */
    len = strlen(head) + strlen(".log") +10;
    if (uvar_outputLabel)
      len += strlen(uvar_outputLabel);

    if ( (fname = (CHAR*)LALCalloc(len,1)) == NULL) {
      ABORT (stat, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
    }
    strcpy (fname, head);
    if (uvar_outputLabel)
      strcat (fname, uvar_outputLabel);
    strcat (fname, ".log");

    if ( (fplog = fopen(fname, "wb" )) == NULL) {
      LogPrintf (LOG_CRITICAL, "Failed to open log-file '%f' for writing.\n\n", fname);
      LALFree (fname);
      ABORT (stat, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
    }

    /* write out a log describing the complete user-input (in cfg-file format) */
    TRY (LALUserVarGetLog(stat->statusPtr, &logstr,  UVAR_LOGFMT_CFGFILE), stat);

    fprintf (fplog, "## LOG-FILE of ComputeFStatistic run\n\n");
    fprintf (fplog, "# User-input:\n");
    fprintf (fplog, "# ----------------------------------------------------------------------\n\n");

    fprintf (fplog, "%s", logstr);
    LALFree (logstr);

    /* append an ident-string defining the exact CVS-version of the code used */
    fprintf (fplog, "\n\n# CVS-versions of executable:\n");
    fprintf (fplog, "# ----------------------------------------------------------------------\n");
    fclose (fplog);
    
    sprintf (command, "ident %s 2> /dev/null | sort -u >> %s", argv[0], fname);
    /* we don't fail here. If system() fails, we assume that */
    /* one of the system-commands was not available, and */
    /* therefore the CVS-versions will not be logged */
    if ( system(command) )
      LogPrintf ( LOG_DEBUG, "\nsystem('%s') returned non-zero status!\n", command );

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

  INITSTATUS(status);
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

  INITSTATUS(status);
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

  /* Free upsampled SFt data */
  if ( uvar_expLALDemod == 3)
    {
      /* Free UpSFTData */
      for (k=0;k<GV.SFTno;k++)
	{
	  /* data */
	  LALFree(UpSFTData[k]->fft->data->data);
	  LALFree(UpSFTData[k]->fft->data);
	  LALFree(UpSFTData[k]->fft);
	  LALFree(UpSFTData[k]);
	} /* for k < SFTno */
      LALFree(UpSFTData);
      UpSFTData = NULL;
    }


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

  if (GV.skyRegionString)
    {
      LALFree ( GV.skyRegionString );
      GV.skyRegionString = NULL;
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

  if ( GV.skyGridFile )
    LALFree ( GV.skyGridFile );

  /* Free ephemeris data */
  LALFree(GV.edat->ephemE);
  LALFree(GV.edat->ephemS);
  LALFree(GV.edat);
  GV.edat = NULL;

  /* free buffer used for fstat.  Its safe to do this because we already did fclose(fpstat) earlier */
  if (fstatbuff) {
    LALFree(fstatbuff);
    fstatbuff = NULL;
  }
  
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
INT4 PrintTopValues(REAL8 TwoFthr, INT4 ReturnMaxN, PulsarDopplerParams searchpos)
{

  INT4 *indexes,i,j,iF,N;
  REAL8 mean=0.0, std=0.0 ,log2val /*=log(2.0)*/;
  
  log2val = medianbias;

  if ( lalDebugLevel > 10)	/* dummy: avoid warnings */
    XLALPrintError (  "%f, %f", searchpos.Alpha, searchpos.Delta);

  for (i=0;i<highFLines->Nclusters;i++){
    N=highFLines->NclustPoints[i];
    for (j=0;j<N;j++){
      iF=highFLines->Iclust[j];
      Fstat.F[iF]=0.0;
    }/*  end j loop over points of i-th cluster  */
  }/*  end i loop over different clusters */

/*    check that ReturnMaxN is sensible */
  if (ReturnMaxN>GV.FreqImax){
    LogPrintf (LOG_CRITICAL,  "PrintTopValues() WARNING: resetting ReturnMaxN=%d to %d\n",
            ReturnMaxN, GV.FreqImax);
    ReturnMaxN=GV.FreqImax;
  }

/*    create an array of indexes */
  if (!(indexes=(INT4 *)LALMalloc(sizeof(INT4)*GV.FreqImax))){
    LogPrintf (LOG_CRITICAL,  "Unable to allocate index array in PrintTopValues()\n");
    return 1;
  }

/*    populate it */
  for (i=0;i<GV.FreqImax;i++)
    indexes[i]=i;

/*   sort array of indexes */
  qsort((void *)indexes, (size_t)GV.FreqImax, sizeof(int), compare);


/*    Normalize */
  TwoFthr*=0.5/log2val;

#ifdef FINDFWHM /* Find the full-width-at-half-maximum */
  {INT4 ic;
  REAL8 leftEdgeFreq=0.0, rightEdgeFreq=10000.0, outlierFWHM=0.0;
  REAL8 leftEdgeTwoF=0.0, rightEdgeTwoF=0.0, maximumTwoF=0.0;

  ic=0;
  while( Fstat.F[ic] <= 0.5 * Fstat.F[indexes[0]] ) {
    ic++;
  } 
  if( ic==0 ) {
    LogPrintf (LOG_CRITICAL,  "Warning: Search frequency band may be too small to cover the outlier.\n");
    leftEdgeTwoF = Fstat.F[ic];
    leftEdgeFreq = GV.spinRange.fkdot[0] + ic*GV.dFreq;
  } else {/* shift slightly downwards to cover the half-maximum point, if possible */
    leftEdgeTwoF = Fstat.F[ic-1];
    leftEdgeFreq = GV.spinRange.fkdot[0] + (ic-1)*GV.dFreq;
  }


  ic=GV.FreqImax-1;
  while( Fstat.F[ic] <= 0.5 * Fstat.F[indexes[0]] ) {
    ic--;
  } 
  if( ic==GV.FreqImax-1 ) {
    LogPrintf (LOG_CRITICAL,  "Warning: Search frequency band may be too small to cover the outlier.\n");
    rightEdgeTwoF = Fstat.F[ic];
    rightEdgeFreq = GV.spinRange.fkdot[0] + ic*GV.dFreq;
  } else { /* shift slightly upwards to cover the half-maximum point, if possible */
    rightEdgeTwoF = Fstat.F[ic+1];
    rightEdgeFreq = GV.spinRange.fkdot[0] + (ic+1)*GV.dFreq;
  }

  maximumTwoF= Fstat.F[indexes[0]];
  outlierFWHM = rightEdgeFreq - leftEdgeFreq;
  fprintf(stdout,"%15.9f %22.6f %15.9f %22.6f %22.6f %15.9f\n",leftEdgeFreq,leftEdgeTwoF,rightEdgeFreq,rightEdgeTwoF,maximumTwoF,outlierFWHM);
  }
#endif

#ifdef FILE_FMAX
  {
    INT4 ntop,err;
    /*    print out the top ones */
    for (ntop=0; ntop<ReturnMaxN; ntop++)
      if (Fstat.F[indexes[ntop]]>TwoFthr){
        err=fprintf(fpmax, "%20.10f %10.8f %10.8f %20.15f\n",
                    GV.spinRange.fkdot[0] + indexes[ntop]*GV.dFreq,
                    searchpos.Alpha, searchpos.Delta, 2.0*log2val*Fstat.F[indexes[ntop]]);
        if (err<=0) {
          LogPrintf (LOG_CRITICAL,  "PrintTopValues couldn't print to Fmax!\n");
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
  mean*=(2.0*log2val);
  std=2.0*log2val*sqrt(std);


/*    Find number of values above threshold (could be made O(log N) */
/*    with binary search! */
  for (i=0; i<GV.FreqImax; i++)
    if (Fstat.F[indexes[i]]<=TwoFthr)
      break;

#ifdef FILE_FMAX_DEBG    
/*    print the output */
  err=fprintf(fpmax,"%10.5f %10.8f %10.8f    %d %10.5f %10.5f %10.5f\n",
	      GV.spinRange.fkdot[0], searchpos.Alpha, searchpos.Delta, 
	      GV.FreqImax-N, mean, std, 2.0*log2val*Fstat.F[indexes[0]]);
#endif
  LALFree(indexes);
#ifdef FILE_FMAX_DEBG    
  if (err<=0) {
    LogPrintf (LOG_CRITICAL,  "PrintTopValues couldn't print to Fmax!\n");
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
  REAL8 THR=10.0;
  
  REAL8 dmp;

  OutliersInput  *outliersInput;
  OutliersParams *outliersParams;
  Outliers       *outliers;
  ClustersInput  *clustersInput;
  ClustersParams *SpClParams;
  Clusters       *SpLines=highFLines;
    
  INT2 smallBlock=1;
  INT4 wings;

  INITSTATUS(stat);
  ATTATCHSTATUSPTR (stat);

  nbins=(UINT4)nbins;

  THR=uvar_Fthreshold;

/* 0.0002 is the max expected width of the F stat curve for signal */
/* with ~ 10 h observation time */

  dmp=0.5+0.0002/ DemodParams->df;
  wings = (INT4) dmp;


#ifdef FILE_FTXT
  /*  file contains freq, PSD, noise floor */
  if( (outfile = fopen("F.txt","wb")) == NULL) {
    ABORT (stat, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
  }
#endif
  /*  file contanis Freq, PSD, noise floor,lines */
#ifdef FILE_FLINES  
  if( (outfile1 = fopen("FLines.txt","wb")) == NULL) {
    ABORT (stat, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
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
    ABORT (stat, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
  }
  outliers->Noutliers=0;

  if ( (outliersParams = (OutliersParams *)LALCalloc(1,sizeof(OutliersParams))) == NULL) {
    ABORT (stat, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
  }
  if ( (outliersInput = (OutliersInput *)LALCalloc(1,sizeof(OutliersInput))) == NULL) {
    ABORT (stat, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
  }
  
  outliersParams->Thr=THR/(2.0*medianbias);
  outliersParams->Floor = FloorF1;
  outliersParams->wings=wings; /*these must be the same as ClustersParams->wings */
  outliersInput->ifmin = (INT4) ((GV.spinRange.fkdot[0] /DemodParams->df)+0.5);
  outliersInput->data = F1;

  /*find values of F above THR and populate outliers with them */
  if (ComputeOutliers(outliersInput, outliersParams, outliers)) {
    ABORT (stat, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
  }

  /*if no outliers were found clean and exit */
   if (outliers->Noutliers == 0){

#ifdef FILE_FTXT
     /*  F.txt file contains Freq, F, noise floor of F   */
     for (i=0;i<nbins;i++){ 
       REAL4 Freq;
       REAL8 r0,r1;
       Freq=GV.fkdot0->data[0] + i*GV.dFreq;
       r0=F1->data[i]*2.0*medianbias;
       r1=FloorF1->data[i]*2.0*medianbias;
       fprintf(outfile,"%f %E %E\n",Freq,r0,r1);
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

     /*      XLALPrintError ( "Nclusters zero \n"); */
     /*      fflush(stderr); */

     goto finished;

   } /* if Noutliers == 0 */
  

   /* if outliers are found get ready to identify clusters of outliers*/
   if ( (SpClParams = (ClustersParams*)LALCalloc(1,sizeof(ClustersParams))) == NULL) {
     ABORT (stat, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
   }
   
   if ( (clustersInput = (ClustersInput *)LALCalloc(1,sizeof(ClustersInput))) == NULL) {
     ABORT (stat, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
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
     REAL8 Freq;
     REAL8 r0,r1,r2;
     j=SpLines->Iclust[i];
     Freq=(GV.fkdot0->data[0]+SpLines->Iclust[i]*GV.dFreq);
     r0=F1->data[j];
     r1=FloorF1->data[j]*2.0*medianbias;
     r2=SpLines->clusters[i]*FloorF1->data[j]*2.0*medianbias;
     fprintf(outfile1,"%20.17f %E %E %E\n",Freq,r0,r1,r2);
   }
#endif
#ifdef FILE_FTXT   
   /*  PSD.txt file contains Freq, PSD, noise floor   */
   for (i=0;i<nbins;i++){ 
     REAL8 Freq;
     REAL8 r0,r1;
     Freq=GV.fkdot0->data[0] + i*GV.dFreq;
     r0=F1->data[i]*2.0*medianbias;
     r1=FloorF1->data[i]*2.0*medianbias;
     fprintf(outfile,"%20.17f %E %E\n",Freq,r0,r1);
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
NormaliseSFTDataRngMdn(LALStatus *stat, INT4 windowSize)
{
#ifdef FILE_SPRNG  
  FILE *outfile;
#endif
  INT4 i,j,m,lpc,il;                         /* loop indices */
  INT4 Ntot,nbins=GV.ifmax-GV.ifmin+1;   /* Number of points in SFT's */
  REAL8Vector *Sp=NULL, *RngMdnSp=NULL;   /* |SFT|^2 and its rngmdn  */
  REAL8 B;                          /* SFT Bandwidth */
  REAL8 deltaT,norm,*N, *Sp1;
  REAL4 xre,xim,xreNorm,ximNorm;

  INITSTATUS(stat);
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

  if( (N = (REAL8 *) LALCalloc(nbins,sizeof(REAL8))) == NULL) {
    ABORT (stat, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
  }
  if( (Sp1 = (REAL8 *) LALCalloc(nbins,sizeof(REAL8))) == NULL) {
    ABORT (stat, COMPUTEFSTAT_EMEM, COMPUTEFSTAT_MSGEMEM);
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
        xre=crealf(SFTData[i]->fft->data->data[j]);
        xim=cimagf(SFTData[i]->fft->data->data[j]);
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
        xre=crealf(SFTData[i]->fft->data->data[j]);
        xim=cimagf(SFTData[i]->fft->data->data[j]);
        xreNorm=N[j]*xre; 
        ximNorm=N[j]*xim; 
        SFTData[i]->fft->data->data[j].realf_FIXME = xreNorm;    
        SFTData[i]->fft->data->data[j].imagf_FIXME = ximNorm;
        Sp1[j]=Sp1[j]+xreNorm*xreNorm+ximNorm*ximNorm;
      }
      
    } /* end loop over SFTs*/

#ifdef FILE_SPRNG  
  if( (outfile = fopen("SpRng.txt","wb")) == NULL) {
    ABORT (stat, COMPUTEFSTAT_ESYS, COMPUTEFSTAT_MSGESYS);
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
  char resolved_name[MAXFILENAMELENGTH];
  if (boinc_resolve_filename(orig_name, resolved_name, sizeof(resolved_name))) {
    LogPrintf (LOG_CRITICAL,   
            "Can't resolve file '%s'\n"
            "If running a non-BOINC test, create [INPUT] or touch [OUTPUT] file\n",
            orig_name);

    boinc_finish(COMPUTEFSTAT_EXIT_BOINCRESOLVE);
  }
  strcpy(orig_name, resolved_name);
  return;
}

void use_boinc_filename1(char **orig_name ) {
  char resolved_name[MAXFILENAMELENGTH];
  if (boinc_resolve_filename(*orig_name, resolved_name, sizeof(resolved_name))) {
    LogPrintf (LOG_CRITICAL,   
            "Can't resolve file '%s'\n"
            "If running a non-BOINC test, create [INPUT] or touch [OUTPUT] file\n",
            *orig_name);
    boinc_finish(COMPUTEFSTAT_EXIT_BOINCRESOLVE);
  }
  LALFree(*orig_name);
  *orig_name = (CHAR*) LALCalloc(strlen(resolved_name)+1,1);
  strcpy(*orig_name, resolved_name);
  return;
}

int globargc=0;
char **globargv=NULL;

void worker(void) {

#if (BOINC_GRAPHICS == 2) 
  if (graphics_lib_handle) {
    if (!(set_search_pos_hook = dlsym(graphics_lib_handle,"set_search_pos"))) {
      LogPrintf (LOG_CRITICAL,   "unable to resolve set_search_pos(): %s\n", dlerror());
      boinc_finish(COMPUTEFSTAT_EXIT_DLOPEN);
    }
    if (!(fraction_done_hook = dlsym(graphics_lib_handle,"fraction_done"))) {
      LogPrintf (LOG_CRITICAL,   "unable to resolve fraction_done(): %s\n", dlerror());
      boinc_finish(COMPUTEFSTAT_EXIT_DLOPEN);
    }
  }
  else
    LogPrintf (LOG_CRITICAL,  "graphics_lib_handle NULL: running without graphics\n");
#endif

#ifndef RUN_POLKA
  int retval=boincmain(globargc,globargv);
#ifdef CLUSTERED_OUTPUT
  Outputfilename=CFstatFilename;
#else
  Outputfilename=FstatFilename;
#endif
#else
  int a1,a2,retval;
  CHAR ckptfname1[MAXFILENAMELENGTH+4];
  CHAR Fstatsfilename1[MAXFILENAMELENGTH+4];

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
#ifdef CLUSTERED_OUTPUT
    strncpy(Fstatsfilename1,CFstatFilename,sizeof(Fstatsfilename1));
#else
    strncpy(Fstatsfilename1,FstatFilename,sizeof(Fstatsfilename1));
#endif
    strncpy(ckptfname1,ckp_fname,sizeof(ckptfname1));
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
      remove (FstatFilename);
      remove (Fstatsfilename1);
    */
  } else {
#ifdef CLUSTERED_OUTPUT
    Outputfilename=CFstatFilename;
#else
    Outputfilename=FstatFilename;
#endif
    remove (ckp_fname);
  }
#endif
#if BOINC_COMPRESS
  /* compress the file if it exists */
  if (uvar_useCompression && !retval) 
    {
      int boinczipret;
      char* zipname = "temp.zip";

      LogPrintf (LOG_DEBUG, "Now boinc_deleting '%s' ... ", zipname);
      boinczipret = boinc_delete_file(zipname);
      if (boinczipret) 
	{
	  LogPrintf (LOG_CRITICAL,  "can't remove old zip file %s. not zipping output.\n",zipname);
	}
      else 
	{
	  LogPrintfVerbatim (LOG_DEBUG, "done.\n");
	  LogPrintf (LOG_DEBUG, "Now boinc_zip'ing '%s' to '%s' ... ", Outputfilename, zipname );
	  boinczipret=boinc_zip(ZIP_IT, zipname , Outputfilename);
	  if (boinczipret) 
	    {
	      LogPrintf (LOG_CRITICAL,   "Error in zipping file %s to temp.zip. Return value %d. not zipping output.\n",
			 Outputfilename, boinczipret);
	    } 
	  else 
	    {
	      LogPrintfVerbatim (LOG_DEBUG, "done.\n");
	      LogPrintf (LOG_DEBUG, "Now boinc_rename'ing '%s' to '%s' ... ", zipname, Outputfilename );
	      boinczipret=boinc_rename(zipname, Outputfilename);
	      if (boinczipret) 
		{
		  LogPrintf (LOG_CRITICAL,   "Error in renaming file %s to %s. boinc_rename() returned %d. not zipping output.\n",
			     zipname, Outputfilename, boinczipret);
		  LogPrintf (LOG_CRITICAL, "System says: %s\n", strerror (errno) );
		}
	      else
		{
		  LogPrintfVerbatim (LOG_DEBUG, " done.\n");
		}
	    } /* if !boinczipret */

	} /* if !boinczipret (deleting) */

    } /* if useCompression && ok */
#endif

  /* report final WU-fpops to client if given */
  if ( LALUserVarWasSet ( &uvar_WUfpops ) ) 
    boinc_ops_cumulative( uvar_WUfpops, 0); /* ignore IOPS */

  boinc_finish(retval);
  return;

} /* worker() */


/*********************************************************************************
 * main() in case of USE_BOINC
 */
int main(int argc, char *argv[])
{
  FILE *fp_debug=NULL;
  int skipsighandler=0;

/* setup windows diagnostics (e.g. redirect stderr into a file!) */
#ifdef _WIN32
  boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED |
                         BOINC_DIAG_HEAPCHECKENABLED |
                         BOINC_DIAG_ARCHIVESTDERR |
                         BOINC_DIAG_REDIRECTSTDERR |
                         BOINC_DIAG_TRACETOSTDERR);
#endif


#define DEBUG_LEVEL_FNAME "CFS_DEBUG_LEVEL"
#define DEBUG_DDD_FNAME   "CFS_DEBUG_DDD"

  LogPrintfVerbatim (LOG_NORMAL, "\n");
  LogPrintf (LOG_NORMAL, "Start of BOINC application '%s'.\n", argv[0]);

  /* see if user has a DEBUG_LEVEL_FNAME file: read integer and set lalDebugLevel */
  if ((fp_debug=fopen("../../" DEBUG_LEVEL_FNAME, "r")) || (fp_debug=fopen("./" DEBUG_LEVEL_FNAME, "r")))
    {
      int read_int;

      LogPrintf (LOG_NORMAL, "Found '%s' file\n", DEBUG_LEVEL_FNAME);
      if ( 1 == fscanf(fp_debug, "%d", &read_int ) ) 
	{
	  LogPrintf (LOG_NORMAL, "...containing int: Setting lalDebugLevel -> %d\n", read_int );
	  lalDebugLevel = read_int;
	}
      else
	{
	  LogPrintf (LOG_NORMAL, "...with no parsable int: Setting lalDebugLevel -> 1\n");
	  lalDebugLevel = 1;
	}
      fclose (fp_debug);

    } /* if DEBUG_LEVEL_FNAME file found */

  
#if defined(__GNUC__)
  /* see if user has created a DEBUG_DDD_FNAME file: turn on debuggin using 'ddd' */
  if ((fp_debug=fopen("../../" DEBUG_DDD_FNAME, "r")) || (fp_debug=fopen("./" DEBUG_DDD_FNAME, "r")) ) 
    {
      char commandstring[256];
      char resolved_name[MAXFILENAMELENGTH];
      char *ptr;
      pid_t process_id=getpid();
      
      fclose(fp_debug);
      LogPrintf ( LOG_NORMAL, "Found '%s' file, trying debugging with 'ddd'\n", DEBUG_DDD_FNAME);
      
      /* see if the path is absolute or has slashes.  If it has
	 slashes, take tail name */
      if ((ptr = strrchr(argv[0], '/'))) {
	ptr++;
      } else {
	ptr = argv[0];
      }
      
      /* if file name is an XML soft link, resolve it */
      if (boinc_resolve_filename(ptr, resolved_name, sizeof(resolved_name)))
	LogPrintf (LOG_NORMAL,  "Unable to boinc_resolve_filename(%s), so no debugging\n", ptr);
      else {
	skipsighandler=1;
	snprintf(commandstring,sizeof(commandstring),"ddd %s %d &", resolved_name ,process_id);
	system(commandstring);
	sleep(20);
      }
    } /* DEBUGGING */
#endif // GNUC
  
  
  globargc=argc;
  globargv=argv;


#ifndef _WIN32
  /* install signal-handler for SIGTERM, SIGINT and SIGABRT(?) 
   * NOTE: it is critical to catch SIGINT, because a user
   * pressing Ctrl-C under boinc should not directly kill the
   * app (which is attached to the same terminal), but the app 
   * should wait for the client to send <quit/> and cleanly exit. 
   */
  boinc_set_signal_handler(SIGTERM, sighandler);
  boinc_set_signal_handler(SIGINT,  sighandler);
  boinc_set_signal_handler(SIGABRT, sighandler);

  /* install signal handler (for ALL threads) for catching
   * Segmentation violations, floating point exceptions, Bus
   * violations and Illegal instructions */
  if ( !skipsighandler )
    {
      boinc_set_signal_handler(SIGSEGV, sighandler);
      boinc_set_signal_handler(SIGFPE,  sighandler);
      boinc_set_signal_handler(SIGILL,  sighandler);
      boinc_set_signal_handler(SIGBUS,  sighandler);
    } /* if !skipsighandler */
#else /* WIN32 */
  signal(SIGTERM, sighandler);
  signal(SIGINT,  sighandler);
  signal(SIGABRT, sighandler);
  if ( !skipsighandler ) {
      signal(SIGSEGV, sighandler);
      signal(SIGFPE,  sighandler);
      signal(SIGILL,  sighandler);
  }
#endif /* WIN32 */
  
#if (BOINC_GRAPHICS == 1)
  set_search_pos_hook = set_search_pos;
  fraction_done_hook = &fraction_done;
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
    LogPrintf (LOG_CRITICAL,  "boinc_init_graphics[_lib]() returned %d. This indicates an error...\n", retval);
    boinc_finish(COMPUTEFSTAT_EXIT_WORKER );
  }
#endif /*  BOINC_GRAPHICS>0 */
    
    boinc_init();
    worker();
    boinc_finish(COMPUTEFSTAT_EXIT_WORKER );
    /* we never get here!! */
    return 0;
}
  
  
#ifdef __GLIBC__
  /* needed to define backtrace() which is glibc specific*/
#include <execinfo.h>
#endif /* __GLIBC__ */

/* signal handlers */
void sighandler(int sig){
#ifdef __GLIBC__
  void *array[64];
  size_t size;
#endif
  static int killcounter = 0;

  /* RP: not sure what this is for. FIXME: better remove?
#ifndef _WIN32
  sigset_t signalset;
  sigemptyset(&signalset);
  sigaddset(&signalset, sig);
  pthread_sigmask(SIG_BLOCK, &signalset, NULL);
#endif
  */

  LALStatus *mystat = &global_status;     /* the only place in the universe where we're
                                           * allowed to use the global_status struct !! */
  /* lets start by ignoring ANY further occurences of this signal
     (hopefully just in THIS thread, if truly implementing POSIX threads */
  LogPrintfVerbatim(LOG_CRITICAL, "\n");
  LogPrintf (LOG_CRITICAL, "APP DEBUG: Application caught signal %d.\n\n", sig );

  /* ignore TERM interrupts once  */
  if ( sig == SIGTERM || sig == SIGINT )
    {
      killcounter ++;

      if ( killcounter >= 4 )
        {
          LogPrintf (LOG_CRITICAL, "APP DEBUG: got 4th kill-signal, guess you mean it. Exiting now\n\n");
          boinc_finish(COMPUTEFSTAT_EXIT_USER);
        }
      else
        return;

    } /* termination signals */

  if (mystat)
    LogPrintf (LOG_CRITICAL,   "Stack trace of LAL functions in worker thread:\n");
  while (mystat) {
    LogPrintf (LOG_CRITICAL,   "%s at line %d of file %s\n", mystat->function, mystat->line, mystat->file);
    if (!(mystat->statusPtr)) {
      const char *p=mystat->statusDescription;
      LogPrintf (LOG_CRITICAL,   "At lowest level status code = %d, description: %s\n", mystat->statusCode, p?p:"NO LAL ERROR REGISTERED");
    }
    mystat=mystat->statusPtr;
  }
  
#ifdef __GLIBC__
  /* now get TRUE stacktrace */
  size = backtrace (array, 64);
  LogPrintf (LOG_CRITICAL,   "Obtained %zd stack frames for this thread.\n", size);
  LogPrintf (LOG_CRITICAL,   "Use gdb command: 'info line *0xADDRESS' to print corresponding line numbers.\n");
  backtrace_symbols_fd(array, size, filenum(stderr));
#endif /* __GLIBC__ */
  /* sleep a few seconds to let the OTHER thread(s) catch the signal too... */
  sleep(5);
  boinc_finish(COMPUTEFSTAT_EXIT_SIGNAL);
  return;
} /* sighandler */
#endif /*USE_BOINC*/


/** Check presence and consistency of checkpoint-file and use to set loopcounter if valid.
 *
 *  The name of the checkpoint-file is FNAME.ckp
 */
void
getCheckpointCounters(LALStatus *stat,		/**< pointer to LALStatus structure */
                      UINT4 *loopcounter,	/**< [out] number of completed loops (refers to main-loop in main()) */
                      UINT4 *checksum,		/**< [out] checksum of file (up the bytecounter bytes) */
                      long *bytecounter,	/**< [out] bytes nominally written to fstats file (for consistency-check) */
                      const CHAR *fstat_fname,	/**< [in] Name of Fstats-file */
                      const CHAR *ckpfn)	/**< [in] check-point file name */
{
  FILE *fp;
  UINT4 lcount;         /* loopcounter */
  long bcount;          /* and bytecounter read from ckp-file */
  long flen;
  char lastnewline='\0';
  UINT4 cksum;           /* as read for .ckp file */
  UINT4 computecksum;    /* as computed from contents of Fstats file */
  int i;
#ifdef DEBUG_CHECKPOINTING
  int savelaldebuglevel=lalDebugLevel;
  lalDebugLevel=1;
#endif
 
  INITSTATUS(stat);
  ASSERT ( fstat_fname, stat, COMPUTEFSTAT_ENULL, COMPUTEFSTAT_MSGENULL);
  ASSERT ( ckpfn, stat, COMPUTEFSTAT_ENULL, COMPUTEFSTAT_MSGENULL);
  ASSERT ( loopcounter, stat, COMPUTEFSTAT_ENULL, COMPUTEFSTAT_MSGENULL);
  ASSERT ( bytecounter, stat, COMPUTEFSTAT_ENULL, COMPUTEFSTAT_MSGENULL);

  /* if anything goes wrong in here: start main-loop from beginning  */
  *loopcounter = 0;     
  *bytecounter = 0;

  
  /* try opening checkpoint-file read-only */
  if (!(fp = fopen(ckpfn, "rb"))) 
    {
      LogPrintf (LOG_NORMAL, "Checkpoint-file '%s' not found.\n", ckpfn);
      RETURN(stat);
    }
  LogPrintf (LOG_NORMAL, "Found checkpoint-file '%s' \n", ckpfn);
  
  /* try reading checkpoint-counters: three INT's loopcounter, checksum, and fstat_bytecounter */
  if ( 4 != fscanf (fp, "%" LAL_UINT4_FORMAT " %" LAL_UINT4_FORMAT " %ld\nDONE%c", &lcount, &cksum, &bcount, &lastnewline) || lastnewline!='\n') 
    {
      LogPrintfVerbatim (LOG_CRITICAL, "Failed to read checkpoint-counters from '%s'!\n", ckpfn);
      fclose(fp);
      RETURN(stat);
    }
  fclose( fp );
  LogPrintf (LOG_DEBUG, "Read checkpoint-counters (%d/%d/%ld)\n", 
	     lcount, cksum, bcount );

  /* checkpoint-file read successfully: check consistency with fstats-file */
  if ( (fp = fopen(fstat_fname, "rb")) == NULL )
    {
      LogPrintf (LOG_CRITICAL, "Opening Fstats file '%s' failed.\n", fstat_fname);
      RETURN(stat);
    }

  /* seek to end of fstats file */
  if ( fseek( fp, 0, SEEK_END) ) 
    {
      LogPrintf (LOG_CRITICAL, "fseek(SEEK_END) on Fstat-file '%s' failed!\n", fstat_fname);
      fclose(fp);
      RETURN(stat);
    }
  flen = ftell(fp);

  /* is bytecounter consistent with length of this file? */
  if ( bcount > flen) 
    {
      LogPrintf (LOG_CRITICAL, "Corrupted Fstat-file '%s': has %ld bytes instead of %ld.\n",
		 fstat_fname, flen, bcount);
      fclose(fp);
      RETURN(stat);
    }
  
  /* if we're using a toplist, read the Fstat-file into the toplist, giving also 
   * the checksum 
   */
  if ( toplist ) 
    {
      LogPrintf (LOG_NORMAL, "Trying to read Fstat-file into toplist ...\n");
      if ( fseek(fp, 0, SEEK_SET) )
	{ 
	  LogPrintf (LOG_CRITICAL, "fseek(SEEK_SET) failed on Fstat-file.\n");
	  fclose(fp);
	  RETURN(stat);
	}

      if( read_fstat_toplist_from_fp( toplist, fp, &computecksum, bcount) < 0 ) 
	{
	  LogPrintf (LOG_CRITICAL, "Couldn't read '%s' into toplist.\n", fstat_fname);
	  fclose(fp);
	  RETURN(stat);
	}
    } /* if toplist */
  else
    {
      /* compute checksum */
      computecksum = 0;
      if ( fseek( fp, 0, SEEK_SET ) )
	{
	  LogPrintf ( LOG_CRITICAL, "fseek(SEEK_SET) failed on Fstats-file '%s'\n", fstat_fname);
	  fclose(fp);
	  RETURN(stat);
	}
      for (i=0; i < bcount; i++) 
	{
	  int onechar = getc(fp);
	  if (onechar == EOF) 
	    {
	      fclose(fp);
	      RETURN(stat);
	    }
	  else
	    computecksum += onechar;
	} /* for i < bcount */
    
    } /* if !toplist */


  if ( computecksum != cksum) 
    {
      LogPrintf ( LOG_CRITICAL, "Fstats file '%s' has incorrect checksum.\n", fstat_fname );
      LogPrintf ( LOG_CRITICAL, "Checkpointed checksum: %d, computed: %d\n", cksum, computecksum );
      fclose(fp);
      RETURN(stat);
    }
  else
    LogPrintf (LOG_NORMAL, "Checksum Ok. Successfully read_fstat_toplist_from_fp()\n");

  
  *loopcounter = lcount;
  *bytecounter = bcount;
  *checksum=cksum;

#ifdef DEBUG_CHECKPOINTING
  lalDebugLevel=savelaldebuglevel;
#endif
  
  fclose(fp);
  RETURN(stat);

} /* getCheckpointCounters() */ 


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


/** Set up the search-grid and prepare DopplerSkyScan for stepping through parameter-space.
 * \note this is a bit ugly as it's using global uvar_ User-input variables.
 */
void
InitSearchGrid ( LALStatus *status, 
		 DopplerSkyScanState *scan, 
		 ConfigVariables *cfg)
{
  DopplerSkyScanInit scanInit = empty_DopplerSkyScanInit; /* init-structure for DopperScanner */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* Prepare input-structure for initialization of DopplerSkyScan
   */
  LogPrintf (LOG_DEBUG, "Setting up template grid ...\n");
  scanInit.metricType = (LALPulsarMetricType) uvar_metricType;
  scanInit.dAlpha = uvar_dAlpha;
  scanInit.dDelta = uvar_dDelta;
  scanInit.gridType = (DopplerGridType) uvar_gridType;
  scanInit.metricMismatch = uvar_metricMismatch;
  scanInit.obsBegin.gpsSeconds = cfg->Ti;
  scanInit.obsBegin.gpsNanoSeconds = 0;
  scanInit.obsDuration = (REAL8) (cfg->Tf - cfg->Ti);
  scanInit.projectMetric = uvar_projectMetric;
  scanInit.Detector = &cfg->Detector;
  scanInit.ephemeris = cfg->edat;         /* used by Ephemeris-based metric */
  scanInit.skyGridFile = cfg->skyGridFile;      /* if applicable */

  scanInit.skyRegionString = cfg->skyRegionString;
  scanInit.Freq = cfg->spinRange.fkdot[0] + cfg->spinRange.fkdotBand[0];

  /* the following call generates a skygrid plus determines dFreq, df1dot,..
   * Currently the complete template-grid is more like a 'foliation' of
   * the parameter-space by skygrids, stacked along f and f1dot, 
   * not a full 4D metric template-grid
   */
  TRY ( InitDopplerSkyScan( status->statusPtr, scan, &scanInit), status); 

  /* ---------- should we write the sky-grid to disk? ---------- */
  if ( uvar_outputSkyGrid ) 
    {
      LogPrintf (LOG_NORMAL,   "Now writing sky-grid into file '%s' ...", uvar_outputSkyGrid);
      TRY (writeSkyGridFile( status->statusPtr, scan->skyGrid, uvar_outputSkyGrid), status);
      LogPrintfVerbatim (LOG_NORMAL, " done.\n\n");
    }

  DETATCHSTATUSPTR (status);
  RETURN(status);
} /* InitSearchGrid() */

int
debug_dump_commandline (int argc,  char *argv[])
{
  int i;

  LogPrintf (LOG_DEBUG, "Commandline:\n");
  for (i=0; i < argc; i++)
    LogPrintfVerbatim ( LOG_DEBUG, "%s ", argv[i]);
  LogPrintfVerbatim (LOG_DEBUG, "\n");
  
  return 0;
} /* debug_dump_commandline() */


/* check if given file is a zip-file by checking the zip-magic header 'PK\003\044'
 * RETURN: 1 = true, 0 = false, -1 = error
 */
int
is_zip_file ( const char *fname )
{
  FILE *fp;
  CHAR zip_magic[] = {'P', 'K', 3, 4 };
  CHAR file_header[4];

  if ( (fp = fopen( fname, "rb")) == NULL )
    {
      LogPrintf (LOG_CRITICAL, "Failed to open '%s' for reading.\n", fname);
      return -1;
    }
	
  if ( 4 != fread ( file_header, sizeof(CHAR), 4, fp ) )
    {
      LogPrintf (LOG_CRITICAL, "Failed to read first 4 bytes from '%s'.\n", fname);
      return -1;
    }

  fclose(fp);

  if ( memcmp ( file_header, zip_magic, 4 ) )
    return 0;	/* false: no zip-file */
  else
    return 1;	/* yep, found magic zip-header */

} /* is_zip_file() */


	    
