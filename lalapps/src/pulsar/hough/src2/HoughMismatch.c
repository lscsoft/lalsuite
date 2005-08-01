/****************************************************
 *
 * Filename: HoughMismatch.c
 *
 * Revision: $Id$
 * 
 * Estimating mismatch of grid used in Hough search
 * 
 ****************************************************/


#include "../src/MCInjectHoughS2.h"


INT4 lalDebugLevel;

/* defaults */
#define EARTHEPHEMERIS "./earth00-04.dat"
#define SUNEPHEMERIS "./sun00-04.dat"
#define MAXFILES 3000 /* maximum number of files to read in a directory */
#define MAXFILENAMELENGTH 256 /* maximum # of characters  of a SFT filename */
#define IFO 2         /*  detector, 1:GEO, 2:LLO, 3:LHO */
#define THRESHOLD 1.6 /* thresold for peak selection, with respect to the */
#define NFSIZE  21 /* n-freq. span of the cylinder, to account for spin-down */
#define BLOCKSRNGMED 101 /* Running median window size */

/* default injected pulsar parameters */
#define F0 250.0          /*  frequency to build the LUT and start search */
#define FDOT 0.0 /* default spindown parameter */
#define ALPHA 0.0		/* center of the sky patch (in radians) */
#define DELTA  (-LAL_PI_2)
#define COSIOTA 0.5
#define PHI0 0.0
#define PSI 0.0
#define H0 (1.0e-23)

/* default file and directory names */
#define SFTDIRECTORY "/home/badkri/L1sfts"
#define FILEOUT "./MismatchOut"      /* prefix file output */
#define TRUE (1==1)
#define FALSE (1==0)

NRCSID (HOUGHMISMATCHC, "$Id$");

int main( int argc, char *argv[]){

  static LALStatus            status;  
  static LALDetector          detector;
  static LIGOTimeGPSVector    timeV;
  static REAL8Cart3CoorVector velV;
  static REAL8Vector          timeDiffV;
  static REAL8Vector          foft;

  /* data about injected signal */
  static PulsarData           pulsarInject;

  /* the template */
  static HoughTemplate        pulsarTemplate;

  FILE  *fpOUT = NULL; /* output file pointer */
  FILE  *fpLog = NULL; /* log file pointer */
  CHAR   *logstr=NULL; /* log string containing user input variables */
  CHAR *fnamelog=NULL; /* name of log file */
  INT4 nfSizeCylinder;

  /* user input variables */
  BOOLEAN uvar_help;
  INT4 uvar_ifo, uvar_blocksRngMed;
  REAL8 uvar_peakThreshold;
  REAL8 uvar_alpha, uvar_delta, uvar_h0, uvar_f0;
  REAL8 uvar_psi, uvar_phi0, uvar_fdot, uvar_cosiota;
  CHAR *uvar_earthEphemeris=NULL;
  CHAR *uvar_sunEphemeris=NULL;
  CHAR *uvar_sftDir=NULL;
  CHAR *uvar_fnameout=NULL;


  /*  set up the default parameters  */
  lalDebugLevel = 0;

  nfSizeCylinder = NFSIZE;
  /* LALDebugLevel must be called before anything else */
  SUB( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  /* set other user input variables */
  uvar_help = FALSE;
  uvar_peakThreshold = THRESHOLD;
  uvar_ifo = IFO;
  uvar_blocksRngMed = BLOCKSRNGMED;

  /* set default pulsar parameters */
  uvar_h0 = H0;
  uvar_alpha = ALPHA;
  uvar_delta = DELTA;
  uvar_f0 =  F0;
  uvar_fdot = FDOT;
  uvar_psi = PSI;
  uvar_cosiota = COSIOTA;
  uvar_phi0 = PHI0;

  /* now set the default filenames */
  uvar_earthEphemeris = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);

  uvar_sftDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sftDir,SFTDIRECTORY);

  uvar_fnameout = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_fnameout, FILEOUT);

  /* register user input variables */
  SUB( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",            &uvar_help),            &status);  
  SUB( LALRegisterINTUserVar(    &status, "ifo",             'i', UVAR_OPTIONAL, "Detector GEO(1) LLO(2) LHO(3)", &uvar_ifo ),            &status);
  SUB( LALRegisterINTUserVar(    &status, "blocksRngMed",    'w', UVAR_OPTIONAL, "RngMed block size",             &uvar_blocksRngMed),    &status);
  SUB( LALRegisterREALUserVar(   &status, "peakThreshold",   't', UVAR_OPTIONAL, "Peak selection threshold",      &uvar_peakThreshold),   &status);
  SUB( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_OPTIONAL, "Earth Ephemeris file",          &uvar_earthEphemeris),  &status);
  SUB( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_OPTIONAL, "Sun Ephemeris file",            &uvar_sunEphemeris),    &status);
  SUB( LALRegisterSTRINGUserVar( &status, "sftDir",          'D', UVAR_OPTIONAL, "SFT Directory",                 &uvar_sftDir),          &status);
  SUB( LALRegisterSTRINGUserVar( &status, "fnameout",        'o', UVAR_OPTIONAL, "Output file prefix",            &uvar_fnameout),        &status);
  SUB( LALRegisterREALUserVar(   &status, "alpha",           'r', UVAR_OPTIONAL, "Right ascension",               &uvar_alpha),           &status);
  SUB( LALRegisterREALUserVar(   &status, "delta",           'l', UVAR_OPTIONAL, "Declination",                   &uvar_delta),           &status);
  SUB( LALRegisterREALUserVar(   &status, "h0",              'm', UVAR_OPTIONAL, "h0 to inject",                  &uvar_h0),              &status);
  SUB( LALRegisterREALUserVar(   &status, "f0",              'f', UVAR_OPTIONAL, "Start search frequency",        &uvar_f0),              &status);
  SUB( LALRegisterREALUserVar(   &status, "psi",             'p', UVAR_OPTIONAL, "Polarization angle",            &uvar_psi),             &status);
  SUB( LALRegisterREALUserVar(   &status, "phi0",            'P', UVAR_OPTIONAL, "Initial phase",                 &uvar_phi0),            &status);
  SUB( LALRegisterREALUserVar(   &status, "cosiota",         'c', UVAR_OPTIONAL, "Cosine of iota",                &uvar_cosiota),         &status);
  SUB( LALRegisterREALUserVar(   &status, "fdot",            'd', UVAR_OPTIONAL, "Spindown parameter",            &uvar_fdot),            &status);

  /* read all command line variables */
  SUB( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 
  
  /* write the log file */
  fnamelog = (CHAR *)LALMalloc( 512*sizeof(CHAR));
  strcpy(fnamelog, uvar_fnameout);
  strcat(fnamelog, "_log");
  /* open the log file for writing */
  if ((fpLog = fopen(fnamelog, "w")) == NULL) {
    fprintf(stderr, "Unable to open file %s for writing\n", fnamelog);
    LALFree(fnamelog);
    exit(1);
  }

  /* get the log string */
  SUB( LALUserVarGetLog(&status, &logstr, UVAR_LOGFMT_CFGFILE), &status);  

  fprintf( fpLog, "## Log file for HoughMismatch\n\n");
  fprintf( fpLog, "# User Input:\n");
  fprintf( fpLog, "#-------------------------------------------\n");
  fprintf( fpLog, logstr);
  LALFree(logstr);

  /* append an ident-string defining the exact CVS-version of the code used */
  {
    CHAR command[1024] = "";
    fprintf (fpLog, "\n\n# CVS-versions of executable:\n");
    fprintf (fpLog, "# -----------------------------------------\n");
    fclose (fpLog);
    
    sprintf (command, "ident %s | sort -u >> %s", argv[0], fnamelog);
    system (command);	/* we don't check this. If it fails, we assume that */
    			/* one of the system-commands was not available, and */
    			/* therefore the CVS-versions will not be logged */

    LALFree(fnamelog); 
  }





  
  INFO( DRIVEHOUGHCOLOR_MSGENORM );
  return DRIVEHOUGHCOLOR_ENORM;

}
