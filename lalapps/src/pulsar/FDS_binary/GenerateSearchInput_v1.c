/*********************************************************************************/
/*    Program to generate an input file for lalapps_ComputeFStatisticBinary      */
/*                                                                               */
/*			           C. Messenger                                  */
/*                                                                               */
/*                         BIRMINGHAM UNIVERISTY -  2004                         */
/*********************************************************************************/

#include "GenerateBinaryMesh_v1.h"

INT4 lalDebugLevel=3;
static LALStatus status;

REAL8 alpha,delta;
REAL8 sma,period,ecc,argp;
INT4 tperisec,tperins;
REAL8 fmin,band,fres;
CHAR datadir[256],efiles[56],basename[256],yr[256],ifo[256],bintempfile[256],outfile[256],fout[256];
BOOLEAN signalflag;
BOOLEAN estimflag;
BOOLEAN binflag;
REAL8 doppler,thresh;
INT4 dterms,window;
REAL8 tspan;
INT4 tstart;

extern char *optarg;
extern int optind, opterr, optopt;
int ReadCommandLine(int argc,char *argv[]);
int OutputConfigFile();
int OutputBinTemplateFile();

int main(int argc,char *argv[]) 
{
 
  if (ReadCommandLine(argc,argv)) return 1;
  
  if (OutputConfigFile()) return 2;

  if (binflag) {
    if (OutputBinTemplateFile()) return 3;
  }
    
  return 0;

}

/*******************************************************************************/

int OutputBinTemplateFile() 
{

  BinaryMeshFileHeader BMFheader;
  FILE *btfp;
  INT4 ntemplate = 1;
  
  btfp=fopen(bintempfile,"w");
  if (btfp==NULL) {
    fprintf(stderr,"Unable to open file %s\n",bintempfile);
    return 1;
  }

  /* setup the header input */
  BMFheader.fmax=fmin+band;
  BMFheader.tspan=tspan;
  BMFheader.tstart.gpsSeconds=tstart;
  BMFheader.tstart.gpsNanoSeconds=0;
  BMFheader.Nfilters=ntemplate;
  BMFheader.mismatch=0.0;
  BMFheader.sma_0=sma;
  BMFheader.sma_MIN=sma;
  BMFheader.sma_MAX=sma;
  BMFheader.tperi_0.gpsSeconds=tperisec;
  BMFheader.tperi_0.gpsNanoSeconds=tperins;
  BMFheader.tperi_MIN.gpsSeconds=tperisec;
  BMFheader.tperi_MIN.gpsNanoSeconds=tperins;
  BMFheader.tperi_MAX.gpsSeconds=tperisec;
  BMFheader.tperi_MAX.gpsNanoSeconds=tperins;
  BMFheader.ecc_MIN=ecc;
  BMFheader.ecc_MAX=ecc;
  BMFheader.argp_MIN=argp;
  BMFheader.argp_MAX=argp;
  BMFheader.period_MIN=period;
  BMFheader.period_MAX=period;
  BMFheader.metric_XX=1.0;     /* set the metric components to 1 */
  BMFheader.metric_XY=1.0;
  BMFheader.metric_YY=1.0;
  sprintf(BMFheader.version,"v1_single");
  sprintf(BMFheader.det,ifo);
  BMFheader.RA=alpha;
  BMFheader.dec=delta;

  if (WriteMeshFileHeader(btfp,&BMFheader)) return 1;

  /* output single filter */
  fprintf(btfp,"%6.12f\t%6.12f\t%d\t%d\t%6.12f\t%6.12f\n",sma,period,tperisec,tperins,ecc,argp);

  fclose(btfp);

  return 0;

}

/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[]) 
{
  INT4 c, errflg = 0;
  char *temp;
  optarg = NULL;
  
   /* Initialize default values */
  alpha=0.0; /* a */
  delta=0.0; /* d */
  sprintf(ifo,"LLO"); /* I */            
  fmin=600.0; /* F */
  band=0.0; /* b */
  fres=0.0; /* r*/

  window=60; /* w */
  doppler=4e-4; /* p */
  dterms=16; /* k */
  estimflag=0; /* m */
  signalflag=1; /* S */
  binflag=0;   /* determines wether we output a binary template file or not */
  thresh=0.0; /* T */
  tstart=0; /* q */
  tspan=0; /* Q */

  sprintf(efiles,"./");/* g */  
  sprintf(yr,"00-04"); /* k */
  sprintf(datadir," "); /* v */
  sprintf(basename," "); /* u */
  sprintf(fout,"Fout.data"); /* o */
  sprintf(bintempfile,"bintemplate.data"); /* t */

  sma=0.0; /* A */
  ecc=0.0; /* e */
  tperisec=0; /* Q */
  tperins=0; /* q */
  period=0.0; /* P */
  argp=0.0; /* A */

  {
    int option_index = 0;
    static struct option long_options[] = {
      {"alpha", required_argument, 0, 'a'},
      {"delta", required_argument, 0, 'd'},
      {"ifo", required_argument, 0, 'I'},
      {"fmin", required_argument, 0, 'f'},
      {"band", required_argument, 0, 'b'},
      {"fres", required_argument, 0, 'r'},
      {"doppler", required_argument, 0, 'p'},
      {"dterms", required_argument, 0, 'k'},
      {"ephem", required_argument, 0, 'E'},
      {"yr", required_argument, 0, 'y'},
      {"datadir", required_argument, 0, 'D'},
      {"basename", required_argument, 0, 'B'},
      {"estim", no_argument, 0, 'm'},
      {"signal", no_argument, 0, 'S'},
      {"thresh", required_argument, 0, 'T'},
      {"bintempfile", required_argument, 0, 't'},
      {"doppler", required_argument, 0, 'p'},
      {"fout", required_argument, 0, 'O'},
      {"window", required_argument, 0, 'w'},
      {"smaxis", required_argument, 0, 'A'},
      {"period", required_argument, 0, 'P'},
      {"tperisec", required_argument, 0, 'X'},
      {"tperinan", required_argument, 0, 'x'},
      {"ecc", required_argument, 0, 'e'},
      {"argp", required_argument, 0, 'u'},
      {"tstart", required_argument, 0, 'q'},
      {"tspan", required_argument, 0, 'Q'},
      {"outfile", required_argument, 0, 'F'},
      {"help", required_argument, 0, 'h'}
    };
   
  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt_long (argc, argv,"ha:d:I:f:b:r:k:E:y:D:B:mST:t:p:O:w:A:P:X:x:q:Q:e:u:",long_options, &option_index)))!=-1)
    switch (c) {
      
    case 'a':
      alpha=atof(optarg);
      break;
    case 'd':
      delta=atof(optarg);
      break;
    case 'I':
      temp=optarg;
      sprintf(ifo,temp);
      break;
    case 'f':
      fmin=atof(optarg);
      break;
    case 'b':
      band=atof(optarg);
      break;
    case 'r':
      fres=atof(optarg);
      break;
    case 'k':
      dterms=atoi(optarg);
      break;
    case 'p':
      doppler=atof(optarg);
      break;
    case 'E':
      temp=optarg;
      sprintf(efiles,temp);
      break;
    case 'y':
      temp=optarg;
      sprintf(yr,temp);
      break;
    case 'D':
      temp=optarg;
      sprintf(datadir,temp);
      break;
    case 'B':
      temp=optarg;
      sprintf(basename,temp);
      break;
    case 'm':
      estimflag=1;
      break;
    case 'S':
      signalflag=1;
      break;
    case 'T':
      thresh=atof(optarg);
      break;
    case 't':
      temp=optarg;
      sprintf(bintempfile,temp);
      break;
    case 'w':
      window=atoi(optarg);
      break;
    case 'A':
      sma=atof(optarg);
      binflag=1;
      break;
    case 'P':
      period=atof(optarg);
      binflag=1;
      break;
    case 'X':
      tperisec=atoi(optarg);
      binflag=1;
      break;
    case 'x':
      tperins=atoi(optarg);
      binflag=1;
      break;
    case 'q':
      tstart=atoi(optarg);
      break;
    case 'Q':
      tspan=atof(optarg);
      break;
    case 'e':
      ecc=atof(optarg);
      binflag=1;
      break;
    case 'u':
      argp=atof(optarg);
      binflag=1;
      break;
    case 'O':
      temp=optarg;
      sprintf(fout,temp);
      break;
    case 'F':
      temp=optarg;
      sprintf(outfile,temp);
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\t--alpha       FLOAT\t Sky position alpha (equatorial coordinates) in radians [DEFAULT=0.0 ]\n");
      fprintf(stdout,"\t--delta       FLOAT\t Sky position delta (equatorial coordinates) in radians [DEFAULT=0.0 ]\n");
      fprintf(stdout,"\t--ifo         STRING\t Detector (LLO,LHO,GEO,VIRGO,TAMA) [DEFAULT=LLO]\n");
      fprintf(stdout,"\t--fmin        FLOAT\t Minimum search frequency in Hz [DEFAULT=0.0] \n");
      fprintf(stdout,"\t--band        FLOAT\t Bandwidth to be searched in Hz [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--fres        FLOAT\t Frequency resolution to be used in Hz [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--ephem       STRING\t Location of ephemeris data [DEFAULT=./]\n");
      fprintf(stdout,"\t--yr          STRING\t Year of ephemeris to obe read [DEFAULT=00-04]\n");
      fprintf(stdout,"\t--datadir     STRING\t Directory containing data to be searhed [DEFAULT=NULL]\n");
      fprintf(stdout,"\t--basename    STRING\t Location and basename of output SFT's [DEFAULT=NULL]\n");
      fprintf(stdout,"\t--estim       BOOLEAN\t Set if estimated parameters required [DEFAULT=FALSE]\n");
      fprintf(stdout,"\t--signal      BOOLEAN\t Set if data has only signal ppresent [DEFAULT=TRUE]\n");
      fprintf(stdout,"\t--window      INTEGER\t Window size to used in running median [DEFAULT=60]\n");
      fprintf(stdout,"\t--dterms      INTEGER\t Number of terms to use in the Dirichelet kernel [DEFAULT=16]\n");
      fprintf(stdout,"\t--bintempfile STRING\t Name of file to contain binary template [DEFAULT=bintemplate.data]\n");
      fprintf(stdout,"\t--thresh      FLOAT\t The threshold set on the F statistic [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--smaxis      FLOAT\t Projected semi-major axis of orbit in seconds [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--period      FLOAT\t Period of orbit in seconds [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--tperisec    INTEGER\t Observed time of periapse passage in seconds [DEFAULT=0] \n");
      fprintf(stdout,"\t--tperinan    INTEGER\t Observed time of periapse passage in nanoseconds [DEFAULT=0]\n");
      fprintf(stdout,"\t--tstart      INTEGER\t The start time of the observation [DEFAULT=0] \n");
      fprintf(stdout,"\t--tspan       FLOAT\t The span of the observation in seconds [DEFAULT=0]\n");
      fprintf(stdout,"\t--ecc         FLOAT\t Orbital eccentricity [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--argp        FLOAT\t Argument of orbital periapse in radians [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--fout        STRING\t Name of file to output Fstat results to [DEFAULT=NULL]\n");
      fprintf(stdout,"\t--outfile     STRING\t Name of file to output input search file [DEFAULT=NULL]\n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }

  }
  
  /* update global variable and return */
  return errflg;
}

/*******************************************************************************/

int OutputConfigFile() 
{
  FILE *fp;

  /*   %strcpy(filename,inDataFilename); */
  fp=fopen(outfile,"w");
  if (fp==NULL) {
    fprintf(stderr,"Unable to open file %s\n",outfile);
    return 1;
  }

  fprintf(fp,"## settings generated by GenSearchInput.c \n");
  fprintf(fp,"## ---------------------------------------------------\n");
  fprintf(fp,"\n");
  fprintf(fp,"## ---------------------------------------------------\n");
  fprintf(fp,"## Designed so that if the binary flag is not set this\n");
  fprintf(fp,"## input file used with lalapps_ComputeFstatisticBinary\n");
  fprintf(fp,"## will produce the same results as\n");
  fprintf(fp,"## lalapps_ComputeFStatistic.\n");
  fprintf(fp,"## ---------------------------------------------------\n");
  fprintf(fp,"\n");
  fprintf(fp,"## ---------------------------------------------------\n");
  fprintf(fp,"## REQUIRED user variables\n");
  fprintf(fp,"## ---------------------------------------------------\n");
  fprintf(fp,"\n");
  fprintf(fp,"DataDir         = %s\t# Directory where SFT's are located\n",datadir);
  fprintf(fp,"EphemDir        = %s\t# Directory where Ephemeris files are located\n",efiles);
  fprintf(fp,"Freq            = %6.12f\t\t# Starting search frequency in Hz\n",fmin);
  fprintf(fp,"\n");
  fprintf(fp,"IFO             = %s\t\t# 0=GEO, 1=LLO, 2=LHO or 3=Roman Bar\n",ifo);
  fprintf(fp,"\n");
  fprintf(fp,"## need to specify exactly ONE of (Alpha,Delta) or skyRegion\n");
  fprintf(fp,"Alpha           = %6.12f\t# Right ascension (equatorial coordinates) in radians\n",alpha);
  fprintf(fp,"Delta           = %6.12f\t# Declination (equatorial coordinates) in radians\n",delta);
  fprintf(fp,"AlphaBand       = 0.0           # Band in alpha (equatorial coordinates) in radians\n");
  fprintf(fp,"DeltaBand       = 0.0           # Band in delta (equatorial coordinates) in radians\n");
  fprintf(fp,"\n");
  fprintf(fp,"#skyRegion      = ( 6.25, 0.5 ), (6.3, 0.5), (6.3, 0.51), (6.25, 0.51) # just an example\n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");
  fprintf(fp,"## -------------------------------------------------------------\n");
  fprintf(fp,"## OPTIONAL user variables, defaults are shown here (if any)\n");
  fprintf(fp,"## -------------------------------------------------------------\n");
  fprintf(fp,"\n");
  fprintf(fp,"## ---------- SEARCH-PARAMETERS:\n");
  fprintf(fp,"\n");
  fprintf(fp,"BaseName        = %s\t# The base name of the input file\n",basename);
  fprintf(fp,"EphemYear       = %s\t\t# Year (or range of years) of ephemeris files\n",yr);
  if (estimflag) {
    fprintf(fp,"EstimSigParam   = true         # Do Signal Parameter Estimation\n");
  }
  else {
    fprintf(fp,"EstimSigParam   = false         # Do Signal Parameter Estimation\n");
  }
  if (signalflag) {
    fprintf(fp,"SignalOnly      = true         # Signal only flag\n");
  }
  else{
    fprintf(fp,"SignalOnly      = false         # Signal only flag\n");
  }
  fprintf(fp,"\n");
  fprintf(fp,"Dterms          = %d\t\t# Number of Dirichlet kernels\n",dterms);
  fprintf(fp,"\n");
  fprintf(fp,"FreqBand        = %6.12f\t# Search frequency band in Hz\n",band);
  fprintf(fp,"dFreq           = %6.12f\t# Frequency resolution in Hz (default: 1/(2*Tsft*Nsft)\n",fres);
  fprintf(fp,"\n");
  fprintf(fp,"#dAlpha         = 0.001         # Resolution in alpha (rad) for gridType=\"flat\"|\"isotropic\"\n");
  fprintf(fp,"                                # in the case of the isotropic-grid this is the value at the equator\n");
  fprintf(fp,"#dDelta         = 0.001         # Resolution in delta (rad) for gridType=\"flat\"|\"isotropic\"\n");
  fprintf(fp,"\n");
  fprintf(fp,"#f1dot          = 0.0           # First spindown parameter f1dot\n");
  fprintf(fp,"#f1dotBand      = 0.0           # Search-band for f1dot\n");
  fprintf(fp,"#df1dot         = 0.0           # Resolution for f1dot (default 1/(2*Tobs*Tsft*Nsft)\n");
  fprintf(fp,"\n");
  fprintf(fp,"Fthreshold      = %lf\t# Signal Set the threshold for selection of 2F\n",thresh);
  fprintf(fp,"\n");
  fprintf(fp,"## ---------- BINARY-SEARCH PARAMETERS:\n");
  fprintf(fp,"\n");
  fprintf(fp,"binary                  = true                          # Binary search flag\n");
  fprintf(fp,"binarytemplatefile      = %s\t# Name and location of binary template file\n",bintempfile);
  fprintf(fp,"dopplermax              = %lf\t# Maximum safe doppler shift factor from binary\n",doppler);
  fprintf(fp,"windowsize              = %d\t\t# Window size used in estimation of noise floor\n",window);
  fprintf(fp,"\n");
  fprintf(fp,"## ---------- TEMPLATE-GRID PARAMETERS:\n");
  fprintf(fp,"\n");
  fprintf(fp,"#gridType       = 0             # Template grid: 0=flat, 1=isotropic, 2=metric, 3=file\n");
  fprintf(fp,"#metricType     = 1             # Metric: 0=none, 1=Ptole-analytic, 2=Ptole-numeric, 3=exact\n");
  fprintf(fp,"#metricMismatch = 0.02          # Maximal mismatch for metric tiling\n");
  fprintf(fp,"#skyGridFile    = franken.grid  # file containing a sky-grid for gridType=\"file\"\n");
  fprintf(fp,"\n");
  fprintf(fp,"## ---------- OUTPUT PARAMETERS\n");
  fprintf(fp,"\n");
  fprintf(fp,"#outputLabel    = label         # Label to be appended to all output file-names\n");
  fprintf(fp,"outputFstat     = %s\t# filename for output of calculated F-statistic field over parameters\n",fout);
  fprintf(fp,"#outputSkyGrid  = southernHemi.grid     # output sky-grid into this file\n");
  fprintf(fp,"#openDX         = false         # make output-files openDX readable (only outputFstat so far)\n");
  fprintf(fp,"\n");   
  fprintf(fp,"## ---------- BOINC-specific PARAMETERS\n");
  fprintf(fp,"\n");
  fprintf(fp,"#mergedSFTFile  = NULL          #    Merged SFT's file to be used with BOINC\n");


  fclose(fp);
  return 0;
  
}
/*******************************************************************************/

