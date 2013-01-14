/* Nested Sampler Using LAL bayesian framework
 (C) John Veitch 2009

 */

#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/Units.h>
#include "LALInspiralMCMC.h"
#include "LALInspiralMCMCUser.h"
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/Random.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALDetectors.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALNoiseModels.h>
#include <lal/Date.h>
#include <lal/LALInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/FrequencySeries.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeSeries.h>
#include <lal/VectorOps.h>
#include <LALAppsVCSInfo.h>
#include <lalapps.h>
#include <fftw3.h>

#include "nest_calc.h"

#define MAXSTR 128
#define TIMESLIDE 10 /* Length of time to slide data to lose coherency */
#define DEBUG 0
#define USAGE "lalapps_inspnest ARGUMENTS [OPTIONS]\n \
Necessary ARGUMENTS:\n \
-o outfile\t:\tOutput samples to outfile\n \
--length duration\t:\tUse duration seconds of data to compute PSD\n \
--Nsegs INT\t:\tNumber of data segments for PSd estimation\n \
-I IFO\t:\tSpecify interferometer, one of H1, H1, L1, V1, or G1\n \
-C STRING\t:\tSpecify reading data from frame channel STRING\n \
-i cachefile\t:\tRead data from LIGO cache file cachefile.\n \
\tif cachefile is LALLIGO, LAL2kLIGO, LALGEO, LALVirgo, LALAdLIGO or LALEGO\n \
\tfake noise will be generated using the approprate noise curve.\n \
\tUse more [... -i FC -I IFO -C Channel] for as many data sources as desired\n \
\n\n\tYou must specify one of the following trigger types\n \
[--XMLfile PATH\t:\tRead SnglInspiralTable from PATH]\n \
[--inj PATH\t:\tRead SimInspiralTable from PATH and perform injection (Use [-F] to fake injection)]\n \
[--end_time GPSTIME\t:\tSpecify end time prior centred at GPSTIME]\n \
 \n\n \
Optional OPTIONS:\n \
[--Nlive INT (1000)\t:\tNumber of live points in nested sampler]\n \
[--Nmcmc INT (100)\t:\tNumber of MCMC points in chain for each sample]\n \
[--Nruns INT (1)\t:\tRun INT parallel samplings of the shrinking distribution\n \
[--seed INT\t:\tSpecify nested sampling random seed, default will use date]\n \
[--dataseed INT\t:\t Seed for faking data]\n \
[-v, --verbose\t:\tProduce statistics while running]\n \
[--GPSstart datastart\t:\tStart PSD estimation from time datastart, will guess if not specified]\n \
[--srate rate (4096)\t:\tDownsample data to rate Hz]\n \
[--pad padding (1s)\t:\tPadding for PSD Tukey window\n \
[--event INT (0)\t:\tUse event INT from Sim or Sngl InspiralTable]\n \
[--mmin FLOAT, --mmax FLOAT\t:\tSpecify min and max prior chirp masses\n \
[--dmin FLOAT (1), --dmax FLOAT (100)\t:\tSpecify min and max prior distances in Mpc\n \
[--approximant STRING (TaylorF2)\t:\tUse a different approximant where STRING is (TaylorF2|TaylorT2|TaylorT3|TaylorT4|AmpCorPPN|IMRPhenomFA|IMRPhenomFB|IMRPhenomFB_NS|IMRPhenomFB_Chi|EOBNR|SpinTaylor|IMRPhenomFB_Chi_low)|IMRPhenomB_Chi_low]\n \
[--amporder INT\t:\tAmplitude order to use, requires --approximant AmpCorPPN]\n \
[--phaseorder INT\t:\tPhase PN order to use, multiply by two, i.e. 3.5PN=7. (Default 4 = 2.0PN)]\n\
[--H1GPSshift FLOAT\t: Specify timeslide in H1]\n \
[--L1GPSshift FLOAT\t: Specify timeslide in L1]\n \
[--V1GPSshift FLOAT\t: Specify timeslide in V1]\n \
[--timeslide\t:\tTimeslide data]\n[--studentt\t:\tuse student-t likelihood function]\n \
[--ra FLOAT --dec FLOAT\t:\tSpecify fixed RA and dec to use (DEGREES)]\n \
[--grb\t:\tuse GRB prior ]\n[--decohere offset\t:\tOffset injection in each IFO]\n \
[--dt FLOAT (0.01)\t:\ttime window (0.01s)]\n \
[--injSNR FLOAT\t:\tScale injection to have network SNR of FLOAT]\n \
[--SNRfac FLOAT\t:\tScale injection SNR by a factor FLOAT]\n \
[--pinparams STRING\t:\tList parameters to be fixed to their injected values (, separated) i.e. --pinparams mchirp,longitude\n \
[--version\t:\tPrint version information and exit]\n \
[--datadump DATA.txt\t:\tOutput frequency domain PSD and data segment to DATA.txt]\n \
[--flow NUM\t:\t:Set low frequency cutoff (default 40Hz)]\n\
[--chimin NUM\t:\tMin value of chi spin parameter]\n\
[--etamin NUM\t:\tMinimum value of eta]\n\
[--chimax NUM\t:\tMax value of chi spin parameter]\n\
[--snrpath PATH\t:\tOutput SNRs to a file in PATH]\n\
[--skyloc\t:\tUse Sky localisation prior for Larry's injections]\n\
[--m1m2\t:\tUse m1 and m2 as mass parameters instead of chirp mass and eta]\n\
\n\n \
Optional PhenSpinTaylorRD_template OPTIONS:\n \
[--onespin_flag INT\t:\tSet S2=(0,0,0) in PhenSpinTaylorRD template waveform]\n \
[--nospin_flag INT\t:\tSet S1=S2=(0,0,0) in PhenSpinTaylorRD template waveform]\n \
[--m_tot_min FLOAT\t:\tSet lower limit on total mass for PhenSpinTaylorRD template waveform. Default is 2 solar masses]\n \
[--m_tot_max FLOAT\t:\tSet upper limit on total mass for PhenSpinTaylorRD template waveform. Default is 35 solar masses]\n \
[--d_min FLOAT\t:\tSet lower limit on source distance for PhenSpinTaylorRD template waveform. Default is 1Mpc.]\n \
[--d_max FLOAT\t:\tSet upper limit on source distance for PhenSpinTaylorRD template waveform. Default is 100Mpc.]\n \
[--s1_mag_min FLOAT\t:\tSet lower limit on spin magnitude of body 1 for PhenSpinTaylorRD template waveform. Default is 0.]\n \
[--s1_mag_max FLOAT\t:\tSet upper limit on spin magnitude of body 1 for PhenSpinTaylorRD template waveform. Default is 1.]\n \
[--s1_theta_min FLOAT\t:\tSet lower limit on spin polar angle for body 1 for PhenSpinTaylorRD template waveform. Default is 0.]\n \
[--s1_theta_max FLOAT\t:\tSet upper limit on spin polar angle for body 1 for PhenSpinTaylorRD template waveform. Default is PI.]\n \
[--s2_mag_min FLOAT\t:\tSet lower limit on spin magnitude of body 2 for PhenSpinTaylorRD template waveform. Default is 0.]\n \
[--s2_mag_max FLOAT\t:\tSet upper limit on spin magnitude of body 2 for PhenSpinTaylorRD template waveform. Default is 1.]\n \
[--s2_theta_min FLOAT\t:\tSet lower limit on spin polar angle for body 2 for PhenSpinTaylorRD template waveform. Default is 0.]\n \
[--s2_theta_max FLOAT\t:\tSet upper limit on spin polar angle for body 2 for PhenSpinTaylorRD template waveform. Default is PI.]\n \
[--s_phi_min FLOAT\t:\tSet lower limit on spin1-spin2 azimuthal angle for PhenSpinTaylorRD template waveform. Default is 0.]\n \
[--s_phi_max FLOAT\t:\tSet upper limit on spin1-spin2 azimuthal angle for PhenSpinTaylorRD template waveform. Default is 2PI.]\n \
[--lat_min FLOAT\t:\tSet lower limit on source plane declination for PhenSpinTaylorRD template waveform. Default is -PI/2.]\n \
[--lat_max FLOAT\t:\tSet upper limit on source plane declination for PhenSpinTaylorRD template waveform. Default is PI/2.]\n \
[--long_min FLOAT\t:\tSet lower limit on source plane RA for PhenSpinTaylorRD template waveform. Default is 0.]\n \
[--long_max FLOAT\t:\tSet upper limit on source plane RA for PhenSpinTaylorRD template waveform. Default is 2PI.]\n \
[--iota_min FLOAT\t:\tSet lower limit on source plane inclination for PhenSpinTaylorRD template waveform. Default is 0.]\n \
[--iota_max FLOAT\t:\tSet upper limit on source plane inclination for PhenSpinTaylorRD template waveform. Default is PI.]\n \
[--help\t:\tPrint this message]\n \
"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

extern CHAR outfile[FILENAME_MAX];
CHAR *datadump=NULL;
extern double etawindow;
extern double timewindow;
double etamin=0.03;
CHAR **CacheFileNames = NULL;
CHAR **ChannelNames = NULL;
CHAR **IFOnames = NULL;
CHAR UserChannel[512];
CHAR **UserChannelNames = NULL;
CHAR *SNRpath = NULL;
int nChannel=0;
UINT4 nIFO=0;
int fakeinj =0;
int skylocprior=0;
REAL8 duration=0;
LIGOTimeGPS datastart;
INT4 SampleRate=0;
REAL8 minFreq=48.0;
REAL4 padding=1.0;
INT4 nSegs=0;
INT4 Nruns=1;
INT4 dataseed=0;
REAL4 fLow=40.0; /* Low-frequency cutoff */
UINT4 Nlive=1000;
CHAR *inputXMLFile;
CHAR *injXMLFile=NULL;
CHAR approx[128]="TaylorF2";
UINT4 event=0;
REAL8 manual_end_time=0;
REAL8 manual_mass_low=2.0;
REAL8 manual_mass_high=35.0;
REAL8 manual_RA=-4200.0;
REAL8 manual_dec=-4200.0;
REAL8 manual_dist_max=100.0;
REAL8 manual_dist_min=1.0;
int Nmcmc = 100;
double injSNR=-1.0;
extern INT4 seed;
int NINJA=0;
int verbose=0;
int timeslides=0;
int specifictimeslides=0;
int studentt=0;
int estimatenoise=1;
int SkyPatch=0;
int FakeFlag=0;
int GRBflag=0;
int SkyLocFlag=0;
REAL8 SNRfac=1.0;
REAL4 H1GPSshift = 0.0, L1GPSshift = 0.0, V1GPSshift = 0.0;
int HighMassFlag=0;
int decohereflag=0;
int onespin_flag=0;
int nospin_flag=0;
REAL8 offset=0.0;
extern const LALUnit strainPerCount;
INT4 ampOrder=0;
INT4 phaseOrder=4;
char *pinned_params=NULL;
UINT4 fLowFlag=0;
UINT4 m1m2Flag=0;

REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length);
int checkParamInList(const char *list, const char *param);

/* variables for NestInitManualPhenSpinRD_manual */
double compmassmin=1.;
double m_tot_min=2.;
double m_tot_max=35.;
double m_tot_min_highmass=30.;
double m_tot_max_highmass=300.;
double d_min=1.;
double d_max=100.;
double s1_mag_min=0.;
double s1_mag_max=0.95;
double s2_mag_min=0.;
double s2_mag_max=0.95;
double s1_theta_min=0.;
double s2_theta_min=0.;
double s_phi_min=0.;
double s1_theta_max=LAL_PI;
double s2_theta_max=LAL_PI;
double s_phi_max=LAL_TWOPI;
double iota_min=0.;
double iota_max=LAL_PI;
double long_min=0.;
double long_max=LAL_TWOPI;
double lat_min=-LAL_PI/2.;
double lat_max=LAL_PI/2.;
double manual_chi_min=-1.;
double manual_chi_max=1.;
int mc_flag=0;
/* */
void NestInitManual(LALMCMCParameter *parameter, void *iT);
void NestInitManualIMRB(LALMCMCParameter *parameter, void *iT);
void NestInitManualIMRBChi(LALMCMCParameter *parameter, void *iT);
void NestInitNINJAManual(LALMCMCParameter *parameter, void *iT);
void NestInitSkyPatch(LALMCMCParameter *parameter, void *iT);
void NestInitGRB(LALMCMCParameter *parameter, void *iT);
void NestInitSkyLoc(LALMCMCParameter *parameter, void *iT);
void NestInitInj(LALMCMCParameter *parameter, void *iT);
void NestInitManualPhenSpinRD(LALMCMCParameter *parameter, void *iT);
void initialise(int argc, char *argv[]);
void PrintSNRsToFile(REAL8* SNRs,SimInspiralTable *inj_table,LALMCMCInput *inputMCMC);

REAL8TimeSeries *readTseries(CHAR *cachefile, CHAR *channel, LIGOTimeGPS start, REAL8 length)
{
	LALStatus status;
	FrCache *cache = NULL;
	FrStream *stream = NULL;
	REAL8TimeSeries *out = NULL;
	memset(&status,0,sizeof(LALStatus));
    fprintf(stdout,"Attempting to open %s at time %lf\n",cachefile,start.gpsSeconds+1e-9*start.gpsNanoSeconds);
	cache  = XLALFrImportCache( cachefile );
	if(cache==NULL) {fprintf(stderr,"ERROR: Unable to import cache file %s\n",cachefile); exit(-1);}
	stream = XLALFrCacheOpen( cache );
	if(stream==NULL) {fprintf(stderr,"ERROR: Unable to open stream from frame cache file\n"); exit(-1);}
	out = XLALFrInputREAL8TimeSeries( stream, channel, &start, length , 0 );
	if(out==NULL) fprintf(stderr,"ERROR: unable to read channel %s from %s at time %i\nCheck the specified data duration is not too long\n",channel,cachefile,start.gpsSeconds);
	LAL_CALL(LALDestroyFrCache(&status,&cache),&status);
	LAL_CALL(LALFrClose(&status,&stream),&status);
	return out;
}

void initialise(int argc, char *argv[]){
	int i;
	int nCache=0; /* records the number of caches */
	int nifo=0;
	double GPS;
	/*	sprintf(outfile,"default.dat"); */
	/* Sets up global variables from the command line */
	static struct option long_options[]=
	{	{"cache",required_argument,0,'i'},
		{"seed",required_argument,0,'z'},
		{"dataseed",required_argument,0,'D'},
		{"GPSstart",required_argument,0,'G'},
		{"length",required_argument,0,'T'},
		{"srate",required_argument,0,'R'},
		{"pad",required_argument,0,'P'},
		{"Nsegs",required_argument,0,'S'},
		{"IFO",required_argument,0,'I'},
		{"Nlive",required_argument,0,'N'},
		{"XMLfile",required_argument,0,'X'},
		{"Nmcmc",required_argument,0,'M'},
		{"Nruns",required_argument,0,'r'},
		{"grb",no_argument,0,'b'},
		{"out",required_argument,0,'o'},
		{"inj",required_argument,0,'j'},
		{"fake",no_argument,0,'F'},
		{"injSNR",required_argument,0,'p'},
		{"dt",required_argument,0,'t'},
		{"event",required_argument,0,'E'},
		{"NINJA",no_argument,0,'n'},
		{"end_time",required_argument,0,'Z'},
		{"mmin",required_argument,0,'m'},
		{"mmax",required_argument,0,'g'},
		{"etamin",required_argument,0,'e'},
		{"verbose",no_argument,0,'v'},
		{"approximant",required_argument,0,'A'},
		{"timeslide",no_argument,0,'L'},
		{"H1GPSshift",required_argument,0,31},
		{"L1GPSshift",required_argument,0,32},
		{"V1GPSshift",required_argument,0,33},
		{"studentt",no_argument,0,'l'},
		{"ra",required_argument,0,'O'},
		{"dec",required_argument,0,'a'},
		{"SNRfac",required_argument,0,14},
		{"channel",required_argument,0,'C'},
		{"highmass",no_argument,0,15},
		{"decohere",required_argument,0,16},
		{"amporder",required_argument,0,17},
		{"phaseorder",required_argument,0,20},
		{"dmin",required_argument,0,18},
		{"dmax",required_argument,0,19},
		{"version",no_argument,0,'V'},
		{"help",no_argument,0,'h'},
		{"pinparams",required_argument,0,21},
		{"datadump",required_argument,0,22},
		{"flow",required_argument,0,23},
		{"nospin",required_argument,0,25},
		{"onespin",required_argument,0,26},
		{"M_min",required_argument,0,40},
		{"M_max",required_argument,0,41},
		{"d_min",required_argument,0,42},
		{"d_max",required_argument,0,43},
		{"s1_mag_min",required_argument,0,44},
		{"s1_mag_max",required_argument,0,45},
		{"s2_mag_min",required_argument,0,46},
		{"s2_mag_max",required_argument,0,47},
		{"s1_theta_min",required_argument,0,48},
		{"s1_theta_max",required_argument,0,49},
		{"s2_theta_min",required_argument,0,50},
		{"s2_theta_max",required_argument,0,51},
		{"s_phi_min",required_argument,0,52},
		{"s_phi_max",required_argument,0,53},
		{"lat_min",required_argument,0,56},
		{"lat_max",required_argument,0,57},
		{"long_min",required_argument,0,58},
		{"long_max",required_argument,0,59},
		{"iota_min",required_argument,0,60},
		{"iota_max",required_argument,0,61},
		{"m_tot_min",required_argument,0,62},
		{"m_tot_max",required_argument,0,63},
		{"chimin",required_argument,0,64}, /* N.B. ASCII codes 65 - 90 and 97-122 are letters */
		{"chimax",required_argument,0,91},
		{"mc_flag",no_argument,0,100},
		{"snrpath",required_argument,0,123},
		{"skyloc",no_argument,0,13},
		{"m1m2",no_argument,0,24},
		{0,0,0,0}};

	if(argc<=1) {fprintf(stderr,USAGE); exit(-1);}
	while((i=getopt_long(argc,argv,"hi:D:G:T:R:g:m:z:P:C:S:I:N:t:X:O:a:M:o:j:e:Z:A:E:nlFVvb",long_options,&i))!=-1){ switch(i) {
		case 64:
			manual_chi_min=atof(optarg);
			break;
		case 24:
			m1m2Flag=1;
			break;
		case 91:
			manual_chi_max=atof(optarg);
			break;
		case 40:
			m_tot_min=atof(optarg);
			break;
		case 41:
			m_tot_max=atof(optarg);
			break;
		case 42:
			d_min=atof(optarg);
			break;
		case 43:
			d_max=atof(optarg);
			break;
		case 44:
			s1_mag_min=atof(optarg);
			if (s1_mag_min<0.) {fprintf(stderr,"ERROR: The minimum spin amplituade is 0, please set --s1_mag_min 0 or more.\n");exit(1);}
			break;
		case 45:
			s1_mag_max=atof(optarg);
            if (s1_mag_max>1.) {fprintf(stderr,"ERROR: The maxmimum spin amplituade is 1, please set --s1_mag_max 1 or less.\n");exit(1);}
			break;
		case 46:
			s2_mag_min=atof(optarg);
			if (s2_mag_min<0.) {fprintf(stderr,"ERROR: The minimum spin amplituade is 0, please set --s2_mag_min 0 or more.\n");exit(1);}
			break;
		case 47:
			s2_mag_max=atof(optarg);
            if (s2_mag_max>1.) {fprintf(stderr,"ERROR: The maxmimum spin amplituade is 1, please set --s2_mag_max 1 or less.\n");exit(1);}			
			break;
		case 48:
			s1_theta_min=atof(optarg);
			if (s1_theta_min<0.) {fprintf(stderr,"ERROR: The minimum spin polar angle is 0, please set --s1_theta_min 0 or more.\n");exit(1);}
			break;
		case 49:
			s1_theta_max=atof(optarg);
            if (s1_theta_max>LAL_PI) {fprintf(stderr,"ERROR: The maximum spin polar angle is pi, please set --s1_theta_max pi or less.\n");exit(1);}
			break;
		case 50:
            s2_theta_min=atof(optarg);
            if (s2_theta_min<0.) {fprintf(stderr,"ERROR: The minimum spin polar angle is 0, please set --s1_theta_min 0 or more.\n");exit(1);}
            break;
		case 51:
            s2_theta_max=atof(optarg);
            if (s2_theta_max>LAL_PI) {fprintf(stderr,"ERROR: The maximum spin polar angle is pi, please set --s1_theta_max pi or less.\n");exit(1);}
			break;
		case 52:
			s_phi_min=atof(optarg);
			if (s_phi_min<0.) {fprintf(stderr,"ERROR: The minimum spin1-spin2 azimuthal angle is 0, please set --s_phi_min 0 or more.\n");exit(1);}
			break;
		case 53:
            s_phi_min=atof(optarg);
            if (s_phi_min>LAL_TWOPI) {fprintf(stderr,"ERROR: The maximum spin1-spin2 azimuthal angle is 2pi, please set --s_phi_max 2pi or less.\n");exit(1);}
            break;
		case 56:
			lat_min=atof(optarg);
			if (lat_min<-LAL_PI/2.) {fprintf(stderr,"ERROR: The minimum latitude is -pi/2, please set --lat_min -pi/2 or more.\n");exit(1);}
			break;
		case 57:
			lat_max=atof(optarg);
			if (lat_max>LAL_PI/2.) {fprintf(stderr,"ERROR: The maximum latitude is pi/2, please set --lat_max pi/2 or less.\n");exit(1);}
			break;
		case 58:
			long_min=atof(optarg);
			if (long_min<0.) {fprintf(stderr,"ERROR: The minimum longitude is 0, please set --long_min 0 or more.\n");exit(1);}
			break;
		case 59:
			long_max=atof(optarg);
			if (long_max>LAL_TWOPI) {fprintf(stderr,"ERROR: The maximum longitude is 2pi, please set --long_max 2pi or less.\n");exit(1);}
			break;
		case 60:
			iota_min=atof(optarg);
			if (iota_min<0.) {fprintf(stderr,"ERROR: The minimum inclination is 0, please set --iota_min 0 or more.\n");exit(1);}
			break;
		case 61:
			iota_max=atof(optarg);
			if (iota_max>LAL_PI) {fprintf(stderr,"ERROR: The maximum inclination is pi, please set --iota_max pi or less.\n");exit(1);} 
			break;
		case 62:
			m_tot_min=atof(optarg);
			if (m_tot_min<0.) {fprintf(stderr,"ERROR: Masses must be positive.\n");exit(1);}
			break;
		case 63:
			m_tot_max=atof(optarg);
			if (m_tot_max>35.) {fprintf(stderr,"Warning: Highmass flag not set.\n");exit(1);}
			break;
		case 25:
			 nospin_flag=atof(optarg);
		   	 break;
		case 26:
			onespin_flag=atof(optarg);
			break;
		case 'h':
			fprintf(stdout,USAGE);
			exit(0);
			break;
		case 21:
			pinned_params=calloc(strlen(optarg)+1 ,sizeof(char));
			memcpy(pinned_params,optarg,strlen(optarg)+1);
			break;
		case 123:
			SNRpath = calloc(strlen(optarg)+1,sizeof(char));
			memcpy(SNRpath,optarg,strlen(optarg)+1);
			break;
		case 'V':
			fprintf(stdout,"LIGO/LSC Bayesian parameter estimation and evidence calculation code\nfor CBC signals, using nested sampling algorithm.\nJohn Veitch <john.veitch@ligo.org>\n");
			XLALOutputVersionString(stderr,0);
			exit(0);
			break;
		case 18:
			manual_dist_min=atof(optarg);
			break;
		case 19:
			manual_dist_max=atof(optarg);
			break;
		case 17:
			ampOrder=atoi(optarg);
			if(ampOrder>5) {fprintf(stderr,"ERROR: The maximum amplitude order is 5, please set --ampOrder 5 or less\n"); exit(1);}
			break;
		case 14:
			SNRfac=atof(optarg);
			break;
		case 31:
			H1GPSshift = atof(optarg);
			specifictimeslides=1;
			break;
		case 32:
			L1GPSshift = atof(optarg);
			specifictimeslides=1;
			break;
		case 33:
			V1GPSshift = atof(optarg);
			specifictimeslides=1;
			break;
		case 16:
			decohereflag=1;
			offset=atof(optarg);
			break;
		case 15:
			HighMassFlag=1;
			break;
		case 20:
			phaseOrder=atoi(optarg);
			break;
		case 'i': /* This type of arragement builds a list of file names for later use */
			if(nCache==0) CacheFileNames=malloc(sizeof(char *));
			else		CacheFileNames=realloc(CacheFileNames,(nCache+1)*sizeof(char *));
			CacheFileNames[nCache]=malloc(strlen(optarg)+1);
			strcpy(CacheFileNames[nCache++],optarg);
			break;
		case 'C':
			if(nChannel==0) UserChannelNames=malloc(sizeof(char *));
			else UserChannelNames=realloc(UserChannelNames,(nChannel+1)*sizeof(char *));
			UserChannelNames[nChannel]=malloc(strlen(optarg)+1);
			strcpy(UserChannelNames[nChannel++],optarg);
			break;
		case 'D':
			dataseed=atoi(optarg);
			break;
		case 'O':
			manual_RA=atof(optarg)*LAL_PI/180.0;
			SkyPatch=1;
			break;
		case 'b':
			GRBflag=1;
			break;
		case 'a':
			manual_dec=atof(optarg)*LAL_PI/180.0;
			SkyPatch=1;
			break;
		case 'A':
			strncpy(approx,optarg,128);
			break;
		case 'l':
			studentt=1;
			break;
		case 'v':
			verbose=1;
			break;
		case 'm':
			manual_mass_low=atof(optarg);
			printf("setting m_low=%e\n",manual_mass_low);
			break;
		case 'g':
			manual_mass_high=atof(optarg);
			printf("setting m_high=%e\n",manual_mass_high);
			break;
		case 't':
			timewindow=atof(optarg);
			break;
		case 'z':
			seed=atoi(optarg);
			break;
		case 'E':
			event=atoi(optarg);
			break;
		case 'p':
			injSNR=atof(optarg);
			break;
		case 'Z':
			manual_end_time=atof(optarg);
			break;
		case 'e':
			etamin=atof(optarg);
			break;
		case 'r':
			Nruns=atoi(optarg);
			break;
		case 'F':
			fakeinj=1;
			break;
		case 'S':
			nSegs=atoi(optarg);
			break;
		case 'M':
			Nmcmc = atof(optarg);
			break;
		case 'j':
			injXMLFile=(CHAR *)malloc(strlen(optarg)+1);
			strcpy(injXMLFile,optarg);
			break;
		case 'X':
			inputXMLFile=(CHAR *)malloc(strlen(optarg)+1);
			strcpy(inputXMLFile,optarg);
			break;
		case 22:
			datadump=(CHAR *)malloc(strlen(optarg)+1);
			strcpy(datadump,optarg);
			break;
		case 'N':
			Nlive=atoi(optarg);
			break;
		case 'I':
			if(nifo==0) {IFOnames=malloc(sizeof(char **)); ChannelNames=malloc(sizeof(char **));}
			else	{IFOnames=realloc(IFOnames,(nifo+1)*sizeof(CHAR **)); ChannelNames=realloc(ChannelNames,(nChannel+1)*sizeof(char **));}
			IFOnames[nifo]=malloc(strlen(optarg)+1);
            printf("strlen(optarg)=%zu, optarg=%s\n",strlen(optarg),optarg);
			ChannelNames[nifo]=malloc(MAXSTR+1);
			/*strcpy(IFOnames[nifo],optarg);*/
            sprintf(IFOnames[nifo],"%s",optarg);
            nifo=nifo+1;
			break;
		case 'o':
			strcpy(outfile,optarg);
			break;
		case 'G':
			GPS=atof(optarg);
			XLALGPSSetREAL8(&datastart,GPS);
			break;
		case 'T':
			duration=atof(optarg);
			break;
		case 'R':
			SampleRate=atoi(optarg);
			break;
		case 'P':
			padding=atof(optarg);
			break;
		case 'n':
			NINJA=1;
			fLow=30.0;
			break;
		case 'L':
			timeslides=1;
			break;
		case 23:
			fLow=atof(optarg);
			fLowFlag=1;
			break;
		case 13:
			skylocprior=1;
			break;
		default:
			fprintf(stdout,USAGE); exit(0);
			break;
	}
	}

	if(inputXMLFile==NULL && injXMLFile==NULL && manual_end_time==0){fprintf(stderr,"Error, you must specify --inj or --XMLfile for trigger list\nOr --end_time, --dt, --Mmin and --Mmax for manual search"); exit(-1);}
	/* Check that the channel/cache combo adds up */
	if(nifo!=nCache || nCache==0) {fprintf(stderr,"Error: You must have equal numbers of IFOs and frame caches, and they must be paired in the correct order!\n");
	exit(-1); }
	if(nChannel>0 && nChannel!=nCache) {fprintf(stderr,"Error: You must specify a channel for each cache file\n"); exit(-1);}
	nIFO=nifo;
	/*	for(i=0;i<nIFO;i++) fprintf(stdout,"%s\t|%s\t| %s\n",IFOnames[i],CacheFileNames[i],ChannelNames[i]); */
	if(Nmcmc==0){fprintf(stderr,"Error: --Nmcmc not specified or zero, use >0\n"); exit(-1);}
	if(SampleRate==0){fprintf(stderr,"Error: --srate not specified. Using 4096 Hz which may NOT be what you want!\n"); SampleRate=4096;}
	if(nSegs==0){fprintf(stderr,"Error: --Nsegs must be greater than 0\n"); exit(-1);}
	if(Nlive<=1){fprintf(stderr,"Error: Nlive must be >1"); exit(-1);}
	if(studentt) estimatenoise=0;
	return;
}

/* =========================== MAIN ==================================== */

int main( int argc, char *argv[])
{
	static LALStatus status;
	LALMCMCParameter **Live = NULL; /* Structure which holds the parameters */
	LALMCMCInput	inputMCMC;
	LALMCMCInput	inputMCMC_N;
	REAL8TimeSeries *RawData;
	UINT4			seglen=0;
	SnglInspiralTable *inputCurrent = NULL;
	SimInspiralTable *injTable = NULL;
	INT4 numTmplts = 0;
	UINT4 i,j;
	REAL8FFTPlan *fwdplan = NULL;
	REAL8FFTPlan *revplan = NULL;
	REAL8Window  *windowplan = NULL;
	INT4 stride=0;
	/* REAL8 strideDur=0.0; */
	REAL8 evidence=0;
	INT4 UNUSED segnum=0;
	RandomParams *randparam=NULL;
	RandomParams *datarandparam=NULL;
	REAL4 TSoffset;
	LIGOTimeGPS realstart,segmentStart;
	REAL8 networkSNR=0.0;
    REAL8 * SNRs=NULL;
    SNRs=calloc(nIFO+1 ,sizeof(REAL8));
    
	lal_errhandler = LAL_ERR_EXIT;
        set_debug_level( "33" );

	seed=0;
	etawindow=1.0;
	timewindow=0.05;
	initialise(argc,argv); /* Get the arguments and act on them */
	if( timeslides && specifictimeslides ){
		fprintf( stderr, "Error: can not use both random and specific timeslides.\n");
		exit( 1 );
	}
	if(inputXMLFile!=NULL){
		/* read in the input file */
		numTmplts = LALSnglInspiralTableFromLIGOLw( &inputCurrent, inputXMLFile, 0, -1);
		if ( numTmplts < 0 )
		{
			fprintf( stderr, "Error: unable to read trigger %i from %s\n", event,inputXMLFile );
			exit( 1 );
		}
		i=0;
		while(i<event) {i++; inputCurrent = inputCurrent->next;}
	}
	REAL8 segDur = duration/(REAL8)nSegs;
	realstart=datastart;
	/* Number of sample in a segment */
	seglen=(UINT4)(segDur*SampleRate);
	/*	seglen=(INT4)pow(2.0,ceil(log2((REAL8)seglen)));*/  /* Make it a power of two for the FFT */
	segDur = seglen/SampleRate;
	nSegs =(INT4)floor(duration/segDur);

	fprintf(stderr,"Choosing %i segments length %i, (%f s)\n",nSegs,seglen,segDur);

	stride = seglen; /* Overlap the padding */
	/* strideDur = stride / SampleRate; */


	if(segDur<=2.0*padding){fprintf(stderr,"ERROR: Seg length %lf s too small for padding %lf s\n",segDur,padding);exit(-1);}
	if(segDur-2.0*padding<6.0){fprintf(stderr,"WARNING: using <6s segments (excl. padding) unadvisable, your current unpadded seglen is %lf s\n",segDur-2.0*padding);}

	int check=0;
	fwdplan = XLALCreateForwardREAL8FFTPlan( seglen, 0 );
	revplan = XLALCreateReverseREAL8FFTPlan( seglen, 0 );
	memset(&inputMCMC,0,sizeof(inputMCMC)); /* CLEAR THE INPUT STRUCTURE! */
	memset(&inputMCMC_N,0,sizeof(inputMCMC_N)); /* CLEAR THE INPUT STRUCTURE! */
	inputMCMC.deltaT=(REAL8 )(1.0/SampleRate);
	inputMCMC.verbose=verbose;
	char strainname[20]="LSC-STRAIN";
	if(NINJA) sprintf(strainname,"STRAIN"); /* Different strain channel name for NINJA */

	/* Make a copy of the detectors list */
	LALDetector *localCachedDetectors=calloc(LAL_NUM_DETECTORS,sizeof(LALDetector));
	for(i=0;i<LAL_NUM_DETECTORS;i++) memcpy(&(localCachedDetectors[i]),&lalCachedDetectors[i],sizeof(LALDetector));

	/* Set up Detector structures */
	for(i=0;i<nIFO;i++){
		if(!strcmp(IFOnames[i],"H1")) {
			inputMCMC.detector[i]=&localCachedDetectors[LALDetectorIndexLHODIFF];
			if(nChannel>0) sprintf(ChannelNames[i],"%s",UserChannelNames[i]);
			else sprintf((ChannelNames[i]),"H1:%s",strainname);
			continue;}
		if(!strcmp(IFOnames[i],"H2")) {
			inputMCMC.detector[i]=&localCachedDetectors[LALDetectorIndexLHODIFF];
			if (nChannel>0) sprintf(ChannelNames[i],"%s",UserChannelNames[i]);
			else sprintf((ChannelNames[i]),"H2:%s",strainname);
			continue;}
		if(!strcmp(IFOnames[i],"LLO")||!strcmp(IFOnames[i],"L1")) {
			inputMCMC.detector[i]=&localCachedDetectors[LALDetectorIndexLLODIFF];
			if (nChannel>0) sprintf(ChannelNames[i],"%s",UserChannelNames[i]);
			else sprintf((ChannelNames[i]),"L1:%s",strainname);
			continue;}
		if(!strcmp(IFOnames[i],"V1")||!strcmp(IFOnames[i],"VIRGO")) {
			inputMCMC.detector[i]=&localCachedDetectors[LALDetectorIndexVIRGODIFF];
			if(!NINJA) sprintf((ChannelNames[i]),"V1:h_16384Hz");
			else sprintf((ChannelNames[i]),"V1:STRAIN");
			if (nChannel>0) sprintf(ChannelNames[i],"%s",UserChannelNames[i]);
			continue;}
		if(!strcmp(IFOnames[i],"GEO")||!strcmp(IFOnames[i],"G1")) {
			inputMCMC.detector[i]=&localCachedDetectors[LALDetectorIndexGEO600DIFF];
			if (nChannel>0) sprintf(ChannelNames[i],"%s",UserChannelNames[i]);
			else sprintf((ChannelNames[i]),"G1:DER_DATA_H");
			continue;}
		inputMCMC_N.detector[i]=inputMCMC.detector[i];
		/*		if(!strcmp(IFOnames[i],"TAMA")||!strcmp(IFOnames[i],"T1")) {inputMCMC.detector[i]=&lalCachedDetectors[LALDetectorIndexTAMA300DIFF]; continue;}*/
		fprintf(stderr,"Unknown interferometer %s. Valid codes: H1 H2 L1 V1 GEO\n",IFOnames[i]); exit(-1);
	}

	inputMCMC.fLow = fLow;
	inputMCMC_N.fLow = inputMCMC.fLow;


	/* Prepare for injections */
	UINT4 Ninj=0;
	PPNParamStruc InjParams;
	LIGOTimeGPS injstart;
	memset(&injstart,0,sizeof(LIGOTimeGPS));
	memset(&InjParams,0,sizeof(PPNParamStruc));
	if(NULL!=injXMLFile) {Ninj=SimInspiralTableFromLIGOLw(&injTable,injXMLFile,0,0);
		if(Ninj<event) {fprintf(stderr,"Error reading event %i from %s\n",event,injXMLFile); exit(-1);}
		i=0;
		while(i<event) {i++; injTable = injTable->next;} /* Select event */
		if (!fLowFlag){
			if(injTable->f_lower>0.0) inputMCMC.fLow = injTable->f_lower;
			else {injTable->f_lower = inputMCMC.fLow;
				fprintf(stderr,"Warning, injection does not specify f_lower, using default %lf\n",inputMCMC.fLow);}
		}
//		InjParams.deltaT=1.0/SampleRate;
//		InjParams.fStartIn=(REAL4)inputMCMC.fLow;
//		memset(&InjectGW,0,sizeof(CoherentGW));
		fprintf(stderr,"Injected event %i:\tMass1: %lf\tMass2: %lf\n\tDistance: %lf Mpc\teta: %lf\n",event,injTable->mass1,injTable->mass2,injTable->distance,injTable->eta);
		/*		memcpy(&(InjParams.epoch),&(injTable->geocent_end_time),sizeof(LIGOTimeGPS)); */
//		Approximant injapprox;
//		fprintf(stderr,"INJ: end time = %lf\n",injTable->geocent_end_time.gpsSeconds + injTable->geocent_end_time.gpsNanoSeconds*1.e-9);
//		LALGetApproximantFromString(&status,injTable->waveform,&injapprox);
//		if(injapprox!=GeneratePPN) {fprintf(stderr,"WARNING!!!!! Not using GeneratePPN approximant may result in offset of the end time!\n");}
//		LALGenerateInspiral(&status,&InjectGW,injTable,&InjParams);
//		if(status.statusCode!=0) {fprintf(stderr,"Error generating injection!!!\n"); REPORTSTATUS(&status); }
		/****************************************************************************************************/
		/********** THIS IS ONLY NECESSARY WHILE THE LALGenerateInspiral and LALInspiralParameterCalc *******/
		/********** GIVE DIFFERENT CHIRP TIMES !                                                      *******/

//		insptemplate.totalMass=InjParams.mTot;
//		insptemplate.eta = InjParams.eta;
//		insptemplate.approximant = TaylorF2;
//		insptemplate.order = LAL_PNORDER_TWO;
//		insptemplate.fLower = inputMCMC.fLow;
//		insptemplate.massChoice = totalMassAndEta;
//		LALInspiralParameterCalc(&status,&insptemplate);
		/*InjParams.tc = insptemplate.tC;*/
//		fprintf(stderr,"GenerateInspiral chirp time=%lf, ParameterCalc chirp time = %lf\n",InjParams.tc,insptemplate.tC);
		/*****************************************************************************************************/

//		injstart = injTable->geocent_end_time;
//		XLALGPSAdd(&injstart, -InjParams.tc); /* makes injstart the time at fLow */
		/*		fprintf(stderr,"start time = %lf\n",injstart.gpsSeconds + injstart.gpsNanoSeconds*1.e-9); */
//		fprintf(stderr,"INJ: Injected wave chirp time: %lf s\n",InjParams.tc);
//		if(InjectGW.h) memcpy(&(InjectGW.h->epoch),&injstart,sizeof(LIGOTimeGPS));
//		if(InjectGW.a) memcpy(&(InjectGW.a->epoch),&injstart,sizeof(LIGOTimeGPS));
//		if(InjectGW.f) memcpy(&(InjectGW.f->epoch),&injstart,sizeof(LIGOTimeGPS));
//		if(InjectGW.phi) memcpy(&(InjectGW.phi->epoch),&injstart,sizeof(LIGOTimeGPS));
//		if(InjectGW.shift) memcpy(&(InjectGW.shift->epoch),&injstart,sizeof(LIGOTimeGPS));
	}

	/* Get the end time of the trigger or injection */
	int ETgpsSeconds,ETgpsNanoseconds;
	if(NULL!=injXMLFile) {
		ETgpsSeconds=injTable->geocent_end_time.gpsSeconds;
		ETgpsNanoseconds=injTable->geocent_end_time.gpsNanoSeconds;}
	else if(NULL!=inputXMLFile) {
		ETgpsSeconds = inputCurrent->end_time.gpsSeconds;
		ETgpsNanoseconds=inputCurrent->end_time.gpsNanoSeconds;}
	else {
		ETgpsNanoseconds = (INT4)1.e9*fmod(manual_end_time,1.0);
		ETgpsSeconds = (INT4) floor(manual_end_time);
	}

	/* If the data segment is not specified, try to guess it from the trigger time */
	
	if(datastart.gpsSeconds==0){
		fprintf(stderr,"You did not specify a data start time, attempt to guess.\nAdjusting GPSstart to %i for trigger %i\n",ETgpsSeconds-(INT4)duration/2,event);
		datastart.gpsSeconds=ETgpsSeconds-(INT4)duration/2;
		datastart.gpsNanoSeconds=0;
		realstart=datastart;
	}

	if(ETgpsSeconds>datastart.gpsSeconds+duration) {fprintf(stderr,"Error, trigger lies outwith data range %i - %i\n",datastart.gpsSeconds,datastart.gpsSeconds+(INT4)duration); exit(-1);}


	/* Read in the data for each IFO */
	for(i=0,j=0;i<nIFO;i++){
		INT4 TrigSegStart,TrigSample;
		inputMCMC.ifoID[i] = IFOnames[i];
		inputMCMC.deltaF = (REAL8)SampleRate/seglen;
		inputMCMC_N.ifoID[i] = inputMCMC.ifoID[i];
		inputMCMC_N.deltaF = inputMCMC.deltaF;
		datastart=realstart; /* Reset the datastart in case it has been slid previously */
		segmentStart = datastart;
		datarandparam=XLALCreateRandomParams(dataseed+(INT2)IFOnames[i][0]+(INT2)IFOnames[i][1]);

		/* Check for synthetic data */
		if(!(strcmp(CacheFileNames[i],"LALLIGO") && strcmp(CacheFileNames[i],"LALVirgo") && strcmp(CacheFileNames[i],"LALGEO") && strcmp(CacheFileNames[i],"LALEGO") && strcmp(CacheFileNames[i],"LALAdLIGO")))
		{
			typedef void (NoiseFunc)(LALStatus *status,REAL8 *psd,REAL8 f);
			NoiseFunc *PSD=NULL;
			FakeFlag=1;
			REAL8 scalefactor=1;
			/* Selection of the noise curve */
			if(!strcmp(CacheFileNames[i],"LALLIGO")) {PSD = &LALLIGOIPsd; scalefactor=9E-46;}
			if(!strcmp(CacheFileNames[i],"LALVirgo")) {PSD = &LALVIRGOPsd; scalefactor=1.0;}
			if(!strcmp(CacheFileNames[i],"LALGEO")) {PSD = &LALGEOPsd; scalefactor=1E-46;}
			if(!strcmp(CacheFileNames[i],"LALEGO")) {PSD = &LALEGOPsd; scalefactor=1.0;}
			if(!strcmp(CacheFileNames[i],"LALAdLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 1E-49;}
			if(!strcmp(CacheFileNames[i],"LAL2kLIGO")) {PSD = &LALAdvLIGOPsd; scalefactor = 36E-46;}
			if(PSD==NULL) {fprintf(stderr,"Error: unknown simulated PSD: %s\n",CacheFileNames[i]); exit(-1);}
			inputMCMC.invspec[i]=(REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("inverse spectrum",&realstart,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);				  inputMCMC_N.invspec[i]=(REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("inverse spectrum",&realstart,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
			/* Create fake data power spectral DENSITY */
			for(j=0;j<inputMCMC.invspec[i]->data->length;j++){ PSD(&status,&(inputMCMC.invspec[i]->data->data[j]),j*inputMCMC.deltaF);}
			/* Allocate buffer for fake freq domain data */
			inputMCMC.stilde[i] = (COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("stilde",&realstart,0.0,inputMCMC.deltaF,&lalDimensionlessUnit,seglen/2 +1);
			inputMCMC_N.stilde[i] = (COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("stilde",&realstart,0.0,inputMCMC_N.deltaF,&lalDimensionlessUnit,seglen/2 +1);
			memcpy(&(inputMCMC.stilde[i]->epoch),&segmentStart,sizeof(LIGOTimeGPS));
			/*			inputMCMC.stilde[i]->epoch = datastart;
			 XLALGPSAdd(&(inputMCMC.stilde[i]->epoch), (REAL8)TrigSegStart/(REAL8)SampleRate);*/
			
			/* Create the fake data */
			for(j=0;j<inputMCMC.invspec[i]->data->length;j++){
				inputMCMC.invspec[i]->data->data[j]=1.0/(scalefactor*inputMCMC.invspec[i]->data->data[j]);
				inputMCMC.stilde[i]->data->data[j].re=XLALNormalDeviate(datarandparam)/(2.0*sqrt(inputMCMC.invspec[i]->data->data[j]*inputMCMC.deltaF));
				inputMCMC.stilde[i]->data->data[j].im=XLALNormalDeviate(datarandparam)/(2.0*sqrt(inputMCMC.invspec[i]->data->data[j]*inputMCMC.deltaF));
			}
		}
		else FakeFlag=0;

		if(timeslides&&!FakeFlag){ /* Set up time slides by randomly offsetting the data */
			LALCreateRandomParams(&status,&randparam,seed);
			LALUniformDeviate(&status,&TSoffset,randparam);
			TSoffset=(TSoffset-0.5)*TIMESLIDE;
			datastart = realstart;
			XLALGPSAdd(&datastart, TSoffset);
			fprintf(stderr,"Slid %s by %f s from %10.10lf to %10.10lf\n",IFOnames[i],TSoffset,realstart.gpsSeconds+1e-9*realstart.gpsNanoSeconds,datastart.gpsSeconds+1e-9*datastart.gpsNanoSeconds);
			XLALDestroyRandomParams(randparam);
		}
		
		if(specifictimeslides && !FakeFlag){ /* Set up time slides by offsetting the data by user defined value */
			if( ( !strcmp(IFOnames[i],"H1") && H1GPSshift != 0.0 ) || ( !strcmp(IFOnames[i],"L1") &&
					L1GPSshift != 0.0 ) || ( !strcmp(IFOnames[i],"V1") && V1GPSshift != 0.0 ) ) {
				if(!strcmp(IFOnames[i],"H1"))
					TSoffset=H1GPSshift;
				else if(!strcmp(IFOnames[i],"L1"))
					TSoffset=L1GPSshift;
				else if(!strcmp(IFOnames[i],"V1"))
					TSoffset=V1GPSshift;
				datastart = realstart;
				XLALGPSAdd(&datastart, TSoffset);
				fprintf(stderr,"Slid %s by %f s from %10.10lf to %10.10lf\n",IFOnames[i],TSoffset,realstart.gpsSeconds+1e-9*realstart.gpsNanoSeconds,datastart.gpsSeconds+1e-9*datastart.gpsNanoSeconds);
			}
		}
		
		TrigSample=(INT4)(SampleRate*(ETgpsSeconds - realstart.gpsSeconds));
		TrigSample+=(INT4)(1e-9*SampleRate*ETgpsNanoseconds - 1e-9*SampleRate*realstart.gpsNanoSeconds);
		/*TrigSegStart=TrigSample+SampleRate*(0.5*(segDur-InjParams.tc)) - seglen; */ /* Centre the injection */
		TrigSegStart=TrigSample+ (2*SampleRate) - seglen; /* Put trigger 2 s before end of segment */
		if(InjParams.tc>segDur) fprintf(stderr,"Warning! Your template is longer than the data segment\n");
		XLALGPSAdd(&segmentStart, (REAL8)TrigSegStart/(REAL8)SampleRate);
		memcpy(&(inputMCMC.epoch),&segmentStart,sizeof(LIGOTimeGPS));
		
		/* set up a Tukey Window */
		if (inputMCMC.window==NULL) inputMCMC.window = windowplan = XLALCreateTukeyREAL8Window( seglen, 0.1*(8.0/(seglen/SampleRate))); /* 0.1 agreed on beta parameter for review with 8s window*/
		/* if (inputMCMC.window==NULL) inputMCMC.window = windowplan = XLALCreateTukeyREAL8Window( seglen,(REAL8)2.0*padding*SampleRate/(REAL8)seglen); */ /* Original window, commented out for review */
		/* Read the data from disk into a vector (RawData) */
		if(!FakeFlag){
			RawData = readTseries(CacheFileNames[i],ChannelNames[i],datastart,duration); /* This reads the raw data from the cache */
			if(RawData==NULL){fprintf(stderr,"Error opening %s in %s\n",ChannelNames[i],CacheFileNames[i]); exit(-1);}
			if(timeslides || specifictimeslides){
				memcpy(&(RawData->epoch),&realstart,sizeof(LIGOTimeGPS));
				memcpy(&datastart,&realstart,sizeof(LIGOTimeGPS));
			}
			/* Resample the time series */
			if(SampleRate) check=XLALResampleREAL8TimeSeries(RawData,1.0/SampleRate);
			if(check) {fprintf(stderr,"check=%d, failed to resample from %lf Hz to %d Hz\n",check,1.0/RawData->deltaT,SampleRate); exit(-1);}
			/* Need to resize the raw data to be an integer multiple of the seglen */
			fprintf(stderr,"Shrinking... (lost %d samples from end)\n",RawData->data->length-(seglen*nSegs));
			RawData=(REAL8TimeSeries *)XLALShrinkREAL8TimeSeries(RawData,(size_t) 0, (size_t) seglen*nSegs);
			/* Estimate the noise PSD */
			if(estimatenoise){ /* Spectrum not used with student-t likelihood */
				/* Set up inverse spectrum structure */
				inputMCMC.invspec[i] = (REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("inverse spectrum",&RawData->epoch,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);
				inputMCMC_N.invspec[i] = (REAL8FrequencySeries *)XLALCreateREAL8FrequencySeries("inverse spectrum",&RawData->epoch,0.0,(REAL8)(SampleRate)/seglen,&lalDimensionlessUnit,seglen/2 +1);

				/* Compute power spectrum */
				if(DEBUG) fprintf(stderr,"Computing power spectrum, seglen %i\n",seglen);
				check=XLALREAL8AverageSpectrumMedian( inputMCMC.invspec[i] ,RawData,(UINT4)seglen,(UINT4)stride,windowplan,fwdplan);
				check|=XLALREAL8SpectrumInvertTruncate( inputMCMC.invspec[i], inputMCMC.fLow, seglen, (seglen-stride)/4, fwdplan, revplan );

				if(check) {fprintf(stderr,"Cannot create spectrum, check=%x\n",check); exit(-1);}
				/* POWER SPECTRUM SHOULD HAVE UNITS OF TIME! */
			}

			/* Set up to read trigger time independently */
			/* Shift the data if requested like for PSD */
			LIGOTimeGPS realsegstart;
			memcpy(&realsegstart,&segmentStart,sizeof(LIGOTimeGPS));
               		if(timeslides&&!FakeFlag){ /* Set up time slides by randomly offsetting the data */
                        	LALCreateRandomParams(&status,&randparam,seed);
                        	LALUniformDeviate(&status,&TSoffset,randparam);
                        	TSoffset=(TSoffset-0.5)*TIMESLIDE;
                        	XLALGPSAdd(&segmentStart, TSoffset);
	                        fprintf(stderr,"Slid %s by %f s from %10.10lf to %10.10lf\n",IFOnames[i],TSoffset,realstart.gpsSeconds+1e-9*realstart.gpsNanoSeconds,datastart.gpsSeconds+1e-9*datastart.gpsNanoSeconds);
        	                XLALDestroyRandomParams(randparam);
			}
			if(specifictimeslides && !FakeFlag){ /* Set up time slides by offsetting the data by user defined value */
                        	if( ( !strcmp(IFOnames[i],"H1") && H1GPSshift != 0.0 ) || ( !strcmp(IFOnames[i],"L1") &&
                                        L1GPSshift != 0.0 ) || ( !strcmp(IFOnames[i],"V1") && V1GPSshift != 0.0 ) ) {
                                if(!strcmp(IFOnames[i],"H1"))
                                        TSoffset=-H1GPSshift;
                                else if(!strcmp(IFOnames[i],"L1"))
                                        TSoffset=-L1GPSshift;
                                else if(!strcmp(IFOnames[i],"V1"))
                                        TSoffset=-V1GPSshift;
                                XLALGPSAdd(&segmentStart, TSoffset);
                                fprintf(stderr,"Slid %s by %f s from %10.10lf to %10.10lf\n",IFOnames[i],TSoffset,realsegstart.gpsSeconds+1e-9*realsegstart.gpsNanoSeconds,segmentStart.gpsSeconds+1e-9*segmentStart.gpsNanoSeconds);
                        	}
                	}
			inputMCMC.segment[i]=readTseries(CacheFileNames[i],ChannelNames[i],segmentStart,(REAL8)seglen/SampleRate);
			/* Copy the real start time into vector */
			memcpy(&segmentStart,&realsegstart,sizeof(LIGOTimeGPS));
			memcpy(&(inputMCMC.segment[i]->epoch),&realsegstart,sizeof(LIGOTimeGPS));

			if(SampleRate) check=XLALResampleREAL8TimeSeries(inputMCMC.segment[i],1.0/SampleRate);

			if(InjParams.tc>segDur-padding) fprintf(stderr,"Warning, flat-top is shorter than injected waveform!\n");
			/* Store the appropriate data in the input structure */

			memcpy(&(inputMCMC.invspec[i]->epoch),&(inputMCMC.segment[i]->epoch),sizeof(LIGOTimeGPS));
			if(DEBUG) fprintf(stderr,"Data segment %d in %s from %f to %f, including padding\n",i,IFOnames[i],((float)TrigSegStart)/((float)SampleRate),((float)(TrigSegStart+seglen))/((float)SampleRate) );

			inputMCMC.stilde[i] = (COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("stilde",&realstart,0.0,inputMCMC.deltaF,&lalDimensionlessUnit,seglen/2 +1);
			inputMCMC_N.stilde[i] = (COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("stilde",&realstart,0.0,inputMCMC_N.deltaF,&lalDimensionlessUnit,seglen/2 +1);
			XLALDestroyREAL8TimeSeries(RawData);

			/* Window and FFT the data */
			XLALDDVectorMultiply(inputMCMC.segment[i]->data,inputMCMC.segment[i]->data,windowplan->data);
			check=XLALREAL8TimeFreqFFT(inputMCMC.stilde[i],inputMCMC.segment[i],fwdplan); /* XLALREAL8TimeFreqFFT multiplies by deltaT */
			for(j=0;j<inputMCMC.stilde[i]->data->length;j++) {
				inputMCMC.stilde[i]->data->data[j].re/=sqrt(windowplan->sumofsquares / windowplan->data->length);
				inputMCMC.stilde[i]->data->data[j].im/=sqrt(windowplan->sumofsquares / windowplan->data->length);
			}
		} /* End if(!FakeFlag) */

		/* Perform injection in time domain */
		if(NULL!=injXMLFile && fakeinj==0) {
			DetectorResponse det;
			REAL8 SNR=0.0;
			LIGOTimeGPS realSegStart;
			memset(&det,0,sizeof(DetectorResponse));
			det.site=inputMCMC.detector[i];
			/* Inject incoherently */
			if(decohereflag){
				memcpy(&realSegStart,&segmentStart,sizeof(realSegStart));
				XLALGPSAdd(&segmentStart,((REAL8) i+1)*offset);
				fprintf(stdout,"Offset injection by %lf s\n",((REAL8) i+1)*offset);
			}
			/* Create a buffer long enough to hold the signal */
			UINT4 bufferlength = (UINT4)(100.0/inputMCMC.deltaT);
			if(bufferlength<seglen) bufferlength=seglen;
			LIGOTimeGPS bufferstart;
			memcpy(&bufferstart,&segmentStart,sizeof(LIGOTimeGPS));
			XLALGPSAdd(&bufferstart,((REAL8)seglen*inputMCMC.deltaT));
			XLALGPSAdd(&bufferstart,-((REAL8)bufferlength*inputMCMC.deltaT));

			REAL4TimeSeries *injWave=(REAL4TimeSeries *)XLALCreateREAL4TimeSeries(IFOnames[i],&(bufferstart),0.0,inputMCMC.deltaT,&lalADCCountUnit,(size_t)bufferlength);

			for (j=0;j<injWave->data->length;j++) injWave->data->data[j]=0.0;
//			LALSimulateCoherentGW(&status,injWave,&InjectGW,&det);
			COMPLEX8FrequencySeries *resp = XLALCreateCOMPLEX8FrequencySeries("response",&bufferstart,0.0,inputMCMC.deltaF,(const LALUnit *)&strainPerCount,seglen);
			for(j=0;j<resp->data->length;j++) {resp->data->data[j].re=(REAL4)1.0; resp->data->data[j].im=0.0;}
			SimInspiralTable this_injection;
			memcpy(&this_injection,injTable,sizeof(SimInspiralTable));
			this_injection.next=NULL;
			LALFindChirpInjectSignals(&status,injWave,&this_injection,resp);
			XLALDestroyCOMPLEX8FrequencySeries(resp);
			printf("Finished InjectSignals\n");
			fprintf(stderr,"Cutting injection buffer from %d to %d\n",bufferlength,seglen);

			TrigSegStart=(INT4)((segmentStart.gpsSeconds-injWave->epoch.gpsSeconds)*SampleRate);
			TrigSegStart+=(INT4)((segmentStart.gpsNanoSeconds - injWave->epoch.gpsNanoSeconds)*1e-9*SampleRate);

			injWave=(REAL4TimeSeries *)XLALCutREAL4TimeSeries(injWave,TrigSegStart,seglen);
			fprintf(stderr,"Cut buffer start time=%lf, segment start time=%lf\n",injWave->epoch.gpsSeconds+1e-9*injWave->epoch.gpsNanoSeconds,inputMCMC.stilde[i]->epoch.gpsSeconds + 1.0e-9*inputMCMC.stilde[i]->epoch.gpsNanoSeconds);
			REPORTSTATUS(&status);
			if(decohereflag) {
				memcpy(&segmentStart,&realSegStart,sizeof(realSegStart));
				memcpy(&(injWave->epoch),&realSegStart,sizeof(realSegStart));
			}
			REAL8TimeSeries *inj8Wave=(REAL8TimeSeries *)XLALCreateREAL8TimeSeries("injection",&segmentStart,0.0,inputMCMC.deltaT,&lalDimensionlessUnit,(size_t)seglen);
			for (j=0;j<injWave->data->length;j++) inj8Wave->data->data[j]=(REAL8)injWave->data->data[j]; /* Move into a REAL8 vector */
			
			
			/* Compute the frequency domain wave for SNR calculation */
			RealFFTPlan *inj_plan = XLALCreateForwardREAL4FFTPlan( seglen, 0 );
			COMPLEX16FrequencySeries *injF = (COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("injFD",&(segmentStart),0.0,inputMCMC.deltaF,&lalDimensionlessUnit,seglen/2 +1);
			/* Window the data */
			REAL4 WinNorm = sqrt(windowplan->sumofsquares/windowplan->data->length);
			for(j=0;j<inj8Wave->data->length;j++) inj8Wave->data->data[j]*=SNRfac*windowplan->data->data[j]/WinNorm;
			XLALREAL8TimeFreqFFT(injF,inj8Wave,fwdplan); /* This calls XLALREAL8TimeFreqFFT which normalises by deltaT */
			
			REPORTSTATUS(&status);
			if(estimatenoise){
				for(j=(UINT4) (inputMCMC.fLow/inputMCMC.invspec[i]->deltaF),SNR=0.0;j<inputMCMC.invspec[i]->data->length;j++){
					SNR+=((REAL8)injF->data->data[j].re)*((REAL8)injF->data->data[j].re)*inputMCMC.invspec[i]->data->data[j];
					SNR+=((REAL8)injF->data->data[j].im)*((REAL8)injF->data->data[j].im)*inputMCMC.invspec[i]->data->data[j];}
				SNR*=4.0*inputMCMC.invspec[i]->deltaF; /* Get units correct - factor of 4 for 1-sided */
			}
			LALDestroyREAL4FFTPlan(&status,&inj_plan);

			networkSNR+=SNR;
            SNRs[i]=sqrt(SNR);
			SNR=sqrt(SNR);

			/* Actually inject the waveform */
			if(!FakeFlag) for(j=0;j<inj8Wave->data->length;j++) inputMCMC.segment[i]->data->data[j]+=(REAL8)inj8Wave->data->data[j];
			for(j=0;j<injF->data->length;j++) {
				inputMCMC_N.stilde[i]->data->data[j].re=inputMCMC.stilde[i]->data->data[j].re;
			  	inputMCMC_N.stilde[i]->data->data[j].im=inputMCMC.stilde[i]->data->data[j].im;
			  	inputMCMC_N.invspec[i]->data->data[j]=inputMCMC.invspec[i]->data->data[j];			
				inputMCMC.stilde[i]->data->data[j].re+=(REAL8)injF->data->data[j].re;
				inputMCMC.stilde[i]->data->data[j].im+=(REAL8)injF->data->data[j].im;
			}
#if DEBUG
			FILE *waveout;
			char wavename[100];
			sprintf(wavename,"wave_%s.dat",IFOnames[i]);
			waveout=fopen(wavename,"w");
			for(j=0;j<injF->data->length;j++) fprintf(waveout,"%10.10lf %10.10e %10.10e\n",j*inputMCMC.deltaF,injF->data->data[j].re,injF->data->data[j].im);
			fclose(waveout);
#endif
			XLALDestroyCOMPLEX16FrequencySeries(injF);

			XLALDestroyREAL4TimeSeries(injWave);
			XLALDestroyREAL8TimeSeries(inj8Wave);

			if(status.statusCode==0) {fprintf(stderr,"Injected signal into %s. SNR=%lf\n",IFOnames[i],SNR);}
			else {fprintf(stderr,"injection failed!!!\n"); REPORTSTATUS(&status); exit(-1);}
		}

	XLALDestroyRandomParams(datarandparam);
	} /* End loop over IFOs */
    
    if(NULL!=injXMLFile && fakeinj==0 && SNRpath!=NULL) {
    /* Print the SNRs in a file */
    PrintSNRsToFile(SNRs,injTable,&inputMCMC);
    }
    
	/* Data is now all in place in the inputMCMC structure for all IFOs and for one trigger */

	if(estimatenoise && DEBUG){
		for(j=0;j<nIFO;j++){
			char filename[100];
			sprintf(filename,"indata_%s.dat",IFOnames[j]);
			FILE *outinit=fopen(filename,"w");
			for(i=0;i<inputMCMC.stilde[j]->data->length;i++) fprintf(outinit,"%e %e %e %e %e %e\n",
						 inputMCMC.stilde[j]->f0 + i*inputMCMC.stilde[0]->deltaF,
						 inputMCMC.stilde[j]->data->data[i].re,
				 		 inputMCMC.stilde[j]->data->data[i].im,
						 1./inputMCMC.invspec[j]->data->data[i],
						 inputMCMC_N.stilde[j]->data->data[i].re,
				 		 inputMCMC_N.stilde[j]->data->data[i].im);	
			fclose(outinit);
		}
	}

	/* Set up the structure */
	inputMCMC.injectionTable = injTable;
	inputMCMC.numberDataStreams = nIFO;
	inputMCMC_N.numberDataStreams = nIFO;
	
	inputMCMC.numPoints = seglen;
	inputMCMC.stride = stride;
	inputMCMC.inspiralTable = inputCurrent;
	inputMCMC.fwdplan = fwdplan;
	inputMCMC.revplan = revplan;
	inputMCMC.numberDraw = Nmcmc;
	inputMCMC.annealingTemp = 0.1;
	/* randparams need to be handled differently from the MCMC code*/
	LALCreateRandomParams(&status,&(inputMCMC.randParams),seed);


	/* Set up the approximant to use in the likelihood function */
	CHAR TT2[]="TaylorT2"; CHAR TT3[]="TaylorT3"; CHAR TT4[]="TaylorT4"; CHAR TF2[]="TaylorF2";
	CHAR BBH[]="IMRPhenomFA"; CHAR BBHSpin1[]="IMRPhenomFB_NS"; CHAR BBHSpin2[]="IMRPhenomFB";
	CHAR BBHSpin3[]="IMRPhenomFB_Chi"; CHAR EBNR[]="EOBNR"; CHAR AMPCOR[]="AmpCorPPN";
	CHAR ST[]="SpinTaylor"; CHAR LowMassIMRFB[]="IMRPhenomFB_Chi_low"; CHAR LowMassIMRB[]="IMRPhenomB_Chi_low";
	CHAR PSTRD[]="PhenSpinTaylorRD";
	/*CHAR PSTRD[]="PhenSpinTaylorRD"; */ /* Commented out until PhenSpin waveforms are in master */
	inputMCMC.approximant = TaylorF2; /* Default */
	if(!strcmp(approx,TF2)) inputMCMC.approximant=TaylorF2;
	else if(!strcmp(approx,TT2)) inputMCMC.approximant=TaylorT2;
	else if(!strcmp(approx,TT3)) inputMCMC.approximant=TaylorT3;
    else if(!strcmp(approx,TT4)) inputMCMC.approximant=TaylorT4;
	else if(!strcmp(approx,BBH)) inputMCMC.approximant=IMRPhenomFA;
    else if(!strcmp(approx,BBHSpin1)) inputMCMC.approximant=IMRPhenomFB;
    else if(!strcmp(approx,BBHSpin2)) inputMCMC.approximant=IMRPhenomFB;
    else if(!strcmp(approx,BBHSpin3)) inputMCMC.approximant=IMRPhenomFB;
    else if(!strcmp(approx,LowMassIMRFB)) inputMCMC.approximant=IMRPhenomFB;
    else if(!strcmp(approx,LowMassIMRB)) inputMCMC.approximant=IMRPhenomB;
    else if(!strcmp(approx,EBNR)) inputMCMC.approximant=EOBNR;
	else if(!strcmp(approx,AMPCOR)) inputMCMC.approximant=AmpCorPPN;
	else if(!strcmp(approx,ST)) inputMCMC.approximant=SpinTaylor;
	else if(!strcmp(approx,PSTRD)) {
	    inputMCMC.approximant=PhenSpinTaylorRD;
	    if (strstr(approx,"inspiralOnly")) {
	      inputMCMC.inspiralOnly=1;
	    }
	    else {
	      inputMCMC.inspiralOnly=0;
	    }
	    if (strstr(approx,"fixedStep")) {
	      inputMCMC.fixedStep=1;
	    }
	    else {
	      inputMCMC.fixedStep=0;
	    }
	    if (strstr(approx,"TotalJ")) {
	      inputMCMC.axisChoice=LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J;
	    }
	    else if (strstr(approx,"OrbitalL")) {
	      inputMCMC.axisChoice=LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L;
	    }
	    else {
	      inputMCMC.axisChoice=LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW;
	    }
	}
	else {fprintf(stderr,"Unknown approximant: %s\n",approx); exit(-1);}

	if(inputMCMC.approximant!=AmpCorPPN && ampOrder!=0){
		fprintf(stderr,"Warning, setting amp order %i but not using AmpCorPPN. Amplitude corrected waveforms will NOT be generated!\n",ampOrder);
	}
	inputMCMC.ampOrder=ampOrder;

	if(phaseOrder>7 && inputMCMC.approximant!=EOBNR)
	{
		fprintf(stderr,"Error: Cannot go above 3.5PN in phase using this template!\n");
		exit(1);
	}
	switch(phaseOrder)
	{
		case 0:
		{
			inputMCMC.phaseOrder=LAL_PNORDER_NEWTONIAN;
			break;
		}
		case 1:
		{
			inputMCMC.phaseOrder=LAL_PNORDER_HALF;
			break;
		}
		case 2:
		{
			inputMCMC.phaseOrder=LAL_PNORDER_ONE;
			break;
		}
		case 3:
		{
			inputMCMC.phaseOrder=LAL_PNORDER_ONE_POINT_FIVE;
			break;
		}
		case 4:
		{
			inputMCMC.phaseOrder=LAL_PNORDER_TWO;
			break;
		}
		case 5:
		{
			inputMCMC.phaseOrder=LAL_PNORDER_TWO_POINT_FIVE;
			break;
		}
		case 6:
		{
			inputMCMC.phaseOrder=LAL_PNORDER_THREE;
			break;
		}
		case 7:
		{
			inputMCMC.phaseOrder=LAL_PNORDER_THREE_POINT_FIVE;
			break;
		}
		case 8:
		{
			inputMCMC.phaseOrder=LAL_PNORDER_PSEUDO_FOUR;
			break;
		}
		default:
			inputMCMC.phaseOrder=LAL_PNORDER_TWO;
	}

	/* Set the initialisation and likelihood functions */
	if(SkyPatch) {inputMCMC.funcInit = NestInitSkyPatch; goto doneinit;}
	if(SkyLocFlag) {inputMCMC.funcInit = NestInitSkyLoc; goto doneinit;}
	if(NULL!=inputXMLFile) inputMCMC.funcInit = NestInit2PN;
	else if(NINJA && NULL==injXMLFile) inputMCMC.funcInit = NestInitNINJAManual;
	else if(NINJA) inputMCMC.funcInit = NestInitInjNINJA;
	else {if(NULL!=injXMLFile) inputMCMC.funcInit = NestInitInj;
	else inputMCMC.funcInit = NestInitManual;}
doneinit:
	if(studentt) inputMCMC.funcLikelihood = MCMCSTLikelihoodMultiCoherentF;
	else inputMCMC.funcLikelihood = MCMCLikelihoodMultiCoherentF;
	if(inputMCMC.approximant==AmpCorPPN) inputMCMC.funcLikelihood = MCMCLikelihoodMultiCoherentAmpCor;
	if(inputMCMC.approximant==EOBNR) {
		inputMCMC.funcLikelihood = MCMCLikelihoodMultiCoherentAmpCor; /* EOBNR uses CoherentGW the same way */
		inputMCMC.phaseOrder = LAL_PNORDER_PSEUDO_FOUR;
	}
		
	inputMCMC.funcPrior = NestPrior;
	if(GRBflag) {inputMCMC.funcPrior = GRBPrior;
		inputMCMC.funcInit = NestInitGRB;
	}
	if(HighMassFlag) inputMCMC.funcPrior = NestPriorHighMass;
	if(skylocprior) inputMCMC.funcPrior = NestPriorSkyLoc;

    if(!strcmp(approx,BBHSpin1)) {
        inputMCMC.funcPrior = NestPriorHighMass;
        inputMCMC.funcLikelihood = MCMCLikelihoodMultiCoherentF;
        inputMCMC.funcInit = NestInitManual;
    }
		/* IMRPhenomFB */
    if(!strcmp(approx,BBHSpin2)) {
        inputMCMC.funcPrior = NestPriorHighMass;
        inputMCMC.funcLikelihood = MCMCLikelihoodMultiCoherentF;
        inputMCMC.funcInit = NestInitManualIMRB;
    }
		/* IMRPhenomFB_Chi */
    if(!strcmp(approx,BBHSpin3)) {
        inputMCMC.funcPrior = NestPriorHighMass;
        inputMCMC.funcLikelihood = MCMCLikelihoodMultiCoherentF;
        inputMCMC.funcInit = NestInitManualIMRBChi;
    }
    if((!strcmp(approx,LowMassIMRB)) || (!strcmp(approx,LowMassIMRFB)) ){
	inputMCMC.funcPrior = NestPrior;
	inputMCMC.funcLikelihood = MCMCLikelihoodMultiCoherentF;
	inputMCMC.funcInit = NestInitManualIMRBChi;
    }
   
    	if(!strcmp(approx,PSTRD)) {
	  inputMCMC.funcPrior = NestPrior;
	  inputMCMC.funcLikelihood = MCMCLikelihoodMultiCoherentF_PhenSpin;
	  inputMCMC.likelihoodPlan = NULL;
	  inputMCMC.funcInit = NestInitManualPhenSpinRD;
	  inputMCMC.Fwfc = XLALCreateREAL4Vector(inputMCMC.numPoints);
	  inputMCMC.Fwfp = XLALCreateREAL4Vector(inputMCMC.numPoints);
          InspiralTemplate template; 
		template.fCutoff=SampleRate/2.-1.;
		template.tSampling=SampleRate;
		template.approximant=PhenSpinTaylorRD;
		template.totalMass=mc2mass1(manual_mass_low,0.25)+mc2mass2(manual_mass_low,0.25);
		
		template.eta=0.24;
		template.massChoice=totalMassAndEta;
		template.fLower=inputMCMC.fLow;
        template.order=inputMCMC.phaseOrder;
        template.nStartPad = 0;
        template.nEndPad =0;
        template.startTime = 0.0;
        template.ieta = 1;

		LALInspiralParameterCalc(&status,&template);
		LALInspiralWaveLength(&status, &inputMCMC.mylength, template);
	        if(inputMCMC.mylength>inputMCMC.numPoints){	
		printf("myl = %d ,total mass=%11.4E\n",inputMCMC.mylength,template.totalMass);
		
		LALCreateForwardREAL4FFTPlan(&status,&inputMCMC.longplan,inputMCMC.mylength,FFTW_ESTIMATE);}
          #if DEBUGMODEL !=0 
          fprintf(stdout,"\n   Creating FFTW plan...\n");
          #endif
          LALCreateForwardREAL4FFTPlan(&status,&inputMCMC.likelihoodPlan,inputMCMC.numPoints,FFTW_ESTIMATE);
          #if DEBUGMODEL !=0
          fprintf(stdout,"    Done.\n");
          #endif
        

	} 
	
	if(HighMassFlag) inputMCMC.funcPrior = NestPriorHighMass;
     
	/* Live is an array of LALMCMCParameter * types */
	Live = (LALMCMCParameter **)LALMalloc(Nlive*sizeof(LALMCMCParameter *));
	for (i=0;i<Nlive;i++) Live[i]=(LALMCMCParameter *)LALMalloc(sizeof(LALMCMCParameter));

	if(networkSNR!=0.0) fprintf(stdout,"Injected signal network SNR= %lf\n",sqrt(networkSNR));

	double ReducedChiSq=0;
	/* variance of dimensionful real part d(f_k) (= variance of imaginary part) is zeta^2 */
	/* zeta^2 = N/(4deltaT) * S(f_k)  (S(f_k) dimensionful one-sided) */

	if(estimatenoise){
		for (i=(int)fLow/inputMCMC.invspec[0]->deltaF;i<inputMCMC.stilde[0]->data->length;i++) ReducedChiSq+=(pow(inputMCMC.stilde[0]->data->data[i].re,2.0)+pow(inputMCMC.stilde[0]->data->data[i].im,2.0))*inputMCMC.invspec[0]->data->data[i];
		ReducedChiSq *= 2.0*inputMCMC.invspec[0]->deltaF/(inputMCMC.stilde[0]->data->length-(fLow/inputMCMC.invspec[0]->deltaF)); /* should be N */
	}
	fprintf(stdout,"reduced chi squared = %e\n",ReducedChiSq);
	fprintf(stdout,"Number of points in F-domain above fLow = %i\n",(int)inputMCMC.stilde[0]->data->length-(int)(fLow/(double)inputMCMC.stilde[0]->deltaF));

	/* Output data if requested */
	if(datadump)
	{
		CHAR dumpfile[FILENAME_MAX];
		for(j=0;j<inputMCMC.numberDataStreams;j++){
			sprintf(dumpfile,"%s_%s.dat",datadump,IFOnames[j]);
			FILE *dataoutfile=fopen(dumpfile,"w");
			for(i=0;i<inputMCMC.stilde[j]->data->length;i++)
			{
				if(estimatenoise)
					fprintf(dataoutfile,"%12.5e  %14.8e  %14.8e  %14.8e  %14.8e  %14.8e\n",(REAL8)i*inputMCMC.invspec[j]->deltaF,1./inputMCMC.invspec[j]->data->data[i],inputMCMC.stilde[j]->data->data[i].re,inputMCMC.stilde[j]->data->data[i].im,inputMCMC_N.stilde[j]->data->data[i].re,inputMCMC_N.stilde[j]->data->data[i].im);
				else
					fprintf(dataoutfile,"%12.5e  %14.8e  %14.8e  %14.8e  %14.8e\n",(REAL8)i*inputMCMC.stilde[j]->deltaF,inputMCMC.stilde[j]->data->data[i].re,inputMCMC.stilde[j]->data->data[i].im,inputMCMC_N.stilde[j]->data->data[i].re,inputMCMC_N.stilde[j]->data->data[i].im);
			}
			fclose(dataoutfile);
		}
	}

			evidence=computeZ(&inputMCMC_N);
	#if DEBUG
	  fprintf(stdout,"Evidence of noise w/o injection: %18.8e\n",evidence);
	#endif
	
	evidence = nestZ(Nruns,Nlive,Live,&inputMCMC);
	fprintf(stdout,"logZ = %lf\n",evidence);

	/* Clean up */
	XLALDestroyREAL8Window(windowplan);
	if(!strcmp(approx,PSTRD)){
	XLALDestroyREAL4Vector(inputMCMC.Fwfc);
	XLALDestroyREAL4Vector(inputMCMC.Fwfp);
	}
	for(i=0;i<nIFO;i++){
		XLALDestroyCOMPLEX16FrequencySeries(inputMCMC.stilde[i]);
		if(estimatenoise) XLALDestroyREAL8FrequencySeries(inputMCMC.invspec[i]);
		XLALDestroyREAL8TimeSeries(inputMCMC.segment[i]);
	}
	return(0);
} /* End main() */



void NestInitManualPhenSpinRD(LALMCMCParameter *parameter, void *iT)
{
  (void)iT;

  double dmin=d_min;
  double dmax=d_max;
 

/*  if ( (manual_mass_high > manual_mass_low) && (manual_mass_low>2.*singleMassMin) ) {
    totalMassMin=manual_mass_low;
    totalMassMax=manual_mass_high;
  }
  else {
    if (HighMassFlag) {
      totalMassMin=m_tot_min_highmass;
      totalMassMax=m_tot_max_highmass;
    }
  }
  */

  double eta=etamin+gsl_rng_uniform(RNG)*(0.25-etamin);

  double logMc;
//  double mcmax=0.435275*totalMassMax;
 // double mcmin = m2mc(m1min,m2min);
 // double mcmax = m2mc(m1maxhalf,m2maxhalf);

  double lMcmin=log(manual_mass_low);
  double lMcmax=log(manual_mass_high);


  parameter->param=NULL;
  parameter->dimension = 0;
  
  logMc=log((gsl_rng_uniform(RNG)*(33.)+2.)*pow(eta,3./5.));

    XLALMCMCAddParam(parameter,"logmc",logMc,lMcmin,lMcmax,0);

    XLALMCMCAddParam(parameter,"eta",eta,etamin,0.25,0);

    XLALMCMCAddParam(parameter,"time",(gsl_rng_uniform(RNG)-0.5)*timewindow +manual_end_time,manual_end_time-0.5*timewindow,manual_end_time+0.5*timewindow,0);
    XLALMCMCAddParam(parameter,"phi", LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
    XLALMCMCAddParam(parameter,"distMpc", (dmax-dmin)*gsl_rng_uniform(RNG)+dmin,dmin,dmax,0);
   
    if ((long_min<=0.)&&(long_max>=2.*LAL_PI))
    XLALMCMCAddParam(parameter,"ra",gsl_rng_uniform(RNG)*2.*LAL_PI,0.,2.*LAL_PI,1);
    else
    XLALMCMCAddParam(parameter,"ra",gsl_rng_uniform(RNG)*(long_max-long_min)+long_min,long_min,long_max,0);

    XLALMCMCAddParam(parameter,"dec", lat_min+gsl_rng_uniform(RNG)*(lat_max-lat_min),lat_min,lat_max,0);
    XLALMCMCAddParam(parameter,"psi",gsl_rng_uniform(RNG)*LAL_PI,0,LAL_PI,1);
    XLALMCMCAddParam(parameter,"iota", gsl_rng_uniform(RNG)*(iota_max-iota_min)+iota_min ,iota_min,iota_max,0);
 
    if(onespin_flag==0 && nospin_flag==0){
    
    double spin1min=s1_mag_min;
    double spin1max=s1_mag_max;

    double spin2min=s2_mag_min;
    double spin2max=s2_mag_max;
    
    double spin1thetamin=s1_theta_min;
    double spin1thetamax=s1_theta_max;
        
    double spin2thetamin=s2_theta_min;
    double spin2thetamax=s2_theta_max;

    double spinphimin=s_phi_min;
    double spinphimax=s_phi_max;


    XLALMCMCAddParam(parameter,"a1",    (spin1max-spin1min)*gsl_rng_uniform(RNG)+spin1min,  spin1min, spin1max, 0);
    XLALMCMCAddParam(parameter,"theta1",(spin1thetamax-spin1thetamin)*gsl_rng_uniform(RNG)+spin1thetamin, spin1thetamin, spin1thetamax, 0);

    XLALMCMCAddParam(parameter,"a2",     (spin2max-spin2min)*gsl_rng_uniform(RNG)+spin2min,  spin2min, spin2max,0);
    XLALMCMCAddParam(parameter,"theta2",(spin2thetamax-spin2thetamin)*gsl_rng_uniform(RNG)+spin2thetamin,  spin2thetamin, spin2thetamax, 0);

    if ((spinphimax>=2.*LAL_PI)&&(spinphimin<=0.))
    XLALMCMCAddParam(parameter,"phi1",  2.*LAL_PI*gsl_rng_uniform(RNG), 0., 2.*LAL_PI,1);
   	else
    XLALMCMCAddParam(parameter,"phi1",  (spinphimax-spinphimin)*gsl_rng_uniform(RNG)+spinphimin,  spinphimin, spinphimax, 0);
			if ((spinphimax>=2.*LAL_PI)&&(spinphimin<=0.))
    XLALMCMCAddParam(parameter,"phi2",  2.*LAL_PI*gsl_rng_uniform(RNG), 0., 2.*LAL_PI,1);
        else
    XLALMCMCAddParam(parameter,"phi2",  (spinphimax-spinphimin)*gsl_rng_uniform(RNG)+spinphimin,  spinphimin, spinphimax, 0);		


			}
	
	else{
		if(onespin_flag==1 && nospin_flag==0){
			double spin2min=s2_mag_min;
			double spin2max=s2_mag_max;
        		double spin2thetamin=s2_theta_min;
			double spin2thetamax=s2_theta_max;
			double spinphimin=s_phi_min;
    			double spinphimax=s_phi_max;
			XLALMCMCAddParam(parameter,"a2",     (spin2max-spin2min)*gsl_rng_uniform(RNG)+spin2min,  spin2min, spin2max,0);
			XLALMCMCAddParam(parameter,"theta2",(spin2thetamax-spin2thetamin)*gsl_rng_uniform(RNG)+spin2thetamin,  spin2thetamin, spin2thetamax, 0);
		
			if ((spinphimax>=2.*LAL_PI)&&(spinphimin<=0.))
    XLALMCMCAddParam(parameter,"phi2",  2.*LAL_PI*gsl_rng_uniform(RNG), 0., 2.*LAL_PI,1);
        else
    XLALMCMCAddParam(parameter,"phi2",  (spinphimax-spinphimin)*gsl_rng_uniform(RNG)+spinphimin,  spinphimin, spinphimax, 0);
		
		}
	}
    return;
}

void NestInitGRB(LALMCMCParameter *parameter, void *iT){
	REAL8 grb_time;
	SimInspiralTable *injTable = (SimInspiralTable *)iT;
	REAL4 localetawin;
	REAL8 mcmin,mcmax,m1min,m1max,m2min,m2max;
	REAL8 deltaLong=0.01;
	REAL8 deltaLat=0.01;
	REAL8 trueLong=0.0,trueLat=0.0;

	parameter->param = NULL;
	parameter->dimension = 0;

	if(iT!=NULL){
		grb_time = (REAL8) injTable->geocent_end_time.gpsSeconds + (REAL8)injTable->geocent_end_time.gpsNanoSeconds *1.0e-9;
		trueLong = (REAL8)injTable->longitude;
		trueLat = (REAL8)injTable->latitude;
	}
	/*else*/   {
		grb_time = manual_end_time;
		if(manual_RA!=-4200.0) trueLong = manual_RA;
		if(manual_dec!=-4200.0) trueLat = manual_dec;
    }
	/*etamin = etamin<0.01?0.01:etamin;*/
	double etamax = 0.25;

	/* GRB priors are below */
	m1min=1.0;
	m1max=3.0;
	m2min=1.0;
	m2max=35.0;

	mcmin = m2mc(m1min,m2min);
	mcmax = m2mc(m1max,m2max);
	etamin = 0.027;

	localetawin=etamax-etamin;
	double lmmin=log(mcmin);
	double lmmax=log(mcmax);
	XLALMCMCAddParam(parameter,"logmc",lmmin+(lmmax-lmmin)*gsl_rng_uniform(RNG),lmmin,lmmax,0);

	/*  XLALMCMCAddParam(parameter,"mchirp",mcmin+(mcmax-mcmin)*gsl_rng_uniform(RNG),mcmin,mcmax,0);*/
	XLALMCMCAddParam(parameter, "eta", gsl_rng_uniform(RNG)*localetawin+etamin , etamin, etamax, 0);
	XLALMCMCAddParam(parameter, "time",             (gsl_rng_uniform(RNG)-0.5)*timewindow + grb_time ,grb_time-0.5*timewindow,grb_time+0.5*timewindow,0);
	XLALMCMCAddParam(parameter, "phi",              LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter, "distMpc", 99.0*gsl_rng_uniform(RNG)+1.0, 1.0, 100.0, 0);

	XLALMCMCAddParam(parameter,"ra",trueLong,trueLong-0.5*deltaLong,trueLong+0.5*deltaLong,-1);
	XLALMCMCAddParam(parameter,"dec",trueLat,trueLat-0.5*deltaLat,trueLat+0.5*deltaLat,-1);

	XLALMCMCAddParam(parameter,"psi",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,1);
	XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);


	return;
}

void NestInitSkyLoc(LALMCMCParameter *parameter, void *iT)
{
	SimInspiralTable *injTable = (SimInspiralTable *) iT;
	parameter->param=NULL;
	parameter->dimension=0;
	double inM1 = injTable->mass1;
	double inM2 = injTable->mass2;
	double inEta = injTable->eta;
	double inTime = injTable->geocent_end_time.gpsSeconds + 1e-9*injTable->geocent_end_time.gpsNanoSeconds;
	double inMc = m2mc(inM1,inM2);
	double deltaM=0.05; double deltaEta=0.01;
	double etaMin=inEta-0.5*deltaEta; double etaMax=inEta+0.5*deltaEta;
	etaMin=etaMin<0.0?0.0:etaMin;
	etaMax=etaMax>0.25?0.25:etaMax;
	deltaEta=etaMax-etaMin;
	double lmmin=log(inMc-deltaM);
	double lmmax=log(inMc+deltaM);
	double m_comp_min=1.0;
	double m_comp_max=15.0;
	
	if(m1m2Flag){
		XLALMCMCAddParam(parameter,"m1",m_comp_min+(m_comp_max-m_comp_min)*gsl_rng_uniform(RNG),m_comp_min,m_comp_max,0);
		XLALMCMCAddParam(parameter,"m2",m_comp_min+(m_comp_max-m_comp_min)*gsl_rng_uniform(RNG),m_comp_min,m_comp_max,0);
	}
	else{
		XLALMCMCAddParam(parameter,"logmc",lmmin+(lmmax-lmmin)*gsl_rng_uniform(RNG),lmmin,lmmax,0);

	/*  XLALMCMCAddParam(parameter,"mchirp",(gsl_rng_uniform(RNG)-0.5)*deltaM + inMc,inMc-0.5*deltaM,inMc+0.5*deltaM,0);*/
		XLALMCMCAddParam(parameter,"eta",(gsl_rng_uniform(RNG))*deltaEta + etaMin,etaMin,etaMax,0);
	}
	XLALMCMCAddParam(parameter,"time",(gsl_rng_uniform(RNG)-0.5)*timewindow+inTime,inTime-0.5*timewindow,inTime+0.5*timewindow,0);
	XLALMCMCAddParam(parameter,"phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter,"logdist", (log(40.0)-log(10.0))*gsl_rng_uniform(RNG)+log(10.0), log(10.0), log(40.0), 0);
	XLALMCMCAddParam(parameter,"ra",LAL_TWOPI*gsl_rng_uniform(RNG),0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter,"dec",LAL_PI*(gsl_rng_uniform(RNG)-0.5),-LAL_PI/2.0,LAL_PI/2.0,0);
	XLALMCMCAddParam(parameter,"psi",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,1);
	XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);
	return;
}

/* FIXME: parameter iT is unused */
void NestInitSkyPatch(LALMCMCParameter *parameter, void UNUSED *iT)
{
	double mcmin,mcmax;
	double deltaLong=0.001;
	double deltaLat=0.001;
	parameter->param=NULL;
	parameter->dimension = 0;
	fprintf(stderr,"Using longitude = %f, latitude = %f\n",manual_RA,manual_dec);
	mcmin=m2mc(manual_mass_low/2.0,manual_mass_low/2.0);
	mcmax=m2mc(manual_mass_high/2.0,manual_mass_high/2.0);

	double lmmin=log(mcmin);
	double lmmax=log(mcmax);
	XLALMCMCAddParam(parameter,"logmc",lmmin+(lmmax-lmmin)*gsl_rng_uniform(RNG),lmmin,lmmax,0);

	/*	XLALMCMCAddParam(parameter,"mchirp",mcmin+(mcmax-mcmin)*gsl_rng_uniform(RNG),mcmin,mcmax,0);*/
	/*	XLALMCMCAddParam(parameter,"mtotal",manual_mass_low+mwin*gsl_rng_uniform(RNG),manual_mass_low,manual_mass_high,0);*/
	XLALMCMCAddParam(parameter,"eta",etamin+gsl_rng_uniform(RNG)*(0.25-etamin),etamin,0.25,0);
	XLALMCMCAddParam(parameter,"time",(gsl_rng_uniform(RNG)-0.5)*timewindow +manual_end_time,manual_end_time-0.5*timewindow,manual_end_time+0.5*timewindow,0);
	XLALMCMCAddParam(parameter,"phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter,"distMpc", 99.0*gsl_rng_uniform(RNG)+1.0, 1.0, 100.0, 0);
	XLALMCMCAddParam(parameter,"ra",manual_RA,manual_RA-0.5*deltaLong,manual_RA+0.5*deltaLong,-1);
	XLALMCMCAddParam(parameter,"dec",manual_dec,manual_dec-0.5*deltaLat,manual_dec+0.5*deltaLat,-1);
	XLALMCMCAddParam(parameter,"psi",LAL_PI*gsl_rng_uniform(RNG), 0, LAL_PI, 1);
	XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);
	return;
}

/* FIXME: parameter iT is unused */
void NestInitManual(LALMCMCParameter *parameter, void UNUSED *iT)
{
	double mcmin,mcmax;
	parameter->param=NULL;
	parameter->dimension = 0;
	mcmin=manual_mass_low*pow(etamin,3./5.);
	mcmax=manual_mass_high*pow(0.25,3./5.);
	double lmmin=log(mcmin);
	double lmmax=log(mcmax);
	double lDmin=log(manual_dist_min);
	double lDmax=log(manual_dist_max);

	XLALMCMCAddParam(parameter,"logmc",lmmin+(lmmax-lmmin)*gsl_rng_uniform(RNG),lmmin,lmmax,0);
	/*	XLALMCMCAddParam(parameter,"mchirp",mcmin+(mcmax-mcmin)*gsl_rng_uniform(RNG),mcmin,mcmax,0);*/
	/*	XLALMCMCAddParam(parameter,"mtotal",manual_mass_low+mwin*gsl_rng_uniform(RNG),manual_mass_low,manual_mass_high,0);*/
	XLALMCMCAddParam(parameter,"eta",etamin+gsl_rng_uniform(RNG)*(0.25-etamin),etamin,0.25,0);
	XLALMCMCAddParam(parameter,"time",(gsl_rng_uniform(RNG)-0.5)*timewindow +manual_end_time,manual_end_time-0.5*timewindow,manual_end_time+0.5*timewindow,0);
	XLALMCMCAddParam(parameter,"phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
/*	XLALMCMCAddParam(parameter,"distMpc", (dmax-dmin)*gsl_rng_uniform(RNG)+dmin,dmin,dmax, 0);*/
	XLALMCMCAddParam(parameter,"logdist",lDmin+gsl_rng_uniform(RNG)*(lDmax-lDmin),lDmin,lDmax,0);
	XLALMCMCAddParam(parameter,"ra",LAL_TWOPI*gsl_rng_uniform(RNG),0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter,"dec",LAL_PI*(gsl_rng_uniform(RNG)-0.5),-LAL_PI/2.0,LAL_PI/2.0,0);
	XLALMCMCAddParam(parameter,"psi",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,1);
	XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);

	return;
}

/* FIXME: parameter iT is unused */
void NestInitManualIMRB(LALMCMCParameter *parameter, void UNUSED *iT)
{
	double mcmin,mcmax;
	parameter->param=NULL;
	parameter->dimension = 0;
	mcmin=manual_mass_low;
	mcmax=manual_mass_high;

    double lmmin=log(mcmin);
	double lmmax=log(mcmax);

    double ldmin=log(manual_dist_min);
    double ldmax=log(manual_dist_max);

    double spin1zmin=-1.;
    double spin1zmax=1.;

    double spin2zmin=-1.;
    double spin2zmax=1.;

	XLALMCMCAddParam(parameter,"logmc",lmmin+(lmmax-lmmin)*gsl_rng_uniform(RNG),lmmin,lmmax,0);
	/*	XLALMCMCAddParam(parameter,"mchirp",mcmin+(mcmax-mcmin)*gsl_rng_uniform(RNG),mcmin,mcmax,0);*/
	/*	XLALMCMCAddParam(parameter,"mtotal",manual_mass_low+mwin*gsl_rng_uniform(RNG),manual_mass_low,manual_mass_high,0);*/
	XLALMCMCAddParam(parameter,"eta",etamin+gsl_rng_uniform(RNG)*(0.25-etamin),etamin,0.25,0);
	XLALMCMCAddParam(parameter,"time",(gsl_rng_uniform(RNG)-0.5)*timewindow +manual_end_time,manual_end_time-0.5*timewindow,manual_end_time+0.5*timewindow,0);
	XLALMCMCAddParam(parameter,"phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
/*	XLALMCMCAddParam(parameter,"distMpc", (dmax-dmin)*gsl_rng_uniform(RNG)+dmin,dmin,dmax, 0);*/
	XLALMCMCAddParam(parameter,"logdist",ldmin+gsl_rng_uniform(RNG)*(ldmax-ldmin),ldmin,ldmax,0);
	XLALMCMCAddParam(parameter,"ra",LAL_TWOPI*gsl_rng_uniform(RNG),0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter,"dec",LAL_PI*(gsl_rng_uniform(RNG)-0.5),-LAL_PI/2.0,LAL_PI/2.0,0);
	XLALMCMCAddParam(parameter,"psi",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,1);
	XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);

    XLALMCMCAddParam(parameter,"spin1z",(spin1zmax-spin1zmin)*gsl_rng_uniform(RNG)+spin1zmin,spin1zmin,spin1zmax,0);
    XLALMCMCAddParam(parameter,"spin2z",(spin2zmax-spin2zmin)*gsl_rng_uniform(RNG)+spin2zmin,spin2zmin,spin2zmax,0);
	return;
}

/* FIXME: parameter iT is unused */
void NestInitManualIMRBChi(LALMCMCParameter *parameter, void UNUSED *iT)
{
	double mcmin,mcmax;
	parameter->param=NULL;
	parameter->dimension = 0;
	mcmin=manual_mass_low;
	mcmax=manual_mass_high;

	double lmmin=log(mcmin);
	double lmmax=log(mcmax);

    double ldmin=log(manual_dist_min);
    double ldmax=log(manual_dist_max);
    double chiSpinmin=manual_chi_min;
    double chiSpinmax=manual_chi_max;

	XLALMCMCAddParam(parameter,"logmc",lmmin+(lmmax-lmmin)*gsl_rng_uniform(RNG),lmmin,lmmax,0);
	/*	XLALMCMCAddParam(parameter,"mchirp",mcmin+(mcmax-mcmin)*gsl_rng_uniform(RNG),mcmin,mcmax,0);*/
	/*	XLALMCMCAddParam(parameter,"mtotal",manual_mass_low+mwin*gsl_rng_uniform(RNG),manual_mass_low,manual_mass_high,0);*/
	XLALMCMCAddParam(parameter,"eta",etamin+gsl_rng_uniform(RNG)*(0.25-etamin),etamin,0.25,0);
	XLALMCMCAddParam(parameter,"time",(gsl_rng_uniform(RNG)-0.5)*timewindow +manual_end_time,manual_end_time-0.5*timewindow,manual_end_time+0.5*timewindow,0);
	XLALMCMCAddParam(parameter,"phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
/*	XLALMCMCAddParam(parameter,"distMpc", (dmax-dmin)*gsl_rng_uniform(RNG)+dmin,dmin,dmax, 0);*/
	XLALMCMCAddParam(parameter,"logdist",ldmin+gsl_rng_uniform(RNG)*(ldmax-ldmin),ldmin,ldmax,0);
	XLALMCMCAddParam(parameter,"ra",LAL_TWOPI*gsl_rng_uniform(RNG),0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter,"dec",LAL_PI*(gsl_rng_uniform(RNG)-0.5),-LAL_PI/2.0,LAL_PI/2.0,0);
	XLALMCMCAddParam(parameter,"psi",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,1);
	XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);

    XLALMCMCAddParam(parameter,"chiSpin",(chiSpinmax-chiSpinmin)*gsl_rng_uniform(RNG)+chiSpinmin,chiSpinmin,chiSpinmax,0);

    return;
}

/* FIXME: parameter iT is unused */
void NestInitNINJAManual(LALMCMCParameter *parameter, void UNUSED *iT){
	REAL8 trg_time,mcmin,mcmax;
	REAL4 localetawin;
	parameter->param = NULL;
	parameter->dimension = 0;
	trg_time = manual_end_time;

	/*double etamin = eta-0.5*etawindow;
	 etamin = etamin<0.01?0.01:etamin;*/
	/*double etamax = eta+0.5*etawindow;
	 etamax = etamax>0.25?0.25:etamax;*/
	double etamax=0.25;
	localetawin=etamax-etamin;
	mcmin=m2mc(25.,25.);
	mcmax=m2mc(75.,75.);
	/*              parameter structure, name of parameter, initial value of parameter, minimum value parameter, maximum value of parameter, wrapped?) */
	XLALMCMCAddParam(parameter,"mchirp",mcmin+(mcmax-mcmin)*gsl_rng_uniform(RNG),mcmin,mcmax,0);
	/*XLALMCMCAddParam(parameter,"mtotal",gsl_rng_uniform(RNG)*100.0+50.0,50.0,150.0,0);*/
	/*XLALMCMCAddParam(parameter,"mtotal",3.0+27.0*gsl_rng_uniform(RNG),3.0,30.0,0);*/
	XLALMCMCAddParam(parameter, "eta", gsl_rng_uniform(RNG)*localetawin+etamin , etamin, etamax, 0);
	XLALMCMCAddParam(parameter, "time",             (gsl_rng_uniform(RNG)-0.5)*timewindow + trg_time ,trg_time-0.5*timewindow,trg_time+0.5*timewindow,0);
	XLALMCMCAddParam(parameter, "phi",              LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter, "distMpc", 499.0*gsl_rng_uniform(RNG)+1.0, 1.0, 500.0, 0);
	XLALMCMCAddParam(parameter,"ra",LAL_TWOPI*gsl_rng_uniform(RNG),0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter,"dec",LAL_PI*(gsl_rng_uniform(RNG)-0.5),-LAL_PI/2.0,LAL_PI/2.0,0);
	XLALMCMCAddParam(parameter,"psi",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,1);
	XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);


	return;
}

void NestInitInj(LALMCMCParameter *parameter, void *iT){
	REAL8 trg_time;
	SimInspiralTable *injTable = (SimInspiralTable *)iT;
	REAL4 UNUSED mtot, UNUSED eta, UNUSED mwindow, localetawin;
	REAL8 UNUSED mc, mcmin, mcmax, lmmin, lmmax;
	parameter->param = NULL;
	parameter->dimension = 0;
	trg_time = (REAL8) injTable->geocent_end_time.gpsSeconds + (REAL8)injTable->geocent_end_time.gpsNanoSeconds *1.0e-9;
	mtot = injTable->mass1 + injTable->mass2;
	eta = injTable->eta;
	mwindow = 0.2;
	/*etamin = etamin<0.01?0.01:etamin;*/
	etamin=0.01;
	double etamax = 0.25;
	mc=m2mc(injTable->mass1,injTable->mass2);
	mcmin=manual_mass_low*pow(etamin,3./5.);
	REAL8 m_comp_min=1., m_comp_max=15.;
	mcmax=manual_mass_high*pow(etamax,3./5.);

	lmmin=log(mcmin);
	lmmax=log(mcmax);
	localetawin=etamax-etamin;
	
	LALMCMCParam *head;
	
	if(m1m2Flag){
		if(checkParamInList(pinned_params,"m1")||checkParamInList(pinned_params,"m2"))
		{XLALMCMCAddParam(parameter,"m1",injTable->mass1,m_comp_min,m_comp_max,0);
		XLALMCMCAddParam(parameter,"m2",injTable->mass2,m_comp_min,m_comp_max,0);
		}
		else {
			XLALMCMCAddParam(parameter,"m1",m_comp_min+(m_comp_max-m_comp_min)*gsl_rng_uniform(RNG),m_comp_min,m_comp_max,0);
			XLALMCMCAddParam(parameter,"m2",m_comp_min+(m_comp_max-m_comp_min)*gsl_rng_uniform(RNG),m_comp_min,m_comp_max,0);
		}
	}
	else{
		if(checkParamInList(pinned_params,"logmc")||checkParamInList(pinned_params,"mchirp"))
			XLALMCMCAddParam(parameter,"logmc",log(injTable->mchirp),lmmin,lmmax,-1);
		else
			XLALMCMCAddParam(parameter,"logmc",lmmin+(lmmax-lmmin)*gsl_rng_uniform(RNG),lmmin,lmmax,0);
		if(checkParamInList(pinned_params,"eta"))
			XLALMCMCAddParam(parameter,"eta",injTable->eta,etamin,etamax,-1);
		else
			XLALMCMCAddParam(parameter, "eta", gsl_rng_uniform(RNG)*localetawin+etamin , etamin, etamax, 0);
		
	}
	if(checkParamInList(pinned_params,"time"))
		XLALMCMCAddParam(parameter,"time",trg_time,trg_time-0.5*timewindow,trg_time+0.5*timewindow,-1);
	else
		XLALMCMCAddParam(parameter, "time",		(gsl_rng_uniform(RNG)-0.5)*timewindow + trg_time,trg_time-0.5*timewindow,trg_time+0.5*timewindow,0);
	
	if(checkParamInList(pinned_params,"phi"))
		XLALMCMCAddParam(parameter,"phi",injTable->coa_phase,0,LAL_TWOPI,-1);
	else
		XLALMCMCAddParam(parameter, "phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
	
	if(checkParamInList(pinned_params,"dist") || checkParamInList(pinned_params,"logdist") || checkParamInList(pinned_params,"distance") || checkParamInList(pinned_params,"logdistance"))
		XLALMCMCAddParam(parameter,"logdist",log(injTable->distance),log(manual_dist_min),log(manual_dist_max),-1);
	else
		XLALMCMCAddParam(parameter,"logdist",(log(manual_dist_max)-log(manual_dist_min))*gsl_rng_uniform(RNG)+log(manual_dist_min) ,log(manual_dist_min),log(manual_dist_max),0);

	if(checkParamInList(pinned_params,"ra")||checkParamInList(pinned_params,"longitude")||checkParamInList(pinned_params,"RA"))
		XLALMCMCAddParam(parameter,"ra",injTable->longitude,0,LAL_TWOPI,-1);
	else
		XLALMCMCAddParam(parameter,"ra",gsl_rng_uniform(RNG)*LAL_TWOPI,0,LAL_TWOPI,1);
	if(checkParamInList(pinned_params,"dec") || checkParamInList(pinned_params,"latitude") || checkParamInList(pinned_params,"dec"))
		XLALMCMCAddParam(parameter,"dec",injTable->latitude,-LAL_PI/2.0,LAL_PI/2.0,-1);
	else
		XLALMCMCAddParam(parameter,"dec", acos(2.0*gsl_rng_uniform(RNG)-1.0)-LAL_PI/2.0,-LAL_PI/2.0,LAL_PI/2.0,0);

	if(checkParamInList(pinned_params,"psi")||checkParamInList(pinned_params,"polarization"))
		XLALMCMCAddParam(parameter,"psi",injTable->polarization,0,LAL_PI,-1);
	else
		XLALMCMCAddParam(parameter,"psi",gsl_rng_uniform(RNG)*LAL_PI,0,LAL_PI,1);
	
	if(checkParamInList(pinned_params,"iota") || checkParamInList(pinned_params,"inclination"))
		XLALMCMCAddParam(parameter,"iota", injTable->inclination, 0, LAL_PI, -1);
	else
		XLALMCMCAddParam(parameter,"iota", acos(2.0*gsl_rng_uniform(RNG)-1.0) ,0,LAL_PI,0);

	for (head=parameter->param;head;head=head->next)
	{
		if(head->core->wrapping==-1)
			fprintf(stdout,"Fixed parameter %s to %lf\n",head->core->name,head->value);
	}

	return;

}

int checkParamInList(const char *list, const char *param)
{
	/* Check for param in comma-seperated list */
	char *post=NULL,*pos=NULL;
	if (list==NULL) return 0;
	if (param==NULL) return 0;

	if(!(pos=strstr(list,param))) return 0;
	
	/* The string is a substring. Check that it is a token */
	/* Check the character before and after */
	if(pos!=list)
		if(*(pos-1)!=',')
			return 0;

	post=&(pos[strlen(param)]);
	if(*post!='\0')
		if(*post!=',')
			return 0;
	return 1;
}

void PrintSNRsToFile(REAL8* SNRs,SimInspiralTable *inj_table,LALMCMCInput *inputMCMC){
/* open the SNR file */
    char SnrName[70];
    char ListOfIFOs[10];
    REAL8 NetSNR=0.0;
    snprintf(ListOfIFOs, 10, " ");

    for (UINT4 det_i=0;det_i<nIFO;det_i++){
         sprintf(ListOfIFOs,"%s%s",ListOfIFOs,inputMCMC->ifoID[det_i]);
        }
    
    sprintf(SnrName,"%s/snr_%s_%10.1f.dat",SNRpath,ListOfIFOs,(REAL8) inj_table->geocent_end_time.gpsSeconds+ (REAL8) inj_table->geocent_end_time.gpsNanoSeconds*1.0e-9);
    FILE * snrout = fopen(SnrName,"w");
    if(!snrout){
	fprintf(stderr,"Unable to open the path %s for writing SNR files\n",SNRpath);
	exit(1);
    }
    for (UINT4 det_i=0;det_i<nIFO;det_i++){
        fprintf(snrout,"%s:\t",inputMCMC->ifoID[det_i]);
        fprintf(snrout,"%4.2f\n",SNRs[det_i]);
        NetSNR+=(SNRs[det_i]*SNRs[det_i]);
    }		
    if (nIFO>1){  fprintf(snrout,"Network:\t");
    fprintf(snrout,"%4.2f\n",sqrt(NetSNR));}
    fprintf(snrout,"\n");
    fclose(snrout);

}
