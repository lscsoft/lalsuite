/*********************************************************************************/
/*         Program to generate an input file for lalapps_makefakedata            */
/*              that consists of randomly generated parameters                   */
/*                                                                               */
/*			           C. Messenger                                  */
/*                                                                               */
/*                         BIRMINGHAM UNIVERISTY -  2004                         */
/*********************************************************************************/

#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <time.h>
#include <getopt.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include "GenerateBinaryMesh_v1.h"
#include <lal/Random.h>

INT4 lalDebugLevel=3;
static LALStatus status;

REAL8 alphaMIN,alphaMAX,sindeltaMIN,sindeltaMAX;
REAL8 smaMIN,smaMAX,periodMIN,periodMAX,eccMIN,eccMAX,argpMIN,argpMAX;
REAL8 tobsMIN,tobsMAX;
REAL8 tstartMIN,tstartMAX;
LIGOTimeGPS tperiMIN,tperiMAX;
CHAR infile[256],outfile[256],stamps[256],ifo[256];
REAL8 phi,psi,cosiota;
REAL8 f0MIN,f0MAX,h0MIN,h0MAX;
REAL8 f1dot,f2dot,f3dot;
REAL8 aplus,across;
INT4 seed;
INT4 nsft,band,fmin;
REAL8 tobs;
INT4 tstart,tsft;
BOOLEAN alphaflag=0;
BOOLEAN deltaflag=0;
BOOLEAN f0flag=0;
BOOLEAN h0flag=0;
BOOLEAN smaflag=0;
BOOLEAN periodflag=0;
BOOLEAN tperiflag=0;
BOOLEAN eccflag=0;
BOOLEAN argpflag=0;
BOOLEAN phiflag=0;
BOOLEAN psiflag=0;
BOOLEAN cosiotaflag=0;
BOOLEAN aplusflag=0;
BOOLEAN acrossflag=0;
BOOLEAN ifoflag=0;
BOOLEAN tobsflag=0;
BOOLEAN tstartflag=0;
BOOLEAN nsftflag=0;
BOOLEAN tsftflag=0;
BOOLEAN bandflag=0;
BOOLEAN fminflag=0;
BOOLEAN allskyflag=0;
REAL8 alpha,delta,sindelta,f0,h0,sma,period,argp,ecc;
LIGOTimeGPS tperiGPS;
INT4 ndet=6;  /* number of detectors to choose from */

extern char *optarg;
extern int optind, opterr, optopt;
int GenRandomParams(); 
int GenTimeStamps();
int OuptutRandomConfigFile();
int MakeSafeBinaryParams(); 
int ReadCommandLine(int argc,char *argv[]);
int OutputRandomConfigFile();
ssize_t getline(char **lineptr, size_t *n, FILE *stream);

int main(int argc,char *argv[]) 
{
  
  if (ReadCommandLine(argc,argv)) return 1;

  if (GenRandomParams()) return 2;
  
  if (tstartflag) {
    if (GenTimeStamps()) return 3;
    }

  if (OutputRandomConfigFile()) return 3;
  

  return 0;

}

/*******************************************************************************/

int GenTimeStamps() 
{

  FILE *fpstamps;
  INT4 j;

  nsft=floor(tobs/tsft);
  tobs=(REAL8)nsft*tsft;

  /* opening the output timestamps file */
  fpstamps=fopen(stamps,"w");
  if (fpstamps==NULL) {
    fprintf(stderr,"Unable to open file %s\n",stamps);
    return 1;
  }
  
  /* write out timestamps file */
  for (j=0;j<nsft;j++) {
    fprintf(fpstamps,"%d\t0\n",tstart+(j*tsft));
  }

  fclose(fpstamps);

  return 0;

}

/*******************************************************************************/

int GenRandomParams() 
{
  
  RandomParams *params=NULL;
  REAL4Vector *vector=NULL;
  LALTimeInterval deltaT;
  REAL8 deltaTFLT;
  UINT4 i,p,k;
  REAL8 phiMIN=0.0;
  REAL8 phiMAX=LAL_TWOPI;
  REAL8 psiMIN=0.0;
  REAL8 psiMAX=LAL_TWOPI;
  REAL8 cosiotaMIN=-1.0;
  REAL8 cosiotaMAX=1.0;
  REAL8 RandDeltaT,rand; 
 

  p=0;
  /* work out how many random variables we will require */
  if (allskyflag) p=p+2;
  if ((f0MAX-f0MIN)!=0.0) p++;
  if ((h0MAX-h0MIN)!=0.0) p++;
  if ((smaMAX-smaMIN)!=0.0) p++;
  if ((periodMAX-periodMIN)!=0.0) p++;
  if ((tperiMAX.gpsSeconds-tperiMIN.gpsSeconds)!=0) p++;
  if ((eccMAX-eccMIN)!=0.0) p++;
  if ((argpMAX-argpMIN)!=0.0) p++;
  if (phiflag) p++;
  if (psiflag) p++;
  if (cosiotaflag) p++;
  if ((tobsMAX-tobsMIN)!=0.0) p++;
  if ((tstartMAX-tstartMIN)!=0) p++;
  if (ifoflag) p++;

  /* create a vector of length p and generate random parameters*/
  if (p>0) {
    LALCreateVector(&status,&vector,p);
    LALCreateRandomParams(&status,&params,seed);
  }

  /* fill vector with random sequence */
  for (i=0;i<vector->length;i++)
    {
      LALUniformDeviate(&status,vector->data+i,params);
    }

  k=0;
  /* calculate the random numbers */ 
  if (alphaflag) {
    alpha=alphaMIN+vector->data[k]*(alphaMAX-alphaMIN);
    k++;
  }
  if (deltaflag) {
    sindelta=sindeltaMIN+vector->data[k]*(sindeltaMAX-sindeltaMIN);
    delta=asin(sindelta);
    k++;
  }
  if (f0flag) {
    f0=f0MIN+vector->data[k]*(f0MAX-f0MIN);
    k++;
  }
  if (h0flag) {
    h0=h0MIN+vector->data[k]*(h0MAX-h0MIN);
    k++;
  }
  if (smaflag) {
    {
      REAL8 a,b,smaPOW;
      a=log10(smaMIN);
      b=log10(smaMAX);
      smaPOW=a+vector->data[k]*(b-a);
      sma=pow(10.0,smaPOW);
      k++;
    }
  }
  if (periodflag) {
    {
      REAL8 a,b,periodPOW;
      a=log10(periodMIN);
      b=log10(periodMAX);
      periodPOW=a+vector->data[k]*(b-a);
      period=pow(10.0,periodPOW);
      k++;
    }
  }
  if (tperiflag) {
    /* find difference between min and max */
    LALDeltaGPS(&status,&deltaT,&tperiMAX,&tperiMIN);
    LALIntervalToFloat(&status,&deltaTFLT,&deltaT);
   
    /* calculate random increment */
    RandDeltaT=vector->data[k]*deltaTFLT;
   
    /* add random increment */
    LALAddFloatToGPS(&status,&tperiGPS,&tperiMIN,RandDeltaT);
   
    k++;
  }
  if (eccflag) {
    /* do logarithmic interval */
    if (eccMIN>0.0) {
      REAL8 a,b,eccPOW;
      a=log10(eccMIN);
      b=log10(eccMAX);
      eccPOW=a+vector->data[k]*(b-a);
      ecc=pow(10.0,eccPOW);
      k++;
    }
    else {
      /* do uniform interval */
      ecc=vector->data[k]*eccMAX;
    }
  }
  if (argpflag) {
    argp=argpMIN+vector->data[k]*(argpMAX-argpMIN);
    k++;
  }
  if (phiflag) {
    phi=phiMIN+vector->data[k]*(phiMAX-phiMIN);
    k++;
  }
  if (psiflag) {
    psi=psiMIN+(REAL8)vector->data[k]*(psiMAX-psiMIN);
    k++;
  }
  if (cosiotaflag) {
    cosiota=cosiotaMIN+vector->data[k]*(cosiotaMAX-cosiotaMIN);
    aplus=0.5*h0*(1.0+cosiota*cosiota);
    across=h0*cosiota;
    k++;
  }
  if (tobsflag) {
    {
      REAL8 a,b,tobsPOW;
      a=log10(tobsMIN);
      b=log10(tobsMAX);
      tobsPOW=a+vector->data[k]*(b-a);
      tobs=pow(10.0,tobsPOW);
      k++;
    }
  }
  if (tstartflag) {
    tstart=(INT4)((REAL8)tstartMIN+vector->data[k]*(REAL8)(tstartMAX-tstartMIN));
    k++;
  }
  if (ifoflag) {
    rand=(REAL8)ndet*vector->data[k];
    if ((rand>=0.0)&&(rand<1.0)) sprintf(ifo,"GEO");
    if ((rand>=1.0)&&(rand<2.0)) sprintf(ifo,"LLO");
    if ((rand>=2.0)&&(rand<3.0)) sprintf(ifo,"LHO"); 
    if ((rand>=3.0)&&(rand<4.0)) sprintf(ifo,"VIRGO"); 
    if ((rand>=4.0)&&(rand<5.0)) sprintf(ifo,"TAMA"); 
    if ((rand>=5.0)&&(rand<6.0)) sprintf(ifo,"CIT"); 
  }
  
  /* clean up memory */
  if (p>0) {
    LALDestroyRandomParams(&status,&params);
    LALDestroyVector(&status,&vector);
  }

  /* calculate number of sfts */
  if ((tsftflag) || (tobsflag)) {
    nsft=floor(tobs/tsft);
    tobs=(REAL8)nsft*tsft;
    nsftflag=1;
  }

  /* calculate min frequency */
  if ((bandflag) || (fminflag)) {
    fmin=floor(f0-((REAL8)band/2.0));
  }

  /* make the binary parameters realistic for the chosen sft length */
  if (MakeSafeBinaryParams()) return 4;

  /* make sure that the periapse passage is before the start time */
  if ((tperiGPS.gpsSeconds>tstart)&&(tstartflag)) 
    {
      REAL8 norb=0.0;
      norb=(REAL8)(tperiGPS.gpsSeconds-tstart+600)/period;
      tperiGPS.gpsSeconds=tperiGPS.gpsSeconds-(INT4)((floor(norb)+1.0)*period);
    }

 

  return 0;
 
}

/*******************************************************************************/

int MakeSafeBinaryParams() 
{

  /* This section just checks wether we are going to get significant */
  /* frequency drift within an SFT and wether we are getting into a  */
  /* relativistic regime. */

  REAL8 maxTsft,projvmax,ratio;
  REAL8 safe=0.25;
  REAL8 relativesafe=0.1;

  /* here we test for frequency drift within an sft */
  maxTsft=safe*sqrt((1.0-ecc)/(f0*sma))*period*(1.0-ecc)/(LAL_TWOPI*(1+ecc));
    /* if things are extreme then reduce the period to make things safe */
  if (maxTsft<(REAL8)tsft) {
    ratio=floor(tsft/maxTsft)+1.0;
      period=period*(2.0*ratio);
    }

  /* here we test for relativistic properties */
  projvmax=LAL_TWOPI*sma*sqrt((1.0+ecc)/(1.0-ecc))/period;
    if (projvmax>relativesafe) {
    period=period*(floor(projvmax/relativesafe)+1.0);
  }


  return 0;

}


/*******************************************************************************/

  int ReadCommandLine(int argc,char *argv[]) 
{
  INT4 c, errflg = 0;
  char *temp;
  optarg = NULL;
  
  /* Initialize default values */
  alphaMIN=0.0; 
  alphaMAX=LAL_TWOPI; 
  sindeltaMIN=-1.0; 
  sindeltaMAX=1.0; 
  allskyflag=0; /* A */

  phiflag=0; /* p */
  psiflag=0; /* P */
  cosiotaflag=0; /* c */
  
  h0MIN=1.0; /* g */
  h0MAX=1.0; /* G */
  h0=1.0;
  f0MIN=600.0; /* f */
  f0MAX=600.0; /* F */
  band=0;
 
  smaMIN=0.0; /* r */
  smaMAX=0.0; /* R */
  eccMIN=0.0; /* e */
  eccMAX=0.0; /* E */
  tperiMIN.gpsSeconds=0; /* t */
  tperiMIN.gpsNanoSeconds=0;
  tperiMAX.gpsSeconds=0; /* T */
  tperiMAX.gpsNanoSeconds=0;
  periodMIN=0.0; /* q */
  periodMAX=0.0; /* Q */
  argpMIN=0.0; /* k */
  argpMAX=0.0; /* K */

  tsft=60; /* v */
  tobsMIN=1e4; /* b */
  tobsMAX=1e4; /* B */
  tobs=1e4;
  tstartMIN=731100000; /* w */
  tstartMAX=731100000; /* W */
  tstart=731100000;
  ifoflag=0; /* I */

  sprintf(infile,"bin.cfg"); /* i */
  sprintf(outfile,"rand.cfg"); /* o */
  /* sprintf(stamps,"stamps.data"); */ /* m */
  seed=0; /* z */

  {
    int option_index = 0;
    static struct option long_options[] = {
                {"allsky", no_argument, 0, 'A'},
		{"phi", no_argument, 0, 'p'},
                {"psi", no_argument, 0, 'P'},
                {"cosiota", no_argument, 0, 'c'},
                {"h0MIN", required_argument, 0, 'g'},
		{"h0MAX", required_argument, 0, 'G'},
                {"f0MIN", required_argument, 0, 'f'},
		{"f0MAX", required_argument, 0, 'F'},
		{"smaMIN", required_argument, 0, 'r'},
		{"smaMAX", required_argument, 0, 'R'},
                {"periodMIN", required_argument, 0, 'q'},
		{"periodMAX", required_argument, 0, 'Q'},
                {"tperiMIN", required_argument, 0, 't'},
                {"tperiMAX", required_argument, 0, 'T'},
                {"eccMIN", required_argument, 0, 'e'},
		{"eccMAX", required_argument, 0, 'E'},
                {"argpMIN", required_argument, 0, 'k'},
		{"argpMAX", required_argument, 0, 'K'},
		{"tsft", required_argument, 0, 'v'},
		{"tobsMIN", required_argument, 0, 'b'},
		{"tobsMAX", required_argument, 0, 'B'},
		{"tstartMIN", required_argument, 0, 'w'},
		{"tstartMAX", required_argument, 0, 'W'},
		{"det", no_argument, 0, 'I'},
		{"band", required_argument, 0, 'n'},
		{"stamps", required_argument, 0, 'm'},
		{"infile", required_argument, 0, 'i'},
		{"outfile", required_argument, 0, 'o'},
		{"seed", required_argument, 0, 'z'},
		{"help", no_argument, 0, 'h'}
    };
  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt_long (argc, argv,"hA:pPcg:G:f:F:r:R:q:Q:t:T:e:E:k:K:b:B:w:W:Iz:n:i:m:o:",long_options, &option_index)))!=-1)
    switch (c) {
      
    case 'A':
      alphaflag=1;
      deltaflag=1;
      allskyflag=1;
      break;
    case 'p':
      phiflag=1;
      break;
    case 'P':
      psiflag=1;
      break;
    case 'c':
      cosiotaflag=1;
      aplusflag=1;
      acrossflag=1;
      break;
    case 'g':
      h0MIN=atof(optarg);
      h0flag=1;
      aplusflag=1;
      acrossflag=1;
      break;
    case 'G':
      h0MAX=atof(optarg);
      h0flag=1;
      aplusflag=1;
      acrossflag=1;
      break;
    case 'f':
      f0MIN=atof(optarg);
      f0flag=1;
      break;
    case 'F':
      f0MAX=atof(optarg);
      f0flag=1;
      break;
    case 'r':
      smaMIN=atof(optarg);
      smaflag=1;
      break;
    case 'R':
      smaMAX=atof(optarg);
      smaflag=1;
      break;
    case 'q':
      periodMIN=atof(optarg);
      periodflag=1;
      break;
    case 'Q':
      periodMAX=atof(optarg);
      periodflag=1;
      break;
    case 't':
      tperiMIN.gpsSeconds=atoi(optarg);
      tperiflag=1;
      break;
    case 'T':
      tperiMAX.gpsSeconds=atoi(optarg);
      tperiflag=1;
      break;
    case 'e':
      eccMIN=atof(optarg);
      eccflag=1;
      break;
    case 'E':
      eccMAX=atof(optarg);
      eccflag=1;
      break;
    case 'k':
      argpMIN=atof(optarg);
      argpflag=1;
      break;
    case 'K':
      argpMAX=atof(optarg);
      argpflag=1;
      break;
    case 'b':
      tobsMIN=atof(optarg);
      tobsflag=1;
      break;
    case 'B':
      tobsMAX=atof(optarg);
      tobsflag=1;
      break;
    case 'w':
      tstartMIN=atoi(optarg);
      tstartflag=1;
      break;
    case 'W':
      tstartMAX=atoi(optarg);
      tstartflag=1;
      break;
    case 'I':
      ifoflag=1;
      break;
    case 'z':
      seed=atoi(optarg);
      break;
    case 'v':
      tsft=atoi(optarg);
      tsftflag=1;
      break;
    case 'n':
      band=atoi(optarg);
      bandflag=1;
      fminflag=1;
      break;
    case 'm':
      temp=optarg;
      sprintf(stamps,temp);
      break; 
    case 'i':
      temp=optarg;
      sprintf(infile,temp);
      break;
    case 'o':
      temp=optarg;
      sprintf(outfile,temp);
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\t--allsky    BOOLEAN\t Set if random sky position required [DEFAULT=false]\n");
      fprintf(stdout,"\t--phi       BOOLEAN\t Set if random initial phase required  [DEFAULT=false]\n");
      fprintf(stdout,"\t--psi       BOOLEAN\t Set if random polarisation angle required  [DEFAULT=false]\n");
      fprintf(stdout,"\t--cosiota   BOOLEAN\t Set if random Cos(iota) required [DEFAULT=false]\n");
      fprintf(stdout,"\t--h0MIN     FLOAT\t MIN Gravitational wave amplitude [DEFAULT=1.0]\n");
      fprintf(stdout,"\t--h0MAX     FLOAT\t MAX Gravitational wave amplitude [DEFAULT=1.0]\n");
      fprintf(stdout,"\t--f0MIN     FLOAT\t MIN Gravitational wave frequency in Hz [DEFAULT=600.0]\n");
      fprintf(stdout,"\t--f0MAX     FLOAT\t MAX Gravitational wave frequency in Hz [DEFAULT=600.0]\n");    
      fprintf(stdout,"\t--smaMIN    FLOAT\t Projected semi-major axis of orbit in seconds [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--smaMAX    FLOAT\t Projected semi-major axis of orbit in seconds [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--periodMIN FLOAT\t MIN Period of orbit in seconds [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--periodMAX FLOAT\t MAX Period of orbit in seconds [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--tperiMIN  INTEGER\t MIN Observed time of periapse passage in seconds [DEFAULT=0] \n");
      fprintf(stdout,"\t--tperiMAX  INTEGER\t MAX Observed time of periapse passage in nanoseconds [DEFAULT=0]\n");
      fprintf(stdout,"\t--eccMIN    FLOAT\t MIN Orbital eccentricity [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--eccMAX    FLOAT\t MAX Orbital eccentricity [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--argpMIN   FLOAT\t MIN Argument of orbital periapse in radians [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--argpMAX   FLOAT\t MAX Argument of orbital periapse in radians [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--seed      INTEGER\t Seed for random number generation [DEFAULT=0 (from clock)]\n");
      fprintf(stdout,"\t--tsft      INTEGER\t The sft length required in seconds [DEFAULT=60]\n");
      fprintf(stdout,"\t--tobsMIN   FLOAT\t MIN value of random observation time [DEFAULT=1e4]\n");
      fprintf(stdout,"\t--tobsMAX   FLOAT\t MAX value of random observation time [DEFAULT=1e4]\n");
      fprintf(stdout,"\t--det       BOOLEAN\t Set if random detector required  [DEFAULT=false]\n");
      fprintf(stdout,"\t--tstartMIN INTEGER\t MIN value of GPS start of observation time in seconds [DEFAULT=731100000]\n");
      fprintf(stdout,"\t--tstartMAX INTEGER\t MAX value of GPS start of observation time in seconds [DEFAULT=731100000]\n");
      fprintf(stdout,"\t--stamps    STRING\t Name of output timstamps file [DEFAULT=NULL]\n"); 
      fprintf(stdout,"\t--infile    STRING\t Name of input configuration file [DEFAULT=NULL]\n");
      fprintf(stdout,"\t--outfile   STRING\t Name of output random configuration file [DEFAULT=NULL]\n");
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

int OutputRandomConfigFile() 
{
  FILE *fpin,*fpout,*fpcla;
  char name[256];
  char *line;
  char temp[256];
  char clargs[2048];
  /* defining the actual strings to identify in the input file */
  char *phi0text="phi0";
  char *psitext="psi";
  char *longitudetext="longitude";
  char *latitudetext="latitude";
  char *ifotext="detector";
  char *nsfttext="nTsft";
  char *tsfttext="Tsft";
  char *tstarttext="startTime";
  char *aPlustext="aPlus";
  char *aCrosstext="aCross";
  char *f0text="f0";
  char *orbitSemiMajorAxistext="orbitSemiMajorAxis";
  char *orbitPeriodtext="orbitPeriod";
  char *orbitEccentricitytext="orbitEccentricity";
  char *orbitArgPeriapsetext="orbitArgPeriapse";
  char *orbitTperiSSBsectext="orbitTperiSSBsec";
  char *orbitTperiSSBnstext="orbitTperiSSBns";
  char *fmintext="fmin";
  char *bandtext="Band";
  size_t len = 0;
  int i;
 
  /* opening the input config file */
  fpin=fopen(infile,"r");
  if (fpin==NULL) {
    fprintf(stderr,"Unable to open file %s\n",infile);
    return 1;
  }
  /* opening the output tconfig file */
  fpout=fopen(outfile,"w");
  if (fpout==NULL) {
    fprintf(stderr,"Unable to open file %s\n",outfile);
    return 1;
  }

  /* He we are simply outputting to file a string of commandline arguments  */
  /* which can be catted into the program lalapps_GenerateSearchInput so    */
  /* that a search of exactly matched parameters can be done */
  fpcla=fopen("clargsfile.out","w");
  if (fpcla==NULL) {
    fprintf(stderr,"Unable to open file\n");
    return 1;
  }
  
  
  i=0;
  /* cycle through the input file and identify the appropriate lines */ 
  while (getline(&line, &len, fpin)!=EOF) 
  {
    /* get the first column string to compare */
    sscanf(line,"%s",name);
    
    /* if any of the following strings are found and are to be replaced the change line */
    if ((strcmp(name,ifotext)==0)&&(ifoflag)) {
      sprintf(line,"detector\t= %s\t\t# Detector: LHO, LLO, VIRGO, GEO, TAMA, CIT, ROME\n",ifo);
      sprintf(temp,"--ifo %s ",ifo); 
      strcat(clargs,temp);
    }
    if ((strcmp(name,nsfttext)==0)&&(nsftflag)) {
      sprintf(line,"nTsft\t\t= %d\t\t# number of SFTs to calculate\n",nsft);
    }
    if ((strcmp(name,fmintext)==0)&&(fminflag)) {
      sprintf(line,"fmin\t\t= %d\t\t# lowest SFT-frequency in Hz\n",fmin);
    }
    if ((strcmp(name,bandtext)==0)&&(bandflag)) {
      sprintf(line,"Band\t\t= %d\t\t# SFT frequency band in Hz\n",band);
    }
    if ((strcmp(name,longitudetext)==0)&&(alphaflag)) {
      sprintf(line,"longitude\t= %6.12f\t# source longitude (in_radians)\n",alpha);
      sprintf(temp,"--alpha %6.12f ",alpha); 
      strcat(clargs,temp);
    }
    if ((strcmp(name,latitudetext)==0)&&(deltaflag)) {
      sprintf(line,"latitude\t= %6.12f\t# source latitude (in radians)\n",delta);
      sprintf(temp,"--delta %6.12f ",delta); 
      strcat(clargs,temp);
    }
    if ((strcmp(name,aPlustext)==0)&&(aplusflag)) {
      sprintf(line,"aPlus\t\t= %6.12f\t# plus-polarization amplitude a_+ (strain)\n",aplus);
    }
    if ((strcmp(name,aCrosstext)==0)&&(acrossflag)) {
      sprintf(line,"aCross\t\t= %6.12f\t# cross-polarization amplitude a_x (strain)\n",across);
    }
    if ((strcmp(name,psitext)==0)&&(psiflag)) {
      sprintf(line,"psi\t\t= %6.12f\t# wave polarization angle Psi\n",psi);
    }
    if ((strcmp(name,phi0text)==0)&&(phiflag)) {
      sprintf(line,"phi0\t\t= %6.12f\t# initial wave-phase phi0 (at reference-time tRef)\n",phi);
    }
    if ((strcmp(name,f0text)==0)&&(f0flag)) {
      sprintf(line,"f0\t\t= %6.12f\t# intrinsic signal frequency f0 (at tRef)\n",f0);
      f0flag=0;
      sprintf(temp,"--fmin %6.12f ",f0); 
      strcat(clargs,temp);
    }
    if ((strcmp(name,tsfttext)==0)&&(tsftflag)) {
      sprintf(line,"Tsft\t\t= %d\t# length of SFTs in seconds\n",tsft);
    }
    if ((strcmp(name,tstarttext)==0)&&(tstartflag)) {
      sprintf(line,"startTime\t= %d\t# GPS start time of (contiguous) output time-series\n",tstart);
    }
    if ((strcmp(name,orbitSemiMajorAxistext)==0)&&(smaflag)) {
      sprintf(line,"orbitSemiMajorAxis\t= %6.12f\t# Projected orbital semi-major axis a in seconds (i.e. a*sin(i)/c)\n",sma);
      sprintf(temp,"--smaxis %6.12f ",sma); 
      strcat(clargs,temp);
    }
    if ((strcmp(name,orbitEccentricitytext)==0)&&(eccflag)) {
      sprintf(line,"orbitEccentricity\t= %6.12f\t# Orbital eccentricity\n",ecc);
       sprintf(temp,"--ecc %6.12f ",ecc); 
      strcat(clargs,temp);
    }
    if ((strcmp(name,orbitTperiSSBsectext)==0)&&(tperiflag)) {
      sprintf(line,"orbitTperiSSBsec\t= %d\t\t# 'observed' (SSB) time of periapsis passage. Seconds.\n",tperiGPS.gpsSeconds);
       sprintf(temp,"--tperisec %d ",tperiGPS.gpsSeconds); 
      strcat(clargs,temp);
    }
    if ((strcmp(name,orbitTperiSSBnstext)==0)&&(tperiflag)) {
      sprintf(line,"orbitTperiSSBns\t\t= %d\t\t# 'observed' (SSB) time of periapsis passage. Nanoseconds.\n",tperiGPS.gpsNanoSeconds);
       sprintf(temp,"--tperinan %d ",tperiGPS.gpsNanoSeconds); 
      strcat(clargs,temp);
    }
    if ((strcmp(name,orbitPeriodtext)==0)&&(periodflag)) {
      sprintf(line,"orbitPeriod\t\t= %6.12f\t# Orbital period (seconds)\n",period);
      sprintf(temp,"--period %6.12f ",period); 
      strcat(clargs,temp);
    }
    if ((strcmp(name,orbitArgPeriapsetext)==0)&&(argpflag)) {
      sprintf(line,"orbitArgPeriapse\t= %6.12f\t# Argument of periapsis (radians)\n",argp); 
      sprintf(temp,"--argp %6.12f ",argp); 
      strcat(clargs,temp);
    }
    
    /* output the appropriate line to the output file */
    fprintf(fpout,"%s",line);

    
  }
  fprintf(fpcla,"%s",clargs);
 
  fclose(fpin);
  fclose(fpout);
  fclose(fpcla);
  return 0;
  
}
/*******************************************************************************/

