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
#include <lal/Random.h>
#include "GenerateBinaryMesh_v1.h"
#include "GenerateRandomBinaryInput_v1.h"


INT4 lalDebugLevel=3;
static LALStatus status;

REAL8 raMIN,raMAX,decMIN,decMAX;
REAL8 smaMIN,smaMAX,periodMIN,periodMAX,eccMIN,eccMAX,argpMIN,argpMAX;
REAL8 tobsMIN,tobsMAX;
REAL8 tstartMIN,tstartMAX;
LIGOTimeGPS tperiMIN,tperiMAX;
CHAR infile[256],primarystampsfile[256],secondarystampsfile[256];
CHAR primarytemplatefile[256],secondarytemplatefile[256];
CHAR primaryoutputfile[256],secondaryoutputfile[256];
CHAR primarynoisedir[256],secondarynoisedir[256];
CHAR primarysftbase[256],secondarysftbase[256];
REAL8 freqMIN,freqMAX,h0MIN,h0MAX;
INT4 primarystart,secondarystart;
INT4 seed;
REAL8 freqband,f_min;
INT4 tstart,tsft;

/* these flags are set from the clargs */
BOOLEAN detflag=0;
BOOLEAN freqflag=0;
BOOLEAN h0flag=0;
BOOLEAN smaflag=0;
BOOLEAN periodflag=0;
BOOLEAN tperiflag=0;
BOOLEAN eccflag=0;
BOOLEAN argpflag=0;
BOOLEAN phiflag=0;
BOOLEAN psiflag=0;
BOOLEAN cosiotaflag=0;
BOOLEAN ifoflag=0;
BOOLEAN tobsflag=0;
BOOLEAN tstartflag=0;
BOOLEAN bintempfileflag=0;
BOOLEAN coflag=0;
BOOLEAN stampsflag=0;
BOOLEAN allskyflag=0;
INT4 ndet=2;  /* number of detectors to choose from */

extern char *optarg;
extern int optind, opterr, optopt;

int GenRandomParams(RandomParameters *); 
int OutputRandomConfigFile(char *,RandomParameters *,RandomParameters *);
int ReadCommandLine(int argc,char *argv[]);
int ReadOrbitalParams(char *,RandomParameters *);
int AddRandomParams(RandomParameters *);
int InitialiseRandParams(RandomParameters *);

int main(int argc,char *argv[]) 
{
  
  RandomParameters primaryrandparams;
  RandomParameters secondaryrandparams;
  RandomParameters *randparams=NULL;

  /* read the clargs */
  if (ReadCommandLine(argc,argv)) return 1;

  /*printf("primary temp file is %s\n",primarytemplatefile);*/

  /* initilise all the random boundary structures */
  if (InitialiseRandParams(&primaryrandparams)) return 2;
  if (InitialiseRandParams(&secondaryrandparams)) return 2;

  /* read in the random parameter boundaries from the mesh files */
  if (bintempfileflag) {
    if (ReadOrbitalParams(primarytemplatefile,&primaryrandparams)) return 2;
    /* primaryrandparams.start=primarystart; */
    LALSnprintf(primaryrandparams.noisedir,256,primarynoisedir);
    LALSnprintf(primaryrandparams.sftbase,256,primarysftbase);
    LALSnprintf(primaryrandparams.stampsfile,256,primarystampsfile);
    randparams=&primaryrandparams;

    /*printf("done primary template read\n");*/

    if (coflag) {
      if (ReadOrbitalParams(secondarytemplatefile,&secondaryrandparams)) return 2;
      /* secondaryrandparams.start=secondarystart; */
      LALSnprintf(secondaryrandparams.noisedir,256,secondarynoisedir);
      LALSnprintf(secondaryrandparams.sftbase,256,secondarysftbase);
      LALSnprintf(secondaryrandparams.stampsfile,256,secondarystampsfile);
      /*printf("done secondary template read\n");*/

      if (primaryrandparams.start<secondaryrandparams.start) {
	randparams=&primaryrandparams;
      }
      else {
	randparams=&secondaryrandparams;
      }
    }
    
  }
  else {
    randparams=&primaryrandparams;
    LALSnprintf(randparams->noisedir,256,primarynoisedir);
    LALSnprintf(randparams->sftbase,256,primarysftbase);
    LALSnprintf(randparams->stampsfile,256,primarystampsfile);
  }

  /*printf("sorted the template boundaries\n");*/

  /* add all the other parameter ranges from clargs */
  if (AddRandomParams(randparams)) return 3;

  /*printf("added the rand cla boundaries\n");*/

  /* generate the random parameters */
  if (GenRandomParams(randparams)) return 2;
  
  /*printf("generated the rand params\n");*/

  /* output the random parameter config files */ 
  if (OutputRandomConfigFile(primaryoutputfile,randparams,&primaryrandparams)) return 3;
  
  /*printf("output to primary\n");*/

  if (coflag) {
    if (OutputRandomConfigFile(secondaryoutputfile,randparams,&secondaryrandparams)) return 3;
  }

  /*printf("output to secondary\n");*/

  return 0;

}

/*******************************************************************************/

int AddRandomParams(RandomParameters *randparams) 
{

  /* This function completes the randparams boundaries using the clargs */

  if (smaflag) {
    if (!bintempfileflag) {
      randparams->sma=(REAL8range *)LALMalloc(sizeof(REAL8range));
    }
    randparams->sma->min=smaMIN;
    randparams->sma->max=smaMAX;
  }
  if (tperiflag) {
    if (!bintempfileflag) randparams->tperi=(LIGOTimeGPSrange *)LALMalloc(sizeof(LIGOTimeGPSrange));
    randparams->tperi->min.gpsSeconds=tperiMIN.gpsSeconds;
    randparams->tperi->max.gpsSeconds=tperiMAX.gpsSeconds;
    randparams->tperi->min.gpsNanoSeconds=tperiMIN.gpsNanoSeconds;
    randparams->tperi->max.gpsNanoSeconds=tperiMAX.gpsNanoSeconds;
  }
  if (periodflag) {
    if (!bintempfileflag) randparams->period=(REAL8range *)LALMalloc(sizeof(REAL8range));
    randparams->period->min=periodMIN;
    randparams->period->max=periodMAX;
  }
  if (eccflag) {
    if (!bintempfileflag) randparams->ecc=(REAL8range *)LALMalloc(sizeof(REAL8range));
    randparams->ecc->min=eccMIN;
    randparams->ecc->max=eccMAX;
  }
  if (argpflag) {
    if (!bintempfileflag) randparams->argp=(REAL8range *)LALMalloc(sizeof(REAL8range));
    randparams->argp->min=argpMIN;
    randparams->argp->max=argpMAX;
  }
  if (freqflag) {
    randparams->freq=(REAL8range *)LALMalloc(sizeof(REAL8range));
    randparams->freq->min=freqMIN;
    randparams->freq->max=freqMAX;
  }
  if (phiflag) {
    randparams->phi=(REAL8range *)LALMalloc(sizeof(REAL8range));
    randparams->phi->min=0.0;
    randparams->phi->max=LAL_TWOPI;
  }
  if (psiflag) {
    randparams->psi=(REAL8range *)LALMalloc(sizeof(REAL8range));
    randparams->psi->min=0.0;
    randparams->psi->max=LAL_TWOPI;
  }
  if (cosiotaflag) {
    randparams->cosiota=(REAL8range *)LALMalloc(sizeof(REAL8range));
    randparams->cosiota->min=-1.0;
    randparams->cosiota->max=1.0;
  }
  if (h0flag) {
    randparams->h0=(REAL8range *)LALMalloc(sizeof(REAL8range));
    randparams->h0->min=h0MIN;
    randparams->h0->max=h0MAX;
  }
  if (allskyflag) {
    randparams->ra=(REAL8range *)LALMalloc(sizeof(REAL8range));
    randparams->ra->min=0.0;
    randparams->ra->max=LAL_TWOPI;
    randparams->dec=(REAL8range *)LALMalloc(sizeof(REAL8range));
    randparams->dec->min=-LAL_PI/2.0;
    randparams->dec->max=LAL_PI/2.0;
  }
  if (detflag) {
    randparams->det=(INT4range *)LALMalloc(sizeof(INT4range));
    randparams->det->min=0;
    randparams->det->max=ndet;
  }

  
  return 0;

}

/*******************************************************************************/

int ReadOrbitalParams(char *templatefile,RandomParameters *randomparams) 
{

  /* this function reads in the orbital parameter space boundaries from the */
  /* input template files.  It also gets the sky location */


  FILE *fp;
  char smaminname[256];
  char smamaxname[256];
  char raname[256];
  char decname[256];
  char tperisecminname[256];
  char tperisecmaxname[256];
  char tperinanominname[256];
  char tperinanomaxname[256];
  char eccminname[256];
  char eccmaxname[256];
  char argpminname[256];
  char argpmaxname[256];
  char periodminname[256];
  char periodmaxname[256];
  BinaryMeshFileHeader BMFheader;
  
  sprintf(smaminname,"Source_projected_orbital_semi_major_axis_MIN_sec");
  sprintf(smamaxname,"Source_projected_orbital_semi_major_axis_MAX_sec");
  sprintf(raname,"Search_right_ascension");
  sprintf(decname,"Search_declination");
  sprintf(tperisecminname,"Source_SSB_periapse_passage_MIN_GPS_sec");
  sprintf(tperisecmaxname,"Source_SSB_periapse_passage_MAX_GPS_sec");
  sprintf(tperinanominname,"Source_SSB_periapse_passage_MIN_GPS_nanosec");
  sprintf(tperinanomaxname,"Source_SSB_periapse_passage_MAX_GPS_nanosec");
  sprintf(eccminname,"Source_orbital_eccentricity_MIN");
  sprintf(eccmaxname,"Source_orbital_eccentricity_MAX");
  sprintf(argpminname,"Source_argument_of_periapse_MIN_rad");
  sprintf(argpmaxname,"Source_argument_of_periapse_MAX_rad");
  sprintf(periodminname,"Source_orbital_period_MIN_sec");
  sprintf(periodmaxname,"Source_orbital_period_MAX_sec");
  

  /* open the primary binary template file */
  fp = fopen(templatefile,"r");
  if (fp==NULL) {
    fprintf(stderr,"ERROR : could not open template file %s\n",templatefile);
    exit(1);
  }

  if (ReadMeshFileHeader(fp,&BMFheader)) return 2;
  
  /* extract the parameters from the header */
  LALSnprintf(randomparams->detector,256,BMFheader.det);

  randomparams->sma=(REAL8range *)LALMalloc(sizeof(REAL8range));
  randomparams->sma->min=BMFheader.sma_MIN;
  randomparams->sma->max=BMFheader.sma_MAX;

  randomparams->tperi=(LIGOTimeGPSrange *)LALMalloc(sizeof(LIGOTimeGPSrange));
  randomparams->tperi->min.gpsSeconds=BMFheader.tperi_MIN.gpsSeconds;
  randomparams->tperi->max.gpsSeconds=BMFheader.tperi_MAX.gpsSeconds;
  randomparams->tperi->min.gpsNanoSeconds=BMFheader.tperi_MIN.gpsNanoSeconds;
  randomparams->tperi->max.gpsNanoSeconds=BMFheader.tperi_MAX.gpsNanoSeconds;
  
  randomparams->ecc=(REAL8range *)LALMalloc(sizeof(REAL8range));
  randomparams->ecc->min=BMFheader.ecc_MIN;
  randomparams->ecc->max=BMFheader.ecc_MAX;

  randomparams->argp=(REAL8range *)LALMalloc(sizeof(REAL8range));
  randomparams->argp->min=BMFheader.argp_MIN;
  randomparams->argp->max=BMFheader.argp_MAX;

  randomparams->period=(REAL8range *)LALMalloc(sizeof(REAL8range));
  randomparams->period->min=BMFheader.period_MIN;
  randomparams->period->max=BMFheader.period_MAX;

  randomparams->ra=(REAL8range *)LALMalloc(sizeof(REAL8range));
  randomparams->ra->min=BMFheader.RA;
  randomparams->ra->max=BMFheader.RA;

  randomparams->dec=(REAL8range *)LALMalloc(sizeof(REAL8range));
  randomparams->dec->min=BMFheader.dec;
  randomparams->dec->max=BMFheader.dec;

  randomparams->start=BMFheader.tstart.gpsSeconds;

  fclose(fp);

  return 0;
  
}

/*******************************************************************************/

int InitialiseRandParams(RandomParameters *random) 
{

  /* this function just initialises the random parameters by pointing them to NULL */

  /* allocate some memory */
  /* rand=(RandomParameters *)LALMalloc(sizeof(RandomParameters)); */

  random->ra=NULL;
  random->dec=NULL;
  random->sma=NULL;
  random->tperi=NULL;
  random->period=NULL;
  random->ecc=NULL;
  random->argp=NULL;
  random->freq=NULL;
  random->phi=NULL;
  random->psi=NULL;
  random->cosiota=NULL;
  random->h0=NULL;
  random->ra=NULL;
  random->dec=NULL;
  random->det=NULL;

  return 0;

}

/*******************************************************************************/

int GenRandomParams(RandomParameters *randparams) 
{
  
  RandomParams *params=NULL;
  REAL4Vector *vector=NULL;
  LALTimeInterval deltaT;
  REAL8 deltaTFLT;
  UINT4 i,p,k;
  REAL8 RandDeltaT; 
  REAL8 sindecMIN,sindecMAX,sindec;

  p=0;
  /* work out how many random variables we will require */
  if (randparams->sma!=NULL) p++;
  if (randparams->tperi!=NULL) p++;
  if (randparams->period!=NULL) p++;
  if (randparams->ecc!=NULL) p++;
  if (randparams->argp!=NULL) {
    p++;
  }
  if (randparams->freq!=NULL) p++;
  if (randparams->h0!=NULL) p++;
  if (randparams->phi!=NULL) p++;
  if (randparams->psi!=NULL) p++;
  if (randparams->cosiota!=NULL) p++;
  if (randparams->ra!=NULL) p++;
  if (randparams->dec!=NULL) p++;
  if (randparams->det!=NULL) p++;

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
  if (randparams->ra!=NULL) {
    randparams->ra->value=randparams->ra->min+vector->data[k]*(randparams->ra->max-randparams->ra->min);
    k++;
  }
  if (randparams->dec!=NULL) {
    sindecMIN=sin(randparams->dec->min);
    sindecMAX=sin(randparams->dec->max);
    sindec=sindecMIN+vector->data[k]*(sindecMAX-sindecMIN);
    randparams->dec->value=asin(sindec);
    k++;
  }
  if (randparams->freq!=NULL) {
    randparams->freq->value=randparams->freq->min+vector->data[k]*(randparams->freq->max-randparams->freq->min);
    k++;
  }
  if (randparams->h0!=NULL) {
    randparams->h0->value=randparams->h0->min+vector->data[k]*(randparams->h0->max-randparams->h0->min);
    k++;
  }
  if (randparams->sma!=NULL) {
    randparams->sma->value=randparams->sma->min+vector->data[k]*(randparams->sma->max-randparams->sma->min);
    k++;
  }
  if (randparams->period!=NULL) {
    randparams->period->value=randparams->period->min+vector->data[k]*(randparams->period->max-randparams->period->min);
    k++;
  }
  if (randparams->tperi!=NULL) {
    /* find difference between min and max */
    LALDeltaGPS(&status,&deltaT,&randparams->tperi->max,&randparams->tperi->min);
    LALIntervalToFloat(&status,&deltaTFLT,&deltaT);
   
    /* calculate random increment */
    RandDeltaT=vector->data[k]*deltaTFLT;
   
    /* add random increment */
    LALAddFloatToGPS(&status,&randparams->tperi->value,&randparams->tperi->min,RandDeltaT);
    k++;
  }
  if (randparams->ecc!=NULL) {
    randparams->ecc->value=randparams->ecc->min+(REAL8)vector->data[k]*(randparams->ecc->max-randparams->ecc->min);
    k++;
  }
  if (randparams->argp!=NULL) {
    randparams->argp->value=randparams->argp->min+(REAL8)vector->data[k]*(randparams->argp->max-randparams->argp->min);
    k++;
  }
  if (randparams->phi!=NULL) {
    randparams->phi->value=randparams->phi->min+(REAL8)vector->data[k]*(randparams->phi->max-randparams->phi->min);
    k++;
  }
  if (randparams->psi!=NULL) {
    randparams->psi->value=randparams->psi->min+(REAL8)vector->data[k]*(randparams->psi->max-randparams->psi->min);
    k++;
  }
  if (randparams->cosiota!=NULL) {
    randparams->cosiota->value=randparams->cosiota->min+(REAL8)vector->data[k]*(randparams->cosiota->max-randparams->cosiota->min);
    k++;
  }
  if (randparams->det!=NULL) {
    randparams->det->value=floor((REAL8)randparams->det->min+(REAL8)vector->data[k]*((REAL8)randparams->det->max-(REAL8)randparams->det->min));
    if ((randparams->det->value<1.0)&&(randparams->det->value>=0.0)) sprintf(randparams->detector,"LLO");
    else if ((randparams->det->value<2.0)&&(randparams->det->value>=1.0)) sprintf(randparams->detector,"LHO");
    k++;
  }
    
  /* clean up memory */
  if (p>0) {
    LALDestroyRandomParams(&status,&params);
    LALDestroyVector(&status,&vector);
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
  raMIN=0.0; 
  raMAX=LAL_TWOPI; 
  decMIN=-LAL_PI/2.0; 
  decMAX=LAL_PI/2.0; 
 
  phiflag=0; /* p */
  psiflag=0; /* P */
  cosiotaflag=0; /* c */
  allskyflag=0;
  detflag=0;

  h0MIN=1.0; /* g */
  h0MAX=1.0; /* G */
  freqMIN=600.0; /* f */
  freqMAX=600.0; /* F */
  freqband=0.0;
  
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
  sprintf(primarytemplatefile," ");
  sprintf(secondarytemplatefile," ");

  primarystart=731100000; /* w */
  secondarystart=731100000; /* W */
 
  coflag=0;
  sprintf(infile," "); /* i */
  sprintf(primaryoutputfile," "); /* o */
  sprintf(secondaryoutputfile," ");
  sprintf(primarystampsfile," ");  /* m */
  sprintf(secondarystampsfile," ");
  sprintf(primarynoisedir," "); 
  sprintf(secondarynoisedir," ");
  sprintf(primarysftbase," ");  
  sprintf(secondarysftbase," ");
  seed=0; /* z */

  {
    int option_index = 0;
    static struct option long_options[] = {
      {"ptemplatefile", required_argument, 0, 'S'},
      {"stemplatefile", required_argument, 0, 's'},
      {"pstart", required_argument, 0, 'A'},
      {"sstart", required_argument, 0, 'a'},
      {"phi", no_argument, 0, 'p'},
      {"psi", no_argument, 0, 'P'},
      {"cosiota", no_argument, 0, 'c'},
      {"allsky", no_argument, 0, 'j'},
      {"det", no_argument, 0, 'J'},
      {"h0MIN", required_argument, 0, 'g'},
      {"h0MAX", required_argument, 0, 'G'},
      {"f0MIN", required_argument, 0, 'f'},
      {"f0MAX", required_argument, 0, 'F'},
      {"fband", required_argument, 0, 'b'},
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
      {"pstampsfile", required_argument, 0, 'm'},
      {"sstampsfile", required_argument, 0, 'M'},
      {"pnoisedir", required_argument, 0, 'n'},
      {"snoisedir", required_argument, 0, 'N'},
      {"psftbase", required_argument, 0, 'd'},
      {"ssftbase", required_argument, 0, 'D'},
      {"infile", required_argument, 0, 'i'},
      {"co", no_argument, 0, 'C'},
      {"poutfile", required_argument, 0, 'o'},
      {"soutfile", required_argument, 0, 'x'},
      {"seed", required_argument, 0, 'z'},
      {"help", no_argument, 0, 'h'}
    };
    /* Scan through list of command line arguments */
    while (!errflg && ((c = getopt_long (argc, argv,"hS:s:A:a:pPcJjg:G:f:F:b:r:R:q:Q:t:T:e:E:k:K:m:M:n:N:d:D:i:Co:x:z:",long_options, &option_index)))!=-1)
      switch (c) {
      case 'S':
	temp=optarg;
	sprintf(primarytemplatefile,temp);
	bintempfileflag=1;
	break;
      case 's':
	temp=optarg;
	sprintf(secondarytemplatefile,temp);
	break;
      case 'A':
	primarystart=atoi(optarg);
	break;
      case 'a':
	secondarystart=atoi(optarg);
	break;
      case 'p':
	phiflag=1;
	break;
      case 'P':
	psiflag=1;
	break;
      case 'c':
	cosiotaflag=1;
	break;
      case 'j':
	allskyflag=1;
	break;
      case 'J':
	detflag=1;
	break;
      case 'g':
	h0MIN=atof(optarg);
	h0flag=1;
	break;
      case 'G':
	h0MAX=atof(optarg);
	h0flag=1;
	break;
      case 'f':
	freqMIN=atof(optarg);
	freqflag=1;
	break;
      case 'F':
	freqMAX=atof(optarg);
	break;
      case 'b':
	freqband=atof(optarg);
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
      case 'z':
	seed=atoi(optarg);
	break;
      case 'C':
	coflag=1;
	break;
      case 'm':
	temp=optarg;
	sprintf(primarystampsfile,temp);
	break; 
      case 'M':
	temp=optarg;
	sprintf(secondarystampsfile,temp);
	break; 
      case 'n':
	temp=optarg;
	sprintf(primarynoisedir,temp);
	break; 
      case 'N':
	temp=optarg;
	sprintf(secondarynoisedir,temp);
	break; 
      case 'd':
	temp=optarg;
	sprintf(primarysftbase,temp);
	break; 
      case 'D':
	temp=optarg;
	sprintf(secondarysftbase,temp);
	break; 
      case 'i':
	temp=optarg;
	sprintf(infile,temp);
	break;
      case 'o':
	temp=optarg;
	sprintf(primaryoutputfile,temp);
	break;
      case 'x':
	temp=optarg;
	sprintf(secondaryoutputfile,temp);
	break;
      case 'h':
	/* print usage/help message */
	fprintf(stdout,"Arguments are:\n");
	fprintf(stdout,"\t--ptemplatefile    STRING\t The name of the source file [DEFAULT=]\n");
	fprintf(stdout,"\t--stemplatefile    STRING\t The name of the source [DEFAULT=]\n");
	fprintf(stdout,"\t--phi              BOOLEAN\t Set if random initial phase required  [DEFAULT=false]\n");
	fprintf(stdout,"\t--psi              BOOLEAN\t Set if random polarisation angle required  [DEFAULT=false]\n");
	fprintf(stdout,"\t--cosiota          BOOLEAN\t Set if random Cos(iota) required [DEFAULT=false]\n");
	fprintf(stdout,"\t--h0MIN            REAL8\t MIN Gravitational wave amplitude [DEFAULT=1.0]\n");
	fprintf(stdout,"\t--h0MAX            REAL8\t MAX Gravitational wave amplitude [DEFAULT=1.0]\n");
	fprintf(stdout,"\t--f0MIN            REAL8\t MIN Gravitational wave frequency in Hz [DEFAULT=600.0]\n");
	fprintf(stdout,"\t--f0MAX            REAL8\t MAX Gravitational wave frequency in Hz [DEFAULT=600.0]\n");    
	fprintf(stdout,"\t--fband            REAL8\t The bnadwidth of data to be genereted [DEFAULT=0.0]\n"); 
	fprintf(stdout,"\t--smaMIN           REAL8\t Projected semi-major axis of orbit in seconds [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--smaMAX           REAL8\t Projected semi-major axis of orbit in seconds [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--periodMIN        REAL8\t MIN Period of orbit in seconds [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--periodMAX        REAL8\t MAX Period of orbit in seconds [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--tperiMIN         INT4\t MIN Observed time of periapse passage in seconds [DEFAULT=0] \n");
	fprintf(stdout,"\t--tperiMAX         INT4\t MAX Observed time of periapse passage in nanoseconds [DEFAULT=0]\n");
	fprintf(stdout,"\t--eccMIN           REAL8\t MIN Orbital eccentricity [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--eccMAX           REAL8\t MAX Orbital eccentricity [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--argpMIN          REAL8\t MIN Argument of orbital periapse in radians [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--argpMAX          REAL8\t MAX Argument of orbital periapse in radians [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--seed             INT4\t Seed for random number generation [DEFAULT=0 (from clock)]\n");
	fprintf(stdout,"\t--infile           STRING\t Name of input configuration file [DEFAULT=NULL]\n");
	fprintf(stdout,"\t--co               BOOLEAN\t Set if coincidence analysis is being performed [DEFAULT=]\n"); 
	fprintf(stdout,"\t--pstampsfile      STRING\t Name of output timstamps file for the primary detector [DEFAULT=NULL]\n"); 
	fprintf(stdout,"\t--sstampsfile      STRING\t Name of output timstamps file for the secondary detector [DEFAULT=NULL]\n"); 
	fprintf(stdout,"\t--pnoisedir        STRING\t Location of the primary data to be injected into [DEFAULT=NULL]\n"); 
	fprintf(stdout,"\t--snoisedir        STRING\t Location of the secondary data to be injected into [DEFAULT=NULL]\n"); 
	fprintf(stdout,"\t--psftbase         STRING\t Location of the primary output dataset + sft basename [DEFAULT=NULL]\n"); 
	fprintf(stdout,"\t--ssftbase         STRING\t Location of the secondary output dataset + sft basename [DEFAULT=NULL]\n"); 
	fprintf(stdout,"\t--poutfile         STRING\t Name of output random configuration file for the primary detector [DEFAULT=NULL]\n");
	fprintf(stdout,"\t--soutfile         STRING\t Name of output random configuration file for the secondary detector [DEFAULT=NULL]\n");
	exit(0);
	break;
      default:
	/* unrecognized option */
	errflg++;
	fprintf(stderr,"Unrecognized option argument %c\n",c);
	exit(1);
	break;
      }

    if ((detflag)&&(coflag)) {
      printf("ERROR : can't select random detector and coincidence analysis !!!\n");
      exit(1);
    }
    
  }
  
  /* update global variable and return */
  return errflg;
}


/*******************************************************************************/

int OutputRandomConfigFile(char *outfile,RandomParameters *random,RandomParameters *rand_det) 
{

  /* this function scans through the original input template config file and */
  /* replaces the appropriate variables with the random values */

  FILE *fpin,*fpout,*fppython;
  REAL8 aplus,across,h0;
  INT4 errcode;
  char name[256],line[256];
  char phi0text[256];
  char psitext[256];
  char longitudetext[256];
  char latitudetext[256];
  char ifotext[256];
  char durationtext[256];
  char tsfttext[256];
  char tstarttext[256];
  char aPlustext[256];
  char aCrosstext[256];
  char f0text[256];
  char fbandtext[256];
  char fmintext[256];
  char orbitSemiMajorAxistext[256];
  char orbitPeriodtext[256];
  char orbitEccentricitytext[256];
  char orbitArgPeriapsetext[256];
  char orbitTperiSSBsectext[256];
  char orbitTperiSSBnstext[256];
  char stampstext[256];
  char sftbasetext[256];
  char noisedirtext[256];
  char reftimetext[256];
  INT4 i;

  /* defining the actual strings to identify in the input file */
  LALSnprintf(phi0text,256,"phi0");
  LALSnprintf(psitext,256,"psi");
  LALSnprintf(longitudetext,256,"longitude");
  LALSnprintf(latitudetext,256,"latitude");
  LALSnprintf(ifotext,256,"detector");
  LALSnprintf(durationtext,256,"duration");
  LALSnprintf(tsfttext,256,"Tsft");
  LALSnprintf(tstarttext,256,"startTime");
  LALSnprintf(reftimetext,256,"refTime");
  LALSnprintf(aPlustext,256,"aPlus");
  LALSnprintf(aCrosstext,256,"aCross");
  LALSnprintf(f0text,256,"f0");
  LALSnprintf(fbandtext,256,"Band");
  LALSnprintf(fmintext,256,"fmin");
  LALSnprintf(stampstext,256,"timestampsFile");
  LALSnprintf(orbitSemiMajorAxistext,256,"orbitSemiMajorAxis");
  LALSnprintf(orbitPeriodtext,256,"orbitPeriod");
  LALSnprintf(orbitEccentricitytext,256,"orbitEccentricity");
  LALSnprintf(orbitArgPeriapsetext,256,"orbitArgPeriapse");
  LALSnprintf(orbitTperiSSBsectext,256,"orbitTperiSSBsec");
  LALSnprintf(orbitTperiSSBnstext,256,"orbitTperiSSBns");
  LALSnprintf(sftbasetext,256,"outSFTbname");
  LALSnprintf(noisedirtext,256,"noiseSFTs");

  /* opening the input config file */
  fpin=fopen(infile,"r");
  if (fpin==NULL) {
    fprintf(stderr,"Unable to open file %s\n",infile);
    return 1;
  }
  /* opening the output config file */
  fpout=fopen(outfile,"w");
  if (fpout==NULL) {
    fprintf(stderr,"Unable to open file %s\n",outfile);
    return 1;
  }

  /* here we are opening another output file that contains the same output */
  /* parameters but in a python readable config file */
  fppython=fopen("randparams_python.data","w");
  if (fppython==NULL) {
    fprintf(stderr,"Unable to open file\n");
    return 1;
  }
  fprintf(fppython,"[randomparameters]\n");
  
  /* define h0 value here */
  if (random->h0==NULL) h0=1.0;
  else h0=random->h0->value;


  i=0;
  /* read in the input file and identify the appropriate lines */ 
  while (fgets(line,255,fpin))  
    {

    /* get the first column string to compare */
    errcode=sscanf(line,"%s",name); 
    /* fprintf(stdout,"code is %d name is #%s#\n",errcode,name); */
       

    if (errcode>0) {
      /* if any of the following strings are found and are to be replaced the change line */
      if (strcmp(name,ifotext)==0) {
	sprintf(line,"detector\t= %s\t\t# Detector: LHO, LLO, VIRGO, GEO, TAMA, CIT, ROME\n",rand_det->detector);
	fprintf(fppython,"det = %s\n",rand_det->detector); 
      }
      if (strcmp(name,reftimetext)==0) {
        sprintf(line,"refTime\t= %d\t\t# ?\n",rand_det->start);
      }
      else if ((strcmp(name,longitudetext)==0)&&(random->ra!=NULL)) {
	sprintf(line,"longitude\t= %6.12f\t# source longitude (in_radians)\n",random->ra->value);
	fprintf(fppython,"ra = %6.12f\n",random->ra->value); 
      }
      else if ((strcmp(name,latitudetext)==0)&&(random->dec!=NULL)) {
	sprintf(line,"latitude\t= %6.12f\t# source latitude (in radians)\n",random->dec->value);
	fprintf(fppython,"declination = %6.12f\n",random->dec->value); 
      }
      else if ((strcmp(name,aPlustext)==0)&&(random->cosiota!=NULL)) {
	aplus=0.5*h0*(1.0+random->cosiota->value*random->cosiota->value);
	sprintf(line,"aPlus\t\t= %6.12e\t# plus-polarization amplitude a_+ (strain)\n",aplus);
	fprintf(fppython,"cosiota = %6.12f\n",random->cosiota->value);
      }
      else if ((strcmp(name,aCrosstext)==0)&&(random->cosiota!=NULL)) {
	across=h0*random->cosiota->value;
	sprintf(line,"aCross\t\t= %6.12e\t# cross-polarization amplitude a_x (strain)\n",across);
      }
      else if ((strcmp(name,psitext)==0)&&(random->psi!=NULL)) {
	sprintf(line,"psi\t\t= %6.12f\t# wave polarization angle Psi\n",random->psi->value);
	fprintf(fppython,"psi = %6.12f\n",random->psi->value);
      }
      else if ((strcmp(name,phi0text)==0)&&(random->phi!=NULL)) {
	sprintf(line,"phi0\t\t= %6.12f\t# initial wave-phase phi0 (at reference-time tRef)\n",random->phi->value);
	fprintf(fppython,"phi = %6.12f\n",random->phi->value);
      }
      else if ((strcmp(name,f0text)==0)&&(random->freq!=NULL)) {
	sprintf(line,"f0\t\t= %6.12f\t# intrinsic signal frequency f0 (at tRef)\n",random->freq->value);
	fprintf(fppython,"f0 = %6.12f\n",random->freq->value); 
      }
      else if ((strcmp(name,fmintext)==0)&&(random->freq!=NULL)) {
	f_min=(REAL8)floor(random->freq->value-((REAL8)freqband/2.0)+0.5);
	sprintf(line,"fmin\t\t= %6.12f\t# Lowest frequency in output SFT (= heterodyning frequency)\n",f_min);
      }
      else if ((strcmp(name,fbandtext)==0)&&(random->freq!=NULL)) {
	sprintf(line,"Band\t\t= %6.12f\t# bandwidth of output SFT in Hz (= 1/2 sampling frequency)\n",freqband);
      }
      else if (strcmp(name,tstarttext)==0) {
	sprintf(line,"startTime\t= %d\t# GPS start time of (contiguous) output time-series\n",rand_det->start);
	fprintf(fppython,"tstart = %d\n",rand_det->start); 
      }
      else if (strcmp(name,stampstext)==0) {
	sprintf(line,"timestampsFile\t= %s\t# Timestamps file\n",rand_det->stampsfile);
      }
      else if (strcmp(name,sftbasetext)==0) {
	sprintf(line,"outSFTbname\t= %s\t# Path and basefilename of output SFT files\n",rand_det->sftbase);
      }
      else if (strcmp(name,noisedirtext)==0) {
	/* note that we put the * on the end to bypass wildcard evaluation problems in the python script */
	sprintf(line,"noiseSFTs\t= %s/*\t# Glob-like pattern specifying noise-SFTs to be added to signal\n",rand_det->noisedir);
      }
      else if ((strcmp(name,orbitSemiMajorAxistext)==0)&&(random->sma!=NULL)) {
	sprintf(line,"orbitSemiMajorAxis\t= %6.12f\t# Projected orbital semi-major axis a in seconds (i.e. a*sin(i)/c)\n",random->sma->value);
	fprintf(fppython,"sma = %6.12f\n",random->sma->value); 
      }
      else if ((strcmp(name,orbitEccentricitytext)==0)&&(random->ecc!=NULL)) {
	sprintf(line,"orbitEccentricity\t= %6.12f\t# Orbital eccentricity\n",random->ecc->value);
	fprintf(fppython,"ecc = %6.12f\n",random->ecc->value); 
      }
      else if ((strcmp(name,orbitTperiSSBsectext)==0)&&(random->tperi!=NULL)) {
	sprintf(line,"orbitTperiSSBsec\t= %d\t\t# 'observed' (SSB) time of periapsis passage. Seconds.\n",random->tperi->value.gpsSeconds);
	fprintf(fppython,"tpsec = %d\n",random->tperi->value.gpsSeconds); 
      }
      else if ((strcmp(name,orbitTperiSSBnstext)==0)&&(random->tperi!=NULL)) {
	sprintf(line,"orbitTperiSSBns\t\t= %d\t\t# 'observed' (SSB) time of periapsis passage. Nanoseconds.\n",random->tperi->value.gpsNanoSeconds);
	fprintf(fppython,"tpnano = %d\n",random->tperi->value.gpsNanoSeconds); 
      }
      else if ((strcmp(name,orbitPeriodtext)==0)&&(random->period!=NULL)) {
	sprintf(line,"orbitPeriod\t\t= %6.12f\t# Orbital period (seconds)\n",random->period->value);
	fprintf(fppython,"period = %6.12f\n",random->period->value); 
      }
      else if ((strcmp(name,orbitArgPeriapsetext)==0)&&(random->argp!=NULL)) {
	sprintf(line,"orbitArgPeriapse\t= %6.12f\t# Argument of periapsis (radians)\n",random->argp->value); 
	fprintf(fppython,"argp = %6.12f\n",random->argp->value); 
      }
      
      /* output the appropriate line to the output file */
      fprintf(fpout,"%s",line);
    }
  }
  
  fclose(fpin);
  fclose(fpout);
  fclose(fppython);

  return 0;
  
}
/*******************************************************************************/

