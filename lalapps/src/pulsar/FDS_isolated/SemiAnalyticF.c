/*********************************************************************************/
/*                 Semi-Analytic calculation of the F-statistic                  */
/*                                                                               */
/*			           X. Siemens                                    */
/*                                                                               */
/*                             UWM - February 2003                               */
/*********************************************************************************/

#include <errno.h>
#include "SemiAnalyticF.h"

LIGOTimeGPS *timestamps=NULL;       /* Time stamps from SFT data */
INT4 lalDebugLevel=3;
static LALStatus status;
AMCoeffs amc;
BOOLEAN stampsflag;

extern char *optarg;
extern int optind, opterr, optopt;

int main(int argc,char *argv[]) 
{

  /* Reads command line arguments into the CommandLineArgs struct. 
     In the absence of command line arguments it sets some defaults */
  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;

  if (stampsflag) {
    if (ReadTimeStamps(CommandLineArgs)) return 2;
  }
  else {
    if (MakeTimeStamps(CommandLineArgs)) return 3;
  }

  if (ComputeF(CommandLineArgs)) return 5;

  if (Freemem()) return 8;

  return 0;

}

/*******************************************************************************/

int MakeTimeStamps(struct CommandLineArgsTag CLA) 
{
  int i;
 
  /* allocate memory for timestamps */
  timestamps=(LIGOTimeGPS *)LALMalloc(CLA.nTsft*sizeof(LIGOTimeGPS)); 
      
  /* generate timetamps */
  for (i=0;i<CLA.nTsft;i++){
    timestamps[i].gpsSeconds=CLA.gpsStart+(int)(i*CLA.tsft);
    timestamps[i].gpsNanoSeconds=0;
  } 
  
  return 0;
  
}
/*******************************************************************************/

int ReadTimeStamps(struct CommandLineArgsTag CLA) 
{
  FILE *fp;
  int i,r;
 
 
  /*   %strcpy(filename,inDataFilename); */
  fp=fopen(CLA.timestamps,"r");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",CLA.timestamps);
    return 1;
  }
  timestamps=(LIGOTimeGPS *)LALMalloc(CLA.nTsft*sizeof(LIGOTimeGPS)); 
      
  
  for (i=0;i<CLA.nTsft;i++){
    r=fscanf(fp,"%d  %d\n", &timestamps[i].gpsSeconds, &timestamps[i].gpsNanoSeconds);
    if ( r !=2 ) {
      fprintf(stderr,"Unable to read datum # %d\n",i);
      fprintf(stderr,"from file %s\n",CLA.timestamps);
      return 1; 
    } 
  } 
  
  fclose(fp);
  return 0;
  
}
/*******************************************************************************/

int ComputeF(struct CommandLineArgsTag CLA)
{

  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  EphemerisData *edat=NULL;          /* Stores earth/sun ephemeris data for barycentering */
  LALDetector Detector;              /* Our detector*/
  EarthState earth;
  AMCoeffsParams *amParams;
  LIGOTimeGPS *midTS=NULL;           /* Time stamps for amplitude modulation coefficients */
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4 leap;

  REAL8 A,B,C,D,alpha,beta,delta,kappa,A1,A2,A3,A4,h0,cosi,psi,phi,To,Sh,F;

  char filenameE[256],filenameS[256];

  INT4 k;
  
  strcpy(filenameE,CLA.efiles);
  strcat(filenameE,"/earth00-04.dat");

  strcpy(filenameS,CLA.efiles);
  strcat(filenameS,"/sun00-04.dat");

  edat=(EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = filenameE;     
  (*edat).ephiles.sunEphemeris = filenameS;         

  LALLeapSecs(&status,&leap,&timestamps[0],&formatAndAcc);
  (*edat).leap=leap; 

  LALInitBarycenter(&status, edat);               /* Reads in ephemeris files */

  if(CLA.detector == 1)
    {
      Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];   
    }
  if(CLA.detector == 2)
    {
      Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];   
    }

  if(CLA.detector  == 3)
    {
        if (CreateDetector(&Detector)) return 5;
    }
/* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.site.location[0]=Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=Detector.location[2]/LAL_C_SI;
  baryinput.alpha=CLA.skyalpha;
  baryinput.delta=CLA.skydelta;
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
  amParams->edat = edat;
  amParams->das->pDetector = &Detector; 
  amParams->das->pSource->equatorialCoords.latitude = CLA.skydelta;
  amParams->das->pSource->equatorialCoords.longitude = CLA.skyalpha;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!!!!!!!!*/
  amParams->leapAcc=formatAndAcc.accuracy;

/* Allocate space for AMCoeffs */
  amc.a = NULL;
  amc.b = NULL;
  LALSCreateVector(&status, &(amc.a), (UINT4)  CLA.nTsft);
  LALSCreateVector(&status, &(amc.b), (UINT4)  CLA.nTsft);

 /* Mid point of each SFT */
   midTS = (LIGOTimeGPS *)LALCalloc(CLA.nTsft,sizeof(LIGOTimeGPS));
   for(k=0; k<CLA.nTsft; k++)
     {
       REAL8 teemp=0.0;
       LALGPStoFloat(&status,&teemp, &(timestamps[k]));
       teemp += 0.5*CLA.tsft;
       LALFloatToGPS(&status,&(midTS[k]), &teemp);
     }
   
   LALComputeAM(&status, &amc, midTS, amParams);

/*    fprintf(stdout,"A = %f\n",amc.A); */
/*    fprintf(stdout,"B = %f\n",amc.B); */
/*    fprintf(stdout,"C = %f\n",amc.C); */
/*    fprintf(stdout,"D = %f\n",amc.C); */

   A=amc.A;
   B=amc.B;
   C=amc.C;
   D=amc.D; 

   cosi=CLA.cosiota;
   phi=CLA.phi;
   psi=CLA.psi;
   h0=CLA.h0;

   A1=h0*( 0.5*(1+pow(cosi,2))*cos(2.0*psi)*cos(2*phi)-cosi*sin(2.0*psi)*sin(2.0*phi));
   A2=h0*( 0.5*(1+pow(cosi,2))*sin(2.0*psi)*cos(2*phi)+cosi*cos(2.0*psi)*sin(2.0*phi));
   A3=h0*(-0.5*(1+pow(cosi,2))*cos(2.0*psi)*sin(2*phi)-cosi*sin(2.0*psi)*cos(2.0*phi));
   A4=h0*(-0.5*(1+pow(cosi,2))*sin(2.0*psi)*sin(2*phi)+cosi*cos(2.0*psi)*cos(2.0*phi));

   alpha = 1.0/2.0*A*A1+1.0/2.0*C*A2;
   beta  = 1.0/2.0*B*A2+1.0/2.0*C*A1;
   delta = 1.0/2.0*A*A3+1.0/2.0*C*A4;
   kappa = 1.0/2.0*B*A4+1.0/2.0*C*A3;

   To=CLA.nTsft*CLA.tsft;

   Sh=pow(CLA.sqrtSh,2);

   F=((B*pow(alpha,2)+A*pow(beta,2)-2*C*alpha*beta)/D+(B*pow(delta,2)+A*pow(kappa,2)-2*C*delta*kappa)/D)*To/Sh;

   fprintf(stdout,"%e\n",F);
/*    fprintf(stdout,"\nSNR = %e\n",sqrt(F)); */
	

   LALFree(midTS);

   LALFree(amParams->das->pSource);
   LALFree(amParams->das);
   LALFree(amParams);
   
   LALFree(edat->ephemE);
   LALFree(edat->ephemS);
   LALFree(edat);

   return 0;
}


/*******************************************************************************/

int CreateDetector(LALDetector *Detector){

/*   LALDetector Detector;  */
  LALFrDetector detector_params;
  LALDetectorType bar;
  LALDetector Detector1;

/*   detector_params=(LALFrDetector )LALMalloc(sizeof(LALFrDetector)); */
 
  bar=LALDETECTORTYPE_CYLBAR;
  strcpy(detector_params.name,"NAUTILUS");
  detector_params.vertexLongitudeRadians=12.67*LAL_PI/180.0;
  detector_params.vertexLatitudeRadians=41.82*LAL_PI/180.0;
  detector_params.vertexElevation=300.0;
  detector_params.xArmAltitudeRadians=0.0;
  detector_params.xArmAzimuthRadians=44.0*LAL_PI/180.0;

  LALCreateDetector(&status,&Detector1,&detector_params,bar);

  *Detector=Detector1;

  return 0;
}

/*******************************************************************************/

   int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{
  INT4 c, errflg = 0;
  optarg = NULL;
  
  /* Initialize default values */
  CLA->skyalpha=0.0;
  CLA->skydelta=0.0;
  CLA->tsft=0.0;
  CLA->detector=0;               
  CLA->nTsft=0;            
  CLA->timestamps=NULL;
  CLA->gpsStart=-1;
  stampsflag=0;
  CLA->efiles=NULL;
  CLA->phi=0.0;
  CLA->psi=0.0;
  CLA->cosiota=0.0;
  CLA->h0=1.0;
  CLA->sqrtSh=1.0;

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"ha:d:E:n:T:S:D:t:Q:Y:i:s:N:"))!=-1))
    switch (c) {
    case 'T':
      /* timestamps file */
      CLA->timestamps=optarg;
      stampsflag=1;
      break;
    case 'S':
      /* start time of observation */
      CLA->gpsStart=atoi(optarg);
      break;
    case 'E':
      /* ephemeris file location */
      CLA->efiles=optarg;
      break;
    case 'a':
      /* sky position alpha */
      CLA->skyalpha=atof(optarg);
      break;
    case 'd':
      /* sky position delta */
      CLA->skydelta=atof(optarg);
      break;
    case 'Q':
      /* phi0 - initial phase */
      CLA->phi=atof(optarg);
      break;
    case 'Y':
      /* psi - polarisation */
      CLA->psi=atof(optarg);
      break;
    case 'i':
      /* iota */
      CLA->cosiota=atof(optarg);
      break;
    case 's':
      /* strain */
      CLA->h0=atof(optarg);
      break;
    case 'N':
      /* sqrt(Sh) */
      CLA->sqrtSh=atof(optarg);
      break;
    case 't':
      /* duration fo an sft */
      CLA->tsft=atof(optarg);
      break;
    case 'n':
      /* number of SFTs */
      CLA->nTsft=atof(optarg);
      break;
    case 'D':
      /* detector number */
      CLA->detector=atof(optarg);
      break;

    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\t-a\tFLOAT\t Sky position alpha (equatorial coordinates) in radians \n");
      fprintf(stdout,"\t-d\tFLOAT\t Sky position delta (equatorial coordinates) in radians \n");
      fprintf(stdout,"\t-Q\tFLOAT\t Initial phase in radians (default 0.0)\n");
      fprintf(stdout,"\t-Y\tFLOAT\t Polarisation in radians (default 0.0)\n");
      fprintf(stdout,"\t-i\tFLOAT\t Cos(iota) (default 0.0)\n");
      fprintf(stdout,"\t-s\tFLOAT\t Strain (default 1.0)\n");
      fprintf(stdout,"\t-N\tFLOAT\t Noise floor (sqrt(Sh)) in 1/sqrt(Hz) (default 1.0)\n");
      fprintf(stdout,"\t**** Only one of the following two needs to be set ****\n");
      fprintf(stdout,"\t-T\tSTRING\t Name of timestamps file \n");
      fprintf(stdout,"\t-S\tINTEGER\t GPS start time of continuous observation\n");
      fprintf(stdout,"\t-t\tFLOAT\t Time baseline of SFTs \n");
      fprintf(stdout,"\t-n\tINTEGER\t Number of data points to read from timestamps file \n");
      fprintf(stdout,"\t-E\tSTRING\t Directory where Ephemeris files are located \n");
      fprintf(stdout,"\t-D\tINTEGER\t Detector number (1=LLO, 2=LHO, 3=Roman Bar) \n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }

  if(CLA->efiles == NULL)
    {
      fprintf(stderr,"No ephemeris data (earth??.dat, sun??.dat) directory specified; input directory with -E option.\n");
      fprintf(stderr,"Try ./lalapps_SemiAnalyticF -h \n");
      return 1;
    }      
  if(CLA->detector == 0)
    {
      fprintf(stderr,"No detector specified; set with -D option.\n");
      fprintf(stderr,"Try ./lalapps_SemiAnalyticF -h \n");
     return 1;
    }      
  if((CLA->gpsStart>=0)&&(CLA->timestamps!=NULL))
    {
      fprintf(stderr,"Both start time and timestamps file specified - just need one !!\n");
      fprintf(stderr,"Try ./lalapps_SemiAnalyticF -h \n");
     return 1;
    }   
  if((CLA->gpsStart<0)&&(CLA->timestamps==NULL))
    {
      fprintf(stderr,"Need to specify gpsStart time or a timestamps file !!\n");
      fprintf(stderr,"Try ./lalapps_SemiAnalyticF -h \n");
     return 1;
    }   


  /* update global variable and return */
  return errflg;
}


/*******************************************************************************/

int Freemem() 
{

  /*Free timestamps*/
  LALFree(timestamps);

  /*Free DemodParams*/
  LALSDestroyVector(&status, &(amc.a));
  LALSDestroyVector(&status, &(amc.b));

  LALCheckMemoryLeaks();
  
  return 0;
}
