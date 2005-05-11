/*********************************************************************************/
/*         Program to read in a 2-D binary orbit template file and find          */
/*                 the closest N templates in parameter space                    */
/*                                                                               */
/*			           C. Messenger                                  */
/*                                                                               */
/*                         BIRMINGHAM UNIVERISTY -  2004                         */
/*********************************************************************************/


#include "GenerateBinaryMesh_v1.h"

INT4 lalDebugLevel=3;
static LALStatus status;
REAL8 sma0,period,ecc,argp;
REAL8 *sma,*dist;
LIGOTimeGPS tperi0,*tperi,TstartDET,TstartSSB,tperiCURRENT;
REAL8 f_max,tspan,ecc,argp,period;
INT4 tperisecMIN,tperisecMAX,tperinsMIN,tperinsMAX,TstartsecMIN,TstartsecMAX;
REAL8 smaMIN,smaMAX,mismatch,sma0FULL;
REAL8 gamXX,gamXY,gamYY;
REAL8 X,Y,X0,Y0; 
REAL8 RA,dec;
REAL8 smaCURRENT;
INT4 det,Nsub,Nfull;
INT4 tperisec0FULL,tperins0FULL;
CHAR ifo[256];
REAL8 RA,dec;
FILE *sbfp,*fbfp;
LIGOTimeGPS TstartSSB;
LIGOTimeGPS TstartDET;
CHAR fullbankfile[256],subbankfile[256],ephemfile[256],yr[256],dmp[256];
LIGOTimeGPS tperiGPS;
EphemerisData *edat=NULL;          /* Stores earth/sun ephemeris data for barycentering */
LALDetector Detector;              /* Our detector*/
EarthState earth;
EmissionTime emit;
LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
INT4 leap;
BinaryMeshFileHeader BMFheader;
INT4 exactflag;

extern char *optarg;
extern int optind, opterr, optopt;

int ReadCommandLine(int argc,char *argv[]);
int ReadFullBank(void);
int OutputSortedDist(void);
int FreeMem(void);
int SetupBaryInput(void);
int GetSSBTime(LIGOTimeGPS *, LIGOTimeGPS *);
int CalculateDistance(REAL8 *, REAL8 *, REAL8 *);
int OutputSortedDist(void);
int compare(const void *, const void *);

int main(int argc,char *argv[]) 
{
 
 
  if (ReadCommandLine(argc,argv)) return 1;
 
  if (ReadFullBank()) return 3;
 
  if (OutputSortedDist()) return 4;
   
  if (FreeMem()) return 5;
  
  return 0;

}
/****************************************************************************/

int GetSSBTime(LIGOTimeGPS *tdet, LIGOTimeGPS *tssb)
{

  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */

 /* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.tgps.gpsSeconds=tdet->gpsSeconds;
  baryinput.tgps.gpsNanoSeconds=tdet->gpsNanoSeconds;
  baryinput.site.location[0]=Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=Detector.location[2]/LAL_C_SI;
  baryinput.alpha=RA;
  baryinput.delta=dec;
  baryinput.dInv=0.e0;

  /* Setup barycentric time at start of observation (first timestamp) */
  LALBarycenterEarth(&status,&earth,tdet,edat);
  LALBarycenter(&status,&emit,&baryinput,&earth);
  tssb->gpsSeconds=(emit.te.gpsSeconds);
  tssb->gpsNanoSeconds=emit.te.gpsNanoSeconds;
  
  return 0;
}

/*******************************************************************************/

int SetupBaryInput(void)
{

  CHAR filenameE[256],filenameS[256];
  FILE *fp;
  
  /* make the full file name/location for the ephemeris files */
  strcpy(filenameE,ephemfile);
  strcat(filenameE,"/earth");
  strcat(filenameE,yr);
  strcat(filenameE,".dat");
  strcpy(filenameS,ephemfile);
  strcat(filenameS,"/sun");
  strcat(filenameS,yr);
  strcat(filenameS,".dat");

  /* open the ephemeris files to check they exist */
  fp=fopen(filenameE,"r");
  if (fp==NULL) 
    {
      fprintf(stderr,"Could not find %s\n",filenameE);
      return 1;
    }
  fclose(fp);
  fp=fopen(filenameS,"r");
  if (fp==NULL) 
    {
      fprintf(stderr,"Could not find %s\n",filenameS);
      return 1;
    }
  fclose(fp);

  /* allocate memory for the ephemeris */
  edat=(EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = filenameE;     
  (*edat).ephiles.sunEphemeris = filenameS;         

  /* set up leap second information */
  LALLeapSecs(&status,&leap,&TstartDET,&formatAndAcc);
  (*edat).leap=leap;

  /* Read in ephemeris files */
  LALInitBarycenter(&status,edat);             

  /* setup chosen detector */
  if(strcmp(ifo,"GEO")!=0) Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if(strcmp(ifo,"LLO")!=0) Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(strcmp(ifo,"LHO")!=0) Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
 
  /* First job is to convert the observation start time from detector time to SSB time */
  if (GetSSBTime(&TstartDET,&TstartSSB)) return 1;

  return 0;

}
            
/*******************************************************************************/

int ReadFullBank(void) 
{

  RTPLocation RTPloc;
  XYLocation XYloc;
  INT4 i;
  REAL8 tperi0FLT;
  REAL8 tperiMINFLT;
  REAL8 tperiMAXFLT;
  REAL8 dummy;
  REAL8 dist_min=999999999.0;
  REAL8 Xmin,Ymin;
 
  /* here we read in the full template file and its header information */
  
  /* open the full bank file */
  fbfp=fopen(fullbankfile,"r");
  if (fbfp==NULL) {
    fprintf(stderr,"Unable to open file %s\n",fullbankfile);
    return 1;
  }

  if (ReadMeshFileHeader(fbfp,&BMFheader)) return 1;

  /* put all the output into the correct variables */
  f_max=BMFheader.f_max;
  tspan=BMFheader.tspan;
  TstartDET.gpsSeconds=BMFheader.tstart.gpsSeconds;
  TstartDET.gpsNanoSeconds=BMFheader.tstart.gpsNanoSeconds;
  Nfull=BMFheader.Nfilters;
  mismatch=BMFheader.mismatch;
  sma0FULL=BMFheader.sma_0;
  smaMIN=BMFheader.sma_MIN;
  smaMAX=BMFheader.sma_MAX;
  tperisec0FULL=BMFheader.tperi_0.gpsSeconds;
  tperins0FULL=BMFheader.tperi_0.gpsNanoSeconds;
  tperisecMIN=BMFheader.tperi_MIN.gpsSeconds;
  tperinsMIN=BMFheader.tperi_MIN.gpsNanoSeconds;
  tperisecMAX=BMFheader.tperi_MAX.gpsSeconds;
  tperinsMAX=BMFheader.tperi_MAX.gpsNanoSeconds;
  ecc=BMFheader.ecc_MIN;
  argp=BMFheader.argp_MIN;
  period=BMFheader.period_MIN;
  gamXX=BMFheader.metric_XX;
  gamXY=BMFheader.metric_XY;
  gamYY=BMFheader.metric_YY;
  RA=BMFheader.RA;
  dec=BMFheader.dec;
  sprintf(ifo,BMFheader.det);

  /* do a quick check that the targeted area lies within the full mesh boundaries */
  if ((sma0<smaMIN)||(sma0>smaMAX)) {
    fprintf(stderr,"ERROR : targetted value of semi major axis outside range of input mesh\n");
    exit(1);
  }
  tperi0FLT=tperi0.gpsSeconds+(1e-9)*tperi0.gpsNanoSeconds;
  tperiMINFLT=tperisecMIN+(1e-9)*tperinsMIN;
  tperiMAXFLT=tperisecMAX+(1e-9)*tperinsMAX;
  /* printf("o is %lf min is %lf max is %lf\n",tperi0FLT,tperiMINFLT,tperiMAXFLT); */
  if ((tperi0FLT<tperiMINFLT)||(tperi0FLT>tperiMAXFLT)) {
    fprintf(stderr,"ERROR : targetted value of periapse passage time outside range of input mesh\n");
    exit(1);
  }
      
  /* get the ssb time of the start of the observation so also do setup baryinput */
  if (SetupBaryInput()) return 1;
  if (GetSSBTime(&TstartDET,&TstartSSB)) return 2;

  /* allocate some memory */
  sma=(REAL8 *)LALMalloc(Nfull*sizeof(REAL8));
  tperi=(LIGOTimeGPS *)LALMalloc(Nfull*sizeof(LIGOTimeGPS));
  dist=(REAL8 *)LALMalloc(Nfull*sizeof(REAL8));
  
  /* fill in some obvious info */
  RTPloc.ecc=0.0;
  RTPloc.argp=0.0;
  RTPloc.period=period;
  RTPloc.tstartSSB.gpsSeconds=TstartSSB.gpsSeconds;
  RTPloc.tstartSSB.gpsNanoSeconds=TstartSSB.gpsNanoSeconds;

  /* calculate X0 and Y0 at centre of p-space */
  RTPloc.sma=sma0;
  RTPloc.tperi.gpsSeconds=tperi0.gpsSeconds;
  RTPloc.tperi.gpsNanoSeconds=tperi0.gpsNanoSeconds;
  
  /* do the parameter conversion */
  if (ConvertRTperitoXY(&RTPloc,&XYloc,&dummy)) return 1; 
  X0=XYloc.X;
  Y0=XYloc.Y;

  /* cycle through the templates and calculate the distance from the central point */ 
  for (i=0;i<Nfull;i++) {
    fscanf(fbfp,"%lf%lf%d%d%lf%lf",&sma[i],&dummy,&tperi[i].gpsSeconds,&tperi[i].gpsNanoSeconds,&dummy,&dummy);
    RTPloc.sma=sma[i];
    RTPloc.tperi.gpsSeconds=tperi[i].gpsSeconds;
    RTPloc.tperi.gpsNanoSeconds=tperi[i].gpsNanoSeconds;
    if (ConvertRTperitoXY(&RTPloc,&XYloc,&dummy)) return 1; 
    if (CalculateDistance(&XYloc.X,&XYloc.Y,&dist[i])) return 2;
  
    if (dist[i]<dist_min) {
      dist_min = dist[i];
      Xmin=XYloc.X;
      Ymin=XYloc.Y;
    }
  }
  
  
  fclose(fbfp);

  return 0;

}

/*******************************************************************************/

int FreeMem(void)
{

  LALFree(dist);
  LALFree(sma);
  LALFree(tperi);

  return 0;

}

/*******************************************************************************/

int CalculateDistance(REAL8 *x,REAL8 *y,REAL8 *ds)
{

  /* this function simply calculates the ditance from the central location using the metric */

  REAL8 distSQ;

  distSQ=gamXX*(*x-X0)*(*x-X0)+gamYY*(*y-Y0)*(*y-Y0)+2.0*gamXY*(*x-X0)*(*y-Y0);
  *ds=sqrt(distSQ);

  return 0;

}

/*******************************************************************************/
/** Sorting function to sort into INCREASING order. Used in SortDistances(). */
int compare(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=dist[*(const int *)ip];
  dj=dist[*(const int *)jp];

  if (di>dj)
    return 1;
  
  if (di==dj)
    return 0;

  return -1;
}

/*******************************************************************************/

INT4 OutputSortedDist(void)
{

  INT4 *indexes,i,nlow;
   
  /*    create an array of indexes */
  if (!(indexes=(INT4 *)LALMalloc(sizeof(INT4)*Nfull))){
    fprintf(stderr,"Unable to allocate index array in CalculateDist()\n");
    return 1;
  }

  /*    populate it */
  for (i=0;i<Nfull;i++)
    indexes[i]=i;

  /*   sort array of indexes */
  qsort((void *)indexes, (size_t)Nfull, sizeof(int), compare);

  

  sbfp=fopen(subbankfile,"w");
  if (sbfp==NULL) {
    fprintf(stderr,"Unable to open file %s\n",subbankfile);
    return 1;
  }

  if (sbfp!=NULL) {

    /* setup the header input */
    BMFheader.f_max=f_max;
    BMFheader.tspan=tspan;
    BMFheader.tstart.gpsSeconds=TstartDET.gpsSeconds;
    BMFheader.tstart.gpsNanoSeconds=TstartDET.gpsNanoSeconds;
    BMFheader.Nfilters=Nsub;
    BMFheader.mismatch=mismatch;
    BMFheader.sma_0=sma0FULL;
    BMFheader.sma_MIN=smaMIN;
    BMFheader.sma_MAX=smaMAX;
    BMFheader.tperi_0.gpsSeconds=tperisec0FULL;
    BMFheader.tperi_0.gpsNanoSeconds=tperins0FULL;
    BMFheader.tperi_MIN.gpsSeconds=tperisecMIN;
    BMFheader.tperi_MIN.gpsNanoSeconds=tperinsMIN;
    BMFheader.tperi_MAX.gpsSeconds=tperisecMAX;
    BMFheader.tperi_MAX.gpsNanoSeconds=tperinsMAX;
    BMFheader.ecc_MIN=ecc;
    BMFheader.ecc_MAX=ecc;
    BMFheader.argp_MIN=argp;
    BMFheader.argp_MAX=argp;
    BMFheader.period_MIN=period;
    BMFheader.period_MAX=period;
    BMFheader.metric_XX=gamXX;
    BMFheader.metric_XY=gamXY;
    BMFheader.metric_YY=gamYY;
    sprintf(BMFheader.version,"v1_sub");
    sprintf(BMFheader.det,ifo);
    BMFheader.RA=RA;
    BMFheader.dec=dec;

  
    /* write header information to the sub mesh */
    if (WriteMeshFileHeader(sbfp,&BMFheader)) return 3;

    if (exactflag==0) {
      
      /*    print out the lowest Nsub ones */
      for (nlow=0; nlow<Nsub; nlow++) {
	
	fprintf(sbfp,"%6.12f %6.12f %d %d %6.12f %6.12f\n",
		sma[indexes[nlow]],period,tperi[indexes[nlow]].gpsSeconds,
		tperi[indexes[nlow]].gpsNanoSeconds,ecc,argp);
	
      }
    }
    else {
      fprintf(sbfp,"%6.12f %6.12f %d %d %6.12f %6.12f\n",
		sma0,period,tperi0.gpsSeconds,
		tperi0.gpsNanoSeconds,ecc,argp);

    }
    
    fclose(sbfp);
  }
  else {
    printf("problem outputting to file %s\n",subbankfile); 
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
  sma0=0.0;  /* a */
  tperi0.gpsSeconds=0; /* T */ 
  tperi0.gpsNanoSeconds=0; /* t */
  Nsub=1; /* N */
  exactflag=0;
 
  sprintf(fullbankfile," "); /* f */
  sprintf(subbankfile," "); /* s */
  sprintf(ephemfile," "); /* E */
  sprintf(yr,"00-04"); /* I */

  {
    int option_index = 0;
    static struct option long_options[] = {
                {"sma", required_argument, 0, 'r'},
		{"tpsec", required_argument, 0, 'T'},
		{"tpnano", required_argument, 0, 't'},
		{"fullbank", required_argument, 0, 'f'},
                {"subbank", required_argument, 0, 's'},
		{"ephem", required_argument, 0, 'E'},
                {"yr", required_argument, 0, 'y'},
		{"Nsub", required_argument, 0, 'N'},
		{"exact", no_argument, 0, 'x'},
		{"help", no_argument, 0, 'h'}
    };
  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt_long (argc, argv,"hr:T:t:f:s:E:y:N:x",long_options, &option_index)))!=-1)
    switch (c) {
    case 'r':
      sma0=atof(optarg);
      break;
    case 'T':
      tperi0.gpsSeconds=atoi(optarg);
      break;
    case 't':
      tperi0.gpsNanoSeconds=atoi(optarg);
      break;
    case 'f':
      temp=optarg;
      sprintf(fullbankfile,temp);
      break; 
    case 's':
      temp=optarg;
      sprintf(subbankfile,temp);
      break;
    case 'E':
      temp=optarg;
      sprintf(ephemfile,temp);
      break;
    case 'y':
      temp=optarg;
      sprintf(yr,temp);
      break;
    case 'N':
      Nsub=atoi(optarg);
      break;
    case 'x':
      exactflag=1;
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\t--sma\tFLOAT\t The central value of the projected semi-major axis in seconds [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--tpsec\tINT\t The central value of time of periapse passage in integer seconds [DEFAULT=0]\n");
      fprintf(stdout,"\t--tpnano\tINT\t The central value of time of periapse passage in integer nanoseconds [DEFAULT=0]\n");
      fprintf(stdout,"\t--fullbank\tSTRING\t The name of the full template bank file [DEFAULT=NULL]\n");
      fprintf(stdout,"\t--subbank\t\tSTRING\t The name of the sub template file [DEFAULT=NULL]\n");
      fprintf(stdout,"\t--ephem\tSTRING\t The location of the ephemeris files [DEFAULT=./]\n");
      fprintf(stdout,"\t--yr\t\tSTRING\t The year of the ephemeris file to use [DEFAULT=00-04]\n");
      fprintf(stdout,"\t--Nsub\t\tINT\t The number of filters required in the sub template file [DEFAULT=1]\n");
      fprintf(stdout,"\t--exact\t\tBOOLEAN\t Set this if you require a single template at the given coords [DEFAULT=0]\n");
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
 
  if (exactflag==1) Nsub=1;
 
  /* update global variable and return */
  return errflg;
}



/*******************************************************************************/

