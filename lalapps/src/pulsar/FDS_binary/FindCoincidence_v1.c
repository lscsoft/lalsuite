/*********************************************************************************/
/* Program to identify a subset of orbital templates in a secondary detector     */
/* given a single template location in the primary detector.                     */
/*                                                                               */
/*			           C. Messenger                                  */
/*                                                                               */
/*                         BIRMINGHAM UNIVERISTY -  2004                         */
/*********************************************************************************/

#include "FindCoincidence_v1.h"
#include "ReadSourceFile_v1.h"
#include "GenerateBinaryMesh_v1.h"

INT4 lalDebugLevel=3;
static LALStatus status;
REAL8 *dist;
LIGOTimeGPS tperi0,*tperi,TstartDET,TstartSSB,tperiCURRENT;
REAL8 f_max,tspan;
INT4 tperisecMIN,tperisecMAX,tperinsMIN,tperinsMAX,TstartsecMIN,TstartsecMAX;
REAL8 smaMIN,smaMAX,mismatch,sma0FULL;
REAL8 gamXX,gamXY,gamYY;
REAL8 smaCURRENT;
INT4 Nsub,Nfull;
INT4 tperisec0FULL,tperins0FULL;
CHAR ifo[256];
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

/* clargs */
char primarybank[256];
char secondarybank[256];
char presultsdir[256];
char sresultsdir[256];
REAL8 f_min, f_max;
REAL8 df;
char coresultsdir[256];
char ephdir[256];
char yr[256];

extern char *optarg;
extern int optind, opterr, optopt;

int ReadCommandLine(int argc,char *argv[]);
int ReadHeader(char *,BinaryMeshFileHeader *);
int FindCoincidence(BinaryMeshFileHeader *,BinaryMeshFileHeader *,REAL8,Result,Results *,CoResults **);
int OutputCoincidence(char *,REAL8,REAL8,CoResults *);
int FreeMem(void);
int SetupBaryInput(char *,char *,char *,LIGOTimeGPS *);
int GetSSBTime(LALDetector *,REAL8 *,REAL8 *,LIGOTimeGPS *, LIGOTimeGPS *);
int CalculateDistance(XYLocation *, XYLocation *, BinaryMeshFileHeader *, REAL8 *);
int CheckCoincidence(RTPLocation *,RTPLocation *,BinaryMeshFileHeader *,BinaryMeshFileHeader *);
int CheckConsistency(BinaryMeshFileHeader *, BinaryMeshFileHeader *);
int ReadResults(char *,REAL8,REAL8,REAL8 *,Results **);
int CalculateSignificance(REAL8,REAL8); 

int main(int argc,char *argv[]) 
{
 
  BinaryMeshFileHeader p_BMFheader;
  BinaryMeshFileHeader s_BMFheader;
  Results *p_results;
  Results *s_results;
  CoResults *co_results;
  REAL8 *dummy_df=NULL;
  INT4 i;

  /* read in the command linie arguments */
  if (ReadCommandLine(argc,argv)) return 1;

  /* read in the header information from the primary template bank */
  if (ReadHeader(primarybank,&p_BMFheader)) return 2;

  /* read in the header information from the secondary template bank */
  if (ReadHeader(secondarybank,&s_BMFheader)) return 2;

  /* read in primary results */
  if (ReadResults(presultsdir,f_min,f_max,dummy_df,&p_results)) return 3;

  /* read in secondary results */
  if (ReadResults(sresultsdir,f_min,f_max,&df,&s_results)) return 3;

  /* loop over each primary result */
  for (i=0;i<p_results->Nresults;i++) {

    /* find the coincident templates in the secondary bank */
    if (FindCoincidence(&p_BMFheader,&s_BMFheader,df,p_results->result[i],s_results,&co_results)) return 3;
    
    /* output the coincident results */
    if (OutputCoincidence(coresultsdir,f_min,f_max,co_results)) return 5;
    
    /* free up any allocated memory */
    if (FreeMem()) return 6;
    
  }
  

  return 0;

}

/****************************************************************************/

int ReadResults(char *resultsdir, REAL8 min_f, REAL8 max_f, REAL8 *deltaf, Results **results)
{

  /* this function reads in files from a given results directory.  It stores only */
  /* those results that lie within a given frequency window. */
  
  INT4 fileno=0;
  FILE *fp;
  char **filelist;
  char command[512];
  glob_t globbuf;
  INT4 i,k;
  INT4 nfiles;
  REAL8 f,ra,dec,sma,p,ecc,argp,mean,std,twoF;
  INT4 tpsec,tpnano,N;

  /* set the correct frequency boundaries */
  if (deltaf!=NULL) {
    min_f=min_f-(*deltaf);
    max_f=max_f+(*deltaf);
  }

  
  /* set up the datadir name */
  strcpy(command, resultsdir);
  strcat(command,"/Fstats_*");
    
  /* set up some glob stuff */
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);
  
  /* check if there are any SFT's in the directory */
  if(globbuf.gl_pathc==0)
    {
      fprintf (stderr,"\nNo SFTs in directory %s ... Exiting.\n", resultsdir);
      exit(1);
    }
  
  /* allocate memory for the pathnames */
  filelist=(char **)LALMalloc(globbuf.gl_pathc*sizeof(char *));
  for (i=0;i<(INT4)globbuf.gl_pathc;i++) filelist[i]=(char *)LALMalloc(256*sizeof(char));
  
  /* read all file names into memory */
  while ((UINT4)fileno < globbuf.gl_pathc) 
    {
      strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
      if (fileno > MAXFILES)
	{
	  fprintf(stderr,"\nToo many files in directory! Exiting... \n");
	  exit(1);
	}
    }
  globfree(&globbuf);

  nfiles=fileno;

  k=0;
  /* this first loop is just to couont how many results we need to store */
  /* open each file and read contents */
  for (i=0;i<nfiles;i++) {
    
    fp=fopen(filelist[i],"r");
    if (fp==NULL) {
      printf("ERROR : could not open file %s, which is strange\n",filelist[i]);
      exit(1);
    }

    /* read in each line */
    while (fscanf(fp,"%lf %lf %lf %lf %lf %d %d %lf %lf %d %lf %lf %lf\n",
		  &f,&ra,&dec,&sma,&p,&tpsec,&tpnano,&ecc,&argp,&N,&mean,&std,&twoF)!=EOF) {
  
      /* if in the correct frequency band then count it */
      if ((f>min_f)&&(f<max_f)) {
  	k++;
      }
      
    }
    
    /* close the current file */
    fclose(fp);
    
  }

  /* define the number of results for this directory in this band */
  N=k;

  /* allocate some memory */
  (*results)=(Results *)LALMalloc(sizeof(Results));
  (*results)->result=(Result *)LALMalloc(N*sizeof(Result)); 

  k=0;
  /* open each file and read contents and actually store them */
  for (i=0;i<nfiles;i++) {
    
    fp=fopen(filelist[i],"r");
    if (fp==NULL) {
      printf("ERROR : could not open file %s, which is strange\n",filelist[i]);
      exit(1);
    }

    /* read in each line */
    while (fscanf(fp,"%lf %lf %lf %lf %lf %d %d %lf %lf %d %lf %lf %lf\n",
		  &f,&ra,&dec,&sma,&p,&tpsec,&tpnano,&ecc,&argp,&N,&mean,&std,&twoF)!=EOF) {
      
      /* if in the correct frequency band then store the result */
      if ((f>min_f)&&(f<max_f)) {
	(*results)->result[k].freq=f;
	(*results)->result[k].RA=ra;
	(*results)->result[k].dec=dec;
	(*results)->result[k].sma=sma;
	(*results)->result[k].period=p;
	(*results)->result[k].tp.gpsSeconds=tpsec;
	(*results)->result[k].tp.gpsNanoSeconds=tpnano;
	(*results)->result[k].ecc=ecc;
	(*results)->result[k].argp=argp;
	(*results)->result[k].twoF=twoF;
	k++;
      }
      
    }

    /* close the current file */
    fclose(fp);

  }

  return 0;

}

/****************************************************************************/

int CheckConsistency(BinaryMeshFileHeader *p_BMFheader, BinaryMeshFileHeader *s_BMFheader)
{

  /* this function compares the header information from both template files */
  /* and checks for consistency */

  /* check for consistency in maximum frequency */
  if (p_BMFheader->f_max!=s_BMFheader->f_max) {
    fprintf(stdout,"ERROR : fmax not equal in both detectors !!\n");
    exit(1);
  }
   
  /* here we check for a zero-range in orbital period for the primary detector */
  if (p_BMFheader->period_MIN!=p_BMFheader->period_MAX) {
    fprintf(stdout,"ERROR : Orbital period range is non-zero in primary detector !!\n");
    exit(1);
  }
  
  /* here we check for a zero-range in orbital period for the secondary detector */
  if (s_BMFheader->period_MIN!=s_BMFheader->period_MAX) {
    fprintf(stdout,"ERROR : Orbital period range is non-zero in secondary detector !!\n");
    exit(1);
  }
  
  /* here we check for consistent orbital periods between detectors  */
  if (p_BMFheader->period_MIN!=s_BMFheader->period_MIN) {
    fprintf(stdout,"ERROR : Primary and Secondary detectors have different orbital periods !!\n");
    exit(1);
  }

  /* here we chack for identical sky positions */
  if (p_BMFheader->RA!=s_BMFheader->RA) {
    fprintf(stdout,"ERROR : RA not equal in both detectors !!\n");
    exit(1);
  }
  if (p_BMFheader->dec!=s_BMFheader->dec) {
    fprintf(stdout,"ERROR : dec not equal in both detectors !!\n");
    exit(1);
  }
  
     /* could add more checks in the future */

  return 0;

}

/****************************************************************************/

int GetSSBTime(LALDetector *Det, REAL8 *skyRA, REAL8 *skydec, LIGOTimeGPS *tdet, LIGOTimeGPS *tssb)
{

  /* this function returns the SSB time given a GPS time, a detector, and a sky location */

  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */

 /* setup the input for the barycentering  */
  baryinput.tgps.gpsSeconds=tdet->gpsSeconds;
  baryinput.tgps.gpsNanoSeconds=tdet->gpsNanoSeconds;
  baryinput.site.location[0]=Det->location[0]/LAL_C_SI;
  baryinput.site.location[1]=Det->location[1]/LAL_C_SI;
  baryinput.site.location[2]=Det->location[2]/LAL_C_SI;
  baryinput.alpha=*skyRA;
  baryinput.delta=*skydec;
  baryinput.dInv=0.e0;

  /* Call the barycentering routines */
  LALBarycenterEarth(&status,&earth,tdet,edat);
  LALBarycenter(&status,&emit,&baryinput,&earth);
  tssb->gpsSeconds=(emit.te.gpsSeconds);
  tssb->gpsNanoSeconds=emit.te.gpsNanoSeconds;
  
  return 0;
}

/*******************************************************************************/

int SetupBaryInput(char *ephfile, char *year, char *detector, LIGOTimeGPS *tstart)
{

  /* this function sets up all the required variables and structures for */
  /* doing barycentering */

  CHAR filenameE[256],filenameS[256];
  FILE *fp;
  
  /* make the full file name/location for the ephemeris files */
  strcpy(filenameE,ephfile);
  strcat(filenameE,"/earth");
  strcat(filenameE,year);
  strcat(filenameE,".dat");
  strcpy(filenameS,ephfile);
  strcat(filenameS,"/sun");
  strcat(filenameS,year);
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
  LALLeapSecs(&status,&leap,tstart,&formatAndAcc);
  (*edat).leap=leap;

  /* Read in ephemeris files */
  LALInitBarycenter(&status,edat);             

  /* setup chosen detector */
  if(strcmp(detector,"GEO")!=0) Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if(strcmp(detector,"LLO")!=0) Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(strcmp(detector,"LHO")!=0) Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
 
  return 0;

}
            
/*******************************************************************************/

int ReadHeader(char *meshfilename, BinaryMeshFileHeader *BMF) 
{

  /* here we read template file header information */

  FILE *fp;
  
  /* open the full bank file */
  fp=fopen(meshfilename,"r");
  if (fp==NULL) {
    fprintf(stderr,"Unable to open file %s\n",meshfilename);
    return 1;
  }
  
  /* call the routine for reading the header into memory */
  if (ReadMeshFileHeader(fp,BMF)) return 1;
    
  return 0;
  
}

/*******************************************************************************/

int FindCoincidence(BinaryMeshFileHeader *p_BMFheader,BinaryMeshFileHeader *s_BMFheader,REAL8 deltaf, Result p_result, Results *s_results, CoResults **co_results)
{

  /* this function cycles over every secondary event and compares it to the current */
  /* primary event to check for coinicidence */

  RTPLocation p_RTPloc;
  RTPLocation s_RTPloc;
  REAL8 min_f;
  REAL8 max_f;
  INT4 k,i;
  REAL8 f_s;
  REAL8 p_sig=0.0;
  REAL8 s_sig=0.0;
  REAL8 co_sig=0.0;
  INT4 MaxCo;
  

   /* now lets check if the primary and secondary headers have consistent information */
  if (CheckConsistency(p_BMFheader,s_BMFheader)) return 2;

  /* get the ssb time of the start of the primary observation so also do setup baryinput */
  if (SetupBaryInput(ephdir,yr,p_BMFheader->det,&(p_BMFheader->tstart))) return 1;
  if (GetSSBTime(&Detector,&(p_BMFheader->RA),&(p_BMFheader->dec),&(p_BMFheader->tstart),&p_RTPloc.tstartSSB)) return 2;
  
  /* get the ssb time of the start of the secondary observation so also do setup baryinput */
  if (SetupBaryInput(ephdir,yr,s_BMFheader->det,&(s_BMFheader->tstart))) return 3;
  if (GetSSBTime(&Detector,&(s_BMFheader->RA),&(s_BMFheader->dec),&(s_BMFheader->tstart),&s_RTPloc.tstartSSB)) return 4;

  /*printf("primary ssbtime is %d %d\n",p_RTPloc.tstartSSB.gpsSeconds,p_RTPloc.tstartSSB.gpsNanoSeconds);
  printf("secondary ssbtime is %d %d\n",s_RTPloc.tstartSSB.gpsSeconds,s_RTPloc.tstartSSB.gpsNanoSeconds);*/

   /* fill in info for central point definied by primary template location */
  p_RTPloc.sma=p_result.sma;
  p_RTPloc.tperi.gpsSeconds=p_result.tp.gpsSeconds;
  p_RTPloc.tperi.gpsNanoSeconds=p_result.tp.gpsNanoSeconds;
  p_RTPloc.ecc=p_result.ecc;
  p_RTPloc.argp=p_result.argp;
  p_RTPloc.period=p_result.period;

  /* set the frequency range */
  min_f=p_result.freq-deltaf;
  max_f=p_result.freq+deltaf;

  k=0;
  /* pre-allocate memory to store coincident locations */
  MaxCo=10000;
  (*co_results)=(CoResults *)LALMalloc(sizeof(CoResults));
  (*co_results)->primary_result=(Result *)LALMalloc(MaxCo*sizeof(Result));
  (*co_results)->secondary_result=(Result *)LALMalloc(MaxCo*sizeof(Result));
  (*co_results)->significance=(Significance *)LALMalloc(MaxCo*sizeof(Significance));

  /* cycle through the secondary events */ 
  for (i=0;i<s_results->Nresults;i++) {

    /* check if the secondary frequency is in the range */
    f_s=s_results->result[i].freq;
    if ((f_s>min_f)&&(f_s<max_f)) {
          
      /* store the current secondary parameters in an RTPloc structure */
      s_RTPloc.sma=s_results->result[i].sma;
      s_RTPloc.tperi.gpsSeconds=s_results->result[i].tp.gpsSeconds;
      s_RTPloc.tperi.gpsNanoSeconds=s_results->result[i].tp.gpsNanoSeconds;
      s_RTPloc.ecc=s_results->result[i].ecc;
      s_RTPloc.argp=s_results->result[i].argp;
      s_RTPloc.period=s_results->result[i].period;

      /* call a routine that checks whether this location is consistent with coincidence */
      if (CheckCoincidence(&p_RTPloc,&s_RTPloc,p_BMFheader,s_BMFheader)) {

	/* calculate the significances of the seperate results */
	if (CalculateSignificance(p_result.twoF,p_sig)) return 3;
	if (CalculateSignificance(s_results->result[i].twoF,s_sig)) return 4;
	co_sig=p_sig+s_sig-(p_sig*s_sig);

	/* if it is consistent then store it */
	(*co_results)->primary_result[k].freq=p_result.freq;
	(*co_results)->primary_result[k].RA=p_result.RA;
	(*co_results)->primary_result[k].dec=p_result.dec;
	(*co_results)->primary_result[k].sma=p_RTPloc.sma;
	(*co_results)->primary_result[k].period=p_RTPloc.period;
	(*co_results)->primary_result[k].tp.gpsSeconds=p_RTPloc.tperi.gpsSeconds;
	(*co_results)->primary_result[k].tp.gpsNanoSeconds=p_RTPloc.tperi.gpsNanoSeconds;
	(*co_results)->primary_result[k].ecc=p_RTPloc.ecc;
	(*co_results)->primary_result[k].argp=p_RTPloc.argp;
	(*co_results)->primary_result[k].twoF=p_result.twoF;
	(*co_results)->significance[k].log10oneminusp_sig=log10(p_sig);
	(*co_results)->secondary_result[k].freq=s_results->result[i].freq;
	(*co_results)->secondary_result[k].RA=s_results->result[i].RA;
	(*co_results)->secondary_result[k].dec=s_results->result[i].dec;
	(*co_results)->secondary_result[k].sma=s_RTPloc.sma;
	(*co_results)->secondary_result[k].period=s_RTPloc.period;
	(*co_results)->secondary_result[k].tp.gpsSeconds=s_RTPloc.tperi.gpsSeconds;
	(*co_results)->secondary_result[k].tp.gpsNanoSeconds=s_RTPloc.tperi.gpsNanoSeconds;
	(*co_results)->secondary_result[k].ecc=s_RTPloc.ecc;
	(*co_results)->secondary_result[k].argp=s_RTPloc.argp;
	(*co_results)->secondary_result[k].freq=s_results->result[i].twoF;
	(*co_results)->significance[k].log10oneminuss_sig=log10(s_sig);
	(*co_results)->significance[k].log10oneminusco_sig=log10(co_sig);
	k++;
      }
    }
  }

  /* this is the number of coincident templates found */
  (*co_results)->Nresults=k;

  return 0;

}

/*******************************************************************************/

int FreeMem()
{

  /* this routine simply frees up memory at the end */

  LALFree(edat);

  return 0;

}


/*******************************************************************************/

int CheckCoincidence(RTPLocation *p_RTPloc,RTPLocation *s_RTPloc,BinaryMeshFileHeader *p_BMFheader,BinaryMeshFileHeader *s_BMFheader)
{

  /* This function calclates whether a location in the secondary detector is a coincidence event */ 
  /* it uses a different method based on intersecting boxes */

  REAL8Array *eigvec=NULL;
  REAL8Vector *eigval=NULL;
  REAL8 p_a,p_b;
  REAL8 s_a,s_b;
  RTPLocation temp_RTPloc;
  XYLocation pins_XYloc;
  XYLocation sins_XYloc;
  REAL8 dummy;
  UINT4Vector *dimlength=NULL;
  INT4 dim=2;
  INT4 i,j,k;
  Corner *p_corner;
  Corner *s_corner;
  INT4 *line;
  REAL8 xp_a,xp_b;
  REAL8 yp_a,yp_b;
  REAL8 xs_a,xs_b;
  REAL8 ys_a,ys_b;
  REAL8 p_grad,s_grad;
  REAL8 p_int,s_int;
  REAL8 p_int_xmin,p_int_xmax;
  REAL8 p_int_ymin,p_int_ymax;
  REAL8 s_int_xmin,s_int_xmax;
  REAL8 s_int_ymin,s_int_ymax;
  INT4 p_true=0;
  INT4 s_true=0;
  REAL8 x_cross;
  REAL8 y_cross;

  /* Lets convert primary detector central location into XY coords in secondary detector */
  /* This means finding the equivelent XY location in the secondary detector, ie. if a signal */
  /* were present at this location in the primary, then this is where it would be in the */
  /* secondary detector */

  /* so we use primary sma and tperi location */
  temp_RTPloc.sma=p_RTPloc->sma;
  temp_RTPloc.tperi.gpsSeconds=p_RTPloc->tperi.gpsSeconds;
  temp_RTPloc.tperi.gpsNanoSeconds=p_RTPloc->tperi.gpsNanoSeconds;
  temp_RTPloc.ecc=0.0;
  temp_RTPloc.argp=0.0;
  temp_RTPloc.period=p_RTPloc->period;

  /* but use secondary detector SSB start time (this is the parameter that causes each detectors */
  /* XY space covering to NOT neccesserily overlap) */ 
  temp_RTPloc.tstartSSB.gpsSeconds=s_RTPloc->tstartSSB.gpsSeconds;
  temp_RTPloc.tstartSSB.gpsNanoSeconds=s_RTPloc->tstartSSB.gpsNanoSeconds;
  
  /* Now we take this RT location and find the XY location in the secondary detector */
  if (ConvertRTperitoXY(&temp_RTPloc,&pins_XYloc,&dummy)) return 1;

  /* Now convert current secondary detector location into XY in secondary detector */
  /* Note that now we are considering points in XY parameter space based on the */
  /* secondary detector search */
  if (ConvertRTperitoXY(s_RTPloc,&sins_XYloc,&dummy)) return 1;
  
  /* Also remember that in XY space the metric is constant and so eigenvectors */
  /* and eigenvalues remain constant for each detector.  Therefore boxes defining */
  /* closest neighbour regions are invariant under translation ie. the primary */
  /* detector eigenvectors/values do not change now that we are dealing with the */
  /* secondary detector */
  
  /* Now need to find the angle that the first eigenvector associated with the primary */
  /* metric makes with the X-axis.  So need to diagonalise the gamma matrix  */
  /* here we find the eigenvalues and eigenvectors */
  LALU4CreateVector(&status,&dimlength,dim);
  dimlength->data[0]=dim;
  dimlength->data[1]=dim;
  LALDCreateArray(&status,&eigvec,dimlength);
  LALDCreateVector(&status,&eigval,DIM);
  LALU4DestroyVector(&status,&dimlength);
  eigvec->data[0]=p_BMFheader->metric_XX;
  eigvec->data[1]=p_BMFheader->metric_XY;
  eigvec->data[2]=p_BMFheader->metric_XY;
  eigvec->data[3]=p_BMFheader->metric_YY;
  LALDSymmetricEigenVectors(&status,eigval,eigvec);

  /* find the length of the eigenvectors defining the size of the box */
  p_a=sqrt(p_BMFheader->mismatch/(2.0*eigval->data[0]));
  p_b=sqrt(p_BMFheader->mismatch/(2.0*eigval->data[1]));

  /* allocate memory for corners of box structures */
  p_corner=(Corner *)LALMalloc(4*sizeof(Corner));
  s_corner=(Corner *)LALMalloc(4*sizeof(Corner));  
 
  /* identify the corners of the primary box */
  p_corner[0].x=pins_XYloc.X+(p_a*eigvec->data[0]+p_b*eigvec->data[1]);
  p_corner[0].y=pins_XYloc.Y+(p_a*eigvec->data[2]+p_b*eigvec->data[3]);
  p_corner[1].x=pins_XYloc.X+(p_a*eigvec->data[0]-p_b*eigvec->data[1]);
  p_corner[1].y=pins_XYloc.Y+(p_a*eigvec->data[2]-p_b*eigvec->data[3]);
  p_corner[2].x=pins_XYloc.X+(-p_a*eigvec->data[0]-p_b*eigvec->data[1]);
  p_corner[2].y=pins_XYloc.Y+(-p_a*eigvec->data[2]-p_b*eigvec->data[3]);
  p_corner[3].x=pins_XYloc.X+(-p_a*eigvec->data[0]+p_b*eigvec->data[1]);
  p_corner[3].y=pins_XYloc.Y+(-p_a*eigvec->data[2]+p_b*eigvec->data[3]);

  /* Now we deal with the eigenvectors for the second detector */
  /* first lets diagonalise the second template banks matrix */
  eigvec->data[0]=s_BMFheader->metric_XX;
  eigvec->data[1]=s_BMFheader->metric_XY;
  eigvec->data[2]=s_BMFheader->metric_XY;
  eigvec->data[3]=s_BMFheader->metric_YY;
  LALDSymmetricEigenVectors(&status,eigval,eigvec);

  /* find the dimensions of the secondary box */
  s_a=sqrt(s_BMFheader->mismatch/(2.0*eigval->data[0]));
  s_b=sqrt(s_BMFheader->mismatch/(2.0*eigval->data[1]));
  
  /* now we position the secondary box centered at each corner of the primary box */ 
  /* for each corner position the secondary box centered on that corner and record */
  /* each of the 4 locations of the secondary box corners */
  /* loop over the primary corners */
  
  s_corner[0].x=sins_XYloc.X+(s_a*eigvec->data[0]+s_b*eigvec->data[1]);
  s_corner[0].y=sins_XYloc.Y+(s_a*eigvec->data[2]+s_b*eigvec->data[3]);
  s_corner[1].x=sins_XYloc.X+(s_a*eigvec->data[0]-s_b*eigvec->data[1]);
  s_corner[1].y=sins_XYloc.Y+(s_a*eigvec->data[2]-s_b*eigvec->data[3]);
  s_corner[2].x=sins_XYloc.X+(-s_a*eigvec->data[0]-s_b*eigvec->data[1]);
  s_corner[2].y=sins_XYloc.Y+(-s_a*eigvec->data[2]-s_b*eigvec->data[3]);
  s_corner[3].x=sins_XYloc.X+(-s_a*eigvec->data[0]+s_b*eigvec->data[1]);
  s_corner[3].y=sins_XYloc.Y+(-s_a*eigvec->data[2]+s_b*eigvec->data[3]);
  
  /* identify the indices for pairs of corners */
  line=(INT4 *)LALMalloc(4*sizeof(INT4));
  line[0]=1;
  line[1]=2;
  line[2]=3;
  line[3]=0;
 
  k=0;
  /* now we loop over each primary box corner */
  for (i=0;i<4;i++) {

    /* identify current line end points */
    xp_a=p_corner[i].x;
    yp_a=p_corner[i].y;
    xp_b=p_corner[line[i]].x;
    yp_b=p_corner[line[i]].y;

    /* calculate the current primary line equation */
    p_grad=(yp_a-yp_b)/(xp_a-xp_b);
    p_int=yp_a-(p_grad*xp_a);

    /* calculate possible intersection region */
    p_int_xmin=xp_a;
    p_int_xmax=xp_b;
    if(p_int_xmin>p_int_xmax) {
      p_int_xmin=xp_b;
      p_int_xmax=xp_a;
    }
    p_int_ymin=yp_a;
    p_int_ymax=yp_b;
    if(p_int_ymin>p_int_ymax) {
      p_int_ymin=yp_b;
      p_int_ymax=yp_a;
    }
    

   

    /* now we loop over each of the 4 lines joing the corners of the current secondary box */
    for (j=1;j<4;j++) {

      /* identify current line end points */
      xs_a=s_corner[j].x;
      ys_a=s_corner[j].y;
      xs_b=s_corner[line[j]].x;
      ys_b=s_corner[line[j]].y;

      /* calculate possible intersection region */
      s_int_xmin=xs_a;
      s_int_xmax=xs_b;
      if(s_int_xmin>s_int_xmax) {
	s_int_xmin=xs_b;
	s_int_xmax=xs_a;
      }
      s_int_ymin=ys_a;
      s_int_ymax=ys_b;
      if(s_int_ymin>s_int_ymax) {
	s_int_ymin=ys_b;
	s_int_ymax=ys_a;
      }

      /* calculate the current primary line equation */
      s_grad=(ys_a-ys_b)/(xs_a-xs_b);
      s_int=ys_a-(s_grad*xs_a);
      
          /* find where the lines intersect */
      x_cross=(p_int-s_int)/(s_grad-p_grad);
      y_cross=(p_grad*x_cross)+p_int;

      p_true=0;
      s_true=0;

      /* if they intesect in the correct region then it is a coincident pair */
      if ((x_cross>=p_int_xmin)&&(x_cross<=p_int_xmax)&&(y_cross>=p_int_ymin)&&(y_cross<=p_int_ymax)) p_true=1;
      if ((x_cross>=s_int_xmin)&&(x_cross<=s_int_xmax)&&(y_cross>=s_int_ymin)&&(y_cross<=s_int_ymax)) s_true=1;
      if ((p_true==1)&&(s_true==1)) {
	printf("grad int p %f %f\n",p_grad,p_int);
	printf("grad int p %f %f\n",s_grad,s_int);
	printf("p min-max %f %f\n",p_int_xmin,p_int_xmax);
	printf("p min-max %f %f\n",p_int_ymin,p_int_ymax);
	printf("s min-max %f %f\n",s_int_xmin,s_int_xmax);
	printf("s min-max %f %f\n",s_int_ymin,s_int_ymax);
	printf("%f %f %f %f %f %f %f %f %f %f\n",xp_a,yp_a,xp_b,yp_b,xs_a,ys_a,xs_b,ys_b,x_cross,y_cross);
      
	return 1;
	
      }
  
    }
  }
   
  LALFree(line);
  /* not a coincident pair */
  return 0;
  

  

}

/*******************************************************************************/

int CalculateDistance(XYLocation *XYlocONE, XYLocation *XYlocTWO, BinaryMeshFileHeader *BMF, REAL8 *ds)
{

  /* this function simply calculates the distance from the central location using the metric */

  REAL8 distSQ;
  REAL8 XX,YY,XY;
  REAL8 X0,Y0;
  REAL8 X,Y;

  /* define some shorter variable names */
  XX=BMF->metric_XX;
  XY=BMF->metric_XY;
  YY=BMF->metric_YY;
  X0=XYlocONE->X;
  Y0=XYlocONE->Y;
  X=XYlocTWO->X;
  Y=XYlocTWO->Y;

  
  /* calculate distance in the space defined by the given metric */
  distSQ=XX*(X-X0)*(X-X0)+YY*(Y-Y0)*(Y-Y0)+2.0*XY*(X-X0)*(Y-Y0);
  
  *ds=sqrt(distSQ);

  return 0;

}

/*******************************************************************************/

INT4 OutputCoincidence(char *outdir,REAL8 min_f,REAL8 max_f,CoResults *co_results)
{

  /* this function outputs the coincidence events to a file in the given directory */
  
  INT4 k;
  FILE *fp;
  char temp[256];
  char outfile[256];


  /* define output file name */
  sprintf(temp,"%s/coincidence_%f-%f.data",outdir,min_f,max_f);
  strcat(outfile,temp);

  /* open output file */
  fp=fopen(outfile,"w");
  if (fp==NULL) {
    fprintf(stderr,"Unable to open file %s\n",outfile);
    return 1;
  }

  /* loop over all the coincidence events */
  for (k=0;k<(*co_results).Nresults;k++) {
    
    fprintf(fp,"%f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f\n",
	    co_results->primary_result[k].freq,
	    co_results->primary_result[k].RA,
	    co_results->primary_result[k].dec,
	    co_results->primary_result[k].sma,
	    co_results->primary_result[k].period,
	    co_results->primary_result[k].tp.gpsSeconds,
	    co_results->primary_result[k].tp.gpsNanoSeconds,
	    co_results->primary_result[k].ecc,
	    co_results->primary_result[k].argp,
	    co_results->primary_result[k].twoF,
	    co_results->significance[k].log10oneminusp_sig,
	    co_results->secondary_result[k].freq,
	    co_results->secondary_result[k].RA,
	    co_results->secondary_result[k].dec,
	    co_results->secondary_result[k].sma,
	    co_results->secondary_result[k].period,
	    co_results->secondary_result[k].tp.gpsSeconds,
	    co_results->secondary_result[k].tp.gpsNanoSeconds,
	    co_results->secondary_result[k].ecc,
	    co_results->secondary_result[k].argp,
	    co_results->secondary_result[k].twoF,
	    co_results->significance[k].log10oneminuss_sig,
	    co_results->significance[k].log10oneminusco_sig);
      
  }
  
  /* close the output file */
  fclose(fp);
  
  return 0;
}

/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[]) 
{
  INT4 c, errflg = 0;
  char *temp;
  optarg = NULL;
  
  /* Initialize default values */    
  sprintf(primarybank," "); /* p */
  sprintf(secondarybank," "); /* s */
  sprintf(presultsdir," ");
  sprintf(sresultsdir," ");
  f_min=0.0;
  f_max=0.0;
  df=0.0;
  sprintf(coresultsdir," "); /* c */
  sprintf(ephdir," "); /* E */
  sprintf(yr,"00-04"); /* I */

  {
    int option_index = 0;
    static struct option long_options[] = {
                {"fmin", required_argument, 0, 'f'},
		{"fmax", required_argument, 0, 'F'},
		{"df", required_argument, 0, 'd'},
		{"pbank", required_argument, 0, 'b'},
		{"sbank", required_argument, 0, 'B'},
                {"presultsdir", required_argument, 0, 'r'},
		{"sresultsdir", required_argument, 0, 'R'},
                {"coresultsdir", required_argument, 0, 'o'},
		{"ephdir", required_argument, 0, 'E'},
		{"yr", required_argument, 0, 'y'},
		{"help", no_argument, 0, 'h'}
    };
  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt_long (argc, argv,"hf:F:d:b:B:r:R:o:E:y:",long_options, &option_index)))!=-1)
    switch (c) {
    case 'f':
      f_min=atof(optarg);
      break;
    case 'F':
      f_max=atof(optarg);
      break;
    case 'd':
      df=atof(optarg);
      break;
    case 'b':
      temp=optarg;
      sprintf(primarybank,temp);
      break; 
    case 'B':
      temp=optarg;
      sprintf(secondarybank,temp);
      break;
    case 'r':
      temp=optarg;
      sprintf(presultsdir,temp);
      break;
    case 'R':
      temp=optarg;
      sprintf(sresultsdir,temp);
      break;
    case 'o':
      temp=optarg;
      sprintf(coresultsdir,temp);
      break;
    case 'E':
      temp=optarg;
      sprintf(ephdir,temp);
      break;
    case 'y':
      temp=optarg;
      sprintf(yr,temp);
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\t--fmin          REAL8\t The minimum primary event frequency [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--fmax          REAL8\t The maximum primary event frequency [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--df            REAL8\t The maximum difference in event frequency between detectors [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--pbank         STRING\t The name of the primary template bank file [DEFAULT=]\n");
      fprintf(stdout,"\t--sbank         STRING\t The name of the secondary template bank file [DEFAULT=]\n");
      fprintf(stdout,"\t--presultsdir   STRING\t The location of the primary Fstat results directory [DEFAULT=]\n");
      fprintf(stdout,"\t--sresultsdir   STRING\t The location of the secondary Fstat results directory [DEFAULT=]\n");
      fprintf(stdout,"\t--coresultsdir  STRING\t The location of the coincidence results directory [DEFAULT=]\n");
      fprintf(stdout,"\t--ephdir        STRING\t The location of the ephemeris files [DEFAULT=]\n");
      fprintf(stdout,"\t--yr            STRING\t The year of the ephemeris file to use [DEFAULT=00-04]\n");
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

int CalculateSignificance(REAL8 twoF,REAL8 sig) 
{

  /* this function calculates the theoretical significance of a twoF value */
 

  sig=(1.0+(twoF/2.0))*exp((-1.0)*twoF/2.0);

  return 0;

}
