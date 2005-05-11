/************************************************************************************/
/*         Program to generate a grid of filters in parameter space for use         */
/*         in a search for signals from sources in circular binary orbits.          */
/*         A metric is calculated, projected and diagonalised.  Filters are         */
/*         then placed within the specified boundaries of the space and output      */
/*         to file ready for a search.                                              */
/*                                                                                  */
/*			           C. Messenger                                     */
/*                                                                                  */
/*                         BIRMINGHAM UNIVERISTY -  2005                            */
/************************************************************************************/

#include "GenerateBinaryMesh_v1.h"
#include "ReadSourceFile_v1.h"
/*#include "PeriapseShift_v1.h"*/
#include "ComputeFStatisticBinary_v2.h"

LIGOTimeGPS tperi_0_new;
LIGOTimeGPS tperi_MIN_new;
LIGOTimeGPS tperi_MAX_new;
LIGOTimeGPS tperi_GLOBAL;
INT4 NORB;
static LALStatus status;

int FreeMem(EphemerisData *);
int GenerateMesh(GlobVar,REAL4VectorSequence **,XYparameterspace *,Metric *);
int SetupPspaceParams(GlobVar,RTparameterspace *,XYparameterspace *);
int GenMetricComp(GlobVar,REAL8,REAL8 *,REAL8 *,Metric *);
int CheckRTBoundary(REAL8 *,LIGOTimeGPS *,RTparameterspace *);
int ConvertMesh(GlobVar,REAL4VectorSequence **,RTMesh *,RTparameterspace *);
int OutputRTMesh(GlobVar *,REAL8,RTMesh *,Metric *);
int ReadDataParams(char *, GlobVar *);
int SetGlobalVariables(CLargs,GlobVar *,binarysource);
int ReadCommandLine(int argc,char *argv[],CLargs *CLA);
int ConvertTperitoPhase(void);
int SetupBaryInput(GlobVar *,EphemerisData **,LALDetector *);
int CheckInput(GlobVar);
int GetSSBTime(GlobVar,LIGOTimeGPS *, LIGOTimeGPS *,EphemerisData **,LALDetector);

int main(int argc, char **argv){
  
  CLargs CLA;
  GlobVar GV;
  XYparameterspace XYspace;
  RTparameterspace RTspace;
  Metric XYMetric;
  RTMesh RTmesh;
  REAL4VectorSequence *XYmesh=NULL;
  binarysource sourceparams;
  EphemerisData *edat=NULL;          /* Stores earth/sun ephemeris data for barycentering */
  LALDetector Detector;              /* Our detector*/
  INT4 j;

  /* read the command line arguments */
  if (ReadCommandLine(argc,argv,&CLA)) return 1;

  /* if we have input a data dir (note that datadir dominates over clargs) */
  if (CLA.datadirflag) {
   
    /* look at the data directory and retrieve some parameters */
       if (ReadDataParams(CLA.datadir,&GV)) return 2;
  
  }
  /* else we use clargs for tstart and tspan */
  else {
    GV.tstart.gpsSeconds=CLA.tstart.gpsSeconds;
    GV.tstart.gpsNanoSeconds=CLA.tstart.gpsNanoSeconds;
    GV.tspan=CLA.tspan;
  }

  /* read the source parameters from the source file */
  if (ReadSource(CLA.sourcefile,CLA.source,&GV.tstart,&sourceparams)) return 2;

  /*printf("read sourcefile\n");*/

  /* Fill in the remaining global variables */
  if (SetGlobalVariables(CLA,&GV,sourceparams)) return 2;

  /*printf("set globvar\n");*/

  /* set up the required things for barycentering */
  if (SetupBaryInput(&GV,&edat,&Detector)) return 2;

  /*printf("setup baryinput\n");*/

  /* First job is to convert the observation start time from detector time to SSB time */
  if (GetSSBTime(GV,&GV.tstart,&GV.tstartSSB,&edat,Detector)) return 1;

  /*printf("got SSB time as %d %d\n",GV.tstartSSB.gpsSeconds,GV.tstartSSB.gpsNanoSeconds);*/
  

  /* setup the parameter space */
  if (SetupPspaceParams(GV,&RTspace,&XYspace)) return 2;
  
  /*printf("setup pspace\n");*/

  /* check the validity of the input */ 
  if (CheckInput(GV)) return 2;

  /*printf("checked validity\n");*/

  /* loop over the number of template banks we are going to generate */
  for (j=0;j<GV.nband;j++) {

    /* calculate the metric components for this max frequency */
    if (GenMetricComp(GV,GV.f_max[j],&(XYspace.X_0),&(XYspace.Y_0),&XYMetric)) return 3;
      
    /*printf("set gen metric comp\n");*/

    /* generate the mesh in XY space */
    if (GenerateMesh(GV,&XYmesh,&XYspace,&XYMetric)) return 4;
    
    /*printf("made XY mesh\n");*/

    /* convert this mesh to RT space */
    if (ConvertMesh(GV,&XYmesh,&RTmesh,&RTspace)) return 5;
   
    /*printf("made RT mesh\n");*/

    /* output the mesh to file */
    if (OutputRTMesh(&GV,GV.f_max[j],&RTmesh,&XYMetric)) return 6;
    
    /*printf("output mesh\n");*/

  }
  
  if (FreeMem(edat)) return 7;
  
  LALCheckMemoryLeaks();

  exit(0);

}

/***********************************************************************************/

int SetGlobalVariables(CLargs CLA,GlobVar *GV,binarysource sourceparams)
{

  /* here we just put the source params and dataparams into the global varaibles */

  INT4 j;
  REAL8 f;

  GV->band=CLA.band;
  GV->RA=sourceparams.skypos.ra;
  GV->dec=sourceparams.skypos.dec;
  GV->period=sourceparams.orbit.period;
  GV->sma_0=sourceparams.orbit.sma;
  GV->tperi_0.gpsSeconds=sourceparams.orbit.tperi.gpsSeconds;
  GV->tperi_0.gpsNanoSeconds=sourceparams.orbit.tperi.gpsNanoSeconds;
  GV->sma_MIN=sourceparams.orbit.sma_min;
  GV->sma_MAX=sourceparams.orbit.sma_max;
  GV->tperi_MIN.gpsSeconds=sourceparams.orbit.tperi_min.gpsSeconds;
  GV->tperi_MIN.gpsNanoSeconds=sourceparams.orbit.tperi_min.gpsNanoSeconds;
  GV->tperi_MAX.gpsSeconds=sourceparams.orbit.tperi_max.gpsSeconds;
  GV->tperi_MAX.gpsNanoSeconds=sourceparams.orbit.tperi_max.gpsNanoSeconds;
  GV->mismatch=CLA.mismatch;
  GV->mismatchedflag=CLA.mismatchedflag;
  GV->exactflag=CLA.exactflag;
  sprintf(GV->ifo,CLA.ifo);
  sprintf(GV->ephemdir,CLA.ephemdir);
  sprintf(GV->yr,CLA.yr);
  sprintf(GV->meshdir,CLA.meshdir);
  sprintf(GV->source,CLA.source);

  /* allocate memory for the number of bands */
  if (sourceparams.freq.f_err[0]!=0.0) {
    GV->nband=(INT4)ceil((sourceparams.freq.f_max[0]-sourceparams.freq.f_min[0])/GV->band);
  }
  else {
    GV->nband=1;
  }
  
  if (sourceparams.freq.nband==2) GV->nband+=(INT4)ceil((sourceparams.freq.f_max[1]-sourceparams.freq.f_min[1])/GV->band);
  GV->f_max=(REAL8 *)LALMalloc(GV->nband*sizeof(REAL8));

  j=0;
  /* if we have two bands */
  if (sourceparams.freq.nband==2) {
    f=sourceparams.freq.f_max[1];

    while (f>sourceparams.freq.f_min[1]) {
    
      GV->f_max[j]=f;
      f=f-GV->band;
      j++;
    }
  }

  /* if we have one or two  bands */
  f=sourceparams.freq.f_max[0];
  GV->f_max[0]=f;
  while (f>sourceparams.freq.f_min[0]) {    
    GV->f_max[j]=f;
    f=f-GV->band;
    j++;
  }

  return 0;

}

/***********************************************************************************/

int ReadDataParams(char *datadir, GlobVar *GV)
{

  /* the main point of this function is to extract the observation span */
  /* and observation start time at the detector site */

  INT4 fileno=0;
  FILE *fp;
  size_t errorcode;
  char **filelist;
  char command[512];
  glob_t globbuf;
  UINT4 i;
  INT4 nfiles;
  

  /* set up the datadir name */
  strcpy(command, datadir);
  strcat(command,"/*");
    
  /* set up some glob stuff */
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);
  
  /* check if there are any SFT's in the directory */
  if(globbuf.gl_pathc==0)
    {
      fprintf (stderr,"\nNo SFTs in directory %s ... Exiting.\n", datadir);
      exit(1);
    }
  
  /* allocate memory for the pathnames */
  filelist=(char **)LALMalloc(globbuf.gl_pathc*sizeof(char *));
  for (i=0;i<globbuf.gl_pathc;i++) filelist[i]=(char *)LALMalloc(256*sizeof(char));
  
  /* read all file names into memory */
  while ((UINT4)fileno < globbuf.gl_pathc) 
    {
      strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
      if (fileno > MAXSFTS)
	{
	  fprintf(stderr,"\nToo many files in directory! Exiting... \n");
	  exit(1);
	}
    }
  globfree(&globbuf);

  nfiles=fileno;
  
  /* open the first file to read header information */
  if (!(fp=fopen(filelist[0],"rb"))) {
	fprintf(stderr,"Weird... %s doesn't exist!\n",filelist[0]);
	return 1;
      }
      
  /* Read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",filelist[0]);
      return 1;
    }
  
  /* Check that the time base is positive */
  if (header.tbase<=0.0)
    {
      fprintf(stderr,"Timebase %f from data file %s non-positive!\n",
	      header.tbase,filelist[fileno]);
      return 3;
    }

  /* read in start time */
  GV->tstart.gpsSeconds=header.gps_sec;
  GV->tstart.gpsNanoSeconds=0;

  /* close the test file */
  fclose(fp);

  /* open up the last file in the list */
   if (!(fp=fopen(filelist[nfiles-1],"rb"))) {
	fprintf(stderr,"Weird... %s doesn't exist!\n",filelist[nfiles-1]);
	return 1;
      }
      
  /* Read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",filelist[fileno]);
      return 1;
    }


  /* read in end time */
  GV->tspan=(header.gps_sec-GV->tstart.gpsSeconds)+header.tbase;
  LALFree(filelist);


  fclose(fp);
  
  return 0;

}

/***********************************************************************************/

int FreeMem(EphemerisData *edat)
{
 
  /* Here we just free up the memory that has yet to be freed */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  return 0;
 
}

/***********************************************************************************/

int GenerateMesh(GlobVar GV,REAL4VectorSequence **XYmesh,XYparameterspace *XYspace,Metric *XYMetric)
{
  
  /* this function lays the mesh by calling some of the Flatmesh routines */

  static FlatMeshParamStruc FMparams;
  UINT4Vector *dimlength=NULL;
  REAL4Array *Matrix=NULL;
  REAL4Array *InvMatrix=NULL;
  REAL4VectorSequence *MatrixVec=NULL;
  REAL4VectorSequence *InvMatrixVec=NULL;
  REAL4VectorSequence *Corners=NULL;
  REAL4Vector *vector=NULL;
  REAL4 *dummy=NULL;
  RandomParams *params=NULL;
  CreateVectorSequenceIn in;
  INT4 seed=0;
  UINT4 i;
  REAL8 a,b;

  /* allocate memory for matrix and inverse in sequence form (for FlatMesh routines) */
  in.length=DIM;
  in.vectorLength=DIM;
  LALSCreateVectorSequence(&status,&MatrixVec,&in);
  LALSCreateVectorSequence(&status,&InvMatrixVec,&in);
  LALSCreateVectorSequence(&status,&Corners,&in);
 
  /* allocate memory for the array forms of matrix and inverse (for matrix routines) */
  LALU4CreateVector(&status,&dimlength,(UINT4)DIM);
  dimlength->data[0]=DIM;
  dimlength->data[1]=DIM;
  LALSCreateArray(&status,&Matrix,dimlength);
  LALSCreateArray(&status,&InvMatrix,dimlength);
  LALU4DestroyVector(&status,&dimlength);

  /* here we correctly fill the matrix array with correct normalisation (FlatMesh.h lal doc) */
  /* also fill it in correct order for use in project routine */
  Matrix->data[0]=2.0*((REAL4)sqrt(GV.mismatch))*XYMetric->eigenvec->data[0]/sqrt(abs((REAL4)DIM*XYMetric->eigenval->data[0]));
  Matrix->data[1]=2.0*((REAL4)sqrt(GV.mismatch))*XYMetric->eigenvec->data[2]/sqrt(abs((REAL4)DIM*XYMetric->eigenval->data[0]));
  Matrix->data[2]=2.0*((REAL4)sqrt(GV.mismatch))*XYMetric->eigenvec->data[1]/sqrt(abs((REAL4)DIM*XYMetric->eigenval->data[1]));
  Matrix->data[3]=2.0*((REAL4)sqrt(GV.mismatch))*XYMetric->eigenvec->data[3]/sqrt(abs((REAL4)DIM*XYMetric->eigenval->data[1]));

  /* copy contents of matrix array into matrix sequence */
  MatrixVec->data[0]=Matrix->data[0];
  MatrixVec->data[1]=Matrix->data[1];
  MatrixVec->data[2]=Matrix->data[2];
  MatrixVec->data[3]=Matrix->data[3];

  /* Invert the matrix (Matrix functions require it in array form) */
  /* it corrupts the contents of Matrix */
  LALSMatrixInverse(&status,dummy,Matrix,InvMatrix);

  /* copy the contents of the Inverse array to the inverse vector sequence */  
  InvMatrixVec->data[0]=InvMatrix->data[0];
  InvMatrixVec->data[1]=InvMatrix->data[1];
  InvMatrixVec->data[2]=InvMatrix->data[2];
  InvMatrixVec->data[3]=InvMatrix->data[3];

  /* destroy the array form of the matrix and inverse */
  LALSDestroyArray(&status,&Matrix);
  LALSDestroyArray(&status,&InvMatrix);

  /* if we are generating a single mismatched template */
  if (GV.mismatchedflag) {

    /* allocate memory for the single template */
    in.length=1;
    in.vectorLength=DIM;
    LALCreateVectorSequence(&status,XYmesh,&in);

    /* generate two random numbers */
    LALCreateVector(&status,&vector,2);
    LALCreateRandomParams(&status,&params,seed);

    /* fill vector with random sequence between -1 -> 1 */
    for (i=0;i<vector->length;i++) {
      LALUniformDeviate(&status,vector->data+i,params);
      vector->data[i]=-1.0 + 2.0*vector->data[i];
    }
    
    /* find the length of the box */
    a=vector->data[0]*sqrt(GV.mismatch/(2.0*XYMetric->eigenval->data[0]));
    b=vector->data[1]*sqrt(GV.mismatch/(2.0*XYMetric->eigenval->data[1]));
    
    (*XYmesh)->length=1;
    (*XYmesh)->data[0]=XYspace->X_0+(a*XYMetric->eigenvec->data[0]+b*XYMetric->eigenvec->data[1]);
    (*XYmesh)->data[1]=XYspace->Y_0+(a*XYMetric->eigenvec->data[2]+b*XYMetric->eigenvec->data[3]);

    /*printf("Xo and Y0 = %f %f\n",XYspace->X_0,XYspace->Y_0);
      printf("mismatched = %f %f\n",(*XYmesh)->data[0],(*XYmesh)->data[1]);*/

  }
  /* else we are generating an exact template */
  else if (GV.exactflag) {

    /* allocate memory for the single template */
    in.length=1;
    in.vectorLength=DIM;
    LALCreateVectorSequence(&status,XYmesh,&in);
     
    (*XYmesh)->length=1;
    (*XYmesh)->data[0]=XYspace->X_0;
    (*XYmesh)->data[1]=XYspace->Y_0;

  }
  /* else we are actually making a mesh */
  else {
 
    /* fill in the corners stucture */
    Corners->data[0]=(REAL4)XYspace->X_MIN;
    Corners->data[1]=(REAL4)XYspace->Y_MIN;
    Corners->data[2]=(REAL4)XYspace->X_MAX;
    Corners->data[3]=(REAL4)XYspace->Y_MAX;
    
    /* fill in the params structure */
    FMparams.matrix=MatrixVec;
    FMparams.matrixInv=InvMatrixVec;
    FMparams.xMin=NULL;
    FMparams.xMax=NULL;
    FMparams.controlPoints = Corners;
    FMparams.intersection = LALRectIntersect;
    
    /* this next bit just copies corners to Xmin and Xmax structures */
    LALSCreateVector(&status, &(FMparams.xMin), (UINT4)DIM );
    LALSCreateVector(&status, &(FMparams.xMax), (UINT4)DIM );
    FMparams.xMin->data[0]=Corners->data[0];
    FMparams.xMin->data[1]=Corners->data[1];
    FMparams.xMax->data[0]=Corners->data[2];
    FMparams.xMax->data[1]=Corners->data[3];

    /* Compute the mesh, and clean up local memory. */
    LALCreateFlatMesh(&status,XYmesh, &FMparams );
    LALSDestroyVector(&status,&(FMparams.xMin));
    LALSDestroyVector(&status,&(FMparams.xMax));
    LALSDestroyVectorSequence(&status,&Corners);

    /* check here if we have actually generated any points in XY space */
    if ((*XYmesh)->length<1) {
      fprintf(stderr,"ERROR : No points have been generated in the XY space\n");
      fprintf(stderr,"        Space may be too small -> Could be a one filter target\n");
      fprintf(stderr,"        Further investigation is required\n");
      fprintf(stderr,"          We will output a single filter in the center of the space\n");

      LALDestroyVectorSequence(&status,XYmesh);
      
      /* allocate memory for the single template */
      in.length=1;
      in.vectorLength=DIM;
      LALCreateVectorSequence(&status,XYmesh,&in);

      /* define the single template */
      (*XYmesh)->length=1;
      (*XYmesh)->data[0]=XYspace->X_0;
      (*XYmesh)->data[1]=XYspace->Y_0;

    }

  }
  
  LALSDestroyVectorSequence(&status,&MatrixVec);
  LALSDestroyVectorSequence(&status,&InvMatrixVec);

  return 0;

}

/***********************************************************************************/

int SetupBaryInput(GlobVar *GV,EphemerisData **edat,LALDetector *Detector)
{
  
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4 leap;
  CHAR filenameE[256],filenameS[256];
  FILE *fp;

  /* make the full file name/location for the ephemeris files */
  strcpy(filenameE,GV->ephemdir);
  strcat(filenameE,"/earth");
  strcat(filenameE,GV->yr);
  strcat(filenameE,".dat");
  strcpy(filenameS,GV->ephemdir);
  strcat(filenameS,"/sun");
  strcat(filenameS,GV->yr);
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
  (*edat)=(EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*(*edat)).ephiles.earthEphemeris = filenameE;     
  (*(*edat)).ephiles.sunEphemeris = filenameS;   

  /* set up leap second information */
  LALLeapSecs(&status,&leap,&GV->tstart,&formatAndAcc);
  (*(*edat)).leap=leap;

  /* Read in ephemeris files */
  LALInitBarycenter(&status,(*edat));             

  /* setup chosen detector */
  if(strcmp(GV->ifo,"GEO")!=0) *Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if(strcmp(GV->ifo,"LLO")!=0) *Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(strcmp(GV->ifo,"LHO")!=0) *Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];

  return 0;

}

/***********************************************************************************/

int GetSSBTime(GlobVar GV,LIGOTimeGPS *tdet, LIGOTimeGPS *tssb,EphemerisData **edat,LALDetector Detector)
{

  EarthState earth;
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */ 
  EmissionTime emit;

/* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.tgps.gpsSeconds=tdet->gpsSeconds;
  baryinput.tgps.gpsNanoSeconds=tdet->gpsNanoSeconds;
  baryinput.site.location[0]=Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=Detector.location[2]/LAL_C_SI;
  baryinput.alpha=GV.RA;
  baryinput.delta=GV.dec;
  baryinput.dInv=0.e0;

  /* Setup barycentric time at start of observation (first timestamp) */
  LALBarycenterEarth(&status,&earth,tdet,(*edat));
  LALBarycenter(&status,&emit,&baryinput,&earth);
  tssb->gpsSeconds=(emit.te.gpsSeconds);
  tssb->gpsNanoSeconds=emit.te.gpsNanoSeconds;
  
  return 0;
}

/*******************************************************************************/

int SetupPspaceParams(GlobVar GV,RTparameterspace *RTspace,XYparameterspace *XYspace)
{

  /* this function takes the boundary information given in sma,tperi coordinates   */
  /* and calculates the new overestimate of the boundary in X,Y space.             */
  /* It also converts the central values of sma,tperi into corresponding central   */
  /* values of X and Y. */
  
  RTPLocation RTPloc;
  XYLocation XYloc;
  REAL8Vector *Xtemp=NULL;
  REAL8Vector *Ytemp=NULL;
  REAL8Vector *alpha_temp=NULL;
  REAL8 alpha_MAX;
  REAL8 alpha_MIN;
  REAL8 dummy;
  INT4 i;
  UINT4 TWODIM;
  REAL8 tperi_err;

  /* here we simply adjust the periapse passage time so that it is within */
  /* a single orbit of the observation start time.  We do not account for */
  /* error propogation, the user must do this by hand and be present in the */
  /* search parameters error ranges. */
  
  /* here we call PeripaseShift to find the most recent periapse passage */
  if (PeriapseShift(GV.tperi_0,&(tperi_0_new),GV.tstartSSB,GV.period,&NORB)) return 2;
  /*if (PeriapseShift(GV.tperi_MIN,&(tperi_MIN_new),GV.tstartSSB,GV.period)) return 2;
    if (PeriapseShift(GV.tperi_MAX,&(tperi_MAX_new),GV.tstartSSB,GV.period)) return 2;*/

  /* try to solve straddling problem where pspace straddles a periapse passage time */
  LALDeltaFloatGPS(&status,&tperi_err,&(GV.tperi_MIN),&(GV.tperi_0));
  LALAddFloatToGPS(&status,&tperi_MIN_new,&tperi_0_new,tperi_err);
  LALDeltaFloatGPS(&status,&tperi_err,&(GV.tperi_MAX),&(GV.tperi_0));
  LALAddFloatToGPS(&status,&tperi_MAX_new,&tperi_0_new,tperi_err);

  /*  printf("old tperi 0 %d %d\n",GV.tperi_0.gpsSeconds,GV.tperi_0.gpsNanoSeconds);
  printf("new tperi 0 %d %d\n",tperi_0_new.gpsSeconds,tperi_0_new.gpsNanoSeconds);
  printf("old tperi min %d %d\n",GV.tperi_MIN.gpsSeconds,GV.tperi_MIN.gpsNanoSeconds);
  printf("new tperi min %d %d\n",tperi_MIN_new.gpsSeconds,tperi_MIN_new.gpsNanoSeconds);
  printf("old tperi max %d %d\n",GV.tperi_MAX.gpsSeconds,GV.tperi_MAX.gpsNanoSeconds);
  printf("new tperi max %d %d\n",tperi_MAX_new.gpsSeconds,tperi_MAX_new.gpsNanoSeconds);*/
 

  /* start filling in the RTspace structure */
  RTspace->sma_0=GV.sma_0;
  RTspace->sma_MIN=GV.sma_MIN;
  RTspace->sma_MAX=GV.sma_MAX;
  RTspace->tperi_0.gpsSeconds=tperi_0_new.gpsSeconds;
  RTspace->tperi_0.gpsNanoSeconds=tperi_0_new.gpsNanoSeconds;
  RTspace->tperi_MIN.gpsSeconds=tperi_MIN_new.gpsSeconds;
  RTspace->tperi_MIN.gpsNanoSeconds=tperi_MIN_new.gpsNanoSeconds;
  RTspace->tperi_MAX.gpsSeconds=tperi_MAX_new.gpsSeconds;
  RTspace->tperi_MAX.gpsNanoSeconds=tperi_MAX_new.gpsNanoSeconds;

  /*printf("RTspace stuff is %d %d - %d %d - %d %d\n",RTspace->tperi_0.gpsSeconds,RTspace->tperi_0.gpsNanoSeconds,RTspace->tperi_MIN.gpsSeconds, RTspace->tperi_MIN.gpsNanoSeconds,RTspace->tperi_MAX.gpsSeconds,RTspace->tperi_MAX.gpsNanoSeconds);
    printf("RTspace stuff is %f %f %f\n",RTspace->sma_0,RTspace->sma_MIN,RTspace->sma_MAX);*/

  /* allocate memory to Xtemp and Ytemp */
  TWODIM=2*(UINT4)DIM;
  LALDCreateVector(&status,&Xtemp,TWODIM);
  LALDCreateVector(&status,&Ytemp,TWODIM);
  LALDCreateVector(&status,&alpha_temp,TWODIM);

  /* fill up RTP location structures */
  RTPloc.period=GV.period;
  RTPloc.ecc=ECC;   /* this is set to zero in the header */
  RTPloc.argp=ARGP; /* this is set to zero in the header */
  RTPloc.tstartSSB.gpsSeconds=GV.tstartSSB.gpsSeconds;
  RTPloc.tstartSSB.gpsNanoSeconds=GV.tstartSSB.gpsNanoSeconds;
  
  /* convert all 4 corners of the RT space to corners in the XY space */
  RTPloc.sma=RTspace->sma_MIN;
  RTPloc.tperi.gpsSeconds=RTspace->tperi_MIN.gpsSeconds;
  RTPloc.tperi.gpsNanoSeconds=RTspace->tperi_MIN.gpsNanoSeconds;
  if (ConvertRTperitoXY(&RTPloc,&XYloc,&alpha_temp->data[0])) return 1;
  Xtemp->data[0]=XYloc.X;
  Ytemp->data[0]=XYloc.Y;
  RTPloc.sma=RTspace->sma_MIN;
  RTPloc.tperi.gpsSeconds=RTspace->tperi_MAX.gpsSeconds;
  RTPloc.tperi.gpsNanoSeconds=RTspace->tperi_MAX.gpsNanoSeconds;
  if (ConvertRTperitoXY(&RTPloc,&XYloc,&alpha_temp->data[1])) return 2;
  Xtemp->data[1]=XYloc.X;
  Ytemp->data[1]=XYloc.Y;						
  RTPloc.sma=RTspace->sma_MAX;
  RTPloc.tperi.gpsSeconds=RTspace->tperi_MIN.gpsSeconds;
  RTPloc.tperi.gpsNanoSeconds=RTspace->tperi_MIN.gpsNanoSeconds;
  if (ConvertRTperitoXY(&RTPloc,&XYloc,&alpha_temp->data[2])) return 3;
  Xtemp->data[2]=XYloc.X;
  Ytemp->data[2]=XYloc.Y;
  RTPloc.sma=RTspace->sma_MAX;
  RTPloc.tperi.gpsSeconds=RTspace->tperi_MAX.gpsSeconds;
  RTPloc.tperi.gpsNanoSeconds=RTspace->tperi_MAX.gpsNanoSeconds;
  if (ConvertRTperitoXY(&RTPloc,&XYloc,&alpha_temp->data[3])) return 3;
  Xtemp->data[3]=XYloc.X;
  Ytemp->data[3]=XYloc.Y;
  
  /* initialise the XYspace bouondaries */
  for (i=0;i<4;i++) {
    XYspace->X_MAX=Xtemp->data[0];
    XYspace->X_MIN=Xtemp->data[0];
    XYspace->Y_MAX=Ytemp->data[0];
    XYspace->Y_MIN=Ytemp->data[0];
    alpha_MIN=alpha_temp->data[0];
    alpha_MAX=alpha_temp->data[0];
  }

  /* loop over the input corners and find the minimum and maximum */
  for (i=0;i<4;i++) {
    if (Xtemp->data[i]>XYspace->X_MAX) XYspace->X_MAX=Xtemp->data[i];
    if (Xtemp->data[i]<XYspace->X_MIN) XYspace->X_MIN=Xtemp->data[i];
    if (Ytemp->data[i]>XYspace->Y_MAX) XYspace->Y_MAX=Ytemp->data[i];
    if (Ytemp->data[i]<XYspace->Y_MIN) XYspace->Y_MIN=Ytemp->data[i];
  }

  /* put initial phases in order (alpha lies in range 0 - 2PI) */
   for (i=0;i<4;i++) {
    if (alpha_temp->data[i]>alpha_MAX) alpha_MAX=alpha_temp->data[i];
    if (alpha_temp->data[i]<alpha_MIN) alpha_MIN=alpha_temp->data[i];
   }
   
   /* clean up somoe temporary memory space */
   LALDDestroyVector(&status,&Xtemp);
   LALDDestroyVector(&status,&Ytemp);
   LALDDestroyVector(&status,&alpha_temp);
   
   /* check if we cross a multiple of PI/2 in which case we have to redefine max and min boundaries */
   if ((alpha_MIN>0.0)&&(alpha_MIN<LAL_PI/2.0)&&(alpha_MAX>3.0*LAL_PI/2.0)&&(alpha_MAX<LAL_TWOPI)) XYspace->X_MAX=RTspace->sma_MAX;
   if ((alpha_MIN>0.0)&&(alpha_MIN<LAL_PI/2.0)&&(alpha_MAX>LAL_PI/2.0)&&(alpha_MAX<LAL_PI)) XYspace->Y_MAX=RTspace->sma_MAX;
   if ((alpha_MIN>LAL_PI/2.0)&&(alpha_MIN<LAL_PI)&&(alpha_MAX>LAL_PI)&&(alpha_MAX<3.0*LAL_PI/2.0)) XYspace->X_MIN=(-1.0)*RTspace->sma_MAX;
   if ((alpha_MIN>LAL_PI)&&(alpha_MIN<3.0*LAL_PI/2.0)&&(alpha_MAX>3.0*LAL_PI/2.0)&&(alpha_MAX<LAL_TWOPI)) XYspace->Y_MIN=(-1.0)*RTspace->sma_MAX;
   
   /* find centre of XY space */
   RTPloc.sma=RTspace->sma_0;
   RTPloc.tperi.gpsSeconds=RTspace->tperi_0.gpsSeconds;
   RTPloc.tperi.gpsNanoSeconds=RTspace->tperi_0.gpsNanoSeconds;
   if (ConvertRTperitoXY(&RTPloc,&XYloc,&dummy)) return 4;
   XYspace->X_0=XYloc.X;
   XYspace->Y_0=XYloc.Y;
   
   return 0;
   
}

/***********************************************************************************/

int GenMetricComp(GlobVar GV,REAL8 f_max,REAL8 *X,REAL8 *Y,Metric *XYMetric)
{
 
  /* This function calculates the metric elements of the g-metric and then          */
  /* projects to find the gamma-metric elements.  It diagonalises the gamma metric  */
  /* and returns eigenvalues and eigenvectors plus other things in a metric         */
  /* structure */
  
  REAL8 dphi[3];
  REAL8 dphisq[3];
  REAL8 dphidphi[3];
  UINT4Vector *dimlength=NULL;
  REAL8Vector *gMetric=NULL;
  REAL8Array *elementtemp=NULL;
  REAL8 T;
  REAL8 w;
  REAL8 f;
  REAL8 pi;
  INT4 i,j;

  /* set up some more concise variables */
  T=GV.tspan;
  f=(REAL8)f_max;
  w=LAL_TWOPI/(GV.period);
  pi=LAL_PI;
  
  /* here we calculate each element of the unprojected metric */
  
  /* <df> */
  dphi[0]=(2.0*pi/T)*((T*T/2.0)+(*X/w)*(cos(w*T)-1.0) - ((*Y)/w)*sin(w*T) );

  /* <dX> */
  dphi[1]=((2.0*pi*f)/(T*w)) * (1.0-cos(w*T));

  /* <dY> */
  dphi[2]=((2.0*pi*f)/(T*w))*sin(w*T);

  /* <dfdf> */
  dphisq[0]=(4.0*pi*pi/T) * ( (T*T*T/3.0) + (2.0*(*X)*T*cos(w*T)/w) - (2.0*(*X)*sin(w*T)/(w*w)) - (2.0*(*Y)*T*sin(w*T)/w) - (2.0*(*Y)*cos(w*T)/(w*w)) + (2.0*(*Y)/(w*w)) + ((*X)*(*X)*T/2.0) - ((*X)*(*X)*sin(2.0*w*T)/(4.0*w)) + ((*Y)*(*Y)*T/2.0) + ((*Y)*(*Y)*sin(2.0*w*T)/(4.0*w)) - ((*X)*(*Y)*cos(2.0*w*T)/(2.0*w)) + ((*X)*(*Y)/(2.0*w)) );

  /* <dXdX> */
  dphisq[1]=(2.0*pi*pi*f*f/T)*(T-(sin(2.0*w*T)/(2.0*w)));

  /* <dYdY> */
  dphisq[2]=(2.0*pi*pi*f*f/T)*(T+(sin(2.0*w*T)/(2.0*w)));

  /* <dfdX> */
  dphidphi[0]=(4.0*pi*pi*f/T)*( (sin(w*T)/(w*w)) - (T*cos(w*T)/w) - ((*X)*T/2.0) + ((*X)*sin(2.0*w*T)/(4.0*w)) + ((*Y)*cos(2.0*w*T)/(4.0*w)) - ((*Y)/(4.0*w)) );

  /* <dXdY> */
  dphidphi[1]=(pi*pi*f*f/(T*w))*(1.0-cos(2.0*w*T));

  /* <dfdY> */
  dphidphi[2]=(4.0*pi*pi*f/T) * ( (T*sin(w*T)/w) + (cos(w*T)/(w*w)) - (1.0/(w*w)) + ((*X)*cos(2.0*w*T)/(4.0*w)) - ((*X)/(4.0*w)) - ((*Y)*T/2.0) - ((*Y)*sin(2.0*w*T)/(4.0*w)) );

  /* allocate some memory for the metric structure */
  LALU4CreateVector(&status,&dimlength,DIM);
  dimlength->data[0]=DIM;
  dimlength->data[1]=DIM;
  XYMetric->element=NULL;
  XYMetric->eigenvec=NULL;
  XYMetric->eigenval=NULL;
  LALDCreateArray(&status,&(XYMetric->element),dimlength);
  LALDCreateArray(&status,&(elementtemp),dimlength);
  LALDCreateArray(&status,&(XYMetric->eigenvec),dimlength);
  dimlength->data[0]=DIM+1;
  dimlength->data[1]=DIM+1;
  LALDCreateVector(&status,&gMetric,(DIM+1)*(DIM+2)/2);
  LALDCreateVector(&status,&(XYMetric->eigenval),DIM);
  LALU4DestroyVector(&status,&dimlength);

  /* now we calculate the elements of the g-metric */
  /* they are stored in the order defined in CoherentMetric.c */
  gMetric->data[0]=dphisq[0]-(dphi[0]*dphi[0]);  /* gff */
  gMetric->data[2]=dphisq[1]-(dphi[1]*dphi[1]);  /* gXX */
  gMetric->data[5]=dphisq[2]-(dphi[2]*dphi[2]);  /* gYY */
  gMetric->data[1]=dphidphi[0]-(dphi[0]*dphi[1]); /* gfX */
  gMetric->data[4]=dphidphi[1]-(dphi[1]*dphi[2]); /* gXY */
  gMetric->data[3]=dphidphi[2]-(dphi[2]*dphi[0]); /* gfY */
  
  /* Now we project out the frequency component */
  LALProjectMetric(&status,gMetric,0);

  /*printf("generated metric gamma elements as %f %f\n",gMetric->data[2],gMetric->data[4]);
    printf("generated metric gamma elements as %f %f\n",gMetric->data[4],gMetric->data[5]);*/

  /* put the now projected metric into the correct location */
  /* stored normally with redundancy because symmetric */
  XYMetric->element->data[0]=gMetric->data[2];
  XYMetric->element->data[1]=gMetric->data[4];
  XYMetric->element->data[2]=gMetric->data[4];
  XYMetric->element->data[3]=gMetric->data[5];

  /*printf("generated metric gamma elements as %f %f\n",XYMetric->element->data[0],XYMetric->element->data[1]);
    printf("generated metric gamma elements as %f %f\n",XYMetric->element->data[2],XYMetric->element->data[3]);*/

  /* copy the information to a temporary variable before hand */
  elementtemp->dimLength[0]=XYMetric->element->dimLength[0];
  elementtemp->dimLength[1]=XYMetric->element->dimLength[1];
  elementtemp->data[0]=XYMetric->element->data[0];
  elementtemp->data[1]=XYMetric->element->data[1];
  elementtemp->data[2]=XYMetric->element->data[2];
  elementtemp->data[3]=XYMetric->element->data[3];
  
  /* calculate the determinant of the matrix */
  LALDMatrixDeterminant(&status,&(XYMetric->determinant),(elementtemp));
 
  /* big check here ! must stop if determinant is NOT positive */
  if (XYMetric->determinant<=0.0) {
    fprintf(stderr,"ERROR : metric determinant < 0.0 -> observation time span probably too low\n");
    exit(1);
  }   

  /* copy the metric information to the eigenvec variable beforehand */
  for (i=0;i<DIM;i++) {
    for (j=0;j<DIM;j++) {
      XYMetric->eigenvec->data[i*DIM+j]=XYMetric->element->data[i*DIM+j];
    }
    XYMetric->eigenvec->dimLength[i]=XYMetric->element->dimLength[i];
  }

  /* here we find the eigenvalues and eigenvectors */
  LALDSymmetricEigenVectors(&status,XYMetric->eigenval,XYMetric->eigenvec);

  /* then put elements back into the metric structure */
  LALDDestroyVector(&status,&gMetric);
  /* LALDDestroyArray(&status,&elementtemp);*/ /* removed this, dont know why, but it stops the thing crashing */

  /* printf("generated metric gamma elements as %f %f\n",XYMetric->element->data[0],XYMetric->element->data[1]);
     printf("generated metric gamma elements as %f %f\n",XYMetric->element->data[2],XYMetric->element->data[3]); */

  return 0;

}

/***********************************************************************************/

int CheckRTBoundary(REAL8 *sma_temp,LIGOTimeGPS *tp_temp,RTparameterspace *RTspace)
{

  /* this routine simply checks wether a point in RT space lies within the space boundaries */

  INT4 result;

  if (*sma_temp<(RTspace->sma_MIN)) return 0;
  if (*sma_temp>(RTspace->sma_MAX)) return 0;
  LALCompareGPS(&status,&result,tp_temp,&(RTspace->tperi_MIN));
  if (result==-1) return 0;
  LALCompareGPS(&status,&result,tp_temp,&(RTspace->tperi_MAX));
  if (result==1) return 0;

  return 1;

}

/***********************************************************************************/

int ConvertMesh(GlobVar GV,REAL4VectorSequence **XYmesh,RTMesh *RTmesh,RTparameterspace *RTspace)
{

  /* this section coverts the XY mesh to a RT mesh and retains only those points */
  /* within the original RT boundary. */
  
  RTPLocation RTPloc;
  XYLocation XYloc;
  LIGOTimeGPS *tperi_vec=NULL;
  LIGOTimeGPS temp;
  REAL8Vector *sma_vec=NULL;
  UINT4 i,j;
  INT4 *tempNORB=NULL;

  /* first allocate some memory for the new mesh (using all temporary local variables */
  if ((*XYmesh)->length<1) {
    printf("ERROR : Zero length XY mesh, strange !!!");
    exit(1);
  }
  RTmesh->length=(*XYmesh)->length;
  RTmesh->sma=NULL;
  LALDCreateVector(&status,&(RTmesh->sma),RTmesh->length);

  /* point the temp pointer to the allocated space */
  sma_vec=RTmesh->sma;
  RTmesh->tperi = (LIGOTimeGPS *)LALCalloc( RTmesh->length , sizeof(LIGOTimeGPS) );

  j=0;
  /* point the temp pointer to the alocated space */
  tperi_vec=RTmesh->tperi;
  
  /* fill in the XY location structure */
  XYloc.period=GV.period;
  XYloc.tstartSSB.gpsSeconds=GV.tstartSSB.gpsSeconds;
  XYloc.tstartSSB.gpsNanoSeconds=GV.tstartSSB.gpsNanoSeconds;
  XYloc.ecc=ECC;
  XYloc.argp=ARGP;
  
 
  /* covert XY to RT mesh and trim */
  for (i=0;i<(*XYmesh)->length;i++) {
    XYloc.X=(REAL8)(*XYmesh)->data[i*2];
    XYloc.Y=(REAL8)(*XYmesh)->data[i*2+1];
    
    /* convert this point in XY space to RT space */
    ConvertXYtoRTperi(&XYloc,&RTPloc);
    /*printf("converted XY %f %f -> RT %f %d %d\n",XYloc.X,XYloc.Y,RTPloc.sma,RTPloc.tperi.gpsSeconds,RTPloc.tperi.gpsNanoSeconds);*/

    
    /* Now check if this RT point lies within the original boundaries */
    /* if it does save it to memory for later output */
    if (CheckRTBoundary(&RTPloc.sma,&RTPloc.tperi,RTspace)==1) {
      
      /* shift the periapse passage time back to original input boundaries */
      if (PeriapseShiftBack(GV.tstartSSB,GV.tperi_0,RTPloc.tperi,&temp,GV.period,NORB)) return 3;
     
      RTPloc.tperi.gpsSeconds=temp.gpsSeconds;
      RTPloc.tperi.gpsNanoSeconds=temp.gpsNanoSeconds; 
      sma_vec->data[j]=RTPloc.sma;
      tperi_vec[j].gpsSeconds=RTPloc.tperi.gpsSeconds;
      tperi_vec[j].gpsNanoSeconds=RTPloc.tperi.gpsNanoSeconds;
      j++;
    }
  }
  
  if (j>0) {
    /* define new trimmed length */
    RTmesh->length=j;
    
    /* clean up the XY mesh memory */    
    LALSDestroyVectorSequence(&status,&(*XYmesh));
    
    /* change the length of the RTmesh array */
    LALDResizeVector(&status,&(RTmesh->sma),RTmesh->length);
    LALRealloc(RTmesh->tperi,RTmesh->length*sizeof(LIGOTimeGPS));
  }
  
  /* else if we are doing a mimatched template */
  else if (GV.mismatchedflag) {
    
    /* define new trimmed length */
    RTmesh->length=1;

    XYloc.X=(REAL8)(*XYmesh)->data[0];
    XYloc.Y=(REAL8)(*XYmesh)->data[1];
    
    /* convert this point in XY space to RT space */
    ConvertXYtoRTperi(&XYloc,&RTPloc);

    /* shift the periapse passage time back to original input boundaries */
    if (PeriapseShiftBack(GV.tstartSSB,GV.tperi_0,RTPloc.tperi,&temp,GV.period,NORB)) return 3;
    
    RTPloc.tperi.gpsSeconds=temp.gpsSeconds;
    RTPloc.tperi.gpsNanoSeconds=temp.gpsNanoSeconds; 
    sma_vec->data[j]=RTPloc.sma;
    tperi_vec[j].gpsSeconds=RTPloc.tperi.gpsSeconds;
    tperi_vec[j].gpsNanoSeconds=RTPloc.tperi.gpsNanoSeconds;

    /* clean up the XY mesh memory */    
    LALSDestroyVectorSequence(&status,&(*XYmesh));
    
    /* change the length of the RTmesh array */
    LALDResizeVector(&status,&(RTmesh->sma),RTmesh->length);
    LALRealloc(RTmesh->tperi,RTmesh->length*sizeof(LIGOTimeGPS));
    
    fprintf(stderr,"WARNING : a randomly mismatched signal has been put in the original parameter space.\n");
  
  }

  /* else we are doing an exact template or just a single template at the center */
  else {

    /* define new trimmed length */
    RTmesh->length=1;
    
    /* define the central template */
    sma_vec->data[0]=GV.sma_0;
    tperi_vec[0].gpsSeconds=GV.tperi_0.gpsSeconds;
    tperi_vec[0].gpsNanoSeconds=GV.tperi_0.gpsNanoSeconds;
    
    /* clean up the XY mesh memory */    
    LALSDestroyVectorSequence(&status,&(*XYmesh));
    
    /* change the length of the RTmesh array */
    LALDResizeVector(&status,&(RTmesh->sma),RTmesh->length);
    LALRealloc(RTmesh->tperi,RTmesh->length*sizeof(LIGOTimeGPS));
    
    fprintf(stderr,"WARNING : none of the points lie in the original parameter space.\n");
    fprintf(stderr,"          This could be a single filter target but futher investigation\n");
    fprintf(stderr,"          is required (edge effects !!)\n");
    fprintf(stderr,"          We will place a single template at the center of the space\n");
  }
  
  return 0;

}

/***********************************************************************************/

 int OutputRTMesh(GlobVar *GV,REAL8 f_max,RTMesh *RTmesh,Metric *XYMetric) 
{

  /* this section simply outputs the final mesh to file including the header information */

  FILE *fp;
  BinaryMeshFileHeader BMFheader;
  UINT4 i;
  char filename[256],ext[256],sourcename[256];

  /* generate file name for this fmax */
  strcpy(filename,GV->meshdir);
  strcpy(sourcename,GV->source);
  sprintf(ext,"/mesh_%s_%s_%.6f.mesh",GV->ifo,sourcename,f_max);
  strcat(filename,ext);

  fp=fopen(filename,"w");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",filename);
  }

  /* setup the header input */
  BMFheader.f_max=(REAL8)f_max;
  BMFheader.tspan=GV->tspan;
  BMFheader.tstart.gpsSeconds=GV->tstart.gpsSeconds;
  BMFheader.tstart.gpsNanoSeconds=GV->tstart.gpsNanoSeconds;
  BMFheader.Nfilters=RTmesh->length;
  BMFheader.mismatch=GV->mismatch;
  BMFheader.sma_0=GV->sma_0;
  BMFheader.sma_MIN=GV->sma_MIN;
  BMFheader.sma_MAX=GV->sma_MAX;
  BMFheader.tperi_0.gpsSeconds=GV->tperi_0.gpsSeconds;
  BMFheader.tperi_0.gpsNanoSeconds=GV->tperi_0.gpsNanoSeconds;
  BMFheader.tperi_MIN.gpsSeconds=GV->tperi_MIN.gpsSeconds;
  BMFheader.tperi_MIN.gpsNanoSeconds=GV->tperi_MIN.gpsNanoSeconds;
  BMFheader.tperi_MAX.gpsSeconds=GV->tperi_MAX.gpsSeconds;
  BMFheader.tperi_MAX.gpsNanoSeconds=GV->tperi_MAX.gpsNanoSeconds;
  BMFheader.ecc_MIN=ECC;
  BMFheader.ecc_MAX=ECC;
  BMFheader.argp_MIN=ARGP;
  BMFheader.argp_MAX=ARGP;
  BMFheader.period_MIN=GV->period;
  BMFheader.period_MAX=GV->period;
  BMFheader.metric_XX=XYMetric->element->data[0];
  BMFheader.metric_XY=XYMetric->element->data[1];
  BMFheader.metric_YY=XYMetric->element->data[3];
  sprintf(BMFheader.version,"v1");
  sprintf(BMFheader.det,GV->ifo);
  BMFheader.RA=GV->RA;
  BMFheader.dec=GV->dec;

  if (WriteMeshFileHeader(fp,&BMFheader)) return 1;

  /* output the filters in the form ready for a search */
  for (i=0;i<RTmesh->length;i++) {
   
    fprintf(fp,"%6.12f %6.12f %d %d %6.12f %6.12f\n", \
	    RTmesh->sma->data[i],GV->period,RTmesh->tperi[i].gpsSeconds, \
	    RTmesh->tperi[i].gpsNanoSeconds,ECC,ARGP);
  }

  fclose(fp);

  return 0;

}

/***********************************************************************************/

 int ReadCommandLine(int argc,char *argv[],CLargs *CLA) 
{
  INT4 c, errflg = 0;
  CHAR *temp;
  optarg = NULL;
  
  /* Initialize default values */
  sprintf(CLA->sourcefile,"sources.data");
  sprintf(CLA->source,"sco-x1");
  sprintf(CLA->datadir," ");
  CLA->tstart.gpsSeconds=0;
  CLA->tstart.gpsNanoSeconds=0;
  CLA->tspan=0.0;
  CLA->mismatch=0.0;
  sprintf(CLA->ephemdir,"./");
  sprintf(CLA->yr,"00-04");
  sprintf(CLA->ifo,"LLO");
  sprintf(CLA->meshdir,"./");
  CLA->datadirflag=0;
  CLA->mismatchedflag=0;
  CLA->exactflag=0;

  {
    int option_index = 0;
    static struct option long_options[] = {
      {"sourcefile", required_argument, 0, 'S'},
      {"source", required_argument, 0, 's'},
      {"datadir", required_argument, 0, 'D'},
      {"tstart", required_argument, 0, 'T'},
      {"tspan", required_argument, 0, 't'},
      {"mismatch", required_argument, 0, 'm'},
      {"band", required_argument, 0, 'b'},
      {"ephdir", required_argument, 0, 'E'},
      {"yr", required_argument, 0, 'y'},
      {"detector", required_argument, 0, 'I'},
      {"meshdir", required_argument, 0, 'o'},
      {"mismatched", no_argument, 0, 'X'},
      {"exact", no_argument, 0, 'x'},
      {"help", no_argument, 0, 'h'}
    };
    /* Scan through list of command line arguments */
    while (!errflg && ((c = getopt_long (argc, argv,"hS:s:D:T:t:m:b:E:y:I:o:Xx",long_options, &option_index)))!=-1)
      switch (c) {
      case 'S':
	temp=optarg;
	sprintf(CLA->sourcefile,temp);
	break;
      case 's':
	temp=optarg;
	sprintf(CLA->source,temp);
	break;
      case 'T':
	CLA->tstart.gpsSeconds=atoi(optarg);
	CLA->tstart.gpsNanoSeconds=0;
	break;
      case 't':
	CLA->tspan=atof(optarg);
	break;
      case 'D':
	temp=optarg;
	sprintf(CLA->datadir,temp);
	CLA->datadirflag=1;
	break;
      case 'm':
	CLA->mismatch=atof(optarg);
	break;
      case 'b':
	CLA->band=atof(optarg);
	break;
      case 'E':
	temp=optarg;
	sprintf(CLA->ephemdir,temp);
	break;
      case 'y':
	temp=optarg;
	sprintf(CLA->yr,temp);
	break;
      case 'I':
	temp=optarg;
	sprintf(CLA->ifo,temp);
	break;
      case 'o':
	temp=optarg;
	sprintf(CLA->meshdir,temp);
	break;
      case 'X':
	CLA->mismatchedflag=1;
	break;
      case 'x':
	CLA->exactflag=1;
	break;	
      case 'h':
	/* print usage/help message */
	fprintf(stdout,"Arguments are:\n");
	fprintf(stdout,"\t--sourcefile  STRING\t Name of the file containing source parameters [DEFAULT=]\n");
	fprintf(stdout,"\t--source      STRING\t Name of source [DEFAULT=]\n");
	fprintf(stdout,"\t--datadir     STRING\t Directory containing the data to be searched [DEFAULT=]\n");
	fprintf(stdout,"\t--tstart      INT4\t The start time of the observation (GPS) [DEFAULT=]\n");
	fprintf(stdout,"\t--tspan       REAL8\t The span of the observation (sec) [DEFAULT=]\n");
	fprintf(stdout,"\t--mismatch    REAL8\t The mismatch to be used [DEFAULT=]\n");
	fprintf(stdout,"\t--band        REAL8\t The size of the bands to be used [DEFAULT=]\n");
	fprintf(stdout,"\t--ephdir      STRING\t Location of ephemeris files earth?.dat and sun?.dat [DEFAULT=NULL]\n");
	fprintf(stdout,"\t--yr          STRING\t Year(s) specifying ephemeris files [DEFAULT=00-04]\n");
	fprintf(stdout,"\t--detector    STRING\t Interferometer being used for the search (LLO,LHO,GEO,TAMA,CIT,VIRGO) [DEFAULT=LLO]\n");
	fprintf(stdout,"\t--meshdir     STRING\t Name of output mesh file [DEFAULT=mesh.out]\n");
	fprintf(stdout,"\t--mismatched  BOOLEAN\t Set this flag if you require a single mismatched filter [DEFAULT=0]\n");
	fprintf(stdout,"\t--exact       BOOLEAN\t Set this flag if you require a single exactly matched filter [DEFAULT=0]\n");
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
  /* Need to add some CLA error checking here */
  
  /* update global variable and return */
  return errflg;
}

/************************************************************************/

int CheckInput(GlobVar GV)
{

  /* this routine takes the CLA inputs and does some basic validity checks */

  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  LALDate beginDate;
  LALDate endDate;
  LIGOTimeGPS beginGPS;
  LIGOTimeGPS endGPS;
  INT4 i;
 
  /* check that each of the f_max values are positive */
  for (i=0;i<GV.nband;i++) {
    if (GV.f_max[i]<0) {
      fprintf(stderr,"ERROR : MAX search frequency must be > 0 \n");
      exit(1);
    }
  }
  /* check that the time span is positive */
  if (GV.tspan<0.0) {
    fprintf(stderr,"ERROR : Observation time span must be > 0 \n");
    exit(1);
  }
  
  /* set up some unix date variables */
  beginDate.unixDate.tm_sec=0;
  beginDate.unixDate.tm_min=0;
  beginDate.unixDate.tm_hour=0;
  beginDate.unixDate.tm_mday=0;
  beginDate.unixDate.tm_mon=0;
  endDate.unixDate.tm_sec=0;
  endDate.unixDate.tm_min=0;
  endDate.unixDate.tm_hour=0;
  endDate.unixDate.tm_mday=0;
  endDate.unixDate.tm_mon=0;

  /* set the start and end times for the given year */
  if (strcmp(GV.yr,"98")==0) {
    beginDate.unixDate.tm_year=98;
    endDate.unixDate.tm_year=99;
  }
  if (strcmp(GV.yr,"99")==0) {
    beginDate.unixDate.tm_year=99;
    endDate.unixDate.tm_year=100;
  }
  if (strcmp(GV.yr,"00")==0) {
    beginDate.unixDate.tm_year=100;
    endDate.unixDate.tm_year=101;
  }
  if (strcmp(GV.yr,"01")==0) {
    beginDate.unixDate.tm_year=101;
    endDate.unixDate.tm_year=102;
  }
  if (strcmp(GV.yr,"02")==0) {
    beginDate.unixDate.tm_year=102;
    endDate.unixDate.tm_year=103;
  }
  if (strcmp(GV.yr,"03")==0) {
    beginDate.unixDate.tm_year=103;
    endDate.unixDate.tm_year=104;
  }
  if (strcmp(GV.yr,"04")==0) {
    beginDate.unixDate.tm_year=104;
    endDate.unixDate.tm_year=105;
  }
  if (strcmp(GV.yr,"05")==0) {
    beginDate.unixDate.tm_year=105;
    endDate.unixDate.tm_year=106;
  }
  if (strcmp(GV.yr,"06")==0) {
    beginDate.unixDate.tm_year=106;
    endDate.unixDate.tm_year=107;
  }
  if (strcmp(GV.yr,"07")==0) {
    beginDate.unixDate.tm_year=107;
    endDate.unixDate.tm_year=108;
  }
  if (strcmp(GV.yr,"08")==0) {
    beginDate.unixDate.tm_year=108;
    endDate.unixDate.tm_year=109;
  }
  if (strcmp(GV.yr,"09")==0) {
    beginDate.unixDate.tm_year=109;
    endDate.unixDate.tm_year=110;
  }
  if (strcmp(GV.yr,"00-04")==0) {
    beginDate.unixDate.tm_year=100;
    endDate.unixDate.tm_year=105;
  }
  if (strcmp(GV.yr,"03-06")==0) {
    beginDate.unixDate.tm_year=100;
    endDate.unixDate.tm_year=106;
  }
  if (strcmp(GV.yr,"05-09")==0) {
    beginDate.unixDate.tm_year=105;
    endDate.unixDate.tm_year=110;
  }
 
  /* convert the beginning and end of the relevant year(s) to a GPS time */
  LALUTCtoGPS(&status,&beginGPS,&beginDate,&(formatAndAcc.accuracy));
  LALUTCtoGPS(&status,&endGPS,&endDate,&(formatAndAcc.accuracy));
 
  /* check that the start time lies within the ephemeris span */
  if ((GV.tstart.gpsSeconds<beginGPS.gpsSeconds)||(GV.tstart.gpsSeconds+(INT4)GV.tspan>endGPS.gpsSeconds)) {
    fprintf(stderr,"Start time (+ observation span) must lie within time of ephemeris file\n");
    fprintf(stderr,"start time = %d beginGPS = %d endGPS = %d tobs = %f\n",GV.tstart.gpsSeconds,beginGPS.gpsSeconds,endGPS.gpsSeconds,GV.tspan);
    exit(1);
  }
  /* check that the RA is sensible */
  if ((GV.RA<0.0)||(GV.RA>LAL_TWOPI)) {
    fprintf(stderr,"ERROR : Source RA must be within range (0 -> 2PI) \n");
    exit(1);
  }
  /* check that the declination is sensible */
  if ((GV.dec<(-0.5)*LAL_PI)||(GV.dec>(0.5)*LAL_PI)) {
    fprintf(stderr,"ERROR : Source dec must be within range (-PI/2 -> PI/2) \n");
    exit(1);
  }
  /* check that the sma is positive */
  if (GV.sma_0<0.0) {
    fprintf(stderr,"ERROR : Central value of Orbital semi-major axis must be > 0 \n");
    exit(1);
  }
  /* check that the period is positive */
  if (GV.period<0.0) {
    fprintf(stderr,"ERROR : Orbital period must be > 0 \n");
    exit(1);
  }
  /* check that the mismatch is positive */
  if (GV.mismatch<0.0) {
    fprintf(stderr,"ERROR : Mismatch must be > 0 \n");
    exit(1);
  }
  /* check that the period is < 1 */
  if (GV.mismatch>1.0) {
    fprintf(stderr,"WARNING : Mismatch should be < 1\n");
    exit(1);
  }
  /* check that the detector name is valid */
  if ((strcmp(GV.ifo,"LLO")!=0)&&(strcmp(GV.ifo,"LHO")!=0)&&(strcmp(GV.ifo,"GEO")!=0)) {
    fprintf(stderr,"Not a known detector name \n");
    exit(1);
  }
  /* check that the periapse passage time lies within the ephemeris span */ 
  if ((GV.tperi_0.gpsSeconds<beginGPS.gpsSeconds)||(GV.tperi_0.gpsSeconds>endGPS.gpsSeconds)) {
    fprintf(stderr,"Central value of periapse passage time must lie within time of ephemeris file\n");
    exit(1);
  }
  /* check that the sma is wwithin its own boundaries */
  if ((GV.sma_MIN>GV.sma_0)||(GV.sma_MAX<GV.sma_0)) {
    fprintf(stderr,"Central value of the orbital semi-major axis not within MIN and MAX range \n");
    exit(1);
  }
  /* check that the periapse passage lies within its own range */
  if ((GV.tperi_MIN.gpsSeconds>GV.tperi_0.gpsSeconds)||(GV.tperi_MAX.gpsSeconds<GV.tperi_0.gpsSeconds)) {
    fprintf(stderr,"Central value of the periapse passage time not within MIN and MAX range \n");
    exit(1);
  }
  /* check that the maximum periapse time is greater than the minimum time */
  if (GV.tperi_MIN.gpsSeconds>GV.tperi_MAX.gpsSeconds) {
    fprintf(stderr,"MIN value of the periapse passage must be < MAX value \n");
    exit(1);
  }

  /* now for a few more physical checks */

  /* check for relativistic speeds in circular orbit */
  if ((GV.sma_MAX*LAL_TWOPI/GV.period)>0.01) {
    fprintf(stderr,"WARNING : MAX parameters indicate system is ~ relativistic !! \n");
  }
 
  /* more checks will be added as they are thought of */
  
  return 0;

}

/*******************************************************************************/
