/************************************************************************************/
/*         Program to generate a grid of filters in parameter space for use         */
/*         in a search for signals from sources in circular binary orbits.          */
/*         A metric is calculated, projected and diagonalised.  Filters are         */
/*         then placed within the specified boundaries of the space and output      */
/*         to file ready for a search.                                              */
/*                                                                                  */
/*			           C. Messenger                                     */
/*                                                                                  */
/*                         BIRMINGHAM UNIVERISTY -  2004                            */
/************************************************************************************/

#include "GenerateBinaryMesh_v1.h"

CLargs CLA;
REAL8 sma_GLOBAL;
LIGOTimeGPS tperi_GLOBAL;
LIGOTimeGPS tstartSSB;
EphemerisData *edat=NULL;          /* Stores earth/sun ephemeris data for barycentering */
LALDetector Detector;              /* Our detector*/
EarthState earth;
EmissionTime emit;
LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
INT4 leap;
static LALStatus status;

int FreeMem();
int GenerateMesh(REAL4VectorSequence **,XYparameterspace *,Metric *);
int SetupPspaceParams(RTparameterspace *,XYparameterspace *);
int GenMetricComp(REAL8 *,REAL8 *,Metric *);
int CheckRTBoundary(REAL8 *,LIGOTimeGPS *,RTparameterspace *);
int ConvertMesh(REAL4VectorSequence **,RTMesh *,RTparameterspace *);
int OutputRTMesh(RTMesh *,Metric *, RTparameterspace *);
int ReadCommandLine(int argc,char *argv[]);
int ConvertTperitoPhase();
int SetupBaryInput();
int CheckInput();
int GetSSBTime(LIGOTimeGPS *, LIGOTimeGPS *);

int main(int argc, char **argv){
  
  XYparameterspace XYspace;
  RTparameterspace RTspace;
  Metric XYMetric;
  RTMesh RTmesh;
  REAL4VectorSequence *XYmesh=NULL;

  if (ReadCommandLine(argc,argv)) return 1;

  if (CheckInput()) return 2;

  if (SetupBaryInput()) return 2;

  if (SetupPspaceParams(&RTspace,&XYspace)) return 2;

  if (GenMetricComp(&(XYspace.X_0),&(XYspace.Y_0),&XYMetric)) return 3;

  if (GenerateMesh(&XYmesh,&XYspace,&XYMetric)) return 4;
  
  if (ConvertMesh(&XYmesh,&RTmesh,&RTspace)) return 5;
    
  if (OutputRTMesh(&RTmesh,&XYMetric,&RTspace)) return 6;
  
  if (FreeMem()) return 7;
  
  LALCheckMemoryLeaks();

  exit(1);

}

/***********************************************************************************/

int FreeMem()
{
 
  /* Here we just free up the memory that has yet to be freed */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  return 0;
 
}

/***********************************************************************************/

int GenerateMesh(REAL4VectorSequence **XYmesh,XYparameterspace *XYspace,Metric *XYMetric)
{
  
  /* this function lays the mesh by calling some of the Flatmesh routines */

  static FlatMeshParamStruc FMparams;
  UINT4Vector *dimlength=NULL;
  REAL4Array *Matrix=NULL;
  REAL4Array *InvMatrix=NULL;
  REAL4VectorSequence *MatrixVec=NULL;
  REAL4VectorSequence *InvMatrixVec=NULL;
  REAL4VectorSequence *Corners=NULL;
  REAL4 *dummy=NULL;
  CreateVectorSequenceIn in;
  
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
  Matrix->data[0]=2.0*((REAL4)sqrt(CLA.mismatch))*XYMetric->eigenvec->data[0]/sqrt(abs((REAL4)DIM*XYMetric->eigenval->data[0]));
  Matrix->data[1]=2.0*((REAL4)sqrt(CLA.mismatch))*XYMetric->eigenvec->data[2]/sqrt(abs((REAL4)DIM*XYMetric->eigenval->data[0]));
  Matrix->data[2]=2.0*((REAL4)sqrt(CLA.mismatch))*XYMetric->eigenvec->data[1]/sqrt(abs((REAL4)DIM*XYMetric->eigenval->data[1]));
  Matrix->data[3]=2.0*((REAL4)sqrt(CLA.mismatch))*XYMetric->eigenvec->data[3]/sqrt(abs((REAL4)DIM*XYMetric->eigenval->data[1]));

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
  LALSDestroyVectorSequence(&status,&MatrixVec);
  LALSDestroyVectorSequence(&status,&InvMatrixVec);

  /* check here if we have actually generated any points in XY space */
  if ((*XYmesh)->length<1) {
    fprintf(stderr,"ERROR : No points have been generated in the XY space\n");
    fprintf(stderr,"        Space may be too small -> Could be a one filter target\n");
    fprintf(stderr,"        Further investigation is required\n");
    exit(1);
  }
  
  return 0;

}

/***********************************************************************************/

int SetupBaryInput()
{

  CHAR filenameE[256],filenameS[256];
  FILE *fp;
  
  /* make the full file name/location for the ephemeris files */
  strcpy(filenameE,CLA.ephemdir);
  strcat(filenameE,"/earth");
  strcat(filenameE,CLA.yr);
  strcat(filenameE,".dat");
  strcpy(filenameS,CLA.ephemdir);
  strcat(filenameS,"/sun");
  strcat(filenameS,CLA.yr);
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
  LALLeapSecs(&status,&leap,&CLA.tstart,&formatAndAcc);
  (*edat).leap=leap;

  /* Read in ephemeris files */
  LALInitBarycenter(&status,edat);             

  /* setup chosen detector */
  if(strcmp(CLA.ifo,"GEO")!=0) Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if(strcmp(CLA.ifo,"LLO")!=0) Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(strcmp(CLA.ifo,"LHO")!=0) Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
 
  /* First job is to convert the observation start time from detector time to SSB time */
  if (GetSSBTime(&CLA.tstart,&tstartSSB)) return 1;

  return 0;

}

/***********************************************************************************/

int GetSSBTime(LIGOTimeGPS *tdet, LIGOTimeGPS *tssb)
{

  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */ 
  EmissionTime emit;

/* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.tgps.gpsSeconds=tdet->gpsSeconds;
  baryinput.tgps.gpsNanoSeconds=tdet->gpsNanoSeconds;
  baryinput.site.location[0]=Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=Detector.location[2]/LAL_C_SI;
  baryinput.alpha=CLA.RA;
  baryinput.delta=CLA.dec;
  baryinput.dInv=0.e0;

  /* Setup barycentric time at start of observation (first timestamp) */
  LALBarycenterEarth(&status,&earth,tdet,edat);
  LALBarycenter(&status,&emit,&baryinput,&earth);
  tssb->gpsSeconds=(emit.te.gpsSeconds);
  tssb->gpsNanoSeconds=emit.te.gpsNanoSeconds;
  
  return 0;
}

/*******************************************************************************/

int SetupPspaceParams(RTparameterspace *RTspace,XYparameterspace *XYspace)
{

  /* this function takes the boundary information given in sma,tperi coordinates   */
  /* and calculates the new overestimate of the boundary in X,Y space.             */
  /* It also converts the central values of sma,tperi into corresponding central   */
  /* values of X and Y. */
  
  LALTimeInterval interval;
  RTPLocation RTPloc;
  XYLocation XYloc;
  REAL8 intervalFLT;
  INT4 NorbCorrection;
  REAL8Vector *Xtemp=NULL;
  REAL8Vector *Ytemp=NULL;
  REAL8Vector *alpha_temp=NULL;
  REAL8 alpha_MAX;
  REAL8 alpha_MIN;
  REAL8 dummy;
  INT4 i;
  UINT4 TWODIM;

  /* check that periapse passage time is within one orbit of the start time */
  /* just a safety thing to make sure that orbital phase errors are calculated */
  /* by hand and correctly (at present) */
  /* we do no error correction here.  We just exit if the periapse time is not within one orbit */
  if (GetSSBTime(&(CLA.tstart),&tstartSSB)) return 1;
  LALDeltaGPS(&status,&interval,&tstartSSB,&(CLA.tperi_0));
  LALIntervalToFloat(&status,&intervalFLT,&interval);
  NorbCorrection=floor(intervalFLT/(CLA.period));
  if (NorbCorrection>0) {
    fprintf(stderr,"ERROR : start time not within one orbit of periapse passage !! \n");
    exit(1);
  }
  
  /* start filling in the RTspace structure */
  RTspace->sma_0=CLA.sma_0;
  RTspace->sma_MIN=CLA.sma_MIN;
  RTspace->sma_MAX=CLA.sma_MAX;
  RTspace->tperi_0.gpsSeconds=CLA.tperi_0.gpsSeconds;
  RTspace->tperi_0.gpsNanoSeconds=CLA.tperi_0.gpsNanoSeconds;
  RTspace->tperi_MIN.gpsSeconds=CLA.tperi_MIN.gpsSeconds;
  RTspace->tperi_MIN.gpsNanoSeconds=CLA.tperi_MIN.gpsNanoSeconds;
  RTspace->tperi_MAX.gpsSeconds=CLA.tperi_MAX.gpsSeconds;
  RTspace->tperi_MAX.gpsNanoSeconds=CLA.tperi_MAX.gpsNanoSeconds;

  /* allocate memory to Xtemp and Ytemp */
  TWODIM=2*(UINT4)DIM;
  LALDCreateVector(&status,&Xtemp,TWODIM);
  LALDCreateVector(&status,&Ytemp,TWODIM);
  LALDCreateVector(&status,&alpha_temp,TWODIM);

  /* fill up RTP location structures */
  RTPloc.period=CLA.period;
  RTPloc.ecc=ECC;   /* this is set to zero in the header */
  RTPloc.argp=ARGP; /* this is set to zero in the header */
  RTPloc.tstartSSB.gpsSeconds=tstartSSB.gpsSeconds;
  RTPloc.tstartSSB.gpsNanoSeconds=tstartSSB.gpsNanoSeconds;
  
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
   
   /* clean up somoe temporarymemory space */
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

  return 0;

}

/***********************************************************************************/

int GenMetricComp(REAL8 *X,REAL8 *Y,Metric *XYMetric)
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
  REAL8 f;
  REAL8 w;
  REAL8 pi;
  INT4 i,j;

  /* set up some more concise variables */
  T=CLA.tspan;
  f=CLA.fmax;
  w=LAL_TWOPI/(CLA.period);
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

  printf("generated metric gamma elements as %lf %lf\n",gMetric->data[2],gMetric->data[4]);
  printf("generated metric gamma elements as %lf %lf\n",gMetric->data[4],gMetric->data[5]);

  /* put the now projected metric into the correct location */
  /* stored normally with redundancy because symmetric */
  XYMetric->element->data[0]=gMetric->data[2];
  XYMetric->element->data[1]=gMetric->data[4];
  XYMetric->element->data[2]=gMetric->data[4];
  XYMetric->element->data[3]=gMetric->data[5];

  printf("generated metric gamma elements as %lf %lf\n",XYMetric->element->data[0],XYMetric->element->data[1]);
  printf("generated metric gamma elements as %lf %lf\n",XYMetric->element->data[2],XYMetric->element->data[3]);

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
  LALDDestroyArray(&status,&elementtemp);

  printf("generated metric gamma elements as %lf %lf\n",XYMetric->element->data[0],XYMetric->element->data[1]);
  printf("generated metric gamma elements as %lf %lf\n",XYMetric->element->data[2],XYMetric->element->data[3]);

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

int ConvertMesh(REAL4VectorSequence **XYmesh,RTMesh *RTmesh,RTparameterspace *RTspace)
{

  /* this section coverts the XY mesh to a RT mesh and retains only those points */
  /* within the original RT boundary. */
  
  RTPLocation RTPloc;
  XYLocation XYloc;
  LIGOTimeGPS *tperi_vec=NULL;
  REAL8Vector *sma_vec=NULL;
  UINT4 i,j;
 
  /* first allocate some memory for the new mesh (using all temporary local variables */
  RTmesh->length=(*XYmesh)->length;
  RTmesh->sma=NULL;
  LALDCreateVector(&status,&(RTmesh->sma),RTmesh->length);

  /* point the temp pointer to the allocated space */
  sma_vec=RTmesh->sma;

  RTmesh->tperi = (LIGOTimeGPS *)LALCalloc( RTmesh->length , sizeof(LIGOTimeGPS) );

  /* point the temp pointer to the alocated space */
  tperi_vec=RTmesh->tperi;
  
  /* fill in the XY location structure */
  XYloc.period=CLA.period;
  XYloc.tstartSSB.gpsSeconds=tstartSSB.gpsSeconds;
  XYloc.tstartSSB.gpsNanoSeconds=tstartSSB.gpsNanoSeconds;
  XYloc.ecc=ECC;
  XYloc.argp=ARGP;

  /* covert XY to RT mesh and trim */
  j=0;
  for (i=0;i<(*XYmesh)->length;i++) {
    XYloc.X=(REAL8)(*XYmesh)->data[i*2];
    XYloc.Y=(REAL8)(*XYmesh)->data[i*2+1];
    /* convert this point in XY space to RT space */
    ConvertXYtoRTperi(&XYloc,&RTPloc);
    /* Now check if this RT point lies within the original boundaries */
    /* if it does save it to memory for later output */
      if (CheckRTBoundary(&RTPloc.sma,&RTPloc.tperi,RTspace)==1) {
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
  else {
    fprintf(stderr,"WARNING : none of the points lie in the original parameter space.\n");
    fprintf(stderr,"          This could be a single filter target but futher investigation\n");
    fprintf(stderr,"          is required (edge effects !!)\n");
  }
  return 0;

}

/***********************************************************************************/

 int OutputRTMesh(RTMesh *RTmesh,Metric *XYMetric, RTparameterspace *RTspace) 
{

  /* this section simply outputs the final mesh to file including the header information */

  FILE *fp;
  BinaryMeshFileHeader BMFheader;
  UINT4 i;

  fp=fopen(CLA.meshfile,"w");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",CLA.meshfile);
  }

  /* setup the header input */
  BMFheader.fmax=CLA.fmax;
  BMFheader.tspan=CLA.tspan;
  BMFheader.tstart.gpsSeconds=CLA.tstart.gpsSeconds;
  BMFheader.tstart.gpsNanoSeconds=CLA.tstart.gpsNanoSeconds;
  BMFheader.Nfilters=RTmesh->length;
  BMFheader.mismatch=CLA.mismatch;
  BMFheader.sma_0=RTspace->sma_0;
  BMFheader.sma_MIN=RTspace->sma_MIN;
  BMFheader.sma_MAX=RTspace->sma_MAX;
  BMFheader.tperi_0.gpsSeconds=RTspace->tperi_0.gpsSeconds;
  BMFheader.tperi_0.gpsNanoSeconds=RTspace->tperi_0.gpsNanoSeconds;
  BMFheader.tperi_MIN.gpsSeconds=RTspace->tperi_MIN.gpsSeconds;
  BMFheader.tperi_MIN.gpsNanoSeconds=RTspace->tperi_MIN.gpsNanoSeconds;
  BMFheader.tperi_MAX.gpsSeconds=RTspace->tperi_MAX.gpsSeconds;
  BMFheader.tperi_MAX.gpsNanoSeconds=RTspace->tperi_MAX.gpsNanoSeconds;
  BMFheader.ecc_MIN=ECC;
  BMFheader.ecc_MAX=ECC;
  BMFheader.argp_MIN=ARGP;
  BMFheader.argp_MAX=ARGP;
  BMFheader.period_MIN=CLA.period;
  BMFheader.period_MAX=CLA.period;
  BMFheader.metric_XX=XYMetric->element->data[0];
  BMFheader.metric_XY=XYMetric->element->data[1];
  BMFheader.metric_YY=XYMetric->element->data[3];
  sprintf(BMFheader.version,"v1");
  sprintf(BMFheader.det,CLA.ifo);
  BMFheader.RA=CLA.RA;
  BMFheader.dec=CLA.dec;

  if (WriteMeshFileHeader(fp,&BMFheader)) return 1;

  /* output the filters in the form ready for a search */
  for (i=0;i<RTmesh->length;i++) {
    fprintf(fp,"%6.12f %6.12f %d %d %6.12f %6.12f\n", \
	    RTmesh->sma->data[i],CLA.period,RTmesh->tperi[i].gpsSeconds, \
	    RTmesh->tperi[i].gpsNanoSeconds,ECC,ARGP);
  }

  fclose(fp);

  return 0;

}

/***********************************************************************************/

 int ReadCommandLine(int argc,char *argv[]) 
{
  INT4 c, errflg = 0;
  CHAR *temp;
  optarg = NULL;
  
  /* Initialize default values */
  CLA.fmax=0.0;
  CLA.tspan=0.0;
  CLA.tstart.gpsSeconds=0;
  CLA.tstart.gpsNanoSeconds=0;
  CLA.RA=0.0;
  CLA.dec=0.0;
  CLA.sma_0=0.0;
  CLA.tperi_0.gpsSeconds=0;
  CLA.tperi_0.gpsNanoSeconds=0;
  CLA.period=0.0;
  CLA.sma_MIN=0.0;
  CLA.sma_MAX=0.0;
  CLA.tperi_MIN.gpsSeconds=0;
  CLA.tperi_MIN.gpsNanoSeconds=0;
  CLA.tperi_MAX.gpsSeconds=0;
  CLA.tperi_MAX.gpsNanoSeconds=0;
  CLA.mismatch=0.0;
  sprintf(CLA.ephemdir,"./");
  sprintf(CLA.yr,"00-04");
  sprintf(CLA.ifo,"LLO");
  sprintf(CLA.meshfile,"mesh.out");

  {
    int option_index = 0;
    static struct option long_options[] = {
      {"fmax", required_argument, 0, 'f'},
      {"tspan", required_argument, 0, 's'},
      {"tstart", required_argument, 0, 'T'},
      {"ra", required_argument, 0, 'a'},
      {"dec", required_argument, 0, 'd'},
      {"sma", required_argument, 0, 'A'},
      {"tperi", required_argument, 0, 'p'},
      {"period", required_argument, 0, 'P'},
      {"smaMIN", required_argument, 0, 'q'},
      {"smaMAX", required_argument, 0, 'Q'},
      {"tperiMIN", required_argument, 0, 'w'},
      {"tperiMAX", required_argument, 0, 'W'},
      {"mismatch", required_argument, 0, 'm'},
      {"ephdir", required_argument, 0, 'E'},
      {"yr", required_argument, 0, 'y'},
      {"ifo", required_argument, 0, 'I'},
      {"meshfile", required_argument, 0, 'o'},
      {"help", no_argument, 0, 'h'}
    };
    /* Scan through list of command line arguments */
    while (!errflg && ((c = getopt_long (argc, argv,"hf:s:T:a:d:A:p:P:q:Q:w:W:m:E:y:I:o:",long_options, &option_index)))!=-1)
      switch (c) {
      case 'f':
	CLA.fmax=atof(optarg);
	break;
      case 's':
	CLA.tspan=atof(optarg);
	break;
      case 'T':
	CLA.tstart.gpsSeconds=atoi(optarg);
	CLA.tstart.gpsNanoSeconds=0;
	break;
      case 'a':
	CLA.RA=atof(optarg);
	break;
      case 'd':
	CLA.dec=atof(optarg);
	break;
      case 'A':
	CLA.sma_0=atof(optarg);
	break;
      case 'p':
	CLA.tperi_0.gpsSeconds=atoi(optarg);
	CLA.tperi_0.gpsNanoSeconds=0;
	break;
      case 'P':
	CLA.period=atof(optarg);
	break;
      case 'q':
	CLA.sma_MIN=atof(optarg);
	break;
      case 'Q':
	CLA.sma_MAX=atof(optarg);
	break;
      case 'w':
	CLA.tperi_MIN.gpsSeconds=atoi(optarg);
	CLA.tperi_MIN.gpsNanoSeconds=0;
	break;
      case 'W':
	CLA.tperi_MAX.gpsSeconds=atoi(optarg);
	CLA.tperi_MAX.gpsNanoSeconds=0;
	break;
      case 'm':
	CLA.mismatch=atof(optarg);
	break;
      case 'E':
	temp=optarg;
	sprintf(CLA.ephemdir,temp);
	break;
      case 'y':
	temp=optarg;
	sprintf(CLA.yr,temp);
	break;
      case 'I':
	temp=optarg;
	sprintf(CLA.ifo,temp);
	break;
      case 'o':
	temp=optarg;
	sprintf(CLA.meshfile,temp);
	break;
      case 'h':
	/* print usage/help message */
	fprintf(stdout,"Arguments are:\n");
	fprintf(stdout,"\t--fmax      REAL8\t Maximum search frequency in Hz [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--tspan     REAL8\t Observation time span in seconds [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--tstart    INT4\t Start time of observation at detector site in GPS seconds [DEFAULT=0]\n");
	fprintf(stdout,"\t--ra        REAL8\t Source right ascension in radians (equatorial) [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--dec       REAL8\t Source declination in radians (equatorial) [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--sma       REAL8\t Central value of projected semi-major axis of orbit in seconds [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--tperi     INT4\t Central value of observed periapse passage measured in the SSB in GPS seconds [DEFAULT=0]\n");
	fprintf(stdout,"\t--period    REAL8\t Orbital period in seconds [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--smaMIN    REAL8\t MINIMUM value of projected semi-major axis in seconds [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--smaMAX    REAL8\t MAXIMUM value of projected semi-major axis in seconds [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--tperiMIN     INT4\t MINIMUM value of observed periapse passage measured in the SSB in GPS seconds [DEFAULT=0]\n");
	fprintf(stdout,"\t--tperiMAX     INT4\t MAXIMUM value of observed periapse passage measured in the SSB in GPS seconds [DEFAULT=0]\n");
	fprintf(stdout,"\t--mismatch  REAL8\t Mismatch required for 2D orbital parameter mesh [DEFAULT=0.0]\n");
	fprintf(stdout,"\t--ephdir    STRING\t Location of ephemeris files earth?.dat and sun?.dat [DEFAULT=NULL]\n");
	fprintf(stdout,"\t--yr        STRING\t Year(s) specifying ephemeris files [DEFAULT=00-04]\n");
	fprintf(stdout,"\t--ifo       STRING\t Interferometer being used for the search (LLO,LHO,GEO,TAMA,CIT,VIRGO) [DEFAULT=LLO]\n");
	fprintf(stdout,"\t--meshfile  STRING\t Name of output mesh file [DEFAULT=mesh.out]\n");
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

int CheckInput()
{

  /* this routine takes the CLA inputs and does some basic validity checks */

  LALDate beginDate;
  LALDate endDate;
  LIGOTimeGPS beginGPS;
  LIGOTimeGPS endGPS;
 
  
  if (CLA.fmax<0.0) {
    fprintf(stderr,"MAX search frequency must be > 0 \n");
    exit(1);
  }
  if (CLA.tspan<0.0) {
    fprintf(stderr,"Observation time span must be > 0 \n");
    exit(1);
  }
  if ((strcmp(CLA.yr,"00")!=0)&&(strcmp(CLA.yr,"01")!=0)&&(strcmp(CLA.yr,"02")!=0)&&(strcmp(CLA.yr,"03")!=0)&&(strcmp(CLA.yr,"00-04")!=0)) {
    fprintf(stderr,"Not a known year for ephemeris file \n");
    exit(1);
  }
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
  if (strcmp(CLA.yr,"00")!=0) {
    beginDate.unixDate.tm_year=100;
    endDate.unixDate.tm_year=101;
  }
  if (strcmp(CLA.yr,"01")!=0) {
    beginDate.unixDate.tm_year=101;
    endDate.unixDate.tm_year=102;
  }
  if (strcmp(CLA.yr,"02")!=0) {
    beginDate.unixDate.tm_year=102;
    endDate.unixDate.tm_year=103;
  }
  if (strcmp(CLA.yr,"03")!=0) {
    beginDate.unixDate.tm_year=103;
    endDate.unixDate.tm_year=104;
  }
  if (strcmp(CLA.yr,"00-04")!=0) {
    beginDate.unixDate.tm_year=100;
    endDate.unixDate.tm_year=105;
  }
  /* convert the beginning and end of the relevant year(s) to a GPS time */
  LALUTCtoGPS(&status,&beginGPS,&beginDate,&(formatAndAcc.accuracy));
  LALUTCtoGPS(&status,&endGPS,&endDate,&(formatAndAcc.accuracy));
  if ((CLA.tstart.gpsSeconds<beginGPS.gpsSeconds)||(CLA.tstart.gpsSeconds+(INT4)CLA.tspan>endGPS.gpsSeconds)) {
    fprintf(stderr,"Start time (+ observation span) must lie within time of ephemeris file\n");
    exit(1);
  }
  if ((CLA.RA<0.0)||(CLA.RA>LAL_TWOPI)) {
    fprintf(stderr,"Source RA must be within range (0 -> 2PI) \n");
    exit(1);
  }
  if ((CLA.dec<(-0.5)*LAL_PI)||(CLA.dec>(0.5)*LAL_PI)) {
    fprintf(stderr,"Source dec must be within range (-PI/2 -> PI/2) \n");
    exit(1);
  }
  if (CLA.sma_0<0.0) {
    fprintf(stderr,"Central value of Orbital semi-major axis must be > 0 \n");
    exit(1);
  }
  if (CLA.period<0.0) {
    fprintf(stderr,"Orbital period must be > 0 \n");
    exit(1);
  }
  if (CLA.mismatch<0.0) {
    fprintf(stderr,"Mismatch must be > 0 \n");
    exit(1);
  }
  if ((strcmp(CLA.ifo,"LLO")!=0)&&(strcmp(CLA.ifo,"LHO")!=0)&&(strcmp(CLA.ifo,"GEO")!=0)) {
    fprintf(stderr,"Not a known detector name \n");
    exit(1);
  }
  if ((CLA.tperi_0.gpsSeconds<beginGPS.gpsSeconds-(INT4)CLA.tspan)||(CLA.tperi_0.gpsSeconds>endGPS.gpsSeconds)) {
    fprintf(stderr,"Central value of periapse passage time (- observation span) must lie within time of ephemeris file\n");
    exit(1);
  }
  if ((CLA.sma_MIN>CLA.sma_0)||(CLA.sma_MAX<CLA.sma_0)) {
    fprintf(stderr,"Central value of the orbital semi-major axis not within MIN and MAX range \n");
    exit(1);
  }
  if ((CLA.tperi_MIN.gpsSeconds>CLA.tperi_0.gpsSeconds)||(CLA.tperi_MAX.gpsSeconds<CLA.tperi_0.gpsSeconds)) {
    fprintf(stderr,"Central value of the periapse passage time not within MIN and MAX range \n");
    exit(1);
  }
  if (CLA.tperi_MIN.gpsSeconds>CLA.tperi_MAX.gpsSeconds) {
    fprintf(stderr,"MIN value of the periapse passage must be < MAX value \n");
    exit(1);
  }

   /* now for a few more physical checks */

  /* check for relativistic speeds in circular orbit */
  if ((CLA.sma_MAX*LAL_TWOPI/CLA.period)>0.01) {
    fprintf(stderr,"WARNING : MAX parameters indicate system is ~ relativistic !! \n");
  }
 
  /* more checks will be added as they are thought of */

  return 0;

}

/*******************************************************************************/
