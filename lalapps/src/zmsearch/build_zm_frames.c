#include<stdio.h>
#include<math.h>
#include<stdarg.h> 
#include <FrameL.h>

#include <lalapps.h>
#include <series.h>
#include <processtable.h>
#include <lalappsfrutils.h>

#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Calibration.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/RealFFT.h>
#include <lal/Interpolate.h>

RCSID( "$Id$" );

#define PLAYC_ENORM  0
#define PLAYC_ESUB   1
#define PLAYC_EARG   2
#define PLAYC_EVAL   3
#define PLAYC_EFILE  4
#define PLAYC_EINPUT 5
#define PLAYC_EMEM   6

#define PLAYC_MSGENORM  "Normal exit"
#define PLAYC_MSGESUB   "Subroutine failed"
#define PLAYC_MSGEARG   "Error parsing arguments"
#define PLAYC_MSGEVAL   "Input argument out of valid range"
#define PLAYC_MSGEFILE  "Could not open file"
#define PLAYC_MSGEINPUT "Error reading file"
#define PLAYC_MSGEMEM   "Out of memory"

/* Usage format string. */       
#define USAGE "Usage: %s [-o outfile] [-i infile]\n"
#define MAX_LENGTH 32768

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                           do                                                                             if (lalDebugLevel & LALERROR ) {                                                                                LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"                              "        %s %s\n", (code), *argv, __FILE__,                                      __LINE__, "test", statement ? statement :                                      "", (msg) );                                                 }                                                                              while(0) 

float **matrix(long nrow, long ncol)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i;
        float **m;

        /* allocate pointers to rows */
        m=(float **) malloc((size_t)((nrow)*sizeof(float*)));
        if (!m) exit(1);

        /* allocate rows and set pointers to them */
        m[0]=(float *) malloc((size_t)((nrow*ncol)*sizeof(float)));
        if (!m[0]) exit(2);

        for(i=1;i<nrow;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

void free_matrix(float **m)
{
  free(m[0]);
  free(m);
}



void swap4(int n,char *here);
int readZM(const char* fname,  float **time,  float **amplitude); 
int readOtt(const char* fname,  float **time,  float **amplitude); 
int snprintf(char *str,size_t size,const char *format, ...);

int main(int argc, char **argv)
{
  static LALStatus stat;                                                      
  INT4 arg=1;                   /* counters                            */

  CHAR *outfile = NULL;           /* name of ascii outfile */
  CHAR *infile  = NULL;          /* name of ilwd outfile */
  CHAR *channel = NULL;          /* name of ilwd outfile */
  FILE *fp;
  int i,j,N[2048],p,l, k, m;
  int useZMWaveforms=1;
  int nfiles=0;
  float n[2048];
  float dt=1.0/16384.0;            /*time diff in ms.=0.2441ms*/
  float unit=1.0;                   /* seconds,  1000.0 would be ms */
  float *time, *amp;
  char filename[30], outputfile[30];   
  float **t, **h;
  float *tmp;
  float x; 
  float **h_intpolat, **t_intpolat;

  float      t_target[5];
  float      h_target[5];
  float      target;
  SInterpolatePar intpar   ={5,t_target,h_target};
  SInterpolateOut intout;

  /* frame output data */
  struct FrFile *frOutFile  = NULL;
  struct FrameH *outFrame   = NULL;

  /* raw input data storage */
  REAL4TimeSeries               zmseries;

  /*******************************************************************
   * PARSE ARGUMENTS (arg stores the current position)               *
   *******************************************************************/

  if (argc <= 1){
    LALPrintError( USAGE, *argv );
    return 0;
  }

  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        outfile = argv[arg++];
      }else{
        ERROR( PLAYC_EARG, PLAYC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PLAYC_EARG;
      }
    }
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        infile = argv[arg++];
      }else{
        ERROR( PLAYC_EARG, PLAYC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PLAYC_EARG;
      }                      
    }
    else if ( !strcmp( argv[arg], "-c" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        channel = argv[arg++];
      }else{
        ERROR( PLAYC_EARG, PLAYC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PLAYC_EARG;
      }                      
    }
    else if ( !strcmp( argv[arg], "-nozm" ) ) {
        arg++;
        useZMWaveforms = 0;
    }
    else if ( !strcmp( argv[arg], "-n" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        nfiles = atoi(argv[arg++]);
      }else{
        ERROR( PLAYC_EARG, PLAYC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PLAYC_EARG;
      }                      
    }
    else if ( !strcmp( argv[arg], "-u" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        unit = atof(argv[arg++]);
        dt = unit * dt;
      }else{
        ERROR( PLAYC_EARG, PLAYC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PLAYC_EARG;
      }                      
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      LALPrintError( USAGE, *argv );
      return PLAYC_EARG;
    }
  } /* End of argument parsing loop. */

  if ( !outfile ){
      LALPrintError( USAGE, *argv );
      return PLAYC_EARG;
  }

  /* allocate some memory */
  t = matrix(nfiles,MAX_LENGTH);
  memset( t[0], 0, nfiles*MAX_LENGTH*sizeof(float) );
  h = matrix(nfiles,MAX_LENGTH);
  memset( h[0], 0, nfiles*MAX_LENGTH*sizeof(float) );


  /******************************************************************
   * Read in the ZM waveforms
   ******************************************************************/
  N[0]=0;
  for(i=0;i<nfiles;i++){

    fprintf(stdout,"Reading in file %s-%d.bin\n", infile, i);
    snprintf(filename,30,"%s-%d.bin", infile, i);

    if ( useZMWaveforms ){
      
      N[i] = readZM(filename,&time,&amp);

      /* print out in same format as Ott et al */
      snprintf(filename,30,"%s-%d.trh", infile, i);
      fp = fopen(filename,"w");
      for(j=0;j<N[i];j++){   
        fprintf(fp, "%e\t%e\t%e\n", unit*(*(time+j)), *(amp+j), 1e-20*(*(amp+j)));
      }
      fclose(fp);
    
    } else {

      N[i] = readOtt(filename,&time,&amp);

    }
    fprintf(stdout,"Read %d lines from file %s-%d.bin\n", N[i],infile, i);

    for(j=0;j<N[i];j++){   
      t[i][j] = unit*(*(time+j));
      h[i][j] = *(amp+j);
    }
    free(time);
    free(amp);
  }

  /*calculate h+ from A*/
  for(i=0;i<nfiles;i++){
    for(j=0;j<N[i];j++){
      /*look into pg.2 of documentation:here angle=90 & R=1cm.*/
      h[i][j] *= (0.273);             
    }
  }


  /************************ find l *****************************************/

  fp = fopen("no.dat","w");
  for(i=0;i<nfiles;i++){
    n[i] = (t[i][N[i]-1]-t[i][0])/dt;
    fprintf(fp,"%f\n",n[i]);
  }
  fclose(fp);

  x=0;
  if(n[1]>=n[0])
    x = n[1];
  else if(n[1]<=n[0])
    x=n[0];
  for(i=2;i<nfiles;i++){
    if(n[i]>=x)
      x = n[i];
    else if(n[i]<=x)
      x = x;
  }
  printf("%f\n",x);

  for(p=0;pow(2,p)<x;++p){
  }

  l=2*pow(2,p);                          /* l is the required arraysize */

  printf("number of points = %d\n",l);
  LALSCreateVector( &stat, &(zmseries.data), l );
  zmseries.deltaT = dt/unit;
  printf("deltaT = %e\n",zmseries.deltaT);
  zmseries.epoch.gpsSeconds = 0;
  zmseries.epoch.gpsNanoSeconds = 0;
  snprintf( zmseries.name, LALNameLength * sizeof(CHAR), "SIM" );


  /******************************interpolate******************************************/

  /* allocate some memory */
  t_intpolat = matrix(nfiles,l);
  h_intpolat = matrix(nfiles,l);

  for(i=0;i<nfiles;i++){
    t_intpolat[i][0]=0;
    h_intpolat[i][0]=0;
  }

  for(i=0;i<nfiles;i++){
    t_intpolat[i][0] = t[i][0];  
    for(j=0;j<(l-1);j++){
      t_intpolat[i][j+1] = t_intpolat[i][j]+dt;
      h_intpolat[i][j+1] = 0;
    }
  }

  for(i=0;i<nfiles;i++){ 
    for(j=0;t_intpolat[i][j]<=t[i][N[i]-1];j++){
      for(k=0;t[i][k]<=t_intpolat[i][j];k++){

      }

      for(m=0;m<4;m++){
        t_target[m]=0;
        h_target[m]=0;
      }

      if(k>=2 && k<=(N[i]-3)){
        for(m=0;m<=4;m++){
          t_target[m]=t[i][k-2+m];
          h_target[m]=h[i][k-2+m];
        }
      }
      else if(k==1){
        for(m=0;m<=4;m++){
          t_target[m]=t[i][k-1+m];
          h_target[m]=h[i][k-1+m];
        }
      }
      else if(k==0){
        for(m=0;m<=4;m++){
          t_target[m]=t[i][k+m];
          h_target[m]=h[i][k+m];
        }
      }
      else if(k==(N[i]-2)){
        for(m=0;m<=4;m++){
          t_target[m]=t[i][k-3+m];
          h_target[m]=h[i][k-3+m];
        }
      }
      else if(k==(N[i]-1)){
        for(m=0;m<=4;m++){
          t_target[m]=t[i][k-4+m];
          h_target[m]=h[i][k-4+m];
        }
      }

      target = t_intpolat[i][j];

      LAL_CALL( LALSPolynomialInterpolation( &stat, &intout, target, &intpar ), 
          &stat);
      h_intpolat[i][j]=intout.y;
    }
  }

  {
    fp = fopen(outfile,"w");
    for(j=0;j<N[1];j++){
      fprintf(fp,"%e %e\n",t_intpolat[1][j], h_intpolat[1][j]);
    }  
    fclose(fp); 
  }

  free_matrix(t);
  free_matrix(h);


  /*************************************************************
   * Build the frames from the interpolated waveforms
   ************************************************************/
  for(i=0;i<nfiles;i++){
    snprintf(outputfile,30,"%s_%0d",channel,i);

    memset( zmseries.data->data, 0, zmseries.data->length*sizeof(REAL4) );
    for(j=0;j<l;j++){
      zmseries.data->data[j]=h_intpolat[i][j];
    }

    outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
        &zmseries, "strain", outputfile );

  }

  /****************************************************************
   * Write out the frame file containing all the different channels
   ****************************************************************/
  frOutFile = FrFileONew( "zm-waveforms.gwf", 0 );
  FrameWrite( outFrame, frOutFile );
  FrFileOEnd( frOutFile );

  return 0;
}


/* function to read in ZM waveforms */
int readZM(const char* fname,  float **time,  float **amplitude)
{
  long n;
  int i;
  FILE *fp;     


  fp = fopen(fname,"r+b");

  fread(&i, sizeof(long), 1, fp);
  fread(&n, sizeof(long), 1, fp);
  fread(&i, sizeof(long), 1, fp);
  fread(&i, sizeof(long), 1, fp);
  swap4(1,(char *)&n);


  (*time) = (float *)malloc( n*sizeof(float));
  fread(*time, sizeof(float), n, fp);
  swap4(n,(char *)(*time));

  (*amplitude) = (float *)malloc( n*sizeof(float));
  fread(*amplitude, sizeof(float), n, fp);
  swap4(n,(char *)(*amplitude));     


  fclose(fp);

  return n;
} 

/* function to read in Ott et al waveforms */
int readOtt(const char* fname,  float **time,  float **amplitude)
{
  int i;
  FILE *fp;     
  CHAR line[1024];
  int nlines=0;
  float tmp,tmpt,tmph;
  float *itime, *iamp;

  fp = fopen(fname,"r");
  /* read data into arrays */
  nlines=0;
  while (1) {
    if (fgets(line,sizeof(line),fp)==NULL) {
      break;
    }
    if (line[0] != '#'){
      sscanf(line,"%e %e %e\n",&tmpt,&tmp,&tmph);
      nlines++;
    }
  }
  fclose(fp);

  itime = (*time) = (float *)malloc( nlines*sizeof(float));
  iamp  = (*amplitude) = (float *)malloc( nlines*sizeof(float));

  fp = fopen(fname,"r");
  /* read data into arrays */
  i=0;
  while (1) {
    if (fgets(line,sizeof(line),fp)==NULL) {
      break;
    }
    if (line[0] != '#'){
      sscanf(line,"%e %e %e\n",&itime[i],&tmp,&iamp[i]);
      i++;
    }
  }
  fclose(fp);

  return nlines;
} 

/* byte swap because the data was written on a machine which is big-endian */
void swap4(int n, char *b1) {
  char *b2,*b3,*b4;
  char temp;
  int i; for (i=0;i<n;i++) {
    b2=b1+1;
    b3=b2+1;
    b4=b3+1; 

    temp=*b1;
    *b1=*b4;
    *b4=temp;
    temp=*b2;
    *b2=*b3;
    *b3=temp;
    b1+=4;
  }
  return;
}








