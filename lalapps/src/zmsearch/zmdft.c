#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
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
#include <lal/Units.h>
#include <lal/FrameStream.h>
#include <lal/ReadNoiseSpectrum.h>

RCSID("$Id$");

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "zmdft"

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
#define USAGE \
"Usage: %s [options]\n\n"\
"  --outfile outfile            Output file for something\n"\
"  --dirname dir                Name of directory with zm-waveforms.gwf\n"\
"  --min-match mu               Minimal match for orthonormalization\n"\
  "\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )  \
do if (lalDebugLevel & LALERROR )      \
{                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
      "%s %s\n", (code), *argv, __FILE__,  __LINE__, \
      PLAYC, statement ? statement :"", (msg) ); \
}while(0)

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

        for(i=1;i<=nrow;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}



int snprintf(char *str,size_t size,const char *format, ...); 
int zmnormalise(int n,int n1,float **amp,float rhosq,float **eamplitude,REAL4FrequencySeries *spectrum);
void zmfft(int n1, float *ampfft, float *ampfftre, float *ampfftim);
void correl(int num, float *re, float *RE, float *im, float *IM, float *corr);
float zminproduct(int num, float *re, float *im, float *amp);

/* ZM waveform parameters */
FrChanIn   channelIn;               /* channnel information                  */
INT4       numPoints      = 4096;    /* number of samples from frames: 0.25s */
INT4       sampleRate     = 16384;   /* sample rate in Hz                    */
CHAR      *dirname        = NULL;    /* name of directory with frame file    */
CHAR      *outFile        = NULL;    /* name of ascii outfile                */

/* some global output flags */
int verbose = 0;

int main(int argc, char **argv)
{
  static LALStatus      stat;
  FrStream             *stream        = NULL;
  FILE *fp;
  int i, j;
  int n, m;
  char filename[30], outputfile[30];
  float Amp[78][4096];
  float *amp[78];
  float *eamplitude[78];
  float minimal_match=0.0;

  /* data storage */
  REAL4TimeSeries            series;
  REAL4FrequencySeries       spectrum;

  /* frame output data */
  struct FrFile *frOutFile  = NULL;
  struct FrameH *outFrame   = NULL;


  /*******************************************************************
   * PARSE ARGUMENTS 
   *******************************************************************/

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,       &verbose,          1 },
    {"dirname",                 required_argument, 0,                'd'},
    {"outfile",                 required_argument, 0,                'o'},
    {"min-match",               required_argument, 0,                'm'},
    {"help",                    no_argument,       0,                'h'}, 
    {"debug-level",             required_argument, 0,                'z'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };
  int switchc;

  /*
   * 
   * initialize things
   *
   */

  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );
  setvbuf( stdout, NULL, _IONBF, 0 );

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    int optarg_len = 0;

    switchc = getopt_long_only( argc, argv, 
        "o:hz:V", long_options, 
	&option_index );

    /* detect the end of the options */
    if ( switchc == -1 )
    {
      break;
    }

    switch ( switchc )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "Error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'd':
        optarg_len = strlen( optarg ) + 1;
        dirname = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( dirname, optarg, optarg_len );
        break;

      case 'o':
        optarg_len = strlen( optarg ) + 1;
        outFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( outFile, optarg, optarg_len );
        break;

      case 'm':
        minimal_match = atof( optarg );
        break;

      case 'h':
        /* help message */
        fprintf( stderr, USAGE , argv[0]);
        exit( 1 );
        break;

      case 'z':
        set_debug_level( optarg );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Compute dimension of vector space spanned by ZM\n" 
            "Saikat Ray-Majumdar and Patrick Brady\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
        break;

      case '?':
        fprintf( stderr, USAGE , argv[0]);
        exit( 1 );
        break;

      default:
        fprintf( stderr, "Error: Unknown error while parsing options\n" );
        fprintf( stderr, USAGE, argv[0] );
        exit( 1 );
    }
  }

  /* check that the frame file can be located */
  if ( !dirname ){
    fprintf(stderr,"Must supply directory name with zm-waveforms.gwf\n");
    fprintf( stderr, USAGE, argv[0] );
    exit( 1 );
  }
  
  /* check that minimal match is supplied  */
  if ( minimal_match <= 0.0 ){
    fprintf(stderr,"Must supply minimal match\n");
    fprintf( stderr, USAGE, argv[0] );
    exit( 1 );
  }
  /* End of argument parsing loop. */


  /*******************************************************
   * step 1:  read in the ZM waveforms from the frame file
   ******************************************************/

  /* create and initialize the time series vector */
  series.data = NULL;
  LAL_CALL( LALCreateVector( &stat, &series.data, numPoints), &stat);
  memset( series.data->data, 0, series.data->length*sizeof(REAL4) );
  series.epoch.gpsSeconds     = 0;
  series.epoch.gpsNanoSeconds = 0;
  series.deltaT = 1.0/((REAL8) sampleRate);
  series.f0 = 0.0;
  series.sampleUnits = lalADCCountUnit;

  /* get the frame stream */
  LAL_CALL( LALFrOpen( &stat, &stream, dirname, "zm-waveforms.gwf" ), &stat);

  /* initialize everything */
  strcpy(series.name,"ZMSIM_ZM_WAVEFORM_0" );
  channelIn.name = series.name;
  channelIn.type = LAL_PROC_CHAN;
  LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream), &stat);
  series.epoch.gpsSeconds     = 0;
  series.epoch.gpsNanoSeconds = 0;

  /* loop over the frame file */
  for(i=0;i<78;i++){

    /* make sure we're at the start of the frame file */
    LAL_CALL( LALFrSeek(&stat, &(series.epoch), stream), &stat);

    /* which waveform are we going to select? */
    snprintf(filename,30,"ZMSIM_ZM_WAVEFORM_%d",i);
    strcpy(series.name, filename);
    channelIn.name = series.name;

    /* get the data */
    LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream), 
        &stat);

    /* put the data into the array that we have for it */
    for(j=0;j<series.data->length;j++){
      /* amp. of the interpolated waveforms*/
      Amp[i][j] = series.data->data[j];  
    }
  }

  if ( verbose )
    fprintf(stdout,"Finished reading in the ZM waveforms\n");


  /****************************************************************
   * Read in the noise power spectrum
   ***************************************************************/
  spectrum.f0 = 0;
  spectrum.deltaF = (REAL4)(sampleRate) / (REAL4)(numPoints);
  spectrum.data = NULL;

  LAL_CALL( LALCreateVector( &stat, &(spectrum.data), numPoints/2+1 ), 
      &stat );

  LAL_CALL( LALReadNoiseSpectrum( &stat, &spectrum, "noise-initial_ligo.dat"), &stat);

  if ( verbose )
    fprintf(stdout,"Finished reading in the noise spectrum\n");



  /**************************** 
   * step 2 
   ****************************/

  for(i=0;i<78;i++){
    amp[i] = Amp[i];
  }

  n=78;
  m=0;

  /** m = no. of basis vectors required**/
  m=zmnormalise(n,numPoints,amp,minimal_match,eamplitude,&spectrum); 
  printf("%f %d\n",minimal_match,m);

  /*************************************************************
   * Build the frames from the interpolated waveforms
   ************************************************************/
  for(i=1;i<m;i++){
    snprintf(outputfile,30,"ZMBASIS_%0d",i);

    for(j=0;j<numPoints;j++){
      series.data->data[j]=eamplitude[i][j];
    }

    outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
        &series, "strain", outputfile );

  }
  LAL_CALL( LALSDestroyVector( &stat, &(series.data) ), &stat);

  /****************************************************************
   * Write out the frame file containing all the different channels
   ****************************************************************/
  frOutFile = FrFileONew( "zm-basis.gwf", 0 );
  FrameWrite( outFrame, frOutFile );
  FrFileOEnd( frOutFile );

  return 0;   
}






/*****************************************************
 *****************************************************/

#define NMAX 100
#define N1MAX 4096

int zmnormalise(int n,int n1,float **amp,float rhosq,float **eamplitude,REAL4FrequencySeries *spectrum)
{
  int i,j,p,l,k,m;
  float *ampfft=NULL, *ampfftre=NULL, *ampfftim=NULL;
  float sum1, sum_ii;
  float **ampre, **ampim;
  float **correlation, **eamp;
  float rhom[N1MAX];
  float *redum[NMAX], *imdum[NMAX], *ampdum[NMAX];
  float e_real[NMAX][N1MAX/2+1], e_imagin[NMAX][N1MAX/2+1];
  float c[NMAX];
  float *rho_max, *shift, *corr;
  float *RE, *IM, *re, *im;
  float *ampshift, *o_re, *o_im;
  float *oamp, *oamp_re, *oamp_im;
  float *dum1, *dum2, *dum3;
  int index,toshift;
  float sumfinal, rho_max_min;

  if ( n > NMAX || n1 > N1MAX )
  {
    fprintf( stderr, "Error:\n" );
    fprintf( stderr, "n  = %d, nmax  = %d\n", n, NMAX );
    fprintf( stderr, "n1 = %d, n1max = %d\n", n1, N1MAX );
    exit( 1 );
  }

  ampre=matrix(NMAX,N1MAX/2+1);
  ampim=matrix(NMAX,N1MAX/2+1);

  /************************
   * Step 3:normalise the waveforms
   * **********************/
  ampfftre=(float*)malloc((n1/2+1)*sizeof(float));
  ampfftim=(float*)malloc((n1/2+1)*sizeof(float));

  for(i=0;i<n;i++){


    ampfft = amp[i];

    sum1=0;

    zmfft(n1,ampfft,ampfftre,ampfftim);  /** this function finds the fft **/

    for(j=0;j<=n1/2;j++){
      spectrum->data->data[j] /= spectrum->data->data[0];
    }
    
    for(j=0;j<=n1/2;j++){
      sum1 += 4.0*(ampfftre[j]*ampfftre[j]+ampfftim[j]*ampfftim[j]) / 
        (spectrum->data->data[j]);
    }
    
    for(j=0;j<n1;j++){
      amp[i][j]=amp[i][j]/sqrt(sum1);
    }
    
    for(j=0;j<=n1/2;j++){
      ampre[i][j]=(ampfftre[j]/sqrt(sum1))/(sqrt(spectrum->data->data[j]));
      ampim[i][j]=(ampfftim[j]/sqrt(sum1))/(sqrt(spectrum->data->data[j]));
    }

  }
  free(ampfftre);
  free(ampfftim);

  /*********************
   * Step 4: assign the dummy pointers 
   *******************************************/

  for(i=0;i<n;i++){
    ampdum[i]=amp[i];
    redum[i]=ampre[i];
    imdum[i]=ampim[i];
  }

  /*********************
   * Step 5: Set e_0
   ******************************************/
  eamp=matrix(NMAX,N1MAX);
  for(j=0;j<n1;j++){
    eamp[0][j]=ampdum[0][j];
  }
  for(j=0;j<=n1/2;j++){
    e_real[0][j]=redum[0][j];
    e_imagin[0][j]=imdum[0][j];
  }

  /*******************
   * Step 6: Orthonormalisation
   * ******************************************/
  correlation=matrix(NMAX,N1MAX);

  for(i=1;i<(n);i++){

    /********************
     * Step 6a
     * ******************************************/
    rho_max=(float*)malloc(n*sizeof(float));
    shift=(float*)malloc(n*sizeof(float));

    for(p=0;p<n;p++){
      rho_max[p]=0;
      shift[p]=0;
    }
    /*******************
     * Step 6b
     * ***************************************/
    for(j=i;j<n;j++){

      for(l=0;l<n1;l++){               /*rhom:rho^2; set rho^2 to 0*/
        rhom[l]=0;
      }

      RE = redum[j];
      IM = imdum[j];

      /********
       * Step 6b-1: find the correlation of the waveform pointed 
       * to by the j-th dummy pointer with the basis vectors already 
       * formed
       ****************************/

      for(k=i-1;k>=0;k--){
        re = e_real[k];
        im = e_imagin[k];

        corr=(float*)malloc(n1*sizeof(float));

        correl(n1,re,RE,im,IM,corr);  /********function*********/

        for(l=0;l<n1;l++){
          correlation[k][l]=0;
        }

        for(l=0;l<n1;l++){   
          correlation[k][l]=corr[l];

        }

        free(corr);
      }

      /*********Step 6b-2******************************/

      for(l=0;l<n1;l++){
        for(k=i-1;k>=0;k--){
          rhom[l]+=(correlation[k][l]*correlation[k][l]);
        }
      }

      /**********Step 6b-3*****************************/

      if(rhom[0]>=rhom[1]){
        rho_max[j]=rhom[0];
        shift[j]=0;
      }
      else if(rhom[0]<rhom[1]){
        rho_max[j]=rhom[1];
        shift[j]=1;
      }
      for(l=2;l<n1;l++){
        if(rho_max[j]>=rhom[l]){
          rho_max[j]=rho_max[j];
          shift[j]=shift[j];
        }
        else if(rho_max[j]<rhom[l]){
          rho_max[j]=rhom[l];
          shift[j]=l;
        }
      }
    }

    /*************Step 6c********************************/

    index=0;
    toshift=0;
    rho_max_min=0;

    if(rho_max[i]>=rho_max[i+1]){
      rho_max_min=rho_max[i+1];
      index=(i+1);
      toshift=shift[i+1];
    }
    else if(rho_max[i]<rho_max[i+1]){
      rho_max_min=rho_max[i];
      index=i;
      toshift=shift[i];
    }
    for(k=i+2;k<n;k++){
      if(rho_max_min>=rho_max[k]){
        rho_max_min=rho_max[k];
        index=k;
        toshift=shift[k];
      }
      else if(rho_max_min<rho_max[k]){
        rho_max_min=rho_max_min;
        index=index;
        toshift=toshift;
      }
    }
    /*******************************************/

    /*	printf("index=%d toshift=%d\n",index, toshift);
        printf("corr_min=%f\n",4.0*rho_max_min);*/


    /**************Step 6d: Check******************************/

    if(rhosq<=rho_max_min && 1.0>=rho_max_min)
      break;

    /**************Step 6e*****************************/

    ampshift=(float*)malloc(n1*sizeof(float));	
    if(toshift<=n1/2){
      for(j=toshift;j<n1;j++){
        ampshift[j]=ampdum[index][j-toshift];
      }
      for(j=0;j<toshift;j++){
        ampshift[j]=0;
      }
    }
    else if(toshift>n1/2){
      toshift=n1-toshift;
      for(j=0;j<(n1-toshift);j++){
        ampshift[j]=ampdum[index][j+toshift];
      }
      for(j=(n1-toshift);j<n1;j++){
        ampshift[j]=0;
      }
    }

    /************Step 6f**************************************/

    for(k=i-1;k>=0;k--){

      o_re=e_real[k];
      o_im=e_imagin[k];

      c[k]=zminproduct(n1,o_re,o_im,ampshift);  /*******function********/
    }

    /************Step 6g: Orthogonalize***************************************/

    oamp=(float*)malloc(n1*sizeof(float));

    for(j=0;j<n1;j++){
      sumfinal =0;
      for(k=i-1;k>=0;k--){
        sumfinal+=(c[k]*eamp[k][j]);
      }
      oamp[j]=ampshift[j]-sumfinal;
    }

    /************Step 6h***********************************/
    oamp_re=(float*)malloc((n1/2+1)*sizeof(float));
    oamp_im=(float*)malloc((n1/2+1)*sizeof(float));


    zmfft(n1,oamp,oamp_re,oamp_im);            /******function******/

    sum_ii=0;	  
    for(j=0;j<=n1/2;j++){
      sum_ii+=4.0*(oamp_re[j]*oamp_re[j]+oamp_im[j]*oamp_im[j]);
    }


    for(j=0;j<n1;j++){
      eamp[i][j]=oamp[j]/sqrt(sum_ii);
    }
    /*  for(j=0;j<n1;j++){	
        eamp[i][j]=oamp[j];
        }*/
    eamplitude[i]=eamp[i];


    for(j=0;j<=n1/2;j++){
      e_real[i][j]=oamp_re[j]/sqrt(sum_ii);
      e_imagin[i][j]=oamp_im[j]/sqrt(sum_ii);
    }
    /*  for(j=0;j<=n1/2;j++){
        e_real[i][j]=oamp_re[j];
        e_imagin[i][j]=oamp_im[j];
        }*/

    /***********************Step 6i***************************************/

    dum1=redum[i];
    redum[i]=redum[index];
    redum[index]=dum1;
    dum2=imdum[i];
    imdum[i]=imdum[index];
    imdum[index]=dum2;
    dum3=ampdum[i];				  
    ampdum[i]=ampdum[index];
    ampdum[index]=dum3;

    /*************************************************************/

    free(rho_max);
    free(shift);
    free(oamp);
    free(oamp_re);
    free(oamp_im);
    free(ampshift);


    /************************************************************/

  }
  /************************return the no. of basis vectors**************/
  m=i-1;

  return m;
}

#undef NMAX
#undef N1MAX







void zmfft(int n1, float *ampfft, float *ampfftre, float *ampfftim)
{
  static LALStatus stat;
  int i;

  RealFFTPlan *pfwd = NULL;   /* FFTW uses a plan to assess best FFT method*/    
  REAL4Vector *hvec = NULL;        /* Fundamental LAL data types */
  COMPLEX8Vector *Hvec = NULL;     /* see package "std" code and */


  LALCreateForwardRealFFTPlan( &stat, &pfwd, n1, 0);

  /* Create an S (float) vector of length "n" to hold time data */
  LALSCreateVector( &stat, &hvec, n1 );

  /* Create C (complex) vector of length n/2+1 to hold FFT */
  LALCCreateVector( &stat, &Hvec, n1/2 + 1 );


  for(i=0;i<n1;i++){
    hvec->data[i]=ampfft[i];
  }

  /* do a forward FFT */  
  LALForwardRealFFT( &stat, Hvec, hvec, pfwd );

  for ( i=0 ; i<=n1/2; i++){
    REAL4 re1=Hvec->data[i].re;
    REAL4 im1=Hvec->data[i].im;
    ampfftre[i] = re1;
    ampfftim[i] = im1;
  }


  LALDestroyRealFFTPlan( &stat, &pfwd );

  LALSDestroyVector( &stat, &hvec );
  LALCDestroyVector( &stat, &Hvec );


  return ;

}


void correl(int num, float *re, float *RE, float *im, float *IM, float *corr)
{

  static LALStatus stat;
  int i;
  float *productre,*productim;



  RealFFTPlan *prev = NULL;   /* This one is for the reverse FFT */
  REAL4Vector *hvec = NULL;   /* Fundamental LAL data types */
  COMPLEX8Vector *Hvec = NULL;     /* see package "std" code and */



  /* Create an FFTW plan for reverse REAL FFT */
  LALCreateReverseRealFFTPlan( &stat, &prev, num, 0);

  /* Create an S (float) vectors of length "num" to hold time data */
  LALSCreateVector( &stat, &hvec, num );

  /* Create C (complex) vectors of length num/2+1 to hold FFT */
  LALCCreateVector( &stat, &Hvec, num/2+1 );


  productre =(float*) malloc((num/2+1)*sizeof(float));
  productim = (float*)malloc((num/2+1)*sizeof(float));
  /* productre[num/2]=0;
     productim[num/2]=0;*/
  for(i=0;i<=num/2;i++){
    productre[i] = re[i]*RE[i] + im[i]*IM[i];
    productim[i] = -IM[i]*re[i] + RE[i]*im[i];
  }

  for(i=0;i<=num/2;i++){
    REAL4 pre = productre[i];
    REAL4 pim = productim[i];
    Hvec->data[i].re = pre;
    Hvec->data[i].im = pim;
  }

  LALReverseRealFFT(&stat, hvec, Hvec, prev);

  for(i=0;i<num;i++){
    corr[i]=2.0*hvec->data[i];
  }

  LALDestroyRealFFTPlan( &stat, &prev );

  /* get rid of the vectors */
  LALSDestroyVector( &stat, &hvec );
  LALCDestroyVector( &stat, &Hvec );

  free(productre);
  free(productim);

  return ;
}


float zminproduct(int num, float *re, float *im, float *amp)
{

  static LALStatus stat;
  int i;
  float *reA,*imA;
  float sum;



  RealFFTPlan *pfwd = NULL;   /* FFTW uses a plan to assess best FFT method*/
  REAL4Vector *hvec = NULL;        /* Fundamental LAL data types */
  COMPLEX8Vector *Hvec = NULL;     /* see package "std" code and */


  /* Create an FFTW plan for forward REAL FFT */
  LALCreateForwardRealFFTPlan( &stat, &pfwd, num, 0);


  /* Create an S (float) vestor of length "n" to hold time data */
  LALSCreateVector( &stat, &hvec, num );


  /* Create C (complex) vector of length n/2+1 to hold FFT */
  LALCCreateVector( &stat, &Hvec, num/2 + 1 );


  for(i=0;i<num;i++){
    hvec->data[i]=amp[i];
  }

  /* do a forward FFT */
  LALForwardRealFFT( &stat, Hvec, hvec, pfwd );

  reA=malloc((num/2+1)*sizeof(float));
  imA=malloc((num/2+1)*sizeof(float));


  for ( i=0 ; i<=num/2; i++){
    REAL4 re1=Hvec->data[i].re;
    REAL4 im1=Hvec->data[i].im;
    reA[i] = re1;
    imA[i] = im1;
  }


  sum=0;

  for(i=0;i<=num/2;i++){
    sum += reA[i]*re[i] + imA[i]*im[i];
  }


  LALDestroyRealFFTPlan( &stat, &pfwd );


  /* get rid of the vectors */
  LALSDestroyVector( &stat, &hvec );
  LALCDestroyVector( &stat, &Hvec );

  free(reA);
  free(imA);

  return 4.0*sum;
}

