#include<stdio.h>
#include<math.h>
#include<stdarg.h> 
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Calibration.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/RealFFT.h>
#include <lal/Interpolate.h>

NRCSID( PLAYC, "$Id$" );

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

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                   do if (lalDebugLevel & LALERROR )                                                 {                                                                                LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"                              "        %s %s\n", (code), *argv, __FILE__,                                      __LINE__, PLAYC, statement ? statement :                                      "", (msg) );                                                 }      while(0)

INT4 lalDebugLevel = LALMSGLVL3;

int snprintf(char *str,size_t size,const char *format, ...); 
int zmnormalise(int n,int n1,float **amp,float rhosq,float **eamplitude);
void zmfft(int n1, float *ampfft, float *ampfftre, float *ampfftim);
void correl(int num, float *re, float *RE, float *im, float *IM, float *corr);
float zminproduct(int num, float *re, float *im, float *amp);

int main(int argc, char **argv)
{
                                                        
  INT4 arg=1;                   /* counters                            */
  
  CHAR *outfile = NULL;           /* name of ascii outfile */
  CHAR *infile = NULL;          /* name of ilwd outfile */
  FILE *fp;
  int i, j;
  int n1;
  int rc, n, m;
  float *c, *d;
  float *cc, *dd;
  char filename[30], outputfile[30];
  float time[78][16385];
  float Amp[78][16385];
  float *amp[78];
  float *eamplitude[78];
  float rhosq;
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
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      LALPrintError( USAGE, *argv );
      return PLAYC_EARG;
    }
  } /* End of argument parsing loop. */

  /******************************* step 1 ****************************************************/

 for(i=0;i<78;i++){
   
    snprintf(filename,30,"zminpol-%d.dat",i);
  
    n1=0;
  
  fp = fopen(filename,"r");
 
  while (fscanf(fp,"%*f %*f")!= EOF)                     /* read in the interpolated waveforms; n2 = 16384 */
                  ++n1;
  
  c = cc = (float*)malloc(n1*sizeof(float));
  d = dd = (float*)malloc(n1*sizeof(float));
   
  rewind (fp);
  while((rc = fscanf(fp,"%f %f", cc++, dd++)) == 2) continue;

  if(rc !=EOF) 
    printf("error reading file");

  fclose(fp);


      for(j=0;j<n1;j++){
	time[i][j] = c[j];
	Amp[i][j] = d[j];  /*time and amp. of the interpolated waveforms*/
    }
  free(c);
  free(d);
   
 }
 

 /**************************** step 2 ******************************************************/

 for(i=0;i<78;i++){
   amp[i] = Amp[i];
 }

 n=78;
 m=0;
 rhosq=0.6;

 
 m=zmnormalise(n,n1,amp,rhosq,eamplitude); /** m = no. of basis vectors required**/
 printf("%f %d\n",rhosq,m);
   

 /****************************step 7:Write the basis vectors to output file******************/


   for(i=1;i<=m;i++){
 snprintf(outputfile,30,"e_amplitude-%d.dat",i);

 fp = fopen(outputfile,"w");
 
 for(j=0;j<n1;j++){
   fprintf(fp,"%f %f\n",time[i][j],eamplitude[i][j]);
 }

 fclose(fp);
 }  



 return 0;   
}
 





/******************************************************************************************
 ******************************************************************************************/

enum { nmax = 100, n1max = 100 };

int zmnormalise(int n,int n1,float **amp,float rhosq,float **eamplitude)
{
  int i,j,p,l,k,m;
  float *ampfft, *ampfftre, *ampfftim;
  float sum1, sum_ii;
  float ampre[nmax][n1max/2+1], ampim[nmax][n1max/2+1];
  float correlation[nmax][n1max], eamp[nmax][n1max];
  float rhom[n1max];
  float *redum[nmax], *imdum[nmax], *ampdum[nmax];
  float e_real[nmax][n1max/2+1], e_imagin[nmax][n1max/2+1];
  float c[nmax];
  float *rho_max, *shift, *corr;
  float *RE, *IM, *re, *im;
  float *ampshift, *o_re, *o_im;
  float *oamp, *oamp_re, *oamp_im;
  float *dum1, *dum2, *dum3;
  int index,toshift;
  float sumfinal, rho_max_min;

  if ( n > nmax || n1 > n1max )
  {
    fprintf( stderr, "Error:\n" );
    fprintf( stderr, "n  = %d, nmax  = %d\n", n, nmax );
    fprintf( stderr, "n1 = %d, n1max = %d\n", n1, n1max );
    exit( 1 );
  }




  /************Step 3:normalise the waveforms**********************/
  for(i=0;i<n;i++){
  
  
      ampfft = amp[i];
 
  ampfftre=(float*)malloc((n1/2+1)*sizeof(float));
  ampfftim=(float*)malloc((n1/2+1)*sizeof(float));

  sum1=0;

  zmfft(n1,ampfft,ampfftre,ampfftim);  /** this function finds the fft **/

  for(j=0;j<=n1/2;j++){
    sum1 += 4.0*(ampfftre[j]*ampfftre[j]+ampfftim[j]*ampfftim[j]);
      }
  for(j=0;j<n1;j++){
    amp[i][j]=amp[i][j]/sqrt(sum1);
  }
  for(j=0;j<=n1/2;j++){
    ampre[i][j]=ampfftre[j]/sqrt(sum1);
    ampim[i][j]=ampfftim[j]/sqrt(sum1);
  }
 
  free(ampfftre);
  free(ampfftim);

  }

  /*********************Step 4: assign the dummy pointers ******************************************/

  for(i=0;i<n;i++){
    ampdum[i]=amp[i];
    redum[i]=ampre[i];
    imdum[i]=ampim[i];
  }

  /*********************Step 5: Set e_0*****************************************/
  for(j=0;j<n1;j++){
  eamp[0][j]=ampdum[0][j];
  }
  for(j=0;j<=n1/2;j++){
  e_real[0][j]=redum[0][j];
  e_imagin[0][j]=imdum[0][j];
  }

  /*******************Step 6: Orthonormalisation******************************************/


  for(i=1;i<(n-5);i++){
  
       /********************Step 6a******************************************/
	rho_max=(float*)malloc(n*sizeof(float));
	shift=(float*)malloc(n*sizeof(float));

	for(p=0;p<n;p++){
	  rho_max[p]=0;
	  shift[p]=0;
	}
	/*******************Step 6b***************************************/
	for(j=i;j<n;j++){

	  for(l=0;l<n1;l++){               /*rhom:rho^2; set rho^2 to 0*/
	    rhom[l]=0;
	  }

	    RE = redum[j];
	    IM = imdum[j];

	              /********Step 6b-1: find the correlation of the waveform pointed to by the 
             j-th dummy pointer with the basis vectors already formed***************************/

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

