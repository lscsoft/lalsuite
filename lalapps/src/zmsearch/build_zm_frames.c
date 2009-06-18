/*
*  Copyright (C) 2007 Patrick Brady, Saikat Ray-Majumder
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
#define USAGE "Usage: %s --infileindex waveform id --channel channel suffix " \
        "--nfiles #files --frame-length framelength(sec:Must be power of 2) " \
        "[--verbose] [--printrawdata] [--printinterpolateddata] [--help]\n"

#define MAX_LENGTH 81920
#define TRUE  1
#define FALSE 0


static float **matrix(long nrow, long ncol)
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

static void free_matrix(float **m)
{
  free(m[0]);
  free(m);
}

/* byte swap because the data was written on a machine which is big-endian */
static void swap4(int n, char *b1) {
  char *b2,*b3,*b4;
  char temp;
  int i; 

  for (i=0;i<n;i++) {
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

/* function to read in ZM waveforms */
static int readZM(const char* fname,  float **time,  float **amplitude, float **ampcross)
{
  long n;
  int i;
  FILE *fp=NULL;     


  fp = fopen(fname,"r+b");

  if(!fp){
    fprintf(stderr,"File %s doesn't exit\n",fname);
    exit(6);
  }

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

  (*ampcross) = (float *)malloc( n*sizeof(float));
  memset( *ampcross, 0, n*sizeof(float) );

  fclose(fp);

  return n;
} 

/* function to read in Ott et al waveforms */
static int readOtt(const char* fname,  float **time,  float **amplitude, float **ampcross)
{
  int i;
  FILE *fp;     
  CHAR line[1024];
  int nlines=0;
  float tmp,tmpt,tmph;
  float *itime, *iamp;

  fp = fopen(fname,"r");

  if(!fp){
    fprintf(stderr,"File %s doesn't exit\n",fname);
    exit(6);
  }

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

  (*ampcross) = (float *)malloc( nlines*sizeof(float));
  memset( *ampcross, 0, nlines*sizeof(float) );

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

static int readWarren(const char* fname,  float **time,  float **amplitude, float **ampcross)
{
  int i;
  FILE *fp;     
  CHAR line[1024];
  int nlines=0;
  float tmp,tmpt,tmph;
  float *itime, *iamp, *iampcross;

  fp = fopen(fname,"r");

  if(!fp){
    fprintf(stderr,"File %s doesn't exist\n",fname);
    exit(6);
  }

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
  iampcross  = (*ampcross) = (float *)malloc( nlines*sizeof(float));

  fp = fopen(fname,"r");
  /* read data into arrays */
  i=0;
  while (1) {
    if (fgets(line,sizeof(line),fp)==NULL) {
      break;
    }
    if (line[0] != '#'){
      sscanf(line,"%e %e %e\n",&itime[i],&iamp[i],&iampcross[i]);
      i++;
    }
  }
  fclose(fp);

  return nlines;
} 

static void build_interpolate ( float **h, float **t, float **h_intpolat, float **t_intpolat, float dt, int l, int *N, int nfiles ){

  static LALStatus stat;

  float      t_target[5];
  float      h_target[5];
  float      target;
  SInterpolatePar intpar   ={5,t_target,h_target};
  SInterpolateOut intout;

  INT4 i,j,k,m;

  for(i=0;i<nfiles;i++){
    for(j=0;j<l;j++){
      t_intpolat[i][j] = j*dt;
      h_intpolat[i][j] = 0;
    }
  }

  for(i=0;i<nfiles;i++){
    for(j=0;t_intpolat[i][j]<t[i][N[i]-1];j++){
      for(k=0;t[i][k]<=t_intpolat[i][j];k++){
      }

      for(m=0;m<5;m++){
        t_target[m]=0;
        h_target[m]=0;
      }

      if(k>=2 && k<=(N[i]-3)){
        for(m=0;m<5;m++){
          t_target[m]=t[i][k-2+m];
          h_target[m]=h[i][k-2+m];
        }
      }
      else if(k==1){
        for(m=0;m<5;m++){
          t_target[m]=t[i][k-1+m];
          h_target[m]=h[i][k-1+m];
        }
      }
      else if(k==0){
        for(m=0;m<5;m++){
          t_target[m]=t[i][k+m];
          h_target[m]=h[i][k+m];
        }
      }
      else if(k==(N[i]-2)){
        for(m=0;m<5;m++){
          t_target[m]=t[i][k-3+m];
          h_target[m]=h[i][k-3+m];
        }
      }
      else if(k==(N[i]-1)){
        for(m=0;m<5;m++){
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
}
 
static void find_tpeak( float **hplus_intpolat, float **hcross_intpolat, float dt, INT4 nfiles, INT4 l, int *tpeak){

  int i,j;

  float Fplus, Fcross;
  float **hsquare, **t;
  float htmp;
  int ttmp;

  t = matrix(nfiles,l);
  hsquare = matrix(nfiles,l);

  Fplus = 1.0;
  Fcross = 0.0;

  for(i=0; i<nfiles;i++){
    for(j=0;j<l;j++){
      t[i][j] = j*dt;
      hsquare[i][j] = Fplus*hplus_intpolat[i][j]*Fplus*hplus_intpolat[i][j] + Fcross*hcross_intpolat[i][j]*Fcross*hcross_intpolat[i][j];
    }
  }

  for(i=0;i<nfiles;i++){
    htmp = hsquare[i][0];
    ttmp = 0;
    for(j=1;j<l;j++){
      if(hsquare[i][j] > htmp){
	htmp = hsquare[i][j];  
	ttmp = j;
      }  
    }
    tpeak[i] = ttmp;
  }

  free_matrix(hsquare);
}

struct options_t {
 	CHAR *infileindex;
	CHAR *channel;
	INT4 nfiles;

	int verbose;
	int printrawdata;
	int printinterpolateddata;

	INT4 samplerate;
	REAL4 lframe ;
	INT4 useZMWaveforms;
	INT4 useOttWaveforms;
	INT4 useWarrenWaveforms;
	INT4 useKuduWaveforms;
};

static void set_option_defaults(struct options_t *options)
{
	options->infileindex = NULL;
	options->channel = NULL;
	options->nfiles = 0;

	options->verbose = FALSE;
	options->printrawdata = FALSE;
	options->printinterpolateddata = FALSE;

	options->samplerate = 0;
	options->lframe = 0;
	options->useZMWaveforms = 0;
	options->useOttWaveforms = 0;
	options->useWarrenWaveforms = 0;
	options->useKuduWaveforms = 0;
}


static void parse_command_line(int argc, char **argv, struct options_t *options)
{
	struct option long_options[] = {
		/* these options set a flag */
		{"verbose",                no_argument,        &options->verbose, TRUE},
		{"printrawdata",           no_argument,        &options->printrawdata, TRUE},
		{"printinterpolateddata",  no_argument,        &options->printinterpolateddata, TRUE},
		/* parameters which determine the output xml file */
		{"infileindex",     required_argument,  NULL,  'a'},
		{"channel",         required_argument,  NULL,  'b'},
		{"nfiles",          required_argument,  NULL,  'c'},
		{"sample-rate",     required_argument,  NULL,  'd'},
		{"frame-length",    required_argument,  NULL,  'e'},
		{"help",            no_argument,        NULL,  'o'}, 
		{NULL, 0, NULL, 0}
	};
	int c;
	int option_index;

	do {
		switch(c = getopt_long(argc, argv, "a:b:c:d:e:", long_options, &option_index)) {
			case -1:
			case 0:
			break;

			case 'a':
			options->infileindex = optarg;
			if ( ! strcmp("zm", optarg))
			  options->useZMWaveforms = 1;
			if ( ! strcmp("ott", optarg))
			  options->useOttWaveforms = 1;
			if ( ! strcmp("warren", optarg))
			  options->useWarrenWaveforms = 1;
			break;

			case 'b':
			/*
			 * output cache file name
			 */
			options->channel = optarg;
			break;	

			case 'c':
			/*
			 * output cache file name
			 */
			options->nfiles = atoi(optarg);
			break;

			case 'd':
			/*
			 * output cache file name
			 */
			options->samplerate = atoi(optarg);
			break;		

			case 'e':
			/*
			 * output cache file name
			 */
			options->lframe = atoi(optarg);
			break;

			case ':':
			case '?':
			case 'o':
			default:
			/*
			 * print usage
			 */
			LALPrintError(USAGE, *argv);
			exit(1);
		}
	} while(c != -1);

	if(optind < argc) {
		fprintf(stderr, "extraneous command line arguments:\n");
		while(optind < argc)
			fprintf(stderr, "%s\n", argv[optind++]);
		exit(2);
	}

	if(!options->infileindex) {
		LALPrintError( "Input file naming index must be specified\n" );
		exit(3);
	}

	if(!options->lframe) {
		LALPrintError( "Frame length in seconds must be specified\n" );
		exit(4);
	}

	if(!(options->useZMWaveforms || options->useOttWaveforms || options->useWarrenWaveforms)) {
		LALPrintError( "Invalid input for --infileindex; Must be one of zm/ott/warren\n" );
		exit(5);
	}
}


/*void swap4(int n,char *here);
  int readZM(const char* fname,  float **time,  float **amplitude); 
  int readOtt(const char* fname,  float **time,  float **amplitude);*/ 

int main(int argc, char **argv)
{
  static LALStatus stat;                                                      

  struct options_t options;       /* command line */

  INT4 nfiles = 0;                /* no. of input files */
  REAL4 dt = 1.0/16384.0;         /* default sample rate */

  CHAR filename[30], outputfile[30];

  FILE *fp;

  int i,j,l,tshift,N[2048];
  int lframe ;
  float n[2048];
  float x; 
  int *tpeakindex;
  int tpeakmax;

  float *time, *ampplus, *ampcross;
   
  float **hplus, **hcross, **t;
  float **hplus_intpolat, **hcross_intpolat, **t_intpolat;
  float **hplus_pad, **hcross_pad, **t_pad;

  /* frame output data */
  struct FrFile *frOutFile  = NULL;
  struct FrameH *outFrame   = NULL;

  /* raw input data storage */
  REAL4TimeSeries               zmseries;

  /*******************************************************************
   * PARSE ARGUMENTS (arg stores the current position)               *
   *******************************************************************/
  set_option_defaults(&options);
  parse_command_line(argc, argv, &options);

  /* set some of the variables */
  nfiles = options.nfiles;               /* no. of input files */
  if ( options.samplerate )
    dt = 1.0/options.samplerate;         /* required sample  rate */
  lframe = (int)(options.lframe*options.samplerate);  /* no.of points in the frame */

  /* allocate some memory */
  t = matrix(nfiles,MAX_LENGTH);
  memset( t[0], 0, nfiles*MAX_LENGTH*sizeof(float) );
  hplus = matrix(nfiles,MAX_LENGTH);
  memset( hplus[0], 0, nfiles*MAX_LENGTH*sizeof(float) );
  hcross = matrix(nfiles,MAX_LENGTH);
  memset( hcross[0], 0, nfiles*MAX_LENGTH*sizeof(float) );

  /******************************************************************
   * Read in the waveforms
   ******************************************************************/
  N[0]=0;
  for(i=0;i<nfiles;i++){
    /* read in Zm waveforms */
    if ( options.useZMWaveforms ){
      if (options.verbose)
	fprintf(stdout,"Reading in file %s-%d.bin\n", options.infileindex, i);
      snprintf(filename,30,"%s-%d.bin", options.infileindex, i);
      N[i] = readZM(filename,&time,&ampplus,&ampcross);

      if (options.printrawdata){
	/* print out in same format as Ott et al */
	snprintf(filename,30,"%s-%d.trh", options.infileindex, i);
	/* write out the ZM waveforms */
	fp = fopen(filename,"w");
	for(j=0;j<N[i];j++){  
	  fprintf(fp, "%e\t%e\t%e\n", *(time+j), *(ampplus+j), 1e-20*(*(ampplus+j)));
	}
	fclose(fp);
      }    
    } 
    /* read in Ott waveforms */
    else if (options.useOttWaveforms ){
      if (options.verbose)
	fprintf(stdout,"Reading in file %s-%d.bin\n", options.infileindex, i);
      snprintf(filename,30,"%s-%d.bin", options.infileindex, i);
      N[i] = readOtt(filename,&time,&ampplus,&ampcross);
    }
    /*read in Warren waveforms*/
    else if (options.useWarrenWaveforms){
      if (options.verbose)
	fprintf(stdout,"Reading in file %s-%d.dat\n", options.infileindex, i);
      snprintf(filename,30,"%s-%d.bin", options.infileindex, i);
      snprintf(filename,30,"%s-%d.dat", options.infileindex, i);
      N[i] = readWarren(filename,&time,&ampplus,&ampcross);

      if (options.printrawdata){
	snprintf(filename,30,"%s-%d-raw.dat", options.infileindex, i);
	/* write out the Warren waveforms */
	fp = fopen(filename,"w");
	for(j=0;j<N[i];j++){  
	  fprintf(fp, "%e %e %e\n", *(time+j), *(ampplus+j), *(ampcross+j));
	}
	fclose(fp);
      }
    }

    if (options.verbose)
      fprintf(stdout,"Read %d lines from file %s-%d\n", N[i],options.infileindex, i);

    /* store the raw time and amplitudes in t & h */
    for(j=0;j<N[i];j++){   
      t[i][j] = *(time+j);
      hplus[i][j] = *(ampplus+j);
      hcross[i][j] = *(ampcross+j);
    }

    free(time);
    free(ampplus);
    free(ampcross);
  }

  /*calculate h+ from A
  for(i=0;i<nfiles;i++){
    for(j=0;j<N[i];j++){
      look into pg.2 of documentation:here angle=90 & R=1cm.
      h[i][j] *= (0.273);             
    }
  }
  */

  /*********************************************************************************
   * Find l: the required array size to store the waveforms and interpolate to fit l
   ********************************************************************************/

  for(i=0;i<nfiles;i++){
    n[i] = (t[i][N[i]-1]-t[i][0])/dt;
  }

  x=n[0];
  for(i=1;i<nfiles;i++){
    if(n[i]>x)
      x = n[i];
  }

  for(i=0;pow(2,i)<x;++i){
  }
  l=pow(2,i);                          /* l is the required arraysize */

  if ( options.verbose )
    fprintf(stdout,"Number of points = %d\n",l);

  /******************************interpolate**************************/

  /* allocate some memory */
  t_intpolat = matrix(nfiles,l);
  memset( t_intpolat[0], 0, nfiles*l*sizeof(float) );
  hplus_intpolat = matrix(nfiles,l);
  memset( hplus_intpolat[0], 0, nfiles*l*sizeof(float) );
  hcross_intpolat = matrix(nfiles,l);
  memset( hcross_intpolat[0], 0, nfiles*l*sizeof(float) );

  build_interpolate(hplus, t, hplus_intpolat, t_intpolat, dt, l, N, nfiles);
  build_interpolate(hcross, t, hcross_intpolat, t_intpolat, dt, l, N, nfiles);
  
  if ( options.printinterpolateddata ){
    /* print out the interpolated data */
    for(i=0;i<nfiles;i++){
      snprintf(filename,30,"%s-%d-interpolated.dat", options.infileindex, i);
      /* write out the ZM waveforms */
      fp = fopen(filename,"w");
      for(j=0;j<l;j++){  
	fprintf(fp, "%e\t%e\t%e\n", t_intpolat[i][j], hplus_intpolat[i][j], hcross_intpolat[i][j]);
      }
      fclose(fp);
    }
  }

  tpeakindex = malloc( nfiles*sizeof(float));
  find_tpeak(hplus_intpolat, hcross_intpolat, dt, nfiles, l, tpeakindex);

  tpeakmax = tpeakindex[0];
  for(i=1;i<nfiles;i++){
    if (tpeakindex[i] > tpeakmax)
      tpeakmax = tpeakindex[i];
  }

  free_matrix(t);
  free_matrix(hplus);
  free_matrix(hcross);
  /*********************************************************************** 
   * Now adjust the time series so that the peak is at the center of the 
   * frame(whose length is specified by the user)
   **********************************************************************/

  if(lframe<2*tpeakmax){
    fprintf(stderr,"Error: Have to increase the frame length. Make it greater than %e secs\n",2*tpeakmax*dt);
    exit(8);
  }

  t_pad = matrix(nfiles,lframe);
  hplus_pad = matrix(nfiles,lframe);
  hcross_pad = matrix(nfiles,lframe);

  tshift = 0;
  for ( i=0; i<nfiles; i++){
    for (j=0; j<lframe; j++){
      t_pad[i][j] = j*dt;
      hplus_pad[i][j] = 0;
      hcross_pad[i][j] = 0;
    }
    tshift = lframe/2 - tpeakindex[i];    
    for (j = 0; j < l; j++){
      hplus_pad[i][j+tshift] = hplus_intpolat[i][j];
      hcross_pad[i][j+tshift] = hcross_intpolat[i][j];
    }
  }

  free_matrix(t_intpolat);
  free_matrix(hplus_intpolat);
  free_matrix(hcross_intpolat);

  if ( options.printinterpolateddata ){
    /* print out the interpolated data */
    for(i=0;i<nfiles;i++){
      snprintf(filename,30,"%s-%d-padded.dat", options.infileindex, i);
      /* write out the waveforms */
      fp = fopen(filename,"w");
      for(j=0;j<lframe;j++){  
	fprintf(fp, "%e\t%e\t%e\n", t_pad[i][j], hplus_pad[i][j], hcross_pad[i][j]);
      }
      fclose(fp);
    }
  }

  /* create the required vector to store the waveform */
  LALSCreateVector( &stat, &(zmseries.data), lframe );
  zmseries.deltaT = dt;
  if ( options.verbose)
    printf("deltaT = %e\n",zmseries.deltaT);
  zmseries.epoch.gpsSeconds = 0;
  zmseries.epoch.gpsNanoSeconds = 0;
  snprintf( zmseries.name, LALNameLength * sizeof(CHAR), "SIM" );

  /*************************************************************
   * Build the frames from the interpolated and padded waveforms
   ************************************************************/
  for(i=0;i<nfiles;i++){
    snprintf(outputfile,30,"plus_%s_%0d",options.infileindex,i);

    memset( zmseries.data->data, 0, zmseries.data->length*sizeof(REAL4) );
    for(j=0;j<lframe;j++){
      zmseries.data->data[j]=hplus_pad[i][j];
    }

    outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
        &zmseries, "strain", outputfile );
  }

  for(i=0;i<nfiles;i++){
    snprintf(outputfile,30,"cross_%s_%0d",options.infileindex,i);

    memset( zmseries.data->data, 0, zmseries.data->length*sizeof(REAL4) );
    for(j=0;j<lframe;j++){
      zmseries.data->data[j]=hcross_pad[i][j];
    }

    outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
        &zmseries, "strain", outputfile );
  }

  free_matrix(t_pad);
  free_matrix(hplus_pad);
  free_matrix(hcross_pad);

  /****************************************************************
   * Write out the frame file containing all the different channels
   ****************************************************************/
  snprintf(filename,30,"HL-SIM_%s-%d-%d.gwf",options.infileindex,zmseries.epoch.gpsSeconds,(int)(zmseries.data->length*zmseries.deltaT));
  frOutFile = FrFileONew( filename, 0 );
  FrameWrite( outFrame, frOutFile );
  FrFileOEnd( frOutFile );

  LALCheckMemoryLeaks();
  return 0;
}










