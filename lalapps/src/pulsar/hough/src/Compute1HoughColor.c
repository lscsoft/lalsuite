/*-----------------------------------------------------------------------
 *
 * File Name: COMPUTE1HOUGHCOLOR.c
 *
 * Authors: Krishnan, B. Sintes, A.M., 
 *
 * Revision: $Id$
 *  To compute the Hough transform for only one point in parameter space
 *-----------------------------------------------------------------------
 */

/* Input shoud be from 
             SFT files and 
	     the same In.data file used for signal injection
             it should also read the times and velocities used in DriveHoughColor
   This code will output only a number
*/

#include "./DriveHoughColor.h"

#ifdef TIMING
#include "./timer/cycle_counter/Intel/GCC/cycle_counter.h"
#endif

/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID (COMPUTE1HOUGHCOLORC, "$Id$");

/* ************************************************************
 * Usage format string. 
 */

#define USAGE "Usage: %s 
 [-d debuglevel] 
        Default value:  lalDebugLevel=1
 [-D directory] 
        Directory where the input SFT data are located. 
	If not set, the program will look into ./data1
 [-I input data file]	
 	This is a string of the input data file use for signal injection.
	It might contain a path.
	If not set, the program will read ./In.data
 [-V time-velocity data file]	
 	This is a string of the time-velocity data file use for signal injection.
	It might contain a path.
	If not set, the program will read ./velocity.data
 [-b search frequency band (in Hz)]
        Bandwith to be analyzed (as in the wide area search code)
        Default: 2.0 Hz
 [-t peak threshold selection ] 
        Threshold relative to the PSD for the selection of peak in the
        time-frequency plane 
        Default: 1.6 (almost optimal)
 [-w running median window size ] 
        To estimate the psd 
 	Default: 135 and 0.696837167218819
 \n"


/* ***************************************************************
 * Constant Declarations.  Default parameters.
 */

INT4 lalDebugLevel=1;
#define MAXFILES 3000 /* maximum number of files to read in a directory */
#define MAXFILENAMELENGTH 256 /* maximum # of characters  of a SFT filename */
#define SFTDIRECTORY "./data1"
#define FILEOUT "./outHM1/HM"      /* prefix file output */
#define FILEINDATA "./In.data"  /* name of the in.data file with inject signal info */
#define FILEVELOCITY "./velocity.data"  /* name: file with time-velocity info */

#define THRESHOLD 1.6 /* thresold for peak selection, with respect to the
                              the averaged power in the search band */
#define F0 250.0          /*  frequency to build the LUT and start search */
#define FBAND 2.0          /* search frequency band  (in Hz) */
#define NFSIZE  0 /* n-freq. span of the cylinder, to account for spin-down search */
#define BLOCKSRNGMED 135 /* Default running median window size */

/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */

typedef struct tagInjectInputData{
  REAL8        Tsft;   /* Tsft_in_sec */
  UINT4        nTsft;  /* number of elements */
  REAL8	       fmin;   /* first_SFT_frequency_in_Hz */
  REAL4        fBand;  /* SFT_freq_band_in_Hz */
  REAL4        sigma;  /* sigma_(std_of_noise.When=0_only_signal_present) */
  REAL4        aPlus;
  REAL4        aCross;
  REAL4        psi;
  REAL8        phi0;
  REAL8        f0;
  REAL8        latitude;   /* of the source in radians */
  REAL8        longitude;  /* of the source in radians */
  REAL8Vector  spindown;   /* SpinOrder and parameters */
  CHAR         timeStampsName[128];  /* timestamps file name */
} InjectInputData;

/*
 *  Functions Declarations (i.e., prototypes). Not declared in DriveHoughColor.h
 */
 
 void ReadInputInjectDataFile(LALStatus           *status,
                             InjectInputData     *injectInData,
                             CHAR                *fname);

	

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv------------------------------------ */
int main(int argc, char *argv[]){

  static LALStatus           status;  /* LALStatus pointer */
  static InjectInputData     injectInData;
  
  static LIGOTimeGPSVector1    timeV;
  static REAL8Cart3CoorVector velV;
  static REAL8Cart3Coor       sourceLocation;
  static REAL8Vector          timeDiffV;
  static REAL8Vector          foft;
  /* ------------------------------------------------------- */
  UINT4   numberCount;
  CHAR  *fnameInData = NULL;
  CHAR  *fnameVelocity = NULL;
  CHAR   *directory = NULL; /* the directory where the SFT  could be */

  CHAR   filelist[MAXFILES][MAXFILENAMELENGTH];

  INT4   mObsCoh, nfSizeCylinder;
  INT8   f0Bin, fLastBin;           /* freq. bin to perform search */
  REAL8  peakThreshold, normalizeThr;
  REAL8  f0, fSearchBand, timeBase;
   
  UINT2  blocksRngMed;

#ifdef TIMING
  unsigned long long start, stop;
#endif
 
 
#ifdef TIMING
   start = realcc();
#endif

  /******************************************************************/
  /*    Set up the default parameters.      */
  /* ****************************************************************/
 
  fSearchBand = FBAND;
  peakThreshold = THRESHOLD;
  nfSizeCylinder = NFSIZE;
  blocksRngMed = BLOCKSRNGMED;
  SUB( RngMedBias( &status, &normalizeThr, blocksRngMed ), &status );
  directory   = SFTDIRECTORY;
  fnameInData = FILEINDATA;
  fnameVelocity = FILEVELOCITY;
  
  /*****************************************************************/
  /*    Parse argument list.  i stores the current position.       */
  /*****************************************************************/
  {  
    INT4  arg;        /* Argument counter */

    arg = 1;
    while ( arg < argc ) {
      /* Parse debuglevel option. */
      if ( !strcmp( argv[arg], "-d" ) ) {
	if ( argc > arg + 1 ) {
	  arg++;
	  lalDebugLevel = atoi( argv[arg++] );
	} else {
	  ERROR( DRIVEHOUGHCOLOR_EARG, DRIVEHOUGHCOLOR_MSGEARG, 0 );
	  LALPrintError( USAGE, *argv );
	  return DRIVEHOUGHCOLOR_EARG;
	}
      }
      /* Parse filename of In.Data file used for signal injection. */
      else if ( !strcmp( argv[arg], "-I" ) ) {
	if ( argc > arg + 1 ) {
	  arg++;
	  fnameInData = argv[arg++];
	} else {
	  ERROR( DRIVEHOUGHCOLOR_EARG, DRIVEHOUGHCOLOR_MSGEARG, 0 );
	  LALPrintError( USAGE, *argv );
	  return DRIVEHOUGHCOLOR_EARG;
	}
      }
       /* Parse filename of Velocity.Data file corresponding to the SFTs. */
      else if ( !strcmp( argv[arg], "-V" ) ) {
	if ( argc > arg + 1 ) {
	  arg++;
	  fnameVelocity = argv[arg++];
	} else {
	  ERROR( DRIVEHOUGHCOLOR_EARG, DRIVEHOUGHCOLOR_MSGEARG, 0 );
	  LALPrintError( USAGE, *argv );
	  return DRIVEHOUGHCOLOR_EARG;
	}
      }
      /* Parse directory SFT path  option. */
      else if ( !strcmp( argv[arg], "-D" ) ) {
	if ( argc > arg + 1 ) {
	  arg++;
	  directory = argv[arg++];
	} else {
	  ERROR( DRIVEHOUGHCOLOR_EARG, DRIVEHOUGHCOLOR_MSGEARG, 0 );
	  LALPrintError( USAGE, *argv );
	  return DRIVEHOUGHCOLOR_EARG;
	}
      }
      /* Parse search frequency band option. */
      else if ( !strcmp( argv[arg], "-b" ) ) {
	if ( argc > arg + 1 ) {
	  arg++;
	  fSearchBand = atof(argv[arg++]);
	} else {
	  ERROR( DRIVEHOUGHCOLOR_EARG, DRIVEHOUGHCOLOR_MSGEARG, 0 );
	  LALPrintError( USAGE, *argv );
	  return DRIVEHOUGHCOLOR_EARG;
	}
      }
      /* Parse peak threshold option. */
      else if ( !strcmp( argv[arg], "-t" ) ) {
	if ( argc > arg + 1 ) {
	  arg++;
	  peakThreshold = atof(argv[arg++]);
	} else {
	  ERROR( DRIVEHOUGHCOLOR_EARG, DRIVEHOUGHCOLOR_MSGEARG, 0 );
	  LALPrintError( USAGE, *argv );
	  return DRIVEHOUGHCOLOR_EARG;
	}
      }
      /* Parse Running median window size. */
      else if ( !strcmp( argv[arg], "-w" ) ) {
	if ( argc > arg + 1 ) {
	  arg++;
	  blocksRngMed = atoi( argv[arg++] );
	  SUB( RngMedBias( &status, &normalizeThr, blocksRngMed ), &status );  
	} else {
	  ERROR( DRIVEHOUGHCOLOR_EARG, DRIVEHOUGHCOLOR_MSGEARG, 0 );
	  LALPrintError( USAGE, *argv );
	  return DRIVEHOUGHCOLOR_EARG;
	}
      }
     /* Unrecognized option. */
      else {
	ERROR( DRIVEHOUGHCOLOR_EARG, DRIVEHOUGHCOLOR_MSGEARG, 0 );
	LALPrintError( USAGE, *argv );
	return DRIVEHOUGHCOLOR_EARG;
      }
    } /* End of argument parsing loop. */
  }
  /******************************************************************/
  /******************************************************************/
 
 
  /******************************************************************/
  /* reading the Indata file */
  /******************************************************************/
 {  
   REAL8UnitPolarCoor sourceLocPolar;
   
   SUB( ReadInputInjectDataFile(&status, &injectInData, fnameInData), &status );
   /* remember to free the memory before the end of the program of
                  injectIndata.spindown.data      */
   f0 = injectInData.f0;
  
   /* note comments in makefakedata_v2 are wrong */
   sourceLocPolar.delta = injectInData.latitude;
   sourceLocPolar.alpha = injectInData.longitude;
  
   sourceLocation.x= cos(sourceLocPolar.delta)* cos(sourceLocPolar.alpha);
   sourceLocation.y= cos(sourceLocPolar.delta)* sin(sourceLocPolar.alpha);
   sourceLocation.z= sin(sourceLocPolar.delta);
 }
  /******************************************************************/
  /* Looking into the SFT data files to get the names and how many there are*/
  /******************************************************************/
  { 
    CHAR     command[256];
    glob_t   globbuf;
    UINT4    j;
     
    strcpy(command, directory);
    strcat(command, "/*SFT*.*");
    
    globbuf.gl_offs = 1;
    glob(command, GLOB_ERR, NULL, &globbuf);
    
    if(globbuf.gl_pathc==0)
      {
	fprintf(stderr,"No SFTs in directory %s ... Exiting.\n", directory);
	return 1;  /* stop the program */
      }
    
    /* we will read up to a certain number of SFT files, but not all 
       if there are too many ! */ 
    mObsCoh = MIN (MAXFILES, globbuf.gl_pathc);
    
    /* Remember to do the following: 
       globfree(&globbuf); after reading the file names. The file names are 
       globbuf.gl_pathv[fileno]   that one can copy into whatever as:
       strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);  */
    
    for (j=0; j < mObsCoh; j++){
      strcpy(filelist[j],globbuf.gl_pathv[j]);
    }
    globfree(&globbuf);	
  }

  /* ****************************************************************/
  /*  Reading the first headerfile of the first SFT  */
  /* ****************************************************************/
  {
    SFTHeader1    header;
    CHAR   *fname = NULL;
    
    fname = filelist[0];
    SUB( ReadSFTbinHeader1( &status, &header, fname ), &status );
   /* SUB( ReadSFTbinHeader1( &status, &header,&(filelist[0]) ), &status ); */
 
    timeBase= header.timeBase; /* Coherent integration time */
  }
  /* ****************************************************************/
  /*  Reading the time & detector-velocity corresponding to each SFT  */
  /* ****************************************************************/
  velV.length = mObsCoh;
  velV.data = NULL;
  velV.data = (REAL8Cart3Coor *)LALMalloc(mObsCoh*sizeof(REAL8Cart3Coor));

  timeV.length = mObsCoh;
  timeV.time = NULL;  
  timeV.time = (LIGOTimeGPS *)LALMalloc(mObsCoh*sizeof(LIGOTimeGPS));

  {
    UINT4   j; 
    FILE   *fp = NULL;
    INT4    r;

    fp = fopen( fnameVelocity, "r");
    if (fp==NULL){
      fprintf(stderr,"Unable to find file %s\n",fnameVelocity);
      return 1;
    }
    /* read data format:  INT4 INT4  REAL8 REAL8 REAL8 */
    for (j=0; j<mObsCoh;++j){
      r= fscanf(fp, "%d %d %lf %lf %lf\n",
                &(timeV.time[j].gpsSeconds),
		&(timeV.time[j].gpsNanoSeconds), 
		&(velV.data[j].x), &(velV.data[j].y),&(velV.data[j].z) );
		
      if ( r !=5 ) {
       fprintf(stderr,"Unable to assign how many SFTs from %s\n",fnameVelocity);
       return 1;
      }   
    }
    fclose(fp);
  }
  /******************************************************************/
  /* compute the time difference relative to startTime for all SFT */
  /******************************************************************/
  timeDiffV.length = mObsCoh;
  timeDiffV.data = NULL; 
  timeDiffV.data = (REAL8 *)LALMalloc(mObsCoh*sizeof(REAL8));

  {   
    REAL8   t0, ts, tn;
    UINT4   j; 

    ts = timeV.time[0].gpsSeconds;
    tn = timeV.time[0].gpsNanoSeconds * 1.00E-9;
    t0=ts+tn;
    timeDiffV.data[0] = 0.0;

    for(j=1; j< velV.length; ++j){
      ts = timeV.time[j].gpsSeconds;
      tn = timeV.time[j].gpsNanoSeconds * 1.00E-9;  
      timeDiffV.data[j] = ts+tn -t0; 
    }  
  }
 
 
  /* ****************************************************************/
  /* Computing the frequency path f(t) = f0(t)* (1+v/c.n)  */
  /* ****************************************************************/

  foft.length = mObsCoh;
  foft.data = NULL;
  foft.data = (REAL8 *)LALMalloc(mObsCoh*sizeof(REAL8));
  
  {
    REAL8   f0new, vcProdn, timeDiffN;
    UINT4   j,i, nspin, factorialN; 
  
    nspin = injectInData.spindown.length;
    for (j=0; j<mObsCoh; ++j){  /* loop for all different time stamps */
      vcProdn = velV.data[j].x * sourceLocation.x
              + velV.data[j].y * sourceLocation.y
              + velV.data[j].z * sourceLocation.z;
      f0new = f0;
      factorialN = 1;
      timeDiffN = timeDiffV.data[j];
      
      for (i=0; i<nspin;++i){ /* loop for spin-down values */
        factorialN *=(i+1);
        f0new += injectInData.spindown.data[i]* timeDiffN / factorialN;
	timeDiffN *= timeDiffN;
      }
      foft.data[j] = f0new * (1.0 +vcProdn);
    }
    
  }
 
 
  /* ****************************************************************/
  /* reading from SFT, times and generating peakgrams  and producing the
      number-count*/
  /* ****************************************************************/

  numberCount=0;
  

  { 
    COMPLEX8SFTData1  sft;
    REAL8PeriodoPSD   periPSD;
    UCHARPeakGram     pg1;
    
    INT4   length, fWings;
    REAL8  threshold;
    UINT4  j, index; 
    CHAR   *fname = NULL;
  
    /* bandwith to be read should account for Doppler effects and 
       possible spin-down-up */

    f0Bin = f0*timeBase;     /* initial frequency to be analyzed */

    length =  fSearchBand*timeBase; 
    fLastBin = f0Bin+length;   /* final frequency to be analyzed */

 /* To accomodate for running median 
 *    fWings =  floor( fLastBin * VTOT +0.5) + nfSizeCylinder;
 */
    fWings =  floor( fLastBin * VTOT +0.5) + nfSizeCylinder + blocksRngMed;
    length = 1 + length + 2*fWings;

    sft.length = length;
    sft.fminBinIndex = f0Bin-fWings;
    sft.data = NULL;
    sft.data = (COMPLEX8 *)LALMalloc(length* sizeof(COMPLEX8));

    periPSD.periodogram.length = length;
    periPSD.periodogram.data = NULL;
    periPSD.periodogram.data = (REAL8 *)LALMalloc(length* sizeof(REAL8));
    periPSD.psd.length = length;
    periPSD.psd.data = NULL;
    periPSD.psd.data = (REAL8 *)LALMalloc(length* sizeof(REAL8));

    threshold = peakThreshold/normalizeThr; 

    pg1.length = length;
    pg1.data = NULL;
    pg1.data = (UCHAR *)LALMalloc(length* sizeof(UCHAR));

    for (j=0; j < mObsCoh; j++)
      {
      fname = filelist[j];
      SUB( ReadCOMPLEX8SFTbinData1( &status, &sft, fname ), &status );
      SUB( COMPLEX8SFT2Periodogram1(&status, &periPSD.periodogram, &sft), &status );
  
      /* for color noise */    
      SUB( LALPeriodo2PSDrng( &status, 
	  &periPSD.psd, &periPSD.periodogram, &blocksRngMed), &status );
       /* SUB( Periodo2PSDrng( &status, 
	  &periPSD.psd, &periPSD.periodogram, &blocksRngMed), &status ); */
          
      SUB( LALSelectPeakColorNoise(&status,&pg1,&threshold,&periPSD), &status); 

      index = floor( foft.data[j]*timeBase -sft.fminBinIndex+0.5); 
      numberCount+=pg1.data[index]; /* adds 0 or 1 to the counter*/

 
      /******************************************************************/
      /* for debugging purposes only */
      /*FILE *fptemp=NULL;
      INT4 counter;
      fptemp = fopen("tempout","w"); 
      for (counter=0; counter < length; counter++)
	{ 
	  fprintf(fptemp, "%d   %g\n", counter+1, periPSD.periodogram.data[counter]); 
	}
	fclose(fptemp); */
      /******************************************************************/

      /* for debugging purposes only */
      /* if (pg1.data[index] == 0) fprintf(stdout, "%d\n", j); */
      /* fprintf(stdout, "%d\n", index); */


      }

    LALFree(sft.data);
    LALFree(periPSD.periodogram.data);
    LALFree(periPSD.psd.data);
    LALFree(pg1.data);
  }
  
  /******************************************************************/
  /* printing result in stdout */
  /******************************************************************/
   /* ? fprintf(stdout, " %20.15f %12.8e\n", numberCount, comand line arguments ?); */
  {
  
    INT4 j;
    fprintf(stdout, "%d  %f  %f  %f  %f  %f", numberCount, injectInData.f0, injectInData.latitude, injectInData.longitude, injectInData.psi, injectInData.phi0);
    /* write spindown parameters if any */
    if (injectInData.spindown.length > 0)
      {
	for (j=0; j<injectInData.spindown.length; j++)
	  fprintf(stdout, " %e", injectInData.spindown.data[j]);
      }
    fprintf(stdout, "\n");
  }
  /******************************************************************/
  /* Free memory and exit */
  /******************************************************************/
  {
    LALFree(timeV.time);
    LALFree(timeDiffV.data);
    LALFree(velV.data);
    LALFree(foft.data);
 
 
    LALCheckMemoryLeaks();
  } 

#ifdef TIMING
  stop = realcc();
  printf(" All: %llu\n", stop-start);
#endif

  INFO( DRIVEHOUGHCOLOR_MSGENORM );
  return DRIVEHOUGHCOLOR_ENORM;
}

/******************************************************************/
/*  estimate psd from periodogram using runing-median  */
/*  this is just a wrapper to Mohanty's rngmed */
/******************************************************************/
void Periodo2PSDrng (LALStatus  *status,
                     REAL8Periodogram1    *psd,
                     REAL8Periodogram1    *peri,
		     UINT2                *blocksRNG){ 

  INT4   length, j;
  UINT2  blockSize, nblocks2;
  REAL8  *data;
  REAL8  *medians;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "Periodo2PSDrng", COMPUTE1HOUGHCOLORC);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (psd,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (peri, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (blocksRNG, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  psd->epoch.gpsSeconds     = peri->epoch.gpsSeconds;
  psd->epoch.gpsNanoSeconds = peri->epoch.gpsNanoSeconds;
  psd->timeBase     = peri->timeBase;
  psd->fminBinIndex = peri->fminBinIndex;
  
  length = peri->length;
  blockSize = *blocksRNG;
  
  ASSERT (length==psd->length,  status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  ASSERT (blockSize,  status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  ASSERT (blockSize<=length, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  
  if (length > 0){
    ASSERT (peri->data,   status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (psd->data,   status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
    
    data=peri->data;
    nblocks2 = floor(blockSize/2.0);
    medians = psd->data+nblocks2;
    
    /* calling Mohanty's function. It is not LAL compliant */
    rngmed(data, length, blockSize, medians);
    
    for (j=0; j<nblocks2 ; ++j){psd->data[j] = medians[0]; }
    for (j=nblocks2+length-blockSize+1; j<length ; ++j){
      psd->data[j] = medians[length-blockSize];
    }     
  }
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* *************************************************************************/
/* this function should read files containing the following:
	60.0    %Tsft_in_sec
	23920   %nTsft
	1283.0  %first_SFT_frequency_in_Hz
	2.0     %SFT_freq_band_in_Hz
	-1.0    %sigma_(std_of_noise.When=0_only_signal_present)
	APLUS   %Aplus
	ACROSS  %Across
	PSI     %psi
	PHI0    %phi0
	FREQVAL %f0
	0.3766960246    %latitude_in_degrees
	5.1471621319    %longitude_in_degrees
	1       %SpinOrder
	-8.6626e-14 %1stsdp
	TS1     %name_of_time-stamps_file
*/
/* *************************************************************************/
void ReadInputInjectDataFile(LALStatus  *status,
                   InjectInputData     *injectInData,
                   CHAR                *fname) {

  FILE      *fp = NULL; 
  CHAR       dmp[128]; 
  INT4       r, msp, i;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "ReadInputInjectDataFile", COMPUTE1HOUGHCOLORC);
  ATTATCHSTATUSPTR (status);

 /*   Make sure the arguments are not NULL: */
  ASSERT (injectInData,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (fname, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  fp = fopen( fname, "r"); 
  ASSERT (fp, status, DRIVEHOUGHCOLOR_EFILE,  DRIVEHOUGHCOLOR_MSGEFILE); 

  /*Tsft REAL8*/
  r=fscanf(fp, "%lf %s\n",&(injectInData->Tsft),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
  
  r=fscanf(fp, "%d %s\n",&(injectInData->nTsft),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
  
  /*fmin REAL8*/
  r=fscanf(fp, "%lf %s\n",&(injectInData->fmin),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
  
  /*Band REAL4*/
  r=fscanf(fp, "%f %s\n",&(injectInData->fBand),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
 
  /*sigma REAL4*/
  r=fscanf(fp, "%f %s\n",&(injectInData->sigma),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
  
  /*Ap REAL4*/
  r=fscanf(fp, "%f %s\n",&(injectInData->aPlus),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
 
  /*Ac REAL4*/
  r=fscanf(fp, "%f %s\n",&(injectInData->aCross),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
 
  /*psi REAL4*/
  r=fscanf(fp, "%f %s\n",&(injectInData->psi),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 

  /*phi0 REAL8*/
  r=fscanf(fp, "%lf %s\n",&(injectInData->phi0),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 

  /*f0 REAL8*/
  r=fscanf(fp, "%lf %s\n",&(injectInData->f0),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
 
  /*alpha - from input file in  radians.*/
  r=fscanf(fp, "%lf %s\n",&(injectInData->latitude),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
 
  /*delta - from input file in radians.*/
  r=fscanf(fp, "%lf %s\n",&(injectInData->longitude),dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 


  /* max spin-down parameter order */
  r=fscanf(fp, "%d %s\n",&msp, dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 

  injectInData->spindown.length = msp;
  injectInData->spindown.data   = NULL;
  
  /* if there are spin down parameters read them in */
  if (msp > 0){
    injectInData->spindown.data =(REAL8 *)LALMalloc(msp*sizeof(REAL8));
    for (i=0;i<msp;i++){
      r=fscanf(fp, "%le %s\n",&(injectInData->spindown.data[i]),dmp);
      ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
    }/*endfor over different spindown orders*/
  }/*endif there were spindown values at all*/
  
  /* timestamps file name */
  r=fscanf(fp, "%s %s\n",injectInData->timeStampsName,dmp);
  ASSERT (r==2, status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG); 
  
  fclose(fp);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}
  
 
  
