
/*-----------------------------------------------------------------------
 *
 * File Name: makefakedata_v1.c
 *
 * Authors: Papa, M.A., 
 *
 * Revision: $Id: extractSFTband.c,v 1.2 2004/06/22 21:03:56 siemens Exp $
 *
 *-----------------------------------------------------------------------
 */

/*
 * 1.  An author and Id block
 */

/************************************ <lalVerbatim file="makefakedataCV">
Author: Papa, M.A. 
$Id: extractSFTband.c,v 1.2 2004/06/22 21:03:56 siemens Exp $
************************************* </lalVerbatim> */

#include <lal/LALStdlib.h>
NRCSID (MAKEFAKEDATAC, "$Id: extractSFTband.c,v 1.2 2004/06/22 21:03:56 siemens Exp $");

/***************************************************/

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <glob.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include "lal/LALVersion.h"

int mycalls=0;
int myclears=0;
void *ptr;

/* For debugging memory allocation */
#if(0)
#define LALMalloc(how)  ptr=LALMalloc(how);mycalls++,fprintf(stderr,"Allocating %d bytes with %d calls tp %p\n",(INT4)(how),mycalls,ptr)
#define LALFree(quantity) LALFree(ptr=quantity); myclears++; fprintf(stderr,"Clearing with %d calls to LALFree at %p\n",myclears,ptr)
#endif

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, MAKEFAKEDATAC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              MAKEFAKEDATAC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( MAKEFAKEDATAC_ESUB, MAKEFAKEDATAC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return MAKEFAKEDATAC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

/***************************/

#define MAXFILES 40000         /* Maximum # of files in a directory */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

/*Default input data file name*/
char *inDataFilename="InExtract.data";
char *basefilename="TEST_SFT";
char *noisedir="./";
char filelist[MAXFILES][MAXFILENAMELENGTH];

/* timebaseline of SFT in sec, band SFT in Hz */
REAL4 Tsft,B;
/* smallest frequency in the band B */
REAL8 fmin;
/* How many SFTS we'll produce*/
INT4 nTsft;

/*How many samples in the original SFT*/
INT4 nsamplesOriginal;

/*This will hold the SFT*/
COMPLEX8Vector *fvec=NULL;

/* Store time stamps from SFT data */
LIGOTimeGPS *timestamps=NULL;


INT4 lalDebugLevel=3;

/* Prototypes for the functions defined in this file */
int make_filelist(LALStatus *);
int read_file(LALStatus *, int argc, char *argv[]);
int prepare_fvec(LALStatus *);
int read_noise(LALStatus *, int iSFT);
int write_SFTS(LALStatus *, int iSFT);
int freemem(LALStatus *);

int main(int argc,char *argv[]) {

  static LALStatus status;
  INT4 iSFT;

  fvec=NULL;
   
  /* read input variables from input data file */
  if (read_file(&status,argc,argv))
    return 1;

  if (make_filelist(&status))
    return 1;

  timestamps=(LIGOTimeGPS *)LALMalloc(nTsft*sizeof(LIGOTimeGPS)); 
 
  for (iSFT=0;iSFT<nTsft;iSFT++){

   /* read data, create fvec vector and populate it with the data */ 
   if (read_noise(&status, iSFT))
      return 1;

   if (write_SFTS(&status, iSFT))
     return 1;
  }
  
  if (freemem(&status))
    return 1;
  
  LALCheckMemoryLeaks(); 

  return 0;
}



/* This routine frees up all the memory */
int freemem(LALStatus* status){


  /* Clean up FFTs of signal and of noise - each a complex8vector */
  LALCDestroyVector(status, &fvec);
 /* Clean up timestamps */
  LALFree(timestamps);

  return 0;
}





int make_filelist(LALStatus* status) {

  INT4 fileno=0;
  char command[256];
  glob_t globbuf;

  strcpy(command,noisedir);
  strcat(command,"*SFT*");
  
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);

  /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */
  while (fileno< globbuf.gl_pathc) 
    {
      strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
      if (fileno > MAXFILES)
	{
	  fprintf(stderr,"Too many files in directory! Exiting... \n");
	  return 1;
	}
    }
  globfree(&globbuf);


  return 0;

}


int read_noise(LALStatus* status, int iSFT) {

  FILE *fp;
  size_t errorcode;
  int imin, len2,i;
  REAL4 norm;

  struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
  } header;
  
  /**********************/


  /* open file and get info from it*/
  fp=fopen(filelist[iSFT],"r");
  /* read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",filelist[iSFT]);
      return 1;
    }

  /* check that data is correct endian order */
  if (header.endian!=1.0)
    {
      fprintf(stderr,"First object in file %s is not (double)1.0!\n",filelist[iSFT]);
      fprintf(stderr,"It could be a file format error (big/little\n");
      fprintf(stderr,"endian) or the file might be corrupted\n\n");
      return 2;
    }

  /*************************************/
  Tsft=header.tbase;
  imin=(INT4)(fmin*Tsft+0.5);
  fmin=imin/Tsft;

  /* check frequency range */
  if ((fmin*Tsft < header.firstfreqindex) ||
      ((fmin*Tsft+(ceil(B*Tsft))) > header.firstfreqindex + header.nsamples)){
    fprintf(stderr,"Frequency band of noise data out of range !\n");
    return 1;
  }

  nsamplesOriginal=header.nsamples;
  timestamps[iSFT].gpsSeconds=header.gps_sec;
  timestamps[iSFT].gpsNanoSeconds=header.gps_nsec;

  /* seek to position */
  if (0 != fseek(fp, (imin - header.firstfreqindex) * 4 * 2, SEEK_CUR)){
    fprintf(stderr,"file too short (could'n fssek to position\n");
    return 1;
  }

  len2=(INT4)(B*Tsft+0.5)+1;
  if (iSFT == 0) {
    LALCCreateVector(status, &fvec, (UINT4)len2);
  }

  norm=(REAL4)(len2)/(REAL4)header.nsamples;

  /* read the data */
  if (1 != fread(fvec->data,len2*2*4,1,fp)) {
      fprintf(stderr,"Could not read the data \n");
    return 1;
  }

  fclose(fp);

  for (i=0;i<len2;i++){

    fvec->data[i].re=norm*fvec->data[i].re;
    fvec->data[i].im=norm*fvec->data[i].im;    
  }


  return 0;
}





int write_SFTS(LALStatus* status, int iSFT){


  FILE *fp;
  REAL4 rpw,ipw;
  char filename[256], filenumber[16];
  int i,errorcode;

  struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
  } header;
  
  /* Open SFT data file */
  strcpy(filename,basefilename);
  sprintf(filenumber,".%05d",iSFT); /*ACHTUNG: used to be iSFT+1 - now starts from .00000*/
  strcat(filename,filenumber);
  fp=fopen(filename,"w");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",filename);
    return 1;
  }

  header.endian=1.0;
  header.gps_sec=timestamps[iSFT].gpsSeconds;
  header.gps_nsec=timestamps[iSFT].gpsNanoSeconds;
  header.tbase=Tsft;
  header.firstfreqindex=(INT4)(fmin*Tsft+0.5);
  header.nsamples=fvec->length; 
  
  /* write header */
  errorcode=fwrite((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1){
    printf("Error in writing header into file!\n");
    return 1;
  }

  for (i=0;i<fvec->length;i++){

    rpw=fvec->data[i].re;
    ipw=fvec->data[i].im;

    errorcode=fwrite((void*)&rpw, sizeof(REAL4),1,fp);  
    if (errorcode!=1){
      printf("Error in writing data into SFT file!\n");
      return 1;
    }
    errorcode=fwrite((void*)&ipw, sizeof(REAL4),1,fp);  
    if (errorcode!=1){
      printf("Error in writing data into SFT file!\n");
      return 1;
    }
        
  }
  
  fclose(fp);
  return 0;  
  
}


int read_file(LALStatus* status, int argc,char *argv[]) {
  
  char filename[256], dmp[128];
  int c, errflg = 0;
  int r;
  FILE *fp;
  extern char *optarg;
  
  /* scan through the list of arguments on the command line 
     and get the input data filename*/
  
  while (!errflg && ((c = getopt(argc, argv,"hn:d:N:b:f:"))!=-1))
    switch (c) {
    case 'n':
      /* Name of output SFT data file */
      basefilename=optarg;
      break;
    case 'd':
      /* Name of input (noise) SFT file */
      noisedir=optarg;
      break;
    case 'N':
      /* Name of input (noise) SFT file */
      nTsft=atof(optarg);
      break;
    case 'b':
      /* Name of input (noise) SFT file */
      B=atof(optarg);
      break;
    case 'f':
      /* Name of input (noise) SFT file */
      fmin=atof(optarg);
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stderr,"Arguments are:\n");
      fprintf(stderr,"\t-n\tCharacter String\t( Name of output data file)\n");
      fprintf(stderr,"\t-d\tCharacter String\t( Name of input data directory)\n");
      fprintf(stderr,"\t-N\tINTEGER\t(Number of SFTs)\n");
      fprintf(stderr,"\t-b\tFLOAT\t(Band to extract)\n");
      fprintf(stderr,"\t-f\tFLOAT\t(Starting frequency)\n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
  
  return errflg;
}





