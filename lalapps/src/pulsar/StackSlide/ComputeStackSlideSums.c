/*********************************************************************************/
/*                                                                               */
/* File: ComputeStackSlideSums.c                                                 */
/* Purpose: Drive StackSlide jobs under Condor and the Grid                      */
/* Author: Mendell, G. and Landry, M.                                            */
/* Started: 03 December 2003                                                     */
/* $Id$         */
/*                                                                               */
/*********************************************************************************/
/* Notes: */
/* BLK = one Block of frequency domain data */
/* STK = one Stack of frequency domain data; Block are coherenty combined to make Stacks. */
/* SUM = Summed Stacks after Sliding */

/* TO DO LIST: */
/* Add handling of status pointer in main function */
/* Get units from input data.  For now these are hard coded based on input data type. */
/* Currently inputDataTypeFlag == 0 = SFTs is the only supported input data type. Need to add types 1 = PSDs, 2 = F-stat */
/* Currently normalization code only works when input data is SFTs */
/* Input SFTs must be calibrated; Calibration currently cannot be done from this code. */
/* Currently normalizationFlag == 1 only type of normalization supported */
/* Check that all allocated memory is deallocated */
/* Code currently assumes params->numSUMsPerParamSpacePt = params->duration/params->tSUM = 1; need to upgrade to more general case */
/* Compute Correct Power Stats, snr, confidence, etc... */
/* Need analysis of SUMs to search for events above threshold and width of events; time analyze SUMs to veto lines vs signals. */
/* Input SFTs are not sorted; thus actualStartTime and actualEndTime are not always correct; this only affects search summary table but should be fixed */
/* Check if endtime is OK, or should endtime - tBKL be used to not get SFTs with end time past gpsStartTime + duration? */
/* For now fix DeltaRA, DeltaDec,and DeltaFDeriv1 to default values when finding mismatch during Monte Carlo; need to make these adjustable */

/* REVISIONS: */
/* 12/03/03 gam; Based code on StackSlide.c from LALWrapper DSO */
/* 01/06/04 gam; Fix bug when deallocating BLKData and params->numBLKs != GV.SFTno. */
/* 01/23/04 gam; Glob for *SFT.* rather than *SFT* to avoid non-SFT files that have SFT in their name. */
/* 01/28/04 gam; Include what to glob for in params->sftDirectory rather than hard coding */
/* 04/15/04 gam; Fix help when no command line arguments */
/* 05/07/04 gam; add alternative to using glob */
/* 05/07/04 gam; Do not check d_type when using readdir since it is often UNKOWN. */
/* 05/11/04 gam; Add RunStackSlideMonteCarloSimulation to software inject signals into the SFTs for Monte Carlo simulations */
/* 05/21/04 gam; Normalize the SFTs from LALSignalToSFTs as per the normalization in makefakedata_v2.c and lalapps/src/pulsar/make_sfts.c. */
/* 05/21/04 gam; Save the input SFTs and Parameter space; run Monte Carlo simulation on the parameter space with random offsets */
/* 05/26/04 gam; Continue work on Monto Carlo code. */
/* 05/26/04 gam; Change finishedSUMs to finishedSUMs; add startSUMs; defaults are TRUE; use to control I/O during Monte Carlo */
/* 05/26/04 gam; Add whichMCSUM = which Monte Carlo SUM; default is -1. */
/* 05/28/04 gam; Use LALUniformDeviate from LAL utilities package (include <lal/Random.h>) to generate random mismatch during Monte Carlo. */
/* 06/01/04 gam; Make sure CHECKSTATUSPTR called after every LAL function call */
/* 06/01/04 gam; pPulsarSignalParams->startTimeGPS and duration should be based on gpsStartTime and duration, not actualStartTime */
/* 06/01/04 gam; For now fix DeltaRA, DeltaDec,and DeltaFDeriv1 to default values when finding mismatch during Monte Carlo */
/* 07/14/04 gam; add option to use function LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs; see LAL inject packge GeneratePulsarSignal.c and .h */
/* 07/19/04 gam; if (params->testFlag & 4) > 0 use LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs for Monte Carlo simulations */
/* 08/02/04 gam; if (params->testFlag & 4) > 0 ComputeSky uses reference time GPSin: params->timeStamps[0].gpsSeconds, params->timeStamps[0].gpsNanoSeconds */
/* 08/02/04 gam; Set pSFTandSignalParams->resTrig = 0 to avoid serious bug in LALFastGeneratePulsarSFTs when using lookup tables (LUTs) for trig calls. */
/* 12/06/04 gam; get params->sampleRate, = effective sample rate, from the SFTs; calculate params->deltaT after reading SFTs. */
/* 12/06/04 gam; add params->gpsEpochStartTimeNan; get gpsEpochStartTime, gpsEpochStartTimeNan, and gpsStartTime from command line; */ 
/* 12/06/04 gam; if (params->testFlag & 8) > 0 use fixed values for psi and cosIota during Monte Carlo simulations */

/*********************************************/
/*                                           */
/* START SECTION: define preprocessor flags  */
/*                                           */
/*********************************************/
/* #define DEBUG_READSFT_CODE */
#define INCLUDE_RUNSTACKSLIDEMONTECARLO_CODE
/* #define DEBUG_MONTECARLOTIMEDOMAIN_DATA */
/* #define PRINT_MONTECARLOTIMEDOMAIN_DATA */
/* #define DEBUG_MONTECARLOSFT_DATA */
/* #define PRINT_MONTECARLOSFT_DATA */
/* #define PRINTCOMPARISON_INPUTVSMONTECARLOSFT_DATA */
/* #define DEBUG_RANDOMTRIALPARAMETERS */
/* #define DEBUG_LALFASTGENERATEPULSARSFTS */
/* #define DEBUG_SETFIXED_RANDVAL */
/* #define PRINT_ONEMONTECARLO_OUTPUTSFT */
/* #define PRINT_MAXPOWERANDBINEACHSFT */
/*********************************************/
/*                                           */
/* END SECTION: define preprocessor flags    */
/*                                           */
/*********************************************/

/*********************************************/
/*                                           */
/* START SECTION: include header files       */
/* (Note that most include statements        */
/*  are in the header files below.)          */
/*                                           */
/*********************************************/
#include "ComputeStackSlideSums.h"
/*********************************************/
/*                                           */
/* END SECTION: include header files         */
/*                                           */
/*********************************************/

NRCSID( COMPUTESTACKSLIDESUMSC, "$Id$");

/*********************************************/
/*                                           */
/* START SECTION: declare global variables   */
/*                                           */
/*********************************************/
GlobalVariables GV;
/*********************************************/
/*                                           */
/* END SECTION: declare global variables     */
/*                                           */
/*********************************************/

/******************************************/
/*                                        */
/* START FUNCTION: main                   */
/*                                        */
/******************************************/
int main(int argc,char *argv[]) 
{ 
  static LALStatus status;
  StackSlideSearchParams *params;  

  if (argc < 2) {
        /* fprintf(stdout, "No command line arguments found.  See README_STACKSLIDE or LAL StackSlide documentation for help\n"); */
        fprintf(stdout, "No command line arguments found. For help see README_ComputeStackSlideSums or check for StackSlide documentation in LAL.\n");
        fflush(stdout);
        return 1;
  }
  
  /* Initialize status */
  status.statusCode = 0;
  status.statusPtr = NULL;
  status.statusDescription = NULL;
    
  /* Allocate memory */
  params = (StackSlideSearchParams *)LALMalloc(sizeof( StackSlideSearchParams ));  
       
  /* Initialize params using the command line arguments */
  StackSlideInitSearch(&status,params,argc,argv);
  INTERNAL_CHECKSTATUS_FROMMAIN(status)  
      
  if (SetGlobalVariables(params)) return 2;
    
  if (params->inputDataTypeFlag == 0) {  
     if (ReadSFTData(params)) return 3; /* Note that for now it is assumed that BLK data is SFT data */
  } else {
     return 3;
  }

  StackSlideConditionData(&status,params);  
  INTERNAL_CHECKSTATUS_FROMMAIN(status)  

  if ( (params->testFlag & 2) > 0 ) {
    /* 05/11/04 gam; Add code to software inject signals into the SFTs for Monte Carlo simulations */    
    RunStackSlideMonteCarloSimulation(&status,params);
    INTERNAL_CHECKSTATUS_FROMMAIN(status);    
  } else {      
    StackSlideApplySearch(&status,params);
    INTERNAL_CHECKSTATUS_FROMMAIN(status)
  }
  
  StackSlideFinalizeSearch(&status,params);
  INTERNAL_CHECKSTATUS_FROMMAIN(status)    
      
  if (Freemem(params)) return 4;

  LALFree(params);  
  
  return 0;

}
/******************************************/
/*                                        */
/* END FUNCTION: main                     */
/*                                        */
/******************************************/

/******************************************/
/*                                        */
/* START FUNCTION: ReadSFTData            */
/*                                        */
/******************************************/
int ReadSFTData(StackSlideSearchParams *params)
{

  /* Based on ReadSFTData in lalapps/src/pulsar/ComputeFStatistic.c written by M.A. Papa and X. Siemens */
  /* Modified this function to read in only number of BLKs requested from SFTs for times requested! */
  
  INT4 blkno=0, fileno=0,ndeltaf,offset; /* Added blkno = index of which BLK we are reading in */
  UINT4 stoptime;
  FILE *fp;
  size_t errorcode;
 
  /* INITSTATUS( status, "ReadSFTData", COMPUTESTACKSLIDESUMSC);
  ATTATCHSTATUSPTR (status); */
  
  /* params->BLKData=(FFT **)LALMalloc(GV.SFTno*sizeof(FFT *));
  params->timeStamps=(LIGOTimeGPS *)LALMalloc(GV.SFTno*sizeof(LIGOTimeGPS)); */
  /* Code copied from StackSlide.c DSO LALConditionData function: */  
  params->BLKData=(FFT **)LALMalloc(params->numBLKs*sizeof(FFT *));
  params->timeStamps=(LIGOTimeGPS *)LALMalloc(params->numBLKs*sizeof(LIGOTimeGPS));
  /* params->timeIntervals =(gpsTimeInterval *)LALCalloc(params->numBLKs, sizeof(gpsTimeInterval)); */
	
  ndeltaf=GV.ifmax-GV.ifmin+1; 
  stoptime = params->gpsStartTimeSec + (UINT4)floor(params->duration);
  
  for (fileno=0;fileno<GV.SFTno;fileno++)
  {
      /* open FIRST file and get info from it*/

      fp=fopen(GV.filelist[fileno],"r");
      if (fp==NULL) 
	{
	  fprintf(stderr,"Weird... %s doesn't exist!\n",GV.filelist[fileno]);
	  return 1;
	}
      /* Read in the header from the file */
      errorcode=fread((void*)&header,sizeof(header),1,fp);
      if (errorcode!=1) 
	{
	  fprintf(stderr,"No header in data file %s\n",GV.filelist[fileno]);
	  return 1;
	}

      /* Check that data is correct endian order */
      if (header.endian!=1.0)
	{
	  fprintf(stderr,"First object in file %s is not (double)1.0!\n",GV.filelist[fileno]);
	  fprintf(stderr,"It could be a file format error (big/little\n");
	  fprintf(stderr,"endian) or the file might be corrupted\n\n");
	  return 2;
	}
    
      /* Check that the time base is positive */
      if (header.tbase<=0.0)
	{
	  fprintf(stderr,"Timebase %f from data file %s non-positive!\n",
		  header.tbase,GV.filelist[fileno]);
	  return 3;
	}
        
      /* Check that are frequency bins needed are in data set */
      
       if (GV.ifmin<header.firstfreqindex || 
	  GV.ifmax>header.firstfreqindex+header.nsamples) 
	{
	  fprintf(stderr,"Freq index range %d->%d not in %d to %d (file %s)\n",
		GV.ifmin,GV.ifmax,header.firstfreqindex,
		  header.firstfreqindex+header.nsamples,GV.filelist[fileno]);
	  return 4;
	}
		
      /* Check that this SFT is within the times requested */			
      if (((UINT4)header.gps_sec >= params->gpsStartTimeSec) && ((UINT4)header.gps_sec < stoptime)) {	
        /* Put time stamps from file into array */
        params->timeStamps[blkno].gpsSeconds=header.gps_sec;
        params->timeStamps[blkno].gpsNanoSeconds=header.gps_nsec;

        /* Move forward in file */
        offset=(GV.ifmin-header.firstfreqindex)*2*sizeof(REAL4);
        errorcode=fseek(fp,offset,SEEK_CUR);
        if (errorcode) 
	  {
	    perror(GV.filelist[fileno]);
	    fprintf(stderr,"Can't get to offset %d in file %s\n",offset,GV.filelist[fileno]);
	    return 5;
	  }

        /* Make data structures */
        /* ndeltaf=GV.ifmax-GV.ifmin+1; */ /* Do just once above the loop */
        params->BLKData[blkno]=(FFT *)LALMalloc(sizeof(FFT));
        params->BLKData[blkno]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
        params->BLKData[blkno]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
        params->BLKData[blkno]->fft->data->data=(COMPLEX8 *)LALMalloc(ndeltaf*sizeof(COMPLEX8));

        /* Fill in actual SFT data, and housekeeping */
        errorcode=fread((void*)(params->BLKData[blkno]->fft->data->data),ndeltaf*sizeof(COMPLEX8),1,fp);
        params->BLKData[blkno]->fft->epoch=params->timeStamps[blkno];
        params->BLKData[blkno]->fft->f0=GV.ifmin*GV.df;
        params->BLKData[blkno]->fft->deltaF=GV.df;
        params->BLKData[blkno]->fft->data->length=ndeltaf;
	      
        blkno++;	
     
    } /* end if ((header.gps_sec >= params->gpsStartTimeSec) && (header.gps_sec < stoptime)) */
    
    fclose(fp);     /* close file */
    
    if (blkno >= params->numBLKs) {
      	params->finishedBLKs = 1;
      	break;
    }

  } /* end for (fileno=0;fileno<GV.SFTno;fileno++) */

  /* Some requested input BLK data is missing */
  if (!params->finishedBLKs) {
     /* ABORT( status, DRIVESTACKSLIDEH_EMISSINGBLKDATA, DRIVESTACKSLIDEH_MSGEMISSINGBLKDATA); */
     fprintf(stderr,"Some requested input BLK data is missing\n");
     return 6;     
  }
  
  /* CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status); */
  
  return 0;  
}
/******************************************/
/*                                        */
/* END FUNCTION: ReadSFTData              */
/*                                        */
/******************************************/

/******************************************/
/*                                        */
/* START FUNCTION: SetGlobalVariables     */
/*                                        */
/******************************************/
int SetGlobalVariables(StackSlideSearchParams *params)
{

  /* char command[256]; */ /* 01/28/04 gam */
  FILE *fp;
  size_t errorcode;
  REAL8 df;                         /* freq resolution */
  INT4 fileno=0;   
  /* glob_t globbuf; */     /* 05/07/04 gam; add alternative to using glob */
  DIR *pSFTDir;             /* 05/07/04 gam; pointer to SFT directory */
  struct dirent *dirEntry;  /* 05/07/04 gam; structure that holds a directory entry */
  INT2 i, j, sftDirNameLen, indSlash, asteriskCount; /* 05/07/04 gam; Used when parsing params->sftDirectory */
  CHAR sftDirectory[256],sftPattern[256];            /* 05/07/04 gam; parse params->sftDirectory to get sftDirectory and sftPattern. */
    
  /* INITSTATUS( status, "ReadSFTData", COMPUTESTACKSLIDESUMSC);
  ATTATCHSTATUSPTR (status); */
  
  /* 05/07/04 gam; START add alternative to using glob */
  sftDirNameLen = strlen(params->sftDirectory);
  if (sftDirNameLen > 256) {
     fprintf(stderr,"Error: params->sftDirectory cannot contain more than 256 characters. \n");
     return 1;
  }  
  indSlash = -1;
  /* Find the last slash in params->sftDirectory */
  for (i=sftDirNameLen-1;i>=0;i--) {
        #ifdef DEBUG_READSFT_CODE
          fprintf(stdout,"char[%i]  = %c\n",i,params->sftDirectory[i]);
          fflush(stdout);
        #endif
        if(params->sftDirectory[i] == '/') {
           indSlash = i;
           break;
        }
  }
  #ifdef DEBUG_READSFT_CODE
        fprintf(stdout,"indSlash = %i\n",indSlash);
        fflush(stdout);
  #endif  
  /* If the last character is not a / then check for a pattern to match */  
  if ( indSlash < (sftDirNameLen-1) ) {
     asteriskCount = 0;
     for (i=indSlash+1;i<sftDirNameLen;i++) {
        if(params->sftDirectory[i] == '*') {
           asteriskCount++;
        }
     }
     #ifdef DEBUG_READSFT_CODE
        fprintf(stdout,"asteriskCount = %i\n",asteriskCount);
        fflush(stdout);
     #endif
     if (asteriskCount > 0) {
       if (indSlash == -1) {
         strcpy(sftDirectory,"./"); /* params->sftDirectory has no slashes but has asterisks; set sftDirectory to current directory */
       } else {
         strncpy(sftDirectory,params->sftDirectory,indSlash+1);
         sftDirectory[indSlash+1] = '\0';
       }
       j = 0;
       for (i=indSlash+1;i<sftDirNameLen;i++) {
          #ifdef DEBUG_READSFT_CODE
            fprintf(stdout,"char[%i]  = %c\n",i,params->sftDirectory[i]);
            fflush(stdout);
          #endif
          if(params->sftDirectory[i] != '*') {
             sftPattern[j] = params->sftDirectory[i];
             j++;   
          }
       }
       sftPattern[j] = '\0';
     } else {
       /* params->sftDirectory did not end in a slash and contained no asterisks after last slash */
       strcpy(sftDirectory,params->sftDirectory);
       strcat(sftDirectory,"/");  /* sftDirectory must end with a slash */
       sftPattern[0] = '\0';  /* No pattern to match given, so match everything. */
     }
  } else {
    strcpy(sftDirectory,params->sftDirectory);
    sftPattern[0] = '\0';  /* No pattern to match given, so match everything. */
  }
  #ifdef DEBUG_READSFT_CODE
     fprintf(stdout,"sftDirectory, sftPattern = %s, %s\n",sftDirectory,sftPattern);
     fflush(stdout);
  #endif  
  pSFTDir = opendir(sftDirectory);
  dirEntry = readdir(pSFTDir);
  while (dirEntry != NULL) {
     /* if ( strstr(dirEntry->d_name,sftPattern) && (dirEntry->d_type == DT_REG) ) */ /* 05/07/04 gam */
     if ( strstr(dirEntry->d_name,sftPattern) ) {
        #ifdef DEBUG_READSFT_CODE
          fprintf(stdout,"d_name, d_type = %s, %i\n",dirEntry->d_name,dirEntry->d_type);
          fflush(stdout);
        #endif
        strcpy(GV.filelist[fileno],sftDirectory);
        strcat(GV.filelist[fileno],dirEntry->d_name);
        fileno++;
        if (fileno > MAXFILES) {
          fprintf(stderr,"Too many files match params->sftDirectory! Exiting... \n");
          return 1;
        }
     }
     dirEntry = readdir(pSFTDir);
  }  
  closedir(pSFTDir);
  if(fileno==0) {
    fprintf(stderr,"No SFTs found in params->sftDirectory %s ... Exiting.\n",params->sftDirectory);
    return 1;
  }
  /* 05/07/04 gam; END add alternative to using glob */

/* 05/07/04 gam; use preprocessor flag to ignor glob code: */
#ifdef NOTHING      
  /* strcpy(command,params->sftDirectory); */  /* 01/28/04 gam */
  /* strcat(command,"/ *SFT*"); */  /* 01/23/04 gam; Glob for *SFT.* rather than *SFT* */
  /* strcat(command,"/ *SFT.*"); */  /* 01/28/04 gam (Also note added space between / and * to get rid of warning about comment inside comment. */

  globbuf.gl_offs = 1;
  /* glob(command, GLOB_ERR, NULL, &globbuf); */ /* 01/28/04 gam; Include what to glob for in params->sftDirectory rather than hard coding */
  glob(params->sftDirectory, GLOB_ERR, NULL, &globbuf);  

  /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */

  if(globbuf.gl_pathc==0)
    {
      /* fprintf(stderr,"No SFTs in params->sftDirectory %s ... Exiting.\n",params->sftDirectory); */ /* 01/28/04 gam */
      fprintf(stderr,"No SFTs match params->sftDirectory %s ... Exiting.\n",params->sftDirectory);
      return 1;
    }

  while ((UINT4)fileno < globbuf.gl_pathc) 
    {
      strcpy(GV.filelist[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
      if (fileno > MAXFILES)
	{
	  /* fprintf(stderr,"Too many files in params->sftDirectory! Exiting... \n"); */ /* 01/28/04 gam */
	  fprintf(stderr,"Too many files match params->sftDirectory! Exiting... \n");
	  return 1;
	}
    }
  globfree(&globbuf);
#endif

  GV.SFTno=fileno; /* remember this is 1 more than the index value */
        
  #ifdef DEBUG_READSFT_CODE
          fprintf(stdout,"The first SFT in GV.filelist is %s\n",GV.filelist[0]);
          fflush(stdout);
  #endif

  /* open FIRST file and get info from it*/
  fp=fopen(GV.filelist[0],"r");
  /* read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",GV.filelist[0]);
      return 1;
    }

  /* check that data is correct endian order */
  if (header.endian!=1.0)
    {
      fprintf(stderr,"First object in file %s is not (double)1.0!\n",GV.filelist[0]);
      fprintf(stderr,"It could be a file format error (big/little\n");
      fprintf(stderr,"endian) or the file might be corrupted\n\n");
      return 2;
    }
  fclose(fp);

  GV.Ti=header.gps_sec;  /* INITIAL TIME */

  /* open LAST file and get info from it*/
  fp=fopen(GV.filelist[fileno-1],"r");
  /* read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",GV.filelist[fileno-1]);
      return 1;
    }
  /* check that data is correct endian order */
  if (header.endian!=1.0)
    {
      fprintf(stderr,"First object in file %s is not (double)1.0!\n",GV.filelist[fileno-1]);
      fprintf(stderr,"It could be a file format error (big/little\n");
      fprintf(stderr,"endian) or the file might be corrupted\n\n");
      return 2;
    }
  fclose(fp);

  GV.Tf=header.gps_sec+header.tbase;  /* FINAL TIME */

  GV.tsft=header.tbase;  /* Time baseline of SFTs */
    
  GV.nsamples=header.nsamples;    /* # of freq. bins */
  
  /* 12/06/04 gam; get params->sampleRate, = effective sample rate, from the SFTs */
  params->sampleRate = 2.0*((REAL8)(header.nsamples))/header.tbase;
  params->deltaT = 1.0/params->sampleRate;

  /* frequency resolution: used only for printing! */
  df=(1.0)/(1.0*header.tbase);
  GV.df=df;
  
  /* Verify meta data is correct */
  if (header.tbase != params->tEffBLK) {
  	/* ABORT( status, DRIVESTACKSLIDEH_EFREQSTEPSIZE, DRIVESTACKSLIDEH_MSGEFREQSTEPSIZE); */
        fprintf(stderr,"Input BLK frequency step size does not agree with that expected\n");
        return 3;	
  }

  /* Set up index of min and max frequency to get from SFTs */  
  GV.ifmin = floor(params->f0BLK*params->tEffBLK + 0.5);
  GV.ifmax = GV.ifmin + params->nBinsPerBLK -1;
  
  /* CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status); */

  return 0;  
}
/******************************************/
/*                                        */
/* END FUNCTION: SetGlobalVariables       */
/*                                        */
/******************************************/

/******************************************/
/*                                        */
/* START FUNCTION: Freemem                */
/*                                        */
/******************************************/
int Freemem(StackSlideSearchParams *params)
{

  INT4 i;

/*  INITSTATUS( status, "Freemem", COMPUTESTACKSLIDESUMSC);
  ATTATCHSTATUSPTR (status); */

  /*Free BLKData*/
  /* 01/06/04 gam; Fix bug when deallocating BLKData and params->numBLKs != GV.SFTno. */
  /* for (i=0;i<GV.SFTno;i++) */ /* 01/06/04 gam */
  for (i=0;i<params->numBLKs;i++)
    {
      LALFree(params->BLKData[i]->fft->data->data);
      LALFree(params->BLKData[i]->fft->data);
      LALFree(params->BLKData[i]->fft);
      LALFree(params->BLKData[i]);
    }
  LALFree(params->BLKData);

  /*Free timeStamps*/
  LALFree(params->timeStamps);

  LALCheckMemoryLeaks();
  
/*  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status); */
  
  return 0;
}
/******************************************/
/*                                        */
/* END FUNCTION: Freemem                  */
/*                                        */
/******************************************/

/*********************************************************/
/*                                                       */
/* START FUNCTION: RunStackSlideMonteCarloSimulation     */
/* Injects fake signals and runs Monte Carlo simulation. */
/*                                                       */
/*********************************************************/
void RunStackSlideMonteCarloSimulation(LALStatus *status, StackSlideSearchParams *params)
{
#ifdef INCLUDE_RUNSTACKSLIDEMONTECARLO_CODE
  INT4    i = 0; /* all purpose index */
  INT4    j = 0; /* all purpose index */
  INT4    k = 0; /* all purpose index */  
  PulsarSignalParams *pPulsarSignalParams = NULL;
  REAL4TimeSeries *signal = NULL;
  LIGOTimeGPS GPSin;            /* reference-time for pulsar parameters at the detector; will convert to SSB! */
  LALDetector cachedDetector;
  REAL8 cosIota;                    /* cosine of inclination angle iota of the source */
  REAL8 h_0 = params->threshold4;   /* Source amplitude */
  SFTParams *pSFTParams = NULL;
  LIGOTimeGPSVector timestamps;
  SFTVector *outputSFTs = NULL;
  /* 07/14/04 gam; next two are needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
  SkyConstAndZeroPsiAMResponse *pSkyConstAndZeroPsiAMResponse;
  SFTandSignalParams *pSFTandSignalParams;
  INT4 iLastSkyPos; /* 07/14/04 gam; keep track of index to the last sky position */
  REAL4 renorm;               /* Need to renormalize SFTs to account for different sample rates */  
  FFT **savBLKData = NULL;    /* 05/21/04 gam; save the input noise SFTs */
  REAL8 **savSkyPosData;      /* 05/26/04 gam; save Sky Position Data */
  REAL8 **savFreqDerivData;   /* 05/26/04 gam; save Frequency Derivative Data */
  INT4 numSkyPosTotal = 0;    /* 05/26/04 gam; total number of Sky positions to cover */
  INT4 numFreqDerivTotal = 0; /* 05/26/04 gam; total number of Frequency evolution models to cover */
  INT4 numParamSpacePts = 0;  /* 05/26/04 gam; total number of points in the parameter space to cover */
  INT4 numSUMsTotal = 0;      /* 05/26/04 gam; total Number of Sums output = numSUMsPerParamPt*numParamSpacePts */
  INT4 numSUMsTotalm1 = -1;                 /* 05/26/04 gam */
  INT4 numFreqDerivIncludingNoSpinDown = 0; /* 05/26/04 gam */
  INT4 kSUM = 0;                            /* 05/26/04 gam */
  INT4 iFreq = 0;                           /* 05/26/04 gam */  
  REAL8 f0SUM = 0.0;                        /* 05/26/04 gam */
  REAL8 bandSUM = 0.0;                      /* 05/26/04 gam */
  INT4 nBinsPerSUM = 0;                     /* 05/26/04 gam */
  INT4 nBinsPerSUMm1 = -1;                  /* 05/26/04 gam */  
  INT4 seed=0; /* 05/28/04 gam; next 14 are for computing random mismatch */
  REAL4 randval;
  RandomParams *randPar=NULL;
  FILE *fpRandom;
  INT4 rndCount;
  REAL8 cosTmpDEC;
  REAL8 tmpDeltaRA;
  /* REAL8 DeltaRA = params->stksldSkyPatchData->deltaRA;
  REAL8 DeltaDec = params->stksldSkyPatchData->deltaDec;
  REAL8 DeltaFDeriv1 = params->deltaFDeriv1; */
  /* 06/01/04 gam; For now fix DeltaRA, DeltaDec,and DeltaFDeriv1 to default values when finding mismatch during Monte Carlo */
  REAL8 DeltaRA = 0.01;
  REAL8 DeltaDec = 0.01;
  REAL8 DeltaFDeriv1 = -1.0e-10;
  REAL8 DeltaFDeriv2 = params->deltaFDeriv2;
  REAL8 DeltaFDeriv3 = params->deltaFDeriv3;
  REAL8 DeltaFDeriv4 = params->deltaFDeriv4;
  REAL8 DeltaFDeriv5 = params->deltaFDeriv5;

  INITSTATUS( status, "RunStackSlideMonteCarloSimulation", COMPUTESTACKSLIDESUMSC );
  ATTATCHSTATUSPTR(status);
  
  /* 05/21/04 gam; save the input noise SFTs */
  savBLKData=(FFT **)LALMalloc(params->numBLKs*sizeof(FFT *));
  for (i=0;i<params->numBLKs;i++)
  {
        savBLKData[i]=(FFT *)LALMalloc(sizeof(FFT));
        savBLKData[i]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
        savBLKData[i]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
        savBLKData[i]->fft->data->data=(COMPLEX8 *)LALMalloc(params->nBinsPerBLK*sizeof(COMPLEX8));

        savBLKData[i]->fft->epoch = params->BLKData[i]->fft->epoch;
        savBLKData[i]->fft->f0 = params->BLKData[i]->fft->f0;
        savBLKData[i]->fft->deltaF = params->BLKData[i]->fft->deltaF;
        savBLKData[i]->fft->data->length = params->BLKData[i]->fft->data->length;
        for(j=0;j<params->nBinsPerBLK;j++) {
           savBLKData[i]->fft->data->data[j].re = params->BLKData[i]->fft->data->data[j].re;
           savBLKData[i]->fft->data->data[j].im = params->BLKData[i]->fft->data->data[j].im;
        }
  }  
 
  /* 05/21/04 gam; save the input skyPosData */
  savSkyPosData=(REAL8 **)LALMalloc(params->numSkyPosTotal*sizeof(REAL8 *));
  for(i=0;i<params->numSkyPosTotal;i++)
  {
        savSkyPosData[i] = (REAL8 *)LALMalloc(2*sizeof(REAL8));
        savSkyPosData[i][0] = params->skyPosData[i][0];
        savSkyPosData[i][1] = params->skyPosData[i][1];
  }
  /* reallocate memory for the skyPosData structure */
  for(i=0;i<params->numSkyPosTotal;i++)
  {
      LALFree(params->skyPosData[i]);
  }
  LALFree(params->skyPosData);
  numSkyPosTotal = params->numSkyPosTotal;  /* sav original number of sky positions */
  params->numSkyPosTotal = 1;               /* Set up to run on one sky position at a time */
  params->skyPosData=(REAL8 **)LALMalloc(params->numSkyPosTotal*sizeof(REAL8 *));
  params->skyPosData[0] = (REAL8 *)LALMalloc(2*sizeof(REAL8));
  
  /* 05/21/04 gam; save the input freqDerivData */
  if (params->numSpinDown > 0) {
    savFreqDerivData=(REAL8 **)LALMalloc(params->numFreqDerivTotal*sizeof(REAL8 *));
    for(i=0;i<params->numFreqDerivTotal;i++)
    {
        savFreqDerivData[i] = (REAL8 *)LALMalloc(params->numSpinDown*sizeof(REAL8));
        for(j=0;j<params->numSpinDown;j++)
        {
          savFreqDerivData[i][j] = params->freqDerivData[i][j];
        }
    }
    /* reallocate memory for the freqDerivData structure */
    for(i=0;i<params->numFreqDerivTotal;i++)
    {
        LALFree(params->freqDerivData[i]);
    }
    LALFree(params->freqDerivData);
    numFreqDerivTotal = params->numFreqDerivTotal;
    params->numFreqDerivTotal = 1;
    params->freqDerivData=(REAL8 **)LALMalloc(params->numFreqDerivTotal*sizeof(REAL8 *));
    params->freqDerivData[0] = (REAL8 *)LALMalloc(params->numSpinDown*sizeof(REAL8));
  }
  numParamSpacePts = params->numParamSpacePts;
  numSUMsTotal = params->numSUMsTotal;
  params->numParamSpacePts = 1;
  params->numSUMsTotal = 1;
  if (numFreqDerivTotal != 0) {
     numFreqDerivIncludingNoSpinDown = numFreqDerivTotal;
  } else {
     numFreqDerivIncludingNoSpinDown = 1;  /* Even if numSpinDown = 0 still need to count case of zero spindown. */
  }
   
  /* Allocate memory for PulsarSignalParams and initialize */
  pPulsarSignalParams = (PulsarSignalParams *)LALMalloc(sizeof(PulsarSignalParams));
  pPulsarSignalParams->pulsar.position.system = COORDINATESYSTEM_EQUATORIAL;
  pPulsarSignalParams->pulsar.spindown = NULL;
  if (params->numSpinDown > 0) {
    LALDCreateVector(status->statusPtr, &(pPulsarSignalParams->pulsar.spindown),((UINT4)params->numSpinDown));
    CHECKSTATUSPTR (status);
  }
  pPulsarSignalParams->orbit = NULL; 
  pPulsarSignalParams->transfer = NULL;  
  /* Set up pulsar site */
  if (strstr(params->IFO, "LHO")) {
       cachedDetector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  } else if (strstr(params->IFO, "LLO")) {
       cachedDetector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  } else if (strstr(params->IFO, "GEO")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  } else if (strstr(params->IFO, "VIRGO")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  } else if (strstr(params->IFO, "TAMA")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  } else {
      /* "Invalid or null IFO" */
      /* ABORT( status, DRIVESTACKSLIDEH_EIFO, DRIVESTACKSLIDEH_MSGEIFO); */
  }    
  pPulsarSignalParams->site = &cachedDetector;     
  pPulsarSignalParams->ephemerides = params->edat;
  /* pPulsarSignalParams->startTimeGPS = params->actualStartTime; */ /* 06/01/04 gam; uses gpsStartTime requested; not actualStartTime */
  pPulsarSignalParams->startTimeGPS.gpsSeconds = (INT4)params->gpsStartTimeSec;
  pPulsarSignalParams->startTimeGPS.gpsNanoSeconds = (INT4)params->gpsStartTimeNan;
  pPulsarSignalParams->duration = (UINT4)params->duration;
  pPulsarSignalParams->samplingRate = (REAL8)ceil(2.0*params->bandBLK); /* Make sampleRate an integer so that T*samplingRate = integer for integer T */
  pPulsarSignalParams->fHeterodyne = params->f0BLK;  
  /* Find the time at the SSB that corresponds to the arrive time at the detector of first data requested. */
  GPSin.gpsSeconds = (INT4)params->gpsEpochStartTimeSec;     /* 12/06/04 gam; GPS epoch Sec that gives reference time in SSB */
  GPSin.gpsNanoSeconds = (INT4)params->gpsEpochStartTimeNan; /* 12/06/04 gam; GPS epoch Nan that gives reference time in SSB */
  
  /* Allocate memory for SFTParams and initialize */
  pSFTParams = (SFTParams *)LALMalloc(sizeof(SFTParams));
  pSFTParams->Tsft = params->tBLK;
  timestamps.length = params->numBLKs;
  timestamps.data = params->timeStamps; 
  pSFTParams->timestamps = &timestamps;
  pSFTParams->noiseSFTs = NULL;  
  
  /* 05/26/04 gam; initialize variables that keep track of which SUM we are working on */
  params->startSUMs = 1;   /* 05/26/04 gam; use to control I/O during Monte Carlo. Default is TRUE. */
  params->finishSUMs = 0;  /* 05/26/04 gam; use to control I/O during Monte Carlo. Default is TRUE. */
  params->whichMCSUM = -1; /* 05/26/04 gam; which SUM the Monte Carlo Simulation is running on. Default is -1 */
  numSUMsTotalm1 = numSUMsTotal - 1;  /* Index of last SUM */
  
  /* 05/26/04 gam; initialize variables that keep track of which frequency we are working on */  
  f0SUM = params->f0SUM;
  bandSUM = params->bandSUM;
  nBinsPerSUM = params->nBinsPerSUM;  
  nBinsPerSUMm1 = nBinsPerSUM - 1;
  params->keepThisNumber = 1;      /* During MC only keep 1 event; this will be the injected event */
  
  /* check whether we are outputing just the loudest events */
  if ( ((params->outputEventFlag & 2) > 0) && (params->thresholdFlag <= 0) ) {
        params->bandSUM = params->dfSUM;
        params->nBinsPerSUM = 1;
  } else {
     /* TO DO: MAYBE SHOULD ABORT ?? */
  }
  
  /* 05/28/04 gam; Initial seed and randPar to use LALUniformDeviate to generate random mismatch during Monte Carlo. */
  fpRandom = fopen("/dev/urandom","r");
  rndCount = fread(&seed, sizeof(INT4),1, fpRandom);
  fclose(fpRandom);
  /* seed = 1234; */ /* Test value */
  LALCreateRandomParams(status->statusPtr, &randPar, seed);
  CHECKSTATUSPTR (status);
  
  /* 07/14/04 gam; allocate memory for structs needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
  if ( (params->testFlag & 4) > 0 ) {
     #ifdef DEBUG_LALFASTGENERATEPULSARSFTS
          fprintf(stdout,"allocate memory for structs needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs \n");
          fflush(stdout);
     #endif
     pSkyConstAndZeroPsiAMResponse = (SkyConstAndZeroPsiAMResponse *)LALMalloc(sizeof(SkyConstAndZeroPsiAMResponse));
     pSkyConstAndZeroPsiAMResponse->skyConst = (REAL8 *)LALMalloc((2*params->numSpinDown*(params->numBLKs+1)+2*params->numBLKs+3)*sizeof(REAL8));
     pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi = (REAL4 *)LALMalloc(params->numBLKs*sizeof(REAL4));
     pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi = (REAL4 *)LALMalloc(params->numBLKs*sizeof(REAL4));
     pSFTandSignalParams = (SFTandSignalParams *)LALMalloc(sizeof(SFTandSignalParams));
     /* create lookup table (LUT) values for doing trig */
     /* pSFTandSignalParams->resTrig = 64; */ /* length sinVal and cosVal; resolution of trig functions = 2pi/resTrig */
     pSFTandSignalParams->resTrig = 0; /* 08/02/04 gam; avoid serious bug when using LUTs for trig calls */
     pSFTandSignalParams->trigArg = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
     pSFTandSignalParams->sinVal  = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
     pSFTandSignalParams->cosVal  = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
     for (k=0; k<=pSFTandSignalParams->resTrig; k++) {
       pSFTandSignalParams->trigArg[k]= ((REAL8)LAL_TWOPI) * ((REAL8)k) / ((REAL8)pSFTandSignalParams->resTrig);
       pSFTandSignalParams->sinVal[k]=sin( pSFTandSignalParams->trigArg[k] );
       pSFTandSignalParams->cosVal[k]=cos( pSFTandSignalParams->trigArg[k] );
     }
     pSFTandSignalParams->pSigParams = pPulsarSignalParams;
     pSFTandSignalParams->pSFTParams = pSFTParams;
     pSFTandSignalParams->nSamples = GV.nsamples;
     pPulsarSignalParams->samplingRate = 2.0*params->bandBLK; /* can set samplingRate to exactly 2.0*bandBLK when using LALFastGeneratePulsarSFTs */
     GPSin.gpsSeconds = params->timeStamps[0].gpsSeconds;         /* 08/02/04 gam */
     GPSin.gpsNanoSeconds = params->timeStamps[0].gpsNanoSeconds; /* 08/02/04 gam */
  } else { 
     pSkyConstAndZeroPsiAMResponse = NULL;
     pSFTandSignalParams = NULL;
  }

  /*********************************************************/
  /*                                                       */
  /* START SECTION: MONTE CARLO LOOP OVER PARAMETER SPACE  */
  /*                                                       */
  /*********************************************************/
  iLastSkyPos = -1; /* 07/14/04 gam initialize */
  for(kSUM=0;kSUM<numSUMsTotal;kSUM++) {
  
    params->whichMCSUM = kSUM; /* keep track in StackSlideApplySearch of which SUM we are injecting into. */  
    
    /* 05/26/04 gam; set up pointer to current sky position and spindown parameters */
    i = kSUM/numFreqDerivIncludingNoSpinDown;   /* index to params->skyPosData for this SUM; */
    j = kSUM % numFreqDerivIncludingNoSpinDown; /* index to params->freqDerivData for this SUM; */
    params->skyPosData[0][0] = savSkyPosData[i][0];
    params->skyPosData[0][1] = savSkyPosData[i][1];        
    LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status); /* 05/28/04 gam */
    #ifdef DEBUG_SETFIXED_RANDVAL
       randval = params->threshold5; /* Temporarily use threshold5 for this; need to add testParameters to commandline. */
    #endif    
    pPulsarSignalParams->pulsar.position.latitude = params->skyPosData[0][1] + (((REAL8)randval) - 0.5)*DeltaDec;
    cosTmpDEC = cos(savSkyPosData[i][1]);
    if (cosTmpDEC != 0.0) {
          tmpDeltaRA = DeltaRA/cosTmpDEC;
    } else {
          tmpDeltaRA = 0.0; /* We are at a celestial pole */
    }
    LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status); /* 05/28/04 gam */
    #ifdef DEBUG_SETFIXED_RANDVAL
       randval = params->threshold5; /* Temporarily use threshold5 for this; need to add testParameters to commandline. */
    #endif    
    pPulsarSignalParams->pulsar.position.longitude = params->skyPosData[0][0] + (((REAL8)randval) - 0.5)*tmpDeltaRA;
    /* pPulsarSignalParams->pulsar.position.longitude = params->skyPosData[0][0];
    pPulsarSignalParams->pulsar.position.latitude = params->skyPosData[0][1]; */
    #ifdef DEBUG_RANDOMTRIALPARAMETERS
        fprintf(stdout,"kSUM = %i, search RA = %23.10e, inject RA = %23.10e \n",kSUM,params->skyPosData[0][0],pPulsarSignalParams->pulsar.position.longitude);
        fprintf(stdout,"kSUM = %i, search DEC = %23.10e, inject DEC = %23.10e \n",kSUM,params->skyPosData[0][1],pPulsarSignalParams->pulsar.position.latitude);
        fflush(stdout);      
    #endif

    if (params->numSpinDown > 0) {
      for(k=0;k<params->numSpinDown;k++) {
        params->freqDerivData[0][k] = savFreqDerivData[j][k];
        LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status); /* 05/28/04 gam */
        #ifdef DEBUG_SETFIXED_RANDVAL
           randval = params->threshold5; /* Temporarily use threshold5 for this; need to add testParameters to commandline. */
        #endif
        if (k == 0) {
          pPulsarSignalParams->pulsar.spindown->data[k] = params->freqDerivData[0][k] + (((REAL8)randval) - 0.5)*DeltaFDeriv1;
        } else if (k == 1) {
          pPulsarSignalParams->pulsar.spindown->data[k] = params->freqDerivData[0][k] + (((REAL8)randval) - 0.5)*DeltaFDeriv2;
        } else if (k == 2) {
          pPulsarSignalParams->pulsar.spindown->data[k] = params->freqDerivData[0][k] + (((REAL8)randval) - 0.5)*DeltaFDeriv3;
        } else if (k == 3) {
          pPulsarSignalParams->pulsar.spindown->data[k] = params->freqDerivData[0][k] + (((REAL8)randval) - 0.5)*DeltaFDeriv4;
        } else if (k == 4) {
          pPulsarSignalParams->pulsar.spindown->data[k] = params->freqDerivData[0][k] + (((REAL8)randval) - 0.5)*DeltaFDeriv5;
        } /* END if (k == 0) ELSE ... */
        /* pPulsarSignalParams->pulsar.spindown->data[k] = params->freqDerivData[0][k]; */
        #ifdef DEBUG_RANDOMTRIALPARAMETERS
          fprintf(stdout,"kSUM = %i, search fDeriv[%i] = %23.10e, inject fDeriv[%i] = %23.10e \n",kSUM,k,params->freqDerivData[0][k],k,pPulsarSignalParams->pulsar.spindown->data[k]);
          fflush(stdout);      
        #endif
      }
    }
    
    /* 07/14/04 gam; one per sky position update the reference time for this sky position */
    if ( iLastSkyPos != i ) {
      LALConvertGPS2SSB(status->statusPtr,&(pPulsarSignalParams->pulsar.tRef), GPSin, pPulsarSignalParams);
      CHECKSTATUSPTR (status);
      
      /* 07/14/04 gam; one per sky position fill in SkyConstAndZeroPsiAMResponse for use with LALFastGeneratePulsarSFTs */
      if ( (params->testFlag & 4) > 0 ) {
        #ifdef DEBUG_LALFASTGENERATEPULSARSFTS
            fprintf(stdout,"about to call LALComputeSkyAndZeroPsiAMResponse \n");
            fflush(stdout);
        #endif
        LALComputeSkyAndZeroPsiAMResponse (status->statusPtr, pSkyConstAndZeroPsiAMResponse, pSFTandSignalParams);
        CHECKSTATUSPTR (status);
      }      
    }    
    iLastSkyPos = i; /* Update the index to the last Sky Position */
    
    /****************************************************/
    /*                                                  */
    /* START SECTION: MONTE CARLO LOOP OVER FREQUENCIES */
    /*                                                  */
    /****************************************************/    
    for(iFreq=0;iFreq<nBinsPerSUM;iFreq++) {
    
      if (iFreq > 0) {
        params->startSUMs = 0;  /* Indicate that SUMs already started */
      }
    
      if ((kSUM == numSUMsTotalm1) && (iFreq == nBinsPerSUMm1)) {
          params->finishSUMs = 1;  /* This is the last injection */
      }

      /* 12/06/04 gam */
      if ( (params->testFlag & 8) > 0 ) {
         pPulsarSignalParams->pulsar.psi = params->orientationAngle;
         cosIota =params->cosInclinationAngle;
      } else {
         /* get random value for psi */
         LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status); /* 05/28/04 gam */
         #ifdef DEBUG_SETFIXED_RANDVAL
            randval = params->threshold5; /* Temporarily use threshold5 for this; need to add testParameters to commandline. */
         #endif
         pPulsarSignalParams->pulsar.psi = (randval - 0.5) * ((REAL4)LAL_PI_2);
         /* pPulsarSignalParams->pulsar.psi = 0.0; */
    
         /* get random value for cosIota */
         LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status); /* 05/28/04 gam */
         #ifdef DEBUG_SETFIXED_RANDVAL
            randval = params->threshold5; /* Temporarily use threshold5 for this; need to add testParameters to commandline. */
         #endif
         cosIota = 2.0*((REAL8)randval) - 1.0;
         /* cosIota = 1.0; */
      }
      
      /* h_0 is fixed equal to params->threshold4 above; get A_+ and A_x from h_0 and random cosIota */
      pPulsarSignalParams->pulsar.aPlus = (REAL4)(0.5*h_0*(1.0 + cosIota*cosIota));
      pPulsarSignalParams->pulsar.aCross = (REAL4)(h_0*cosIota);

      /* get random value for phi0 */
      LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status); /* 05/28/04 gam */
      #ifdef DEBUG_SETFIXED_RANDVAL
         randval = params->threshold5; /* Temporarily use threshold5 for this; need to add testParameters to commandline. */
      #endif
      pPulsarSignalParams->pulsar.phi0 = ((REAL8)randval) * ((REAL8)LAL_TWOPI);
      /* pPulsarSignalParams->pulsar.phi0 = 0.0; */

      pPulsarSignalParams->pulsar.f0 = f0SUM + iFreq*params->dfSUM; /* We add in the mismatch after adjusting params->f0SUM */
      /* check whether we are outputing just the loudest events */
      if ( ((params->outputEventFlag & 2) > 0) && (params->thresholdFlag <= 0) ) {
        params->f0SUM = pPulsarSignalParams->pulsar.f0;
      } else {
        /* TO DO: MAYBE SHOULD ABORT ABOVE AND ALWAY SET params->f0SUM as above ?? */
      }
      /* Now add mismatch to pPulsarSignalParams->pulsar.f0 */
      LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status); /* 05/28/04 gam */
      #ifdef DEBUG_SETFIXED_RANDVAL
         randval = params->threshold5; /* Temporarily use threshold5 for this; need to add testParameters to commandline. */
      #endif      
      pPulsarSignalParams->pulsar.f0 = pPulsarSignalParams->pulsar.f0 + (((REAL8)randval) - 0.5)*params->dfSUM;
      
      #ifdef DEBUG_RANDOMTRIALPARAMETERS
        fprintf(stdout,"iFreq = %i, inject h_0 = %23.10e \n",iFreq,h_0);
        fprintf(stdout,"iFreq = %i, inject cosIota = %23.10e, A_+ = %23.10e, A_x = %23.10e \n",iFreq,cosIota,pPulsarSignalParams->pulsar.aPlus,pPulsarSignalParams->pulsar.aCross);
        fprintf(stdout,"iFreq = %i, inject psi = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.psi);
        fprintf(stdout,"iFreq = %i, inject phi0 = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.phi0);
        fprintf(stdout,"iFreq = %i, search f0 = %23.10e, inject f0 = %23.10e \n",iFreq,f0SUM + iFreq*params->dfSUM,pPulsarSignalParams->pulsar.f0);
        fflush(stdout);
      #endif
      #ifdef DEBUG_MONTECARLOTIMEDOMAIN_DATA
        fprintf(stdout,"iFreq = %i, params->f0SUM = %23.10e, pPulsarSignalParams->pulsar.f0 = %23.10e \n",iFreq,params->f0SUM,pPulsarSignalParams->pulsar.f0);
        fprintf(stdout,"nBinsPerSUM = %i, params->nBinsPerSUM = %i\n",nBinsPerSUM,params->nBinsPerSUM);
        fflush(stdout);	
      #endif
      
      /* 07/14/04 gam; check if using LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
      if ( (params->testFlag & 4) > 0 ) {
        #ifdef DEBUG_LALFASTGENERATEPULSARSFTS
          fprintf(stdout,"About to call LALFastGeneratePulsarSFTs \n");
          fflush(stdout);
        #endif
        /* 07/14/04 gam; use SkyConstAndZeroPsiAMResponse from LALComputeSkyAndZeroPsiAMResponse and SFTandSignalParams to generate SFTs fast. */
        LALFastGeneratePulsarSFTs (status->statusPtr, &outputSFTs, pSkyConstAndZeroPsiAMResponse, pSFTandSignalParams);
        CHECKSTATUSPTR (status);
        
        #ifdef PRINT_ONEMONTECARLO_OUTPUTSFT
          if (params->whichMCSUM == 0 && iFreq == 0) {
            REAL4  fPlus;
            REAL4  fCross;
            i=0; /* index of which outputSFT to output */
            fPlus = pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi[i]*cos(pPulsarSignalParams->pulsar.psi) + pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi[i]*sin(pPulsarSignalParams->pulsar.psi);
            fCross = pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi[i]*cos(pPulsarSignalParams->pulsar.psi) - pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi[i]*sin(pPulsarSignalParams->pulsar.psi);
            fprintf(stdout,"iFreq = %i, inject h_0 = %23.10e \n",iFreq,h_0);
            fprintf(stdout,"iFreq = %i, inject cosIota = %23.10e, A_+ = %23.10e, A_x = %23.10e \n",iFreq,cosIota,pPulsarSignalParams->pulsar.aPlus,pPulsarSignalParams->pulsar.aCross);
            fprintf(stdout,"iFreq = %i, inject psi = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.psi);
            fprintf(stdout,"iFreq = %i, fPlus, fCross = %23.10e,  %23.10e \n",iFreq,fPlus,fCross);
            fprintf(stdout,"iFreq = %i, inject phi0 = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.phi0);
            fprintf(stdout,"iFreq = %i, search f0 = %23.10e, inject f0 = %23.10e \n",iFreq,f0SUM + iFreq*params->dfSUM,pPulsarSignalParams->pulsar.f0);
            fprintf(stdout,"outputSFTs->data[%i].data->length = %i \n",i,outputSFTs->data[i].data->length);
            fflush(stdout);
            for(j=0;j<params->nBinsPerBLK;j++) {
               /* fprintf(stdout,"%i %g %g \n",j,outputSFTs->data[i].data->data[j].re,outputSFTs->data[i].data->data[j].im); */
               fprintf(stdout,"%i %g \n",j,outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im);
               fflush(stdout);
            }
          }
        #endif
        
        #ifdef PRINT_MAXPOWERANDBINEACHSFT
         {
           INT4 jMaxPwr;
           REAL4 maxPwr;
           REAL4 pwr;
           for(i=0;i<params->numBLKs;i++) {
             maxPwr = 0.0;
             jMaxPwr = -1;
             for(j=0;j<params->nBinsPerBLK;j++) {
               pwr = outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im;
               if (pwr > maxPwr) {
                  maxPwr = pwr;
                  jMaxPwr = j;
               }
             }
             fprintf(stdout,"Max power for SFT %i is in bin %i = %g \n",i,jMaxPwr,maxPwr);
             fflush(stdout);
           }
         }
        #endif

        /* Add outputSFTs with injected signal to input noise SFTs; no renorm should be needed */
        for(i=0;i<params->numBLKs;i++) {
          for(j=0;j<params->nBinsPerBLK;j++) {
             #ifdef PRINTCOMPARISON_INPUTVSMONTECARLOSFT_DATA
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].re, outputSFTs->data[%i].data->data[%i].re = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[j].re,outputSFTs->data[i].data->data[j].re);
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].im, outputSFTs->data[%i].data->data[%i].im = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[j].im,outputSFTs->data[i].data->data[j].im);
               fflush(stdout);
             #endif
             params->BLKData[i]->fft->data->data[j].re = savBLKData[i]->fft->data->data[j].re + outputSFTs->data[i].data->data[j].re;
             params->BLKData[i]->fft->data->data[j].im = savBLKData[i]->fft->data->data[j].im + outputSFTs->data[i].data->data[j].im;
          }
        }
        /* next is else for if ( (params->testFlag & 4) > 0 ) else... */
      } else {
        #ifdef PRINT_ONEMONTECARLO_OUTPUTSFT
          if (params->whichMCSUM == 0 && iFreq == 0) {
            /* To compare with LALFastGeneratePulsarSFTs, which uses ComputeSky, use first time stamp as reference time */
            GPSin.gpsSeconds = params->timeStamps[0].gpsSeconds;
            GPSin.gpsNanoSeconds = params->timeStamps[0].gpsNanoSeconds;
            LALConvertGPS2SSB(status->statusPtr,&(pPulsarSignalParams->pulsar.tRef), GPSin, pPulsarSignalParams); CHECKSTATUSPTR (status);
          }
        #endif

        signal = NULL; /* Call LALGeneratePulsarSignal to generate signal */
        LALGeneratePulsarSignal(status->statusPtr, &signal, pPulsarSignalParams);
        CHECKSTATUSPTR (status);

        #ifdef DEBUG_MONTECARLOTIMEDOMAIN_DATA
          fprintf(stdout,"signal->deltaT = %23.10e \n",signal->deltaT);
          fprintf(stdout,"signal->epoch.gpsSeconds = %i \n",signal->epoch.gpsSeconds);
          fprintf(stdout,"signal->epoch.gpsNanoSeconds = %i \n",signal->epoch.gpsNanoSeconds);
          fprintf(stdout,"pPulsarSignalParams->duration*pPulsarSignalParams->samplingRate, signal->data->length = %g %i\n",pPulsarSignalParams->duration*pPulsarSignalParams->samplingRate,signal->data->length);
          fflush(stdout);
        #endif
        #ifdef PRINT_MONTECARLOTIMEDOMAIN_DATA
          for(i=0;i<signal->data->length;i++)  {
            fprintf(stdout,"signal->data->data[%i] = %g\n",i,signal->data->data[i]);
            fflush(stdout);
          }
        #endif

        outputSFTs = NULL; /* Call LALSignalToSFTs to generate outputSFTs with injection data */
        LALSignalToSFTs(status->statusPtr, &outputSFTs, signal, pSFTParams);
        CHECKSTATUSPTR (status);

        #ifdef PRINT_ONEMONTECARLO_OUTPUTSFT
           if (params->whichMCSUM == 0 && iFreq == 0) {
              REAL4 realHetPhaseCorr;
              REAL4 imagHetPhaseCorr;
              REAL8 deltaTHeterodyne;
              i=0; /* index of which outputSFT to output */
              fprintf(stdout,"iFreq = %i, inject h_0 = %23.10e \n",iFreq,h_0);
              fprintf(stdout,"iFreq = %i, inject cosIota = %23.10e, A_+ = %23.10e, A_x = %23.10e \n",iFreq,cosIota,pPulsarSignalParams->pulsar.aPlus,pPulsarSignalParams->pulsar.aCross);
              fprintf(stdout,"iFreq = %i, inject psi = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.psi);
              fprintf(stdout,"iFreq = %i, inject phi0 = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.phi0);
              fprintf(stdout,"iFreq = %i, search f0 = %23.10e, inject f0 = %23.10e \n",iFreq,f0SUM + iFreq*params->dfSUM,pPulsarSignalParams->pulsar.f0);
              fprintf(stdout,"outputSFTs->data[%i].data->length = %i \n",i,outputSFTs->data[i].data->length);
              /* deltaTHeterodyne = the time from the refererence time to the start of an SFT */
              GPSin.gpsSeconds = params->timeStamps[i].gpsSeconds;
              GPSin.gpsNanoSeconds = params->timeStamps[i].gpsNanoSeconds;
              TRY (LALDeltaFloatGPS (status->statusPtr, &deltaTHeterodyne, &GPSin, &(pPulsarSignalParams->pulsar.tRef)), status);
              /* deltaTHeterodyne = 0.0; */ /* UNCOMMENT TO TURN OFF CORRECTION OF PHASE */
              fprintf(stdout,"t_h = %23.10e\n",deltaTHeterodyne);
              fflush(stdout);
              realHetPhaseCorr = cos(LAL_TWOPI*pPulsarSignalParams->fHeterodyne*deltaTHeterodyne);
              imagHetPhaseCorr = sin(LAL_TWOPI*pPulsarSignalParams->fHeterodyne*deltaTHeterodyne);
              renorm = ((REAL4)GV.nsamples)/((REAL4)(outputSFTs->data[i].data->length - 1));
              for(j=0;j<params->nBinsPerBLK;j++) {
                /* fprintf(stdout,"%i %g %g \n",j,renorm*outputSFTs->data[i].data->data[j].re*realHetPhaseCorr - renorm*outputSFTs->data[i].data->data[j].im*imagHetPhaseCorr,renorm*outputSFTs->data[i].data->data[j].re*imagHetPhaseCorr + renorm*outputSFTs->data[i].data->data[j].im*realHetPhaseCorr); */
                fprintf(stdout,"%i %g \n",j,renorm*renorm*outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + renorm*renorm*outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im);
                fflush(stdout);
              }
           }
        #endif

        #ifdef PRINT_MAXPOWERANDBINEACHSFT
         {
           INT4 jMaxPwr;
           REAL4 maxPwr;
           REAL4 pwr;
           for(i=0;i<params->numBLKs;i++) {
             maxPwr = 0.0;
             jMaxPwr = -1;
             renorm = ((REAL4)GV.nsamples)/((REAL4)(outputSFTs->data[i].data->length - 1));
             for(j=0;j<params->nBinsPerBLK;j++) {
               pwr = outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im;
               pwr = renorm*renorm*pwr;
               if (pwr > maxPwr) {
                  maxPwr = pwr;
                  jMaxPwr = j;
               }
             }
             fprintf(stdout,"Max power for SFT %i is in bin %i = %g \n",i,jMaxPwr,maxPwr);
             fflush(stdout);
           }
         }
        #endif
        
        /* Add outputSFTs with injected signal to input noise SFTs; renorm is needed. */
        /* 05/21/04 gam; Normalize the SFTs from LALSignalToSFTs as per the normalization in makefakedata_v2.c and lalapps/src/pulsar/make_sfts.c. */
        renorm = ((REAL4)GV.nsamples)/((REAL4)(outputSFTs->data[0].data->length - 1)); /* 05/21/04 gam; should be the same for all outputSFTs; note minus 1 is needed */
        for(i=0;i<params->numBLKs;i++) {
          /* renorm = ((REAL4)GV.nsamples)/((REAL4)outputSFTs->data[i].data->length); */ /* 05/21/04 gam; do once correctly above */ /* Should be the same for all SFTs, but just in case recompute. */
          #ifdef DEBUG_MONTECARLOSFT_DATA  
            fprintf(stdout,"Mutiplying outputSFTs->data[%i] with renorm = %g \n",i,renorm);
            fflush(stdout);
          #endif  
          for(j=0;j<params->nBinsPerBLK;j++) {
             #ifdef PRINTCOMPARISON_INPUTVSMONTECARLOSFT_DATA
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].re, outputSFTs->data[%i].data->data[%i].re = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[j].re,renorm*outputSFTs->data[i].data->data[j].re);
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].im, outputSFTs->data[%i].data->data[%i].im = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[j].im,renorm*outputSFTs->data[i].data->data[j].im);
               fflush(stdout);
             #endif
             /* 05/21/04 gam; changed next two lines so that params->BLKData is sum of saved noise SFTs and SFTs with injection data. */
             params->BLKData[i]->fft->data->data[j].re = savBLKData[i]->fft->data->data[j].re + renorm*outputSFTs->data[i].data->data[j].re;
             params->BLKData[i]->fft->data->data[j].im = savBLKData[i]->fft->data->data[j].im + renorm*outputSFTs->data[i].data->data[j].im;
          }
        }
      } /* END if ( (params->testFlag & 4) > 0 ) else */

      #ifdef DEBUG_MONTECARLOSFT_DATA  
        fprintf(stdout,"params->numBLKs, outputSFTs->length = %i, %i \n",params->numBLKs,outputSFTs->length);
        fflush(stdout);
        for(i=0;i<outputSFTs->length;i++) {
          fprintf(stdout,"params->f0BLK, outputSFTs->data[%i].f0 = %g, %g\n",i,params->f0BLK,outputSFTs->data[i].f0);
          fprintf(stdout,"params->tBLK, 1.0/outputSFTs->data[%i].deltaF = %g, %g\n",i,params->tBLK,1.0/outputSFTs->data[i].deltaF);
          fprintf(stdout,"params->timeStamps[%i].gpsSeconds, outputSFTs->data[%i].epoch.gpsSeconds = %i, %i\n",i,i,params->timeStamps[i].gpsSeconds,outputSFTs->data[i].epoch.gpsSeconds);
          fprintf(stdout,"params->timeStamps[%i].gpsNanoSeconds, outputSFTs->data[%i].epoch.gpsNanoSeconds = %i, %i\n",i,i,params->timeStamps[i].gpsNanoSeconds,outputSFTs->data[i].epoch.gpsNanoSeconds);
          fprintf(stdout,"params->nBinsPerBLK, outputSFTs->data[%i].data->length = %i, %i\n",i,params->nBinsPerBLK,outputSFTs->data[i].data->length);
          fflush(stdout);
        }
      #endif
      #ifdef PRINT_MONTECARLOSFT_DATA
        for(i=0;i<outputSFTs->length;i++) {
          for(j=0;j<outputSFTs->data[i].data->length;j++) {
            fprintf(stdout,"outputSFTs->data[%i].data->data[%i].re, outputSFTs->data[%i].data->data[%i].im = %g, %g\n",i,j,i,j,outputSFTs->data[i].data->data[j].re,outputSFTs->data[i].data->data[j].im);
            fflush(stdout);
          }
        }
      #endif
      
      /* 07/14/04 gam; check if using LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
      if ( !((params->testFlag & 4) > 0) ) {
        LALDestroySFTVector(status->statusPtr, &outputSFTs);
        CHECKSTATUSPTR (status); /* 06/01/04 gam */
      
        LALFree(signal->data->data);
        LALFree(signal->data);
        LALFree(signal);
      }
      
      /* 07/14/04 gam; Example code */  /*
      SkyConstAndZeroPsiAMResponse *pSkyConstAndZeroPsiAMResponse;
      SFTandSignalParams *pSFTandSignalParams;

      pSkyConstAndZeroPsiAMResponse = (SkyConstAndZeroPsiAMResponse *)LALMalloc(sizeof(SkyConstAndZeroPsiAMResponse));
      pSkyConstAndZeroPsiAMResponse->skyConst = (REAL8 *)LALMalloc((2*params->numSpinDown*(params->numBLKs+1)+2*params->numBLKs+3)*sizeof(REAL8));
      pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi = (REAL4 *)LALMalloc(params->numBLKs*sizeof(REAL4));
      pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi = (REAL4 *)LALMalloc(params->numBLKs*sizeof(REAL4));
      pSFTandSignalParams = (SFTandSignalParams *)LALMalloc(sizeof(SFTandSignalParams));
      pSFTandSignalParams->resTrig = 64;
      pSFTandSignalParams->trigArg = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
      pSFTandSignalParams->sinVal  = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
      pSFTandSignalParams->cosVal  = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
      for (k=0; k<=pSFTandSignalParams->resTrig; k++) {
       pSFTandSignalParams->trigArg[k]= ((REAL8)LAL_TWOPI) * ((REAL8)k) / ((REAL8)pSFTandSignalParams->resTrig);
       pSFTandSignalParams->sinVal[k]=sin( pSFTandSignalParams->trigArg[k] );
       pSFTandSignalParams->cosVal[k]=cos( pSFTandSignalParams->trigArg[k] );
      }
      pSFTandSignalParams->pSigParams = pPulsarSignalParams;
      pSFTandSignalParams->pSFTParams = pSFTParams;
      pSFTandSignalParams->nSamples = GV.nsamples;
           
      LALComputeSkyAndZeroPsiAMResponse (status->statusPtr, pSkyConstAndZeroPsiAMResponse, pSFTandSignalParams);
      CHECKSTATUSPTR (status);
      LALFastGeneratePulsarSFTs (status->statusPtr, &outputSFTs, pSkyConstAndZeroPsiAMResponse, pSFTandSignalParams);
      CHECKSTATUSPTR (status);
      
      LALFree(pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi);
      LALFree(pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi);
      LALFree(pSkyConstAndZeroPsiAMResponse->skyConst);
      LALFree(pSkyConstAndZeroPsiAMResponse);
      LALFree(pSFTandSignalParams->trigArg);
      LALFree(pSFTandSignalParams->sinVal);
      LALFree(pSFTandSignalParams->cosVal);
      LALFree(pSFTandSignalParams);
      */

      StackSlideApplySearch(status->statusPtr,params);
      CHECKSTATUSPTR (status);
        
    } /* END for(iFreq=0;iFreq<nBinsPerSUM;iFreq++) */
    /****************************************************/
    /*                                                  */
    /* END SECTION: MONTE CARLO LOOP OVER FREQUENCIES   */
    /*                                                  */
    /****************************************************/
  } /* END for(kSUM=0;kSUM<numSUMsTotal;kSUM++) */
  /*********************************************************/
  /*                                                       */
  /* END SECTION: MONTE CARLO LOOP OVER PARAMETER SPACE    */
  /*                                                       */
  /*********************************************************/

  LALDestroyRandomParams(status->statusPtr, &randPar); /* 05/28/04 gam */
  CHECKSTATUSPTR (status);
  
  LALFree(pSFTParams);
      
  if (params->numSpinDown > 0) {
    LALDDestroyVector(status->statusPtr, &(pPulsarSignalParams->pulsar.spindown));
    CHECKSTATUSPTR (status);
  }
  LALFree(pPulsarSignalParams);

  /* 07/14/04 gam; deallocate memory for structs needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
  if ( (params->testFlag & 4) > 0 ) {
    #ifdef DEBUG_LALFASTGENERATEPULSARSFTS
      fprintf(stdout,"deallocate memory for structs needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs \n");
      fflush(stdout);
    #endif
    LALDestroySFTVector(status->statusPtr, &outputSFTs);
    CHECKSTATUSPTR (status);
    LALFree(pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi);
    LALFree(pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi);
    LALFree(pSkyConstAndZeroPsiAMResponse->skyConst);
    LALFree(pSkyConstAndZeroPsiAMResponse);
    LALFree(pSFTandSignalParams->trigArg);
    LALFree(pSFTandSignalParams->sinVal);
    LALFree(pSFTandSignalParams->cosVal);
    LALFree(pSFTandSignalParams);
  }
          
  /* 05/21/04 gam; free memory of saved data */
  for (i=0;i<params->numBLKs;i++) {
      LALFree(savBLKData[i]->fft->data->data);
      LALFree(savBLKData[i]->fft->data);
      LALFree(savBLKData[i]->fft);
      LALFree(savBLKData[i]);
  }
  LALFree(savBLKData);
  
  /* deallocate memory for the savSkyPosData structure */
  for(i=0;i<numSkyPosTotal;i++)
  {
      LALFree(savSkyPosData[i]);
  }
  LALFree(savSkyPosData);
  
  if (params->numSpinDown > 0) {
    /* deallocate memory for the savFreqDerivData structure */
    for(i=0;i<numFreqDerivTotal;i++)
    {
        LALFree(savFreqDerivData[i]);
    }
    LALFree(savFreqDerivData);
  }
  
  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);
#endif  
}
/*********************************************************/
/*                                                       */
/* END FUNCTION: RunStackSlideMonteCarloSimulation       */
/* Injects fake signals and runs Monte Carlo simulation. */
/*                                                       */
/*********************************************************/
