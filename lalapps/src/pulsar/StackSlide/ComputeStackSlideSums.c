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

/* REVISIONS: */
/* 12/03/03 gam; Based code on StackSlide.c from LALWrapper DSO */
/* 01/06/04 gam; Fix bug when deallocating BLKData and params->numBLKs != GV.SFTno. */
/* 01/23/04 gam; Glob for *SFT.* rather than *SFT* to avoid non-SFT files that have SFT in their name. */
/* 01/28/04 gam; Include what to glob for in params->sftDirectory rather than hard coding */
/* 04/15/04 gam; Fix help when no command line arguments */

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
    
  StackSlideApplySearch(&status,params);
  INTERNAL_CHECKSTATUS_FROMMAIN(status)    
  
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
  glob_t globbuf;
  
  /* INITSTATUS( status, "ReadSFTData", COMPUTESTACKSLIDESUMSC);
  ATTATCHSTATUSPTR (status); */

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

  GV.SFTno=fileno; /* remember this is 1 more than the index value */

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
