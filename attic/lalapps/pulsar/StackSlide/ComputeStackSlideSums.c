/*
*  Copyright (C) 2007 Gregory Mendell, Virginia Re
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

/**
 * \file
 * \ingroup pulsarApps
 * \author G.Mendell, M.Landry
 * \brief
 * Drive StackSlide jobs under Condor and the Grid
 */

/*********************************************************************************/
/*                                                                               */
/* File: ComputeStackSlideSums.c                                                 */
/* Purpose: Drive StackSlide jobs under Condor and the Grid                      */
/* Author: Mendell, G. and Landry, M.                                            */
/* Started: 03 December 2003                                                     */
/*                                                                               */
/*********************************************************************************/
/* Notes: */
/* BLK = one Block of frequency domain data */
/* STK = one Stack of frequency domain data; one or more Blocks are combined to make Stacks. */
/* SUM = Sum of Stacks after Sliding */

/* REVISIONS: */
/* 12/03/03 gam; Based code on StackSlide.c from LALWrapper DSO */
/* 01/06/04 gam; Fix bug when deallocating BLKData and params->numBLKs != GV.SFTno. */
/* 01/23/04 gam; Glob for *SFT.* rather than *SFT* to avoid non-SFT files that have SFT in their name. */
/* 01/28/04 gam; Include what to glob for in params->sftDirectory rather than hard coding */
/* 04/15/04 gam; Fix help when no command line arguments */
/* 05/07/04 gam; add alternative to using glob */
/* 05/07/04 gam; Do not check d_type when using readdir since it is often UNKOWN. */
/* 04/12/05 gam; Change default Monte Carlo to use RunStackSlideIsolatedMonteCarloSimulation in StackSlideIsolated.c */
/* 05/25/05 gam; remove obsolete inputDataTypeFlag and other obsolete code and comments */
/* 08/31/05 gam; Sort the list of SFT files: GV.filelist. The ideas is that sorting by filename */
/*               == sorting by GPS time for properly named SFT files! */
/* 08/31/05 gam; For clarity, change the condition when reading in BLKs of data from until blkno >= params->numBLKs to until blkno == params->numBLKs. */
/* 08/31/05 gam; Fix a few misleading comments; remove old comments; add new ones where clarification is needed */
/* 08/31/05 gam; Change stoptime to require input data to exist strickly in the interval [params->gpsStartTimeSec, params->gpsStartTimeSec + params->duration) */

/*********************************************/
/*                                           */
/* START SECTION: define preprocessor flags  */
/*                                           */
/*********************************************/
/* #define DEBUG_READSFT_CODE */
/* #define DEBUG_FILELIST_SORTING */
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
  StackSlideSearchParams *params;  /* This is the main container for all the parameters used to run a search */

  if (argc < 2) {
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

  if (ReadSFTData(params)) return 3; /* Note that for now BLK data is always SFT data */

  StackSlideConditionData(&status,params);  
  INTERNAL_CHECKSTATUS_FROMMAIN(status)  

  if ( (params->testFlag & 2) > 0 ) {
    /* 05/11/04 gam; Add code to software inject signals into the SFTs for Monte Carlo simulations */    
    /* RunStackSlideMonteCarloSimulation(&status,params); */ /* 04/12/05 gam; moved to StackSlideIsolated.c */

	 	  
/*!05/09/07 */	  if ( params->binaryFlag == 1 ){
/*!*/        RunStackSlideBinaryMonteCarloSimulation(&status,params,GV.nsamples);
/*!*/        INTERNAL_CHECKSTATUS_FROMMAIN(status);    

/*!*/               }else{
              RunStackSlideIsolatedMonteCarloSimulation(&status,params,GV.nsamples);
              INTERNAL_CHECKSTATUS_FROMMAIN(status);    
/*!*/	        }
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
 
  /* Allocate memory for BLK data */  
  params->BLKData=(FFT **)LALMalloc(params->numBLKs*sizeof(FFT *));
  params->timeStamps=(LIGOTimeGPS *)LALMalloc(params->numBLKs*sizeof(LIGOTimeGPS));
  /* params->timeIntervals =(gpsTimeInterval *)LALCalloc(params->numBLKs, sizeof(gpsTimeInterval)); */

  ndeltaf=GV.ifmax-GV.ifmin+1; 
  /* stoptime = params->gpsStartTimeSec + (UINT4)floor(params->duration); */ /* 08/31/05 gam */
  stoptime = params->gpsStartTimeSec + ((UINT4)floor(params->duration)) - ((UINT4)ceil(params->tBLK));

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
      /* if (((UINT4)header.gps_sec >= params->gpsStartTimeSec) && ((UINT4)header.gps_sec < stoptime)) */ /* 08/31/05 gam */
      if (((UINT4)header.gps_sec >= params->gpsStartTimeSec) && ((UINT4)header.gps_sec <= stoptime)) {
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
     
    } /* end if ((header.gps_sec >= params->gpsStartTimeSec) && (header.gps_sec <= stoptime)) */

    fclose(fp);     /* close file */

    /* if (blkno >= params->numBLKs) */ /* 08/31/05 gam; this will not change anything, but make it clear blkno is precisly == params->numBLKs */
    if (blkno == params->numBLKs) {
        params->finishedBLKs = 1;
        break;
    }

  } /* end for (fileno=0;fileno<GV.SFTno;fileno++) */

  /* Some requested input BLK data is missing */
  if (!params->finishedBLKs) {
     fprintf(stderr,"Some requested input BLK data is missing\n");
     return 6;     
  }

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

  /* 08/31/05 gam */  
  /***********************************************************/
  /*                                                         */
  /* START SECTION: sort the list of SFT files: GV.filelist. */
  /*                                                         */
  /***********************************************************/
   #ifdef DEBUG_FILELIST_SORTING  
     for (fileno=0;fileno<GV.SFTno;fileno++) {
       fprintf(stdout,"Before sorting GV.filelist[%i] = %s\n",fileno,GV.filelist[fileno]);
       fflush(stdout);
     }
     fileno = GV.SFTno; /* restore value */
   #endif
  
   qsort(GV.filelist,GV.SFTno,sizeof(GV.filelist[0]), (__compar_fn_t)strcmp); /* 08/31/05 gam; sort the filelist! */

   #ifdef DEBUG_FILELIST_SORTING  
     for (fileno=0;fileno<GV.SFTno;fileno++) {
       fprintf(stdout,"After sorting GV.filelist[%i] = %s\n",fileno,GV.filelist[fileno]);
       fflush(stdout);
     }
     fileno = GV.SFTno; /* restore value */
   #endif
  /***********************************************************/
  /*                                                         */
  /* END SECTION: sort the list of SFT files: GV.filelist.   */
  /*                                                         */
  /***********************************************************/

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

  return 0;
}
/******************************************/
/*                                        */
/* END FUNCTION: Freemem                  */
/*                                        */
/******************************************/
