/*
*  Copyright (C) 2007 Chris Messenger
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
 * \author Chris Messenger
 * \brief
 * This code has been written to quantify the sensitivity of a particular detector
 * given a specific sky position.  It is used to identify the most sensitive
 * stretch of data available from a large dataset.
 */

/************************************************************************************/
/* This code has been written to quantify the sensitivity of a particular detector  */
/* given a specific sky position.  It is used to identify the most sensitive        */
/* stretch of data available from a large dataset.                                  */
/*                                                                                  */
/*			           C. Messenger                                     */
/*                                                                                  */
/*                         BIRMINGHAM UNIVERISTY -  2004                            */
/************************************************************************************/


#define LAL_USE_OLD_COMPLEX_STRUCTS
#include "CalculateSensitivity_v1.h" 

int ReadSource(char *,char *,LIGOTimeGPS *,binarysource *); 
int ReadData(char *,binarysource *,dataset *,REAL8,INT4,FFT ***);
int CalculateSh(FFT ***,dataset *,REAL8 *);
int FreeSFTs(dataset *,FFT ***);   
int GetAB(EphemerisData *,char *,dataset *,binarysource *,sensresults *); 
int CalculateSensitivity(sensitivityparams *,binarysource *,sensresults **);
int OutputResult(char *,sensitivityparams *,sensresults *);
int SetupParams(sensitivityparams *,dataset *);
int SelectSFTs(FFT ***,FFT ***,dataset *,dataset *);
int ReadEphemeris(char *,char *,EphemerisData **);

char outdir[256],datadir[256],ephdir[256],yr[256],det[256],sourcename[256],sourcefile[256];
INT4 start;
INT4 end;
REAL8 tspan;
INT4 tstep;
INT4 maxsftno;
static LALStatus status;

int ReadCommandLine(int argc,char *argv[]);

int main(int argc, char **argv){
  
  FFT **SFTData;
  FFT **tempSFTData;
  dataset fulldataparams;
  sensresults *results=NULL;
  sensitivityparams sensparams;
  LIGOTimeGPS *obsstart=NULL;
  binarysource sourceparams;
  EphemerisData *edat=NULL;
  INT4 i,j;

   /* read the command line arguments */
  if (ReadCommandLine(argc,argv)) return 1;
  
  /* read the required source parameters from the input source file */
  if (ReadSource(sourcefile,sourcename,obsstart,&sourceparams)) return 1; 

  /* set up some structure inputs */
  fulldataparams.start.gpsSeconds=start;
  fulldataparams.start.gpsNanoSeconds=0;
  fulldataparams.end.gpsSeconds=end;
  fulldataparams.end.gpsNanoSeconds=0;
  sensparams.tspan=tspan;
  sensparams.overlap=tstep;
  
  /* start a loop over the number of bands */
  for (i=0;i<sourceparams.freq.nband;i++) {
    
    /*printf("*** Reading in full dataset\n");*/

    /* read the sub data into memory */
    if (ReadData(datadir,&sourceparams,&fulldataparams,sensparams.tspan,i,&SFTData)) return 2; 

    /*printf("*** Setting up parameters\n");*/

    /* set up sensitivity parameters ie number of chunks */
    if (SetupParams(&sensparams,&fulldataparams)) return 3;

    /* allocate results space */ /* sloppy */
    if (i==0) results=(sensresults *)LALMalloc(sensparams.nchunks*sizeof(sensresults));
      
    /*printf("*** Calculating noise floor sensitivity\n");*/

    /* initialise the noise floor results */
    for (j=0;j<sensparams.nchunks;j++) {
      results[j].ShAV[0]=0.0;
      results[j].ShAV[1]=0.0;
    }

    /* start a loop over the time chunks */
    for (j=0;j<sensparams.nchunks;j++) {

      /* if there are more than 1 stamps in this chunk */
      if (sensparams.dataparams[j].sftno>1) {

	/* copy required SFT's to a sub area */
	if (SelectSFTs(&SFTData,&tempSFTData,&sensparams.dataparams[j],&fulldataparams)) return 3;

	/* calculate the noise floor */
	if (CalculateSh(&tempSFTData,&sensparams.dataparams[j],&results[j].ShAV[i])) return 3; 

	/* free memory */
	if (FreeSFTs(&sensparams.dataparams[j],&tempSFTData)) return 7; 
      
      }

    }

    /*printf("*** Freeing the full dataset\n"); */

    /* free the full SFT data */
    if (FreeSFTs(&fulldataparams,&SFTData)) return 8;
        
    /* freeing the full data set parameters structure */
    LALFree(fulldataparams.stamps);
    
  }

  /*printf("*** Calculating sky postion sensitivity\n");*/

  /* read in epehemeris info */
  if (ReadEphemeris(ephdir,yr,&edat)) return 9;

  /* loop over the chunks again to calulate A and B */
  for (j=0;j<sensparams.nchunks;j++) {
      
    /* if we have more than 1 SFT in this chunk */
    if (sensparams.dataparams[j].sftno>1) { 

      /* calculate the values of A and B */
      if (GetAB(edat,det,&sensparams.dataparams[j],&sourceparams,&results[j])) return 4;      
    
    }
    
    /* freeing the timestamps in sensparams structure */
    LALFree(sensparams.dataparams[j].stamps);

  }

  /* freeing the ephemeris data */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat); 

  /*printf("*** Calculating overall sensitivity and outputting to file\n");*/

  /* evaluate the sensitivity of all the chunks for this source */
  if (CalculateSensitivity(&sensparams,&sourceparams,&results)) return 5;

  /* output the results */
  if (OutputResult(outdir,&sensparams,results)) return 6;
    
  exit(0);
  
}
/******************************************************************************/

int SelectSFTs(FFT ***SFTData,FFT ***tempSFTData,dataset *dataparams,dataset *fulldataparams)
{

  /* this function extracts SFTs from a larger set of SFT's */

  INT4 i,j,k,s;

  /* allocate memory */
  (*tempSFTData)=(FFT **)LALMalloc(dataparams->sftno*sizeof(FFT *));

    k=0;
  /* loop over all the SFT's */
  for (i=0;i<fulldataparams->sftno;i++) {
        
    /* loop over all the required stamps */
    for (j=0;j<dataparams->sftno;j++) {
  
    
      /* if the stamps match */
      if ((*SFTData)[i]->fft->epoch.gpsSeconds==dataparams->stamps[j].gpsSeconds) {
       
	/* allocate the memory for the temporary SFTs */
	(*tempSFTData)[k]=(FFT *)LALMalloc(sizeof(FFT));
	(*tempSFTData)[k]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries)); 
	(*tempSFTData)[k]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
	(*tempSFTData)[k]->fft->data->data=(COMPLEX8 *)LALMalloc((*SFTData)[i]->fft->data->length*sizeof(COMPLEX8));

	/* loop over the SFT data and copy it */
	for (s=0;s<(INT4)(*SFTData)[i]->fft->data->length;s++) {
	  (*tempSFTData)[k]->fft->data->data[s].realf_FIXME=crealf((*SFTData)[i]->fft->data->data[s]);
	  (*tempSFTData)[k]->fft->data->data[s].im=(*SFTData)[i]->fft->data->data[s].im;
	}

	/* increment the loop counter */
	k++;
	
      }
      
    }

  }
  
  return 0;
  
}


/******************************************************************************/

int SetupParams(sensitivityparams *sensparams,dataset *fulldataparams)
{

  /* this function sets up the required parameters for splitting the timestamps up */

  INT4 i,j,k;
  INT4 maxstamps;
  
  /*printf("tspan is %f\n",sensparams->tspan);
  printf("tsft is %d\n",fulldataparams->tsft);
  printf("fullstart id %d\n",fulldataparams->start.gpsSeconds);
  printf("fullend is %d\n",fulldataparams->end.gpsSeconds+(INT4)sensparams->tspan);*/

  maxstamps=(INT4)((REAL8)sensparams->tspan/(REAL8)fulldataparams->tsft);

  /*printf("maxstamps is %d\n",maxstamps);
    printf("overlap is %d\n",sensparams->overlap);*/

  /* calculate the number of chunks */
  
  sensparams->fullstart.gpsSeconds=fulldataparams->start.gpsSeconds;
  sensparams->fullstart.gpsNanoSeconds=fulldataparams->start.gpsNanoSeconds;
  sensparams->fullend.gpsSeconds=fulldataparams->end.gpsSeconds;
  sensparams->fullend.gpsNanoSeconds=fulldataparams->end.gpsNanoSeconds;
  sensparams->overlap=tstep;
  sensparams->nchunks=(INT4)((REAL8)(fulldataparams->end.gpsSeconds-fulldataparams->start.gpsSeconds)/(REAL8)sensparams->overlap);
  
  
  /* allocate memory for the chunks */
  sensparams->dataparams=(dataset *)LALMalloc(sensparams->nchunks*sizeof(dataset));

  /* loop over the chunks and define the boundaries of each chunk */
  for (i=0;i<sensparams->nchunks;i++) {
    
    sensparams->dataparams[i].start.gpsSeconds=fulldataparams->start.gpsSeconds+(INT4)(i*sensparams->overlap);
    sensparams->dataparams[i].end.gpsSeconds=sensparams->dataparams[i].start.gpsSeconds+sensparams->tspan;
    sensparams->dataparams[i].start.gpsNanoSeconds=0;
    sensparams->dataparams[i].end.gpsNanoSeconds=0;

    /* allocate memory for the stamps */
    sensparams->dataparams[i].stamps=(LIGOTimeGPS *)LALCalloc(maxstamps,sizeof(LIGOTimeGPS));

    k=0;
    /* loop over all the timestamps and copy the required SFTs */
    for (j=0;j<fulldataparams->sftno;j++) {
      
      /* if the stamp is in the range */
      if ((fulldataparams->stamps[j].gpsSeconds>=sensparams->dataparams[i].start.gpsSeconds)&&(fulldataparams->stamps[j].gpsSeconds<sensparams->dataparams[i].end.gpsSeconds)) {
    
	/* record the stamp */
	sensparams->dataparams[i].stamps[k].gpsSeconds=fulldataparams->stamps[j].gpsSeconds;
	sensparams->dataparams[i].stamps[k].gpsNanoSeconds=fulldataparams->stamps[j].gpsNanoSeconds;
	
	k++;
      }

      
      
    }
    
    /* fill in correct start and end info */
    sensparams->dataparams[i].sftno=k;
    if (k!=0) {
      sensparams->dataparams[i].start.gpsSeconds=sensparams->dataparams[i].stamps[0].gpsSeconds;
      sensparams->dataparams[i].start.gpsNanoSeconds=0;
      sensparams->dataparams[i].end.gpsSeconds=sensparams->dataparams[i].stamps[k-1].gpsSeconds+fulldataparams->tsft;
      sensparams->dataparams[i].end.gpsNanoSeconds=0;
      sensparams->dataparams[i].tsft=fulldataparams->tsft;
      sensparams->dataparams[i].nbins=fulldataparams->nbins;
      sensparams->dataparams[i].nsamples=fulldataparams->nsamples;
      sensparams->dataparams[i].tobs=fulldataparams->tsft*sensparams->dataparams[i].sftno;
    }
    else {
      sensparams->dataparams[i].tobs=0;
    }

  }


  return 0;

}
  
/******************************************************************************/
int OutputResult(char *outputdir,sensitivityparams *sensparams,sensresults *results) 
{

  /* this function simply outputs the result to file and standard out */

  FILE *fp;
  INT4 i;
  char outputfile[512];
  
  strcpy(outputfile,outputdir);
  sprintf(outputfile,"%s/sensitivity_%s_%s_%d-%d.data",outputdir,det,sourcename,start,end);
 
  fp=fopen(outputfile,"w");
  if (fp==NULL) {
    printf("ERROR : could not open output file %s\n",outputfile);
    exit(1);
  }

  fprintf(fp,"# Sensitivity results\n");
  fprintf(fp,"# -------------------------------------------------------------------------------- #\n");
  fprintf(fp,"# sourcefile = %s\n",sourcefile);
  fprintf(fp,"# source = %s\n",sourcename);
  fprintf(fp,"# full_start_time = %d\n",sensparams->fullstart.gpsSeconds);
  fprintf(fp,"# full_end_time = %d\n",sensparams->fullend.gpsSeconds);
  fprintf(fp,"# -------------------------------------------------------------------------------- #\n");
  fprintf(fp,"# tstart\ttend\t\tSh_1\t\tSh_2\t\tSh_wtd\t\tA\t\tB\t\tT_obs\t\tQ\n#\n");
      
  for (i=0;i<sensparams->nchunks;i++) {
    
    if (sensparams->dataparams[i].sftno>1) {
      
      fprintf(fp,"%d\t%d\t%e\t%e\t%e\t%f\t%f\t%f\t%e\n",
	      sensparams->dataparams[i].start.gpsSeconds,sensparams->dataparams[i].end.gpsSeconds,
	      results[i].ShAV[0],results[i].ShAV[1],results[i].ShAVweight,
	      results[i].A,results[i].B,
	      sensparams->dataparams[i].tobs,
	      results[i].Q);
    }
    
  }
 
  fclose(fp);
  
  return 0;
  
}

/******************************************************************************/

int CalculateSensitivity(sensitivityparams *sensparams,binarysource *sourceparams,sensresults **results)
{

  /* this function simply takes the calculated values of A, B, Sh, and T   */
  /* and combines them to give a quantity proportional to the detectable   */
  /* GW amplitude h0.  It weights each Sh by the size of its band.  It     */
  /* does this for all results */

  REAL8 band;
  REAL8 bandtot=0.0;
  INT4 i,j;

  /* printf("results A and B are %f %f\n",(*results[0]).A,(*results[0]).B);
     printf("Sh is %e %e\n",(*results[0]).ShAV[0],(*results[0]).ShAV[1]);*/

  
  /* loop over the number of chunks */
  for (j=0;j<sensparams->nchunks;j++) {

    if (sensparams->dataparams[j].tobs>sensparams->dataparams[j].tsft) {

      /* calculate the weighted ShAV */
      bandtot=0.0;
      (*results)[j].ShAVweight=0.0;
      for (i=0;i<sourceparams->freq.nband;i++) {
	band=sourceparams->freq.f_max[i]-sourceparams->freq.f_min[i];
	(*results)[j].ShAVweight+=band*(*results)[j].ShAV[i];
	bandtot+=band;
      }
      
      (*results)[j].ShAVweight/=bandtot;
      (*results)[j].Q=sqrt(5.0*(*results)[j].ShAVweight/(((*results)[j].A+(*results)[j].B)*sensparams->dataparams[j].tobs));

    }
    else {
      (*results)[j].ShAVweight=0.0;
      (*results)[j].Q=0.0;
    }

  
  }
 

  return 0;

}

/******************************************************************************/

int ReadData(char *datadirectory, binarysource *sourceparams, dataset *dataparams,REAL8 t_span,INT4 bandno,FFT ***SFTData) 
{

  INT4 filenum=0,offset;
  FILE *fp;
  size_t errorcode;
  UINT4 ndeltaf=0;
  char **filelist;
  char command[512];
  glob_t globbuf;
  UINT4 i;
  INT4 ifmin,ifmax;
  REAL8 f_min,f_max;
  INT4 fullstart;
  INT4 fullend;
  INT4 nfiles;
  INT4 flag=0;

  /* extract some variables from the sourceparams structure */
  f_min=sourceparams->freq.f_min[bandno];
  f_max=sourceparams->freq.f_max[bandno];

  /* set up the datadir name */
  strcpy(command, datadirectory);
  strcat(command,"/*");
    
  /* set up some glob stuff */
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);
  
  /* check if there are any SFT's in the directory */
  if(globbuf.gl_pathc==0)
    {
      fprintf (stderr,"\nNo SFTs in directory %s ... Exiting.\n", datadirectory);
      exit(1);
    }
  
  /* allocate memory for the pathnames */
  filelist=(char **)LALMalloc(globbuf.gl_pathc*sizeof(char *));
  for (i=0;i<globbuf.gl_pathc;i++) filelist[i]=(char *)LALMalloc(256*sizeof(char));
  
  /* read all file names into memory */
  while ((UINT4)filenum < globbuf.gl_pathc) 
    {
      strcpy(filelist[filenum],globbuf.gl_pathv[filenum]);
      filenum++;
      if (filenum > MAXFILES)
	{
	  fprintf(stderr,"\nToo many files in directory! Exiting... \n");
	  exit(1);
	}
    }
  globfree(&globbuf);

  nfiles=filenum;

  /* open the first file to read header information */
  if (!(fp=fopen(filelist[0],"rb"))) {
	fprintf(stderr,"Weird... %s doesn't exist!\n",filelist[0]);
	return 1;
      }
      
  /* Read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",filelist[0]);
      return 1;
    }
  
  /* Check that the time base is positive */
  if (header.tbase<=0.0)
    {
      fprintf(stderr,"Timebase %f from data file %s non-positive!\n",
	      header.tbase,filelist[filenum]);
      return 3;
    }

  /* read in start time */
  fullstart=header.gps_sec;

  /* close the test file */
  fclose(fp);

  /* define SFT time, nsamples, and max SFT number  */
  dataparams->tsft=header.tbase;
  dataparams->nsamples=header.nsamples;

  /* open up the last file in the list */
   if (!(fp=fopen(filelist[nfiles-1],"rb"))) {
	fprintf(stderr,"Weird... %s doesn't exist!\n",filelist[nfiles-1]);
	return 1;
      }
      
  /* Read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",filelist[filenum]);
      return 1;
    }

  /* read in end time */
  fullend=header.gps_sec+dataparams->tsft;


  fclose(fp);

  maxsftno=(INT4)((REAL8)(fullend-fullstart+t_span+(REAL8)dataparams->tsft)/(REAL8)dataparams->tsft);

  /* check that requested start time is within the dataset */
  if (((REAL8)(dataparams->start.gpsSeconds)+t_span)<fullstart) {
    fprintf(stderr,"ERROR : requested start time + t_span is before the dataset start\n");
    exit(1);
  }
  if ((REAL8)(dataparams->start.gpsSeconds)>fullend) {
    fprintf(stderr,"ERROR : requested start time is after the dataset end\n");
    exit(1);
  }

  /* allocate maximium memory for the SFT's and timestamps */ 
  (*SFTData)=(FFT **)LALMalloc(maxsftno*sizeof(FFT *));
  dataparams->stamps=(LIGOTimeGPS *)LALMalloc(maxsftno*sizeof(LIGOTimeGPS));

  /* loop over all files in the given data directory */
  dataparams->sftno=0;
  filenum=0;
  flag=0;
  while ((filenum<nfiles)&&(flag==0)) {
    
    /* open the SFT file */
    if (!(fp=fopen(filelist[filenum],"rb"))) {
      fprintf(stderr,"Weird... %s doesn't exist!\n",filelist[filenum]);
      return 1;
    }
      
    /* Read in the header from the file */
    errorcode=fread((void*)&header,sizeof(header),1,fp);
    if (errorcode!=1) 
      {
	fprintf(stderr,"No header in data file %s\n",filelist[filenum]);
	return 1;
      }
        
    /* if the data is within the requested observation window */
    if ((header.gps_sec>=dataparams->start.gpsSeconds)&&((header.gps_sec)<dataparams->end.gpsSeconds+t_span)) {

        /* Check that data is correct endian order */
      if (header.endian!=1.0)
	{
	  fprintf(stderr,"First object in file %s is not (double)1.0!\n",filelist[filenum]);
	  fprintf(stderr,"It could be a file format error (big/little\n");
	  fprintf(stderr,"endian) or the file might be corrupted\n\n");
	  return 2;
	}
      
      /* Check that the time base is positive */
      if (header.tbase<=0.0)
	{
	  fprintf(stderr,"Timebase %f from data file %s non-positive!\n",
		  header.tbase,filelist[filenum]);
	  return 3;
	}
      
      /* define indexes for SFT access */
      ifmax=ceil(f_max*dataparams->tsft);
      ifmin=floor(f_min*dataparams->tsft);
      
      /* Check that are frequency bins needed are in data set */
      if (ifmin<header.firstfreqindex || 
	  ifmax>header.firstfreqindex+header.nsamples) 
	{
	  fprintf(stderr,"Freq index range %d->%d not in %d to %d (file %s)\n",
		  ifmin,ifmax,header.firstfreqindex,
		  header.firstfreqindex+header.nsamples,filelist[filenum]);
	  return 4;
	}
      
      /*printf("found an SFT with time %d\n",header.gps_sec);*/
    
      /* Put time stamps from file into array */
      dataparams->stamps[dataparams->sftno].gpsSeconds = header.gps_sec;
      dataparams->stamps[dataparams->sftno].gpsNanoSeconds = header.gps_nsec;
      
      /* Move forward in file */
      offset=(ifmin-header.firstfreqindex)*2*sizeof(REAL4);
      errorcode=fseek(fp,offset,SEEK_CUR);
      if (errorcode) 
	{
	  perror(filelist[filenum]);
	  fprintf(stderr,"Can't get to offset %d in file %s\n",offset,filelist[filenum]);
	  return 5;
	}
      
      /* Make data structures */
      ndeltaf=ifmax-ifmin+1;
      (*SFTData)[dataparams->sftno]=(FFT *)LALMalloc(sizeof(FFT));
      (*SFTData)[dataparams->sftno]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries)); 
      (*SFTData)[dataparams->sftno]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
      (*SFTData)[dataparams->sftno]->fft->data->data=(COMPLEX8 *)LALMalloc(ndeltaf*sizeof(COMPLEX8));
	
      /* Fill in actual SFT data, and housekeeping */
      errorcode=fread((void*)((*SFTData)[dataparams->sftno]->fft->data->data), sizeof(COMPLEX8), ndeltaf, fp);
      if (errorcode!=ndeltaf){
	perror(filelist[filenum]);
	fprintf(stderr, "The SFT data was truncated.  Only read %zu not %d complex floats\n", errorcode, ndeltaf);
	return 6;
      }
	
      /* fill in extra SFT info */
      (*SFTData)[dataparams->sftno]->fft->epoch=dataparams->stamps[dataparams->sftno];
      (*SFTData)[dataparams->sftno]->fft->f0 = ifmin / dataparams->tsft;
      (*SFTData)[dataparams->sftno]->fft->deltaF = 1.0 / dataparams->tsft;
      (*SFTData)[dataparams->sftno]->fft->data->length = ndeltaf;
	
      /* increment true sft number */
      dataparams->sftno++;
      
      /* ending if data in observation window */ 
    }

    filenum++;

    /* close the file */
    fclose(fp);     
      
    }

  /* if we have found no SFTs ini the range */
  if (dataparams->sftno==0) {
    fprintf(stderr,"ERROR : No SFT's in directory %s between %d and %d !!! Exiting nicely\n",datadirectory,dataparams->start.gpsSeconds,dataparams->end.gpsSeconds+(INT4)t_span);
    exit(0);
  }

  
  /* define true start and end times and observation time */
  /*dataparams->start.gpsSeconds=dataparams->stamps[0].gpsSeconds;
  dataparams->start.gpsNanoSeconds=dataparams->stamps[0].gpsNanoSeconds;
  dataparams->end.gpsSeconds=dataparams->stamps[dataparams->sftno-1].gpsSeconds+dataparams->tsft;
  dataparams->end.gpsNanoSeconds=dataparams->stamps[dataparams->sftno-1].gpsNanoSeconds;
  dataparams->tobs=(REAL8)(dataparams->tsft*dataparams->sftno);*/
  dataparams->nbins=ndeltaf;

  for (i=0;i<(UINT4)nfiles;i++) LALFree(filelist[i]);
  LALFree(filelist);
   
  XLALRealloc((*SFTData),dataparams->sftno*sizeof(FFT *));
  XLALRealloc(dataparams->stamps,dataparams->sftno*sizeof(LIGOTimeGPS));


  return 0;  
}


/*********************************************************************************************/

int CalculateSh(FFT ***SFTData, dataset *dataparams,REAL8 *ShAV)
{

  INT4 k,j,i;                         /* loop indices */  
  REAL8 B;                        /* SFT Bandwidth */
  REAL8 deltaT,N;
  REAL8 *Sh;


  /*printf("tsft is %d\n",dataparams->tsft);
  printf("tobs is %fn",dataparams->tobs);
  printf("nbinsis %d\n",dataparams->nbins);
  printf("nsamples%d\n",dataparams->nsamples);
  printf("some data is now %e\n",(*SFTData[0])->fft->data->data[0].re);*/

   /* allocate memory */
  Sh=(REAL8 *)LALMalloc(dataparams->nbins*sizeof(REAL8));
  for (i=0;i<dataparams->nbins;i++) Sh[i]=0.0;

  /* loop over frequency bins */
  for (j=0;j<dataparams->nbins;j++) { 
    /* loop over each SFT */
    for (k=0;k<dataparams->sftno;k++) {
      Sh[j]=Sh[j]+
	(REAL8)crealf((*SFTData)[k]->fft->data->data[j]) * (REAL8)crealf((*SFTData)[k]->fft->data->data[j])+
	(REAL8)(*SFTData)[k]->fft->data->data[j].im * (REAL8)(*SFTData)[k]->fft->data->data[j].im;
    }
  }

  /* loop over the power and divide by Nsft to get the averege */
  for (j=0;j<dataparams->nbins;j++) { 
    Sh[j]=Sh[j]/(REAL8)dataparams->sftno;
  }

  N=2.0*(REAL8)dataparams->nsamples;
  B=(REAL8)dataparams->nsamples/dataparams->tsft;
  deltaT=1.0/(2.0*B);
  (*ShAV)=0.0;

   /* loop over the average power and correctly normalise as described by Xavier */
  for (j=0;j<dataparams->nbins;j++) { 
    Sh[j]=Sh[j]*2.0*deltaT/(REAL8)N;
    (*ShAV)+=Sh[j];
  }
  (*ShAV)/=(REAL8)dataparams->nbins;
    
  LALFree(Sh);

  return 0;
}

/*******************************************************************************/

int FreeSFTs(dataset *dataparams, FFT ***SFTData)
{

  INT4 k;
 
  /* here we free up the SFT memory */
  for (k=0;k<dataparams->sftno;k++)
    {
      /*printf("freeing #%d SFT\n",k);*/
      LALFree((*SFTData)[k]->fft->data->data);
      LALFree((*SFTData)[k]->fft->data);
      LALFree((*SFTData)[k]->fft);
      LALFree((*SFTData)[k]);
    }
  LALFree(*SFTData);

  return 0;
 
}


/*******************************************************************************/

int GetAB(EphemerisData *edat, char *detector, dataset *dataparams,binarysource *sourceparams,sensresults *results)
{

  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  LALDetector Detector;              /* Our detector*/
  EarthState earth;
  AMCoeffsParams *amParams;
  LIGOTimeGPS *midTS=NULL;           /* Time stamps for amplitude modulation coefficients */
  INT4 k;
  AMCoeffs amc;

  if(strcmp(detector,"GEO")) Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if(strcmp(detector,"LLO")) Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(strcmp(detector,"LHO")) Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
 
 
/* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.site.location[0]=Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=Detector.location[2]/LAL_C_SI;
  baryinput.alpha=sourceparams->skypos.ra;
  baryinput.delta=sourceparams->skypos.dec;
  baryinput.dInv=0.e0;

 
/* amParams structure to compute a(t) and b(t) */

/* Allocate space for amParams stucture */
/* Here, amParams->das is the Detector and Source info */
  amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
/* Fill up AMCoeffsParams structure */
  amParams->baryinput = &baryinput;
  amParams->earth = &earth; 
  amParams->edat = edat;
  amParams->das->pDetector = &Detector; 
  amParams->das->pSource->equatorialCoords.latitude = sourceparams->skypos.ra;
  amParams->das->pSource->equatorialCoords.longitude = sourceparams->skypos.dec;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!!!!!!!!*/

  /* Allocate space for AMCoeffs */
  amc.a = NULL;
  amc.b = NULL;
  LALSCreateVector(&status, &(amc.a), (UINT4) (dataparams->sftno-1));
  LALSCreateVector(&status, &(amc.b), (UINT4) (dataparams->sftno-1));

  /* Mid point of each SFT */
  midTS = (LIGOTimeGPS *)LALCalloc((dataparams->sftno-1),sizeof(LIGOTimeGPS));
  
  for(k=0; k<(dataparams->sftno-1); k++)
    { 
     

      midTS[k].gpsSeconds=(dataparams->stamps[k]).gpsSeconds+(int)(dataparams->tsft/2.0);
      midTS[k].gpsNanoSeconds=0;
    }
  
 
  LALComputeAM(&status, &amc, midTS, amParams); 

 
  results->A=amc.A;
  results->B=amc.B;

  LALFree(midTS);


  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);


  return 0;
}
/***********************************************************************************/

int ReadEphemeris(char *ephemdir,char *year,EphemerisData **edat) 
{
  
  char filenameE[256],filenameS[256];
  FILE *fp;
  
  strcpy(filenameE,ephemdir);
  strcat(filenameE,"/earth");
  strcat(filenameE,year);
  strcat(filenameE,".dat");
  
  strcpy(filenameS,ephemdir);
  strcat(filenameS,"/sun");
  strcat(filenameS,year);
  strcat(filenameS,".dat");

  /* *** Make sure the e-files are really there *** */
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
  /* ********************************************** */

  (*edat)=(EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*(*edat)).ephiles.earthEphemeris = filenameE;     
  (*(*edat)).ephiles.sunEphemeris = filenameS;    

  LALInitBarycenter(&status, *edat);               /* Reads in ephemeris files */

  return 0;

}

/***********************************************************************************/

 int ReadCommandLine(int argc,char *argv[]) 
{
  INT4 c, errflg = 0;
  CHAR *temp;
  optarg = NULL;
  
  /* Initialize default values */
  sprintf(datadir," ");
  sprintf(outdir," ");
  sprintf(sourcefile," ");
  sprintf(sourcename," ");
  start=0;
  end=0;
  tspan=0.0;
  tstep=0;
  sprintf(det,"LLO");
  sprintf(ephdir," ");
  sprintf(yr,"00-04");
  
  {
    int option_index = 0;
    static struct option long_options[] = {
      {"datadir", required_argument, 0, 'D'},
      {"outdir", required_argument, 0, 'O'},
      {"sourcefile", required_argument, 0, 'S'},
      {"source", required_argument, 0, 's'},
      {"start", required_argument, 0, 'a'},
      {"end", required_argument, 0, 'b'},
      {"tspan", required_argument, 0, 't'},
      {"tstep", required_argument, 0, 'T'},
      {"ephdir", required_argument, 0, 'E'},
      {"yr", required_argument, 0, 'y'},
      {"det", required_argument, 0, 'I'},
      {"help", no_argument, 0, 'h'}
    };
    /* Scan through list of command line arguments */
    while (!errflg && ((c = getopt_long (argc, argv,"hD:O:S:s:a:b:t:T:E:y:I:",long_options, &option_index)))!=-1)
      switch (c) {
      case 'D':
	temp=optarg;
	sprintf(datadir,"%s",temp);
	break;
      case 'O':
	temp=optarg;
	sprintf(outdir,"%s",temp);
	break;
      case 'S':
	temp=optarg;
	sprintf(sourcefile,"%s",temp);
	break;
      case 's':
	temp=optarg;
	sprintf(sourcename,"%s",temp);
	break;
      case 'a':
	start=atoi(optarg);
	break;
      case 'b':
	end=atoi(optarg);
	break;
      case 't':
	tspan=atof(optarg);
	break;
      case 'T':
	tstep=atoi(optarg);
	break;
      case 'E':
	temp=optarg;
	sprintf(ephdir,"%s",temp);
	break;
      case 'y':
	temp=optarg;
	sprintf(yr,"%s",temp);
	break;
      case 'I':
	temp=optarg;
	sprintf(det,"%s",temp);
	break;
      case 'h':
	/* print usage/help message */
	fprintf(stdout,"Arguments are:\n");
	fprintf(stdout,"\t--datadir     STRING\t Name of the directory where SFT's are stored [DEFAULT= ]\n");
	fprintf(stdout,"\t--sourcefile  STRING\t Location of the sourcefile containing the source parameters [DEFAULT= ]\n");
	fprintf(stdout,"\t--source      STRING\t Name of the source [DEFAULT= ]\n");
	fprintf(stdout,"\t--start       INT4\t GPS start time for full analysis [DEFAULT = datadir start]\n");
	fprintf(stdout,"\t--end         INT4\t GPS end time for full analysis [DEFAULT = datadir end]\n");
	fprintf(stdout,"\t--tspan       INT4\t Maximum observation span [DEFAULT = 0.0]\n");
	fprintf(stdout,"\t--tstep       INT4\t Step size [DEFAULT = 0.0]\n");
	fprintf(stdout,"\t--ephdir      STRING\t Location of ephemeris files earth?.dat and sun?.dat [DEFAULT=NULL]\n");
	fprintf(stdout,"\t--yr          STRING\t Year(s) specifying ephemeris files [DEFAULT=00-04]\n");
	fprintf(stdout,"\t--det         STRING\t Detector being used for the search (LLO,LHO,GEO,TAMA,CIT,VIRGO) [DEFAULT=LLO]\n");
	fprintf(stdout,"\t--outdir      STRING\t Name of output directory containing sensitivity info [DEFAULT=]\n");
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

