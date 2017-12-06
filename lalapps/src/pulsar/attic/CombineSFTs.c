/*
*  Copyright (C) 2007 Xavier Siemens
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
 * \author X. Siemens
 * \brief SFT combining code
 */

/*********************************************************************************/
/*                               SFT combining code                              */
/*                                                                               */
/*			             X. Siemens                                  */
/*                                                                               */
/*                                UWM - July  2003                               */
/*********************************************************************************/

#include "CombineSFTs.h"
char filelist[MAXFILES][MAXFILENAMELENGTH];
REAL8 sTsft,lTsft,B,OrigB;
int filenumber=0;
FFT **SFTData=NULL;                 /* SFT Data for LALDemod */
REAL8 *sinVal,*cosVal;        /*LUT values computed by the routine do_trig_lut*/
LIGOTimeGPS timestamps;       /* Time stamps from SFT data */

int ifmin,ifmax,if0,if1;

int main(int argc,char *argv[]) 
{
  int i=0,t1,t2,k;
  FILE *fsft1,*fsft2;
  size_t ecode;

  /* Reads command line arguments into the CommandLineArgs struct. 
     In the absence of command line arguments it sets some defaults */
  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;

  /* Sets global variables and sticks them into the GlobalVariables struct */  
  if (CreateFileList(CommandLineArgs)) return 2;

  if (AllocateMem(CommandLineArgs)) return 3;


  while(i<filenumber)
    {
      
      fsft1=fopen(filelist[i],"r");
      /* read in the header from the file */
      ecode=fread((void*)&header,sizeof(header),1,fsft1);
      if (ecode!=1) 
	{
	  fprintf(stderr,"No header in data file %s\n",filelist[i]);
	  return 1;
	}
      /* check that data is correct endian order */
      if (header.endian!=1.0)
	{
	  fprintf(stderr,"First object in file %s is not (double)1.0!\n",filelist[i]);
	  fprintf(stderr,"It could be a file format error (big/little\n");
	  fprintf(stderr,"endian) or the file might be corrupted\n\n");
	  return 2;
	}
      fclose(fsft1);
      t1=header.gps_sec;

      timestamps.gpsSeconds=header.gps_sec;
      timestamps.gpsNanoSeconds=header.gps_nsec;

      fsft2=fopen(filelist[i+CommandLineArgs.number-1],"r");
      /* read in the header from the file */
      ecode=fread((void*)&header,sizeof(header),1,fsft2);
      if (ecode!=1) 
	{
	  fprintf(stderr,"No header in data file %s\n",filelist[i]);
	  return 1;
	}
      /* check that data is correct endian order */
      if (header.endian!=1.0)
	{
	  fprintf(stderr,"First object in file %s is not (double)1.0!\n",filelist[i]);
	  fprintf(stderr,"It could be a file format error (big/little\n");
	  fprintf(stderr,"endian) or the file might be corrupted\n\n");
	  return 2;
	}
      fclose(fsft2);
      t2=header.gps_sec+sTsft;

/*       if( (int)((t2-t1)/lTsft +0.5) ==   1) */
/* 	{ */

	  fprintf(stdout,"Combining SFTs %d thru %d\n",i,i+CommandLineArgs.number-1);
	  /* Have CommandLineArgs.number consecutive SFTs */
	  if (ReadSFTs(CommandLineArgs,i)) return 3;
	  
	  if (CSFTs(CommandLineArgs)) return 4;
	  	  
	  i += CommandLineArgs.number;
/* 	} */
/*       else i++; */
    }


  /* Free the SFT data  */
  for (k=0;k<CommandLineArgs.number;k++)
    {
      LALFree(SFTData[k]->fft->data->data);
      LALFree(SFTData[k]->fft->data);
      LALFree(SFTData[k]->fft);
      LALFree(SFTData[k]);
    }
  LALFree(SFTData);
  
  LALFree(sinVal);
  LALFree(cosVal);

  LALCheckMemoryLeaks();

  return 0;

}
/*******************************************************************************/
int AllocateMem(struct CommandLineArgsTag CLA)
{ 
  int ndeltaf,k;
  INT4  res;                    /*resolution of the argument of the trig functions: 2pi/res.*/

  SFTData=(FFT **)LALMalloc(CLA.number*sizeof(FFT *));

  ndeltaf=ifmax-ifmin+1;

  for (k=0;k<CLA.number;k++)
    {
      SFTData[k]=(FFT *)LALMalloc(sizeof(FFT));
      SFTData[k]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
      SFTData[k]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
      SFTData[k]->fft->data->data=(COMPLEX8 *)LALMalloc(ndeltaf*sizeof(COMPLEX8));
    }

  /* res=10*(params->mCohSFT); */
  /* This size LUT gives errors ~ 10^-7 with a three-term Taylor series */
   res=64;
   sinVal=(REAL8 *)LALMalloc((res+1)*sizeof(REAL8));
   cosVal=(REAL8 *)LALMalloc((res+1)*sizeof(REAL8)); 
   for (k=0; k<=res; k++){
     sinVal[k]=sin((LAL_TWOPI*k)/res);
     cosVal[k]=cos((LAL_TWOPI*k)/res);
   }

  return 0;

}
/*******************************************************************************/
int CSFTs(struct CommandLineArgsTag CLA)
{ 
  INT4 alpha,m;                 /* loop indices */
  REAL8	xTemp;	                /* temp variable for phase model */
  REAL8	deltaF;	                /* width of SFT band */
  INT4	k, k1;	                /* defining the sum over which is calculated */
  REAL8	x;		        /* local variable for holding x */
  REAL8	realXP, imagXP; 	/* temp variables used in computation of */
  REAL8	realP, imagP;	        /* real and imaginary parts of P, see CVS */
  INT4	nDeltaF;	        /* number of frequency bins per SFT band */
  INT4	sftIndex;	        /* more temp variables */
  REAL8	y;		        /* local variable for holding y */
  REAL8 realQ, imagQ;

  INT4 index;

  COMPLEX16 llSFT;
  COMPLEX8 lSFT;

  REAL8 f;
  FILE *fp;
  char filename[256], fnuber[16];
  int errorcode;

  /* variable redefinitions for code readability */
  deltaF=(*SFTData)->fft->deltaF;
  nDeltaF=(*SFTData)->fft->data->length;

  /* Open SFT data file */
  strcpy(filename,CLA.outputdirectory);
  sprintf(fnuber,"CSFT.%09d",timestamps.gpsSeconds); /*ACHTUNG: used to be iSFT+1 - now starts from .00000*/
  strcat(filename,fnuber);
  fp=fopen(filename,"w");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",filename);
    return 1;
  }

  header.endian=1.0;
  header.gps_sec=timestamps.gpsSeconds;
  header.gps_nsec=timestamps.gpsNanoSeconds;
  header.tbase=lTsft;
  header.firstfreqindex=(INT4)(if0*lTsft/sTsft+0.5);
  header.nsamples=CLA.number*(if1-if0);
  
  /* write header */
  errorcode=fwrite((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1){
    printf("Error in writing header into file!\n");
    return 1;
  }


  // Loop over frequencies to be demodulated
  for(m = 0 ; m <= CLA.number*(if1-if0)  ; m++ )
  {
    llSFT.re =0.0;
    llSFT.im =0.0;

    f=if0*deltaF+m*deltaF/CLA.number;

    // Loop over SFTs that contribute to F-stat for a given frequency
    for(alpha=0;alpha<CLA.number;alpha++)
      {
	REAL8 tsin, tcos, tempFreq;
	COMPLEX8 *Xalpha=SFTData[alpha]->fft->data->data;
	xTemp=(REAL8)if0+(REAL8)m/(REAL8)CLA.number;
	realXP=0.0;
	imagXP=0.0;
	      	/* find correct index into LUT -- pick closest point */
	tempFreq=xTemp-(INT4)xTemp;
	index=(INT4)(tempFreq*64+0.5); //just like res above
	      
	{
	  REAL8 d=LAL_TWOPI*(tempFreq-(REAL8)index/64.0);//just like res above
	  REAL8 d2=0.5*d*d;
	  REAL8 ts=sinVal[index];
	  REAL8 tc=cosVal[index];
		
	  tsin=ts+d*tc-d2*ts;
	  tcos=tc-d*ts-d2*tc-1.0;
	}

        tempFreq=LAL_TWOPI*(tempFreq+CLA.Dterms-1);
        k1=(INT4)xTemp-CLA.Dterms+1;
        /* Loop over terms in dirichlet Kernel */
        for(k=0;k<2*CLA.Dterms;k++)
	  {
	    COMPLEX8 Xalpha_k;
	    x=tempFreq-LAL_TWOPI*(REAL8)k;
	    realP=tsin/x;
	    imagP=tcos/x;

	    /* If x is small we need correct x->0 limit of Dirichlet kernel */
	    if(fabs(x) < 0.000001) 
	      {
		realP=1.0;
		imagP=0.0;
	      }	 
 
	    sftIndex=k1+k-ifmin;

	    /* these four lines compute P*xtilde */
	    Xalpha_k=Xalpha[sftIndex];
	    realXP += Xalpha_k.re*realP;
	    realXP -= Xalpha_k.im*imagP;
	    imagXP += Xalpha_k.re*imagP;
	    imagXP += Xalpha_k.im*realP;
	  }
      
	y=-LAL_TWOPI*(alpha)*(if0+(REAL8)m/(REAL8)CLA.number);

	      //    f=if0*deltaF+m*deltaF/CLA.number;

/*              y=-LAL_TWOPI*(alpha+1)*((REAL8)beta*oneOverMCohSFT+(REAL8)l); */

	realQ = cos(y);
	imagQ = sin(y);

	/* implementation of amplitude demodulation */
	{
	  REAL8 realQXP = realXP*realQ-imagXP*imagQ;
	  REAL8 imagQXP = realXP*imagQ+imagXP*realQ;
	  llSFT.re += realQXP;
	  llSFT.im += imagQXP;
	}
      }      

    /* Adjust normalisation to compensate for change in B; this way we can compute things just with header */

     lSFT.re = B/OrigB*llSFT.re; 
     lSFT.im = B/OrigB*llSFT.im; 
    
    errorcode=fwrite((void*)&lSFT.re, sizeof(REAL4),1,fp);  
    errorcode=fwrite((void*)&lSFT.im, sizeof(REAL4),1,fp);  

  }

  fclose(fp);

  return 0;

}



/*******************************************************************************/



int ReadSFTs(struct CommandLineArgsTag CLA, int i)
{
  INT4 fn=0,offset;
  FILE *fp;
  size_t errorcode;

  for (fn=0;fn<CLA.number;fn++)
    {
      /* open FIRST file and get info from it*/

      fp=fopen(filelist[fn+i],"r");
      if (fp==NULL) 
        {
          fprintf(stderr,"Weird... %s doesn't exist!\n",filelist[fn+i]);
          return 1;
        }
      /* Read in the header from the file */
      errorcode=fread((void*)&header,sizeof(header),1,fp);
      if (errorcode!=1) 
        {
          fprintf(stderr,"No header in data file %s\n",filelist[fn+i]);
          return 1;
        }

      /* Check that data is correct endian order */
      if (header.endian!=1.0)
        {
          fprintf(stderr,"First object in file %s is not (double)1.0!\n",filelist[fn+i]);
          fprintf(stderr,"It could be a file format error (big/little\n");
          fprintf(stderr,"endian) or the file might be corrupted\n\n");
          return 2;
        }
      /* Check that the time base is positive */
      if (header.tbase<=0.0)
        {
          fprintf(stderr,"Timebase %f from data file %s non-positive!\n",
                  header.tbase,filelist[fn+i]);
          return 3;
        }
        

      /* Check that are frequency bins needed are in data set */
      if ( ifmin < header.firstfreqindex || 
           ifmax > header.firstfreqindex+header.nsamples) 
        {
          fprintf(stderr,"Freq index range %d -> %d  not in %d to %d (file %s)\n",
                  ifmin,ifmax, header.firstfreqindex,
                  header.firstfreqindex+header.nsamples,filelist[fn+i]);
          return 4;
        }

      /* Move forward in file */
      offset=(ifmin-header.firstfreqindex)*2*sizeof(REAL4);
      errorcode=fseek(fp,offset,SEEK_CUR);
      if (errorcode) 
        {
          perror(filelist[fn+i]);
          fprintf(stderr,"Can't get to offset %d in file %s\n",offset,filelist[fn+i]);
          return 5;
        }

    /* Make data structures */

    /* Fill in actual SFT data, and housekeeping */
    errorcode=fread((void*)(SFTData[fn]->fft->data->data),(ifmax-ifmin+1)*sizeof(COMPLEX8),1,fp);
    SFTData[fn]->fft->f0=ifmin/sTsft;
    SFTData[fn]->fft->deltaF=1.0/sTsft;
    SFTData[fn]->fft->data->length=ifmax-ifmin+1;

    fclose(fp);     /* close file */
    }



  return 0;  
}
/*******************************************************************************/

int CreateFileList(struct CommandLineArgsTag CLA)
{

  char command[256];
  glob_t globbuf;
  FILE *fp;
  size_t errorcode;
  float minimumF,maximumF;

  strcpy(command,CLA.inputdirectory);
  strcat(command,"/*SFT*");
  
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);

  /* read in file names */

  if(globbuf.gl_pathc==0)
    {
      fprintf(stderr,"No SFTs in directory %s ... Exiting.\n",CLA.inputdirectory);
      return 1;
    }

  while (filenumber< globbuf.gl_pathc) 
    {
      strcpy(filelist[filenumber],globbuf.gl_pathv[filenumber]);
      filenumber++;
      if (filenumber > MAXFILES)
        {
          fprintf(stderr,"Too many files in directory! Exiting... \n");
          return 1;
        }
    }
  fprintf(stdout,"%d files in directory\n",filenumber);
  globfree(&globbuf);

  /* open FIRST file and get info from it*/
  fp=fopen(filelist[0],"r");
  /* read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",filelist[0]);
      return 1;
    }

  /* check that data is correct endian order */
  if (header.endian!=1.0)
    {
      fprintf(stderr,"First object in file %s is not (double)1.0!\n",filelist[0]);
      fprintf(stderr,"It could be a file format error (big/little\n");
      fprintf(stderr,"endian) or the file might be corrupted\n\n");
      return 2;
    }
  fclose(fp);

  sTsft=header.tbase;  // Time baseline of short SFTs
  lTsft=CLA.number*header.tbase;  // Time baseline of long SFTs

  OrigB=((float)header.nsamples) / sTsft; // Original bandwidth of SFTs

  minimumF=150.0; //150.0
  maximumF=160.0; //1300.0  is what we agreed with map

  ifmin=floor(minimumF*sTsft) - CLA.Dterms;
  ifmax=ceil(maximumF*sTsft) + CLA.Dterms;

  if0=floor(minimumF*sTsft);
  if1=ceil(maximumF*sTsft);

  B=((float) CLA.number*(if1-if0))/(CLA.number*sTsft);

  return 0;  
}

/*******************************************************************************/

   int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{
  INT4 c, errflg = 0;
  optarg = NULL;
  
  /* Initialize default values */
  CLA->Dterms=16;
  CLA->inputdirectory="";
  CLA->outputdirectory="";
  CLA->number=30;

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"ht:d:i:o:n:"))!=-1))
    switch (c) {
    case 't':
      /* frequency bandwidth */
      CLA->Dterms=atof(optarg);
      break;
    case 'i':
      /* starting observation time */
      CLA->inputdirectory=optarg;
      break;
    case 'o':
      /* starting observation time */
      CLA->outputdirectory=optarg;
      break;
    case 'n':
      /* starting observation time */
      CLA->number=atof(optarg);
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\t-t\tINTEGER\t Number of terms to keep in Dirichlet kernel sum (default 16)\n");
      fprintf(stdout,"\t-i\tSTRING\t Directory where input SFT's are located (not set by default) \n");
      fprintf(stdout,"\t-o\tSTRING\t Directory where input SFT's will be located (not set by default) \n");
      fprintf(stdout,"\t-n\tINTEGER\t Number of input SFTs that make one output SFT (30 by default) \n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
  if(CLA->inputdirectory == "")
    {
      fprintf(stderr,"No SFT directory specified; input directory with -D option.\n");
      fprintf(stderr,"Try ./CombineSFTs -h \n");
      return 1;
    }      
  if(CLA->outputdirectory == "")
    {
      fprintf(stderr,"No output directory specified; input directory with -D option.\n");
      fprintf(stderr,"Try ./CombineSFTs -h \n");
      return 1;
    }      

  /* update global variable and return */
  return errflg;
}


/*******************************************************************************/
