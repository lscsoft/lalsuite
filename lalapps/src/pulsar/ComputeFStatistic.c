/*********************************************************************************/
/*                    F-statistic generation code for known pulsars              */
/*                                                                               */
/*			      M.A. Papa and X. Siemens                           */
/*                                                                               */
/*                 Albert Einstein Institute/UWM - September 2002                */
/*********************************************************************************/

#include <errno.h>
#include "ComputeFStatistic.h"

FFT **SFTData=NULL;                 /* SFT Data for LALDemod */
DemodPar *DemodParams  = NULL;      /* Demodulation parameters for LALDemod */
LIGOTimeGPS *timestamps=NULL;       /* Time stamps from SFT data */
double *F;
INT4 lalDebugLevel=0,i;
static LALStatus status;
AMCoeffs amc;
GlobalVariables GV;
REAL4 MeanOneOverSh=0.0;

int main(int argc,char *argv[]) 
{

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;

  if (SetGlobalVariables(CommandLineArgs)) return 2;

  if (ReadSFTData()) return 3;

  if (NormaliseSFTData()) return 4;

  if (CreateDemodParams(CommandLineArgs)) return 5;

  LALDemod(&status,F,SFTData,DemodParams);

  for(i=0;i < GV.imax ;i++)
    {
      fprintf(stdout,"%20.10f %20.15f\n",GV.startingdemodfreq+i*GV.demodfreqres,F[i]);
    }

  if (Freemem()) return 8;

  return 0;

}

/*******************************************************************************/

int CreateDemodParams(struct CommandLineArgsTag CLA)
{
  CSParams *csParams  = NULL;        /* ComputeSky parameters */
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  EphemerisData *edat=NULL;          /* Stores earth/sun ephemeris data for barycentering */
  LALDetector Detector;              /* Our detector*/
  EarthState earth;
  EmissionTime emit;
  AMCoeffsParams *amParams;
  LIGOTimeGPS *midTS=NULL;           /* Time stamps for amplitude modulation coefficients */
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4 leap;


  char filenameE[256],filenameS[256];

  INT4 i,k;
  FILE *fp;
  
  strcpy(filenameE,CLA.efiles);
  strcat(filenameE,"/earth02.dat");

  strcpy(filenameS,CLA.efiles);
  strcat(filenameS,"/sun02.dat");

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

  edat=(EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = filenameE;     
  (*edat).ephiles.sunEphemeris = filenameS;         

  LALLeapSecs(&status,&leap,&timestamps[0],&formatAndAcc);
  (*edat).leap=leap;


  LALInitBarycenter(&status, edat);               /* Reads in ephemeris files */

  if(CLA.IFO == 0) Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if(CLA.IFO == 1) Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(CLA.IFO == 2) Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
  if(CLA.IFO == 3)
    {
        if (CreateDetector(&Detector)) return 5;
    }
 
/* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.site.location[0]=Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=Detector.location[2]/LAL_C_SI;
  baryinput.alpha=CLA.skyalpha;
  baryinput.delta=CLA.skydelta;
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
  amParams->das->pSource->equatorialCoords.latitude = CLA.skydelta;
  amParams->das->pSource->equatorialCoords.longitude = CLA.skyalpha;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!!!!!!!!*/
  amParams->leapAcc=formatAndAcc.accuracy;


/* Allocate space for AMCoeffs */
  amc.a = NULL;
  amc.b = NULL;
  LALSCreateVector(&status, &(amc.a), (UINT4) GV.SFTno);
  LALSCreateVector(&status, &(amc.b), (UINT4) GV.SFTno);

 /* Mid point of each SFT */
   midTS = (LIGOTimeGPS *)LALCalloc(GV.SFTno,sizeof(LIGOTimeGPS));
   for(k=0; k<GV.SFTno; k++)
     { 
       REAL8 teemp=0.0;
       LALGPStoFloat(&status,&teemp, &(timestamps[k]));
       teemp += 0.5*GV.tsft;
       LALFloatToGPS(&status,&(midTS[k]), &teemp);
     }
   
   LALComputeAM(&status, &amc, midTS, amParams); 

   /*   LALComputeAM(&status, &amc, midTS, amParams);  */

/* ComputeSky stuff*/
  csParams=(CSParams *)LALMalloc(sizeof(CSParams));
  csParams->tGPS=timestamps;  
  csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
  csParams->mObsSFT=GV.SFTno;     /* Changed this from GV.mobssft !!!!!! */
  csParams->tSFT=GV.tsft;
  csParams->edat=edat;
  csParams->baryinput=&baryinput;
  csParams->spinDwnOrder=CLA.spinDwnOrder;
  csParams->skyPos[0]=CLA.skyalpha;
  csParams->skyPos[1]=CLA.skydelta;
  csParams->earth = &earth;
  csParams->emit = &emit;

/* Finally, DemodParams */

  /* Allocate DemodParams structure */
  DemodParams=(DemodPar *)LALMalloc(sizeof(DemodPar));
  
  /* space for sky constants */
/* Based on maximum index for array of as and bs sky constants as from ComputeSky.c */
  k=2*(GV.SFTno-1)*(CLA.spinDwnOrder+1)+2*CLA.spinDwnOrder+2; 
  DemodParams->skyConst=(REAL8 *)LALMalloc(k*sizeof(REAL8));
  DemodParams->amcoe=&amc;
  DemodParams->spinDwnOrder=CLA.spinDwnOrder;
  DemodParams->SFTno=GV.SFTno;

  DemodParams->f0=GV.startingdemodfreq;
  DemodParams->imax=GV.imax;
  DemodParams->df=GV.demodfreqres;

  DemodParams->Dterms=GV.Dterms;
  DemodParams->ifmin=GV.ifmin;

  ComputeSky(&status,DemodParams->skyConst,0,csParams);        /* compute the */
						               /* "sky-constants" A and B */
  /* space for spin down params */
  DemodParams->spinDwn=(REAL8 *)LALMalloc(CLA.spinDwnOrder*sizeof(REAL8));
  /* Spindown parameters for  */
  for (i=0;i<CLA.spinDwnOrder;i++)
    {
      DemodParams->spinDwn[i]=CLA.spinParams[i]; 
    }  

  LALFree(midTS);

  LALFree(csParams->skyPos);
  LALFree(csParams);

  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);

  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  return 0;
}

/*******************************************************************************/

int NormaliseSFTData()
{
  INT4 i,j;                         /* loop indices */
  INT4 nbins=GV.ifmax-GV.ifmin+1;   /* Number of points in SFT's */
  REAL4 SFTsqav;                  /* Average of Square of SFT */
  REAL8 B;                        /* SFT Bandwidth */
  REAL8 deltaT,N;


  /* loop over each SFTs */
  for (i=0;i<GV.SFTno;i++)         
    {
      
      SFTsqav=0.0;
      /* loop over SFT data to estimate noise */
      for (j=0;j<nbins;j++)               
	  {
	    SFTsqav=SFTsqav+
	      SFTData[i]->fft->data->data[j].re * SFTData[i]->fft->data->data[j].re+
	      SFTData[i]->fft->data->data[j].im * SFTData[i]->fft->data->data[j].im;
	  }
      SFTsqav=SFTsqav/(1.0*nbins);              /* Actual average of Square of SFT */
      MeanOneOverSh=2.0*GV.nsamples*GV.nsamples/(SFTsqav*GV.tsft)+MeanOneOverSh;      

      N=1.0/sqrt(2.0*SFTsqav);

/* signal only case */  
      if(GV.noise == 1)
	{
	  B=(1.0*GV.nsamples)/(1.0*GV.tsft);
	  deltaT=1.0/(2.0*B);
	  N=deltaT/sqrt(GV.tsft);
	}

  /* loop over SFT data to Normalise it*/
      for (j=0;j<nbins;j++)               
	{
	  SFTData[i]->fft->data->data[j].re = N*SFTData[i]->fft->data->data[j].re; 
	  SFTData[i]->fft->data->data[j].im = N*SFTData[i]->fft->data->data[j].im;
	}
    } 
	
   MeanOneOverSh=MeanOneOverSh/(1.0*GV.SFTno); 
  return 0;
}

/*******************************************************************************/

int ReadSFTData(void)
{
  INT4 fileno=0,ndeltaf,offset;
  FILE *fp;
  size_t errorcode;

  SFTData=(FFT **)LALMalloc(GV.SFTno*sizeof(FFT *));
  timestamps=(LIGOTimeGPS *)LALMalloc(GV.SFTno*sizeof(LIGOTimeGPS));

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
      /* Put time stamps from file into array */
      timestamps[fileno].gpsSeconds=header.gps_sec;
      timestamps[fileno].gpsNanoSeconds=header.gps_nsec;

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
    ndeltaf=GV.ifmax-GV.ifmin+1;
    SFTData[fileno]=(FFT *)LALMalloc(sizeof(FFT));
    SFTData[fileno]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
    SFTData[fileno]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
    SFTData[fileno]->fft->data->data=(COMPLEX8 *)LALMalloc(ndeltaf*sizeof(COMPLEX8));

    /* Fill in actual SFT data, and housekeeping */
    errorcode=fread((void*)(SFTData[fileno]->fft->data->data),ndeltaf*sizeof(COMPLEX8),1,fp);
    SFTData[fileno]->fft->epoch=timestamps[fileno];
    SFTData[fileno]->fft->f0=GV.ifmin*GV.df;
    SFTData[fileno]->fft->deltaF=GV.df;
    SFTData[fileno]->fft->data->length=ndeltaf;

    fclose(fp);     /* close file */

    }
  return 0;  
}

/*******************************************************************************/

int SetGlobalVariables(struct CommandLineArgsTag CLA)
{

  char command[256];
  FILE *fp;
  size_t errorcode;
  REAL8 df;                         /* freq resolution */
  INT4 fileno=0;   
  glob_t globbuf;

  strcpy(command,CLA.directory);
  strcat(command,"/*SFT*");
  
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);

  /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */

  if(globbuf.gl_pathc==0)
    {
      fprintf(stderr,"No SFTs in directory %s ... Exiting.\n",CLA.directory);
      return 1;
    }

  while ((UINT4)fileno < globbuf.gl_pathc) 
    {
      strcpy(GV.filelist[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
      if (fileno > MAXFILES)
	{
	  fprintf(stderr,"Too many files in directory! Exiting... \n");
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
    
  /* variables for starting demodulation frequency, band and resolution */
  GV.startingdemodfreq=CLA.startingdemodfreq;
  GV.demodband=CLA.demodband;
  GV.demodfreqres=CLA.demodfreqres;


  /* if user has not input demodulation frequency resolution; set to 1/Tobs */
  if(CLA.demodfreqres == 0.0 ) GV.demodfreqres=1.0/(header.tbase*GV.SFTno);

  GV.imax=(int)(GV.demodband/GV.demodfreqres+.5)+1;  /*Number of frequency values to calculate F for */

  GV.nsamples=header.nsamples;    /* # of freq. bins */

  /* frequency resolution: used only for printing! */
  df=(1.0)/(1.0*header.tbase);
  GV.df=df;

  GV.ifmax=ceil((1.0+DOPPLERMAX)*(GV.startingdemodfreq+GV.demodband)*GV.tsft)+CLA.Dterms;
  GV.ifmin=floor((1.0-DOPPLERMAX)*GV.startingdemodfreq*GV.tsft)-CLA.Dterms;

  GV.Dterms=CLA.Dterms;

  GV.noise=CLA.noise;

  /* Tell the user what we have arrived at */
/*    fprintf(stdout,"\n");    */
/*    fprintf(stdout,"SFT time baseline:                  %f min\n",header.tbase/60.0);  */
/*    fprintf(stdout,"SFT freq resolution:                %f Hz\n",df);  */
/*    fprintf(stdout,"Starting search frequency:          %f Hz\n",GV.startingdemodfreq);  */
/*    fprintf(stdout,"Demodulation frequency band:        %f Hz\n",GV.demodband);  */
/*    fprintf(stdout,"# of SFT in a DeFT:                 %f\n",ceil((1.0*(GV.Tf - GV.Ti))/header.tbase));  */
/*    fprintf(stdout,"Actual # of SFTs:                   %d\n",GV.SFTno);  */
/*    fprintf(stdout,"==> DeFT baseline:                  %f hours\n",(GV.Tf - GV.Ti)/3600.0);  */
/*    fprintf(stdout,"# of points taken from each SFT:    %d %d\n",GV.ifmax,GV.ifmin);  */
/*    fprintf(stdout,"# of points in SFT:    %d\n",GV.nsamples);  */


   /* allocate F-statistic array */

  F=(double *)LALMalloc(GV.imax*sizeof(double));


  return 0;  
}

/*******************************************************************************/

   int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{
  INT4 c, errflg = 0;
  optarg = NULL;
  
  /* Initialize default values */
  CLA->Dterms=16;
  CLA->startingdemodfreq=0.0;
  CLA->skyalpha=0.0;
  CLA->skydelta=0.0;
  CLA->spinDwnOrder=1;               
  CLA->spinParams[0]=0.0;            
  CLA->spinParams[1]=0.0;            
  CLA->spinParams[2]=0.0;            
  CLA->spinParams[3]=0.0;            
  CLA->directory="";
  CLA->efiles="";
  CLA->IFO=-1;
  CLA->noise=0;
  CLA->demodband=0.0;
  CLA->demodfreqres=0.0;

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"ht:f:a:d:D:E:s:1:2:3:4:I:n:b:r:"))!=-1))
    switch (c) {
    case 't':
      /* frequency bandwidth */
      CLA->Dterms=atof(optarg);
      break;
    case 'I':
      /* frequency bandwidth */
      CLA->IFO=atof(optarg);
      break;
    case 'n':
      /* frequency bandwidth */
      CLA->noise=atof(optarg);
      break;
    case 'f':
      /* first search frequency */
      CLA->startingdemodfreq=atof(optarg);
      break;
    case 'b':
      /* first search frequency */
      CLA->demodband=atof(optarg);
      break;
    case 'r':
      /* first search frequency */
      CLA->demodfreqres=atof(optarg);
      break;
    case 'D':
      /* starting observation time */
      CLA->directory=optarg;
      break;
    case 'E':
      /* starting observation time */
      CLA->efiles=optarg;
      break;
    case 'a':
      /* sky position alpha */
      CLA->skyalpha=atof(optarg);
      break;
    case 'd':
      /* sky position delta */
      CLA->skydelta=atof(optarg);
      break;
    case 's':
      /* Spin down order */
      CLA->spinDwnOrder=atof(optarg);
      break;
      if(CLA->spinDwnOrder > 4)  /* THIS DOESN'T WORK, WEIRD!!!!!*/
      {
	fprintf(stderr,"Maximum Spindown order is 4.");
	errflg++;
	exit(1);
	break;
      }
    case '1':
      /* Spin down order */
      CLA->spinParams[0]=atof(optarg);
      break;
    case '2':
      /* Spin down order */
      CLA->spinParams[1]=atof(optarg);
      break;
    case '3':
      /* Spin down order */
      CLA->spinParams[2]=atof(optarg);
      break;
    case '4':
      /* Spin down order */
      CLA->spinParams[3]=atof(optarg);
      break;

    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\t-f\tFLOAT\t Starting search frequency in Hz (not set by default)\n");
      fprintf(stdout,"\t-b\tFLOAT\t Demodulation frequency band in Hz (set=0.0 by default)\n");
      fprintf(stdout,"\t-r\tFLOAT\t Demodulation frequency resolution in Hz (set to 1/Tobs by default)\n");
      fprintf(stdout,"\t-t\tINTEGER\t Number of terms to keep in Dirichlet kernel sum (default 16)\n");
      fprintf(stdout,"\t-a\tFLOAT\t Sky position alpha (equatorial coordinates) in radians (default 0.0)\n");
      fprintf(stdout,"\t-d\tFLOAT\t Sky position delta (equatorial coordinates) in radians (default 0.0)\n");
      fprintf(stdout,"\t-D\tSTRING\t Directory where SFT's are located (not set by default) \n");
      fprintf(stdout,"\t-E\tSTRING\t Directory where Ephemeris files are located (not set by default) \n");
      fprintf(stdout,"\t-I\tSTRING\t Detector; must be set to 0=GEO, 1=LLO, 2=LHO or 3=Roman Bar (not set by default) \n");
      fprintf(stdout,"\t-n\tSTRING\t Noise only; must be set to 0=Regular data, 1=Signal only data (default 0) \n");
      fprintf(stdout,"\t-s\tINTEGER\t Spindown order (default 1) \n");
      fprintf(stdout,"\t-1\tFLOAT\t 1st spindown parameter (default 0.0)\n");
      fprintf(stdout,"\t-2\tFLOAT\t 2nd spindown parameter (default 0.0)\n");
      fprintf(stdout,"\t-3\tFLOAT\t 3rd spindown parameter (default 0.0)\n");
      fprintf(stdout,"\t-4\tFLOAT\t 4th spindown parameter (default 0.0)\n");

      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
  if(CLA->IFO == -1)
    {
      fprintf(stderr,"No IFO specified; input with -I option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      return 1;
    }      
  if(CLA->directory == "")
    {
      fprintf(stderr,"No SFT directory specified; input directory with -D option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      return 1;
    }      
  if(CLA->efiles == "")
    {
      fprintf(stderr,"No ephemeris data (earth??.dat, sun??.dat) directory specified; input directory with -E option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      return 1;
    }      
  if(CLA->startingdemodfreq == 0.0)
    {
      fprintf(stderr,"No search frequency specified; set with -f option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
     return 1;
    }      

  /* update global variable and return */
  return errflg;
}


/*******************************************************************************/

int CreateDetector(LALDetector *Detector){

/*   LALDetector Detector;  */
  LALFrDetector detector_params;
  LALDetectorType bar;
  LALDetector Detector1;

/*   detector_params=(LALFrDetector )LALMalloc(sizeof(LALFrDetector)); */
 
  bar=LALDETECTORTYPE_CYLBAR;
  strcpy(detector_params.name,"NAUTILUS");
  detector_params.vertexLongitudeRadians=12.67*LAL_PI/180.0;
  detector_params.vertexLatitudeRadians=41.82*LAL_PI/180.0;
  detector_params.vertexElevation=300.0;
  detector_params.xArmAltitudeRadians=0.0;
  detector_params.xArmAzimuthRadians=44.0*LAL_PI/180.0;

  LALCreateDetector(&status,&Detector1,&detector_params,bar);

  *Detector=Detector1;

  return 0;
}

/*******************************************************************************/

int Freemem() 
{

  INT4 i;

  /*Free SFTData*/

  for (i=0;i<GV.SFTno;i++)
    {
      LALFree(SFTData[i]->fft->data->data);
      LALFree(SFTData[i]->fft->data);
      LALFree(SFTData[i]->fft);
      LALFree(SFTData[i]);
    }
  LALFree(SFTData);

  /*Free timestamps*/
  LALFree(timestamps);

  LALFree(F);

  /*Free DemodParams*/

  LALSDestroyVector(&status, &(DemodParams->amcoe->a));
  LALSDestroyVector(&status, &(DemodParams->amcoe->b));

  LALFree(DemodParams->skyConst);
  LALFree(DemodParams->spinDwn);
  LALFree(DemodParams);

  LALCheckMemoryLeaks();
  
  return 0;
}
