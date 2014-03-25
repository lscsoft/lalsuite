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
 * \brief SFT Calibration Code
 */

/*********************************************************************************/
/*                                SFT Calibration Code                           */
/*                                                                               */
/*                                     X. Siemens                                */
/*                                                                               */
/*     (with much help from Bruce Allen, Patrick Brady and Jolien Creighton)     */
/*                                                                               */
/*                                UWM - November 2002                            */
/*********************************************************************************/

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <lal/LALDatatypes.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALDemod.h>

#define MAXFILES 70000          /* Maximum # of files in a directory */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename*/

#define MAXLINESRS   500000     /* Maximum # of lines in a Response or Sensing file */
#define MAXLINESF    100000     /* Maximum # of lines in a Factors file */

struct CommandLineArgsTag
{
  char *directory;
  char *caldirectory;
  char *run;
  char *IFO;
  char *outputdirectory;
} CommandLineArgs;

struct headertag
{
  REAL8 endian;
  INT4  gps_sec;
  INT4  gps_nsec;
  REAL8 tbase;
  INT4  firstfreqindex;
  INT4  nsamples;
} header;

typedef struct ResponseFunctionTag
{
  REAL4 Frequency[MAXLINESRS];
  REAL4 Magnitude[MAXLINESRS];
  REAL4 Phase[MAXLINESRS];
  REAL4 re[MAXLINESRS];
  REAL4 im[MAXLINESRS];
} Response;

typedef struct SensingFunctionTag
{
  REAL4 Frequency[MAXLINESRS];
  REAL4 Magnitude[MAXLINESRS];
  REAL4 Phase[MAXLINESRS];
  REAL4 re[MAXLINESRS];
  REAL4 im[MAXLINESRS];
} Sensing;

struct FactorsTag
{
  INT4  t[MAXLINESF];      /* in GPS seconds */
  REAL4 alpha[MAXLINESF];
  REAL4 alpha_beta[MAXLINESF];
} Factors;

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int ReadSFTDirectory(struct CommandLineArgsTag CLA);
int ReadCalibrationFiles(struct CommandLineArgsTag CLA);
int ComputeInitialRSFunctions(void);
int CalibrateSfts(struct CommandLineArgsTag CLA);
int Freemem(void);


COMPLEX8 tmpa, tmpb, tmpc;
REAL4 tmpx, tmpy;

INT4 SFTno,RealSFTno;
REAL4 *p,*pC;
char filelist[MAXFILES][MAXFILENAMELENGTH];
Sensing Sraw,So;
Response Rraw,Ro;

int Freemem(void);

extern char *optarg;
extern int optind, opterr, optopt;

int main(int argc,char *argv[])
{

  /* Reads command line arguments into the CommandLineArgs struct.
     In the absence of command line arguments it sets some defaults */
  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;

  /* Reads in SFT directory for filenames and total number of files */
  fprintf(stderr,"\n");
  fprintf(stderr,"Reading SFT directory:                ");
  if (ReadSFTDirectory(CommandLineArgs)) return 2;
  fprintf(stderr," Done\n");

  /* Reads Calibration Files */
  fprintf(stderr,"Reading Calibration Files:            ");
  if (ReadCalibrationFiles(CommandLineArgs)) return 3;
  fprintf(stderr," Done\n");

  /* Computes the initial response function for the frequency bins of the SFTs */
  fprintf(stderr,"Computing initial response function:  ");
  if (ComputeInitialRSFunctions()) return 4;
  fprintf(stderr," Done\n");

  /* Calibrates the SFTs */
  fprintf(stderr,"Calibrating %05d SFTs:                ",SFTno);
  if (CalibrateSfts(CommandLineArgs)) return 4;
  fprintf(stderr," Done\n");

  fprintf(stderr,"   (Info: Out of %d SFTs had calibration information for %d SFTs)\n",
          SFTno,RealSFTno);

  /* Free memory*/
  fprintf(stderr,"Freeing allocated memory:             ");
  if (Freemem()) return 8;
  fprintf(stderr," Done\n \n");

  return 0;

}


/*******************************************************************************/
int CalibrateSfts(struct CommandLineArgsTag CLA)
{
  FILE *fp,*fpo;
  INT4 i,j,t,k=0;
  REAL4 dt,a,b,alpha,alpha_beta;
  COMPLEX8 R,C,H;
  size_t errorcode;
  char filename[256],filenumber[16];

  RealSFTno=SFTno;

  for (i=0;i<SFTno;i++)       /* Loop over SFTs */
    {
      fp=fopen(filelist[i],"r");
      if (fp==NULL)
        {
          fprintf(stderr,"Weird... %s doesn't exist!\n",filelist[i]);
          return 1;
        }
      /* Read in the header from the file */
      errorcode=fread((void*)&header,sizeof(header),1,fp);
      if (errorcode!=1)
        {
          fprintf(stderr,"No header in data file %s\n",filelist[i]);
          return 1;
        }

      /* Check that data is correct endian order */
      if (header.endian!=1.0)
        {
          fprintf(stderr,"First object in file %s is not (double)1.0!\n",filelist[i]);
          fprintf(stderr,"It could be a file format error (big/little\n");
          fprintf(stderr,"endian) or the file might be corrupted\n\n");
          return 2;
        }

      /* Check that the time base is positive */
      if (header.tbase<=0.0)
        {
          fprintf(stderr,"Timebase %f from data file %s non-positive!\n",
                  header.tbase,filelist[i]);
          return 3;
        }

      t=header.gps_sec;

      /* If the SFT time is smaller than the 1st time of the calibration
       then there is no calibration info for that file and it should be skipped */
      if(t < Factors.t[0] )
        {
          RealSFTno=RealSFTno-1;
          fclose(fp);
          fprintf(stderr,"Threw out %s\n",filelist[i]);
          fprintf(stderr,"Time of SFT smaller than first factors file time\n");
          continue;
        }

      /* Advance in time until sft time is smaller than Factors time */
      while ( k < MAXLINESF-1 && t > Factors.t[k])
        {
          k++;
        }

      /* Since now Factors.t[k-1] < t =< Factors.t[k] ||*/
      /* can compute a weighted average of the factors at i-1 and i */

      /* check both bounds!!!!! */
      if(t<Factors.t[k-1] || t > Factors.t[k])
        {
          fprintf(stderr,"Time of SFT does not lie between two lines in factors file!\n");
          fprintf(stderr,"Time of SFT is: %d; Factors times are %d and %d\n",t,Factors.t[k-1],Factors.t[k]);
          return 1;
        }


      /* But first, if any factors are = 0 we don't know what the calibration is
         and we should skip the SFT */
      if(Factors.alpha_beta[k]==0 && Factors.alpha_beta[k-1]==0)
        {
              RealSFTno=RealSFTno-1;
              fclose(fp);
              fprintf(stderr,"Threw out %s\n",filelist[i]);
    fprintf(stderr,"Factor %d=%f  Factor %d=%f\n",k,Factors.alpha_beta[k],k-1,Factors.alpha_beta[k-1]);
              continue;
        }
      if(Factors.alpha_beta[k]==0 || Factors.alpha_beta[k-1]==0)
        {
          if(t != Factors.t[k-1] && t != Factors.t[k])
            {
              RealSFTno=RealSFTno-1;
              fclose(fp);
              fprintf(stderr,"Threw out %s\n",filelist[i]);
    fprintf(stderr,"Factor %d=%f  Factor %d=%f\n",k,Factors.alpha_beta[k],k-1,Factors.alpha_beta[k-1]);
              continue;
            }
        }

      dt=(REAL4)(Factors.t[k]-Factors.t[k-1]);

      if(dt <= 0)
        {
          fprintf(stderr,"Time interval smaller or equal to zero, dt= %f !\n",dt);
          return 1;
        }


     /* If the time interval between two points in the file is > 60 that means there's no cal info */
      if(dt > 60.0)
        {
          if(t != Factors.t[k-1] && t != Factors.t[k])
            {
              RealSFTno=RealSFTno-1;
              fclose(fp);
              fprintf(stderr,"Threw out %s\n",filelist[i]);
              fprintf(stderr,"Time interval grater than 60s dt= %f !\n",dt);
              fprintf(stderr,"%d    %d\n",Factors.t[k],Factors.t[k-1]);
              continue;
            }
        }

      a=(REAL4)(t-Factors.t[k-1]) / dt;
      if (a>1.0) a=1.0;
      b=1.0-a;

/*       fprintf(stderr,"----> Calibrating %s\n",filelist[i]); */

/*       fprintf(stderr,"a=%f  b=%f\n",a,b); */

      alpha=a*Factors.alpha[k]+b*Factors.alpha[k-1];
      alpha_beta=a*Factors.alpha_beta[k]+b*Factors.alpha_beta[k-1];

  /*     fprintf(stderr,"alpha=%f  alpha_beta=%f\n",alpha,alpha_beta); */


      if(alpha == 0.0  ||  alpha_beta == 0.0)
        {
              RealSFTno=RealSFTno-1;
              fclose(fp);
              fprintf(stderr,"Threw out %s\n",filelist[i]);
              fprintf(stderr,"alpha=%f alpha_beta=%f!\n",alpha,alpha_beta);
              fprintf(stderr,"%d    %d\n",Factors.t[k],Factors.t[k-1]);
              continue;
        }

      if(alpha < 0.3 )
        {
              RealSFTno=RealSFTno-1;
              fclose(fp);
              fprintf(stderr,"Threw out %s\n",filelist[i]);
              fprintf(stderr,"alpha=%f too small!\n",alpha);
              fprintf(stderr,"%d    %d\n",Factors.t[k],Factors.t[k-1]);
              continue;
        }

      if(alpha > 2.0 )
        {
              RealSFTno=RealSFTno-1;
              fclose(fp);
              fprintf(stderr,"Threw out %s\n",filelist[i]);
              fprintf(stderr,"alpha=%f too large!\n",alpha);
              fprintf(stderr,"%d    %d\n",Factors.t[k],Factors.t[k-1]);
              continue;
        }

 /*      fprintf(stderr,"%d    %d\n",Factors.t[k],Factors.t[k-1]); */
/*       fprintf(stderr,"%f    %f\n",Factors.alpha_beta[k],Factors.alpha_beta[k-1]); */

      strcpy(filename,CLA.outputdirectory);
      strcat(filename,"/CAL_SFT.");
      sprintf(filenumber,"%09d",header.gps_sec);
      strcat(filename,filenumber);

      fpo=fopen(filename,"w");
      if (fpo==NULL)
        {
          fprintf(stderr,"Problems opening file %s !\n",filename);
          return 1;
        }

      /* first thing is to write the header into it */
      errorcode=fwrite((void*)&header,sizeof(header),1,fpo);
      if (errorcode!=1)
        {
        printf("Error in writing header into file %s!\n",filename);
        return 1;
        }

      errorcode=fread((void*)p,2*header.nsamples*sizeof(REAL4),1,fp);
      if (errorcode!=1)
        {
          printf("Dang! Error reading data in SFT file %s!\n",filename);
          return 1;
        }
      fclose(fp);

      /* then we read points into from original sft file, calibrate them and write them
         into the new file */

       /* Loop over frequency bins in each SFT       */
       for (j=0;j<header.nsamples;j++)
         {
           int jre=2*j;
           int jim=jre+1;

           R = crectf( Ro.re[j], Ro.im[j] );

           C = crectf( So.re[j], So.im[j] );

           /* compute the reference open loop function H0 */
           H = (C * R) - 1.0;

           /* update the open loop function */
           H *= ((REAL4) alpha_beta);

           /* update the sensing function */
           C *= ((REAL4) alpha);

           /* compute the updated response function */
           H = crectf( crealf(H) + 1.0, cimagf(H) );
           R = H / C;

           /* the jth elements of p and pC are th real parts and the (j+1)th the imaginary part */
           pC[jre]=crealf(R)*p[jre]-cimagf(R)*p[jim];
           pC[jim]=crealf(R)*p[jim]+cimagf(R)*p[jre];
       }

     errorcode=fwrite((void*)pC,2*header.nsamples*sizeof(REAL4),1,fpo);
      if (errorcode!=1){
        printf("Error in writing data into SFT file!\n");
        return 1;
      }
      fclose(fpo);
    }

  return 0;

}

/*******************************************************************************/

int ComputeInitialRSFunctions(void)
{
  FILE *fp;
  INT4 i,j;
  REAL4 f,df,a,b;
  size_t errorcode;

  /* open FIRST file and get info from it*/

  fp=fopen(filelist[0],"r");
  if (fp==NULL)
    {
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

  /* Check that data is correct endian order */
  if (header.endian!=1.0)
    {
      fprintf(stderr,"First object in file %s is not (double)1.0!\n",filelist[0]);
      fprintf(stderr,"It could be a file format error (big/little\n");
      fprintf(stderr,"endian) or the file might be corrupted\n\n");
      return 2;
    }

  /* Check that the time base is positive */
  if (header.tbase<=0.0)
    {
      fprintf(stderr,"Timebase %f from data file %s non-positive!\n",
              header.tbase,filelist[0]);
      return 3;
    }
  fclose(fp);

  /* Allocate pointers for SFT data -- to be used later */
  p=(REAL4 *)LALMalloc(2*header.nsamples*sizeof(REAL4));
  pC=(REAL4 *)LALMalloc(2*header.nsamples*sizeof(REAL4));

  i=0;

  for (j=0;j<header.nsamples;j++)
    {
      f = header.firstfreqindex/header.tbase + j/header.tbase;

      while (i < MAXLINESRS-1 && f > Rraw.Frequency[i]) i++;

      if(i == MAXLINESRS-1 && f > Rraw.Frequency[i])
        {
          fprintf(stderr,"No calibration info for frequency %f!\n",f);
          return 1;
        }
      /* checks order */

      /* Since now Rraw.Frequency[i-1] < f =< Rraw.Frequency[i] ||*/
      /* can compute a weighted average of the raw frequencies at i-1 and i */

      /* check both bounds!!!!! */
      if(f < Rraw.Frequency[i-1] || f > Rraw.Frequency[i])
        {
          fprintf(stderr,"Frequency %f in SFT does not lie between two lines in Response file!\n",f);
          return 1;
        }

      /* If the frequencies are closely spaced this may create dangerous floating point errors */
      df=Rraw.Frequency[i]-Rraw.Frequency[i-1];

      a=(f-Rraw.Frequency[i-1])/df;
      if (a>1.0) a=1.0;
      b=1.0-a;

      Ro.Frequency[j]=f;
      Ro.Magnitude[j]=a*Rraw.Magnitude[i]+b*Rraw.Magnitude[i-1];
      Ro.Phase[j]=a*Rraw.Phase[i]+b*Rraw.Phase[i-1];

      Ro.re[j]=a*Rraw.re[i]+b*Rraw.re[i-1];
      Ro.im[j]=a*Rraw.im[i]+b*Rraw.im[i-1];


      So.Frequency[j]=f;
      So.Magnitude[j]=a*Sraw.Magnitude[i]+b*Sraw.Magnitude[i-1];
      So.Phase[j]=a*Sraw.Phase[i]+b*Sraw.Phase[i-1];

      So.re[j]=a*Sraw.re[i]+b*Sraw.re[i-1];
      So.im[j]=a*Sraw.im[i]+b*Sraw.im[i-1];
    }

  return 0;
}

/*******************************************************************************/

int ReadCalibrationFiles(struct CommandLineArgsTag CLA)
{
  char SensingFile[256],ResponseFile[256],FactorsFile[256];
  char line[256];

  INT4 i;
  FILE *fpS,*fpR,*fpF;

 strcpy(SensingFile,CLA.caldirectory);
 strcat(SensingFile,"/"); strcat(SensingFile,CLA.run);
 strcat(SensingFile,"-"); strcat(SensingFile,CLA.IFO);
 strcat(SensingFile,"-CAL-CAV_GAIN.txt");

 strcpy(ResponseFile,CLA.caldirectory);
 strcat(ResponseFile,"/"); strcat(ResponseFile,CLA.run);
 strcat(ResponseFile,"-"); strcat(ResponseFile,CLA.IFO);
 strcat(ResponseFile,"-CAL-RESPONSE.txt");

 strcpy(FactorsFile,CLA.caldirectory);
 strcat(FactorsFile,"/");strcat(FactorsFile,CLA.run);
 strcat(FactorsFile,"-");strcat(FactorsFile,CLA.IFO);
 strcat(FactorsFile,"-CAL-FACTORS.txt");


 /* This is kinda messy... Unfortunately there's no good way of doing this */

 /* ------ Open and read Sensing file ------ */
 i=0;
 fpS=fopen(SensingFile,"r");
 if (fpS==NULL)
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",SensingFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpS))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i >= MAXLINESRS)
       {
         fprintf(stderr,"Too many lines in file %s! Exiting... \n", SensingFile);
         return 1;
       }
     sscanf(line,"%e %e %e",&Sraw.Frequency[i],&Sraw.Magnitude[i],&Sraw.Phase[i]);
     Sraw.re[i]=Sraw.Magnitude[i]*cos(Sraw.Phase[i]);
     Sraw.im[i]=Sraw.Magnitude[i]*sin(Sraw.Phase[i]);

     i++;
   }
 fclose(fpS);
 /* -- close Sensing file -- */

 /* ------ Open and read Response File ------ */
 i=0;
 fpR=fopen(ResponseFile,"r");
 if (fpR==NULL)
   {
     fprintf(stderr,"Weird... %s doesn't exist!\n",ResponseFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpR))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i >= MAXLINESRS)
       {
         fprintf(stderr,"Too many lines in file %s! Exiting... \n", ResponseFile);
         return 1;
       }
     sscanf(line,"%e %e %e",&Rraw.Frequency[i],&Rraw.Magnitude[i],&Rraw.Phase[i]);
     Rraw.re[i]=Rraw.Magnitude[i]*cos(Rraw.Phase[i]);
     Rraw.im[i]=Rraw.Magnitude[i]*sin(Rraw.Phase[i]);
     i++;
   }
 fclose(fpR);
 /* -- close Response file -- */

 /* ------ Open and read Factors file ------ */
 i=0;
 fpF=fopen(FactorsFile,"r");
 if (fpF==NULL)
   {
     fprintf(stderr,"Weird... %s doesn't exist!\n",FactorsFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpR))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i >= MAXLINESF)
       {
         fprintf(stderr,"Too many lines in file %s! Exiting... \n", FactorsFile);
         return 1;
       }
     sscanf(line,"%d %e %e",&Factors.t[i],&Factors.alpha_beta[i],&Factors.alpha[i]);
     i++;
   }
 fclose(fpF);
 /* -- close Factors file -- */

  return 0;
}

/*******************************************************************************/

int ReadSFTDirectory(struct CommandLineArgsTag CLA)
{
  char command[256];
  INT4 filenum=0;
  glob_t globbuf;


  /* check out what's in SFT directory */
  strcpy(command,CLA.directory);
  strcat(command,"/*SFT*");
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR|GLOB_MARK, NULL, &globbuf);

  /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */
  while (filenum < (int) globbuf.gl_pathc)
    {
      strcpy(filelist[filenum],globbuf.gl_pathv[filenum]);
      filenum++;
      if (filenum > MAXFILES)
        {
          fprintf(stderr,"Too many files in directory! Exiting... \n");
          return 1;
        }
    }
  globfree(&globbuf);

  SFTno=filenum;  /* Global variable that keeps track of no of SFTs */

  return 0;
}


/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA)
{
  INT4 c, errflg = 0;
  optarg = NULL;

  /* Initialize default values */
  CLA->directory=NULL;
  CLA->run=NULL;
  CLA->IFO=NULL;
  CLA->caldirectory=NULL;
  CLA->outputdirectory=NULL;

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"hb:D:r:I:C:o:"))!=-1))
    switch (c) {
    case 'D':
      /* SFT directory */
      CLA->directory=optarg;
      break;
    case 'C':
      /* calibration files directory */
      CLA->caldirectory=optarg;
      break;
    case 'r':
      /* run */
      CLA->run=optarg;
      break;
    case 'I':
      /* interferometer */
      CLA->IFO=optarg;
      break;
    case 'o':
      /* calibrated sft output directory */
      CLA->outputdirectory=optarg;
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stderr,"Arguments are:\n");
      fprintf(stderr,"\t-D\tSTRING\t(Directory where SFTs are located)\n");
      fprintf(stderr,"\t-C\tSTRING\t(Directory where calibration files are located)\n");
      fprintf(stderr,"\t-r\tSTRING\t(Run (ie: E7, S1, S2, etc.)\n");
      fprintf(stderr,"\t-I\tSTRING\t(Interferometer (ie: L1, H1, H2 etc.)\n");
      fprintf(stderr,"\t-o\tSTRING\t(Ouput directory for calibrated SFTs)\n");
      fprintf(stderr,"(eg: ./CalibrateSFTs -D ../KnownPulsarDemo/data/ -C ./ -r S1 -I L1 -o ./) \n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }

  if(CLA->directory == NULL)
    {
      fprintf(stderr,"No directory specified; input directory with -D option.\n");
      fprintf(stderr,"For help type ./CalibrateSFTs -h \n");
      return 1;
    }
  if(CLA->caldirectory == NULL)
    {
      fprintf(stderr,"No calibration directory specified; input directory with -C option.\n");
      fprintf(stderr,"For help type ./CalibrateSFTs -h \n");
      return 1;
    }
  if(CLA->outputdirectory == NULL)
    {
      fprintf(stderr,"No output directory specified; input directory with -o option.\n");
      fprintf(stderr,"For help type ./CalibrateSFTs -h \n");
      return 1;
    }
  if(CLA->run == NULL)
    {
      fprintf(stderr,"No run specified; input run with -r option.\n");
      fprintf(stderr,"For help type ./CalibrateSFTs -h \n");
      return 1;
    }
  if(CLA->IFO == NULL)
    {
      fprintf(stderr,"No interferometer specified; input interferometer with -I option.\n");
      fprintf(stderr,"For help type ./CalibrateSFTs -h \n");
      return 1;
    }

  return errflg;
}


/*******************************************************************************/


int Freemem(void)
{


  LALFree(pC);
  LALFree(p);

  LALCheckMemoryLeaks();

  return 0;
}
