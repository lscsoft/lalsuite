/*********************************************************************************/
/*                                 Noise vs. Time Code                           */
/*                                                                               */
/*			               X. Siemens                                */
/*                                                                               */
/*                                UWM - December 2002                            */
/*********************************************************************************/

#include "FreqAverager.h"
INT4 SFTno,RealSFTno;                /* Global variables that keep track of no of SFTs */
INT4 lalDebugLevel=3;
INT4 ifmin,band;
REAL4 *p;
REAL8 *po;
char filelist[MAXFILES][MAXFILENAMELENGTH];
REAL8 N,deltaT;

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

  fprintf(stderr,"Computing power spectral density:    ");
  if (ComputePSD(CommandLineArgs)) return 2;
  fprintf(stderr," Done\n");  

  /* Free memory*/
  fprintf(stderr,"Freeing allocated memory:             ");  
  if (Freemem()) return 8;
  fprintf(stderr," Done\n \n");  

  return 0;

}

/*******************************************************************************/

int ComputePSD(struct CommandLineArgsTag CLA)
{
  FILE *fp,*fpo;
  INT4 i,j=0,offset,ndeltaf;
  size_t errorcode;
  double Sh;
  char filename[256];

  strcpy(filename,CLA.outputfile);
  fpo=fopen(filename,"w");
  if (fpo==NULL) 
    {
      fprintf(stderr,"Could not open %s!\n",filename);
      return 1;
    }

  for (i=0;i<SFTno;i++)       /* Loop over SFTs          */
    {

      fprintf(stderr,"In file %d of %d\n",i,SFTno);

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
      
      offset=(ifmin-header.firstfreqindex)*2*sizeof(REAL4);
      errorcode=fseek(fp,offset,SEEK_CUR);
      ndeltaf=band+1;

      errorcode=fread((void*)p,2*ndeltaf*sizeof(REAL4),1,fp);  
      if (errorcode!=1)
	{
	  printf("Dang! Error reading data in SFT file %s!\n",filelist[i]);
	  return 1;
	}
      fclose(fp);

      /*Loop over frequency bins in each SFT      */
      for (j=0;j<ndeltaf;j++)
	 {
	   int jre=2*j;
	   int jim=jre+1;

	   po[i]=po[i]+(p[jre]*p[jre]+p[jim]*p[jim])/((REAL4) ndeltaf);
	 }
      Sh=2.0*deltaT/N * po[i];
      fprintf(fpo,"%d  %d  %d  %20.19e\n",i,header.gps_sec,header.gps_nsec,sqrt(Sh));
    }

  /* write output file */
      
  fclose(fpo);
  
  
  return 0;

}

/*******************************************************************************/

int ReadSFTDirectory(struct CommandLineArgsTag CLA)
{
  char command[256];
  size_t errorcode;
  FILE *fp;
  INT4 fileno=0,j;
  glob_t globbuf;


  /* check out what's in SFT directory */
  strcpy(command,CLA.directory);
  strcat(command,"/*SFT*");
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR|GLOB_MARK, NULL, &globbuf);

  /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */
  while (fileno < (int)globbuf.gl_pathc) 
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

  SFTno=fileno;  /* Global variable that keeps track of no of SFTs */


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
  po=(REAL8 *)LALMalloc(SFTno*sizeof(REAL8));

  for (j=0;j<SFTno;j++)
    {
      po[j]=0.0;   
    }
 
 deltaT=1.0/CLA.s;  /* deltaT is the time resolution of the original data */
 N=header.tbase*CLA.s; /* the number of time data points in one sft  */

 if(CLA.s == 0.0)
   {
     deltaT=header.tbase/(2.0*header.nsamples); /* deltaT is the time resolution of the original data */
     N=2.0*header.nsamples; /* the number of time data points in one sft */
   }

 fprintf(stdout,"deltaT=%f, N=%f\n",deltaT,N);


  ifmin=floor(CLA.f0*header.tbase);
  band=ceil(CLA.b*header.tbase);

      if (ifmin<header.firstfreqindex || 
	  ifmin+band > header.firstfreqindex+header.nsamples) 
	{
	  fprintf(stderr,"Freq index range %d->%d not in %d to %d\n",
		ifmin,ifmin+band,header.firstfreqindex,
		  header.firstfreqindex+header.nsamples);
	  return 4;
	}

  return 0;
}


/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{
  INT4 c, errflg = 0;
  optarg = NULL;
  
  /* Initialize default values */
  CLA->directory="";
  CLA->outputfile="";
  CLA->f0=0.0;
  CLA->b=0.0;
  CLA->s=0.0;

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"hb:D:r:I:C:o:f:b:s:"))!=-1))
    switch (c) {
    case 'D':
      /* SFT directory */
      CLA->directory=optarg;
      break;
    case 'o':
      /* calibrated sft output directory */
      CLA->outputfile=optarg;
      break;
    case 'f':
      /* calibrated sft output directory */
      CLA->f0=atof(optarg);
      break;
    case 'b':
      /* calibrated sft output directory */
      CLA->b=atof(optarg);
      break;
    case 's':
      /* sampling rate */
      CLA->s=atof(optarg);
      break;
   case 'h':
      /* print usage/help message */
      fprintf(stderr,"Arguments are:\n");
      fprintf(stderr,"\t-D\tSTRING\t(Directory where SFTs are located)\n");
      fprintf(stderr,"\t-o\tSTRING\t(Ascii file to output to)\n");
      fprintf(stderr,"\t-f\tFLOAT\t(Starting Frequency in Hz; default is 0.0 Hz)\n");
      fprintf(stderr,"\t-b\tFLOAT\t(Frequenzy band in Hz)\n");
      fprintf(stderr,"\t-s\tFLOAT\t(Sampling Rate in Hz; if = 0 then header info is used)\n");
      fprintf(stderr,"(eg: ./psd -D ../KnownPulsarDemo/data/ -o strain.txt -f 1010.7 -b 3.0) \n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }

  if(CLA->b == 0.0)
    {
      fprintf(stderr,"No frequency band specified; input band with -b option.\n");
      fprintf(stderr,"For help type ./psd -h \n");
      return 1;
    }      
  if(CLA->directory == "")
    {
      fprintf(stderr,"No directory specified; input directory with -D option.\n");
      fprintf(stderr,"For help type ./psd -h \n");
      return 1;
    }      
  if(CLA->outputfile == "")
    {
      fprintf(stderr,"No output directory specified; input directory with -o option.\n");
      fprintf(stderr,"For help type ./psd -h \n");
      return 1;
    }      
  return errflg;
}


/*******************************************************************************/


int Freemem()
{

  LALFree(po);
  LALFree(p);

  LALCheckMemoryLeaks();
  
  return 0;
}
