/*********************************************************************************/
/*                                 In.data file maker                            */
/*                                                                               */
/*			               X. Siemens                                */
/*          (takes in an Fstats file and a ParamMLE file and constructs In.data) */
/*                                                                               */
/*                                  UWM - May  2004                            */
/*********************************************************************************/

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <getopt.h>

int ReadCommandLine(int argc,char *argv[]);

char *fstatsfile,*parammlefile,*indatafile,*tsfile;
double Tsft,fstart,band;
int nSFT;

int main(int argc,char *argv[]) 
{
  int i;
  FILE *fp1;
  char line[256];
  double f0,alpha,delta,Ap,Ac,psi,phi0,dmp;

  if (ReadCommandLine(argc,argv)) return 1;

 /* ------ Open and read fstats file ------ */
 i=0;
 fp1=fopen(fstatsfile,"r");
 if (fp1==NULL) 
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",fstatsfile);
     return 1;
   }
 while(fgets(line,sizeof(line),fp1) && i == 0)
   {
     if (i > 1)
       {
	 fprintf(stderr,"Too many lines in file %s! Exiting... \n",fstatsfile);
	 return 1;
       }
     if(i==0)
       sscanf(line,"%le %le %le %le  %le %le %le",&f0,&alpha,&delta,&dmp,&dmp,&dmp,&dmp);
     i++;
   }
 fclose(fp1);     
 /* -- close fstats file -- */

 /* ------ Open and read parammle file ------ */
 i=0;
 fp1=fopen(parammlefile,"r");
 if (fp1==NULL) 
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",parammlefile);
     return 1;
   }
 while(fgets(line,sizeof(line),fp1))
   {
     if (i > 1)
       {
	 fprintf(stderr,"Too many lines in file %s! Exiting... \n",parammlefile);
	 return 1;
       }
     if(i==0)
       sscanf(line,"%le %le %le %le %le %le",&dmp,&dmp,&Ap,&Ac,&psi,&phi0);
     i++;
   }
 fclose(fp1);     
 /* -- close fstats file -- */

 /* ------ Open and write In.data file ------ */
 fp1=fopen(indatafile,"w");
 if (fp1==NULL) 
   {
     fprintf(stderr,"That's weird... could not open %s!\n",indatafile);
     return 1;
   }

 fprintf(fp1,"%le ##TSFT\n",Tsft);
 fprintf(fp1,"%d ##nSFT\n",nSFT);
 fprintf(fp1,"%le ##fstart\n",fstart);
 fprintf(fp1,"%le ##band\n",band);
 fprintf(fp1,"0.0 ##sigma\n");
 fprintf(fp1,"%le ##Ap\n",Ap);
 fprintf(fp1,"%le ##Ac\n",Ac);
 fprintf(fp1,"%le ##psi\n",psi);
 fprintf(fp1,"%le ##phi0\n",phi0);
 fprintf(fp1,"%le ##f0\n",f0);
 fprintf(fp1,"%le ##delta\n",delta);
 fprintf(fp1,"%le ##alpha\n",alpha);
 fprintf(fp1,"0 ##spinorder\n");
 fprintf(fp1,"%s ##timestampsfile\n",tsfile);


 fclose(fp1);     
 
  return 0;

}


/*******************************************************************************/


int ReadCommandLine(int argc,char *argv[]) 
{
  int errflg = 0;
  int c; 
  int option_index = 0;

  const char *optstring = "hf:p:t:o:l:n:s:b:";
  struct option long_options[] =
    {
      {"fstatsfile", 		required_argument, 0, 	'f'},
      {"param-mle", 		required_argument, 0, 	'p'},
      {"time-stamps", 		required_argument, 0, 	't'},
      {"outputfile", 		required_argument, 0, 	'o'},
      {"Tsft", 		        required_argument, 0, 	'l'},
      {"nSFT", 		        required_argument, 0, 	'n'},
      {"fstart", 		required_argument, 0, 	's'},
      {"band", 		        required_argument, 0, 	'b'},
      {"help", 			no_argument, 0, 	'h'},
      {0, 0, 0, 0}
    };

  
  fstatsfile=NULL;
  parammlefile=NULL;
  indatafile=NULL;
  tsfile=NULL;

  Tsft=0.0;
  nSFT=0;
  fstart=-1.0;
  band=0.0;

  /* Scan through list of command line arguments */
  while (1)
    {
      c = getopt_long(argc, argv, optstring, long_options, &option_index);
      
      if (c == -1) 
	break;

      switch (c) {
      case 'f':
	/* SFT directory */
	fstatsfile=optarg;
	break;
      case 'p':
	/* calibration files directory */
	parammlefile=optarg;
	break;
      case 'o':
	/* calibration files directory */
	indatafile=optarg;
	break;
      case 't':
	/* Spin down order */
	tsfile=optarg;
	break;
      case 'l':
	/* Spin down order */
	Tsft=atof(optarg);
	break;
      case 's':
	/* Spin down order */
	fstart=atof(optarg);
	break;
      case 'n':
	/* Spin down order */
	nSFT=atof(optarg);
	break;
      case 'b':
	/* Spin down order */
	band=atof(optarg);
	break;
      case 'h':
	/* print usage/help message */
	fprintf(stderr,"Arguments are (defaults):\n");
	fprintf(stderr,"\t--fstatsfile\t(-f)\tSTRING\tName of Fstats file\n");
	fprintf(stderr,"\t--param-mle\t(-p)\tSTRING\tName ParamMLE file\n");
	fprintf(stderr,"\t--outputfile\t(-o)\tSTRING\tName of ouput In.data file for makefakedata\n");
	fprintf(stderr,"\t--time-stamps\t(-t)\tSTRING\tName of time stamps file\n");
	fprintf(stderr,"\t--Tsft\t(-l)\tFLOAT\tTime baseline of sfts\n");
	fprintf(stderr,"\t--nSFT\t(-n)\tINTEGER\tNumber of sfts\n");
	fprintf(stderr,"\t--fstart\t(-s)\tFLOAT\tStarting freq. of sfts\n");
	fprintf(stderr,"\t--band\t(-b)\tFLOAT\tFreq. band of sfts\n");
	fprintf(stderr,"\t--help\t(-h)\t \tThis message\n");
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
  if(fstatsfile==NULL)
    {
      fprintf(stderr,"No Fstats file specified; input with -f option.\n");
      fprintf(stderr,"For help type ./makeInvetofile -h \n");
      return 1;
    }
  if(parammlefile==NULL)
    {
      fprintf(stderr,"No ParamMLE file specified; input with -p option.\n");
      fprintf(stderr,"For help type ./makeInvetofile -h \n");
      return 1;
    }
  if(indatafile==NULL)
    {
      fprintf(stderr,"No output In.data file specified; input with -o option.\n");
      fprintf(stderr,"For help type ./makeInvetofile -h \n");
      return 1;
    }
  if(tsfile==NULL)
    {
      fprintf(stderr,"No time-stamps file specified; input with -t option.\n");
      fprintf(stderr,"For help type ./makeInvetofile -h \n");
      return 1;
    }

  return errflg;
}
