/*********************************************************************************/
/*                     polka - the pulsar koinzidenz analysis code               */
/*                                                                               */
/*			               X. Siemens                                */
/*                   (takes in two Fstats file to look for coincidence)          */
/*                                                                               */
/*                                  UWM - March  2004                            */
/*********************************************************************************/

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <getopt.h>
#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>

#define MAXCANDIDATES   5000000     /* Maximum # of allowed candidates */
#define MAXCOINC   5000000     /* Maximum # of allowed coincident candidates */

struct CommandLineArgsTag 
{
  char *FstatsFile1; /* Names of Fstat files to be read in */
  char *FstatsFile2;
  char *FstatsFile3; /* Names of Fstat files to be read in to compute false alarm */
  char *FstatsFile4;
  char *OutputFile;
  REAL8 Deltaf;      /* Size of coincidence window in Hz */
  REAL8 DeltaAlpha;  /* Size of coincidence window in radians */
  REAL8 DeltaDelta;  /* Size of coincidence window in radians */
  REAL8 fmin;        /* Minimum frequency of candidate in first IFO */
  REAL8 fmax;        /* Maximum frequency of candidate in first IFO */
} CommandLineArgs;

typedef struct CandidateTag 
{
  REAL8 f[MAXCANDIDATES];        /* Frequency */
  REAL8 Alpha[MAXCANDIDATES];    /* longitude */
  REAL8 Delta[MAXCANDIDATES];    /* latitude */
  REAL8 F[MAXCANDIDATES];        /* Maximum value of F for the cluster */
  REAL8 fa[MAXCANDIDATES];       /* false alarm probability for that candidate */
} Candidate;

typedef struct CoincidentCandidateTag 
{
  REAL8 f1[MAXCOINC],f2[MAXCOINC];            /* Frequency */
  REAL8 Alpha1[MAXCOINC],Alpha2[MAXCOINC];    /* longitude */
  REAL8 Delta1[MAXCOINC],Delta2[MAXCOINC];    /* latitude */
  REAL8 F1[MAXCOINC],F2[MAXCOINC];        /* Maximum value of F for the cluster */
  REAL8 fa[MAXCOINC],fa1[MAXCOINC],fa2[MAXCOINC];       /* false alarm probability for that candidate */
} CoincidentCandidate;

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int ReadCandidateFiles(struct CommandLineArgsTag CLA);
int compare1F(const void *ip, const void *jp);
int compare2F(const void *ip, const void *jp);
int compare3F(const void *ip, const void *jp);
int compare4F(const void *ip, const void *jp);
int compare2f(const void *ip, const void *jp);
int compareCCfa(const void *ip, const void *jp);
void locate(double xx[], int n, double x, int *j, int *indices);

int NCands1,NCands2,NCCands,NCands3,NCands4;       /* Global variables that keep track of no of candidates */
INT4 lalDebugLevel=0;
Candidate C1,C2,C3,C4; /* Candidate structures */
CoincidentCandidate CC;

int main(int argc,char *argv[]) 
{
  INT4 *indices1F=NULL,*indices2f=NULL,*indices2F=NULL,*indicesCCfa=NULL,*indices3F=NULL,*indices4F=NULL;
  REAL8 MaxAngularDistance;
  int i,k;
  FILE *fpOut;
 
  /* Reads command line arguments */
  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;

  /* Reads in candidare files */
  if (ReadCandidateFiles(CommandLineArgs)) return 2;

  /* create arrays of indices */
  if (!(indices1F=(INT4 *)LALMalloc(sizeof(INT4)*NCands1))){
    fprintf(stderr,"Unable to allocate index array in main\n");
    return 1;
  }
  if (!(indices2F=(INT4 *)LALMalloc(sizeof(INT4)*NCands2))){
    fprintf(stderr,"Unable to allocate index array in main\n");
    return 1;
  }
  if (!(indices2f=(INT4 *)LALMalloc(sizeof(INT4)*NCands2))){
    fprintf(stderr,"Unable to allocate index array in main\n");
    return 1;
  }

  /* populate arrays of indices */
  for (i=0;i<NCands1;i++) indices1F[i]=i;
  for (i=0;i<NCands2;i++) indices2F[i]=i;
  for (i=0;i<NCands2;i++) indices2f[i]=i;

  /* sort arrays of indices in DECREASING order*/
  qsort((void *)indices1F, (size_t)NCands1, sizeof(int), compare1F);
  qsort((void *)indices2F, (size_t)NCands2, sizeof(int), compare2F);
  qsort((void *)indices2f, (size_t)NCands2, sizeof(int), compare2f);

  if(CommandLineArgs.FstatsFile3 != NULL)
    {
      if (!(indices3F=(INT4 *)LALMalloc(sizeof(INT4)*NCands3))){
	fprintf(stderr,"Unable to allocate index array in main\n");
	return 1;
      }
      for (i=0;i<NCands3;i++) indices3F[i]=i;
      qsort((void *)indices3F, (size_t)NCands3, sizeof(int), compare3F);
    }      
  if(CommandLineArgs.FstatsFile4 != NULL)
    {
      if (!(indices4F=(INT4 *)LALMalloc(sizeof(INT4)*NCands4))){
	fprintf(stderr,"Unable to allocate index array in main\n");
	return 1;
      }
      for (i=0;i<NCands4;i++) indices4F[i]=i;
      qsort((void *)indices4F, (size_t)NCands4, sizeof(int), compare4F);
    }      

  k=0; /* kounts koinzident events */
  MaxAngularDistance=sqrt(pow(CommandLineArgs.DeltaAlpha,2)+pow(CommandLineArgs.DeltaDelta,2))+1e-8;

  /* go through list */
  for (i=0; i < NCands1; i++)
    {
      REAL8 f1min,f1max,difff;
      REAL8 f1,Alpha1,Delta1,F1;
      int if2min,if2max,f;

      /* Minimum and maximum frequencies acceptable for coincidence */
      f1=C1.f[indices1F[i]];
      /* if candidate frequency does not lie within bounds specified by user go to next in list */
      if(f1 < CommandLineArgs.fmin || f1 > CommandLineArgs.fmax) continue;

      f1min=f1-CommandLineArgs.Deltaf;
      f1max=f1+CommandLineArgs.Deltaf;

      /* Find nearest index to f1min and f1max; function explained below */
      locate(C2.f,NCands2,f1min,&if2min,indices2f);
      locate(C2.f,NCands2,f1max,&if2max,indices2f);

      /* alpha */
      Alpha1=C1.Alpha[indices1F[i]];
      /* delta */
      Delta1=C1.Delta[indices1F[i]];
      /* F */
      F1=C1.F[indices1F[i]];
      
      for (f=if2max; f <= if2min; f++)
	{
	  REAL8 Alpha2=C2.Alpha[indices2f[f]],Delta2=C2.Delta[indices2f[f]];
	  REAL8 n1[3],n2[3],AngularDistance;
	  
	  n1[0]=cos(Alpha1)*cos(Delta1);
	  n1[1]=sin(Alpha1)*cos(Delta1);
	  n1[2]=sin(Delta1);

	  n2[0]=cos(Alpha2)*cos(Delta2);
	  n2[1]=sin(Alpha2)*cos(Delta2);
	  n2[2]=sin(Delta2);

	  AngularDistance=acos(n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2]);
	  difff=fabs(f1 - C2.f[indices2f[f]]);

	  /* check difference in frequencies because we're not guaranteed 
	     sufficient closeness at the edges of array */
	  if ( difff <= CommandLineArgs.Deltaf) 
	    {
	      if ( AngularDistance <= MaxAngularDistance )
		{	
		  int j;
		      
		  CC.f1[k]=f1;
		  CC.Alpha1[k]=Alpha1;
		  CC.Delta1[k]=Delta1;
		  CC.F1[k]=F1;

		  if(CommandLineArgs.FstatsFile3 != NULL)
		    {
		      locate(C3.F,NCands3,CC.F1[k],&j,indices3F);
		      CC.fa1[k]=(REAL8)(j+1)/(REAL8)NCands3;
		    }
		  else CC.fa1[k]=(REAL8)(i+1)/(REAL8)NCands1;

		  CC.f2[k]=C2.f[indices2f[f]];
		  CC.Alpha2[k]=C2.Alpha[indices2f[f]];
		  CC.Delta2[k]=C2.Delta[indices2f[f]];
		  CC.F2[k]=C2.F[indices2f[f]];

		  if(CommandLineArgs.FstatsFile4 != NULL)
		    {
		      locate(C4.F,NCands4,CC.F2[k],&j,indices4F);
		      CC.fa2[k]=(REAL8)(j+1)/(REAL8)NCands4;
		    }else{
		      locate(C2.F,NCands2,CC.F2[k],&j,indices2F);
		      CC.fa2[k]=(REAL8)(j+1)/(REAL8)NCands2;
		    }

		  CC.fa[k]=CC.fa1[k]*CC.fa2[k];

		  k++;
		  
		}
	    }
	}
      /* next candidate for 1st ifo */
    }     

  NCCands=k;
  /* allocate space */
  if (!(indicesCCfa=(INT4 *)LALMalloc(sizeof(INT4)*NCCands))){
    fprintf(stderr,"Unable to allocate index array in main\n");
    return 1;
  }
  for (i=0;i<NCCands;i++) indicesCCfa[i]=i;
  /* sort in increasing probability of joint false alarm */
  qsort((void *)indicesCCfa, (size_t)NCCands, sizeof(int), compareCCfa);

  /* open and write the file */
  fpOut=fopen(CommandLineArgs.OutputFile,"w"); 	 
  for (i=0;i<NCCands;i++) 
    {
      k=indicesCCfa[i];
      fprintf(fpOut,"%1.15le %le %le %le %le %1.15le %le %le %le %le %le\n",
	      CC.f1[k],CC.Alpha1[k],CC.Delta1[k],
	      CC.F1[k],CC.fa1[k],
	      CC.f2[k],CC.Alpha2[k],CC.Delta2[k],
	      CC.F2[k],CC.fa2[k],CC.fa[k]);
    }
  fclose(fpOut);

  LALFree(indices1F);
  LALFree(indices2F);
  LALFree(indices2f);
  LALFree(indicesCCfa);

  return 0;

}

/*******************************************************************************/

void locate(double xx[], int n, double x, int *j, int *indices) 
     /* locates x in array of xx */
{ 

 int ju,jm,jl; 


 if( x <= xx[indices[n-1]] ) 
   {
     *j=n-1;
     return;
   }
 if( x >= xx[indices[0]] ) 
   {
     *j=0;
     return;
   }

 jl=0;  
 ju=n;

 while (ju-jl > 1) 
   {
     jm=(ju+jl)/2;
     if ( x <= xx[indices[jm]] ) 
       jl=jm;
     else ju=jm; 
   }
 
 *j=jl;

}


/*******************************************************************************/

/* Sorting function to sort 1st candidate into DECREASING order of F */
int compare1F(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=C1.F[*(int *)ip];
  dj=C1.F[*(int *)jp];

  if (di<dj)
    return 1;
  
  if (di==dj)
    return (ip < jp);

  return -1;
}
/*******************************************************************************/

/* Sorting function to sort 1st candidate into DECREASING order of F */
int compare2F(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=C2.F[*(int *)ip];
  dj=C2.F[*(int *)jp];

  if (di<dj)
    return 1;
  
  if (di==dj)
    return (ip < jp);

  return -1;
}

/*******************************************************************************/

/* Sorting function to sort 1st candidate into DECREASING order of F */
int compare3F(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=C3.F[*(int *)ip];
  dj=C3.F[*(int *)jp];

  if (di<dj)
    return 1;
  
  if (di==dj)
    return (ip < jp);

  return -1;
}
/*******************************************************************************/

/* Sorting function to sort 1st candidate into DECREASING order of F */
int compare4F(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=C4.F[*(int *)ip];
  dj=C4.F[*(int *)jp];

  if (di<dj)
    return 1;
  
  if (di==dj)
    return (ip < jp);

  return -1;
}

/*******************************************************************************/

/* Sorting function to sort second candidate list into DECREASING order of f */
int compare2f(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=C2.f[*(int *)ip];
  dj=C2.f[*(int *)jp];

  if (di<dj)
    return 1;
  
  if (di==dj)
    return (ip < jp);

  return -1;
}
/*******************************************************************************/

/* Sorting function to sort second candidate list into increasing order of fa */
int compareCCfa(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=CC.fa[*(int *)ip];
  dj=CC.fa[*(int *)jp];

  if (di<dj)
    return -1;
  
  if (di==dj)
    return (ip > jp);

  return 1;
}


/*******************************************************************************/

int ReadCandidateFiles(struct CommandLineArgsTag CLA)
{
  INT4 i;
  FILE *fp;
  char line[256];
  REAL8 dmp;


 /* This is kinda messy... Unfortunately there's no good way of doing this */
 
 /* ------ Open and read 1st candidate file ------ */
 i=0;
 fp=fopen(CLA.FstatsFile1,"r");
 if (fp==NULL) 
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.FstatsFile1);
     return 1;
   }
 while(fgets(line,sizeof(line),fp))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i >= MAXCANDIDATES)
       {
	 fprintf(stderr,"Too many lines in file %s! Exiting... \n",CLA.FstatsFile1);
	 return 1;
       }
     sscanf(line,"%le %le %le %le %le %le %le",&C1.f[i],&C1.Alpha[i],&C1.Delta[i],
	    &dmp,&dmp,&dmp,&C1.F[i]);
     i++;
   }
 NCands1=i;
 fclose(fp);     
 /* -- close 1st candidate file -- */

 /* ------ Open and read 2nd candidate file ------ */
 i=0;
 fp=fopen(CLA.FstatsFile2,"r");
 if (fp==NULL) 
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.FstatsFile2);
     return 1;
   }
 while(fgets(line,sizeof(line),fp))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i >= MAXCANDIDATES)
       {
	 fprintf(stderr,"Too many lines in file %s! Exiting... \n",CLA.FstatsFile2);
	 return 1;
       }
     sscanf(line,"%le %le %le %le %le %le %le",&C2.f[i],&C2.Alpha[i],&C2.Delta[i],
	    &dmp,&dmp,&dmp,&C2.F[i]);
     i++;
   }
 NCands2=i;
 fclose(fp);     
 /* -- close 1st candidate file -- */

 if(CLA.FstatsFile3 != NULL)
   {
     /* ------ Open and read 3rd candidate file ------ */
     i=0;
     fp=fopen(CLA.FstatsFile3,"r");
     if (fp==NULL) 
       {
	 fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.FstatsFile3);
	 return 1;
       }
     while(fgets(line,sizeof(line),fp))
       {
	 if(*line == '#') continue;
	 if(*line == '%') continue;
	 if (i >= MAXCANDIDATES)
	   {
	     fprintf(stderr,"Too many lines in file %s! Exiting... \n",CLA.FstatsFile3);
	     return 1;
	   }
	 sscanf(line,"%le %le %le %le %le %le %le",&C3.f[i],&C3.Alpha[i],&C3.Delta[i],
		&dmp,&dmp,&dmp,&C3.F[i]);
	 i++;
       }
     NCands3=i;
     fclose(fp);     
    }      

 if(CLA.FstatsFile4 != NULL)
   {
     /* ------ Open and read 4th candidate file ------ */
     i=0;
     fp=fopen(CLA.FstatsFile4,"r");
     if (fp==NULL) 
       {
	 fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.FstatsFile4);
	 return 1;
       }
     while(fgets(line,sizeof(line),fp))
       {
	 if(*line == '#') continue;
	 if(*line == '%') continue;
	 if (i >= MAXCANDIDATES)
	   {
	     fprintf(stderr,"Too many lines in file %s! Exiting... \n",CLA.FstatsFile4);
	     return 1;
	   }
	 sscanf(line,"%le %le %le %le %le %le %le",&C4.f[i],&C4.Alpha[i],&C4.Delta[i],
		&dmp,&dmp,&dmp,&C4.F[i]);
	 i++;
       }
     NCands4=i;
     fclose(fp);     
    }      


  return 0;
}

/*******************************************************************************/


int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{
  INT2 errflg = 0;
  INT4 c; 
  INT4 option_index = 0;

  const char *optstring = "h1:2:3:4:f:a:d:m:M:o:";
  struct option long_options[] =
    {
      {"fstatsfile1", 		required_argument, 0, 	'1'},
      {"fstatsfile2", 		required_argument, 0, 	'2'},
      {"frequency-window", 	required_argument, 0, 	'f'},
      {"delta-window", 		required_argument, 0, 	'd'},
      {"alpha-window", 		required_argument, 0, 	'a'},
      {"fmin",   		required_argument, 0, 	's'},
      {"fmax",   		required_argument, 0, 	'e'},
      {"outputfile", 		required_argument, 0, 	'o'},
      {"help", 			no_argument, 0, 	'h'},
      {0, 0, 0, 0}
    };

  /* Initialize default values */
  CLA->FstatsFile1=NULL;
  CLA->FstatsFile2=NULL;
  CLA->FstatsFile3=NULL;
  CLA->FstatsFile4=NULL;
  CLA->OutputFile=NULL;
  CLA->Deltaf=0.0;
  CLA->DeltaAlpha=0;
  CLA->DeltaDelta=0;
  CLA->fmin=0;
  CLA->fmax=0;


  /* Scan through list of command line arguments */
  while (1)
    {
      c = getopt_long(argc, argv, optstring, long_options, &option_index);
      
      if (c == -1) 
	break;

      switch (c) {
      case '1':
	/* SFT directory */
	CLA->FstatsFile1=optarg;
	break;
      case '2':
	/* calibration files directory */
	CLA->FstatsFile2=optarg;
	break;
      case '3':
	/* SFT directory */
	CLA->FstatsFile3=optarg;
	break;
      case '4':
	/* calibration files directory */
	CLA->FstatsFile4=optarg;
	break;
      case 'o':
	/* calibration files directory */
	CLA->OutputFile=optarg;
	break;
      case 'f':
	/* Spin down order */
	CLA->Deltaf=atof(optarg);
	break;
      case 'a':
	/* Spin down order */
	CLA->DeltaAlpha=atof(optarg);
	break;
      case 's':
	/* Spin down order */
	CLA->fmin=atof(optarg);
	break;
      case 'e':
	/* Spin down order */
	CLA->fmax=atof(optarg);
	break;
      case 'd':
	/* Spin down order */
	CLA->DeltaDelta=atof(optarg);
	break;
      case 'h':
	/* print usage/help message */
	fprintf(stderr,"Arguments are (defaults):\n");
	fprintf(stderr,"\t--fstatsfile1 (-1)\tSTRING\tFirst candidates Fstats file\n");
	fprintf(stderr,"\t--fstatsfile1 (-1)\tSTRING\tFirst candidates Fstats file\n");
	fprintf(stderr,"\t--fstatsfile3 (-3)\tSTRING\tFstats used to compute false alarm for -1\n");
	fprintf(stderr,"\t--fstatsfile4 (-4)\tSTRING\tFstats used to compute false alarm for -2\n");
	fprintf(stderr,"\t--outputfile  (-o)\tSTRING\tName of ouput candidates file\n");
	fprintf(stderr,"\t--frequency-window (-f)\tFLOAT\tFrequency window in Hz (0.0)\n");
	fprintf(stderr,"\t--alpha-window (-a)\tFLOAT\tAlpha window in radians (0.0)\n");
	fprintf(stderr,"\t--delta-window (-d)\tFLOAT\tDelta window in radians (0.0)\n");
	fprintf(stderr,"\t--fmin (-s)\tFLOAT\t Minimum frequency of candidate in 1st IFO\n");
	fprintf(stderr,"\t--fmax (-e)\tFLOAT\t Maximum frequency of candidate in 1st IFO\n");
	fprintf(stderr,"\t--help        (-h)\tFLOAT\tThis message\n");
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

  if(CLA->FstatsFile1 == NULL)
    {
      fprintf(stderr,"No 1st candidates file specified; input with -1 option.\n");
      fprintf(stderr,"For help type ./polka -h \n");
      return 1;
    }      
  if(CLA->FstatsFile2 == NULL)
    {
      fprintf(stderr,"No 2nd candidates file specified; input with -2 option.\n");
      fprintf(stderr,"For help type ./polka -h \n");
      return 1;
    }      
  if(CLA->OutputFile == NULL)
    {
      fprintf(stderr,"No ouput filename specified; input with -o option.\n");
      fprintf(stderr,"For help type ./polka -h \n");
      return 1;
    }      

  if(CLA->fmin == 0.0)
    {
      fprintf(stderr,"No minimum frequency specified.\n");
      fprintf(stderr,"For help type ./polka -h \n");
      return 1;
    }      

  if(CLA->fmax == 0.0)
    {
      fprintf(stderr,"No maximum frequency specified.\n");
      fprintf(stderr,"For help type ./polka -h \n");
      return 1;
    }      

  return errflg;
}


/* Explanation of locate

the function locate returns the lower value of the two indices of the array xx
that bound the value given x.

consider xx to be in descending order:

fmax                                      fmin
   |     |     |     |     |     |     |     |
j: 0     1     2     3     4     5     6     7

            ^           ^            ^
	 f2+Df	        f2        f2-Df

locate will return if2max=1 and if2min=5

*/
