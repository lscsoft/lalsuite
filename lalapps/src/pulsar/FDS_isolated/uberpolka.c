/*********************************************************************************/
/*       uberpolka - the pulsar koinzidenz analysis code for einstein@home       */
/*                                                                               */
/*                                     X. Siemens                                */
/*                   (takes in two Fstats file to look for coincidence)          */
/*                    modified by Bernd Machenschalk for Einstein@Home           */
/*                                                                               */
/*                                  UWM - January  2005                          */
/*********************************************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/LALStatusMacros.h>
#include <lal/ConfigFile.h>

#include <lalapps.h>

#include "getopt.h"

RCSID ("$Id$");

/* some error codes and messages */
#define POLKAC_ENULL            1
#define POLKAC_ENONULL          2
#define POLKAC_ESYS             3
#define POLKAC_EINVALIDFSTATS   4
#define POLKAC_EMEM             5

#define POLKAC_MSGENULL         "Arguments contained an unexpected null pointer"
#define POLKAC_MSGENONULL       "Input pointer was not NULL"
#define POLKAC_MSGESYS          "System call failed (probably file IO"
#define POLKAC_MSGEINVALIDFSTATS "Invalid Fstats file"
#define POLKAC_MSGEMEM          "Sorry, ran out of memory... bye."

#ifndef USE_BOINC
#define USE_BOINC 0
#endif

#if USE_BOINC
/* BOINC includes */
#include "boinc_api.h"
#include "filesys.h"
/* alias fopen - this will take care of architecture-specific problem handling */
#define fopen boinc_fopen
/* this global variable communicates the output filename to the calling routine in CFS */
extern CHAR *Outputfilename;
#endif

struct PolkaCommandLineArgsTag 
{
  char *FstatsFile1; /* Names of Fstat files to be read in */
  char *FstatsFile2;
  char *OutputFile;
  REAL8 Deltaf;      /* Size of coincidence window in Hz */
  REAL8 DeltaAlpha;  /* Size of coincidence window in radians */
  REAL8 DeltaDelta;  /* Size of coincidence window in radians */
  REAL8 fmin;        /* Minimum frequency of candidate in first IFO */
  REAL8 fmax;        /* Maximum frequency of candidate in first IFO */
  UINT4 EAH;         /* Einstein at home flag for alternative output */ 
} PolkaCommandLineArgs;

/* This structure contains the indices corresponding to the 
coarse frequency and sky bins */
typedef struct CandINDICESTag
{
  REAL8 f;        /* Frequency */
  REAL8 Alpha;    /* longitude */
  REAL8 Delta;    /* latitude */
  REAL8 F;        /* Maximum value of F for the cluster */
  INT4 iFreq;
  INT4 iDelta;
  INT4 iCand;
  INT4 iCandSorted;
} CandINDICES;

typedef struct CandidateTag 
{
  REAL8 *f;        /* Frequency */
  REAL8 *Alpha;    /* longitude */
  REAL8 *Delta;    /* latitude */
  REAL8 *F;        /* Maximum value of F for the cluster */
  REAL8 *fa;       /* false alarm probability for that candidate */
  UINT4 *Ctag;     /* tag for candidate if it's been found in coincidence */
  INT4  *CtagCounter;     /* contains the cumulative sum of coincident candidates so far */
  UINT4 length;    /* number of candidates in list */
} CandidateList;


typedef struct CoincidentPairsTag 
{
  UINT4 c1;             /* number in Fstats file that corresponds to first member of pair */
  UINT4 c2;             /* number in Fstats file that corresponds to second member of pair */
  REAL8 fa;             /* joint false alarm for that pair */
} CoincidentPairs;

int ReadCommandLine(int argc,char *argv[],struct PolkaCommandLineArgsTag *CLA);
int ReadCandidateFiles(struct PolkaCommandLineArgsTag CLA);
int ReadOneCandidateFile(CandidateList *CList, const char *fname);
int compareCIStructs(const void *ip, const void *jp);
int compareCPfa(const void *ip, const void *jp);
int FineCoincidenceTest(CandINDICES c1, CandINDICES c2, struct PolkaCommandLineArgsTag CLA);
int OutputCoincidences(struct PolkaCommandLineArgsTag CLA);

extern INT4 lalDebugLevel;

CandidateList CList1, CList2; /* treat up to 4 candidate files */
CandINDICES *SortedC1,*SortedC2;
CoincidentPairs *CP;

INT4 numCoincidences=0;
REAL8 MaxAngularDistance;

#ifndef FALSE
#define FALSE (1==0)
#endif
#ifndef TRUE
#define TRUE  (1==1)
#endif

/* main() mapped to polka() if using boinc */
#if USE_BOINC
int polka(int argc,char *argv[])
#else
int main(int argc,char *argv[]) 
#endif
{
  UINT4 i;
  lalDebugLevel = 0;

  /* Reads command line arguments */
  if (ReadCommandLine(argc,argv,&PolkaCommandLineArgs)) return 1;

  MaxAngularDistance=sqrt(pow(PolkaCommandLineArgs.DeltaAlpha,2)+pow(PolkaCommandLineArgs.DeltaDelta,2))+1e-8;

  /* Reads in candidare files */
  if (ReadCandidateFiles(PolkaCommandLineArgs)) return 2;

  if (CList1.length != 0 && CList2.length != 0 )
    {
      /* create arrays of candidates to be sorted */
      if (!(SortedC1=(CandINDICES *)LALMalloc(sizeof(CandINDICES) * CList1.length))){
        fprintf(stderr,"Unable to allocate index1 array in main\n");
        return 1;
      }
      if (!(SortedC2=(CandINDICES *)LALMalloc(sizeof(CandINDICES) * CList2.length))){
        fprintf(stderr,"Unable to allocate index2 array in main\n");
        return 1;
      }

      /* Initialise arrays of sorted candidates to use for bsearch */
      for (i=0;i<CList1.length;i++) 
        {
          SortedC1[i].f=CList1.f[i];
          SortedC1[i].Delta=CList1.Delta[i];
          SortedC1[i].Alpha=CList1.Alpha[i];
          SortedC1[i].F=CList1.F[i];
          SortedC1[i].iFreq=(INT4) (CList1.f[i]/(PolkaCommandLineArgs.Deltaf));
          SortedC1[i].iDelta=(INT4)(CList1.Delta[i]/(PolkaCommandLineArgs.DeltaDelta));
	  SortedC1[i].iCand=i;
        }
      for (i=0;i<CList2.length;i++) 
        {
          SortedC2[i].f=CList2.f[i];
          SortedC2[i].Delta=CList2.Delta[i];
          SortedC2[i].Alpha=CList2.Alpha[i];
          SortedC2[i].F=CList2.F[i];
          SortedC2[i].iFreq=(INT4) (CList2.f[i]/(PolkaCommandLineArgs.Deltaf));
          SortedC2[i].iDelta=(INT4)(CList2.Delta[i]/(PolkaCommandLineArgs.DeltaDelta));
	  SortedC2[i].iCand=i;
        }
      
      /* sort arrays of candidates */
      qsort(SortedC1, (size_t)CList1.length, sizeof(CandINDICES), compareCIStructs);
      qsort(SortedC2, (size_t)CList2.length, sizeof(CandINDICES), compareCIStructs);

      for (i=0;i<CList1.length;i++) SortedC1[i].iCandSorted=i;
      for (i=0;i<CList2.length;i++) SortedC2[i].iCandSorted=i;
      
      /* loop iover candidates in first array */
      for (i=0; i < CList1.length; i++)
        {
#if USE_BOINC
          /* make sure the cpu time is updated */ 
          if (boinc_time_to_checkpoint())
            boinc_checkpoint_completed();
#endif
          if  (SortedC1[i].f >= PolkaCommandLineArgs.fmin && SortedC1[i].f <= PolkaCommandLineArgs.fmax)
            {
              int iFreq2, iDelta2;

              /* Loop to run bsearch etc. on all surrounding boxes */
              for (iFreq2=SortedC1[i].iFreq-1; iFreq2 <= SortedC1[i].iFreq+1; iFreq2++)
                {
		  for (iDelta2=SortedC1[i].iDelta-1; iDelta2 <= SortedC1[i].iDelta+1; iDelta2++)
		    {
		      CandINDICES *p, can;
		      
		      can.iFreq=iFreq2;
		      can.iDelta=iDelta2;
                          
		      p=bsearch(&can,SortedC2,(size_t)CList2.length, sizeof(CandINDICES),compareCIStructs);
		      
		      if (p != NULL)
			{
			  /* Now we've found at least one candidate */
			  /* we need to move to the right edge (without segfaulting!) */
			  
			  while ( p->iCandSorted > 0 && !compareCIStructs(p, p-1) )
			    p--;
			  
                          /* Now p points to first coincident event in the second list */
			  
                          /* Now loop over candidates found in the second list and do the fine coincidence test */
                          if(FineCoincidenceTest(SortedC1[i],*p, PolkaCommandLineArgs)) return 3;
                          while ( p->iCandSorted <  (int)CList2.length-1 &&  !compareCIStructs(p, p+1) )
                            { 
                              p++;
                              if(FineCoincidenceTest(SortedC1[i],*p, PolkaCommandLineArgs)) return 3;
                            }

			}/* check that besearch was non-null */
		    } /* loop over deltas */
		}/* loop over frequencies */    
	    } /* check that frequency lies between two input bounds */
	}/* loop over 1st candidate list */
      
      LALFree(SortedC1);
      LALFree(SortedC2);
      
    }

  /* freeing a CList is a bit tedious, so we use a macro */
#define freeCList(x) do { LALFree((x).f); LALFree((x).Alpha); LALFree((x).Delta); LALFree((x).F); LALFree((x).fa); LALFree((x).Ctag);LALFree((x).CtagCounter);} while(0)
  
  
  /* Ouput candidates */
  if (OutputCoincidences( PolkaCommandLineArgs )) return 4;
  
  if(CList1.length) freeCList(CList1);
  if(CList2.length) freeCList(CList2);

  LALCheckMemoryLeaks(); 

  return 0;
 
}

/*******************************************************************************/

int OutputCoincidences(struct PolkaCommandLineArgsTag CLA)
{
  INT4 *indicesCCfa=NULL;
  INT4 i;
  FILE *fpOut;
#if USE_BOINC
  static char resolved_filename[256];
#endif

  /* allocate space */
  if (numCoincidences != 0){ 
    if (!(indicesCCfa=(INT4 *)LALMalloc(sizeof(INT4) * numCoincidences))){
      fprintf(stderr,"Unable to allocate index array in main\n");
      return 1;
    }
  }

  for (i=0; i < numCoincidences; i++) indicesCCfa[i]=i;
  qsort((void *)indicesCCfa, (size_t)numCoincidences, sizeof(int), compareCPfa);

  /* open and write the file */
#if USE_BOINC
  if (boinc_resolve_filename(CLA.OutputFile, resolved_filename, sizeof(resolved_filename))) {
    fprintf(stderr,
            "Can't resolve file \"%s\"\n"
            "If running a non-BOINC test, create [INPUT] or touch [OUTPUT] file\n",
            CLA.OutputFile);
    boinc_finish(2);
  }
  fpOut=fopen(resolved_filename,"w");
#else
  fpOut=fopen(CLA.OutputFile,"w");       
#endif
  if (!CLA.EAH)
    {
      /* sort in increasing probability of joint false alarm */
      for (i=0; i < numCoincidences; i++) 
        {
          UINT4 k = indicesCCfa[i];  /* print out ordered by joint significance */
          UINT4 k1 = CP[k].c1;
          UINT4 k2 = CP[k].c2;
          fprintf(fpOut,"%1.15le %le %le %le %le %1.15le %le %le %le %le %le\n",
                  CList1.f[k1], CList1.Alpha[k1], CList1.Delta[k1], 
                  CList1.F[k1], CList1.fa[k1],
                  CList2.f[k2], CList2.Alpha[k2], CList2.Delta[k2], 
                  CList2.F[k2], CList2.fa[k2],CList1.fa[k1]*CList2.fa[k2]);
        }
    }else{
      fprintf(fpOut,"%%1\n");
      {
        int k=-1;    
        for (i=0; i < (int)CList1.length; i++)
          {
            if (CList1.Ctag[i]) 
              {
                k++;
                fprintf(fpOut,"%16.12f %10.8f %10.8f %20.17f\n",
                        CList1.f[i],CList1.Alpha[i],CList1.Delta[i],CList1.F[i]);
                CList1.CtagCounter[i]=k;
              }
          }
      }
      fprintf(fpOut,"%%2\n");
      {
        int k=-1;
        for (i=0; i < (int)CList2.length; i++)
          {
            if (CList2.Ctag[i]) 
              {
                k++;
                fprintf(fpOut,"%16.12f %10.8f %10.8f %20.17f\n",
                        CList2.f[i],CList2.Alpha[i],CList2.Delta[i],CList2.F[i]);  
                CList2.CtagCounter[i]=k;
              }    
          }
      }
      fprintf(fpOut,"%%coincidences\n");
      /* sort in increasing probability of joint false alarm */
      qsort((void *)indicesCCfa, (size_t)numCoincidences, sizeof(int), compareCPfa);
      
      for (i=0; i < numCoincidences; i++) 
        {
          UINT4 k = indicesCCfa[i];  /* print out ordered by joint significance */
          fprintf(fpOut,"%d %d %le\n",
                  CList1.CtagCounter[CP[k].c1],CList2.CtagCounter[CP[k].c2],CP[k].fa);
        }
    }
  /* write end marker */
  fprintf(fpOut,"%%DONE\n");    
#if USE_BOINC
  /* make the output filename known to teh boinc main() routine in ComputeFStatistic */
  Outputfilename=resolved_filename;
#endif
  fclose(fpOut);

  if (numCoincidences != 0){ 
    LALFree ( CP );
    LALFree(indicesCCfa);
  }

  return 0;
}

/*******************************************************************************/
#define BLOCKSIZE 16384

int FineCoincidenceTest(CandINDICES c1, CandINDICES c2, struct PolkaCommandLineArgsTag CLA)
{
  
  REAL8 f1=c1.f,f2=c2.f, F1=c1.F, F2=c2.F;
  REAL8 Alpha1=c1.Alpha,Delta1=c1.Delta;
  REAL8 Alpha2=c2.Alpha,Delta2=c2.Delta;
  REAL8 n1[3],n2[3],AngularDistance;
  REAL8 cosAngularDistance, difff;

  n1[0]=cos(Alpha1)*cos(Delta1);
  n1[1]=sin(Alpha1)*cos(Delta1);
  n1[2]=sin(Delta1);
  
  n2[0]=cos(Alpha2)*cos(Delta2);
  n2[1]=sin(Alpha2)*cos(Delta2);
  n2[2]=sin(Delta2);
 
  cosAngularDistance=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
  if (cosAngularDistance  >  1.0) cosAngularDistance =  1.0;
  if (cosAngularDistance  < -1.0) cosAngularDistance = -1.0;
        
  AngularDistance=acos(cosAngularDistance);
  
  difff=fabs(f1 - f2);
              
  /* check difference in frequencies because we're not guaranteed 
     sufficient closeness at the edges of array */
  if ( difff <= CLA.Deltaf) 
    {
      if ( AngularDistance <= MaxAngularDistance )
        {       
          CoincidentPairs *thisCP;

          /* tag the candidates that have been found in coincidence */
          CList1.Ctag[c1.iCand]=1;
          CList2.Ctag[c2.iCand]=1;
                  
          /* seems we found a coincident candidate: let's make space for it to be stored */
          
          if (numCoincidences % BLOCKSIZE == 0)
            {
              if ( (CP = LALRealloc ( CP, (numCoincidences + BLOCKSIZE) * sizeof(CoincidentPairs) )) == NULL) {
                fprintf (stderr, "Error: polka ran out of memory allocating %d coincident candidate (CP).\n",numCoincidences);
                return 1;
              }
            }

          numCoincidences ++;
          
          thisCP = &(CP[ numCoincidences - 1]); /* point to current new coincidences */
          
          CList1.fa[c1.iCand]=(1+F1/2)*exp(-F1/2);
          CList2.fa[c2.iCand]=(1+F2/2)*exp(-F2/2);
          
          thisCP->c1=c1.iCand;
          thisCP->c2=c2.iCand;
          thisCP->fa=CList1.fa[c1.iCand]*CList2.fa[c2.iCand];

        }
    }
  return 0;
} 

/*******************************************************************************/

/* Sorting function to sort candidate indices INCREASING order of f, delta, alpha */
int compareCIStructs(const void *a, const void *b)
{
  const CandINDICES *ip = a;
  const CandINDICES *jp = b;
  INT4 ifreq1,ifreq2;

  ifreq1=ip->iFreq;
  ifreq2=jp->iFreq;

  if (ifreq1 < ifreq2)
    return -1;
  
  if (ifreq1 == ifreq2)
    {
      INT4 iDelta1, iDelta2;

      iDelta1=ip->iDelta;
      iDelta2=jp->iDelta;
      
      if (iDelta1 < iDelta2)
        return -1;

      if (iDelta1 > iDelta2)
        return 1;
  
      if (iDelta1 == iDelta2)           
	return 0;
    }
  return 1;
}

/*******************************************************************************/

/* Sorting function to sort second candidate list into increasing order of fa */
int compareCPfa(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=CP[*(const int *)ip].fa;
  dj=CP[*(const int *)jp].fa;

  if (di<dj)
    return -1;
  
  if (di==dj)
    return (ip > jp);

  return 1;
}

/*******************************************************************************/

int ReadCandidateFiles(struct PolkaCommandLineArgsTag CLA)
{

  if (ReadOneCandidateFile ( &CList1, CLA.FstatsFile1)) return 1;

  if (ReadOneCandidateFile ( &CList2, CLA.FstatsFile2)) return 1;

  return 0;

} /* ReadCandidateFiles() */


/*******************************************************************************/

#define DONE_MARKER "%DONE"

/* read and parse the given candidate 'Fstats'-file fname into the candidate-list CList */
int  ReadOneCandidateFile (CandidateList *CList, const char *fname)
{
  UINT4 i;
  UINT4 numlines;
  REAL8 dmp;
  char line1[256],line2[5];
  CandidateList cands;
  FILE *fp;
  INT4 read;

  /* ------ Open and count candidates file ------ */
  i=0;
  fp=fopen(fname,"r");
  if (fp==NULL) 
    {
      fprintf(stderr,"File %s doesn't exist!\n",fname);
      return 1;
    }
  while(fgets(line1,sizeof(line1),fp))
    {
      strncpy(line2,line1,5);
      i++;
    }
  numlines=i;
  fclose(fp);     
  /* -- close candidate file -- */

  if ( numlines == 0) 
    {
      LALPrintError ("ERROR: File '%s' is empty and is not properly terminated by '%s' marker!\n\n", fname, DONE_MARKER);
      return 1;
    }

  /* check validity of this Fstats-file */
  if ( strcmp(line2, DONE_MARKER ) ) 
    {
      LALPrintError ("ERROR: File '%s' is not properly terminated by '%s' marker!\n\n", fname, DONE_MARKER);
      return 1;
    }
  else
    numlines --;        /* avoid stepping on DONE-marker */
  
  /* reserve memory for fstats-file contents */
  if (numlines > 0)
    {
      cands.f     = LALCalloc (numlines, sizeof(REAL8));
      cands.Alpha = LALCalloc (numlines, sizeof(REAL8));
      cands.Delta = LALCalloc (numlines, sizeof(REAL8));
      cands.F     = LALCalloc (numlines, sizeof(REAL8));
      cands.fa    = LALCalloc (numlines, sizeof(REAL8));
      cands.Ctag  = LALCalloc (numlines, sizeof(UINT4));
      cands.CtagCounter  = LALCalloc (numlines, sizeof(INT4));

      if ( !cands.f || !cands.Alpha || !cands.Delta || !cands.F || !cands.fa || !cands.Ctag || !cands.CtagCounter )
        {
          LALPrintError ("Could not allocate memory for candidate file %s\n\n", fname);
          return 1;
        }
    }

  /* ------ Open and count candidates file ------ */
  i=0;
  fp=fopen(fname,"r");
  if (fp==NULL) 
    {
      fprintf(stderr,"File %s doesn't exist!\n",fname);
      return 1;
    }
  while(fgets(line1,sizeof(line1),fp) && i < numlines )
    {

      cands.Ctag[i]=0;
      cands.CtagCounter[i]=-1;

      read = sscanf (line1, 
                     "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
                     " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT, 
                     &(cands.f[i]), &(cands.Alpha[i]), &(cands.Delta[i]), &dmp, &dmp, &dmp, &(cands.F[i]) );

      if ( read != 7 )
        {
          LALPrintError ("Failed to parse line %d in file '%s' \n", i+1, fname);
          LALFree (cands.f);
          LALFree (cands.Alpha);
          LALFree (cands.Delta);
          LALFree (cands.F);
          LALFree (cands.fa);
          LALFree (cands.Ctag);
          LALFree (cands.CtagCounter);
          return 1;
        }

      i++;
    }
  fclose(fp);     
  /* -- close 1st candidate file -- */
  
  /* return final candidate-list */
  CList->length = numlines;
  CList->f      = cands.f;
  CList->Alpha  = cands.Alpha;
  CList->Delta  = cands.Delta;
  CList->F      = cands.F;
  CList->fa     = cands.fa;
  CList->Ctag   = cands.Ctag;
  CList->CtagCounter   = cands.CtagCounter;

  return 0;

} /* ReadOneCandidateFile() */

/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct PolkaCommandLineArgsTag *CLA) 
{
  INT2 errflg = 0;
  INT4 c; 
  INT4 option_index = 0;

  const char *optstring = "h1:2:f:a:d:m:M:o:s:e:b";
  struct option long_options[] =
    {
      {"fstatsfile1",           required_argument, 0,   '1'},
      {"fstatsfile2",           required_argument, 0,   '2'},
      {"frequency-window",      required_argument, 0,   'f'},
      {"delta-window",          required_argument, 0,   'd'},
      {"alpha-window",          required_argument, 0,   'a'},
      {"fmin",                  required_argument, 0,   's'},
      {"fmax",                  required_argument, 0,   'e'},
      {"outputfile",            required_argument, 0,   'o'},
      {"EAHoutput",             no_argument, 0,         'b'},
      {"help",                  no_argument, 0,         'h'},
      {0, 0, 0, 0}
    };

  /* Initialize default values */
  CLA->FstatsFile1=NULL;
  CLA->FstatsFile2=NULL;
  CLA->OutputFile=NULL;
  CLA->Deltaf=0.0;
  CLA->DeltaAlpha=0;
  CLA->DeltaDelta=0;
  CLA->fmin=0;
  CLA->fmax=0;
  CLA->EAH=0;

  /* reset gnu getopt */
  optind = 0;

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
      case 'b':
        /* Spin down order */
        CLA->EAH=1;
        break;
      case 'h':
        /* print usage/help message */
        fprintf(stderr,"Arguments are (defaults):\n");
        fprintf(stderr,"\t--fstatsfile1 (-1)\tSTRING\tFirst candidates Fstats file\n");
        fprintf(stderr,"\t--fstatsfile2 (-2)\tSTRING\tSecond candidates Fstats file\n");
        fprintf(stderr,"\t--outputfile  (-o)\tSTRING\tName of ouput candidates file\n");
        fprintf(stderr,"\t--frequency-window (-f)\tFLOAT\tFrequency window in Hz (0.0)\n");
        fprintf(stderr,"\t--alpha-window (-a)\tFLOAT\tAlpha window in radians (0.0)\n");
        fprintf(stderr,"\t--delta-window (-d)\tFLOAT\tDelta window in radians (0.0)\n");
        fprintf(stderr,"\t--fmin (-s)\tFLOAT\t Minimum frequency of candidate in 1st IFO\n");
        fprintf(stderr,"\t--fmax (-e)\tFLOAT\t Maximum frequency of candidate in 1st IFO\n");
        fprintf(stderr,"\t--EAHoutput (-b)\tFLAG\t Einstein at home output flag. \n");
        fprintf(stderr,"\t--help        (-h)\t\tThis message\n");
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

