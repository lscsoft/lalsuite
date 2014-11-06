/*
*  Copyright (C) 2007 Bruce Allen, Bernd Machenschalk, Xavier Siemens
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
 * \author Xavier Siemens,  Bruce Allen,  Bernd Machenschalk
 * \brief takes in two Fstats file to look for coincidence
 */

/*********************************************************************************/
/*       uberpolka - the pulsar koinzidenz analysis code for einstein@home       */
/*                                                                               */
/*                   Xavier Siemens,  Bruce Allen,  Bernd Machenschalk           */
/*                                                                               */
/*                   (takes in two Fstats file to look for coincidence)          */
/*                                                                               */
/*                                  UWM - January  2005                          */
/*********************************************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#if HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <math.h>

#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/LALStatusMacros.h>
#include <lal/ConfigFile.h>

#include <lalapps.h>

#include "getopt.h"

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
/* communicating the progress to the graphics thread */
extern double *fraction_done_hook;
/* define XLALPrintError locally again, otherwise the stderr redirection
   doesn't seem to work on the Mac */
int myPrintError( const char *fmt, ... )
{
  int n;
  va_list ap;
  va_start( ap, fmt );
  n = vfprintf( stderr, fmt, ap );
  va_end( ap );
  return n;
}
#else
#define myPrintError XLALPrintError
#endif

/* this is defined in C99 and *should* be in math.h.  Long term
   protect this with a HAVE_FINITE */
#ifdef _MSC_VER
#include <float.h>
#define finite _finite
#endif

#define UBERPOLKA_EXIT_ERRCLINE 31
#define UBERPOLKA_EXIT_READCND  32
#define UBERPOLKA_EXIT_FCTEST   33
#define UBERPOLKA_EXIT_OUTFAIL  34

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
typedef struct CandidateListTag
{
  REAL8 f;           /* Frequency */
  REAL4 Alpha;       /* longitude -> REAL4*/
  REAL4 Delta;       /* latitude -> REAL4 */
  REAL4 F;           /* Maximum value of F for the cluster -> REAL4*/
  REAL4 lfa;         /* log of false alarm probability for that candidate ->REAL4*/
  INT4 iCand;
  INT4 CtagCounter; /* contains the cumulative sum of coincident candidates so far */
  INT4  iFreq;       /* INT2 , delete? */
  INT2 iDelta;      /* -157-157 -> INT2, delete?  */
  UINT4  Ctag;        /* tag for candidate if it's been found in coincidence, just a Bit, maybe coded as CtagC-2 */
} CandidateList; /* ~ Fstat lines */ 

typedef struct CoincidentPairsTag 
{
  INT4 c1;             /* number in Fstats file that corresponds to first member of pair */
  INT4 c2;             /* number in Fstats file that corresponds to second member of pair */
  /* REAL8 fa; */       /* joint false alarm for that pair */
  REAL4 lfa;            /* log of joint false alarm for that pair */
} CoincidentPairs; /* ~ coninc */

int ReadCommandLine(int argc,char *argv[],struct PolkaCommandLineArgsTag *CLA);
int ReadCandidateFiles(struct PolkaCommandLineArgsTag CLA);
int ReadOneCandidateFile(CandidateList **CList, const char *fname);
int compareCIStructs(const void *ip, const void *jp);
int compareCdaf(const void *ip, const void *jp);
int compareCPfa(const void *ip, const void *jp);
int FineCoincidenceTest(CandidateList c1, CandidateList c2, struct PolkaCommandLineArgsTag CLA);
int OutputCoincidences(struct PolkaCommandLineArgsTag CLA);


CandidateList *SortedC1,*SortedC2;
CoincidentPairs *CP;

UINT4  CLength=0,CLength1=0,CLength2=0;

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
#if USE_BOINC
  REAL8 local_fraction_done;
#endif


  /* Reads command line arguments */
  if (ReadCommandLine(argc,argv,&PolkaCommandLineArgs)) {
    fprintf(stderr,"ReadCommandLine failed\n");
    return UBERPOLKA_EXIT_ERRCLINE;
  }

  MaxAngularDistance=sqrt(pow(PolkaCommandLineArgs.DeltaAlpha,2)+pow(PolkaCommandLineArgs.DeltaDelta,2))+1e-8;

  /* Reads in candidare files, set CLength1 and CLength2 */
  if (ReadCandidateFiles(PolkaCommandLineArgs)) {
    fprintf(stderr,"ReadCandidateFiles failed\n");
    return UBERPOLKA_EXIT_READCND;
  }

  if (CLength1 != 0 && CLength2 != 0 )
    {
      /* Initialise arrays of sorted candidates to use for bsearch */
      for (i=0;i<CLength1;i++)
        {
          SortedC1[i].iFreq=(INT4) (SortedC1[i].f/(PolkaCommandLineArgs.Deltaf));
          SortedC1[i].iDelta=(INT4)(SortedC1[i].Delta/(PolkaCommandLineArgs.DeltaDelta));
         }
      for (i=0;i<CLength2;i++)
        {
          SortedC2[i].iFreq=(INT4) (SortedC2[i].f/(PolkaCommandLineArgs.Deltaf));
          SortedC2[i].iDelta=(INT4)(SortedC2[i].Delta/(PolkaCommandLineArgs.DeltaDelta));
        }

      /* sort arrays of candidates */
      qsort(SortedC1, (size_t)CLength1, sizeof(CandidateList), compareCIStructs);
      qsort(SortedC2, (size_t)CLength2, sizeof(CandidateList), compareCIStructs);

      for (i=0;i<CLength1;i++) SortedC1[i].iCand=i;
      for (i=0;i<CLength2;i++) SortedC2[i].iCand=i;
      
      /* loop iover candidates in first array */
      for (i=0; i < CLength1; i++)
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
                      CandidateList *p, can;
                      
                      can.iFreq=iFreq2;
                      can.iDelta=iDelta2;
                          
                      p=bsearch(&can,SortedC2,(size_t)CLength2, sizeof(CandidateList),compareCIStructs);
                      
                      if (p != NULL)
                        {
                          /* Now we've found at least one candidate */
                          /* we need to move to the right edge (without segfaulting!) */
                          
                          while ( p->iCand > 0 && !compareCIStructs(p, p-1) )
                            p--;
                          
                          /* Now p points to first coincident event in the second list */
                          
                          /* Now loop over candidates found in the second list and do the fine coincidence test */
                          if(FineCoincidenceTest(SortedC1[i],*p, PolkaCommandLineArgs)) {
                            fprintf(stderr,"FineCoincidenceTest failed\n");
                            return UBERPOLKA_EXIT_FCTEST;
                          }
                          while ( (int)p->iCand <  (int)CLength2-1 &&  !compareCIStructs(p, p+1) )
                            { 
                              p++;
                              if(FineCoincidenceTest(SortedC1[i],*p, PolkaCommandLineArgs)) {
                                fprintf(stderr,"FineCoincidenceTest failed in loop\n");
                                return UBERPOLKA_EXIT_FCTEST;
                              }
                            }

                        }/* check that besearch was non-null */
                    } /* loop over deltas */
                }/* loop over frequencies */    
            } /* check that frequency lies between two input bounds */
#if USE_BOINC
          local_fraction_done = 0.99 + 0.01 * (double)i / (double)CLength1;
          /* update progress, the last % is reserved for polka */
          boinc_fraction_done(local_fraction_done);
          /* pass variable externally to graphics routines */
          if (fraction_done_hook != NULL)
            *fraction_done_hook = local_fraction_done;
#endif
        }/* loop over 1st candidate list */      
    }/* check that we have candidates in both files */
  
  /* Ouput candidates */
  if (OutputCoincidences( PolkaCommandLineArgs )) {
    fprintf(stderr,"OutputCoincidences failed\n");
    return UBERPOLKA_EXIT_OUTFAIL;
  }
  
  if(CLength1) LALFree(SortedC1);
  if(CLength2) LALFree(SortedC2);

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
      myPrintError("Unable to allocate index array in main\n");
      return 1;
    }
  }

  for (i=0; i < numCoincidences; i++) indicesCCfa[i]=i;
  qsort((void *)indicesCCfa, (size_t)numCoincidences, sizeof(int), compareCPfa);

  /* open and write the file */
#if USE_BOINC
  if (boinc_resolve_filename(CLA.OutputFile, resolved_filename, sizeof(resolved_filename))) {
    myPrintError(
            "Can't resolve file \"%s\"\n"
            "If running a non-BOINC test, create [INPUT] or touch [OUTPUT] file\n",
            CLA.OutputFile);
    boinc_finish(2);
  }
  fpOut=fopen(resolved_filename,"wb");
  if (!fpOut){
    myPrintError("Unable to open output file \"%s\" for writing\n", resolved_filename);
    return 1;
  }
#else
  fpOut=fopen(CLA.OutputFile,"wb");
  if (!fpOut){
    myPrintError("Unable to open output file \"%s\" for writing\n", CLA.OutputFile);
    return 1;
  }
#endif
  if (!CLA.EAH)
    {
      /* sort in increasing probability of joint false alarm */
      for (i=0; i < numCoincidences; i++) 
        {
          UINT4 k = indicesCCfa[i];  /* print out ordered by joint significance */
          UINT4 k1 = CP[k].c1;
          UINT4 k2 = CP[k].c2;
          fprintf(fpOut,"%1.15e %e %e %e %e %1.15e %e %e %e %e %e\n",
                  SortedC1[k1].f, SortedC1[k1].Alpha, SortedC1[k1].Delta, 
                  SortedC1[k1].F, exp(SortedC1[k1].lfa),
                  SortedC2[k2].f, SortedC2[k2].Alpha, SortedC2[k2].Delta, 
                  SortedC2[k2].F, exp(SortedC2[k2].lfa),
                                  exp(SortedC1[k1].lfa)*exp(SortedC2[k2].lfa));
        }
    }else{
      /* sort by delta, alpha, f (just like Fstats file) */
      qsort(SortedC1, (size_t)CLength1, sizeof(CandidateList), compareCdaf);
      qsort(SortedC2, (size_t)CLength2, sizeof(CandidateList), compareCdaf);

      fprintf(fpOut,"%%1\n");
      {
        int k=-1;    
        for (i=0; i < (int)CLength1; i++)
          {
            if (SortedC1[i].Ctag) 
              {
                k++;
                fprintf(fpOut,"%16.12f %10.8f %10.8f %20.17f\n",
                        SortedC1[i].f,SortedC1[i].Alpha,SortedC1[i].Delta,SortedC1[i].F);
                SortedC1[i].CtagCounter=k;
              }
          }
      }
      fprintf(fpOut,"%%2\n");
      {
        int k=-1;
        for (i=0; i < (int)CLength2; i++)
          {
            if (SortedC2[i].Ctag) 
              {
                k++;
                fprintf(fpOut,"%16.12f %10.8f %10.8f %20.17f\n",
                        SortedC2[i].f,SortedC2[i].Alpha,SortedC2[i].Delta,SortedC2[i].F);  
                SortedC2[i].CtagCounter=k;
              }    
          }
      }
      /* sort arrays of candidates back to the order they were in */
      qsort(SortedC1, (size_t)CLength1, sizeof(CandidateList), compareCIStructs);
      qsort(SortedC2, (size_t)CLength2, sizeof(CandidateList), compareCIStructs);

      fprintf(fpOut,"%%coincidences\n");
      /* sort in increasing probability of joint false alarm */
      qsort((void *)indicesCCfa, (size_t)numCoincidences, sizeof(int), compareCPfa);
      
      for (i=0; i < numCoincidences; i++) 
        {
          UINT4 k = indicesCCfa[i];  /* print out ordered by joint significance */
          fprintf(fpOut,"%d %d %e\n",
                  SortedC1[CP[k].c1].CtagCounter,SortedC2[CP[k].c2].CtagCounter,exp(CP[k].lfa));
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

int FineCoincidenceTest(CandidateList c1, CandidateList c2, struct PolkaCommandLineArgsTag CLA)
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
          SortedC1[c1.iCand].Ctag=1;
          SortedC2[c2.iCand].Ctag=1;
                  
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
          /*
          SortedC1[c1.iCand].fa=(1+F1/2)*exp(-F1/2);
          SortedC2[c2.iCand].fa=(1+F2/2)*exp(-F2/2);

          log(a*exp(b)) == log(exp(log(a))*exp(b)) == log(exp(log(a)+b)) == log(a)+b
          */
          SortedC1[c1.iCand].lfa=log(1+F1/2)-F1/2;
          SortedC2[c2.iCand].lfa=log(1+F2/2)-F2/2;
          
          thisCP->c1=c1.iCand;
          thisCP->c2=c2.iCand;
          thisCP->lfa=SortedC1[c1.iCand].lfa+SortedC2[c2.iCand].lfa;

        }
    }
  return 0;
} 

/*******************************************************************************/

/* Sorting function to sort candidate indices INCREASING order of f, delta, alpha */
int compareCIStructs(const void *a, const void *b)
{
  const CandidateList *ip = a;
  const CandidateList *jp = b;
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

/* Sorting function to sort candidate indices INCREASING order of f, delta, alpha */
int compareCdaf(const void *a, const void *b)
{
  const CandidateList *ip = a;
  const CandidateList *jp = b;

  REAL8 Delta1, Delta2;

  Delta1=ip->Delta;
  Delta2=jp->Delta;

  if (Delta1 < Delta2)
    return -1;
  
  if (Delta1 == Delta2)
    {
      REAL8 Alpha1, Alpha2;

      Alpha1=ip->Alpha;
      Alpha2=jp->Alpha;
          
      if (Alpha1 < Alpha2)
        return -1;
  
      if (Alpha1 > Alpha2)
        return 1;
  
      if (Alpha1 == Alpha2)
        {
          REAL8 freq1,freq2;

          freq1=ip->f;
          freq2=jp->f;

          if (freq1 < freq2)
            return -1;

          if (freq1 > freq2)
            return 1;

          /* This should never ever happen */
          if (freq1 == freq2)
            return 0;
        }
    }

  return 1;

}

/*******************************************************************************/

/* Sorting function to sort second candidate list into increasing order of fa */
int compareCPfa(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=CP[*(const int *)ip].lfa;
  dj=CP[*(const int *)jp].lfa;

  if (di<dj)
    return -1;
  
  if (di==dj)
    return (ip > jp);

  return 1;
}

/*******************************************************************************/

int ReadCandidateFiles(struct PolkaCommandLineArgsTag CLA)
{
  
  if (ReadOneCandidateFile ( &SortedC1, CLA.FstatsFile1)) return 1;
  CLength1=CLength;

  if (ReadOneCandidateFile ( &SortedC2, CLA.FstatsFile2)) return 1;
  CLength2=CLength;

  return 0;

} /* ReadCandidateFiles() */


/*******************************************************************************/

#define DONE_MARKER "%DONE\n"

/* read and parse the given candidate 'Fstats'-file fname into the candidate-list CList */
int  ReadOneCandidateFile (CandidateList **CList, const char *fname)
{
  UINT4 i;
  INT4 j, uplim;
  UINT4 numlines;
  REAL8 epsilon=1e-5;
  char line1[256];
  FILE *fp;
  INT4 read;
  UINT4 checksum=0;
  UINT4 bytecount=0;

  /* ------ Open and count candidates file ------ */
  i=0;
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      myPrintError("File %s doesn't exist!\n",fname);
      return 1;
    }
  while(fgets(line1,sizeof(line1),fp)) {
    unsigned int k;
    size_t len=strlen(line1);

    /* check that each line ends with a newline char (no overflow of
       line1 or null chars read) */
    if (!len || line1[len-1] != '\n') {
      myPrintError(
              "Line %d of file %s is too long or has no NEWLINE.  First 255 chars are:\n%s\n",
              i+1, fname, line1);
      fclose(fp);
      return 1;
    }

    /* increment line counter */
    i++;

    /* maintain a running checksum and byte count */
    bytecount+=len;
    for (k=0; k<len; k++)
      checksum+=(int)line1[k];
  }
  numlines=i;
  /* -- close candidate file -- */
  fclose(fp);     

  if ( numlines == 0) 
    {
      myPrintError ("ERROR: File '%s' has no lines so is not properly terminated by: %s", fname, DONE_MARKER);
      return 1;
    }

  /* output a record of the running checksun amd byte count */
  myPrintError( "%s: bytecount %" LAL_UINT4_FORMAT " checksum %" LAL_UINT4_FORMAT "\n", fname, bytecount, checksum);

  /* check validity of this Fstats-file */
  if ( strcmp(line1, DONE_MARKER ) ) 
    {
      myPrintError ("ERROR: File '%s' is not properly terminated by: %sbut has %s instead", fname, DONE_MARKER, line1);
      return 1;
    }
  else
    numlines --;        /* avoid stepping on DONE-marker */

  CLength=numlines;
  
  /* reserve memory for fstats-file contents */
  if (numlines > 0)
    {
      *CList = (CandidateList *)LALMalloc (numlines*sizeof(CandidateList));
      if ( !CList )
        {
          myPrintError ("Could not allocate memory for candidate file %s\n\n", fname);
          return 1;
        }
    }

  /* ------ Open and count candidates file ------ */
  i=0;
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      myPrintError("fopen(%s) failed!\n", fname);
      LALFree ((*CList));
      return 1;
    }
  while(i < numlines && fgets(line1,sizeof(line1),fp))
    {
      REAL8 dmp1, dmp2, dmp3;
      char newline='\0';
      CandidateList *cl=&(*CList)[i];

      if (strlen(line1)==0 || line1[strlen(line1)-1] != '\n') {
        myPrintError(
                "Line %d of file %s is too long or has no NEWLINE.  First 255 chars are:\n%s\n",
                i+1, fname, line1);
        LALFree ((*CList));
        fclose(fp);
        return 1;
      }
      
      cl->Ctag=0;
      cl->CtagCounter=-1;

      read = sscanf (line1, 
                     "%" LAL_REAL8_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL8_FORMAT 
                     " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL4_FORMAT "%c", 
                     &(cl->f), &(cl->Alpha), &(cl->Delta), &dmp1, &dmp2, &dmp3, &(cl->F), &newline );

      /* check that values that are read in are sensible */
      if (
          cl->f < 0.0                        ||
          cl->F < 0.0                        ||
          cl->Alpha <         0.0 - epsilon  ||
          cl->Alpha >   LAL_TWOPI + epsilon  ||
          cl->Delta < -0.5*LAL_PI - epsilon  ||
          cl->Delta >  0.5*LAL_PI + epsilon  ||                                                                 
          !isfinite(cl->f)                   ||
          !isfinite(cl->Alpha)               ||
          !isfinite(cl->Delta)               ||
          !isfinite(dmp1)                    ||
          !isfinite(dmp2)                    ||
          !isfinite(dmp3)                    ||
          !isfinite(cl->F)
          ) {
          myPrintError(
                  "Line %d of file %s has invalid values.\n"
                  "First 255 chars are:\n"
                  "%s\n"
                  "1st and 7th field should be positive.\n" 
                  "2nd field should lie between 0 and %1.15f.\n" 
                  "3rd field should lie between %1.15f and %1.15f.\n"
                  "All fields should be finite\n",
                  i+1, fname, line1, (double)LAL_TWOPI, (double)-LAL_PI/2.0, (double)LAL_PI/2.0);
          LALFree ((*CList));
          fclose(fp);
          return 1;
        }
           
           

      /* check that the FIRST character following the Fstat value is a
         newline.  Note deliberate LACK OF WHITE SPACE char before %c
         above */
      if (newline != '\n') {
        myPrintError(
                "Line %d of file %s had extra chars after F value and before newline.\n"
                "First 255 chars are:\n"
                "%s\n",
                i+1, fname, line1);
        LALFree ((*CList));
        fclose(fp);
        return 1;
      }

      /* check that we read 7 quantities with exactly the right format */
      if ( read != 8 )
        {
          myPrintError ("Found %d not %d values on line %d in file '%s'\n"
                         "Line in question is\n%s",
                         read, 8, i+1, fname, line1);               
          LALFree ((*CList));
          fclose(fp);
          return 1;
        }

      i++;
    }
  /* check that we read ALL lines! */
  if (i != numlines) {
    myPrintError(
            "Read of file %s terminated after %d line but numlines=%d\n",
            fname, i, numlines);
    LALFree((*CList));
    fclose(fp);
    return 1;
  }

  /* read final line with %DONE\n marker */
  if (!fgets(line1, sizeof(line1), fp)) {
    myPrintError(
            "Failed to find marker line of file %s\n",
            fname);
    LALFree((*CList));
    fclose(fp);
    return 1;
  }

  /* check for %DONE\n marker */
  if (strcmp(line1, DONE_MARKER)) {
    myPrintError(
            "Failed to parse marker: 'final' line of file %s contained %s not %s",
            fname, line1, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    return 1;
  }

  /* check that we are now at the end-of-file */
  if (fgetc(fp) != EOF) {
    myPrintError(
            "File %s did not terminate after %s",
            fname, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    return 1;
  }

  /* -- close candidate file -- */
  fclose(fp);     

  /* Finally, check that candidates are in the correct order */
  uplim = -1+(int)CLength;
  for (j=0; j<uplim; j++)
    {
      if(compareCdaf(*CList+j,*CList+j+1) != -1)
        {
          myPrintError( "Candidates in line %d and %d in Fstats file %s are not in the correct order. Exiting.\n", 
                  j+1,j+2, fname);
          LALFree ((*CList));
          return 1;
        }
    }

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
        myPrintError("Arguments are (defaults):\n");
        myPrintError("\t--fstatsfile1 (-1)\tSTRING\tFirst candidates Fstats file\n");
        myPrintError("\t--fstatsfile2 (-2)\tSTRING\tSecond candidates Fstats file\n");
        myPrintError("\t--outputfile  (-o)\tSTRING\tName of ouput candidates file\n");
        myPrintError("\t--frequency-window (-f)\tFLOAT\tFrequency window in Hz (0.0)\n");
        myPrintError("\t--alpha-window (-a)\tFLOAT\tAlpha window in radians (0.0)\n");
        myPrintError("\t--delta-window (-d)\tFLOAT\tDelta window in radians (0.0)\n");
        myPrintError("\t--fmin (-s)\tFLOAT\t Minimum frequency of candidate in 1st IFO\n");
        myPrintError("\t--fmax (-e)\tFLOAT\t Maximum frequency of candidate in 1st IFO\n");
        myPrintError("\t--EAHoutput (-b)\tFLAG\t Einstein at home output flag. \n");
        myPrintError("\t--help        (-h)\t\tThis message\n");
        exit(0);
        break;
      default:
        /* unrecognized option */
        errflg++;
        myPrintError("Unrecognized option argument %c\n",c);
        exit(1);
        break;
      }
    }

  if(CLA->FstatsFile1 == NULL)
    {
      myPrintError("No 1st candidates file specified; input with -1 option.\n");
      myPrintError("For help type %s -h\n", argv[0]);
      return 1;
    }      
  if(CLA->FstatsFile2 == NULL)
    {
      myPrintError("No 2nd candidates file specified; input with -2 option.\n");
      myPrintError("For help type %s -h\n", argv[0]);
      return 1;
    }      
  if(CLA->OutputFile == NULL)
    {
      myPrintError("No ouput filename specified; input with -o option.\n");
      myPrintError("For help type %s -h\n", argv[0]);
      return 1;
    }      

  if(CLA->fmin == 0.0)
    {
      myPrintError("No minimum frequency specified.\n");
      myPrintError("For help type %s -h\n", argv[0]);
      return 1;
    }      

  if(CLA->fmax == 0.0)
    {
      myPrintError("No maximum frequency specified.\n");
      myPrintError("For help type %s -h\n", argv[0]);
      return 1;
    }      

  return errflg;
}
