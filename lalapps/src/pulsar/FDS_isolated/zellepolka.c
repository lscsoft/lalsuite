/*********************************************************************************/
/*       uberpolka - the pulsar koinzidenz analysis code for einstein@home       */
/*                                                                               */
/*                   Xavier Siemens,  Bruce Allen,  Bernd Machenschalk           */
/*                                                                               */
/*                   (takes in one Fstats file to look for coincidence)          */
/*                                                                               */
/*                                  UWM - January  2005                          */
/*********************************************************************************/




/* ----------------------------------------------------------------------------- */
/* defines */
#ifndef FALSE
#define FALSE (1==0)
#endif
#ifndef TRUE
#define TRUE  (1==1)
#endif

#define DONE_MARKER "%DONE\n"

#ifndef USE_BOINC
#define USE_BOINC 0
#endif


/* ----------------------------------------------------------------------------- */
/* file includes */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "getopt.h"
#include <math.h>

#if HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/LALStatusMacros.h>
#include <lal/ConfigFile.h>
#include <lal/UserInput.h>

#include <lalapps.h>


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
#endif

/* this is defined in C99 and *should* be in math.h.  Long term
   protect this with a HAVE_FINITE */
#ifdef _MSC_VER
#include <float.h>
#define finite _finite
#else
int finite(double);
#endif


/* ----------------------------------------------------------------------------- */
/* some error codes and messages */
#define POLKAC_ENULL            1
#define POLKAC_ENONULL          2
#define POLKAC_ESYS             3
#define POLKAC_EINVALIDFSTATS   4
#define POLKAC_EMEM             5
#define POLKAC_ENORMAL          6

#define POLKAC_MSGENULL         "Arguments contained an unexpected null pointer"
#define POLKAC_MSGENONULL       "Input pointer was not NULL"
#define POLKAC_MSGESYS          "System call failed (probably file IO"
#define POLKAC_MSGEINVALIDFSTATS "Invalid Fstats file"
#define POLKAC_MSGEMEM          "Sorry, ran out of memory... bye."

#define UBERPOLKA_EXIT_ERRCLINE 31
#define UBERPOLKA_EXIT_READCND  32
#define UBERPOLKA_EXIT_FCTEST   33
#define UBERPOLKA_EXIT_OUTFAIL  34

#define POLKA_EXIT_OK 0




/* ----------------------------------------------------------------------------- */
/* structures */
struct PolkaCommandLineArgsTag 
{
  char *FstatsFile; /* Names of Fstat files to be read in */
  char *OutputFile;
  REAL8 Deltaf;      /* Size of coincidence window in Hz */
  REAL8 DeltaAlpha;  /* Size of coincidence window in radians */
  REAL8 DeltaDelta;  /* Size of coincidence window in radians */
  REAL8 Shiftf;      /* Parallel shift of frequency in Hz of cell */
  REAL8 ShiftAlpha;  /* Parallel shift of Alpha in Hz of cell */
  REAL8 ShiftDelta;  /* Parallel shift of Delta in Hz of cell */
  UINT4 EAH;         /* Einstein at home flag for alternative output */ 
} PolkaCommandLineArgs;

/* This structure contains the indices corresponding to the 
coarse frequency and sky bins */
typedef struct CandidateListTag
{
  UINT4 iCand;       /* Candidate id: unique with in this program.  */
  REAL8 f;           /* Frequency of the candidate */
  REAL8 Alpha;       /* right ascension of the candidate */
  REAL8 Delta;       /* declination  of the candidate */
  REAL8 F;           /* Maximum value of F for the cluster */
  INT4 FileID;       /* File ID to specify from which file the candidate under consideration originaly comes. */
  INT4  iFreq;       /* Frequency index */
  INT4 iDelta;       /* Declination index. This can be negative. */
  INT4 iAlpha;       /* Right ascension index */
  REAL8 lfa;         /* minus log of false alarm probability for that candidate */
} CandidateList;     /* ~ Fstat lines */ 

struct int4_linked_list {
  INT4 data;
  struct int4_linked_list *next;
}; 

typedef struct CellDataTag
{
  INT4 iFreq;          /* Frequency of the cell */
  INT4 iDelta;         /* Declination of the cell */
  INT4 iAlpha;         /* Right ascension of the cell */
  INT4 nCand;          /* number of the data in this cell. */
  REAL8 significance;  /* minus log of joint false alarm of the candidates in this cell. */
  struct int4_linked_list *CandID;  /* linked structure that has candidate id-s of the candidates in this cell. */
  REAL8 Freq;
  REAL8 Alpha;
  REAL8 Delta;
} CellData;


/* ----------------------------------------------------------------------------- */
/* Function declarelations */
void ReadCommandLineArgs(LALStatus *stat, int argc,char *argv[], struct PolkaCommandLineArgsTag *CLA); 
int ReadCandidateFiles(struct PolkaCommandLineArgsTag CLA);
int ReadOneCandidateFile(CandidateList **CList, const char *fname);
int compareCIStructs(const void *ip, const void *jp);
int rintcompare(INT4 *idata1, INT4 *idata2, size_t s); /* compare two INT4 arrays of size s.*/
int rfloatcompare(REAL8 *rdata1, REAL8 *rdata2, size_t s); /* compare two REAL8 arrays of size s.*/
void delete_int4_linked_list(struct int4_linked_list *list_ptr);
struct int4_linked_list *add_int4_data(struct int4_linked_list *list_ptr, INT4 *data);
void get_info_of_the_cell( CellData *cd, struct int4_linked_list *list_ptr );
int compareSignificance(const void *a, const void *b);
void FreeMemory(LALStatus *stat, struct PolkaCommandLineArgsTag *CLA, CellData *cell);

/*
int ReadCommandLine(int argc,char *argv[],struct PolkaCommandLineArgsTag *CLA);
*/

/* ----------------------------------------------------------------------------- */
/* Global variables. */
LALStatus global_status;
extern INT4 lalDebugLevel;
extern int vrbflg;

CandidateList *SortedC;
UINT4  CLength=0;




RCSID ("$Id$");


/* ------------------------------------------------------------------------------------------*/
/* Code starts here.                                                                         */
/* ------------------------------------------------------------------------------------------*/
/* main() mapped to polka() if using boinc */
#if USE_BOINC
int polka(int argc,char *argv[])
#else
int main(int argc,char *argv[]) 
#endif
{
  CellData *cell = NULL;
  UINT4 icell, icand;
#if USE_BOINC
  REAL8 local_fraction_done;
#endif

  LALStatus *stat = &global_status;
  lalDebugLevel = 0 ;  
  vrbflg = 1;   /* verbose error-messages */
  LAL_CALL (LALGetDebugLevel(stat, argc, argv, 'v'), stat);

  /* Reads command line arguments */
#if 0
  if (ReadCommandLine(argc,argv,&PolkaCommandLineArgs)) {
    fprintf(stderr,"ReadCommandLine failed\n");
    return UBERPOLKA_EXIT_ERRCLINE;
  }
#endif

  LAL_CALL( ReadCommandLineArgs( stat, argc,argv, &PolkaCommandLineArgs ), stat); 

  /* Reads in candidare files, set CLength */
  if (ReadCandidateFiles(PolkaCommandLineArgs)) {
    fprintf(stderr,"ReadCandidateFiles failed\n");
    return UBERPOLKA_EXIT_READCND;
  }

  if (CLength == 0  )
    {
      LALPrintError("CLength = %ud!\n",CLength);
      exit(1);
    }/* check that we have candidates. */

  cell = (CellData *) LALCalloc(sizeof(CellData),CLength);
  for(icell=0;icell<CLength;icell++) {
    cell[icell].CandID = (struct int4_linked_list *) LALCalloc(sizeof(struct int4_linked_list),1);
    cell[icell].CandID->next = NULL;
    cell[icell].iFreq = 0;
    cell[icell].iDelta = 0;
    cell[icell].iAlpha = 0;
    cell[icell].nCand = 0;
    cell[icell].Freq = 0.0;
    cell[icell].Delta = 0.0;
    cell[icell].Alpha = 0.0;
    cell[icell].significance = 0;
  }
  
  
  /* Initialise arrays of sorted candidates to use for bsearch */
  for (icand=0;icand<CLength;icand++)
    {
      SortedC[icand].iFreq=(INT4) (SortedC[icand].f/(PolkaCommandLineArgs.Deltaf) + PolkaCommandLineArgs.Shiftf  );
      SortedC[icand].iDelta=(INT4)(SortedC[icand].Delta/(PolkaCommandLineArgs.DeltaDelta)  + PolkaCommandLineArgs.ShiftDelta );
      SortedC[icand].iAlpha=(INT4)(SortedC[icand].Alpha*cos(SortedC[icand].Delta)/(PolkaCommandLineArgs.DeltaAlpha)  + PolkaCommandLineArgs.ShiftAlpha  );
    }
  
  /* sort arrays of candidates */
  qsort(SortedC, (size_t)CLength, sizeof(CandidateList), compareCIStructs);
  
  for (icand=0;icand<CLength;icand++) SortedC[icand].iCand=icand;


  /* initialization */
  icell = 0;
  icand = 0;
  cell[icell].iFreq = SortedC[icand].iFreq;
  cell[icell].iDelta = SortedC[icand].iDelta;
  cell[icell].iAlpha = SortedC[icand].iAlpha;
  cell[icell].CandID->data = SortedC[icand].iCand; 
  cell[icell].nCand = 1;

  /* ---------------------------------------------------------------------------------------------------------------*/      
  /* loop over candidates  */
  icell = 0;
  for (icand=1; icand < CLength; icand++)
    {
#if USE_BOINC
      /* make sure the cpu time is updated */ 
      if (boinc_time_to_checkpoint())
	boinc_checkpoint_completed();
#endif
      

      if( SortedC[icand].iFreq  == cell[icell].iFreq  && 
	  SortedC[icand].iDelta == cell[icell].iDelta &&
	  SortedC[icand].iAlpha == cell[icell].iAlpha ) 
	{ 
	  /* This candidate is in this cell. */
	  INT4 lastFileIDinThisCell = SortedC[cell[icell].CandID->data].FileID;
	  if( SortedC[icand].FileID != lastFileIDinThisCell ) 
	    {
	      /* This candidate has a different file id from the candidates in this cell. */
	      cell[icell].CandID = add_int4_data( cell[icell].CandID, &(SortedC[icand].iCand) );
	      cell[icell].nCand += 1;
	    } 
	  else 
	    {
	      /* This candidate has the same file id to one of candidates in this cell. */
	      /* Because the array is already sorted in the DECREASING ORDER OF F, 
		 we do nothing here. */
	    }
	} 
      else 
	{	  
	  /* This candidate is outside of this cell. */
	  icell++;
	  cell[icell].iFreq = SortedC[icand].iFreq;
	  cell[icell].iDelta = SortedC[icand].iDelta;
	  cell[icell].iAlpha = SortedC[icand].iAlpha;
	  cell[icell].CandID->data = SortedC[icand].iCand; 
	  cell[icell].nCand = 1;
	}
           
#if USE_BOINC
      local_fraction_done = 0.99 + 0.01 * (double)i / (double)CLength;
      /* update progress, the last % is reserved for polka */
      boinc_fraction_done(local_fraction_done);
      /* pass variable externally to graphics routines */
      if (fraction_done_hook != NULL)
	*fraction_done_hook = local_fraction_done;
#endif
    }/* loop over candidate list */      
  /* ---------------------------------------------------------------------------------------------------------------*/      



  /* ---------------------------------------------------------------------------------------------------------------*/      

  /* sort arrays of candidates */
  { 
    FILE *fp = NULL;
    UINT4 ncell=icell+1;

    for(icell=0;icell<ncell;icell++) {
      get_info_of_the_cell( &cell[icell], cell[icell].CandID );
    }

    qsort(cell, (size_t)ncell, sizeof(CellData), compareSignificance);

    fp = fopen(PolkaCommandLineArgs.OutputFile,"w");
    for(icell=0;icell<ncell;icell++) {
      fprintf(fp,"%" LAL_REAL4_FORMAT "\t%" LAL_REAL4_FORMAT "\t%" LAL_REAL4_FORMAT "\t" 
	      "%10d %10d %10d %10d %22.12f\n",
	      cell[icell].Freq,cell[icell].Delta,cell[icell].Alpha,
	      cell[icell].iFreq,cell[icell].iDelta,cell[icell].iAlpha,
	      cell[icell].nCand,
	      cell[icell].significance);
    }
    fclose(fp);
  }


  LAL_CALL( FreeMemory(stat, &PolkaCommandLineArgs, cell), stat);

  LALCheckMemoryLeaks(); 

  return 0;
 
}

/*******************************************************************************/
/* Free memory */
void FreeMemory(LALStatus *stat, struct PolkaCommandLineArgsTag *CLA, CellData *cell)
{
  INITSTATUS( stat, "FreeMemory", rcsid );
  ATTATCHSTATUSPTR (stat);

  UINT4 icell;

  LALFree(CLA->FstatsFile);
  LALFree(CLA->OutputFile);

  LALFree(SortedC);

  for(icell=0;icell<CLength;icell++) {
    delete_int4_linked_list( cell[icell].CandID );
  }
  LALFree(cell);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
}



/*******************************************************************************/
/* add data to linked structure */
struct int4_linked_list *add_int4_data(struct int4_linked_list *list_ptr, INT4 *data)
{
  struct int4_linked_list *p = NULL;
  p = (struct int4_linked_list *) LALMalloc(sizeof(struct int4_linked_list));
  p->data = *(data);
  p->next = list_ptr;
  list_ptr = p;
  return list_ptr;
}


/* delete data to linked structure */
void delete_int4_linked_list(struct int4_linked_list *list_ptr)
{
  struct int4_linked_list *q;
  while( list_ptr !=NULL ) {
    q = list_ptr->next;
    LALFree( list_ptr );
    list_ptr = q;
  }
  return;
}


/* get info of this cell. */
void get_info_of_the_cell( CellData *cd, struct int4_linked_list *list_ptr )
{
  INT4 idx;
  struct int4_linked_list *p;
  p = list_ptr;
  while( p !=NULL ) {
    idx = p->data;
    cd->significance += SortedC[idx].lfa;
    cd->Alpha += SortedC[idx].Alpha;
    cd->Delta += SortedC[idx].Delta;
    cd->Freq += SortedC[idx].f;
    p = p->next;
  }

  cd->Alpha /= cd->nCand;
  cd->Delta /= cd->nCand;
  cd->Freq  /= cd->nCand;
  return;
}



/*******************************************************************************/
/* Sorting function to sort candidate indices INCREASING order of f, delta, alpha, and 
   DECREASING ORDER OF F STATISTIC. */
int compareCIStructs(const void *a, const void *b)
{
  const CandidateList *ip = a;
  const CandidateList *jp = b;
  int res;
  INT4 ap[4],bp[4];

  ap[0]=ip->iFreq;
  ap[1]=ip->iDelta;
  ap[2]=ip->iAlpha;
  ap[3]=ip->FileID;

  bp[0]=jp->iFreq;
  bp[1]=jp->iDelta;
  bp[2]=jp->iAlpha;
  bp[3]=jp->FileID;

  res = rintcompare( ap,  bp, 4);
  if( res == 0 ) {
    REAL8 F1, F2;
    F1=ip->F;
    F2=jp->F;
    /* I put F1 and F2 inversely, because I would like to get decreasingly-ordered set. */ 
    res = rfloatcompare( &F2,  &F1, 1);
  } 
  return res;
} /* int compareCIStructs() */




int compareSignificance(const void *a, const void *b)
{
  const CellData *ip = a;
  const CellData *jp = b;
  int res;

  REAL8 F1, F2;
  F1=ip->significance;
  F2=jp->significance;
  /* I put F1 and F2 inversely, because I would like to get decreasingly-ordered set. */ 
  res = rfloatcompare( &F2,  &F1, 1);
  return res;
} /* int compareCIStructs() */






int rfloatcompare(REAL8 *ap, REAL8 *bp, size_t n) {
  if( (*ap) == (*bp) ) { 
    if ( n > 1 ){  
      return rfloatcompare( ap+1, bp+1, n-1 );
    } else {
      return 0;
    }
  }
  if ( (*ap) < (*bp) ) 
    return -1;    
  return 1;
} /* int rfloatcompare() */


int rintcompare(INT4 *ap, INT4 *bp, size_t n) {
  if( (*ap) == (*bp) ) { 
    if ( n > 1 ){  
      return rintcompare( ap+1, bp+1, n-1 );
    } else {
      return 0;
    }
  }
  if ( (*ap) < (*bp) ) 
    return -1;    
  return 1;
} /* int rintcompare() */



/*******************************************************************************/
/* We would use glob in the future to read different files ? */
int ReadCandidateFiles(struct PolkaCommandLineArgsTag CLA)
{
  if (ReadOneCandidateFile ( &SortedC, CLA.FstatsFile)) return 1;

  return 0;

} /* ReadCandidateFiles() */


/*******************************************************************************/



/* read and parse the given candidate 'Fstats'-file fname into the candidate-list CList */
int  ReadOneCandidateFile (CandidateList **CList, const char *fname)
{
  UINT4 i;
  /*  INT4 j, uplim; */
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
      LALPrintError("File %s doesn't exist!\n",fname);
      return 1;
    }
  while(fgets(line1,sizeof(line1),fp)) {
    unsigned int k;
    size_t len=strlen(line1);

    /* check that each line ends with a newline char (no overflow of
       line1 or null chars read) */
    if (!len || line1[len-1] != '\n') {
      LALPrintError(
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
      LALPrintError ("ERROR: File '%s' has no lines so is not properly terminated by: %s", fname, DONE_MARKER);
      return 1;
    }

  /* output a record of the running checksun amd byte count */
  LALPrintError( "%s: bytecount %" LAL_UINT4_FORMAT " checksum %" LAL_UINT4_FORMAT "\n", fname, bytecount, checksum);

  /* check validity of this Fstats-file */
  if ( strcmp(line1, DONE_MARKER ) ) 
    {
      LALPrintError ("ERROR: File '%s' is not properly terminated by: %sbut has %s instead", fname, DONE_MARKER, line1);
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
          LALPrintError ("Could not allocate memory for candidate file %s\n\n", fname);
          return 1;
        }
    }

  /* ------ Open and count candidates file ------ */
  i=0;
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      LALPrintError("fopen(%s) failed!\n", fname);
      LALFree ((*CList));
      return 1;
    }
  while(i < numlines && fgets(line1,sizeof(line1),fp))
    {
      char newline='\0';
      CandidateList *cl=&(*CList)[i];

      if (strlen(line1)==0 || line1[strlen(line1)-1] != '\n') {
        LALPrintError(
                "Line %d of file %s is too long or has no NEWLINE.  First 255 chars are:\n%s\n",
                i+1, fname, line1);
        LALFree ((*CList));
        fclose(fp);
        return 1;
      }
      
      read = sscanf (line1, 
                     "%" LAL_INT4_FORMAT "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
                     "%c", 
                     &(cl->FileID), &(cl->f), &(cl->Alpha), &(cl->Delta), &(cl->F), &newline );

      /* check that values that are read in are sensible */
      if (
          cl->FileID < 0                        ||
          cl->f < 0.0                        ||
          cl->F < 0.0                        ||
          cl->Alpha <         0.0 - epsilon  ||
          cl->Alpha >   LAL_TWOPI + epsilon  ||
          cl->Delta < -0.5*LAL_PI - epsilon  ||
          cl->Delta >  0.5*LAL_PI + epsilon  ||
	  !finite(cl->FileID)                     ||                                                                 
          !finite(cl->f)                     ||
          !finite(cl->Alpha)                 ||
          !finite(cl->Delta)                 ||
          !finite(cl->F)
          ) {
          LALPrintError(
                  "Line %d of file %s has invalid values.\n"
                  "First 255 chars are:\n"
                  "%s\n"
                  "1st and 4th field should be positive.\n" 
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
        LALPrintError(
                "Line %d of file %s had extra chars after F value and before newline.\n"
                "First 255 chars are:\n"
                "%s\n",
                i+1, fname, line1);
        LALFree ((*CList));
        fclose(fp);
        return 1;
      }

      /* check that we read 6 quantities with exactly the right format */
      if ( read != 6 )
        {
          LALPrintError ("Found %d not %d values on line %d in file '%s'\n"
                         "Line in question is\n%s",
                         read, 6, i+1, fname, line1);               
          LALFree ((*CList));
          fclose(fp);
          return 1;
        }


      cl->lfa = cl->F - log( 1.0 + cl->F );

      i++;
    } /*  end of main while loop */
  /* check that we read ALL lines! */
  if (i != numlines) {
    LALPrintError(
            "Read of file %s terminated after %d line but numlines=%d\n",
            fname, i, numlines);
    LALFree((*CList));
    fclose(fp);
    return 1;
  }

  /* read final line with %DONE\n marker */
  if (!fgets(line1, sizeof(line1), fp)) {
    LALPrintError(
            "Failed to find marker line of file %s\n",
            fname);
    LALFree((*CList));
    fclose(fp);
    return 1;
  }

  /* check for %DONE\n marker */
  if (strcmp(line1, DONE_MARKER)) {
    LALPrintError(
            "Failed to parse marker: 'final' line of file %s contained %s not %s",
            fname, line1, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    return 1;
  }

  /* check that we are now at the end-of-file */
  if (fgetc(fp) != EOF) {
    LALPrintError(
            "File %s did not terminate after %s",
            fname, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    return 1;
  }

  /* -- close candidate file -- */
  fclose(fp);     

  return 0;

} /* ReadOneCandidateFile() */



/*******************************************************************************/
#if 0
int ReadCommandLine(int argc,char *argv[],struct PolkaCommandLineArgsTag *CLA) 
{
  INT2 errflg = 0;
  INT4 c; 
  INT4 option_index = 0;

  const char *optstring = "hI:f:a:d:F:D:A:o:E";
  struct option long_options[] =
    {
      {"input-data-file",       required_argument, 0,   'I'},
      {"frequency-window",      required_argument, 0,   'f'},
      {"delta-window",          required_argument, 0,   'd'},
      {"alpha-window",          required_argument, 0,   'a'},
      {"frequency-shift",       required_argument, 0,   'F'},
      {"delta-shift",           required_argument, 0,   'D'},
      {"alpha-shift",           required_argument, 0,   'A'},
      {"outputfile",            required_argument, 0,   'o'},
      {"EAHoutput",             no_argument, 0,         'E'},
      {"help",                  no_argument, 0,         'h'},
      {0, 0, 0, 0}
    };

  /* Initialize default values */
  CLA->FstatsFile=NULL;
  CLA->OutputFile=NULL;
  CLA->Deltaf=0.0;
  CLA->DeltaAlpha=0;
  CLA->DeltaDelta=0;

  CLA->Shiftf=0.0;
  CLA->ShiftAlpha=0;
  CLA->ShiftDelta=0;

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
      case 'I':
        /* SFT directory */
        CLA->FstatsFile=optarg;
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
      case 'd':
        /* Spin down order */
        CLA->DeltaDelta=atof(optarg);
        break;
      case 'F':
        /* Spin down order */
        CLA->Shiftf=atof(optarg);
        break;
      case 'A':
        /* Spin down order */
        CLA->ShiftAlpha=atof(optarg);
        break;
      case 'D':
        /* Spin down order */
        CLA->ShiftDelta=atof(optarg);
        break;
      case 'E':
        /* Spin down order */
        CLA->EAH=1;
        break;
      case 'h':
        /* print usage/help message */
        LALPrintError("Arguments are (defaults):\n");
        LALPrintError("\t--fstatsfile (-I)\tSTRING\tFirst candidates Fstats file\n");
        LALPrintError("\t--outputfile  (-o)\tSTRING\tName of ouput candidates file\n");
        LALPrintError("\t--frequency-window (-f)\tFLOAT\tFrequency window in Hz (0.0)\n");
        LALPrintError("\t--alpha-window (-a)\tFLOAT\tAlpha window in radians (0.0)\n");
        LALPrintError("\t--delta-window (-d)\tFLOAT\tDelta window in radians (0.0)\n");
        LALPrintError("\t--frequency-shift (-F)\tFLOAT\tFrequency window in Hz (0.0)\n");
        LALPrintError("\t--alpha-shift (-A)\tFLOAT\tAlpha window in radians (0.0)\n");
        LALPrintError("\t--delta-shift (-D)\tFLOAT\tDelta window in radians (0.0)\n");
        LALPrintError("\t--EAHoutput (-E)\tFLAG\t Einstein at home output flag. \n");
        LALPrintError("\t--help        (-h)\t\tThis message\n");
        exit(0);
        break;
      default:
        /* unrecognized option */
        errflg++;
        LALPrintError("Unrecognized option argument %c\n",c);
        exit(1);
        break;
      }
    }

  if(CLA->FstatsFile == NULL)
    {
      LALPrintError("No candidates file specified; input with -I option.\n");
      LALPrintError("For help type %s -h\n", argv[0]);
      return 1;
    }      
  if(CLA->OutputFile == NULL)
    {
      LALPrintError("No ouput filename specified; input with -o option.\n");
      LALPrintError("For help type %s -h\n", argv[0]);
      return 1;
    }      


  return errflg;
}
#endif


void ReadCommandLineArgs(LALStatus *stat, int argc,char *argv[], struct PolkaCommandLineArgsTag *CLA) 
{
  INITSTATUS( stat, "ReadCommandLineArgs", rcsid );
  ATTATCHSTATUSPTR (stat);

  ASSERT( CLA != NULL, stat, POLKAC_ENULL, POLKAC_MSGENULL);

  CHAR* uvar_InputData;
  CHAR* uvar_OutputData;
  REAL8 uvar_FreqWindow;
  REAL8 uvar_AlphaWindow;
  REAL8 uvar_DeltaWindow;
  REAL8 uvar_FreqShift;
  REAL8 uvar_AlphaShift;
  REAL8 uvar_DeltaShift;
  BOOLEAN uvar_help;
  BOOLEAN uvar_EAHoutput;


  uvar_help = 0;
  uvar_EAHoutput = 0;

  uvar_InputData = NULL;
  uvar_OutputData = NULL;

  uvar_FreqWindow = 0.0;
  uvar_AlphaWindow = 0.0;
  uvar_DeltaWindow = 0.0;

  uvar_FreqShift = 0.0;
  uvar_AlphaShift = 0.0;
  uvar_DeltaShift = 0.0;

  /* register all our user-variables */
  LALregBOOLUserVar(stat,       help,           'h', UVAR_HELP,     "Print this message"); 
  LALregBOOLUserVar(stat,       EAHoutput,      'E', UVAR_OPTIONAL, "Einstein at Home output"); 

  LALregSTRINGUserVar(stat,     InputData,      'I', UVAR_REQUIRED, "Input candidates Fstats file.");
  LALregSTRINGUserVar(stat,     OutputData,     'o', UVAR_REQUIRED, "Ouput candidates file name");

  LALregREALUserVar(stat,       FreqWindow,     'f', UVAR_REQUIRED, "Frequency window in Hz");
  LALregREALUserVar(stat,       AlphaWindow,    'a', UVAR_REQUIRED, "Right ascension window in radians");
  LALregREALUserVar(stat,       DeltaWindow,    'd', UVAR_REQUIRED, "Declination window in radians");
  LALregREALUserVar(stat,       FreqShift,      'F', UVAR_OPTIONAL, "Frequency shift in FreqWindow");
  LALregREALUserVar(stat,       AlphaShift,     'A', UVAR_OPTIONAL, "Right ascension shift in AlphaWindow");
  LALregREALUserVar(stat,       DeltaShift,     'D', UVAR_OPTIONAL, "Declination shift in DeltaWindow");

  TRY (LALUserVarReadAllInput(stat->statusPtr,argc,argv),stat); 

  if (uvar_help) {	/* if help was requested, we're done here */
    LALDestroyUserVars(stat->statusPtr);
    exit(POLKA_EXIT_OK);
  }


  CLA->FstatsFile = (CHAR *) LALMalloc(strlen(uvar_InputData)+1);
  if(CLA->FstatsFile == NULL)
    {
      LALPrintError("No candidates file specified; input with -I option.\n");
      LALPrintError("For help type %s -h\n", argv[0]);
      RETURN (stat);
    }      
  CLA->OutputFile = (CHAR *) LALMalloc(strlen(uvar_OutputData)+1);
  if(CLA->OutputFile == NULL)
    {
      LALFree(CLA->FstatsFile); 
      LALPrintError("No ouput filename specified; input with -o option.\n");
      LALPrintError("For help type %s -h\n", argv[0]);
      RETURN (stat);
    }      

  strcpy(CLA->FstatsFile,uvar_InputData);
  strcpy(CLA->OutputFile,uvar_OutputData);

  CLA->Deltaf = uvar_FreqWindow;
  CLA->DeltaAlpha = uvar_AlphaWindow;
  CLA->DeltaDelta = uvar_DeltaWindow;

  CLA->Shiftf = uvar_FreqShift;
  CLA->ShiftAlpha = uvar_AlphaShift;
  CLA->ShiftDelta = uvar_DeltaShift;

  CLA->EAH = uvar_EAHoutput;

  LALDestroyUserVars(stat->statusPtr);
  BEGINFAIL(stat) {
    LALFree(CLA->FstatsFile);
    LALFree(CLA->OutputFile);
  } ENDFAIL(stat);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
}

