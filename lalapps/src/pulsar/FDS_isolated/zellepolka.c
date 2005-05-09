/********************************************************************************************/
/*       zellepolka - the pulsar koinzidenz analysis code for einstein@home postprocess     */
/*                                                                                          */
/*               Xavier Siemens,  Bruce Allen,  Bernd Machenschalk,  Yousuke Itoh           */
/*                                                                                          */
/*                   (takes in one Fstats files to look for coincidence)                    */
/*                                                                                          */
/*                                  UWM - April  2005                                       */
/*                             Based on uberpolka written by                                */
/*                        Xavier Siemens, Bruce Allen, Bernd Machenschalk                   */
/********************************************************************************************/



/* ----------------------------------------------------------------------------- */
/* defines */
#ifndef FALSE
#define FALSE (1==0)
#endif
#ifndef TRUE
#define TRUE  (1==1)
#endif

#define DONE_MARKER "%DONE\n"
/* maximum depth of a linked structure. */
#define LINKEDSTR_MAX_DEPTH 1024 

#ifndef USE_BOINC
#define USE_BOINC 0
#endif



/* ----------------------------------------------------------------------------- */
/* file includes */
#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "getopt.h"
#include <math.h>

#include <unistd.h>


/* 
   To use unzip, you need to have unzip-5.5x from, say,  the InfoZip webpage, 
   and readzipfile_util.h and .c. from yousuke. 
#define USE_UNZIP  
*/




#ifdef USE_UNZIP
#include "unzip.h"
#include "readzipfile_util.h"
#endif

#ifdef HAVE_GLOB_H
#include <glob.h>
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
#define POLKAC_EUNZIP           7
#define POLKAC_EGLOB           8

#define POLKAC_MSGENULL         "Arguments contained an unexpected null pointer"
#define POLKAC_MSGENONULL       "Input pointer was not NULL"
#define POLKAC_MSGESYS          "System call failed (probably file IO"
#define POLKAC_MSGEINVALIDFSTATS "Invalid Fstats file"
#define POLKAC_MSGEMEM          "Sorry, ran out of memory... bye."
#define POLKAC_MSGEUNZIP          "Cannot use unzip."
#define POLKAC_MSGEGLOB          "Cannot use glob."


#define POLKA_EXIT_OK 0
#define POLKA_EXIT_ERR     31
#define POLKA_EXIT_READCND  32
#define POLKA_EXIT_FCTEST   33
#define POLKA_EXIT_OUTFAIL  34


/* ----------------------------------------------------------------------------- */
/* structures */
typedef struct PolkaConfigVarsTag 
{
  CHAR *FstatsFile;  /**<  Names of Fstat files to be read in */
  CHAR *OutputFile;  /**<  Names of output file */
  CHAR *InputDir;    /**<  Directory name of input files */
  CHAR *BaseName;    /**<  Base name of input files */
  CHAR **Filelist;   /**<  Array of filenames to load Fstats file from */
  UINT4 NFiles;      /**<  Number of input files read */
  INT4 Nthr;        /**<  Show exective results of cells with numbers of coincidence above Nthr. */
  REAL4 Sthr;        /**<  Show exective results of cells with significance above Sthr. */
  BOOLEAN AutoOut; 
  REAL8 Deltaf;      /**<  Size of coincidence window in Hz */
  REAL8 DeltaAlpha;  /**<  Size of coincidence window in radians */
  REAL8 DeltaDelta;  /**<  Size of coincidence window in radians */
  REAL8 Shiftf;      /**<  Parallel shift of frequency in Hz of cell */
  REAL8 ShiftAlpha;  /**<  Parallel shift of Alpha in Hz of cell */
  REAL8 ShiftDelta;  /**<  Parallel shift of Delta in Hz of cell */
} PolkaConfigVariables;

/* This structure contains the indices corresponding to the 
coarse frequency and sky bins */
typedef struct CandidateListTag
{
  UINT4 iCand;       /**<  Candidate id: unique with in this program.  */
  REAL8 f;           /**<  Frequency of the candidate */
  REAL8 Alpha;       /**<  right ascension of the candidate */
  REAL8 Delta;       /**<  declination  of the candidate */
  REAL8 TwoF;           /**<  Maximum value of F for the cluster */
  INT4 FileID;       /**<  File ID to specify from which file the candidate under consideration originaly comes. */
  INT4  iFreq;       /**<  Frequency index */
  INT4 iDelta;       /**<  Declination index. This can be negative. */
  INT4 iAlpha;       /**<  Right ascension index */
} CandidateList;     /**<  ~ Fstat lines */ 

struct int4_linked_list {
  INT4 data;
  struct int4_linked_list *next;
}; 

typedef struct CellDataTag
{
  INT4 iFreq;          /**<  Frequency of the cell */
  INT4 iDelta;         /**<  Declination of the cell */
  INT4 iAlpha;         /**<  Right ascension of the cell */
  INT4 nCand;          /**<  number of the data in this cell. */
  REAL8 significance;  /**<  minus log of joint false alarm of the candidates in this cell. */
  struct int4_linked_list *CandID;  /**<  linked structure that has candidate id-s of the candidates in this cell. */
  REAL8 Freq;
  REAL8 Alpha;
  REAL8 Delta;
} CellData;


/* ----------------------------------------------------------------------------- */
/* Function declarelations */
void ReadCommandLineArgs(LALStatus *stat, int argc,char *argv[], struct PolkaConfigVarsTag *CLA); 
void ReadCandidateFiles(LALStatus *stat, CandidateList **Clist, struct PolkaConfigVarsTag *CLA, UINT4 *datalen);
void ReadOneCandidateFile(LALStatus *stat, CandidateList **CList, const char *fname, UINT4 *datalen);
void PrepareCells( LALStatus *stat, CellData **cell, UINT4 datalen );
int compareCandidates(const void *ip, const void *jp);
int compareSignificance(const void *a, const void *b);
int rintcompare(INT4 *idata1, INT4 *idata2, size_t s); /* compare two INT4 arrays of size s.*/
int rfloatcompare(REAL8 *rdata1, REAL8 *rdata2, size_t s); /* compare two REAL8 arrays of size s.*/
void add_int4_data(LALStatus *stat, struct int4_linked_list **list_ptr, INT4 *data);
void delete_int4_linked_list(struct int4_linked_list *list_ptr);
void get_info_of_the_cell( CellData *cd, CandidateList *CList);
void PrintResult(LALStatus *stat, struct PolkaConfigVarsTag *CLA, CellData *cell, UINT4 *ncell, CandidateList *CList);
void FreeMemory(LALStatus *stat, struct PolkaConfigVarsTag *CLA, CellData *cell, CandidateList *CList, UINT4 datalen);
void FreeConfigVars(LALStatus *stat, struct PolkaConfigVarsTag *CLA );
void GetFilesListInThisDir(LALStatus *stat, CHAR *directory, CHAR *basename, CHAR ***filelist, UINT4 *nfiles );
void print_Fstat_of_the_cell( FILE *fp, CellData *cd, CandidateList *CList, INT4 icell_start, INT4 icell_end, REAL8 sig_thr, REAL8 ncand_thr );
int compareNCandidate(const void *a, const void *b);
void print_info_of_the_cell( FILE *fp, CellData *cd, INT4 icell_start, INT4 icell_end, REAL8 sig_thr, REAL8 ncand_thr);


#ifdef USE_UNZIP
void ReadCandidateListFromZipFile (LALStatus *stat, CandidateList **CList, CHAR *fname, UINT4 *candlen, INT4 *FileID);
#endif

/* ----------------------------------------------------------------------------- */
/* Global variables. */
LALStatus global_status;
extern INT4 lalDebugLevel;
extern int vrbflg;

RCSID ("$Id$");



/* ------------------------------------------------------------------------------------------*/
/* Code starts here.                                                                         */
/* ------------------------------------------------------------------------------------------*/
/* ########################################################################################## */
/* main() mapped to polka() if using boinc */
#if USE_BOINC
int polka(int argc,char *argv[])
#else
int main(int argc,char *argv[]) 
#endif
{
  LALStatus *stat = &global_status;
  lalDebugLevel = 0 ;  
  vrbflg = 1;   /* verbose error-messages */

  UINT4  CLength=0;
  CandidateList *SortedC = NULL;
  CellData *cell = NULL;
  UINT4 icell, icand, ncell;
#if USE_BOINC
  REAL8 local_fraction_done;
#endif

  PolkaConfigVariables PolkaConfigVars;


  LAL_CALL (LALGetDebugLevel(stat, argc, argv, 'v'), stat);

  /* Reads command line arguments */
  LAL_CALL( ReadCommandLineArgs( stat, argc,argv, &PolkaConfigVars ), stat); 

  /* Reads in candidare files, set CLength */
  LAL_CALL( ReadCandidateFiles(stat, &SortedC, &PolkaConfigVars, &CLength), stat);

  /* Prepare cells. */
  LAL_CALL( PrepareCells( stat, &cell, CLength ), stat);  



  /* --------------------------------------------------------------------------------*/      
  /* initialization */
  /* Initialise arrays of sorted candidates. */
  for (icand=0;icand<CLength;icand++)
    {
      SortedC[icand].iFreq=(INT4) (SortedC[icand].f/(PolkaConfigVars.Deltaf) + PolkaConfigVars.Shiftf  );
      SortedC[icand].iDelta=(INT4)(SortedC[icand].Delta/(PolkaConfigVars.DeltaDelta)  + PolkaConfigVars.ShiftDelta );
      SortedC[icand].iAlpha=(INT4)(SortedC[icand].Alpha*cos(SortedC[icand].Delta)/(PolkaConfigVars.DeltaAlpha)  + PolkaConfigVars.ShiftAlpha  );
      SortedC[icand].iCand=icand; /* Keep the original ordering before sort to refer the orignal data later. */
    }

  /* sort arrays of candidates */
  qsort(SortedC, (size_t)CLength, sizeof(CandidateList), compareCandidates);
  

  /* Initialise the first cell by the first candidate. */
  icell = 0;
  icand = 0;
  cell[icell].iFreq = SortedC[icand].iFreq;
  cell[icell].iDelta = SortedC[icand].iDelta;
  cell[icell].iAlpha = SortedC[icand].iAlpha;
  cell[icell].CandID->data = icand; 
  cell[icell].nCand = 1;

  /* ------------------------------------------------------------------------------*/      
  /* main loop over candidates  */
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
	      add_int4_data( stat, &(cell[icell].CandID), &(icand) );
	      cell[icell].nCand += 1;
	    } 
	  else 
	    {
	      /* This candidate has the same file id to one of candidates in this cell. */
	      /* Because the array is already sorted in the DECREASING ORDER OF F, 
		 we do nothing here. */
	    } /* if( SortedC[icand].FileID != lastFileIDinThisCell ) */
	} /*  if( SortedC[icand].iFreq  == cell[icell].iFreq  && .. ) */ 
      else 
	{	  
	  /* This candidate is outside of this cell. */
	  icell++;
	  cell[icell].iFreq = SortedC[icand].iFreq;
	  cell[icell].iDelta = SortedC[icand].iDelta;
	  cell[icell].iAlpha = SortedC[icand].iAlpha;
	  cell[icell].CandID->data = icand;
	  cell[icell].nCand = 1;
	} /*  if( SortedC[icand].iFreq  == cell[icell].iFreq  && .. ) */ 

           
#if USE_BOINC
      local_fraction_done = 0.99 + 0.01 * (double)i / (double)CLength;
      /* update progress, the last % is reserved for polka */
      boinc_fraction_done(local_fraction_done);
      /* pass variable externally to graphics routines */
      if (fraction_done_hook != NULL)
	*fraction_done_hook = local_fraction_done;
#endif
    } /* for (icand=1; icand < CLength; icand++): loop over candidate list */      
  /* ---------------------------------------------------------------------------------------- */      


  /* Get the information in each cell. */
  ncell=icell+1; /* number of the cells in which more than or at least one candidate exists. */
  for(icell=0;icell<ncell;icell++) {
    get_info_of_the_cell( &cell[icell], SortedC);
  }  


  /* -----------------------------------------------------------------------------------------*/      
  /* Output results */
  LAL_CALL( PrintResult( stat, &PolkaConfigVars, cell, &ncell, SortedC),stat );



  /* -----------------------------------------------------------------------------------------*/      
  /* Clean-up */
  LAL_CALL( FreeMemory(stat, &PolkaConfigVars, cell, SortedC, CLength), stat);

  LALCheckMemoryLeaks(); 

  return 0;
 
} /* main() */


/* ########################################################################################## */
/* ------------------------------------------------------------------------------*/      
/* Initialize the code: allocate memory, set initial values.                     */
/* ------------------------------------------------------------------------------*/      
void PrepareCells( LALStatus *stat, CellData **cell, UINT4 CLength )
{
  INITSTATUS( stat, "InitCode", rcsid );
  ATTATCHSTATUSPTR (stat);

  UINT4 icell, ncell;
  INT4 errflg = 0;

  *cell = (CellData *) LALCalloc(sizeof(CellData),CLength);
  if( *cell == NULL ) {
    ABORT (stat, POLKAC_EMEM, POLKAC_MSGEMEM);
  }

  for(icell=0;icell<CLength;icell++) {
    (*cell)[icell].CandID = NULL;
    (*cell)[icell].CandID = (struct int4_linked_list *) LALCalloc(sizeof(struct int4_linked_list),1);
    if( (*cell)[icell].CandID == NULL ) {
      errflg = 1;
      break;
    }
    (*cell)[icell].CandID->next = NULL;
    (*cell)[icell].iFreq = 0;
    (*cell)[icell].iDelta = 0;
    (*cell)[icell].iAlpha = 0;
    (*cell)[icell].nCand = 0;
    (*cell)[icell].Freq = 0.0;
    (*cell)[icell].Delta = 0.0;
    (*cell)[icell].Alpha = 0.0;
    (*cell)[icell].significance = 0;
  }


  if( errflg != 0 ) {
    ncell = icell;
    for(icell=0;icell<ncell;icell++) {
      LALFree( (*cell)[icell].CandID );
    }
    ABORT (stat, POLKAC_EMEM, POLKAC_MSGEMEM);
  }


  DETATCHSTATUSPTR (stat);
  RETURN (stat);
} /* PrepareCells() */


/* ########################################################################################## */
/* Output results */
void PrintResult(LALStatus *stat, struct PolkaConfigVarsTag *CLA, CellData *cell, UINT4 *ncell, CandidateList *CList)
{
  INITSTATUS( stat, "PrintResult", rcsid );
  UINT4 icell;
  CHAR fnameSigTime[]="polka_significant_outlier_2FofTime"; /* Time variation of 2F of some significant outliers. */
  CHAR fnameSigCell[]="polka_significant_outlier_CellData"; /* Cell information of some significant outliers*/
  CHAR fnameCoiTime[]="polka_coincident_outlier_2FofTime";  /* Time variation of 2F of some coincident outliers. */
  CHAR fnameCoiCell[]="polka_coincident_outlier_CellData";  /* Cell information of some coincident outliers*/
  FILE *fp = NULL, *fpSigTime = NULL, *fpSigCell = NULL, *fpCoiTime = NULL, *fpCoiCell = NULL;
  INT4 *count;
  INT4 nc, nmax,idxmax = 0;
  REAL4 Sigmax = 0.0;

  /* ------------------------------------------------------------- */
  /* First Sort arrays of candidates based on number of candidate. */ 
  qsort(cell, (size_t) (*ncell), sizeof(CellData), compareNCandidate);


  nmax = cell[0].nCand; /* This is the number of the maximum coincidences. */

  if( (count = (INT4 *) LALCalloc( (size_t) (nmax + 1), sizeof(INT4))) == NULL ) {
    LALPrintError("Could not allocate Memory! \n");
    ABORT (stat, POLKAC_EMEM, POLKAC_MSGEMEM);
  }
  


  /* ------------------------------------------------------------- */
  /* Print out to the user-specified output file all the information in all the cell. 
     This file can be too huge to be tractable.*/
  if( (fp = fopen(CLA->OutputFile,"w")) == NULL ) 
    {
      LALPrintError("\n Cannot open file %s\n",CLA->OutputFile); 
      exit(POLKA_EXIT_ERR);
    }
  /* output for all the cells */
  print_info_of_the_cell( fp, cell, 0,(*ncell),0,0);
  fclose(fp);



  /* number counts and find the most significant event. */
  for(icell=0;icell<(*ncell);icell++) {
    nc=cell[icell].nCand;
    count[nc] += 1;
    if( Sigmax < cell[icell].significance) {
      Sigmax = cell[icell].significance;
      idxmax = icell;
    }
  }


  /* ------------------------------------------------------------- */
  /* output summary table. */
  if(lalDebugLevel < 3 ) {
    fprintf(stderr,"%% Maximly significant cell : freq [Hz]\tdec [rad]\tra [rad]  # [events]   Sig" "\n");
    fprintf(stderr, "%%\t\t\t     ");
    print_info_of_the_cell( stderr, cell, idxmax,idxmax+1,0,0);
    fprintf(stderr,"%% Maximly coincident cell  : freq [Hz]\tdec [rad]\tra [rad]  # [events]   Sig" "\n");
    fprintf(stderr, "%%\t\t\t     ");
    print_info_of_the_cell( stderr, cell, 0,1,0,0);

    nmax = cell[0].nCand;
    fprintf(stderr,"%% # of coincidences: \n");
    for(nc=0;nc<=nmax;nc++) {
      fprintf(stderr,"%7d",nc);
    }
    fprintf(stderr,"\n");
    fprintf(stderr,"%% # of cells       : \n");
    for(nc=0;nc<=nmax;nc++) { 
      fprintf(stdout,"%7d",count[nc]);
    }
    fprintf(stdout,"\n");
  }
  LALFree( count );


 
  if( CLA->AutoOut || cell[0].nCand >= CLA->Nthr ) 
    {
      if( (fpCoiCell = fopen(fnameCoiCell,"w")) == NULL || (fpCoiTime = fopen(fnameCoiTime,"w")) == NULL )
	{ 
	  LALPrintError("\n Cannot open file %s or %s\n",fnameCoiCell,fnameCoiTime); 
	  exit(POLKA_EXIT_ERR); 
	}
    }

  if( CLA->AutoOut || Sigmax > CLA->Sthr ) 
    {
      if( (fpSigCell = fopen(fnameSigCell,"w")) == NULL || (fpSigTime = fopen(fnameSigTime,"w")) == NULL )
	{ 
	  LALPrintError("\n Cannot open file %s or %s\n",fnameSigCell,fnameSigTime); 
	  exit(POLKA_EXIT_ERR); 
	}
    }



  /* ------------------------------------------------------------- */
  if( CLA->AutoOut ) 
    {  /* Output the info of the most significant and the most coincident event. */

      /* Output the info of the most coincident event. */
      /* Information of the cell. */
      print_info_of_the_cell( fpCoiCell, cell, 0,1,0,0);
      /* Print F stat from each file contributed to this cell. */
      print_Fstat_of_the_cell( fpCoiTime, cell, CList, 0,1,0,0 );

      /* Output the info of the most significant event. */
      /* Information of the cell. */
      print_info_of_the_cell( fpSigCell, cell, idxmax,idxmax+1,0,0);
      /* Print F stat from each file contributed to this cell. */
      print_Fstat_of_the_cell( fpSigTime, cell, CList, idxmax,idxmax+1,0,0 );

    } /* if( CLA->AutoOut ) */ 
  else 
    {
      /* output only on outliers larger than Nthr on number of coincidences and Sthr on significance.*/
      if( cell[0].nCand >= CLA->Nthr ) 
	{
	  
	  /* Information of the cell. */
	  print_info_of_the_cell( fpCoiCell, cell, 0, 0, 0, CLA->Nthr);
	  /* Print F stat from each file contributed to this cell. */
	  print_Fstat_of_the_cell( fpCoiTime, cell, CList, 0, 0, 0, CLA->Nthr );

	} /* if( cell[0].nCand > CLA->Nthr ) */
      
      
      /* ------------------------------------------------------------- */
      /* Second Sort arrays of candidates based on significance, if necessary. */ 
      /* output only on outliers */
      if( Sigmax > CLA->Sthr ) 
	{
	  qsort(cell, (size_t) (*ncell), sizeof(CellData), compareSignificance);

	  /* Information of the cell. */
	  print_info_of_the_cell( fpSigCell, cell, 0, 0, CLA->Sthr, 0);
	  /* Print F stat from each file contributed to this cell. */
	  print_Fstat_of_the_cell( fpSigTime, cell, CList, 0, 0, CLA->Sthr, 0 );

	} /* if( cell[0].significance > CLA->Sthr ) */
    } /* else of if( CLA->AutoOut ) */


  if( CLA->AutoOut || cell[0].nCand >= CLA->Nthr ) {
    fclose(fpCoiTime);
    fclose(fpCoiCell);
  }

  if( CLA->AutoOut || Sigmax > CLA->Sthr ) { 
    fclose(fpSigTime);
    fclose(fpSigCell);
  }


  RETURN (stat);
} /* PrintResult() */


/* ########################################################################################## */
/* Print_info_of_the_cell() */
void print_info_of_the_cell( FILE *fp, CellData *cd, INT4 icell_start, INT4 icell_end, REAL8 sig_thr, REAL8 ncand_thr)
{
  INT4 icell;

  icell = icell_start;
  while( icell < icell_end && 
	 cd[icell].significance > sig_thr && 
	 cd[icell].nCand > ncand_thr ) 
    {
      fprintf(fp,"%" LAL_REAL4_FORMAT "\t%" LAL_REAL4_FORMAT "\t%" LAL_REAL4_FORMAT "\t" 
	      "%" LAL_INT4_FORMAT "\t" "%" LAL_REAL4_FORMAT "\n",
	      cd[icell].Freq, cd[icell].Delta, cd[icell].Alpha,
	      cd[icell].nCand,
	      cd[icell].significance);
      icell++;
    }
}




/* ########################################################################################## */
/* Free memory */
void FreeMemory(LALStatus *stat, struct PolkaConfigVarsTag *CLA, CellData *cell, CandidateList *CList, UINT4 CLength)
{
  INITSTATUS( stat, "FreeMemory", rcsid );
  ATTATCHSTATUSPTR (stat);

  UINT4 icell;

  FreeConfigVars( stat->statusPtr, CLA );

  if( CList != NULL ) LALFree(CList);

  /* FIX (?) ME:  
     This part takes really long, when lalDebugLevel = 3. I do not know why.*/
  if( cell != NULL ) {
    for(icell=0;icell<CLength;icell++) {
      delete_int4_linked_list( cell[icell].CandID );
    }
    LALFree(cell);
  }
    
  DETATCHSTATUSPTR (stat);
  RETURN (stat);
} /* FreeMemory */


/* ########################################################################################## */
void FreeConfigVars(LALStatus *stat, struct PolkaConfigVarsTag *CLA )
{
  INITSTATUS( stat, "FreeConfigVars", rcsid );
  ATTATCHSTATUSPTR (stat);

  UINT4 k;

  if( CLA->FstatsFile != NULL ) LALFree(CLA->FstatsFile);
  if( CLA->OutputFile != NULL ) LALFree(CLA->OutputFile);
  if( CLA->InputDir != NULL ) LALFree(CLA->InputDir);
  if( CLA->BaseName != NULL ) LALFree(CLA->BaseName);


  if( (CLA->InputDir != NULL) && (CLA->BaseName != NULL) ) 
    { /* We have used glob and allocated mem for Filelist.*/
      for (k=0;k<CLA->NFiles;k++)
	{
	  if(CLA->Filelist[k] != NULL ) 
	    LALFree (CLA->Filelist[k]);
	} 
      LALFree (CLA->Filelist);
    }

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
} /* FreeCOnfigVars */


/* ########################################################################################## */
/* add data to linked structure */
void add_int4_data(LALStatus *stat, struct int4_linked_list **list_ptr, INT4 *data)
{
  INITSTATUS( stat, "add_int4_data", rcsid );

  struct int4_linked_list *p = NULL;
  p = (struct int4_linked_list *) LALMalloc(sizeof(struct int4_linked_list));
  if( p == NULL ) {
    LALPrintError("Could not allocate Memory! \n");
    ABORT (stat, POLKAC_EMEM, POLKAC_MSGEMEM);
  }
  p->data = *(data);
  p->next = *list_ptr;
  *list_ptr = p;

  RETURN (stat);
} /* void add_int4_data() */


/* ########################################################################################## */
/* delete data to linked structure */
void delete_int4_linked_list( struct int4_linked_list *list_ptr )
{
  INT4 ic;
  struct int4_linked_list *q;
  ic = 0;
  while( list_ptr !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) {  
    q = list_ptr->next;
    LALFree( list_ptr );
    list_ptr = q;
    ic++;
  }
  if( ic >  LINKEDSTR_MAX_DEPTH ) {
    LALPrintError("Maximum depth of linked structure reached!");
    exit(POLKA_EXIT_ERR);
  }
  return;
}

/* ########################################################################################## */
/* get info of this cell. */
void get_info_of_the_cell( CellData *cd, CandidateList *CList )
{
  INT4 idx, ic;
  REAL8 lfa;
  struct int4_linked_list *p;
  p = cd->CandID;

  ic = 0;
  while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) { 
    idx = p->data;
    lfa = CList[idx].TwoF/2.0 - log(1.0 + CList[idx].TwoF/2.0);
    cd->significance += lfa;
    cd->Alpha += CList[idx].Alpha;
    cd->Delta += CList[idx].Delta;
    cd->Freq += CList[idx].f;
    p = p->next;
    ic++;
  }

  if( ic >  LINKEDSTR_MAX_DEPTH ) {
    LALPrintError("Maximum depth of linked structure reached!");
    exit(POLKA_EXIT_ERR);
  }

  cd->Alpha /= cd->nCand;
  cd->Delta /= cd->nCand;
  cd->Freq  /= cd->nCand;
  return;
}



/* ########################################################################################## */
/* print F stat. */
void print_Fstat_of_the_cell( FILE *fp, CellData *cd, CandidateList *CList, INT4 icell_start, INT4 icell_end, REAL8 sig_thr, REAL8 ncand_thr )
{
  INT4 idx, ic, icell;
  struct int4_linked_list *p;


  icell = icell_start;
  while( icell < icell_end && 
	 cd[icell].significance > sig_thr && 
	 cd[icell].nCand > ncand_thr ) 
    {

      p = cd[icell].CandID;
      
      ic = 0;
      while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) { 
	idx = p->data;
	fprintf(fp,"%" LAL_INT4_FORMAT "\t%" LAL_INT4_FORMAT "\t%" LAL_REAL4_FORMAT "\n", 
		icell, CList[idx].FileID, CList[idx].TwoF );
	p = p->next;
	ic++;
      } /*   while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) {  */

      if( ic >  LINKEDSTR_MAX_DEPTH ) {
	LALPrintError("Maximum depth of linked structure reached!");
	exit(POLKA_EXIT_ERR);
      }

      icell++;
    } /*   while( icell < icell_end && ...  */

  return;
}


#if 0
/* ########################################################################################## */
/* print F stat. */
void print_info_of_the_cell( FILE *fp, CellData *cd, INT4 icell_start, INT4 icell_end, REAL8 sig_thr, REAL8 ncand_thr)
{
  INT4 idx, ic;
  struct int4_linked_list *p;
  p = cd->CandID;

  ic = 0;
  while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) { 
    idx = p->data;
    fprintf(fp,"%" LAL_INT4_FORMAT "\t%" LAL_INT4_FORMAT "\t%" LAL_REAL4_FORMAT "\n", 
	    outlier_id, CList[idx].FileID, CList[idx].TwoF );
    p = p->next;
    ic++;
  }

  if( ic >  LINKEDSTR_MAX_DEPTH ) {
    LALPrintError("Maximum depth of linked structure reached!");
    exit(POLKA_EXIT_ERR);
  }

  return;
}
#endif



/* ########################################################################################## */
/* Sorting function to sort candidate indices INCREASING order of f, delta, alpha, and 
   DECREASING ORDER OF F STATISTIC. */
int compareCandidates(const void *a, const void *b)
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
    F1=ip->TwoF;
    F2=jp->TwoF;
    /* I put F1 and F2 inversely, because I would like to get decreasingly-ordered set. */ 
    res = rfloatcompare( &F2,  &F1, 1);
  } 
  return res;
} /* int compareCandidates() */



/* ########################################################################################## */
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
  if( res == 0 ) {
    INT4 n1, n2;
    n1=ip->nCand;
    n2=jp->nCand;
    /* I put n1 and n2 inversely, because I would like to get decreasingly-ordered set. */ 
    res = rintcompare( &n2,  &n1, 1);
  } 


  return res;
} /* int compareSignificance() */



/* ########################################################################################## */
int compareNCandidate(const void *a, const void *b)
{
  const CellData *ip = a;
  const CellData *jp = b;
  int res;

  UINT4 n1, n2;
  n1=ip->nCand;
  n2=jp->nCand;
  /* I put n1 and n2 inversely, because I would like to get decreasingly-ordered set. */ 
  res = rintcompare( &n2,  &n1, 1);
  if( res == 0 ) {
    REAL8 F1, F2;
    F1=ip->significance;
    F2=jp->significance;
    /* I put F1 and F2 inversely, because I would like to get decreasingly-ordered set. */ 
    res = rfloatcompare( &F2,  &F1, 1);
  } 

  return res;
} /* int compareNCandidate() */



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


/* ########################################################################################## */
/* We would use glob in the future to read different files ? */
void ReadCandidateFiles(LALStatus *stat, CandidateList **CList, struct PolkaConfigVarsTag *CLA, UINT4 *clen)
{
  INITSTATUS( stat, "ReadCandidateFiles", rcsid );
  ATTATCHSTATUSPTR (stat);

  UINT4 kc;


  if( (CLA->InputDir != NULL) && (CLA->BaseName != NULL) ) 
    {
      TRY( GetFilesListInThisDir( stat->statusPtr, CLA->InputDir, CLA->BaseName, &(CLA->Filelist), &(CLA->NFiles) ), stat );
      
      *clen = 0;     /* We first have to set the candidate list length zero. */
      *CList = NULL; /* We first have to nullify the list. */
      for (kc=0;kc<CLA->NFiles;kc++)
	{
	  if( lalDebugLevel > 1 ) {
	    fprintf(stderr,"%s, %d\n",CLA->Filelist[kc],CLA->NFiles);
	  }
#ifdef USE_UNZIP
	  {INT4 FileID = 2*kc; /* the factor 2 because we have 2 sections in each file. */
	  ReadCandidateListFromZipFile( stat->statusPtr, CList, CLA->Filelist[kc], clen, &FileID);
	  }
#endif
	} 

    } /* if( (CLA->InputDir != NULL) && (CLA->BaseName != NULL) )  */
  else if ( CLA->FstatsFile != NULL ) 
    {
      TRY( ReadOneCandidateFile(stat->statusPtr, CList, CLA->FstatsFile, clen), stat );
      /* The last file is from last file.*/
      CLA->NFiles = (*CList)[*clen-1].FileID;
    } /* if( (CLA->InputDir != NULL) && (CLA->BaseName != NULL) )  */
  else 
    { /* We should not be here. */
      LALPrintError("\nYou have to specify either input data directory or input data file.\n");
      exit(POLKA_EXIT_ERR);;
    }

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
} /* ReadCandidateFiles() */



/* ########################################################################################## */
void GetFilesListInThisDir(LALStatus *stat, CHAR *directory, CHAR *basename, CHAR ***filelist, UINT4 *nfiles )
{
#ifdef HAVE_GLOB_H   
  CHAR command[512];
  UINT4 fileno=0;
  glob_t globbuf;
#endif

  INITSTATUS (stat, "GetFilesListInThisDir", rcsid);
  ATTATCHSTATUSPTR (stat);

#ifndef HAVE_GLOB_H   
  LALPrintError("Cannot use GetFilesListInThisDir() without glob.");
  ABORT( stat, POLKAC_EGLOB, POLKAC_MSGEGLOB);
#endif

  strcpy(command, directory);
  strcat(command,"/*");
  
  strcat(command, basename);
  strcat(command,"*");

#ifdef HAVE_GLOB_H
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);
  
  /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */
  
  if(globbuf.gl_pathc==0)
    {
      LALPrintError ("\nNo Input files in directory %s ... Exiting.\n\n", directory);
      ABORT (stat, POLKAC_ESYS, POLKAC_MSGESYS);
    }

  /* prepare memory for all filenames */
  *filelist = NULL;
  if ( ( *filelist = (CHAR**)LALCalloc(globbuf.gl_pathc, sizeof(CHAR*))) == NULL) {
    ABORT (stat, POLKAC_EMEM, POLKAC_MSGEMEM);
  }
  while ( fileno < (UINT4) globbuf.gl_pathc) 
    {
      (*filelist)[fileno] = NULL;
      if ( ((*filelist)[fileno] = (CHAR*)LALCalloc(1, strlen(globbuf.gl_pathv[fileno])+1)) == NULL) {
	ABORT (stat, POLKAC_EMEM, POLKAC_MSGEMEM);
      }
      strcpy((*filelist)[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
    }
  globfree(&globbuf);

  *nfiles = fileno; /* remember this is 1 more than the index value */
#endif

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
}




/* ########################################################################################## */
/* read and parse the given candidate 'Fstats'-file fname into the candidate-list CList */
#ifdef USE_UNZIP
/*
TODO:
Check if *CList is either NULL or the memory of which is previously allocated by alloc() or the kind.
(how?).
*/
void  
ReadCandidateListFromZipFile (LALStatus *stat, CandidateList **CList, CHAR *fname, UINT4 *candlen, INT4 *FileID)
{
  FILE *fp;
  const UINT4 max_num_candidates = 8000000; /* maximum tractable number of candidate events. */
  UINT4 numlines;
  INT4 nread;
  REAL8 epsilon=1e-5;
  UINT4 ic;
  INT4 length; /* length of file */
  CHAR *line, *endp; /* pointers to start and end of line */
  INT4 section = 0;    /* 0: non-POLKA, 1,2: IFO sections,
			 3: coincidence, 4: end-of-file */
  UINT4 nlines[2] = {0,0}; /* number of events for each IFO */
  const INT4 MAX_SECS = 4;

  UzpBuffer uzpbuff;

  INITSTATUS (stat, "ReadCandidateListFromZipFile", rcsid);
  ATTATCHSTATUSPTR (stat);


  /* check if file exists.*/
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      LALPrintError("File %s doesn't exist!\n",fname);
      ABORT( stat, POLKAC_ESYS, POLKAC_MSGESYS ); 
     }
  fclose(fp);

  /* Check if the memory to be allocated is not huge 
     (say, < 512MB. sizeof(CandidateList) ~ 60B. 512/60 = 8000000). */
  if( *candlen > max_num_candidates ) {
    LALPrintError("\nMaximum number of candidate events reached.\n");
    LALPrintError("\nWe have %u events while the maximum allowed number of events is %u.\n",*candlen,max_num_candidates);
    ABORT( stat, POLKAC_ESYS, POLKAC_MSGESYS ); 
  }



  uzpbuff.strptr = NULL;

  /* ------------------------------------------------------------------------- */
  /*  Open and count the size of the candidates file */
  /* Read into buffer.  If this fails, we can't proceed. */
  if ( getfile( &uzpbuff, fname )  < 0 ) {
    if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
    LALPrintError("Cannot read file %s . \n",fname);
    ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  length = uzpbuff.strlength;
  line = uzpbuff.strptr;
  if ( !line || length == 0 || *line == '\0' ) {
    if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
    LALPrintError ("Unknown format of the file  %s.\n\n", fname);
    ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* ------------------------------------------------------------------------- */
  /* Check for correct ending tag.  If this isn't present, it is
     safest not to proceed (it greatly simplifies error trapping). */
  line += length;
  if ( ( length < 8 || strncmp( line - 8, "\n%DONE\r\n", 8 ) ) &&
       ( length < 7 || strncmp( line - 7, "\n%DONE\n", 7 ) ) ) {
    free(uzpbuff.strptr);
    LALPrintError("File %s does not end with the DONE_MARKER. \n",fname);
    ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* ------------------------------------------------------------------------- */
  /* Start reading file data line-by-line and count the number of candidate events. */
  for ( line = uzpbuff.strptr; section < MAX_SECS;
	*endp = '\n', line = endp + 1 ) {

    /* Find end of line.  Previous endchecking assures us we will not
       fall off of the end of the file. */
    endp = line;
    while ( *endp != '\n' )
      endp++;

    /* Check for POLKA section divisions or EOF marker. */
    if ( !strncmp( line, "%1", 2 ) ) {
      section = 1;
      continue;
    } else if ( !strncmp( line, "%2", 2 ) ) {
      if( section != 1 ) { /* We should have section 1 before 2. */
	if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
	LALPrintError("Unknown format file %s.",fname);
	ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }
      section = 2;
      continue;
    } else if ( !strncmp( line, "%coincidence", 12 ) ) {
      if( section != 2 ) { /* We should have section 2 before 3. */
	if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
	LALPrintError("Unknown format file %s.",fname);
	ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }
      section = 3;
      break; /* We are not interested in the section 3 here. */
    }  /*   if ( !strncmp( line, "%1", 2 ) ) {*/
 

    /* Do non-POLKA checks: */
    if ( section == 0 ) 
      {
	LALPrintError("Unknown format file %s.",fname);
	if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
	ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }
    /* Do POLKA IFO-section checks */
    else if ( section == 1 || section == 2 ) 
      {
	nlines[section-1] += 1;
      }
    /* Do POLKA coincidence-section checks. */
    else 
      { /* we should not be here */
	if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
	LALPrintError("Unknown format file %s.",fname);
	ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      } /*     if ( section == 0 )  */


    /* Done reading this line. */
  } /*   for ( line = uzpbuff.strptr; section < MAX_SECS; ... */
  /* ------------------------------------------------------------------------- */

  numlines = nlines[0] + nlines[1]; /* we have two sections. */

  if( numlines == 0 ) { /* This file is empty. Go to the next file.*/
    if( lalDebugLevel > 1 ) {
      LALPrintError( "No candidate events in the file %s\n\n", fname);
    }
    free(uzpbuff.strptr);
    DETATCHSTATUSPTR (stat);
    RETURN (stat);
  } 

  /* ------------------------------------------------------------------------- */
  /* reserve memory for fstats-file contents */
  if ( numlines > 0) 
    { 
      CandidateList *tmp;
      tmp = (CandidateList *)LALRealloc (*CList, ( *candlen + numlines )*sizeof(CandidateList));
      if ( !tmp ) 
	{ 
	  if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
	  LALPrintError("Could not allocate memory for candidate file %s\n\n", fname);
	  ABORT (stat, POLKAC_EMEM, POLKAC_MSGEMEM);
	}
      *CList = tmp;
    }


  /* ------------------------------------------------------------------------- */
  /* Start reading file data line-by-line. */
  section = 0;
  ic = (*candlen);
  for ( line = uzpbuff.strptr; section < MAX_SECS;
	*endp = '\n', line = endp + 1 ) {


    /* Find end of line.  Previous endchecking assures us we will not
       fall off of the end of the file. */
    endp = line;
    while ( *endp != '\n' )
      endp++;
    *endp = '\0'; 
    /* Replace *endp = '\n' by '\0' makes it easy to read file line by line. */

    /* Check for POLKA section divisions or EOF marker. */
    if ( !strncmp( line, "%1", 2 ) ) {
      section = 1;
      continue;
    } else if ( !strncmp( line, "%2", 2 ) ) {
      section = 2;
      continue;
    } else if ( !strncmp( line, "%coincidence", 12 ) ) {
      section = 3;
      break; /* We are not interested in the section 3 here. */
    } 
 

    if ( section == 1 || section == 2 ) 
      {
	CandidateList *cl=&(*CList)[ic];
	ic++;

	nread = sscanf (line, 
			"%lf %lf %lf %lf", 
			&(cl->f), &(cl->Alpha), &(cl->Delta), &(cl->TwoF) );
	cl->FileID = (*FileID) + section - 1; /* section can be either 1 or 2. */


	/* ------------------------------------------------------------------------- */
	/* check that values that are read in are sensible */
	if (
	    cl->FileID < 0                        ||
	    cl->f < 0.0                        ||
	    cl->TwoF < 0.0                        ||
	    cl->Alpha <         0.0 - epsilon  ||
	    cl->Alpha >   LAL_TWOPI + epsilon  ||
	    cl->Delta < -0.5*LAL_PI - epsilon  ||
	    cl->Delta >  0.5*LAL_PI + epsilon  ||
	    !finite(cl->FileID)                     ||                                                                 
	    !finite(cl->f)                     ||
	    !finite(cl->Alpha)                 ||
	    !finite(cl->Delta)                 ||
	    !finite(cl->TwoF)
	    ) {
	  LALPrintError(
			"Line %d of file %s has invalid values.\n"
			"First 255 chars are:\n"
			"%s\n"
			"1st and 4th field should be positive.\n" 
			"2nd field should lie between 0 and %1.15f.\n" 
			"3rd field should lie between %1.15f and %1.15f.\n"
			"All fields should be finite\n",
			ic+1, fname, line, (double)LAL_TWOPI, (double)-LAL_PI/2.0, (double)LAL_PI/2.0);
	  LALFree ((*CList));
	  free( uzpbuff.strptr );
	  ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
	} /* end of the check of the range of the values.*/
      } /*     if ( section == 1 || section == 2 )  */
    /* Do POLKA coincidence-section checks. */
    else 
      { /* we should not be here */
	LALPrintError("Unknown format file %s.",fname);
	if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
	ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      } /* if ( section == 1 || section == 2 )  */
    
    /* Done reading this line and filling CList. */

           

    /* ------------------------------------------------------------------------- */
    /* check that we read 4 quantities with exactly the right format */
    if ( nread != 4 )
      {
	LALPrintError ("Found %d not %d values on line %d in file '%s'\n"
		       "Line in question is\n%s",
		       nread, 4, ic+1, fname, line);               
	LALFree ((*CList));
	free( uzpbuff.strptr );
	ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }

  } /*   for ( line = uzpbuff.strptr; section < MAX_SECS; ... ) */

  free( uzpbuff.strptr ); /* uzpbuff is allocated by malloc() in getfile(). It is user's responsibility to free this. */


  if (ic != (*candlen) + numlines ) {
    LALPrintError(
            "Read of file %s terminated after %d line but numlines=%d\n",
            fname, ic, numlines);
    LALFree((*CList));
    ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  if ( section != 3 ) {
    LALPrintError(
            "Read of file %s terminated not by coincidence section but %s\n",
            fname, line);
    LALFree((*CList));
    ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }


  (*candlen) += numlines; /* total number of candidate so far */


  DETATCHSTATUSPTR (stat);
  RETURN (stat);


} /* void  ReadCandidateListFromZipFile () */
#endif /* #ifdef USE_UNZIP */




/* ########################################################################################## */
/* read and parse the given candidate 'Fstats'-file fname into the candidate-list CList */
void  ReadOneCandidateFile (LALStatus *stat, CandidateList **CList, const CHAR *fname, UINT4 *candlen)
{
  UINT4 i;
  UINT4 numlines;
  REAL8 epsilon=1e-5;
  char line1[256];
  FILE *fp;
  INT4 nread;
  UINT4 checksum=0;
  UINT4 bytecount=0;


  INITSTATUS( stat, "ReadOneCandidateFile", rcsid );
  ATTATCHSTATUSPTR (stat);

  /* ------ Open and count candidates file ------ */
  i=0;
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      LALPrintError("File %s doesn't exist!\n",fname);
      ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
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
      ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
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
      ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
    }

  /* output a record of the running checksun amd byte count */
  LALPrintError( "%% %s: bytecount %" LAL_UINT4_FORMAT " checksum %" LAL_UINT4_FORMAT "\n", fname, bytecount, checksum);

  /* check validity of this Fstats-file */
  if ( strcmp(line1, DONE_MARKER ) ) 
    {
      LALPrintError ("ERROR: File '%s' is not properly terminated by: %sbut has %s instead", fname, DONE_MARKER, line1);
      ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
    }
  else
    numlines --;        /* avoid stepping on DONE-marker */

  *candlen=numlines;


  if (*candlen <= 0  )
    {
      LALPrintError("candidate length = %ud!\n",*candlen);
      exit(POLKA_EXIT_ERR);;
    }/* check that we have candidates. */


  
  /* reserve memory for fstats-file contents */
  if (numlines > 0) 
    { 
      *CList = (CandidateList *)LALMalloc (numlines*sizeof(CandidateList));
      if ( !CList ) 
        { 
          LALPrintError ("Could not allocate memory for candidate file %s\n\n", fname);
	  ABORT (stat, POLKAC_EMEM, POLKAC_MSGEMEM);
        }
    }

  /* ------ Open and count candidates file ------ */
  i=0;
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      LALPrintError("fopen(%s) failed!\n", fname);
      LALFree ((*CList));
      ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
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
	ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }
      
      nread = sscanf (line1, 
                     "%" LAL_INT4_FORMAT "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
                     "%c", 
                     &(cl->FileID), &(cl->f), &(cl->Alpha), &(cl->Delta), &(cl->TwoF), &newline );

      /* check that values that are read in are sensible */
      if (
          cl->FileID < 0                        ||
          cl->f < 0.0                        ||
          cl->TwoF < 0.0                        ||
          cl->Alpha <         0.0 - epsilon  ||
          cl->Alpha >   LAL_TWOPI + epsilon  ||
          cl->Delta < -0.5*LAL_PI - epsilon  ||
          cl->Delta >  0.5*LAL_PI + epsilon  ||
	  !finite(cl->FileID)                     ||                                                                 
          !finite(cl->f)                     ||
          !finite(cl->Alpha)                 ||
          !finite(cl->Delta)                 ||
          !finite(cl->TwoF)
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
	  ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
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
	ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }

      /* check that we read 6 quantities with exactly the right format */
      if ( nread != 6 )
        {
          LALPrintError ("Found %d not %d values on line %d in file '%s'\n"
                         "Line in question is\n%s",
                         nread, 6, i+1, fname, line1);               
          LALFree ((*CList));
          fclose(fp);
	  ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
        }



      i++;
    } /*  end of main while loop */
  /* check that we read ALL lines! */
  if (i != numlines) {
    LALPrintError(
            "Read of file %s terminated after %d line but numlines=%d\n",
            fname, i, numlines);
    LALFree((*CList));
    fclose(fp);
    ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* read final line with %DONE\n marker */
  if (!fgets(line1, sizeof(line1), fp)) {
    LALPrintError(
            "Failed to find marker line of file %s\n",
            fname);
    LALFree((*CList));
    fclose(fp);
    ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* check for %DONE\n marker */
  if (strcmp(line1, DONE_MARKER)) {
    LALPrintError(
            "Failed to parse marker: 'final' line of file %s contained %s not %s",
            fname, line1, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* check that we are now at the end-of-file */
  if (fgetc(fp) != EOF) {
    LALPrintError(
            "File %s did not terminate after %s",
            fname, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    ABORT (stat, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* -- close candidate file -- */
  fclose(fp);     

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* ReadOneCandidateFile() */


/* ########################################################################################## */
void ReadCommandLineArgs(LALStatus *stat, int argc,char *argv[], struct PolkaConfigVarsTag *CLA) 
{
  INITSTATUS( stat, "ReadCommandLineArgs", rcsid );
  ATTATCHSTATUSPTR (stat);

  ASSERT( CLA != NULL, stat, POLKAC_ENULL, POLKAC_MSGENULL);


  CHAR* uvar_InputData;
  CHAR* uvar_OutputData;

  CHAR* uvar_InputDirectory;
  CHAR* uvar_BaseName;

  BOOLEAN uvar_AutoOut;
  INT4 uvar_Nthr;      
  REAL8 uvar_Sthr;      

  REAL8 uvar_FreqWindow;
  REAL8 uvar_AlphaWindow;
  REAL8 uvar_DeltaWindow;
  REAL8 uvar_FreqShift;
  REAL8 uvar_AlphaShift;
  REAL8 uvar_DeltaShift;
  BOOLEAN uvar_help;


  const CHAR BNAME[] = "Test";

  uvar_AutoOut = 0;
  uvar_help = 0;

  uvar_InputData = NULL;
  uvar_OutputData = NULL;

  uvar_InputDirectory = NULL;
  uvar_BaseName = (CHAR*)LALCalloc (1, strlen(BNAME)+1);
  strcpy (uvar_BaseName, BNAME);

  /* The following numbers are arbitrary. */
  uvar_Nthr = 65536;     
  uvar_Sthr = 1.0e5;      

  uvar_FreqWindow = 0.0;
  uvar_AlphaWindow = 0.0;
  uvar_DeltaWindow = 0.0;

  uvar_FreqShift = 0.0;
  uvar_AlphaShift = 0.0;
  uvar_DeltaShift = 0.0;



  /* register all our user-variables */
  LALregBOOLUserVar(stat,       help,           'h', UVAR_HELP,     "Print this message"); 

  LALregSTRINGUserVar(stat,     OutputData,     'o', UVAR_REQUIRED, "Ouput candidates file name");

  LALregSTRINGUserVar(stat,     InputData,      'I', UVAR_OPTIONAL, "Input candidates Fstats file.");
  LALregSTRINGUserVar(stat,     InputDirectory, 'i', UVAR_OPTIONAL,"Input candidates Fstats files directory.");
  LALregSTRINGUserVar(stat,     BaseName,       'b', UVAR_OPTIONAL,"BaseName of the Input Fstats files");

  LALregINTUserVar(stat,        Nthr,            0,  UVAR_OPTIONAL, "Threshold on number of coincidence");
  LALregREALUserVar(stat,       Sthr,            0,  UVAR_OPTIONAL, "Threshold on significance.");
  LALregBOOLUserVar(stat,       AutoOut,         0,  UVAR_OPTIONAL, "Set Nthr and Sthr to print most significant cell only."); 

  LALregREALUserVar(stat,       FreqWindow,     'f', UVAR_REQUIRED, "Frequency window in Hz");
  LALregREALUserVar(stat,       AlphaWindow,    'a', UVAR_REQUIRED, "Right ascension window in radians");
  LALregREALUserVar(stat,       DeltaWindow,    'd', UVAR_REQUIRED, "Declination window in radians");
  LALregREALUserVar(stat,       FreqShift,      'F', UVAR_OPTIONAL, "Frequency shift in FreqWindow");
  LALregREALUserVar(stat,       AlphaShift,     'A', UVAR_OPTIONAL, "Right ascension shift in AlphaWindow");
  LALregREALUserVar(stat,       DeltaShift,     'D', UVAR_OPTIONAL, "Declination shift in DeltaWindow");

  TRY (LALUserVarReadAllInput(stat->statusPtr,argc,argv),stat); 


  if (uvar_help) {	/* if help was requested, we're done here */
    LALPrintError("%s\n",rcsid);
    fflush(stderr);
    LALDestroyUserVars(stat->statusPtr);
    exit(POLKA_EXIT_OK);
  }


  if( LALUserVarWasSet (&uvar_InputData) && 
      LALUserVarWasSet (&uvar_InputDirectory) ) {
    LALPrintError("\nCannot set both of InputData and InputDirectory\n");
    exit(POLKA_EXIT_OK);
  }

  if( (!LALUserVarWasSet (&uvar_InputData)) && 
      (!LALUserVarWasSet (&uvar_InputDirectory)) ) {
    LALPrintError("\nPlease set either InputData and InputDirectory\n");
    exit(POLKA_EXIT_OK);
  }


#ifndef USE_UNZIP
  if( LALUserVarWasSet (&uvar_InputDirectory) ) {
    LALPrintError("\nCannot use -i option without enabling unzip.\n");
    ABORT( stat, POLKAC_EUNZIP, POLKAC_MSGEUNZIP );
  }
#endif


  CLA->FstatsFile = NULL;
  CLA->OutputFile = NULL;
  CLA->InputDir = NULL;
  CLA->BaseName = NULL;


  if( LALUserVarWasSet (&uvar_InputData) ) {
    CLA->FstatsFile = (CHAR *) LALMalloc(strlen(uvar_InputData)+1);
    if(CLA->FstatsFile == NULL)
      {
	LALPrintError("No candidates file specified; input with -I option.\n");
	LALPrintError("For help type %s -h\n", argv[0]);
	RETURN (stat);
      }      
    strcpy(CLA->FstatsFile,uvar_InputData);
  }

  CLA->OutputFile = (CHAR *) LALMalloc(strlen(uvar_OutputData)+1);
  if(CLA->OutputFile == NULL)
    {
      FreeConfigVars( stat->statusPtr, CLA );
      ABORT (stat, POLKAC_EMEM, POLKAC_MSGEMEM);
    }      

  strcpy(CLA->OutputFile,uvar_OutputData);

  if( LALUserVarWasSet (&uvar_InputDirectory) ) {
#ifndef HAVE_GLOB_H   
    LALPrintError("Sorry, but you cannot use this feature without glob.h.\n");
    exit(POLKA_EXIT_ERR);
#endif
    CLA->InputDir = (CHAR *) LALMalloc(strlen(uvar_InputDirectory)+1);
    if(CLA->InputDir == NULL)
      {
	FreeConfigVars( stat->statusPtr, CLA );
	ABORT (stat, POLKAC_EMEM, POLKAC_MSGEMEM);
      }          
    strcpy(CLA->InputDir,uvar_InputDirectory);
  }



  CLA->BaseName = (CHAR *) LALMalloc(strlen(uvar_BaseName)+1);
  if(CLA->BaseName == NULL)
    {
      FreeConfigVars( stat->statusPtr, CLA );
      RETURN (stat);
    }          
  strcpy(CLA->BaseName,uvar_BaseName);


  CLA->AutoOut = uvar_AutoOut;
  CLA->Nthr = uvar_Nthr;
  CLA->Sthr = uvar_Sthr;


  CLA->Deltaf = uvar_FreqWindow;
  CLA->DeltaAlpha = uvar_AlphaWindow;
  CLA->DeltaDelta = uvar_DeltaWindow;

  CLA->Shiftf = uvar_FreqShift;
  CLA->ShiftAlpha = uvar_AlphaShift;
  CLA->ShiftDelta = uvar_DeltaShift;

  LALDestroyUserVars(stat->statusPtr);
  BEGINFAIL(stat) {
    LALFree(CLA->FstatsFile);
    LALFree(CLA->OutputFile);
  } ENDFAIL(stat);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
} /* void ReadCommandLineArgs()  */

