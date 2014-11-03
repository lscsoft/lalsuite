/*
*  Copyright (C) 2007 Badri Krishnan,  Holger Pletsch, Reinhard Prix, Yousuke Itoh
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

/********************************************************************************************/
/*      zellepolka - the pulsar coincidence analysis code for einstein at home postprocess  */
/*                                                                                          */
/*      Xavier Siemens,  Bruce Allen,  Bernd Machenschalk,  Yousuke Itoh,  Holger Pletsch   */
/*                                                                                          */
/*                                                                                          */
/*                                  UWM - April  2005                                       */
/*                             Based on uberpolka written by                                */
/*                        Xavier Siemens, Bruce Allen, Bernd Machenschalk                   */
/********************************************************************************************/

/*! 
   @file
   @author Xavier Siemens,  Bruce Allen,  Bernd Machenschalk,  Yousuke Itoh , Holger Pletsch
   @brief The pulsar coincidence analysis code for einstein at home post-processing --- counting
    number of events in cells construcetd in 4D parameters space.

    \par Inputs and outputs
<ul>
<li> This code takes one single file generated from the EaH zipped result files by a python code combiner_v4.py,
to be found at:
http://www.lsc-group.phys.uwm.edu/cgi-bin/cvs/viewcvs.cgi/einsteinathome/CFS/post_processing/combiner_v4.py?cvsroot=lscsoft

<li> This code outputs one output file whose name is specified by the user, five files with default names, and 
some information on stdout.
<ol>
<li> The one output file contains all the information (f,a,d,f1dot,ncandidate,sig) in all the cells. This can be huge.
<li> One file outputs time variation of 2F of some significant outliers.
<li> One file outputs cell information of some significant outliers
<li> One file outputs time variation of 2F of some coincident outliers.
<li> One file outputs cell information of some coincident outliers
<li> One file outputs the cell info of the maximum coincident event over each Frequency cell but all over the sky and spin-down.
<li> Outputs summary table on stdout. This table shows how many cells have how many counts. 
     (e.g. 2 cells have 10 coincidences.)
<li> Outputs most coincident cell with its events.
</ol>
</ul>

\par Algorithm
<ol>
<li> First construct a grid in four dimensional parameters space (frequency, right ascension, declination and f1dot).
   A varying cell width in declination is used according to the metric-grid used in the E\@H S4 search
   (the actual implementation uses a Gaussian declination-model, with maximum around the equator),
   a uniform cell-grid in the frequency and spin-down.
   We can change the cell-grid spacing of the each parameter, and shift the cell-grid as a whole for 
  each parameter.  <br> 
<li>  Then we count how many candidate events are found in each cell.
  Even if a file has more than one event in a cell, we say we have one event from that file, 
  and take the largest F statistic event among the events in the cell from that file.
  Therefore, the greatest number counts in a cell is the total number of data-stretches from
  which the result files may have originated.
</ol>

*/



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

#define ADDITIONAL_MEM 32768

/* 
This can be defined for testing purposes only. So it should not be defined 
If defined, it does:
 + output entire information about ALL cells
 + output candidate events in cells which have less than 17 coincidences
 + randomly set candidate events' 2F-values

#define SKYTESTMODE  
*/

/*                                                                                                                                                                      
   To use unzip, you need to have unzip-5.5x from, say,  the InfoZip webpage,                                                                                                       
   and readzipfile_util.h and .c. from yousuke.                                                                                                                  
   ( not needed any more in S4 , will be done by combiner_v4.py )                                                                                                                   
#define USE_UNZIP                                                                                                                                                                     
*/

#define EPSEDGE 1e-6

/* ----------------------------------------------------------------------------- */
/* file includes */
#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "getopt.h"
#include <math.h>
#include <unistd.h>


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
#include <locale.h>

/* this is defined in C99 and *should* be in math.h.  Long term
   protect this with a HAVE_FINITE */
#ifdef _MSC_VER
#include <float.h>
#define finite _finite
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
#define POLKAC_EGLOB            8

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
/*!
Configuration variables have three categories. 

Cell parameters: Define dimensions of cells and center of the cells.
@param Deltaf     REAL8 Size of coincidence window in Hz (= Size of cell in frequency)
@param DeltaAlpha REAL8 Size of coincidence window in radians (at the equator delta=0)
@param DeltaDelta REAL8 Size of coincidence window in radians (at the equator delta=0)
@param DeltaF1dot REAL8 Size of coincidence window of spindown d/dt f (= Size of cell in f1dot)
@param Shiftf     REAL8 Parallel shift of frequency in Hz of cell 
@param ShiftAlpha REAL8 Parallel shift of Alpha in Hz of cell
@param ShiftDelta REAL8 Parallel shift of Delta in Hz of cell 
@param ShiftF1dot REAL8 Parallel shift of F1dot of cell
@param Kappa      REAL8 Tuning parameter for declination window

Input: Control the way to read a file or files in a directory
@param *FstatsFile CHAR Names of Fstat files to be read in.
@param *InputDir   CHAR Directory name of input files 
@param *BaseName   CHAR Base name of input files 
@param **Filelist  CHAR Array of filenames to load Fstats file from 
@param NFiles;     UINT4 Number of input files read
@param TwoFthr;    REAL8 Threshold for 2F values

Output: Control the way to return results
@param OutputFile  CHAR Names of output file
@param Nthr        INT4 Show exective results of cells with numbers of coincidence above Nthr.
@param Sthr        REAL4 Show exective results of cells with significance above Sthr.
@param AutoOut     BOOLEAN If set, output the info of the most significant and the most coincident event.
*/

typedef struct PolkaConfigVarsTag 
{
  REAL8 TwoFthr;     /*  Threshold for TwoF values */
  REAL8 Deltaf;      /*  Size of coincidence window in Hz */
  REAL8 DeltaF1dot;  /*  Size of coincidence window of spindown */
  REAL8 DeltaAlpha;  /*  Size of coincidence window in radians (at equator) */
  REAL8 DeltaDelta;  /*  Size of coincidence window in radians (at equator) */
  REAL8 Kappa;       /*  Tuning parameter for declination window */
  REAL8 Shiftf;      /*  Parallel shift of frequency of cell */
  REAL8 ShiftAlpha;  /*  Parallel shift of Alpha of cell */
  REAL8 ShiftDelta;  /*  Parallel shift of Delta of cell */
  REAL8 ShiftF1dot;  /*  Parallel shift of F1dot spindown of cell */
  CHAR *FstatsFile;  /*  Names of Fstat files to be read in */
  CHAR *OutputFile;  /*  Names of output file */
  CHAR *InputDir;    /*  Directory name of input files */
  CHAR *BaseName;    /*  Base name of input files */
  CHAR *EahRun;      /*  E\@H identifying run label */
  CHAR **Filelist;   /*  Array of filenames to load Fstats file from */
  UINT4 NFiles;      /*  Number of input files read */
  INT4 Nthr;         /*  Show exective results of cells with numbers of coincidence above Nthr. */
  INT4 CellGrid;
  INT4 MinCellCoin;   /*  Output cells having Coin's > MinCellCoin into separate files. */
  REAL4 Sthr;        /*  Show exective results of cells with significance above Sthr. */
  BOOLEAN AutoOut;
  BOOLEAN UseUnzip;
} PolkaConfigVars;


/* This structure contains the indices corresponding to the 
coarse frequency and sky bins */
/*!
CandidateList

@param iCand         Candidate id: unique within this program.  
@param f             Frequency of the candidate 
@param Alpha         right ascension of the candidate
@param Delta         declination  of the candidate 
@param F1dot         spindown (d/dt f) of the candidate
@param TwoF          2F of this candidate event
@param FileID        File ID to specify from which file the candidate under consideration originaly comes. 
@param iFreq         Frequency index of this candidate event
@param iDelta        Declination index. This can be negative. 
@param iAlpha        Right ascension index of this candidate event
@param iF1dot        Spindwon index of this candidate event
*/
#pragma pack(1) 
typedef struct CandidateListTag
{
  REAL8 f;           /*  Frequency of the candidate */
  REAL8 F1dot;       /*  spindown (d/dt f) of the candidate */
  REAL4 Alpha;       /*  right ascension of the candidate */
  REAL4 Delta;       /*  declination  of the candidate */
  REAL4 TwoF;        /*  Maximum value of F for the cluster */
  INT4 FileID;       /*  File ID to specify from which file the candidate under consideration originaly comes. */
  INT4 iFreq;        /*  Frequency index */
  INT4 iDelta;       /*  Declination index. This can be negative. */
  INT4 iAlpha;       /*  Right ascension index */
  INT4 iF1dot;       /*  Spindown index */
  INT8 iCand;       /*  Candidate id: unique with in this program.  */
} CandidateList;     /*   Fstat lines */ 
#pragma pack(0)


/*!
Liked list containing one INT8 data..

@param data INT8
@param next int8_linked_list*
*/
struct int8_linked_list {
  INT8 data;
  struct int8_linked_list *next;
}; 


/*!
Structure containg data in cells

@param Freq          REAL8 Frequency index of the cell 
@param Alpha         REAL4 Right ascension index of the cell
@param Delta         REAL4 Declination index of the cell
@param F1dot         REAL8 Spindown index of the cell 
@param iFreq         INT4 Frequency index of this candidate event
@param iDelta        INT4 Declination index. This can be negative. 
@param iAlpha        INT4 Right ascension index of this candidate event
@param iF1dot        INT4 Spindown index of this candidate event
@param significance  REAL8 minus log of joint false alarm of the candidates in this cell
@param nCand;        INT4 number of the events in this cell
@param CandID        int8_linked_list* linked structure that has candidate id-s of the candidates in this cell
*/

#pragma pack(1) 
typedef struct CellDataTag
{
  REAL8 Freq;          /*  Frequency index of the cell */
  REAL8 F1dot;          /*  Spindown index of the cell */
  REAL4 Alpha;         /*  Right ascension index of the cell */
  REAL4 Delta;         /*  Declination index of the cell */
  REAL4 significance;  /*  minus log of joint false alarm of the candidates in this cell. */
  INT4 iFreq;          /*  Frequency index of this candidate event */
  INT4 iDelta;         /*  Declination index of this candidate event */
  INT4 iAlpha;         /*  Right ascension index of this candidate event */
  INT4 iF1dot;         /*  Spindown index of this candidate event */
  INT4 nCand;          /*  number of the events in this cell. */
  struct int8_linked_list *CandID;  /* linked structure that has candidate id-s of the candidates in this cell. */
} CellData;
#pragma pack(0)


/* ----------------------------------------------------------------------------- */
/* Function declarelations */
void ReadCommandLineArgs( LALStatus *, INT4 argc, CHAR *argv[], PolkaConfigVars *CLA ); 
void GetFilesListInThisDir(LALStatus *, const CHAR *directory, const CHAR *basename, CHAR ***filelist, UINT4 *nfiles );
void ReadCandidateFiles( LALStatus *, CandidateList **Clist, PolkaConfigVars *CLA, INT8 *datalen );
void ReadOneCandidateFile( LALStatus *, CandidateList **CList, const CHAR *fname, INT8 *datalen, const REAL8 myFthr );
void ReadOneCandidateFileV2( LALStatus *lalStatus, CandidateList **CList, const CHAR *fname, INT8 *candlen );

#ifdef USE_UNZIP
void ReadCandidateListFromZipFile (LALStatus *, CandidateList **CList, CHAR *fname, INT8 *candlen, const INT4 *FileID);
#endif

void RePrepareCells( LALStatus *, CellData **cell, const INT8 CLength , const INT8 iposition);
void PrepareCells( LALStatus *, CellData **cell, const INT8 datalen );

void add_int8_data(LALStatus *, struct int8_linked_list **list_ptr, const INT8 *data);
void delete_int8_linked_list( LALStatus *, struct int8_linked_list *list_ptr);

void get_info_of_the_cell( LALStatus *, CellData *cd, const CandidateList *CList);

int compareREAL8arrays(const REAL8 *ap, const REAL8 *bp, size_t n);
int compareINT4arrays(const INT4 *ap, const INT4 *bp, size_t n);
int compareINT8arrays(const INT8 *ap, const INT8 *bp, size_t n);
int compareSignificances(const void *a, const void *b);

void PrintResult( LALStatus *, const PolkaConfigVars *CLA, CellData *cell, const INT8 *ncell, 
		  CandidateList *CList, const INT4 cellgridnum, INT8 CellListi[]);

void print_Fstat_of_the_cell( LALStatus *, FILE *fp, const CellData *cd, const CandidateList *CList, const INT8 icell_start, 
			      const INT8 icell_end, const REAL8 sig_thr, const REAL8 ncand_thr );
void print_info_of_the_cell( LALStatus *lalStatus, FILE *fp, const CellData *cd, const INT8 icell_start, 
			     const INT8 icell_end, const REAL8 sig_thr, const REAL8 ncand_thr);
void print_info_of_cell_and_ifo_S4R2a( LALStatus *, FILE *fp, const CellData *cd, const CandidateList *CList, const INT8 icell_start,
				 const INT8 icell_end, const REAL8 sig_thr, const REAL8 ncand_thr );
void print_info_of_cell_and_ifo_S5R1a( LALStatus *, FILE *fp, const CellData *cd, const CandidateList *CList, const INT8 icell_start,
				       const INT8 icell_end, const REAL8 sig_thr, const REAL8 ncand_thr );

void print_cand_of_most_coin_cell( LALStatus *lalStatus, CellData *cd, const CandidateList *CList);

void FreeMemory(LALStatus *, PolkaConfigVars *CLA, CandidateList *CList);
void FreeMemoryCellsOnly(LALStatus *, CellData *cell, const INT8 datalen);
void FreeConfigVars(LALStatus *, PolkaConfigVars *CLA );

void sortCandidates(INT8 *data, INT8 N);
void sortCandidates2(INT8 *data, INT8 left, INT8 right);
void sortFreqCells2(INT8 *data, INT8 left, INT8 right);


/* ----------------------------------------------------------------------------- */
/* Global Variables */

CandidateList *SortedC = NULL; /* need to be global to access them in compare functions for qsort */
CellData *global_cell = NULL;

/*
ListForSort *SortedCList = NULL;
ListForSort *cellList = NULL;
*/

/*! @param global_status LALStatus Used to initialize LALStatus lalStatus. */
LALStatus global_status;
/*! @param lalDebugLevel INT4 Control debugging behaviours. Defined in lalapps.h */
/*! @param vrbflg        INT4 Control debugging messages. Defined in lalapps.h */
extern INT4 vrbflg;

/* ------------------------------------------------------------------------------------------*/
/* Code starts here.                                                                         */
/* ------------------------------------------------------------------------------------------*/
/* ########################################################################################## */
/*!
  Main function

  @param[in] argc   INT4
  @param[in] argv[] CHAR*
  @return    return 0 on normal exit.  
*/
int main(INT4 argc,CHAR *argv[]) 
{
  PolkaConfigVars PCV; /* Configuration variables */
  REAL8 *LookupDelta1 = NULL;
  REAL8 *LookupDelta2 = NULL;
  REAL8 DeltaWin=0;
  REAL8 DeltaBorder=0;
  REAL8 d0=0,a0=0,k0=0;
  INT4  NumDeltaWins, idb;
  INT4  iDeltaMax1=0;
  INT4  iDeltaMax2=0;
  INT4  DeltaDeltaStep=0;
  INT8  CLength=0; /* Length of the candidate-events list */
  INT8  *SortedCListi = NULL; /* List of candidate-events INT-labels */
  INT8  *cellListi = NULL;  /* List of cell INT-labels */
  INT8  icell, icand, ncell;
  INT4  cc1,cc2,cc3,cc4,bb1,bb2,bb3,bb4; /* cell-grid shifting variables */
  INT4  selectGrid; /* denotes one of the 16 shifted cell-grid */
  INT8  sizecells; /* Length of the cell list */
  
  LALStatus *lalStatus = &global_status;
  setlocale(LC_ALL, "C");
  vrbflg = 1;   /* verbose error-messages */

#ifdef SKYTESTMODE
  /* Set the seed */
  FILE *devrandom = NULL;
  INT4 errorcode = 0;
  INT4 seed = 0;
  if (!(devrandom=fopen("/dev/urandom","r")))  {
    fprintf(stderr,"Unable to open device /dev/urandom\n");
  }
  errorcode=fread((void*)&seed,sizeof(INT4),1,devrandom);
  if (errorcode!=1)  {
    fprintf(stderr, "Error reading /dev/random file!\n");
  }
  fclose(devrandom);
  srand( seed );
#endif

  /* Reads command line arguments */
  LAL_CALL( ReadCommandLineArgs( lalStatus, argc,argv, &PCV ), lalStatus); 

  /* Reads in candidare files, set CLength */
  LAL_CALL( ReadCandidateFiles(lalStatus, &SortedC, &PCV, &CLength), lalStatus);

  /* More efficent sorting in memory by using an integer list */
  SortedCListi = (INT8 *) LALCalloc(CLength, sizeof(INT8) );
  if( SortedCListi == NULL ) {
    return POLKAC_EMEM;
  }


  /* --------------------------------------------------------------------------------*/      
  
  if( !strcmp(PCV.EahRun,"S4R2a") || !strcmp(PCV.EahRun,"S5R1a") ) {

    d0=PCV.DeltaDelta;
    a0=PCV.DeltaAlpha;
    k0=PCV.Kappa;
  
    /* Construct the 1st lookup table for the declination coincidence windows */
    DeltaWin = d0;
    DeltaBorder = DeltaWin / 2; /* half width of first cell */
    NumDeltaWins = 1;
    /* Count cells in declination direction */
    while ( DeltaBorder < (LAL_PI*0.5) ) {
      DeltaWin = (a0 / (fabs(sin(fabs(DeltaBorder) - (k0 * a0)))));
      if (DeltaWin > d0) {
	DeltaWin = d0;
      }
      DeltaBorder = DeltaBorder + DeltaWin;
      NumDeltaWins++;
    }
    /* Allocate memory for the declination cell lookup table */
    LookupDelta1 = (REAL8 *) LALCalloc(NumDeltaWins+1, sizeof(REAL8) );
    if( LookupDelta1 == NULL ) {
      return POLKAC_EMEM;
    }
    /* Finally calclate the 1st declination lookup table */
    idb=0;
    DeltaWin = d0;
    DeltaBorder = DeltaWin / 2;
    LookupDelta1[idb] = DeltaBorder;
    while ( DeltaBorder < (LAL_PI*0.5) ) {
      idb++;
      DeltaWin = (a0 / (fabs(sin(fabs(DeltaBorder) - (k0 * a0)))));
      if (DeltaWin > d0) {
	DeltaWin = d0;
      }
      DeltaBorder = DeltaBorder + DeltaWin;
      LookupDelta1[idb] = DeltaBorder;
    }
    iDeltaMax1 = idb;
    
    
    /* Construct the 2nd lookup table for the declination coincidence windows */
    DeltaWin = d0;
    DeltaBorder = DeltaWin; /* split at equator */
    NumDeltaWins = 1;
    /* Count cells in declination direction */
    while ( DeltaBorder < (LAL_PI*0.5) ) {
      DeltaWin = (a0 / (fabs(sin(fabs(DeltaBorder) - (k0 * a0)))));
      if (DeltaWin > d0) {
	DeltaWin = d0;
      }
      DeltaBorder = DeltaBorder + DeltaWin;
      NumDeltaWins++;
    }
    /* Allocate memory for the declination cell 2nd lookup table */
    LookupDelta2 = (REAL8 *) LALCalloc(NumDeltaWins+1, sizeof(REAL8) );
    if( LookupDelta2 == NULL ) {
      return POLKAC_EMEM;
    }
    /* Finally calclate the 2nd declination lookup table */
    idb=0;
    DeltaWin = d0;
    DeltaBorder = DeltaWin;
    LookupDelta2[idb] = DeltaBorder;
    while ( DeltaBorder < (LAL_PI*0.5) ) {
      idb++;
      DeltaWin = (a0 / (fabs(sin(fabs(DeltaBorder) - (k0 * a0)))));
      if (DeltaWin > d0) {
	DeltaWin = d0;
      }
      DeltaBorder = DeltaBorder + DeltaWin;
      LookupDelta2[idb] = DeltaBorder;
    }
    iDeltaMax2 = idb;
  
  } /* if( !strcmp(CLA->EahRun,"S4R2a") || !strcmp(CLA->EahRun,"S5R1a") )  */
    

  /*--------------------------------------------------------------------------*/
  /* Loops for the 2^4 = 16 cell-grid shifts along the 4 different dimensions */
  selectGrid = 0;
  bb1=0; bb2=0; bb3=0; bb4=0;
  cc1=0; cc2=0; cc3=0; cc4=0;
  
  for (bb1=0; bb1<2; bb1++) {
    for (bb2=0; bb2<2; bb2++) {
      for (bb3=0; bb3<2; bb3++) {
	for (bb4=0; bb4<2; bb4++) {

#ifdef SKYTESTMODE
	  if (bb1 ==0 && bb4 == 0) {
#endif
	  /*if (selectGrid == PCV.CellGrid) {*/
	    cc1=bb1; /* shift in f */
	    cc2=bb2; /* shift in alpha */
	    cc3=bb3; /* shift in delta */
	    cc4=bb4; /* shift in f1dot */
	  /*}*/
	    
	    fprintf(stdout,"\n%% Selected CellGrid: %d %d %d %d\n", cc1,cc2,cc3,cc4);

	     /* Prepare cells. */
	    global_cell = NULL;
	    sizecells = ADDITIONAL_MEM;
	    LAL_CALL( PrepareCells( lalStatus, &global_cell, sizecells ), lalStatus);  

	    /* Assigning four indices to each candidate event */
	    for (icand=0;icand<CLength;icand++) 
	      {
		
		/* Assign the FREQUENCY index to the candidate event */
                SortedC[icand].iFreq = floor( ((SortedC[icand].f) / (PCV.Deltaf)) + (cc1 * 0.5)  );

		/* Assign the SKY indices to the candidate event */
	
		if( !strcmp(PCV.EahRun,"S4R2a") || !strcmp(PCV.EahRun,"S5R1a") ) 
		  {
		    /* Check if Alpha has a correct value and reset if necessary */
		    if( SortedC[icand].Alpha < 0.0 ) {
		      SortedC[icand].Alpha = SortedC[icand].Alpha + LAL_TWOPI;
		    }
		    if( SortedC[icand].Alpha > LAL_TWOPI ) {
		      SortedC[icand].Alpha = SortedC[icand].Alpha - LAL_TWOPI;
		    }
		    
		    /* Assign the ALPHA index to the candidate event */
		    SortedC[icand].iAlpha = floor( ((SortedC[icand].Alpha) * (cos(SortedC[icand].Delta)) / (PCV.DeltaAlpha))  + (cc2 * 0.5) );
		  
		    /* Assign the DELTA index to the candidate event */
		    DeltaDeltaStep=0;
		    if (cc3 == 0) { /* unshifted cell-grid */
		      while ( LookupDelta1[DeltaDeltaStep] < fabs(SortedC[icand].Delta) ) {
			DeltaDeltaStep++;
			if (DeltaDeltaStep >= iDeltaMax1)
			  break;
		      }
		      if ( SortedC[icand].Delta < 0 ) {
			SortedC[icand].iDelta = -(INT4)(2 * DeltaDeltaStep);
		      }
		      else {
			SortedC[icand].iDelta = (INT4)(2 * DeltaDeltaStep);
		      }
		    }
		    else { /* shifted cell-grid */
		      while ( LookupDelta2[DeltaDeltaStep] < fabs(SortedC[icand].Delta) ) { 
			DeltaDeltaStep++;
			if (DeltaDeltaStep >= iDeltaMax2)
			  break;
		      }
		      if ( SortedC[icand].Delta < 0 ) {
			SortedC[icand].iDelta = -(INT4)((2 * DeltaDeltaStep) + 1);
		      }
		      else {
			SortedC[icand].iDelta = (INT4)((2 * DeltaDeltaStep) + 1);
		      }
		    }
		  }
		    
		else 
		  {
		    /* Assign the ALPHA index to the candidate event */
		    SortedC[icand].iAlpha = floor( ((SortedC[icand].Alpha) * (cos(SortedC[icand].Delta)) / (PCV.DeltaAlpha))  + (cc2 * 0.5) );

		    /* Assign the DELTA index to the candidate event */  
		    SortedC[icand].iDelta = floor( ((SortedC[icand].Delta) / (PCV.DeltaDelta))  + (cc3 * 0.5) );
		  }

		/* Assign the F1DOT index to the candidate event */
		SortedC[icand].iF1dot = floor( ((SortedC[icand].F1dot) / (PCV.DeltaF1dot))  + (cc4 * 0.5)  );
		
		/* Keep the original ordering before sort to refer the orignal data later. */
		SortedC[icand].iCand=icand;

		/* List for sorting purposes, gives speed while rearranging items in memory */
		SortedCListi[icand]=icand;
	    }

#if 0
	    /* Bruce's sorting function */
	    sortCandidates(SortedCListi, CLength);
#endif

#if 1
	    /* Holger's sorting function */
	    sortCandidates2(SortedCListi, 0, CLength-1);
	    fprintf(stdout,"%% Sorting of candidates finished.\n");
#endif
	    /* Initialise the first cell by the first candidate. */
	    icell = 0;
	    icand = 0;
	    global_cell[icell].iFreq = SortedC[SortedCListi[icand]].iFreq;
	    global_cell[icell].iDelta = SortedC[SortedCListi[icand]].iDelta;
	    global_cell[icell].iAlpha = SortedC[SortedCListi[icand]].iAlpha;
	    global_cell[icell].iF1dot = SortedC[SortedCListi[icand]].iF1dot;
	    global_cell[icell].CandID->data = SortedC[SortedCListi[icand]].iCand; 
	    global_cell[icell].nCand = 1;
	    	    
	    /* ------------------------------------------------------------------------------*/      
	    /* main loop over candidates  */
	    icell = 0;
	    for (icand=1; icand < CLength; icand++)
	      {
		/* Skip candidate events with 2F values below the threshold of TwoFthr. */
		if ( SortedC[SortedCListi[icand]].TwoF > PCV.TwoFthr ) 
		  {
		    if( SortedC[SortedCListi[icand]].iFreq  == global_cell[icell].iFreq  && 
			SortedC[SortedCListi[icand]].iDelta == global_cell[icell].iDelta &&
			SortedC[SortedCListi[icand]].iAlpha == global_cell[icell].iAlpha &&
			SortedC[SortedCListi[icand]].iF1dot == global_cell[icell].iF1dot ) 
		      {
			/* This candidate is in this cell. */
			INT4 lastFileIDinThisCell = SortedC[global_cell[icell].CandID->data].FileID;
			if( SortedC[SortedCListi[icand]].FileID != lastFileIDinThisCell ) 
			  { 
			    /* This candidate has a different file id from the candidates in this cell. */
			    LAL_CALL( add_int8_data( lalStatus, &(global_cell[icell].CandID), &(SortedC[SortedCListi[icand]].iCand) ), lalStatus );
			    global_cell[icell].nCand += 1;
			  }  
			else  
			  { 
			    /* This candidate has the same file id to one of candidates in this cell. */ 
			    /* 	       Because the array is already sorted in the DECREASING ORDER OF 2F,  */
			    /* 		 we do nothing here. */
			  }  /*if( SortedC[icand].FileID != lastFileIDinThisCell ) */
		      } /*  if( SortedC[icand].iFreq  == cell[icell].iFreq  && .. ) */ 
		    else 
		      {	  
			/* This candidate is outside of this cell. */
			icell++;
			
			/* Re-allocate Memory for more cells */
			if( icell >= sizecells ) {
			  sizecells = sizecells + ADDITIONAL_MEM;
			  LAL_CALL( RePrepareCells(lalStatus, &global_cell, sizecells, icell), lalStatus);
			}
			
			global_cell[icell].iFreq = SortedC[SortedCListi[icand]].iFreq;
			global_cell[icell].iDelta = SortedC[SortedCListi[icand]].iDelta;
			global_cell[icell].iAlpha = SortedC[SortedCListi[icand]].iAlpha;
			global_cell[icell].iF1dot = SortedC[SortedCListi[icand]].iF1dot;
			global_cell[icell].CandID->data = SortedC[SortedCListi[icand]].iCand;
			global_cell[icell].nCand = 1;
		
		      } /*  if( SortedC[icand].iFreq  == cell[icell].iFreq  && .. ) */ 
		    
		  } /* if ( SortedC[icand].TwoF > PCV.TwoFthr ) */
		
	      } /* for (icand=1; icand < CLength; icand++): loop over candidate list */      
	    
	    /* ---------------------------------------------------------------------------------------- */      
	    
	    
	    /* Get the information in each cell. */
	    ncell=icell+1; /* number of the cells in which more than or at least one candidate exists. */
	    /* cellList = (ListForSort *) LALCalloc( ncell, sizeof(ListForSort) ); */
	    cellListi = (INT8 *) LALCalloc(ncell, sizeof(INT8) );

	    for(icell=0;icell<ncell;icell++) {
	      LAL_CALL( get_info_of_the_cell( lalStatus, &global_cell[icell], SortedC), lalStatus);
	      cellListi[icell] = icell;
	    }  
	 
	    

	    /* -----------------------------------------------------------------------------------------*/      
	    /* Output results */
	    LAL_CALL( PrintResult( lalStatus, &PCV, global_cell, &ncell, SortedC, selectGrid, cellListi),lalStatus );
	    
	    
	    /* clean memory for cells */
	    LAL_CALL( FreeMemoryCellsOnly(lalStatus, global_cell, sizecells), lalStatus);
	    LALCheckMemoryLeaks(); 
	    
	    selectGrid++;

#ifdef SKYTESTMODE
	  }
#endif
	}
      }
    }
  }
	    

  /* -----------------------------------------------------------------------------------------*/      
  /* Clean-up */
  LAL_CALL( FreeMemory(lalStatus, &PCV, SortedC), lalStatus);
  LALCheckMemoryLeaks(); 

  
  return(POLKA_EXIT_OK);
 
} /* main() */



/* ########################################################################################## */
/* ------------------------------------------------------------------------------*/      
/* Initialize the code: allocate memory, set initial values.                     */
/* ------------------------------------------------------------------------------*/      
/*!
  Allocate memory for the cells.
  This function initialize the celldata variable.

  @param[in,out] lalStatus LALStatus*
  @param[out]    cell      CellData** CellData structure to be initialized
  @param[in]     CLength   INT8      Number of the cells
*/
void PrepareCells( LALStatus *lalStatus, CellData **cell, const INT8 CLength )
{
  INT8 icell, ncell;
  INT4 errflg = 0;

  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  ASSERT( *cell == NULL, lalStatus, POLKAC_ENONULL, POLKAC_MSGENONULL);

  *cell = (CellData *) LALCalloc( CLength, sizeof(CellData) );
  if( *cell == NULL ) {
    ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
  }

  for(icell=0;icell<CLength;icell++) {
    (*cell)[icell].CandID = NULL;
    (*cell)[icell].CandID = (struct int8_linked_list *) LALCalloc( 1, sizeof(struct int8_linked_list) );
    if( (*cell)[icell].CandID == NULL ) {
      errflg = 1;
      break;
    }
    (*cell)[icell].CandID->next = NULL;
    (*cell)[icell].iFreq = 0;
    (*cell)[icell].iDelta = 0;
    (*cell)[icell].iAlpha = 0;
    (*cell)[icell].iF1dot = 0;
    (*cell)[icell].nCand = 0;
    (*cell)[icell].Freq = 0.0;
    (*cell)[icell].Delta = 0.0;
    (*cell)[icell].Alpha = 0.0;
    (*cell)[icell].F1dot = 0.0;
    (*cell)[icell].significance = 0;
  }


  if( errflg != 0 ) {
    ncell = icell;
    for(icell=0;icell<ncell;icell++) {
      LALFree( (*cell)[icell].CandID );
    }
    ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
  }


  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
} /* PrepareCells() */




/* ########################################################################################## */
/* ------------------------------------------------------------------------------*/      
/* Re-allocate memory, set initial values.                     */
/* ------------------------------------------------------------------------------*/      
/*!
  Re-Allocate memory for the cells and initialize the additional celldata variables.
*/
 
void RePrepareCells( LALStatus *lalStatus,	/**< LALStatus* pointer */
                     CellData **cell,		/**< CellData structure to be initialized */
                     const INT8 CLength,	/**< Number of the cells */
                     const INT8 iposition	/**< FIXME: !TO BE DOCUMENTED! */
                     )
{
  INT8 icell, ncell;
  INT4 errflg = 0;
  CellData *tmp;
  
  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  /*  ASSERT( *cell == NULL, lalStatus, POLKAC_ENONULL, POLKAC_MSGENONULL); */


  tmp = (CellData *)LALRealloc (*cell, ( CLength * sizeof(CellData)) );
  if ( !tmp ) 
    { 
      XLALPrintError("Could not re-allocate memory for cells \n\n");
      ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
    }
  *cell = tmp;

  if( *cell == NULL ) {
    ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
  }

  for(icell=iposition;icell<CLength;icell++) {
    (*cell)[icell].CandID = NULL;
    (*cell)[icell].CandID = (struct int8_linked_list *) LALCalloc( 1, sizeof(struct int8_linked_list) );
    if( (*cell)[icell].CandID == NULL ) {
      errflg = 1;
      break;
    }
    (*cell)[icell].CandID->next = NULL;
    (*cell)[icell].iFreq = 0;
    (*cell)[icell].iDelta = 0;
    (*cell)[icell].iAlpha = 0;
    (*cell)[icell].iF1dot = 0;
    (*cell)[icell].nCand = 0;
    (*cell)[icell].Freq = 0.0;
    (*cell)[icell].Delta = 0.0;
    (*cell)[icell].Alpha = 0.0;
    (*cell)[icell].F1dot = 0.0;
    (*cell)[icell].significance = 0;
  }


  if( errflg != 0 ) {
    ncell = icell;
    for(icell=iposition;icell<ncell;icell++) {
      LALFree( (*cell)[icell].CandID );
    }
    ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
  }


  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
} /* RePrepareCells() */







/* ########################################################################################## */
/*! 
  Output results 
*/
void PrintResult(LALStatus *lalStatus, 		/**< LALStatus pointer */
		 const PolkaConfigVars *CLA, 	/**< PolkaConfigVars* */
		 CellData *cell, 		/**< CellData* */
		 const INT8 *ncell, 		/**< Number of the cells */
		 CandidateList *CList, 		/**< CandidateList */
		 const INT4 cellgridnum, 	/**< FIXME: !TO BE DOCUMENTED! */
		 INT8 CellListi[])		/**< FIXME: !TO BE DOCUMENTED! */
{
  INT8 icell;
  CHAR fnameSigTime[256]; /* Time variation of 2F of some significant outliers. */
  CHAR fnameSigCell[256]; /* Cell information of some significant outliers*/
  CHAR fnameCoiTime[256];  /* Time variation of 2F of some coincident outliers. */
  CHAR fnameCoiCell[256];  /* Cell information of some coincident outliers*/
  /* The cell info of the maximum coincident event over each Frequency cell but all over the sky.*/
  CHAR fnameMaxOverSky[256]; 
  CHAR cgn[256];  /* Cell-grid number*/
  CHAR fnameAllCells[256]; /* Output file to write all cell information to */
  CHAR fnameMinCellCoin[256];

  FILE *fp = NULL, *fpSigTime = NULL, *fpSigCell = NULL, *fpCoiTime = NULL, *fpCoiCell = NULL, *fpMinCellCoin = NULL;
  INT8 *count;
  INT8 nc, nmax,idxmax = 0;
  INT8 idxmaxcoin = 0;
  REAL4 Sigmax = 0.0;
  
  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  ASSERT( cell != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( CLA != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( CList != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);

  sprintf(fnameSigTime,"polka_significant_outlier_2FofTime_%d", cellgridnum);
  sprintf(fnameSigCell,"polka_significant_outlier_CellData_%d", cellgridnum);
  sprintf(fnameCoiTime,"polka_coincident_outlier_2FofTime_%d", cellgridnum);
  sprintf(fnameCoiCell,"polka_coincident_outlier_CellData_%d", cellgridnum);
  sprintf(fnameMaxOverSky,"polka_maxcoincident_over_each_freqcell_and_allsky_%d", cellgridnum);
  sprintf(cgn,"_%02d_cg", cellgridnum);
  strcpy(fnameAllCells,CLA->OutputFile);
  strcat(fnameAllCells,cgn);
  /* ------------------------------------------------------------- */

  
  /* This is the number of the maximum coincidences. */
  nmax=-1;

  /* elimininated qsort to find maximum by: */

  for(icell=0;icell<(*ncell);icell++) {
    if(nmax <= cell[CellListi[icell]].nCand){
      
      if (nmax == cell[CellListi[icell]].nCand){
	if (cell[idxmaxcoin].significance < cell[CellListi[icell]].significance){
	  nmax=cell[CellListi[icell]].nCand;
	  idxmaxcoin=CellListi[icell];
	}
      }

      else {
	nmax=cell[CellListi[icell]].nCand;
	idxmaxcoin=CellListi[icell];
      }

    }
  }


 
  /* Allocate memory for the histogram of coincidences */
  if( (count = (INT8 *) LALCalloc( (size_t) (nmax + 1), sizeof(INT8))) == NULL ) {
    XLALPrintError("Could not allocate Memory! \n");
    ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
  }
  


#ifdef SKYTESTMODE
  /* ------------------------------------------------------------- */
  /* Print out to the user-specified output file all the information in all the cell. 
     This file can be too huge to be tractable.*/
  if( (fp = fopen(fnameAllCells,"w")) == NULL ) 
    {
      XLALPrintError("\n Cannot open file %s\n",fnameAllCells); 
      ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
    }
  /* output for all the cells */
  print_info_of_the_cell( lalStatus->statusPtr, fp, cell, 0,(*ncell),0,0);
  BEGINFAIL(lalStatus) {fclose(fp);} ENDFAIL(lalStatus);

#endif


  /* compute the histogram of coincidences and find the most significant cell */
  for(icell=0;icell<(*ncell);icell++) {
    nc=cell[CellListi[icell]].nCand;
    count[nc] += 1;
    if( Sigmax < cell[CellListi[icell]].significance) {
      Sigmax = cell[CellListi[icell]].significance;
      idxmax = CellListi[icell];
    }
  }


  /* ------------------------------------------------------------- */
  /* output the summary table including the histogram of coincidences and 
     candidate events of most coincident cell */
  if(lalDebugLevel < 3 ) {
    fprintf(stdout,"%% Most significant cell : freq [Hz]\tdec [rad]\tra [rad]  \tF1dot \t\t   #[events]\tSig" "\n");
    fprintf(stdout, "%%\t\t\t   ");
    TRY( print_info_of_the_cell( lalStatus->statusPtr, stdout, cell, idxmax,idxmax+1,0,0), lalStatus);
    fprintf(stdout,"%% Most coincident cell  : freq [Hz]\tdec [rad]\tra [rad]  \tF1dot \t\t   #[events]\tSig" "\n");
    fprintf(stdout, "%%\t\t\t   ");
    TRY( print_info_of_the_cell( lalStatus->statusPtr, stdout, cell, idxmaxcoin,idxmaxcoin+1,0,0), lalStatus);

    fprintf(stdout,"%% # of coincidences: \n");
    for(nc=0;nc<=nmax;nc++) {
      fprintf(stdout,"%" LAL_INT8_FORMAT ,nc);
    }

    fprintf(stdout,"\n");
    fprintf(stdout,"%% # of cells       : \n");
    for(nc=0;nc<=nmax;nc++) { 
      fprintf(stdout, "%" LAL_INT8_FORMAT,count[nc]);
    }
    
    fprintf(stdout,"\n%%\n%% Candidate events of most coincident cell : \n%% data-seg \tfreq [Hz]\tdec [rad]\tra [rad]  \tF1dot[Hz/s]\t\t2F" "\n");
    TRY( print_cand_of_most_coin_cell( lalStatus->statusPtr, &cell[idxmaxcoin], CList), lalStatus);

  
  }
  LALFree( count );

 
  if( CLA->AutoOut || cell[idxmaxcoin].nCand >= CLA->Nthr ) 
    {
      if( (fpCoiCell = fopen(fnameCoiCell,"w")) == NULL || (fpCoiTime = fopen(fnameCoiTime,"w")) == NULL )
	{ 
	  XLALPrintError("\n Cannot open file %s or %s\n",fnameCoiCell,fnameCoiTime); 
	  exit(POLKA_EXIT_ERR); 
	}
    }

  if( CLA->AutoOut || Sigmax > CLA->Sthr ) 
    {
      if( (fpSigCell = fopen(fnameSigCell,"w")) == NULL || (fpSigTime = fopen(fnameSigTime,"w")) == NULL )
	{ 
	  XLALPrintError("\n Cannot open file %s or %s\n",fnameSigCell,fnameSigTime); 
	  exit(POLKA_EXIT_ERR); 
	}
    }



  /* ------------------------------------------------------------- */
  if( CLA->AutoOut ) 
    {  /* Output the info of the most significant and the most coincident event. */

      /* Output the info of the most coincident event. */
      /* Information of the cell. */
      print_info_of_the_cell( lalStatus->statusPtr, fpCoiCell, cell, idxmaxcoin,idxmaxcoin+1,0,0);
      BEGINFAIL(lalStatus) {fclose(fpCoiCell);} ENDFAIL(lalStatus);
      /* Print F stat from each file contributed to this cell. */
      print_Fstat_of_the_cell( lalStatus->statusPtr, fpCoiTime, cell, CList,  idxmaxcoin,idxmaxcoin+1,0,0 );
      BEGINFAIL(lalStatus) {fclose(fpCoiTime);} ENDFAIL(lalStatus);
    
      /* Output the info of the most significant event. */
      /* Information of the cell. */
      print_info_of_the_cell( lalStatus->statusPtr, fpSigCell, cell, idxmax,idxmax+1,0,0);
      BEGINFAIL(lalStatus) {fclose(fpSigCell);} ENDFAIL(lalStatus);
      /* Print F stat from each file contributed to this cell. */
      print_Fstat_of_the_cell( lalStatus->statusPtr, fpSigTime, cell, CList, idxmax,idxmax+1,0,0 );
      BEGINFAIL(lalStatus) {fclose(fpSigTime);} ENDFAIL(lalStatus);


    } /* if( CLA->AutoOut ) */ 
  else 
    {
      /* output only on outliers larger than Nthr on number of coincidences and Sthr on significance.*/
      if( cell[CellListi[0]].nCand >= CLA->Nthr ) 
	{
	  
	  /* Information of the cell. */
	  print_info_of_the_cell( lalStatus->statusPtr, fpCoiCell, cell, 0, 0, 0, CLA->Nthr);
	  BEGINFAIL(lalStatus) {fclose(fpCoiCell);} ENDFAIL(lalStatus);
	  /* Print F stat from each file contributed to this cell. */
	  print_Fstat_of_the_cell( lalStatus->statusPtr, fpCoiTime, cell, CList, 0, 0, 0, CLA->Nthr );
	  BEGINFAIL(lalStatus) {fclose(fpCoiTime);} ENDFAIL(lalStatus);

	} /* if( cell[0].nCand > CLA->Nthr ) */
      
      
      /* ------------------------------------------------------------- */
      /* Second Sort arrays of candidates based on significance, if necessary. */ 
      /* output only on outliers */
      if( Sigmax > CLA->Sthr ) 
	{
	  qsort(cell, (size_t) (*ncell), sizeof(CellData), compareSignificances);

	  /* Information of the cell. */
	  print_info_of_the_cell( lalStatus->statusPtr, fpSigCell, cell, 0, 0, CLA->Sthr, 0);
	  BEGINFAIL(lalStatus) {fclose(fpSigCell);} ENDFAIL(lalStatus);
	  /* Print F stat from each file contributed to this cell. */
	  print_Fstat_of_the_cell( lalStatus->statusPtr, fpSigTime, cell, CList, 0, 0, CLA->Sthr, 0 );
	  BEGINFAIL(lalStatus) {fclose(fpSigTime);} ENDFAIL(lalStatus);

	} /* if( cell[0].significance > CLA->Sthr ) */
    } /* else of if( CLA->AutoOut ) */


  if( CLA->AutoOut || cell[idxmaxcoin].nCand >= CLA->Nthr ) {
    fclose(fpCoiTime);
    fclose(fpCoiCell);
  }

  if( CLA->AutoOut || Sigmax > CLA->Sthr ) { 
    fclose(fpSigTime);
    fclose(fpSigCell);
  }


  /* ---------------------------------------------------------------------------------- */
  /* Output the maximum coincident event over each frequency cell and over all the sky. */

  /* faster sorting. */
  sortFreqCells2(CellListi, 0, (*ncell)-1 );
  fprintf(stdout,"%% Sorting of cells finished.\n");
  
  {
    INT8 prev_iFreq = -1;

    if( ( fp = fopen(fnameMaxOverSky,"w") ) == NULL ) {
      { 
	XLALPrintError("\n Cannot open file %s or %s\n",fnameCoiCell,fnameCoiTime); 
	exit(POLKA_EXIT_ERR); 
      }
    }

    if( !strcmp(CLA->EahRun,"S5R1a") ) {
      for( icell=0; icell<(*ncell); icell++ ) {
        if( cell[CellListi[icell]].iFreq != prev_iFreq ) {
          print_info_of_cell_and_ifo_S5R1a( lalStatus->statusPtr, fp, cell, CList, CellListi[icell], CellListi[icell]+1, 0, 0);
          BEGINFAIL(lalStatus) {fclose(fp);} ENDFAIL(lalStatus);
        }
        prev_iFreq = cell[CellListi[icell]].iFreq;
      }
    }
    else {
      if( !strcmp(CLA->EahRun,"S4R2a") ) {
	for( icell=0; icell<(*ncell); icell++ ) {
	  if( cell[CellListi[icell]].iFreq != prev_iFreq ) {
	    print_info_of_cell_and_ifo_S4R2a( lalStatus->statusPtr, fp, cell, CList, CellListi[icell], CellListi[icell]+1, 0, 0);
	    BEGINFAIL(lalStatus) {fclose(fp);} ENDFAIL(lalStatus);
	  }
	  prev_iFreq = cell[CellListi[icell]].iFreq;
	}
      }
      else {

	/* ---------------------------------------------------------------------------------- */
	/* Ouput individual candidate events of each cells with more than 'MinCellCoin' coincidences */
	if( CLA->MinCellCoin > 0 ) {
	  for( icell=0; icell<(*ncell); icell++ ) {
            if( cell[CellListi[icell]].iFreq != prev_iFreq ) {
              print_info_of_the_cell( lalStatus->statusPtr, fp, cell, CellListi[icell], CellListi[icell]+1, 0, 0);
              BEGINFAIL(lalStatus) {fclose(fp);} ENDFAIL(lalStatus);

	      if( cell[CellListi[icell]].nCand >= CLA->MinCellCoin )
		{
		  sprintf(fnameMinCellCoin,"ZP_G%02d_COIN%02d_CID%" LAL_INT8_FORMAT "_.dat", cellgridnum, cell[CellListi[icell]].nCand, CellListi[icell]);

		  if( (fpMinCellCoin = fopen(fnameMinCellCoin,"w")) == NULL )
		    {
		      XLALPrintError("\n Cannot open file %s or %s\n",fnameCoiCell,fnameCoiTime);
		      exit(POLKA_EXIT_ERR);
		    }
		  print_Fstat_of_the_cell( lalStatus->statusPtr, fpMinCellCoin, cell, CList, CellListi[icell], CellListi[icell]+1, 0, 0 );
		  fclose(fpMinCellCoin);
		}
	      }
            prev_iFreq = cell[CellListi[icell]].iFreq;
          }

	}
	else {
	  for( icell=0; icell<(*ncell); icell++ ) {
	    if( cell[CellListi[icell]].iFreq != prev_iFreq ) {
	      print_info_of_the_cell( lalStatus->statusPtr, fp, cell, CellListi[icell], CellListi[icell]+1, 0, 0);
	      BEGINFAIL(lalStatus) {fclose(fp);} ENDFAIL(lalStatus);
	      }
	    prev_iFreq = cell[CellListi[icell]].iFreq;
	  }

	}
      }
    }
    fclose(fp);
  }
  

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
} /* PrintResult() */




/* ########################################################################################## */
/*! 
  Print_info_of_the_cell

  This function basically shows the information of the outliers, where 
  those are in the parameters space and how coincident those are.
 
  Print out into FILE \b fp the infos of the cells 
  \li whose indices are between \b icell_start and \b icell_end, and 
  \li in which numbers of the events are above \b ncand_thr, and 
  \li in which significances are above \b sig_thr.
*/
void print_info_of_the_cell( LALStatus *lalStatus, 
			     FILE *fp, 
			     const CellData *cd, 
			     const INT8 icell_start, 
			     const INT8 icell_end, 
			     const REAL8 sig_thr, 
			     const REAL8 ncand_thr )
{
  INT8 icell;

  INITSTATUS(lalStatus);
  ASSERT( cd != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);

  icell = icell_start;
  while( icell < icell_end && 
	 cd[icell].significance > sig_thr && 
	 cd[icell].nCand > ncand_thr ) 
    {
      fprintf(fp,"%.6f\t%.6f\t%.6f\t%g\t\t%d\t%.6f\n", 
	      cd[icell].Freq, cd[icell].Delta, cd[icell].Alpha, cd[icell].F1dot, cd[icell].nCand, cd[icell].significance);
      icell++;
    }


  RETURN (lalStatus);
} /* void print_info_of_the_cell() */




/* ########################################################################################## */
/*!
  Free memory 

  Free Configuration variables \b CLA, CellData variable \b cell, CandidateList var \b CList.
*/
void FreeMemory( LALStatus *lalStatus,	/**< LALStatus*  pointer */ 
                 PolkaConfigVars *CLA, 	/**< configuration variables structure */
                 CandidateList *CList) /**< CandidateList structure */
{
  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  FreeConfigVars( lalStatus->statusPtr, CLA );

  if( CList != NULL ) LALFree(CList);

  /* This is now being done by FreeMemoryCellsOnly ! ############
  if( cell != NULL ) {
    for(icell=0;icell<datalen;icell++) {
      TRY( delete_int8_linked_list( lalStatus->statusPtr, cell[icell].CandID ), lalStatus );
    }
    LALFree(cell);
  }
  */

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
} /* FreeMemory */


/* ########################################################################################## */
/*!
  Free memory 

  Free Configuration variables \b CLA, CellData variable \b cell, CandidateList var \b CList.
*/
void FreeMemoryCellsOnly( LALStatus *lalStatus, /**< LALStatus*  pointer */
			  CellData *cell, 	/**< CellData structure */
			  const INT8 datalen)	/**< Number of the cells */
{
  INT8 icell;

  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  /* FIX (?) ME:  
     This part takes really long, when lalDebugLevel = 3. I do not know why.*/
  if( cell != NULL ) {
    for(icell=0;icell<datalen;icell++) {
      TRY( delete_int8_linked_list( lalStatus->statusPtr, cell[icell].CandID ), lalStatus );
    }
    LALFree(cell);
  }
    
  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
} /* FreeMemory */


/* ########################################################################################## */
/*!
  Free Configuration variables \b CLA.

  @param[in,out] lalStatus LALStatus* 
  @param[in]     CLA       PolkaConfigVars* configuration variables structure
*/
void FreeConfigVars(LALStatus *lalStatus, PolkaConfigVars *CLA )
{
  UINT4 k;

  INITSTATUS(lalStatus);

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

  RETURN (lalStatus);
} /* FreeCOnfigVars */


/* ########################################################################################## */
/*!
  add data to linked structure 

  @param[in,out] lalStatus   LALStatus* 
  @param[in,out] list_ptr    int8_linked_list**
  @param[in]     data        INT8*
*/
void add_int8_data(LALStatus *lalStatus, struct int8_linked_list **list_ptr, const INT8 *data)
{
  struct int8_linked_list *p = NULL;

  INITSTATUS(lalStatus);

  p = (struct int8_linked_list *) LALMalloc(sizeof(struct int8_linked_list));
  if( p == NULL ) {
    XLALPrintError("Could not allocate Memory! \n");
    ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
  }
  p->data = *(data);
  p->next = *list_ptr;
  *list_ptr = p;

  RETURN (lalStatus);
} /* void add_int8_data() */


/* ########################################################################################## */
/*!
  delete a linked structure 

  @param[in,out] lalStatus   LALStatus* 
  @param[in]     list_ptr    int8_linked_list*
*/
void delete_int8_linked_list( LALStatus *lalStatus, struct int8_linked_list *list_ptr )
{
  INT8 ic;
  struct int8_linked_list *q;

  INITSTATUS(lalStatus);

  ic = 0;
  while( list_ptr !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) {  
    q = list_ptr->next;
    LALFree( list_ptr );
    list_ptr = q;
    ic++;
  }
  if( ic >  LINKEDSTR_MAX_DEPTH ) {
    XLALPrintError("Maximum depth of linked structure reached!");
    exit(POLKA_EXIT_ERR);
  }

  RETURN (lalStatus);
} /* void delete_int8_linked_list() */


/* ########################################################################################## */
/*!
  get info of this cell. 

  We have indices of the candidate events contained in each cell 
  before the call of this function. This function computes the joint 
  significance, average alpha, average delta, average spin-down 
  and average frequency of the events in each cell and stores them 
  into a CellData structure.

  @param[in,out] lalStatus   LALStatus* 
  @param[in,out] cd          CellData*
  @param[in]     CList       CandidateList*  
*/
void get_info_of_the_cell( LALStatus *lalStatus, CellData *cd, const CandidateList *CList )
{
  INT8 idx, ic;
  REAL8 lfa;
  struct int8_linked_list *p;

  INITSTATUS(lalStatus);
  ASSERT( cd != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( CList != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);

  p = cd->CandID;

  ic = 0;
  while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) { 
    idx = p->data;
    lfa = CList[idx].TwoF/2.0 - log(1.0 + CList[idx].TwoF/2.0);
    cd->significance += lfa;
    cd->F1dot += CList[idx].F1dot;
    cd->Alpha += CList[idx].Alpha;
    cd->Delta += CList[idx].Delta;
    cd->Freq += CList[idx].f;

#ifdef SKYTESTMODE
    /* Testing for sky-cells having less than 17 coincidences */
    if( cd->nCand < 17 ) {
      if( CList[idx].Alpha > LAL_TWOPI ) {
	fprintf(stdout,"skypointx %d %.6f %.6f\n",cd->nCand,CList[idx].Alpha - LAL_TWOPI,CList[idx].Delta);
      }
      else {
	if( CList[idx].Alpha < 0 ) {
	  fprintf(stdout,"skypointx %d %.6f %.6f\n",cd->nCand,CList[idx].Alpha + LAL_TWOPI,CList[idx].Delta);
	}
	else {
	  fprintf(stdout,"skypointx %d %.6f %.6f\n",cd->nCand,CList[idx].Alpha,CList[idx].Delta);
	}
      }
    }
#endif

    p = p->next;
    ic++;
  }
  if( ic >  LINKEDSTR_MAX_DEPTH ) {
    XLALPrintError("Maximum depth of linked structure reached!");
    exit(POLKA_EXIT_ERR);
  }

  cd->F1dot /= cd->nCand;
  cd->Alpha /= cd->nCand;
  cd->Delta /= cd->nCand;
  cd->Freq  /= cd->nCand;
  
  if( cd->Alpha > LAL_TWOPI ) {
    cd->Alpha = cd->Alpha - LAL_TWOPI;  
  }
  
  if( cd->Alpha < 0 ) {
    cd->Alpha = cd->Alpha + LAL_TWOPI;
  }

  RETURN (lalStatus);
} /* void get_info_of_the_cell() */



/* ########################################################################################## */
/*!
  print candidates of most coincident cell. 

  We have indices of the candidate events contained in each cell 
  before the call of this function. This function returns all
  the candidates belonging to the most coincident cell.

  @param[in,out] lalStatus   LALStatus* 
  @param[in,out] cd          CellData*
  @param[in]     CList       CandidateList*  
*/

void print_cand_of_most_coin_cell( LALStatus *lalStatus, CellData *cd, const CandidateList *CList )
{
  INT8 idx, ic;
  REAL8 AlphaTmp=0;
  struct int8_linked_list *p;

  INITSTATUS(lalStatus);
  ASSERT( cd != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( CList != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);  

  p = cd->CandID;

  ic = 0;
  while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) { 
    idx = p->data;
    AlphaTmp = CList[idx].Alpha;
    if( CList[idx].Alpha > LAL_TWOPI ) {
      AlphaTmp = CList[idx].Alpha - LAL_TWOPI;
    }
    if( CList[idx].Alpha < 0 ) {
      AlphaTmp = CList[idx].Alpha + LAL_TWOPI;
    }
    fprintf(stdout,"  %d\t\t%.6f\t%.6f\t%.6f\t%g\t\t%.6f\n", 
	    CList[idx].FileID, CList[idx].f, CList[idx].Delta, AlphaTmp, CList[idx].F1dot, CList[idx].TwoF);
    
    p = p->next;
    ic++;
  }

  if( ic >  LINKEDSTR_MAX_DEPTH ) {
    XLALPrintError("Maximum depth of linked structure reached!");
    exit(POLKA_EXIT_ERR);
  }

  RETURN (lalStatus);
} /* void print_cand_of_most_coin_cell() */



/* ########################################################################################## */
/*!
  print F stat. 

  Print out into FILE \b fp the F stats of the cells 
  \li whose indices are between \b icell_start and \b icell_end, and 
  \li in which numbers of the events are above \b ncand_thr, and 
  \li in which significances are above \b sig_thr.

  This function basically shows how F statistics variers from file to file.

  @param[in,out] lalStatus   LALStatus* 
  @param[in]     fp          FILE*
  @param[in]     cd          CellData*
  @param[in]     CList       CandidateList*  
  @param[in]     icell_start INT8  Starting index of a cell
  @param[in]     icell_end   INT8  Ending index of a cell
  @param[in]     sig_thr     REAL8 Threshold on significance of the candidate events 
  above which results will be printed out.
  @param[in]     ncand_thr   REAL8 Threshold on number of the candidate events 
  above which results will be printed out.
*/
void print_Fstat_of_the_cell( LALStatus *lalStatus, 
			      FILE *fp, 
			      const CellData *cd, 
			      const CandidateList *CList, 
			      const INT8 icell_start, 
			      const INT8 icell_end, 
			      const REAL8 sig_thr, 
			      const REAL8 ncand_thr )
{
  INT8 idx, ic, icell;
  struct int8_linked_list *p;

  INITSTATUS(lalStatus);
  ASSERT( cd != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( CList != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);

  icell = icell_start;
  while( icell < icell_end && 
	 cd[icell].significance > sig_thr && 
	 cd[icell].nCand > ncand_thr ) 
    {

      p = cd[icell].CandID;
      
      ic = 0;
      while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) { 
	idx = p->data;

	fprintf(fp,"%" LAL_INT8_FORMAT "\t%d\t%.6f\t%.6f\t%.6f\t%g\t%.6f\n", 
		icell, CList[idx].FileID, CList[idx].f, CList[idx].Alpha, CList[idx].Delta, CList[idx].F1dot, CList[idx].TwoF );

	p = p->next;
	ic++;
      } /*   while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) {  */

      if( ic >  LINKEDSTR_MAX_DEPTH ) {
	XLALPrintError("Maximum depth of linked structure reached!");
	exit(POLKA_EXIT_ERR);
      }

      icell++;
    } /*   while( icell < icell_end && ...  */

  RETURN (lalStatus);
} /* void print_Fstat_of_the_cell( ) */


/* ########################################################################################## */
/*!
  Sorting function to sort cells indices in the 
  DECREASING ORDER OF a significanc and a number of candidate events in cells.

  Compare two cells in terms of the significance.
  If those are the same, then compare them in terms of the number of the candidates.
  If we use qsort, we will have cells ordered as 
  cell[0], cell[1], ....
  where cell[0] has a largher significance (or if it is equal to that of cell[1]
  then a larger number of the candidate events) than that of cell[1].


  @param[in] a CellData* to be compared. 
  @param[in] b CellData* to be compared. 
  @return If a<b, return -1, if a==b return 0, otherwise return 1. 
*/
int compareSignificances(const void *a, const void *b)
{
  const CellData *ip = a;
  const CellData *jp = b;
  int res;

  REAL8 F1, F2;
  F1=ip->significance;
  F2=jp->significance;
  /* I put F1 and F2 inversely, because I would like to get decreasingly-ordered set. */ 
  res = compareREAL8arrays( &F2,  &F1, 1);
  if( res == 0 ) {
    INT8 n1, n2;
    n1=ip->nCand;
    n2=jp->nCand;
    /* I put n1 and n2 inversely, because I would like to get decreasingly-ordered set. */ 
    res = compareINT8arrays( &n2,  &n1, 1);
  } 


  return res;
} /* int compareSignificances() */




/* ############################################################################### */

/*!
  Compare two REAL8 arrays of the same size \b n.
  First compare ap[0] and bp[0]. If ap[0] < bp[0], then 
  return -1. If ap[0] > bp[0], then return 1. If 
  ap[0] == bp[0], then compare ap[1] with bp[1]. Do the 
  same untill we reach the stage where we compare ap[n-1] 
  with bp[n-1]. If ap[n-1]==bp[n-1], then return 0.

  @param[in] ap REAL8 array to be compared
  @param[in] bp REAL8 array to be compared
  @param[in] n  Size of the array
  @return If ap<bp, return -1, if ap==bp return 0, otherwise return 1. 
*/
int compareREAL8arrays(const REAL8 *ap, const REAL8 *bp, size_t n) 
{
  if( (*ap) == (*bp) ) { 
    if ( n > 1 ){  
      return compareREAL8arrays( ap+1, bp+1, n-1 );
    } else {
      return 0;
    }
  }
  if ( (*ap) < (*bp) ) 
    return -1;    
  return 1;
} /* int compareREAL8arrays() */



/*!
  Compare two INT4 arrays of the same size \b n.
  First compare ap[0] and bp[0]. If ap[0] < bp[0], then
  return -1. If ap[0] > bp[0], then return 1. If
  ap[0] == bp[0], then compare ap[1] with bp[1]. Do the
  same untill we reach the stage where we compare ap[n-1]
  with bp[n-1]. If ap[n-1]==bp[n-1], then return 0.

  @param[in] ap INT4 array to be compared
  @param[in] bp INT4 array to be compared
  @param[in] n  Size of the array
  @return If ap<bp, return -1, if ap==bp return 0, otherwise return 1.
*/
int compareINT4arrays(const INT4 *ap, const INT4 *bp, size_t n)
{
  if( (*ap) == (*bp) ) {
    if ( n > 1 ){
      return compareINT4arrays( ap+1, bp+1, n-1 );
    } else {
      return 0;
    }
  }
  if ( (*ap) < (*bp) )
    return -1;
  return 1;
} /* int compareINT4arrays() */


/*!
  Compare two INT8 arrays of the same size \b n.
  First compare ap[0] and bp[0]. If ap[0] < bp[0], then 
  return -1. If ap[0] > bp[0], then return 1. If 
  ap[0] == bp[0], then compare ap[1] with bp[1]. Do the 
  same untill we reach the stage where we compare ap[n-1] 
  with bp[n-1]. If ap[n-1]==bp[n-1], then return 0.

  @param[in] ap INT8 array to be compared
  @param[in] bp INT8 array to be compared
  @param[in] n  Size of the array
  @return If ap<bp, return -1, if ap==bp return 0, otherwise return 1. 
*/
int compareINT8arrays(const INT8 *ap, const INT8 *bp, size_t n) 
{
  if( (*ap) == (*bp) ) { 
    if ( n > 1 ){  
      return compareINT8arrays( ap+1, bp+1, n-1 );
    } else {
      return 0;
    }
  }
  if ( (*ap) < (*bp) ) 
    return -1;    
  return 1;
} /* int compareINT8arrays() */


/* ########################################################################################## */
/*!
  Read Candidate File(s) and store the events into CandidateList str \b CList.

  If an input directory (\b CLA->InputDir ) is specified, call \b GetFilesListInThisDir() 
  to find the list of the files, and 
  then call \b ReadCandidateListFromZipFile() to fill CanidateList structure \b CList.

  If an input file (\b CLA->FstatsFile ) is specified, call ReadOneCandidateFile() to fill 
  CanidateList structure \b CList;

  @param[in,out] lalStatus LALStatus* 
  @param[out]    CList     CandidateList**  Candidate events struecture to be filled
  @param[in,out] CLA       PolkaConfigVars* Configuration variables structure
  @param[out]    clen      INT8*           The total number of the candidate events in the files.
*/
void 
ReadCandidateFiles(LALStatus *lalStatus, 
		   CandidateList **CList, 
		   PolkaConfigVars *CLA, 
		   INT8 *clen)
{
  UINT4 kc;
  /* UINT4 *CLenFthr;*/
  /* REAL8 percentage = 0; */

  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  ASSERT( CLA != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( *CList == NULL, lalStatus, POLKAC_ENONULL, POLKAC_MSGENONULL);

  if( (CLA->InputDir != NULL) && (CLA->BaseName != NULL) ) 
    {
      CLA->Filelist = NULL;
      TRY( GetFilesListInThisDir( lalStatus->statusPtr, 
				  CLA->InputDir, 
				  CLA->BaseName, 
				  &(CLA->Filelist), 
				  &(CLA->NFiles) ), 
	   lalStatus );
      
      *clen = 0;     /* We first have to set the candidate list length zero. */
      /**CLenFthr = 0;*/
      *CList = NULL; /* We first have to nullify the list. */
      for (kc=0;kc<CLA->NFiles;kc++)
	{
	  if( lalDebugLevel > 1 ) {
	    fprintf(stdout,"%s\n",CLA->Filelist[kc]);
	  }

	  if( CLA->UseUnzip ) 
	    {
#ifdef USE_UNZIP
	      {INT4 FileID = 2*kc; /* the factor 2 because we have 2 sections in each file. */
	      TRY( ReadCandidateListFromZipFile( lalStatus->statusPtr, 
						 CList, 
						 CLA->Filelist[kc], 
						 clen, 
						 &FileID), 
		   lalStatus);
	      }
#endif
	    } /* if( CLA->UseUnzip ) */
	  else 
	    {	      
	      TRY( ReadOneCandidateFileV2( lalStatus->statusPtr, 
					   CList, 
					   CLA->Filelist[kc], 
					   clen ), 
		   lalStatus);	      
	    }
	} /*       for (kc=0;kc<CLA->NFiles;kc++) */

    } /* if( (CLA->InputDir != NULL) && (CLA->BaseName != NULL) )  */
  else if ( CLA->FstatsFile != NULL ) 
    {
      *CList = NULL;
      TRY( ReadOneCandidateFile(lalStatus->statusPtr, CList, CLA->FstatsFile, clen, CLA->TwoFthr ), lalStatus );
      /* The last file is from last file.*/
      CLA->NFiles = (*CList)[*clen-1].FileID;
    } /* if( (CLA->InputDir != NULL) && (CLA->BaseName != NULL) )  */
  else 
    { /* We should not be here. */
      XLALPrintError("\nYou have to specify either input data directory or input data file.\n");
      exit(POLKA_EXIT_ERR);;
    }

  /* 
     percentage = ( (REAL8) *CLenFthr / *clen ) * 100.0;

     fprintf(stdout,"\n%%Number of the candidate events in this file/directory = %u.\n%% --- Threshold for 2F: %.3f\t Number of candidates kept: %u  or  %.3f%% --- \n",*clen, CLA->TwoFthr, *CLenFthr, percentage);
  */

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
} /* ReadCandidateFiles() */





/* ########################################################################################## */
/*!
  Get the list of the files which has the base name \b basename in 
  \b directory and store the list into \b filelist. 
  Count the number of the files \b nfiles. This function checks 
  \li if HAVE_GLOB_H has been defined. If not, the function aborts.

  @param[in,out]  lalStatus LALStatus*
  @param[in]      directory CHAR*   Directory name for which files list will be made.
  @param[in]      basename  CHAR*   The base name of the files to be listed. 
  @param[out]     filelist  CHAR*** The list of the files will be stored in this strucutre.
  @param[out]     nfiles    UINT4*  The number of the files which has the basename 
  \b basename in the \b directory.
*/
void 
GetFilesListInThisDir( LALStatus *lalStatus, 
		       const CHAR *directory, 
		       const CHAR *basename, 
		       CHAR ***filelist, 
		       UINT4 *nfiles )
{
#ifdef HAVE_GLOB_H   
  CHAR *command = NULL;
  UINT4 filenum=0;
  glob_t globbuf;
#endif

  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);
  ASSERT( directory != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( basename != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( *filelist == NULL, lalStatus, POLKAC_ENONULL, POLKAC_MSGENONULL);

#ifndef HAVE_GLOB_H   
  XLALPrintError("Cannot use GetFilesListInThisDir() without glob.");
  ABORT( lalStatus, POLKAC_EGLOB, POLKAC_MSGEGLOB);
#endif

  command = (CHAR*) LALMalloc( strlen(directory) + strlen("/*") + strlen(basename) + strlen("*") + 1 );
  if( command == NULL ) {
    ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
  }

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
      XLALPrintError ("\nNo Input files in directory %s ... Exiting.\n\n", directory);
      ABORT (lalStatus, POLKAC_ESYS, POLKAC_MSGESYS);
    }

  /* prepare memory for all filenames */
  *filelist = NULL;
  if ( ( *filelist = (CHAR**)LALCalloc(globbuf.gl_pathc, sizeof(CHAR*))) == NULL) {
    ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
  }
  while ( filenum < (UINT4) globbuf.gl_pathc) 
    {
      (*filelist)[filenum] = NULL;
      if ( ((*filelist)[filenum] = (CHAR*)LALCalloc(1, strlen(globbuf.gl_pathv[filenum])+1)) == NULL) {
	ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
      }
      strcpy((*filelist)[filenum],globbuf.gl_pathv[filenum]);
      filenum++;
    }
  globfree(&globbuf);

  *nfiles = filenum; /* remember this is 1 more than the index value */
#endif

  LALFree(command);

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
}




/* ########################################################################################## */
#ifdef USE_UNZIP
/*
TODO:
Check if *CList is either NULL or the memory of which is previously allocated by alloc() or the kind.
(how?).
*/
/*!
  Read the given zipped candidate 'Fstats'-file \b fname and append the events in the file to 
  the candidate-list \b CList. 
  This function is invoked only when \b USE_UNZIP is defined.
  The function aborts almost all the cases when the checks below failed.
  This function checks 
  \li if the file \b fname is readable. 
  \li if the number of the candidate event is smaller than the hardcoded number 8000000.
  \li if the file \b fname has the correct ending tag "%DONE".
  \li if the file has sections %1, %2 and %coincidence in this order.
  \li if the ranges of the values in the file are sensible.
  \li if the number of each row of the file is correct.
  \li if we could read all the events in the file.

  @param[in,out] lalStatus LALStatus* 
  @param[in,out] CList     CandidateList**  CandidateList str to be appended
  @param[in]     fname     CHAR* the name of the file to be read
  @param[in,out] candlen   INT8* total number of the candidate events so far. 
  This will be updated after reading the file \b fname. 
  @param[in]     FileID    INT4* The \b FileID of the file to be read. Assign a \b FildID 
  to each event and record which file a certain event comes from.
*/
void  
ReadCandidateListFromZipFile( LALStatus *lalStatus, 
			      CandidateList **CList, 
			      CHAR *fname, 
			      INT8 *candlen, 
			      const INT4 *FileID )
{
  FILE *fp;
  const INT8 max_num_candidates = 8000000; /* maximum tractable number of candidate events. */
  INT8 numlines;
  INT8 nread;
  REAL8 epsilon=1e-5;
  INT8 ic;
  INT8 length; /* length of file */
  CHAR *line, *endp; /* pointers to start and end of line */
  INT4 section = 0;    /* 0: non-POLKA, 1,2: IFO sections,
			 3: coincidence, 4: end-of-file */
  INT8 nlines[2] = {0,0}; /* number of events for each IFO */
  const INT4 MAX_SECS = 4;

  UzpBuffer uzpbuff;

  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);
  ASSERT( fname != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);

  /* check if file exists.*/
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      XLALPrintError("File %s doesn't exist!\n",fname);
      ABORT( lalStatus, POLKAC_ESYS, POLKAC_MSGESYS ); 
     }
  fclose(fp);

  /* Check if the memory to be allocated is not huge 
     (say, < 512MB. sizeof(CandidateList) ~ 60B. 512/60 = 8000000). */
  if( *candlen > max_num_candidates ) {
    XLALPrintError("\nMaximum number of candidate events reached.\n");
    XLALPrintError("\nWe have %u events while the maximum allowed number of events is %u.\n",*candlen,max_num_candidates);
    ABORT( lalStatus, POLKAC_ESYS, POLKAC_MSGESYS ); 
  }



  uzpbuff.strptr = NULL;

  /* ------------------------------------------------------------------------- */
  /*  Open and count the size of the candidates file */
  /* Read into buffer.  If this fails, we can't proceed. */
  if ( getfile( &uzpbuff, fname )  < 0 ) {
    if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
    XLALPrintError("Cannot read file %s . \n",fname);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  length = uzpbuff.strlength;
  line = uzpbuff.strptr;
  if ( !line || length == 0 || *line == '\0' ) {
    if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
    XLALPrintError ("Unknown format of the file  %s.\n\n", fname);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* ------------------------------------------------------------------------- */
  /* Check for correct ending tag.  If this isn't present, it is
     safest not to proceed (it greatly simplifies error trapping). */
  line += length;
  if ( ( length < 8 || strncmp( line - 8, "\n%DONE\r\n", 8 ) ) &&
       ( length < 7 || strncmp( line - 7, "\n%DONE\n", 7 ) ) ) {
    free(uzpbuff.strptr);
    XLALPrintError("File %s does not end with the DONE_MARKER. \n",fname);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
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
	XLALPrintError("Unknown format file %s.",fname);
	ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }
      section = 2;
      continue;
    } else if ( !strncmp( line, "%coincidence", 12 ) ) {
      if( section != 2 ) { /* We should have section 2 before 3. */
	if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
	XLALPrintError("Unknown format file %s.",fname);
	ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }
      section = 3;
      break; /* We are not interested in the section 3 here. */
    }  /*   if ( !strncmp( line, "%1", 2 ) ) {*/
 

    /* Do non-POLKA checks: */
    if ( section == 0 ) 
      {
	XLALPrintError("Unknown format file %s.",fname);
	if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
	ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
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
	XLALPrintError("Unknown format file %s.",fname);
	ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      } /*     if ( section == 0 )  */


    /* Done reading this line. */
  } /*   for ( line = uzpbuff.strptr; section < MAX_SECS; ... */
  /* ------------------------------------------------------------------------- */

  numlines = nlines[0] + nlines[1]; /* we have two sections. */

  if( numlines == 0 ) { /* This file is empty. Go to the next file.*/
    if( lalDebugLevel > 1 ) {
      XLALPrintError( "No candidate events in the file %s\n\n", fname);
    }
    free(uzpbuff.strptr);
    DETATCHSTATUSPTR (lalStatus);
    RETURN (lalStatus);
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
	  XLALPrintError("Could not allocate memory for candidate file %s\n\n", fname);
	  ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
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
            cl->FileID < 0                     ||
            cl->f < 0.0                        ||
            cl->TwoF < 0.0                     ||
            cl->Alpha <         0.0 - epsilon  ||
            cl->Alpha >   LAL_TWOPI + epsilon  ||
            cl->Delta < -0.5*LAL_PI - epsilon  ||
            cl->Delta >  0.5*LAL_PI + epsilon  ||
            !isfinite(cl->FileID)              ||
            !isfinite(cl->f)                   ||
            !isfinite(cl->Alpha)               ||
            !isfinite(cl->Delta)               ||
            !isfinite(cl->TwoF)
	    ) {
	  XLALPrintError(
			"Line %d of file %s has invalid values.\n"
			"First 255 chars are:\n"
			"%s\n"
			"1st and 4th field should be positive.\n" 
			"2nd field should lie between 0 and %1.15f.\n" 
			"3rd field should lie between %1.15f and %1.15f.\n"
			"All fields should be finite\n",
			ic+1, fname, line, (REAL8)LAL_TWOPI, (REAL8)-LAL_PI/2.0, (REAL8)LAL_PI/2.0);
	  LALFree ((*CList));
	  free( uzpbuff.strptr );
	  ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
	} /* end of the check of the range of the values.*/
      } /*     if ( section == 1 || section == 2 )  */
    /* Do POLKA coincidence-section checks. */
    else 
      { /* we should not be here */
	XLALPrintError("Unknown format file %s.",fname);
	if( uzpbuff.strptr != NULL ) free(uzpbuff.strptr);
	ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      } /* if ( section == 1 || section == 2 )  */
    
    /* Done reading this line and filling CList. */

           

    /* ------------------------------------------------------------------------- */
    /* check that we read 4 quantities with exactly the right format */
    if ( nread != 4 )
      {
	XLALPrintError ("Found %d not %d values on line %d in file '%s'\n"
		       "Line in question is\n%s",
		       nread, 4, ic+1, fname, line);               
	LALFree ((*CList));
	free( uzpbuff.strptr );
	ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }

  } /*   for ( line = uzpbuff.strptr; section < MAX_SECS; ... ) */

  free( uzpbuff.strptr ); /* uzpbuff is allocated by malloc() in getfile(). It is user's responsibility to free this. */


  if (ic != (*candlen) + numlines ) {
    XLALPrintError(
            "Read of file %s terminated after %d line but numlines=%d\n",
            fname, ic, numlines);
    LALFree((*CList));
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  if ( section != 3 ) {
    XLALPrintError(
            "Read of file %s terminated not by coincidence section but %s\n",
            fname, line);
    LALFree((*CList));
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }


  (*candlen) += numlines; /* total number of candidate so far */


  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);


} /* void  ReadCandidateListFromZipFile () */
#endif /* #ifdef USE_UNZIP */



/* ########################################################################################## */
/*!
  Read one candidate-events file and fill CandidateList structure \b CList. 
  Count the number of the candidate events and fill it in \b candlen.
  The function aborts almost all the cases when the checks below failed.
  This function checks 
  \li if the file \b fname is readable. 
  \li if the file \b fname has the correct ending tag "%DONE".
  \li if the ranges of the values in the file are sensible.
  \li if the number of each row of the file is correct.
  \li if we could read all the events in the file.

  This function prints the bytecounts and the checksum of the file \b fname.

  @param[in,out] lalStatus LALStatus* 
  @param[out]    CList     CandidateList** CandidateList str to be filled in this code 
  @param[in]     fname     CHAR* the name of the file to be read
  @param[in,out]    candlen   INT8* total number of the candidate events
*/
void  
ReadOneCandidateFileV2( LALStatus *lalStatus, 
		      CandidateList **CList, 
		      const CHAR *fname, 
		      INT8 *candlen )
{
  INT8 i;
  INT8 numlines;
  REAL8 epsilon=1e-5;
  CHAR line1[256];
  FILE *fp;
  INT8 nread;
  INT8 checksum=0;
  INT8 bytecount=0;


  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);
  ASSERT( fname != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);


  /* ------ Open and count candidates file ------ */
  i=0;
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      XLALPrintError("File %s doesn't exist!\n",fname);
      ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
     }
  while(fgets(line1,sizeof(line1),fp)) {
    UINT8 k;
    size_t len=strlen(line1);

    /* check that each line ends with a newline char (no overflow of
       line1 or null chars read) */
    if (!len || line1[len-1] != '\n') {
      XLALPrintError(
              "Line %d of file %s is too long or has no NEWLINE.  First 255 chars are:\n%s\n",
              i+1, fname, line1);
      fclose(fp);
      ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
     }

    /* increment line counter */
    i++;

    /* maintain a running checksum and byte count */
    bytecount+=len;
    for (k=0; k<len; k++)
      checksum+=(INT8)line1[k];
  }
  numlines=i;
  /* -- close candidate file -- */
  fclose(fp);     

  if ( numlines == 0) 
    {
      XLALPrintError ("ERROR: File '%s' has no lines so is not properly terminated by: %s", fname, DONE_MARKER);
      ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
    }

  /* output a record of the running checksun and byte count */
  fprintf(stdout, "%% %s: bytecount %" LAL_INT8_FORMAT " checksum %" LAL_INT8_FORMAT "\n", fname, bytecount, checksum);

  /* check validity of this Fstats-file */
  if ( strcmp(line1, DONE_MARKER ) ) 
    {
      XLALPrintError ("ERROR: File '%s' is not properly terminated by: %sbut has %s instead", fname, DONE_MARKER, line1);
      ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
    }
  else
    numlines --;        /* avoid stepping on DONE-marker */



  if( numlines == 0 ) { /* This file is empty. Go to the next file.*/
    if( lalDebugLevel > 1 ) {
      XLALPrintError( "No candidate events in the file %s\n\n", fname);
    }
    DETATCHSTATUSPTR (lalStatus);
    RETURN (lalStatus);
  } 

#if 0 /* Do we need to check this for INT8? */
  if ( numlines < 0  )
    {
      XLALPrintError("candidate length overflow (or indeed negative) = %ud!\n",numlines);
      exit(POLKA_EXIT_ERR);
    }/* check that we have candidates. */
#endif


  /* ------------------------------------------------------------------------- */
  /* reserve memory for fstats-file contents */
  if ( numlines > 0) 
    { 
      CandidateList *tmp;
      tmp = (CandidateList *)LALRealloc (*CList, ( *candlen + numlines )*sizeof(CandidateList));
      if ( !tmp ) 
	{ 
	  XLALPrintError("Could not allocate memory for candidate file %s\n\n", fname);
	  ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
	}
      *CList = tmp;
    }

  

  /* ------ Open and count candidates file ------ */
  i=0; /* append the new candidate events to the existing list. */
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      XLALPrintError("fopen(%s) failed!\n", fname);
      LALFree ((*CList));
      ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
    }
  while(i < numlines && fgets(line1,sizeof(line1),fp))
    {
      CHAR newline='\0';
      CandidateList *cl=&(*CList)[i+(*candlen)];

      if (strlen(line1)==0 || line1[strlen(line1)-1] != '\n') {
        XLALPrintError(
                "Line %d of file %s is too long or has no NEWLINE.  First 255 chars are:\n%s\n",
                i+1, fname, line1);
        LALFree ((*CList));
        fclose(fp);
	ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }
      
      nread = sscanf (line1, 
                     "%" LAL_INT4_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL8_FORMAT 
                     " %" LAL_REAL4_FORMAT "%c", 
                     &(cl->FileID), &(cl->f), &(cl->Alpha), &(cl->Delta), &(cl->F1dot), &(cl->TwoF), &newline );

      /* check that values that are read in are sensible */
      if (
          cl->FileID < 0                     ||
          cl->f < 0.0                        ||
          cl->TwoF < 0.0                     ||
          cl->Alpha <         0.0 - epsilon  ||
          cl->Alpha >   LAL_TWOPI + epsilon  ||
          cl->Delta < -0.5*LAL_PI - epsilon  ||
          cl->Delta >  0.5*LAL_PI + epsilon  ||
          !isfinite(cl->FileID)              ||
          !isfinite(cl->f)                   ||
          !isfinite(cl->Alpha)               ||
          !isfinite(cl->F1dot)               ||
          !isfinite(cl->Delta)               ||
          !isfinite(cl->TwoF)
          ) {
          XLALPrintError(
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
	  ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }
           
           

      /* check that the FIRST character following the Fstat value is a
         newline.  Note deliberate LACK OF WHITE SPACE char before %c
         above */
      if (newline != '\n') {
        XLALPrintError(
                "Line %d of file %s had extra chars after F value and before newline.\n"
                "First 255 chars are:\n"
                "%s\n",
                i+1, fname, line1);
        LALFree ((*CList));
        fclose(fp);
	ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }

      /* check that we read 7 quantities with exactly the right format */
      if ( nread != 7 )
        {
          XLALPrintError ("Found %d not %d values on line %d in file '%s'\n"
                         "Line in question is\n%s",
                         nread, 7, i+1, fname, line1);               
          LALFree ((*CList));
          fclose(fp);
	  ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
        }



      i++;
    } /*  end of main while loop */
  /* check that we read ALL lines! */
  if (i != numlines) {
    XLALPrintError(
            "Read of file %s terminated after %d line but numlines=%d\n",
            fname, i, numlines);
    LALFree((*CList));
    fclose(fp);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* read final line with %DONE\n marker */
  if (!fgets(line1, sizeof(line1), fp)) {
    XLALPrintError(
            "Failed to find marker line of file %s\n",
            fname);
    LALFree((*CList));
    fclose(fp);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* check for %DONE\n marker */
  if (strcmp(line1, DONE_MARKER)) {
    XLALPrintError(
            "Failed to parse marker: 'final' line of file %s contained %s not %s",
            fname, line1, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* check that we are now at the end-of-file */
  if (fgetc(fp) != EOF) {
    XLALPrintError(
            "File %s did not terminate after %s",
            fname, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* -- close candidate file -- */
  fclose(fp);     


  (*candlen) += numlines; /* total number of candidate so far */


  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);

} /* ReadOneCandidateFileV2() */





/* ########################################################################################## */
/*!
  Read one candidate-events file and fill CandidateList structure \b CList. 
  Count the number of the candidate events and fill it in \b candlen.
  The function aborts almost all the cases when the checks below failed.
  This function checks 
  \li if the file \b fname is readable. 
  \li if the file \b fname has the correct ending tag "%DONE".
  \li if the ranges of the values in the file are sensible.
  \li if the number of each row of the file is correct.
  \li if we could read all the events in the file.

  This function prints the bytecounts and the checksum of the file \b fname.
*/
void
ReadOneCandidateFile( LALStatus *lalStatus, 	/**< LALStatus pointer */
		      CandidateList **CList, 	/**< CandidateList str to be filled in this code  */
		      const CHAR *fname, 	/**< the name of the file to be read */
		      INT8 *candlen, 		/**< total number of the candidate events */
                      const REAL8 myFthr	/**< FIXME: !TO BE DOCUMENTED! */
                      )
{
  INT8 i;
  INT8 numlines;
  REAL8 epsilon=1e-5;
  CHAR line1[256];
  FILE *fp;
  INT8 nread;
  INT8 checksum=0;
  INT8 bytecount=0;
  INT8 numlinesFthr=0;


  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);
  ASSERT( fname != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( *CList == NULL, lalStatus, POLKAC_ENONULL, POLKAC_MSGENONULL);

  /* ------ Open and count candidates file ------ */
  i=0;
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      XLALPrintError("File %s doesn't exist!\n",fname);
      ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
     }

  while(fgets(line1,sizeof(line1),fp) != NULL) {
    UINT8 k;
    size_t len=strlen(line1);

    /* check that each line ends with a newline char (no overflow of
       line1 or null chars read) */
    if (!len || line1[len-1] != '\n') {
      XLALPrintError(
              "Line %d of file %s is too long or has no NEWLINE.  First 255 chars are:\n%s\n",
              i+1, fname, line1);
      fclose(fp);
      ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
     }

    /* increment line counter */
    i++;
    
    /* maintain a running checksum and byte count */
    bytecount+=len;
    for (k=0; k<len; k++)
      checksum+=(INT8)line1[k];
  }
  numlines=i;
  /* -- close candidate file -- */
  fclose(fp);     

  if ( numlines == 0) 
    {
      XLALPrintError ("ERROR: File '%s' has no lines so is not properly terminated by: %s", fname, DONE_MARKER);
      ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
    }

  /* output a record of the running checksun and byte count */
  fprintf(stdout, "%% %s: bytecount %" LAL_INT8_FORMAT " checksum %" LAL_INT8_FORMAT "\n", fname, bytecount, checksum);

  /* check validity of this Fstats-file */
  if ( strcmp(line1, DONE_MARKER ) ) 
    {
      XLALPrintError ("ERROR: File '%s' is not properly terminated by: %sbut has %s instead", fname, DONE_MARKER, line1);
      ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
    }
  else
    numlines --;        /* avoid stepping on DONE-marker */

  *candlen=numlines;

#if 0 /* Do we need to check this? */
  if (*candlen <= 0  )
    {
      XLALPrintError("candidate length = %ud!\n",*candlen);
      exit(POLKA_EXIT_ERR);;
    }/* check that we have candidates. */
#endif

  
  /* reserve memory for fstats-file contents */
  if (numlines > 0) 
    { 
      *CList = (CandidateList *)LALMalloc (numlines*sizeof(CandidateList));
      if ( !CList ) 
        { 
          XLALPrintError ("Could not allocate memory for candidate file %s\n\n", fname);
	  ABORT (lalStatus, POLKAC_EMEM, POLKAC_MSGEMEM);
        }
    }

  /* ------ Open and count candidates file ------ */
  i=0;
  fp=fopen(fname,"rb");
  if (fp==NULL) 
    {
      XLALPrintError("fopen(%s) failed!\n", fname);
      LALFree ((*CList));
      ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
    }
  while(i < numlines && fgets(line1,sizeof(line1),fp))
    {
      CHAR newline='\0';
      CandidateList *cl=&(*CList)[i];

      if (strlen(line1)==0 || line1[strlen(line1)-1] != '\n') {
        XLALPrintError(
                "Line %d of file %s is too long or has no NEWLINE.  First 255 chars are:\n%s\n",
                i+1, fname, line1);
        LALFree ((*CList));
        fclose(fp);
	ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }
      
      nread = sscanf (line1, 
                     "%" LAL_INT4_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL8_FORMAT 
                     " %" LAL_REAL4_FORMAT "%c", 
                     &(cl->FileID), &(cl->f), &(cl->Alpha), &(cl->Delta), &(cl->F1dot), &(cl->TwoF), &newline );

           
      /* find number of candidates that are above the 2F threshold. */
      if ( cl->TwoF > myFthr ) {
	numlinesFthr++;
      }
      
#ifdef SKYTESTMODE
      /* Pick a random values for the 2F-values, so that for cells covering two or more sky-points from
       the same data segment, each time pick a different sky-point from with in the cell. */
      cl->TwoF = (((float)rand()/RAND_MAX)*999)+1;
      if ( cl->TwoF <= 0.0 ) {
        cl->TwoF = EPSEDGE;
      }
#endif

      /* check that values that are read in are sensible */
      if (
          cl->FileID < 0                     ||
          cl->f < 0.0                        ||
          cl->TwoF < 0.0                     ||
          cl->Alpha <         0.0 - epsilon  ||
          cl->Alpha >   LAL_TWOPI + epsilon  ||
          cl->Delta < -0.5*LAL_PI - epsilon  ||
          cl->Delta >  0.5*LAL_PI + epsilon  ||
          !isfinite(cl->FileID)              ||
          !isfinite(cl->f)                   ||
          !isfinite(cl->Alpha)               ||
          !isfinite(cl->Delta)               ||
          !isfinite(cl->F1dot)               ||
          !isfinite(cl->TwoF)
          ) {
          XLALPrintError(
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
	  ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }
           
           

      /* check that the FIRST character following the Fstat value is a
         newline.  Note deliberate LACK OF WHITE SPACE char before %c
         above */
      if (newline != '\n') {
        XLALPrintError(
                "Line %d of file %s had extra chars after F value and before newline.\n"
                "First 255 chars are:\n"
                "%s\n",
                i+1, fname, line1);
        LALFree ((*CList));
        fclose(fp);
	ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
      }

      /* check that we read 7 quantities with exactly the right format */
      if ( nread != 7 )
        {
          XLALPrintError ("Found %d not %d values on line %d in file '%s'\n"
                         "Line in question is\n%s",
                         nread, 7, i+1, fname, line1);               
          LALFree ((*CList));
          fclose(fp);
	  ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
        }

      i++;
    } /*  end of main while loop */

  /* number of candidates above the 2F threshold. */
  /* *candilenFthr = numlinesFthr; */

  /* check that we read ALL lines! */
  if (i != numlines) {
    XLALPrintError(
            "Read of file %s terminated after %d line but numlines=%d\n",
            fname, i, numlines);
    LALFree((*CList));
    fclose(fp);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* read final line with %DONE\n marker */
  if (!fgets(line1, sizeof(line1), fp)) {
    XLALPrintError(
            "Failed to find marker line of file %s\n",
            fname);
    LALFree((*CList));
    fclose(fp);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* check for %DONE\n marker */
  if (strcmp(line1, DONE_MARKER)) {
    XLALPrintError(
            "Failed to parse marker: 'final' line of file %s contained %s not %s",
            fname, line1, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* check that we are now at the end-of-file */
  if (fgetc(fp) != EOF) {
    XLALPrintError(
            "File %s did not terminate after %s",
            fname, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    ABORT (lalStatus, POLKAC_EINVALIDFSTATS, POLKAC_MSGEINVALIDFSTATS);
  }

  /* -- close candidate file -- */
  fclose(fp);     

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);

} /* ReadOneCandidateFile() */







/* ########################################################################################## */
/*!
  Read command line arguments and fill a PolkaConfigVars structure \b CLA. 
  Almost all failures in this code invoke an exit of the code.
  This function checks 
  \li if a user has glob and specify \b input data dir.
  \li if a user specified either \b input data file or \b input data dir but not both.
  \li if a user did not define USE_UNZIP but specify \b input data dir.

  @param[in,out] lalStatus LALStatus* 
  @param[in]     argc      INT4  
  @param[in]     argv[]    CHAR* 
  @param[out]    CLA       PolkaConfigVars* Configuration variables
*/
void 
ReadCommandLineArgs( LALStatus *lalStatus, 
		     INT4 argc, 
		     CHAR *argv[], 
		     PolkaConfigVars *CLA ) 
{

  CHAR* uvar_InputData;
  CHAR* uvar_OutputData;

  CHAR* uvar_InputDirectory;
  CHAR* uvar_BaseName;

  CHAR* uvar_EahRun;

  BOOLEAN uvar_AutoOut;
  INT4 uvar_Nthr;  
  INT4 uvar_MinCellCoin;
  INT4 uvar_CellGrid;
  REAL8 uvar_Sthr;      
  REAL8 uvar_TwoFthr;

  REAL8 uvar_FreqWindow;
  REAL8 uvar_AlphaWindow;
  REAL8 uvar_DeltaWindow;
  REAL8 uvar_Kappa;
  REAL8 uvar_F1dotWindow;
  REAL8 uvar_FreqShift;
  REAL8 uvar_AlphaShift;
  REAL8 uvar_DeltaShift;
  REAL8 uvar_F1dotShift;
  BOOLEAN uvar_help;
  BOOLEAN uvar_UseUnzip;

  const CHAR BNAME[] = "Test";

  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  ASSERT( CLA != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);


  uvar_AutoOut = 0;
  uvar_help = 0;

  uvar_InputData = NULL;
  uvar_OutputData = NULL;

  uvar_InputDirectory = NULL;

  uvar_BaseName = (CHAR*)LALCalloc (1, strlen(BNAME)+1);
  strcpy (uvar_BaseName, BNAME);

  uvar_EahRun = (CHAR*)LALCalloc (1, strlen(BNAME)+1);
  strcpy (uvar_EahRun, BNAME);

  /* The following numbers are arbitrary. */
  uvar_MinCellCoin = 0;
  uvar_Nthr = 65536;     
  uvar_Sthr = 1.0e5; 
  uvar_TwoFthr = 0.0;

  uvar_FreqWindow = 0.0;
  uvar_AlphaWindow = 0.0;
  uvar_DeltaWindow = 0.0;
  uvar_F1dotWindow = 0.0;
  uvar_Kappa = 1;
  uvar_CellGrid = 0;
  
  uvar_FreqShift = 0.0;
  uvar_AlphaShift = 0.0;
  uvar_DeltaShift = 0.0;
  uvar_F1dotShift = 0.0;

  uvar_UseUnzip = 0;


  /* register all our user-variables */
  LALregBOOLUserVar(lalStatus,       help,           'h', UVAR_HELP,     "Print this message"); 
  LALregBOOLUserVar(lalStatus,       UseUnzip,       'z', UVAR_OPTIONAL, "Use Unzip"); 

  LALregSTRINGUserVar(lalStatus,     OutputData,     'o', UVAR_REQUIRED, "Ouput candidates file name");

  LALregSTRINGUserVar(lalStatus,     InputData,      'I', UVAR_OPTIONAL, "Input candidates Fstats file.");
  LALregSTRINGUserVar(lalStatus,     InputDirectory, 'i', UVAR_OPTIONAL, "Input candidates Fstats files directory.");
  LALregSTRINGUserVar(lalStatus,     BaseName,       'b', UVAR_OPTIONAL, "BaseName of the Input Fstats files");

  LALregINTUserVar(lalStatus,        Nthr,            0,  UVAR_OPTIONAL, "Threshold on number of coincidences");
  LALregREALUserVar(lalStatus,       Sthr,            0,  UVAR_OPTIONAL, "Threshold on significance.");
  LALregBOOLUserVar(lalStatus,       AutoOut,         0,  UVAR_OPTIONAL, "Set Nthr and Sthr to print most significant cell only."); 

  LALregREALUserVar(lalStatus,       TwoFthr,         0,  UVAR_OPTIONAL, "Threshold on TwoF values for candidates.");
  LALregINTUserVar(lalStatus,        CellGrid,       'g', UVAR_OPTIONAL, "Select CellGrid to use ( from 0 to 15 ).");

  LALregREALUserVar(lalStatus,       FreqWindow,     'f', UVAR_REQUIRED, "Frequency window in Hz");
  LALregREALUserVar(lalStatus,       F1dotWindow,    's', UVAR_REQUIRED, "First Spindown parameter window");
  LALregREALUserVar(lalStatus,       AlphaWindow,    'a', UVAR_REQUIRED, "Right Ascension window in radians");
  LALregREALUserVar(lalStatus,       DeltaWindow,    'd', UVAR_REQUIRED, "Declination window in radians");
  LALregREALUserVar(lalStatus,       Kappa,          'k', UVAR_OPTIONAL, "Tuning parameter for declination window");

  LALregREALUserVar(lalStatus,       FreqShift,      'F', UVAR_OPTIONAL, "Frequency shift in FreqWindow");
  LALregREALUserVar(lalStatus,       F1dotShift,     'S', UVAR_OPTIONAL, "First Spindown shift in F1dotWindow");
  LALregREALUserVar(lalStatus,       AlphaShift,     'A', UVAR_OPTIONAL, "Right Ascension shift in AlphaWindow");
  LALregREALUserVar(lalStatus,       DeltaShift,     'D', UVAR_OPTIONAL, "Declination shift in DeltaWindow");

  LALregSTRINGUserVar(lalStatus,     EahRun,         'r', UVAR_OPTIONAL, "E@H identifying run label for ifo split-up (S4R2a, S5R1a, Nautilus)");
  LALregINTUserVar(lalStatus,        MinCellCoin,     0,  UVAR_OPTIONAL, "Output all cells  with Coin's > MinCellCoin into separate files");

  TRY (LALUserVarReadAllInput(lalStatus->statusPtr,argc,argv),lalStatus); 


  if (uvar_help) {	/* if help was requested, we're done here */
    XLALPrintError("%s\n","$Id$");
    fflush(stderr);
    LALDestroyUserVars(lalStatus->statusPtr);
    exit(POLKA_EXIT_OK);
  }


  if( LALUserVarWasSet (&uvar_InputData) && 
      LALUserVarWasSet (&uvar_InputDirectory) ) {
    XLALPrintError("\nCannot set both of InputData and InputDirectory\n");
    exit(POLKA_EXIT_ERR);
  }

  if( (!LALUserVarWasSet (&uvar_InputData)) && 
      (!LALUserVarWasSet (&uvar_InputDirectory)) ) {
    XLALPrintError("\nPlease set either InputData and InputDirectory\n");
    exit(POLKA_EXIT_ERR);
  }


  if( uvar_UseUnzip ) {
#ifndef USE_UNZIP
    XLALPrintError("\n unzip can be used only when compiling with unzip enabled.\n");
    exit(POLKA_EXIT_ERR);
#endif
  }

  CLA->UseUnzip = uvar_UseUnzip;
  CLA->CellGrid = uvar_CellGrid;

  CLA->FstatsFile = NULL;
  CLA->OutputFile = NULL;
  CLA->InputDir = NULL;
  CLA->BaseName = NULL;
  CLA->EahRun = NULL;

  if( LALUserVarWasSet (&uvar_InputData) ) {
    CLA->FstatsFile = (CHAR *) LALMalloc(strlen(uvar_InputData)+1);
    if(CLA->FstatsFile == NULL)
      {
	XLALPrintError("No candidates file specified; input with -I option.\n");
	XLALPrintError("For help type %s -h\n", argv[0]);
	exit(POLKA_EXIT_ERR);
      }      
    strcpy(CLA->FstatsFile,uvar_InputData);
  }

  CLA->OutputFile = (CHAR *) LALMalloc(strlen(uvar_OutputData)+1);
  if(CLA->OutputFile == NULL)
    {
      TRY( FreeConfigVars( lalStatus->statusPtr, CLA ), lalStatus);
      exit(POLKA_EXIT_ERR);
    }      

  strcpy(CLA->OutputFile,uvar_OutputData);

  if( LALUserVarWasSet (&uvar_InputDirectory) ) {
#ifndef HAVE_GLOB_H   
    XLALPrintError("Sorry, but you cannot use this feature without glob.h.\n");
    exit(POLKA_EXIT_ERR);
#endif
    CLA->InputDir = (CHAR *) LALMalloc(strlen(uvar_InputDirectory)+1);
    if(CLA->InputDir == NULL)
      {
	TRY( FreeConfigVars( lalStatus->statusPtr, CLA ), lalStatus);
	exit(POLKA_EXIT_ERR);
      }          
    strcpy(CLA->InputDir,uvar_InputDirectory);
  }



  CLA->BaseName = (CHAR *) LALMalloc(strlen(uvar_BaseName)+1);
  if(CLA->BaseName == NULL)
    {
      TRY( FreeConfigVars( lalStatus->statusPtr, CLA ), lalStatus);
      exit(POLKA_EXIT_ERR);
    }          
  strcpy(CLA->BaseName,uvar_BaseName);
  
  CLA->EahRun = (CHAR *) LALMalloc(strlen(uvar_EahRun)+1);
  if(CLA->EahRun == NULL)
    {
      TRY( FreeConfigVars( lalStatus->statusPtr, CLA ), lalStatus);
      exit(POLKA_EXIT_ERR);
    }
  strcpy(CLA->EahRun,uvar_EahRun);


  CLA->AutoOut = uvar_AutoOut;
  CLA->MinCellCoin = uvar_MinCellCoin;
  CLA->Nthr = uvar_Nthr;
  CLA->Sthr = uvar_Sthr;
  CLA->TwoFthr = uvar_TwoFthr;

  CLA->Deltaf = uvar_FreqWindow;
  CLA->DeltaAlpha = uvar_AlphaWindow;
  CLA->DeltaDelta = uvar_DeltaWindow;
  CLA->DeltaF1dot = uvar_F1dotWindow;
  CLA->Kappa = uvar_Kappa;

  CLA->Shiftf = uvar_FreqShift;
  CLA->ShiftAlpha = uvar_AlphaShift;
  CLA->ShiftDelta = uvar_DeltaShift;
  CLA->ShiftF1dot = uvar_F1dotShift;

  LALDestroyUserVars(lalStatus->statusPtr);
  BEGINFAIL(lalStatus) {
    LALFree(CLA->FstatsFile);
    LALFree(CLA->OutputFile);
  } ENDFAIL(lalStatus);

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
} /* void ReadCommandLineArgs()  */



#if 0
/* ################################################################### */
/*
  Bruce's  sorting function to replace qsort. 
*/

void sortCandidates(INT8 *data, INT8 N)
{
  INT8 i, j, z;
  INT8 v, t;
  
  if(N<=1) return;
  
  /* 
     Partition elements
     BRUCE: first swap data[0] and data[rand(0,N-1)]
  */
  z=(INT8) ((N-1) * (rand() / (RAND_MAX+1.0)));
  t = data[0]; 
  data[0] = data[z]; 
  data[z] = t;


  v = data[0];
  i = 0;
  j = N;

  for(;;)
  {
    while(compareCandidatesList(&v, &data[++i]) && (i < N) ) { }  /*  BRUCE: replace data[]<v with compare function */
    while(compareCandidatesList(&data[--j], &v) ) { }   /* BRUCE: replace data[]>v with compare function */
    if(i >= j) break;
    t = data[i]; 
    data[i] = data[j]; 
    data[j] = t;
  }

  t = data[i-1]; 
  data[i-1] = data[0]; 
  data[0] = t;
  
 
  sortCandidates(data, i-1);
  sortCandidates(data+i, N-i);
}
#endif

/* ################################################################### */
/* ################################################################### */
/*
  Holger's sorting function to replace qsort. 

  Sorting function to sort cell indices in the INCREASING order of f, delta, alpha, FileID and 
  DECREASING ORDER OF a significance.
*/

void sortCandidates2(INT8 *data, INT8 left, INT8 right)
{
  INT8 i, last, t, rzahl=0;
  
  if (left >= right) 
    return;

  rzahl = (INT8) ( rand() % (right-left) );
  t =  data[left+rzahl];
  data[left+rzahl] = data[left];
  data[left] = t;
  
  last=left;

  for (i = left + 1; i <= right; i++)
  {
    if(SortedC[data[i]].iFreq < SortedC[data[left]].iFreq) {  
      last++;
      t =  data[last]; 
      data[last] = data[i]; 
      data[i] = t;
    }
    else if (SortedC[data[i]].iFreq == SortedC[data[left]].iFreq) {
      if (SortedC[data[i]].iDelta < SortedC[data[left]].iDelta) {
	last++;
	t =  data[last]; 
	data[last] = data[i]; 
	data[i] = t;
      }
      else if (SortedC[data[i]].iDelta == SortedC[data[left]].iDelta) {
	if (SortedC[data[i]].iAlpha < SortedC[data[left]].iAlpha) {
	  last++;
	  t =  data[last]; 
	  data[last] = data[i]; 
	  data[i] = t;
	}
	else if (SortedC[data[i]].iAlpha == SortedC[data[left]].iAlpha) {
	  if (SortedC[data[i]].iF1dot < SortedC[data[left]].iF1dot) {
	    last++;
	    t =  data[last]; 
	    data[last] = data[i]; 
	    data[i] = t;
	  }
	  else if (SortedC[data[i]].iF1dot == SortedC[data[left]].iF1dot) {
	    if (SortedC[data[i]].FileID < SortedC[data[left]].FileID) {
	      last++;
	      t =  data[last]; 
	      data[last] = data[i]; 
	      data[i] = t;
	    }
	    else if (SortedC[data[i]].FileID == SortedC[data[left]].FileID) {
	      if (SortedC[data[i]].TwoF > SortedC[data[left]].TwoF) {
		last++;
		t =  data[last]; 
		data[last] = data[i]; 
		data[i] = t;
	      }
	    }
	  }
	}
      }
    }
  }

  t = data[left]; 
  data[left] = data[last]; 
  data[last] = t;

  sortCandidates2(data, left, last-1);
  sortCandidates2(data, last+1, right);
}


/* ################################################################### */
/* ################################################################### */
/*
  Holger's sorting function to replace qsort. 
*/

void sortFreqCells2(INT8 *data, INT8 left, INT8 right)
{
  INT8 i, last, t, rzahl=0;
  
  if (left >= right) 
    return;

  rzahl = (INT8) ( rand() % (right-left) );
  t =  data[left+rzahl];
  data[left+rzahl] = data[left];
  data[left] = t;
  
  last=left;

  for (i = left + 1; i <= right; i++)
  {
    if(global_cell[data[i]].iFreq < global_cell[data[left]].iFreq)  {
      last++;
      t =  data[last]; 
      data[last] = data[i]; 
      data[i] = t;
    }
    else if (global_cell[data[i]].iFreq == global_cell[data[left]].iFreq) {
      if(global_cell[data[i]].nCand > global_cell[data[left]].nCand)  {
	last++;
	t =  data[last]; 
	data[last] = data[i]; 
	data[i] = t;
      }
      else if (global_cell[data[i]].nCand == global_cell[data[left]].nCand) {
	 if (global_cell[data[i]].significance > global_cell[data[left]].significance)  {
	   last++;
	   t =  data[last]; 
	   data[last] = data[i]; 
	   data[i] = t;
	 }
      }
    }
  } 

  t = data[left]; 
  data[left] = data[last]; 
  data[last] = t;

  sortFreqCells2(data, left, last-1);
  sortFreqCells2(data, last+1, right);
}
/* ########################################################################################## */



/* ########################################################################################## */
void print_info_of_cell_and_ifo_S4R2a( LALStatus *lalStatus,
                              FILE *fp,
                              const CellData *cd,
                              const CandidateList *CList,
                              const INT8 icell_start,
                              const INT8 icell_end,
                              const REAL8 sig_thr,
                              const REAL8 ncand_thr )
{
  INT8 idx, ic, icell;
  INT4 cH1,cL1;
  struct int8_linked_list *p;

  INITSTATUS(lalStatus);
  ASSERT( cd != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( CList != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);

  cH1 = 0;
  cL1 = 0;

  icell = icell_start;
  while( icell < icell_end &&
         cd[icell].significance > sig_thr &&
         cd[icell].nCand > ncand_thr )
    {

      cH1 = 0;
      cL1 = 0;
      p = cd[icell].CandID;
      ic = 0;
      while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) {
        idx = p->data;

	switch( CList[idx].FileID ) 
	  {
	  case 6537:
	  case 6497:
	  case 5828:
	  case 6120:
	  case 5955:
	  case 5613:
	  case 6126:
	  case 5946:
	  case 6130:
	  case 5515:
	    cH1++;
	    break;
	  case 6341:
	  case 6102:
	  case 5813:
	  case 5783:
	  case 5538:
	  case 6514:
	  case 5653:
	    cL1++;
	    break;
	  }
	
        p = p->next;
        ic++;
      } /*   while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) {  */

      if( ic >  LINKEDSTR_MAX_DEPTH ) {
        XLALPrintError("Maximum depth of linked structure reached!");
        exit(POLKA_EXIT_ERR);
      }

      if( cd[icell].nCand != (cH1+cL1) ) {
        XLALPrintError("Split-up of number of coincidences among detectors incorrect!");
        exit(POLKA_EXIT_ERR);
      }


      fprintf(fp,"%.6f\t%.6f\t%.6f\t%g\t\t%d\t%.6f\t%d\t%d\n",\
	      cd[icell].Freq, cd[icell].Delta, cd[icell].Alpha, cd[icell].F1dot, cd[icell].nCand, cd[icell].significance,cH1,cL1);
      icell++;

    } /*   while( icell < icell_end && ...  */

  RETURN (lalStatus);
} /* void print_info_of_cell_and_ifo_S4R2a( ) */


/* ########################################################################################## */
void print_info_of_cell_and_ifo_S5R1a( LALStatus *lalStatus,
				       FILE *fp,
				       const CellData *cd,
				       const CandidateList *CList,
				       const INT8 icell_start,
				       const INT8 icell_end,
				       const REAL8 sig_thr,
				       const REAL8 ncand_thr )
{
  INT8 idx, ic, icell;
  INT4 cH1,cL1;
  struct int8_linked_list *p;

  INITSTATUS(lalStatus);
  ASSERT( cd != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);
  ASSERT( CList != NULL, lalStatus, POLKAC_ENULL, POLKAC_MSGENULL);

  cH1 = 0;
  cL1 = 0;

  icell = icell_start;
  while( icell < icell_end &&
         cd[icell].significance > sig_thr &&
         cd[icell].nCand > ncand_thr )
    {

      cH1 = 0;
      cL1 = 0;
      p = cd[icell].CandID;
      ic = 0;
      while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) {
        idx = p->data;

        switch( CList[idx].FileID )
          {
	  case 5800:
	  case 6065:
	  case 6066:
	  case 5923:
	  case 5728:
	  case 5684:
	  case 5883:
	  case 5719:
	  case 6074:
	  case 5909:
	  case 6069:
	  case 5726:
	  case 6780:
	  case 7049:
	  case 6026:
	  case 5702:
	  case 5905:
	  case 6502:
	  case 6493:
	  case 5790:
	  case 5949:
	  case 5768:
	  case 0:
          case 1:
          case 2:
          case 3:
          case 4:
          case 5:
          case 6:
          case 7:
          case 8:
          case 9:
          case 10:
          case 11:
          case 12:
          case 13:
          case 14:
          case 15:
          case 16:
          case 17:
          case 18:
          case 19:
          case 20:
          case 21:
	    cH1++;
            break;

	  case 6284:
	  case 5843:
	  case 6500:
	  case 5925:
	  case 6744:
	  case 5835:
          case 22:
          case 23:
          case 24:
          case 25:
          case 26:
          case 27:
            cL1++;
            break;
          }

        p = p->next;
        ic++;
      } /*   while( p !=NULL && ic <= LINKEDSTR_MAX_DEPTH ) {  */
      
      if( ic >  LINKEDSTR_MAX_DEPTH ) {
	XLALPrintError("Maximum depth of linked structure reached!");
	exit(POLKA_EXIT_ERR);
      }
      
      if( cd[icell].nCand != (cH1+cL1) ) {
	XLALPrintError("Split-up of number of coincidences among detectors incorrect!");
	exit(POLKA_EXIT_ERR);
	}
      
      fprintf(fp,"%.6f\t%.6f\t%.6f\t%g\t\t%d\t%.6f\t%d\t%d\n",\
	      cd[icell].Freq, cd[icell].Delta, cd[icell].Alpha, cd[icell].F1dot, cd[icell].nCand, cd[icell].significance,cH1,cL1);
      icell++;
      
    } /*   while( icell < icell_end && ...  */
  
  RETURN (lalStatus);
} /* void print_info_of_cell_and_ifo_S5R1a( ) */

