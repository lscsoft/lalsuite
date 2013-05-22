/*
 *  Copyright (C) 2007  Holger Pletsch
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
/*      FreqFiducialTime - shifting frequency parameters of Einstein at Home result files   */
/*              to a fixed fiducial time for post-processing coincidence analysis           */
/*                                                                                          */
/*                              Holger Pletsch, UWM - March 2006                            */ 
/********************************************************************************************/

/* This code simply shifts all frequency parameters of a combined result file generated
   by combiner_v2.py, to a fixed fiducial GPS time for later coincidence analysis. 
   Note: The code makes use of the current Einstein at Home setup file.
   (i.e. 'CFS_S4R2_setup.h')  
*/


/***************************************************************************************************************************/
/*   INSTRUCTIONS for building "FreqFiducialTime" under LALApps                                                            */
/*                                                                                                                         */
/*   Files from Einstein@Home Workunit Generator that are needed from:                                                     */
/*   http://www.lsc-group.phys.uwm.edu/cgi-bin/cvs/viewcvs.cgi/einsteinathome/CFS/workunit_generator/?cvsroot=lscsoft      */
/*                                                                                                                         */
/*   - CFS_S4R2_setup.C (version 1.4, 2006/07/21)                                                                          */
/*   - WU_generator_misc.C (version 1.2, 2006/06/08)                                                                       */
/*   - WU_generator_daemon.h (version 1.6, 2006/07/21)                                                                     */
/*                                                                                                                         */
/*   Add the following lines to Makefile.am:                                                                               */
/*                                                                                                                         */
/*   INSTANCE = a                                                                                                          */
/*   warnLevel = -Wall -W -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -fno-common */
/*   WUDIR = ${LALAPPS_SRC}/pulsar/FDS_isolated                                                                            */
/*   CFLAGS := $(CFLAGS) $(CXXFLAGS) -g -O3 $(warnLevel)                                                                   */
/*                                                                                                                         */
/*   bin_PROGRAMS = ... FreqFiducialTime                                                                                   */
/*                                                                                                                         */
/*   CFS_S4R2_setup_CFLAGS := $(CFLAGS) -DDEFAULT_INSTANCE=\"$(INSTANCE)\"                                                 */
/*   CFS_S4R2_setup_SOURCES = CFS_S4R2_setup.C WU_generator_daemon.h                                                       */
/*   WU_generator_misc_CFLAGS := $(CFLAGS) $(CXXFLAGS) -g $(warnLevel)                                                     */
/*   WU_generator_misc_SOURCES = WU_generator_misc.C WU_generator_daemon.h                                                 */
/*   FreqFiducialTime_IFLAGS = -I$(WUDIR)                                                                                  */
/*   FreqFiducialTime_CFLAGS = $(CXXFLAGS) -Wall -g                                                                        */
/*   FreqFiducialTime_SOURCES = FreqFiducialTime.c CFS_S4R2_setup.C WU_generator_misc.C WU_generator_daemon.h              */
/*                                                                                                                         */
/***************************************************************************************************************************/

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


/* ----------------------------------------------------------------------------- */
/* file includes */
#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "getopt.h"
#include <math.h>


#include <unistd.h>

#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/LALStatusMacros.h>
#include <lal/ConfigFile.h>
#include <lal/UserInput.h>

#include <lalapps.h>

#include "WU_generator_daemon.h"  /* current Einstein at Home setup */


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
#define FIDUCIALC_ENULL            1
#define FIDUCIALC_ENONULL          2
#define FIDUCIALC_ESYS             3
#define FIDUCIALC_EINVALIDFSTATS   4
#define FIDUCIALC_EMEM             5
#define FIDUCIALC_ENORMAL          6


#define FIDUCIALC_MSGENULL         "Arguments contained an unexpected null pointer"
#define FIDUCIALC_MSGENONULL       "Input pointer was not NULL"
#define FIDUCIALC_MSGESYS          "System call failed (probably file IO"
#define FIDUCIALC_MSGEINVALIDFSTATS "Invalid Fstats file"
#define FIDUCIALC_MSGEMEM          "Sorry, ran out of memory... bye."


#define FIDUCIAL_EXIT_OK 0
#define FIDUCIAL_EXIT_ERR     31
#define FIDUCIAL_EXIT_READCND  32
#define FIDUCIAL_EXIT_FCTEST   33
#define FIDUCIAL_EXIT_OUTFAIL  34


/* ----------------------------------------------------------------------------- */
/* structures */

typedef struct FiducialTimeConfigVarsTag 
{
  REAL8 ThrTwoF;
  INT4 InNumLines;
  INT4 FiducialTime; /*  The fiducial GPS time */
  CHAR *OutputFile;  /*  Name of output file */
  CHAR *InputFile;   /*  Name of input file (combined result file produced by combiner_v2.py */
} FiducialTimeConfigVars;

typedef struct CandidateListTag
{

  REAL8 f;           /*  Frequency of the candidate */
  REAL8 Alpha;       /*  right ascension of the candidate */
  REAL8 Delta;       /*  declination  of the candidate */
  REAL8 F1dot;       /*  spindown (d/dt f) of the candidate */
  REAL8 TwoF;        /*  Maximum value of F for the cluster */
  INT4 FileID;       /*  File ID to specify from which file the candidate under consideration originally came. */
  INT4 DataStretch;
  CHAR resultfname[256];  /*  Name of the particular result file where values originally came from. */
} CandidateList;     




/* ----------------------------------------------------------------------------- */
/* Function declarelations */
void ReadCommandLineArgs( LALStatus *, INT4 argc, CHAR *argv[], FiducialTimeConfigVars *CLA ); 
void ReadCombinedFile( LALStatus *lalStatus, CandidateList **CList, FiducialTimeConfigVars *CLA, long *candlen );

void ComputeFiducialTimeFrequency( LALStatus *,	FiducialTimeConfigVars *CLA, CandidateList *CList, INT4 candlen );

void PrintResultFile( LALStatus *, const FiducialTimeConfigVars *CLA, CandidateList *CList, long candlen );

void FreeMemory( LALStatus *, FiducialTimeConfigVars *CLA, CandidateList *CList, const UINT4 datalen );
void FreeConfigVars( LALStatus *, FiducialTimeConfigVars *CLA );




/* ----------------------------------------------------------------------------- */
/* Global Variables */
/*! @param global_status LALStatus Used to initialize LALStatus lalStatus. */
LALStatus global_status;
/*! @param lalDebugLevel INT4 Control debugging behaviours. Defined in lalapps.h */
/*! @param vrbflg        INT4 Control debugging messages. Defined in lalapps.h */
extern INT4 vrbflg;

const REAL8 FIXED_FIDUCIAL_TIME = 793555944; /* here e.g. GPS startTime of first data stretch in S4 is chosen */

/* ------------------------------------------------------------------------------------------*/
/* Code starts here.                                                                         */
/* ------------------------------------------------------------------------------------------*/
/* ########################################################################################## */
/*!
  Main function  
*/

int main(INT4 argc,CHAR *argv[]) 
{
  LALStatus *lalStatus = &global_status;
  long CLength=0;
  CandidateList *AllmyC = NULL;
  FiducialTimeConfigVars FTCV;


  vrbflg = 1;   /* verbose error-messages */

  /* Reads command line arguments */
  LAL_CALL( ReadCommandLineArgs( lalStatus, argc, argv, &FTCV ), lalStatus); 

  /* Reads in combined candidare file, set CLength */
  LAL_CALL( ReadCombinedFile(lalStatus, &AllmyC, &FTCV, &CLength), lalStatus);

  /* -----------------------------------------------------------------------------------------*/      
  /* Compute shifting of frequency parameters */
  LAL_CALL( ComputeFiducialTimeFrequency( lalStatus, &FTCV, AllmyC, CLength),lalStatus );
 
  /* -----------------------------------------------------------------------------------------*/      
  /* Output result file */
  LAL_CALL( PrintResultFile( lalStatus, &FTCV, AllmyC, CLength),lalStatus );

  /* -----------------------------------------------------------------------------------------*/      
  /* Clean-up */
  LAL_CALL( FreeMemory(lalStatus, &FTCV, AllmyC, CLength), lalStatus );

  LALCheckMemoryLeaks(); 

  return(FIDUCIAL_EXIT_OK);
 
} /* main() */


/* ########################################################################################## */

void PrintResultFile(LALStatus *lalStatus, const FiducialTimeConfigVars *CLA, CandidateList *CList, long candlen)
{
  
  INT4 iindex;
  FILE *fp = NULL;
  INT4 *count;
  INT4 nmax = 0;
  const CHAR *fname;
  
  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  ASSERT( CLA != NULL, lalStatus, FIDUCIALC_ENULL, FIDUCIALC_MSGENULL);
  ASSERT( CList != NULL, lalStatus, FIDUCIALC_ENULL, FIDUCIALC_MSGENULL);
  
  if( (count = (INT4 *) LALCalloc( (size_t) (nmax + 1), sizeof(INT4))) == NULL ) {
    XLALPrintError("Could not allocate Memory! \n");
    ABORT (lalStatus, FIDUCIALC_EMEM, FIDUCIALC_MSGEMEM);
  }
  
  fname=CLA->OutputFile;  


  /* ------------------------------------------------------------- */
  /* Print out to the user-specified output file.*/
  if(strcmp(fname,"-")==0){
    fp=stdout;
  }
  else {
    if( (fp = fopen(fname,"w")) == NULL ) 
      {
	XLALPrintError("\n Cannot open output file %s\n",CLA->OutputFile); 
	ABORT (lalStatus, FIDUCIALC_EMEM, FIDUCIALC_MSGEMEM);
      }
  }

  /* output lines */
  /*INITSTATUS(lalStatus); */
 
  iindex=0;
  
  while(iindex < candlen) 
    {
      fprintf(fp,"%" LAL_INT4_FORMAT " %.13g %.7g %.7g %.5g %.6g\n",
	      CList[iindex].DataStretch, 
	      CList[iindex].f, 
	      CList[iindex].Alpha, 
	      CList[iindex].Delta, 
	      CList[iindex].F1dot, 
	      CList[iindex].TwoF );
      iindex++;
    }

  fprintf(fp, "%s", DONE_MARKER);
  
  BEGINFAIL(lalStatus) {fclose(fp);} ENDFAIL(lalStatus);

  LALFree( count );

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
} /* PrintResult() */



/* ########################################################################################## */

void FreeMemory( LALStatus *lalStatus, 
	    FiducialTimeConfigVars *CLA, 
	    CandidateList *CList, 
	    const UINT4 CLength)
{
  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  FreeConfigVars( lalStatus->statusPtr, CLA );

  if( CList != NULL ) LALFree(CList);

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
} /* FreeMemory */


/* ########################################################################################## */

void FreeConfigVars(LALStatus *lalStatus, FiducialTimeConfigVars *CLA )
{
  INITSTATUS(lalStatus);

  if( CLA->OutputFile != NULL ) LALFree(CLA->OutputFile);
  if( CLA->InputFile != NULL ) LALFree(CLA->InputFile);

  RETURN (lalStatus);
} /* FreeCOnfigVars */


/* ########################################################################################## */

void ReadCombinedFile( LALStatus *lalStatus, 
		      CandidateList **CList, 
		      FiducialTimeConfigVars *CLA, 
		      long *candlen )
{
  long i,jj;
  INT4  numlines;
  REAL8 epsilon=1e-5;
  CHAR line1[256];
  FILE *fp;
  long nread;
  UINT4 checksum=0;
  UINT4 bytecount=0;
  INT4 sizelist=16384;
  const CHAR *fname;
        
  fname = CLA->InputFile;
 
  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);
  ASSERT( fname != NULL, lalStatus, FIDUCIALC_ENULL, FIDUCIALC_MSGENULL);
  ASSERT( *CList == NULL, lalStatus, FIDUCIALC_ENONULL, FIDUCIALC_MSGENONULL);

  /* ------ Open and count candidates file ------ */
  i=0;

  if(strcmp(fname,"-")==0){
    fp=stdin;
    }
  else {
    fp=fopen(fname,"rb");
    if (fp==NULL) 
      {
	XLALPrintError("File %s doesn't exist!\n",fname);
	ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
      }
  }

  
  if(strcmp(fname,"-")==0){
    numlines = CLA->InNumLines;
#if 0
    while(gets(line1)) {
      UINT4 k;
      size_t len=strlen(line1);
      
      /* increment line counter */
      i++;
      
      /* maintain a running checksum and byte count */
      bytecount+=len;
      for (k=0; k<len; k++)
	checksum+=(INT4)line1[k];
    }
#endif
  }
  
  else{
    while( fgets(line1,sizeof(line1),fp) ) {
      UINT4 k;
      size_t len=strlen(line1);
      
      /* check that each line ends with a newline char (no overflow of
	 line1 or null chars read) */
      if (!len || line1[len-1] != '\n') {
	XLALPrintError(
		      "Line %d of file %s is too long or has no NEWLINE.  First 255 chars are:\n%s\n",
		      i+1, fname, line1);
	fclose(fp);
	ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
      }
      
      /* increment line counter */
      i++;
      
      /* maintain a running checksum and byte count */
      bytecount+=len;
      for (k=0; k<len; k++)
	checksum+=(INT4)line1[k];
    }
    /* -- close candidate file -- */
    fclose(fp);  
    numlines=i;
  }
    
 
  
  if ( numlines == 0) 
    {
      XLALPrintError ("ERROR: File '%s' has no lines so is not properly terminated by: %s", fname, DONE_MARKER);
      ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
    }

  if(strcmp(fname,"-") != 0){
    /* output a record of the running checksun amd byte count */
    XLALPrintError( "%% %s: bytecount %" LAL_UINT4_FORMAT " checksum %" LAL_UINT4_FORMAT "\n", fname, bytecount, checksum);
    
    /* check validity of this Fstats-file */
    if ( strcmp(strcat(line1,"\n"), DONE_MARKER ) ) 
      {
	XLALPrintError ("ERROR: File '%s' is not properly terminated by: %sbut has %s instead", fname, DONE_MARKER, line1);
	ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
      }
    else
      numlines --;        /* avoid stepping on DONE-marker */

    *candlen=numlines;
  }


#if 0 /* Do we need to check this? */
  if (*candlen <= 0  )
    {
      XLALPrintError("candidate length = %ud!\n",*candlen);
      exit(FIDUCIAL_EXIT_ERR);;
    }/* check that we have candidates. */
#endif

  
  /* start reserving memory for fstats-file contents */
  if (numlines > 0) 
    { 
      *CList = (CandidateList *)LALMalloc (sizelist*sizeof(CandidateList));
      if ( !CList ) 
        { 
          XLALPrintError ("Could not allocate memory for candidate file %s\n\n", fname);
	  ABORT (lalStatus, FIDUCIALC_EMEM, FIDUCIALC_MSGEMEM);
        }
    }

  /* ------ Open and count candidates file ------ */
  i=0;
  if(strcmp(fname,"-")==0){
    fp=stdin;
  }
  else {
    fp=fopen(fname,"rb");
    if (fp==NULL) 
      {
	XLALPrintError("fopen(%s) failed!\n", fname);
	LALFree ((*CList));
	ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
      }
  }

  
  jj=0;

  /* printf("Number of lines to read: %d\n", numlines); */
    
  if(strcmp(fname,"-")==0){
    while(i < numlines && gets(line1) )
      {
       	CHAR newline='\0';
	
	if (jj >= sizelist) {
	  sizelist = sizelist + 16384;
	  CandidateList *tmp;
	  tmp = (CandidateList *)LALRealloc(*CList, (sizelist * sizeof(CandidateList)) );
	  if ( !tmp ){
	    XLALPrintError("couldnot re-allocate memory for candidate list \n\n");
	    ABORT (lalStatus, FIDUCIALC_EMEM, FIDUCIALC_MSGEMEM);
	  }
	  *CList = tmp;
	}
	
	CandidateList *cl=&(*CList)[jj];
	
	nread = sscanf (line1,
			"%s %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT, 
			&(cl->resultfname), &(cl->f), &(cl->Alpha), &(cl->Delta), &(cl->F1dot), &(cl->TwoF) );

	if(cl->TwoF >= CLA->ThrTwoF) {
	  /* check that values that are read in are sensible, 
	     (result file names will be checked later, when getting 
	     the search parameters from Einstein at Home  setup library) */
	  if (
	      cl->f < 0.0                        ||
	      cl->TwoF < 0.0                     ||
	      cl->Alpha <         0.0 - epsilon  ||
	      cl->Alpha >   LAL_TWOPI + epsilon  ||
	      cl->Delta < -0.5*LAL_PI - epsilon  ||
	      cl->Delta >  0.5*LAL_PI + epsilon  ||
	      !finite(cl->FileID)                ||                                                                 
	      !finite(cl->f)                     ||
	      !finite(cl->Alpha)                 ||
	      !finite(cl->Delta)                 ||
	      !finite(cl->F1dot)                 ||
	      !finite(cl->TwoF)
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
	    ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
	  }
	  
          

	  /* check that we read 6 quantities with exactly the right format */
	  if ( nread != 6 )
	    {
	      XLALPrintError ("Found %d not %d values on line %d in file '%s'\n"
			     "Line in question is\n%s",
			     nread, 6, i+1, fname, line1);               
	      LALFree ((*CList));
	      ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
	    }
	  
	  jj++;
	  
	} /* if(cl->TwoF >= CLA->ThrTwoF) */
	
	i++;
      } /*  end of main while loop */
  }
  else{
    while(i < numlines && fgets(line1,sizeof(line1),fp) )
      {
	CHAR newline='\0';
	
	if (jj >= sizelist) {
	  sizelist = sizelist + 16384;
	  CandidateList *tmp;
	  tmp = (CandidateList *)LALRealloc(*CList, (sizelist * sizeof(CandidateList)) );
	  if ( !tmp ){
	    XLALPrintError("couldnot re-allocate memory for candidate list \n\n");
	    ABORT (lalStatus, FIDUCIALC_EMEM, FIDUCIALC_MSGEMEM);
	  }
	  *CList = tmp;
	}
	
	CandidateList *cl=&(*CList)[jj];
	
	if (strlen(line1)==0 || line1[strlen(line1)-1] != '\n') {
	  XLALPrintError(
			"Line %d of file %s is too long or has no NEWLINE.  First 255 chars are:\n%s\n",
			i+1, fname, line1);
	  LALFree ((*CList));
	  fclose(fp);
	  ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
	}
	
	nread = sscanf (line1,
			"%s %" LAL_INT4_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" 
			LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT "%c", 
			&(cl->resultfname), &(cl->FileID), &(cl->f), &(cl->Alpha), &(cl->Delta), &(cl->F1dot), &(cl->TwoF), &newline );
	if(cl->TwoF >= CLA->ThrTwoF) {
	  /* check that values that are read in are sensible, 
	     (result file names will be checked later, when getting 
	     the search parameters from Einstein at Home  setup library) */
	  if (
	      cl->FileID < 0                     ||
	      cl->f < 0.0                        ||
	      cl->TwoF < 0.0                     ||
	      cl->Alpha <         0.0 - epsilon  ||
	      cl->Alpha >   LAL_TWOPI + epsilon  ||
	      cl->Delta < -0.5*LAL_PI - epsilon  ||
	      cl->Delta >  0.5*LAL_PI + epsilon  ||
	      !finite(cl->FileID)                ||                                                                 
	      !finite(cl->f)                     ||
	      !finite(cl->Alpha)                 ||
	      !finite(cl->Delta)                 ||
	      !finite(cl->F1dot)                 ||
	      !finite(cl->TwoF)
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
	    ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
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
	    ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
	  }
	  
	  /* check that we read 8 quantities with exactly the right format */
	  if ( nread != 8 )
	    {
	      XLALPrintError ("Found %d not %d values on line %d in file '%s'\n"
			     "Line in question is\n%s",
			     nread, 8, i+1, fname, line1);               
	      LALFree ((*CList));
	      fclose(fp);
	      ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
	    }
	  
	  jj++;
	  
	} /* if(cl->TwoF >= CLA->ThrTwoF) */
	
	i++;
      } /*  end of main while loop */
  }


  *candlen=jj-1;

  /* check that we read ALL lines! */
  if (i != numlines) {
    XLALPrintError(
            "Reading of file %s terminated after %d line but numlines=%d\n",
            fname, i, numlines);
    LALFree((*CList));
    fclose(fp);
    ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
  }

  /* read final line with %DONE\n marker */
  if (!fgets(line1, sizeof(line1), fp)) {
    XLALPrintError(
            "Failed to find marker line of file %s\n",
            fname);
    LALFree((*CList));
    fclose(fp);
    ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
  }

  /* check for %DONE\n marker */
  if (strcmp(line1, DONE_MARKER)) {
    XLALPrintError(
            "Failed to parse marker: 'final' line of file %s contained %s not %s",
            fname, line1, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
  }

  /* check that we are now at the end-of-file */
  if (fgetc(fp) != EOF) {
    XLALPrintError(
            "File %s did not terminate after %s",
            fname, DONE_MARKER);
    LALFree ((*CList));
    fclose(fp);
    ABORT (lalStatus, FIDUCIALC_EINVALIDFSTATS, FIDUCIALC_MSGEINVALIDFSTATS);
  }

  /* -- close candidate file -- */
  fclose(fp);     

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);

} /* ReadCombinedFile() */




/* ########################################################################################## */

void ReadCommandLineArgs( LALStatus *lalStatus, 
		     INT4 argc, 
		     CHAR *argv[], 
		     FiducialTimeConfigVars *CLA ) 
{

  CHAR* uvar_InputData;
  CHAR* uvar_OutputData;
  BOOLEAN uvar_help;
  REAL8 uvar_ThrTwoF;
  INT4 uvar_InNumLines;
  INT4 uvar_FiducialTime;

  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  ASSERT( CLA != NULL, lalStatus, FIDUCIALC_ENULL, FIDUCIALC_MSGENULL);

  uvar_help = 0;
  uvar_InputData = NULL;
  uvar_OutputData = NULL;
  uvar_ThrTwoF = -1.0;
  uvar_FiducialTime = FIXED_FIDUCIAL_TIME;
  
  /* register all our user-variables */
  LALregBOOLUserVar(lalStatus,       help,           'h', UVAR_HELP,     "Print this message"); 
  
  LALregINTUserVar(lalStatus,        FiducialTime,   'f', UVAR_OPTIONAL, "Fiducial GPS time");
  LALregSTRINGUserVar(lalStatus,     OutputData,     'o', UVAR_OPTIONAL, "Ouput file name");
  LALregSTRINGUserVar(lalStatus,     InputData,      'i', UVAR_OPTIONAL, "Input file name");
  LALregREALUserVar(lalStatus,       ThrTwoF,        't', UVAR_OPTIONAL, "Threshold on values of 2F");
  LALregINTUserVar(lalStatus,        InNumLines,     'l', UVAR_OPTIONAL, "Number of lines of input");


  TRY (LALUserVarReadAllInput(lalStatus->statusPtr,argc,argv),lalStatus); 


  if (uvar_help) {	/* if help was requested, we're done here */
    XLALPrintError("%s\n","$Id$");
    fflush(stderr);
    LALDestroyUserVars(lalStatus->statusPtr);
    exit(FIDUCIAL_EXIT_OK);
  }

  CLA->FiducialTime = 0;
  CLA->OutputFile = NULL;
  CLA->InputFile = NULL;
  CLA->ThrTwoF = -1.0;
  CLA->InNumLines = 0;


  if(uvar_OutputData != ""){
    if(uvar_OutputData !="-"){
      CLA->OutputFile = (CHAR *) LALMalloc(strlen(uvar_OutputData)+1);
      if(CLA->OutputFile == NULL) {
	TRY( FreeConfigVars( lalStatus->statusPtr, CLA ), lalStatus);
	XLALPrintError("Output file '%s' is incorrect. \n\n",uvar_OutputData);
	exit(FIDUCIAL_EXIT_ERR);
      }
    }
    strcpy(CLA->OutputFile,uvar_OutputData);
  }
  
  else 
    { 
      TRY( FreeConfigVars( lalStatus->statusPtr, CLA ), lalStatus);
      XLALPrintError("Output file '%s' is incorrect. \n\n",uvar_OutputData);
      exit(FIDUCIAL_EXIT_ERR);
    }

  
  if(uvar_InputData != ""){
    if(uvar_InputData !="-"){
      CLA->InputFile = (CHAR *) LALMalloc(strlen(uvar_InputData)+1);
      if(CLA->InputFile == NULL){
	TRY( FreeConfigVars( lalStatus->statusPtr, CLA ), lalStatus);
	XLALPrintError("Input file '%s' is incorrect. \n\n",uvar_InputData);
	exit(FIDUCIAL_EXIT_ERR);
      }
    }
    strcpy(CLA->InputFile,uvar_InputData);
  }
  else
    {
      TRY( FreeConfigVars( lalStatus->statusPtr, CLA ), lalStatus);
      XLALPrintError("Input file '%s' is incorrect. \n\n",uvar_InputData);
      exit(FIDUCIAL_EXIT_ERR);
    }

  if(uvar_FiducialTime > 0){
    CLA->FiducialTime = uvar_FiducialTime;
  }
  else    {
    TRY( FreeConfigVars( lalStatus->statusPtr, CLA ), lalStatus);
    XLALPrintError("The fiducial GPS time %d is incorrect. \n\n",uvar_FiducialTime);
    exit(FIDUCIAL_EXIT_ERR);
  }
  CLA->InNumLines = uvar_InNumLines;
  CLA->ThrTwoF = uvar_ThrTwoF;

  LALDestroyUserVars(lalStatus->statusPtr);
  BEGINFAIL(lalStatus) {
    LALFree(CLA->InputFile);
    LALFree(CLA->OutputFile);
  } ENDFAIL(lalStatus);

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);
} /* void ReadCommandLineArgs()  */



/* ########################################################################################## */


/*----------------------------------------------------- */
void ComputeFiducialTimeFrequency( LALStatus *lalStatus,
				   FiducialTimeConfigVars *CLA,
				   CandidateList *CList, 
				   INT4 candlen)
{
  REAL8 f_CFS;
  REAL8 F1dot_CFS;
  REAL8 f_fiducial;
  REAL8 deltaT;
  INT4 iindex;
  INT4 segtmp;
  WU_search_params_t wparams;


  INITSTATUS(lalStatus);
  ATTATCHSTATUSPTR (lalStatus);

  f_CFS=0;
  f_fiducial=0;
  F1dot_CFS=0;
  iindex=0;
  deltaT=0;
  
  /* Find data segment from resultname */
  findSearchParams4Result( CList[0].resultfname, &wparams );

  /* Compute Data-StretchID */
  segtmp = (INT4)(wparams.endTime / (wparams.endTime - wparams.startTime));
  
  switch( segtmp )
    {
      /* H1 ----------- */

    case 6537:
      segtmp=0;
      break;
      
    case 6497:
      segtmp=1;
      break;
      
    case 5828:
      segtmp=2;
      break;
      
    case 6120:
      segtmp=3;
      break;
      
    case 5955:
      segtmp=4;
      break;
      
    case 5613:
      segtmp=5;
      break;

    case 6126:
      segtmp=6;
      break;
      
    case 5946:
      segtmp=7;
      break;
      
    case 6130:
      segtmp=8;
      break;

    case 5515:
      segtmp=9;
      break;

      /* L1 ----------- */
      
    case 6341:
      segtmp=10;
      break;

    case 6102:
      segtmp=11;
      break;

    case 5813:
      segtmp=12;
      break;

    case 5783:
      segtmp=13;
      break;

    case 5538:
      segtmp=14;
      break;

    case 6514:
      segtmp=15;
      break;

    case 5653:
      segtmp=16;
      break;

    }
  
  while( iindex < candlen )
    {
      f_CFS = CList[iindex].f;
      F1dot_CFS = CList[iindex].F1dot;
      
      /* Get search parameters from Einstein at Home setup library */
      findSearchParams4Result( CList[iindex].resultfname, &wparams );

      /* Compute Data-StretchID */
      CList[iindex].DataStretch = segtmp;
      
      /* Fixed fiducial time = e.g. GPS time of first SFT in S4 */
      deltaT = wparams.startTime - CLA->FiducialTime;

      /* Compute new frequency values at fixed fiducial time */
      f_fiducial = f_CFS - (F1dot_CFS * deltaT);
     
      /* Replace f values by the new ones, that all refer to the same */
      CList[iindex].f = f_fiducial;
 
      f_CFS = 0;
      f_fiducial = 0;
      F1dot_CFS = 0;
      deltaT=0;

      iindex++;

    } /* while( iindex < candlen)  */

  DETATCHSTATUSPTR (lalStatus);
  RETURN (lalStatus);

} /* ComputeFiducialTimeFrequencies () */
