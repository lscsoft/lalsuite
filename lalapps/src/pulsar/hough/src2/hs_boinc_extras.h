/*
*  Copyright (C) 2007 Bernd Machenschalk
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

/* Extras for BOINC compilation of HierarchicalSearch
*/

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>

/* BOINC includes */
#include "filesys.h"

#include <lal/LALError.h>
#include <lal/LALRCSID.h>
#include "HoughFStatToplist.h"

#define HSBOINCEXTRASHRCSID "$Id$"

#define EAH_LOGLEVEL 1        /* LOG_DEBUG */
#define EAH_LALDEBUGLEVEL 33  /* DebugLevel = 1, but without time-consuming memory debugging */

/* linking proper functions to the hooks in HierarchicalSearch.c */

/* use a local copy of ComputeFStatFreqBand() and related functions for E@H-specific optimizations */
#ifdef EAH_CUDA
#include "HierarchicalSearch.h"
#include <cuda_fstat.h>
#define COMPUTEFSTATHOUGHMAP LocalComputeFstatHoughMap
#define REARRANGE_SFT_DATA   cuda_prepare_sfts (&stackMultiSFT, nStacks, &fstatVector)
#define COMPUTEFSTATFREQBAND(a,b,c,d,e,f,g) cuda_ComputeFStatFreqBand(a,b,c,d,e,f,g,k)
#else
#define REARRANGE_SFT_DATA
#ifndef EAH_OPTIMIZATION
#define COMPUTEFSTATHOUGHMAP ComputeFstatHoughMap
#define COMPUTEFSTATFREQBAND ComputeFStatFreqBand
#else
#define COMPUTEFSTATHOUGHMAP LocalComputeFstatHoughMap
#define COMPUTEFSTATFREQBAND LocalComputeFStatFreqBand

#include "HierarchicalSearch.h"

extern void
LocalComputeFStatFreqBand ( LALStatus *status, 
                            REAL8FrequencySeries *FstatVector,
                            const PulsarDopplerParams *doppler,
                            const MultiSFTVector *multiSFTs, 
                            const MultiNoiseWeights *multiWeights,
                            const MultiDetectorStateSeries *multiDetStates,
                            const ComputeFParams *params);
extern void
LocalComputeFstatHoughMap ( LALStatus *status,
			    SemiCohCandidateList  *out,   /* output candidates */
			    HOUGHPeakGramVector *pgV, /* peakgram vector */
			    SemiCoherentParams *params);
#endif
#endif

#define SHOW_PROGRESS show_progress
#define fopen boinc_fopen

#ifndef HS_CHECKPOINTING
#define HS_CHECKPOINTING 1
#endif

#define INSERT_INTO_HOUGHFSTAT_TOPLIST insert_into_houghFStat_toplist

#if (HS_CHECKPOINTING)
#define GET_CHECKPOINT(toplist,total,count,outputname,cptname)\
  { int ret = init_and_read_checkpoint(toplist,total,count,outputname,cptname);\
    if(ret < 0) {\
      LogPrintf(LOG_CRITICAL, HIERARCHICALSEARCH_MSGECHECKPT " (%d)\n",ret);\
      return(HIERARCHICALSEARCH_ECHECKPT);\
    } else if (ret == 2) {\
      return(HIERARCHICALSEARCH_ENORM);\
    }\
  }
#define SET_CHECKPOINT set_checkpoint()

#else
#define SET_CHECKPOINT
#define GET_CHECKPOINT(toplist,total,count,outputname,cptname) *total=0;
#endif

#ifdef  __cplusplus
extern "C" {
#endif

extern int global_cpu_type;

extern LALStatus *global_status;

/* function prototypes, they are defined in boinc_extras.c */

/** show progress of the App.
    NOTE: This also sets the count & total (skypos) for checkpointing */
extern void show_progress(REAL8 rac,   REAL8 dec,
			  REAL8 count, REAL8 total,
			  REAL8 freq,  REAL8 fband);

/** inits checkpointing for the toplist and reads the last checkpoint if present
    This expects all passed variables (toplist, total, count) to be already
    initialized. In case of an error, the toplist is cleared and the count
    is set to 0, total is effectivly ignored.
    If *cptname (name of the checkpoint file) is NULL,
    the name is constructed by appending ".cpt" to the output filename.

    The function returns
    0 if no checkpoint could be found,
   -1 if a checkpoint was found but it couldn't be read,
   -2 if an error occured (out of memory),
    1 if a checkpoint was found and previous output could be read
    2 if nothing to do (previously written output file was found)
*/
extern int init_and_read_checkpoint(toplist_t*toplist, UINT4*count,
                                     UINT4 total, char*outputname, char*cptname);

/** actually writes a checkpoint only if it's boinc_time_to_checkpoint() */
extern void set_checkpoint(void);

/** writes the toplist to the final (ASCII) output file */
extern void write_and_close_checkpointed_file (void);

/** LALApps error handler for BOINC */
extern int BOINC_LAL_ErrHand (LALStatus*, const char*, const char*, const int, volatile const char*);

/** attach gdb to the running process; for debugging. */
extern void attach_gdb(void);

/** play with floating-point exceptions */
extern void enable_floating_point_exceptions(void);

/** generate a segfault (for testing purposes) */
extern int segfault (void);

/** the main() function of HierarchicalSerach.c becomes the extern MAIN(),
    the actual main() function of the BOINC App is defined in boinc_extras.c
*/
extern int MAIN(int,char**);

#ifdef  __cplusplus
}
#endif
