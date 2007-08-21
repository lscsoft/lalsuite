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
   Author: Bernd Machenschalk
*/

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>

/* BOINC includes */
#include "filesys.h"

#include <lal/LALError.h>
#include <lal/LALRCSID.h>
NRCSID(HSBOINCEXTRASHRCSID,"$Id$");
#include "FstatToplist.h"

#define EAH_LOGLEVEL 1        /* LOG_DEBUG */
#define EAH_LALDEBUGLEVEL 33  /* DebugLevel = 1, but without time-consuming memory debugging */

/* linking proper functions to the hooks in HierarchicalSearch.c */

/* use a local copy of ComputeFStatFreqBand() and related functions for E@H-specific optimizations */
#ifndef EAH_OPTIMIZATION
#define COMPUTEFSTATFREQBAND ComputeFStatFreqBand
#else
#define COMPUTEFSTATFREQBAND LocalComputeFStatFreqBand

extern int
LocalXLALComputeFaFb ( Fcomponents *FaFb,
		  const SFTVector *sfts, 
		  const PulsarSpins fkdot,
		  const SSBtimes *tSSB,
		  const AMCoeffs *amcoe,
		  const ComputeFParams *params);

extern void
LocalComputeFStat ( LALStatus *, Fcomponents *Fstat, 
		    const PulsarDopplerParams *doppler,
		    const MultiSFTVector *multiSFTs,
		    const MultiNoiseWeights *multiWeights,
		    const MultiDetectorStateSeries *multiDetStates,
		    const ComputeFParams *params,
		    ComputeFBuffer *cfBuffer );

extern void
LocalComputeFStatFreqBand ( LALStatus *status, 
			    REAL8FrequencySeries *FstatVector,
			    const PulsarDopplerParams *doppler,
			    const MultiSFTVector *multiSFTs, 
			    const MultiNoiseWeights *multiWeights,
			    const MultiDetectorStateSeries *multiDetStates,
			    const ComputeFParams *params);
#endif

#define SHOW_PROGRESS show_progress
#define fopen boinc_fopen

#ifndef HS_CHECKPOINTING
#define HS_CHECKPOINTING 1
#endif

#if (HS_CHECKPOINTING)
#define GET_CHECKPOINT(toplist,total,count,outputname,cptname)\
  { int ret = init_and_read_checkpoint(toplist,total,count,outputname,cptname);\
    if(ret < 0) {\
      fprintf(stderr, HIERARCHICALSEARCH_MSGCHECKPT "\n");\
      return(HIERARCHICALSEARCH_ECHECKPT);\
    } else if (ret == 2) {\
      return(HIERARCHICALSEARCH_ENORM);\
    }\
  }
#define SET_CHECKPOINT set_checkpoint()
#define INSERT_INTO_FSTAT_TOPLIST add_checkpoint_candidate

#else
#define SET_CHECKPOINT
#define GET_CHECKPOINT(toplist,total,count,outputname,cptname) *total=0;
#define INSERT_INTO_FSTAT_TOPLIST insert_into_fstat_toplist
#endif

#ifdef  __cplusplus
extern "C" {
#endif

extern LALStatus *global_status;

/* function prototypes, they are defined in boinc_extras.c */

/** allows the App to register another output file to be put into the
    zip archive that is sent back to the server */
extern void register_output_file(char*filename);

/** show progress of the App.
    NOTE: This also set the count & total (skypos) for checkpointing */
extern void show_progress(double rac, double dec, UINT4 count, UINT4 total);

/** inits checkpointing for the toplist and reads the last checkpoint if present
    This expects all passed variables (toplist, total, count) to be already initialized.
    The variables are modified only if a previous checkpoint was found.
    If *cptname (name of the checkpoint file) is NULL,
    the name is constructed by appending ".cpt" to the output filename.
    The FILE* should be the one that checpointed_fopen() above has returned.

    The function returns
    0 if no checkpoint could be found,
   -1 if a checkpoint was found but it or the previous output couldn't be read,
   -2 if an error occured (out of memory),
    1 if a checkpoint was found and previous output could be read
    2 a previously written end marker was detected
*/
extern int init_and_read_checkpoint(toplist_t*toplist, UINT4*count,
				     UINT4 total, char*outputname, char*cptname);

/** This corresponds to insert_into_fstat_toplist().
    It inserts a candidate into the toplist, updates the file
    and "compacts" it if necessary (i.e. bytes > maxsize).
    NOTE that the toplist parameter is just a dummy to make the interface
         compatible to insert_into_fstat_toplist(). The operation is
         actually performed on the toplist passed to the least recent call
         of init_and_read_checkpoint(), which, however, should be the same
         in all reasonable cases. */
extern int add_candidate_and_checkpoint (toplist_t*toplist, FstatOutputEntry cand);

/** replacement for add_candidate_and_checkpoint() */
extern int add_checkpoint_candidate (toplist_t*toplist, FstatOutputEntry cand);

/** actually writes a checkpoint only if it's "boinc time to checkpoint"
    and compacts the output file if necessary */
extern void set_checkpoint(void);

/** does the final (compact) write of the file and cleans up checkpointing stuff
    The checkpoint file remains there in case the App gets inteerupted afterwards
    but befor boinc_finish was called */
extern void write_and_close_checkpointed_file (void);

/** LALApps error handler for BOINC */
int BOINC_LAL_ErrHand (LALStatus*, const char*, const char*, const int, volatile const char*);


/** the main() function of HierarchicalSerach.c becomes the extern MAIN(),
    the real main() function of the BOINC App is defined in boinc_extras.c
*/
extern int MAIN(int,char**);

#ifdef  __cplusplus
}
#endif

