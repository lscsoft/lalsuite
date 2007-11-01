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

#define INSERT_INTO_FSTAT_TOPLIST insert_into_fstat_toplist

#if (HS_CHECKPOINTING)
#define GET_CHECKPOINT(toplist,total,count,outputname,cptname)\
  { int ret = init_and_read_checkpoint(toplist,total,count,outputname,cptname);\
    if(ret < 0) {\
      LogPrintf(LOG_CRITICAL, HIERARCHICALSEARCH_MSGCHECKPT " (%d)\n",ret);\
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

extern LALStatus *global_status;

/* function prototypes, they are defined in boinc_extras.c */

/** allows the App to register another output file to be put into the
    zip archive that is sent back to the server */
extern void register_output_file(char*filename);

/** show progress of the App.
    NOTE: This also set the count & total (skypos) for checkpointing */
extern void show_progress(double rac, double dec, UINT4 count, UINT4 total);

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

/** actually writes a checkpoint only if it's "boinc time to checkpoint"
    and compacts the output file if necessary */
extern void set_checkpoint(void);

/** does the final (compact) write of the file and cleans up checkpointing stuff
    The checkpoint file remains there in case the App gets inteerupted afterwards
    but befor boinc_finish was called */
extern void write_and_close_checkpointed_file (void);

/** LALApps error handler for BOINC */
int BOINC_LAL_ErrHand (LALStatus*, const char*, const char*, const int, volatile const char*);

/** attach gdb to the running process; for debugging. */
void attach_gdb(void);

/** play with FPU status and control word
    (BSD doesn't seem to have C99 fenv.h etc.) */
typedef UINT2 fpuw_t;
typedef UINT4 ssew_t;
extern void  set_fpu_control_word(const fpuw_t word);
extern fpuw_t get_fpu_control_word(void);
extern fpuw_t get_fpu_status(void);
extern ssew_t get_sse_control_status(void);
extern void set_sse_control_status(const ssew_t cword);

/* constants in FPU status word and control word mask */
#define FPU_STATUS_INVALID      (1<<0)
#define FPU_STATUS_DENORMALIZED (1<<1)
#define FPU_STATUS_ZERO_DIVIDE  (1<<2)
#define FPU_STATUS_OVERFLOW     (1<<3)
#define FPU_STATUS_UNDERFLOW    (1<<4)
#define FPU_STATUS_PRECISION    (1<<5)
#define FPU_STATUS_STACK_FAULT  (1<<6)
#define FPU_STATUS_ERROR_SUMM   (1<<7)
#define FPU_STATUS_COND_0       (1<<8)
#define FPU_STATUS_COND_1       (1<<9)
#define FPU_STATUS_COND_2       (1<<10)
#define FPU_STATUS_COND_3       (1<<14)
/* for SSE, status and control information is in the same register
   the status bits 0-5 are identical to the FPU status bits,
   the exception mask bits follow */
#define SSE_MASK_INVALID        (1<<7)
#define SSE_MASK_DENORMALIZED   (1<<8)
#define SSE_MASK_ZERO_DIVIDE    (1<<9)
#define SSE_MASK_OVERFLOW       (1<<10)
#define SSE_MASK_UNDERFLOW      (1<<11)
#define SSE_MASK_PRECISION      (1<<12)


/** the main() function of HierarchicalSerach.c becomes the extern MAIN(),
    the real main() function of the BOINC App is defined in boinc_extras.c
*/
extern int MAIN(int,char**);

#ifdef  __cplusplus
}
#endif
