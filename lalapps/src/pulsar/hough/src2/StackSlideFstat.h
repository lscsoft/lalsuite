/*  
 *  Copyright (C) 2005 Gregory Mendell
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
 * 
 */



/**
 * \file StackSlideFstat.h
 * \brief Header file for StackSlideFstat.c 
 * \author Gregory Mendell
 * \date $Date$
 * 
 ****/


#ifndef _STACKSLIDEFSTAT_H
#define _STACKSLIDEFSTAT_H

/* standard includes */
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h> 

/* lal includes */
#include <lal/UserInput.h>
#include <lal/LALStdlib.h>
#include <lal/PulsarDataTypes.h>
#include <lal/SFTfileIO.h>
#include <lal/AVFactories.h>
#include <lal/RngMedBias.h>
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Velocity.h>
#include <lal/LALDemod.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/Date.h>
#include <lal/LALHough.h> 
#include <lal/NormalizeSFTRngMed.h>
#include <lal/ComputeFstat.h>
#include <lal/Statistics.h>
#include <lal/GeneratePulsarSignal.h> 
#include <lal/LogPrintf.h>
#include "./HierarchicalSearch.h"
/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif


/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID( STACKSLIDEFSTATH, "$Id$" );

/******************************************************
 *  Error codes and messages.
 */
 
#define STACKSLIDEFSTAT_ENORM 0
#define STACKSLIDEFSTAT_ESUB  1
#define STACKSLIDEFSTAT_EARG  2
#define STACKSLIDEFSTAT_EBAD  3
#define STACKSLIDEFSTAT_EFILE 4
#define STACKSLIDEFSTAT_ENULL 5
#define STACKSLIDEFSTAT_EVAL  6
#define STACKSLIDEFSTAT_ENONULL 7

#define STACKSLIDEFSTAT_MSGENORM "Normal exit"
#define STACKSLIDEFSTAT_MSGESUB  "Subroutine failed"
#define STACKSLIDEFSTAT_MSGEARG  "Error parsing arguments"
#define STACKSLIDEFSTAT_MSGEBAD  "Bad argument values"
#define STACKSLIDEFSTAT_MSGEFILE "Could not create output file"
#define STACKSLIDEFSTAT_MSGENULL "Null pointer"
#define STACKSLIDEFSTAT_MSGEVAL "Invalid value"
#define STACKSLIDEFSTAT_MSGENONULL "Pointer not null"


/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */

/* prototypes */

/* Function that stackslides a vector of Fstat frequency series or any REAL8FrequencySeriesVector. */
void StackSlideVecF(LALStatus *status, 
                    SemiCohCandidateList *out, 
                    REAL4FrequencySeriesVector *vecF, 
                    SemiCoherentParams *params);


void StackSlideVecF_HoughMode(LALStatus *status,
			      SemiCohCandidateList  *out,       
			      REAL4FrequencySeriesVector *vecF, 
			      SemiCoherentParams *params);   

/* Calculate f(t) using the master equation given by Eq. 6.18 in gr-qc/0407001 */
/* Returns f(t) in outputPoint.fkdot[0] */
void LALappsFindFreqFromMasterEquation(LALStatus *status, 
                                       PulsarDopplerParams *outputPoint,
                                       PulsarDopplerParams *inputPoint,
                                       REAL8 *vel,
                                       REAL8 deltaT,
                                       UINT2 numSpindown);

/* Get StackSlide candidates using a fixed threshold */
void GetStackSlideCandidates_threshold(LALStatus *status, 
                                       SemiCohCandidateList *out,
                                       REAL8FrequencySeries *stackslideSum, 
                                       PulsarDopplerParams *outputPoint,
                                       PulsarDopplerParams *outputPointUnc,
                                       REAL8 threshold);

/* Get StackSlide candidates as a toplist */
void GetStackSlideCandidates_toplist(LALStatus *status,
                                     toplist_t *list,
                                     REAL8FrequencySeries *stackslideSum,
                                     PulsarDopplerParams *outputPoint,
                                     PulsarDopplerParams *outputPointUnc);

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _STACKSLIDEFSTAT_H */
