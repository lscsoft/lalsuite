/*
 *  Copyright (C) 2007 Badri Krishnan
 *  Copyright (C) 2008 Christine Chung, Badri Krishnan and John Whelan
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
 * \author Christine Chung, Badri Krishnan, John Whelan
 * \date 2008
 * \file 
 * \ingroup pulsar
 * \brief Header-file for LAL routines for CW cross-correlation searches
 *
 * $Id$
 *
 */
 
/*
 *   Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _PULSARCROSSCORR_H
#define _PULSARCROSSCORR_H

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if HAVE_GLOB_H
#include <glob.h>
#endif
#include <time.h>
#include <errno.h> 

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/LALDatatypes.h>
#include <lal/LALHough.h>
#include <lal/RngMedBias.h>
#include <lal/LALRunningMedian.h>
#include <lal/Velocity.h>
#include <lal/Statistics.h>
#include <lal/ComputeFstat.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SFTClean.h>
#include <gsl/gsl_cdf.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>


/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif


/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID (PULSARCROSSCORRH, "$Id$");

/******************************************************
 *  Error codes and messages.
 */
 
#define PULSARCROSSCORR_ENULL 1
#define PULSARCROSSCORR_ENONULL 2
#define PULSARCROSSCORR_EVAL 3

#define PULSARCROSSCORR_MSGENULL "Null pointer"
#define PULSARCROSSCORR_MSGENONULL "Non-null pointer"
#define PULSARCROSSCORR_MSGEVAL "Invalid value"

/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */

  typedef enum
  { SAME,
    DIFFERENT,
    ALL
  } DetChoice;


  /** struct holding info about skypoints */
  typedef struct tagSkyPatchesInfo{
    UINT4 numSkyPatches;
    REAL8 *alpha;
    REAL8 *delta;
    REAL8 *alphaSize;
    REAL8 *deltaSize;
  } SkyPatchesInfo;

  typedef struct tagSFTDetectorInfo{
    COMPLEX8FrequencySeries *sft;
    REAL8 vDetector[3];
    REAL8 rDetector[3];
    REAL8 a;
    REAL8 b;
  } SFTDetectorInfo;

  /* define structs to hold combinations of F's and A's */
  typedef struct tagCrossCorrAmps {
    REAL8 Aplussq;
    REAL8 Acrosssq;
    REAL8 AplusAcross;
  } CrossCorrAmps;

  typedef struct tagCrossCorrBeamFn{
    REAL8 a;
    REAL8 b;
  } CrossCorrBeamFn;

  typedef struct tagSFTListElement {
    SFTtype sft;
    struct SFTListElement *nextSFT;
  } SFTListElement;

  typedef struct tagPSDListElement {
    REAL8FrequencySeries psd;
    struct PSDListElement *nextPSD;
  } PSDListElement;

  typedef struct tagREAL8ListElement {
    REAL8 val;
    struct REAL8ListElement *nextVal;
  } REAL8ListElement;

  typedef struct {
    CrossCorrBeamFn beamfn;
    struct CrossCorrBeamFnListElement *nextBeamfn;
  } CrossCorrBeamFnListElement;
/*
 *  Functions Declarations (i.e., prototypes).
 */

void LALCreateSFTPairsIndicesFrom2SFTvectors(LALStatus          *status,
					     INT4VectorSequence **out,
					     SFTListElement     *in,
					     REAL8	        lag,
					     INT4		listLength,
					     DetChoice 		detChoice,
					     BOOLEAN	        autoCorrelate);

void LALCorrelateSingleSFTPair(LALStatus                *status,
			       COMPLEX16                *out,
			       COMPLEX8FrequencySeries  *sft1,
			       COMPLEX8FrequencySeries  *sft2,
			       REAL8FrequencySeries     *psd1,
			       REAL8FrequencySeries     *psd2,
			       REAL8                    freq1,
			       REAL8                    freq2);

void LALGetSignalFrequencyInSFT(LALStatus                *status,
				REAL8                    *out,
				LIGOTimeGPS		 *epoch,
				PulsarDopplerParams      *dopp,
				REAL8Vector              *vel);

void LALGetSignalPhaseInSFT(LALStatus               *status,
			    REAL8                   *out,
			    LIGOTimeGPS		    *epoch,
			    PulsarDopplerParams     *dopp,
			    REAL8Vector             *pos);

void LALCalculateSigmaAlphaSq(LALStatus            *status,
			      REAL8                *out,
			      REAL8                freq1,
			      REAL8                freq2,
			      REAL8FrequencySeries *psd1,
			      REAL8FrequencySeries *psd2);

void LALCalculateAveUalpha(LALStatus *status,
			COMPLEX16 *out,
			REAL8     phiI,
			REAL8     phiJ,
			CrossCorrBeamFn beamfnsI,
			CrossCorrBeamFn beamfnsJ,
			REAL8     sigmasq);

void LALCalculateUalpha(LALStatus *status,
			COMPLEX16 *out,
			CrossCorrAmps amplitudes,
			REAL8     phiI,
			REAL8     phiJ,
			CrossCorrBeamFn beamfnsI,
			CrossCorrBeamFn beamfnsJ,
			REAL8     sigmasq,
			REAL8     *psi);

void LALCalculateCrossCorrPower(LALStatus       *status,
				REAL8	        *out,
				COMPLEX16Vector *yalpha,
				COMPLEX16Vector *ualpha);

void LALNormaliseCrossCorrPower(LALStatus        *status,
				REAL8		 *out,
				COMPLEX16Vector  *ualpha,
				REAL8Vector      *sigmaAlphasq);

/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _PULSARCROSSCORR_H */
