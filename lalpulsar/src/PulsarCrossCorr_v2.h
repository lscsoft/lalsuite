/*
 *  Copyright (C) 2012, 2013 John Whelan, Shane Larson and Badri Krishnan
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
#ifndef _PULSARCROSSCORRV2_H
#define _PULSARCROSSCORRV2_H

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup PulsarCrossCorr_v2_h Header PulsarCrossCorr_v2.h
 * \ingroup pkg_pulsarCrossCorr
 * \author John Whelan, Shane Larson, Badri Krishnan
 * \date 2012
 * \brief Header-file for XLAL routines for v2 CW cross-correlation searches
 *
 */
/*@{*/

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

/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */

/** Index to refer to an SFT given a set of SFTs from several different detectors */
  typedef struct tagSFTIndex {
    UINT4 detInd; /**< index of detector in list */
    UINT4 sftInd; /**< index of SFT in list for this detector */
  } SFTIndex;

/** List of SFT indices */
  typedef struct tagSFTIndexList {
    UINT4    length; /**< number of SFTs */
    SFTIndex *data; /**< array of SFT indices */
  } SFTIndexList;

/** Index to refer to a pair of SFTs */
  typedef struct tagSFTPairIndex {
#if 0
    SFTIndex sftInd1; /**< index of 1st SFT in pair */
    SFTIndex sftInd2; /**< index of 2nd SFT in pair */
#endif
    UINT4 sftNum[2]; /**< ordinal numbers of first and second SFTs */
  } SFTPairIndex;

/** List of SFT pair indices */
  typedef struct tagSFTPairIndexList {
    UINT4    length; /**< number of SFT Pairs */
    SFTPairIndex *data; /**< array of SFT Pair indices */
  } SFTPairIndexList;

/*
 *  Functions Declarations (i.e., prototypes).
 */

int XLALGetDopplerShiftedFrequencyInfo
  (
   REAL8Vector         *shiftedFreqs,
   UINT4Vector         *lowestBins,
   REAL8Vector         *kappaValues,
   UINT4               numBins,
   PulsarDopplerParams *dopp,
   SFTIndexList        *sfts,
   MultiSSBtimes       *multiTimes,
   REAL8               Tsft
  )
;

int XLALCreateSFTIndexListFromMultiSFTVect
  (
   SFTIndexList        **indexList,
   MultiSFTVector      *sfts
  )
;

int XLALCreateSFTPairIndexList
  (
   SFTPairIndexList  **pairIndexList,
   SFTIndexList       *indexList,
   MultiSFTVector     *sfts,
   REAL8               maxLag,
   BOOLEAN             inclAutoCorr
  )
;

int XLALCalculateCrossCorrSigmaUnshifted
  (
   REAL8Vector      **sigma_alpha,
   SFTPairIndexList  *pairIndexList,
   SFTIndexList      *indexList,
   MultiPSDVector    *psds,
   REAL8              freq,
   REAL8              Tsft       
  )
;

/*@}*/

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _PULSARCROSSCORR_H */
