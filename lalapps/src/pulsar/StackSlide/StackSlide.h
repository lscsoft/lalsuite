/*
*  Copyright (C) 2007 Gregory Mendell, Virginia Re
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

/* REVISIONS: */
/* 04/26/04 gam; Change LALStackSlide to StackSlide and LALSTACKSLIDE to STACKSLIDE for initial entry to LALapps. */
/* 06/05/04 gam; Add gpsStartTimeSec and gpsStartTimeNan to StackSlideSkyParams; set these to epoch that gives T0 at SSB. */
/* 12/03/04 gam; Clean up indentation; remove extraneous or obsolete comments. */
/* 12/03/04 gam; Add parameter: BOOLEAN divideSUMsByNumSTKs. */
/* 04/12/05 gam; Simplify StackSlideParams struct; change REAL8 **freqDerivData to REAL8 *freqDerivData; */
/* 05/06/05 gam; Add function SumStacks with just creates a SUM from the STKs without sliding */

/**
 * \author M. Landry, G. Mendell
 * \file
 *
 * \brief Computes frequency model, slide stacks accordingly, and sum them.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/StackSlide.h>
 * \endcode
 * This is a short summary...
 */

#ifndef _STACKSLIDE_H
#define _STACKSLIDE_H

#include <lal/LALStdlib.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALBarycenter.h>
 #include <lal/Date.h>
#include <lal/FindRoot.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Author-defined error codes */

/**\name Error Codes */ /*@{*/
#define STACKSLIDEH_ENULL 1
#define STACKSLIDEH_ENNUL 2
#define STACKSLIDEH_ENEGA 4
#define STACKSLIDEH_MSGENULL "Null Pointer"
#define STACKSLIDEH_MSGENNUL "Non-Null Pointer"
#define STACKSLIDEH_MSGENEGA "Bad Negative Value"
#define STACKSLIDECOMPUTESKYH_ENULL 6
#define STACKSLIDECOMPUTESKYH_ENNUL 8
#define STACKSLIDECOMPUTESKYH_ENEGA 10
#define STACKSLIDECOMPUTESKYH_MSGENULL "Null Pointer in StackSlideComputeSky"
#define STACKSLIDECOMPUTESKYH_MSGENNUL "Non-Null Pointer in StackSlideComputeSky"
#define STACKSLIDECOMPUTESKYH_MSGENEGA "Bad Negative Value in StackSlideComputeSky"

#define COMPUTESKYBINARYH_ENULL 1
#define COMPUTESKYBINARYH_ENNUL 2
#define COMPUTESKYBINARYH_ERANG 3
#define COMPUTESKYBINARYH_ENEGA 4
#define COMPUTESKYBINARYH_MSGENULL "Null Pointer"
#define COMPUTESKYBINARYH_MSGENNUL "Non-Null Pointer"
#define COMPUTESKYBINARYH_MSGERANG "Input parameter out of range"
#define COMPUTESKYBINARYH_MSGENEGA "Bad Negative Value"
/*@}*/

#define ACC 1e-9

/**
 * \heading{Structures}
 *
 * \code
 * struct StackSlideParams
 * \endcode
 * \c StackSlideParams
 * This structure contains the parameters for the <tt>StackSlide()</tt> routine.  The parameters are:
 *
 * \code
 * typedef struct
 * tagStackSlideParams
 * {
 *   REAL8 **skyPosData;
 *   REAL8 **freqDerivData;
 *   REAL8 *ParamsSMA;
 *   REAL8 *ParamsTperi;
 *   INT4 numSkyPosTotal;
 *   INT4 numFreqDerivTotal;
 *   BOOLEAN divideSUMsByNumSTKs;
 *   REAL8 f0STK;
 *   REAL8 refFreq;
 *   REAL8 f0SUM;
 *   REAL8 tSTK;
 *   REAL8 tSUM;
 *   INT4  nBinsPerSUM;
 *   INT4  numSTKs;
 *   INT2 binaryFlag;
 *   REAL8 dfSUM;
 *   UINT4 gpsStartTimeSec;
 *   UINT4 gpsStartTimeNan;
 *   INT4 numSpinDown;
 *   EphemerisData *edat;
 *   LIGOTimeGPS *timeStamps;
 *   BarycenterInput *baryinput;
 *   REAL8 		SemiMajorAxis;
 *   REAL8           OrbitalPeriod;
 *   REAL8           OrbitalEccentricity;
 *   REAL8           ArgPeriapse;
 *   UINT4 TperiapseSSBSec;
 *   UINT4 TperiapseSSBNanoSec;
 *   REAL8 deltaSMA;
 *   REAL8 SMAcentral;
 *   INT4 iFreqDeriv;
 *   REAL8 LoudestEvent;
 *   REAL8 peakFreq;
 *   INT4 nMaxSMA;
 * }
 * StackSlideParams;
 * \endcode
 */

/* 04/12/05 gam; Simplify StackSlideParams struct; change REAL8 **freqDerivData to REAL8 *freqDerivData; */
typedef struct
tagStackSlideParams
{
        REAL8 *freqDerivData;
        BOOLEAN divideSUMsByNumSTKs;
        REAL8 f0STK;
        REAL8 f0SUM;
        REAL8 tSTK;
        REAL8 tSUM;
        INT4  nBinsPerSUM;
        INT4  numSTKs;
        REAL8 dfSUM;
        UINT4 gpsStartTimeSec;
        UINT4 gpsStartTimeNan;
        INT4 numSpinDown;
}
StackSlideParams;

typedef struct
tagStackSlideSkyParams
{
	INT8		spinDwnOrder;	/* max spindown parameter order */
	INT8		mObsSFT;	/* number of coherent timescales */
	REAL8		tSFT;		/* timescale of SFT */
	LIGOTimeGPS	*tGPS;		/* GPS time of 1st data sample of each SFT */
	UINT4 gpsStartTimeSec;          /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
	UINT4 gpsStartTimeNan;          /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
	REAL8 		*skyPos; 	/* array of sky positions */
        REAL8 		SemiMajorAxis;  /* orbital radius of binary (in sec) */
        REAL8           OrbitalPeriod;         /* Period of binary (in sec) */
        REAL8           OrbitalEccentricity;   /* Orbital eccentricy */
        REAL8           ArgPeriapse;    /* Argument of Periapse */
        LIGOTimeGPS     TperiapseSSB;   /* Instance of periapse passage measured in the SSB frame */
	UINT4 gpsTperiapseSSBSec;
	BarycenterInput *baryinput;
	EmissionTime *emit;
	EarthState *earth;
	EphemerisData *edat;
	REAL8 dInv;
       }
StackSlideSkyParams;

typedef struct
tagTdotsAndDeltaTs
{
	REAL8  *vecTDots;    /* 1-d array of (dT/dt)'s for frequency calculation */
	REAL8  **vecDeltaTs; /* 2-d array of (T-T_0)'s for frequency calculation */
}
TdotsAndDeltaTs;

void StackSlide(	LALStatus *status,
			REAL4FrequencySeries **SUMData,
			REAL4FrequencySeries **STKData,
			TdotsAndDeltaTs *pTdotsAndDeltaTs,
			StackSlideParams *params);

/* 05/06/05 gam; Add function SumStacks with just creates a SUM from the STKs without sliding */
/*********************************************************************************/
/*              START function: SumStacks                                        */
/*********************************************************************************/
void SumStacks( 	LALStatus *status,
			REAL4FrequencySeries **SUMData,
			REAL4FrequencySeries **STKData,
			StackSlideParams *params);

void StackSlideComputeSky (LALStatus 	*status,
			TdotsAndDeltaTs 	*pTdotsAndDeltaTs,
			StackSlideSkyParams 	*params);

void StackSlideComputeSkyBinary (LALStatus 	*status,
			TdotsAndDeltaTs 	*pTdotsAndDeltaTs,
			INT8 		iSkyCoh,
			StackSlideSkyParams 	*params);

#ifdef __cplusplus
}
#endif

#endif /* _STACKSLIDE_H */
