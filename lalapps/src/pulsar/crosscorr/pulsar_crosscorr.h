/*
 *  Copyright (C) 2007 Badri Krishnan
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
 * \brief Header for CW cross-correlation search
 *
 * $Id: pulsar_crosscorr.h,v 1.7 2009/03/13 00:43:04 cchung Exp $
 *
 */
 
/*
 *   Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _PULSAR_CROSSCORR_H
#define _PULSAR_CROSSCORR_H

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
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
#include <lalapps.h>
#include <gsl/gsl_cdf.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/PulsarCrossCorr.h>


/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif


/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID (PULSAR_CROSSCORRH, "$Id: pulsar_crosscorr.h,v 1.7 2009/03/13 00:43:04 cchung Exp $");

/******************************************************
 *  Error codes and messages.
 */
 
#define PULSAR_CROSSCORR_ENORM 0
#define PULSAR_CROSSCORR_ESUB  1
#define PULSAR_CROSSCORR_EARG  2
#define PULSAR_CROSSCORR_EBAD  3
#define PULSAR_CROSSCORR_EFILE 4
#define PULSAR_CROSSCORR_EDIR 5
#define PULSAR_CROSSCORR_ENULL 6
#define PULSAR_CROSSCORR_ENONULL 7
#define PULSAR_CROSSCORR_EVAL 8
#define PULSAR_CROSSCORR_EMEM 9

#define PULSAR_CROSSCORR_MSGENORM "Normal exit"
#define PULSAR_CROSSCORR_MSGESUB  "Subroutine failed"
#define PULSAR_CROSSCORR_MSGEARG  "Error parsing arguments"
#define PULSAR_CROSSCORR_MSGEBAD  "Bad argument values"
#define PULSAR_CROSSCORR_MSGEFILE "Could not create output file"
#define PULSAR_CROSSCORR_MSGEDIR  "Could not create directory"
#define PULSAR_CROSSCORR_MSGENULL "Null pointer"
#define PULSAR_CROSSCORR_MSGENONULL "Non-null pointer"
#define PULSAR_CROSSCORR_MSGEVAL "Invalid value"
#define PULSAR_CROSSCORR_MSGEMEM "Out of memory"

#define PIXELFACTOR  2 


/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */

/*
 *  Functions Declarations (i.e., prototypes).
 */

void SetUpRadiometerSkyPatches(LALStatus *status,
			       SkyPatchesInfo *out,
			       CHAR *skyFileName,
			       CHAR *skyRegion,
			       REAL8 dAlpha,
			       REAL8 dDelta);

void InitDoppParams(LALStatus *status,
 		    REAL8Vector *fdots,
		    PulsarDopplerParams *thisPoint,
		    LIGOTimeGPS refTime,
  		    REAL8 f_current,
 		    REAL8 q1_current,
		    REAL8 q2_current,
	 	    REAL8 n_current,
		    REAL8 fdot_current,
		    REAL8 fddot_current);


void GetBeamInfo(LALStatus *status, 
		 CrossCorrBeamFnListElement *beamHead, 
		 SFTListElement *sftHead, 
		 REAL8ListElement *freqHead,
		 REAL8ListElement *phaseHead, 
		 SkyPosition skypos, 
		 EphemerisData *edat, 
		 PulsarDopplerParams *thisPoint);

void CalculateFdots (LALStatus *status,
		     REAL8Vector *fdots,
		     REAL8 f0,
		     REAL8 q1,
		     REAL8 q2,
		     REAL8 n);

void CopySFTFromCatalog(LALStatus *status,
		   	SFTCatalog *catalog,
			SFTVector **sft,
			REAL8 fMin,
			REAL8 fMax,
			INT4 sftindex);

void AddSFTtoList(LALStatus *status,
		  SFTListElement **sftHead,
		  SFTListElement **sftTail,
		  SFTtype *sft);

void AddPSDtoList(LALStatus *status,
		  PSDListElement **psdHead,
		  PSDListElement **psdTail,
		  INT4 length);

void AddREAL8toList(LALStatus *status,
		    REAL8ListElement **head,
		    REAL8ListElement **tail);

void AddBeamFntoList(LALStatus *status,
		     CrossCorrBeamFnListElement **beamHead,
		     CrossCorrBeamFnListElement **beamTail);

void DeleteSFTHead (LALStatus *status, 
		    SFTListElement **sftHead);

void DeletePSDHead (LALStatus *status, 
		    PSDListElement **psdHead);

void DeleteREAL8Head (LALStatus *status,
		      REAL8ListElement **head);

void DeleteBeamFnHead (LALStatus *status,
		       CrossCorrBeamFnListElement **beamHead);

		     
/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _PULSAR_CROSSCORR_H */
