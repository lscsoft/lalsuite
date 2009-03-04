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
 * $Id$
 *
 */
 
/*
 *   Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _PULSAR_CROSSCORR_V0_H
#define _PULSAR_CROSSCORR_V0_H

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
#include <lal/PulsarCrossCorr_v0.h>


/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif


/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID (PULSAR_CROSSCORRH, "$Id$");

/******************************************************
 *  Error codes and messages.
 */
 
#define PULSAR_CROSSCORR_V0_ENORM 0
#define PULSAR_CROSSCORR_V0_ESUB  1
#define PULSAR_CROSSCORR_V0_EARG  2
#define PULSAR_CROSSCORR_V0_EBAD  3
#define PULSAR_CROSSCORR_V0_EFILE 4
#define PULSAR_CROSSCORR_V0_EDIR 4
#define PULSAR_CROSSCORR_V0_ENULL 5
#define PULSAR_CROSSCORR_V0_ENONULL 6
#define PULSAR_CROSSCORR_V0_EVAL 5
#define PULSAR_CROSSCORR_V0_EMEM 14

#define PULSAR_CROSSCORR_V0_MSGENORM "Normal exit"
#define PULSAR_CROSSCORR_V0_MSGESUB  "Subroutine failed"
#define PULSAR_CROSSCORR_V0_MSGEARG  "Error parsing arguments"
#define PULSAR_CROSSCORR_V0_MSGEBAD  "Bad argument values"
#define PULSAR_CROSSCORR_V0_MSGEFILE "Could not create output file"
#define PULSAR_CROSSCORR_V0_MSGEDIR  "Could not create directory"
#define PULSAR_CROSSCORR_V0_MSGENULL "Null pointer"
#define PULSAR_CROSSCORR_V0_MSGENONULL "Non-null pointer"
#define PULSAR_CROSSCORR_V0_MSGEVAL "Invalid value"
#define PULSAR_CROSSCORR_V0_MSGEMEM "Out of memory"

#define PIXELFACTOR  2 


/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */


/*
 *  Functions Declarations (i.e., prototypes).
 */

void SetUpRadiometerSkyPatches_v0(LALStatus *status,
			       SkyPatchesInfo_v0 *out,
			       CHAR *skyFileName,
			       CHAR *skyRegion,
			       REAL8 dAlpha,
			       REAL8 dDelta);

/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _PULSAR_CROSSCORR_V0_H */
