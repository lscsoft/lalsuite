/*
 * Copyright (C) 2006 S.Fairhurst, B. Krishnan, L.Santamaria
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

#ifndef _NRWAVEINJECT_H
#define _NRWAVEINJECT_H

/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/NRWaveIO.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/SphericalHarmonics.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

/** \addtogroup NRWaveInject_h
 * \author S.Fairhurst, B. Krishnan, L.Santamaria
 *
 * \brief Module for generating h(t) from Numrel waveforms
 *
 */
/*@{*/

/** \name Error codes */
/*@{*/
#define NRWAVEINJECT_ENULL    1	/**< Null pointer */
#define NRWAVEINJECT_EFILE    2	/**< Error in file-IO */
#define NRWAVEINJECT_ENONULL  3	/**< Not a Null pointer */
#define NRWAVEINJECT_ENOMEM   4	/**< Memory ellocation error */
#define NRWAVEINJECT_EVAL     5	/**< Invalid value */
#define NRWAVEINJECT_EFORMAT  6	/**< Meta data file format incorrect */
/*@}*/

/** \cond DONT_DOXYGEN */
#define NRWAVEINJECT_MSGENULL   "Null pointer"
#define NRWAVEINJECT_MSGEFILE   "Error in file-IO"
#define NRWAVEINJECT_MSGENONULL "Not a Null pointer"
#define NRWAVEINJECT_MSGENOMEM  "Memory ellocation error"
#define NRWAVEINJECT_MSGEVAL    "Invalid value"
#define NRWAVEINJECT_MSGEFORMAT "Meta data file format incorrect"
/** \endcond */

#define NINJA_MIN_MODE 2
#define NINJA_MAX_MODE 5

/** enum for list of numrel groups */
typedef enum{
  NINJA_GROUP_AEI = 0,
  NINJA_GROUP_CIT,
  NINJA_GROUP_LSU,
  NINJA_GROUP_JENA,
  NINJA_GROUP_RIT,
  NINJA_GROUP_CORNELL,
  NINJA_GROUP_PSU,
  NINJA_GROUP_FAU,
  NINJA_GROUP_UTB,
  NINJA_GROUP_UIUC,
  NINJA_GROUP_PRINCETON,
  NINJA_GROUP_LAST
} NumRelGroup;


REAL4TimeVectorSeries *
XLALSumStrain(
    REAL4TimeVectorSeries *tempstrain,
    REAL4TimeVectorSeries *strain);

REAL8TimeVectorSeries *
XLALSumStrainREAL8(
    REAL8TimeVectorSeries *tempstrain,
    REAL8TimeVectorSeries *strain);

/* REAL4TimeVectorSeries * */
INT4
XLALOrientNRWave(
    REAL4TimeVectorSeries *strain,
    UINT4                  modeL,
    INT4                   modeM,
    REAL4                  inclination,
    REAL4                  coa_phase);

REAL8TimeVectorSeries *
XLALOrientNRWaveREAL8(
    REAL8TimeVectorSeries *strain,
    UINT4                  modeL,
    INT4                   modeM,
    REAL4                  inclination,
    REAL4                  coa_phase);

void
XLALOrientNRWaveTimeSeriesREAL8(
    REAL8TimeSeries        *plus,
    REAL8TimeSeries        *cross,
    UINT4                  modeL,
    INT4                   modeM,
    REAL4                  inclination,
    REAL4                  coa_phase);

REAL4TimeSeries *
XLALCalculateNRStrain(
    REAL4TimeVectorSeries *strain,
    SimInspiralTable      *thisInj,
    const CHAR            *ifo,
    INT4                   sampleRate);

REAL4TimeSeries *
XLALInterpolateNRWave( REAL4TimeSeries *in,
                       INT4      sampleRate);

REAL8TimeSeries *
XLALInterpolateNRWaveREAL8( REAL8TimeSeries *in,
                       INT4      sampleRate);

INT4
XLALFindNRFile( NRWaveMetaData *out,
                NRWaveCatalog *nrCatalog,
                const SimInspiralTable  *inj,
                INT4  modeL,
                INT4  modeM);

void LALInjectStrainGW( LALStatus *status,
                        REAL4TimeSeries *injData,
                        REAL4TimeVectorSeries *strain,
                        SimInspiralTable *thisInj,
                        CHAR *ifo,
                        REAL8 dynRange);

void LALInjectStrainGWREAL8( LALStatus                 *status,
                             REAL8TimeSeries           *injData,
                             REAL8TimeVectorSeries     *strain,
                             SimInspiralTable          *thisInj,
                             CHAR                      *ifo,
                             REAL8                     dynRange);


INT4
XLALFindNRCoalescenceTime(REAL8 *tc,
                          const REAL4TimeVectorSeries *in);

INT4
XLALFindNRCoalescenceTimeFromhoft(REAL8 *tc,
                                  const REAL4TimeSeries *in);

INT4
XLALFindNRCoalescenceTimeREAL8(REAL8 *tc,
                               const REAL8TimeSeries *in);

INT4
XLALFindNRCoalescencePlusCrossREAL8(REAL8 *tc,
                                    const REAL8TimeSeries *plus,
                                    const REAL8TimeSeries *cross);


/** channel name for nr data in frame file */
CHAR* XLALGetNinjaChannelName(const CHAR *polarisation, /**< either plus or cross */
                              UINT4 l, /**< azimuthal mode index */
                              INT4 m /**< polar mode index */);

NumRelGroup XLALParseNumRelGroupName( CHAR *name);


/*@}*/

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif           /* Close double-include protection _NRWAVEINJECT_H */

