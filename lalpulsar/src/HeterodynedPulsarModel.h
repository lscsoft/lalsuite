/*
 * Copyright (C) 2018 Matthew Pitkin
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

#ifndef _HETERODYNEDPULSARMODEL_H
#define _HETERODYNEDPULSARMODEL_H

#include <math.h>

#include <lal/LALError.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/ReadPulsarParFile.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/SFTutils.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <gsl/gsl_sf_gamma.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \ingroup lalpulsar
 * \author Matthew Pitkin
 * \date 2018
 *
 * \brief Pulsar signal-generation routines for heterodyned data.
 *
 */

/*@{*/


typedef struct tagDetResponseTimeLookupTable{
  LIGOTimeGPS t0;         /**< GPS time epoch for the look-up table */
  UINT4 ntimebins;        /**< number of time bins in the look-up table */
  LALDetector *det;       /**< detector */
  REAL8 alpha;            /**< the source right ascension in radians */
  REAL8 delta;            /**< the source declination in radians */
  REAL8 psi;              /**< the polarisation angle in radians */
  REAL8Vector *fplus;     /**< tensor plus polarisation response */
  REAL8Vector *fcross;    /**< tensor cross polarisation response */
  REAL8Vector *fx;        /**< vector "x" polarisation response */
  REAL8Vector *fy;        /**< vector "y" polarisation response */
  REAL8Vector *fb;        /**< scalar breathing mode polarisation response */
  REAL8Vector *fl;        /**< scalar longitudinal mode polarisation response */
}DetResponseTimeLookupTable;


/* ---------- Function prototypes ---------- */

REAL8Vector *XLALHeterodynedPulsarPhaseDifference( PulsarParameters *params,
                                                   const LIGOTimeGPSVector *datatimes,
                                                   REAL8 freqfactor,
                                                   REAL8Vector *ssbdts,
                                                   UINT4 updateSSBDelay,
                                                   REAL8Vector *bsbdts,
                                                   UINT4 updateBSBDelay,
                                                   const LALDetector *detector,
                                                   const EphemerisData *ephem,
                                                   const TimeCorrectionData *tdat,
                                                   TimeCorrectionType ttype );

REAL8Vector *XLALHeterodynedPulsarGetSSBDelay( PulsarParameters *pars,
                                               const LIGOTimeGPSVector *datatimes,
                                               const LALDetector *detector,
                                               const EphemerisData *ephem,
                                               const TimeCorrectionData *tdat,
                                               TimeCorrectionType ttype );

REAL8Vector *XLALHeterodynedPulsarGetBSBDelay( PulsarParameters *pars,
                                               const LIGOTimeGPSVector *datatimes,
                                               const REAL8Vector *dts,
                                               const EphemerisData *edat );

void XLALGetEarthPosVel( EarthState *earth,
                         const EphemerisData *edat,
                         const LIGOTimeGPS *tGPS );

COMPLEX16TimeSeries* XLALHeterodynedPulsarGetAmplitudeModel( PulsarParameters *pars,
                                                             REAL8 freqfactor,
                                                             UINT4 varyphase,
                                                             UINT4 useroq,
                                                             UINT4 nonGR,
                                                             const LIGOTimeGPSVector *timestamps,
                                                             const DetResponseTimeLookupTable *resp );

COMPLEX16TimeSeries* XLALHeterodynedPulsarGetModel( PulsarParameters *pars,
                                                    REAL8 freqfactor,
                                                    UINT4 varyphase,
                                                    UINT4 useroq,
                                                    UINT4 nonGR,
                                                    const LIGOTimeGPSVector *timestamps,
                                                    REAL8Vector *hetssbdelays,
                                                    UINT4 updateSSBDelay,
                                                    REAL8Vector *hetbsbdelays,
                                                    UINT4 updateBSBDelay,
                                                    const DetResponseTimeLookupTable *resp,
                                                    const EphemerisData *ephem,
                                                    const TimeCorrectionData *tdat,
                                                    TimeCorrectionType ttype );

DetResponseTimeLookupTable* XLALDetResponseLookupTable( REAL8 t0,
                                                        const LALDetector *det,
                                                        REAL8 alpha,
                                                        REAL8 delta,
                                                        UINT4 timeSteps,
                                                        REAL8 avedt );

void XLALDestroyDetResponseTimeLookupTable( DetResponseTimeLookupTable* resp );

void XLALPulsarSourceToWaveformParams( PulsarParameters *params );


/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif