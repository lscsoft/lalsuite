/*
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef _LIGOMETADATAINSPIRALUTILS_H
#define _LIGOMETADATAINSPIRALUTILS_H

#ifdef  __cplusplus
extern "C" {
#endif

#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/Random.h>
#include <lal/SkyCoordinates.h>

/**
 * \defgroup LIGOMetadataInspiralUtils_h Header LIGOMetadataInspiralUtils.h
 * \ingroup lalinspiral_UNCLASSIFIED
 * \brief unknown
 */
/** @{ */

/** @} */ /* end:LIGOMetadataInspiralUtils_h */

/* ---------- inspiral specific functions ---------- */

/* sngl inspiral */
void
LALFreeSnglInspiral (
    LALStatus          *status,
    SnglInspiralTable **eventHead
    );

int
XLALFreeSnglInspiral (
    SnglInspiralTable **eventHead
    );

void
LALSortSnglInspiral (
    LALStatus          *status,
    SnglInspiralTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    );

SnglInspiralTable *
XLALSortSnglInspiral (
    SnglInspiralTable  *eventHead,
    int(*comparfunc)    (const void *, const void *)
    );

int
LALCompareSnglInspiralByTime (
    const void *a,
    const void *b
    );

SnglInspiralTable *
XLALTimeCutSingleInspiral(
    SnglInspiralTable          *eventList,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    );

SnglInspiralTable *
XLALIfoCutSingleInspiral(
    SnglInspiralTable         **eventHead,
    char                       *ifo
    );

INT4
XLALCountSnglInspiral(
    SnglInspiralTable *head
    );

SnglInspiralTable *
XLALMassCut(
    SnglInspiralTable          *eventHead,
    const char                 *massCut,
    REAL4                       massRangeLow,
    REAL4                       massRangeHigh,
    REAL4                       mass2RangeLow,
    REAL4                       mass2RangeHigh
    );

/* sim inspiral */

void
LALGalacticInspiralParamsToSimInspiralTable(
    LALStatus                  *status,
    SimInspiralTable           *output,
    GalacticInspiralParamStruc *input,
    RandomParams               *params
    );

void
LALInspiralSiteTimeAndDist(
    LALStatus         *status,
    SimInspiralTable  *output,
    LALDetector       *detector,
    LIGOTimeGPS       *endTime,
    REAL4             *effDist,
    SkyPosition       *skyPos
    );

void
LALPopulateSimInspiralSiteInfo(
    LALStatus                  *status,
    SimInspiralTable           *output
    );

void
XLALSortSimInspiral(
    SimInspiralTable **head,
    int (*comparefunc)(const SimInspiralTable * const *,
      const SimInspiralTable * const *)
    );

int
XLALCompareSimInspiralByGeocentEndTime(
	const SimInspiralTable * const *a,
	const SimInspiralTable * const *b
    );

int
XLALFreeSimInspiral (
    SimInspiralTable **eventHead
    );

long
XLALSimInspiralAssignIDs (
    SimInspiralTable *head,
    long process_id,
    long simulation_id
    );

int
XLALSimInspiralChirpMassCut(
    SimInspiralTable   **eventHead,
    REAL4                minChirpMass,
    REAL4                maxChirpMass
    );

int
XLALSimInspiralCompMassCut(
    SimInspiralTable   **eventHead,
    REAL4                minCompMass,
    REAL4                maxCompMass,
    REAL4                minCompMass2,
    REAL4                maxCompMass2
    );

int
XLALSimInspiralTotalMassCut(
    SimInspiralTable   **eventHead,
    REAL4                minTotalMass,
    REAL4                maxTotalMass
    );

#ifdef  __cplusplus
}
#endif

#endif /* _LIGOMETADATAINSPIRALUTILS_H */
