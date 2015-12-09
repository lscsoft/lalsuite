/*
*  Copyright (C) 2010, 2012 Evan Goetz
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

#ifndef __ANTENNA_H__
#define __ANTENNA_H__

#include <lal/VectorMath.h>
#include <lal/DetectorStates.h>

INT4 CompBinShifts(INT4Vector *output, const SSBtimes *ssbTimes, const REAL8 freq, const REAL8 Tsft, const REAL4 dopplerMultiplier);
INT4 CompAntennaPatternWeights(REAL4VectorAligned *output, const SkyPosition skypos, const REAL8 t0, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const BOOLEAN linPolOn, const REAL8 polAngle, const LALDetector det);
INT4 CompAntennaPatternWeights2(REAL4VectorAligned *output, const SkyPosition skypos, const LIGOTimeGPSVector *timestamps, const LALDetector det, const REAL8 *cosi, const REAL8 *psi);
INT4 CompAntennaVelocity(REAL4VectorAligned *output, const SkyPosition skypos, const REAL8 t0, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const LALDetector det, EphemerisData *edat);

REAL4 CompDetectorDeltaVmax(const REAL8 t0, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const LALDetector det, EphemerisData *edat);
REAL4 CompDetectorVmax(const REAL8 t0, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const LALDetector det, EphemerisData *edat);
REAL4 CompDetectorVmax2(const LIGOTimeGPSVector *timestamps, const LALDetector det, EphemerisData *edat);

#endif

