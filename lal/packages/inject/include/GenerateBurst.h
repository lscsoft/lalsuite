/*
*  Copyright (C) 2007 Jolien Creighton, Patrick Brady, Saikat Ray-Majumder
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


#ifndef _GENERATEBURST_H
#define _GENERATEBURST_H


#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <gsl/gsl_rng.h>


/* FIXME:  which of these are still needed? */
#include <lal/LALStdlib.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>
#include <lal/LIGOMetadataTables.h>


#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID(GENERATEBURSTH, "$Id$");


#define GENERATEBURSTH_EMEM 3
#define GENERATEBURSTH_MSGEMEM "Out of memory"


typedef enum {
	sineGaussian,
	Gaussian,
	Ringdown,
	Ringup
} SimBurstType;


typedef struct tagBurstParamStruc {
	REAL8 deltaT;	/* requested sampling interval (s) */
	CoordinateSystem system;	/* coordinate system to assume for simBurst */
} BurstParamStruc;


/*
 * ============================================================================
 *
 *                            Function Prototypes
 *
 * ============================================================================
 */


int XLALBandAndTimeLimitedWhiteNoiseBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 duration,
	REAL8 frequency,
	REAL8 bandwidth,
	REAL8 int_hdot_squared,
	REAL8 delta_t,
	gsl_rng *rng
);


int XLALGenerateStringCusp(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 amplitude,
	REAL8 f_high,
	REAL8 delta_t
);


void LALBurstInjectSignals(
	LALStatus *status,
	REAL4TimeSeries *series,
	SimBurstTable *injections,
	COMPLEX8FrequencySeries *resp,
	INT4 calType
);


#ifdef  __cplusplus
#pragma {
}
#endif


#endif /* _GENERATEBURST_H */
