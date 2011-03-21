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


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>


#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID(GENERATEBURSTH, "$Id$");


/*
 * ============================================================================
 *
 *                            Function Prototypes
 *
 * ============================================================================
 */


int XLALGenerateSimBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	const SimBurst *sim_burst,
	double delta_t
);


int XLALBurstInjectSignals(
	REAL8TimeSeries *h,
	const SimBurst *sim_burst,
	const COMPLEX16FrequencySeries *response
);

int XLALBurstInjectHNullSignals(
	REAL8TimeSeries *series,
	const SimBurst *sim_burst
);

#ifdef  __cplusplus
#pragma {
}
#endif


#endif /* _GENERATEBURST_H */
