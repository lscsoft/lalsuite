/*
 * Copyright (C) 2008 J. Creighton, T. Creighton
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
#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/SkyCoordinates.h>

REAL8TimeSeries * XLALSimDetectorStrainREAL8TimeSeries( REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, SkyPosition *position, REAL8 psi, LALDetector *detector );

REAL8TimeSeries * XLALSimInjectionREAL8TimeSeries( REAL8TimeSeries *h, LIGOTimeGPS *start, REAL8 deltaT, UINT4 length, COMPLEX16FrequencySeries *response );
