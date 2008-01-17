/*
 * Copyright (C) 2008 J. Creighton
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

#include <math.h>
#include <lal/LALSimBurst.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Date.h>

/* note: eccentricity and polarization are angles in waveform hrss space */
/* plus is always a cosine gaussian; cross is always a sine gaussian */
int XLALSimBurstSineGaussian( REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LIGOTimeGPS *epoch, REAL8 deltaT, REAL8 Q, REAL8 f0, REAL8 hrss, REAL8 eccentricity, REAL8 polarization )
{
	LIGOTimeGPS start;
	REAL8 duration = 100.0 * Q / f0; /* CHECKME: long enough ??? */
	UINT4 length = floor( 0.5 + duration / deltaT );
	REAL8 a;
	REAL8 b;
	REAL8 cgrss;
	REAL8 sgrss;
	REAL8 hplusrss;
	REAL8 hcrossrss;

	REAL8 h0plus;
	REAL8 h0cross;
	UINT4 j;

	/* semimajor and semiminor axes of waveform ellipsoid */
	a = hrss / sqrt( 2.0 - eccentricity * eccentricity );
	b = a * sqrt( 1.0 - eccentricity * eccentricity );
	
	/* rss of plus and cross polarizations */
	hplusrss  = a * cos( polarization ) - b * sin( polarization );
	hcrossrss = b * cos( polarization ) + a * sin( polarization );

	/* rss of unit amplitude cosine- and sine-gaussian waveforms */
	/* see: K. Riles, LIGO-T040055-00.pdf */
	cgrss = sqrt( (Q/(4.0*f0*sqrt(LAL_PI))) * (1.0 + exp(-Q*Q)) );
	sgrss = sqrt( (Q/(4.0*f0*sqrt(LAL_PI))) * (1.0 - exp(-Q*Q)) );

	/* "peak" amplitudes of plus and cross */
	h0plus  = hplusrss / cgrss;
	h0cross = hplusrss / sgrss;

	/* update length to be even, correct duration, and shift start time */
	length = length % 2 ? length + 1 : length;
	duration = length * deltaT;
	start = *epoch;
	XLALGPSAdd( &start, -0.5 * duration );
	
	*hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &start, 0.0, deltaT, &lalStrainUnit, length );
	*hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &start, 0.0, deltaT, &lalStrainUnit, length );

	for ( j = 0; j < length; ++j ) {
		REAL4 t   = j * deltaT - 0.5 * duration;
		REAL4 phi = LAL_TWOPI * f0 * t;
		REAL4 fac = exp( -0.5*phi*phi/(Q*Q) );
		(*hplus)->data->data[j]  = h0plus * fac * cos( phi );
		(*hcross)->data->data[j] = h0cross * fac * sin( phi );
	}

	return 0;
}
