/*
 *  Copyright (C) 2012 Chris Pankow
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
 *
 * \author Chris Pankow
 *
 * \file
 *
 * \brief Testing SphHarmTimeSeries code
 *
 * */


#include "lal/LALSimInspiral.h"
#include "lal/TimeSeries.h"
#include "lal/Date.h"
#include "lal/Units.h"

int main(void){
		lalDebugLevel=7;
		// Empty time series -- technically works, but doesn't make a lot
		// of sense
		SphHarmTimeSeries *ts = XLALSphHarmTimeSeriesAddMode( NULL, NULL, 0, 0 );

		// blow it away
		XLALDestroySphHarmTimeSeries( ts );
		ts = NULL;

		LIGOTimeGPS epoch = {0,0};
		int l;
		int m;

		COMPLEX16TimeSeries *h_lm;
		for( l=0; l<3; l++ ){
			for( m=-l; m<=l; m++ ){
				h_lm = XLALCreateCOMPLEX16TimeSeries(
						"test hlm",
						&(epoch),
						0,
						1.0/16384,
						&(lalStrainUnit),
						100
					);
				ts = XLALSphHarmTimeSeriesAddMode( ts, h_lm, l, m );
				// time series makes a duplicate of the input, so this isn't
				// needed.
				XLALDestroyCOMPLEX16TimeSeries(h_lm);
			}
		}
		printf( "%u major modes added.", XLALSphHarmTimeSeriesGetMaxL( ts )+1 );

		h_lm = XLALCreateCOMPLEX16TimeSeries(
				"test hlm",
				&(epoch),
				0,
				1.0/16384,
				&(lalStrainUnit),
				100
			);

		// Overwrite a component mode
		SphHarmTimeSeries **check = &ts;
		ts = XLALSphHarmTimeSeriesAddMode( ts, h_lm, 2, -2 );
		XLALDestroyCOMPLEX16TimeSeries(h_lm);
		if( check != &ts ){
			fprintf( stderr, "Returned structure does not match after mode overwrite." );
			return 1;
		}

		XLALDestroySphHarmTimeSeries( ts );

		LALCheckMemoryLeaks();

		return 0;
}
