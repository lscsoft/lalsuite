/*
 *  Copyright (C) 2012 Chris Pankow, R. O'Shaughnessy
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
 * \author Chris Pankow
 *
 * \file
 *
 * \brief Functions to take an arbitrary waveform time series and impose the
 * effects of causing the viewing angle to precess about a cone of L around J.
 * The cone currently has a constant opening angle.
 *
 */

#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/LALDatatypes.h>
#include <lal/SphericalHarmonics.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

int XLALSimInspiralPrecessionRotateModes(
				SphHarmTimeSeries* h_lm,
				REAL8TimeSeries* alpha,
				REAL8TimeSeries* beta,
				REAL8TimeSeries* gam
);

int XLALSimInspiralConstantPrecessionConeWaveformModes(
				SphHarmTimeSeries** h_lm_tmp,
				double precess_freq,
				double a,
				double phi_precess,
				double alpha_0,
				double beta_0
);

int XLALSimInspiralConstantPrecessionConeWaveform(
				REAL8TimeSeries** hp,
				REAL8TimeSeries** hx,
				SphHarmTimeSeries* h_lm,
				double precess_freq,
				double a,
				double phi_precess,
				double alpha_0,
				double beta_0
);

int XLALSimInspiralOrientationMatrixForL2(
				REAL8 mtx[3][3], 
				COMPLEX16 h22,
				COMPLEX16 h2m2,
				COMPLEX16 h21,
				COMPLEX16 h2m1,
				COMPLEX16 h20
);

int XLALSimInspiralOrientationMatrixDirection(
				REAL8 vec[3],
				REAL8 mtx[3][3]
);
