/*
 * Copyright (C) 2015 Jolien Creighton
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

#include <lal/LALStdlib.h>
#include <lal/LALSimSphHarmMode.h>
#include <lal/SphericalHarmonics.h>
#include <lal/TimeSeries.h>
#include "check_series_macros.h"


/**
 * @addtogroup LALSimSphHarmMode_h
 * @{
 */

/**
 * Multiplies a mode h(l,m) by a spin-2 weighted spherical harmonic
 * to obtain hplus - i hcross, which is added to the time series.
 *
 * Implements the sum of a single term of Eq. (11) of:
 * Lawrence E. Kidder, \"Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit\", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 *
 * If sym is non-zero, symmetrically add the m and -m terms assuming
 * that \f$h(l,-m) = (-1)^l h(l,m)*\f$; see Eq. (78) ibid.
 */
int XLALSimAddMode(
		REAL8TimeSeries *hplus,      /**< +-polarization waveform */
	       	REAL8TimeSeries *hcross,     /**< x-polarization waveform */
	       	COMPLEX16TimeSeries *hmode,  /**< complex mode h(l,m) */
	       	REAL8 theta,                 /**< polar angle (rad) */
	       	REAL8 phi,                   /**< azimuthal angle (rad) */
	       	int l,                       /**< mode number l */
	       	int m,                       /**< mode number m */
	       	int sym                      /**< flag to add -m mode too */
		)
{
	COMPLEX16 Y;
	UINT4 j;

	LAL_CHECK_VALID_SERIES(hmode, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hplus, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hcross, XLAL_FAILURE);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hmode, XLAL_FAILURE);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hcross, hmode, XLAL_FAILURE);

	Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
	for ( j = 0; j < hmode->data->length; ++j ) {
		COMPLEX16 hpc;
		hpc = Y * hmode->data->data[j];
		hplus->data->data[j] += creal(hpc);
		hcross->data->data[j] += -cimag(hpc);
	}
	if ( sym ) { /* equatorial symmetry: add in -m mode */
		Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m);
		if ( l % 2 ) /* l is odd */
			Y = -Y;
		for ( j = 0; j < hmode->data->length; ++j ) {
			COMPLEX16 hpc;
			hpc = Y * conj(hmode->data->data[j]);
			hplus->data->data[j] += creal(hpc);
			hcross->data->data[j] += -cimag(hpc);
		}
	}
	return 0;
}

/**
 * For all valid TimeSeries contained within hmode structure,
 * multiplies a mode h(l,m) by a spin-2 weighted spherical harmonic
 * to obtain hplus - i hcross, which is added to the time series.
 *
 * Implements the sum of Eq. (11) of:
 * Lawrence E. Kidder, \"Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit\", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
int XLALSimAddModeFromModes(
		    REAL8TimeSeries *hplus,      /**< +-polarization waveform */
	       	REAL8TimeSeries *hcross,     /**< x-polarization waveform */
	       	SphHarmTimeSeries *hmode,    /**< complex modes h(l,m) */
	       	REAL8 theta,                 /**< polar angle (rad) */
	       	REAL8 phi                    /**< azimuthal angle (rad) */
		)
{
    SphHarmTimeSeries* this = hmode;

    while ( this ) {
        if ( !this->tdata ) {
            this = this->next;
            continue;
        }

        XLALSimAddMode(hplus, hcross, hmode->mode, theta, phi, this->l, this->m, 1);
        this = this->next;
    }
	return 0;
}

/**
 * Returns the h+, hx waveforms constructed from all valid TimeSeries 
 * contained within hmode structure. 
 *
 * @sa XLALSimAddModeFromModes() and XLALSimAddMode()
 */
int XLALSimNewTimeSeriesFromModes(
		REAL8TimeSeries **hplus,     /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross,    /**< x-polarization waveform */
	       	SphHarmTimeSeries *hmode,    /**< complex modes h(l,m) */
	       	REAL8 theta,                 /**< polar angle (rad) */
	       	REAL8 phi                    /**< azimuthal angle (rad) */
		)
{

    if (!hmode) {
        XLALPrintError("NULL mode structure passed.\n");
        XLAL_ERROR(XLAL_EINVAL);
    }
    if (*hplus || *hcross) {
        XLALPrintError("hplus and hcross time series must be NULL.\n");
        XLAL_ERROR(XLAL_EINVAL);
    }

    *hplus = XLALCreateREAL8TimeSeries("hplus", &(hmode->mode->epoch),
                hmode->mode->f0, hmode->mode->deltaT, &lalStrainUnit,
                hmode->mode->data->length);
    *hcross = XLALCreateREAL8TimeSeries("hplus", &(hmode->mode->epoch),
                hmode->mode->f0, hmode->mode->deltaT, &lalStrainUnit,
                hmode->mode->data->length);
    memset((*hplus)->data->data, 0, (*hplus)->data->length*sizeof(REAL8));
    memset((*hcross)->data->data, 0, (*hcross)->data->length*sizeof(REAL8));

    XLALSimAddModeFromModes(*hplus, *hcross, hmode, theta, phi);

	return 0;
}

/** @} */
