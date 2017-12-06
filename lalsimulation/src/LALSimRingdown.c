/*
 * Copyright (C) 2011 J. Clark
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

#include <complex.h>

#include <lal/LALSimRingdown.h>

#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>
#include <lal/Date.h>

#define EPS LAL_REAL4_EPS

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


/**
 * Computes the waveform for the ringdown of a Newtonian neutron star
 * quasinormal mode (l,m).
 *
 */
int XLALSimRingdown(
    REAL8TimeSeries **hplus,	/**< plus-polarization waveform [returned] */
    REAL8TimeSeries **hcross,	/**< cross-polarization waveform [returned] */
	double delta_t,			/**< sampling interval (s) */
	double Amp,				/**< initial intrinsic amplitude of ringdown */
	double omega0,			/* f-mode oscillation frequency (rad) */
	double quality,			/* quality factor = pi*decay*frequency*/
	double phi0,			/**< initial phase of ringdown (rad) */
    double theta,		/**< inclination of source's spin axis (rad) */
    double azimuth,		/**< azimuthal angle */
    int l,		        /**< polar mode number */
    int m    	        /**< azimuthal mode number */
)
{


	LIGOTimeGPS epoch;
    int length;

	unsigned j; 		  		/* iterator for waveform timestamps */

	/* Out-of-loop calculations */
	COMPLEX16 comega0;			/* Complex frequency comega0 = sqrt(-1) x omega0 */
	COMPLEX16 expPhase;			/* exponential form of phase term */
	COMPLEX16 comega_tau_dt; 	/* complex exponent with time-dependence for ring-down */

	double numfac=0.0;
	double AmpPlus=0.0;		/* Amplitude factor for plus polarisation    */
	double AmpCross=0.0;	/* Amplitude factor for cross polarisation   */
    //double h0Plus;              /* Peak amplitude for plus (includes phase)  */
    //double h0Cross;             /* Peak amplitude for cross (includes phase) */

	double cosTheta;
	cosTheta=cos(theta);

    /* Get tau from quality factor and frequency */
    /* quality = LAL_PI*f0*tau */
    /*         = omega0/2 * tau */
    /* tau = quality/(omega/2) */
    double tau = 2.0*quality/omega0;

    /* length of the injection time series is 5 * the decay time, rounded to
     * the nearest odd integer */
	length = (int) floor(5.0 * tau / delta_t / 2.0);
	length = 2 * length + 1;

	/* the First sample is t = 0 */
	XLALGPSSetREAL8(&epoch, 0.0);

	/* allocate the time series */
	*hplus = XLALCreateREAL8TimeSeries("ringdown +", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("ringdown x", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	if(!*hplus || !*hcross) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}
	
	/* Expressions for spherical harmonics */
	if ( l==2 && m==2 ) {
		numfac = sqrt(2.0*LAL_PI/15.0);
		AmpPlus = -1.0*numfac*(1.0 + cosTheta*cosTheta) * Amp;
		AmpCross = numfac*2.0*cosTheta * Amp;
	}

	/* Time-savers */
	comega0 = I * omega0;
	comega_tau_dt =  (comega0 + 1.0/tau) * -1.0*delta_t;
	expPhase = cexp(I * -(phi0+azimuth));

    /* Actual peak amplitudes */
    //h0Plus  = AmpPlus*creal(expPhase);
    //h0Cross = AmpCross*cimag(expPhase);

    /* populate */
    for (j = 0; j < (*hplus)->data->length; ++j)
    {
        COMPLEX16 h;
        /* h = Ylm * a_lm * exp[-t(1/tau + i*omega0)] * exp(-i*(phi0+azimuth)) */
        h = cexp(comega_tau_dt * (j-floor(0.5*(*hplus)->data->length))) * expPhase;
        (*hplus)->data->data[j]  = AmpPlus*creal(h);
        (*hcross)->data->data[j] = AmpCross*cimag(h);
    }

	return 0;
}


