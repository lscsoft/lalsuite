/*
 *  Copyright (C) 2008, 2009 Jordi Burguet-Castell, Xavier Siemens
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
#include <lal/ComputeDataQualityVector.h>

#include <math.h>  /* to use isnan() and isinf() */


/**
 * Compute the Data Quality Vector as defined in
 * https://www.lsc-group.phys.uwm.edu/daswg/wiki/S6OnlineGroup/CalibratedData
 *
 * Copied for reference:
 *
 * State vector channel: IFO:DMT-STATE_VECTOR
 *
 * The data type is float(!) so they have to be first converted to int.
 * Then, bits 0-4 define the states. Bits 5-15 are always set to 1,
 * such that the word 0xffff means science mode:
 *
 * bit0=SCI (operator set to go to science mode)
 * bit1=CON conlog unsets this bit is non-harmless epics changes
 * bit2=UP (set by locking scripts)
 * bit3=!INJ Injections unset this bit
 * bit4=EXC Unauthorized excitations cause this bit to be unset
 *
 * Channel Name: IFO:DMT-DATA_QUALITY_VECTOR where IFO is one of (H1, H2, L1)
 * Sample Rate: 1 Hz, for any state to be true, it must be true every
 * sample of h(t) within the second
 * Type: INT4
 *
 * The quality channel will use the following bitmask definition
 *
 * SCIENCE       1  // SV_SCIENCE & LIGHT
 * INJECTION     2  // Injection: same as statevector
 * UP            4  // SV_UP & LIGHT
 * CALIBRATED    8  // SV_UP & LIGHT & (not TRANSIENT)
 * BADGAMMA     16  // Calibration is bad (outside 0.8 < gamma < 1.2)
 * LIGHT        32  // Light in the arms ok
 * MISSING      64  // Indication that data was dropped in DMT
 *
 * All the arrays must have been previously allocated.
 * r_x is the number of x_value(s) that are in a single DQ sample (=
 * x.length/n_dq). So it is the "rate" of x with respect to DQ samples.
 *
 * Other special meanings:
 * t_bad_left:  time (in s) of the last NOT-UP event in the Data Quality
 * t_bad_right: time (in s) of the next NOT-UP event in the Data Quality
 * wings:       duration (in s) of the wings used for calibration
 *
 * If t_bad_left < wings then it doesn't matter how much less it
 * is. Same thing for t_bad_right > wings.
 */
int XLALComputeDQ(REAL4* sv_data, int r_sv,
                  REAL4* lax_data, REAL4* lay_data, int r_light,
                  COMPLEX16* gamma_data, int r_gamma,
                  int t_bad_left, int t_bad_right, int wings,
                  int missing,
                  int* dq_data, int n_dq)  /* output */
{
    int i, j;                    /* counters */
    float sum_x, sum_y;
    int light;                   /* is there light in the arms? */
    int sv_science, sv_injection, sv_up;  /* state vector info */
    int science, injection, up;   /* DQ (not identical to the sv_*) */
    int calibrated;
    int badgamma;
    int dq_value;

    /* Fill one by one the contents of the DQ vector */
    for (i = 0; i < n_dq; i++) {
        /* light */
        sum_x = 0;
        sum_y = 0;
        for (j = 0; j < r_light; j++) {
            sum_x += lax_data[i*r_light + j];
            sum_y += lay_data[i*r_light + j];
        }

        light = (sum_x/r_light > 100 && sum_y/r_light > 100);
        /* "is the mean higher than 100 for both arms?" */

        /* science, injection, up (stuff coming from the state vector) */
        sv_science = 1;    /* in science mode */
        sv_injection = 0;  /* with no injection going on */
        sv_up = 1;         /* and IFO is up */
        for (j = 0; j < r_sv; j++) {  /* go over SV samples in a DQ sample */
            int s = (int) sv_data[i*r_sv + j];  /* convert from float to int */
            if ((s & (1 << 0)) == 0)  sv_science = 0;
            if ((s & (1 << 3)) == 0)  sv_injection = 1;
            if ((s & (1 << 2)) == 0)  sv_up = 0;
        }

        science = sv_science && light;
        injection = sv_injection;
        up = sv_up && light;  /* these are the definitions for the DQ vector */


        /* calibrated */
        /* Because we will have to compute UP for the Data Quality
         * everywhere before anything, to know if something funny
         * happens within a "wings" distance from this data, we will
         * leave the computation of the calibrated flag for next
         * loop.
         */
    /*  calibrated = up && (! transient);  */

        /* badgamma */
        badgamma = 0;

        for (j = 0; j < r_gamma; j++) {
            REAL8 re = creal(gamma_data[i*r_gamma + j]);
            if (re < 0.8 || re > 1.2 || isnan(re) || isinf(re))
                badgamma = 1;
        }

        /* data quality */
        dq_value = 0;
        if (science)    dq_value += (1 << 0);
        if (injection)  dq_value += (1 << 1);
        if (up)         dq_value += (1 << 2);
   /*   if (calibrated) dq_value += (1 << 3);  */  /* we'll do that later */
        if (badgamma)   dq_value += (1 << 4);
        if (light)      dq_value += (1 << 5);
        if (missing)    dq_value += (1 << 6);  /* directly from the argument */

        dq_data[i] = dq_value;
    }

    /* Now look for the transients and fill the "calibrated" bit. */
    for (i = 0; i < n_dq; i++) {
        calibrated = 1;

        if (i - wings < t_bad_left || i + wings > n_dq + t_bad_right)
            calibrated = 0;

        for (j = 0; j < wings; j++) {
            int pos = i - j;
            if (pos > 0) {
                if ((dq_data[pos] & (1 << 2)) == 0)  /* if not up */
                    calibrated = 0;
            }
            pos = i + j;  /* note that this includes dq_data[i] having UP=1 */
            if (pos < n_dq) {
                if ((dq_data[pos] & (1 << 2)) == 0)  /* if not up */
                    calibrated = 0;
            }
        }

        if (calibrated) dq_data[i] += (1 << 3);
        /* we checked that we were DQ up all the time in +/- wings
         * seconds, so that also includes sv_up && light */
    }

    return 0;
}
