/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
 * Copyright (C) 2012 Leo Singer
 */


#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include <lal/TriggerInterpolation.h>

int main(__attribute__((unused)) int argc, __attribute__((unused)) char **argv)
{
    int result;
    double tmax;
    COMPLEX16 ymax;

    CubicSplineTriggerInterpolant *interp = XLALCreateCubicSplineTriggerInterpolant(2);
    if (!interp)
        exit(EXIT_FAILURE);

    {
        /* Concave-up, increasing series: maximum should occur at upper sample boundary */
        const COMPLEX16 y[] = {0, 1, 8, 27, 256};

        result = XLALCOMPLEX16ApplyCubicSplineTriggerInterpolant(interp, &tmax, &ymax, &y[2]);
        if (result)
            exit(EXIT_FAILURE);

        if (tmax != 1)
            exit(EXIT_FAILURE);
        if (creal(ymax) != 27)
            exit(EXIT_FAILURE);
        if (cimag(ymax) != 0)
            exit(EXIT_FAILURE);
    }

    {
        /* Concave-up, decreasing series: maximum should occur at lower sample boundary */
        const COMPLEX16 y[] = {256, 27, 8, 1, 0};

        result = XLALCOMPLEX16ApplyCubicSplineTriggerInterpolant(interp, &tmax, &ymax, &y[2]);
        if (result)
            exit(EXIT_FAILURE);

        if (tmax != -1)
            exit(EXIT_FAILURE);
        if (creal(ymax) != 27)
            exit(EXIT_FAILURE);
        if (cimag(ymax) != 0)
            exit(EXIT_FAILURE);
    }

    {
        /* Test some random sample data. */
        const COMPLEX16 y[] = {0, 1.23412285-1.56122275*I, -2.06157986+0.78298714*I, -0.52811449+0.47859189*I, 0.42450302+2.06877714*I};

        result = XLALCOMPLEX16ApplyCubicSplineTriggerInterpolant(interp, &tmax, &ymax, &y[2]);
        if (result)
            exit(EXIT_FAILURE);

        if (fabs(0.108106552865 - tmax) > 1e-8)
            exit(EXIT_FAILURE);
        if (fabs(-2.10041938439 - creal(ymax)) > 1e-8)
            exit(EXIT_FAILURE);
        if (fabs(0.85409051094 - cimag(ymax)) > 1e-8)
            exit(EXIT_FAILURE);
    }

    XLALDestroyCubicSplineTriggerInterpolant(interp);
    exit(EXIT_SUCCESS);
}
