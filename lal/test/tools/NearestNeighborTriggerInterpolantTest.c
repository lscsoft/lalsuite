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
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 *
 * Copyright (C) 2012 Leo Singer
 */


#include <complex.h>
#include <stdlib.h>

#include <lal/TriggerInterpolation.h>

int main(__attribute__ ((unused)) int argc, __attribute__ ((unused)) char **argv)
{
    int result;
    double tmax;

    NearestNeighborTriggerInterpolant *interp = XLALCreateNearestNeighborTriggerInterpolant(0);
    if (!interp)
        exit(EXIT_FAILURE);

    {
        const COMPLEX16 y[] = {1+2*I};
        COMPLEX16 ymax;

        result = XLALCOMPLEX16ApplyNearestNeighborTriggerInterpolant(interp, &tmax, &ymax, y);
        if (result)
            exit(EXIT_FAILURE);

        if (tmax != 0)
            exit(EXIT_FAILURE);
        if (creal(ymax) != 1)
            exit(EXIT_FAILURE);
        if (cimag(ymax) != 2)
            exit(EXIT_FAILURE);
    }

    {
        const COMPLEX8 y[] = {1+2*I};
        COMPLEX8 ymax;

        result = XLALCOMPLEX8ApplyNearestNeighborTriggerInterpolant(interp, &tmax, &ymax, y);
        if (result)
            exit(EXIT_FAILURE);

        if (tmax != 0)
            exit(EXIT_FAILURE);
        if (crealf(ymax) != 1)
            exit(EXIT_FAILURE);
        if (cimagf(ymax) != 2)
            exit(EXIT_FAILURE);
    }

    {
        const REAL8 y[] = {1};
        REAL8 ymax;

        result = XLALREAL8ApplyNearestNeighborTriggerInterpolant(interp, &tmax, &ymax, y);
        if (result)
            exit(EXIT_FAILURE);

        if (tmax != 0)
            exit(EXIT_FAILURE);
        if (ymax != 1)
            exit(EXIT_FAILURE);
    }

    {
        const REAL4 y[] = {1};
        REAL4 ymax;

        result = XLALREAL4ApplyNearestNeighborTriggerInterpolant(interp, &tmax, &ymax, y);
        if (result)
            exit(EXIT_FAILURE);

        if (tmax != 0)
            exit(EXIT_FAILURE);
        if (ymax != 1)
            exit(EXIT_FAILURE);
    }

    XLALDestroyNearestNeighborTriggerInterpolant(interp);
    exit(EXIT_SUCCESS);
}
