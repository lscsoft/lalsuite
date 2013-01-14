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
#include <stdlib.h>

#include <lal/TriggerInterpolation.h>

int main(__attribute__((unused)) int argc, __attribute__((unused)) char **argv)
{
    int result;
    double tmax;
    COMPLEX16 ymax;
    const COMPLEX16 y[] = {1+2*I};

    NearestNeighborTriggerInterpolant *interp = XLALCreateNearestNeighborTriggerInterpolant(0);
    if (!interp)
        exit(EXIT_FAILURE);

    result = XLALApplyNearestNeighborTriggerInterpolant(interp, &tmax, &ymax, y);
    if (result)
        exit(EXIT_FAILURE);

    if (tmax != 0)
        exit(EXIT_FAILURE);
    if (creal(ymax) != 1)
        exit(EXIT_FAILURE);
    if (cimag(ymax) != 2)
        exit(EXIT_FAILURE);

    XLALDestroyNearestNeighborTriggerInterpolant(interp);
    exit(EXIT_SUCCESS);
}
