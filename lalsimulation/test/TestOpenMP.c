/*
 *  Copyright (C) 2014 Leo Singer
 *
 *  Check that XLALSimInspiralTaylorF2 gives identical answers no matter how
 *  many OpenMP threads are used.
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


#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>

#include <omp.h>
#include <stdio.h>
#include <string.h>


/* Return 1 if two COMPLEX16FrequencySeries differ,
 * or 0 if they are identical. Print a message to stderr
 * describing the first field that differs. */
static int series_differ(
    COMPLEX16FrequencySeries *a,
    COMPLEX16FrequencySeries *b
) {
    int ret = 1;

    if (strcmp(a->name, b->name) != 0)
        fputs("name differs", stderr);
    else if (XLALGPSCmp(&a->epoch, &b->epoch) != 0)
        fputs("epoch differs", stderr);
    else if (a->f0 != b->f0)
        fputs("f0 differs", stderr);
    else if (a->deltaF != b->deltaF)
        fputs("deltaF differs", stderr);
    else if (XLALUnitCompare(&a->sampleUnits, &b->sampleUnits) != 0)
        fputs("sampleUnits differs", stderr);
    else if (a->data->length != b->data->length)
        fputs("length differs", stderr);
    else if (memcmp(a->data->data, b->data->data,
                    a->data->length * sizeof(a->data->data[0])) != 0)
        fputs("data differs", stderr);
    else
        ret = 0;

    return ret;
}


int main (int argc, char **argv) {
    int num_threads;
    COMPLEX16FrequencySeries *base_htilde = NULL;

    /* Ignore unused parameters. */
    (void)argc;
    (void)argv;

    /* Check that using 2-8 threads gives an answer that is identical
     * to using 1 thread. */
    for (num_threads = 1; num_threads < 8; num_threads++)
    {
        COMPLEX16FrequencySeries *htilde = NULL;
        omp_set_num_threads(num_threads);

        XLALSimInspiralTaylorF2(
            &htilde, 0, 1. / 2048,
            1.4 * LAL_MSUN_SI, 1.4 * LAL_MSUN_SI,
            0, 0, 10, 0., 2048, 1e6 * LAL_PC_SI,
            0, 0,
            LAL_SIM_INSPIRAL_SPIN_ORDER_0PN,
            LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN,
            LAL_PNORDER_THREE_POINT_FIVE,
            LAL_PNORDER_THREE_POINT_FIVE);

        if (num_threads == 1)
            base_htilde = htilde;
        else if (series_differ(base_htilde, htilde))
            return 1;
        else
            XLALDestroyCOMPLEX16FrequencySeries(htilde);
    }

    return 0;
}
