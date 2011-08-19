/*
*  Copyright (C) 2011 Nickolas Fotopoulos
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
#include <stdlib.h>

#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/XLALError.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*
 * main
 */
int main (int UNUSED argc , char UNUSED **argv) {
    FILE *f;
    ssize_t i;
    REAL8 deltaF = 0.25;
    LIGOTimeGPS tRef = {0., 0.};
    COMPLEX16FrequencySeries *htilde = NULL;
    COMPLEX16 *dataPtr;

    /* fail hard */
    XLALSetErrorHandler(XLALAbortErrorHandler);

    /* generate waveform */
    XLALSimIMRPhenomAGenerateFD(&htilde, &tRef, 0., 10., deltaF, 5. * LAL_MSUN_SI, 5.1 * LAL_MSUN_SI, 10., 2000., 1e8 * LAL_PC_SI);

    /* dump file */
    f = fopen("test.dat", "w");
    fprintf(f, "# f htilde.re htilde.im\n");
    dataPtr = htilde->data->data;
    for (i=0; i < htilde->data->length; i++)
      fprintf(f, "%e %e %e\n", i * deltaF, dataPtr[i].re, dataPtr[i].im);
    fclose(f);
    XLALDestroyCOMPLEX16FrequencySeries(htilde);
    return 0;
}
