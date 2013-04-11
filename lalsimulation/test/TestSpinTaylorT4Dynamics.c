/*
 *  Copyright (C) 2011 Evan Ochsner
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
#include <time.h>

#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/LALSimInspiral.h>
#include <lal/XLALError.h>
#include <lal/LALAdaptiveRungeKutta4.h>

#define UNUSED(expr) do { (void)(expr); } while (0)

int main (int argc , char **argv) 
{
    FILE *f;
    int i, len;
    REAL8 lnhatx, lnhaty, lnhatz, e1x, e1y, e1z;
    REAL8TimeSeries *V = NULL;
    REAL8TimeSeries *Phi = NULL;
    REAL8TimeSeries *S1x = NULL;
    REAL8TimeSeries *S1y = NULL;
    REAL8TimeSeries *S1z = NULL;
    REAL8TimeSeries *S2x = NULL;
    REAL8TimeSeries *S2y = NULL;
    REAL8TimeSeries *S2z = NULL;
    REAL8TimeSeries *LNhatx = NULL;
    REAL8TimeSeries *LNhaty = NULL;
    REAL8TimeSeries *LNhatz = NULL;
    REAL8TimeSeries *E1x = NULL;
    REAL8TimeSeries *E1y = NULL;
    REAL8TimeSeries *E1z = NULL;
    REAL8 m1 = 12.0 * LAL_MSUN_SI;
    REAL8 m2 = 3.0 * LAL_MSUN_SI;
    REAL8 s1x = 0.7;
    REAL8 s1y = 0.2;
    REAL8 s1z = 0.1;
    REAL8 s2x = 0.1;
    REAL8 s2y = 0.8;
    REAL8 s2z = 0.1;
    REAL8 inclination = 1.;
    REAL8 lambda1 = 500.;
    REAL8 lambda2 = 500.;
    REAL8 quadparam1 = 1.;
    REAL8 quadparam2 = 1.;
    LALSimInspiralSpinOrder spinO = -1; // All spin effects
    LALSimInspiralTidalOrder tideO = -1; // All tidal effects
    REAL8 deltaT = 1. / 16384.;
    REAL8 fStart = 40.;
    REAL8 fEnd = 0.;
    INT4 phaseO = 7;
    Approximant approx = SpinTaylorT4;

    lnhatx = sin(inclination);
    lnhaty = 0.;
    lnhatz = cos(inclination);
    e1x = cos(inclination);
    e1y = 0.;
    e1z = -sin(inclination);

    UNUSED(argc);
    UNUSED(argv);

    // Test SpinTaylorT4
    XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi, &S1x, &S1y,
            &S1z, &S2x, &S2y, &S2z, &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z,
            deltaT, m1, m2, fStart, fEnd, s1x, s1y, s1z, s2x, s2y,
            s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
            quadparam1, quadparam2, spinO, tideO, phaseO, approx);
    len = V->data->length;
    f = fopen("ST4-dynamics.dat", "w");
    for(i = 0; i < len; i++)
    {
        fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                i * deltaT, Phi->data->data[i], pow( V->data->data[i], 3.),
                LNhatx->data->data[i], LNhaty->data->data[i],
                LNhatz->data->data[i], S1x->data->data[i], S1y->data->data[i],
                S1z->data->data[i], S2x->data->data[i], S2y->data->data[i],
                S2z->data->data[i], E1x->data->data[i],
                E1y->data->data[i], E1z->data->data[i]);
    }
    fclose(f);

    // Test SpinTaylorT2
    approx = SpinTaylorT2;
    XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi, &S1x, &S1y,
            &S1z, &S2x, &S2y, &S2z, &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z,
            deltaT, m1, m2, fStart, fEnd, s1x, s1y, s1z, s2x, s2y,
            s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
            quadparam1, quadparam2, spinO, tideO, phaseO, approx);
    len = V->data->length;
    f = fopen("ST2-dynamics.dat", "w");
    for(i = 0; i < len; i++)
    {
        fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                i * deltaT, Phi->data->data[i], pow( V->data->data[i], 3.),
                LNhatx->data->data[i], LNhaty->data->data[i],
                LNhatz->data->data[i], S1x->data->data[i], S1y->data->data[i],
                S1z->data->data[i], S2x->data->data[i], S2y->data->data[i],
                S2z->data->data[i], E1x->data->data[i],
                E1y->data->data[i], E1z->data->data[i]);
    }
    fclose(f);

    /* Destroy vectors of dynamical variables, check for errors then exit */
    XLALDestroyREAL8TimeSeries(V);
    XLALDestroyREAL8TimeSeries(Phi);
    XLALDestroyREAL8TimeSeries(S1x);
    XLALDestroyREAL8TimeSeries(S1y);
    XLALDestroyREAL8TimeSeries(S1z);
    XLALDestroyREAL8TimeSeries(S2x);
    XLALDestroyREAL8TimeSeries(S2y);
    XLALDestroyREAL8TimeSeries(S2z);
    XLALDestroyREAL8TimeSeries(LNhatx);
    XLALDestroyREAL8TimeSeries(LNhaty);
    XLALDestroyREAL8TimeSeries(LNhatz);
    XLALDestroyREAL8TimeSeries(E1x);
    XLALDestroyREAL8TimeSeries(E1y);
    XLALDestroyREAL8TimeSeries(E1z);

    return 0;
}
