/*
 *  Copyright (C) 2013 Evan Ochsner
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
 *
 * \author Evan Ochsner
 *
 * \file
 *
 * \brief Check ChooseTD/FDWaveformFromCache is consistent with ChooseWaveoform
 *
 * */

#include <math.h>
#include <time.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>

int main(void) {
    clock_t s1, e1, s2, e2;
    double diff1, diff2;
    unsigned int i;
    REAL8 plusdiff, crossdiff, temp;
    LALSimInspiralWaveformCache *cache = XLALCreateSimInspiralWaveformCache();
    REAL8TimeSeries *hplus = NULL;
    REAL8TimeSeries *hcross = NULL;
    REAL8TimeSeries *hplusC = NULL;
    REAL8TimeSeries *hcrossC = NULL;
    REAL8 m1 = 10. * LAL_MSUN_SI, m2 = 10 * LAL_MSUN_SI;
    REAL8 s1x = 0., s1y = 0., s1z = 0., s2x = 0., s2y = 0., s2z = 0.;
    REAL8 f_min = 40., f_ref = 0., lambda1 = 0., lambda2 = 0.;
    REAL8 phiref = 0., dt = 1./16384.;
    int ret, phaseO = 7, ampO = 0;
    Approximant approx = SEOBNRv1;
    REAL8 inc1 = 0.2, inc2 = 1.3;
    REAL8 dist1 = 1.e6 * LAL_PC_SI, dist2 = 2.e6 * LAL_PC_SI;

    // Generate waveform via usual ChooseWaveform path
    s1 = clock();
    ret = XLALSimInspiralChooseTDWaveform(&hplus, &hcross, phiref, dt, m1, m2,
        s1x, s1y, s1z, s2x, s2y, s2z, f_min, f_ref, dist1, inc1,
        lambda1, lambda2, NULL, NULL, ampO, phaseO, approx);
    e1 = clock();
    diff1 = (double) (e1 - s1) / CLOCKS_PER_SEC;
    if( ret != XLAL_SUCCESS )
        XLAL_ERROR(XLAL_EFUNC);

    // Generate waveform via FromCache - will call ChooseWaveform this 1st time
    s2 = clock();
    ret = XLALSimInspiralChooseTDWaveformFromCache(&hplusC, &hcrossC, phiref,
        dt, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, f_min, f_ref, dist1, inc1,
        lambda1, lambda2, NULL, NULL, ampO, phaseO, approx, cache);
    e2 = clock();
    diff2 = (double) (e2 - s2) / CLOCKS_PER_SEC;
    if( ret != XLAL_SUCCESS )
        XLAL_ERROR(XLAL_EFUNC);

    // Find level of agreement
    plusdiff = crossdiff = 0.;
    for(i=0; i < hplus->data->length; i++)
    {
        temp = abs(hplus->data->data[i] - hplusC->data->data[i]);
        if(temp > plusdiff) plusdiff = temp;
        temp = abs(hcross->data->data[i] - hcrossC->data->data[i]);
        if(temp > crossdiff) crossdiff = temp;
    }
    printf("Comparing waveforms from ChooseTDWaveform and ChooseTDWaveformFromCache\n");
    printf("when both must be generated from scratch...\n");
    printf("ChooseTDWaveform took %f seconds\n", diff1);
    printf("ChooseTDWaveformFromCache took %f seconds\n", diff2);
    printf("Largest difference in plus polarization is: %.16g\n", plusdiff);
    printf("Largest difference in cross polarization is: %.16g\n", crossdiff);

    // Generate another waveform via ChooseWaveform path
    s1 = clock();
    ret = XLALSimInspiralChooseTDWaveform(&hplus, &hcross, phiref, dt, m1, m2,
        s1x, s1y, s1z, s2x, s2y, s2z, f_min, f_ref, dist2, inc2,
        lambda1, lambda2, NULL, NULL, ampO, phaseO, approx);
    e1 = clock();
    diff1 = (double) (e1 - s1) / CLOCKS_PER_SEC;
    if( ret != XLAL_SUCCESS )
        XLAL_ERROR(XLAL_EFUNC);

    // Generate waveform via FromCache - will transform previous waveform
    s2 = clock();
    ret = XLALSimInspiralChooseTDWaveformFromCache(&hplusC, &hcrossC, phiref,
        dt, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, f_min, f_ref, dist2, inc2,
        lambda1, lambda2, NULL, NULL, ampO, phaseO, approx, cache);
    e2 = clock();
    diff2 = (double) (e2 - s2) / CLOCKS_PER_SEC;
    if( ret != XLAL_SUCCESS )
        XLAL_ERROR(XLAL_EFUNC);

    // Find level of agreement
    plusdiff = crossdiff = 0.;
    for(i=0; i < hplus->data->length; i++)
    {
        temp = abs(hplus->data->data[i] - hplusC->data->data[i]);
        if(temp > plusdiff) plusdiff = temp;
        temp = abs(hcross->data->data[i] - hcrossC->data->data[i]);
        if(temp > crossdiff) crossdiff = temp;
    }
    printf("Comparing waveforms from ChooseTDWaveform and ChooseTDWaveformFromCache\n");
    printf("when the latter is cached and transformed...\n");
    printf("ChooseTDWaveform took %f seconds\n", diff1);
    printf("ChooseTDWaveformFromCache took %f seconds\n", diff2);
    printf("Largest difference in plus polarization is: %.16g\n", plusdiff);
    printf("Largest difference in cross polarization is: %.16g\n", crossdiff);

    return 0;
}
