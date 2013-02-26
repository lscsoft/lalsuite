/*
*  Copyright (C) 2013 Alex Nitz
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
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/XLALError.h>
#include <math.h>
#include <stdio.h>

#include <lal/LALInspiralPyCBCTemplate.h>

void XLALInspiralPyCBCTemplatePhase (COMPLEX8Vector* htilde, REAL4Vector* sincos_look,
                    REAL4Vector* f13, REAL4Vector* logv_look,
                    int kmin, int  phase_order,
                    float delta_f, float piM, float pfaN,
                    float pfa2, float pfa3, float pfa4, float pfa5, float pfl5,
                    float pfa6, float pfl6, float pfa7, float tC, float v0){

    float dp = LAL_TWOPI / (sincos_look->length);
    float piM13 = cbrtf(piM);
    float logpiM13 = log(piM13);
    float logv0 = log(v0);
    float log4 = log(4);

    if ( !htilde )
        XLAL_ERROR_VOID(XLAL_EFAULT, "XLALInspiralPyCBCTemplatePhase");

    for (unsigned int i=0; i< htilde->length; i++){
        int index = i + kmin;
        const float f = index * delta_f;
        const float v =  piM13 * f13->data[index];
        const float logv = logv_look->data[index] * 1.0/3.0 + logpiM13;
        const float v5 = v * v * v * v * v;
        float phasing = 0.;
        float shft = -LAL_TWOPI * tC;

        switch (phase_order)
        {
            case -1:
            case 7:
                phasing = pfa7 * v;
            case 6:
                phasing = (phasing + pfa6 + pfl6 * (logv + log4) ) * v;
            case 5:
                phasing = (phasing + pfa5 + pfl5 * (logv - logv0) ) * v;
            case 4:
                phasing = (phasing + pfa4) * v;
            case 3:
                phasing = (phasing + pfa3) * v;
            case 2:
                phasing = (phasing + pfa2) * v * v;
            case 0:
                phasing += 1.;
                break;
            default:
                break;
        }

        phasing *= pfaN / v5;
        phasing += shft * f + LAL_PI_4;

        float sphase = phasing - (int) (phasing / LAL_TWOPI) * LAL_TWOPI;
        float cphase = (phasing + LAL_PI_2) - (int) ((phasing + LAL_PI_2 ) / LAL_TWOPI) * LAL_TWOPI;

        int sindex = sphase / dp;
        int cindex = cphase / dp;

        float pcos = sincos_look->data[cindex];
        float psin = sincos_look->data[sindex];
        htilde->data[i] = (pcos - psin*I);
    }
}
