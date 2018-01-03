/*
*  Copyright (C) 2007 Anand Sengupta
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

#include <lal/LALStdlib.h>
#include <lal/LALNoiseModelsInspiral.h>
#include <lal/LALConstants.h>

void LALEstimateEffectiveDistance (
        LALStatus               *status,
        InspiralTemplate        param,
        REAL8                   df,
        REAL8Vector             *psd,
        REAL8                   lal_nm_snr,
        REAL8                   *effDistance
        )
{
    REAL8   powerNorm, totalMass, flso, distanceNorm;
    REAL8   f, ins_amp;
    INT4    i;
    REAL8   dynRange = 1.0;

    INITSTATUS(status);
    ATTATCHSTATUSPTR (status);

    powerNorm = 0.;
    totalMass = param.totalMass*LAL_MTSUN_SI;

    flso      = 1.L/(pow(6.L,1.5L)*totalMass*LAL_PI);

    /* Integrate over frequency - we can't start from i=0 (f=0) because f^-7/3 will become infinity */
    for (i=1; i<(INT4)psd->length; i++)  {
        f = i*df;
        if (f > flso) break;
        if (psd->data[i]) powerNorm += pow(f,(-7./3.))/psd->data[i];
    }

    /* I am not sure if there should be a factor of 2.0 here inside the sqrt ()
       i.e distanceNorm = sqrt(2.0*powerNorm * df);
       Multiplying by 2.0 makes dist agree with HW injections */
    distanceNorm = 2.*sqrt(2.*powerNorm * df);

    ins_amp = (LAL_MTSUN_SI * LAL_C_SI / (1.0e6 *  LAL_PC_SI))
            * sqrt( 5.0*param.mu / 96.0 )
            * ( pow( param.totalMass/(LAL_PI*LAL_PI) , 0.33333 ) / pow(LAL_MTSUN_SI, 1.0 / 6.0) ) ;

    distanceNorm *= (ins_amp * sqrt(dynRange));

    /* We need to calculate randIn.SignalAmp = distanceNorm / deff (in
     * Mpc)*/
    (*effDistance) = (distanceNorm / lal_nm_snr);

    /* Normal exit */
    DETATCHSTATUSPTR (status);
    RETURN (status);

}
