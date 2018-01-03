/*
 *  Copyright (C) 2012 Duncan M. Macleod
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
#include <stdio.h>
#include <lal/LALDetCharRange.h>
#include <lal/LALConstants.h>

REAL4 
XLALREAL4InspiralRange(
    const REAL4FrequencySeries *spectrum,
    REAL4 snr,
    REAL4 mass_1,
    REAL4 mass_2,
    REAL4 f_min
)
{
    UINT4 i;
    REAL4 f;

    /* calculate chirp mass */
    REAL4 mchirp = pow(mass_1 * mass_2, 0.6) / pow(mass_1 + mass_2, 0.4)\
                   * LAL_MSUN_SI;

    /* calculate prefactor of integrand */
    REAL4 amp = (5 * pow(LAL_C_SI, 1./3.) *\
                     pow(mchirp * LAL_G_SI / pow(LAL_C_SI, 2), 5./3.) *\
                     pow(1.77, 2)) /
                (96 * pow(LAL_PI, 4./3.) * pow(snr, 2));

    /* calculate isco */
    REAL4 f_isco = pow(LAL_C_SI, 3) /
                   (LAL_G_SI * LAL_MSUN_SI * pow(6, 1.5) *
                    LAL_PI * (mass_1 + mass_2));

    /* integrate */
    REAL4 integral = 0;
    REAL4 Mpc = LAL_PC_SI*1e6;
    for (i=0; i < spectrum->data->length; ++i)
    {
        f = spectrum->f0 + (i * spectrum->deltaF);
        if ((f_min <= f) & (f < f_isco))
        {
            integral += (amp * spectrum->deltaF *\
                            pow(f, -7./3.) / spectrum->data->data[i]) /\
                        pow(Mpc, 2);
        }
    }

    REAL4 inspiral_range = pow(integral, 0.5);

    return inspiral_range;
}

REAL8 
XLALREAL8InspiralRange(
    const REAL8FrequencySeries *spectrum,
    REAL8 snr,
    REAL8 mass_1,
    REAL8 mass_2,
    REAL8 f_min
)
{
    UINT4 i;

    /* calculate chirp mass */
    REAL8 mchirp = pow(mass_1 * mass_2, 0.6) / pow(mass_1 + mass_2, 0.4)\
                   * LAL_MSUN_SI;

    /* calculate prefactor of integrand */
    REAL8 amp = (5 * pow(LAL_C_SI, 1./3.) *\
                     pow(mchirp * LAL_G_SI / pow(LAL_C_SI, 2), 5./3.) *\
                     pow(1.77, 2)) /
                (96 * pow(LAL_PI, 4./3.) * pow(snr, 2));

    /* calculate isco */
    REAL8 f_isco = pow(LAL_C_SI, 3) /
                   (LAL_G_SI * LAL_MSUN_SI * pow(6, 1.5) *
                    LAL_PI * (mass_1 + mass_2));

    /* integrate */
    REAL8 f, integral = 0;
    REAL8 Mpc = LAL_PC_SI*1e6;
    for (i=0; i < spectrum->data->length; ++i)
    {
        f = spectrum->f0 + (i * spectrum->deltaF);
        if ((f_min <= f) & (f < f_isco))
        {
            integral += (amp * spectrum->deltaF *\
                            pow(f, -7./3.) / spectrum->data->data[i]) /\
                        pow(Mpc, 2);
        }
    }

    REAL8 inspiral_range = pow(integral, 0.5);

    return inspiral_range;
}
