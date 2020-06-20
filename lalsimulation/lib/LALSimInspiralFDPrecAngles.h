#ifndef _LALSIM_INS_FD_PREC_ANGLES
#define _LALSIM_INS_FD_PREC_ANGLES

/*
 * Copyright (C) 2017 Katerina Chatziioannou, Sebastian Khan
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

typedef struct tagvector
{
    double x;
    double y;
    double z;
} vector;

typedef struct tagsystemquantites
{
    double onethird;
    double constants_u[6];
    double constants_phiz[6];
    double constants_zeta[6];
    double constants_L[6];
    double phiz_0, zeta_0, constant_of_S;
    double c_1, Ssqave, sqrtSsqave, Seff, c1_2, nu_2, nu_4, c_1_over_nu, S0_norm;
    double S1_norm_2, S2_norm_2;
    double dot1, dot2, dot12, dot1n, dot2n;
    double deltam_over_M, nu, q;
} sysq;

#endif /* _LALSIM_INS_FD_PREC_ANGLES */
