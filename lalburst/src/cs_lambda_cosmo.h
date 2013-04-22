/*
*  Copyright (C) 2007 Xavier Siemens
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

/*********************************************************************************/
/*     Cosmological functions for cosmic string burst computation (header)       */
/*                                                                               */
/*                  Jolien Creighton, Irit Maor, Xavier Siemens                  */
/*                                                                               */
/*                         UWM/Caltech - September 2006                          */
/*********************************************************************************/

#ifndef CS_LAMBDA_COSMO_H
#define CS_LAMBDA_COSMO_H

#include <lal/cs_cosmo.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#define LAMBDA_Z_EQ 5440.0		/* redshift of matter-radiation equilib: 10^3.94 */
#define LAMBDA_H_0 2.27e-18		/* Hubble constant (s^-1) */
#define LAMBDA_OMEGA_M 0.279
#define LAMBDA_OMEGA_R 8.5e-5

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* CS_LAMBDA_COSMO_H */
