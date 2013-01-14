/*
 * Copyright (C) 2007 S.Fairhurst, B. Krishnan, L.Santamaria, C. Robinson, 
 * C. Pankow
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


/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Factorial.h>

#ifndef _SPHERICALHARMONICS_H
#define _SPHERICALHARMONICS_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

/**
 * \addtogroup SphericalHarmonics_h
 * \author S.Fairhurst, B. Krishnan, L.Santamaria, C. Robinson, C. Pankow
 *
 * \brief Library of Spin-weighted Spherical Harmonic functions and an 
 * implementation of the Wigner-D matrices
 *
 *
 */
/*@{*/
COMPLEX16 XLALSpinWeightedSphericalHarmonic( REAL8 theta, REAL8 phi, int s, int l, int m );
int XLALScalarSphericalHarmonic( COMPLEX16 *y, UINT4 l, INT4  m, REAL8 theta, REAL8 phi );
INT4 XLALSphHarm ( COMPLEX16 *out, UINT4   L, INT4 M, REAL4 theta, REAL4   phi );
double XLALJacobiPolynomial( int n, int alpha, int beta, double x );
double XLALWignerdMatrix( int l, int mp, int m, double beta );
COMPLEX16 XLALWignerDMatrix( int l, int mp, int m, double alpha, double beta, double gam );
/*@}*/


#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif /* _SPHERICALHARMONICS_H */
