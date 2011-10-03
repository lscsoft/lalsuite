/*
 * Copyright (C) 2007 S.Fairhurst, B. Krishnan, L.Santamaria, C. Robinson
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
 * \defgroup SphericalHarmonics Spin-weighted Spherical Harmonics
 * \ingroup support
 * \author S.Fairhurst, B. Krishnan, L.Santamaria, C. Robinson
 *
 * \brief Library of spherical harmonic functions
 *

 *
 */



/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

#ifndef _SPHERICALHARMONICS_H
#define _SPHERICALHARMONICS_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( SPHERICALHARMONICSH, "$Id$");

/**
 * Computes the (s)Y(l,m) spin-weighted spherical harmonic.
 *
 * From somewhere ....
 *
 * See also:
 * Implements Equations (II.9)-(II.13) of
 * D. A. Brown, S. Fairhurst, B. Krishnan, R. A. Mercer, R. K. Kopparapu,
 * L. Santamaria, and J. T. Whelan,
 * "Data formats for numerical relativity waves",
 * arXiv:0709.0093v1 (2007).
 *
 * Currently only supports s=-2, l=2,3,4,5,6,7,8 modes.
 */
COMPLEX16 XLALSpinWeightedSphericalHarmonic(
                                   REAL8 theta,  /**< polar angle (rad) */
                                   REAL8 phi,    /**< azimuthal angle (rad) */
                                   int s,        /**< spin weight */
                                   int l,        /**< mode number l */
                                   int m         /**< mode number m */
    );

/**
 * Computes the scalar spherical harmonic \f$ Y_{lm}(\theta, \phi) \f$.
 */
int
XLALScalarSphericalHarmonic(
                         COMPLEX16 *y, /**< output */
                         UINT4 l,      /**< value of l */
                         INT4  m,      /**< value of m */
                         REAL8 theta,  /**< angle theta */
                         REAL8 phi     /**< angle phi */
                         );

/**
 * Computes the spin 2 weighted spherical harmonic. This function is now
 * deprecated and will be removed soon. All calls should be replaced with
 * calls to XLALSpinWeightedSphHarm.
 */
INT4 XLALSphHarm ( COMPLEX16 *out, /**< output */
                   UINT4   L,      /**< value of L */
                   INT4 M,         /**< value of M */
                   REAL4 theta,    /**< angle with respect to the z axis */
                   REAL4   phi     /**< angle with respect to the x axis */
                   );

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif /* _SPHERICALHARMONICS_H */
