/*
*  Copyright (C) 2011 Craig Robinson, Enrico Barausse
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
 * \author Craig Robinson, Collin Capano, Yi Pan
 *
 * \ file
 *
 * \brief Functions for producing EOB waveforms for
 * spinning binaries, as described in Barausse and Buonanno PRD 81, 084024 (2010).
 */

#ifndef _LALSIMIMREOBSPINAUXFUNCS_C
#define _LALSIMIMREOBSPINAUXFUNCS_C

#include <complex.h>
#include "LALSimIMREOBNRv2.h"

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

static int XLALSimIMRSpinEOBCalculateSigmaKerr(
                                   REAL8Vector *sigmaKerr,
                                   REAL8        mass1,
                                   REAL8        mass2,
                                   REAL8Vector *s1,
                                   REAL8Vector *s2 );

static int XLALSimIMRSpinEOBCalculateSigmaStar(
                                   REAL8Vector *sigmaStar,
                                   REAL8        mass1,
                                   REAL8        mass2,
                                   REAL8Vector *s1,
                                   REAL8Vector *s2 );

/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

/**
 * Function to calculate normalized spin of the deformed-Kerr background in SEOBNRv1.
 * Eq. 5.2 of Barausse and Buonanno PRD 81, 084024 (2010).
 */
static int XLALSimIMRSpinEOBCalculateSigmaKerr( 
                  REAL8Vector *sigmaKerr, /**<< OUTPUT, normalized (to total mass) spin of deformed-Kerr */
                  REAL8        mass1,     /**<< mass 1 */
                  REAL8        mass2,     /**<< mass 2 */
                  REAL8Vector *s1,        /**<< spin vector 1 */
                  REAL8Vector *s2         /**<< spin vector 2 */)
{

  UINT4 i;
  REAL8 totalMass = mass1 + mass2;

  for ( i = 0; i < 3; i++ )
  {
    sigmaKerr->data[i] = (s1->data[i] + s2->data[i])/(totalMass*totalMass);
  }

  return XLAL_SUCCESS;
}

/**
 * Function to calculate normalized spin of the test particle in SEOBNRv1.
 * Eq. 5.3 of Barausse and Buonanno PRD 81, 084024 (2010).
 */
static int XLALSimIMRSpinEOBCalculateSigmaStar( 
                  REAL8Vector *sigmaStar, /**<< OUTPUT, normalized (to total mass) spin of test particle */
                  REAL8        mass1,     /**<< mass 1 */
                  REAL8        mass2,     /**<< mass 2 */
                  REAL8Vector *s1,        /**<< spin vector 1 */
                  REAL8Vector *s2         /**<< spin vector 2 */)
{

  UINT4 i;
  REAL8 totalMass = mass1 + mass2;

  for ( i = 0; i < 3; i++ )
  {
    sigmaStar->data[i] = (mass2/mass1 * s1->data[i] + mass1/mass2 * s2->data[i])/(totalMass*totalMass);
  }

  return XLAL_SUCCESS;
}


#endif /* _LALSIMIMREOBSPINAUXFUNCS_C */
