/*
*  Copyright (C) 2007 Duncan Brown
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpPTF.h
 *
 * Author: Brown, D. A., and Fazi, D.
 *
 *-----------------------------------------------------------------------
 */

#ifndef _FINDCHIRPPTFH_H
#define _FINDCHIRPPTFH_H

#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInspiral.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \defgroup FindChirpPTF_h Header FindChirpPTF.h
 * \ingroup lalinspiral_findchirp
 * \author Brown, D. A., and Fazi, D.
 *
 * \brief Provides structures and functions to filter interferometer data using the
 * physical template family.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/FindChirpPTF.h>
 * \endcode
 *
 */
/** @{ */

REAL4Vector*
XLALPTFOmegaPNCoeffsOrbital(
    REAL4                m1,
    REAL4                m2
    );

REAL4Vector*
XLALPTFOmegaPNCoeffsSpin(
    REAL4                m1,
    REAL4                m2,
    REAL4                chi1,
    REAL4                chi2,
    REAL4                Q1,
    REAL4                Q2
    );

REAL4Vector*
XLALPTFOmegaPNCoeffsEnergy(
    REAL4                m1,
    REAL4                m2,
    REAL4                chi1,
    REAL4                chi2,
    REAL4                Q1,
    REAL4                Q2
    );

INT4
XLALFindChirpPTFWaveform(
    REAL4Vector         *PTFphi,
    REAL4Vector         *PTFomega_2_3,
    REAL4VectorSequence *PTFe1,
    REAL4VectorSequence *PTFe2,
    InspiralTemplate    *InspTmplt,
    REAL8                deltaT
    );

/** @} */ /* end:FindChirpPTF.h */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FINDCHIRPPTFH_H */
