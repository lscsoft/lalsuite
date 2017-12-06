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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
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

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup FindChirpPTF_h
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
/*@{*/

void
LALFindChirpPTFTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *InspTmplt,
    FindChirpTmpltParams       *params
    );

void
LALFindChirpPTFNormalize(
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    FindChirpSegment           *fcSeg,
    FindChirpDataParams        *params
    );

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

INT4
XLALPTFWaveformDerivatives(
              REAL8  t,
        const REAL8  y[14],
              REAL8  dydt[14],
              void  *params );

void
LALFindChirpPTFFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );

/*@}*/ /* end:FindChirpPTF.h */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FINDCHIRPPTFH_H */
