//
// Copyright (C) 2014 Reinhard Prix
// Copyright (C) 2012, 2013, 2014 David Keitel, Bernd Machenschalk, Reinhard Prix, Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

#ifndef _ESTIMATEAMPLITUDEPARAMS_H
#define _ESTIMATEAMPLITUDEPARAMS_H

#include <lal/LALStdlib.h>
#include <lal/PulsarDataTypes.h>
#include <lal/LALComputeAM.h>
#include <lal/LALComputeAM.h>
#include <lal/SSBtimes.h>

#ifdef  __cplusplus
extern "C" {
#endif

///
/// \defgroup EstimateAmplitudeParams_h Header EstimateAmplitudeParams.h
/// \ingroup lalpulsar_coh
/// \authors Reinhard Prix
///
/// \brief Functions to estimate amplitude parameters and convert between different parametrizations.
/// * ### Synopsis ###
///
/// \code
/// #include <lal/EstimateAmplitudeParams.h>
/// \endcode
///

// @{

// ---------- API function prototypes ----------
int XLALEstimatePulsarAmplitudeParams ( PulsarCandidate *pulsarParams, const LIGOTimeGPS* FaFb_refTime,
                                        const COMPLEX8 Fa, const COMPLEX8 Fb, const AntennaPatternMatrix *Mmunu );

int XLALAmplitudeParams2Vect ( PulsarAmplitudeVect A_Mu, const PulsarAmplitudeParams Amp );
int XLALAmplitudeVect2Params( PulsarAmplitudeParams *Amp, const PulsarAmplitudeVect A_Mu );

// @}

#ifdef  __cplusplus
}
#endif

#endif // _ESTIMATEAMPLITUDEPARAMS_H
