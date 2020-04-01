/*
* Copyright (C) 2019 Cecilio García Quirós, Geraint Pratten
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

#ifndef _LALSIM_IMR_PHENOMXPHM_H
#define _LALSIM_IMR_PHENOMXPHM_H


#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <complex.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/LALSimInspiral.h>


#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomX_precession.h"


/* Multibanding grid for precessing angles */
INT4 XLALSimIMRPhenomXPHMMultibandingGrid(REAL8Sequence **coarseFreqs, UINT4 ell, UINT4 emmprime, IMRPhenomXWaveformStruct *pWF, LALDict *lalParams);

#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMXPHM_H */
