#ifndef _LALSIM_IMR_PHENOMX_H
#define _LALSIM_IMR_PHENOMX_H

/*
 * Copyright (C) 2018 Geraint Pratten
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
 * \author Geraint Pratten
 *
 */

#ifdef __cplusplus
extern "C" {
#endif


/* CONSTANTS */
/* Dimensionless frequency (Mf) at which define the end of the waveform */
#define f_CUT 0.3

#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"

#include "LALSimIMRPhenomX_precession.h"

/* Decleration for internal function to generate aligned-spin, 22 only IMRPhenomXAS waveform */
int IMRPhenomXASGenerateFD(
  COMPLEX16FrequencySeries **htilde22,
  const REAL8Sequence *freqs,
  IMRPhenomXWaveformStruct *pWF,
  LALDict *lalParams
);


int IMRPhenomXCheckForUniformFrequencies(REAL8Sequence *frequencies,REAL8 df);

#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMX_H */
