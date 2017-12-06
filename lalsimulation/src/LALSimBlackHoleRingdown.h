/*
 * Copyright (C) 2011 J. Creighton
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

#ifndef _LALSIMBLACKHOLERINGDOWN_H
#define _LALSIMBLACKHOLERINGDOWN_H

#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/**
 * @author Jolien Creighton
 * @addtogroup LALSimBlackHoleRingdown_h Header LALSimBlackHoleRingdown.h
 * @ingroup lalsimulation_inspiral
 * @brief Routines to generate black hole ringdown waveforms.
 * @details
 * These routines generate black hole quasinormal modes, spin-weighted
 * spheroidal harmonics, and ringdown gravitational waveforms.
 */

/* LOW-LEVEL ROUTINES (USE LEAVER'S CONVENSIONS G = c = 2M = 1) */

int XLALSimBlackHoleRingdownModeEigenvaluesLeaver(COMPLEX16 *A, COMPLEX16 *omega, double a, int l, int m, int s
);
COMPLEX16 XLALSimBlackHoleRingdownSpheroidalWaveFunctionLeaver(double mu, double a, int l, int m, int s, COMPLEX16 A, COMPLEX16 omega);


/* HIGH-LEVEL ROUTINES */

int XLALSimBlackHoleRingdownMode(double *frequency, double *quality, double mass, double dimensionless_spin, int l, int m, int s);
COMPLEX16 XLALSimBlackHoleRingdownSpheroidalWaveFunction(double theta, double dimensionless_spin, int l, int m, int s);
int XLALSimBlackHoleRingdown(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const LIGOTimeGPS *t0, double phi0, double deltaT, double mass, double dimensionless_spin, double fractional_mass_loss, double distance, double inclination, int l, int m);
INT4 XLALSimIMREOBFinalMassSpin(REAL8 *finalMass, REAL8 *finalSpin, const REAL8 mass1, const REAL8 mass2, const REAL8 spin1[3], const REAL8 spin2[3], Approximant approximant);
INT4 XLALSimIMREOBGenerateQNMFreqV2(COMPLEX16Vector *modefreqs, const REAL8 mass1, const REAL8 mass2, const REAL8 spin1[3], const REAL8 spin2[3], UINT4 l, INT4 m, UINT4 nmodes, Approximant approximant);
INT4 XLALSimIMREOBGenerateQNMFreqV2fromFinal(COMPLEX16Vector *modefreqs, const REAL8 finalMass, const REAL8 finalSpin, UINT4 l, INT4 m, UINT4 nmodes);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMBLACKHOLERINGDOWN_H */
