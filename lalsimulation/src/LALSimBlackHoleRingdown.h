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

REAL8 XLALSimRadiusKerrISCO ( REAL8 a );
REAL8 XLALSimEnergyKerrISCO ( REAL8 rISCO );
REAL8 XLALSimAngMomKerrISCO ( REAL8 rISCO );

/* Constants entering the final mass formulas of SEOBNRv2,4 */
/* See http://arxiv.org/pdf/1206.3803.pdf */
static const REAL8 p0coeff = 0.04826;
static const REAL8 p1coeff = 0.01559;
static const REAL8 p2coeff = 0.00485;
/* See http://arxiv.org/pdf/0904.2577.pdf */
static const REAL8 t0coeff = -2.8904;
static const REAL8 t2coeff = -3.5171;
static const REAL8 t3coeff = 2.5763;
static const REAL8 s4coeff = -0.1229;
static const REAL8 s5coeff = 0.4537;
/* See https://dcc.ligo.org/T1400476 */
static const REAL8 s9coeff = 2.763032781169752;
static const REAL8 s8coeff = -2.6081232221537394;
static const REAL8 s7coeff = 1.2657111864932808;
static const REAL8 s6coeff = -0.7835007857591175;
static const REAL8 s5v2coeff = -0.3264724801557159;
static const REAL8 s4v2coeff = -0.27506210736300474;
static const REAL8 t0v2coeff = -2.649826989941522;
static const REAL8 t3v2coeff = 3.910637513328723;
static const REAL8 t2v2coeff = -3.850983155206041;

/* Constants entering the final spin formulas of SEOBNRv4 */
/* Table I of https://arxiv.org/pdf/1605.01938v2.pdf */
static const REAL8 k00 = -5.977230835551017; // Solving Eq.(11) of https://arxiv.org/pdf/1605.01938v2.pdf
static const REAL8 k01 = 3.39221;
static const REAL8 k02 = 4.48865;
static const REAL8 k03 = -5.77101;
static const REAL8 k04 = -13.0459;
static const REAL8 k10 = 35.1278;
static const REAL8 k11 = -72.9336;
static const REAL8 k12 = -86.0036;
static const REAL8 k13 = 93.7371;
static const REAL8 k14 = 200.975;
static const REAL8 k20 = - 146.822;
static const REAL8 k21 = 387.184;
static const REAL8 k22 = 447.009;
static const REAL8 k23 = -467.383;
static const REAL8 k24 = -884.339;
static const REAL8 k30 = 223.911;
static const REAL8 k31 = -648.502;
static const REAL8 k32 = -697.177;
static const REAL8 k33 = 753.738;
static const REAL8 k34 = 1166.89;

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMBLACKHOLERINGDOWN_H */
