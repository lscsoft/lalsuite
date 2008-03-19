/*
 * Copyright (C) 2008 J. Creighton, S. Fairhurst, B. Krishnan, L. Santamaria
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
#include <lal/LALDatatypes.h>

#define LAL_PN_MODE_L_MAX 3 


/**
 * Computes the (s)Y(l,m) spin-weighted spherical harmonic.
 *
 * From somewhere ...
 *
 * See also:
 * Equations (II.9)-(II.13) of
 * D. A. Brown, S. Fairhurst, B. Krishnan, R. A. Mercer, R. K. Kopparapu,
 * L. Santamaria, and J. T. Whelan,
 * "Data formats for numerical relativity waves",
 * arXiv:0709.0093v1 (2007).
 *
 * Currently only supports s=-2, l=2,3,4,5 modes.
 *
 */
COMPLEX16 XLALSpinWeightedSphericalHarmonic(REAL8 theta, REAL8 phi, int s, int l, int m);

/**
 * Multiplies a mode h(l,m) by a spin-2 weighted spherical harmonic
 * to obtain hplus - i hcross, which is added to the time series.
 *
 * Implements the sum of a single term of Eq. (11) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 *
 * If sym is non-zero, symmetrically add the m and -m terms assuming
 * that h(l,-m) = (-1)^l h(l,m)*; see Eq. (78) ibid.
 */
int XLALSimAddMode(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, COMPLEX16TimeSeries *hmode, REAL8 theta, REAL8 phi, int l, int m, int sym);

/**
 * Computes the rate of increase of the orbital frequency for a post-Newtonian
 * inspiral.  This function returns dx/dt rather than the true angular
 * acceleration.
 *
 * Implements Equation (6) of
 * Yi Pan, Alessandra Buonanno, Yanbei Chen, and Michele Vallisneri,
 * "A physical template family for gravitational waves from precessing
 * binaries of spinning compact objects: Application to single-spin binaries"
 * arXiv:gr-qc/0310034v3 (2007).
 *
 * Note: this equation is actually dx/dt rather than (domega/dt)/(omega)^2
 * so the leading coefficient is different.  Also, this function applies
 * for non-spinning objects.
 *
 * Compare the overall coefficient, with nu=1/4, to Equation (45) of
 * Michael Boyle, Duncan A. Brown, Lawrence E. Kidder, Abdul H. Mroue, 
 * Harald P. Pfeiﬀer, Mark A. Scheel, Gregory B. Cook, and Saul A. Teukolsky 
 * "High-accuracy comparison of numerical relativity simulations with
 * post-Newtonian expansions"
 * arXiv:0710.0158v1 (2007).
 */
REAL8 XLALSimInspiralPNAngularAcceleration(REAL8 x, REAL8 m1, REAL8 m2, int O);

/**
 * Computes the orbital angular velocity from the quantity x.
 * This is from the definition of x.
 *
 * Implements Equation (46) of
 * Michael Boyle, Duncan A. Brown, Lawrence E. Kidder, Abdul H. Mroue, 
 * Harald P. Pfeiﬀer, Mark A. Scheel, Gregory B. Cook, and Saul A. Teukolsky 
 * "High-accuracy comparison of numerical relativity simulations with
 * post-Newtonian expansions"
 * arXiv:0710.0158v1 (2007).
 */
REAL8 XLALSimInspiralPNAngularVelocity(REAL8 x, REAL8 m1, REAL8 m2);

/**
 * Computes the orbital energy at a fixed frequency and pN order.
 *
 * Implements Equation (152) of
 * Luc Blanchet,
 * "Gravitational Radiation from Post-Newtonian Sources and Inspiralling
 * Compact Binaries",
 * http://www.livingreviews.org/lrr-2006-4/index.html
 *
 * This is the same as Equation (10) (where the spin of the objects
 * is zero) of:
 * Yi Pan, Alessandra Buonanno, Yanbei Chen, and Michele Vallisneri,
 * "A physical template family for gravitational waves from precessing
 * binaries of spinning compact objects: Application to single-spin binaries"
 * arXiv:gr-qc/0310034v3 (2007). 
 * Note: this equation is actually dx/dt rather than (domega/dt)/(omega)^2
 * so the leading coefficient is different.
 */
REAL8 XLALSimInspiralPNEnergy(REAL8 x, REAL8 m1, REAL8 m2, int O);

/**
 * Evolves a post-Newtonian orbit using the Taylor T4 method.
 *
 * See:
 * Michael Boyle, Duncan A. Brown, Lawrence E. Kidder, Abdul H. Mroue, 
 * Harald P. Pfeiﬀer, Mark A. Scheel, Gregory B. Cook, and Saul A. Teukolsky 
 * "High-accuracy comparison of numerical relativity simulations with
 * post-Newtonian expansions"
 * arXiv:0710.0158v1 (2007).
 */
int XLALSimInspiralPNEvolveOrbitTaylorT4(REAL8TimeSeries **x, REAL8TimeSeries **phi, LIGOTimeGPS *tc, REAL8 phic, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fmin, int O);

/**
 * Computes h(l,m) modes of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (79)-(116) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16 XLALSimInspiralPNMode22(REAL8 x, REAL8 phi, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16 XLALSimInspiralPNMode21(REAL8 x, REAL8 phi, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16 XLALSimInspiralPNMode33(REAL8 x, REAL8 phi, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16 XLALSimInspiralPNMode32(REAL8 x, REAL8 phi, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16 XLALSimInspiralPNMode31(REAL8 x, REAL8 phi, REAL8 m1, REAL8 m2, REAL8 r, int O);

/**
 * Computes h(l,m) mode timeseries of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * See Eqns. (79)-(116) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(REAL8TimeSeries *x, REAL8TimeSeries *phi, REAL8 m1, REAL8 m2, REAL8 r, int O, int l, int m);

/**
 * Given an orbit evolution phasing, construct the waveform h+ and hx.
 *
 * Implements Equation (11) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
int XLALSimInspiralPNPolarizationWaveforms(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8TimeSeries *x, REAL8TimeSeries *phi, REAL8 m1, REAL8 m2, REAL8 r, REAL8 i, int O);

/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 * 
 * This routine allows the user to specify different pN orders
 * for phasing calcuation vs. amplitude calculations.
 */
int XLALSimInspiralPNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LIGOTimeGPS *tc, REAL8 phic, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fmin, REAL8 r, REAL8 i, int amplitudeO, int phaseO);

/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 * 
 * This routine uses the same pN order for phasing and amplitude
 * (unless the order is -1 in which case the highest available
 * order is used for both of these -- which might not be the same).
 */
int XLALSimInspiralPN(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LIGOTimeGPS *tc, REAL8 phic, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fmin, REAL8 r, REAL8 i, int O);

/**
 * Driver routine to compute the restricted post-Newtonian inspiral waveform.
 * 
 * This routine computes the phasing to the specified order, but
 * only computes the amplitudes to the Newtonian (quadrupole) order.
 */
int XLALSimInspiralPNRestricted(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LIGOTimeGPS *tc, REAL8 phic, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fmin, REAL8 r, REAL8 i, int O);
