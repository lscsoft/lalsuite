/*
 * Copyright (C) 2011 E. Ochsner
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
#include <lal/LALSimInspiral.h>
#include <lal/LALAdaptiveRungeKutta4.h>

NRCSID(LALSIMINSPIRALSPINTAYLORT4H, "$Id$");

/**
 * This function evolves the orbital equations for a precessing binary using 
 * the "TaylorT4" approximant for solving the orbital dynamics 
 * (see arXiv:0907.0700 for a review of the various PN approximants).
 *
 * It returns time series of the "orbital velocity", orbital phase, 
 * and components for both individual spin vectors, the "Newtonian"
 * orbital angular momentum (which defines the instantaneous plane)
 * and "E1", a basis vector in the instantaneous orbital plane.
 * Note that LNhat and E1 completely specify the instantaneous orbital plane.
 * It also returns the time and phase of the final time step
 *
 * FIXME: Do we want tc, phic or tStart, phiStart or both or something else?
 *
 * For input, the function takes the two masses, the initial orbital phase, 
 * Values of S1, S2, LNhat, E1 vectors at starting time,
 * the desired time step size, the starting GW frequency, 
 * and PN order at which to evolve the phase,
 * 
 * NOTE: All vectors are given in the so-called "radiation frame", 
 * where the direction of propagation is the z-axis, the principal "+" 
 * polarization axis is the x-axis, and the y-axis is given by the RH rule.
 * You must give the initial values in this frame, and the time series of the
 * vector components will also be returned in this frame
 */
int XLALSimInspiralPNEvolveOrbitSpinTaylorT4(
	REAL8TimeSeries **V,      /**< post-Newtonian parameter [returned]*/
	REAL8TimeSeries **Phi,    /**< orbital phase            [returned]*/
	REAL8TimeSeries **S1x,    /**< Spin1 vector x component [returned]*/
	REAL8TimeSeries **S1y,    /**< "    "    "  y component [returned]*/
	REAL8TimeSeries **S1z,    /**< "    "    "  z component [returned]*/
	REAL8TimeSeries **S2x,    /**< Spin2 vector x component [returned]*/
	REAL8TimeSeries **S2y,    /**< "    "    "  y component [returned]*/
	REAL8TimeSeries **S2z,    /**< "    "    "  z component [returned]*/
	REAL8TimeSeries **LNhatx, /**< unit orbital ang. mom. x [returned]*/
	REAL8TimeSeries **LNhaty, /**< "    "    "  y component [returned]*/
	REAL8TimeSeries **LNhatz, /**< "    "    "  z component [returned]*/
	REAL8TimeSeries **E1x,    /**< orb. plane basis vector x[returned]*/
	REAL8TimeSeries **E1y,    /**< "    "    "  y component [returned]*/
	REAL8TimeSeries **E1z,    /**< "    "    "  z component [returned]*/
	REAL8 m1,                 /**< mass of companion 1 */
	REAL8 m2,                 /**< mass of companion 2 */
	LIGOTimeGPS *tStart,      /**< start time of output vectors */
	REAL8 phiStart,           /**< orbital phase at initial time */
	REAL8 s1x,                /**< initial value of S1x */
	REAL8 s1y,                /**< initial value of S1y */
	REAL8 s1z,                /**< initial value of S1z */
	REAL8 s2x,                /**< initial value of S2x */
	REAL8 s2y,                /**< initial value of S2y */
	REAL8 s2z,                /**< initial value of S2z */
	REAL8 lnhatx,             /**< initial value of LNhatx */
	REAL8 lnhaty,             /**< initial value of LNhaty */
	REAL8 lnhatz,             /**< initial value of LNhatz */
	REAL8 e1x,                /**< initial value of E1x */
	REAL8 e1y,                /**< initial value of E1y */
	REAL8 e1z,                /**< initial value of E1z */
	REAL8 deltaT,             /**< sampling interval (s) */
	REAL8 fStart,             /**< start frequency */
	LALSpinInteraction spinFlags,  /**< flags to control spin effects */
	INT4 phaseO               /**< twice post-Newtonian order */
	);

