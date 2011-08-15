/*
*  Copyright (C) 2011 Stas Babak, David Churches, Drew Keppel, Duncan Brown, Jolien Creighton, David McKechan, Patrick Brady, Peter Shawhan, Reinhard Prix, B.S. Sathyaprakash, Anand Sengupta, Craig Robinson , Sean Seader, Thomas Cokelaer, Riccardo Sturani,  Laszlo Vereb
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

typedef struct
tagexpnCoeffsTaylorT3 {
   int ieta;
   /* Taylor expansion coefficents in phi(t)*/
   REAL8 ptaN, pta2, pta3, pta4, pta5, pta6, pta7, ptl6;
   /* Taylor expansion coefficents in f(t)*/
   REAL8 ftaN, fta2, fta3, fta4, fta5, fta6, fta7, ftl6;

   /* sampling interval*/
   REAL8 samplinginterval;
   /* symmetric mass ratio, total mass, component masses*/
   REAL8 eta, totalmass, m1, m2;
   /* unknown 3PN parameters, euler constant*/
   REAL8 lambda, theta, EulerC;

   /* initial and final values of frequency, time, velocity; lso
    values of velocity and frequency; final phase.*/
   REAL8 f0, fn, t0, tn, v0, vn, vf, vlso, flso, phiC;

   /* last stable orbit and pole defined by various Taylor and P-approximants*/
   REAL8 vlsoT0, vlsoT2, vlsoT4, vlsoT6;
}  expnCoeffsTaylorT3;


typedef REAL8 (SimInspiralPhasing3)(
   REAL8 td,
   expnCoeffsTaylorT3 *ak);

typedef REAL8 (SimInspiralFrequency3)(
   REAL8 td,
   expnCoeffsTaylorT3 *ak);

typedef struct
tagexpnFuncTaylorT3
{
   SimInspiralPhasing3 *phasing3;
   SimInspiralFrequency3 *frequency3;
} expnFuncTaylorT3;

REAL8 XLALSimInspiralChirpLength(
		REAL8 m1,		/**< mass of companion 1 */
		REAL8 m2,		/**< mass of companion 2 */
		REAL8 f_min,		/**< start frequency */
		int O			/**< twice post-Newtonian order */
		);

int XLALSimInspiralTaylorT3Setup(
	expnCoeffsTaylorT3 *ak,	/**< coefficients for TaylorT3 evolution [returned] */
	expnFuncTaylorT3 *f,	/**< functions for TaylorT3 evolution [returned] */
       	REAL8 deltaT,		/**< sampling interval */
	REAL8 m1,		/**< mass of companion 1 */
	REAL8 m2,		/**< mass of companion 2 */
	REAL8 f_min,		/**< start frequency */
	int O			/**< twice post-Newtonian order */
	);

