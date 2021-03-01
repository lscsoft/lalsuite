/*
*  Copyright (C) 2007 David Churches, B.S. Sathyaprakash, Drew Keppel
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

/**
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief The code \ref LALInspiralPhasing3.c calculates the phase the waveform
 * from an inspiralling binary system as a function of time up to second post-Nowtonian
 * order.
 *
 * <tt>LALInspiralPhasing3()</tt>
 * <ul>
 * <li> \c phase: Output, the phase of the wave at the current epoch.</li>
 * <li> \c td: Input, the PN expansion coefficients of phase \f$\phi^t_k\f$ as a function
 * of time (cf. \ref table_flux "this table".</li>
 * <li> \c ak: Input containing PN expansion coefficients.</li>
 * </ul>
 *
 * ### Description ###
 *
 * The phase of the inspiral wave corresponding to the \c Approximant \c TaylorT2
 * as in \eqref{eq_InspiralWavePhase3}.
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * None.
 *
 * ### Notes ###
 *
 * None.
 *
 */

#include <math.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALInspiral.h>
#include <lal/XLALError.h>


REAL8
XLALInspiralPhasing3_0PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta5 = pow(td,-0.625);
  phase = (ak->ptaN/theta5);

  return phase;
}



REAL8
XLALInspiralPhasing3_2PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta5 = theta2*theta2*theta;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2);

  return phase;
}



REAL8
XLALInspiralPhasing3_3PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3,theta5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta5 = theta2*theta3;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3);

  return phase;
}



REAL8
XLALInspiralPhasing3_4PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4);

  return phase;
}



REAL8
XLALInspiralPhasing3_5PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4
         + ak->pta5 * log(td/ak->tn) * theta5);

  return phase;
}



REAL8
XLALInspiralPhasing3_6PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4
         + ak->pta5*log(td/ak->tn)*theta5
         +(ak->ptl6*log(td/256.) + ak->pta6)*theta6);

  return phase;
}



REAL8
XLALInspiralPhasing3_7PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6,theta7;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;
  theta7 = theta6*theta;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4
         + ak->pta5*log(td/ak->tn)*theta5
         +(ak->ptl6*log(td/256.) + ak->pta6)*theta6
         + ak->pta7*theta7);

  return phase;
}
