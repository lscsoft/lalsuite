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
 * \brief The code \ref LALInspiralPhasing2.c calculates the phase of an inspiral
 * waveform as a function of the
 * instantaneous frequency of the wave, up to \f$2^{nd}\f$ post--Newtonian order.
 *
 * ### Prototypes ###
 *
 * <tt>LALInspiralPhasing2()</tt>
 * <ul>
 * <li> \c phase: Output, the phase of the wave at the current epoch.</li>
 * <li> \c v: Input, the PN expansion parameter at the current epoch.</li>
 * <li> \c ak: Input containing PN expansion coefficients.</li>
 * </ul>
 *
 * ### Description ###
 *
 * The phase of the inspiral wave corresponding to the \c Approximant #TaylorT2
 * as in \eqref{eq_InspiralWavePhase2} (<tt>correct equation?</tt>)
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
XLALInspiralPhasing2_0PN (
   REAL8       v,
   expnCoeffs *ak
   )
{
  REAL8 v5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  v5 = pow(v,5.);
  phase = ak->phiC
         + ak->pvaN / v5;

  return phase;
}



REAL8
XLALInspiralPhasing2_2PN (
   REAL8       v,
   expnCoeffs *ak
   )
{
  REAL8 v2,v5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  v2 = v*v;
  v5 = v2*v2*v;
  phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2);

  return phase;
}



REAL8
XLALInspiralPhasing2_3PN (
   REAL8       v,
   expnCoeffs *ak
   )
{
  REAL8 v2,v3,v5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  v2 = v*v;
  v3 = v2*v;
  v5 = v3*v2;
  phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2
         + ak->pva3 * v3);

  return phase;
}



REAL8
XLALInspiralPhasing2_4PN (
   REAL8       v,
   expnCoeffs *ak
   )
{
  REAL8 v2,v3,v4,v5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2
         + ak->pva3 * v3
         + ak->pva4 * v4);

  return phase;
}



REAL8
XLALInspiralPhasing2_5PN (
   REAL8       v,
   expnCoeffs *ak
   )
{
  REAL8 v2,v3,v4,v5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2
         + ak->pva3 * v3
         + ak->pva4 * v4
         + ak->pva5 * log(v/ak->vlso) * v5);

  return phase;
}



REAL8
XLALInspiralPhasing2_6PN (
   REAL8       v,
   expnCoeffs *ak
   )
{
  REAL8 v2,v3,v4,v5,v6;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2
         + ak->pva3 * v3
         + ak->pva4 * v4
         + ak->pva5 * log(v/ak->vlso) * v5
         + (ak->pva6 + ak->pvl6*log(4*v)) * v6);

  return phase;
}



REAL8
XLALInspiralPhasing2_7PN (
   REAL8       v,
   expnCoeffs *ak
   )
{
  REAL8 v2,v3,v4,v5,v6,v7;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  v7 = v6*v;
  phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2
         + ak->pva3 * v3
         + ak->pva4 * v4
         + ak->pva5 * log(v/ak->vlso) * v5
         + (ak->pva6 + ak->pvl6*log(4*v)) * v6
         + ak->pva7 * v7);

  return phase;
}
