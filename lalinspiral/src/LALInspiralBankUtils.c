/*
*  Copyright (C) 2007 Thomas Cokelaer
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
 * \defgroup LALInspiralBankUtils_c Module LALInspiralBankUtils.c
 * \ingroup LALInspiralBank_h
 * \author Cokelaer Thomas
 *
 * \brief NONE
 *
 * ### Description ###
 *
 * In a parameter space defined by \f$m_1\f$ and \f$m_2\f$, or equivalently, \f$M=m_1+m_2\f$ and \f$\eta=\frac{m_1 m_2}{M^2}\f$, the conversion
 * to chirp-time parameter such as \f$\tau_0\f$ and \f$\tau_3\f$ si quite common. In particular, it is interesting to get the value of
 * \f$\tau_3\f$ when only \f$\tau_0\f$ is known, and a constraint on the masses exists (e.g., \f$m_1=m_2\f$ or one of the mass equals mMin or mMax.
 * This modules contains a few functions to perform these conversion.
 *
 * ### Algorithm ###
 *
 * We know that
 * \anchor eq_tau0a \f{equation}{
 * \tau_0 = \frac{A_0}{\eta} M^{-5/2},
 * \tag{eq_tau0a}
 * \f}
 * and
 * \f{equation}{
 * \tau_3 = \frac{A_3}{\eta} M^{-2/3},
 * \f}
 * where
 * \f{equation}{
 * A_0 = \frac{5}{256 (\pi *f_L)^{8/3}},
 * \f}
 * and
 * \f{equation}{
 * A_3 = \frac{\pi}{8 (\pi *f_L)^{5/3}},
 * \f}
 *
 * Therefore, it is straightforward to express \f$\tau_3\f$ as a function of \f$\tau_0\f$ amd \f$\eta\f$:
 * \anchor eq_tau3b \f{equation}{
 * \tau_3 = \frac{A3}{\eta} \left( \frac{\tau_0 \eta}{ A_0} \right)^{2/5}
 * \tag{eq_tau3b}
 * \f}
 * if \f$\eta=0.25\f$ on the equal-mass line, then
 * \anchor eq_tau3a \f{equation}{
 * \tau_3 = 4 A3 \left( \frac{\tau_0}{ 4 A_0} \right)^{2/5}
 * \tag{eq_tau3a}
 * \f}
 *
 * Equation\eqref{eq_tau3b} returns \f$\tau_3\f$ given in \f$M, \eta\f$ and \f$f_L\f$ and is defined
 * in\c XLALInspiralTau3FromNonEqualMass().
 *
 * Equation\eqref{eq_tau3a} returns tau3 in the particular case \f$m_1=m_2\f$, given
 * \f$\tau_0\f$ only, and is defined in \c XLALInspiralTau3FromTau0AndEqualMassLine().
 *
 * Equation\eqref{eq_tau0a} returns \f$tau_0\f$ given \f$M, \eta\f$ and \f$f_L\f$, and is defined
 * \c XLALInspiralTau0FromMEta().
 *
 * Finally, \c XLALInspiralMFromTau0AndNonEqualMass() returns \f$M\f$ when \f$\tau_0\f$ is known
 * and a constraint exists on one of the individual mass (e.g., \f$m_1=\textrm{mMax}\f$ or
 * \f$m_1=\textrm{mMin}\f$). This functions requires a little more algebra and is used in the
 * HybridHexagonal placement. The \ref LALInspiralHybridHexagonalBank_c describes this algebra.
 *
 */
/*@{*/

#include <stdio.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FindRoot.h>
#include <lal/LALInspiralBank.h>
#include <math.h>


/** \see See \ref LALInspiralBankUtils_c for documentation */
REAL4
XLALInspiralTau3FromTau0AndEqualMassLine(
    REAL4               tau0,
    REAL4               fL
    )
{
  REAL4 A0, A3, tau3=0;


  A0 = (5.0 / 256.0) * pow(LAL_PI * fL, (-8.0/3.0));
  A3  = LAL_PI / (8.0 * pow(LAL_PI * fL, (5.0/3.0)));

  tau3 = 4 * A3 * pow(tau0/4/A0, 2./5.);

  return tau3;
}


/** \see See \ref LALInspiralBankUtils_c for documentation */
REAL4
XLALInspiralTau3FromNonEqualMass(
    REAL4              	M,
    REAL4 		eta,
    REAL4		fL
 )
{
  REAL4 A3;
  REAL4 tau3 = 0;

  A3  = LAL_PI / (8.0 * pow(LAL_PI*fL, (5.0/3.0)));
  tau3 = A3 * pow(M * LAL_MTSUN_SI, -2.0/3.0) / eta;

  return tau3;
}

/** \see See \ref LALInspiralBankUtils_c for documentation */
REAL4
XLALInspiralTau0FromMEta(
    REAL4              	M,
    REAL4 		eta,
    REAL4		fL
 )
{

/* This function returns tau3, computed from M and eta*/

  REAL4 A0;
  REAL4 tau0 = 0;

  A0 = (5.0 / 256.0) * pow( LAL_PI * fL, (-8.0/3.0));
  tau0 = A0 * pow(M*LAL_MTSUN_SI, -5.0/3.0) / eta;

  return tau0;
}

/** \see See \ref LALInspiralBankUtils_c for documentation */
REAL8
XLALInspiralMFromTau0AndNonEqualMass(
  REAL8 tau0,
  REAL8 extremMass,
  REAL8 fL)
{
  REAL8 result, A0, p, q, x;

  A0 = (5.0 / 256.0) * pow( LAL_PI * fL, (-8.0/3.0));

  /* from tau0, and M, we can get a poylomial expression where M is the
  unknowm of the form x^3+px+q =0 where x = M^(1/3) and p and q as follows :*/
  p = -A0/tau0/extremMass/LAL_MTSUN_SI;
  q = -extremMass * LAL_MTSUN_SI;

  x = pow((-q/2-0.5*sqrt((27*q*q + 4*p*p*p)/27)), 1./3.);
  x += pow((-q/2+0.5*sqrt((27*q*q + 4*p*p*p)/27)), 1./3.);

  /* This is a real solution and M is simply */
  result = x*x*x/LAL_MTSUN_SI;

  return result;
}
/*@}*/
