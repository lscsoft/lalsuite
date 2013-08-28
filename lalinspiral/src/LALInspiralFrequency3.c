/*
*  Copyright (C) 2007 David Churches, B.S. Sathyaprakash, Duncan Brown, Drew Keppel
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
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief The code \ref LALInspiralFrequency3.c calculates the frequency the
 * waveform from an inspiralling binary system as a function of time up to 3.5
 * post-Nowtonian order.
 *
 * \heading{Prototypes}
 *
 * <tt>LALInspiralFrequency3()</tt>
 * <ul>
 * <li> \c frequency: Output containing the inspiral waveform.</li>
 * <li> \c td: Input containing PN expansion coefficients \f$F_k\f$ (cf. Table\tableref{table_flux})
 * of frequency as a function of time.</li>
 * <li> \c ak: Input containing all PN expansion coefficients.</li>
 * </ul>
 *
 * \heading{Description}
 *
 * This module computes the instantaneous frequency of an inspiral wave using
 * \f{equation}{
 * F(t) = F_N(\theta) \sum F_k \theta^k,
 * \f}
 * where the expansion coefficients \f$F_k,\f$ Newtonian value \f$F_N\f$ and the
 * time-variable \f$\theta\f$ are defined in Table\tableref{table_flux}.
 *
 * \heading{Algorithm}
 *
 * \heading{Uses}
 * None.
 *
 * \heading{Notes}
 * The frequency evolution defined by post-Newtonian expansion is not monotonic.
 * Indeed, the equations become highly inaccurate close to the last stable orbit (lso)
 * and breakdown at or slightly after lso, and the frequency begins to decrease at later times.
 * It turns out that the evolution is monotonic at least up to lso.
 *
 */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

void
LALInspiralFrequency3_0PN (
   LALStatus  *status,
   REAL8      *frequency,
   REAL8      td,
   expnCoeffs *ak
   )
{
  XLALPrintDeprecationWarning("LALInspiralFrequency3_0PN", "XLALInspiralFrequency3_0PN");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT(frequency, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  *frequency = XLALInspiralFrequency3_0PN(td, ak);
  if (XLAL_IS_REAL8_FAIL_NAN(*frequency))
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



REAL8
XLALInspiralFrequency3_0PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta3;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta3 = theta*theta*theta;

  frequency = theta3*ak->ftaN;

  return frequency;
}



void
LALInspiralFrequency3_2PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{
  XLALPrintDeprecationWarning("LALInspiralFrequency3_2PN", "XLALInspiralFrequency3_2PN");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT(frequency, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  *frequency = XLALInspiralFrequency3_2PN(td, ak);
  if (XLAL_IS_REAL8_FAIL_NAN(*frequency))
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



REAL8
XLALInspiralFrequency3_2PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2);

  return frequency;
}



void
LALInspiralFrequency3_3PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{
  XLALPrintDeprecationWarning("LALInspiralFrequency3_3PN", "XLALInspiralFrequency3_3PN");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT(frequency, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  *frequency = XLALInspiralFrequency3_3PN(td, ak);
  if (XLAL_IS_REAL8_FAIL_NAN(*frequency))
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



REAL8
XLALInspiralFrequency3_3PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3);

  return frequency;
}



void
LALInspiralFrequency3_4PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{
  XLALPrintDeprecationWarning("LALInspiralFrequency3_4PN", "XLALInspiralFrequency3_4PN");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT(frequency, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  *frequency = XLALInspiralFrequency3_4PN(td, ak);
  if (XLAL_IS_REAL8_FAIL_NAN(*frequency))
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



REAL8
XLALInspiralFrequency3_4PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3,theta4;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4);

  return frequency;
}



void
LALInspiralFrequency3_5PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{
  XLALPrintDeprecationWarning("LALInspiralFrequency3_5PN", "XLALInspiralFrequency3_5PN");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT(frequency, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  *frequency = XLALInspiralFrequency3_5PN(td, ak);
  if (XLAL_IS_REAL8_FAIL_NAN(*frequency))
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



REAL8
XLALInspiralFrequency3_5PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5);

  return frequency;
}



void
LALInspiralFrequency3_6PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{
  XLALPrintDeprecationWarning("LALInspiralFrequency3_6PN", "XLALInspiralFrequency3_6PN");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT(frequency, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  *frequency = XLALInspiralFrequency3_6PN(td, ak);
  if (XLAL_IS_REAL8_FAIL_NAN(*frequency))
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



REAL8
XLALInspiralFrequency3_6PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5
             + (ak->fta6 + ak->ftl6*log(td))*theta6);

  return frequency;
}



void
LALInspiralFrequency3_7PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{
  XLALPrintDeprecationWarning("LALInspiralFrequency3_7PN", "XLALInspiralFrequency3_7PN");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT(frequency, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  *frequency = XLALInspiralFrequency3_7PN(td, ak);
  if (XLAL_IS_REAL8_FAIL_NAN(*frequency))
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



REAL8
XLALInspiralFrequency3_7PN (
   REAL8       td,
   expnCoeffs *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6,theta7;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;
  theta7 = theta6*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5
             + (ak->fta6 + ak->ftl6*log(td))*theta6
             + ak->fta7*theta7);

  return frequency;
}
