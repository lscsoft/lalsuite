/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Craig Robinson
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

#include <lal/LALInspiralBank.h>

/** \ingroup LALInspiralBank_h
 * \brief This function takes as input \f$\tau_{0}\f$, \f$\tau_{3}\f$ and \f$f_a\f$ (the lower frequency of the detectors
 * sensitivity), it then calculates \f$m\f$ (the total mass of the binary), \f$\eta\f$ (the
 * symmetric mass ratio) and the individual mass of the compact objects.
 * \author Churches, D. K.
 *
 * We start with the definition of the chirp times \f$\tau_{0}\f$ and \f$\tau_{3}\f$,
 * \f{equation}{
 * \tau_{0} = \frac{5}{256 (\pi f_{a} )^{8/3} m^{5/3} \eta}
 * \f}
 * and
 * \f{equation}{
 * \tau_{3} = \frac{1}{8 (\pi^{2} f_{a}^{5} )^{1/3} m^{2/3} \eta}
 * \f}
 * These equations may be inverted to yield
 * \f{equation}{
 * m = \frac{5}{32 \pi^{2} f_{a}} \frac{\tau_{3}}{\tau_{0}}
 * \f}
 * and
 * \f{equation}{
 * \eta = \left( \frac{2 \pi^{2}}{25 f_{a}^{3}} \frac{\tau_{0}^{2}}{\tau_{3}^{5}}
 * \right)^{1/3}\f}
 * The individual masses may be calculated as follows.  We have
 *
 * \anchor mass1 \f{equation}{
 * m = m_{1} + m_{2}
 * \tag{mass1}
 * \f}
 * and
 * \anchor eta1 \f{equation}{
 * \eta = \frac{m_{1} m_{2}}{(m_{1} + m_{2})^{2}}
 * \tag{eta1}
 * \f}
 * From Eq.\eqref{mass1} we may eliminate either \f$m_{1}\f$ or \f$m_{2}\f$,
 * \f{equation}{
 * m_{1} = m - m_{2}
 * \f}
 * This may be substituted into Eq.\eqref{eta1} to give
 * \f{equation}{
 * \eta = \frac{(m - m_{2}) m_{2}}{\left[ (m - m{2}) + m_{2} \right]^{2}}
 * \f}
 * which may be re--arranged to give
 * \f{equation}{
 * m_{2}^{2} - m m_{2} + \eta m^{2} = 0
 * \f}
 * i.e.                                          \
 * \f{equation}{
 * m_{2} = \frac{ m \pm \sqrt{m^{2}(1 - 4 \eta) }}{2}
 * \f}
 * Therefore, since we know that \f$\eta \leq 1/4\f$, real roots are guaranteed.
 * If we had eliminated \f$m_{2}\f$ rather than \f$m_{1}\f$ then we would have arrived at an identical
 * expression for
 * \f$m_{1}\f$, and so of one object has mass
 * \f{equation}{
 * m_{1} = \frac{m + \sqrt{m^{2}(1-4 \eta)}}{2}
 * \f}
 * then the other object must have mass
 * \f{equation}{
 * m_{2} = \frac{m - \sqrt{m^{2}(1-4 \eta)}}{2}
 * \f}
 */
void LALInspiralComputeParams(LALStatus            *status,	/**< LAL status pointer */
                              InspiralTemplate     *pars,	/**< [out] inspiral wave parameter structure */
                              InspiralBankParams   bankParams,	/**< [in] the parameters of the template bank */
                              InspiralCoarseBankIn coarseIn	/**< [in] input parameters specifying the coarse bank */
                              )
{


  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);
  ASSERT (pars,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

  pars->fLower = coarseIn.fLower;

  ASSERT (bankParams.x0 > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT (bankParams.x1 > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
/*ASSERT (bankParams.x1 < bankParams.x0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);*/
  ASSERT ((INT4)coarseIn.space >= 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT ((INT4)coarseIn.space <= 1, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

  switch (coarseIn.space) {
     case Tau0Tau2:
        pars->t0 = bankParams.x0;
        pars->t2 = bankParams.x1;
        pars->massChoice = t02;
        break;
     case Tau0Tau3:
        pars->t0 = bankParams.x0;
        pars->t3 = bankParams.x1;
        pars->massChoice = t03;
        break;
     default:
	if (lalDebugLevel&LALINFO)
	{
		LALPrintError("LALInspiralComputeParams: No choice for parameter space");
	}
  }
  LALInspiralParameterCalc(status->statusPtr, pars);
  pars->fCutoff = 1.L/(pow(6.L,1.5L) * LAL_PI * pars->totalMass * LAL_MTSUN_SI);
  CHECKSTATUSPTR(status);

  DETATCHSTATUSPTR(status);
  RETURN (status);
}
