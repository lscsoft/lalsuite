/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, Peter Shawhan, B.S. Sathyaprakash, Craig Robinson , Thomas Cokelaer
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
#include <stdio.h>

/**
\ingroup LALInspiralBank_h
\author Churches, D. K. and Sathyaprakash, B.S.
\brief Function which checks whether or not a pair of parameter values are consistent with the search space.

Module which checks whether or not a pair of parameter
values \f$\tau_{0}\f$ and \f$\tau_{2(3)}\f$ correspond to
a user specified range of component masses <tt>(mMin,mMax)</tt> OR to a
minimum value of the component masses \c mMin and maximum total
mass <tt>MMax.</tt> In the first case chirptimes satisfying the
constraint \c mMin\f$\le m_1, m_2 \le\f$\c mMax are accepted
as valid systems. In the second cases chirptimes satisfying the
constraint \c mMin\f$\le m_1, m_2,\f$ and \c MMax\f$\le m=m_1+m_2\f$
are treated as valid.

\heading{Description}

We start with the definition of the chirp times \f$\tau_{0}\f$ and \f$\tau_{3}\f$,
\f{equation}{
\tau_{0} = \frac{5}{256 (\pi f_{a} )^{8/3} m^{5/3} \eta}
\f}
and
\f{equation}{
\tau_{3} = \frac{1}{8 (\pi^{2} f_{a}^{5} )^{1/3} m^{2/3} \eta}
\f}
 These equations may be inverted to yield
\f{equation}{
m = \frac{5}{32 \pi^{2} f_{a}} \frac{\tau_{3}}{\tau_{0}}
\f}
and
\f{equation}{
\eta = \left( \frac{2 \pi^{2}}{25 f_{a}^{3}} \frac{\tau_{0}^{2}}{\tau_{3}^{1/3}}
\right)^{5}\f}

The individual masses may be calculated as follows.  We have
\anchor eq_mass \f{equation}{
m = m_{1} + m_{2}
\tag{eq_mass}
\f}
and
\anchor eq_eta \f{equation}{
\eta = \frac{m_{1} m_{2}}{(m_{1} + m_{2})^{2}}
\tag{eq_eta}
\f}
From Eq.\eqref{eq_mass} we may eliminate either \f$m_{1}\f$ or \f$m_{2}\f$,
\f{equation}{
m_{1} = m - m_{2}
\f}
This may be substituted into Eq.\eqref{eq_eta} to give
\f{equation}{
\eta = \frac{(m - m_{2}) m_{2}}{\left[ (m - m_{2}) + m_{2} \right]^{2}}
= \frac{(m - m_{2}) m_{2}}{m^{2}}
\f}
which may be re--arranged to give
\f{equation}{
m_{2}^{2} - m m_{2} + \eta m^{2} = 0,
\f}
i.e.
\f{equation}{
m_{2} = \frac{ m \pm \sqrt{m^{2}(1 - 4 \eta) }}{2}
\f}
Therefore, since we know that \f$\eta \leq 1/4\f$, real roots are guaranteed.
If we had eliminated \f$m_{2}\f$ rather than \f$m_{1}\f$ then we would have arrived at an identical
expression for
\f$m_{1}\f$, and so of one object has mass
\f{equation}{
m_{1} = \frac{m + \sqrt{m^{2}(1-4 \eta)}}{2}
\f}
then the other object must have mass
\f{equation}{
m_{2} = \frac{m - \sqrt{m^{2}(1-4 \eta)}}{2}
\f}
This function is also given \c mMin and \c MMax as inputs, which it may
use to calculate the minimum value of \f$\eta\f$ which is possible with those inputs,
\f{equation}{
\eta_{min} = \mathtt{ \frac{mMin(MMax - mMin)}{MMax^{2}} }
\f}

To recap, the function calculates \f$m\f$, \f$\eta\f$, \f$\eta_{min}\f$ and \f$m_{1,2}\f$.
It then checks that
\f{equation}{
\eta_{min} \leq \eta \leq 1/4
\f}
and that
\f{equation}{
m_{1} \geq \mathtt{mMin}
\f}
and
\f{equation}{
m_{2} \geq \mathtt{mMin}
\f}
 */
void LALInspiralValidParams(
    LALStatus            *status,	/**< LAL status pointer */
    INT4                 *valid,	/**< [out] 0 means invalid template, 1 means valid */
    InspiralBankParams   bankParams,	/**< [in] Input */
    InspiralCoarseBankIn coarseIn	/**< [in] Input */
    )

{

  InspiralTemplate *Pars=NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.fLower > 0.L, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  *valid = 0;

  if ( bankParams.x0 <=0 || bankParams.x1 <=0 )
  {
    LALInfo( status, "x0 or x1 are less than or equal to zero" );
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }

  Pars = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );
  if ( ! Pars )
  {
    ABORT (status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
  }

  /* First set the chirp times of Pars to be as in bankParams */
  Pars->t0 = bankParams.x0;
  Pars->fLower = coarseIn.fLower;
  switch ( coarseIn.space )
  {
    case Tau0Tau2:
      Pars->t2 = bankParams.x1;
      Pars->massChoice = t02;
      break;
    case Tau0Tau3:
      Pars->t3 = bankParams.x1;
      Pars->massChoice = t03;
      break;
    default:
      ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
  }

  /* Compute all the parameters, including masses, */
  /* corresponding to (t0,t2/t3)                   */
  LALInspiralParameterCalc( status->statusPtr, Pars );
  CHECKSTATUSPTR( status );

  /* If the masses are in the correct range accept as valid parameters */
  switch (coarseIn.massRange)
  {
    case MinComponentMassMaxTotalMass:
      if (
          Pars->mass1 >= coarseIn.mMin &&
          Pars->mass2 >= coarseIn.mMin &&
          Pars->totalMass <= coarseIn.MMax &&
          Pars->eta <= 0.25 &&
          Pars->eta >= coarseIn.etamin
         )
      {
        *valid = 1;
      }
      break;

    case MinMaxComponentMass:
      if (
          Pars->mass1 >= coarseIn.mMin &&
          Pars->mass2 >= coarseIn.mMin &&
          Pars->mass1 <= coarseIn.mMax &&
          Pars->mass2 <= coarseIn.mMax &&
          Pars->eta <= 0.25 &&
          Pars->eta >= coarseIn.etamin
         )
      {
        *valid = 1;
      }
      break;

    case MinMaxComponentTotalMass:
      if (
          Pars->mass1 >= coarseIn.mMin &&
          Pars->mass2 >= coarseIn.mMin &&
          Pars->totalMass <= coarseIn.MMax &&
          Pars->totalMass >= coarseIn.MMin &&
          Pars->eta <= 0.25 &&
          Pars->eta >= coarseIn.etamin
         )
      {
        *valid = 1;
      }
      break;

    default:
      ABORT(status, 999, "Invalid choice for enum InspiralBankMassRange");
  }

  LALFree( Pars );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
