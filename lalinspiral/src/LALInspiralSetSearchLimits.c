/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, Peter Shawhan, Craig Robinson , Thomas Cokelaer
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
#include <lal/LALStdlib.h>

/**
 * \ingroup LALInspiralBank_h
 * \brief Function which calculates the minimum and maximum values of \f$\tau_{0}\f$ and \f$\tau_{3}\f$.
 * \author Churches, D. K.
 *
 * This Function calculates the minimum and maximum values of \f$\tau_{0}\f$ and \f$\tau_{3}\f$
 * as determined by the total mass of the binary \f$m\f$ and the symmetric
 * mass ratio \f$\eta\f$.  The function also calulates the coordinates of the
 * first template in the bank. These coordinates are \f$\tau_{0}=\tau_{0min}\f$,
 * \f$\tau_{3}=\tau_{3min}\f$.
 *
 * ### Description ###
 *
 * We start with the definition of the chirp times \f$\tau_{0}\f$ and \f$\tau_{3}\f$,
 * \f{equation}{
 * \tau_{0} = \frac{5}{256 (\pi f_{a} )^{8/3} m^{5/3} \eta}
 * \f}
 *
 * and
 *
 * \f{equation}{
 * \tau_{3} = \frac{1}{8 (\pi^{2} f_{a}^{5} )^{1/3} m^{2/3} \eta}
 * \f}
 *
 * \f$\tau_{0}\f$ is minimised when \f$\eta=1/4\f$ and \f$\mathtt{m=MMax}\f$.
 * \f$\tau_{0}\f$ is maximised when \f$\eta=1/4\f$ and \f$\mathtt{m=2mMin}\f$.
 * \f$\tau_{3}\f$ is minimised when \f$\eta=1/4\f$ and \f$\mathtt{m=MMax}\f$.
 * \f$\tau_{3}\f$ is maximised when
 * \f$\eta=\mathtt{ mMin(MMax-mMin)/MMax^{2} }\f$.
 */
void
LALInspiralSetSearchLimits (
    LALStatus            *status,	/**< LAL status pointer */
    InspiralBankParams   *bankParams,	/**< [out] containing the boundary of search, current lattice point, etc. */
    InspiralCoarseBankIn  coarseIn	/**< [in] specifies the parameters of the search space */
    )

{
   InspiralTemplate *Pars1=NULL, *Pars2=NULL, *Pars3=NULL, *Pars4=NULL;

   INITSTATUS(status);
   ATTATCHSTATUSPTR( status );

   ASSERT( bankParams, status,
       LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
   ASSERT( coarseIn.space == Tau0Tau2 || coarseIn.space == Tau0Tau3, status,
       LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
   ASSERT( coarseIn.mMin > 0, status,
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
   ASSERT( coarseIn.MMax >= 2. * coarseIn.mMin, status,
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
   ASSERT( coarseIn.mmCoarse > 0., status,
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
   ASSERT( coarseIn.mmCoarse < 1., status,
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
   ASSERT( coarseIn.fLower > 0., status,
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
   ASSERT( coarseIn.tSampling > 0., status,
       LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

   Pars1 = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );
   Pars2 = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );
   Pars3 = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );
   Pars4 = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );

   if ( ! Pars1 || ! Pars2 || ! Pars3 || !Pars4 )
   {
     ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
   }

   /* Initiate three parameter vectors consistent with the coarseIn structure */
   LALInspiralSetParams(status->statusPtr, Pars1, coarseIn);
   CHECKSTATUSPTR(status);
   LALInspiralSetParams(status->statusPtr, Pars2, coarseIn);
   CHECKSTATUSPTR(status);
   LALInspiralSetParams(status->statusPtr, Pars3, coarseIn);
   CHECKSTATUSPTR(status);
   LALInspiralSetParams(status->statusPtr, Pars4, coarseIn);
   CHECKSTATUSPTR(status);

   Pars1->massChoice = Pars2->massChoice = Pars3->massChoice = m1Andm2;
   Pars4->massChoice = m1Andm2;

   /* Calculate the value of the parameters at the three corners */
   /* of the search space                                        */
   Pars1->mass1 = Pars1->mass2 = coarseIn.MMax/2.;
   LALInspiralParameterCalc( status->statusPtr, Pars1 );
   CHECKSTATUSPTR( status );

   if ( coarseIn.massRange == MinMaxComponentTotalMass )
   {
     Pars2->mass1 = Pars2->mass2 = coarseIn.MMin/2.;
   }
   else
   {
     Pars2->mass1 = Pars2->mass2 = coarseIn.mMin;
   }
   LALInspiralParameterCalc( status->statusPtr, Pars2 );
   CHECKSTATUSPTR( status );

   Pars3->mass1 = coarseIn.mMin;
   Pars3->mass2 = coarseIn.MMax - coarseIn.mMin;
   LALInspiralParameterCalc( status->statusPtr, Pars3 );
   CHECKSTATUSPTR( status );

   if ( coarseIn.massRange == MinMaxComponentTotalMass )
   {
     Pars4->mass1 = coarseIn.mMin;
     Pars4->mass2 = coarseIn.MMin - coarseIn.mMin;
     LALInspiralParameterCalc( status->statusPtr, Pars4 );
     CHECKSTATUSPTR( status );
   }
   else
   {
     Pars4->t0 = 0.0;
   }

   /* Find the minimum and maximum values of the parameters and set     */
   /* the search space.  (The minimum values of chirp times are those   */
   /* corresponding to m1 = m2 = MMax/2, i.e., Pars1 structure.         */
   bankParams->x0 = bankParams->x0Min = Pars1->t0;
   bankParams->x0Max = (Pars2->t0 > Pars4->t0) ? Pars2->t0 : Pars4->t0;

   switch ( coarseIn.space )
   {
     case Tau0Tau2:
       bankParams->x1 = bankParams->x1Min = Pars1->t2;
       bankParams->x1Max = (Pars2->t2 > Pars3->t2) ? Pars2->t2 : Pars3->t2;
       break;

     case Tau0Tau3:
       bankParams->x1 = bankParams->x1Min = Pars1->t3;
       bankParams->x1Max = (Pars2->t3 > Pars3->t3) ? Pars2->t3 : Pars3->t3;
       break;

     default:
       ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
   }

   LALFree( Pars1 );
   LALFree( Pars2 );
   LALFree( Pars3 );
   LALFree( Pars4 );

   DETATCHSTATUSPTR( status );
   RETURN( status );
}
