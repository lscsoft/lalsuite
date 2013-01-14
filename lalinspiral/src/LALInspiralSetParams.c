/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, Thomas Cokelaer
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

#include <stdio.h>
#include <lal/LALInspiralBank.h>

/** \ingroup LALInspiralBank_h
 * \brief A routine that fills an \c ::InspiralTemplate structure based on the values in the \c ::InspiralCoarseBankIn structure.
 * \author Churches, D. K and Sathyaprakash, B.S.
 *
 * This function takes as an input a structure of type
 * \c ::InspiralCoarseBankIn and it fills up the
 * elements of a structure of type \c ::InspiralTemplate.
 * The function sets the fields
 * \c massChoice, \c ieta, \c signalAmplitude,
 * \c tSampling, \c fLower, \c fCutoff,
 * \c order,
 * \c approximant, \c nStartPad, \c nEndPad.
 */
void LALInspiralSetParams(LALStatus            *status,		/**< LAL status pointer */
                          InspiralTemplate     *tempPars,	/**< [out] Output */
                          InspiralCoarseBankIn coarseIn		/**< [in] Input */
                          )
{

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);
   ASSERT (tempPars,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT ((INT4)coarseIn.space >= 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT ((INT4)coarseIn.space <= 1, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

   switch (coarseIn.space) {
      case Tau0Tau2:
         tempPars->massChoice = t02;
         break;
      case Tau0Tau3:
         tempPars->massChoice = t03;
         break;
      default: /* JC: DEFAULT CASE ADDED HERE */
         ABORT( status, 9999, "Default case in switch." );
   }

   tempPars->ieta = 1;
   tempPars->signalAmplitude = 1.;
   tempPars->tSampling = coarseIn.tSampling;
   tempPars->fLower = coarseIn.fLower;
   tempPars->fCutoff = coarseIn.fUpper;
   tempPars->order = coarseIn.order;
   tempPars->approximant = coarseIn.approximant;
   tempPars->nStartPad = 0;
   tempPars->nEndPad = 0;
   tempPars->Theta = 0.;
   tempPars->OmegaS = 0.;
   tempPars->startTime = 0.;
   tempPars->startPhase = 0.;
   tempPars->spin1[0] = tempPars->spin1[1] = tempPars->spin1[2] = 0.;
   tempPars->spin2[0] = tempPars->spin2[1] = tempPars->spin2[2] = 0.;
   tempPars->inclination = 0.;
   tempPars->eccentricity = 0.;

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
