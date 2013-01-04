/*
*  Copyright (C) 2007 B.S. Sathyaprakash
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
 * \brief Function to find the longest template in a template bank.
 * \author Sathyaprakash, B.S.
 *
 * Given the parameters of a template bank find the longest template
 * in the bank. This is done by looking at the duration for which
 * a signal corresponding to smallest masses lasts. One simply calls
 * the \c LALInspiralWaveLength code for a system consisting
 * of two stars each of mass <tt>mMin.</tt>
 */
void
LALInspiralLongestTemplateInBank
   (
   LALStatus            *status,
   UINT4                *templateLength,
   InspiralCoarseBankIn *coarseIn
   )
{

   InspiralTemplate param;
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);
   ASSERT (coarseIn,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   param.startTime = 0.0;
   param.startPhase = 0.0;
   param.nStartPad = 0;
   param.nEndPad = 0;
   param.ieta = 1;
   param.Theta = 0.;
   param.OmegaS = 0.;
   param.mass1 = coarseIn->mMin;
   param.mass2 = coarseIn->mMin;
   param.fLower = coarseIn->fLower;
   param.fCutoff = coarseIn->fUpper;
   param.tSampling = coarseIn->tSampling;
   param.signalAmplitude = 1.;
   param.order = coarseIn->order;
   param.approximant = coarseIn->approximant;
   param.massChoice = m1Andm2;
   LALInspiralParameterCalc (status->statusPtr, &param);
   CHECKSTATUSPTR(status);
   *templateLength = 0;
   LALInspiralWaveLength (status->statusPtr, templateLength, param);
   CHECKSTATUSPTR(status);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
