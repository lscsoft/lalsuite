/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
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

/**
 * \ingroup LALInspiralBank_h
 * \author Sathyaprakash, B.S. and Churches, D. K.
 * \brief Function to create a fine grid of templates.
 *
 * The fine grid algorithm is a very simple algorithm that computes a uniform
 * grid of templates around a given coordinate point -- which can in particular be
 * a coarse grid point -- from a knowledge of the metric at the coordinate point
 * and the coarse and fine grid minimal matches, \f$D\tau_{0,3}\f$ and
 * \f$d\tau_{0,3},\f$ respectively. Since \f$D\tau\f$ is not necessarily an
 * integral multiple of \f$d\tau\f$ the rectangular fine grid about the point
 * in question will be larger than required. The algorithm chooses templates
 * \e symmetrically about the given coarse grid point. It does so
 * by laying a rectangular lattice of templates with spacings
 * \f$d\tau_0\f$ and \f$d\tau_3,\f$ in the rectangular region defined by
 * \f$(\tau_0 - \Delta \tau_0, \tau_3 - \Delta \tau_3),\f$
 * \f$(\tau_0 + \Delta \tau_0, \tau_3 - \Delta \tau_3),\f$
 * \f$(\tau_0 + \Delta \tau_0, \tau_3 + \Delta \tau_3)\f$ and
 * \f$(\tau_0 - \Delta \tau_0, \tau_3 + \Delta \tau_3),\f$
 * where
 * \f[\Delta\tau_0 = d\tau_0 \left [ \frac{D\tau_0}{2d\tau_0} \right ], \f]
 * and for any \f$x\f$, \f$[x]\f$ denotes the smallest integer greater than or
 * equal to \f$x\f$.
 * \image html  LALInspiralBankHfine.png "Fig.[fig_fine]: Algorithm sketching the construction of a rectangular fine grid around a given coordinate point"
 * \image latex LALInspiralBankHfine.pdf "Algorithm sketching the construction of a rectangular fine grid around a given coordinate point" width=4.5in
 * The algorithm takes as input a structure of type
 * \c InspiralFineBankIn and returns a <tt>pointer-to-a-pointer</tt>
 * of type \c InspiralTemplateList as well as the number of fine grid
 * templates \c int around the lattice point in question.
 *
 * The spacing between fine grid templates is chosen
 * to be a constant determined by the metric at the coarse grid point; for
 * example,
 * \f[d\tau_0 = \sqrt{\frac{2 (1 - MM_\textrm{Fine})}{g_{00}} }.\f]
 * Only those grid points that are within the parameter space boundary, or
 * have the vertices of the ambiguity rectangle inside the parameter
 * space, are kept and others are discarded.
 *
 * ### Algorithm ###
 *
 * The Fine grid algorithm works as follows:
 * <tt>
 * <ul>
 * <li> From input structure extract coordinates of the grid point \f$(\tau_0^G, \tau_3^G).\f$
 * <li> Compute coarse and fine grid spacings \f$(D\tau_0, D\tau_3)\f$ and \f$(d\tau_0, d\tau_3)\f$
 * <li>Compute half-sides of the <i>smallest</i> symmetric rectangle about \f$(\tau_0, \tau_3)\f$:
 * <ul>
 * <li> \f$\Delta\tau_0 =  d\tau_0 \,\mathrm{ceil}[D\tau_0/(2d\tau_0)],\f$ \f$\Delta\tau_3 =  d\tau_3 \,\mathrm{ceil}[D\tau_3/(2d\tau_3)],\f$
 * </ul>
 * <li> Begin at \f$\tau_3 = \tau_3^G - \Delta \tau_3,\f$
 * <li> do while (\f$\tau_3 <= \tau_3^G+\Delta \tau_3\f$)<br>
 * {<br>
 * <ul>
 * <li> Begin at \f$\tau_0 = \tau_0^G - \Delta \tau_0,\f$
 * <li> do while (\f$\tau_0 <= \tau_0^G+\Delta \tau_0\f$)<br>
 * {<br>
 * <ul>
 * <li> if (\f$(\tau_0, \tau_3)\f$ is inside the parameter space)<br>
 * {<br>
 * <ul>
 * <li> Add (\f$\tau_0, \tau_3\f$) to InspiralTemplateList
 * <li> numTemplates++
 * </ul>
 * }<br>
 * <li> Increment \f$\tau_0:\f$ \f$\tau_0 = \tau_0 + d\tau_0\f$
 * </ul>
 * }<br>
 * <li>Increment \f$\tau_3:\f$ \f$\tau_3 = \tau_3 + d\tau_3\f$
 * </ul>
 * }
 * </ul>
 * </tt>
 */
void LALInspiralCreateFineBank(LALStatus            *status,	/**< LAL status pointer */
                               InspiralTemplateList **outlist,	/**< [out] containing an array of template bank parameters */
                               INT4                 *nlist,	/**< [out] the number of fine bank templates around a given coarse-mesh point */
                               InspiralFineBankIn   fineIn	/**< [in] the parameters required to find the fine bank */
                               )
{

  REAL8 x0, x1, Dx0, Dx1, dx0, dx1, x0FineMin, x1FineMin;
  INT4  i, j, validPars, bins0, bins1;
  InspiralTemplate   *tempPars=NULL;
  InspiralBankParams *bankPars=NULL;


  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);
  ASSERT ((INT4)fineIn.coarseIn.space>=0,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT ((INT4)fineIn.coarseIn.space<=1,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT ((REAL8)fineIn.templateList.params.t0 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

  /* set the number of fine templates generated to zero      */
  /* otherwise the LALCalloc will fail horribly, since there */
  /* is nothing that guarantees that *nlist is a reasonable  */
  /* number when the function is called                      */
  *nlist = 0;

  tempPars = (InspiralTemplate *) LALCalloc(1, sizeof(InspiralTemplate));
  bankPars = (InspiralBankParams *) LALCalloc(1, sizeof(InspiralBankParams));
  *tempPars = fineIn.templateList.params;
  switch (fineIn.coarseIn.space) {
    case Tau0Tau2:
      ASSERT (fineIn.templateList.params.t2 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
      bankPars->x0 = fineIn.templateList.params.t0;
      bankPars->x1 = fineIn.templateList.params.t2;
      break;
    case Tau0Tau3:
      ASSERT (fineIn.templateList.params.t3 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
      bankPars->x0 = fineIn.templateList.params.t0;
      bankPars->x1 = fineIn.templateList.params.t3;
      break;
    default: /* JC: DEFAULT CASE ADDED HERE */
      ABORT( status, 9999, "Default case in switch." );
  }

  LALInspiralUpdateParams(status->statusPtr,bankPars,fineIn.templateList.metric,fineIn.coarseIn.mmCoarse);
  CHECKSTATUSPTR(status);
  x0 = bankPars->x0;
  x1 = bankPars->x1;
  Dx0 = bankPars->dx0;
  Dx1 = bankPars->dx1;

  LALInspiralUpdateParams(status->statusPtr,bankPars,fineIn.templateList.metric,fineIn.coarseIn.mmFine);
  CHECKSTATUSPTR(status);
  dx0 = bankPars->dx0;
  dx1 = bankPars->dx1;

  bins0 = (int)(Dx0/dx0) + 1;
  bins1 = (int)(Dx1/dx1) + 1;

  x0FineMin = x0 - (float) bins0/2. * dx0;
  x1FineMin = x1 - (float) bins1/2. * dx1;

  bankPars->x1 = x1FineMin;
  for(i=0; i<=bins1; i++) {
     bankPars->x0 = x0FineMin;
     for(j=0; j<=bins0; j++) {
       LALInspiralValidTemplate(status->statusPtr, &validPars, *bankPars, fineIn.coarseIn);
       CHECKSTATUSPTR(status);
       if (validPars) {
         LALInspiralComputeParams(status->statusPtr, tempPars, *bankPars, fineIn.coarseIn);
         CHECKSTATUSPTR(status);
/*
    On failure realloc() returns a NULL to outlist, hence there is
    no need to explicitly free the outlist
*/
         if (!(*outlist = (InspiralTemplateList*)
            LALRealloc(*outlist, sizeof(InspiralTemplateList)*(*nlist+1)))) {
            ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
	    outlist = NULL;
         }
         memset( *outlist + *nlist, 0, sizeof(InspiralTemplateList) );
         (*outlist)[*nlist].params = *tempPars;
         ++(*nlist);
       }
       bankPars->x0+=bankPars->dx0;
     }
     bankPars->x1+=bankPars->dx1;
  }
  if (tempPars!=NULL) LALFree(tempPars);
  if (bankPars!=NULL) LALFree(bankPars);
  DETATCHSTATUSPTR(status);
  RETURN (status);
}
