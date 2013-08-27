/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Steven Caudill, Anand Sengupta, Craig Robinson , Thomas Cokelaer, Chris Van Den Broeck
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

/** \defgroup LALInspiralCreateCoarseBank_c Module LALInspiralCreateCoarseBank.c
 * \ingroup LALInspiralBank_h
   \author Churches, D. K and Sathyaprakash, B.S.
   \brief Functions to create a coarse grid of templates.

The coarse grid algorithm works in two stages:
After computing the minimum and maximum chirp-times corresponding to the
search space: \f$(\tau_0^\mathrm{min}, \tau_0^\mathrm{max}),\f$
\f$(\tau_2^\mathrm{min}, \tau_2^\mathrm{max})\f$ (or [Note: In what follows
we will only mention \f$\tau_3\f$; however, the algorithm is itself valid,
and has been implemented, in the case of \f$(\tau_0,\tau_2)\f$ too. However,
we recommend that the space \f$\tau_0\f$-\f$\tau_3\f$ be used.]
\f$(\tau_3^\mathrm{min}, \tau_3^\mathrm{max})\f$) the algorithm
<ol>
<li> chooses a lattice of templates along the equal mass curve and then</li>

<li> lays a rectangular grid in the rectangular region defined by
the minimum and maximum values of the chirp-times and retain only
if (a) the point lies in the parameter space, OR (b) one of the
vertices defined by the rectangle lies in the parameter space.</li>
</ol>

\heading{Templates along the equal mass curve}
The algorithm works in two
stages: In the first stage, templates are built along the equal
mass (that is, \f$\eta=1/4\f$) curve starting from the minimum value
of the Newtonian chirp-time and stopping at its maximum value.
Given the \f$n\f$\ th template at \f$O\f$ with parameters \f$(\tau_0^{(n)},\tau_3^{(n)}),\f$
and given also the distance between templates in our preferred coordinates
\f$(D\tau_0^{(n)},D\tau_3^{(n)}),\f$
consider lines \f$\tau_0 = \tau_0^{(n)} + D\tau_0^{(n)}\f$
(\f$QA\f$ of Fig.\figref{fig_equal_mass_algo}) and
\f$\tau_3 = \tau_3^{(n)} + D\tau_3^{(n)}\f$
(\f$PB\f$ of Fig.\figref{fig_equal_mass_algo}).

\image html  LALInspiralBankHequalmass.png "Fig.[fig_equal_mass_algo]: Algorithm sketching the placement of templates along eta=1/4 curve"
\image latex LALInspiralBankHequalmass.pdf "Algorithm sketching the placement of templates along eta=1/4 curve" width=4.0in

The template next to
\f$(\tau_0^{(n)},\tau_3^{(n)}),\f$ on the equal mass curve, must lie
either along \f$PB\f$ or along \f$QA\f$ (cf. Fig.\figref{fig_equal_mass_algo} in order
that all the signals that may lie on \f$OAB\f$
are spanned by at least one of the two templates.  Clearly, if we were
to place the \f$(n+1)\f$\ th template at \f$B,\f$ some of the signals won't
have the required minimal match. However, placing the \f$(n+1)\f$\ th template at
\f$A\f$ suffices our requirement.
(Note, however, that there is
no guarantee that this will always work; it works only if the curve
along which templates are being laid is a slowly varying function.)
To locate the \f$(n+1)\f$\ th template we
compute the following pairs of coordinates:
\f{eqnarray}{
\tau_0^{(n+1)} = \tau_0^{(n)} + D\tau_0^{(n)}, \ \
\tau_3^{(n+1)} =  4A_3 \left ( \frac{\tau_0^{(n+1)}}{4A_0} \right )^{2/5}
\nonumber \\
\tau_3^{(n+1)} = \tau_3^{(n)} + D\tau_3^{(n)}, \ \
\tau_0^{(n+1)} =  4A_0 \left ( \frac{\tau_3^{(n+1)}}{4A_3} \right )^{5/2},
\f}
where
\f{equation}{
A_0=\frac{5}{256 (\pi f_0)^{8/3}}, \ \ A_3=\frac{\pi}{8 (\pi f_0)^{5/3}}.
\f}
Of the two pairs, the required pair is the one that is closer to the
starting point \f$(\tau_0^{(n)},\tau_3^{(n)}).\f$

\heading{Templates in the rest of the parameter space}
In the second stage, the algorithm begins again at the point
\f$(\tau_0^\mathrm{min}, \tau_3^\mathrm{min}),\f$
corresponding distance between templates
\f$(D\tau_0^\mathrm{min}, D\tau_3^\mathrm{min}),\f$ and chooses a rectangular lattice
of templates in the rectangular region defined by
\f$(\tau_0^\mathrm{min}, \tau_3^\mathrm{min})\f$
\f$(\tau_0^\mathrm{max}, \tau_3^\mathrm{min})\f$
\f$(\tau_0^\mathrm{max}, \tau_3^\mathrm{max})\f$  and
\f$(\tau_0^\mathrm{min}, \tau_3^\mathrm{max})\f$.
The implementation of the algorithm along the equal mass curve and
in a rectangular lattice in the rest of the parameter space is shown
plotted in Fig.\figref{fig_coarse}, where the templates
chosen are represented as points.

\image html  LALInspiralBankHCoarse2.png "Fig.[fig_coarse]: Algorithm sketching the construction of a rectangular lattice of templates"
\image latex LALInspiralBankHCoarse2.pdf "Algorithm sketching the construction of a rectangular lattice of templates" width=4.5in

\heading{Algorithm}

The algorithm to lay templates along the equal-mass curve is as follows:
<tt>
<ul>
<li> Begin at \f$\tau_0 = \tau_0^\mathrm{min}\f$
<li> do while \f$(\tau_0 < \tau_0^\mathrm{max})\f$<br>
   {
   <ul>
   <li> \f$\tau_0^A = \tau_0 + D\tau_0, \ \ \tau_3^A = 4A_3 \left ( {\tau_0^A}/{4A_0} \right )^{2/5}\f$
   <li> \f$\tau_3^B = \tau_3 + D\tau_3, \ \ \tau_0^B = 4A_0 \left ( {\tau_3^B}/{4A_3} \right )^{5/2}\f$
   <li> if (\f$(\tau_0^A,\tau_3^A)\f$ is closer \f$(\tau_0,\tau_3)\f$ than \f$(\tau_0^B,\tau_3^B)\f$)<br>
        {<br>
           \f$\tau_0 = \tau_0^A, \tau_3 = \tau_3^A\f$<br>
        }<br>
   <li> else <br>
        {<br>
        \f$\tau_0 = \tau_0^B, \tau_3 = \tau_3^B\f$ <br>
        }<br>
   <li> Add \f$(\tau_0, \tau_3)\f$ to InspiralTemplateList
   <li> numTemplates++
   <li> Compute metric at \f$(\tau_0, \tau_3)\f$
   <li> Compute distance between templates at  new point: \f$(D\tau_0, D\tau_3)\f$
   </ul>
}
</ul>
</tt>

The algorithm to lay templates in the rest of the parameter space
is as follows:
<tt>
<ul>
<li> Begin at \f$\tau_0 = \tau_0^\mathrm{min}, \tau_3 = \tau_3^\mathrm{min}\f$
<li> Compute metric at \f$(\tau_0, \tau_3)\f$
<li> Compute distance between templates at  new point: \f$(D\tau_0, D\tau_3)\f$
<li> Add \f$(\tau_0, \tau_3)\f$ to InspiralTemplateList
<li> numTemplates++
<li> do while (\f$\tau_3 <= \tau_3^\mathrm{max}\f$)<br>
   {<br>
   <ul>
   <li> do while (\f$\tau_0 <= \tau_0^\mathrm{max}\f$)<br>
       {<br>
       <ul>
       <li> if (\f$(\tau_0, \tau_3)\f$ is inside the parameter space)<br>
            {<br>
            <ul>
            <li> Compute metric at (\f$\tau_0, \tau_3\f$)
            <li> Compute distance between templates at  new point: (\f$D\tau_0, D\tau_3\f$)
            <li> Add (\f$\tau_0, \tau_3\f$) to InspiralTemplateList
            <li> numTemplates++
            </ul>
            }<br>
       <li> Increment \f$\tau_0:\f$ \f$\tau_0 = \tau_0 + D\tau_0\f$<br>
       </ul>
       }<br>
   <li> Increment \f$\tau_3:\f$ \f$\tau_3 = \tau_3 + D\tau_3\f$
   <li> Get next template along \f$\tau_3=\mathrm{const.}\f$: \f$(\tau_0, \tau_3)\f$<br>
   </ul>
   }<br>
</ul>
</tt>
*/
/*@{*/


#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdio.h>
#include <lal/FindRoot.h>


/** \see See \ref LALInspiralCreateCoarseBank_c for documentation */
void
LALInspiralCreateCoarseBank(
    LALStatus            *status,	/**< LAL-status pointer */
    InspiralTemplateList **list,	/**< [out] an array containing the template bank parameters */
    INT4                 *nlist,	/**< [out] the number of templates found by the function */
    InspiralCoarseBankIn coarseIn	/**< [in] specifies the search space, range of masses, etc */
    )
{
  INT4 i;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.shf.data, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( coarseIn.shf.data->data, status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL );
  ASSERT( coarseIn.mmCoarse > 0.L, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.mmCoarse < 1.L, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.fLower > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.tSampling > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.tSampling >= 2.*coarseIn.fUpper, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  switch( coarseIn.approximant )
  {
    case BCV:
      ASSERT( coarseIn.space == Psi0Psi3, status,
          LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
      LALInspiralCreateBCVBank( status->statusPtr, list, nlist, coarseIn );
      CHECKSTATUSPTR( status );
      break;

    case AmpCorPPN:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case TaylorF1:
    case TaylorF2:
    case Eccentricity:
    case PadeT1:
    case PadeF1:
    case EOB:
    case EOBNR:
    case EOBNRv2:
    case IMRPhenomA:
    case IMRPhenomB:
    case TaylorEt:
    case TaylorN:
    case FindChirpPTF:
      ASSERT( coarseIn.space == Tau0Tau2 || coarseIn.space == Tau0Tau3, status,
          LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );

      /* Thomas:: we can use either a square placement or an hexagonal
       * placement. The hexagonal placement is along the eigenvalues but
       * not the square one.*/
      if (coarseIn.gridSpacing == Hexagonal){
	LALInspiralCreatePNCoarseBankHexa( status->statusPtr, list, nlist, coarseIn );
	CHECKSTATUSPTR( status );
      }
      else if (coarseIn.gridSpacing == HybridHexagonal){
	LALInspiralCreatePNCoarseBankHybridHexa( status->statusPtr, list, nlist, coarseIn );
	CHECKSTATUSPTR( status );
      }
      else if (coarseIn.gridSpacing == SquareNotOriented){
	LALInspiralCreatePNCoarseBank( status->statusPtr, list, nlist, coarseIn );
	CHECKSTATUSPTR( status );
      }
      else {
        ABORT( status, LALINSPIRALBANKH_EGRIDSPACING, LALINSPIRALBANKH_MSGEGRIDSPACING );
      }

      /* Anand:: Nudge the templates only if using max-total-mass cut */
      if ( coarseIn.massRange == MinComponentMassMaxTotalMass  ||
             coarseIn.massRange == MinMaxComponentTotalMass )
      {
          LALNudgeTemplatesToConstantTotalMassLine( status->statusPtr, list, (*nlist), coarseIn);
          CHECKSTATUSPTR( status );
      }

      break;

    default:
      ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
      break;
  }

  /* record the minimal match of the bank in the template and   */
  /* set up the tmplts as a linked list so that it can be       */
  /* manipulated easily by findchirp                            */
  for ( i = 0; i < *nlist - 1 ; ++i )
  {
    (*list)[i].params.minMatch = (REAL4) coarseIn.mmCoarse;
    (*list)[i].params.next     = &((*list)[i+1].params);
    (*list)[i].params.fine     = NULL;
  }
  (*list)[i].params.minMatch = (REAL4) coarseIn.mmCoarse;
  (*list)[i].params.next = NULL;
  (*list)[i].params.fine = NULL;

  DETATCHSTATUSPTR(status);
  RETURN (status);
}

/** Anand: 26 October 2006
 * This function nudges the templates in the list to
 * the (max-total-mass = constant) line.
 * This is done only for those templates whose total
 * mass exceeds the desired max-total-mass value. The
 * templates are nudged along the metric eigen direction
 * until they lie on the said line.
 */
void
LALNudgeTemplatesToConstantTotalMassLine(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 nlist,
    InspiralCoarseBankIn coarseIn
    )
{
  InspiralTemplate      *tempPars=NULL;
  InspiralMetric        *metric=NULL;
  InspiralMomentsEtc    moments;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* If there are no templates, return now */
  if ( nlist <= 0 )
  {
      LALWarning( status, "number of templates is <= 0 ! " );

      DETATCHSTATUSPTR(status);
      RETURN (status);
  }

  /* Allocate memory (only required to calculate noise moments) */
  tempPars = (InspiralTemplate *)
          LALCalloc( 1, sizeof(InspiralTemplate) );
  metric   = (InspiralMetric *)
          LALCalloc( 1, sizeof(InspiralMetric) );

  /* Init the tempPars */
  LALInspiralSetParams( status->statusPtr, tempPars, coarseIn );
  CHECKSTATUSPTR( status );

  tempPars->totalMass  = coarseIn.MMax;
  tempPars->eta        = 0.25;
  tempPars->ieta       = 1.L;
  tempPars->fLower     = coarseIn.fLower;
  tempPars->massChoice = totalMassAndEta;
  LALInspiralParameterCalc( status->statusPtr, tempPars );
  CHECKSTATUSPTR( status );

  /* Get the moments of the PSD required in the computation of the metric */
  LALGetInspiralMoments( status->statusPtr, &moments, &coarseIn.shf, tempPars );
  CHECKSTATUSPTR( status );

  /* Loop over template list and nudge the templates if required */
  {
      INT4   i;
      REAL4  P, Q, M, C, ms; /*, t0, t3;*/

      M = coarseIn.MMax*LAL_MTSUN_SI;
      P = (5./256.)*pow( (LAL_PI*coarseIn.fLower), -8./3. ) ;
      Q = (LAL_PI/8.)*pow( (LAL_PI*coarseIn.fLower), -5./3. ) ;

      for (i=0; i < nlist; i++)
      {
          /* If the totalMass of this template exceeds max-total-mass
           * then nudge along the metric eigen-direction.
           */
          if ( (*list)[i].params.totalMass > coarseIn.MMax )
          {
             ms = tan( LAL_PI/2. + (*list)[i].metric.theta );
             C  = (*list)[i].params.t3 - ms*((*list)[i].params.t0);

             /* Calculate the new co-ordinates in tau0-tau3 space */
             (*list)[i].params.t3  = C / ( 1. -  (ms*P/(M*Q)) );
             (*list)[i].params.t0  = P*(*list)[i].params.t3/(M*Q);

             /* Calculate the other parameters */
             LALInspiralParameterCalc( status->statusPtr, &(*list)[i].params );
             CHECKSTATUSPTR( status );

             /* Check that the new point has not gone down below the
              * equal mass line. If it has, set it to m1=m2=coarseIn.MMax/2.0
              */
             if ( (*list)[i].params.eta > 0.25L )
             {
                 InputMasses originalMassChoice = (*list)[i].params.massChoice;

                 (*list)[i].params.totalMass = coarseIn.MMax ;
                 (*list)[i].params.eta       = 0.25L;
                 (*list)[i].params.massChoice = totalMassAndEta;

                 LALInspiralParameterCalc( status->statusPtr, &(*list)[i].params );
                 CHECKSTATUSPTR( status );

                 /* Reset the massChoice to whatever it was */
                 (*list)[i].params.massChoice = originalMassChoice;
             }

             /* Recalculate the metric at this new point */
             LALInspiralComputeMetric( status->statusPtr, &((*list)[i].metric),
                     &((*list)[i].params), &moments );
             CHECKSTATUSPTR( status );
          }
      }/* Loop over templates */
  }

  /* Clean up */
  LALFree( tempPars );
  LALFree( metric );

  /* Normal exit */
  DETATCHSTATUSPTR(status);
  RETURN (status);
}

/** \see See \ref LALInspiralCreateCoarseBank_c for documentation */
void
LALInspiralCreatePNCoarseBank(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    )
{
  InspiralBankParams bankPars, bankParsOld;
  InspiralTemplate *tempPars;
  InspiralMetric metric;
  InspiralMomentsEtc moments;
  INT4 validPars;
  REAL8 x01, x02, x11, x12, dist1, dist2, ndx1, ndx2, a25;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.mMin > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.mMax > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.MMax >= 2.*coarseIn.mMin, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  ndx1 = 0.0;
  ndx2 = 0.0;
  a25 = 0.0;

  /* Number of templates is nlist */
  *nlist = 0;

  /* Set the elements of the metric and tempPars structures in  */
  /* conformity with the coarseIn structure                     */
  if ( !
      (tempPars = (InspiralTemplate *)LALCalloc( 1, sizeof(InspiralTemplate) ))
      )
  {
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }
  metric.space = coarseIn.space;
  LALInspiralSetParams( status->statusPtr, tempPars, coarseIn );
  CHECKSTATUSPTR( status );

  /* Identify the boundary of search and parameters for the     */
  /* first lattice point                                        */
  LALInspiralSetSearchLimits( status->statusPtr, &bankPars, coarseIn );
  CHECKSTATUSPTR( status );
  tempPars->totalMass = coarseIn.MMax;
  tempPars->eta = 0.25;
  tempPars->ieta = 1.L;
  tempPars->fLower = coarseIn.fLower;
  tempPars->massChoice = totalMassAndEta;
  LALInspiralParameterCalc( status->statusPtr, tempPars );
  CHECKSTATUSPTR( status );

  /* Get the moments of the PSD integrand and other parameters */
  /* required in the computation of the metric                 */
  LALGetInspiralMoments( status->statusPtr, &moments, &coarseIn.shf, tempPars );
  CHECKSTATUSPTR( status );

  /* compute the metric at this point, update bankPars and add */
  /* the params to the list                                    */
  LALInspiralComputeMetric( status->statusPtr, &metric, tempPars, &moments );
  CHECKSTATUSPTR( status );
  LALInspiralUpdateParams( status->statusPtr, &bankPars, metric,
      coarseIn.mmCoarse );
  CHECKSTATUSPTR( status );

  /* add the first template to the template list */
  *list = (InspiralTemplateList*)
    LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist + 1) );
  if ( ! *list )
  {
    LALFree( tempPars );
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }
  memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );

  (*list)[*nlist].ID = *nlist;
  (*list)[*nlist].params = *tempPars;
  (*list)[*nlist].metric = metric;
  ++(*nlist);

  /* First lay templates along the equal mass curve; i.e. eta=1/4.      */
  /* Choose the constant and the index converting the chirp times to    */
  /* one another along the curve depending on whether the templates     */
  /* are laid along the tau0-tau2 or tau0-tau3 space                    */
  switch ( coarseIn.space )
  {
    case Tau0Tau2:
      ndx1 = 0.6L;
      ndx2 = 1.L/ndx1;
      a25 = pow(64.L/5.L, ndx1)*(2435.L/8064.L)/pow(LAL_PI*coarseIn.fLower,.4L);
      break;

    case Tau0Tau3:
      a25 = LAL_PI_2 * pow(64.L/5.L, .4L)/pow(LAL_PI * coarseIn.fLower, .6L);
      ndx1 = 0.4L;
      ndx2 = 2.5L;
      break;

    case Psi0Psi3:
    case PTFIntrinsic:
    case PTFFull:
      ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
      break;
  }

  bankParsOld = bankPars;

  while ( bankPars.x0 < bankPars.x0Max )
  {
    x01 = bankPars.x0 + bankPars.dx0;
    x11 = a25 * pow(x01,ndx1);
    x12 = bankPars.x1 + bankPars.dx1;
    x02 = pow(x12/a25,ndx2);
    dist1 = pow(bankPars.x0 - x01,2.L) + pow(bankPars.x1 - x11, 2.L);
    dist2 = pow(bankPars.x0 - x02,2.L) + pow(bankPars.x1 - x12, 2.L);
    if ( dist1 < dist2 )
    {
      bankPars.x0 = x01;
      bankPars.x1 = x11;
    }
    else
    {
      bankPars.x0 = x02;
      bankPars.x1 = x12;
    }

    /* If this is a valid point add it to our list */
    LALInspiralValidTemplate( status->statusPtr,
        &validPars, bankPars, coarseIn );
    CHECKSTATUSPTR( status );

    if ( validPars )
    {
      LALInspiralComputeParams( status->statusPtr,
          tempPars, bankPars, coarseIn);
      CHECKSTATUSPTR( status );
      LALInspiralComputeMetric( status->statusPtr,
          &metric, tempPars, &moments );
      CHECKSTATUSPTR( status );
      LALInspiralUpdateParams( status->statusPtr,
          &bankPars, metric, coarseIn.mmCoarse );
      CHECKSTATUSPTR( status );

      *list = (InspiralTemplateList *)
        LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist + 1) );
      if ( ! *list )
      {
        LALFree( tempPars );
        ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
      memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );

      (*list)[*nlist].ID = *nlist;
      (*list)[*nlist].params = *tempPars;
      (*list)[*nlist].metric = metric;
      ++(*nlist);
    }
  }

  /* Begin with the parameters found at the first lattice point */
  bankPars = bankParsOld;

  /* Loop along x1 and x0 coordinates until maximum values are reached */
  while ( bankPars.x1 <= bankPars.x1Max )
  {
    /* step along the tau0 axis until the boundary is reached */
    while ( bankPars.x0 <= bankPars.x0Max )
    {
      /* If this is a valid point add it to our list */
      LALInspiralValidTemplate( status->statusPtr,
          &validPars, bankPars, coarseIn );
      CHECKSTATUSPTR( status );

      if ( validPars )
      {
        LALInspiralComputeParams( status->statusPtr,
            tempPars, bankPars, coarseIn );
        CHECKSTATUSPTR( status );
        LALInspiralComputeMetric( status->statusPtr,
            &metric, tempPars, &moments );
        CHECKSTATUSPTR( status );
        LALInspiralUpdateParams( status->statusPtr,
            &bankPars, metric, coarseIn.mmCoarse );
        CHECKSTATUSPTR( status );

        *list = (InspiralTemplateList *)
          LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist + 1) );
        if ( ! *list )
        {
          LALFree( tempPars );
          ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
        }
        memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );

        (*list)[*nlist].ID = *nlist;
        (*list)[*nlist].params = *tempPars;
        (*list)[*nlist].metric = metric;
        ++(*nlist);
      }

      bankPars.x0 += bankPars.dx0;
    }
    bankPars = bankParsOld;
    bankPars.x1 += bankPars.dx1;

    /* Find the t0 coordinate of the next template close to the t2/t3 axis */
    LALInspiralNextTemplate( status->statusPtr, &bankPars, metric );
    CHECKSTATUSPTR( status );

    /* Hop along t0-axis until t0 is inside the region of interest or quit */
    LALInspiralValidTemplate( status->statusPtr,
        &validPars, bankPars, coarseIn );
    CHECKSTATUSPTR( status );
    while ( validPars == 0 && bankPars.x0 < bankPars.x0Max )
    {
      bankPars.x0 += bankPars.dx0;
      LALInspiralValidTemplate( status->statusPtr,
          &validPars, bankPars, coarseIn );
      CHECKSTATUSPTR( status );
    }
    bankParsOld = bankPars;
  }
  LALFree( tempPars );

  DETATCHSTATUSPTR( status );
  RETURN ( status );
}
/*@}*/
