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

/*  <lalVerbatim file="LALInspiralCreateCoarseBankCV">
Author: Churches, D. K and Sathyaprakash, B.S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralCreateCoarseBank.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralCreateCoarseBankCP}
\idx{LALInspiralCreateCoarseBank()}

\begin{itemize}
   \item \texttt{list,} Output, an array containing the template bank parameters
   \item \texttt{nlist,} Output, the number of templates found by the function
   \item \texttt{coarseIn,} Input, specifies the search space, range of masses, etc.
\end{itemize}

The coarse grid algorithm works in two stages:
After computing the minimum and maximum chirp-times corresponding to the
search space: $(\tau_0^{\rm min}, \tau_0^{\rm max}),$
$(\tau_2^{\rm min}, \tau_2^{\rm max})$ (or\footnote{In what follows
we will only mention $\tau_3$; however, the algorithm is itself valid,
and has been implemented, in the case of $(\tau_0,\tau_2)$ too. However,
we recommend that the space $\tau_0$-$\tau_3$ be used.}
$(\tau_3^{\rm min}, \tau_3^{\rm max})$) the algorithm

\begin{enumerate}
\item chooses a lattice of templates along the equal mass curve and then

\item lays a rectangular grid in the rectangular region defined by
the minimum and maximum values of the chirp-times and retain only
if (a) the point lies in the parameter space, OR (b) one of the
vertices defined by the rectangle lies in the parameter space.
\end{enumerate}

\subsubsection*{Description}
\paragraph*{Templates along the equal mass curve}
The algorithm works in two
stages: In the first stage, templates are built along the equal
mass (that is, $\eta=1/4$) curve starting from the minimum value
of the Newtonian chirp-time and stopping at its maximum value.
Given the $n$~th template at $O$ with parameters $(\tau_0^{(n)},\tau_3^{(n)}),$
and given also the distance between templates in our preferred coordinates
$(D\tau_0^{(n)},D\tau_3^{(n)}),$
consider lines $\tau_0 = \tau_0^{(n)} + D\tau_0^{(n)}$
($QA$ of Fig.~\ref{fig:equal mass algo}) and
$\tau_3 = \tau_3^{(n)} + D\tau_3^{(n)}$
($PB$ of Fig.~\ref{fig:equal mass algo}).
\begin{figure}[h]
\begin{center}
\includegraphics[angle=-90,width=4.0 true in]{LALInspiralBankHequalmass}
\end{center}
\caption{Algorithm sketching the placement of templates along $\eta=1/4$ curve.}
\label{fig:equal mass algo}
\end{figure}
The template next to
$(\tau_0^{(n)},\tau_3^{(n)}),$ on the equal mass curve, must lie
either along $PB$ or along $QA$ (cf. Fig.~\ref{fig:equal mass algo}) in order
that all the signals that may lie on $OAB$
are spanned by at least one of the two templates.  Clearly, if we were
to place the $(n+1)$~th template at $B,$ some of the signals won't
have the required minimal match. However, placing the $(n+1)$~th template at
$A$ suffices our requirement.
(Note, however, that there is
no guarantee that this will always work; it works only if the curve
along which templates are being laid is a slowly varying function.)
To locate the $(n+1)$~th template we
compute the following pairs of coordinates:
\begin{eqnarray}
\tau_0^{(n+1)} = \tau_0^{(n)} + D\tau_0^{(n)}, \ \
\tau_3^{(n+1)} =  4A_3 \left ( \frac{\tau_0^{(n+1)}}{4A_0} \right )^{2/5}
\nonumber \\
\tau_3^{(n+1)} = \tau_3^{(n)} + D\tau_3^{(n)}, \ \
\tau_0^{(n+1)} =  4A_0 \left ( \frac{\tau_3^{(n+1)}}{4A_3} \right )^{5/2},
\end{eqnarray}
where
\begin{equation}
A_0=\frac{5}{256 (\pi f_0)^{8/3}}, \ \ A_3=\frac{\pi}{8 (\pi f_0)^{5/3}}.
\end{equation}
Of the two pairs, the required pair is the one that is closer to the
starting point $(\tau_0^{(n)},\tau_3^{(n)}).$

\paragraph*{Templates in the rest of the parameter space}
In the second stage, the algorithm begins again at the point
$(\tau_0^{\rm min}, \tau_3^{\rm min}),$
corresponding distance between templates
$(D\tau_0^{\rm min}, D\tau_3^{\rm min}),$ and chooses a rectangular lattice
of templates in the rectangular region defined by
$(\tau_0^{\rm min}, \tau_3^{\rm min})$
$(\tau_0^{\rm max}, \tau_3^{\rm min})$
$(\tau_0^{\rm max}, \tau_3^{\rm max})$  and
$(\tau_0^{\rm min}, \tau_3^{\rm max})$.
The implementation of the algorithm along the equal mass curve and
in a rectangular lattice in the rest of the parameter space is shown
plotted in Fig.~\ref{fig:coarse}, where the templates
chosen are represented as points.
\begin{figure}[h]
\begin{center}
\includegraphics[angle=-90,width=4.5 true in]{LALInspiralBankHCoarse2}
\caption{Algorithm sketching the construction of a rectangular lattice of templates.}
\label{fig:coarse}
\end{center}
\end{figure}


\subsubsection*{Algorithm}

The algorithm to lay templates along the equal-mass curve is as follows:
\begin{obeylines}
\texttt{
\hskip 1 true cm Begin at $\tau_0 = \tau_0^{\rm min}$
\hskip 1 true cm do while $(\tau_0 < \tau_0^{\rm max})$
\hskip 1 true cm \{
\hskip 2 true cm $\tau_0^A = \tau_0 + D\tau_0, \ \ \tau_3^A = 4A_3 \left ( {\tau_0^A}/{4A_0} \right )^{2/5}$
\hskip 2 true cm $\tau_3^B = \tau_3 + D\tau_3, \ \ \tau_0^B = 4A_0 \left ( {\tau_3^B}/{4A_3} \right )^{5/2}$
\hskip 2 true cm if ($(\tau_0^A,\tau_3^A)$ is closer $(\tau_0,\tau_3)$ than $(\tau_0^B,\tau_3^B)$)
\hskip 2 true cm \{
\hskip 3 true cm $\tau_0 = \tau_0^A, \tau_3 = \tau_3^A$
\hskip 2 true cm \}
\hskip 2 true cm else
\hskip 2 true cm \{
\hskip 3 true cm $\tau_0 = \tau_0^B, \tau_3 = \tau_3^B$
\hskip 2 true cm \}
\hskip 2 true cm Add $(\tau_0, \tau_3)$ to InspiralTemplateList
\hskip 2 true cm numTemplates++
\hskip 2 true cm Compute metric at $(\tau_0, \tau_3)$
\hskip 2 true cm Compute distance between templates at  new point: $(D\tau_0, D\tau_3)$
\hskip 1 true cm \}
}
\end{obeylines}

The algorithm to lay templates in the rest of the parameter space
is as follows:
\begin{obeylines}
\texttt{
\hskip 1 true cm Begin at $\tau_0 = \tau_0^{\rm min}, \tau_3 = \tau_3^{\rm min}$
\hskip 1 true cm Compute metric at $(\tau_0, \tau_3)$
\hskip 1 true cm Compute distance between templates at  new point: $(D\tau_0, D\tau_3)$
\hskip 1 true cm Add $(\tau_0, \tau_3)$ to InspiralTemplateList
\hskip 1 true cm numTemplates++
\hskip 1 true cm do while ($\tau_3 <= \tau_3^{\rm max}$)
\hskip 1 true cm \{
\hskip 2 true cm do while ($\tau_0 <= \tau_0^{\rm max}$)
\hskip 2 true cm \{
\hskip 3 true cm if ($(\tau_0, \tau_3)$ is inside the parameter space)
\hskip 3 true cm \{
\hskip 4 true cm Compute metric at ($\tau_0, \tau_3$)
\hskip 4 true cm Compute distance between templates at  new point: ($D\tau_0, D\tau_3$)
\hskip 4 true cm Add ($\tau_0, \tau_3$) to InspiralTemplateList
\hskip 4 true cm numTemplates++
\hskip 3 true cm \}
\hskip 3 true cm Increment $\tau_0:$ $\tau_0 = \tau_0 + D\tau_0$
\hskip 2 true cm \}
\hskip 2 true cm Increment $\tau_3:$ $\tau_3 = \tau_3 + D\tau_3$
\hskip 2 true cm Get next template along $\tau_3={\rm const.}$: $(\tau_0, \tau_3)$
\hskip 1 true cm \}
}
\end{obeylines}


\subsubsection*{Uses}
\begin{verbatim}
LALInspiralNextTemplate()
LALInspiralParameterCalc()
LALInspiralSetParams()
LALInspiralSetSearchLimits()
LALInspiralComputeParams()
LALInspiralComputeMetric()
LALInspiralUpdateParams()
LALInspiralValidTemplate()
\end{verbatim}

\subsubsection*{Notes}
\clearpage



\input{LALInspiralCreateFlatBankCP}
\idx{LALInspiralCreateFlatBank()}
\begin{itemize}
   \item \texttt{list,} Output, an array containing the template bank parameters
   \item \texttt{bankParams,} Input. It is necessary and sufficient to input
   the eigenvalues of the metric and the angle between the $x_0$ axis and the
   semi-major axis of the ambiguity ellipse, that is,
   \texttt{bankParams.metric.g00, bankParams.metric.g11, bankParams.metric.theta,}
   the minimal match, \texttt{bankParams.minimalMatch} and the range of the two
   coordinates over which templates must be chosen:
({\tt bankParams->x0Min}, {\tt bankParams->x0Max}) and
({\tt bankParams->x1Min}, {\tt bankParams->x1Max}).
\end{itemize}
The code expects {\tt list->vectorLength=2} and allocates just the
requisite amount of memory to {\tt list} and returns the number
of grid points in {\tt list->length.} The data points {\tt list->data[2j],}
{\tt j=1,2,\ldots, list->length,} contain the $x_0$-coordinates of the grid
and data points {\tt list->data[2j+1],} contain the $x_1$-coordinates
of the grid.

\subsubsection*{Description}
Given the {\tt metric} and the {\tt minimalMatch} this routine calls
{\tt bank/LALInspiralUpdateParams} to get the spacings in user coordinates (which are
not necessarily the eigen-directions) and lays a uniform grid of templates in
the range specified in ({\tt bankParams->x0Min}, {\tt bankParams->x0Max}) and
({\tt bankParams->x1Min}, {\tt bankParams->x1Max}).
\subsubsection*{Algorithm}
The algorithm to lay templates is as follows: Given the increments $Dx_0$ and
$Dx_1$ found from calling {\tt bank/LALInspiralUpdateParams} lay a rectangular
grid in the space of $(x_0, x_1).$
\begin{obeylines}
\texttt{
\hskip 1 true cm $x_1 = x_1^{\min}$
\hskip 1 true cm do while ($x_1 <= x_1^{\rm max}$)
\hskip 1 true cm \{
\hskip 2 true cm $x_0 = x_0^{\min}$
\hskip 2 true cm do while ($x_0 <= x_0^{\rm max}$)
\hskip 2 true cm \{
\hskip 3 true cm Add ($x_0, x_1$) to list
\hskip 3 true cm numTemplates++
\hskip 3 true cm Increment $x_0:$ $x_0 = x_0 + Dx_0$
\hskip 2 true cm \}
\hskip 2 true cm Increment $x_1:$ $x_1 = x_1 + Dx_1$
\hskip 1 true cm \}
}
\end{obeylines}

\subsubsection*{Uses}
\begin{verbatim}
LALInspiralUpdateParams()
LALRalloc()
\end{verbatim}

\subsubsection*{Notes}



\vspace{0.1in}
\vfill{\footnotesize\input{LALInspiralCreateCoarseBankCV}}
</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdio.h>
#include <lal/FindRoot.h>



NRCSID(LALINSPIRALCREATECOARSEBANKC, "$Id$");


/*  <lalVerbatim file="LALInspiralCreateCoarseBankCP"> */
void
LALInspiralCreateCoarseBank(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    )
{  /*  </lalVerbatim>  */
  INT4 i;

  INITSTATUS( status,
      "LALInspiralCreateCoarseBank", LALINSPIRALCREATECOARSEBANKC );
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

/* Anand:: 26 October 2006
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

  INITSTATUS( status, "LALNudgeTemplatesToConstantTotalMassLine",
      LALINSPIRALCREATECOARSEBANKC );
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

  INITSTATUS( status, "LALInspiralCreateCoarseBank",
      LALINSPIRALCREATECOARSEBANKC );
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


