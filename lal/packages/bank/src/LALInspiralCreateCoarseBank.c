/*  <lalVerbatim file="LALInspiralCreateCoarseBankCV">
Author: Churches, D. K and Sathyaprakash, B.S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralCreateCoarseBank.c}}

Module to create a coarse grid of templates.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralCreateCoarseBankCP}
\index{\verb&LALInspiralCreateCoarseBank()&}

\subsubsection*{Description}
The function first lays a set of templates along the $\eta=1/4$ curve
and then a grid in the $\tau_0$--$\tau_2$ space or $\tau_0$--$\tau_3$ space,
depending on the choice of \texttt{coarseIn.space=Tau0Tau2} or \texttt{Tau0Tau3}

\subsubsection*{Algorithm}


\subsubsection*{Uses}
\begin{verbatim}
LALInspiralNextTemplate()
LALInspiralSetParams()
LALInspiralSetSearchLimits()
LALInspiralComputeParams()
LALInspiralComputeMetric()
LALInspiralUpdateParams()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralCreateCoarseBankCV}}

</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALInspiralBank.h>


NRCSID (LALINSPIRALCREATECOARSEBANKC, "Id: $");

/*  <lalVerbatim file="LALInspiralCreateCoarseBankCP"> */

void LALInspiralCreateCoarseBank(LALStatus            *status, 
                                 InspiralTemplateList **list, 
                                 INT4                 *nlist,
                                 InspiralCoarseBankIn coarseIn) 
{  /*  </lalVerbatim>  */

   static InspiralBankParams bankPars, bankParsOld;
   static InspiralTemplate tempPars;
   static InspiralMetric metric;
   static INT4 validPars, i, pass;
   static REAL8 a25, x01, x02, x11, x12, dist1, dist2, ndx1, ndx2;

   INITSTATUS (status, "LALInspiralCreateCoarseBank", LALINSPIRALCREATECOARSEBANKC);
   ATTATCHSTATUSPTR(status);

   ASSERT (coarseIn.NoisePsd, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (coarseIn.mMin > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (coarseIn.MMax > 2.*coarseIn.mMin, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (coarseIn.mmCoarse > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (coarseIn.fLower > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (coarseIn.tSampling > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (coarseIn.tSampling >= 2.*coarseIn.fUpper, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (coarseIn.space >= 0, status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);
   ASSERT (coarseIn.space <= 1, status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);

/* Number of templates is nlist */
   *nlist = i = 0;

/* Set the elemets of the metric and tempPars structures in 
   conformity with the coarseIn structure */

   metric.NoisePsd = coarseIn.NoisePsd;
   metric.space = coarseIn.space;
   metric.iflso = coarseIn.iflso;
   LALInspiralSetParams(status->statusPtr, &tempPars, coarseIn);
   CHECKSTATUSPTR(status);

/* Identify the boundary of search and parameters for the first lattice point */
   LALInspiralSetSearchLimits(status->statusPtr, &bankPars, coarseIn);
   CHECKSTATUSPTR(status);
   LALInspiralComputeParams(status->statusPtr, &tempPars, bankPars, coarseIn);
   CHECKSTATUSPTR(status);

/* Compute the metric at this point, update bankPars and add it to the list */
   pass = 1;
   LALInspiralComputeMetric(status->statusPtr, &metric, tempPars, pass);
   CHECKSTATUSPTR(status);
   pass = coarseIn.iflso;
   LALInspiralUpdateParams(status->statusPtr, &bankPars, metric, coarseIn.mmCoarse);
   CHECKSTATUSPTR(status);
   *list = (InspiralTemplateList*) 
      LALRealloc(*list, sizeof(InspiralTemplateList)*(i+1));
   (*list)[i].ID = i; (*list)[i].params = tempPars; (*list)[i].metric = metric; ++i; 

/* First lay templates along the equal mass curve; i.e. eta=1/4. 
   Choose the constant and the index converting the chirp times to
   one another along the curve depending on whether the templates
   are laid along the tau0-tau2 or tau0-tau3 parameters
*/

   switch (coarseIn.space) {
      case Tau0Tau2:
         ndx1 = 0.6;
         ndx2 = 1./ndx1;
         a25 = pow(64./5., ndx1)*(2435./8064.)/pow(LAL_PI*coarseIn.fLower,.4);
         break;
      case Tau0Tau3:
         a25 = LAL_PI_2 * pow(64./5., .4)/pow(LAL_PI * coarseIn.fLower, .6);
         ndx1 = 0.4;
         ndx2 = 2.5;
         break;
   }
 
   bankParsOld = bankPars;
   while (bankPars.x0 < bankPars.x0Max) {
      x01 = bankPars.x0 + bankPars.dx0;
      x11 = a25 * pow(x01,ndx1);
      x12 = bankPars.x1 + bankPars.dx1;
      x02 = pow(x12/a25,ndx2);
      dist1 = pow(bankPars.x0 - x01,2.) + pow(bankPars.x1 - x11, 2.);
      dist2 = pow(bankPars.x0 - x02,2.) + pow(bankPars.x1 - x12, 2.);
      if (dist1 < dist2) {
         bankPars.x0 = x01;
         bankPars.x1 = x11;
      } else {
         bankPars.x0 = x02;
         bankPars.x1 = x12;
      }
/* If this is a valid point add it to our list */
   LALInspiralValidParams(status->statusPtr, &validPars, bankPars, coarseIn); 
   CHECKSTATUSPTR(status);
      if (validPars) {
         LALInspiralComputeParams(status->statusPtr, &tempPars, bankPars, coarseIn);
         CHECKSTATUSPTR(status);
         LALInspiralComputeMetric(status->statusPtr, &metric, tempPars, pass);
         CHECKSTATUSPTR(status);
         LALInspiralUpdateParams(status->statusPtr, &bankPars, metric, coarseIn.mmCoarse);
         CHECKSTATUSPTR(status);
         *list = (InspiralTemplateList*) 
            LALRealloc(*list, sizeof(InspiralTemplateList)*(i+1));
         (*list)[i].ID = i; (*list)[i].params = tempPars; (*list)[i].metric = metric; ++i; 
      }
   }

/* Begin with the parameters found at the first lattice point */
   bankPars = bankParsOld;
                      
/* Loop along x1 and x0 coordinates until maximum values are reached */
   while ( (bankPars.x1 += bankPars.dx1) <= bankPars.x1Max) {

/* Find the t0 coordinate of the next template close to the t2/t3 axis */
      LALInspiralNextTemplate(status->statusPtr, &bankPars, metric);
      CHECKSTATUSPTR(status);

/* Hop along t0-axis until t0 is inside the region of interest or quit */
      LALInspiralValidParams(status->statusPtr, &validPars, bankPars, coarseIn);
      CHECKSTATUSPTR(status);
      while (validPars==0 && bankPars.x0 < bankPars.x0Max) {
         bankPars.x0 += bankPars.dx0;
         LALInspiralValidParams(status->statusPtr, &validPars, bankPars, coarseIn);
         CHECKSTATUSPTR(status);
      }

/* If this is a valid point, compute the parameter and metric at 
   this point; add it to our list*/
      if (validPars) {
         LALInspiralComputeParams(status->statusPtr, &tempPars, bankPars, coarseIn); 
         CHECKSTATUSPTR(status);
         LALInspiralComputeMetric(status->statusPtr, &metric, tempPars, pass);  
         CHECKSTATUSPTR(status);
         LALInspiralUpdateParams(status->statusPtr, &bankPars, metric, coarseIn.mmCoarse);
         CHECKSTATUSPTR(status);
         *list = (InspiralTemplateList*) 
            LALRealloc(*list, sizeof(InspiralTemplateList)*(i+1));
         (*list)[i].ID = i; (*list)[i].params = tempPars; (*list)[i].metric = metric; ++i; 
      }

/* step along the tau0 axis until the boundary is reached */
      while ((bankPars.x0 += bankPars.dx0) <= bankPars.x0Max) {

/* If this is a valid point add it to our list */
         LALInspiralValidParams(status->statusPtr, &validPars, bankPars, coarseIn); 
         CHECKSTATUSPTR(status);
         if (validPars) {
            LALInspiralComputeParams(status->statusPtr, &tempPars, bankPars, coarseIn);
            CHECKSTATUSPTR(status);
            LALInspiralComputeMetric(status->statusPtr, &metric, tempPars, pass);
            CHECKSTATUSPTR(status);
            LALInspiralUpdateParams(status->statusPtr, &bankPars, metric, coarseIn.mmCoarse);
            CHECKSTATUSPTR(status);
            *list = (InspiralTemplateList*) 
               LALRealloc(*list, sizeof(InspiralTemplateList)*(i+1));
            (*list)[i].ID = i; (*list)[i].params = tempPars; (*list)[i].metric = metric; ++i; 
         }
      }
   } 
   *nlist = i;
   DETATCHSTATUSPTR(status);
   RETURN (status);
}

