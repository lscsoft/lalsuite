/*  <lalVerbatim file="LALInspiralCreateFineBankCV">
Author: Sathyaprakash, B.S. and Churches, D. K. 
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralCreateFineBank.c}}

Function to create a fine grid of templates.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralCreateFineBankCP}
\idx{LALInspiralCreateFineBank()}
\begin{itemize}
   \item \texttt{outlist,} Output, containing an array of template bank parameters 
   \item \texttt{nlist,} Output, the number of fine bank templates around a given coarse-mesh point
   \item \texttt{fineIn,} Input, the parameters required to find the fine bank
\end{itemize}

\subsubsection*{Description}

The fine grid algorithm is a very simple algorithm that computes a uniform
grid of templates around a given coordinate point -- which can in particular be 
a coarse grid point -- from a knowledge of the metric at the coordinate point
and the coarse and fine grid minimal matches, $D\tau_{0,3}$ and
$d\tau_{0,3},$ respectively. Since $D\tau$ is not necessarily an
integral multiple of $d\tau$ the rectangular fine grid about the point
in question will be larger than required. The algorithm chooses templates
{\it symmetrically} about the given coarse grid point. It does so
by laying a rectangular lattice of templates with spacings 
$d\tau_0$ and $d\tau_3,$ in the rectangular region defined by 
$(\tau_0 - \Delta \tau_0, \tau_3 - \Delta \tau_3),$
$(\tau_0 + \Delta \tau_0, \tau_3 - \Delta \tau_3),$
$(\tau_0 + \Delta \tau_0, \tau_3 + \Delta \tau_3)$ and
$(\tau_0 - \Delta \tau_0, \tau_3 + \Delta \tau_3),$
where 
$$\Delta\tau_0 = d\tau_0 \left [ \frac{D\tau_0}{2d\tau_0} \right ], $$
and for any $x$, $[x]$ denotes the smallest integer greater than or
equal to $x$.
\begin{figure}[h]
\centering\includegraphics[angle=-90,width=4.5 true in]{LALInspiralBankHfine}
\caption{Algorithm sketching the construction of a rectangular fine grid around a given coordinate point.}
\label{fig:fine}
\end{figure}
The algorithm takes as input a structure of type 
\texttt{InspiralFineBankIn} and returns a \texttt{pointer-to-a-pointer} 
of type \texttt{InspiralTemplateList} as well as the number of fine grid
templates \texttt {int} around the lattice point in question.

The spacing between fine grid templates is chosen
to be a constant determined by the metric at the coarse grid point; for
example, 
$$d\tau_0 = \sqrt{\frac{2 (1 - MM_{\rm Fine})}{g_{00}} }.$$
Only those grid points that are within the parameter space boundary, or
have the vertices of the ambiguity rectangle inside the parameter
space, are kept and others are discarded.
 
\subsubsection*{Algorithm}

The Fine grid algorithm works as follows:
\begin{obeylines}
\texttt{
\hskip 1 true cm From input structure extract coordinates of the grid point $(\tau_0^G, \tau_3^G).$
\hskip 1 true cm Compute coarse and fine grid spacings $(D\tau_0, D\tau_3)$ and $(d\tau_0, d\tau_3)$
\hskip 1 true cm Compute half-sides of the {\it smallest} symmetric rectangle about $(\tau_0, \tau_3)$: 
\hskip 2 true cm $\Delta\tau_0 =  d\tau_0 {\it ceil}[D\tau_0/(2d\tau_0)],$ $\Delta\tau_3 =  d\tau_3 {\it ceil}[D\tau_3/(2d\tau_3)],$
\hskip 1 true cm Begin at $\tau_3 = \tau_3^G - \Delta \tau_3,$ 
\hskip 1 true cm do while ($\tau_3 <= \tau_3^G+\Delta \tau_3$)
\hskip 1 true cm \{
\hskip 2 true cm Begin at $\tau_0 = \tau_0^G - \Delta \tau_0,$ 
\hskip 2 true cm do while ($\tau_0 <= \tau_0^G+\Delta \tau_0$)
\hskip 2 true cm \{
\hskip 3 true cm if ($(\tau_0, \tau_3)$ is inside the parameter space)
\hskip 3 true cm \{
\hskip 4 true cm Add ($\tau_0, \tau_3$) to InspiralTemplateList
\hskip 4 true cm numTemplates++
\hskip 3 true cm \}
\hskip 3 true cm Increment $\tau_0:$ $\tau_0 = \tau_0 + d\tau_0$ 
\hskip 2 true cm \}
\hskip 2 true cm Increment $\tau_3:$ $\tau_3 = \tau_3 + d\tau_3$ 
\hskip 1 true cm \}
}
\end{obeylines}

\subsubsection*{Uses}
\begin{verbatim}
LALInspiralComputeParams()
LALInspiralUpdateParams()
LALInspiralValidTemplate()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralCreateFineBankCV}}

</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALInspiralBank.h>

NRCSID (LALINSPIRALCREATEFINEBANKC, "$Id$");

/*  <lalVerbatim file="LALInspiralCreateFineBankCP"> */

void LALInspiralCreateFineBank(LALStatus            *status,  
                               InspiralTemplateList **outlist,
                               INT4                 *nlist,
                               InspiralFineBankIn   fineIn)
{ /* </lalVerbatim> */
 
  REAL8 x0, x1, Dx0, Dx1, dx0, dx1, x0FineMin, x1FineMin, x0FineMax, x1FineMax;     
  INT4  i, j, validPars, bins0, bins1; 
  InspiralTemplate   *tempPars=NULL;  
  InspiralBankParams *bankPars=NULL;


  INITSTATUS (status, "LALInspiralCreateFineBank", LALINSPIRALCREATEFINEBANKC);
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
  x0FineMax = x0 + (float) bins0/2. * dx0;
  x1FineMin = x1 - (float) bins1/2. * dx1;
  x1FineMax = x1 + (float) bins1/2. * dx1;

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

