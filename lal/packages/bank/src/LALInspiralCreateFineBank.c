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
\index{\verb&LALInspiralCreateFineBank()&}

\subsubsection*{Description}

A set of fine grid templates is chosen around a given coarse grid point
$(\tau_0, \tau_{2(3)})$ in a rectangular region of size 
$D\tau_0/2 \times D\tau_{2(3)}/2,$ where $D\tau_{0(2,3)}$ is the 
coarse grid spacing. 

\subsubsection*{Algorithm}

The spacing between fine grid templates is chosen
to be a constant determined by the metric at the coarse grid point; for
example, 
$$d\tau_0 = \sqrt{\frac{2 (1 - MM_{\rm Fine})}{g_{00}} }.$$
Only those grid points that are within the parameter space boundary are kept
and others are discarded.
 
\subsubsection*{Uses}
\begin{verbatim}
LALInspiralComputeParams()
LALInspiralUpdateParams()
LALInspiralValidParams()
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
 
  REAL8 x0FineMin, x1FineMin, x0FineMax, x1FineMax;     
  INT4  i, validPars; 
  static InspiralTemplate   tempPars;  
  static InspiralBankParams bankPars;


  INITSTATUS (status, "LALInspiralCreateFineBank", LALINSPIRALCREATEFINEBANKC);
  ATTATCHSTATUSPTR(status);
  ASSERT (fineIn.coarseIn.space>=0,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT (fineIn.coarseIn.space<=1,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT (fineIn.templateList.params.t0 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

  tempPars = fineIn.templateList.params;
  switch (fineIn.coarseIn.space) {
    case Tau0Tau2:
      ASSERT (fineIn.templateList.params.t2 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
      bankPars.x0 = fineIn.templateList.params.t0;
      bankPars.x1 = fineIn.templateList.params.t2;
      break;
    case Tau0Tau3:
      ASSERT (fineIn.templateList.params.t3 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
      bankPars.x0 = fineIn.templateList.params.t0;
      bankPars.x1 = fineIn.templateList.params.t3;
      break;
  }

  LALInspiralUpdateParams(status->statusPtr,&bankPars,fineIn.templateList.metric,fineIn.coarseIn.mmCoarse); 

  CHECKSTATUSPTR(status);

  x0FineMin = bankPars.x0 - bankPars.dx0/2.;
  x0FineMax = bankPars.x0 + bankPars.dx0/2.;
  x1FineMin = bankPars.x1 - bankPars.dx1/2.;
  x1FineMax = bankPars.x1 + bankPars.dx1/2.;

  LALInspiralUpdateParams(status->statusPtr,&bankPars,fineIn.templateList.metric,fineIn.coarseIn.mmFine); 
  CHECKSTATUSPTR(status);

  i=0;
    for(bankPars.x1=x1FineMin; bankPars.x1<=x1FineMax; bankPars.x1+=bankPars.dx1) {  
    for(bankPars.x0=x0FineMin; bankPars.x0<=x0FineMax; bankPars.x0+=bankPars.dx0) { 
	LALInspiralValidParams(status->statusPtr, &validPars, bankPars, fineIn.coarseIn);
      CHECKSTATUSPTR(status);
       if (validPars) {
          LALInspiralComputeParams(status->statusPtr, &tempPars, bankPars, fineIn.coarseIn);
          CHECKSTATUSPTR(status);
          *outlist = (InspiralTemplateList*) 
             LALRealloc(*outlist, sizeof(InspiralTemplateList)*(i+1));
          (*outlist)[i].params = tempPars;
          ++i; 
       }
    }   
    }
   *nlist = i;
   DETATCHSTATUSPTR(status);
   RETURN (status);
}
