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

\vfill{\footnotesize\input{LALInspiralCreateCoarseBankCV}}

</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

static void
GetInspiralMoments (
		LALStatus            *status,
		InspiralMomentsEtc   *moments,
		REAL8FrequencySeries *psd,
		InspiralTemplate     *params );

NRCSID (LALINSPIRALCREATECOARSEBANKC, "Id: $");

/*  <lalVerbatim file="LALInspiralCreateCoarseBankCP"> */

void 
LALInspiralCreateCoarseBank(
		LALStatus            *status, 
		InspiralTemplateList **list, 
		INT4                 *nlist,
		InspiralCoarseBankIn coarseIn) 
{  /*  </lalVerbatim>  */

  InspiralBankParams bankPars, bankParsOld;
  InspiralTemplate *tempPars;
  InspiralMetric metric;
  InspiralMomentsEtc moments;
  INT4 validPars, i;
  REAL8 x01, x02, x11, x12, dist1, dist2, ndx1, ndx2, a25;

  INITSTATUS (status, "LALInspiralCreateCoarseBank", LALINSPIRALCREATECOARSEBANKC);
  ATTATCHSTATUSPTR(status);
  
  ASSERT (coarseIn.shf.data, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT (coarseIn.shf.data->data, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT (coarseIn.mMin > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT (coarseIn.mMax > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT (coarseIn.MMax >= 2.*coarseIn.mMin, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT (coarseIn.mmCoarse > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT (coarseIn.fLower > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT (coarseIn.tSampling > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT (coarseIn.tSampling >= 2.*coarseIn.fUpper, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT ((INT4)coarseIn.space >= 0, status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);
  ASSERT ((INT4)coarseIn.space <= 1, status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);
  
  ndx1 = 0.0;
  ndx2 = 0.0;
  a25 = 0.0;
/* Number of templates is nlist */
  *nlist = 0;

  /* Set the elements of the metric and tempPars structures in 
   conformity with the coarseIn structure
  */ 
  if (! (tempPars = (InspiralTemplate *)LALMalloc(sizeof(InspiralTemplate))) )
  {
	  ABORT (status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
  }
  metric.space = coarseIn.space;
  LALInspiralSetParams(status->statusPtr, tempPars, coarseIn);
  CHECKSTATUSPTR(status);

  /* Identify the boundary of search and parameters for the first lattice point */

  LALInspiralSetSearchLimits(status->statusPtr, &bankPars, coarseIn);
  CHECKSTATUSPTR(status);
  tempPars->totalMass = coarseIn.MMax;
  tempPars->eta = 0.25;
  tempPars->ieta = 1.L;
  tempPars->fLower = coarseIn.fLower;
  tempPars->massChoice = totalMassAndEta; 
  LALInspiralParameterCalc(status->statusPtr, tempPars);
  CHECKSTATUSPTR(status);

  /* Get the moments of the PSD integrand and other parameters
   required in the computation of the metric
  */
  GetInspiralMoments (status->statusPtr, &moments, &coarseIn.shf, tempPars);
  CHECKSTATUSPTR(status);


  /* compute the metric at this point, update bankPars and add the params to the list */
	  
  LALInspiralComputeMetric(status->statusPtr, &metric, tempPars, &moments);
  CHECKSTATUSPTR(status);
	  
  LALInspiralUpdateParams(status->statusPtr, &bankPars, metric, coarseIn.mmCoarse);
  /* printf("%e %e %e dx0=%e dx1=%e\n", metric.g00, metric.g11, metric.theta, bankPars.dx0, bankPars.dx1); */
  CHECKSTATUSPTR(status);
  
  /* Let's add the first template to the template list
   On failure realloc() returns a NULL to list, hence there is
   no need to explicitly free list 
  */
  if (!(*list = (InspiralTemplateList*) 
			  LALRealloc(*list, sizeof(InspiralTemplateList)*(*nlist+1)))) {
	  LALFree(tempPars);
	  ABORT (status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
  }
	  
  (*list)[*nlist].ID = *nlist; 
  (*list)[*nlist].params = *tempPars; 
  (*list)[*nlist].metric = metric; 
  ++(*nlist); 

  /* First lay templates along the equal mass curve; i.e. eta=1/4. 
   Choose the constant and the index converting the chirp times to
   one another along the curve depending on whether the templates
   are laid along the tau0-tau2 or tau0-tau3 space
  */

  switch (coarseIn.space) {
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
  }
  
  bankParsOld = bankPars;
  while (bankPars.x0 < bankPars.x0Max) {
    x01 = bankPars.x0 + bankPars.dx0;
    x11 = a25 * pow(x01,ndx1);
    x12 = bankPars.x1 + bankPars.dx1;
    x02 = pow(x12/a25,ndx2);
    dist1 = pow(bankPars.x0 - x01,2.L) + pow(bankPars.x1 - x11, 2.L);
    dist2 = pow(bankPars.x0 - x02,2.L) + pow(bankPars.x1 - x12, 2.L);
    if (dist1 < dist2) {
      bankPars.x0 = x01;
      bankPars.x1 = x11;
    } else {
      bankPars.x0 = x02;
      bankPars.x1 = x12;
    }
    /* If this is a valid point add it to our list */

    LALInspiralValidTemplate(status->statusPtr, &validPars, bankPars, coarseIn); 
    CHECKSTATUSPTR(status);
    if (validPars) {
      LALInspiralComputeParams(status->statusPtr, tempPars, bankPars, coarseIn);
      CHECKSTATUSPTR(status);
      LALInspiralComputeMetric(status->statusPtr, &metric, tempPars, &moments);
      CHECKSTATUSPTR(status);
      LALInspiralUpdateParams(status->statusPtr, &bankPars, metric, coarseIn.mmCoarse);
      /* printf("%e %e %e dx0=%e dx1=%e\n", metric.g00, metric.g11, metric.theta, bankPars.dx0, bankPars.dx1); */
      CHECKSTATUSPTR(status);
      /*
       On failure realloc() returns a NULL to outlist, hence there is
       no need to explicitly free the outlist 
      */
      if (!(*list = (InspiralTemplateList*) 
			      LALRealloc(*list, sizeof(InspiralTemplateList)*(*nlist+1)))) 
      {
	      LALFree(tempPars);
	      ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
      }
      (*list)[*nlist].ID = *nlist; 
      (*list)[*nlist].params = *tempPars; 
      (*list)[*nlist].metric = metric; 
      ++(*nlist); 
    }
  }
  
  /* Begin with the parameters found at the first lattice point */

  bankPars = bankParsOld;
  
  /* Loop along x1 and x0 coordinates until maximum values are reached */
  while (bankPars.x1 <= bankPars.x1Max) 
  {

    /* step along the tau0 axis until the boundary is reached */

    while (bankPars.x0 <= bankPars.x0Max) 
    {
    /* If this is a valid point add it to our list */

      LALInspiralValidTemplate(status->statusPtr, &validPars, bankPars, coarseIn); 
      CHECKSTATUSPTR(status);
      if (validPars) 
      {
	LALInspiralComputeParams(status->statusPtr, tempPars, bankPars, coarseIn);
	CHECKSTATUSPTR(status);
	LALInspiralComputeMetric(status->statusPtr, &metric, tempPars, &moments);
	CHECKSTATUSPTR(status);
	LALInspiralUpdateParams(status->statusPtr, &bankPars, metric, coarseIn.mmCoarse);
	/* printf("%e %e %e dx0=%e dx1=%e\n", metric.g00, metric.g11, metric.theta, bankPars.dx0, bankPars.dx1); */
	CHECKSTATUSPTR(status);
	if (!(*list = (InspiralTemplateList*) 
				LALRealloc(*list, sizeof(InspiralTemplateList)*(*nlist+1)))) 
        {
		LALFree(tempPars);
		ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
	}
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

    LALInspiralNextTemplate(status->statusPtr, &bankPars, metric);
    CHECKSTATUSPTR(status);

    /* Hop along t0-axis until t0 is inside the region of interest or quit */

    LALInspiralValidTemplate(status->statusPtr, &validPars, bankPars, coarseIn);
    CHECKSTATUSPTR(status);
    while (validPars==0 && bankPars.x0 < bankPars.x0Max) {
      bankPars.x0 += bankPars.dx0;
      LALInspiralValidTemplate(status->statusPtr, &validPars, bankPars, coarseIn);
      CHECKSTATUSPTR(status);
    }
    bankParsOld = bankPars;
  } 
  LALFree(tempPars);

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









static void
GetInspiralMoments (
		LALStatus            *status,
		InspiralMomentsEtc   *moments,
		REAL8FrequencySeries *psd,
		InspiralTemplate     *params )
{

   UINT4 k;
   InspiralMomentsIn in;

   INITSTATUS (status, "LALInspiralCreateCoarseBank", LALINSPIRALCREATECOARSEBANKC);
   ATTATCHSTATUSPTR(status);
  
   ASSERT (params, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (params->fLower>0, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (moments, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (psd, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   moments->a01 = 3.L/5.L;
   moments->a21 = 11.L * LAL_PI/12.L;
   moments->a22 = 743.L/2016.L * pow(25.L/(2.L*LAL_PI*LAL_PI), 1.L/3.L);
   moments->a31 = -3.L/2.L;
   moments->a41 = 617.L * LAL_PI * LAL_PI / 384.L;
   moments->a42 = 5429.L/5376.L * pow ( 25.L * LAL_PI/2.L, 1.L/3.L);
   moments->a43 = 1.5293365L/1.0838016L * pow(5.L/(4.L*pow(LAL_PI,4.L)), 1.L/3.L);
   
   /* setup the input structure needed in the computation of the moments */

   in.shf = psd;
   in.shf->f0 /= params->fLower;
   in.shf->deltaF /= params->fLower;
   in.xmin = params->fLower/params->fLower;
   in.xmax = params->fCutoff/params->fLower;
	   
   /* First compute the norm */

   in.norm = 1.L;
   in.ndx = 7.L/3.L; 
   LALInspiralMoments(status->statusPtr, &moments->j[7], in); 
   CHECKSTATUSPTR(status);
   in.norm = moments->j[7];

   if (lalDebugLevel & LALINFO)
   {
	   fprintf (stderr, "a01=%e a21=%e a22=%e a31=%e a41=%e a42=%e a43=%e \n", 
			   moments->a01, moments->a21, moments->a22, moments->a31, 
			   moments->a41, moments->a42, moments->a43);
   
	   fprintf(stderr, "j7=%e\n", moments->j[7]);
   }

   /* Normalised moments of the noise PSD from 1/3 to 17/3. */

   for (k=1; k<=17; k++)
   {
	   in.ndx = (REAL8) k /3.L; 
	   LALInspiralMoments(status->statusPtr,&moments->j[k],in);  
	   CHECKSTATUSPTR(status);
	   if (lalDebugLevel==1) fprintf(stderr, "j%1i=%e\n", k,moments->j[k]);
   }
   in.shf->deltaF *= params->fLower;
   in.shf->f0 *= params->fLower;
  
   DETATCHSTATUSPTR(status);
   RETURN (status);
}















/*  <lalVerbatim file="LALInspiralCreateCoarseBankCP"> */
void 
LALInspiralCreateFlatBank(
		LALStatus            *status, 
		REAL4VectorSequence  *list, 
		InspiralBankParams   *bankParams) 
{  /*  </lalVerbatim>  */

  InspiralMetric *metric; 
  REAL8 minimalMatch; 
  REAL8 x0, x1;
  UINT4 nlist = 0;

  INITSTATUS (status, "LALInspiralCreateCoarseBank", LALINSPIRALCREATECOARSEBANKC);
  ATTATCHSTATUSPTR(status);
  /* From the knowledge of the metric and the minimal match find the constant
   * increments bankParams->dx0 and bankParmams->dx1
   */
  metric = bankParams->metric;
  minimalMatch = bankParams->minimalMatch;
  LALInspiralUpdateParams(status->statusPtr, bankParams, *metric, minimalMatch);

  /* Construct the template bank */

  for (x1 = bankParams->x1Min; x1 <= bankParams->x1Max; x1 += bankParams->dx1)
  {
	  for (x0 = bankParams->x0Min; x0 <= bankParams->x0Max; x0 += bankParams->dx0)
	  {
		  UINT4 ndx = 2*nlist;
		  list->data = (REAL4 *)LALRealloc( list->data, (ndx+2)*sizeof(REAL4) );
		  if (!list->data)
		  {
			  ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
		  }
		  list->data[ndx] = x0;
		  list->data[ndx+1] = x1;
		  ++nlist; 
	  }
  }

  list->length = nlist;
  DETATCHSTATUSPTR(status);
  RETURN (status);
}
