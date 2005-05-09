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


\subsection{Module \texttt{LALInspiralCreateBCVBank.c}}
Lay a flat grid of BCV templates in the user specified range
of the parameters $(\psi_0, \psi_3)$ in {\tt coarseIn} structure
(see below).
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralCreateBCVBankCP}
\idx{LALInspiralCreateBCVBank()}
\begin{itemize}
   \item \texttt{list,} Output, an array containing the template bank parameters.
   \item \texttt{nlist,} Output, the number of templates in bank.
\end{itemize}

\subsubsection*{Description}
Given the range of the parameters $(\psi_0, \psi_3),$  
number of templates in the {\tt fCut} direction,
{\it minimalMatch}, noise spectral density, upper and
lower frequency cutoffs (all in the input structure {\tt coarseIn})
this routine outputs the list of templates in the BCV bank
for the parameters $(\psi_0, \psi_3, f_{\rm cut}).$  
\subsubsection*{Algorithm}
A flat signal manifold is assumed and templates are laid
uniform in the three dimensions.  See below for an explanation
of how templates are chosen in the {\tt fcut} direction.

\subsubsection*{Uses}
\begin{verbatim}
LALInspiralUpdateParams()
LALRalloc()
\end{verbatim}

\subsubsection*{Notes}
\clearpage


\subsection{Module \texttt{LALInspiralCreateFlatBank.c}}
Lay a flat grid of templates in the user defined $(x_0, x_1)$ 
coordinates and range.
\subsubsection*{Prototypes}
\vspace{0.1in}
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
\clearpage

\subsection{Module \texttt{LALInspiralBCVFcutBank.c}}
Given a grid of templates with distinct values of $(\psi_0, \psi_3)$
this routine returns a new grid in which every template has {\tt numFcutTemplates}
partners differing from one another in the ending frequency {\tt fendBCV.}
A call to this function should be preceeded by a call to \texttt {LALInspiralCreateFlatBank.c},
or a similar function, that gives a grid in $(\psi_0, \psi_3)$ space.
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralBCVFcutBankCP}
\idx{LALInspiralBCVFcutBank()}
\begin{itemize}
   \item \texttt{list,} Output/Input, an array initially containing the template 
   bank with the values of {\tt list[j]->psi0, list[j]->psi3, list[j]->fLower,} specified,
   is replaced on return with a re-sized array specifying also {\tt list->fFinal.}
   \item \texttt{Nlist,} Output/Input, the number of templates in the Input bank is
   replaced by the number of templates in the output bank.
   \item \texttt{numFcutTemplates,} Input, the largest number of templates for the
   parameter $f_{cut}$ of BCV.
\end{itemize}

\subsubsection*{Description}

A lattice of templates for BCV models should include,
in addition to the values of $(\psi_0, \psi_3)$ 
a range of $f_{\rm cut}$ -- the cutoff frequency. 
The right approach would be
to compute the metric in the three-dimensional space of 
$(\psi_0, \psi_3, f_{\rm cut})$ and to choose templates as
dictated by the metric. However, analytic computation of the
metric has not been easy. Therefore, it has become necessary
(at least for the time being) to make alternate choice of
the cutoff frequencies. 

In this routine we implement a simple
choice based on physical grounds: The post-Newtonian models
predict an ending frequency that is larger than, but close to,
the Schwarzschild last-stable orbit frequency
$f_{\rm lso} = (6^{3/2} \pi M )^{-1}$ where $M$ is the total mass,
while the effective one-body model has an ending frequency close
to the light-ring, whose Schwarzschild value is 
$f_{\rm lr} = (3^{3/2} \pi M )^{-1}.$ It is necessary to know
the total mass of the system in both cases.  However, not all
pairs of $(\psi_0, \psi_3)$ can be inverted to get a positive
$M$ but only when $\psi_0 > 0$ and $\psi_3 < 0.$ Even then
it is not guaranteed that the symmetric mass ratio will be
less than $1/4,$ a necessary condition so that the component
masses are found to be real. However, we do not demand that the
symmetric mass ratio is less than a quarter. If the total mass 
is non-negative then we compute the $(f_{\rm lso}, f_{\rm lr})$
and choose a user specified {\tt numFcutTemplates} number of 
templates with their cutoff frequency {\tt list->fFinal} defined
uniformly spaced in the range $[f_{\rm lso},\ f_{\rm lr}].$

Furthermore, this routine discards all templates for which
either the mass is not defined or, when defined, {\tt list->fFinal} is
smaller than the user defined lower frequency cutoff or larger
than the Nyquist frequency of templates.
Thus, the number of templates returned by this routine could 
be larger or fewer than the input number of templates.

\subsubsection*{Algorithm}
Given $(\psi_0, \psi_3)$ one can solve for $(M, \eta)$ using:
\begin{equation}
M = \frac{-\psi_3}{16 \pi^2 \psi_0},\ \ \eta = \frac{3}{128 \psi_0 (\pi M)^{5/3}}.
\end{equation}
Given the total mass compute the last stable orbit and light-ring frequencies using
\begin{equation}
f_{\rm lso} = (6^{3/2} \pi M)^{-1},\ \  f_{\rm lr} = (3^{3/2} \pi M)^{-1}.
\end{equation}
Divide the range $(f_{\rm lso}, f_{\rm lr})$ so as to have $n_{\rm cut}={\tt numFcutTemplates}$
templates over this range: 
\begin{equation}
df = f_{\rm lr} \frac {\left ( 1 - 2^{-3/2} \right ) }{ (n_{\rm cut} -1) }.
\end{equation}
Next, choose templates at $f_k = f_{\rm lr} - k \times df,$ where $k=0, \ldots, n_{\rm cut}-1.$
Note that by definition $f_0 = f_{\rm lr}$ and $f_{n_{\rm cut}-1} = f_{\rm lso};$
there are exatly $n_{\rm cut}$ templates in the range $(f_{\rm lso}, f_{\rm lr}).$
We discard a template if either $M$ is not defined or if $f_{\rm cut}$ is smaller
than the lower frequency cutoff specified in  \texttt{list[j]->fLower.}

\subsubsection*{Uses}
\begin{verbatim}
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

static void 
PSItoMasses (
    InspiralTemplate *params, 
    UINT4 *valid,
    REAL4 highGM
);

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

    case TaylorT1: 
    case TaylorT2: 
    case TaylorT3: 
    case TaylorF1: 
    case TaylorF2: 
    case PadeT1: 
    case PadeF1: 
    case EOB: 
      ASSERT( coarseIn.space == Tau0Tau2 || coarseIn.space == Tau0Tau3, status,
          LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
      LALInspiralCreatePNCoarseBank( status->statusPtr, list, nlist, coarseIn );
      CHECKSTATUSPTR( status );
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


/*  <lalVerbatim file="LALInspiralCreateBCVBankCP"> */
void 
LALInspiralCreateBCVBank (
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    ) 
/*  </lalVerbatim>  */
{  
  INT4 j;
  INT4 nlistOld;
  static InspiralBankParams bankParams;
  static InspiralMetric metric;
  static InspiralTemplate params;
  static CreateVectorSequenceIn in; 
  static REAL4VectorSequence *tempList=NULL;

  INITSTATUS( status, "LALInspiralCreateBCVBank", 
      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.psi0Min > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.psi0Max > coarseIn.psi0Min, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.psi3Min < 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.psi3Max > coarseIn.psi3Min, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.LowGM < coarseIn.HighGM, status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  params.fLower = coarseIn.fLower;
  params.fCutoff = coarseIn.fUpper;
  params.alpha = coarseIn.alpha;

  LALInspiralComputeMetricBCV( status->statusPtr, 
      &metric, &coarseIn.shf, &params );
  CHECKSTATUSPTR( status );

  if ( lalDebugLevel & LALINFO ) 
  {
    REAL8 dx0 = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/metric.g00 );
    REAL8 dx1 = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/metric.g11 );
    LALPrintError( "G00=%e G01=%e G11=%e\n", 
        metric.G00, metric.G01, metric.G11 );
    LALPrintError( "g00=%e g11=%e theta=%e\n", 
        metric.g00, metric.g11, metric.theta );
    LALPrintError( "dp0=%e dp1=%e\n", dx0, dx1 );
  }

  bankParams.metric = &metric;
  bankParams.minimalMatch = coarseIn.mmCoarse;
  bankParams.x0Min = coarseIn.psi0Min;
  bankParams.x0Max = coarseIn.psi0Max;
  bankParams.x1Min = coarseIn.psi3Min;
  bankParams.x1Max = coarseIn.psi3Max;
  
  in.length = 1;
  in.vectorLength = 2;
  LALSCreateVectorSequence( status->statusPtr, &tempList, &in );
  CHECKSTATUSPTR( status );

  /*  LALInspiralCreateFlatBankS3( status->statusPtr, tempList, &bankParams , coarseIn);*/
  LALInspiralCreateFlatBank( status->statusPtr, tempList, &bankParams);
  CHECKSTATUSPTR( status );

  *nlist = tempList->length;
  *list = (InspiralTemplateList *) 
    LALCalloc( *nlist, sizeof(InspiralTemplateList) );
  if ( ! *list )
  {
    ABORT (status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
  }

  for ( j = 0; j < *nlist; ++j )
  {
    /* Retain only those templates that have meaningful chirptimes:*/
    (*list)[j].params.psi0 = (REAL8) tempList->data[2*j];
    (*list)[j].params.psi3 = (REAL8) tempList->data[2*j+1];
    (*list)[j].params.fLower = params.fLower;
    (*list)[j].params.nStartPad = 0;
    (*list)[j].params.nEndPad = 0;
    (*list)[j].params.tSampling= coarseIn.tSampling;
    (*list)[j].params.distance =  1.;
    (*list)[j].params.signalAmplitude= 1.;
    (*list)[j].params.approximant= BCV;
    (*list)[j].params.massChoice= psi0Andpsi3;
    (*list)[j].params.order= twoPN;
    (*list)[j].metric = metric;
    (*list)[j].params.alpha = coarseIn.alpha;
  }

  nlistOld = *nlist;
  /* If coarseIn.lowGM == - 1 then LowGM is  unphysical. Hence, we use a 
   * Regular grid in cutoff frequency which is independant of LowGM or HighGM
   * and which lays between Flower and Fsampling/2
   *
   * If coarseIn.LowGM > 0 it means that we want to use a bank where
   * the cutoff frequency has some physical values betwwen Low and
   * HighGM 
   * option -2 is used to try hexagonal grid and so on. In that case lowGM is
   * equal to 3 */
  if (coarseIn.LowGM  == -2)
    {
      LALInspiralBCVBankFcutS3( status->statusPtr, 
				list, nlist, coarseIn);
      CHECKSTATUSPTR( status );
    }
  else  if (coarseIn.LowGM  == -1)
  {
	  LALInspiralBCVRegularFcutBank( status->statusPtr, 
	      list, nlist, coarseIn);
	  CHECKSTATUSPTR( status );
  }
  else
  {
	  LALInspiralBCVFcutBank( status->statusPtr, 
	      list, nlist, coarseIn);
	  CHECKSTATUSPTR( status );
  }

  if ( lalDebugLevel & LALINFO ) 
  {
    LALPrintError( 
        "Templates before %d and after %d calling LALInspiralBCVBank\n", 
        nlistOld, *nlist );
  }

  LALSDestroyVectorSequence( status->statusPtr, &tempList );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="LALInspiralCreateFlatBankCP"> */
void 
LALInspiralCreateFlatBank (
    LALStatus            *status, 
    REAL4VectorSequence  *list, 
    InspiralBankParams   *bankParams
    )
/* </lalVerbatim> */
{  
  InspiralMetric *metric; 
  REAL8 minimalMatch; 
  REAL8 x0, x1;
  UINT4 nlist = 0;

  INITSTATUS( status, "LALInspiralCreateFlatBank", 
      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  /* From the knowledge of the metric and the minimal match find the */
  /* constant increments bankParams->dx0 and bankParmams->dx1        */
  metric = bankParams->metric;
  minimalMatch = bankParams->minimalMatch;
  LALInspiralUpdateParams( status->statusPtr, 
      bankParams, *metric, minimalMatch );
  CHECKSTATUSPTR( status );

  /* Construct the template bank */
  for (x1 = bankParams->x1Min; x1 <= bankParams->x1Max; x1 += bankParams->dx1)
  {
    for (x0 = bankParams->x0Min; x0 <= bankParams->x0Max; x0 += bankParams->dx0)
    {
      UINT4 ndx = 2 * nlist;
      list->data = (REAL4 *) LALRealloc( list->data, (ndx+2) * sizeof(REAL4) );
      if ( !list->data )
      {
        ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
      }
      list->data[ndx] = x0;
      list->data[ndx + 1] = x1;
      ++nlist; 
    }
  }

  list->length = nlist;

  DETATCHSTATUSPTR(status);
  RETURN (status);
}


/*  <lalVerbatim file="LALInspiralBCVFcutBankCP"> */
void 
LALInspiralBCVFcutBank (
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                *NList, 
    InspiralCoarseBankIn coarseIn
    ) 
/*  </lalVerbatim>  */
{  
  UINT4 nf, nlist, j, ndx;
  REAL8 frac, fendBCV;
  REAL4 LowGM, HighGM;

  INITSTATUS( status, "LALInspiralBCVFcutBank", LALINSPIRALCREATECOARSEBANKC );

  nf = coarseIn.numFcutTemplates;
  ndx = nlist = *NList;

  LowGM = coarseIn.LowGM;
  HighGM = coarseIn.HighGM;

  /* if we have only one layer, we don't need HighGM. 
   * And default value for LowGM is  3GM*/
  if ( nf == 1 )
  {
    frac = 1;
  }
  else
  {
    frac = (1.L - 1.L/pow(HighGM/3., 1.5L)) / (nf-1.L);
  }
  
  for ( j = 0; j < nlist; ++j )
  {
    UINT4 valid = 0;
    
    PSItoMasses( &((*list)[j].params), &valid , LowGM);
    
    if ( valid )
    {
      UINT4 i;
      REAL8 fMax; 

      fMax = (*list)[j].params.fFinal;

      for ( i = 0; i < nf; ++i )
      {
	fendBCV = fMax * (1.L - (REAL8) i * frac);

        if ( (*list)[j].params.tSampling <= 0 )
        {
          ABORT( status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
        }
        if ( fendBCV > (*list)[j].params.fLower && 
            fendBCV < (*list)[j].params.tSampling / 2.0 )
        {
          ++ndx;

	    *list = (InspiralTemplateList *) 
            LALRealloc( *list, ndx * sizeof(InspiralTemplateList) );
          if ( ! *list )
          {
            ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
          }
          memset( *list + ndx - 1, 0, sizeof(InspiralTemplate) );
          (*list)[ndx-1] = (*list)[j];
          (*list)[ndx-1].params.fFinal = fendBCV;
          (*list)[ndx-1].metric = (*list)[0].metric;
          (*list)[ndx-1].nLayer = i;

        }
      }
    }
  }
  for ( j = nlist; j < ndx; ++j )
  {
    (*list)[j-nlist] = (*list)[j];
  }

  *NList = ndx - nlist;
  *list = LALRealloc( *list, *NList * sizeof(InspiralTemplateList) );
  if ( ! *list )
  {
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }

  RETURN( status );
}


static void
PSItoMasses (
    InspiralTemplate *params,
    UINT4            *valid,
    REAL4             HighGM
    )
{
  if ( params->psi0 <= 0.L || params->psi3 >= 0.L )
  {
    *valid = 0;
  }
  else
  {
    REAL8 totalMass; 
    REAL8 eta; 
    REAL8 eightBy3 = 8.L/3.L;
    REAL8 twoBy3=2.L/3.L;
    REAL8 fiveBy3 = 5.L/3.L;

    params->totalMass = -params->psi3/(16.L*LAL_PI * LAL_PI * params->psi0);
    eta = params->eta = 
      3.L/(128.L * params->psi0 * pow(LAL_PI*params->totalMass, fiveBy3));
    totalMass = params->totalMass;
    params->fFinal = 1.L/( LAL_PI * pow(HighGM, 1.5L) * params->totalMass );
    params->totalMass /= LAL_MTSUN_SI;
    *valid = 1;

#if 0
    if (params->eta > 0.25L) 
    {
      *valid = 0;
    }
    else
    {
      LALInspiralParameterCalc( status->statusPtr, params );
      CHECKSTATUSPTR( status );
      *valid = 1;
    }
#endif

    params->t0 = 5.0 / ( 256.0*eta*pow(totalMass, fiveby3) * 
        pow(LAL_PI * params->fLower, eightBy3));
    params->t3 = LAL_PI/(8.0*eta*pow(totalMass, twoBy3) * 
        pow(LAL_PI * params->fLower, fiveBy3));
  }
}






void 
LALInspiralBCVBankFcutS3 (
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                *NList, 
    InspiralCoarseBankIn coarseIn
    ) 
/*  </lalVerbatim>  */
{  
  UINT4 nlist, j, ndx;
  REAL4 frac;
  REAL4 LowGM, HighGM;
  REAL4 fendBCV;
  INT4  nf;
  

  INITSTATUS( status, "LALInspiralBCVBankFcutS3", LALINSPIRALCREATECOARSEBANKC );

  nf    = coarseIn.numFcutTemplates;

  ndx   = nlist = *NList;

  LowGM         = coarseIn.LowGM  = 3.;
  HighGM        = coarseIn.HighGM;


  for ( j = 0; j < nlist; ++j )
  {
    UINT4 valid = 0;
    
    LALEmpiricalPSI2MassesConversion( &((*list)[j].params), &valid , LowGM);   
    
    if (valid)
      {
	UINT4 i;
	REAL8 fMax; 

	fMax = (*list)[j].params.fFinal; 
	if ( (*list)[j].params.tSampling <= 0 )
	  {
	    ABORT( status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
	  }
        /* the user might request only one layer */	
	if (nf == 1)
	    frac = 1;
	else 
	    frac = (1.L - 1.L/pow(HighGM/3., 1.5L)) / (nf-1.L);
	 
        /* sometimes since fMin is greater than the nyquist frequency, there
         * is no template generated. This is not acceptable. We need at least
         * one frequency at the nyquist frequency. */
        if (((fMax * (1.L - (REAL4) (nf-1) * frac)) >= (*list)[j].params.tSampling/2.0)) 
         {
           fMax = (*list)[j].params.tSampling/2.0 - 1. ;
           frac = -1;
         }
	
        for (i=0; i<nf; i++)
        {
          fendBCV = fMax * (1.L - (REAL4) i * frac);


	    /* we search for valid expression of fendBCV and populate the bank */
	    if ( fendBCV > (*list)[j].params.fLower * 1.5 && 
		 fendBCV < (*list)[j].params.tSampling / 2.0 )
	      {
		
		++ndx;
		
		*list = (InspiralTemplateList *) 
		  LALRealloc( *list, ndx * sizeof(InspiralTemplateList) );
		if ( ! *list )
		  {
		    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
		  }
		memset( *list + ndx - 1, 0, sizeof(InspiralTemplate) );
		(*list)[ndx-1] = (*list)[j];
		(*list)[ndx-1].params.fFinal = fendBCV;
		(*list)[ndx-1].metric = (*list)[0].metric;
		(*list)[ndx-1].nLayer = i;
		
	      }
	    
	  }
	
      }
  }
  
  
  for ( j = nlist; j < ndx; ++j )
  {
    (*list)[j-nlist] = (*list)[j];
  }

  *NList = ndx - nlist;
  *list = LALRealloc( *list, *NList * sizeof(InspiralTemplateList) );
  if ( ! *list )
  {
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }


  RETURN( status );
}



/*  <lalVerbatim file="LALInspiralBCVRegularFcutBankCP"> */
void 
LALInspiralBCVRegularFcutBank (
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                *NList, 
    InspiralCoarseBankIn coarseIn
    ) 
/*  </lalVerbatim>  */
{  
  /* no restriction of  physical masses. 
   * And regular layer of templates in the Frequency dimension */
  UINT4 i,nf, nlist, j, ndx;
  REAL8 fendBCV;

  INITSTATUS( status, "LALInspiralBCVFcutBank", LALINSPIRALCREATECOARSEBANKC );

  nf = coarseIn.numFcutTemplates;
  ndx = nlist = *NList;
  
  for ( j = 0; j < nlist; ++j )
  {     
      for ( i = 1; i <=nf; ++i )
      {
	fendBCV = (*list)[j].params.fLower 
		+ i * ((*list)[j].params.tSampling/2.0 - (*list)[j].params.fLower) / nf ;
        ++ndx;

	*list = (InspiralTemplateList *) 
        LALRealloc( *list, ndx * sizeof(InspiralTemplateList) );
        if ( ! *list )
          {
            ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
          }
         memset( *list + ndx - 1, 0, sizeof(InspiralTemplate) );
         (*list)[ndx-1] = (*list)[j];
         (*list)[ndx-1].params.fFinal = fendBCV;
         (*list)[ndx-1].metric = (*list)[0].metric;
         (*list)[ndx-1].nLayer = i;
      }
  }
    
  
  for ( j = nlist; j < ndx; ++j )
  {
    (*list)[j-nlist] = (*list)[j];
  }

 *NList = ndx - nlist;
  *list = LALRealloc( *list, *NList * sizeof(InspiralTemplateList) );
  if ( ! *list )
  {
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }

  RETURN( status );
}







void
LALEmpiricalPSI2MassesConversion (
    InspiralTemplate *params,
    UINT4            *valid,
    REAL4             lightring
    )
{
  if ( params->psi0 <= 0.L || params->psi3 >= 0.L )
    {
      *valid = 0;
    }
  else
    {
      params->totalMass = -params->psi3/(16.L*LAL_PI * LAL_PI * params->psi0);
      
      params->totalMass = params->totalMass * 2  ; /* The factor 2 is purely empiricail and 
						      comes from simulaitons. ?It seems indeed
						      tyhat the relation between psi0 and psi3 
						      which gives the total mass is not really
						      suitable. Ususally, the total mass is 
						      twice as much as the estimated one.
						   */
      params->fFinal = 1.L/( LAL_PI * pow(lightring, 1.5L) * params->totalMass );
      params->totalMass /= LAL_MTSUN_SI;

    *valid = 1;

  }
}



/* <lalVerbatim file="LALInspiralCreateFlatBankCP"> */
void 
LALInspiralCreateFlatBankS3 (
    LALStatus            *status, 
    REAL4VectorSequence  *list, 
    InspiralBankParams   *bankParams,
    InspiralCoarseBankIn coarseIn
    )
/* </lalVerbatim> */
{  
  InspiralMetric *metric; 
  REAL8 minimalMatch; 
  REAL8 x0, x1, dx1, dx0, x, y;
  UINT4 nlist = 0;
  INT4 layer = 1;

  INITSTATUS( status, "LALInspiralCreateFlatBankS3", 
      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  /* From the knowledge of the metric and the minimal match 
     find the constant increments bankParams->dx0 and 
     bankParmams->dx1        */
  metric = bankParams->metric;
  minimalMatch = bankParams->minimalMatch;

  switch (coarseIn.gridType){
  case  OrientedHexagonal:
    dx0 = sqrt(2.L * (1.L - minimalMatch)/metric->g00 );
    dx1 = sqrt(2.L * (1.L - minimalMatch)/metric->g11 );
    dx0 *=3./2./sqrt(2.);
    dx1 *=sqrt(3./2.);
    break;
  case OrientedSquare:
    dx0 = sqrt(2.L * (1.L - minimalMatch)/metric->g00 );
    dx1 = sqrt(2.L * (1.L - minimalMatch)/metric->g11 );
    break;
  case  Hexagonal:
    LALInspiralUpdateParams( status->statusPtr, 
			     bankParams, *metric, minimalMatch );
    CHECKSTATUSPTR( status );
    dx0 = bankParams->dx0 * 3./2./sqrt(2.);
    dx1 = bankParams->dx1 * sqrt(3./2.);
    break;

  case  Square:
    LALInspiralUpdateParams( status->statusPtr, 
			     bankParams, *metric, minimalMatch );
    CHECKSTATUSPTR( status );
    dx0 = bankParams->dx0;
    dx1 = bankParams->dx1;
    break;
  }


  
  switch (coarseIn.gridType){
  case OrientedHexagonal:
  case Hexagonal:
    
    /* x1==psi3 and x0==psi0 */
    for (x1 = bankParams->x1Min -1e6;  x1 <= bankParams->x1Max + 1e6; x1 += dx1)
      {
	layer++;
	for (x0 = bankParams->x0Min - 1e6 +dx0/2.*(layer%2); x0 <= bankParams->x0Max+1e6; x0 += dx0 )
	  {
	    UINT4 ndx = 2 * nlist;
	
	    if ( coarseIn.gridType == OrientedHexagonal) 
	      {
	    
		x =  x0 *cos(metric->theta) + sin(metric->theta)* x1;
		y =  x0 *sin(metric->theta) - cos(metric->theta)* x1;
	      }
	    else
	      {
		x = x0;
		y = x1;
	      }
	    
	    if ( (x > bankParams->x0Min -dx0/2.) && (y < bankParams->x1Max + dx1/2.) && 
		 (x < bankParams->x0Max +dx0/2.) && (y > bankParams->x1Min - dx1/2.))
	      {
		list->data = (REAL4 *) LALRealloc( list->data, (ndx+2) * sizeof(REAL4) );
		if ( !list->data )
		  {
		    ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
		  }
		list->data[ndx] = x;
		list->data[ndx + 1] = y;
		    ++nlist; 
	      }
	  }
      
      }
      break;
  
  case  OrientedSquare:
  case  Square:

    /* !! dx1 and dx0 are computed in a different way de[pending on the 
       value of BANKGRId */
    for (x1 = bankParams->x1Min -1e6;  x1 <= bankParams->x1Max + 1e6; x1 += dx1)
   {
	
	for (x0 = bankParams->x0Min - 1e6 ; x0 <= bankParams->x0Max+1e6; x0 += dx0 )

	  {
	    UINT4 ndx = 2 * nlist; 

	    if (coarseIn.gridType == OrientedSquare)
	      {
		x =  x0 *cos(metric->theta) + sin(metric->theta)* x1 ;
		y =  x0 *sin(metric->theta) - cos(metric->theta)* x1;
	      }
	    else if (coarseIn.gridType == Square)
	      {
		x = x0;
		y = x1;
	      }
	    if ( (x > bankParams->x0Min - dx0/2.) && (y < bankParams->x1Max + dx1/2.) && 
		 (x < bankParams->x0Max + dx0/2.) && (y > bankParams->x1Min - dx1/2.))
	    
	      {
		list->data = (REAL4 *) LALRealloc( list->data, (ndx+2) * sizeof(REAL4) );
		if ( !list->data )
		  {
		    ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
		  }
		list->data[ndx] = x;
		list->data[ndx + 1] = y;
		++nlist; 
	      }
	  }
      }
    break;
  }
  



  list->length = nlist;

  DETATCHSTATUSPTR(status);
  RETURN (status);
}
