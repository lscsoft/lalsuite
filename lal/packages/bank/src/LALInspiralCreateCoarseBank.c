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
#include <lal/FindRoot.h>


/* Thomas:: temporary definition for SPA hexagonal grid. */
static REAL4  A0;
static REAL4  A3;

typedef struct{
REAL4 ct;
REAL4 b;
}
PRIN;

static void LALSPAF(LALStatus *status,  REAL4 *result, REAL4 x, void *t3);
/* end comments ::Thomas */

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

      /* Thomas:: we can use either a square placement or an hexagonal
       * placement. The hexagonal placement is along the eigenvalues but 
       * not the square one.*/
      if (coarseIn.gridSpacing == Hexagonal){
	LALInspiralCreatePNCoarseBankHexa( status->statusPtr, list, nlist, coarseIn );
	CHECKSTATUSPTR( status );
      }
      else if (coarseIn.gridSpacing == SquareNotOriented){
	LALInspiralCreatePNCoarseBank( status->statusPtr, list, nlist, coarseIn );
	CHECKSTATUSPTR( status );
      }
      else {
        ABORT( status, LALINSPIRALBANKH_EGRIDSPACING, LALINSPIRALBANKH_MSGEGRIDSPACING );
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


void 
LALInspiralCreatePNCoarseBankHexa(
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    ) 
{  

  InspiralBankParams    bankPars;
  InspiralTemplate      *tempPars;
  InspiralMomentsEtc    moments;
  INT4                  i;
  InspiralCell          *cells;
  REAL4                 piFl;
  HexaGridParam         gridParam;
  CellEvolution         cellEvolution;
  INT4 firstId=0;
  CellList *cellList=NULL;

  INITSTATUS( status, "LALInspiralCreateCoarseBank", 
      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );


  ASSERT( coarseIn.mMin > 0., status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.mMax > 0., status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.MMax >= 2.*coarseIn.mMin, status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  /* Set the elements of the metric and tempPars structures in  */
  /* conformity with the coarseIn structure                     */ 
  if ( !(tempPars = (InspiralTemplate *) 
                LALCalloc( 1, sizeof(InspiralTemplate)))) {
    LALFree(tempPars);
    LALFree(cells);
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }

  
  LALInspiralSetParams( status->statusPtr, tempPars, coarseIn );
  CHECKSTATUSPTR( status );
  
  /* Identify the boundary of search and parameters for the     */
  /* first lattice point                                        */
  LALInspiralSetSearchLimits( status->statusPtr, &bankPars, coarseIn );
  CHECKSTATUSPTR( status );
  
  tempPars->totalMass   = coarseIn.MMax;
  tempPars->eta         = 0.25;
  tempPars->ieta        = 1.L;
  tempPars->fLower      = coarseIn.fLower;
  tempPars->massChoice  = m1Andm2;
  tempPars->mass1       = coarseIn.mMin;
  tempPars->mass2       = coarseIn.mMax;
  
  LALInspiralParameterCalc( status->statusPtr, tempPars );
  CHECKSTATUSPTR( status );
  
  /* Get the moments of the PSD integrand and other parameters */
  /* required in the computation of the metric  once for all.   */
  LALGetInspiralMoments( status->statusPtr, &moments, &coarseIn.shf, tempPars );
  CHECKSTATUSPTR( status );
  
  /* Allocate memory for one cell */
  cells = (InspiralCell*)
    LALCalloc(1,   sizeof(InspiralCell) );

  /*define gridParam*/
  gridParam.mm = coarseIn.mmCoarse;
  gridParam.x0Min     = bankPars.x0Min;
  gridParam.x0Max     = bankPars.x0Max;
  gridParam.x1Min     = bankPars.x1Min;
  gridParam.x1Max     = bankPars.x1Max;
  gridParam.mMin      = coarseIn.mMin;
  gridParam.mMax      = coarseIn.mMax;
  gridParam.etaMin    = coarseIn.etamin;
  gridParam.space     = coarseIn.space;

  cellEvolution.nTemplate = 1;
  cellEvolution.nTemplateMax = 1;
  cellEvolution.fertile = 0;

  /* initialise that first cell */

  tempPars->massChoice  = t03;
  cells[0].t0           = tempPars->t0;
  cells[0].t3           = tempPars->t3;

  /* some aliases */
  piFl  = LAL_PI * tempPars->fLower;
  A0    = 5. / pow(piFl, 8./3.) / 256.;
  A3    = LAL_PI / pow(piFl, 5./3.)/8.;


  LALCellInit(status->statusPtr, 
	      &cells, firstId, 
	      &moments, tempPars,
	      &gridParam, &cellEvolution, 
	      &cellList);
  CHECKSTATUSPTR( status );




  {
    INT4        k, kk;


    INT4       *list=NULL;
    CellList *ptr=NULL;
    INT4 length;


    while (cellEvolution.fertile) {
      length = Length(cellList);

      if (list!=NULL) 
	free(list);

      if (! (list =  LALMalloc(length*sizeof(INT4))))
      {
	ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
      ptr = cellList;

      for (k=0; k< length; k++){
	list[k]=ptr->id;	
	ptr = ptr->next;

      }      

      for (kk = 0; kk < length; kk++) 
	{	
	
	k = list[kk];

	if ( cells[k].status == Fertile) {

          LALPopulateCell(status->statusPtr, &moments, &cells,
              k,  tempPars, &gridParam, &cellEvolution, &cellList);          
	  CHECKSTATUSPTR( status );         	 	  


        }
      }
    }




  }

  if (cellList != NULL)
    printf("wierd behaviour here\n");



  *nlist = cellEvolution.nTemplate;

  {
    INT4 k ;
    INT4 length;
    length = cellEvolution.nTemplate;

    for (k=0; k<length; k++)
      {  
	REAL4  a,b, x0, tempA3;
	SFindRootIn input;
	INT4 valid;

	PRIN  prin;
	
	tempA3              = pow(A3, -5./2.)/pow(0.25,-1.5);
	tempPars->t0        = cells[k].t0;
	tempPars->t3        = cells[k].t3;
	
	if(cells[k].RectPosition[0] == Below ) {
	  
	  a = tan(cells[k].metric.theta);
	  b = cells[k].t3 - a * cells[k].t0;
	  
	  input.function = LALSPAF;
	  input.xmin = cells[k].t3-1e-3;
	  input.xmax = 1000;
	  input.xacc = 1e-6;
	  
	  prin.ct = a * A0 * tempA3;
	  prin.b = b;
	  
	  LALSBisectionFindRoot(status->statusPtr,&x0, &input, (void *)&prin);
	  CHECKSTATUSPTR( status );         
	  
	  tempPars->t3 = x0 + 1e-3; 
	  tempPars->t0 = (tempPars->t3 - b)/a;
	  if (tempPars->t0 > 0) {
	    LALInspiralParameterCalc(status->statusPtr, tempPars);
	    CHECKSTATUSPTR( status );         		  
	  }
	
	  cells[k].t0  = tempPars->t0;
	  cells[k].t3  = tempPars->t3;    
	  
	  /* update its position values */
	  valid = 1;
	  GetPositionRectangle(status->statusPtr, &cells, k,  tempPars , 
			       &gridParam, 
			       &cellEvolution, 
			       &cellList, 
			       &valid);
	  
	  {
	    INT4 above=0, below=0, in=0, out=0;
	    switch (cells[k].RectPosition[1]){
	    case In:    in    +=1; break;
	    case Below: below +=1; break;
	    case Above: above +=1; break;
	    case Out:   out   +=1; break;
	    }
	    switch (cells[k].RectPosition[2]){
	    case In:    in    +=1; break;
	    case Below: below +=1; break;
	    case Above: above +=1; break;
	    case Out:   out   +=1; break;
	    }
	    switch (cells[k].RectPosition[3]){
	    case In:    in    +=1; break;
	    case Below: below +=1; break;
	    case Above: above +=1; break;
	    case Out:   out   +=1; break;
	    }
	    switch (cells[k].RectPosition[4]){
	    case In:    in    +=1; break;
	    case Below: below +=1; break;
	    case Above: above +=1; break;
	    case Out:   out   +=1; break;
	    }
	    
	    if (above == 2 && cells[k].position == In){
	      cells[cells[k].child[0]].position = Out;
	    }
	    
	  }
		
	} 
		
      }
  }

  for (i=0; i<cellEvolution.nTemplate; i++) {
    if (cells[i].position == In ) {
      *nlist = *nlist +1; 
    }
  }

    /* allocate appropriate memory */

  
  *list = (InspiralTemplateList*) 
    LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist+1) );
  if ( ! *list )
      {
	LALFree( tempPars );
	ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
  memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );
  
  
  {
    *nlist = 0 ;
    for (i=0; i<cellEvolution.nTemplate; i++) {
      	if (cells[i].position == In) {
          tempPars->t0  = cells[i].t0;
          tempPars->t3  = cells[i].t3;
	    
          (*list)[*nlist].ID            = *nlist; 
          (*list)[*nlist].params        = *tempPars; 
          (*list)[*nlist].metric        = cells[i].metric; 
          ++(*nlist); 
        }
    }
  }
  
    
  LALFree( cells );
  LALFree( tempPars );
  
  DETATCHSTATUSPTR( status );
  RETURN ( status );
}
  




void
LALPopulateCell(LALStatus               *status,
		InspiralMomentsEtc      *moments,
		InspiralCell            **cell, 
		INT4                     headId,
		InspiralTemplate        *paramsIn,
		HexaGridParam           *gridParam,
		CellEvolution           *cellEvolution, 
		CellList **cellList
		)
{
  REAL4 dx0, dx1,  newt0, newt3;  
  INT4 i, id1, id2;
  REAL4 theta, ctheta,stheta;
  INT4 offSpring;
  INT4 it;
  INT4 add=0;

  INITSTATUS( status, "LALPopulateCell", 
	      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  /* aliases to get the characteristics of the parent template, that we refer
   * to its ID (headId) */  
  
  dx0           = (*cell)[headId].dx0;
  dx1           = (*cell)[headId].dx1;
  theta         = (*cell)[headId].metric.theta;
  ctheta        = cos(theta);
  stheta        = sin(theta);
  offSpring     = cellEvolution->nTemplate;

  
  /* Around the parent, the offspring can be at most 6 (hexagonal grid). 
   * By default the child are unset. If so it is created and have the 
   * properties of its parents. However, a child migh have been created 
   * earlier. In that case, we do not do anything.  */  
  it = 0 ; 

  for (i = 0; i < 6; i++) {
    if ((*cell)[headId].child[i] == -1) {
      add++;
      /* reallocate memory by set of 1000 cells if needed*/
      if ( (offSpring+add)>cellEvolution->nTemplateMax){
        *cell = (InspiralCell*) 
          LALRealloc( *cell, sizeof(InspiralCell) * (cellEvolution->nTemplateMax + 1000) );
        if ( ! cell ) {
          ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
        }
        cellEvolution->nTemplateMax +=  1000;
      }
      
      /* creates the child connection if needed. A child heritates the
       * properties of its parent */
   
      switch ( i ){
      case 0:
	newt0   = dx0 ;
	newt3   = 0 ;
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 *ctheta + stheta* newt3;
	(*cell)[offSpring + it].t3   += newt0 *stheta - ctheta* newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      case 1:
	newt0   =   dx0/2. ;
	newt3   =   -dx1 *sqrt(3./2) ;
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      case 2:
	newt0   =  -dx0/2 ;
	newt3   =  -dx1 *sqrt(3./2);
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      case 3:
	newt0   = -dx0 ;
	newt3   = 0;
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      case 4:
	newt0   =  -dx0/2. ;
	newt3   =  dx1 *sqrt(3./2);
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      case 5:
	newt0   = dx0/2. ;
	newt3   = dx1 *sqrt(3./2);
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList);
	break;
      }      
      
      /* Now, tricky part, if a child has been creating, he must have a
       * connection with its parents and vice-versa.  */
      if ((*cell)[offSpring + it].child[(i+3)%6] == -1){
	(*cell)[offSpring + it].child[(i+3)%6] = (*cell)[headId].ID;
	(*cell)[headId].child[i] = offSpring+it;
      }
      /* a new cell index */
      it += 1;
    }
  }
  
  /* how many new cells have been created ? */
  cellEvolution->nTemplate += it;


  /* Here, the parent has its 6 children set; he become sterile. */
  (*cell)[headId].status = Sterile;
  (cellEvolution->fertile)=cellEvolution->fertile-1;
  Delete(cellList, headId);


  /* what shall we do with that parent. Is he valid ? inside the space,
   * outside since eta>0.25 but close to the boundary .... */  
   	
  {
    if ((*cell)[headId].RectPosition[0]==Above && (*cell)[headId].in==1)
      {
	(*cell)[headId].RectPosition[0]=Out;
      }
  }
  

  
  /* propagate  connection annexe au freres pour eviter redondance */  
  for (i=0; i<6; i++){/* for each child*/
    id1 = (*cell)[headId].child[i%6];
    id2 = (*cell)[headId].child[(i+1)%6];
    (*cell)[id1].child[(i+2)%6] = (*cell)[id2].ID;
    (*cell)[id2].child[(i+4+1)%6] = (*cell)[id1].ID;   
  }
  
  
  /* enfin trouver position[0] (In/out)? of the children. */
  for (i=0; i<6; i++){/* for each child find position[0]*/
    id1 = (*cell)[headId].child[i%6];

    if ((*cell)[id1].status == Fertile) {
      LALSPAValidPosition(status->statusPtr, cell, id1, 
			  moments, cellEvolution, cellList);
      CHECKSTATUSPTR( status );
	  
      if ((*cell)[id1].position != In ) {
        if ((*cell)[id1].status == Fertile) {
          (*cell)[id1].status= Sterile;
          cellEvolution->fertile=cellEvolution->fertile-1;
	  Delete(cellList, id1);
        }
      }
    }
  }
    
  DETATCHSTATUSPTR( status );
  RETURN ( status );
}



void 
LALCellInit(    LALStatus               *status,
                InspiralCell            **cell, 
                INT4                    id,
                InspiralMomentsEtc      *moments, 
                InspiralTemplate        *paramsIn, 
		HexaGridParam           *gridParam, 
		CellEvolution           *cellEvolution,
		CellList **cellList)
{
  
  INT4          i;
  INT4 valid;   
  INITSTATUS( status, "LALCellInit", 
	      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );
  
  /* a new cell is created; by default it can create new children, 
     therefore it is fertile */
  cellEvolution->fertile = cellEvolution->fertile + 1;;
  (*cell)[id].status = Fertile;  
  append(cellList, id);


  /* all of whom are unset and do not have any id set*/
  for (i = 0; i < 6; i++) {
    (*cell)[id].child[i] = -1;
  } 
 
  /* filled some values related to the space */
  (*cell)[id].ID        = id;  
  (*cell)[id].position  = In;
  (*cell)[id].metric.space = gridParam->space;


  /* before ant further computation, check that t0, t3 is positive.*/
  if ((*cell)[id].t0 > 0 && (*cell)[id].t3 > 0){
    /* Get the metric at the position of the cell */ 
    paramsIn->t0 = (*cell)[id].t0;
    paramsIn->t3 = (*cell)[id].t3;

    LALInspiralComputeMetric( status->statusPtr, 
			      &((*cell)[id].metric),
			      paramsIn,
			      moments);
    CHECKSTATUSPTR( status );
  
    /* let us store the dx0 and dx3 at that point. */
    (*cell)[id].dx0 = sqrt(2.L * (1.L - gridParam->mm)/(*cell)[id].metric.g00 );
    (*cell)[id].dx1 = sqrt(2.L * (1.L - gridParam->mm)/(*cell)[id].metric.g11 );

    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[0]), paramsIn, gridParam);
    CHECKSTATUSPTR( status );

    /* if outside, this is a sterile cell which can not propagate */  
    if ((*cell)[id].RectPosition[0] == Out) {
      (*cell)[id].position      = Out;
      for (i = 0; i < 5; i++){
        (*cell)[id].RectPosition[i] = Out;
      }
      (*cell)[id].status = Sterile;
      (cellEvolution->fertile)=cellEvolution->fertile-1;
      Delete(cellList, id);
      
      
      DETATCHSTATUSPTR(status);
      RETURN(status);
    }
    else{
      valid = 1;
      GetPositionRectangle(status->statusPtr, &(*cell), id,  paramsIn , 
			   gridParam, cellEvolution, &(*cellList), &valid);
    }
  }
  else{/* if t0 or t3 < 0 , this is not a valid cell*/
    valid = 0;   
  }

  
  if (valid == 0){
    for (i=0; i<5; i++){(*cell)[id].RectPosition[i] = Out;}
    (*cell)[id].position = Out;
    (*cell)[id].status = Sterile;
    (cellEvolution->fertile)=cellEvolution->fertile-1;
    Delete(cellList, id);
  }


  DETATCHSTATUSPTR(status);
  RETURN(status);
}




void
GetPositionRectangle(LALStatus *status, 
		     InspiralCell **cell,
		     INT4 id,
		     InspiralTemplate *params, 
		     HexaGridParam *gridParam, 
		     CellEvolution *cellEvolution, 
		     CellList **cellList, 
		     INT4 *valid)
{

  RectangleIn   RectIn;
  RectangleOut  RectOut;
  InspiralTemplate paramsIn;



  INITSTATUS( status, "GetPosition", 
	      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  RectIn.x0    = params->t0;
  RectIn.y0    = params->t3;
  RectIn.dx    = (*cell)[id].dx0 ;
  RectIn.dy    = (*cell)[id].dx1 ;
  RectIn.theta = (*cell)[id].metric.theta;
  
  LALRectangleVertices(status->statusPtr, &RectOut, &RectIn);
  CHECKSTATUSPTR( status );

  paramsIn = *params;

  paramsIn.t0 = RectOut.x1;
  paramsIn.t3 = RectOut.y1;
  
  if (RectOut.x1>0 && RectOut.y1>0){
    LALFindPosition(status->statusPtr,(*cell)[id].dx0, (*cell)[id].dx1, 
		    &((*cell)[id].RectPosition[1]), 
		    &paramsIn, 
		    gridParam);    

    CHECKSTATUSPTR( status );
  }
  else {
    *valid = 0; 
    DETATCHSTATUSPTR(status);
    RETURN(status);     
  }
  
  paramsIn.t0 = RectOut.x2;
  paramsIn.t3 = RectOut.y2;
  if (RectOut.x2>0 && RectOut.y2>0){
    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[2]), &paramsIn, gridParam);
    CHECKSTATUSPTR( status );
    
  }
  else
    {
      *valid = 0;
      DETATCHSTATUSPTR(status);
      RETURN(status);     
    }
  
  paramsIn.t0 = RectOut.x3;
  paramsIn.t3 = RectOut.y3; 
  if (RectOut.x3>0 && RectOut.y3>0){
    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[3]), &paramsIn, gridParam);
    CHECKSTATUSPTR( status );
  }
  else
    {
      *valid = 0 ;
      DETATCHSTATUSPTR(status);
      RETURN(status);     
    }
  
  paramsIn.t0 = RectOut.x4;
  paramsIn.t3 = RectOut.y4;
  if (RectOut.x4>0 && RectOut.y4>0){
    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[4]), &paramsIn, gridParam); 
    CHECKSTATUSPTR( status );
  }
  else
    {
      *valid = 0;
      DETATCHSTATUSPTR(status);
      RETURN(status);     
    }
    
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

  if (coarseIn.LowGM  == -2)
  {
    LALInspiralCreateFlatBankS3( status->statusPtr, tempList, &bankParams , coarseIn);
    CHECKSTATUSPTR( status );
  }
  else 
  {
    LALInspiralCreateFlatBank( status->statusPtr, tempList, &bankParams);
    CHECKSTATUSPTR( status );
  }

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
         * one frequency at the nyquist frequency otherwise low masses
         * systems are missed. */
        if (((fMax * (1.L - (REAL4) (nf-1) * frac)) >= (*list)[j].params.tSampling/2.0)) 
         {
           fMax = (*list)[j].params.tSampling/2.0 - 1. ;
           frac = -1;
         }
         
        /*Similarly, for high masses. */
        /*if (((fMax * (1.L - (REAL4) (nf-1) * frac)) <= (*list)[j].params.fLower * 1.5)) 
         {
           fMax = (*list)[j].params.fLower * 1.5 ;           
         }
        */
        for (i=0; i<nf; i++)
        {
          fendBCV = fMax * (1.L - (REAL4) i * frac);


	    /* we search for valid expression of fendBCV and populate the bank */
	    if ( fendBCV >= (*list)[j].params.fLower * 1.5 && 
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
      params->totalMass = params->totalMass * 2.  ; /* The factor 2 is purely empiricail and 
						      comes from simulaitons. ?It seems indeed
						      tyhat the relation between psi0 and psi3 
						      which gives the total mass is not really
						      suitable. Ususally, the total mass is 
						      twice as much as the estimated one.
						   */
      params->fFinal = 1.L/( LAL_PI * pow(lightring, 1.5L) * params->totalMass );
      params->totalMass /= LAL_MTSUN_SI; /* it it used later ? */

    *valid = 1;

  }
}



/* <lalVerbatim file="LALInspiralCreateFlatBankS3CP"> */
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
  REAL8 x0, x1, dx1=0, dx0=0, x, y;
  UINT4 nlist = 0;
  INT4 layer  = 1;
  INT4 valid = -1;
  INITSTATUS( status, "LALInspiralCreateFlatBankS3", 
      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  /* From the knowledge of the metric and the minimal match 
     find the constant increments bankParams->dx0 and 
     bankParmams->dx1        */
  metric = bankParams->metric;
  minimalMatch = bankParams->minimalMatch;

  switch (coarseIn.gridSpacing){
  case Hexagonal:
    dx0 = sqrt(2.L * (1.L - minimalMatch)/metric->g00 );
    dx1 = sqrt(2.L * (1.L - minimalMatch)/metric->g11 );
    dx0 *=3./2./sqrt(2.);
    dx1 *=sqrt(3./2.);
    break;
  case Square:
    dx0 = sqrt(2.L * (1.L - minimalMatch)/metric->g00 );
    dx1 = sqrt(2.L * (1.L - minimalMatch)/metric->g11 );
    break;
  case  HexagonalNotOriented:
    LALInspiralUpdateParams( status->statusPtr, 
			     bankParams, *metric, minimalMatch );
    CHECKSTATUSPTR( status );
    dx0 = bankParams->dx0 * 3./2./sqrt(2.);
    dx1 = bankParams->dx1 * sqrt(3./2.);
    break;

  case  SquareNotOriented:
    LALInspiralUpdateParams( status->statusPtr, 
			     bankParams, *metric, minimalMatch );
    CHECKSTATUSPTR( status );
    dx0 = bankParams->dx0;
    dx1 = bankParams->dx1;
    break;
  }


  
  switch (coarseIn.gridSpacing){
  case Hexagonal:
  case HexagonalNotOriented:
    
    /* x1==psi3 and x0==psi0 */
    for (x1 = bankParams->x1Min -1e6;  x1 <= bankParams->x1Max + 1e6; x1 += dx1)
      {
	layer++;
	for (x0 = bankParams->x0Min - 1e6 +dx0/2.*(layer%2); x0 <= bankParams->x0Max+1e6; x0 += dx0 )
	  {
	    UINT4 ndx = 2 * nlist;
	
	    if ( coarseIn.gridSpacing == Hexagonal) 
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
                LALExcludeTemplate(status->statusPtr, &valid, bankParams, x, y);
                if (valid)
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
      }
      break;
  
  case  Square:
  case  SquareNotOriented:

    /* !! dx1 and dx0 are computed in a different way de[pending on the 
       value of BANKGRId */
    for (x1 = bankParams->x1Min -1e6;  x1 <= bankParams->x1Max + 1e6; x1 += dx1)
      {
	
	for (x0 = bankParams->x0Min - 1e6 ; x0 <= bankParams->x0Max+1e6; x0 += dx0 )

	  {
	    UINT4 ndx = 2 * nlist; 

	    if (coarseIn.gridSpacing == Square)
	      {
		x =  x0 *cos(metric->theta) + sin(metric->theta)* x1 ;
		y =  x0 *sin(metric->theta) - cos(metric->theta)* x1;
	      }
	    else if (coarseIn.gridSpacing == SquareNotOriented)
	      {
		x = x0;
		y = x1;
	      }
	    if ( (x > bankParams->x0Min - dx0/2.) && (y < bankParams->x1Max + dx1/2.) && 
		 (x < bankParams->x0Max + dx0/2.) && (y > bankParams->x1Min - dx1/2.))
	    
	      {
		LALExcludeTemplate(status->statusPtr, &valid, bankParams, x, y);
                if (valid)
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
   }
    break;
  }
  



  list->length = nlist;

  DETATCHSTATUSPTR(status);
  RETURN (status);
}



void
LALExcludeTemplate(
    LALStatus            *status, 
    INT4                 *valid, 
    InspiralBankParams   *bankParams,
    REAL4                 x,
    REAL4                 y)
{
  REAL4 psi0Int = 300000.;
  REAL4 psi3Int = -3000.;
  


  INITSTATUS( status, "LALExcludeTemplate", 
      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );
 

  if (x > psi0Int && y < psi3Int )
  {
    *valid = 0 ;
  }
  else 
  {
    *valid = 1;
  }

  DETATCHSTATUSPTR(status);
  RETURN (status);
}





void
LALSPAValidPosition(LALStatus *status, 
		    InspiralCell **cell,
		    INT4 id1,
		    InspiralMomentsEtc *moments, 
		    CellEvolution *cellEvolution, 
		    CellList **cellList
		    )
{
  INT4 below=0, in=0, out=0, above=0;
  
  INITSTATUS( status, "LALSPAFindPosition", 
	      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );


  switch ((*cell)[id1].RectPosition[1]){
  case In:    in    +=1; break;
  case Below: below +=1; break;
  case Above: above +=1; break;
  case Out:   out   +=1; break;
  }
  switch ((*cell)[id1].RectPosition[2]){
  case In:    in    +=1; break;
  case Below: below +=1; break;
  case Above: above +=1; break;
  case Out:   out   +=1; break;
  }
  switch ((*cell)[id1].RectPosition[3]){
  case In:    in    +=1; break;
  case Below: below +=1; break;
  case Above: above +=1; break;
  case Out:   out   +=1; break;
  }
  switch ((*cell)[id1].RectPosition[4]){
  case In:    in    +=1; break;
  case Below: below +=1; break;
  case Above: above +=1; break;
  case Out:   out   +=1; break;
  }
  switch ((*cell)[id1].RectPosition[0]){
  case In:    in    +=1; break;
  case Below: below +=1; break;
  case Above: above +=1; break;
  case Out:   out   +=1; break;
  }

  (*cell)[id1].in = in;
  
  if ((*cell)[id1].RectPosition[0]==In)
    {
      (*cell)[id1].position = In;
      if ((*cell)[id1].status == Sterile)
	{
	  (*cell)[id1].status = Fertile;
	  (cellEvolution->fertile)=cellEvolution->fertile+1;; 
	  append(cellList, id1);
	}
      DETATCHSTATUSPTR(status);
      RETURN(status);
  }




  if ( above == 5){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);

      }
  }
  else if ( below == 5){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);

      }
  }  
  else if ( out == 5){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);

      }
  }
  else if (in >= 1){
    (*cell)[id1].position = In;
    if ((*cell)[id1].status == Sterile)
      {
	(*cell)[id1].status = Fertile;
	(cellEvolution->fertile)=cellEvolution->fertile+1;
	append(cellList, id1);

  
      }

  }
  else if (above+below >= 5){
    if(out==1){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);

      }
    }
    else
    {
    (*cell)[id1].position = In;
    if ((*cell)[id1].status == Sterile)
      {
	(*cell)[id1].status = Fertile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);


      }
    }  
  }
  else{
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);

      }
  }

	    

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


void
LALFindPosition(LALStatus               *status, 
		REAL4                   dx0, 
		REAL4                   dx1,
		Position                *position, 
		InspiralTemplate        *paramsIn,
		HexaGridParam           *gridParam
)

{
  REAL8 mint3;  
  REAL4   eta, totalMass,ieta, oneby4, tiny, piFl;

  INITSTATUS( status, "LALFindPosition", 
	      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );


  ieta 	        = 1.;
  oneby4 	= 1./4.;
  tiny 	        = 1.e-10;
  piFl 	        = LAL_PI * paramsIn->fLower;
  
  ASSERT(paramsIn->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(paramsIn->t3 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  
  /* given t0, t3 we get the totalMass and eta. 
     We do not need to call ParameterCalc again and again here. */
  totalMass     = A0 * paramsIn->t3/(A3 * paramsIn->t0);
  eta           = A0/(paramsIn->t0 * pow(totalMass, fiveby3));
  
  /* be sure eta is inside the space if it is suppose to be */
  if (eta > oneby4) {
    eta-=tiny;
  }
  
  /* let us fill the param strucutre now : eta, Mass, mass1, mass2*/
  paramsIn->eta = eta;
  totalMass     = paramsIn->totalMass = totalMass/LAL_MTSUN_SI;
  if (eta <= oneby4) {
    paramsIn->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
    paramsIn->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
  }
  
  /* does t3 positive*/
  
  if ((paramsIn->t3-dx1)<0){ 
    mint3 = 0;
  }
  else{
    mint3 = paramsIn->t3-dx1;
  }
  
  if ( 
      (paramsIn->t0 <gridParam->x0Min - dx0/2)
      ||(paramsIn->t0 >gridParam->x0Max + dx0/2) 
          || (paramsIn->t3 <= mint3))
    {
      *position = Out;
      DETATCHSTATUSPTR(status);
      RETURN(status);
    }   

  if (
      paramsIn->mass1 >= gridParam->mMin &&
      paramsIn->mass2 >= gridParam->mMin &&
      paramsIn->mass1 <= gridParam->mMax &&
      paramsIn->mass2 <= gridParam->mMax &&
      paramsIn->eta <= 0.25 && 
      paramsIn->eta >= gridParam->etaMin
      ) 
    {
      *position = In;
    }
  else
    if (paramsIn->eta > .25){
      *position = Below; 
    }
    else{
      *position = Above;
    }    
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}





static void LALSPAF(LALStatus *status,  
		    REAL4 *result, 
		    REAL4 t3, 
		    void *param)
{
  REAL4 ct, b;
  PRIN *prin;

  INITSTATUS( status, "LALSPAF", 
	      LALINSPIRALCREATECOARSEBANKC );
  ATTATCHSTATUSPTR( status );

  prin = (PRIN *)param;
  ct = prin->ct;
  b  = prin->b;


  *result = ct*pow(t3,5./2.) - t3 + b;

  DETATCHSTATUSPTR( status );  
  RETURN(status);
}



void
print_list(CellList *head)
{
  if (head==NULL){
    printf("\n");
  }
  else { 
    printf(" %d", head->id);
    print_list(head->next);
  }
}


int Length(CellList *list)
{

  int count = 0; 

  while (list!=NULL){
    count++;
    list = list->next;
  }
  return count;
}


void append(CellList **headRef, INT4 id)
{
  CellList *current;


  if ((current = malloc(sizeof(*current))) == NULL) {
    {
      printf("Error with malloc\n");
      exit(0);
    }
  }

  current->id = id;
  current->next = *headRef;
  *headRef = current;

}




void DeleteList(CellList **headRef)
{
  CellList *tmp;
  
  while (headRef !=NULL){
    tmp = (*headRef)->next;
    free(headRef);
    (*headRef) = tmp;
  }
  
}



void Delete(CellList **headRef, INT4 id)
{
  CellList *ptr  = NULL;
  CellList *prev = NULL;


  for (ptr = *headRef; ptr != NULL; ptr= ptr->next){
    if (id == ptr->id)
      break;

    prev = ptr;
  }

  if (ptr==NULL)
    return ;

  if (prev!=NULL){
    prev->next = ptr->next;
  }
  else{
    *headRef = ptr->next;
  }

  /* free the data here if needed such as free(ptr->id); */
  free(ptr);

  return ;


 
}
