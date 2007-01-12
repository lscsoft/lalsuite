/*  <lalVerbatim file="LALInspiralBCVBankCV">
Author: Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralBCVBank.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralBCVBankCP}
\idx{LALInspiralBCVBank()}


\subsubsection*{Prototypes}
\vspace{0.1in}
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
</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdio.h>
#include <lal/FindRoot.h>

static void 
PSItoMasses (
    InspiralTemplate *params, 
    UINT4 *valid,
    REAL4 highGM
);

NRCSID(LALINSPIRALBCVBANKC, "$Id$");

/*  <lalVerbatim file="LALInspiralBCVBankCP"> */
void 
LALInspiralCreateBCVBank (
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    ) 
{  /*  </lalVerbatim>  */
  INT4 j;
  INT4 nlistOld;
  static InspiralBankParams bankParams;
  static InspiralMetric metric;
  static InspiralTemplate params;
  static CreateVectorSequenceIn in; 
  static REAL4VectorSequence *tempList=NULL;

  INITSTATUS( status, "LALInspiralCreateBCVBank", 
      LALINSPIRALBCVBANKC );
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

  if (coarseIn.gridSpacing  != S2BCV)
  {
    LALInspiralCreateFlatBankS3S4( status->statusPtr, tempList, &bankParams , coarseIn);
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
   *  */
  if (coarseIn.gridSpacing != S2BCV)
    {
      LALInspiralBCVBankFcutS3S4( status->statusPtr, 
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
      LALINSPIRALBCVBANKC );
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

  INITSTATUS( status, "LALInspiralBCVFcutBank", LALINSPIRALBCVBANKC );

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






/*  <lalVerbatim file="LALInspiralBCVFcutBankCP"> */
void 
LALInspiralBCVBankFcutS3S4 (
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
  

  INITSTATUS( status, "LALInspiralBCVBankFcutS3S4", LALINSPIRALBCVBANKC );

  nf    = coarseIn.numFcutTemplates;

  ndx   = nlist = *NList;

  LowGM         =  3.;
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

  INITSTATUS( status, "LALInspiralBCVFcutBank", LALINSPIRALBCVBANKC );

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
      params->totalMass = params->totalMass * 2.  ; /* The factor 2 is purely empirical and 
						       comes from simulations. It seems indeed
						       that the relation between psi0 and psi3 
						       which gives the total mass is not really
						       suitable. Usually, the total mass is 
						       twice as much as the estimated one.
						    */
      params->fFinal = 1.L/( LAL_PI * pow(lightring, 1.5L) * params->totalMass );
      params->totalMass /= LAL_MTSUN_SI; 
      
      *valid = 1;
    }
}



/* <lalVerbatim file="LALInspiralCreateFlatBankS3S4CP"> */
void 
LALInspiralCreateFlatBankS3S4 (
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
  REAL4 xp[8] = {3000, 40000, 100000, 300000, 550000, 550000, 250000,3000};
  REAL4 yp[8] = {0, -3000, -3000, -2000, -1500, -300, -300, 0};
  INT4 npol = 8;

  
  INITSTATUS( status, "LALInspiralCreateFlatBankS3S4", 
      LALINSPIRALBCVBANKC );
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
		if (coarseIn.insidePolygon == True)
		{
                  LALInsidePolygon(status->statusPtr, xp, yp, npol, x, y, &valid);
		}
                else
		{
		  LALExcludeTemplate(status->statusPtr, &valid, bankParams, x, y);
		}



                if (valid == 1)
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
		if (coarseIn.insidePolygon == True) {
		  LALInsidePolygon(status->statusPtr, xp, yp, npol, x, y, &valid);
                }
		else
		{
		  LALExcludeTemplate(status->statusPtr, &valid, bankParams, x, y);
		}
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


/* Thomas: 31 Aug 2006. This function is redundant with the polygon fit. 
It was design for BBH and therefore had tight boundary. For a more general 
purpose, I extend the range to generous values
 */
void
LALExcludeTemplate(
    LALStatus            *status, 
    INT4                 *valid, 
    InspiralBankParams   *bankParams,
    REAL4                 x,
    REAL4                 y)
{
  REAL4 psi0Int = 2500000.;
  REAL4 psi3Int = -10000.;

  INITSTATUS( status, "LALExcludeTemplate", 
      LALINSPIRALBCVBANKC );
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
