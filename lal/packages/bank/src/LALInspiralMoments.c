/*  <lalVerbatim file="LALInspiralMomentsCV">
Authors: Brown, D. A., and Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

#if 0
<lalLaTeX>

\subsection{Module \texttt{LALInspiralMoments.c}}

Module to calculate the moment of the noise power spectral density. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralMomentsCP}
\idx{LALInspiralMoments()}

\subsubsection*{Description}

The moments of the noise curve are defined as
\begin{equation}
I(q)  \equiv S_{h}(f_{0}) \int^{f_{c}/f_{0}}_{f_{s}/f_{0}}
\frac{x^{-q/3}}{s_{h}(x f_{0})} \, dx \,.  
\end{equation}
Because in practice we will always divide one of these moments by another, we
do not need to include the $S_{h}(f_{0})$ term, which always cancels.
This function calculates the integral
\begin{equation}
I = \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-q/3}}{s_{h}(x f_{0})} \, dx \,.
\end{equation} 
It then divides this quantity by a normalisation constant which has been
passed to the function. In the case of calculating the components of the
metric for the signal manifold for the purpose of generating a template bank,
this constant is given by $I(7)$, because of the definition of the quantity
\begin{equation}
J(q) \equiv \frac{I(q)}{I(7)} \,.
\end{equation}

\subsubsection*{Algorithm}
Given the exponent \texttt{pars.ndx} and limits of integration
\texttt{pars.xmin} and \texttt{pars.xmax} this function returns the moment of
the power spectral density specified by by the frequency series
\texttt{pars.shf} according to
\begin{equation}
\mathtt{moment} = \int_{\mathtt{xmin}}^{\mathtt{xmax}} 
\frac{x^{-\mathtt{ndx}}}{S_h(x)}\, dx \, .
\end{equation}

\subsubsection*{Uses}
\begin{verbatim}
LALInspiralMomentsIntegrand
LALDRombergIntegrate
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralMomentsCV}}

</lalLaTeX>
#endif

#include <lal/LALInspiralBank.h>
#include <lal/Integrate.h>

NRCSID(LALINSPIRALMOMENTSC, "$Id$");

/*  <lalVerbatim file="LALInspiralMomentsCP"> */

void LALInspiralMoments(LALStatus         *status,
                        REAL8             *moment,
                        InspiralMomentsIn pars)
{ /* </lalVerbatim> */

   REAL8 fMax;
   DIntegrateIn In;
   void *funcparams = (void *) &pars;

   INITSTATUS (status, "LALInspiralMoments", LALINSPIRALMOMENTSC);
   ATTATCHSTATUSPTR(status);

   ASSERT (pars.shf, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (pars.shf->data, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (pars.shf->data->data, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   /* make sure that the minimum and maximum of the integral are within */
   /* the frequency series                                              */
   fMax = pars.shf->f0 + (REAL8) pars.shf->data->length * 
     pars.shf->deltaF;
   if ( pars.xmin < pars.shf->f0 || pars.xmax > fMax )
   {
     ABORT( status, 999, "limits outside range of frequency series" );
   }

   In.function = LALInspiralMomentsIntegrand;
   In.xmin = pars.xmin;
   In.xmax = pars.xmax;
   In.type = ClosedInterval;

   LALDRombergIntegrate (status->statusPtr, moment, &In, funcparams);
   CHECKSTATUSPTR(status);

   *moment /= pars.norm;
   DETATCHSTATUSPTR(status);
   RETURN (status);
}
