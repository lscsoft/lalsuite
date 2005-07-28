/* <lalVerbatim file="LALInspiralValidParamsCV">
Author: Churches, D. K. and Sathyaprakash, B.S.
$Id$
</lalVerbatim>  */


/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralValidParams.c}}
Module which checks whether or not a pair of parameter 
values are consistent with the search space.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralValidParamsCP}
\idx{LALInspiralValidParams()}
\begin{itemize}
   \item \texttt{valid,} Output, 0 means invalid template, 1 means valid
   \item \texttt{bankParams,} Input
   \item \texttt{coarseIn,} Input
\end{itemize}

Module which checks whether or not a pair of parameter 
values $\tau_{0}$ and $\tau_{2(3)}$ correspond to
a user specified range of component masses \texttt{(mMin,mMax)} OR to a
minimum value of the component masses \texttt{mMin} and maximum total
mass \texttt{MMax.} In the first case chirptimes satisfying the
constraint \texttt{mMin}~$\le m_1, m_2 \le$~\texttt{mMax} are accepted
as valid systems. In the second cases chirptimes satisfying the
constraint \texttt{mMin}~$\le m_1, m_2,$ and \texttt{MMax}$\le m=m_1+m_2$
are treated as valid. 

\subsubsection*{Description}

We start with the definition of the chirp times $\tau_{0}$ and $\tau_{3}$,
\begin{equation}
\tau_{0} = \frac{5}{256 (\pi f_{a} )^{8/3} m^{5/3} \eta}
\end{equation}
and
\begin{equation}
\tau_{3} = \frac{1}{8 (\pi^{2} f_{a}^{5} )^{1/3} m^{2/3} \eta}
\end{equation}
 These equations may be inverted to yield
\begin{equation}
m = \frac{5}{32 \pi^{2} f_{a}} \frac{\tau_{3}}{\tau_{0}}
\end{equation}
and
\begin{equation}
\eta = \left( \frac{2 \pi^{2}}{25 f_{a}^{3}} \frac{\tau_{0}^{2}}{\tau_{3}^{1/3}}
\right)^{5}\end{equation}

The individual masses may be calculated as follows.  We have
\begin{equation}
m = m_{1} + m_{2}
\label{mass}
\end{equation}
and
\begin{equation}
\eta = \frac{m_{1} m_{2}}{(m_{1} + m_{2})^{2}}
\label{eta}
\end{equation}
From Eq.(\ref{mass}) we may eliminate either $m_{1}$ or $m_{2}$,
\begin{equation}
m_{1} = m - m_{2}
\end{equation}
This may be substituted into Eq.(\ref{eta}) to give
\begin{equation}
\eta = \frac{(m - m_{2}) m_{2}}{\left[ (m - m_{2}) + m_{2} \right]^{2}}
     = \frac{(m - m_{2}) m_{2}}{m^{2}}
\end{equation}
which may be re--arranged to give
\begin{equation}
m_{2}^{2} - m m_{2} + \eta m^{2} = 0,
\end{equation}
i.e.\
\begin{equation}
m_{2} = \frac{ m \pm \sqrt{m^{2}(1 - 4 \eta) }}{2}
\end{equation}
Therefore, since we know that $\eta \leq 1/4$, real roots are guaranteed.
If we had eliminated $m_{2}$ rather than $m_{1}$ then we would have arrived at an identical
expression for
$m_{1}$, and so of one object has mass
\begin{equation}
m_{1} = \frac{m + \sqrt{m^{2}(1-4 \eta)}}{2}
\end{equation}
then the other object must have mass
\begin{equation}
m_{2} = \frac{m - \sqrt{m^{2}(1-4 \eta)}}{2}
\end{equation}
This function is also given \texttt{mMin} and \texttt{MMax} as inputs, which it may 
use to calculate the minimum value of $\eta$ which is possible with those inputs,
\begin{equation}
\eta_{min} = \mathtt{ \frac{mMin(MMax - mMin)}{MMax^{2}} }
\end{equation}

To recap, the function calculates $m$, $\eta$, $\eta_{min}$ and $m_{1,2}$.
It then checks that
\begin{equation}
\eta_{min} \leq \eta \leq 1/4
\end{equation}
and that
\begin{equation}
m_{1} \geq \mathtt{mMin}
\end{equation}
and
\begin{equation}
m_{2} \geq \mathtt{mMin}
\end{equation}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralValidParamsCV}}

</lalLaTeX>  */


#include <lal/LALInspiralBank.h>
#include <stdio.h>

NRCSID (LALINSPIRALVALIDPARAMSC, "$Id$");

/*  <lalVerbatim file="LALInspiralValidParamsCP">  */

void LALInspiralValidParams(
    LALStatus            *status,
    INT4                 *valid,
    InspiralBankParams   bankParams, 
    InspiralCoarseBankIn coarseIn
    )
/* </lalVerbatim> */
{

  InspiralTemplate *Pars=NULL;

  INITSTATUS( status, "LALInspiralValidParams", LALINSPIRALVALIDPARAMSC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( coarseIn.fLower > 0.L, status, 
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  *valid = 0;

  if ( bankParams.x0 <=0 || bankParams.x1 <=0 )
  {
    LALInfo( status, "x0 or x1 are less than or equal to zero" );
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }

  Pars = (InspiralTemplate *) LALCalloc( 1, sizeof(InspiralTemplate) );
  if ( ! Pars )
  {
    ABORT (status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
  }

  /* First set the chirp times of Pars to be as in bankParams */
  Pars->t0 = bankParams.x0;
  Pars->fLower = coarseIn.fLower;
  switch ( coarseIn.space ) 
  {
    case Tau0Tau2:
      Pars->t2 = bankParams.x1;
      Pars->massChoice = t02;
      break;
    case Tau0Tau3:
      Pars->t3 = bankParams.x1;
      Pars->massChoice = t03;
      break;
    default:
      ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
  }

  /* Compute all the parameters, including masses, */
  /* corresponding to (t0,t2/t3)                   */
  LALInspiralParameterCalc( status->statusPtr, Pars );
  CHECKSTATUSPTR( status );

  /* If the masses are in the correct range accept as valid parameters */
  switch (coarseIn.massRange) 
  {
    case MinComponentMassMaxTotalMass:
      if (
          Pars->mass1 >= coarseIn.mMin &&
          Pars->mass2 >= coarseIn.mMin &&
          Pars->totalMass <= coarseIn.MMax &&
          Pars->eta <= 0.25 && 
          Pars->eta >= coarseIn.etamin
         ) 
      {
        *valid = 1;
      }
      break;
    
    case MinMaxComponentMass:
      if (
          Pars->mass1 >= coarseIn.mMin &&
          Pars->mass2 >= coarseIn.mMin &&
          Pars->mass1 <= coarseIn.mMax &&
          Pars->mass2 <= coarseIn.mMax &&
          Pars->eta <= 0.25 && 
          Pars->eta >= coarseIn.etamin
         ) 
      {
        *valid = 1;
      }
      break;
      
    default:
      ABORT(status, 999, "Invalid choice for enum InspiralBankMassRange");
  }

  LALFree( Pars );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
