/* <lalVerbatim file="LALInspiralValidTemplateCV">
Author: Churches, D. K.
$Id$
</lalVerbatim>  */


/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralValidTemplate.c}}

Module which checks whether or not a pair of parameter 
values $\tau_{0}$ and $\tau_{2(3)}$ correspond to
physical values for the masses of a binary and 
their symmetric mass ratio $\eta$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralValidTemplateCP}
\index{\verb&LALInspiralValidTemplate()&}

\subsubsection*{Description}

We start with the definition of the chirp times $\tau_{0}$ and $\tau_{3}$,
\begin{equation}
\tau_{0} = \frac{5}{256 (\pi f_{a} )^{8/3} m^{5/3} \eta}
\end{equation}

and

\begin{equation}
\tau_{3} = \frac{1}{8 (\pi^{2} f_{a}^{3} )^{1/3} m^{2/3} \eta}
\end{equation}
 These equations may be inverted to yield
\begin{equation}
m = \frac{5}{32 \pi^{2} f_{a}} \frac{\tau_{3}}{\tau_{0}}
\end{equation}

and

\begin{equation}
\eta = \left( \frac{2 \pi^{2}}{25 f_{a}^{3}} \frac{\tau_{0}^{2}}{\tau_{3}^{5}}
\right)^{5}\end{equation}

The individual masses may be calculated as follows.

We have

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
\eta = \frac{(m - m_{2}) m_{2}}{\left[ (m - m{2}) + m_{2} \right]^{2}}
\end{equation}
 
which may be re--arranged to give
 
\begin{equation}
m_{2}^{2} - m m_{2} + \eta m^{2} = 0
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

This function is also given \texttt{mMin} and \texttt{MMax} as inputs, which it may use to calculate the
minimum value of $\eta$ which is possible with those inputs,
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

\vfill{\footnotesize\input{LALInspiralValidTemplateCV}}

</lalLaTeX>  */


#include <lal/LALInspiralBank.h>
#include <stdio.h>

NRCSID (LALINSPIRALVALIDTEMPLATEC, "$Id$");

/*  <lalVerbatim file="LALInspiralValidTemplateCP">  */

void 
LALInspiralValidTemplate(
  LALStatus            *status,
  INT4                 *valid,
  InspiralBankParams   bankParams, 
  InspiralCoarseBankIn coarseIn)
{  /*  </lalVerbatim>  */


  INITSTATUS (status, "LALInspiralValidTemplate", LALINSPIRALVALIDTEMPLATEC);
  ATTATCHSTATUSPTR(status);
  ASSERT (bankParams.x0 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT (bankParams.x1 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  ASSERT (coarseIn.fLower > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

  *valid = 0;
/* 
  We have a valid template either if the template itself, or one
  of the vertices of the 'ambiguity rectangle', is in the region of
  interest
*/

  LALInspiralValidParams(status->statusPtr, valid, bankParams, coarseIn); 
  CHECKSTATUSPTR(status);

  if (*valid == 1) 
  {
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }

  bankParams.x1 = bankParams.x1 - bankParams.dx1/2.;
  LALInspiralValidParams(status->statusPtr, valid, bankParams, coarseIn); 
  CHECKSTATUSPTR(status);
  if (*valid == 1) 
  {
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }

/*
  bankParams.x0 = bankParams.x0 - 2.*bankParams.dx0;
  LALInspiralValidParams(status->statusPtr, valid, bankParams, coarseIn); 
  CHECKSTATUSPTR(status);
  if (*valid == 1) 
  {
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
  bankParams.x1 = bankParams.x1 + bankParams.dx1;
  LALInspiralValidParams(status->statusPtr, valid, bankParams, coarseIn); 
  CHECKSTATUSPTR(status);
  if (*valid == 1) 
  {
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
  bankParams.x1 = bankParams.x1 - 2.*bankParams.dx1;
  LALInspiralValidParams(status->statusPtr, valid, bankParams, coarseIn); 
  CHECKSTATUSPTR(status);
  if (*valid == 1) 
  {
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
*/

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
