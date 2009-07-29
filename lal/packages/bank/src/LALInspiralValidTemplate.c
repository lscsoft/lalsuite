/*
*  Copyright (C) 2007 David Churches, Peter Shawhan, Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* <lalVerbatim file="LALInspiralValidTemplateCV">
Author: Churches, D. K. and Sathyaprakash, B.S.
$Id$
</lalVerbatim>  */


/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralValidTemplate.c}}

Module which checks whether or not a given template should
be kept in the template list.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralValidTemplateCP}
\index{\verb&LALInspiralValidTemplate()&}
\begin{itemize}
   \item \texttt{valid,} Output, 0 means invalid template, 1 means valid
   \item \texttt{bankParams,} Input
   \item \texttt{coarseIn,} Input
\end{itemize}

\subsubsection*{Description}
Given the parameter values $\tau_{0}$ and $\tau_{2(3)}$ this code
checks to see if they correspond to physical values of the masses of
a binary and their symmetric mass ratio $\eta.$ The parameter values
will be accepted as valid parameters {\em even though} they
may not lie within the search space but their span does, as described below.
At the moment the code allows extra templates only
in the positive-$\tau_{2(3)}$ direction only. We have found
that placing templates in other directions is redundant.


\subsubsection*{Algorithm}
%-- The following is not an accurate description of the current algorithm:
%Compute the coordinates at the corners of the ambiguity rectangle.
%Accept if at least one of those points corresponds to the search space.

Consider the point $(\tau_0,\tau_{2(3)})$ describing the template, and
also a point at $(\tau_0,\tau_{2(3)}$$-$$\mbox{\texttt{bankParams.dx1/2}})$ ,
{\it i.e.}\ displaced in the negative $\tau_{2(3)}$ direction.
\texttt{bankParams.dx1} is calculated from the metric and corresponds
to the vertical spacing between the horizontal rows of templates being
considered.
Accept the template if at least one of those points is within the search space.

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


  INITSTATUS( status, "LALInspiralValidTemplate", LALINSPIRALVALIDTEMPLATEC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.fLower > 0, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

  *valid = 0;
  if ( bankParams.x0 <=0 || bankParams.x1 <=0 )
  {
    LALInfo( status, "x0 or x1 is less than or equal to zero" );
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }

  /* We have a valid template either if the template itself, or one     */
  /* of the vertices of the 'ambiguity rectangle', is in the region of  */
  /* interest                                                           */

  LALInspiralValidParams( status->statusPtr, valid, bankParams, coarseIn );
  CHECKSTATUSPTR( status );

  if ( *valid == 1 )
  {
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }

  bankParams.x1 = bankParams.x1 - bankParams.dx1/2.;

  LALInspiralValidParams( status->statusPtr, valid, bankParams, coarseIn );
  CHECKSTATUSPTR( status );

  if ( *valid == 1 )
  {
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }

#if 0
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
#endif

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
