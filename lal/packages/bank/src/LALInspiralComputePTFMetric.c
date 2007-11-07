/*
*  Copyright (C) 2007 Yi Pan, Duncan Brown
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

/*  <lalVerbatim file="LALInspiralComputePTFMetricCV">
Author: Yi Pan, Duncan Brown
$Id$
</lalVerbatim>  */

#if 0
<lalLaTeX>
\subsection{Module \texttt{LALInspiralComputePTFMetric.c}}

Module to compute the components of the metric which is used to describe
distances on Physical Template Family signal manifold.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralComputePTFIntrinsicMetricCP}
\idx{LALInspiralComputePTFIntrinsicMetric()}
\begin{itemize}
   \item \texttt{metric,} Output, the metric at the lattice point defined by \texttt{params}
   \item \texttt{psd,} Input, the power spectral density of the data
   \item \texttt{params,} Input, the parameters where metric must be computed
   in the computation of the metric.
\end{itemize}

\input{LALInspiralComputePTFFullMetricCP}
\idx{LALInspiralComputePTFFullMetric()}
\begin{itemize}
   \item \texttt{metric,} Output, the metric at the lattice point defined by \texttt{params}
   \item \texttt{psd,} Input, the power spectral density of the data
   \item \texttt{params,} Input, the parameters where metric must be computed
   in the computation of the metric.
\end{itemize}

\subsubsection*{Description}
We calculate the components of the metric using the procedure outlined 
by Yi.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\begin{verbatim}
LALMalloc
LALFree
\end{verbatim}

\subsubsection*{Notes}
 
\vfill{\footnotesize\input{LALInspiralComputePTFMetricCV}}

</lalLaTeX>
#endif

#include <stdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>

/* <lalVerbatim file="LALInspiralComputePTFIntrinsicMetricCP">  */
void XLALInspiralComputePTFIntrinsticMetric (
    InspiralMetric             *metric,
    REAL8FrequencySeries       *psd,
    InspiralTemplate           *params
    )
/* </lalVerbatim> */
{

}

/* <lalVerbatim file="LALInspiralComputePTFFullMetricCP">  */
void XLALInspiralComputePTFFullMetric (
    InspiralMetric             *metric,
    REAL8FrequencySeries       *psd,
    InspiralTemplate           *params
    )
/* </lalVerbatim> */
{

}

