/************************************ <lalVerbatim file="LALMomentCV">
Author: Tibbits, M M
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{LALMoment.c}}
\label{s:LALMoment.c}

Routine to compute various moments of data.

\subsubsection*{Prototypes}
\input{LALMomentCP}

\subsubsection*{Description}
The data structure passed in is either a REAL8 or a REAL4 Sequence.  The only parameter is which moment to calculate.
The function the sums the data, calculates the average and then it returns the average for the first moment, it returns the variance for the second moment, and it returns the n-th moment about the mean for higher order moments.

\subsubsection*{Algorithm}
\begin{itemize}
\item \textit{Find the mean (here referred to as $ \overline{x} $).}
\item \textit{Sum, over all the elements, the quantity: \((x[k] - \overline{x})^{n}. \)}
\item \textit{Divide the sum just made by N-1. Call it moment-n}
\item \textit{If n is greater than 2:}
\begin{itemize}
\item \textit{Sum, over all the elements, the quantity: \((x[k] - \overline{x})^{n}. \)}
\item \textit{Divide the sum just made by N. Call it moment-n}
\end{itemize}
\item \textit{Return moment-n}
\end{itemize}

\subsubsection*{Uses}

Determination of a specific moment of a set of data.

\subsubsection*{Notes}

\begin{itemize}
\item \textit{Moments less than two are not allowed.}
\item \textit{The result structure must be Non-NULL when passed in.}
\item \textit{The function assumes that the length member of the data passed in is correct.}
\end{itemize}

\vfill{\footnotesize\input{LALMomentCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALMoment.h>

NRCSID( LALMOMENTC, "$Id$");

define(`TYPECODE',`D')
include(`LALMomentBase.m4')

define(`TYPECODE',`S')
include(`LALMomentBase.m4')
