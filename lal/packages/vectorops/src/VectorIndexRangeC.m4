/* -*- C -*- */
/******** <lalVerbatim file="VectorIndexRangeCV"> ********
Author: David Chin <dwchin@umich.edu> +1-734-709-9119
$Id$
********* </lalVerbatim> ********/

/******************************* <lalLaTeX file="VectorIndexRangeC">
\subsection{Module \texttt{VectorIndexRange.c}}
\label{s:VectorIndexRange.c}


\subsection*{Description}

\subsection*{Algorithms}

\subsection*{Uses}

\begin{itemize}
\item \texttt{LAL*CreateVector}
\end{itemize}

\subsection*{Notes}

\subsubsection*{Prototypes}

 *************************************************** </lalLaTeX> */

#include <math.h>
#include "VectorIndexRange.h"

NRCSID( VECTORINDEXRANGEC, "$Id$" );

define(`TYPECODE',`CHAR')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`I2')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`I4')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`I8')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`U2')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`U4')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`U8')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`S')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`D')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`C')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`Z')
include(`VectorIndexRangeBaseC.m4')

/******************************* <lalLaTeX file="VectorIndexRangeC">

\vfill{\footnotesize\input{VectorIndexRangeCV}}

 *************************************************** </lalLaTeX> */
