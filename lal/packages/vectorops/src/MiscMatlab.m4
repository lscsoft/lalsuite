/******** <lalVerbatim file="MiscMatlabCV"> ********
Author: Tibbits, M. M.
$Id$
********* </lalVerbatim> ********/

/******************************* <lalLaTeX file="MiscMatlabC">
\subsection{Module \texttt{MiscMatlab.c}}
\label{s:MiscMatlab.c}

This file reproduces the last few matlab functions that we needed for our purposes.
It creates useable forms of cumsum, sum, max, and finally an implemenation of the
array addressing in matlab.  Matlab has an easy of inverting a vector, (end: -1: 1)
and the final function, FlipVector returns a result vector that has been flipped in
that same manner.

\subsection*{Description}

This file reproduces the last few matlab functions that we needed for our purposes.
It creates useable forms of cumsum, sum, max, and finally an implemenation of the
array addressing in matlab.  Matlab has an easy of inverting a vector, (end: -1: 1)
and the final function, FlipVector returns a result vector that has been flipped in
that same manner.

\subsection*{Algorithms}

The algorithms are the same as in matlab.  Flip vector was discussed above.  Sum
takes the sum of all of the elements in a vector.  Cum sum takes an input vector:

vector input[25];
vector output[25];

output[0] = input[0];
output[1] = input[0] + input[1];
output[2] = input[0] + input[1] + input[2];

etc

\subsection*{Uses}

\begin{itemize}
\item \texttt{LALDCreateVector}
\end{itemize}

\subsection*{Notes}

At the current time none of the operations have been specified for neither the
complex datatypes nor the unsigned datatypes.

Also, the prototypes are out of order as I have used m4 to create all of the
functions from one codebase.

\subsubsection*{Prototypes}

 *************************************************** </lalLaTeX> */

#include <math.h>
#include "Matrix.h"

NRCSID( MATLABMATRIXSUMC, "$Id$");

define(`TYPECODE',`D')
include(`MiscMatlabBase.m4')

define(`TYPECODE',`S')
include(`MiscMatlabBase.m4')

define(`TYPECODE',`I2')
include(`MiscMatlabBase.m4')

define(`TYPECODE',`I4')
include(`MiscMatlabBase.m4')

define(`TYPECODE',`I8')
include(`MiscMatlabBase.m4')
