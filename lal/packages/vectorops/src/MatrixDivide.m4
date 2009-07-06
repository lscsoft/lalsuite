/******** <lalVerbatim file="MatrixDivideCV"> ********
Author: Tibbits, M. M.
$Id$
********* </lalVerbatim> ********/

/******************************* <lalLaTeX file="MatrixDivideC">
\subsection{Module \texttt{MatrixDivide.c}}
\label{s:MatrixDivide.c}

This file is dedicated to reproducing the matlab function " ./ ".  This file
has several declarations of the same function taking all forms of available
input.  This being said, I have yet to script the complex actions and their
counterparts.

\subsection*{Description}

This file is to help make the conversion from Matlab to c much earier.
In this file, we have created all of the versions of ./ that we plan on
using.

\subsection*{Algorithms}

The algorithm is the same as it is in matlab.  The dot in front of an operator
in matlab signifies that if either or both of the operands are vectors, then
the operation will be carried out member by member.  For instance

vector a[25];
vector b[25];
vector c[25];

c = a ./ b;

The result of this is:

c[0] =	a[0] /	b[0];
c[1] =	a[1] /	b[1];
.	.	.
.	.	.
.	.	.

etc.

\subsection*{Uses}

\begin{itemize}
\item \texttt{LALDCreateVector}
\end{itemize}

\subsection*{Notes}

At the current time none of the operations have been specified for neither the
complex datatypes nor the unsigned datatypes.

\subsubsection*{Prototypes}

 *************************************************** </lalLaTeX> */


#include <math.h>
#include "Matrix.h"

NRCSID( MATLABMATRIXDIVC, "$Id$");

define(`TYPECODE',`D')
define(`TYPECODE2',`D')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`D')
define(`TYPECODE2',`S')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`D')
define(`TYPECODE2',`I2')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`D')
define(`TYPECODE2',`I4')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`D')
define(`TYPECODE2',`I8')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`S')
define(`TYPECODE2',`S')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`S')
define(`TYPECODE2',`I2')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`S')
define(`TYPECODE2',`I4')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`S')
define(`TYPECODE2',`I8')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`I2')
define(`TYPECODE2',`I2')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`I4')
define(`TYPECODE2',`I2')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`I4')
define(`TYPECODE2',`I4')
include(`MatrixDivideBase.m4')


define(`TYPECODE',`I8')
define(`TYPECODE2',`I2')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`I8')
define(`TYPECODE2',`I4')
include(`MatrixDivideBase.m4')

define(`TYPECODE',`I8')
define(`TYPECODE2',`I8')
include(`MatrixDivideBase.m4')
